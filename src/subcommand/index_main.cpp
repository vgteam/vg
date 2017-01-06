// index.cpp: define the "vg index" subcommand, which makes xg, GCSA2, and RocksDB indexes

#include <omp.h>
#include <unistd.h>
#include <getopt.h>

#include <string>
#include <vector>
#include <regex>

#include "subcommand.hpp"

#include "../vg.hpp"
#include "../index.hpp"
#include "../stream.hpp"
#include "../vg_set.hpp"
#include "../utility.hpp"
#include "../path_index.hpp"

#include "gcsa/gcsa.h"
#include "gcsa/algorithms.h"

using namespace std;
using namespace vg;
using namespace vg::subcommand;

void help_index(char** argv) {
    cerr << "usage: " << argv[0] << " index [options] <graph1.vg> [graph2.vg ...]" << endl
         << "Creates an index on the specified graph or graphs. All graphs indexed must " << endl
         << "already be in a joint ID space, and the graph containing the highest-ID node " << endl
         << "must come first." << endl
         << "xg options:" << endl
         << "    -x, --xg-name FILE     use this file to store a succinct, queryable version of" << endl
         << "                           the graph(s) (effectively replaces rocksdb)" << endl
         << "    -v, --vcf-phasing FILE import phasing blocks from the given VCF file as threads" << endl
         << "    -r, --rename V=P       rename contig V in the VCFs to path P in the graph (may repeat)" << endl
         << "    -T, --store-threads    use gPBWT to store the embedded paths as threads" << endl
         << "gcsa options:" << endl
         << "    -g, --gcsa-out FILE    output a GCSA2 index instead of a rocksdb index" << endl
         << "    -i, --dbg-in FILE      optionally use deBruijn graph encoded in FILE rather than an input VG (multiple allowed" << endl
         << "    -k, --kmer-size N      index kmers of size N in the graph" << endl
         << "    -X, --doubling-steps N use this number of doubling steps for GCSA2 construction" << endl
         << "    -Z, --size-limit N     limit of memory to use for GCSA2 construction in gigabytes" << endl
         << "    -O, --path-only        only index the kmers in paths embedded in the graph" << endl
         << "    -F, --forward-only     omit the reverse complement of the graph from indexing" << endl
         << "    -e, --edge-max N       only consider paths which make edge choices at <= this many points" << endl
         << "    -j, --kmer-stride N    step distance between succesive kmers in paths (default 1)" << endl
         << "    -d, --db-name PATH     create rocksdb in PATH directory (default: <graph>.index/)" << endl
         << "                           or GCSA2 index in PATH file (default: <graph>" << gcsa::GCSA::EXTENSION << ")" << endl
         << "                           (this is required if you are using multiple graphs files)" << endl
         << "    -t, --threads N        number of threads to use" << endl
         << "    -p, --progress         show progress" << endl
         << "    -V, --verify-index     validate the GCSA2 index using the input kmers (important for testing)" << endl
         << "rocksdb options (ignored with -g):" << endl
         << "    -s, --store-graph      store graph as xg" << endl
         << "    -m, --store-mappings   input is .gam format, store the mappings in alignments by node" << endl
         << "    -a, --store-alignments input is .gam format, store the alignments by node" << endl
         << "    -A, --dump-alignments  graph contains alignments, output them in sorted order" << endl
         << "    -N, --node-alignments  input is (ideally, sorted) .gam format, cross reference nodes by alignment traversals" << endl
         << "    -P, --prune KB         remove kmer entries which use more than KB kilobytes" << endl
         << "    -n, --allow-negs       don't filter out relative negative positions of kmers" << endl
         << "    -D, --dump             print the contents of the db to stdout" << endl
         << "    -M, --metadata         describe aspects of the db stored in metadata" << endl
         << "    -L, --path-layout      describes the path layout of the graph" << endl
         << "    -S, --set-kmer         assert that the kmer size (-k) is in the db" << endl
        //<< "    -b, --tmp-db-base S    use this base name for temporary indexes" << endl
         << "    -C, --compact          compact the index into a single level (improves performance)" << endl
         << "    -o, --discard-overlaps if phasing vcf calls alts at overlapping variants, call all but the first one as ref" << endl;
}
#define debug
int main_index(int argc, char** argv) {

    if (argc == 2) {
        help_index(argv);
        return 1;
    }

    string rocksdb_name;
    string gcsa_name;
    string xg_name;
    // Where should we import haplotype phasing paths from, if anywhere?
    string vcf_name;
    // This maps from graph path name (FASTA name) to VCF contig name
    map<string,string> path_to_vcf;
    vector<string> dbg_names;
    int kmer_size = 0;
    bool path_only = false;
    int edge_max = 0;
    int kmer_stride = 1;
    int prune_kb = -1;
    bool store_graph = false;
    bool dump_index = false;
    bool describe_index = false;
    bool show_progress = false;
    bool set_kmer_size = false;
    bool path_layout = false;
    bool store_alignments = false;
    bool store_node_alignments = false;
    bool store_mappings = false;
    bool allow_negs = false;
    bool compact = false;
    bool dump_alignments = false;
    int doubling_steps = 3;
    bool verify_index = false;
    bool forward_only = false;
    size_t size_limit = 200; // in gigabytes
    bool store_threads = false; // use gPBWT to store paths
    bool discard_overlaps = false;

    int c;
    optind = 2; // force optind past command positional argument
    while (true) {
        static struct option long_options[] =
        {
            //{"verbose", no_argument,       &verbose_flag, 1},
            {"db-name", required_argument, 0, 'd'},
            {"kmer-size", required_argument, 0, 'k'},
            {"doubling-steps", required_argument, 0, 'X'},
            {"edge-max", required_argument, 0, 'e'},
            {"kmer-stride", required_argument, 0, 'j'},
            {"store-graph", no_argument, 0, 's'},
            {"store-alignments", no_argument, 0, 'a'},
            {"dump-alignments", no_argument, 0, 'A'},
            {"store-mappings", no_argument, 0, 'm'},
            {"dump", no_argument, 0, 'D'},
            {"metadata", no_argument, 0, 'M'},
            {"set-kmer", no_argument, 0, 'S'},
            {"threads", required_argument, 0, 't'},
            {"progress",  no_argument, 0, 'p'},
            {"prune",  required_argument, 0, 'P'},
            {"path-layout", no_argument, 0, 'L'},
            {"compact", no_argument, 0, 'C'},
            {"allow-negs", no_argument, 0, 'n'},
            {"gcsa-name", required_argument, 0, 'g'},
            {"xg-name", required_argument, 0, 'x'},
            {"vcf-phasing", required_argument, 0, 'v'},
            {"rename", required_argument, 0, 'r'},
            {"verify-index", no_argument, 0, 'V'},
            {"forward-only", no_argument, 0, 'F'},
            {"size-limit", no_argument, 0, 'Z'},
            {"path-only", no_argument, 0, 'O'},
            {"store-threads", no_argument, 0, 'T'},
            {"node-alignments", no_argument, 0, 'N'},
            {"dbg-in", required_argument, 0, 'i'},
            {"discard-overlaps", no_argument, 0, 'o'},
            {0, 0, 0, 0}
        };

        int option_index = 0;
        c = getopt_long (argc, argv, "d:k:j:pDshMt:b:e:SP:LmaCnAg:X:x:v:r:VFZ:Oi:TNo",
                long_options, &option_index);

        // Detect the end of the options.
        if (c == -1)
            break;

        switch (c)
        {
        case 'd':
            rocksdb_name = optarg;
            break;

        case 'x':
            xg_name = optarg;
            break;

        case 'v':
            vcf_name = optarg;
            break;
            
        case 'r':
            {
                // Parse the rename old=new
                string key_value(optarg);
                auto found = key_value.find('=');
                if (found == string::npos || found == 0 || found + 1 == key_value.size()) {
                    cerr << "error:[vg construct] could not parse rename " << key_value << endl;
                    exit(1);
                }
                // Parse out the two parts
                string vcf_contig = key_value.substr(0, found);
                string graph_contig = key_value.substr(found + 1);
                // Add the name mapping
                path_to_vcf[graph_contig] = vcf_contig;
            }
            break;

        case 'P':
            prune_kb = atoi(optarg);
            break;

        case 'k':
            kmer_size = atoi(optarg);
            break;


        case 'O':
            path_only = true;
            break;

        case 'e':
            edge_max = atoi(optarg);
            break;


        case 'j':
            kmer_stride = atoi(optarg);
            break;

        case 'p':
            show_progress = true;
            break;

        case 'D':
            dump_index = true;
            break;

        case 'M':
            describe_index = true;
            break;

        case 'L':
            path_layout = true;
            break;

        case 'S':
            set_kmer_size = true;
            break;

        case 's':
            store_graph = true;
            break;

        case 'a':
            store_alignments = true;
            break;

        case 'A':
            dump_alignments = true;
            break;

        case 'm':
            store_mappings = true;
            break;

        case 'n':
            allow_negs = true;
            break;

        case 'C':
            compact = true;
            break;

        case 't':
            omp_set_num_threads(atoi(optarg));
            break;

        case 'g':
            gcsa_name = optarg;
            break;

        case 'V':
            verify_index = true;
            break;
        case 'i':
            dbg_names.push_back(optarg);
            break;
        case 'F':
            forward_only = true;
            break;

        case 'X':
            doubling_steps = atoi(optarg);
            break;

        case 'Z':
            size_limit = atoi(optarg);
            break;

        case 'T':
            store_threads = true;
            break;

        case 'o':
            discard_overlaps = true;
            break;

        case 'N':
            store_node_alignments = true;
            break;

        case 'h':
        case '?':
            help_index(argv);
            exit(1);
            break;

        default:
            abort ();
        }
    }

    if (edge_max == 0) edge_max = kmer_size + 1;

    vector<string> file_names;
    while (optind < argc) {
        string file_name = get_input_file_name(optind, argc, argv);
        file_names.push_back(file_name);
    }
    
    if (file_names.size() <= 0 && dbg_names.empty()){
        //cerr << "No graph provided for indexing. Please provide a .vg file or GCSA2-format deBruijn graph to index." << endl;
        //return 1;
    }

    if (kmer_size == 0 && !gcsa_name.empty() && dbg_names.empty()) {
        // gcsa doesn't do anything if we tell it a kmer size of 0.
        cerr << "error:[vg index] kmer size for GCSA2 index must be >0" << endl;
        return 1;
    }

    if (kmer_size < 0) {
        cerr << "error:[vg index] kmer size cannot be negative" << endl;
        return 1;
    }

    if (kmer_stride <= 0) {
        // kmer strides of 0 (or negative) are silly.
        cerr << "error:[vg index] kmer stride must be positive and nonzero" << endl;
        return 1;
    }

    if (!xg_name.empty()) {
        // We need to build an xg index

        // We'll fill this with the opened VCF file if we need one.
        vcflib::VariantCallFile variant_file;

        if (!vcf_name.empty()) {
            // There's a VCF we should load haplotype info from

            variant_file.open(vcf_name);
            if (!variant_file.is_open()) {
                cerr << "error:[vg index] could not open " << vcf_name << endl;
                return 1;
            } else if (show_progress) {
                cerr << "Opened variant file " << vcf_name << endl;
            }

        }

        // We want to siphon off the "_alt_<variant>_<number>" paths from "vg
        // construct -a" and not index them, and use them for creating haplotype
        // threads.
        // TODO: a better way to store path metadata
        map<string, Path> alt_paths;
        // This is matched against the entire string.
        regex is_alt("_alt_.+_[0-9]+");

        if (file_names.empty()) {
            // VGset or something segfaults when we feed it no graphs.
            cerr << "error:[vg index] at least one graph is required to build an xg index" << endl;
            return 1;
        }

        // store the graphs
        VGset graphs(file_names);
        // Turn into an XG index, except for the alt paths which we pull out and load into RAM instead.
        xg::XG index;
        graphs.to_xg(index, store_threads, is_alt, alt_paths);

        if (show_progress) {
            cerr << "Built base XG index" << endl;
        }

        // We're going to collect all the phase threads as XG threads (which
        // aren't huge like Protobuf Paths), and then insert them all into xg in
        // a batch, for speed. This will take a lot of memory (although not as
        // much as a real vg::Paths index or vector<Path> would)
        vector<xg::XG::thread_t> all_phase_threads;

        if (variant_file.is_open()) {
            // Now go through and add the varaints.

            // How many phases are there?
            size_t num_samples = variant_file.sampleNames.size();
            // And how many phases?
            size_t num_phases = num_samples * 2;

            for (size_t path_rank = 1; path_rank <= index.max_path_rank(); path_rank++) {
                // Find all the reference paths and loop over them. We'll just
                // assume paths that don't start with "_" might appear in the
                // VCF. We need to use the xg path functions, since we didn't
                // load up the whole vg graph.

                // What path is this?
                string path_name = index.path_name(path_rank);
                
                // Convert to VCF space if applicable
                string vcf_contig_name = path_to_vcf.count(path_name) ? path_to_vcf[path_name] : path_name;
                
                if (show_progress) {
                    cerr << "Processing path " << path_name << " as VCF contig " << vcf_contig_name << endl;
                }

                // We already know it's not a variant's alt, since those were
                // removed, so it might be a primary contig.

                // How many bases is it?
                size_t path_length = index.path_length(path_name);
                
                // We're going to extract it and index it, so we don't keep
                // making queries against it for every sample.
                PathIndex path_index(index.path(path_name));

                // Allocate some threads to store phase threads
                vector<xg::XG::thread_t> active_phase_threads{num_phases};
                // We need to remember how many paths of a particular phase have
                // already been generated.
                vector<int> saved_phase_paths(num_phases, 0);

                // What's the first reference position after the last variant?
                vector<size_t> nonvariant_starts(num_phases, 0);

                // Completed ones just get dumped into the index
                auto finish_phase = [&](size_t phase_number) {
                    // We have finished a phase (because an unphased variant
                    // came up or we ran out of variants); dump it into the
                    // index under a name and make a new Path for that phase.

                    // Find where this path is in our vector
                    xg::XG::thread_t& to_save = active_phase_threads[phase_number];

                    if (to_save.size() > 0) {
                        // Only actually do anything if we put in some mappings.

                        // Count this thread from this phase as being saved.
                        saved_phase_paths[phase_number]++;

                        // We don't tie threads from a pahse together in the
                        // index yet.

                        // Copy the thread over to our batch that we GPBWT all
                        // at once, exploiting the fact that VCF-derived graphs
                        // are DAGs.
                        all_phase_threads.push_back(to_save);

                        // Clear it out for re-use
                        to_save.clear();
                    }
                };

                // We need a way to convert Mappings to ThreadMappings
                // TODO: add a converting constructor?
                auto mapping_to_thread_mapping = [](const Mapping& mapping) {
                    xg::XG::ThreadMapping thread_mapping;
                    thread_mapping.node_id = mapping.position().node_id();
                    thread_mapping.is_reverse = mapping.position().is_reverse();
                    return thread_mapping;
                };
                // And NodeSides to thread mappings
                auto node_side_to_thread_mapping = [](const NodeSide& side) {
                    xg::XG::ThreadMapping thread_mapping;
                    thread_mapping.node_id = side.node;
                    thread_mapping.is_reverse = side.is_end;
                    return thread_mapping;
                };

                // We need a way to dump mappings into pahse threads. The
                // mapping edits and rank and offset info will be ignored; the
                // Mapping just represents an oriented node traversal.
                auto append_mapping = [&](size_t phase_number, const xg::XG::ThreadMapping& mapping) {
                    // Find the path to add to
                    xg::XG::thread_t& to_extend = active_phase_threads[phase_number];

                    // See if the edge we need to follow exists
                    if (to_extend.size() > 0) {
                        // If there's a previous mapping, go find it
                        const xg::XG::ThreadMapping& previous = to_extend[to_extend.size() - 1];

                        // Break out the IDs and flags we need to check for the edge
                        const int64_t& last_node = previous.node_id;
                        const bool& last_from_start = previous.is_reverse;

                        const int64_t& new_node = mapping.node_id;
                        const bool& new_to_end = mapping.is_reverse;

                        if (!index.has_edge(xg::make_edge(last_node, last_from_start, new_node, new_to_end))) {
                            // We can't have a thread take this edge (or an
                            // equivalent). Split and emit the current mappings
                            // and start a new path.
#ifdef debug
                            cerr << "warning:[vg index] phase " << phase_number << " wants edge "
                                << last_node << (last_from_start ? "L" : "R") << " - "
                                << new_node << (new_to_end ? "R" : "L")
                                << " which does not exist. Splitting!" << endl;
#endif
                            finish_phase(phase_number);
                        }
                    }

                    // Add the ThreadMapping
                    active_phase_threads[phase_number].push_back(mapping);
                };

                // We need an easy way to append any reference mappings from the
                // last variant up until a certain position (which may be past
                // the end of the entire reference path).
                auto append_reference_mappings_until = [&](size_t phase_number, size_t end) {
                    // We need to look and add in the mappings to the
                    // intervening reference nodes from the last variant, if
                    // any. For which we need access to the last variant's past-
                    // the-end reference position.
                    size_t ref_pos = nonvariant_starts[phase_number];
                    
                    // Get an iterator to the next node visit to add
                    PathIndex::iterator next_to_add = path_index.find_position(ref_pos);
                    
                    while(ref_pos < end && next_to_add != path_index.end()) {
                        // While there is intervening reference
                        // sequence, add it to our phase.

                        // What node side is the node that covers here?
                        NodeSide ref_side = next_to_add->second;

                        // Stick it in the phase path
                        append_mapping(phase_number, node_side_to_thread_mapping(ref_side));

                        // Advance to what's after that mapping, pulling node
                        // length from the path index
                        ref_pos += path_index.node_length(next_to_add);
                        
                        // Budge the iterator over so we don't need to do
                        // another tree query.
                        next_to_add++;
                    }
                    nonvariant_starts[phase_number] = ref_pos;
                };

                // We also have another function to handle each variant as it comes in.
                auto handle_variant = [&](vcflib::Variant& variant) {
                    // So we have a variant

                    // Grab its id, or make one by hashing stuff if it doesn't
                    // have an ID.
                    string var_name = make_variant_id(variant);
    
                    // We have alt paths like _alt_<var_name>_0 ...
                    // _alt_<var_name>_n. Up to one of them may be missing, in
                    // which case it represents a 0-length path that's just the
                    // edge from the node before the variable part of the
                    // variant to the node after.
                    
                    // If we take the ref allele when the ref path is missing,
                    // we don't care! We'll make mappings through here when we
                    // hit the next nonreference variant or the end of the
                    // contig and add the reference matches.
                    
                    // If we take an allele that's present, we go up through the
                    // end of the ref node that's before it, and then visit the
                    // allele.
                    
                    // If we take an alt allele when its path is missing, we go
                    // up to the end of the ref node before it, and then mark us
                    // as complete through there plus the length of the nodes
                    // along the ref path for the variant (which must be
                    // nonempty).

                    for (int sample_number = 0; sample_number < num_samples; sample_number++) {
                        // For each sample

                        // What sample is it?
                        string& sample_name = variant_file.sampleNames[sample_number];

                        // Parse it out and see if it's phased.
                        string genotype = variant.getGenotype(sample_name);

                        // Find the phasing bar
                        auto bar_pos = genotype.find('|');

                        if (bar_pos == string::npos || bar_pos == 0 || bar_pos + 1 >= genotype.size()) {
                            // If it isn't phased, or we otherwise don't like
                            // it, we need to break phasing paths.
                            for (int phase_offset = 0; phase_offset < 2; phase_offset++) {
                                // For each of the two phases for the sample
                                
                                // Remember where the end of the last variant was
                                auto cursor = nonvariant_starts[sample_number * 2 + phase_offset];
                                
                                // Make the phase thread reference up to the
                                // start of this variant. Doesn't have to be
                                // into the variable region.
                                append_reference_mappings_until(sample_number * 2 + phase_offset, variant.position);
                            
                                // Finish the phase thread and start a new one
                                finish_phase(sample_number * 2 + phase_offset);
                                
                                // Walk the cursor back so we repeat the
                                // reference segment, which we need to do in
                                // order to properly handle zero-length alleles
                                // at the ends of phase blocks.
                                nonvariant_starts[sample_number * 2 + phase_offset] = cursor;
                                
                                // TODO: we still can't handle deletions
                                // adjacent to SNPs where phasing gets lost. We
                                // have to have intervening reference bases. But
                                // that's a defect of the data model.
                            }
                        }
                        
                        // If it is phased, parse out the two alleles and handle
                        // each separately.
                        vector<string> alt_indices({genotype.substr(0, bar_pos),
                                genotype.substr(bar_pos + 1)});

                        for (int phase_offset = 0; phase_offset < 2; phase_offset++) {
                            // Handle each phase and its alt
                            string& alt_string = alt_indices[phase_offset];
                            
                            if (alt_string == ".") {
                                // This is a missing data call. Skip it. TODO:
                                // that means we'll just treat it like a
                                // reference call, when really we should break
                                // the haplotype here and not touch either alt.
                                // But if the reference path is empty, and we
                                // don't have a handy alt, we don't necessarily
                                // know where the site actually *is*, so we
                                // can't break phasing. What we should really do
                                // is iterate alt numbers until we find one, but
                                // that's going to be slow.
                                continue;
                            }
                            
                            // Otherwise it must be a proper number reference.
                            // Parse it.
                            int alt_index = stoi(alt_string);

                            if (alt_index != 0) {
                                // If this sample doesn't take the reference
                                // path at this variant, we need to actually go
                                // through it and not just call
                                // append_reference_mappings_until
                            
                                // We need to fill this in with the first
                                // reference position covered by the ref allele
                                // of this site, as actually represented in the
                                // path for the ref alt (i.e. after clipping
                                // fixed bases off the start and end in the
                                // VCF). This is the base after the insertion
                                // for pure insertions.
                                size_t first_ref_base = 0;
                                
                                // We need to look for the ref path for this variant
                                string ref_path_name = "_alt_" + var_name + "_0";
                                auto ref_path_iter = alt_paths.find(ref_path_name);
                                
                                // We also need to look for the path for this alt of this
                                // variant. 
                                string alt_path_name = "_alt_" + var_name + "_" + to_string(alt_index);
                                auto alt_path_iter = alt_paths.find(alt_path_name);
                                
                                
                                if (ref_path_iter != alt_paths.end() && ref_path_iter->second.mapping_size() != 0) {
                                    // We have the ref path so we can just look at its first node
                                    auto first_ref_node = ref_path_iter->second.mapping(0).position().node_id();
                                    
                                    // Find the first place it starts in the ref path
                                    first_ref_base = path_index.by_id.at(first_ref_node).first;
                                } else if (alt_path_iter != alt_paths.end() && alt_path_iter->second.mapping_size() != 0)  {
                                    // We have an alt path, so we can look at
                                    // the ref node before it and go one after
                                    // its end
                                    
                                    // Find the first node in the alt
                                    auto first_alt_id = alt_path_iter->second.mapping(0).position().node_id();
                                    bool first_alt_orientation = alt_path_iter->second.mapping(0).position().is_reverse();
                                    
                                    // Get all the edges coming in to it
                                    auto left_edges = (first_alt_orientation ? index.edges_on_end(first_alt_id) :
                                        index.edges_on_start(first_alt_id));
                                        
                                    // We need to fill in the ref to past the
                                    // end of the latest reference node that can
                                    // come before this alt.
                                    first_ref_base = 0;
                                    for (auto& edge : left_edges) {
                                        // For every edge, see what other node it attaches to
                                        auto other_id = (edge.from() == first_alt_id ? edge.to() : edge.from());
                                        if (other_id == first_alt_id) {
                                            // Skip self loops
                                            continue;
                                        }

                                        if (path_index.by_id.count(other_id) == 0) {
                                            // Skip nodes that aren't in the reference path
                                            continue;
                                        }

                                        // Look up where the node starts in the reference
                                        auto start = path_index.by_id.at(other_id).first;
                                        // There plus the length of the node will be the first ref base in our site
                                        first_ref_base = max(first_ref_base, start + index.node_length(other_id));
                                        
                                        // TODO: handling of cases where the alt
                                        // connects to multiple reference nodes
                                        // in different orientations.
                                    }
                                } else {
                                    // We lack both the ref and the alt path.
                                    // This site must have been skipped during
                                    // construction.
                                    cerr << "warning:[vg index] Alt and ref paths for " << var_name 
                                        << " at " << variant.sequenceName << ":" << variant.position
                                        << " missing/empty! Was variant skipped during construction?" << endl;
                                    continue;
                                }
                                
                                // Now we know the first ref base in our ref
                                // allele. What's the past-the-end base after we
                                // go through our ref allele.
                                size_t last_ref_base = first_ref_base;
                                if (ref_path_iter != alt_paths.end()) {
                                    for (size_t i = 0; i < ref_path_iter->second.mapping_size(); i++) {
                                        // Scoot it along with the length of
                                        // every node on our reference allele
                                        // path.
                                        last_ref_base += index.node_length(
                                            ref_path_iter->second.mapping(i).position().node_id());
                                    }
                                }
                            
                                if ((nonvariant_starts[sample_number * 2 + phase_offset] <= first_ref_base) ||
                                    !discard_overlaps) {
                                    
                                    // We need reference mappings from the last
                                    // variant up until the first actually
                                    // variable ref base in this site
                                    append_reference_mappings_until(sample_number * 2 + phase_offset, first_ref_base);

                                    for (size_t i = 0; (alt_path_iter != alt_paths.end() &&
                                        i < alt_path_iter->second.mapping_size()); i++) {
                                        // Then blit mappings from the alt over to the phase thread
                                        append_mapping(sample_number * 2 + phase_offset,
                                            mapping_to_thread_mapping(alt_path_iter->second.mapping(i)));
                                    }

                                    // Say we've accounted for the reference on
                                    // this path through the end of the variable
                                    // region, which we have.
                                    nonvariant_starts[sample_number * 2 + phase_offset] = last_ref_base;
                                }
                            }
                        }

                        // Now we have processed both phasings for this sample.
                    }
                };

                // Look for variants only on this path
                variant_file.setRegion(vcf_contig_name);

                // Set up progress bar
                ProgressBar* progress = nullptr;
                // Message needs to last as long as the bar itself.
                string progress_message = "loading variants for " + vcf_contig_name;
                if (show_progress) {
                    progress = new ProgressBar(path_length, progress_message.c_str());
                    progress->Progressed(0);
                }

                // Allocate a place to store actual variants
                vcflib::Variant var(variant_file);

                // How many variants have we done?
                size_t variants_processed = 0;
                while (variant_file.is_open() && variant_file.getNextVariant(var) && var.sequenceName == vcf_contig_name) {
                    // this ... maybe we should remove it as for when we have calls against N
                    bool isDNA = allATGC(var.ref);
                    for (vector<string>::iterator a = var.alt.begin(); a != var.alt.end(); ++a) {
                        if (!allATGC(*a)) isDNA = false;
                    }
                    // only work with DNA sequences
                    if (!isDNA) {
                        continue;
                    }
                    
                    var.position -= 1; // convert to 0-based
                    
                    // Handle the variant
                    handle_variant(var);


                    if (variants_processed++ % 1000 == 0 && progress != nullptr) {
                        // Say we made progress
                        progress->Progressed(var.position);
                    }
                }

                if (variants_processed > 0) {
                    // There were actually some variants on this path. We only
                    // want to actually have samples traverse the path if there
                    // were variants on it.

                    // Now finish up all the threads
                    for (size_t i = 0; i < num_phases; i++) {
                        // Each thread runs out until the end of the reference path
                        append_reference_mappings_until(i, path_length);

                        // And then we save all the threads
                        finish_phase(i);
                    }
                }

                if (progress != nullptr) {
                    // Throw out our progress bar
                    delete progress;
                    cerr << endl;
                    if (show_progress) {
                        cerr << "Processed " << variants_processed << " variants" << endl;
                    }
                }

            }

            if (show_progress) {
                cerr << "Inserting all phase threads into DAG..." << endl;
            }

            // Now insert all the threads in a batch into the known-DAG VCF-
            // derived graph.
            index.insert_threads_into_dag(all_phase_threads);
            all_phase_threads.clear();

        }

        if (show_progress) {
            cerr << "Saving index to disk..." << endl;
        }

        // save the xg version to the file name we've been given
        ofstream db_out(xg_name);
        index.serialize(db_out);
        db_out.close();
    }

    if (!gcsa_name.empty()) {
        // We need to make a gcsa index.

        // Configure GCSA2 verbosity so it doesn't spit out loads of extra info
        if (!show_progress) gcsa::Verbosity::set(gcsa::Verbosity::SILENT);

        // Load up the graphs
        vector<string> tmpfiles;
        if (dbg_names.empty()) {
            VGset graphs(file_names);
            graphs.show_progress = show_progress;
            // Go get the kmers of the correct size
            tmpfiles = graphs.write_gcsa_kmers_binary(kmer_size, path_only, forward_only);
        } else {
            tmpfiles = dbg_names;
        }
        // Make the index with the kmers
        gcsa::InputGraph input_graph(tmpfiles, true);
        gcsa::ConstructionParameters params;
        params.setSteps(doubling_steps);
        params.setLimit(size_limit);

        // build the GCSA index
        gcsa::GCSA* gcsa_index = new gcsa::GCSA(input_graph, params);

        // build the LCP array
        string lcp_name = gcsa_name + ".lcp";
        gcsa::LCPArray* lcp_array = new gcsa::LCPArray(input_graph, params);

        if (verify_index) {
            //cerr << "verifying index" << endl;
            if (!gcsa::verifyIndex(*gcsa_index, lcp_array, input_graph)) {
                cerr << "[vg::main]: GCSA2 index verification failed" << endl;
            }
        }

        // clean up input graph temp files
        if (dbg_names.empty()) {
            for (auto& tfn : tmpfiles) {
                remove(tfn.c_str());
            }
        }

        // Save the GCSA2 index
        sdsl::store_to_file(*gcsa_index, gcsa_name);
        delete gcsa_index;

        // Save the LCP array
        sdsl::store_to_file(*lcp_array, lcp_name);
        delete lcp_array;

    }

    if (!rocksdb_name.empty()) {

        Index index;

        if (compact) {
            index.open_for_write(rocksdb_name);
            index.compact();
            index.flush();
            index.close();
        }

        // todo, switch to xg for graph storage
        // index should write and load index/xg or such
        // then a handful of functions used in main.cpp and mapper.cpp need to be rewritten to use the xg index
        if (store_graph && file_names.size() > 0) {
            index.open_for_write(rocksdb_name);
            VGset graphs(file_names);
            graphs.show_progress = show_progress;
            graphs.store_in_index(index);
            //index.flush();
            //index.close();
            // reopen to index paths
            // this requires the index to be queryable
            //index.open_for_write(db_name);
            graphs.store_paths_in_index(index);
            index.compact();
            index.flush();
            index.close();
        }

        if (store_node_alignments && file_names.size() > 0) {
            index.open_for_bulk_load(rocksdb_name);
            int64_t aln_idx = 0;
            function<void(Alignment&)> lambda = [&index,&aln_idx](Alignment& aln) {
                index.cross_alignment(aln_idx++, aln);
            };
            for (auto& file_name : file_names) {
                get_input_file(file_name, [&](istream& in) {
                    stream::for_each_parallel(in, lambda);
                });
            }
            index.flush();
            index.close();
        }

        if (store_alignments && file_names.size() > 0) {
            index.open_for_bulk_load(rocksdb_name);
            function<void(Alignment&)> lambda = [&index](Alignment& aln) {
                index.put_alignment(aln);
            };
            for (auto& file_name : file_names) {
                get_input_file(file_name, [&](istream& in) {
                    stream::for_each_parallel(in, lambda);
                });
            }
            index.flush();
            index.close();
        }

        if (dump_alignments) {
            vector<Alignment> output_buf;
            index.open_read_only(rocksdb_name);
            auto lambda = [&output_buf](const Alignment& aln) {
                output_buf.push_back(aln);
                stream::write_buffered(cout, output_buf, 100);
            };
            index.for_each_alignment(lambda);
            stream::write_buffered(cout, output_buf, 0);
            index.close();
        }

        if (store_mappings && file_names.size() > 0) {
            index.open_for_bulk_load(rocksdb_name);
            function<void(Alignment&)> lambda = [&index](Alignment& aln) {
                const Path& path = aln.path();
                for (int i = 0; i < path.mapping_size(); ++i) {
                    index.put_mapping(path.mapping(i));
                }
            };
            for (auto& file_name : file_names) {
                get_input_file(file_name, [&](istream& in) {
                    stream::for_each_parallel(in, lambda);
                });
            }
            index.flush();
            index.close();
        }

        if (kmer_size != 0 && file_names.size() > 0) {
            index.open_for_bulk_load(rocksdb_name);
            VGset graphs(file_names);
            graphs.show_progress = show_progress;
            graphs.index_kmers(index, kmer_size, path_only, edge_max, kmer_stride, allow_negs);
            index.flush();
            index.close();
            // forces compaction
            index.open_for_write(rocksdb_name);
            index.flush();
            index.compact();
            index.close();
        }

        if (prune_kb >= 0) {
            if (show_progress) {
                cerr << "pruning kmers > " << prune_kb << " on disk from " << rocksdb_name << endl;
            }
            index.open_for_write(rocksdb_name);
            index.prune_kmers(prune_kb);
            index.compact();
            index.close();
        }

        if (set_kmer_size) {
            assert(kmer_size != 0);
            index.open_for_write(rocksdb_name);
            index.remember_kmer_size(kmer_size);
            index.close();
        }

        if (dump_index) {
            index.open_read_only(rocksdb_name);
            index.dump(cout);
            index.close();
        }

        if (describe_index) {
            index.open_read_only(rocksdb_name);
            set<int> kmer_sizes = index.stored_kmer_sizes();
            cout << "kmer sizes: ";
            for (auto kmer_size : kmer_sizes) {
                cout << kmer_size << " ";
            }
            cout << endl;
            index.close();
        }

        if (path_layout) {
            index.open_read_only(rocksdb_name);
            //index.path_layout();
            map<string, int64_t> path_by_id = index.paths_by_id();
            map<string, pair<pair<int64_t, bool>, pair<int64_t, bool>>> layout;
            map<string, int64_t> length;
            index.path_layout(layout, length);
            for (auto& p : layout) {
                // Negate IDs for backward nodes
                cout << p.first << " " << p.second.first.first * (p.second.first.second ? -1 : 1) << " "
                    << p.second.second.first * (p.second.second.second ? -1 : 1) << " " << length[p.first] << endl;
            }
            index.close();
        }
    }

    return 0;

}

// Register subcommand
static Subcommand vg_construct("index", "index graphs or alignments for random access or mapping", main_index);
