// construct.cpp: define the "vg construct" subcommand.

#include <omp.h>
#include <unistd.h>
#include <getopt.h>
#include <memory>

#include "subcommand.hpp"

#include <vg/io/stream.hpp>
#include "../constructor.hpp"
#include "../msa_converter.hpp"
#include "../region.hpp"

#include <bdsg/hash_graph.hpp>

using namespace std;
using namespace vg;
using namespace vg::subcommand;

void help_construct(char** argv) {
    cerr << "usage: " << argv[0] << " construct [options] >new.vg" << endl
         << "options:" << endl
         << "construct from a reference and variant calls:" << endl
         << "    -r, --reference FILE   input FASTA reference (may repeat)" << endl
         << "    -v, --vcf FILE         input VCF (may repeat)" << endl
         << "    -n, --rename V=F       match contig V in the VCFs to contig F in the FASTAs (may repeat)" << endl
         << "    -a, --alt-paths        save paths for alts of variants by SHA1 hash" << endl
         << "    -A, --alt-paths-plain  save paths for alts of variants by variant ID if possible, otherwise SHA1" << endl
         << "                           (IDs must be unique across all input VCFs)" << endl
         << "    -R, --region REGION    specify a VCF contig name or 1-based inclusive region (may repeat, if on different contigs)" << endl
         << "    -C, --region-is-chrom  don't attempt to parse the regions (use when the reference" << endl
         << "                           sequence name could be inadvertently parsed as a region)" << endl
         << "    -z, --region-size N    variants per region to parallelize (default: 1024)" << endl
         << "    -t, --threads N        use N threads to construct graph (defaults to numCPUs)" << endl
         << "    -S, --handle-sv        include structural variants in construction of graph." << endl
         << "    -I, --insertions FILE  a FASTA file containing insertion sequences "<< endl
         << "                           (referred to in VCF) to add to graph." << endl
         << "    -f, --flat-alts N      don't chop up alternate alleles from input VCF" << endl
         << "    -l, --parse-max N      don't chop up alternate alleles from input VCF longer than N (default: 100)" << endl
         << "    -i, --no-trim-indels   don't remove the 1bp reference base from alt alleles of indels." << endl
         << "    -N, --in-memory        construct the entire graph in memory before outputting it." <<endl
         << "construct from a multiple sequence alignment:" << endl
         << "    -M, --msa FILE         input multiple sequence alignment" << endl
         << "    -F, --msa-format       format of the MSA file (options: fasta, clustal; default fasta)" << endl
         << "    -d, --drop-msa-paths   don't add paths for the MSA sequences into the graph" << endl
         << "shared construction options:" << endl
         << "    -m, --node-max N       limit the maximum allowable node sequence size (default: 32)" << endl
         << "                           nodes greater than this threshold will be divided" << endl
         << "                           Note: nodes larger than ~1024 bp can't be GCSA2-indexed" << endl
         << "    -p, --progress         show progress" << endl;

}

int main_construct(int argc, char** argv) {

    if (argc == 2) {
        help_construct(argv);
        return 1;
    }

    // Make a constructor to fill in
    Constructor constructor;

    // We also parse some arguments separately.
    vector<string> fasta_filenames;
    vector<string> vcf_filenames;
    vector<string> insertion_filenames;
    vector<string> regions;
    bool region_is_chrom = false;
    string msa_filename;
    int max_node_size = 32;
    bool keep_paths = true;
    bool construct_in_memory = false;
    string msa_format = "fasta";
    bool show_progress = false;

    int c;
    optind = 2; // force optind past command positional argument
    while (true) {
        static struct option long_options[] =
            {
                /* These options set a flag. */
                //{"verbose", no_argument,       &verbose_flag, 1},
                {"vcf", required_argument, 0, 'v'},
                {"reference", required_argument, 0, 'r'},
                {"msa", required_argument, 0, 'M'},
                {"msa-format", required_argument, 0, 'F'},
                {"drop-msa-paths", no_argument, 0, 'd'},
                {"rename", required_argument, 0, 'n'},
                {"alt-paths", no_argument, 0, 'a'},
                {"alt-paths-plain", no_argument, 0, 'A'},
                {"handle-sv", no_argument, 0, 'S'},
                {"insertions", required_argument, 0, 'I'},
                {"progress",  no_argument, 0, 'p'},
                {"region-size", required_argument, 0, 'z'},
                {"threads", required_argument, 0, 't'},
                {"region", required_argument, 0, 'R'},
                {"region-is-chrom", no_argument, 0, 'C'},
                {"node-max", required_argument, 0, 'm'},
                {"flat-alts", no_argument, 0, 'f'},
                {"parse-max", required_argument, 0, 'l'},
                {"no-trim-indels", no_argument, 0, 'i'},
                {"in-memory", no_argument, 0, 'N'},
                {0, 0, 0, 0}
            };

        int option_index = 0;
        c = getopt_long (argc, argv, "v:r:n:ph?z:t:R:m:aACfl:SI:M:dF:iN",
                         long_options, &option_index);

        /* Detect the end of the options. */
        if (c == -1)
            break;

        switch (c)
        {
        case 'v':
            vcf_filenames.push_back(optarg);
            break;

        case 'M':
            msa_filename = optarg;
            break;
            
        case 'F':
            msa_format = optarg;
            break;
            
        case 'd':
            keep_paths = false;
            break;

        case 'i':
            constructor.trim_indels = false;
            break;

        case 'N':
            construct_in_memory = true;
            break;

        case 'r':
            fasta_filenames.push_back(optarg);
            break;

        case 'S':
            constructor.do_svs = true;
            break;

        case 'I':
            insertion_filenames.push_back(optarg);
            break;

            
        case 'n':
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
                string fasta_contig = key_value.substr(found + 1);
                // Add the name mapping
                constructor.add_name_mapping(vcf_contig, fasta_contig);
            }
            break;

        case 'a':
            constructor.alt_paths = true;
            break;

        case 'A':
            constructor.alt_paths = true;
            constructor.sha1_variant_name = false;
            break;

        case 'p':
            show_progress = true;
            break;

        case 'z':
            constructor.vars_per_chunk = parse<int>(optarg);
            break;

        case 'R':
            regions.push_back(optarg);
            break;

        case 'C':
            region_is_chrom = true;
            break;

        case 't':
            omp_set_num_threads(parse<int>(optarg));
            break;

        case 'm':
            max_node_size = parse<int>(optarg);
            break;

        case 'f':
            constructor.flat = true;
            break;
            
        case 'l':
            constructor.max_parsed_variant_size = parse<size_t>(optarg);
            break;

        case 'h':
        case '?':
            /* getopt_long already printed an error message. */
            help_construct(argv);
            exit(1);
            break;

        default:
            throw runtime_error("Not implemented: " + to_string(c));
        }
    }
    
    if (max_node_size == 0) {
        // Make sure we can actually make nodes
        cerr << "error:[vg construct] max node size cannot be 0" << endl;
        exit(1);
    }
    
    if (!msa_filename.empty() && !fasta_filenames.empty()) {
        cerr << "error:[vg construct] cannot construct from a reference/VCF and an MSA simultaneously" << endl;
        exit(1);
    }
    
    if (!fasta_filenames.empty()) {
        // Actually use the Constructor.
        // TODO: If we aren't always going to use the Constructor, refactor the subcommand to not always create and configure it.

        // Copy shared parameters into the constructor
        constructor.max_node_size = max_node_size;
        constructor.show_progress = show_progress;

        unordered_set<string> used_region_contigs; 
        for (auto& region : regions) {
            // We want to limit to one or more region
            if (!region_is_chrom) {
                // We are allowed to parse the region.
                // Break out sequence name and region bounds
                string seq_name;
                int64_t start_pos = -1, stop_pos = -1;
                parse_region(region,
                             seq_name,
                             start_pos,
                             stop_pos);
                             
                if (used_region_contigs.count(seq_name)) {
                    cerr << "error:[vg construct] cannot construct multiple regions of " << seq_name << endl;
                    exit(1);
                }
                used_region_contigs.insert(seq_name);
                
                if (start_pos > 0 && stop_pos > 0) {
                    // These are 0-based, so if both are nonzero we got a real set of coordinates
                    if (constructor.show_progress) {
                        cerr << "Restricting to " << seq_name << " from " << start_pos << " to " << stop_pos << endl;
                    }
                    constructor.allowed_vcf_names.insert(seq_name);
                    // Make sure to correct the coordinates to 0-based exclusive-end, from 1-based inclusive-end
                    constructor.allowed_vcf_regions[seq_name] = make_pair(start_pos - 1, stop_pos);
                } else if (start_pos < 0 && stop_pos < 0) {
                    // We just got a name
                    cerr << "Restricting to " << seq_name << " from 1 to end" << endl;
                    constructor.allowed_vcf_names.insert(seq_name);
                } else {
                    // This doesn't make sense. Does it have like one coordinate?
                    cerr << "error:[vg construct] could not parse " << region << endl;
                    exit(1);
                }
            } else {
                // We have been told not to parse the region
                cerr << "Restricting to " << region << " from 1 to end" << endl;
                constructor.allowed_vcf_names.insert(region);
            }
        }
        
        
        if (fasta_filenames.empty()) {
            cerr << "error:[vg construct] a reference is required for graph construction" << endl;
            return 1;
        }
        if (insertion_filenames.size() > 1){
            cerr << "Error: only one insertion file may be provided." << endl;
            exit(1);
        }
        
        if (construct_in_memory) {
            // Build the whole thing into memory
            bdsg::HashGraph constructed;
            constructor.construct_graph(fasta_filenames, vcf_filenames, insertion_filenames, &constructed);
            constructed.serialize(cout);
        } else {
            // Make an emitter that serializes the actual Graph objects, with buffering.
            // But just serialize one graph at a time in each group.
            // Make sure to compress the output.
            vg::io::ProtobufEmitter<Graph> emitter(cout, true, 1);

            // We need a callback to handle pieces of graph as they are produced.
            auto callback = [&](Graph& big_chunk) {
                // Sort the nodes by ID so that the serialized chunks come out in sorted order
                // TODO: We still interleave chunks from different threads working on different contigs
                std::sort(big_chunk.mutable_node()->begin(), big_chunk.mutable_node()->end(), [](const Node& a, const Node& b) -> bool {
                    // Return true if a comes before b
                    return a.id() < b.id();
                });
            
                // We don't validate the chunk because its end node may be held
                // back for the next chunk, while edges and path mappings for it
                // still live in this chunk. Also, we no longer create a VG to
                // re-chunk the chunk (because we can now handle chunks up to about
                // 1 GB serialized), and the VG class has the validator.
                
                // One thread at a time can write to the emitter and the output stream
    #pragma omp critical (emitter)
                emitter.write_copy(big_chunk); 
            };
            
            // Construct the graph.
            constructor.construct_graph(fasta_filenames, vcf_filenames, insertion_filenames, callback);

            // The output will be flushed when the ProtobufEmitter we use in the callback goes away.
            // Don't add an extra EOF marker or anything.
            
            // NB: If you worry about "still reachable but possibly lost" warnings in valgrind,
            // this would free all the memory used by protobuf:
            //ShutdownProtobufLibrary();
        }
    }
    else if (!msa_filename.empty()) {
        
        ifstream msa_file(msa_filename);
        if (!msa_file) {
            cerr << "error:[vg construct] could not open MSA file " << msa_filename << endl;
            exit(1);
        }
        
        MSAConverter msa_converter;
        msa_converter.show_progress = show_progress;
        
        msa_converter.load_alignments(msa_file, msa_format);
        VG msa_graph = msa_converter.make_graph(keep_paths, max_node_size);
        
        msa_graph.serialize_to_ostream(cout);
    }
    else {
        cerr << "error:[vg construct] a reference or an MSA is required for construct" << endl;
        exit(1);
    }

    return 0;
}

// Register subcommand
static Subcommand vg_construct("construct", "graph construction", PIPELINE, 2, main_construct);

