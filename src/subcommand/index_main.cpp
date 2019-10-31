// index.cpp: define the "vg index" subcommand, which makes xg, GCSA2, GBWT, and RocksDB indexes

#include <omp.h>
#include <unistd.h>
#include <getopt.h>

#include <random>
#include <string>
#include <vector>

#include "subcommand.hpp"

#include "../vg.hpp"
#include "../index.hpp"
#include "xg.hpp"
#include <vg/io/stream.hpp>
#include <vg/io/vpkg.hpp>
#include "../stream_index.hpp"
#include "../vg_set.hpp"
#include "../utility.hpp"
#include "../region.hpp"
#include "../snarls.hpp"
#include "../min_distance.hpp"
#include "../source_sink_overlay.hpp"
#include "../gbwt_helper.hpp"

#include <gcsa/gcsa.h>
#include <gcsa/algorithms.h>
#include <gbwt/variants.h>

using namespace std;
using namespace vg;
using namespace vg::subcommand;

void help_index(char** argv) {
    cerr << "usage: " << argv[0] << " index [options] <graph1.vg> [graph2.vg ...]" << endl
         << "Creates an index on the specified graph or graphs. All graphs indexed must " << endl
         << "already be in a joint ID space." << endl
         << "general options:" << endl
         << "    -b, --temp-dir DIR     use DIR for temporary files" << endl
         << "    -t, --threads N        number of threads to use" << endl
         << "    -p, --progress         show progress" << endl
         << "xg options:" << endl
         << "    -x, --xg-name FILE     use this file to store a succinct, queryable version of the graph(s), or read for GCSA or distance indexing" << endl
         << "    -L, --xg-alts          include alt paths in xg" << endl
         << "gbwt options:" << endl
         << "    -v, --vcf-phasing FILE generate threads from the haplotypes in the VCF file FILE" << endl
         << "    -W, --ignore-missing   don't warn when variants in the VCF are missing from the graph; silently skip them" << endl
         << "    -e, --parse-only FILE  store the VCF parsing with prefix FILE without generating threads" << endl
         << "    -T, --store-threads    generate threads from the embedded paths" << endl
         << "    -M, --store-gam FILE   generate threads from the alignments in FILE (many allowed)" << endl
         << "    -G, --gbwt-name FILE   store the threads as GBWT in FILE" << endl
         << "    -H, --write-haps FILE  store the threads as sequences in FILE" << endl
         << "    -z, --actual-phasing   do not make unphased homozygous genotypes phased"<< endl
         << "    -P, --force-phasing    replace unphased genotypes with randomly phased ones" << endl
         << "    -o, --discard-overlaps skip overlapping alternate alleles if the overlap cannot be resolved" << endl
         << "    -B, --batch-size N     number of samples per batch (default 200)" << endl
         << "    -u, --buffer-size N    GBWT construction buffer size in millions of nodes (default 100)" << endl
         << "    -n, --id-interval N    store haplotype ids at one out of N positions (default 1024)" << endl
         << "    -R, --range X..Y       process samples X to Y (inclusive)" << endl
         << "    -r, --rename V=P       rename contig V in the VCFs to path P in the graph (may repeat)" << endl
         << "    -I, --region C:S-E     operate on only the given 1-based region of the given VCF contig (may repeat)" << endl
         << "    -E, --exclude SAMPLE   exclude any samples with the given name from haplotype indexing" << endl
         << "gcsa options:" << endl
         << "    -g, --gcsa-out FILE    output a GCSA2 index to the given file" << endl
         << "    -i, --dbg-in FILE      use kmers from FILE instead of input VG (may repeat)" << endl
         << "    -f, --mapping FILE     use this node mapping in GCSA2 construction" << endl
         << "    -k, --kmer-size N      index kmers of size N in the graph (default " << gcsa::Key::MAX_LENGTH << ")" << endl
         << "    -X, --doubling-steps N use this number of doubling steps for GCSA2 construction (default " << gcsa::ConstructionParameters::DOUBLING_STEPS << ")" << endl
         << "    -Z, --size-limit N     limit temporary disk space usage to N gigabytes (default " << gcsa::ConstructionParameters::SIZE_LIMIT << ")" << endl
         << "    -V, --verify-index     validate the GCSA2 index using the input kmers (important for testing)" << endl
         << "gam indexing options:" << endl
         << "    -l, --index-sorted-gam input is sorted .gam format alignments, store a GAI index of the sorted GAM in INPUT.gam.gai" << endl
         << "vg in-place indexing options:" << endl
         << "    --index-sorted-vg      input is ID-sorted .vg format graph chunks, store a VGI index of the sorted vg in INPUT.vg.vgi" << endl
         << "rocksdb options:" << endl
         << "    -d, --db-name  <X>     store the RocksDB index in <X>" << endl
         << "    -m, --store-mappings   input is .gam format, store the mappings in alignments by node" << endl
         << "    -a, --store-alignments input is .gam format, store the alignments by node" << endl
         << "    -A, --dump-alignments  graph contains alignments, output them in sorted order" << endl
         << "    -N, --node-alignments  input is (ideally, sorted) .gam format," << endl
         << "                           cross reference nodes by alignment traversals" << endl
         << "    -D, --dump             print the contents of the db to stdout" << endl
         << "    -C, --compact          compact the index into a single level (improves performance)" << endl
         << "snarl distance index options" << endl
         << "    -s  --snarl-name FILE  load snarls from FILE" << endl
         << "    -j  --dist-name FILE   use this file to store a snarl-based distance index" << endl
         << "    -w  --max_dist N       cap beyond which the maximum distance is no longer accurate. If this is not included or is 0, don't build maximum distance index" << endl;
}

// Convert Path to a GBWT path.
gbwt::vector_type path_to_gbwt(const Path& path) {
    gbwt::vector_type result(path.mapping_size());
    for (size_t i = 0; i < result.size(); i++) {
        result[i] = mapping_to_gbwt(path.mapping(i));
    }
    return result;
}

// Find all predecessor nodes of the path, ignoring self-loops.
gbwt::vector_type predecessors(const xg::XG& xg_index, const Path& path) {
    gbwt::vector_type result;
    if (path.mapping_size() == 0) {
        return result;
    }

    vg::id_t first_node = path.mapping(0).position().node_id();
    bool is_reverse = path.mapping(0).position().is_reverse();
    
#ifdef debug
    cerr << "Look for predecessors of node " << first_node << " " << is_reverse << " which is first in alt path" << endl;
#endif

    xg_index.follow_edges(xg_index.get_handle(first_node), !is_reverse, [&] (handle_t next) {
            if (xg_index.get_id(next) != first_node) {
                result.push_back(gbwt::Node::encode(xg_index.get_id(next), xg_index.get_is_reverse(next)));
            }
        });

    return result;
}

std::vector<std::string> parseGenotypes(const std::string& vcf_line, size_t num_samples);

int main_index(int argc, char** argv) {

    if (argc == 2) {
        help_index(argv);
        return 1;
    }

    #define OPT_BUILD_VGI_INDEX 1000

    // Which indexes to build.
    bool build_xg = false, build_gbwt = false, write_threads = false, build_gcsa = false, build_rocksdb = false, build_dist = false;

    // Files we should read.
    string vcf_name, mapping_name;
    vector<string> dbg_names;

    // Files we should write.
    string xg_name, gbwt_name, parse_name, threads_name, gcsa_name, rocksdb_name, dist_name, snarl_name;

    // General
    bool show_progress = false;

    // GBWT
    bool warn_on_missing_variants = true;
    size_t found_missing_variants = 0; // Track the number of variants in the phasing VCF that aren't found in the graph
    size_t max_missing_variant_warnings = 10; // Only report up to this many of them
    bool index_haplotypes = false, index_paths = false, index_gam = false;
    bool parse_only = false;
    vector<string> gam_file_names;
    bool phase_homozygous = true, force_phasing = false, discard_overlaps = false;
    size_t samples_in_batch = 200;
    size_t gbwt_buffer_size = gbwt::DynamicGBWT::INSERT_BATCH_SIZE / gbwt::MILLION; // Millions of nodes.
    size_t id_interval = gbwt::DynamicGBWT::SAMPLE_INTERVAL;
    std::pair<size_t, size_t> sample_range(0, ~(size_t)0); // The semiopen range of samples to process.
    map<string, string> path_to_vcf; // Path name conversion from --rename.
    map<string, pair<size_t, size_t>> regions; // Region restrictions for contigs, in VCF name space, as 0-based exclusive-end ranges.
    unordered_set<string> excluded_samples; // Excluded sample names from --exclude.

    // GCSA
    gcsa::size_type kmer_size = gcsa::Key::MAX_LENGTH;
    gcsa::ConstructionParameters params;
    bool verify_gcsa = false;
    
    // Gam index (GAI)
    bool build_gai_index = false;
    
    // VG in-place index (VGI)
    bool build_vgi_index = false;

    // RocksDB
    bool dump_index = false;
    bool store_alignments = false;
    bool store_node_alignments = false;
    bool store_mappings = false;
    bool dump_alignments = false;

    //Distance index
    int cap = -1;
    bool include_maximum = false;

    // Unused?
    bool compact = false;

    // Include alt paths in xg
    bool xg_alts = false;

    int c;
    optind = 2; // force optind past command positional argument
    while (true) {
        static struct option long_options[] =
        {
            // General
            {"temp-dir", required_argument, 0, 'b'},
            {"threads", required_argument, 0, 't'},
            {"progress",  no_argument, 0, 'p'},

            // XG
            {"xg-name", required_argument, 0, 'x'},
            {"thread-db", required_argument, 0, 'F'},
            {"xg-alts", no_argument, 0, 'L'},

            // GBWT
            {"vcf-phasing", required_argument, 0, 'v'},
            {"ignore-missing", no_argument, 0, 'W'},
            {"parse-only", required_argument, 0, 'e'},
            {"store-threads", no_argument, 0, 'T'},
            {"store-gam", required_argument, 0, 'M'},
            {"gbwt-name", required_argument, 0, 'G'},
            {"write-haps", required_argument, 0, 'H'},
            {"actual-phasing", no_argument, 0, 'z'},
            {"force-phasing", no_argument, 0, 'P'},
            {"discard-overlaps", no_argument, 0, 'o'},
            {"batch-size", required_argument, 0, 'B'},
            {"buffer-size", required_argument, 0, 'u'},
            {"id-interval", required_argument, 0, 'n'},
            {"range", required_argument, 0, 'R'},
            {"rename", required_argument, 0, 'r'},
            {"region", required_argument, 0, 'I'},
            {"exclude", required_argument, 0, 'E'},

            // GCSA
            {"gcsa-out", required_argument, 0, 'g'},
            {"dbg-in", required_argument, 0, 'i'},
            {"mapping", required_argument, 0, 'f'},
            {"kmer-size", required_argument, 0, 'k'},
            {"doubling-steps", required_argument, 0, 'X'},
            {"size-limit", required_argument, 0, 'Z'},
            {"verify-index", no_argument, 0, 'V'},
            
            // GAM index (GAI)
            {"index-sorted-gam", no_argument, 0, 'l'},
            
            // VG in-place index (VGI)
            {"index-sorted-vg", no_argument, 0, OPT_BUILD_VGI_INDEX},

            // RocksDB
            {"db-name", required_argument, 0, 'd'},
            {"store-mappings", no_argument, 0, 'm'},
            {"store-alignments", no_argument, 0, 'a'},
            {"dump-alignments", no_argument, 0, 'A'},
            {"node-alignments", no_argument, 0, 'N'},
            {"dump", no_argument, 0, 'D'},
            {"compact", no_argument, 0, 'C'},

            //Snarl distance index
            {"snarl-name", required_argument, 0, 's'},
            {"dist-name", required_argument, 0, 'j'},
            {"max-dist", required_argument, 0, 'w'},
            {0, 0, 0, 0}
        };

        int option_index = 0;
        c = getopt_long (argc, argv, "b:t:px:Lv:We:TM:G:H:zPoB:u:n:R:r:I:E:g:i:f:k:X:Z:Vld:maANDCs:j:w:h",
                long_options, &option_index);

        // Detect the end of the options.
        if (c == -1)
            break;

        switch (c)
        {
        // General
        case 'b':
            temp_file::set_dir(optarg);
            break;
        case 't':
            omp_set_num_threads(parse<int>(optarg));
            break;
        case 'p':
            show_progress = true;
            break;

        // XG
        case 'x':
            build_xg = true;
            xg_name = optarg;
            break;
        case 'L':
            xg_alts = true;
            break;

        // GBWT
        case 'v':
            index_haplotypes = true;
            build_xg = true;
            vcf_name = optarg;
            break;
        case 'W':
            warn_on_missing_variants = false;
            break;
        case 'e':
            parse_only = true;
            parse_name = optarg;
            break;
        case 'T':
            index_paths = true;
            build_xg = true;
            break;
        case 'M':
            index_gam = true;
            build_gbwt = true;
            gam_file_names.push_back(optarg);
            break;
        case 'G':
            build_gbwt = true;
            gbwt_name = optarg;
            break;
        case 'H':
            write_threads = true;
            threads_name = optarg;
            break;
        case 'z':
            phase_homozygous = false;
            break;
        case 'P':
            force_phasing = true;
            break;
        case 'o':
            discard_overlaps = true;
            break;
        case 'B':
            samples_in_batch = std::max(parse<size_t>(optarg), 1ul);
            break;
        case 'u':
            gbwt_buffer_size = std::max(parse<size_t>(optarg), 1ul);
            break;
        case 'n':
            id_interval = parse<size_t>(optarg);
            break;
        case 'R':
            {
                // Parse first..last
                string temp(optarg);
                size_t found = temp.find("..");
                if(found == string::npos || found == 0 || found + 2 == temp.size()) {
                    cerr << "error: [vg index] could not parse range " << temp << endl;
                    exit(1);
                }
                sample_range.first = parse<size_t>(temp.substr(0, found));
                sample_range.second = parse<size_t>(temp.substr(found + 2)) + 1;
            }
            break;
        case 'r':
            {
                // Parse the rename old=new
                string key_value(optarg);
                auto found = key_value.find('=');
                if (found == string::npos || found == 0 || found + 1 == key_value.size()) {
                    cerr << "error: [vg index] could not parse rename " << key_value << endl;
                    exit(1);
                }
                // Parse out the two parts
                string vcf_contig = key_value.substr(0, found);
                string graph_contig = key_value.substr(found + 1);
                // Add the name mapping
                path_to_vcf[graph_contig] = vcf_contig;
            }
            break;
        case 'I':
            {
                // We want to parse this region specifier
                string region(optarg);
                
                Region parsed;
                parse_region(region, parsed);
                if (parsed.start <= 0 || parsed.end <= 0) {
                    // We need both range bounds, and we can't accept 0 since input is 1-based.
                    cerr << "error: [vg index] could not parse 1-based region " << optarg << endl;
                }
                
                // Make sure to correct the coordinates to 0-based exclusive-end, from 1-based inclusive-end
                regions[parsed.seq] = make_pair((size_t) (parsed.start - 1), (size_t) parsed.end);
            }
            break;
        case 'E':
            excluded_samples.insert(optarg);
            break;

        // GCSA
        case 'g':
            build_gcsa = true;
            gcsa_name = optarg;
            break;
        case 'i':
            dbg_names.push_back(optarg);
            break;
        case 'f':
            mapping_name = optarg;
            break;
        case 'k':
            kmer_size = std::max(parse<size_t>(optarg), 1ul);
            break;
        case 'X':
            params.setSteps(parse<size_t>(optarg));
            break;
        case 'Z':
            params.setLimit(parse<size_t>(optarg));
            break;
        case 'V':
            verify_gcsa = true;
            break;
            
        // Gam index (GAI)
        case 'l':
            build_gai_index = true;
            break;
            
        // VGI index
        case OPT_BUILD_VGI_INDEX:
            build_vgi_index = true;
            break;

        // RocksDB
        case 'd':
            build_rocksdb = true;
            rocksdb_name = optarg;
            break;
        case 'm':
            store_mappings = true;
            break;
        case 'a':
            store_alignments = true;
            break;
        case 'A':
            dump_alignments = true;
            break;
        case 'N':
            store_node_alignments = true;
            break;
        case 'D':
            dump_index = true;
            break;

        case 'C':
            compact = true;
            break;

        //Snarl distance index
        case 's':
            snarl_name = optarg;
            break;
        case 'j':
            build_dist = true;
            dist_name = optarg;
            break;
        case 'w':
            build_dist = true;
            cap = parse<int>(optarg);
            include_maximum = true;
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

    vector<string> file_names;
    while (optind < argc) {
        string file_name = get_input_file_name(optind, argc, argv);
        file_names.push_back(file_name);
    }

    if (xg_name.empty() && gbwt_name.empty() && parse_name.empty() && threads_name.empty() &&
        gcsa_name.empty() && rocksdb_name.empty() && !build_gai_index && !build_vgi_index && dist_name.empty()) {
        cerr << "error: [vg index] index type not specified" << endl;
        return 1;
    }

    if ((build_gbwt || write_threads) && !(index_haplotypes || index_paths || index_gam)) {
        cerr << "error: [vg index] cannot build GBWT without threads" << endl;
        return 1;
    }

    if (parse_only && (!index_haplotypes || index_paths || index_gam)) {
        cerr << "error: [vg index] --parse-only works with --vcf-phasing only" << endl;
        return 1;
    }

    if ((index_haplotypes || index_paths || index_gam) && !(build_gbwt || parse_only || write_threads)) {
        cerr << "error: [vg index] no output format specified for the threads" << endl;
        return 1;
    }

    if (index_gam && (index_haplotypes || index_paths)) {
        cerr << "error: [vg index] GAM threads are incompatible with haplotype/path threads" << endl;
        return 1;
    }

    if (file_names.size() <= 0 && dbg_names.empty()){
        //cerr << "No graph provided for indexing. Please provide a .vg file or GCSA2-format deBruijn graph to index." << endl;
        //return 1;
    }
    
    if (file_names.size() != 1 && build_gai_index) {
        cerr << "error: [vg index] can only index exactly one sorted GAM file at a time" << endl;
        return 1;
    }
    
    if (file_names.size() != 1 && build_vgi_index) {
        cerr << "error: [vg index] can only index exactly one sorted VG file at a time" << endl;
        return 1;
    }
    
    if (file_names.size() > 1 && build_dist) {
        // Allow zero filenames for the index-from-xg mode
        cerr << "error: [vg index] can only create one distance index at a time" << endl;
        return 1;
    }
    
    if (build_gcsa && kmer_size > gcsa::Key::MAX_LENGTH) {
        cerr << "error: [vg index] GCSA2 cannot index with kmer size greater than " << gcsa::Key::MAX_LENGTH << endl;
        return 1;
    }
    
    if (build_xg && build_gcsa && file_names.empty()) {
        // Really we want to build a GCSA by *reading* and XG
        build_xg = false;
        // We'll continue in the build_gcsa section
    }
    if (build_dist && file_names.empty()) {
        //If we want to build the distance index from the xg
        build_xg = false;
    }

    // Build XG
    xg::XG* xg_index = new xg::XG();
    map<string, Path> alt_paths;
    if (build_xg) {
        if (file_names.empty()) {
            // VGset or something segfaults when we feed it no graphs.
            cerr << "error: [vg index] at least one graph is required to build an xg index" << endl;
            return 1;
        }
        VGset graphs(file_names);
        graphs.to_xg(*xg_index, xg_alts ? [](const string&) {return false;} : Paths::is_alt, index_haplotypes ? &alt_paths : nullptr);
        if (index_haplotypes && xg_alts) {
            assert(alt_paths.empty());
            // they weren't filtered in the above: re-extract here
            xg_index->for_each_path_handle([&](path_handle_t path_handle) {
                    string path_name = xg_index->get_path_name(path_handle);
                    if (Paths::is_alt(path_name)) {
                        alt_paths[path_name] = path_from_path_handle(*xg_index, path_handle);
                    }
                });
        }
        if (show_progress) {
            cerr << "Built base XG index" << endl;
        }
    }

#ifdef debug
    cerr << "Alt paths:" << endl;
    for (auto& kv : alt_paths) {
        cerr << kv.first << ": " << kv.second.mapping_size() << " entries" << endl;
    }
#endif

    // Generate threads
    if (index_haplotypes || index_paths || index_gam) {

        // Use the same temp directory as VG.
        gbwt::TempFile::setDirectory(temp_file::get_dir());

        // GBWT metadata.
        std::vector<std::string> sample_names, contig_names;
        size_t haplotype_count = 0;
        size_t true_sample_offset = 0; // Id of the first VCF sample.

        // Determine node id width.
        size_t id_width;
        if (!index_gam) {
            id_width = gbwt::bit_length(gbwt::Node::encode(xg_index->max_node_id(), true));
        } else { // indexing a GAM
            if (show_progress) {
                cerr << "Finding maximum node id in GAM..." << endl;
            }
            vg::id_t max_id = 0;
            size_t alignments_in_gam = 0;
            function<void(Alignment&)> lambda = [&](Alignment& aln) {
                gbwt::vector_type buffer;
                for (auto& m : aln.path().mapping()) {
                    max_id = max(m.position().node_id(), max_id);
                }
                alignments_in_gam++;
            };
            for (auto& file_name : gam_file_names) {
                get_input_file(file_name, [&](istream& in) {
                    vg::io::for_each(in, lambda);
                });
            }
            id_width = gbwt::bit_length(gbwt::Node::encode(max_id, true));
            sample_names.reserve(alignments_in_gam); // We store alignment names as sample names.
        }
        
        if (show_progress) {
            cerr << "Node id width: " << id_width << endl;
        }

        // Do we build GBWT?
        gbwt::GBWTBuilder* gbwt_builder = 0;
        if (build_gbwt) {
            if (show_progress) {
                cerr << "GBWT parameters: buffer size " << gbwt_buffer_size << ", id interval " << id_interval << endl;
            }
            gbwt::Verbosity::set(gbwt::Verbosity::SILENT);  // Make the construction thread silent.
            gbwt_builder = new gbwt::GBWTBuilder(id_width, gbwt_buffer_size * gbwt::MILLION, id_interval);
            gbwt_builder->index.addMetadata();
        }

        // Do we write threads?
        gbwt::text_buffer_type binary_file;
        if (write_threads) {
            if (show_progress) { cerr << "Writing the threads to " << threads_name << endl; }
            binary_file = gbwt::text_buffer_type(threads_name, std::ios::out, gbwt::MEGABYTE, id_width);
        }

        // Store a thread.
        auto store_thread = [&](const gbwt::vector_type& to_save) {
            if (build_gbwt) {
                gbwt_builder->insert(to_save, true); // Insert in both orientations.
            }
            if (write_threads) {
                for (auto node : to_save) { binary_file.push_back(node); }
                binary_file.push_back(gbwt::ENDMARKER);
            }
        };

        // Store the name of the thread in GBWT metadata.
        auto store_thread_name = [&](gbwt::size_type sample, gbwt::size_type contig, gbwt::size_type phase, gbwt::size_type count) {
            if (build_gbwt) {
                gbwt_builder->index.metadata.addPath({
                    static_cast<gbwt::PathName::path_name_type>(sample),
                    static_cast<gbwt::PathName::path_name_type>(contig),
                    static_cast<gbwt::PathName::path_name_type>(phase),
                    static_cast<gbwt::PathName::path_name_type>(count)
                });
            }
        };

        // Store contig names.
        if (index_paths || index_haplotypes) {
            xg_index->for_each_path_handle([&](path_handle_t path_handle) {
                    string path_name = xg_index->get_path_name(path_handle);
                    if (!Paths::is_alt(path_name)) {
                        contig_names.push_back(path_name);
                    }
                });
        }
        // Convert paths to threads
        if (index_paths) {
            if (show_progress) {
                cerr << "Converting paths to threads..." << endl;
            }
            int path_rank = 0;
            xg_index->for_each_path_handle([&](path_handle_t path_handle) {
                    ++path_rank;
                    if (xg_index->is_empty(path_handle) || Paths::is_alt(xg_index->get_path_name(path_handle))) {
                        return;
                    }
                    gbwt::vector_type buffer;
                    buffer.reserve(xg_index->get_step_count(path_handle));
                    for (handle_t handle : xg_index->scan_path(path_handle)) {
                        buffer.push_back(gbwt::Node::encode(xg_index->get_id(handle), xg_index->get_is_reverse(handle)));
                    }
                    store_thread(buffer);
                    store_thread_name(true_sample_offset, path_rank - 1, 0, 0);
                });

            // GBWT metadata: We assume that the XG index contains the reference paths.
            sample_names.emplace_back("ref");
            haplotype_count++;
            true_sample_offset++;
        }

        // Index GAM, using alignment names as sample names.
        if (index_gam) {
            if (show_progress) {
                cerr << "Converting GAM to threads..." << endl;
            }
            function<void(Alignment&)> lambda = [&](Alignment& aln) {
                gbwt::vector_type buffer;
                for (auto& m : aln.path().mapping()) {
                    buffer.push_back(mapping_to_gbwt(m));
                }
                store_thread(buffer);
                store_thread_name(sample_names.size(), 0, 0, 0);
                sample_names.emplace_back(aln.name());
                haplotype_count++;
                true_sample_offset++;
            };
            for (auto& file_name : gam_file_names) {
                get_input_file(file_name, [&](istream& in) {
                    vg::io::for_each(in, lambda);
                });
            }
        }

        // Generate haplotypes
        if (index_haplotypes) {
            size_t total_variants_processed = 0;
            vcflib::VariantCallFile variant_file;
            variant_file.parseSamples = false; // vcflib parsing is very slow if there are many samples.
            variant_file.open(vcf_name);
            if (!variant_file.is_open()) {
                cerr << "error: [vg index] could not open " << vcf_name << endl;
                return 1;
            } else if (show_progress) {
                cerr << "Opened variant file " << vcf_name << endl;
            }
            std::mt19937 rng(0xDEADBEEF);
            std::uniform_int_distribution<std::mt19937::result_type> random_bit(0, 1);

            // How many samples are there?
            size_t num_samples = variant_file.sampleNames.size();
            if (num_samples == 0) {
                cerr << "error: [vg index] the variant file does not contain phasings" << endl;
                return 1;
            }

            // Determine the samples we want to index.
            sample_range.second = std::min(sample_range.second, num_samples);
            sample_names.insert(sample_names.end(),
                                variant_file.sampleNames.begin() + sample_range.first,
                                variant_file.sampleNames.begin() + sample_range.second);
            haplotype_count += 2 * (sample_range.second - sample_range.first);  // Assuming a diploid genome
            if (show_progress) {
                cerr << "Haplotype generation parameters:" << endl;
                cerr << "- Samples " << sample_range.first << " to " << (sample_range.second - 1) << endl;
                cerr << "- Batch size " << samples_in_batch << endl;
                if (phase_homozygous) {
                    cerr << "- Phase homozygous genotypes" << endl;
                }
                if (force_phasing) {
                    cerr << "- Force phasing" << endl;
                }
                if (discard_overlaps) {
                    cerr << "- Discard overlaps" << endl;
                }
            }

            // Process each VCF contig corresponding to an XG path.
            vector<path_handle_t> path_handles;
            // 1st pass: scan for all non-alt paths (they are handled separately)
            xg_index->for_each_path_handle([&](path_handle_t path_handle) {
                    if (!alt_paths.count(xg_index->get_path_name(path_handle))) {
                        path_handles.push_back(path_handle);
                    }
                });
            size_t max_path_rank = path_handles.size();
            for (size_t path_rank = 1; path_rank <= max_path_rank; path_rank++) {
                string path_name = xg_index->get_path_name(path_handles[path_rank - 1]);
                if (Paths::is_alt(path_name)) {
                    continue;
                }
                string vcf_contig_name = path_to_vcf.count(path_name) ? path_to_vcf[path_name] : path_name;
                if (show_progress) {
                    cerr << "Processing path " << path_name << " as VCF contig " << vcf_contig_name << endl;
                }
                string parse_file = parse_name + '_' + vcf_contig_name;

                // Structures to parse the VCF file into.
                gbwt::VariantPaths variants(xg_index->get_step_count(xg_index->get_path_handle(path_name)));
                variants.setSampleNames(sample_names);
                variants.setContigName(path_name);
                std::vector<gbwt::PhasingInformation> phasings;

                // Add the reference to VariantPaths.
                for (handle_t handle : xg_index->scan_path(path_handles[path_rank - 1])) {
                    variants.appendToReference(gbwt::Node::encode(xg_index->get_id(handle), xg_index->get_is_reverse(handle)));
                }
                variants.indexReference();

                // Create a PhasingInformation for each batch.
                for (size_t batch_start = sample_range.first; batch_start < sample_range.second; batch_start += samples_in_batch) {
                    if (parse_only) {
                        // Use a permanent file.
                        phasings.emplace_back(parse_file, batch_start, std::min(samples_in_batch, sample_range.second - batch_start));
                        variants.addFile(phasings.back().name(), phasings.back().offset(), phasings.back().size());
                    } else {
                        // Use a temporary file.
                        phasings.emplace_back(batch_start, std::min(samples_in_batch, sample_range.second - batch_start));
                    }
                }

                // Set the VCF region or process the entire contig.
                if (regions.count(vcf_contig_name)) {
                    auto region = regions[vcf_contig_name];
                    if (show_progress) {
                        cerr << "- Setting region " << region.first << " to " << region.second << endl;
                    }
                    variant_file.setRegion(vcf_contig_name, region.first, region.second);
                } else {
                    variant_file.setRegion(vcf_contig_name);
                }

                // Parse the variants and the phasings.
                vcflib::Variant var(variant_file);
                size_t variants_processed = 0;
                std::vector<bool> was_diploid(sample_range.second, true); // Was the sample diploid at the previous site?
                while (variant_file.is_open() && variant_file.getNextVariant(var) && var.sequenceName == vcf_contig_name) {
                    // Skip variants with non-DNA sequence, as they are not included in the graph.
                    bool isDNA = allATGC(var.ref);
                    for (vector<string>::iterator a = var.alt.begin(); a != var.alt.end(); ++a) {
                        if (!allATGC(*a)) isDNA = false;
                    }
                    if (!isDNA) {
                        continue;
                    }

                    // Determine the reference nodes for the current variant and create a variant site.
                    // If the variant is not an insertion, there should be a path for the ref allele.
                    
                    std::string var_name = make_variant_id(var);
                    std::string ref_path_name = "_alt_" + var_name + "_0";
                    auto ref_path_iter = alt_paths.find(ref_path_name);
                    gbwt::vector_type ref_path;
                    size_t ref_pos = variants.invalid_position();
                    if (ref_path_iter != alt_paths.end() && ref_path_iter->second.mapping_size() != 0) {
                        ref_path = path_to_gbwt(ref_path_iter->second);
                        ref_pos = variants.firstOccurrence(ref_path.front());
                        if (ref_pos == variants.invalid_position()) {
                            cerr << "warning: [vg index] invalid ref path for " << var_name << " at "
                                 << var.sequenceName << ":" << var.position << endl;
                            continue;
                        }
                    } else {
                        // Try using alt paths instead.
                        bool found = false;
                        for (size_t alt_index = 1; alt_index < var.alleles.size(); alt_index++) {
                            std::string alt_path_name = "_alt_" + var_name + "_" + to_string(alt_index);
                            size_t candidate_pos = 0;
                            bool candidate_found = false;
                            auto alt_path_iter = alt_paths.find(alt_path_name);
                            if (alt_path_iter != alt_paths.end()) {
                                gbwt::vector_type pred_nodes = predecessors(*xg_index, alt_path_iter->second);
                                for (auto node : pred_nodes) {
                                    size_t pred_pos = variants.firstOccurrence(node);
                                    if (pred_pos != variants.invalid_position()) {
#ifdef debug
                                        cerr << "Found predecessor node " << gbwt::Node::id(node) << " " << gbwt::Node::is_reverse(node)
                                             << " occurring at valid pos " << pred_pos << endl;
#endif
                                        candidate_pos = std::max(candidate_pos, pred_pos + 1);
                                        candidate_found = true;
                                        found = true;
                                    }
                                }
                                // For each alternate allele, find the rightmost reference node among
                                // its predecessors. If multiple alleles have candidates for the
                                // reference position, choose the leftmost one.
                                if (candidate_found) {
                                    ref_pos = std::min(ref_pos, candidate_pos);
                                }
                            }
                        }
                        if (!found) {
                            // This variant from the VCF is just not in the graph

                            found_missing_variants++;

                            if (warn_on_missing_variants) {
                                if (found_missing_variants <= max_missing_variant_warnings) {
                                    // The user might not know it. Warn them in case they mixed up their VCFs.
                                    cerr << "warning: [vg index] alt and ref paths for " << var_name
                                         << " at " << var.sequenceName << ":" << var.position
                                         << " missing/empty! Was the variant skipped during construction?" << endl;
                                    if (found_missing_variants == max_missing_variant_warnings) {
                                        cerr << "warning: [vg index] suppressing further missing variant warnings" << endl;
                                    }
                                }
                            }

                            // Skip this variant and move on to the next as if it never appeared.
                            continue;
                        }
                    }
                    variants.addSite(ref_pos, ref_pos + ref_path.size());

                    // Add alternate alleles to the site.
                    for (size_t alt_index = 1; alt_index < var.alleles.size(); alt_index++) {
                        std::string alt_path_name = "_alt_" + var_name + "_" + to_string(alt_index);
                        auto alt_path_iter = alt_paths.find(alt_path_name);
                        if (alt_path_iter != alt_paths.end()) {
                            variants.addAllele(path_to_gbwt(alt_path_iter->second));
                        } else {
                            variants.addAllele(ref_path);
                        }
                    }

                    // Store the phasings in PhasingInformation structures.
                    std::vector<std::string> genotypes = parseGenotypes(var.originalLine, num_samples);
                    for (size_t batch = 0; batch < phasings.size(); batch++) {
                        std::vector<gbwt::Phasing> current_phasings;
                        for (size_t sample = phasings[batch].offset(); sample < phasings[batch].limit(); sample++) {
                            current_phasings.emplace_back(genotypes[sample], was_diploid[sample], phase_homozygous);
                            was_diploid[sample] = current_phasings.back().diploid;
                            if(force_phasing) {
                                current_phasings.back().forcePhased([&]() {
                                        return random_bit(rng);
                                    });
                            }
                        }
                        phasings[batch].append(current_phasings);
                    }
                    variants_processed++;
                } // End of variants.
                if (show_progress) {
                    cerr << "- Parsed " << variants_processed << " variants" << endl;
                    size_t phasing_bytes = 0;
                    for (size_t batch = 0; batch < phasings.size(); batch++) {
                        phasing_bytes += phasings[batch].bytes();
                    }
                    cerr << "- Phasing information: " << gbwt::inMegabytes(phasing_bytes) << " MB" << endl;
                }

                // Save memory:
                // - Delete the alt paths if we no longer need them.
                // - Delete the XG index if we no longer need it.
                // - Close the phasings files.
                if (path_rank == max_path_rank) {
                    alt_paths.clear();
                    if (xg_name.empty()) {
                        delete xg_index;
                        xg_index = nullptr;
                    }
                }
                for (size_t batch = 0; batch < phasings.size(); batch++) {
                    phasings[batch].close();
                }

                // Save the VCF parse or generate the haplotypes.
                if (parse_only) {
                    if (!sdsl::store_to_file(variants, parse_file)) {
                        cerr << "error: [vg index] cannot write parse file " << parse_file << endl;
                        return 1;
                    }
                } else {
                    for (size_t batch = 0; batch < phasings.size(); batch++) {
                        gbwt::generateHaplotypes(variants, phasings[batch],
                                                 [&](gbwt::size_type sample) -> bool {
                                                     return (excluded_samples.find(variant_file.sampleNames[sample]) == excluded_samples.end());
                                                 },
                                                 [&](const gbwt::Haplotype& haplotype) {
                                                     store_thread(haplotype.path);
                                                     store_thread_name(haplotype.sample + true_sample_offset - sample_range.first,
                                                                       path_rank - 1,
                                                                       haplotype.phase,
                                                                       haplotype.count);
                                                 },
                                                 [&](gbwt::size_type, gbwt::size_type) -> bool {
                                                     return discard_overlaps;
                                                 });
                        if (show_progress) {
                            cerr << "- Processed samples " << phasings[batch].offset() << " to " << (phasings[batch].offset() + phasings[batch].size() - 1) << endl;
                        }
                    }
                } // End of haplotype generation for the current contig.
            
                // Record the number of variants we saw on this contig
                total_variants_processed += variants_processed;
            
            } // End of contigs.
            
        if (warn_on_missing_variants && found_missing_variants > 0) {
            cerr << "warning: [vg index] Found " << found_missing_variants << "/" << total_variants_processed
                 << " variants in phasing VCF but not in graph! Do your graph and VCF match?" << endl;
        }
    } // End of haplotypes.
        
    // Write the threads to disk.
    alt_paths.clear();
    if (!parse_only) {
        if (build_gbwt) {
            gbwt_builder->finish();
            gbwt_builder->index.metadata.setSamples(sample_names);
            gbwt_builder->index.metadata.setHaplotypes(haplotype_count);
            gbwt_builder->index.metadata.setContigs(contig_names);
            if (show_progress) {
                cerr << "GBWT metadata: "; gbwt::operator<<(cerr, gbwt_builder->index.metadata); cerr << endl;
                cerr << "Saving GBWT to disk..." << endl;
            }
                
                // Save encapsulated in a VPKG
                vg::io::VPKG::save(gbwt_builder->index, gbwt_name);
                
                delete gbwt_builder; gbwt_builder = nullptr;
            }
            if (write_threads) {
                binary_file.close();
            }
        }
    } // End of thread indexing.

    // Save XG
    if (build_xg && !xg_name.empty()) {
        if (show_progress) {
            cerr << "Saving XG index to disk..." << endl;
        }
        // Save encapsulated in a VPKG
        vg::io::VPKG::save(*xg_index, xg_name); 
    }
    delete xg_index; xg_index = nullptr;

    // Build GCSA
    if (build_gcsa) {

        // Configure GCSA2 verbosity so it doesn't spit out loads of extra info
        if (!show_progress) {
            gcsa::Verbosity::set(gcsa::Verbosity::SILENT);
        }

        // Use the same temp directory as VG.
        gcsa::TempFile::setDirectory(temp_file::get_dir());

        double start = gcsa::readTimer();

        // Generate temporary kmer files
        bool delete_kmer_files = false;
        if (dbg_names.empty()) {
            if (show_progress) {
                cerr << "Generating kmer files..." << endl;
            }
            
            if (!file_names.empty()) {
                // Get the kmers from a VGset.
                VGset graphs(file_names);
                size_t kmer_bytes = params.getLimitBytes();
                dbg_names = graphs.write_gcsa_kmers_binary(kmer_size, kmer_bytes);
                params.reduceLimit(kmer_bytes);
                delete_kmer_files = true;
            } else if (!xg_name.empty()) {
                // Get the kmers from an XG
                
                get_input_file(xg_name, [&](istream& xg_stream) {
                    // Load the XG
                    auto xg = vg::io::VPKG::load_one<xg::XG>(xg_stream);
                
                    // Make an overlay on it to add source and sink nodes
                    // TODO: Don't use this directly; unify this code with VGset's code.
                    SourceSinkOverlay overlay(xg.get(), kmer_size);
                    
                    // Get the size limit
                    size_t kmer_bytes = params.getLimitBytes();
                    
                    // Write just the one kmer temp file
                    dbg_names.push_back(write_gcsa_kmers_to_tmpfile(overlay, kmer_size, kmer_bytes,
                        overlay.get_id(overlay.get_source_handle()),
                        overlay.get_id(overlay.get_sink_handle())));
                        
                    // Feed back into the size limit
                    params.reduceLimit(kmer_bytes);
                    delete_kmer_files = true;
                
                });
            } else {
                cerr << "error: [vg index] cannot generate GCSA index without either a vg or an xg" << endl;
                exit(1);
            }
        }

        // Build the index
        if (show_progress) {
            cerr << "Building the GCSA2 index..." << endl;
        }
        gcsa::InputGraph input_graph(dbg_names, true, gcsa::Alphabet(), mapping_name);
        gcsa::GCSA gcsa_index(input_graph, params);
        gcsa::LCPArray lcp_array(input_graph, params);
        if (show_progress) {
            double seconds = gcsa::readTimer() - start;
            cerr << "GCSA2 index built in " << seconds << " seconds, "
                 << gcsa::inGigabytes(gcsa::memoryUsage()) << " GB" << endl;
            cerr << "I/O volume: " << gcsa::inGigabytes(gcsa::readVolume()) << " GB read, "
                 << gcsa::inGigabytes(gcsa::writeVolume()) << " GB write" << endl;
        }

        // Save the indexes
        if (show_progress) {
            cerr << "Saving the index to disk..." << endl;
        }
        vg::io::VPKG::save(gcsa_index, gcsa_name);
        vg::io::VPKG::save(lcp_array, gcsa_name + ".lcp");

        // Verify the index
        if (verify_gcsa) {
            if (show_progress) {
                cerr << "Verifying the index..." << endl;
            }
            if (!gcsa::verifyIndex(gcsa_index, &lcp_array, input_graph)) {
                cerr << "warning: [vg index] GCSA2 index verification failed" << endl;
            }
        }

        // Delete the temporary kmer files
        if (delete_kmer_files) {
            for (auto& filename : dbg_names) {
                temp_file::remove(filename);
            }
        }
    }
    
    if (build_gai_index) {
        // Index a sorted GAM file.
        
        get_input_file(file_names.at(0), [&](istream& in) {
            // Grab the input GAM stream and wrap it in a cursor
            vg::io::ProtobufIterator<Alignment> cursor(in);
            
            // Index the file
            StreamIndex<Alignment> index;
            index.index(cursor);
 
            // Save the GAM index in the appropriate place.
            // TODO: Do we really like this enforced naming convention just beacuse samtools does it?
            ofstream index_out(file_names.at(0) + ".gai");
            if (!index_out.good()) {
                cerr << "error: [vg index] could not open " << file_names.at(0) << ".gai" << endl;
                exit(1);
            }
            index.save(index_out);
        });
    }
    
    if (build_vgi_index) {
        // Index an ID-sorted VG file.
        get_input_file(file_names.at(0), [&](istream& in) {
            // Grab the input VG stream and wrap it in a cursor
            vg::io::ProtobufIterator<Graph> cursor(in);
            
            // Index the file
            StreamIndex<Graph> index;
            index.index(cursor);
 
            // Save the index in the appropriate place.
            // TODO: Do we really like this enforced naming convention just beacuse samtools does it?
            ofstream index_out(file_names.at(0) + ".vgi");
            if (!index_out.good()) {
                cerr << "error: [vg index] could not open " << file_names.at(0) << ".vgi" << endl;
                exit(1);
            }
            index.save(index_out);
        });
        
    }

    if (build_rocksdb) {

        Index index;

        if (compact) {
            index.open_for_write(rocksdb_name);
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
                    vg::io::for_each(in, lambda);
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
                    vg::io::for_each(in, lambda);
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
                vg::io::write_buffered(cout, output_buf, 100);
            };
            index.for_each_alignment(lambda);
            vg::io::write_buffered(cout, output_buf, 0);
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
                    vg::io::for_each(in, lambda);
                });
            }
            index.flush();
            index.close();
        }

        if (dump_index) {
            index.open_read_only(rocksdb_name);
            index.dump(cout);
            index.close();
        }

    }


    //Build new snarl-based minimum distance index
    if (build_dist) {
        if (file_names.empty() && xg_name.empty()) {
            cerr << "error: [vg index] one graph is required to build a distance index" << endl;
            return 1;
        } else if (file_names.size() > 1 || (file_names.size() == 1 && !xg_name.empty())) {
            cerr << "error: [vg index] only one graph at a time can be used to build a distance index" << endl;
        } else if (dist_name.empty()) {
            cerr << "error: [vg index] distance index requires an output file" << endl;
            return 1;
        } else if (snarl_name.empty()) {
            cerr << "error: [vg index] distance index requires a snarl file" << endl;
            return 1;
            
        } else {
            //Get snarl manager
            ifstream snarl_stream(snarl_name);
            if (!snarl_stream) {
                cerr << "error: [vg index] cannot open Snarls file" << endl;
                exit(1);
            }
            SnarlManager* snarl_manager = new SnarlManager(snarl_stream);
            snarl_stream.close();

            //Get graph and build dist index
            if (file_names.empty() && !xg_name.empty()) {
                // We were given a -x specifically to read as XG
                
                ifstream xg_stream(xg_name);
                auto xg = vg::io::VPKG::load_one<xg::XG>(xg_stream);

                // Create the MinimumDistanceIndex
                MinimumDistanceIndex di(xg.get(), snarl_manager);
                // Save the completed DistanceIndex
                ofstream ostream(dist_name);
                di.serialize(ostream);

            } else {
                // We were given a graph generically
                auto graph = vg::io::VPKG::load_one<handlegraph::HandleGraph>(file_names.at(0));
    
                // Create the MinimumDistanceIndex
                MinimumDistanceIndex di(graph.get(), snarl_manager);
                // Save the completed DistanceIndex
                ofstream ostream(dist_name);
                di.serialize(ostream);
//                vg::io::VPKG::save(di, dist_name);
            }
          
            
        }

    }
    if (show_progress) {
        cerr << "Memory usage: " << gcsa::inGigabytes(gcsa::memoryUsage()) << " GB" << endl;
    }
    return 0;
}

std::vector<std::string> parseGenotypes(const std::string& vcf_line, size_t num_samples) {
    std::vector<std::string> result;

    // The 9th tab-separated field should start with "GT".
    size_t offset = 0;
    for (int i = 0; i < 8; i++) {
        size_t pos = vcf_line.find('\t', offset);
        if (pos == std::string::npos) {
            std::cerr << "error: [vg index] VCF line does not contain genotype information" << std::endl;
            std::exit(EXIT_FAILURE);
        }
        offset = pos + 1;
    }
    if (vcf_line.substr(offset, 2) != "GT") {
        std::cerr << "error: [vg index] VCF line does not contain genotype information" << std::endl;
        std::exit(EXIT_FAILURE);
    }

    // Genotype strings are the first colon-separated fields in the 10th+ tab-separated fields.
    offset = vcf_line.find('\t', offset);
    while (offset != std::string::npos && offset + 1 < vcf_line.length()) {
        offset++;
        size_t pos = vcf_line.find_first_of("\t:", offset);
        if (pos == std::string::npos) {
            pos = vcf_line.length();
        }
        result.emplace_back(vcf_line.substr(offset, pos - offset));
        offset = vcf_line.find('\t', offset);
    }

    if (result.size() != num_samples) {
        std::cerr << "error: [vg index] expected " << num_samples << " samples, got " << result.size() << std::endl;
        std::exit(EXIT_FAILURE);
    }

    return result;
}

// Register subcommand
static Subcommand vg_construct("index", "index graphs or alignments for random access or mapping", PIPELINE, 2, main_index);
