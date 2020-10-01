/** \file gbwt_main.cpp
 *
 * Defines the "vg gbwt" subcommand, which wraps up access for commands we'd otherwise find
 * in the gbwt submodule.  */

#include <omp.h>
#include <unistd.h>
#include <getopt.h>

#include <list>
#include <fstream>

#include "subcommand.hpp"

#include "../gbwt_helper.hpp"
#include <vg/io/vpkg.hpp>
#include <bdsg/overlays/overlay_helper.hpp>

#include <gbwt/fast_locate.h>
#include <gbwtgraph/gbwtgraph.h>
#include <gbwtgraph/path_cover.h>

using namespace std;
using namespace vg;
using namespace vg::subcommand;

#include <unistd.h>

enum merge_mode { merge_none, merge_insert, merge_fast };
enum path_cover_mode { path_cover_none, path_cover_augment, path_cover_local, path_cover_greedy };
enum index_type { index_none, index_compressed, index_dynamic };

void load_gbwt(const std::string& filename, gbwt::GBWT& index, bool show_progress);
void load_gbwt(const std::string& filename, gbwt::DynamicGBWT& index, bool show_progress);

void get_compressed(gbwt::GBWT& compressed_index, gbwt::DynamicGBWT& dynamic_index, index_type& in_use, const std::string& filename, bool show_progress);
void get_dynamic(gbwt::GBWT& compressed_index, gbwt::DynamicGBWT& dynamic_index, index_type& in_use, const std::string& filename, bool show_progress);

void get_graph(std::unique_ptr<HandleGraph>& graph, bool& in_use, const std::string& filename, bool show_progress);
void clear_graph(std::unique_ptr<HandleGraph>& graph, bool& in_use);

void help_gbwt(char** argv) {
    std::cerr << "usage: " << argv[0] << " gbwt [options] [args]" << std::endl;
    std::cerr << std::endl;
    std::cerr << "Manipulate GBWTs. Each input arg represents one input GBWT." << std::endl;
    std::cerr << std::endl;
    std::cerr << "General options:" << std::endl;
    std::cerr << "    -x, --xg-name FILE      read the graph from FILE" << std::endl;
    std::cerr << "    -o, --output FILE       write output GBWT to FILE" << std::endl;
    std::cerr << "    -p, --progress          show progress and statistics" << std::endl;
    std::cerr << std::endl;
    std::cerr << "Step 1: Merge multiple input GBWTs (requires -o; use deps/gbwt/merge_gbwt for more options):" << std::endl;
    std::cerr << "    -m, --merge             use the insertion algorithm" << std::endl;
    std::cerr << "    -f, --fast              fast merging algorithm (node ids must not overlap)" << std::endl;
    std::cerr << std::endl;
    std::cerr << "Step 2: Remove samples (requires -o and one input GBWT):" << std::endl;
    std::cerr << "    -R, --remove-sample X   remove sample X from the index" << std::endl;
    std::cerr << std::endl;
    std::cerr << "Step 3: Path cover (requires -o, -x, and one of { -a, -l, -P }):" << std::endl;
    std::cerr << "    -a, --augment-gbwt      add a path cover of missing components (one input GBWT)" << std::endl;
    std::cerr << "    -l, --local-haplotypes  sample local haplotypes (one input GBWT)" << std::endl;
    std::cerr << "    -P, --path-cover        build a greedy path cover (no input GBWTs)" << std::endl;
    std::cerr << "    -n, --num-paths N       find N paths per component (default " << gbwtgraph::PATH_COVER_DEFAULT_N << ")" << std::endl;
    std::cerr << "    -k, --context-length N  use N-node contexts (default " << gbwtgraph::PATH_COVER_DEFAULT_K << ")" << std::endl;
    std::cerr << "    -b, --buffer-size N     GBWT construction buffer size in millions of nodes (default " << (gbwt::DynamicGBWT::INSERT_BATCH_SIZE / gbwt::MILLION) << ")" << std::endl;
    std::cerr << "    -i, --id-interval N     store path ids at one out of N positions (default " << gbwt::DynamicGBWT::SAMPLE_INTERVAL << ")" << std::endl;
    std::cerr << std::endl;
    std::cerr << "Step 4: GBWTGraph construction (requires -x and one input GBWT):" << std::endl;
    std::cerr << "    -g, --graph-name FILE   build GBWTGraph and store it in FILE" << std::endl;
    std::cerr << std::endl;
    std::cerr << "Step 5: R-index construction (one input GBWT):" << std::endl;
    std::cerr << "    -r, --r-index FILE      build an r-index and store it in FILE" << std::endl;
    std::cerr << "    -t, --threads N         use N extraction threads (default " << omp_get_max_threads() << ")" << std::endl;
    std::cerr << std::endl;
    std::cerr << "Step 6: Metadata (one input GBWT; use deps/gbwt/metadata_tool to modify):" << std::endl;
    std::cerr << "    -M, --metadata          print basic metadata" << std::endl;
    std::cerr << "    -C, --contigs           print the number of contigs" << std::endl;
    std::cerr << "    -H, --haplotypes        print the number of haplotypes" << std::endl;
    std::cerr << "    -S, --samples           print the number of samples" << std::endl;
    std::cerr << "    -L, --list-names        list contig/sample names (use with -C or -S)" << std::endl;
    std::cerr << "    -T, --thread-names      list thread names" << std::endl;
    std::cerr << std::endl;
    std::cerr << "Step 7: Threads (one input GBWT):" << std::endl;
    std::cerr << "    -c, --count-threads     print the number of threads" << std::endl;
    std::cerr << "    -e, --extract FILE      extract threads in SDSL format to FILE" << std::endl;
    std::cerr << std::endl;
}


int main_gbwt(int argc, char** argv)
{
    if (argc == 2) {
        help_gbwt(argv);
        std::exit(EXIT_FAILURE);
    }

    // Requirements and modes.
    bool produces_one_gbwt = false;
    merge_mode merge = merge_none;
    path_cover_mode path_cover = path_cover_none;
    bool metadata_mode = false, thread_mode = false;

    // Other parameters and flags.
    bool show_progress = false;
    bool count_threads = false;
    bool metadata = false, contigs = false, haplotypes = false, samples = false, list_names = false, thread_names = false;
    size_t num_paths = gbwtgraph::PATH_COVER_DEFAULT_N, context_length = gbwtgraph::PATH_COVER_DEFAULT_K;
    size_t buffer_size = gbwt::DynamicGBWT::INSERT_BATCH_SIZE;
    size_t id_interval = gbwt::DynamicGBWT::SAMPLE_INTERVAL;

    // File/sample names.
    string gbwt_output, thread_output;
    string graph_output, xg_name;
    std::string r_index_name;
    std::string to_remove;

    int c;
    optind = 2; // force optind past command positional argument
    while (true) {
        static struct option long_options[] =
            {
                // General
                { "xg-name", required_argument, 0, 'x' },
                { "output", required_argument, 0, 'o' },
                { "progress",  no_argument, 0, 'p' },

                // Merging
                { "merge", no_argument, 0, 'm' },
                { "fast", no_argument, 0, 'f' },

                // Remove sample
                { "remove-sample", required_argument, 0, 'R' },

                // Path cover
                { "augment-gbwt", no_argument, 0, 'a' },
                { "local-haplotypes", no_argument, 0, 'l' },
                { "path-cover", no_argument, 0, 'P' },
                { "num-paths", required_argument, 0, 'n' },
                { "context-length", required_argument, 0, 'k' },
                { "buffer-size", required_argument, 0, 'b' },
                { "id-interval", required_argument, 0, 'i' },

                // GBWTGraph
                { "graph-name", required_argument, 0, 'g' },

                // R-index
                { "r-index", required_argument, 0, 'r' },
                { "threads", required_argument, 0, 't' },

                // Metadata
                { "metadata", no_argument, 0, 'M' },
                { "contigs", no_argument, 0, 'C' },
                { "haplotypes", no_argument, 0, 'H' },
                { "samples", no_argument, 0, 'S' },
                { "list-names", no_argument, 0, 'L' },
                { "thread-names", no_argument, 0, 'T' },

                // Threads
                { "count-threads", no_argument, 0, 'c' },
                { "extract", required_argument, 0, 'e' },

                { "help", no_argument, 0, 'h' },
                { 0, 0, 0, 0 }
            };

        int option_index = 0;
        c = getopt_long(argc, argv, "x:o:pmfR:alPn:k:b:i:g:r:t:MCHSLTce:h?", long_options, &option_index);

        /* Detect the end of the options. */
        if (c == -1)
            break;

        switch (c)
        {
        // General
        case 'x':
            xg_name = optarg;
            break;
        case 'o':
            gbwt_output = optarg;
            break;
        case 'p':
            show_progress = true;
            break;

        // Merging
        case 'm':
            merge = merge_insert;
            produces_one_gbwt = true;
            break;
        case 'f':
            merge = merge_fast;
            produces_one_gbwt = true;
            break;

        // Remove sample
        case 'R':
            to_remove = optarg;
            break;

        // Path cover
        case 'a':
            assert(path_cover == path_cover_none);
            path_cover = path_cover_augment;
            break;
        case 'l':
            assert(path_cover == path_cover_none);
            path_cover = path_cover_local;
            break;
        case 'P':
            assert(path_cover == path_cover_none);
            path_cover = path_cover_greedy;
            produces_one_gbwt = true;
            break;
        case 'n':
            num_paths = parse<size_t>(optarg);
            break;
        case 'k':
            context_length = parse<size_t>(optarg);
            break;
        case 'b':
            buffer_size = parse<size_t>(optarg) * gbwt::MILLION;
            break;
        case 'i':
            id_interval = parse<size_t>(optarg);
            break;

        // GBWTGraph
        case 'g':
            graph_output = optarg;
            break;

        // Build r-index
        case 'r':
            r_index_name = optarg;
            break;
        case 't':
            omp_set_num_threads(parse<int>(optarg));
            break;

        // Metadata
        case 'M':
            metadata = true;
            metadata_mode = true;
            break;
        case 'C':
            contigs = true;
            metadata_mode = true;
            break;
        case 'H':
            haplotypes = true;
            metadata_mode = true;
            break;
        case 'S':
            samples = true;
            metadata_mode = true;
            break;
        case 'L':
            list_names = true;
            metadata_mode = true;
            break;
        case 'T':
            thread_names = true;
            metadata_mode = true;
            break;

        // Threads
        case 'c':
            count_threads = true;
            thread_mode = true;
            break;
        case 'e':
            thread_output = optarg;
            thread_mode = true;
            break;

        case 'h':
        case '?':
            /* getopt_long already printed an error message. */
            help_gbwt(argv);
            std::exit(EXIT_FAILURE);
            break;

        default:
            std::exit(EXIT_FAILURE);
        }
    }
    int input_gbwts = argc - optind;
    std::string gbwt_name = (input_gbwts == 0 ? "" : argv[optind]);

    // Sanity checks.
    if (merge != merge_none) {
        if (input_gbwts < 2 || gbwt_output.empty()) {
            std::cerr << "error: [vg gbwt]: merging requires multiple input GBWTs and output GBWT" << std::endl;
            std::exit(EXIT_FAILURE);
        }
    }
    if (!to_remove.empty()) {
        if (input_gbwts != 1 || gbwt_output.empty()) {
            std::cerr << "error: [vg gbwt]: removing a sample requires one input GBWT and output GBWT" << std::endl;
            std::exit(EXIT_FAILURE);
        }
    }
    if (path_cover != path_cover_none) {
        if (gbwt_output.empty() || xg_name.empty()) {
            std::cerr << "error: [vg gbwt]: path cover options require output GBWT and a graph" << std::endl;
            std::exit(EXIT_FAILURE);
        }
        if (path_cover == path_cover_greedy && input_gbwts != 0) {
            std::cerr << "error: [vg gbwt]: greedy path cover does not use input GBWTs" << std::endl;
            std::exit(EXIT_FAILURE);
        }
        if ((path_cover == path_cover_local || path_cover == path_cover_augment) && !(input_gbwts == 1 || produces_one_gbwt)) {
            std::cerr << "error: [vg gbwt]: path cover options -a and -l require one input GBWT" << std::endl;
            std::exit(EXIT_FAILURE);
        }
        if (num_paths == 0) {
            std::cerr << "error: [vg gbwt] number of paths must be non-zero for path cover" << std::endl;
            std::exit(EXIT_FAILURE);
        }
        if (context_length < gbwtgraph::PATH_COVER_MIN_K) {
            std::cerr << "error: [vg gbwt] context length must be at least " << gbwtgraph::PATH_COVER_MIN_K << " for path cover" << std::endl;
            std::exit(EXIT_FAILURE);
        }
        if (buffer_size == 0) {
            std::cerr << "error: [vg gbwt] GBWT construction buffer size cannot be 0" << std::endl;
            std::exit(EXIT_FAILURE);
        }
    }
    if (!graph_output.empty()) {
        if (xg_name.empty() || !(input_gbwts == 1 || produces_one_gbwt)) {
            std::cerr << "error: [vg gbwt]: GBWTGraph construction requires a graph and one input GBWT" << std::endl;
            std::exit(EXIT_FAILURE);
        }
    }
    if (!r_index_name.empty()) {
        if (!(input_gbwts == 1 || produces_one_gbwt)) {
            std::cerr << "error: [vg gbwt]: r-index construction requires one input GBWT" << std::endl;
            std::exit(EXIT_FAILURE);
        }
    }
    if (metadata_mode) {
        if (!(input_gbwts == 1 || produces_one_gbwt)) {
            std::cerr << "error: [vg gbwt]: metadata operations one input GBWT" << std::endl;
            std::exit(EXIT_FAILURE);
        }
    }
    if (thread_mode) {
        if (!(input_gbwts == 1 || produces_one_gbwt)) {
            std::cerr << "error: [vg gbwt]: thread operations one input GBWT" << std::endl;
            std::exit(EXIT_FAILURE);
        }
    }


    // Let GBWT operate silently and use the same temporary directory as vg.
    gbwt::Verbosity::set(gbwt::Verbosity::SILENT);
    gbwt::TempFile::setDirectory(temp_file::get_dir());

    // This is the data we are using.
    gbwt::GBWT compressed_index;
    gbwt::DynamicGBWT dynamic_index;
    index_type in_use = index_none;
    std::unique_ptr<HandleGraph> input_graph;
    bool graph_in_use = false;


    // Merge multiple input GBWTs.
    if (merge != merge_none) {
        double start = gbwt::readTimer();
        if (show_progress) {
            std::string algo_name = (merge == merge_fast ? "fast" : "insertion");
            std::cerr << "Merging " << input_gbwts << " input GBWTs (" << algo_name << " algorithm)" << std::endl;
        }
        if (merge == merge_fast) {
            std::vector<gbwt::GBWT> indexes(input_gbwts);
            for (int i = optind; i < argc; i++) {
                load_gbwt(argv[i], indexes[i - optind], show_progress);
            }
            if (show_progress) {
                std::cerr << "Merging the GBWTs" << std::endl;
            }
            compressed_index = gbwt::GBWT(indexes);
            in_use = index_compressed;
        }
        else if (merge == merge_insert) {
            load_gbwt(gbwt_name, dynamic_index, show_progress);
            in_use = index_dynamic;
            for (int curr = optind + 1; curr < argc; curr++) {
                gbwt::GBWT next;
                load_gbwt(argv[curr], next, show_progress);
                if (next.size() > 2 * dynamic_index.size()) {
                    std::cerr << "warning: [vg gbwt] merging " << argv[curr] << " into a substantially smaller index" << std::endl;
                    std::cerr << "warning: [vg gbwt] merging would be faster in another order" << std::endl;
                }
                if (show_progress) {
                    std::cerr << "Inserting " << next.sequences() << " sequences of total length " << next.size() << std::endl;
                }
                dynamic_index.merge(next);
            }
        }
        if (show_progress) {
            double seconds = gbwt::readTimer() - start;
            std::cerr << "GBWTs merged in " << seconds << " seconds, " << gbwt::inGigabytes(gbwt::memoryUsage()) << " GiB" << std::endl;
            std::cerr << std::endl;
        }
    }


    // Remove a sample from the GBWT.
    if (!to_remove.empty()) {
        double start = gbwt::readTimer();
        if (show_progress) {
            std::cerr << "Removing sample " << to_remove << " from the index" << std::endl;
        }
        get_dynamic(compressed_index, dynamic_index, in_use, gbwt_name, show_progress);
        if (!(dynamic_index.hasMetadata() && dynamic_index.metadata.hasPathNames() && dynamic_index.metadata.hasSampleNames())) {
            std::cerr << "error: [vg gbwt] the index does not contain metadata with thread and sample names" << std::endl;
            std::exit(EXIT_FAILURE);
        }
        gbwt::size_type sample_id = dynamic_index.metadata.sample(to_remove);
        std::vector<gbwt::size_type> path_ids = dynamic_index.metadata.removeSample(sample_id);
        if (path_ids.empty()) {
            std::cerr << "error: [vg gbwt] the index does not contain sample " << to_remove << std::endl;
            std::exit(EXIT_FAILURE);
        }
        if (show_progress) {
            std::cerr << "Removing " << path_ids.size() << " threads" << std::endl;
        }
        size_t foo = dynamic_index.remove(path_ids);
        if (show_progress) {
            double seconds = gbwt::readTimer() - start;
            std::cerr << "Sample removed in " << seconds << " seconds, " << gbwt::inGigabytes(gbwt::memoryUsage()) << " GiB" << std::endl;
            std::cerr << std::endl;
        }
    }


    // Path cover options.
    if (path_cover != path_cover_none) {
        double start = gbwt::readTimer();
        if (show_progress) {
            std::cerr << "Finding a " << num_paths << "-path cover with context length " << context_length << std::endl;
        }
        get_graph(input_graph, graph_in_use, xg_name, show_progress);
        if (path_cover == path_cover_greedy) {
            if (show_progress) {
                std::cerr << "Algorithm: greedy" << std::endl;
            }
            compressed_index = gbwtgraph::path_cover_gbwt(*input_graph, num_paths, context_length, buffer_size, id_interval, show_progress);
            in_use = index_compressed;
        } else if (path_cover == path_cover_augment) {
            if (show_progress) {
                std::cerr << "Algorithm: augment" << std::endl;
            }
            get_dynamic(compressed_index, dynamic_index, in_use, gbwt_name, show_progress);
            gbwtgraph::augment_gbwt(*input_graph, dynamic_index, num_paths, context_length, buffer_size, id_interval, show_progress);
        } else {
            if (show_progress) {
                std::cerr << "Algorithm: local haplotypes" << std::endl;
            }
            get_compressed(compressed_index, dynamic_index, in_use, gbwt_name, show_progress);
            compressed_index = gbwtgraph::local_haplotypes(*input_graph, compressed_index, num_paths, context_length, buffer_size, id_interval, show_progress);
        }
        if (show_progress) {
            double seconds = gbwt::readTimer() - start;
            std::cerr << "Path cover built in " << seconds << " seconds, " << gbwt::inGigabytes(gbwt::memoryUsage()) << " GiB" << std::endl;
            std::cerr << std::endl;
        }
    }


    // Now we can serialize the GBWT.
    if (!gbwt_output.empty() && in_use != index_none) {
        double start = gbwt::readTimer();
        if (show_progress) {
            std::cerr << "Serializing the GBWT to " << gbwt_output << std::endl;
        }
        if (in_use == index_compressed) {
            vg::io::VPKG::save(compressed_index, gbwt_output);
        } else if (in_use == index_dynamic) {
            vg::io::VPKG::save(dynamic_index, gbwt_output);
        }
        if (show_progress) {
            double seconds = gbwt::readTimer() - start;
            std::cerr << "GBWT serialized in " << seconds << " seconds, " << gbwt::inGigabytes(gbwt::memoryUsage()) << " GiB" << std::endl;
            std::cerr << std::endl;
        }
    }


    // GBWTGraph construction.
    if (!graph_output.empty()) {
        double start = gbwt::readTimer();
        if (show_progress) {
            std::cerr << "Building GBWTGraph" << std::endl;
        }
        get_graph(input_graph, graph_in_use, xg_name, show_progress);
        get_compressed(compressed_index, dynamic_index, in_use, gbwt_name, show_progress);
        if (show_progress) {
            std::cerr << "Starting the construction" << std::endl;
        }
        gbwtgraph::GBWTGraph graph(compressed_index, *input_graph);
        if (show_progress) {
            std::cerr << "Serializing GBWTGraph to " << graph_output << std::endl;
        }
        vg::io::VPKG::save(graph, graph_output);
        if (show_progress) {
            double seconds = gbwt::readTimer() - start;
            std::cerr << "GBWTGraph built in " << seconds << " seconds, " << gbwt::inGigabytes(gbwt::memoryUsage()) << " GiB" << std::endl;
            std::cerr << std::endl;
        }
    }


    // We no longer need the graph.
    clear_graph(input_graph, graph_in_use);


    // R-index construction.
    if (!r_index_name.empty()) {
        double start = gbwt::readTimer();
        if (show_progress) {
            std::cerr << "Building r-index" << std::endl;
        }
        get_compressed(compressed_index, dynamic_index, in_use, gbwt_name, show_progress);
        if (show_progress) {
            std::cerr << "Starting the construction" << std::endl;
        }
        gbwt::FastLocate r_index(compressed_index);
        if (show_progress) {
            std::cerr << "Serializing the r-index to " << r_index_name << std::endl;
        }
        vg::io::VPKG::save(r_index, r_index_name);
        if (show_progress) {
            double seconds = gbwt::readTimer() - start;
            std::cerr << "R-index built in " << seconds << " seconds, " << gbwt::inGigabytes(gbwt::memoryUsage()) << " GiB" << std::endl;
            std::cerr << std::endl;
        }
    }


    // Metadata options.
    if (metadata_mode) {
        get_compressed(compressed_index, dynamic_index, in_use, gbwt_name, show_progress);
        if (!compressed_index.hasMetadata()) {
            std::cerr << "error: [vg gbwt] the GBWT does not contain metadata" << std::endl;
            std::exit(EXIT_FAILURE);
        }
        if (metadata) {
            gbwt::operator<<(std::cout, compressed_index.metadata) << std::endl;
        }
        if (contigs) {
            if (list_names) {
                if (compressed_index.metadata.hasContigNames()) {
                    for (size_t i = 0; i < compressed_index.metadata.contigs(); i++) {
                        std::cout << compressed_index.metadata.contig(i) << std::endl;
                    }
                } else {
                    std::cerr << "error: [vg gbwt] the metadata does not contain contig names" << std::endl;
                    std::exit(EXIT_FAILURE);
                }
            } else {
                std::cout << compressed_index.metadata.contigs() << std::endl;
            }
        }
        if (haplotypes) {
            std::cout << compressed_index.metadata.haplotypes() << std::endl;
        }
        if (samples) {
            if (list_names) {
                if (compressed_index.metadata.hasSampleNames()) {
                    for (size_t i = 0; i < compressed_index.metadata.samples(); i++) {
                        std::cout << compressed_index.metadata.sample(i) << std::endl;
                    }
                } else {
                    std::cerr << "error: [vg gbwt] the metadata does not contain sample names" << std::endl;
                    std::exit(EXIT_FAILURE);
                }
            } else {
                std::cout << compressed_index.metadata.samples() << std::endl;
            }
        }
        if (thread_names) {
            if (compressed_index.metadata.hasPathNames()) {
                for (size_t i = 0; i < compressed_index.metadata.paths(); i++) {
                    std::cout << thread_name(compressed_index, i) << std::endl;
                }
            } else {
                std::cerr << "error: [vg gbwt] the metadata does not contain thread names" << std::endl;
            }
        }
    }

    // Thread options.
    if (thread_mode) {
        // Extract threads in SDSL format.
        if (!thread_output.empty()) {
            double start = gbwt::readTimer();
            if (show_progress) {
                std::cerr << "Extracting threads to " << thread_output << std::endl;
            }
            get_compressed(compressed_index, dynamic_index, in_use, gbwt_name, show_progress);
            if (show_progress) {
                std::cerr << "Starting the extraction" << std::endl;
            }
            gbwt::size_type node_width = gbwt::bit_length(compressed_index.sigma() - 1);
            gbwt::text_buffer_type out(thread_output, std::ios::out, gbwt::MEGABYTE, node_width);
            for (gbwt::size_type id = 0; id < compressed_index.sequences(); id += 2) { // Ignore reverse complements.
                gbwt::vector_type sequence = compressed_index.extract(id);
                for (auto node : sequence) {
                    out.push_back(node);
                }
                out.push_back(gbwt::ENDMARKER);
            }
            out.close();
            if (show_progress) {
                double seconds = gbwt::readTimer() - start;
                std::cerr << "Threads extracted in " << seconds << " seconds, " << gbwt::inGigabytes(gbwt::memoryUsage()) << " GiB" << std::endl;
                std::cerr << std::endl;
            }
        }

        // There are two sequences for each thread.
        if (count_threads) {
            get_compressed(compressed_index, dynamic_index, in_use, gbwt_name, show_progress);
            std::cout << (compressed_index.sequences() / 2) << std::endl;
        }
    }

    return 0;
}


//----------------------------------------------------------------------------

// Utility functions

void load_gbwt(const std::string& filename, gbwt::GBWT& index, bool show_progress) {
    if (show_progress) {
        std::cerr << "Loading " << filename << std::endl;
    }
    std::unique_ptr<gbwt::GBWT> loaded = vg::io::VPKG::load_one<gbwt::GBWT>(filename);
    if (loaded.get() == nullptr) {
        std::cerr << "error: [vg gbwt] could not load compressed GBWT " << filename << std::endl;
        std::exit(EXIT_FAILURE);
    }
    index = std::move(*loaded);
}

void load_gbwt(const std::string& filename, gbwt::DynamicGBWT& index, bool show_progress) {
    if (show_progress) {
        std::cerr << "Loading " << filename << std::endl;
    }
    std::unique_ptr<gbwt::DynamicGBWT> loaded = vg::io::VPKG::load_one<gbwt::DynamicGBWT>(filename);
    if (loaded.get() == nullptr) {
        std::cerr << "error: [vg gbwt] could not load dynamic GBWT " << filename << std::endl;
        std::exit(EXIT_FAILURE);
    }
    index = std::move(*loaded);
}

void get_compressed(gbwt::GBWT& compressed_index, gbwt::DynamicGBWT& dynamic_index, index_type& in_use, const std::string& filename, bool show_progress) {
    if (in_use == index_compressed) {
        return;
    } else if (in_use == index_dynamic) {
        if (show_progress) {
            std::cerr << "Converting dynamic GBWT into compressed GBWT" << std::endl;
        }
        compressed_index = gbwt::GBWT(dynamic_index);
        dynamic_index = gbwt::DynamicGBWT();
        in_use = index_compressed;
    } else {
        load_gbwt(filename, compressed_index, show_progress);
        in_use = index_compressed;
    }
}

void get_dynamic(gbwt::GBWT& compressed_index, gbwt::DynamicGBWT& dynamic_index, index_type& in_use, const std::string& filename, bool show_progress) {
    if (in_use == index_dynamic) {
        return;
    } else if (in_use == index_compressed) {
        if (show_progress) {
            std::cerr << "Converting compressed GBWT into dynamic GBWT" << std::endl;
        }
        dynamic_index = gbwt::DynamicGBWT(compressed_index);
        compressed_index = gbwt::GBWT();
        in_use = index_dynamic;
    } else {
        load_gbwt(filename, dynamic_index, show_progress);
        in_use = index_dynamic;
    }
}

void get_graph(std::unique_ptr<HandleGraph>& graph, bool& in_use, const std::string& filename, bool show_progress) {
    if (in_use) {
        return;
    } else {
        if (show_progress) {
            std::cerr << "Loading input graph from " << filename << std::endl;
        }
        graph = vg::io::VPKG::load_one<HandleGraph>(filename);
        if (graph == nullptr) {
            std::cerr << "error: [vg gbwt] could not load graph " << filename << std::endl;
            std::exit(EXIT_FAILURE);
        }
        in_use = true;
    }
}

void clear_graph(std::unique_ptr<HandleGraph>& graph, bool& in_use) {
    if (!in_use) {
        return;
    } else {
        graph.reset();
        in_use = false;
    }
}

//----------------------------------------------------------------------------


// Register subcommand
static Subcommand vg_gbwt("gbwt", "Manipulate GBWTs", main_gbwt);

