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

enum operation_mode { operation_none, operation_merge, operation_graph, operation_remove, operation_r_index, operation_other };
enum path_cover_mode { path_cover_none, path_cover_augment, path_cover_local, path_cover_greedy };

void help_gbwt(char** argv) {
    std::cerr << "usage: " << argv[0] << " [options] [args]" << std::endl;
    std::cerr << std::endl;
    std::cerr << "Manipulate GBWTs." << std::endl;
    std::cerr << std::endl;
    std::cerr << "General options:" << std::endl;
    std::cerr << "    -o, --output FILE       write output GBWT to FILE" << std::endl;
    std::cerr << "    -p, --progress          show progress and statistics" << std::endl;
    std::cerr << std::endl;
    std::cerr << "Merging (requires -o; use deps/gbwt/merge_gbwt for more options):" << std::endl;
    std::cerr << "    -m, --merge             merge the GBWT files from the input args and write to output" << std::endl;
    std::cerr << "    -f, --fast              fast merging algorithm (node ids must not overlap; implies -m)" << std::endl;
    std::cerr << std::endl;
    std::cerr << "Threads (one input GBWT):" << std::endl;
    std::cerr << "    -c, --count-threads     print the number of threads" << std::endl;
    std::cerr << "    -e, --extract FILE      extract threads in SDSL format to FILE" << std::endl;
    std::cerr << std::endl;
    std::cerr << "GBWTGraph construction (based on one input GBWT or path cover GBWT):" << std::endl;
    std::cerr << "    -g, --graph-name FILE   build GBWTGraph and serialize it to FILE (required)" << std::endl;
    std::cerr << "    -x, --xg-name FILE      use the node sequences from the graph in FILE (required)" << std::endl;
    std::cerr << std::endl;
    std::cerr << "Path cover (for GBWTGraph construction; requires -o and one of { -a, -l, -P }):" << std::endl;
    std::cerr << "    -a, --augment-gbwt      add a path cover of missing components (one input GBWT)" << std::endl;
    std::cerr << "    -l, --local-haplotypes  sample local haplotypes (one input GBWT)" << std::endl;
    std::cerr << "    -P, --path-cover        build a greedy path cover (no input GBWTs)" << std::endl;
    std::cerr << "    -n, --num-paths N       find N paths per component (default " << gbwtgraph::PATH_COVER_DEFAULT_N << ")" << std::endl;
    std::cerr << "    -k, --context-length N  use N-node contexts (default " << gbwtgraph::PATH_COVER_DEFAULT_K << ")" << std::endl;
    std::cerr << "    -b, --buffer-size N     GBWT construction buffer size in millions of nodes (default " << (gbwt::DynamicGBWT::INSERT_BATCH_SIZE / gbwt::MILLION) << ")" << std::endl;
    std::cerr << "    -i, --id-interval N     store path ids at one out of N positions (default " << gbwt::DynamicGBWT::SAMPLE_INTERVAL << ")" << std::endl;
    std::cerr << std::endl;
    std::cerr << "Metadata (one input GBWT; use deps/gbwt/metadata_tool to modify):" << std::endl;
    std::cerr << "    -M, --metadata          print basic metadata" << std::endl;
    std::cerr << "    -C, --contigs           print the number of contigs" << std::endl;
    std::cerr << "    -H, --haplotypes        print the number of haplotypes" << std::endl;
    std::cerr << "    -S, --samples           print the number of samples" << std::endl;
    std::cerr << "    -L, --list-names        list contig/sample names (use with -C or -S)" << std::endl;
    std::cerr << "    -T, --thread-names      list thread names" << std::endl;
    std::cerr << std::endl;
    std::cerr << "Remove sample (one input GBWT; requires -o):" << std::endl;
    std::cerr << "    -R, --remove-sample X   remove sample X from the index" << std::endl;
    std::cerr << std::endl;
    std::cerr << "R-index construction (one input GBWT):" << std::endl;
    std::cerr << "    -r, --r-index FILE      build an r-index and store it in FILE" << std::endl;
    std::cerr << "    -t, --threads N         use N extraction threads (default " << omp_get_max_threads() << ")" << std::endl;
    std::cerr << std::endl;
}


int main_gbwt(int argc, char** argv)
{
    if (argc == 2) {
        help_gbwt(argv);
        std::exit(EXIT_FAILURE);
    }

    // Operation modes.
    operation_mode mode = operation_none;
    path_cover_mode path_cover = path_cover_none;

    bool fast_merging = false;
    bool show_progress = false;
    bool count_threads = false;
    bool metadata = false, contigs = false, haplotypes = false, samples = false, list_names = false, thread_names = false;
    size_t num_paths = gbwtgraph::PATH_COVER_DEFAULT_N, context_length = gbwtgraph::PATH_COVER_DEFAULT_K;
    size_t buffer_size = gbwt::DynamicGBWT::INSERT_BATCH_SIZE;
    size_t id_interval = gbwt::DynamicGBWT::SAMPLE_INTERVAL;
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
                { "output", required_argument, 0, 'o' },
                { "progress",  no_argument, 0, 'p' },

                // Merging
                { "merge", no_argument, 0, 'm' },
                { "fast", no_argument, 0, 'f' },

                // Threads
                { "count-threads", no_argument, 0, 'c' },
                { "extract", required_argument, 0, 'e' },

                // GBWTGraph
                { "graph-name", required_argument, 0, 'g' },
                { "xg-name", required_argument, 0, 'x' },

                // Path cover
                { "augment-gbwt", no_argument, 0, 'a' },
                { "local-haplotypes", no_argument, 0, 'l' },
                { "path-cover", no_argument, 0, 'P' },
                { "num-paths", required_argument, 0, 'n' },
                { "context-length", required_argument, 0, 'k' },
                { "buffer-size", required_argument, 0, 'b' },
                { "id-interval", required_argument, 0, 'i' },

                // Metadata
                { "metadata", no_argument, 0, 'M' },
                { "contigs", no_argument, 0, 'C' },
                { "haplotypes", no_argument, 0, 'H' },
                { "samples", no_argument, 0, 'S' },
                { "list-names", no_argument, 0, 'L' },
                { "thread-names", no_argument, 0, 'T' },

                // Remove sample
                { "remove-sample", required_argument, 0, 'R' },

                // R-index
                { "r-index", required_argument, 0, 'r' },
                { "threads", required_argument, 0, 't' },

                { "help", no_argument, 0, 'h' },
                { 0, 0, 0, 0 }
            };

        int option_index = 0;
        c = getopt_long(argc, argv, "o:pmfce:g:x:alPn:k:b:i:MCHSLTR:r:t:h?", long_options, &option_index);

        /* Detect the end of the options. */
        if (c == -1)
            break;

        switch (c)
        {
        // General
        case 'o':
            gbwt_output = optarg;
            break;
        case 'p':
            show_progress = true;
            break;

        // Merging
        case 'm':
            mode = operation_merge;
            break;
        case 'f':
            mode = operation_merge;
            fast_merging = true;
            break;

        // Threads
        case 'c':
            mode = operation_other;
            count_threads = true;
            break;
        case 'e':
            mode = operation_other;
            thread_output = optarg;
            break;

        // GBWTGraph
        case 'g':
            mode = operation_graph;
            graph_output = optarg;
            break;
        case 'x':
            xg_name = optarg;
            break;

        // Path cover
        case 'a':
            path_cover = path_cover_augment;
            break;
        case 'l':
            path_cover = path_cover_local;
            break;
        case 'P':
            path_cover = path_cover_greedy;
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

        // Metadata
        case 'M':
            mode = operation_other;
            metadata = true;
            break;
        case 'C':
            mode = operation_other;
            contigs = true;
            break;
        case 'H':
            mode = operation_other;
            haplotypes = true;
            break;
        case 'S':
            mode = operation_other;
            samples = true;
            break;
        case 'L':
            list_names = true;
            break;
        case 'T':
            mode = operation_other;
            thread_names = true;
            break;

        // Remove sample
        case 'R':
            mode = operation_remove;
            to_remove = optarg;
            break;

        // Build r-index
        case 'r':
            mode = operation_r_index;
            r_index_name = optarg;
            break;
        case 't':
            omp_set_num_threads(parse<int>(optarg));
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

    // Sanity checks.
    if (mode == operation_none) {
        help_gbwt(argv);
        std::exit(EXIT_FAILURE);
    }

    // Let GBWT operate silently.
    gbwt::Verbosity::set(gbwt::Verbosity::SILENT);


    // Merging options.
    if (mode == operation_merge) {

        // Ugly hack here. GBWT prints to stdout, and we want to direct it to stderr.
        std::streambuf* cout_buf = std::cout.rdbuf();
        std::cout.rdbuf(cerr.rdbuf());

        size_t input_files = argc - optind;
        size_t total_inserted = 0;
        if (input_files <= 1) {
            std::cerr << "error: [vg gbwt] at least two input gbwt files required to merge" << std::endl;
            std::exit(EXIT_FAILURE);
        }
        if (gbwt_output.empty()) {
            std::cerr << "error: [vg gbwt] output file must be specified with -o" << std::endl;
            std::exit(EXIT_FAILURE);
        }
        if (show_progress) {
            gbwt::printHeader("Algorithm"); std::cout << (fast_merging ? "fast" : "insert") << std::endl;
            gbwt::printHeader("Input files"); std::cout << input_files << std::endl;
            gbwt::printHeader("Output name"); std::cout << gbwt_output << std::endl;
            std::cout << std::endl;
        }

        double start = gbwt::readTimer();
        if(fast_merging)
        {
            vector<gbwt::GBWT> indexes(argc - optind);
            for(int i = optind; i < argc; i++)
            {
                string input_name = argv[i];
                
                // Try loading the GBWT
                std::unique_ptr<gbwt::GBWT> loaded = vg::io::VPKG::load_one<gbwt::GBWT>(input_name);
                if (loaded.get() == nullptr) {
                    std::cerr << "error: [vg gbwt] could not load GBWT " << input_name << std::endl;
                    std::exit(EXIT_FAILURE);
                }
                
                // Move out of the std::unique_ptr and into the vector that the GBWT library needs.
                indexes[i - optind] = std::move(*loaded);
                
                if (show_progress) {
                    gbwt::printStatistics(indexes[i - optind], input_name);
                }
                total_inserted += indexes[i - optind].size();
            }
            
            // Merge the GBWT.
            gbwt::GBWT merged(indexes);
            
            // Save to a file in VPKG-encapsulated format.
            vg::io::VPKG::save(merged, gbwt_output);
            if (show_progress) {
                gbwt::printStatistics(merged, gbwt_output);
            }
        }
        else
        {
            std::unique_ptr<gbwt::DynamicGBWT> index;
            {
                // Try to load the first GBWT as a dynamic one.
                std::string input_name = argv[optind];
                index = vg::io::VPKG::load_one<gbwt::DynamicGBWT>(input_name);
                if (index.get() == nullptr) {
                    std::cerr << "error: [vg gbwt] could not load dynamic GBWT " << input_name << std::endl;
                    std::exit(EXIT_FAILURE);
                }
                if (show_progress) {
                    gbwt::printStatistics(*index, input_name);
                }
            }
            for (int curr = optind + 1; curr < argc; curr++)
            {
                std::string input_name = argv[curr];
                std::unique_ptr<gbwt::GBWT> next = vg::io::VPKG::load_one<gbwt::GBWT>(input_name);
                if (next.get() == nullptr) {
                    std::cerr << "error: [vg gbwt]: could not load GBWT " << input_name << std::endl;
                    std::exit(EXIT_FAILURE);
                }
                if (show_progress) {
                    gbwt::printStatistics(*next, input_name);
                }
                if (next->size() > 2 * index->size()) {
                    std::cerr << "warning: [vg gbwt] merging " << input_name << " into a substantially smaller index" << std::endl;
                    std::cerr << "warning: [vg gbwt] merging would be faster in another order" << std::endl;
                }
                index->merge(*next);
                total_inserted += next->size();
            }
            
            // Save to a file in VPKG-encapsulated format.
            vg::io::VPKG::save(*index, gbwt_output);
            if (show_progress) { 
                gbwt::printStatistics(*index, gbwt_output);
            }
        }
        double seconds = gbwt::readTimer() - start;

        if (show_progress) {
            std::cout << "Inserted " << total_inserted << " nodes in " << seconds << " seconds ("
                      << (total_inserted / seconds) << " nodes/second)" << std::endl;
            std::cout << "Memory usage " << gbwt::inGigabytes(gbwt::memoryUsage()) << " GB" << std::endl;
            std::cout << std::endl;
        }

        // Revert the hack and exit.
        std::cout.rdbuf(cout_buf);
        std::exit(EXIT_SUCCESS);
    }


    // GBWTGraph construction.
    if (mode == operation_graph) {
        if (path_cover != path_cover_none) {
            if (path_cover == path_cover_local && optind + 1 != argc) {
                std::cerr << "error: [vg gbwt] no input GBWT given for finding local haplotypes" << std::endl;
                std::exit(EXIT_FAILURE);
            }
            if (path_cover == path_cover_augment && optind + 1 != argc) {
                std::cerr << "error: [vg gbwt] no input GBWT to augment" << std::endl;
                std::exit(EXIT_FAILURE);
            }
            if (path_cover == path_cover_greedy && optind != argc) {
                std::cerr << "error: [vg gbwt] input GBWTs are not used with greedy path cover" << std::endl;
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
            if (gbwt_output.empty()) {
                std::cerr << "error: [vg gbwt] output file not specified for path cover GBWT" << std::endl;
                std::exit(EXIT_FAILURE);
            }
        } else {
            if (optind + 1 != argc) {
                std::cerr << "error: [vg gbwt] no input GBWT given for GBWTGraph construction" << std::endl;
                std::exit(EXIT_FAILURE);
            }
        }
        if (xg_name.empty()) {
            std::cerr << "error: [vg gbwt] GBWTGraph construction requires an input graph" << std::endl;
            std::exit(EXIT_FAILURE);
        }

        // Load the input graph.
        if (show_progress) {
            std::cerr << "Loading input graph " << xg_name << std::endl;
        }
        std::unique_ptr<HandleGraph> handle_graph = vg::io::VPKG::load_one<HandleGraph>(xg_name);
        if (handle_graph == nullptr) {
            std::cerr << "error: [vg gbwt] could not load graph " << xg_name << std::endl;
            std::exit(EXIT_FAILURE);
        }

        // Load or build GBWT.
        std::unique_ptr<gbwt::GBWT> loaded_gbwt(nullptr);
        gbwt::GBWT generated_gbwt;
        if (path_cover == path_cover_greedy) {
            if (show_progress) {
                std::cerr << "Finding " << num_paths << "-path cover with context length " << context_length << std::endl;
            }
            double start = gbwt::readTimer();
            generated_gbwt = gbwtgraph::path_cover_gbwt(*handle_graph, num_paths, context_length, buffer_size, id_interval, show_progress);
            if (show_progress) {
                double seconds = gbwt::readTimer() - start;
                std::cerr << "GBWT built in " << seconds << " seconds, " << gbwt::inGigabytes(gbwt::memoryUsage()) << " GiB" << std::endl;
                std::cerr << "Serializing path cover GBWT to " << gbwt_output << std::endl;
            }
            vg::io::VPKG::save(generated_gbwt, gbwt_output);
        } else if (path_cover == path_cover_augment) {
            if (show_progress) {
                std::cerr << "Loading GBWT " << argv[optind] << std::endl;
            }
            std::unique_ptr<gbwt::DynamicGBWT> input_gbwt = vg::io::VPKG::load_one<gbwt::DynamicGBWT>(argv[optind]);
            if (input_gbwt.get() == nullptr) {
                std::cerr << "error: [vg gbwt] could not load GBWT " << argv[optind] << std::endl;
                std::exit(EXIT_FAILURE);
            }
            if (show_progress) {
                std::cerr << "Finding " << num_paths << "-path cover with context length " << context_length << " for missing components" << std::endl;
            }
            double start = gbwt::readTimer();
            gbwtgraph::augment_gbwt(*handle_graph, *input_gbwt, num_paths, context_length, buffer_size, id_interval, show_progress);
            if (show_progress) {
                double seconds = gbwt::readTimer() - start;
                std::cerr << "GBWT augmented in " << seconds << " seconds, " << gbwt::inGigabytes(gbwt::memoryUsage()) << " GiB" << std::endl;
                std::cerr << "Serializing augmented GBWT to " << gbwt_output << std::endl;
            }
            generated_gbwt = gbwt::GBWT(*input_gbwt);
            vg::io::VPKG::save(generated_gbwt, gbwt_output);
        } else {
            if (show_progress) {
                std::cerr << "Loading GBWT " << argv[optind] << std::endl;
            }
            loaded_gbwt = vg::io::VPKG::load_one<gbwt::GBWT>(argv[optind]);
            if (loaded_gbwt.get() == nullptr) {
                std::cerr << "error: [vg gbwt] could not load GBWT " << argv[optind] << std::endl;
                std::exit(EXIT_FAILURE);
            }
            if (path_cover == path_cover_local) {
                if (show_progress) {
                    std::cerr << "Finding " << num_paths << "-path cover with context length " << context_length << std::endl;
                }
                double start = gbwt::readTimer();
                generated_gbwt = gbwtgraph::local_haplotypes(*handle_graph, *loaded_gbwt, num_paths, context_length, buffer_size, id_interval, show_progress);
                if (show_progress) {
                    double seconds = gbwt::readTimer() - start;
                    std::cerr << "GBWT built in " << seconds << " seconds, " << gbwt::inGigabytes(gbwt::memoryUsage()) << " GiB" << std::endl;
                    std::cerr << "Serializing local haplotypes to " << gbwt_output << std::endl;
                }
                vg::io::VPKG::save(generated_gbwt, gbwt_output);
                loaded_gbwt.reset();
            }
        }
        const gbwt::GBWT& selected_gbwt = (path_cover == path_cover_none ? *loaded_gbwt : generated_gbwt);

        if (show_progress) {
            std::cerr << "Building GBWTGraph" << std::endl;
        }
        gbwtgraph::GBWTGraph graph(selected_gbwt, *handle_graph);
        if (show_progress) {
            std::cerr << "Serializing GBWTGraph to " << graph_output << std::endl;
        }
        vg::io::VPKG::save(graph, graph_output);
        std::exit(EXIT_SUCCESS);
    }


    // Remove a sample from the GBWT.
    if (mode == operation_remove) {
        if (optind + 1 != argc) {
            std::cerr << "error: [vg gbwt] one input GBWT is required for removing a sample" << std::endl;
            std::exit(EXIT_FAILURE);
        }
        if (gbwt_output.empty()) {
            std::cerr << "error: [vg gbwt] output file not specified for removing a sample" << std::endl;
            std::exit(EXIT_FAILURE);
        }
        std::unique_ptr<gbwt::DynamicGBWT> index = vg::io::VPKG::load_one<gbwt::DynamicGBWT>(argv[optind]);
        if (index.get() == nullptr) {
            std::cerr << "error: [vg gbwt] could not load dynamic GBWT " << argv[optind] << std::endl;
            std::exit(EXIT_FAILURE);
        }
        if (!(index->hasMetadata() && index->metadata.hasPathNames() && index->metadata.hasSampleNames())) {
            std::cerr << "error: [vg gbwt] the index does not contain metadata with thread and sample names" << std::endl;
            std::exit(EXIT_FAILURE);
        }
        gbwt::size_type sample_id = index->metadata.sample(to_remove);
        std::vector<gbwt::size_type> path_ids = index->metadata.removeSample(sample_id);
        if (path_ids.empty()) {
            std::cerr << "error: [vg gbwt] the index does not contain sample " << to_remove << std::endl;
            std::exit(EXIT_FAILURE);
        }
        size_t foo = index->remove(path_ids);
        std::string output = (gbwt_output.empty() ? argv[optind] : gbwt_output);
        vg::io::VPKG::save(*index, output);
        std::exit(EXIT_SUCCESS);
    }


    // R-index construction.
    if (mode == operation_r_index) {
        if (show_progress) {
            std::cerr << "Loading GBWT index " << argv[optind] << std::endl;
        }
        std::unique_ptr<gbwt::GBWT> index = vg::io::VPKG::load_one<gbwt::GBWT>(argv[optind]);
        if (show_progress) {
            std::cerr << "Building the r-index" << std::endl;
        }
        gbwt::FastLocate r_index(*index);
        if (show_progress) {
            std::cerr << "Saving the r-index to " << r_index_name << std::endl;
        }
        vg::io::VPKG::save(r_index, r_index_name);
    }


    // Other options.
    if (mode == operation_other) {
        if (optind + 1 != argc) {
            std::cerr << "error: [vg gbwt] selected options require one input GBWT" << std::endl;
            std::exit(EXIT_FAILURE);
        }
        std::unique_ptr<gbwt::GBWT> index = vg::io::VPKG::load_one<gbwt::GBWT>(argv[optind]);
        if (index.get() == nullptr) {
            std::cerr << "error: [vg gbwt] could not load GBWT " << argv[optind] << std::endl;
            std::exit(EXIT_FAILURE);
        }

        // Extract threads in SDSL format.
        if (!thread_output.empty()) {
            gbwt::size_type node_width = gbwt::bit_length(index->sigma() - 1);
            gbwt::text_buffer_type out(thread_output, std::ios::out, gbwt::MEGABYTE, node_width);
            for (gbwt::size_type id = 0; id < index->sequences(); id += 2) { // Ignore reverse complements.
                gbwt::vector_type sequence = index->extract(id);
                for (auto node : sequence) {
                    out.push_back(node);
                }
                out.push_back(gbwt::ENDMARKER);
            }
            out.close();
        }

        // There are two sequences for each thread.
        if (count_threads) {
            std::cout << (index->sequences() / 2) << std::endl;
        }

        if (metadata) {
            if (index->hasMetadata()) {
                gbwt::operator<<(std::cout, index->metadata) << std::endl;
            }
        }

        if (contigs) {
            if (index->hasMetadata()) {
                if (list_names) {
                    if (index->metadata.hasContigNames()) {
                        for (size_t i = 0; i < index->metadata.contigs(); i++) {
                            std::cout << index->metadata.contig(i) << std::endl;
                        }
                    } else {
                        std::cerr << "error: [vg gbwt] the metadata does not contain contig names" << std::endl;
                    }
                } else {
                    std::cout << index->metadata.contigs() << std::endl;
                }
            } else {
                std::cout << -1 << std::endl;
            }
        }

        if (haplotypes) {
            if (index->hasMetadata()) {
                std::cout << index->metadata.haplotypes() << std::endl;
            } else {
                std::cout << -1 << std::endl;
            }
        }

        if (samples) {
            if (index->hasMetadata()) {
                if (list_names) {
                    if (index->metadata.hasSampleNames()) {
                        for (size_t i = 0; i < index->metadata.samples(); i++) {
                            std::cout << index->metadata.sample(i) << std::endl;
                        }
                    } else {
                        std::cerr << "error: [vg gbwt] the metadata does not contain sample names" << std::endl;
                    }
                } else {
                    std::cout << index->metadata.samples() << std::endl;
                }
            } else {
                std::cout << -1 << std::endl;
            }
        }

        if (thread_names) {
            if (index->hasMetadata() && index->metadata.hasPathNames()) {
                for (size_t i = 0; i < index->metadata.paths(); i++) {
                    std::cout << thread_name(*index, i) << std::endl;
                }
            } else {
                std::cerr << "error: [vg gbwt] the metadata does not contain thread names" << std::endl;
            }
        }
    }

    return 0;
}


// Register subcommand
static Subcommand vg_gbwt("gbwt", "Manipulate GBWTs", main_gbwt);

