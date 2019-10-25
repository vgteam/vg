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

#include "../xg.hpp"
#include "../gbwt_helper.hpp"
#include <vg/io/vpkg.hpp>
#include <bdsg/overlay_helper.hpp>

#include <gbwtgraph/gbwtgraph.h>
#include <gbwtgraph/path_cover.h>

using namespace std;
using namespace vg;
using namespace vg::subcommand;

#include <unistd.h>


void help_gbwt(char** argv) {
    std::cerr << "usage: " << argv[0] << " [options] [args]" << std::endl;
    std::cerr << std::endl;
    std::cerr << "Manipulate GBWTs." << std::endl;
    std::cerr << std::endl;
    std::cerr << "General options:" << std::endl;
    std::cerr << "    -o, --output X          write output GBWT to X (required with -m, -f, and -P)" << std::endl;
    std::cerr << "    -p, --progress          show progress and statistics" << std::endl;
    std::cerr << std::endl;
    std::cerr << "Merging (requires -o; use deps/gbwt/merge_gbwt for more options):" << std::endl;
    std::cerr << "    -m, --merge             merge the GBWT files from the input args and write to output" << std::endl;
    std::cerr << "    -f, --fast              fast merging algorithm (node ids must not overlap; implies -m)" << std::endl;
    std::cerr << std::endl;
    std::cerr << "Threads (one GBWT file as an input arg):" << std::endl;
    std::cerr << "    -c, --count-threads     print the number of threads" << std::endl;
    std::cerr << "    -e, --extract FILE      extract threads in SDSL format to FILE" << std::endl;
    std::cerr << std::endl;
    std::cerr << "GBWTGraph construction (0 or 1 GBWT files as input args):" << std::endl;
    std::cerr << "    -g, --graph-name FILE   build GBWTGraph and serialize it to FILE (requires -x)" << std::endl;
    std::cerr << "    -x, --xg-name FILE      use the node sequences from the graph in FILE" << std::endl;
    std::cerr << "    -l, --local-haplotypes  sample local haplotypes from the input (requires -o)" << std::endl;
    std::cerr << "    -P, --path-cover        build GBWT from a greedy path cover (requires -o)" << std::endl;
    std::cerr << "    -n, --num-paths N       find N paths per component in -l or -P (default" << gbwtgraph::PATH_COVER_DEFAULT_N << ")" << std::endl;
    std::cerr << "    -k, --context-length N  use N-node contexts in -l or -P (default " << gbwtgraph::PATH_COVER_DEFAULT_K << ")" << std::endl;
    std::cerr << "    -b, --buffer-size N     GBWT construction buffer size in millions of nodes (default " << (gbwt::DynamicGBWT::INSERT_BATCH_SIZE / gbwt::MILLION) << ")" << std::endl;
    std::cerr << "    -i, --id-interval N     store path ids at one out of N positions (default " << gbwt::DynamicGBWT::SAMPLE_INTERVAL << ")" << std::endl;
    std::cerr << std::endl;
    std::cerr << "Metadata (one GBWT file as an input arg; use deps/gbwt/metadata_tool to modify):" << std::endl;
    std::cerr << "    -M, --metadata          print basic metadata" << std::endl;
    std::cerr << "    -C, --contigs           print the number of contigs" << std::endl;
    std::cerr << "    -H, --haplotypes        print the number of haplotypes" << std::endl;
    std::cerr << "    -S, --samples           print the number of samples" << std::endl;
    std::cerr << "    -L, --list-names        list contig/sample names (use with -C or -S)" << std::endl;
    std::cerr << "    -T, --thread-names      list thread names" << std::endl;
    std::cerr << "    -R, --remove-sample X   remove sample X from the index (overwrites input if -o is not used)" << std::endl;
    std::cerr << std::endl;
}


int main_gbwt(int argc, char** argv)
{
    if (argc == 2) {
        help_gbwt(argv);
        std::exit(EXIT_FAILURE);
    }

    // Operation modes. Only one of these can be active.
    bool merge = false, build_graph = false, remove_sample = false, other_options = false;

    bool fast_merging = false;
    bool show_progress = false;
    bool count_threads = false;
    bool metadata = false, contigs = false, haplotypes = false, samples = false, list_names = false, thread_names = false;
    bool local_haplotypes = false, path_cover = false;
    size_t num_paths = gbwtgraph::PATH_COVER_DEFAULT_N, context_length = gbwtgraph::PATH_COVER_DEFAULT_K;
    size_t buffer_size = gbwt::DynamicGBWT::INSERT_BATCH_SIZE;
    size_t id_interval = gbwt::DynamicGBWT::SAMPLE_INTERVAL;
    string gbwt_output, thread_output;
    string graph_output, xg_name;
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
                { "list_names", no_argument, 0, 'L' },
                { "thread_names", no_argument, 0, 'T' },
                { "remove-sample", required_argument, 0, 'R' },

                { "help", no_argument, 0, 'h' },
                { 0, 0, 0, 0 }
            };

        int option_index = 0;
        c = getopt_long(argc, argv, "o:pmfce:g:x:lPn:k:b:i:MCHSLTR:h?", long_options, &option_index);

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
            merge = true;
            break;
        case 'f':
            fast_merging = true;
            merge = true;
            break;

        // Threads
        case 'c':
            count_threads = true;
            other_options = true;
            break;
        case 'e':
            thread_output = optarg;
            other_options = true;
            break;

        // GBWTGraph
        case 'g':
            graph_output = optarg;
            build_graph = true;
            break;
        case 'x':
            xg_name = optarg;
            break;
        case 'l':
            local_haplotypes = true;
            break;
        case 'P':
            path_cover = true;
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
            metadata = true;
            other_options = true;
            break;
        case 'C':
            contigs = true;
            other_options = true;
            break;
        case 'H':
            haplotypes = true;
            other_options = true;
            break;
        case 'S':
            samples = true;
            other_options = true;
            break;
        case 'L':
            list_names = true;
            break;
        case 'T':
            thread_names = true;
            other_options = true;
            break;
        case 'R':
            to_remove = optarg;
            remove_sample = true;
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
    size_t active_modes = 0;
    if (merge) {
        active_modes++;
    }
    if (build_graph) {
        active_modes++;
    }
    if (remove_sample) {
        active_modes++;
    }
    if (other_options) {
        active_modes++;
    }
    if (active_modes > 1) {
        std::cerr << "error: [vg gbwt] merging, graph construction, removing samples, and other options are mutually exclusive" << std::endl;
        std::exit(EXIT_FAILURE);
    }
    if (active_modes == 0) {
        help_gbwt(argv);
        std::exit(EXIT_FAILURE);
    }

    // Let GBWT operate silently.
    gbwt::Verbosity::set(gbwt::Verbosity::SILENT);


    // Merging options.
    if (merge) {

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
                string input_name = argv[optind];
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
                string input_name = argv[curr];
                std::unique_ptr<gbwt::GBWT> next = vg::io::VPKG::load_one<gbwt::GBWT>(input_name);
                if (next.get() == nullptr) {
                    std::cerr << "error: [vg gbwt]: could not load GBWT " << input_name << std::endl;
                    std::exit(EXIT_FAILURE);
                }
                if (show_progress) {
                    gbwt::printStatistics(*next, input_name);
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
    if (build_graph) {
        if (local_haplotypes || path_cover) {
            if (local_haplotypes && path_cover) {
                std::cerr << "error: [vg gbwt] cannot build both local haplotypes and path cover" << std::endl;
                std::exit(EXIT_FAILURE);
            }
            if (local_haplotypes && optind + 1 != argc) {
                std::cerr << "error: [vg gbwt] no input GBWT given for finding local haplotypes" << std::endl;
                std::exit(EXIT_FAILURE);
            }
            if (path_cover && optind != argc) {
                std::cerr << "error: [vg gbwt] input GBWTs are not used with path cover" << std::endl;
                std::exit(EXIT_FAILURE);
            }
            if (num_paths == 0) {
                std::cerr << "error: [vg gbwt] number of paths must be non-zero for -l and -P" << std::endl;
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
        if (path_cover) {
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
        } else {
            if (show_progress) {
                std::cerr << "Loading GBWT " << argv[optind] << std::endl;
            }
            loaded_gbwt = vg::io::VPKG::load_one<gbwt::GBWT>(argv[optind]);
            if (loaded_gbwt.get() == nullptr) {
                std::cerr << "error: [vg gbwt]: could not load GBWT " << argv[optind] << std::endl;
                std::exit(EXIT_FAILURE);
            }
            if (local_haplotypes) {
                if (show_progress) {
                    std::cerr << "Building temporary GBWTGraph" << std::endl;
                }
                gbwtgraph::GBWTGraph temp_graph(*loaded_gbwt, *handle_graph);
                if (show_progress) {
                    std::cerr << "Finding " << num_paths << "-path cover with context length " << context_length << std::endl;
                }
                double start = gbwt::readTimer();
                generated_gbwt = gbwtgraph::local_haplotypes(temp_graph, num_paths, context_length, buffer_size, id_interval, show_progress);
                if (show_progress) {
                    double seconds = gbwt::readTimer() - start;
                    std::cerr << "GBWT built in " << seconds << " seconds, " << gbwt::inGigabytes(gbwt::memoryUsage()) << " GiB" << std::endl;
                    std::cerr << "Serializing local haplotypes to " << gbwt_output << std::endl;
                }
                vg::io::VPKG::save(generated_gbwt, gbwt_output);
                loaded_gbwt.reset();
            }
        }
        const gbwt::GBWT& selected_gbwt = ((path_cover || local_haplotypes) ? generated_gbwt : *loaded_gbwt);

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
    if (remove_sample) {
        if (optind + 1 != argc) {
            std::cerr << "error: [vg gbwt] one input GBWT is required for removing a sample" << std::endl;
            std::exit(EXIT_FAILURE);
        }
        std::unique_ptr<gbwt::DynamicGBWT> index = vg::io::VPKG::load_one<gbwt::DynamicGBWT>(argv[optind]);
        if (index.get() == nullptr) {
            std::cerr << "error: [vg gbwt]: could not load dynamic GBWT " << argv[optind] << std::endl;
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


    // Other options.
    if (other_options) {
        if (optind + 1 != argc) {
            std::cerr << "error: [vg gbwt] selected options require one input GBWT" << std::endl;
            std::exit(EXIT_FAILURE);
        }
        std::unique_ptr<gbwt::GBWT> index = vg::io::VPKG::load_one<gbwt::GBWT>(argv[optind]);
        if (index.get() == nullptr) {
            std::cerr << "error: [vg gbwt]: could not load GBWT " << argv[optind] << std::endl;
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

