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
    std::cerr << "Merging (use deps/gbwt/merge_gbwt for more options):" << std::endl;
    std::cerr << "    -m, --merge            merge the GBWT files from the input args and write to output" << std::endl;
    std::cerr << "    -o, --output X         write output GBWT to X (required)" << std::endl;
    std::cerr << "    -f, --fast             fast merging algorithm (node ids must not overlap; implies -m)" << std::endl;
    std::cerr << "    -p, --progress         show progress and statistics" << std::endl;
    std::cerr << "Threads (one GBWT file as an input arg):" << std::endl;
    std::cerr << "    -c, --count-threads    print the number of threads" << std::endl;
    std::cerr << "    -e, --extract FILE     extract threads in SDSL format to FILE" << std::endl;
    std::cerr << "GBWTGraph construction (0 or 1 GBWT files as input args):" << std::endl;
    std::cerr << "    -g, --graph-name FILE  build a GBWT graph and serialize it to FILE (requires -x)" << std::endl;
    std::cerr << "    -x, --xg-name FILE     use the node sequences from the graph in FILE" << std::endl;
    // FIXME path cover
    std::cerr << "Metadata (one GBWT file as an input arg; use deps/gbwt/metadata_tool to modify):" << std::endl;
    std::cerr << "    -M, --metadata         print basic metadata" << std::endl;
    std::cerr << "    -C, --contigs          print the number of contigs" << std::endl;
    std::cerr << "    -H, --haplotypes       print the number of haplotypes" << std::endl;
    std::cerr << "    -S, --samples          print the number of samples" << std::endl;
    std::cerr << "    -L, --list-names       list contig/sample names (use with -C or -S)" << std::endl;
    std::cerr << "    -T, --thread-names     list thread names" << std::endl;
    std::cerr << "    -R, --remove-sample X  remove sample X from the index (use -o to change output)" << std::endl;
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
    string gbwt_output, thread_output;
    string graph_output, xg_name;
    std::string to_remove;

    int c;
    optind = 2; // force optind past command positional argument    
    while (true) {
        static struct option long_options[] =
            {
                // Merging
                {"merge", no_argument, 0, 'm'},
                {"output", required_argument, 0, 'o'},
                {"fast", no_argument, 0, 'f'},
                {"progress",  no_argument, 0, 'p'},

                // Threads
                {"count-threads", no_argument, 0, 'c'},
                {"extract", required_argument, 0, 'e'},

                // GBWTGraph
                {"graph-name", required_argument, 0, 'g'},
                {"xg-name", required_argument, 0, 'x'},

                // Metadata
                {"metadata", no_argument, 0, 'M'},
                {"contigs", no_argument, 0, 'C'},
                {"haplotypes", no_argument, 0, 'H'},
                {"samples", no_argument, 0, 'S'},
                {"list_names", no_argument, 0, 'L'},
                {"thread_names", no_argument, 0, 'T'},
                {"remove-sample", required_argument, 0, 'R'},

                {"help", no_argument, 0, 'h'},
                {0, 0, 0, 0}
            };

        int option_index = 0;
        c = getopt_long(argc, argv, "mo:fpce:g:x:MCHSLTR:h?", long_options, &option_index);

        /* Detect the end of the options. */
        if (c == -1)
            break;

        switch (c)
        {
        // Merging
        case 'm':
            merge = true;
            break;
        case 'o':
            gbwt_output = optarg;
            break;
        case 'f':
            fast_merging = true;
            merge = true;
            break;
        case 'p':
            show_progress = true;
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
                unique_ptr<gbwt::GBWT> loaded = vg::io::VPKG::load_one<gbwt::GBWT>(input_name);
                if (loaded.get() == nullptr) {
                    std::cerr << "error: [vg gbwt] could not load GBWT " << input_name << std::endl;
                    std::exit(EXIT_FAILURE);
                }
                
                // Move out of the unique_ptr and into the vector that the GBWT library needs.
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
            unique_ptr<gbwt::DynamicGBWT> index;
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
                unique_ptr<gbwt::GBWT> next = vg::io::VPKG::load_one<gbwt::GBWT>(input_name);
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
        if (optind + 1 != argc) {
            std::cerr << "error: [vg gbwt] GBWTGraph construction requires one input GBWT" << std::endl;
            std::exit(EXIT_FAILURE);
        }
        if (xg_name.empty()) {
            std::cerr << "error: [vg gbwt] GBWTGraph construction requires an input graph" << std::endl;
            std::exit(EXIT_FAILURE);
        }
        unique_ptr<gbwt::GBWT> index = vg::io::VPKG::load_one<gbwt::GBWT>(argv[optind]);
        if (index.get() == nullptr) {
            std::cerr << "error: [vg gbwt]: could not load GBWT " << argv[optind] << std::endl;
            std::exit(EXIT_FAILURE);
        }
        unique_ptr<HandleGraph> handle_graph = vg::io::VPKG::load_one<HandleGraph>(xg_name);
        if (handle_graph == nullptr) {
            std::cerr << "error: [vg gbwt] could not load graph " << xg_name << std::endl;
            std::exit(EXIT_FAILURE);
        }
        gbwtgraph::GBWTGraph graph(*index, *handle_graph);
        vg::io::VPKG::save(graph, graph_output);
        std::exit(EXIT_SUCCESS);
    }


    // Remove a sample from the GBWT.
    if (remove_sample) {
        if (optind + 1 != argc) {
            std::cerr << "error: [vg gbwt] one input GBWT is required for removing a sample" << std::endl;
            std::exit(EXIT_FAILURE);
        }
        unique_ptr<gbwt::DynamicGBWT> index = vg::io::VPKG::load_one<gbwt::DynamicGBWT>(argv[optind]);
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
        unique_ptr<gbwt::GBWT> index = vg::io::VPKG::load_one<gbwt::GBWT>(argv[optind]);
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

