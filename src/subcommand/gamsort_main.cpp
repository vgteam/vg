/** \file gamsort_main.cpp
 *
 * Defines the "vg gamsort" subcommand for sorting and shuffling GAM and GAF files.
 */

#include "subcommand.hpp"

#include "../alignment.hpp"
#include "../gaf_sorter.hpp"
#include "../stream_index.hpp"
#include "../stream_sorter.hpp"

#include <vg/io/stream.hpp>
#include <getopt.h>

using namespace vg;
using namespace vg::subcommand;

//------------------------------------------------------------------------------

// We limit the max threads, and only allow thread count to be lowered, to
// prevent tcmalloc from giving each thread a very large heap for many
// threads. On my machine we can keep about 4 threads busy.
constexpr static size_t GAM_MAX_THREADS = 4;

// gaf_sorter defaults to 1 thread. If we assume that the input and the output
// are bgzip-compressed, we should be able to saturate 5 worker threads.
// Because intermediate merges are independent, we can use more threads.
// The final single-threaded merge can use 5 bgzip threads.

void help_gamsort(char **argv)
{
    std::cerr << "usage: " << argv[0] << " " << argv[1] << " [options] input > output" << std::endl;
    std::cerr << std::endl;
    std::cerr << "Sort a GAM/GAF file, or index a sorted GAM file." << std::endl;
    std::cerr << std::endl;
    std::cerr << "General options:" << std::endl;
    std::cerr << "    -p, --progress          show progress" << std::endl;
    std::cerr << "    -s, --shuffle           shuffle reads by hash" << std::endl;
    std::cerr << "    -t, --threads N         use N worker threads (default: " << GAM_MAX_THREADS << " for GAM, " << GAFSorterParameters::THREADS << " for GAF)" << std::endl;
    std::cerr << std::endl;
    std::cerr << "GAM sorting options:" << std::endl;
    std::cerr << "    -i, --index FILE        produce an index of the sorted GAM file" << std::endl;
    std::cerr << "    -d, --dumb-sort         use naive sorting algorithm (no tmp files, faster for small GAMs)" << std::endl;
    std::cerr << std::endl;
    std::cerr << "GAF sorting options:" << std::endl;
    std::cerr << "    -G, --gaf-input         input is a GAF file" << std::endl;
    std::cerr << "    -c, --chunk-size N      number of reads per chunk (default: " << GAFSorterParameters::RECORDS_PER_FILE << ")" << std::endl;
    std::cerr << "    -m, --merge-width N     number of files to merge at once (default: " << GAFSorterParameters::FILES_PER_MERGE << ")" << std::endl;
    std::cerr << "    -S, --stable            use stable sorting" << std::endl;
    std::cerr << std::endl;
}

//------------------------------------------------------------------------------

int main_gamsort(int argc, char **argv)
{
    // General options.
    string input_format = "GAM";

    // GAM sorting options.
    size_t num_threads = GAM_MAX_THREADS;
    string index_filename;
    bool easy_sort = false;
    bool shuffle = false;
    bool show_progress = false;

    // GAF sorting options.
    GAFSorterParameters gaf_params;

    int c;
    optind = 2; // force optind past command positional argument
    while (true)
    {
        static struct option long_options[] = {
            { "progress", no_argument, 0, 'p' },
            { "shuffle", no_argument, 0, 's' },
            { "threads", required_argument, 0, 't' },
            { "index", required_argument, 0, 'i' },
            { "dumb-sort", no_argument, 0, 'd' },
            { "gaf-input", no_argument, 0, 'G' },
            { "chunk-size", required_argument, 0, 'c' },
            { "merge-width", required_argument, 0, 'm' },
            { "stable", no_argument, 0, 'S' },
            { "help", no_argument, 0, 'h' },
            { 0, 0, 0, 0 }
        };
        int option_index = 0;
        c = getopt_long(argc, argv, "pst:i:dGc:m:Sh", long_options, &option_index);

        // Detect the end of the options.
        if (c == -1)
            break;

        switch (c)
        {
        // General options.
        case 'p':
            show_progress = true;
            gaf_params.progress = true;
            break;
        case 's':
            shuffle = true;
            gaf_params.key_type = GAFSorterRecord::key_hash;
            break;
        case 't':
            {
                size_t parsed = std::max(parse<size_t>(optarg), size_t(1));
                num_threads = std::min(parsed, num_threads);
                gaf_params.threads = parsed;
            }
            break;

        // GAM sorting options.
        case 'i':
            index_filename = optarg;
            break;
        case 'd':
            easy_sort = true;
            break;

        // GAF sorting options.
        case 'G':
            input_format = "GAF";
            break;
        case 'c':
            gaf_params.records_per_file = parse<size_t>(optarg);
            break;
        case 'm':
            gaf_params.files_per_merge = parse<size_t>(optarg);
            break;
        case 'S':
            gaf_params.stable = true;
            break;

        case 'h':
        case '?':
        default:
            help_gamsort(argv);
            exit(1);
        }
    }

    if (argc < 3){
        help_gamsort(argv);
        exit(1);
    }
    
    omp_set_num_threads(num_threads);

    if (input_format == "GAM") {
        if (shuffle && !index_filename.empty()) {
            cerr << "[vg gamsort] Indexing is not allowed when shuffling GAM files." << endl;
            exit(1);
        }
        get_input_file(optind, argc, argv, [&](istream& gam_in) {

            GAMSorter gs(shuffle ? GAMSorter::Order::RANDOM : GAMSorter::Order::BY_GRAPH_POSITION, show_progress);

            // Do a normal GAMSorter sort
            std::unique_ptr<GAMIndex> index;
        
            if (!index_filename.empty()) {
                // Make an index
                index = unique_ptr<GAMIndex>(new GAMIndex());
            }
        
            if (easy_sort) {
                // Sort in a single pass in memory
                gs.easy_sort(gam_in, std::cout, index.get());
            } else {
                // Sort using fan-in-limited temp file merging
                gs.stream_sort(gam_in, std::cout, index.get());
            }
        
            if (index.get() != nullptr) {
                // Save the index
                std::ofstream index_out(index_filename);
                index->save(index_out);
            }
        });
    } else if (input_format == "GAF") {
        std::string input_file = get_input_file_name(optind, argc, argv);
        std::string output_file = "-"; // output to stdout
        bool success = sort_gaf(input_file, output_file, gaf_params);
        if (!success) {
            std::exit(EXIT_FAILURE);
        }
    }

    return 0;
}

static Subcommand vg_gamsort("gamsort", "Sort a GAM/GAF file or index a sorted GAM file.", main_gamsort);
