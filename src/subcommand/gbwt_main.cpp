/** \file gbwt_main.cpp
 *
 * Defines the "vg gbwt" subcommand, which wraps up access for commands we'd otherwise find
 * in the gbwt submodule.  This is handy if we want to be able to merge gbwt's without adding */

#include <omp.h>
#include <unistd.h>
#include <getopt.h>

#include <list>
#include <fstream>

#include "subcommand.hpp"

#include "../vg.hpp"

using namespace std;
using namespace vg;
using namespace vg::subcommand;

#include <unistd.h>

#include <gbwt/dynamic_gbwt.h>

using namespace gbwt;
using namespace std;

void help_gbwt(char** argv) {
    cerr << "usage: " << argv[0] << " add [options] [args]" << endl
         << "Manipuate GBWTs." << endl
         << "merging:" << endl
         << "    -m, --merge            Merge the GBWT files from the input args and write to output" << endl
         << "    -o, --output X         Use X as the base name for output (required)" << endl
         << "    -b, --batches N        Use batches of N sequences for merging (default: "
         << DynamicGBWT::MERGE_BATCH_SIZE << ")" << endl
         << "    -f, --fast             Fast merging algorithm (node ids must not overlap)" << endl
         << "    -p, --progress         Show progress and statistics" << endl;
}


int main_gbwt(int argc, char** argv)
{
    if (argc == 2) {
        help_gbwt(argv);
        return 1;
    }

    bool merge = false;
    size_type batch_size = DynamicGBWT::MERGE_BATCH_SIZE;
    bool fast_merging = false;
    bool show_progress = false;
    string output;

    int c;
    optind = 2; // force optind past command positional argument    
    while (true) {
        static struct option long_options[] =
            {
                {"merge", no_argument, 0, 'm'},
                {"output", required_argument, 0, 'o'},
                {"batches", required_argument, 0, 'b'},
                {"fast", no_argument, 0, 'f'},
                {"progress",  no_argument, 0, 'p'},                
                {"help", no_argument, 0, 'h'},
                {0, 0, 0, 0}
            };

        int option_index = 0;
        c = getopt_long(argc, argv, "mo:b:fph?", long_options, &option_index);

        /* Detect the end of the options. */
        if (c == -1)
            break;

        switch (c)
        {

        case 'm':
            merge = true;
            break;

        case 'o':
            output = optarg;
            break;
            
        case 'b':
            batch_size = std::stoul(optarg);
            break;

        case 'f':
            fast_merging = true;
            break;

        case 'p':
            show_progress = true;
            break;

        case 'h':
        case '?':
            /* getopt_long already printed an error message. */
            help_gbwt(argv);
            exit(1);
            break;

        default:
            abort ();
        }
    }

    if (merge) {
        size_type input_files = argc - optind;
        size_type total_inserted = 0;
        if (input_files <= 1) {
            cerr << "[vg gbwt] error: At least two input gbwt files required to merge" << endl;
            exit(1);
        }
        if (output.empty()) {
            cerr << "[vg gbwt] error: Output file must be specified with -o" << endl;
        }
        if (show_progress) {
            printHeader("Algorithm"); std::cout << (fast_merging ? "fast" : "insert") << std::endl;
            printHeader("Input files"); std::cout << input_files << std::endl;
            printHeader("Output name"); std::cout << output << std::endl;
            if(!fast_merging) { printHeader("Batch size"); std::cout << batch_size << std::endl; }
            std::cout << std::endl;
        }

        double start = readTimer();

        if(fast_merging)
        {
            std::vector<GBWT> indexes(argc - optind);
            for(int i = optind; i < argc; i++)
            {
                std::string input_name = argv[i];
                sdsl::load_from_file(indexes[i - optind], input_name + GBWT::EXTENSION);
                if (show_progress) {
                    printStatistics(indexes[i - optind], input_name);
                }
                total_inserted += indexes[i - optind].size();
            }
            GBWT merged(indexes);
            sdsl::store_to_file(merged, output + GBWT::EXTENSION);
            if (show_progress) {
                printStatistics(merged, output);
            }
        }
        else
        {
            DynamicGBWT index;
            {
                std::string input_name = argv[optind];
                sdsl::load_from_file(index, input_name + DynamicGBWT::EXTENSION);
                if (show_progress) {
                    printStatistics(index, input_name);
                }
                optind++;
            }
            while(optind < argc)
            {
                std::string input_name = argv[optind];
                GBWT next;
                sdsl::load_from_file(next, input_name + GBWT::EXTENSION);
                if (show_progress) {
                    printStatistics(next, input_name);
                }
                index.merge(next, batch_size);
                total_inserted += next.size();
                optind++;
            }
            sdsl::store_to_file(index, output + DynamicGBWT::EXTENSION);
            if (show_progress) { 
                printStatistics(index, output);
            }
        }
        
        double seconds = readTimer() - start;

        if (show_progress) {

            std::cout << "Inserted " << total_inserted << " nodes in " << seconds << " seconds ("
                      << (total_inserted / seconds) << " nodes/second)" << std::endl;
            std::cout << "Memory usage " << inGigabytes(memoryUsage()) << " GB" << std::endl;
            std::cout << std::endl;
        }
        
        return 0;
    }

    return 0;
}


// Register subcommand
static Subcommand vg_gbwt("gbwt", "Manipuate GBWTs", main_gbwt);

