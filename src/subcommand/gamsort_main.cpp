#include "../stream_sorter.hpp"
#include <vg/io/stream.hpp>
#include "../stream_index.hpp"
#include <getopt.h>
#include "subcommand.hpp"

/**
* GAM sort main
*/

using namespace std;
using namespace vg;
using namespace vg::subcommand;
void help_gamsort(char **argv)
{
    cerr << "gamsort: sort a GAM file, or index a sorted GAM file" << endl
         << "Usage: " << argv[1] << " [Options] gamfile" << endl
         << "Options:" << endl
         << "  -i / --index FILE       produce an index of the sorted GAM file" << endl
         << "  -d / --dumb-sort        use naive sorting algorithm (no tmp files, faster for small GAMs)" << endl
         << "  -p / --progress         Show progress." << endl
         << "  -t / --threads          Use the specified number of threads." << endl
         << endl;
}

int main_gamsort(int argc, char **argv)
{
    string index_filename;
    bool easy_sort = false;
    bool show_progress = false;
    // We limit the max threads, and only allow thread count to be lowered, to
    // prevent tcmalloc from giving each thread a very large heap for many
    // threads.
    // On my machine we can keep about 4 threads busy.
    size_t num_threads = 4;
    int c;
    optind = 2; // force optind past command positional argument
    while (true)
    {
        static struct option long_options[] =
            {
                {"index", required_argument, 0, 'i'},
                {"dumb-sort", no_argument, 0, 'd'},
                {"rocks", required_argument, 0, 'r'},
                {"progress", no_argument, 0, 'p'},
                {"threads", required_argument, 0, 't'},
                {0, 0, 0, 0}};
        int option_index = 0;
        c = getopt_long(argc, argv, "i:dhpt:",
                        long_options, &option_index);

        // Detect the end of the options.
        if (c == -1)
            break;

        switch (c)
        {
        case 'i':
            index_filename = optarg;
            break;
        case 'd':
            easy_sort = true;
            break;
        case 'p':
            show_progress = true;
            break;
        case 't':
            num_threads = min(parse<size_t>(optarg), num_threads);
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

    get_input_file(optind, argc, argv, [&](istream& gam_in) {

        GAMSorter gs(show_progress);

        // Do a normal GAMSorter sort
        unique_ptr<GAMIndex> index;
        
        if (!index_filename.empty()) {
            // Make an index
            index = unique_ptr<GAMIndex>(new GAMIndex());
        }
        
        if (easy_sort) {
            // Sort in a single pass in memory
            gs.easy_sort(gam_in, cout, index.get());
        } else {
            // Sort using fan-in-limited temp file merging
            gs.stream_sort(gam_in, cout, index.get());
        }
        
        if (index.get() != nullptr) {
            // Save the index
            ofstream index_out(index_filename);
            index->save(index_out);
        }
    });

    return 0;
}

static Subcommand vg_gamsort("gamsort", "Sort a GAM file or index a sorted GAM file.", main_gamsort);
