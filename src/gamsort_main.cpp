#include "GAMSorter.hpp"
#include "stream.hpp"
#include "getopt.h"

/**
* GAM sort main
*/

void gamsort_help(char** argv){
    cerr << "gamsort: sort a GAM file (or index it) without Rocksdb" << endl
    << "Usage: " << argv[1] << " [Options] gamfile" << endl
    << "Options:" << endl
    << "  -i / --index <INDEX>    produce a node-to-alignment index"
    << "  -d / --dumb-sort      use naive sorting algorithm (no tmp files, faster for small GAMs)"
    << endl;
}

int gam_sort_main(int argc, char** argv){
    string gamfile;
    string index_file;
    bool dumb_sort = false;
    bool is_paired = false;
    int c;
    optind = 2; // force optind past command positional argument
    while (true) {
        static struct option long_options[] =
        {
            {"index", required_argument, 0, 'i'},
            {"dumb-sort", no_argument, 0, 'd'},
            {"paired", no_argument, 0, 'p'},
            {0, 0, 0, 0}

        };
        int option_index = 0;
        c = getopt_long (argc, argv, "i:dh",
                long_options, &option_index);

        // Detect the end of the options.
        if (c == -1)
            break;

        switch (c)
        {
            case 'i':
                index_file = optarg;
                break;
            case 'd':
                dumb_sort = true;
            case 'h':
            case '?':
            default:
                help_srpe(argv);
                abort();
        }

    }

    gamfile = argv[optind];

    GAMSorter gs;

    if (dumb_sort){
        gs.dumb_sort(gamfile);
    }

    if (!index_file.empty()){
        gs.write_index(gamfile, gamfile + ".gai", false);
    }

    return 1;
}