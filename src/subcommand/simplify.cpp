// simplify.cpp: define the "vg simplify" subcommand, which removes small variation

#include <omp.h>
#include <unistd.h>
#include <getopt.h>

#include "subcommand.hpp"

#include "../vg.hpp"
// This provides the CactusSiteFinder
#include "../genotypekit.hpp"


using namespace std;
using namespace vg;
using namespace vg::subcommand;

void help_simplify(char** argv) {
    cerr << "usage: " << argv[0] << " simplify [options] <old.vg >new.vg" << endl
         << "options:" << endl
         << "    -m, --min-size N      remove bubbles with fewer than N bases involved (default: 10)" << endl
         << "    -p, --progress        show progress" << endl
         << "    -t, --threads N       use N threads to construct graph (defaults to numCPUs)" << endl;
}

int main_simplify(int argc, char** argv) {

    if (argc == 2) {
        help_simplify(argv);
        return 1;
    }


    size_t min_size = 10;
    bool show_progress = false;

    int c;
    while (true) {
        static struct option long_options[] =
            {
                {"min-size", required_argument, 0, 'm'},
                {"progress",  no_argument, 0, 'p'},
                {"threads", required_argument, 0, 't'},
                {"help", no_argument, 0, 'h'},
                {0, 0, 0, 0}
            };

        int option_index = 0;
        c = getopt_long (argc, argv, "m:pt:h?",
                         long_options, &option_index);

        /* Detect the end of the options. */
        if (c == -1)
            break;

        switch (c)
        {

        case 'm':
            min_size = atoi(optarg);
            break;

        case 'p':
            show_progress = true;
            break;

        case 't':
            omp_set_num_threads(atoi(optarg));
            break;

        case 'h':
        case '?':
            /* getopt_long already printed an error message. */
            help_simplify(argv);
            exit(1);
            break;

        default:
            abort ();

        }
    }

    // NB: If you worry about "still reachable but possibly lost" warnings in valgrind,
    // this would free all the memory used by protobuf:
    //ShutdownProtobufLibrary();

    return 0;
}

// Register subcommand
static Subcommand vg_construct("simplify", "graph simplification", main_simplify);

