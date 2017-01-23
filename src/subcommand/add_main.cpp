/** \file add_main.cpp
 *
 * Defines the "vg add" subcommand, which adds in variation from a VCF to an
 * existing graph.
 */

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

void help_add(char** argv) {
    cerr << "usage: " << argv[0] << " add [options] old.vg >new.vg" << endl
         << "options:" << endl
         << "    -v, --vcf FILE         add in variants from the given VCF file" << endl
         << "    -p, --progress         show progress" << endl
         << "    -t, --threads N        use N threads (defaults to numCPUs)" << endl;
}

int main_add(int argc, char** argv) {

    if (argc == 2) {
        help_add(argv);
        return 1;
    }

    
    // We can have one or more VCFs
    vector<string> vcf_filenames;
    bool show_progress = false;

    int c;
    optind = 2; // force optind past command positional argument
    while (true) {
        static struct option long_options[] =
            {
                {"vcf", required_argument, 0, 'v'},
                {"progress",  no_argument, 0, 'p'},
                {"threads", required_argument, 0, 't'},
                {"help", no_argument, 0, 'h'},
                {0, 0, 0, 0}
            };

        int option_index = 0;
        c = getopt_long (argc, argv, "v:pt:h?",
                         long_options, &option_index);

        /* Detect the end of the options. */
        if (c == -1)
            break;

        switch (c)
        {

        case 'v':
            vcf_filenames.push_back(optarg);

        case 'p':
            show_progress = true;
            break;
            
        case 't':
            omp_set_num_threads(atoi(optarg));
            break;

        case 'h':
        case '?':
            /* getopt_long already printed an error message. */
            help_add(argv);
            exit(1);
            break;

        default:
            abort ();

        }
    }
    
    // Load the graph
    VG* graph;
    get_input_file(optind, argc, argv, [&](istream& in) {
        graph = new VG(in, show_progress);
    });
    
    if (graph == nullptr) {
        cerr << "error:[vg add]: Could not load graph" << endl;
        exit(1);
    }
    
    // TODO: all the actual work
    
    delete graph;

    // NB: If you worry about "still reachable but possibly lost" warnings in valgrind,
    // this would free all the memory used by protobuf:
    //ShutdownProtobufLibrary();

    return 0;
}

// Register subcommand
static Subcommand vg_add("add", "add variants from a VCF to a graph", main_add);

