/** \file concat_main.cpp
 *
 * Defines the "vg concat" subcommand
 */


#include <omp.h>
#include <unistd.h>
#include <getopt.h>

#include <iostream>

#include "subcommand.hpp"

#include "../vg.hpp"

using namespace std;
using namespace vg;
using namespace vg::subcommand;

void help_concat(char** argv) {
    cerr << "usage: " << argv[0] << " concat [options] <graph1.vg> [graph2.vg ...] >merged.vg" << endl
        << "Concatenates graphs in order by adding edges from the tail nodes of the" << endl
        << "predecessor to the head nodes of the following graph. Node IDs are" << endl
        << "compacted, so care should be taken if consistent IDs are required." << endl;
}

int main_concat(int argc, char** argv) {

    if (argc == 2) {
        help_concat(argv);
        return 1;
    }

    int c;
    optind = 2; // force optind past command positional argument
    while (true) {
        static struct option long_options[] =
        {
            {"help", no_argument, 0, 'h'},
            {0, 0, 0, 0}
        };

        int option_index = 0;
        c = getopt_long (argc, argv, "h",
                long_options, &option_index);

        // Detect the end of the options.
        if (c == -1)
            break;

        switch (c)
        {
            case 'h':
            case '?':
                help_concat(argv);
                exit(1);
                break;

            default:
                abort ();
        }
    }

    list<VG*> graphs;

    while (optind < argc) {
        VG* graph;
        get_input_file(optind, argc, argv, [&](istream& in) {
            graph = new VG(in);
        });
        graphs.push_back(graph);
    }

    VG merged;
    for (list<VG*>::iterator g = graphs.begin(); g != graphs.end(); ++g) {
        merged.append(**g);
    }

    // output
    merged.serialize_to_ostream(std::cout);

    return 0;
}

// Register subcommand
static Subcommand vg_concat("concat", "concatenate graphs tail-to-head", main_concat);

