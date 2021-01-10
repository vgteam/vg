/** \file join_main.cpp
 *
 * Defines the "vg join" subcommand, which attaches graphs together.
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

void help_join(char** argv) {
    cerr << "usage: " << argv[0] << " join [options] <graph1.vg> [graph2.vg ...] >joined.vg" << endl
        << "Joins graphs and sub-graphs into a single variant graph by connecting their" << endl
        << "heads to a single root node with sequence 'N'." << endl
        << "Assumes a single id namespace for all graphs to join." << endl;
}

int main_join(int argc, char** argv) {

    if (argc == 2) {
        help_join(argv);
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
                help_join(argv);
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

    VG joined;
    for (list<VG*>::iterator g = graphs.begin(); g != graphs.end(); ++g) {
        // Stick all the graphs together, complaining if they use the same node IDs (since they probably shouldn't).
        joined.extend(**g, true);
    }

    // combine all subgraphs
    joined.join_heads();

    // output
    joined.serialize_to_ostream(std::cout);

    return 0;
}

// Register subcommand
static Subcommand vg_join("join", "combine graphs via a new head", DEPRECATED, main_join);

