/** \file validate_main.cpp
 *
 * Defines the "vg validate" subcommand
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

void help_validate(char** argv) {
    cerr << "usage: " << argv[0] << " validate [options] graph" << endl
        << "Validate the graph." << endl
        << endl
        << "options:" << endl
        << "    default: check all aspects of the graph, if options are specified do only those" << endl
        << "    -n, --nodes    verify that we have the expected number of nodes" << endl
        << "    -e, --edges    verify that the graph contains all nodes that are referred to by edges" << endl
        << "    -p, --paths    verify that contiguous path segments are connected by edges" << endl
        << "    -o, --orphans  verify that all nodes have edges" << endl;
}

int main_validate(int argc, char** argv) {

    if (argc <= 2) {
        help_validate(argv);
        return 1;
    }

    bool check_nodes = false;
    bool check_edges = false;
    bool check_orphans = false;
    bool check_paths = false;

    int c;
    optind = 2; // force optind past command positional argument
    while (true) {
        static struct option long_options[] =
        {
            {"help", no_argument, 0, 'h'},
            {"nodes", no_argument, 0, 'n'},
            {"edges", no_argument, 0, 'e'},
            {"paths", no_argument, 0, 'o'},
            {"orphans", no_argument, 0, 'p'},
            {0, 0, 0, 0}
        };

        int option_index = 0;
        c = getopt_long (argc, argv, "hneop",
                long_options, &option_index);

        // Detect the end of the options.
        if (c == -1)
            break;

        switch (c)
        {

            case 'n':
                check_nodes = true;
                break;

            case 'e':
                check_edges = true;
                break;

            case 'o':
                check_orphans = true;
                break;

            case 'p':
                check_paths = true;
                break;

            case 'h':
            case '?':
                help_validate(argv);
                exit(1);
                break;

            default:
                abort ();
        }
    }

    VG* graph;
    get_input_file(optind, argc, argv, [&](istream& in) {
        graph = new VG(in);
    });

    // if we chose a specific subset, do just them
    if (check_nodes || check_edges || check_orphans || check_paths) {
        if (graph->is_valid(check_nodes, check_edges, check_orphans, check_paths)) {
            return 0;
        } else {
            return 1;
        }
        // otherwise do everything
    } else if (graph->is_valid()) {
        return 0;
    } else {
        return 1;
    }
}

// Register subcommand
static Subcommand vg_validate("validate", "validate the semantics of a graph", main_validate);

