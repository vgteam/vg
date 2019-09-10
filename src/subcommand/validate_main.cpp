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
#include "../xg.hpp"
#include "../alignment.hpp"
#include <vg/io/vpkg.hpp>

using namespace std;
using namespace vg;
using namespace vg::subcommand;

void help_validate(char** argv) {
    cerr << "usage: " << argv[0] << " validate [options] [graph]" << endl
         << "Validate the graph." << endl
         << endl
         << "options:" << endl
         << "    default: check all aspects of the graph, if options are specified do only those" << endl
         << "    -n, --nodes     verify that we have the expected number of nodes" << endl
         << "    -e, --edges     verify that the graph contains all nodes that are referred to by edges" << endl
         << "    -p, --paths     verify that contiguous path segments are connected by edges" << endl
         << "    -o, --orphans   verify that all nodes have edges" << endl
         << "    -a, --gam FILE  verify that edits in the alignment fit on nodes in the graph" << endl
         << "    -x, --xg FILE   index to use for -a" << endl;
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
    string xg_path;
    string gam_path;

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
            {"gam", required_argument, 0, 'a'},
            {"xg", required_argument, 0, 'x'},
            {0, 0, 0, 0}
        };

        int option_index = 0;
        c = getopt_long (argc, argv, "hneopa:x:",
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

            case 'a':
                gam_path = optarg;
                break;

            case 'x':
                xg_path= optarg;
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

    if (!gam_path.empty() || !xg_path.empty()) {
        // GAM validation is its entirely own thing
        if (xg_path.empty()) {
            cerr << "error:[vg validate] xg index (-x) required with (-a)" << endl;
            return 1;
        } else if (gam_path.empty()) {
            cerr << "error:[vg validate] gam alignment (-a) required with (-x)" << endl;
            return 1;
        } else if (check_nodes || check_edges || check_orphans || check_paths) {
            cerr << "error:[vg validate] -n, -e -o, -p cannot be used with -a and -x" << endl;
            return 1;
        }
        ifstream in(xg_path.c_str());
        unique_ptr<PathPositionHandleGraph> xindex = vg::io::VPKG::load_one<PathPositionHandleGraph>(in);
        in.close();
        get_input_file(gam_path, [&](istream& in) {
                vg::io::for_each<Alignment>(in, [&](Alignment& aln) {
                        if (!alignment_is_valid(aln, xindex.get())) {
                            exit(1);
                        }
                    });
            });
        return 0;
    } else {

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
}

// Register subcommand
static Subcommand vg_validate("validate", "validate the semantics of a graph or gam", DEVELOPMENT, main_validate);

