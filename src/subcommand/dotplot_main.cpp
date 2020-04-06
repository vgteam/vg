/** \file dotplot_main.cpp
 *
 * Defines the "vg dotplot" subcommand, which renders dotplots from an xg index relating the embedded paths
 */


#include <omp.h>
#include <unistd.h>
#include <getopt.h>

#include <iostream>

#include "subcommand.hpp"

#include "../vg.hpp"
#include "../xg.hpp"
#include <vg/io/vpkg.hpp>
#include <bdsg/overlays/overlay_helper.hpp>
#include "../position.hpp"

#include "algorithms/nearest_offsets_in_paths.hpp"

using namespace std;
using namespace vg;
using namespace vg::subcommand;

void help_dotplot(char** argv) {
    cerr << "usage: " << argv[0] << " dotplot [options]" << endl
         << "options:" << endl
         << "  input:" << endl
         << "    -x, --xg FILE         use the graph or the XG index FILE" << endl;
    //<< "  output:" << endl;
}

int main_dotplot(int argc, char** argv) {

    if (argc == 2) {
        help_dotplot(argv);
        return 1;
    }

    string xg_file;

    int c;
    optind = 2; // force optind past command positional argument
    while (true) {
        static struct option long_options[] =
        {
            {"xg", required_argument, 0, 'x'},
            {0, 0, 0, 0}
        };

        int option_index = 0;
        c = getopt_long (argc, argv, "hx:",
                long_options, &option_index);

        // Detect the end of the options.
        if (c == -1)
            break;

        switch (c)
        {

        case 'x':
            xg_file = optarg;
            break;

        case 'h':
        case '?':
            help_dotplot(argv);
            exit(1);
            break;

        default:
            abort ();
        }
    }

    if (xg_file.empty()) {
        cerr << "[vg dotplot] Error: an xg index is required" << endl;
        exit(1);
    } else {
        unique_ptr<PathHandleGraph> path_handle_graph = vg::io::VPKG::load_one<PathHandleGraph>(xg_file);
        bdsg::PathPositionOverlayHelper overlay_helper;
        PathPositionHandleGraph* xindex = overlay_helper.apply(path_handle_graph.get());

        cout << "query.name" << "\t"
             << "query.pos" << "\t"
             << "orientation" << "\t"
             << "target.name" << "\t"
             << "target.pos" << endl;
        xindex->for_each_handle([&](const handle_t& h) {
                vg::id_t id = xindex->get_id(h);
                for (size_t i = 0; i < xindex->get_length(xindex->get_handle(id)); ++i) {
                    pos_t p = make_pos_t(id, false, i);
                    map<string, vector<pair<size_t, bool> > > offsets = algorithms::offsets_in_paths(&(*xindex), p);
                    // cross the offsets in output
                    for (auto& o : offsets) {
                        auto& name1 = o.first;
                        for (auto& t : o.second) {
                            for (auto& m : offsets) {
                                auto& name2 = m.first;
                                for (auto& w : m.second) {
                                    cout << name1 << "\t"
                                         << t.first << "\t"
                                         << (t.second == w.second ? "+" : "-") << "\t"
                                         << name2 << "\t"
                                         << w.first << endl;
                                }
                            }
                        }
                    }
                }
            });
    }
    
    return 0;

}

// Register subcommand
static Subcommand vg_dotplot("dotplot", "generate the dotplot matrix from the embedded paths in an xg index", main_dotplot);

