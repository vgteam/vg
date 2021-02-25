/** \file concat_main.cpp
 *
 * Defines the "vg concat" subcommand
 */


#include <omp.h>
#include <unistd.h>
#include <getopt.h>

#include <iostream>

#include "subcommand.hpp"
#include "../option.hpp"
#include "../xg.hpp"
#include "../vg.hpp"
#include <vg/io/stream.hpp>
#include <vg/io/vpkg.hpp>
#include "../io/save_handle_graph.hpp"
#include <handlegraph/mutable_path_mutable_handle_graph.hpp>

using namespace std;
using namespace vg;
using namespace vg::subcommand;

void help_concat(char** argv) {
    cerr << "usage: " << argv[0] << " concat [options] <graph1.vg> [graph2.vg ...] >merged.vg" << endl
         << "Concatenates graphs in order by adding edges from the tail nodes of the" << endl
         << "predecessor to the head nodes of the following graph.  If node ID spaces overlap "
         << "between graphs, they will be resolved (as in vg ids -j)" << endl
         << endl
         << "Options:" << endl
         << "    -p, --only-join-paths         Only add edges necessary to join up appended paths (as opposed between all heads/tails)" << endl
         << endl;
}

int main_concat(int argc, char** argv) {

    if (argc == 2) {
        help_concat(argv);
        return 1;
    }

    bool only_join_paths = false;
    
    int c;
    optind = 2; // force optind past command positional argument
    while (true) {
        static struct option long_options[] =
        {
            {"help", no_argument, 0, 'h'},
            {"only-join-paths", no_argument, 0, 'p'},
            {0, 0, 0, 0}
        };

        int option_index = 0;
        c = getopt_long (argc, argv, "hp",
                long_options, &option_index);

        // Detect the end of the options.
        if (c == -1)
            break;

        switch (c)
        {
        case 'p':
            only_join_paths = true;
            break;
        case 'h':
        case '?':
            help_concat(argv);
            exit(1);
            break;
        default:
            abort ();
        }
    }

    unique_ptr<MutablePathMutableHandleGraph> first_graph;
    get_input_file(optind, argc, argv, [&](istream& in) {
            first_graph = vg::io::VPKG::load_one<MutablePathMutableHandleGraph>(in);
        });
    int64_t max_node_id = first_graph->max_node_id();

    while (optind < argc) {

        unique_ptr<MutablePathMutableHandleGraph> graph;
        get_input_file(optind, argc, argv, [&](istream& in) {
                graph = vg::io::VPKG::load_one<MutablePathMutableHandleGraph>(in);
            });

        // join the id spaces if necessary
        int64_t delta = max_node_id - graph->min_node_id();
        if (delta >= 0) {
            graph->increment_node_ids(delta + 1);
        }
        max_node_id = graph->max_node_id();

        handlealgs::append_path_handle_graph(graph.get(), first_graph.get(), only_join_paths);
    }

    // Serialize the graph using VPKG.
    vg::io::save_handle_graph(first_graph.get(), cout);

    return 0;
}

// Register subcommand
static Subcommand vg_concat("concat", "concatenate graphs tail-to-head", DEPRECATED, main_concat);

