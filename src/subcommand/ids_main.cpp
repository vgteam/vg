/** \file ids_main.cpp
 *
 * Defines the "vg ids" subcommand, which modifies node IDs.
 */


#include <omp.h>
#include <unistd.h>
#include <getopt.h>

#include <iostream>

#include "subcommand.hpp"

#include "../vg.hpp"
#include "../vg_set.hpp"
#include "../algorithms/topological_sort.hpp"

#include <gcsa/support.h>

using namespace std;
using namespace vg;
using namespace vg::subcommand;

void help_ids(char** argv) {
    cerr << "usage: " << argv[0] << " ids [options] <graph1.vg> [graph2.vg ...] >new.vg" << endl
        << "options:" << endl
        << "    -c, --compact        minimize the space of integers used by the ids" << endl
        << "    -i, --increment N    increase ids by N" << endl
        << "    -d, --decrement N    decrease ids by N" << endl
        << "    -j, --join           make a joint id space for all the graphs that are supplied" << endl
        << "                         by iterating through the supplied graphs and incrementing" << endl
        << "                         their ids to be non-conflicting (modifies original files)" << endl
        << "    -m, --mapping FILE   create an empty node mapping for vg prune" << endl
        << "    -s, --sort           assign new node IDs in (generalized) topological sort order" << endl;
}

int main_ids(int argc, char** argv) {

    if (argc == 2) {
        help_ids(argv);
        return 1;
    }

    bool join = false;
    bool compact = false;
    bool sort = false;
    int64_t increment = 0;
    int64_t decrement = 0;
    std::string mapping_name;

    int c;
    optind = 2; // force optind past command positional argument
    while (true) {
        static struct option long_options[] =
        {
            {"compact", no_argument, 0, 'c'},
            {"increment", required_argument, 0, 'i'},
            {"decrement", required_argument, 0, 'd'},
            {"join", no_argument, 0, 'j'},
            {"mapping", required_argument, 0, 'm'},
            {"sort", no_argument, 0, 's'},
            {"help", no_argument, 0, 'h'},
            {0, 0, 0, 0}
        };

        int option_index = 0;
        c = getopt_long (argc, argv, "hci:d:jm:s",
                long_options, &option_index);

        // Detect the end of the options.
        if (c == -1)
            break;

        switch (c)
        {
            case 'c':
                compact = true;
                break;

            case 'i':
                increment = parse<int>(optarg);
                break;

            case 'd':
                decrement = parse<int>(optarg);
                break;

            case 'j':
                join = true;
                break;

            case 'm':
                mapping_name = optarg;
                break;

            case 's':
                sort = true;
                break;

            case 'h':
            case '?':
                help_ids(argv);
                exit(1);
                break;

            default:
                abort ();
        }
    }

    if (!join && mapping_name.empty()) {
        VG* graph;
        get_input_file(optind, argc, argv, [&](istream& in) {
            graph = new VG(in);
        });

        if (sort) {
            // Set up the nodes so we go through them in topological order
            algorithms::sort(graph);
        }

        if (compact || sort) {
            // Compact only, or compact to re-assign IDs after sort
            graph->compact_ids();
        }

        if (increment != 0) {
            graph->increment_node_ids(increment);
        }

        if (decrement != 0) {
            graph->decrement_node_ids(decrement);
        }

        graph->serialize_to_ostream(std::cout);
        delete graph;
    } else {

        vector<string> graph_file_names;
        while (optind < argc) {
            string file_name = get_input_file_name(optind, argc, argv);
            graph_file_names.push_back(file_name);
        }

        VGset graphs(graph_file_names);
        vg::id_t max_node_id = (join ? graphs.merge_id_space() : graphs.max_node_id());
        if (!mapping_name.empty()) {
            gcsa::NodeMapping mapping(max_node_id + 1);
            std::ofstream out(mapping_name, std::ios_base::binary);
            if (!out) {
                std::cerr << "[vg ids]: cannot create node mapping file " << mapping_name << std::endl;
            } else {
                mapping.serialize(out);
                out.close();
            }
        }
    }

    return 0;

}

// Register subcommand
static Subcommand vg_ids("ids", "manipulate node ids", TOOLKIT, main_ids);

