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
#include <vg/io/stream.hpp>
#include <vg/io/vpkg.hpp>
#include <handlegraph/mutable_path_mutable_handle_graph.hpp>
#include "bdsg/packed_graph.hpp"
#include "bdsg/hash_graph.hpp"
#include <bdsg/overlays/overlay_helper.hpp>
#include "../io/save_handle_graph.hpp"
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
        unique_ptr<MutablePathMutableHandleGraph> graph;
        string graph_filename = get_input_file_name(optind, argc, argv);
        graph = vg::io::VPKG::load_one<MutablePathMutableHandleGraph>(graph_filename);            
            
        if (sort || compact) {
            // We need to reassign IDs
            hash_map<nid_t, nid_t> new_ids;
            
            if (compact && !sort) {
                // We are compacting, but do not need to topologically sort
                
                // Loop over all the nodes in the graph's order and assign them new IDs in ID order.
                // This is slower than it needs to be, but gets us nice results even on graphs that don't preserve node order.
                // TODO: counting all the nodes may be an O(1) scan of the graph to save some vector copies.
                vector<nid_t> all_ids;
                all_ids.reserve(graph->get_node_count());
                graph->for_each_handle([&](const handle_t& h) {
                    all_ids.emplace_back(graph->get_id(h));
                });
                std::sort(all_ids.begin(), all_ids.end());
                
                // Now invert the vector's mapping
                new_ids.reserve(all_ids.size());
                for (nid_t i = 1; i < all_ids.size() + 1; i++) {
                    new_ids[all_ids[i - 1]] = i;
                }
            } else {
                // We are sorting to assign IDs, which inherently compacts.
                
                // We only need to sort the ID numbers, not the graph's iteration order (if any).
                auto handle_order = handlealgs::topological_order(graph.get());
                
                // Now invert the order's mapping
                new_ids.reserve(handle_order.size());
                for (nid_t i = 1; i < handle_order.size() + 1; i++) {
                    new_ids[graph->get_id(handle_order[i - 1])] = i;
                }
            }
            
            // Now assign the IDs. If we find any e.g. dangling paths or
            // edges with no nodes we will crash.
            graph->reassign_node_ids([&](const nid_t& old_id) {
                return new_ids.at(old_id);
            });
        }

        if (increment != 0) {
            graph->increment_node_ids(increment);
        }

        if (decrement != 0) {
            graph->increment_node_ids(-increment);
        }

        vg::io::save_handle_graph(graph.get(), cout);
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

