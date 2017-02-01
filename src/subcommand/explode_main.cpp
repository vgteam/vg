/**
 * \file explode_main.cpp: break a graph into connected components
 */

#include <omp.h>
#include <unistd.h>
#include <getopt.h>

#include "subcommand.hpp"

#include "../vg.hpp"
#include "../stream.hpp"
#include "../utility.hpp"


using namespace std;
using namespace vg;
using namespace vg::subcommand;

void help_explode(char** argv) {
    cerr << "usage: " << argv[0] << " explode [options] source.vg part_dir" << endl
         << "Breaks a graph into connected components in their own files in the given directory" << endl
         << endl
         << "options:" << endl
         << "general:" << endl
         << "    -t, --threads N          for tasks that can be done in parallel, use this many threads [1]" << endl
         << "    -h, --help" << endl;
}

int main_explode(int argc, char** argv) {

    if (argc == 2) {
        help_explode(argv);
        return 1;
    }

    int threads = 1;
    
    int c;
    optind = 2; // force optind past command positional argument
    while (true) {
        static struct option long_options[] =
        {
            {"help", no_argument, 0, 'h'},
            {"threads", required_argument, 0, 't'},
            {0, 0, 0, 0}
        };

        int option_index = 0;
        c = getopt_long (argc, argv, "ht:",
                long_options, &option_index);


        // Detect the end of the options.
        if (c == -1)
            break;

        switch (c)
        {

        case 't':
            omp_set_num_threads(atoi(optarg));
            break;

        case 'h':
        case '?':
        default:
            help_explode(argv);
            exit(1);
            break;
        }
    }

    VG* graph;
    get_input_file(optind, argc, argv, [&](istream& in) {
        graph = new VG(in);
    });

    // Grab the directory name to put stuff in
    string output_dir = get_input_file_name(optind, argc, argv);
    
    // Make sure it exists if it's not a directory already
    mkdir(output_dir.c_str(), 0644);
    // Ignore failure
    
    // Now we explode the VG
    
    // Count through the components we build
    size_t component_index = 0;
    
    // Track all the nodes we've already assigned to subgraphs
    set<Node*> used;
    
    graph->for_each_node([&](Node* start) {
        if (!used.count(start)) {
            // It's a new connected component!
            VG component;
            graph->for_each_connected_node(start, [&](Node* n) {
                // Mark this connected node as used in a component.
                used.insert(n);
                
                // Copy node over
                component.create_node(n->sequence(), n->id());
                
                // Copy over its edges
                for (auto* e : graph->edges_of(n)) {
                    component.add_edge(*e);
                }
                
                // Copy paths over
                for (auto& path : graph->paths.get_node_mapping(n)) {
                    for (auto& m : path.second) {
                        component.paths.append_mapping(path.first, *m);
                    }
                }
            });
            
            component.sync_paths();
            
            // Save the component
            string filename = output_dir + "/component" + to_string(component_index) + ".vg";
            component.serialize_to_file(filename);
            
            component_index++;
        }
    });
    
    if (graph != nullptr) {
        delete graph;
    }
    
    return 0;
}

// Register subcommand
static Subcommand vg_explode("explode", "split graph into connected components", main_explode);


