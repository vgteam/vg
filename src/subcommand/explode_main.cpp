/**
 * \file explode_main.cpp: break a graph into connected components
 */

#include <omp.h>
#include <unistd.h>
#include <getopt.h>

#include "subcommand.hpp"

#include "../vg.hpp"
#include <vg/io/stream.hpp>
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
            omp_set_num_threads(parse<int>(optarg));
            break;

        case 'h':
        case '?':
        default:
            help_explode(argv);
            exit(1);
            break;
        }
    }

    cerr << "vg explode is deprecated.  Please use \"vg chunk -C source.vg -b part_dir/component\" for same* functionality as \"vg explode source.vg part_dir\"" << endl
         << " * (unlike explode, the output directory must already exist when running chunk, though)" << endl;
    return 1;

    VG* graph;
    get_input_file(optind, argc, argv, [&](istream& in) {
        graph = new VG(in);
    });

    // Grab the directory name to put stuff in
    string output_dir = get_output_file_name(optind, argc, argv);
    
    // Make sure it exists if it's not a directory already
    mkdir(output_dir.c_str(), 0755);
    // Ignore failure
    
    // Now we explode the VG
    
    // Count through the components we build
    size_t component_index = 0;
    
    // Track all the nodes we've already assigned to subgraphs
    unordered_set<handle_t> used;
    
    graph->for_each_handle([&](const handle_t& start) {
        if (!used.count(start)) {
            // It's a new connected component!
            // TODO: this could be replaced by any handle graph type now
            VG component;
            
            // We want to track the path names in each component
            set<path_handle_t> paths;
            
            deque<handle_t> queue{start};
            
            // Mark this connected node as used in a component.
            used.insert(start);
            
            while (!queue.empty()) {
                
                handle_t handle = queue.front();
                queue.pop_front();
                
                // Copy node over
                handle_t new_handle = component.create_handle(graph->get_sequence(handle),
                                                              graph->get_id(handle));
                
                // Copy over its edges and queue the next handles
                graph->follow_edges(handle, false, [&](const handle_t& next) {
                    if (component.has_node(graph->get_id(next))) {
                        component.create_edge(new_handle, component.get_handle(graph->get_id(next),
                                                                               graph->get_is_reverse(next)));
                    }
                    if (!used.count(next)) {
                        queue.push_back(next);
                        used.insert(next);
                    }
                });
                graph->follow_edges(handle, true, [&](const handle_t& prev) {
                    if (component.has_node(graph->get_id(prev))) {
                        component.create_edge(component.get_handle(graph->get_id(prev),
                                                                   graph->get_is_reverse(prev)), new_handle);
                    }
                    if (!used.count(prev)) {
                        queue.push_back(prev);
                        used.insert(prev);
                    }
                });
                
                // Record paths
                graph->for_each_step_on_handle(handle, [&](const step_handle_t& step) {
                    paths.insert(graph->get_path_handle_of_step(step));
                });
            }
            
            // Copy the paths over
            for (path_handle_t path_handle : paths) {
                path_handle_t new_path_handle = component.create_path_handle(graph->get_path_name(path_handle),
                                                                             graph->get_is_circular(path_handle));
                for (handle_t handle : graph->scan_path(path_handle)) {
                    component.append_step(new_path_handle, component.get_handle(graph->get_id(handle),
                                                                                graph->get_is_reverse(handle)));
                }
            }
            
            // Save the component
            string filename = output_dir + "/component" + to_string(component_index) + ".vg";
            
            // Now report what paths went into the component in parseable TSV
            cout << filename;
            for (auto& path_handle : paths) {
                cout << "\t" << graph->get_path_name(path_handle);
            }
            cout << endl;
            
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
static Subcommand vg_explode("explode", "split graph into connected components", DEPRECATED, main_explode);


