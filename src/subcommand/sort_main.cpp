/** \file sort_main.cpp
 *
 * Defines the "vg sort" subcommand, which sorts graph nodes.
 */


#include <omp.h>
#include <unistd.h>
#include <getopt.h>

#include <iostream>

#include "subcommand.hpp"

#include "../vg.hpp"
#include "../stream_index.hpp"
#include "../flow_sort.hpp"
#include "../algorithms/gfa_to_handle.hpp"
#include "../algorithms/id_sort.hpp"
#include <vg/io/vpkg.hpp>

using namespace std;
using namespace vg;
using namespace vg::subcommand;

void help_sort(char** argv){
    cerr << "usage: " << argv[0] << " sort [options] > sorted.vg " << endl
         << "options: " << endl
         << "    -a, --algorithm NAME   sort by the given algorithm (eades, max-flow, id, or topo; default id)" << endl
         << "    -g, --gfa              input in GFA format" << endl
         << "    -r, --ref              reference name, for eades and max-flow algorithms; makes -a default to max-flow" << endl
         << "    -w, --without-grooming no grooming mode for eades" << endl
         << "    -I, --index-to FILE    produce an index of an id-sorted vg file to the given filename" << endl
         << endl;
}

int main_sort(int argc, char *argv[]) {

    // What should we sort the graph by?
    string algorithm;
    
    // Default input format is vg, but we can also read GFA
    bool gfa_input = false;
    
    string reference_name;
    bool without_grooming = false;
    string sorted_index_filename;
    
    int c;
    optind = 2; // force optind past command positional argument
    while (true) {
        static struct option long_options[] =
            {
                {"algorithm", required_argument, 0, 'a'},
                {"gfa", no_argument, 0, 'g'},
                {"ref", required_argument, 0, 'r'},
                {"without-grooming", no_argument, 0, 'w'},
                {"index-to", no_argument, 0, 'I'},
                {0, 0, 0, 0}
            };

        int option_index = 0;
        c = getopt_long (argc, argv, "a:gr:wI:",
                         long_options, &option_index);

        /* Detect the end of the options. */
        if (c == -1)
            break;

        switch (c)
        {
        case 'a':
            algorithm = optarg;
            break;
        case 'g':
            gfa_input = true;
            break;
        case 'r':
            reference_name = optarg;
            if (algorithm.empty()) {
                algorithm = "max-flow";
            }
            break;
        case 'w':
            without_grooming = true;
            break;
        case 'I':
            sorted_index_filename = optarg;
            break;
        case 'h':
        case '?':
            /* getopt_long already printed an error message. */
            help_sort(argv);
            exit(1);
            break;
        default:
            abort();
        }
    }
    
    if (algorithm.empty()) {
        // Set the default algorithm
        algorithm = "id";
    }
    
    // Validate the algorithm selection and option combination
    if (algorithm == "id" || algorithm == "topo") {
        if (!reference_name.empty()) {
            cerr << "error[vg sort]: Reference name not used with " << algorithm << " sort algorithm" << endl;
            exit(1);
        }
        if (without_grooming) {
            cerr << "error[vg sort]: Not sensible to turn off grooming with " << algorithm << " sort algorithm" << endl;
            exit(1);
        }
    } else if (algorithm == "max-flow" || algorithm == "eades") {
        if (reference_name.empty()) {
            cerr << "error[vg sort]: Reference name required with " << algorithm << " sort algorithm" << endl;
            exit(1);
        }
    } else {
        cerr << "error[vg sort]: Unrecognized sort algorithm " << algorithm << endl;
        exit(1);
    }
    if (!sorted_index_filename.empty() && algorithm != "id") {
        cerr << "error[vg sort]: Sorted VG index can only be produced when sorting by ID" << endl;
        exit(1);
    }
    
    // With the input graph file
    string filename = get_input_file_name(optind, argc, argv);
    
    // We will load it into this graph
    std::unique_ptr<VG> graph;
    
    if (gfa_input) {
        // Read as GFA
        graph.reset(new VG());
        try {
            algorithms::gfa_to_path_handle_graph(filename, graph.get());
        } catch(algorithms::GFAFormatError& e) {
            // GFA loading has failed because the file is invalid
            cerr << e.what() << endl;
            exit(1);
        } catch(ios_base::failure& e) {
            // GFA loading has failed because the file couldn't be read
            cerr << e.what() << endl;
            exit(1);
        }
    } else {
        // Read as Handle Graph and copy into VG           
        unique_ptr<handlegraph::MutablePathMutableHandleGraph> handle_graph =
            vg::io::VPKG::load_one<handlegraph::MutablePathMutableHandleGraph>(filename);
        VG* vg_graph = dynamic_cast<vg::VG*>(graph.get());
        if (vg_graph != nullptr) {
            graph.reset(vg_graph);
            handle_graph.release();
        } else {
            // Copy instead.
            vg_graph = new vg::VG();
            handlealgs::copy_path_handle_graph(handle_graph.get(), vg_graph);
            // Give the unique_ptr ownership and delete the graph we loaded.
            graph.reset(vg_graph);
            // Make sure the paths are all synced up
            vg_graph->paths.to_graph(vg_graph->graph);
        }
    }
    
    // Now sort the graph
    
    if (algorithm == "max-flow" || algorithm == "eades") {
        // Do max flow sort or Eades algorithm sort
        FlowSort flow_sort(*graph.get());
        if (algorithm == "eades") {
            // Use Eades algorithm
            flow_sort.fast_linear_sort(reference_name, !without_grooming);
        } else {
            // Use max flow
            flow_sort.max_flow_sort(reference_name);
        }
    } else if (algorithm == "id") {
        // Sort by ID
        graph.get()->id_sort();
    } else if (algorithm == "topo") {
        // Sort topologically
        graph.get()->sort();
    } else {
        throw runtime_error("Unimplemented sort algorithm: " + algorithm);
    }
    
    // We have an optional index, which will outlive our emitter
    unique_ptr<StreamIndex<Graph>> index;
    if (!sorted_index_filename.empty()) {
        // Make an index we can use later for graph random access
        index = unique_ptr<StreamIndex<Graph>>(new StreamIndex<Graph>());
    }
    
    // Maintain our own group buffer at a higher scope than the emitter.
    vector<Graph> group_buffer;
    
    {
        // Make our own emitter for serialization
        vg::io::ProtobufEmitter<Graph> emitter(std::cout);
        
        if (index) {
            emitter.on_message([&](const Graph& g) {
                // Copy every graph that is emitted.
                // TODO: Just compute indexing stats instead.
                group_buffer.push_back(g);
            });
        
            emitter.on_group([&](int64_t start_vo, int64_t past_end_vo) {
                // On every group, tell the index to record the group stats, and clear the buffer.
                index->add_group(group_buffer, start_vo, past_end_vo);
                group_buffer.clear();
            });
        }
        
        // Save the sorted graph to the emitter
        graph->serialize_to_emitter(emitter);
    }
    
    if (index) {
        // Now save out the index
        ofstream index_out(sorted_index_filename);
        index->save(index_out);
    }
        
    return 0;
}

// Register subcommand
static Subcommand vg_sort("sort", "sort variant graph by various algorithms", DEPRECATED, main_sort);

