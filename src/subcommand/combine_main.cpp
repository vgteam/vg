/** \file combine_main.cpp
 *
 * Defines the "vg combine" subcommand
 */


#include <omp.h>
#include <unistd.h>
#include <getopt.h>

#include <iostream>

#include <vg/io/vpkg.hpp>
#include <vg/io/message_iterator.hpp>
#include <vg/io/blocked_gzip_input_stream.hpp>

#include "subcommand.hpp"

#include "../handle.hpp"
#include "../vg.hpp"
#include "../io/save_handle_graph.hpp"

using namespace std;
using namespace vg;
using namespace vg::subcommand;

void help_combine(char** argv) {
    cerr << "usage: " << argv[0] << " combine [options] <graph1.vg> [graph2.vg ...] >merged.vg" << endl
         << "Combines one or more graphs into a single file, regardless of input format." << endl
         << "Node IDs will be modified as needed to resolve conflicts (in same manner as vg ids -j)." << endl
         << endl
         << "Options:" << endl
         << "    -c, --cat-proto       Merge graphs by converting each to Protobuf (if not already) and catting the results."
         << "                          Node IDs not modified [DEPRECATED]" << endl
         << "    -p, --connect-paths   Add edges necessary to connect paths with the same name present in different graphs." << endl
         << "                          ex: If path x is present in graphs N-1 and N, then an edge connecting the last node of x in N-1 " << endl
         << "                          and the first node of x in N will be added." << endl;
}

static int cat_proto_graphs(int argc, char** argv);

int main_combine(int argc, char** argv) {

    if (argc == 2) {
        help_combine(argv);
        return 1;
    }

    bool connect_paths = false;
    bool cat_proto = false;
    
    int c;
    optind = 2; // force optind past command positional argument
    while (true) {
        static struct option long_options[] =
        {
            {"help", no_argument, 0, 'h'},
            {"connect-paths", no_argument, 0, 'p'},
            {"cat-proto", no_argument, 0, 'c'},
            {0, 0, 0, 0}
        };

        int option_index = 0;
        c = getopt_long (argc, argv, "hpc",
                long_options, &option_index);

        // Detect the end of the options.
        if (c == -1)
            break;

        switch (c)
        {
        case 'p':
            connect_paths = true;
            break;
        case 'c':
            cat_proto = true;
            break;
        case 'h':
        case '?':
            help_combine(argv);
            exit(1);
            break;

        default:
            abort ();
        }
    }

    if (cat_proto) {
        if (connect_paths) 
        cerr << "warning [vg combine]: --cat-proto/-c option is deprecated and will be removed in a future version of vg." << endl;
        return cat_proto_graphs(argc, argv);
    }
    
    unique_ptr<MutablePathMutableHandleGraph> first_graph;
    string first_graph_filename = get_input_file_name(optind, argc, argv);
    first_graph = vg::io::VPKG::load_one<MutablePathMutableHandleGraph>(first_graph_filename);
    int64_t max_node_id = first_graph->max_node_id();

    while (optind < argc) {

        unique_ptr<MutablePathMutableHandleGraph> graph;
        string graph_filename = get_input_file_name(optind, argc, argv);
        graph = vg::io::VPKG::load_one<MutablePathMutableHandleGraph>(graph_filename);

        // join the id spaces if necessary
        int64_t delta = max_node_id - graph->min_node_id();
        if (delta >= 0) {
            graph->increment_node_ids(delta + 1);
        }
        max_node_id = graph->max_node_id();

        if (connect_paths) {
            handlealgs::append_path_handle_graph(graph.get(), first_graph.get(), true);
        } else {
            graph->for_each_path_handle([&](path_handle_t path_handle) {
                    string path_name = graph->get_path_name(path_handle);
                    if (first_graph->has_path(path_name)) {
                        cerr << "error [vg combine]: Paths with name \"" << path_name << "\" found in multiple input graphs. If they are consecutive subpath ranges, they can be connected by using the -p option." << endl;
                        exit(1);
                    }
                });
            handlealgs::copy_path_handle_graph(graph.get(), first_graph.get());
        }        
    }

    // Serialize the graph using VPKG.
    vg::io::save_handle_graph(first_graph.get(), cout);

    return 0;
}

// Register subcommand
static Subcommand vg_combine("combine", "merge multiple graph files together", main_combine);


// This is the original vg combine logic, which itself mimics using "cat" to join up protobuf files
// Since it relies on the Protobuf format itself, particular the ability to stream together chunks that
// would otherwise be invalid individually, it is probably never going to be ported to the handle graph
// api, which is why it's been relegated to the deprecated bin
int cat_proto_graphs(int argc, char** argv) {
    
    while (optind < argc) {
        get_input_file(optind, argc, argv, [&](istream& in) {
            // We're producing output in uncompressed, "VG"-type-tagged, VPKG Protobuf format.
            // We will check if this file is uncompressed or compressed VG-type-tagged data.
            
            if (vg::io::BlockedGzipInputStream::SmellsLikeGzip(in)) {
                // It is compressed.
                
                // Save our start position
                auto start = in.tellg();
                
                {
                    // Try decompressing.
                    vg::io::BlockedGzipInputStream decompressed(in);
                    if (decompressed.IsBGZF() && vg::io::MessageIterator::sniff_tag(decompressed) == "VG") {
                        // We have Blocked GZIP which we can potentially just forward.
                        // It looks like compressed VG Protobuf data.
                        
                        // Decompress it all to stdout, using the ZeroCopyInputStream API.
                        char* buffer = nullptr;
                        int bytes = 0;
                        while (cout && decompressed.Next((const void**) &buffer, &bytes)) {
                            // Each time we get bytes, write them to stdout.
                            cout.write(buffer, bytes);
                        }
                        
                        if (!cout) {
                            cerr << "error [vg combine]: Could not write decompressed data to output stream." << endl;
                            exit(1);
                        }
                        
                        // Do the next input file
                        return;
                    }
                }
                
                // We may have hit EOF.
                in.clear();
                
                // If we get here, it wasn't compressed VG Protobuf.
                // So we need to go back to the start of the file, since the decompressor read some.
                in.seekg(start);
                
            } else if (vg::io::MessageIterator::sniff_tag(in) == "VG") {
                // It isn't compressed, but it looks like uncompressed VG Protobuf.
                // Send the uncompressed data to stdout.
                cout << in.rdbuf();
                
                if (!cout) {
                    cerr << "error [vg combine]: Could not write raw data to output stream." << endl;
                    exit(1);
                }
                
                // Do the next input file
                return;
            }
            
            // If we get here, it isn't compressed or uncompressed VG protobuf.
            // Read it as a PathHandleGraph
            unique_ptr<PathHandleGraph> graph = vg::io::VPKG::load_one<PathHandleGraph>(in);
            
            // Convert to vg::VG
            VG* vg_graph = dynamic_cast<vg::VG*>(graph.get());
            if (vg_graph == nullptr) {
                vg_graph = new vg::VG();
                handlealgs::copy_path_handle_graph(graph.get(), vg_graph);
                // Give the unique_ptr ownership and delete the graph we loaded.
                graph.reset(vg_graph);
                // Make sure the paths are all synced up
                vg_graph->paths.to_graph(vg_graph->graph);
            }
            
            {
                // Save to stdout, uncompressed
                vg::io::ProtobufEmitter<Graph> emitter(cout, false);
                vg_graph->serialize_to_emitter(emitter);
                // Make sure the emitter goes away and writes before we check on the stream.
            }
            
            if (!cout) {
                cerr << "error [vg combine]: Could not write converted graph to output stream." << endl;
                exit(1);
            }
            
        });
    }
    return 0;
}
