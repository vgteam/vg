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
#include "../algorithms/copy_graph.hpp"
#include "../vg.hpp"

using namespace std;
using namespace vg;
using namespace vg::subcommand;

void help_combine(char** argv) {
    cerr << "usage: " << argv[0] << " combine [options] <graph1.vg> [graph2.vg ...] >merged.vg" << endl
        << "Combines one or more graphs into a single file, regardless of input format." << endl
        << "Handling of duplicate nodes or edges is undefined. Graphs must be normal, seekable files." << endl;
}

int main_combine(int argc, char** argv) {

    if (argc == 2) {
        help_combine(argv);
        return 1;
    }

    int c;
    optind = 2; // force optind past command positional argument
    while (true) {
        static struct option long_options[] =
        {
            {"help", no_argument, 0, 'h'},
            {0, 0, 0, 0}
        };

        int option_index = 0;
        c = getopt_long (argc, argv, "h",
                long_options, &option_index);

        // Detect the end of the options.
        if (c == -1)
            break;

        switch (c)
        {
            case 'h':
            case '?':
                help_combine(argv);
                exit(1);
                break;

            default:
                abort ();
        }
    }

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
                algorithms::copy_path_handle_graph(graph.get(), vg_graph);
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

// Register subcommand
static Subcommand vg_combine("combine", "merge multiple graph files together", main_combine);

