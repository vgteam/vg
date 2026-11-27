/**
 * \file register_loader_gfaz.cpp
 */

#include "register_loader_gfaz.hpp"

#include <GFAz/core/codec/serialization.hpp>
#include <GFAz/handlegraph_adapter/gfaz_handle_graph.hpp>
#include <vg/io/registry.hpp>

#include "../algorithms/gfaz_to_handle.hpp"
#include "../gfa.hpp"
#include "../handle.hpp"

#include "save_handle_graph.hpp"

#include <cstdlib>
#include <iostream>

namespace vg {
namespace io {

using namespace std;
using namespace vg::io;

void register_loader_gfaz() {

    std::uint32_t magic_number = gfaz::GFAZ_MAGIC;
    std::string magic_string(reinterpret_cast<char*>(&magic_number), sizeof(magic_number));

    // Graph types that only need to be read are served straight from the
    // container: nothing is converted, and paths stay compressed and are
    // decoded a block at a time.
    //
    // VG_GFAZ_CONVERT=1 is a debugging switch, like VG_FULL_TRACEBACK: it loads
    // a GFAHandleGraph for every type instead, to tell whether a problem comes
    // from reading in place or from the file. It is an environment variable
    // because loaders are registered before any options are parsed.
    const char* convert_env = getenv("VG_GFAZ_CONVERT");
    bool always_convert = convert_env != nullptr && string(convert_env) != "0";

    if (!always_convert) {
        Registry::register_bare_loader_with_magic_and_filename<gfaz::GfazHandleGraph, PathPositionHandleGraph, PathHandleGraph, HandleGraph>("GFAZView", magic_string, [](istream& input, const string& filename) -> void* {
            try {
                gfaz::GfazHandleGraphOptions options;
                // Use vg's thread count instead of GFAz's own default.
                options.num_threads = -1;
                if (!filename.empty() && filename != "-") {
                    // Loading by name also picks up a saved path index.
                    return (void*) gfaz::GfazHandleGraph::load(filename, options).release();
                }
                return (void*) new gfaz::GfazHandleGraph(gfaz::deserialize_compressed_data(input), options);
            } catch (std::exception& e) {
                cerr << "error[register_loader_gfaz] GFAZ ";
                if (!filename.empty() && filename != "-") {
                    cerr << "file " << filename;
                } else {
                    cerr << "stream";
                }
                cerr << " cannot be loaded: " << e.what() << endl;
                exit(1);
            }
        });
    }

    auto convert_loader = [](istream& input, const string& filename) -> void* {
        GFAHandleGraph* gfa_graph = new GFAHandleGraph();
        try {
            // Load the graph
            algorithms::GFAzParser parser;
            algorithms::attach_parser(parser, &gfa_graph->gfa_id_space);
            algorithms::attach_parser(parser, gfa_graph);
            parser.parse(input);
            gfa_graph->gfa_id_space.invert_translation();
        } catch (algorithms::GFAFormatError& e) {
            cerr << "error[register_loader_gfaz] GFAZ ";
            if (!filename.empty() && filename != "-") {
                cerr << "file " << filename;
            } else {
                cerr << "stream";
            }
            cerr << " cannot be loaded: " << e.what() << endl;
            exit(1);
        }
        return (void*) gfa_graph;
    };

    if (always_convert) {
        Registry::register_bare_loader_with_magic_and_filename<GFAHandleGraph, MutablePathDeletableHandleGraph, MutablePathMutableHandleGraph, MutableHandleGraph, PathHandleGraph, HandleGraph>("GFAZ", magic_string, convert_loader);
    } else {
        // Types that get modified still need a real mutable graph.
        Registry::register_bare_loader_with_magic_and_filename<GFAHandleGraph, MutablePathDeletableHandleGraph, MutablePathMutableHandleGraph, MutableHandleGraph>("GFAZ", magic_string, convert_loader);
    }
}

}
}
