/**
 * \file register_loader_gfaz.cpp
 */

#include "register_loader_gfaz.hpp"

#include <GFAz/core/codec/serialization.hpp>
#include <vg/io/registry.hpp>

#include "../algorithms/gfaz_to_handle.hpp"
#include "../gfa.hpp"
#include "../handle.hpp"

#include "save_handle_graph.hpp"

#include <iostream>

namespace vg {
namespace io {

using namespace std;
using namespace vg::io;

void register_loader_gfaz() {

    std::uint32_t magic_number = gfaz::GFAZ_MAGIC;
    std::string magic_string(reinterpret_cast<char*>(&magic_number), sizeof(magic_number));

    Registry::register_bare_loader_with_magic_and_filename<GFAHandleGraph, MutablePathDeletableHandleGraph, MutablePathMutableHandleGraph, MutableHandleGraph, PathHandleGraph, HandleGraph>("GFAZ", magic_string, [](istream& input, const string& filename) -> void* {
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
    });
}

}
}
