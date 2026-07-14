/**
 * \file register_loader_gfaz.cpp
 */

#include "register_loader_gfaz.hpp"

#include <GFAz/core/codec/serialization.hpp>
#include <vg/io/registry.hpp>

#include "../algorithms/gfaz_to_handle.hpp"
#include "../gfa.hpp"
#include "../handle.hpp"
#include "../utility.hpp"

#include "save_handle_graph.hpp"

#include <cstdio>
#include <fstream>
#include <iostream>
#include <unistd.h>

namespace vg {
namespace io {

using namespace std;
using namespace vg::io;

void register_loader_gfaz() {

    std::uint32_t magic_number = gfaz::GFAZ_MAGIC;
    std::string magic_string(reinterpret_cast<char*>(&magic_number), sizeof(magic_number));

    Registry::register_bare_loader_with_magic_and_filename<GFAHandleGraph, MutablePathDeletableHandleGraph, MutablePathMutableHandleGraph, MutableHandleGraph, PathHandleGraph, HandleGraph>("GFAZ", magic_string, [](istream& input, const string& filename) -> void* {
        GFAHandleGraph* gfa_graph = new GFAHandleGraph();
        string temp_name;
        try {
            string load_name = filename;
            if (load_name.empty() || load_name == "-") {
                // Dump the whole stream to a temporary file
                // TODO: When GFAz lets us load from a stream we can stop doing this.
                // TODO: When we get a stream, is the filename meant to be empty or "-"? One of those is wrong!
                temp_name = temp_file::create("gfaz-load");
                ofstream temp_stream(temp_name, ios::binary);
                temp_stream << input.rdbuf();
                temp_stream.close();
                load_name = temp_name;
            }
            // Load the graph
            algorithms::GFAzParser parser;
            algorithms::attach_parser(parser, &gfa_graph->gfa_id_space);
            algorithms::attach_parser(parser, gfa_graph);
            parser.parse(load_name);
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
        if (!temp_name.empty()) {
            // Clean up the temporary file
            unlink(temp_name.c_str());
        }
        return (void*) gfa_graph;
    });
}

}
}
