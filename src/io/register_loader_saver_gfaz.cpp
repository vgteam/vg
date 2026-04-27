/**
 * \file register_loader_saver_gfaz.cpp
 */

#include <vg/io/registry.hpp>
#include "register_loader_saver_gfaz.hpp"
#include "algorithms/gfa_to_handle.hpp"

#include "gfa.hpp"

#include "handle.hpp"
#include "save_handle_graph.hpp"
#include "../utility.hpp"

#include <cstdio>
#include <fstream>
#include <iostream>
#include <unistd.h>

namespace vg {
namespace io {

using namespace std;
using namespace vg::io;

static bool sniff_gfaz(istream& stream) {
    if (!stream) {
        return false;
    }
    char buffer[sizeof(uint32_t)];
    size_t buffer_used = 0;
    while (stream.peek() != EOF && buffer_used < sizeof(buffer)) {
        buffer[buffer_used++] = (char) stream.get();
    }
    for (size_t i = 0; i < buffer_used; ++i) {
        stream.unget();
        if (!stream) {
            throw runtime_error("Ungetting failed after " + to_string(i) + " characters");
        }
    }
    if (!stream) {
        stream.clear();
    }
    return algorithms::GFAzParser::has_magic(buffer, buffer_used);
}

void register_loader_saver_gfaz() {
    Registry::register_bare_loader_saver_with_header_check<GFAHandleGraph, MutablePathDeletableHandleGraph, MutablePathMutableHandleGraph, MutableHandleGraph, PathHandleGraph, HandleGraph>("GFAZ", sniff_gfaz, [](istream& input, const string& filename) -> void* {
        GFAHandleGraph* gfa_graph = new GFAHandleGraph();
        string temp_name;
        try {
            string load_name = filename;
            if (load_name.empty() || load_name == "-") {
                temp_name = temp_file::create("gfaz-load");
                ofstream temp_stream(temp_name, ios::binary);
                temp_stream << input.rdbuf();
                temp_stream.close();
                load_name = temp_name;
            }
            algorithms::GFAzParser parser;
            parser.external_id_map = &gfa_graph->gfa_id_space;
            algorithms::parser_to_path_handle_graph(parser, gfa_graph);
            parser.parse(load_name);
            gfa_graph->gfa_id_space.invert_translation();
        } catch (algorithms::GFAFormatError& e) {
            cerr << "error[register_loader_saver_gfaz] GFAZ ";
            if (!filename.empty() && filename != "-") {
                cerr << "file " << filename;
            } else {
                cerr << "stream";
            }
            cerr << " cannot be loaded: " << e.what() << endl;
            exit(1);
        }
        if (!temp_name.empty()) {
            unlink(temp_name.c_str());
        }
        return (void*) gfa_graph;
    }, [](const void* gfa_graph_void, ostream& output) {
        assert(gfa_graph_void != nullptr);
        graph_to_gfa((const GFAHandleGraph*) gfa_graph_void, output);
    });
}

}
}
