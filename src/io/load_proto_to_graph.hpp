#ifndef VG_IO_LOAD_PROTO_TO_GRAPH_HPP_INCLUDED
#define VG_IO_LOAD_PROTO_TO_GRAPH_HPP_INCLUDED

/**
 * \file load_proto_to_graph.hpp
 * Read VPKG-packaged VG Protobuf into any MutablePathMutableHandleGraph.
 */

#include "../handle.hpp"
#include <vg/io/vpkg.hpp>
#include <vg/io/registry.hpp>

namespace vg {

namespace io {

using namespace std;
using namespace vg;

/**
 * Read all string messages supplied by the given message sender as Protobuf
 * Graph objects, and create the specified graph in the destination graph.
 *
 * Paths need to be cached until the end for ranks to be respected.
 */
void load_proto_to_graph(const vg::io::message_sender_function_t& for_each_message, vg::MutablePathMutableHandleGraph* destination);

}

}

#endif
