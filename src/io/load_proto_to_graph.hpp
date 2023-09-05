#ifndef VG_IO_LOAD_PROTO_TO_GRAPH_HPP_INCLUDED
#define VG_IO_LOAD_PROTO_TO_GRAPH_HPP_INCLUDED

/**
 * \file load_proto_to_graph.hpp
 * Read VG Protobuf into any MutablePathMutableHandleGraph.
 * Also useful for converting streams of Protobuf Graph objects into a MutablePathMutableHandleGraph.
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
void load_proto_to_graph(vg::MutablePathMutableHandleGraph* destination, const vg::io::message_sender_function_t& for_each_message);

/**
 * Call the given function with a callback which it can call with a series of
 * Protobuf Graph objects, possibly in multiple threads. The Protobuf Graph
 * objects may have dangling edges.
 *
 * Resolves all the dangling edges and writes all the graph data into the given
 * MutablePathMutableHandleGraph, with the destination graph being protected
 * from concurrent modification.
 */
void load_proto_to_graph(vg::MutablePathMutableHandleGraph* destination, const function<void(const function<void(Graph&)>&)>& chunk_sender);

}

}

#endif
