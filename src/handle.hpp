#ifndef VG_HANDLE_HPP_INCLUDED
#define VG_HANDLE_HPP_INCLUDED

/** \file 
 * One stop shop for libhandlegraph types and things we need to work with them.
 */

#include <handlegraph/handle_graph.hpp>
#include <handlegraph/mutable_handle_graph.hpp>
#include <handlegraph/path_handle_graph.hpp>
#include <handlegraph/path_position_handle_graph.hpp>
#include <handlegraph/mutable_path_handle_graph.hpp>
#include <handlegraph/mutable_path_mutable_handle_graph.hpp>
#include <handlegraph/deletable_handle_graph.hpp>
#include <handlegraph/mutable_path_deletable_handle_graph.hpp>
#include <handlegraph/util.hpp>

#include "hash_map.hpp"
#include <vg/vg.pb.h>
#include "types.hpp"

namespace vg {

using namespace std;

// Import all the handle stuff into the vg namespace for transition purposes.
using handle_t = handlegraph::handle_t;
using path_handle_t = handlegraph::path_handle_t;
using step_handle_t = handlegraph::step_handle_t;
using edge_t = handlegraph::edge_t;

using HandleGraph = handlegraph::HandleGraph;
using MutableHandleGraph = handlegraph::MutableHandleGraph;
using PathHandleGraph = handlegraph::PathHandleGraph;
using PathPositionHandleGraph = handlegraph::PathPositionHandleGraph;
using MutablePathHandleGraph = handlegraph::MutablePathHandleGraph;
using MutablePathMutableHandleGraph = handlegraph::MutablePathMutableHandleGraph;
using DeletableHandleGraph = handlegraph::DeletableHandleGraph;
using MutablePathDeletableHandleGraph = handlegraph::MutablePathDeletableHandleGraph;
using SerializableHandleGraph = handlegraph::SerializableHandleGraph;

/**
 * Define wang hashes for handles.
 */
template<>
struct wang_hash<handle_t> {
    size_t operator()(const handlegraph::handle_t& handle) const {
        return wang_hash<std::int64_t>()(handlegraph::as_integer(handle));
    }
};
    

/**
 * This is the interface for a graph that represents a transformation of some underlying
 * HandleGraph where every node in the overlay corresponds to a node in the underlying
 * graph, but where more than one node in the overlay can map to the same underlying node.
 */
class ExpandingOverlayGraph : public HandleGraph {

public:
    
    /**
     * Returns the handle in the underlying graph that corresponds to a handle in the
     * overlay
     */
    virtual handle_t get_underlying_handle(const handle_t& handle) const = 0;
};

}

#endif
