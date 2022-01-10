#ifndef VG_HANDLE_HPP_INCLUDED
#define VG_HANDLE_HPP_INCLUDED

/** \file 
 * One stop shop for libhandlegraph types and things we need to work with them.
 */

#include <handlegraph/util.hpp>

#include <handlegraph/handle_graph.hpp>
#include <handlegraph/mutable_handle_graph.hpp>
#include <handlegraph/path_handle_graph.hpp>
#include <handlegraph/path_position_handle_graph.hpp>
#include <handlegraph/mutable_path_handle_graph.hpp>
#include <handlegraph/mutable_path_mutable_handle_graph.hpp>
#include <handlegraph/deletable_handle_graph.hpp>
#include <handlegraph/mutable_path_deletable_handle_graph.hpp>
#include <handlegraph/named_node_back_translation.hpp>

#include <handlegraph/expanding_overlay_graph.hpp>
#include <handlegraph/serializable_handle_graph.hpp>

#include <handlegraph/algorithms/append_graph.hpp>
#include <handlegraph/algorithms/apply_orientations.hpp>
#include <handlegraph/algorithms/are_equivalent.hpp>
#include <handlegraph/algorithms/copy_graph.hpp>
#include <handlegraph/algorithms/count_walks.hpp>
#include <handlegraph/algorithms/dagify.hpp>
#include <handlegraph/algorithms/dijkstra.hpp>
#include <handlegraph/algorithms/eades_algorithm.hpp>
#include <handlegraph/algorithms/extend.hpp>
#include <handlegraph/algorithms/find_shortest_paths.hpp>
#include <handlegraph/algorithms/find_tips.hpp>
#include <handlegraph/algorithms/is_acyclic.hpp>
#include <handlegraph/algorithms/reverse_complement.hpp>
#include <handlegraph/algorithms/split_strands.hpp>
#include <handlegraph/algorithms/strongly_connected_components.hpp>
#include <handlegraph/algorithms/topological_sort.hpp>
#include <handlegraph/algorithms/chop.hpp>
#include <handlegraph/algorithms/weakly_connected_components.hpp>


#include "hash_map.hpp"
#include <vg/vg.pb.h>
#include "types.hpp"

namespace vg {

using namespace std;

namespace handlealgs = handlegraph::algorithms;

// Import all the handle stuff into the vg namespace for transition purposes.
using handle_t = handlegraph::handle_t;
using nid_t = handlegraph::nid_t;
using path_handle_t = handlegraph::path_handle_t;
using step_handle_t = handlegraph::step_handle_t;
using edge_t = handlegraph::edge_t;
using oriented_node_range_t = handlegraph::oriented_node_range_t;

using HandleGraph = handlegraph::HandleGraph;
using RankedHandleGraph = handlegraph::RankedHandleGraph;
using MutableHandleGraph = handlegraph::MutableHandleGraph;
using PathHandleGraph = handlegraph::PathHandleGraph;
using PathPositionHandleGraph = handlegraph::PathPositionHandleGraph;
using MutablePathHandleGraph = handlegraph::MutablePathHandleGraph;
using MutablePathMutableHandleGraph = handlegraph::MutablePathMutableHandleGraph;
using DeletableHandleGraph = handlegraph::DeletableHandleGraph;
using MutablePathDeletableHandleGraph = handlegraph::MutablePathDeletableHandleGraph;
using SerializableHandleGraph = handlegraph::SerializableHandleGraph;
using VectorizableHandleGraph = handlegraph::VectorizableHandleGraph;
using NamedNodeBackTranslation = handlegraph::NamedNodeBackTranslation;

/**
 * Define wang hashes for handles.
 */
template<>
struct wang_hash<handle_t> {
    size_t operator()(const handlegraph::handle_t& handle) const {
        return wang_hash<std::int64_t>()(handlegraph::as_integer(handle));
    }
};

template<>
struct wang_hash<path_handle_t> {
    size_t operator()(const handlegraph::path_handle_t& handle) const {
        return wang_hash<std::int64_t>()(handlegraph::as_integer(handle));
    }
};
    
using handlegraph::ExpandingOverlayGraph;
}

#endif
