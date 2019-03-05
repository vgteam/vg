#ifndef VG_HANDLE_HPP_INCLUDED
#define VG_HANDLE_HPP_INCLUDED

/** \file 
 * One stop shop for libhandlegraph types and things we need to work with them.
 */

#include <handlegraph/handle_graph.hpp>
#include <handlegraph/mutable_handle_graph.hpp>
#include <handlegraph/path_handle_graph.hpp>
#include <handlegraph/mutable_path_handle_graph.hpp>
#include <handlegraph/mutable_path_mutable_handle_graph.hpp>
#include <handlegraph/deletable_handle_graph.hpp>
#include <handlegraph/mutable_path_deletable_handle_graph.hpp>
#include <handlegraph/util.hpp>

#include "hash_map.hpp"
#include "vg.pb.h"
#include "types.hpp"

namespace vg {

using namespace std;

// Import all the handle stuff into the vg namespace for transition purposes.
using handle_t = handlegraph::handle_t;
using path_handle_t = handlegraph::path_handle_t;
using occurrence_handle_t = handlegraph::occurrence_handle_t;
using edge_t = handlegraph::edge_t;

using HandleGraph = handlegraph::HandleGraph;
using MutableHandleGraph = handlegraph::MutableHandleGraph;
using PathHandleGraph = handlegraph::PathHandleGraph;
using MutablePathHandleGraph = handlegraph::MutablePathHandleGraph;
using MutablePathMutableHandleGraph = handlegraph::MutablePathMutableHandleGraph;
using DeletableHandleGraph = handlegraph::DeletableHandleGraph;
using MutablePathDeletableHandleGraph = handlegraph::MutablePathDeletableHandleGraph;

/**
 * Define wang hashes for handles.
 */
template<>
struct wang_hash<handle_t> {
    size_t operator()(const handlegraph::handle_t& handle) const {
        return wang_hash<std::int64_t>()(handlegraph::as_integer(handle));
    }
};

}

#endif
