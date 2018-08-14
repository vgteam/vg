#ifndef VG_ALGORITHMS_IS_SINGLE_STRANDED_HPP_INCLUDED
#define VG_ALGORITHMS_IS_SINGLE_STRANDED_HPP_INCLUDED

/**
 * \file single_stranded.hpp
 *
 * Defines algorithms for deciding if a graph contains reversing edges.
 */

#include "../handle.hpp"
#include "apply_bulk_modifications.hpp"

#include <unordered_set>
#include <unordered_map>
#include <vector>

namespace vg {
namespace algorithms {

using namespace std;

    /// Returns true if the graph contains no reversing edges (i.e. edges that connected
    /// the locally forward orientation of a node to the locally reverse orientation of
    /// of another node).
    bool is_single_stranded(const HandleGraph* graph);

    /// Returns a vector of handles where the orientation of each handle indicates an
    /// orientation that could be used to convert the graph into a single-stranded graph.
    /// That is, if all of the reverse handles (or all of the forward handles) were swapped
    /// in orientation, the graph would contain no reversing edges. Returns an empty vector
    /// if there is no such combination of node orientations (also if graph has no nodes).
    vector<handle_t> single_stranded_orientation(const HandleGraph* graph);
        
    /// Finds a set of node orientations that can be applied so that there are no
    /// reversing edges (i.e. every edge connects a locally forward node traversal
    /// to another locally forward orientation). If no such combination of orientations
    /// exists, produces an error and exits. Returns a set of the node IDs for nodes that
    /// were swapped in orientation. Potentially invalidates any existing handles.
    unordered_set<id_t> make_single_stranded(MutableHandleGraph* graph);

}
}

#endif
