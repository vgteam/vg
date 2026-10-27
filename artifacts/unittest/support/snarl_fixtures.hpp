#ifndef VG_UNITTEST_SNARL_FIXTURES_HPP_INCLUDED
#define VG_UNITTEST_SNARL_FIXTURES_HPP_INCLUDED

/** \file
 * Conversions from the `Snarl` and `SnarlTraversal` that test fixtures are written in to the bounds
 * and walks the read-likelihood caller takes.
 */

#include <vector>

#include "decomposition_sites.hpp"
#include "handle.hpp"
#include "traversal_finder.hpp"
#include <vg/vg.pb.h>

namespace vg {
namespace unittest {

using namespace std;

/// Get the bounds of `snarl` in `graph`.
inline SiteBounds bounds_of(const HandleGraph& graph, const Snarl& snarl) {
    return SiteBounds{graph.get_handle(snarl.start().node_id(), snarl.start().backward()),
                      graph.get_handle(snarl.end().node_id(), snarl.end().backward())};
}

/// Convert a SnarlTraversal of node visits to a walk.
inline Traversal walk_of(const HandleGraph& graph, const SnarlTraversal& trav) {
    Traversal walk;
    walk.reserve(trav.visit_size());
    for (const Visit& visit : trav.visit()) {
        walk.push_back(graph.get_handle(visit.node_id(), visit.backward()));
    }
    return walk;
}

/// Convert each SnarlTraversal of node visits to a walk.
inline vector<Traversal> walks_of(const HandleGraph& graph, const vector<SnarlTraversal>& travs) {
    vector<Traversal> walks;
    walks.reserve(travs.size());
    for (const SnarlTraversal& trav : travs) {
        walks.push_back(walk_of(graph, trav));
    }
    return walks;
}

}
}

#endif
