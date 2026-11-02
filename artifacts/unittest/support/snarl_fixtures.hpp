#ifndef VG_UNITTEST_SNARL_FIXTURES_HPP_INCLUDED
#define VG_UNITTEST_SNARL_FIXTURES_HPP_INCLUDED

/** \file
 * Conversions from the `Snarl` and `SnarlTraversal` that test fixtures are written in to the bounds
 * and walks the read-likelihood caller takes.
 */

#include <vector>

#include "../../decomposition_sites.hpp"
#include "../../../src/handle.hpp"
#include "../../materialized_site.hpp"
#include "../../site_values.hpp"
#include "../../../src/traversal_finder.hpp"
#include <vg/vg.pb.h>

namespace vg {
namespace unittest {

using namespace std;

/// Get the nodes holding the start and end bounds of `snarl` in `graph`.
inline pair<handle_t, handle_t> bounds_of(const HandleGraph& graph, const Snarl& snarl) {
    return {graph.get_handle(snarl.start().node_id(), snarl.start().backward()),
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

/// Write out `snarl` with the alleles `alleles`.
inline MaterializedSite materialized_of(const HandleGraph& graph, const Snarl& snarl,
                                        const CandidateAlleles& alleles) {
    const handle_t start = graph.get_handle(snarl.start().node_id(), snarl.start().backward());
    const handle_t end = graph.get_handle(snarl.end().node_id(), snarl.end().backward());
    return MaterializedSite{start, graph.get_length(start), end, 0, alleles};
}

/// Write out the edit site `edit`, which runs from its first base to the end of its walk.
inline MaterializedSite materialized_of(const HandleGraph& graph, const EditSite& edit) {
    const Traversal& walk = edit.alleles.walks.front();
    return MaterializedSite{walk.front(), edit.start_offset, walk.back(),
                            graph.get_length(walk.back()), edit.alleles};
}

}
}

#endif
