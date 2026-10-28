#ifndef VG_UNITTEST_SNARL_FIXTURES_HPP_INCLUDED
#define VG_UNITTEST_SNARL_FIXTURES_HPP_INCLUDED

/** \file
 * Conversions from the `Snarl` and `SnarlTraversal` that test fixtures are written in to the bounds
 * and walks the read-likelihood caller takes.
 */

#include <vector>

#include "decomposition_sites.hpp"
#include "handle.hpp"
#include "materialized_site.hpp"
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

/// Write out `snarl` with the alleles `alleles`.
inline MaterializedSite materialized_of(const HandleGraph& graph, const Snarl& snarl,
                                        const CandidateAlleles& alleles) {
    const handle_t start = graph.get_handle(snarl.start().node_id(), snarl.start().backward());
    const handle_t end = graph.get_handle(snarl.end().node_id(), snarl.end().backward());
    return MaterializedSite{start, graph.get_length(start), end, 0, alleles};
}

/// Get the alleles of the edit site `edit`.
inline CandidateAlleles alleles_of(const EditSite& edit) {
    return CandidateAlleles{edit.travs, 0, edit.sequences};
}

/// Write out the edit site `edit`, with its alleles `alleles`, which runs from its first base to
/// the end of its walk.
inline MaterializedSite materialized_of(const HandleGraph& graph, const EditSite& edit,
                                        const CandidateAlleles& alleles) {
    return MaterializedSite{edit.bounds.start, bases_before_site(graph, edit.bounds),
                            edit.bounds.end, graph.get_length(edit.bounds.end), alleles};
}

}
}

#endif
