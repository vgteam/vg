#ifndef VG_UNITTEST_DECOMPOSITION_PAIR_HPP_INCLUDED
#define VG_UNITTEST_DECOMPOSITION_PAIR_HPP_INCLUDED

/**
 * \file decomposition_pair.hpp
 * Provides DecompositionPair, which builds both snarl decomposition
 * implementations of one graph from one run of a snarl finder.
 */

#include <vector>

#include "snarls.hpp"
#include "snarl_distance_index.hpp"
#include "snarl_manager_decomposition.hpp"
#include "snarl_decomposition_fuzzer.hpp"

namespace vg {
namespace unittest {

/**
 * A SnarlManager, its SnarlManagerDecomposition, and a SnarlDistanceIndex,
 * all built from the same decomposition of a graph. The finder's events are
 * recorded once and replayed into both, so a finder that answers differently
 * each time, like SnarlDecompositionFuzzer, still gives both the same
 * decomposition. The graph must outlive the pair.
 */
class DecompositionPair {
public:
    /// Record the finder's decomposition of the graph and build both
    /// implementations from it.
    DecompositionPair(const HandleGraph& graph, const HandleGraphSnarlFinder& finder);

    DecompositionPair(const DecompositionPair& other) = delete;
    DecompositionPair& operator=(const DecompositionPair& other) = delete;

    /// Return true if the decomposition holds anything the two
    /// implementations are known to present differently: a snarl with no
    /// nodes that is not an ultrabubble, a unary snarl, or a cyclic chain.
    bool has_known_difference() const;

    /// The graph decomposed.
    const HandleGraph& graph;
    /// The finder's events, as recorded.
    std::vector<DecompositionEvent> events;
    /// The finished manager built from the events.
    SnarlManager manager;
    /// The manager seen through the decomposition interface.
    SnarlManagerDecomposition adapter;
    /// The distance index built from the events.
    SnarlDistanceIndex distance_index;
};

}
}

#endif
