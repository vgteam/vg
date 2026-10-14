#ifndef VG_UNITTEST_DECOMPOSITION_PAIR_HPP_INCLUDED
#define VG_UNITTEST_DECOMPOSITION_PAIR_HPP_INCLUDED

/**
 * \file decomposition_pair.hpp
 * Provides DecompositionPair, which builds both snarl decomposition
 * implementations of one graph from one run of a snarl finder, and tools for
 * describing what the two implementations present in a way that can be
 * compared between them.
 */

#include <random>
#include <string>
#include <vector>

#include <bdsg/hash_graph.hpp>

#include "snarls.hpp"
#include "snarl_distance_index.hpp"
#include "snarl_manager_decomposition.hpp"
#include "snarl_decomposition_fuzzer.hpp"

namespace vg::unittest {

/**
 * A SnarlManager, a SnarlManagerDecomposition view of it, and a
 * SnarlDistanceIndex, all built from the same decomposition of a graph. The
 * finder's events are recorded once and replayed into the manager and the
 * distance index, so a finder that answers differently each time, like
 * SnarlDecompositionFuzzer, still gives both the same decomposition. The
 * graph must outlive the pair.
 */
class DecompositionPair {
public:
    /// Record the finder's decomposition of the graph and build both
    /// implementations from it.
    DecompositionPair(const HandleGraph& graph_to_decompose, const HandleGraphSnarlFinder& finder);

    /// Build both implementations from a recorded decomposition of the
    /// graph.
    DecompositionPair(const HandleGraph& graph, std::vector<DecompositionEvent> recorded);

    DecompositionPair(const DecompositionPair& other) = delete;
    DecompositionPair& operator=(const DecompositionPair& other) = delete;

    /// Return true if the decomposition holds anything that
    /// SnarlManagerDecomposition presents differently from a
    /// SnarlDistanceIndex: a snarl with no nodes that is not an ultrabubble,
    /// a snarl whose start and end are the same node, or a cyclic chain.
    bool has_known_difference() const;

    /// The graph decomposed.
    const HandleGraph& graph;
    /// The finder's events, as recorded.
    std::vector<DecompositionEvent> events;
    /// The finished manager built from the events.
    SnarlManager manager;
    /// The manager seen through the decomposition interface.
    SnarlManagerDecomposition view;
    /// The distance index built from the events.
    SnarlDistanceIndex distance_index;
};

/**
 * A small random graph, with each node stored reverse-complemented with
 * probability 1/2 so that snarls and chains meet nodes in both
 * orientations, and a DecompositionPair built from a
 * SnarlDecompositionFuzzer's decomposition of it.
 */
class RandomDecomposition {
public:
    /// Make the graph and decompose it, with the fuzzer reversing each chain
    /// with the given probability.
    RandomDecomposition(std::default_random_engine& generator, double chain_flip_probability);

    RandomDecomposition(const RandomDecomposition& other) = delete;
    RandomDecomposition& operator=(const RandomDecomposition& other) = delete;

    bdsg::HashGraph graph;
    DecompositionPair pair;
};

/// Describe a handle as its node ID followed by + if forward or - if
/// reverse.
std::string handle_string(const HandleGraph& graph, const handle_t& handle);

/// Describe a snarl or chain by its bounds, giving the same description
/// whichever way it is read.
std::string bounds_key(const handlegraph::SnarlDecomposition& decomposition, const HandleGraph& graph,
                       const net_handle_t& net);

/// Describe a list of items, giving the same description whether the list is
/// read forward or backward.
std::string unoriented_list(const std::vector<std::string>& items);

}

#endif
