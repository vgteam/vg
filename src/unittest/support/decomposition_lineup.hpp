#ifndef VG_UNITTEST_DECOMPOSITION_LINEUP_HPP_INCLUDED
#define VG_UNITTEST_DECOMPOSITION_LINEUP_HPP_INCLUDED

/**
 * \file decomposition_lineup.hpp
 * Test inputs for SnarlDecomposition implementations: random graphs, random
 * decompositions of them, and a line-up of every implementation built from
 * one decomposition.
 */

#include <deque>
#include <functional>
#include <memory>
#include <random>
#include <string>
#include <vector>

#include <bdsg/hash_graph.hpp>
#include <handlegraph/snarl_decomposition.hpp>

#include "snarls.hpp"
#include "snarl_decomposition_fuzzer.hpp"

namespace vg::unittest {

using handlegraph::SnarlDecomposition;

/**
 * A line-up of every SnarlDecomposition implementation, each built
 * from the same recorded decomposition of a graph.
 *
 * Tests that run one check on each implementation loop over the
 * line-up with for_each(), so a new implementation needs only a new entry in
 * the constructor.
 */
class DecompositionLineup {
public:
    /**
     * Build each implementation from a recorded decomposition.
     *
     * @param graph The graph that was decomposed. The implementations keep
     * references to it, so it must outlive the line-up.
     * @param events The events of a snarl finder's decomposition of the
     * graph, as capture_events() records them.
     */
    DecompositionLineup(const HandleGraph& graph, const std::vector<DecompositionEvent>& events);

    /**
     * Build each implementation from a snarl finder's decomposition.
     *
     * Run the finder once and record its events, so that every
     * implementation gets the same decomposition from a finder that
     * decomposes differently on each run.
     *
     * @param graph As for the other constructor.
     * @param finder The snarl finder to run.
     */
    DecompositionLineup(const HandleGraph& graph, const HandleGraphSnarlFinder& finder);

    DecompositionLineup(const DecompositionLineup& other) = delete;
    DecompositionLineup& operator=(const DecompositionLineup& other) = delete;

    /**
     * Call a function on each implementation, with the name of its
     * class.
     *
     * @post Implementations come in the same order on each call, and in each
     * line-up.
     */
    void for_each(const std::function<void(const std::string&, const SnarlDecomposition&)>& iteratee) const;

private:
    /// One implementation, with the name of its class.
    struct Entry {
        std::string name;
        std::unique_ptr<SnarlDecomposition> decomposition;
    };

    /// The SnarlManagers that the SnarlManagerDecomposition entries view.
    /// Each view keeps a reference to its manager, and adding to the deque
    /// never moves a manager.
    std::deque<SnarlManager> managers;
    /// The implementations, in the order the constructor builds them.
    std::vector<Entry> entries;
};

/**
 * Make a small random graph, with each node stored
 * reverse-complemented with probability 1/2.
 *
 * A node stored reverse-complemented has a reverse handle where the
 * original graph had a forward one. So chains in the result have some nodes
 * reversed in parent.
 */
bdsg::HashGraph random_flipped_graph(std::default_random_engine& generator);

/**
 * Record a decomposition of a graph, with each chain reversed at
 * random.
 *
 * Run an IntegratedSnarlFinder through a SnarlDecompositionFuzzer, which
 * reverses each chain with the given probability.
 *
 * @return The finder's events.
 */
std::vector<DecompositionEvent> fuzzed_decomposition(const HandleGraph& graph, double chain_flip_probability,
                                                     std::default_random_engine& generator);

}

#endif
