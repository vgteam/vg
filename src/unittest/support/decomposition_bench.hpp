#ifndef VG_UNITTEST_DECOMPOSITION_BENCH_HPP_INCLUDED
#define VG_UNITTEST_DECOMPOSITION_BENCH_HPP_INCLUDED

/**
 * \file decomposition_bench.hpp
 * Provides DecompositionBench, which builds two SnarlDecomposition
 * implementations of one graph from one decomposition of it, and functions
 * that make random graphs and decompositions for it.
 */

#include <random>
#include <vector>

#include <bdsg/hash_graph.hpp>

#include "snarls.hpp"
#include "snarl_distance_index.hpp"
#include "snarl_manager_decomposition.hpp"
#include "snarl_decomposition_fuzzer.hpp"

namespace vg::unittest {

/**
 * @brief A test bench that runs a snarl finder on a graph, and builds a
 * SnarlManagerDecomposition and a SnarlDistanceIndex from that one run.
 *
 * @details The bench records the snarl finder's events and replays them into
 * each implementation, so a finder that decomposes the graph differently on
 * each run, such as SnarlDecompositionFuzzer, still gives both the same
 * decomposition.
 */
class DecompositionBench {
public:
    /**
     * @brief Decompose a graph with a snarl finder.
     *
     * @param graph The graph to decompose. The bench keeps a reference to
     * it, so it must outlive the bench.
     * @param finder The snarl finder, which the constructor runs.
     */
    DecompositionBench(const HandleGraph& graph, const HandleGraphSnarlFinder& finder);

    /**
     * @brief Build the implementations from a recorded decomposition.
     *
     * @param graph The graph that was decomposed. The bench keeps a
     * reference to it, so it must outlive the bench.
     * @param events The snarl finder's events, from capture_events().
     */
    DecompositionBench(const HandleGraph& graph, const std::vector<DecompositionEvent>& events);

    DecompositionBench(const DecompositionBench& other) = delete;
    DecompositionBench& operator=(const DecompositionBench& other) = delete;

    /// @brief Get the graph that was decomposed.
    const HandleGraph& graph() const;
    /// @brief Get the SnarlManagerDecomposition built from the decomposition.
    const SnarlManagerDecomposition& manager_view() const;
    /// @brief Get the SnarlDistanceIndex built from the same decomposition.
    const SnarlDistanceIndex& distance_index() const;

private:
    /// @brief The graph that was decomposed.
    const HandleGraph& decomposed_graph;
    /// @brief The snarls that view presents, which must outlive it.
    SnarlManager snarl_manager;
    /// @brief The view of snarl_manager.
    SnarlManagerDecomposition view;
    /// @brief The distance index built from the events.
    SnarlDistanceIndex index;
};

/**
 * @brief Make a small random graph, with each node stored
 * reverse-complemented with probability 1/2.
 *
 * @details Storing nodes in both orientations makes snarls and chains meet
 * them in both orientations.
 */
bdsg::HashGraph random_flipped_graph(std::default_random_engine& generator);

/**
 * @brief Record a decomposition of a graph in which chains are reversed at
 * random.
 *
 * @details Runs an IntegratedSnarlFinder through a SnarlDecompositionFuzzer,
 * which reverses each chain with the given probability.
 *
 * @return The finder's events, for a DecompositionBench.
 */
std::vector<DecompositionEvent> fuzzed_decomposition(const HandleGraph& graph, double chain_flip_probability,
                                                     std::default_random_engine& generator);

}

#endif
