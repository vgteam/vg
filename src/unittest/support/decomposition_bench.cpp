/**
 * \file decomposition_bench.cpp
 * Implements DecompositionBench and the random decompositions for it.
 */

#include "decomposition_bench.hpp"

#include "../../integrated_snarl_finder.hpp"
#include "../../vg.hpp"
#include "random_graph.hpp"
#include "randomly_flipped_nodes.hpp"

namespace vg::unittest {

DecompositionBench::DecompositionBench(const HandleGraph& graph, const HandleGraphSnarlFinder& finder) :
    DecompositionBench(graph, capture_events(finder, graph)) {
}

DecompositionBench::DecompositionBench(const HandleGraph& graph, const std::vector<DecompositionEvent>& events) :
    decomposed_graph(graph),
    snarl_manager(ReplaySnarlFinder(&graph, events).find_snarls()),
    view(snarl_manager, graph) {

    ReplaySnarlFinder replay(&graph, events);
    fill_in_distance_index(&index, &graph, &replay);
}

const HandleGraph& DecompositionBench::graph() const {
    return decomposed_graph;
}

const SnarlManagerDecomposition& DecompositionBench::manager_view() const {
    return view;
}

const SnarlDistanceIndex& DecompositionBench::distance_index() const {
    return index;
}

bdsg::HashGraph random_flipped_graph(std::default_random_engine& generator) {
    VG base_graph;
    small_random_graph(generator, &base_graph);
    return randomly_flipped_nodes(base_graph, 0.5, generator);
}

std::vector<DecompositionEvent> fuzzed_decomposition(const HandleGraph& graph, double chain_flip_probability,
                                                     std::default_random_engine& generator) {
    IntegratedSnarlFinder base_finder(graph);
    SnarlDecompositionFuzzer finder(&graph, &base_finder, chain_flip_probability, generator);
    return capture_events(finder, graph);
}

}
