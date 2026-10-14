/**
 * \file decomposition_lineup.cpp
 * Builds the SnarlDecomposition implementations that tests compare, and the
 * random inputs they are tested on.
 */

#include "decomposition_lineup.hpp"

#include "../../integrated_snarl_finder.hpp"
#include "../../snarl_distance_index.hpp"
#include "../../snarl_manager_decomposition.hpp"
#include "../../vg.hpp"
#include "random_graph.hpp"
#include "randomly_flipped_nodes.hpp"

namespace vg::unittest {

DecompositionLineup::DecompositionLineup(const HandleGraph& graph, const HandleGraphSnarlFinder& finder) :
    DecompositionLineup(graph, capture_events(finder, graph)) {
}

DecompositionLineup::DecompositionLineup(const HandleGraph& graph, const std::vector<DecompositionEvent>& events) {
    managers.push_back(ReplaySnarlFinder(&graph, events).find_snarls());
    entries.push_back({"SnarlManagerDecomposition",
                       std::make_unique<SnarlManagerDecomposition>(managers.back(), graph)});

    auto index = std::make_unique<SnarlDistanceIndex>();
    ReplaySnarlFinder replay(&graph, events);
    fill_in_distance_index(index.get(), &graph, &replay);
    entries.push_back({"SnarlDistanceIndex", std::move(index)});
}

void DecompositionLineup::for_each(
    const std::function<void(const std::string&, const SnarlDecomposition&)>& iteratee) const {
    for (const Entry& entry : entries) {
        iteratee(entry.name, *entry.decomposition);
    }
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
