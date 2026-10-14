/**
 * \file decomposition_pair.cpp
 * Implements DecompositionPair and the tools for describing decompositions.
 */

#include "decomposition_pair.hpp"

#include <algorithm>
#include <format>

#include "../../integrated_snarl_finder.hpp"
#include "../../vg.hpp"
#include "random_graph.hpp"
#include "randomly_flipped_nodes.hpp"

namespace vg::unittest {

DecompositionPair::DecompositionPair(const HandleGraph& graph_to_decompose, const HandleGraphSnarlFinder& finder) :
    DecompositionPair(graph_to_decompose, capture_events(finder, graph_to_decompose)) {
}

DecompositionPair::DecompositionPair(const HandleGraph& graph, std::vector<DecompositionEvent> recorded) :
    graph(graph),
    events(std::move(recorded)),
    manager(ReplaySnarlFinder(&graph, events).find_snarls()),
    view(manager, graph) {

    ReplaySnarlFinder replay(&graph, events);
    fill_in_distance_index(&distance_index, &graph, &replay);
}

bool DecompositionPair::has_known_difference() const {
    bool found = false;
    manager.for_each_snarl_preorder([this, &found](const Snarl* snarl) {
        bool one_node = snarl->start().node_id() == snarl->end().node_id();
        bool empty_non_ultrabubble = snarl->type() != ULTRABUBBLE && manager.is_leaf(snarl) &&
            manager.shallow_contents(snarl, graph, false).first.empty();
        found = found || one_node || empty_non_ultrabubble;
    });
    manager.for_each_chain([&found](const Chain* chain) {
        found = found || get_start_of(*chain) == get_end_of(*chain);
    });
    return found;
}

/// Make a small random graph, storing each node reverse-complemented with
/// probability 1/2.
static bdsg::HashGraph random_flipped_graph(std::default_random_engine& generator) {
    VG base_graph;
    small_random_graph(generator, &base_graph);
    return randomly_flipped_nodes(base_graph, 0.5, generator);
}

/// Record a decomposition of the graph with each chain flipped with the
/// given probability.
static std::vector<DecompositionEvent> fuzzed_events(const HandleGraph& graph, double chain_flip_probability,
                                                     std::default_random_engine& generator) {
    IntegratedSnarlFinder base_finder(graph);
    SnarlDecompositionFuzzer finder(&graph, &base_finder, chain_flip_probability, generator);
    return capture_events(finder, graph);
}

RandomDecomposition::RandomDecomposition(std::default_random_engine& generator, double chain_flip_probability) :
    graph(random_flipped_graph(generator)),
    pair(graph, fuzzed_events(graph, chain_flip_probability, generator)) {
}

std::string handle_string(const HandleGraph& graph, const handle_t& handle) {
    return std::format("{}{}", graph.get_id(handle), graph.get_is_reverse(handle) ? "-" : "+");
}

std::string bounds_key(const handlegraph::SnarlDecomposition& decomposition, const HandleGraph& graph,
                       const net_handle_t& net) {
    handle_t start = decomposition.get_handle(decomposition.get_bound(net, false, true), &graph);
    handle_t end = decomposition.get_handle(decomposition.get_bound(net, true, false), &graph);
    return std::min(handle_string(graph, start) + "/" + handle_string(graph, end),
                    handle_string(graph, graph.flip(end)) + "/" + handle_string(graph, graph.flip(start)));
}

std::string unoriented_list(const std::vector<std::string>& items) {
    std::string forward;
    std::string backward;
    for (size_t i = 0; i < items.size(); i++) {
        forward += items[i] + ",";
        backward += items[items.size() - 1 - i] + ",";
    }
    return std::min(forward, backward);
}

}
