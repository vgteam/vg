/// \file site_traversal_finder.cpp
///
/// Unit tests for SiteTraversalFinder: each finder's handle-pair version
/// agrees with its Snarl version.

#include "catch.hpp"
#include "../traversal_finder.hpp"
#include "../integrated_snarl_finder.hpp"
#include "../vg.hpp"
#include "support/random_graph.hpp"
#include "support/randomness.hpp"

#include <gbwt/dynamic_gbwt.h>

#include <format>

namespace vg::unittest {

/// Convert a SnarlTraversal to handles in the graph.
static Traversal to_traversal(const HandleGraph& graph, const SnarlTraversal& snarl_traversal) {
    Traversal traversal;
    for (const Visit& visit : snarl_traversal.visit()) {
        traversal.push_back(graph.get_handle(visit.node_id(), visit.backward()));
    }
    return traversal;
}

/// Add a path that walks from the given handle along random edges, for at
/// most the given number of steps.
template<typename URNG>
static void add_random_walk(MutablePathHandleGraph& graph, const string& name, handle_t here, size_t max_steps,
                            URNG& generator) {
    path_handle_t path = graph.create_path_handle(name);
    graph.append_step(path, here);
    for (size_t step = 1; step < max_steps; step++) {
        vector<handle_t> next;
        graph.follow_edges(here, false, [&next](const handle_t& handle) {
            next.push_back(handle);
        });
        if (next.empty()) {
            break;
        }
        here = next[std::uniform_int_distribution<size_t>(0, next.size() - 1)(generator)];
        graph.append_step(path, here);
    }
}

/// Encode a handle as a GBWT node.
static gbwt::vector_type::value_type gbwt_node(const HandleGraph& graph, const handle_t& handle) {
    return static_cast<gbwt::vector_type::value_type>(gbwt::Node::encode(graph.get_id(handle),
                                                                         graph.get_is_reverse(handle)));
}

/// Store every path in the graph in a GBWT, in both orientations.
static gbwt::GBWT paths_gbwt(const PathHandleGraph& graph) {
    gbwt::vector_type text;
    graph.for_each_path_handle([&graph, &text](const path_handle_t& path) {
        gbwt::vector_type thread;
        graph.for_each_step_in_path(path, [&graph, &thread](const step_handle_t& step) {
            thread.push_back(gbwt_node(graph, graph.get_handle_of_step(step)));
        });
        text.insert(text.end(), thread.begin(), thread.end());
        text.push_back(static_cast<gbwt::vector_type::value_type>(gbwt::ENDMARKER));
        for (auto it = thread.rbegin(); it != thread.rend(); ++it) {
            text.push_back(static_cast<gbwt::vector_type::value_type>(gbwt::Node::reverse(*it)));
        }
        text.push_back(static_cast<gbwt::vector_type::value_type>(gbwt::ENDMARKER));
    });
    gbwt::DynamicGBWT dynamic_gbwt;
    dynamic_gbwt.insert(text, true);
    return gbwt::GBWT(dynamic_gbwt);
}

/// Check that a finder's handle-pair version agrees with its Snarl version
/// on the given site, and return how many traversals it found.
template<typename Finder>
static size_t check_versions_agree(Finder& finder, const HandleGraph& graph, const Snarl& site) {
    handle_t start = graph.get_handle(site.start().node_id(), site.start().backward());
    handle_t end = graph.get_handle(site.end().node_id(), site.end().backward());
    vector<SnarlTraversal> by_snarl = finder.find_traversals(site);
    vector<Traversal> by_handles = finder.find_traversals(start, end);
    REQUIRE(by_handles.size() == by_snarl.size());
    for (size_t i = 0; i < by_handles.size(); i++) {
        REQUIRE(by_handles[i] == to_traversal(graph, by_snarl[i]));
        REQUIRE(!by_handles[i].empty());
        REQUIRE(by_handles[i].front() == start);
        REQUIRE(by_handles[i].back() == end);
    }
    return by_handles.size();
}

/// Check that a FlowTraversalFinder's weighted handle-pair version agrees
/// with its weighted Snarl version on the given site.
static void check_weighted_versions_agree(FlowTraversalFinder& finder, const HandleGraph& graph, const Snarl& site) {
    handle_t start = graph.get_handle(site.start().node_id(), site.start().backward());
    handle_t end = graph.get_handle(site.end().node_id(), site.end().backward());
    auto [snarl_traversals, snarl_weights] = finder.find_weighted_traversals(site);
    auto [handle_traversals, handle_weights] = finder.find_weighted_traversals(start, end);
    REQUIRE(handle_weights == snarl_weights);
    REQUIRE(handle_traversals.size() == snarl_traversals.size());
    for (size_t i = 0; i < handle_traversals.size(); i++) {
        REQUIRE(handle_traversals[i] == to_traversal(graph, snarl_traversals[i]));
    }
}

/// Get a Snarl holding only the bounds of the given snarl, read forward or
/// reversed.
static Snarl site_of(const Snarl& snarl, bool reversed) {
    Snarl site;
    *site.mutable_start() = reversed ? reverse(snarl.end()) : snarl.start();
    *site.mutable_end() = reversed ? reverse(snarl.start()) : snarl.end();
    return site;
}

TEST_CASE("Handle-pair and Snarl versions of each SiteTraversalFinder agree", "[site_traversal_finder]") {

    std::default_random_engine generator(test_seed_source());

    size_t path_found = 0;
    size_t gbwt_found = 0;
    size_t flow_found = 0;
    for (size_t repeat = 0; repeat < 30; repeat++) {
        VG graph;
        small_random_graph(generator, &graph);

        // Add haplotypes beside the reference path that random_graph made.
        handle_t first = graph.get_handle_of_step(graph.path_begin(graph.get_path_handle("path")));
        for (size_t i = 0; i < 3; i++) {
            add_random_walk(graph, std::format("walk{}", i), first, 2 * graph.get_node_count(), generator);
        }
        gbwt::GBWT gbwt_index = paths_gbwt(graph);

        IntegratedSnarlFinder snarl_finder(graph);
        SnarlManager manager = snarl_finder.find_snarls();

        auto node_weight = [&graph](handle_t handle) {
            return static_cast<double>(graph.get_id(handle) % 5 + 1);
        };
        auto edge_weight = [](edge_t) {
            return 1.0;
        };
        PathTraversalFinder path_finder(graph);
        GBWTTraversalFinder gbwt_finder(graph, gbwt_index);
        FlowTraversalFinder flow_finder(graph, 4, node_weight, edge_weight);

        manager.for_each_snarl_preorder([&path_finder, &gbwt_finder, &flow_finder, &graph,
                                         &path_found, &gbwt_found, &flow_found](const Snarl* snarl) {
            for (bool reversed : {false, true}) {
                Snarl site = site_of(*snarl, reversed);
                path_found += check_versions_agree(path_finder, graph, site);
                gbwt_found += check_versions_agree(gbwt_finder, graph, site);
                flow_found += check_versions_agree(flow_finder, graph, site);
                check_weighted_versions_agree(flow_finder, graph, site);
            }
        });
    }
    REQUIRE(path_found > 0);
    REQUIRE(gbwt_found > 0);
    REQUIRE(flow_found > 0);
}

}
