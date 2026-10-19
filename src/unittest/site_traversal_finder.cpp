/// \file site_traversal_finder.cpp
///
/// Unit tests for the traversal finders that implement both TraversalFinder
/// and SiteTraversalFinder.
///
/// PathTraversalFinder, GBWTTraversalFinder and FlowTraversalFinder each find
/// the walks through a site given either a Snarl, through TraversalFinder, or
/// the two bound handles of the site, through SiteTraversalFinder. Both
/// interfaces must give the same walks for the same site.

#include "catch.hpp"
#include "../traversal_finder.hpp"
#include "../integrated_snarl_finder.hpp"
#include "../vg.hpp"
#include "../gbwt_helper.hpp"
#include "support/random_graph.hpp"
#include "support/randomness.hpp"

#include <format>

namespace vg::unittest {

/// Get the handles of the nodes that a SnarlTraversal visits, in order.
static Traversal to_traversal(const HandleGraph& graph, const SnarlTraversal& snarl_traversal) {
    Traversal traversal;
    for (const Visit& visit : snarl_traversal.visit()) {
        traversal.push_back(graph.get_handle(visit.node_id(), visit.backward()));
    }
    return traversal;
}

/**
 * Add a path that walks from a handle along randomly chosen edges.
 *
 * The path stops after max_steps steps, or earlier at a dead end.
 */
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

/// Store each path in the graph in a GBWT, in both orientations.
static gbwt::GBWT paths_gbwt(const PathHandleGraph& graph) {
    vector<gbwt::vector_type> paths;
    graph.for_each_path_handle([&graph, &paths](const path_handle_t& path) {
        paths.push_back(extract_as_gbwt_path(graph, graph.get_path_name(path)));
    });
    return get_gbwt(paths);
}

/**
 * Require that a finder gives the same walks through a site from both of its
 * interfaces.
 *
 * The walks found from the bound handles must match those found from the
 * Snarl, in the same order. Each walk must start with the start handle of
 * the site and end with its end handle.
 *
 * @return The number of walks found.
 */
template<typename Finder>
static size_t require_interfaces_agree(Finder& finder, const HandleGraph& graph, const Snarl& site) {
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

/**
 * Require that a FlowTraversalFinder gives the same weighted walks through a
 * site from both of its interfaces.
 *
 * The walks and flows found from the bound handles must match those found
 * from the Snarl, in the same order.
 */
static void require_weighted_interfaces_agree(FlowTraversalFinder& finder, const HandleGraph& graph, const Snarl& site) {
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

/**
 * Make a Snarl with the bounds of a snarl and no other fields.
 *
 * When reversed is true, the new start is the end of the snarl and the new
 * end is its start, each flipped, so that finders search the snarl in
 * reverse.
 */
static Snarl site_of(const Snarl& snarl, bool reversed) {
    Snarl site;
    *site.mutable_start() = reversed ? reverse(snarl.end()) : snarl.start();
    *site.mutable_end() = reversed ? reverse(snarl.start()) : snarl.end();
    return site;
}

// Each graph is a small random graph with the path that random_graph() makes
// and three random walks as more paths, all also held in a GBWT. Every snarl is searched forward and
// in reverse. Each finder must find the same walks from the bound handles as
// from the Snarl, and must find some walk in some snarl.
TEST_CASE("Traversal finders find the same traversals from a site's bound handles as from its Snarl",
          "[site_traversal_finder]") {

    std::default_random_engine generator(test_seed_source());

    size_t path_found = 0;
    size_t gbwt_found = 0;
    size_t flow_found = 0;
    for (size_t repeat = 0; repeat < 30; repeat++) {
        VG graph;
        small_random_graph(generator, &graph);

        // Start the random walks where the path from random_graph() starts.
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
        // Vary edge weights too, because the width of a walk is its smallest
        // node or edge weight, and equal edge weights would make all widths
        // the same.
        auto edge_weight = [&graph](edge_t edge) {
            return static_cast<double>((graph.get_id(edge.first) + graph.get_id(edge.second)) % 3 + 1);
        };
        PathTraversalFinder path_finder(graph);
        GBWTTraversalFinder gbwt_finder(graph, gbwt_index);
        FlowTraversalFinder flow_finder(graph, 4, node_weight, edge_weight);

        manager.for_each_snarl_preorder([&path_finder, &gbwt_finder, &flow_finder, &graph,
                                         &path_found, &gbwt_found, &flow_found](const Snarl* snarl) {
            for (bool reversed : {false, true}) {
                Snarl site = site_of(*snarl, reversed);
                path_found += require_interfaces_agree(path_finder, graph, site);
                gbwt_found += require_interfaces_agree(gbwt_finder, graph, site);
                flow_found += require_interfaces_agree(flow_finder, graph, site);
                require_weighted_interfaces_agree(flow_finder, graph, site);
            }
        });
    }
    REQUIRE(path_found > 0);
    REQUIRE(gbwt_found > 0);
    REQUIRE(flow_found > 0);
}

}
