/// \file site_traversal_finder.cpp
///
/// Unit tests for the traversal finders that implement both TraversalFinder
/// and SiteTraversalFinder.
///
/// PathTraversalFinder, GBWTTraversalFinder and FlowTraversalFinder each find
/// a site's traversals given either a Snarl, through TraversalFinder, or the
/// site's two bound handles, through SiteTraversalFinder. Both interfaces
/// must give the same traversals for the same site.

#include "catch.hpp"
#include "../traversal_finder.hpp"
#include "../integrated_snarl_finder.hpp"
#include "../vg.hpp"
#include "../gbwt_helper.hpp"
#include "support/random_graph.hpp"
#include "support/randomness.hpp"

#include <format>

namespace vg::unittest {

/// @brief Convert a SnarlTraversal to the graph handles it visits.
static Traversal to_traversal(const HandleGraph& graph, const SnarlTraversal& snarl_traversal) {
    Traversal traversal;
    for (const Visit& visit : snarl_traversal.visit()) {
        traversal.push_back(graph.get_handle(visit.node_id(), visit.backward()));
    }
    return traversal;
}

/**
 * @brief Add a path that walks from a handle along randomly chosen edges.
 *
 * @details The path stops after the given number of steps, or earlier at a
 * dead end.
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

/// @brief Store every path in the graph in a GBWT, in both orientations.
static gbwt::GBWT paths_gbwt(const PathHandleGraph& graph) {
    vector<gbwt::vector_type> paths;
    graph.for_each_path_handle([&graph, &paths](const path_handle_t& path) {
        paths.push_back(extract_as_gbwt_path(graph, graph.get_path_name(path)));
    });
    return get_gbwt(paths);
}

/**
 * @brief Require that a finder gives the same traversals of a site through
 * both of its interfaces.
 *
 * @post Finding traversals of the site's two bound handles, through
 * SiteTraversalFinder, gives the same traversals in the same order as
 * finding traversals of the Snarl, through TraversalFinder.
 * @post Each traversal is nonempty, starts with the site's start handle and
 * ends with its end handle.
 *
 * @return The number of traversals found.
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
 * @brief Require that a FlowTraversalFinder gives the same weighted
 * traversals of a site through both of its interfaces.
 *
 * @post The weighted traversals found from the site's two bound handles are
 * the same, with the same weights in the same order, as those found from the
 * Snarl.
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

/// @brief Make a Snarl holding only a snarl's bounds, read the snarl's way or
/// in reverse.
static Snarl site_of(const Snarl& snarl, bool reversed) {
    Snarl site;
    *site.mutable_start() = reversed ? reverse(snarl.end()) : snarl.start();
    *site.mutable_end() = reversed ? reverse(snarl.start()) : snarl.end();
    return site;
}

TEST_CASE("Traversal finders find the same traversals from a site's bound handles as from its Snarl",
          "[site_traversal_finder]") {

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
