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

namespace vg {
namespace unittest {

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
        graph.follow_edges(here, false, [&](const handle_t& handle) {
            next.push_back(handle);
        });
        if (next.empty()) {
            break;
        }
        here = next[std::uniform_int_distribution<size_t>(0, next.size() - 1)(generator)];
        graph.append_step(path, here);
    }
}

TEST_CASE("Handle-pair and Snarl versions of each SiteTraversalFinder agree", "[site_traversal_finder]") {

    std::default_random_engine generator(test_seed_source());

    // How many traversals each finder found, over all graphs
    vector<size_t> found(3, 0);
    for (size_t repeat = 0; repeat < 30; repeat++) {
        size_t bases = std::uniform_int_distribution<size_t>(50, 300)(generator);
        size_t variant_bases = std::uniform_int_distribution<size_t>(1, bases / 20)(generator);
        size_t variant_count = std::uniform_int_distribution<size_t>(1, bases / 30)(generator);
        VG graph;
        random_graph(bases, variant_bases, variant_count, &graph);

        // Add haplotypes beside the reference path that random_graph made.
        handle_t first = graph.get_handle_of_step(graph.path_begin(graph.get_path_handle("path")));
        for (size_t i = 0; i < 3; i++) {
            add_random_walk(graph, "walk" + to_string(i), first, 2 * graph.get_node_count(), generator);
        }

        // Store every path in a GBWT, in both orientations.
        gbwt::vector_type text;
        graph.for_each_path_handle([&](const path_handle_t& path) {
            vector<gbwt::node_type> thread;
            graph.for_each_step_in_path(path, [&](const step_handle_t& step) {
                handle_t handle = graph.get_handle_of_step(step);
                thread.push_back(gbwt::Node::encode(graph.get_id(handle), graph.get_is_reverse(handle)));
            });
            text.insert(text.end(), thread.begin(), thread.end());
            text.push_back(gbwt::ENDMARKER);
            for (auto it = thread.rbegin(); it != thread.rend(); ++it) {
                text.push_back(gbwt::Node::reverse(*it));
            }
            text.push_back(gbwt::ENDMARKER);
        });
        gbwt::DynamicGBWT dynamic_gbwt;
        dynamic_gbwt.insert(text, true);
        gbwt::GBWT gbwt_index(dynamic_gbwt);

        IntegratedSnarlFinder snarl_finder(graph);
        SnarlManager manager = snarl_finder.find_snarls();

        auto node_weight = [&](handle_t handle) {
            return (double) (graph.get_id(handle) % 5 + 1);
        };
        auto edge_weight = [&](edge_t edge) {
            return 1.0;
        };
        PathTraversalFinder path_finder(graph);
        GBWTTraversalFinder gbwt_finder(graph, gbwt_index);
        FlowTraversalFinder flow_finder(graph, 4, node_weight, edge_weight);
        vector<TraversalFinder*> finders {&path_finder, &gbwt_finder, &flow_finder};

        manager.for_each_snarl_preorder([&](const Snarl* snarl) {
            for (bool reversed : {false, true}) {
                // The finders only read the site's bounds.
                Snarl site;
                *site.mutable_start() = reversed ? reverse(snarl->end()) : snarl->start();
                *site.mutable_end() = reversed ? reverse(snarl->start()) : snarl->end();
                handle_t start = graph.get_handle(site.start().node_id(), site.start().backward());
                handle_t end = graph.get_handle(site.end().node_id(), site.end().backward());

                for (size_t i = 0; i < finders.size(); i++) {
                    vector<SnarlTraversal> by_snarl = finders[i]->find_traversals(site);
                    SiteTraversalFinder& site_finder = *finders[i];
                    vector<Traversal> by_handles = site_finder.find_traversals(start, end);
                    REQUIRE(by_handles.size() == by_snarl.size());
                    for (size_t j = 0; j < by_handles.size(); j++) {
                        REQUIRE(by_handles[j] == to_traversal(graph, by_snarl[j]));
                        REQUIRE(!by_handles[j].empty());
                        REQUIRE(by_handles[j].front() == start);
                        REQUIRE(by_handles[j].back() == end);
                    }
                    found[i] += by_handles.size();
                }

                pair<vector<SnarlTraversal>, vector<double>> weighted_by_snarl = flow_finder.find_weighted_traversals(site);
                pair<vector<Traversal>, vector<double>> weighted_by_handles = flow_finder.find_weighted_traversals(start, end);
                REQUIRE(weighted_by_handles.second == weighted_by_snarl.second);
                REQUIRE(weighted_by_handles.first.size() == weighted_by_snarl.first.size());
                for (size_t j = 0; j < weighted_by_handles.first.size(); j++) {
                    REQUIRE(weighted_by_handles.first[j] == to_traversal(graph, weighted_by_snarl.first[j]));
                }
            }
        });
    }
    for (size_t count : found) {
        REQUIRE(count > 0);
    }
}

}
}
