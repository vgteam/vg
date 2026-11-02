/// \file decomposition_sites.cpp
///
/// Unit tests for the decomposition_sites functions, run on each
/// SnarlDecomposition implementation built from the same decomposition.

#include "../../src/unittest/catch.hpp"
#include "support/site_tree_queries.hpp"
#include "../decomposition_sites.hpp"
#include "../../src/handle.hpp"
#include "../../src/vg.hpp"
#include "../../src/integrated_snarl_finder.hpp"
#include "../../src/snarl_manager_decomposition.hpp"
#include "../../src/unittest/support/decomposition_lineup.hpp"
#include "../../src/unittest/support/snarl_decomposition_fuzzer.hpp"
#include "../../src/unittest/support/random_graph.hpp"
#include "../../src/unittest/support/randomly_flipped_nodes.hpp"
#include "../../src/unittest/support/randomness.hpp"

#include <bdsg/hash_graph.hpp>

#include <algorithm>
#include <map>

namespace vg {
namespace unittest {

static string handle_string(const HandleGraph& graph, const handle_t& handle) {
    return to_string(graph.get_id(handle)) + (graph.get_is_reverse(handle) ? "-" : "+");
}

/// Describe a pair of bounds the same way whichever way they are read.
static string bounds_key(const HandleGraph& graph, const handle_t& start, const handle_t& end) {
    return min(handle_string(graph, start) + "/" + handle_string(graph, end),
               handle_string(graph, graph.flip(end)) + "/" + handle_string(graph, graph.flip(start)));
}

/// Describe a site or chain by its bounds, whichever way it is read.
static string site_key(const SnarlDecomposition& decomposition, const HandleGraph& graph, const net_handle_t& net) {
    const auto [start, end] = site_bounds(decomposition, graph, net);
    return bounds_key(graph, start, end);
}

/// Describe every site through the decomposition_sites functions, by site
/// key, ignoring orientation, and check the laws they obey on the way.
static map<string, string> describe_sites(const SnarlDecomposition& decomposition, const HandleGraph& graph) {
    map<string, string> descriptions;
    function<void(const net_handle_t&)> describe = [&](const net_handle_t& site) {
        const auto [start, end] = site_bounds(decomposition, graph, site);
        const auto [flipped_start, flipped_end] =
            site_bounds(decomposition, graph, decomposition.flip(site));
        REQUIRE(flipped_start == graph.flip(end));
        REQUIRE(flipped_end == graph.flip(start));

        net_handle_t parent = parent_site(decomposition, site);
        string description = "depth " + to_string(depth(decomposition, site)) +
            (is_leaf(decomposition, site) ? " leaf" : "") +
            " in " + (decomposition.is_root(parent) ? "root" : site_key(decomposition, graph, parent));
        vector<string> chains;
        for (const net_handle_t& chain : child_chains(decomposition, site)) {
            vector<string> children;
            for (const net_handle_t& child : child_sites(decomposition, chain)) {
                REQUIRE(site_key(decomposition, graph, parent_site(decomposition, child)) ==
                        site_key(decomposition, graph, site));
                REQUIRE(depth(decomposition, child) == depth(decomposition, site) + 1);
                children.push_back(site_key(decomposition, graph, child));
                describe(child);
            }
            string forward;
            string backward;
            for (size_t i = 0; i < children.size(); i++) {
                forward += children[i] + ",";
                backward += children[children.size() - 1 - i] + ",";
            }
            chains.push_back(site_key(decomposition, graph, chain) + "(" + min(forward, backward) + ")");
        }
        REQUIRE(chains.empty() == is_leaf(decomposition, site));
        sort(chains.begin(), chains.end());
        for (const string& chain : chains) {
            description += " " + chain;
        }
        descriptions[site_key(decomposition, graph, site)] = description;
    };
    vector<string> top_level;
    for (const net_handle_t& site : top_level_sites(decomposition)) {
        REQUIRE(depth(decomposition, site) == 0);
        REQUIRE(decomposition.is_root(parent_site(decomposition, site)));
        top_level.push_back(site_key(decomposition, graph, site));
        describe(site);
    }
    sort(top_level.begin(), top_level.end());
    for (const string& key : top_level) {
        descriptions["top level"] += key + " ";
    }
    return descriptions;
}

/// Check the functions on a view of a manager against the manager's own
/// answers.
static void check_against_manager(const SnarlManager& manager, const SnarlManagerDecomposition& view,
                                  const HandleGraph& graph) {
    manager.for_each_snarl_preorder([&](const Snarl* snarl) {
        if (manager.is_trivial(snarl, graph)) {
            return;
        }
        net_handle_t site = view.get_snarl_net(snarl);
        // The view orients the snarl forward along its chain, which may be
        // the other way from the manager's Snarl.
        REQUIRE(site_key(view, graph, site) ==
                bounds_key(graph, graph.get_handle(snarl->start().node_id(), snarl->start().backward()),
                           graph.get_handle(snarl->end().node_id(), snarl->end().backward())));
        REQUIRE(is_leaf(view, site) == manager.all_children_trivial(snarl, graph));
        size_t enclosing = 0;
        for (const Snarl* parent = manager.parent_of(snarl); parent != nullptr; parent = manager.parent_of(parent)) {
            enclosing++;
        }
        REQUIRE(depth(view, site) == enclosing);
        net_handle_t parent = parent_site(view, site);
        if (manager.parent_of(snarl) == nullptr) {
            REQUIRE(view.is_root(parent));
        } else {
            REQUIRE(view.get_snarl(parent) == manager.parent_of(snarl));
        }
    });
}

/// Return true if a manager holds a snarl with no nodes that is not an
/// ultrabubble. SnarlManagerDecomposition shows such a snarl as a site, and
/// SnarlDistanceIndex leaves out every snarl with no children.
static bool has_known_difference(const SnarlManager& manager, const HandleGraph& graph) {
    bool found = false;
    manager.for_each_snarl_preorder([&](const Snarl* snarl) {
        found = found || (snarl->type() != ULTRABUBBLE && manager.is_leaf(snarl) &&
                          manager.shallow_contents(snarl, graph, false).first.empty());
    });
    return found;
}

/**
 * Check the decomposition_sites functions on each implementation built from
 * one decomposition.
 *
 * @return The descriptions of the sites from SnarlManagerDecomposition.
 */
static map<string, string> check_sites(const HandleGraph& graph, const vector<DecompositionEvent>& events,
                                       bool compare_implementations) {
    SnarlManager manager = ReplaySnarlFinder(&graph, events).find_snarls();
    SnarlManagerDecomposition view(manager, graph);
    check_against_manager(manager, view, graph);
    map<string, string> sites = describe_sites(view, graph);
    if (compare_implementations && !has_known_difference(manager, graph)) {
        DecompositionLineup(graph, events).for_each([&](const string& name, const SnarlDecomposition& decomposition) {
            INFO("Implementation: " << name);
            REQUIRE(describe_sites(decomposition, graph) == sites);
        });
    }
    return sites;
}

TEST_CASE("decomposition_sites gives the same sites on both implementations", "[decomposition_sites]") {

    SECTION("Three nested snarls") {
        VG graph;
        Node* n1 = graph.create_node("GCA");
        Node* n2 = graph.create_node("T");
        Node* n3 = graph.create_node("G");
        Node* n4 = graph.create_node("CTGA");
        Node* n5 = graph.create_node("GCA");
        Node* n6 = graph.create_node("T");
        Node* n7 = graph.create_node("G");
        Node* n8 = graph.create_node("CTGA");
        graph.create_edge(n1, n2);
        graph.create_edge(n1, n8);
        graph.create_edge(n2, n3);
        graph.create_edge(n2, n6);
        graph.create_edge(n3, n4);
        graph.create_edge(n3, n5);
        graph.create_edge(n4, n5);
        graph.create_edge(n5, n7);
        graph.create_edge(n6, n7);
        graph.create_edge(n7, n8);
        IntegratedSnarlFinder finder(graph);
        map<string, string> sites = check_sites(graph, capture_events(finder, graph), true);
        REQUIRE(sites.at("top level") == "1+/8+ ");
        REQUIRE(sites.at("1+/8+") == "depth 0 in root 2+/7+(2+/7+,)");
        REQUIRE(sites.at("2+/7+") == "depth 1 in 1+/8+ 3+/5+(3+/5+,)");
        REQUIRE(sites.at("3+/5+") == "depth 2 leaf in 2+/7+");
    }

    SECTION("Random graphs") {
        std::default_random_engine generator(test_seed_source());
        for (double chain_flip_probability : {0.0, 0.5}) {
            size_t compared = 0;
            for (size_t repeat = 0; repeat < 100; repeat++) {
                size_t bases = std::uniform_int_distribution<size_t>(50, 300)(generator);
                size_t variant_bases = std::uniform_int_distribution<size_t>(1, bases / 20)(generator);
                size_t variant_count = std::uniform_int_distribution<size_t>(1, bases / 30)(generator);
                VG base_graph;
                random_graph(bases, variant_bases, variant_count, &base_graph);
                bdsg::HashGraph graph = randomly_flipped_nodes(base_graph, 0.5, generator);
                IntegratedSnarlFinder base_finder(graph);
                SnarlDecompositionFuzzer finder(&graph, &base_finder, chain_flip_probability, generator);
                vector<DecompositionEvent> events = capture_events(finder, graph);
                check_sites(graph, events, true);
                SnarlManager manager = ReplaySnarlFinder(&graph, events).find_snarls();
                if (!has_known_difference(manager, graph)) {
                    compared++;
                }
            }
            // Most random graphs have nothing the implementations present differently.
            REQUIRE(compared >= 50);
        }
    }
}

}
}
