/// \file decomposition_sites.cpp
///
/// Unit tests for the decomposition_sites functions, run on a
/// SnarlManagerDecomposition and a SnarlDistanceIndex built from the same
/// decomposition.

#include "catch.hpp"
#include "../decomposition_sites.hpp"
#include "../handle.hpp"
#include "../vg.hpp"
#include "../integrated_snarl_finder.hpp"
#include "support/decomposition_pair.hpp"
#include "support/randomness.hpp"

#include <algorithm>
#include <format>
#include <map>

namespace vg::unittest {

using SiteDescriptions = map<string, string, less<>>;

/// Check that each traversal of a site or chain that does not start or end
/// at a tip has bounds that, read backward, are the bounds of its flip.
static void check_bounds_laws(const SnarlDecomposition& decomposition, const HandleGraph& graph,
                              const net_handle_t& site) {
    decomposition.for_each_traversal(site, [&decomposition, &graph](const net_handle_t& traversal) {
        if (decomposition.starts_at_tip(traversal) || decomposition.ends_at_tip(traversal)) {
            return;
        }
        SiteBounds bounds = site_bounds(decomposition, graph, traversal);
        SiteBounds flipped = site_bounds(decomposition, graph, decomposition.flip(traversal));
        REQUIRE(flipped.start == graph.flip(bounds.end));
        REQUIRE(flipped.end == graph.flip(bounds.start));
    });
}

/// Check the laws a site obeys with its traversals and with its child sites,
/// and then check them for its child sites.
static void check_site_laws(const SnarlDecomposition& decomposition, const HandleGraph& graph,
                            const net_handle_t& site) {
    check_bounds_laws(decomposition, graph, site);
    for (const net_handle_t& chain : child_chains(decomposition, site)) {
        check_bounds_laws(decomposition, graph, chain);
        for (const net_handle_t& child : child_sites(decomposition, chain)) {
            REQUIRE(bounds_key(decomposition, graph, parent_site(decomposition, child)) ==
                    bounds_key(decomposition, graph, site));
            REQUIRE(depth(decomposition, child) == depth(decomposition, site) + 1);
            check_site_laws(decomposition, graph, child);
        }
    }
}

/// Check the laws every site obeys.
static void check_laws(const SnarlDecomposition& decomposition, const HandleGraph& graph) {
    for (const net_handle_t& site : top_level_sites(decomposition)) {
        REQUIRE(depth(decomposition, site) == 0);
        REQUIRE(decomposition.is_root(parent_site(decomposition, site)));
        check_site_laws(decomposition, graph, site);
    }
}

/// Describe a site's depth, whether it is a leaf, its parent, and its child
/// chains with their sites, ignoring orientation.
static string describe_site(const SnarlDecomposition& decomposition, const HandleGraph& graph,
                            const net_handle_t& site) {
    net_handle_t parent = parent_site(decomposition, site);
    string description = std::format("depth {}{} in {}", depth(decomposition, site),
                                is_leaf(decomposition, site) ? " leaf" : "",
                                decomposition.is_root(parent) ? "root" : bounds_key(decomposition, graph, parent));
    vector<string> chains;
    for (const net_handle_t& chain : child_chains(decomposition, site)) {
        vector<string> children;
        for (const net_handle_t& child : child_sites(decomposition, chain)) {
            children.push_back(bounds_key(decomposition, graph, child));
        }
        chains.push_back(std::format("{}({})", bounds_key(decomposition, graph, chain), unoriented_list(children)));
    }
    std::ranges::sort(chains);
    for (const string& chain : chains) {
        description += " " + chain;
    }
    return description;
}

/// Describe a site and every site nested in it, by bounds key.
static void describe_nested_sites(const SnarlDecomposition& decomposition, const HandleGraph& graph,
                                  const net_handle_t& site, SiteDescriptions& descriptions) {
    descriptions[bounds_key(decomposition, graph, site)] = describe_site(decomposition, graph, site);
    for (const net_handle_t& chain : child_chains(decomposition, site)) {
        for (const net_handle_t& child : child_sites(decomposition, chain)) {
            describe_nested_sites(decomposition, graph, child, descriptions);
        }
    }
}

/// Describe every site by bounds key, and list the top-level sites under the
/// key "top level".
static SiteDescriptions describe_sites(const SnarlDecomposition& decomposition, const HandleGraph& graph) {
    SiteDescriptions descriptions;
    vector<string> top_level;
    for (const net_handle_t& site : top_level_sites(decomposition)) {
        top_level.push_back(bounds_key(decomposition, graph, site));
        describe_nested_sites(decomposition, graph, site, descriptions);
    }
    std::ranges::sort(top_level);
    for (const string& key : top_level) {
        descriptions["top level"] += key + " ";
    }
    return descriptions;
}

/// Count the snarls enclosing a snarl.
static size_t ancestor_count(const SnarlManager& manager, const Snarl* snarl) {
    size_t count = 0;
    for (const Snarl* parent = manager.parent_of(snarl); parent != nullptr; parent = manager.parent_of(parent)) {
        count++;
    }
    return count;
}

/// Check the functions on the view for one of the manager's snarls
/// against the manager's own answers.
static void check_site_against_manager(const DecompositionPair& pair, const Snarl* snarl) {
    net_handle_t site = pair.view.get_snarl_net(snarl);
    SiteBounds bounds = site_bounds(pair.view, pair.graph, site);
    REQUIRE(bounds.start == pair.graph.get_handle(snarl->start().node_id(), snarl->start().backward()));
    REQUIRE(bounds.end == pair.graph.get_handle(snarl->end().node_id(), snarl->end().backward()));
    REQUIRE(is_leaf(pair.view, site) == pair.manager.all_children_trivial(snarl, pair.graph));
    REQUIRE(depth(pair.view, site) == ancestor_count(pair.manager, snarl));
    net_handle_t parent = parent_site(pair.view, site);
    if (pair.manager.parent_of(snarl) == nullptr) {
        REQUIRE(pair.view.is_root(parent));
    } else {
        REQUIRE(pair.view.get_snarl(parent) == pair.manager.parent_of(snarl));
    }
}

/// Check the functions on the view against the manager's own answers, for
/// every snarl the view shows.
static void check_against_manager(const DecompositionPair& pair) {
    pair.manager.for_each_snarl_preorder([&pair](const Snarl* snarl) {
        if (!pair.manager.is_trivial(snarl, pair.graph)) {
            check_site_against_manager(pair, snarl);
        }
    });
}

/// Check the functions on the view against the manager, and against the
/// distance index if the pair has no known difference. Return 1 if the
/// distance index was compared and 0 otherwise.
static size_t check_random_decomposition(const DecompositionPair& pair) {
    check_against_manager(pair);
    check_laws(pair.view, pair.graph);
    if (pair.has_known_difference()) {
        return 0;
    }
    check_laws(pair.distance_index, pair.graph);
    REQUIRE(describe_sites(pair.view, pair.graph) == describe_sites(pair.distance_index, pair.graph));
    return 1;
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
        DecompositionPair pair(graph, finder);
        check_against_manager(pair);
        check_laws(pair.view, graph);
        check_laws(pair.distance_index, graph);
        SiteDescriptions sites = describe_sites(pair.view, graph);
        REQUIRE(sites == describe_sites(pair.distance_index, graph));
        REQUIRE(sites.at("top level") == "1+/8+ ");
        REQUIRE(sites.at("1+/8+") == "depth 0 in root 2+/7+(2+/7+,)");
        REQUIRE(sites.at("2+/7+") == "depth 1 in 1+/8+ 3+/5+(3+/5+,)");
        REQUIRE(sites.at("3+/5+") == "depth 2 leaf in 2+/7+");
    }

    SECTION("A snarl that can be turned around in") {
        VG graph;
        Node* n1 = graph.create_node("GCA");
        Node* n2 = graph.create_node("T");
        Node* n3 = graph.create_node("G");
        graph.create_edge(n1, n2);
        graph.create_edge(n2, n3);
        // Turn around on the end of node 2
        graph.create_edge(n2, n2, false, true);
        IntegratedSnarlFinder finder(graph);
        DecompositionPair pair(graph, finder);
        vector<net_handle_t> sites = top_level_sites(pair.view);
        REQUIRE(sites.size() == 1);
        size_t turnarounds = 0;
        pair.view.for_each_traversal(sites.front(), [&pair, &turnarounds](const net_handle_t& traversal) {
            if (pair.view.starts_at(traversal) == pair.view.ends_at(traversal)) {
                turnarounds++;
            }
        });
        REQUIRE(turnarounds > 0);
        check_laws(pair.view, graph);
    }

    SECTION("Random graphs") {
        std::default_random_engine generator(test_seed_source());
        for (double chain_flip_probability : {0.0, 0.5}) {
            size_t compared = 0;
            for (size_t repeat = 0; repeat < 100; repeat++) {
                RandomDecomposition random(generator, chain_flip_probability);
                compared += check_random_decomposition(random.pair);
            }
            // Require that at least 50 of the 100 graphs had no known difference,
            // and so were compared with the distance index.
            REQUIRE(compared >= 50);
        }
    }
}

}
