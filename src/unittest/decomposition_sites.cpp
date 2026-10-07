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
#include "support/random_graph.hpp"
#include "support/randomly_flipped_nodes.hpp"
#include "support/randomness.hpp"

#include <bdsg/hash_graph.hpp>

#include <algorithm>
#include <map>

namespace vg {
namespace unittest {

static string handle_string(const HandleGraph& graph, const handle_t& handle) {
    return to_string(graph.get_id(handle)) + (graph.get_is_reverse(handle) ? "-" : "+");
}

/// Describe a site or chain by its bounds, whichever way it is read.
static string site_key(const SnarlDecomposition& decomposition, const HandleGraph& graph, const net_handle_t& net) {
    SiteBounds bounds = site_bounds(decomposition, graph, net);
    return min(handle_string(graph, bounds.start) + "/" + handle_string(graph, bounds.end),
               handle_string(graph, graph.flip(bounds.end)) + "/" + handle_string(graph, graph.flip(bounds.start)));
}

/// Describe every site through the decomposition_sites functions, by site
/// key, ignoring orientation, and check the laws they obey on the way.
static map<string, string> describe_sites(const SnarlDecomposition& decomposition, const HandleGraph& graph) {
    map<string, string> descriptions;
    function<void(const net_handle_t&)> describe = [&](const net_handle_t& site) {
        SiteBounds bounds = site_bounds(decomposition, graph, site);
        SiteBounds flipped = site_bounds(decomposition, graph, decomposition.flip(site));
        REQUIRE(flipped.start == graph.flip(bounds.end));
        REQUIRE(flipped.end == graph.flip(bounds.start));

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

/// Check the functions on the adapter against the manager's own answers.
static void check_against_manager(const DecompositionPair& pair) {
    pair.manager.for_each_snarl_preorder([&](const Snarl* snarl) {
        if (pair.manager.is_trivial(snarl, pair.graph)) {
            return;
        }
        net_handle_t site = pair.adapter.get_snarl_net(snarl);
        SiteBounds bounds = site_bounds(pair.adapter, pair.graph, site);
        REQUIRE(bounds.start == pair.graph.get_handle(snarl->start().node_id(), snarl->start().backward()));
        REQUIRE(bounds.end == pair.graph.get_handle(snarl->end().node_id(), snarl->end().backward()));
        REQUIRE(is_leaf(pair.adapter, site) == pair.manager.all_children_trivial(snarl, pair.graph));
        size_t enclosing = 0;
        for (const Snarl* parent = pair.manager.parent_of(snarl); parent != nullptr; parent = pair.manager.parent_of(parent)) {
            enclosing++;
        }
        REQUIRE(depth(pair.adapter, site) == enclosing);
        net_handle_t parent = parent_site(pair.adapter, site);
        if (pair.manager.parent_of(snarl) == nullptr) {
            REQUIRE(pair.adapter.is_root(parent));
        } else {
            REQUIRE(pair.adapter.get_snarl(parent) == pair.manager.parent_of(snarl));
        }
    });
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
        map<string, string> sites = describe_sites(pair.adapter, graph);
        REQUIRE(sites == describe_sites(pair.distance_index, graph));
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
                DecompositionPair pair(graph, finder);

                check_against_manager(pair);
                map<string, string> adapter_sites = describe_sites(pair.adapter, graph);
                if (!pair.has_known_difference()) {
                    REQUIRE(adapter_sites == describe_sites(pair.distance_index, graph));
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
