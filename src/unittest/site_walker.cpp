/// \file site_walker.cpp
///
/// Unit tests for SiteWalker and the site values it hands out, run on a SnarlManagerDecomposition
/// and a SnarlDistanceIndex built from the same decomposition.

#include "catch.hpp"
#include "../decomposition_sites.hpp"
#include "../handle.hpp"
#include "../integrated_snarl_finder.hpp"
#include "../site_values.hpp"
#include "../site_walker.hpp"
#include "../vg.hpp"
#include "support/decomposition_pair.hpp"
#include "support/random_graph.hpp"
#include "support/randomly_flipped_nodes.hpp"
#include "support/randomness.hpp"

#include <bdsg/hash_graph.hpp>

#include <algorithm>
#include <map>
#include <mutex>
#include <set>

namespace vg {
namespace unittest {

using namespace vg::multipass;

using std::to_string;
using way_in_t = pair<nid_t, bool>;

static string handle_text(const HandleGraph& graph, const handle_t& handle) {
    return to_string(graph.get_id(handle)) + (graph.get_is_reverse(handle) ? "-" : "+");
}

/// Name a site or chain by its bounds, whichever way round it is read.
static string bounds_key(const HandleGraph& graph, const SiteBounds& bounds) {
    return min(handle_text(graph, bounds.start) + "/" + handle_text(graph, bounds.end),
               handle_text(graph, graph.flip(bounds.end)) + "/"
                   + handle_text(graph, graph.flip(bounds.start)));
}

/// Describe what a site holds, independent of how either implementation orients its chains:
/// each way in to a child site, with the node IDs of the chain it is in, smaller first.
static string describe_children(const HandleGraph& graph, const SiteChildren& children) {
    vector<string> parts;
    for (const auto& entry : children.entries) {
        const ChildChain& chain = children.chains[entry.second];
        nid_t a = graph.get_id(chain.start);
        nid_t b = graph.get_id(chain.end);
        parts.push_back(to_string(entry.first.first) + (entry.first.second ? "-" : "+") + ">"
                        + to_string(min(a, b)) + ":" + to_string(max(a, b)));
    }
    sort(parts.begin(), parts.end());
    string out = children.known ? "known" : "unknown";
    for (const string& part : parts) {
        out += " " + part;
    }
    return out;
}

/// Walk `decomposition` and describe every site visited, by key: how many times it was visited,
/// its depth, the site it is in and what it holds. `answer` is what each call returns.
static map<string, string> walk_sites(const SnarlDecomposition& decomposition,
                                      const HandleGraph& graph,
                                      GraphCaller::RecurseType recurse_type, size_t window,
                                      const function<bool(const SiteView&)>& answer) {
    // What each visit saw, checked once the walk is done, since the calls run on several threads.
    struct Seen {
        string key;
        size_t enclosing;
        size_t depth;
        string enclosing_key;
        string parent_key;
        string children;
    };
    vector<Seen> visits;
    std::mutex lock;
    SiteWalker walker(decomposition, graph);
    walker.walk(recurse_type, window, false, [&](const SiteView& site) {
        const net_handle_t parent = parent_site(decomposition, site.net);
        Seen visit{
            bounds_key(graph, site.bounds), site.enclosing.size(), depth(decomposition, site.net),
            site.enclosing.empty() ? string("root") : bounds_key(graph, site.enclosing[0]),
            decomposition.is_root(parent) ? string("root")
                                          : bounds_key(graph, site_bounds(decomposition, graph,
                                                                          parent)),
            describe_children(graph, site_children(decomposition, graph, site.net, site.bounds))};
        // The view's bounds are the site's.
        {
            lock_guard<std::mutex> guard(lock);
            visits.push_back(std::move(visit));
        }
        return answer(site);
    });
    map<string, string> seen;
    map<string, size_t> times;
    for (const Seen& visit : visits) {
        // The enclosing sites are the site's ancestors.
        REQUIRE(visit.enclosing == visit.depth);
        REQUIRE(visit.enclosing_key == visit.parent_key);
        seen[visit.key] = "depth " + to_string(visit.depth) + " in " + visit.parent_key + ", "
                          + visit.children;
        ++times[visit.key];
    }
    for (auto& kv : seen) {
        kv.second = to_string(times.at(kv.first)) + "x " + kv.second;
    }
    return seen;
}

/// What a walk that queues every site's children must visit: every site the manager does not
/// consider trivial, once each.
static set<string> manager_sites(const DecompositionPair& pair) {
    set<string> sites;
    pair.manager.for_each_snarl_preorder([&](const Snarl* snarl) {
        if (!pair.manager.is_trivial(snarl, pair.graph)) {
            sites.insert(bounds_key(pair.graph, bounds_of(pair.graph, *snarl)));
        }
    });
    return sites;
}

/// The children of `snarl` as the manager's own lookups give them: for each way in to a child
/// snarl that the manager maps back to that child, the chain of the child.
static string manager_children(const DecompositionPair& pair, const Snarl* snarl) {
    SiteChildren children;
    children.known = true;
    for (const Snarl* child : pair.manager.children_of(snarl)) {
        const Chain* chain = pair.manager.chain_of(child);
        const bool in_chain = chain != nullptr && !chain->empty();
        const Visit first = in_chain ? get_start_of(*chain) : child->start();
        const Visit last = in_chain ? get_end_of(*chain) : child->end();
        children.chains.push_back(ChildChain{pair.graph.get_handle(first.node_id(), first.backward()),
                                             pair.graph.get_handle(last.node_id(), last.backward())});
        for (const way_in_t& way_in : {way_in_t(child->start().node_id(), child->start().backward()),
                                     way_in_t(child->end().node_id(), !child->end().backward())}) {
            if (pair.manager.into_which_snarl(way_in.first, way_in.second) == child) {
                children.entries.emplace_back(way_in, children.chains.size() - 1);
            }
        }
    }
    return describe_children(pair.graph, children);
}

/// Check the walk over `pair` in every mode, and against the other implementation unless the
/// decomposition holds something the two present differently.
static void check_walks(const DecompositionPair& pair, size_t window) {
    auto all = [](const SiteView&) { return true; };
    auto none = [](const SiteView&) { return false; };

    map<string, string> adapter_sites = walk_sites(pair.adapter, pair.graph,
                                                   GraphCaller::RecurseAlways, window, all);
    // Every site the manager has, each once.
    set<string> keys;
    for (const auto& kv : adapter_sites) {
        keys.insert(kv.first);
        REQUIRE(kv.second.substr(0, 3) == "1x ");
    }
    REQUIRE(keys == manager_sites(pair));
    // What each site holds, as the manager's own lookups give it.
    pair.manager.for_each_snarl_preorder([&](const Snarl* snarl) {
        if (pair.manager.is_trivial(snarl, pair.graph)) {
            return;
        }
        const string key = bounds_key(pair.graph, bounds_of(pair.graph, *snarl));
        const string& description = adapter_sites.at(key);
        REQUIRE(description.substr(description.find(", ") + 2) == manager_children(pair, snarl));
    });

    // A walk that queues the children of failed calls, where every call fails, visits the same
    // sites; one that queues none visits only the top-level sites.
    REQUIRE(walk_sites(pair.adapter, pair.graph, GraphCaller::RecurseOnFail, window, none)
            == adapter_sites);
    map<string, string> top_level = walk_sites(pair.adapter, pair.graph,
                                               GraphCaller::RecurseNever, window, all);
    for (const auto& kv : top_level) {
        REQUIRE(kv.second.find("depth 0 ") != string::npos);
    }
    size_t roots = 0;
    for (const auto& kv : adapter_sites) {
        roots += kv.second.find("depth 0 ") != string::npos ? 1 : 0;
    }
    REQUIRE(top_level.size() == roots);

    if (!pair.has_known_difference()) {
        REQUIRE(walk_sites(pair.distance_index, pair.graph, GraphCaller::RecurseAlways, window,
                           all) == adapter_sites);
    }
}

/// Describe where each node of the graph lies: the chain it is in and the sites that hold it,
/// innermost first, each named whichever way round it is read.
static map<nid_t, string> node_placements(const SnarlDecomposition& decomposition,
                                          const HandleGraph& graph) {
    map<nid_t, string> out;
    graph.for_each_handle([&](const handle_t& handle) {
        const nid_t id = graph.get_id(handle);
        const NodePlacement place = placement_of_node(decomposition, graph, id);
        string text = "chain " + bounds_key(graph, place.chain) + " in";
        if (place.enclosing.empty()) {
            REQUIRE(decomposition.is_root(place.site));
            text += " root";
        } else {
            // The site handle and the innermost bounds name the same site.
            REQUIRE(!decomposition.is_root(place.site));
            REQUIRE(bounds_key(graph, oriented_bounds(decomposition, graph, place.site))
                    == bounds_key(graph, place.enclosing.front()));
        }
        for (const SiteBounds& site : place.enclosing) {
            text += " " + bounds_key(graph, site);
        }
        out[id] = text;
    });
    return out;
}

TEST_CASE("placement_of_node gives each node's chain and sites on both implementations",
          "[site_walker]") {

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
        const map<nid_t, string> placed = node_placements(pair.adapter, graph);
        REQUIRE(placed == node_placements(pair.distance_index, graph));
        // A node that bounds a site is in that site's chain, held by the chain's parent.
        REQUIRE(placed.at(1) == "chain 1+/8+ in root");
        REQUIRE(placed.at(8) == "chain 1+/8+ in root");
        REQUIRE(placed.at(2) == "chain 2+/7+ in 1+/8+");
        REQUIRE(placed.at(7) == "chain 2+/7+ in 1+/8+");
        REQUIRE(placed.at(3) == "chain 3+/5+ in 2+/7+ 1+/8+");
        REQUIRE(placed.at(5) == "chain 3+/5+ in 2+/7+ 1+/8+");
        // A node that bounds no site is a chain of its own, inside the innermost site.
        REQUIRE(placed.at(4) == "chain 4+/4+ in 3+/5+ 2+/7+ 1+/8+");
        REQUIRE(placed.at(6) == "chain 6+/6+ in 2+/7+ 1+/8+");

        // A chain of one node has a key, and it differs from its neighbours'.
        for (const SnarlDecomposition* decomposition :
             {(const SnarlDecomposition*)&pair.adapter,
              (const SnarlDecomposition*)&pair.distance_index}) {
            const ChildChain four = placement_of_node(*decomposition, graph, 4).chain;
            REQUIRE(graph.get_id(four.start) == 4);
            REQUIRE(graph.get_id(four.end) == 4);
            const size_t expected_key = ((size_t)4 * 1000003) ^ (size_t)4;
            REQUIRE(chain_key_of(graph, four) == expected_key);
            REQUIRE(chain_key_of(graph, four)
                    != chain_key_of(graph, placement_of_node(*decomposition, graph, 3).chain));
            REQUIRE(chain_key_of(graph, placement_of_node(*decomposition, graph, 6).chain)
                    != chain_key_of(graph, placement_of_node(*decomposition, graph, 2).chain));
        }
    }

    SECTION("Random graphs") {
        std::default_random_engine generator(test_seed_source());
        size_t compared = 0;
        size_t nested = 0;
        for (double chain_flip_probability : {0.0, 0.5}) {
            for (size_t repeat = 0; repeat < 50; repeat++) {
                size_t bases = std::uniform_int_distribution<size_t>(50, 300)(generator);
                size_t variant_bases = std::uniform_int_distribution<size_t>(1, bases / 20)(generator);
                size_t variant_count = std::uniform_int_distribution<size_t>(1, bases / 30)(generator);
                VG base_graph;
                random_graph(bases, variant_bases, variant_count, &base_graph);
                bdsg::HashGraph graph = randomly_flipped_nodes(base_graph, 0.5, generator);
                IntegratedSnarlFinder base_finder(graph);
                SnarlDecompositionFuzzer finder(&graph, &base_finder, chain_flip_probability,
                                                generator);
                DecompositionPair pair(graph, finder);
                const map<nid_t, string> placed = node_placements(pair.adapter, graph);
                REQUIRE(placed.size() == graph.get_node_count());
                if (!pair.has_known_difference()) {
                    REQUIRE(placed == node_placements(pair.distance_index, graph));
                    ++compared;
                }
                for (const auto& entry : placed) {
                    nested += entry.second.find(" root") == string::npos ? 1 : 0;
                }
            }
        }
        // The comparison is not vacuous.
        REQUIRE(compared > 10);
        REQUIRE(nested > 0);
    }
}

TEST_CASE("SiteWalker visits the same sites on both implementations", "[site_walker]") {

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
        check_walks(pair, 0);
        check_walks(pair, 4);
        map<string, string> sites = walk_sites(pair.adapter, graph, GraphCaller::RecurseAlways, 0,
                                               [](const SiteView&) { return true; });
        REQUIRE(sites.size() == 3);
        REQUIRE(sites.at("1+/8+") == "1x depth 0 in root, known 2+>2:7 7->2:7");
        REQUIRE(sites.at("2+/7+") == "1x depth 1 in 1+/8+, known 3+>3:5 5->3:5");
        REQUIRE(sites.at("3+/5+") == "1x depth 2 in 2+/7+, known");
    }

    SECTION("Random graphs") {
        std::default_random_engine generator(test_seed_source());
        for (double chain_flip_probability : {0.0, 0.5}) {
            for (size_t repeat = 0; repeat < 50; repeat++) {
                size_t bases = std::uniform_int_distribution<size_t>(50, 300)(generator);
                size_t variant_bases = std::uniform_int_distribution<size_t>(1, bases / 20)(generator);
                size_t variant_count = std::uniform_int_distribution<size_t>(1, bases / 30)(generator);
                VG base_graph;
                random_graph(bases, variant_bases, variant_count, &base_graph);
                bdsg::HashGraph graph = randomly_flipped_nodes(base_graph, 0.5, generator);
                IntegratedSnarlFinder base_finder(graph);
                SnarlDecompositionFuzzer finder(&graph, &base_finder, chain_flip_probability,
                                                generator);
                DecompositionPair pair(graph, finder);
                check_walks(pair, repeat % 2 == 0 ? 0 : 16);
            }
        }
    }
}

}
}
