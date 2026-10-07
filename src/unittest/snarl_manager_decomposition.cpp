/// \file snarl_manager_decomposition.cpp
///
/// Unit tests for SnarlManagerDecomposition, against a SnarlDistanceIndex
/// built from the same decomposition.

#include "catch.hpp"
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
#include <set>

namespace vg {
namespace unittest {

using endpoint_t = handlegraph::SnarlDecomposition::endpoint_t;

static string handle_string(const HandleGraph& graph, const handle_t& handle) {
    return to_string(graph.get_id(handle)) + (graph.get_is_reverse(handle) ? "-" : "+");
}

/// Describe a snarl or chain by its bounds, whichever way it is read.
static string bounds_key(const handlegraph::SnarlDecomposition& decomposition, const HandleGraph& graph,
                         const net_handle_t& net) {
    handle_t start = decomposition.get_handle(decomposition.get_bound(net, false, true), &graph);
    handle_t end = decomposition.get_handle(decomposition.get_bound(net, true, false), &graph);
    return min(handle_string(graph, start) + "/" + handle_string(graph, end),
               handle_string(graph, graph.flip(end)) + "/" + handle_string(graph, graph.flip(start)));
}

/// Describe the tree under a net handle, ignoring orientation: snarls by
/// their bounds, with their child chains in sorted order, and chains by
/// their children in order or reversed, whichever sorts first.
static string tree_description(const handlegraph::SnarlDecomposition& decomposition, const HandleGraph& graph,
                               const net_handle_t& net) {
    if (decomposition.is_node(net)) {
        return "n" + to_string(graph.get_id(decomposition.get_handle(net, &graph)));
    }
    vector<string> children;
    decomposition.for_each_child(net, [&](const net_handle_t& child) {
        children.push_back(tree_description(decomposition, graph, child));
    });
    if (decomposition.is_chain(net)) {
        string forward;
        string backward;
        for (size_t i = 0; i < children.size(); i++) {
            forward += children[i] + ",";
            backward += children[children.size() - 1 - i] + ",";
        }
        return "[" + min(forward, backward) + "]";
    }
    sort(children.begin(), children.end());
    string description = decomposition.is_root(net) ? "root{" : "snarl " + bounds_key(decomposition, graph, net) + "{";
    for (const string& child : children) {
        description += child + ";";
    }
    return description + "}";
}

/// Describe a traversal by the graph handles it is entered and left
/// through, so that traversals of the same thing match across
/// implementations that orient it differently.
static string traversal_key(const handlegraph::SnarlDecomposition& decomposition, const HandleGraph& graph,
                            const net_handle_t& net) {
    if (decomposition.is_root(net)) {
        return "root";
    }
    if (decomposition.is_node(net)) {
        return "node " + handle_string(graph, decomposition.get_handle(net, &graph));
    }
    if (decomposition.is_sentinel(net)) {
        return "bound of " + bounds_key(decomposition, graph, decomposition.get_parent(net)) + " " +
            handle_string(graph, decomposition.get_handle(net, &graph));
    }
    handle_t start = decomposition.get_handle(decomposition.get_bound(net, false, true), &graph);
    handle_t end = decomposition.get_handle(decomposition.get_bound(net, true, false), &graph);
    auto side = [&](endpoint_t endpoint, bool entering) {
        if (endpoint == handlegraph::SnarlDecomposition::TIP) {
            return string("tip");
        }
        bool at_start = (endpoint == handlegraph::SnarlDecomposition::START);
        handle_t handle = at_start ? start : end;
        // The start bound reads in and the end bound reads out.
        return handle_string(graph, entering == at_start ? handle : graph.flip(handle));
    };
    return string(decomposition.is_snarl(net) ? "snarl " : "chain ") + side(decomposition.starts_at(net), true) +
        ">" + side(decomposition.ends_at(net), false);
}

/// Collect every realizable traversal that does not start or end at an
/// internal tip, and every bound sentinel, under a net handle, by traversal
/// key.
static void collect_traversals(const handlegraph::SnarlDecomposition& decomposition, const HandleGraph& graph,
                               const net_handle_t& net, map<string, net_handle_t>& traversals) {
    if (!decomposition.is_root(net)) {
        decomposition.for_each_traversal(net, [&](const net_handle_t& traversal) {
            if (!decomposition.starts_at_tip(traversal) && !decomposition.ends_at_tip(traversal)) {
                traversals[traversal_key(decomposition, graph, traversal)] = traversal;
            }
        });
    }
    if (decomposition.is_snarl(net)) {
        for (bool get_end : {false, true}) {
            for (bool face_in : {false, true}) {
                net_handle_t bound = decomposition.get_bound(net, get_end, face_in);
                traversals[traversal_key(decomposition, graph, bound)] = bound;
            }
        }
    }
    if (!decomposition.is_node(net)) {
        decomposition.for_each_child(net, [&](const net_handle_t& child) {
            collect_traversals(decomposition, graph, child, traversals);
        });
    }
}

/// Check the interface's laws on the adapter: flipping, endpoints, parents
/// of children and of bounds.
static void check_laws(const SnarlManagerDecomposition& adapter, const net_handle_t& parent) {
    adapter.for_each_child(parent, [&](const net_handle_t& child) {
        REQUIRE(adapter.flip(adapter.flip(child)) == child);
        REQUIRE(adapter.starts_at(adapter.flip(child)) == adapter.ends_at(child));
        REQUIRE(adapter.ends_at(adapter.flip(child)) == adapter.starts_at(child));
        REQUIRE(adapter.canonical(adapter.flip(child)) == adapter.canonical(child));
        if (adapter.is_chain(parent)) {
            // Children come out along the chain, so the chain comes back the same way.
            REQUIRE(adapter.get_parent(child) == parent);
        } else {
            REQUIRE(adapter.canonical(adapter.get_parent(child)) == adapter.canonical(parent));
        }
        if (adapter.is_snarl(child)) {
            for (bool get_end : {false, true}) {
                for (bool face_in : {false, true}) {
                    net_handle_t bound = adapter.get_bound(child, get_end, face_in);
                    REQUIRE(adapter.is_sentinel(bound));
                    REQUIRE(adapter.starts_at(bound) == (get_end ? handlegraph::SnarlDecomposition::END
                                                                 : handlegraph::SnarlDecomposition::START));
                    REQUIRE(adapter.starts_at(adapter.flip(bound)) == adapter.starts_at(bound));
                    REQUIRE(adapter.flip(adapter.flip(bound)) == bound);
                    REQUIRE(adapter.canonical(adapter.get_parent(bound)) == adapter.canonical(child));
                }
            }
        }
        if (!adapter.is_node(child)) {
            check_laws(adapter, child);
        }
    });
}

/// Check that every snarl the adapter shows is a nontrivial snarl the
/// manager owns, with the same bounds, and that every such snarl is shown.
static void check_against_manager(const DecompositionPair& pair) {
    const SnarlManagerDecomposition& adapter = pair.adapter;
    size_t shown = 0;
    function<void(const net_handle_t&)> visit = [&](const net_handle_t& net) {
        if (adapter.is_snarl(net)) {
            const Snarl* snarl = adapter.get_snarl(net);
            REQUIRE(!pair.manager.is_trivial(snarl, pair.graph));
            REQUIRE(adapter.get_snarl_net(snarl) == net);
            REQUIRE(adapter.get_handle(adapter.get_bound(net, false, true), &pair.graph) ==
                    pair.graph.get_handle(snarl->start().node_id(), snarl->start().backward()));
            REQUIRE(adapter.get_handle(adapter.get_bound(net, true, false), &pair.graph) ==
                    pair.graph.get_handle(snarl->end().node_id(), snarl->end().backward()));
            shown++;
        }
        if (!adapter.is_node(net)) {
            adapter.for_each_child(net, visit);
        }
    };
    visit(adapter.get_root());
    size_t nontrivial = 0;
    pair.manager.for_each_snarl_preorder([&](const Snarl* snarl) {
        if (!pair.manager.is_trivial(snarl, pair.graph)) {
            nontrivial++;
        }
    });
    REQUIRE(shown == nontrivial);
}

/// Check that walking the net graph, and asking for parents, give the same
/// answers on both implementations, matching traversals by their keys.
/// Parents produced start-to-end may face either way. Going left from a node
/// to the next node in a chain, the distance index reads that node the way
/// the first node is read, where everywhere else both read what they reach
/// the way they travel, so nodes reached going left are compared without
/// orientation.
static void check_net_edges(const DecompositionPair& pair) {
    map<string, net_handle_t> adapter_traversals;
    map<string, net_handle_t> index_traversals;
    collect_traversals(pair.adapter, pair.graph, pair.adapter.get_root(), adapter_traversals);
    collect_traversals(pair.distance_index, pair.graph, pair.distance_index.get_root(), index_traversals);
    vector<string> adapter_keys;
    vector<string> index_keys;
    for (auto& kv : adapter_traversals) {
        adapter_keys.push_back(kv.first);
    }
    for (auto& kv : index_traversals) {
        index_keys.push_back(kv.first);
    }
    REQUIRE(adapter_keys == index_keys);

    auto neighbours = [&](const handlegraph::SnarlDecomposition& decomposition, const net_handle_t& net, bool go_left) {
        set<string> found;
        decomposition.follow_net_edges(net, &pair.graph, go_left, [&](const net_handle_t& next) {
            if (go_left && decomposition.is_node(next)) {
                found.insert("node " + to_string(pair.graph.get_id(decomposition.get_handle(next, &pair.graph))));
            } else {
                found.insert(traversal_key(decomposition, pair.graph, next));
            }
        });
        return found;
    };
    for (auto& kv : adapter_traversals) {
        const net_handle_t& here = kv.second;
        const net_handle_t& there = index_traversals.at(kv.first);
        for (bool go_left : {false, true}) {
            INFO("following " << (go_left ? "left" : "right") << " from " << kv.first);
            REQUIRE(neighbours(pair.adapter, here, go_left) == neighbours(pair.distance_index, there, go_left));
        }
        net_handle_t parent = pair.adapter.get_parent(here);
        string index_parent = traversal_key(pair.distance_index, pair.graph, pair.distance_index.get_parent(there));
        bool directed = pair.adapter.starts_at(here) != pair.adapter.ends_at(here) &&
            !pair.adapter.starts_at_tip(here) && !pair.adapter.ends_at_tip(here);
        if (pair.adapter.is_chain(parent) && directed) {
            REQUIRE(traversal_key(pair.adapter, pair.graph, parent) == index_parent);
        } else {
            REQUIRE((traversal_key(pair.adapter, pair.graph, parent) == index_parent ||
                     traversal_key(pair.adapter, pair.graph, pair.adapter.flip(parent)) == index_parent));
        }
    }
}

/// Run every comparison on one graph.
static void check_pair(const DecompositionPair& pair) {
    check_laws(pair.adapter, pair.adapter.get_root());
    check_against_manager(pair);
    REQUIRE(tree_description(pair.adapter, pair.graph, pair.adapter.get_root()) ==
            tree_description(pair.distance_index, pair.graph, pair.distance_index.get_root()));
    check_net_edges(pair);
}

TEST_CASE("SnarlManagerDecomposition matches a distance index on hand-built graphs",
          "[snarl_manager_decomposition]") {

    SECTION("A node and a chain connected in the root") {
        VG graph;
        Node* n1 = graph.create_node("GCA");
        Node* n2 = graph.create_node("T");
        Node* n3 = graph.create_node("G");
        Node* n4 = graph.create_node("CTGA");
        Node* n5 = graph.create_node("GCA");
        graph.create_edge(n1, n2);
        graph.create_edge(n1, n1);
        graph.create_edge(n2, n3);
        graph.create_edge(n2, n4);
        graph.create_edge(n3, n4);
        graph.create_edge(n4, n5);
        IntegratedSnarlFinder finder(graph);
        DecompositionPair pair(graph, finder);
        check_pair(pair);
        REQUIRE(tree_description(pair.adapter, graph, pair.adapter.get_root()) ==
                "root{[n1,];[n2,snarl 2+/4+{[n3,];},n4,n5,];}");
    }

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
        check_pair(pair);
        REQUIRE(tree_description(pair.adapter, graph, pair.adapter.get_root()) ==
                "root{[n1,snarl 1+/8+{[n2,snarl 2+/7+{[n3,snarl 3+/5+{[n4,];},n5,];[n6,];},n7,];},n8,];}");
    }

    SECTION("A chain with trivial snarls") {
        VG graph;
        Node* n1 = graph.create_node("GCA");
        Node* n2 = graph.create_node("T");
        Node* n3 = graph.create_node("G");
        Node* n4 = graph.create_node("CTGA");
        Node* n5 = graph.create_node("GCA");
        Node* n6 = graph.create_node("T");
        Node* n7 = graph.create_node("G");
        Node* n8 = graph.create_node("CTGA");
        Node* n9 = graph.create_node("GCA");
        Node* n10 = graph.create_node("T");
        graph.create_edge(n1, n2);
        graph.create_edge(n2, n3);
        graph.create_edge(n2, n4);
        graph.create_edge(n3, n5);
        graph.create_edge(n4, n5);
        graph.create_edge(n5, n6);
        graph.create_edge(n6, n7);
        graph.create_edge(n7, n8);
        graph.create_edge(n7, n9);
        graph.create_edge(n8, n9);
        graph.create_edge(n9, n10);
        IntegratedSnarlFinder finder(graph);
        DecompositionPair pair(graph, finder);
        check_pair(pair);
        REQUIRE(tree_description(pair.adapter, graph, pair.adapter.get_root()) ==
                "root{[n1,n2,snarl 2+/5+{[n3,];[n4,];},n5,n6,n7,snarl 7+/9+{[n8,];},n9,n10,];}");
    }
}

TEST_CASE("SnarlManagerDecomposition matches a distance index on random graphs",
          "[snarl_manager_decomposition]") {

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

            if (pair.has_known_difference()) {
                check_laws(pair.adapter, pair.adapter.get_root());
                check_against_manager(pair);
            } else {
                check_pair(pair);
                compared++;
            }
        }
        // Most random graphs have nothing the implementations present differently.
        REQUIRE(compared >= 50);
    }
}

TEST_CASE("SnarlManagerDecomposition presents the known differences", "[snarl_manager_decomposition]") {

    SECTION("A snarl with no nodes that is not an ultrabubble is shown") {
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
        REQUIRE(pair.has_known_difference());
        check_laws(pair.adapter, pair.adapter.get_root());
        check_against_manager(pair);
        REQUIRE(tree_description(pair.adapter, graph, pair.adapter.get_root()) == "root{[n1,n2,snarl 2+/3+{},n3,];}");
        REQUIRE(tree_description(pair.distance_index, graph, pair.distance_index.get_root()) == "root{[n1,n2,n3,];}");
    }

    SECTION("A cyclic chain starts and ends at the same node") {
        VG graph;
        Node* n1 = graph.create_node("GCA");
        Node* n2 = graph.create_node("T");
        Node* n3 = graph.create_node("G");
        Node* n4 = graph.create_node("GG");
        graph.create_edge(n1, n2);
        graph.create_edge(n1, n3);
        graph.create_edge(n2, n4);
        graph.create_edge(n3, n4);
        graph.create_edge(n4, n1);
        IntegratedSnarlFinder finder(graph);
        DecompositionPair pair(graph, finder);
        REQUIRE(pair.has_known_difference());
        check_laws(pair.adapter, pair.adapter.get_root());
        check_against_manager(pair);
        vector<net_handle_t> chains;
        pair.adapter.for_each_child(pair.adapter.get_root(), [&](const net_handle_t& chain) {
            chains.push_back(chain);
        });
        REQUIRE(chains.size() == 1);
        REQUIRE(pair.adapter.get_bound(chains[0], false, true) == pair.adapter.get_bound(chains[0], true, false));
        // The node the chain starts at is listed once, and the walk wraps around.
        vector<net_handle_t> children;
        pair.adapter.for_each_child(chains[0], [&](const net_handle_t& child) {
            children.push_back(child);
        });
        REQUIRE(children.size() == 3);
        REQUIRE(children.front() == pair.adapter.get_bound(chains[0], false, true));
        size_t found = 0;
        pair.adapter.follow_net_edges(children.back(), &graph, false, [&](const net_handle_t& next) {
            REQUIRE(next == children.front());
            found++;
        });
        REQUIRE(found == 1);
    }

    SECTION("A unary snarl has the same node as both bounds") {
        bdsg::HashGraph graph;
        handle_t h1 = graph.create_handle("GCA", 1);
        handle_t h2 = graph.create_handle("T", 2);
        handle_t h3 = graph.create_handle("G", 3);
        graph.create_edge(h1, h2);
        graph.create_edge(h2, h1);
        graph.create_edge(h1, h3);
        graph.create_edge(h3, h1);
        // A cyclic chain of one snarl, from node 1 back around to node 1
        Snarl unary;
        unary.mutable_start()->set_node_id(1);
        unary.mutable_end()->set_node_id(1);
        unary.set_type(UNARY);
        unary.set_start_end_reachable(true);
        vector<Snarl> snarls {unary};
        SnarlManager manager(snarls.begin(), snarls.end());
        SnarlManagerDecomposition adapter(manager, graph);
        REQUIRE(manager.num_snarls() == 1);
        const Snarl* snarl = manager.top_level_snarls().front();
        REQUIRE(snarl->type() == UNARY);
        net_handle_t net = adapter.get_snarl_net(snarl);
        REQUIRE(graph.get_id(adapter.get_handle(adapter.get_bound(net, false, true), &graph)) == 1);
        REQUIRE(graph.get_id(adapter.get_handle(adapter.get_bound(net, true, false), &graph)) == 1);
        check_laws(adapter, adapter.get_root());
        net_handle_t chain = adapter.get_parent(net);
        REQUIRE(adapter.get_bound(chain, false, true) == adapter.get_bound(chain, true, false));
        size_t children = 0;
        adapter.for_each_child(net, [&](const net_handle_t& child) {
            REQUIRE(adapter.is_chain(child));
            children++;
        });
        REQUIRE(children == 2);
    }

    SECTION("Only questions about internal tips throw") {
        VG graph;
        Node* n1 = graph.create_node("GCA");
        Node* n2 = graph.create_node("T");
        Node* n3 = graph.create_node("G");
        graph.create_edge(n1, n2);
        graph.create_edge(n1, n3);
        graph.create_edge(n2, n3);
        IntegratedSnarlFinder finder(graph);
        DecompositionPair pair(graph, finder);
        net_handle_t root = pair.adapter.get_root();
        REQUIRE_THROWS(pair.adapter.for_each_tippy_child(root, [&](const net_handle_t& child) {}));
        REQUIRE_THROWS(pair.adapter.for_each_traversal_start(root, [&](const net_handle_t& child) {}));
    }
}

}
}
