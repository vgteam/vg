/// \file snarl_manager_decomposition.cpp
///
/// Unit tests for SnarlManagerDecomposition, against a SnarlDistanceIndex
/// built from the same decomposition.

#include "catch.hpp"
#include "../handle.hpp"
#include "../vg.hpp"
#include "../integrated_snarl_finder.hpp"
#include "support/decomposition_pair.hpp"
#include "support/randomness.hpp"

#include <bdsg/hash_graph.hpp>

#include <algorithm>
#include <format>
#include <map>
#include <set>

namespace vg::unittest {

using handlegraph::SnarlDecomposition;
using endpoint_t = SnarlDecomposition::endpoint_t;
using TraversalsByKey = map<string, net_handle_t, less<>>;

/// Describe the tree under a net handle, ignoring orientation: snarls by
/// their bounds, with their child chains in sorted order, and chains by
/// their children in order or reversed, whichever sorts first.
static string tree_description(const SnarlDecomposition& decomposition, const HandleGraph& graph,
                               const net_handle_t& net) {
    if (decomposition.is_node(net)) {
        return std::format("n{}", graph.get_id(decomposition.get_handle(net, &graph)));
    }
    vector<string> children;
    decomposition.for_each_child(net, [&decomposition, &graph, &children](const net_handle_t& child) {
        children.push_back(tree_description(decomposition, graph, child));
    });
    if (decomposition.is_chain(net)) {
        return "[" + unoriented_list(children) + "]";
    }
    std::ranges::sort(children);
    string description = decomposition.is_root(net) ? "root{" : "snarl " + bounds_key(decomposition, graph, net) + "{";
    for (const string& child : children) {
        description += child + ";";
    }
    return description + "}";
}

/// Describe where a traversal of a snarl or chain enters or leaves it: the
/// graph handle reading in or out there, or "tip".
static string traversal_side(const SnarlDecomposition& decomposition, const HandleGraph& graph,
                             const net_handle_t& net, endpoint_t endpoint, bool entering) {
    if (endpoint == SnarlDecomposition::TIP) {
        return "tip";
    }
    bool at_start = (endpoint == SnarlDecomposition::START);
    // The start bound reads in and the end bound reads out.
    handle_t handle = decomposition.get_handle(decomposition.get_bound(net, !at_start, at_start), &graph);
    return handle_string(graph, entering == at_start ? handle : graph.flip(handle));
}

/// Describe a traversal by the graph handles it is entered and left
/// through, so that traversals of the same thing match across
/// implementations that orient it differently.
static string traversal_key(const SnarlDecomposition& decomposition, const HandleGraph& graph,
                            const net_handle_t& net) {
    if (decomposition.is_root(net)) {
        return "root";
    }
    if (decomposition.is_node(net)) {
        return "node " + handle_string(graph, decomposition.get_handle(net, &graph));
    }
    if (decomposition.is_sentinel(net)) {
        return std::format("bound of {} {}", bounds_key(decomposition, graph, decomposition.get_parent(net)),
                           handle_string(graph, decomposition.get_handle(net, &graph)));
    }
    return std::format("{} {}>{}", decomposition.is_snarl(net) ? "snarl" : "chain",
                       traversal_side(decomposition, graph, net, decomposition.starts_at(net), true),
                       traversal_side(decomposition, graph, net, decomposition.ends_at(net), false));
}

/// Collect every bound sentinel of a snarl, by traversal key.
static void collect_sentinels(const SnarlDecomposition& decomposition, const HandleGraph& graph,
                              const net_handle_t& snarl, TraversalsByKey& traversals) {
    for (bool get_end : {false, true}) {
        for (bool face_in : {false, true}) {
            net_handle_t bound = decomposition.get_bound(snarl, get_end, face_in);
            traversals[traversal_key(decomposition, graph, bound)] = bound;
        }
    }
}

/// Collect every realizable traversal that does not start or end at an
/// internal tip, and every bound sentinel, under a net handle, by traversal
/// key.
static void collect_traversals(const SnarlDecomposition& decomposition, const HandleGraph& graph,
                               const net_handle_t& net, TraversalsByKey& traversals) {
    if (!decomposition.is_root(net)) {
        decomposition.for_each_traversal(net, [&decomposition, &graph, &traversals](const net_handle_t& traversal) {
            if (!decomposition.starts_at_tip(traversal) && !decomposition.ends_at_tip(traversal)) {
                traversals[traversal_key(decomposition, graph, traversal)] = traversal;
            }
        });
    }
    if (decomposition.is_snarl(net)) {
        collect_sentinels(decomposition, graph, net, traversals);
    }
    if (!decomposition.is_node(net)) {
        decomposition.for_each_child(net, [&decomposition, &graph, &traversals](const net_handle_t& child) {
            collect_traversals(decomposition, graph, child, traversals);
        });
    }
}

/// Check the interface's laws on the bound sentinels of a snarl.
static void check_bound_laws(const SnarlManagerDecomposition& view, const net_handle_t& snarl) {
    for (bool get_end : {false, true}) {
        for (bool face_in : {false, true}) {
            net_handle_t bound = view.get_bound(snarl, get_end, face_in);
            REQUIRE(view.is_sentinel(bound));
            REQUIRE(view.starts_at(bound) == (get_end ? SnarlDecomposition::END : SnarlDecomposition::START));
            REQUIRE(view.starts_at(view.flip(bound)) == view.starts_at(bound));
            REQUIRE(view.flip(view.flip(bound)) == bound);
            // The interface wants a snarl parent start-to-end.
            REQUIRE(view.get_parent(bound) == view.get_snarl_net(view.get_snarl(snarl)));
        }
    }
    REQUIRE(view.get_parent_traversal(view.get_bound(snarl, false, true), view.get_bound(snarl, true, false)) ==
            view.get_snarl_net(view.get_snarl(snarl)));
}

/// Check the interface's laws on a child: flipping, endpoints, and its
/// parent.
static void check_child_laws(const SnarlManagerDecomposition& view, const net_handle_t& parent,
                             const net_handle_t& child) {
    REQUIRE(view.flip(view.flip(child)) == child);
    REQUIRE(view.starts_at(view.flip(child)) == view.ends_at(child));
    REQUIRE(view.ends_at(view.flip(child)) == view.starts_at(child));
    REQUIRE(view.canonical(view.flip(child)) == view.canonical(child));
    if (view.is_chain(parent)) {
        // The chain was listed start-to-end, so its children read along it,
        // and each one's parent is the chain read the same way.
        REQUIRE(view.get_parent(child) == parent);
    } else {
        REQUIRE(view.canonical(view.get_parent(child)) == view.canonical(parent));
    }
}

/// Check the child and bound laws on the view for everything under a net
/// handle.
static void check_laws(const SnarlManagerDecomposition& view, const net_handle_t& parent) {
    view.for_each_child(parent, [&view, &parent](const net_handle_t& child) {
        check_child_laws(view, parent, child);
        if (view.is_snarl(child)) {
            check_bound_laws(view, child);
        }
        if (!view.is_node(child)) {
            check_laws(view, child);
        }
    });
}

/// Count the snarls the view shows under a net handle, checking that each
/// is a nontrivial snarl the manager owns, with the same bounds.
static size_t count_shown_snarls(const DecompositionPair& pair, const net_handle_t& net) {
    const SnarlManagerDecomposition& view = pair.view;
    size_t shown = 0;
    if (view.is_snarl(net)) {
        const Snarl* snarl = view.get_snarl(net);
        REQUIRE(!pair.manager.is_trivial(snarl, pair.graph));
        REQUIRE(view.get_snarl_net(snarl) == net);
        REQUIRE(view.get_handle(view.get_bound(net, false, true), &pair.graph) ==
                pair.graph.get_handle(snarl->start().node_id(), snarl->start().backward()));
        REQUIRE(view.get_handle(view.get_bound(net, true, false), &pair.graph) ==
                pair.graph.get_handle(snarl->end().node_id(), snarl->end().backward()));
        shown++;
    }
    if (!view.is_node(net)) {
        view.for_each_child(net, [&pair, &shown](const net_handle_t& child) {
            shown += count_shown_snarls(pair, child);
        });
    }
    return shown;
}

/// Check that every snarl the view shows is a nontrivial snarl the
/// manager owns, with the same bounds, and that every such snarl is shown.
static void check_against_manager(const DecompositionPair& pair) {
    size_t nontrivial = 0;
    pair.manager.for_each_snarl_preorder([&pair, &nontrivial](const Snarl* snarl) {
        if (!pair.manager.is_trivial(snarl, pair.graph)) {
            nontrivial++;
        }
    });
    REQUIRE(count_shown_snarls(pair, pair.view.get_root()) == nontrivial);
}

/// Describe what follow_net_edges() reaches from a traversal, starting with
/// the given description of where it starts. Nodes reached going left from a
/// node are described without orientation, because the distance index reads
/// them the way the starting node is read, while everything else is read the
/// way we travel.
static string neighbour_report(const SnarlDecomposition& decomposition, const HandleGraph& graph,
                               const string& start_key, const net_handle_t& net, bool go_left) {
    set<string, less<>> found;
    decomposition.follow_net_edges(net, &graph, go_left,
                                   [&decomposition, &graph, &found, &net, go_left](const net_handle_t& next) {
        if (go_left && decomposition.is_node(net) && decomposition.is_node(next)) {
            found.insert(std::format("node {}", graph.get_id(decomposition.get_handle(next, &graph))));
        } else {
            found.insert(traversal_key(decomposition, graph, next));
        }
    });
    string report = std::format("following {} from {}:", go_left ? "left" : "right", start_key);
    for (const string& key : found) {
        report += std::format(" {};", key);
    }
    return report;
}

/// Check that the parent of a traversal is the same on both implementations.
/// The implementations orient a parent the same way only when it is a chain
/// and the child is start-to-end or end-to-start, so other parents may match
/// either way round.
static void check_parent(const DecompositionPair& pair, const net_handle_t& here, const net_handle_t& there) {
    const SnarlManagerDecomposition& view = pair.view;
    net_handle_t parent = view.get_parent(here);
    string index_parent = traversal_key(pair.distance_index, pair.graph, pair.distance_index.get_parent(there));
    bool directed = view.starts_at(here) != view.ends_at(here) &&
        !view.starts_at_tip(here) && !view.ends_at_tip(here);
    if (view.is_chain(parent) && directed) {
        REQUIRE(traversal_key(view, pair.graph, parent) == index_parent);
    } else {
        REQUIRE((traversal_key(view, pair.graph, parent) == index_parent ||
                 traversal_key(view, pair.graph, view.flip(parent)) == index_parent));
    }
}

/// Check that both implementations have the same traversals, that walking
/// the net graph from each gives the same answers, and that each has the
/// same parent, matching traversals by their keys.
static void check_net_edges(const DecompositionPair& pair) {
    TraversalsByKey view_traversals;
    TraversalsByKey index_traversals;
    collect_traversals(pair.view, pair.graph, pair.view.get_root(), view_traversals);
    collect_traversals(pair.distance_index, pair.graph, pair.distance_index.get_root(), index_traversals);
    vector<string> view_keys;
    vector<string> index_keys;
    for (const auto& [key, traversal] : view_traversals) {
        view_keys.push_back(key);
    }
    for (const auto& [key, traversal] : index_traversals) {
        index_keys.push_back(key);
    }
    REQUIRE(view_keys == index_keys);

    for (const auto& [key, here] : view_traversals) {
        const net_handle_t& there = index_traversals.at(key);
        for (bool go_left : {false, true}) {
            REQUIRE(neighbour_report(pair.view, pair.graph, key, here, go_left) ==
                    neighbour_report(pair.distance_index, pair.graph, key, there, go_left));
        }
        check_parent(pair, here, there);
    }
}

/// Run every comparison on one graph.
static void check_pair(const DecompositionPair& pair) {
    check_laws(pair.view, pair.view.get_root());
    check_against_manager(pair);
    REQUIRE(tree_description(pair.view, pair.graph, pair.view.get_root()) ==
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
        REQUIRE(tree_description(pair.view, graph, pair.view.get_root()) ==
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
        REQUIRE(tree_description(pair.view, graph, pair.view.get_root()) ==
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
        REQUIRE(tree_description(pair.view, graph, pair.view.get_root()) ==
                "root{[n1,n2,snarl 2+/5+{[n3,];[n4,];},n5,n6,n7,snarl 7+/9+{[n8,];},n9,n10,];}");
    }
}

TEST_CASE("SnarlManagerDecomposition matches a distance index on random graphs",
          "[snarl_manager_decomposition]") {

    std::default_random_engine generator(test_seed_source());

    for (double chain_flip_probability : {0.0, 0.5}) {
        size_t compared = 0;
        for (size_t repeat = 0; repeat < 100; repeat++) {
            RandomDecomposition random(generator, chain_flip_probability);
            const DecompositionPair& pair = random.pair;

            if (pair.has_known_difference()) {
                check_laws(pair.view, pair.view.get_root());
                check_against_manager(pair);
            } else {
                check_pair(pair);
                compared++;
            }
        }
        // Require that at least 50 of the 100 graphs had no known difference,
        // and so were compared with the distance index.
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
        check_laws(pair.view, pair.view.get_root());
        check_against_manager(pair);
        REQUIRE(tree_description(pair.view, graph, pair.view.get_root()) == "root{[n1,n2,snarl 2+/3+{},n3,];}");
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
        check_laws(pair.view, pair.view.get_root());
        check_against_manager(pair);
        vector<net_handle_t> chains;
        pair.view.for_each_child(pair.view.get_root(), [&chains](const net_handle_t& chain) {
            chains.push_back(chain);
        });
        REQUIRE(chains.size() == 1);
        REQUIRE(pair.view.get_bound(chains[0], false, true) == pair.view.get_bound(chains[0], true, false));
        // The node the chain starts at is listed once, and the walk wraps around.
        vector<net_handle_t> children;
        pair.view.for_each_child(chains[0], [&children](const net_handle_t& child) {
            children.push_back(child);
        });
        REQUIRE(children.size() == 3);
        REQUIRE(children.front() == pair.view.get_bound(chains[0], false, true));
        size_t found = 0;
        pair.view.follow_net_edges(children.back(), &graph, false, [&children, &found](const net_handle_t& next) {
            REQUIRE(next == children.front());
            found++;
        });
        REQUIRE(found == 1);
        // Going left from the first node wraps around to the last node.
        found = 0;
        pair.view.follow_net_edges(children.front(), &graph, true, [&pair, &children, &found](const net_handle_t& next) {
            REQUIRE(next == pair.view.flip(children.back()));
            found++;
        });
        REQUIRE(found == 1);
        // Read backward, the chain still lists its bound node first.
        vector<net_handle_t> backward_children;
        pair.view.for_each_child(pair.view.flip(chains[0]), [&backward_children](const net_handle_t& child) {
            backward_children.push_back(child);
        });
        REQUIRE(backward_children.size() == 3);
        REQUIRE(backward_children[0] == pair.view.flip(children[0]));
        REQUIRE(backward_children[1] == pair.view.flip(children[2]));
        REQUIRE(backward_children[2] == pair.view.flip(children[1]));
        // The chain's bounds give back the chain.
        REQUIRE(pair.view.get_parent_traversal(pair.view.get_bound(chains[0], false, true),
                                               pair.view.get_bound(chains[0], true, false)) == chains[0]);
        // Every edge on the node at a cyclic chain's bounds leads into one of
        // its snarls, so nothing outside the chain can be reached from it.
        for (const net_handle_t& traversal : {chains[0], pair.view.flip(chains[0])}) {
            for (bool go_left : {false, true}) {
                REQUIRE(pair.view.follow_net_edges(traversal, &graph, go_left, [](const net_handle_t&) {
                    return false;
                }));
            }
        }
    }

    SECTION("A snarl from a node around into the same node has that node as both bounds") {
        bdsg::HashGraph graph;
        handle_t h1 = graph.create_handle("GCA", 1);
        handle_t h2 = graph.create_handle("T", 2);
        handle_t h3 = graph.create_handle("G", 3);
        graph.create_edge(h1, h2);
        graph.create_edge(h2, h1);
        graph.create_edge(h1, h3);
        graph.create_edge(h3, h1);
        // A snarl that reads out of node 1 and around back into it, which
        // makes a cyclic chain by itself
        Snarl loop;
        loop.mutable_start()->set_node_id(1);
        loop.mutable_end()->set_node_id(1);
        loop.set_start_end_reachable(true);
        vector<Snarl> snarls {loop};
        SnarlManager manager(snarls.begin(), snarls.end());
        SnarlManagerDecomposition view(manager, graph);
        REQUIRE(manager.num_snarls() == 1);
        const Snarl* snarl = manager.top_level_snarls().front();
        net_handle_t net = view.get_snarl_net(snarl);
        REQUIRE(graph.get_id(view.get_handle(view.get_bound(net, false, true), &graph)) == 1);
        REQUIRE(graph.get_id(view.get_handle(view.get_bound(net, true, false), &graph)) == 1);
        check_laws(view, view.get_root());
        net_handle_t chain = view.get_parent(net);
        REQUIRE(view.get_bound(chain, false, true) == view.get_bound(chain, true, false));
        size_t children = 0;
        view.for_each_child(net, [&view, &children](const net_handle_t& child) {
            REQUIRE(view.is_chain(child));
            children++;
        });
        REQUIRE(children == 2);
    }

    SECTION("A unary snarl is refused") {
        bdsg::HashGraph graph;
        handle_t h1 = graph.create_handle("GCA", 1);
        handle_t h2 = graph.create_handle("T", 2);
        graph.create_edge(h1, h2);
        graph.create_edge(h2, graph.flip(h1));
        // A unary snarl leaves node 1 through its end side and comes back in
        // through the same side.
        Snarl unary;
        unary.mutable_start()->set_node_id(1);
        unary.mutable_end()->set_node_id(1);
        unary.mutable_end()->set_backward(true);
        unary.set_type(UNARY);
        unary.set_start_end_reachable(true);
        vector<Snarl> snarls {unary};
        SnarlManager manager(snarls.begin(), snarls.end());
        auto make_view = [&manager, &graph]() {
            SnarlManagerDecomposition view(manager, graph);
        };
        REQUIRE_THROWS_AS(make_view(), std::invalid_argument);
    }

    SECTION("Traversal bounds in different trivial chains are refused") {
        VG graph;
        Node* n1 = graph.create_node("GCA");
        Node* n2 = graph.create_node("T");
        Node* n3 = graph.create_node("G");
        Node* n4 = graph.create_node("CTGA");
        graph.create_edge(n1, n2);
        graph.create_edge(n1, n3);
        graph.create_edge(n2, n4);
        graph.create_edge(n3, n4);
        IntegratedSnarlFinder finder(graph);
        DecompositionPair pair(graph, finder);
        net_handle_t node2 = pair.view.get_net(graph.get_handle(n2->id()), &graph);
        net_handle_t node3 = pair.view.get_net(graph.get_handle(n3->id()), &graph);
        REQUIRE(pair.view.is_chain(pair.view.get_parent_traversal(node2, node2)));
        REQUIRE_THROWS_AS(pair.view.get_parent_traversal(node2, node3), std::invalid_argument);
    }

    SECTION("Questions about internal tips throw") {
        VG graph;
        Node* n1 = graph.create_node("GCA");
        Node* n2 = graph.create_node("T");
        Node* n3 = graph.create_node("G");
        graph.create_edge(n1, n2);
        graph.create_edge(n1, n3);
        graph.create_edge(n2, n3);
        IntegratedSnarlFinder finder(graph);
        DecompositionPair pair(graph, finder);
        net_handle_t root = pair.view.get_root();
        REQUIRE_THROWS(pair.view.for_each_tippy_child(root, [](const net_handle_t&) {}));
        REQUIRE_THROWS(pair.view.for_each_traversal_start(root, [](const net_handle_t&) {}));
    }
}

}
