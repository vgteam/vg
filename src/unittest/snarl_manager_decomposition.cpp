/// \file snarl_manager_decomposition.cpp
///
/// Unit tests for SnarlManagerDecomposition, which must present the same
/// snarl tree as a SnarlDistanceIndex built from the same decomposition.

#include "catch.hpp"
#include "../handle.hpp"
#include "../vg.hpp"
#include "../integrated_snarl_finder.hpp"
#include "support/decomposition_bench.hpp"
#include "support/decomposition_comparison.hpp"
#include "support/randomness.hpp"

#include <bdsg/hash_graph.hpp>

#include <string>
#include <vector>

namespace vg::unittest {

using endpoint_t = SnarlDecomposition::endpoint_t;

/// @brief List the children of a net graph element, as for_each_child()
/// gives them.
static vector<net_handle_t> children_of(const SnarlDecomposition& decomposition, const net_handle_t& net) {
    vector<net_handle_t> children;
    decomposition.for_each_child(net, [&children](const net_handle_t& child) {
        children.push_back(child);
    });
    return children;
}

/// @brief Label the children of a net graph element, in the order
/// for_each_child() gives them.
static vector<string> child_labels(const SnarlDecomposition& decomposition, const HandleGraph& graph,
                                   const net_handle_t& net) {
    TraversalLabeler labeler(decomposition, graph);
    vector<string> labels;
    for (const net_handle_t& child : children_of(decomposition, net)) {
        labels.push_back(labeler.label(child));
    }
    return labels;
}

/// @brief Label the realizable traversals of a net graph element, in the
/// order for_each_traversal() gives them.
static vector<string> traversal_labels(const SnarlDecomposition& decomposition, const HandleGraph& graph,
                                       const net_handle_t& net) {
    TraversalLabeler labeler(decomposition, graph);
    vector<string> labels;
    decomposition.for_each_traversal(net, [&labeler, &labels](const net_handle_t& traversal) {
        labels.push_back(labeler.label(traversal));
    });
    return labels;
}

/**
 * @brief Require the SnarlDecomposition interface's laws about one traversal.
 *
 * @post Flipping the traversal twice gives it back.
 * @post Flipping a traversal other than a sentinel swaps where it starts and
 * where it ends.
 * @post The traversal and its flip have the same canonical traversal.
 */
static void require_traversal_laws(const SnarlDecomposition& decomposition, const net_handle_t& traversal) {
    REQUIRE(decomposition.flip(decomposition.flip(traversal)) == traversal);
    if (!decomposition.is_sentinel(traversal)) {
        REQUIRE(decomposition.starts_at(decomposition.flip(traversal)) == decomposition.ends_at(traversal));
        REQUIRE(decomposition.ends_at(decomposition.flip(traversal)) == decomposition.starts_at(traversal));
    }
    REQUIRE(decomposition.canonical(decomposition.flip(traversal)) == decomposition.canonical(traversal));
}

/**
 * @brief Require the SnarlDecomposition interface's laws about a snarl's or
 * chain's bounds.
 *
 * @post A snarl's bounds are sentinels. Each one starts at the endpoint of
 * its bound, keeps it when flipped, and has the snarl, start-to-end, as its
 * parent.
 * @post If the snarl or chain can be crossed start-to-end,
 * get_parent_traversal() of the start bound facing in and the end bound
 * facing out gives the snarl or chain start-to-end.
 */
static void require_bound_laws(const SnarlDecomposition& decomposition, const net_handle_t& net) {
    net_handle_t forward = decomposition.flip(decomposition.flip(net));
    decomposition.for_each_traversal(net, [&decomposition, &forward](const net_handle_t& traversal) {
        if (decomposition.starts_at(traversal) == SnarlDecomposition::START &&
            decomposition.ends_at(traversal) == SnarlDecomposition::END) {
            forward = traversal;
        }
    });
    if (decomposition.starts_at(forward) == SnarlDecomposition::START &&
        decomposition.ends_at(forward) == SnarlDecomposition::END) {
        REQUIRE(decomposition.get_parent_traversal(decomposition.get_bound(net, false, true),
                                                   decomposition.get_bound(net, true, false)) == forward);
    }
    if (!decomposition.is_snarl(net)) {
        return;
    }
    // The interface wants a snarl parent start-to-end.
    forward = decomposition.get_parent(decomposition.get_bound(net, false, true));
    REQUIRE(decomposition.starts_at(forward) == SnarlDecomposition::START);
    REQUIRE(decomposition.ends_at(forward) == SnarlDecomposition::END);
    REQUIRE(decomposition.canonical(forward) == decomposition.canonical(net));
    for (bool get_end : {false, true}) {
        for (bool face_in : {false, true}) {
            net_handle_t bound = decomposition.get_bound(net, get_end, face_in);
            REQUIRE(decomposition.is_sentinel(bound));
            REQUIRE(decomposition.starts_at(bound) == (get_end ? SnarlDecomposition::END : SnarlDecomposition::START));
            REQUIRE(decomposition.starts_at(decomposition.flip(bound)) == decomposition.starts_at(bound));
            REQUIRE(decomposition.get_parent(bound) == forward);
        }
    }
}

/**
 * @brief Require the SnarlDecomposition interface's laws under a net graph
 * element.
 *
 * @details For every child under the element, and every realizable traversal
 * of each, requires the laws of require_traversal_laws() and
 * require_bound_laws(), and these:
 *
 * @post A child of a start-to-end chain has that chain, start-to-end, as its
 * parent.
 * @post A child of a snarl or the root has a parent with the same canonical
 * traversal as the snarl or the root.
 * @post A child's canonical traversal is one of its realizable traversals,
 * if it has any.
 * @post A chain read end-to-start lists the children it lists start-to-end
 * in reverse order, each flipped, except that a chain whose two bounds are
 * the same node still lists that node first.
 */
static void require_interface_laws(const SnarlDecomposition& decomposition, const HandleGraph& graph,
                                   const net_handle_t& parent) {
    for (const net_handle_t& child : children_of(decomposition, parent)) {
        if (decomposition.is_chain(parent)) {
            REQUIRE(decomposition.get_parent(child) == parent);
        } else {
            REQUIRE(decomposition.canonical(decomposition.get_parent(child)) == decomposition.canonical(parent));
        }
        bool has_traversals = false;
        bool canonical_is_realizable = false;
        decomposition.for_each_traversal(child, [&decomposition, &child, &has_traversals,
                                                 &canonical_is_realizable](const net_handle_t& traversal) {
            require_traversal_laws(decomposition, traversal);
            has_traversals = true;
            canonical_is_realizable = canonical_is_realizable || traversal == decomposition.canonical(child);
        });
        REQUIRE(canonical_is_realizable == has_traversals);
        if (decomposition.is_snarl(child) || decomposition.is_chain(child)) {
            require_bound_laws(decomposition, child);
        }
        if (decomposition.is_chain(child)) {
            vector<net_handle_t> expected = children_of(decomposition, child);
            std::ranges::reverse(expected);
            for (net_handle_t& grandchild : expected) {
                grandchild = decomposition.flip(grandchild);
            }
            if (decomposition.get_handle(decomposition.get_bound(child, false, true), &graph) ==
                decomposition.get_handle(decomposition.get_bound(child, true, false), &graph)) {
                std::ranges::rotate(expected, expected.end() - 1);
            }
            REQUIRE(children_of(decomposition, decomposition.flip(child)) == expected);
        }
        if (!decomposition.is_node(child)) {
            require_interface_laws(decomposition, graph, child);
        }
    }
}

/// @brief Require the interface's laws on the view, and that the view
/// presents the same tree as the distance index.
static void require_same_as_distance_index(const DecompositionBench& bench) {
    require_interface_laws(bench.manager_view(), bench.graph(), bench.manager_view().get_root());
    require_same_tree(bench.distance_index(), bench.manager_view(), bench.graph());
}

TEST_CASE("SnarlManagerDecomposition presents the same tree as a distance index",
          "[snarl_manager_decomposition]") {

    SECTION("A node and a chain connected in the root are presented the same way") {
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
        DecompositionBench bench(graph, finder);
        require_same_as_distance_index(bench);
    }

    SECTION("Three nested snarls are presented the same way") {
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
        DecompositionBench bench(graph, finder);
        require_same_as_distance_index(bench);
    }

    SECTION("A snarl with no nodes and an edge that turns around at a bound is hidden in both") {
        VG graph;
        Node* n1 = graph.create_node("GCA");
        Node* n2 = graph.create_node("T");
        Node* n3 = graph.create_node("G");
        graph.create_edge(n1, n2);
        graph.create_edge(n2, n3);
        // Turn around on the end of node 2
        graph.create_edge(n2, n2, false, true);
        IntegratedSnarlFinder finder(graph);
        DecompositionBench bench(graph, finder);
        require_same_as_distance_index(bench);
    }

    SECTION("A cyclic chain is presented the same way") {
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
        DecompositionBench bench(graph, finder);
        require_same_as_distance_index(bench);
    }

    SECTION("A snarl holding a tip node is presented the same way") {
        VG graph;
        Node* n1 = graph.create_node("GCA");
        Node* n2 = graph.create_node("T");
        Node* n3 = graph.create_node("G");
        Node* n4 = graph.create_node("CTGA");
        Node* n5 = graph.create_node("A");
        graph.create_edge(n1, n2);
        graph.create_edge(n1, n3);
        graph.create_edge(n2, n4);
        graph.create_edge(n3, n4);
        // Node 5 hangs off node 3 and leads nowhere.
        graph.create_edge(n3, n5);
        IntegratedSnarlFinder finder(graph);
        DecompositionBench bench(graph, finder);
        require_same_as_distance_index(bench);
    }

    SECTION("Random graphs are presented the same way") {
        std::default_random_engine generator(test_seed_source());
        for (double chain_flip_probability : {0.0, 0.5}) {
            for (size_t repeat = 0; repeat < 100; repeat++) {
                bdsg::HashGraph graph = random_flipped_graph(generator);
                DecompositionBench bench(graph, fuzzed_decomposition(graph, chain_flip_probability, generator));
                require_same_as_distance_index(bench);
            }
        }
    }
}

TEST_CASE("SnarlManagerDecomposition shows each of the manager's snarls that has children, with its bounds",
          "[snarl_manager_decomposition]") {
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
    SnarlManager manager = finder.find_snarls();
    SnarlManagerDecomposition view(manager, graph);

    vector<net_handle_t> root_chains = children_of(view, view.get_root());
    REQUIRE(child_labels(view, graph, view.get_root()) == vector<string> {"chain 1+>8+"});
    REQUIRE(child_labels(view, graph, root_chains[0]) == vector<string> {"node 1+ SE", "snarl 1+>8+", "node 8+ SE"});
    net_handle_t outer = children_of(view, root_chains[0])[1];
    REQUIRE(child_labels(view, graph, outer) == vector<string> {"chain 2+>7+"});
    net_handle_t middle = children_of(view, children_of(view, outer)[0])[1];
    REQUIRE(child_labels(view, graph, middle) == vector<string> {"chain 3+>5+", "chain 6+>6+"});

    manager.for_each_snarl_preorder([&view, &graph](const Snarl* snarl) {
        net_handle_t net = view.get_snarl_net(snarl);
        REQUIRE(view.get_snarl(net) == snarl);
        REQUIRE(view.get_handle(view.get_bound(net, false, true), &graph) ==
                graph.get_handle(snarl->start().node_id(), snarl->start().backward()));
        REQUIRE(view.get_handle(view.get_bound(net, true, false), &graph) ==
                graph.get_handle(snarl->end().node_id(), snarl->end().backward()));
    });
}

TEST_CASE("SnarlManagerDecomposition hides snarls with no children", "[snarl_manager_decomposition]") {
    VG graph;
    Node* n1 = graph.create_node("GCA");
    Node* n2 = graph.create_node("T");
    Node* n3 = graph.create_node("G");
    Node* n4 = graph.create_node("CTGA");
    Node* n5 = graph.create_node("GCA");
    graph.create_edge(n1, n2);
    graph.create_edge(n2, n3);
    graph.create_edge(n2, n4);
    graph.create_edge(n3, n4);
    graph.create_edge(n4, n5);
    IntegratedSnarlFinder finder(graph);
    SnarlManager manager = finder.find_snarls();
    SnarlManagerDecomposition view(manager, graph);

    vector<net_handle_t> root_chains = children_of(view, view.get_root());
    REQUIRE(root_chains.size() == 1);
    REQUIRE(child_labels(view, graph, root_chains[0]) ==
            vector<string> {"node 1+ SE", "node 2+ SE", "snarl 2+>4+", "node 4+ SE", "node 5+ SE"});

    // The snarls from node 1 to node 2 and from node 4 to node 5 hold no
    // nodes, so they are hidden.
    size_t hidden = 0;
    manager.for_each_snarl_preorder([&manager, &view, &graph, &hidden](const Snarl* snarl) {
        if (manager.shallow_contents(snarl, graph, false).first.empty()) {
            REQUIRE_THROWS_AS(view.get_snarl_net(snarl), std::invalid_argument);
            hidden++;
        }
    });
    REQUIRE(hidden == 2);

    // Going right from a chain node reaches its neighbour read the way the
    // walk travels. Going left reaches it read the way the traversal reads.
    vector<net_handle_t> chain_children = children_of(view, root_chains[0]);
    TraversalLabeler labeler(view, graph);
    vector<string> reached;
    view.follow_net_edges(chain_children[0], &graph, false, [&labeler, &reached](const net_handle_t& next) {
        reached.push_back(labeler.label(next));
    });
    REQUIRE(reached == vector<string> {"node 2+ SE"});
    reached.clear();
    view.follow_net_edges(chain_children[1], &graph, true, [&labeler, &reached](const net_handle_t& next) {
        reached.push_back(labeler.label(next));
    });
    REQUIRE(reached == vector<string> {"node 1+ SE"});
}

TEST_CASE("A cyclic chain lists its shared node once and leads back around to it",
          "[snarl_manager_decomposition]") {
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
    DecompositionBench bench(graph, finder);
    const SnarlManagerDecomposition& view = bench.manager_view();

    vector<net_handle_t> chains = children_of(view, view.get_root());
    REQUIRE(chains.size() == 1);
    const net_handle_t& chain = chains[0];
    vector<net_handle_t> children = children_of(view, chain);
    REQUIRE(children.size() == 3);
    REQUIRE(children.front() == view.get_bound(chain, false, true));
    // The end bound is the same node, read the same way.
    handle_t shared = view.get_handle(view.get_bound(chain, false, true), &graph);
    REQUIRE(view.get_handle(view.get_bound(chain, true, false), &graph) == shared);

    // Leaving the last child along the chain reaches the shared node, read
    // the same way.
    vector<handle_t> reached;
    view.follow_net_edges(children.back(), &graph, false, [&view, &graph, &reached](const net_handle_t& next) {
        REQUIRE(view.is_node(next));
        reached.push_back(view.get_handle(next, &graph));
    });
    REQUIRE(reached == vector<handle_t> {shared});

    // Read in reverse, the chain still lists its shared node first.
    vector<net_handle_t> reverse_children = children_of(view, view.flip(chain));
    REQUIRE(reverse_children == vector<net_handle_t> {view.flip(children[0]), view.flip(children[2]),
                                                      view.flip(children[1])});

    // Every edge on the shared node leads into one of the chain's snarls, so
    // nothing outside the chain can be reached from it.
    for (const net_handle_t& traversal : {chain, view.flip(chain)}) {
        for (bool go_left : {false, true}) {
            REQUIRE(view.follow_net_edges(traversal, &graph, go_left, [](const net_handle_t&) {
                return false;
            }));
        }
    }
}

TEST_CASE("A snarl whose start and end are the same node in the same orientation makes a cyclic chain by itself",
          "[snarl_manager_decomposition]") {
    bdsg::HashGraph graph;
    handle_t h1 = graph.create_handle("GCA", 1);
    handle_t h2 = graph.create_handle("T", 2);
    handle_t h3 = graph.create_handle("G", 3);
    graph.create_edge(h1, h2);
    graph.create_edge(h2, h1);
    graph.create_edge(h1, h3);
    graph.create_edge(h3, h1);
    // The snarl reads out of node 1's end and back into node 1's start.
    Snarl loop;
    loop.mutable_start()->set_node_id(1);
    loop.mutable_end()->set_node_id(1);
    loop.set_start_end_reachable(true);
    vector<Snarl> snarls {loop};
    SnarlManager manager(snarls.begin(), snarls.end());
    SnarlManagerDecomposition view(manager, graph);
    require_interface_laws(view, graph, view.get_root());

    net_handle_t net = view.get_snarl_net(manager.top_level_snarls().front());
    REQUIRE(view.get_handle(view.get_bound(net, false, true), &graph) == h1);
    REQUIRE(view.get_handle(view.get_bound(net, true, false), &graph) == h1);
    REQUIRE(child_labels(view, graph, net) == vector<string> {"chain 2+>2+", "chain 3+>3+"});
    net_handle_t chain = view.get_parent(net);
    REQUIRE(child_labels(view, graph, chain) == vector<string> {"node 1+ SE", "snarl 1+>1+"});
}

TEST_CASE("SnarlManagerDecomposition refuses a unary snarl", "[snarl_manager_decomposition]") {
    bdsg::HashGraph graph;
    handle_t h1 = graph.create_handle("GCA", 1);
    handle_t h2 = graph.create_handle("T", 2);
    graph.create_edge(h1, h2);
    graph.create_edge(h2, graph.flip(h1));
    // The snarl reads out of node 1's end and back into the same side.
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

TEST_CASE("get_parent_traversal refuses free node bounds that make no realizable traversal of one trivial chain",
          "[snarl_manager_decomposition]") {
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
    DecompositionBench bench(graph, finder);
    const SnarlManagerDecomposition& view = bench.manager_view();
    net_handle_t node2 = view.get_net(graph.get_handle(n2->id()), &graph);
    net_handle_t node3 = view.get_net(graph.get_handle(n3->id()), &graph);
    REQUIRE(view.is_chain(view.get_parent_traversal(node2, node2)));
    REQUIRE_THROWS_AS(view.get_parent_traversal(node2, node3), std::invalid_argument);
    // A trivial chain cannot be entered and left at the same end.
    REQUIRE_THROWS_AS(view.get_parent_traversal(node2, view.flip(node2)), std::invalid_argument);
}

TEST_CASE("A snarl holding a tip node reached from its start can be traversed from its start to the tip",
          "[snarl_manager_decomposition]") {
    VG graph;
    Node* n1 = graph.create_node("GCA");
    Node* n2 = graph.create_node("T");
    Node* n3 = graph.create_node("G");
    Node* n4 = graph.create_node("CTGA");
    Node* n5 = graph.create_node("A");
    graph.create_edge(n1, n2);
    graph.create_edge(n1, n3);
    graph.create_edge(n2, n4);
    graph.create_edge(n3, n4);
    // Node 5 hangs off node 3, reached going forward from node 1.
    graph.create_edge(n3, n5);
    IntegratedSnarlFinder finder(graph);
    DecompositionBench bench(graph, finder);
    const SnarlManagerDecomposition& view = bench.manager_view();

    vector<net_handle_t> root_chains = children_of(view, view.get_root());
    REQUIRE(root_chains.size() == 1);
    net_handle_t snarl = children_of(view, root_chains[0])[1];
    REQUIRE(view.is_snarl(snarl));
    REQUIRE(traversal_labels(view, graph, snarl) ==
            vector<string> {"snarl 1+>4+", "snarl 1+>tip", "snarl 4->1-", "snarl tip>1-"});
    REQUIRE(traversal_labels(view, graph, root_chains[0]) ==
            vector<string> {"chain 1+>4+", "chain 1+>tip", "chain 4->1-", "chain tip>1-"});

    // Traversals to the tip can be made from the bounds and the children
    // that hold the tip.
    net_handle_t tip_chain = view.get_net(graph.get_handle(n5->id()), &graph);
    tip_chain = view.get_parent(tip_chain);
    TraversalLabeler labeler(view, graph);
    REQUIRE(labeler.label(view.get_parent_traversal(view.get_bound(snarl, false, true), tip_chain)) ==
            "snarl 1+>tip");
    REQUIRE(labeler.label(view.get_parent_traversal(view.get_bound(root_chains[0], false, true), snarl)) ==
            "chain 1+>tip");
    REQUIRE_THROWS_AS(view.get_parent_traversal(view.get_bound(snarl, true, true), tip_chain),
                      std::invalid_argument);

    // A walk that stops at the tip goes no further along the chain.
    view.for_each_traversal(snarl, [&view, &graph](const net_handle_t& traversal) {
        for (bool go_left : {false, true}) {
            bool leaves_through_tip = go_left ? view.starts_at_tip(traversal) : view.ends_at_tip(traversal);
            if (leaves_through_tip) {
                REQUIRE(view.follow_net_edges(traversal, &graph, go_left, [](const net_handle_t&) {
                    return false;
                }));
            }
        }
    });
}

}
