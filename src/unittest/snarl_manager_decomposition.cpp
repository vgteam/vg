/// \file snarl_manager_decomposition.cpp
///
/// Unit tests for SnarlManagerDecomposition, and for the SnarlDecomposition
/// interface's laws on each implementation.

#include "catch.hpp"
#include "../handle.hpp"
#include "../vg.hpp"
#include "../integrated_snarl_finder.hpp"
#include "../snarl_manager_decomposition.hpp"
#include "support/decomposition_comparison.hpp"
#include "support/decomposition_lineup.hpp"
#include "support/randomness.hpp"

#include <bdsg/hash_graph.hpp>

#include <algorithm>
#include <functional>
#include <string>
#include <vector>

namespace vg::unittest {

using endpoint_t = SnarlDecomposition::endpoint_t;

/// List the children of a net graph element, as for_each_child() gives them.
static vector<net_handle_t> children_of(const SnarlDecomposition& decomposition, const net_handle_t& net) {
    vector<net_handle_t> children;
    decomposition.for_each_child(net, [&children](const net_handle_t& child) {
        children.push_back(child);
    });
    return children;
}

/// Name the children of a net graph element with element_name(), in the
/// order for_each_child() gives them.
static vector<string> child_names(const SnarlDecomposition& decomposition, const HandleGraph& graph,
                                  const net_handle_t& net) {
    vector<string> names;
    for (const net_handle_t& child : children_of(decomposition, net)) {
        names.push_back(element_name(decomposition, graph, child));
    }
    return names;
}

/// Write the realizable traversals of a net graph element, in the order
/// for_each_traversal() gives them, each as two letters for the endpoints
/// where it starts and ends: S for START, E for END and T for TIP.
static vector<string> traversal_kinds(const SnarlDecomposition& decomposition, const net_handle_t& net) {
    auto letter = [](endpoint_t endpoint) {
        return endpoint == SnarlDecomposition::START ? 'S' : endpoint == SnarlDecomposition::END ? 'E' : 'T';
    };
    vector<string> kinds;
    decomposition.for_each_traversal(net, [&](const net_handle_t& traversal) {
        kinds.push_back({letter(decomposition.starts_at(traversal)), letter(decomposition.ends_at(traversal))});
    });
    return kinds;
}

/// Describe a traversal for a failure message: its element's name, or a
/// sentinel's graph handle, and the endpoints where it starts and ends.
static string describe(const SnarlDecomposition& decomposition, const HandleGraph& graph,
                       const net_handle_t& traversal) {
    auto letter = [](endpoint_t endpoint) {
        return endpoint == SnarlDecomposition::START ? 'S' : endpoint == SnarlDecomposition::END ? 'E' : 'T';
    };
    string element;
    if (decomposition.is_sentinel(traversal)) {
        handle_t handle = decomposition.get_handle(traversal, &graph);
        element = "sentinel " + std::to_string(graph.get_id(handle)) + (graph.get_is_reverse(handle) ? "-" : "+");
    } else if (decomposition.is_root(traversal)) {
        element = "root";
    } else {
        element = element_name(decomposition, graph, traversal);
    }
    return element + " " + letter(decomposition.starts_at(traversal)) + letter(decomposition.ends_at(traversal));
}

/// List the traversals that follow_net_edges() gives.
static vector<net_handle_t> neighbours(const SnarlDecomposition& decomposition, const HandleGraph& graph,
                                      const net_handle_t& here, bool go_left) {
    vector<net_handle_t> reached;
    decomposition.follow_net_edges(here, &graph, go_left, [&reached](const net_handle_t& next) {
        reached.push_back(next);
    });
    return reached;
}

/**
 * @brief Require the SnarlDecomposition interface's laws about one traversal.
 *
 * @post Flipping the traversal twice gives it back.
 * @post Flipping a traversal other than a sentinel swaps where it starts and
 * where it ends.
 * @post The traversal and its flip have the same canonical traversal.
 * @post Going left from the traversal gives exactly the traversals from
 * which going right gives it.
 */
static void require_traversal_laws(const SnarlDecomposition& decomposition, const HandleGraph& graph,
                                   const net_handle_t& traversal) {
    REQUIRE(decomposition.flip(decomposition.flip(traversal)) == traversal);
    if (!decomposition.is_sentinel(traversal)) {
        REQUIRE(decomposition.starts_at(decomposition.flip(traversal)) == decomposition.ends_at(traversal));
        REQUIRE(decomposition.ends_at(decomposition.flip(traversal)) == decomposition.starts_at(traversal));
    }
    REQUIRE(decomposition.canonical(decomposition.flip(traversal)) == decomposition.canonical(traversal));
    INFO("Traversal: " << describe(decomposition, graph, traversal));
    for (const net_handle_t& left : neighbours(decomposition, graph, traversal, true)) {
        INFO("Going left reaches: " << describe(decomposition, graph, left));
        vector<net_handle_t> right = neighbours(decomposition, graph, left, false);
        REQUIRE(std::ranges::find(right, traversal) != right.end());
    }
    for (const net_handle_t& right : neighbours(decomposition, graph, traversal, false)) {
        INFO("Going right reaches: " << describe(decomposition, graph, right));
        vector<net_handle_t> left = neighbours(decomposition, graph, right, true);
        REQUIRE(std::ranges::find(left, traversal) != left.end());
    }
}

/**
 * @brief Require the SnarlDecomposition interface's laws about a snarl's or
 * chain's bounds.
 *
 * @post A snarl's bounds are sentinels. Each one starts at the endpoint of
 * its bound, keeps it when flipped, and has the snarl, start-to-end, as its
 * parent.
 * @post A chain's bounds are nodes, each traversed start-to-end or
 * end-to-start.
 * @post If the snarl or chain can be crossed start-to-end,
 * get_parent_traversal() of the start bound facing in and the end bound
 * facing out gives the snarl or chain start-to-end.
 */
static void require_bound_laws(const SnarlDecomposition& decomposition, const net_handle_t& net) {
    vector<string> kinds = traversal_kinds(decomposition, net);
    if (std::ranges::find(kinds, "SE") != kinds.end()) {
        net_handle_t traversal = decomposition.get_parent_traversal(decomposition.get_bound(net, false, true),
                                                                    decomposition.get_bound(net, true, false));
        REQUIRE(decomposition.starts_at(traversal) == SnarlDecomposition::START);
        REQUIRE(decomposition.ends_at(traversal) == SnarlDecomposition::END);
        REQUIRE(decomposition.canonical(traversal) == decomposition.canonical(net));
    }
    for (bool get_end : {false, true}) {
        for (bool face_in : {false, true}) {
            net_handle_t bound = decomposition.get_bound(net, get_end, face_in);
            if (decomposition.is_chain(net)) {
                REQUIRE(decomposition.is_node(bound));
                REQUIRE(decomposition.starts_at(bound) != decomposition.ends_at(bound));
                continue;
            }
            REQUIRE(decomposition.is_sentinel(bound));
            REQUIRE(decomposition.starts_at(bound) == (get_end ? SnarlDecomposition::END : SnarlDecomposition::START));
            REQUIRE(decomposition.starts_at(decomposition.flip(bound)) == decomposition.starts_at(bound));
            net_handle_t parent = decomposition.get_parent(bound);
            REQUIRE(decomposition.starts_at(parent) == SnarlDecomposition::START);
            REQUIRE(decomposition.ends_at(parent) == SnarlDecomposition::END);
            REQUIRE(decomposition.canonical(parent) == decomposition.canonical(net));
        }
    }
}

/**
 * @brief Require the SnarlDecomposition interface's laws under a net graph
 * element.
 *
 * @details For each child under the element, and each realizable traversal
 * of each, requires the laws of require_traversal_laws() and
 * require_bound_laws(), and these:
 *
 * @post A child of a start-to-end chain has that chain, start-to-end, as its
 * parent.
 * @post A child of a snarl or the root has a parent with the same canonical
 * traversal as the snarl or the root.
 * @post A child's canonical traversal is one of its realizable traversals,
 * if it has any.
 * @post A chain traversed end-to-start lists the children it lists
 * start-to-end in reverse order, each flipped, except that a chain whose two
 * bounds are the same node still lists that node first.
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
        decomposition.for_each_traversal(child, [&](const net_handle_t& traversal) {
            require_traversal_laws(decomposition, graph, traversal);
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

/// Make the small graphs that the law and comparison tests share, each named
/// for the structure it has.
static vector<pair<string, std::function<void(VG&)>>> example_graphs() {
    return {
        {"a node and a chain connected in the root", [](VG& graph) {
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
        }},
        {"three nested snarls", [](VG& graph) {
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
        }},
        {"a snarl with no children and an edge that turns around at a bound", [](VG& graph) {
            Node* n1 = graph.create_node("GCA");
            Node* n2 = graph.create_node("T");
            Node* n3 = graph.create_node("G");
            graph.create_edge(n1, n2);
            graph.create_edge(n2, n3);
            graph.create_edge(n2, n2, false, true);
        }},
        {"a cyclic chain", [](VG& graph) {
            Node* n1 = graph.create_node("GCA");
            Node* n2 = graph.create_node("T");
            Node* n3 = graph.create_node("G");
            Node* n4 = graph.create_node("GG");
            graph.create_edge(n1, n2);
            graph.create_edge(n1, n3);
            graph.create_edge(n2, n4);
            graph.create_edge(n3, n4);
            graph.create_edge(n4, n1);
        }},
        {"a snarl holding a tip node", [](VG& graph) {
            Node* n1 = graph.create_node("GCA");
            Node* n2 = graph.create_node("T");
            Node* n3 = graph.create_node("G");
            Node* n4 = graph.create_node("CTGA");
            Node* n5 = graph.create_node("A");
            graph.create_edge(n1, n2);
            graph.create_edge(n1, n3);
            graph.create_edge(n2, n4);
            graph.create_edge(n3, n4);
            graph.create_edge(n3, n5);
        }}
    };
}

/// Require that every implementation in a line-up presents the same tree
/// shape as the first.
static void require_same_tree_shapes(const DecompositionLineup& lineup, const HandleGraph& graph) {
    optional<TreeShape> first;
    lineup.for_each([&](const string& name, const SnarlDecomposition& decomposition) {
        INFO("Implementation: " << name);
        TreeShape shape = tree_shape(decomposition, graph);
        if (!first) {
            first = shape;
        } else {
            REQUIRE(shape == *first);
        }
    });
}

TEST_CASE("Every SnarlDecomposition implementation presents the same tree shape",
          "[snarl_manager_decomposition]") {
    for (const auto& [description, build] : example_graphs()) {
        SECTION("The same shape for " + description) {
            VG graph;
            build(graph);
            IntegratedSnarlFinder finder(graph);
            require_same_tree_shapes(DecompositionLineup(graph, finder), graph);
        }
    }

    SECTION("The same shapes for random graphs") {
        std::default_random_engine generator(test_seed_source());
        for (double chain_flip_probability : {0.0, 0.5}) {
            for (size_t repeat = 0; repeat < 100; repeat++) {
                bdsg::HashGraph graph = random_flipped_graph(generator);
                DecompositionLineup lineup(graph, fuzzed_decomposition(graph, chain_flip_probability, generator));
                require_same_tree_shapes(lineup, graph);
            }
        }
    }
}

/// Check the interface's laws on every implementation of one name in a
/// line-up made for each of the example graphs and some random graphs.
static void require_laws_on_implementation(const string& implementation) {
    auto check = [&implementation](const DecompositionLineup& lineup, const HandleGraph& graph) {
        lineup.for_each([&](const string& name, const SnarlDecomposition& decomposition) {
            if (name == implementation) {
                require_interface_laws(decomposition, graph, decomposition.get_root());
            }
        });
    };
    for (const auto& [description, build] : example_graphs()) {
        INFO("Graph: " << description);
        VG graph;
        build(graph);
        IntegratedSnarlFinder finder(graph);
        check(DecompositionLineup(graph, finder), graph);
    }
    std::default_random_engine generator(test_seed_source());
    for (double chain_flip_probability : {0.0, 0.5}) {
        for (size_t repeat = 0; repeat < 100; repeat++) {
            bdsg::HashGraph graph = random_flipped_graph(generator);
            check(DecompositionLineup(graph, fuzzed_decomposition(graph, chain_flip_probability, generator)), graph);
        }
    }
}

TEST_CASE("SnarlManagerDecomposition keeps the SnarlDecomposition interface's laws",
          "[snarl_manager_decomposition]") {
    require_laws_on_implementation("SnarlManagerDecomposition");
}

TEST_CASE("SnarlDistanceIndex keeps the SnarlDecomposition interface's laws",
          "[snarl_manager_decomposition][snarl_distance_interface][!mayfail]") {
    require_laws_on_implementation("SnarlDistanceIndex");
}

TEST_CASE("SnarlManagerDecomposition shows each of the manager's snarls, oriented forward along its chain",
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
    std::default_random_engine generator(test_seed_source());
    // Reverse some chains, so that some snarls are backward in their chains.
    IntegratedSnarlFinder base_finder(graph);
    SnarlDecompositionFuzzer finder(&graph, &base_finder, 0.5, generator);
    SnarlManager manager = finder.find_snarls();
    SnarlManagerDecomposition view(manager, graph);

    vector<net_handle_t> root_chains = children_of(view, view.get_root());
    REQUIRE(child_names(view, graph, view.get_root()) == vector<string> {"chain 1+/8+"});
    REQUIRE(child_names(view, graph, root_chains[0]).size() == 3);

    manager.for_each_snarl_preorder([&](const Snarl* snarl) {
        net_handle_t net = view.get_snarl_net(snarl);
        REQUIRE(view.get_snarl(net) == snarl);
        bool backward = (*manager.chain_of(snarl))[manager.chain_rank_of(snarl)].second;
        handle_t start = graph.get_handle(snarl->start().node_id(), snarl->start().backward());
        handle_t end = graph.get_handle(snarl->end().node_id(), snarl->end().backward());
        REQUIRE(view.get_handle(view.get_bound(net, false, true), &graph) == (backward ? graph.flip(end) : start));
        REQUIRE(view.get_handle(view.get_bound(net, true, false), &graph) == (backward ? graph.flip(start) : end));
        // The snarl's parent chain, start-to-end, lists the snarl start-to-end.
        vector<net_handle_t> siblings = children_of(view, view.get_parent(net));
        REQUIRE(std::ranges::find(siblings, net) != siblings.end());
    });
}

TEST_CASE("SnarlManagerDecomposition shows a snarl with no children between its bounds",
          "[snarl_manager_decomposition]") {
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
    REQUIRE(child_names(view, graph, root_chains[0]) ==
            vector<string> {"node 1", "snarl 1+/2+", "node 2", "snarl 2+/4+", "node 4", "snarl 4+/5+", "node 5"});
}

TEST_CASE("A chain can turn around at a self-loop in a snarl with no children", "[snarl_manager_decomposition]") {
    // The chain is node 1, the snarl from 1 to 4, node 4, the snarl from 4 to
    // 5, and node 5. The snarl from 4 to 5 has no children, and an edge from
    // the end of 4 back into the end of 4.
    bdsg::HashGraph graph;
    handle_t n1 = graph.create_handle("A");
    handle_t n2 = graph.create_handle("C");
    handle_t n3 = graph.create_handle("G");
    handle_t n4 = graph.create_handle("T");
    handle_t n5 = graph.create_handle("A");
    graph.create_edge(n1, n2);
    graph.create_edge(n1, n3);
    graph.create_edge(n2, n4);
    graph.create_edge(n3, n4);
    graph.create_edge(n4, n5);
    graph.create_edge(n4, graph.flip(n4));
    IntegratedSnarlFinder finder(graph);
    SnarlManager manager = finder.find_snarls();
    SnarlManagerDecomposition view(manager, graph);

    net_handle_t chain = view.get_parent(view.get_net(n1, &graph));
    REQUIRE(view.starts_at(chain) == SnarlDecomposition::START);
    // A walk from the chain's start crosses the first snarl, turns around in
    // the second, and leaves by the chain's start. A walk from the chain's
    // end enters node 4 on the side away from the self-loop, so it cannot
    // turn around.
    REQUIRE(traversal_kinds(view, chain) == vector<string> {"SS", "SE", "ES"});
}

TEST_CASE("A cyclic chain lists its bound node once and leads back around to it",
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
    SnarlManager manager = finder.find_snarls();
    SnarlManagerDecomposition view(manager, graph);

    vector<net_handle_t> chains = children_of(view, view.get_root());
    REQUIRE(chains.size() == 1);
    const net_handle_t& chain = chains[0];
    vector<net_handle_t> children = children_of(view, chain);
    // The chain is node 1, the snarl from 1 to 4, node 4, and the snarl from
    // 4 back to 1, which has no children.
    REQUIRE(children.size() == 4);
    // Both bounds are the first child, traversed the way the chain goes.
    REQUIRE(view.get_bound(chain, false, true) == children.front());
    REQUIRE(view.get_bound(chain, true, false) == children.front());

    // Leaving the last child forward along the chain, a walk next visits the
    // first child.
    REQUIRE(neighbours(view, graph, children.back(), false) == vector<net_handle_t> {children.front()});

    // Traversed in reverse, the chain still lists its bound node first.
    REQUIRE(children_of(view, view.flip(chain)) ==
            vector<net_handle_t> {view.flip(children[0]), view.flip(children[3]), view.flip(children[2]),
                                  view.flip(children[1])});

    // Each edge on the bound node leads into one of the chain's snarls, so
    // nothing outside the chain is next to it.
    for (const net_handle_t& traversal : {chain, view.flip(chain)}) {
        for (bool go_left : {false, true}) {
            REQUIRE(neighbours(view, graph, traversal, go_left).empty());
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
    // The snarl leaves node 1 through its end and comes back in through its
    // start.
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
    REQUIRE(child_names(view, graph, net) == vector<string> {"chain 2+/2+", "chain 3+/3+"});
    REQUIRE(child_names(view, graph, view.get_parent(net)) == vector<string> {"node 1", "snarl 1+/1+"});
}

TEST_CASE("SnarlManagerDecomposition refuses a unary snarl", "[snarl_manager_decomposition]") {
    bdsg::HashGraph graph;
    handle_t h1 = graph.create_handle("GCA", 1);
    handle_t h2 = graph.create_handle("T", 2);
    graph.create_edge(h1, h2);
    graph.create_edge(h2, graph.flip(h1));
    // The snarl leaves node 1 through its end and comes back in through the
    // same side.
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
    SnarlManager manager = finder.find_snarls();
    SnarlManagerDecomposition view(manager, graph);
    net_handle_t node2 = view.get_net(graph.get_handle(n2->id()), &graph);
    net_handle_t node3 = view.get_net(graph.get_handle(n3->id()), &graph);
    REQUIRE(view.is_chain(view.get_parent_traversal(node2, node2)));
    REQUIRE_THROWS_AS(view.get_parent_traversal(node2, node3), std::invalid_argument);
    // A trivial chain cannot be entered and left at the same end.
    REQUIRE_THROWS_AS(view.get_parent_traversal(node2, view.flip(node2)), std::invalid_argument);
}

TEST_CASE("A snarl holding a tip node that its start reaches can be traversed from its start to the tip",
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
    // Node 5 hangs off the end of node 3, and its own end is a dead end.
    graph.create_edge(n3, n5);
    IntegratedSnarlFinder finder(graph);
    SnarlManager manager = finder.find_snarls();
    SnarlManagerDecomposition view(manager, graph);

    vector<net_handle_t> root_chains = children_of(view, view.get_root());
    REQUIRE(root_chains.size() == 1);
    const net_handle_t& chain = root_chains[0];
    net_handle_t snarl = children_of(view, chain)[1];
    REQUIRE(view.is_snarl(snarl));
    REQUIRE(traversal_kinds(view, snarl) == vector<string> {"SE", "ST", "ES", "TS"});
    REQUIRE(traversal_kinds(view, chain) == vector<string> {"SE", "ST", "ES", "TS"});

    SECTION("The snarl's tippy child is node 5's trivial chain, traversed away from the dead end") {
        net_handle_t tip_chain = view.get_parent(view.get_net(graph.get_handle(n5->id(), true), &graph));
        vector<net_handle_t> tippy;
        view.for_each_tippy_child(snarl, [&tippy](const net_handle_t& child) {
            tippy.push_back(child);
        });
        REQUIRE(tippy == vector<net_handle_t> {tip_chain});
        net_handle_t from_tip = view.get_parent_traversal(tip_chain, view.get_bound(snarl, false, false));
        REQUIRE(view.starts_at(from_tip) == SnarlDecomposition::TIP);
        REQUIRE(view.ends_at(from_tip) == SnarlDecomposition::START);
    }

    SECTION("Traversals to the tip can be made from the bounds and the children that hold the tip") {
        net_handle_t tip_chain = view.get_parent(view.get_net(graph.get_handle(n5->id()), &graph));
        net_handle_t to_tip = view.get_parent_traversal(view.get_bound(snarl, false, true), tip_chain);
        REQUIRE(to_tip == view.get_parent_traversal(view.get_bound(snarl, false, true), view.flip(tip_chain)));
        REQUIRE(view.starts_at(to_tip) == SnarlDecomposition::START);
        REQUIRE(view.ends_at(to_tip) == SnarlDecomposition::TIP);
        net_handle_t chain_to_tip = view.get_parent_traversal(view.get_bound(chain, false, true), snarl);
        REQUIRE(view.starts_at(chain_to_tip) == SnarlDecomposition::START);
        REQUIRE(view.ends_at(chain_to_tip) == SnarlDecomposition::TIP);
        REQUIRE_THROWS_AS(view.get_parent_traversal(view.get_bound(snarl, true, true), tip_chain),
                          std::invalid_argument);
    }

    SECTION("A walk that stops at the tip goes no further along the chain") {
        view.for_each_traversal(snarl, [&](const net_handle_t& traversal) {
            if (view.ends_at_tip(traversal)) {
                REQUIRE(neighbours(view, graph, traversal, false).empty());
            }
            if (view.starts_at_tip(traversal)) {
                REQUIRE(neighbours(view, graph, traversal, true).empty());
            }
        });
    }
}

TEST_CASE("A snarl can reach a tip inside a child chain", "[snarl_manager_decomposition]") {
    // The snarl from 1 to 7 holds the chain of 2, the snarl from 2 to 5, and
    // 5. In that inner snarl, node 6 hangs off 3 with nothing after it. The
    // test builds the snarls by hand, because a snarl finder would end the
    // root chain at node 6.
    bdsg::HashGraph graph;
    handle_t n1 = graph.create_handle("A", 1);
    handle_t n2 = graph.create_handle("C", 2);
    handle_t n3 = graph.create_handle("G", 3);
    handle_t n4 = graph.create_handle("T", 4);
    handle_t n5 = graph.create_handle("A", 5);
    handle_t n6 = graph.create_handle("C", 6);
    handle_t n7 = graph.create_handle("G", 7);
    graph.create_edge(n1, n2);
    graph.create_edge(n2, n3);
    graph.create_edge(n2, n4);
    graph.create_edge(n3, n5);
    graph.create_edge(n4, n5);
    graph.create_edge(n3, n6);
    graph.create_edge(n5, n7);
    graph.create_edge(n1, n7);
    Snarl outer;
    outer.mutable_start()->set_node_id(1);
    outer.mutable_end()->set_node_id(7);
    outer.set_start_end_reachable(true);
    Snarl inner;
    inner.mutable_start()->set_node_id(2);
    inner.mutable_end()->set_node_id(5);
    inner.set_start_end_reachable(true);
    *inner.mutable_parent()->mutable_start() = outer.start();
    *inner.mutable_parent()->mutable_end() = outer.end();
    vector<Snarl> snarls {outer, inner};
    SnarlManager manager(snarls.begin(), snarls.end());
    SnarlManagerDecomposition view(manager, graph);

    net_handle_t inner_chain = view.get_parent(view.get_net(n2, &graph));
    REQUIRE(child_names(view, graph, inner_chain) == vector<string> {"node 2", "snarl 2+/5+", "node 5"});
    REQUIRE(traversal_kinds(view, inner_chain) == vector<string> {"SE", "ST", "ES", "TS"});
    net_handle_t outer_snarl = view.get_parent(inner_chain);
    REQUIRE(traversal_kinds(view, outer_snarl) == vector<string> {"SE", "ST", "ES", "TS"});
    require_interface_laws(view, graph, view.get_root());
}

}
