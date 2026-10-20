/// \file snarl_manager_decomposition.cpp
///
/// Unit tests for SnarlManagerDecomposition, and tests that check each
/// SnarlDecomposition implementation against the interface's laws and
/// against the other implementations.

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
#include <set>
#include <functional>
#include <string>
#include <vector>

namespace vg::unittest {

using endpoint_t = SnarlDecomposition::endpoint_t;

/// Name the children of a net graph element or the root with element_name(),
/// in the order for_each_child() gives them.
static vector<string> child_names(const SnarlDecomposition& decomposition, const HandleGraph& graph,
                                  const net_handle_t& net) {
    vector<string> names;
    decomposition.for_each_child(net, [&](const net_handle_t& child) {
        names.push_back(element_name(decomposition, graph, child));
    });
    return names;
}

/// Write an endpoint as S for START, E for END or T for TIP.
static char endpoint_letter(endpoint_t endpoint) {
    return endpoint == SnarlDecomposition::START ? 'S' : endpoint == SnarlDecomposition::END ? 'E' : 'T';
}

/// Write each realizable traversal of a net graph element as the
/// endpoint_letter() of its start and of its end, in the order
/// for_each_traversal() gives them.
static vector<string> traversal_kinds(const SnarlDecomposition& decomposition, const net_handle_t& net) {
    vector<string> kinds;
    decomposition.for_each_traversal(net, [&](const net_handle_t& traversal) {
        kinds.push_back({endpoint_letter(decomposition.starts_at(traversal)),
                         endpoint_letter(decomposition.ends_at(traversal))});
    });
    return kinds;
}

/**
 * Describe a traversal for a failure message.
 *
 * The description is the element_name() of the net graph element, then the
 * endpoint_letter() of the start and of the end of the traversal. A root is
 * named `root`, and a sentinel is named `sentinel` and its graph handle.
 */
static string describe(const SnarlDecomposition& decomposition, const HandleGraph& graph,
                       const net_handle_t& traversal) {
    string element;
    if (decomposition.is_sentinel(traversal)) {
        handle_t handle = decomposition.get_handle(traversal, &graph);
        element = "sentinel " + std::to_string(graph.get_id(handle)) + (graph.get_is_reverse(handle) ? "-" : "+");
    } else if (decomposition.is_root(traversal)) {
        element = "root";
    } else {
        element = element_name(decomposition, graph, traversal);
    }
    return element + " " + endpoint_letter(decomposition.starts_at(traversal)) +
        endpoint_letter(decomposition.ends_at(traversal));
}

/// Require the SnarlDecomposition interface's laws about one traversal.
static void require_traversal_laws(const SnarlDecomposition& decomposition, const HandleGraph& graph,
                                   const net_handle_t& traversal) {
    // Flipping the traversal twice gives it back.
    REQUIRE(decomposition.flip(decomposition.flip(traversal)) == traversal);
    if (!decomposition.is_sentinel(traversal)) {
        // Flipping a traversal of a net graph element swaps where it starts
        // and where it ends.
        REQUIRE(decomposition.starts_at(decomposition.flip(traversal)) == decomposition.ends_at(traversal));
        REQUIRE(decomposition.ends_at(decomposition.flip(traversal)) == decomposition.starts_at(traversal));
    }
    // The traversal and its flip have the same canonical traversal.
    REQUIRE(decomposition.canonical(decomposition.flip(traversal)) == decomposition.canonical(traversal));
    INFO("Traversal: " << describe(decomposition, graph, traversal));
    for (bool go_left : {true, false}) {
        decomposition.follow_net_edges(traversal, &graph, go_left, [&](const net_handle_t& next) {
            INFO("Going " << (go_left ? "left" : "right") << " gives: " << describe(decomposition, graph, next));
            bool comes_back = false;
            decomposition.follow_net_edges(next, &graph, !go_left, [&](const net_handle_t& back) {
                comes_back = comes_back || back == traversal;
            });
            // Going the other way from each traversal that going one way
            // gives leads back to this traversal.
            REQUIRE(comes_back);
        });
    }
}

/**
 * Require the SnarlDecomposition interface's laws about a snarl's or
 * chain's bounds.
 *
 * @param net A traversal of a snarl or chain.
 */
static void require_bound_laws(const SnarlDecomposition& decomposition, const net_handle_t& net) {
    vector<string> kinds = traversal_kinds(decomposition, net);
    if (std::ranges::find(kinds, "SE") != kinds.end()) {
        net_handle_t traversal = decomposition.get_parent_traversal(decomposition.get_bound(net, false, true),
                                                                    decomposition.get_bound(net, true, false));
        // The start bound facing in and the end bound facing out make the
        // snarl or chain START_END.
        REQUIRE(decomposition.starts_at(traversal) == SnarlDecomposition::START);
        REQUIRE(decomposition.ends_at(traversal) == SnarlDecomposition::END);
        REQUIRE(decomposition.canonical(traversal) == decomposition.canonical(net));
    }
    for (bool get_end : {false, true}) {
        for (bool face_in : {false, true}) {
            net_handle_t bound = decomposition.get_bound(net, get_end, face_in);
            if (decomposition.is_chain(net)) {
                // A chain's bounds are node traversals that are START_END or
                // END_START.
                REQUIRE(decomposition.is_node(bound));
                REQUIRE(decomposition.starts_at(bound) != decomposition.ends_at(bound));
                continue;
            }
            // A snarl's bounds are sentinels.
            REQUIRE(decomposition.is_sentinel(bound));
            // A sentinel of the start bound starts at START, and a sentinel of
            // the end bound at END, facing in or out.
            REQUIRE(decomposition.starts_at(bound) == (get_end ? SnarlDecomposition::END : SnarlDecomposition::START));
            REQUIRE(decomposition.starts_at(decomposition.flip(bound)) == decomposition.starts_at(bound));
            net_handle_t parent = decomposition.get_parent(bound);
            // The parent of each sentinel is its snarl START_END.
            REQUIRE(decomposition.starts_at(parent) == SnarlDecomposition::START);
            REQUIRE(decomposition.ends_at(parent) == SnarlDecomposition::END);
            REQUIRE(decomposition.canonical(parent) == decomposition.canonical(net));
        }
    }
}

/**
 * Require the SnarlDecomposition interface's laws on everything under a snarl,
 * chain or root.
 *
 * Check each child, and everything under each child that is a snarl or chain.
 *
 * @param parent The root, or a traversal of a snarl or chain.
 */
static void require_interface_laws(const SnarlDecomposition& decomposition, const HandleGraph& graph,
                                   const net_handle_t& parent) {
    decomposition.for_each_child(parent, [&](const net_handle_t& child) {
        if (decomposition.is_chain(parent)) {
            // A chain's child has as its parent the chain traversal that
            // listed it.
            REQUIRE(decomposition.get_parent(child) == parent);
        } else {
            // A snarl's or the root's child has as its parent some traversal
            // of the snarl or root.
            REQUIRE(decomposition.canonical(decomposition.get_parent(child)) == decomposition.canonical(parent));
        }
        bool has_traversals = false;
        bool canonical_is_realizable = false;
        decomposition.for_each_traversal(child, [&](const net_handle_t& traversal) {
            require_traversal_laws(decomposition, graph, traversal);
            has_traversals = true;
            canonical_is_realizable = canonical_is_realizable || traversal == decomposition.canonical(child);
        });
        // The canonical traversal is realizable, if any traversal is.
        REQUIRE(canonical_is_realizable == has_traversals);
        if (decomposition.is_snarl(child) || decomposition.is_chain(child)) {
            require_bound_laws(decomposition, child);
        }
        if (decomposition.is_chain(child)) {
            vector<net_handle_t> expected;
            decomposition.for_each_child(child, [&](const net_handle_t& grandchild) {
                expected.push_back(decomposition.flip(grandchild));
            });
            std::ranges::reverse(expected);
            if (decomposition.get_handle(decomposition.get_bound(child, false, true), &graph) ==
                decomposition.get_handle(decomposition.get_bound(child, true, false), &graph)) {
                std::ranges::rotate(expected, expected.end() - 1);
            }
            vector<net_handle_t> flipped_children;
            decomposition.for_each_child(decomposition.flip(child), [&](const net_handle_t& grandchild) {
                flipped_children.push_back(grandchild);
            });
            // The flipped chain lists the chain's children in reverse order,
            // each flipped. A looping chain still lists its boundary node
            // first.
            REQUIRE(flipped_children == expected);
        }
        if (!decomposition.is_node(child)) {
            require_interface_laws(decomposition, graph, child);
        }
    });
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
        {"a looping chain", [](VG& graph) {
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
        {"a snarl with a tip", [](VG& graph) {
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
        }},
        {"a snarl with a walk from one tip to another", [](VG& graph) {
            // The end of 3 and the end of 5 both enter the end of 6, so a walk
            // from the start of 5 can go through 6 to the start of 6.
            Node* n1 = graph.create_node("GCA");
            Node* n2 = graph.create_node("T");
            Node* n3 = graph.create_node("G");
            Node* n4 = graph.create_node("CTGA");
            Node* n5 = graph.create_node("A");
            Node* n6 = graph.create_node("C");
            graph.create_edge(n1, n2);
            graph.create_edge(n1, n3);
            graph.create_edge(n2, n4);
            graph.create_edge(n3, n4);
            graph.create_edge(n3, n6, false, true);
            graph.create_edge(n5, n6, false, true);
        }},
        {"a looping chain with a turnaround and a tip on the far side", [](VG& graph) {
            // A walk from the start of 1 can turn around through the end of 2
            // and go around the loop, entering the snarl from 1 to 4 at its
            // end. Only that way can it reach the tip at the start of 5.
            Node* n1 = graph.create_node("GCA");
            Node* n2 = graph.create_node("T");
            Node* n3 = graph.create_node("G");
            Node* n4 = graph.create_node("GG");
            Node* n5 = graph.create_node("A");
            graph.create_edge(n1, n2);
            graph.create_edge(n1, n3);
            graph.create_edge(n2, n4);
            graph.create_edge(n3, n4);
            graph.create_edge(n4, n1);
            graph.create_edge(n1, n2, false, true);
            graph.create_edge(n5, n3);
        }},
        {"a chain that a walk cannot cross with a tip past it", [](VG& graph) {
            // Inside the snarl from 1 to 7, the chain from 2 to 4 holds a snarl
            // that no walk crosses: 3 hangs off 2 and 5 hangs off 4. Node 6
            // hangs off the end of 4, and the end of 6 is a tip.
            Node* n1 = graph.create_node("GCA");
            Node* n2 = graph.create_node("T");
            Node* n3 = graph.create_node("G");
            Node* n4 = graph.create_node("CTGA");
            Node* n5 = graph.create_node("A");
            Node* n6 = graph.create_node("C");
            Node* n7 = graph.create_node("G");
            graph.create_edge(n1, n2);
            graph.create_edge(n2, n3);
            graph.create_edge(n3, n3, false, true);
            graph.create_edge(n5, n4);
            graph.create_edge(n5, n5, true, false);
            graph.create_edge(n4, n6);
            graph.create_edge(n4, n7);
            graph.create_edge(n1, n7);
        }}
    };
}

/// Require that tree_shape() of each implementation in a line-up is the
/// tree_shape() of the first.
static void require_same_tree_shapes(const DecompositionLineup& lineup, const HandleGraph& graph) {
    optional<TreeShape> first;
    lineup.for_each([&](const string& name, const SnarlDecomposition& decomposition) {
        INFO("Implementation: " << name);
        TreeShape shape = tree_shape(decomposition, graph);
        if (!first) {
            first = shape;
            return;
        }
        // Compare entry by entry, so that a failure shows which differs.
        set<string> keys;
        for (const auto& [key, value] : shape) {
            keys.insert(key);
        }
        for (const auto& [key, value] : *first) {
            keys.insert(key);
        }
        for (const string& key : keys) {
            INFO("Entry: " << key);
            REQUIRE(shape[key] == (*first)[key]);
        }
    });
}

TEST_CASE("tree_shape() is the same for each SnarlDecomposition implementation",
          "[snarl_manager_decomposition]") {
    for (const auto& [description, build] : example_graphs()) {
        SECTION("tree_shape() is the same for " + description) {
            VG graph;
            build(graph);
            IntegratedSnarlFinder finder(graph);
            require_same_tree_shapes(DecompositionLineup(graph, finder), graph);
        }
    }

    SECTION("tree_shape() is the same for random graphs") {
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

/**
 * Require the SnarlDecomposition interface's laws on one implementation.
 *
 * Each implementation gets a test case of its own, so a failure names the
 * implementation that broke a law. The example graphs hold the structures
 * that the laws treat specially, such as looping chains and tips. The random
 * graphs have nodes and snarls reversed in parent.
 *
 * @param implementation The class name that DecompositionLineup::for_each()
 * gives the implementation.
 */
static void require_laws_on_implementation(const string& implementation) {
    auto check = [&implementation](const DecompositionLineup& lineup, const HandleGraph& graph) {
        bool found = false;
        lineup.for_each([&](const string& name, const SnarlDecomposition& decomposition) {
            if (name == implementation) {
                found = true;
                require_interface_laws(decomposition, graph, decomposition.get_root());
            }
        });
        // The line-up has an implementation with the name.
        REQUIRE(found);
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
          "[snarl_manager_decomposition][snarl_distance_interface]") {
    require_laws_on_implementation("SnarlDistanceIndex");
}

TEST_CASE("SnarlManagerDecomposition get_snarl_net() gives each manager snarl oriented forward along its chain",
          "[snarl_manager_decomposition]") {
    // The root chain is node 1, the snarl from 1 to 8, and node 8. The snarl
    // from 1 to 8 holds the chain of node 2, the snarl from 2 to 7, and node 7.
    // The snarl from 2 to 7 holds node 6 and the chain of node 3, the snarl
    // from 3 to 5, and node 5.
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
    // Reverse some chains, so that some snarls are reversed in parent.
    IntegratedSnarlFinder base_finder(graph);
    SnarlDecompositionFuzzer finder(&graph, &base_finder, 0.5, generator);
    SnarlManager manager = finder.find_snarls();
    SnarlManagerDecomposition view(manager, graph);

    REQUIRE(child_names(view, graph, view.get_root()) == vector<string> {"chain 1+/8+"});
    view.for_each_child(view.get_root(), [&](const net_handle_t& root_chain) {
        REQUIRE(child_names(view, graph, root_chain).size() == 3);
    });

    manager.for_each_snarl_preorder([&](const Snarl* snarl) {
        net_handle_t net = view.get_snarl_net(snarl);
        REQUIRE(view.get_snarl(net) == snarl);
        bool is_reverse = (*manager.chain_of(snarl))[manager.chain_rank_of(snarl)].second;
        handle_t start = graph.get_handle(snarl->start().node_id(), snarl->start().backward());
        handle_t end = graph.get_handle(snarl->end().node_id(), snarl->end().backward());
        REQUIRE(view.get_handle(view.get_bound(net, false, true), &graph) == (is_reverse ? graph.flip(end) : start));
        REQUIRE(view.get_handle(view.get_bound(net, true, false), &graph) == (is_reverse ? graph.flip(start) : end));
        bool listed = false;
        view.for_each_child(view.get_parent(net), [&](const net_handle_t& sibling) {
            listed = listed || sibling == net;
        });
        // for_each_child() of the parent chain of the snarl lists the same
        // snarl traversal.
        REQUIRE(listed);
    });
}

TEST_CASE("SnarlManagerDecomposition for_each_child() of a chain lists a childless snarl between the snarl's boundary nodes",
          "[snarl_manager_decomposition]") {
    // The chain is nodes 1, 2, 4 and 5, with a snarl between each two. The
    // snarl from 2 to 4 holds node 3, and the other two are trivial snarls.
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

    size_t root_chain_count = 0;
    view.for_each_child(view.get_root(), [&](const net_handle_t& root_chain) {
        root_chain_count++;
        REQUIRE(child_names(view, graph, root_chain) ==
                vector<string> {"node 1", "snarl 1+/2+", "node 2", "snarl 2+/4+", "node 4", "snarl 4+/5+", "node 5"});
    });
    REQUIRE(root_chain_count == 1);
}

TEST_CASE("SnarlManagerDecomposition for_each_traversal() gives START_START for a chain with a self-loop in a childless snarl",
          "[snarl_manager_decomposition]") {
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
    // end enters node 4 through its end and leaves through its start. The
    // self-loop can only be taken leaving 4 through its end, so this walk
    // cannot turn around.
    REQUIRE(traversal_kinds(view, chain) == vector<string> {"SS", "SE", "ES"});
}

TEST_CASE("SnarlManagerDecomposition for_each_child() of a looping chain lists the boundary node once and follow_net_edges() of the last child gives the boundary node",
          "[snarl_manager_decomposition]") {
    // The chain is node 1, the snarl from 1 to 4, node 4, and the snarl from
    // 4 back to 1, which has no children.
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

    vector<net_handle_t> chains;
    view.for_each_child(view.get_root(), [&chains](const net_handle_t& chain) {
        chains.push_back(chain);
    });
    REQUIRE(chains.size() == 1);
    const net_handle_t& chain = chains[0];
    vector<net_handle_t> children;
    view.for_each_child(chain, [&children](const net_handle_t& child) {
        children.push_back(child);
    });
    REQUIRE(children.size() == 4);
    // Both bounds are the first child, oriented forward along the chain.
    REQUIRE(view.get_bound(chain, false, true) == children.front());
    REQUIRE(view.get_bound(chain, true, false) == children.front());

    vector<net_handle_t> after_last;
    view.follow_net_edges(children.back(), &graph, false, [&after_last](const net_handle_t& next) {
        after_last.push_back(next);
    });
    // Leaving the last child forward along the chain, a walk next visits the
    // first child.
    REQUIRE(after_last == vector<net_handle_t> {children.front()});

    vector<net_handle_t> flipped_children;
    view.for_each_child(view.flip(chain), [&flipped_children](const net_handle_t& child) {
        flipped_children.push_back(child);
    });
    // Traversed in reverse, the chain still lists its boundary node first.
    REQUIRE(flipped_children ==
            vector<net_handle_t> {view.flip(children[0]), view.flip(children[3]), view.flip(children[2]),
                                  view.flip(children[1])});

    // Each edge on the boundary node is in one of the snarls of the chain, so
    // follow_net_edges() of the chain gives nothing in either direction.
    for (const net_handle_t& traversal : {chain, view.flip(chain)}) {
        for (bool go_left : {false, true}) {
            bool gave_any = false;
            view.follow_net_edges(traversal, &graph, go_left, [&gave_any](const net_handle_t&) {
                gave_any = true;
            });
            REQUIRE(!gave_any);
        }
    }
}

TEST_CASE("SnarlManagerDecomposition puts a snarl whose start() and end() Visits are equal in a looping chain of its own",
          "[snarl_manager_decomposition]") {
    // Nodes 2 and 3 each make a loop from the end of node 1 back to its start.
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

TEST_CASE("SnarlManagerDecomposition constructor throws invalid_argument for a unary snarl",
          "[snarl_manager_decomposition]") {
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

TEST_CASE("SnarlManagerDecomposition get_parent_traversal() throws invalid_argument for nodes of two trivial chains or an unrealizable traversal",
          "[snarl_manager_decomposition]") {
    // The chain is node 1, the snarl from 1 to 4, and node 4. The snarl holds
    // nodes 2 and 3.
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

TEST_CASE("SnarlManagerDecomposition get_parent_traversal() throws invalid_argument for sentinels facing the wrong way",
          "[snarl_manager_decomposition]") {
    // The chain is node 1, the snarl from 1 to 4, and node 4. The snarl holds
    // nodes 2 and 3.
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
    net_handle_t snarl = view.get_parent(view.get_parent(view.get_net(graph.get_handle(n2->id()), &graph)));
    REQUIRE(view.is_snarl(snarl));
    net_handle_t start_in = view.get_bound(snarl, false, true);
    net_handle_t end_out = view.get_bound(snarl, true, false);
    REQUIRE(view.get_parent_traversal(start_in, end_out) == view.get_snarl_net(view.get_snarl(snarl)));
    REQUIRE_THROWS_AS(view.get_parent_traversal(view.flip(start_in), end_out), std::invalid_argument);
    REQUIRE_THROWS_AS(view.get_parent_traversal(start_in, view.flip(end_out)), std::invalid_argument);
}

TEST_CASE("SnarlManagerDecomposition for_each_traversal() gives START_TIP and TIP_START for a snarl with a walk from its start to a tip",
          "[snarl_manager_decomposition]") {
    // The chain is node 1, the snarl from 1 to 4, and node 4. The snarl holds
    // nodes 2, 3 and 5.
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
    // An edge joins the end of node 3 to node 5, and the end of node 5 is a
    // tip.
    graph.create_edge(n3, n5);
    IntegratedSnarlFinder finder(graph);
    SnarlManager manager = finder.find_snarls();
    SnarlManagerDecomposition view(manager, graph);

    vector<net_handle_t> root_chains;
    view.for_each_child(view.get_root(), [&root_chains](const net_handle_t& root_chain) {
        root_chains.push_back(root_chain);
    });
    REQUIRE(root_chains.size() == 1);
    const net_handle_t& chain = root_chains[0];
    vector<net_handle_t> chain_children;
    view.for_each_child(chain, [&chain_children](const net_handle_t& child) {
        chain_children.push_back(child);
    });
    net_handle_t snarl = chain_children[1];
    REQUIRE(view.is_snarl(snarl));
    REQUIRE(traversal_kinds(view, snarl) == vector<string> {"SE", "ST", "ES", "TS"});
    REQUIRE(traversal_kinds(view, chain) == vector<string> {"SE", "ST", "ES", "TS"});

    SECTION("for_each_tippy_child() of the snarl gives node 5's trivial chain END_START") {
        net_handle_t tip_chain = view.get_parent(view.get_net(graph.get_handle(n5->id(), true), &graph));
        vector<net_handle_t> tippy;
        view.for_each_tippy_child(snarl, [&tippy](const net_handle_t& child) {
            tippy.push_back(child);
        });
        REQUIRE(tippy == vector<net_handle_t> {tip_chain});
        REQUIRE(view.starts_at(tip_chain) == SnarlDecomposition::END);
        REQUIRE(view.ends_at(tip_chain) == SnarlDecomposition::START);
        net_handle_t from_tip = view.get_parent_traversal(tip_chain, view.get_bound(snarl, false, false));
        REQUIRE(view.starts_at(from_tip) == SnarlDecomposition::TIP);
        REQUIRE(view.ends_at(from_tip) == SnarlDecomposition::START);
    }

    SECTION("get_parent_traversal() gives START_TIP from the start bound and node 5's trivial chain") {
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

    SECTION("follow_net_edges() from the TIP end of a snarl traversal gives nothing") {
        view.for_each_traversal(snarl, [&](const net_handle_t& traversal) {
            for (bool go_left : {false, true}) {
                if (go_left ? view.starts_at_tip(traversal) : view.ends_at_tip(traversal)) {
                    bool gave_any = false;
                    view.follow_net_edges(traversal, &graph, go_left, [&gave_any](const net_handle_t&) {
                        gave_any = true;
                    });
                    REQUIRE(!gave_any);
                }
            }
        });
    }
}

TEST_CASE("SnarlManagerDecomposition for_each_traversal() gives START_TIP for a snarl with a walk from its start to a tip in a child chain",
          "[snarl_manager_decomposition]") {
    // The snarl from 1 to 7 holds the chain of node 2, the snarl from 2 to 5,
    // and node 5. In the snarl from 2 to 5, an edge joins the end of 3 to node
    // 6, and the end of 6 is a tip. The test builds the snarls by hand,
    // because a snarl finder may make node 6 an end of the root chain.
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
