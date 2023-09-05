/** \file
 *
 * Unit tests for the PhaseUnfolder, which replaces pruned regions of the
 * graph with paths supported by XG paths or GBWT threads.
 */

#include <iostream>
#include <map>

#include <omp.h>

#include <gbwt/dynamic_gbwt.h>

#include "../phase_unfolder.hpp"
#include "vg/io/json2pb.h"
#include "xg.hpp"

#include "catch.hpp"

namespace vg {
namespace unittest {

void check_unfolded_nodes(VG& vg_graph,
                          const xg::XG& xg_index,
                          const PhaseUnfolder& unfolder,
                          const std::set<vg::id_t>& expected_nodes,
                          const std::multiset<vg::id_t>& corresponding_nodes) {

    REQUIRE(vg_graph.get_node_count() == expected_nodes.size());

    SECTION("the set of nodes should be the expected one") {
        std::set<vg::id_t> found_nodes;
        vg_graph.for_each_node([&](Node* node) {
            found_nodes.insert(node->id());
        });
        REQUIRE(found_nodes == expected_nodes);
    }

    SECTION("the duplicated nodes should correspond to the right original ones") {
        std::multiset<vg::id_t> found_nodes;
        vg_graph.for_each_node([&](Node* node) {
            found_nodes.insert(unfolder.get_mapping(node->id()));
        });
        REQUIRE(found_nodes == corresponding_nodes);
    }

    SECTION("the duplicated nodes should have the original sequences") {
        bool ok = true;
        vg_graph.for_each_node([&](Node* node) {
            ok &= (node->sequence() == xg_index.get_sequence(xg_index.get_handle(unfolder.get_mapping(node->id()))));
        });
        REQUIRE(ok);
    }
}

void check_unfolded_edges(VG& vg_graph,
                          const PhaseUnfolder& unfolder,
                          std::multiset<std::pair<vg::id_t, vg::id_t>> corresponding_edges,
                          std::set<std::vector<vg::id_t>> expected_paths,
                          vg::id_t next_id) {

    REQUIRE(vg_graph.edge_count() == corresponding_edges.size());

    SECTION("the duplicated edges must correspond to the right original edges") {
        std::multiset<std::pair<vg::id_t, vg::id_t>> found_edges;
        vg_graph.for_each_edge([&](Edge* edge) {
            vg::id_t from_id = unfolder.get_mapping(edge->from());
            vg::id_t to_id = unfolder.get_mapping(edge->to());
            if (edge->from_start() && edge->to_end()) {
                std::swap(from_id, to_id);
            }
            found_edges.emplace(from_id, to_id);
        });
        REQUIRE(found_edges == corresponding_edges);
    }

    SECTION("the duplicated edges must form the right paths") {
        for (vg::id_t node = next_id; node <= vg_graph.max_node_id(); node++) {
            std::vector<Edge*> edges = vg_graph.edges_of(vg_graph.get_node(node));
            std::set<vg::id_t> predecessors, successors;
            for (Edge* edge : edges) {
                vg::id_t from_id = edge->from();
                vg::id_t to_id = edge->to();
                if (edge->from_start() && edge->to_end()) {
                    std::swap(from_id, to_id);
                }
                if (to_id == node) {
                    predecessors.insert(from_id);
                } else {
                    successors.insert(to_id);
                }
            }
            for (vg::id_t pred : predecessors) {
                for (vg::id_t succ : successors) {
                    std::vector<vg::id_t> path { unfolder.get_mapping(pred), unfolder.get_mapping(node), unfolder.get_mapping(succ) };
                    auto iter = expected_paths.find(path);
                    REQUIRE(iter != expected_paths.end());
                    expected_paths.erase(iter);
                }
            }
        }
        REQUIRE(expected_paths.empty());
    }
}

void prune_and_unfold(VG& vg_graph, const std::set<vg::id_t>& to_remove, PhaseUnfolder& unfolder) {
    for (vg::id_t node : to_remove) {
        vg_graph.destroy_node(node);
    }
    unfolder.unfold(vg_graph);

    SECTION("the unfolded graph should contain all indexed paths") {
        omp_set_num_threads(1);
        size_t failures = unfolder.verify_paths(vg_graph);
        REQUIRE(failures == 0);
    }
}

void prune_and_restore(VG& vg_graph, const std::set<vg::id_t>& to_remove, PhaseUnfolder& unfolder) {
    for (vg::id_t node : to_remove) {
        vg_graph.destroy_node(node);
    }
    unfolder.restore_paths(vg_graph);

    SECTION("the restored graph should contain all indexed paths") {
        omp_set_num_threads(1);
        size_t failures = unfolder.verify_paths(vg_graph);
        REQUIRE(failures == 0);
    }
}

/*
  A toy graph for GA(T|GGG)TA(C|A)A with some additional edges. The idea is to remove nodes
  3, 4, 7, 8, and 9. There are two components in the complement graph, and in one of them
  all unfolded paths and threads are maximal. The maximal paths should all end in the same
  duplicate of node 9.
*/
const std::string unfolder_graph = R"(
{
    "node": [
        {"id": 1, "sequence": "G"},
        {"id": 2, "sequence": "A"},
        {"id": 3, "sequence": "T"},
        {"id": 4, "sequence": "GGG"},
        {"id": 5, "sequence": "T"},
        {"id": 6, "sequence": "A"},
        {"id": 7, "sequence": "C"},
        {"id": 8, "sequence": "A"},
        {"id": 9, "sequence": "A"}
    ],
    "edge": [
        {"from": 1, "to": 2},
        {"from": 1, "to": 4},
        {"from": 1, "to": 6},
        {"from": 2, "to": 3},
        {"from": 2, "to": 4},
        {"from": 3, "to": 5},
        {"from": 4, "to": 5},
        {"from": 5, "to": 6},
        {"from": 6, "to": 7},
        {"from": 6, "to": 8},
        {"from": 7, "to": 9},
        {"from": 8, "to": 9}
    ]
}
)";

// The above graph with a reference path.
const std::string unfolder_graph_path = R"(
{
    "node": [
        {"id": 1, "sequence": "G"},
        {"id": 2, "sequence": "A"},
        {"id": 3, "sequence": "T"},
        {"id": 4, "sequence": "GGG"},
        {"id": 5, "sequence": "T"},
        {"id": 6, "sequence": "A"},
        {"id": 7, "sequence": "C"},
        {"id": 8, "sequence": "A"},
        {"id": 9, "sequence": "A"}
    ],
    "edge": [
        {"from": 1, "to": 2},
        {"from": 1, "to": 4},
        {"from": 1, "to": 6},
        {"from": 2, "to": 3},
        {"from": 2, "to": 4},
        {"from": 3, "to": 5},
        {"from": 4, "to": 5},
        {"from": 5, "to": 6},
        {"from": 6, "to": 7},
        {"from": 6, "to": 8},
        {"from": 7, "to": 9},
        {"from": 8, "to": 9}
    ],
    "path": [
        {"name": "hint", "mapping": [
            {"position": {"node_id": 1}, "rank" : 1 },
            {"position": {"node_id": 2}, "rank" : 2 },
            {"position": {"node_id": 3}, "rank" : 3 },
            {"position": {"node_id": 5}, "rank" : 4 },
            {"position": {"node_id": 6}, "rank" : 5 },
            {"position": {"node_id": 7}, "rank" : 6 },
            {"position": {"node_id": 9}, "rank" : 7 }
        ]}
    ]
}
)";

TEST_CASE("PhaseUnfolder can unfold XG paths", "[phaseunfolder][indexing]") {

    // Build an XG index with a path.
    Graph graph_with_path;
    json2pb(graph_with_path, unfolder_graph_path.c_str(), unfolder_graph_path.size());
    xg::XG xg_index;
    xg_index.from_path_handle_graph(VG(graph_with_path));

    // Build an empty GBWT index.
    gbwt::GBWT gbwt_index;

    // Build a PhaseUnfolder with duplicate node ids starting from 10.
    vg::id_t next_id = 10;
    PhaseUnfolder unfolder(xg_index, gbwt_index, next_id);

    // Build a VG graph.
    VG vg_graph;
    Graph temp_graph;
    json2pb(temp_graph, unfolder_graph.c_str(), unfolder_graph.size());
    vg_graph.merge(temp_graph);

    // Remove branching regions from the VG graph, including the last node,
    // but keep the edge (1, 6) in the graph.
    std::set<vg::id_t> to_remove { 3, 4, 7, 8, 9 };
    prune_and_unfold(vg_graph, to_remove, unfolder);

    // Check the nodes.
    SECTION("the unfolded graph should have 7 nodes") {
        std::set<vg::id_t> expected_nodes { 1, 2, 5, 6, 10, 11, 12 };
        std::multiset<vg::id_t> corresponding_nodes { 1, 2, 3, 5, 6, 7, 9 };
        check_unfolded_nodes(vg_graph, xg_index, unfolder, expected_nodes, corresponding_nodes);
    }

    // Check the edges.
    SECTION("the unfolded graph should have 7 edges") {
        std::multiset<std::pair<vg::id_t, vg::id_t>> corresponding_edges {
            { 1, 2 }, { 1, 6 }, { 2, 3 }, { 3, 5 }, { 5, 6 }, { 6, 7 }, { 7, 9 }
        };
        std::set<std::vector<vg::id_t>> expected_paths {
            { 2, 3, 5 }, { 6, 7, 9 }
        };
        check_unfolded_edges(vg_graph, unfolder, corresponding_edges, expected_paths, next_id);
    }
}

TEST_CASE("PhaseUnfolder can restore XG paths", "[phaseunfolder][indexing]") {

    // Build an XG index with a path.
    Graph graph_with_path;
    json2pb(graph_with_path, unfolder_graph_path.c_str(), unfolder_graph_path.size());
    xg::XG xg_index;
    xg_index.from_path_handle_graph(VG(graph_with_path));

    // Build an empty GBWT index.
    gbwt::GBWT gbwt_index;

    // Build a PhaseUnfolder with duplicate node ids starting from 10.
    vg::id_t next_id = 10;
    PhaseUnfolder unfolder(xg_index, gbwt_index, next_id);

    // Build a VG graph.
    VG vg_graph;
    Graph temp_graph;
    json2pb(temp_graph, unfolder_graph.c_str(), unfolder_graph.size());
    vg_graph.merge(temp_graph);

    // Remove branching regions from the VG graph, including the last node,
    // but keep the edge (1, 6) in the graph.
    std::set<vg::id_t> to_remove { 3, 4, 7, 8, 9 };
    prune_and_restore(vg_graph, to_remove, unfolder);

    // Check the nodes.
    SECTION("the unfolded graph should have 7 nodes") {
        std::set<vg::id_t> expected_nodes { 1, 2, 3, 5, 6, 7, 9 };
        std::multiset<vg::id_t> corresponding_nodes { 1, 2, 3, 5, 6, 7, 9 };
        check_unfolded_nodes(vg_graph, xg_index, unfolder, expected_nodes, corresponding_nodes);
    }

    // Check the edges.
    SECTION("the unfolded graph should have 7 edges") {
        std::multiset<std::pair<vg::id_t, vg::id_t>> corresponding_edges {
            { 1, 2 }, { 1, 6 }, { 2, 3 }, { 3, 5 }, { 5, 6 }, { 6, 7 }, { 7, 9 }
        };
        std::set<std::vector<vg::id_t>> expected_paths {
        };
        check_unfolded_edges(vg_graph, unfolder, corresponding_edges, expected_paths, next_id);
    }
}

TEST_CASE("PhaseUnfolder can unfold GBWT threads", "[phaseunfolder][indexing]") {

    // Build an XG index without a path.
    Graph graph_without_path;
    json2pb(graph_without_path, unfolder_graph.c_str(), unfolder_graph.size());
    xg::XG xg_index;
    xg_index.from_path_handle_graph(VG(graph_without_path));

    // Build a GBWT with three threads including a duplicate. We want to have
    // only one instance of short_path unfolded, but we want separate copies
    // of node 4, as it is on both alt_path and short_path.
    gbwt::vector_type alt_path {
        static_cast<gbwt::vector_type::value_type>(gbwt::Node::encode(1, false)),
        static_cast<gbwt::vector_type::value_type>(gbwt::Node::encode(2, false)),
        static_cast<gbwt::vector_type::value_type>(gbwt::Node::encode(4, false)),
        static_cast<gbwt::vector_type::value_type>(gbwt::Node::encode(5, false)),
        static_cast<gbwt::vector_type::value_type>(gbwt::Node::encode(6, false)),
        static_cast<gbwt::vector_type::value_type>(gbwt::Node::encode(8, false)),
        static_cast<gbwt::vector_type::value_type>(gbwt::Node::encode(9, false))
    };
    gbwt::vector_type short_path {
        static_cast<gbwt::vector_type::value_type>(gbwt::Node::encode(1, false)),
        static_cast<gbwt::vector_type::value_type>(gbwt::Node::encode(4, false)),
        static_cast<gbwt::vector_type::value_type>(gbwt::Node::encode(5, false)),
        static_cast<gbwt::vector_type::value_type>(gbwt::Node::encode(6, false)),
        static_cast<gbwt::vector_type::value_type>(gbwt::Node::encode(7, false)),
        static_cast<gbwt::vector_type::value_type>(gbwt::Node::encode(9, false))
    };
    std::vector<gbwt::vector_type> gbwt_threads {
        short_path, alt_path, short_path
    };
    gbwt::GBWT gbwt_index = get_gbwt(gbwt_threads);

    // Build a PhaseUnfolder with duplicate node ids starting from 10.
    vg::id_t next_id = 10;
    PhaseUnfolder unfolder(xg_index, gbwt_index, next_id);

    // Build a VG graph.
    VG vg_graph;
    Graph temp_graph;
    json2pb(temp_graph, unfolder_graph.c_str(), unfolder_graph.size());
    vg_graph.merge(temp_graph);

    // Remove branching regions from the VG graph, including the last node,
    // but keep the edge (1, 6) in the graph.
    std::set<vg::id_t> to_remove { 3, 4, 7, 8, 9 };
    prune_and_unfold(vg_graph, to_remove, unfolder);

    // Check the nodes.
    SECTION("the unfolded graph should have 9 nodes") {
        std::set<vg::id_t> expected_nodes { 1, 2, 5, 6, 10, 11, 12, 13, 14 };
        std::multiset<vg::id_t> corresponding_nodes { 1, 2, 4, 4, 5, 6, 7, 8, 9 };
        check_unfolded_nodes(vg_graph, xg_index, unfolder, expected_nodes, corresponding_nodes);
    }

    // Check the edges.
    SECTION("the unfolded graph should have 11 edges") {
        std::multiset<std::pair<vg::id_t, vg::id_t>> corresponding_edges {
            { 1, 2 }, { 1, 4 }, { 1, 6 }, { 2, 4 }, { 4, 5 }, { 4, 5 }, { 5, 6 }, { 6, 7 }, { 6, 8 }, { 7, 9 }, { 8, 9 }
        };
        std::set<std::vector<vg::id_t>> expected_paths {
            { 1, 4, 5 }, { 2, 4, 5 }, { 6, 7, 9 }, { 6, 8, 9 }
        };
        check_unfolded_edges(vg_graph, unfolder, corresponding_edges, expected_paths, next_id);
    }
}

TEST_CASE("PhaseUnfolder can unfold both XG paths and GBWT threads", "[phaseunfolder][indexing]") {

    // Build an XG index with a path.
    Graph graph_with_path;
    json2pb(graph_with_path, unfolder_graph_path.c_str(), unfolder_graph_path.size());
    xg::XG xg_index;
    xg_index.from_path_handle_graph(VG(graph_with_path));

    // Build a GBWT with three threads including a duplicate. We want to have
    // only one instance of short_path unfolded, but we want separate copies
    // of node 4, as it is on both alt_path and short_path.
    gbwt::vector_type alt_path {
        static_cast<gbwt::vector_type::value_type>(gbwt::Node::encode(1, false)),
        static_cast<gbwt::vector_type::value_type>(gbwt::Node::encode(2, false)),
        static_cast<gbwt::vector_type::value_type>(gbwt::Node::encode(4, false)),
        static_cast<gbwt::vector_type::value_type>(gbwt::Node::encode(5, false)),
        static_cast<gbwt::vector_type::value_type>(gbwt::Node::encode(6, false)),
        static_cast<gbwt::vector_type::value_type>(gbwt::Node::encode(8, false)),
        static_cast<gbwt::vector_type::value_type>(gbwt::Node::encode(9, false))
    };
    gbwt::vector_type short_path {
        static_cast<gbwt::vector_type::value_type>(gbwt::Node::encode(1, false)),
        static_cast<gbwt::vector_type::value_type>(gbwt::Node::encode(4, false)),
        static_cast<gbwt::vector_type::value_type>(gbwt::Node::encode(5, false)),
        static_cast<gbwt::vector_type::value_type>(gbwt::Node::encode(6, false)),
        static_cast<gbwt::vector_type::value_type>(gbwt::Node::encode(7, false)),
        static_cast<gbwt::vector_type::value_type>(gbwt::Node::encode(9, false))
    };
    std::vector<gbwt::vector_type> gbwt_threads {
        short_path, alt_path, short_path
    };
    gbwt::GBWT gbwt_index = get_gbwt(gbwt_threads);

    // Build a PhaseUnfolder with duplicate node ids starting from 10.
    vg::id_t next_id = 10;
    PhaseUnfolder unfolder(xg_index, gbwt_index, next_id);

    // Build a VG graph.
    VG vg_graph;
    Graph temp_graph;
    json2pb(temp_graph, unfolder_graph.c_str(), unfolder_graph.size());
    vg_graph.merge(temp_graph);

    // Remove branching regions from the VG graph, including the last node,
    // but keep the edge (1, 6) in the graph.
    std::set<vg::id_t> to_remove { 3, 4, 7, 8, 9 };
    prune_and_unfold(vg_graph, to_remove, unfolder);

    // Check the nodes.
    SECTION("the unfolded graph should have 10 nodes") {
        std::set<vg::id_t> expected_nodes { 1, 2, 5, 6, 10, 11, 12, 13, 14, 15 };
        std::multiset<vg::id_t> corresponding_nodes { 1, 2, 3, 4, 4, 5, 6, 7, 8, 9 };
        check_unfolded_nodes(vg_graph, xg_index, unfolder, expected_nodes, corresponding_nodes);
    }

    // Check the edges.
    SECTION("the unfolded graph should have 13 edges") {
        std::multiset<std::pair<vg::id_t, vg::id_t>> corresponding_edges {
            { 1, 2 }, { 1, 4 }, { 1, 6 }, { 2, 3 }, { 2, 4 }, { 3, 5 }, { 4, 5 }, { 4, 5 }, { 5, 6 }, { 6, 7 }, { 6, 8 }, { 7, 9 }, { 8, 9 }
        };
        std::set<std::vector<vg::id_t>> expected_paths {
            { 1, 4, 5 }, { 2, 3, 5 }, { 2, 4, 5 }, { 6, 7, 9 }, { 6, 8, 9 }
        };
        check_unfolded_edges(vg_graph, unfolder, corresponding_edges, expected_paths, next_id);
    }
}

/*
  A toy graph for GAT(A|T)ACA. The idea is to remove nodes 3 to 6, which form a
  bubble. All paths through the bubble share prefixes 2 -> 3 and suffixes 6 -> 7,
  which should be merged. We can also try extending a short haplotype with the
  reference.
*/
const std::string unfolder_graph_simple = R"(
{
    "node": [
        {"id": 1, "sequence": "G"},
        {"id": 2, "sequence": "A"},
        {"id": 3, "sequence": "T"},
        {"id": 4, "sequence": "A"},
        {"id": 5, "sequence": "T"},
        {"id": 6, "sequence": "A"},
        {"id": 7, "sequence": "C"},
        {"id": 8, "sequence": "A"}
    ],
    "edge": [
        {"from": 1, "to": 2},
        {"from": 2, "to": 3},
        {"from": 3, "to": 4},
        {"from": 3, "to": 5},
        {"from": 4, "to": 6},
        {"from": 5, "to": 6},
        {"from": 6, "to": 7},
        {"from": 7, "to": 8}
    ]
}
)";

// The above graph with a reference path.
const std::string unfolder_graph_simple_path = R"(
{
    "node": [
        {"id": 1, "sequence": "G"},
        {"id": 2, "sequence": "A"},
        {"id": 3, "sequence": "T"},
        {"id": 4, "sequence": "A"},
        {"id": 5, "sequence": "T"},
        {"id": 6, "sequence": "A"},
        {"id": 7, "sequence": "C"},
        {"id": 8, "sequence": "A"}
    ],
    "edge": [
        {"from": 1, "to": 2},
        {"from": 2, "to": 3},
        {"from": 3, "to": 4},
        {"from": 3, "to": 5},
        {"from": 4, "to": 6},
        {"from": 5, "to": 6},
        {"from": 6, "to": 7},
        {"from": 7, "to": 8}
    ],
    "path": [
        {"name": "hint", "mapping": [
            {"position": {"node_id": 1}, "rank" : 1 },
            {"position": {"node_id": 2}, "rank" : 2 },
            {"position": {"node_id": 3}, "rank" : 3 },
            {"position": {"node_id": 5}, "rank" : 4 },
            {"position": {"node_id": 6}, "rank" : 5 },
            {"position": {"node_id": 7}, "rank" : 6 },
            {"position": {"node_id": 8}, "rank" : 7 }
        ]}
    ]
}
)";

TEST_CASE("PhaseUnfolder can merge shared prefixes and suffixes", "[phaseunfolder][indexing]") {

    // Build an XG index.
    Graph simple_graph;
    json2pb(simple_graph, unfolder_graph_simple.c_str(), unfolder_graph_simple.size());
    xg::XG xg_index;
    xg_index.from_path_handle_graph(VG(simple_graph));

    // Build a GBWT with both possible threads.
    gbwt::vector_type upper_path {
        static_cast<gbwt::vector_type::value_type>(gbwt::Node::encode(1, false)),
        static_cast<gbwt::vector_type::value_type>(gbwt::Node::encode(2, false)),
        static_cast<gbwt::vector_type::value_type>(gbwt::Node::encode(3, false)),
        static_cast<gbwt::vector_type::value_type>(gbwt::Node::encode(4, false)),
        static_cast<gbwt::vector_type::value_type>(gbwt::Node::encode(6, false)),
        static_cast<gbwt::vector_type::value_type>(gbwt::Node::encode(7, false)),
        static_cast<gbwt::vector_type::value_type>(gbwt::Node::encode(8, false))
    };
    gbwt::vector_type lower_path {
        static_cast<gbwt::vector_type::value_type>(gbwt::Node::encode(1, false)),
        static_cast<gbwt::vector_type::value_type>(gbwt::Node::encode(2, false)),
        static_cast<gbwt::vector_type::value_type>(gbwt::Node::encode(3, false)),
        static_cast<gbwt::vector_type::value_type>(gbwt::Node::encode(5, false)),
        static_cast<gbwt::vector_type::value_type>(gbwt::Node::encode(6, false)),
        static_cast<gbwt::vector_type::value_type>(gbwt::Node::encode(7, false)),
        static_cast<gbwt::vector_type::value_type>(gbwt::Node::encode(8, false))
    };
    std::vector<gbwt::vector_type> gbwt_threads {
        upper_path, lower_path
    };
    gbwt::GBWT gbwt_index = get_gbwt(gbwt_threads);

    // Build a PhaseUnfolder with duplicate node ids starting from 9.
    vg::id_t next_id = 9;
    PhaseUnfolder unfolder(xg_index, gbwt_index, next_id);

    // Build a VG graph.
    VG vg_graph;
    Graph temp_graph;
    json2pb(temp_graph, unfolder_graph_simple.c_str(), unfolder_graph_simple.size());
    vg_graph.merge(temp_graph);

    // Remove the bubble, including its endpoints.
    std::set<vg::id_t> to_remove { 3, 4, 5, 6 };
    prune_and_unfold(vg_graph, to_remove, unfolder);

    // Check the nodes.
    SECTION("the unfolded graph should have 8 nodes") {
        std::set<vg::id_t> expected_nodes { 1, 2, 7, 8, 9, 10, 11, 12 };
        std::multiset<vg::id_t> corresponding_nodes { 1, 2, 3, 4, 5, 6, 7, 8 };
        check_unfolded_nodes(vg_graph, xg_index, unfolder, expected_nodes, corresponding_nodes);
    }

    // Check the edges.
    SECTION("the unfolded graph should have 8 edges") {
        std::multiset<std::pair<vg::id_t, vg::id_t>> corresponding_edges {
            { 1, 2 }, { 2, 3 }, { 3, 4 }, { 3, 5 }, { 4, 6 }, { 5, 6 }, { 6, 7 }, { 7, 8 }
        };
        std::set<std::vector<vg::id_t>> expected_paths {
            { 2, 3, 4 }, { 2, 3, 5 }, { 3, 4, 6 }, { 3, 5, 6 }, { 4, 6, 7 }, { 5, 6, 7 }
        };
        check_unfolded_edges(vg_graph, unfolder, corresponding_edges, expected_paths, next_id);
    }
}

TEST_CASE("PhaseUnfolder can extend short threads", "[phaseunfolder][indexing]") {

    // Build an XG index.
    Graph simple_graph_with_path;
    json2pb(simple_graph_with_path, unfolder_graph_simple_path.c_str(), unfolder_graph_simple_path.size());
    xg::XG xg_index;
    xg_index.from_path_handle_graph(VG(simple_graph_with_path));

    // Build a GBWT for the fragment that is different from the reference.
    gbwt::vector_type short_fragment {
        static_cast<gbwt::vector_type::value_type>(gbwt::Node::encode(3, false)),
        static_cast<gbwt::vector_type::value_type>(gbwt::Node::encode(4, false)),
        static_cast<gbwt::vector_type::value_type>(gbwt::Node::encode(6, false))
    };
    std::vector<gbwt::vector_type> gbwt_threads { short_fragment };
    gbwt::GBWT gbwt_index = get_gbwt(gbwt_threads);

    // Build a PhaseUnfolder with duplicate node ids starting from 9.
    vg::id_t next_id = 9;
    PhaseUnfolder unfolder(xg_index, gbwt_index, next_id);

    // Build a VG graph.
    VG vg_graph;
    Graph temp_graph;
    json2pb(temp_graph, unfolder_graph_simple.c_str(), unfolder_graph_simple.size());
    vg_graph.merge(temp_graph);

    // Remove the bubble, including its endpoints.
    std::set<vg::id_t> to_remove { 3, 4, 5, 6 };
    prune_and_unfold(vg_graph, to_remove, unfolder);

    // Check the nodes.
    SECTION("the unfolded graph should have 8 nodes") {
        std::set<vg::id_t> expected_nodes { 1, 2, 7, 8, 9, 10, 11, 12 };
        std::multiset<vg::id_t> corresponding_nodes { 1, 2, 3, 4, 5, 6, 7, 8 };
        check_unfolded_nodes(vg_graph, xg_index, unfolder, expected_nodes, corresponding_nodes);
    }

    // Check the edges.
    SECTION("the unfolded graph should have 8 edges") {
        std::multiset<std::pair<vg::id_t, vg::id_t>> corresponding_edges {
            { 1, 2 }, { 2, 3 }, { 3, 4 }, { 3, 5 }, { 4, 6 }, { 5, 6 }, { 6, 7 }, { 7, 8 }
        };
        std::set<std::vector<vg::id_t>> expected_paths {
            { 2, 3, 4 }, { 2, 3, 5 }, { 3, 4, 6 }, { 3, 5, 6 }, { 4, 6, 7 }, { 5, 6, 7 }
        };
        check_unfolded_edges(vg_graph, unfolder, corresponding_edges, expected_paths, next_id);
    }
}

}
}
