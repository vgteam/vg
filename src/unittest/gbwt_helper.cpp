/** \file
 *
 * Unit tests for gbwt_helper.cpp, which provides utility functions using GBWT.
 */

#include "../gbwt_helper.hpp"
#include "../json2pb.h"

#include "catch.hpp"

#include <set>
#include <vector>

#include <omp.h>

namespace vg {
namespace unittest {

namespace {

/*
  A toy graph for GA(T|GGG)TA(C|A)A with some additional edges.
*/
const std::string gbwt_helper_graph = R"(
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

} // anonymous namespace

// Build a GBWT with three threads including a duplicate.
gbwt::GBWT build_gbwt_index() {
    std::vector<gbwt::vector_type> gbwt_threads {
        short_path, alt_path, short_path
    };

    return get_gbwt(gbwt_threads);
}

TEST_CASE("GBWTGraph works correctly", "[gbwt_helper]") {

    // Build an XG index.
    Graph graph;
    json2pb(graph, gbwt_helper_graph.c_str(), gbwt_helper_graph.size());
    xg::XG xg_index(graph);

    // Build a GBWT with three threads including a duplicate.
    gbwt::GBWT gbwt_index = build_gbwt_index();

   // Build a GBWT-backed graph.
    GBWTGraph gbwt_graph(gbwt_index, xg_index);

    // Verify node ids.
    std::set<id_t> correct_nodes {
        static_cast<id_t>(1),
        static_cast<id_t>(2),
        static_cast<id_t>(4),
        static_cast<id_t>(5),
        static_cast<id_t>(6),
        static_cast<id_t>(7),
        static_cast<id_t>(8),
        static_cast<id_t>(9)
    };
    SECTION("node id space is correct") {
        REQUIRE(gbwt_graph.node_size() == correct_nodes.size());
        REQUIRE(gbwt_graph.min_node_id() == *(correct_nodes.begin()));
        REQUIRE(gbwt_graph.max_node_id() == *(correct_nodes.rbegin()));
        for (id_t id = gbwt_graph.min_node_id(); id <= gbwt_graph.max_node_id(); id++) {
            bool should_exist = (correct_nodes.find(id) != correct_nodes.end());
            REQUIRE(gbwt_graph.has_node(id) == should_exist);
        }
    }

    // Conversions between handles and (id, orientation) pairs.
    SECTION("handles encode node ids and orientations") {
        for (id_t id : correct_nodes) {
            handle_t forward_handle = gbwt_graph.get_handle(id, false);
            handle_t reverse_handle = gbwt_graph.get_handle(id, true);
            REQUIRE(gbwt_graph.get_id(forward_handle) == id);
            REQUIRE(!(gbwt_graph.get_is_reverse(forward_handle)));
            REQUIRE(gbwt_graph.get_id(reverse_handle) == id);
            REQUIRE(gbwt_graph.get_is_reverse(reverse_handle));
            handle_t flipped_fw = gbwt_graph.flip(forward_handle);
            handle_t flipped_rev = gbwt_graph.flip(reverse_handle);
            REQUIRE(as_integer(forward_handle) != as_integer(reverse_handle));
            REQUIRE(as_integer(flipped_fw) == as_integer(reverse_handle));
            REQUIRE(as_integer(flipped_rev) == as_integer(forward_handle));
        }
    }

    // Node sequences and sequence lengths.
    SECTION("node sequences match those in XG") {
        for (id_t id : correct_nodes) {
            handle_t gbwt_fw = gbwt_graph.get_handle(id, false);
            handle_t gbwt_rev = gbwt_graph.get_handle(id, true);
            handle_t xg_fw = xg_index.get_handle(id, false);
            handle_t xg_rev = xg_index.get_handle(id, true);
            REQUIRE(gbwt_graph.get_length(gbwt_fw) == xg_index.get_length(xg_fw));
            REQUIRE(gbwt_graph.get_sequence(gbwt_fw) == xg_index.get_sequence(xg_fw));
            REQUIRE(gbwt_graph.get_length(gbwt_rev) == xg_index.get_length(xg_rev));
            REQUIRE(gbwt_graph.get_sequence(gbwt_rev) == xg_index.get_sequence(xg_rev));
        }
    }

    // Verify edges.
    typedef std::pair<gbwt::node_type, gbwt::node_type> gbwt_edge;
    std::set<gbwt_edge> correct_edges {
        { gbwt::Node::encode(1, false), gbwt::Node::encode(2, false) },
        { gbwt::Node::encode(1, false), gbwt::Node::encode(4, false) },
        { gbwt::Node::encode(2, false), gbwt::Node::encode(4, false) },
        { gbwt::Node::encode(4, false), gbwt::Node::encode(5, false) },
        { gbwt::Node::encode(5, false), gbwt::Node::encode(6, false) },
        { gbwt::Node::encode(6, false), gbwt::Node::encode(7, false) },
        { gbwt::Node::encode(6, false), gbwt::Node::encode(8, false) },
        { gbwt::Node::encode(7, false), gbwt::Node::encode(9, false) },
        { gbwt::Node::encode(8, false), gbwt::Node::encode(9, false) }
    };
    std::set<gbwt_edge> reverse_edges;
    for (gbwt_edge edge : correct_edges) {
        reverse_edges.insert(gbwt_edge(gbwt::Node::reverse(edge.second), gbwt::Node::reverse(edge.first)));
    }
    SECTION("graph contains correct edges") {
        std::set<gbwt_edge> fw_succ, fw_pred, rev_succ, rev_pred;
        for (id_t id : correct_nodes) {
            handle_t forward_handle = gbwt_graph.get_handle(id, false);
            handle_t reverse_handle = gbwt_graph.get_handle(id, true);
            gbwt_graph.follow_edges(forward_handle, false, [&](const handle_t& handle) {
                fw_succ.insert(gbwt_edge(GBWTGraph::handle_to_node(forward_handle), GBWTGraph::handle_to_node(handle)));
            });
            gbwt_graph.follow_edges(forward_handle, true, [&](const handle_t& handle) {
                fw_pred.insert(gbwt_edge(GBWTGraph::handle_to_node(handle), GBWTGraph::handle_to_node(forward_handle)));
            });
            gbwt_graph.follow_edges(reverse_handle, false, [&](const handle_t& handle) {
                rev_succ.insert(gbwt_edge(GBWTGraph::handle_to_node(reverse_handle), GBWTGraph::handle_to_node(handle)));
            });
            gbwt_graph.follow_edges(reverse_handle, true, [&](const handle_t& handle) {
                rev_pred.insert(gbwt_edge(GBWTGraph::handle_to_node(handle), GBWTGraph::handle_to_node(reverse_handle)));
            });
        }
        REQUIRE(fw_succ == correct_edges);
        REQUIRE(fw_pred == correct_edges);
        REQUIRE(rev_succ == reverse_edges);
        REQUIRE(rev_pred == reverse_edges);
    }

    // for_each_handle() finds all nodes.
    SECTION("iterate over all forward handles") {
        std::vector<handle_t> found_handles, parallel_handles;
        gbwt_graph.for_each_handle([&](const handle_t& handle) {
            found_handles.push_back(handle);
        }, false);
        int old_thread_count = omp_get_max_threads();
        omp_set_num_threads(2);
        gbwt_graph.for_each_handle([&](const handle_t& handle) {
#pragma omp critical
            {
                parallel_handles.push_back(handle);
            }
        }, false);
        omp_set_num_threads(old_thread_count);
        REQUIRE(found_handles.size() == correct_nodes.size());
        REQUIRE(parallel_handles.size() == correct_nodes.size());
        for (handle_t& handle : found_handles) {
            REQUIRE(correct_nodes.find(gbwt_graph.get_id(handle)) != correct_nodes.end());
            REQUIRE(!(gbwt_graph.get_is_reverse(handle)));
        }
        for (handle_t& handle : parallel_handles) {
            REQUIRE(correct_nodes.find(gbwt_graph.get_id(handle)) != correct_nodes.end());
            REQUIRE(!(gbwt_graph.get_is_reverse(handle)));
        }
    }

    // Extract GBWT paths using forward traversal.
    std::set<gbwt::vector_type> correct_paths { alt_path, short_path };
    SECTION("haplotype-consistent forward traversal") {
        typedef std::pair<gbwt::SearchState, gbwt::vector_type> state_type;
        std::vector<gbwt::vector_type> found_paths;
        std::stack<state_type> states;
        handle_t first_node = gbwt_graph.get_handle(1, false);
        states.push({ gbwt_graph.get_state(first_node),
                      { static_cast<gbwt::vector_type::value_type>(GBWTGraph::handle_to_node(first_node)) } });
        while (!states.empty()) {
            state_type curr = states.top();
            states.pop();
            bool extend_success = false;
            gbwt_graph.follow_edges(curr.first, [&](const gbwt::SearchState& next_search) -> bool {
                if (!next_search.empty()) {
                    extend_success = true;
                    gbwt::vector_type next_path = curr.second;
                    next_path.push_back(next_search.node);
                    states.push( { next_search, next_path });
                }
                return true;
            });
            if (!extend_success) {
                found_paths.push_back(curr.second);
            }
        }
        REQUIRE(found_paths.size() == correct_paths.size());
        for (const gbwt::vector_type& path : found_paths) {
            REQUIRE(correct_paths.find(path) != correct_paths.end());
        }
    }

    // Extract GBWT paths using bidirectional traversal.
    SECTION("haplotype-consistent bidirectional traversal") {
        typedef std::pair<gbwt::BidirectionalState, gbwt::vector_type> state_type;
        std::vector<gbwt::vector_type> found_paths;
        std::stack<state_type> fw_states, rev_states;
        handle_t first_node = gbwt_graph.get_handle(4, false);
        fw_states.push({ gbwt_graph.get_bd_state(first_node),
                       { static_cast<gbwt::vector_type::value_type>(GBWTGraph::handle_to_node(first_node)) } });
        while (!fw_states.empty()) {
            state_type curr = fw_states.top();
            fw_states.pop();
            bool extend_success = false;
            gbwt_graph.follow_edges(curr.first, false, [&](const gbwt::BidirectionalState& next_search) -> bool {
                if (!next_search.empty()) {
                    extend_success = true;
                    gbwt::vector_type next_path = curr.second;
                    next_path.push_back(next_search.forward.node);
                    fw_states.push( { next_search, next_path });
                }
                return true;
            });
            if (!extend_success) {
                rev_states.push(curr);
            }
        }
        while (!rev_states.empty()) {
            state_type curr = rev_states.top();
            rev_states.pop();
            bool extend_success = false;
            gbwt_graph.follow_edges(curr.first, true, [&](const gbwt::BidirectionalState& next_search) -> bool {
                if (!next_search.empty()) {
                    extend_success = true;
                    gbwt::vector_type next_path {
                        static_cast<gbwt::vector_type::value_type>(gbwt::Node::reverse(next_search.backward.node))
                    };
                    next_path.insert(next_path.end(), curr.second.begin(), curr.second.end());
                    rev_states.push( { next_search, next_path });
                }
                return true;
            });
            if (!extend_success) {
                found_paths.push_back(curr.second);
            }
        }
        REQUIRE(found_paths.size() == correct_paths.size());
        for (const gbwt::vector_type& path : found_paths) {
            REQUIRE(correct_paths.find(path) != correct_paths.end());
        }
    }
}

TEST_CASE("for_each_kmer() finds the correct kmers", "[gbwt_helper]") {

    // Build an XG index.
    Graph graph;
    json2pb(graph, gbwt_helper_graph.c_str(), gbwt_helper_graph.size());
    xg::XG xg_index(graph);

    // Build a GBWT with three threads including a duplicate.
    gbwt::GBWT gbwt_index = build_gbwt_index();

    // These are the 3-mers the traversal should find.
    typedef std::pair<std::vector<std::pair<pos_t, size_t>>, std::string> kmer_type;
    std::set<kmer_type> correct_kmers {
        // alt_path forward
        { { { make_pos_t(1, false, 0), static_cast<size_t>(1) },
            { make_pos_t(2, false, 0), static_cast<size_t>(1) },
            { make_pos_t(4, false, 0), static_cast<size_t>(1) } }, "GAG" },
        { { { make_pos_t(2, false, 0), static_cast<size_t>(1) },
            { make_pos_t(4, false, 0), static_cast<size_t>(2) } }, "AGG" },
        { { { make_pos_t(4, false, 0), static_cast<size_t>(3) } }, "GGG" },
        { { { make_pos_t(4, false, 1), static_cast<size_t>(2) },
            { make_pos_t(5, false, 0), static_cast<size_t>(1) } }, "GGT" },
        { { { make_pos_t(4, false, 2), static_cast<size_t>(1) },
            { make_pos_t(5, false, 0), static_cast<size_t>(1) },
            { make_pos_t(6, false, 0), static_cast<size_t>(1) } }, "GTA" },
        { { { make_pos_t(5, false, 0), static_cast<size_t>(1) },
            { make_pos_t(6, false, 0), static_cast<size_t>(1) },
            { make_pos_t(8, false, 0), static_cast<size_t>(1) } }, "TAA" },
        { { { make_pos_t(6, false, 0), static_cast<size_t>(1) },
            { make_pos_t(8, false, 0), static_cast<size_t>(1) },
            { make_pos_t(9, false, 0), static_cast<size_t>(1) } }, "AAA" },

        // alt_path reverse
        { { { make_pos_t(9, true, 0), static_cast<size_t>(1) },
            { make_pos_t(8, true, 0), static_cast<size_t>(1) },
            { make_pos_t(6, true, 0), static_cast<size_t>(1) } }, "TTT" },
        { { { make_pos_t(8, true, 0), static_cast<size_t>(1) },
            { make_pos_t(6, true, 0), static_cast<size_t>(1) },
            { make_pos_t(5, true, 0), static_cast<size_t>(1) } }, "TTA" },
        { { { make_pos_t(6, true, 0), static_cast<size_t>(1) },
            { make_pos_t(5, true, 0), static_cast<size_t>(1) },
            { make_pos_t(4, true, 0), static_cast<size_t>(1) } }, "TAC" },
        { { { make_pos_t(5, true, 0), static_cast<size_t>(1) },
            { make_pos_t(4, true, 0), static_cast<size_t>(2) } }, "ACC" },
        { { { make_pos_t(4, true, 0), static_cast<size_t>(3) } }, "CCC" },
        { { { make_pos_t(4, true, 1), static_cast<size_t>(2) },
            { make_pos_t(2, true, 0), static_cast<size_t>(1) } }, "CCT" },
        { { { make_pos_t(4, true, 2), static_cast<size_t>(1) },
            { make_pos_t(2, true, 0), static_cast<size_t>(1) },
            { make_pos_t(1, true, 0), static_cast<size_t>(1) } }, "CTC" },

        // short_path forward
        { { { make_pos_t(1, false, 0), static_cast<size_t>(1) },
            { make_pos_t(4, false, 0), static_cast<size_t>(2) } }, "GGG" },
        { { { make_pos_t(4, false, 0), static_cast<size_t>(3) } }, "GGG" },
        { { { make_pos_t(4, false, 1), static_cast<size_t>(2) },
            { make_pos_t(5, false, 0), static_cast<size_t>(1) } }, "GGT" },
        { { { make_pos_t(4, false, 2), static_cast<size_t>(1) },
            { make_pos_t(5, false, 0), static_cast<size_t>(1) },
            { make_pos_t(6, false, 0), static_cast<size_t>(1) } }, "GTA" },
        { { { make_pos_t(5, false, 0), static_cast<size_t>(1) },
            { make_pos_t(6, false, 0), static_cast<size_t>(1) },
            { make_pos_t(7, false, 0), static_cast<size_t>(1) } }, "TAC" },
        { { { make_pos_t(6, false, 0), static_cast<size_t>(1) },
            { make_pos_t(7, false, 0), static_cast<size_t>(1) },
            { make_pos_t(9, false, 0), static_cast<size_t>(1) } }, "ACA" },

        // short_path reverse
        { { { make_pos_t(9, true, 0), static_cast<size_t>(1) },
            { make_pos_t(7, true, 0), static_cast<size_t>(1) },
            { make_pos_t(6, true, 0), static_cast<size_t>(1) } }, "TGT" },
        { { { make_pos_t(7, true, 0), static_cast<size_t>(1) },
            { make_pos_t(6, true, 0), static_cast<size_t>(1) },
            { make_pos_t(5, true, 0), static_cast<size_t>(1) } }, "GTA" },
        { { { make_pos_t(6, true, 0), static_cast<size_t>(1) },
            { make_pos_t(5, true, 0), static_cast<size_t>(1) },
            { make_pos_t(4, true, 0), static_cast<size_t>(1) } }, "TAC" },
        { { { make_pos_t(5, true, 0), static_cast<size_t>(1) },
            { make_pos_t(4, true, 0), static_cast<size_t>(2) } }, "ACC" },
        { { { make_pos_t(4, true, 0), static_cast<size_t>(3) } }, "CCC" },
        { { { make_pos_t(4, true, 1), static_cast<size_t>(2) },
            { make_pos_t(1, true, 0), static_cast<size_t>(1) } }, "CCC" }
    };

    // Extract all haplotype-consistent 3-mers.
    std::set<kmer_type> found_kmers;
    auto lambda = [&found_kmers](const std::vector<std::pair<pos_t, size_t>>& traversal, const std::string& seq) {
        found_kmers.insert(kmer_type(traversal, seq));
    };
    for_each_kmer(xg_index, gbwt_index, 3, lambda, false);

    // Check the number of kmers.
    SECTION("the traversal should have found the correct number of kmers") {
        REQUIRE(found_kmers.size() == correct_kmers.size());
    }

    // Check that each kmer was found.
    SECTION("the traversal should have found the correct kmers") {
        for (const kmer_type& kmer : correct_kmers) {
            REQUIRE(found_kmers.find(kmer) != found_kmers.end());
        }
    }
}

TEST_CASE("for_each_window() finds the correct windows with GBWT", "[gbwt_helper]") {

    // Build an XG index.
    Graph graph;
    json2pb(graph, gbwt_helper_graph.c_str(), gbwt_helper_graph.size());
    xg::XG xg_index(graph);

    // Build a GBWT with three threads including a duplicate.
    gbwt::GBWT gbwt_index = build_gbwt_index();

    // These are the windows the traversal should find.
    typedef std::pair<std::vector<std::pair<pos_t, size_t>>, std::string> kmer_type;
    std::set<kmer_type> correct_kmers {
        // alt_path forward
        { { { make_pos_t(1, false, 0), static_cast<size_t>(1) },
            { make_pos_t(2, false, 0), static_cast<size_t>(1) },
            { make_pos_t(4, false, 0), static_cast<size_t>(1) } }, "GAG" },
        { { { make_pos_t(2, false, 0), static_cast<size_t>(1) },
            { make_pos_t(4, false, 0), static_cast<size_t>(2) } }, "AGG" },
        { { { make_pos_t(4, false, 0), static_cast<size_t>(3) },
            { make_pos_t(5, false, 0), static_cast<size_t>(1) },
            { make_pos_t(6, false, 0), static_cast<size_t>(1) } }, "GGGTA" },
        { { { make_pos_t(5, false, 0), static_cast<size_t>(1) },
            { make_pos_t(6, false, 0), static_cast<size_t>(1) },
            { make_pos_t(8, false, 0), static_cast<size_t>(1) } }, "TAA" },
        { { { make_pos_t(6, false, 0), static_cast<size_t>(1) },
            { make_pos_t(8, false, 0), static_cast<size_t>(1) },
            { make_pos_t(9, false, 0), static_cast<size_t>(1) } }, "AAA" },

        // alt_path reverse
        { { { make_pos_t(9, true, 0), static_cast<size_t>(1) },
            { make_pos_t(8, true, 0), static_cast<size_t>(1) },
            { make_pos_t(6, true, 0), static_cast<size_t>(1) } }, "TTT" },
        { { { make_pos_t(8, true, 0), static_cast<size_t>(1) },
            { make_pos_t(6, true, 0), static_cast<size_t>(1) },
            { make_pos_t(5, true, 0), static_cast<size_t>(1) } }, "TTA" },
        { { { make_pos_t(6, true, 0), static_cast<size_t>(1) },
            { make_pos_t(5, true, 0), static_cast<size_t>(1) },
            { make_pos_t(4, true, 0), static_cast<size_t>(1) } }, "TAC" },
        { { { make_pos_t(5, true, 0), static_cast<size_t>(1) },
            { make_pos_t(4, true, 0), static_cast<size_t>(2) } }, "ACC" },
        { { { make_pos_t(4, true, 0), static_cast<size_t>(3) },
            { make_pos_t(2, true, 0), static_cast<size_t>(1) },
            { make_pos_t(1, true, 0), static_cast<size_t>(1) } }, "CCCTC" },

        // short_path forward
        { { { make_pos_t(1, false, 0), static_cast<size_t>(1) },
            { make_pos_t(4, false, 0), static_cast<size_t>(2) } }, "GGG" },
        { { { make_pos_t(4, false, 0), static_cast<size_t>(3) },
            { make_pos_t(5, false, 0), static_cast<size_t>(1) },
            { make_pos_t(6, false, 0), static_cast<size_t>(1) } }, "GGGTA" },
        { { { make_pos_t(5, false, 0), static_cast<size_t>(1) },
            { make_pos_t(6, false, 0), static_cast<size_t>(1) },
            { make_pos_t(7, false, 0), static_cast<size_t>(1) } }, "TAC" },
        { { { make_pos_t(6, false, 0), static_cast<size_t>(1) },
            { make_pos_t(7, false, 0), static_cast<size_t>(1) },
            { make_pos_t(9, false, 0), static_cast<size_t>(1) } }, "ACA" },

        // short_path reverse
        { { { make_pos_t(9, true, 0), static_cast<size_t>(1) },
            { make_pos_t(7, true, 0), static_cast<size_t>(1) },
            { make_pos_t(6, true, 0), static_cast<size_t>(1) } }, "TGT" },
        { { { make_pos_t(7, true, 0), static_cast<size_t>(1) },
            { make_pos_t(6, true, 0), static_cast<size_t>(1) },
            { make_pos_t(5, true, 0), static_cast<size_t>(1) } }, "GTA" },
        { { { make_pos_t(6, true, 0), static_cast<size_t>(1) },
            { make_pos_t(5, true, 0), static_cast<size_t>(1) },
            { make_pos_t(4, true, 0), static_cast<size_t>(1) } }, "TAC" },
        { { { make_pos_t(5, true, 0), static_cast<size_t>(1) },
            { make_pos_t(4, true, 0), static_cast<size_t>(2) } }, "ACC" },
        { { { make_pos_t(4, true, 0), static_cast<size_t>(3) },
            { make_pos_t(1, true, 0), static_cast<size_t>(1) } }, "CCCC" }
    };

    // Extract all haplotype-consistent windows of length 3.
    std::set<kmer_type> found_kmers;
    auto lambda = [&found_kmers](const std::vector<std::pair<pos_t, size_t>>& traversal, const std::string& seq) {
        found_kmers.insert(kmer_type(traversal, seq));
    };
    for_each_window(xg_index, gbwt_index, 3, lambda, false);

    // Check the number of kmers.
    SECTION("the traversal should have found the correct number of kmers") {
        REQUIRE(found_kmers.size() == correct_kmers.size());
    }

    // Check that each kmer was found.
    SECTION("the traversal should have found the correct kmers") {
        for (const kmer_type& kmer : correct_kmers) {
            REQUIRE(found_kmers.find(kmer) != found_kmers.end());
        }
    }
}

TEST_CASE("for_each_window() finds the correct windows without GBWT", "[gbwt_helper]") {

    // Build an XG index.
    Graph graph;
    json2pb(graph, gbwt_helper_graph.c_str(), gbwt_helper_graph.size());
    xg::XG xg_index(graph);

    // These are the windows the traversal should find.
    typedef std::pair<std::vector<std::pair<pos_t, size_t>>, std::string> kmer_type;
    std::set<kmer_type> correct_kmers {
        // forward
        { { { make_pos_t(1, false, 0), static_cast<size_t>(1) },
            { make_pos_t(2, false, 0), static_cast<size_t>(1) },
            { make_pos_t(3, false, 0), static_cast<size_t>(1) } }, "GAT" },
        { { { make_pos_t(1, false, 0), static_cast<size_t>(1) },
            { make_pos_t(2, false, 0), static_cast<size_t>(1) },
            { make_pos_t(4, false, 0), static_cast<size_t>(1) } }, "GAG" },
        { { { make_pos_t(1, false, 0), static_cast<size_t>(1) },
            { make_pos_t(4, false, 0), static_cast<size_t>(2) } }, "GGG" },
        { { { make_pos_t(1, false, 0), static_cast<size_t>(1) },
            { make_pos_t(6, false, 0), static_cast<size_t>(1) },
            { make_pos_t(7, false, 0), static_cast<size_t>(1) } }, "GAC" },
        { { { make_pos_t(1, false, 0), static_cast<size_t>(1) },
            { make_pos_t(6, false, 0), static_cast<size_t>(1) },
            { make_pos_t(8, false, 0), static_cast<size_t>(1) } }, "GAA" },
        { { { make_pos_t(2, false, 0), static_cast<size_t>(1) },
            { make_pos_t(3, false, 0), static_cast<size_t>(1) },
            { make_pos_t(5, false, 0), static_cast<size_t>(1) } }, "ATT" },
        { { { make_pos_t(2, false, 0), static_cast<size_t>(1) },
            { make_pos_t(4, false, 0), static_cast<size_t>(2) } }, "AGG" },
        { { { make_pos_t(3, false, 0), static_cast<size_t>(1) },
            { make_pos_t(5, false, 0), static_cast<size_t>(1) },
            { make_pos_t(6, false, 0), static_cast<size_t>(1) } }, "TTA" },
        { { { make_pos_t(4, false, 0), static_cast<size_t>(3) },
            { make_pos_t(5, false, 0), static_cast<size_t>(1) },
            { make_pos_t(6, false, 0), static_cast<size_t>(1) } }, "GGGTA" },
        { { { make_pos_t(5, false, 0), static_cast<size_t>(1) },
            { make_pos_t(6, false, 0), static_cast<size_t>(1) },
            { make_pos_t(7, false, 0), static_cast<size_t>(1) } }, "TAC" },
        { { { make_pos_t(5, false, 0), static_cast<size_t>(1) },
            { make_pos_t(6, false, 0), static_cast<size_t>(1) },
            { make_pos_t(8, false, 0), static_cast<size_t>(1) } }, "TAA" },
        { { { make_pos_t(6, false, 0), static_cast<size_t>(1) },
            { make_pos_t(7, false, 0), static_cast<size_t>(1) },
            { make_pos_t(9, false, 0), static_cast<size_t>(1) } }, "ACA" },
        { { { make_pos_t(6, false, 0), static_cast<size_t>(1) },
            { make_pos_t(8, false, 0), static_cast<size_t>(1) },
            { make_pos_t(9, false, 0), static_cast<size_t>(1) } }, "AAA" },

        // reverse
        { { { make_pos_t(9, true, 0), static_cast<size_t>(1) },
            { make_pos_t(8, true, 0), static_cast<size_t>(1) },
            { make_pos_t(6, true, 0), static_cast<size_t>(1) } }, "TTT" },
        { { { make_pos_t(9, true, 0), static_cast<size_t>(1) },
            { make_pos_t(7, true, 0), static_cast<size_t>(1) },
            { make_pos_t(6, true, 0), static_cast<size_t>(1) } }, "TGT" },
        { { { make_pos_t(8, true, 0), static_cast<size_t>(1) },
            { make_pos_t(6, true, 0), static_cast<size_t>(1) },
            { make_pos_t(5, true, 0), static_cast<size_t>(1) } }, "TTA" },
        { { { make_pos_t(8, true, 0), static_cast<size_t>(1) },
            { make_pos_t(6, true, 0), static_cast<size_t>(1) },
            { make_pos_t(1, true, 0), static_cast<size_t>(1) } }, "TTC" },
        { { { make_pos_t(7, true, 0), static_cast<size_t>(1) },
            { make_pos_t(6, true, 0), static_cast<size_t>(1) },
            { make_pos_t(5, true, 0), static_cast<size_t>(1) } }, "GTA" },
        { { { make_pos_t(7, true, 0), static_cast<size_t>(1) },
            { make_pos_t(6, true, 0), static_cast<size_t>(1) },
            { make_pos_t(1, true, 0), static_cast<size_t>(1) } }, "GTC" },
        { { { make_pos_t(6, true, 0), static_cast<size_t>(1) },
            { make_pos_t(5, true, 0), static_cast<size_t>(1) },
            { make_pos_t(4, true, 0), static_cast<size_t>(1) } }, "TAC" },
        { { { make_pos_t(6, true, 0), static_cast<size_t>(1) },
            { make_pos_t(5, true, 0), static_cast<size_t>(1) },
            { make_pos_t(3, true, 0), static_cast<size_t>(1) } }, "TAA" },
        { { { make_pos_t(5, true, 0), static_cast<size_t>(1) },
            { make_pos_t(4, true, 0), static_cast<size_t>(2) } }, "ACC" },
        { { { make_pos_t(5, true, 0), static_cast<size_t>(1) },
            { make_pos_t(3, true, 0), static_cast<size_t>(1) },
            { make_pos_t(2, true, 0), static_cast<size_t>(1) } }, "AAT" },
        { { { make_pos_t(4, true, 0), static_cast<size_t>(3) },
            { make_pos_t(2, true, 0), static_cast<size_t>(1) },
            { make_pos_t(1, true, 0), static_cast<size_t>(1) } }, "CCCTC" },
        { { { make_pos_t(4, true, 0), static_cast<size_t>(3) },
            { make_pos_t(1, true, 0), static_cast<size_t>(1) } }, "CCCC" },
        { { { make_pos_t(3, true, 0), static_cast<size_t>(1) },
            { make_pos_t(2, true, 0), static_cast<size_t>(1) },
            { make_pos_t(1, true, 0), static_cast<size_t>(1) } }, "ATC" }
    };

    // Extract all haplotype-consistent windows of length 3.
    std::set<kmer_type> found_kmers;
    auto lambda = [&found_kmers](const std::vector<std::pair<pos_t, size_t>>& traversal, const std::string& seq) {
        found_kmers.insert(kmer_type(traversal, seq));
    };
    for_each_window(xg_index, 3, lambda, false);

    // Check the number of kmers.
    SECTION("the traversal should have found the correct number of kmers") {
        REQUIRE(found_kmers.size() == correct_kmers.size());
    }

    // Check that each kmer was found.
    SECTION("the traversal should have found the correct kmers") {
        for (const kmer_type& kmer : correct_kmers) {
            REQUIRE(found_kmers.find(kmer) != found_kmers.end());
        }
    }
}

}
}
