/** \file
 *
 * Unit tests for various graph algorithms.
 */

// These are the algorithms we want to test.
#include "../algorithms/extract_subchain.hpp"

// Other includes.
#include <bdsg/hash_graph.hpp>

#include <random>

#include "catch.hpp"

namespace vg {

namespace unittest {

//------------------------------------------------------------------------------

namespace {

// Generates a random graph with the given number of snarls, each containing the given number of nodes.
// There is one shared node between each pair of adjacent snarls.
// Also returns the list of ids for the nodes between the snarls.
// If a snarl is bounded by nodes u and w, then every node v with u < v < w is in the snarl.
std::pair<bdsg::HashGraph, std::vector<nid_t>> random_graph(size_t snarls, size_t snarl_size, std::uint32_t seed = 0x12345678) {
    std::vector<std::string> labels { "A", "C", "G", "T" };
    bdsg::HashGraph graph;
    std::vector<nid_t> shared_nodes;
    auto create_node = [&](nid_t id, bool shared) {
        graph.create_handle(labels[id % labels.size()], id);
        if (shared) {
            shared_nodes.push_back(id);
        }
    };

    create_node(1, true);
    std::mt19937 rng(seed);
    nid_t next = 2;
    for (size_t i = 0; i < snarls; i++) {
        for (nid_t id = next; id < next + snarl_size; id++) {
            create_node(id, false);
        }
        create_node(next + snarl_size, true);

        // Chain of even offsets.
        nid_t prev = next - 1;
        for (nid_t id = next; id < next + snarl_size; id += 2) {
            graph.create_edge(graph.get_handle(prev), graph.get_handle(id));
            prev = id;
        }
        graph.create_edge(graph.get_handle(prev), graph.get_handle(next + snarl_size));

        // Chain of odd offsets.
        prev = next - 1;
        for (nid_t id = next + 1; id < next + snarl_size; id += 2) {
            graph.create_edge(graph.get_handle(prev), graph.get_handle(id));
            prev = id;
        }
        graph.create_edge(graph.get_handle(prev), graph.get_handle(next + snarl_size));

        // Every other edge has a 5% chance of existing in random orientation.
        for (nid_t from = next - 1; from <= next + snarl_size; from++) {
            for (nid_t to = from + 1; to <= next + snarl_size; to++) {
                if (rng() % 20 == 0) {
                    bool from_is_reverse = (from == next - 1 ? false : rng() % 2 == 0);
                    bool to_is_reverse = (to == next + snarl_size ? false : rng() % 2 == 0);
                    graph.create_edge(graph.get_handle(from, from_is_reverse), graph.get_handle(to, to_is_reverse));
                }
            }
        }

        next += snarl_size + 1;
    }

    return std::make_pair(graph, shared_nodes);
}

} // anonymous namespace

//------------------------------------------------------------------------------

TEST_CASE("Extract subchain", "[graph_algorithms]") {
    bdsg::HashGraph graph;
    std::vector<nid_t> shared_nodes;
    std::tie(graph, shared_nodes) = random_graph(10, 10);

    auto check_subchain = [&](const hash_set<nid_t>& nodes, nid_t min, nid_t max) {
        size_t expected_size = max - min + 1;
        REQUIRE(nodes.size() == expected_size);
        for (nid_t id : nodes) {
            REQUIRE(min <= id);
            REQUIRE(id <= max);
        }
    };

    SECTION("Empty subchain") {
        nid_t start = 1, end = 1;
        hash_set<nid_t> subchain = extract_subchain(graph, graph.get_handle(start), graph.get_handle(end));
        check_subchain(subchain, start, end);
    }

    SECTION("Empty subchain, reverse") {
        nid_t start = 1, end = 1;
        handle_t from = graph.get_handle(end, true);
        handle_t to = graph.get_handle(start, true);
        hash_set<nid_t> subchain = extract_subchain(graph, from, to);
        check_subchain(subchain, start, end);
    }

    SECTION("Single snarl") {
        for (size_t i = 0; i + 1 < shared_nodes.size(); i++) {
            nid_t start = shared_nodes[i], end = shared_nodes[i + 1];
            hash_set<nid_t> subchain = extract_subchain(graph, graph.get_handle(start), graph.get_handle(end));
            check_subchain(subchain, start, end);
        }
    }

    SECTION("Single snarl, reverse") {
        for (size_t i = 0; i + 1 < shared_nodes.size(); i++) {
            nid_t start = shared_nodes[i], end = shared_nodes[i + 1];
            handle_t from = graph.get_handle(end, true);
            handle_t to = graph.get_handle(start, true);
            hash_set<nid_t> subchain = extract_subchain(graph, from, to);
            check_subchain(subchain, start, end);
        }
    }

    SECTION("Multiple snarls") {
        for (size_t i = 0; i < shared_nodes.size(); i++) {
            for (size_t j = i + 1; j < shared_nodes.size(); j++) {
                nid_t start = shared_nodes[i], end = shared_nodes[j];
                hash_set<nid_t> subchain = extract_subchain(graph, graph.get_handle(start), graph.get_handle(end));
                check_subchain(subchain, start, end);
            }
        }
    }

    SECTION("Multiple snarls, reverse") {
        for (size_t i = 0; i < shared_nodes.size(); i++) {
            for (size_t j = i + 1; j < shared_nodes.size(); j++) {
                nid_t start = shared_nodes[i], end = shared_nodes[j];
                handle_t from = graph.get_handle(end, true);
                handle_t to = graph.get_handle(start, true);
                hash_set<nid_t> subchain = extract_subchain(graph, from, to);
                check_subchain(subchain, start, end);
            }
        }
    }
}

//------------------------------------------------------------------------------

} // namespace unittest

} // namespace vg
