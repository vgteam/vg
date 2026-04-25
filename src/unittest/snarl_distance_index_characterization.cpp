// Characterization tests for snarl_distance_index.
// These tests lock down the serialized byte layout (via FNV-1a hash) of the
// distance index for six canonical graphs.  Hash constants are captured from
// the first run on unmodified code and must remain identical after the file
// split in PR 1.

#include "../snarl_distance_index.hpp"
#include "../integrated_snarl_finder.hpp"
#include <bdsg/hash_graph.hpp>
#include "catch.hpp"
#include "../path.hpp"

namespace vg {
namespace unittest {

static uint64_t fnv1a(const std::vector<uint8_t>& data) {
    uint64_t h = 14695981039346656037ULL;
    for (uint8_t b : data) { h ^= b; h *= 1099511628211ULL; }
    return h;
}

static std::vector<uint8_t> serialize_index(const SnarlDistanceIndex& idx) {
    std::vector<uint8_t> buf;
    idx.serialize([&](const void* p, size_t n) {
        const uint8_t* bytes = static_cast<const uint8_t*>(p);
        buf.insert(buf.end(), bytes, bytes + n);
    });
    return buf;
}

// Walk the snarl tree rooted at handle and collect all snarls into out.
static void collect_snarls(const SnarlDistanceIndex& idx, const net_handle_t& handle,
                            std::vector<net_handle_t>& out) {
    if (idx.is_snarl(handle)) {
        out.push_back(handle);
    }
    // Nodes and sentinels have no snarl-tree children; recursing into them throws.
    if (idx.is_node(handle) || idx.is_sentinel(handle)) {
        return;
    }
    idx.for_each_child(handle, [&](const net_handle_t& child) -> bool {
        collect_snarls(idx, child, out);
        return true;
    });
}

// ─── Fixture 1: linear chain ─────────────────────────────────────────────────
// h1 → h2 → h3 → h4 → h5   (no bubbles, just one chain)

TEST_CASE("Characterization: linear chain", "[snarl_characterization]") {
    bdsg::HashGraph graph;
    handle_t h1 = graph.create_handle("A");
    handle_t h2 = graph.create_handle("A");
    handle_t h3 = graph.create_handle("A");
    handle_t h4 = graph.create_handle("A");
    handle_t h5 = graph.create_handle("A");
    graph.create_edge(h1, h2);
    graph.create_edge(h2, h3);
    graph.create_edge(h3, h4);
    graph.create_edge(h4, h5);

    IntegratedSnarlFinder finder(graph);
    SnarlDistanceIndex idx;
    fill_in_distance_index(&idx, &graph, &finder);

    SECTION("serialization hash") {
        auto buf = serialize_index(idx);
        constexpr size_t EXPECTED_SIZE = 1024;
        constexpr uint64_t EXPECTED_HASH = 4461810471415873827ULL;
        REQUIRE(buf.size() == EXPECTED_SIZE);
        REQUIRE(fnv1a(buf) == EXPECTED_HASH);
    }

    SECTION("no non-trivial snarls") {
        std::vector<net_handle_t> snarls;
        collect_snarls(idx, idx.get_root(), snarls);
        // A linear chain has no non-trivial snarls (root snarls don't count as
        // regular/irregular in the sense of check_regularity).
        for (const net_handle_t& s : snarls) {
            // Root snarls wrapping a component are allowed; skip them.
            if (!idx.is_root_snarl(s)) {
                // No internal snarls expected in a plain linear chain.
                REQUIRE(false);
            }
        }
    }

    SECTION("subgraph in distance range") {
        // From the middle node h3 look forward; nodes 2 steps away should be h5.
        std::unordered_set<nid_t> sub;
        path_handle_t ph = graph.create_path_handle("path_linear");
        graph.append_step(ph, h3);
        Path path = path_from_path_handle(graph, ph);
        subgraph_in_distance_range(idx, path, &graph, 2, 3, sub, true);
        REQUIRE(sub.count(graph.get_id(h5)));
        REQUIRE(!sub.count(graph.get_id(h1)));
        REQUIRE(!sub.count(graph.get_id(h2)));
    }
}

// ─── Fixture 2: simple bubble ─────────────────────────────────────────────────
// h1 → h2 → h4
// h1 → h3 → h4
// One snarl (h1,h4) with two single-node children → should be regular.

TEST_CASE("Characterization: simple bubble", "[snarl_characterization]") {
    bdsg::HashGraph graph;
    handle_t h1 = graph.create_handle("A");
    handle_t h2 = graph.create_handle("A");
    handle_t h3 = graph.create_handle("A");
    handle_t h4 = graph.create_handle("A");
    graph.create_edge(h1, h2);
    graph.create_edge(h1, h3);
    graph.create_edge(h2, h4);
    graph.create_edge(h3, h4);

    IntegratedSnarlFinder finder(graph);
    SnarlDistanceIndex idx;
    fill_in_distance_index(&idx, &graph, &finder);

    SECTION("serialization hash") {
        auto buf = serialize_index(idx);
        constexpr size_t EXPECTED_SIZE = 1024;
        constexpr uint64_t EXPECTED_HASH = 10070957726680237483ULL;
        REQUIRE(buf.size() == EXPECTED_SIZE);
        REQUIRE(fnv1a(buf) == EXPECTED_HASH);
    }

    SECTION("snarl is regular") {
        std::vector<net_handle_t> snarls;
        collect_snarls(idx, idx.get_root(), snarls);
        bool found_internal = false;
        for (const net_handle_t& s : snarls) {
            if (!idx.is_root_snarl(s)) {
                found_internal = true;
                REQUIRE(idx.is_regular_snarl(s));
            }
        }
        REQUIRE(found_internal);
    }

    SECTION("subgraph in distance range") {
        // From h1 look forward; h2 and h3 are 1 step away, h4 is 2 steps away.
        std::unordered_set<nid_t> sub;
        path_handle_t ph = graph.create_path_handle("path_bubble");
        graph.append_step(ph, h1);
        Path path = path_from_path_handle(graph, ph);
        subgraph_in_distance_range(idx, path, &graph, 1, 2, sub, true);
        REQUIRE(sub.count(graph.get_id(h2)));
        REQUIRE(sub.count(graph.get_id(h3)));
        REQUIRE(sub.count(graph.get_id(h4)));
    }
}

// ─── Fixture 3: nested chain with loop ───────────────────────────────────────
// h1 → h2 → h3 → h4 → h5
// h2 → flip(h2)  (self-loop, allows reversing at h2)
// h3 → h5        (shortcut creating a nested snarl)

TEST_CASE("Characterization: nested chain with loop", "[snarl_characterization]") {
    bdsg::HashGraph graph;
    handle_t h1 = graph.create_handle("A");
    handle_t h2 = graph.create_handle("A");
    handle_t h3 = graph.create_handle("A");
    handle_t h4 = graph.create_handle("A");
    handle_t h5 = graph.create_handle("A");
    graph.create_edge(h1, h2);
    graph.create_edge(h2, h3);
    graph.create_edge(h3, h4);
    graph.create_edge(h4, h5);
    graph.create_edge(h2, graph.flip(h2)); // self-loop
    graph.create_edge(h3, h5);            // shortcut

    IntegratedSnarlFinder finder(graph);
    SnarlDistanceIndex idx;
    fill_in_distance_index(&idx, &graph, &finder);

    SECTION("serialization hash") {
        auto buf = serialize_index(idx);
        constexpr size_t EXPECTED_SIZE = 1024;
        constexpr uint64_t EXPECTED_HASH = 16246149163740101819ULL;
        REQUIRE(buf.size() == EXPECTED_SIZE);
        REQUIRE(fnv1a(buf) == EXPECTED_HASH);
    }

    SECTION("snarls exist") {
        std::vector<net_handle_t> snarls;
        collect_snarls(idx, idx.get_root(), snarls);
        // There should be at least one non-root snarl due to the shortcut h3→h5.
        bool found = false;
        for (const net_handle_t& s : snarls) {
            if (!idx.is_root_snarl(s)) { found = true; break; }
        }
        REQUIRE(found);
    }

    SECTION("subgraph in distance range") {
        std::unordered_set<nid_t> sub;
        path_handle_t ph = graph.create_path_handle("path_nested");
        graph.append_step(ph, h1);
        Path path = path_from_path_handle(graph, ph);
        subgraph_in_distance_range(idx, path, &graph, 2, 3, sub, true);
        // h3 is at distance 2 (through h2 then h3) and h5 is reachable via the shortcut.
        REQUIRE(sub.count(graph.get_id(h3)));
    }
}

// ─── Fixture 4: multi-component root ─────────────────────────────────────────
// Component 1: h1 → h2 → h3
// Component 2: h4 → h5 → h6   (no edges between components)

TEST_CASE("Characterization: multi-component root", "[snarl_characterization]") {
    bdsg::HashGraph graph;
    handle_t h1 = graph.create_handle("A");
    handle_t h2 = graph.create_handle("A");
    handle_t h3 = graph.create_handle("A");
    handle_t h4 = graph.create_handle("A");
    handle_t h5 = graph.create_handle("A");
    handle_t h6 = graph.create_handle("A");
    graph.create_edge(h1, h2);
    graph.create_edge(h2, h3);
    graph.create_edge(h4, h5);
    graph.create_edge(h5, h6);

    IntegratedSnarlFinder finder(graph);
    SnarlDistanceIndex idx;
    fill_in_distance_index(&idx, &graph, &finder);

    SECTION("serialization hash") {
        auto buf = serialize_index(idx);
        constexpr size_t EXPECTED_SIZE = 1024;
        constexpr uint64_t EXPECTED_HASH = 13763592152412395439ULL;
        REQUIRE(buf.size() == EXPECTED_SIZE);
        REQUIRE(fnv1a(buf) == EXPECTED_HASH);
    }

    SECTION("connected component count") {
        REQUIRE(idx.connected_component_count() == 2);
    }

    SECTION("subgraph in distance range") {
        // From h1 look forward: h3 is 2 steps away.
        std::unordered_set<nid_t> sub;
        path_handle_t ph = graph.create_path_handle("path_multicomp");
        graph.append_step(ph, h1);
        Path path = path_from_path_handle(graph, ph);
        subgraph_in_distance_range(idx, path, &graph, 2, 3, sub, true);
        REQUIRE(sub.count(graph.get_id(h3)));
        // h4/h5/h6 are in a different component and not reachable.
        REQUIRE(!sub.count(graph.get_id(h4)));
        REQUIRE(!sub.count(graph.get_id(h5)));
        REQUIRE(!sub.count(graph.get_id(h6)));
    }
}

// ─── Fixture 5: oversized snarl ──────────────────────────────────────────────
// h1 → h2 → h6
// h1 → h3 → h6
// h1 → h4 → h6
// h1 → h5 → h6
// Snarl (h1,h6) has 4 internal children; with size_limit=3 → oversized.

TEST_CASE("Characterization: oversized snarl", "[snarl_characterization]") {
    bdsg::HashGraph graph;
    handle_t h1 = graph.create_handle("A");
    handle_t h2 = graph.create_handle("A");
    handle_t h3 = graph.create_handle("A");
    handle_t h4 = graph.create_handle("A");
    handle_t h5 = graph.create_handle("A");
    handle_t h6 = graph.create_handle("A");
    graph.create_edge(h1, h2);
    graph.create_edge(h1, h3);
    graph.create_edge(h1, h4);
    graph.create_edge(h1, h5);
    graph.create_edge(h2, h6);
    graph.create_edge(h3, h6);
    graph.create_edge(h4, h6);
    graph.create_edge(h5, h6);

    IntegratedSnarlFinder finder(graph);
    SnarlDistanceIndex idx;
    fill_in_distance_index(&idx, &graph, &finder, /*size_limit=*/3, false, /*silence_warnings=*/true);

    SECTION("serialization size") {
        // Hub-label content is non-deterministic (contraction hierarchy uses
        // hash-based graph structures), so we only lock down the byte count.
        auto buf = serialize_index(idx);
        constexpr size_t EXPECTED_SIZE = 7168;
        REQUIRE(buf.size() == EXPECTED_SIZE);
    }

    SECTION("oversized snarl exists") {
        std::vector<net_handle_t> snarls;
        collect_snarls(idx, idx.get_root(), snarls);
        bool found_oversized = false;
        for (const net_handle_t& s : snarls) {
            if (!idx.is_root_snarl(s) && idx.is_oversized_snarl(s)) {
                found_oversized = true;
            }
        }
        REQUIRE(found_oversized);
    }

    SECTION("subgraph in distance range") {
        // From h1 look forward; h2,h3,h4,h5 are 1 step away, h6 is 2 steps.
        std::unordered_set<nid_t> sub;
        path_handle_t ph = graph.create_path_handle("path_oversized");
        graph.append_step(ph, h1);
        Path path = path_from_path_handle(graph, ph);
        subgraph_in_distance_range(idx, path, &graph, 1, 2, sub, true);
        REQUIRE(sub.count(graph.get_id(h2)));
        REQUIRE(sub.count(graph.get_id(h3)));
        REQUIRE(sub.count(graph.get_id(h4)));
        REQUIRE(sub.count(graph.get_id(h5)));
        REQUIRE(sub.count(graph.get_id(h6)));
    }
}

// ─── Fixture 6: irregular snarl ──────────────────────────────────────────────
// h1 → h2 → h4
// h1 → h3 → h4
// h2 → h3       (cross-edge between children → snarl is not regular)

TEST_CASE("Characterization: irregular snarl", "[snarl_characterization]") {
    bdsg::HashGraph graph;
    handle_t h1 = graph.create_handle("A");
    handle_t h2 = graph.create_handle("A");
    handle_t h3 = graph.create_handle("A");
    handle_t h4 = graph.create_handle("A");
    graph.create_edge(h1, h2);
    graph.create_edge(h1, h3);
    graph.create_edge(h2, h4);
    graph.create_edge(h3, h4);
    graph.create_edge(h2, h3); // cross-edge

    IntegratedSnarlFinder finder(graph);
    SnarlDistanceIndex idx;
    fill_in_distance_index(&idx, &graph, &finder);

    SECTION("serialization hash") {
        auto buf = serialize_index(idx);
        constexpr size_t EXPECTED_SIZE = 1024;
        constexpr uint64_t EXPECTED_HASH = 14645746962011564342ULL;
        REQUIRE(buf.size() == EXPECTED_SIZE);
        REQUIRE(fnv1a(buf) == EXPECTED_HASH);
    }

    SECTION("snarl is not regular") {
        std::vector<net_handle_t> snarls;
        collect_snarls(idx, idx.get_root(), snarls);
        bool found_irregular = false;
        for (const net_handle_t& s : snarls) {
            if (!idx.is_root_snarl(s) && !idx.is_regular_snarl(s)) {
                found_irregular = true;
            }
        }
        REQUIRE(found_irregular);
    }

    SECTION("subgraph in distance range") {
        // From h1, h2 and h3 are 1 step away, h4 is 2 steps.
        std::unordered_set<nid_t> sub;
        path_handle_t ph = graph.create_path_handle("path_irregular");
        graph.append_step(ph, h1);
        Path path = path_from_path_handle(graph, ph);
        subgraph_in_distance_range(idx, path, &graph, 1, 2, sub, true);
        REQUIRE(sub.count(graph.get_id(h2)));
        REQUIRE(sub.count(graph.get_id(h3)));
        REQUIRE(sub.count(graph.get_id(h4)));
        REQUIRE(!sub.count(graph.get_id(h1)));
    }
}

} // namespace unittest
} // namespace vg
