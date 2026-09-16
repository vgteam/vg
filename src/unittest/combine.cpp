/// \file combine.cpp
///
/// Unit tests for GraphCombiner's internals.
///
/// GraphCombiner::combine() is covered end to end by test/t/51_vg_combine.t.
/// These tests drive the steps underneath it, whose intermediate results -- the ID
/// shifts, the trimmed boundary, the merged coordinates -- the command line cannot see.
///
/// Error branches are absent on purpose: logger.error() exits the process, which would
/// take the whole test binary with it. They live in test/t/51_vg_combine.t.

#include "catch.hpp"
#include "../combine.hpp"

#include "bdsg/hash_graph.hpp"
#include "bdsg/overlays/path_position_overlays.hpp"

namespace vg {
namespace unittest {

using namespace std;
using handlegraph::nid_t;
using handlegraph::offset_t;
using handlegraph::PathMetadata;
using handlegraph::PathSense;
using handlegraph::subrange_t;

// A child class to expose the internals for testing.
class TestGraphCombiner : public GraphCombiner {
public:
    using GraphCombiner::GraphCombiner;
    using GraphCombiner::Chunk;
    using GraphCombiner::SeamState;
    using GraphCombiner::RefKey;
    using GraphCombiner::GroupKey;
    using GraphCombiner::compute_path_length;
    using GraphCombiner::path_start_offset;
    using GraphCombiner::renumber_out_of_way;
    using GraphCombiner::group_key_of;
    using GraphCombiner::path_prefix;
    using GraphCombiner::path_suffix;
    using GraphCombiner::check_no_split_paths;
    using GraphCombiner::describe_reference;
    using GraphCombiner::sort_by_reference_offset;
    using GraphCombiner::check_overlaps_agree;
    using GraphCombiner::check_all_paths_start_at;
    using GraphCombiner::check_all_paths_end_at;
    using GraphCombiner::trim_and_shift;
    using GraphCombiner::ingest;
    using GraphCombiner::share_nodes_by_id;
    using GraphCombiner::copy_paths_checked;
    using GraphCombiner::connect_seam;
    using GraphCombiner::merge_path_fragments;
    using GraphCombiner::dest;
    using GraphCombiner::max_node_id;
    using GraphCombiner::seams;
};

/// Spell out the sequence a path visits, so tests can assert on what a path
/// means rather than on which nodes happen to carry it.
static string spell_path(const PathHandleGraph& graph, const string& name) {
    string sequence;
    path_handle_t path = graph.get_path_handle(name);
    for (handle_t h : graph.scan_path(path)) {
        sequence += graph.get_sequence(h);
    }
    return sequence;
}

/// Make a chunk holding an empty HashGraph.
static TestGraphCombiner::Chunk make_chunk(const string& name) {
    TestGraphCombiner::Chunk chunk;
    chunk.name = name;
    chunk.graph = make_unique<bdsg::HashGraph>();
    return chunk;
}

/// Create a path with full metadata and walk it over the given handles.
static path_handle_t add_path(MutablePathHandleGraph& graph, const PathSense& sense,
                              const string& sample, const string& locus, size_t haplotype,
                              size_t phase_block, const subrange_t& subrange,
                              const vector<handle_t>& steps, bool is_circular = false) {
    path_handle_t path = graph.create_path(sense, sample, locus, haplotype, phase_block,
                                           subrange, is_circular);
    for (handle_t h : steps) {
        graph.append_step(path, h);
    }
    return path;
}

/// Create a REFERENCE-sense path, the sense --seam trim positions a chunk by.
static path_handle_t add_ref_path(MutablePathHandleGraph& graph, const string& sample,
                                  const string& locus, const subrange_t& subrange,
                                  const vector<handle_t>& steps) {
    return add_path(graph, PathSense::REFERENCE, sample, locus, 0,
                    PathMetadata::NO_PHASE_BLOCK, subrange, steps);
}

/// Collect a graph's path names, sorted so tests do not depend on iteration order.
static vector<string> sorted_path_names(const PathHandleGraph& graph) {
    vector<string> names;
    graph.for_each_path_handle([&](const path_handle_t& p) {
        names.push_back(graph.get_path_name(p));
    });
    std::sort(names.begin(), names.end());
    return names;
}

TEST_CASE("compute_path_length adds up what a path visits", "[combine]") {
    bdsg::HashGraph graph;
    handle_t n1 = graph.create_handle("AAAA");
    handle_t n2 = graph.create_handle("CC");
    handle_t n3 = graph.create_handle("GGGGG");
    graph.create_edge(n1, n2);
    graph.create_edge(n2, n3);
    add_ref_path(graph, "GRCh38", "chr1", subrange_t{0, 11}, {n1, n2, n3});

    SECTION("a graph without positions is scanned step by step") {
        // HashGraph is not a PathPositionHandleGraph, so this takes the scan_path branch.
        path_handle_t path = graph.get_path_handle("GRCh38#0#chr1[0-11]");
        REQUIRE(TestGraphCombiner::compute_path_length(&graph, path) == 11);
    }

    SECTION("a positional graph is asked directly, and agrees") {
        bdsg::PositionOverlay positional(&graph);
        path_handle_t path = positional.get_path_handle("GRCh38#0#chr1[0-11]");
        REQUIRE(TestGraphCombiner::compute_path_length(&positional, path) == 11);
    }

    SECTION("an empty path has no length") {
        graph.create_path_handle("empty");
        path_handle_t path = graph.get_path_handle("empty");
        REQUIRE(TestGraphCombiner::compute_path_length(&graph, path) == 0);
    }

    SECTION("a step's orientation does not change how much it contributes") {
        handle_t n4 = graph.create_handle("TTT");
        path_handle_t path = graph.create_path_handle("flipped");
        graph.append_step(path, graph.flip(n4));
        REQUIRE(TestGraphCombiner::compute_path_length(&graph, path) == 3);
    }
}

TEST_CASE("path_start_offset reads a start offset off path metadata", "[combine]") {
    SECTION("a subrange gives the offset outright") {
        REQUIRE(TestGraphCombiner::path_start_offset(subrange_t{100, 200}, 7, true) == 100);
        REQUIRE(TestGraphCombiner::path_start_offset(subrange_t{100, 200}, 7, false) == 100);
    }

    SECTION("a subrange wins over a phase block even when blocks are offsets") {
        REQUIRE(TestGraphCombiner::path_start_offset(subrange_t{100, PathMetadata::NO_END_POSITION},
                                                     7, true) == 100);
    }

    SECTION("without a subrange, -P reads the phase block as the offset") {
        REQUIRE(TestGraphCombiner::path_start_offset(PathMetadata::NO_SUBRANGE, 7, true) == 7);
    }

    SECTION("without -P a phase block is an opaque identifier, not a position") {
        REQUIRE(TestGraphCombiner::path_start_offset(PathMetadata::NO_SUBRANGE, 7, false) == 0);
    }

    SECTION("no subrange and no phase block starts at zero") {
        REQUIRE(TestGraphCombiner::path_start_offset(PathMetadata::NO_SUBRANGE,
                                                     PathMetadata::NO_PHASE_BLOCK, true) == 0);
        REQUIRE(TestGraphCombiner::path_start_offset(PathMetadata::NO_SUBRANGE,
                                                     PathMetadata::NO_PHASE_BLOCK, false) == 0);
    }
}

TEST_CASE("renumber_out_of_way shifts node IDs so they do not collide with IDs already used", "[combine]") {
    SECTION("IDs already above the mark are left alone") {
        bdsg::HashGraph graph;
        graph.create_handle("A", 10);
        graph.create_handle("C", 12);
        REQUIRE(TestGraphCombiner::renumber_out_of_way(&graph, 3) == 0);
        REQUIRE(graph.min_node_id() == 10);
        REQUIRE(graph.max_node_id() == 12);
    }

    SECTION("overlapping IDs are shifted clear, and stay clear") {
        bdsg::HashGraph graph;
        graph.create_handle("A", 1);
        graph.create_handle("C", 5);
        int64_t shift = TestGraphCombiner::renumber_out_of_way(&graph, 10);
        REQUIRE(shift == 10);
        REQUIRE(graph.min_node_id() == 11);
        REQUIRE(graph.max_node_id() == 15);
    }

    SECTION("a graph starting exactly on the mark still moves off it") {
        // delta is 0 here, which is not negative, so the shift is 1: IDs have to end
        // up strictly above the mark, not on it.
        bdsg::HashGraph graph;
        graph.create_handle("A", 7);
        REQUIRE(TestGraphCombiner::renumber_out_of_way(&graph, 7) == 1);
        REQUIRE(graph.min_node_id() == 8);
    }

    SECTION("the shift the caller gets back is the one that was applied") {
        bdsg::HashGraph graph;
        graph.create_handle("A", 4);
        graph.create_handle("C", 9);
        int64_t shift = TestGraphCombiner::renumber_out_of_way(&graph, 20);
        REQUIRE(graph.get_id(graph.get_handle(4 + shift)) == 4 + shift);
        REQUIRE(graph.has_node(9 + shift));
        REQUIRE(!graph.has_node(4));
    }
}

TEST_CASE("group_key_of decides which paths are pieces of the same thing", "[combine]") {
    bdsg::HashGraph graph;
    handle_t n1 = graph.create_handle("AAAA");

    add_path(graph, PathSense::HAPLOTYPE, "HG002", "chr1", 1, 0,
             PathMetadata::NO_SUBRANGE, {n1});
    add_path(graph, PathSense::HAPLOTYPE, "HG002", "chr1", 1, 10,
             PathMetadata::NO_SUBRANGE, {n1});

    path_handle_t block_0 = graph.get_path_handle("HG002#1#chr1#0");
    path_handle_t block_10 = graph.get_path_handle("HG002#1#chr1#10");

    SECTION("phase blocks tell paths apart by default") {
        REQUIRE(TestGraphCombiner::group_key_of(graph, block_0, false)
                != TestGraphCombiner::group_key_of(graph, block_10, false));
    }

    SECTION("-P folds the phase block out, so the two become one group") {
        REQUIRE(TestGraphCombiner::group_key_of(graph, block_0, true)
                == TestGraphCombiner::group_key_of(graph, block_10, true));
        // Folded out means replaced by the sentinel, not kept as some other value.
        REQUIRE(std::get<5>(TestGraphCombiner::group_key_of(graph, block_0, true))
                == PathMetadata::NO_PHASE_BLOCK);
    }

    SECTION("a subrange is not part of the key, so fragments group together") {
        add_path(graph, PathSense::REFERENCE, "GRCh38", "chr1", 0,
                 PathMetadata::NO_PHASE_BLOCK, subrange_t{0, 4}, {n1});
        add_path(graph, PathSense::REFERENCE, "GRCh38", "chr1", 0,
                 PathMetadata::NO_PHASE_BLOCK, subrange_t{4, 8}, {n1});
        REQUIRE(TestGraphCombiner::group_key_of(graph, graph.get_path_handle("GRCh38#0#chr1[0-4]"), false)
                == TestGraphCombiner::group_key_of(graph, graph.get_path_handle("GRCh38#0#chr1[4-8]"), false));
    }

    SECTION("each identity field on its own is enough to split the group") {
        add_path(graph, PathSense::HAPLOTYPE, "HG002", "chr1", 2, 0,
                 PathMetadata::NO_SUBRANGE, {n1});           // other haplotype
        add_path(graph, PathSense::HAPLOTYPE, "HG002", "chr2", 1, 0,
                 PathMetadata::NO_SUBRANGE, {n1});           // other locus
        add_path(graph, PathSense::HAPLOTYPE, "HG005", "chr1", 1, 0,
                 PathMetadata::NO_SUBRANGE, {n1});           // other sample
        add_path(graph, PathSense::REFERENCE, "HG002", "chr1", 1,
                 PathMetadata::NO_PHASE_BLOCK, subrange_t{0, 4}, {n1});  // other sense

        auto key = [&](const string& name) {
            return TestGraphCombiner::group_key_of(graph, graph.get_path_handle(name), false);
        };
        REQUIRE(key("HG002#1#chr1#0") != key("HG002#2#chr1#0"));
        REQUIRE(key("HG002#1#chr1#0") != key("HG002#1#chr2#0"));
        REQUIRE(key("HG002#1#chr1#0") != key("HG005#1#chr1#0"));
        REQUIRE(key("HG002#1#chr1#0") != key("HG002#1#chr1[0-4]"));
    }
}

TEST_CASE("path_prefix and path_suffix read the ends of a path", "[combine]") {
    bdsg::HashGraph graph;
    handle_t n1 = graph.create_handle("AAAA");
    handle_t n2 = graph.create_handle("CC");
    handle_t n3 = graph.create_handle("GGGTT");
    graph.create_edge(n1, n2);
    graph.create_edge(n2, n3);
    path_handle_t path = add_ref_path(graph, "GRCh38", "chr1", subrange_t{0, 11},
                                      {n1, n2, n3});

    SECTION("a span inside the first or last node stops there") {
        REQUIRE(TestGraphCombiner::path_prefix(graph, path, 3) == "AAA");
        REQUIRE(TestGraphCombiner::path_suffix(graph, path, 3) == "GTT");
    }

    SECTION("a span crossing node boundaries keeps reading") {
        REQUIRE(TestGraphCombiner::path_prefix(graph, path, 7) == "AAAACCG");
        REQUIRE(TestGraphCombiner::path_suffix(graph, path, 7) == "CCGGGTT");
    }

    SECTION("asking for more than the path holds yields all of it") {
        REQUIRE(TestGraphCombiner::path_prefix(graph, path, 20) == "AAAACCGGGTT");
        REQUIRE(TestGraphCombiner::path_suffix(graph, path, 20) == "AAAACCGGGTT");
    }

    SECTION("asking for nothing yields nothing, even on an empty path") {
        path_handle_t empty = add_ref_path(graph, "GRCh38", "chr2", subrange_t{0, 0}, {});
        REQUIRE(TestGraphCombiner::path_prefix(graph, path, 0) == "");
        REQUIRE(TestGraphCombiner::path_suffix(graph, path, 0) == "");
        REQUIRE(TestGraphCombiner::path_prefix(graph, empty, 5) == "");
        REQUIRE(TestGraphCombiner::path_suffix(graph, empty, 5) == "");
    }

    SECTION("a reverse step reads in the path's orientation") {
        handle_t n4 = graph.create_handle("TTGCA");
        graph.create_edge(n3, graph.flip(n4));
        path_handle_t flipped = add_ref_path(graph, "GRCh38", "chr3", subrange_t{0, 16},
                                             {n1, n2, n3, graph.flip(n4)});
        REQUIRE(TestGraphCombiner::path_suffix(graph, flipped, 6) == "TTGCAA");
    }
}

TEST_CASE("check_no_split_paths passes graphs whose paths are each whole", "[combine]") {
    // Paths that merely look alike are not pieces of one path.
    // Rejection exits the process, so it is tested in test/t/51_vg_combine.t.
    TestGraphCombiner combiner(SeamPolicy::RENUMBER, false, false);

    SECTION("a single path is not a split path") {
        TestGraphCombiner::Chunk chunk = make_chunk("one.pg");
        handle_t n1 = chunk.graph->create_handle("AAAA");
        add_ref_path(*chunk.graph, "GRCh38", "chr1", subrange_t{0, 4}, {n1});
        combiner.check_no_split_paths(chunk);
        REQUIRE(chunk.graph->get_path_count() == 1);
    }

    SECTION("paths differing only by haplotype are separate paths") {
        TestGraphCombiner::Chunk chunk = make_chunk("haps.pg");
        handle_t n1 = chunk.graph->create_handle("AAAA");
        add_path(*chunk.graph, PathSense::HAPLOTYPE, "HG002", "chr1", 1, 0,
                 subrange_t{0, 4}, {n1});
        add_path(*chunk.graph, PathSense::HAPLOTYPE, "HG002", "chr1", 2, 0,
                 subrange_t{0, 4}, {n1});
        combiner.check_no_split_paths(chunk);
        REQUIRE(chunk.graph->get_path_count() == 2);
    }

    SECTION("paths differing only by locus are separate paths") {
        TestGraphCombiner::Chunk chunk = make_chunk("loci.pg");
        handle_t n1 = chunk.graph->create_handle("AAAA");
        add_ref_path(*chunk.graph, "GRCh38", "chr1", subrange_t{0, 4}, {n1});
        add_ref_path(*chunk.graph, "GRCh38", "chr2", subrange_t{0, 4}, {n1});
        combiner.check_no_split_paths(chunk);
        REQUIRE(chunk.graph->get_path_count() == 2);
    }

    SECTION("phase blocks are separate paths until -P says otherwise") {
        TestGraphCombiner::Chunk chunk = make_chunk("blocks.pg");
        handle_t n1 = chunk.graph->create_handle("AAAA");
        add_path(*chunk.graph, PathSense::HAPLOTYPE, "HG002", "chr1", 1, 0,
                 PathMetadata::NO_SUBRANGE, {n1});
        add_path(*chunk.graph, PathSense::HAPLOTYPE, "HG002", "chr1", 1, 10,
                 PathMetadata::NO_SUBRANGE, {n1});
        combiner.check_no_split_paths(chunk);
        REQUIRE(chunk.graph->get_path_count() == 2);
    }
}

TEST_CASE("describe_reference reads a chunk's position off its REFERENCE path", "[combine]") {

    SECTION("every field is filled in from the path's metadata") {
        TestGraphCombiner combiner(SeamPolicy::COORD_TRIM, false, false);
        TestGraphCombiner::Chunk chunk = make_chunk("left.pg");
        handle_t n1 = chunk.graph->create_handle("AAAA");
        handle_t n2 = chunk.graph->create_handle("CC");
        handle_t n3 = chunk.graph->create_handle("GGGGG");
        chunk.graph->create_edge(n1, n2);
        chunk.graph->create_edge(n2, n3);
        add_ref_path(*chunk.graph, "GRCh38", "chr1", subrange_t{100, 111}, {n1, n2, n3});
        // A haplotype alongside it must not be mistaken for the reference.
        add_path(*chunk.graph, PathSense::HAPLOTYPE, "HG002", "chr1", 1, 0,
                 subrange_t{100, 111}, {n1, n2, n3});

        combiner.describe_reference(chunk);

        REQUIRE(chunk.ref_name == "GRCh38#0#chr1[100-111]");
        REQUIRE(chunk.ref_key == TestGraphCombiner::RefKey{"GRCh38", "chr1", 0, false});
        REQUIRE(chunk.ref_start_offset == 100);
        // The end is the start plus the path's own length, not the subrange's end.
        REQUIRE(chunk.ref_end_offset == 111);
        REQUIRE(chunk.graph->get_id(chunk.ref_start_handle) == chunk.graph->get_id(n1));
        REQUIRE(chunk.graph->get_id(chunk.ref_end_handle) == chunk.graph->get_id(n3));
        REQUIRE(!chunk.graph->get_is_reverse(chunk.ref_start_handle));
        REQUIRE(!chunk.graph->get_is_reverse(chunk.ref_end_handle));
    }

    SECTION("the end offset follows the sequence, not the subrange's claim") {
        // A subrange claiming 100-999 over 11bp of sequence still ends at 111.
        TestGraphCombiner combiner(SeamPolicy::COORD_TRIM, false, false);
        TestGraphCombiner::Chunk chunk = make_chunk("claim.pg");
        handle_t n1 = chunk.graph->create_handle("AAAAAAAAAAA");
        add_ref_path(*chunk.graph, "GRCh38", "chr1", subrange_t{100, 999}, {n1});

        combiner.describe_reference(chunk);
        REQUIRE(chunk.ref_start_offset == 100);
        REQUIRE(chunk.ref_end_offset == 111);
    }

    SECTION("a single-node chunk starts and ends on the same node") {
        TestGraphCombiner combiner(SeamPolicy::COORD_TRIM, false, false);
        TestGraphCombiner::Chunk chunk = make_chunk("solo.pg");
        handle_t n1 = chunk.graph->create_handle("AAAA");
        add_ref_path(*chunk.graph, "GRCh38", "chr1", subrange_t{0, 4}, {n1});

        combiner.describe_reference(chunk);
        REQUIRE(chunk.ref_start_handle == chunk.ref_end_handle);
        REQUIRE(chunk.ref_start_offset == 0);
        REQUIRE(chunk.ref_end_offset == 4);
    }

    SECTION("a reference with no subrange starts at zero, with or without -P") {
        // A REFERENCE-sense path may not carry a phase block, so -P has nothing to
        // read here and the chunk can only be placed at the origin.
        for (bool phase_block_is_offset : {false, true}) {
            TestGraphCombiner combiner(SeamPolicy::COORD_TRIM, false, phase_block_is_offset);
            TestGraphCombiner::Chunk chunk = make_chunk("bare.pg");
            handle_t n1 = chunk.graph->create_handle("AAAA");
            add_ref_path(*chunk.graph, "GRCh38", "chr1", PathMetadata::NO_SUBRANGE, {n1});

            combiner.describe_reference(chunk);
            REQUIRE(chunk.ref_start_offset == 0);
            REQUIRE(chunk.ref_end_offset == 4);
        }
    }

    SECTION("a circular reference is a different reference from a linear one") {
        TestGraphCombiner combiner(SeamPolicy::COORD_TRIM, false, false);
        TestGraphCombiner::Chunk chunk = make_chunk("circle.pg");
        handle_t n1 = chunk.graph->create_handle("AAAA");
        add_path(*chunk.graph, PathSense::REFERENCE, "GRCh38", "chrM", 0,
                 PathMetadata::NO_PHASE_BLOCK, subrange_t{0, 4}, {n1}, true);

        combiner.describe_reference(chunk);
        REQUIRE(chunk.ref_key == TestGraphCombiner::RefKey{"GRCh38", "chrM", 0, true});
    }
}

TEST_CASE("sort_by_reference_offset groups references and orders each by offset", "[combine]") {
    TestGraphCombiner combiner(SeamPolicy::COORD_TRIM, false, false);

    // Build chunks by hand: the sort reads only ref_key and ref_start_offset, so
    // these need no graphs.
    auto positioned = [](const string& name, const string& locus, offset_t start) {
        TestGraphCombiner::Chunk chunk;
        chunk.name = name;
        chunk.ref_key = TestGraphCombiner::RefKey{"GRCh38", locus, 0, false};
        chunk.ref_start_offset = start;
        return chunk;
    };
    // Same, for chunks whose reference differs in something other than locus.
    auto keyed = [](const string& name, const TestGraphCombiner::RefKey& ref_key,
                    offset_t start) {
        TestGraphCombiner::Chunk chunk;
        chunk.name = name;
        chunk.ref_key = ref_key;
        chunk.ref_start_offset = start;
        return chunk;
    };

    SECTION("chunks of one reference come out in ascending offset order") {
        vector<TestGraphCombiner::Chunk> chunks;
        chunks.push_back(positioned("c", "chr1", 200));
        chunks.push_back(positioned("a", "chr1", 0));
        chunks.push_back(positioned("b", "chr1", 100));

        combiner.sort_by_reference_offset(chunks);
        REQUIRE(chunks[0].name == "a");
        REQUIRE(chunks[1].name == "b");
        REQUIRE(chunks[2].name == "c");
    }

    SECTION("references keep first-seen order, and their chunks stay contiguous") {
        // chr2 is seen first, so it stays first even though chr1 sorts earlier.
        vector<TestGraphCombiner::Chunk> chunks;
        chunks.push_back(positioned("y2", "chr2", 100));
        chunks.push_back(positioned("x2", "chr1", 100));
        chunks.push_back(positioned("y1", "chr2", 0));
        chunks.push_back(positioned("x1", "chr1", 0));

        combiner.sort_by_reference_offset(chunks);
        REQUIRE(chunks[0].name == "y1");
        REQUIRE(chunks[1].name == "y2");
        REQUIRE(chunks[2].name == "x1");
        REQUIRE(chunks[3].name == "x2");
    }

    SECTION("equal offsets keep the order they arrived in") {
        vector<TestGraphCombiner::Chunk> chunks;
        chunks.push_back(positioned("first", "chr1", 50));
        chunks.push_back(positioned("second", "chr1", 50));
        chunks.push_back(positioned("third", "chr1", 50));

        combiner.sort_by_reference_offset(chunks);
        REQUIRE(chunks[0].name == "first");
        REQUIRE(chunks[1].name == "second");
        REQUIRE(chunks[2].name == "third");
    }

    SECTION("haplotype and circularity make separate references of one locus") {
        // All four share sample and locus, so only the other two RefKey fields can
        // tell the three references apart.
        vector<TestGraphCombiner::Chunk> chunks;
        chunks.push_back(keyed("hap2", {"GRCh38", "chr1", 2, false}, 50));
        chunks.push_back(positioned("hap0_late", "chr1", 100));
        chunks.push_back(positioned("hap0_early", "chr1", 0));
        chunks.push_back(keyed("circular", {"GRCh38", "chr1", 0, true}, 50));

        combiner.sort_by_reference_offset(chunks);
        // References lead in first-seen order; hap0's two chunks follow in offset order.
        REQUIRE(chunks[0].name == "hap2");
        REQUIRE(chunks[1].name == "hap0_early");
        REQUIRE(chunks[2].name == "hap0_late");
        REQUIRE(chunks[3].name == "circular");
    }
}

TEST_CASE("check_all_paths_start_at and check_all_paths_end_at pass chunks with flush boundaries", "[combine]") {
    // Legitimate chunks are not turned away.
    // Rejection exits the process, so it is tested in test/t/51_vg_combine.t.
    TestGraphCombiner combiner(SeamPolicy::COORD_TRIM, false, false);

    TestGraphCombiner::Chunk chunk = make_chunk("flush.pg");
    handle_t n1 = chunk.graph->create_handle("AAAA");
    handle_t n2 = chunk.graph->create_handle("C");
    handle_t n3 = chunk.graph->create_handle("G");
    handle_t n4 = chunk.graph->create_handle("TTTTT");
    chunk.graph->create_edge(n1, n2);
    chunk.graph->create_edge(n1, n3);
    chunk.graph->create_edge(n2, n4);
    chunk.graph->create_edge(n3, n4);
    // Three paths through a bubble: they take different alleles but share both ends.
    add_ref_path(*chunk.graph, "GRCh38", "chr1", subrange_t{0, 10}, {n1, n2, n4});
    add_path(*chunk.graph, PathSense::HAPLOTYPE, "HG002", "chr1", 1, 0,
             subrange_t{0, 10}, {n1, n3, n4});
    add_path(*chunk.graph, PathSense::HAPLOTYPE, "HG002", "chr1", 2, 0,
             subrange_t{0, 10}, {n1, n2, n4});
    chunk.ref_start_handle = n1;
    chunk.ref_end_handle = n4;

    SECTION("all paths beginning on the boundary node are accepted") {
        combiner.check_all_paths_start_at(chunk);
        REQUIRE(chunk.graph->get_path_count() == 3);
    }

    SECTION("all paths ending on the boundary node are accepted") {
        combiner.check_all_paths_end_at(chunk);
        REQUIRE(chunk.graph->get_path_count() == 3);
    }

    SECTION("a boundary the paths cross in reverse is matched by orientation too") {
        // The check compares handles, so a reverse boundary needs reverse steps.
        TestGraphCombiner::Chunk flipped = make_chunk("flipped.pg");
        handle_t a = flipped.graph->create_handle("AAAA");
        handle_t b = flipped.graph->create_handle("CCC");
        flipped.graph->create_edge(flipped.graph->flip(a), b);
        add_ref_path(*flipped.graph, "GRCh38", "chr1", subrange_t{0, 7},
                     {flipped.graph->flip(a), b});
        flipped.ref_start_handle = flipped.graph->flip(a);
        flipped.ref_end_handle = b;

        combiner.check_all_paths_start_at(flipped);
        combiner.check_all_paths_end_at(flipped);
        REQUIRE(flipped.graph->get_path_count() == 1);
    }
}

TEST_CASE("trim_and_shift cuts the overlap away and moves coordinates with it", "[combine]") {

    SECTION("no overlap and no fuse leaves the chunk exactly as it was") {
        TestGraphCombiner combiner(SeamPolicy::COORD_TRIM, false, false);
        TestGraphCombiner::Chunk chunk = make_chunk("clean.pg");
        handle_t n1 = chunk.graph->create_handle("AAAA");
        handle_t n2 = chunk.graph->create_handle("CCC");
        chunk.graph->create_edge(n1, n2);
        add_ref_path(*chunk.graph, "GRCh38", "chr1", subrange_t{100, 107}, {n1, n2});
        chunk.ref_start_handle = n1;
        nid_t before = chunk.graph->get_id(n1);

        nid_t after = combiner.trim_and_shift(chunk, 0, 5);

        REQUIRE(after == before);
        REQUIRE(chunk.graph->get_node_count() == 2);
        REQUIRE(sorted_path_names(*chunk.graph) == vector<string>{"GRCh38#0#chr1[100-107]"});
    }

    SECTION("no overlap but fusing pulls every coordinate back by the welded node") {
        // Welding buries the left chunk's last node inside this chunk's first, so
        // this chunk's coordinates move back by that node's length.
        TestGraphCombiner combiner(SeamPolicy::COORD_TRIM, true, false);
        TestGraphCombiner::Chunk chunk = make_chunk("weld.pg");
        handle_t n1 = chunk.graph->create_handle("AAAA");
        handle_t n2 = chunk.graph->create_handle("CCC");
        chunk.graph->create_edge(n1, n2);
        add_ref_path(*chunk.graph, "GRCh38", "chr1", subrange_t{100, 107}, {n1, n2});
        chunk.ref_start_handle = n1;

        combiner.trim_and_shift(chunk, 0, 5);

        REQUIRE(chunk.graph->get_node_count() == 2);
        REQUIRE(sorted_path_names(*chunk.graph) == vector<string>{"GRCh38#0#chr1[95-107]"});
    }

    SECTION("an overlap splits the first node, drops the front, and shifts coordinates") {
        TestGraphCombiner combiner(SeamPolicy::COORD_TRIM, false, false);
        TestGraphCombiner::Chunk chunk = make_chunk("overlap.pg");
        handle_t n1 = chunk.graph->create_handle("AAAACC");
        handle_t n2 = chunk.graph->create_handle("GGG");
        chunk.graph->create_edge(n1, n2);
        add_ref_path(*chunk.graph, "GRCh38", "chr1", subrange_t{100, 109}, {n1, n2});
        add_path(*chunk.graph, PathSense::HAPLOTYPE, "HG002", "chr1", 1, 0,
                 subrange_t{100, 109}, {n1, n2});
        chunk.ref_start_handle = n1;

        nid_t new_start = combiner.trim_and_shift(chunk, 4, 5);

        // The overlapping "AAAA" is gone; "CC" remains and is the new start.
        REQUIRE(chunk.graph->get_node_count() == 2);
        REQUIRE(chunk.graph->get_sequence(chunk.graph->get_handle(new_start)) == "CC");
        REQUIRE(spell_path(*chunk.graph, "GRCh38#0#chr1[104-109]") == "CCGGG");
        REQUIRE(spell_path(*chunk.graph, "HG002#1#chr1#0[104-109]") == "CCGGG");
    }

    SECTION("an overlap of one base still leaves the rest of the node behind") {
        TestGraphCombiner combiner(SeamPolicy::COORD_TRIM, false, false);
        TestGraphCombiner::Chunk chunk = make_chunk("thin.pg");
        handle_t n1 = chunk.graph->create_handle("AC");
        add_ref_path(*chunk.graph, "GRCh38", "chr1", subrange_t{100, 102}, {n1});
        chunk.ref_start_handle = n1;

        nid_t new_start = combiner.trim_and_shift(chunk, 1, 5);

        REQUIRE(chunk.graph->get_node_count() == 1);
        REQUIRE(chunk.graph->get_sequence(chunk.graph->get_handle(new_start)) == "C");
        REQUIRE(spell_path(*chunk.graph, "GRCh38#0#chr1[101-102]") == "C");
    }

    SECTION("with -P the phase block moves instead of a subrange") {
        TestGraphCombiner combiner(SeamPolicy::COORD_TRIM, false, true);
        TestGraphCombiner::Chunk chunk = make_chunk("pblock.pg");
        handle_t n1 = chunk.graph->create_handle("AAAACC");
        add_path(*chunk.graph, PathSense::HAPLOTYPE, "HG002", "chr1", 1, 100,
                 PathMetadata::NO_SUBRANGE, {n1});
        chunk.ref_start_handle = n1;

        combiner.trim_and_shift(chunk, 4, 5);

        REQUIRE(sorted_path_names(*chunk.graph) == vector<string>{"HG002#1#chr1#104"});
    }

    SECTION("a path with no coordinates at all is given a subrange to carry the shift") {
        TestGraphCombiner combiner(SeamPolicy::COORD_TRIM, false, false);
        TestGraphCombiner::Chunk chunk = make_chunk("bare.pg");
        handle_t n1 = chunk.graph->create_handle("AAAACC");
        add_path(*chunk.graph, PathSense::GENERIC, PathMetadata::NO_SAMPLE_NAME, "scaffold",
                 PathMetadata::NO_HAPLOTYPE, PathMetadata::NO_PHASE_BLOCK,
                 PathMetadata::NO_SUBRANGE, {n1});
        chunk.ref_start_handle = n1;

        combiner.trim_and_shift(chunk, 4, 0);

        REQUIRE(sorted_path_names(*chunk.graph) == vector<string>{"scaffold[4]"});
    }

    SECTION("a shift that would run past zero stops at zero") {
        // Fusing pulls coordinates back by 50, but the path only starts at 10.
        TestGraphCombiner combiner(SeamPolicy::COORD_TRIM, true, false);
        TestGraphCombiner::Chunk chunk = make_chunk("clamp.pg");
        handle_t n1 = chunk.graph->create_handle("AAAA");
        add_ref_path(*chunk.graph, "GRCh38", "chr1", subrange_t{10, 14}, {n1});
        chunk.ref_start_handle = n1;

        combiner.trim_and_shift(chunk, 0, 50);

        REQUIRE(sorted_path_names(*chunk.graph) == vector<string>{"GRCh38#0#chr1[0-14]"});
    }

    SECTION("trimming keeps each path on its own alleles through a bubble") {
        TestGraphCombiner combiner(SeamPolicy::COORD_TRIM, false, false);
        TestGraphCombiner::Chunk chunk = make_chunk("bubble.pg");
        handle_t n1 = chunk.graph->create_handle("AAAACC");
        handle_t n2 = chunk.graph->create_handle("G");
        handle_t n3 = chunk.graph->create_handle("T");
        handle_t n4 = chunk.graph->create_handle("TTTTT");
        chunk.graph->create_edge(n1, n2);
        chunk.graph->create_edge(n1, n3);
        chunk.graph->create_edge(n2, n4);
        chunk.graph->create_edge(n3, n4);
        add_ref_path(*chunk.graph, "GRCh38", "chr1", subrange_t{100, 112}, {n1, n2, n4});
        add_path(*chunk.graph, PathSense::HAPLOTYPE, "HG002", "chr1", 1, 0,
                 subrange_t{100, 112}, {n1, n3, n4});
        chunk.ref_start_handle = n1;

        combiner.trim_and_shift(chunk, 4, 0);

        REQUIRE(spell_path(*chunk.graph, "GRCh38#0#chr1[104-112]") == "CCGTTTTT");
        REQUIRE(spell_path(*chunk.graph, "HG002#1#chr1#0[104-112]") == "CCTTTTTT");
    }
}

TEST_CASE("ingest moves a chunk into the accumulator", "[combine]") {

    SECTION("the first chunk becomes the accumulator untouched") {
        TestGraphCombiner combiner(SeamPolicy::RENUMBER, false, false);
        TestGraphCombiner::Chunk chunk = make_chunk("first.pg");
        chunk.graph->create_handle("AAAA", 5);
        chunk.graph->create_handle("CCC", 9);

        int64_t shift = combiner.ingest(chunk);

        REQUIRE(shift == 0);
        REQUIRE(combiner.max_node_id == 9);
        REQUIRE(combiner.dest != nullptr);
        REQUIRE(combiner.dest->has_node(5));
        REQUIRE(combiner.dest->has_node(9));
        // The graph was moved out of the chunk, not copied.
        REQUIRE(chunk.graph == nullptr);
    }

    SECTION("renumber lifts a second chunk's IDs clear of the first") {
        TestGraphCombiner combiner(SeamPolicy::RENUMBER, false, false);
        TestGraphCombiner::Chunk left = make_chunk("left.pg");
        handle_t l1 = left.graph->create_handle("AAAA", 1);
        add_ref_path(*left.graph, "GRCh38", "chr1", subrange_t{0, 4}, {l1});
        combiner.ingest(left);

        TestGraphCombiner::Chunk right = make_chunk("right.pg");
        handle_t r1 = right.graph->create_handle("CCC", 1);
        handle_t r2 = right.graph->create_handle("GG", 2);
        right.graph->create_edge(r1, r2);
        add_ref_path(*right.graph, "GRCh38", "chr1", subrange_t{4, 9}, {r1, r2});

        int64_t shift = combiner.ingest(right);

        REQUIRE(shift > 0);
        // Both chunks' nodes are present and none was overwritten.
        REQUIRE(combiner.dest->get_node_count() == 3);
        REQUIRE(combiner.dest->get_sequence(combiner.dest->get_handle(1)) == "AAAA");
        REQUIRE(combiner.dest->get_sequence(combiner.dest->get_handle(1 + shift)) == "CCC");
        REQUIRE(combiner.dest->has_edge(combiner.dest->get_handle(1 + shift),
                                        combiner.dest->get_handle(2 + shift)));
        REQUIRE(combiner.max_node_id == 2 + shift);
        REQUIRE(sorted_path_names(*combiner.dest)
                == vector<string>{"GRCh38#0#chr1[0-4]", "GRCh38#0#chr1[4-9]"});
        // The source is released rather than held until combine finishes.
        REQUIRE(right.graph == nullptr);
    }

    SECTION("shared keeps IDs as they are, so nothing is shifted") {
        TestGraphCombiner combiner(SeamPolicy::SHARED_IDS, false, false);
        TestGraphCombiner::Chunk left = make_chunk("left.pg");
        handle_t l1 = left.graph->create_handle("AAAA", 1);
        handle_t l2 = left.graph->create_handle("CCC", 2);
        left.graph->create_edge(l1, l2);
        add_ref_path(*left.graph, "GRCh38", "chr1", subrange_t{0, 7}, {l1, l2});
        combiner.ingest(left);

        TestGraphCombiner::Chunk right = make_chunk("right.pg");
        handle_t r2 = right.graph->create_handle("CCC", 2);   // the shared node
        handle_t r3 = right.graph->create_handle("GG", 3);
        right.graph->create_edge(r2, r3);
        add_ref_path(*right.graph, "GRCh38", "chr1", subrange_t{4, 9}, {r2, r3});

        int64_t shift = combiner.ingest(right);

        REQUIRE(shift == 0);
        REQUIRE(combiner.dest->get_node_count() == 3);
        REQUIRE(combiner.max_node_id == 3);
        REQUIRE(spell_path(*combiner.dest, "GRCh38#0#chr1[4-9]") == "CCCGG");
    }
}

TEST_CASE("share_nodes_by_id treats a repeated ID as the same node", "[combine]") {
    TestGraphCombiner combiner(SeamPolicy::SHARED_IDS, false, false);
    combiner.dest = make_unique<bdsg::HashGraph>();
    handle_t d1 = combiner.dest->create_handle("AAAA", 1);
    handle_t d2 = combiner.dest->create_handle("CCC", 2);
    combiner.dest->create_edge(d1, d2);

    SECTION("an ID already present is reused, and a new one is created") {
        TestGraphCombiner::Chunk chunk = make_chunk("more.pg");
        chunk.graph->create_handle("CCC", 2);      // already in dest
        chunk.graph->create_handle("GG", 3);       // new

        combiner.share_nodes_by_id(chunk);

        REQUIRE(combiner.dest->get_node_count() == 3);
        REQUIRE(combiner.dest->get_sequence(combiner.dest->get_handle(2)) == "CCC");
        REQUIRE(combiner.dest->get_sequence(combiner.dest->get_handle(3)) == "GG");
    }

    SECTION("an edge that is already there is not created twice") {
        TestGraphCombiner::Chunk chunk = make_chunk("dup_edge.pg");
        handle_t c1 = chunk.graph->create_handle("AAAA", 1);
        handle_t c2 = chunk.graph->create_handle("CCC", 2);
        chunk.graph->create_edge(c1, c2);

        combiner.share_nodes_by_id(chunk);

        size_t edges = 0;
        combiner.dest->for_each_edge([&](const edge_t&) { edges++; });
        REQUIRE(edges == 1);
        REQUIRE(combiner.dest->get_node_count() == 2);
    }

    SECTION("an edge's orientation survives the copy") {
        TestGraphCombiner::Chunk chunk = make_chunk("flip_edge.pg");
        handle_t c2 = chunk.graph->create_handle("CCC", 2);
        handle_t c3 = chunk.graph->create_handle("GG", 3);
        // 2 to the reverse of 3: the copy must not quietly straighten this out.
        chunk.graph->create_edge(c2, chunk.graph->flip(c3));

        combiner.share_nodes_by_id(chunk);

        REQUIRE(combiner.dest->has_edge(combiner.dest->get_handle(2),
                                        combiner.dest->get_handle(3, true)));
        REQUIRE(!combiner.dest->has_edge(combiner.dest->get_handle(2),
                                         combiner.dest->get_handle(3, false)));
    }

    SECTION("nodes arriving on their own, with no edges, still land") {
        TestGraphCombiner::Chunk chunk = make_chunk("island.pg");
        chunk.graph->create_handle("TTTT", 7);

        combiner.share_nodes_by_id(chunk);

        REQUIRE(combiner.dest->has_node(7));
        REQUIRE(combiner.dest->get_sequence(combiner.dest->get_handle(7)) == "TTTT");
    }
}

TEST_CASE("copy_paths_checked carries a path's metadata across with it", "[combine]") {
    TestGraphCombiner combiner(SeamPolicy::SHARED_IDS, false, false);
    combiner.dest = make_unique<bdsg::HashGraph>();
    handle_t d1 = combiner.dest->create_handle("AAAA", 1);
    handle_t d2 = combiner.dest->create_handle("CCC", 2);
    combiner.dest->create_edge(d1, d2);

    TestGraphCombiner::Chunk chunk = make_chunk("paths.pg");
    handle_t c1 = chunk.graph->create_handle("AAAA", 1);
    handle_t c2 = chunk.graph->create_handle("CCC", 2);
    chunk.graph->create_edge(c1, c2);
    add_ref_path(*chunk.graph, "GRCh38", "chr1", subrange_t{0, 7}, {c1, c2});
    add_path(*chunk.graph, PathSense::HAPLOTYPE, "HG002", "chr1", 1, 3,
             PathMetadata::NO_SUBRANGE, {c1, c2});
    add_path(*chunk.graph, PathSense::GENERIC, PathMetadata::NO_SAMPLE_NAME, "scaffold",
             PathMetadata::NO_HAPLOTYPE, PathMetadata::NO_PHASE_BLOCK,
             PathMetadata::NO_SUBRANGE, {c2}, true);

    combiner.copy_paths_checked(chunk);

    SECTION("every path arrives, under its own name") {
        REQUIRE(sorted_path_names(*combiner.dest)
                == vector<string>{"GRCh38#0#chr1[0-7]", "HG002#1#chr1#3", "scaffold"});
    }

    SECTION("the metadata fields come across one for one") {
        path_handle_t hap = combiner.dest->get_path_handle("HG002#1#chr1#3");
        REQUIRE(combiner.dest->get_sense(hap) == PathSense::HAPLOTYPE);
        REQUIRE(combiner.dest->get_sample_name(hap) == "HG002");
        REQUIRE(combiner.dest->get_locus_name(hap) == "chr1");
        REQUIRE(combiner.dest->get_haplotype(hap) == 1);
        REQUIRE(combiner.dest->get_phase_block(hap) == 3);

        path_handle_t ref = combiner.dest->get_path_handle("GRCh38#0#chr1[0-7]");
        REQUIRE(combiner.dest->get_sense(ref) == PathSense::REFERENCE);
        REQUIRE(combiner.dest->get_subrange(ref) == subrange_t{0, 7});

        REQUIRE(combiner.dest->get_is_circular(combiner.dest->get_path_handle("scaffold")));
    }

    SECTION("the steps land in order, on dest's own nodes") {
        REQUIRE(spell_path(*combiner.dest, "GRCh38#0#chr1[0-7]") == "AAAACCC");
        REQUIRE(combiner.dest->get_step_count(
                    combiner.dest->get_path_handle("GRCh38#0#chr1[0-7]")) == 2);
    }
}

TEST_CASE("connect_seam joins the running reference to the next chunk", "[combine]") {

    SECTION("without fusing, the two are joined by an edge") {
        TestGraphCombiner combiner(SeamPolicy::COORD_TRIM, false, false);
        combiner.dest = make_unique<bdsg::HashGraph>();
        handle_t left = combiner.dest->create_handle("AAAA", 1);
        handle_t start = combiner.dest->create_handle("CCC", 2);
        handle_t end = combiner.dest->create_handle("GG", 3);
        combiner.dest->create_edge(start, end);

        TestGraphCombiner::SeamState seam{combiner.dest->get_id(left), false, 4, "left.pg"};
        handle_t new_end = combiner.connect_seam(seam, "right.pg", start, end);

        REQUIRE(combiner.dest->has_edge(left, start));
        REQUIRE(combiner.dest->get_node_count() == 3);
        // Nothing was welded, so the chunk's own end is still the reference's end.
        REQUIRE(new_end == end);
    }

    SECTION("an edge that is already there is not added again") {
        TestGraphCombiner combiner(SeamPolicy::COORD_TRIM, false, false);
        combiner.dest = make_unique<bdsg::HashGraph>();
        handle_t left = combiner.dest->create_handle("AAAA", 1);
        handle_t start = combiner.dest->create_handle("CCC", 2);
        combiner.dest->create_edge(left, start);

        TestGraphCombiner::SeamState seam{combiner.dest->get_id(left), false, 4, "left.pg"};
        combiner.connect_seam(seam, "right.pg", start, start);

        size_t edges = 0;
        combiner.dest->for_each_edge([&](const edge_t&) { edges++; });
        REQUIRE(edges == 1);
    }

    SECTION("fusing welds the two boundary nodes into one") {
        TestGraphCombiner combiner(SeamPolicy::COORD_TRIM, true, false);
        combiner.dest = make_unique<bdsg::HashGraph>();
        handle_t before = combiner.dest->create_handle("TT", 1);
        handle_t left = combiner.dest->create_handle("AAAA", 2);
        handle_t start = combiner.dest->create_handle("CCC", 3);
        handle_t end = combiner.dest->create_handle("GG", 4);
        combiner.dest->create_edge(before, left);
        combiner.dest->create_edge(start, end);
        add_ref_path(*combiner.dest, "GRCh38", "chr1", subrange_t{0, 6}, {before, left});
        combiner.max_node_id = 4;

        TestGraphCombiner::SeamState seam{combiner.dest->get_id(left), false, 6, "left.pg"};
        handle_t new_end = combiner.connect_seam(seam, "right.pg", start, end);

        REQUIRE(combiner.dest->get_node_count() == 3);
        REQUIRE(spell_path(*combiner.dest, "GRCh38#0#chr1[0-6]") == "TTAAAACCC");
        REQUIRE(new_end == end);
        REQUIRE(combiner.max_node_id >= combiner.dest->max_node_id());
    }

    SECTION("fusing a single-node chunk returns the weld, not the stale handle") {
        // When the chunk's reference ends on the very node being welded, the old
        // handle is destroyed by the fuse; the caller must get the new one back.
        TestGraphCombiner combiner(SeamPolicy::COORD_TRIM, true, false);
        combiner.dest = make_unique<bdsg::HashGraph>();
        handle_t left = combiner.dest->create_handle("AAAA", 1);
        handle_t start = combiner.dest->create_handle("CCC", 2);
        combiner.max_node_id = 2;

        TestGraphCombiner::SeamState seam{combiner.dest->get_id(left), false, 4, "left.pg"};
        handle_t new_end = combiner.connect_seam(seam, "right.pg", start, start);

        REQUIRE(combiner.dest->get_node_count() == 1);
        REQUIRE(combiner.dest->has_node(combiner.dest->get_id(new_end)));
        REQUIRE(combiner.dest->get_sequence(new_end) == "AAAACCC");
        REQUIRE(combiner.max_node_id >= combiner.dest->max_node_id());
    }
}

TEST_CASE("merge_path_fragments joins pieces of one path back together", "[combine]") {

    SECTION("two abutting fragments become one path spanning both") {
        TestGraphCombiner combiner(SeamPolicy::RENUMBER, false, false);
        combiner.dest = make_unique<bdsg::HashGraph>();
        handle_t n1 = combiner.dest->create_handle("AAAA");
        handle_t n2 = combiner.dest->create_handle("CCC");
        add_ref_path(*combiner.dest, "GRCh38", "chr1", subrange_t{0, 4}, {n1});
        add_ref_path(*combiner.dest, "GRCh38", "chr1", subrange_t{4, 7}, {n2});

        combiner.merge_path_fragments();

        REQUIRE(sorted_path_names(*combiner.dest) == vector<string>{"GRCh38#0#chr1[0-7]"});
        REQUIRE(spell_path(*combiner.dest, "GRCh38#0#chr1[0-7]") == "AAAACCC");
        // Nothing joined the two before, so the merge has to say they abut.
        REQUIRE(combiner.dest->has_edge(n1, n2));
    }

    SECTION("a lone fragment is left alone, and keeps its name") {
        TestGraphCombiner combiner(SeamPolicy::RENUMBER, false, false);
        combiner.dest = make_unique<bdsg::HashGraph>();
        handle_t n1 = combiner.dest->create_handle("AAAA");
        add_ref_path(*combiner.dest, "GRCh38", "chr1", subrange_t{0, 4}, {n1});

        combiner.merge_path_fragments();

        REQUIRE(sorted_path_names(*combiner.dest) == vector<string>{"GRCh38#0#chr1[0-4]"});
    }

    SECTION("fragments sharing a run of steps at the seam do not repeat it") {
        TestGraphCombiner combiner(SeamPolicy::SHARED_IDS, false, false);
        combiner.dest = make_unique<bdsg::HashGraph>();
        handle_t n1 = combiner.dest->create_handle("AAAA");
        handle_t n2 = combiner.dest->create_handle("CC");
        handle_t n3 = combiner.dest->create_handle("GGG");
        combiner.dest->create_edge(n1, n2);
        combiner.dest->create_edge(n2, n3);
        // Both fragments cross n2: it must appear once in the merge, not twice.
        add_ref_path(*combiner.dest, "GRCh38", "chr1", subrange_t{0, 6}, {n1, n2});
        add_ref_path(*combiner.dest, "GRCh38", "chr1", subrange_t{4, 9}, {n2, n3});

        combiner.merge_path_fragments();

        REQUIRE(spell_path(*combiner.dest, "GRCh38#0#chr1[0-9]") == "AAAACCGGG");
        REQUIRE(combiner.dest->get_step_count(
                    combiner.dest->get_path_handle("GRCh38#0#chr1[0-9]")) == 3);
    }

    SECTION("fragments of different paths do not bleed into each other") {
        TestGraphCombiner combiner(SeamPolicy::RENUMBER, false, false);
        combiner.dest = make_unique<bdsg::HashGraph>();
        handle_t x1 = combiner.dest->create_handle("AAAA");
        handle_t x2 = combiner.dest->create_handle("CCC");
        handle_t y1 = combiner.dest->create_handle("TTTT");
        handle_t y2 = combiner.dest->create_handle("GGG");
        add_ref_path(*combiner.dest, "GRCh38", "chr1", subrange_t{0, 4}, {x1});
        add_ref_path(*combiner.dest, "GRCh38", "chr1", subrange_t{4, 7}, {x2});
        add_ref_path(*combiner.dest, "GRCh38", "chr2", subrange_t{0, 4}, {y1});
        add_ref_path(*combiner.dest, "GRCh38", "chr2", subrange_t{4, 7}, {y2});

        combiner.merge_path_fragments();

        REQUIRE(sorted_path_names(*combiner.dest)
                == vector<string>{"GRCh38#0#chr1[0-7]", "GRCh38#0#chr2[0-7]"});
        REQUIRE(spell_path(*combiner.dest, "GRCh38#0#chr1[0-7]") == "AAAACCC");
        REQUIRE(spell_path(*combiner.dest, "GRCh38#0#chr2[0-7]") == "TTTTGGG");
        REQUIRE(!combiner.dest->has_edge(x2, y1));
    }

    SECTION("fragments arriving out of order are put back in coordinate order") {
        TestGraphCombiner combiner(SeamPolicy::RENUMBER, false, false);
        combiner.dest = make_unique<bdsg::HashGraph>();
        handle_t n3 = combiner.dest->create_handle("GGG");
        handle_t n1 = combiner.dest->create_handle("AAAA");
        handle_t n2 = combiner.dest->create_handle("CC");
        add_ref_path(*combiner.dest, "GRCh38", "chr1", subrange_t{6, 9}, {n3});
        add_ref_path(*combiner.dest, "GRCh38", "chr1", subrange_t{0, 4}, {n1});
        add_ref_path(*combiner.dest, "GRCh38", "chr1", subrange_t{4, 6}, {n2});

        combiner.merge_path_fragments();

        REQUIRE(spell_path(*combiner.dest, "GRCh38#0#chr1[0-9]") == "AAAACCGGG");
    }

    SECTION("a merged HAPLOTYPE path keeps the phase block it needs") {
        TestGraphCombiner combiner(SeamPolicy::RENUMBER, false, false);
        combiner.dest = make_unique<bdsg::HashGraph>();
        handle_t n1 = combiner.dest->create_handle("AAAA");
        handle_t n2 = combiner.dest->create_handle("CCC");
        add_path(*combiner.dest, PathSense::HAPLOTYPE, "HG002", "chr22", 2, 4,
                 subrange_t{0, 4}, {n1});
        add_path(*combiner.dest, PathSense::HAPLOTYPE, "HG002", "chr22", 2, 4,
                 subrange_t{4, 7}, {n2});

        combiner.merge_path_fragments();

        REQUIRE(sorted_path_names(*combiner.dest) == vector<string>{"HG002#2#chr22#4[0-7]"});
        path_handle_t merged = combiner.dest->get_path_handle("HG002#2#chr22#4[0-7]");
        REQUIRE(combiner.dest->get_phase_block(merged) == 4);
        REQUIRE(combiner.dest->get_sense(merged) == PathSense::HAPLOTYPE);
    }

    SECTION("with -P, phase blocks position the fragments and fold into a subrange") {
        TestGraphCombiner combiner(SeamPolicy::RENUMBER, false, true);
        combiner.dest = make_unique<bdsg::HashGraph>();
        handle_t n1 = combiner.dest->create_handle("AAAA");
        handle_t n2 = combiner.dest->create_handle("CCC");
        add_path(*combiner.dest, PathSense::HAPLOTYPE, "HG002", "chr1", 1, 0,
                 PathMetadata::NO_SUBRANGE, {n1});
        add_path(*combiner.dest, PathSense::HAPLOTYPE, "HG002", "chr1", 1, 4,
                 PathMetadata::NO_SUBRANGE, {n2});

        combiner.merge_path_fragments();

        REQUIRE(sorted_path_names(*combiner.dest) == vector<string>{"HG002#1#chr1#0[0-7]"});
        REQUIRE(spell_path(*combiner.dest, "HG002#1#chr1#0[0-7]") == "AAAACCC");
    }

    SECTION("a circular path stays circular through the merge") {
        TestGraphCombiner combiner(SeamPolicy::RENUMBER, false, false);
        combiner.dest = make_unique<bdsg::HashGraph>();
        handle_t n1 = combiner.dest->create_handle("AAAA");
        handle_t n2 = combiner.dest->create_handle("CCC");
        add_path(*combiner.dest, PathSense::REFERENCE, "GRCh38", "chrM", 0,
                 PathMetadata::NO_PHASE_BLOCK, subrange_t{0, 4}, {n1}, true);
        add_path(*combiner.dest, PathSense::REFERENCE, "GRCh38", "chrM", 0,
                 PathMetadata::NO_PHASE_BLOCK, subrange_t{4, 7}, {n2}, true);

        combiner.merge_path_fragments();

        REQUIRE(sorted_path_names(*combiner.dest) == vector<string>{"GRCh38#0#chrM[0-7]"});
        REQUIRE(combiner.dest->get_is_circular(
                    combiner.dest->get_path_handle("GRCh38#0#chrM[0-7]")));
    }

    SECTION("three fragments merge in one pass") {
        TestGraphCombiner combiner(SeamPolicy::RENUMBER, false, false);
        combiner.dest = make_unique<bdsg::HashGraph>();
        handle_t n1 = combiner.dest->create_handle("AAAA");
        handle_t n2 = combiner.dest->create_handle("CC");
        handle_t n3 = combiner.dest->create_handle("GGG");
        add_ref_path(*combiner.dest, "GRCh38", "chr1", subrange_t{0, 4}, {n1});
        add_ref_path(*combiner.dest, "GRCh38", "chr1", subrange_t{4, 6}, {n2});
        add_ref_path(*combiner.dest, "GRCh38", "chr1", subrange_t{6, 9}, {n3});

        combiner.merge_path_fragments();

        REQUIRE(sorted_path_names(*combiner.dest) == vector<string>{"GRCh38#0#chr1[0-9]"});
        REQUIRE(spell_path(*combiner.dest, "GRCh38#0#chr1[0-9]") == "AAAACCGGG");
        REQUIRE(combiner.dest->has_edge(n1, n2));
        REQUIRE(combiner.dest->has_edge(n2, n3));
    }
}

}
}
