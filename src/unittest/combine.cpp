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
    TestGraphCombiner(SeamPolicy seam_policy, bool fuse, bool phase_block_is_offset)
        : GraphCombiner(seam_policy, fuse, phase_block_is_offset, "test") {
    }
    using GraphCombiner::Chunk;
    using GraphCombiner::SeamState;
    using GraphCombiner::BasePathKey;
    using GraphCombiner::GroupKey;
    using GraphCombiner::path_start_offset;
    using GraphCombiner::shift_ids_above;
    using GraphCombiner::group_key_of;
    using GraphCombiner::path_prefix;
    using GraphCombiner::path_suffix;
    using GraphCombiner::check_no_split_paths;
    using GraphCombiner::check_subrange_matches_length;
    using GraphCombiner::describe_reference;
    using GraphCombiner::sort_by_base_path_offset;
    using GraphCombiner::check_overlaps_agree;
    using GraphCombiner::check_all_paths_start_at_reference_first_node;
    using GraphCombiner::check_all_paths_end_at_reference_last_node;
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

/// Create a path with the given metadata, that visits the given steps.
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

/// Create a REFERENCE-sense path.
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

/// Build a chunk by hand on a locus of GRCh38, without a graph; sorting doesn't look at it.
static TestGraphCombiner::Chunk chunk_on_locus(const string& name, const string& locus,
                                               offset_t start) {
    TestGraphCombiner::Chunk chunk;
    chunk.name = name;
    chunk.ref.base_path_key = TestGraphCombiner::BasePathKey{"GRCh38", locus, 0};
    chunk.ref.start_offset = start;
    return chunk;
}

/// Same, for chunks that differ in sample or haplotype.
static TestGraphCombiner::Chunk chunk_on_base_path(const string& name,
                                                   const TestGraphCombiner::BasePathKey& base_path_key,
                                                   offset_t start) {
    TestGraphCombiner::Chunk chunk;
    chunk.name = name;
    chunk.ref.base_path_key = base_path_key;
    chunk.ref.start_offset = start;
    return chunk;
}

TEST_CASE("path_start_offset reads the offset from a subrange", "[combine]") {
    TestGraphCombiner offsets(SeamPolicy::RENUMBER, false, true);
    TestGraphCombiner not_offsets(SeamPolicy::RENUMBER, false, false);
    REQUIRE(offsets.path_start_offset(subrange_t{100, 200}, 7) == 100);
    REQUIRE(not_offsets.path_start_offset(subrange_t{100, 200}, 7) == 100);
}

TEST_CASE("path_start_offset prefers a subrange over a phase block when phase blocks are offsets", "[combine]") {
    TestGraphCombiner offsets(SeamPolicy::RENUMBER, false, true);
    REQUIRE(offsets.path_start_offset(subrange_t{100, PathMetadata::NO_END_POSITION},
                                      7) == 100);
}

TEST_CASE("path_start_offset reads the phase block when there is no subrange and phase blocks are offsets", "[combine]") {
    TestGraphCombiner offsets(SeamPolicy::RENUMBER, false, true);
    REQUIRE(offsets.path_start_offset(PathMetadata::NO_SUBRANGE, 7) == 7);
}

TEST_CASE("path_start_offset ignores the phase block when phase blocks are not offsets", "[combine]") {
    TestGraphCombiner not_offsets(SeamPolicy::RENUMBER, false, false);
    REQUIRE(not_offsets.path_start_offset(PathMetadata::NO_SUBRANGE, 7) == 0);
}

TEST_CASE("path_start_offset returns zero with no subrange and no phase block", "[combine]") {
    TestGraphCombiner offsets(SeamPolicy::RENUMBER, false, true);
    TestGraphCombiner not_offsets(SeamPolicy::RENUMBER, false, false);
    REQUIRE(offsets.path_start_offset(PathMetadata::NO_SUBRANGE,
                                      PathMetadata::NO_PHASE_BLOCK) == 0);
    REQUIRE(not_offsets.path_start_offset(PathMetadata::NO_SUBRANGE,
                                          PathMetadata::NO_PHASE_BLOCK) == 0);
}

TEST_CASE("shift_ids_above leaves IDs already above the given ID alone", "[combine]") {
    bdsg::HashGraph graph;
    graph.create_handle("A", 10);
    graph.create_handle("C", 12);
    REQUIRE(TestGraphCombiner::shift_ids_above(&graph, 3) == 0);
    REQUIRE(graph.min_node_id() == 10);
    REQUIRE(graph.max_node_id() == 12);
}

TEST_CASE("shift_ids_above shifts overlapping IDs above the given ID", "[combine]") {
    bdsg::HashGraph graph;
    graph.create_handle("A", 1);
    graph.create_handle("C", 5);
    int64_t shift = TestGraphCombiner::shift_ids_above(&graph, 10);
    REQUIRE(shift == 10);
    REQUIRE(graph.min_node_id() == 11);
    REQUIRE(graph.max_node_id() == 15);
}

TEST_CASE("shift_ids_above moves IDs that start exactly at the given ID above it", "[combine]") {
    // delta is 0 here, which is not negative, so the shift is 1: IDs have to end
    // up strictly above the mark, not on it.
    bdsg::HashGraph graph;
    graph.create_handle("A", 7);
    REQUIRE(TestGraphCombiner::shift_ids_above(&graph, 7) == 1);
    REQUIRE(graph.min_node_id() == 8);
}

TEST_CASE("shift_ids_above returns the shift it applied", "[combine]") {
    bdsg::HashGraph graph;
    graph.create_handle("A", 4);
    graph.create_handle("C", 9);
    int64_t shift = TestGraphCombiner::shift_ids_above(&graph, 20);
    REQUIRE(graph.get_id(graph.get_handle(4 + shift)) == 4 + shift);
    REQUIRE(graph.has_node(9 + shift));
    REQUIRE(!graph.has_node(4));
}

TEST_CASE("group_key_of groups pieces of one path together", "[combine]") {
    bdsg::HashGraph graph;
    handle_t n1 = graph.create_handle("AAAA");

    add_path(graph, PathSense::HAPLOTYPE, "HG002", "chr1", 1, 0,
             PathMetadata::NO_SUBRANGE, {n1});
    add_path(graph, PathSense::HAPLOTYPE, "HG002", "chr1", 1, 10,
             PathMetadata::NO_SUBRANGE, {n1});

    path_handle_t block_0 = graph.get_path_handle("HG002#1#chr1#0");
    path_handle_t block_10 = graph.get_path_handle("HG002#1#chr1#10");

    TestGraphCombiner offsets(SeamPolicy::RENUMBER, false, true);
    TestGraphCombiner not_offsets(SeamPolicy::RENUMBER, false, false);

    SECTION("phase blocks tell paths apart by default") {
        REQUIRE(not_offsets.group_key_of(graph, block_0)
                != not_offsets.group_key_of(graph, block_10));
    }

    SECTION("when phase blocks are offsets the phase block is left out of the key") {
        REQUIRE(offsets.group_key_of(graph, block_0)
                == offsets.group_key_of(graph, block_10));
        // Left out means replaced by NO_PHASE_BLOCK, not kept as some other value.
        REQUIRE(std::get<4>(offsets.group_key_of(graph, block_0))
                == PathMetadata::NO_PHASE_BLOCK);
    }

    SECTION("fragments with different subranges share a key") {
        add_path(graph, PathSense::REFERENCE, "GRCh38", "chr1", 0,
                 PathMetadata::NO_PHASE_BLOCK, subrange_t{0, 4}, {n1});
        add_path(graph, PathSense::REFERENCE, "GRCh38", "chr1", 0,
                 PathMetadata::NO_PHASE_BLOCK, subrange_t{4, 8}, {n1});
        REQUIRE(not_offsets.group_key_of(graph, graph.get_path_handle("GRCh38#0#chr1[0-4]"))
                == not_offsets.group_key_of(graph, graph.get_path_handle("GRCh38#0#chr1[4-8]")));
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
            return not_offsets.group_key_of(graph, graph.get_path_handle(name));
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

    SECTION("asking for zero bases or reading an empty path yields nothing") {
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

TEST_CASE("check_no_split_paths accepts a graph with a single path", "[combine]") {
    TestGraphCombiner combiner(SeamPolicy::RENUMBER, false, false);
    TestGraphCombiner::Chunk chunk = make_chunk("one.pg");
    handle_t n1 = chunk.graph->create_handle("AAAA");
    add_ref_path(*chunk.graph, "GRCh38", "chr1", subrange_t{0, 4}, {n1});
    combiner.check_no_split_paths(chunk);
    REQUIRE(chunk.graph->get_path_count() == 1);
}

TEST_CASE("check_no_split_paths accepts paths that differ only by haplotype", "[combine]") {
    TestGraphCombiner combiner(SeamPolicy::RENUMBER, false, false);
    TestGraphCombiner::Chunk chunk = make_chunk("haps.pg");
    handle_t n1 = chunk.graph->create_handle("AAAA");
    add_path(*chunk.graph, PathSense::HAPLOTYPE, "HG002", "chr1", 1, 0,
             subrange_t{0, 4}, {n1});
    add_path(*chunk.graph, PathSense::HAPLOTYPE, "HG002", "chr1", 2, 0,
             subrange_t{0, 4}, {n1});
    combiner.check_no_split_paths(chunk);
    REQUIRE(chunk.graph->get_path_count() == 2);
}

TEST_CASE("check_no_split_paths accepts paths that differ only by locus", "[combine]") {
    TestGraphCombiner combiner(SeamPolicy::RENUMBER, false, false);
    TestGraphCombiner::Chunk chunk = make_chunk("loci.pg");
    handle_t n1 = chunk.graph->create_handle("AAAA");
    add_ref_path(*chunk.graph, "GRCh38", "chr1", subrange_t{0, 4}, {n1});
    add_ref_path(*chunk.graph, "GRCh38", "chr2", subrange_t{0, 4}, {n1});
    combiner.check_no_split_paths(chunk);
    REQUIRE(chunk.graph->get_path_count() == 2);
}

TEST_CASE("check_no_split_paths accepts paths that differ only by phase block when phase blocks are not offsets", "[combine]") {
    TestGraphCombiner combiner(SeamPolicy::RENUMBER, false, false);
    TestGraphCombiner::Chunk chunk = make_chunk("blocks.pg");
    handle_t n1 = chunk.graph->create_handle("AAAA");
    add_path(*chunk.graph, PathSense::HAPLOTYPE, "HG002", "chr1", 1, 0,
             PathMetadata::NO_SUBRANGE, {n1});
    add_path(*chunk.graph, PathSense::HAPLOTYPE, "HG002", "chr1", 1, 10,
             PathMetadata::NO_SUBRANGE, {n1});
    combiner.check_no_split_paths(chunk);
    REQUIRE(chunk.graph->get_path_count() == 2);
}

TEST_CASE("check_no_split_paths rejects two pieces of one path in one graph", "[combine]") {
    TestGraphCombiner combiner(SeamPolicy::RENUMBER, false, false);
    TestGraphCombiner::Chunk chunk = make_chunk("split.pg");
    handle_t n1 = chunk.graph->create_handle("AAAA");
    handle_t n2 = chunk.graph->create_handle("CCC");
    add_ref_path(*chunk.graph, "GRCh38", "chr1", subrange_t{0, 4}, {n1});
    add_ref_path(*chunk.graph, "GRCh38", "chr1", subrange_t{4, 7}, {n2});
    REQUIRE_THROWS_AS(combiner.check_no_split_paths(chunk), InputConflictError);
    REQUIRE_THROWS_WITH(combiner.check_no_split_paths(chunk), Catch::Contains("split.pg"));
}

TEST_CASE("check_no_split_paths rejects two phase blocks of one path when phase blocks are offsets", "[combine]") {
    TestGraphCombiner combiner(SeamPolicy::RENUMBER, false, true);
    TestGraphCombiner::Chunk chunk = make_chunk("blocks.pg");
    handle_t n1 = chunk.graph->create_handle("AAAA");
    add_path(*chunk.graph, PathSense::HAPLOTYPE, "HG002", "chr1", 1, 0,
             PathMetadata::NO_SUBRANGE, {n1});
    add_path(*chunk.graph, PathSense::HAPLOTYPE, "HG002", "chr1", 1, 10,
             PathMetadata::NO_SUBRANGE, {n1});
    REQUIRE_THROWS_AS(combiner.check_no_split_paths(chunk), InputConflictError);
    REQUIRE_THROWS_WITH(combiner.check_no_split_paths(chunk), Catch::Contains("blocks.pg"));
}

TEST_CASE("check_subrange_matches_length rejects a subrange end that disagrees with the sequence length", "[combine]") {
    // A subrange claiming 100-999 over 11bp of sequence.
    TestGraphCombiner::Chunk chunk = make_chunk("claim.pg");
    handle_t n1 = chunk.graph->create_handle("AAAAAAAAAAA");
    add_ref_path(*chunk.graph, "GRCh38", "chr1", subrange_t{100, 999}, {n1});

    REQUIRE_THROWS_AS(TestGraphCombiner::check_subrange_matches_length(chunk),
                      PathLengthMismatchError);
    REQUIRE_THROWS_WITH(TestGraphCombiner::check_subrange_matches_length(chunk),
                        Catch::Contains("GRCh38#0#chr1[100-999]"));
}

TEST_CASE("describe_reference fills in every field from the reference path's metadata", "[combine]") {
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

    chunk.ref = combiner.describe_reference(chunk);

    REQUIRE(chunk.ref.path_name == "GRCh38#0#chr1[100-111]");
    REQUIRE(chunk.ref.base_path_key == TestGraphCombiner::BasePathKey{"GRCh38", "chr1", 0});
    REQUIRE(chunk.ref.start_offset == 100);
    REQUIRE(chunk.ref.end_offset == 111);
    REQUIRE(chunk.graph->get_id(chunk.ref.start_handle) == chunk.graph->get_id(n1));
    REQUIRE(chunk.graph->get_id(chunk.ref.end_handle) == chunk.graph->get_id(n3));
    REQUIRE(!chunk.graph->get_is_reverse(chunk.ref.start_handle));
    REQUIRE(!chunk.graph->get_is_reverse(chunk.ref.end_handle));
}

TEST_CASE("describe_reference gives a single-node chunk the same start and end node", "[combine]") {
    TestGraphCombiner combiner(SeamPolicy::COORD_TRIM, false, false);
    TestGraphCombiner::Chunk chunk = make_chunk("solo.pg");
    handle_t n1 = chunk.graph->create_handle("AAAA");
    add_ref_path(*chunk.graph, "GRCh38", "chr1", subrange_t{0, 4}, {n1});

    chunk.ref = combiner.describe_reference(chunk);
    REQUIRE(chunk.ref.start_handle == chunk.ref.end_handle);
    REQUIRE(chunk.ref.start_offset == 0);
    REQUIRE(chunk.ref.end_offset == 4);
}

TEST_CASE("describe_reference starts a reference path with no subrange at zero whether or not phase blocks are offsets", "[combine]") {
    // A REFERENCE-sense path may not carry a phase block, so with no subrange there
    // is no offset to read and the chunk can only be placed at offset 0.
    for (bool phase_block_is_offset : {false, true}) {
        TestGraphCombiner combiner(SeamPolicy::COORD_TRIM, false, phase_block_is_offset);
        TestGraphCombiner::Chunk chunk = make_chunk("bare.pg");
        handle_t n1 = chunk.graph->create_handle("AAAA");
        add_ref_path(*chunk.graph, "GRCh38", "chr1", PathMetadata::NO_SUBRANGE, {n1});

        chunk.ref = combiner.describe_reference(chunk);
        REQUIRE(chunk.ref.start_offset == 0);
        REQUIRE(chunk.ref.end_offset == 4);
    }
}

TEST_CASE("describe_reference leaves circularity out of the base path key", "[combine]") {
    TestGraphCombiner combiner(SeamPolicy::COORD_TRIM, false, false);
    TestGraphCombiner::Chunk chunk = make_chunk("circle.pg");
    handle_t n1 = chunk.graph->create_handle("AAAA");
    add_path(*chunk.graph, PathSense::REFERENCE, "GRCh38", "chrM", 0,
             PathMetadata::NO_PHASE_BLOCK, subrange_t{0, 4}, {n1}, true);

    chunk.ref = combiner.describe_reference(chunk);
    REQUIRE(chunk.ref.base_path_key == TestGraphCombiner::BasePathKey{"GRCh38", "chrM", 0});
    REQUIRE(chunk.ref.is_circular);
}

TEST_CASE("describe_reference rejects an input with no reference path", "[combine]") {
    TestGraphCombiner combiner(SeamPolicy::COORD_TRIM, false, false);
    TestGraphCombiner::Chunk chunk = make_chunk("noref.pg");
    add_path(*chunk.graph, PathSense::GENERIC, PathMetadata::NO_SAMPLE_NAME, "plain_contig",
             PathMetadata::NO_HAPLOTYPE, PathMetadata::NO_PHASE_BLOCK,
             PathMetadata::NO_SUBRANGE, {chunk.graph->create_handle("AAAA")});
    REQUIRE_THROWS_AS(combiner.describe_reference(chunk), ReferencePathError);
    REQUIRE_THROWS_WITH(combiner.describe_reference(chunk), Catch::Contains("noref.pg"));
}

TEST_CASE("describe_reference rejects an input with two reference paths", "[combine]") {
    TestGraphCombiner combiner(SeamPolicy::COORD_TRIM, false, false);
    TestGraphCombiner::Chunk chunk = make_chunk("tworef.pg");
    handle_t n1 = chunk.graph->create_handle("AAAA");
    add_ref_path(*chunk.graph, "GRCh38", "chr1", subrange_t{0, 4}, {n1});
    add_ref_path(*chunk.graph, "GRCh38", "chr2", subrange_t{0, 4}, {n1});
    REQUIRE_THROWS_AS(combiner.describe_reference(chunk), ReferencePathError);
    REQUIRE_THROWS_WITH(combiner.describe_reference(chunk), Catch::Contains("tworef.pg"));
}

TEST_CASE("describe_reference rejects a reference path that visits its first node in reverse", "[combine]") {
    TestGraphCombiner combiner(SeamPolicy::COORD_TRIM, false, false);
    TestGraphCombiner::Chunk chunk = make_chunk("rev.pg");
    handle_t n1 = chunk.graph->create_handle("AAAA");
    handle_t n2 = chunk.graph->create_handle("CCC");
    chunk.graph->create_edge(chunk.graph->flip(n1), n2);
    add_ref_path(*chunk.graph, "GRCh38", "chr1", subrange_t{0, 7}, {chunk.graph->flip(n1), n2});
    REQUIRE_THROWS_AS(combiner.describe_reference(chunk), ReferencePathError);
    REQUIRE_THROWS_WITH(combiner.describe_reference(chunk), Catch::Contains("GRCh38#0#chr1[0-7]"));
}

TEST_CASE("sort_by_base_path_offset sorts chunks of one base path by offset", "[combine]") {
    TestGraphCombiner combiner(SeamPolicy::COORD_TRIM, false, false);
    vector<TestGraphCombiner::Chunk> chunks;
    chunks.push_back(chunk_on_locus("c", "chr1", 200));
    chunks.push_back(chunk_on_locus("a", "chr1", 0));
    chunks.push_back(chunk_on_locus("b", "chr1", 100));

    combiner.sort_by_base_path_offset(chunks);
    REQUIRE(chunks[0].name == "a");
    REQUIRE(chunks[1].name == "b");
    REQUIRE(chunks[2].name == "c");
}

TEST_CASE("sort_by_base_path_offset keeps base paths in first-seen order with their chunks together", "[combine]") {
    TestGraphCombiner combiner(SeamPolicy::COORD_TRIM, false, false);
    // chr2 is seen first, so it stays first even though chr1 sorts earlier.
    vector<TestGraphCombiner::Chunk> chunks;
    chunks.push_back(chunk_on_locus("y2", "chr2", 100));
    chunks.push_back(chunk_on_locus("x2", "chr1", 100));
    chunks.push_back(chunk_on_locus("y1", "chr2", 0));
    chunks.push_back(chunk_on_locus("x1", "chr1", 0));

    combiner.sort_by_base_path_offset(chunks);
    REQUIRE(chunks[0].name == "y1");
    REQUIRE(chunks[1].name == "y2");
    REQUIRE(chunks[2].name == "x1");
    REQUIRE(chunks[3].name == "x2");
}

TEST_CASE("sort_by_base_path_offset keeps chunks with equal offsets in input order", "[combine]") {
    TestGraphCombiner combiner(SeamPolicy::COORD_TRIM, false, false);
    vector<TestGraphCombiner::Chunk> chunks;
    chunks.push_back(chunk_on_locus("first", "chr1", 50));
    chunks.push_back(chunk_on_locus("second", "chr1", 50));
    chunks.push_back(chunk_on_locus("third", "chr1", 50));

    combiner.sort_by_base_path_offset(chunks);
    REQUIRE(chunks[0].name == "first");
    REQUIRE(chunks[1].name == "second");
    REQUIRE(chunks[2].name == "third");
}

TEST_CASE("sort_by_base_path_offset treats haplotypes of one locus as separate base paths", "[combine]") {
    TestGraphCombiner combiner(SeamPolicy::COORD_TRIM, false, false);
    // All three share sample and locus, so only haplotype can tell the two base
    // paths apart.
    vector<TestGraphCombiner::Chunk> chunks;
    chunks.push_back(chunk_on_base_path("hap2", {"GRCh38", "chr1", 2}, 50));
    chunks.push_back(chunk_on_locus("hap0_late", "chr1", 100));
    chunks.push_back(chunk_on_locus("hap0_early", "chr1", 0));

    combiner.sort_by_base_path_offset(chunks);
    // Base paths lead in first-seen order; hap0's two chunks follow in offset order.
    REQUIRE(chunks[0].name == "hap2");
    REQUIRE(chunks[1].name == "hap0_early");
    REQUIRE(chunks[2].name == "hap0_late");
}

TEST_CASE("check_overlaps_agree rejects neighbouring chunks that spell different sequences over their overlap", "[combine]") {
    TestGraphCombiner combiner(SeamPolicy::COORD_TRIM, false, false);
    vector<TestGraphCombiner::Chunk> chunks;
    chunks.push_back(make_chunk("left.pg"));
    add_ref_path(*chunks[0].graph, "GRCh38", "chr1", subrange_t{0, 10},
                 {chunks[0].graph->create_handle("AAAAAAAAAA")});
    chunks[0].ref = combiner.describe_reference(chunks[0]);
    chunks.push_back(make_chunk("right.pg"));
    handle_t r1 = chunks[1].graph->create_handle("AACCC");
    handle_t r2 = chunks[1].graph->create_handle("GGGGG");
    chunks[1].graph->create_edge(r1, r2);
    add_ref_path(*chunks[1].graph, "GRCh38", "chr1", subrange_t{6, 16}, {r1, r2});
    chunks[1].ref = combiner.describe_reference(chunks[1]);

    // The overlap is [6, 10): the left chunk has AAAA there and the right chunk AACC.
    REQUIRE_THROWS_AS(combiner.check_overlaps_agree(chunks), ChunkOverlapError);
    REQUIRE_THROWS_WITH(combiner.check_overlaps_agree(chunks), Catch::Contains("at offset 8"));
}

TEST_CASE("check_all_paths_start_at_reference_first_node and check_all_paths_end_at_reference_last_node pass chunks whose paths all start and end there", "[combine]") {
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
    chunk.ref.start_handle = n1;
    chunk.ref.end_handle = n4;

    SECTION("all paths beginning on the boundary node are accepted") {
        combiner.check_all_paths_start_at_reference_first_node(chunk);
        REQUIRE(chunk.graph->get_path_count() == 3);
    }

    SECTION("all paths ending on the boundary node are accepted") {
        combiner.check_all_paths_end_at_reference_last_node(chunk);
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
        flipped.ref.start_handle = flipped.graph->flip(a);
        flipped.ref.end_handle = b;

        combiner.check_all_paths_start_at_reference_first_node(flipped);
        combiner.check_all_paths_end_at_reference_last_node(flipped);
        REQUIRE(flipped.graph->get_path_count() == 1);
    }
}

TEST_CASE("check_all_paths_start_at_reference_first_node rejects a path running against the reference path", "[combine]") {
    TestGraphCombiner combiner(SeamPolicy::COORD_TRIM, false, false);
    TestGraphCombiner::Chunk chunk = make_chunk("reverse.pg");
    handle_t a = chunk.graph->create_handle("AAAA");
    handle_t b = chunk.graph->create_handle("CCC");
    chunk.graph->create_edge(a, b);
    add_ref_path(*chunk.graph, "GRCh38", "chr1", subrange_t{0, 7}, {a, b});
    add_path(*chunk.graph, PathSense::HAPLOTYPE, "HG002", "chr1", 1, 0,
             subrange_t{0, 7}, {chunk.graph->flip(b), chunk.graph->flip(a)});
    chunk.ref.start_handle = a;
    chunk.ref.end_handle = b;
    REQUIRE_THROWS_AS(combiner.check_all_paths_start_at_reference_first_node(chunk), ChunkOverlapError);
}

TEST_CASE("check_all_paths_start_at_reference_first_node rejects a path that starts on a different node", "[combine]") {
    TestGraphCombiner combiner(SeamPolicy::COORD_TRIM, false, false);
    TestGraphCombiner::Chunk chunk = make_chunk("other_start.pg");
    handle_t n1 = chunk.graph->create_handle("AAAA");
    handle_t n2 = chunk.graph->create_handle("CCC");
    handle_t n3 = chunk.graph->create_handle("T");
    chunk.graph->create_edge(n1, n2);
    chunk.graph->create_edge(n3, n2);
    add_ref_path(*chunk.graph, "GRCh38", "chr1", subrange_t{0, 7}, {n1, n2});
    add_path(*chunk.graph, PathSense::HAPLOTYPE, "HG002", "chr1", 1, 0,
             subrange_t{3, 7}, {n3, n2});
    chunk.ref.start_handle = n1;
    chunk.ref.end_handle = n2;
    REQUIRE_THROWS_AS(combiner.check_all_paths_start_at_reference_first_node(chunk), ChunkOverlapError);
    REQUIRE_THROWS_WITH(combiner.check_all_paths_start_at_reference_first_node(chunk),
                        Catch::Contains("HG002#1#chr1#0[3-7]"));
}

TEST_CASE("check_all_paths_end_at_reference_last_node rejects a path that ends on a different node", "[combine]") {
    TestGraphCombiner combiner(SeamPolicy::COORD_TRIM, false, false);
    TestGraphCombiner::Chunk chunk = make_chunk("other_end.pg");
    handle_t n1 = chunk.graph->create_handle("AAAA");
    handle_t n2 = chunk.graph->create_handle("CCC");
    handle_t n3 = chunk.graph->create_handle("T");
    chunk.graph->create_edge(n1, n2);
    chunk.graph->create_edge(n1, n3);
    add_ref_path(*chunk.graph, "GRCh38", "chr1", subrange_t{0, 7}, {n1, n2});
    add_path(*chunk.graph, PathSense::HAPLOTYPE, "HG002", "chr1", 1, 0,
             subrange_t{0, 5}, {n1, n3});
    chunk.ref.start_handle = n1;
    chunk.ref.end_handle = n2;
    REQUIRE_THROWS_AS(combiner.check_all_paths_end_at_reference_last_node(chunk), ChunkOverlapError);
    REQUIRE_THROWS_WITH(combiner.check_all_paths_end_at_reference_last_node(chunk),
                        Catch::Contains("HG002#1#chr1#0[0-5]"));
}

TEST_CASE("trim_and_shift with no overlap and no fusing leaves the chunk unchanged", "[combine]") {
    TestGraphCombiner combiner(SeamPolicy::COORD_TRIM, false, false);
    TestGraphCombiner::Chunk chunk = make_chunk("clean.pg");
    handle_t n1 = chunk.graph->create_handle("AAAA");
    handle_t n2 = chunk.graph->create_handle("CCC");
    chunk.graph->create_edge(n1, n2);
    add_ref_path(*chunk.graph, "GRCh38", "chr1", subrange_t{100, 107}, {n1, n2});
    chunk.ref.start_handle = n1;
    nid_t before = chunk.graph->get_id(n1);

    nid_t after = combiner.trim_and_shift(chunk, 0, 5);

    REQUIRE(after == before);
    REQUIRE(chunk.graph->get_node_count() == 2);
    REQUIRE(sorted_path_names(*chunk.graph) == vector<string>{"GRCh38#0#chr1[100-107]"});
}

TEST_CASE("trim_and_shift with fusing and no overlap moves coordinates back by the left boundary node length", "[combine]") {
    // Fusing merges the left chunk's last node into this chunk's first, so
    // this chunk's coordinates move back by that node's length.
    TestGraphCombiner combiner(SeamPolicy::COORD_TRIM, true, false);
    TestGraphCombiner::Chunk chunk = make_chunk("fuse.pg");
    handle_t n1 = chunk.graph->create_handle("AAAA");
    handle_t n2 = chunk.graph->create_handle("CCC");
    chunk.graph->create_edge(n1, n2);
    add_ref_path(*chunk.graph, "GRCh38", "chr1", subrange_t{100, 107}, {n1, n2});
    chunk.ref.start_handle = n1;

    combiner.trim_and_shift(chunk, 0, 5);

    REQUIRE(chunk.graph->get_node_count() == 2);
    REQUIRE(sorted_path_names(*chunk.graph) == vector<string>{"GRCh38#0#chr1[95-107]"});
}

TEST_CASE("trim_and_shift splits the first node at the overlap and shifts coordinates", "[combine]") {
    TestGraphCombiner combiner(SeamPolicy::COORD_TRIM, false, false);
    TestGraphCombiner::Chunk chunk = make_chunk("overlap.pg");
    handle_t n1 = chunk.graph->create_handle("AAAACC");
    handle_t n2 = chunk.graph->create_handle("GGG");
    chunk.graph->create_edge(n1, n2);
    add_ref_path(*chunk.graph, "GRCh38", "chr1", subrange_t{100, 109}, {n1, n2});
    add_path(*chunk.graph, PathSense::HAPLOTYPE, "HG002", "chr1", 1, 0,
             subrange_t{100, 109}, {n1, n2});
    chunk.ref.start_handle = n1;

    nid_t new_start = combiner.trim_and_shift(chunk, 4, 5);

    // The overlapping "AAAA" is gone; "CC" remains and is the new start.
    REQUIRE(chunk.graph->get_node_count() == 2);
    REQUIRE(chunk.graph->get_sequence(chunk.graph->get_handle(new_start)) == "CC");
    REQUIRE(spell_path(*chunk.graph, "GRCh38#0#chr1[104-109]") == "CCGGG");
    REQUIRE(spell_path(*chunk.graph, "HG002#1#chr1#0[104-109]") == "CCGGG");
}

TEST_CASE("trim_and_shift with a one-base overlap removes only the first base of the node", "[combine]") {
    TestGraphCombiner combiner(SeamPolicy::COORD_TRIM, false, false);
    TestGraphCombiner::Chunk chunk = make_chunk("thin.pg");
    handle_t n1 = chunk.graph->create_handle("AC");
    add_ref_path(*chunk.graph, "GRCh38", "chr1", subrange_t{100, 102}, {n1});
    chunk.ref.start_handle = n1;

    nid_t new_start = combiner.trim_and_shift(chunk, 1, 5);

    REQUIRE(chunk.graph->get_node_count() == 1);
    REQUIRE(chunk.graph->get_sequence(chunk.graph->get_handle(new_start)) == "C");
    REQUIRE(spell_path(*chunk.graph, "GRCh38#0#chr1[101-102]") == "C");
}

TEST_CASE("trim_and_shift shifts the phase block instead of a subrange when phase blocks are offsets", "[combine]") {
    TestGraphCombiner combiner(SeamPolicy::COORD_TRIM, false, true);
    TestGraphCombiner::Chunk chunk = make_chunk("pblock.pg");
    handle_t n1 = chunk.graph->create_handle("AAAACC");
    add_path(*chunk.graph, PathSense::HAPLOTYPE, "HG002", "chr1", 1, 100,
             PathMetadata::NO_SUBRANGE, {n1});
    chunk.ref.start_handle = n1;

    combiner.trim_and_shift(chunk, 4, 5);

    REQUIRE(sorted_path_names(*chunk.graph) == vector<string>{"HG002#1#chr1#104"});
}

TEST_CASE("trim_and_shift adds a subrange starting at the trimmed length to a path with no subrange", "[combine]") {
    TestGraphCombiner combiner(SeamPolicy::COORD_TRIM, false, false);
    TestGraphCombiner::Chunk chunk = make_chunk("bare.pg");
    handle_t n1 = chunk.graph->create_handle("AAAACC");
    add_path(*chunk.graph, PathSense::GENERIC, PathMetadata::NO_SAMPLE_NAME, "scaffold",
             PathMetadata::NO_HAPLOTYPE, PathMetadata::NO_PHASE_BLOCK,
             PathMetadata::NO_SUBRANGE, {n1});
    chunk.ref.start_handle = n1;

    combiner.trim_and_shift(chunk, 4, 0);

    REQUIRE(sorted_path_names(*chunk.graph) == vector<string>{"scaffold[4]"});
}

TEST_CASE("trim_and_shift rejects a shift that would make a start coordinate negative", "[combine]") {
    // Fusing moves coordinates back by 50, but the path only starts at 10.
    TestGraphCombiner combiner(SeamPolicy::COORD_TRIM, true, false);
    TestGraphCombiner::Chunk chunk = make_chunk("negative.pg");
    handle_t n1 = chunk.graph->create_handle("AAAA");
    add_ref_path(*chunk.graph, "GRCh38", "chr1", subrange_t{10, 14}, {n1});
    chunk.ref.start_handle = n1;

    REQUIRE_THROWS_AS(combiner.trim_and_shift(chunk, 0, 50), ChunkOverlapError);
    REQUIRE_THROWS_WITH(combiner.trim_and_shift(chunk, 0, 50),
                        Catch::Contains("GRCh38#0#chr1[10-14]"));
}

TEST_CASE("trim_and_shift keeps each path on its own alleles through a bubble", "[combine]") {
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
    chunk.ref.start_handle = n1;

    combiner.trim_and_shift(chunk, 4, 0);

    REQUIRE(spell_path(*chunk.graph, "GRCh38#0#chr1[104-112]") == "CCGTTTTT");
    REQUIRE(spell_path(*chunk.graph, "HG002#1#chr1#0[104-112]") == "CCTTTTTT");
}

TEST_CASE("trim_and_shift rejects an overlap longer than the first node", "[combine]") {
    TestGraphCombiner combiner(SeamPolicy::COORD_TRIM, false, false);
    TestGraphCombiner::Chunk chunk = make_chunk("short.pg");
    handle_t n1 = chunk.graph->create_handle("AA");
    add_ref_path(*chunk.graph, "GRCh38", "chr1", subrange_t{0, 2}, {n1});
    chunk.ref.start_handle = n1;
    REQUIRE_THROWS_AS(combiner.trim_and_shift(chunk, 3, 0), ChunkOverlapError);
    REQUIRE_THROWS_WITH(combiner.trim_and_shift(chunk, 3, 0), Catch::Contains("short.pg"));
}

TEST_CASE("trim_and_shift rejects an overlap equal to the first node", "[combine]") {
    TestGraphCombiner combiner(SeamPolicy::COORD_TRIM, false, false);
    TestGraphCombiner::Chunk chunk = make_chunk("exact.pg");
    handle_t n1 = chunk.graph->create_handle("AA");
    add_ref_path(*chunk.graph, "GRCh38", "chr1", subrange_t{0, 2}, {n1});
    chunk.ref.start_handle = n1;
    REQUIRE_THROWS_AS(combiner.trim_and_shift(chunk, 2, 0), ChunkOverlapError);
    REQUIRE_THROWS_WITH(combiner.trim_and_shift(chunk, 2, 0), Catch::Contains("exact.pg"));
}

TEST_CASE("trim_and_shift rejects an edge into the part of the first node it would cut off", "[combine]") {
    TestGraphCombiner combiner(SeamPolicy::COORD_TRIM, false, false);
    TestGraphCombiner::Chunk chunk = make_chunk("dangle.pg");
    handle_t other = chunk.graph->create_handle("TT");
    handle_t n1 = chunk.graph->create_handle("AAAACC");
    chunk.graph->create_edge(other, n1);
    add_ref_path(*chunk.graph, "GRCh38", "chr1", subrange_t{0, 6}, {n1});
    chunk.ref.start_handle = n1;
    REQUIRE_THROWS_AS(combiner.trim_and_shift(chunk, 4, 0), ChunkOverlapError);
}

TEST_CASE("ingest makes the first input dest unchanged", "[combine]") {
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

TEST_CASE("ingest shifts a second input's IDs above the first when renumbering", "[combine]") {
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

TEST_CASE("ingest does not shift node IDs when nodes are shared by ID", "[combine]") {
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

TEST_CASE("share_nodes_by_id reuses an ID already in dest and creates a new one", "[combine]") {
    TestGraphCombiner combiner(SeamPolicy::SHARED_IDS, false, false);
    combiner.dest = make_unique<bdsg::HashGraph>();
    handle_t d1 = combiner.dest->create_handle("AAAA", 1);
    handle_t d2 = combiner.dest->create_handle("CCC", 2);
    combiner.dest->create_edge(d1, d2);

    TestGraphCombiner::Chunk chunk = make_chunk("more.pg");
    chunk.graph->create_handle("CCC", 2);      // already in dest
    chunk.graph->create_handle("GG", 3);       // new

    combiner.share_nodes_by_id(chunk);

    REQUIRE(combiner.dest->get_node_count() == 3);
    REQUIRE(combiner.dest->get_sequence(combiner.dest->get_handle(2)) == "CCC");
    REQUIRE(combiner.dest->get_sequence(combiner.dest->get_handle(3)) == "GG");
}

TEST_CASE("share_nodes_by_id copies nodes with no edges into the resulting graph", "[combine]") {
    TestGraphCombiner combiner(SeamPolicy::SHARED_IDS, false, false);
    combiner.dest = make_unique<bdsg::HashGraph>();
    TestGraphCombiner::Chunk chunk = make_chunk("island.pg");
    chunk.graph->create_handle("TTTT", 7);

    combiner.share_nodes_by_id(chunk);

    REQUIRE(combiner.dest->has_node(7));
    REQUIRE(combiner.dest->get_sequence(combiner.dest->get_handle(7)) == "TTTT");
}

TEST_CASE("share_nodes_by_id keeps the orientation of every kind of edge", "[combine]") {
    for (bool from_rev : {false, true}) {
        for (bool to_rev : {false, true}) {
            TestGraphCombiner combiner(SeamPolicy::SHARED_IDS, false, false);
            combiner.dest = make_unique<bdsg::HashGraph>();
            TestGraphCombiner::Chunk chunk = make_chunk("edge.pg");
            chunk.graph->create_handle("AAAA", 1);
            chunk.graph->create_handle("CCC", 2);
            chunk.graph->create_edge(chunk.graph->get_handle(1, from_rev),
                                     chunk.graph->get_handle(2, to_rev));

            combiner.share_nodes_by_id(chunk);

            REQUIRE(combiner.dest->get_edge_count() == 1);
            REQUIRE(combiner.dest->has_edge(combiner.dest->get_handle(1, from_rev),
                                            combiner.dest->get_handle(2, to_rev)));
        }
    }
}

TEST_CASE("share_nodes_by_id rejects a repeated ID with a different sequence", "[combine]") {
    TestGraphCombiner combiner(SeamPolicy::SHARED_IDS, false, false);
    combiner.dest = make_unique<bdsg::HashGraph>();
    combiner.dest->create_handle("CCC", 2);
    TestGraphCombiner::Chunk chunk = make_chunk("clash.pg");
    chunk.graph->create_handle("GGG", 2);
    REQUIRE_THROWS_AS(combiner.share_nodes_by_id(chunk), InputConflictError);
    REQUIRE_THROWS_WITH(combiner.share_nodes_by_id(chunk), Catch::Contains("Node 2"));
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

    SECTION("every path is copied under its own name") {
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

    SECTION("the steps are copied in order onto dest's nodes") {
        REQUIRE(spell_path(*combiner.dest, "GRCh38#0#chr1[0-7]") == "AAAACCC");
        REQUIRE(combiner.dest->get_step_count(
                    combiner.dest->get_path_handle("GRCh38#0#chr1[0-7]")) == 2);
    }
}

TEST_CASE("copy_paths_checked rejects a path name that is already in dest", "[combine]") {
    TestGraphCombiner combiner(SeamPolicy::RENUMBER, false, false);
    combiner.dest = make_unique<bdsg::HashGraph>();
    add_ref_path(*combiner.dest, "GRCh38", "chr1", subrange_t{0, 4},
                 {combiner.dest->create_handle("AAAA")});
    TestGraphCombiner::Chunk chunk = make_chunk("dup.pg");
    add_ref_path(*chunk.graph, "GRCh38", "chr1", subrange_t{0, 4},
                 {chunk.graph->create_handle("AAAA")});
    REQUIRE_THROWS_AS(combiner.copy_paths_checked(chunk), InputConflictError);
    REQUIRE_THROWS_WITH(combiner.copy_paths_checked(chunk), Catch::Contains("GRCh38#0#chr1[0-4]"));
}

TEST_CASE("connect_seam without fusing joins the chunk to the base path with an edge", "[combine]") {
    TestGraphCombiner combiner(SeamPolicy::COORD_TRIM, false, false);
    combiner.dest = make_unique<bdsg::HashGraph>();
    handle_t left = combiner.dest->create_handle("AAAA", 1);
    handle_t start = combiner.dest->create_handle("CCC", 2);
    handle_t end = combiner.dest->create_handle("GG", 3);
    combiner.dest->create_edge(start, end);

    TestGraphCombiner::SeamState base_path_end{combiner.dest->get_id(left), false, 4, "left.pg"};
    handle_t new_end = combiner.connect_seam(base_path_end, start, end, "right.pg");

    REQUIRE(combiner.dest->has_edge(left, start));
    REQUIRE(combiner.dest->get_node_count() == 3);
    // Nothing was fused, so the chunk's own end is still the reference's end.
    REQUIRE(new_end == end);
}

TEST_CASE("connect_seam with fusing merges the two boundary nodes into one", "[combine]") {
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

    TestGraphCombiner::SeamState base_path_end{combiner.dest->get_id(left), false, 6, "left.pg"};
    handle_t new_end = combiner.connect_seam(base_path_end, start, end, "right.pg");

    REQUIRE(combiner.dest->get_node_count() == 3);
    REQUIRE(spell_path(*combiner.dest, "GRCh38#0#chr1[0-6]") == "TTAAAACCC");
    REQUIRE(new_end == end);
    REQUIRE(combiner.max_node_id >= combiner.dest->max_node_id());
}

TEST_CASE("connect_seam fusing a single-node chunk returns the handle to the fused node", "[combine]") {
    // When the chunk's reference ends on the very node being fused, the old
    // handle is destroyed by the fuse; the caller must get the new one back.
    TestGraphCombiner combiner(SeamPolicy::COORD_TRIM, true, false);
    combiner.dest = make_unique<bdsg::HashGraph>();
    handle_t left = combiner.dest->create_handle("AAAA", 1);
    handle_t start = combiner.dest->create_handle("CCC", 2);
    combiner.max_node_id = 2;

    TestGraphCombiner::SeamState base_path_end{combiner.dest->get_id(left), false, 4, "left.pg"};
    handle_t new_end = combiner.connect_seam(base_path_end, start, start, "right.pg");

    REQUIRE(combiner.dest->get_node_count() == 1);
    REQUIRE(combiner.dest->has_node(combiner.dest->get_id(new_end)));
    REQUIRE(combiner.dest->get_sequence(new_end) == "AAAACCC");
    REQUIRE(combiner.max_node_id >= combiner.dest->max_node_id());
}

TEST_CASE("connect_seam rejects fusing a boundary node that has an edge on the side being fused", "[combine]") {
    TestGraphCombiner combiner(SeamPolicy::COORD_TRIM, true, false);
    combiner.dest = make_unique<bdsg::HashGraph>();
    handle_t left = combiner.dest->create_handle("AAAA", 1);
    handle_t start = combiner.dest->create_handle("CCC", 2);
    handle_t other = combiner.dest->create_handle("GG", 3);
    combiner.dest->create_edge(left, other);
    combiner.max_node_id = 3;

    TestGraphCombiner::SeamState base_path_end{combiner.dest->get_id(left), false, 4, "left.pg"};
    REQUIRE_THROWS_AS(combiner.connect_seam(base_path_end, start, start, "right.pg"), FuseError);
    REQUIRE_THROWS_WITH(combiner.connect_seam(base_path_end, start, start, "right.pg"),
                        Catch::Contains("left.pg"));
}

TEST_CASE("merge_path_fragments merges two abutting fragments into one path spanning both", "[combine]") {
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

TEST_CASE("merge_path_fragments keeps the name of a lone fragment", "[combine]") {
    TestGraphCombiner combiner(SeamPolicy::RENUMBER, false, false);
    combiner.dest = make_unique<bdsg::HashGraph>();
    handle_t n1 = combiner.dest->create_handle("AAAA");
    add_ref_path(*combiner.dest, "GRCh38", "chr1", subrange_t{0, 4}, {n1});

    combiner.merge_path_fragments();

    REQUIRE(sorted_path_names(*combiner.dest) == vector<string>{"GRCh38#0#chr1[0-4]"});
}

TEST_CASE("merge_path_fragments keeps the name of a lone path with no subrange", "[combine]") {
    TestGraphCombiner combiner(SeamPolicy::RENUMBER, false, false);
    combiner.dest = make_unique<bdsg::HashGraph>();
    handle_t n1 = combiner.dest->create_handle("AAAACCCC");
    add_ref_path(*combiner.dest, "GRCh38", "chr1", PathMetadata::NO_SUBRANGE, {n1});

    combiner.merge_path_fragments();

    REQUIRE(sorted_path_names(*combiner.dest) == vector<string>{"GRCh38#0#chr1"});
}

TEST_CASE("merge_path_fragments does not repeat steps that two fragments share", "[combine]") {
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

TEST_CASE("merge_path_fragments does not merge fragments of different paths", "[combine]") {
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

TEST_CASE("merge_path_fragments merges three fragments given out of order in coordinate order", "[combine]") {
    TestGraphCombiner combiner(SeamPolicy::RENUMBER, false, false);
    combiner.dest = make_unique<bdsg::HashGraph>();
    handle_t n3 = combiner.dest->create_handle("GGG");
    handle_t n1 = combiner.dest->create_handle("AAAA");
    handle_t n2 = combiner.dest->create_handle("CC");
    add_ref_path(*combiner.dest, "GRCh38", "chr1", subrange_t{6, 9}, {n3});
    add_ref_path(*combiner.dest, "GRCh38", "chr1", subrange_t{0, 4}, {n1});
    add_ref_path(*combiner.dest, "GRCh38", "chr1", subrange_t{4, 6}, {n2});

    combiner.merge_path_fragments();

    REQUIRE(sorted_path_names(*combiner.dest) == vector<string>{"GRCh38#0#chr1[0-9]"});
    REQUIRE(spell_path(*combiner.dest, "GRCh38#0#chr1[0-9]") == "AAAACCGGG");
    REQUIRE(combiner.dest->has_edge(n1, n2));
    REQUIRE(combiner.dest->has_edge(n2, n3));
}

TEST_CASE("merge_path_fragments keeps the phase block of a merged haplotype path", "[combine]") {
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

TEST_CASE("merge_path_fragments positions fragments by phase block when phase blocks are offsets", "[combine]") {
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

TEST_CASE("merge_path_fragments keeps a circular path circular", "[combine]") {
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

TEST_CASE("merge_path_fragments rejects fragments that disagree on circularity", "[combine]") {
    TestGraphCombiner combiner(SeamPolicy::RENUMBER, false, false);
    combiner.dest = make_unique<bdsg::HashGraph>();
    handle_t n1 = combiner.dest->create_handle("AAAA");
    handle_t n2 = combiner.dest->create_handle("CCC");
    add_path(*combiner.dest, PathSense::REFERENCE, "GRCh38", "chrM", 0,
             PathMetadata::NO_PHASE_BLOCK, subrange_t{0, 4}, {n1}, true);
    add_path(*combiner.dest, PathSense::REFERENCE, "GRCh38", "chrM", 0,
             PathMetadata::NO_PHASE_BLOCK, subrange_t{4, 7}, {n2}, false);
    REQUIRE_THROWS_AS(combiner.merge_path_fragments(), InputConflictError);
    REQUIRE_THROWS_WITH(combiner.merge_path_fragments(), Catch::Contains("chrM"));
}

TEST_CASE("merge_path_fragments rejects a gap between fragments", "[combine]") {
    TestGraphCombiner combiner(SeamPolicy::RENUMBER, false, false);
    combiner.dest = make_unique<bdsg::HashGraph>();
    add_ref_path(*combiner.dest, "GRCh38", "chr1", subrange_t{0, 4},
                 {combiner.dest->create_handle("AAAA")});
    add_ref_path(*combiner.dest, "GRCh38", "chr1", subrange_t{6, 9},
                 {combiner.dest->create_handle("CCC")});
    REQUIRE_THROWS_AS(combiner.merge_path_fragments(), FragmentMergeError);
    REQUIRE_THROWS_WITH(combiner.merge_path_fragments(), Catch::Contains("Gap of 2 bp"));
}

TEST_CASE("merge_path_fragments rejects phase blocks that do not line up as offsets", "[combine]") {
    // Blocks 0 and 2 overlap, so this is not a gap, but 0 + 4bp should be followed by block 4.
    TestGraphCombiner combiner(SeamPolicy::RENUMBER, false, true);
    combiner.dest = make_unique<bdsg::HashGraph>();
    add_path(*combiner.dest, PathSense::HAPLOTYPE, "HG002", "chr1", 1, 0,
             PathMetadata::NO_SUBRANGE, {combiner.dest->create_handle("AAAA")});
    add_path(*combiner.dest, PathSense::HAPLOTYPE, "HG002", "chr1", 1, 2,
             PathMetadata::NO_SUBRANGE, {combiner.dest->create_handle("CCCC")});
    REQUIRE_THROWS_AS(combiner.merge_path_fragments(), FragmentMergeError);
    REQUIRE_THROWS_WITH(combiner.merge_path_fragments(), Catch::Contains("HG002#1#chr1#0"));
}

TEST_CASE("combine rejects chunks of one base path that disagree on circularity", "[combine]") {
    vector<GraphCombiner::Input> inputs(2);
    inputs[0].name = "left.pg";
    inputs[0].graph = make_unique<bdsg::HashGraph>();
    add_path(*inputs[0].graph, PathSense::REFERENCE, "GRCh38", "chrM", 0,
             PathMetadata::NO_PHASE_BLOCK, subrange_t{0, 4},
             {inputs[0].graph->create_handle("AAAA")}, true);
    inputs[1].name = "right.pg";
    inputs[1].graph = make_unique<bdsg::HashGraph>();
    add_path(*inputs[1].graph, PathSense::REFERENCE, "GRCh38", "chrM", 0,
             PathMetadata::NO_PHASE_BLOCK, subrange_t{4, 7},
             {inputs[1].graph->create_handle("CCC")}, false);

    TestGraphCombiner combiner(SeamPolicy::COORD_TRIM, false, false);
    REQUIRE_THROWS_AS(combiner.combine(std::move(inputs)), InputConflictError);
}

}
}
