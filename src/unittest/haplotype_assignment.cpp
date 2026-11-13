#include "../haplotype_assignment.hpp"

#include "../gbwt_helper.hpp"
#include "../integrated_snarl_finder.hpp"
#include "../snarl_distance_index.hpp"
#include "../vg.hpp"

#include "catch.hpp"

namespace vg {
namespace unittest {

namespace {

VG make_tiny_bubble_graph() {
    // 1 -> 2 -> 4 -> 5
    //  \-> 3 -/
    VG graph;
    Node* n1 = graph.create_node("A");
    Node* n2 = graph.create_node("C");
    Node* n3 = graph.create_node("G");
    Node* n4 = graph.create_node("T");
    Node* n5 = graph.create_node("A");
    graph.create_edge(n1, n2);
    graph.create_edge(n1, n3);
    graph.create_edge(n2, n4);
    graph.create_edge(n3, n4);
    graph.create_edge(n4, n5);
    return graph;
}

VG make_nested_bubble_graph() {
    // Node IDs. From n8 on, a node's ID is one less than its variable's number.
    //   10 -> 1 -> 2 -> 7 -> 9 -> 11
    //          \-> 3 -> (4|5) -> 6 -> 8 -/
    // The snarl between 3 and 6 is nested in the snarl between 1 and 9. The
    // flanking nodes 10 and 11 are long, so the top-level chain runs from 10
    // to 11. Without them the snarl finder makes the cycle through 1 and 9 the
    // top-level chain, and 1 and 9 become tips.
    VG graph;
    Node* n1 = graph.create_node("A");
    Node* n2 = graph.create_node("C");
    Node* n3 = graph.create_node("G");
    Node* n4 = graph.create_node("T");
    Node* n5 = graph.create_node("A");
    Node* n6 = graph.create_node("C");
    Node* n8 = graph.create_node("G");
    Node* n9 = graph.create_node("T");
    Node* n10 = graph.create_node("A");

    graph.create_edge(n1, n2);
    graph.create_edge(n1, n3);
    graph.create_edge(n2, n8);
    graph.create_edge(n8, n10);

    graph.create_edge(n3, n4);
    graph.create_edge(n3, n5);
    graph.create_edge(n4, n6);
    graph.create_edge(n5, n6);
    graph.create_edge(n6, n9);
    graph.create_edge(n9, n10);

    Node* left_flank = graph.create_node("ACGTACGTACGTACGTACGT");
    Node* right_flank = graph.create_node("TGCATGCATGCATGCATGCA");
    graph.create_edge(left_flank, n1);
    graph.create_edge(n10, right_flank);

    return graph;
}

gbwt::vector_type gbwt_path(std::initializer_list<nid_t> nodes) {
    gbwt::vector_type path;
    path.reserve(nodes.size());
    for (nid_t node_id : nodes) {
        path.push_back(gbwt::Node::encode(node_id, false));
    }
    return path;
}

Alignment make_alignment(std::initializer_list<nid_t> nodes) {
    Alignment aln;
    Path* path = aln.mutable_path();
    for (nid_t node_id : nodes) {
        Mapping* mapping = path->add_mapping();
        Position* pos = mapping->mutable_position();
        pos->set_node_id(node_id);
        pos->set_is_reverse(false);
        pos->set_offset(0);
    }
    return aln;
}

} // namespace

TEST_CASE("HaplotypeAssigner finds segments and snarl boundaries on a tiny bubble graph",
          "[giraffe][haplotype_assignment]") {

    VG graph = make_tiny_bubble_graph();

    // Three haplotypes. h0 and h2 take the same path.
    //   h0: 1-2-4-5
    //   h1: 1-3-4-5
    //   h2: 1-2-4-5
    std::vector<gbwt::vector_type> threads = {
        gbwt_path({1, 2, 4, 5}),
        gbwt_path({1, 3, 4, 5}),
        gbwt_path({1, 2, 4, 5})
    };
    gbwt::GBWT gbwt_index = get_gbwt(threads);
    gbwtgraph::GBZ gbz(std::move(gbwt_index), graph, nullptr);

    IntegratedSnarlFinder snarl_finder(graph);
    SnarlDistanceIndex dist_index;
    fill_in_distance_index(&dist_index, &graph, &snarl_finder);

    HaplotypeAssigner assigner(gbz, dist_index);

    SECTION("A path that two haplotypes follow is one segment with both") {
        Alignment aln = make_alignment({1, 2, 4, 5});
        auto segments = assigner.segments(aln);

        REQUIRE(segments.size() == 1);
        REQUIRE(segments[0].first_mapping == 0);
        REQUIRE(segments[0].past_last_mapping == 4);
        REQUIRE(segments[0].visit_count == 2);
        REQUIRE(segments[0].seq_ids.size() == 2);
        REQUIRE(segments[0].seq_ids[0] != segments[0].seq_ids[1]);
    }

    SECTION("The other branch has its single haplotype") {
        Alignment aln = make_alignment({1, 3, 4, 5});
        auto segments = assigner.segments(aln);

        REQUIRE(segments.size() == 1);
        REQUIRE(segments[0].visit_count == 1);
        REQUIRE(segments[0].seq_ids.size() == 1);
    }

    SECTION("A path inside a snarl gets the snarl's boundary nodes") {
        Alignment aln = make_alignment({2});
        auto segments = assigner.segments(aln);
        SnarlBoundaries boundaries = assigner.snarl_boundaries(aln);

        REQUIRE(segments.size() == 1);
        REQUIRE(segments[0].visit_count == 2);
        REQUIRE(boundaries.left.value_or(0) == 1);
        REQUIRE(boundaries.right.value_or(0) == 4);
    }

    SECTION("A path that no snarl contains gets no boundary nodes") {
        // 4 is the snarl's end boundary node, so it is outside the snarl.
        Alignment aln = make_alignment({2, 4});
        auto segments = assigner.segments(aln);
        SnarlBoundaries boundaries = assigner.snarl_boundaries(aln);

        REQUIRE(segments.size() == 1);
        REQUIRE(segments[0].visit_count == 2);
        REQUIRE_FALSE(boundaries.left.has_value());
        REQUIRE_FALSE(boundaries.right.has_value());
    }

    SECTION("The cap shortens seq_ids below visit_count") {
        assigner.max_haplotypes_to_report = 1;
        Alignment aln = make_alignment({1, 2, 4, 5});
        auto segments = assigner.segments(aln);

        REQUIRE(segments.size() == 1);
        REQUIRE(segments[0].visit_count == 2);
        REQUIRE(segments[0].seq_ids.size() == 1);
    }

    SECTION("A haplotype that ends inside the snarl is still in the segment") {
        // h3 visits 1 and 2 and stops there, before the snarl's end.
        std::vector<gbwt::vector_type> with_short = {
            gbwt_path({1, 2, 4, 5}),
            gbwt_path({1, 3, 4, 5}),
            gbwt_path({1, 2, 4, 5}),
            gbwt_path({1, 2})
        };
        gbwtgraph::GBZ gbz_with_short(get_gbwt(with_short), graph, nullptr);
        HaplotypeAssigner assigner_with_short(gbz_with_short, dist_index);

        Alignment aln = make_alignment({2});
        SnarlBoundaries boundaries = assigner_with_short.snarl_boundaries(aln);

        // The walk to the right goes on without h3. The segment has all three.
        REQUIRE(assigner_with_short.segments(aln)[0].visit_count == 3);
        REQUIRE(boundaries.left.value_or(0) == 1);
        REQUIRE(boundaries.right.value_or(0) == 4);
    }

    SECTION("An alignment with no path has no segments and no boundary nodes") {
        Alignment aln;
        SnarlBoundaries boundaries = assigner.snarl_boundaries(aln);

        REQUIRE(assigner.segments(aln).empty());
        REQUIRE_FALSE(boundaries.left.has_value());
        REQUIRE_FALSE(boundaries.right.has_value());
    }
}

TEST_CASE("HaplotypeAssigner finds snarl boundaries in nested snarls",
          "[giraffe][haplotype_assignment][nested]") {

    VG graph = make_nested_bubble_graph();

    // h0 takes the short branch. h1 and h3 take the long branch through 4, and
    // h2 takes it through 5.
    std::vector<gbwt::vector_type> threads = {
        gbwt_path({1, 2, 7, 9}),
        gbwt_path({1, 3, 4, 6, 8, 9}),
        gbwt_path({1, 3, 5, 6, 8, 9}),
        gbwt_path({1, 3, 4, 6, 8, 9})
    };
    gbwt::GBWT gbwt_index = get_gbwt(threads);
    gbwtgraph::GBZ gbz(std::move(gbwt_index), graph, nullptr);

    IntegratedSnarlFinder snarl_finder(graph);
    SnarlDistanceIndex dist_index;
    fill_in_distance_index(&dist_index, &graph, &snarl_finder);

    HaplotypeAssigner assigner(gbz, dist_index);

    SECTION("A path through the inner snarl has the haplotypes of its allele") {
        Alignment aln = make_alignment({3, 4, 6, 8});
        auto segments = assigner.segments(aln);

        REQUIRE(segments.size() == 1);
        REQUIRE(segments[0].visit_count == 2);
        REQUIRE(segments[0].seq_ids.size() == 2);
    }

    SECTION("A node of the inner snarl gets the inner snarl's boundary nodes") {
        Alignment aln = make_alignment({4});
        SnarlBoundaries boundaries = assigner.snarl_boundaries(aln);

        REQUIRE(assigner.segments(aln)[0].visit_count == 2);
        REQUIRE(boundaries.left.value_or(0) == 3);
        REQUIRE(boundaries.right.value_or(0) == 6);
    }

    SECTION("The walk stops where the haplotypes come from different nodes") {
        // Haplotypes reach 6 from both 4 and 5, so the left walk stops without
        // a boundary node.
        Alignment aln = make_alignment({6});
        SnarlBoundaries boundaries = assigner.snarl_boundaries(aln);

        REQUIRE(assigner.segments(aln)[0].visit_count == 3);
        REQUIRE_FALSE(boundaries.left.has_value());
        REQUIRE(boundaries.right.value_or(0) == 9);
    }

    SECTION("A path that leaves the inner snarl gets the outer snarl's boundary nodes") {
        // 4 is in both snarls, and 6 only in the outer one.
        Alignment aln = make_alignment({4, 6});
        SnarlBoundaries boundaries = assigner.snarl_boundaries(aln);

        REQUIRE(assigner.segments(aln)[0].visit_count == 2);
        REQUIRE(boundaries.left.value_or(0) == 1);
        REQUIRE(boundaries.right.value_or(0) == 9);
    }

    SECTION("The walk goes through a child snarl") {
        // Without h2, every haplotype at 3 goes on through 4, inside the inner
        // snarl, and the walk goes on to the outer snarl's end.
        std::vector<gbwt::vector_type> through_4 = {
            gbwt_path({1, 2, 7, 9}),
            gbwt_path({1, 3, 4, 6, 8, 9}),
            gbwt_path({1, 3, 4, 6, 8, 9})
        };
        gbwtgraph::GBZ gbz_through_4(get_gbwt(through_4), graph, nullptr);
        HaplotypeAssigner assigner_through_4(gbz_through_4, dist_index);

        Alignment aln = make_alignment({3});
        SnarlBoundaries boundaries = assigner_through_4.snarl_boundaries(aln);

        REQUIRE(assigner_through_4.segments(aln)[0].visit_count == 2);
        REQUIRE(boundaries.left.value_or(0) == 1);
        REQUIRE(boundaries.right.value_or(0) == 9);
    }

    SECTION("A path that no haplotype follows end to end has two segments and no boundary nodes") {
        // No haplotype goes from 3 straight to 6.
        Alignment aln = make_alignment({3, 6});
        auto segments = assigner.segments(aln);
        SnarlBoundaries boundaries = assigner.snarl_boundaries(aln);

        REQUIRE(segments.size() == 2);
        REQUIRE(segments[0].past_last_mapping == 1);
        REQUIRE(segments[1].first_mapping == 1);
        REQUIRE_FALSE(boundaries.left.has_value());
        REQUIRE_FALSE(boundaries.right.has_value());
    }
}

TEST_CASE("haplotype_name names a haplotype by sample, phase and contig, and a generic path by its contig",
          "[giraffe][haplotype_assignment]") {

    gbwt::Metadata metadata;
    metadata.setSamples(std::vector<std::string>{gbwtgraph::GENERIC_PATH_SAMPLE_NAME, "HG002"});
    metadata.setContigs(std::vector<std::string>{"chrM"});
    // A generic path, then two fragments of one haplotype.
    metadata.addPath(0, 0, 0, 0);
    metadata.addPath(1, 0, 2, 0);
    metadata.addPath(1, 0, 2, 1);

    REQUIRE(haplotype_name(metadata, 0) == "chrM");
    REQUIRE(haplotype_name(metadata, 1) == "HG002#2#chrM");
    REQUIRE(haplotype_name(metadata, 2) == "HG002#2#chrM");
}

} // namespace unittest
} // namespace vg
