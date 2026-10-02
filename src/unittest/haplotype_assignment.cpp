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
    // Outer bubble topology:
    //   1 -> 2 -> 7 -> 9
    //    \-> 3 -> (4|5) -> 6 -> 8 -> 9
    // Inner bubble is (4|5) between 3 and 6, nested inside the long outer branch.
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

TEST_CASE("HaplotypeAssigner returns expected haplotypes on a tiny bubble graph",
          "[giraffe][haplotype_assignment]") {

    VG graph = make_tiny_bubble_graph();

    // Build a tiny GBWT with three haplotypes:
    //   h0: 1-2-4-5
    //   h1: 1-3-4-5
    //   h2: 1-2-4-5 (duplicate path, distinct sequence id)
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

    SECTION("Exact branch path reports exact known haplotype multiplicity") {
        Alignment aln = make_alignment({1, 2, 4, 5});
        auto result = assigner.assign(aln, false);

        REQUIRE(result.candidate_count == 2);
        REQUIRE(result.seq_ids.size() == 2);
        REQUIRE(result.truncated == false);
        REQUIRE(result.left_boundary_node == 0);
        REQUIRE(result.right_boundary_node == 0);
        REQUIRE(result.seq_ids[0] != result.seq_ids[1]);
    }

    SECTION("Other branch path reports its single haplotype") {
        Alignment aln = make_alignment({1, 3, 4, 5});
        auto result = assigner.assign(aln, false);

        REQUIRE(result.candidate_count == 1);
        REQUIRE(result.seq_ids.size() == 1);
        REQUIRE(result.truncated == false);
    }

    SECTION("Extension mode preserves or narrows haplotype set and reports boundaries") {
        Alignment aln = make_alignment({2, 4});
        auto no_extend = assigner.assign(aln, false);
        auto extended = assigner.assign(aln, true);

        REQUIRE(no_extend.candidate_count == 2);
        REQUIRE(extended.candidate_count <= no_extend.candidate_count);
        REQUIRE(extended.left_boundary_node == 1);
        REQUIRE(extended.right_boundary_node == 5);
    }

    SECTION("Truncation flag and list length are consistent with cap") {
        assigner.max_haplotypes_to_report = 1;
        Alignment aln = make_alignment({1, 2, 4, 5});
        auto result = assigner.assign(aln, false);

        REQUIRE(result.candidate_count == 2);
        REQUIRE(result.seq_ids.size() == 1);
        REQUIRE(result.truncated == true);
    }
}

TEST_CASE("HaplotypeAssigner handles nested snarls and extension edge cases",
          "[giraffe][haplotype_assignment][nested]") {

    VG graph = make_nested_bubble_graph();

    // h0: short outer branch
    // h1/h3: long branch via inner allele 4
    // h2: long branch via inner allele 5
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

    SECTION("Inner-allele path resolves to expected multiplicity") {
        Alignment aln = make_alignment({3, 4, 6, 8});
        auto result = assigner.assign(aln, false);

        REQUIRE(result.candidate_count == 2);
        REQUIRE(result.seq_ids.size() == 2);
        REQUIRE(result.truncated == false);
    }

    SECTION("Nested snarl extension finds inner boundaries for inner node") {
        Alignment aln = make_alignment({4});
        auto result = assigner.assign(aln, true);

        REQUIRE(result.candidate_count == 2);
        REQUIRE(result.left_boundary_node == 3);
        REQUIRE(result.right_boundary_node == 6);
    }

    SECTION("Extension stops conservatively at branching predecessor set") {
        // Node 6 has two valid predecessors (4 and 5) across different
        // haplotypes; left extension should stop without boundary.
        Alignment aln = make_alignment({6});
        auto result = assigner.assign(aln, true);

        REQUIRE(result.candidate_count == 3);
        REQUIRE(result.left_boundary_node == 0);
        REQUIRE(result.right_boundary_node == 9);
    }

    SECTION("Extension is monotonic on candidate count") {
        Alignment aln = make_alignment({3, 6});
        auto no_extend = assigner.assign(aln, false);
        auto extended = assigner.assign(aln, true);

        REQUIRE(extended.candidate_count <= no_extend.candidate_count);
    }
}

} // namespace unittest
} // namespace vg
