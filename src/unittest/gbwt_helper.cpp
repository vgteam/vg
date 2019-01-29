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

// Build a GBWT with three threads including a duplicate.
gbwt::GBWT build_gbwt_index() {
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

    return get_gbwt(gbwt_threads);
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
    int old_threads = omp_get_max_threads();
    omp_set_num_threads(1);
    auto lambda = [&found_kmers](const GBWTTraversal& window) {
        found_kmers.insert(kmer_type(window.traversal, window.seq));
    };
    for_each_kmer(xg_index, gbwt_index, 3, lambda);
    omp_set_num_threads(old_threads);

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

TEST_CASE("for_each_window() finds the correct windows", "[gbwt_helper]") {

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
    int old_threads = omp_get_max_threads();
    omp_set_num_threads(1);
    auto lambda = [&found_kmers](const GBWTTraversal& window) {
        found_kmers.insert(kmer_type(window.traversal, window.seq));
    };
    for_each_window(xg_index, gbwt_index, 3, lambda);
    omp_set_num_threads(old_threads);

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
