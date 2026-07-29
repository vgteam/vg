/** \file
 *
 * Unit tests for utility classes and functions for working with GBWT, GBWTGraph, and GCSA2.
 */

#include "../gbwt_helper.hpp"

#include <gbwtgraph/utils.h>

#include "catch.hpp"



namespace vg {

namespace unittest {

//------------------------------------------------------------------------------

namespace {

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

gbwt::vector_type empty_path;

gbwt::vector_type sample_1_a {
    static_cast<gbwt::vector_type::value_type>(gbwt::Node::encode(11, false)),
    static_cast<gbwt::vector_type::value_type>(gbwt::Node::encode(12, false)),
    static_cast<gbwt::vector_type::value_type>(gbwt::Node::encode(14, false)),
    static_cast<gbwt::vector_type::value_type>(gbwt::Node::encode(15, false)),
    static_cast<gbwt::vector_type::value_type>(gbwt::Node::encode(17, false)),
};

gbwt::vector_type sample_2_a {
    static_cast<gbwt::vector_type::value_type>(gbwt::Node::encode(11, false)),
    static_cast<gbwt::vector_type::value_type>(gbwt::Node::encode(13, false)),
    static_cast<gbwt::vector_type::value_type>(gbwt::Node::encode(14, false)),
    static_cast<gbwt::vector_type::value_type>(gbwt::Node::encode(16, false)),
    static_cast<gbwt::vector_type::value_type>(gbwt::Node::encode(17, false)),
};

gbwt::vector_type sample_1_b {
    static_cast<gbwt::vector_type::value_type>(gbwt::Node::encode(21, false)),
    static_cast<gbwt::vector_type::value_type>(gbwt::Node::encode(22, false)),
    static_cast<gbwt::vector_type::value_type>(gbwt::Node::encode(24, false)),
    static_cast<gbwt::vector_type::value_type>(gbwt::Node::encode(23, true)),
    static_cast<gbwt::vector_type::value_type>(gbwt::Node::encode(21, true)),
};

gbwt::vector_type sample_2_b {
    static_cast<gbwt::vector_type::value_type>(gbwt::Node::encode(21, false)),
    static_cast<gbwt::vector_type::value_type>(gbwt::Node::encode(22, false)),
    static_cast<gbwt::vector_type::value_type>(gbwt::Node::encode(24, false)),
    static_cast<gbwt::vector_type::value_type>(gbwt::Node::encode(25, false)),
};

void check_paths(const gbwt::GBWT& index, const std::vector<gbwt::vector_type>& truth) {
    REQUIRE(index.sequences() == 2 * truth.size());
    for (gbwt::size_type i = 0; i < index.sequences(); i += 2) {
        gbwt::vector_type path = index.extract(i);
        REQUIRE(path == truth[i / 2]);
    }
}

} // anonymous namespace

//------------------------------------------------------------------------------

TEST_CASE("GBWT reconstruction", "[index_helpers]") {
    std::vector<gbwt::vector_type> source {
        short_path, alt_path, empty_path, short_path,
    };
    gbwt::GBWT original = get_gbwt(source);

    SECTION("simple replacements") {
        std::vector<RebuildJob::mapping_type> mappings {
            { { 8 }, { } }, // delete 4
            { { 12, 14 }, { 14 } }, // delete 6 if followed by 7
            { { 12, 16 }, { 12, 20, 16 } }, // visit 10 between 6 and 8
        };
        std::vector<gbwt::vector_type> truth {
            { 2, 10, 14, 18 },
            { 2, 4, 10, 12, 20, 16, 18 },
            { },
            { 2, 10, 14, 18 },
        };
        gbwt::GBWT index = rebuild_gbwt(original, mappings);
        check_paths(index, truth);
    }

    SECTION("reverse replacements") {
        std::vector<RebuildJob::mapping_type> mappings {
            { { 9 }, { } }, // delete 4
            { { 15, 13 }, { 15 } }, // delete 6 if followed by 7
            { { 17, 13 }, { 17, 21, 13 } }, // visit 10 between 6 and 8
        };
        std::vector<gbwt::vector_type> truth {
            { 2, 10, 14, 18 },
            { 2, 4, 10, 12, 20, 16, 18 },
            { },
            { 2, 10, 14, 18 },
        };
        gbwt::GBWT index = rebuild_gbwt(original, mappings);
        check_paths(index, truth);
    }

    SECTION("replacements with context") {
        std::vector<RebuildJob::mapping_type> mappings {
            { { 8 }, { 6, 8 } }, // add 3 before 4; do it only once
            { { 8 }, { 8, 24 } }, // add 12 after 4; this does not happen because 4 was already consumed
            { { 10, 12 }, { 22, 12 } }, // replace 5 with 11 if followed by 6
            { { 12, 16 }, { 12, 20, 16 } }, // visit 10 between 6 and 8; this works because 6 was not consumed
        };
        std::vector<gbwt::vector_type> truth {
            { 2, 6, 8, 22, 12, 14, 18 },
            { 2, 4, 6, 8, 22, 12, 20, 16, 18 },
            { },
            { 2, 6, 8, 22, 12, 14, 18 },
        };
        gbwt::GBWT index = rebuild_gbwt(original, mappings);
        check_paths(index, truth);
    }

    SECTION("impossible replacements") {
        std::vector<RebuildJob::mapping_type> mappings {
            { { 6 }, { } }, // delete 3
            { { 4, 10 }, { 10 } }, // delete 2 if followed by 5
            { { 18, 20 }, { 18, 10, 20 } }, // visit 5 between 9 and 10
        };
        gbwt::GBWT index = rebuild_gbwt(original, mappings);
        check_paths(index, source);
    }
}

TEST_CASE("Multiple rebuild_gbwt jobs", "[index_helpers]") {
    // We order the threads by (phase, contig), but rebuild_gbwt() will
    // reorder them by contig (with stable sorting).
    std::vector<gbwt::vector_type> source {
        sample_1_a, sample_1_b, sample_2_a, sample_2_b,
    };
    gbwt::GBWT index = get_gbwt(source);
    index.addMetadata();
    index.metadata.setSamples(1);
    index.metadata.setHaplotypes(2);
    index.metadata.setContigs(2);
    index.metadata.addPath(0, 0, 0, 0);
    index.metadata.addPath(0, 1, 0, 0);
    index.metadata.addPath(0, 0, 1, 0);
    index.metadata.addPath(0, 1, 1, 0);

    std::vector<RebuildJob> jobs {
        {
            {
                { { 24, 28, 30 }, { 24, 30 } }, // remove 14 in context 12 15.
                { { 30 }, { 30, 32 } }, // visit 16 after 15
            },
            7
        },
        {
            {
                { { 48, 47 }, { 46, 47 } }, // replace 24 with 23 if followed by reverse 23
                { { 48, 50 }, { 49, 50 } }, // flip 24 if followed by 25
            },
            5
        },
    };
    std::unordered_map<nid_t, size_t> node_to_job {
        { 11, 0 }, { 12, 0 }, { 13, 0 }, { 14, 0 }, { 15, 0 }, { 16, 0 }, { 17, 0 },
        { 21, 1 }, { 22, 1 }, { 23, 1 }, { 24, 1 }, { 25, 1 },
    };
    std::vector<gbwt::vector_type> truth {
        { 22, 24, 30, 32, 34 },
        { 22, 26, 28, 32, 34 },
        { 42, 44, 46, 47, 43 },
        { 42, 44, 49, 50 },
    };
    gbwt::Metadata true_metadata;
    true_metadata.setSamples(1);
    true_metadata.setHaplotypes(2);
    true_metadata.setContigs(2);
    true_metadata.addPath(0, 0, 0, 0);
    true_metadata.addPath(0, 0, 1, 0);
    true_metadata.addPath(0, 1, 0, 0);
    true_metadata.addPath(0, 1, 1, 0);

    SECTION("single-threaded") {
        RebuildParameters parameters;
        gbwt::GBWT rebuilt = rebuild_gbwt(index, jobs, node_to_job, parameters);
        check_paths(rebuilt, truth);
        REQUIRE(rebuilt.metadata == true_metadata);
    }

    SECTION("multi-threaded") {
        RebuildParameters parameters;
        parameters.num_jobs = 2;
        gbwt::GBWT rebuilt = rebuild_gbwt(index, jobs, node_to_job, parameters);
        check_paths(rebuilt, truth);
        REQUIRE(rebuilt.metadata == true_metadata);
    }
}

//------------------------------------------------------------------------------

namespace {

// A short haplotype fragment on contig 0 (nodes 11, 12, 14).
gbwt::vector_type wrap_hap_a {
    static_cast<gbwt::vector_type::value_type>(gbwt::Node::encode(11, false)),
    static_cast<gbwt::vector_type::value_type>(gbwt::Node::encode(12, false)),
    static_cast<gbwt::vector_type::value_type>(gbwt::Node::encode(14, false)),
};

// A short haplotype fragment on contig 1 (nodes 21, 22, 24).
gbwt::vector_type wrap_hap_b {
    static_cast<gbwt::vector_type::value_type>(gbwt::Node::encode(21, false)),
    static_cast<gbwt::vector_type::value_type>(gbwt::Node::encode(22, false)),
    static_cast<gbwt::vector_type::value_type>(gbwt::Node::encode(24, false)),
};

// Extract the forward path for a given path id from a bidirectional GBWT.
gbwt::vector_type extract_forward(const gbwt::GBWT& index, gbwt::size_type path_id) {
    return index.extract(gbwt::Path::encode(path_id, false));
}

} // anonymous namespace

TEST_CASE("double_origin_fragment doubles a path", "[index_helpers]") {
    SECTION("non-empty path") {
        gbwt::vector_type path = wrap_hap_a;
        double_origin_fragment(path);
        gbwt::vector_type expected = wrap_hap_a;
        expected.insert(expected.end(), wrap_hap_a.begin(), wrap_hap_a.end());
        REQUIRE(path == expected);
    }

    SECTION("empty path is unchanged") {
        gbwt::vector_type path;
        double_origin_fragment(path);
        REQUIRE(path.empty());
    }
}

TEST_CASE("wrap_haplotype_paths doubles origin fragments", "[index_helpers]") {
    // Two contigs, each with a single haplotype (sample "sampleA"), plus a
    // generic reference path on contig 0 that must never be wrapped.
    std::vector<gbwt::vector_type> source {
        wrap_hap_a, // path 0: haplotype on contig 0, count 0
        wrap_hap_b, // path 1: haplotype on contig 1, count 0
        wrap_hap_a, // path 2: generic reference on contig 0, count 0
    };
    gbwt::GBWT index = get_gbwt(source);
    index.addMetadata();
    index.metadata.setSamples(std::vector<std::string>{ "sampleA", gbwtgraph::GENERIC_PATH_SAMPLE_NAME });
    index.metadata.setContigs(std::vector<std::string>{ "chr1", "chrM" });
    index.metadata.setHaplotypes(1);
    index.metadata.addPath(0, 0, 0, 0); // sampleA, chr1, phase 0, count 0
    index.metadata.addPath(0, 1, 0, 0); // sampleA, chrM, phase 0, count 0
    index.metadata.addPath(1, 0, 0, 0); // generic,  chr1, phase 0, count 0

    SECTION("wrap one contig") {
        std::unordered_set<std::string> contigs { "chrM" };
        gbwt::GBWT wrapped = wrap_haplotype_paths(index, contigs);

        // chr1 haplotype untouched.
        REQUIRE(extract_forward(wrapped, 0) == wrap_hap_a);
        // chrM haplotype doubled.
        gbwt::vector_type expected_b = wrap_hap_b;
        expected_b.insert(expected_b.end(), wrap_hap_b.begin(), wrap_hap_b.end());
        REQUIRE(extract_forward(wrapped, 1) == expected_b);
        // Generic reference on chr1 untouched.
        REQUIRE(extract_forward(wrapped, 2) == wrap_hap_a);
        // Metadata preserved.
        REQUIRE(wrapped.metadata == index.metadata);
    }

    SECTION("wrapping a generic reference contig does not affect the reference") {
        // chr1 has both a haplotype and a generic path; wrapping chr1 should
        // double only the haplotype, not the generic reference path.
        std::unordered_set<std::string> contigs { "chr1" };
        gbwt::GBWT wrapped = wrap_haplotype_paths(index, contigs);

        gbwt::vector_type expected_a = wrap_hap_a;
        expected_a.insert(expected_a.end(), wrap_hap_a.begin(), wrap_hap_a.end());
        REQUIRE(extract_forward(wrapped, 0) == expected_a); // haplotype doubled
        REQUIRE(extract_forward(wrapped, 1) == wrap_hap_b); // chrM untouched
        REQUIRE(extract_forward(wrapped, 2) == wrap_hap_a); // generic untouched
    }

    SECTION("empty contig set is a no-op copy") {
        std::unordered_set<std::string> contigs;
        gbwt::GBWT wrapped = wrap_haplotype_paths(index, contigs);
        REQUIRE(extract_forward(wrapped, 0) == wrap_hap_a);
        REQUIRE(extract_forward(wrapped, 1) == wrap_hap_b);
        REQUIRE(extract_forward(wrapped, 2) == wrap_hap_a);
    }
}

TEST_CASE("wrap_haplotype_paths preserves reference samples", "[index_helpers]") {
    std::vector<gbwt::vector_type> source {
        wrap_hap_a, // path 0: reference sample on contig 0
        wrap_hap_b, // path 1: haplotype on contig 1
    };
    gbwt::GBWT index = get_gbwt(source);
    index.addMetadata();
    index.metadata.setSamples(std::vector<std::string>{ "refSample", "sampleA" });
    index.metadata.setContigs(std::vector<std::string>{ "chr1", "chrM" });
    index.metadata.setHaplotypes(2);
    index.metadata.addPath(0, 0, 0, 0); // refSample, chr1
    index.metadata.addPath(1, 1, 0, 0); // sampleA,   chrM
    index.tags.set(gbwtgraph::REFERENCE_SAMPLE_LIST_GBWT_TAG, "refSample");

    std::unordered_set<std::string> contigs { "chr1", "chrM" };
    gbwt::GBWT wrapped = wrap_haplotype_paths(index, contigs);

    // The reference-sample path on chr1 has REFERENCE sense, so it must not be
    // doubled even though chr1 is named.
    REQUIRE(extract_forward(wrapped, 0) == wrap_hap_a);
    // The haplotype on chrM is doubled.
    gbwt::vector_type expected_b = wrap_hap_b;
    expected_b.insert(expected_b.end(), wrap_hap_b.begin(), wrap_hap_b.end());
    REQUIRE(extract_forward(wrapped, 1) == expected_b);
    // Reference-samples tag preserved.
    REQUIRE(wrapped.tags.get(gbwtgraph::REFERENCE_SAMPLE_LIST_GBWT_TAG) == "refSample");
}

TEST_CASE("wrap_haplotype_paths rejects a missing origin fragment", "[index_helpers]") {
    // A single haplotype fragment on contig chrM whose count is 1 (its origin,
    // count 0, was truncated out of the graph). Wrapping must throw.
    std::vector<gbwt::vector_type> source { wrap_hap_b };
    gbwt::GBWT index = get_gbwt(source);
    index.addMetadata();
    index.metadata.setSamples(std::vector<std::string>{ "sampleA" });
    index.metadata.setContigs(std::vector<std::string>{ "chrM" });
    index.metadata.setHaplotypes(1);
    index.metadata.addPath(0, 0, 0, 1); // count 1: no origin fragment

    std::unordered_set<std::string> contigs { "chrM" };
    REQUIRE_THROWS_AS(wrap_haplotype_paths(index, contigs), std::runtime_error);
}

//------------------------------------------------------------------------------

}
}
