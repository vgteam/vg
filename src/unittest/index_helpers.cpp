/** \file
 *
 * Unit tests for utility classes and functions for working with GBWT, GBWTGraph, and GCSA2.
 */

#include "../gbwt_helper.hpp"

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
        short_path, alt_path, short_path
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

TEST_CASE("Parallel GBWT reconstruction", "[index_helpers]") {
    std::vector<gbwt::vector_type> source {
        sample_1_a, sample_1_b, sample_2_a, sample_2_b,
    };
    gbwt::GBWT index = get_gbwt(source);

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

    SECTION("single-threaded") {
        RebuildParameters parameters;
        gbwt::GBWT rebuilt = rebuild_gbwt(index, jobs, node_to_job, parameters);
        check_paths(rebuilt, truth);
    }

    SECTION("multi-threaded") {
        RebuildParameters parameters;
        parameters.num_jobs = 2;
        gbwt::GBWT rebuilt = rebuild_gbwt(index, jobs, node_to_job, parameters);
        check_paths(rebuilt, truth);
    }
}

//------------------------------------------------------------------------------

}
}
