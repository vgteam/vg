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

    SECTION("simple replacements") {
        std::vector<gbwt::vector_type> source {
            short_path, alt_path, short_path
        };
        gbwt::GBWT index = get_gbwt(source);
        std::vector<std::pair<gbwt::vector_type, gbwt::vector_type>> mappings {
            { { 8 }, { } }, // delete 4
            { { 12, 14 }, { 14 } }, // delete 6 if followed by 7
            { { 12, 16 }, { 12, 20, 16 } }, // visit 10 between 6 and 8
        };
        std::vector<gbwt::vector_type> truth {
            { 2, 10, 14, 18 },
            { 2, 4, 10, 12, 20, 16, 18 },
            { 2, 10, 14, 18 },
        };
        index = rebuild_gbwt(index, mappings);
        check_paths(index, truth);
    }

    SECTION("reverse replacements") {
        std::vector<gbwt::vector_type> source {
            short_path, alt_path, short_path
        };
        gbwt::GBWT index = get_gbwt(source);
        std::vector<std::pair<gbwt::vector_type, gbwt::vector_type>> mappings {
            { { 9 }, { } }, // delete 4
            { { 15, 13 }, { 15 } }, // delete 6 if followed by 7
            { { 17, 13 }, { 17, 21, 13 } }, // visit 10 between 6 and 8
        };
        std::vector<gbwt::vector_type> truth {
            { 2, 10, 14, 18 },
            { 2, 4, 10, 12, 20, 16, 18 },
            { 2, 10, 14, 18 },
        };
        index = rebuild_gbwt(index, mappings);
        check_paths(index, truth);
    }

    SECTION("impossible replacements") {
        std::vector<gbwt::vector_type> source {
            short_path, alt_path, short_path
        };
        gbwt::GBWT index = get_gbwt(source);
        std::vector<std::pair<gbwt::vector_type, gbwt::vector_type>> mappings {
            { { 6 }, { } }, // delete 3
            { { 4, 10 }, { 10 } }, // delete 2 if followed by 5
            { { 18, 20 }, { 18, 10, 20 } }, // visit 5 between 9 and 10
        };
        index = rebuild_gbwt(index, mappings);
        check_paths(index, source);
    }
}

//------------------------------------------------------------------------------

}
}
