#include "catch.hpp"
#include <stdio.h>
#include <iostream>
#include "../min_width_int_vector.hpp"

namespace vg{
namespace unittest{
using namespace std;

    TEST_CASE("Array of ints added one at a time", "[minint]") {
        SECTION ("[0]") {
            min_width_int_vector_t minint_vector (1);
            minint_vector.push_back(0);
            REQUIRE(minint_vector.size() == 1);
            REQUIRE(minint_vector.at(0) == 0);
        }
        SECTION ("[1]") {
            min_width_int_vector_t minint_vector (1);
            minint_vector.push_back(1);
            REQUIRE(minint_vector.size() == 1);
            REQUIRE(minint_vector.at(0) == 1);
        }
        SECTION ("[1, 2]") {
            min_width_int_vector_t minint_vector(2);
            minint_vector.push_back(1);
            minint_vector.push_back(2);
            REQUIRE(minint_vector.size() == 2);
            REQUIRE(minint_vector.at(0) == 1);
            REQUIRE(minint_vector.at(1) == 2);
        }
        SECTION ("more values") {
            vector<size_t> values {1, 3243, 123634, 53454, 0};
            min_width_int_vector_t minint_vector(1+(size_t)std::floor(std::log2(123634)));
            for (auto& x : values) {
               minint_vector.push_back(x); 
            }
            assert(minint_vector.size() == values.size());
            for (size_t i = 0 ; i < values.size() ; i++) {
                assert(minint_vector.at(i) == values[i]);
            }
        }
    }
    TEST_CASE("Array of ints from vector", "[minint]") {
        SECTION ("[0]") {
            vector<size_t> original {0};
            min_width_int_vector_t minint_vector;
            minint_vector.from_vector(original);
            REQUIRE(minint_vector.size() == 1);
            REQUIRE(minint_vector.at(0) == 0);
            REQUIRE(minint_vector.get_bit_width() == 1);
        }
        SECTION ("[1]") {
            vector<size_t> original {1};
            min_width_int_vector_t minint_vector;
            minint_vector.from_vector(original);
            REQUIRE(minint_vector.size() == 1);
            REQUIRE(minint_vector.at(0) == 1);
            REQUIRE(minint_vector.get_bit_width() == 1);
        }
        SECTION ("[1, 2]") {
            vector<size_t> original {1, 2};
            min_width_int_vector_t minint_vector;
            minint_vector.from_vector(original);

            REQUIRE(minint_vector.size() == 2);
            REQUIRE(minint_vector.at(0) == 1);
            REQUIRE(minint_vector.at(1) == 2);
            REQUIRE(minint_vector.get_bit_width() == 2);
        }
        SECTION ("more values") {
            vector<size_t> values {1, 3243, 123634, 53454, 0};
            min_width_int_vector_t minint_vector (3);
            minint_vector.from_vector(values, 123634);
            REQUIRE(minint_vector.get_bit_width() == 1+(size_t)std::floor(std::log2(123634)));
            assert(minint_vector.size() == values.size());
            for (size_t i = 0 ; i < values.size() ; i++) {
                assert(minint_vector.at(i) == values[i]);
            }
        }
        SECTION ("more values without bitwidth") {
            vector<size_t> values {1, 3243, 123634, 53454, 0};
            min_width_int_vector_t minint_vector;
            minint_vector.from_vector(values);
            assert(minint_vector.size() == values.size());
            for (size_t i = 0 ; i < values.size() ; i++) {
                assert(minint_vector.at(i) == values[i]);
            }
            REQUIRE(minint_vector.get_bit_width() == 1+(size_t)std::floor(std::log2(123634)));
        }
    }
}
}
