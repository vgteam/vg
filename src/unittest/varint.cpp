#include "catch.hpp"
#include <stdio.h>
#include <iostream>
#include "../varint.hpp"

namespace vg{
namespace unittest{
using namespace std;

    TEST_CASE("Array of ints", "[varint]") {
        SECTION ("[0]") {
            varint_vector_t varint_vector;
            varint_vector.add_value(0);
            pair<size_t, size_t> value_and_index = varint_vector.get_value_and_next_index(0);
            REQUIRE(value_and_index.first == 0);
            REQUIRE(value_and_index.second == 1);
        }
        SECTION ("[1]") {
            varint_vector_t varint_vector;
            varint_vector.add_value(1);
            pair<size_t, size_t> value_and_index = varint_vector.get_value_and_next_index(0);
            REQUIRE(value_and_index.first == 1);
            REQUIRE(value_and_index.second == 1);
        }
        SECTION ("[1, 2]") {
            varint_vector_t varint_vector;
            varint_vector.add_value(1);
            varint_vector.add_value(2);
            pair<size_t, size_t> value_and_index = varint_vector.get_value_and_next_index(0);
            REQUIRE(value_and_index.first == 1);
            REQUIRE(value_and_index.second == 1);
            value_and_index = varint_vector.get_value_and_next_index(1);
            REQUIRE(value_and_index.first == 2);
            REQUIRE(value_and_index.second == 2);
        }
        SECTION ("more values") {
            cerr << endl;
            vector<size_t> values {1, 56435345, 23423, 5, 123498275, 0, 213, 14253452324, std::numeric_limits<size_t>::max(), 0, 23123241234234, std::numeric_limits<size_t>::max()-1};
            varint_vector_t varint_vector;
            for (auto& x : values) {
               varint_vector.add_value(x); 
            }
            cerr << endl;
            size_t index = 0;//index in the varint vector
            size_t i = 0; //index in values
            while (i < values.size()) {
                pair<size_t, size_t> value_and_index = varint_vector.get_value_and_next_index(index);
                REQUIRE(value_and_index.first == values[i]);
                cerr << value_and_index.first << endl;
                index = value_and_index.second;
                i++;
            }
            REQUIRE(i == values.size());
        }
    }
}
}
