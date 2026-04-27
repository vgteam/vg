/// \file longest_overlap.cpp
///  
/// Unit tests for the longest overlap function
///

#include <iostream>
#include <string>
#include "../longest_overlap.hpp"
#include "../utility.hpp"
#include "catch.hpp"


namespace vg {
namespace unittest {
using namespace std;

TEST_CASE("Longest overlap algorithm works on curated inputs", "[string][z]") {
    
    SECTION("Simple case") {
        string str1 = "ACG";
        string str2 = "CGT";
        REQUIRE(longest_overlap(str1, str2) == 2);
        REQUIRE(longest_overlap(str2, str1) == 0);
    }

    SECTION("Empty cases") {
        string str1 = "ACG";
        string str2 = "";
        REQUIRE(longest_overlap(str1, str2) == 0);
        REQUIRE(longest_overlap(str2, str1) == 0);
        REQUIRE(longest_overlap(str2, str2) == 0);
    }

    SECTION("Random cases") {
        for (size_t i = 0; i < 500; ++i) {
            string str1 = pseudo_random_sequence(10, 2 * i);
            string str2 = pseudo_random_sequence(10, 2 * i + 1);
            size_t l = 0;
            for (size_t i = min(str1.size(), str2.size()); i != 0; --i) {
                if (str1.substr(str1.size() - i, i) == str2.substr(0, i)) {
                    l = i;
                    break;
                }
            }
            REQUIRE(longest_overlap(str1, str2) == l);
        }

    }
    
}
   
}
}
        
