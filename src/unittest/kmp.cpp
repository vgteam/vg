/// \file annotation.cpp
///  
/// Unit tests for the annotations system on Alignments and MultipathAlignments
///

#include <iostream>
#include <string>
#include "../kmp.hpp"
#include "../utility.hpp"
#include "catch.hpp"


namespace vg {
namespace unittest {
using namespace std;

TEST_CASE("Knuth-Morris-Pratt implementation produces correct results", "[kmp]") {
    
    // examples taken from https://www.youtube.com/watch?v=4jY57Ehc14Y
    vector<pair<string, string>> searches;
    searches.emplace_back("ONIONIONSPL", "ONIONS");
    searches.emplace_back("QABCGABCABCD", "ABCD");
    searches.emplace_back("AAAAAAAAD", "AAAD");
    searches.emplace_back("TRAILTRAIN", "TRAIN");
    searches.emplace_back("ABCABC", "XYZ");
    searches.emplace_back("AABAAAC", "AAB");
    searches.emplace_back("AABAAAC", "AAC");
    for (size_t i = 0; i < 100; ++i) {
        searches.emplace_back(pseudo_random_sequence(20, i),
                              pseudo_random_sequence(5, 230897 * i + 98712));
        searches.emplace_back(searches.back().first, searches.back().first.substr(6, 10));
    }
    
    for (auto& search : searches) {
        
        string text, pattern;
        tie(text, pattern) = search;
        
        auto table = make_prefix_suffix_table(pattern.c_str(), pattern.size());
        size_t result = kmp_search(text.c_str(), text.size(),
                                   pattern.c_str(), pattern.size(),
                                   table);
        
        REQUIRE(result == text.find(pattern));
    }
    
}
   
}
}
        
