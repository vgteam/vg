/// \file mem.cpp
///  
/// unit tests for MEMs and their clustering algorithms
///

#include <iostream>
#include "vg/io/json2pb.h"
#include <vg/vg.pb.h>
#include "../mem.hpp"
#include "../cluster.hpp"
#include "catch.hpp"

namespace vg {
namespace unittest {
    
TEST_CASE( "ShuffledPairs loops over all pairs", "[mem]" ) {
    
    SECTION( "ShuffledPairs finds the right number of pairs for 0 items" ) {
        
        ShuffledPairs pairs(0);
        
        set<pair<size_t, size_t>> found;
        
        for (auto pair : pairs) {
            found.insert(pair);
            REQUIRE(pair.first < pair.second);
        }
        
        REQUIRE(found.size() == 0);
    }
    
    SECTION( "ShuffledPairs finds the right number of pairs for 1 items" ) {
        
        ShuffledPairs pairs(1);
        
        set<pair<size_t, size_t>> found;
        
        for (auto pair : pairs) {
            found.insert(pair);
            REQUIRE(pair.first < pair.second);
        }
        
        REQUIRE(found.size() == 0);
    }
    
    SECTION( "ShuffledPairs finds the right number of pairs for 2 items" ) {
        
        ShuffledPairs pairs(2);
        
        set<pair<size_t, size_t>> found;
        
        for (auto pair : pairs) {
            found.insert(pair);
            REQUIRE(pair.first < pair.second);
        }
        
        REQUIRE(found.size() == 1);
    }
    
    SECTION( "ShuffledPairs finds the right number of pairs for 3 items" ) {
        
        ShuffledPairs pairs(3);
        
        set<pair<size_t, size_t>> found;
        
        for (auto pair : pairs) {
            found.insert(pair);
            REQUIRE(pair.first < pair.second);
        }
        
        REQUIRE(found.size() == 3);
    }
    
    SECTION( "ShuffledPairs finds the right number of pairs for 4 items" ) {
        
        ShuffledPairs pairs(4);
        
        set<pair<size_t, size_t>> found;
        
        for (auto pair : pairs) {
            found.insert(pair);
            REQUIRE(pair.first < pair.second);
        }
        
        REQUIRE(found.size() == 6);
    }
    
    SECTION( "ShuffledPairs finds the right number of pairs for 5 items" ) {
        
        ShuffledPairs pairs(5);
        
        set<pair<size_t, size_t>> found;
        
        for (auto pair : pairs) {
            found.insert(pair);
            REQUIRE(pair.first < pair.second);
        }
        
        REQUIRE(found.size() == 5 * 4 / 2);
    }
    
    SECTION( "ShuffledPairs finds the right number of pairs for 6 items" ) {
        
        ShuffledPairs pairs(6);
        
        set<pair<size_t, size_t>> found;
        
        for (auto pair : pairs) {
            found.insert(pair);
            REQUIRE(pair.first < pair.second);
        }
        
        REQUIRE(found.size() == 6 * 5 / 2);
    }
    
    SECTION( "ShuffledPairs finds the right number of pairs for 8 items" ) {
        
        ShuffledPairs pairs(8);
        
        set<pair<size_t, size_t>> found;
        
        for (auto pair : pairs) {
            found.insert(pair);
            REQUIRE(pair.first < pair.second);
        }
        
        REQUIRE(found.size() == 8 * 7 / 2);
    }
    
    SECTION( "ShuffledPairs finds the right number of pairs for 9 items" ) {
        
        ShuffledPairs pairs(9);
        
        set<pair<size_t, size_t>> found;
        
        for (auto pair : pairs) {
            REQUIRE(found.count(pair) == 0);
            found.insert(pair);
            REQUIRE(pair.first < pair.second);
        }
        
        REQUIRE(found.size() == 9 * 8 / 2);
    }
    
    SECTION( "ShuffledPairs finds the right number of pairs for 10 items" ) {
        
        ShuffledPairs pairs(10);
        
        set<pair<size_t, size_t>> found;
        
        for (auto pair : pairs) {
            REQUIRE(found.count(pair) == 0);
            found.insert(pair);
            REQUIRE(pair.first < pair.second);
        }
        
        REQUIRE(found.size() == 10 * 9 / 2);
    }
    
    SECTION( "ShuffledPairs finds the right number of pairs for 100 items" ) {
        
        ShuffledPairs pairs(100);
        
        set<pair<size_t, size_t>> found;
        
        for (auto pair : pairs) {
            REQUIRE(found.count(pair) == 0);
            found.insert(pair);
            REQUIRE(pair.first < pair.second);
        }
        
        REQUIRE(found.size() == 100 * 99 / 2);
    }
    
}
}
}
