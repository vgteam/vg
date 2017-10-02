/**
 * \file 
 * unittest/handle.cpp: test cases for the implementations of the HandleGraph class.
 */

#include "catch.hpp"

#include "../handle.hpp"
#include "../json2pb.h"

#include <iostream>
#include <limits>
#include <vector>

namespace vg {
namespace unittest {

using namespace std;

TEST_CASE( "Handle utility functions work", "[handle]" ) {

    SECTION("Handles work like ints") {
        
        SECTION("Handles are int-sized") {
            REQUIRE(sizeof(handle_t) == sizeof(int64_t));
        }
        
        SECTION("Handles can hold a range of positive integer values") {
            for (int64_t i = 0; i < 100; i++) {
                REQUIRE(as_integer(as_handle(i)) == i);
            }      
            
            REQUIRE(as_integer(as_handle(numeric_limits<int64_t>::max())) == numeric_limits<int64_t>::max());
            
        }
        
    }
    
    SECTION("Handle equality works") {
        vector<handle_t> handles;
        
        for (size_t i = 0; i < 100; i++) {
            handles.push_back(as_handle(i));
        }
        
        for (size_t i = 0; i < handles.size(); i++) {
            for (size_t j = 0; j < handles.size(); j++) {
                if (i == j) {
                    REQUIRE(handles[i] == handles[j]);
                    REQUIRE(! (handles[i] != handles[j]));
                } else {
                    REQUIRE(handles[i] != handles[j]);
                    REQUIRE(! (handles[i] == handles[j]));
                }
            }
        }
    }
    
}

}
}
