///
/// \file position.cpp
///  
/// Unit tests for Position and pos_t manipulation
///

#include "catch.hpp"
#include "../position.hpp"
#include "../vg.pb.h"

namespace vg {
namespace unittest {

TEST_CASE( "Position can be reversed", "[position]" ) {
    Position forward;
    forward.set_node_id(123);
    forward.set_is_reverse(false);
    
    SECTION( "At start of node" ) {
    
        forward.set_offset(0);

        SECTION( "On single-base node" ) {
            Position reversed = reverse(forward, 1);
            
            // Should come out as the first base on the reverse strand.
            REQUIRE(reversed.node_id() == forward.node_id());
            REQUIRE(reversed.is_reverse() == true);
            REQUIRE(reversed.offset() == 0);
            
        }
        
        SECTION( "On small even node" ) {
            Position reversed = reverse(forward, 2);
            
            // Should come out as the last base on the reverse strand.
            REQUIRE(reversed.node_id() == forward.node_id());
            REQUIRE(reversed.is_reverse() == true);
            REQUIRE(reversed.offset() == 1);
        }
        
        SECTION( "On small odd node" ) {
            Position reversed = reverse(forward, 3);
            
            // Should come out as the last base on the reverse strand.
            REQUIRE(reversed.node_id() == forward.node_id());
            REQUIRE(reversed.is_reverse() == true);
            REQUIRE(reversed.offset() == 2);
        }
        
        SECTION( "On large even node" ) {
            Position reversed = reverse(forward, 300);
            
            // Should come out as the last base on the reverse strand.
            REQUIRE(reversed.node_id() == forward.node_id());
            REQUIRE(reversed.is_reverse() == true);
            REQUIRE(reversed.offset() == 299);
        }
        
        SECTION( "On large odd node" ) {
            Position reversed = reverse(forward, 301);
            
            // Should come out as the last base on the reverse strand.
            REQUIRE(reversed.node_id() == forward.node_id());
            REQUIRE(reversed.is_reverse() == true);
            REQUIRE(reversed.offset() == 300);
        }
        
    }
    
    SECTION( "In the middle" ) {
    
        forward.set_offset(5);
        
        SECTION( "On even node, before center" ) {
            Position reversed = reverse(forward, 12);
            
            //           11
            // 012345678901
            //      *
            // 11          
            // 109876543210
            
            // Should come out as base 6
            REQUIRE(reversed.node_id() == forward.node_id());
            REQUIRE(reversed.is_reverse() == true);
            REQUIRE(reversed.offset() == 6);
            
        }

        SECTION( "On even node, past center" ) {
            Position reversed = reverse(forward, 10);
            
            // 0123456789
            //      *
            // 9876543210
            
            // Should come out as base 4
            REQUIRE(reversed.node_id() == forward.node_id());
            REQUIRE(reversed.is_reverse() == true);
            REQUIRE(reversed.offset() == 4);
            
        }
        
        SECTION( "On odd node, at center" ) {
            Position reversed = reverse(forward, 11);
            
            //           1
            // 01234567890
            //      *
            // 09876543210
            // 1
            
            // Should come out as base 5
            REQUIRE(reversed.node_id() == forward.node_id());
            REQUIRE(reversed.is_reverse() == true);
            REQUIRE(reversed.offset() == 5);
            
        }
    }
}

}
}
