///
/// \file position.cpp
///  
/// Unit tests for Position and pos_t manipulation
///

#include "catch.hpp"
#include "../position.hpp"
#include <vg/vg.pb.h>

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
            
            // 0 1
            //  X
            // 1 0
            
            // Should come out as after the last base on the reverse strand.
            REQUIRE(reversed.node_id() == forward.node_id());
            REQUIRE(reversed.is_reverse() == true);
            REQUIRE(reversed.offset() == 1);
            
        }
        
        SECTION( "On small even node" ) {
            Position reversed = reverse(forward, 2);
            
            // Should come out as after the last base on the reverse strand.
            REQUIRE(reversed.node_id() == forward.node_id());
            REQUIRE(reversed.is_reverse() == true);
            REQUIRE(reversed.offset() == 2);
        }
        
        SECTION( "On small odd node" ) {
            Position reversed = reverse(forward, 3);
            
            // Should come out as after the last base on the reverse strand.
            REQUIRE(reversed.node_id() == forward.node_id());
            REQUIRE(reversed.is_reverse() == true);
            REQUIRE(reversed.offset() == 3);
        }
        
        SECTION( "On large even node" ) {
            Position reversed = reverse(forward, 300);
            
            // Should come out as after the last base on the reverse strand.
            REQUIRE(reversed.node_id() == forward.node_id());
            REQUIRE(reversed.is_reverse() == true);
            REQUIRE(reversed.offset() == 300);
        }
        
        SECTION( "On large odd node" ) {
            Position reversed = reverse(forward, 301);
            
            // Should come out as after the last base on the reverse strand.
            REQUIRE(reversed.node_id() == forward.node_id());
            REQUIRE(reversed.is_reverse() == true);
            REQUIRE(reversed.offset() == 301);
        }
        
    }
    
    SECTION( "In the middle" ) {
    
        forward.set_offset(5);
        
        SECTION( "On even node, before center" ) {
            Position reversed = reverse(forward, 12);
            
            //                     1 1 1
            // 0 1 2 3 4 5 6 7 8 9 0 1 2
            //  X X X X X*X X X X X X X   
            // 1 1 1          
            // 2 1 0 9 8 7 6 5 4 3 2 1 0
            
            // Should come out as offset 7
            REQUIRE(reversed.node_id() == forward.node_id());
            REQUIRE(reversed.is_reverse() == true);
            REQUIRE(reversed.offset() == 7);
            
        }

        SECTION( "On even node, past center" ) {
            Position reversed = reverse(forward, 10);
            
            //                     1
            // 0 1 2 3 4 5 6 7 8 9 0
            //           *
            // 1 9 8 7 6 5 4 3 2 1 0
            // 0
            
            // Should come out as offset 5
            REQUIRE(reversed.node_id() == forward.node_id());
            REQUIRE(reversed.is_reverse() == true);
            REQUIRE(reversed.offset() == 5);
            
        }
        
        SECTION( "On odd node, at center" ) {
            Position reversed = reverse(forward, 11);
            
            // Should come out as offset 6
            REQUIRE(reversed.node_id() == forward.node_id());
            REQUIRE(reversed.is_reverse() == true);
            REQUIRE(reversed.offset() == 6);
            
        }
    }
}

}
}
