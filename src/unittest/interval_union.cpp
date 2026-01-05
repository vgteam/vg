/// \file interval_union.cpp
///  
/// unit tests for IntervalUnion
///

#include "../interval_union.hpp"
#include "catch.hpp"

#include "bdsg/hash_graph.hpp"

namespace vg {
namespace unittest {

    TEST_CASE("IntervalUnion produces correct results on a small example", "[interval]") {
        
        IntervalUnion iunion;

        REQUIRE(iunion.total_size() == 0);
        REQUIRE(iunion.component_size() == 0);
        REQUIRE(iunion.overlap(50, 100) == 0);

        // (50, 100)
        iunion.add(50, 100);
        REQUIRE(iunion.total_size() == 50);
        REQUIRE(iunion.component_size() == 1);
        REQUIRE(iunion.overlap(0, 40) == 0);
        REQUIRE(iunion.overlap(0, 50) == 0);
        REQUIRE(iunion.overlap(100, 150) == 0);
        REQUIRE(iunion.overlap(110, 150) == 0);
        REQUIRE(iunion.overlap(50, 100) == 50);
        REQUIRE(iunion.overlap(75, 150) == 25);
        REQUIRE(iunion.overlap(25, 75) == 25);
        REQUIRE(iunion.overlap(60, 70) == 10);

        // (50, 100)
        iunion.add(60, 80);
        REQUIRE(iunion.total_size() == 50);
        REQUIRE(iunion.component_size() == 1);
        REQUIRE(iunion.overlap(0, 40) == 0);
        REQUIRE(iunion.overlap(0, 50) == 0);
        REQUIRE(iunion.overlap(100, 150) == 0);
        REQUIRE(iunion.overlap(110, 150) == 0);
        REQUIRE(iunion.overlap(50, 100) == 50);
        REQUIRE(iunion.overlap(75, 150) == 25);
        REQUIRE(iunion.overlap(25, 75) == 25);
        REQUIRE(iunion.overlap(60, 70) == 10);
       
        // (50, 100), (200, 250)
        iunion.add(200, 250);
        REQUIRE(iunion.total_size() == 100);
        REQUIRE(iunion.component_size() == 2);
        REQUIRE(iunion.overlap(50, 250) == 100);
        REQUIRE(iunion.overlap(40, 260) == 100);
        REQUIRE(iunion.overlap(75, 225) == 50);
        REQUIRE(iunion.overlap(40, 225) == 75);
        REQUIRE(iunion.overlap(70, 270) == 80);
        REQUIRE(iunion.overlap(20, 75) == 25);
        REQUIRE(iunion.overlap(75, 120) == 25);
        REQUIRE(iunion.overlap(225, 270) == 25);
        REQUIRE(iunion.overlap(20, 110) == 50);
        REQUIRE(iunion.overlap(190, 260) == 50);

        // (40, 100), (200, 250)
        iunion.add(40, 60);
        REQUIRE(iunion.total_size() == 110);
        REQUIRE(iunion.overlap(30, 260) == 110);
        REQUIRE(iunion.overlap(30, 100) == 60);
        REQUIRE(iunion.overlap(50, 100) == 50);
        REQUIRE(iunion.overlap(60, 110) == 40);

        // (40, 110), (200, 250)
        iunion.add(50, 110);
        REQUIRE(iunion.total_size() == 120);
        REQUIRE(iunion.component_size() == 2);
        REQUIRE(iunion.overlap(30, 260) == 120);
        REQUIRE(iunion.overlap(30, 100) == 60);
        REQUIRE(iunion.overlap(60, 120) == 50);

        // (40, 250)
        iunion.add(100, 210);
        REQUIRE(iunion.total_size() == 210);
        REQUIRE(iunion.component_size() == 1);
        REQUIRE(iunion.overlap(30, 260) == 210);
        REQUIRE(iunion.overlap(30, 50) == 10);
        REQUIRE(iunion.overlap(200, 300) == 50);
        REQUIRE(iunion.overlap(80, 130) == 50);
        REQUIRE(iunion.overlap(0, 10) == 0);

        // (10, 20), (40, 250)
        iunion.add(10, 20);
        REQUIRE(iunion.total_size() == 220);
        REQUIRE(iunion.component_size() == 2);
        REQUIRE(iunion.overlap(0, 30) == 10);
        REQUIRE(iunion.overlap(10, 100) == 70);

        // (10, 20), (40, 250), (280, 300)
        iunion.add(280, 300);
        REQUIRE(iunion.total_size() == 240);
        REQUIRE(iunion.component_size() == 3);
        REQUIRE(iunion.overlap(15, 295) == 230);


        // (10, 300)
        iunion.add(20, 280);
        REQUIRE(iunion.total_size() == 290);
        REQUIRE(iunion.component_size() == 1);
        REQUIRE(iunion.overlap(0, 30) == 20);
        REQUIRE(iunion.overlap(0, 10) == 0);
        REQUIRE(iunion.overlap(280, 310) == 20);
        REQUIRE(iunion.overlap(300, 310) == 0);
        REQUIRE(iunion.overlap(150, 200) == 50);
    }
}
}
