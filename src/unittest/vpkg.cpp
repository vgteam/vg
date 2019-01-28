///
/// \file vpkg.cpp
///  
/// Tests for VPKG packaging format
///


#include "catch.hpp"

#include "../stream/vpkg.hpp"
#include <gcsa/gcsa.h>
#include <sstream>
#include <tuple>
#include <memory>

namespace vg {
namespace unittest {
using namespace std;

TEST_CASE("We can serialize and re-read an empty GCSA", "[vpkg][gcsa]") {

    gcsa::GCSA empty_index;
    
    // Make sure we can save the empty index to a stream at all.
    stringstream teststream;
    empty_index.serialize(teststream);
    
    stringstream ss;
    
    stream::VPKG::save(empty_index, ss);
    
    // There should be some data
    REQUIRE(ss.str().size() != 0);

    tuple<unique_ptr<gcsa::GCSA>> loaded = stream::VPKG::load_all<gcsa::GCSA>(ss);
    
    // We should get something allocated
    REQUIRE(get<0>(loaded).get() != nullptr);

}

}

}

