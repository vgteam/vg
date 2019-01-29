///
/// \file vpkg.cpp
///  
/// Tests for VPKG packaging format
///


#include "catch.hpp"

#include "../stream/vpkg.hpp"
#include "../xg.hpp"
#include "../json2pb.h"
#include <gcsa/gcsa.h>
#include <sstream>
#include <tuple>
#include <memory>

namespace vg {
namespace unittest {
using namespace std;

TEST_CASE("We can serialize and re-read an empty GCSA", "[vpkg][gcsa]") {

    gcsa::GCSA empty_index;
    
    stringstream ss;
    
    stream::VPKG::save(empty_index, ss);
    
    // There should be some data
    REQUIRE(ss.str().size() != 0);

    tuple<unique_ptr<gcsa::GCSA>> loaded = stream::VPKG::load_all<gcsa::GCSA>(ss);
    
    // We should get something allocated
    REQUIRE(get<0>(loaded).get() != nullptr);

}

TEST_CASE("We can read and write XG", "[vpkg][xg]") {

    string graph_json = R"(
    {"node":[{"id":1,"sequence":"GATT"},
    {"id":2,"sequence":"ACA"}],
    "edge":[{"to":2,"from":1}]}
    )";
    
    // Load the JSON
    Graph proto_graph;
    json2pb(proto_graph, graph_json.c_str(), graph_json.size());
    
    // Build the xg index
    xg::XG xg_index(proto_graph);

    stringstream ss;
    
    SECTION("We can read from a VPKG-wrapped stream") {
        stream::VPKG::save(xg_index, ss);
        
        // There should be some data
        REQUIRE(ss.str().size() != 0);
        
        unique_ptr<xg::XG> loaded = stream::VPKG::load_one<xg::XG>(ss);
        
        // Make sure we got something
        REQUIRE(loaded.get() != nullptr);
        
        // Make sure it is the thing we saved
        REQUIRE(loaded->node_size() == 2);
        REQUIRE(loaded->get_sequence(loaded->get_handle(1, false)) == "GATT");
    }
    
    SECTION("We can read from a bare stream") {
        xg_index.serialize(ss);
        
        // There should be some data
        REQUIRE(ss.str().size() != 0);
        
        unique_ptr<xg::XG> loaded = stream::VPKG::load_one<xg::XG>(ss);
        
        // Make sure we got something
        REQUIRE(loaded.get() != nullptr);
        
        // Make sure it is the thing we saved
        REQUIRE(loaded->node_size() == 2);
        REQUIRE(loaded->get_sequence(loaded->get_handle(1, false)) == "GATT");
    }
}

}

}

