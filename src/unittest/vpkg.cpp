///
/// \file vpkg.cpp
///  
/// Tests for VPKG packaging format
///


#include "catch.hpp"

#include <vg/io/vpkg.hpp>
#include "../xg.hpp"
#include "../seed_clusterer.hpp"
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
    
    vg::io::VPKG::save(empty_index, ss);
    
    // There should be some data
    REQUIRE(ss.str().size() != 0);

    tuple<unique_ptr<gcsa::GCSA>> loaded = vg::io::VPKG::load_all<gcsa::GCSA>(ss);
    
    // We should get something allocated
    REQUIRE(get<0>(loaded).get() != nullptr);

}

TEST_CASE("We can read and write XG", "[vpkg][handlegraph][xg]") {

    string graph_json = R"(
    {"node":[{"id":1,"sequence":"GATT"},
    {"id":2,"sequence":"ACA"}],
    "edge":[{"to":2,"from":1}]}
    )";
    
    // Load the JSON
    Graph proto_graph;
    json2pb(proto_graph, graph_json.c_str(), graph_json.size());
    
    // Build the xg index
    XG xg_index(proto_graph);

    stringstream ss;
    
    SECTION("We can read from a VPKG-wrapped stream as XG") {
        vg::io::VPKG::save(xg_index, ss);
        
        // There should be some data
        REQUIRE(ss.str().size() != 0);
        
        unique_ptr<XG> loaded = vg::io::VPKG::load_one<XG>(ss);
        
        // Make sure we got something
        REQUIRE(loaded.get() != nullptr);
        
        // Make sure it is the thing we saved
        REQUIRE(loaded->get_node_count() == 2);
        REQUIRE(loaded->get_sequence(loaded->get_handle(1, false)) == "GATT");
    }
    
    SECTION("We can read from a VPKG-wrapped stream as a HandleGraph") {
        vg::io::VPKG::save(xg_index, ss);
        
        // There should be some data
        REQUIRE(ss.str().size() != 0);
        
        unique_ptr<HandleGraph> loaded = vg::io::VPKG::load_one<HandleGraph>(ss);
        
        // Make sure we got something
        REQUIRE(loaded.get() != nullptr);
        
        // Make sure it is the thing we saved
        REQUIRE(loaded->get_node_count() == 2);
        REQUIRE(loaded->get_sequence(loaded->get_handle(1, false)) == "GATT");
    }
    
    SECTION("We can read from a bare stream") {
        xg_index.serialize(ss);
        
        // There should be some data
        REQUIRE(ss.str().size() != 0);
        
        unique_ptr<XG> loaded = vg::io::VPKG::load_one<XG>(ss);
        
        // Make sure we got something
        REQUIRE(loaded.get() != nullptr);
        
        // Make sure it is the thing we saved
        REQUIRE(loaded->get_node_count() == 2);
        REQUIRE(loaded->get_sequence(loaded->get_handle(1, false)) == "GATT");
    }
}

TEST_CASE("We can read a VG from an empty file", "[vpkg][vg][empty]") {
    stringstream ss;
    unique_ptr<vg::VG> loaded = vg::io::VPKG::try_load_one<vg::VG>(ss);
    
    // It should be empty but exist because this type is default constructible
    REQUIRE(std::is_default_constructible<vg::VG>::value);
    REQUIRE(loaded.get() != nullptr);
    REQUIRE(loaded->get_node_count() == 0);
}

TEST_CASE("We cannot read a SnarlSeedClusterer from an empty file", "[vpkg][snarlseedclusterer][empty]") {
    stringstream ss;
    unique_ptr<SnarlSeedClusterer> loaded = vg::io::VPKG::try_load_one<SnarlSeedClusterer>(ss);
    
    // It should be null because this type is not default constructible
    REQUIRE(!std::is_default_constructible<SnarlSeedClusterer>::value);
    REQUIRE(loaded.get() == nullptr);
}

TEST_CASE("We cannot read a base HandleGraph from an empty file", "[vpkg][handlegraph][empty]") {
    stringstream ss;
    unique_ptr<HandleGraph> loaded = vg::io::VPKG::try_load_one<HandleGraph>(ss);
    
    // It should be null because this type is not default constructible (because it is abstract)
    REQUIRE(!std::is_default_constructible<HandleGraph>::value);
    REQUIRE(loaded.get() == nullptr);
}

TEST_CASE("We can read VG from a VPKG-wrapped stream as a HandleGraph", "[vpkg][handlegraph][vg]") {
    string graph_json = R"(
    {"node":[{"id":1,"sequence":"GATT"},
    {"id":2,"sequence":"ACA"}],
    "edge":[{"to":2,"from":1}]}
    )";
    
    // Load the JSON
    Graph proto_graph;
    json2pb(proto_graph, graph_json.c_str(), graph_json.size());
    
    // Build the VG
    vg::VG vg_graph(proto_graph);
    
    // Save it
    stringstream ss;
    vg::io::VPKG::save(vg_graph, ss);
    
    // There should be some data
    REQUIRE(ss.str().size() != 0);
    
    unique_ptr<HandleGraph> loaded = vg::io::VPKG::load_one<HandleGraph>(ss);
    
    // Make sure we got something
    REQUIRE(loaded.get() != nullptr);
    
    // Make sure it is the thing we saved
    REQUIRE(loaded->get_node_count() == 2);
    REQUIRE(loaded->get_sequence(loaded->get_handle(1, false)) == "GATT");
}

TEST_CASE("We can read an empty VG as a HandleGraph", "[vpkg][handlegraph][vg][empty]") {
    string graph_json = "{}";
    
    // Load the JSON
    Graph proto_graph;
    json2pb(proto_graph, graph_json.c_str(), graph_json.size());
    
    // Build the VG
    vg::VG vg_graph(proto_graph);
    
    // Save it
    stringstream ss;
    vg::io::VPKG::save(vg_graph, ss);
    
    // There should be some data
    REQUIRE(ss.str().size() != 0);
    
    unique_ptr<HandleGraph> loaded = vg::io::VPKG::load_one<HandleGraph>(ss);
    
    // Make sure we got something
    REQUIRE(loaded.get() != nullptr);
    
    // Make sure it is the thing we saved
    REQUIRE(loaded->get_node_count() == 0);
}

}

}

