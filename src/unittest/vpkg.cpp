///
/// \file vpkg.cpp
///  
/// Tests for VPKG packaging format
///


#include "catch.hpp"

#include <vg/io/vpkg.hpp>
#include <bdsg/hash_graph.hpp>
#include <bdsg/packed_graph.hpp>
#include "xg.hpp"
#include "../vg.hpp"
#include "../snarl_seed_clusterer.hpp"
#include "vg/io/json2pb.h"
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

    tuple<unique_ptr<gcsa::GCSA>> loaded = vg::io::VPKG::try_load_first<gcsa::GCSA>(ss);
    
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
    xg::XG xg_index;
    xg_index.from_path_handle_graph(VG(proto_graph));

    stringstream ss;
    
    SECTION("We can read from a VPKG-wrapped stream as XG") {
        vg::io::VPKG::save(xg_index, ss);
        
        // There should be some data
        REQUIRE(ss.str().size() != 0);
        
        unique_ptr<xg::XG> loaded = vg::io::VPKG::load_one<xg::XG>(ss);
        
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
        
        unique_ptr<xg::XG> loaded = vg::io::VPKG::load_one<xg::XG>(ss);
        
        // Make sure we got something
        REQUIRE(loaded.get() != nullptr);
        
        // Make sure it is the thing we saved
        REQUIRE(loaded->get_node_count() == 2);
        REQUIRE(loaded->get_sequence(loaded->get_handle(1, false)) == "GATT");
    }
}

TEST_CASE("We cannot read a vg::VG from an empty file", "[vpkg][vg][empty]") {
    stringstream ss;
    unique_ptr<vg::VG> loaded = vg::io::VPKG::try_load_one<vg::VG>(ss);
    
    REQUIRE(loaded.get() == nullptr);
}

TEST_CASE("We cannot read a HandleGraph from an empty file", "[vpkg][handlegraph][empty]") {
    stringstream ss;
    unique_ptr<bdsg::HandleGraph> loaded = vg::io::VPKG::try_load_one<bdsg::HandleGraph>(ss);
    
    // We used to be able to load vg::VG graphs from empty files but we can't
    // really know that when loading a handle graph and picking the type.
    REQUIRE(loaded.get() == nullptr);
}

TEST_CASE("We cannot read a SnarlDistanceIndexClusterer from an empty file", "[vpkg][snarlseedclusterer][empty]") {
    stringstream ss;
    unique_ptr<SnarlDistanceIndexClusterer> loaded = vg::io::VPKG::try_load_one<SnarlDistanceIndexClusterer>(ss);
    
    // It should be null because this type is not default constructible
    REQUIRE(!std::is_default_constructible<SnarlDistanceIndexClusterer>::value);
    REQUIRE(loaded.get() == nullptr);
}

TEST_CASE("We cannot read a base HandleGraph from an empty file", "[vpkg][handlegraph][empty]") {
    stringstream ss;
    unique_ptr<HandleGraph> loaded = vg::io::VPKG::try_load_one<HandleGraph>(ss);
    
    // It should be null because this type is not default constructible (because it is abstract)
    REQUIRE(!std::is_default_constructible<HandleGraph>::value);
    REQUIRE(loaded.get() == nullptr);
}

TEST_CASE("We can read VG from a VPKG-wrapped stream as a VG", "[vpkg][handlegraph][vg]") {
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
    
    unique_ptr<VG> loaded = vg::io::VPKG::load_one<VG>(ss);
    
    // Make sure we got something
    REQUIRE(loaded.get() != nullptr);
    
    // Make sure it is the thing we saved
    REQUIRE(loaded->get_node_count() == 2);
    REQUIRE(loaded->get_sequence(loaded->get_handle(1, false)) == "GATT");
}

TEST_CASE("We can read VG from a VPKG-wrapped stream as a HandleGraph which is a HashGraph", "[vpkg][handlegraph][vg]") {
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
    
    // We should have a HashGraph, actually, and not a vg::VG
    REQUIRE(dynamic_cast<vg::VG*>(loaded.get()) == nullptr);
    REQUIRE(dynamic_cast<bdsg::HashGraph*>(loaded.get()) != nullptr);
    
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

TEST_CASE("We prefer to read a graph as the first provided type that matches", "[vpkg]") {
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
    
    // We use PackedGraph as a distractor here; HashGraph and VG can load from
    // vg Protobuf but PackedGraph can't as of this writing.
    auto loaded = vg::io::VPKG::try_load_first<bdsg::PackedGraph, HandleGraph, vg::VG>(ss);
    
    // Make sure we got the right thing
    REQUIRE(get<0>(loaded).get() == nullptr);
    REQUIRE(get<1>(loaded).get() != nullptr);
    REQUIRE(get<2>(loaded).get() == nullptr);
    
    // Make sure it is the thing we saved
    REQUIRE(get<1>(loaded)->get_node_count() == 2);
    REQUIRE(get<1>(loaded)->get_sequence(get<1>(loaded)->get_handle(1, false)) == "GATT");
}

TEST_CASE("We can send more than 2 GB of data through the bare-to-encapsulated conversion", "[vpkg]") {
    size_t DATA_SIZE = 3L * 1024 * 1024 * 1024;
    
    std::streamsize to_write = DATA_SIZE;
    // Make sure we don't overflow the amount we are supposed to be writing
    REQUIRE((size_t)to_write == (size_t)DATA_SIZE);
    
    // Map a giant buffer, since we can't just ask for a 3 GB string
    void* data = mmap(nullptr, DATA_SIZE, PROT_READ, MAP_PRIVATE | MAP_ANONYMOUS, -1, 0);
    REQUIRE(data != MAP_FAILED);
    
    std::atomic<size_t> bytes_seen;
    std::atomic<size_t> messages_seen;
    vg::io::message_consumer_function_t emit_message = [&](const string& message) {
        // Count the data we got
        size_t message_number = messages_seen++;
#ifdef debug
        #pragma omp critical (cerr)
        std::cerr << "Message " << message_number << " of size " << message.size() << " bytes observed." << std::endl;
#endif
        bytes_seen += message.size();
    };
    
    auto explain_stream = [](const ostream& o) {
        #pragma omp critical (cerr)
        cerr << "good: " << o.good() << " bad: " << o.bad() << " fail: " << o.fail() << " eof: " << o.eof() << endl;
    };
    
    vg::io::with_function_calling_stream(emit_message, [&](ostream& out) {
        // Write one giant chunk all at once
        
        // First turn on exceptions so that if something goes wrong we stop right away.
        out.exceptions(std::ifstream::badbit | std::ifstream::failbit);
        
#ifdef debug
        explain_stream(out);
        #pragma omp critical (cerr)
        std::cerr << "Writing " << to_write << " bytes at " << out.tellp() << "..." << std::endl;
#endif

        out.write((const char*)data, to_write);
        
#ifdef debug
        #pragma omp critical (cerr)
        std::cerr << "Writing over. Stream is now at " << out.tellp() << "." << std::endl;
        explain_stream(out);
#endif
    });
    
    // Unmap the giant buffer
    munmap(data, DATA_SIZE);
    
    // Make sure it all got there
    REQUIRE(bytes_seen == DATA_SIZE);

}

}

}

