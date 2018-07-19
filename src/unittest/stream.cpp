/// \file stream.cpp
///  
/// Unit tests for stream functions

#include "../stream.hpp"

#include "vg.pb.h"

#include "catch.hpp"

#include <sstream>
#include <iostream>

namespace vg {
namespace unittest {
using namespace std;

TEST_CASE("Protobuf messages that are all default can be stored and retrieved", "[stream]") {
    stringstream datastream;
    
    // Write one empty message
    REQUIRE(stream::write<Graph>(datastream, 1, [](size_t i) {
       return Graph();
    }));
    
    // Look for it
    int seen = 0;
    stream::for_each<Graph>(datastream, [&](const Graph& item) {
        seen++;
    });
    
    // Make sure it comes back
    REQUIRE(seen == 1);
}

TEST_CASE("Protobuf messages can be written and read back", "[stream]") {
    stringstream datastream;
    
    // Define some functions to make and check fake Protobuf objects
    using message_t = Position;
    
    auto get_message = [&](size_t index) {
        message_t item;
        item.set_node_id(index);
#ifdef debug
        cerr << "Made item " << index << endl;
#endif
        return item;
    };
    
    size_t index_expected = 0;
    auto check_message = [&](const message_t& item) {
#ifdef debug
        cerr << "Read item " << item.node_id() << endl;
#endif
        REQUIRE(item.node_id() == index_expected);
        index_expected++;
    };
    
    // Serialize some objects
    REQUIRE(stream::write<message_t>(datastream, 10, get_message));
    
#ifdef debug
    // Dump the compressed data
    auto data = datastream.str();
    for (size_t i = 0; i < data.size(); i++) {
        ios state(nullptr);
        state.copyfmt(cerr);
        cerr << setfill('0') << setw(2) << hex << (int)(uint8_t)data[i] << " ";
        if (i % 8 == 7) {
            cerr << endl;
        }
        cerr.copyfmt(state);
    }
    cerr << endl;
#endif
   
    // Read them back
    stream::for_each<message_t>(datastream, check_message);

}

TEST_CASE("Multiple write calls work correctly on the same stream", "[stream]") {
    stringstream datastream;

    // Define some functions to make and check fake Protobuf objects
    using message_t = Position;
    
    size_t index_to_make = 0;
    auto get_message = [&](size_t index) {
        message_t item;
        item.set_node_id(index_to_make);
        index_to_make++;
        return item;
    };
    
    size_t index_expected = 0;
    auto check_message = [&](const message_t& item) {
        REQUIRE(item.node_id() == index_expected);
        index_expected++;
    };
    
    for (size_t i = 0; i < 10; i++) {
        // Serialize some objects
        REQUIRE(stream::write<message_t>(datastream, 1, get_message));
    }
    
    // Read them back
    stream::for_each<message_t>(datastream, check_message);

}

TEST_CASE("Single auto-chunking write calls work correctly", "[stream]") {
    stringstream datastream;

    // Define some functions to make and check fake Protobuf objects
    using message_t = Graph;
    
    id_t id_to_make = 0;
    auto get_message = [&](size_t start, size_t count) {
        message_t item;
        
        for(size_t i = 0; i < count; i++) {
            // Put nodes with those IDs in
            Node* added = item.add_node();
            added->set_id(id_to_make);
            id_to_make++;
#ifdef debug
            cerr << "Emit item " << added->id() << endl;
#endif
        }
        
        return item;
    };
    
    id_t id_expected = 0;
    auto check_message = [&](const message_t& item) {
        for (size_t i = 0; i < item.node_size(); i++) {
            // Make sure they come out in the same order
#ifdef debug
            cerr << "Found item " << item.node(i).id() << endl;
#endif
            REQUIRE(item.node(i).id() == id_expected);
            id_expected++;
        }
    };
    
    // Serialize some objects with dynamic chunking
    REQUIRE(stream::write<message_t>(datastream, 10, 1, get_message));
    
    // Read them back
    stream::for_each<message_t>(datastream, check_message);
    
    // Make sure we saw them all
    REQUIRE(id_expected == 10);

}

TEST_CASE("Multiple auto-chunking write calls work correctly on the same stream", "[stream]") {
    stringstream datastream;

    // Define some functions to make and check fake Protobuf objects
    using message_t = Graph;
    
    id_t id_to_make = 0;
    auto get_message = [&](size_t start, size_t count) {
        message_t item;
        
        for(size_t i = 0; i < count; i++) {
            // Put nodes with those IDs in
            Node* added = item.add_node();
            added->set_id(id_to_make);
            id_to_make++;
#ifdef debug
            cerr << "Emit item " << added->id() << endl;
#endif
        }
        
        return item;
    };
    
    id_t id_expected = 0;
    auto check_message = [&](const message_t& item) {
        for (size_t i = 0; i < item.node_size(); i++) {
            // Make sure they come out in the same order
#ifdef debug
            cerr << "Found item " << item.node(i).id() << endl;
#endif
            REQUIRE(item.node(i).id() == id_expected);
            id_expected++;
        }
    };
    
    for (size_t i = 0; i < 3; i++) {
        // Serialize some objects with dynamic chunking
        REQUIRE(stream::write<message_t>(datastream, 10, 1, get_message));
    }
    
    // Read them back
    stream::for_each<message_t>(datastream, check_message);
    
    // Make sure we saw them all
    REQUIRE(id_expected == 30);

}

}

}
