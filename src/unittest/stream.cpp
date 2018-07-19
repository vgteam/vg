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
        cerr << "Made item " << index << endl;
        return item;
    };
    
    size_t index_expected = 0;
    auto check_message = [&](const message_t& item) {
        cerr << "Read item " << item.node_id() << endl;
        REQUIRE(item.node_id() == index_expected);
        index_expected++;
    };
    
    // Serialize some objects
    REQUIRE(stream::write<message_t>(datastream, 10, get_message));
    
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
   
    // Read them back
    stream::for_each<message_t>(datastream, check_message);

}


}

}
