/// \file alignment.cpp
///  
/// unit tests for Alignments and their utility functions
///

#include <iostream>
#include <string>
#include <unordered_map>
#include "vg/io/json2pb.h"
#include <vg/vg.pb.h>
#include "../algorithms/back_translate.hpp"
#include "catch.hpp"

namespace vg {
namespace unittest {
using namespace std;

class MockBackTranslation : public NamedNodeBackTranslation {
public:
    
    unordered_map<oriented_node_range_t, vector<oriented_node_range_t>> translation;
    unordered_map<nid_t, string> node_names;
    
    /**
     * Translate the given range of bases on the given orientation of the given
     * node in the current graph, to zero or more ranges on orientations of
     * nodes in some prior graph.
     */
    std::vector<oriented_node_range_t> translate_back(const oriented_node_range_t& range) const {
        return translation.at(range);
    }
    
    /**
     * Get the name of a node in the graph that translate_back() translates
     * into, given its number.
     */
    std::string get_back_graph_node_name(const nid_t& back_node_id) const {
        return node_names.at(back_node_id);
    }


};

TEST_CASE("A Path can be back-translated", "[algorithms][back_translate]") {

    // Make a path
    string path_string = R"(
        {
            "mapping": [
                {
                    "position": {"node_id": 1},
                    "edit": [
                        {"from_length": 1, "to_length": 1},
                        {"from_length": 1}
                    ]
                },
                {
                    "position": {"node_id": 2},
                    "edit": [
                        {"from_length": 1},
                        {"from_length": 1, "to_length": 1}
                    ]
                }
            ]
        }
    )";
    Path p;
    json2pb(p, path_string.c_str(), path_string.size());
    
    // Define a translation back to a named node space
    MockBackTranslation trans;
    trans.translation[oriented_node_range_t(1, false, 0, 2)] = {oriented_node_range_t(1, false, 5, 2)};
    trans.translation[oriented_node_range_t(2, false, 0, 2)] = {oriented_node_range_t(1, false, 7, 2)};
    trans.node_names[1] = "TheNode";
    
    // Translate the path
    vg::algorithms::back_translate_in_place(&trans, p);

    // Check the result
    REQUIRE(p.mapping_size() == 1);
    REQUIRE(p.mapping(0).edit_size() == 3);
    REQUIRE(p.mapping(0).position().node_id() == 0);
    REQUIRE(p.mapping(0).position().name() == "TheNode");
    REQUIRE(p.mapping(0).position().offset() == 5);
}


}
}
