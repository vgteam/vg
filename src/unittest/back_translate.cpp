/// \file alignment.cpp
///  
/// unit tests for Alignments and their utility functions
///

#include <iostream>
#include <string>
#include <unordered_map>
#include "vg/io/json2pb.h"
#include <vg/vg.pb.h>
#include <vg/io/alignment_io.hpp>
#include <bdsg/hash_graph.hpp>
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

TEST_CASE("An Alignment can be back-translated while converting to GAF", "[algorithms][back_translate]") {

    string alignment_string = R"(
        {
            "name": "francine",
            "mapping_quality": 30,
            "sequence": "GATTACA",
            "path": {"mapping": [
                {
                    "position": {"node_id": 2, "offset": 2},
                    "edit": [
                        {"from_length": 1, "to_length": 1},
                        {"from_length": 1}
                    ]
                },
                {
                    "position": {"node_id": 3},
                    "edit": [
                        {"from_length": 1, "to_length": 1},
                        {"to_length": 1, "sequence": "T"},
                        {"from_length": 1, "to_length": 1}
                    ]
                },
                {
                    "position": {"node_id": 4},
                    "edit": [
                        {"from_length": 1, "to_length": 1},
                        {"from_length": 2},
                        {"from_length": 2, "to_length": 2}
                    ]
                }
            ]}
        }
    )";
    
    Alignment a;
    json2pb(a, alignment_string.c_str(), alignment_string.size());
    
    bdsg::HashGraph g;
    
    handle_t h1 = g.create_handle("G");
    handle_t h2 = g.create_handle("GGGG");
    handle_t h3 = g.create_handle("AT");
    handle_t h4 = g.create_handle("ACACAAA");
    handle_t h5 = g.create_handle("A");
    
    g.create_edge(h1, h2);
    g.create_edge(h2, h3);
    g.create_edge(h3, h4);
    g.create_edge(h4, h5);

    // Define a translation back to a named node space
    MockBackTranslation trans;
    trans.translation[oriented_node_range_t(2, false, 2, 2)] = {oriented_node_range_t(1, false, 3, 2)};
    trans.translation[oriented_node_range_t(3, false, 0, 2)] = {oriented_node_range_t(1, false, 5, 2)};
    trans.translation[oriented_node_range_t(4, false, 0, 5)] = {oriented_node_range_t(2, false, 0, 5)};
    // Make sure to include reverse strand translation info for segment length sniffing.
    trans.translation[oriented_node_range_t(2, true, 0, 0)] = {oriented_node_range_t(1, true, 2, 0)};
    trans.translation[oriented_node_range_t(3, true, 0, 0)] = {oriented_node_range_t(1, true, 0, 0)};
    trans.translation[oriented_node_range_t(4, true, 0, 0)] = {oriented_node_range_t(2, true, 1, 0)};
    trans.node_names[1] = "FirstSegment";
    trans.node_names[2] = "SecondSegment";
    
    SECTION("Translating GAF generation produces the right GAF") {
        auto node_space = vg::io::alignment_to_gaf(g, a);
        stringstream s;
        s << node_space;
        
        // See column definitions at https://github.com/lh3/gfatools/blob/master/doc/rGFA.md#the-graph-alignment-format-gaf
        // Alignment block length is longest involved sequence.
        // Note that we combine adjacent duplicate operations in cs across node boundaries.
        // Note that end position is 0-based inclusive
        REQUIRE(s.str() == "francine\t7\t0\t7\t+\t>2>3>4\t13\t2\t10\t6\t8\t30\tcs:Z::1-G:1+T:2-CA:2");
    }
    
    SECTION("Translating GAF generation with a translation produces the right GAF") {
        auto back_translated = vg::io::alignment_to_gaf(g, a, &trans);
        stringstream s;
        s << back_translated;
        
        REQUIRE(s.str() == "francine\t7\t0\t7\t+\t>FirstSegment>SecondSegment\t15\t3\t11\t6\t8\t30\tcs:Z::1-G:1+T:2-CA:2");
    }

}


}
}
