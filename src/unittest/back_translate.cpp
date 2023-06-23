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
#include "alignment.hpp"

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

TEST_CASE("A Path can be back-translated properly when it goes around cycles", "[algorithms][back_translate]") {

    // Make a path
    string path_string = R"(
        {
            "mapping": [
                {
                    "position": {"node_id": 1},
                    "edit": [
                        {"from_length": 10, "to_length": 10}
                    ]
                },
                {
                    "position": {"node_id": 1},
                    "edit": [
                        {"from_length": 10, "to_length": 10}
                    ]
                },
                {
                    "position": {"node_id": 1},
                    "edit": [
                        {"from_length": 10, "to_length": 10}
                    ]
                },
                {
                    "position": {"node_id": 1, "is_reverse": true},
                    "edit": [
                        {"from_length": 10, "to_length": 10}
                    ]
                },
                {
                    "position": {"node_id": 2},
                    "edit": [
                        {"from_length": 5, "to_length": 5}
                    ]
                },
                {
                    "position": {"node_id": 3},
                    "edit": [
                        {"from_length": 5, "to_length": 5}
                    ]
                },
                {
                    "position": {"node_id": 2},
                    "edit": [
                        {"from_length": 5, "to_length": 5}
                    ]
                },
                {
                    "position": {"node_id": 3},
                    "edit": [
                        {"from_length": 5, "to_length": 5}
                    ]
                },
                {
                    "position": {"node_id": 2},
                    "edit": [
                        {"from_length": 5, "to_length": 5}
                    ]
                }
            ]
        }
    )";
    Path p;
    json2pb(p, path_string.c_str(), path_string.size());
    
    // Define a translation back to a named node space.
    // We have node 1 just be as is and nodes 2 and 3 come from two halves of segment 2
    MockBackTranslation trans;
    trans.translation[oriented_node_range_t(1, false, 0, 10)] = {oriented_node_range_t(1, false, 0, 10)};
    trans.translation[oriented_node_range_t(1, true, 0, 10)] = {oriented_node_range_t(1, true, 0, 10)};
    trans.translation[oriented_node_range_t(2, false, 0, 5)] = {oriented_node_range_t(2, false, 0, 5)};
    trans.translation[oriented_node_range_t(3, false, 0, 5)] = {oriented_node_range_t(2, false, 5, 5)};
    trans.node_names[1] = "A";
    trans.node_names[2] = "B";
    
    // Translate the path
    vg::algorithms::back_translate_in_place(&trans, p);

    // Check the result
    REQUIRE(p.mapping_size() == 7);
    REQUIRE(p.mapping(0).edit_size() == 1);
    REQUIRE(p.mapping(0).edit(0).from_length() == 10);
    REQUIRE(p.mapping(0).position().name() == "A");
    REQUIRE(p.mapping(0).position().offset() == 0);
    REQUIRE(p.mapping(0).position().is_reverse() == false);
    
    REQUIRE(p.mapping(1).edit_size() == 1);
    REQUIRE(p.mapping(1).edit(0).from_length() == 10);
    REQUIRE(p.mapping(1).position().name() == "A");
    REQUIRE(p.mapping(1).position().offset() == 0);
    REQUIRE(p.mapping(1).position().is_reverse() == false);
    
    REQUIRE(p.mapping(2).edit_size() == 1);
    REQUIRE(p.mapping(2).edit(0).from_length() == 10);
    REQUIRE(p.mapping(2).position().name() == "A");
    REQUIRE(p.mapping(2).position().offset() == 0);
    REQUIRE(p.mapping(2).position().is_reverse() == false);
    
    REQUIRE(p.mapping(3).edit_size() == 1);
    REQUIRE(p.mapping(3).edit(0).from_length() == 10);
    REQUIRE(p.mapping(3).position().name() == "A");
    REQUIRE(p.mapping(3).position().offset() == 0);
    REQUIRE(p.mapping(3).position().is_reverse() == true);
    
    REQUIRE(p.mapping(4).edit_size() == 1);
    REQUIRE(p.mapping(4).edit(0).from_length() == 10);
    REQUIRE(p.mapping(4).position().name() == "B");
    REQUIRE(p.mapping(4).position().offset() == 0);
    REQUIRE(p.mapping(4).position().is_reverse() == false);
    
    REQUIRE(p.mapping(5).edit_size() == 1);
    REQUIRE(p.mapping(5).edit(0).from_length() == 10);
    REQUIRE(p.mapping(5).position().name() == "B");
    REQUIRE(p.mapping(5).position().offset() == 0);
    REQUIRE(p.mapping(5).position().is_reverse() == false);
    
    REQUIRE(p.mapping(6).edit_size() == 1);
    REQUIRE(p.mapping(6).edit(0).from_length() == 5);
    REQUIRE(p.mapping(6).position().name() == "B");
    REQUIRE(p.mapping(6).position().offset() == 0);
    REQUIRE(p.mapping(6).position().is_reverse() == false);
}

TEST_CASE("An Alignment can be back-translated while converting to GAF", "[algorithms][back_translate]") {

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

    REQUIRE(alignment_from_length(a) == 9);
    REQUIRE(alignment_to_length(a) == 7);
    size_t block_length = std::max(alignment_from_length(a), alignment_to_length(a));

    SECTION("Translating GAF generation produces the right GAF") {
        auto node_space = vg::io::alignment_to_gaf(g, a);
        stringstream s;
        s << node_space;
        
        // See column definitions at https://github.com/lh3/gfatools/blob/master/doc/rGFA.md#the-graph-alignment-format-gaf
        // Alignment block length is longest involved sequence.
        // Note that we combine adjacent duplicate operations in cs across node boundaries.
        // Note that end position is 0-based inclusive
        // Path length is 13 bp, and we run [2 to 11) (9 bases)
        REQUIRE(s.str() == "francine\t7\t0\t7\t+\t>2>3>4\t13\t2\t11\t6\t" + std::to_string(block_length) + "\t30\tcs:Z::1-G:1+T:2-CA:2");
    }
    
    SECTION("Translating GAF generation with a translation produces the right GAF") {
        auto back_translated = vg::io::alignment_to_gaf(g, a, &trans);
        stringstream s;
        s << back_translated;
        
        REQUIRE(s.str() == "francine\t7\t0\t7\t+\t>FirstSegment>SecondSegment\t15\t3\t12\t6\t" + std::to_string(block_length) + "\t30\tcs:Z::1-G:1+T:2-CA:2");
    }
}

TEST_CASE("A reverse-strand Alignment can be back-translated while converting to GAF", "[algorithms][back_translate]") {

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
    trans.translation[oriented_node_range_t(3, true, 1, 1)] = {oriented_node_range_t(1, true, 1, 1)};
    trans.translation[oriented_node_range_t(2, true, 0, 2)] = {oriented_node_range_t(1, true, 2, 2)};
    // Make sure to include reverse strand translation info for segment length sniffing.
    trans.translation[oriented_node_range_t(2, false, 0, 0)] = {oriented_node_range_t(1, false, 1, 0)};
    trans.translation[oriented_node_range_t(3, false, 0, 0)] = {oriented_node_range_t(1, false, 5, 0)};
    trans.node_names[1] = "FirstSegment";
    trans.node_names[2] = "SecondSegment";


    string alignment_string = R"(
        {
            "name": "steve",
            "mapping_quality": 30,
            "sequence": "TCC",
            "path": {"mapping": [
                {
                    "position": {"node_id": 3, "offset": 1, "is_reverse": true},
                    "edit": [
                        {"from_length": 1, "to_length": 1}
                    ]
                },
                {
                    "position": {"node_id": 2, "is_reverse": true},
                    "edit": [
                        {"from_length": 2, "to_length": 2}
                    ]
                }
            ]}
        }
    )";
    
    // We take up 3 bases.
    // We leave 1 base of node 3 on the left, and 2 bases of node 4 on the right.
    // We leave 1 base of segment 1 on the left and 3 bases of segment 1 on the right.

    Alignment a;
    json2pb(a, alignment_string.c_str(), alignment_string.size());

    REQUIRE(alignment_from_length(a) == 3);
    REQUIRE(alignment_to_length(a) == 3);
    size_t block_length = std::max(alignment_from_length(a), alignment_to_length(a));
    
    auto back_translated = vg::io::alignment_to_gaf(g, a, &trans);
    stringstream s;
    s << back_translated;
    
    // Should be mapped to a length-7 path, at offset 1 to 4, 3 matches in a length 3 block.
    REQUIRE(s.str() == "steve\t3\t0\t3\t+\t<FirstSegment\t7\t1\t4\t3\t" + std::to_string(block_length) + "\t30\tcs:Z::3");
}

    

}
}
