/**
 * unittest/genotypekit.cpp: test cases for genotypekit modular genotyper pieces
 */

#include "catch.hpp"
#include "../bubbles.hpp"
#include "../vg.hpp"
#include "../json2pb.h"

namespace vg {
namespace unittest {

TEST_CASE("bubbles can be found with Cactus", "[bubbles]") {
    
    // Build a toy graph
    const string graph_json = R"(
    
    {
        "node": [
            {"id": 1, "sequence": "G"},
            {"id": 2, "sequence": "A"},
            {"id": 3, "sequence": "T"},
            {"id": 4, "sequence": "GGG"},
            {"id": 5, "sequence": "T"},
            {"id": 6, "sequence": "A"},
            {"id": 7, "sequence": "C"},
            {"id": 8, "sequence": "A"},
            {"id": 9, "sequence": "A"}
        ],
        "edge": [
            {"from": 1, "to": 2},
            {"from": 1, "to": 6},
            {"from": 2, "to": 3},
            {"from": 2, "to": 4},
            {"from": 3, "to": 5},
            {"from": 4, "to": 5},
            {"from": 5, "to": 6},
            {"from": 6, "to": 7},
            {"from": 6, "to": 8},
            {"from": 7, "to": 9},
            {"from": 8, "to": 9}
            
        ],
        "path": [
            {"name": "hint", "mapping": [
                {"position": {"node_id": 1}, "rank" : 1 },
                {"position": {"node_id": 6}, "rank" : 2 },
                {"position": {"node_id": 8}, "rank" : 3 },
                {"position": {"node_id": 9}, "rank" : 4 }
            ]}
        ]
    }
    
    )";
    
    // Make an actual graph
    VG graph;
    Graph chunk;
    json2pb(chunk, graph_json.c_str(), graph_json.size());
    graph.merge(chunk);
    
    BubbleTree* bubble_tree = ultrabubble_tree(graph);
    
    bubble_tree->for_each_preorder([&](const TreeNode<Bubble>* node) {
        cerr << "Found bubble " << node->v.start << " to " << node->v.end << " containing ";
        for (auto& id : node->v.contents) {
            cerr << id << " ";
        }
        cerr << endl;
    });
    
    SECTION("Tree contains a root node that is a fake global bubble") {
        REQUIRE(bubble_tree->root != nullptr);
        REQUIRE(bubble_tree->root->v.start.node == 0);
        REQUIRE(bubble_tree->root->v.end.node == 0);
        
        SECTION("Root node has 2 child bubbles") {
            REQUIRE(bubble_tree->root->children.size() == 2);
            
            TreeNode<Bubble>* child1 = bubble_tree->root->children[0];
            TreeNode<Bubble>* child2 = bubble_tree->root->children[1];
            
            SECTION("First child is from 1 end to 6 start") {
                REQUIRE(child1->v.start.node == 1);
                REQUIRE(child1->v.start.is_end == true);
                REQUIRE(child1->v.end.node == 6);
                REQUIRE(child1->v.end.is_end == false);
                
                SECTION("First child has a child from 2 end to 5 start") {
                    REQUIRE(child1->children.size() == 1);
                    
                    TreeNode<Bubble>* subchild = child1->children[0];
                    
                    REQUIRE(subchild->v.start.node == 2);
                    REQUIRE(subchild->v.start.is_end == true);
                    REQUIRE(subchild->v.end.node == 5);
                    REQUIRE(subchild->v.end.is_end == false);
                    
                    SECTION("Subchild child has no children") {
                        REQUIRE(subchild->children.size() == 0);
                    }
                        
                }
                
            }
            
            SECTION("Second child is from 6 end to 9 start") {
                REQUIRE(child2->v.start.node == 6);
                REQUIRE(child2->v.start.is_end == true);
                REQUIRE(child2->v.end.node == 9);
                REQUIRE(child2->v.end.is_end == false);
                
                SECTION("First child has no children") {
                    REQUIRE(child2->children.size() == 0);
                }
            }
            
        }
        
    }
    
    delete bubble_tree;
    
}

}
}
