/**
 * unittest/genotypekit.cpp: test cases for genotypekit modular genotyper pieces
 */

#include "catch.hpp"
#include "../bubbles.hpp"
#include "../vg.hpp"
#include "../json2pb.h"

namespace vg {
namespace unittest {

TEST_CASE("bubbles can be found", "[bubbles]") {
    
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
    graph.extend(chunk);
    
    // We need to see the path.
    REQUIRE(graph.paths.size() == 1);
    
    BubbleTree* bubble_tree = ultrabubble_tree(graph);
    
#ifdef debug
    bubble_tree->for_each_preorder([&](const TreeNode<Bubble>* node) {
        cerr << "Found bubble " << node->v.start << " to " << node->v.end << " containing ";
        for (auto& id : node->v.contents) {
            cerr << id << " ";
        }
        cerr << endl;
    });
#endif
    
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
                    
                    SECTION("Subchild has no children") {
                        REQUIRE(subchild->children.size() == 0);
                    }
                        
                }
                
            }
            
            SECTION("Second child is from 6 end to 9 start") {
                REQUIRE(child2->v.start.node == 6);
                REQUIRE(child2->v.start.is_end == true);
                REQUIRE(child2->v.end.node == 9);
                REQUIRE(child2->v.end.is_end == false);
                
                SECTION("Second child has no children") {
                    REQUIRE(child2->children.size() == 0);
                }
            }
            
        }
        
    }
    
    delete bubble_tree;
    
}

TEST_CASE("bubbles can be found in graphs with only heads", "[bubbles]") {
    
    // Build a toy graph
    const string graph_json = R"(
    
    {
        "node": [
            {"id": 1, "sequence": "G"},
            {"id": 2, "sequence": "A"}
        ],
        "edge": [
            {"from": 1, "to": 2, "to_end": true}
            
        ],
        "path": [
            {"name": "hint", "mapping": [
                {"position": {"node_id": 1}, "rank" : 1 },
                {"position": {"node_id": 2, "is_reverse": true}, "rank" : 2 }
            ]}
        ]
    }
    
    )";
    
    // Make an actual graph
    VG graph;
    Graph chunk;
    json2pb(chunk, graph_json.c_str(), graph_json.size());
    graph.extend(chunk);
    
    BubbleTree* bubble_tree = ultrabubble_tree(graph);

#ifdef debug
    bubble_tree->for_each_preorder([&](const TreeNode<Bubble>* node) {
        cerr << "Found bubble " << node->v.start << " to " << node->v.end << " containing ";
        for (auto& id : node->v.contents) {
            cerr << id << " ";
        }
        cerr << endl;
    });
#endif
    
    SECTION("Tree contains a root node that is a fake global bubble") {
        REQUIRE(bubble_tree->root != nullptr);
        REQUIRE(bubble_tree->root->v.start.node == 0);
        REQUIRE(bubble_tree->root->v.end.node == 0);
        
        SECTION("Root node has 1 child bubble") {
            REQUIRE(bubble_tree->root->children.size() == 1);
            
            TreeNode<Bubble>* child1 = bubble_tree->root->children[0];
            
            SECTION("First child is from 1 end to 2 end") {
                REQUIRE(child1->v.start.node == 1);
                REQUIRE(child1->v.start.is_end == true);
                REQUIRE(child1->v.end.node == 2);
                REQUIRE(child1->v.end.is_end == true);
                
                SECTION("First child has no children") {
                    REQUIRE(child1->children.size() == 0);
                }
                
            }
            
        }
        
    }
    
    delete bubble_tree;
    
}


TEST_CASE("bubbles can be found in bigger graphs with only heads", "[bubbles]") {
    
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
            {"id": 9, "sequence": "T"}
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
            {"from": 7, "to": 9, "to_end": true},
            {"from": 8, "to": 9, "to_end": true}
            
        ],
        "path": [
            {"name": "hint", "mapping": [
                {"position": {"node_id": 1}, "rank" : 1 },
                {"position": {"node_id": 6}, "rank" : 2 },
                {"position": {"node_id": 8}, "rank" : 3 },
                {"position": {"node_id": 9, "is_reverse": true}, "rank" : 4 }
            ]}
        ]
    }
    
    )";
    
    // Make an actual graph
    VG graph;
    Graph chunk;
    json2pb(chunk, graph_json.c_str(), graph_json.size());
    graph.extend(chunk);
    
    BubbleTree* bubble_tree = ultrabubble_tree(graph);
    
#ifdef debug
    bubble_tree->for_each_preorder([&](const TreeNode<Bubble>* node) {
        cerr << "Found bubble " << node->v.start << " to " << node->v.end << " containing ";
        for (auto& id : node->v.contents) {
            cerr << id << " ";
        }
        cerr << endl;
    });
#endif
    
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
                    
                    SECTION("Subchild has no children") {
                        REQUIRE(subchild->children.size() == 0);
                    }
                        
                }
                
            }
            
            SECTION("Second child is from 6 end to 9 end") {
                REQUIRE(child2->v.start.node == 6);
                REQUIRE(child2->v.start.is_end == true);
                REQUIRE(child2->v.end.node == 9);
                REQUIRE(child2->v.end.is_end == true);
                
                SECTION("Second child has no children") {
                    REQUIRE(child2->children.size() == 0);
                }
            }
            
        }
        
    }
    
    delete bubble_tree;
    
}

TEST_CASE("bubbles can be found in graphs with only tails", "[bubbles]") {
    
    // Build a toy graph
    const string graph_json = R"(
    
    {
        "node": [
            {"id": 1, "sequence": "C"},
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
            {"from": 1, "to": 2, "from_start": true},
            {"from": 1, "to": 6, "from_start": true},
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
                {"position": {"node_id": 1, "is_reverse": true}, "rank" : 1 },
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
    graph.extend(chunk);
    
    BubbleTree* bubble_tree = ultrabubble_tree(graph);
    
#ifdef debug
    bubble_tree->for_each_preorder([&](const TreeNode<Bubble>* node) {
        cerr << "Found bubble " << node->v.start << " to " << node->v.end << " containing ";
        for (auto& id : node->v.contents) {
            cerr << id << " ";
        }
        cerr << endl;
    });
#endif
    
    SECTION("Tree contains a root node that is a fake global bubble") {
        REQUIRE(bubble_tree->root != nullptr);
        REQUIRE(bubble_tree->root->v.start.node == 0);
        REQUIRE(bubble_tree->root->v.end.node == 0);
        
        SECTION("Root node has 2 child bubbles") {
            REQUIRE(bubble_tree->root->children.size() == 2);
            
            TreeNode<Bubble>* child1 = bubble_tree->root->children[0];
            TreeNode<Bubble>* child2 = bubble_tree->root->children[1];
            
            SECTION("First child is from 1 start to 6 start") {
                REQUIRE(child1->v.start.node == 1);
                REQUIRE(child1->v.start.is_end == false);
                REQUIRE(child1->v.end.node == 6);
                REQUIRE(child1->v.end.is_end == false);
                
                SECTION("First child should have all the contained nodes in its contents, including contents of its children") {
                    REQUIRE(child1->v.contents.size() == 6);
                };
                
                SECTION("First child has a child from 2 end to 5 start") {
                    REQUIRE(child1->children.size() == 1);
                    
                    TreeNode<Bubble>* subchild = child1->children[0];
                    
                    REQUIRE(subchild->v.start.node == 2);
                    REQUIRE(subchild->v.start.is_end == true);
                    REQUIRE(subchild->v.end.node == 5);
                    REQUIRE(subchild->v.end.is_end == false);
                    
                    SECTION("Subchild has no children") {
                        REQUIRE(subchild->children.size() == 0);
                    }
                        
                }
                
            }
            
            SECTION("Second child is from 6 end to 9 start") {
                REQUIRE(child2->v.start.node == 6);
                REQUIRE(child2->v.start.is_end == true);
                REQUIRE(child2->v.end.node == 9);
                REQUIRE(child2->v.end.is_end == false);
                
                SECTION("Second child has no children") {
                    REQUIRE(child2->children.size() == 0);
                }
            }
            
        }
        
    }
    
    delete bubble_tree;
    
}

TEST_CASE("bubbles can be found when heads cannot reach tails", "[bubbles]") {
    
    // Build a toy graph
    // Looks like:
    //
    //    1\
    //  2---3
    //  v\4 v
    //
    // The head is 1, the tail is 4, 2 and 3 have self loops, and the graph is connected.
    // But the head cannot reach the tail.
    const string graph_json = R"(
    
    {
        "node": [
            {"id": 1, "sequence": "A"},
            {"id": 2, "sequence": "A"},
            {"id": 3, "sequence": "A"},
            {"id": 4, "sequence": "A"}
        ],
        "edge": [
            {"from": 1, "to": 3},
            {"from": 2, "to": 4},
            {"from": 2, "to": 3},
            {"from": 2, "to": 2},
            {"from": 3, "to": 3}            
        ]
    }
    
    )";
    
    // Make an actual graph
    VG graph;
    Graph chunk;
    json2pb(chunk, graph_json.c_str(), graph_json.size());
    graph.extend(chunk);
    
    BubbleTree* bubble_tree = ultrabubble_tree(graph);
    
#ifdef debug
    bubble_tree->for_each_preorder([&](const TreeNode<Bubble>* node) {
        cerr << "Found bubble " << node->v.start << " to " << node->v.end << " containing ";
        for (auto& id : node->v.contents) {
            cerr << id << " ";
        }
        cerr << endl;
    });
#endif
    
    SECTION("Tree contains a root node that is a fake global bubble") {
        REQUIRE(bubble_tree->root != nullptr);
        REQUIRE(bubble_tree->root->v.start.node == 0);
        REQUIRE(bubble_tree->root->v.end.node == 0);
        
        SECTION("Root should have one child actual bubble") {
            REQUIRE(bubble_tree->root->children.size() == 1);
            
            TreeNode<Bubble>* child1 = bubble_tree->root->children[0];
            
            SECTION("The child should contain all the nodes as contents") {
                REQUIRE(child1->v.contents.size() == 4);
            }
            
            // TODO: When unary snarls are exposed, make sure we have unary children.
        }
    }
    
    delete bubble_tree;
    
}

TEST_CASE("bubbles can be found in a graph with no heads or tails", "[bubbles]") {
    
    // Build a toy graph
    const string graph_json = R"(
    
    {
        "node": [
            {"id": 1, "sequence": "G"},
            {"id": 2, "sequence": "A"}
        ],
        "edge": [
            {"from": 1, "to": 2},
            {"from": 2, "to": 1}
            
        ],
        "path": [
            {"name": "hint", "mapping": [
                {"position": {"node_id": 1}, "rank" : 1 },
                {"position": {"node_id": 2}, "rank" : 2 }
            ]}
        ]
    }
    
    )";
    
    // Make an actual graph
    VG graph;
    Graph chunk;
    json2pb(chunk, graph_json.c_str(), graph_json.size());
    graph.extend(chunk);
    
    BubbleTree* bubble_tree = ultrabubble_tree(graph);
    
#ifdef debug
    bubble_tree->for_each_preorder([&](const TreeNode<Bubble>* node) {
        cerr << "Found bubble " << node->v.start << " to " << node->v.end << " containing ";
        for (auto& id : node->v.contents) {
            cerr << id << " ";
        }
        cerr << endl;
    });
#endif
    
    SECTION("Tree contains a root node that is a fake global bubble") {
        REQUIRE(bubble_tree->root != nullptr);
        REQUIRE(bubble_tree->root->v.start.node == 0);
        REQUIRE(bubble_tree->root->v.end.node == 0);
        
        SECTION("Root node has 2 child bubbles") {
            REQUIRE(bubble_tree->root->children.size() == 2);
            
            TreeNode<Bubble>* child1 = bubble_tree->root->children[0];
            TreeNode<Bubble>* child2 = bubble_tree->root->children[1];
            
            SECTION("First child is from 1 start/end to 2 end/start") {
                REQUIRE(child1->v.start.node == 1);
                REQUIRE(child1->v.end.node == 2);
                REQUIRE(child1->v.start.is_end != child1->v.end.is_end);
                
                SECTION("First child has no children") {
                    REQUIRE(child1->children.size() == 0);
                }
                
            }
            
            SECTION("Second child is from 1 start/end to 2 end/start") {
                REQUIRE(child2->v.start.node == 1);
                REQUIRE(child2->v.end.node == 2);
                REQUIRE(child1->v.start.is_end != child1->v.end.is_end);
                
                SECTION("Second child has no children") {
                    REQUIRE(child2->children.size() == 0);
                }
                
            }
            
        }
        
    }
    
    delete bubble_tree;
    
}

TEST_CASE("bubbles can be found in a graph with no heads or tails or paths", "[bubbles]") {
    
    // Build a toy graph
    const string graph_json = R"(
    
    {
        "node": [
            {"id": 1, "sequence": "G"},
            {"id": 2, "sequence": "A"}
        ],
        "edge": [
            {"from": 1, "to": 2},
            {"from": 2, "to": 1}
            
        ]
    }
    
    )";
    
    // Make an actual graph
    VG graph;
    Graph chunk;
    json2pb(chunk, graph_json.c_str(), graph_json.size());
    graph.extend(chunk);
    
    BubbleTree* bubble_tree = ultrabubble_tree(graph);
    
#ifdef debug
    bubble_tree->for_each_preorder([&](const TreeNode<Bubble>* node) {
        cerr << "Found bubble " << node->v.start << " to " << node->v.end << " containing ";
        for (auto& id : node->v.contents) {
            cerr << id << " ";
        }
        cerr << endl;
    });
#endif
    
    SECTION("Tree contains a root node that is a fake global bubble") {
        REQUIRE(bubble_tree->root != nullptr);
        REQUIRE(bubble_tree->root->v.start.node == 0);
        REQUIRE(bubble_tree->root->v.end.node == 0);
        
        SECTION("Root node has 2 child bubbles") {
            REQUIRE(bubble_tree->root->children.size() == 2);
            
            // We should have 2 bubbles, one looking around the cycle in each direction.
            // TODO: can't really say much about its contents.
            
        }
        
    }
    
    delete bubble_tree;
    
}

TEST_CASE("bubbles are created based on most distant connected tips", "[bubbles]") {
    
    // Build a toy graph
    // Looks like:
    //
    //       1\
    //  5--2---3--6
    //      \4  
    //
    // The head is 1, the tail is 4, 2 and 3 have self loops, and the graph is connected.
    // But the head cannot reach the tail.
    const string graph_json = R"(
    
    {
        "node": [
            {"id": 1, "sequence": "A"},
            {"id": 2, "sequence": "A"},
            {"id": 3, "sequence": "A"},
            {"id": 4, "sequence": "A"},
            {"id": 5, "sequence": "A"},
            {"id": 6, "sequence": "A"}
        ],
        "edge": [
            {"from": 1, "to": 3},
            {"from": 2, "to": 4},
            {"from": 2, "to": 3},
            {"from": 5, "to": 2},
            {"from": 3, "to": 6}
        ]
    }
    
    )";
    
    // Make an actual graph
    VG graph;
    Graph chunk;
    json2pb(chunk, graph_json.c_str(), graph_json.size());
    graph.extend(chunk);
    
    BubbleTree* bubble_tree = ultrabubble_tree(graph);
    
#ifdef debug
    bubble_tree->for_each_preorder([&](const TreeNode<Bubble>* node) {
        cerr << "Found bubble " << node->v.start << " to " << node->v.end << " containing ";
        for (auto& id : node->v.contents) {
            cerr << id << " ";
        }
        cerr << endl;
    });
#endif
    
    SECTION("Tree contains a root node that is a fake global bubble") {
        REQUIRE(bubble_tree->root != nullptr);
        REQUIRE(bubble_tree->root->v.start.node == 0);
        REQUIRE(bubble_tree->root->v.end.node == 0);
        
        SECTION("Root should have three child actual bubbles") {
            REQUIRE(bubble_tree->root->children.size() == 3);
            
            TreeNode<Bubble>* child1 = bubble_tree->root->children[0];
            TreeNode<Bubble>* child2 = bubble_tree->root->children[1];
            TreeNode<Bubble>* child3 = bubble_tree->root->children[2];
            
            SECTION("First child is trivail snarl from 2 start to 5 end") {
                // Not 5 end to 2 start because we seem to be sorting by node ID.
                REQUIRE(child1->v.start.node == 2);
                REQUIRE(child1->v.start.is_end == false);
                REQUIRE(child1->v.end.node == 5);
                REQUIRE(child1->v.end.is_end == true);
                
            }
            
            SECTION("Second child is tip-containing snarl from 2 end to 3 start") {
                REQUIRE(child2->v.start.node == 2);
                REQUIRE(child2->v.start.is_end == true);
                REQUIRE(child2->v.end.node == 3);
                REQUIRE(child2->v.end.is_end == false);
                
            }
            
            SECTION("Third child is trivial snarl from 3 end to 6 start") {
                REQUIRE(child3->v.start.node == 3);
                REQUIRE(child3->v.start.is_end == true);
                REQUIRE(child3->v.end.node == 6);
                REQUIRE(child3->v.end.is_end == false);
                
            }
            
        }
        
        
        
    }
    
    delete bubble_tree;
    
}


}
}
