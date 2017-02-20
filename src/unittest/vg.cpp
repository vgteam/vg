/**
 * unittest/vg.cpp: test cases for vg::VG methods
 */

#include "catch.hpp"
#include "vg.hpp"

namespace vg {
namespace unittest {

using namespace std;

// Turn a JSON string into a VG graph
VG string_to_graph(const string& json) {
    VG graph;
    Graph chunk;
    json2pb(chunk, json.c_str(), json.size());
    graph.merge(chunk);
    
    return graph;
}

TEST_CASE("is_acyclic() should return whether the graph is acyclic", "[vg][cycles]") {
    
    SECTION("a tiny DAG should be acyclic") {
        const string graph_json = R"(
        
        {
            "node": [
                {"id": 1, "sequence": "G"},
                {"id": 2, "sequence": "A"}
            ],
            "edge": [
                {"from": 1, "to": 2}
            ]
        }
    
        )";
        
        VG graph = string_to_graph(graph_json);
        
        REQUIRE(graph.is_acyclic() == true);
    }
    
    SECTION("a tiny cyclic graph should be cyclic") {
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
        
        VG graph = string_to_graph(graph_json);
        
        REQUIRE(graph.is_acyclic() == false);
    }
    
    SECTION("a tiny cyclic graph using from_start and to_end should be cyclic") {
        const string graph_json = R"(
        
        {
            "node": [
                {"id": 1, "sequence": "G"},
                {"id": 2, "sequence": "A"}
            ],
            "edge": [
                {"from": 1, "to": 2},
                {"from": 1, "to": 2, "from_start": true, "to_end": true}
            ]
        }
    
        )";
        
        VG graph = string_to_graph(graph_json);
        
        REQUIRE(graph.is_acyclic() == false);
    }
    
    SECTION("a tiny cyclic graph using from_start and to_end the other way should be cyclic") {
        const string graph_json = R"(
        
        {
            "node": [
                {"id": 1, "sequence": "G"},
                {"id": 2, "sequence": "A"}
            ],
            "edge": [
                {"from": 2, "to": 1},
                {"from": 2, "to": 1, "from_start": true, "to_end": true}
            ]
        }
    
        )";
        
        VG graph = string_to_graph(graph_json);
        
        REQUIRE(graph.is_acyclic() == false);
    }
    
    SECTION("a nontrivial DAG should be acyclic") {
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
                
            ]
        }
    
        )";
        
        VG graph = string_to_graph(graph_json);
        
        REQUIRE(graph.is_acyclic() == true);
    }
}

TEST_CASE("unfold() should properly unfold a graph out to the requested length", "[vg][unfold]") {

    SECTION("unfolding across a single reversing edge should double the entire graph") {

        const string graph_json = R"(
        {
          "node": [
            {"sequence": "CACACACACACTGAGATACCATCTCCCGTCGGTCAGAATGGCCATTTGTCCAAAGT","id": 1},
            {"sequence": "CCAACATTTACATGCTGACAAGGCAGCGGAGGAAAGCAAACACTG","id": 2},
            {"sequence": "C","id": 3},
            {"sequence": "T","id": 4},
            {"sequence": "TACACTCTTGGAGGGAA","id": 5},
            {"sequence": "T","id": 6},
            {"sequence": "C","id": 7},
            {"sequence": "AAAAACTAG","id": 8},
            {"sequence": "AGTTGCAT","id": 9},
            {"sequence": "TTCTCTGATGATGAG","id": 10},
            {"sequence": "TGATGTTGAGGGTTTTTTTTGTCT","id": 11},
            {"sequence": "ATTGGTCACTTGTACATCTTATTTTTACAA","id": 12},
            {"sequence":"GAACGTCTTCATGTCTTTTACCCATTGTTTGATTGGCTTATTTGTCTTTCTGCCTCTTGATTTTTTTAAGGTCACGTTAGAGTCTGGCTATTAGGTCTTGGTCAGAAGCATAGTTTGGGAACATTTTCTCCCATCCTGTAGGCTATGTTTTTACTGTGCTGGTGATTTATTTTGCTGTCTGGCAGTGTTTTAGTTTCTTAGGCCCAACTTGTCCATTTTGGTTTTTCTTGCCGTTGCTCTTAGGGACTAAGTTATTTAAATTCTTTACCAAAGCCCATGTTGAGAAAGGTATTTCCTAGTTTTTCTTGTAGGACTTGTATAGTTTGTAGTCTTCTGCTGAAATCTTTCATTCACCTTGAGTTAGTTTTTGCATCTTGTGAGAGGTAAGGCTGTAGTGCTGTTCATCTGCAAGTGGCTAGACACTATCCCAGTGCCATTCATTGCACAGTGAGCCCTTTCCACATCTGAATTTTGGTGTTTCAAAGGTCAGATGGTTGTGGGTATGTGGGGTTGCTTCTGGGTTTCCTATTCTGTCTAGGTGGAGGTGGGTCTGTAGCTTTTCTGTGAAATTTGAAATTGGAGAGTGTGGTCCATCTGACATCGTCTGAATCTCCCAGCCTGGCTTTGGAGTTTCAGAGCATTTTGTGGTCCCATGGGAATTTTAGCATTCATTGTTTCTTCACATTGCTTTCAAAAAACAAAATCGACCCCATCCTAAAGGTGTACAGGTAGGGGGTAGAGTGGAATATTTGGATCCATGC","id": 13}
          ],
          "edge": [
            {"from": 1,"to": 9,"from_start": true},
            {"from": 1,"to": 2},
            {"from": 2,"to": 3},
            {"from": 2,"to": 4},
            {"from": 3, "to": 5},
            {"from": 4,"to": 5},
            {"from": 5,"to": 6},
            {"from": 5,"to": 7},
            {"from": 6,"to": 8},
            {"from": 7,"to": 8},
            {"from": 9,"to": 10},
            {"from": 10,"to": 11},
            {"from": 11,"to": 12},
            {"from": 12,"to": 13}
          ]
        }
        )";
        
        VG graph = string_to_graph(graph_json);
        
        map<id_t, pair<id_t, bool> > node_translation;
        VG completely_unfolded = graph.unfold(10000, node_translation);
        
        REQUIRE(completely_unfolded.size() == graph.size() * 2);
    }

}

TEST_CASE("expand_context_by_length() should respect barriers", "[vg][context]") {

    const string graph_json = R"(
    {
      "node": [
        {"sequence": "CCATTTGTCCAAAGT","id": 1},
        {"sequence": "AAGCAAACACTG","id": 2},
        {"sequence": "C","id": 3},
        {"sequence": "T","id": 4},
        {"sequence": "TACACTCTTGGAGGGAA","id": 5},
        {"sequence": "T","id": 6},
        {"sequence": "C","id": 7},
        {"sequence": "AAAAACTAG","id": 8},
        {"sequence": "AGTTGCAT","id": 9},
        {"sequence": "TTCTCTGATGATGAG","id": 10},
        {"sequence": "TGATGTTGAGGGTTTTTTTTGTCT","id": 11},
        {"sequence": "ATTGGTCACTTGTACATCTTATTTTTACAA","id": 12},
        {"sequence":"GAACGTTT", "id": 13}
      ],
      "edge": [
        {"from": 1,"to": 9,"from_start": true},
        {"from": 1,"to": 2},
        {"from": 2,"to": 3},
        {"from": 2,"to": 4},
        {"from": 3, "to": 5},
        {"from": 4,"to": 5},
        {"from": 5,"to": 6},
        {"from": 5,"to": 7},
        {"from": 6,"to": 8},
        {"from": 7,"to": 8},
        {"from": 9,"to": 10},
        {"from": 10,"to": 11},
        {"from": 11,"to": 12},
        {"from": 12,"to": 13}
      ]
    }
    )";
    
    VG graph = string_to_graph(graph_json);

    SECTION("barriers on either end of the seed node should stop anything being extracted") {

        VG context;
        context.add_node(*graph.get_node(3));
        graph.expand_context_by_length(context, 1000, false, true, {NodeSide(3, false), NodeSide(3, true)});
        
        REQUIRE(context.size() == 1);
    }
    
    SECTION("barriers should stop edges being formed") {

        VG context;
        context.add_node(*graph.get_node(3));
        context.add_node(*graph.get_node(4));
        // Note that we wouldn't get any edges between 3 and 4, if there were
        // any, because context expansion sees no edges between seed nodes.
        graph.expand_context_by_length(context, 1000, false, true, {NodeSide(3, false), NodeSide(3, true)});
        
        SECTION("node 4 should have both attached edges") {
            REQUIRE(context.has_edge(NodeSide(4, false), NodeSide(2, true)) == true);
            REQUIRE(context.has_edge(NodeSide(4, true), NodeSide(5, false)) == true);
        }
        
        SECTION("node 3 should have no atatched edges") {
            REQUIRE(context.has_edge(NodeSide(3, false), NodeSide(2, true)) == false);
            REQUIRE(context.has_edge(NodeSide(3, true), NodeSide(5, false)) == false);
        }
        
        
    }

}

TEST_CASE("bluntify() should resolve overlaps", "[vg][bluntify]") {
    
    SECTION("an overlap across a normal edge should be resolved") {
        const string graph_json = R"(
        
        {
            "node": [
                {"id": 1, "sequence": "GAA"},
                {"id": 2, "sequence": "AAT"}
            ],
            "edge": [
                {"from": 1, "to": 2, "overlap": 2}
            ]
        }
    
        )";
        
        VG graph = string_to_graph(graph_json);
        graph.bluntify();
        
        SECTION("the bluntified graph should have 3 nodes") {
            REQUIRE(graph.node_count() == 3);
            
            SECTION("their sequences should be G, AA, and T") {
                Node* g_node = nullptr;
                Node* aa_node = nullptr;
                Node* t_node = nullptr;
                
                graph.for_each_node([&](Node* n) {
                    if (n->sequence() == "G") {
                        g_node = n;
                    } else if (n->sequence() == "AA") {
                        aa_node = n;
                    } else if (n->sequence() == "T") {
                        t_node = n;
                    }
                });
                
                REQUIRE(g_node != nullptr);
                REQUIRE(aa_node != nullptr);
                REQUIRE(t_node != nullptr);
                
                SECTION("the right side of the G node should be connected to the left side of the AA node") {
                    REQUIRE(graph.has_edge(NodeSide(g_node->id(), true), NodeSide(aa_node->id(), false)));
                }
                
                SECTION("the right side of the AA node should be connected to the left side of the T node") {
                    REQUIRE(graph.has_edge(NodeSide(aa_node->id(), true), NodeSide(t_node->id(), false)));
                }
            }
        }
        
        SECTION("the bluntified graph should have 2 edges") {
            REQUIRE(graph.edge_count() == 2);
            
            SECTION("no edge should have overlap") {
                graph.for_each_edge([&](Edge* e) {
                    REQUIRE(e->overlap() == 0);
                });
            }
        }
        
    }
    
    SECTION("an overlap across a doubly-reversing edge should be resolved") {
        const string graph_json = R"(
        
        {
            "node": [
                {"id": 1, "sequence": "TTC"},
                {"id": 2, "sequence": "ATT"}
            ],
            "edge": [
                {"from": 1, "to": 2, "overlap": 2, "from_start": true, "to_end": true}
            ]
        }
    
        )";
        
        VG graph = string_to_graph(graph_json);
        graph.bluntify();
        
        SECTION("the bluntified graph should have 3 nodes") {
            REQUIRE(graph.node_count() == 3);
            
            SECTION("their sequences should be C, TT, and A") {
                Node* c_node = nullptr;
                Node* tt_node = nullptr;
                Node* a_node = nullptr;
                
                graph.for_each_node([&](Node* n) {
                    if (n->sequence() == "C") {
                        c_node = n;
                    } else if (n->sequence() == "TT") {
                        tt_node = n;
                    } else if (n->sequence() == "A") {
                        a_node = n;
                    }
                });
                
                REQUIRE(c_node != nullptr);
                REQUIRE(tt_node != nullptr);
                REQUIRE(a_node != nullptr);
                
                SECTION("the right side of the TT node should be connected to the left side of the C node") {
                    REQUIRE(graph.has_edge(NodeSide(c_node->id(), false), NodeSide(tt_node->id(), true)));
                }
                
                SECTION("the right side of the A node should be connected to the left side of the TT node") {
                    REQUIRE(graph.has_edge(NodeSide(tt_node->id(), false), NodeSide(a_node->id(), true)));
                }
            }
        }
        
        SECTION("the bluntified graph should have 2 edges") {
            REQUIRE(graph.edge_count() == 2);
            
            SECTION("no edge should have overlap") {
                graph.for_each_edge([&](Edge* e) {
                    REQUIRE(e->overlap() == 0);
                });
            }
        }
        
    }
    
    SECTION("an overlap across a reversing edge should be resolved") {
        const string graph_json = R"(
        
        {
            "node": [
                {"id": 1, "sequence": "TTC"},
                {"id": 2, "sequence": "AAT"}
            ],
            "edge": [
                {"from": 1, "to": 2, "overlap": 2, "from_start": true}
            ]
        }
    
        )";
        
        VG graph = string_to_graph(graph_json);
        graph.bluntify();
        
        SECTION("the bluntified graph should have 3 nodes") {
            REQUIRE(graph.node_count() == 3);
            
            SECTION("their sequences should be C, TT or AA, and T") {
                Node* c_node = nullptr;
                Node* middle_node = nullptr;
                bool is_tt;
                Node* t_node = nullptr;
                
                graph.for_each_node([&](Node* n) {
                    if (n->sequence() == "C") {
                        c_node = n;
                    } else if (n->sequence() == "TT") {
                        middle_node = n;
                        is_tt = true;
                    } else if (n->sequence() == "AA") {
                        middle_node = n;
                        is_tt = false;
                    } else if (n->sequence() == "T") {
                        t_node = n;
                    }
                });
                
                REQUIRE(c_node != nullptr);
                REQUIRE(middle_node != nullptr);
                REQUIRE(t_node != nullptr);
                
                SECTION("the left side of the C node should be connected to the right/left side of the TT/AA node") {
                    REQUIRE(graph.has_edge(NodeSide(c_node->id(), false), NodeSide(middle_node->id(), is_tt)));
                }
                
                SECTION("the left/right side of the TT/AA node should be connected to the left side of the T node") {
                    REQUIRE(graph.has_edge(NodeSide(middle_node->id(), !is_tt), NodeSide(t_node->id(), false)));
                }
            }
        }
        
        SECTION("the bluntified graph should have 2 edges") {
            REQUIRE(graph.edge_count() == 2);
            
            SECTION("no edge should have overlap") {
                graph.for_each_edge([&](Edge* e) {
                    REQUIRE(e->overlap() == 0);
                });
            }
        }
        
    }
    
    SECTION("extraneous overlap should be pruned back") {
        const string graph_json = R"(
        
        {
            "node": [
                {"id": 1, "sequence": "GAA"},
                {"id": 2, "sequence": "AAT"}
            ],
            "edge": [
                {"from": 1, "to": 2, "overlap": 3}
            ]
        }
    
        )";
        
        VG graph = string_to_graph(graph_json);
        graph.bluntify();
        
        SECTION("the bluntified graph should have 3 nodes") {
            REQUIRE(graph.node_count() == 3);
            
            SECTION("their sequences should be G, AA, and T") {
                Node* g_node = nullptr;
                Node* aa_node = nullptr;
                Node* t_node = nullptr;
                
                graph.for_each_node([&](Node* n) {
                    if (n->sequence() == "G") {
                        g_node = n;
                    } else if (n->sequence() == "AA") {
                        aa_node = n;
                    } else if (n->sequence() == "T") {
                        t_node = n;
                    }
                });
                
                REQUIRE(g_node != nullptr);
                REQUIRE(aa_node != nullptr);
                REQUIRE(t_node != nullptr);
                
                SECTION("the right side of the G node should be connected to the left side of the AA node") {
                    REQUIRE(graph.has_edge(NodeSide(g_node->id(), true), NodeSide(aa_node->id(), false)));
                }
                
                SECTION("the right side of the AA node should be connected to the left side of the T node") {
                    REQUIRE(graph.has_edge(NodeSide(aa_node->id(), true), NodeSide(t_node->id(), false)));
                }
            }
        }
        
        SECTION("the bluntified graph should have 2 edges") {
            REQUIRE(graph.edge_count() == 2);
            
            SECTION("no edge should have overlap") {
                graph.for_each_edge([&](Edge* e) {
                    REQUIRE(e->overlap() == 0);
                });
            }
        }
        
    }
    
    SECTION("overlaps should be able to overlap in the middle") {
        const string graph_json = R"(
        
        {
            "node": [
                {"id": 1, "sequence": "GAA"},
                {"id": 2, "sequence": "AA"},
                {"id": 3, "sequence": "AAT"}
            ],
            "edge": [
                {"from": 1, "to": 2, "overlap": 2},
                {"from": 2, "to": 3, "overlap": 2}
            ]
        }
    
        )";
        
        VG graph = string_to_graph(graph_json);
        graph.bluntify();
        
        SECTION("the bluntified graph should have 3 nodes") {
            REQUIRE(graph.node_count() == 3);
            
            SECTION("their sequences should be G, AA, and T") {
                Node* g_node = nullptr;
                Node* aa_node = nullptr;
                Node* t_node = nullptr;
                
                graph.for_each_node([&](Node* n) {
                    if (n->sequence() == "G") {
                        g_node = n;
                    } else if (n->sequence() == "AA") {
                        aa_node = n;
                    } else if (n->sequence() == "T") {
                        t_node = n;
                    }
                });
                
                REQUIRE(g_node != nullptr);
                REQUIRE(aa_node != nullptr);
                REQUIRE(t_node != nullptr);
                
                SECTION("the right side of the G node should be connected to the left side of the AA node") {
                    REQUIRE(graph.has_edge(NodeSide(g_node->id(), true), NodeSide(aa_node->id(), false)));
                }
                
                SECTION("the right side of the AA node should be connected to the left side of the T node") {
                    REQUIRE(graph.has_edge(NodeSide(aa_node->id(), true), NodeSide(t_node->id(), false)));
                }
            }
        }
        
        SECTION("the bluntified graph should have 2 edges") {
            REQUIRE(graph.edge_count() == 2);
            
            SECTION("no edge should have overlap") {
                graph.for_each_edge([&](Edge* e) {
                    REQUIRE(e->overlap() == 0);
                });
            }
        }
        
    }
    
    SECTION("overlaps should be able to overlap in the middle across reversing edges") {
        const string graph_json = R"(
        
        {
            "node": [
                {"id": 1, "sequence": "TTC"},
                {"id": 2, "sequence": "AA"},
                {"id": 3, "sequence": "AAT"}
            ],
            "edge": [
                {"from": 1, "to": 2, "overlap": 2, "from_start": true},
                {"from": 2, "to": 3, "overlap": 2}
            ]
        }
    
        )";
        
        VG graph = string_to_graph(graph_json);
        graph.bluntify();
        
        SECTION("the bluntified graph should have 3 nodes") {
            REQUIRE(graph.node_count() == 3);
            
            SECTION("their sequences should be C, TT or AA, and T") {
                Node* c_node = nullptr;
                Node* middle_node = nullptr;
                bool is_tt;
                Node* t_node = nullptr;
                
                graph.for_each_node([&](Node* n) {
                    if (n->sequence() == "C") {
                        c_node = n;
                    } else if (n->sequence() == "TT") {
                        middle_node = n;
                        is_tt = true;
                    } else if (n->sequence() == "AA") {
                        middle_node = n;
                        is_tt = false;
                    } else if (n->sequence() == "T") {
                        t_node = n;
                    }
                });
                
                REQUIRE(c_node != nullptr);
                REQUIRE(middle_node != nullptr);
                REQUIRE(t_node != nullptr);
                
                SECTION("the left side of the C node should be connected to the right/left side of the TT/AA node") {
                    REQUIRE(graph.has_edge(NodeSide(c_node->id(), false), NodeSide(middle_node->id(), is_tt)));
                }
                
                SECTION("the left/right side of the TT/AA node should be connected to the left side of the T node") {
                    REQUIRE(graph.has_edge(NodeSide(middle_node->id(), !is_tt), NodeSide(t_node->id(), false)));
                }
            }
        }
        
        SECTION("the bluntified graph should have 2 edges") {
            REQUIRE(graph.edge_count() == 2);
            
            SECTION("no edge should have overlap") {
                graph.for_each_edge([&](Edge* e) {
                    REQUIRE(e->overlap() == 0);
                });
            }
        }
        
    }
    
    SECTION("non-overlapping edges should be preserved") {
        const string graph_json = R"(
        
        {
            "node": [
                {"id": 1, "sequence": "CAAAA"},
                {"id": 2, "sequence": "AAAT"},
                {"id": 3, "sequence": "GGG"},
                {"id": 4, "sequence": "CC"}
            ],
            "edge": [
                {"from": 1, "to": 2, "overlap": 3},
                {"from": 3, "to": 1},
                {"from": 2, "to": 4} 
            ]
        }
    
        )";
        
        VG graph = string_to_graph(graph_json);
        graph.bluntify();
        graph.unchop();
        
        SECTION("the unchopped bluntified graph should have one node") {
            REQUIRE(graph.node_count() == 1);
            
            Node* the_node = nullptr;
            graph.for_each_node([&](Node* n) {
                the_node = n;
            });
            
            REQUIRE(the_node != nullptr);
            
            SECTION("that node should be GGGCAAAATCC") {
                REQUIRE(the_node->sequence() == "GGGCAAAATCC");
            }
        }
        
        SECTION("the unchopped bluntified graph has no edges") {
            REQUIRE(graph.edge_count() == 0);
        }
        
    }
}

}
}
