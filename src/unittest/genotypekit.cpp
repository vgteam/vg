/**
 * unittest/genotypekit.cpp: test cases for genotypekit modular genotyper pieces and traversal finders
 */

#include "catch.hpp"
#include "../genotypekit.hpp"
#include "../snarls.hpp"
#include "../cactus_snarl_finder.hpp"
#include "../integrated_snarl_finder.hpp"
#include "../traversal_finder.hpp"
#include "xg.hpp"
#include "../haplotype_extracter.hpp"

namespace Catch {

// Make it so we can see node pointers
std::string toString(vg::Node* const& value) {
  if(value != nullptr) {
    return "&<Node " + std::to_string(value->id()) + ">";
  } else {
    return "&<null node>";
  }
}
    
// Make it so we can see edge pointers
std::string toString(vg::Edge* const& value) {
  if(value != nullptr) {
    auto sides = vg::NodeSide::pair_from_edge(value);
    std::stringstream stream;
    stream << "&<Edge " << sides.first << " -> " << sides.second << ">";
    return stream.str();
  } else {
    return "&<null edge>";
  }
}

    
// And so that we can see sets of things
template<typename Item> struct StringMaker<std::set<Item>> {
  static std::string convert(std::set<Item> const& value) {
    std::stringstream stream;
    stream << "{";
            
            
    for(auto i = value.begin(); i != value.end(); ++i) {
      if(i != value.begin()) {
        // We need the separator
        stream << ", ";
      }
                
      stream << toString(*i);
    }
            
    stream << "}";
    return stream.str();
  }
};
    
}

namespace vg {
namespace unittest {

TEST_CASE("sites can be found with Cactus", "[genotype]") {
    
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
    
  // Make a CactusSnarlFinder
  unique_ptr<SnarlFinder> finder(new CactusSnarlFinder(graph));
    
  SECTION("CactusSnarlFinder should find two top-level sites") {
        
    SnarlManager manager = finder->find_snarls();
        
    auto sites = manager.top_level_snarls();
        
    REQUIRE(sites.size() == 2);
        
    // Order them
    const Snarl* site_1 = sites[0]->start().node_id() > sites[1]->start().node_id() ? sites[1] : sites[0];
    const Snarl* site_2 = sites[0]->start().node_id() > sites[1]->start().node_id() ? sites[0] : sites[1];
        
    SECTION("the first site should be 1 fwd to 6 fwd") {
      REQUIRE(site_1->start().node_id() == 1);
      REQUIRE(site_1->start().backward() == false);
      REQUIRE(site_1->end().node_id() == 6);
      REQUIRE(site_1->end().backward() == false);
            
      SECTION("and should contain exactly nodes 1 through 6") {
        auto nodes = manager.deep_contents(site_1, graph, true).first;
        set<Node*> correct{graph.get_node(1), graph.get_node(2),
             graph.get_node(3), graph.get_node(4),
             graph.get_node(5), graph.get_node(6)};
                
        REQUIRE(nodes.size() == correct.size());
        for (id_t node_id : nodes) {
            REQUIRE(correct.count(graph.get_node(node_id)));
        }
      }
            
      SECTION("and should contain exactly edges 1->6, 1->2, 2->3, 2->4, 3->5, 4->5, 5->6") {
          auto edges = manager.deep_contents(site_1, graph, true).second;
          set<Edge*> correct{
              graph.get_edge(NodeSide(1, true), NodeSide(6)),
                  graph.get_edge(NodeSide(1, true), NodeSide(2)),
                  graph.get_edge(NodeSide(2, true), NodeSide(3)),
                  graph.get_edge(NodeSide(2, true), NodeSide(4)),
                  graph.get_edge(NodeSide(3, true), NodeSide(5)),
                  graph.get_edge(NodeSide(4, true), NodeSide(5)),
                  graph.get_edge(NodeSide(5, true), NodeSide(6))
                  };
                
          REQUIRE(edges.size() == correct.size());
          for (const edge_t& edge_handle : edges) {
              REQUIRE(correct.count(graph.get_edge(NodeTraversal(graph.get_node(graph.get_id(edge_handle.first)),
                                                                 graph.get_is_reverse(edge_handle.first)),
                                                   NodeTraversal(graph.get_node(graph.get_id(edge_handle.second)),
                                                                 graph.get_is_reverse(edge_handle.second)))));
          }
      }
            
      SECTION("and should contain one child") {
        REQUIRE(manager.children_of(site_1).size() == 1);
                
        auto& child = manager.children_of(site_1)[0];
                
        SECTION("that child should be 2 fwd to 5 fwd") {
          REQUIRE(child->start().node_id() == 2);
          REQUIRE(child->start().backward() == false);
          REQUIRE(child->end().node_id() == 5);
          REQUIRE(child->end().backward() == false);
        }
      }
            
    }
        
    SECTION("the second site should be 6 fwd to 9 fwd") {
      REQUIRE(site_2->start().node_id() == 6);
      REQUIRE(site_2->start().backward() == false);
      REQUIRE(site_2->end().node_id() == 9);
      REQUIRE(site_2->end().backward() == false);
            
      SECTION("and should contain no children") {
        REQUIRE(manager.children_of(site_2).size() == 0);
      }
    }
        
        
  }
    
}

TEST_CASE("sites can be found with the IntegratedSnarlFinder", "[genotype][integrated-snarl-finder]") {
    
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
    
  // Make an IntegratedSnarlFinder
  unique_ptr<SnarlFinder> finder(new IntegratedSnarlFinder(graph));
    
  SECTION("IntegratedSnarlFinder should find two top-level sites") {
        
    SnarlManager manager = finder->find_snarls();
        
    auto sites = manager.top_level_snarls();
        
    REQUIRE(sites.size() == 2);
        
    // Order them
    const Snarl* site_1 = sites[0]->start().node_id() > sites[1]->start().node_id() ? sites[1] : sites[0];
    const Snarl* site_2 = sites[0]->start().node_id() > sites[1]->start().node_id() ? sites[0] : sites[1];
        
    SECTION("the first site should be 1 fwd to 6 fwd") {
      REQUIRE(site_1->start().node_id() == 1);
      REQUIRE(site_1->start().backward() == false);
      REQUIRE(site_1->end().node_id() == 6);
      REQUIRE(site_1->end().backward() == false);
            
      SECTION("and should contain exactly nodes 1 through 6") {
        auto nodes = manager.deep_contents(site_1, graph, true).first;
        set<Node*> correct{graph.get_node(1), graph.get_node(2),
             graph.get_node(3), graph.get_node(4),
             graph.get_node(5), graph.get_node(6)};
                
        REQUIRE(nodes.size() == correct.size());
        for (id_t node_id : nodes) {
            REQUIRE(correct.count(graph.get_node(node_id)));
        }
      }
            
      SECTION("and should contain exactly edges 1->6, 1->2, 2->3, 2->4, 3->5, 4->5, 5->6") {
          auto edges = manager.deep_contents(site_1, graph, true).second;
          set<Edge*> correct{
              graph.get_edge(NodeSide(1, true), NodeSide(6)),
                  graph.get_edge(NodeSide(1, true), NodeSide(2)),
                  graph.get_edge(NodeSide(2, true), NodeSide(3)),
                  graph.get_edge(NodeSide(2, true), NodeSide(4)),
                  graph.get_edge(NodeSide(3, true), NodeSide(5)),
                  graph.get_edge(NodeSide(4, true), NodeSide(5)),
                  graph.get_edge(NodeSide(5, true), NodeSide(6))
                  };
                
          REQUIRE(edges.size() == correct.size());
          for (const edge_t& edge_handle : edges) {
              REQUIRE(correct.count(graph.get_edge(NodeTraversal(graph.get_node(graph.get_id(edge_handle.first)),
                                                                 graph.get_is_reverse(edge_handle.first)),
                                                   NodeTraversal(graph.get_node(graph.get_id(edge_handle.second)),
                                                                 graph.get_is_reverse(edge_handle.second)))));
          }
      }
            
      SECTION("and should contain one child") {
        REQUIRE(manager.children_of(site_1).size() == 1);
                
        auto& child = manager.children_of(site_1)[0];
                
        SECTION("that child should be 2 fwd to 5 fwd") {
          REQUIRE(child->start().node_id() == 2);
          REQUIRE(child->start().backward() == false);
          REQUIRE(child->end().node_id() == 5);
          REQUIRE(child->end().backward() == false);
        }
      }
            
    }
        
    SECTION("the second site should be 6 fwd to 9 fwd") {
      REQUIRE(site_2->start().node_id() == 6);
      REQUIRE(site_2->start().backward() == false);
      REQUIRE(site_2->end().node_id() == 9);
      REQUIRE(site_2->end().backward() == false);
            
      SECTION("and should contain no children") {
        REQUIRE(manager.children_of(site_2).size() == 0);
      }
    }
        
        
  }
}

TEST_CASE("IntegratedSnarlFinder works when cactus graph contains back-to-back cycles along root path", "[genotype][integrated-snarl-finder]") {
    
    // Build a toy graph
    const string graph_json = R"(

    {
        "node": [
            {"id": 1, "sequence": "GGGGG"},
            {"id": 2, "sequence": "A"},
            {"id": 3, "sequence": "C"},
            {"id": 4, "sequence": "G"},
            {"id": 5, "sequence": "T"},
            {"id": 6, "sequence": "AAAAA"}
        ],
        "edge": [
            {"from": 1, "to": 2},
            {"from": 1, "to": 3},
            {"from": 2, "to": 4},
            {"from": 2, "to": 5},
            {"from": 3, "to": 4},
            {"from": 3, "to": 5},
            {"from": 4, "to": 6},
            {"from": 5, "to": 6}
            
        ]
    }

    )";

    // Make an actual graph
    VG graph;
    Graph chunk;
    json2pb(chunk, graph_json.c_str(), graph_json.size());
    graph.merge(chunk);

    // Make an IntegratedSnarlFinder
    unique_ptr<SnarlFinder> finder(new IntegratedSnarlFinder(graph));

    SnarlManager manager = finder->find_snarls();
        
    auto sites = manager.top_level_snarls();
        
    // There should just be 1 top snarl, with an internal adjacency component
    REQUIRE(sites.size() == 1);
}

TEST_CASE("IntegratedSnarlFinder works on an all bridge edge Y graph with specific numbering", "[genotype][integrated-snarl-finder]") {
    
    // Build a toy graph
    const string graph_json = R"(
    {"node":[{"id":"2","sequence":"G"},{"id":"3","sequence":"G"},{"id":"4","sequence":"G"},{"id":"5","sequence":"G"},{"id":"6","sequence":"G"},{"id":"11","sequence":"G"}],
    "edge":[{"from":"2","to":"3"},{"from":"3","to":"6"},{"from":"4","to":"5"},{"from":"5","to":"6"},{"from":"6","to":"11"}]}    
    )";

    // Make an actual graph
    VG graph;
    Graph chunk;
    json2pb(chunk, graph_json.c_str(), graph_json.size());
    graph.merge(chunk);

    // Make an IntegratedSnarlFinder
    unique_ptr<SnarlFinder> finder(new IntegratedSnarlFinder(graph));

    SnarlManager manager = finder->find_snarls();
        
    auto sites = manager.top_level_snarls();
        
    // There should be 3 snarls in a chain, and 1 nested snarl
    REQUIRE(sites.size() == 3);
    REQUIRE(manager.num_snarls() == 4);
    
    // We don't care which pair of Y branches it roots along.
}

TEST_CASE("IntegratedSnarlFinder roots correctly an all bridge edge Y graph with winning longest path", "[genotype][integrated-snarl-finder]") {
    
    // Build a toy graph
    const string graph_json = R"(
    {"node":[{"id":"2","sequence":"G"},{"id":"3","sequence":"G"},{"id":"4","sequence":"GG"},{"id":"5","sequence":"G"},{"id":"6","sequence":"G"},{"id":"11","sequence":"GG"}],
    "edge":[{"from":"2","to":"3"},{"from":"3","to":"6"},{"from":"4","to":"5"},{"from":"5","to":"6"},{"from":"6","to":"11"}]}    
    )";

    // Make an actual graph
    VG graph;
    Graph chunk;
    json2pb(chunk, graph_json.c_str(), graph_json.size());
    graph.merge(chunk);

    // Make an IntegratedSnarlFinder
    unique_ptr<SnarlFinder> finder(new IntegratedSnarlFinder(graph));

    SnarlManager manager = finder->find_snarls();
        
    auto sites = manager.top_level_snarls();
        
    // There should be 3 snarls in a chain, and 1 nested snarl
    REQUIRE(sites.size() == 3);
    REQUIRE(manager.num_snarls() == 4);
    
    // We don't care which pair of Y branches it roots along.
    
    // Top snarls should have the right bounds
    REQUIRE(sites[0]->start().node_id() == 6);
    REQUIRE(sites[0]->start().backward() == false);
    REQUIRE(sites[0]->end().node_id() == 11);
    REQUIRE(sites[0]->end().backward() == false);
    // These are coming out backward vs chain order
    REQUIRE(manager.chain_rank_of(sites[0]) == 2);
    
    REQUIRE(sites[1]->start().node_id() == 5);
    REQUIRE(sites[1]->start().backward() == false);
    REQUIRE(sites[1]->end().node_id() == 6);
    REQUIRE(sites[1]->end().backward() == false);
    REQUIRE(manager.chain_rank_of(sites[1]) == 1);
    
    REQUIRE(sites[2]->start().node_id() == 4);
    REQUIRE(sites[2]->start().backward() == false);
    REQUIRE(sites[2]->end().node_id() == 5);
    REQUIRE(sites[2]->end().backward() == false);
    REQUIRE(manager.chain_rank_of(sites[2]) == 0);
    
}

TEST_CASE("IntegratedSnarlFinder works when cactus graph contains longer back-to-back cycles along root path", "[genotype][integrated-snarl-finder]") {
    
    // Build a toy graph
    const string graph_json = R"(

    {
        "node": [
            {"id": 1, "sequence": "GGGGG"},
            {"id": 2, "sequence": "A"},
            {"id": 21, "sequence": "A"},
            {"id": 22, "sequence": "A"},
            {"id": 3, "sequence": "C"},
            {"id": 31, "sequence": "C"},
            {"id": 32, "sequence": "C"},
            {"id": 4, "sequence": "G"},
            {"id": 5, "sequence": "T"},
            {"id": 6, "sequence": "AAAAA"}
        ],
        "edge": [
            {"from": 1, "to": 2},
            {"from": 1, "to": 3},
            {"from": 2, "to": 21},
            {"from": 21, "to": 22},
            {"from": 22, "to": 4},
            {"from": 22, "to": 5},
            {"from": 3, "to": 31},
            {"from": 31, "to": 32},
            {"from": 32, "to": 4},
            {"from": 32, "to": 5},
            {"from": 4, "to": 6},
            {"from": 5, "to": 6}
            
        ]
    }

    )";

    // Make an actual graph
    VG graph;
    Graph chunk;
    json2pb(chunk, graph_json.c_str(), graph_json.size());
    graph.merge(chunk);

    // Make an IntegratedSnarlFinder
    unique_ptr<SnarlFinder> finder(new IntegratedSnarlFinder(graph));

    SnarlManager manager = finder->find_snarls();
        
    auto sites = manager.top_level_snarls();
        
    // There should just be 1 top snarl, with an internal adjacency component
    REQUIRE(sites.size() == 1);
}

TEST_CASE("IntegratedSnarlFinder works on a complex bundle-y region with a nested snarl", "[genotype][integrated-snarl-finder]") {
    
    // Build a toy graph
    const string graph_json = R"(
        {"edge": [{"from": "129672", "to": "129673"},
                  {"from": "129662", "to": "129663"}, 
                  {"from": "129662", "to": "129664"}, 
                  {"from": "129664", "to": "129665"}, 
                  {"from": "129664", "to": "129666"}, 
                  {"from": "129666", "to": "129668"}, 
                  {"from": "129666", "to": "129669"}, 
                  {"from": "129666", "to": "129667"}, 
                  {"from": "129667", "to": "129668"}, 
                  {"from": "129667", "to": "129669"}, 
                  {"from": "129669", "to": "129670"}, 
                  {"from": "129669", "to": "129673"}, 
                  {"from": "129671", "to": "129672"}, 
                  {"from": "129668", "to": "129670"}, 
                  {"from": "129668", "to": "129673"}, 
                  {"from": "129665", "to": "129668"}, 
                  {"from": "129665", "to": "129669"}, 
                  {"from": "129665", "to": "129667"}, 
                  {"from": "129670", "to": "129671"}, 
                  {"from": "129670", "to": "129672"}, 
                  {"from": "129663", "to": "129665"}, 
                  {"from": "129663", "to": "129666"}], 
        "node": [{"id": "129672", "sequence": "AT"}, 
                 {"id": "129662", "sequence": "CAGGTCAAACTGTGAT"}, 
                 {"id": "129664", "sequence": "T"}, 
                 {"id": "129666", "sequence": "T"}, 
                 {"id": "129667", "sequence": "G"}, 
                 {"id": "129669", "sequence": "G"}, 
                 {"id": "129671", "sequence": "T"}, 
                 {"id": "129668", "sequence": "A"}, 
                 {"id": "129665", "sequence": "A"}, 
                 {"id": "129670", "sequence": "A"}, 
                 {"id": "129673", "sequence": "ATATATATATACTTATTGTAAAAATCTTTAGA"}, 
                 {"id": "129663", "sequence": "G"}]}
    )";

    // Make an actual graph
    VG graph;
    Graph chunk;
    json2pb(chunk, graph_json.c_str(), graph_json.size());
    graph.merge(chunk);

    // Make an IntegratedSnarlFinder
    unique_ptr<SnarlFinder> finder(new IntegratedSnarlFinder(graph));

    SnarlManager manager = finder->find_snarls();
        
    auto sites = manager.top_level_snarls();
        
    // There should just be 1 top snarl.
    REQUIRE(sites.size() == 1);
    
    // It should have the right bounds
    REQUIRE(sites[0]->start().node_id() == 129662);
    REQUIRE(sites[0]->start().backward() == false);
    REQUIRE(sites[0]->end().node_id() == 129673);
    REQUIRE(sites[0]->end().backward() == false);
           
    // It should have one child snarl
    REQUIRE(manager.children_of(sites[0]).size() == 1);
    auto& child = manager.children_of(sites[0])[0];
    
    // And the child should have the right bounds
    REQUIRE(child->start().node_id() == 129670);
    REQUIRE(child->start().backward() == false);
    REQUIRE(child->end().node_id() == 129672);
    REQUIRE(child->end().backward() == false);
}

TEST_CASE("CactusSnarlFinder safely handles a single node graph", "[genotype][cactus-snarl-finder]") {
    
  // Build a toy graph
  const string graph_json = R"(
    
    {
        "node": [
            {"id": 1, "sequence": "GATTACA"}
        ]
    }
    
    )";
    
  // Make an actual graph
  VG graph;
  Graph chunk;
  json2pb(chunk, graph_json.c_str(), graph_json.size());
  graph.merge(chunk);
    
  // Make a CactusSnarlFinder
  unique_ptr<SnarlFinder> finder(new CactusSnarlFinder(graph));
    
  SECTION("CactusSnarlFinder returns empty snarl manager instead of throwing or crashing") {
    REQUIRE(finder->find_snarls().num_snarls() == 0);
  }
  
}

TEST_CASE("IntegratedSnarlFinder safely handles a completely empty graph", "[genotype][integrated-snarl-finder]") {
    
    // Build a toy graph
    const string graph_json = "{}";

    // Make an actual graph
    VG graph;
    Graph chunk;
    json2pb(chunk, graph_json.c_str(), graph_json.size());
    graph.merge(chunk);

    // Make a IntegratedSnarlFinder
    unique_ptr<SnarlFinder> finder(new IntegratedSnarlFinder(graph));

    // There should be no snarls but the root implicit snarl.
    REQUIRE(finder->find_snarls().num_snarls() == 0);
}

TEST_CASE("IntegratedSnarlFinder safely handles a single node graph", "[genotype][integrated-snarl-finder]") {
    
    // Build a toy graph
    const string graph_json = R"(

    {
        "node": [
            {"id": 1, "sequence": "GATTACA"}
        ]
    }

    )";

    // Make an actual graph
    VG graph;
    Graph chunk;
    json2pb(chunk, graph_json.c_str(), graph_json.size());
    graph.merge(chunk);

    // Make an IntegratedSnarlFinder
    unique_ptr<SnarlFinder> finder(new IntegratedSnarlFinder(graph));

    // There should be no snarls but the root implicit snarl.
    REQUIRE(finder->find_snarls().num_snarls() == 0);
}

TEST_CASE("IntegratedSnarlFinder produces all the correct types of single-node chains", "[genotype][integrated-snarl-finder]") {
    
    // Build a toy graph
    const string graph_json = R"(

    {
        "node": [
            {"id": 1, "sequence": "GATTACA"},
            {"id": 2, "sequence": "GATT"},
            {"id": 3, "sequence": "ACA"},
            {"id": 4, "sequence": "T"},
            {"id": 5, "sequence": "A"}
        ], "edge": [
            {"from": 2, "to": 3},
            {"from": 2, "to": 4},
            {"from": 4, "to": 3},
            {"from": 2, "to": 5}
        ]
    }

    )";

    // Make an actual graph
    VG graph;
    Graph chunk;
    json2pb(chunk, graph_json.c_str(), graph_json.size());
    graph.merge(chunk);

    // Make an IntegratedSnarlFinder
    IntegratedSnarlFinder finder(graph);
    
    // We will log evberything the snarl finder shows us.
    unordered_set<handle_t> chain_starts;
    unordered_set<handle_t> chain_ends;
    unordered_set<handle_t> snarl_starts;
    unordered_set<handle_t> snarl_ends;
    
    // Chain visits are specified to happen for single-node chains
    finder.traverse_decomposition([&](const handle_t& chain_start) {
        // Deal with the start of a chain
#ifdef debug
        std::cerr << "Chain start: " << graph.get_id(chain_start) << (graph.get_is_reverse(chain_start) ? "R" : "F") << std::endl;
#endif
        chain_starts.insert(chain_start);
    }, [&](const handle_t& chain_end){
        // Deal with the end of a chain
#ifdef debug
        std::cerr << "Chain end: " << graph.get_id(chain_end) << (graph.get_is_reverse(chain_end) ? "R" : "F") << std::endl;
#endif
        chain_ends.insert(chain_end);
    }, [&](const handle_t& snarl_start){
        // Deal with the start of a snarl
#ifdef debug
        std::cerr << "Snarl start: " << graph.get_id(snarl_start) << (graph.get_is_reverse(snarl_start) ? "R" : "F") << std::endl;
#endif
        snarl_starts.insert(snarl_start);
    }, [&](const handle_t& snarl_end) {
        // Deal with the end of a snarl
#ifdef debug
        std::cerr << "Snarl end: " << graph.get_id(snarl_end) << (graph.get_is_reverse(snarl_end) ? "R" : "F") << std::endl;
#endif
        snarl_ends.insert(snarl_end);
    });
    
    // The trivial snarl should be found
    REQUIRE((snarl_starts.count(graph.get_handle(2)) || snarl_starts.count(graph.get_handle(3, true))));
    REQUIRE((snarl_ends.count(graph.get_handle(3)) || snarl_ends.count(graph.get_handle(2, true))));
    
    // It should also be found as contained in a chain
    REQUIRE((chain_starts.count(graph.get_handle(2)) || chain_starts.count(graph.get_handle(3, true))));
    REQUIRE((chain_ends.count(graph.get_handle(3)) || chain_ends.count(graph.get_handle(2, true))));
    
    // The single-node chain should be found for the single-node component
    REQUIRE((chain_starts.count(graph.get_handle(1)) || chain_starts.count(graph.get_handle(1, true))));
    REQUIRE((chain_ends.count(graph.get_handle(1)) || chain_ends.count(graph.get_handle(1, true))));
    
    // The single-node chain should be found for the contained node in the snarl
    REQUIRE((chain_starts.count(graph.get_handle(4)) || chain_starts.count(graph.get_handle(4, true))));
    REQUIRE((chain_ends.count(graph.get_handle(4)) || chain_ends.count(graph.get_handle(4, true))));
    
    // The single-node chain should be found for the dangling bridge in the snarl
    REQUIRE((chain_starts.count(graph.get_handle(5)) || chain_starts.count(graph.get_handle(5, true))));
    REQUIRE((chain_ends.count(graph.get_handle(5)) || chain_ends.count(graph.get_handle(5, true))));
}

TEST_CASE("IntegratedSnarlFinder safely handles a path when forced to root at one end", "[genotype][integrated-snarl-finder]") {
    
    // Build a toy graph
    const string graph_json = R"(

    {
        "node": [
            {"id": 1, "sequence": "GATTACA"},
            {"id": 2, "sequence": "CATTAG"},
            {"id": 3, "sequence": "A"},
            {"id": 4, "sequence": "A"}
        ], "edge": [
            {"from": 1, "to": 2},
            {"from": 2, "to": 3},
            {"from": 3, "to": 4},
            {"from": 4, "to": 3}
        ]
    }

    )";

    // Make an actual graph
    VG graph;
    Graph chunk;
    json2pb(chunk, graph_json.c_str(), graph_json.size());
    graph.merge(chunk);

    // Make a CactusSnarlFinder
    unique_ptr<SnarlFinder> finder(new IntegratedSnarlFinder(graph));

    // There should be a snarl along the body and a snarl for the cycle at the tip.
    REQUIRE(finder->find_snarls().num_snarls() == 2);
}

TEST_CASE("IntegratedSnarlFinder safely handles a single node connected component in a larger graph", "[genotype][integrated-snarl-finder]") {
    
    // Build a toy graph
    const string graph_json = R"(

    {
        "node": [
            {"id": 1, "sequence": "GATTACA"},
            {"id": 2, "sequence": "GAT"},
            {"id": 3, "sequence": "TACA"}
        ], "edge": [
            {"from": 2, "to": 3}
        ]
    }

    )";

    // Make an actual graph
    VG graph;
    Graph chunk;
    json2pb(chunk, graph_json.c_str(), graph_json.size());
    graph.merge(chunk);

    // Make an IntegratedSnarlFinder
    unique_ptr<SnarlFinder> finder(new IntegratedSnarlFinder(graph));
    
    auto manager = finder->find_snarls();
    
    // We expect the one trivial snarl only.
    REQUIRE(manager.num_snarls() == 1);
    
    auto sites = manager.top_level_snarls();
        
    // There should just be 1 top snarl.
    REQUIRE(sites.size() == 1);
    
    // It should have the right bounds
    REQUIRE(sites[0]->start().node_id() == 2);
    REQUIRE(sites[0]->start().backward() == false);
    REQUIRE(sites[0]->end().node_id() == 3);
    REQUIRE(sites[0]->end().backward() == false);
}

TEST_CASE("IntegratedSnarlFinder safely handles a single node cycle", "[genotype][integrated-snarl-finder]") {
    
    // Build a toy graph
    const string graph_json = R"(

    {
        "node": [
            {"id": 1, "sequence": "GATTACA"}
        ], "edge": [
            {"from": 1, "to": 1}
        ]
    }

    )";

    // Make an actual graph
    VG graph;
    Graph chunk;
    json2pb(chunk, graph_json.c_str(), graph_json.size());
    graph.merge(chunk);

    // Make an IntegratedSnarlFinder
    unique_ptr<SnarlFinder> finder(new IntegratedSnarlFinder(graph));
    
    auto manager = finder->find_snarls();
    
    // Like other fully-connected cases we regard this as contents of the root
    // snarl. So no snarls come out.
    REQUIRE(manager.num_snarls() == 0);
}

TEST_CASE("IntegratedSnarlFinder safely handles a totally connected graph", "[genotype][integrated-snarl-finder]") {
    
    // Build a toy graph
    const string graph_json = R"(

    {
        "node": [
            {"id": 1, "sequence": "GAT"},
            {"id": 2, "sequence": "TA"},
            {"id": 3, "sequence": "CA"}
        ], "edge": [
            {"from": 1, "to": 2},
            {"from": 1, "to": 3},
            {"from": 2, "to": 3},
            {"from": 2, "to": 1},
            {"from": 3, "to": 1},
            {"from": 3, "to": 2}
        ]
    }

    )";

    // Make an actual graph
    VG graph;
    Graph chunk;
    json2pb(chunk, graph_json.c_str(), graph_json.size());
    graph.merge(chunk);

    // Make an IntegratedSnarlFinder
    unique_ptr<SnarlFinder> finder(new IntegratedSnarlFinder(graph));
    
    auto manager = finder->find_snarls();
    
    // There should be no snarls but the root snarl where everything is just contents.
    REQUIRE(manager.num_snarls() == 0);
    
}

TEST_CASE("IntegratedSnarlFinder prefers to root at a bridge edge path in a tie", "[genotype][integrated-snarl-finder]") {
    
    // Build a toy graph
    const string graph_json = R"(

    {
        "node": [
            {"id": 1, "sequence": "GATT"},
            {"id": 2, "sequence": "GAT"},
            {"id": 3, "sequence": "TACA"},
            {"id": 4, "sequence": "ACA"}
        ], "edge": [
            {"from": 2, "to": 3},
            {"from": 3, "to": 2},
            {"from": 2, "to": 1},
            {"from": 1, "to": 4}
        ]
    }

    )";

    // Make an actual graph
    VG graph;
    Graph chunk;
    json2pb(chunk, graph_json.c_str(), graph_json.size());
    graph.merge(chunk);

    // Make an IntegratedSnarlFinder
    unique_ptr<SnarlFinder> finder(new IntegratedSnarlFinder(graph));
    
    auto manager = finder->find_snarls();
    
    // There should be 2 total snarls.
    REQUIRE(manager.num_snarls() == 2);
    
    auto sites = manager.top_level_snarls();
    
    // Both should be top-level, with connectivity in the root snarl
    REQUIRE(sites.size() == 2);
    
    // The cycle snarl should be set up so the bridge path snarl is outside it.
    REQUIRE(sites[0]->start().node_id() == 2);
    REQUIRE(sites[0]->start().backward() == true);
    REQUIRE(sites[0]->end().node_id() == 3);
    REQUIRE(sites[0]->end().backward() == true);
    
    // The bridge chain snarl should have the right bounds
    REQUIRE(sites[1]->start().node_id() == 1);
    REQUIRE(sites[1]->start().backward() == false);
    REQUIRE(sites[1]->end().node_id() == 4);
    REQUIRE(sites[1]->end().backward() == false);
}

TEST_CASE("IntegratedSnarlFinder prefers to root at a cycle that is 1 bp longer", "[genotype][integrated-snarl-finder]") {
    
    // Build a toy graph
    const string graph_json = R"(

    {
        "node": [
            {"id": 1, "sequence": "GATT"},
            {"id": 2, "sequence": "GAT"},
            {"id": 3, "sequence": "TAACA"},
            {"id": 4, "sequence": "ACA"}
        ], "edge": [
            {"from": 2, "to": 3},
            {"from": 3, "to": 2},
            {"from": 2, "to": 1},
            {"from": 1, "to": 4}
        ]
    }

    )";

    // Make an actual graph
    VG graph;
    Graph chunk;
    json2pb(chunk, graph_json.c_str(), graph_json.size());
    graph.merge(chunk);

    // Make an IntegratedSnarlFinder
    unique_ptr<SnarlFinder> finder(new IntegratedSnarlFinder(graph));
    
    auto manager = finder->find_snarls();
    
    // There should be 3 total snarls
    REQUIRE(manager.num_snarls() == 3);
    
    auto sites = manager.top_level_snarls();
    
    // 2 of them should be top-level
    REQUIRE(sites.size() == 2);
    
    // They should both be around the cycle
    
    REQUIRE(sites[0]->start().node_id() == 2);
    REQUIRE(sites[0]->start().backward() == false);
    REQUIRE(sites[0]->end().node_id() == 3);
    REQUIRE(sites[0]->end().backward() == false);
    
    REQUIRE(sites[1]->start().node_id() == 3);
    REQUIRE(sites[1]->start().backward() == false);
    REQUIRE(sites[1]->end().node_id() == 2);
    REQUIRE(sites[1]->end().backward() == false);
}

TEST_CASE("IntegratedSnarlFinder sees tips as disqualifying ultrabubbles", "[genotype][integrated-snarl-finder]") {
    
    // Build a toy graph
    const string graph_json = R"(

    {
        "node": [
            {"id": 1, "sequence": "GGGGG"},
            {"id": 2, "sequence": "A"},
            {"id": 3, "sequence": "C"},
            {"id": 4, "sequence": "G"},
            {"id": 5, "sequence": "T"},
            {"id": 6, "sequence": "AAAAA"}
        ],
        "edge": [
            {"from": 1, "to": 2},
            {"from": 2, "to": 5},
            {"from": 5, "to": 6},
            {"from": 3, "to": 4},
            {"from": 4, "to": 5}
        ]
    }

    )";

    // Make an actual graph
    VG graph;
    Graph chunk;
    json2pb(chunk, graph_json.c_str(), graph_json.size());
    graph.merge(chunk);

    // Make an IntegratedSnarlFinder
    unique_ptr<SnarlFinder> finder(new IntegratedSnarlFinder(graph));

    SnarlManager manager = finder->find_snarls();
        
    auto sites = manager.top_level_snarls();
        
    // There should just be 3 top level snarls, two trivial and the one we care about.
    REQUIRE(sites.size() == 3);
    
    // Grab the snarl we care about
    const Snarl* snarl = manager.into_which_snarl(2, false);
    
    REQUIRE(snarl != nullptr);
    
    // The snarl shouldn't be an ultrabubble because it has tips.
    REQUIRE(snarl->type() == UNCLASSIFIED);
    
    // Grab some other trivial snarl
    const Snarl* snarl2 = manager.into_which_snarl(1, false);
    
    // Make sure it is still an ultrabubble
    REQUIRE(snarl2 != nullptr);
    REQUIRE(snarl2->type() == ULTRABUBBLE);
}

TEST_CASE("CactusSnarlFinder throws an error instead of crashing when the graph has no edges", "[genotype][cactus-snarl-finder]") {
    
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
        ]
    }
    
    )";
    
  // Make an actual graph
  VG graph;
  Graph chunk;
  json2pb(chunk, graph_json.c_str(), graph_json.size());
  graph.merge(chunk);
    
  // Make a CactusSnarlFinder
  unique_ptr<SnarlFinder> finder(new CactusSnarlFinder(graph));
    
  SECTION("CactusSnarlFinder should fail gracefully") {
    REQUIRE_THROWS(finder->find_snarls());
  }
  
}

TEST_CASE("fixed priors can be assigned to genotypes", "[genotype]") {
    
  GenotypePriorCalculator* calculator = new FixedGenotypePriorCalculator();
    
  Genotype het;
  het.add_allele(0);
  het.add_allele(1);
    
  Genotype hom_alt;
  hom_alt.add_allele(1);
  hom_alt.add_allele(1);
    
  Genotype hom_ref;
  hom_ref.add_allele(0);
  hom_ref.add_allele(0);
    
  SECTION("homozygote priors should be equal") {
    REQUIRE(calculator->calculate_log_prior(hom_alt) == calculator->calculate_log_prior(hom_ref));
  }
    
  SECTION("homozygotes should be more likely than heterozygotes") {
    REQUIRE(calculator->calculate_log_prior(het) < calculator->calculate_log_prior(hom_ref));
    REQUIRE(calculator->calculate_log_prior(het) < calculator->calculate_log_prior(hom_alt));
  }
    
  SECTION("haploid genotypes should have nonzero prior") {
    Genotype haploid;
    haploid.add_allele(5);
    REQUIRE(calculator->calculate_log_prior(haploid) > prob_to_logprob(0));
  }
    
  SECTION("zero-ploid genotypes should have nonzero prior") {
    Genotype empty;
    REQUIRE(calculator->calculate_log_prior(empty) > prob_to_logprob(0));
  }
    
  SECTION("polyploid genotypes should have nonzero prior") {
    Genotype polyploid;
    for(int i = 0; i < 100; i++) {
      polyploid.add_allele(i);
    }
    REQUIRE(calculator->calculate_log_prior(polyploid) > prob_to_logprob(0));
  }
    
  delete calculator;
}

TEST_CASE("TrivialTraversalFinder can find traversals", "[genotype]") {
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
    
  // Make a site
  Snarl site;
  site.mutable_start()->set_node_id(2);
  site.mutable_end()->set_node_id(5);
  site.set_type(ULTRABUBBLE);
  site.set_start_end_reachable(true);
  site.set_directed_acyclic_net_graph(true);
    
  // Make the TraversalFinder
  TraversalFinder* finder = new TrivialTraversalFinder(graph);
    
  SECTION("at least one path can be found") {
    auto site_traversals = finder->find_traversals(site);
        
    REQUIRE(!site_traversals.empty());
        
    SECTION("the path must visit 3 node to span the site") {
      REQUIRE(site_traversals.front().visit_size() == 3);
            
      SECTION("the site must follow one of the two paths in the correct orientation") {
        bool followed_path_1 = site_traversals.front().visit(1).node_id() == 3 &&
           site_traversals.front().visit(1).backward() == false;
                
        bool followed_path_2 = site_traversals.front().visit(1).node_id() == 4 &&
           site_traversals.front().visit(1).backward() == false;
                
        REQUIRE(followed_path_1 != followed_path_2);
      }
    }
  }
    
  delete finder;
    

}
    
TEST_CASE("ExhaustiveTraversalFinder finds all paths on a bubble with an inverted node", "[genotype]") {
  VG graph;
  Node* n0 = graph.create_node("A");
  Node* n1 = graph.create_node("C");
  Node* n2 = graph.create_node("G");
  Node* n3 = graph.create_node("T");
  graph.create_edge(n0, n1, false, true);
  graph.create_edge(n0, n2, false, false);
  graph.create_edge(n3, n1, true, false);
  graph.create_edge(n2, n3, false, false);
    
  Snarl site;
  site.mutable_start()->set_node_id(n0->id());
  site.mutable_end()->set_node_id(n3->id());
  site.set_type(ULTRABUBBLE);
    
  list<Snarl> snarls{site};
    
  SnarlManager manager(snarls.begin(), snarls.end());
    
  ExhaustiveTraversalFinder finder(graph, manager);
    
  bool found_trav_1 = false;
  bool found_trav_2 = false;
    
  for (auto snarl : manager.top_level_snarls()) {
    auto travs = finder.find_traversals(*snarl);
    for (SnarlTraversal trav : travs) {
      if (trav.visit_size() == 3) {
        if (trav.visit(1).node_id() == n1->id() && trav.visit(1).backward()) {
          found_trav_1 = true;
        }
        else if (trav.visit(1).node_id() == n2->id() && !trav.visit(1).backward()) {
          found_trav_2 = true;
        }
      }
    }
  }
    
  REQUIRE(found_trav_1);
  REQUIRE(found_trav_2);
}

TEST_CASE("CactusSnarlFinder can differentiate ultrabubbles from snarls", "[genotype][cactus-snarl-finder]") {

    SECTION("Directed cycle does not count as ultrabubble") {
        // Build a toy graph
        const string graph_json = R"(
        {
        "node": [
            {"id": 1, "sequence": "AGGAGGG"},
            {"id": 2, "sequence": "A"},
            {"id": 3, "sequence": "T"},
            {"id": 4, "sequence": "GGG"},
            {"id": 5, "sequence": "GTTTGTG"}
        ],
        "edge": [
            {"from": 1, "to": 2},
            {"from": 2, "to": 3},
            {"from": 3, "to": 4},
            {"from": 1, "to": 4, "to_end": true},
            {"from": 4, "to": 5}
        ]
        }
        )";
        
        // Make an actual graph
        VG graph;
        Graph chunk;
        json2pb(chunk, graph_json.c_str(), graph_json.size());
        graph.merge(chunk);

        // Find the snarls
        CactusSnarlFinder cubs(graph);
        SnarlManager snarl_manager = cubs.find_snarls();
        const vector<const Snarl*>& snarl_roots = snarl_manager.top_level_snarls();

        // Make sure we have one non-ultrabubble start at 1-5
        REQUIRE(snarl_roots.size() == 1);
        const Snarl* snarl = snarl_roots[0];
        int64_t start = snarl->start().node_id();
        int64_t end = snarl->end().node_id();
        if (start > end) {
          std::swap(start, end);
        }
        REQUIRE((start == 1 && end == 5) == true);
        REQUIRE(snarl->type() != ULTRABUBBLE);
    }

    SECTION("Ultrabubble flagged as ultrabubble") {

        // Build a toy graph
        const string graph_json = R"(
        {
        "node": [
            {"id": 1, "sequence": "G"},
            {"id": 2, "sequence": "A"},
            {"id": 3, "sequence": "T"},
            {"id": 4, "sequence": "GGG"},
            {"id": 5, "sequence": "GT"},
            {"id": 6, "sequence": "GT"}
        ],
        "edge": [
            {"from": 1, "to": 2, "to_end": true},
            {"from": 2, "to": 4, "from_start": true},
            {"from": 2, "to": 5, "from_start": true},
            {"from": 4, "to": 6, "to_end": true},
            {"from": 3, "to": 1, "from_start": true, "to_end": true},
            {"from": 4, "to": 3, "from_start": true, "to_end": true},
            {"from": 5, "to": 3, "from_start": true, "to_end": true},
            {"from": 6, "to": 5, "to_end": true},
            {"from": 1, "to": 6, "from_start": true}
        ]
        }
        )";
    
        // Make an actual graph
        VG graph;
        Graph chunk;
        json2pb(chunk, graph_json.c_str(), graph_json.size());
        graph.merge(chunk);

        // Find the snarls
        CactusSnarlFinder cubs(graph);
        SnarlManager snarl_manager = cubs.find_snarls();
        const vector<const Snarl*>& snarl_roots = snarl_manager.top_level_snarls();

        // Make sure we have one ultrabubble from 1 forward to 6 reverse, and
        // another ultrabubble closing the cycle from 6 reverse to 1 forward.
        REQUIRE(snarl_roots.size() == 2);
        const Snarl* snarl1 = snarl_roots[0];
        const Snarl* snarl2 = snarl_roots[1];
            
        if (snarl1->start().node_id() > snarl1->end().node_id()) {
          // Flip it around so it goes from lower to higher numbers.
          snarl_manager.flip(snarl1);
        }
            
        if (snarl2->start().node_id() > snarl2->end().node_id()) {
          // Flip it around so it goes from lower to higher numbers.
          snarl_manager.flip(snarl2);
        }
            
        if (snarl1->start().node_id() > snarl2->start().node_id() ||
            (snarl1->start().node_id() == snarl2->start().node_id() &&
             snarl1->start().backward() > snarl2->start().backward())) {
          // Make sure to order them deterministically
          std::swap(snarl1, snarl2);
        }
            
        REQUIRE(snarl1->start().node_id() == 1);
        REQUIRE(snarl1->start().backward() == false);
        REQUIRE(snarl1->end().node_id() == 6);
        REQUIRE(snarl1->end().backward() == true);
        REQUIRE(snarl1->type() == ULTRABUBBLE);
            
        REQUIRE(snarl2->start().node_id() == 1);
        REQUIRE(snarl2->start().backward() == true);
        REQUIRE(snarl2->end().node_id() == 6);
        REQUIRE(snarl2->end().backward() == false);
        REQUIRE(snarl2->type() == ULTRABUBBLE);
    }   

}

TEST_CASE("IntegratedSnarlFinder can differentiate ultrabubbles from snarls", "[genotype][integrated-snarl-finder]") {

    SECTION("Directed cycle does not count as ultrabubble") {
        // Build a toy graph
        // Make sure start and end nodes are long enough to suggest the rooting we want.
        const string graph_json = R"(
        {
        "node": [
            {"id": 1, "sequence": "AGGAGGG"},
            {"id": 2, "sequence": "A"},
            {"id": 3, "sequence": "T"},
            {"id": 4, "sequence": "GGG"},
            {"id": 5, "sequence": "GTTTGTG"}
        ],
        "edge": [
            {"from": 1, "to": 2},
            {"from": 2, "to": 3},
            {"from": 3, "to": 4},
            {"from": 1, "to": 4, "to_end": true},
            {"from": 4, "to": 5}
        ]
        }
        )";
        
        // Make an actual graph
        VG graph;
        Graph chunk;
        json2pb(chunk, graph_json.c_str(), graph_json.size());
        graph.merge(chunk);

        // Find the snarls
        IntegratedSnarlFinder cubs(graph);
        SnarlManager snarl_manager = cubs.find_snarls();
        const vector<const Snarl*>& snarl_roots = snarl_manager.top_level_snarls();

        // Make sure we have one non-ultrabubble start at 1-5
        REQUIRE(snarl_roots.size() == 1);
        const Snarl* snarl = snarl_roots[0];
        int64_t start = snarl->start().node_id();
        int64_t end = snarl->end().node_id();
        if (start > end) {
          std::swap(start, end);
        }
        REQUIRE((start == 1 && end == 5) == true);
        REQUIRE(snarl->type() != ULTRABUBBLE);
    }

    SECTION("Ultrabubble flagged as ultrabubble") {

        // Build a toy graph
        const string graph_json = R"(
        {
        "node": [
            {"id": 1, "sequence": "G"},
            {"id": 2, "sequence": "A"},
            {"id": 3, "sequence": "T"},
            {"id": 4, "sequence": "GGG"},
            {"id": 5, "sequence": "GT"},
            {"id": 6, "sequence": "GT"}
        ],
        "edge": [
            {"from": 1, "to": 2, "to_end": true},
            {"from": 2, "to": 4, "from_start": true},
            {"from": 2, "to": 5, "from_start": true},
            {"from": 4, "to": 6, "to_end": true},
            {"from": 3, "to": 1, "from_start": true, "to_end": true},
            {"from": 4, "to": 3, "from_start": true, "to_end": true},
            {"from": 5, "to": 3, "from_start": true, "to_end": true},
            {"from": 6, "to": 5, "to_end": true},
            {"from": 1, "to": 6, "from_start": true}
        ]
        }
        )";
    
        // Make an actual graph
        VG graph;
        Graph chunk;
        json2pb(chunk, graph_json.c_str(), graph_json.size());
        graph.merge(chunk);

        // Find the snarls
        IntegratedSnarlFinder cubs(graph);
        SnarlManager snarl_manager = cubs.find_snarls();
        const vector<const Snarl*>& snarl_roots = snarl_manager.top_level_snarls();

        // Make sure we have one ultrabubble from 1 forward to 6 reverse, and
        // another ultrabubble closing the cycle from 6 reverse to 1 forward.
        REQUIRE(snarl_roots.size() == 2);
        const Snarl* snarl1 = snarl_roots[0];
        const Snarl* snarl2 = snarl_roots[1];
            
        if (snarl1->start().node_id() > snarl1->end().node_id()) {
          // Flip it around so it goes from lower to higher numbers.
          snarl_manager.flip(snarl1);
        }
            
        if (snarl2->start().node_id() > snarl2->end().node_id()) {
          // Flip it around so it goes from lower to higher numbers.
          snarl_manager.flip(snarl2);
        }
            
        if (snarl1->start().node_id() > snarl2->start().node_id() ||
            (snarl1->start().node_id() == snarl2->start().node_id() &&
             snarl1->start().backward() > snarl2->start().backward())) {
          // Make sure to order them deterministically
          std::swap(snarl1, snarl2);
        }
            
        REQUIRE(snarl1->start().node_id() == 1);
        REQUIRE(snarl1->start().backward() == false);
        REQUIRE(snarl1->end().node_id() == 6);
        REQUIRE(snarl1->end().backward() == true);
        REQUIRE(snarl1->type() == ULTRABUBBLE);
            
        REQUIRE(snarl2->start().node_id() == 1);
        REQUIRE(snarl2->start().backward() == true);
        REQUIRE(snarl2->end().node_id() == 6);
        REQUIRE(snarl2->end().backward() == false);
        REQUIRE(snarl2->type() == ULTRABUBBLE);
    }   

}

TEST_CASE("RepresentativeTraversalFinder finds traversals correctly", "[genotype][representativetraversalfinder]") {

  SECTION("Traversal-finding should work on a substitution inside a deletion") {
    
    // Build a toy graph
    // Should be a substitution (3 or 4 between 2 and 5) nested inside a 1 to 6 deletion
    const string graph_json = R"(
    
        {
        "node": [
            {"id": 1, "sequence": "G"},
            {"id": 2, "sequence": "A"},
            {"id": 3, "sequence": "T"},
            {"id": 4, "sequence": "GGG"},
            {"id": 5, "sequence": "GT"},
            {"id": 6, "sequence": "GT"}
        ],
        "edge": [
            {"from": 1, "to": 2},
            {"from": 2, "to": 3},
            {"from": 2, "to": 4},
            {"from": 3, "to": 5},
            {"from": 4, "to": 5},
            {"from": 5, "to": 6},
            {"from": 1, "to": 6}
        ]
        }
        )";
    
    // Make an actual graph
    VG graph;
    Graph chunk;
    json2pb(chunk, graph_json.c_str(), graph_json.size());
    graph.merge(chunk);

    // Find the snarls
    CactusSnarlFinder cubs(graph);
    SnarlManager snarl_manager = cubs.find_snarls();
    const vector<const Snarl*>& snarl_roots = snarl_manager.top_level_snarls();

    // We should have a single root snarl
    REQUIRE(snarl_roots.size() == 1);
    const Snarl* parent = snarl_roots[0];
        
    REQUIRE(parent->type() == ULTRABUBBLE);

    auto children = snarl_manager.children_of(parent);
    REQUIRE(children.size() == 1);
    const Snarl* child = children[0];
        
    REQUIRE(child->type() == ULTRABUBBLE);
        
    // We need an AugmentedGraph wraping the graph to use the
    // RepresentativeTraversalFinder
    SupportAugmentedGraph augmented;
    augmented.graph = graph;
        
    // Make a RepresentativeTraversalFinder
    RepresentativeTraversalFinder finder(augmented.graph, snarl_manager, 100, 100, 1000);
        
    SECTION("there should be two traversals of the substitution") {
                        
      vector<SnarlTraversal> traversals = finder.find_traversals(*child);
            
      REQUIRE(traversals.size() == 2);
            
      SECTION("all the traversals should be size 1, with just the middle node") {
        // TODO: this will have to change when we start to support traversals of non-ultrabubbles.
        for (auto traversal : traversals) {
          REQUIRE(traversal.visit_size() == 3);
                    
          // Make sure the middle is a node and not a child snarl
          REQUIRE(traversal.visit(1).node_id() != 0);
        }
      }
            
      SECTION("one should hit node 3") {
        bool found = false;
        for (auto traversal : traversals) {
          for (size_t i = 0; i < traversal.visit_size(); i++) {
            // Check every visit on every traversal
            if (traversal.visit(i).node_id() == 3) {
              found = true;
            }
          }
        }
        REQUIRE(found);
      }
            
      SECTION("one should hit node 4") {
        bool found = false;
        for (auto traversal : traversals) {
          for (size_t i = 0; i < traversal.visit_size(); i++) {
            // Check every visit on every traversal
            if (traversal.visit(i).node_id() == 4) {
              found = true;
            }
          }
        }
        REQUIRE(found);
      }
                
    }
        
    SECTION("there should be two traversals of the parent deletion") {
            
      vector<SnarlTraversal> traversals = finder.find_traversals(*parent);

      REQUIRE(traversals.size() == 2);
            
      SECTION("one should be empty") {
        bool found = false;
        for (auto traversal : traversals) {
          if (traversal.visit_size() == 2) {
            found = true;
          }
        }
        REQUIRE(found);
      }
            
      SECTION("one should visit just the child site") {
        bool found = false;
        for (auto traversal : traversals) {
          if (traversal.visit_size() == 3) {
            auto& visit = traversal.visit(1);
                        
            if(visit.node_id() == 0 && visit.snarl().start() == child->start() &&
               visit.snarl().end() == child->end()) {
              found = true;
            }
          }
        }
        REQUIRE(found);
      }
                
    }
  }

}

TEST_CASE("RepresentativeTraversalFinder finds traversals of simple inversions", "[genotype][representativetraversalfinder]") {


    // Build a toy graph
    // Should be a substitution (3 or 4 between 2 and 5) nested inside a 1 to 6 deletion
    const string graph_json = R"(
    {
    "node": [
        {"id": 1, "sequence": "G"},
        {"id": 2, "sequence": "ATTA"},
        {"id": 3, "sequence": "CA"}
    ],
    "edge": [
        {"from": 1, "to": 2},
        {"from": 2, "to": 3},
        {"from": 1, "to": 2, "to_end": true},
        {"from": 2, "from_start": true, "to": 3}
    ]
    }
    )";

    // Make an actual graph
    VG graph;
    Graph chunk;
    json2pb(chunk, graph_json.c_str(), graph_json.size());
    graph.merge(chunk);

    // Find the snarls
    CactusSnarlFinder cubs(graph);
    SnarlManager snarl_manager = cubs.find_snarls();
    const vector<const Snarl*>& snarl_roots = snarl_manager.top_level_snarls();

    // We should have a single root snarl that is the inversion
    REQUIRE(snarl_roots.size() == 1);
    const Snarl* root = snarl_roots[0];
        
    // It should *not* be an ultrabubble
    REQUIRE(root->type() != ULTRABUBBLE);
    // But it should be connected through
    REQUIRE(root->start_end_reachable());

    // We need an AugmentedGraph wraping the graph to use the
    // RepresentativeTraversalFinder
    SupportAugmentedGraph augmented;
    augmented.graph = graph;
        
    // Make a RepresentativeTraversalFinder
    RepresentativeTraversalFinder finder(augmented.graph, snarl_manager, 100, 100, 1000);
        
    vector<SnarlTraversal> traversals = finder.find_traversals(*root);
    
    // There should be inverted and normal traversals
    REQUIRE(traversals.size() == 2);
    
    bool found_inverted = false;
    bool found_normal = false;
    
    for (auto& trav : traversals) {
        REQUIRE(trav.visit_size() == 3);
        REQUIRE(trav.visit(0).node_id() == 1);
        REQUIRE(trav.visit(1).node_id() == 2);
        REQUIRE(trav.visit(2).node_id() == 3);
        
        if (trav.visit(1).backward()) {
            found_inverted = true;
        } else {
            found_normal = true;
        }
    }
    
    REQUIRE(found_inverted);
    REQUIRE(found_normal);
}

TEST_CASE("GBWTTraversalFinder finds traversals for GBWT threads", "[genotype][gbwttraversalfinder]") {

    // (copied from haplotypes.cpp)
    
    // This graph is the start of xy2 from test/small
    string graph_json = R"({"node": [{"id": 1, "sequence": "CAAATAAGGCTT"}, {"id": 2, "sequence": "G"}, {"id": 3, "sequence": "GGAAATTTTC"}, {"id": 4, "sequence": "C"}, {"id": 5, "sequence": "TGGAGTTCTATTATATTCC"}, {"id": 6, "sequence": "G"}, {"id": 7, "sequence": "A"}, {"id": 8, "sequence": "ACTCTCTGGTTCCTG"}, {"id": 9, "sequence": "A"}, {"id": 10, "sequence": "G"}, {"id": 11, "sequence": "TGCTATGTGTAACTAGTAATGGTAATGGATATGTTGGGCTTTTTTCTTTGATTTATTTGAAGTGACGTTTGACAATCTATCACTAGGGGTAATGTGGGGAAATGGAAAGAATACAAGATTTGGAGCCA"}], "edge": [{"from": 1, "to": 2}, {"from": 1, "to": 3}, {"from": 2, "to": 3}, {"from": 3, "to": 4}, {"from": 3, "to": 5}, {"from": 4, "to": 5}, {"from": 5, "to": 6}, {"from": 5, "to": 7}, {"from": 6, "to": 8}, {"from": 7, "to": 8}, {"from": 8, "to": 9}, {"from": 8, "to": 10}, {"from": 9, "to": 11}, {"from": 10, "to": 11}]})";
  
    // Load the JSON
    vg::Graph proto_graph;
    json2pb(proto_graph, graph_json.c_str(), graph_json.size());
    // Build the xg index
    xg::XG xg_index;
    xg_index.from_path_handle_graph(vg::VG(proto_graph));
    
    gbwt::Verbosity::set(gbwt::Verbosity::SILENT);
  
    // Populate a map from node ID to encoded GBWT node.
    map<vg::id_t, gbwt::vector_type::value_type> tm;
    map<vg::id_t, gbwt::vector_type::value_type> tr;
    xg_index.for_each_handle([&](const vg::handle_t& here) {
            tm[xg_index.get_id(here)] = gbwt::Node::encode(xg_index.get_id(here), false);
            tr[xg_index.get_id(here)] = gbwt::Node::encode(xg_index.get_id(here), true);
        });
  
    // Make the threads
    gbwt::vector_type thread0 = {tm[1], tm[3], tm[4], tm[5], tm[7], tm[8], tm[10], tm[11], static_cast<gbwt::vector_type::value_type>(gbwt::ENDMARKER)};
    gbwt::vector_type thread1 = {tm[1], tm[2], tm[3], tm[4], tm[5], tm[6], tm[8], tm[9], tm[11], static_cast<gbwt::vector_type::value_type>(gbwt::ENDMARKER)};
    gbwt::vector_type thread2 = {tm[2], tm[3], tm[4], tm[5], tm[6], tm[8], tm[9], tm[11], static_cast<gbwt::vector_type::value_type>(gbwt::ENDMARKER)};
    gbwt::vector_type thread3 = {tm[1], tm[3], tm[4], tm[5], tm[6], tm[8], static_cast<gbwt::vector_type::value_type>(gbwt::ENDMARKER)};
    
    gbwt::vector_type thread1r = {tr[11], tr[9], tr[8], tr[6], tr[5], tr[4], tr[3], tr[2], tr[1], static_cast<gbwt::vector_type::value_type>(gbwt::ENDMARKER)};
  
    vector<gbwt::vector_type > haplotypes_to_add = {
        thread0,
        thread1,
        thread2,
        thread3,
        thread1r
    };

    function<SnarlTraversal(gbwt::vector_type&)> gv_to_st = [&](gbwt::vector_type& gv) {
        SnarlTraversal st;
        for (int i = 0; i < gv.size() - 1; ++i) {
            Visit* v= st.add_visit();
            v->set_node_id(gbwt::Node::id(gv[i]));
            v->set_backward(gbwt::Node::is_reverse(gv[i]));
        }
        return st;
    };

    // a fake snarl spanning all the variants
    Snarl one_eleven;
    one_eleven.mutable_start()->set_node_id(1);
    one_eleven.mutable_start()->set_backward(false);
    one_eleven.mutable_end()->set_node_id(11);
    one_eleven.mutable_end()->set_backward(false);

    Snarl eleven_one;
    eleven_one.mutable_start()->set_node_id(11);
    eleven_one.mutable_start()->set_backward(true);
    eleven_one.mutable_end()->set_node_id(1);
    eleven_one.mutable_end()->set_backward(true);

    // just a single forward thread
    gbwt::DynamicGBWT gbwt_d_1f;
    gbwt_d_1f.insert(thread1);
    gbwt::GBWT gbwt_1f(gbwt_d_1f);
    GBWTTraversalFinder trav_finder_1f(xg_index, gbwt_1f);
    vector<SnarlTraversal> travs_1f = trav_finder_1f.find_traversals(one_eleven);
    REQUIRE(travs_1f.size() == 1);
    REQUIRE(travs_1f[0].visit_size() == thread1.size() - 1);
    REQUIRE(travs_1f[0] == gv_to_st(thread1));

    vector<SnarlTraversal> travs_1f1 = trav_finder_1f.find_traversals(eleven_one);
    REQUIRE(travs_1f1.size() == 1);
    SnarlTraversal flip_1f;
    for (int i = travs_1f1[0].visit_size() - 1; i >= 0; --i) {
        Visit* visit = flip_1f.add_visit();
        *visit = travs_1f1[0].visit(i);
        visit->set_backward(!visit->backward());
    }
    REQUIRE(flip_1f == travs_1f[0]);

    // just a single reverse thread
    gbwt::DynamicGBWT gbwt_d_1r;
    gbwt_d_1r.insert(thread1);
    gbwt::GBWT gbwt_1r(gbwt_d_1r);
    GBWTTraversalFinder trav_finder_1r(xg_index, gbwt_1r);
    vector<SnarlTraversal> travs_1r = trav_finder_1r.find_traversals(one_eleven);
    REQUIRE(travs_1r.size() == 1);
    REQUIRE(travs_1r[0].visit_size() == thread1r.size() - 1);

    vector<SnarlTraversal> travs_1r1 = trav_finder_1r.find_traversals(eleven_one);
    REQUIRE(travs_1r1.size() == 1);
    REQUIRE(travs_1r1[0] == gv_to_st(thread1r));
    SnarlTraversal flip_1r;
    for (int i = travs_1r1[0].visit_size() - 1; i >= 0; --i) {
        Visit* visit = flip_1r.add_visit();
        *visit = travs_1r1[0].visit(i);
        visit->set_backward(!visit->backward());
    }
    REQUIRE(flip_1r == travs_1r[0]);

    // forward and reverse version of the same thread.  we only want one back
    gbwt::DynamicGBWT gbwt_d_1fr;
    gbwt_d_1fr.insert(thread1);
    gbwt_d_1fr.insert(thread1r);
    gbwt::GBWT gbwt_1fr(gbwt_d_1fr);
    GBWTTraversalFinder trav_finder_1fr(xg_index, gbwt_1fr);
    vector<SnarlTraversal> travs_1fr = trav_finder_1fr.find_traversals(one_eleven);
    REQUIRE(travs_1fr.size() == 1);
    REQUIRE(travs_1fr[0].visit_size() == thread1.size() - 1);
    REQUIRE(travs_1f[0] == gv_to_st(thread1));

    // all the threads
    gbwt::DynamicGBWT gbwt_d_all;
    for (auto& thread : haplotypes_to_add) {
        gbwt_d_all.insert(thread);
    }
    gbwt::GBWT gbwt_all(gbwt_d_all);
    GBWTTraversalFinder trav_finder_all(xg_index, gbwt_all);
    vector<SnarlTraversal> travs_all = trav_finder_all.find_traversals(one_eleven);
    // only 2 unique spanning traversals
    REQUIRE(travs_all.size() == 2);
    // canonicalize order for comparison
    if (travs_all[0].visit_size() > travs_all[1].visit_size()) {
        std::swap(travs_all[0], travs_all[1]);
    }
    REQUIRE(gv_to_st(thread0) == travs_all[0]);
    REQUIRE(gv_to_st(thread1) == travs_all[1]);
}

}
}
