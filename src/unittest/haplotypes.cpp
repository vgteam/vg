//
//  xg.cpp
//  
// Tests for haplotype score functions
//

#include "catch.hpp"
#include "haplotypes.hpp"
#include <iostream>

namespace vg {
    namespace unittest {
        using namespace std;

xg::XG::ThreadMapping make_xg_mapping(int64_t node_id) {
  xg::XG::ThreadMapping to_return;
  to_return.node_id = node_id;
  to_return.is_reverse = false;
  return to_return;
}

TEST_CASE("We can score haplotypes", "[haplo-score]") {
  //     / [3] \   / [6] \         ###################
  //  [1]       [4]-------[7]      # baby test graph #
  //     \ [2] /   \ [5] /         ###################
  string graph_json = R"(
  {"node":[{"id":1,"sequence":"GATT"},
  {"id":2,"sequence":"ACA"},
  {"id":3,"sequence":"TT"},
  {"id":4,"sequence":"ACA"},
  {"id":5,"sequence":"TTAG"},
  {"id":6,"sequence":"GG"},
  {"id":7,"sequence":"ATTACA"}],
  "edge":[{"to":2,"from":1},
  {"to":3,"from":1},
  {"to":4,"from":2},
  {"to":4,"from":3},
  {"to":5,"from":4},
  {"to":6,"from":4},
  {"to":7,"from":4},
  {"to":7,"from":5},
  {"to":7,"from":6}
  ]}
  )";
  
  vector<xg::XG::ThreadMapping> tm = {
    make_xg_mapping(1), // dummy for 1-based indexing
    make_xg_mapping(1),
    make_xg_mapping(2),
    make_xg_mapping(3),
    make_xg_mapping(4),
    make_xg_mapping(5),
    make_xg_mapping(6),
    make_xg_mapping(7)
  };
  
  thread_t t1_2_4_5_7 = {tm[1], tm[2], tm[4], tm[5], tm[7]};
  thread_t t1_3_4_5_7 = {tm[1], tm[3], tm[4], tm[5], tm[7]};
  thread_t t1_2_4_6_7 = {tm[1], tm[2], tm[4], tm[6], tm[7]};
  thread_t t1_3_4_6_7 = {tm[1], tm[3], tm[4], tm[6], tm[7]};
  thread_t t1_2_4_d_7 = {tm[1], tm[2], tm[4], tm[7]};
  thread_t t1_3_4_d_7 = {tm[1], tm[3], tm[4], tm[7]};
  
  thread_t query = {tm[1], tm[2], tm[4]};
  
  vector<thread_t> haplotypes_to_add = {
    t1_2_4_5_7,
    t1_2_4_5_7,
    t1_2_4_5_7,
    t1_2_4_5_7,
    t1_2_4_5_7,
    t1_3_4_5_7,
    t1_3_4_5_7,
    t1_2_4_5_7,
    t1_2_4_d_7,
    t1_2_4_d_7,
    t1_3_4_6_7,
    t1_3_4_6_7
  };
  
  vector<string> haplotype_names = {
    "t0",
    "t1",
    "t2",
    "t3",
    "t4",
    "t5",
    "t6",
    "t7",
    "t8",
    "t9",
    "t10",
    "t11"
  };
      
  // Load the JSON
  Graph proto_graph;
  json2pb(proto_graph, graph_json.c_str(), graph_json.size());
  // Build the xg index
  xg::XG xg_index(proto_graph);
  xg_index.insert_threads_into_dag(haplotypes_to_add, haplotype_names);
  // build a score-penalty memo: recombination penalty 9, population size 12
  RRMemo memo(9, 12);
  
  // initial node
  hDP_graph_accessor ga1(xg_index, tm[1], memo);
  haplo_DP_column hdpc(ga1);
  REQUIRE(hdpc.get_sizes()[0] == 12);
  REQUIRE(hdpc.current_sum() <= 0);
  double last_sum = hdpc.current_sum();
  
  // second node
  hDP_graph_accessor ga12(xg_index, tm[1], tm[2], memo);
  hdpc.extend(ga12);
  REQUIRE(hdpc.get_sizes().size() == 1);
  REQUIRE(hdpc.current_sum() < last_sum);
  REQUIRE(hdpc.get_sizes()[0] == 8);
  last_sum = hdpc.current_sum();
  
  // third node
  hDP_graph_accessor ga24(xg_index, tm[2], tm[4], memo);
  hdpc.extend(ga24);
  REQUIRE(hdpc.get_sizes().size() == 2);
  REQUIRE(hdpc.get_sizes()[0] == 4);
  REQUIRE(hdpc.get_sizes()[1] == 8);
  REQUIRE(hdpc.current_sum() > last_sum);
  
  last_sum = hdpc.current_sum();

  REQUIRE(haplo_DP::score(query, xg_index, memo).second);
  REQUIRE(fabs(haplo_DP::score(query, xg_index, memo).first - last_sum) < 0.000001);
  
  // recognize nonexistent edge
  hDP_graph_accessor ga14(xg_index, tm[1], tm[4], memo);
  REQUIRE(!ga14.has_edge());
  
  thread_t bad_query = {tm[1], tm[4], tm[2]};
  REQUIRE(!(haplo_DP::score(bad_query, xg_index, memo).second));
}
}
}