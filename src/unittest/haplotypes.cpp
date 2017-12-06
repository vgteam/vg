//
// Tests for haplotype score functions
//

#include "catch.hpp"
#include "haplotypes.hpp"

namespace unittest {
using namespace std;

xg::XG::ThreadMapping make_xg_mapping(int64_t node_id) {
  xg::XG::ThreadMapping to_return;
  to_return.node_id = node_id;
  to_return.is_reverse = false;
  return to_return;
}

TEST_CASE("We can score haplotypes using gPBWT", "[haplo-score][xg]") {
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
  
  vector<xg::XG::ThreadMapping> t1_2_4_5_7 = {tm[1], tm[2], tm[4], tm[5], tm[7]};
  vector<xg::XG::ThreadMapping> t1_3_4_5_7 = {tm[1], tm[3], tm[4], tm[5], tm[7]};
  vector<xg::XG::ThreadMapping> t1_2_4_6_7 = {tm[1], tm[2], tm[4], tm[6], tm[7]};
  vector<xg::XG::ThreadMapping> t1_3_4_6_7 = {tm[1], tm[3], tm[4], tm[6], tm[7]};
  vector<xg::XG::ThreadMapping> t1_2_4_d_7 = {tm[1], tm[2], tm[4], tm[7]};
  vector<xg::XG::ThreadMapping> t1_3_4_d_7 = {tm[1], tm[3], tm[4], tm[7]};
  
  vector<xg::XG::ThreadMapping> query = {tm[1], tm[2], tm[4]};
  
  vector<vector<xg::XG::ThreadMapping> > haplotypes_to_add = {
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
  vg::Graph proto_graph;
  json2pb(proto_graph, graph_json.c_str(), graph_json.size());
  // Build the xg index
  xg::XG xg_index(proto_graph);
  xg_index.insert_threads_into_dag(haplotypes_to_add, haplotype_names);
  // build a score-penalty memo: recombination penalty 9, population size 12
  haplo::haploMath::RRMemo memo(9, 12);
  
  // initial node
  haplo::hDP_graph_accessor ga1(xg_index, tm[1], memo);
  haplo::haplo_DP_column hdpc(ga1);
  REQUIRE(hdpc.get_sizes()[0] == 12);
  REQUIRE(hdpc.current_sum() <= 0);
  double last_sum = hdpc.current_sum();
  
  // second node
  haplo::hDP_graph_accessor ga12(xg_index, tm[1], tm[2], memo);
  hdpc.extend(ga12);
  REQUIRE(hdpc.get_sizes().size() == 1);
  REQUIRE(hdpc.current_sum() < last_sum);
  REQUIRE(hdpc.get_sizes()[0] == 8);
  last_sum = hdpc.current_sum();
  
  // third node
  haplo::hDP_graph_accessor ga24(xg_index, tm[2], tm[4], memo);
  hdpc.extend(ga24);
  REQUIRE(hdpc.get_sizes().size() == 2);
  REQUIRE(hdpc.get_sizes()[0] == 4);
  REQUIRE(hdpc.get_sizes()[1] == 8);
  REQUIRE(hdpc.current_sum() > last_sum);
  
  last_sum = hdpc.current_sum();
  
  pair<double, bool> result_from_thread = haplo::haplo_DP::score(query, xg_index, memo);  
  REQUIRE(result_from_thread.second);
  REQUIRE(fabs(result_from_thread.first - last_sum) < 0.000001);
  
  // recognize nonexistent edge
  haplo::hDP_graph_accessor ga14(xg_index, tm[1], tm[4], memo);
  REQUIRE(!ga14.has_edge());
  
  vector<xg::XG::ThreadMapping> bad_query = {tm[1], tm[4], tm[2]};
  REQUIRE(!(haplo::haplo_DP::score(bad_query, xg_index, memo).second));
}

TEST_CASE("We can score haplotypes using GBWT", "[haplo-score][gbwt]") {
  gbwt::DynamicGBWT gbwt_index;
  
  vector<gbwt::node_type> tm = {
    gbwt::Node::encode(1, false), // dummy for 1-based indexing
    gbwt::Node::encode(1, false),
    gbwt::Node::encode(2, false),
    gbwt::Node::encode(3, false),
    gbwt::Node::encode(4, false),
    gbwt::Node::encode(5, false),
    gbwt::Node::encode(6, false),
    gbwt::Node::encode(7, false)
  };
  
  vector<gbwt::node_type> t1_2_4_5_7 = {tm[1], tm[2], tm[4], tm[5], tm[7], gbwt::ENDMARKER};
  vector<gbwt::node_type> t1_3_4_5_7 = {tm[1], tm[3], tm[4], tm[5], tm[7], gbwt::ENDMARKER};
  vector<gbwt::node_type> t1_2_4_6_7 = {tm[1], tm[2], tm[4], tm[6], tm[7], gbwt::ENDMARKER};
  vector<gbwt::node_type> t1_3_4_6_7 = {tm[1], tm[3], tm[4], tm[6], tm[7], gbwt::ENDMARKER};
  vector<gbwt::node_type> t1_2_4_d_7 = {tm[1], tm[2], tm[4], tm[7], gbwt::ENDMARKER};
  vector<gbwt::node_type> t1_3_4_d_7 = {tm[1], tm[3], tm[4], tm[7], gbwt::ENDMARKER};
  
  vector<vector<gbwt::node_type> > haplotypes_to_add = {
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
  
  vector<size_t> node_lengths = {
    0, // dummy for 1-based indexing
    4,
    3,
    2,
    3,
    4,
    2,
    6
  };
  
  for(size_t i = 0; i < haplotypes_to_add.size(); i++) {
    gbwt_index.insert(haplotypes_to_add[i]);
  }
  
  REQUIRE(gbwt_index.nodeSize(tm[1]) == 12);
  REQUIRE(gbwt_index.nodeSize(tm[2]) == 8);
  REQUIRE(gbwt_index.nodeSize(tm[1]) == 12);
  
  haplo::haploMath::RRMemo memo(9, 12);
  
  // initial node
  haplo::hDP_gbwt_graph_accessor<gbwt::DynamicGBWT> ga1(gbwt_index, tm[1], node_lengths[1], memo);
  haplo::haplo_DP_column hdpc(ga1);
  REQUIRE(hdpc.get_sizes()[0] == 12);
  REQUIRE(hdpc.current_sum() <= 0);
  double last_sum = hdpc.current_sum();
  
  // second node
  haplo::hDP_gbwt_graph_accessor<gbwt::DynamicGBWT> ga12(gbwt_index, tm[1], tm[2], node_lengths[2], memo);
  hdpc.extend(ga12);
  REQUIRE(hdpc.get_sizes().size() == 1);
  REQUIRE(hdpc.current_sum() < last_sum);
  REQUIRE(hdpc.get_sizes()[0] == 8);
  last_sum = hdpc.current_sum();

  // third node
  haplo::hDP_gbwt_graph_accessor<gbwt::DynamicGBWT> ga24(gbwt_index, tm[2], tm[4], node_lengths[4], memo);
  hdpc.extend(ga24);
  REQUIRE(hdpc.get_sizes().size() == 2);
  REQUIRE(hdpc.get_sizes()[0] == 4);
  REQUIRE(hdpc.get_sizes()[1] == 8);
  REQUIRE(hdpc.current_sum() > last_sum);
  
  last_sum = hdpc.current_sum();
  
  vector<gbwt::node_type> query_nodes = {tm[1], tm[2], tm[4]};
  vector<size_t> query_node_lengths = {node_lengths[1], node_lengths[2], node_lengths[4]};
  haplo::gbwt_thread_t query(query_nodes, query_node_lengths);
  
  pair<double, bool> result_from_thread = haplo::haplo_DP::score(query, gbwt_index, memo);  
  REQUIRE(fabs(result_from_thread.first - last_sum) < 0.000001);
}
}