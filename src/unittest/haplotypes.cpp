//
// Tests for haplotype score functions
//

#include "catch.hpp"
#include "haplotypes.hpp"

#include <numeric>

namespace unittest {
using namespace std;

xg::XG::ThreadMapping make_xg_mapping(int64_t node_id) {
  xg::XG::ThreadMapping to_return;
  to_return.node_id = node_id;
  to_return.is_reverse = false;
  return to_return;
}

//TODO test paths that aren't in the index

TEST_CASE("We can score haplotypes using gPBWT", "[haplo-score][xg]") {
  //     / [3] \   / [6] \         ###################
  //  [1]       [4]-------[7]      # baby test graph #
  //    \\ [2] /   \ [5] /         ###################
  //    \__[8]
  string graph_json = R"(
  {"node":[{"id":1,"sequence":"GATT"},
  {"id":2,"sequence":"ACA"},
  {"id":3,"sequence":"TT"},
  {"id":4,"sequence":"ACA"},
  {"id":5,"sequence":"TTAG"},
  {"id":6,"sequence":"GG"},
  {"id":7,"sequence":"ATTACA"},
  {"id":8,"sequence":"AAAA"}],
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
    make_xg_mapping(7),
    make_xg_mapping(8)
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
  haplotypes_to_add.clear();
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
  
  vector<xg::XG::ThreadMapping> missing_edge = {tm[1], tm[4]};  
  REQUIRE(haplo::haplo_DP::score(missing_edge, xg_index, memo).second);
  
  vector<xg::XG::ThreadMapping> empty_node = {tm[1], tm[8]};
  REQUIRE(!(haplo::haplo_DP::score(empty_node, xg_index, memo).second));
}

TEST_CASE("We can score haplotypes using GBWT", "[haplo-score][gbwt]") {

  gbwt::Verbosity::set(gbwt::Verbosity::SILENT);
  gbwt::DynamicGBWT* gbwt_index = new gbwt::DynamicGBWT;
  
  vector<gbwt::node_type> tm = {
    gbwt::Node::encode(1, false), // dummy for 1-based indexing
    gbwt::Node::encode(1, false),
    gbwt::Node::encode(2, false),
    gbwt::Node::encode(3, false),
    gbwt::Node::encode(4, false),
    gbwt::Node::encode(5, false),
    gbwt::Node::encode(6, false),
    gbwt::Node::encode(7, false),
    gbwt::Node::encode(8, false),
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
    6,
    3
  };
  
  for(size_t i = 0; i < haplotypes_to_add.size(); i++) {
    gbwt_index->insert(haplotypes_to_add[i]);
  }
  
  REQUIRE(gbwt_index->nodeSize(tm[1]) == 12);
  REQUIRE(gbwt_index->nodeSize(tm[2]) == 8);
  REQUIRE(gbwt_index->nodeSize(tm[1]) == 12);
  
  haplo::haploMath::RRMemo memo(9, 12);
  
  // initial node
  haplo::hDP_gbwt_graph_accessor<gbwt::DynamicGBWT> ga1(*gbwt_index, tm[1], node_lengths[1], memo);
  haplo::haplo_DP_column hdpc(ga1);
  REQUIRE(hdpc.get_sizes()[0] == 12);
  REQUIRE(hdpc.current_sum() <= 0);
  double last_sum = hdpc.current_sum();
  
  // second node
  haplo::hDP_gbwt_graph_accessor<gbwt::DynamicGBWT> ga12(*gbwt_index, tm[1], tm[2], node_lengths[2], memo);
  hdpc.extend(ga12);
  REQUIRE(hdpc.get_sizes().size() == 1);
  REQUIRE(hdpc.current_sum() < last_sum);
  REQUIRE(hdpc.get_sizes()[0] == 8);
  last_sum = hdpc.current_sum();
  
  // third node
  haplo::hDP_gbwt_graph_accessor<gbwt::DynamicGBWT> ga24(*gbwt_index, tm[2], tm[4], node_lengths[4], memo);
  hdpc.extend(ga24);
  REQUIRE(hdpc.get_sizes().size() == 2);
  REQUIRE(hdpc.get_sizes()[0] == 4);
  REQUIRE(hdpc.get_sizes()[1] == 8);
  REQUIRE(hdpc.current_sum() > last_sum);
  
  last_sum = hdpc.current_sum();
  
  vector<gbwt::node_type> query_nodes = {tm[1], tm[2], tm[4]};
  vector<size_t> query_node_lengths = {node_lengths[1], node_lengths[2], node_lengths[4]};
  haplo::gbwt_thread_t query(query_nodes, query_node_lengths);
  
  pair<double, bool> result_from_thread = haplo::haplo_DP::score(query, *gbwt_index, memo);  
  REQUIRE(fabs(result_from_thread.first - last_sum) < 0.000001);
  
  // recognize nonexistent edge
  haplo::hDP_gbwt_graph_accessor<gbwt::DynamicGBWT> ga14(*gbwt_index, tm[1], tm[4], 3, memo);
  REQUIRE(!ga14.has_edge());
  
  query_nodes = {tm[1], tm[4]};
  query_node_lengths = {node_lengths[1], node_lengths[4]};
  haplo::gbwt_thread_t missing_edge(query_nodes, query_node_lengths);
  
  REQUIRE(haplo::haplo_DP::score(missing_edge, *gbwt_index, memo).second);
  
  query_nodes = {tm[1], tm[8]};
  query_node_lengths = {node_lengths[1], node_lengths[8]};
  haplo::gbwt_thread_t empty_node(query_nodes, query_node_lengths);
  REQUIRE(!(haplo::haplo_DP::score(empty_node, *gbwt_index, memo).second));
  delete gbwt_index;
}

TEST_CASE("We can recognize a required crossover", "[hapo-score][gbwt]") {
  // This graph is the start of xy2 from test/small
  string graph_json = R"({"node": [{"id": 1, "sequence": "CAAATAAGGCTT"}, {"id": 2, "sequence": "G"}, {"id": 3, "sequence": "GGAAATTTTC"}, {"id": 4, "sequence": "C"}, {"id": 5, "sequence": "TGGAGTTCTATTATATTCC"}, {"id": 6, "sequence": "G"}, {"id": 7, "sequence": "A"}, {"id": 8, "sequence": "ACTCTCTGGTTCCTG"}, {"id": 9, "sequence": "A"}, {"id": 10, "sequence": "G"}, {"id": 11, "sequence": "TGCTATGTGTAACTAGTAATGGTAATGGATATGTTGGGCTTTTTTCTTTGATTTATTTGAAGTGACGTTTGACAATCTATCACTAGGGGTAATGTGGGGAAATGGAAAGAATACAAGATTTGGAGCCA"}], "edge": [{"from": 1, "to": 2}, {"from": 1, "to": 3}, {"from": 2, "to": 3}, {"from": 3, "to": 4}, {"from": 3, "to": 5}, {"from": 4, "to": 5}, {"from": 5, "to": 6}, {"from": 5, "to": 7}, {"from": 6, "to": 8}, {"from": 7, "to": 8}, {"from": 8, "to": 9}, {"from": 8, "to": 10}, {"from": 9, "to": 11}, {"from": 10, "to": 11}]})";
  
  // Load the JSON
  vg::Graph proto_graph;
  json2pb(proto_graph, graph_json.c_str(), graph_json.size());
  // Build the xg index
  xg::XG xg_index(proto_graph);
    
  gbwt::Verbosity::set(gbwt::Verbosity::SILENT);
  gbwt::DynamicGBWT* gbwt_index = new gbwt::DynamicGBWT;
  
  // Populate a map from node ID to encoded GBWT node.
  map<vg::id_t, gbwt::node_type> tm;
  // And a map of node lengths
  map<vg::id_t, size_t> node_lengths;
  xg_index.for_each_handle([&](const vg::handle_t& here) {
    tm[xg_index.get_id(here)] = gbwt::Node::encode(xg_index.get_id(here), false);
    node_lengths[xg_index.get_id(here)] = xg_index.get_length(here);
  });
  
  // Make the two threads
  vector<gbwt::node_type> thread0 = {tm[1], tm[3], tm[4], tm[5], tm[7], tm[8], tm[10], tm[11], gbwt::ENDMARKER};
  vector<gbwt::node_type> thread1 = {tm[1], tm[2], tm[3], tm[4], tm[5], tm[6], tm[8], tm[9], tm[11], gbwt::ENDMARKER};
  
  vector<vector<gbwt::node_type> > haplotypes_to_add = {
    thread0,
    thread1
  };
  
  for(size_t i = 0; i < haplotypes_to_add.size(); i++) {
    gbwt_index->insert(haplotypes_to_add[i]);
  }
  
  REQUIRE(gbwt_index->nodeSize(tm[1]) == 2);
  REQUIRE(gbwt_index->nodeSize(tm[7]) == 1);
  REQUIRE(gbwt_index->nodeSize(tm[10]) == 1);
  
  haplo::haploMath::RRMemo memo(9, 2);
  
  // Now we trace a haplotype that should match
  vector<gbwt::node_type> should_match = {tm[1], tm[3], tm[4], tm[5], tm[7], tm[8], tm[10], tm[11]};
  //And one that should need a crossover
  vector<gbwt::node_type> should_crossover = {tm[1], tm[3], tm[4], tm[5], tm[7], tm[8], tm[9], tm[11]};
  
  // initial node (same for both)
  size_t i = 0;
  auto here = should_match[i];
  haplo::hDP_gbwt_graph_accessor<gbwt::DynamicGBWT> ga(*gbwt_index, here, node_lengths[gbwt::Node::id(here)], memo);
  // We have to run two columns in parallel since there is no deep copy and copies share mutable state
  haplo::haplo_DP_column hdpc0(ga);
  haplo::haplo_DP_column hdpc1(ga);
  auto sizes = hdpc0.get_sizes();
  REQUIRE(sizes.size() == 1);
  REQUIRE(accumulate(sizes.begin(), sizes.end(), 0) == 2);
  REQUIRE(hdpc0.current_sum() <= 0);
  
  
  // Run them out until they differ
  for (i++; should_match[i] == should_crossover[i]; i++) {
    auto last = here;
    here = should_match[i];
    haplo::hDP_gbwt_graph_accessor<gbwt::DynamicGBWT> ga(*gbwt_index, last, here, node_lengths[gbwt::Node::id(here)], memo);
    hdpc0.extend(ga);
    hdpc1.extend(ga);
    auto sizes = hdpc0.get_sizes();
    // In some iterations we will match two haplotypes but in some we will match only one.
    REQUIRE(sizes.size() <= 2);
    REQUIRE(accumulate(sizes.begin(), sizes.end(), 0) <= 2);
  }
  
  // Now do different things for each
  auto last = here;
  
  // The first one should match something
  here = should_match[i];
  haplo::hDP_gbwt_graph_accessor<gbwt::DynamicGBWT> ga0(*gbwt_index, last, here, node_lengths[gbwt::Node::id(here)], memo);
  hdpc0.extend(ga0);
  sizes = hdpc0.get_sizes();
  REQUIRE(sizes.size() == 1);
  REQUIRE(accumulate(sizes.begin(), sizes.end(), 0) == 1);
  auto last_sum0 = hdpc0.current_sum();
  
  // The second one should match nothing and force a crossover
  here = should_crossover[i];
  haplo::hDP_gbwt_graph_accessor<gbwt::DynamicGBWT> ga1(*gbwt_index, last, here, node_lengths[gbwt::Node::id(here)], memo);
  hdpc1.extend(ga1);
  sizes = hdpc1.get_sizes();
  REQUIRE(sizes.size() == 1);
  REQUIRE(accumulate(sizes.begin(), sizes.end(), 0) == 1);
  auto last_sum1 = hdpc1.current_sum();
  
  // Crossing over has to be worse
  REQUIRE(last_sum1 < last_sum0);
    
    
}

}
