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

using thread_t = vector<xg::XG::ThreadMapping>;

vg::Path path_from_thread_t(thread_t& t) {
	vg::Path toReturn;
	int rank = 1;
	for(int i = 0; i < t.size(); i++) {
		vg::Mapping* mapping = toReturn.add_mapping();

    // Set up the position
    mapping->mutable_position()->set_node_id(t[i].node_id);
    mapping->mutable_position()->set_is_reverse(t[i].is_reverse);

    // Set the rank
    mapping->set_rank(rank++);
  }
  // We're done making the path
  return toReturn;
}

TEST_CASE("We can represent appropriate graphs according to linear reference", "[slls][bubble-finding]") {

  string graph_json = R"(
  {"node":[{"id":1,"sequence":"AAA"},
  {"id":2,"sequence":"C"},
  {"id":3,"sequence":"T"},
  {"id":4,"sequence":"AAA"},
  {"id":5,"sequence":"C"},
  {"id":6,"sequence":"C"},
  {"id":7,"sequence":"T"},
  {"id":8,"sequence":"AAA"},
  {"id":9,"sequence":"G"},
  {"id":10,"sequence":"AAA"},
  {"id":11,"sequence":"C"},
  {"id":12,"sequence":"T"},
  {"id":13,"sequence":"AAA"},
  {"id":14,"sequence":"C"},
  {"id":15,"sequence":"G"},
  {"id":16,"sequence":"TG"},
  {"id":17,"sequence":"AAA"}],
  "edge":[{"to":2,"from":1},
  {"to":3,"from":1},
  {"to":4,"from":2},
  {"to":4,"from":3},
  {"to":5,"from":4},
  {"to":6,"from":4},
  {"to":7,"from":6},
  {"to":8,"from":5},
  {"to":8,"from":7},
  {"to":9,"from":8},
  {"to":10,"from":9},
  {"to":10,"from":8},
  {"to":11,"from":10},
  {"to":12,"from":10},
  {"to":13,"from":10},
  {"to":13,"from":11},
  {"to":13,"from":12},
  {"to":14,"from":13},
  {"to":15,"from":13},
  {"to":16,"from":13},
  {"to":17,"from":14},
  {"to":17,"from":15},
  {"to":17,"from":16}],
  "path":[{"name":"0","mapping":[{"position":{"node_id":1},"rank":1},
  {"position":{"node_id":2},"rank":2},
  {"position":{"node_id":4},"rank":3},
  {"position":{"node_id":5},"rank":4},
  {"position":{"node_id":8},"rank":5},
  {"position":{"node_id":9},"rank":6},
  {"position":{"node_id":10},"rank":7},
  {"position":{"node_id":11},"rank":8},
  {"position":{"node_id":13},"rank":9},
  {"position":{"node_id":15},"rank":10},
  {"position":{"node_id":17},"rank":11}]}]})";

  vector<xg::XG::ThreadMapping> tm = {
    make_xg_mapping(1),
    make_xg_mapping(1),
    make_xg_mapping(2),
    make_xg_mapping(3),
    make_xg_mapping(4),
    make_xg_mapping(5),
    make_xg_mapping(6),
    make_xg_mapping(7),
    make_xg_mapping(8),
    make_xg_mapping(9),
    make_xg_mapping(10),
    make_xg_mapping(11),
    make_xg_mapping(12),
    make_xg_mapping(13),
    make_xg_mapping(14),
    make_xg_mapping(15),
    make_xg_mapping(16),
    make_xg_mapping(17)
  };
  
  vg::Graph proto_graph;
  json2pb(proto_graph, graph_json.c_str(), graph_json.size());
  // Build the xg index
  xg::XG xg_index(proto_graph);

  thread_t test_thread = {tm[1], tm[2], tm[4], tm[5], tm[8], tm[9], tm[10], tm[13], tm[15], tm[17]};
  vg::Path test_path = path_from_thread_t(test_thread);

  string disjoint_test_file = "disjoint_test.slls";
  ofstream dj_slls_out;
  dj_slls_out.open(disjoint_test_file, ios::out | ios::trunc);
  dj_slls_out << "0\t23\t2\n7\t19\n3\n0\t0\t0\t3\t0\n0\t1\t2\n0\t0\t0\t3\t0\n0\t1\t2";
  dj_slls_out.close();  
  ifstream dj_slls_in;
  dj_slls_in.open(disjoint_test_file, ios::in);
  haplo::linear_haplo_DP disjoint_lin_DP(dj_slls_in, -3, -2, xg_index, 0);
  dj_slls_in.close();
  remove(disjoint_test_file.c_str());

  string matching_test_file = "matching_test.slls";
  ofstream slls_out;
  slls_out.open(matching_test_file, ios::out | ios::trunc);
  slls_out << "0\t23\t2\n7\t19\n3\n0\t0\t0\t3\t0\n0\t1\t2\n0\t0\t0\t3\t0\n0\t1\t2";
  slls_out.close();  
  ifstream slls_in;
  slls_in.open(matching_test_file, ios::in);
  haplo::linear_haplo_DP matching_lin_DP(slls_in, -3, -2, xg_index, 0);
  slls_in.close();
  remove(matching_test_file.c_str());

  vector<int64_t> found_potential_snps = disjoint_lin_DP.potential_snps(test_path);
  vector<int64_t> found_potential_deletions = disjoint_lin_DP.potential_deletions(test_path);

  REQUIRE(disjoint_lin_DP.is_snv(2));
  REQUIRE(disjoint_lin_DP.is_snv(12));
  REQUIRE(disjoint_lin_DP.is_snv(9));
  REQUIRE(!disjoint_lin_DP.is_snv(5));
  REQUIRE(!disjoint_lin_DP.is_snv(15));
  REQUIRE(!disjoint_lin_DP.path_to_input_haplotype(test_path)->is_valid());
  REQUIRE(!disjoint_lin_DP.path_to_input_haplotype(test_path)->is_valid());
}

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
}