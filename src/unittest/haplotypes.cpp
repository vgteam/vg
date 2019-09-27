//
// Tests for haplotype score functions
//

#include "catch.hpp"
#include "haplotypes.hpp"
#include "xg.hpp"
#include "vg.hpp"

#include <numeric>

namespace unittest {
using namespace std;

using thread_t = vector<gbwt::node_type>;

vg::Path path_from_thread_t(thread_t& t) {
    vg::Path toReturn;
    int rank = 1;
    for (int i = 0; i < t.size(); i++) {
        vg::Mapping* mapping = toReturn.add_mapping();
        
        // Set up the position
        mapping->mutable_position()->set_node_id(gbwt::Node::id(t[i]));
        mapping->mutable_position()->set_is_reverse(gbwt::Node::is_reverse(t[i]));
        
        // Set the rank
        mapping->set_rank(rank++);
    }
    // We're done making the path
    return toReturn;
}

TEST_CASE("We can represent appropriate graphs according to linear reference", "[slls][bubble-finding]") {
  vector<gbwt::node_type> tm = {
    gbwt::Node::encode(1, false),
    gbwt::Node::encode(1, false),
    gbwt::Node::encode(2, false),
    gbwt::Node::encode(3, false),
    gbwt::Node::encode(4, false),
    gbwt::Node::encode(5, false)
  };
  
  // POSITIVE SNVs
  
  string SNP_graph_json = R"(
  {"node":[
    {"id":1,"sequence":"AAA"},
    {"id":2,"sequence":"C"},
    {"id":3,"sequence":"T"},
    {"id":4,"sequence":"AAA"}
  ],
  "edge":[
    {"to":2,"from":1},
    {"to":3,"from":1},
    {"to":4,"from":2},
    {"to":4,"from":3}
  ],
  "path":[
    {"name":"reference","mapping":[
      {"position":{"node_id":1},"rank":1},
      {"position":{"node_id":2},"rank":2},
      {"position":{"node_id":4},"rank":3}
    ]}
  ]}
  )";
  
  thread_t SNP_thread = {tm[1], tm[3], tm[4]};
    
  string del_graph_json = R"(
  {"node":[
    {"id":1,"sequence":"AAA"},
    {"id":2,"sequence":"C"},
    {"id":3,"sequence":"T"},
    {"id":4,"sequence":"AAA"}
  ],
  "edge":[
    {"to":2,"from":1},
    {"to":3,"from":1},
    {"to":4,"from":1},
    {"to":4,"from":2},
    {"to":4,"from":3}
  ],
  "path":[
    {"name":"reference","mapping":[
      {"position":{"node_id":1},"rank":1},
      {"position":{"node_id":2},"rank":2},
      {"position":{"node_id":4},"rank":3}
    ]}
  ]}
  )";
  
  thread_t del_ref_thread = {tm[1], tm[2], tm[4]};
  thread_t del_thread = {tm[1], tm[4]};
  
  vg::Graph SNP_proto_graph;
  json2pb(SNP_proto_graph, SNP_graph_json.c_str(), SNP_graph_json.size());
  // Build the xg index
  xg::XG SNP_xg_index;
  SNP_xg_index.from_path_handle_graph(vg::VG(SNP_proto_graph));
  vg::path_handle_t SNP_ref_path_handle = SNP_xg_index.get_path_handle("reference");
  
  vg::Graph del_proto_graph;
  json2pb(del_proto_graph, del_graph_json.c_str(), del_graph_json.size());
  // Build the xg index
  xg::XG del_xg_index;
  del_xg_index.from_path_handle_graph(vg::VG(del_proto_graph));
  vg::path_handle_t del_ref_path_handle = del_xg_index.get_path_handle("reference");
  
  // NEGATIVE SNVs
  
  string long_graph_json = R"(
  {"node":[
    {"id":1,"sequence":"AAA"},
    {"id":2,"sequence":"C"},
    {"id":3,"sequence":"CT"},
    {"id":4,"sequence":"AAA"}
  ],
  "edge":[
    {"to":2,"from":1},
    {"to":3,"from":1},
    {"to":4,"from":2},
    {"to":4,"from":3}
  ],
  "path":[
    {"name":"reference","mapping":[
      {"position":{"node_id":1},"rank":1},
      {"position":{"node_id":2},"rank":2},
      {"position":{"node_id":4},"rank":3}
    ]}
  ]}
  )";
  
  thread_t long_thread = {tm[1], tm[2], tm[4]};
      
  string double_graph_json = R"(
  {"node":[
    {"id":1,"sequence":"AAA"},
    {"id":2,"sequence":"C"},
    {"id":3,"sequence":"G"},
    {"id":5,"sequence":"T"},
    {"id":4,"sequence":"AAA"}
  ],
  "edge":[
    {"to":2,"from":1},
    {"to":3,"from":1},
    {"to":4,"from":2},
    {"to":5,"from":3},
    {"to":4,"from":5}
  ],
  "path":[
    {"name":"reference","mapping":[
      {"position":{"node_id":1},"rank":1},
      {"position":{"node_id":2},"rank":2},
      {"position":{"node_id":4},"rank":3}
    ]}
  ]}
  )";
  
  thread_t double_thread = {tm[1], tm[2], tm[4]};

  vg::Graph long_proto_graph;
  json2pb(long_proto_graph, long_graph_json.c_str(), long_graph_json.size());
  // Build the xg index
  xg::XG long_xg_index;
  long_xg_index.from_path_handle_graph(vg::VG(long_proto_graph));
  vg::path_handle_t long_ref_path_handle = long_xg_index.get_path_handle("reference");
  
  vg::Graph double_proto_graph;
  json2pb(double_proto_graph, double_graph_json.c_str(), double_graph_json.size());
  // Build the xg index
  xg::XG double_xg_index;
  double_xg_index.from_path_handle_graph(vg::VG(double_proto_graph));
  vg::path_handle_t double_ref_path_handle = double_xg_index.get_path_handle("reference");

  string matching_test_file = "matching_test.slls";
  ofstream slls_out;
  slls_out.open(matching_test_file, ios::out | ios::trunc);
  slls_out << "0\t7\t1\n3\n3\n2\t1\t0\t0\t0\n2";
  slls_out.close();  
  ifstream slls_in;
  slls_in.open(matching_test_file, ios::in);
  haplo::linear_haplo_structure SNP_lin_DP(slls_in, -3, -2, SNP_xg_index, SNP_ref_path_handle);
  slls_in.close();
  slls_in.open(matching_test_file, ios::in);
  haplo::linear_haplo_structure del_lin_DP(slls_in, -3, -2, del_xg_index, del_ref_path_handle);
  slls_in.close();
  slls_in.open(matching_test_file, ios::in);
  haplo::linear_haplo_structure long_lin_DP(slls_in, -3, -2, long_xg_index, long_ref_path_handle);
  slls_in.close();
  slls_in.open(matching_test_file, ios::in);
  haplo::linear_haplo_structure double_lin_DP(slls_in, -3, -2, double_xg_index, double_ref_path_handle);
  slls_in.close();
  slls_in.open(matching_test_file, ios::in);
  haplo::linear_haplo_structure match_lin_DP(slls_in, -3, -2, SNP_xg_index, SNP_ref_path_handle);
  slls_in.close();
  remove(matching_test_file.c_str());
  
  string nonmatching_test_file = "nonmatching_test.slls";
  ofstream nslls_out;
  nslls_out.open(nonmatching_test_file, ios::out | ios::trunc);
  nslls_out << "0\t7\t1\n5\n3\n2\t1\t0\t0\t0\n2";
  nslls_out.close();  
  ifstream nslls_in;
  nslls_in.open(nonmatching_test_file, ios::in);
  haplo::linear_haplo_structure nonmatch_lin_DP(nslls_in, -3, -2, SNP_xg_index, SNP_ref_path_handle);
  nslls_in.close();
  remove(nonmatching_test_file.c_str());

  vg::Path SNP_path = path_from_thread_t(SNP_thread);
  vg::Path del_path = path_from_thread_t(del_thread);
  vg::Path del_ref_path = path_from_thread_t(del_ref_thread);
  vg::Path long_path = path_from_thread_t(long_thread);
  vg::Path double_path = path_from_thread_t(double_thread);
  SECTION( "atomic node-type checking" ) {
    REQUIRE(SNP_lin_DP.is_solitary_ref(1) == true);
    REQUIRE(del_lin_DP.is_solitary_ref(1) == true);
    REQUIRE(del_lin_DP.is_solitary_ref(4) == true);
    REQUIRE(SNP_lin_DP.is_solitary_ref(4) == true);
    REQUIRE(SNP_lin_DP.is_solitary_ref(2) == false);
    REQUIRE(SNP_lin_DP.is_snv(2) == true);
    REQUIRE(SNP_lin_DP.is_snv(3) == true);
    REQUIRE(SNP_lin_DP.is_snv(3) == true);
    REQUIRE(long_lin_DP.is_snv(3) == false);
    REQUIRE(double_lin_DP.is_snv(3) == false);
    REQUIRE(SNP_lin_DP.get_type(1) == haplo::linear_haplo_structure::ref_span);
    REQUIRE(SNP_lin_DP.get_type(2) == haplo::linear_haplo_structure::snv);
    REQUIRE(long_lin_DP.get_type(3) == haplo::linear_haplo_structure::invalid);
    REQUIRE(double_lin_DP.get_type(3) == haplo::linear_haplo_structure::invalid);
  }
  SECTION( "graph with SNP" ) {
    bool passed = true;
    haplo::linear_haplo_structure::SNVvector SNP_snvs = SNP_lin_DP.SNVs(SNP_path);
    try{
      haplo::linear_haplo_structure::SNVvector SNP_snvs = SNP_lin_DP.SNVs(SNP_path);
    } catch(haplo::linear_haplo_structure::linearUnrepresentable& e) {
      passed = false;
    }
    REQUIRE(passed);
  }
  SECTION( "graph with deletion" ) {
    bool passed = true;
    haplo::linear_haplo_structure::SNVvector del_snvs = del_lin_DP.SNVs(del_path);
    try{
      haplo::linear_haplo_structure::SNVvector SNP_snvs = SNP_lin_DP.SNVs(SNP_path);
    } catch(haplo::linear_haplo_structure::linearUnrepresentable& e) {
      passed = false;
    }
    REQUIRE(passed);
  }
  SECTION( "graph with too long polymorphism" ) {
    bool passed = true;
    try{
      haplo::linear_haplo_structure::SNVvector long_snvs = long_lin_DP.SNVs(long_path);
    } catch(haplo::linear_haplo_structure::linearUnrepresentable& e) {
      passed = false;
    }
    REQUIRE(!passed);
  }
  SECTION( "graph with 2-node path in superbubble" ) {
    bool passed = true;
    try{
      haplo::linear_haplo_structure::SNVvector double_snvs = double_lin_DP.SNVs(double_path);
    } catch(haplo::linear_haplo_structure::linearUnrepresentable& e) {
      passed = false;
    }
    REQUIRE(!passed);
  }
}
//TODO test paths that aren't in the index

TEST_CASE("We can score haplotypes using GBWT", "[haplo-score][gbwt]") {

  gbwt::Verbosity::set(gbwt::Verbosity::SILENT);
  gbwt::DynamicGBWT* gbwt_index = new gbwt::DynamicGBWT;
  
  gbwt::vector_type tm = {
    static_cast<gbwt::vector_type::value_type>(gbwt::Node::encode(1, false)), // dummy for 1-based indexing
    static_cast<gbwt::vector_type::value_type>(gbwt::Node::encode(1, false)),
    static_cast<gbwt::vector_type::value_type>(gbwt::Node::encode(2, false)),
    static_cast<gbwt::vector_type::value_type>(gbwt::Node::encode(3, false)),
    static_cast<gbwt::vector_type::value_type>(gbwt::Node::encode(4, false)),
    static_cast<gbwt::vector_type::value_type>(gbwt::Node::encode(5, false)),
    static_cast<gbwt::vector_type::value_type>(gbwt::Node::encode(6, false)),
    static_cast<gbwt::vector_type::value_type>(gbwt::Node::encode(7, false)),
    static_cast<gbwt::vector_type::value_type>(gbwt::Node::encode(8, false)),
  };
  
  gbwt::vector_type t1_2_4_5_7 = {tm[1], tm[2], tm[4], tm[5], tm[7], static_cast<gbwt::vector_type::value_type>(gbwt::ENDMARKER)};
  gbwt::vector_type t1_3_4_5_7 = {tm[1], tm[3], tm[4], tm[5], tm[7], static_cast<gbwt::vector_type::value_type>(gbwt::ENDMARKER)};
  gbwt::vector_type t1_2_4_6_7 = {tm[1], tm[2], tm[4], tm[6], tm[7], static_cast<gbwt::vector_type::value_type>(gbwt::ENDMARKER)};
  gbwt::vector_type t1_3_4_6_7 = {tm[1], tm[3], tm[4], tm[6], tm[7], static_cast<gbwt::vector_type::value_type>(gbwt::ENDMARKER)};
  gbwt::vector_type t1_2_4_d_7 = {tm[1], tm[2], tm[4], tm[7], static_cast<gbwt::vector_type::value_type>(gbwt::ENDMARKER)};
  gbwt::vector_type t1_3_4_d_7 = {tm[1], tm[3], tm[4], tm[7], static_cast<gbwt::vector_type::value_type>(gbwt::ENDMARKER)};
  
  vector<gbwt::vector_type> haplotypes_to_add = {
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
  
  gbwt::vector_type query_nodes = {tm[1], tm[2], tm[4]};
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
  xg::XG xg_index;
  xg_index.from_path_handle_graph(vg::VG(proto_graph));
    
  gbwt::Verbosity::set(gbwt::Verbosity::SILENT);
  gbwt::DynamicGBWT* gbwt_index = new gbwt::DynamicGBWT;
  
  // Populate a map from node ID to encoded GBWT node.
  map<vg::id_t, gbwt::vector_type::value_type> tm;
  // And a map of node lengths
  map<vg::id_t, size_t> node_lengths;
  xg_index.for_each_handle([&](const vg::handle_t& here) {
    tm[xg_index.get_id(here)] = gbwt::Node::encode(xg_index.get_id(here), false);
    node_lengths[xg_index.get_id(here)] = xg_index.get_length(here);
  });
  
  // Make the two threads
  gbwt::vector_type thread0 = {tm[1], tm[3], tm[4], tm[5], tm[7], tm[8], tm[10], tm[11], static_cast<gbwt::vector_type::value_type>(gbwt::ENDMARKER)};
  gbwt::vector_type thread1 = {tm[1], tm[2], tm[3], tm[4], tm[5], tm[6], tm[8], tm[9], tm[11], static_cast<gbwt::vector_type::value_type>(gbwt::ENDMARKER)};
  
  vector<gbwt::vector_type > haplotypes_to_add = {
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
  gbwt::vector_type should_match = {tm[1], tm[3], tm[4], tm[5], tm[7], tm[8], tm[10], tm[11]};
  //And one that should need a crossover
  gbwt::vector_type should_crossover = {tm[1], tm[3], tm[4], tm[5], tm[7], tm[8], tm[9], tm[11]};
  
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
