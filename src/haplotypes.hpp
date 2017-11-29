#ifndef HAPLOTYPE_FWD_ALG_H
#define HAPLOTYPE_FWD_ALG_H

#include <cmath>
#include <vector>

#include "vg.pb.h"
#include "xg.hpp"

using namespace std;
using thread_t = vector<xg::XG::ThreadMapping>;

namespace haploMath{
  double logsum(double a, double b);
  double logdiff(double a, double b);
  double int_weighted_sum(vector<double> values, vector<int64_t> counts);
  double int_weighted_sum(double* values, int64_t* counts, size_t n_entries);
}

// -----------------------------------------------------------------------------
//  RRMemo
//  Created by Jordan Eizenga on 6/21/16
// -----------------------------------------------------------------------------
//  Stores memoized values of constants and scaling factors which are used in
//  the forward matrix extension
// -----------------------------------------------------------------------------
struct RRMemo {
private:
  std::vector<double> S_multipliers;
  double T_multiplier;

  std::vector< std::vector<double> > S;
  std::vector<double> T;

  double rho;
  double exp_rho;

  double continue_probability;

  double S_value(int height, int width);
  double T_value(int width);

  std::vector<double> logS_bases;

public:
  size_t population_size;

  RRMemo(double recombination_penalty, size_t population_size);

  double recombination_penalty();
  double log_recombination_penalty();
  double cont_probability();

  double log_continue_factor(int64_t totwidth);
  double continue_factor(int64_t totwidth);

  double rr_diff(int height, int width);
  double rr_same(int height, int width);
  double rr_adj(int width);
  double rr_all(int height, int width);

  double logT_base;
  double logT(int width);
  double logS(int height, int width);
  double logRRDiff(int height, int width);
};

// -----------------------------------------------------------------------------

struct int_itvl_t{
  int64_t bottom;
  int64_t top;
  int64_t size() const;
  bool empty() const;
  static bool nondisjoint(const int_itvl_t& A, const int_itvl_t& B);
  static bool disjoint(const int_itvl_t& A, const int_itvl_t& B);
  // NB that the set of int_itvl_t's is closed under intersections but not under
  // unions or differences
  static int_itvl_t intersection(const int_itvl_t& A, const int_itvl_t& B);
};

// -----------------------------------------------------------------------------

struct haplo_DP_edge_memo {
private:  
  vector<vg::Edge> in;
  vector<vg::Edge> out;
public:
  haplo_DP_edge_memo();                      // for constructing null edge_memos
  haplo_DP_edge_memo(xg::XG& graph, 
                     xg::XG::ThreadMapping last_node, 
                     xg::XG::ThreadMapping node);
  const vector<vg::Edge>& edges_in() const;
  const vector<vg::Edge>& edges_out() const;
};

// -----------------------------------------------------------------------------

struct hDP_graph_accessor {
  const xg::XG::ThreadMapping old_node;
  const xg::XG::ThreadMapping new_node;
  const haplo_DP_edge_memo edges;
  const xg::XG& graph;
  RRMemo& memo;
  
  // accessor for noninitial nodes in a haplotype
  hDP_graph_accessor(xg::XG& graph, 
                     xg::XG::ThreadMapping old_node, 
                     xg::XG::ThreadMapping new_node, 
                     RRMemo& memo);
  // accessor for initial node in a haplotype                   
  // old_node and edge-vectors are null; do not use to extend nonempty states
  hDP_graph_accessor(xg::XG& graph, 
                     xg::XG::ThreadMapping new_node,
                     RRMemo& memo);
                     
  int64_t new_side() const;
  int64_t new_height() const;
  int64_t old_height() const;
  int64_t new_length() const;
  
  bool has_edge() const;
};

// -----------------------------------------------------------------------------

// struct hDP_gbwt_graph_accessor {
//   
// };

// -----------------------------------------------------------------------------

struct haplo_DP_rectangle{
private:
  int64_t inner_value;
  xg::XG::ThreadSearchState state;
  int64_t previous_index = -1;
public:
  double R;
  void extend(hDP_graph_accessor& ga);
  void false_extend(hDP_graph_accessor& ga, int_itvl_t delta);
  int64_t I() const;
  int64_t interval_size() const;
  void set_prev_idx(int64_t index);
  int64_t prev_idx() const;
  bool is_new() const;
  void calculate_I(int64_t succ_o_val);
};

// -----------------------------------------------------------------------------

struct haplo_DP_column {
private:
  vector<double> previous_values;
  vector<int64_t> previous_sizes;
  vector<haplo_DP_rectangle*> entries;
  double previous_sum;
  double sum;
  double length;
  void binary_extend_intervals(hDP_graph_accessor& ga, 
                               int_itvl_t indices, 
                               int_itvl_t ss_deltas, 
                               int_itvl_t state_sizes);
  void standard_extend(hDP_graph_accessor& ga);
  void update_inner_values();
  void update_score_vector(RRMemo& memo);
  double previous_R(size_t i) const;
public:
  haplo_DP_column(hDP_graph_accessor& ga);
  void extend(hDP_graph_accessor& ga);
  vector<int64_t> get_sizes() const;
  vector<double> get_scores() const;
  double current_sum() const;
};

thread_t path_to_thread_t(const vg::Path& path);

// -----------------------------------------------------------------------------

typedef pair<double, bool> haplo_score_type;

struct haplo_DP {
private:
  haplo_DP_column DP_column;
public:
  haplo_DP(hDP_graph_accessor& ga);
  haplo_DP_column* get_current_column();
  void print(ostream& out) const;
  
  static haplo_score_type score(const vg::Path& path, xg::XG& graph, RRMemo& memo);
  static haplo_score_type score(const thread_t& thread, xg::XG& graph, RRMemo& memo);
};

#endif