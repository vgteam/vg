#ifndef HAPLOTYPE_ENUMERATOR_H
#define HAPLOTYPE_ENUMERATOR_H

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <iostream>
#include <vector>

#include "vg.pb.h"
#include "xg.hpp"

using namespace std;
using thread_t = vector<xg::XG::ThreadMapping>;

double logsum(double a, double b);
double logdiff(double a, double b);

//  RRMemo functions
//  Created by Jordan Eizenga on 6/21/16.
struct RRMemo {
private:
  void initialize(double recombination_penalty);

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
  int population_size = 5008;

  RRMemo();
  RRMemo(double recombination_penalty);
  ~RRMemo(void);

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

  double logSN(vector<double> logRs, vector<int> Is);
};

class rectangle {
private:
  // We don't use these yet (we're using relative indices instead) but they will
  // be used in edit-propsal
  int a_index;
public:
  ~rectangle(void) {};
  xg::XG::ThreadSearchState state;

  // Switched to indices from pointers as quick fix for messy bug caused by
  // reindexing on removal of empty rectangles... TODO? fix this
  int prev = -1;
  int next = -1;
  int J = 0;
  int I = 0;
  double R = 0;
  double logR = 0;

  // Computes J at next_id for the strip corresponding to state
  // NB that this also calls rectangle::extend
  int get_next_J(xg::XG::ThreadMapping next_node, xg::XG& graph);
  int get_next_J(thread_t& extension, xg::XG& graph);
  int get_next_J(xg::XG::ThreadMapping next_node, xg::XG& graph, vector<vg::Edge>& edges_in, vector<vg::Edge>& edges_out);
  int get_next_J(thread_t& extension, xg::XG& graph, vector<vg::Edge>& edges_in, vector<vg::Edge>& edges_out);

  // Extends the gPBWT search state by node
  void extend(xg::XG::ThreadMapping next_node, xg::XG& graph);
  void extend(xg::XG::ThreadMapping next_node, xg::XG& graph, vector<vg::Edge>& edges_in, vector<vg::Edge>& edges_out);


  // Updates the ThreadSearchState if it can be inferred from the surrounding rectangles
  void simple_extend(thread_t& extension, xg::XG& graph, int delta_start, int delta_end);
  void simple_extend(xg::XG::ThreadMapping next_node, xg::XG& graph, int delta_start, int delta_end);
};

// A cross-section is a column of rectangles S^a_b, a <= b. Each "rectangle" in
// the sense of recomb-rectangle functions is a whole cross_section
struct cross_section {
private:
  xg::XG::ThreadMapping node;
  int b_index;
public:
  cross_section(int64_t new_height,int b,xg::XG::ThreadMapping new_node);
  ~cross_section(void) {};
  vector<rectangle> S;
  int height; // height (in consistent thread_ts)
  int width = 1; // width (in base pairs)
  inline xg::XG::ThreadMapping get_node();
  inline bool has_squashed_nodes();
  inline xg::XG::ThreadMapping get_last_node();
  // Which nodes were skipped?
  thread_t bridge;
  cross_section cs_shell();
};

// A haplo_d indexes |A| + 1 columns of rectangles S^*_b according in A-order
class haplo_d {
public:
  rectangle empty_rect;
  vector<cross_section> cs;
  int64_t tot_width = 0;
  haplo_d();
  haplo_d(const thread_t& t, xg::XG& graph);
  ~haplo_d(void) {};
  // calculate_Is() needs to be called for the cross_sections have I values in
  // their rectangles. The haplo_d constructor only builds the most recent (in
  // terms of node history) 1 or 2 rectangles at each node
  // DEPRECATED
  void calculate_Is(xg::XG& graph);
  // IN FAVOR OF
  void log_calculate_Is(xg::XG& graph);
  void seeded_log_calculate_Is(xg::XG& graph);
  void binaryI(xg::XG& graph, thread_t extension, int b, int atop, int abottom, int deltaItop, int deltaIbottom, int Jtop, int Jbottom, int level);
  void binaryI(xg::XG& graph, thread_t extension, int b, int atop, int abottom, int deltaJtop, int deltaJbottom, int Jtop, int Jbottom, int level, vector<vg::Edge>& edges_in, vector<vg::Edge>& edges_out);

  void print(ostream& stream);
  void print_detailed(ostream& stream);
  void print_detailed_searchstates(ostream& stream);
  void print_graphical(ostream& stream);
  void print_neighbours(ostream& stream);

  // accessing substructure for probability calculations:
  inline double prev_R(int b, int a);
  inline double prev_logR(int b, int a);
  inline int prev_I(int b, int a);
  vector<double> prev_logRs(int b);
  vector<int> prev_Is(int b);
  vector<double> current_logRs(int b);
  vector<int> current_Is(int b);

  double probability(RRMemo& memo);
  double log_probability(RRMemo& memo);

  pair<vector<pair<double,int>>,vector<pair<double,int>>> identify_structures(double leave_threshold,
          double return_threshold, int timeout, xg::XG& graph);

  vector<rectangle*> trace_strip(int a, int offset, int distance);
  vector<rectangle*> trace_strip(int offset);

  void build_start(xg::XG::ThreadMapping node, xg::XG& graph);

  void initialize_skeleton(thread_t& t, xg::XG& graph);
  void initialize_skeleton(thread_t& t, int start, cross_section& prevAs, xg::XG& graph);
  void initialize_skeleton(thread_t& t, pair<int,int> interval, cross_section& prevAs, xg::XG& graph);

  inline bool has_joining_node(int index);
  rectangle* joining_node(int index);
  rectangle* last_continuing(int index);
};

haplo_d recombine_arms(haplo_d& left, haplo_d& right, int left_cut, int right_join, xg::XG& graph);

// making thread_t's to operate on
thread_t path_to_thread_t(vg::Path& path);
thread_t extract_thread(xg::XG& index, xg::XG::ThreadMapping node, int64_t offset, int64_t max_length);

bool check_for_edges(int64_t old_node_id, bool old_node_is_reverse, int64_t new_node_id, bool new_node_is_reverse, xg::XG& index);
void logRR_tests(double recombination_penalty);

int find_node(thread_t& t, xg::XG::ThreadMapping node, int hint);
int find_node(thread_t& t, xg::XG::ThreadMapping node);
int find_node(haplo_d& h, xg::XG::ThreadMapping node, int hint);


#endif
