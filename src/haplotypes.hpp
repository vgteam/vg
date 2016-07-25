#ifndef HAPLOTYPE_ENUMERATOR_H
#define HAPLOTYPE_ENUMERATOR_H

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <iostream>
#include <vector>

#include "xg.hpp"

using namespace std;
using namespace xg;

//  RRMemo functions
//  Created by Jordan Eizenga on 6/21/16.
struct RRMemo {
private:

  std::vector<double> S_multipliers;
  double T_multiplier;

  std::vector< std::vector<double> > S;
  std::vector<double> T;

  double rho;
  double exp_rho;

  double S_value(int height, int width);
  double T_value(int width);

public:
  RRMemo(double recombination_penalty);
  ~RRMemo(void);

  double recombination_penalty();

  double rr_diff(int height, int width);
  double rr_same(int height, int width);
  double rr_adj(int width);
  double rr_all(int height, int width);
};

class rectangle {
private:
  XG::ThreadSearchState state;
  // We don't use these yet (we're using relative indices instead) but they will
  // be used in edit-propsal
  int a_index;
public:
  ~rectangle(void) {};
  // Pointer to the rectangle in the same strip in the previous cross-section
  rectangle* prev = nullptr;
  int J = 0;
  int I = 0;
  double R = 0;
  // Computes J at next_id for the strip corresponding to state
  // NB that this also calls rectangle::extend
  int get_next_J(XG::ThreadMapping next_node, XG& graph);
  // Extends state by node next_id
  void extend(XG::ThreadMapping next_node, XG& graph);
};

// A cross-section is a column of rectangles S^a_b, a <= b. Each "rectangle" in
// the sense of recomb-rectangle functions is a whole cross_section
struct cross_section {
private:
  XG::ThreadMapping node;
  int b_index;
public:
  cross_section(int64_t new_height,int b,XG::ThreadMapping new_node);
  ~cross_section(void) {};
  vector<rectangle> S;
  int height;
  int width = 1;
  inline XG::ThreadMapping get_node();
};

using thread_t = vector<XG::XG::ThreadMapping>;

// A haplo_d indexes |A| + 1 columns of rectangles S^*_b according in A-order
class haplo_d {
public:
  rectangle empty_rect;
  vector<cross_section> cs;
  haplo_d(const thread_t& t, XG& graph);
  ~haplo_d(void) {};
  // calculate_Is() needs to be called before the cross_sections have I values in
  // their rectangles. The haplo_d constructor only builds the most recent (in
  // terms of node history) 1 or 2 rectangles at each node
  void calculate_Is(XG& graph);
  double probability(double recombination_penalty);
};

#endif
