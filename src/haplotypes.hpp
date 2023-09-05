#ifndef HAPLOTYPE_FWD_ALG_H
#define HAPLOTYPE_FWD_ALG_H

#include <cmath>
#include <vector>
#include <iostream>
#include <unordered_set>

#include <vg/vg.pb.h>
#include "handle.hpp"

#include <gbwt/gbwt.h>
#include <gbwt/dynamic_gbwt.h>
#include <sublinearLS/reference.hpp>
#include <sublinearLS/penalty_set.hpp>
#include <sublinearLS/input_haplotype.hpp>
#include <sublinearLS/probability.hpp>


using namespace std;

////////////////////////////////////////////////////////////////////////////////
// GENERAL USE
////////////////////////////////////////////////////////////////////////////////
// 
// A. Construct the following shared objects
//    1. An index, either a
//        i.  gbwt::GBWT
//        ii. gbwt::DynamicGBWT
//    2. An appropriate ScoreProvider implementation that will use the index.
//       It is also responsible for determining the population size from its index, if able.
//       It can also implement incremental haplotype search, because we need that
//       functionality in places where the haplotype index is abstracted as a ScoreProvider.
//    3. A memo for shared values used in calculations; a
//             haplo::haploMath::RRMemo, which takes the parameters
//                    i.    double -log(recombination probability)
//                    ii.   size_t population size
//
// B. Then, on a per-query-path basis (ie a vg::Path), call
//      provider->score(const vg::Path&, haploMath::RRMemo& memo)
//    
//      this will and return an output
//        pair<double, bool> where
//             arg 1  double  log(calculate probability)
//             arg 2  bool    whether the path is valid with respect to the index
//                            in terms of whether all edges in the path exist in
//                            the index
//
//
// Internally, these will generally call a template specialization of
//    haplo::haplo_DP::score(const vg::Path&, indexType& index, haploMath::RRMemo& memo)
//    which takes in inputs
//        i.   const vg::Path& Path
//        ii.  indexType& index where indexType is one of
//             a. gbwt::GBWT
//             b. gbwt::DynamicGBWT
//        iii. haploMath::RRMemo
//
////////////////////////////////////////////////////////////////////////////////

namespace haplo {

// If this global is set, warn the user when scoring fails
extern bool warn_on_score_fail;

namespace haploMath{
  double logsum(double a, double b);
  double logdiff(double a, double b);
  double int_weighted_sum(vector<double> values, vector<int64_t> counts);
  double int_weighted_sum(double* values, int64_t* counts, size_t n_entries);

  // ---------------------------------------------------------------------------
  //  RRMemo
  //  Created by Jordan Eizenga on 6/21/16
  // ---------------------------------------------------------------------------
  //  Stores memoized values of constants and scaling factors which are used in
  //  the forward matrix extension
  // ---------------------------------------------------------------------------
  //
  //
  struct RRMemo {
  private:
    // LINEAR SPACE CONSTANTS --------------------------------------------------
    size_t population_size;
    double exp_rho;                          // lin space recombination penalty
    
    // LOG SPACE CONSTANTS -----------------------------------------------------
    double rho;                              // log space recombination penalty
    double log_continue_probability;         // 
    std::vector<double> logS_bases;

  public:
    RRMemo(double recombination_penalty, size_t population_size);

    double log_population_size();
    double log_recombination_penalty();
    double log_continue_factor(int64_t totwidth);
    
    double logT_base;
    double logT(int width);
    double logS(int height, int width);
    double logRRDiff(int height, int width);
  };
}

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

struct gbwt_thread_t {
private:
  gbwt::vector_type nodes;
  vector<size_t> node_lengths;
public:
  gbwt_thread_t();
  gbwt_thread_t(const gbwt::vector_type& nodes, const vector<size_t>& node_lengths);
  void push_back(gbwt::node_type node, size_t node_length);
  gbwt::vector_type::value_type& operator[](size_t i);
  const gbwt::vector_type::value_type& operator[](size_t i) const;
  gbwt::vector_type::value_type& back();
  const gbwt::vector_type::value_type& back() const;
  size_t nodelength(size_t i) const;
  size_t size() const;
};

gbwt_thread_t path_to_gbwt_thread_t(const vg::Path& path);

template<class GBWTType>
class hDP_gbwt_graph_accessor {
public:
  const gbwt::node_type old_node;
  const gbwt::node_type new_node;
  const GBWTType& graph;
  haploMath::RRMemo& memo;
  size_t length;
  
  // accessor for noninitial nodes in a haplotype
  hDP_gbwt_graph_accessor(GBWTType& graph, 
                          gbwt::node_type old_node, 
                          gbwt::node_type new_node,
                          size_t new_length, 
                          haploMath::RRMemo& memo);
  // accessor for initial node in a haplotype                   
  // old_node and edge-vectors are null; do not use to extend nonempty states
  hDP_gbwt_graph_accessor(GBWTType& graph, 
                          gbwt::node_type new_node,
                          size_t new_length,
                          haploMath::RRMemo& memo);
                     
  size_t new_length() const;
  size_t new_side() const;
  size_t new_height() const;
  bool inclusive_interval() const { return true; }
  
  bool has_edge() const;  
  void print(ostream& output_stream) const;
};

// -----------------------------------------------------------------------------

struct haplo_DP_rectangle{
private:
  typedef pair<size_t, size_t> gen_range_t;
  typedef size_t gen_flat_node;
  
  int64_t inner_value;
  gen_range_t state = make_pair(1, 0);
  gen_flat_node flat_node;
  int64_t previous_index = -1;
  bool int_is_inc;
public:
  haplo_DP_rectangle();
  haplo_DP_rectangle(bool inclusive_interval);
  double R;
  void set_offset(int offset);
  template<class accessorType>
  void false_extend(accessorType& ga, int_itvl_t delta);
  template<class GBWTType>
  void extend(hDP_gbwt_graph_accessor<GBWTType>& ga);
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
  vector<shared_ptr<haplo_DP_rectangle>> entries;
  double previous_sum;
  double sum;
  double length;
  template<class accessorType>
  void binary_extend_intervals(accessorType& ga, 
                               int_itvl_t indices, 
                               int_itvl_t ss_deltas, 
                               int_itvl_t state_sizes);
  template<class accessorType>
  void standard_extend(accessorType& ga);
  void update_inner_values();
  void update_score_vector(haploMath::RRMemo& memo);
  double previous_R(size_t i) const;
public:
  template<class accessorType>
  haplo_DP_column(accessorType& ga);
  ~haplo_DP_column();
  template<class accessorType>
  void extend(accessorType& ga);
  vector<int64_t> get_sizes() const;
  vector<double> get_scores() const;
  double current_sum() const;
  void print(ostream& out) const;
  bool is_empty() const;
};

//------------------------------------------------------------------------------
// Outward-facing
//------------------------------------------------------------------------------
// return type for score function:
// -  haplo_score_type.first       double      log_haplotype_likelihood
// -  haplo_score_type.second      bool        true iff all edges in path/thread
//                                             exist in the index            
typedef pair<double, bool> haplo_score_type;
//------------------------------------------------------------------------------

struct haplo_DP {
private:
  haplo_DP_column DP_column;
public:
//------------------------------------------------------------------------------
// API functions
  template<class GBWTType>
  static haplo_score_type score(const vg::Path& path, GBWTType& graph, haploMath::RRMemo& memo);
//------------------------------------------------------------------------------

// public member functions which are not part of the API
  template<class accessorType>
  haplo_DP(accessorType& ga);
  haplo_DP_column* get_current_column();
  template<class GBWTType>
  static haplo_score_type score(const gbwt_thread_t& thread, GBWTType& graph, haploMath::RRMemo& memo);
};

//------------------------------------------------------------------------------

struct linear_haplo_structure{
private:

  siteIndex* index = nullptr;
  haplotypeCohort* cohort = nullptr;
  penaltySet* penalties = nullptr;
  const vg::PathPositionHandleGraph& graph;
  vg::path_handle_t ref_path_handle;
public:
  typedef enum nodeType{
    ref_span,
    snv,
    invalid
  } nodeType;
  
  class linearUnrepresentable : public runtime_error {
    using runtime_error::runtime_error;
  };
  
  struct SNVvector{
  private:
    // Positions of SNVs
    vector<size_t> ref_positions;
    // Allele values at those positions
    vector<alleleValue> alleles;
  public:
    void push_back(alleleValue allele, size_t ref_pos, bool deletion);
    size_t ref_position(size_t i) const;
    alleleValue allele(size_t i) const;
    size_t size() const;
  };
  
  /// Make a new linear_haplo_structure with the given indexes, mutation and recombination scoring parameters, and reference path in the graph.
  /// Penalties *must* be negative, and ought to be something like -9*2.3 mutation and -6*2.3 recombination.
  linear_haplo_structure(istream& slls_index, double log_mut_penalty, double log_recomb_penalty,
                         const vg::PathPositionHandleGraph& graph, vg::path_handle_t ref_path_handle);
  ~linear_haplo_structure();
  haplo_score_type score(const vg::Path& path) const;
  
  inputHaplotype* path_to_input_haplotype(const vg::Path& path) const;
  
  // Get the position of a node in the reference path, assuming the reference
  // path is acyclic. Enforces the constraint that all the nodes in the path
  // are locally forward.
  size_t position_assuming_acyclic(int64_t node_id) const;
  nodeType get_type(int64_t node_id) const;
  // Get the node ID of the node at the given mapping index in the path.
  int64_t path_mapping_node_id(const vg::Path& path, size_t i) const;
  // Get the offset of the path in the node at the given mapping index in the path.
  size_t path_mapping_offset(const vg::Path& path, size_t i) const;
  int64_t get_SNP_ref_position(size_t node_id) const;
  SNVvector SNVs(const vg::Path& path) const;
  bool is_snv(int64_t node_id) const;
  bool is_solitary_ref(int64_t node_id) const;
  bool sn_deletion_between_ref(int64_t left, int64_t right) const;
  int64_t get_ref_following(int64_t node_id) const;
  alleleValue get_SNV_allele(int64_t node_id) const;
};


/// Incremental haplotype search range type used for ScoreProvider's
/// incremental search API. Default constructed, represents an empty or
/// un-started search. Supports an empty() and a length(). Copyable.
/// TODO: There's some overlap with the graph accessors here, but I don't
/// understand them enough to work out exactly what it is and eliminate it.
/// TODO: This should become a real type (base class and implementations
/// wrapping implementation-specific data) when we get any other
/// implementations.
using IncrementalSearchState = gbwt::SearchState;

/// Interface abstracting over the various ways of generating haplotype scores.
/// You probably want the implementations: GBWTScoreProvider, LinearScoreProvider
/// TODO: Const-ify the indexes used
class ScoreProvider {
public:
  /// Score the given path usign the given memo
  virtual pair<double, bool> score(const vg::Path&, haploMath::RRMemo& memo) = 0;
  /// Return the haplotype count (number of expected haplotypes that agree with
  /// a path that is fixed in the population) that should be used to construct
  /// the memo to pass to score, or -1 if the indexes backing the ScoreProvider
  /// do not make this information available.
  virtual int64_t get_haplotype_count() const;
  
  // We have optional support for incremental search. TODO: We need a search
  // state abstraction that encompasses GBWT search state and other search
  // state implementations.
  
  /// Return true if this ScoreProvider supports incremental search for
  /// counting haplotypes.
  virtual bool has_incremental_search() const;
  
  /// Start a new search state with the node visit described by the given
  /// Position, if incremental search is supported.
  virtual IncrementalSearchState incremental_find(const vg::Position& pos) const;
  
  /// Extend the given search state with the node visit described by the given
  /// Position, if incremental search is supported.
  virtual IncrementalSearchState incremental_extend(const IncrementalSearchState& state, const vg::Position& pos) const;
  
  virtual ~ScoreProvider() = default;
};

/// Score haplotypes using a GBWT haplotype database (normal or dynamic)
template<class GBWTType>
class GBWTScoreProvider : public ScoreProvider {
public:
  GBWTScoreProvider(GBWTType& index);
  pair<double, bool> score(const vg::Path&, haploMath::RRMemo& memo);
  
  int64_t get_haplotype_count() const;
  
  bool has_incremental_search() const;
  IncrementalSearchState incremental_find(const vg::Position& pos) const;
  IncrementalSearchState incremental_extend(const IncrementalSearchState& state, const vg::Position& pos) const;
private:
  GBWTType& index;
};

/// Score haplotypes using a linear_haplo_structure
class LinearScoreProvider : public ScoreProvider {
public:
  LinearScoreProvider(const linear_haplo_structure& index);
  pair<double, bool> score(const vg::Path&, haploMath::RRMemo& memo);
private:
  const linear_haplo_structure& index;
};


//------------------------------------------------------------------------------
// template implementations
//------------------------------------------------------------------------------

template<class GBWTType>
hDP_gbwt_graph_accessor<GBWTType>::hDP_gbwt_graph_accessor(GBWTType& graph, 
                                                 gbwt::node_type new_node,
                                                 size_t new_length, 
                                                 haploMath::RRMemo& memo) :
  graph(graph), length(new_length),
  old_node(gbwt::invalid_node()), new_node(new_node), memo(memo) {
  
  // Nothing will work well if we have a node that has no record at all in the GBWT.
  assert(graph.contains(new_node));
    
}

template<class GBWTType>
hDP_gbwt_graph_accessor<GBWTType>::hDP_gbwt_graph_accessor(GBWTType& graph, 
                                                 gbwt::node_type old_node, 
                                                 gbwt::node_type new_node,
                                                 size_t new_length, 
                                                 haploMath::RRMemo& memo) :
  graph(graph), length(new_length),
  old_node(old_node), new_node(new_node), memo(memo) {
  
  // Nothing will work well if we have a node that has no record at all in the GBWT.
  assert(graph.contains(new_node));
  if (old_node != gbwt::invalid_node()) {
    assert(graph.contains(old_node));
  }
    
}

template<class GBWTType>
size_t hDP_gbwt_graph_accessor<GBWTType>::new_length() const {
  return length;
}

template<class GBWTType>
size_t hDP_gbwt_graph_accessor<GBWTType>::new_side() const {
  return new_node;
}

template<class GBWTType>
size_t hDP_gbwt_graph_accessor<GBWTType>::new_height() const {
  return graph.nodeSize(new_node);
}

template<class GBWTType>
void hDP_gbwt_graph_accessor<GBWTType>::print(ostream& output_stream) const {
  output_stream << "From node: ID " << gbwt::Node::id(old_node) << " is_reverse " << gbwt::Node::is_reverse(old_node) << " ; To node: ID " << gbwt::Node::id(new_node) << " is_reverse " << gbwt::Node::is_reverse(new_node);
  output_stream << " ; Reference haplotypes visiting To Node: " << new_height() << endl;
}

template<class GBWTType>
bool hDP_gbwt_graph_accessor<GBWTType>::has_edge() const {
  return graph.hasEdge(old_node, new_node);
}

//------------------------------------------------------------------------------

template<class GBWTType>
void haplo_DP_rectangle::extend(hDP_gbwt_graph_accessor<GBWTType>& ga) {
  if(previous_index == -1) {
    state = make_pair(0, ga.graph.nodeSize(ga.new_node) - 1);
  } else {
    state = ga.graph.LF(ga.old_node, state, ga.new_node);
  }
  flat_node = ga.new_node;
  inner_value = -1;
}

template<class accessorType>
void haplo_DP_rectangle::false_extend(accessorType& ga, 
                                      int_itvl_t delta) {
  flat_node = ga.new_side();
  state.first -= delta.bottom;
  state.second -= delta.top;
}

//------------------------------------------------------------------------------

template<class accessorType>
haplo_DP_column::haplo_DP_column(accessorType& ga) {
  haplo_DP_rectangle* first_rectangle = new haplo_DP_rectangle(ga.inclusive_interval());
  assert(first_rectangle != nullptr);
  entries.push_back(shared_ptr<haplo_DP_rectangle>(first_rectangle));
  first_rectangle->extend(ga);
  update_inner_values();
  update_score_vector(ga.memo);
}

template<class accessorType>
void haplo_DP_column::standard_extend(accessorType& ga) {
  previous_values = get_scores();
  previous_sizes = get_sizes();
  haplo_DP_rectangle* new_rectangle = new haplo_DP_rectangle(ga.inclusive_interval());
  assert(new_rectangle != nullptr);
  new_rectangle->extend(ga);
  decltype(entries) new_entries;
  new_entries.push_back(shared_ptr<haplo_DP_rectangle>(new_rectangle));
  for(size_t i = 0; i < entries.size(); i++) {
    // extend candidate
    auto candidate = entries.at(i);
    candidate->set_prev_idx(i);
    candidate->extend(ga);
    // check if the last rectangle added is nonempty
    if(candidate->interval_size() == new_entries.back()->interval_size()) {
      new_entries.pop_back();
    }
    if(candidate->interval_size() != 0) {
      new_entries.push_back(candidate);
    } else {
      break;
    }
  }
  entries = new_entries;
  
  if(is_empty()) {
    cerr << "In haplotype scoring at ";
    ga.print(cerr);
    throw runtime_error("Haplotype scoring calculation produced empty rectangle-vector");
  }
  
  update_inner_values();
}

template<class accessorType>
void haplo_DP_column::extend(accessorType& ga) {
  standard_extend(ga);
  length = (double)(ga.new_length());
  update_score_vector(ga.memo);
}

//------------------------------------------------------------------------------

template<class accessorType>
haplo_DP::haplo_DP(accessorType& ga) : DP_column(ga) {
  
}

template<class GBWTType>
haplo_score_type haplo_DP::score(const vg::Path& path, GBWTType& graph, haploMath::RRMemo& memo) {
  return score(path_to_gbwt_thread_t(path), graph, memo);
}

template<class GBWTType>
haplo_score_type haplo_DP::score(const gbwt_thread_t& thread, GBWTType& graph, haploMath::RRMemo& memo) {
  if (thread.size() == 0) {
    if (warn_on_score_fail) {
      cerr << "[WARNING] Path is empty and cannot be scored" << endl;
      cerr << "Cannot compute a meaningful haplotype likelihood score" << endl;
    }
    return pair<double, bool>(nan(""), false);
  }
  if (!graph.contains(thread[0])) {
    // We start on a node that has no haplotype index entry
    if (warn_on_score_fail) {
      cerr << "[WARNING] Path starts outside of haplotype index and cannot be scored" << endl;
      cerr << "Cannot compute a meaningful haplotype likelihood score" << endl;
    }
    return pair<double, bool>(nan(""), false);
  }
  
  hDP_gbwt_graph_accessor<GBWTType> ga_i(graph, thread[0], thread.nodelength(0), memo);
  haplo_DP hdp(ga_i);
  if(ga_i.new_height() == 0) {
    if (warn_on_score_fail) {
      cerr << "[WARNING] Initial node in path is visited by 0 reference haplotypes" << endl;
      cerr << "Cannot compute a meaningful haplotype likelihood score" << endl;
      ga_i.print(cerr);
    }
    return pair<double, bool>(nan(""), false);
  }
#ifdef debug
  cerr << "After entry 0 (" << gbwt::Node::id(thread[0]) << ") height: " << ga_i.new_height() << " intervals: ";
  for (auto& interval : hdp.DP_column.get_sizes()) {
    cerr << interval << " ";
  }
  cerr << "score: " << hdp.DP_column.current_sum() << endl;
#endif
  for(size_t i = 1; i < thread.size(); i++) {
    if (!graph.contains(thread[i])) {
      if (warn_on_score_fail) { 
        cerr << "[WARNING] Node " << i + 1 << " in path leaves haplotype index and cannot be scored" << endl;
        cerr << "Cannot compute a meaningful haplotype likelihood score" << endl;
      }
      return pair<double, bool>(nan(""), false);
    }
    
    hDP_gbwt_graph_accessor<GBWTType> ga(graph, thread[i-1], thread[i], thread.nodelength(i), memo);
    if(ga.new_height() == 0) {
      if (warn_on_score_fail) {
        cerr << "[WARNING] Node " << i + 1 << " in path is visited by 0 reference haplotypes" << endl;
        cerr << "Cannot compute a meaningful haplotype likelihood score" << endl;
        ga.print(cerr);
      }
      return pair<double, bool>(nan(""), false);
    } else {
      hdp.DP_column.extend(ga);
    }
#ifdef debug
    cerr << "After entry " << i << "/" << thread.size() << " (" << gbwt::Node::id(thread[i]) << ") height: " << ga.new_height() << " intervals: ";
    for (auto& interval : hdp.DP_column.get_sizes()) {
      cerr << interval << " ";
    }
    cerr << "score: " << hdp.DP_column.current_sum() << endl;
#endif
  }
  return pair<double, bool>(hdp.DP_column.current_sum(), true);
}

//------------------------------------------------------------------------------

template<class GBWTType>
GBWTScoreProvider<GBWTType>::GBWTScoreProvider(GBWTType& index) : index(index) {
  // Nothing to do!
}

template<class GBWTType>
pair<double, bool> GBWTScoreProvider<GBWTType>::score(const vg::Path& path, haploMath::RRMemo& memo) {
  return haplo_DP::score(path, index, memo);
}

template<class GBWTType>
int64_t GBWTScoreProvider<GBWTType>::get_haplotype_count() const {
  if (!index.hasMetadata()) {
    // No metadata available
    return -1;
  }
  
  // TODO: Does this haplotype count have the same expected-count-for-fixed-path semantics that we want?
  // Or does it count fragments of haplotypes?
  return index.metadata.haplotypes();
}

template<class GBWTType>
bool GBWTScoreProvider<GBWTType>::has_incremental_search() const {
  // We are going to implement incremental search.
  return true;
}

template<class GBWTType>
IncrementalSearchState GBWTScoreProvider<GBWTType>::incremental_find(const vg::Position& pos) const {
  return index.find(gbwt::Node::encode(pos.node_id(), pos.is_reverse()));
}

template<class GBWTType>
IncrementalSearchState GBWTScoreProvider<GBWTType>::incremental_extend(const IncrementalSearchState& state, const vg::Position& pos) const {
  return index.extend(state, gbwt::Node::encode(pos.node_id(), pos.is_reverse()));
}

} // namespace haplo

#endif
