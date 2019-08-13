#include "haplotypes.hpp"
#include "path.hpp"
#include "position.hpp"

namespace haplo {

// By default, should we warn when haplotype scoring fails?
bool warn_on_score_fail = false;

/*******************************************************************************
hDP_gbwt_graph_accessor
*******************************************************************************/

gbwt_thread_t::gbwt_thread_t() {
}

gbwt_thread_t::gbwt_thread_t(const gbwt::vector_type& nodes, const vector<size_t>& node_lengths) : nodes(nodes), node_lengths(node_lengths) {
}

void gbwt_thread_t::push_back(gbwt::node_type node, size_t node_length) {
  nodes.push_back(node);
  node_lengths.push_back(node_length);
}

gbwt::vector_type::value_type& gbwt_thread_t::operator[](size_t i) {
  return nodes[i];  
}

const gbwt::vector_type::value_type& gbwt_thread_t::operator[](size_t i) const {
  return nodes[i];  
}

gbwt::vector_type::value_type& gbwt_thread_t::back() {
  return nodes.back();
}

const gbwt::vector_type::value_type& gbwt_thread_t::back() const {
  return nodes.back();
}

size_t gbwt_thread_t::nodelength(size_t i) const {
  return node_lengths[i];  
}

size_t gbwt_thread_t::size() const {
  return nodes.size();
}

gbwt_thread_t path_to_gbwt_thread_t(const vg::Path& path) {
  gbwt_thread_t t;
  for(size_t i = 0; i < path.mapping_size(); i++) {
    vg::Mapping mapping = path.mapping(i);
    auto pos = mapping.position();
    gbwt::node_type n = gbwt::Node::encode(pos.node_id(), pos.is_reverse());
    size_t node_length = vg::mapping_from_length(mapping);
    t.push_back(n, node_length);
  }
  return t;
}

/*******************************************************************************
haplo_DP_rectangle
*******************************************************************************/

haplo_DP_rectangle::haplo_DP_rectangle() {
  
}

haplo_DP_rectangle::haplo_DP_rectangle(bool inclusive_interval) : int_is_inc(inclusive_interval) {
  
}

void haplo_DP_rectangle::calculate_I(int64_t succ_o_val) {
  inner_value = interval_size() - succ_o_val;
}

int64_t haplo_DP_rectangle::I() const {
  return inner_value;
}

int64_t haplo_DP_rectangle::interval_size() const {
  return state.second - state.first + int_is_inc;
}

void haplo_DP_rectangle::set_prev_idx(int64_t index) {
  previous_index = index;
}

int64_t haplo_DP_rectangle::prev_idx() const {
  return previous_index;
}

bool haplo_DP_rectangle::is_new() const {
  return (previous_index == -1);
}

/*******************************************************************************
int_itvl_t
*******************************************************************************/

int64_t int_itvl_t::size() const {
  return top - bottom;
}

bool int_itvl_t::empty() const {
  return top == bottom;
}

bool int_itvl_t::nondisjoint(const int_itvl_t& A, const int_itvl_t& B) {
  return (B.bottom <= A.top) && (A.bottom <= B.top);
}

bool int_itvl_t::disjoint(const int_itvl_t& A, const int_itvl_t& B) {
  return !nondisjoint(A, B);
}

int_itvl_t int_itvl_t::intersection(const int_itvl_t& A, const int_itvl_t& B) {
  int_itvl_t I;
  if(int_itvl_t::nondisjoint(A, B)) {
    I.top = A.top > B.top ? B.top : A.top;
    I.bottom = A.bottom > B.bottom ? B.bottom : A.bottom;
  } else {
    I.top = I.bottom = 0;
  }
  return I;
}

/*******************************************************************************
haplo_DP_column
*******************************************************************************/

haplo_DP_column::~haplo_DP_column() {
}

void haplo_DP_column::update_inner_values() {
  for(size_t i = 0; i + 1 < entries.size(); i++) {
    assert(entries[i].get() != nullptr);
    assert(entries[i+1].get() != nullptr);
    entries[i]->calculate_I(entries[i+1]->interval_size());
  }
  if(!entries.empty()) {
    (entries.back())->calculate_I(0);
  }
}

// void haplo_DP_column::binary_extend_intervals(hDP_graph_accessor& ga, int_itvl_t indices, int_itvl_t ss_deltas, int_itvl_t state_sizes) {
//   if(indices.size() <= 1) {
//     return;
//   } else if(ss_deltas.size() == 0) {
//     int_itvl_t correction = prev - curr range start
//     for(size_t i = indices.bottom + 1; i < indices.top; i++) {
//       
//       rect.false_extend(ga, correction);
//       
//     }
//   } else {
//     int64_t mid_index = indices.bottom + indices.size()/2;
//     
//   }
// }

void haplo_DP_column::update_score_vector(haploMath::RRMemo& memo) {
  assert(!entries.empty());
  auto r_0 = entries.at(0);
  assert(r_0.get() != nullptr);
  if(entries.size() == 1 && entries.at(0)->prev_idx() == -1) {
    r_0->R = -memo.log_population_size();
    sum = r_0->R + log(r_0->interval_size());
    previous_values = {r_0->R};
    previous_sizes = {r_0->interval_size()};
    return;
  } else {
    previous_sum = sum;
    int64_t offset = (int64_t)(entries.at(0)->is_new());
    // TODO optimize by using arrays
    vector<double> continuing_Rs(entries.size() - offset);
    vector<int64_t> continuing_counts(entries.size() - offset);
    for(size_t i = offset; i < entries.size(); i++) {
      assert(entries[i].get() != nullptr);
      continuing_Rs.at(i - offset) = previous_values[entries[i]->prev_idx()];
      continuing_counts.at(i - offset) = entries[i]->I();
    }
    
    double logpS1S2RRS = previous_sum + 
                         memo.log_recombination_penalty() + 
                         memo.logS(r_0->interval_size(), length);
    double logS1 = haploMath::int_weighted_sum(continuing_Rs, continuing_counts);
    double logS1RRD = logS1 + memo.logRRDiff(r_0->interval_size(), length);
    size_t i = 0;
    if(r_0->prev_idx() == -1) {
      r_0->R = logpS1S2RRS;
      i = 1;
    }
    if(length == 1) {
      for(; i < entries.size(); i++) {
        assert(entries[i].get() != nullptr);
        double logLHS = memo.logT_base +
                        previous_R(i) +
                        memo.logT(length);
        entries[i]->R = haploMath::logsum(logLHS, logpS1S2RRS);
      }
    } else {
      for(; i < entries.size(); i++) {
        assert(entries[i].get() != nullptr);
        double logLHS = memo.logT_base +
                        haploMath::logsum(logS1RRD, previous_R(i) + memo.logT(length));
        entries[i]->R = haploMath::logsum(logLHS, logpS1S2RRS);
      }
    }
  }
  previous_values = get_scores();
  previous_sizes = get_sizes();
  sum = haploMath::int_weighted_sum(previous_values, previous_sizes);
}

double haplo_DP_column::previous_R(size_t i) const {
  assert(entries.at(i).get() != nullptr);
  return previous_values[(entries.at(i))->prev_idx()];
}

vector<double> haplo_DP_column::get_scores() const {
  vector<double> to_return;
  for(size_t i = 0; i < entries.size(); i++) {
    assert(entries[i].get() != nullptr);
    to_return.push_back(entries[i]->R);
  }
  return to_return;
}

vector<int64_t> haplo_DP_column::get_sizes() const {
  vector<int64_t> to_return;
  for(size_t i = 0; i < entries.size(); i++) {
    assert(entries[i].get() != nullptr);
    to_return.push_back(entries[i]->I());
  }
  return to_return;
}

double haplo_DP_column::current_sum() const {
  return sum;
}

void haplo_DP_column::print(ostream& out) const {
  for(size_t i = 0; i < get_sizes().size(); i++) {
    out << "[";
    for(size_t j = 0; j < get_sizes().size() - i - 1; j++) {
      out << "  ";
    }
    assert(entries.at(i).get() != nullptr);
    out << entries.at(i)->I() << "] : " << entries.at(i)->interval_size() << endl;
  }
}

bool haplo_DP_column::is_empty() const {
  return entries.size() == 0;
}

haplo_DP_column* haplo_DP::get_current_column() {
  return &DP_column;
}

/*******************************************************************************
linear_haplo_structure
*******************************************************************************/

linear_haplo_structure::nodeType linear_haplo_structure::get_type(int64_t node_id) const {
  if(is_solitary_ref(node_id)) {
    return ref_span;
  }
  if(is_snv(node_id)) {
    return snv;
  }
  return invalid;
}

int64_t linear_haplo_structure::path_mapping_node_id(const vg::Path& path, size_t i) const {
  vg::Mapping this_mapping = path.mapping(i);
  auto last_pos = this_mapping.position();
  return last_pos.node_id();
}

size_t linear_haplo_structure::path_mapping_offset(const vg::Path& path, size_t i) const {
  vg::Mapping this_mapping = path.mapping(i);
  auto last_pos = this_mapping.position();
  return last_pos.offset();
}

int64_t linear_haplo_structure::get_SNP_ref_position(size_t node_id) const {
    // walk to the left and the neighbor there
    vg::handle_t lnbr;
    bool found_lnbr = !graph.follow_edges(graph.get_handle(node_id), true,
                                          [&](const vg::handle_t& prev) {
        lnbr = prev;
        return false;
                                          });
    if (!found_lnbr) {
        throw runtime_error("SNP at node ID " + to_string(node_id) + " does not have neighbors that can be used to find reference path " + graph.get_path_name(ref_path_handle));
    }
    
    // walk back to the right and get the position of the allele that's on the
    // reference path
    int64_t ref_pos;
    bool found_ref_pos = !graph.follow_edges(lnbr, false,
                                             [&](const vg::handle_t& next) {
        return graph.for_each_step_on_handle(next, [&](const vg::step_handle_t& step) {
            if (graph.get_path_handle_of_step(step) == ref_path_handle) {
                ref_pos = graph.get_position_of_step(step);
                return false;
            }
            return true;
        });
    });
    
    if (!found_ref_pos) {
        throw runtime_error("SNP at node ID " + to_string(node_id) + " is not adjacent to the reference path " + graph.get_path_name(ref_path_handle));
    }
    
    return ref_pos;
}

void linear_haplo_structure::SNVvector::push_back(alleleValue allele, size_t ref_pos, bool deletion) {
  ref_positions.push_back(ref_pos);
  if(deletion) {
    alleles.push_back(gap);
  } else {
    alleles.push_back(allele);
  }
}

alleleValue linear_haplo_structure::get_SNV_allele(int64_t node_id) const {
  char allele_char = graph.get_sequence(graph.get_handle(node_id)).at(0);
  return allele::from_char(allele_char);
}

size_t linear_haplo_structure::SNVvector::ref_position(size_t i) const {
  return ref_positions.at(i);
}

alleleValue linear_haplo_structure::SNVvector::allele(size_t i) const {
  return alleles.at(i);
}

size_t linear_haplo_structure::SNVvector::size() const {
  return alleles.size();
}

bool linear_haplo_structure::sn_deletion_between_ref(int64_t left, int64_t right) const {
  int64_t gap = position_assuming_acyclic(right) - position_assuming_acyclic(left) - graph.get_length(graph.get_handle(left));
  if(gap == 0) {
    return false;
  } else if(gap == 1) {
    return true;
  } else {
    throw linearUnrepresentable("indel not of length 1");
  }
}

int64_t linear_haplo_structure::get_ref_following(int64_t node_id) const {
    // walk to the right and get all nodes on the reference path
    vector<vg::handle_t> refs;
    graph.follow_edges(graph.get_handle(node_id), false, [&](const vg::handle_t& next) {
        graph.for_each_step_on_handle(next, [&](const vg::step_handle_t& step) {
            if (graph.get_path_handle_of_step(step) == ref_path_handle) {
                refs.push_back(graph.get_handle_of_step(step));
                return false;
            }
            return true;
        });
    });
        
    if (refs.empty()) {
        throw runtime_error("SNP at node ID " + to_string(node_id) + " does not have a following node on the reference path " + graph.get_path_name(ref_path_handle));
    }
    
    size_t smallest = numeric_limits<size_t>::max();
    vg::handle_t node_at_smallest;
    for(size_t i = 0; i < refs.size(); i++) {
        auto pos = position_assuming_acyclic(graph.get_id(refs[i]));
        if(pos < smallest) {
            smallest = pos;
            node_at_smallest = refs[i];
        }
    }
    return graph.get_id(node_at_smallest);
}

linear_haplo_structure::SNVvector linear_haplo_structure::SNVs(const vg::Path& path) const {
  SNVvector to_return;
  if(path.mapping_size() < 1) {
    return to_return;
  }
  
  int64_t last_node = path_mapping_node_id(path, 0);
  nodeType last_type = get_type(last_node);
  size_t last_pos;
  if(last_type == ref_span) {
    // Find where the last_node (and not the path) starts in the reference
    last_pos = position_assuming_acyclic(last_node);
  } else if(last_type == snv) {
    last_pos = get_SNP_ref_position(last_node);
    to_return.push_back(get_SNV_allele(last_node), last_pos, false);
  } else {
    throw linearUnrepresentable("not an SNV");
  }

  for(size_t i = 1; i < path.mapping_size(); i++) {
    int64_t this_node = path_mapping_node_id(path, i);
    // Only the first mapping can have an offset.
    assert(path_mapping_offset(path, i) == 0);
    nodeType this_type = get_type(this_node);
    size_t this_pos;

    if(this_type == invalid) {
      throw linearUnrepresentable("not an SNV");
    } else if(this_type == snv) {
      this_pos = get_SNP_ref_position(this_node);
      if(this_pos != last_pos + graph.get_length(graph.get_handle(last_node))) {
        throw linearUnrepresentable("indel immediately before SNV");
      }
      to_return.push_back(get_SNV_allele(this_node), this_pos, false);
    } else {
      this_pos = position_assuming_acyclic(this_node);
      if(last_type == snv) {
        if(this_pos != last_pos + graph.get_length(graph.get_handle(last_node))) {
          throw linearUnrepresentable("indel immediately after SNV");
        }
      } else {
        if(sn_deletion_between_ref(last_node, this_node)) {
          int64_t ref_node = get_ref_following(last_node);
          to_return.push_back(gap, position_assuming_acyclic(ref_node), true);
        }
      }
    }
    
    last_node = this_node;
    last_type = this_type;
    last_pos = this_pos;
  }
  return to_return;
}

size_t linear_haplo_structure::position_assuming_acyclic(int64_t node_id) const {
    
    // check occurrences o this node on paths
    size_t pos;
    bool found_pos = !graph.for_each_step_on_handle(graph.get_handle(node_id), [&](const vg::step_handle_t& step) {
        
        // get the pos if the path matches
        if (graph.get_path_handle_of_step(step) == ref_path_handle) {
            pos = graph.get_position_of_step(step);
            return false;
        }
        return true;
    });
    
  if (!found_pos) {
      throw runtime_error("requested position-in-path of node " + to_string(node_id) + " not in path " + graph.get_path_name(ref_path_handle));
  }
    
  return pos;
}


bool linear_haplo_structure::is_solitary_ref(int64_t node_id) const {
    vg::handle_t handle = graph.get_handle(node_id);
    bool on_ref = !graph.for_each_step_on_handle(handle, [&](const vg::step_handle_t& step) {
        return graph.get_path_handle_of_step(step) != ref_path_handle;
    });
    
    if (!on_ref) {
        return false;
    }
    
    bool is_deletion_neighbour = true;
    graph.follow_edges(handle, true, [&](const vg::handle_t& prev) {
        if (graph.get_degree(prev, false) != 1) {
            graph.follow_edges(prev, false, [&](const vg::handle_t& next) {
                if (next != handle) {
                    size_t rr_count = 0;
                    graph.follow_edges(next, false, [&](const vg::handle_t& next_next) {
                        rr_count++;
                        if (next_next != handle || rr_count > 1) {
                            is_deletion_neighbour = false;
                        }
                    });
                }
            });
        }
    });
    
    graph.follow_edges(handle, false, [&](const vg::handle_t& next) {
        if (graph.get_degree(next, true) != 1) {
            graph.follow_edges(next, true, [&](const vg::handle_t& prev) {
                if (prev != handle) {
                    size_t ll_count = 0;
                    graph.follow_edges(prev, true, [&](const vg::handle_t& prev_prev) {
                        ll_count++;
                        if (prev_prev != handle || ll_count > 1) {
                            is_deletion_neighbour = false;
                        }
                    });
                }
            });
        }
    });
    
    if (!is_deletion_neighbour) {
        return false;
    }
    return true;
}

bool linear_haplo_structure::is_snv(int64_t node_id) const {
    // has only one left and one right neighbour
    vg::handle_t lnbr, rnbr;
    size_t lnbr_count = 0, rnbr_count = 0;
    
    vg::handle_t handle = graph.get_handle(node_id);
    graph.follow_edges(handle, true, [&](const vg::handle_t& prev) {
        lnbr = prev;
        ++lnbr_count;
    });
    graph.follow_edges(handle, false, [&](const vg::handle_t& next) {
        rnbr = next;
        ++rnbr_count;
    });
    
    if (lnbr_count != 1 || rnbr_count != 1) {
        // has too many or too few neighbours to be an SNV
        return false;
    }
    
    unordered_set<vg::handle_t> from_lnbr;
    bool all_snv = graph.follow_edges(lnbr, false, [&](const vg::handle_t& next) {
        from_lnbr.insert(next);
        return graph.get_length(next) == 1;
    });
    
    if (!all_snv) {
        // some alleles are not SNVs
        return false;
    }
    
    size_t from_rnbr_count = 0;
    bool all_match = graph.follow_edges(rnbr, true, [&](const vg::handle_t& prev) {
        ++from_rnbr_count;
        return (bool) from_lnbr.count(prev);
    });
    
    if (!all_match || from_rnbr_count != from_lnbr.size()) {
        // we didn't find all of the neighbors in common
        return false;
    }
    
    return true;
}

inputHaplotype* linear_haplo_structure::path_to_input_haplotype(const vg::Path& path) const {
  if(path.mapping_size() == 0) {
    return new inputHaplotype();
  }
  
  SNVvector SNV_candidates;
  
  try {
    SNV_candidates = SNVs(path);
  } catch(linearUnrepresentable& e) {
    return new inputHaplotype();
  }
  
  // Determine the start position of the path in the reference.
  size_t start;
  // This is complicated by the fact that the path may start on an SNV alt
  // allele node which does not occur on the reference path.
  int64_t start_node = path_mapping_node_id(path, 0);
  nodeType start_type = get_type(start_node);
  if(start_type == ref_span) {
    // We start on a reference node.
    start = position_assuming_acyclic(start_node) + path_mapping_offset(path, 0);
  } else if(start_type == snv) {
    // We start on an alt
    start = get_SNP_ref_position(start_node) + path_mapping_offset(path, 0);
  } else {
    // We start at some other weird place.
    return new inputHaplotype();
  }
 
  size_t length = 0;
  for(size_t i = 0; i < path.mapping_size(); i++) {
    vg::Mapping mapping = path.mapping(i);
    length += vg::mapping_from_length(mapping);
  }  
  
  if(SNV_candidates.size() == 0) {
    // We need no sites, but one before-first site, after-last-site novel SNP
    // count entry. TODO: actually count novel SNPs?
    return new inputHaplotype(vector<alleleValue>(0), vector<size_t>(1, 0), index, start, length);
  }
  
  vector<size_t> positions;
  vector<alleleValue> alleles;
  
  for(size_t i = 0; i < SNV_candidates.size(); i++) {
    if(index->is_site(SNV_candidates.ref_position(i))) {
      positions.push_back(SNV_candidates.ref_position(i));
      alleles.push_back(SNV_candidates.allele(i));
    }
  }

  for(size_t i = 1; i < positions.size(); i++) {
    if(positions[i] <= positions[i - 1]) {
      return new inputHaplotype();
    }
  }
  
  size_t last_i = std::numeric_limits<size_t>::max();
  for(size_t i = 0; i < SNV_candidates.size(); i++) {
    if(!index->is_site(SNV_candidates.ref_position(i))) {
      // This is a single-base indel or something else that looks like a SNP but isn't.
      // We can skip it like we do with all non-SNP variation.
      continue;
    }
    
    // If we know it's int he index, get its index in the index
    size_t this_i = index->get_site_index(SNV_candidates.ref_position(i));
    if(last_i != std::numeric_limits<size_t>::max() && this_i != last_i + 1) {
      // The last SNP existed and this SNP does not come after it.
      // All the SNPs have to be consecutive in the index.
      return new inputHaplotype();
    } else {
      // Either this was the first SNP or it followed right after the previous one in the sites list.
      last_i = this_i;
    }
  }
  
  // With actual alleles, we need one novel SNP count between each pair of
  // alleles, plus one before the first and one after the last. 
  inputHaplotype* to_return = new inputHaplotype(alleles, vector<size_t>(alleles.size() + 1, 0), index, start, length);
  return to_return;
}

linear_haplo_structure::linear_haplo_structure(istream& slls_index, double log_mut_penalty,
                                               double log_recomb_penalty,
                                               const vg::PathPositionHandleGraph& graph,
                                               vg::path_handle_t ref_path_handle) : graph(graph), ref_path_handle(ref_path_handle) {
  
  if (log_mut_penalty > 0) {
    throw runtime_error("log mutation penalty must be negative");
  }
  
  if (log_recomb_penalty > 0) {
    throw runtime_error("log recombination penalty must be negative");
  }
  
  index = new siteIndex(slls_index);
  cohort = new haplotypeCohort(slls_index, index);
  penalties = new penaltySet(log_recomb_penalty, log_mut_penalty, cohort->get_n_haplotypes());
}

linear_haplo_structure::~linear_haplo_structure() {
  delete index;
  delete cohort;
  delete penalties;
}

haplo_score_type linear_haplo_structure::score(const vg::Path& path) const {
  inputHaplotype* observed = path_to_input_haplotype(path);
  haplo_score_type to_return;
  if(observed->is_valid()) {
    fastFwdAlgState observed_state(index, penalties, cohort);
    double result = observed_state.calculate_probability(observed);
    to_return = haplo_score_type(result, true);
  } else {
    to_return = haplo_score_type(nan(""), false);
  }  
  delete observed;
  return(to_return);
}

/*******************************************************************************
ScoreProvider
*******************************************************************************/

int64_t ScoreProvider::get_haplotype_count() const {
  // By default, say that we don't know the haplotype count.
  return -1;
}

bool ScoreProvider::has_incremental_search() const {
  // By default, say that we lack incremental search support.
  return false;
}

IncrementalSearchState ScoreProvider::incremental_find(const vg::Position& pos) const {
  throw runtime_error("Incremental search not implemented");
}

IncrementalSearchState ScoreProvider::incremental_extend(const IncrementalSearchState& state, const vg::Position& pos) const {
  throw runtime_error("Incremental search not implemented");
}

/*******************************************************************************
LinearScoreProvider
*******************************************************************************/

LinearScoreProvider::LinearScoreProvider(const linear_haplo_structure& index) : index(index) {
  // Nothing to do!
}

pair<double, bool> LinearScoreProvider::score(const vg::Path& path, haploMath::RRMemo& memo) {
  // Memo is ignored; all penalties come from the index itself.
  auto scored = index.score(path);
  
  if (scored.second) {
      // Yohei says there should never be NANs when it worked
      assert(!std::isnan(scored.first));
  }
  
  return scored;
}

/*******************************************************************************
math functions
*******************************************************************************/

namespace haploMath {

RRMemo::RRMemo(double recombination_penalty, size_t population_size) : 
    population_size(population_size) {
  
  rho = -recombination_penalty - log(population_size - 1);
  exp_rho = exp(rho);
  assert(exp_rho < 1);
  
  log_continue_probability = log1p(- exp_rho * (population_size - 1));

  // log versions
  logT_base = log1p(-exp_rho);
  // Populate the tabel out to twice the haplotype count.
  // In regions between unphased variants, we can have twice as many hits as real haplotypes in the index.
  for(int i = 0; i < population_size * 2; i++) {
    logS_bases.push_back(log1p(i*exp_rho));
  }
}

double logdiff(double a, double b) {
  if(b > a) {
    double c = a;
    a = b;
    b = c;
  }
  return a + log1p(-exp(b - a));
}

double logsum(double a, double b) {
  if(b > a) {
    double c = a;
    a = b;
    b = c;
  }
  return a + log1p(exp(b - a));
}

double int_weighted_sum(double* values, int64_t* counts, size_t n_values) {
  if(n_values == 0) {
    return 0;
  } else if(n_values == 1) {
    return values[0] + log(counts[0]);
  } else {
    double max_summand = values[0] + log(counts[0]);
    int max_index = 0;
    vector<double> summands;
    for(int i = 0; i < n_values; i++){
      summands.push_back(values[i] + log(counts[i]));
      if(summands.back() > max_summand) {
        max_summand = summands.back();
        max_index = i;
      }
    }
    double sum = 0;
    for(int i = 0; i < summands.size(); i++) {
      if(i != max_index) {
        sum += exp(summands[i]-max_summand);
      }
    }
    return max_summand + log1p(sum);
  }
}


double int_weighted_sum(vector<double> values, vector<int64_t> counts) {
  if(values.size() == 0) {
    return 0;
  } else if(values.size() == 1) {
    return values[0] + log(counts[0]);
  } else {
    double max_summand = values[0] + log(counts[0]);
    int max_index = 0;
    vector<double> summands;
    for(int i = 0; i < values.size(); i++){
      summands.push_back(values[i] + log(counts[i]));
      if(summands.back() > max_summand) {
        max_summand = summands.back();
        max_index = i;
      }
    }
    double sum = 0;
    for(int i = 0; i < summands.size(); i++) {
      if(i != max_index) {
        sum += exp(summands[i]-max_summand);
      }
    }
    return max_summand + log1p(sum);
  }
}

double RRMemo::logT(int width) {
  return (width-1)*logT_base; //logT_base = log(1 - exp_rho)
}

double RRMemo::logS(int height, int width) {
  if (height <= logS_bases.size()) {
    // Fulfil from lookup table
    return (width-1)*logS_bases[height-1]; //logS_base = log(1 + i*exp_rho)
  } else {
    // We must have a cycle or something; we have *way* more hits than haplotypes.
    // Uncommon; just recompute the logS base as we do in the constructor.
    return (width-1)*log1p((height-1)*exp_rho);
  }
}

double RRMemo::logRRDiff(int height, int width) {
  return haploMath::logdiff(logS(height,width),logT(width)) - log(height);
}

double RRMemo::log_continue_factor(int64_t totwidth) {
  return totwidth * log_continue_probability;
}

double RRMemo::log_recombination_penalty() {
  return rho;
}

double RRMemo::log_population_size() {
  return log(population_size);
}

} // namespace haploMath

} // namespace haplo
