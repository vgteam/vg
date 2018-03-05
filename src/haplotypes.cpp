#include "haplotypes.hpp"
#include "path.hpp"
#include "position.hpp"

namespace haplo {

// By default, should we warn when haplotype scoring fails?
bool warn_on_score_fail = false;

/*******************************************************************************
haplo_DP_edge_memo
*******************************************************************************/

haplo_DP_edge_memo::haplo_DP_edge_memo() : 
  in(vector<vg::Edge>(0)), out(vector<vg::Edge>(0)) {
  
}

haplo_DP_edge_memo::haplo_DP_edge_memo(xg::XG& graph, 
                                       xg::XG::ThreadMapping last_node,
                                       xg::XG::ThreadMapping node) {
  if(has_edge(graph, last_node, node)) {
    out = last_node.is_reverse ? 
         graph.edges_on_start(last_node.node_id) : 
         graph.edges_on_end(last_node.node_id);
    in = node.is_reverse ? 
        graph.edges_on_end(node.node_id) :
        graph.edges_on_start(node.node_id);
  } else {
    out = vector<vg::Edge>(0);
    in = vector<vg::Edge>(0);
  }
}

const vector<vg::Edge>& haplo_DP_edge_memo::edges_in() const {
  return in;
}

const vector<vg::Edge>& haplo_DP_edge_memo::edges_out() const {
  return out;
}

bool haplo_DP_edge_memo::is_null() const {
  return out.size() == 0;
}

bool haplo_DP_edge_memo::has_edge(xg::XG& graph, xg::XG::ThreadMapping old_node, xg::XG::ThreadMapping new_node) {
  vg::Edge edge_taken = xg::make_edge(old_node.node_id, old_node.is_reverse, new_node.node_id, new_node.is_reverse);

  bool edge_found = false;

  const vector<vg::Edge>& edges = old_node.is_reverse ? graph.edges_on_start(old_node.node_id) : 
                                                   graph.edges_on_end(old_node.node_id);
  
  for(auto& edge : edges) {
    // Look at every edge in order.
    if(xg::edges_equivalent(edge, edge_taken)) {
      // If we found the edge we're taking, break.
      edge_found = true;
      break;
    }
  }
  return edge_found;  
}

/*******************************************************************************
hDP_graph_accessor
*******************************************************************************/

hDP_graph_accessor::hDP_graph_accessor(xg::XG& graph, 
                                       xg::XG::ThreadMapping new_node, 
                                       haploMath::RRMemo& memo) :
  graph(graph), edges(haplo_DP_edge_memo()), 
  old_node(xg::XG::ThreadMapping()), new_node(new_node), memo(memo) {
    
}

hDP_graph_accessor::hDP_graph_accessor(xg::XG& graph, 
                                       xg::XG::ThreadMapping old_node, 
                                       xg::XG::ThreadMapping new_node, 
                                       haploMath::RRMemo& memo) :
  graph(graph), old_node(old_node), new_node(new_node), memo(memo),
  edges(haplo_DP_edge_memo(graph, old_node, new_node)) {
    
}

bool hDP_graph_accessor::has_edge() const {
  return !edges.is_null();
}

int64_t hDP_graph_accessor::new_side() const {
  return graph.id_to_rank(new_node.node_id) * 2 + new_node.is_reverse;
}

int64_t hDP_graph_accessor::new_height() const {
  return graph.node_height(new_node);
}

int64_t hDP_graph_accessor::old_height() const {
  return graph.node_height(old_node);
}

int64_t hDP_graph_accessor::new_length() const {
  return graph.node_length(new_node.node_id);
}

void hDP_graph_accessor::print(ostream& output_stream) const {
  output_stream << "From node: ID " << old_node.node_id << " is_reverse " << old_node.is_reverse << " ; To node: ID " << old_node.node_id << " is_reverse " << new_node.is_reverse << " ; Reference haplotypes visiting To Node: " << new_height() << endl;
}

/*******************************************************************************
hDP_gbwt_graph_accessor
*******************************************************************************/

gbwt_thread_t::gbwt_thread_t() {
}

gbwt_thread_t::gbwt_thread_t(const vector<gbwt::node_type>& nodes, const vector<size_t>& node_lengths) : nodes(nodes), node_lengths(node_lengths) {
}

void gbwt_thread_t::push_back(gbwt::node_type node, size_t node_length) {
  nodes.push_back(node);
  node_lengths.push_back(node_length);
}

gbwt::node_type& gbwt_thread_t::operator[](size_t i) {
  return nodes[i];  
}

const gbwt::node_type& gbwt_thread_t::operator[](size_t i) const {
  return nodes[i];  
}

gbwt::node_type& gbwt_thread_t::back() {
  return nodes.back();
}

const gbwt::node_type& gbwt_thread_t::back() const {
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

void haplo_DP_rectangle::extend(hDP_graph_accessor& ga) {
  int64_t new_side = ga.new_side();
  if(previous_index == -1) {
    // We're extending an empty state
    state.first = 0;
    state.second = ga.new_height();
  } else if(!ga.edges.is_null()) {
    state.first = ga.graph.where_to(flat_node, 
                                    state.first, 
                                    new_side, 
                                    ga.edges.edges_in(), 
                                    ga.edges.edges_out());
    state.second = ga.graph.where_to(flat_node, 
                                    state.second, 
                                    new_side, 
                                    ga.edges.edges_in(), 
                                    ga.edges.edges_out());
  } else {
    // gPBWT can't extend across an edge it doesn't know about; don't try
    state.first = 0;
    state.second = 0;
  }
  flat_node = new_side;
  inner_value = -1;
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
  auto r_0 = entries.at(0);
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
      for(i; i < entries.size(); i++) {
        double logLHS = memo.logT_base +
                        previous_R(i) +
                        memo.logT(length);
        entries[i]->R = haploMath::logsum(logLHS, logpS1S2RRS);
      }
    } else {
      for(i; i < entries.size(); i++) {
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
  return previous_values[(entries.at(i))->prev_idx()];
}

vector<double> haplo_DP_column::get_scores() const {
  vector<double> to_return;
  for(size_t i = 0; i < entries.size(); i++) {
    to_return.push_back(entries[i]->R);
  }
  return to_return;
}

vector<int64_t> haplo_DP_column::get_sizes() const {
  vector<int64_t> to_return;
  for(size_t i = 0; i < entries.size(); i++) {
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
    out << entries.at(i)->I() << "] : " << entries.at(i)->interval_size() << endl;
  }
}

bool haplo_DP_column::is_empty() const {
  return entries.size() == 0;
}

/*******************************************************************************
haplo_DP
*******************************************************************************/
haplo_score_type haplo_DP::score(const vg::Path& path, xg::XG& graph, haploMath::RRMemo& memo) {
  return score(path_to_thread_t(path), graph, memo);
}

haplo_score_type haplo_DP::score(const thread_t& thread, xg::XG& graph, haploMath::RRMemo& memo) {
  hDP_graph_accessor ga_i(graph, thread[0], memo);
  haplo_DP hdp(ga_i);
#ifdef debug
  cerr << "After entry 0 (" << thread[0].node_id << ") height: " << ga_i.new_height() << " score: " << hdp.DP_column.current_sum() << endl;
#endif
  if(ga_i.new_height() == 0) {
    if (warn_on_score_fail) {
      cerr << "[WARNING] Initial node in path is visited by 0 reference haplotypes" << endl;
      cerr << "Cannot compute a meaningful haplotype likelihood score" << endl;
      ga_i.print(cerr);
    }
    return pair<double, bool>(nan(""), false);
  }
  for(size_t i = 1; i < thread.size(); i++) {
    hDP_graph_accessor ga(graph, thread[i-1], thread[i], memo);
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
  cerr << "After entry " << i << " (" << thread[i].node_id << ") height: " << ga.new_height() << " score: " << hdp.DP_column.current_sum() << endl;
#endif
  }
  return pair<double, bool>(hdp.DP_column.current_sum(), true);
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

int64_t linear_haplo_structure::get_SNP_ref_position(size_t node_id) const {
  vector<vg::Edge> lnbr_edges = xg_index.edges_on_start(node_id);
  int64_t lnbr = lnbr_edges[0].from();
  vector<vg::Edge> SNP_allele_edges = xg_index.edges_on_end(lnbr);
  for(size_t i = 0; i < SNP_allele_edges.size(); i++) {
    if(xg_index.path_contains_node(xg_index.path_name(xg_ref_rank), SNP_allele_edges[i].to())) {
      return position_assuming_acyclic(SNP_allele_edges[i].to());
    }
  }
  throw runtime_error("no ref allele at SNP");
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
  char allele_char = xg_index.node_sequence(node_id).at(0);
  return allele::from_char(allele_char);
}

size_t linear_haplo_structure::SNVvector::ref_position(size_t i) const {
  return ref_positions.at(i);
}

alleleValue linear_haplo_structure::SNVvector::allele(size_t i) const {
  return alleles.at(i);
}

size_t linear_haplo_structure::SNVvector::size() const {
  return nodes.size();
}

bool linear_haplo_structure::sn_deletion_between_ref(int64_t left, int64_t right) const {
  int64_t gap = position_assuming_acyclic(right) - position_assuming_acyclic(left) - xg_index.node_length(left);
  if(gap == 0) {
    return false;
  } else if(gap == 1) {
    return true;
  } else {
    throw linearUnrepresentable("indel not of length 1");
  }
}

int64_t linear_haplo_structure::get_ref_following(int64_t node_id) const {
  vector<vg::Edge> r_edges = xg_index.edges_on_end(node_id);
  vector<int64_t> refs;
  for(size_t i = 0; i < r_edges.size(); i++) {
    if(xg_index.path_contains_node(xg_index.path_name(xg_ref_rank), r_edges[i].to())) {
      refs.push_back(r_edges[i].to());
    }
  }
  if(refs.size() == 0) {
    throw runtime_error("no ref node following");
  }
  size_t smallest = SIZE_MAX;
  int64_t node = refs[0];
  for(size_t i = 0; i < refs.size(); i++) {
    if(position_assuming_acyclic(refs[i]) < smallest) {
      smallest = position_assuming_acyclic(refs[i]);
      node = refs[i];
    }
  }
  return node;
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
    last_pos = position_assuming_acyclic(last_node);
  } else if(last_type == snv) {
    last_pos = get_SNP_ref_position(last_node);
    to_return.push_back(get_SNV_allele(last_node), last_pos, false);
  } else {
    throw linearUnrepresentable("not an SNV");
  }

  for(size_t i = 1; i < path.mapping_size(); i++) {
    int64_t this_node = path_mapping_node_id(path, i);
    nodeType this_type = get_type(this_node);
    size_t this_pos;

    if(this_type == invalid) {
      throw linearUnrepresentable("not an SNV");
    } else if(this_type == snv) {
      this_pos = get_SNP_ref_position(this_node);
      if(this_pos != last_pos + xg_index.node_length(last_node)) {
        throw linearUnrepresentable("indel immediately before SNV");
      }
      to_return.push_back(get_SNV_allele(this_node), this_pos, false);
    } else {
      this_pos = position_assuming_acyclic(this_node);
      if(last_type == snv) {
        if(this_pos != last_pos + xg_index.node_length(last_node)) {
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
  if(!xg_index.path_contains_node(xg_index.path_name(xg_ref_rank), node_id)) {
    throw runtime_error("requested position-in-path of node not in path");
  }
  return xg_index.position_in_path(node_id, xg_ref_rank).at(0);
}


bool linear_haplo_structure::is_solitary_ref(int64_t node_id) const {
  if(!xg_index.path_contains_node(xg_index.path_name(xg_ref_rank), node_id)) {
    return false;
  }
  
  vector<vg::Edge> l_edges = xg_index.edges_on_start(node_id);
  vector<vg::Edge> r_edges = xg_index.edges_on_end(node_id);
  for(size_t i = 0; i < l_edges.size(); i++) {
    if(xg_index.edges_on_end(l_edges[i].from()).size() != 1) {
      bool is_deletion_neighbour = true;
      vector<vg::Edge> neighbour_r_edges = xg_index.edges_on_end(l_edges[i].from());
      for(size_t i = 0; i < neighbour_r_edges.size(); i++) {
        if(neighbour_r_edges[i].to() != node_id) {
          vector<vg::Edge> neighbour_rr_edges = xg_index.edges_on_end(neighbour_r_edges[i].to());
          if(!(neighbour_rr_edges.size() == 1 && neighbour_rr_edges[0].to() == node_id)) {
            is_deletion_neighbour = false;
          }
        }
      }
      if(!is_deletion_neighbour) {
        return false;
      }
    }
  }
  for(size_t i = 0; i < r_edges.size(); i++) {
    if(xg_index.edges_on_start(r_edges[i].to()).size() != 1) {
      bool is_deletion_neighbour = true;
      vector<vg::Edge> neighbour_l_edges = xg_index.edges_on_start(r_edges[i].to());
      for(size_t i = 0; i < neighbour_l_edges.size(); i++) {
        if(neighbour_l_edges[i].from() != node_id) {
          vector<vg::Edge> neighbour_ll_edges = xg_index.edges_on_start(neighbour_l_edges[i].from());
          if(!(neighbour_ll_edges.size() == 1 && neighbour_ll_edges[0].from() == node_id)) {
            is_deletion_neighbour = false;
          }
        }
      }
      if(!is_deletion_neighbour) {
        return false;
      }
    }
  }
  return true;
}

bool linear_haplo_structure::is_snv(int64_t node_id) const {  
  // has only one left and one right neighbour
  int64_t lnbr;
  int64_t rnbr;
  
  vector<vg::Edge> l_edges = xg_index.edges_on_start(node_id);
  vector<vg::Edge> r_edges = xg_index.edges_on_end(node_id);

  if(l_edges.size() == 1 && r_edges.size() == 1) {
    lnbr = l_edges[0].from();
    rnbr = r_edges[0].to();
  } else {
    // has too many or too few neighbours to be an SNV
    return false;
  }
  
  vector<vg::Edge> lnbr_edges = xg_index.edges_on_end(lnbr);
  vector<vg::Edge> rnbr_edges = xg_index.edges_on_start(rnbr);
  if(lnbr_edges.size() != rnbr_edges.size()) {
    // neigbours must have a node not in common
    return false;
  }  
  
  vector<int64_t> r_of_lnbr(lnbr_edges.size());
  vector<int64_t> l_of_rnbr(rnbr_edges.size());
  for(size_t i = 0; i < lnbr_edges.size(); i++) {
    r_of_lnbr[i] = lnbr_edges[i].to();
  }
  for(size_t i = 0; i < rnbr_edges.size(); i++) {
    l_of_rnbr[i] = rnbr_edges[i].from();
  }
  for(size_t i = 0; i < r_of_lnbr.size(); i++) {
    if(xg_index.node_length(r_of_lnbr[i]) != 1) {
      if(r_of_lnbr[i] != rnbr) {
        return false;
      }
    }
  }
    
  // given guarantee that edge_sets contain no duplicates, check that there is an injective map 
  // from r_of_lnbr to l_of_rnbr. If so, must be bijection as are finite sets of same size
  for(size_t i = 0; i < r_of_lnbr.size(); i++) {
    if(r_of_lnbr[i] == rnbr) {
      r_of_lnbr[i] = -1;
    } else {
      for(size_t j = 0; j < l_of_rnbr.size(); j++) {
        if(r_of_lnbr[i] == l_of_rnbr[j]) {
          r_of_lnbr[i] = -1;
        }
      }
    }
  }
  for(size_t i = 0; i < r_of_lnbr.size(); i++) {
    if(r_of_lnbr[i] != -1) {
      return false;
    }
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
  
  size_t start = position_assuming_acyclic(path_mapping_node_id(path, 0));
  size_t length = 0;
  for(size_t i = 0; i < path.mapping_size(); i++) {
    vg::Mapping mapping = path.mapping(i);
    length += vg::mapping_from_length(mapping);
  }  
  
  if(SNV_candidates.size() == 0) {
    return new inputHaplotype(vector<alleleValue>(0), vector<size_t>(0), index, start, length);
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
  
  size_t last_i = index->get_site_index(SNV_candidates.ref_position(0));
  for(size_t i = 1; i < SNV_candidates.size(); i++) {
    size_t this_i = index->get_site_index(SNV_candidates.ref_position(i)); 
    if(this_i != last_i + 1) {
      return new inputHaplotype();
    } else {
      last_i = this_i;
    }
  }
  
  inputHaplotype* to_return = new inputHaplotype(alleles, vector<size_t>(positions.size(), 0), index, start, length);
  return to_return;
}

linear_haplo_structure::linear_haplo_structure(istream& slls_index, double log_mut_penalty, double log_recomb_penalty, xg::XG& xg_index, size_t xg_ref_rank) : xg_index(xg_index), xg_ref_rank(xg_ref_rank) {
  if(xg_ref_rank > xg_index.max_path_rank()) {
    throw runtime_error("reference path rank out of bounds");
  }
  index = new siteIndex(slls_index);
  cohort = new haplotypeCohort(slls_index, index);
  penalties = new penaltySet(log_recomb_penalty, log_mut_penalty, cohort->get_n_haplotypes());
}

linear_haplo_structure::~linear_haplo_structure() {
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
XGScoreProvider
*******************************************************************************/

XGScoreProvider::XGScoreProvider(xg::XG& index) : index(index) {
  // Nothing to do!
}

pair<double, bool> XGScoreProvider::score(const vg::Path& path, haploMath::RRMemo& memo) {
  return haplo_DP::score(path, index, memo);
}

/*******************************************************************************
LinearScoreProvider
*******************************************************************************/

LinearScoreProvider::LinearScoreProvider(const linear_haplo_structure& index) : index(index) {
  // Nothing to do!
}

pair<double, bool> LinearScoreProvider::score(const vg::Path& path, haploMath::RRMemo& memo) {
  // Memo is ignored; all penalties come from the index itself.
  return index.score(path);
}

/*******************************************************************************
path conversion
*******************************************************************************/

thread_t path_to_thread_t(const vg::Path& path) {
  thread_t t;
  for(size_t i = 0; i < path.mapping_size(); i++) {
    vg::Mapping mapping = path.mapping(i);
    auto pos = mapping.position();
    xg::XG::ThreadMapping m = {pos.node_id(), pos.is_reverse()};
    t.push_back(m);
  }
  return t;
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
  for(int i = 0; i < population_size; i++) {
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
  return (width-1)*logS_bases[height-1]; //logS_base = log(1 + i*exp_rho)
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
