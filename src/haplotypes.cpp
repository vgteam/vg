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
  }
  return pair<double, bool>(hdp.DP_column.current_sum(), true);
}

haplo_DP_column* haplo_DP::get_current_column() {
  return &DP_column;
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
