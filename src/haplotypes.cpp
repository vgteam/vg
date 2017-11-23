#include "haplotypes.hpp"
#include "path.hpp"
#include "position.hpp"

/*******************************************************************************
haplo_DP_edge_memo
*******************************************************************************/

haplo_DP_edge_memo::haplo_DP_edge_memo() : 
  in(vector<vg::Edge>(0)), out(vector<vg::Edge>(0)) {
  
}

haplo_DP_edge_memo::haplo_DP_edge_memo(xg::XG& graph, 
                                       xg::XG::ThreadMapping last_node,
                                       xg::XG::ThreadMapping node) {
   out = last_node.is_reverse ? 
         graph.edges_on_start(last_node.node_id) : 
         graph.edges_on_end(last_node.node_id);
   in = node.is_reverse ? 
        graph.edges_on_end(node.node_id) :
        graph.edges_on_start(node.node_id);
}

const vector<vg::Edge>& haplo_DP_edge_memo::edges_in() const {
  return in;
}

const vector<vg::Edge>& haplo_DP_edge_memo::edges_out() const {
  return out;
}

/*******************************************************************************
hDP_graph_accessor
*******************************************************************************/

hDP_graph_accessor::hDP_graph_accessor(xg::XG& graph, 
                                       xg::XG::ThreadMapping new_node, 
                                       RRMemo& memo) :
  graph(graph), edges(haplo_DP_edge_memo()), 
  old_node(xg::XG::ThreadMapping()), new_node(new_node), memo(memo) {
    
}

hDP_graph_accessor::hDP_graph_accessor(xg::XG& graph, 
                                       xg::XG::ThreadMapping old_node, 
                                       xg::XG::ThreadMapping new_node, 
                                       RRMemo& memo) :
  graph(graph), edges(haplo_DP_edge_memo(graph, old_node, new_node)), 
  old_node(old_node), new_node(new_node), memo(memo) {
    
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


/*******************************************************************************
haplo_DP_rectangle
*******************************************************************************/

void haplo_DP_rectangle::extend(hDP_graph_accessor& ga) {
  int64_t new_side = ga.new_side();
  if(state.current_side == 0) {
    // We're extending an empty state
    state.range_start = 0;
    state.range_end = ga.new_height();
  } else {
    state.range_start = ga.graph.where_to(state.current_side, 
                                          state.range_start, 
                                          new_side, 
                                          ga.edges.edges_in(), 
                                          ga.edges.edges_out());
    state.range_end = ga.graph.where_to(state.current_side, 
                                        state.range_end, 
                                        new_side, 
                                        ga.edges.edges_in(), 
                                        ga.edges.edges_out());
  }
  state.current_side = new_side;
  inner_value = -1;
}

void haplo_DP_rectangle::false_extend(hDP_graph_accessor& ga, 
                                      int_itvl_t delta) {
  state.current_side = ga.new_side();
  state.range_start -= delta.bottom;
  state.range_end -= delta.top;
}

void haplo_DP_rectangle::calculate_I(int64_t succ_o_val) {
  inner_value = interval_size() - succ_o_val;
}

int64_t haplo_DP_rectangle::I() const {
  return inner_value;
}

int64_t haplo_DP_rectangle::interval_size() const {
  return state.range_end - state.range_start;
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
  //TODO double-check this
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

haplo_DP_column::haplo_DP_column(hDP_graph_accessor& ga) {
  haplo_DP_rectangle* first_rectangle = new haplo_DP_rectangle();
  entries = {first_rectangle};
  first_rectangle->extend(ga);
  update_inner_values();
  update_score_vector(ga.memo);
}

void haplo_DP_column::standard_extend(hDP_graph_accessor& ga) {
  previous_values = get_scores();
  previous_sizes = get_sizes();
  haplo_DP_rectangle* new_rectangle = new haplo_DP_rectangle();
  new_rectangle->extend(ga);
  vector<haplo_DP_rectangle*> new_entries = {new_rectangle};
  for(size_t i = 0; i < entries.size(); i++) {
    // extend candidate
    haplo_DP_rectangle* candidate = entries[i];
    candidate->extend(ga);
    candidate->set_prev_idx(i);
    // check if the last rectangle added is nonempty
    if(candidate->interval_size() == new_entries.back()->interval_size()) {
      haplo_DP_rectangle* to_delete = new_entries.back();
      new_entries.pop_back();
      if(to_delete->prev_idx() != -1) {
        entries[to_delete->prev_idx()] = nullptr;
      }
      delete to_delete;
    }
    if(candidate->interval_size() != 0) {
      new_entries.push_back(candidate);
    } else {
      break;
    }
  }
  entries = new_entries;
  update_inner_values();
}

void haplo_DP_column::update_inner_values() {
  for(size_t i = 0; i < entries.size() - 1; i++) {
    entries[i]->calculate_I(entries[i+1]->interval_size());
  }
  (entries.back())->calculate_I(0);
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

void haplo_DP_column::extend(hDP_graph_accessor& ga) {
  standard_extend(ga);
  length = (double)(ga.new_length());
  update_score_vector(ga.memo);
}

void haplo_DP_column::update_score_vector(RRMemo& memo) {
  haplo_DP_rectangle* r_0 = entries[0];
  if(entries.size() == 1 && entries[0]->prev_idx() == -1) {
    r_0->R = -log(memo.population_size);
    sum = r_0->R + log(r_0->interval_size());
    previous_values = {r_0->R};
    previous_sizes = {r_0->interval_size()};
    return;
  } else {
    previous_sum = sum;
    int64_t offset = (int64_t)(entries[0]->is_new());
    // TODO optimize by using arrays
    vector<double> continuing_Rs(entries.size() - offset);
    vector<int64_t> continuing_counts(entries.size() - offset);
    for(size_t i = offset; i < entries.size(); i++) {
      continuing_Rs[i - offset] = previous_values[entries[i]->prev_idx()];
      continuing_counts[i] = entries[i]->I();
    }
    
    double logpS1S2RRS = previous_sum + 
                         memo.log_recombination_penalty() + 
                         memo.logS(r_0->interval_size(), length);
    double logS1 = int_weighted_sum(continuing_Rs, continuing_counts);
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
        entries[i]->R = logsum(logLHS, logpS1S2RRS);
      }
    } else {
      for(i; i < entries.size(); i++) {
        double logLHS = memo.logT_base +
                        logsum(logS1RRD, previous_R(i) + memo.logT(length));
        entries[i]->R = logsum(logLHS, logpS1S2RRS);
      }
    }
  }
  previous_values = get_scores();
  previous_sizes = get_sizes();
  sum = int_weighted_sum(previous_values, previous_sizes);
}

double haplo_DP_column::previous_R(size_t i) const {
  return previous_values[(entries[i])->prev_idx()];
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

/*******************************************************************************
haplo_DP
*******************************************************************************/

haplo_DP::haplo_DP(hDP_graph_accessor& ga) : DP_column(haplo_DP_column(ga)) {
  
}

double haplo_DP::score(const vg::Path& path, xg::XG& graph, RRMemo& memo) {
  return score(path_to_thread_t(path), graph, memo);
}

double haplo_DP::score(const thread_t& thread, xg::XG& graph, RRMemo& memo) {
  hDP_graph_accessor ga_i(graph, thread[0], thread[0], memo);
  haplo_DP hdp(ga_i);
  for(size_t i = 1; i < thread.size(); i++) {
    hDP_graph_accessor ga(graph, thread[i-1], thread[i], memo);
    hdp.DP_column.extend(ga);
  }
  return hdp.DP_column.current_sum();
}

void haplo_DP::print(ostream& out) const {
  //TODO
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

RRMemo::RRMemo(double recombination_penalty, size_t population_size) : 
    rho(recombination_penalty), population_size(population_size) {
  exp_rho = exp(-rho);
  assert(exp_rho < 1);
  continue_probability = 1 - population_size * exp_rho;
  exp_rho = exp_rho/continue_probability;
  rho = -log(exp_rho);
  S.push_back(std::vector<double>(1, 1.0));
  S_multipliers.push_back(1.0);
  T.push_back(1.0);
  T_multiplier = 1.0 - exp_rho;

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

double int_weighted_sum(vector<double> values, vector<int64_t> counts) {
  if(values.size() == 1) {
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

  return logdiff(logS(height,width),logT(width)) - log(height);
}

double RRMemo::log_continue_factor(int64_t totwidth) {
  return totwidth * log(continue_probability);
}

double RRMemo::continue_factor(int64_t totwidth) {
  return exp(log_continue_factor(totwidth));
}

double RRMemo::recombination_penalty() {
  return exp_rho;
}

double RRMemo::log_recombination_penalty() {
  return -rho;
}

double RRMemo::cont_probability() {
  return continue_probability;
}

double RRMemo::S_value(int height, int width) {

  while (S.size() < height) {
    S_multipliers.push_back(S_multipliers[S_multipliers.size() - 1] + exp_rho);
    S.push_back(std::vector<double>(1, 1.0));
  }
  std::vector<double>& S_row = S[height - 1];
  double S_multiplier = S_multipliers[height - 1];

  while (S_row.size() < width) {
    S_row.push_back(S_row[S_row.size() - 1] * S_multiplier);
  }

  return S_row[width - 1];
}

double RRMemo::T_value(int width) {

  while (T.size() < width) {
    T.push_back(T[T.size() - 1] * T_multiplier);
  }

  return T[width - 1];
}

double RRMemo::rr_diff(int height, int width) {

  if (height < 1 || width < 1) {
    cerr << "error:[RRMemo] height and width of recombination rectangle must be >= 1" << endl;
  }

  return (S_value(height, width) - T_value(width)) / height;
}

double RRMemo::rr_same(int height, int width) {

  if (height < 1 || width < 1) {
    cerr << "error:[RRMemo] height and width of recombination rectangle must be >= 1" << endl;
  }

  double T_val = T_value(width);
  return (S_value(height, width) - T_val) / height + T_val;
}

double RRMemo::rr_adj(int width) {

  if (width < 1) {
    cerr << "error:[RRMemo] height and width of recombination rectangle must be >= 1" << endl;
  }

  return T_value(width);
}

double RRMemo::rr_all(int height, int width) {

  if (height < 1 || width < 1) {
    cerr << "error:[RRMemo] height and width of recombination rectangle must be >= 1" << endl;
  }

  return S_value(height, width);
}

// Unmemoized implementations of same polynomials

double rr_diff(int height, int width, double recombination_penalty) {
  double exp_rho = exp(-recombination_penalty);
  return (pow(1.0 + (height - 1.0) * exp_rho, width - 1.0) - pow(1.0 - exp_rho, width - 1.0)) / height;
}

double rr_same(int height, int width, double recombination_penalty) {
  double exp_rho = exp(-recombination_penalty);
  double T_val = pow(1.0 - exp_rho, width - 1.0);
  return (pow(1.0 + (height - 1.0) * exp_rho, width - 1.0) - T_val) / height + T_val;
}

double rr_adj(int width, double recombination_penalty) {
  return pow(1.0 - exp(-recombination_penalty), width - 1.0);
}

double rr_all(int height, int width, double recombination_penalty) {
  double exp_rho = exp(-recombination_penalty);
  return exp_rho * pow(1.0 + (height - 1.0) * exp_rho, width - 1.0);
}