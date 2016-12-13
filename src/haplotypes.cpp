#include <iostream>
#include "haplotypes.hpp"
#include "path.hpp"
#include "position.hpp"

using namespace std;
using namespace xg;

RRMemo::RRMemo() {
  initialize(9.0);
}

void RRMemo::initialize(double recombination_penalty) {
  rho = recombination_penalty;
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

RRMemo::RRMemo(double recombination_penalty)  {
  initialize(recombination_penalty);
}

RRMemo::~RRMemo(void) {

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

double RRMemo::logSN(vector<double> logRs, vector<int> Is) {
  if(logRs.size() == 1) {
    return logRs[0] + log(Is[0]);
  } else {
    double max_summand = logRs[0] + log(Is[0]);
    int max_index = 0;
    vector<double> summands;
    for(int i = 0; i < logRs.size(); i++){
      summands.push_back(logRs[i] + log(Is[i]));
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

// Rectangular decomposition building:

cross_section::cross_section(int64_t height, int b, XG::ThreadMapping node) : height(height), b_index(b), node(node) {

}

void rectangle::extend(XG::ThreadMapping next_node, XG& graph, vector<vg::Edge>& edges_in, vector<vg::Edge>& edges_out) {
  int64_t next_side = graph.id_to_rank(next_node.node_id) * 2 + next_node.is_reverse;
  if(state.current_side == 0) {
    // We're extending an empty state
    state.range_start = 0;
    state.range_end = graph.node_height(next_node);
  } else {
    // bool edge_exists = check_for_edges(graph.rank_to_id(state.current_side / 2),state.current_side % 2,
    //        next_node.node_id, next_node.is_reverse, graph);
    // Else, look at where the path goes to and apply the where_to function to
    // shrink the range down.
    //if(edge_exists) {
      state.range_start = graph.where_to(state.current_side, state.range_start, next_side, edges_in, edges_out);
      state.range_end = graph.where_to(state.current_side, state.range_end, next_side, edges_in, edges_out);
    //} else {
    //  state.range_end = state.range_start;
    //}
  }
  state.current_side = next_side;
}

int rectangle::get_next_J(XG::ThreadMapping next_node, XG& graph, vector<vg::Edge>& edges_in, vector<vg::Edge>& edges_out) {
  extend(next_node, graph, edges_in, edges_out);
  J = state.count();
  return state.count();
}

int rectangle::get_next_J(thread_t& extension, XG& graph, vector<vg::Edge>& edges_in, vector<vg::Edge>& edges_out) {
  if(extension.size() == 1) {
    return get_next_J(extension.back(), graph, edges_in, edges_out);
  } else {
    xg::XG::ThreadMapping second_last_node = extension.end()[-2];
    state.current_side = graph.id_to_rank(second_last_node.node_id) * 2 + second_last_node.is_reverse;
    extend(extension.back(), graph, edges_in, edges_out);
    J = state.count();
    return state.count();
  }
}

void rectangle::extend(XG::ThreadMapping next_node, XG& graph) {
  int64_t next_side = graph.id_to_rank(next_node.node_id) * 2 + next_node.is_reverse;
  if(state.current_side == 0) {
    // We're extending an empty state
    state.range_start = 0;
    state.range_end = graph.node_height(next_node);
  } else {
    // bool edge_exists = check_for_edges(graph.rank_to_id(state.current_side / 2),state.current_side % 2,
    //        next_node.node_id, next_node.is_reverse, graph);
    // Else, look at where the path goes to and apply the where_to function to
    // shrink the range down.
    //if(edge_exists) {
      state.range_start = graph.where_to(state.current_side, state.range_start, next_side);
      state.range_end = graph.where_to(state.current_side, state.range_end, next_side);
    //} else {
    //  state.range_end = state.range_start;
    //}
  }
  state.current_side = next_side;
}

void rectangle::simple_extend(thread_t& extension, xg::XG& graph, int delta_start = 0, int delta_end = 0) {
  if(extension.size() > 0) {
    xg::XG::ThreadMapping next_node = extension.back();
    int64_t next_side = graph.id_to_rank(next_node.node_id) * 2 + next_node.is_reverse;
    state.current_side = next_side;
  }
  state.range_start -= delta_start;
  state.range_end -= delta_end;
}

void rectangle::simple_extend(xg::XG::ThreadMapping next_node, xg::XG& graph, int delta_start = 0, int delta_end = 0) {
  int64_t next_side = graph.id_to_rank(next_node.node_id) * 2 + next_node.is_reverse;
  state.current_side = next_side;
  state.range_start -= delta_start;
  state.range_end -= delta_end;
}

int rectangle::get_next_J(XG::ThreadMapping next_node, XG& graph) {
  extend(next_node, graph);
  J = state.count();
  return state.count();
}

int rectangle::get_next_J(thread_t& extension, XG& graph) {
  if(extension.size() == 1) {
    return get_next_J(extension.back(), graph);
  } else {
    xg::XG::ThreadMapping second_last_node = extension.end()[-2];
    state.current_side = graph.id_to_rank(second_last_node.node_id) * 2 + second_last_node.is_reverse;
    extend(extension.back(), graph);
    J = state.count();
    return state.count();
  }
}

inline xg::XG::ThreadMapping cross_section::get_node() {
  return node;
}

inline bool cross_section::has_squashed_nodes() {
  return (bridge.size() != 0);
}

inline xg::XG::ThreadMapping cross_section::get_last_node() {
  if(has_squashed_nodes()) {
    return bridge.back();
  } else {
    return get_node();
  }
}

bool check_for_edges(int64_t old_node_id, bool old_node_is_reverse, int64_t new_node_id,
          bool new_node_is_reverse, xg::XG& index) {
  // What edge are we following
  Edge edge_taken = make_edge(old_node_id, old_node_is_reverse, new_node_id, new_node_is_reverse);

  // Make sure we find it
  bool edge_found = false;

  vector<Edge> edges = new_node_is_reverse ? index.edges_on_end(new_node_id) : index.edges_on_start(new_node_id);

  for(auto& edge : edges) {
    // Look at every edge in order.
    if(edges_equivalent(edge, edge_taken)) {
      // If we found the edge we're taking, break.
      edge_found = true;
      break;
    }
  }
  if(edge_found == false) {  cerr << "did not find edge between" << old_node_id << " and " << new_node_id << endl;}
  return edge_found;
}

bool check_if_thread_t_broken(const thread_t& t, XG& graph) {
  bool broken = false;
  XG::ThreadMapping current_node = t[0];
  for(int i = 1; i < t.size(); i++) {
    XG::ThreadMapping next_node = t[i];
    bool edge_exists = check_for_edges(current_node.node_id, current_node.is_reverse,
            next_node.node_id, next_node.is_reverse, graph);
    if(!edge_exists) {
      broken = true;
      break;
    }
    current_node = next_node;
  }
  return broken;
}

pair<vector<pair<double,int>>,vector<pair<double,int>>> haplo_d::identify_structures(double leave_threshold,
        double return_threshold, int timeout, xg::XG& graph) {
  vector<pair<double,int>> recombinations;
  vector<pair<double,int>> private_polys;
  // timer ensures that you don't identify the "same" complex variant twice
  int timer = 0;
  for(int b = 1; b < cs.size()-1; b++) {
    if(timer < timeout) {
      timer += cs[b].width;
    } else {
      // does a lot of thread leave?
      int leavers = cs[b-1].S[0].J - cs[b].S[cs[b-1].S[0].next].J;
      double leave_ratio = leavers / cs[b-1].height;
      int starters = 0;
      if(cs[b].S[0].prev == -1) {
        starters = cs[b].S[0].I;
      }
      double return_ratio = 0;
      if(leavers > 0) {
        return_ratio = starters / leavers;
      }
      if(leave_ratio > leave_threshold) {
        if(return_ratio > return_threshold) {
          // returning "defines" recombination structurally--therefore we count it instead of leaving
          recombinations.push_back(make_pair(return_ratio,leavers));
        } else {
          private_polys.push_back(make_pair(leave_ratio,leavers));
        }
        // reset the timer
        timer = cs[b].width;
      }
    }
  }
  return make_pair(recombinations,private_polys);
}

haplo_d::haplo_d() {

}

// haplo_d::haplo_d(const thread_t& t, xg::XG& graph) {
//   initialize_skeleton(t, graph);
// }

haplo_d::haplo_d(const thread_t& t, XG& graph) {
  rectangle rect;
  rect.J = rect.get_next_J(t[0],graph);
  // At the leftmost node there is only one strip, so I = J
  rect.I = rect.J;
  int last_height = rect.J;
  cs.push_back(cross_section(rect.J,0,t[0]));
  cs.back().S.push_back(rect);
  int width = 0;
  int new_height;
  bool add_rectangle;
  bool add_A;
  for(int i = 1; i < t.size(); i++) {
    // Count the number of base pairs since the last entry or exit node
    width += graph.node_length(t[i-1].node_id);
    new_height = graph.node_height(t[i]);
    if(cs.back().S.size() != 0) {
      rect = cs.back().S[0];
      rect.J = rect.get_next_J(t[i],graph); // step this strip forward
      // Did any threads leave?
      if(last_height > rect.J) {
        add_A = 1;
      }
      // Are there any threads here which didn't come from the previous node?
      if(rect.J < new_height) {
        add_rectangle = 1;
        add_A = 1;
      }
      // This is an entry or exit node, add a cross-section to the vector of
      // cross-sections (which corresponds to the "A" set in the theory doc)
      if(add_A) {
        cs.back().width = width;
        width = 0;
        cs.push_back(cross_section(new_height,i,t[i]));
      } else {
        // This isn't a node where anything leaves or joins, let's skip over it
        cs.back().bridge.push_back(t[i]);
        for (size_t a = 0; a < cs.back().S.size(); a++) {
          cs.back().S[a].extend(t[i],graph);
        }
      }
      // This is an entry node; we also need a new rectangle corresponding to the
      // new strip. We need to do this *before* we populate since cross_sections
      // arrange rectangles newest -> oldest
      // NB that add_rectangle implies add_A
      if(add_rectangle) {
        rectangle new_rect;
        new_rect.extend(t[i],graph);
        new_rect.J = new_height;
        cs.back().height = new_rect.J;
        cs.back().S.push_back(new_rect);
        cs.back().S.back().I = new_rect.J - rect.J;
      }
      if(add_A) {
        int b = cs.size()-1;
        if(rect.J > 0) {
          cs[b].S.push_back(rect);
          cs[b].S.back().prev = 0;
          cs[b-1].S[0].next = cs[b].S.size()-1;
        }
      }
      last_height = new_height;
      add_A = 0;
      add_rectangle = 0;
    } else {
      cs.back().width = width;
      width = 0;
      cs.push_back(cross_section(new_height,i,t[i]));
      if(new_height > 0) {
        rectangle new_rect;
        new_rect.extend(t[i],graph);
        new_rect.J = new_height;
        cs.back().height = new_rect.J;
        cs.back().S.push_back(new_rect);
        cs.back().S.back().I = new_rect.J - rect.J;
      }
    }
  }
  if(cs.size() == 1) {
    cs.back().width = width;
  }
  cs.back().width += graph.node_length(t.back().node_id) - 1;
  for(int i = 0; i < cs.size(); i++) {
    tot_width += cs[i].width;
  }
}

void haplo_d::calculate_Is(XG& graph) {
  // node 0 was done in the haplo_d constructor; start at node 1
  for(int b = 1; b < cs.size(); b++) {
    // make sure that there is at least one rectangle here
    if(cs[b].S.size() != 0) {
      // get side and orientation of the next element in our query thread_t
      XG::ThreadMapping next_node = cs[b].get_node();
      // if J = 0 for a rectangle, then J must be 0 for all older rectangles
      bool nonempty_J = (cs[b].S.back().J > 0);
      if(nonempty_J) {
        int new_J;
        // start at a = 1 since the haplo_d initializer handles a = 0
        for(int a = 1; a < cs[b-1].S.size(); a++) {
          rectangle new_rect = cs[b-1].S[a];
          if (cs[b-1].bridge.size() != 0 &&
                new_rect.state.current_side != graph.id_to_rank(cs[b-1].bridge.back().node_id) * 2 + cs[b-1].bridge.back().is_reverse) {
            for(int i = 0; i < cs[b-1].bridge.size(); i++) {
              new_rect.extend(cs[b-1].bridge[i],graph);
            }
          }
          new_J = new_rect.get_next_J(next_node,graph);
          new_rect.J = new_J;
          if(new_J != 0) {
            if(new_J == cs[b].S.back().J) {
              // There was nothing in the last rectangle, so let's remove it
              cs[b].S.pop_back();
            } else {
              // There were threads in the last rectangle, so let's count them
              cs[b].S.back().I = cs[b].S.back().J - new_J;
              cs[b-1].S[a].next = cs[b].S.size();
            }
            // Let's assume for now that there is something in this rectangle;
            // if it turns out at the next step that we're wrong, then we'll
            // remove it
            cs[b].S.push_back(new_rect);
            cs[b].S.back().prev = a;
          } else {
            // don't try to add more rectangles, we know that J = 0 from here on
            break;
          }
        }
      } else {
        // this shouldn't be here
        cs[b].S.pop_back();
      }
      cs[b].S.back().I = cs[b].S.back().J;
    }
  }
}

// Querying haplo_d structure:

thread_t path_to_thread_t(Path& path) {
  thread_t t;
  for(size_t i = 0; i < path.mapping_size(); i++) {
    Mapping mapping = path.mapping(i);
    auto pos = mapping.position();
    XG::ThreadMapping m = {pos.node_id(), pos.is_reverse()};
    t.push_back(m);
  }
  return t;
}

// Calculating probabilitiees

inline double haplo_d::prev_R(int b, int a) {
  if(cs[b].S[a].prev == -1) {
    return 0;
  } else {
    return cs[b-1].S[cs[b].S[a].prev].R;
  }
}

inline double haplo_d::prev_logR(int b, int a) {
  return cs[b-1].S[cs[b].S[a].prev].logR;
}

inline int haplo_d::prev_I(int b, int a) {
  if(cs[b].S[a].prev == -1) {
    return 0;
  } else {
    return cs[b-1].S[cs[b].S[a].prev].I;
  }
}

vector<double> haplo_d::prev_logRs(int b) {
  vector<double> returnRs;
  for(int i = 0; i < cs[b].S.size(); i++) {
    if(cs[b].S[i].prev != -1) {
      returnRs.push_back(prev_logR(b,i));
    }
  }
  return returnRs;
}

vector<int> haplo_d::prev_Is(int b) {
  vector<int> returnIs;
  for(int i = 0; i < cs[b].S.size(); i++) {
    if(cs[b].S[i].prev != -1) {
      returnIs.push_back(prev_I(b,i));
    }
  }
  return returnIs;
}

vector<double> haplo_d::current_logRs(int b) {
  vector<double> returnRs;
  for(int i = 0; i < cs[b].S.size(); i++) {
    returnRs.push_back(cs[b].S[i].logR);
  }
  return returnRs;
}

vector<int> haplo_d::current_Is(int b) {
  vector<int> returnIs;
  for(int i = 0; i < cs[b].S.size(); i++) {
    returnIs.push_back(cs[b].S[i].I);
  }
  return returnIs;
}


double haplo_d::probability(RRMemo& memo) {
  // Make sure that we have nonempty rectangles through which to propagate our
  // calculation
  for (int b = 0; b < cs.size(); b++) {
    if(cs[b].height == 0) {
      return 0;
    }
  }
  // defined same as in writeup
  double S1 = 0;
  double S1S2 = 0;
  // compute R for the first interval (which has no predecessor)
  // NB: we are always working at the left edge of a cross_section
  cerr << "starting probability calculation. Tot length " << tot_width << ", |A| "<< cs.size() << endl;
  cs[0].S[0].R = memo.rr_all(cs[0].height,cs[0].width);
  for (int b = 1; b < cs.size(); b++) {
    S1 = 0;
    S1S2 = 0;
    if(cs[b].width != 1) {
      for(int a = 0; a < cs[b].S.size(); a++) {
        // N.B. that R's are r^a_b's rather that R^a_b's. Thus the I factor
        S1 += (prev_R(b,a)) * (prev_I(b,a));
      }
    }
    for(int a = 0; a < cs[b-1].S.size(); a++) {
      S1S2 += cs[b-1].S[a].R * cs[b-1].S[a].I;
    }
    // calculate contributions from all continuing strips
    for(int a = 0; a < cs[b].S.size(); a++) {
      if(cs[b].S[a].prev != -1) {
        cs[b].S[a].R =
        (1 - memo.recombination_penalty()) * ((S1 * memo.rr_diff(cs[b].height, cs[b].width)) + (prev_R(b,a) * memo.rr_adj(cs[b].width)) +
        (memo.recombination_penalty() * S1S2 * memo.rr_all(cs[b].height,cs[b].width)));
        cs[b].S[a].R = cs[b].S[a].R * memo.continue_factor(cs[b].width);
      } else {
        cs[b].S[a].R = (memo.recombination_penalty() * S1S2 * memo.rr_all(cs[b].height,cs[b].width));
        cs[b].S[a].R = cs[b].S[a].R * memo.continue_factor(cs[b].width);
      }
    }
  }
  double total_probability_haplotype = 0;
  for(int a = 0; a < cs.back().S.size(); a++) {
    total_probability_haplotype += cs.back().S[a].R * cs.back().S[a].I;
  }
  return total_probability_haplotype;
}

double haplo_d::log_probability(RRMemo& memo) {
  // Make sure that we have nonempty rectangles through which to propagate our
  // calculation
  for (int b = 0; b < cs.size(); b++) {
    if(cs[b].height == 0 || cs[b].S.size() == 0) {
      return nan("");
    }
  }
  // defined same as in writeup
  double logS1 = 0;
  double logS1S2 = 0;
  double logpS1S2RRS = 0;
  double logS1RRD = 0;
  double logLHS = 0;
  // compute R for the first interval (which has no predecessor)
  // we are always working at the left edge of a cross_section
  cerr << "starting log-probability calculation. Tot length " << tot_width << ", |A| " << cs.size() << endl;
  cs[0].S[0].logR = memo.logS(cs[0].height,cs[0].width);
  for (int b = 1; b < cs.size(); b++) {
    logS1S2 = memo.logSN(current_logRs(b-1),current_Is(b-1));
    logpS1S2RRS = logS1S2 + memo.log_recombination_penalty() + memo.logS(cs[b].height,cs[b].width);
    if(prev_logRs(b).size() == 0) {
      return nan("");
    }
    logS1 = memo.logSN(prev_logRs(b),prev_Is(b));
    logS1RRD = logS1 + memo.logRRDiff(cs[b].height,cs[b].width);
    // calculate contributions from all continuing strips
    for(int a = 0; a < cs[b].S.size(); a++) {
      if(cs[b].S[a].prev != -1) {
        if(cs[b].width == 1) {
          logLHS = memo.logT_base + prev_logR(b,a) + memo.logT(cs[b].width);
        } else {
          logLHS = memo.logT_base + logsum(logS1RRD,prev_logR(b,a) + memo.logT(cs[b].width));
        }
        cs[b].S[a].logR = logsum(logLHS,logpS1S2RRS);
      } else {
        cs[b].S[a].logR = logpS1S2RRS;
      }
    }
  }
  double total_probability_haplotype = memo.logSN(current_logRs(cs.size()-1),current_Is(cs.size()-1));
  return total_probability_haplotype;
}

void logRR_tests(double recombination_penalty) {
  RRMemo memo = RRMemo(recombination_penalty);
  vector<int> Is = {4000,1000,500,10,5};
  vector<int> lengths = {1000,20,5,3,1};
  for(int i = 0; i < lengths.size(); i++) {

    cerr << "logsum(" << Is[i] << "," << lengths[i] << ") = " << logsum(log(Is[i]),log(lengths[i])) << " | ";
    cerr << "log(" << Is[i] << "+" << lengths[i] << ") = " << log(Is[i] + lengths[i]) << endl;
  }
  for(int i = 0; i < lengths.size(); i++) {
    cerr << "logdiff(" << Is[i] << "," << lengths[i] << ") = " << logdiff(log(Is[i]),log(lengths[i])) << " | ";
    cerr << "log(" << Is[i] << "-" << lengths[i] << ") = " << log(Is[i] - lengths[i]) << endl;
  }

  for(int i = 0; i < lengths.size(); i++) {
    for(int j = 0; j < Is.size(); j++) {
      cerr << "logDiff(" << Is[j] <<","<< lengths[i] << ") = " << memo.logRRDiff(Is[j],lengths[i]) << " | ";
      cerr << "log[Diff(" << Is[j] <<","<< lengths[i] << ")] = " << log(memo.rr_diff(Is[j],lengths[i])) << endl;
    }
  }
  vector<double> Rs = {0.00001,0.000002,0.000000003,0.0000000000004,0.000000000000000000007};
  vector<double> logRs;
  for(int i = 0; i < Rs.size(); i++) {
    logRs.push_back(log(Rs[i]));
    cerr << "exp(logR) = " << logRs[i] << " | " << "R = " << Rs[i] << endl;
  }
  cerr << "logSN = " << memo.logSN(logRs,Is) << " | " << "log[SN] = ";
  double result = 0;
  double intermediate = 0;
  for(int i = 0; i < Rs.size(); i++) {
    intermediate = Rs[i]*Is[i];
    result += intermediate;
  }
  cerr << log(result) << endl;
  cerr << "------------" << endl;
}

bool RR_tests(void) {
  // RRMemo tests
  double epsilon = 0.0000001;

  double memo_val;
  double direct_val;

  for (double rho = 1.0; rho < 5.0; rho += 1.0) {

    RRMemo memo(rho);

    for (int c = 1; c < 10; c++) {
      for (int n = 1; n < 10; n++) {

        memo_val = memo.rr_diff(n, c);
        direct_val = rr_diff(n, c, rho);

        if (fabs(memo_val - direct_val) > epsilon) {
          cerr << "FAIL: rr_diff, n = " << n << ", c = " << c << ", rho = " << rho
          << ", memo value = " << memo_val << ", direct value = " << direct_val << endl;
          exit(1);
        }

        memo_val = memo.rr_same(n, c);
        direct_val = rr_same(n, c, rho);

        if (fabs(memo_val - direct_val) > epsilon) {
          cerr << "FAIL: rr_same, n = " << n << ", c = " << c << ", rho = " << rho
          << ", memo value = " << memo_val << ", direct value = " << direct_val << endl;
          exit(1);
        }

        memo_val = memo.rr_all(n, c);
        direct_val = rr_all(n, c, rho);

        if (fabs(memo_val - direct_val) > epsilon) {
          cerr << "FAIL: rr_all, n = " << n << ", c = " << c << ", rho = " << rho
          << ", memo value = " << memo_val << ", direct value = " << direct_val << endl;
          exit(1);
        }
      }

      memo_val = memo.rr_adj(c);
      direct_val = rr_adj(c, rho);

      if (fabs(memo_val - direct_val) > epsilon) {
        cerr << "FAIL: rr_adj, c = " << c << ", rho = " << rho
        << ", memo value = " << memo_val << ", direct value = " << direct_val << endl;
        exit(1);
      }
    }
  }

  cerr << "RR tests passed!" << endl;
  return true;
}

void haplo_d::print(ostream& stream) {
  for(int a = 0; a < cs.size(); a++) {
    int joining = cs[a].S[0].prev == -1 ? cs[a].S[0].I : 0;
    int left = a == 0 ? 0 : cs[a-1].S[0].J - cs[a].S[cs[a-1].S[0].next].J;
    stream  << "[" << a << "]\t"      // position
            << cs[a].width << "\t"      // sequence length
            << cs[a].S[0].J << "\t"   // consistent haplotypes
            << cs[a].S.size() << "\t" // A_curr
            << joining << "\t"
            << left << endl;
  }
}

void haplo_d::print_detailed(ostream& stream) {
  for(int a = 0; a < cs.size(); a++) {
    stream  << "[" << a << "]\t"      // position
            << cs[a].get_node().node_id << "\t"
            << cs[a].width << "\t"      // sequence length
            << cs[a].S[0].J << "\t"   // consistent haplotypes
            << cs[a].S.size() << "\t";// A_curr
    for(int b = 0; b < cs[a].S.size(); b++) {
      stream << cs[a].S[b].J << "\t";
    }
    stream << endl;
  }
}

void haplo_d::print_detailed_searchstates(ostream& stream) {
  for(int a = 0; a < cs.size(); a++) {
    stream  << "[" << a << "]\t"      // position
            << cs[a].get_node().node_id << "\t"
            << cs[a].width << "\t"      // sequence length
            << cs[a].S[0].J << "\t"   // consistent haplotypes
            << cs[a].S.size() << "\t";// A_curr
    for(int b = 0; b < cs[a].S.size(); b++) {
      // stream << cs[a].S[b].J << "\t";
      // stream << cs[a].S[b].state.current_side << ":[";
      stream << "[";
      stream << cs[a].S[b].state.range_start << ",";
      stream << cs[a].S[b].state.range_end << "] (" << cs[a].S[b].prev << "," << cs[a].S[b].next << ")\t";
    }
    stream << endl;
  }
}

void haplo_d::print_graphical(ostream& stream) {
  for(int a = 0; a < cs.size(); a++) {

  }
}

void haplo_d::print_neighbours(ostream& stream) {
  for(int i = 0; i < cs.size(); i++) {
    stream << cs[i].get_node().node_id<< ", J's: \t";
    for(int j = 0; j < cs[i].S.size(); j++) {
      stream << cs[i].S[j].J << "\t\t";
    }
    stream << endl;
    stream << cs[i].get_node().node_id<< ", prev: \t";
    for(int j = 0; j < cs[i].S.size(); j++) {
      stream << cs[i].S[j].prev << "\t";
    }
    stream << endl;
    stream << cs[i].get_node().node_id<< ", next: \t";
    for(int j = 0; j < cs[i].S.size(); j++) {
      stream << cs[i].S[j].next << "\t";
    }
    stream << endl;
  }
}

thread_t extract_thread(xg::XG& index, xg::XG::ThreadMapping node, int64_t offset = 0, int64_t max_length = 0) {
  thread_t path;
  int64_t side = (node.node_id)*2 + node.is_reverse;
  bool continue_search = true;
  while(continue_search) {
    xg::XG::ThreadMapping m = {index.rank_to_id(side / 2), (bool) (side % 2)};
    path.push_back(m);
    // Work out where we go:
    // What edge of the available edges do we take?
    int64_t edge_index = index.bs_get(side, offset);
    assert(edge_index != index.BS_SEPARATOR);
    if(edge_index == index.BS_NULL) {
        // Path ends here.
        break;
    } else {
        // If we've got an edge, convert to an actual edge index
        edge_index -= 2;
    }
    // We also should not have negative edges.
    assert(edge_index >= 0);
    // Look at the edges we could have taken next
    vector<Edge> edges_out = side % 2 ? index.edges_on_start(index.rank_to_id(side / 2)) :
          index.edges_on_end(index.rank_to_id(side / 2));
    assert(edge_index < edges_out.size());
    Edge& taken = edges_out[edge_index];
    // Follow the edge
    int64_t other_node = taken.from() == index.rank_to_id(side / 2) ? taken.to() :
          taken.from();
    bool other_orientation = (side % 2) != taken.from_start() != taken.to_end();
    // Get the side
    int64_t other_side = index.id_to_rank(other_node) * 2 + other_orientation;
    // Go there with where_to
    offset = index.where_to(side, offset, other_side);
    side = other_side;
    // Continue the process from this new side
    if(max_length != 0) {
      if(path.size() >= max_length) {
        continue_search = false;
      }
    }
  }
  return path;
}

void haplo_d::log_calculate_Is(xg::XG& graph) {
  // Things which were calculated in the constructor:
  // -- A's
  // -- J for the top continuing and any new rectangle
  // -- I for any new rectangle
  // ostream& stream = cout;
  vector<Edge> edges_out;
  vector<Edge> edges_in;
  for(int b = 1; b < cs.size(); b++) {
    // xg::XG::ThreadMapping lastnode;
    // if(cs[b-1].bridge.size() == 0) {
    //   lastnode = cs[b-1].get_node();
    // } else {
    //   lastnode = cs[b-1].bridge.back();
    // }
    edges_out = cs[b-1].get_last_node().is_reverse ? graph.edges_on_start(cs[b-1].get_last_node().node_id) : graph.edges_on_end(cs[b-1].get_last_node().node_id);
    // edges_out = lastnode.is_reverse ? graph.edges_on_start(lastnode.node_id) : graph.edges_on_end(lastnode.node_id);
    edges_in = cs[b].get_node().is_reverse ? graph.edges_on_end(cs[b].get_node().node_id) : graph.edges_on_start(cs[b].get_node().node_id);
    vector<rectangle>& prevAs = cs[b-1].S;
    vector<rectangle>& currAs = cs[b].S;
    bool new_threads = (prevAs[0].next == 1);
    // make sure that there is at least one rectangle here
    if(prevAs.size() == 0) {
      cerr << "[vg haplo error] no consistent haplotypes at node " << cs[b-1].get_node().node_id << endl;
    } else if(prevAs.size() == 1) {
      currAs.back().I = currAs.back().J;
      // currAs has size at most 2
      if(currAs.size() == 2) {
        currAs[0].I = currAs[0].J - currAs[1].J;
      }
    } else if(prevAs.size() >= 2) {
      // We're going to have to extend, so let's grab the next node
      XG::ThreadMapping next_node = cs[b].get_node();
      // Let's also grab the nodes which we'll skip over between this and the last node
      thread_t extension = cs[b-1].bridge;
      extension.push_back(next_node);
      // if J = 0 for a rectangle, then J must be 0 for all older rectangles
      if(currAs.back().J == 0) {
        currAs.pop_back();
      } else {
        int deltaJ = prevAs[0].J - currAs[prevAs[0].next].J;
        if(deltaJ == 0) {
          // cerr << b << ", deltaJ = 0" << endl;
          currAs[prevAs[0].next].I = prevAs[0].I;
          int delta_start = prevAs[0].state.range_start - currAs[prevAs[0].next].state.range_start;
          int delta_end = prevAs[0].state.range_end - currAs[prevAs[0].next].state.range_end;
          for(int a = 1; a < prevAs.size(); a++) {
            rectangle new_rect = prevAs[a];
            new_rect.simple_extend(extension, graph, delta_start, delta_end);
            new_rect.prev = a;
            currAs.push_back(new_rect);
            prevAs[a].next = currAs.size()-1;
          }
        } else {
          // binaryI(XG&, thread_t, b, atop, abott,    dJtop, dJbott, Jtop,               Jbott, indent level)
          binaryI(graph, extension, b, 0, prevAs.size(), deltaJ, 0, currAs[prevAs[0].next].J, 0, 0, edges_in, edges_out);
          for(int a = 0; a < currAs.size() - 1; a++) {
            currAs[a].I = currAs[a].J - currAs[a+1].J;
          }
          currAs.back().I = currAs.back().J;
        }
      }
    }
  }
}

void haplo_d::seeded_log_calculate_Is(xg::XG& graph) {
  // Things which were calculated in the constructor:
  // -- A's
  // -- J for the top continuing and any new rectangle
  // -- I for any new rectangle
  // ostream& stream = cout;
  vector<Edge> edges_out;
  vector<Edge> edges_in;
  for(int b = 1; b < cs.size(); b++) {
    vector<rectangle>& prevAs = cs[b-1].S;
    vector<rectangle>& currAs = cs[b].S;
    // xg::XG::ThreadMapping lastnode;
    // if(cs[b-1].bridge.size() == 0) {
    //   lastnode = cs[b-1].get_node();
    // } else {
    //   lastnode = cs[b-1].bridge.back();
    // }
    // edges_out = lastnode.is_reverse ? graph.edges_on_start(lastnode.node_id) : graph.edges_on_end(lastnode.node_id);
    edges_out = cs[b-1].get_last_node().is_reverse ? graph.edges_on_start(cs[b-1].get_last_node().node_id) : graph.edges_on_end(cs[b-1].get_last_node().node_id);
    edges_in = cs[b].get_node().is_reverse ? graph.edges_on_end(cs[b].get_node().node_id) : graph.edges_on_start(cs[b].get_node().node_id);
    bool new_threads = (prevAs[0].next == 1);
    // make sure that there is at least one rectangle here
    if(prevAs.size() == 0) {
      cerr << "[vg haplo error] no consistent haplotypes at node " << cs[b-1].get_node().node_id << endl;
    } else if(prevAs.size() == 1) {
      currAs.back().I = currAs.back().J;
      // currAs has size at most 2
      if(currAs.size() == 2) {
        currAs[0].I = currAs[0].J - currAs[1].J;
      }
    } else if(prevAs.size() >= 2) {
      // We're going to have to extend, so let's grab the next node
      XG::ThreadMapping next_node = cs[b].get_node();
      // Let's also grab the nodes which we'll skip over between this and the last node
      thread_t extension = cs[b-1].bridge;
      extension.push_back(next_node);
      // if J = 0 for a rectangle, then J must be 0 for all older rectangles
      if(currAs.back().J == 0) {
        currAs.pop_back();
      } else {
        int deltaJ = prevAs[0].J - currAs[prevAs[0].next].J;
        if(deltaJ == 0) {
          // cerr << b << ", deltaJ = 0" << endl;
          currAs[prevAs[0].next].I = prevAs[0].I;
          int delta_start = prevAs[0].state.range_start - currAs[prevAs[0].next].state.range_start;
          int delta_end = prevAs[0].state.range_end - currAs[prevAs[0].next].state.range_end;
          for(int a = 1; a < prevAs.size(); a++) {
            rectangle new_rect = prevAs[a];
            new_rect.simple_extend(extension, graph, delta_start, delta_end);
            new_rect.prev = a;
            currAs.push_back(new_rect);
            prevAs[a].next = currAs.size()-1;
          }
        } else {
          vector<int> previously_big;
          int big_cutoff = 400;
          for(int i = 1; i < prevAs.size(); i++) {
            if(prevAs[i].I >= big_cutoff) {
              previously_big.push_back(i);
            }
          }
          // cerr << "made big list, it's size " << previously_big.size() << endl;
          vector<rectangle> big_rectangles;
          vector<int> big_deltas;
          vector<int> big_Js;
          for(int i = 0; i < previously_big.size(); i++) {
            big_rectangles.push_back(prevAs[previously_big[i]]);
            int Jbig = big_rectangles.back().get_next_J(extension,graph, edges_in, edges_out);
            // cerr << Jbig << "\t" << flush;
            big_Js.push_back(Jbig);
            big_deltas.push_back(prevAs[previously_big[i]].J - Jbig);
            if(Jbig == 0) {
              break;
            }
          }
          // cerr << endl;
          // cerr << "collected attributes of big rectangles, " << big_Js.size() << " are nonempty" << endl;
          if(big_Js.size() > 0) {
            int Aabove = 0;
            int Jabove = currAs[prevAs[0].next].J;
            int dJabove = prevAs[0].J - currAs[prevAs[0].next].J;
            for(int i = 0; i < previously_big.size(); i++) {
              if(big_Js[i] == Jabove) {
                // all rectangles between are actually empty
                prevAs[currAs.back().prev].next = -1;
                currAs.pop_back();
              } else {
                binaryI(graph, extension, b, Aabove, previously_big[i], dJabove, big_deltas[i], Jabove, big_Js[i], 0, edges_in, edges_out);
              }
              Aabove = previously_big[i];
              Jabove = big_Js[i];
              dJabove = big_deltas[i];
              if(big_Js[i] == 0) {
                // Don't build smaller rectangles
                // Don't add this rectangle
                break;
              } else {
                big_rectangles[i].prev = previously_big[i];
                currAs.push_back(big_rectangles[i]);
                prevAs[previously_big[i]].next = currAs.size()-1;
              }
            }
            if(big_Js.back() != 0) {
              binaryI(graph, extension, b, previously_big.back(), prevAs.size(), big_deltas.back(), 0, big_Js.back(), 0, 0, edges_in, edges_out);
            }
          } else {
            binaryI(graph, extension, b, 0, prevAs.size(), deltaJ, 0, currAs[prevAs[0].next].J, 0, 0, edges_in, edges_out);
          }
          for(int a = 0; a < currAs.size() - 1; a++) {
            currAs[a].I = currAs[a].J - currAs[a+1].J;
          }
          currAs.back().I = currAs.back().J;
        }
      }
    }
  }
}

void haplo_d::binaryI(xg::XG& graph, thread_t extension, int b, int atop, int abottom, int deltaJtop, int deltaJbottom, int Jtop, int Jbottom, int level, vector<vg::Edge>& edges_in, vector<vg::Edge>& edges_out) {
  vector<rectangle>& prevAs = cs[b-1].S;
  vector<rectangle>& currAs = cs[b].S;
  // for(int i = 0; i < level; i++) {
  //   cerr << "  ";
  // }
  // cerr << "called b-I," << b << "\t b = [" << atop << "," << abottom << "]\tJ = [" << Jtop << "," << Jbottom << "]" << endl;

  if(abottom <= atop + 1) {
    // You've reached max recursion depth
    return;
  } else if(deltaJtop == deltaJbottom) {
    // cerr << "simple_extend" << endl;
    // The nexting behavior of the J-counts ensures that there are no changes to
    // thread membership in the interval we're evaluating here. We can safely
    // extend everythings
    int delta_start = prevAs[atop].state.range_start - currAs[prevAs[atop].next].state.range_start;
    int delta_end = prevAs[atop].state.range_end - currAs[prevAs[atop].next].state.range_end;
    for(int i = atop + 1; i < abottom; i++) {
      rectangle rect = prevAs[i];
      rect.simple_extend(extension, graph, delta_start, delta_end);
      rect.prev = i;
      rect.J = rect.state.count();
      currAs.push_back(rect);
      prevAs[i].next = currAs.size()-1;
    }
  } else {
    // At least one thread diverges within this interval; let's find it
    int mid = atop + (abottom-atop)/2;
    rectangle rect_mid = prevAs[mid];
    int Jmid = rect_mid.get_next_J(extension, graph, edges_in, edges_out);
    int deltaJmid = prevAs[mid].J - Jmid;
    if(Jmid == Jtop) {
      // all rectangles between are actually empty
      prevAs[currAs.back().prev].next = -1;
      currAs.pop_back();
    } else {
      binaryI(graph, extension, b, atop, mid, deltaJtop, deltaJmid, Jtop, Jmid, level+1, edges_in, edges_out);
    }
    if(Jmid == 0) {
      // Don't build smaller rectangles
      // Don't add this rectangle
      return;
    } else {
      rect_mid.prev = mid;
      currAs.push_back(rect_mid);
      prevAs[mid].next = currAs.size()-1;
      if(Jbottom == Jmid) {
        prevAs[currAs.back().prev].next = -1;
        currAs.pop_back();
      } else {
        binaryI(graph, extension, b, mid, abottom, deltaJmid, deltaJbottom, Jmid, Jbottom, level+1, edges_in, edges_out);
      }
    }
  }
}

void haplo_d::binaryI(xg::XG& graph, thread_t extension, int b, int atop, int abottom, int deltaJtop, int deltaJbottom, int Jtop, int Jbottom, int level) {
  vector<rectangle>& prevAs = cs[b-1].S;
  vector<rectangle>& currAs = cs[b].S;
  // for(int i = 0; i < level; i++) {
  //   cerr << "  ";
  // }
  // cerr << "called b-I," << b << "\t b = [" << atop << "," << abottom << "]\tJ = [" << Jtop << "," << Jbottom << "]" << endl;

  if(abottom <= atop + 1) {
    // You've reached max recursion depth
    return;
  } else if(deltaJtop == deltaJbottom) {
    // cerr << "simple_extend" << endl;
    // The nexting behavior of the J-counts ensures that there are no changes to
    // thread membership in the interval we're evaluating here. We can safely
    // extend everythings
    int delta_start = prevAs[atop].state.range_start - currAs[prevAs[atop].next].state.range_start;
    int delta_end = prevAs[atop].state.range_end - currAs[prevAs[atop].next].state.range_end;
    for(int i = atop + 1; i < abottom; i++) {
      rectangle rect = prevAs[i];
      rect.simple_extend(extension, graph, delta_start, delta_end);
      rect.prev = i;
      rect.J = rect.state.count();
      currAs.push_back(rect);
      prevAs[i].next = currAs.size()-1;
    }
  } else {
    // At least one thread diverges within this interval; let's find it
    int mid = atop + (abottom-atop)/2;
    rectangle rect_mid = prevAs[mid];
    int Jmid = rect_mid.get_next_J(extension,graph);
    int deltaJmid = prevAs[mid].J - Jmid;
    if(Jmid == Jtop) {
      // all rectangles between are actually empty
      prevAs[currAs.back().prev].next = -1;
      currAs.pop_back();
    } else {
      binaryI(graph, extension, b, atop, mid, deltaJtop, deltaJmid, Jtop, Jmid, level+1);
    }
    if(Jmid == 0) {
      // Don't build smaller rectangles
      // Don't add this rectangle
      return;
    } else {
      rect_mid.prev = mid;
      currAs.push_back(rect_mid);
      prevAs[mid].next = currAs.size()-1;
      if(Jbottom == Jmid) {
        prevAs[currAs.back().prev].next = -1;
        currAs.pop_back();
      } else {
        binaryI(graph, extension, b, mid, abottom, deltaJmid, deltaJbottom, Jmid, Jbottom, level+1);
      }
    }
  }
}

vector<rectangle*> haplo_d::trace_strip(int offset) {
  return trace_strip(0, offset, -1);
}

vector<rectangle*> haplo_d::trace_strip(int a, int offset, int distance) {
  vector<rectangle*> to_return;
  assert(a < cs.size() && a >= 0);
  assert(offset < cs[a].S.size());
  if(offset < 0) {
    to_return.push_back(&cs[a].S.back());
  } else {
    to_return.push_back(&cs[a].S[offset]);
  }
  if(distance < 0 || distance > cs.size() - a) {
    distance = cs.size() - a;
  }
  for(int i = a + 1; i < distance + a; i++) {
    int to_add_index = to_return.back()->next;
    if(to_add_index < 0) {
      return to_return;
    } else {
      to_return.push_back(&cs[i].S[to_add_index]);
    }
  }
}

void haplo_d::build_start(xg::XG::ThreadMapping node, xg::XG& graph) {
  assert(cs.size() == 0);
  rectangle rect;
  rect.J = rect.get_next_J(node,graph);
  // At the leftmost node there is only one strip, so I = J
  rect.I = rect.J;
  int last_height = rect.J;
  cs.push_back(cross_section(rect.J,0,node));
  cs.back().S.push_back(rect);
}

void haplo_d::initialize_skeleton(thread_t& t, xg::XG& graph) {
  build_start(t[0],graph);
  initialize_skeleton(t, 1, cs[0],graph);
}

void haplo_d::initialize_skeleton(thread_t& t, int start, cross_section& prevAs, xg::XG& graph) {
  assert(cs.size() == 0);
  initialize_skeleton(t, make_pair(start, t.size()-1), prevAs, graph);
}

//TODO: translate this
void haplo_d::initialize_skeleton(thread_t& t, pair<int,int> interval, cross_section& prevAs, xg::XG& graph) {
  rectangle rect;
  int new_height;
  int last_height = prevAs.S[0].J;
  bool add_rectangle;
  bool add_A;
  //TODO: fix this
  int width = 0;
  for(int i = interval.first; i <= interval.second; i++) {
    // Count the number of base pairs since the last entry or exit node
    width += graph.node_length(t[i-1].node_id);
    new_height = graph.node_height(t[i]);
    if(cs.back().S.size() != 0) {
      if(i == interval.first) {
        prevAs.S[0];
      } else {
        rect = cs.back().S[0];
      }
      rect.J = rect.get_next_J(t[i],graph); // step this strip forward
      // Did any threads leave?
      if(last_height > rect.J) {
        add_A = 1;
      }
      // Are there any threads here which didn't come from the previous node?
      if(rect.J < new_height) {
        add_rectangle = 1;
        add_A = 1;
      }
      // This is an entry or exit node, add a cross-section to the vector of
      // cross-sections (which corresponds to the "A" set in the theory doc)
      if(add_A) {
        cs.back().width = width;
        width = 0;
        cs.push_back(cross_section(new_height,i,t[i]));
      } else {
        // This isn't a node where anything leaves or joins, let's skip over it
        cs.back().bridge.push_back(t[i]);
        for (size_t a = 0; a < cs.back().S.size(); a++) {
          cs.back().S[a].extend(t[i],graph);
        }
      }
      // This is an entry node; we also need a new rectangle corresponding to the
      // new strip. We need to do this *before* we populate since cross_sections
      // arrange rectangles newest -> oldest
      // NB that add_rectangle implies add_A
      if(add_rectangle) {
        rectangle new_rect;
        new_rect.extend(t[i],graph);
        new_rect.J = new_height;
        cs.back().height = new_rect.J;
        cs.back().S.push_back(new_rect);
        cs.back().S.back().I = new_rect.J - rect.J;
      }
      if(add_A) {
        int b = cs.size()-1;
        if(rect.J > 0) {
          cs[b].S.push_back(rect);
          cs[b].S.back().prev = 0;
          cs[b-1].S[0].next = cs[b].S.size()-1;
        }
      }
      last_height = new_height;
      add_A = 0;
      add_rectangle = 0;
    } else {
      cs.back().width = width;
      width = 0;
      cs.push_back(cross_section(new_height,i,t[i]));
      if(new_height > 0) {
        rectangle new_rect;
        new_rect.extend(t[i],graph);
        new_rect.J = new_height;
        cs.back().height = new_rect.J;
        cs.back().S.push_back(new_rect);
        cs.back().S.back().I = new_rect.J - rect.J;
      }
    }
  }
  if(cs.size() == 1) {
    cs.back().width = width;
  }
  cs.back().width += graph.node_length(t.back().node_id) - 1;
  for(int i = 0; i < cs.size(); i++) {
    tot_width += cs[i].width;
  }
}

inline bool haplo_d::has_joining_node(int index) {
  if(cs[index].S[0].prev == -1) {
    return true;
  } else {
    return false;
  }
}

rectangle* haplo_d::joining_node(int index) {
  if(!has_joining_node(index)) {
    return nullptr;
  } else {
    return &cs[index].S[0];
  }
}

rectangle* haplo_d::last_continuing(int index) {
  if(!has_joining_node(index)) {
    return &cs[index].S[0];
  } else {
    if(cs[index].S.size() > 1) {
      return &cs[index].S[1];
    } else {
      return nullptr;
    }
  }
}

cross_section cross_section::cs_shell() {
  cross_section to_return(height, b_index, node);
  to_return.width = width;
  to_return.bridge = bridge;
  return to_return;
}

haplo_d recombine_arms(haplo_d& left, haplo_d& right, int left_cut, int right_join, xg::XG& graph) {
  haplo_d to_return;
  if(!right.has_joining_node(right_join)) {
    return to_return;
  } else {
    vector<rectangle*> boundary = right.trace_strip(right_join, 0, -1);
    rectangle rect = left.cs[left_cut].S[0];
    thread_t extension = left.cs[left_cut].bridge;
    int lastJ = rect.J;
    vector<int> boundaryDeltas;
    vector<int> boundaryJs;
    for(int i = 0; i < boundary.size(); i++) {
      to_return.cs.push_back(right.cs[right_join + i].cs_shell());
      extension.push_back(to_return.cs[i].get_node());
      int new_J = rect.get_next_J(extension,graph);
      boundaryDeltas.push_back(new_J - lastJ);
      if(new_J > 0) {
        boundaryJs.push_back(new_J);
      } else {
        break;
      }
      if(boundary[i]->J - new_J > 0) {
        rectangle joiners;
        if(i > 0) {
          joiners.prev = 0;
          to_return.cs[i-1].S[0].next = 0;
        } else {
          joiners.prev = -1;
        }
        joiners.J = boundary[i]->J;
        joiners.I = boundary[i]->J - new_J;
        to_return.cs[i].S.push_back(joiners);
      }
      if(new_J > 0) {
        rectangle continuing;
        continuing.J = new_J;
        if(i > 0) {
          continuing.prev = to_return.cs[i-1].S.size() - 1;
          to_return.cs[i-1].S.back().next = to_return.cs[i].S.size();
        }
        to_return.cs[i].S.push_back(continuing);
      }
      lastJ = new_J;
      extension = left.cs[left_cut].bridge;
    }
    for(int i = 0; i < left.cs[left_cut].S.size(); i++) {
      //build first column
    }
    vector<Edge> edges_out;
    vector<Edge> edges_in;
    for(int i = 1; i < boundaryJs.size(); i++) {
      edges_out = to_return.cs[i-1].get_last_node().is_reverse ? graph.edges_on_start(to_return.cs[i-1].get_last_node().node_id) : graph.edges_on_end(to_return.cs[i-1].get_last_node().node_id);
      edges_in = to_return.cs[i].get_node().is_reverse ? graph.edges_on_end(to_return.cs[i].get_node().node_id) : graph.edges_on_start(to_return.cs[i].get_node().node_id);
      to_return.binaryI(graph, extension, i, to_return.cs[i].S.back().prev, to_return.cs[i-1].S.size(), boundaryDeltas[i], 0, boundaryJs[i], 0, 0, edges_in, edges_out);
      for(int a = 0; a < to_return.cs[i].S.size() - 1; a++) {
        to_return.cs[i].S[a].I = to_return.cs[i].S[a].J - to_return.cs[i].S[a+1].J;
      }
      to_return.cs[i].S.back().I = to_return.cs[i].S.back().J;
    }
  }
}

int find_node(thread_t& t, xg::XG::ThreadMapping node) {
  int hint = (t.size() - 1)/2;
  return find_node(t, node, hint);
}

int find_node(haplo_d& h, xg::XG::ThreadMapping node, int hint) {
  thread_t t;
  for(int i = 0; i < h.cs.size(); i++) {
    t.push_back(h.cs[i].get_node());
  }
  return find_node(t, node, hint);
}

int find_node(thread_t& t, xg::XG::ThreadMapping node, int hint) {
  if(hint > t.size()) {
    hint = (t.size() - 1)/2;
  }
  if(t[hint].node_id == node.node_id && t[hint].is_reverse == node.is_reverse) {
    return hint;
  } else {
    int above = t.size() - hint;
    int bound = max(above,hint);
    int search_up = hint;
    int search_down = hint;
    for(int i = 1; i < bound; i++) {
      if(search_up < t.size() - 1) {
        search_up++;
        if(t[search_up].node_id == node.node_id && t[search_up].is_reverse == node.is_reverse) {
          return search_up;
        }
      }
      if(search_down > 0) {
        search_down--;
        if(t[search_down].node_id == node.node_id && t[search_down].is_reverse == node.is_reverse) {
          return search_down;
        }
      }
    }
  }
  // wasn't found!
  return -1;
}
