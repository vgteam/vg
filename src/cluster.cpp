#include <string>
#include <algorithm>
#include <utility>
#include <cstring>

#include "cluster.hpp"
#include "algorithms/subgraph.hpp"
#include "algorithms/extract_containing_graph.hpp"
#include "utility.hpp"

//#define debug_mem_clusterer

using namespace std;
using namespace structures;

namespace vg {
    
//size_t OrientedDistanceClusterer::PRUNE_COUNTER = 0;
//size_t OrientedDistanceClusterer::CLUSTER_TOTAL = 0;
//size_t OrientedDistanceClusterer::MEM_FILTER_COUNTER = 0;
//size_t OrientedDistanceClusterer::MEM_TOTAL = 0;
//size_t OrientedDistanceClusterer::SPLIT_ATTEMPT_COUNTER = 0;
//size_t OrientedDistanceClusterer::SUCCESSFUL_SPLIT_ATTEMPT_COUNTER = 0;
//size_t OrientedDistanceClusterer::PRE_SPLIT_CLUSTER_COUNTER = 0;
//size_t OrientedDistanceClusterer::POST_SPLIT_CLUSTER_COUNTER = 0;
    
MEMChainModel::MEMChainModel(
    const vector<size_t>& aln_lengths,
    const vector<vector<MaximalExactMatch> >& matches,
    const function<int64_t(pos_t)>& approx_position,
    const function<unordered_map<path_handle_t, vector<pair<size_t, bool> > >(pos_t)>& path_position,
    const function<double(const MaximalExactMatch&, const MaximalExactMatch&)>& transition_weight,
    int band_width,
    int position_depth,
    int max_connections) {
    // store the MEMs in the model
    int frag_n = 0;
    for (auto& fragment : matches) {
        ++frag_n;
        for (auto& mem : fragment) {
            // copy the MEM for each specific hit in the base graph
            // and add it in as a vertex
            for (auto& node : mem.nodes) {
                //model.emplace_back();
                //auto m = model.back();
                auto pos = make_pos_t(node);
                MEMChainModelVertex m;
                m.mem = mem;
                m.weight = mem.length();
                m.prev = nullptr;
                m.score = 0;
                m.mem.positions = path_position(pos);
                m.mem.positions[handlegraph::as_path_handle(0)].push_back(make_pair(approx_position(pos), is_rev(pos)));
                m.mem.nodes.clear();
                m.mem.nodes.push_back(node);
                m.mem.fragment = frag_n;
                m.mem.match_count = mem.match_count;
                model.push_back(m);
            }
        }
    }
    // index the model with the positions
    for (vector<MEMChainModelVertex>::iterator v = model.begin(); v != model.end(); ++v) {
        for (auto& chr : v->mem.positions) {
            auto& chr_positions = positions[chr.first];
            for (auto& pos : chr.second) {
                chr_positions[pos.first].push_back(v);
            }
        }
    }
    // sort the vertexes at each approx position by their matches and trim
    for (auto& chr : positions) {
        for (auto& p : chr.second) {
            auto& pos = p.second;
            std::sort(pos.begin(), pos.end(), [](const vector<MEMChainModelVertex>::iterator& v1,
                                                 const vector<MEMChainModelVertex>::iterator& v2) {
                          return v1->mem.length() > v2->mem.length();
                      });
            pos.resize(min(pos.size(), (size_t)position_depth));
        }
    }
    // now build up the model using the positional bandwidth
    set<pair<vector<MEMChainModelVertex>::iterator, vector<MEMChainModelVertex>::iterator> > seen;
    for (unordered_map<path_handle_t, map<int64_t, vector<vector<MEMChainModelVertex>::iterator> > >::iterator c = positions.begin(); c != positions.end(); ++
c) {
        for (map<int64_t, vector<vector<MEMChainModelVertex>::iterator> >::iterator p = c->second.begin(); p != c->second.end(); ++p) {
            for (auto& v1 : p->second) {
                // For each vertex...
                if (redundant_vertexes.count(v1)) continue;
                // ...that isn't redundant
                auto q = p;
                while (++q != c->second.end() && abs(p->first - q->first) < band_width) {
                    for (auto& v2 : q->second) {
                        // For each other vertex...
                    
                        if (redundant_vertexes.count(v2)) continue;
                        // ...that isn't redudnant
                    
                        // if this is an allowable transition, run the weighting function on it
                        if (!seen.count(make_pair(v1, v2))
                            && v1->next_cost.size() < max_connections
                            && v2->prev_cost.size() < max_connections) {
                            // There are not too many connections yet
                            seen.insert(make_pair(v1, v2));
                            if (v1->mem.fragment < v2->mem.fragment
                                || (v1->mem.fragment == v2->mem.fragment && v1->mem.begin < v2->mem.begin)) {
                                // Transition is allowable because the first comes before the second
                            
                                double weight = transition_weight(v1->mem, v2->mem);
                                if (weight > -std::numeric_limits<double>::max()) {
                                    v1->next_cost.push_back(make_pair(&*v2, weight));
                                    v2->prev_cost.push_back(make_pair(&*v1, weight));
                                }
                            } else if (v1->mem.fragment > v2->mem.fragment
                                       || (v1->mem.fragment == v2->mem.fragment && v1->mem.begin > v2->mem.begin)) {
                                // Really we want to think about the transition going the other way
                            
                                double weight = transition_weight(v2->mem, v1->mem);
                                if (weight > -std::numeric_limits<double>::max()) {
                                    v2->next_cost.push_back(make_pair(&*v1, weight));
                                    v1->prev_cost.push_back(make_pair(&*v2, weight));
                                }
                            }
                        }
                    }
                }
            }
        }
    }
}

void MEMChainModel::score(const unordered_set<MEMChainModelVertex*>& exclude) {
    // propagate the scores in the model
    for (auto& m : model) {
        // score is equal to the max inbound + mem.weight
        if (exclude.count(&m)) continue; // skip if vertex was whole cluster
        m.score = m.weight;
        for (auto& p : m.prev_cost) {
            if (p.first == nullptr) continue; // this transition is masked out
            double proposal = m.weight + p.second + p.first->score;
            if (proposal > m.score) {
                m.prev = p.first;
                m.score = proposal;
            }
        }
    }
}

MEMChainModelVertex* MEMChainModel::max_vertex(void) {
    MEMChainModelVertex* maxv = nullptr;
    for (auto& m : model) {
        if (maxv == nullptr || m.score > maxv->score) {
            maxv = &m;
        }
    }
    return maxv;
}

void MEMChainModel::clear_scores(void) {
    for (auto& m : model) {
        m.score = 0;
        m.prev = nullptr;
    }
}

vector<vector<MaximalExactMatch> > MEMChainModel::traceback(int alt_alns, bool paired, bool debug) {
    vector<vector<MaximalExactMatch> > traces;
    traces.reserve(alt_alns); // avoid reallocs so we can refer to pointers to the traces
    unordered_set<MEMChainModelVertex*> exclude;
    for (auto& v : redundant_vertexes) exclude.insert(&*v);
    for (int i = 0; i < alt_alns; ++i) {
        // score the model, accounting for excluded traces
        clear_scores();
        score(exclude);
#ifdef debug_mapper
#pragma omp critical
        {
            if (debug) {
                cerr << "MEMChainModel::traceback " << i << endl;
                display(cerr);
            }
        }
#endif
        vector<MEMChainModelVertex*> vertex_trace;
        {
            // find the maximum score
            auto* vertex = max_vertex();
            // check if we've exhausted our MEMs
            if (vertex == nullptr || vertex->score == 0) break;
#ifdef debug_mapper
#pragma omp critical
            {
                if (debug) cerr << "maximum score " << vertex->mem.sequence() << " " << vertex << ":" << vertex->score << endl;
            }
#endif
            // make trace
            while (vertex != nullptr) {
                vertex_trace.push_back(vertex);
                if (vertex->prev != nullptr) {
                    vertex = vertex->prev;
                } else {
                    break;
                }
            }
        }
        // if we have a singular match or reads are not paired, record not to use it again
        if (paired && vertex_trace.size() == 1) {
            exclude.insert(vertex_trace.front());
        }
        // fill this out when we're paired to help mask out in-fragment transitions
        unordered_set<MEMChainModelVertex*> chain_members;
        if (paired) for (auto v : vertex_trace) chain_members.insert(v);
        traces.emplace_back();
        auto& mem_trace = traces.back();
        for (auto v = vertex_trace.rbegin(); v != vertex_trace.rend(); ++v) {
            auto& vertex = **v;
            if (!paired) exclude.insert(&vertex);
            if (v != vertex_trace.rbegin()) {
                auto y = v - 1;
                MEMChainModelVertex* prev = *y;
                // mask out used transitions
                for (auto& p : vertex.prev_cost) {
                    if (p.first == prev) {
                        p.first = nullptr;
                    } else if (paired && p.first != nullptr
                               && p.first->mem.fragment != vertex.mem.fragment
                               && chain_members.count(p.first)) {
                        p.first = nullptr;
                    }
                }
            }
            mem_trace.push_back(vertex.mem);
        }
        //display_dot(cerr, vertex_trace); // for debugging
    }
    return traces;
}

// show model
void MEMChainModel::display(ostream& out) {
    for (auto& vertex : model) {
        out << vertex.mem.sequence() << ":" << vertex.mem.fragment << " " << &vertex << ":" << vertex.score << "@";
        for (auto& node : vertex.mem.nodes) {
            id_t id = gcsa::Node::id(node);
            size_t offset = gcsa::Node::offset(node);
            bool is_rev = gcsa::Node::rc(node);
            out << id << (is_rev ? "-" : "+") << ":" << offset << " ";
        }
        out << "prev: ";
        for (auto& p : vertex.prev_cost) {
            auto& next = p.first;
            if (p.first == nullptr) continue;
            out << p.first << ":" << p.second << "@";
            for (auto& node : next->mem.nodes) {
                id_t id = gcsa::Node::id(node);
                size_t offset = gcsa::Node::offset(node);
                bool is_rev = gcsa::Node::rc(node);
                out << id << (is_rev ? "-" : "+") << ":" << offset << " ";
            }
            out << " ; ";
        }
        out << " next: ";
        for (auto& p : vertex.next_cost) {
            auto& next = p.first;
            if (p.first == nullptr) continue;
            out << p.first << ":" << p.second << "@";
            for (auto& node : next->mem.nodes) {
                id_t id = gcsa::Node::id(node);
                size_t offset = gcsa::Node::offset(node);
                bool is_rev = gcsa::Node::rc(node);
                out << id << (is_rev ? "-" : "+") << ":" << offset << " ";
            }
            out << " ; ";
        }
        out << endl;
    }
}

void MEMChainModel::display_dot(ostream& out, vector<MEMChainModelVertex*> vertex_trace) {
    map<MEMChainModelVertex*, int> vertex_ids;
    int i = 0;
    for (auto& vertex : model) {
        vertex_ids[&vertex] = ++i;
    }
    map<MEMChainModelVertex*,int> in_trace;
    i = 0;
    for (auto& v : vertex_trace) {
        in_trace[v] = ++i;
    }
    out << "digraph memchain {" << endl;
    for (auto& vertex : model) {
        out << vertex_ids[&vertex]
            << " [label=\"id:" << vertex_ids[&vertex]
            << " seq:" << vertex.mem.sequence()
            << " score:" << vertex.score
            << " pos:[";
        for (auto& p : vertex.mem.positions) {
            for (auto& o : p.second) {
                out << handlegraph::as_integer(p.first) << ":" << o.first << ":" << (o.second?"-":"+") << ",";
            }
        }
        out << "]\"";
        if (in_trace.find(&vertex) != in_trace.end()) {
            out << ",color=red";
        }
        out << ",shape=box];" << endl;
        /*
        for (auto& p : vertex.prev_cost) {
            out << vertex_ids[p.first] << " -> " << vertex_ids[&vertex] << " [label=\"" << p.second << "\"];" << endl;
        }
        */
        for (auto& p : vertex.next_cost) {
            //out << in_trace[&vertex] << " " << in_trace[p.first] << endl;
            out << vertex_ids[&vertex] << " -> " << vertex_ids[p.first] << " [label=\"" << p.second << "\"";
            if (in_trace.find(&vertex) != in_trace.end() && in_trace.find(p.first) != in_trace.end() &&
                in_trace[&vertex] - 1 == in_trace[p.first]) {
                out << ",color=red";
            }
            out << "];" << endl;
        }
    }
    out << "}" << endl;
}

ShuffledPairs::ShuffledPairs(size_t num_items) : num_items(num_items), num_pairs(num_items * num_items), larger_prime(1), primitive_root(1) {
    
    // Find a prime that is at least as large as but at most a constant factor
    // larger than the number of pairs using the pre-computed list
    // If there are 0 or 1 pairs, we let the number stay at 1 (which is not prime) because that
    // turns out to be convenient
    for (size_t i = 0; larger_prime < num_pairs; i++) {
        larger_prime = spaced_primes[i];
        primitive_root = primitive_roots_of_unity[i];
    }
}

ShuffledPairs::iterator ShuffledPairs::begin() const {
    return iterator(*this, 0);
}

ShuffledPairs::iterator ShuffledPairs::end() const {
    return iterator(*this, larger_prime - 1);
}

ShuffledPairs::iterator::iterator(const ShuffledPairs& iteratee, size_t start_at) : iteratee(iteratee), permutation_idx(start_at), permuted(1) {
    
    if (permutation_idx + 1 >= iteratee.larger_prime) {
        // Don't run the dereference if we're at the end already. Handles the empty case with nothing to pair.
        return;
    }
    
    // deterministically generate pseudo-shuffled pairs in constant time per pair, adapted from
    // https://stackoverflow.com/questions/10054732/create-a-random-permutation-of-1-n-in-constant-space
    
    // See the pair we would return
    pair<size_t, size_t> returned = *(*this);
    while (permutation_idx < iteratee.larger_prime - 1 &&
           (returned.first >= returned.second || returned.first >= iteratee.num_items)) {
        
        // Advance until it's valid or we hit the end.
        permutation_idx++;
        permuted = (permuted * iteratee.primitive_root) % iteratee.larger_prime;
        returned = *(*this);
    }
}

ShuffledPairs::iterator& ShuffledPairs::iterator::operator++() {
    // Advance the permutation index
    permutation_idx++;
    permuted = (permuted * iteratee.primitive_root) % iteratee.larger_prime;
    
    // See the pair we would return
    pair<size_t, size_t> returned = *(*this);
    while (permutation_idx < iteratee.larger_prime - 1 &&
           (returned.first >= returned.second || returned.first >= iteratee.num_items)) {
        
        // Advance until it's valid or we hit the end.
        permutation_idx++;
        permuted = (permuted * iteratee.primitive_root) % iteratee.larger_prime;
        returned = *(*this);
    }
    
    return *this;
}

pair<size_t, size_t> ShuffledPairs::iterator::operator*() const {
    return pair<size_t, size_t>(permuted / iteratee.num_items, permuted % iteratee.num_items);
}

bool ShuffledPairs::iterator::operator==(const iterator& other) const {
    return permutation_idx == other.permutation_idx;
}

bool ShuffledPairs::iterator::operator!=(const iterator& other) const {
    return !(*this == other);
}
    
int32_t MEMClusterer::estimate_edge_score(const MaximalExactMatch* mem_1, const MaximalExactMatch* mem_2,
                                          int64_t graph_dist, const GSSWAligner* aligner) const {
    
    // the length of the sequence in between the MEMs (can be negative if they overlap)
    int64_t between_length = mem_2->begin - mem_1->end;
    
    if (between_length < 0) {
        // the MEMs overlap, but this can occur in some insertions and deletions
        // because the SMEM algorithm is "greedy" in taking up as much of the read
        // as possible
        // we can check if this happened directly, but it's expensive
        // so for now we just give it the benefit of the doubt but adjust the edge
        // score so that the matches don't get double counted
        
        int64_t extra_dist = abs(graph_dist - between_length);
        
        return aligner->match * between_length
               - (extra_dist ? (extra_dist - 1) * aligner->gap_extension + aligner->gap_open : 0);
    }
    else {
        int64_t gap_length = abs(between_length - graph_dist);
        // the read length in between the MEMs is the same as the distance, suggesting a pure mismatch
        return gap_length ? -((gap_length - 1) * aligner->gap_extension + aligner->gap_open) : 0;
    }
}

void MEMClusterer::deduplicate_cluster_pairs(vector<pair<pair<size_t, size_t>, int64_t>>& cluster_pairs,
                                             int64_t optimal_separation) {
    
    // sort so that pairs with same clusters are adjacent
    sort(cluster_pairs.begin(), cluster_pairs.end());
    
#ifdef debug_mem_clusterer
    cerr << "pairs before deduplicating:" << endl;
    for (const auto& pair_record : cluster_pairs) {
        cerr << pair_record.first.first << ", " << pair_record.first.second << ": " << pair_record.second << endl;
    }
    cerr << "target separation " << optimal_separation << endl;
#endif
    
    size_t removed_so_far = 0;
    
    for (size_t i = 0; i < cluster_pairs.size();) {
        // find the range of values that have the same pair of indices
        size_t range_end = i + 1;
        while (range_end < cluster_pairs.size() ? cluster_pairs[i].first == cluster_pairs[range_end].first : false) {
            range_end++;
        }
        
        // find the pair that is closest to the middle of the target interval
        int64_t best_separation = cluster_pairs[i].second;
        size_t best_idx = i;
        for (size_t j = i + 1; j < range_end; j++) {
            if (abs(cluster_pairs[j].second - optimal_separation) < abs(best_separation - optimal_separation)) {
                best_separation = cluster_pairs[j].second;
                best_idx = j;
            }
        }
        
        // move the best pair with these indices into the part of the vector we will keep
        cluster_pairs[i - removed_so_far] = cluster_pairs[best_idx];
        
        // we remove the entire interval except for one
        removed_so_far += range_end - i - 1;
        i = range_end;
    }
    
    // trim off the end of the vector, which now contains arbitrary values
    cluster_pairs.resize(cluster_pairs.size() - removed_so_far);
    
#ifdef debug_mem_clusterer
    cerr << "pairs after deduplicating:" << endl;
    for (const auto& pair_record : cluster_pairs) {
        cerr << pair_record.first.first << ", " << pair_record.first.second << ": " << pair_record.second << endl;
    }
#endif
}
    
MEMClusterer::HitGraph::HitGraph(const vector<MaximalExactMatch>& mems, const Alignment& alignment,
                                 const GSSWAligner* aligner, size_t min_mem_length, bool track_components,
                                 const match_fanouts_t* fanouts) :
    track_components(track_components), components(0, false)
{
    // there generally will be at least as many nodes as MEMs, so we can speed up the reallocation
    nodes.reserve(mems.size());
    
    for (const MaximalExactMatch& mem : mems) {
        
        //#pragma omp atomic
        //        MEM_TOTAL += mem.nodes.size();
        
        if (mem.length() < min_mem_length) {
#ifdef debug_mem_clusterer
            cerr << "skipping short MEM " << mem << endl;
#endif
            
            //#pragma omp atomic
            //            MEM_FILTER_COUNTER += mem.nodes.size();
            continue;
        }
        
        int32_t mem_score = aligner->score_exact_match(mem.begin, mem.end,
                                                       alignment.quality().begin() + (mem.begin - alignment.sequence().begin()));
        
        // adjust the score downward for fan-out mismatches
        if (fanouts && fanouts->count(&mem)) {
            for (const auto& fanout : fanouts->at(&mem)) {
                mem_score += (aligner->score_mismatch(fanout.first, fanout.first + 1,
                                                      alignment.quality().begin() + (fanout.first - alignment.sequence().begin()))
                              - aligner->score_exact_match(fanout.first, fanout.first + 1,
                                                           alignment.quality().begin() + (fanout.first - alignment.sequence().begin())));
            }
        }
        
#ifdef debug_mem_clusterer
        cerr << "adding nodes for MEM " << mem << " with score " << mem_score << endl;
        if (fanouts && fanouts->count(&mem)) {
            cerr << "fanouts:" << endl;
            for (auto fanout : fanouts->at(&mem)) {
                cerr << "\t" << (fanout.first - mem.begin) << ": " << *fanout.first << " -> " << fanout.second << endl;
            }
        }
        cerr << "locations:" << endl;
#endif
        for (gcsa::node_type mem_hit : mem.nodes) {
            nodes.emplace_back(mem, make_pos_t(mem_hit), mem_score);
#ifdef debug_mem_clusterer
            cerr << "\t" << nodes.size() - 1 << ": " << make_pos_t(mem_hit) << endl;
#endif
        }
    }
    
    // init the component tracker
    if (track_components) {
        components = UnionFind(nodes.size(), false);
    }
}

void MEMClusterer::HitGraph::add_edge(size_t from, size_t to, int32_t weight, int64_t distance) {
    nodes[from].edges_from.emplace_back(to, weight, distance);
    nodes[to].edges_to.emplace_back(from, weight, distance);
    
    if (track_components) {
        components.union_groups(from, to);
    }
}
    
void MEMClusterer::HitGraph::connected_components(vector<vector<size_t>>& components_out) const {
    
    components_out.clear();
    vector<bool> enqueued(nodes.size());
    
    // check each node in turn to find new components
    for (size_t i = 0; i < nodes.size(); i++) {
        if (enqueued[i]) {
            // we've already found this node from some component
            continue;
        }
        
        // this node belongs to a component we haven't found yet, use DFS to find the rest
        vector<size_t> stack {i};
        enqueued[i] = true;
        components_out.emplace_back(1, i);
        
        while (!stack.empty()) {
            
            const HitNode& node = nodes[stack.back()];
            stack.pop_back();
            
            // search in both forward and backward directions
            
            for (const HitEdge& edge : node.edges_from) {
                
                if (!enqueued[edge.to_idx]) {
                    stack.push_back(edge.to_idx);
                    enqueued[edge.to_idx] = true;
                    components_out.back().push_back(edge.to_idx);
                }
            }
            
            for (const HitEdge& edge : node.edges_to) {
                
                if (!enqueued[edge.to_idx]) {
                    stack.push_back(edge.to_idx);
                    enqueued[edge.to_idx] = true;
                    components_out.back().push_back(edge.to_idx);
                }
            }
        }
    }
}
    
void MEMClusterer::HitGraph::prune_low_scoring_edges(vector<vector<size_t>>& components, size_t component_idx, double score_factor) {
    
    vector<size_t>& component = components[component_idx];
    
    // get the topological order within this component (expressed in indexes into the component vector)
    vector<size_t> component_order;
    component_topological_order(component, component_order);
    
#ifdef debug_mem_clusterer
    cerr << "doing backwards DP" << endl;
#endif
    
    vector<int32_t> backwards_dp_score(component.size());
    unordered_map<size_t, size_t> node_idx_to_component_idx;
    for (size_t i = 0; i < component.size(); i++) {
        backwards_dp_score[i] = nodes[component[i]].score;
        node_idx_to_component_idx[component[i]] = i;
    }
    
    // do dynamic programming backwards within this component
    for (int64_t i = component_order.size() - 1; i >= 0; i--) {
        size_t idx = component_order[i];
        size_t node_idx = component[idx];
        for (HitEdge& edge : nodes[node_idx].edges_to) {
            size_t local_to_idx = node_idx_to_component_idx[edge.to_idx];
            int32_t dp_score = backwards_dp_score[idx] + edge.weight;
            if (dp_score > backwards_dp_score[local_to_idx]) {
                backwards_dp_score[local_to_idx] = dp_score;
            }
        }
    }
    
#ifdef debug_mem_clusterer
    cerr << "backwards dp scores:" << endl;
    for (size_t i = 0; i < component.size(); i++) {
        cerr << "\t" << component[i] << ": " << backwards_dp_score[i] << endl;
    }
#endif
    
    // the minimum score we will require each edge to be a part of
    int32_t min_score = *max_element(backwards_dp_score.begin(), backwards_dp_score.end()) * score_factor;
    
#ifdef debug_mem_clusterer
    cerr << "looking for edges with max score less than " << min_score << endl;
#endif
    
    for (size_t i = 0; i < component.size(); i++) {
        size_t node_idx = component[i];
        HitNode& node = nodes[node_idx];
        for (size_t j = 0; j < node.edges_from.size(); ) {
            HitEdge& edge = node.edges_from[j];
            
            // don't remove edges that look nearly perfect (helps keep redundant sub-MEMs in the cluster with
            // their parent so that they can be removed later)
            if (abs((edge.distance + (node.mem->end - node.mem->begin))
                    - (nodes[edge.to_idx].mem->begin - node.mem->begin)) <= 1) {
#ifdef debug_mem_clusterer
                cerr << "preserving edge because distance looks good" << endl;
#endif
                j++;
                continue;
            }
            
            // the forward-backward score of this edge
            int32_t edge_score = node.dp_score + edge.weight + backwards_dp_score[node_idx_to_component_idx[edge.to_idx]];
            
            // is the max score across this edge too low?
            if (edge_score < min_score) {
                
#ifdef debug_mem_clusterer
                cerr << "removing edge " << node_idx << "->" << edge.to_idx << " with weight " << edge.weight << " and max score " << edge_score << endl;
#endif
                
                // remove it's reverse counterpart
                HitNode& dest_node = nodes[edge.to_idx];
                for (size_t k = 0; k < dest_node.edges_to.size(); k++) {
                    if (dest_node.edges_to[k].to_idx == node_idx) {
#ifdef debug_mem_clusterer
                        cerr << "removing bwd edge " << edge.to_idx << "->" << dest_node.edges_to[k].to_idx << " with weight " << dest_node.edges_to[k].weight << " and max score " << edge_score << endl;
#endif
                        dest_node.edges_to[k] = dest_node.edges_to.back();
                        dest_node.edges_to.pop_back();
                        break;
                    }
                }
                
                // remove the edge
                node.edges_from[j] = node.edges_from.back();
                node.edges_from.pop_back();
            }
            else {
                j++;
            }
        }
    }
    
#ifdef debug_mem_clusterer
    cerr << "reidentifying connected components" << endl;
#endif
    
    // use DFS to identify the connected components again
    vector<vector<size_t>> new_components;
    
    vector<bool> enqueued(component.size(), false);
    for (size_t i = 0; i < component.size(); i++) {
        if (enqueued[i]) {
            continue;
        }
        new_components.emplace_back();
        vector<size_t> stack(1, component[i]);
        enqueued[i] = true;
        while (!stack.empty()) {
            size_t node_idx = stack.back();
            stack.pop_back();
            
            new_components.back().push_back(node_idx);
            
            for (HitEdge& edge : nodes[node_idx].edges_from) {
                size_t local_idx = node_idx_to_component_idx[edge.to_idx];
                if (!enqueued[local_idx]) {
                    stack.push_back(edge.to_idx);
                    enqueued[local_idx] = true;
                }
            }
            
            for (HitEdge& edge : nodes[node_idx].edges_to) {
                size_t local_idx = node_idx_to_component_idx[edge.to_idx];
                if (!enqueued[local_idx]) {
                    stack.push_back(edge.to_idx);
                    enqueued[local_idx] = true;
                }
            }
        }
    }
    
    // did we break this connected component into multiple connected components?
    if (new_components.size() > 1) {
#ifdef debug_mem_clusterer
        stringstream strm;
        strm << "splitting cluster:" << endl;
        for (auto& comp : new_components) {
            for (size_t i : comp) {
                strm << "\t" << i << " " << nodes[i].mem->sequence() << " " << nodes[i].start_pos << endl;
            }
            strm << endl;
        }
        cerr << strm.str();
#endif
        // the the original component
        components[component_idx] = move(new_components[0]);
        // add the remaining to the end
        for (size_t i = 1; i < new_components.size(); i++) {
            components.emplace_back(move(new_components[i]));
        }
    }
}
  
size_t MEMClusterer::HitGraph::median_mem_coverage(const vector<size_t>& component, const Alignment& aln) const {
    
    // express the MEMs as intervals along the read sequence
    vector<pair<int64_t, int64_t>> mem_intervals;
    for (size_t node_idx : component) {
        mem_intervals.emplace_back(nodes[node_idx].mem->begin - aln.sequence().begin(), nodes[node_idx].mem->end - aln.sequence().begin());
    }
    
    // put the intervals in order by starting index and then descending by length
    sort(mem_intervals.begin(), mem_intervals.end(), [](const pair<int64_t, int64_t>& a, const pair<int64_t, int64_t>& b) {
        return a.first < b.first || (a.first == b.first && a.second > b.second);
    });
    
#ifdef debug_median_algorithm
    cerr << "intervals:" << endl;
    for (const auto& interval : mem_intervals) {
        cerr << "\t[" << interval.first << ", " << interval.second << ")" << endl;
    }
#endif
    
    unordered_map<int64_t, int64_t> coverage_count;
    
    // a pointer to the read index we're currently at
    int64_t at = 0;
    // to keep track of how many intervals cover the current segment
    int64_t depth = 0;
    
    // we can keep track of the SMEM we're in by checking whether we've passed its final index
    pair<int64_t, int64_t> curr_smem(0, 0);
    // and the number of hits of this SMEM we've seen
    int64_t curr_smem_hit_count = 0;
    // we will skip one copy of each sub-MEM (heurstically assuming it's redundant with the parent)
    // per copy of the SMEM
    unordered_map<pair<int64_t, int64_t>, int64_t> skipped_sub_mems;
    
    // the sort order ensures we will encounter the interval starts in order, we use a priority queue
    // to also ensure that we will encounter their ends in order
    priority_queue<int64_t, vector<int64_t>, greater<int64_t>> ends;
    
    for (size_t i = 0; i < mem_intervals.size(); i++) {
        pair<int64_t, int64_t>& interval = mem_intervals[i];
        
#ifdef debug_median_algorithm
        cerr << "iter for interval [" << interval.first << ", " << interval.second << "), starting at " << at << endl;
#endif
        
        if (interval.second > curr_smem.second) {
            // we're in a MEM that covers distinct sequence from the current SMEM, so this is
            // a new SMEM (because of sort order)
            curr_smem = interval;
            curr_smem_hit_count = 1;
#ifdef debug_median_algorithm
            cerr << "\tthis is a new SMEM" << endl;
#endif
        }
        else if (interval == curr_smem) {
            // this is another hit of the same SMEM, increase the count
            curr_smem_hit_count++;
#ifdef debug_median_algorithm
            cerr << "\tthis is a repeat of the current SMEM" << endl;
#endif
        }
        else if (skipped_sub_mems[interval] < curr_smem_hit_count) {
            // we're in a MEM that covers a strict subinterval of the current SMEM, so skip
            // one sub-MEM per hit of the SMEM on the assumption that it's redundant
            skipped_sub_mems[interval]++;
#ifdef debug_median_algorithm
            cerr << "\tthis is a sub-MEM we must skip" << endl;
#endif
            continue;
        }
        
        // add the coverage of any segments that come before the start of this interval
        while (ends.empty() ? false : ends.top() <= interval.first) {
#ifdef debug_median_algorithm
            cerr << "\ttraversing interval end at " << ends.top() << " adding " << ends.top() - at << " to depth " << depth << endl;
#endif
            coverage_count[depth] += ends.top() - at;
            at = ends.top();
            ends.pop();
            
            // an interval is leaving scope, decrement the depth
            depth--;
        }
        
        // if there's an initial interval of 0 depth, we ignore it (helps with read-end effects from sequencers)
        if (at > 0 || depth > 0) {
#ifdef debug_median_algorithm
            cerr << "\ttraversing pre-interval segment staring from " << at << " adding " << interval.first - at << " to depth " << depth << endl;
#endif
            coverage_count[depth] += interval.first - at;
        }
#ifdef debug_median_algorithm
        else {
            cerr << "\tskipping an initial segment from " << at << " to " << interval.first << " with depth " << depth << endl;
        }
#endif
        
        
        at = interval.first;
        // an interval is entering scope, increment the depth
        depth++;
        ends.push(interval.second);
    }
    
    // run through the rest of the ends
    while (!ends.empty()) {
#ifdef debug_median_algorithm
        cerr << "\ttraversing interval end at " << ends.top() << " adding " << ends.top() - at << " to depth " << depth << endl;
#endif
        coverage_count[depth] += ends.top() - at;
        at = ends.top();
        ends.pop();
        
        // an interval is leaving scope, decrement the depth
        depth--;
    }
    
    // NOTE: we used to count the final interval of depth 0 here, but now we ignore 0-depth terminal intervals
    // because it seems to help with the read-end effects of sequencers (which can lead to match dropout)
    //coverage_count[0] += aln.sequence().size() - at;
    
    // convert it into a CDF over read coverage
    vector<pair<int64_t, int64_t>> cumul_coverage_count(coverage_count.begin(), coverage_count.end());
    sort(cumul_coverage_count.begin(), cumul_coverage_count.end());
    
#ifdef debug_median_algorithm
    cerr << "\tcoverage distr is: " ;
    for (const auto& record : cumul_coverage_count) {
        cerr << record.first << ":" << record.second << "  ";
    }
    cerr << endl;
#endif
    
    int64_t cumul = 0;
    for (pair<int64_t, int64_t>& coverage_record : cumul_coverage_count) {
        coverage_record.second += cumul;
        cumul = coverage_record.second;
    }
    
    // bisect to find the median
    int64_t target = aln.sequence().size() / 2;
    if (target <= cumul_coverage_count[0].second) {
        return cumul_coverage_count[0].first;
    }
    int64_t low = 0;
    int64_t hi = cumul_coverage_count.size() - 1;
    int64_t mid;
    while (hi > low + 1) {
        mid = (hi + low) / 2;
        
        if (target <= cumul_coverage_count[mid].second) {
            hi = mid;
        }
        else {
            low = mid;
        }
    }
#ifdef debug_median_algorithm
    cerr << "\tmedian is " << cumul_coverage_count[hi].first << endl;
#endif
    return cumul_coverage_count[hi].first;
}
    
void MEMClusterer::HitGraph::perform_dp() {
    
    for (HitNode& node : nodes) {
        // as in local alignment, minimum score is the score of node itself
        node.dp_score = node.score;
    }
    
#ifdef debug_mem_clusterer
    cerr << "computing topological order for clustering DP" << endl;
#endif
    
    vector<size_t> order;
    topological_order(order);
    
    for (size_t i : order) {
        HitNode& node = nodes[i];
#ifdef debug_mem_clusterer
        cerr << "at node " << i << " with DP score " << node.dp_score << " and node score " << node.score << endl;
#endif
        // for each edge out of this node
        for (HitEdge& edge : node.edges_from) {
            
            // check if the path through the node out of this edge increase score of target node
            HitNode& target_node = nodes[edge.to_idx];
            int32_t extend_score = node.dp_score + edge.weight + target_node.score;
            if (extend_score > target_node.dp_score) {
#ifdef debug_mem_clusterer
                cerr << "extending DP to node " << edge.to_idx << " with score " << extend_score << endl;
#endif
                target_node.dp_score = extend_score;
            }
        }
    }
}

void MEMClusterer::HitGraph::topological_order(vector<size_t>& order_out) const {
    
    // initialize return value
    order_out.clear();
    order_out.resize(nodes.size());
    size_t order_idx = nodes.size() - 1;
    
    // initialize iteration structures
    vector<bool> enqueued(nodes.size(), false);
    vector<size_t> edge_index(nodes.size(), 0);
    vector<size_t> stack;
    
    // iterate through starting nodes
    for (size_t init_node_idx = 0; init_node_idx < nodes.size(); init_node_idx++) {
        if (enqueued[init_node_idx]) {
            continue;
        }
        // navigate through graph with DFS
        stack.push_back(init_node_idx);
        enqueued[init_node_idx] = true;
        while (!stack.empty()) {
            size_t node_idx = stack.back();
            size_t& edge_idx = edge_index[node_idx];
            if (edge_idx < nodes[node_idx].edges_from.size()) {
                size_t target_idx = nodes[node_idx].edges_from[edge_idx].to_idx;
                if (enqueued[target_idx]) {
                    edge_index[node_idx]++;
                }
                else {
                    stack.push_back(target_idx);
                    enqueued[target_idx] = true;
                }
            }
            else {
                // add to topological order in reverse finishing order
                stack.pop_back();
                order_out[order_idx] = node_idx;
                order_idx--;
            }
        }
    }
}

void MEMClusterer::HitGraph::component_topological_order(const vector<size_t>& component,
                                                         vector<size_t>& order_out) const {
    // initialize return value
    order_out.clear();
    order_out.resize(component.size());
    
    vector<size_t> in_degree(component.size());
    vector<size_t> stack;
    unordered_map<size_t, size_t> node_idx_to_component_idx;
    for (size_t i = 0; i < component.size(); i++) {
        in_degree[i] = nodes[component[i]].edges_to.size();
        if (in_degree[i] == 0) {
            stack.push_back(i);
        }
        node_idx_to_component_idx[component[i]] = i;
    }
    
    size_t order_idx = 0;
    while (!stack.empty()) {
        size_t i = stack.back();
        stack.pop_back();
        
        order_out[order_idx] = i;
        order_idx++;
        
        for (const HitEdge& edge : nodes[component[i]].edges_from) {
            size_t j = node_idx_to_component_idx[edge.to_idx];
            in_degree[j]--;
            if (in_degree[j] == 0) {
                stack.push_back(j);
            }
        }
    }
}

void MEMClusterer::HitGraph::identify_sources_and_sinks(vector<size_t>& sources_out,
                                                        vector<size_t>& sinks_out) const {
    
    sources_out.clear();
    sinks_out.clear();
    
    vector<bool> is_source(nodes.size(), true);
    
    for (size_t i = 0; i < nodes.size(); i++) {
        if (nodes[i].edges_from.empty()) {
            sinks_out.push_back(i);
        }
        
        for (const HitEdge& edge : nodes[i].edges_from) {
            is_source[edge.to_idx] = false;
        }
    }
    
    for (size_t i = 0; i < nodes.size(); i++) {
        if (is_source[i]) {
            sources_out.push_back(i);
        }
    }
}

vector<MEMClusterer::cluster_t> MEMClusterer::HitGraph::clusters(const Alignment& alignment,
                                                                 const GSSWAligner* aligner,
                                                                 int32_t max_qual_score,
                                                                 int32_t log_likelihood_approx_factor,
                                                                 size_t min_median_mem_coverage_for_split,
                                                                 double suboptimal_edge_pruning_factor,
                                                                 double cluster_multiplicity_diff) {
    
    vector<cluster_t> to_return;
    if (nodes.size() == 0) {
        // this should only happen if we have filtered out all MEMs, so there are none to cluster
        return to_return;
    }
    
#ifdef debug_mem_clusterer
    cerr << "performing approximate DP across MEMs" << endl;
#endif
    perform_dp();
    
#ifdef debug_mem_clusterer
    cerr << "finding top tracebacks within connected components" << endl;
#endif
    // find the weakly connected components, which should correspond to mappings
    vector<vector<size_t>> components;
    connected_components(components);
    
#ifdef debug_mem_clusterer
    cerr << "traceback returns the following components: " << endl;
    for (size_t i = 0; i < components.size(); i++)  {
        vector<size_t>& component = components[i];
        cerr << "\tcomponent " << i << ":" << endl;
        for (size_t idx : component) {
            cerr << "\t\t" << idx << " " << nodes[idx].start_pos << " ";
            for (auto iter = nodes[idx].mem->begin; iter != nodes[idx].mem->end; iter++) {
                cerr << *iter;
            }
            cerr << endl;
        }
    }
#endif
    
    if (min_median_mem_coverage_for_split) {
#ifdef debug_mem_clusterer
        cerr << "looking for high coverage clusters to split" << endl;
#endif
        size_t num_original_components = components.size();
        for (size_t i = 0; i < num_original_components; i++) {
#ifdef debug_mem_clusterer
            cerr << "component " << i << " has median coverage " << median_mem_coverage(components[i], alignment) << endl;
#endif
            size_t curr_num_components = components.size();
            if (median_mem_coverage(components[i], alignment) >= min_median_mem_coverage_for_split) {
                //#pragma omp atomic
                //                SPLIT_ATTEMPT_COUNTER++;
#ifdef debug_mem_clusterer
                cerr << "attempting to prune and split cluster" << endl;
#endif
                
                prune_low_scoring_edges(components, i, suboptimal_edge_pruning_factor);
                
                if (components.size() > curr_num_components) {
                    //#pragma omp atomic
                    //                    SUCCESSFUL_SPLIT_ATTEMPT_COUNTER++;
                }
            }
        }
#ifdef debug_mem_clusterer
        vector<vector<size_t>> current_components;
        connected_components(current_components);
        cerr << "after splitting, from " << num_original_components << " to " << current_components.size() << " connected components" << endl;
#endif
        //#pragma omp atomic
        //        PRE_SPLIT_CLUSTER_COUNTER += num_original_components;
        //#pragma omp atomic
        //        POST_SPLIT_CLUSTER_COUNTER += components.size();
    }
    
    
    // find the node with the highest DP score in each connected component
    // each record is a pair of (score lower bound, node index)
    vector<pair<int32_t, vector<size_t>>> component_traceback_ends(components.size(),
                                                                   make_pair(numeric_limits<int32_t>::min(), vector<size_t>()));
    for (size_t i = 0; i < components.size(); i++) {
        vector<size_t>& component = components[i];
        pair<int32_t, vector<size_t>>& traceback_end = component_traceback_ends[i];
        for (size_t j = 0; j < component.size(); j++) {
            int32_t dp_score = nodes[component[j]].dp_score;
            if (dp_score > traceback_end.first) {
                // this is better than all previous scores, so throw anything we have away
                traceback_end.first = dp_score;
                traceback_end.second.clear();
                traceback_end.second.push_back(component[j]);
            }
            else if (dp_score == traceback_end.first) {
                // this is equivalent to the current best, so hold onto both
                traceback_end.second.push_back(component[j]);
            }
        }
    }
    //#pragma omp atomic
    //    CLUSTER_TOTAL += component_traceback_ends.size();
    
    std::make_heap(component_traceback_ends.begin(), component_traceback_ends.end());
    
    // estimate the minimum score a cluster must obtain to even affect the mapping quality
    // TODO: this approximation could break down sometimes, need to look into it
    int32_t top_score = component_traceback_ends.front().first;
    int32_t suboptimal_score_cutoff = top_score - log_likelihood_approx_factor * aligner->mapping_quality_score_diff(max_qual_score);
    // keep track of the scores of the clusters we take off the heap
    vector<int32_t> returned_cluster_scores;
    while (!component_traceback_ends.empty()) {
        // get the next highest scoring traceback end(s)
        auto traceback_end = component_traceback_ends.front();
        
#ifdef debug_mem_clusterer
        cerr << "checking traceback of component starting at " << traceback_end.second.front() << endl;
#endif
        // if this cluster does not look like it even affect the mapping quality of the top scoring
        // cluster, don't bother forming it
        if (traceback_end.first < suboptimal_score_cutoff) {
#ifdef debug_mem_clusterer
            cerr << "skipping rest of components on account of low score of " << traceback_end.first << " compared to max score " << top_score << " and cutoff " << suboptimal_score_cutoff << endl;
#endif
            
            //#pragma omp atomic
            //            PRUNE_COUNTER += component_traceback_ends.size() + 1;
            break;
        }
        
        // we're going to add this cluster to the return vector, take it off the heap
        std::pop_heap(component_traceback_ends.begin(), component_traceback_ends.end());
        component_traceback_ends.pop_back();
        
        // get the index of the node
        vector<size_t>& trace_stack = traceback_end.second;
        
        // traceback all optimal paths in this connected component
        
        // keep track of which indexes have already been added to the stack
        unordered_set<size_t> stacked{trace_stack.begin(), trace_stack.end()};
        
        while (!trace_stack.empty()) {
            size_t trace_idx = trace_stack.back();
            trace_stack.pop_back();
#ifdef debug_mem_clusterer
            cerr << "\ttracing back from " << trace_idx << " with DP score " << nodes[trace_idx].dp_score << " and node score " << nodes[trace_idx].score << endl;
#endif
            
            int32_t target_source_score = nodes[trace_idx].dp_score - nodes[trace_idx].score;
            for (HitEdge& edge : nodes[trace_idx].edges_to) {
#ifdef debug_mem_clusterer
                cerr << "\t\ttrace from " << edge.to_idx << " would have score " << nodes[edge.to_idx].dp_score + edge.weight + nodes[trace_idx].score << endl;
#endif
                if (nodes[edge.to_idx].dp_score + edge.weight == target_source_score && !stacked.count(edge.to_idx)) {
                    trace_stack.push_back(edge.to_idx);
                    stacked.insert(edge.to_idx);
#ifdef debug_mem_clusterer
                    cerr << "\t\tidentifying this as a proper traceback that we have not yet traced" << endl;
#endif
                }
            }
        }
        
        // make a cluster
        to_return.emplace_back();
        auto& cluster = to_return.back();
        for (size_t traced_idx : stacked) {
            HitNode& node = nodes[traced_idx];
            cluster.first.emplace_back(node.mem, node.start_pos);
        }
        // it starts with multiplicity 1 be default
        cluster.second = 1.0;
        // keep track of its score for further multiplicity calculations
        returned_cluster_scores.push_back(traceback_end.first);
        
        // put the cluster in order by read position
        sort(cluster.first.begin(), cluster.first.end(), [](const hit_t& hit_1, const hit_t& hit_2) {
            return hit_1.first->begin < hit_2.first->begin ||
            (hit_1.first->begin == hit_2.first->begin && hit_1.first->end < hit_2.first->end);
        });
    }
    
    // find out how many of the remaining clusters had similar score to the final
    // ones we're returning
    int32_t tail_equiv_diff = round(aligner->mapping_quality_score_diff(cluster_multiplicity_diff));
    int32_t min_tail_score = returned_cluster_scores.back() - tail_equiv_diff;
    int64_t num_tail_cutoff = 0;
    while (!component_traceback_ends.empty() &&
           component_traceback_ends.front().first >= min_tail_score) {
        // count it and remove it
        ++num_tail_cutoff;
        std::pop_heap(component_traceback_ends.begin(), component_traceback_ends.end());
        component_traceback_ends.pop_back();
    }
    if (num_tail_cutoff > 0) {
        // find out how many of the clusters we're returning also have similar score
        int32_t max_tail_score = returned_cluster_scores.back() + tail_equiv_diff;
        int64_t max_tail_idx = to_return.size() - 1;
        while (max_tail_idx > 0 && returned_cluster_scores[max_tail_idx - 1] <= max_tail_score) {
            --max_tail_idx;
        }
        
        // assign the corresponding multiplicity to all of clusters we're returning with similar scores
        double cluster_multiplicity = (double(to_return.size() - max_tail_idx + num_tail_cutoff)
                                       / double(to_return.size() - max_tail_idx));
        for (int64_t i = max_tail_idx; i < to_return.size(); ++i) {
            to_return[i].second = cluster_multiplicity;
        }
    }
    
    return to_return;
}
    
vector<MEMClusterer::cluster_t> MEMClusterer::clusters(const Alignment& alignment,
                                                       const vector<MaximalExactMatch>& mems,
                                                       const GSSWAligner* aligner,
                                                       size_t min_mem_length ,
                                                       int32_t max_qual_score,
                                                       int32_t log_likelihood_approx_factor,
                                                       size_t min_median_mem_coverage_for_split,
                                                       double suboptimal_edge_pruning_factor,
                                                       double cluster_multiplicity_diff,
                                                       const match_fanouts_t* fanouts) {
    
    HitGraph hit_graph = make_hit_graph(alignment, mems, aligner, min_mem_length, fanouts);
    return hit_graph.clusters(alignment, aligner, max_qual_score, log_likelihood_approx_factor,
                              min_median_mem_coverage_for_split, suboptimal_edge_pruning_factor,
                              cluster_multiplicity_diff);
    
}

MEMClusterer::HitGraph NullClusterer::make_hit_graph(const Alignment& alignment, const vector<MaximalExactMatch>& mems,
                                                     const GSSWAligner* aligner, size_t min_mem_length,
                                                     const match_fanouts_t* fanouts) {
    // intialize the hit nodes, but do not add any edges, and ignore the min mem length
    return HitGraph(mems, alignment, aligner, 1, fanouts);
}

vector<pair<pair<size_t, size_t>, int64_t>>  NullClusterer::pair_clusters(const Alignment& alignment_1,
                                                                          const Alignment& alignment_2,
                                                                          const vector<cluster_t*>& left_clusters,
                                                                          const vector<cluster_t*>& right_clusters,
                                                                          const vector<pair<size_t, size_t>>& left_alt_cluster_anchors,
                                                                          const vector<pair<size_t, size_t>>& right_alt_cluster_anchors,
                                                                          int64_t optimal_separation,
                                                                          int64_t max_deviation) {
    // do not cluster pairs.
    return vector<pair<pair<size_t, size_t>, int64_t>>();
}
    
PathOrientedDistanceMeasurer::PathOrientedDistanceMeasurer(const PathPositionHandleGraph* graph,
                                                           const PathComponentIndex* path_component_index) :
    graph(graph), path_component_index(path_component_index) {
    
}
    
int64_t PathOrientedDistanceMeasurer::oriented_distance(const pos_t& pos_1, const pos_t& pos_2) {
    
    /*
     * STEP 1: TRAVERSE OUTWARD FROM BOTH POSITIONS TO FIND A SHARED PATH STRAND
     */
    
    // maps of oriented paths to (handle, oriented distance) tuples
    unordered_map<pair<path_handle_t, bool>, vector<pair<step_handle_t, int64_t>>> path_strand_dists_1;
    unordered_map<pair<path_handle_t, bool>, vector<pair<step_handle_t, int64_t>>> path_strand_dists_2;
    
    unordered_set<pair<path_handle_t, bool>> shared_path_strands;
    
    // ensure that the paths of the start nodes are added, even if their ends are too far away
    // from the positions for the search to explore
    // TODO: this leaves the ambiguity that a node might occur multiple times on the same path, in which case
    // the tie for closest traversal to the path is broken arbitrarily
    handle_t handle_1 = graph->get_handle(id(pos_1), is_rev(pos_1));
    handle_t handle_2 = graph->get_handle(id(pos_2), is_rev(pos_2));
    
    for (const step_handle_t& step : graph->steps_of_handle(handle_1)) {
        pair<path_handle_t, bool> path_occurrence(graph->get_path_handle_of_step(step),
                                                  graph->get_handle_of_step(step) != handle_1);
        path_strand_dists_1[path_occurrence].emplace_back(step, -((int64_t) offset(pos_1)));
        
#ifdef debug_algorithms
        cerr << "[PathDistance] first position " << id(pos_1) << "[" << offset(pos_1) << "]" << (is_rev(pos_1) ? "-" : "+") << " has an initial path occurrence on " << as_integer(graph->get_path_handle_of_step(step)) << (graph->get_handle_of_step(step) != handle_1 ? "-" : "+") << endl;
#endif
    }
    
    for (const step_handle_t& step : graph->steps_of_handle(handle_2)) {
        pair<path_handle_t, bool> path_occurrence(graph->get_path_handle_of_step(step),
                                                  graph->get_handle_of_step(step) != handle_2);
        path_strand_dists_2[path_occurrence].emplace_back(step, -((int64_t) offset(pos_2)));
        
#ifdef debug_algorithms
        cerr << "[PathDistance] second position " << id(pos_2) << "[" << offset(pos_2) << "]" << (is_rev(pos_2) ? "-" : "+") << " has an initial path occurrence on " << as_integer(graph->get_path_handle_of_step(step)) << (graph->get_handle_of_step(step) != handle_2 ? "-" : "+") << endl;
#endif
        
        if (path_strand_dists_1.count(path_occurrence)) {
#ifdef debug_algorithms
            cerr << "[PathDistance] this occurrence is on a shared path" << endl;
#endif
            shared_path_strands.insert(path_occurrence);
        }
        
    }
    
    
    // if we already found shared paths on the start nodes, don't search anymore
    if (shared_path_strands.empty() && max_walk > 0) {
#ifdef debug_algorithms
        cerr << "[PathDistance] no shared paths detected, beginning traversals" << endl;
#endif
        
        // priority queues over traversals
        // distance is measure at the end of the node, so it's actually the distance
        // to the next nodes we will traverse to
        // there is a separate queue for each of the positions
        RankPairingHeap<pair<handle_t, bool>, int64_t, greater<int64_t>> queue_1, queue_2;
        
        queue_1.push_or_reprioritize(make_pair(handle_1, true), offset(pos_1) - graph->get_length(handle_1));
        queue_1.push_or_reprioritize(make_pair(handle_1, false), -offset(pos_1));
        queue_2.push_or_reprioritize(make_pair(handle_2, true), offset(pos_2) - graph->get_length(handle_2));
        queue_2.push_or_reprioritize(make_pair(handle_2, false), -offset(pos_2));
        
        while (!(queue_1.empty() && queue_2.empty()) && shared_path_strands.empty()) {
            
#ifdef debug_algorithms
            cerr << "[PathDistance] choosing queue for next traversal" << endl;
#endif
            // we'll use whichever queue has the shortest traversal so far
            auto curr_queue = &queue_1;
            auto curr_path_strand_dists = &path_strand_dists_1;
            auto other_path_strand_dists = &path_strand_dists_2;
            if (queue_1.empty() ? true : (queue_2.empty() ? false : queue_1.top().second > queue_2.top().second)) {
                curr_queue = &queue_2;
                std::swap(curr_path_strand_dists, other_path_strand_dists);
            }
            
            auto trav = curr_queue->top();
            curr_queue->pop();
            
#ifdef debug_algorithms
            cerr << "[PathDistance] traversing " << graph->get_id(trav.first.first) << (graph->get_is_reverse(trav.first.first) ? "-" : "+") << " in " << (trav.first.second ? "leftward" : "rightward") << " direction at distance " << trav.second << endl;
#endif
            
            // don't look any further if the next closest traversal is beyond the maximum distance
            if (trav.second > (int64_t) max_walk) {
                break;
            }
            
            int64_t dist = trav.second + graph->get_length(trav.first.first);
            
            if (!(curr_queue == &queue_1 && trav.first.first == handle_1)
                && !(curr_queue == &queue_2 && trav.first.first == handle_2)) {
                // this is not one of the start positions, so it might have new paths on it
                for (const step_handle_t& step : graph->steps_of_handle(trav.first.first)) {
        
                    pair<path_handle_t, bool> path_occurrence(graph->get_path_handle_of_step(step),
                                                              graph->get_handle_of_step(step) != trav.first.first);
                    
#ifdef debug_algorithms
                    cerr << "\ttrav is on path " << as_integer(path_occurrence.first) << " in " << (path_occurrence.second ? "reverse" : "forward") << " orientation" << endl;
#endif
                    
                    if (!curr_path_strand_dists->count(path_occurrence)) {
                        // record the oriented distance to the forward beginning of the node, relative to the start traversal
                        (*curr_path_strand_dists)[path_occurrence].emplace_back(step, trav.first.second ? -dist : trav.second);
                        // have we found nodes that share a path yet?
                        if (other_path_strand_dists->count(path_occurrence)) {
                            shared_path_strands.insert(path_occurrence);
                        }
                    }
                }
            }
            
            graph->follow_edges(trav.first.first, trav.first.second, [&](const handle_t& next) {
#ifdef debug_algorithms
                cerr << "\tfollowing edge to " << graph->get_id(next) << (graph->get_is_reverse(next) ? "-" : "+") << " at dist " << dist << endl;
#endif
                curr_queue->push_or_reprioritize(make_pair(next, trav.first.second), dist);
            });
        }
    }
    
#ifdef debug_algorithms
    cerr << "[PathDistance] found a shared path or exhausted search distance" << endl;
#endif
    
    /*
     * STEP 2: COMPUTE THE MINIMUM DISTANCE ALONG ANY SHARED PATH STRANDS DISCOVERED
     */
    
    // we will look for minimum absolute distance, so set it to the max to begin
    int64_t approx_dist = std::numeric_limits<int64_t>::max();
    for (const pair<path_handle_t, bool>& oriented_path : shared_path_strands) {
        
#ifdef debug_algorithms
        cerr << "[PathDistance] estimating distance with shared path " << as_integer(oriented_path.first) << (oriented_path.second ? "-" : "+") << endl;
#endif
        
        for (const pair<step_handle_t, int64_t>& node_trav_1 : path_strand_dists_1[oriented_path]) {
            for (const pair<step_handle_t, int64_t>& node_trav_2 : path_strand_dists_2[oriented_path]) {
                
                // the net distance searched between the two points to get to nodes on the path
                int64_t relative_offset = node_trav_1.second - node_trav_2.second;
                
#ifdef debug_algorithms
                cerr << "[PathDistance] search offset adds up to " << relative_offset << endl;
#endif
                
                // add in the interval along the path
                if (oriented_path.second) {
                    // the interval is on the reverse strand, so we need measure from the end of the node,
                    // which is also the start of the next node
                    relative_offset += (graph->get_position_of_step(graph->get_next_step(node_trav_1.first))
                                        - graph->get_position_of_step(graph->get_next_step(node_trav_2.first)));
                }
                else {
                    relative_offset += (graph->get_position_of_step(node_trav_2.first)
                                        - graph->get_position_of_step(node_trav_1.first));
                }
                                        
#ifdef debug_algorithms
                cerr << "[PathDistance] estimating distance on path " << as_integer(oriented_path.first) << (oriented_path.second ? "-" : "+") << " at " << relative_offset << endl;
#endif
                
                // find the minimum absolute distance, but retain signing
                if (abs(relative_offset) < abs(approx_dist)) {
                    approx_dist = relative_offset;
                }
            }
        }
    }
    
#ifdef debug_algorithms
    cerr << "[PathDistance] minimum distance is estimated at " << approx_dist << endl;
#endif
    
    return approx_dist;
}
    
vector<vector<size_t>> PathOrientedDistanceMeasurer::get_buckets(const function<pos_t(size_t)>& get_position, size_t num_items) {
#ifdef debug_mem_clusterer
    cerr << "using paths to bucket distance comparisons" << endl;
#endif
    
    // the return value
    vector<vector<size_t>> buckets;
    
    // we will associate each path strand with the index of a bucket
    unordered_map<pair<path_handle_t, bool>, size_t> bucket_of_path_strand;
    
    // we will also keep track of any hits that were not on any path
    vector<size_t> non_path_hits;
    
    for (size_t i = 0; i < num_items; i++) {
        pos_t pos = get_position(i);
#ifdef debug_mem_clusterer
        cerr << "adding position " << pos << " to memo" << endl;
#endif
        
        // iterate over the path steps that this node is on
        bool on_path = false;
        for (const step_handle_t& step : graph->steps_of_handle(graph->get_handle(id(pos)))) {
            on_path = true;
            
            // key indicating a path and a strand
            pair<path_handle_t, bool> key(graph->get_path_handle_of_step(step),
                                          graph->get_is_reverse(graph->get_handle_of_step(step)) != is_rev(pos));
            
            size_t bucket;
            if (!bucket_of_path_strand.count(key)) {
                // add a new bucket
                bucket_of_path_strand[key] = buckets.size();
                bucket = buckets.size();
                buckets.emplace_back();
            }
            else {
                // access the old bucket
                bucket = bucket_of_path_strand[key];
            }
            buckets[bucket].push_back(i);
        }
        
        if (!on_path) {
            // record that this hit was not on any paths
            non_path_hits.push_back(i);
        }
    }
    
    // check the nearest nodes to each non-path hit to see if we can use them to bucket the item
    for (size_t non_path_hit : non_path_hits) {
        
        pos_t pos = get_position(non_path_hit);
        
        handle_t handle = graph->get_handle(id(pos), is_rev(pos));
        size_t right_dist = graph->get_length(handle) - offset(pos);
        size_t trav_dist = min(offset(pos), right_dist);
        if (trav_dist <= max_walk) {
            // we want to consider neighbors out this far according to our walk parameter
            
            graph->follow_edges(handle, offset(pos) < right_dist, [&](const handle_t& neighbor) {
                // check whether this neighbor is on any paths
                for (const step_handle_t& step : graph->steps_of_handle(neighbor)) {
                    
                    // key indicating a path and a strand
                    pair<path_handle_t, bool> key(graph->get_path_handle_of_step(step),
                                                  graph->get_is_reverse(graph->get_handle_of_step(step)) != is_rev(pos));
                    
                    size_t bucket;
                    if (!bucket_of_path_strand.count(key)) {
                        // add a new bucket
                        bucket_of_path_strand[key] = buckets.size();
                        bucket = buckets.size();
                        buckets.emplace_back();
                    }
                    else {
                        // access the old bucket
                        bucket = bucket_of_path_strand[key];
                    }
                    buckets[bucket].push_back(non_path_hit);
                    // we can stop after this bucketing
                    return false;
                }
                return true;
            });
        }
    }
    
    return buckets;
}

vector<pair<size_t, size_t>> PathOrientedDistanceMeasurer::exclude_merges(vector<vector<size_t>>& current_groups,
                                                                          const function<pos_t(size_t)>& get_position){
    
    
    // the pairs that we are going to exclude
    vector<pair<size_t, size_t>> excludes;
    
    if (!path_component_index) {
#ifdef debug_mem_clusterer
        cerr << "no path component index, skipping process of excluding merges" << endl;
#endif
        return excludes;
    }
    
#ifdef debug_mem_clusterer
    cerr << "using path component index to exclude strand merges" << endl;
#endif
    
    // use the component path set index to exclude some distance measurements between groups we can tell are on separate
    // strands a priori
    // TODO: I wonder if there's a way to do this without the quadratic loop (although it's quadratic in number of connected
    // components, so probably not all that bad)

#ifdef debug_mem_clusterer
    cerr << "groups: " << endl;
    for (auto& group : current_groups) {
        for (size_t idx : group) {
            cerr << idx << " ";
        }
        cerr << endl;
    }
#endif
    
    // returns the path and a bool indicating whether the search was successful
    function<pair<path_handle_t, bool>(const vector<size_t>&)> find_path_of_group = [&](const vector<size_t>& group) {
        // try to find a member of the group that is on a path
        for (size_t i : group) {
            handle_t handle = graph->get_handle(id(get_position(i)));
            for (const step_handle_t& step : graph->steps_of_handle(handle)) {
                return make_pair(graph->get_path_handle_of_step(step), true);
            }
            
        }
        // try to find a member whose neighbor is on a path
        for (size_t i : group) {
            pos_t pos = get_position(i);
            handle_t handle = graph->get_handle(id(pos));
            size_t right_dist = graph->get_length(handle) - offset(pos);
            size_t trav_dist = min(offset(pos), right_dist);
            if (trav_dist <= max_walk) {
                path_handle_t result;
                bool not_found = graph->follow_edges(handle, offset(pos) < right_dist, [&](const handle_t& neighbor) {
                    for (const step_handle_t& step : graph->steps_of_handle(neighbor)) {
                        result = graph->get_path_handle_of_step(step);
                        return false;
                    }
                    return true;
                });
                if (!not_found) {
                    return make_pair(result, true);
                }
            }
        }
        
        // we ran through every hit and did not find a path
        return make_pair(handlegraph::as_path_handle(0), false);
    };
    
    for (size_t i = 1; i < current_groups.size(); i++) {
        pair<path_handle_t, bool> i_path = find_path_of_group(current_groups[i]);
        
        if (!i_path.second) {
            continue;
        }
        
        for (size_t j = 0; j < i; j++) {
            pair<path_handle_t, bool> j_path = find_path_of_group(current_groups[j]);
            if (!j_path.second) {
                continue;
            }
            
            // we can exclude any hits that are on separate connected components
            if (!path_component_index->paths_on_same_component(i_path.first, j_path.first)) {
                excludes.emplace_back(i, j);
            }
        }
    }
    
    return excludes;
}
    
SnarlOrientedDistanceMeasurer::SnarlOrientedDistanceMeasurer(SnarlDistanceIndex* distance_index) : distance_index(distance_index) {
    
    // nothing to do
}

int64_t SnarlOrientedDistanceMeasurer::oriented_distance(const pos_t& pos_1, const pos_t& pos_2) {
    
#ifdef debug_mem_clusterer
    cerr << "measuring distance between " << pos_1 << " and " << pos_2 << endl;
#endif
    
    size_t forward_dist = minimum_distance(*distance_index, pos_1, pos_2);
    size_t backward_dist = minimum_distance(*distance_index, pos_2, pos_1);
    
    // -1 is the sentinel returned by the distance index if the distance is not measurable
    if (forward_dist == std::numeric_limits<size_t>::max() && backward_dist == std::numeric_limits<size_t>::max()) {
        // convert to the sentinel used by this interface
        return numeric_limits<int64_t>::max();
    }
    else if (forward_dist == std::numeric_limits<size_t>::max()) {
        return -(int64_t)backward_dist;
    }
    else if (backward_dist == std::numeric_limits<size_t>::max()) {
        return forward_dist;
    }
    else {
        return forward_dist < backward_dist ? forward_dist : -(int64_t)backward_dist;
    }
    
}

vector<vector<size_t>> SnarlOrientedDistanceMeasurer::get_buckets(const function<pos_t(size_t)>& get_position, size_t num_items) {
    // we don't do bucketed distance measurements with this method, return it empty
    return vector<vector<size_t>>();
}

vector<pair<size_t, size_t>> SnarlOrientedDistanceMeasurer::exclude_merges(vector<vector<size_t>>& current_groups,
                                                                           const function<pos_t(size_t)>& get_position) {
    // we don't do merge exclusion with this method, return it empty
    return vector<pair<size_t, size_t>>();
}

   
    
    
MEMClusterer::HitGraph OrientedDistanceClusterer::make_hit_graph(const Alignment& alignment, const vector<MaximalExactMatch>& mems,
                                                                 const GSSWAligner* aligner, size_t min_mem_length,
                                                                 const match_fanouts_t* fanouts) {
    
    HitGraph hit_graph(mems, alignment, aligner, min_mem_length, false, fanouts);
    
    // Get all the distances between nodes, in a forrest of unrooted trees of
    // nodes that we know are on a consistent strand.
    unordered_map<pair<size_t, size_t>, int64_t> recorded_finite_dists = get_on_strand_distance_tree(hit_graph.nodes.size(),
                                                                                                     [&](size_t node_number) {
                                                                                                         return hit_graph.nodes[node_number].start_pos;
                                                                                                     },
                                                                                                     [&](size_t node_number) {
                                                                                                         return 0;
                                                                                                     });
    
    // Flatten the trees to maps of relative position by node ID.
    vector<unordered_map<size_t, int64_t>> strand_relative_position = flatten_distance_tree(hit_graph.nodes.size(), recorded_finite_dists);
    
#ifdef debug_mem_clusterer
    for (const auto& strand : strand_relative_position) {
        cerr << "strand reconstruction: "  << endl;
        vector<size_t> order;
        for (const auto& record : strand) {
            order.push_back(record.first);
        }
        sort(order.begin(), order.end(), [&](size_t a, size_t b) {return strand.at(a) < strand.at(b);});
        for (const auto i : order) {
            int64_t strand_pos = strand.at(i);
            cerr << "\t" << i << ":\t" << strand_pos << "\t" << hit_graph.nodes[i].mem->sequence() << endl;
        }
    }
#endif
    
    // now we use the strand clusters and the estimated distances to make the DAG for the
    // approximate MEM alignment
    
    int64_t gap_open_score = aligner->gap_open;
    int64_t gap_extension_score = aligner->gap_extension;
    
    int64_t forward_gap_length = min<int64_t>(aligner->longest_detectable_gap(alignment), max_gap) + max_expected_dist_approx_error;
    for (const unordered_map<size_t, int64_t>& relative_pos : strand_relative_position) {
        
        // sort the nodes by relative position
        vector<pair<int64_t, size_t>> sorted_pos;
        for (const pair<size_t, int64_t>& pos_record : relative_pos) {
            sorted_pos.emplace_back(pos_record.second, pos_record.first);
        }
        std::sort(sorted_pos.begin(), sorted_pos.end());
        
        // find edges within each strand cluster by first identifying the interval of MEMs that meets
        // the graph distance constrant for each MEM and then checking for read colinearity and the
        // reverse distance constraint
        int64_t last_idx = sorted_pos.size() - 1;
        int64_t low = 0, hi = 0;
        for (int64_t i = 0; i < sorted_pos.size(); i++) {
            
            int64_t strand_pos = sorted_pos[i].first;
            size_t pivot_idx = sorted_pos[i].second;
            HitNode& pivot = hit_graph.nodes[pivot_idx];
            int64_t pivot_length = pivot.mem->end - pivot.mem->begin;
            int64_t suffix_length = alignment.sequence().end() - pivot.mem->begin;
            
            // the limits of how far away we might detect edges to add to the clustering graph
            int64_t target_low_pos = strand_pos - max_expected_dist_approx_error;
            int64_t target_hi_pos = strand_pos + suffix_length + forward_gap_length;
            
            // move the lower boundary of the search interval to the lowest value inside the
            // the target interval
            if (sorted_pos[low].first > target_low_pos) {
                while (low > 0 ? sorted_pos[low - 1].first > target_low_pos : false) {
                    low--;
                }
            }
            else {
                while (low < sorted_pos.size() ? sorted_pos[low].first < target_low_pos : false) {
                    low++;
                }
            }
            
            // move the upper boundary of the search interval to the highest value inside the
            // the target interval
            if (sorted_pos[hi].first > target_hi_pos) {
                while (hi > 0 ?  sorted_pos[hi].first > target_hi_pos : false) {
                    hi--;
                }
            }
            else {
                while (hi < last_idx ? sorted_pos[hi + 1].first <= target_hi_pos : false) {
                    hi++;
                }
            }
            
#ifdef debug_mem_clusterer
            cerr << "checking for possible edges from " << sorted_pos[i].second << " to MEMs between " << sorted_pos[low].first << "(" << sorted_pos[low].second << ") and " << sorted_pos[hi].first << "(" << sorted_pos[hi].second << "), which is inside the interval (" << target_low_pos << ", " << target_hi_pos << ")" << endl;
#endif
            
            for (int64_t j = low; j <= hi; j++) {
                // don't make self edges
                if (i == j) {
                    continue;
                }
                
                int64_t next_idx = sorted_pos[j].second;
                HitNode& next = hit_graph.nodes[next_idx];
                
                // the estimated distance between the end of the pivot and the start of the next MEM in the graph
                int64_t graph_dist = sorted_pos[j].first - strand_pos - pivot_length;
                
                if (next.mem->begin >= pivot.mem->begin && next.mem->end <= pivot.mem->end
                    && abs((sorted_pos[j].first - strand_pos) - (next.mem->begin - pivot.mem->begin)) <= 1) {
                    // this looks like a redundant sub-MEM
                    
                    // we add a dummy edge, but only to connect the nodes' components and join the clusters,
                    // not to actually use in dynamic programming (given arbitrary low weight that should not
                    // cause overflow)
                    hit_graph.add_edge(pivot_idx, next_idx, numeric_limits<int32_t>::lowest() / 2, graph_dist);
                    
                    continue;
                }
                else if (next.mem->begin <= pivot.mem->begin || next.mem->end <= pivot.mem->end) {
                    // these MEMs cannot be colinear along the read
                    
                    // note: we allow one of the start/end positions to be the same here even though they can't
                    // techinically overlap because it tends to soak up redundant sub-MEMs into the same connected
                    // component so that they don't get their own cluster
                    
                    continue;
                }
                
                // add the edge in
                int32_t edge_score = estimate_edge_score(pivot.mem, next.mem, graph_dist, aligner);
                hit_graph.add_edge(pivot_idx, next_idx, edge_score, graph_dist);
                
#ifdef debug_mem_clusterer
                cerr << "adding edge to MEM " << sorted_pos[j].first << "(" << sorted_pos[j].second << ") with weight " << edge_score << endl;
#endif
            }
        }
    }
    
    return hit_graph;
}

OrientedDistanceClusterer::OrientedDistanceClusterer(OrientedDistanceMeasurer& distance_measurer,
                                                     size_t max_expected_dist_approx_error)
    : distance_measurer(distance_measurer), max_expected_dist_approx_error(max_expected_dist_approx_error) {
        
}

unordered_map<pair<size_t, size_t>, int64_t> OrientedDistanceClusterer::get_on_strand_distance_tree(size_t num_items,
                                                                                                    const function<pos_t(size_t)>& get_position,
                                                                                                    const function<int64_t(size_t)>& get_offset) {
    
    // for recording the distance of any pair that we check with a finite distance
    unordered_map<pair<size_t, size_t>, int64_t> recorded_finite_dists;
    
    // for recording the number of times elements of a strand cluster have been compared
    // and found an infinite distance
    map<pair<size_t, size_t>, size_t> num_infinite_dists;
    
    // we use a union find to keep track of which MEMs have been identified as being on the same strand
    UnionFind component_union_find(num_items);
    
    size_t num_possible_merges_remaining = (num_items * (num_items - 1)) / 2;
    
    int64_t max_failed_distance_probes = 2;
    
    // an initial pass that only looks at easily identifiable buckets
    extend_dist_tree_by_buckets(get_position, get_offset, num_items, recorded_finite_dists,
                                component_union_find, num_possible_merges_remaining);
    
    // another initial pass that tries to identify groups that cannot be merged
    exclude_dist_tree_merges(get_position, num_infinite_dists, component_union_find,
                             num_possible_merges_remaining, max_failed_distance_probes);
    
    // TODO: permutations that try to assign singletons
    
    // a second pass that measures distances between randomly selected pairs
    size_t nlogn = ceil(num_items * log(num_items));
    extend_dist_tree_by_permutations(get_position, get_offset, num_items, max_failed_distance_probes, nlogn,
                                     recorded_finite_dists, num_infinite_dists, component_union_find, num_possible_merges_remaining);
    
    return recorded_finite_dists;
}
    
void OrientedDistanceClusterer::exclude_dist_tree_merges(const function<pos_t(size_t)>& get_position,
                                                         map<pair<size_t, size_t>, size_t>& num_infinite_dists,
                                                         UnionFind& component_union_find,
                                                         size_t& num_possible_merges_remaining,
                                                         int64_t max_failed_distance_probes) {
    
    // the current set of groups after bucketed merging
    vector<vector<size_t>> current_groups = component_union_find.all_groups();

    // pairs of groups that we can easily identify as being on separate connected components
    vector<pair<size_t, size_t>> excluded_pairs =  distance_measurer.exclude_merges(current_groups, get_position);
    
    // mark these pairs as unmergeable and update the accounting accordingly
    for (const auto& excluded_pair : excluded_pairs) {
        size_t group_1 = component_union_find.find_group(current_groups[excluded_pair.first].front());
        size_t group_2 = component_union_find.find_group(current_groups[excluded_pair.second].front());
        
        num_infinite_dists[make_pair(group_1, group_2)] = max_failed_distance_probes + 1;
        num_infinite_dists[make_pair(group_2, group_1)] = max_failed_distance_probes + 1;
        
        num_possible_merges_remaining -= component_union_find.group_size(group_1) * component_union_find.group_size(group_2);
    }
}
    
void OrientedDistanceClusterer::extend_dist_tree_by_buckets(const function<pos_t(size_t)>& get_position,
                                                            const function<int64_t(size_t)>& get_offset,
                                                            size_t num_items,
                                                            unordered_map<pair<size_t, size_t>, int64_t>& recorded_finite_dists,
                                                            UnionFind& component_union_find,
                                                            size_t& num_possible_merges_remaining) {
    
    vector<vector<size_t>> buckets = distance_measurer.get_buckets(get_position, num_items);
    
    // Ensure a deterministic, system independent ordering
    for (vector<size_t>& bucket : buckets) {
        sort(bucket.begin(), bucket.end());
    }
    sort(buckets.begin(), buckets.end());
    
    // use the path strands to bucket distance measurements
    for (vector<size_t>& bucket : buckets) {

        for (size_t i = 1; i < bucket.size(); i++) {
            size_t prev = bucket[i - 1];
            size_t here = bucket[i];
            
            // have these items already been identified as on the same strand?
            if (component_union_find.find_group(prev) == component_union_find.find_group(here)) {
                continue;
            }
            
            // estimate the distance
            pos_t pos_prev = get_position(prev);
            pos_t pos_here = get_position(here);
            
#ifdef debug_mem_clusterer
            cerr << "measuring distance between " << prev << " at " << pos_prev << " and " << here << " at " << pos_here << endl;
#endif
            
            int64_t dist = distance_measurer.oriented_distance(pos_prev, pos_here);
            
            
            // did we get a successful estimation?
            if (dist == numeric_limits<int64_t>::max()) {
#ifdef debug_mem_clusterer
                cerr << "they don't appear to be have a measurable distance, skipping" << endl;
#endif
                continue;
            }
            
            // add the fixed offset from the hit position
            dist += get_offset(here) - get_offset(prev);
            
#ifdef debug_mem_clusterer
            cerr << "recording distance at " << dist << endl;
#endif
            
            // merge them into a strand cluster
            recorded_finite_dists[make_pair(prev, here)] = dist;
            num_possible_merges_remaining -= component_union_find.group_size(prev) * component_union_find.group_size(here);
            component_union_find.union_groups(prev, here);
        }
    }
}
    
void OrientedDistanceClusterer::extend_dist_tree_by_permutations(const function<pos_t(size_t)>& get_position,
                                                                 const function<int64_t(size_t)>& get_offset,
                                                                 size_t num_items,
                                                                 int64_t max_failed_distance_probes,
                                                                 size_t decrement_frequency,
                                                                 unordered_map<pair<size_t, size_t>, int64_t>& recorded_finite_dists,
                                                                 map<pair<size_t, size_t>, size_t>& num_infinite_dists,
                                                                 UnionFind& component_union_find,
                                                                 size_t& num_possible_merges_remaining) {
    
    // We want to run through all possible pairsets of node numbers in a permuted order.
    ShuffledPairs shuffled_pairs(num_items);
    auto current_pair = shuffled_pairs.begin();
    size_t pairs_checked = 0;
    
    // a simulated annealing parameter loosely inspired by the cutoff for an Erdos-Renyi random graph
    // to be connected with probability approaching 1
    size_t current_max_num_probes = max_failed_distance_probes;
    
    while (num_possible_merges_remaining > 0 && current_pair != shuffled_pairs.end() && current_max_num_probes > 0) {
        // slowly lower the number of distances we need to check before we believe that two clusters are on
        // separate strands
#ifdef debug_mem_clusterer
        cerr << "checked " << pairs_checked << " pairs with max probes " << current_max_num_probes << ", decrement frequency " << decrement_frequency << ", merges remaining " << num_possible_merges_remaining << endl;
#endif
        
        if (pairs_checked % decrement_frequency == 0 && pairs_checked != 0) {
            current_max_num_probes--;
#ifdef debug_mem_clusterer
            cerr << "reducing the max number of probes to " << current_max_num_probes << endl;
#endif
            for (const pair<pair<size_t, size_t>, size_t>& inf_dist_record : num_infinite_dists) {
                // break symmetry so we don't repeat the operation twice
                if (inf_dist_record.first.first < inf_dist_record.first.second && inf_dist_record.second == current_max_num_probes) {
                    // this merge just fell below the new maximum number of distance probes
                    size_t strand_size_1 = component_union_find.group_size(inf_dist_record.first.first);
                    size_t strand_size_2 = component_union_find.group_size(inf_dist_record.first.second);
                    num_possible_merges_remaining -= strand_size_1 * strand_size_2;
#ifdef debug_mem_clusterer
                    cerr << "after reduction, the total number of probes between strand " << inf_dist_record.first.first << " and " << inf_dist_record.first.second <<  " is above max, reducing possible merges by " << strand_size_1 * strand_size_2 << " to " << num_possible_merges_remaining << endl;
#endif
                }
            }
        }
        
        
        pair<size_t, size_t> node_pair = *current_pair;
        ++current_pair;
        
        pairs_checked++;
        
        size_t strand_1 = component_union_find.find_group(node_pair.first);
        size_t strand_2 = component_union_find.find_group(node_pair.second);
        
#ifdef debug_mem_clusterer
        cerr << "checking MEMs " << node_pair.first << " and " << node_pair.second << " in cluster " << strand_1 << " and " << strand_2 << endl;
#endif
        
        if (strand_1 == strand_2) {
            // these are already identified as on the same strand, don't need to do it again
#ifdef debug_mem_clusterer
            cerr << "already on same strand" << endl;
#endif
            continue;
        }
        
        auto num_failed_probes = num_infinite_dists.find(make_pair(strand_1, strand_2));
        if (num_failed_probes == num_infinite_dists.end() ? false : num_failed_probes->second >= current_max_num_probes) {
            // we've already checked multiple distances between these strand clusters and
            // none have returned a finite distance, so we conclude that they are in fact
            // on separate clusters and decline to check any more distances
#ifdef debug_mem_clusterer
            cerr << "already have checked distance above maximum number of probes" << endl;
#endif
            continue;
        }
        
        const pos_t& pos_1 = get_position(node_pair.first);
        const pos_t& pos_2 = get_position(node_pair.second);
        
        int64_t oriented_dist = distance_measurer.oriented_distance(pos_1, pos_2);
        
#ifdef debug_mem_clusterer
        cerr << "distance between " << pos_1 << " and " << pos_2 << " estimated at " << oriented_dist << endl;
#endif
        
        if (oriented_dist == std::numeric_limits<int64_t>::max()) {
            // distance is estimated at infinity, so these are either on different strands
            // or the path heuristic failed to find a shared path
            
            if (num_failed_probes == num_infinite_dists.end()) {
                num_failed_probes = num_infinite_dists.insert(pair<pair<size_t, size_t>, size_t>(make_pair(strand_1, strand_2), 1)).first;
                num_infinite_dists[make_pair(strand_2, strand_1)] = 1;
            }
            else {
                num_failed_probes->second++;
                num_infinite_dists[make_pair(strand_2, strand_1)]++;
            }
            
            
            // this infinite distance pushed the count over the maximum number of probes, so remove
            // these merges from the pool of potential merges remaining
            if (num_failed_probes->second >= current_max_num_probes) {
                size_t strand_size_1 = component_union_find.group_size(strand_1);
                size_t strand_size_2 = component_union_find.group_size(strand_2);
                
                num_possible_merges_remaining -= strand_size_1 * strand_size_2;
                
#ifdef debug_mem_clusterer
                cerr << "number of probes " << num_failed_probes->second << " crossed max threshold of " << current_max_num_probes << ", reducing possible merges by " << strand_size_1 * strand_size_2 << " to " << num_possible_merges_remaining << endl;
#endif
            }
        }
        else {
            // the distance is finite, so merge the strand clusters
            
            // add the fixed offset of the hit from the start position
            oriented_dist += get_offset(node_pair.second) - get_offset(node_pair.first);
            
            recorded_finite_dists[node_pair] = oriented_dist;
            
            size_t strand_size_1 = component_union_find.group_size(strand_1);
            size_t strand_size_2 = component_union_find.group_size(strand_2);
            
            component_union_find.union_groups(node_pair.first, node_pair.second);
            
            // remove these from the pool of remaining merges
            num_possible_merges_remaining -= strand_size_1 * strand_size_2;
            
            size_t strand_retaining = component_union_find.find_group(node_pair.first);
            size_t strand_removing = strand_retaining == strand_1 ? strand_2 : strand_1;
            
#ifdef debug_mem_clusterer
            cerr << "probe triggered group merge, reducing possible merges by " << strand_size_1 * strand_size_2 << " to " << num_possible_merges_remaining << " and retaining strand " << strand_retaining << endl;
#endif
            
            // get the ranges in the counter for failed distance probe records for both of the strands
            auto removing_iter = num_infinite_dists.lower_bound(make_pair(strand_removing, 0));
            auto removing_end = num_infinite_dists.upper_bound(make_pair(strand_removing, numeric_limits<size_t>::max()));
            auto retaining_iter = num_infinite_dists.lower_bound(make_pair(strand_retaining, 0));
            auto retaining_end = num_infinite_dists.upper_bound(make_pair(strand_retaining, numeric_limits<size_t>::max()));
            
            vector<pair<size_t, size_t>> unseen_comparisons;
            while (removing_iter != removing_end && retaining_iter != retaining_end) {
                if (removing_iter->first.second == retaining_iter->first.second) {
                    // both the removing and the retaining strand cluster have failed probes against this cluster so
                    // we need to combine the records
                    
                    // check if we've already marked some of these merges as off limits
                    bool retaining_already_blocked = retaining_iter->second >= current_max_num_probes;
                    bool removing_already_blocked = removing_iter->second >= current_max_num_probes;
                    
                    // add the counts together
                    retaining_iter->second += removing_iter->second;
                    num_infinite_dists[make_pair(retaining_iter->first.second, strand_retaining)] += removing_iter->second;
                    
                    // update the number of possible merges remaining
                    if (retaining_already_blocked && !removing_already_blocked) {
                        num_possible_merges_remaining -= (strand_retaining == strand_1 ? strand_size_2 : strand_size_1) * component_union_find.group_size(removing_iter->first.second);
                        
#ifdef debug_mem_clusterer
                        cerr << "after merge, the total number of probes against strand " << removing_iter->first.second << " increased to " << retaining_iter->second << ", above current max of " << current_max_num_probes << ", but the retaining strand is already blocked, reducing possible merges by " << (strand_retaining == strand_1 ? strand_size_2 : strand_size_1) * component_union_find.group_size(removing_iter->first.second) << " to " << num_possible_merges_remaining << endl;
#endif
                    }
                    else if (removing_already_blocked && !retaining_already_blocked) {
                        num_possible_merges_remaining -= (strand_retaining == strand_1 ? strand_size_1 : strand_size_2) * component_union_find.group_size(removing_iter->first.second);
                        
#ifdef debug_mem_clusterer
                        cerr << "after merge, the total number of probes against strand " << removing_iter->first.second << " increased to " << retaining_iter->second << ", above current max of " << current_max_num_probes << ", but the removing strand is already blocked, reducing possible merges by " << (strand_retaining == strand_1 ? strand_size_1 : strand_size_2) * component_union_find.group_size(removing_iter->first.second) << " to " << num_possible_merges_remaining << endl;
#endif
                    }
                    else if (!retaining_already_blocked && !removing_already_blocked && retaining_iter->second >= current_max_num_probes) {
                        num_possible_merges_remaining -= (strand_size_1 + strand_size_2) * component_union_find.group_size(removing_iter->first.second);
                        
#ifdef debug_mem_clusterer
                        cerr << "after merge, the total number of probes against strand " << removing_iter->first.second << " increased to " << retaining_iter->second << ", above current max of " << current_max_num_probes << ", reducing possible merges by " << (strand_size_1 + strand_size_2) * component_union_find.group_size(removing_iter->first.second) << " to " << num_possible_merges_remaining << endl;
#endif
                        
                    }
                    removing_iter++;
                    retaining_iter++;
                }
                else if (removing_iter->first.second < retaining_iter->first.second) {
                    // the strand being removed has probes against this strand cluster, but the strand being
                    // retained does not, mark this and save it for later so that we don't invalidate the range
                    unseen_comparisons.emplace_back(removing_iter->first.second, removing_iter->second);
                    removing_iter++;
                }
                else {
                    // the strand being retained has probes against this strand cluster, but the strand being
                    // removed does not, check if we need to add the removing strand to the remaining merges
                    // counter
                    if (retaining_iter->second >= current_max_num_probes) {
                        num_possible_merges_remaining -= (strand_retaining == strand_1 ? strand_size_2 : strand_size_1) * component_union_find.group_size(retaining_iter->first.second);
                        
#ifdef debug_mem_clusterer
                        cerr << "after merge, the total number of probes against strand " << retaining_iter->first.second << " increased to " << retaining_iter->second << ", above current max of " << current_max_num_probes << ", but the retaining strand is already blocked, reducing possible merges by " << (strand_retaining == strand_1 ? strand_size_2 : strand_size_1) * component_union_find.group_size(retaining_iter->first.second) << " to " << num_possible_merges_remaining << endl;
#endif
                    }
                    retaining_iter++;
                }
            }
            
            // finish off either range
            while (removing_iter != removing_end) {
                unseen_comparisons.emplace_back(removing_iter->first.second, removing_iter->second);
                removing_iter++;
            }
            while (retaining_iter != retaining_end) {
                if (retaining_iter->second >= current_max_num_probes) {
                    num_possible_merges_remaining -= (strand_retaining == strand_1 ? strand_size_2 : strand_size_1) * component_union_find.group_size(retaining_iter->first.second);
                    
#ifdef debug_mem_clusterer
                    cerr << "after merge, the total number of probes against strand " << retaining_iter->first.second << " increased to " << retaining_iter->second << ", above current max of " << current_max_num_probes << ", but the retaining strand is already blocked, reducing possible merges by " << (strand_retaining == strand_1 ? strand_size_2 : strand_size_1) * component_union_find.group_size(retaining_iter->first.second) << " to " << num_possible_merges_remaining << endl;
#endif
                }
                retaining_iter++;
            }
            
            
            // add the probes between the removing strands and clusters that had never been compared to the retaining strand
            for (const pair<size_t, size_t>& unseen_comparison : unseen_comparisons) {
                num_infinite_dists[make_pair(unseen_comparison.first, strand_retaining)] = unseen_comparison.second;
                num_infinite_dists[make_pair(strand_retaining, unseen_comparison.first)] = unseen_comparison.second;
                
                if (unseen_comparison.second >= current_max_num_probes) {
                    num_possible_merges_remaining -= (strand_retaining == strand_1 ? strand_size_1 : strand_size_2) * component_union_find.group_size(unseen_comparison.first);
                    
#ifdef debug_mem_clusterer
                    cerr << "after merge, the total number of probes against strand " << unseen_comparison.first << " increased to " << unseen_comparison.second << ", above current max of " << current_max_num_probes << ", but the removing strand is already blocked, reducing possible merges by " << (strand_retaining == strand_1 ? strand_size_2 : strand_size_1) * component_union_find.group_size(unseen_comparison.first) << " to " << num_possible_merges_remaining << endl;
#endif
                }
            }
            
            // find the range containing the records with the removing strand again (it may have changed since we
            // altered the map)
            removing_iter = num_infinite_dists.lower_bound(make_pair(strand_removing, 0));
            removing_end = num_infinite_dists.upper_bound(make_pair(strand_removing, numeric_limits<size_t>::max()));
            if (removing_iter != removing_end) {
                // move the end so that it is an inclusive range
                removing_end--;
                
                // erase the range
                if (removing_iter == removing_end) {
                    if (removing_iter->first.first != removing_iter->first.second) {
                        num_infinite_dists.erase(make_pair(removing_iter->first.second, removing_iter->first.first));
                    }
                    num_infinite_dists.erase(removing_iter);
                }
                else {
                    // erase the previous position on each iteration so that we don't invalidate the iterator before
                    // we use it to move to the next position
                    auto removing_iter_prev = removing_iter;
                    removing_iter++;
                    while (removing_iter != removing_end) {
                        if (removing_iter_prev->first.first != removing_iter_prev->first.second) {
                            num_infinite_dists.erase(make_pair(removing_iter_prev->first.second, removing_iter_prev->first.first));
                        }
                        num_infinite_dists.erase(removing_iter_prev);
                        removing_iter_prev = removing_iter;
                        removing_iter++;
                    }
                    if (removing_iter_prev->first.first != removing_iter_prev->first.second) {
                        num_infinite_dists.erase(make_pair(removing_iter_prev->first.second, removing_iter_prev->first.first));
                    }
                    num_infinite_dists.erase(removing_iter_prev);
                    if (removing_iter->first.first != removing_iter->first.second) {
                        num_infinite_dists.erase(make_pair(removing_iter->first.second, removing_iter->first.first));
                    }
                    num_infinite_dists.erase(removing_iter);
                }
            }
        }
    }
}

vector<unordered_map<size_t, int64_t>> OrientedDistanceClusterer::flatten_distance_tree(size_t num_items,
                                                                                        const unordered_map<pair<size_t, size_t>, int64_t>& recorded_finite_dists) {
    
#ifdef debug_mem_clusterer
    cerr << "constructing strand distance tree from " << num_items << " distances records:" << endl;
    for (const auto& record : recorded_finite_dists) {
        cerr << "\t" << record.first.first << "->" << record.first.second << ": " << record.second << endl;
    }
#endif
    
    // build the graph of relative distances in adjacency list representation
    // by construction each strand cluster will be an undirected, unrooted tree
    vector<vector<size_t>> strand_distance_tree(num_items);
    for (const auto& dist_record : recorded_finite_dists) {
        strand_distance_tree[dist_record.first.first].push_back(dist_record.first.second);
        strand_distance_tree[dist_record.first.second].push_back(dist_record.first.first);
    }
    
    // now approximate the relative positions along the strand by traversing each tree and
    // treating the distances we estimated as transitive
    vector<unordered_map<size_t, int64_t>> strand_relative_position;
    vector<bool> processed(num_items, false);
    for (size_t i = 0; i < num_items; i++) {
        if (processed[i]) {
            continue;
        }
        
#ifdef debug_mem_clusterer
        cerr << "beginning a distance tree traversal at item " << i << endl;
#endif
        strand_relative_position.emplace_back();
        unordered_map<size_t, int64_t>& relative_pos = strand_relative_position.back();
        
        // arbitrarily make this node the 0 point
        relative_pos[i] = 0;
        processed[i] = true;
        
        // traverse the strand's tree with DFS
        list<size_t> queue{i};
        while (!queue.empty()) {
            size_t curr = queue.back();
            queue.pop_back();
            
            int64_t curr_pos = relative_pos[curr];
            
            for (size_t next : strand_distance_tree[curr]) {
                if (processed[next]) {
                    continue;
                }
                
                // invert the sign of the distance if we originally measured it in the other order
                int64_t dist = recorded_finite_dists.count(make_pair(curr, next)) ?
                recorded_finite_dists.at(make_pair(curr, next)) :
                -recorded_finite_dists.at(make_pair(next, curr));
                
                // find the position relative to the previous node we just traversed
                relative_pos[next] = curr_pos + dist;
                processed[next] = true;
                
                queue.push_back(next);
            }
        }
    }
    
    return strand_relative_position;
}

vector<pair<size_t, size_t>> OrientedDistanceClusterer::compute_tail_mem_coverage(const Alignment& alignment,
                                                                                  const vector<MaximalExactMatch>& mems) {
    
    // include an index for the past-the-last position on the read
    vector<pair<size_t, size_t>> mem_tail_coverage(alignment.sequence().size() + 1);
    
    if (mems.empty()) {
        return mem_tail_coverage;
    }
    
    // convert the MEMs to the read interval they cover
    vector<pair<int64_t, int64_t>> mem_intervals;
    mem_intervals.reserve(mems.size());
    for (int64_t i = 0; i < mems.size(); i++) {
        if (!mems[i].nodes.empty()) {
            mem_intervals.emplace_back(mems[i].begin - alignment.sequence().begin(),
                                       mems[i].end - alignment.sequence().begin());
        }
    }
    
    // ensure that the intervals are sorted lexicographically
    if (!std::is_sorted(mem_intervals.begin(), mem_intervals.end())) {
        std::sort(mem_intervals.begin(), mem_intervals.end());
    }
    
    // find number of SMEM beginnings strictly to the left of each position
    
    int64_t last_mem_idx = mem_intervals.size() - 1;
    int64_t mem_idx = 0;
    size_t smem_count = 0;
    
    // iterate through any sub-MEMs contained in the SMEM that share its start position
    int64_t curr_mem_begin = mem_intervals[mem_idx].first;
    int64_t curr_mem_end = mem_intervals[mem_idx].second;
    while (mem_idx < last_mem_idx ? mem_intervals[mem_idx + 1].first == curr_mem_begin : false) {
        mem_idx++;
    }
    for (int64_t i = 0; i < mem_tail_coverage.size(); i++) {
        
        mem_tail_coverage[i].first = smem_count;
        
        // are we encountering the start of another SMEM
        if (mem_idx < mem_intervals.size() ? i == mem_intervals[mem_idx].first : false) {
            smem_count++;
            // iterate to the next MEM that contains some new sequence
            curr_mem_end = mem_intervals[mem_idx].second;
            mem_idx++;
            while (mem_idx < mems.size() ? mem_intervals[mem_idx].second <= curr_mem_end : false) {
                mem_idx++;
            }
            // iterate through any sub-MEMs contained in the SMEM that share its start position
            curr_mem_begin = mem_intervals[mem_idx].first;
            while (mem_idx < last_mem_idx ? mem_intervals[mem_idx + 1].first == curr_mem_begin : false) {
                mem_idx++;
            }
        }
    }
    
    // now use insertion sort to switch the lexicographic ordering
    for (int64_t i = 1; i < mem_intervals.size(); i++) {
        int64_t j = i;
        while (mem_intervals[j].second < mem_intervals[j - 1].second ||
               (mem_intervals[j].second == mem_intervals[j - 1].second && mem_intervals[j].first < mem_intervals[j - 1].first)) {
            std::swap(mem_intervals[j], mem_intervals[j - 1]);
            j--;
            if (j == 0) {
                break;
            }
        }
    }
    
#ifdef debug_mem_clusterer
    cerr << "reversed lexicographic ordering of intervals" << endl;
    for (auto interval : mem_intervals) {
        cerr << "\t" << interval.first << " " << interval.second << endl;
    }
#endif
    
    // find number of SMEM ends strictly to the right of each position
    
    mem_idx = last_mem_idx;
    smem_count = 0;
    
    // iterate through any sub-MEMs contained in the SMEM that share its end position
    curr_mem_begin = mem_intervals[mem_idx].first;
    curr_mem_end = mem_intervals[mem_idx].second;
    while (mem_idx > 0 ? mem_intervals[mem_idx - 1].second == curr_mem_end : false) {
        mem_idx--;
    }
    
    for (int64_t i = mem_tail_coverage.size() - 1; i >= 0; i--) {
        
        mem_tail_coverage[i].second = smem_count;
        
        if (mem_idx >= 0 ? i == mem_intervals[mem_idx].second : false) {
            smem_count++;
            // iterate to the next MEM that contains some new sequence
            curr_mem_begin = mem_intervals[mem_idx].first;
            mem_idx--;
            while (mem_idx >= 0 ? mem_intervals[mem_idx].first >= curr_mem_begin : false) {
                mem_idx--;
            }
            // iterate through any sub-MEMs contained in the SMEM that share its end position
            curr_mem_end = mem_intervals[mem_idx].second;
            while (mem_idx > 0 ? mem_intervals[mem_idx - 1].second == curr_mem_end : false) {
                mem_idx--;
            }
        }
    }
    
#ifdef debug_mem_clusterer
    cerr << "computed left MEM coverage" << endl;
    for (auto pos : mem_tail_coverage) {
        cerr << pos.first << " ";
    }
    cerr << endl;
    cerr << "computed right MEM coverage" << endl;
    for (auto pos : mem_tail_coverage) {
        cerr << pos.second << " ";
    }
    cerr << endl;
#endif
    
    return mem_tail_coverage;
}

vector<pair<pair<size_t, size_t>, int64_t>> OrientedDistanceClusterer::pair_clusters(const Alignment& alignment_1,
                                                                                     const Alignment& alignment_2,
                                                                                     const vector<cluster_t*>& left_clusters,
                                                                                     const vector<cluster_t*>& right_clusters,
                                                                                     const vector<pair<size_t, size_t>>& left_alt_cluster_anchors,
                                                                                     const vector<pair<size_t, size_t>>& right_alt_cluster_anchors,
                                                                                     int64_t optimal_separation,
                                                                                     int64_t max_deviation) {
    
#ifdef debug_mem_clusterer
    cerr << "beginning clustering of MEM cluster pairs for " << left_clusters.size() << " left clusters and " << right_clusters.size() << " right clusters" << endl;
    cerr << "looking for pairs with separation within of " << max_deviation << " from " << optimal_separation << endl;
#endif
    
    // We will fill this in with all sufficiently close pairs of clusters from different reads.
    vector<pair<pair<size_t, size_t>, int64_t>> to_return;
    
    // We think of the clusters as a single linear ordering, with our clusters coming first.
    size_t total_clusters = left_clusters.size() + right_clusters.size();
    
    // We also allow for different MEM hits to serve as the position we assign the cluster to
    size_t total_alt_anchors = left_alt_cluster_anchors.size() + right_alt_cluster_anchors.size();
    
    // The total number of positions we will allow for clustering
    size_t total_cluster_positions = total_clusters + total_alt_anchors;
    
    // Compute distance trees for sets of clusters that are distance-able on consistent strands.
    unordered_map<pair<size_t, size_t>, int64_t> distance_tree = get_on_strand_distance_tree(total_cluster_positions,
         [&](size_t cluster_num) {
             // Assumes the clusters are nonempty.
             if (cluster_num < left_clusters.size()) {
                 // Grab the pos_t for the first hit in the cluster, which is sorted to be the largest one.
                 return left_clusters[cluster_num]->first.front().second;
             }
             else if (cluster_num < total_clusters) {
                 // Grab the pos_t for the largest hit from the other cluster
                 return right_clusters[cluster_num - left_clusters.size()]->first.front().second;
             }
             else if (cluster_num < total_clusters + left_alt_cluster_anchors.size()) {
                 // Grab a lower pos_t in the list of hits according to the alt anchor
                 const pair<size_t, size_t>& alt_anchor = left_alt_cluster_anchors[cluster_num - total_clusters];
                 return left_clusters[alt_anchor.first]->first.at(alt_anchor.second).second;
             }
             else {
                 // Grab an alternate pos_t for a right cluster
                 const pair<size_t, size_t>& alt_anchor = right_alt_cluster_anchors[cluster_num - total_clusters - left_alt_cluster_anchors.size()];
                 return right_clusters[alt_anchor.first]->first.at(alt_anchor.second).second;
             }
             
         },
         [&](size_t cluster_num) {
             // Give the offset of the position we chose to either the start or end of the read
             if (cluster_num < left_clusters.size()) {
                 return alignment_1.sequence().begin() - left_clusters[cluster_num]->first.front().first->begin;
             }
             else if (cluster_num < total_clusters) {
                 return alignment_2.sequence().end() - right_clusters[cluster_num - left_clusters.size()]->first.front().first->begin;
             }
             else if (cluster_num < total_clusters + left_alt_cluster_anchors.size()) {
                 const pair<size_t, size_t>& alt_anchor = left_alt_cluster_anchors[cluster_num - total_clusters];
                 return alignment_1.sequence().begin() - left_clusters[alt_anchor.first]->first.at(alt_anchor.second).first->begin;
             }
             else {
                 const pair<size_t, size_t>& alt_anchor = right_alt_cluster_anchors[cluster_num - total_clusters - left_alt_cluster_anchors.size()];
                 return alignment_2.sequence().end() - right_clusters[alt_anchor.first]->first.at(alt_anchor.second).first->begin;
             }
         });
    
    // Flatten the distance tree to a set of linear spaces, one per tree.
    vector<unordered_map<size_t, int64_t>> linear_spaces = flatten_distance_tree(total_cluster_positions, distance_tree);
    
#ifdef debug_mem_clusterer
    for (const auto& strand : linear_spaces) {
        cerr << "strand reconstruction: "  << endl;
        for (const auto& record : strand) {
            if (record.first < left_clusters.size()) {
                cerr << "\t" << record.first << " left: " << record.second << "\t" << left_clusters[record.first]->first.front().second << endl;
            }
            else if (record.first < total_clusters) {
                cerr << "\t" << record.first - left_clusters.size() << " right: " << record.second << "\t" << right_clusters[record.first - left_clusters.size()]->first.front().second << endl;
            }
            else if (record.first < total_clusters + left_alt_cluster_anchors.size()) {
                const pair<size_t, size_t>& alt_anchor = left_alt_cluster_anchors[record.first - total_clusters];
                cerr << "\t" << alt_anchor.first << "(alt " << alt_anchor.second << ") left: " << record.second << "\t" << left_clusters[alt_anchor.first]->first.front().second << endl;
            }
            else {
                const pair<size_t, size_t>& alt_anchor = right_alt_cluster_anchors[record.first - total_clusters - left_alt_cluster_anchors.size()];
                cerr << "\t" << alt_anchor.first << "(alt " << alt_anchor.second << ") right: " << record.second << "\t" << right_clusters[alt_anchor.first]->first.front().second << endl;
            }

        }
    }
#endif
    
    // choose bounds based on whether we're measuring stranded distances
    int64_t max_inter_cluster_distance = optimal_separation + max_deviation;
    int64_t min_inter_cluster_distance = optimal_separation - max_deviation;
    
    for (const unordered_map<size_t, int64_t>& linear_space : linear_spaces) {
        // For each linear space
        
        // The linear space may run forward or reverse relative to our read.
        
        // This will hold pairs of relative position and cluster number
        vector<pair<int64_t, size_t>> sorted_pos;
        for (auto& cluster_and_pos : linear_space) {
            // Flip each pair around and put it in the list to sort.
            sorted_pos.emplace_back(cluster_and_pos.second, cluster_and_pos.first);
        }
        // Sort the list ascending by the first item (relative position)
        std::sort(sorted_pos.begin(), sorted_pos.end());
        
        // Now scan for opposing pairs within the distance limit.
        // TODO: this is going to be O(n^2) in the number of clusters in range.
        // Note: but only if there are a lot of clusters within the range, if the
        // clusters are distributed sparsely it will be approximately linear
        
        // Keep a cursor to the start of the window and the end of the window.
        // When adding each new thing to the window, eject anything too far
        // behind it, then compare it to everything that is left.
        size_t window_start = 0;
        size_t window_last = 0;
        
        for (size_t i = 0; i < sorted_pos.size(); i++) {
            // we're looking for left to right connections, so don't start from the right
            if ((sorted_pos[i].second >= left_clusters.size() && sorted_pos[i].second < total_clusters)
                || (sorted_pos[i].second >= total_clusters + left_alt_cluster_anchors.size())) {
                continue;
            }
            
            size_t left_idx = (sorted_pos[i].second < total_clusters ?
                               sorted_pos[i].second : left_alt_cluster_anchors[sorted_pos[i].second - total_clusters].first);
            
            // the interval of linearized coordinates we want to form pairs to
            int64_t coord_interval_start = sorted_pos[i].first + min_inter_cluster_distance;
            int64_t coord_interval_end = sorted_pos[i].first + max_inter_cluster_distance;
            
#ifdef debug_mem_clusterer
            if (sorted_pos[i].second < total_clusters) {
                cerr << "looking for clusters consistent with cluster that starts with " << left_clusters[sorted_pos[i].second]->first.front().second << " at relative position " << sorted_pos[i].first << " in coordinate window " << coord_interval_start << ":" << coord_interval_end << endl;
            }
            else {
                const pair<size_t, size_t>& alt_anchor = left_alt_cluster_anchors[sorted_pos[i].second - total_clusters];
                cerr << "looking for clusters consistent with (alt) cluster that starts with " << left_clusters[alt_anchor.first]->first.front().second << " at relative position " << sorted_pos[i].first << " in coordinate window " << coord_interval_start << ":" << coord_interval_end << endl;
            }
#endif
            
            // move the window bounds forward until it's inside the coordinate interval
            while (window_start < sorted_pos.size() ? sorted_pos[window_start].first < coord_interval_start : false) {
                window_start++;
#ifdef debug_mem_clusterer
                if (window_start == sorted_pos.size()) {
                    cerr << "window is beyond the end of the clusters" << endl;
                }
                else {
                    cerr << "moving window start to relative position " << sorted_pos[window_start].first << endl;
                }
#endif
            }
            while (window_last + 1 < sorted_pos.size() ? sorted_pos[window_last + 1].first < coord_interval_end : false) {
                window_last++;
#ifdef debug_mem_clusterer
                cerr << "moving window end to relative position " << sorted_pos[window_last - 1].first << endl;
#endif
            }
            
            // add each pair of clusters that's from the two read ends to the return value
            for (size_t j = window_start; j <= window_last; j++) {
                if (sorted_pos[j].second < left_clusters.size()
                    || (sorted_pos[j].second >= total_clusters && sorted_pos[j].second < total_clusters + left_alt_cluster_anchors.size())) {
#ifdef debug_mem_clusterer
                    size_t idx = sorted_pos[j].second < total_clusters ? sorted_pos[j].second : sorted_pos[j].second - total_clusters;
                    cerr << "cluster at relative position " << sorted_pos[idx].first << " is from the same end, skipping" << endl;
#endif
                    continue;
                }
                
                size_t right_idx;
                if (sorted_pos[j].second < total_clusters) {
                    right_idx = sorted_pos[j].second - left_clusters.size();
                }
                else {
                    right_idx = right_alt_cluster_anchors[sorted_pos[j].second - total_clusters - left_alt_cluster_anchors.size()].first;
                }
                
#ifdef debug_mem_clusterer
                cerr << "adding pair (" << left_idx << ", " << right_idx << ") with cluster relative position " << sorted_pos[j].first << " starting with " << right_clusters[right_idx]->first.front().second << endl;
#endif
                
                to_return.emplace_back(make_pair(left_idx, right_idx),
                                       sorted_pos[j].first - sorted_pos[i].first);
            }
        }
    }
    
    if (!left_alt_cluster_anchors.empty() || !right_alt_cluster_anchors.empty()) {
        // get rid of extra copies of pairs due to alternate anchor positions
        deduplicate_cluster_pairs(to_return, optimal_separation);
    }
    
    return to_return;
}

SnarlMinDistance::SnarlMinDistance(SnarlDistanceIndex& distance_index) : distance_index(distance_index) {
    // nothing else to do
}

int64_t SnarlMinDistance::operator()(const pos_t& pos_1, const pos_t& pos_2) {
    size_t distance = minimum_distance(distance_index, pos_1, pos_2);
    return distance == std::numeric_limits<size_t>::max() ? -1 : (int64_t)distance;
}

TipAnchoredMaxDistance::TipAnchoredMaxDistance(SnarlDistanceIndex& distance_index) : distance_index(distance_index) {
    // nothing else to do
}

int64_t TipAnchoredMaxDistance::operator()(const pos_t& pos_1, const pos_t& pos_2) {
    return maximum_distance(distance_index, pos_1, pos_2);
}

TargetValueSearch::TargetValueSearch(const HandleGraph& handle_graph,
                                     DistanceHeuristic* upper_bound_heuristic,
                                     DistanceHeuristic* lower_bound_heuristic) :
    handle_graph(handle_graph), upper_bound_heuristic(upper_bound_heuristic), lower_bound_heuristic(lower_bound_heuristic) {
    // nothing else to do
}

bool TargetValueSearch::tv_path_exists(const pos_t& pos_1, const pos_t& pos_2, int64_t target_value, int64_t tolerance) {
    return !tv_path(pos_1, pos_2, target_value, tolerance).empty();
}
    
vector<handle_t> TargetValueSearch::tv_path(const pos_t& pos_1, const pos_t& pos_2, int64_t target_value, int64_t tolerance) {

    //TODO: Doesn't work for cyclic graphs since max dist returns cap, not infinity
    bool exact_min = true;//TODO: Put this somewhere else. True if the min heuristic is exact
    DistanceHeuristic& min_distance = *lower_bound_heuristic;
    DistanceHeuristic& max_distance = *upper_bound_heuristic;
    
    int64_t offset_1 = offset(pos_1);
    int64_t offset_2 = offset(pos_2);

    //map each node to the target values from that node that are out of the min
    //and max bounds for the node but within tolerance of the target
    //Only keep the smaller/larger value that is closest to the target - maximum shorter distance
    hash_map<pair<id_t, bool>, int64_t> node_to_target_longer;//Too long
    hash_map<pair<id_t, bool>, int64_t> node_to_target_shorter;//Too short 

    //Path that is closest to the target and difference from target
    pair<int64_t, pair<pair<id_t, bool>, int64_t>> next_best 
                               (-1, make_pair(make_pair(0, false), -1));

    //Best that is too long - use when min heuristic finds actual minimum 
    //difference between target and best dist, node, and target from that node
    //difference between actual best path and target, the node itself, and the 
    // target from that node
    pair<int64_t, pair<pair<id_t, bool>, int64_t>> best_long 
                               (-1, make_pair(make_pair(0, false), -1));


    //map each node and target for node to the node+target leading to it
    //TODO: Maybe better to map each node to the actual path and remove old ones as necessary
    hash_map<pair<pair<id_t, bool>, int64_t>, pair<pair<id_t, bool>, int64_t>>
             node_to_path; 
    
    //reachable node 
    vector<pair<pair<id_t, bool>, int64_t>> next_nodes; //node and target
    next_nodes.push_back(make_pair(make_pair(id(pos_1), is_rev(pos_1)),
                                        target_value + offset_1));

    handle_t h = handle_graph.get_handle(id(pos_1), is_rev(pos_1));




    //TODO: maybe move this somewhere else
    auto get_min_path = [&](pair<pair<id_t, bool>, int64_t> node) 
       -> vector<handle_t> {
        /* Assuming that the path from node to pos_2 is the best path,
         * find the path from pos_1 to pos_2 that passes through node with the
         * given target value
         */ 

        //Get the path from pos_1 to node
        list<handle_t> result;
        auto prev = node_to_path.find(node);
        handle_t curr_handle = handle_graph.get_handle(node.first.first, 
                                                  node.first.second);
        result.push_front(curr_handle);
        while (prev != node_to_path.end()) {
            pair<pair<id_t, bool>, int64_t> prev_n = prev->second;
            curr_handle = handle_graph.get_handle(
                                 prev_n.first.first, prev_n.first.second);
            result.push_front(curr_handle);
            prev = node_to_path.find(prev_n);
            
        }

        vector<handle_t> path (result.begin(), result.end());
        //Path contains handles from pos_1 to node 

        //Get the path from node to pos_2
        pos_t curr_pos = make_pos_t(node.first.first, node.first.second, 0);
        int64_t dist = min_distance(curr_pos, pos_2);
        while (id(curr_pos) != id(pos_2) || is_rev(curr_pos) != is_rev(pos_2)){

            handle_t handle = handle_graph.get_handle(id(curr_pos), 
                                                       is_rev(curr_pos));
            auto try_next = [&](const handle_t& h)-> bool {
                curr_pos = make_pos_t(handle_graph.get_id(h), 
                                      handle_graph.get_is_reverse(h), 0);

                int64_t node_len = handle_graph.get_length(handle); 
                if (min_distance(curr_pos, pos_2) + node_len == dist) {
                    //If this node is on a minimum path
                    dist = dist - node_len; 
                    path.push_back(h);
                    return false;
                } else {
                    return true;
                } 

            };
            
            handle_graph.follow_edges(handle, false, try_next);

        }   
        return path;
         
    };

    int64_t min = min_distance(pos_1, pos_2);
    if (min == -1 || target_value + tolerance < min) {
        // The positions are too far apart, or are unreachable
        return vector<handle_t>();
    }
    
    int64_t max = max_distance(pos_1, pos_2);
    if (target_value - tolerance > max) {
        // The positions are too close together
        return vector<handle_t>();
    }

    //////////// Phase 1 of tsv search: get target for each reachable node
    while (next_nodes.size() != 0) {
        //Traverse graph in DFS order, find the target at each node

        pair<pair<id_t, bool>, int64_t> next = next_nodes.back();
        next_nodes.pop_back();
        pair<id_t, bool> curr_node = next.first;
        int64_t curr_target = next.second;

  
        if (curr_node.first == id(pos_2) && curr_node.second == is_rev(pos_2)) {
            //If this node is the end node
            if (curr_target == offset_2) {
                //If perfect path
                list<handle_t> result;
                auto prev = node_to_path.find(make_pair(curr_node, curr_target));
                handle_t handle = handle_graph.get_handle(curr_node.first, curr_node.second);
                result.push_front(handle);
                while (prev != node_to_path.end()) {
                    pair<pair<id_t, bool>, int64_t> prev_n = prev->second;
                    handle = handle_graph.get_handle(
                                 prev_n.first.first, prev_n.first.second);
                    result.push_front(handle);
                    prev = node_to_path.find(prev_n);
            
                }
                return vector<handle_t>(result.begin(), result.end());
            } else {
                int64_t diff = abs(curr_target-offset_2 );
                if (next_best.first == -1 || diff < next_best.first) {
                    next_best.first = diff;
                    next_best.second = make_pair(curr_node, curr_target);
                } 
            }
        }

        //If this is any other node or the target was not hit
 
        handle_t curr_handle = handle_graph.get_handle(curr_node.first, curr_node.second);           
        int64_t new_target = curr_target - handle_graph.get_length(curr_handle);
 
        vector<handle_t> best_path;//Use this if the best path can be found using min distance

        auto add_next = [&](const handle_t& h)-> bool {
            //For each adjacent node, add it to next nodes if end node is 
            //reachable with target 

            id_t id = handle_graph.get_id(h);
            bool rev = handle_graph.get_is_reverse(h);
            pos_t new_pos = make_pos_t(id, rev, 0);

            int64_t min_dist = min_distance(new_pos, pos_2); 
            int64_t max_dist = max_distance(new_pos, pos_2); 
            int64_t lower_target = std::max((int64_t)0, (new_target - tolerance));
            int64_t upper_target = new_target + tolerance;  
 
            if (exact_min && min_dist != -1 && min_dist == new_target) {
                //If the minimum path is the best path
                node_to_path[make_pair(make_pair(id, rev), new_target)]=
                             make_pair(curr_node, curr_target);
                best_path = get_min_path(make_pair(make_pair(id, rev), new_target)); 
                return false;
            }

            if (min_dist != -1 && 
                    min_dist <= new_target && new_target <= max_dist) {
                //If the target is within the distance bounds

                auto prev = node_to_path.find(make_pair(make_pair(id, rev), 
                                                        new_target));
                if (prev == node_to_path.end()) {
                    //If this node hasn't been seen before 
                    node_to_path[make_pair(make_pair(id, rev), new_target)]=
                             make_pair(curr_node, curr_target);
                    next_nodes.emplace_back(make_pair(id, rev), new_target);
                } 

            } else if (min_dist != -1 && 
                       ((lower_target <= min_dist && min_dist <= upper_target) ||
                        (lower_target <= max_dist && max_dist <= upper_target))){

                //If no path will hit the target but there are paths 
                //within tolerance, then save for later
                //TODO: Could take a shortcut if we assume that the min dist is actual min dist

                auto prev_max_target = node_to_target_shorter.find(
                                                           make_pair(id, rev));
                auto prev_min_target = node_to_target_longer.find(
                                                           make_pair(id, rev));
                if (min_dist >= new_target) {
                    //All paths too long - want to minimize distance
                        
                    if (exact_min && (best_long.first == -1 || 
                                    min_dist - new_target < best_long.first)) {
                        //If the min heuristic is exact, only save one longer 
                        //path
                        //If the min from here is better than previous longer path
                        best_long.first = min_dist - new_target;
                        best_long.second = make_pair(make_pair(id, rev), 
                                                                    new_target);
                        node_to_path[make_pair(make_pair(id, rev), 
                               new_target)] = make_pair(curr_node, curr_target);

                    } else if (!exact_min && 
                             (prev_min_target == node_to_target_longer.end() ||
                                     new_target < prev_min_target->second)) {
                        //Target is better (smaller)than previous from this node
                        node_to_target_longer.erase(make_pair(id, rev));
                        node_to_target_longer.emplace(make_pair(id, rev), 
                                                                    new_target);
                        node_to_path[make_pair(make_pair(id, rev), new_target)]
                                         = make_pair(curr_node, curr_target);
                    }
                } else if (max_dist <= new_target && 
                            (prev_max_target == node_to_target_shorter.end() ||
                             new_target > prev_max_target->second)){
                    //All paths too short;
                    node_to_target_shorter.erase(make_pair(id, rev));
                    node_to_target_shorter.emplace(
                                                make_pair(id, rev), new_target);
                    node_to_path[make_pair(make_pair(id, rev), 
                         new_target)] =  make_pair(curr_node, curr_target);
                }
            }
            return true;
        };
        if (!handle_graph.follow_edges(curr_handle, false, add_next)){

            return best_path;

        }
        
    }
    return tv_phase2(pos_1, pos_2, target_value, tolerance, node_to_target_shorter, node_to_target_longer, best_long, next_best, node_to_path);
}
 
vector<handle_t> TargetValueSearch::tv_phase2(const pos_t& pos_1, const pos_t& pos_2, int64_t target_value, int64_t tolerance,
                                              hash_map<pair<id_t, bool>, int64_t>& node_to_target_shorter,
                                              hash_map<pair<id_t, bool>, int64_t>& node_to_target_longer,
                                              pair<int64_t, pair<pair<id_t, bool>, int64_t>>& best_long,
                                              pair<int64_t, pair<pair<id_t, bool>, int64_t>>& next_best,
                                              hash_map<pair<pair<id_t, bool>, int64_t>, pair<pair<id_t, bool>, int64_t>>& node_to_path) {
//TODO: Any path that has been found is probably still pretty good, could return it here 

    DistanceHeuristic& min_distance = *lower_bound_heuristic;
    DistanceHeuristic& max_distance = *upper_bound_heuristic;
    int64_t offset_1 = offset(pos_1);
    int64_t offset_2 = offset(pos_2);
    auto get_min_path = [&](pair<pair<id_t, bool>, int64_t> node) 
       -> vector<handle_t> {
        /* Assuming that the path from node to pos_2 is the best path,
         * find the path from pos_1 to pos_2 that passes through node with the
         * given target value
         */ 

        //Get the path from pos_1 to node
        list<handle_t> result;
        auto prev = node_to_path.find(node);
        handle_t curr_handle = handle_graph.get_handle(node.first.first, 
                                                  node.first.second);
        result.push_front(curr_handle);
        while (prev != node_to_path.end()) {
            pair<pair<id_t, bool>, int64_t> prev_n = prev->second;
            curr_handle = handle_graph.get_handle(
                                 prev_n.first.first, prev_n.first.second);
            result.push_front(curr_handle);
            prev = node_to_path.find(prev_n);
            
        }

        vector<handle_t> path (result.begin(), result.end());
        //Path contains handles from pos_1 to node 

        //Get the path from node to pos_2
        pos_t curr_pos = make_pos_t(node.first.first, node.first.second, 0);
        int64_t dist = min_distance(curr_pos, pos_2);
        while (id(curr_pos) != id(pos_2) || is_rev(curr_pos) != is_rev(pos_2)){

            handle_t handle = handle_graph.get_handle(id(curr_pos), 
                                                       is_rev(curr_pos));
            auto try_next = [&](const handle_t& h)-> bool {
                curr_pos = make_pos_t(handle_graph.get_id(h), 
                                      handle_graph.get_is_reverse(h), 0);

                int64_t node_len = handle_graph.get_length(handle); 
                if (min_distance(curr_pos, pos_2) + node_len == dist) {
                    //If this node is on a minimum path
                    dist = dist - node_len; 
                    path.push_back(h);
                    return false;
                } else {
                    return true;
                } 

            };
            
            handle_graph.follow_edges(handle, false, try_next);

        }   
        return path;
         
    };
    ///////// Phase 2
    //If there is no perfect path, look for ones still within tolerance
    auto cmp = [] (pair<pair<pair<id_t, bool>, int64_t>, int64_t> x,
                          pair<pair<pair<id_t, bool>, int64_t>, int64_t> y) {
        //Comparison function for priority queue
        return (x.second > y.second);
    };
    priority_queue<pair<pair<pair<id_t, bool>, int64_t>, int64_t>,
                  vector<pair<pair<pair<id_t, bool>, int64_t>, int64_t>>,
                  decltype(cmp)> reachable(cmp);

    //Put all nodes into 
    for (auto it : node_to_target_shorter) {
       pair<id_t, bool> node = it.first;
       int64_t target = it.second;
       pos_t pos = make_pos_t(node.first, node.second, 0);
       int64_t diff = target - max_distance(pos, pos_2) ; 
       reachable.push(make_pair(make_pair(node, target), diff));

    }
    for (auto it : node_to_target_longer) {
       pair<id_t, bool> node = it.first;
       int64_t target = it.second;
       pos_t pos = make_pos_t(node.first, node.second, 0);
       int64_t diff = min_distance(pos, pos_2) - target; 
       reachable.push(make_pair(make_pair(node, target), diff));

    } 

    while (reachable.size() != 0) {
        //Continue A* search of nodes that cannot reach pos_2 with target length

        pair<pair<pair<id_t, bool>, int64_t>,int64_t> next = reachable.top();
        reachable.pop();
        pair<id_t, bool> curr_node = next.first.first;
        int64_t curr_target = next.first.second;

        handle_t curr_handle = handle_graph.get_handle(curr_node.first, curr_node.second);
        pair<pair<id_t, bool>, int64_t> prev_node (curr_node, curr_target);
  
        if (curr_node.first == id(pos_2) && curr_node.second == is_rev(pos_2)) {
            //If this node is the end node
            int64_t diff = abs( curr_target - offset_2);
            if (next_best.first == -1 || diff < next_best.first) {
                next_best.first = diff;
                next_best.second = prev_node;
            }
        } else {

            //If this is any other node
            //TODO: Should be able to traverse the start node twice if this is a cyclic graph
            
            int64_t new_target = curr_target -
                                          handle_graph.get_length(curr_handle);
            auto add_next = [&](const handle_t& h)-> bool {
                id_t id = handle_graph.get_id(h);
                bool rev = handle_graph.get_is_reverse(h);
                pos_t new_pos = make_pos_t(id, rev, 0);


                int64_t min_dist = min_distance(new_pos, pos_2); 
                int64_t max_dist = max_distance(new_pos, pos_2); 

                if (min_dist != -1) {
                    auto prev_max_target = node_to_target_shorter.find(
                                                           make_pair(id, rev));
                    auto prev_min_target = node_to_target_longer.find(
                                                           make_pair(id, rev));
                    if (min_dist >= new_target) {
                    //If paths are too long
                        if ( min_dist <= new_target + tolerance && 
                             (prev_min_target == node_to_target_longer.end() ||
                                     new_target < prev_min_target->second)) {
                            //If this target is better than last one
                            node_to_target_longer.erase(make_pair(id, rev));
                            node_to_target_longer.emplace(make_pair(id, rev),
                                                                   new_target);
                            node_to_path[make_pair(make_pair(id, rev), 
                                                      new_target)] = prev_node;
                            int64_t diff = min_dist - new_target;
                            reachable.push(make_pair(make_pair(
                                       make_pair(id, rev), new_target), diff));
                        }
                     
                    } else if (max_dist <= new_target){
                        //All paths too short
                        if ( max_dist >= new_target - tolerance && 
                            (prev_max_target == node_to_target_shorter.end() ||
                                 new_target > prev_max_target->second)){

                            node_to_target_shorter.erase(make_pair(id, rev));
                            node_to_target_shorter.emplace(make_pair(id, rev), 
                                                                    new_target);
                            node_to_path[make_pair(make_pair(id, rev), 
                                                       new_target)] = prev_node;
                            int64_t diff = new_target - max_dist;
                            reachable.push(make_pair(make_pair(make_pair(
                                                   id, rev),new_target), diff));
                        }
                    } else {
                        //Target is within bounds again
                        //TODO: Maybe keep track of whether the path is too long or too short
                        auto prev = node_to_path.find(make_pair(make_pair(
                                                         id, rev), new_target));
                        if (prev == node_to_path.end()) {
                            //If this node hasn't been seen before 
                            node_to_path[make_pair(make_pair(id, rev), 
                                                       new_target)]= prev_node;
                            reachable.push(make_pair(make_pair(make_pair(
                                                    id, rev),  new_target), 0));
                        } 
                    }
                }

                return true;
            };
            handle_graph.follow_edges(curr_handle, false, add_next);
        }
    }

    if (best_long.first != -1 && best_long.first <= tolerance &&
                (next_best.first == -1 || 
                  best_long.first <= next_best.first)) {
        //Get path for the best that is longer than the target

        return get_min_path(best_long.second);

    } else if (next_best.first != -1 && next_best.first <= tolerance) {
        
        //Backtrack to get path
        list<handle_t> result;
        auto prev = node_to_path.find(next_best.second);
        handle_t handle = handle_graph.get_handle(next_best.second.first.first,
                                                 next_best.second.first.second);
        result.push_front(handle);
        while (prev != node_to_path.end()) {
            pair<pair<id_t, bool>, int64_t> prev_node = prev->second;
            handle = handle_graph.get_handle(prev_node.first.first, 
                                            prev_node.first.second); 
            result.push_front(handle);
            prev = node_to_path.find(prev_node);
            
        }

        return vector<handle_t>(result.begin(), result.end());
    } else {

        return vector<handle_t>();
    }
}

int64_t TargetValueSearch::tv_path_length(const pos_t& pos_1, const pos_t& pos_2, int64_t target_value, int64_t tolerance) {

    vector<handle_t> path = tv_path(pos_1, pos_2, target_value, tolerance);
    if (path.empty()) {
        return numeric_limits<int64_t>::max();
    }
    else {
        // TODO: we should move tv_path into an internal function that also returns length,
        // there shouldn't be any reason to recompute it here!
        int64_t length = offset(pos_2) - offset(pos_1);
        for (size_t i = 0, end = path.size() - 1; i < end; i++) {
            length += handle_graph.get_length(path[i]);
        }
        return length;
    }
}
    
TVSClusterer::TVSClusterer(const HandleGraph* handle_graph, SnarlDistanceIndex* distance_index) :
      tvs(*handle_graph, new TipAnchoredMaxDistance(*distance_index), new SnarlMinDistance(*distance_index))   {
    
    // nothing else to do
}
    
MEMClusterer::HitGraph TVSClusterer::make_hit_graph(const Alignment& alignment, const vector<MaximalExactMatch>& mems,
                                                    const GSSWAligner* aligner, size_t min_mem_length,
                                                    const match_fanouts_t* fanouts) {
    
    
    // intialize with nodes
    HitGraph hit_graph(mems, alignment, aligner, min_mem_length, false, fanouts);
    
    // assumes that MEMs are given in lexicographic order by read interval
    for (size_t i = 0; i < hit_graph.nodes.size(); i++) {
        HitNode& hit_node_1 = hit_graph.nodes[i];
        
        for (size_t j = i + 1; j < hit_graph.nodes.size(); j++){
            
            HitNode& hit_node_2 = hit_graph.nodes[j];
            
            if (hit_node_2.mem->begin <= hit_node_1.mem->begin
                && hit_node_2.mem->end <= hit_node_1.mem->end) {
                // this node is at the same place or earlier in the read, so they can't be colinear
                
#ifdef debug_mem_clusterer
                cerr << "nodes " << i << " (" << hit_node_1.start_pos << ") and " << j << " (" << hit_node_2.start_pos << ") are not read colinear" << endl;
#endif
                continue;
            }
            
            // how far apart do we expect them to be based on the read?
            int64_t read_separation = hit_node_2.mem->begin - hit_node_1.mem->begin;
            
            // how long of an insert/deletion could we detect based on the scoring parameters?
            size_t longest_gap = min<int64_t>(min(aligner->longest_detectable_gap(alignment, hit_node_1.mem->end),
                                                  aligner->longest_detectable_gap(alignment, hit_node_2.mem->begin)),
                                              max_gap);
            
#ifdef debug_mem_clusterer
            cerr << "estimating distance between " << i << " (pos " << hit_node_1.start_pos << ") and " << j << " (pos " << hit_node_2.start_pos << ") with target " << read_separation << " and tolerance " << longest_gap << endl;
#endif
            
            // how close can we get to the expected distance, restricting to detectable edits
            int64_t tv_len = tvs.tv_path_length(hit_node_1.start_pos, hit_node_2.start_pos, read_separation, longest_gap);
            
#ifdef debug_mem_clusterer
            cerr << "estimate distance at " << tv_len << endl;
#endif
            
            if (tv_len == read_separation
                && ((hit_node_2.mem->begin >= hit_node_1.mem->begin && hit_node_2.mem->end < hit_node_1.mem->end)
                    || (hit_node_2.mem->begin > hit_node_1.mem->begin && hit_node_2.mem->end <= hit_node_1.mem->end))) {
                // this has the appearance of being a redundant hit of a sub-MEM, which we don't want to form
                // a separate cluster
                
                // we add a dummy edge, but only to connect the nodes' components and join the clusters,
                // not to actually use in dynamic programming (given arbitrary low weight that should not
                // cause overflow)
                hit_graph.add_edge(i, j, numeric_limits<int32_t>::lowest() / 2, tv_len);
            }
            else if (tv_len != numeric_limits<int64_t>::max()
                     && hit_node_2.mem->begin >= hit_node_1.mem->begin
                     && hit_node_2.mem->end >= hit_node_1.mem->end) {
                // there's a path within in the limit, and these hits are read colinear
                
                // the distance from the end of the first hit to the beginning of the next
                int64_t graph_dist = tv_len - (hit_node_1.mem->end - hit_node_1.mem->begin);
                
                // add the corresponding edge
                hit_graph.add_edge(i, j, estimate_edge_score(hit_node_1.mem, hit_node_2.mem, graph_dist, aligner), graph_dist);
                
            }
        }
    }
    
    return hit_graph;
}
    
vector<pair<pair<size_t, size_t>, int64_t>> TVSClusterer::pair_clusters(const Alignment& alignment_1,
                                                                        const Alignment& alignment_2,
                                                                        const vector<cluster_t*>& left_clusters,
                                                                        const vector<cluster_t*>& right_clusters,
                                                                        const vector<pair<size_t, size_t>>& left_alt_cluster_anchors,
                                                                        const vector<pair<size_t, size_t>>& right_alt_cluster_anchors,
                                                                        int64_t optimal_separation,
                                                                        int64_t max_deviation) {
    
#ifdef debug_mem_clusterer
    cerr << "clustering pairs of clusters" << endl;
#endif
    
    vector<pair<pair<size_t, size_t>, int64_t>> to_return;
    
    for (size_t i = 0, i_end = left_clusters.size() + left_alt_cluster_anchors.size(); i < i_end; i++) {
        
        // choose the appropriate left cluster and assign it a position
        size_t left_clust_idx;
        hit_t left_clust_hit;
        if (i < left_clusters.size()) {
            left_clust_idx = i;
            left_clust_hit = left_clusters[i]->first.front();
        }
        else {
            auto& alt_anchor = left_alt_cluster_anchors[i - left_clusters.size()];
            left_clust_idx = alt_anchor.first;
            left_clust_hit = left_clusters[left_clust_idx]->first.at(alt_anchor.second);
        }
        
        for (size_t j = 0, j_end  = right_clusters.size() + right_alt_cluster_anchors.size(); j < j_end; j++) {
            
            // choose the appropriate right cluster and assign it a position
            size_t right_clust_idx;
            hit_t right_clust_hit;
            if (j < right_clusters.size()) {
                right_clust_idx = j;
                right_clust_hit = right_clusters[j]->first.front();
            }
            else {
                auto& alt_anchor = right_alt_cluster_anchors[j - right_clusters.size()];
                right_clust_idx = alt_anchor.first;
                right_clust_hit = right_clusters[right_clust_idx]->first.at(alt_anchor.second);
            }
            
            // adjust the target value by how far away we are from the ends of the fragment
            int64_t left_clip = left_clust_hit.first->begin - alignment_1.sequence().begin();
            int64_t right_clip = alignment_2.sequence().end() - right_clust_hit.first->begin;
            int64_t target_separation = optimal_separation - left_clip - right_clip;
            
#ifdef debug_mem_clusterer
            cerr << "measuring distance between cluster " << left_clust_idx << " (" << left_clust_hit.second << ") and " << right_clust_idx << " (" << right_clust_hit.second << ") with target of " << target_separation << " and max deviation " << max_deviation << endl;
#endif
            
            // find the closest distance to this in the path
            int64_t tv_dist = tvs.tv_path_length(left_clust_hit.second, right_clust_hit.second,
                                                 target_separation, max_deviation);
            
#ifdef debug_mem_clusterer
            cerr << "estimate distance at " << tv_dist << endl;
#endif
            
            if (tv_dist != numeric_limits<int64_t>::max()) {
                // we found a suitable path, add it to the return vector
                to_return.emplace_back(make_pair(left_clust_idx, right_clust_idx),
                                       tv_dist + left_clip + right_clip);
                
            }
        }
    }
    
    if (!left_alt_cluster_anchors.empty() || !right_alt_cluster_anchors.empty()) {
        // get rid of extra copies of pairs due to alternate anchor positions
        deduplicate_cluster_pairs(to_return, optimal_separation);
    }
    
    return to_return;
}

MinDistanceClusterer::MinDistanceClusterer(SnarlDistanceIndex* distance_index) : distance_index(distance_index) {
    // nothing to do
}
    
vector<pair<pair<size_t, size_t>, int64_t>> MinDistanceClusterer::pair_clusters(const Alignment& alignment_1,
                                                                                const Alignment& alignment_2,
                                                                                const vector<cluster_t*>& left_clusters,
                                                                                const vector<cluster_t*>& right_clusters,
                                                                                const vector<pair<size_t, size_t>>& left_alt_cluster_anchors,
                                                                                const vector<pair<size_t, size_t>>& right_alt_cluster_anchors,
                                                                                int64_t optimal_separation,
                                                                                int64_t max_deviation) {
#ifdef debug_mem_clusterer
    cerr << "clustering pairs of clusters" << endl;
#endif
    
    vector<pair<pair<size_t, size_t>, int64_t>> to_return;
    
    for (size_t i = 0, i_end = left_clusters.size() + left_alt_cluster_anchors.size(); i < i_end; i++) {
        
        // choose the appropriate left cluster and assign it a position
        size_t left_clust_idx;
        hit_t left_clust_hit;
        if (i < left_clusters.size()) {
            left_clust_idx = i;
            left_clust_hit = left_clusters[i]->first.front();
        }
        else {
            auto& alt_anchor = left_alt_cluster_anchors[i - left_clusters.size()];
            left_clust_idx = alt_anchor.first;
            left_clust_hit = left_clusters[left_clust_idx]->first.at(alt_anchor.second);
        }
        
        for (size_t j = 0, j_end  = right_clusters.size() + right_alt_cluster_anchors.size(); j < j_end; j++) {
            
            // choose the appropriate right cluster and assign it a position
            size_t right_clust_idx;
            hit_t right_clust_hit;
            if (j < right_clusters.size()) {
                right_clust_idx = j;
                right_clust_hit = right_clusters[j]->first.front();
            }
            else {
                auto& alt_anchor = right_alt_cluster_anchors[j - right_clusters.size()];
                right_clust_idx = alt_anchor.first;
                right_clust_hit = right_clusters[right_clust_idx]->first.at(alt_anchor.second);
            }
            
#ifdef debug_mem_clusterer
            cerr << "measuring distance between cluster " << left_clust_idx << " (" << left_clust_hit.second << ") and " << right_clust_idx << " (" << right_clust_hit.second << ") with target of " << optimal_separation << " and max deviation " << max_deviation << endl;
#endif
            
            // what is the minimum distance between these hits?
            int64_t min_dist = minimum_distance(*distance_index, left_clust_hit.second, right_clust_hit.second);
            if (min_dist == std::numeric_limits<size_t>::max()) {
                size_t rev_min_dist = minimum_distance(*distance_index, left_clust_hit.second, right_clust_hit.second);
                if (rev_min_dist == std::numeric_limits<size_t>::max()) {
                    // these are not reachable, don't make a pair
                    continue;
                }
                else {
                    // this is reachable by traversing backwards, give it negative distance
                    min_dist = -(int64_t)rev_min_dist;
                }
            }
            
            
#ifdef debug_mem_clusterer
            cerr << "estimate distance at " << min_dist << endl;
#endif
            
            // adjust the distance by how far away we are from the ends of the fragment
            int64_t left_clip = left_clust_hit.first->begin - alignment_1.sequence().begin();
            int64_t right_clip = alignment_2.sequence().end() - right_clust_hit.first->begin;
            int64_t adjusted_dist = min_dist + left_clip + right_clip;
            
            if (adjusted_dist >= optimal_separation - max_deviation
                && adjusted_dist <= optimal_separation + max_deviation) {
                // we found a suitable path, add it to the return vector
                to_return.emplace_back(make_pair(left_clust_idx, right_clust_idx), adjusted_dist);
                
            }
        }
    }
    
    if (!left_alt_cluster_anchors.empty() || !right_alt_cluster_anchors.empty()) {
        // get rid of extra copies of pairs due to alternate anchor positions
        deduplicate_cluster_pairs(to_return, optimal_separation);
    }
    
    return to_return;
}

MEMClusterer::HitGraph MinDistanceClusterer::make_hit_graph(const Alignment& alignment,
                                                            const vector<MaximalExactMatch>& mems,
                                                            const GSSWAligner* aligner,
                                                            size_t min_mem_length,
                                                            const match_fanouts_t* fanouts) {
    
    // intialize with nodes
    HitGraph hit_graph(mems, alignment, aligner, min_mem_length, false, fanouts);
    
    // assumes that MEMs are given in lexicographic order by read interval
    for (size_t i = 0, j_begin = 1; i < hit_graph.nodes.size(); ++i) {
        
        HitNode& hit_node_1 = hit_graph.nodes[i];
        
        // start either at the first un-equal position or at the next position
        j_begin = max(j_begin, i + 1);
        
        // skip measuring to any additional hits of the same read interval
        while (j_begin < hit_graph.nodes.size() &&
               hit_graph.nodes[j_begin].mem->begin == hit_node_1.mem->begin &&
               hit_graph.nodes[j_begin].mem->end == hit_node_1.mem->end) {
            
            // this node is at the same place in the read, so they can't be colinear
#ifdef debug_mem_clusterer
            cerr << "nodes " << i << " (" << hit_node_1.start_pos << ") and " << j_begin << " (" << hit_graph.nodes[j_begin].start_pos << ") are not read colinear" << endl;
#endif
            ++j_begin;
        }
        
        for (size_t j = j_begin; j < hit_graph.nodes.size(); ++j){
            
            HitNode& hit_node_2 = hit_graph.nodes[j];
            
            // what is the minimum distance between these hits?
            int64_t min_dist = minimum_distance(*distance_index, hit_node_1.start_pos, hit_node_2.start_pos);
            if (min_dist == std::numeric_limits<size_t>::max()) {
                // these are not reachable, don't make an edge
                continue;
            }
            
            // how far apart do we expect them to be based on the read?
            int64_t read_separation = hit_node_2.mem->begin - hit_node_1.mem->begin;
            
            // how long of an insert/deletion could we detect based on the scoring parameters?
            size_t longest_gap = min<int64_t>(min(aligner->longest_detectable_gap(alignment, hit_node_1.mem->end),
                                                  aligner->longest_detectable_gap(alignment, hit_node_2.mem->begin)),
                                              max_gap);
            
            // is it possible that an alignment containing both could be detected with local alignment?
            if (abs(read_separation - min_dist) > longest_gap) {
                continue;
            }
            
            if (min_dist == read_separation && hit_node_2.mem->begin >= hit_node_1.mem->begin && hit_node_2.mem->end <= hit_node_1.mem->end) {
                // this has the appearance of being a redundant hit of a sub-MEM, which we don't want to form
                // a separate cluster
                
                // we add a dummy edge, but only to connect the nodes' components and join the clusters,
                // not to actually use in dynamic programming (given arbitrary low weight that should not
                // cause overflow)
                hit_graph.add_edge(i, j, numeric_limits<int32_t>::lowest() / 2, min_dist);
            }
            else if (hit_node_2.mem->begin >= hit_node_1.mem->begin
                     && hit_node_2.mem->end >= hit_node_1.mem->end) {
                // there's a path within in the limit, and these hits are read colinear
                
                // the distance from the end of the first hit to the beginning of the next
                int64_t graph_dist = min_dist - (hit_node_1.mem->end - hit_node_1.mem->begin);
                
                // add the corresponding edge
                hit_graph.add_edge(i, j, estimate_edge_score(hit_node_1.mem, hit_node_2.mem, graph_dist, aligner), graph_dist);
                
            }
        }
    }
    
    return hit_graph;
}

GreedyMinDistanceClusterer::GreedyMinDistanceClusterer(SnarlDistanceIndex* distance_index) : MinDistanceClusterer(distance_index) {
    // nothing else to do
}

MEMClusterer::HitGraph GreedyMinDistanceClusterer::make_hit_graph(const Alignment& alignment, const vector<MaximalExactMatch>& mems,
                                                                  const GSSWAligner* aligner, size_t min_mem_length,
                                                                  const match_fanouts_t* fanouts) {
    
    // init the hit graph's nodes
    HitGraph hit_graph(mems, alignment, aligner, min_mem_length, false, fanouts);
    
    // we will initialize this with the next backward and forward comparisons for each hit node
    vector<pair<int64_t, int64_t>> next_comparisons;
    next_comparisons.reserve(2 * hit_graph.nodes.size());
    
    // assumes that MEMs are given in lexicographic order by read interval
    for (size_t i = 0, j = 1; i < hit_graph.nodes.size(); ++i) {
        
        // where we expect to find the next match on the read
        auto target = hit_graph.nodes[i].mem->end + expected_separation;
        
        // start either at the previous target or one past i
        j = max(j, i + 1);
        
        // move forward until we are past the target
        while (j < hit_graph.nodes.size() && hit_graph.nodes[j].mem->begin < target) {
            ++j;
        }
        
        // move backward until the next index is past the target
        while (j > i + 1 && hit_graph.nodes[j - 1].mem->begin >= target) {
            --j;
        }
        
        // the next backwards comparison
        if (j > i + 1 && target - hit_graph.nodes[j - 1].mem->begin >= min_separation) {
            next_comparisons.emplace_back(i, j - 1);
        }
        // the backwards comparison
        if (j < hit_graph.nodes.size() && target - hit_graph.nodes[j].mem->begin <= max_separation) {
            next_comparisons.emplace_back(i, j);
        }
    }
    
    // the iteration order is starting from distances near the expected distance first, but also
    // favoring forward distances by some pre-determined multiplier
    auto priority_cmp = [&](const pair<size_t, size_t>& a, const pair<size_t, size_t>& b) {
        int64_t a_dist = (hit_graph.nodes[a.second].mem->begin - hit_graph.nodes[a.first].mem->end) - expected_separation;
        int64_t b_dist = (hit_graph.nodes[b.second].mem->begin - hit_graph.nodes[b.first].mem->end) - expected_separation;
        return (a_dist < 0 ? -a_dist * forward_multiplier : a_dist) > (b_dist < 0 ? -b_dist * forward_multiplier : b_dist);
    };
    
    // establish the initial heap ordering
    make_heap(next_comparisons.begin(), next_comparisons.end(), priority_cmp);
    
    // we will block off seeds as they become incorporated into clusters
    // pairs indicate whether a node is blocked for edges (into, out of) it
    vector<pair<bool, bool>> blocked(hit_graph.nodes.size(), pair<bool, bool>(false, false));
    
    // iterate through the comparisons
    while (!next_comparisons.empty())  {
        
        pop_heap(next_comparisons.begin(), next_comparisons.end(), priority_cmp);
        auto comparison = next_comparisons.back();
        next_comparisons.pop_back();
        
#ifdef debug_mem_clusterer
        cerr << "greedy cluster comparing:" << endl;
        cerr << "\t" << comparison.first << ": " << hit_graph.nodes[comparison.first].start_pos << " " << hit_graph.nodes[comparison.first].mem->sequence() << endl;
        cerr << "\t" << comparison.second << ": " << hit_graph.nodes[comparison.second].start_pos << " " << hit_graph.nodes[comparison.second].mem->sequence() << endl;
#endif
        
        if (blocked[comparison.first].second) {
            // we've already greedily accumulated out of this match, so we don't
            // need to look for more connections from it
            
            // TODO: this is actually not correct, we should try to add the next one here...
            continue;
        }
        
        auto& hit_node_1 = hit_graph.nodes[comparison.first];
        auto& hit_node_2 = hit_graph.nodes[comparison.second];
        
        int64_t read_dist = hit_node_2.mem->begin - hit_node_1.mem->end;
        
        if (!blocked[comparison.second].first) {
            
            // what is the minimum distance between these hits?
            size_t min_dist = minimum_distance(*distance_index, hit_node_1.start_pos, hit_node_2.start_pos);
            
#ifdef debug_mem_clusterer
            cerr << "read dist: " << read_dist << ", min dist: " << min_dist << ", graph dist: " << min_dist - (hit_node_1.mem->end - hit_node_1.mem->begin) << endl;
#endif
            
            // TODO: i'm ignoring sub-matches here because it's intended to be used with the stripped
            // algorithm. that might come back to haunt me later
            
            if (min_dist != std::numeric_limits<size_t>::max()) {
                // we were able to measure a distance
                
                // how long of an insert/deletion could we detect based on the scoring parameters?
                int64_t longest_gap = min<int64_t>(min(aligner->longest_detectable_gap(alignment, hit_node_1.mem->end),
                                                       aligner->longest_detectable_gap(alignment, hit_node_2.mem->begin)),
                                                   max_gap);
                
                // the distance from the end of the first hit to the beginning of the next
                int64_t graph_dist = (int64_t)min_dist - (hit_node_1.mem->end - hit_node_1.mem->begin);
                
                // is it possible that an alignment containing both could be detected with local alignment?
                if (abs(read_dist - graph_dist) < longest_gap) {
                    // there's a path within in the limit
                    
#ifdef debug_mem_clusterer
                    cerr << "found hit edge" << endl;
#endif
                    
                    // add the corresponding edge
                    hit_graph.add_edge(comparison.first, comparison.second,
                                       estimate_edge_score(hit_node_1.mem, hit_node_2.mem, graph_dist, aligner),
                                       graph_dist);
                    
                    // we won't look for any more connections involving this end of these two
                    blocked[comparison.first].second = true;
                    blocked[comparison.second].first = true;
                }
            }
        }
        
        if (!blocked[comparison.first].second) {
            
            // TODO: this is actually not correct, we should try to continue until finding an unblocked
            // node
            
            // we didn't just block off connections out of this match,
            // so we can queue up the next one
            if (read_dist >= 0 && comparison.second + 1 < hit_graph.nodes.size() &&
                hit_graph.nodes[comparison.second + 1].mem->begin - hit_node_1.mem->end <= max_separation) {
                // the next in the forward direction
                next_comparisons.emplace_back(comparison.first, comparison.second + 1);
                push_heap(next_comparisons.begin(), next_comparisons.end(), priority_cmp);
            }
            else if (read_dist < 0 && hit_graph.nodes[comparison.second - 1].mem->begin > hit_node_1.mem->begin &&
                     hit_graph.nodes[comparison.second - 1].mem->begin - hit_node_1.mem->end >= min_separation) {
                // the next in the backward direction, requiring read colinearity
                next_comparisons.emplace_back(comparison.first, comparison.second - 1);
                push_heap(next_comparisons.begin(), next_comparisons.end(), priority_cmp);
            }
        }
    }
    
    return hit_graph;
}

ComponentMinDistanceClusterer::ComponentMinDistanceClusterer(SnarlDistanceIndex* distance_index) : MinDistanceClusterer(distance_index) {
    // nothing else to do
}

MEMClusterer::HitGraph ComponentMinDistanceClusterer::make_hit_graph(const Alignment& alignment, const vector<MaximalExactMatch>& mems,
                                                                     const GSSWAligner* aligner, size_t min_mem_length,
                                                                     const match_fanouts_t* fanouts) {
    
    // init the hit graph's nodes
    HitGraph hit_graph(mems, alignment, aligner, min_mem_length, false, fanouts);
    
    // shim the hit graph nodes into the seed clusterer algorithm interface
    vector<SnarlDistanceIndexClusterer::Seed> positions(hit_graph.nodes.size());
    for (size_t i = 0; i < hit_graph.nodes.size(); ++i)  {
        positions[i].pos = hit_graph.nodes[i].start_pos;
    }
 
    typedef SnarlDistanceIndexClusterer::Cluster Cluster;
    SnarlDistanceIndexClusterer seed_clusterer(*distance_index);
    // TODO: magic number, want enough space for the max gap and the inter-seed distance but how to do this in
    // a principled way?
    std::vector<Cluster> distance_components = seed_clusterer.cluster_seeds(positions, 2 * max_gap);
    
    // these components are returned by the structures::UnionFind::all_groups() method, which
    // always returns them in sorted order, so we can assume that they are still lexicographically
    // ordered internally
    for (Cluster& cluster : distance_components) {
        std::vector<size_t>& component = cluster.seeds;
#ifdef debug_mem_clusterer
        cerr << "looking edges in distance component containing:" << endl;
        for (size_t i : component) {
            cerr << "\t" << i << " " << hit_graph.nodes[i].start_pos << " " << hit_graph.nodes[i].mem->sequence() << endl;
        }
#endif
        
        for (size_t i = 0, j_begin = 1; i < component.size(); ++i) {
            
            HitNode& hit_node_1 = hit_graph.nodes[component[i]];
            auto from = hit_node_1.mem->begin;
            
            // start either at the previous target or one past i
            j_begin = max(j_begin, i + 1);
            
            // move forward until we are within the window
            while (j_begin < component.size() &&
                   hit_graph.nodes[component[j_begin]].mem->begin - from < min_read_separation) {
                ++j_begin;
            }
            
            // move backward until we are just within the window
            while (j_begin > i + 1 &&
                   hit_graph.nodes[component[j_begin - 1]].mem->begin - from >= min_read_separation) {
                --j_begin;
            }
            
            int64_t connections_made = 0;
            
            for (size_t j = j_begin;
                 j < component.size() && hit_graph.nodes[component[j]].mem->begin - from <= max_gap; ++j) {
                
                HitNode& hit_node_2 = hit_graph.nodes[component[j]];
                
                // TODO: this code is getting repetitive, i should probably factor it into a MinDistanceClusterer method
                
                int64_t min_dist = minimum_distance(*distance_index, hit_node_1.start_pos, hit_node_2.start_pos);
                if (min_dist != std::numeric_limits<size_t>::max()) {
                    // how long of an insert/deletion could we detect based on the scoring parameters?
                    int64_t longest_gap = min<int64_t>(min(aligner->longest_detectable_gap(alignment, hit_node_1.mem->end),
                                                           aligner->longest_detectable_gap(alignment, hit_node_2.mem->begin)),
                                                       max_gap);
                    
                    // the distance from the end of the first hit to the beginning of the next
                    int64_t graph_dist = (int64_t)min_dist - (hit_node_1.mem->end - hit_node_1.mem->begin);
                    
                    // the distance between the seeds on the read
                    int64_t read_dist = hit_node_2.mem->begin - hit_node_1.mem->end;
                    
                    if (min_dist == hit_node_2.mem->begin - hit_node_1.mem->begin &&
                        ((hit_node_2.mem->begin >= hit_node_1.mem->begin && hit_node_2.mem->end <= hit_node_1.mem->end)
                         || (hit_node_1.mem->begin >= hit_node_2.mem->begin && hit_node_1.mem->end <= hit_node_2.mem->end))) {
                        // this has the appearance of being a redundant hit of a sub-MEM, which we don't want to form
                        // a separate cluster
                        
#ifdef debug_mem_clusterer
                        cerr << "adding dummy hit edge edge " << component[i] << " -> " << component[j] << " to join componenents" << endl;
#endif
                        
                        // we add a dummy edge, but only to connect the nodes' components and join the clusters,
                        // not to actually use in dynamic programming (given arbitrary low weight that should not
                        // cause overflow)
                        hit_graph.add_edge(component[i], component[j], numeric_limits<int32_t>::lowest() / 2, graph_dist);
                    }
                    else if (abs(read_dist - graph_dist) < longest_gap) {
                        // there's a path within in the limit
                        
#ifdef debug_mem_clusterer
                        cerr << "adding hit edge " << component[i] << " -> " << component[j] << ", read dist " << read_dist << ", graph dist " << graph_dist << endl;
#endif
                        
                        // add the corresponding edge
                        hit_graph.add_edge(component[i], component[j],
                                           estimate_edge_score(hit_node_1.mem, hit_node_2.mem, graph_dist, aligner),
                                           graph_dist);
                        
                        // check if we've made enough connnections to stop early (usually the very next
                        // match along the read is the true "next" match in the cluster)
                        ++connections_made;
                        if (early_stop_number && connections_made >= early_stop_number) {
                            break;
                        }
                    }
                }
            }
        }
    }
    
    return hit_graph;
}
    
// collect node starts to build out graph
vector<pair<gcsa::node_type, size_t> > mem_node_start_positions(const HandleGraph& graph, const vg::MaximalExactMatch& mem) {
    // walk the match, getting all the nodes that it touches
    string mem_seq = mem.sequence();
    vector<pair<gcsa::node_type, size_t> > positions;
    hash_set<pair<gcsa::node_type, size_t> > seen_pos;
    //hash_set<handle_t> handles;
    hash_set<pair<gcsa::node_type, size_t> > next; // handles and offsets of the exact match set we need to extend
    auto start_pos = mem.nodes.front();
    next.insert(make_pair(start_pos, 0));
    while (!next.empty()) {
        hash_set<pair<gcsa::node_type, size_t> > todo;
        for (auto& h : next) {
            auto& pos = h.first;
            size_t query_offset = h.second;
            // check if we match each node in next
            auto handle = graph.get_handle(gcsa::Node::id(pos), gcsa::Node::rc(pos));
            string h_seq = graph.get_sequence(handle);
            size_t mem_todo = mem_seq.size() - query_offset;
            size_t overlap = min((size_t)mem_todo, (size_t)(h_seq.size()-gcsa::Node::offset(pos)));
            /*
            cerr << pos << " " << mem_todo << " " << overlap << endl
                 << mem_seq.substr(query_offset, overlap) << endl
                 << h_seq.substr(offset(pos), overlap) << endl;
            */
            // if we do, insert into nodes
            if (mem_seq.substr(query_offset, overlap) == h_seq.substr(gcsa::Node::offset(pos), overlap)) {
                if (!seen_pos.count(h)) {
                    seen_pos.insert(h);
                    auto q = h;
                    q.second = mem_todo - overlap;
                    positions.push_back(q);
                }
                // if we continue past this node, insert our next nodes into nexts
                if (mem_todo - overlap > 0) {
                    size_t new_off = query_offset + overlap;
                    graph.follow_edges(handle, false, [&](const handle_t& next) {
                            todo.insert(make_pair(gcsa::Node::encode(graph.get_id(next), 0, graph.get_is_reverse(next)), new_off));
                            return true;
                        });
                }
            } else {
                // still store at least this node and the remainder
                // but record that we have more left
                if (!seen_pos.count(h)) {
                    seen_pos.insert(h);
                    auto q = h;
                    q.second = mem_todo;
                    positions.push_back(q);
                }
            }
        }
        next = todo;
    }
    // ensure positions are sorted by remainder
    std::sort(positions.begin(), positions.end(),
              [](const pair<gcsa::node_type, size_t>& a, const pair<gcsa::node_type, size_t>& b) {
                  return a.second > b.second;
              });
    return positions;
}

bdsg::HashGraph cluster_subgraph_containing(const HandleGraph& base, const Alignment& aln, const vector<vg::MaximalExactMatch>& cluster, const GSSWAligner* aligner) {
    vector<pos_t> positions;
    vector<size_t> forward_max_dist;
    vector<size_t> backward_max_dist;
    positions.reserve(cluster.size());
    forward_max_dist.reserve(cluster.size());
    backward_max_dist.reserve(cluster.size());
    for (auto& mem : cluster) {
        // get the start position of the MEM
        positions.push_back(make_pos_t(mem.nodes.front()));
        // search far enough away to get any hit detectable without soft clipping
        forward_max_dist.push_back(aligner->longest_detectable_gap(aln, mem.end)
                                   + (aln.sequence().end() - mem.begin));
        backward_max_dist.push_back(aligner->longest_detectable_gap(aln, mem.begin)
                                    + (mem.begin - aln.sequence().begin()));
    }
    auto cluster_graph = bdsg::HashGraph();
    algorithms::extract_containing_graph(&base, &cluster_graph, positions, forward_max_dist, backward_max_dist);
    return cluster_graph;
}

bdsg::HashGraph cluster_subgraph_walk(const HandleGraph& base, const Alignment& aln, const vector<vg::MaximalExactMatch>& mems, double expansion) {
    assert(mems.size());
    auto& start_mem = mems.front();
    auto start_pos = make_pos_t(start_mem.nodes.front());
    auto rev_start_pos = reverse(start_pos, base.get_length(base.get_handle(id(start_pos))));
    // Even if the MEM is right up against the start of the read, it may not be
    // part of the best alignment. Make sure to have some padding.
    // TODO: how much padding?
    bdsg::HashGraph graph;
    int inside_padding = max(1, (int)aln.sequence().size()/16);
    int end_padding = max(8, (int)aln.sequence().size()/8);
    int get_before = end_padding + (int)(expansion * (int)(start_mem.begin - aln.sequence().begin()));
    if (get_before) {
        //algorithms::extract_context(base, graph, base.get_handle(id(rev_start_pos), is_rev(rev_start_pos)), offset(rev_start_pos), get_before, false, true);
    }
    //cerr << "======================================================" << endl;
    for (int i = 0; i < mems.size(); ++i) {
        auto& mem = mems[i];
        //cerr << mem << endl;
        vector<pair<gcsa::node_type, size_t> > match_positions = mem_node_start_positions(base, mem);
        if (!match_positions.size()) {
            // TODO XXX is MEM merging causing this to occur?
            match_positions.push_back(make_pair(mem.nodes.front(), mem.length()));
        }
        for (auto& p : match_positions) {
            handle_t h = base.get_handle(gcsa::Node::id(p.first), gcsa::Node::rc(p.first));
            algorithms::extract_context(base, graph, h, gcsa::Node::offset(p.first), 0);
        }
        // extend after the last match node with the expansion
        auto& p = match_positions.back();
        auto& pos = p.first;
        int mem_remainder = p.second;
        //cerr << p.first << " " << p.second << endl;
        handle_t h = base.get_handle(gcsa::Node::id(pos));
        int get_after = //base.get_length(h);
            (i+1 == mems.size() ?
               end_padding +
               expansion * ((int)(aln.sequence().end() - mem.end) + mem_remainder)
               :
               inside_padding +
               expansion * ((int)(mems[i+1].begin - mem.end) + mem_remainder));
        if (get_after > 0) {
            algorithms::extract_context(base, graph, h, gcsa::Node::offset(pos), get_after, true, false);
        }
    }
    //algorithms::expand_subgraph_by_steps(base, graph, 0);
    algorithms::expand_subgraph_to_length(base, graph, aln.sequence().size() * expansion, false);
    return graph;
}

}
