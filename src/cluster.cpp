#include <string>
#include <algorithm>
#include <utility>
#include <cstring>

#include "cluster.hpp"

//#define debug_od_clusterer

namespace vg {
    
using namespace std;

MEMChainModel::MEMChainModel(
    const vector<size_t>& aln_lengths,
    const vector<vector<MaximalExactMatch> >& matches,
    const function<int64_t(pos_t)>& approx_position,
    const function<map<string, vector<size_t> >(pos_t)>& path_position,
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
                m.mem.positions[""].push_back(approx_position(pos));
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
                chr_positions[pos].push_back(v);
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
    // for each vertex merge if we go equivalently forward in the positional space and forward in the read to the next position
    // scan forward
    for (map<string, map<int64_t, vector<vector<MEMChainModelVertex>::iterator> > >::iterator c = positions.begin(); c != positions.end(); ++c) {
        for (map<int64_t, vector<vector<MEMChainModelVertex>::iterator> >::iterator p = c->second.begin(); p != c->second.end(); ++p) {
            for (auto& v1 : p->second) {
                if (redundant_vertexes.count(v1)) continue;
                auto q = p;
                while (++q != c->second.end() && abs(p->first - q->first) < band_width) {
                    for (auto& v2 : q->second) {
                        if (redundant_vertexes.count(v2)) continue;
                        if (mems_overlap(v1->mem, v2->mem)
                            && abs(v2->mem.begin - v1->mem.begin) == abs(q->first - p->first)) {
                            if (v2->mem.length() < v1->mem.length()) {
                                redundant_vertexes.insert(v2);
                                if (v2->mem.end > v1->mem.end) {
                                    v1->weight += v2->mem.end - v1->mem.end;
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    // scan reverse
    for (map<string, map<int64_t, vector<vector<MEMChainModelVertex>::iterator> > >::iterator c = positions.begin(); c != positions.end(); ++c) {
        for (map<int64_t, vector<vector<MEMChainModelVertex>::iterator> >::reverse_iterator p = c->second.rbegin(); p != c->second.rend(); ++p) {
            for (auto& v1 : p->second) {
                if (redundant_vertexes.count(v1)) continue;
                auto q = p;
                while (++q != c->second.rend() && abs(p->first - q->first) < band_width) {
                    for (auto& v2 : q->second) {
                        if (redundant_vertexes.count(v2)) continue;
                        if (mems_overlap(v1->mem, v2->mem)
                            && abs(v2->mem.begin - v1->mem.begin) == abs(p->first - q->first)) {
                            if (v2->mem.length() < v1->mem.length()) {
                                redundant_vertexes.insert(v2);
                                if (v2->mem.end > v1->mem.end) {
                                    v1->weight += v2->mem.end - v1->mem.end;
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    // now build up the model using the positional bandwidth
    set<pair<vector<MEMChainModelVertex>::iterator, vector<MEMChainModelVertex>::iterator> > seen;
    for (map<string, map<int64_t, vector<vector<MEMChainModelVertex>::iterator> > >::iterator c = positions.begin(); c != positions.end(); ++c) {
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
                                || v1->mem.fragment == v2->mem.fragment && v1->mem.begin < v2->mem.begin) {
                                // Transition is allowable because the first comes before the second
                            
                                double weight = transition_weight(v1->mem, v2->mem);
                                if (weight > -std::numeric_limits<double>::max()) {
                                    v1->next_cost.push_back(make_pair(&*v2, weight));
                                    v2->prev_cost.push_back(make_pair(&*v1, weight));
                                }
                            } else if (v1->mem.fragment > v2->mem.fragment
                                       || v1->mem.fragment == v2->mem.fragment && v1->mem.begin > v2->mem.begin) {
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

void MEMChainModel::score(const set<MEMChainModelVertex*>& exclude) {
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
    set<MEMChainModelVertex*> exclude;
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
        set<MEMChainModelVertex*> chain_members;
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

OrientedDistanceClusterer::OrientedDistanceClusterer(const Alignment& alignment,
                                                     const vector<MaximalExactMatch>& mems,
                                                     const QualAdjAligner& aligner,
                                                     xg::XG* xgindex,
                                                     size_t max_expected_dist_approx_error,
                                                     size_t min_mem_length,
                                                     node_occurrence_on_paths_memo_t* node_path_memo,
                                                     handle_memo_t* handle_memo) :
    OrientedDistanceClusterer(alignment, mems, nullptr, &aligner, xgindex, max_expected_dist_approx_error,
                              min_mem_length, node_path_memo, handle_memo) {
    // nothing else to do
}

OrientedDistanceClusterer::OrientedDistanceClusterer(const Alignment& alignment,
                                                     const vector<MaximalExactMatch>& mems,
                                                     const Aligner& aligner,
                                                     xg::XG* xgindex,
                                                     size_t max_expected_dist_approx_error,
                                                     size_t min_mem_length,
                                                     node_occurrence_on_paths_memo_t* node_path_memo,
                                                     handle_memo_t* handle_memo) :
    OrientedDistanceClusterer(alignment, mems, &aligner, nullptr, xgindex, max_expected_dist_approx_error,
                              min_mem_length, node_path_memo, handle_memo) {
    // nothing else to do
}

OrientedDistanceClusterer::OrientedDistanceClusterer(const Alignment& alignment,
                                                     const vector<MaximalExactMatch>& mems,
                                                     const Aligner* aligner,
                                                     const QualAdjAligner* qual_adj_aligner,
                                                     xg::XG* xgindex,
                                                     size_t max_expected_dist_approx_error,
                                                     size_t min_mem_length,
                                                     node_occurrence_on_paths_memo_t* node_path_memo,
                                                     handle_memo_t* handle_memo) : aligner(aligner), qual_adj_aligner(qual_adj_aligner) {
    
    // there generally will be at least as many nodes as MEMs, so we can speed up the reallocation
    nodes.reserve(mems.size());
    
    for (const MaximalExactMatch& mem : mems) {
        
        //#pragma omp atomic
        //        MEM_TOTAL += mem.nodes.size();
        
        if (mem.length() < min_mem_length) {
#ifdef debug_od_clusterer
            cerr << "skipping short MEM " << mem << endl;
#endif
            //#pragma omp atomic
            //            MEM_FILTER_COUNTER += mem.nodes.size();
            continue;
        }
        
        // calculate the longest gaps we could detect to the left and right of this MEM
        int32_t mem_score;
        if (aligner) {
            mem_score = aligner->score_exact_match(mem.begin, mem.end);
        }
        else {
            mem_score = qual_adj_aligner->score_exact_match(mem.begin, mem.end, alignment.quality().begin()
                                                            + (mem.begin - alignment.sequence().begin()));
        }
        
#ifdef debug_od_clusterer
        cerr << "adding nodes for MEM " << mem << endl;
#endif
        for (gcsa::node_type mem_hit : mem.nodes) {
            nodes.emplace_back(mem, make_pos_t(mem_hit), mem_score);
#ifdef debug_od_clusterer
            cerr << "\t" << nodes.size() - 1 << ": " << make_pos_t(mem_hit) << endl;
#endif
        }
    }
    
    // Get all the distances between nodes, in a forrest of unrooted trees of
    // nodes that we know are on a consistent strand.
    unordered_map<pair<size_t, size_t>, int64_t> recorded_finite_dists = get_on_strand_distance_tree(nodes.size(), xgindex,
                                                                                                     [&](size_t node_number) {
                                                                                                         return nodes[node_number].start_pos;
                                                                                                     },
                                                                                                     [&](size_t node_number) {
                                                                                                         return 0;
                                                                                                     },
                                                                                                     node_path_memo,
                                                                                                     handle_memo);
    
    // Flatten the trees to maps of relative position by node ID.
    vector<unordered_map<size_t, int64_t>> strand_relative_position = flatten_distance_tree(nodes.size(), recorded_finite_dists);
    
#ifdef debug_od_clusterer
    for (const auto& strand : strand_relative_position) {
        cerr << "strand reconstruction: "  << endl;
        for (const auto& record : strand) {
            cerr << "\t" << record.first << ": " << record.second << "\t" << nodes[record.first].mem->sequence() << endl;
        }
    }
#endif
    
    // now we use the strand clusters and the estimated distances to make the DAG for the
    // approximate MEM alignment
    
    int64_t match_score, mismatch_score, gap_open_score, gap_extension_score, max_gap;
    if (aligner) {
        match_score = aligner->match;
        gap_open_score = aligner->gap_open;
        gap_extension_score = aligner->gap_extension;
        max_gap = aligner->longest_detectable_gap(alignment);
    }
    else {
        match_score = qual_adj_aligner->match;
        gap_open_score = qual_adj_aligner->gap_open;
        gap_extension_score = qual_adj_aligner->gap_extension;
        max_gap = qual_adj_aligner->longest_detectable_gap(alignment);
    }
    
    int64_t forward_gap_length = max_gap + max_expected_dist_approx_error;
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
        int64_t low = 0, hi = last_idx;
        for (int64_t i = 0; i < sorted_pos.size(); i++) {
            
            int64_t strand_pos = sorted_pos[i].first;
            size_t pivot_idx = sorted_pos[i].second;
            ODNode& pivot = nodes[pivot_idx];
            int64_t pivot_length = pivot.mem->end - pivot.mem->begin;
            int64_t suffix_length = alignment.sequence().end() - pivot.mem->end;
            
            // the limits of how far away we might detect edges to add to the clustering graph
            int64_t target_low_pos = strand_pos - max_expected_dist_approx_error;
            int64_t target_hi_pos = strand_pos + suffix_length + forward_gap_length;
            
            // move the lower boundary of the search interval to the lowest value inside the
            // the target interval
            while (low < sorted_pos.size() ? sorted_pos[low].first < target_low_pos : false) {
                low++;
            }
            
            // move the upper boundary of the search interval to the highest value inside the
            // the target interval (this one can move in either direction because pivots are
            // different lengths)
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
            
#ifdef debug_od_clusterer
            cerr << "checking for possible edges from " << sorted_pos[i].second << " to MEMs between " << sorted_pos[low].first << "(" << sorted_pos[low].second << ") and " << sorted_pos[hi].first << "(" << sorted_pos[hi].second << "), which is inside the interval (" << target_low_pos << ", " << target_hi_pos << ")" << endl;
#endif
            
            for (int64_t j = low; j <= hi; j++) {
                int64_t next_idx = sorted_pos[j].second;
                ODNode& next = nodes[next_idx];
                
                if (next.mem->begin < pivot.mem->begin || next.mem->end < pivot.mem->end
                    || (next.mem->begin == pivot.mem->begin && next.mem->end == pivot.mem->end)) {
                    // these MEMs cannot be colinear along the read (also filters out j == i)
                    
                    // note: we allow one of the start/end positions to be the same here even though they can't
                    // techinically overlap because it tends to soak up redundant sub-MEMs into the same connected
                    // component so that they don't get their own cluster
                    
                    continue;
                }
                
                // the length of the sequence in between the MEMs (can be negative if they overlap)
                int64_t between_length = next.mem->begin - pivot.mem->end;
                // the estimated distance between the end of the pivot and the start of the next MEM in the graph
                int64_t graph_dist = sorted_pos[j].first - strand_pos - pivot_length;
                
                int32_t edge_score;
                if (between_length < 0) {
                    // the MEMs overlap, but this can occur in some insertions and deletions
                    // because the SMEM algorithm is "greedy" in taking up as much of the read
                    // as possible
                    // we can check if this happened directly, but it's expensive
                    // so for now we just give it the benefit of the doubt but adjust the edge
                    // score so that the matches don't get double counted
                    
                    int64_t extra_dist = abs(graph_dist - between_length);
                    
                    edge_score = match_score * between_length
                                 - (extra_dist ? (extra_dist - 1) * gap_extension_score + gap_open_score : 0);
                }
                else {
                    int64_t gap_length = abs(between_length - graph_dist);
                    // the read length in between the MEMs is the same as the distance, suggesting a pure mismatch
                    edge_score = gap_length ? -((gap_length - 1) * gap_extension_score + gap_open_score) : 0;
                }
                
#ifdef debug_od_clusterer
                cerr << "adding edge to MEM " << sorted_pos[j].first << "(" << sorted_pos[j].second << ") with weight " << edge_score << endl;
#endif
                
                // add the edges in
                pivot.edges_from.emplace_back(next_idx, edge_score);
                next.edges_to.emplace_back(pivot_idx, edge_score);
            }
        }
    }
}

unordered_map<pair<size_t, size_t>, int64_t> OrientedDistanceClusterer::get_on_strand_distance_tree(size_t num_items, xg::XG* xgindex,
                                                                                                    const function<pos_t(size_t)>& get_position,
                                                                                                    const function<int64_t(size_t)>& get_offset,
                                                                                                    node_occurrence_on_paths_memo_t* node_path_memo,
                                                                                                    handle_memo_t* handle_memo) {
    
    // for recording the distance of any pair that we check with a finite distance
    unordered_map<pair<size_t, size_t>, int64_t> recorded_finite_dists;
    
    // we use a union find to keep track of which MEMs have been identified as being on the same strand
    UnionFind component_union_find(num_items);
    
    size_t num_possible_merges_remaining = (num_items * (num_items - 1)) / 2;
    
    // an initial pass that only looks at nodes on path
    extend_dist_tree_by_path_buckets(num_possible_merges_remaining,component_union_find, recorded_finite_dists,
                                     num_items, xgindex, get_position, get_offset, node_path_memo, handle_memo);
    
    // TODO: permutations that try to assign singletons
    
    // a second pass that tries fill in the tree by traversing to the nearest shared path
    size_t nlogn = ceil(num_items * log(num_items));
    extend_dist_tree_by_permutations(2, 50, nlogn, num_possible_merges_remaining, component_union_find, recorded_finite_dists,
                                     num_items, xgindex, get_position, get_offset, node_path_memo, handle_memo);
    
    return recorded_finite_dists;
}
    
void OrientedDistanceClusterer::extend_dist_tree_by_path_buckets(size_t& num_possible_merges_remaining,
                                                                 UnionFind& component_union_find,
                                                                 unordered_map<pair<size_t, size_t>, int64_t>& recorded_finite_dists,
                                                                 size_t num_items,
                                                                 xg::XG* xgindex,
                                                                 const function<pos_t(size_t)>& get_position,
                                                                 const function<int64_t(size_t)>& get_offset,
                                                                 node_occurrence_on_paths_memo_t* paths_of_node_memo,
                                                                 handle_memo_t* handle_memo) {
    if (!paths_of_node_memo) {
        return;
    }
    
    for (size_t i = 0; i < num_items; i++) {
        pos_t pos = get_position(i);
        if (!paths_of_node_memo->count(id(pos))) {
            (*paths_of_node_memo)[id(pos)] = xgindex->oriented_paths_of_node(id(pos));
        }
    }
    
    // identify singletons and reverse the memo so that it tells us which hits occur on a strand of a path
    unordered_map<pair<size_t, bool>, vector<size_t>> items_on_path_strand;
    vector<size_t> singletons;
    for (size_t i = 0; i < num_items; i++) {
        pos_t pos = get_position(i);
        auto& path_records = paths_of_node_memo->at(id(pos));
        if (path_records.empty()) {
            singletons.push_back(i);
        }
        else {
            for (pair<size_t, vector<pair<size_t, bool>>>& path_record : path_records) {
                for (pair<size_t, bool>& node_occurence : path_record.second) {
                    items_on_path_strand[make_pair(path_record.first, node_occurence.second != is_rev(pos))].push_back(i);
                }
            }
        }
    }
    
    // check the nearest nodes to each singleton to see if we can use it to bucket the item
    for (size_t i : singletons) {
        pos_t pos = get_position(i);
        handle_t handle = xgindex->memoized_get_handle(id(pos), is_rev(pos), handle_memo);
        size_t right_dist = xgindex->get_length(handle) - offset(pos);
        size_t trav_dist = min(offset(pos), right_dist);
        // TODO: magic number (matches the distance used in the permutations step)
        if (trav_dist <= 50) {
            bool go_left = offset(pos) < right_dist;
            function<bool(const handle_t&)> bucket_using_neighbors = [&](const handle_t& handle) {
                id_t neighbor_id = xgindex->get_id(handle);
                bool neighbor_rev = xgindex->get_is_reverse(handle);
                if (!paths_of_node_memo->count(neighbor_id)) {
                    (*paths_of_node_memo)[neighbor_id] = xgindex->oriented_paths_of_node(neighbor_id);
                }
                auto& path_records = paths_of_node_memo->at(neighbor_id);
                for (pair<size_t, vector<pair<size_t, bool>>>& path_record : path_records) {
                    for (pair<size_t, bool>& node_occurence : path_record.second) {
                        items_on_path_strand[make_pair(path_record.first, node_occurence.second != neighbor_rev)].push_back(i);
                    }
                }
                return true;
            };
            
            xgindex->follow_edges(handle, go_left, bucket_using_neighbors);
        }
    }
    
    
    // make sure the items are unique with each list of hits and generate a system-independent ordering over strands
    vector<pair<size_t, bool>> buckets;
    buckets.reserve(items_on_path_strand.size());
    for (pair<const pair<size_t, bool>, vector<size_t>>& strand_bucket : items_on_path_strand) {
        sort(strand_bucket.second.begin(), strand_bucket.second.end());
        auto new_end = unique(strand_bucket.second.begin(), strand_bucket.second.end());
        strand_bucket.second.resize(new_end - strand_bucket.second.begin());
        buckets.push_back(strand_bucket.first);
    }
    sort(buckets.begin(), buckets.end());
    
    // use the path strands to bucket distance measurements
    for (pair<size_t, bool>& strand_bucket : buckets) {
#ifdef debug_od_clusterer
        cerr << "doing a bucketed comparison of items on the path ranked " << strand_bucket.first << ", strand " << (strand_bucket.second ? "-" : "+") << endl;
#endif
        vector<size_t>& bucket = items_on_path_strand[strand_bucket];
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
            
#ifdef debug_od_clusterer
            cerr << "measuring distance between " << prev << " at " << pos_prev << " and " << here << " at " << pos_here << endl;
#endif
            
            int64_t dist = xgindex->closest_shared_path_oriented_distance(id(pos_prev), offset(pos_prev), is_rev(pos_prev),
                                                                          id(pos_here), offset(pos_here), is_rev(pos_here),
                                                                          0, paths_of_node_memo, handle_memo);
            
            
            // did we get a successful estimation?
            if (dist == numeric_limits<int64_t>::max()) {
#ifdef debug_od_clusterer
                cerr << "they don't appear to be on the same path, skipping" << endl;
#endif
                continue;
            }
            
            // add the fixed offset from the hit position
            dist += get_offset(here) - get_offset(prev);
            
#ifdef debug_od_clusterer
            cerr << "recording distance at " << dist << endl;
#endif
            
            // merge them into a strand cluster
            recorded_finite_dists[make_pair(prev, here)] = dist;
            num_possible_merges_remaining -= component_union_find.group_size(prev) * component_union_find.group_size(here);
            component_union_find.union_groups(prev, here);
        }
    }
}
    
void OrientedDistanceClusterer::extend_dist_tree_by_permutations(int64_t max_failed_distance_probes,
                                                                 int64_t max_search_distance_to_path,
                                                                 size_t decrement_frequency,
                                                                 size_t& num_possible_merges_remaining,
                                                                 UnionFind& component_union_find,
                                                                 unordered_map<pair<size_t, size_t>, int64_t>& recorded_finite_dists,
                                                                 size_t num_items,
                                                                 xg::XG* xgindex,
                                                                 const function<pos_t(size_t)>& get_position,
                                                                 const function<int64_t(size_t)>& get_offset,
                                                                 node_occurrence_on_paths_memo_t* paths_of_node_memo,
                                                                 handle_memo_t* handle_memo) {
    
    // for recording the number of times elements of a strand cluster have been compared
    // and found an infinite distance
    map<pair<size_t, size_t>, size_t> num_infinite_dists;
    
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
#ifdef debug_od_clusterer
        size_t direct_merges_remaining = 0;
        vector<vector<size_t>> groups = component_union_find.all_groups();
        for (size_t i = 1; i < groups.size(); i++) {
            for (size_t j = 0; j < i; j++) {
                size_t strand_1 = component_union_find.find_group(groups[i].front());
                size_t strand_2 = component_union_find.find_group(groups[j].front());
                
                if (num_infinite_dists.count(make_pair(strand_1, strand_2))) {
                    if (num_infinite_dists[make_pair(strand_1, strand_2)] >= current_max_num_probes) {
                        continue;
                    }
                }
                direct_merges_remaining += groups[i].size() * groups[j].size();
            }
        }
        cerr << "checked " << pairs_checked << " pairs with max probes " << current_max_num_probes << ", decrement frequency " << decrement_frequency << ", directly calculated merges " << direct_merges_remaining << " and maintained merges " << num_possible_merges_remaining << endl;
#endif
        
        if (pairs_checked % decrement_frequency == 0 && pairs_checked != 0) {
            current_max_num_probes--;
#ifdef debug_od_clusterer
            cerr << "reducing the max number of probes to " << current_max_num_probes << endl;
#endif
            for (const pair<pair<size_t, size_t>, size_t>& inf_dist_record : num_infinite_dists) {
                // break symmetry so we don't repeat the operation twice
                if (inf_dist_record.first.first < inf_dist_record.first.second && inf_dist_record.second == current_max_num_probes) {
                    // this merge just fell below the new maximum number of distance probes
                    size_t strand_size_1 = component_union_find.group_size(inf_dist_record.first.first);
                    size_t strand_size_2 = component_union_find.group_size(inf_dist_record.first.second);
                    num_possible_merges_remaining -= strand_size_1 * strand_size_2;
#ifdef debug_od_clusterer
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
        
#ifdef debug_od_clusterer
        cerr << "checking MEMs " << node_pair.first << " and " << node_pair.second << " in cluster " << strand_1 << " and " << strand_2 << endl;
#endif
        
        if (strand_1 == strand_2) {
            // these are already identified as on the same strand, don't need to do it again
#ifdef debug_od_clusterer
            cerr << "already on same strand" << endl;
#endif
            continue;
        }
        
        auto num_failed_probes = num_infinite_dists.find(make_pair(strand_1, strand_2));
        if (num_failed_probes == num_infinite_dists.end() ? false : num_failed_probes->second >= current_max_num_probes) {
            // we've already checked multiple distances between these strand clusters and
            // none have returned a finite distance, so we conclude that they are in fact
            // on separate clusters and decline to check any more distances
#ifdef debug_od_clusterer
            cerr << "already have checked distance above maximum number of probes" << endl;
#endif
            continue;
        }
        
        const pos_t& pos_1 = get_position(node_pair.first);
        const pos_t& pos_2 = get_position(node_pair.second);
        
        int64_t oriented_dist = xgindex->closest_shared_path_oriented_distance(id(pos_1), offset(pos_1), is_rev(pos_1),
                                                                               id(pos_2), offset(pos_2), is_rev(pos_2),
                                                                               max_search_distance_to_path, paths_of_node_memo,
                                                                               handle_memo);
        
#ifdef debug_od_clusterer
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
                
#ifdef debug_od_clusterer
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
            
#ifdef debug_od_clusterer
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
                        
#ifdef debug_od_clusterer
                        cerr << "after merge, the total number of probes against strand " << removing_iter->first.second << " increased to " << retaining_iter->second << ", above current max of " << current_max_num_probes << ", but the retaining strand is already blocked, reducing possible merges by " << (strand_retaining == strand_1 ? strand_size_2 : strand_size_1) * component_union_find.group_size(removing_iter->first.second) << " to " << num_possible_merges_remaining << endl;
#endif
                    }
                    else if (removing_already_blocked && !retaining_already_blocked) {
                        num_possible_merges_remaining -= (strand_retaining == strand_1 ? strand_size_1 : strand_size_2) * component_union_find.group_size(removing_iter->first.second);
                        
#ifdef debug_od_clusterer
                        cerr << "after merge, the total number of probes against strand " << removing_iter->first.second << " increased to " << retaining_iter->second << ", above current max of " << current_max_num_probes << ", but the removing strand is already blocked, reducing possible merges by " << (strand_retaining == strand_1 ? strand_size_1 : strand_size_2) * component_union_find.group_size(removing_iter->first.second) << " to " << num_possible_merges_remaining << endl;
#endif
                    }
                    else if (!retaining_already_blocked && !removing_already_blocked && retaining_iter->second >= current_max_num_probes) {
                        num_possible_merges_remaining -= (strand_size_1 + strand_size_2) * component_union_find.group_size(removing_iter->first.second);
                        
#ifdef debug_od_clusterer
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
                        
#ifdef debug_od_clusterer
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
                    
#ifdef debug_od_clusterer
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
                    
#ifdef debug_od_clusterer
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
    
#ifdef debug_od_clusterer
    cerr << "constructing strand distance tree" << endl;
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
        
#ifdef debug_od_clusterer
        cerr << "beginning a distance tree traversal at MEM " << i << endl;
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
    
#ifdef debug_od_clusterer
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
    
#ifdef debug_od_clusterer
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

void OrientedDistanceClusterer::topological_order(vector<size_t>& order_out) {
    
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

void OrientedDistanceClusterer::identify_sources_and_sinks(vector<size_t>& sources_out,
                                                           vector<size_t>& sinks_out) {
    
    sources_out.clear();
    sinks_out.clear();
    
    vector<bool> is_source(nodes.size(), true);
    
    for (size_t i = 0; i < nodes.size(); i++) {
        if (nodes[i].edges_from.empty()) {
            sinks_out.push_back(i);
        }
        
        for (ODEdge& edge : nodes[i].edges_from) {
            is_source[edge.to_idx] = false;
        }
    }
    
    for (size_t i = 0; i < nodes.size(); i++) {
        if (is_source[i]) {
            sources_out.push_back(i);
        }
    }
}

void OrientedDistanceClusterer::connected_components(vector<vector<size_t>>& components_out) {
    
    components_out.clear();
    vector<bool> enqueued(nodes.size());
    
    // check each node in turn to find new components
    for (size_t dfs_start_idx = 0; dfs_start_idx < nodes.size(); dfs_start_idx++) {
        if (enqueued[dfs_start_idx]) {
            // we've already found this node from some component
            continue;
        }
        
        // this node belongs to a component we haven't found yet, use DFS to find the rest
        vector<size_t> stack {dfs_start_idx};
        enqueued[dfs_start_idx] = true;
        components_out.emplace_back(1, dfs_start_idx);
        
        while (!stack.empty()) {
            
            ODNode& node = nodes[stack.back()];
            stack.pop_back();
            
            // search in both forward and backward directions
            
            for (ODEdge& edge : node.edges_from) {
                
                if (!enqueued[edge.to_idx]) {
                    stack.push_back(edge.to_idx);
                    enqueued[edge.to_idx] = true;
                    components_out.back().push_back(edge.to_idx);
                }
            }
            
            for (ODEdge& edge : node.edges_to) {
                
                if (!enqueued[edge.to_idx]) {
                    stack.push_back(edge.to_idx);
                    enqueued[edge.to_idx] = true;
                    components_out.back().push_back(edge.to_idx);
                }
            }
        }
    }
}

void OrientedDistanceClusterer::perform_dp() {
    
    for (ODNode& node : nodes) {
        // as in local alignment, minimum score is the score of node itself
        node.dp_score = node.score;
    }
    
#ifdef debug_od_clusterer
    cerr << "computing topological order for clustering DP" << endl;
#endif
    
    vector<size_t> order;
    topological_order(order);
    
    for (size_t i : order) {
        ODNode& node = nodes[i];
#ifdef debug_od_clusterer
        cerr << "at node " << i << " with DP score " << node.dp_score << " and node score " << node.score << endl;
#endif
        // for each edge out of this node
        for (ODEdge& edge : node.edges_from) {
            
            // check if the path through the node out of this edge increase score of target node
            ODNode& target_node = nodes[edge.to_idx];
            int32_t extend_score = node.dp_score + edge.weight + target_node.score;
            if (extend_score > target_node.dp_score) {
#ifdef debug_od_clusterer
                cerr << "extending DP to node " << edge.to_idx << " with score " << extend_score << endl;
#endif
                target_node.dp_score = extend_score;
            }
        }
    }
}

vector<OrientedDistanceClusterer::cluster_t> OrientedDistanceClusterer::clusters(int32_t max_qual_score,
                                                                                 int32_t log_likelihood_approx_factor) {
    
    vector<vector<pair<const MaximalExactMatch*, pos_t>>> to_return;
    if (nodes.size() == 0) {
        // this should only happen if we have filtered out all MEMs, so there are none to cluster
        return to_return;
    }
    
#ifdef debug_od_clusterer
    cerr << "performing approximate DP across MEMs" << endl;
#endif
    perform_dp();
    
#ifdef debug_od_clusterer
    cerr << "finding top tracebacks within connected components" << endl;
#endif
    // find the weakly connected components, which should correspond to mappings
    vector<vector<size_t>> components;
    connected_components(components);
    
#ifdef debug_od_clusterer
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
    
    // find the node with the highest DP score in each connected component
    // each record is a pair of (score lower bound, node index)
    vector<pair<int32_t, size_t>> component_traceback_ends(components.size(),
                                                           pair<int32_t, size_t>(numeric_limits<int32_t>::min(), 0));
    for (size_t i = 0; i < components.size(); i++) {
        vector<size_t>& component = components[i];
        pair<int32_t, size_t>& traceback_end = component_traceback_ends[i];
        for (size_t j = 0; j < component.size(); j++) {
            int32_t dp_score = nodes[component[j]].dp_score;
            if (dp_score > traceback_end.first) {
                traceback_end.first = dp_score;
                traceback_end.second = component[j];
            }
        }
    }
    //#pragma omp atomic
    //    CLUSTER_TOTAL += component_traceback_ends.size();
    
    std::make_heap(component_traceback_ends.begin(), component_traceback_ends.end());
    
    // estimate the minimum score a cluster must obtain to even affect the mapping quality
    // TODO: this approximation could break down sometimes, need to look into it
    int32_t top_score = component_traceback_ends.front().first;
    const BaseAligner* base_aligner = aligner ? (BaseAligner*) aligner : (BaseAligner*) qual_adj_aligner;
    int32_t suboptimal_score_cutoff = top_score - log_likelihood_approx_factor * base_aligner->mapping_quality_score_diff(max_qual_score);
    
    while (!component_traceback_ends.empty()) {
        // get the next highest scoring traceback end
        auto traceback_end = component_traceback_ends.front();
        std::pop_heap(component_traceback_ends.begin(), component_traceback_ends.end());
        component_traceback_ends.pop_back();
        
        // get the index of the node
        size_t trace_idx = traceback_end.second;
        
#ifdef debug_od_clusterer
        cerr << "checking traceback of component starting at " << traceback_end.second << endl;
#endif
        // if this cluster does not look like it even affect the mapping quality of the top scoring
        // cluster, don't bother forming it
        if (traceback_end.first < suboptimal_score_cutoff) {
#ifdef debug_od_clusterer
            cerr << "skipping rest of components on account of low score of " << traceback_end.first << " compared to max score " << top_score << " and cutoff " << suboptimal_score_cutoff << endl;
#endif
            
            //#pragma omp atomic
            //            PRUNE_COUNTER += component_traceback_ends.size() + 1;
            break;
        }
        
        // traceback until hitting a node that has its own score (indicates beginning of a local alignment)
        vector<size_t> trace{trace_idx};
        while (nodes[trace_idx].dp_score > nodes[trace_idx].score) {
            int32_t target_source_score = nodes[trace_idx].dp_score - nodes[trace_idx].score;
            for (ODEdge& edge : nodes[trace_idx].edges_to) {
                if (nodes[edge.to_idx].dp_score + edge.weight == target_source_score) {
                    trace_idx = edge.to_idx;
                    trace.push_back(trace_idx);
                    break;
                }
            }
        }
        
        // make a cluster
        to_return.emplace_back();
        auto& cluster = to_return.back();
        for (auto iter = trace.rbegin(); iter != trace.rend(); iter++) {
            ODNode& node = nodes[*iter];
            cluster.emplace_back(node.mem, node.start_pos);
        }
    }
    
    return std::move(to_return);
}

vector<pair<pair<size_t, size_t>, int64_t>> OrientedDistanceClusterer::pair_clusters(const Alignment& alignment_1,
                                                                                     const Alignment& alignment_2,
                                                                                     const vector<cluster_t*>& left_clusters,
                                                                                     const vector<cluster_t*>& right_clusters,
                                                                                     xg::XG* xgindex,
                                                                                     int64_t min_inter_cluster_distance,
                                                                                     int64_t max_inter_cluster_distance,
                                                                                     node_occurrence_on_paths_memo_t* node_path_memo,
                                                                                     handle_memo_t* handle_memo) {
    
#ifdef debug_od_clusterer
    cerr << "beginning clustering of MEM cluster pairs" << endl;
#endif
    
    // We will fill this in with all sufficiently close pairs of clusters from different reads.
    vector<pair<pair<size_t, size_t>, int64_t>> to_return;
    
    // We think of the clusters as a single linear ordering, with our clusters coming first.
    size_t total_clusters = left_clusters.size() + right_clusters.size();
    
    // Compute distance trees for sets of clusters that are distance-able on consistent strands.
    unordered_map<pair<size_t, size_t>, int64_t> distance_tree = get_on_strand_distance_tree(total_clusters, xgindex,
         [&](size_t cluster_num) {
             // Assumes the clusters are nonempty.
             if (cluster_num < left_clusters.size()) {
                 // Grab the pos_t for the first hit in the cluster, which is sorted to be the largest one.
                 return left_clusters[cluster_num]->front().second;
             } else {
                 // Grab the pos_t for the largest hit from the other cluster
                 return right_clusters[cluster_num - left_clusters.size()]->front().second;
             }
         },
         [&](size_t cluster_num) {
             // Give the offset of the position we chose to the start of the read
             if (cluster_num < left_clusters.size()) {
                 return alignment_1.sequence().begin() - left_clusters[cluster_num]->front().first->begin;
             } else {
                 return alignment_2.sequence().begin() - right_clusters[cluster_num - left_clusters.size()]->front().first->begin;
             }
         },
         node_path_memo, handle_memo);
    
    // Flatten the distance tree to a set of linear spaces, one per tree.
    vector<unordered_map<size_t, int64_t>> linear_spaces = flatten_distance_tree(total_clusters, distance_tree);
    
#ifdef debug_od_clusterer
    for (const auto& strand : linear_spaces) {
        cerr << "strand reconstruction: "  << endl;
        for (const auto& record : strand) {
            if (record.first < left_clusters.size()) {
                cerr << "\t" << record.first << " left: " << record.second << "\t" << left_clusters[record.first]->front().second << endl;
            }
            else {
                cerr << "\t" << record.first - left_clusters.size() << " right: " << record.second << "\t" << right_clusters[record.first - left_clusters.size()]->front().second << endl;
            }

        }
    }
#endif
    
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
            if (sorted_pos[i].second >= left_clusters.size()) {
                continue;
            }
            
            // the interval of linearized coordinates we want to form pairs to
            int64_t coord_interval_start = sorted_pos[i].first + min_inter_cluster_distance;
            int64_t coord_interval_end = sorted_pos[i].first + max_inter_cluster_distance;
            
#ifdef debug_od_clusterer
            cerr << "looking for clusters consistent with cluster that starts with " << left_clusters[sorted_pos[i].second]->front().second << " at relative position " << sorted_pos[i].first << " in coordinate window " << coord_interval_start << ":" << coord_interval_end << endl;
#endif
            
            // move the window bounds forward until it's inside the coordinate interval
            while (window_start < sorted_pos.size() ? sorted_pos[window_start].first < coord_interval_start : false) {
                window_start++;
#ifdef debug_od_clusterer
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
#ifdef debug_od_clusterer
                cerr << "moving window end to relative position " << sorted_pos[window_last - 1].first << endl;
#endif
            }
            
            // add each pair of clusters that's from the two read ends to the return value
            for (size_t j = window_start; j <= window_last; j++) {
                if (sorted_pos[j].second >= left_clusters.size()) {
#ifdef debug_od_clusterer
                    cerr << "adding pair with cluster relative position " << sorted_pos[j].first << " starting with " << right_clusters[sorted_pos[j].second - left_clusters.size()]->front().second << endl;
#endif
                    to_return.emplace_back(make_pair(sorted_pos[i].second,
                                                     sorted_pos[j].second - left_clusters.size()),
                                           sorted_pos[j].first - sorted_pos[i].first);
                }
#ifdef debug_od_clusterer
                else {
                    cerr << "cluster at relative position " << sorted_pos[j].first << " is from the same end, skipping" << endl;
                }
#endif
            }
        }
    }
    
    return to_return;
}

Graph cluster_subgraph(const xg::XG& xg, const Alignment& aln, const vector<vg::MaximalExactMatch>& mems, double expansion) {
    assert(mems.size());
    auto& start_mem = mems.front();
    auto start_pos = make_pos_t(start_mem.nodes.front());
    auto rev_start_pos = reverse(start_pos, xg.node_length(id(start_pos)));
    // Even if the MEM is right up against the start of the read, it may not be
    // part of the best alignment. Make sure to have some padding.
    // TODO: how much padding?
    Graph graph;
    int padding = 1;
    int get_before = padding + (int)(expansion * (int)(start_mem.begin - aln.sequence().begin()));
    if (get_before) {
        graph.MergeFrom(xg.graph_context_id(rev_start_pos, get_before));
    }
    for (int i = 0; i < mems.size(); ++i) {
        auto& mem = mems[i];
        auto pos = make_pos_t(mem.nodes.front());
        int get_after = padding + (i+1 == mems.size() ?
                                   expansion * (int)(aln.sequence().end() - mem.begin)
                                   : expansion * max(mem.length(), (int)(mems[i+1].end - mem.begin)));
        graph.MergeFrom(xg.graph_context_id(pos, get_after));
    }
    sort_by_id_dedup_and_clean(graph);
    return graph;
}

}











