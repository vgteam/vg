
#include <string>
#include <algorithm>

#include "mem.hpp"

namespace vg {

using namespace std;

// construct the sequence of the MEM; useful in debugging
string MaximalExactMatch::sequence(void) const {
    string seq; //seq.resize(end-begin);
    string::const_iterator c = begin;
    while (c != end) seq += *c++;
    return seq;
}
    
// length of the MEM
int MaximalExactMatch::length(void) const {
    return end - begin;
}

// counts Ns in the MEM
size_t MaximalExactMatch::count_Ns(void) const {
    return std::count(begin, end, 'N');
}


ostream& operator<<(ostream& out, const MaximalExactMatch& mem) {
    size_t len = mem.begin - mem.end;
    out << mem.sequence() << ":";
    for (auto& node : mem.nodes) {
        id_t id = gcsa::Node::id(node);
        size_t offset = gcsa::Node::offset(node);
        bool is_rev = gcsa::Node::rc(node);
        out << id << (is_rev ? "-" : "+") << ":" << offset << ",";
    }
    return out;
}

bool operator==(const MaximalExactMatch& m1, const MaximalExactMatch& m2) {
    return m1.begin == m2.begin && m1.end == m2.end && m1.nodes == m2.nodes;
}

bool operator<(const MaximalExactMatch& m1, const MaximalExactMatch& m2) {
    return m1.begin < m2.begin && m1.end < m2.end && m1.nodes < m2.nodes;
}

const string mems_to_json(const vector<MaximalExactMatch>& mems) {
    stringstream s;
    s << "[";
    size_t j = 0;
    for (auto& mem : mems) {
        s << "[\"";
        s << mem.sequence();
        s << "\",[";
        size_t i = 0;
        for (auto& node : mem.nodes) {
            s << "\"" << gcsa::Node::decode(node) << "\"";
            if (++i < mem.nodes.size()) s << ",";
        }
        s << "]]";
        if (++j < mems.size()) s << ",";
    }
    s << "]";
    return s.str();
}

// rank the clusters by the number of unique read bases they cover
int cluster_coverage(const vector<MaximalExactMatch>& cluster) {
    set<string::const_iterator> seen;
    for (auto& mem : cluster) {
        string::const_iterator c = mem.begin;
        while (c != mem.end) seen.insert(c++);
    }
    return seen.size();
}

bool mems_overlap(const MaximalExactMatch& mem1,
                  const MaximalExactMatch& mem2) {
    // we overlap if we are not completely separated
    return mem1.fragment == mem2.fragment
        && !(mem1.end <= mem2.begin
             || mem2.end <= mem1.begin);
}

int mems_overlap_length(const MaximalExactMatch& mem1,
                        const MaximalExactMatch& mem2) {
    if (!mems_overlap(mem1, mem2)) {
        return 0;
    } else {
        if (mem1.begin < mem2.begin) {
            if (mem1.end < mem2.end) {
                return mem1.end - mem2.begin;
            } else {
                return mem1.end - mem1.begin;
            }
        } else {
            if (mem2.end < mem1.end) {
                return mem2.end - mem1.begin;
            } else {
                return mem2.end - mem2.begin;
            }
        }
    }
}

bool clusters_overlap(const vector<MaximalExactMatch>& cluster1,
                      const vector<MaximalExactMatch>& cluster2) {
    for (auto& mem1 : cluster1) {
        for (auto& mem2 : cluster2) {
            if (mems_overlap(mem1, mem2)) {
                return true;
            }
        }
    }
    return false;
}



MEMChainModel::MEMChainModel(
    const vector<size_t>& aln_lengths,
    const vector<vector<MaximalExactMatch> >& matches,
    const function<int(pos_t)>& approx_position,
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
                MEMChainModelVertex m;
                m.mem = mem;
                m.weight = mem.length();
                m.prev = nullptr;
                m.score = 0;
                m.approx_position = approx_position(make_pos_t(node));
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
        approx_positions[v->approx_position].push_back(v);
    }
    // sort the vertexes at each approx position by their matches and trim
    for (auto& pos : approx_positions) {
        std::sort(pos.second.begin(), pos.second.end(), [](const vector<MEMChainModelVertex>::iterator& v1,
                                                           const vector<MEMChainModelVertex>::iterator& v2) {
                      return v1->mem.length() > v2->mem.length();
                  });
        if (pos.second.size() > position_depth) {
            for (int i = position_depth; i < pos.second.size(); ++i) {
                redundant_vertexes.insert(pos.second[i]);
            }
        }
        pos.second.resize(min(pos.second.size(), (size_t)position_depth));
    }
    // for each vertex merge if we go equivalently forward in the positional space and forward in the read to the next position
    // scan forward
    for (map<int, vector<vector<MEMChainModelVertex>::iterator> >::iterator p = approx_positions.begin();
         p != approx_positions.end(); ++p) {
        for (auto& v1 : p->second) {
            if (redundant_vertexes.count(v1)) continue;
            auto q = p;
            while (++q != approx_positions.end() && abs(p->first - q->first) < band_width) {
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
    // scan reverse
    for (map<int, vector<vector<MEMChainModelVertex>::iterator> >::reverse_iterator p = approx_positions.rbegin();
         p != approx_positions.rend(); ++p) {
        for (auto& v1 : p->second) {
            if (redundant_vertexes.count(v1)) continue;
            auto q = p;
            while (++q != approx_positions.rend() && abs(p->first - q->first) < band_width) {
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
    // now build up the model using the positional bandwidth
    for (map<int, vector<vector<MEMChainModelVertex>::iterator> >::iterator p = approx_positions.begin();
         p != approx_positions.end(); ++p) {
        // look bandwidth before and bandwidth after in the approx positions
        // after
        for (auto& v1 : p->second) {
            if (redundant_vertexes.count(v1)) continue;
            auto q = p;
            while (++q != approx_positions.end() && abs(p->first - q->first) < band_width) {
                for (auto& v2 : q->second) {
                    if (redundant_vertexes.count(v2)) continue;
                    // if this is an allowable transition, run the weighting function on it
                    if (v1->next_cost.size() < max_connections
                        && v2->prev_cost.size() < max_connections) {
                        if (v1->mem.fragment < v2->mem.fragment
                            || v1->mem.fragment == v2->mem.fragment && v1->mem.begin < v2->mem.begin) {
                            double weight = transition_weight(v1->mem, v2->mem);
                            if (weight > -std::numeric_limits<double>::max()) {
                                v1->next_cost.push_back(make_pair(&*v2, weight));
                                v2->prev_cost.push_back(make_pair(&*v1, weight));
                            }
                        } else if (v1->mem.fragment > v2->mem.fragment
                                   || v1->mem.fragment == v2->mem.fragment && v1->mem.begin > v2->mem.begin) {
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

OrientedDistanceClusterer::OrientedDistanceClusterer(const Alignment& alignment,
                                       const vector<MaximalExactMatch>& mems,
                                       const QualAdjAligner& aligner,
                                       xg::XG* xgindex,
                                       size_t max_expected_dist_approx_error) : aligner(aligner) {
    
    // the maximum graph distances to the right and left sides detectable from each node
    vector<pair<size_t, size_t>> maximum_detectable_gaps;
    
    // there generally will be at least as many nodes as MEMs, so we can speed up the reallocation
    nodes.reserve(mems.size());
    maximum_detectable_gaps.reserve(mems.size());
    
    for (const MaximalExactMatch& mem : mems) {
        
        // calculate the longest gaps we could detect to the left and right of this MEM
        pair<size_t, size_t> max_gaps(aligner.longest_detectable_gap(alignment, mem.begin),
                                      aligner.longest_detectable_gap(alignment, mem.end));
        
        
        int32_t mem_score = aligner.score_exact_match(mem.begin, mem.end,
                                                      alignment.quality().begin() + (mem.begin - alignment.sequence().begin()));
        
        for (gcsa::node_type mem_hit : mem.nodes) {
            nodes.emplace_back(mem, make_pos_t(mem_hit), mem_score);
            maximum_detectable_gaps.push_back(max_gaps);
        }
    }
    
    // now we use the distance approximation to cluster the MEM hits according to the strand
    // they fall on using the oriented distance estimation function
    
    // for recording the distance of any pair that we check with a finite distance
    unordered_map<pair<size_t, size_t>, int64_t> recorded_finite_dists;
    // for recording the number of times elements of a strand cluster have been compared
    // and found an infinite distance
    map<pair<size_t, size_t>, size_t> num_infinite_dists;
    
    // deterministically generate pseudo-shuffled pairs in constant time per pair, adapted from
    // https://stackoverflow.com/questions/1866684/algorithm-to-print-out-a-shuffled-list-in-place-and-with-o1-memory
    size_t num_pairs = nodes.size() * nodes.size();
    assert(num_pairs < 4294967296); // 2^32, else range_max can't get large enough
    size_t permutation_idx = 0;
    size_t range_max = 1;
    while (range_max < num_pairs) {
        range_max *= 2;
    }
    auto get_next_pair = [&]() {
        pair<size_t, size_t> to_return;
        do {
            size_t permuted = ((permutation_idx ^ (range_max - 1)) ^ (permutation_idx << 6) + 0x9e3779b9) & (range_max - 1);
            to_return.first = permuted / nodes.size();
            to_return.second = permuted % nodes.size();
            permutation_idx++;
        } while (to_return.first >= to_return.second || to_return.first >= nodes.size());
        return to_return;
    };
    
    // we use a union find to keep track of which MEMs have been identified as being on the same strand
    UnionFind union_find(nodes.size());
    
    size_t num_possible_merges_remaining = (nodes.size() * (nodes.size() - 1)) / 2;
    size_t pairs_checked = 0;
    
    // a simulated annealing parameter loosely inspired by the cutoff for an Erdos Renyi random graph
    // to be connected with probability approaching 1
    size_t current_max_num_probes = 2 * ((size_t) ceil(log(nodes.size())));
    
    while (num_possible_merges_remaining > 0 && permutation_idx < range_max && current_max_num_probes > 0) {
        // slowly lower the number of distances we need to check before we believe that two clusters are on
        // separate strands
#ifdef debug_multipath_mapper
        size_t direct_merges_remaining = 0;
        vector<vector<size_t>> groups = union_find.all_groups();
        for (size_t i = 1; i < groups.size(); i++) {
            for (size_t j = 0; j < i; j++) {
                size_t strand_1 = union_find.find_group(groups[i].front());
                size_t strand_2 = union_find.find_group(groups[j].front());
                
                if (num_infinite_dists.count(make_pair(strand_1, strand_2))) {
                    if (num_infinite_dists[make_pair(strand_1, strand_2)] >= current_max_num_probes) {
                        continue;
                    }
                }
                direct_merges_remaining += groups[i].size() * groups[j].size();
            }
        }
        cerr << "at permutation index " << permutation_idx << " with directly calculated merges " << direct_merges_remaining << " and maintained merges " << num_possible_merges_remaining << endl;
#endif
        
        if (pairs_checked % nodes.size() == 0 && pairs_checked != 0) {
            current_max_num_probes--;
#ifdef debug_multipath_mapper
            cerr << "reducing the max number of probes to " << current_max_num_probes << endl;
#endif
            for (const pair<pair<size_t, size_t>, size_t>& inf_dist_record : num_infinite_dists) {
                // break symmetry so we don't repeat the operation twice
                if (inf_dist_record.first.first < inf_dist_record.first.second && inf_dist_record.second == current_max_num_probes) {
                    // this merge just fell below the new maximum number of distance probes
                    size_t strand_size_1 = union_find.group_size(inf_dist_record.first.first);
                    size_t strand_size_2 = union_find.group_size(inf_dist_record.first.second);
                    num_possible_merges_remaining -= strand_size_1 * strand_size_2;
#ifdef debug_multipath_mapper
                    cerr << "after reduction, the total number of probes between strand " << inf_dist_record.first.first << " and " << inf_dist_record.first.second <<  " is above max, reducing possible merges by " << strand_size_1 * strand_size_2 << " to " << num_possible_merges_remaining << endl;
#endif
                }
            }
        }
        
        
        pair<size_t, size_t> node_pair = get_next_pair();
    
        pairs_checked++;
        
        size_t strand_1 = union_find.find_group(node_pair.first);
        size_t strand_2 = union_find.find_group(node_pair.second);
        
#ifdef debug_multipath_mapper
        cerr << "checking MEMs " << node_pair.first << " and " << node_pair.second << endl;
#endif

        if (strand_1 == strand_2) {
            // these are already identified as on the same strand, don't need to do it again
#ifdef debug_multipath_mapper
            cerr << "already on same strand" << endl;
#endif
            continue;
        }
        
        auto num_failed_probes = num_infinite_dists.find(make_pair(strand_1, strand_2));
        if (num_failed_probes == num_infinite_dists.end() ? false : num_failed_probes->second >= current_max_num_probes) {
            // we've already checked multiple distances between these strand clusters and
            // none have returned a finite distance, so we conclude that they are in fact
            // on separate clusters and decline to check any more distances
#ifdef debug_multipath_mapper
            cerr << "already have checked distance above maximum number of probes" << endl;
#endif
            continue;
        }
        
        pos_t& pos_1 = nodes[node_pair.first].start_pos;
        pos_t& pos_2 = nodes[node_pair.second].start_pos;
        
        int64_t oriented_dist = xgindex->closest_shared_path_oriented_distance(id(pos_1), offset(pos_1), is_rev(pos_1),
                                                                               id(pos_2), offset(pos_2), is_rev(pos_2));
        
#ifdef debug_multipath_mapper
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
                size_t strand_size_1 = union_find.group_size(strand_1);
                size_t strand_size_2 = union_find.group_size(strand_2);
                
                num_possible_merges_remaining -= strand_size_1 * strand_size_2;
                
#ifdef debug_multipath_mapper
                cerr << "number of probes " << num_failed_probes->second << " crossed max threshold of " << current_max_num_probes << ", reducing possible merges by " << strand_size_1 * strand_size_2 << " to " << num_possible_merges_remaining << endl;
#endif
            }
        }
        else {
            // the distance is finite, so merge the strand clusters
            
            recorded_finite_dists[node_pair] = oriented_dist;
            
            size_t strand_size_1 = union_find.group_size(strand_1);
            size_t strand_size_2 = union_find.group_size(strand_2);
            
            union_find.union_groups(node_pair.first, node_pair.second);
            
            // remove these from the pool of remaining merges
            num_possible_merges_remaining -= strand_size_1 * strand_size_2;
            
            size_t strand_retaining = union_find.find_group(node_pair.first);
            size_t strand_removing = strand_retaining == strand_1 ? strand_2 : strand_1;
            
#ifdef debug_multipath_mapper
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
                        num_possible_merges_remaining -= (strand_retaining == strand_1 ? strand_size_2 : strand_size_1) * union_find.group_size(removing_iter->first.second);
                        
#ifdef debug_multipath_mapper
                        cerr << "after merge, the total number of probes against strand " << removing_iter->first.second << " increased to " << retaining_iter->second << ", above current max of " << current_max_num_probes << ", but the retaining strand is already blocked, reducing possible merges by " << (strand_retaining == strand_1 ? strand_size_2 : strand_size_1) * union_find.group_size(removing_iter->first.second) << " to " << num_possible_merges_remaining << endl;
#endif
                    }
                    else if (removing_already_blocked && !retaining_already_blocked) {
                        num_possible_merges_remaining -= (strand_retaining == strand_1 ? strand_size_1 : strand_size_2) * union_find.group_size(removing_iter->first.second);
                        
#ifdef debug_multipath_mapper
                        cerr << "after merge, the total number of probes against strand " << removing_iter->first.second << " increased to " << retaining_iter->second << ", above current max of " << current_max_num_probes << ", but the removing strand is already blocked, reducing possible merges by " << (strand_retaining == strand_1 ? strand_size_1 : strand_size_2) * union_find.group_size(removing_iter->first.second) << " to " << num_possible_merges_remaining << endl;
#endif
                    }
                    else if (!retaining_already_blocked && !removing_already_blocked && retaining_iter->second >= current_max_num_probes) {
                        num_possible_merges_remaining -= (strand_size_1 + strand_size_2) * union_find.group_size(removing_iter->first.second);
                        
#ifdef debug_multipath_mapper
                        cerr << "after merge, the total number of probes against strand " << removing_iter->first.second << " increased to " << retaining_iter->second << ", above current max of " << current_max_num_probes << ", reducing possible merges by " << (strand_size_1 + strand_size_2) * union_find.group_size(removing_iter->first.second) << " to " << num_possible_merges_remaining << endl;
#endif
                        
                    }
                    removing_iter++;
                    retaining_iter++;
                }
                else if (removing_iter->first.second < retaining_iter->first.second) {
                    // the strand being removed has probes against this strand cluster, but the strand being
                    // retained does not, mark this and save it for later so that we don't invalidate the range
                    //
                    unseen_comparisons.emplace_back(removing_iter->first.second, removing_iter->second);
                    removing_iter++;
                }
                else {
                    // the strand being retained has probes against this strand cluster, but the strand being
                    // removed does not, check if we need to add the removing strand to the remaining merges
                    // counter
                    if (retaining_iter->second >= current_max_num_probes) {
                        num_possible_merges_remaining -= (strand_retaining == strand_1 ? strand_size_2 : strand_size_1) * union_find.group_size(retaining_iter->first.second);
                        
#ifdef debug_multipath_mapper
                        cerr << "after merge, the total number of probes against strand " << retaining_iter->first.second << " increased to " << retaining_iter->second << ", above current max of " << current_max_num_probes << ", but the retaining strand is already blocked, reducing possible merges by " << (strand_retaining == strand_1 ? strand_size_2 : strand_size_1) * union_find.group_size(retaining_iter->first.second) << " to " << num_possible_merges_remaining << endl;
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
                    num_possible_merges_remaining -= (strand_retaining == strand_1 ? strand_size_2 : strand_size_1) * union_find.group_size(retaining_iter->first.second);
                    
#ifdef debug_multipath_mapper
                    cerr << "after merge, the total number of probes against strand " << retaining_iter->first.second << " increased to " << retaining_iter->second << ", above current max of " << current_max_num_probes << ", but the retaining strand is already blocked, reducing possible merges by " << (strand_retaining == strand_1 ? strand_size_2 : strand_size_1) * union_find.group_size(retaining_iter->first.second) << " to " << num_possible_merges_remaining << endl;
#endif
                }
                retaining_iter++;
            }
            
            // add the probes between the removing strands and clusters that had never been compared to the retaining strand
            for (const pair<size_t, size_t>& unseen_comparison : unseen_comparisons) {
                num_infinite_dists[make_pair(unseen_comparison.first, strand_retaining)] = unseen_comparison.second;
                num_infinite_dists[make_pair(strand_retaining, unseen_comparison.first)] = unseen_comparison.second;
                
                if (unseen_comparison.second >= current_max_num_probes) {
                    num_possible_merges_remaining -= (strand_retaining == strand_1 ? strand_size_2 : strand_size_1) * union_find.group_size(unseen_comparison.first);
                    
#ifdef debug_multipath_mapper
                    cerr << "after merge, the total number of probes against strand " << unseen_comparison.first << " increased to " << unseen_comparison.second << ", above current max of " << current_max_num_probes << ", but the removing strand is already blocked, reducing possible merges by " << (strand_retaining == strand_1 ? strand_size_2 : strand_size_1) * union_find.group_size(unseen_comparison.first) << " to " << num_possible_merges_remaining << endl;
#endif
                }
            }
            
            // find the range containing the records with the removing strand again (it may have changed since we
            // altered the map)
            removing_iter = num_infinite_dists.lower_bound(make_pair(strand_removing, 0));
            removing_end = num_infinite_dists.upper_bound(make_pair(strand_removing, numeric_limits<size_t>::max()));
            if (removing_iter != removing_end) {
                // move the end so that it is an inclusive range (prevents the "end" value from changing if we insert
                // something between the last and the past-the-last positions)
                removing_end--;
                
                // erase the range
                if (removing_iter == removing_end) {
                    num_infinite_dists.erase(make_pair(removing_iter->first.second, removing_iter->first.first));
                    num_infinite_dists.erase(removing_iter);
                }
                else {
                    // erase the previous position on each iteration so that we don't invalidate the iterator before
                    // we use it to move to the next position
                    auto removing_iter_prev = removing_iter;
                    removing_iter++;
                    while (removing_iter != removing_end) {
                        num_infinite_dists.erase(make_pair(removing_iter_prev->first.second, removing_iter_prev->first.first));
                        num_infinite_dists.erase(removing_iter_prev);
                        removing_iter_prev = removing_iter;
                        removing_iter++;
                    }
                    num_infinite_dists.erase(make_pair(removing_iter_prev->first.second, removing_iter_prev->first.first));
                    num_infinite_dists.erase(removing_iter_prev);
                    num_infinite_dists.erase(make_pair(removing_iter->first.second, removing_iter->first.first));
                    num_infinite_dists.erase(removing_iter);
                }
            }
        }
    }
    
#ifdef debug_multipath_mapper
    cerr << "constructing strand distance tree" << endl;
#endif
    
    // build the graph of relative distances in adjacency list representation
    // by construction each strand cluster will be an undirected, unrooted tree
    vector<vector<size_t>> strand_distance_tree(nodes.size());
    for (const auto& dist_record : recorded_finite_dists) {
        strand_distance_tree[dist_record.first.first].push_back(dist_record.first.second);
        strand_distance_tree[dist_record.first.second].push_back(dist_record.first.first);
    }
    
    // now approximate the relative positions along the strand by traversing each tree and
    // treating the distances we estimated as transitive
    vector<unordered_map<size_t, int64_t>> strand_relative_position;
    vector<bool> processed(nodes.size(), false);
    for (size_t i = 0; i < nodes.size(); i++) {
        if (processed[i]) {
            continue;
        }
        
#ifdef debug_multipath_mapper
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
                               recorded_finite_dists[make_pair(curr, next)] :
                               -recorded_finite_dists[make_pair(next, curr)];
                
                // find the position relative to the previous node we just traversed
                relative_pos[next] = curr_pos + dist;
                processed[next] = true;
                
                queue.push_back(next);
            }
        }
#ifdef debug_multipath_mapper
        cerr << "strand reconstruction: "  << endl;
        for (const auto& record : relative_pos) {
            cerr << "\t" << record.first << ": " << record.second << endl;
        }
#endif
    }
    
    
    // now we use the strand clusters and the estimated distances to make the DAG for the
    // approximate MEM alignment
    
    int64_t allowance = max_expected_dist_approx_error;
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
            
            // the limits of how far away we might detect edges to add to the clustering graph
            int64_t target_low_pos = strand_pos - allowance;
            int64_t target_hi_pos = strand_pos + pivot_length + maximum_detectable_gaps[pivot_idx].second + allowance;
            
            // move the lower boundary of the search interval to the lowest value inside the
            // the target interval
            while (sorted_pos[low].first < target_low_pos) {
                low++;
            }
            
            // move the upper boundary of the search interval to the highest value inside the
            // the target interval (this one can move in either direction because the maximum
            // detectable gap changes)
            if (sorted_pos[hi].first > target_hi_pos) {
                while (sorted_pos[hi].first > target_hi_pos) {
                    hi--;
                }
            }
            else {
                while (hi == last_idx ? false : sorted_pos[hi + 1].first <= target_hi_pos) {
                    hi++;
                }
            }
            
#ifdef debug_multipath_mapper
            cerr << "checking for possible edges from " << sorted_pos[i].second << " to MEMs between " << sorted_pos[low].first << "(" << sorted_pos[low].second << ") and " << sorted_pos[hi].first << "(" << sorted_pos[hi].second << ")" << endl;
#endif
            
            for (int64_t j = low; j <= hi; j++) {
                int64_t next_idx = sorted_pos[j].second;
                ODNode& next = nodes[next_idx];
                
                if (next.mem->begin <= pivot.mem->begin || next.mem->end <= pivot.mem->end) {
                    // the MEMs cannot be colinear along the read (also filters out j == i)
                    continue;
                }
                
#ifdef debug_multipath_mapper
                cerr << "adding edge to MEM " << sorted_pos[j].first << "(" << sorted_pos[j].second << ")" << endl;
#endif
                
                // the length of the sequence in between the MEMs (can be negative if they overlap)
                int64_t between_length = next.mem->begin - pivot.mem->end;
                // the estimated distance between the end of the pivot and the start of the next MEM in the graph
                int64_t graph_dist = max<int64_t>(0, sorted_pos[j].first - strand_pos - pivot_length);
                // the discrepancy between the graph distance and the read distance
                int64_t gap_length = abs(graph_dist - between_length);
                
                if (gap_length > maximum_detectable_gaps[next_idx].first + allowance) {
                    // the gap between the MEMs is too long to be believable from the next node
                    continue;
                }
                
                int32_t edge_score;
                if (between_length < 0) {
                    // the MEMs overlap, but this can occur in some insertions and deletions
                    // because the SMEM algorithm is "greedy" in taking up as much of the read
                    // as possible
                    // we can check if this happened directly, but it's expensive
                    // so for now we just give it the benefit of the doubt but adjust the edge
                    // score so that the matches don't get double counted
                    
                    int64_t extra_dist = max<int64_t>(0, gap_length);
                    
                    edge_score = -aligner.match * between_length
                                 + (extra_dist ? -(extra_dist - 1) * aligner.gap_extension - aligner.gap_open : 0);
                }
                else if (between_length > graph_dist) {
                    // the read length in between the MEMs is longer than the distance, suggesting a read insert
                    edge_score = -aligner.mismatch * graph_dist - (gap_length - 1) * aligner.gap_extension
                                 - aligner.gap_open;
                }
                else if (between_length < graph_dist) {
                    // the read length in between the MEMs is shorter than the distance, suggesting a read deletion
                    edge_score = -aligner.mismatch * between_length - (gap_length - 1) * aligner.gap_extension
                                 - aligner.gap_open;
                }
                else {
                    // the read length in between the MEMs is the same as the distance, suggesting a pure mismatch
                    edge_score = -aligner.mismatch * between_length;
                }
                
                // add the edges in
                pivot.edges_from.emplace_back(next_idx, edge_score);
                next.edges_to.emplace_back(pivot_idx, edge_score);
            }
        }
    }
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
    
    // as in local alignment, minimum score is the score of node itself
    for (size_t i = 0; i < nodes.size(); i++) {
        nodes[i].dp_score = nodes[i].score;
    }
    
#ifdef debug_multipath_mapper
    cerr << "computing topological order for clustering DP" << endl;
#endif
    
    vector<size_t> order;
    topological_order(order);
    
    for (size_t i : order) {
        ODNode& node = nodes[i];
#ifdef debug_multipath_mapper
        cerr << "at node " << i << " with DP score " << node.dp_score << " and node score " << node.score << endl;
#endif
        // for each edge out of this node
        for (ODEdge& edge : node.edges_from) {
            
            // check if the path through the node out of this edge increase score of target node
            ODNode& target_node = nodes[edge.to_idx];
            int32_t extend_score = node.dp_score + edge.weight + target_node.score;
            if (extend_score > target_node.dp_score) {
#ifdef debug_multipath_mapper
                cerr << "extending DP to node " << edge.to_idx << " with score " << extend_score << endl;
#endif
                target_node.dp_score = extend_score;
            }
        }
    }
}

vector<vector<pair<const MaximalExactMatch*, pos_t>>> OrientedDistanceClusterer::clusters(int32_t max_qual_score) {
    
    vector<vector<pair<const MaximalExactMatch*, pos_t>>> to_return;
    if (nodes.size() == 0) {
        // this should only happen if we have filtered out all MEMs, so there are none to cluster
        return to_return;
    }
    
#ifdef debug_multipath_mapper
    cerr << "performing approximate DP across MEMs" << endl;
#endif
    perform_dp();
    
#ifdef debug_multipath_mapper
    cerr << "finding top tracebacks within connected components" << endl;
#endif
    // find the weakly connected components, which should correspond to mappings
    vector<vector<size_t>> components;
    connected_components(components);
    
    // find the node with the highest DP score in each connected component
    // each record is a pair of (score, node index)
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
    
    std::make_heap(component_traceback_ends.begin(), component_traceback_ends.end());
    
    int32_t top_score = component_traceback_ends.front().first;
    
    while (!component_traceback_ends.empty()) {
        // get the next highest scoring traceback end
        auto traceback_end = component_traceback_ends.front();
        std::pop_heap(component_traceback_ends.begin(), component_traceback_ends.end());
        component_traceback_ends.pop_back();
        
        // get the index of the node
        size_t trace_idx = traceback_end.second;
        
#ifdef debug_multipath_mapper
        cerr << "checking traceback of component starting at " << traceback_end.second << endl;
#endif
        // if this cluster does not look like it even affect the mapping quality of the top scoring
        // cluster, don't bother forming it
        // TODO: this approximation could break down sometimes, need to look into it
        // TODO: is there a way to make the aligner do this? I don't like having this functionality outside of it
        if (4.3429448190325183 * aligner.log_base * (top_score - traceback_end.first) > max_qual_score ) {
#ifdef debug_multipath_mapper
            cerr << "skipping rest of components on account of low score of " << traceback_end.first << " compared to max score " << top_score << endl;
#endif
            break;
        }
        
        // traceback until hitting a node that has its own score (indicates beginning of a local alignment)
        vector<size_t> trace{trace_idx};
        while (nodes[trace_idx].dp_score > nodes[trace_idx].score) {
            int32_t target_source_score = nodes[trace_idx].dp_score - nodes[trace_idx].score;
            cerr << "trace idx " << trace_idx << " DP score " << nodes[trace_idx].dp_score << " node score " << nodes[trace_idx].score << " source score " << target_source_score << endl;
            for (ODEdge& edge : nodes[trace_idx].edges_to) {
                cerr << "\tscore along edge to " << edge.to_idx << " is " <<  nodes[edge.to_idx].dp_score + edge.weight << endl;
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
    
    return to_return;
}


}
