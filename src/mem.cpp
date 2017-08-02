
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

}
