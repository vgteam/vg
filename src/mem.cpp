
#include <string>
#include <algorithm>
#include <cstring>
#include <sstream>

#include "mem.hpp"

namespace vg {
    
using namespace std;

// construct the sequence of the MEM; useful in debugging
string MaximalExactMatch::sequence(void) const {
    return string(begin, end);
}
    
// length of the MEM
int MaximalExactMatch::length(void) const {
    return end - begin;
}

// counts Ns in the MEM
size_t MaximalExactMatch::count_Ns(void) const {
    return std::count(begin, end, 'N');
}

size_t MaximalExactMatch::filter_hits_to(int limit) {
    int keep_every = ceil((double)nodes.size()/(double)limit);
    int i = 0;
    nodes.erase(std::remove_if(nodes.begin(), nodes.end(), [&i,&keep_every](const gcsa::node_type& n) { return i++%keep_every==0; }),
                nodes.end());
    return nodes.size();
}

pair<int64_t, int64_t> mem_min_oriented_distances(const MaximalExactMatch& m1, const MaximalExactMatch& m2) {
    return min_oriented_distances(m1.positions, m2.positions);
}

bool operator==(const MaximalExactMatch& m1, const MaximalExactMatch& m2) {
    return m1.begin == m2.begin && m1.end == m2.end && m1.nodes == m2.nodes;
}

bool operator<(const MaximalExactMatch& m1, const MaximalExactMatch& m2) {
    return m1.begin < m2.begin && m1.end < m2.end && m1.nodes < m2.nodes;
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

vector<string::const_iterator> cluster_cover(const vector<MaximalExactMatch>& cluster) {
    vector<string::const_iterator> seen;
    for (auto& mem : cluster) {
        string::const_iterator c = mem.begin;
        while (c != mem.end) seen.push_back(c++);
    }
    std::sort(seen.begin(), seen.end());
    seen.erase(unique(seen.begin(), seen.end()), seen.end());
    return seen;
}

// rank the clusters by the number of unique read bases they cover
int cluster_coverage(const vector<MaximalExactMatch>& cluster) {
    return cluster_cover(cluster).size();
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

bool clusters_overlap_in_read(const vector<MaximalExactMatch>& cluster1,
                              const vector<MaximalExactMatch>& cluster2) {
    return clusters_overlap_length(cluster1, cluster2) > 0;
}

int clusters_overlap_length(const vector<MaximalExactMatch>& cluster1,
                            const vector<MaximalExactMatch>& cluster2) {
    vector<string::const_iterator> cov1 = cluster_cover(cluster1);
    vector<string::const_iterator> cov2 = cluster_cover(cluster2);
    vector<string::const_iterator> both;
    std::set_intersection(cov1.begin(), cov1.end(),
                          cov2.begin(), cov2.end(),
                          std::back_inserter(both));
    return both.size();
}

vector<pos_t> cluster_nodes(const vector<MaximalExactMatch>& cluster) {
    vector<pos_t> nodes;
    for (auto& mem : cluster) {
        for (auto& node : mem.nodes) {
            nodes.push_back(make_pos_t(node));
            auto& pos = nodes.back();
            get_offset(pos) = 0;
        }
    }
    sort(nodes.begin(), nodes.end());
    return nodes;
}

bool clusters_overlap_in_graph(const vector<MaximalExactMatch>& cluster1,
                               const vector<MaximalExactMatch>& cluster2) {
    vector<pos_t> pos1 = cluster_nodes(cluster1);
    vector<pos_t> pos2 = cluster_nodes(cluster2);
    vector<pos_t> comm;
    set_intersection(pos1.begin(), pos1.end(),
                     pos2.begin(), pos2.end(),
                     std::back_inserter(comm));
    return comm.size() > 0;
}

vector<MaximalExactMatch> translate_mems(const vector<MaximalExactMatch>& mems,
                                         const unordered_map<id_t, pair<id_t, bool> >& trans) {
    // invert the translation
    unordered_map<id_t, vector<pair<id_t, bool>>> inv_trans;
    for (auto& t : trans) {
        id_t new_id = t.first;
        id_t old_id = t.second.first;
        bool flip = t.second.second;
        inv_trans[old_id].push_back(make_pair(new_id, flip));
    }
    vector<MaximalExactMatch> trans_mems;
    for (auto& mem : mems) {
        auto new_mem = mem;
        new_mem.nodes.clear();
        for (auto& node : mem.nodes) {
            id_t node_id = gcsa::Node::id(node);
            for (auto& p : inv_trans[node_id]) {
                id_t new_id = p.first;
                size_t node_offset = gcsa::Node::offset(node);
                bool new_rc = (p.second ? !gcsa::Node::rc(node) : gcsa::Node::rc(node));
                new_mem.nodes.push_back(gcsa::Node::encode(new_id, node_offset, new_rc));
            }
        }
        trans_mems.push_back(new_mem);
    }
    return trans_mems;
}

}
