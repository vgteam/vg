
#include <string>
#include <algorithm>
#include <cstring>
#include <sstream>

#include "mem.hpp"

namespace vg {
    
    
//size_t OrientedDistanceClusterer::PRUNE_COUNTER = 0;
//size_t OrientedDistanceClusterer::CLUSTER_TOTAL = 0;
//size_t OrientedDistanceClusterer::MEM_FILTER_COUNTER = 0;
//size_t OrientedDistanceClusterer::MEM_TOTAL = 0;
    
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

int64_t mem_min_distance(const MaximalExactMatch& m1, const MaximalExactMatch& m2) {
    int64_t distance = std::numeric_limits<int64_t>::max();
    for (auto& seq : m1.positions) {
        auto& name = seq.first;
        auto f = m2.positions.find(name);
        if (f != m2.positions.end()) {
            auto& pos1 = seq.second;
            auto& pos2 = f->second;
            for (auto& p1 : pos1) {
                for (auto& p2 : pos2) {
                    distance = min(abs((int64_t)p1 - (int64_t)p2), distance);
                }
            }
        }
    }
    return distance;
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

bool clusters_overlap_in_read(const vector<MaximalExactMatch>& cluster1,
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
}











