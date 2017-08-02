#ifndef VG_MEM_HPP_INCLUDED
#define VG_MEM_HPP_INCLUDED

#include "gcsa/gcsa.h"
#include "gcsa/lcp.h"
#include "position.hpp"
#include "gssw_aligner.hpp"

#include <functional>
#include <string>
#include <vector>
#include <map>


/**
 * \file mem.hpp
 *
 * Maximal Exact Matches (MEMs) and the chaining and clustering tools to work
 * with them.
 */

namespace vg {

using namespace std;

class MaximalExactMatch {

public:

    //const string* source;
    string::const_iterator begin;
    string::const_iterator end;
    gcsa::range_type range;
    size_t match_count;
    int fragment;
    bool primary; // if not a sub-MEM
    std::vector<gcsa::node_type> nodes;
    map<string, vector<size_t> > positions;
    
    MaximalExactMatch(string::const_iterator b,
                      string::const_iterator e,
                      gcsa::range_type r,
                      size_t m = 0)
        : begin(b), end(e), range(r), match_count(m) { }

    // construct the sequence of the MEM; useful in debugging
    string sequence(void) const;
    // get the length of the MEM
    int length(void) const;
    // uses a callback to fill out the MEM positions
    void fill_positions(const function<map<string, vector<size_t>>(gcsa::node_type)>& node_positions_in_paths);
    // tells if the MEM contains an N
    size_t count_Ns(void) const;

    friend bool operator==(const MaximalExactMatch& m1, const MaximalExactMatch& m2);
    friend bool operator<(const MaximalExactMatch& m1, const MaximalExactMatch& m2);
    friend ostream& operator<<(ostream& out, const MaximalExactMatch& m);

    MaximalExactMatch(void) = default;                                      // Copy constructor
    MaximalExactMatch(const MaximalExactMatch&) = default;               // Copy constructor
    MaximalExactMatch(MaximalExactMatch&&) = default;                    // Move constructor
    MaximalExactMatch& operator=(const MaximalExactMatch&) & = default;  // MaximalExactMatchopy assignment operator
    MaximalExactMatch& operator=(MaximalExactMatch&&) & = default;       // Move assignment operator
    //virtual ~MaximalExactMatch() { }                     // Destructor
};

const string mems_to_json(const vector<MaximalExactMatch>& mems);

// helper for computing the number of bases in the query covered by a cluster
int cluster_coverage(const vector<MaximalExactMatch>& cluster);
// helper to tell if mems are ovelapping
bool mems_overlap(const MaximalExactMatch& mem1,
                  const MaximalExactMatch& mem2);
// distance of overlap, or 0 if there is no overlap
int mems_overlap_length(const MaximalExactMatch& mem1,
                        const MaximalExactMatch& mem2);
// helper to tell if clusters have any overlap
bool clusters_overlap(const vector<MaximalExactMatch>& cluster1,
                      const vector<MaximalExactMatch>& cluster2);

class MEMChainModelVertex {
public:
    MaximalExactMatch mem;
    vector<pair<MEMChainModelVertex*, double> > next_cost; // for forward
    vector<pair<MEMChainModelVertex*, double> > prev_cost; // for backward
    double weight;
    double score;
    int approx_position;
    MEMChainModelVertex* prev;
    MEMChainModelVertex(void) = default;                                      // Copy constructor
    MEMChainModelVertex(const MEMChainModelVertex&) = default;               // Copy constructor
    MEMChainModelVertex(MEMChainModelVertex&&) = default;                    // Move constructor
    MEMChainModelVertex& operator=(const MEMChainModelVertex&) & = default;  // MEMChainModelVertexopy assignment operator
    MEMChainModelVertex& operator=(MEMChainModelVertex&&) & = default;       // Move assignment operator
    virtual ~MEMChainModelVertex() { }                     // Destructor
};

class MEMChainModel {
public:
    vector<MEMChainModelVertex> model;
    map<int, vector<vector<MEMChainModelVertex>::iterator> > approx_positions;
    set<vector<MEMChainModelVertex>::iterator> redundant_vertexes;
    MEMChainModel(
        const vector<size_t>& aln_lengths,
        const vector<vector<MaximalExactMatch> >& matches,
        const function<int(pos_t)>& approx_position,
        const function<double(const MaximalExactMatch&, const MaximalExactMatch&)>& transition_weight,
        int band_width = 10,
        int position_depth = 1,
        int max_connections = 10);
    void score(const set<MEMChainModelVertex*>& exclude);
    MEMChainModelVertex* max_vertex(void);
    vector<vector<MaximalExactMatch> > traceback(int alt_alns, bool paired, bool debug);
    void display(ostream& out);
    void clear_scores(void);
};

class OrientedDistanceClusterer {
public:
    OrientedDistanceClusterer(const Alignment& alignment,
                              const vector<MaximalExactMatch>& mems,
                              const QualAdjAligner& aligner,
                              xg::XG* xgindex,
                              size_t max_expected_dist_approx_error = 8);
    
    /// Returns a vector of clusters. Each cluster is represented a vector of MEM hits. Each hit
    /// contains a pointer to the original MEM and the position of that particular hit in the graph.
    vector<vector<pair<const MaximalExactMatch*, pos_t>>> clusters(int32_t max_qual_score = 60);
    
private:
    class ODNode;
    class ODEdge;
    struct DPScoreComparator;
    
    /// Fills input vectors with indices of source and sink nodes
    void identify_sources_and_sinks(vector<size_t>& sources_out, vector<size_t>& sinks_out);
    
    /// Identify weakly connected components in the graph
    void connected_components(vector<vector<size_t>>& components_out);
    
    /// Fills the input vector with the indices of a topological sort
    void topological_order(vector<size_t>& order_out);
    
    /// Perform dynamic programming and store scores in nodes
    void perform_dp();
    
    vector<ODNode> nodes;
    
    const QualAdjAligner& aligner;
};

class OrientedDistanceClusterer::ODNode {
public:
    ODNode(const MaximalExactMatch& mem, pos_t start_pos, int32_t score) :
            mem(&mem), start_pos(start_pos), score(score) {}
    ODNode() = default;
    ~ODNode() = default;
    
    MaximalExactMatch const* mem;
    
    /// Position of GCSA hit in the graph
    pos_t start_pos;
    /// Score of the exact match this node represents
    int32_t score;
    
    /// Score during dynamic programming
    int32_t dp_score;
    
    /// Edges from this node that are colinear with the read
    vector<ODEdge> edges_from;
    
    /// Edges to this node that are colinear with the read
    vector<ODEdge> edges_to;
};

class OrientedDistanceClusterer::ODEdge {
public:
    ODEdge(size_t to_idx, int32_t weight) :
            to_idx(to_idx), weight(weight) {}
    ODEdge() = default;
    ~ODEdge() = default;
    
    /// Index of the node that the edge points to
    size_t to_idx;
    
    /// Weight for dynamic programming
    int32_t weight;
};

struct OrientedDistanceClusterer::DPScoreComparator {
private:
    const vector<ODNode>& nodes;
public:
    DPScoreComparator() = delete;
    DPScoreComparator(const vector<ODNode>& nodes) : nodes(nodes) {}
    ~DPScoreComparator() {}
    inline bool operator()(const size_t i, const size_t j) {
        return nodes[i].dp_score < nodes[j].dp_score;
    }
};

}

#endif
