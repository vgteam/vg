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

/**
 * Iterate over pairsets of integers in a pseudorandom but deterministic order.
 * We use the same permutation every time for a given number of items to pair
 * up.
 */
class ShuffledPairs {
public:
    /**
     * Make a new iterable pairing up the given number of items.
     */
    ShuffledPairs(size_t num_items);
    
    /**
     * Actual iterator class.
     */
    class iterator {
    public:
        /**
         * Advance to the next pair.
         */
        iterator& operator++();
        
        /**
         * Get the pair pointed to.
         */
        pair<size_t, size_t> operator*() const;
        
        /**
         * see if two iterators are equal.
         */
        bool operator==(const iterator& other) const;
        
        /**
         * see if two iterators are not equal.
         */
        bool operator!=(const iterator& other) const;
        
        friend class ShuffledPairs;
        
        // Default copy constructor and assignment operator.
        iterator(const iterator& other) = default;
        iterator& operator=(const iterator& other) = default;
        
    private:
        // Which permutation does this iterator mean?
        size_t permutation_idx = 0;
        // What are we iterating over?
        const ShuffledPairs& iteratee;
        
        // Make an iterator. Only the friend parent can do it.
        iterator(const ShuffledPairs& iteratee, size_t start_at);
    };
    
    using const_iterator = iterator;
    
    /**
     * Get an iterator to the first pair.
     */
    iterator begin() const;
    /**
     * Get an iterator to the past-the-end pair.
     */
    iterator end() const;
    
private:
    // How many items are we working with to pair up?
    size_t num_items;
    // How many pairs do we make?
    size_t num_pairs;
    // What's the power of 2 with 1 more bit than needed to represent each pair?
    size_t range_max;
};

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
bool clusters_overlap_in_read(const vector<MaximalExactMatch>& cluster1,
                              const vector<MaximalExactMatch>& cluster2);
bool clusters_overlap_in_graph(const vector<MaximalExactMatch>& cluster1,
                               const vector<MaximalExactMatch>& cluster2);
vector<pos_t> cluster_nodes(const vector<MaximalExactMatch>& cluster);

class MEMChainModelVertex {
public:
    MaximalExactMatch mem;
    vector<pair<MEMChainModelVertex*, double> > next_cost; // for forward
    vector<pair<MEMChainModelVertex*, double> > prev_cost; // for backward
    double weight;
    double score;
    int64_t approx_position;
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
    map<int64_t, vector<vector<MEMChainModelVertex>::iterator> > approx_positions;
    set<vector<MEMChainModelVertex>::iterator> redundant_vertexes;
    MEMChainModel(
        const vector<size_t>& aln_lengths,
        const vector<vector<MaximalExactMatch> >& matches,
        const function<int(pos_t)>& approx_position,
        const function<double(const MaximalExactMatch&, const MaximalExactMatch&)>& transition_weight,
        int band_width = 10,
        int position_depth = 1,
        int max_connections = 20);
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
    
    /// Each hit contains a pointer to the original MEM and the position of that
    /// particular hit in the graph.
    using hit_t = pair<const MaximalExactMatch*, pos_t>;
    
    /// Each cluster is a vector of hits.
    using cluster_t = vector<hit_t>;
                              
    /// Returns a vector of clusters. Each cluster is represented a vector of MEM hits. Each hit
    /// contains a pointer to the original MEM and the position of that particular hit in the graph.
    vector<cluster_t> clusters(int32_t max_qual_score = 60);
    
    /**
     * Returns a vvector of pairs of clusters from this clusterer and the other
     * clusterer that are within max_inter_cluster_distance of each other, and
     * are on opposing strands (as would be expected for read pairs).
     */
    vector<pair<cluster_t, cluster_t>> paired_clusters(OrientedDistanceClusterer& other, xg::XG* xgindex,
        size_t max_inter_cluster_distance, int32_t max_qual_score = 60);
    
private:
    class ODNode;
    class ODEdge;
    struct DPScoreComparator;
    
    /**
     * Given a certain number of items, and a callback to get each item's
     * position, build a distance forrest with trees for items that we can
     * verify are on the same strand of the same molecule.
     *
     * We use the distance approximation to cluster the MEM hits according to
     * the strand they fall on using the oriented distance estimation function
     * in xg.
     *
     * Returns a map from item pair (lower number first) to distance (which may
     * be negative) from the first to the second along the items' forward
     * strand.
     */
    unordered_map<pair<size_t, size_t>, int64_t> get_on_strand_distance_tree(size_t num_items, xg::XG* xgindex,
        const function<pos_t(size_t)>& get_position);
        
    /**
     * Given a number of nodes, and a map from node pair to signed relative
     * distance on a consistent strand (defining a forrest of trees, as
     * generated by get_on_strand_distance_tree()), flatten all the trees.
     *
     * Returns a vector of maps from node ID to relative position in linear
     * space, one map per input tree.
     *
     * Assumes all the distances are transitive, even though this isn't quite
     * true in graph space.
     */
    vector<unordered_map<size_t, int64_t>> flatten_distance_tree(size_t num_items,
        const unordered_map<pair<size_t, size_t>, int64_t>& recorded_finite_dists);
    
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
