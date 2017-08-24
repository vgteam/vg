#ifndef VG_MEM_HPP_INCLUDED
#define VG_MEM_HPP_INCLUDED

#include "gcsa/gcsa.h"
#include "gcsa/lcp.h"
#include "position.hpp"
#include "gssw_aligner.hpp"
#include "utility.hpp"

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

// Prime numbers spaced at approximately logarithmic intervals
static constexpr size_t spaced_primes[62] = {2ull, 5ull, 13ull, 29ull, 53ull, 127ull, 227ull, 487ull, 967ull, 2039ull, 4093ull, 8191ull, 16381ull, 32749ull, 65521ull, 131071ull, 262139ull, 524287ull, 1048573ull, 2097143ull, 4194301ull, 8388593ull, 16777213ull, 33554393ull, 67108859ull, 134217689ull, 268435399ull, 536870909ull, 1073741789ull, 2147483647ull, 4294967291ull, 8589934583ull, 17179869143ull, 34359738337ull, 68719476731ull, 137438953447ull, 274877906899ull, 549755813881ull, 1099511627689ull, 2199023255531ull, 4398046511093ull, 8796093022151ull, 17592186044399ull, 35184372088777ull, 70368744177643ull, 140737488355213ull, 281474976710597ull, 562949953421231ull, 1125899906842597ull, 2251799813685119ull, 4503599627370449ull, 9007199254740881ull, 18014398509481951ull, 36028797018963913ull, 72057594037927931ull, 144115188075855859ull, 288230376151711717ull, 576460752303423433ull, 1152921504606846883ull, 2305843009213693951ull, 4611686018427387847ull, 9223372036854775783ull};

// Precomputed primitive roots of unity paired with these primes (chosen randomly from the 20 smallest roots)
static constexpr size_t primitive_roots_of_unity[62] = {1ull, 3ull, 2ull, 21ull, 27ull, 56ull, 17ull, 45ull, 40ull, 28ull, 69ull, 70ull, 40ull, 31ull, 119ull, 75ull, 42ull, 61ull, 60ull, 46ull, 21ull, 13ull, 39ull, 13ull, 29ull, 15ull, 29ull, 32ull, 37ull, 73ull, 56ull, 45ull, 13ull, 90ull, 51ull, 12ull, 32ull, 11ull, 39ull, 24ull, 8ull, 39ull, 7ull, 51ull, 38ull, 67ull, 2ull, 34ull, 62ull, 19ull, 13ull, 30ull, 12ull, 45ull, 31ull, 57ull, 6ull, 57ull, 3ull, 37ull, 68ull, 54ull};

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
        // What is the ordinal value of this element in the permutation?
        size_t permutation_idx;
        // Which permutation does this iterator mean?
        size_t permuted;
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
    // A prime number that is larger than the number of pairs
    size_t larger_prime;
    // A primitive root of unity for the prime
    size_t primitive_root;
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
    
    OrientedDistanceClusterer(const Alignment& alignment,
                              const vector<MaximalExactMatch>& mems,
                              const Aligner& aligner,
                              xg::XG* xgindex,
                              size_t max_expected_dist_approx_error = 8);
    
    /// Each hit contains a pointer to the original MEM and the position of that
    /// particular hit in the graph.
    using hit_t = pair<const MaximalExactMatch*, pos_t>;
    
    /// Each cluster is a vector of hits.
    using cluster_t = vector<hit_t>;
                              
    /// Returns a vector of clusters. Each cluster is represented a vector of MEM hits. Each hit
    /// contains a pointer to the original MEM and the position of that particular hit in the graph.
    vector<cluster_t> clusters(int32_t max_qual_score = 60, int32_t log_likelihood_approx_factor = 0);
    
    /**
     * Given two vectors of clusters, an xg index, and a maximum separation,
     * returns a vector of pairs of cluster numbers (one in each vector) that
     * are in opposite orientations (as would be expected of read pairs) and
     * within the specified distance.
     */
    static vector<pair<size_t, size_t>> pair_clusters(const vector<cluster_t*>& our_clusters,
        const vector<cluster_t*>& their_clusters, xg::XG* xgindex,
        size_t max_inter_cluster_distance);
    
    //static size_t PRUNE_COUNTER;
    //static size_t CLUSTER_TOTAL;
    
private:
    class ODNode;
    class ODEdge;
    struct DPScoreComparator;
    
    /// Internal constructor that public constructors filter into
    OrientedDistanceClusterer(const Alignment& alignment,
                              const vector<MaximalExactMatch>& mems,
                              const Aligner* aligner,
                              const QualAdjAligner* qual_adj_aligner,
                              xg::XG* xgindex,
                              size_t max_expected_dist_approx_error);
    
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
    static unordered_map<pair<size_t, size_t>, int64_t> get_on_strand_distance_tree(size_t num_items, xg::XG* xgindex,
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
    static vector<unordered_map<size_t, int64_t>> flatten_distance_tree(size_t num_items,
        const unordered_map<pair<size_t, size_t>, int64_t>& recorded_finite_dists);
    
    /// Fills the tail score fields in nodes
    vector<pair<size_t, size_t>> compute_tail_mem_coverage(const Alignment& alignment,
                                                           const vector<MaximalExactMatch>& mems);
    
    /// Fills input vectors with indices of source and sink nodes
    void identify_sources_and_sinks(vector<size_t>& sources_out, vector<size_t>& sinks_out);
    
    /// Identify weakly connected components in the graph
    void connected_components(vector<vector<size_t>>& components_out);
    
    /// Fills the input vector with the indices of a topological sort
    void topological_order(vector<size_t>& order_out);
    
    /// Perform dynamic programming and store scores in nodes
    void perform_dp();
    
    vector<ODNode> nodes;
    
    const Aligner* aligner;
    const QualAdjAligner* qual_adj_aligner;
};

class OrientedDistanceClusterer::ODNode {
public:
    ODNode(const MaximalExactMatch& mem, pos_t start_pos, int32_t score) :
            mem(&mem), start_pos(start_pos), score(score) {}
    ODNode() = default;
    ~ODNode() = default;
    
    const MaximalExactMatch* mem;
    
    /// Position of GCSA hit in the graph
    pos_t start_pos;
    
    /// Score of the exact match this node represents
    int32_t score;
    
    /// Lower bound on score used in dynamic programming
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
