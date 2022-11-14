#ifndef VG_CLUSTER_HPP_INCLUDED
#define VG_CLUSTER_HPP_INCLUDED

#include <gcsa/gcsa.h>
#include <gcsa/lcp.h>
#include <structures/union_find.hpp>
#include "position.hpp"
#include "aligner.hpp"
#include "mem.hpp"
#include "handle.hpp"
#include "snarl_distance_index.hpp"
#include "snarl_seed_clusterer.hpp"
#include "path_component_index.hpp"
#include "bdsg/hash_graph.hpp"

#include <functional>
#include <string>
#include <vector>
#include <map>


/**
 * \file cluster.hpp
 *
 * Chaining and clustering tools to work with maximal exact matches.
 */

namespace vg {
    
using namespace std;
using namespace structures;

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
        
        // TODO: This gets implicitly deleted and generates warning because of the const reference
        // member variable
        //iterator& operator=(const iterator& other) = default;
        
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

class MEMChainModelVertex {
public:
    MaximalExactMatch mem;
    vector<pair<MEMChainModelVertex*, double> > next_cost; // for forward
    vector<pair<MEMChainModelVertex*, double> > prev_cost; // for backward
    double weight;
    double score;
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
    unordered_map<path_handle_t, map<int64_t, vector<vector<MEMChainModelVertex>::iterator> > > positions;
    set<vector<MEMChainModelVertex>::iterator> redundant_vertexes;
    MEMChainModel(
        const vector<size_t>& aln_lengths,
        const vector<vector<MaximalExactMatch> >& matches,
        const function<int64_t(pos_t)>& approx_position,
        const function<unordered_map<path_handle_t, vector<pair<size_t, bool> > >(pos_t)>& path_position,
        const function<double(const MaximalExactMatch&, const MaximalExactMatch&)>& transition_weight,
        int band_width = 10,
        int position_depth = 1,
        int max_connections = 20);
    void score(const unordered_set<MEMChainModelVertex*>& exclude);
    MEMChainModelVertex* max_vertex(void);
    vector<vector<MaximalExactMatch> > traceback(int alt_alns, bool paired, bool debug);
    void display(ostream& out);
    void display_dot(ostream& out, vector<MEMChainModelVertex*> vertex_trace);
    void clear_scores(void);
};

/*
 * A base class to hold some shared methods and data types between the TVS,
 * oriented distance, and minimum distance clusterers.
 */
class MEMClusterer {
public:
    MEMClusterer() = default;
    virtual ~MEMClusterer() = default;
    
    /// Each hit contains a pointer to the original MEM and the position of that
    /// particular hit in the graph.
    using hit_t = pair<const MaximalExactMatch*, pos_t>;
    
    /// Each cluster is a vector of hits and a paired multiplicity
    using cluster_t = pair<vector<hit_t>, double>;
    
    /// Represents the mismatches that were allowed in "MEMs" from the fanout
    /// match algorithm
    using match_fanouts_t = unordered_map<const MaximalExactMatch*, deque<pair<string::const_iterator, char>>>;
    
    /// Returns a vector of clusters. Each cluster is represented a vector of MEM hits. Each hit
    /// contains a pointer to the original MEM and the position of that particular hit in the graph.
    vector<cluster_t> clusters(const Alignment& alignment,
                               const vector<MaximalExactMatch>& mems,
                               const GSSWAligner* Aligner,
                               size_t min_mem_length = 1,
                               int32_t max_qual_score = 60,
                               int32_t log_likelihood_approx_factor = 0,
                               size_t min_median_mem_coverage_for_split = 0,
                               double suboptimal_edge_pruning_factor = .75,
                               double cluster_multiplicity_diff = 10.0,
                               const match_fanouts_t* fanouts = nullptr);
    
    /**
     * Given two vectors of clusters and bounds on the distance between clusters,
     * returns a vector of pairs of cluster numbers (one in each vector) matched with the estimated
     * distance.
     *
     * Clusters are assumed to be located at the position of the first MEM hit they contain. Optionally,
     * additional MEMs may be identied as possible anchors for the cluster. Additional anchors are
     * provided as pairs of (cluster index, MEM index within cluster). Only one result will be returned
     * per pair of clusters regardless of how many alternate anchors are given.
     */
    virtual vector<pair<pair<size_t, size_t>, int64_t>> pair_clusters(const Alignment& alignment_1,
                                                                      const Alignment& alignment_2,
                                                                      const vector<cluster_t*>& left_clusters,
                                                                      const vector<cluster_t*>& right_clusters,
                                                                      const vector<pair<size_t, size_t>>& left_alt_cluster_anchors,
                                                                      const vector<pair<size_t, size_t>>& right_alt_cluster_anchors,
                                                                      int64_t optimal_separation,
                                                                      int64_t max_deviation) = 0;
    
    /// The largest discrepency we will allow between the read-implied distances and the estimated  gap distance
    int64_t max_gap = numeric_limits<int64_t>::max();
    
protected:
    
    class HitNode;
    class HitEdge;
    class HitGraph;
    class DPScoreComparator;
    
    /// Initializes a hit graph and adds edges to it, this must be implemented by any inheriting
    /// class
    virtual HitGraph make_hit_graph(const Alignment& alignment, const vector<MaximalExactMatch>& mems,
                                    const GSSWAligner* aligner, size_t min_mem_length,
                                    const match_fanouts_t* fanouts) = 0;
    
    /// Once the distance between two hits has been estimated, estimate the score of the hit graph edge
    /// connecting them
    int32_t estimate_edge_score(const MaximalExactMatch* mem_1, const MaximalExactMatch* mem_2, int64_t graph_dist,
                                const GSSWAligner* aligner) const;
    
    /// Sorts cluster pairs and removes copies of the same cluster pair, choosing only the one whose distance
    /// is closest to the optimal separation
    void deduplicate_cluster_pairs(vector<pair<pair<size_t, size_t>, int64_t>>& cluster_pairs, int64_t optimal_separation);
};
    
class MEMClusterer::HitGraph {
public:
    
    /// Initializes nodes in the hit graph, but does not add edges
    HitGraph(const vector<MaximalExactMatch>& mems, const Alignment& alignment, const GSSWAligner* aligner,
             size_t min_mem_length = 1, bool track_components = false, const match_fanouts_t* fanouts = nullptr);
    
    /// Add an edge
    void add_edge(size_t from, size_t to, int32_t weight, int64_t distance);
    
    /// Returns the top scoring connected components
    vector<cluster_t> clusters(const Alignment& alignment,
                               const GSSWAligner* aligner,
                               int32_t max_qual_score,
                               int32_t log_likelihood_approx_factor,
                               size_t min_median_mem_coverage_for_split,
                               double suboptimal_edge_pruning_factor,
                               double cluster_multiplicity_diff);
    
    vector<HitNode> nodes;
    
private:
    
    /// Identify weakly connected components in the graph
    void connected_components(vector<vector<size_t>>& components_out) const;
    
    /// Prune edges that are not on any traceback that scores highly compared to the best score in the component,
    /// splits up the components (adding some to the end of the vector) if doing so splits a component
    void prune_low_scoring_edges(vector<vector<size_t>>& components, size_t component_idx, double score_factor);
    
    /// Perform dynamic programming and store scores in nodes
    void perform_dp();
    
    /// Fills input vectors with indices of source and sink nodes
    void identify_sources_and_sinks(vector<size_t>& sources_out, vector<size_t>& sinks_out) const;
    
    /// Fills the input vector with the indices of a topological sort
    void topological_order(vector<size_t>& order_out) const;
    
    /// Computes the topological order of
    void component_topological_order(const vector<size_t>& component, vector<size_t>& order_out) const;
    
    /// Returns the median coverage of bases in the reads by bases in the cluster, attempts to remove apparent
    /// redundant sub-MEMs
    size_t median_mem_coverage(const vector<size_t>& component, const Alignment& aln) const;
    
    /// Should we actively keep track of connected components?
    bool track_components;
    
    /// Keeps track of the connected components
    UnionFind components;
};
    
class MEMClusterer::HitNode {
public:
    HitNode(const MaximalExactMatch& mem, pos_t start_pos, int32_t score) : mem(&mem), start_pos(start_pos), score(score) { }
    HitNode() = default;
    ~HitNode() = default;
    
    const MaximalExactMatch* mem;
    
    /// Position of GCSA hit in the graph
    pos_t start_pos;
    
    /// Score of the exact match this node represents
    int32_t score;
    
    /// Score used in dynamic programming
    int32_t dp_score;
    
    /// Edges from this node that are colinear with the read
    vector<HitEdge> edges_from;
    
    /// Edges to this node that are colinear with the read
    vector<HitEdge> edges_to;
};

class MEMClusterer::HitEdge {
public:
    HitEdge(size_t to_idx, int32_t weight, int64_t distance) : to_idx(to_idx), weight(weight), distance(distance) {}
    HitEdge() = default;
    ~HitEdge() = default;
    
    /// Index of the node that the edge points to
    size_t to_idx;
    
    /// Weight for dynamic programming
    int32_t weight;
    
    /// Estimated distance
    int64_t distance;
};

struct MEMClusterer::DPScoreComparator {
private:
    const vector<HitNode>& nodes;
public:
    DPScoreComparator() = delete;
    DPScoreComparator(const vector<HitNode>& nodes) : nodes(nodes) {}
    ~DPScoreComparator() {}
    inline bool operator()(const size_t i, const size_t j) {
        return nodes[i].dp_score < nodes[j].dp_score;
    }
};
        
/*
 * A clustering implementation that actually doesn't do any clustering
 */
class NullClusterer : public MEMClusterer {
public:
    NullClusterer() = default;
    virtual ~NullClusterer() = default;

    /// Concrete implementation of virtual method from MEMClusterer
    vector<pair<pair<size_t, size_t>, int64_t>> pair_clusters(const Alignment& alignment_1,
                                                              const Alignment& alignment_2,
                                                              const vector<cluster_t*>& left_clusters,
                                                              const vector<cluster_t*>& right_clusters,
                                                              const vector<pair<size_t, size_t>>& left_alt_cluster_anchors,
                                                              const vector<pair<size_t, size_t>>& right_alt_cluster_anchors,
                                                              int64_t optimal_separation,
                                                              int64_t max_deviation);

protected:

    /// Concrete implementation of virtual method from MEMClusterer
    /// Note: ignores the min_mem_length parameter
    HitGraph make_hit_graph(const Alignment& alignment, const vector<MaximalExactMatch>& mems, const GSSWAligner* aligner,
                            size_t min_mem_length, const match_fanouts_t* fanouts);
};
    
    
/*
 * An abstract class that provides distances to the oriented distance clusterer
 */
class OrientedDistanceMeasurer {
public:
    virtual ~OrientedDistanceMeasurer() = default;
    
    /// Returns a signed distance, where positive indicates that pos_2 is to the right
    /// of pos_1, and negative indicates to the left. If the distance is infinite or
    /// can't be determined, returns numeric_limits<int64_t>::max().
    virtual int64_t oriented_distance(const pos_t& pos_1, const pos_t& pos_2) = 0;
    
    /// Return a vector of groups that we believe will have finite distances under this metric,
    /// can be empty.
    virtual vector<vector<size_t>> get_buckets(const function<pos_t(size_t)>& get_position, size_t num_items) = 0;
    
    /// Return a vector of pairs of groups (referred to by indexes in the current_groups vector)
    /// that cannot have finite distances between them (typically because they are on separate components).
    virtual vector<pair<size_t, size_t>> exclude_merges(vector<vector<size_t>>& current_groups,
                                                        const function<pos_t(size_t)>& get_position) = 0;
};

/*
 * A distance function that uses an a graph's embedded paths to measure distances, either in a stranded
 * or unstranded manner.
 */
class PathOrientedDistanceMeasurer : public OrientedDistanceMeasurer {

public:
    
    /// Construct a distance service to measures distance along paths in this graph. Optionally
    /// measures all distances on the forward strand of the paths.
    PathOrientedDistanceMeasurer(const PathPositionHandleGraph* graph,
                                 const PathComponentIndex* path_component_index = nullptr);
    
    /// Default desctructor
    ~PathOrientedDistanceMeasurer() = default;
    
    /// Returns a signed distance, where positive indicates that pos_2 is to the right
    /// of pos_1, and negative indicates to the left. If the distance is infinite or
    /// can't be determined, returns numeric_limits<int64_t>::max().
    int64_t oriented_distance(const pos_t& pos_1, const pos_t& pos_2);
    
    /// Return a vector of groups that we believe will have finite distances under this metric,
    /// can be empty.
    vector<vector<size_t>> get_buckets(const function<pos_t(size_t)>& get_position, size_t num_items);
    
    /// Return a vector of pairs of groups (referred to by indexes in the current_groups vector)
    /// that cannot have finite distances between them (typically because they are on separate components).
    vector<pair<size_t, size_t>> exclude_merges(vector<vector<size_t>>& current_groups,
                                                const function<pos_t(size_t)>& get_position);
    
    /// The maximum distance we will walk trying to find a shared path
    size_t max_walk = 50;
    
private:
    
    const PathPositionHandleGraph* graph = nullptr;
    const PathComponentIndex* path_component_index = nullptr;

};
    
/*
 * A distance function that the minimum distance function provided by the Snarl-based
 * distance index
 */
class SnarlOrientedDistanceMeasurer : public OrientedDistanceMeasurer {

public:
    // Construct a distance service to measures distance as the minimum distance in the graph
    SnarlOrientedDistanceMeasurer(SnarlDistanceIndex* distance_index);
    
    /// Default desctructor
    ~SnarlOrientedDistanceMeasurer() = default;
    
    /// Returns a signed distance, where positive indicates that pos_2 is to the right
    /// of pos_1, and negative indicates to the left. If the distance is infinite or
    /// can't be determined, returns numeric_limits<int64_t>::max().
    int64_t oriented_distance(const pos_t& pos_1, const pos_t& pos_2);
    
    /// Return a vector of groups that we believe will have finite distances under this metric,
    /// can be empty.
    vector<vector<size_t>> get_buckets(const function<pos_t(size_t)>& get_position, size_t num_items);
    
    /// Return a vector of pairs of groups (referred to by indexes in the current_groups vector)
    /// that cannot have finite distances between them (typically because they are on separate components).
    vector<pair<size_t, size_t>> exclude_merges(vector<vector<size_t>>& current_groups,
                                                const function<pos_t(size_t)>& get_position);
    
private:
    
    SnarlDistanceIndex* distance_index = nullptr;
};
    
class OrientedDistanceClusterer : public MEMClusterer {
public:
    
    /// Constructor
    OrientedDistanceClusterer(OrientedDistanceMeasurer& distance_measurer,
                              size_t max_expected_dist_approx_error = 8);
    
    /// Concrete implementation of virtual method from MEMClusterer
    vector<pair<pair<size_t, size_t>, int64_t>> pair_clusters(const Alignment& alignment_1,
                                                              const Alignment& alignment_2,
                                                              const vector<cluster_t*>& left_clusters,
                                                              const vector<cluster_t*>& right_clusters,
                                                              const vector<pair<size_t, size_t>>& left_alt_cluster_anchors,
                                                              const vector<pair<size_t, size_t>>& right_alt_cluster_anchors,
                                                              int64_t optimal_separation,
                                                              int64_t max_deviation);
    
    //static size_t PRUNE_COUNTER;
    //static size_t CLUSTER_TOTAL;
    //static size_t MEM_FILTER_COUNTER;
    //static size_t MEM_TOTAL;
    //static size_t PRE_SPLIT_CLUSTER_COUNTER;
    //static size_t SPLIT_ATTEMPT_COUNTER;
    //static size_t SUCCESSFUL_SPLIT_ATTEMPT_COUNTER;
    //static size_t POST_SPLIT_CLUSTER_COUNTER;
    
protected:

    /**
     * Given a certain number of items, and a callback to get each item's
     * position, and a callback to a fixed offset from that position
     * build a distance forest with trees for items that we can
     * verify are on the same strand of the same molecule.
     *
     * We use the distance approximation to cluster the MEM hits according to
     * the strand they fall on using the oriented distance estimation function.
     *
     * Returns a map from item pair (lower number first) to distance (which may
     * be negative) from the first to the second along the items' forward
     * strand.
     */
    unordered_map<pair<size_t, size_t>, int64_t> get_on_strand_distance_tree(size_t num_items,
                                                                             const function<pos_t(size_t)>& get_position,
                                                                             const function<int64_t(size_t)>& get_offset);
    
    /**
     * Adds edges into the distance tree by estimating the distance between pairs
     * generated by a high entropy deterministic permutation
     */
    void extend_dist_tree_by_permutations(const function<pos_t(size_t)>& get_position,
                                          const function<int64_t(size_t)>& get_offset,
                                          size_t num_items,
                                          int64_t max_failed_distance_probes,
                                          size_t decrement_frequency,
                                          unordered_map<pair<size_t, size_t>, int64_t>& recorded_finite_dists,
                                          map<pair<size_t, size_t>, size_t>& num_infinite_dists,
                                          UnionFind& component_union_find,
                                          size_t& num_possible_merges_remaining);
    
    /**
     * Adds edges into the distance tree by estimating the distance only between pairs
     * of items that can be easily identified as having a finite distance (e.g. by sharing
     * a path)
     */
    void extend_dist_tree_by_buckets(const function<pos_t(size_t)>& get_position,
                                     const function<int64_t(size_t)>& get_offset,
                                     size_t num_items,
                                     unordered_map<pair<size_t, size_t>, int64_t>& recorded_finite_dists,
                                     UnionFind& component_union_find,
                                     size_t& num_possible_merges_remaining);
    
    /**
     * Automatically blocks off merges in the distance tree between groups that can be inferred
     * to be on separate components
     */
    void exclude_dist_tree_merges(const function<pos_t(size_t)>& get_position,
                                  map<pair<size_t, size_t>, size_t>& num_infinite_dists,
                                  UnionFind& component_union_find,
                                  size_t& num_possible_merges_remaining,
                                  int64_t max_failed_distance_probes);
    
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
    
    /// Returns a vector containing the number of SMEM beginnings to the left and the number of SMEM
    /// endings to the right of each read position
    vector<pair<size_t, size_t>> compute_tail_mem_coverage(const Alignment& alignment,
                                                           const vector<MaximalExactMatch>& mems);
    
    /// Concrete implementation of virtual method from MEMClusterer
    HitGraph make_hit_graph(const Alignment& alignment, const vector<MaximalExactMatch>& mems, const GSSWAligner* aligner,
                            size_t min_mem_length, const match_fanouts_t* fanouts);
    
    OrientedDistanceMeasurer& distance_measurer;
    size_t max_expected_dist_approx_error;
    bool unstranded;
    
};
    
/*
 * An abtract class that provides a heuristic distance between two positions. The semantics are
 * unspecified for what manner of "approximate" the heuristic is (upperbound, lowerbound, etc.).
 */
class DistanceHeuristic {
public:
    virtual ~DistanceHeuristic() = default;
    
    virtual int64_t operator()(const pos_t& pos_1, const pos_t& pos_2) = 0;
};
    
/*
 * An exact computation of the minimum distance between two positions using the snarl
 * decomposition
 */
class SnarlMinDistance : public DistanceHeuristic {
public:
    SnarlMinDistance() = delete;
    SnarlMinDistance(SnarlDistanceIndex& distance_index);
    ~SnarlMinDistance() = default;
    
    int64_t operator()(const pos_t& pos_1, const pos_t& pos_2);
private:
    SnarlDistanceIndex& distance_index;
};

/*
 * An upperbound on the distance between two positions computed using the distance
 * between those positions and tips. Strict upperbound in DAGs, only an upperbound
 * among a subset of paths in cyclic graphs (as distance is unbounded above).
 */
class TipAnchoredMaxDistance : public DistanceHeuristic {
public:
    TipAnchoredMaxDistance() = delete;
    TipAnchoredMaxDistance(SnarlDistanceIndex& distance_index);
    ~TipAnchoredMaxDistance() = default;
    
    int64_t operator()(const pos_t& pos_1, const pos_t& pos_2);
private:
    SnarlDistanceIndex& distance_index;
};

/*
 * Implements the heuristic solution to the Target Value Search problem described
 * in Kuhn, et al. (2008).
 */
class TargetValueSearch {
public:
    TargetValueSearch() = delete;
    TargetValueSearch(const HandleGraph& handle_graph,
                      DistanceHeuristic* upper_bound_heuristic,
                      DistanceHeuristic* lower_bound_heuristic);
    ~TargetValueSearch() = default;
    
    /// Does a path exist from pos_1 to pos_2 with length within the tolerance from the target value?
    bool tv_path_exists(const pos_t& pos_1, const pos_t& pos_2, int64_t target_value, int64_t tolerance);
    
    /// Returns the length of path from pos_1 to pos_2 with length closest to the target value. If there
    /// is no such path within the tolerance of the target value, returns numeric_limits<int64_t>::max().
    int64_t tv_path_length(const pos_t& pos_1, const pos_t& pos_2, int64_t target_value, int64_t tolerance);
    
    /// Returns a path from pos_1 to pos_2 with length closest to the target value. If there is no such
    /// path within the tolerance of the target value, returns an empty vector.
    vector<handle_t> tv_path(const pos_t& pos_1, const pos_t& pos_2, int64_t target_value, int64_t tolerance);
    
protected:
    
    vector<handle_t> tv_phase2(const pos_t& pos_1,
                               const pos_t& pos_2,
                               int64_t target_value,
                               int64_t tolerance,
                               hash_map<pair<id_t, bool>,int64_t>& node_to_target_shorter,
                               hash_map<pair<id_t, bool>, int64_t>& node_to_target_longer,
                               pair<int64_t, pair<pair<id_t, bool>,int64_t>>& best_lng,
                               pair<int64_t, pair<pair<id_t, bool>, int64_t>>& next_best,
                               hash_map<pair<pair<id_t, bool>, int64_t>, pair<pair<id_t, bool>, int64_t>>& node_to_path);
    
    const HandleGraph& handle_graph;
    unique_ptr<DistanceHeuristic> upper_bound_heuristic;
    unique_ptr<DistanceHeuristic> lower_bound_heuristic;
};
    
/*
 * A MEM clusterer built around the Target Value Search problem.
 */
class TVSClusterer : public MEMClusterer {
public:
    TVSClusterer(const HandleGraph* handle_graph, SnarlDistanceIndex* distance_index);
    ~TVSClusterer() = default;
    
    /// Concrete implementation of virtual method from MEMClusterer
    vector<pair<pair<size_t, size_t>, int64_t>> pair_clusters(const Alignment& alignment_1,
                                                              const Alignment& alignment_2,
                                                              const vector<cluster_t*>& left_clusters,
                                                              const vector<cluster_t*>& right_clusters,
                                                              const vector<pair<size_t, size_t>>& left_alt_cluster_anchors,
                                                              const vector<pair<size_t, size_t>>& right_alt_cluster_anchors,
                                                              int64_t optimal_separation,
                                                              int64_t max_deviation);
    
protected:
    
    /// Concrete implementation of virtual method from MEMClusterer
    HitGraph make_hit_graph(const Alignment& alignment, const vector<MaximalExactMatch>& mems, const GSSWAligner* aligner,
                            size_t min_mem_length, const match_fanouts_t* fanouts);
    
    TargetValueSearch tvs;
};

/*
 * A MEM clusterer based on finding the minimum distance between all pairs of seeds or clusters
 */
class MinDistanceClusterer : public MEMClusterer {
public:
    MinDistanceClusterer(SnarlDistanceIndex* distance_index);
    virtual ~MinDistanceClusterer() = default;
    
    /// Concrete implementation of virtual method from MEMClusterer
    vector<pair<pair<size_t, size_t>, int64_t>> pair_clusters(const Alignment& alignment_1,
                                                              const Alignment& alignment_2,
                                                              const vector<cluster_t*>& left_clusters,
                                                              const vector<cluster_t*>& right_clusters,
                                                              const vector<pair<size_t, size_t>>& left_alt_cluster_anchors,
                                                              const vector<pair<size_t, size_t>>& right_alt_cluster_anchors,
                                                              int64_t optimal_separation,
                                                              int64_t max_deviation);
    
protected:
    
    /// Concrete implementation of virtual method from MEMClusterer
    virtual HitGraph make_hit_graph(const Alignment& alignment, const vector<MaximalExactMatch>& mems, const GSSWAligner* aligner,
                                    size_t min_mem_length, const match_fanouts_t* fanouts);
    
    const HandleGraph* handle_graph;
    SnarlDistanceIndex* distance_index;
};

/*
 * A version of the MinDistanceClusterer that greedily agglomerates seeds into connected components
 * based on minimum distance, iterating over pairs in a sensible order
 */
class GreedyMinDistanceClusterer : public MinDistanceClusterer {
public:
    GreedyMinDistanceClusterer(SnarlDistanceIndex* distance_index);
    ~GreedyMinDistanceClusterer() = default;
    
protected:
    
    /// Concrete implementation of virtual method from MEMClusterer, overides the inherited one from MinDistanceClusterer
    HitGraph make_hit_graph(const Alignment& alignment, const vector<MaximalExactMatch>& mems, const GSSWAligner* aligner,
                            size_t min_mem_length, const match_fanouts_t* fanouts);
    
    
    /// How far apart do we expect the seeds to be on the read?
    const int64_t expected_separation = 20;
    
    /// How more bases would we search forward to find the next seed before we think
    /// it's worth searching 1 base backward?
    const int64_t forward_multiplier = 3;
    
    /// Minimum distance between two seeds on the read
    const int64_t min_separation = -10;
    
    /// Maximum distance between two seeds on the read
    const int64_t max_separation = 250;
    
};

/*
 * A version of the MinDistanceClusterer that uses the SeedClusterer to partition reads
 * into nearby clusters and only measures distances within clusters
 */
class ComponentMinDistanceClusterer : public MinDistanceClusterer {
public:
    ComponentMinDistanceClusterer(SnarlDistanceIndex* distance_index);
    ~ComponentMinDistanceClusterer() = default;
    
protected:
    
    /// Concrete implementation of virtual method from MEMClusterer, overides the inherited one from MinDistanceClusterer
    HitGraph make_hit_graph(const Alignment& alignment, const vector<MaximalExactMatch>& mems, const GSSWAligner* aligner,
                            size_t min_mem_length, const match_fanouts_t* fanouts);
    
    
    /// Minimum distance between two seeds on the read
    const int64_t min_read_separation = 0;
    
    /// The number of connections from one hit in a component to another that we will consider (0 for no maximum)
    const int64_t early_stop_number = 2;
};

/// get the handles that a mem covers
vector<pair<gcsa::node_type, size_t> > mem_node_start_positions(const HandleGraph& graph, const vg::MaximalExactMatch& mem);
/// return a containing subgraph connecting the mems
bdsg::HashGraph cluster_subgraph_containing(const HandleGraph& base, const Alignment& aln, const vector<vg::MaximalExactMatch>& cluster, const GSSWAligner* aligner);
/// return a subgraph for a cluster of MEMs from the given alignment
/// use walking to get the hits
bdsg::HashGraph cluster_subgraph_walk(const HandleGraph& base, const Alignment& aln, const vector<vg::MaximalExactMatch>& mems, double expansion);

}

#endif
