#ifndef VG_ALGORITHMS_READ_COORDINATE_RMQST_HPP_INCLUDED
#define VG_ALGORITHMS_READ_COORDINATE_RMQST_HPP_INCLUDED

#include <cstddef>
#include <vector>

namespace vg {
namespace algorithms {

// Offline segment tree implementation of ReadCoordinateRMQ.
//
// All coordinates must be supplied at construction time for coordinate
// compression. Scores are then inserted online (one per anchor, left to right).
//
// Each segment tree node keeps the top-RMQ_K candidates by score_with_penalty
// for its coordinate sub-range. A query over [min_coord, max_coord] decomposes
// into O(log N) nodes and merges their heaps, costing O(RMQ_K · log N) total —
// versus O(N) per query for the AVL approach when all anchors fall in range
// (e.g. chrM with a huge max_lookback_bases).
class ReadCoordinateRMQ {
public:
    // all_coordinates: read_end values of every anchor in the current
    // chaining group (used only for coordinate compression).
    ReadCoordinateRMQ(double gap_extend,
                      const std::vector<size_t>& all_coordinates);
    ~ReadCoordinateRMQ() = default;

    void insert(size_t coordinate, int score, size_t item_index);

    std::vector<size_t> query_top_k_candidates(size_t max_coordinate,
                                               size_t max_lookback_distance,
                                               int k_candidates) const;

private:
    // Buffer per node. Must be >= k_candidates passed to query_top_k_candidates.
    static constexpr int RMQ_K = 64;

    struct Candidate {
        int  score_with_penalty;
        int  score;
        size_t item_index;
    };

    // Min-heap comparator: smallest score_with_penalty at top (= worst in top-K).
    struct CmpMin {
        bool operator()(const Candidate& a, const Candidate& b) const {
            return a.score_with_penalty > b.score_with_penalty;
        }
    };

    double              gap_extend;
    std::vector<size_t> sorted_coords;  // sorted, deduplicated coordinates
    int                 seg_n;          // power of 2 >= sorted_coords.size()

    // 1-indexed segment tree: root = 1, children of i are 2i and 2i+1,
    // leaves are [seg_n, 2*seg_n - 1].
    // Flat layout: node i occupies seg_flat[i*RMQ_K .. i*RMQ_K+RMQ_K-1],
    // seg_sizes[i] tracks how many candidates are stored (heap front = worst).
    std::vector<Candidate> seg_flat;   // 2 * seg_n * RMQ_K entries
    std::vector<int>       seg_sizes;  // 2 * seg_n counters

    // Insert c into a node's min-heap, evicting the worst if at capacity.
    void push_to_node(int node_idx, const Candidate& c);
};

} // namespace algorithms
} // namespace vg

#endif
