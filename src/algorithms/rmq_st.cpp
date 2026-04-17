#include "rmq_st.hpp"
#include <algorithm>
#include <queue>

namespace vg {
namespace algorithms {

ReadCoordinateRMQ::ReadCoordinateRMQ(double gap_extend,
                                     const std::vector<size_t>& all_coordinates)
    : gap_extend(gap_extend) {

    sorted_coords = all_coordinates;
    std::sort(sorted_coords.begin(), sorted_coords.end());
    sorted_coords.erase(
        std::unique(sorted_coords.begin(), sorted_coords.end()),
        sorted_coords.end());

    seg_n = 1;
    while (seg_n < (int)sorted_coords.size()) seg_n <<= 1;

    // 2 * seg_n nodes, 1-indexed (index 0 unused).
    // Flat layout: node i occupies seg_flat[i*RMQ_K .. i*RMQ_K+RMQ_K-1].
    seg_flat.resize(2 * seg_n * RMQ_K);
    seg_sizes.assign(2 * seg_n, 0);
}

void ReadCoordinateRMQ::push_to_node(int node_idx, const Candidate& c) {
    Candidate* heap = seg_flat.data() + node_idx * RMQ_K;
    int& sz = seg_sizes[node_idx];
    if (sz < RMQ_K) {
        heap[sz++] = c;
        std::push_heap(heap, heap + sz, CmpMin{});
    } else if (c.score_with_penalty > heap[0].score_with_penalty) {
        std::pop_heap(heap, heap + sz, CmpMin{});
        heap[sz - 1] = c;
        std::push_heap(heap, heap + sz, CmpMin{});
    }
}

void ReadCoordinateRMQ::insert(size_t coordinate, int score, size_t item_index) {
    const int swp = score + (int)(gap_extend * (double)coordinate);
    const Candidate c{swp, score, item_index};

    const int pos = (int)(
        std::lower_bound(sorted_coords.begin(), sorted_coords.end(), coordinate)
        - sorted_coords.begin());

    // Walk from leaf to root: each ancestor's heap covers this coordinate,
    // so it is a candidate for any query whose range includes it.
    for (int i = seg_n + pos; i >= 1; i >>= 1) {
        push_to_node(i, c);
    }
}

std::vector<size_t> ReadCoordinateRMQ::query_top_k_candidates(
    size_t max_coordinate, size_t max_lookback_distance, int k_candidates) const {

    if (sorted_coords.empty() || k_candidates <= 0) return {};

    const size_t min_coord = (max_coordinate > max_lookback_distance)
                             ? max_coordinate - max_lookback_distance : 0;

    // Compress the query range.
    const int ql = (int)(
        std::lower_bound(sorted_coords.begin(), sorted_coords.end(), min_coord)
        - sorted_coords.begin());
    const int qr = (int)(
        std::upper_bound(sorted_coords.begin(), sorted_coords.end(), max_coordinate)
        - sorted_coords.begin()) - 1;

    if (ql > qr || ql >= (int)sorted_coords.size()) return {};

    // Min-heap of size k_candidates: tracks the running top-k across nodes.
    std::priority_queue<Candidate, std::vector<Candidate>, CmpMin> result;

    // Standard iterative segment tree decomposition of [ql, qr].
    // The decomposed nodes are disjoint in coordinate space, so each inserted
    // item appears in exactly one of them — no deduplication needed.
    //
    // Key property for chrM: when [ql, qr] spans the entire coordinate range,
    // the loop resolves to the single root node, costing O(RMQ_K) instead of O(N).
    for (int l = seg_n + ql, r = seg_n + qr + 1; l < r; l >>= 1, r >>= 1) {
        auto process = [&](int node_idx) {
            const Candidate* heap = seg_flat.data() + node_idx * RMQ_K;
            const int sz = seg_sizes[node_idx];
            for (int i = 0; i < sz; i++) {
                const Candidate& c = heap[i];
                if ((int)result.size() < k_candidates) {
                    result.push(c);
                } else if (c.score_with_penalty > result.top().score_with_penalty) {
                    result.pop();
                    result.push(c);
                }
            }
        };
        if (l & 1) process(l++);
        if (r & 1) process(--r);
    }

    std::vector<size_t> out;
    out.reserve(result.size());
    while (!result.empty()) {
        out.push_back(result.top().item_index);
        result.pop();
    }
    std::reverse(out.begin(), out.end());
    return out;
}

} // namespace algorithms
} // namespace vg
