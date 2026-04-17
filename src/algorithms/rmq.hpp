#ifndef VG_ALGORITHMS_READ_COORDINATE_RMQ_HPP_INCLUDED
#define VG_ALGORITHMS_READ_COORDINATE_RMQ_HPP_INCLUDED

#include <cstddef>
#include <cstdint>
#include <queue>
#include <vector>

namespace vg {
namespace algorithms {

class ReadCoordinateRMQAVL {
public:
  ReadCoordinateRMQAVL(int gap_extend);
  ~ReadCoordinateRMQAVL();

  void reserve(size_t n);
  void insert(size_t linear_coordinate, int partial_score, size_t item_index);

  // max_nodes_to_visit: hard cap on AVL nodes visited per query (0 =
  // unlimited). Analogous to lookback_item_hard_cap in chain_items_dp.
  std::vector<size_t>
  query_top_k_candidates(size_t current_coordinate,
                         size_t max_lookback_distance, int k_candidates,
                         size_t max_nodes_to_visit = 0) const;

private:
  struct Node {
    size_t coordinate;
    int score;
    int score_with_penalty;
    size_t item_index;

    int max_subtree_score;
    int height;

    int32_t left;
    int32_t right;

    Node(size_t coord, int s, size_t idx, int gap_extend);
  };

  std::vector<Node> nodes;
  int32_t root = -1;
  int gap_extend;

  int get_height(int32_t idx) const;
  int get_max_subtree_score(int32_t idx) const;
  void update_node_metadata(int32_t idx);

  int32_t rotate_right(int32_t y);
  int32_t rotate_left(int32_t x);
  int get_balance_factor(int32_t idx) const;

  int32_t insert_recursive(int32_t node_idx, size_t coordinate, int score,
                           size_t item_index);

  struct Candidate {
    int score;
    int score_with_penalty;
    size_t item_index;
    bool operator>(const Candidate &other) const {
      return score_with_penalty > other.score_with_penalty;
    }
  };

  void
  query_iterative(size_t min_coordinate, size_t max_coordinate,
                  int k_candidates, size_t max_nodes_to_visit,
                  std::priority_queue<Candidate, std::vector<Candidate>,
                                      std::greater<Candidate>> &min_heap) const;
};

} // namespace algorithms
} // namespace vg

#endif
