#ifndef VG_ALGORITHMS_READ_COORDINATE_RMQ_HPP_INCLUDED
#define VG_ALGORITHMS_READ_COORDINATE_RMQ_HPP_INCLUDED

#include <cstddef>
#include <memory>
#include <queue>
#include <vector>

namespace vg {
namespace algorithms {

class ReadCoordinateRMQ {
public:
  ReadCoordinateRMQ(int gap_extend);
  ~ReadCoordinateRMQ();

  void insert(size_t linear_coordinate, int partial_score, size_t item_index, size_t chain_coordinate);

  std::vector<size_t> query_top_k_candidates(size_t current_coordinate,
                                             size_t max_lookback_distance,
                                             int k_candidates) const;

private:
  struct Node {
    size_t coordinate;
    int score;
    int score_with_penalty;
    size_t item_index;

    int max_subtree_score;
    int height;

    std::unique_ptr<Node> left;
    std::unique_ptr<Node> right;

    Node(size_t coord, int s, size_t idx, int gap_extend, size_t chain_coord);
  };

  std::unique_ptr<Node> root;
  int gap_extend;

  int get_height(const std::unique_ptr<Node> &node) const;
  int get_max_subtree_score(const std::unique_ptr<Node> &node) const;
  void update_node_metadata(std::unique_ptr<Node> &node);

  std::unique_ptr<Node> rotate_right(std::unique_ptr<Node> y);
  std::unique_ptr<Node> rotate_left(std::unique_ptr<Node> x);
  int get_balance_factor(const std::unique_ptr<Node> &node) const;

  std::unique_ptr<Node> insert_recursive(std::unique_ptr<Node> node,
                                         size_t coordinate, int score,
                                         size_t item_index, size_t chain_coordinate);

  struct Candidate {
    int score;
    int score_with_penalty;
    size_t item_index;
    bool operator>(const Candidate &other) const {
      return score_with_penalty > other.score_with_penalty;
    }
  };

  void
  query_recursive(const std::unique_ptr<Node> &node, size_t min_coordinate,
                  size_t max_coordinate, int k_candidates,
                  std::priority_queue<Candidate, std::vector<Candidate>,
                                      std::greater<Candidate>> &min_heap) const;
};

} // namespace algorithms
} // namespace vg

#endif