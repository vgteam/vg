#include "rmq.hpp"
#include <algorithm>
#include <limits>

namespace vg {
namespace algorithms {

using namespace std;

ReadCoordinateRMQAVL::Node::Node(size_t coord, int s, size_t idx, int gap_extend)
    : coordinate(coord), score(s),
      score_with_penalty(s + (int)(gap_extend * coord)), item_index(idx),
      max_subtree_score(score_with_penalty), height(1), left(-1), right(-1) {}

ReadCoordinateRMQAVL::ReadCoordinateRMQAVL(int gap_extend)
    : root(-1), gap_extend(gap_extend) {}

ReadCoordinateRMQAVL::~ReadCoordinateRMQAVL() {}

void ReadCoordinateRMQAVL::reserve(size_t n) { nodes.reserve(n); }

void ReadCoordinateRMQAVL::insert(size_t linear_coordinate, int partial_score,
                               size_t item_index) {
  root = insert_recursive(root, linear_coordinate, partial_score, item_index);
}

std::vector<size_t> ReadCoordinateRMQAVL::query_top_k_candidates(
    size_t current_coordinate, size_t max_lookback_distance, int k_candidates,
    size_t max_nodes_to_visit) const {
  if (k_candidates <= 0 || root == -1) {
    return {};
  }

  size_t min_coordinate = (current_coordinate > max_lookback_distance)
                              ? (current_coordinate - max_lookback_distance)
                              : 0;
  size_t max_coordinate = current_coordinate;

  std::priority_queue<Candidate, std::vector<Candidate>,
                      std::greater<Candidate>>
      min_heap;
  query_iterative(min_coordinate, max_coordinate, k_candidates,
                  max_nodes_to_visit, min_heap);

  std::vector<size_t> results;
  results.reserve(min_heap.size());
  while (!min_heap.empty()) {
    results.push_back(min_heap.top().item_index);
    min_heap.pop();
  }
  std::reverse(results.begin(), results.end());
  return results;
}

int ReadCoordinateRMQAVL::get_height(int32_t idx) const {
  return idx >= 0 ? nodes[idx].height : 0;
}

int ReadCoordinateRMQAVL::get_max_subtree_score(int32_t idx) const {
  return idx >= 0 ? nodes[idx].max_subtree_score
                  : std::numeric_limits<int>::min();
}

void ReadCoordinateRMQAVL::update_node_metadata(int32_t idx) {
  if (idx >= 0) {
    Node &node = nodes[idx];
    node.height = 1 + std::max(get_height(node.left), get_height(node.right));
    node.max_subtree_score = node.score_with_penalty;
    if (node.left >= 0) {
      node.max_subtree_score =
          std::max(node.max_subtree_score, nodes[node.left].max_subtree_score);
    }
    if (node.right >= 0) {
      node.max_subtree_score =
          std::max(node.max_subtree_score, nodes[node.right].max_subtree_score);
    }
  }
}

int32_t ReadCoordinateRMQAVL::rotate_right(int32_t y_idx) {
  int32_t x_idx = nodes[y_idx].left;
  nodes[y_idx].left = nodes[x_idx].right;
  update_node_metadata(y_idx);
  nodes[x_idx].right = y_idx;
  update_node_metadata(x_idx);
  return x_idx;
}

int32_t ReadCoordinateRMQAVL::rotate_left(int32_t x_idx) {
  int32_t y_idx = nodes[x_idx].right;
  nodes[x_idx].right = nodes[y_idx].left;
  update_node_metadata(x_idx);
  nodes[y_idx].left = x_idx;
  update_node_metadata(y_idx);
  return y_idx;
}

int ReadCoordinateRMQAVL::get_balance_factor(int32_t idx) const {
  return idx >= 0 ? (get_height(nodes[idx].left) - get_height(nodes[idx].right))
                  : 0;
}

int32_t ReadCoordinateRMQAVL::insert_recursive(int32_t node_idx, size_t coordinate,
                                            int score, size_t item_index) {
  if (node_idx < 0) {
    int32_t new_idx = (int32_t)nodes.size();
    nodes.emplace_back(coordinate, score, item_index, gap_extend);
    return new_idx;
  }

  // Split read and write to guard against reallocation from the recursive
  // emplace_back.
  if (coordinate < nodes[node_idx].coordinate) {
    int32_t new_left =
        insert_recursive(nodes[node_idx].left, coordinate, score, item_index);
    nodes[node_idx].left = new_left;
  } else {
    // Duplicate coordinates go to the right.
    int32_t new_right =
        insert_recursive(nodes[node_idx].right, coordinate, score, item_index);
    nodes[node_idx].right = new_right;
  }

  update_node_metadata(node_idx);

  int balance = get_balance_factor(node_idx);

  // Left Left Case
  if (balance > 1 && coordinate < nodes[nodes[node_idx].left].coordinate) {
    return rotate_right(node_idx);
  }

  // Right Right Case
  if (balance < -1 && coordinate >= nodes[nodes[node_idx].right].coordinate) {
    return rotate_left(node_idx);
  }

  // Left Right Case
  if (balance > 1 && coordinate >= nodes[nodes[node_idx].left].coordinate) {
    nodes[node_idx].left = rotate_left(nodes[node_idx].left);
    return rotate_right(node_idx);
  }

  // Right Left Case
  if (balance < -1 && coordinate < nodes[nodes[node_idx].right].coordinate) {
    nodes[node_idx].right = rotate_right(nodes[node_idx].right);
    return rotate_left(node_idx);
  }

  return node_idx;
}

void ReadCoordinateRMQAVL::query_iterative(
    size_t min_coordinate, size_t max_coordinate, int k_candidates,
    size_t max_nodes_to_visit,
    std::priority_queue<Candidate, std::vector<Candidate>,
                        std::greater<Candidate>> &min_heap) const {
  std::vector<int32_t> stack;
  stack.reserve(64);
  if (root >= 0) {
    stack.push_back(root);
  }

  size_t nodes_visited = 0;

  while (!stack.empty()) {
    int32_t node_idx = stack.back();
    stack.pop_back();
    nodes_visited++;

    // Score pruning before loading the full node.
    if (min_heap.size() >= (size_t)k_candidates &&
        nodes[node_idx].max_subtree_score <=
            min_heap.top().score_with_penalty) {
      continue;
    }

    const Node &node = nodes[node_idx];

    // Process current node if coordinate is in range.
    if (node.coordinate >= min_coordinate &&
        node.coordinate <= max_coordinate) {
      if (min_heap.size() < (size_t)k_candidates) {
        min_heap.push({node.score, node.score_with_penalty, node.item_index});
      } else if (node.score_with_penalty > min_heap.top().score_with_penalty) {
        min_heap.pop();
        min_heap.push({node.score, node.score_with_penalty, node.item_index});
      }
    }

    bool visit_left = node.left >= 0 && node.coordinate > min_coordinate;
    bool visit_right = node.right >= 0 && node.coordinate <= max_coordinate;

    // Push lower-score subtree first so higher-score is popped first,
    // filling the heap faster and raising the pruning threshold sooner.
    if (visit_left && visit_right) {
      if (nodes[node.left].max_subtree_score >=
          nodes[node.right].max_subtree_score) {
        stack.push_back(node.right);
        stack.push_back(node.left);
      } else {
        stack.push_back(node.left);
        stack.push_back(node.right);
      }
    } else if (visit_left) {
      stack.push_back(node.left);
    } else if (visit_right) {
      stack.push_back(node.right);
    }
  }

}

} // namespace algorithms
} // namespace vg
