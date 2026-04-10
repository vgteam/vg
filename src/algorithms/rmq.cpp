#include "rmq.hpp"
#include <algorithm>
#include <limits>

namespace vg {
namespace algorithms {

ReadCoordinateRMQ::Node::Node(size_t coord, int s, size_t idx, int gap_extend, size_t chain_coord)
    : coordinate(coord), score(s),
      score_with_penalty(s + (int)(0.5 * gap_extend * (coord + chain_coord))),
      item_index(idx), max_subtree_score(score_with_penalty), height(1) {}

ReadCoordinateRMQ::ReadCoordinateRMQ(int gap_extend)
    : root(nullptr), gap_extend(gap_extend) {}

ReadCoordinateRMQ::~ReadCoordinateRMQ() {}

void ReadCoordinateRMQ::insert(size_t linear_coordinate, int partial_score,
                               size_t item_index, size_t chain_coordinate) {
  root = insert_recursive(std::move(root), linear_coordinate, partial_score,
                          item_index, chain_coordinate);
}

std::vector<size_t>
ReadCoordinateRMQ::query_top_k_candidates(size_t current_coordinate,
                                          size_t max_lookback_distance,
                                          int k_candidates) const {
  if (k_candidates <= 0 || !root) {
    return {};
  }

  size_t min_coordinate = (current_coordinate > max_lookback_distance)
                              ? (current_coordinate - max_lookback_distance)
                              : 0;
  size_t max_coordinate = current_coordinate;

  std::priority_queue<Candidate, std::vector<Candidate>,
                      std::greater<Candidate>>
      min_heap;
  query_recursive(root, min_coordinate, max_coordinate, k_candidates, min_heap);

  std::vector<size_t> results;
  results.reserve(min_heap.size());
  while (!min_heap.empty()) {
    results.push_back(min_heap.top().item_index);
    min_heap.pop();
  }
  std::reverse(results.begin(), results.end());
  return results;
}

int ReadCoordinateRMQ::get_height(const std::unique_ptr<Node> &node) const {
  return node ? node->height : 0;
}

int ReadCoordinateRMQ::get_max_subtree_score(
    const std::unique_ptr<Node> &node) const {
  return node ? node->max_subtree_score : std::numeric_limits<int>::min();
}

void ReadCoordinateRMQ::update_node_metadata(std::unique_ptr<Node> &node) {
  if (node) {
    node->height =
        1 + std::max(get_height(node->left), get_height(node->right));
    node->max_subtree_score = node->score_with_penalty;
    if (node->left) {
      node->max_subtree_score =
          std::max(node->max_subtree_score, node->left->max_subtree_score);
    }
    if (node->right) {
      node->max_subtree_score =
          std::max(node->max_subtree_score, node->right->max_subtree_score);
    }
  }
}

std::unique_ptr<ReadCoordinateRMQ::Node>
ReadCoordinateRMQ::rotate_right(std::unique_ptr<Node> y) {
  std::unique_ptr<Node> x = std::move(y->left);
  y->left = std::move(x->right);
  update_node_metadata(y);
  x->right = std::move(y);
  update_node_metadata(x);
  return x;
}

std::unique_ptr<ReadCoordinateRMQ::Node>
ReadCoordinateRMQ::rotate_left(std::unique_ptr<Node> x) {
  std::unique_ptr<Node> y = std::move(x->right);
  x->right = std::move(y->left);
  update_node_metadata(x);
  y->left = std::move(x);
  update_node_metadata(y);
  return y;
}


int ReadCoordinateRMQ::get_balance_factor(
    const std::unique_ptr<Node> &node) const {
  return node ? (get_height(node->left) - get_height(node->right)) : 0;
}

std::unique_ptr<ReadCoordinateRMQ::Node>
ReadCoordinateRMQ::insert_recursive(std::unique_ptr<Node> node,
                                    size_t coordinate, int score,
                                    size_t item_index, size_t chain_coordinate) {
  if (!node) {
    return std::unique_ptr<Node>(new Node(coordinate, score, item_index, gap_extend, chain_coordinate));
  }

  if (coordinate < node->coordinate) {
    node->left =
        insert_recursive(std::move(node->left), coordinate, score, item_index, chain_coordinate);
  } else {
    // We allow duplicate coordinates, they go to the right.
    node->right =
        insert_recursive(std::move(node->right), coordinate, score, item_index, chain_coordinate);
  }

  update_node_metadata(node);

  int balance = get_balance_factor(node);

  // Left Left Case
  if (balance > 1 && coordinate < node->left->coordinate) {
    return rotate_right(std::move(node));
  }

  // Right Right Case
  if (balance < -1 && coordinate >= node->right->coordinate) {
    return rotate_left(std::move(node));
  }

  // Left Right Case
  if (balance > 1 && coordinate >= node->left->coordinate) {
    node->left = rotate_left(std::move(node->left));
    return rotate_right(std::move(node));
  }

  // Right Left Case
  if (balance < -1 && coordinate < node->right->coordinate) {
    node->right = rotate_right(std::move(node->right));
    return rotate_left(std::move(node));
  }

  return node;
}

void ReadCoordinateRMQ::query_recursive(
    const std::unique_ptr<Node> &node, size_t min_coordinate,
    size_t max_coordinate, int k_candidates,
    std::priority_queue<Candidate, std::vector<Candidate>,
                        std::greater<Candidate>> &min_heap) const {
  if (!node) {
    return;
  }

  // Pruning: if the max possible score in this subtree is already worse than
  // the worst candidate we have, skip.
  if (min_heap.size() >= (size_t)k_candidates &&
      node->max_subtree_score <= min_heap.top().score_with_penalty) {
    return;
  }

  // Check if the current node is in range.
  if (node->coordinate >= min_coordinate &&
      node->coordinate <= max_coordinate) {
    if (min_heap.size() < (size_t)k_candidates) {
      min_heap.push({node->score, node->score_with_penalty, node->item_index});
    } else if (node->score_with_penalty > min_heap.top().score_with_penalty) {
      min_heap.pop();
      min_heap.push({node->score, node->score_with_penalty, node->item_index});
    }
  }

  // If searching left is potentially useful
  if (node->left && node->coordinate > min_coordinate) {
    query_recursive(node->left, min_coordinate, max_coordinate, k_candidates,
                    min_heap);
  }

  // If searching right is potentially useful
  if (node->right && node->coordinate <= max_coordinate) {
    query_recursive(node->right, min_coordinate, max_coordinate, k_candidates,
                    min_heap);
  }
}

} // namespace algorithms
} // namespace vg

