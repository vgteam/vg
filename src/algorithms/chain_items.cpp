/**
 * \file chain_items.cpp
 * Non-template function implementations for chaining pieces of a read-to-graph alignment.
 */


#include "chain_items.hpp"

namespace vg {
namespace algorithms {

using namespace std;

ostream& operator<<(ostream& out, const traced_score_t& value) {
    if (score_traits<traced_score_t>::source(value) == score_traits<traced_score_t>::nowhere()) {
        return out << score_traits<traced_score_t>::score(value) << " from nowhere";
    }
    return out << score_traits<traced_score_t>::score(value) << " from #" << score_traits<traced_score_t>::source(value);
}

int MatchAssumingChainingScorer::score_transition(size_t read_distance,
                                                  const size_t* graph_distance,
                                                  size_t max_gap_length) const {
                             
    
    if (read_distance == numeric_limits<size_t>::max()) {
        // Overlap in read, so not allowed.
        return std::numeric_limits<int>::min();
    }
    
    if (!graph_distance) {
        // Assume graph distance equals read distance.
        graph_distance = &read_distance;
    }
    
    if (*graph_distance == numeric_limits<size_t>::max()) {
        // No graph connection
        return std::numeric_limits<int>::min();
    }
    
    // Decide how much length changed
    size_t indel_length = (read_distance > *graph_distance) ? read_distance - *graph_distance : *graph_distance - read_distance;
    
    if (indel_length > max_gap_length) {
        // Don't allow an indel this long
        return std::numeric_limits<int>::min();
    }
    
    // Decide how much length stayed the same, and assume it is matches.
    size_t common_length = std::min(read_distance, *graph_distance);
    int score = common_length * this->match;
    
    if (indel_length > 0) {
        score -= this->gap_open;
        if (indel_length > 1) {
            score -= this->gap_extension * (indel_length - 1);
        }
    }
    return score;
}

int IndelOnlyChainingScorer::score_transition(size_t read_distance,
                                              const size_t* graph_distance,
                                              size_t max_gap_length) const {
                             
    
    if (read_distance == numeric_limits<size_t>::max()) {
        // Overlap in read, so not allowed.
        return std::numeric_limits<int>::min();
    }
    
    if (!graph_distance) {
        // Assume graph distance equals read distance.
        graph_distance = &read_distance;
    }
    
    if (*graph_distance == numeric_limits<size_t>::max()) {
        // No graph connection
        return std::numeric_limits<int>::min();
    }
    
    // Decide how much length changed
    size_t indel_length = (read_distance > *graph_distance) ? read_distance - *graph_distance : *graph_distance - read_distance;
    
    if (indel_length > max_gap_length) {
        // Don't allow an indel this long
        return std::numeric_limits<int>::min();
    }
    
    // Then charge for that indel
    int score = 0;
    if (indel_length > 0) {
        score -= this->gap_open;
        if (indel_length > 1) {
            score -= this->gap_extension * (indel_length - 1);
        }
    }
    return score;
}

}
}
