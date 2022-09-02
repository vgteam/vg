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

std::unique_ptr<LookbackStrategy::Problem> FlatLimitLookbackStrategy::setup_problem(size_t item_count) const {
     std::unique_ptr<LookbackStrategy::Problem> to_return;
     to_return.reset(new FlatLimitLookbackStrategy::Problem(*this));
     return std::move(to_return);
}

FlatLimitLookbackStrategy::Problem::Problem(const FlatLimitLookbackStrategy& parent) : strategy(parent) {
    // Nothing to do!
}

void FlatLimitLookbackStrategy::Problem::advance() {
    // Reset the per-destination counters.
    this->lookback_items_used = 0;
    this->lookback_good_items_used = 0;
    this->lookback_reachable_items_used = 0;
    this->best_achieved_score = 0;
}

LookbackStrategy::verdict_t FlatLimitLookbackStrategy::Problem::should_check(size_t item_a, size_t item_b, size_t read_distance, int item_a_score) {
    if (read_distance > this->strategy.lookback_bases) {
        // This is too far in the read. Everything else will be further, so stop.
        return LookbackStrategy::STOP;
    }
    
    // Count this item as a lookback item.
    // TODO: Harmonize terminology so this makes sense.
    this->lookback_items_used++;
    if (this->lookback_items_used > this->strategy.lookback_items) {
        // This is enough items already, so stop.
        return LookbackStrategy::STOP;
    }
    
    if (this->lookback_reachable_items_used >= this->strategy.lookback_reachable_items) {
        // We already found enough reachable previous items, so stop.
        return LookbackStrategy::STOP;
    }
    
    if (this->lookback_good_items_used >= this->strategy.lookback_good_items) {
        // We already found enough good previous items, so stop.
        return LookbackStrategy::STOP;
    }
    
    // Otherwise, do this item.
    return LookbackStrategy::CHECK;
}

void FlatLimitLookbackStrategy::Problem::did_check(size_t item_a, size_t item_b, size_t read_distance, const size_t* graph_distance, int transition_score, int achieved_score) {
    if (transition_score != std::numeric_limits<int>::min()) {
        // We considered the transition to be usable, so remember that this was a reachable item.
        this->lookback_reachable_items_used++;
    }
    if (achieved_score > this->best_achieved_score) {
        // This is a new best so count it as a "good" item.
        this->best_achieved_score = achieved_score;
        this->lookback_good_items_used++;
    }
}

}
}
