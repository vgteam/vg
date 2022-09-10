/**
 * \file chain_items.cpp
 * Non-template function implementations for chaining pieces of a read-to-graph alignment.
 */


#include "chain_items.hpp"

namespace vg {
namespace algorithms {

using namespace std;

#define debug_lookback

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

void LookbackStrategy::Problem::place_in_graph(size_t item, const pos_t& graph_pos, const SnarlDistanceIndex& distance_index, const HandleGraph& graph) {
    // Default implementation: no-op.
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

std::unique_ptr<LookbackStrategy::Problem> ExponentialLookbackStrategy::setup_problem(size_t item_count) const {
     std::unique_ptr<LookbackStrategy::Problem> to_return;
     to_return.reset(new ExponentialLookbackStrategy::Problem(*this));
     return std::move(to_return);
}

ExponentialLookbackStrategy::Problem::Problem(const ExponentialLookbackStrategy& parent) : strategy(parent) {
    // Nothing to do!
}

void ExponentialLookbackStrategy::Problem::advance() {
    // Reset the per-destination state.
    this->limit = this->strategy.initial_search_bases;
    this->best_transition_found = std::numeric_limits<int>::min();
    this->best_achieved_score = std::numeric_limits<int>::min();
    this->good_score_found = false;
}

LookbackStrategy::verdict_t ExponentialLookbackStrategy::Problem::should_check(size_t item_a, size_t item_b, size_t read_distance, int item_a_score) {
    if (read_distance > this->strategy.lookback_bases) {
        // This is further in the read than the real hard limit.
        return LookbackStrategy::STOP;
    } else if (read_distance > this->limit) {
        // This is further than we wanted to look.
        if (this->good_score_found) {
            // And we have something good already, so impose the limit
            return LookbackStrategy::STOP;
        } else {
            // But we still haven't found anything good, so raise the limit.
            this->limit *= this->strategy.scale_factor;
#ifdef debug_lookback
            std::cerr << "\t\tIncreased lookback limit to " << this->limit << " bp because best transition is " << this->best_transition_found << " so we now allow " << (this->strategy.min_good_transition_score_per_base * this->limit) << std::endl;
#endif
            return LookbackStrategy::CHECK;
        }
    } else {
        // We are close enough to not hit any limits
        return LookbackStrategy::CHECK;
    }
}

void ExponentialLookbackStrategy::Problem::did_check(size_t item_a, size_t item_b, size_t read_distance, const size_t* graph_distance, int transition_score, int achieved_score) {
    if (!graph_distance) {
        graph_distance = &read_distance;
    }
    this->best_transition_found = std::max(this->best_transition_found, transition_score);
    this->best_achieved_score = std::max(this->best_achieved_score, achieved_score);
    if (achieved_score > 0 && this->best_transition_found >= this->strategy.min_good_transition_score_per_base * std::max(read_distance, *graph_distance)) {
        // We found a jump that looks plausible given how far we have searched, so we can stop searching way past here.
        this->good_score_found = true;
    }
}


std::unique_ptr<LookbackStrategy::Problem> BucketLookbackStrategy::setup_problem(size_t item_count) const {
     std::unique_ptr<LookbackStrategy::Problem> to_return;
     to_return.reset(new BucketLookbackStrategy::Problem(*this));
     return std::move(to_return);
}

BucketLookbackStrategy::Problem::Problem(const BucketLookbackStrategy& parent) : strategy(parent) {
    // Nothing to do!
}

void BucketLookbackStrategy::Problem::advance() {
    // Reset the per-destination state.
    this->limit = this->strategy.initial_search_bases;
    this->best_transition_found = std::numeric_limits<int>::min();
    this->best_achieved_score = std::numeric_limits<int>::min();
    this->good_score_found = false;
}

void BucketLookbackStrategy::Problem::place_in_graph(size_t item, const pos_t& graph_pos, const SnarlDistanceIndex& distance_index, const HandleGraph& graph) {
    // We need to assign this item to a bucket
    for (size_t bucket = 0; bucket < bucket_heads.size(); bucket++) {
        auto& head = bucket_heads[bucket];
        size_t distance = distance_index.minimum_distance(
            id(head), is_rev(head), offset(head),
            id(graph_pos), is_rev(graph_pos), offset(graph_pos),
            false, &graph);
        if (distance != std::numeric_limits<size_t>::max() && distance < this->strategy.bucket_limit) {
            // This is the right bucket!
            bucket_sizes[bucket]++;
#ifdef debug_lookback
            std::cerr << "\t\tItem " << item << " is " << distance << " away from the head of bucket " << bucket << " so put it there to get size " << bucket_sizes[bucket] << std::endl;
#endif
            item_to_head.emplace(item, bucket);
            item_to_coordinate.emplace(item, distance);
            return;
        }
    }
    // If it's not reachable from the head of any existing bucket, give it a new bucket
#ifdef debug_lookback
    std::cerr << "\t\tItem " << item << " must start a new bucket " << bucket_heads.size() << std::endl;
#endif
    item_to_head.emplace(item, bucket_heads.size());
    item_to_coordinate.emplace(item, 0);
    bucket_heads.emplace_back(graph_pos);
    bucket_sizes.emplace_back(1);
}

LookbackStrategy::verdict_t BucketLookbackStrategy::Problem::should_check(size_t item_a, size_t item_b, size_t read_distance, int item_a_score) {
    if (read_distance > this->strategy.lookback_bases) {
        // This is too far in the read. Everything else will be further, so stop.
        return LookbackStrategy::STOP;
    } else if (read_distance > this->limit) {
        // This is further than we wanted to look.
        if (this->good_score_found) {
            // And we have something good already, so impose the limit
            return LookbackStrategy::STOP;
        } else {
            // But we still haven't found anything good, so raise the limit.
            this->limit *= this->strategy.scale_factor;
#ifdef debug_lookback
            std::cerr << "\t\tIncreased lookback limit to " << this->limit << " bp because best transition is " << this->best_transition_found << " so we now allow " << (this->strategy.min_good_transition_score_per_base * this->limit) << std::endl;
#endif
        }
    }
    // We might fall through after increasing the limit
    if (!bucket_heads.empty() && item_to_head.at(item_a) != item_to_head.at(item_b)) {
        // These items are in different buckets.
#ifdef debug_lookback
        std::cerr << "\t\tItem " << item_a << " is in bucket " << item_to_head.at(item_a) << " size " << bucket_sizes.at(item_to_head.at(item_a)) << " but item " << item_b << " is in bucket " << item_to_head.at(item_b) << " size " << bucket_sizes.at(item_to_head.at(item_b)) << std::endl;
#endif
        return LookbackStrategy::SKIP;
    } else {
        if (!bucket_heads.empty()) {
#ifdef debug_lookback
            std::cerr << "\t\tBoth items are in bucket " << item_to_head.at(item_a) << " size " << bucket_sizes.at(item_to_head.at(item_a)) << std::endl;
#endif
        }
        // Otherwise they are in the same bucket
        
        auto& coordinate_a = item_to_coordinate.at(item_a);
        auto& coordinate_b = item_to_coordinate.at(item_b);
        
#ifdef debug_lookback
        std::cerr << "\t\tBucket coordinates: " << coordinate_a << " vs " << coordinate_b << std::endl;
#endif
        
        if (coordinate_b > coordinate_a) {
            // We can infer a minimum minimum distance from the distances to the head.
            size_t min_min_distance = coordinate_b - coordinate_a;
            if (min_min_distance > read_distance) {
                // And since they must be further apart in the graph than the read, we can check how much further
                size_t inferred_indel = min_min_distance - read_distance;
#ifdef debug_lookback
                std::cerr << "\t\tItems must be at least " << inferred_indel << " further apart in graph than in read" << std::endl;
#endif
                if (inferred_indel > this->strategy.max_inferred_indel) {
                    return LookbackStrategy::SKIP;
                }
            }
        } else {
            // We can still look at the difference in bucket space but it doesn't really mean anything.
            size_t bucket_coordinate_difference = coordinate_a - coordinate_b;
            if (bucket_coordinate_difference > this->strategy.max_suspicious_bucket_coordinate_difference) {
#ifdef debug_lookback
                std::cerr << "\t\tItems seem suspiciously different distances from bucket head" << std::endl;
#endif
                return LookbackStrategy::SKIP;
            }
        }
        
        return LookbackStrategy::CHECK;
    }
}

void BucketLookbackStrategy::Problem::did_check(size_t item_a, size_t item_b, size_t read_distance, const size_t* graph_distance, int transition_score, int achieved_score) {
    if (!graph_distance) {
        graph_distance = &read_distance;
    }
    this->best_transition_found = std::max(this->best_transition_found, transition_score);
    this->best_achieved_score = std::max(this->best_achieved_score, achieved_score);
    if (achieved_score > 0 && this->best_transition_found >= this->strategy.min_good_transition_score_per_base * std::max(read_distance, *graph_distance)) {
        // We found a jump that looks plausible given how far we have searched, so we can stop searching way past here.
        this->good_score_found = true;
    }
}

}
}
