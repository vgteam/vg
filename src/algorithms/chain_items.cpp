/**
 * \file chain_items.cpp
 * Non-template function implementations for chaining pieces of a read-to-graph alignment.
 */


#include "chain_items.hpp"

#include <handlegraph/algorithms/dijkstra.hpp>

namespace vg {
namespace algorithms {

using namespace std;

ostream& operator<<(ostream& out, const TracedScore& value) {
    if (value.source == TracedScore::nowhere()) {
        return out << value.score << " from nowhere";
    }
    return out << value.score << " from #" << value.source;
}


void TracedScore::max_in(const vector<TracedScore>& options, size_t option_number) {
    auto& option = options[option_number];
    if (option.score > this->score || this->source == nowhere()) {
        // This is the new winner.
        this->score = option.score;
        this->source = option_number;
    }
}

TracedScore TracedScore::score_from(const vector<TracedScore>& options, size_t option_number) {
    TracedScore got = options[option_number];
    got.source = option_number;
    return got;
}

TracedScore TracedScore::add_points(int adjustment) const {
    return {this->score + adjustment, this->source};
}

}
}
