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
    if (TracedScore::source(value) == TracedScore::nowhere()) {
        return out << TracedScore::score(value) << " from nowhere";
    }
    return out << TracedScore::score(value) << " from #" << TracedScore::source(value);
}


void TracedScore::max_in(TracedScore& dest, const vector<TracedScore>& options, size_t option_number) {
    auto& option = options[option_number];
    if (score(option) > score(dest) || source(dest) == nowhere()) {
        // This is the new winner.
        score(dest) = score(option);
        source(dest) = option_number;
    }
}

TracedScore TracedScore::score_from(const vector<TracedScore>& options, size_t option_number) {
    TracedScore got = options[option_number];
    source(got) = option_number;
    return got;
}

TracedScore TracedScore::add_points(const TracedScore& s, int adjustment) {
    return annotate(score(s) + adjustment, source(s));
}

}
}
