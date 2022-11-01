/**
 * \file chain_items.cpp
 * Non-template function implementations for chaining pieces of a read-to-graph alignment.
 */


#include "chain_items.hpp"

#include <handlegraph/algorithms/dijkstra.hpp>

namespace vg {
namespace algorithms {

using namespace std;

ostream& operator<<(ostream& out, const traced_score_t& value) {
    if (score_traits<traced_score_t>::source(value) == score_traits<traced_score_t>::nowhere()) {
        return out << score_traits<traced_score_t>::score(value) << " from nowhere";
    }
    return out << score_traits<traced_score_t>::score(value) << " from #" << score_traits<traced_score_t>::source(value);
}

}
}
