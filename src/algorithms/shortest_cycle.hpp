#ifndef VG_ALGORITHMS_SHORTEST_CYCLE_HPP_INCLUDED
#define VG_ALGORITHMS_SHORTEST_CYCLE_HPP_INCLUDED

#include <unordered_map>

#include "../handle.hpp"

#include "structures/rank_pairing_heap.hpp"

namespace vg {
namespace algorithms {

using namespace std;
    
    /// Returns the length of the shortest cycle in the entire graph, or
    /// numeric_limits<size_t>::max() if no cycle exists.
    size_t shortest_cycle_length(const HandleGraph* graph);

    /// Returns the length of the shortest cycle containing the source node, or
    /// numeric_limits<size_t>::max() if no cycle exists.
    size_t shortest_cycle_length(const HandleGraph* graph, const handle_t& source);
}
}

#endif
