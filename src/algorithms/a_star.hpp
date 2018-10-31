#ifndef VG_ALGORITHMS_A_STAR_HPP_INCLUDED
#define VG_ALGORITHMS_A_STAR_HPP_INCLUDED

/**
 * \file a_star.hpp
 *
 * Defines an implementation of the A* algorithm.
 */

#include "../handle.hpp"
#include "topological_sort.hpp"

#include <unordered_map>
#include <vector>

namespace vg {
namespace algorithms {

using namespace std;

    /// Implements the A* heuristic-guided search algorithm. Returns the path from pos_1 to
    /// pos_2 that is either minimal or maximal length, according to the parameters. Allows
    /// an extremal distance beyond which the algorithm will cease looking for paths (this
    /// should be a large value when looking for minimal paths and a small value when looking
    /// for maximum paths). If there is no path between the positions, or none within the
    /// extremal length, an empty vector will be returned.
    template<class DistHeuristic>
    vector<handle_t> a_star(const HandleGraph* graph,
                            const pos_t& pos_1, const pos_t& pos_2,
                            const DistHeuristic& dist_heuristic,
                            bool find_min = true,
                            int64_t extremal_distance = numeric_limits<int64_t>::max());
}
}

#endif
