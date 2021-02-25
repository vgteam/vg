// SPDX-FileCopyrightText: 2014 Erik Garrison
//
// SPDX-License-Identifier: MIT

#ifndef VG_ALGORITHMS_SORTED_ID_RANGES_HPP_INCLUDED
#define VG_ALGORITHMS_SORTED_ID_RANGES_HPP_INCLUDED

/**
 * \file sorted_id_ranges.hpp
 *
 * Defines an algorithm to get the ranges of IDs covered by a HandleGraph.
 */

#include <vector>
#include <utility>

#include "../handle.hpp"

namespace vg {
namespace algorithms {

using namespace std;

/// Get a sorted list of inclusive ranges of IDs used in the given HandleGraph.
vector<pair<id_t, id_t>> sorted_id_ranges(const HandleGraph* graph);


}
}

#endif
