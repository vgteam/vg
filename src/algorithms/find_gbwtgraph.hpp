#ifndef VG_ALGORITHMS_FIND_GBWTGRAPH_HPP_INCLUDED
#define VG_ALGORITHMS_FIND_GBWTGRAPH_HPP_INCLUDED

/**
 * \file find_gbwtgraph.hpp
 *
 * Defines an algorithm for finding the GBWTGraph associated with a handle graph, if any.
 */

#include "../handle.hpp"
#include <gbwtgraph/gbwtgraph.h>

namespace vg {
namespace algorithms {
using namespace std;

/**
 * Find the GBWTGraph that is part of the given handle graph, if any exists.
 * Works on GBWTGraphs and GBZGraphs.
 * Returns null if no such GBWTGraph exists.
 */
const gbwtgraph::GBWTGraph* find_gbwtgraph(const HandleGraph* graph);

}
}

#endif
