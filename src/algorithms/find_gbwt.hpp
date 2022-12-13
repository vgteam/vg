#ifndef VG_ALGORITHMS_FIND_GBWT_HPP_INCLUDED
#define VG_ALGORITHMS_FIND_GBWT_HPP_INCLUDED

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
 * Find the GBWT that is part of the given handle graph, if any exists.
 * Works on GBWTGraphs and GBZGraphs.
 * Returns null if no such GBWT exists.
 */
const gbwt::GBWT* find_gbwt(const HandleGraph* graph);

/**
 * Find a GBWT either by getting it from the given graph or loading it from the
 * given filename into the given unique_ptr.
 */
const gbwt::GBWT* find_gbwt(const HandleGraph* graph, std::unique_ptr<gbwt::GBWT>& holder, const std::string& filename);

}
}

#endif
