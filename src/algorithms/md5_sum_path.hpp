#ifndef VG_ALGORITHMS_MD5_SUM_PATH_HPP_INCLUDED
#define VG_ALGORITHMS_MD5_SUM_PATH_HPP_INCLUDED
/**
 * \file md5_sum_path.hpp
 * Algorithms for MD5-summing paths
 */

#include "handle.hpp"
#include <string>
#include <utility>

namespace vg {
namespace algorithms {

/// Get the Samtools-compatible hex MD5 sum and the length of a path in a graph, in one pass.
std::pair<std::string, size_t> md5_sum_path_with_length(const PathHandleGraph& graph, const path_handle_t& path);

}
}

#endif


