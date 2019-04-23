#ifndef VG_ALGORITHMS_REMOVE_HIGH_DEGREE_HPP_INCLUDED
#define VG_ALGORITHMS_REMOVE_HIGH_DEGREE_HPP_INCLUDED

/**
 * \file remove_high_degree.hpp
 *
 * Defines a process that removes high-degree nodes from a graph
 */

#include <vg/vg.pb.h>

#include "../handle.hpp"
#include <vector>

namespace vg {
namespace algorithms {

using namespace std;

void remove_high_degree_nodes(DeletableHandleGraph& g, int max_degree);

}
}

#endif
