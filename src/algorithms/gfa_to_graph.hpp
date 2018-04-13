#ifndef VG_ALGORITHMS_GFA_TO_GRAPH_HPP_INCLUDED
#define VG_ALGORITHMS_GFA_TO_GRAPH_HPP_INCLUDED

/**
 * \file gfa_to_graph.hpp
 *
 * Defines an algorithm for converting from GFA, including non-perfect-match
 * edge overlaps and edges that specify containment of one node in another, to
 * a blunt-ended MutableHandleGraph.
 */

#include "../handle.hpp"

#include <gfakluge.hpp>

namespace vg {
namespace algorithms {

using namespace std;
using namespace gfak;


/// Import the given GFA file into the given (empty) MutableHandleGraph. If
/// only_perfect_match is set, only completely-M-operation CIGARs for edge
/// overlaps will be used, and sequence differences will be resolved
/// arbitrarily. Otherwise, CIGAR strings will be respected, and mismatches in
/// CIGAR-M sequences will for bubbles.
void gfa_to_graph(GFAKluge& gfa, MutableHandleGraph* graph, bool only_perfect_match = false);


}
}

#endif
