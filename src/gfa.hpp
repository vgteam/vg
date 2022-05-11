#ifndef VG_GFA_HPP_INCLUDED
#define VG_GFA_HPP_INCLUDED

/**
 * \file gfa.hpp
 *
 * Defines GFA I/O algorithms for PathHandleGraphs graphs.
 *
 * Includes an algorithm for converting from GFA, including non-perfect-match
 * edge overlaps and edges that specify containment of one node in another, to
 * a blunt-ended VG.
 */

#include "handle.hpp"

namespace vg {

using namespace std;

/// Export the given VG graph to the given GFA file.
/// Express paths mentioned in rgfa_paths as rGFA.
/// If rgfa_pline is set, also express them as dedicated lines.
/// If use_w_lines is set, reference and haplotype paths will use W lines instead of P lines.
void graph_to_gfa(const PathHandleGraph* graph, ostream& out,
                  const set<string>& rgfa_paths = {},
                  bool rgfa_pline = false,
                  bool use_w_lines = true);

}

#endif
