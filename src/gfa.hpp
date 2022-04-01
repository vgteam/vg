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
void graph_to_gfa(const PathHandleGraph* graph, ostream& out,
                  const set<string>& rgfa_paths = {},
                  bool rgfa_pline = false,
                  bool use_w_lines = false);

}

#endif
