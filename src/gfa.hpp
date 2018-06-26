#ifndef VG_GFA_HPP_INCLUDED
#define VG_GFA_HPP_INCLUDED

/**
 * \file gfa.hpp
 *
 * Defines GFA I/O algorithms for VG graphs.
 *
 * Includes an algorithm for converting from GFA, including non-perfect-match
 * edge overlaps and edges that specify containment of one node in another, to
 * a blunt-ended VG.
 */

#include "vg.hpp"

#include <gfakluge.hpp>

namespace vg {

using namespace std;
using namespace gfak;


/// Import the given GFA file into the given (empty) VG. If
/// only_perfect_match is set, only completely-M-operation CIGARs for edge
/// overlaps will be used, and sequence differences will be resolved
/// arbitrarily. Otherwise, CIGAR strings will be respected, and mismatches in
/// CIGAR-M sequences will for bubbles.
void gfa_to_graph(istream& in, VG* graph, bool only_perfect_match = false);
void gfa_to_graph(GFAKluge& gfa, VG* graph, bool only_perfect_match = false);

/// Export the given VG graph to the given GFA file.
void graph_to_gfa(const VG* graph, ostream& out);
void graph_to_gfa(const VG* graph, GFAKluge& gfa);


}

#endif
