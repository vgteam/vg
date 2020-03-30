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
#include <functional>
#include "vg.hpp"
#include "tinygfa.hpp"

namespace vg {

using namespace std;

/**
 * Import the given GFA file into the given (empty) VG. If only_perfect_match
 * is set, only completely-M-operation CIGARs of length <= node length for
 * edge overlaps will be used, and sequence differences will be resolved
 * arbitrarily. Otherwise, CIGAR strings will be respected, files containing
 * alingments using more bases than the sequences will be rejected, and
 * mismatches in CIGAR-M sequences will form bubbles.
 *
 * Returns true if the import was successful. Returns false if the import
 * failed because the GFA file is invalid. Throws an error if the import failed
 * because of an apparent bug in the import code, or if the GFA tries to do
 * something that might be technically valid but which we don't know how to
 * interpret.
 */
bool gfa_to_graph(istream& in, VG* graph, bool only_perfect_match = false);

/// Export the given VG graph to the given GFA file.
void graph_to_gfa(const VG* graph, ostream& out);


}

#endif
