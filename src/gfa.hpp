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
#include <unordered_map>
#include "vg.hpp"
#include "tinygfa.hpp"

namespace vg {

using namespace std;

/// Export the given VG graph to the given GFA file.
void graph_to_gfa(const VG* graph, ostream& out);


}

#endif
