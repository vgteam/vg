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

/**
 * Serialize the given graph to the given output in GFA format.
 *
 * This effectively calls write_gfa_header(), write_gfa_nodes(),
 * write_gfa_edges(), and write_gfa_paths() with the given arguments.
 * See the documentation of those functions for further details.
 */
void graph_to_gfa(
    const PathHandleGraph* graph, ostream& out,
    const set<string>& rgfa_paths = {}, bool include_rgfa_paths = false,
    bool use_w_lines = true
);

/**
 * Write a GFA header line of a GFA 1.1 file for the given graph.
 *
 * Reference samples will be listed in the RS tag. Reference sense paths
 * included in rgfa_paths will not be included, unless include_rgfa_paths is
 * set.
 *
 * @note Pass write_gfa_paths() the same rgfa_paths and include_rgfa_paths
 * arguments, so that the RS tag matches the paths written.
 */
void write_gfa_header(
    const PathHandleGraph* graph, ostream& out,
    const set<string>& rgfa_paths = {}, bool include_rgfa_paths = false
);

/**
 * Write a GFA segment line for each node of the given graph.
 *
 * This also writes rGFA tags for nodes visited by paths named in rgfa_paths.
 *
 * @throws std::runtime_error if a path named in rgfa_paths visits a node in
 * reverse orientation.
 */
void write_gfa_nodes(
    const PathHandleGraph* graph, ostream& out,
    const set<string>& rgfa_paths = {}
);

/**
 * Write a GFA link line for each edge of the given graph.
 *
 * Each edge is written from a forward handle if either of its directions
 * allows that. If both directions start from a reverse handle, the edge is
 * written from the lower node ID. This keeps the output stable for tests.
 */
void write_gfa_edges(const HandleGraph* graph, ostream& out);

/**
 * Write the paths of the given graph as GFA lines, sorted by name.
 *
 * Paths included in rgfa_paths are skipped, unless include_rgfa_paths is set.
 * Each path is written either as a path (P) or a walk (W) line, depending
 * on use_w_lines. Generic paths and paths with a nonzero phase block are
 * always written as P lines.
 *
 * @note Pass write_gfa_header() the same rgfa_paths and rgfa_pline, so that
 * its RS tag matches the paths written.
 */
void write_gfa_paths(
    const PathHandleGraph* graph, ostream& out,
    const set<string>& rgfa_paths = {}, bool include_rgfa_paths = false,
    bool use_w_lines = true
);

}

#endif
