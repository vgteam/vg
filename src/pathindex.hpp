#ifndef VG_PATHINDEX_H
#define VG_PATHINDEX_H

/** \file
 *
 * Provides an index for indexing an individual path for fast random access.
 * Stores all the mappings uncompressed in memory.
 *
 * Used for the reference path during VCF creation and interpretation.
 */
 
#include <map>
#include <utility>
#include <string>

#include "vg.hpp"
#include "xg.hpp"

namespace vg {

/**
 * Holds indexes of the reference in a graph: position to node, node to position
 * and orientation, and the full reference string.
 */
struct PathIndex {
    /// Index from node ID to first position on the reference string and
    /// orientation it occurs there.
    std::map<int64_t, std::pair<size_t, bool>> by_id;
    
    /// Index from start position on the reference to the side of the node that
    /// begins there. If it is a right side, the node occurs on the path in a
    /// reverse orientation.
    std::map<size_t, vg::NodeSide> by_start;
    
    /// The actual sequence of the path, if desired.
    std::string sequence;
    
    /// Index just a path
    PathIndex(const Path& path);
    
    /// Index a path and pull sequence from a VG graph.
    PathIndex(const Path& path, VG& vg);
    
    /// Index a path and pull sequence from an XG index.
    PathIndex(const Path& path, const xg::XG& vg);
    
    /// Make a PathIndex from a path in a graph
    PathIndex(VG& vg, const string& ref_path_name, bool extract_sequence = false);
    
    /// Make a PathIndex from a path in an indexed graph
    PathIndex(const xg::XG& index, const string& ref_path_name, bool extract_sequence = false);
    
    // Find what node and orientation covers a position
    NodeSide at_position(size_t position);
};

}
 
#endif
