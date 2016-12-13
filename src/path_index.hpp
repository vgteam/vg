#ifndef VG_PATH_INDEX_H
#define VG_PATH_INDEX_H

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
 * and orientation, and the full reference string. Also knows about the lengths
 * of nodes on the path, and lets you iterate back and forth over it.
 */
struct PathIndex {
    /// Index from node ID to first position on the reference string and
    /// orientation it occurs there.
    std::map<int64_t, std::pair<size_t, bool>> by_id;
    
    /// Index from start position on the reference to the side of the node that
    /// begins there. If it is a right side, the node occurs on the path in a
    /// reverse orientation.
    std::map<size_t, vg::NodeSide> by_start;
    
    /// This, combined with by_start, gets us the length of every node on the
    /// indexed path.
    size_t last_node_length;
    
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
    
    /// Find what node and orientation covers a position. The position must not
    /// be greater than the path length.
    NodeSide at_position(size_t position) const;
    
    /// We keep iterators to node occurrences along the ref path.
    using iterator = std::map<size_t, vg::NodeSide>::const_iterator;
    
    /// Get the iterator to the first node occurrence on the indexed path.
    iterator begin() const;
    /// Get the iterator to the last node occurrence on the indexed path.
    iterator end() const;
    
    /// Find the iterator at the given position along the ref path. The position
    /// must not be greater than the path length.
    iterator find_position(size_t position) const;
    
    /// Get the length of the node occurrence on the path represented by this
    /// iterator.
    size_t node_length(const iterator& here) const;
};

}
 
#endif
