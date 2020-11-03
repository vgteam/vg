#ifndef VG_PATH_INDEX_HPP_INCLUDED
#define VG_PATH_INDEX_HPP_INCLUDED

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
#include "handle.hpp"

namespace vg {

using namespace std;

/**
 * Holds indexes of the reference in a graph: position to node, node to position
 * and orientation, and the full reference string. Also knows about the lengths
 * of nodes on the path, and lets you iterate back and forth over it.
 */
struct PathIndex {
    /// Index from node ID to first position on the reference string and
    /// orientation it occurs there.
    map<int64_t, pair<size_t, bool>> by_id;
    
    /// Index from start position on the reference to the side of the node that
    /// begins there. If it is a right side, the node occurs on the path in a
    /// reverse orientation.
    map<size_t, NodeSide> by_start;
    
    /// The actual sequence of the path, if desired.
    std::string sequence;
    
    /// Index from Mapping pointers in a VG Paths object to their actual
    /// positions along their paths. Pointers may dangle if the vg graph
    /// changes the path.
    map<const mapping_t*, size_t> mapping_positions;
    
    /// Index just a path
    PathIndex(const Path& path);
    
    /// Index just a list of mappings
    PathIndex(const list<mapping_t>& mappings);
    
    /// Index a list of mappings embedded in the given vg's Paths object, and
    /// pull sequence from the given vg.
    PathIndex(const list<mapping_t>& mappings, VG& vg);
    
    /// Index a path and pull sequence from a graph.
    PathIndex(const Path& path, const HandleGraph& graph);
    
    /// Make a PathIndex from a path in a graph
    PathIndex(VG& vg, const string& path_name, bool extract_sequence = false);
    
    /// Make a PathIndex from a path in an indexed graph
    PathIndex(const PathHandleGraph& graph, const string& path_name, bool extract_sequence = false);
    
    /// Rebuild the mapping positions map by tracing all the paths in the given
    /// graph. TODO: We ought to move this functionality to the Paths object and
    /// make it use a good datastructure instead of brute force.
    void update_mapping_positions(VG& vg, const string& path_name);
    
    /// Find what node and orientation covers a position. The position must not
    /// be greater than the path length.
    NodeSide at_position(size_t position) const;

    /// Check whether a node is on the reference path.
    bool path_contains_node(int64_t node_id) const ;
    
    /// Check whether a node is on the reference path in a given path-relative orientation.
    bool path_contains_node_in_orientation(int64_t node_id, bool is_reverse) const;
    
    /// Return two flags for if the path contains the given node in forward and
    /// reverse orientation.
    pair<bool, bool> get_contained_orientations(int64_t node_id) const;
    
    /// We keep iterators to node occurrences along the ref path.
    using iterator = map<size_t, vg::NodeSide>::const_iterator;
    
    /// Find the first occurrence of the given node in the given orientation
    iterator find_in_orientation(int64_t node_id, bool is_reverse) const;
    
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
    
    /// Given an end-exclusive range on the path, round outward to the nearest
    /// node boundary positions.
    pair<size_t, size_t> round_outward(size_t start, size_t past_end) const;
    
    /**
     * Update the index to reflect the changes described by a Translation.
     * References to nodes along the "from" path are changed to references to
     * nodes along the "to" path. The translation must contain two paths of
     * equal length, containing only matches. The translation must only divide
     * nodes; it may not join nodes together. The translation must fully account
     * for each old node that it touches (it can't translate only part of a
     * node). The translation may not re-use the ID from one original node for a
     * piece of a different original node. All the Mappings in the Translation
     * must have Edits.
     */
    void apply_translation(const Translation& translation);
    
    /**
     * Update the index to reflect the changes described by the given collection
     * of Translations. These translations are expected to be in the format
     * produced by VG::edit() which is one to Mapping per translation. The
     * vector may include both forward and reverse versions of each to node, and
     * may also include translations mapping nodes that did not change to
     * themselves.
     */
    void apply_translations(const vector<Translation>& translations);
    
protected:

    /// This, combined with by_start, gets us the length of every node on the
    /// indexed path.
    size_t last_node_length;
    
    /// This holds all the places that a particular node occurs, in order.
    /// TODO: use this to replace by_id
    map<id_t, vector<iterator>> node_occurrences;
    
    /// Convert a Translation that partitions old nodes into a map from old node
    /// ID to the Mappings that replace it in its forward orientation.
    map<id_t, vector<Mapping>> parse_translation(const Translation& translation);
    
    /// Given an iterator into by_start, replace the occurrence of the node
    /// there with occurrences of the nodes given in the vector of mappings,
    /// which partition the forward strand of the node being replaced.
    void replace_occurrence(iterator to_replace, const vector<Mapping>& replacements);
    
    
    
};

}
 
#endif
