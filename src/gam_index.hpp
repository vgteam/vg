#ifndef VG_GAM_INDEX_HPP_INCLUDED
#define VG_GAM_INDEX_HPP_INCLUDED

/**
 * \file gam_index.hpp
 * Contains the GAMIndex class, which allows retrieving reads from a (blocked) BAM file for certain queries.
 */
 
#include <iostream>
#include <vector>
#include <set>
#include <unordered_map>
#include <type_traits>

#include "types.hpp"
#include "vg.pb.h"
#include "stream.hpp"

namespace vg {

using namespace std;

/**
 * 
 * An index for a node-ID-sorted GAM file. Reads are sorted by lowest visited
 * node ID, then by highest visited node ID.
 *
 * Works on a BAI-like concept of bins partitioning node ID space.
 *
 * GAM files are serialized as count-prefixed groups of reads, which are the
 * smallest unit that can be deserialized.
 *
 * Every *group* of reads gets assigned to a bin which is the longest bin that
 * completely contains the ID range used in the group.
 *
 * We define *runs* of adjacent groups of reads which have the same bin, which
 * are the basic subject of the index.
 *
 * We then store an index from bin to the virtual offset ranges (start and
 * past-the-end), in order, of runs that are assigned to the bin.
 *
 * We also have a BAI-style linear index, mapping from tiling windows in node
 * ID space to the lowest virtual offset of a group that overlaps the window.
 *
 * The bin structure is that we partition all of node ID space into bins of
 * power-of-2 size, starting with size 2 nodes. We number the bins such that 0
 * is the whole-ID-space bin, divided into 1 and 2, then into 3, 4, 5, and 6,
 * and so on.
 *
 * The tiling windows are just the node IDs down-shifted by a few bits.
 *
 * Unmapped reads are considered to visit node ID 0.
 *
 */
class GAMIndex {
public:
    GAMIndex() = default;
    
    // Methods that actually go get reads for you are going to need a cursor on an open, seekable GAM file.
    using cursor_t = stream::ProtobufIterator<Alignment>;
    
    // Bins are identified of unsigned integers of the same width as node IDs.
    using bin_t = make_unsigned<id_t>::type;
    
    // So are windows, but we give them their own semantic type
    using window_t = make_unsigned<id_t>::type;
    
    /// Load a GAMIndex from a file.
    /// File holds the index, not the GAM.
    void load(istream& from);
    
    /// Save a GAMIndex to a file.
    void save(ostream& to) const;
    
    
    ///////////////////
    // Top-level Alignment-based interface
    ///////////////////
    
    /// Call the given callback with all Alignments in the index that visit the given node.
    void find(cursor_t& cursor, id_t node_id, const function<void(const Alignment&)>& handle_result) const;
    
    /// Call the given callback with all Alignments in the index that visit a node in the given range.
    void find(cursor_t& cursor, id_t first_id, id_t past_last_id, const function<void(const Alignment&)>& handle_result) const;
    
    /// Given a cursor at the beginning of a sorted, readable file, index the file.
    void index(cursor_t& cursor);
    
    /// Add a group articulated as a vector of alignments, between the given virtual offsets.
    /// Must be called in virtual offset order for successive groups.
    void add_group(const vector<Alignment>& alns, int64_t virtual_start, int64_t virtual_past_end);
    
    ///////////////////
    // Lower-level virtual-offset-based interface
    ///////////////////
    
    /// Find all the ranges of run virtual offsets to check for reads visiting the given node ID.
    /// Trims ranges by the linear index.
    vector<pair<int64_t, int64_t>> find(id_t node_id) const;
    
    /// Add a group into the index, based on its minimum and maximum used node
    /// IDs. Must be called for all groups in virtual offset order.
    void add_group(id_t min_id, id_t max_id, int64_t virtual_start, int64_t virtual_past_end);
    
    ///////////////////
    // Lowest-level functions for thinking about bins and windows.
    ///////////////////
    
    /// Compute the bins, from most to least specific, that a node ID occurs in.
    static vector<bin_t> bins_of_id(id_t id);
    
    /// Get the most specific bin that contains both of the given node IDs.
    static bin_t common_bin(id_t a, id_t b);
    
    /// Get the linear index window that the given node ID falls in.
    static window_t window_of_id(id_t id);
    
    
protected:
    
    // How many bits of a node ID do we truncate to get its linear index window?
    const static size_t WINDOW_SHIFT = 8;
    
    /// Maps from bin number to all the ranges of virtual offsets, in order, for runs that land in the given bin.
    /// A run lands in a bin if that bin is the most specific bin that includes both its lowest and highest nodes it uses.
    unordered_map<bin_t, vector<pair<int64_t, int64_t>>> bin_to_ranges;
    
    /// Maps from linear index window to the virtual offset of the first group
    /// in that window. If you are looking for reads that visit a node, they
    /// can't possibly occur in a group before the first offset stored for the
    /// node's window. TODO: Should we make this a vector instead and hope
    /// nobody uses high/sparse node IDs?
    map<window_t, int64_t> window_to_start;
    
};

}
 
#endif


