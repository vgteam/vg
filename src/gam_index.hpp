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
 * You will get non-contiguous virtual offset ranges for a node ID range when
 * some reads run into the range from the left, then reads that start later
 * don't, and then reads that start even later do again.
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
 * Unmapped reads are considered to visit node ID 0. The maximum and minimum
 * id_t values are used as sentinels, so they can't be real nodes.
 *
 * All find operations are thread-safe with respect to each other. Simultaneous
 * adds or finds and ads are prohibited.
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
    
    // Like the XG we support versioning.
    
    /// What's the maximum GAM index version number we can read with this code?
    const static uint32_t MAX_INPUT_VERSION = 1;
    /// What's the version we serialize?
    const static uint32_t OUTPUT_VERSION = 1;
    /// What magic value do we embed in the compressed gam index data?
    const static string MAGIC_BYTES;
    
    
    ///////////////////
    // Top-level Alignment-based interface
    ///////////////////
    
    /// Call the given callback with all Alignments in the index that visit the given node.
    void find(cursor_t& cursor, id_t node_id, const function<void(const Alignment&)> handle_result) const;
    
    /// Call the given callback with all Alignments in the index that visit a node in the given inclusive range.
    void find(cursor_t& cursor, id_t min_node, id_t max_node, const function<void(const Alignment&)> handle_result) const;
    
    /// Call the given callback with all the Alignments in the index that visit
    /// a node in any of the given sorted, coalesced inclusive ranges.
    /// Emits each alignment at most once.
    /// If only_fully_contained is set, only Alignments where *all* the mappings are to nodes in one of the ranges will match.
    void find(cursor_t& cursor, const vector<pair<id_t, id_t>>& ranges, const function<void(const Alignment&)> handle_result,
        bool only_fully_contained = false) const;
    
    /// Given a cursor at the beginning of a sorted, readable file, index the file.
    void index(cursor_t& cursor);
    
    /// Add a group articulated as a vector of alignments, between the given virtual offsets.
    /// Must be called in virtual offset order for successive groups.
    void add_group(const vector<Alignment>& alns, int64_t virtual_start, int64_t virtual_past_end);
    
    ///////////////////
    // Lower-level virtual-offset-based interface
    ///////////////////
    
    // Note that retrieving all the runs overlapping a node ID or node ID range
    // isn't possible. We can use the index to look up addresses to start at,
    // but the only way to know when to stop scanning groups is when you find a
    // group in the file with a minimum node ID that is too large. Then you
    // know to jump to the next start address.
    
    /// Find all the ranges of run virtual offsets from the first position that
    /// might be relevant for the given node ID to the ends of all the bins it
    /// is in. Trims ranges by the linear index on the low end, and returns a
    /// series of potentially abutting but non-overlapping virtual offset
    /// ranges. Does not stop early (because it has no access to the actual
    /// reads to tell when it should stop looking at runs in a bin). So you
    /// will get ranges covering all runs in a bin that follow the runs you are
    /// interested in as well.
    vector<pair<int64_t, int64_t>> find(id_t node_id) const;
    
    /// Find all the ranges of run virtual offsets to check for reads visiting
    /// the given node ID. Relies on a scanning callback, which will be called
    /// repeatedly with the start and past-the-end virtual offsets of runs
    /// which may contain groups touching the given node ID. When called, the
    /// callback should scan the run and return either true if it wants the
    /// next run, or false if it encountered a group with an out-of-range start
    /// and wants to stop iteration. Runs will be emitted in order, and
    /// truncated on the left to either the appropriate lower bound from the
    /// linear index, or the past-the-end of the previous run scanned.
    void find(id_t node_id, const function<bool(int64_t, int64_t)> scan_callback) const;
    
    // Find all the ranges of run virtual offsets to check for reads visiting
    /// the given inclusive node ID range. Relies on a scanning callback, which
    /// will be called repeatedly with the start and past-the-end virtual
    /// offsets of runs which may contain groups touching the given node ID.
    /// When called, the callback should scan the run and return either true if
    /// it wants the next run, or false if it encountered a group with an
    /// out-of-range start and wants to stop iteration. Runs will be emitted in
    /// order, and truncated on the left to either the appropriate lower bound
    /// from the linear index, or the past-the-end of the previous run scanned.
    void find(id_t min_node, id_t max_node, const function<bool(int64_t, int64_t)> scan_callback) const;
    
    /// Add a group into the index, based on its minimum and maximum
    /// (inclusive) used node IDs. Must be called for all groups in virtual
    /// offset order.
    void add_group(id_t min_id, id_t max_id, int64_t virtual_start, int64_t virtual_past_end);
    
    ///////////////////
    // Lowest-level functions for thinking about bins and windows.
    ///////////////////
    
    /// Compute the bins, from most to least specific, that a node ID occurs in.
    static vector<bin_t> bins_of_id(id_t id);
    
    /// Compute the bins, from most to least specific, that any of the node IDs
    /// in the given inclusive range occur in. There may be multiple bins at a
    /// given level of specificity; they will appear in numerical order.
    static vector<bin_t> bins_of_range(id_t min_id, id_t max_id);
    
    /// Get the most specific bin that contains both of the given node IDs.
    static bin_t common_bin(id_t a, id_t b);
    
    /// Get the linear index window that the given node ID falls in. The window
    /// range for a group is its min nodes' window through its max node's
    /// window.
    static window_t window_of_id(id_t id);
    
    
protected:
    
    // How many bits of a node ID do we truncate to get its linear index window?
    const static size_t WINDOW_SHIFT = 8;
    
    /// Maps from bin number to all the ranges of virtual offsets, in order, for runs that land in the given bin.
    /// A run lands in a bin if that bin is the most specific bin that includes both its lowest and highest nodes it uses.
    unordered_map<bin_t, vector<pair<int64_t, int64_t>>> bin_to_ranges;
    
    /// Maps from linear index window to the virtual offset of the first group
    /// that overlaps that window (taking the group as a min-to-max node
    /// range). If you are looking for reads that visit a node, they can't
    /// possibly occur in a group before the first offset stored for the node's
    /// window (or any greater window). TODO: Should we make this a vector
    /// instead and hope nobody uses high/sparse node IDs?
    map<window_t, int64_t> window_to_start;
    
    /// What was the minimum node ID of the last group added?
    /// If this isn't strictly increasing, we're trying to idnex data that is not sorted.
    id_t last_group_min_id = numeric_limits<id_t>::min();
    
};

}
 
#endif


