#ifndef VG_STREAM_INDEX_HPP_INCLUDED
#define VG_STREAM_INDEX_HPP_INCLUDED

/**
 * \file stream_index.hpp
 * Contains the StreamIndex template, which allows lookup by relevant node ID in sorted stream.hpp-formatted files.
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
 * An index for a node-ID-sorted stream.hpp-formatted file, such as GAM or VG.
 *
 * Works on a BAI-like concept of bins partitioning node ID space.
 *
 * Files are serialized as count-prefixed groups of Protobuf messages. These
 * groups are the smallest unit that can be deserialized.
 *
 * Every *group* of messages gets assigned to a bin which is the longest bin
 * that completely contains the ID range used in the group.
 *
 * We define *runs* of adjacent groups which have the same bin, which are the
 * basic subject of the index.
 *
 * We then store an index from bin to the virtual offset ranges (start and
 * past-the-end), in order, of runs that are assigned to the bin.
 *
 * You will get non-contiguous virtual offset ranges for a node ID range when
 * some messages run into the range from the left, then messages that start
 * later don't, and then messages that start even later do again.
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
 * Messages that use no nodes (i.e. unmapped reads) are considered to visit
 * node ID 0. The maximum and minimum id_t values are used as sentinels, so
 * they can't be real nodes.
 *
 * All find operations are thread-safe with respect to each other. Simultaneous
 * adds or finds and adds are prohibited.
 *
 * Most of the basic index API doesn't depend on the message type. So we put it
 * in this base class and inherit form it in templates that provide the
 * high-level interface in terms of message instances.
 */
class StreamIndexBase {
public:
    StreamIndexBase() = default;
    
    // Bins are identified of unsigned integers of the same width as node IDs.
    using bin_t = make_unsigned<id_t>::type;
    
    // So are windows, but we give them their own semantic type
    using window_t = make_unsigned<id_t>::type;
    
    /// Load an index from a file.
    /// File holds the index, not the actual data being indexed.
    /// Index file format doesn't care what type of message is being indexed.
    void load(istream& from);
    
    /// Save an index to a file.
    void save(ostream& to) const;
    
    // Like the XG we support versioning.
    
    /// What's the maximum index version number we can read with this code?
    const static uint32_t MAX_INPUT_VERSION = 1;
    /// What's the version we serialize?
    const static uint32_t OUTPUT_VERSION = 1;
    /// What magic value do we embed in the compressed index data?
    /// TODO: Make this depend on type of message being indexed so we can't mix up index files.
    const static string MAGIC_BYTES;
    
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
    
    /// Return true if the given ID is in any of the sorted, coalesced, inclusive ranges in the vector, and false otherwise.
    /// TODO: Is repeated binary search on the ranges going to be better than an unordered_set of all the individual IDs?
    static bool is_in_range(const vector<pair<id_t, id_t>>& ranges, id_t id);

};

/**
 * An index that provides a higher-level API in terms of the actual messages
 * being indexed. This is the main entry point for users in most cases.
 *
 * All find operations are thread-safe with respect to each other. Simultaneous
 * adds or finds and adds are prohibited.
 *
 */
template<typename Message>
class StreamIndex : public StreamIndexBase {
public:
    StreamIndex() = default;
    
    // Methods that actually go get messages for you are going to need a cursor on an open, seekable data file.
    using cursor_t = stream::ProtobufIterator<Message>;
    
    ///////////////////
    // Top-level message-based interface
    ///////////////////
    
    /// Call the given callback with all messages in the index that visit the given node.
    void find(cursor_t& cursor, id_t node_id, const function<void(const Message&)> handle_result) const;
    
    /// Call the given callback with all messages in the index that visit a node in the given inclusive range.
    void find(cursor_t& cursor, id_t min_node, id_t max_node, const function<void(const Message&)> handle_result) const;
    
    /// Call the given callback with all the messages in the index that visit
    /// a node in any of the given sorted, coalesced inclusive ranges.
    /// Emits each message at most once.
    /// If only_fully_contained is set, only messages where *all* the involved nodes are in one of the ranges will match.
    void find(cursor_t& cursor, const vector<pair<id_t, id_t>>& ranges, const function<void(const Message&)> handle_result,
        bool only_fully_contained = false) const;
    
    /// Given a cursor at the beginning of a sorted, readable file, index the file.
    void index(cursor_t& cursor);
    
    /// Add a group articulated as a vector of messages, between the given virtual offsets.
    /// Must be called in virtual offset order for successive groups.
    void add_group(const vector<Message>& msgs, int64_t virtual_start, int64_t virtual_past_end);
    
    // Unhide overloads from the base
    using StreamIndexBase::find;
    using StreamIndexBase::add_group;
    
protected:
    
    /// Call the given iteratee for each node ID relevant to the given message.
    /// If the iteratee returns false, stop iteration.
    /// Calls the iteratee with 0 if and only if there are no node IDs relevant to the message.
    /// Must be specialized and implemented for any message type for which the StreamIndex is to be instantiated.
    void for_each_id(const Message& msg, const function<bool(const id_t&)> iteratee) const;
    
};

/// Define a GAM index as a stream index over a stream of Alignments
using GAMIndex = StreamIndex<Alignment>;


////////////
// Template Implementations
////////////

template<typename Message>
auto StreamIndex<Message>::find(cursor_t& cursor, id_t min_node, id_t max_node,
    const function<void(const Message&)> handle_result) const -> void {
    
    find(cursor, vector<pair<id_t, id_t>>{{min_node, max_node}}, handle_result);
    
}

template<typename Message>
auto StreamIndex<Message>::find(cursor_t& cursor, id_t node_id, const function<void(const Message&)> handle_result) const -> void {
    find(cursor, node_id, node_id, std::move(handle_result));
}

template<typename Message>
auto StreamIndex<Message>::find(cursor_t& cursor, const vector<pair<id_t, id_t>>& ranges,
    const function<void(const Message&)> handle_result, bool only_fully_contained) const -> void {
    
#ifdef debug
    cerr << "Begin a find query on ranges:" << endl;
    for (auto& range : ranges) {
        cerr << "\t" << range.first << "-" << range.second << endl;
    }
#endif
    
    // We need seek support
    assert(cursor.tell_raw() != -1);
    
    // Because a node in a later range may appear earlier in the file than a
    // node in an earlier range (but in a high-in-the-hierarchy bin), in
    // general we need to jump around in the file. TODO: Use a processed_up_to
    // counter to constrain us to one sweep in the only_fully_contained case.
    
    // To prevent us from scanning groups multiple times over, we keep a map
    // from already-processed group start VO to the VO of the next group (or
    // EOF). We can ride down chains in this map whenever we hit somewhere we
    // have already been, instead of actually re-reading anything.
    unordered_map<int64_t, int64_t> next_unprocessed;
    
    // We access it with this accessor function. It returns the given address
    // if the group there has not been read, or the next unprocessed VO (or EOF
    // VO) if it has.
    auto get_next_unprocessed = [&](int64_t currently_at) {
        // If we have to chain through multiple VOs to find the final one, we store them here.
        vector<int64_t> chain;
        
#ifdef debug
        cerr << "Find next unprocessed group after " << currently_at << endl;
#endif
        
        auto found = next_unprocessed.find(currently_at);
        while(found != next_unprocessed.end()) {
            // We have a place to go.
            
            // Remember this place as a place that needs to go to the final place we find.
            chain.push_back(currently_at);
            
#ifdef debug
            cerr << currently_at << " chains to " << found->second << endl;
#endif
            
            // Advance to the place we found.
            currently_at = found->second;
            found = next_unprocessed.find(currently_at);
        }
        
        // Now we hit the end. Save the final answer back to the map for
        // everything but the last item, so we never need to scan it again
        for (size_t i = 0; i + 1 < chain.size(); i++) {
            next_unprocessed[chain[i]] = currently_at;
        }
        
#ifdef debug
        cerr << "It is " << currently_at << endl;
#endif
        
        return currently_at;
    };
    
    // And this accessor marks a group as processed
    auto mark_processed = [&](int64_t start_vo, int64_t past_end_vo) {
    
#ifdef debug
        cerr << "Mark group " << start_vo << " to " << past_end_vo << " as processed" << endl;
#endif
    
        next_unprocessed[start_vo] = past_end_vo;
    };
    
    for (auto& range : ranges) {
        // For each range of IDs to look up
        
#ifdef debug
        cerr << "Look up range " << range.first << "-" << range.second << endl;
#endif
        
        find(range.first, range.second, [&](int64_t start_vo, int64_t past_end_vo) -> bool {
            // For each matching range of virtual offsets in the index
            
#ifdef debug
            cerr << "Look at VOs " << start_vo << "-" << past_end_vo << endl;
#endif
            
            // Warp the start past any already-processed groups we know about
            start_vo = get_next_unprocessed(start_vo);
            if (start_vo >= past_end_vo) {
                // Skip this whole range and look at the next one
                
#ifdef debug
                cerr << "The VO range has already been processed." << endl;
#endif
                
                return true;
            }
            
            // Now the range starts with a group we have never seen before.
            
            // Seek the cursor, even if we are already at the group in question.
            // TODO: We don't have a good way to tell if we are at the beginning of a group or not.
           
#ifdef debug
            cerr << "Seek cursor to " << start_vo << endl;
#endif
                
            cursor.seek_group(start_vo);
            
            // We need to track each group we encounter, so we can tell when an
            // entire group is past the top end of the ID range we are
            // currently looking up.
            int64_t group_vo = cursor.tell_group();
            id_t group_min_id = numeric_limits<id_t>::max();
            while (cursor.has_next() && cursor.tell_group() < past_end_vo) {
                // Read each message until we find a group that starts out of range
                
                // Which group is this message in?
                auto message_group_vo = cursor.tell_group();
                
                if (message_group_vo != group_vo) {
                    // We finished the previous group.
                    
#ifdef debug
                    cerr << "Finished group " << group_vo << endl;
#endif

                    // Record the group as processed
                    mark_processed(group_vo, message_group_vo);
                    
                    if (group_min_id != numeric_limits<id_t>::max() && group_min_id > range.second) {
                        // Everything in the (non-empty) previous group was too high. We don't care about this group; our iteration is over.
                        
#ifdef debug
                        cerr << "Group was out of bounds for its range with min id " << group_min_id << " > " << range.second << endl;
                        cerr << "Move on to next range" << endl;
#endif
                        
                        // Stop early. Don't finish this run and don't look at the next runs for this query range.
                        return false;
                    }
                    
                    // Otherwise we need to start a new group
                    group_min_id = numeric_limits<id_t>::max();
                    
                    // Zip the group VO ahead to the next unprocessed group (which may be here, or at EOF)
                    group_vo = get_next_unprocessed(message_group_vo);
                    if (group_vo != message_group_vo) {
                        // We want to go to a different group next.
                        if (group_vo >= past_end_vo) {
                            // But it's out of range for this range. Don't go there.
#ifdef debug
                            cerr << "Next unprocessed VO is out of range" << endl;
#endif
                            break;
                        } else {
                            // Seek there and restart the loop to see if we found anything good.
#ifdef debug
                            cerr << "Seek to next unprocessed VO at " << group_vo << endl;
#endif
                            cursor.seek_group(group_vo);
                            continue;
                        }
                    } else {
                        // Otherwise, we are continuing with this group we just found.
#ifdef debug
                        cerr << "Next unprocessed VO is right here." << endl;
#endif
                    }
                }
                
                // Filter the message by the query and yield it if it matches
                const auto& message = *cursor;
                bool message_match = false;
                
                if (message.path().mapping_size() == 0) {
                    // This read is unmapped, so count it as node 0.
                    group_min_id = min(group_min_id, (id_t)0);
                    if (is_in_range(ranges, 0)) {
                        // We want unmapped reads.
                        message_match = true;
                    }
                } else {
                    // The read has mappings
                    for (const auto& mapping : message.path().mapping()) {
                        // Look at each node that is visited
                        auto visited = mapping.position().node_id();
                        group_min_id = min(group_min_id, visited);
                        if (is_in_range(ranges, visited)) {
                            // We want this node.
                            message_match = true;
                            if (!only_fully_contained) {
                                // All we care about is that any of the nodes match
                                break;
                            }
                        } else if (only_fully_contained) {
                            // We need *all* of the nodes to match, and this one didn't.
                            message_match = false;
                            break;
                        }
                    }
                }
                
                if (message_match) {
                    // This message is one that matches the query. Yield it.
                    handle_result(message);
                }
                
                // Look for the next message
                cursor.get_next();
                
            }
           
            if (group_vo < past_end_vo) {
                // We finished a final group, from group_vo to past_end_vo
#ifdef debug
                cerr << "Finished last group " << group_vo << endl;
#endif

                // Mark it finished
                mark_processed(group_vo, past_end_vo);

                if (group_min_id != numeric_limits<id_t>::max() && group_min_id > range.second) {
                    // If the (non-empty) last group had all its node IDs past the max
                    // node ID, we know nothing after it can possibly match, so stop
                    // iteration.
                    
#ifdef debug
                    cerr << "Group was out of bounds with min id " << group_min_id << " > " << range.second << endl;
                    cerr << "Move on to next range" << endl;
#endif
                
                    return false;
                }
            }
            
            // Otherwise, the last group we looked at was not yet out of bounds, so get another range to look at, if one exists.
            return true;
        });
        
    }
}

template<typename Message>
auto StreamIndex<Message>::index(cursor_t& cursor) -> void {
    // Keep track of what group we are in 
    int64_t group_vo = cursor.tell_group();
    // And load all its messages
    vector<Message> group;
    
    // We need to have seek support
    assert(group_vo != -1);
    
    while (cursor.has_next()) {
        // For each message
        
        // Work out what group it is in
        int64_t message_group_vo = cursor.tell_group();
        
        if (message_group_vo != group_vo) {
            // This is the start of a new group
            
            // Record the old group as being up to here
            add_group(group, group_vo, message_group_vo);
            
            // Set up for the new group
            group.clear();
            group_vo = message_group_vo;
        }
        
        // Add the message to the group and move on
        group.emplace_back(std::move(cursor.take()));
    }
    
    if (!group.empty()) {
        // Record the final group. Use wherever the cursor landed at the end as its final virtual offset.
        add_group(group, group_vo, cursor.tell_raw());
    }
}

template<typename Message>
auto StreamIndex<Message>::add_group(const vector<Message>& msgs, int64_t virtual_start, int64_t virtual_past_end) -> void {
    // Find the min and max ID visited by any of the messages
    id_t min_id = numeric_limits<id_t>::max();
    id_t max_id = numeric_limits<id_t>::min();
    
    for (auto& msg : msgs) {
        // For each message
        if (msg.path().mapping_size() == 0) {
            // The read is unmapped, so it belongs to node ID 0
            min_id = min(min_id, (id_t)0);
            max_id = max(max_id, (id_t)0);
        } else {
            for (auto& mapping : msg.path().mapping()) {
                // For each mapping in it, min/max in the ID
                auto id = mapping.position().node_id();
                min_id = min(min_id, id);
                max_id = max(max_id, id);
            }
        }
    }
    
    add_group(min_id, max_id, virtual_start, virtual_past_end);
}

}

#endif


