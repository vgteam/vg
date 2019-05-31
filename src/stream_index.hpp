#ifndef VG_STREAM_INDEX_HPP_INCLUDED
#define VG_STREAM_INDEX_HPP_INCLUDED

/**
 * \file stream_index.hpp
 * Contains the StreamIndex template, which allows lookup by relevant node ID in sorted VPKG-formatted files.
 */
 
#include <iostream>
#include <vector>
#include <set>
#include <unordered_map>
#include <type_traits>

#include "types.hpp"
#include <vg/vg.pb.h>
#include <vg/io/protobuf_iterator.hpp>
#include "scanner.hpp"

namespace vg {

using namespace std;

// For efficiently storing the bins that are populated in the relatively sparse
// bin space that the index uses, we use a compact prefix trie (radix tree) on
// the node ID bit string prefixes that correspond to the bins.

// To make that work we in turn need a bit string type that is easily
// convertible to/from 64-bit numbers.

/**
 * Represents a string of up to 64 bits.
 */
class BitString {
public:
    
    /// Make a new BitString representing the low length bits of the given number
    BitString(uint64_t bits, size_t length);
    
    /// Make an empty BitString
    BitString();
    
    // Copyable and movable
    BitString& operator=(const BitString& other) = default;
    BitString& operator=(BitString&& other) = default;
    BitString(const BitString& other) = default;
    BitString(BitString&& other) = default;
    
    /// Convert the BitString back to a number
    uint64_t to_number() const;
    
    /// Get a suffix of a BitString by dropping the specified number of bits.
    /// If we drop all the bits or more, get an empty BitString.
    BitString drop_prefix(size_t prefix_length) const;
    
    /// Split into a prefix of the given length and a suffix of the rest
    pair<BitString, BitString> split(size_t prefix_length) const;
    
    /// Determine if two BitStrings are equal
    bool operator==(const BitString& other) const;
    
    /// Determine if two BitStrings are unequal
    bool operator!=(const BitString& other) const;
    
    /// Get the length of the longest common prefix (index of the first
    /// mismatching bit) between this BitString and another.
    size_t common_prefix_length(const BitString& other) const;
    
    /// Return true if one BitString is a prefix of the other, or if this BitString has the 0 at the first differing bit.
    bool at_or_before(const BitString& other) const;
    
    /// Return true if one BitString is a prefix of the other, or if the other BitString has the 0 at the first differing bit.
    bool at_or_after(const BitString& other) const;
    
    /// Peek at the top bit and see if it is a 1 (true) or 0 (false).
    /// Empty bit strings get false.
    bool peek() const;
    
    /// Get the length of the BitString
    size_t length() const;
    
    /// Return true if the bit string is empty
    bool empty() const;
    
protected:
    /// Holds the actual bits of the bit string.
    /// We store the bits aligned to the left to make finding the first mismatch easier.
    /// All bits below the last used bit are 0.
    uint64_t bits;
    
    /// Holds the number of bits that are used.
    uint8_t bit_length;
    
    /// How many total bits are possible?
    const static size_t TOTAL_BITS = numeric_limits<uint64_t>::digits;
};

/// Allow BitStrings to be printed for debugging
ostream& operator<<(ostream& out, const BitString& bs);

/**
 * Represents a radix tree keyed by/internally using BitStrings.
 * Each item has a BitString as a key, and items are stored in a trie/prefix tree.
 * Each node has at most one item and at most two children.
 * Movable but not copyable.
 * TODO: Implement copy.
 */
template<typename Item>
class BitStringTree {
public:

    /// Insert the given item under the given key
    void insert(const BitString& key, const Item& value);
    
    /// Enumerate items whose keys match prefixes of the given key, in order from most specific to least specific.
    /// Returns false if stopped early.
    bool traverse_up(const BitString& key, const function<bool(const Item&)>& iteratee) const;
    
    /// Enumerate items which have the given low and high bit strings as a
    /// prefix, or which fall between them, as an in-order traversal of the
    /// tree. Returns false if stopped early.
    bool traverse_in_order(const BitString& low, const BitString& high, const function<bool(const Item&)>& iteratee) const;

protected:
    struct TreeNode {
        /// Each TreeNode represents a prefix over its parent
        BitString prefix;
        /// Each TreeNode holds a 0 child and a 1 child, at most
        /// Their prefixes say what full prefixes they actually correspond to
        unique_ptr<TreeNode> children[2];
        /// Each TreeNode also may hold an item record
        Item content;
        /// This is set if we actually have one
        bool has_content = false;
        
        /// Insert the given item at or under this node.
        /// Its key must have had our prefix already removed from it.
        void insert(const BitString& key, const Item& value);
        
        /// Search down to the node that corresponds to the given key, or where
        /// it would be. Call the iteratee for that node and every parent, from
        /// bottom to top, until the iteratee returns false. If not found, do
        /// not call the iteratee. Retruns true if the iteratee did not ask to
        /// stop, and false otherwise. The key will have already had this
        /// node's prefix removed.
        bool traverse_up(const BitString& key, const function<bool(const Item&)>& iteratee) const;
        
        /// Iterate over elements in the tree with an in-order traversal
        /// between the two given keys, inclusive. Low and high have already
        /// had this node's prefix removed.
        bool traverse_in_order(const BitString& low, const BitString& high, const function<bool(const Item&)>& iteratee) const;
        
        // Note that depth is bounded so we don't need recursion-breaking in our destructor
        
    };
    
    /// The root node has an empty prefix
    TreeNode root;
    
};

/**
 * An index for a node-ID-sorted VPKG-formatted Protobuf file, such as GAM or VG.
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
 * then into 7, 8, 9, 10, 11, 12, 13, and 14, and so on.
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
    /// linear index, or the past-the-end of the previous run scanned (which
    /// should be moot, because runs should not overlap in the index).
    void find(id_t node_id, const function<bool(int64_t, int64_t)> scan_callback) const;
    
    /// Find all the ranges of run virtual offsets to check for reads visiting
    /// the given inclusive node ID range. Relies on a scanning callback, which
    /// will be called repeatedly with the start and past-the-end virtual
    /// offsets of runs which may contain groups touching the given node ID.
    /// When called, the callback should scan the run and return either true if
    /// it wants the next run, or false if it encountered a group with an
    /// out-of-range start and wants to stop iteration. Runs will be emitted in
    /// order, and truncated on the left to either the appropriate lower bound
    /// from the linear index, or the past-the-end of the previous run scanned.
    void find(id_t min_node, id_t max_node, const function<bool(int64_t, int64_t)> scan_callback) const;
    
    /// Iterate over ranges of virtual offsets from the end of the file to the
    /// start. The ranges to *not* necessarily correspond to runs. The ending
    /// VO of the first range iterated may be numeric_limits<int64_t>::max().
    /// The start VO of each range is guaranteed to be a valid VO of a group.
    /// The past-end VO may be the start of the previously iterated range.
    /// We have no scan_forward because you can just do that with a cursor.
    /// Stops when the callback returns false.
    void scan_backward(const function<bool(int64_t, int64_t)> scan_callback) const;
    
    /// Add a group into the index, based on its minimum and maximum
    /// (inclusive) used node IDs. Must be called for all groups in virtual
    /// offset order.
    void add_group(id_t min_id, id_t max_id, int64_t virtual_start, int64_t virtual_past_end);
    
    ///////////////////
    // Lowest-level functions for thinking about bins and windows.
    ///////////////////
    
    /// Get the ID prefix bits corresponding to a bin
    static BitString bin_to_prefix(bin_t bin);
    
    /// Get the given ID as a bit string
    static BitString id_to_prefix(id_t id);
    
    /// Get the most specific bin that contains both of the given node IDs.
    static bin_t common_bin(id_t a, id_t b);
    
    /// Get the linear index window that the given node ID falls in. The window
    /// range for a group is its min nodes' window through its max node's
    /// window.
    static window_t window_of_id(id_t id);
    
    /// Iterate over the *populated* bins in the index, in in-order bin tree
    /// traversal order, that any of the node IDs in the given inclusive range
    /// occur in. Returns false if asked to stop.
    bool used_bins_of_range(id_t min_id, id_t max_id, const function<bool(bin_t)>& iteratee) const;
    
    
protected:
    // How many bits of a node ID do we truncate to get its linear index window?
    const static size_t WINDOW_SHIFT = 8;
    
    /// Maps from bin number to all the ranges of virtual offsets, in order, for runs that land in the given bin.
    /// A run lands in a bin if that bin is the most specific bin that includes both its lowest and highest nodes it uses.
    unordered_map<bin_t, vector<pair<int64_t, int64_t>>> bin_to_ranges;
    
    /// Maps from the bit string representing the prefix of node IDs that a bin
    /// matches to the bin's bin number. Only contains entries for nonempty
    /// bins.
    BitStringTree<bin_t> bins_by_id_prefix;
    
    /// Maps from linear index window to the virtual offset of the first group
    /// that overlaps that window (taking the group as a min-to-max node
    /// range). If you are looking for reads that visit a node, they can't
    /// possibly occur in a group before the first offset stored for the node's
    /// window (or any greater window). TODO: Should we make this a vector
    /// instead and hope nobody uses high/sparse node IDs?
    map<window_t, int64_t> window_to_start;
    
    /// What was the minimum node ID of the last group added?
    /// If this isn't strictly increasing, we're trying to index data that is not sorted.
    id_t last_group_min_id = numeric_limits<id_t>::min();
    
    /// Return true if the given ID is in any of the sorted, coalesced, inclusive ranges in the vector, and false otherwise.
    /// TODO: Is repeated binary search on the ranges going to be better than an unordered_set of all the individual IDs?
    static bool is_in_range(const vector<pair<id_t, id_t>>& ranges, id_t id);
    
private:
    // Not copyable because we contain pointers.
    StreamIndexBase(const StreamIndexBase& other) = delete;
    StreamIndexBase& operator=(const StreamIndexBase& other) = delete;

};

/** An index that provides a higher-level API in terms of the actual messages
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
    using cursor_t = vg::io::ProtobufIterator<Message>;
    
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
    /// IDs may repeat.
    /// If the iteratee returns false, stop iteration.
    /// Calls the iteratee with 0 only if there are no relevant node IDs *or* the message is relevant to queries for unplaced content.
    void for_each_id(const Message& msg, const function<bool(const id_t&)> iteratee) const;
    
};

/// Define a GAM index as a stream index over a stream of Alignments
using GAMIndex = StreamIndex<Alignment>;


////////////
// Template Implementations
////////////

template<typename Item>
void BitStringTree<Item>::insert(const BitString& key, const Item& value) {
    root.insert(key, value);
}

template<typename Item>
bool BitStringTree<Item>::traverse_up(const BitString& key, const function<bool(const Item&)>& iteratee) const {
    return root.traverse_up(key, iteratee);
}

template<typename Item>
bool BitStringTree<Item>::traverse_in_order(const BitString& low, const BitString& high, const function<bool(const Item&)>& iteratee) const {
    return root.traverse_in_order(low, high, iteratee);
}

template<typename Item>
void BitStringTree<Item>::TreeNode::insert(const BitString& key, const Item& value) {

#ifdef debug
    cerr << "Inserting key " << key << " under " << this << " with value " << value << endl;
#endif
    
    if (key.empty()) {
        // It goes here.
        // We can only take one item
#ifdef debug
        cerr << "Item belongs here" << endl;
#endif
        assert(!has_content);
        content = value;
        has_content = true;
    } else {
        // Get the first bit of its prefix
        bool lead_bit = key.peek();
        
        if (!children[lead_bit]) {
            // We need to make a new child to hold this item
            children[lead_bit] = unique_ptr<TreeNode>(new TreeNode());
            // Populate it
            children[lead_bit]->prefix = key;
            children[lead_bit]->content = value;
            children[lead_bit]->has_content = true;
#ifdef debug
            cerr << "Item belongs in new child " << children[lead_bit].get() << " with lead bit "
                << lead_bit << " and prefix " << children[lead_bit]->prefix << endl;
#endif
        } else {
            // We already have a child in this direction.
            
#ifdef debug
            cerr << "Item belongs on branch with existing child " << children[lead_bit].get() << " with lead bit "
                << lead_bit << " and prefix " << children[lead_bit]->prefix << endl;
#endif
            
            // See where the key diverges from our child's key
            auto breakpoint = children[lead_bit]->prefix.common_prefix_length(key);
            
            if (breakpoint >= children[lead_bit]->prefix.length()) {
                // This key belongs inside the child node we have
                
#ifdef debug
                cerr << "Key " << key << " is a prefix of " << children[lead_bit]->prefix << " so item lives at or under child" << endl;
#endif
                
                // Insert recursively
                children[lead_bit]->insert(key.drop_prefix(breakpoint), value);
            } else {
                // The item to be added diverges somewhere along the branch to the child.
                // And it isn't at the 0th bit because we organized our children by bit 0.
                
#ifdef debug
                cerr << "Key " << key << " matches " << children[lead_bit]->prefix << " up through " << breakpoint << endl;
#endif
                
                // Break up the child's key
                auto prefix_parts = children[lead_bit]->prefix.split(breakpoint);
                
                // Create a new node to sit at the split point and wire it in
                unique_ptr<TreeNode> new_child(new TreeNode());
                new_child->prefix = prefix_parts.first;
                children[lead_bit]->prefix = prefix_parts.second;
                new_child->children[prefix_parts.second.peek()] = move(children[lead_bit]);
                children[lead_bit] = move(new_child);
                
#ifdef debug
                cerr << "Added new node " << children[lead_bit].get() << " with shared prefix " << children[lead_bit]->prefix
                    << " and make its " << prefix_parts.second.peek() << " child the old child with prefix " << prefix_parts.second << endl;
#endif
                
                // Recursively insert either at the breakpoint node or under it in the other child slot
                children[lead_bit]->insert(key.drop_prefix(breakpoint), value);
            }
        }
    }
}

template<typename Item>
bool BitStringTree<Item>::TreeNode::traverse_up(const BitString& key, const function<bool(const Item&)>& iteratee) const {
    if (key.empty()) {
        // We are the item being sought
        if (has_content) {
            // We actually have an item, so send it
            return iteratee(content);
        } else {
            // No item so it can't ask to stop
            return true;
        }
    } else {
        // The item must belong to a child slot
        // Get the first bit of its prefix
        bool lead_bit = key.peek();
        
        if (children[lead_bit]) {
            // We have a child that would be responsible for the key
            
            // But how long is the match
            auto breakpoint = children[lead_bit]->prefix.common_prefix_length(key);
            
            if (breakpoint >= children[lead_bit]->prefix.length()) {
                // This key belongs inside the child.
                // Search recursively
                if (children[lead_bit]->traverse_up(key.drop_prefix(breakpoint), iteratee)) {
                    // The child returned true
                    if (has_content) {
                        // We also have an item, so send it
                        return iteratee(content);
                    } else {
                        // We have no item, so just pass up the true
                        return true;
                    }
                } else {
                    // The iteratee has already stopped.
                    return false;
                }
            } else {
                // The key branches off before the child.
                // So we have to process ourselves as the bottom node
                if (has_content) {
                    // We actually have an item, so send it
                    return iteratee(content);
                } else {
                    // No item so it can't ask to stop
                    return true;
                }
            }
        } else {
            // No child exists that is responsible for the key.
            // We have no results under us, so just do us.
            if (has_content) {
                // We have an item, so send it
                return iteratee(content);
            } else {
                // We have no item, so just pass up the true
                return true;
            }
        }
    }
}

template<typename Item>
bool BitStringTree<Item>::TreeNode::traverse_in_order(const BitString& low, const BitString& high,
    const function<bool(const Item&)>& iteratee) const {
    
    // This is actually pretty easy. We know we're in the range. We just need
    // to figure out if our children are in the range, and if so, call them.
    
    // We use special bit string comparison, so our bounded range is all the
    // stuff that we don't have bits showing it comes before the low or after
    // the high.
    
#ifdef debug
    cerr << "Arrived at node " << this << " with range " << low << " - " << high << endl;
#endif
    
    /// Define a function to process each child.
    /// Returns false if we stop early
    auto do_child = [&](bool child_index) -> bool {
        // Grab a reference to the child's unique_ptr
        auto& child = children[child_index];
        
        if (child) {
            // The child exists. Grab its prefix
            auto& child_prefix = child->prefix;
            
#ifdef debug
            cerr << "Child " << child.get() << " has prefix " << child_prefix << endl; 
#endif
            
            if (low.at_or_before(child_prefix) && high.at_or_after(child_prefix)) {
                // The child is in the range.
                
                // But we need to work out which of the range bounds are
                // outside the region which the child is responsible for and
                // not pass them on.
                BitString child_low = child_prefix.at_or_before(low) ? low.drop_prefix(child_prefix.length()) : BitString();
                BitString child_high = child_prefix.at_or_after(high) ? high.drop_prefix(child_prefix.length()) : BitString();
                
                return child->traverse_in_order(child_low, child_high, iteratee);
            }
        }
        
        // If we get here we didn't need to visit the child at all.
        return true;
    };
   
    // Do the left child, then us if we have an item, then the right child.
    return do_child(false) && ((has_content && iteratee(content)) || !has_content) && do_child(true);
}

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
    assert(cursor.tell_group() != -1);
    
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
            while (cursor.has_current() && cursor.tell_group() < past_end_vo) {
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
                
                for_each_id(message, [&](const id_t& found) {
                    // For each ID touched by the message
                    
                    // Min it in, keeping 0 as the sentinel for no nodes touched.
                    group_min_id = min(group_min_id, found);
                    if (is_in_range(ranges, found)) {
                        // We want this node (or unplaced messages like this one).
                        message_match = true;
                        if (!only_fully_contained) {
                            // All we care about is that any of the nodes match.
                            // we know enough to keep this message.
                            return false;
                        }
                    } else if (only_fully_contained) {
                        // We need *all* of the nodes to match, and this one didn't.
                        message_match = false;
                        // We know enough to discard this message.
                        return false;
                    }
                    // Keep looking
                    return true;
                });
                
                if (message_match) {
                    // This message is one that matches the query. Yield it.
                    handle_result(message);
                }
                
                // Look for the next message
                cursor.advance();
                
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
    
    while (cursor.has_current()) {
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
        add_group(group, group_vo, cursor.tell_group());
    }
}

template<typename Message>
auto StreamIndex<Message>::add_group(const vector<Message>& msgs, int64_t virtual_start, int64_t virtual_past_end) -> void {
    // Find the min and max ID visited by any of the messages
    id_t min_id = numeric_limits<id_t>::max();
    id_t max_id = numeric_limits<id_t>::min();
    
    for (auto& msg : msgs) {
        // For each message
        for_each_id(msg, [&](const id_t& found) {
            // For each ID touched by the message
            
            // Min and max in the ID, keeping 0 to represent no mappings.
            min_id = min(min_id, found);
            max_id = max(max_id, found);
            
            // Don't stop early
            return true;
        });
    }
    
    add_group(min_id, max_id, virtual_start, virtual_past_end);
}

template<typename Message>
auto StreamIndex<Message>::for_each_id(const Message& msg, const function<bool(const id_t&)> iteratee) const -> void {
    // Visit all the IDs.
    // Zeros will come out if the message is empty (unplaced) or has an unplaced child message.
    // Duplicates will come out but that is fine.
    IDScanner<Message>::scan(msg, iteratee);
}

}

#endif


