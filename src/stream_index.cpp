#include "stream_index.hpp"

#include <iostream>
#include <queue>

namespace vg {

using namespace std;

BitString::BitString(uint64_t bits, size_t length) : bits(bits), bit_length(length) {
    // Make sure we didn't overflow
    assert(length <= numeric_limits<uint8_t>::max());
    
    if (bit_length == 0) {
        // Just remove all the bits
        this->bits = 0;
    } else {
        // We won't shift out all the bits so it is safe to shift
        this->bits = this->bits << (TOTAL_BITS - bit_length);
    }
}

BitString::BitString() : BitString(0, 0) {
    // Nothing to do
}

auto BitString::to_number() const -> uint64_t {
    auto total_bits = CHAR_BIT * sizeof(decltype(this->bits));
    
    if (bit_length == 0) {
        // We have no used bits
        return 0;
    } else {
        // We won't shift out all the bits so it is safe to shift
        return bits >> (TOTAL_BITS - bit_length);
    }
}

auto BitString::drop_prefix(size_t prefix_length) const -> BitString {
    // Don't let it be too long
    assert(prefix_length <= bit_length);
    
    if (prefix_length == bit_length) {
        // We are losing all our bits
        return BitString();
    }
    
    // Otherwise shift off the dropped bits and return that
    BitString other;
    other.bits = bits << prefix_length;
    other.bit_length = bit_length - prefix_length;
    
    return other;
}

auto BitString::split(size_t prefix_length) const -> pair<BitString, BitString> {
    // Make the prefix
    auto non_prefix_bits = TOTAL_BITS - prefix_length;
    BitString prefix;
    if (non_prefix_bits < TOTAL_BITS) {
        // The prefix will be nonempty
        prefix.bits = (bits >> non_prefix_bits) << non_prefix_bits;
        prefix.bit_length = prefix_length;
    }
    
    // Return the prefix and the rest
    return make_pair(prefix, drop_prefix(prefix_length));
}

auto BitString::operator==(const BitString& other) const -> bool {
    // There's nothing fancy here because we demand all unused bits in the storage are 0.
    return (bits == other.bits && bit_length == other.bit_length);
}

auto BitString::operator!=(const BitString& other) const -> bool {
    // There's nothing fancy here because we demand all unused bits in the storage are 0.
    return (bits != other.bits || bit_length != other.bit_length);
}

auto BitString::common_prefix_length(const BitString& other) const -> size_t {
   // Make a mask where identical bits are 0
   auto mask = bits ^ other.bits;
   // Count the identical bits (count unset leading bits) with a compiler builtin
   size_t identical_bits = __builtin_clzll(mask);
   // Return the number of matching bits before the first mismatch, or the length of the shorter BitString
   return min(min(identical_bits, (size_t) bit_length), (size_t) other.bit_length);
}

auto BitString::peek() const -> bool {
    if (bit_length == 0) {
        return false;
    }
    
    // If we aren't empty, get the high bit
    return ((uint64_t)1 << 63) & bits;
}

auto BitString::length() const -> size_t {
    return bit_length;
}

auto BitString::empty() const -> bool {
    return bit_length == 0;
}

const string StreamIndexBase::MAGIC_BYTES = "GAI!";

auto StreamIndexBase::bin_to_prefix(bin_t bin) -> BitString {

    // The bin is an offset with the low n bits set, plus an n-bit index.
    // We have to figure out the value of the all-1s offset.
    bin_t bin_offset;
    // And how many bits are used (n)
    size_t used_bits;
    
    // The offset will have as many 1s as the used bits in the number if the
    // number is all 1s after some leading 0s. Otherwise the offset will have 1
    // fewer bits than the number.
    
    // Count the leading 0 bits
    size_t leading_zeros = __builtin_clzll(bin);
    
    if (leading_zeros == numeric_limits<bin_t>::digits) {
        // All zero bin is the root and gets the empty prefix
        // Don't even bother with the rest of the logic
        return BitString();
    } else if (leading_zeros == numeric_limits<bin_t>::digits - 1) {
        // It must be bin 1, which has a 1-bit offset
        bin_offset = 1;
        used_bits = 1;
    } else {
        // Make an all-1s value as wide as our bin number
        bin_offset = ((bin_t)~0) >> leading_zeros;
        
        // Estimate the used bits; the bits that aren't 0s are 1s in the offset.
        used_bits = numeric_limits<bin_t>::digits - leading_zeros;
        
        if (bin_offset > bin) {
            // We need 1 fewer bits for the offset.
            bin_offset = bin_offset >> 1;
            // Correct the used bits in the bin index
            used_bits--;
        }
    }
    
    /// Subtract out the offset to get the bin index, and use it as a prefix
    /// with the appropriate number of bits based on what the offset was.
    return BitString(bin - bin_offset, used_bits);
}

auto StreamIndexBase::bins_of_id(id_t id) -> vector<bin_t> {
    vector<bin_t> to_return;
    
    // How many bits are in the number? We get a bin per bit.
    auto number_bits = CHAR_BIT * sizeof(bin_t);
    
    // We need to keep in mind that shifting *all* the bits out of a number in
    // one go is undefined behavior.
    
    // Bin number consists of an index which is a prefix of the number
    bin_t bin_index = ((bin_t)id) >> 1;
    // And an offset to keep from colliding with other bins.
    bin_t bin_offset = ((bin_t)~0) >> 1;
    
    for (int i = 0; i < number_bits; i++) {
        
#ifdef debug
        cerr << hex << id << " " << i << " " << bin_index << " " << bin_offset << " " << (bin_index + bin_offset) << endl;
#endif
    
        to_return.push_back(bin_index + bin_offset);
        bin_index = bin_index >> 1;
        bin_offset = bin_offset >> 1;
    }
    
    return to_return;
}

auto StreamIndexBase::bins_of_range(id_t min_id, id_t max_id) -> vector<bin_t> {
    // We can just get the two bin vectors for the ending bin, and generate an inclusive range of bins at each level.
    
    vector<bin_t> to_return;
    
    auto min_bins = bins_of_id(min_id);
    auto max_bins = bins_of_id(max_id);
    
    assert(min_bins.size() == max_bins.size());
    
    for (size_t i = 0; i < min_bins.size(); i++) {
        // For each specificity level
        for (bin_t bin = min_bins[i]; bin != max_bins[i] + 1; bin++) {
            // For each bin in the inclusive range, emit it
            to_return.push_back(bin);
        }
    }
    
    return to_return;
    
}

auto StreamIndexBase::bins_of_range(id_t min_id, id_t max_id, const function<bool(bin_t)>& iteratee) -> bool {

    // Set this to false if the iteratee asks to stop
    bool keep_going = true;

    // How many bits are in the number? We get a bin per bit.
    auto number_bits = CHAR_BIT * sizeof(bin_t);
    
    // We need to keep in mind that shifting *all* the bits out of a number in
    // one go is undefined behavior.
    
    // Bin number consists of an index which is a prefix of the number
    bin_t min_bin_index = ((bin_t)min_id) >> 1;
    bin_t max_bin_index = ((bin_t)max_id) >> 1;
    // And an offset to keep from colliding with other bins.
    bin_t bin_offset = ((bin_t)~0) >> 1;
    
    for (int i = 0; keep_going && i < number_bits; i++) {
        // For each level
        
        for (bin_t bin_index = min_bin_index; keep_going && bin_index <= max_bin_index; bin_index++) {
            // For each bin in the range of bins at this level
            
            // Compute its actual offset bin number and feed it to the iteratee
            keep_going &= iteratee(bin_index + bin_offset);
        }
        
        // Ascend a level
        min_bin_index = min_bin_index >> 1;
        max_bin_index = max_bin_index >> 1;
        bin_offset = bin_offset >> 1;
    }
    
    return keep_going;
}
    
auto StreamIndexBase::common_bin(id_t a, id_t b) -> bin_t {
    // Convert to unsigned numbers
    bin_t a_bin = a;
    bin_t b_bin = b;
    
    // Define the offset for the bin
    bin_t offset = ((bin_t)~0);
    
    // We're just going to pop off bits until we find the common prefix.
    // Always pop off one bit, even if we are binning a number and itself.
    // TODO: Find a faster way to do this with the appropriate instruction intrinsics.
    do {
        a_bin = a_bin >> 1;
        b_bin = b_bin >> 1;
        offset = offset >> 1;
    } while(a_bin != b_bin);
    return a_bin + offset;
}

auto StreamIndexBase::window_of_id(id_t id) -> window_t  {
    return id >> WINDOW_SHIFT;
}

auto StreamIndexBase::add_group(id_t min_id, id_t max_id, int64_t virtual_start, int64_t virtual_past_end) -> void {
    
    if (min_id < last_group_min_id) {
        // Someone is trying to index an unsorted GAM.
        // This is probably user error, so complain appropriately:
        cerr << "error [vg::GAMIndex]: GAM data being indexed is not sorted. Sort with vg gamsort." << endl;
        exit(1);
    }
    last_group_min_id = min_id;
    
    // Find the bin for the run
    bin_t bin = common_bin(min_id, max_id);
    
#ifdef debug
    cerr << "Group spanning " << min_id << "-" << max_id << " at "
        << virtual_start << "-" << virtual_past_end << " lands in bin " << bin << endl;
#endif
    
    // Find the existing ranges in the bin.
    // We know the previous one, if present, must end at or before this one's start.
    auto& ranges = bin_to_ranges[bin];
    
    if (ranges.empty()) {
        // We just made a new vector of bin ranges. Remember it.
        bins_by_id_prefix.insert(bin_to_prefix(bin), &ranges);
    }
    
    if (!ranges.empty() && ranges.back().second == virtual_start) {
        // We fit right after the last range.
        ranges.back().second = virtual_past_end;
#ifdef debug
        cerr << "Extend existing range to " << ranges.back().first << "-" << ranges.back().second << endl;
#endif
    } else {
        // We need a new range
        bin_to_ranges[bin].emplace_back(virtual_start, virtual_past_end);
    }
    
    for (window_t w = window_of_id(min_id); w <= window_of_id(max_id); w++) {
        // For each window that this group overlaps
        
        if (!window_to_start.count(w)) {
            // If it is the first group we encounter in the window, it must also be the earliest-staring group in the window.
            
            // This is the earliest virtual offset to overlap that window
            window_to_start[w] = virtual_start;
            
#ifdef debug
            cerr << "Start window " << w << endl;
#endif
        }
    }
}

auto StreamIndexBase::find(id_t node_id) const -> vector<pair<int64_t, int64_t>> {
    vector<pair<int64_t, int64_t>> to_return;
    
    find(node_id, [&](int64_t run_start, int64_t run_past_end) -> bool {
        // For each run we find, remember it
        to_return.emplace_back(run_start, run_past_end);
        // Keep getting runs until we run out in the index. We can't actually scan the data.
        return true;
    });
    
    return to_return;
}

auto StreamIndexBase::find(id_t node_id, const function<bool(int64_t, int64_t)> scan_callback) const -> void {
    // Look for a single-node inclusive range
    find(node_id, node_id, std::move(scan_callback));
}

auto StreamIndexBase::find(id_t min_node, id_t max_node, const function<bool(int64_t, int64_t)> scan_callback) const -> void {
    
#ifdef debug
    cerr << "Query for node range " << min_node << "-" << max_node << endl;
#endif
    
    // Find the window that gives us a lower bound on the virtual offset we
    // need to be at to find things that touch this node ID.
    window_t min_window = window_of_id(min_node);
    window_t max_window = window_of_id(max_node);
    
#ifdef debug
    cerr << "Looking for first filled window of " << min_window << "-" << max_window << endl;
#endif
    
    // Find the minimum virtual offset we need to consider
    int64_t min_vo = 0;
    // It will be for the first occupied window at or after the min window but not greater than the max window.
    auto found = window_to_start.lower_bound(min_window);
    if (found != window_to_start.end() && found->first <= max_window) {
        // Some groups overlapped this window, and they started here.
        min_vo = found->second;
        
#ifdef debug
        cerr << "First occupied window is " << found->first << " at offset " << min_vo << endl;
#endif
    } else {
        // No groups overlapped any window within the range, so don't iterate anything.
        
#ifdef debug
        cerr << "No windows occupied; range is empty" << endl;
#endif
        
        return;
    }
    
    // This will hold bins that actually have vectors in the index, as
    // iterators to them in the index.
    vector<decltype(bin_to_ranges)::const_iterator> used_bins;
    
    // Loop over the bins we have to deal with
    bins_of_range(min_node, max_node, [&](bin_t bin_number) -> bool {
        // See if the bin has any entries
        auto found = bin_to_ranges.find(bin_number);
        if (found != bin_to_ranges.end()) {
            // If it does, keep it
            used_bins.push_back(found);
        }
        
        // TODO: Is there a way we can just process all the bins one at a time,
        // instead of merging across them?
        
        return true;
    });
    
    // Define a cursor type within a bin.
    // A cursor has an iterator to a pair of start and past-end VOs that occur in a bin.
    // But we also need to be able to look at the cursor and know if it is done
    // So we also keep the index in used_bins that it belongs to
    using bin_cursor_t = pair<vector<pair<int64_t, int64_t>>::const_iterator, size_t>;
    
    // As we iterate we will be interested in the smallest cursor (i.e. the cursor to the earliest-starting run not already iterated.)
    // Define a way to get that
    // Since the queue puts the "greatest" element at the top, our "less" is really a "greater" to turn it around.
    auto greater_for_cursors = [&](const bin_cursor_t& a, const bin_cursor_t& b) {
        return a.first->first > b.first->first;
    };
    
    // We want a priority queue of cursors, so we can easily find the next one to use
    priority_queue<bin_cursor_t, vector<bin_cursor_t>, decltype(greater_for_cursors)> cursor_queue(greater_for_cursors);
    
    // Set up a cursor in each bin
    // TODO: Could we do one cursor per specificity level instead? That would be faster.
    for (size_t i = 0; i < used_bins.size(); i++) {
        // Look at the start of the bin
        bin_cursor_t cursor = make_pair(used_bins[i]->second.begin(), i);
        
        while(cursor.first != used_bins[cursor.second]->second.end() && cursor.first->second < min_vo) {
            // Skip any runs that end before the window VO
            cursor.first++;
        }
        
        if (cursor.first != used_bins[cursor.second]->second.end()) {
            // If there are still runs in the bin, use them.
            cursor_queue.push(cursor);
        }
        
#ifdef debug
        cerr << "Bin " << used_bins[i]->first << " overlaps the query and has runs after the window min VO" << endl;
#endif
    }
    
    bool keep_going = true;
    
    while (keep_going && !cursor_queue.empty()) {
        // Loop until the user asks us to stop or we run out of things to give them.
    
#ifdef debug
        cerr << "Find earliest-starting run we haven't used yet in any bin" << endl;
#endif
            
        // Pull off the top element
        bin_cursor_t top = cursor_queue.top();
        cursor_queue.pop();
        
        // The bin windows are proper runs, so they can't overlap.
        // So after we deal with this run, we won't have to adjust any other bin cursors.
        
#ifdef debug
        cerr << "Found run " << top.first->first << "-" << top.first->second << endl;
#endif
        
        // Call the callback with the range max(min_vo from the window, that run's start) to that run's end.
        keep_going &= scan_callback(max(min_vo, top.first->first), top.first->second);
        
        if (!keep_going) {
            // The user is done with runs. They must have found a group that has an out-of-range minimum node ID.
            // We are done!
            return;
        }
        
        // Move up the min VO to the past the end of the run we just did.
        // TODO: We shouldn't need to do this since the runs can't overlap
        min_vo = top.first->second;
        
        // Advance what was the top iterator
        top.first++;
        if (top.first != used_bins[top.second]->second.end()) {
            // We haven't yet hit the end of this bin, so the cursor is eligible to be used again
            cursor_queue.push(top);
        } else {
#ifdef debug
            cerr << "Found bin index " << top.second << " is exhausted" << endl;
#endif
        }
    }
}

auto StreamIndexBase::scan_backward(const function<bool(int64_t, int64_t)> scan_callback) const -> void {
    // Remember the previous range's start VO, to be the next range's past-end VO.
    int64_t prev_vo = numeric_limits<int64_t>::max();
    
    for(auto rit = window_to_start.rbegin(); rit != window_to_start.rend(); ++rit) {
        
        // Go over the window offsets we have stored in reverse order.
        // We can use them as handy valid pointers to groups that are easy to go over in reverse order.
        
        if (!scan_callback(rit->second, prev_vo)) {
            // The iteratee is done
            return;
        }
        
        // Remember the start VO to be the next end
        prev_vo = rit->second;
    }
}

/// Return true if the given ID is in any of the sorted, coalesced, inclusive ranges in the vector, and false otherwise.
/// TODO: Is repeated binary search on the ranges going to be better than an unordered_set of all the individual IDs?
auto StreamIndexBase::is_in_range(const vector<pair<id_t, id_t>>& ranges, id_t id) -> bool {
    // Use a binary search
    size_t left = 0;
    size_t past_right = ranges.size();
    
    while (past_right >= left + 1) {
        // We have a nonempty interval
        
        // Find the middle
        size_t center = (left + past_right) / 2;
        assert(center < ranges.size());
        
        // Look at the range there
        auto& range = ranges[center];
        
        if (id < range.first) {
            // If we're before it, go left
            past_right = center;
        } else if (id > range.second) {
            // If we're after it, go right
            left = center + 1;
        } else {
            // If we're in it, return true
            return true;
        }
    }
    
    // If we get here, it wasn't in any range
    return false;
    
}

auto StreamIndexBase::save(ostream& to) const -> void {
    // We aren't going to save as Protobuf messages; we're going to save as a bunch of varints.
    
    // Format is
    // Magic bytes
    // Index version (varint32)
    // Bin count (varint64)
    // For each bin:
    // Bin number (varint64)
    // Run count (varint64)
    // For each run:
    // Start (varint64)
    // Past-end (varint64)
    // And then window count (varint64)
    // And for each window:
    // Window number (varint64)
    // Window start (varint64)
    
    // All the integers are Protobuf variable-length values.
    // The result is gzip-compressed.
    
    ::google::protobuf::io::OstreamOutputStream raw_out(&to);
    ::google::protobuf::io::GzipOutputStream gzip_out(&raw_out);
    ::google::protobuf::io::CodedOutputStream coded_out(&gzip_out);
    
    // Save the magic bytes
    coded_out.WriteRaw((void*)MAGIC_BYTES.c_str(), MAGIC_BYTES.size());
    
    // Save the version
    coded_out.WriteVarint32(OUTPUT_VERSION);
    
    // Save the bin count
    coded_out.WriteVarint64(bin_to_ranges.size());
    for (auto& kv : bin_to_ranges) {
        // For each bin, save the number
        coded_out.WriteVarint64(kv.first);
        // And the number of runs
        coded_out.WriteVarint64(kv.second.size());
        
        for (auto& run : kv.second) {
            // For each run, write the VO range
            coded_out.WriteVarint64(run.first);
            coded_out.WriteVarint64(run.second);
        }
    }
    
    // Save the window count
    coded_out.WriteVarint64(window_to_start.size());
    for (auto& kv : window_to_start) {
        // Save each window's number and start
        coded_out.WriteVarint64(kv.first);
        coded_out.WriteVarint64(kv.second);
    }
    
}

auto StreamIndexBase::load(istream& from) -> void {

    ::google::protobuf::io::IstreamInputStream raw_in(&from);
    ::google::protobuf::io::GzipInputStream gzip_in(&raw_in);
    

    bin_to_ranges.clear();
    window_to_start.clear();
    
    // Define an error handling function
    auto handle = [](bool ok) {
        if (!ok) throw std::runtime_error("GAMIndex::load detected corrupt index file");
    };
    
    // Look for the magic value
    
    // First read a bit of data
    char* buffer;
    int buffer_size = 0;
    while (buffer_size == 0) {
        // We must retry until we get some data, accoridng to the ZeroCopyInputStream spec
        handle(gzip_in.Next((const void**)&buffer, &buffer_size));
    }
    
    // TODO: In theory, we might have arbitrarily small buffers given to us.
    // We assume that the buffers are always big enough to actually peek the magic value and back up.
    assert(buffer_size >= MAGIC_BYTES.size());
    
    // We will fill this in with the version if we find it
    uint32_t input_version = 0;
    
    // Check to see if the magic bytes are there
    if (std::equal(MAGIC_BYTES.begin(), MAGIC_BYTES.end(), buffer)) {
        // We found the magic bytes! We know this is a versioned GAM index file.
        
        // Roll back to just after them
        gzip_in.BackUp(buffer_size - MAGIC_BYTES.size());
        
        // Read the input version
        {
            ::google::protobuf::io::CodedInputStream coded_in(&gzip_in);
            handle(coded_in.ReadVarint32(&input_version));
        }
    } else {
        // No magic bytes means input version 0
        // Roll back everything
        gzip_in.BackUp(buffer_size);
    }
    
    if (input_version > MAX_INPUT_VERSION) {
        throw std::runtime_error("GAMIndex::load can understand only up to index version " + to_string(MAX_INPUT_VERSION) +
            " and file is version " + to_string(input_version));
    }
    
    switch (input_version) {
    case 0:
    case 1:
        // Read the number of bins that are used
        uint64_t bin_count;
        {
            // TODO: To avoid hitting the coded input stream's byte limit (why is
            // it even at this level?) we destory and recreate it for every
            // semantic group.
            ::google::protobuf::io::CodedInputStream coded_in(&gzip_in);
            handle(coded_in.ReadVarint64(&bin_count));
        }
        
        for (size_t i = 0; i < bin_count; i++) {
            // Read the bin number and run count for each bin
            uint64_t bin_number;
            uint64_t run_count;
            {
                ::google::protobuf::io::CodedInputStream coded_in(&gzip_in);
                handle(coded_in.ReadVarint64(&bin_number));
                handle(coded_in.ReadVarint64(&run_count));
            }
            
            // Create the empty bin
            auto& runs = bin_to_ranges[bin_number];
            
            // Remember it in the bins by ID prefix index
            bins_by_id_prefix.insert(bin_to_prefix(bin_number), &runs);
            
            for (size_t j = 0; j < run_count; j++) {
                // Load each run
                uint64_t run_start;
                uint64_t run_end;
                
                {
                    ::google::protobuf::io::CodedInputStream coded_in(&gzip_in);
                    handle(coded_in.ReadVarint64(&run_start));
                    handle(coded_in.ReadVarint64(&run_end));
                }
                
                runs.emplace_back(run_start, run_end);
                
            }
            
        }
        
        
        // Now count the number of windows
        uint64_t window_count;
        {
            ::google::protobuf::io::CodedInputStream coded_in(&gzip_in);
            handle(coded_in.ReadVarint64(&window_count));
        }
        
        for (size_t i = 0; i < window_count; i++) {
            // Load each window
            uint64_t window_number;
            uint64_t window_start;
            
            {
                ::google::protobuf::io::CodedInputStream coded_in(&gzip_in);
                handle(coded_in.ReadVarint64(&window_number));
                handle(coded_in.ReadVarint64(&window_start));
            }
            
            window_to_start[window_number] = window_start;
        }
        break;
    default:
        throw std::runtime_error("Unimplemented GAM index version " + to_string(input_version));
    }

}

///////////
// Template specializations required for the available indexable types
///////////

template<>
auto StreamIndex<Alignment>::for_each_id(const Alignment& msg, const function<bool(const id_t&)> iteratee) const -> void {
    if (msg.path().mapping_size() == 0) {
        // This read is unmapped, so count it as node 0.
        iteratee(0);
    } else {
        // The read has mappings
        for (const auto& mapping : msg.path().mapping()) {
            // For each mapping
            if(!iteratee(mapping.position().node_id())) {
                // We showed it to the iteratee and it decided to stop.
                break;
            }
        }
    }
}

template<>
auto StreamIndex<Graph>::for_each_id(const Graph& msg, const function<bool(const id_t&)> iteratee) const -> void {
    
    // We use this to deduplicate our output.
    // TODO: We assume that node IDs are unique in the graph's nodes as an optimization.
    unordered_set<id_t> touched;
    
    for (auto& node : msg.node()) {
        // Show the iteratee all the nodes
        if (!iteratee(node.id())) {
            return;
        }
        
        // Make sure to count them as touched
        touched.insert(node.id());
    }
    
    for (auto& edge : msg.edge()) {
        // Show the iteratee the endpoints of all the edges
        
        if (!touched.count(edge.from())) {
            // The start isn't yet announced
            
            if (!iteratee(edge.from())) {
                // Stop early if asked
                return;
            }
            
            touched.insert(edge.from());
            
        }
        
        if (!touched.count(edge.to())) {
            // The end isn't yet announced
            
            if (!iteratee(edge.to())) {
                // Stop early if asked
                return;
            }
            
            touched.insert(edge.to());
        }
    }
    
    for (auto& path : msg.path()) {
        for (auto& mapping : path.mapping()) {
            // Show the iteratee all the nodes visited by paths
            const auto& id = mapping.position().node_id();
            
            if (id == 0) {
                // Skip magically unplaced mappings
                continue;
            }
            
            if (touched.count(id)) {
                // Skip things we already announced
                continue;
            }
            
            if (!iteratee(id)) {
                // Stop early if asked
                return;
            }
            
            touched.insert(id);
        }
        
    }
    
    if (touched.empty()) {
        // We didn't announce any nodes. Announce the unplaced sentinel.
        iteratee(0);
    }
    
}

}
