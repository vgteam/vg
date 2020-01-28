#include "gfa.hpp"
#include <gfakluge.hpp>

// Use sonLib pinch graphs
#include <stPinchGraphs.h>

#include <structures/union_find.hpp>

namespace vg {

using namespace std;
using namespace gfak;

// We augment the pinch graph API with functions that can work as if segments are in blocks even when they aren't.

/// Get the segment's orientation in its block (false=backward, true=forward), or true (forward) if there is no block
bool stPinchSegment_getBlockOrientationSafe(stPinchSegment* segment) {
    stPinchBlock* block = stPinchSegment_getBlock(segment);
    return (block != nullptr) ? stPinchSegment_getBlockOrientation(segment) : true;
}

/// Represents a translation from GFA node name string to pinch thread name number.
/// Tries to translate numerical node names as themselves, to the extent possible.
class GFAToPinchTranslator {
private:
    /// Map from string name to numerical name number
    unordered_map<string, int64_t> forward;
    /// Map from numerical name number back to string name
    unordered_map<int64_t, string> backward;
    /// What is the next unused name we can assign?
    int64_t next_unused = 1;
public:
    /// Translate from GFA name to pinch thread name
    int64_t translate(const string& name);
    /// Translate back from pinch thread name to GFA name
    const string& untranslate(const int64_t& name);
};

int64_t GFAToPinchTranslator::translate(const string& name) {
    // Look up the name
    auto found = forward.find(name);
    
    if (found != forward.end()) {
        // We have a translation already. Use it.
        return found->second;
    }
    
    // Otherwise we need to make a translation.
    // Try to find an unused number to use.
    // To start with, we have no clue what to use (0).
    int64_t assigned = 0;
    if (is_number(name)) {
        // There's a preferred number for this string. Try to use it.
        assigned = std::stol(name);
    }
    
    if (assigned <= 0 || backward.count(assigned)) {
         // We need to find an unused number.
         // next_unused is always guaranteed to be unused.
         assigned = next_unused;
         next_unused++;
    }
    
    if (assigned >= next_unused) {
        // If we read in a node ID from the GFA, this can happen.
        // Budge out the assignment cursor past any numbers yet mentioned in the GFA.
        // This is guaranteed to be past the largest assigned name, and therefore unused.
        next_unused = assigned + 1;
    }
    
    // Save the assigned numeric name
    forward[name] = assigned;
    // Mark it used and record the back translation
    backward[assigned] = name;
    // Return it
    return assigned;
}

const string& GFAToPinchTranslator::untranslate(const int64_t& name) {
    // If it was translated, it must have a reverse entry
    return backward.at(name);
}

/// Represents a translation from pinch thread segments' blocks to vg node IDs.
/// Tries to use pinch thread name numbers as IDs for single-thread blocks, and otherwise assigns IDs.
class PinchToVGTranslator {
private:
    /// Map from block or segment pointer to node ID
    unordered_map<void*, id_t> block_to_id;
    /// Track assigned numeric names
    unordered_set<id_t> used;
    /// What is the next unused name we can assign?
    id_t next_unused = 1;
public:
    /// Translate from pinch thread segment's block to node ID
    id_t translate(stPinchSegment* segment);
};

id_t PinchToVGTranslator::translate(stPinchSegment* segment) {
    // Work out what pointer will represent the segment's block
    void* block = (void*) stPinchSegment_getBlock(segment);
    if (block == nullptr) {
        // No block was found. The segment represents itself
        block = (void*) segment;
    }
    
    // Look up the block
    auto found = block_to_id.find(block);
    
    if (found != block_to_id.end()) {
        // We have a translation already. Use it.
        return found->second;
    }
    
    // Otherwise we need to make a translation. Try to find an unused number to
    // use. To start with, we will use the segment's thread number. If it is
    // the only segment in the block, and the block the only block for it, then
    // it is the right hint. Otherwise, it doesn't matter what we hint, because
    // we get to assign IDs arbitrarily.
    id_t assigned = (id_t) stPinchSegment_getName(segment);
    
    if (assigned <= 0 || used.count(assigned)) {
         // We need to find an unused number.
         // next_unused is always guaranteed to be unused.
         assigned = next_unused;
         next_unused++;
    }
    
    if (assigned >= next_unused) {
        // Budge out the assignment cursor past any numbers yet mentioned in the pinch names.
        // This is guaranteed to be past the largest assigned ID, and therefore unused.
        next_unused = assigned + 1;
    }
    
    // Save the assigned ID
    block_to_id[block] = assigned;
    // Mark it used
    used.insert(assigned);
    // Return it
    return assigned;
}

bool gfa_to_graph(istream& in, VG* graph, bool only_perfect_match) {

    // Things That Are Things
    // The GFA has "sequences" and "links" and "paths"
    // Each link's CIGAR is an alignment of the start of the second sequence to the end of the first
    // The pinch graph has "threads", "adjacencies", and "blocks"
    // The vg graph has "nodes" and "edges" and "paths" (again)
    
    // We have two layers of name/ID translation
    // One is from string GFA node name to pinch thread "name" number.
    // One is from pinch block pointer to vg graph node ID
    // For a GFA with no overlap and numeric names we want them to be invisible.
    // Otherwise in some places we will have to rework things.
    

    // New algorithm
    
    // If we are doing only perfect match:
    // Make a pinch thread set
    // Every sequence goes in
    // For each link
    // Discard the link if it is not all-M in any overlap it has
    // Pinch according to the link
    // Convert the pinch blocks to graph nodes
    // Sequence is always supplied by the first block member, whatever that is
    // Links with no overlap bypass the pinch graph entirely and become edges directly by fiddly lookup
    // Then create the paths by finding the each visited GFA sequence's thread in the pinch graph and translating via the blocks to vg space
    // When overlaps happen we know the overlap length in the second sequence (which arbitrarily loses)
    // And so we know what offset in the second sequence to start at
    // And if we didn't unify anything in the pinch graph we know a vg node will start/end there
    
    // If we are allowing non-perfect matches
    // Make the sequences into pinch graph paths the same way
    // Read the links' CIGAR strings and merge based on them
    // Except don't merge mismatched bases
    // Then resolve the sequences of pinch blocks the same way (since all of a block will be identical)
    // Then for the paths, when the path goes over an overlap edge, work out what offset in the next GFA sequence it ought to start at
    // And add a mapping to that vg node
    
    // So overall the algorithm is:
    // Make a pinch thread set
    // Every sequence goes in
    // For each link
    //  If doing perfect match mode, discard if not all M
    //  Compute CIGAR length in each sequence
    //  Pinch according to the CIGAR, only merging mismatches if in perfect match mode
    // Convert the pinch blocks to graph nodes
    // Sequence is always supplied by the first block member, whatever that is
    // Links with no overlap bypass the pinch graph entirely and become edges directly by fiddly lookup
    // Then create the paths by finding the each visited GFA sequence's thread in the pinch graph and translating via the blocks to vg space
    // When overlaps happen we know the overlap length in the second sequence from the link CIGAR
    // And so we know what offset in the second sequence to start at
    // If any link was discarded, also discard the path (and warn).
    
    // So let's start
    
    // Parse the GFA
    GFAKluge gg;
    gg.parse_gfa_file(in);
    // This maps from GFA sequence name to GFA sequence record
    map<string, sequence_elem, custom_key> gfa_sequences = gg.get_name_to_seq();
    // This maps from GFA sequence name to the GFA links for which it is the source
    map<string, vector<edge_elem> > gfa_links = gg.get_seq_to_edges();
    // This maps from path name to GFA path record
    map<string, path_elem> gfa_paths = gg.get_name_to_path();
    
    // Make a pinch thread set
    unique_ptr<stPinchThreadSet, decltype(&stPinchThreadSet_destruct)> pinch(stPinchThreadSet_construct(), &stPinchThreadSet_destruct);

    // Make a translator to convert from GFA string names to numeric pinch thread names
    GFAToPinchTranslator gfa_to_pinch;
    
    for(auto& name_and_record : gfa_sequences) {
        // For each GFA sequence record by string name
        auto& name = name_and_record.first;
        auto& record = name_and_record.second;
        
        // Assign it a numeric pinch thread name
        auto pinch_name = gfa_to_pinch.translate(name);
        
        // Add the thread to the pinch thread set
        stPinchThreadSet_addThread(pinch.get(), pinch_name, 0, record.sequence.size());
    }
    
    // As we go through the links, we need to remmeber links which are straight abutments of sequences.
    // These won't be processed by the pinch graph; we need to store them and process them later.
    // We store them in pinch thread name terms, as source, sink, source reverse, sink reverse.
    vector<tuple<int64_t, int64_t, bool, bool>> abut_links;
    
    // We also need to remember how much of the destination segment you should
    // skip due to overlap when reading across each edge in each direction. We
    // need this for interpreting GFA path records. We measure out to the end
    // of the alignment from the GFA edge, and don't "fix" things like terminal
    // insert/leading delete operations where the alignment really could just
    // be shorter. We keep entries for all non-discarded edges, even those with
    // no overlap.
    unordered_map<tuple<int64_t, int64_t, bool, bool>, size_t> link_skips;
    
    // Finally, we need to keep track of the dangling ends of alignments, which
    // need to be "tucked in". If an overlapping block starts or ends with
    // something other than a match, there will be a dangling node. It is in
    // correspondence with a position in the other sequence, but it doesn't get
    // pinched into it. We want to link it to the right next/previous node,
    // without actually pinching it in.
    
    // Because tucks can result in transative collapses and the creation of a
    // bunch of edges that would be hard to predict, we're going to use a
    // union-find over tucked sides and he sides they tuck into (i.e. have the
    // sme edges as). We can't make that union-find until we know everything
    // that needs to be involved, so we're going to store a list of pairs of
    // (thread name, base, is_left) tripples that get merged.
    using tuck_side_t = tuple<int64_t, size_t, bool>;
    vector<pair<tuck_side_t, tuck_side_t>> tucks;
    
    for (auto& name_and_links : gfa_links) {
        // For each set of links, by source node name
        auto& name = name_and_links.first;
        auto& links = name_and_links.second;
        
        // Find the source pinch thread
        auto source_pinch_name = gfa_to_pinch.translate(name);
        auto source_thread = stPinchThreadSet_getThread(pinch.get(), source_pinch_name);
        
        // And the source sequence string
        auto& source_string = gfa_sequences[name].sequence;
        // And sequence length
        auto source_sequence_length = stPinchThread_getLength(stPinchThreadSet_getThread(pinch.get(), source_pinch_name));
        
        for (auto& link : links) {
            // For each link on this source node
            
            // Get the CIGAR alignment for the link
            auto cigar = vcflib::splitCigar(link.alignment);
            
            if (only_perfect_match) {
                // We only care about all-match/mismatch CIGARs. Does this one qualify?
                bool is_all_match_mismatch = true;
                for (auto& elem : cigar) {
                    if (elem.second != "M" && elem.second != "=" && elem.second != "X") {
                        // This operation isn't any of the match/mismatch operations
                        is_all_match_mismatch = false;
                        break;
                    }
                }
                
                if (!is_all_match_mismatch) {
                    // This CIGAR has other operations in it and we want to discard the link because of it
                    continue;
                }
            }
            
            // Now we know we need to do the link and process the CIGAR.
            
            // Get the CIGAR's length in the source and the sink.
            // TODO: Is it really true that source is reference and sink is query?
            size_t source_alignment_length = 0;
            size_t sink_alignment_length = 0;
            for (auto& elem : cigar) {
                // Parse each CIGAR element
                auto& length = elem.first;
                assert(!elem.second.empty());
                switch(elem.second[0]) {
                case 'M':
                case '=':
                case 'X':
                    // Equal-length in both
                    source_alignment_length += length;
                    sink_alignment_length += length;
                    break;
                case 'I':
                    // Insert = more sink
                    sink_alignment_length += length;
                    break;
                case 'D':
                    // Deletion = more source
                    source_alignment_length += length;
                    break;
                case 'S':
                    // Soft clip = extra sink?
                    // Probably shouldn't be allowed.
                    throw runtime_error("GFA CIGAR contains a soft-clip operation; semantics unclear");
                    break;
                case 'H':
                    // Hard clip = extra sink also, but even weirder than a soft clip.
                    throw runtime_error("GFA CIGAR contains a hard-clip operation; semantics unclear");
                    break;
                default:
                    // This is an invalid operation; the GFA is invalid.
                    cerr << "error:[gfa_to_graph]: GFA CIGAR invalid: " << elem.second << " operation in " << link.to_string_2() << endl;
                    return false;
                }
            }
            
            // Work out what thread the link is to
            auto sink_pinch_name = gfa_to_pinch.translate(link.sink_name);
            auto sink_thread = stPinchThreadSet_getThread(pinch.get(), sink_pinch_name);
            
            // Get the orientations
            bool source_backward = !link.source_orientation_forward;
            bool sink_backward = !link.sink_orientation_forward;
            
            // Record the link skip distances from the alignment, for interpreting paths.
            // When traversion source to sink, skip the alignment's length in the sink
            link_skips[make_tuple(source_pinch_name, sink_pinch_name, source_backward, sink_backward)] = sink_alignment_length;
            // When traversing sink to source, skip the alignment's length in the source.
            link_skips[make_tuple(sink_pinch_name, source_pinch_name, !sink_backward, !source_backward)] = source_alignment_length;
            
#ifdef debug
            cerr << "Found edge " << link.source_name << " = " << source_pinch_name << (source_backward ? 'L' : 'R')
                << " -> " << link.sink_name << " = " << sink_pinch_name << (sink_backward ? 'R' : 'L') << endl;
            cerr << "Skips: " << sink_alignment_length << " forward, " << source_alignment_length << " reverse" << endl; 
#endif
            
            if (source_alignment_length == 0 || sink_alignment_length == 0) {
                // This link is just an end-to-end abutment with no overlap.
                // It can't be sent through the pinch graph; we store it separately.
                // We also have no need to do any tucks.
                abut_links.push_back(make_tuple(source_pinch_name, sink_pinch_name, source_backward, sink_backward));
                
                // TODO: What exactly are the semantics on e.g. an all-insert
                // alignment? Is it really just that the ends abut? Or do we
                // connect in to the end of the insert?
                
                // Skip the link CIGAR execution
                continue;
            }
            
            // Find the sink sequence data
            auto& sink_string = gfa_sequences[link.sink_name].sequence;
            // And sequence length
            auto sink_sequence_length = stPinchThread_getLength(stPinchThreadSet_getThread(pinch.get(), sink_pinch_name));
            
            if (source_alignment_length > source_sequence_length) {
                // The GFA file is invalid and specifies using more bases than are present.
                if (only_perfect_match) {
                    // Be tolerant and reject just this edge.
                    continue;
                }
                
                // Otherwise complain to the user
                cerr << "error:[gfa_to_graph]: GFA file contains a link " << link.source_name << " " << (source_backward ? 'L' : 'R')
                    << " to " << link.sink_name  << " " << (sink_backward ? 'R' : 'L')
                    << " that tries to consume more source sequence (" << source_alignment_length
                    << ") than is present (" << source_sequence_length << ")" << endl;
                    
                return false;
            }
            
            if (sink_alignment_length > sink_sequence_length) {
                // The GFA file is invalid and specifies using more bases than are present.
                if (only_perfect_match) {
                    // Be tolerant and reject just this edge.
                    continue;
                }
                
                // Otherwise complain to the user
                cerr << "error:[gfa_to_graph]: GFA file contains a link " << link.source_name << " " << (source_backward ? 'L' : 'R')
                    << " to " << link.sink_name  << " " << (sink_backward ? 'R' : 'L')
                    << " that tries to consume more sink sequence (" << sink_alignment_length
                    << ") than is present (" << sink_sequence_length << ")" << endl;
                    
                return false;
            }
            
            // Set up some cursors in each node's sequence that go the right
            // direction, based on orientations. Cursors start at the first
            // base in the CIGAR, which may be past the end/before the
            // beginning on the source if the CIGAR is 0 length.
            int64_t source_cursor = source_backward ? (-1 + source_alignment_length) : (source_sequence_length - source_alignment_length);
            int64_t source_motion = source_backward ? -1 : 1;
            int64_t sink_cursor = sink_backward ? (sink_sequence_length - 1) : 0;
            int64_t sink_motion = sink_backward ? -1 : 1;
            
            // Decide if we are pinching in agreeing orientations
            bool pinch_same_strand = (source_backward == sink_backward);
            
            // Interpret the CIGAR string and perform pinches.
            
            for (size_t cigar_index = 0; cigar_index < cigar.size(); cigar_index++) {
                // For each cigar operation
                auto& elem = cigar[cigar_index];
                
                if (elem.first == 0) {
                    // Skip 0-length operations
                    continue;
                }
                
                if (cigar_index != 0 && (elem.second == "D" && cigar[cigar_index - 1].second == "I" ||
                    elem.second == "I" && cigar[cigar_index - 1].second == "D")) {
                    // We found adjacent inserts and deletions, which throws
                    // off our dangling tip tucking algorithm and which
                    // shouldn't happen anyway in a well-behaved alignment.
                    // TODO: accomodate this somehow.
                    throw runtime_error("GFA importer cannot (yet) handle adjacent insertions and deletions.");
                }
                
                // Decompose each operation into a series of suboperations.
                // This gives us an opportunity to replace M operations with =/X
                vector<pair<size_t, string>> suboperations;
                
                if (elem.second == "M" && !only_perfect_match) {
                    // This is an M operation that needs to be decomposed into = and X
                    for (size_t i = 0; i < elem.first; i++) {
                        // For each position along the M operation
                        
                        // Find the next character in the source
                        auto source_char = source_string.at(source_cursor + source_motion * i);
                        if (source_backward) {
                            source_char = reverse_complement(source_char);
                        }
                        // And the sink
                        auto sink_char = sink_string.at(sink_cursor + sink_motion * i);
                        if (sink_backward) {
                            sink_char = reverse_complement(sink_char);
                        }
                        // Work out what kind of operation we need for this pairing of bases.
                        // TODO: Handle Ns specially?
                        string opcode = (source_char == sink_char) ? "=" : "X";
                        
                        if (!suboperations.empty() && suboperations.back().second == opcode) {
                            // We can accumulate onto the existing suboperation of this type
                            suboperations.back().first++;
                        } else {
                            // We need a new suboperation of this type
                            suboperations.push_back(make_pair(1, opcode));
                        }
                    }
                } else {
                    // This operation can be passed through as-is
                    suboperations.push_back(elem);
                }
               
                for (auto subelem_index = 0; subelem_index < suboperations.size(); subelem_index++) {
                    // For each suboperation
                    auto& subelem = suboperations[subelem_index];
                    
                    // Get its length
                    auto& length = subelem.first;
                    
                    // Compute source and sink lengths
                    size_t source_length = 0;
                    size_t sink_length = 0;
                    assert(!subelem.second.empty());
                    switch(subelem.second[0]) {
                    case 'M':
                    case '=':
                    case 'X':
                        source_length = length;
                        // Fall through
                    case 'I':
                        sink_length = length;
                        break;
                    case 'D':
                        source_length = length;
                        break;
                    default:
                        // We should have already checked for weird operations.
                        throw runtime_error("Invalid operation " + subelem.second + " in pre-screened CIGAR");
                    }
                    
                    // Work out the sequence-local start of the region in each sequence that it may apply to, which depends on orientation.
                    int64_t source_region_start = source_backward ? (source_cursor - source_length + 1) : source_cursor;
                    int64_t sink_region_start = sink_backward ? (sink_cursor - sink_length + 1) : sink_cursor;
                    
                    // And also the end positions (last in region, on the other side of the start if the region is empty)
                    int64_t source_region_last = source_backward ? source_cursor : (source_cursor + source_length - 1);
                    int64_t sink_region_last = sink_backward ? sink_cursor : (sink_cursor + sink_length - 1);
                    
                    // Note thatb these are just lowest/highest coordinate in
                    // thread space, not corresponding to each other to to the
                    // start/end of the operation.
                    
#ifdef debug
                    cerr << "Suboperation " << subelem.first << subelem.second << " runs "
                        << source_region_start << " through " << source_region_last << " in " << source_pinch_name
                        << " and " << sink_region_start << " through " << sink_region_last << " in " << sink_pinch_name << endl;
#endif
                    
                    // We need to know when to wire in dangling source/sink bits. 
                    bool is_first_subelement = (cigar_index == 0 && subelem_index == 0);
                    bool is_last_subelement = (cigar_index == cigar.size() - 1 && subelem_index == suboperations.size() - 1);
                    
                    assert(!subelem.second.empty());
                    switch(subelem.second[0]) {
                    case 'M':
                        if (only_perfect_match) {
                            // The whole match can be merged
                            stPinchThread_pinch(source_thread, sink_thread, source_region_start, sink_region_start, length, pinch_same_strand);
                        } else {
                            // If we aren't in always_perfect_match mode this should have become =/X
                            throw runtime_error("Encountered unparsed M operation");
                        }
                        break;
                    case '=':
                        // Always pinch.
                        // TODO: should we check sequence equality?
                        stPinchThread_pinch(source_thread, sink_thread, source_region_start, sink_region_start, length, pinch_same_strand);
                        break;
                    case 'X':
                        // Only pinch if we are forcing matches (in which case this was X originally)
                        if (only_perfect_match) {
                            stPinchThread_pinch(source_thread, sink_thread, source_region_start, sink_region_start, length, pinch_same_strand);
                        } else {
                            if (is_first_subelement) {
                                // We dangled the sink and potentially created
                                // a tip. We need to remember the pinch we
                                // would have made, so we can wire up the
                                // dangling end to what we would have wired it
                                // to if we had pinched it.
                                
                                // Describe the end that is dangling as (thread name, base, is_left)
                                tuck_side_t dangling_end = make_tuple(sink_pinch_name,
                                    sink_backward ? stPinchThread_getLength(sink_thread) - 1 : 0, !sink_backward);
                                // Describe what it would have merged with.
                                tuck_side_t tuck_into = make_tuple(source_pinch_name, source_cursor, !source_backward);
                                // Queue up the edge exchange
                                tucks.emplace_back(dangling_end, tuck_into);
                                
#ifdef debug    
                                cerr << "\tSuboperation requires tucking "
                                    << get<0>(dangling_end) << ":" << get<1>(dangling_end) << (get<2>(dangling_end) ? 'L' : 'R')
                                    << " in with " <<  get<0>(tuck_into) << ":" << get<1>(tuck_into) << (get<2>(tuck_into) ? 'L' : 'R') << endl;
#endif
                                
                                // Splits to allow tucks will be done later.
                            }
                            
                            if (is_last_subelement) {
                                // We dangled the source as well, on the right of the overlap
                                
                                // Describe the end that is dangling as (thread name, base, is_left)
                                tuck_side_t dangling_end = make_tuple(source_pinch_name,
                                    source_backward ? 0 : stPinchThread_getLength(source_thread) - 1, source_backward);
                                // Describe what it would have merged with.
                                tuck_side_t tuck_into = make_tuple(sink_pinch_name, sink_cursor + sink_motion * (sink_length - 1), sink_backward);
                                // Queue up the edge exchange
                                tucks.emplace_back(dangling_end, tuck_into);
                                
#ifdef debug    
                                cerr << "\tSuboperation requires tucking "
                                    << get<0>(dangling_end) << ":" << get<1>(dangling_end) << (get<2>(dangling_end) ? 'L' : 'R')
                                    << " in with " <<  get<0>(tuck_into) << ":" << get<1>(tuck_into) << (get<2>(tuck_into) ? 'L' : 'R') << endl;
#endif
                                
                                // Splits to allow tucks will be done later.
                            }
                        }
                        break;
                    case 'I':
                        // We don't need to do any pinching
                        if (is_first_subelement) {
                            // We dangled the sink on the left of the overlap
                            
                            // Describe the end that is dangling as (thread name, base, is_left)
                            tuck_side_t dangling_end = make_tuple(sink_pinch_name,
                                sink_backward ? stPinchThread_getLength(sink_thread) - 1 : 0, !sink_backward);
                            // Describe what it would have merged with.
                            tuck_side_t tuck_into = make_tuple(source_pinch_name, source_cursor, !source_backward);
                            // Queue up the edge exchange
                            tucks.emplace_back(dangling_end, tuck_into);
                            
#ifdef debug    
                            cerr << "\tSuboperation requires tucking "
                                << get<0>(dangling_end) << ":" << get<1>(dangling_end) << (get<2>(dangling_end) ? 'L' : 'R')
                                << " in with " <<  get<0>(tuck_into) << ":" << get<1>(tuck_into) << (get<2>(tuck_into) ? 'L' : 'R') << endl;
#endif
                            
                            // Splits to allow tucks will be done later.
                        }
                        break;
                    case 'D':
                        // No pinching
                        if (is_last_subelement) {
                            // We dangled the source on the right of the overlap
                            
                            // Describe the end that is dangling as (thread name, base, is_left)
                            tuck_side_t dangling_end = make_tuple(source_pinch_name,
                                source_backward ? 0 : stPinchThread_getLength(source_thread) - 1, source_backward);
                            // Describe what it would have merged with.
                            // We have to back out sink motion since the cursor is right now at the after-the-deletion position.
                            tuck_side_t tuck_into = make_tuple(sink_pinch_name, sink_cursor - sink_motion, sink_backward);
                            // Queue up the edge exchange
                            tucks.emplace_back(dangling_end, tuck_into);
                            
#ifdef debug    
                            cerr << "\tSuboperation requires tucking "
                                << get<0>(dangling_end) << ":" << get<1>(dangling_end) << (get<2>(dangling_end) ? 'L' : 'R')
                                << " in with " <<  get<0>(tuck_into) << ":" << get<1>(tuck_into) << (get<2>(tuck_into) ? 'L' : 'R') << endl;
#endif
                            
                        }
                        break;
                    default:
                        // We should have already checked for weird operations twice now.
                        throw runtime_error("Invalid operation " + subelem.second + " in pre-screened CIGAR");
                    }
                    
                    // Advance the cursors
                    sink_cursor += sink_motion * sink_length;
                    source_cursor += source_motion * source_length;
                }
            }
        }
    }

    // Now all the pinches have been made
    
    // Also we know all the tucks that need to be union-finded.
    
    // Make a set of them
    unordered_set<tuck_side_t> tuck_side_set;
    for (auto& to_union : tucks) {
        tuck_side_set.insert(to_union.first);
        tuck_side_set.insert(to_union.second);
    }
    
    for (auto& side : tuck_side_set) {
        // Don't forget to split the threads to expose the sides we need to copy edges between
        auto& thread_name = get<0>(side);
        auto& position = get<1>(side);
        auto& is_left = get<2>(side);
        
#ifdef debug
        cerr << "Need to expose thread " << thread_name << ":" << position << " side " << (is_left ? 'L' : 'R') << endl;
#endif
        
        // Find the thread
        stPinchThread* thread = stPinchThreadSet_getThread(pinch.get(), thread_name);
        
        if (is_left && position == 0) {
            // Don't try and ask to split at -1; this side is already exposed
            continue;
        }
        
        // Split the thread to expose this side
        stPinchThread_split(thread, is_left ? position - 1 : position);
    }
    
    // Everything else should wait until after graph nodes are made
    
#ifdef debug
    {
        size_t segment_count = 0;
        auto segment_iter = stPinchThreadSet_getSegmentIt(pinch.get());
        stPinchSegment* segment = stPinchThreadSetSegmentIt_getNext(&segment_iter);
        for (; segment != nullptr; segment = stPinchThreadSetSegmentIt_getNext(&segment_iter)) {
            segment_count++;
        }
        
        cerr << "Total pinch segments: " << segment_count << endl;
    }
#endif

    // Convert the pinch blocks into vg nodes
    
    // We use this translator to translate block pointers to node IDs
    PinchToVGTranslator pinch_to_vg;
    
    {
        auto segment_iter = stPinchThreadSet_getSegmentIt(pinch.get());
        stPinchSegment* segment = stPinchThreadSetSegmentIt_getNext(&segment_iter);
        for (; segment != nullptr; segment = stPinchThreadSetSegmentIt_getNext(&segment_iter)) {
            
            // For each segment in the pinch set (including those not in blocks)
            
            // Generate or assign the node ID for its block or itself
            id_t node_id = pinch_to_vg.translate(segment);
            
            if (graph->has_node(node_id)) {
                // We already did this graph node, from another segment in the block.
                continue;
            }
            
            // Get the segment's thread name, + strand offset, length, and orientation
            auto segment_thread_name = stPinchSegment_getName(segment);
            auto segment_start = stPinchSegment_getStart(segment);
            auto segment_length = stPinchSegment_getLength(segment);
            auto segment_backward = !stPinchSegment_getBlockOrientationSafe(segment);
            
            // Go find the source DNA for the GFA sequence that created the thread
            const string& thread_sequence = gfa_sequences[gfa_to_pinch.untranslate(segment_thread_name)].sequence;
            
            // Compute the sequence for the vg node we will make
            string node_sequence = thread_sequence.substr(segment_start, segment_length);
            if (segment_backward) {
                node_sequence = reverse_complement(node_sequence);
            }
            
            // Make the node in the graph
            graph->create_node(node_sequence, node_id);
        }
    }
    
    {
        // Add edges from pinch graph adjacencies
        auto thread_iter = stPinchThreadSet_getIt(pinch.get());
        stPinchThread* thread = stPinchThreadSetIt_getNext(&thread_iter);
        for (; thread != nullptr; thread = stPinchThreadSetIt_getNext(&thread_iter)) {
        
            // For each thread in the pinch thread set
            // Start at the beginning
            stPinchSegment* here = stPinchThread_getFirst(thread);
            assert(here != nullptr);
            
            // Look one segment ahead
            stPinchSegment* next = stPinchSegment_get3Prime(here);
            
            while (next != nullptr) {
                // For each pinch graph connection from here to next
                
                // Get the nodes we are connecting
                id_t here_id = pinch_to_vg.translate(here);
                id_t next_id = pinch_to_vg.translate(next);
                
                // Get their orientations
                bool here_backward = !stPinchSegment_getBlockOrientationSafe(here);
                bool next_backward = !stPinchSegment_getBlockOrientationSafe(next);
                
                // Make the edge if not present already
                Edge* e = graph->create_edge(here_id, next_id, here_backward, next_backward);
                
#ifdef debug
                cerr << "Created pinch graph edge " << pb2json(*e) << endl;
#endif
                
                // Advance right
                here = next;
                next = stPinchSegment_get3Prime(here);
            }
        }
    }
    
    // Add edges from abut_links
    for (auto& abutment : abut_links) {
        // Unpack each abutment record
        auto& source_name = get<0>(abutment);
        auto& sink_name = get<1>(abutment);
        auto& source_backward = get<2>(abutment);
        auto& sink_backward = get<3>(abutment);
        
        // Get the threads by name
        stPinchThread* source_thread = stPinchThreadSet_getThread(pinch.get(), source_name);
        stPinchThread* sink_thread = stPinchThreadSet_getThread(pinch.get(), sink_name);
        
        // Find the segment of the source that is relevant.
        // If the source sequence is forward, it is the last one, but if it is backward, it is the first one.
        stPinchSegment* source_segment = source_backward ? stPinchThread_getFirst(source_thread) : stPinchThread_getLast(source_thread);
        // And conversely for the sink
        stPinchSegment* sink_segment = sink_backward ? stPinchThread_getLast(sink_thread) : stPinchThread_getFirst(sink_thread);
    
        // Get the node IDs to connect
        id_t from_id = pinch_to_vg.translate(source_segment);
        id_t to_id = pinch_to_vg.translate(sink_segment);
        
        // Figure out the orientation of each node. We take whether the segemnt
        // is backward in its node, and flip it if the segment itself is
        // visited backward.
        bool from_start = (!stPinchSegment_getBlockOrientationSafe(source_segment) != source_backward);
        bool to_end = (!stPinchSegment_getBlockOrientationSafe(sink_segment) != sink_backward);
        
        // Make the edge
        Edge* e = graph->create_edge(from_id, to_id, from_start, to_end);
        
#ifdef debug
        cerr << "Created abutment edge " << pb2json(*e) << endl;
#endif
    }
    
    // Copy edges to equivalent nodes caused by tucking
    
#ifdef debug
    cerr << "Before tucks:" << endl;
    cerr << pb2json(graph->graph) << endl;
#endif


    // OK now we need to do the tucking
    
    // Make a vector of tuck sides to give them all indexes
    vector<tuck_side_t> tuck_side_vector(tuck_side_set.begin(), tuck_side_set.end());
    tuck_side_set.clear();
    
    // Get a map to invert the vector
    unordered_map<tuck_side_t, size_t> tuck_side_to_index;
    for (size_t i = 0; i < tuck_side_vector.size(); i++) {
        tuck_side_to_index[tuck_side_vector[i]] = i;
    }
    
    // Create the union-find
    structures::UnionFind tuck_merger(tuck_side_vector.size());
    
    for (auto& to_union : tucks) {
        // Do all the unioning
        auto side1 = tuck_side_to_index[to_union.first];
        auto side2 = tuck_side_to_index[to_union.second];
        tuck_merger.union_groups(side1, side2);
    }
    // We don't need these anymore now that we have the union-find filled in.
    tuck_side_to_index.clear();
    tucks.clear();
    
    // Convert all the tuck sides to outward-pointing handles
    vector<handle_t> tuck_handles(tuck_side_vector.size());
    // Also invert the handle vector so we can get handle indexes
    unordered_map<handle_t, size_t> handle_to_index;
    for (size_t i = 0; i < tuck_side_vector.size(); i++) {
        // Unpack each side
        auto& side = tuck_side_vector[i];
        auto& thread_name = get<0>(side);
        auto& position = get<1>(side);
        auto& is_left = get<2>(side);
        
        // Find the thread
        stPinchThread* thread = stPinchThreadSet_getThread(pinch.get(), thread_name);
        
        // Find the segment
        stPinchSegment* segment = stPinchThread_getSegment(thread, position);
        
        // Find the node and segment-relative orientation
        id_t node = pinch_to_vg.translate(segment);
        bool segment_reverse_in_node = !stPinchSegment_getBlockOrientationSafe(segment);
        
        // Make the handle
        tuck_handles[i] = graph->get_handle(node, is_left != segment_reverse_in_node);
        // And store its index under it
        handle_to_index[tuck_handles[i]] = i;
    }
    tuck_side_vector.clear();
    
    // Keep a set, for each union-find group, of all the destination handles
    // that everything in the group has to have an edge to.
    unordered_map<size_t, unordered_set<handle_t>> destinations;
    
    for (size_t i = 0; i < tuck_handles.size(); i++) {
        // For each tuck side handle
        auto& handle = tuck_handles[i];
        
        // Find the destination set for the union-find group this tuck side handle belongs to
        auto& my_destinations = destinations[tuck_merger.find_group(i)];
        
        graph->follow_edges(handle, false, [&](const handle_t& next) {
            // Loop over every handle it goes to (facing inward)
        
            // Add each handle to the destinations set for the union-find group that owns this tuck side
            my_destinations.insert(next);
            
            // Check and see if this place we go itself belongs to a union-find
            // group (making sure to make it face us first)
            auto found = handle_to_index.find(graph->flip(next));
            
            if (found != handle_to_index.end()) {
                // It does, so we have to pretend it is also all the other things in its group
                for(size_t other_index : tuck_merger.group(found->second)) {
                    // Go get all the handles that tuck in with it
                    handle_t& also_tucked = tuck_handles.at(other_index);
                
                    // And add them all as destinations too
                    my_destinations.insert(graph->flip(also_tucked));
                }
            }
        });
    }
    
    for (size_t i = 0; i < tuck_handles.size(); i++) {
        // For each tuck side handle again
        auto& handle = tuck_handles[i];
        
        // Look up all the destinations for its union-find group
        auto& my_destinations = destinations[tuck_merger.find_group(i)];
        
        for (auto& destination : my_destinations) {
            // Connect it to each of them
            graph->create_edge(handle, destination);
        }
        
    }
    
    // Now all the nodes and edges exist.
    
    // Process the GFA paths
    
    for (auto& name_and_path : gfa_paths) {
        // For each path record by name
        auto& name = name_and_path.first;
        auto& path = name_and_path.second;
        
#ifdef debug
        cerr << "Import path " << name << endl;
        cerr << path.to_string() << endl;
#endif
        
        // Create each path
        graph->paths.create_path(name);
        
        if (path.segment_names.size() == 0) {
            // Empty paths need nothing else.
            continue;
        }
        
        // Start the path with visits to the entirety of the first thread it traces
        {
            // Find the thread to visit
            int64_t thread_name = gfa_to_pinch.translate(path.segment_names[0]);
            // Determine if it is visited backward
            bool thread_backward = !path.orientations[0];
            
            // Get the actual thread
            stPinchThread* thread = stPinchThreadSet_getThread(pinch.get(), thread_name);
            
            // Get the starting end appropriate to the orientation
            stPinchSegment* segment = thread_backward ? stPinchThread_getLast(thread) : stPinchThread_getFirst(thread);
            
#ifdef debug
            cerr << "\tBegin at " << path.segment_names[0]
                << " = " << thread_name << (thread_backward ? 'L' : 'R') << endl;
#endif
            
            while (segment != nullptr) {
                // Look up the node
                id_t node = pinch_to_vg.translate(segment);
                // Compute its visit orientation
                bool node_backward = (!stPinchSegment_getBlockOrientationSafe(segment) != thread_backward);
                
#ifdef debug 
                cerr << "\t\tPath starts with " << stPinchSegment_getLength(segment)
                    << " bases on node " << node << " orientation " << (node_backward ? "rev" : "fwd") << endl;
#endif
                
                // Visit it
                graph->paths.append_mapping(name, node, node_backward, stPinchSegment_getLength(segment), 0);
                // Advance to the next segment
                segment = thread_backward ? stPinchSegment_get5Prime(segment) : stPinchSegment_get3Prime(segment);
            }
        }
        
        // If we find a nonexistent/skipped link we need to abort the entire path.
        bool abort_path = false;
        for (size_t i = 1; i < path.segment_names.size() && !abort_path; i++) {
            // For each subsequent GFA path visit (which becomes a thread)
            
            // Find the thread to visit
            int64_t thread_name = gfa_to_pinch.translate(path.segment_names[i]);
            // Determine if it is visited backward
            bool thread_backward = !path.orientations[i];
            
            // And the previous thread
            int64_t prev_thread_name = gfa_to_pinch.translate(path.segment_names[i - 1]);
            bool prev_thread_backward = !path.orientations[i - 1];
            
#ifdef debug
            cerr << "\tCross edge " << path.segment_names[i - 1]
                << " = " << prev_thread_name << (prev_thread_backward ? 'L' : 'R') 
                << " to " << path.segment_names[i] << " = " << thread_name << (thread_backward ? 'R' : 'L') << endl;
#endif
            
            // Work out how much of this thread the previous thread ate, by looking at the overlaps on the links
            auto overlap_to = link_skips.find(tie(prev_thread_name, thread_name, prev_thread_backward, thread_backward));
            if (overlap_to == link_skips.end()) {
                // This thread crosses a link that isn't there.
                // We want to get rid of the path entirely since we can't represent it.
                
                if (only_perfect_match) {
                    // Edge may have been removed for having a bad alignment.
                    cerr << "warning [gfa_to_graph]: path " << name << ": edge " << path.segment_names[i - 1]
                        << " = " << prev_thread_name << (prev_thread_backward ? 'L' : 'R') 
                        << " to " << path.segment_names[i] << " = " << thread_name << (thread_backward ? 'R' : 'L')
                        << " is not present. It may have been removed due to having a bad alignment. Discarding path!" << endl;
                } else {
                    // All the edges should be there. This is an error in the GFA.
                    stringstream msg;
                    msg << "error [gfa_to_graph]: path " << name << ": edge " << path.segment_names[i - 1]
                        << " = " << prev_thread_name << (prev_thread_backward ? 'L' : 'R') 
                        << " to " << path.segment_names[i] << " = " << thread_name << (thread_backward ? 'R' : 'L')
                        << " is not present. The GFA file is malformed!";
                    throw runtime_error(msg.str());
                }
                
                // Remove the path and skip out on adding the rest of it
                graph->paths.remove_path(name);
                abort_path = true;
                break;
            }
            
            // TODO: We don't check for *vg graph edges* that aren't present, only *links*.
        
            // Start at the near end of the next thread
            stPinchThread* thread = stPinchThreadSet_getThread(pinch.get(), thread_name);
            stPinchSegment* segment = thread_backward ? stPinchThread_getLast(thread) : stPinchThread_getFirst(thread);
            
            // Skip segments until we have accounted for the overlap.
            size_t overlap_skipped = 0;
            while (overlap_skipped < overlap_to->second && segment != nullptr) {
                overlap_skipped += stPinchSegment_getLength(segment);
#ifdef debug 
                cerr << "\t\tSkip overlap of " << stPinchSegment_getLength(segment)
                    << " from segment for node " << pinch_to_vg.translate(segment) << endl;
#endif
                segment = thread_backward ? stPinchSegment_get5Prime(segment) : stPinchSegment_get3Prime(segment);
            }
            
            // We should always reach the overlap at a segment boundary
            assert(overlap_skipped == overlap_to->second);
            
            // Continue adding segments to the path as nodes from there until the end of the thread.
            while (segment != nullptr) {
                // Look up the node
                id_t node = pinch_to_vg.translate(segment);
                // Compute its visit orientation
                bool node_backward = (!stPinchSegment_getBlockOrientationSafe(segment) != thread_backward);
                
#ifdef debug 
                cerr << "\t\tPath follows " << stPinchSegment_getLength(segment)
                    << " bases on node " << node << " orientation " << (node_backward ? "rev" : "fwd") << endl;
#endif
                
                // Visit it
                graph->paths.append_mapping(name, node, node_backward, stPinchSegment_getLength(segment), 0);
                // Advance to the next segment
                segment = thread_backward ? stPinchSegment_get5Prime(segment) : stPinchSegment_get3Prime(segment);
            }
        }
    }
    
    // Save the paths to the graph
    graph->paths.rebuild_mapping_aux();
    graph->paths.to_graph(graph->graph);
    
    // Now the graph is done!
    // TODO: validate graph and paths and assign path mapping ranks
    
    return true;
    
}

void graph_to_gfa(const VG* graph, ostream& out) {
  GFAKluge gg;
  gg.set_version(1.0);
  for (auto h : gg.get_header()){
    out << h.second.to_string();
  }

    // TODO moving to GFAKluge
    // problem: protobuf longs don't easily go to strings....
    
    graph->for_each_node([&](const Node* n) {
        sequence_elem s_elem;
        // Fill seq element for a node
        s_elem.name = to_string(n->id());
        s_elem.sequence = n->sequence();
        out << s_elem.to_string_1() << endl;
        //gg.add_sequence(s_elem);
    });
    
    auto& pathmap = graph->paths._paths;
    for (auto p : pathmap){
        path_elem p_elem;
        p_elem.name = p.first;
        for (auto m : p.second){
            p_elem.segment_names.push_back( std::to_string(m.node_id()) );
            p_elem.orientations.push_back( !m.is_reverse() );
            const Node* n = graph->get_node( m.node_id() );
            stringstream cigaro;
            //cigaro << n->sequence().size() << (p.mapping(m_ind.position().is_reverse()) ? "M" : "M");
            cigaro << n->sequence().size() << (m.is_reverse() ? "M" : "M");
            p_elem.overlaps.push_back( cigaro.str() );
        }
        out << p_elem.to_string_1() << endl;
        //gg.add_path(p_elem.name, p_elem);
    }

    graph->for_each_edge([&](const Edge* e) {
        edge_elem ee;
        ee.type = 1;
        ee.source_name = to_string(e->from());
        ee.sink_name = to_string(e->to());
        ee.source_orientation_forward = ! e->from_start();
        ee.sink_orientation_forward =  ! e->to_end();
        ee.alignment = std::to_string(e->overlap()) + "M";
        
        if (e->from_start() && (e->to_end() || e->to() < e->from())) {
            // Canonicalize edges to be + orientation first if possible, and
            // then low-ID to high-ID if possible, for testability. This edge
            // needs to flip.
            
            // Swap the nodes
            std::swap(ee.source_name, ee.sink_name);
            // Swap the orientations
            std::swap(ee.source_orientation_forward, ee.sink_orientation_forward);
            // Reverse the orientations
            ee.source_orientation_forward = !ee.source_orientation_forward;
            ee.sink_orientation_forward = !ee.sink_orientation_forward;
        }
        
        out << ee.to_string_1() << endl;;
        //gg.add_edge(ee.source_name, ee);
        //link_elem l;
        //l.source_name = to_string(e->from());
        //l.sink_name = to_string(e->to());
        //l.source_orientation_forward = ! e->from_start();
        //l.sink_orientation_forward =  ! e->to_end();
        //l.cigar = std::to_string(e->overlap()) + "M";
        //gg.add_link(l.source_name, l);
    });
    //gg.output_to_stream(cout);
}

}
