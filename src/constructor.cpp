/*
 * constructor.cpp: contains implementations for vg construction functions.
 */


#include "vg.hpp"
#include "constructor.hpp"

#include <set>
#include <tuple>
#include <list>
#include <algorithm>
#include <memory>

namespace vg {

using namespace std;

vcflib::Variant* VcfBuffer::get() {
    if(has_buffer) {
        // We have a variant
        return &buffer;
    } else {
        // No next variant loaded.
        return nullptr;
    }
}

void VcfBuffer::handle_buffer() {
    // The variant in the buffer has been dealt with
    assert(has_buffer);
    has_buffer = false;
}

void VcfBuffer::fill_buffer() {
    if(file != nullptr && file->is_open() && !has_buffer && safe_to_get) {
        // Put a new variant in the buffer if we have a file and the buffer was empty.
        has_buffer = safe_to_get = file->getNextVariant(buffer);
        if(has_buffer) {
            // Convert to 0-based positions.
            // TODO: refactor to use vcflib zeroBasedPosition()...
            buffer.position -= 1;
        }
#ifdef debug
        cerr << "Variant in buffer: " << buffer << endl;
#endif
    }
}

bool VcfBuffer::has_tabix() {
    return file && file->usingTabix;
}

bool VcfBuffer::set_region(const string& contig, int64_t start, int64_t end) {
    if(!has_tabix()) {
        // Nothing to seek in (or no file)
        return false;
    }
    
    if(start != -1 && end != -1) {
        // We have a start and end
        return file->setRegion(contig, start, end);
    } else {
        // Just seek to the whole chromosome
        return file->setRegion(contig);
    }
}

VcfBuffer::VcfBuffer(vcflib::VariantCallFile* file) : file(file) {
    // Our buffer needs to know about the VCF file it is reading from, because
    // it cares about the sample names. If it's not associated properely, we
    // can't getNextVariant into it.
    if (file) {
        // But only do it if we actually have a real file
        buffer.setVariantCallFile(file);
    }
}

ConstructedChunk Constructor::construct_chunk(string reference_sequence, string reference_path_name,
    vector<vcflib::Variant> variants, size_t chunk_offset) const {
    
    // Construct a chunk for this sequence with these variants.
    ConstructedChunk to_return;
    
    // We need as path to add reference mappings to
    Path* ref_path = to_return.graph.add_path();
    ref_path->set_name(reference_path_name);
    
    // We use this to keep track of what the next unused base, if any, in the
    // reference is.
    size_t reference_cursor = 0;
    
    // We use this to number nodes. All chunks number nodes starting from 1.
    id_t next_id = 1;
    
    // We remember nodes ending at these reference positions that haven't yet
    // all been wired up yet. These are on-the-end and not past-the-end
    // positions, so they are start + length - 1.
    map<size_t, set<id_t>> nodes_ending_at;
    // And nodes starting at these reference positions that haven't yet all been
    // wired up. 
    map<size_t, set<id_t>> nodes_starting_at;
    
    // We don't want to wire inserts to each other, so we have a set of all
    // insert endpoints.
    set<id_t> inserts;
    
    // Here we remember deletions that end at paritcular positions in the
    // reference, which are the positions of the last deleted bases. We map from
    // last deleted base to last non-deleted base before the deletion, so we can
    // go look up nodes ending there. Note that these can map to -1.
    map<size_t, set<int64_t>> deletions_ending_at;
    
    // We also need to track all points at which deletions start, so we can
    // search for the next one when deciding where to break the reference.
    set<int64_t> deletion_starts;
    
    // We use this to get the next variant
    auto next_variant = variants.begin();
    
    // We're going to clump overlapping variants together.
    vector<vcflib::Variant*> clump;
    // And we need to rember the highest past-the-end base of anything in the
    // clump, to catch all the overlaps.
    size_t clump_end = 0;
    
    
    // We have a utility function to tack a full length perfect match onto a
    // path. We need the node so we can get its length.
    auto add_match = [](Path* path, Node* node) {
        // Make a mapping for it
        auto* mapping = path->add_mapping();
        mapping->mutable_position()->set_node_id(node->id());
        // Make it a perfect match explicitly
        auto* match_edit = mapping->add_edit();
        match_edit->set_from_length(node->sequence().size());
        match_edit->set_to_length(node->sequence().size());
    };
    
    // Given a string, turn it into nodes of the max node size or smaller, and
    // add them to the graph. Return pointers to the resulting nodes in the
    // graph, in order.
    auto create_nodes = [&](const string& sequence) -> vector<Node*> {
        
        // How big should each node try to be?
        size_t piece_size;
        
        if(greedy_pieces) {
            // Make pieces as big as possible.
            piece_size = max_node_size;
        } else {
        
            // Let's try to divide into evenly-sized pieces
            size_t piece_count = sequence.size() / max_node_size + 1;
            piece_size = sequence.size() / piece_count;
            // Remember we may have a partial piece at the end.
            
            // TODO: we're rounding down and tend to have 1-base tiny pieces at
            // the end.
        }
    
        // We'll fill this in with created nodes
        vector<Node*> created;
        
        // We keep a cursor to the next non-made-into-a-node base
        size_t cursor = 0;
        
        while (cursor < sequence.size()) {
            // There's still sequence to do, so bite off a piece
            size_t next_node_size = std::min(piece_size, sequence.size() - cursor);
            string node_sequence = sequence.substr(cursor, next_node_size);
            
            // Make a node
            auto* node = to_return.graph.add_node();
            node->set_id(next_id++);
            node->set_sequence(node_sequence);
            
            // Tack it on the end of the list
            created.push_back(node);
            
            // Advance the cursor since we made this node
            cursor += next_node_size;
        }
        
        return created;
    };
    
    // We have a function to emit reference nodes from wherever the current
    // cursor position is up to the given position, advancing the cursor. The
    // target position must be <= the length of the reference. This function
    // adds the nodes to the starting and ending position indexes.
    auto add_reference_nodes_until = [&](size_t target_position) {
        
        // Make new nodes for all the sequence we want to add
        auto new_nodes = create_nodes(reference_sequence.substr(reference_cursor, target_position - reference_cursor));
    
        // Remember the total node length we have scanned through
        size_t seen_bases = 0;
    
        for (Node* node : new_nodes) {
            // Add matches on the reference path for all the new nodes
            add_match(ref_path, node);
            
            // Remember where it starts and ends along the reference path
            nodes_starting_at[reference_cursor + seen_bases].insert(node->id());
            nodes_ending_at[reference_cursor + seen_bases + node->sequence().size() - 1].insert(node->id());
            
            // Remember how long that node was so we place the next one right.
            seen_bases += node->sequence().size();
        }
    
        // Advance the cursor
        reference_cursor = target_position;
    };
    
    while (next_variant != variants.end() || !clump.empty()) {
        // While there are more variants, or while we have the last clump to do...
    
        // Group variants into clumps of overlapping variants.
        if (clump.empty() || 
            (next_variant != variants.end() && clump_end > next_variant->position - chunk_offset)) {
            
            // Either there are no variants in the clump, or this variant
            // overlaps the clump. It belongs in the clump
            clump.push_back(&(*next_variant));
            // It may make the clump longer and necessitate adding more variants.
            clump_end = max(clump_end, next_variant->position + next_variant->ref.size() - chunk_offset);
            
            // Try the variant after that
            next_variant++;
        } else {
            // The next variant doesn't belong in this clump.
            // Handle the clump.
            
            // Parse all the variants into VariantAllele edits
            
            // This holds a map from Variant pointer to a vector of lists
            // of VariantAllele edits, one list per non-ref allele of the
            // variant.
            map<vcflib::Variant*, vector<list<vcflib::VariantAllele>>> parsed_clump;
            
            // This determines the order we will process variants in. We use it
            // to sort variants in a clump by hash for the purposes of assigning
            // IDs.
            map<string, vcflib::Variant*> variants_by_name;
            
            // This holds the min and max values for starts and ends of edits
            // not removed from the clump.
            int64_t first_edit_start = numeric_limits<int64_t>::max();
            int64_t last_edit_end = -1;
            
            for (vcflib::Variant* variant : clump) {
            
                // Check the variant's reference sequence to catch bad VCF/FASTA pairings
                auto expected_ref = reference_sequence.substr(variant->position - chunk_offset, variant->ref.size());
            
                if(variant->ref != expected_ref) {
                    // TODO: report error to caller somehow
                    cerr << "error:[vg::Constructor] Variant/reference sequence mismatch: " << variant->ref
                        << " vs " << expected_ref << "; do your VCF and FASTA coordinates match?"<< endl;
                    cerr << "Variant: " << *variant << endl;
                    exit(1);
                }
            
                // Name the variant and place it in the order that we'll
                // actually construct nodes in (see utility.hpp)
                string variant_name = make_variant_id(*variant);
                variants_by_name[variant_name] = variant;
            
                // We need to parse the variant into alts, each of which is a
                // series of VariantAllele edits. This holds the full alt allele
                // string and the edits needed to make it. The VariantAlleles
                // completely cover the alt, and some of them may be perfect
                // matches to stretches of reference sequence. Note that the
                // reference allele of the variant won't appear here.
                map<string, vector<vcflib::VariantAllele>> alternates = flat ? variant->flatAlternates() :
                    variant->parsedAlternates();

                for (auto& kv : alternates) {                
                    // For each alt in the variant
                    
                    if (kv.first == variant->ref) {
                        // Skip the ref, because we can't make any ref nodes
                        // until all the edits for the clump are known.
                        continue;
                    }
                    
                    // With 0 being the first non-ref allele, which alt are we?
                    // Copy the string out of the map
                    string alt_string = kv.first;
                    // Then look it up
                    size_t alt_index = variant->getAltAlleleIndex(alt_string);
                    
                    if (alt_index >= parsed_clump[variant].size()) {
                        // Make sure we have enough room to store the VariantAlleles for this alt.
                        parsed_clump[variant].resize(alt_index + 1);
                    }
                    
                    // Find the list of edits for this alt
                    auto& alt_parts = parsed_clump[variant][alt_index];
                    
#ifdef debug
                    cerr << "Non-ref allele " << alt_index << endl;
#endif
                    
                    // Copy all the VariantAlleles into the list
                    alt_parts.assign(kv.second.begin(), kv.second.end());
                    
                    while (!alt_parts.empty() && alt_parts.front().ref == alt_parts.front().alt) {
                        // Drop leading ref matches
#ifdef debug
                        cerr << "\tDrop " << alt_parts.front().ref << " -> " << alt_parts.front().alt
                            << " @ " << alt_parts.front().position - chunk_offset << endl;
#endif
                        alt_parts.pop_front();
                    }
                    
                    while (!alt_parts.empty() && alt_parts.back().ref == alt_parts.back().alt) {
                        // Drop trailing ref matches
#ifdef debug
                        cerr << "\tDrop " << alt_parts.back().ref << " -> " << alt_parts.back().alt
                            << " @ " << alt_parts.back().position - chunk_offset << endl; 
#endif
                        alt_parts.pop_back();
                    }
                    
#ifdef debug
                    for (auto& edit : alt_parts) {
                        cerr << "\tKept " << edit.ref << " -> " << edit.alt << " @ " << edit.position - chunk_offset << endl; 
                    }
#endif
                    
                    if (!alt_parts.empty()) {
                        // If this alt's interior non-ref portion exists, see if
                        // it extends beyond those of other alts in the clump.
                        first_edit_start = min(first_edit_start, (int64_t) (alt_parts.front().position - chunk_offset));
                        last_edit_end = max(last_edit_end,
                            (int64_t) (alt_parts.back().position - chunk_offset + alt_parts.back().ref.size() - 1));
                    }
                }
            }
            
            // We have to have some non-ref material, even if it occupies 0
            // reference space.
            assert(last_edit_end != -1);
            assert(first_edit_start != numeric_limits<int64_t>::max());
            
#ifdef debug
            cerr << "Edits run between " << first_edit_start << " and " << last_edit_end << endl;
#endif
            
            // Create ref nodes from the end of the last clump (where the cursor
            // is) to the start of this clump's interior non-ref content.
            add_reference_nodes_until(first_edit_start);
        
            // This keeps track of edits that already have nodes, consisting of
            // a ref position, a ref sequence, and an alt sequence. It maps to a
            // vector of pointers to the nodes created, which are owned by the
            // graph. The vector is to handle the case where an edit is too long
            // to be all one node, according to our max node length, and is
            // always nonempty.
            map<tuple<long, string, string>, vector<Node*>> created_nodes;
        
            // This holds on to variant ref paths, which we can't actually fill
            // in until all the variants in the clump have had their non-ref
            // paths done.
            map<vcflib::Variant*, Path*> variant_ref_paths;
        
            for (auto& kv : variants_by_name) {
                // For each variant in the clump, sorted by name
                auto& variant_name = kv.first;
                auto* variant = kv.second;
            
                // Declare its ref path straight away.
                // We fill inthe ref paths after we make all the nodes for the edits.
                variant_ref_paths[variant] = to_return.graph.add_path();
                variant_ref_paths[variant]->set_name("_alt_" + variant_name + "_0");
            
                for (size_t alt_index = 0; alt_index < parsed_clump[variant].size(); alt_index++) {                
                    // For each non-ref alt in the parsed variant
                    
                    // Name the alt after the number that this allele has.
                    // We have to bump the allele index because the first alt is 0.
                    string alt_name = "_alt_" + variant_name + "_" + to_string(alt_index + 1);
                    
                    // There should be a path named after it.
                    Path* alt_path = to_return.graph.add_path();
                    alt_path->set_name(alt_name);
                    
                    for (vcflib::VariantAllele& edit : parsed_clump[variant][alt_index]) {
                        // For each VariantAllele used by the alt
                        
#ifdef debug
                        cerr << "Apply " << edit.ref << " -> " << edit.alt << " @ " << edit.position << endl;
#endif
                        
                        if (edit.alt != "") {
                            // This is a visit to a node for the alt
                        
                            // We need a key to see if a node has been made for this edit already
                            auto key = make_tuple(edit.position - chunk_offset, edit.ref, edit.alt);
                            
                            if (created_nodes.count(key) == 0) {
                                // We don't have a run of nodes for this edit, so make one.
                                vector<Node*> node_run = create_nodes(edit.alt);
                                
                                // Compute where the edit starts and ends in local chunk coordinates
                                auto edit_start = edit.position - chunk_offset;
                                auto edit_end = edit.position - chunk_offset + edit.ref.size() - 1;
                                
#ifdef debug
                                cerr << "Created nodes running " << edit_start << " to " << edit_end << endl;
#endif
                                
                                // Remember where the first one starts and the last one ends, for wiring up later.
                                nodes_starting_at[edit_start].insert(node_run.front()->id());
                                nodes_ending_at[edit_end].insert(node_run.back()->id());
                                
                                // Save it in case any other alts also have this edit.
                                created_nodes[key] = node_run;
                                
                                if (edit.ref == "") {
                                    // This is an insert, so mark its ends as
                                    // such, so they don't connect to other
                                    // insert ends.
                                    inserts.insert(node_run.front()->id());
                                    inserts.insert(node_run.back()->id());
                                    
#ifdef debug
                                    cerr << "Nodes are insert" << endl;
#endif
                                }
                            } else {
#ifdef debug
                                cerr << "Found existing nodes" << endl;
#endif
                            }
                            
                            
                            for (Node* node : created_nodes[key]) {
                                // Add a visit to each node we created/found in
                                // order to the path for this alt of this
                                // variant.
                                add_match(alt_path, node);
                            }
                            
                        } else if (edit.ref != "") {
                            // It's a deletion (and not a weird ""->"" edit).
                            
                            // Add an entry to the deletion arcs
                            
                            // What is the past-the-end position (first non-deleted)
                            size_t arc_end = edit.position - chunk_offset + edit.ref.size();
                            // What is the before-the-beginning position (last non-deleted, may be -1)
                            int64_t arc_start = (int64_t) edit.position - chunk_offset - 1;
                            
#ifdef debug
                            cerr << "Ensure deletion arc " << arc_start << " to " << arc_end << endl;
#endif
                            
                            // Add the arc (if it doesn't exist). We only index
                            // arcs from the end, because we'll make them when
                            // looking for predecessors of nodes. TODO: could we
                            // make special handling of deletions that go to -1
                            // more efficient?
                            deletions_ending_at[arc_end].insert(arc_start);
                            
                            // Remember that an arc comes from this base
                            deletion_starts.insert(arc_start);
                        }
                        
                    }
                    
                }
            }
            
            // Then after you finish all the alts, add the ref nodes that
            // haven't already been created, breaking wherever something can
            // come in or out.
            
            // We need a function to work that out
            auto next_breakpoint_after = [&](size_t position) -> size_t {
                // This returns the position of the base to the left of the next
                // required breakpoint within this clump, after the given
                // position, given created nodes and deletions that already
                // exist.
                
                // If nothing else, we're going to break at the end of the last
                // edit in the clump.
                size_t to_return = last_edit_end;
                
                // See if any nodes are registered as starting after our
                // position. They'll all start before the end of the clump, and
                // we don't care if they start at our position since that
                // breakpoint already happened.
                auto next_start_iter = nodes_starting_at.upper_bound(position);
                
                if(next_start_iter != nodes_starting_at.end()) {
                    // If we found something, walk back where the breakpoint
                    // needs to be so we break before that node starts.
                    to_return = min(to_return, next_start_iter->first - 1);
                }
                
                // See if any nodes are registered as ending at or after our
                // position. We do care if they end at our position, since that
                // means we need to break right here.
                auto next_end_iter = nodes_ending_at.lower_bound(position);
                
                if(next_end_iter != nodes_ending_at.end()) {
                    // If we found something, we need to break where that node
                    // ends.
                    to_return = min(to_return, next_end_iter->first );
                }
                
                // See if any deletions are registered as ending after here.
                // Deletions break the reference before their past-the-end base,
                // so we don't care about deletions ending here exactly.
                auto deletion_end_iter = deletions_ending_at.upper_bound(position);
                
                if(deletion_end_iter != deletions_ending_at.end()) {
                    // If we found something, walk back where the breakpoint
                    // needs to be so we break before the node after the
                    // deletion starts.
                    to_return = min(to_return, deletion_end_iter->first - 1);
                }
                
                // See if any deletions are known to start at or after this
                // base. We care about exact hits now, because deletions break
                // after the base they start at.
                auto deletion_start_iter = deletion_starts.lower_bound(position);
                // We don't need to worry about -1s here. They won't be found
                // with lower_bound on a size_t.
                
                if(deletion_start_iter != deletion_starts.end()) {
                    // If we found something, walk back where the breakpoint
                    // needs to be so we break at the position the deletion
                    // needs to leave from.
                    to_return = min(to_return, (size_t)*deletion_start_iter);
                }
                
                return to_return;
                
            };
            
            // Note that in some cases (i.e. pure inserts) we may not need any ref nodes at all.
            
            while (reference_cursor < last_edit_end + 1) {
                // Until we hot the end
                
                // Find where the next node run must end to attach to stuff
                size_t next_end = next_breakpoint_after(reference_cursor);
            
                // We need to have a reference node/run of nodes (which may have
                // already been created by a reference match) between where the
                // cursor is and where the next breakpoint has to be.
                // This is the sequence it should have.
                string run_sequence = reference_sequence.substr(reference_cursor, next_end - reference_cursor + 1);
                
                // We need a key to see if a node (run) has been made for this sequece already
                auto key = make_tuple(reference_cursor, run_sequence, run_sequence);
                
                if (created_nodes.count(key) == 0) {
                    // We don't have a run of ref nodes up to the next break, so make one
                    vector<Node*> node_run = create_nodes(run_sequence);
                    
                    // Remember where the first one starts and the last one ends, for wiring up later.
                    nodes_starting_at[reference_cursor].insert(node_run.front()->id());
                    nodes_ending_at[next_end].insert(node_run.back()->id());
                    
#ifdef debug
                    cerr << "Created reference nodes running " << reference_cursor << " to " << next_end << endl;
#endif
                    
                    // Save it in case any other alts also have this edit.
                    created_nodes[key] = node_run;
                }
                
                for (Node* node : created_nodes[key]) {
                    // Add a reference visit to each node we created/found
                    add_match(ref_path, node);
                    
                    for (vcflib::Variant* variant : clump) {
                        // For each variant we might also be part of the ref allele of
                        if (reference_cursor >= variant->position - chunk_offset &&
                            reference_cursor < variant->position - chunk_offset + variant->ref.size()) {
                            // If this run of nodes starts within the varaint...
                            // (We know if it starts in the variant it has to
                            // end in the variant, because the variant ends with
                            // a node break)
                            
                            // Add a match along the variant's ref allele path
                            add_match(variant_ref_paths[variant], node);
                        }
                    }
                }
                
                // Advance the reference cursor to after this run of reference nodes
                reference_cursor = next_end + 1;
                
                // Keep going until we have created reference nodes through to
                // the end of the clump's interior edits.
            }
            
            // Now we have gotten through all the places where nodes start, before the end of the clump.
            
            // Now the clump is handled
            clump.clear();
            clump_end = 0;
            // On the next loop we'll grab the next variant for the next clump.
        }
    }

    // Create reference path nodes and mappings after the last clump.
    add_reference_nodes_until(reference_sequence.size());
    
    
    // Create all the edges
    for (auto& kv : nodes_starting_at) {
        if (kv.first == 0) {
            // These are the nodes that abut the left edge of the chunk. Add
            // each of these nodes to the set of left-edge-abuting nodes.
            for (auto& node_id : kv.second) {
                to_return.left_ends.insert(node_id);
            }
        } else {
            // These are nodes that start somewhere else.
            for (auto& right_node : kv.second) {
                // For every node that could occur here
            
                for (auto& left_node : nodes_ending_at[kv.first - 1]) {
                    // For every node that could come before these nodes
                    
                    if (inserts.count(left_node) && inserts.count(right_node)) {
                        // Don't connect two inserts at the same position (or an insert to itself).
#ifdef debug
                        cerr << "Skip insert-insert edge " << left_node << " -> " << right_node << endl;
#endif
                        continue;
                    }
                    
#ifdef debug
                    cerr << "Add normal edge " << left_node << " -> " << right_node << endl;
#endif
                    
                    // Emit an edge
                    auto* edge = to_return.graph.add_edge();
                    edge->set_from(left_node);
                    edge->set_to(right_node);
                }
                
                for (auto& deletion_start : deletions_ending_at[kv.first]) {
                    // For everywhere a deletion can start that comes to here
                    
                    if (deletion_start == -1) {
                        // This node ought to be exposed on the left actually.
                        to_return.left_ends.insert(right_node);
                    
                    } else {
                        // The deletion doesn't go all the way to the left edge
                        // but actually starts at a place where there are nodes.
                        
                        for (auto& left_node : nodes_ending_at[deletion_start]) {
                            // For every node that the deletion could start with
                            
#ifdef debug
                            cerr << "Add deletion edge " << left_node << " -> " << right_node << endl;
#endif
                            
                            // Emit an edge
                            auto* edge = to_return.graph.add_edge();
                            edge->set_from(left_node);
                            edge->set_to(right_node);
                        }
                    }
                }
            }
        }
    }
    
    for(auto& node_id : nodes_ending_at[reference_sequence.size() - 1]) {
        // Add each node that ends at the end of the chunk to the set of such nodes
        to_return.right_ends.insert(node_id);
    }
    
    for(auto& deletion_start : deletions_ending_at[reference_sequence.size()]) {
        // Also add in nodes at the starts of deletions that go to the end of the chunk
        
        if(deletion_start == -1) {
            // Note that we don't handle completely spanning deletions. But
            // those can't be articulated in VCF anyway because alts can't be
            // empty.
            continue;
        }
        
        for (auto& node_id : nodes_ending_at[deletion_start]) {
            // For every node that the deletion could start with
            // Expose it on the right of the graph
            to_return.right_ends.insert(node_id);
        }
    }
    
    // Remember to tell the caller how many IDs we used
    to_return.max_id = next_id - 1;
    
    return to_return;
}

void Constructor::add_name_mapping(const string& vcf_name, const string& fasta_name) {
    // Fill in both one-way maps.
    // TODO: C++ doesn't have a 2-way map right?
    vcf_to_fasta_renames[vcf_name] = fasta_name;
    fasta_to_vcf_renames[fasta_name] = vcf_name;
}

string Constructor::vcf_to_fasta(const string& vcf_name) const {
    return vcf_to_fasta_renames.count(vcf_name) ? vcf_to_fasta_renames.at(vcf_name) : vcf_name;
}

string Constructor::fasta_to_vcf(const string& fasta_name) const {
    return fasta_to_vcf_renames.count(fasta_name) ? fasta_to_vcf_renames.at(fasta_name) : fasta_name;
}

void Constructor::construct_graph(string vcf_contig, FastaReference& reference, VcfBuffer& variant_source,
    function<void(Graph&)> callback) {
    
    
    // Our caller will set up indexing. We just work with the buffered source that we pull variants from.

    // What sequence are we looking for in the fasta? The one we were passed, unless it was renamed.
    string reference_contig = vcf_to_fasta(vcf_contig);

    // At what offset in the reference sequence do we start?
    size_t leading_offset;
    // At what position in the reference sequence do we stop (past-the-end)?
    size_t reference_end;
    
    if (allowed_vcf_regions.count(vcf_contig)) {
        // Only look at the region we were asked for. We will only find variants
        // *completely* contained in this region! Partially-overlapping variants
        // will be discarded!
        leading_offset = allowed_vcf_regions[vcf_contig].first;
        reference_end = allowed_vcf_regions[vcf_contig].second;
    } else {
        // Look at the whole contig
        leading_offset = 0;
        reference_end = reference.sequenceLength(reference_contig);
    }
    
    // Scan through variants until we find one that is on this contig and in this region.
    // If we're using an index, we ought to already be at the right place.
    variant_source.fill_buffer();
    while(variant_source.get() && (variant_source.get()->sequenceName != vcf_contig ||
        variant_source.get()->position < leading_offset ||
        variant_source.get()->position + variant_source.get()->ref.size() > reference_end)) {
        // This variant comes before our region
        
        // Discard variants that come out that are before our region
        variant_source.handle_buffer();
        variant_source.fill_buffer();
    }
    
    // Now we're on the variants we actually want.

    // This is where the next chunk will start in the reference sequence.
    size_t chunk_start = leading_offset;
    
    // We maintain a growing list of variants that will go into a chunk. They
    // are all positioned relative to chunk_start.
    vector<vcflib::Variant> chunk_variants;
    // And we track the largest past-the-end position of all the variants
    size_t chunk_end = 0;
    
    // For chunk wiring, we need to remember the nodes exposed on the end of the
    // previous chunk.
    set<id_t> exposed_nodes;
    
    // We'll also need to bump chunk IDs out of the way. What's the max ID used
    // in previous chunks?
    id_t max_id = 0;
    
    // Whenever a chunk ends with a single node, we separate it out and buffer
    // it here, because we may need to glue it together with subsequent leading
    // nodes that were broken by a chunk boundary.
    Node last_node_buffer;
    
    // Sometimes we need to emit single node reference chunks gluing things
    // together
    auto emit_reference_node = [&](Node& node) {
    
        // Make a single node chunk for the node
        Graph chunk;
        *(chunk.add_node()) = node;
        
        // It needs a primary path mapping.
        Path* path = chunk.add_path();
        path->set_name(reference_contig);
        Mapping* mapping = path->add_mapping();
        mapping->mutable_position()->set_node_id(node.id());
        Edit* edit = mapping->add_edit();
        edit->set_from_length(node.sequence().size());
        edit->set_to_length(node.sequence().size());
        
        // Emit this chunk we were holding back.
        callback(chunk);
    };
    
    // When a chunk gets constructed, we'll call this handler, which will wire
    // it up to the previous chunk, if any, and then call the callback we're
    // supposed to send our graphs out through.
    // Modifies the chunk in place.
    auto wire_and_emit = [&](ConstructedChunk& chunk) {
        // When each chunk comes back:
        
        if (chunk.left_ends.size() == 1 && last_node_buffer.id() != 0) {
            // We have a last node from the last chunk that we want to glom onto
            // this chunk.
            
            // Grab the first node, which we know must be the single source for
            // the chunk.
            Node* mutable_first_node = chunk.graph.mutable_node(0);
            assert(chunk.left_ends.count(mutable_first_node->id()) == 1);
            
            // Combine the sequences for the two nodes
            string combined_sequence = last_node_buffer.sequence() + mutable_first_node->sequence();
            
            if (combined_sequence.size() <= max_node_size) {
                // We can fit both nodes into one node.
                mutable_first_node->set_sequence(combined_sequence);
                
                // We can re-use the ID from the last node, which we discard.
                // Edges to it will get rerouted to the first node. And we know
                // it can't have any mappings except the primary path.
                max_id--;
                
                // We don't need any edges to it, either.
                exposed_nodes.clear();
                
                // Clear the buffer since we moved its sequence and ID into the
                // graph.
                last_node_buffer = Node();
            } else {
                // We need to keep two nodes. Reapportion the sequence between
                // them according to our division algorithm. TODO: can sometimes
                // differ from old construct behavior, but this way will be
                // better.
                size_t piece_size = greedy_pieces ? max_node_size : ((combined_sequence.size() + 1) / 2);
                last_node_buffer.set_sequence(combined_sequence.substr(0, piece_size));
                mutable_first_node->set_sequence(combined_sequence.substr(piece_size));
                
                // Emit the buffered node as a chunk
                emit_reference_node(last_node_buffer);
                // Clear it
                last_node_buffer = Node();
            }

            // Update the mapping lengths on the mutable first node.
            // First we find the primary path
            Path* path = chunk.graph.mutable_path(0);
            assert(path->name() == reference_contig);
            // Then the first mapping
            Mapping* mapping = path->mutable_mapping(0);
            assert(mapping->position().node_id() == mutable_first_node->id());
            assert(mapping->edit_size() == 1);
            // Then the only edit
            Edit* edit = mapping->mutable_edit(0);
            // Correct its length
            edit->set_from_length(mutable_first_node->sequence().size());
            edit->set_to_length(mutable_first_node->sequence().size());
        } else if (last_node_buffer.id() != 0) {
            // There's no single leading node on this next chunk, but we still
            // have a single trailing node to emit.
            
            // Emit it
            emit_reference_node(last_node_buffer);
            // Clear it
            last_node_buffer = Node();
        }
        
        if (chunk.right_ends.size() == 1) {
            // We need to pull out the last node in the chunk. Note that it may
            // also be the first node in the chunk...
            
            // We know it's the last node in the graph
            last_node_buffer = chunk.graph.node(chunk.graph.node_size() - 1);
            
            
            assert(chunk.right_ends.count(last_node_buffer.id()));
            
            // Remove it
            chunk.graph.mutable_node()->RemoveLast();
            
            // Find the primary path
            Path* path = chunk.graph.mutable_path(0);
            assert(path->name() == reference_contig);
            // Then drop last mapping, which has to be to this node
            assert(path->mapping_size() > 0);
            assert(path->mapping(path->mapping_size() - 1).position().node_id() == last_node_buffer.id());
            path->mutable_mapping()->RemoveLast();
            
            // Update its ID separately, since it's no longer in the graph.
            last_node_buffer.set_id(last_node_buffer.id() + max_id);
        }
        
        // Up all the IDs in the graph
        // TODO: this is repeating code that vg::VG has...
        for (size_t i = 0; i < chunk.graph.node_size(); i++) {
            // For each node
            auto* node = chunk.graph.mutable_node(i);
            // Bump the node ID
            node->set_id(node->id() + max_id);
        }
        for (size_t i = 0; i < chunk.graph.edge_size(); i++) {
            // For each edge
            auto* edge = chunk.graph.mutable_edge(i);
            // Bump the edge end IDs
            edge->set_from(edge->from() + max_id);
            edge->set_to(edge->to() + max_id);
        }
        for (size_t i = 0; i < chunk.graph.path_size(); i++) {
            // For each path
            auto* path = chunk.graph.mutable_path(i);
            for (size_t j = 0; j < path->mapping_size(); j++) {
                // For each mapping in the path
                auto* mapping = path->mutable_mapping(j);
                
                // Bump the ID for the mapping's position
                mapping->mutable_position()->set_node_id(mapping->position().node_id() + max_id);
            }
        }
        
        // If there was a previous ConstructedChunk, wire up the edges between them
        for (auto& from_id : exposed_nodes) {
            // For every dangling end in the last chunk
            
            for (auto& to_id : chunk.left_ends) {
                // For every node in the new chunk we can wire it to
                
                // Make the edge in the new chunk
                Edge* new_edge = chunk.graph.add_edge();
                new_edge->set_from(from_id);
                // Make sure to correct the number in the to set.
                new_edge->set_to(to_id + max_id);
            }
        }
        
        // Save the right-side ends from this chunk for the next one, if any
        exposed_nodes.clear();
        for (auto& from_id : chunk.right_ends) {
            // Make sure to correct each ID
            exposed_nodes.insert(from_id + max_id);
        }

        // Remember the new max id, accounting for all the IDs used by this
        // chunk.
        max_id += chunk.max_id;
        
        // Emit the chunk's graph via the callback
        callback(chunk.graph);
    };
    
    while (variant_source.get() && variant_source.get()->sequenceName == vcf_contig &&
        variant_source.get()->position >= leading_offset &&
        variant_source.get()->position + variant_source.get()->ref.size() <= reference_end) {
    
        // While we have variants we want to include
        
        bool variant_acceptable = true;
        for (string& alt : variant_source.get()->alt) {
            // Validate each alt of the variant
            if(!allATGC(alt)) {
                // It may be a symbolic allele or something. Skip this variant.
                cerr << "warning:[vg::Constructor] Unsupported variant allele \"" << alt << "\"; Skipping variant!" << endl;
                variant_acceptable = false;
                break;
            }
        }
        if (!variant_acceptable) {
            // Skip variants that have symbolic alleles or other nonsense we can't parse.
            variant_source.handle_buffer();
            variant_source.fill_buffer();
        } else if (!chunk_variants.empty() && chunk_end > variant_source.get()->position) {
            // If the chunk is nonempty and this variant overlaps what's in there, put it in too and try the next.
            // TODO: this is a lot like the clumping code...
            
            // Add it in
            chunk_variants.push_back(*(variant_source.get()));
            // Expand out how big the chunk needs to be, so we can get other overlapping variants.
            chunk_end = max(chunk_end, chunk_variants.back().position + chunk_variants.back().ref.size());
            
            // Try the next variant
            variant_source.handle_buffer();
            variant_source.fill_buffer();
            
        } else if(chunk_variants.size() < vars_per_chunk && variant_source.get()->position < chunk_start + bases_per_chunk) {
            // Otherwise if this variant is close enough and the chunk isn't too big yet, put it in and try the next.
            
            // TODO: unify with above code?
            
            // Add it in
            chunk_variants.push_back(*(variant_source.get()));
            // Expand out how big the chunk needs to be, so we can get other overlapping variants.
            chunk_end = max(chunk_end, chunk_variants.back().position + chunk_variants.back().ref.size());
            
            // Try the next variant
            variant_source.handle_buffer();
            variant_source.fill_buffer();

        } else {
            // This variant shouldn't go in this chunk.
        
            // Finish the chunk to a point before the next variant, before the
            // end of the reference, before the max chunk size, and after the
            // last variant the chunk contains.
            chunk_end = max(chunk_end,
                min((size_t ) variant_source.get()->position,
                    min((size_t) reference_end,
                        (size_t) (chunk_start + bases_per_chunk))));
            
            // Get the ref sequence we need
            auto chunk_ref = reference.getSubSequence(reference_contig, chunk_start, chunk_end - chunk_start);
            
            // Call the construction
            auto result = construct_chunk(chunk_ref, reference_contig, chunk_variants, chunk_start);
            
            // Wire up and emit the chunk graph
            wire_and_emit(result);
            
            // Set up a new chunk
            chunk_start = chunk_end;
            chunk_end = 0;
            chunk_variants.clear();
            
            // Loop again on the same variant.
        }
    }
    
    // We ran out of variants, so finish this chunk and all the others after it
    // without looking for variants.
    // TODO: unify with above loop?
    while (chunk_start < reference_end) {
        // We haven't finished the whole reference
    
        // Make the chunk as long as it can be
        chunk_end = max(chunk_end,
                    min((size_t) reference_end,
                        (size_t) (chunk_start + bases_per_chunk)));
    
        // Get the ref sequence we need
        auto chunk_ref = reference.getSubSequence(reference_contig, chunk_start, chunk_end - chunk_start);
        
        // Call the construction
        auto result = construct_chunk(chunk_ref, reference_contig, chunk_variants, chunk_start);
        
        // Wire up and emit the chunk graph
        wire_and_emit(result);
        
        // Set up a new chunk
        chunk_start = chunk_end;
        chunk_end = 0;
        chunk_variants.clear();
    }
    
    // All the chunks have been wired and emitted. Now emit the very last node, if any
    emit_reference_node(last_node_buffer);
    
}

void Constructor::construct_graph(const vector<FastaReference*>& references,
    const vector<vcflib::VariantCallFile*>& variant_files,
    function<void(Graph&)> callback) {

    // Make a map from contig name to fasta reference containing it.
    map<string, FastaReference*> reference_for;
    for (auto* reference : references) {
        // For every FASTA reference, make sure it has an index
        assert(reference->index);
        for (auto& kv : *(reference->index)) {
            // For every sequence name and index entry, point to this reference
            reference_for[kv.first] = reference;
        }
    }
    
    // Make VcfBuffers on all the variant files.
    vector<unique_ptr<VcfBuffer>> buffers;
    for (auto* vcf : variant_files) {
        // Every VCF gets a buffer wrapped around it.
        
        if (!vcf->is_open()) {
            // Except those that didn't open.
            continue;
        }
        
        // These will all get destructed when the vector goes away.
        buffers.emplace_back(new VcfBuffer(vcf));
    }
    
    if (!allowed_vcf_names.empty()) {
        // If we have a set of contigs to do, do those directly.
        
        for (string vcf_name : allowed_vcf_names) {
            // For each VCF contig, get the FASTA name
            string fasta_name = vcf_to_fasta(vcf_name);
            // Also the FASTA reference that has that sequence
            assert(reference_for.count(fasta_name));
            FastaReference* reference = reference_for[fasta_name];
        
            // We'll set this to true if we actually find the VCF that contains
            // the variants for this sequence.
            bool found_region = false;
        
            for (auto& buffer : buffers) {
                // For each VCF we are going to read
                if(!buffer->has_tabix()) {
                    // Die if we don't have indexes for everyone.
                    // TODO: report errors to caller instead.
                    cerr << "[vg::Constructor] Error: all VCFs must be indexed when restricting to a region" << endl;
                    exit(1);
                }
                
                // Try seeking to the right contig/region
                if (allowed_vcf_regions.count(vcf_name)) {
                    // Seek to just that region (0-based)
                    found_region = buffer->set_region(vcf_name, allowed_vcf_regions[vcf_name].first,
                        allowed_vcf_regions[vcf_name].second);
                } else {
                    // Seek to just the whole contig
                    found_region = buffer->set_region(vcf_name);
                }
                
                if (found_region) {
                    // This buffer is the one!
                    // Construct the graph for this contig with the FASTA and the VCF.
                    construct_graph(vcf_name, *reference, *buffer, callback);
                    break;
                }
            }
            
            if (!found_region) {
                // None of the VCFs include variants on this sequence.
                // Just build the graph for this sequence with no varaints.
                VcfBuffer empty(nullptr);
                construct_graph(vcf_name, *reference, empty, callback);
            }
        }
    } else {
        // If we have no set of contigs
        
        // Keep track of the contigs we have constructed, by VCF name
        set<string> constructed;
        
        for (auto& buffer : buffers) {
            // Go through all the VCFs
            // TODO: do this in parallel
            
            // Peek at the first variant and see its contig
            buffer->fill_buffer();
            while(buffer->get()) {
                // While there are still variants in the file
                // See what contig the next varianmt is on.
                string vcf_contig = buffer->get()->sequenceName;
                
                // Decide what FASTA contig that is and make sure we have it
                string fasta_contig = vcf_to_fasta(vcf_contig);
                assert(reference_for.count(fasta_contig));
                auto* reference = reference_for[fasta_contig];
                
                // Construct on it with the appropriate FastaReference for that contig
                construct_graph(vcf_contig, *reference, *buffer, callback);
                // Remember we did this one
                constructed.insert(vcf_contig);
                
                // After we're done constructing, scan until VCF EOF or a new contig comes up
                buffer->fill_buffer();
                while (buffer->get() && buffer->get()->sequenceName == vcf_contig) {
                    // Discard anything left on the same contig, since it must be
                    // out of our desired interval for that contig.
                    buffer->handle_buffer();
                    buffer->fill_buffer();
                }
            }
        }
        
        // Then for all the FASTA contigs that didn't appear in the VCFs,
        // construct them with no variants.
        
        for (auto& kv : reference_for) {
            // For every FASTA contig (and the reference that holds it)
            auto& fasta_contig = kv.first;
            FastaReference* reference = kv.second;
            
            // Convert the name to VCF space
            auto vcf_contig = fasta_to_vcf(fasta_contig);
            
            if (constructed.count(vcf_contig)) {
                // Skip contigs we already did in the VCF
                continue;
            }
            
            // Construct all the contigs we didn't do yet with no varaints.
            VcfBuffer empty(nullptr);
            construct_graph(vcf_contig, *reference, empty, callback);
        }
        
        // Now we've constructed everything we can. We're done!
        
        
    }
        
}

// Implementations of VG functions. TODO: refactor out of VG class

VG::VG(vcflib::VariantCallFile& variantCallFile,
       FastaReference& reference,
       string& target_region,
       bool target_is_chrom,
       int vars_per_region,
       int max_node_size,
       bool flat_input_vcf,
       bool load_phasing_paths,
       bool load_variant_alt_paths,
       bool showprog,
       set<string>* allowed_variants) {

    init();

    omp_set_dynamic(1); // use dynamic scheduling

    show_progress = showprog;

    map<string, VG*> refseq_graph;

    vector<string> targets;
    if (!target_region.empty()) {
        targets.push_back(target_region);
    } else {
        for (vector<string>::iterator r = reference.index->sequenceNames.begin();
             r != reference.index->sequenceNames.end(); ++r) {
            targets.push_back(*r);
        }
    }

    // How many phase paths do we want to load?
    size_t num_phasings = load_phasing_paths ? variantCallFile.sampleNames.size() * 2 : 0;
    // We'll later split these where you would have to take an edge that doesn't exist.

    // to scale up, we have to avoid big string memcpys
    // this could be accomplished by some deep surgery on the construction routines
    // however, that could be a silly thing to do,
    // because why break something that's conceptually clear
    // and anyway, we want to break the works into chunks
    //
    // there is a development that could be important
    // our chunk size isn't going to reach into the range where we'll have issues (>several megs)
    // so we'll run this for regions of moderate size, scaling up in the case that we run into a big deletion
    //

    for (vector<string>::iterator t = targets.begin(); t != targets.end(); ++t) {

        //string& seq_name = *t;
        string seq_name;
        string target = *t;
        int start_pos = 0, stop_pos = 0;
        // nasty hack for handling single regions
        if (!target_is_chrom) {
            parse_region(target,
                         seq_name,
                         start_pos,
                         stop_pos);
            if (stop_pos > 0) {
                if (variantCallFile.is_open()) {
                    variantCallFile.setRegion(seq_name, start_pos, stop_pos);
                }
            } else {
                if (variantCallFile.is_open()) {
                    variantCallFile.setRegion(seq_name);
                }
                stop_pos = reference.sequenceLength(seq_name);
            }
        } else {
            // the user said the target is just a sequence name
            // and is unsafe to parse as it may contain ':' or '-'
            // for example "gi|568815592:29791752-29792749"
            if (variantCallFile.is_open()) {
                variantCallFile.setRegion(target);
            }
            stop_pos = reference.sequenceLength(target);
            seq_name = target;
        }
        vcflib::Variant var(variantCallFile);

        vector<vcflib::Variant>* region = NULL;

        // convert from 1-based input to 0-based internal format
        // and handle the case where we are already doing the whole chromosome
        id_t start = start_pos ? start_pos - 1 : 0;
        id_t end = start;

        create_progress("loading variants for " + target, stop_pos-start_pos);
        // get records
        vector<vcflib::Variant> records;

        // This is going to hold the alleles that occur at certain reference
        // positions, in addition to the reference allele. We keep them ordered
        // so we can refer to them by number.
        map<long,vector<vcflib::VariantAllele> > alleles;

        // This is going to hold, for each position, allele combination, a
        // vector of bools marking which phases of which samples visit that
        // allele. Each sample is stored at (sample number * 2) for phase 0 and
        // (sample number * 2 + 1) for phase 1. The reference may not always get
        // an allele, but if anything is reference it will show up as an
        // overlapping allele elsewhere.
        map<pair<long, int>, vector<bool>> phase_visits;

        // This is going to hold visits to VariantAlleles by the reference and
        // nonreference alts of variants. We map from VariantAllele index and
        // number to a list of the variant ID and alt number pairs that use the
        // VariantAllele.
        map<pair<long, int>, vector<pair<string, int>>> variant_alts;

        // We don't want to load all the vcf records into memory at once, since
        // the vcflib internal data structures are big compared to the info we
        // need.
        int64_t variant_chunk_size = 1000;

        auto parse_loaded_variants = [&]() {
            // Parse the variants we have loaded, and clear them out, so we can
            // go back and load a new batch of variants.

            // decompose records into alleles with offsets against our target
            // sequence Dump the collections of alleles (which are ref, alt
            // pairs) into the alleles map. Populate the phase visit map if
            // we're loading phasing paths, and the variant alt path map if
            // we're loading variant alts.
            vcf_records_to_alleles(records, alleles,
                load_phasing_paths ? &phase_visits : nullptr,
                load_variant_alt_paths ? &variant_alts : nullptr,
                flat_input_vcf);
            records.clear(); // clean up
        };

        int64_t i = 0;
        while (variantCallFile.is_open() && variantCallFile.getNextVariant(var)) {
            // this ... maybe we should remove it as for when we have calls against N
            bool isDNA = allATGC(var.ref);
            for (vector<string>::iterator a = var.alt.begin(); a != var.alt.end(); ++a) {
                if (!allATGC(*a)) isDNA = false;
            }
            // only work with DNA sequences
            if (isDNA) {
                string vrepr = var.vrepr();
                var.position -= 1; // convert to 0-based
                if (allowed_variants == nullptr
                    || allowed_variants->count(vrepr)) {
                    records.push_back(var);
                }
            }
            if (++i % 1000 == 0) update_progress(var.position-start_pos);
            // Periodically parse the records down to what we need and throw away the rest.
            if (i % variant_chunk_size == 0) parse_loaded_variants();
        }
        // Finish up any remaining unparsed variants
        parse_loaded_variants();

        destroy_progress();

        // store our construction plans
        deque<Plan*> construction;
        // so we can check which graphs we can safely append
        set<VG*> graph_completed;
        // we add and remove from graph_completed, so track count for logging
        int graphs_completed = 0;
        int final_completed = -1; // hm
        // the construction queue
        list<VG*> graphq;
        int graphq_size = 0; // for efficiency
        // ^^^^ (we need to insert/remove things in the middle of the list,
        // but we also need to be able to quickly determine its size)
        // for tracking progress through the chromosome
        map<VG*, unsigned long> graph_end;

        create_progress("planning construction", stop_pos-start_pos);
        // break into chunks
        int chunk_start = start;
        bool invariant_graph = alleles.empty();
        while (invariant_graph || !alleles.empty()) {
            invariant_graph = false;
            map<long, vector<vcflib::VariantAllele> > new_alleles;
            map<pair<long, int>, vector<bool>> new_phase_visits;
            map<pair<long, int>, vector<pair<string, int>>> new_variant_alts;
            // our start position is the "offset" we should subtract from the
            // alleles and the phase visits for correct construction
            //chunk_start = (!chunk_start ? 0 : alleles.begin()->first);
            int chunk_end = chunk_start;
            bool clean_end = true;
            for (int i = 0; (i < vars_per_region || !clean_end) && !alleles.empty(); ++i) {
                auto pos = alleles.begin()->first - chunk_start;
                chunk_end = max(chunk_end, (int)alleles.begin()->first);
                auto& pos_alleles = alleles.begin()->second;
                // apply offset when adding to the new alleles
                auto& curr_pos = new_alleles[pos];
                for (int j = 0; j < pos_alleles.size(); j++) {
                    // Go through every allele that occurs at this position, and
                    // update it to the offset position in new_alleles
                    auto& allele = pos_alleles[j];

                    // We'll clone and modify it.
                    auto new_allele = allele;
                    int ref_end = new_allele.ref.size() + new_allele.position;
                    // look through the alleles to see if there is a longer chunk
                    if (ref_end > chunk_end) {
                        chunk_end = ref_end;
                    }
                    new_allele.position = pos;
                    // Copy the modified allele over.
                    // No need to deduplicate.
                    curr_pos.push_back(new_allele);

                    // Also handle any visits to this allele
                    // We need the key, consisting of the old position and the allele number there.
                    auto old_allele_key = make_pair(alleles.begin()->first, j);
                    // Make the new key
                    auto new_allele_key = make_pair(pos, j);
                    if(phase_visits.count(old_allele_key)) {
                        // We have some usages of this allele for phase paths. We need to move them over.

                        // Move over the value and insert into the new map. See <http://stackoverflow.com/a/14816487/402891>
                        // TODO: would it be clearer with the braces instead?
                        new_phase_visits.insert(make_pair(new_allele_key, std::move(phase_visits.at(old_allele_key))));

                        // Now we've emptied out/made-undefined the old vector,
                        // so we probably should drop it from the old map.
                        phase_visits.erase(old_allele_key);
                    }

                    if(variant_alts.count(old_allele_key)) {
                        // We have some usages of this allele by variant alts. We need to move them over.

                        // Do a move operation
                        new_variant_alts.insert(make_pair(new_allele_key, std::move(variant_alts.at(old_allele_key))));
                        // Delete the olkd entry (just so we don't keep it around wasting time/space/being unspecified)
                        variant_alts.erase(old_allele_key);
                    }
                }
                alleles.erase(alleles.begin());
                // TODO here we need to see if we are neighboring another variant
                // and if we are, keep constructing
                if (alleles.begin()->first <= chunk_end) {
                    clean_end = false;
                } else {
                    clean_end = true;
                }
            }
            // record end position, use target end in the case that we are at the end
            if (alleles.empty()) chunk_end = stop_pos;

            // we set the head graph to be this one, so we aren't obligated to copy the result into this object
            // make a construction plan
            Plan* plan = new Plan(graphq.empty() && targets.size() == 1 ? this : new VG,
                                  std::move(new_alleles),
                                  std::move(new_phase_visits),
                                  std::move(new_variant_alts),
                                  reference.getSubSequence(seq_name,
                                                           chunk_start,
                                                           chunk_end - chunk_start),
                                  seq_name);
            chunk_start = chunk_end;
#pragma omp critical (graphq)
            {
                graphq.push_back(plan->graph);
                construction.push_back(plan);
                if (show_progress) graph_end[plan->graph] = chunk_end;
                update_progress(chunk_end);
            }
        }
#ifdef debug
        cerr << omp_get_thread_num() << ": graphq size " << graphq.size() << endl;
#endif
        graphq_size = graphq.size();
        destroy_progress();

        // this system is not entirely general
        // there will be a problem when the regions of overlapping deletions become too large
        // then the inter-dependence of each region will make parallel construction in this way difficult
        // because the chunks will get too large

        // use this function to merge graphs both during and after the construction iteration
        auto merge_first_two_completed_graphs =
            [this, start_pos, &graph_completed, &graphq, &graphq_size, &graph_end, &final_completed](void) {
            // find the first two consecutive graphs which are completed
            VG* first = NULL;
            VG* second = NULL;
//#pragma omp critical (cerr)
//            cerr << omp_get_thread_num() << ": merging" << endl;
#pragma omp critical (graphq)
            {
                auto itp = graphq.begin(); // previous
                auto itn = itp; if (itp != graphq.end()) ++itn; // next
                // scan the graphq to find consecutive entries that are both completed
                while (itp != itn // there is > 1 entry
                       && itn != graphq.end() // we aren't yet at the end
                       && !(graph_completed.count(*itp) // the two we're looking at aren't completed
                            && graph_completed.count(*itn))) {
                    ++itp; ++itn;
                }

                if (itn != graphq.end()) {
                    // we have two consecutive graphs to merge!
                    first = *itp;
                    second = *itn;
                    // unset graph completed for both
                    graph_completed.erase(first);
                    graph_completed.erase(second);
                    graphq.erase(itn);
                    --graphq_size;
                }
            }

            if (first && second) {
                // combine graphs
                first->append(*second);
#pragma omp critical (graphq)
                {
                    if (final_completed != -1) update_progress(final_completed++);
                    graph_completed.insert(first);
                    graph_end.erase(second);
                }
                delete second;
            }
        };

        create_progress("constructing graph", construction.size());

        // (in parallel) construct each component of the graph
#pragma omp parallel for
        for (int i = 0; i < construction.size(); ++i) {

            int tid = omp_get_thread_num();
            Plan* plan = construction.at(i);
#ifdef debug
#pragma omp critical (cerr)
            cerr << tid << ": " << "constructing graph " << plan->graph << " over "
                 << plan->alleles.size() << " variants in " <<plan->seq.size() << "bp "
                 << plan->name << endl;
#endif

            // Make the piece of graph, passing along the number of sample phases if we're making phase paths.
            plan->graph->from_alleles(plan->alleles,
                                      plan->phase_visits,
                                      num_phasings,
                                      plan->variant_alts,
                                      plan->seq,
                                      plan->name);

            // Break up the nodes ourselves
            if(max_node_size > 0) {
                plan->graph->dice_nodes(max_node_size);
            }

#pragma omp critical (graphq)
            {
                update_progress(++graphs_completed);
                graph_completed.insert(plan->graph);
#ifdef debug
#pragma omp critical (cerr)
                cerr << tid << ": " << "constructed graph " << plan->graph << endl;
#endif
            }
            // clean up
            delete plan;

            // concatenate chunks of the result graph together
            merge_first_two_completed_graphs();

        }
        destroy_progress();

        // merge remaining graphs
        final_completed = 0;
        create_progress("merging remaining graphs", graphq.size());
#pragma omp parallel
        {
            bool more_to_merge = true;
            while (more_to_merge) {
                merge_first_two_completed_graphs();
                usleep(10);
#pragma omp critical (graphq)
                more_to_merge = graphq_size > 1;
            }
        }
        destroy_progress();

        // parallel end
        // finalize target

        // our target graph should be the only entry in the graphq
        assert(graphq.size() == 1);
        VG* target_graph = graphq.front();

        // store it in our results
        refseq_graph[target] = target_graph;

        create_progress("joining graphs", target_graph->size());
        // clean up "null" nodes that are used for maintaining structure between temporary subgraphs
        target_graph->remove_null_nodes_forwarding_edges();
        destroy_progress();

        // then use topological sorting and re-compression of the id space to make sure that
        create_progress("topologically sorting", target_graph->size());
        target_graph->sort();
        destroy_progress();

        create_progress("compacting ids", target_graph->size());
        // we get identical graphs no matter what the region size is
        target_graph->compact_ids();
        destroy_progress();

    }

    // hack for efficiency when constructing over a single chromosome
    if (refseq_graph.size() == 1) {
        // *this = *refseq_graph[targets.front()];
        // we have already done this because the first graph in the queue is this
    } else {
        // where we have multiple targets
        for (vector<string>::iterator t = targets.begin(); t != targets.end(); ++t) {
            // merge the variants into one graph
            VG& g = *refseq_graph[*t];
            combine(g);
        }
    }
    // rebuild the mapping ranks now that we've combined everything
    paths.clear_mapping_ranks();
    paths.rebuild_mapping_aux();

    if(load_phasing_paths) {
        // Trace through all the phase paths, and, where they take edges that
        // don't exist, break them. TODO: we still might get spurious phasing
        // through a deletion where the two pahsed bits but up against each
        // other.

        create_progress("dividing phasing paths", num_phasings);
        for(size_t i = 0; i < num_phasings; i++) {
            // What's the path we want to trace?
            string original_path_name = "_phase" + to_string(i);

            list<Mapping>& path_mappings = paths.get_path(original_path_name);

            // What section of this phasing do we want to be outputting?
            size_t subpath = 0;
            // Make a name for it
            string subpath_name = "_phase" + to_string(i) + "_" + to_string(subpath);

            // For each mapping, we want to be able to look at the previous
            // mapping.
            list<Mapping>::iterator prev_mapping = path_mappings.end();
            for(list<Mapping>::iterator mapping = path_mappings.begin(); mapping != path_mappings.end(); ++mapping) {
                // For each mapping in the path
                if(prev_mapping != path_mappings.end()) {
                    // We have the previous mapping and this one

                    // Make the two sides of nodes that should be connected.
                    auto s1 = NodeSide(prev_mapping->position().node_id(),
                        (prev_mapping->position().is_reverse() ? false : true));
                    auto s2 = NodeSide(mapping->position().node_id(),
                        (mapping->position().is_reverse() ? true : false));
                    // check that we always have an edge between the two nodes in the correct direction
                    if (!has_edge(s1, s2)) {
                        // We need to split onto a new subpath;
                        subpath++;
                        subpath_name = "_phase" + to_string(i) + "_" + to_string(subpath);
                    }
                }

                // Now we just drop this node onto the current subpath
                paths.append_mapping(subpath_name, *mapping);

                // Save this mapping as the prev one
                prev_mapping = mapping;
            }

            // Now delete the original full phase path.
            // This invalidates the path_mappings reference!!!
            // We use the variant that actually unthreads the path from the indexes and doesn't erase and rebuild them.
            paths.remove_path(original_path_name);

            update_progress(i);
        }
        destroy_progress();


    }

    std::function<bool(string)> all_upper = [](string s){
        //GO until [size() - 1 ] to avoid the newline char
        for (int i = 0; i < s.size() - 1; i++){
            if (!isupper(s[i])){
                return false;
            }
        }
        return true;
    };

    for_each_node([&](Node* node) {
            if (!all_upper(node->sequence())){
                cerr << "WARNING: Lower case letters found during construction" << endl;
                cerr << "Sequences may not map to this graph." << endl;
                cerr << pb2json(*node) << endl;
            }
        });

}

// construct from VCF records
// --------------------------
// algorithm
// maintain a core reference path upon which we add new variants as they come
// addition procedure is the following
// find reference node overlapping our start position
// if it is already the end of a node, add the new node
// if it is not the end of a node, break it, insert edges from old->new
// go to end position of alt allele (could be the same position)
// if it already has a break, just point to the next node in line
// if it is not broken, break it and point to the next node
// add new node for alt alleles, connect to start and end node in reference path
// store the ref mapping as a property of the edges and nodes (this allows deletion edges and insertion subpaths)
//

void VG::vcf_records_to_alleles(vector<vcflib::Variant>& records,
                                map<long, vector<vcflib::VariantAllele> >& altp,
                                map<pair<long, int>, vector<bool>>* phase_visits,
                                map<pair<long, int>, vector<pair<string, int>>>* alt_allele_visits,
                                bool flat_input_vcf) {



#ifdef debug
    cerr << "Processing " << records.size() << " vcf records..." << endl;
#endif

    for (int i = 0; i < records.size(); ++i) {
        vcflib::Variant& var = records.at(i);

        // What name should we use for the variant? We need to make sure it is
        // unique, even if there are multiple variant records at the same
        // position in the VCF. Also, we don't necessarily have every variant in
        // the VCF in our records vector.
        string var_name = make_variant_id(var);

        // decompose to alts
        // This holds a map from alt or ref allele sequence to a series of VariantAlleles describing an alignment.
        map<string, vector<vcflib::VariantAllele> > alternates
            = (flat_input_vcf ? var.flatAlternates() : var.parsedAlternates());

        if(!alternates.count(var.ref)) {
            // Ref is missing, as can happen with flat construction.
            // Stick the ref in, because we need to have ref.
            alternates[var.ref].push_back(vcflib::VariantAllele(var.ref, var.ref, var.position));
        }

        // This holds a map from alt index (0 for ref) to the phase sets
        // visiting it as a bool vector. No bit vector means no visits.
        map<int, vector<bool>> alt_usages;

        if(phase_visits != nullptr) {

            // Parse out what alleles each sample uses in its phase sets at this
            // VCF record.

            // Get all the sample names in order.
            auto& sample_names = var.vcf->sampleNames;

            for(int64_t j = 0; j < sample_names.size(); j++) {
                // For every sample, see if at this variant it uses this
                // allele in one or both phase sets.

                // Grab the genotypes
                string genotype = var.getGenotype(sample_names[j]);

                // Find the phasing bar
                auto bar_pos = genotype.find('|');

                if(bar_pos == string::npos || bar_pos == 0 || bar_pos + 1 >= genotype.size()) {
                    // Not phased here, or otherwise invalid
                    continue;
                }

                if(genotype.substr(0, bar_pos) == "." || genotype.substr(bar_pos + 1) == ".") {
                    // This site is uncalled
                    continue;
                }

                // Parse out the two alt indexes.
                // TODO: complain if there are more.
                int alt1index = stoi(genotype.substr(0, bar_pos));
                int alt2index = stoi(genotype.substr(bar_pos + 1));

                if(!alt_usages.count(alt1index)) {
                    // Make a new bit vector for the alt visited by 1
                    alt_usages[alt1index] = vector<bool>(var.getNumSamples() * 2, false);
                }
                // First phase of this phase set visits here.
                alt_usages[alt1index][j * 2] = true;

                if(!alt_usages.count(alt2index)) {
                    // Make a new bit vector for the alt visited by 2
                    alt_usages[alt2index] = vector<bool>(var.getNumSamples() * 2, false);
                }
                // Second phase of this phase set visits here.
                alt_usages[alt2index][j * 2 + 1] = true;
            }
        }

        for (auto& alleles : alternates) {

            // We'll point this to a vector flagging all the phase visits to
            // this alt (which may be the ref alt), if we want to record those.
            vector<bool>* visits = nullptr;

            // What alt number is this alt? (0 for ref)
            // -1 for nothing needs to visit it and we don't care.
            int alt_number = -1;

#ifdef debug
            cerr << "Considering alt " << alleles.first << " at " << var.position << endl;
            cerr << var << endl;
#endif

            if(phase_visits != nullptr || alt_allele_visits != nullptr) {
                // We actually have visits to look for. We need to know what
                // alt number we have here.

                // We need to copy out the alt sequence to appease the vcflib API
                string alt_sequence = alleles.first;

                // What alt number are we looking at
                if(alt_sequence == var.ref) {
                    // This is the ref allele
                    alt_number = 0;
                } else {
                    // This is an alternate allele
                    alt_number = var.getAltAlleleIndex(alt_sequence) + 1;
                }

#ifdef debug
                cerr << "Alt is number " << alt_number << endl;
#endif

                if(alt_usages.count(alt_number)) {
                    // Something did indeed visit. Point the pointer at the
                    // vector describing what visited.
                    visits = &alt_usages[alt_number];
                }
            }

            for (auto& allele : alleles.second) {
                // For each of the alignment bubbles or matches, add it in as something we'll need for the graph.
                // These may overlap between alleles, and not every allele will have one at all positions.
                // In general it has to be that way, because the alleles themselves can overlap.

                // TODO: we need these to be unique but also ordered by addition
                // order. For now we just check all previous entries before
                // adding and suffer being n^2 in vcf alts per variant. We
                // should use some kind of addition-ordered set.
                int found_at = -1;
                for(int j = 0; j < altp[allele.position].size(); j++) {
                    if(altp[allele.position][j].ref == allele.ref && altp[allele.position][j].alt == allele.alt) {
                        // TODO: no equality for VariantAlleles for some reason.
                        // We already have it at this index
                        found_at = j;
                        break;
                    }
                }
                if(found_at == -1) {
                    // We need to tack this on at the end.
                    found_at = altp[allele.position].size();
                    // Add the bubble made by this part of this alt at this
                    // position.
                    altp[allele.position].push_back(allele);
                }

                if(visits != nullptr && phase_visits != nullptr) {
                    // We have to record a phase visit

                    // What position, allele index pair are we visiting when we
                    // visit this alt?
                    auto visited = make_pair(allele.position, found_at);

                    if(!phase_visits->count(visited)) {
                        // Make sure we have a vector for visits to this allele, not
                        // just this alt. It needs an entry for each phase of each sample.
                        (*phase_visits)[visited] = vector<bool>(var.getNumSamples() * 2, false);
                    }

                    for(size_t j = 0; j < visits->size(); j++) {
                        // We need to toggle on all the phase sets that visited
                        // this alt as using this allele at this position.
                        if(visits->at(j) && !(*phase_visits)[visited].at(j)) {
                            // The bit needs to be set, because all the phases
                            // visiting this alt visit this allele that appears
                            // in it.
                            (*phase_visits)[visited][j] = true;
                        }

                    }
                }

                if(alt_allele_visits != nullptr && alt_number != -1) {
                    // We have to record a visit of this alt of this variant to
                    // this VariantAllele bubble/reference patch.

                    // What position, allele index pair are we visiting when we
                    // visit this alt?
                    auto visited = make_pair(allele.position, found_at);

#ifdef debug
                    cerr << var_name << " alt " << alt_number << " visits allele #" << found_at
                        << " at position " << allele.position << " of " << allele.ref << " -> " << allele.alt << endl;
#endif

                    // Say we visit this allele as part of this alt of this variant.
                    (*alt_allele_visits)[visited].push_back(make_pair(var_name, alt_number));
                }

            }
        }
    }
}

void VG::from_alleles(const map<long, vector<vcflib::VariantAllele> >& altp,
                      const map<pair<long, int>, vector<bool>>& visits,
                      size_t num_phasings,
                      const map<pair<long, int>, vector<pair<string, int>>>& variant_alts,
                      string& seq,
                      string& name) {

    //init();
    this->name = name;

    int tid = omp_get_thread_num();

#ifdef debug
#pragma omp critical (cerr)
    {
        cerr << tid << ": in from_alleles" << endl;
        cerr << tid << ": with " << altp.size() << " vars" << endl;
        cerr << tid << ": and " << num_phasings << " phasings" << endl;
        cerr << tid << ": and " << visits.size() << " phasing visits" << endl;
        cerr << tid << ": and " << variant_alts.size() << " variant alt visits" << endl;
        cerr << tid << ": and " << seq.size() << "bp" << endl;
        if(seq.size() < 100) cerr << seq << endl;
    }
#endif


    // maintains the path of the seq in the graph
    map<long, id_t> seq_node_ids;
    // track the last nodes so that we can connect everything
    // completely when variants occur in succession
    map<long, set<Node*> > nodes_by_end_position;
    map<long, set<Node*> > nodes_by_start_position;


    Node* seq_node = create_node(seq);
    // This path represents the primary path in this region of the graph. We
    // store it as a map for now, and add it in in the real Paths structure
    // later.
    seq_node_ids[0] = seq_node->id();

    // TODO: dice nodes now so we can work only with small ref nodes?
    // But what if we then had a divided middle node?

    // We can't reasonably track visits to the "previous" bunch of alleles
    // because they may really overlap this bunch of alleles and not be properly
    // previous, path-wise. We'll just assume all the phasings visit all the
    // non-variable nodes, and then break things up later. TODO: won't this
    // artificially merge paths if we have an unphased deletion or something?

    // Where did the last variant end? If it's right before this one starts,
    // there might not be an intervening node.
    long last_variant_end = -1;

    for (auto& va : altp) {

        const vector<vcflib::VariantAllele>& alleles = va.second;

        // if alleles are empty, we just cut at this point. TODO: this should
        // never happen with the node size enforcement refactoring.
        if (alleles.empty()) {
            Node* l = NULL; Node* r = NULL;
            divide_path(seq_node_ids, va.first, l, r);
        }


        // If all the alleles here are perfect reference matches, and no
        // variants visit them, we'll have nothing to do.
        bool all_perfect_matches = true;
        for(auto& allele : alleles) {
            if(allele.ref != allele.alt) {
                all_perfect_matches = false;
                break;
            }
        }

        // Are all the alleles here clear of visits by variants?
        bool no_variant_visits = true;

        for (size_t allele_number = 0; allele_number < alleles.size(); allele_number++) {
            if(variant_alts.count(make_pair(va.first, allele_number))) {
                no_variant_visits = false;
                break;
            }
        }

        if(all_perfect_matches && no_variant_visits) {
            // No need to break anything here.

#ifdef debug
#pragma omp critical (cerr)
            {
                cerr << tid << ": Skipping entire allele site at " << va.first << endl;
            }
#endif

            continue;
        }

        // We also need to sort the allele numbers by the lengths of their
        // alleles' reference sequences, to properly handle inserts followed by
        // matches.
        vector<int> allele_numbers_by_ref_length(alleles.size());
        // Fill with sequentially increasing integers.
        // Sometimes the STL actually *does* have the function you want.
        iota(allele_numbers_by_ref_length.begin(), allele_numbers_by_ref_length.end(), 0);

        // Sort the allele numbers by reference length, ascending
        std::sort(allele_numbers_by_ref_length.begin(), allele_numbers_by_ref_length.end(),
            [&](const int& a, const int& b) -> bool {
            // Sort alleles with shorter ref sequences first.
            return alleles[a].ref.size() < alleles[b].ref.size();
        });

#ifdef debug
#pragma omp critical (cerr)
                {
                    cerr << tid << ": Processing " << allele_numbers_by_ref_length.size() << " alleles at " << va.first << endl;
                }
#endif

        // Is this allele the first one processed? Because the first one
        // processed gets to handle adding mappings to the intervening sequence
        // from the previous allele to here.
        bool first_allele_processed = true;

        for (size_t allele_number : allele_numbers_by_ref_length) {
            // Go through all the alleles with their numbers, in order of
            // increasing reference sequence length (so inserts come first)
            auto& allele = alleles[allele_number];

            auto allele_key = make_pair(va.first, allele_number);

            // 0/1 based conversion happens in offset
            long allele_start_pos = allele.position;
            long allele_end_pos = allele_start_pos + allele.ref.size();
            // for ordering, set insertion start position at +1
            // otherwise insertions at the same position will loop infinitely
            //if (allele_start_pos == allele_end_pos) allele_end_pos++;

            if(allele.ref == allele.alt && !visits.count(allele_key) && !variant_alts.count(allele_key)) {
                // This is a ref-only allele with no visits or usages in
                // alleles, which means we don't actually need any cuts if the
                // allele is not visited. If other alleles here are visited,
                // we'll get cuts from them.

#ifdef debug
#pragma omp critical (cerr)
                {
                    cerr << tid << ": Skipping variant at " << allele_start_pos
                         << " allele " << allele.ref << " -> " << allele.alt << endl;
                }
#endif

                continue;
            }

#ifdef debug
#pragma omp critical (cerr)
            {
                cerr << tid << ": Handling variant at " << allele_start_pos
                     << " allele " << allele.ref << " -> " << allele.alt << endl;
            }
#endif

            if (allele_start_pos == 0) {
                // ensures that we can handle variation at first position
                // (important when aligning)
                Node* root = create_node("");
                seq_node_ids[-1] = root->id();
                nodes_by_start_position[-1].insert(root);
                nodes_by_end_position[0].insert(root);
            }



            // We grab all the nodes involved in this allele: before, being
            // replaced, and after.
            Node* left_seq_node = nullptr;
            std::list<Node*> middle_seq_nodes;
            Node* right_seq_node = nullptr;

            // make one cut at the ref-path relative start of the allele, if it
            // hasn't been cut there already. Grab the nodes on either side of
            // that cut.
            divide_path(seq_node_ids,
                        allele_start_pos,
                        left_seq_node,
                        right_seq_node);

            // if the ref portion of the allele is not empty, then we may need
            // to make another cut. If so, we'll have some middle nodes.
            if (!allele.ref.empty()) {
                Node* last_middle_node = nullptr;
                divide_path(seq_node_ids,
                            allele_end_pos,
                            last_middle_node,
                            right_seq_node);


                // Now find all the middle nodes between left_seq_node and
                // last_middle_node along the primary path.

                // Find the node starting at or before, and including,
                // allele_end_pos.
                map<long, id_t>::iterator target = seq_node_ids.upper_bound(allele_end_pos);
                --target;

                // That should be the node to the right of the variant
                assert(target->second == right_seq_node->id());

                // Everything left of there, stopping (exclusive) at
                // left_seq_node if set, should be middle nodes.

                while(target != seq_node_ids.begin()) {
                    // Don't use the first node we start with, and do use the
                    // begin node.
                    target--;
                    if(left_seq_node != nullptr && target->second == left_seq_node->id()) {
                        // Don't put the left node in as a middle node
                        break;
                    }

                    // If we get here we want to take this node as a middle node
                    middle_seq_nodes.push_front(get_node(target->second));
                }

                // There need to be some nodes in the list when we're done.
                // Otherwise something has gone wrong.
                assert(middle_seq_nodes.size() > 0);
            }

            // What nodes actually represent the alt allele?
            std::list<Node*> alt_nodes;
            // create a new alt node and connect the pieces from before
            if (!allele.alt.empty() && !allele.ref.empty()) {
                //cerr << "both alt and ref have sequence" << endl;

                if (allele.ref == allele.alt) {
                    // We don't really need to make a new run of nodes, just use
                    // the existing one. We still needed to cut here, though,
                    // because we can't have only a ref-matching allele at a
                    // place with alleles; there must be some other different
                    // alleles here.
                    alt_nodes = middle_seq_nodes;
                } else {
                    // We need a new node for this sequence
                    Node* alt_node = create_node(allele.alt);
                    create_edge(left_seq_node, alt_node);
                    create_edge(alt_node, right_seq_node);
                    alt_nodes.push_back(alt_node);
                }

                // The ref and alt nodes may be the same, but neither will be an
                // empty list.
                nodes_by_end_position[allele_end_pos].insert(alt_nodes.back());
                nodes_by_end_position[allele_end_pos].insert(middle_seq_nodes.back());
                //nodes_by_end_position[allele_start_pos].insert(left_seq_node);
                nodes_by_start_position[allele_start_pos].insert(alt_nodes.front());
                nodes_by_start_position[allele_start_pos].insert(middle_seq_nodes.front());

            } else if (!allele.alt.empty()) { // insertion

                // Make a single node to represent the inserted sequence
                Node* alt_node = create_node(allele.alt);
                create_edge(left_seq_node, alt_node);
                create_edge(alt_node, right_seq_node);
                alt_nodes.push_back(alt_node);

                // We know the alt nodes list isn't empty.
                // We'rr immediately pulling the node out of the list again for consistency.
                nodes_by_end_position[allele_end_pos].insert(alt_nodes.back());
                nodes_by_end_position[allele_end_pos].insert(left_seq_node);
                nodes_by_start_position[allele_start_pos].insert(alt_nodes.front());

            } else {// otherwise, we have a deletion, or the empty reference alt of an insertion.

                // No alt nodes should be present
                create_edge(left_seq_node, right_seq_node);
                nodes_by_end_position[allele_end_pos].insert(left_seq_node);
                nodes_by_start_position[allele_start_pos].insert(left_seq_node);

            }

#ifdef debug
#pragma omp critical (cerr)
            {
                if (left_seq_node) cerr << tid << ": left_ref " << left_seq_node->id()
                                        << " "
                                        << (left_seq_node->sequence().size() < 100 ? left_seq_node->sequence() : "...")
                                        << endl;
                for(Node* middle_seq_node : middle_seq_nodes) {
                    cerr << tid << ": middle_ref " << middle_seq_node->id()
                         << " " << middle_seq_node->sequence() << endl;
                }
                for(Node* alt_node : alt_nodes) {
                    cerr << tid << ": alt_node " << alt_node->id()
                                << " " << alt_node->sequence() << endl;
                }
                if (right_seq_node) cerr << tid << ": right_ref " << right_seq_node->id()
                                         << " "
                                         << (right_seq_node->sequence().size() < 100 ? right_seq_node->sequence() : "...")
                                         << endl;
            }
#endif

            // How much intervening space is there between this set of alleles'
            // start and the last one's end?
            long intervening_space = allele.position - last_variant_end;
            if(first_allele_processed && num_phasings > 0 && left_seq_node && intervening_space > 0) {
                // On the first pass through, if we are doing phasings, we make
                // all of them visit the left node. We know the left node will
                // be the same on subsequent passes for other alleles starting
                // here, and we only want to make these left node visits once.

                // However, we can only do this if there actually is a node
                // between the last set of alleles and here.

                // TODO: what if some of these phasings aren't actually phased
                // here? We'll need to break up their paths to just have some
                // ref matching paths between variants where they aren't
                // phased...

                for(size_t i = 0; i < num_phasings; i++) {
                    // Everything uses this node to our left, which won't be
                    // broken again.
                    paths.append_mapping("_phase" + to_string(i), left_seq_node->id());
                }

                // The next allele won't be the first one actually processed.
                first_allele_processed = false;
            }
            if(!alt_nodes.empty() && visits.count(allele_key)) {
                // At least one phased path visits this allele, and we have some
                // nodes to path it through.

                // Get the vector of bools for that phasings visit
                auto& visit_vector = visits.at(allele_key);

                for(size_t i = 0; i < visit_vector.size(); i++) {
                    // For each phasing
                    if(visit_vector[i]) {
                        // If we visited this allele, say we did. TODO: use a
                        // nice rank/select thing here to make this not have to
                        // be a huge loop.

                        string phase_name = "_phase" + to_string(i);

                        for(Node* alt_node : alt_nodes) {
                            // Problem: we may have visited other alleles that also used some of these nodes.
                            // Solution: only add on the mappings for new nodes.

                            // TODO: this assumes we'll not encounter
                            // contradictory alleles, only things like "both
                            // shorter ref match and longer ref match are
                            // visited".

                            if(!paths.get_node_mapping(alt_node).count(phase_name)) {
                                // This node has not yet been visited on this path.
                                paths.append_mapping(phase_name, alt_node->id());
                            }
                        }
                    }
                }
            }

            if(variant_alts.count(allele_key)) {

                for(auto name_and_alt : variant_alts.at(allele_key)) {
                    // For each of the alts using this allele, put mappings for this path
                    string path_name = "_alt_" + name_and_alt.first + "_" + to_string(name_and_alt.second);

                    if(!alt_nodes.empty()) {
                        // This allele has some physical presence and is used by some
                        // variants.

                        for(auto alt_node : alt_nodes) {
                            // Put a mapping on each alt node

                            // TODO: assert that there's an edge from the
                            // previous mapping's node (if any) to this one's
                            // node.

                            paths.append_mapping(path_name, alt_node->id());
                        }

#ifdef debug
                        cerr << "Path " << path_name << " uses these alts" << endl;
#endif

                    } else {
                        // TODO: alts that are deletions don't always have nodes
                        // on both sides to visit. Either anchor your VCF
                        // deletions at both ends, or rely on the presence of
                        // mappings to other alleles (allele 0) in this variant
                        // but not this allele to indicate the deletion of
                        // nodes.

#ifdef debug
                        cerr << "Path " << path_name << " would use these alts if there were any" << endl;
#endif
                    }
                }
            }

            if (allele_end_pos == seq.size()) {
                // ensures that we can handle variation at last position (important when aligning)
                Node* end = create_node("");
                seq_node_ids[allele_end_pos] = end->id();
                // for consistency, this should be handled below in the start/end connections
                if (alt_nodes.size() > 0) {
                    create_edge(alt_nodes.back(), end);
                }
                if (middle_seq_nodes.size() > 0) {
                    create_edge(middle_seq_nodes.back(), end);
                }
            }

            //print_edges();
            /*
            if (!is_valid()) {
                cerr << "graph is invalid after variant " << *a << endl;
                std::ofstream out("fail.vg");
                serialize_to_ostream(out);
                out.close();
                exit(1);
            }
            */

        }

        // Now we need to connect up all the extra deges between variant alleles
        // that abut each other.
        map<long, set<Node*> >::iterator ep
            = nodes_by_end_position.find(va.first);
        map<long, set<Node*> >::iterator sp
            = nodes_by_start_position.find(va.first);
        if (ep != nodes_by_end_position.end()
            && sp != nodes_by_start_position.end()) {
            set<Node*>& previous_nodes = ep->second;
            set<Node*>& current_nodes = sp->second;
            for (set<Node*>::iterator n = previous_nodes.begin();
                 n != previous_nodes.end(); ++n) {
                for (set<Node*>::iterator m = current_nodes.begin();
                     m != current_nodes.end(); ++m) {
                    if (node_index.find(*n) != node_index.end()
                        && node_index.find(*m) != node_index.end()
                        && !(previous_nodes.count(*n) && current_nodes.count(*n)
                             && previous_nodes.count(*m) && current_nodes.count(*m))
                        ) {
#ifdef deubg
                        cerr tid << ": connecting previous "
                                 << (*n)->id() << " @end=" << ep->first << " to current "
                                 << (*m)->id() << " @start=" << sp->first << endl;
#endif
                        create_edge(*n, *m);
                    }
                }
            }
        }

        // clean up previous
        while (!nodes_by_end_position.empty() && nodes_by_end_position.begin()->first < va.first) {
            nodes_by_end_position.erase(nodes_by_end_position.begin()->first);
        }

        while (!nodes_by_start_position.empty() && nodes_by_start_position.begin()->first < va.first) {
            nodes_by_start_position.erase(nodes_by_start_position.begin()->first);
        }

        // Now we just have to update where our end was, so the next group of
        // alleles knows if there was any intervening sequence.
        // The (past the) end position is equal to the number of bases not yet used.
        last_variant_end = seq.size() - get_node((*seq_node_ids.rbegin()).second)->sequence().size();

    }

    // Now we're done breaking nodes. This means the node holding the end of the
    // reference sequence can finally be given its mappings for phasings, if
    // applicable.
    if(num_phasings > 0) {
        // What's the last node on the reference path?
        auto last_node_id = (*seq_node_ids.rbegin()).second;
        for(size_t i = 0; i < num_phasings; i++) {
            // Everything visits this last reference node
            paths.append_mapping("_phase" + to_string(i), last_node_id);
        }
    }

    // Put the mapping to the primary path in the graph
    for (auto& p : seq_node_ids) {
        paths.append_mapping(name, p.second);
    }
    // and set the mapping edits
    force_path_match();

    sort();
    compact_ids();

}


}


