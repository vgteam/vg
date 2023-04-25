/**
 * \file
 * constructor.cpp: contains implementations for vg construction functions.
 */
 
#include "constructor.hpp"
#include "utility.hpp"
#include "crash.hpp"
#include "io/load_proto_to_graph.hpp"

#include <IntervalTree.h>

#include <cstdlib>
#include <set>
#include <tuple>
#include <list>
#include <algorithm>
#include <memory>

//#define debug

namespace vg {

    using namespace std;

    void Constructor::trim_to_variable(vector<list<vcflib::VariantAllele>>& parsed_alleles) {

#ifdef debug
        cerr << "Before trimming to variable region:" << endl;
        for (auto& allele : parsed_alleles) {
            cerr << "Allele:" << endl;

            for (auto& edit : allele) {
                cerr << "\tHave " << edit.ref << " -> " << edit.alt << " @ " << edit.position << endl; 
            }
        }
#endif

        // Return the length of perfect match common to all alleles at the left or
        // the right (based on the front parameter). Only looks at the first or last
        // edit in each allele
        auto get_match_count = [&](bool front) -> size_t { 
            // Start with the max possible value
            size_t match_count = numeric_limits<size_t>::max();
            for (auto& allele : parsed_alleles) {
                // Go throught the alleles
                if (allele.empty()) {
                    // No shared ref match possible with an empty allele
                    return 0;
                }
                // Find the edit on the appropriate edge of this alt
                auto& edit = front ? allele.front() : allele.back();
                if (edit.ref != edit.alt) {
                    // This alt has a non-match edit on its edge
                    return 0;
                }

                // Otherwise we have a leading or trailing match so mix in its length
                match_count = min(match_count, edit.ref.size());
            }
            if (match_count == numeric_limits<size_t>::max()) {
                // If nobody exists we have 0 shared matches.
                return 0;
            }
            return match_count;
        };

        for(size_t front_match_count = get_match_count(true); front_match_count > 0; front_match_count = get_match_count(true)) {
            // While we have shared matches at the front
            
#ifdef debug
            cerr << "Edits at the front share " << front_match_count << " match bases and need to be trimmed down" << endl;
#endif
            
            for (auto& allele : parsed_alleles) {
                // Trim each allele
                if (allele.front().ref.size() > front_match_count) {
                    // This perfect match needs to be made shorter
                    auto new_match_string = allele.front().ref.substr(front_match_count);
#ifdef debug
                    cerr << "Trim " << allele.front().ref << " to " << new_match_string
                        << " @ " << allele.front().position << endl;
#endif
                    allele.front().ref = new_match_string;
                    allele.front().alt = new_match_string;
                    
                    // Since we're trimming off the front we need to bump the position up.
                    allele.front().position += front_match_count;
                } else {
                    // This perfect match can be completely eliminated
#ifdef debug
                    cerr << "Drop " << allele.front().ref << " -> " << allele.front().alt
                        << " @ " << allele.front().position << endl;
#endif
                    allele.pop_front();
                }
            }
        }

        for(size_t back_match_count = get_match_count(false); back_match_count > 0; back_match_count = get_match_count(false)) {
            // While we have shared matches at the back
            
#ifdef debug
            cerr << "Edits at the back share " << back_match_count << " match bases and need to be trimmed down" << endl;
#endif
            
            for (auto& allele : parsed_alleles) {
                // Trim each allele
                if (allele.back().ref.size() > back_match_count) {
                    // This perfect match needs to be made shorter
                    auto new_match_string = allele.back().ref.substr(back_match_count);
#ifdef debug
                    cerr << "Trim " << allele.back().ref << " to " << new_match_string
                        << " @ " << allele.back().position << endl;
#endif
                    allele.back().ref = new_match_string;
                    allele.back().alt = new_match_string;
                } else {
                    // This perfect match can be completely eliminated
                    allele.pop_back();

#ifdef debug
                    cerr << "Drop " << allele.back().ref << " -> " << allele.back().alt
                        << " @ " << allele.back().position << endl; 
#endif

                }
            }
        }

#ifdef debug
        cerr << "After trimming to variable region:" << endl;
        for (auto& allele : parsed_alleles) {
            cerr << "Allele: " << endl;
            for (auto& edit : allele) {
                cerr << "\tKept " << edit.ref << " -> " << edit.alt << " @ " << edit.position << endl; 
            }
        }
#endif

    }

    void Constructor::condense_edits(list<vcflib::VariantAllele>& parsed_allele) {
        for(auto i = parsed_allele.begin(); i != parsed_allele.end(); ++i) {
            // Scan through the edits in the alt
            if (i->ref == i->alt) {
                // We can merge into this edit
                auto next = i;
                ++next;

                // We'll use a string stream to generate the combined string
                stringstream combined;
                combined << i->ref;

                while (next != parsed_allele.end() && next->ref == next->alt) {
                    // Glom up all the subsequent matches and remove their nodes.
                    combined << next->ref;
                    next = parsed_allele.erase(next);
                }

                // Put the finished string into the node that led the run
                i->ref = combined.str();
                i->alt = combined.str();
            }
        }
    }

    pair<int64_t, int64_t> Constructor::get_symbolic_bounds(vcflib::Variant var) {
        // TODO: We assume that the variant actually has at least one symbolic alt allele like <INS>.
        // If that is the case, the base at POS must be an anchoring, unmodified base.
        // But you can also have SV tags on something like a CCATG->G right-anchored deletion as long as
        // none of the alleles are symbolic.
        // We really should be calling this on variants that *were* symbolic before canonicalization.
    
        // Move the start 1 base right to account for the required anchor base.
        // This may make us start after the end.
        int64_t start = (int64_t) var.zeroBasedPosition() + 1;
        int64_t end = var.getMaxReferencePos();
        
        return std::make_pair(start, end);
    }


    pair<int64_t, int64_t> Constructor::get_bounds(const vector<list<vcflib::VariantAllele>>& trimmed_variant) {

        // We track the variable site bounds through all the alts
        int64_t variable_start = numeric_limits<int64_t>::max();
        int64_t variable_stop = -1;

        for (auto& trimmed_parts : trimmed_variant) {
            // For every variable core of an alt (which may be empty)
            if (!trimmed_parts.empty()) {
                // We have at least one valid non-match edit on this alt. Expand the range.
                variable_start = min(variable_start, (int64_t) trimmed_parts.front().position - 1);
                variable_stop = max(variable_stop, (int64_t) (trimmed_parts.back().position - 1 + trimmed_parts.back().ref.size() - 1));
            }
        }

        #ifdef debug
        cerr << "Edits for variant run " << variable_start << " through " << variable_stop
            << " ref length " << (variable_stop - variable_start + 1) << endl;
        #endif

        return make_pair(variable_start, variable_stop);
    }
    
    bool Constructor::sanitize_sequence_in_place(string& sequence, const string* sequence_name, size_t sequence_start_offset, const vcflib::Variant* variant) const {
        
        bool made_change = false;
        
        // Make sure the input sequence is upper-case
        string uppercase_sequence = toUppercase(sequence);
        if (uppercase_sequence != sequence) {
            // We had to make a change
            if (warn_on_lowercase) {
                if (variant) {
                    // We are warning about a variant (alt)
                    if (!lowercase_warned_alt) {
                        #pragma omp critical (cerr)
                        {
                            cerr << "warning:[vg::Constructor] Lowercase characters found in "
                                 << "variant; coercing to uppercase:\n" << *const_cast<vcflib::Variant*>(variant) << endl;
                            lowercase_warned_alt = true;
                        }
                    }
                } else {
                    // What sequence are we complaining about?
                    string name_to_warn = sequence_name ? *sequence_name : "DNA sequence";
                    #pragma omp critical (cerr)
                    {
                        // Note that the pragma also protects this mutable map that we update
                        if (!lowercase_warned_sequences.count(name_to_warn)) {
                            // We haven't warned about this sequence yet
                            cerr << "warning:[vg::Constructor] Lowercase characters found in "
                                << name_to_warn << "; coercing to uppercase." << endl;
                            lowercase_warned_sequences.insert(name_to_warn);
                        }    
                    }
                }
            }
            // Replace the original
            sequence = std::move(uppercase_sequence);
            made_change = true;
        }
        
        // Make sure all IUPAC codes are Ns
        string n_sequence = allAmbiguousToN(sequence);
        if (n_sequence != sequence) {
            // We had to make a change
            if (warn_on_ambiguous) {
                if (variant) {
                    // We are warning about a variant (alt).
                    // TODO: We used to always bail for IUPAC codes in a
                    // variant allele; do we really want to not?
                    #pragma omp critical (cerr)
                    {
                        cerr << "warning:[vg::Constructor] Unsupported IUPAC ambiguity codes found in "
                             << "variant; coercing to N:\n" << *const_cast<vcflib::Variant*>(variant) << endl;
                    }
                } else {
                    // What sequence are we complaining about?
                    string name_to_warn = sequence_name ? *sequence_name : "DNA sequence";
                    #pragma omp critical (cerr)
                    {
                        // Note that the pragma also protects this mutable map that we update
                        if (!ambiguous_warned_sequences.count(name_to_warn)) {
                            // We haven't warned about this sequence yet
                            cerr << "warning:[vg::Constructor] Unsupported IUPAC ambiguity codes found in "
                                << name_to_warn << "; coercing to N." << endl;
                            ambiguous_warned_sequences.insert(name_to_warn);
                        }    
                    }
                }
            }
            // Replace the original
            sequence = std::move(n_sequence);
            made_change = true;
        }
        
        // TODO: this is like the forth scan of the whole string we do; can we
        // condense this all into one pass?
        if (!allATGCN(sequence)) {
            // We don't know what to do with gaps, and we want to catch
            // complete garbage.
            
            // We would like an example.
            auto it = sequence.begin();
            while (it != sequence.end() && (*it == 'A' || *it == 'T' || *it == 'G' || *it == 'C' || *it == 'N')) {
                ++it;
            }
            
            #pragma omp critical (cerr)
            {
                cerr << "error:[vg::Constructor] unacceptable character ";
                if (it != sequence.end()) {
                    cerr << "\"" << *it << "\" ";
                }
                cerr << "found in ";
                if (sequence_name) {
                    cerr << *sequence_name;
                } else {
                    cerr << "DNA sequence";
                }
                if (it != sequence.end()) {
                    cerr << " at index " << (it - sequence.begin() + sequence_start_offset);
                }
                if (variant) {
                    cerr << " in variant:\n" << *const_cast<vcflib::Variant*>(variant);
                } else {
                    cerr << ".";
                }
                cerr << endl;
                exit(1);
            }
        }
        
        return made_change;
    }

    ConstructedChunk Constructor::construct_chunk(string reference_sequence, string reference_path_name,
        vector<vcflib::Variant> variants, size_t chunk_offset) const {
        
        std::stringstream status_stream;
        status_stream << "constructing chunk " << reference_path_name << ":" << chunk_offset << " length " << reference_sequence.size();

        set_crash_context(status_stream.str());

        #ifdef debug
        cerr << status_stream.str() << endl;
        #endif

        // Make sure the input sequence is upper-case and all IUPAC codes are Ns
        sanitize_sequence_in_place(reference_sequence, &reference_path_name);
        
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
        
        // We also keep separate maps for reference nodes only, for tracing
        // back through inversions. Since when we trace through an inversion we
        // need to visit every node in a run, we don't just care about the
        // bounding IDs. So we store entire copies of the runs. But since the
        // inversions always go backward, we only need them by their end.
        // These are also on-the-end and not past-the-end positions, matching
        // nodes_ending_at.
        map<size_t, vector<Node*>> ref_runs_by_end;
        
        // We don't want to wire inserts to each other, so we have a set of all
        // insert endpoints.
        set<id_t> inserts;

        // We need to wire up our inversions super specially. These hold the
        // end positions for each inversion anchored at the key, for
        // inversions_starting, or visa-versa, for inversions_ending. In other
        // words, the start position is exclusive and the end position is
        // inclusive.
        map<size_t, set<size_t>> inversions_starting;
        map<size_t, set<size_t>> inversions_ending;
        
        // Here we remember deletions that end at paritcular positions in the
        // reference, which are the positions of the last deleted bases. We map from
        // last deleted base to last non-deleted base before the deletion, so we can
        // go look up nodes ending there. Note that these can map to -1.
        map<size_t, set<int64_t>> deletions_ending_at;

        // We also need to track all points at which deletions start, so we can
        // search for the next one when deciding where to break the reference.
        // We store the last NON-deleted base; the deletion arc attaches to the
        // right side of this base! This base is *not* deleted.
        set<int64_t> deletion_starts;

        // We use this to get the next variant
        auto next_variant = variants.begin();

        // We're going to clump overlapping variants together.
        vector<vcflib::Variant*> clump;
        // And we need to rember the highest past-the-end base of anything in the
        // clump, to catch all the overlaps.
        size_t clump_end = 0;

        // We use this to remember path ranks. It will initialize to 0 for new
        // paths.
        map<Path*, size_t> max_rank;

        // We have a utility function to tack a full length perfect match onto a
        // path. We need the node so we can get its length.
        // Automatically fills in rank, starting from 1.
        auto add_match = [&](Path* path, Node* node, bool is_reverse = false) {
            #ifdef debug
            cerr << "Add node " << node->id() << " orientation " << is_reverse
                << " length " << node->sequence().size() << " to path " << path->name() << endl;
            #endif
        
            // Make a mapping for it
            auto* mapping = path->add_mapping();
            mapping->mutable_position()->set_node_id(node->id());
            mapping->mutable_position()->set_is_reverse(is_reverse);

            // Set the rank to the next available rank in the path.
            mapping->set_rank(++max_rank[path]);

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
                size_t piece_count = sequence.size() / max_node_size;
                if(piece_count > 1) {
                    piece_size = min(max_node_size, max(sequence.size() / piece_count, (size_t) 1));
                } else {
                    piece_size = max_node_size;
                }

                // Remember we may have a partial piece at the end.
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

                if (!created.empty()) {
                    // We need to link to the previous node in this run of sequence
                    auto* edge = to_return.graph.add_edge();
                    edge->set_from(created.back()->id());
                    edge->set_to(node->id());
                }

                // Tack the new node on the end of the list
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

            #ifdef debug
            cerr << "Create reference from cursor at " << reference_cursor << " out up to before "
                << target_position << "/" << reference_sequence.size() << endl;
            #endif

            // Don't get out of the chunk
            if (target_position > reference_sequence.size()) {
                #pragma omp critical (cerr)
                cerr << "error:[vg::Constructor] On " << reference_path_name
                     << ", attempted to add reference nodes until position " << target_position
                     << " but reference is only " << reference_sequence.size() << " long" << endl;
                exit(1);
            }
            if (reference_cursor > reference_sequence.size()) {
                #pragma omp critical (cerr)
                cerr << "error:[vg::Constructor] On " << reference_path_name
                     << ", reference cursor is at " << reference_cursor
                     << " but reference is only " << reference_sequence.size() << " long" << endl;
                exit(1);
            }
            
            if (target_position < reference_cursor) {
                // TODO: should this ever happen? Should we be asked to go backward?
                #ifdef debug
                cerr << "Nothing to do! Already sufficient reference!" << endl;
                #endif
                return;
            }

            // Make new nodes for all the sequence we want to add
            auto new_nodes = create_nodes(reference_sequence.substr(reference_cursor, target_position - reference_cursor));

            // Remember the total node length we have scanned through
            size_t seen_bases = 0;

            if (!new_nodes.empty()) {
                // Add the start node to the starting at map.
                // Interior nodes break at locations that are arbitrary, and we know
                // nothign wants to attach to them there. Plus they're already linked to
                // each other and we shouldn't link them again.
                // Remember where it starts and ends along the reference path
                nodes_starting_at[reference_cursor].insert(new_nodes.front()->id());

                for (Node* node : new_nodes) {
                    // Add matches on the reference path for all the new nodes
                    add_match(ref_path, node);

                    // Remember how long that node was so we place the next one right.
                    seen_bases += node->sequence().size();
                }

                // Add the end node to the ending at map.
                nodes_ending_at[reference_cursor + seen_bases - 1].insert(new_nodes.back()->id());
                
                // Save the whole run for inversion tracing
                #ifdef debug
                cerr << "Create ref run ending at " << reference_cursor + seen_bases - 1 << endl;
                #endif
                ref_runs_by_end[reference_cursor + seen_bases - 1] = std::move(new_nodes);

            }

            // Advance the cursor
            reference_cursor = target_position;
            
            #ifdef debug
            cerr << "Advanced reference cursor for next unmade base to " << reference_cursor << "/" << reference_sequence.size() << endl;
            #endif
            
            if (reference_cursor > reference_sequence.size()) {
                #pragma omp critical (cerr)
                cerr << "error:[vg::Constructor] On " << reference_path_name
                     << ", after adding reference nodes, reference cursor is at " << reference_cursor
                     << " but reference is only " << reference_sequence.size() << " long" << endl;
                exit(1);
            }
        };

        while (next_variant != variants.end() || !clump.empty()) {
            // While there are more variants, or while we have the last clump to do...

            // Group variants into clumps of overlapping variants.
            if (clump.empty() || 
                (next_variant != variants.end() && clump_end > next_variant->zeroBasedPosition() - chunk_offset)) {

                // Either there are no variants in the clump, or this variant
                // overlaps the clump. It belongs in the clump
                clump.push_back(&(*next_variant));
                // It may make the clump longer and necessitate adding more variants.
                // TODO: make sure long SVs don't fall outside chunk
                clump_end = max(clump_end, next_variant->zeroBasedPosition() + next_variant->ref.size() - chunk_offset);

                // Try the variant after that
                next_variant++;
            } else {
                // The next variant doesn't belong in this clump.
                // Handle the clump.
                
                #ifdef debug
                cerr << "Handling clump of " << clump.size() << " variants up to " << (clump_end + chunk_offset) << endl;
                #endif

                // Parse all the variants into VariantAllele edits

                // This holds a map from Variant pointer to a vector of lists
                // of VariantAllele edits, one list per non-ref allele of the
                // variant.
                map<vcflib::Variant*, vector<list<vcflib::VariantAllele>>> parsed_clump;

                // This determines the order we will process variants in. We use it
                // to sort variants in a clump by hash for the purposes of assigning
                // IDs.
                map<string, vcflib::Variant*> variants_by_name;

                // This holds the min and max values for the starts and ends of
                // edits in each variant that are actual change-making edits. These
                // are in chunk coordinates. They are only populated if a variant
                // has a variable region. Equal start and end indicate a 1-base region.
                vector<IntervalTree<int64_t, vcflib::Variant*>::interval> variable_intervals;

                // This holds the min and max values for starts and ends of edits
                // not removed from the clump. These are in chunk coordinates.
                int64_t first_edit_start = numeric_limits<int64_t>::max();
                int64_t last_edit_end = -1;

                // We'll fill this with any variants that should be ignored,
                // out of the clump. This is better than erasing out of a
                // vector.
                set<vcflib::Variant*> skipped;
                
                for (size_t var_num = 0; var_num < clump.size(); var_num++) {
                    // For each variant in the clump
                    vcflib::Variant* variant = clump[var_num];
#ifdef debug
                    cerr << "Handling clump variant " << var_num << "/" << clump.size() << " @ " << variant->zeroBasedPosition() << endl;
#endif
                
                    // Since we make the fasta reference uppercase, we do the VCF too (otherwise vcflib gets mad).
                    // We set this if we modify the variant and vcflib needs to reindex it.
                    bool reindex = false;
                    // We set this if we skipped the variant
                    bool skip_variant = false;
                    for (size_t i = 0; i < variant->alt.size(); i++) {
                        auto& alt = variant->alt[i];
                        // Process all the alts and not the ref
                        if (alt == "*") {
                            // This is a newer VCF feature we don't support,
                            // but not a broken file.
                            #pragma omp critical (cerr)
                            {
                                cerr << "warning:[vg::Constructor] Unsupported allele \"*\" found in "
                                     << "variant, skipping variant:\n" << *variant << endl;
                            }
                            skipped.insert(variant);
                            skip_variant = true;
                            break;
                        }
                        // Sanitize the alt of Ns and lower case characters,
                        // and ensure what remains is something we can use, not
                        // a symbolic SV.
                        bool modified = sanitize_sequence_in_place(alt, nullptr, 0, variant);
                        if (modified) {
                            // Also update the copy in alleles
                            variant->alleles[i + 1] = alt;
                        }
                        reindex |= modified;
                    }
                    if (skip_variant) {
                        // Move to the next variant
                        continue;
                    }
                    // Also process the reference, but blame problems on the reference
                    if (sanitize_sequence_in_place(variant->ref, &reference_path_name, variant->zeroBasedPosition())) {
                        // Also update the copy in alleles
                        variant->alleles[0] = variant->ref;
                        reindex = true;
                    }
                    if (reindex) {
                        // Redo the indexing
                        variant->updateAlleleIndexes();
                    }

                    // Check the variant's reference sequence to catch bad VCF/FASTA pairings
                    auto expected_ref = reference_sequence.substr(variant->zeroBasedPosition() - chunk_offset, variant->ref.size());
                    if(variant->ref != expected_ref) {
                    // TODO: report error to caller somehow
                        #pragma omp critical (cerr)
                        cerr << "error:[vg::Constructor] Variant/reference sequence mismatch: " << variant->ref
                            << " vs pos: " << variant->position << ": " << expected_ref << "; do your VCF and FASTA coordinates match?"<< endl
                            << "Variant: " << *variant << endl;
                            cerr << "zero ind: " << variant->zeroBasedPosition() << " 1-indexed: " << variant->position << endl;
                        exit(1);
                    }
                    
                    // No variants should still be symbolic at this point.
                    // Either we canonicalized them into base-level sequence, or we rejected them when making the clump.
                    // If they had IUPAC codes in them we should have fixed that already too.
                    if (variant->isSymbolicSV()) {
                        #pragma omp critical (cerr)
                        {
                            cerr << "error:[vg::Constructor] On " << reference_path_name << " @ " << variant->zeroBasedPosition()
                                 << ", variant appears to be a symbolic SV, but all variants should have already been converted to explicit sequence edits." << endl;
                            cerr << "error:[vg::Constructor] Offending variant: " << *variant << endl;
                        }
                        exit(1);
                    }
                    // If variants have SVTYPE set, though, we will still use that info instead of the base-level sequence.

                    // Name the variant and place it in the order that we'll
                    // actually construct nodes in (see utility.hpp)
                    string variant_name = make_variant_id(*variant);
                    if (variants_by_name.count(variant_name)) {
                        // Some VCFs may include multiple variants at the same
                        // position with the same ref and alt. We will only take the
                        // first one.
                        #pragma omp critical (cerr)
                        cerr << "warning:[vg::Constructor] Skipping duplicate variant with hash " << variant_name
                            << " at " << variant->sequenceName << ":" << variant->position << endl;
                        skipped.insert(variant);
                        continue;
                    }

                    variants_by_name[variant_name] = variant;

                    // We need to parse the variant into alts, each of which is a
                    // series of VariantAllele edits. This holds the full alt allele
                    // string and the edits needed to make it. The VariantAlleles
                    // completely cover the alt, and some of them may be perfect
                    // matches to stretches of reference sequence. Note that the
                    // reference allele of the variant won't appear here.

                    map<string, vector<vcflib::VariantAllele>> alternates;
                    
                    // Decide if we should parse (i.e. align) the variant alts.
                    // We don';t want to do it if the alignment would be too big.
                    bool can_parse = !flat;
                    if (can_parse) {
                        if (variant->isSymbolicSV()) {
                            // All the variants are probably canonicalized by now and
                            // probably should not look like SVs. And vcflib is
                            // smart enough to give us flat alts when we ask to
                            // parse something that still does look like an SV.
                            // But we still probably want the flat alt
                            // postprocessing, so go with flat alts here.
                            can_parse = false;
                        } else {
                            // We (no longer?) have symbolic alleles, so just
                            // bail if any allele is too long. Only ref and alt
                            // fields will be filled in with sequence; alleles
                            // field may still be symbolic.
                            if (variant->ref.size() > max_parsed_variant_size) {
                                // Ref is too long. Handle as flat.
                                can_parse = false;
                            } else {
                                for (auto& a : variant->alt) {
                                    if (a.size() > max_parsed_variant_size) {
                                        // This alt is too long. Handle as flat.
                                        can_parse = false;
                                        break;
                                    }
                                }
                            }
                        }
                    }
                    
                    
                    if (can_parse) {
                        // Do alignments to parse the alleles
                        alternates = variant->parsedAlternates();
                    } else {
                        alternates = variant->flatAlternates();
                        // if we can, remove the 1bp "standard" base that's added at the beginning of indels
                        if (this->trim_indels){
                            for (auto& v : alternates) {
                                for (auto& a : v.second) {
                                    if (a.ref[0] == a.alt[0]) {
                                        a.ref = a.ref.substr(1);
                                        a.alt = a.alt.substr(1);
                                        ++a.position;
                                    }
                                }
                            }
                        }
                    }
                    
                    // Get the variable bounds in VCF space for all the trimmed alts of this variant
                    // Note: we still want bounds for SVs, we just have to get them differently
                    std::pair<int64_t, int64_t> bounds;
                    
                    if (!variant->canonical){
                        // The variant did not have to be canonicalized.
                        // We will process the variant as a normal variant, based on its ref and alt sequences.
                        
                        for (auto &kv : alternates) {
                            // For each alt in the variant

                            if (kv.first == variant->ref)
                            {
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
                            auto &alt_parts = parsed_clump[variant][alt_index];

                            #ifdef debug
                            cerr << "Non-ref allele " << alt_index << endl;
                            #endif

                            // Copy all the VariantAlleles into the list
                            alt_parts.assign(kv.second.begin(), kv.second.end());

                            // Condense adjacent perfect match edits, so we only break
                            // matching nodes when necessary.
                            condense_edits(alt_parts);
                        }
                        // Trim the alts down to the variant's (possibly empty) variable
                        // region
                        trim_to_variable(parsed_clump[variant]);
                        
                        // The bounds are determined from that
                        bounds = get_bounds(parsed_clump[variant]);
                    } else {
                        // We are going to make the edit specified by the SV
                        // tags, instead of whatever edits are implied by
                        // aligning the ref and alt sequences.
                        #ifdef debug
                        cerr << "Use SV tags to define edit for: " << *variant << endl;
                        #endif
                        // For now, only permit one allele for SVs
                        // in the future, we'll build out VCF lib to fix this.
                        // TODO build out vcflib to fix this.
                        // TODO: Warn if we are lopping anything off!
                        
                        // We need to make sure parsed_clump[variant] has an entry for each allele.
                        // But the contents won't matter since this is an SV.
                        parsed_clump[variant].resize(variant->alt.size());
                        
                        // The bounds are determined symbolically
                        bounds = get_symbolic_bounds(*variant);
                    }
                
                    if (bounds.first != numeric_limits<int64_t>::max() || bounds.second != -1) {
                        // There's a (possibly 0-length) variable region
                        bounds.first -= chunk_offset;
                        bounds.second -= chunk_offset;
                        
                        if (alt_paths && bounds.second >= bounds.first) {
                            // The variant covers a nonempty part of the
                            // reference, and we will need to find it for alt
                            // path generation.

                            // Save the bounds for making reference node path visits
                            // inside the ref allele of the variable region.
                            variable_intervals.push_back(IntervalTree<int64_t, vcflib::Variant*>::interval(bounds.first, bounds.second, variant));
                        }

                        #ifdef debug
                        if (bounds.first < first_edit_start) {
                            cerr << "Expanded first_edit_start to " << bounds.first << " with " << *variant << endl;
                        }
                        if (bounds.second > last_edit_end) {
                            cerr << "Expanded last_edit_end to " << bounds.second << " with " << *variant << endl;
                        }
                        #endif

                        // Expand bounds for the variable region of the chunk as a whole
                        first_edit_start = min(first_edit_start, bounds.first);
                        last_edit_end = max(last_edit_end, bounds.second);
                    } else {
                        // There's no variation in this variant.
                        #pragma omp critical (cerr)
                        cerr << "warning:[vg::Constructor] Skipping variant with no sequence changes at "
                            << variant->sequenceName << ":" << variant->position << endl;
                        skipped.insert(variant);}
                }
                
                if (skipped.size() == clump.size()) {
                    // We skipped all the variants in the clump. Kick back up
                    // to clump building.
                    clump.clear();
                    clump_end = 0;
                    continue;
                }

                // Otherwise, we have to have some non-ref material in the
                // clump, even if it occupies 0 reference space.
                if (first_edit_start == numeric_limits<int64_t>::max() || last_edit_end == -1) {
                    // Sometimes we still see this, so make a report of the offending variants.
                    #pragma omp critical (cerr)
                    {
                        cerr << "error:[vg::Constructor] got improperly bounded region " << first_edit_start << " to " << last_edit_end << " for edits of clump of " << clump.size() << " variants, of which " << skipped.size() << " were skipped." << endl;
                        for (vcflib::Variant* v : clump) {
                            if (!skipped.count(v)) {
                                cerr << "Unskipped variant: " << *v << endl;
                            }
                        }
                    }
                    exit(1);
                }

                #ifdef debug
                cerr << "edits run between " << first_edit_start << " and " << last_edit_end << endl;
                #endif

                // Index the variants in the clump by the reference region they overlap
                IntervalTree<int64_t, vcflib::Variant*> variable_interval_tree(std::move(variable_intervals));

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
                unordered_map<vcflib::Variant*, Path*> variant_ref_paths;
                
                // This holds alt Path pointers and the inversions (start, end)
                // that they need to trace through in their inverted
                // orientations. They can't be traced until the corresponding
                // reference nodes have been made. This can be resolved at the
                // clump level.
                vector<tuple<Path*, size_t, size_t>> inversion_trace_queue;

                for (auto& kv : variants_by_name) {
                    // For each variant in the clump, sorted by name
                    auto& variant_name = kv.first;
                    auto* variant = kv.second;
                    
                    #ifdef debug
                    cerr << "Process variant " << variant_name << " @ " << variant->zeroBasedPosition()
                        << " with " << parsed_clump[variant].size() << " alts" << endl;
                    #endif

                    if (alt_paths) {
                        // Declare its ref path straight away.
                        // We fill in the ref paths after we make all the nodes for the edits.
                        variant_ref_paths[variant] = to_return.graph.add_path();
                        variant_ref_paths[variant]->set_name("_alt_" + variant_name + "_0");
                    }

                    for (size_t alt_index = 0; alt_index < parsed_clump[variant].size(); alt_index++) {                
                        // For each non-ref alt in the parsed variant

                        // Name the alt after the number that this allele has.
                        // We have to bump the allele index because the first alt is 0.
                        string alt_name = "_alt_" + variant_name + "_" + to_string(alt_index + 1);

                        // There should be a path named after it.
                        Path* alt_path = nullptr;
                        if (alt_paths) {
                            alt_path = to_return.graph.add_path();
                            alt_path->set_name(alt_name);
                        }

                        // SV HAX
                        if (this->do_svs && variant->hasSVTags() && variant->canonical) {
                            // This is an SV
                            
                            #ifdef debug
                            cerr << "Process alt " << (alt_index + 1) << " of variant " << variant_name << " as an SV" << endl;
                            #endif
                            
                            string sv_type = variant->info.at("SVTYPE").at(0);

                            if (sv_type == "INS") {
                            
                                // For an insertion, the created node will
                                // start at the next base after the base used
                                // to position the variant, just like we would
                                // do for a non-SV insertion with the same POS
                                // anchored on that base.
                                auto e_start = variant->zeroBasedPosition() - chunk_offset + 1;
                                
                                // The insertion will "end" at the position
                                // *before* that, so the things it is getting
                                // inserted before can link up with it. TODO:
                                // Respect the END tag if it says something
                                // different.
                                auto e_end = e_start - 1;
                            
                                auto inserted_sequence = variant->info.at("SEQ").at(alt_index);
                            
                                // Identify our created node run with a key, in
                                // case it exists already somehow.
                                auto key = make_tuple(e_start, "", inserted_sequence);

                                if (created_nodes.count(key) == 0) {
                                    // Create insertion sequence nodes
                                    vector<Node*> node_run = create_nodes(inserted_sequence);

                                    #ifdef debug
                                    cerr << "Inserted node " << node_run.front()->id() << " starts at " << e_start << endl;
                                    cerr << "Inserted node " << node_run.back()->id() << " ends at " << e_end << endl;
                                    #endif
                                    
                                    nodes_starting_at[e_start].insert(node_run.front()->id());
                                    nodes_ending_at[e_end].insert(node_run.back()->id());

                                    inserts.insert(node_run.front()->id());
                                    inserts.insert(node_run.back()->id());

                                    created_nodes[key] = node_run;

                                    if (alt_paths) {
                                        for (Node* node : created_nodes[key]) {
                                            // Add a visit to each node we created/found in
                                            // order to the path for this alt of this
                                            // variant.
                                            add_match(alt_path, node);
                                        }
                                    }
                                }
                            } else if (sv_type == "DEL") {
                                // Deletions also start after the base used
                                // to anchor them, so you can keep the same
                                // POS between SV and non-SV
                                // representations. The END is inclusive.
                                int64_t arc_start = (int64_t) variant->zeroBasedPosition() - chunk_offset; 
                                size_t arc_end = std::stol(variant->info.at("END").at(alt_index)) - chunk_offset - 1;
                               
                                #ifdef debug
                                cerr << "Deletion arc runs " << arc_start << " to " << arc_end << endl;
                                #endif

                                deletions_ending_at[arc_end].insert(arc_start);
                                deletion_starts.insert(arc_start);
                                
                                // No alt path mappings necessary for the
                                // deletion; the existence of the empty
                                // path is sufficient.
                                
                            } else if (sv_type == "INV"){
                                // Handle inversions
                                // We only need reference nodes, plus two arcs
                                // one from the inverted sequence's beginning to the sequence following
                                // its last node and
                                // one from the end of the sequence preceding the inversion to the back 
                                // of the inverted sequence's last node.
                                
                                // Inversions also require a left anchoring base, according to the spec.
                                //
                                // "If any of the ALT alleles is a symbolic
                                // allele (an angle-bracketed ID String “<ID>”)
                                // then the padding base is required and POS
                                // denotes the coordinate of the base preceding
                                // the polymorphism."
                                //
                                
                                // The END is still inclusive.
                                
                                // Start at the anchoring base
                                int64_t inv_start = (int64_t) variant->zeroBasedPosition() - chunk_offset;
                                size_t inv_end = std::stol(variant->info.at("END").at(alt_index)) - chunk_offset - 1;

                                #ifdef debug
                                cerr << "Inversion arcs connect right of " << inv_start << " and right of " << inv_end << endl;
                                #endif

                                inversions_starting[inv_start].insert(inv_end);
                                inversions_ending[inv_end].insert(inv_start);
                                
                                if (alt_paths) {
                                    // We need to make alt path entries through this inverted sequence, backward.
                                    // But we don't have the reference nodes created yet.
                                    // So we queue them up
                                    inversion_trace_queue.emplace_back(alt_path, inv_start, inv_end);
                                }
                            } else {
                                // Unknown or unsupported SV type
                                cerr << "warning:[vg::Constructor]: unrecognized SV type " << sv_type << endl;
                            }
                        } else {
                            // This is not an SV
                            #ifdef debug
                            cerr << "Process alt " << (alt_index + 1) << " of variant " << variant_name << " as an ordinary variant" << endl;
                            #endif
                        
                            for (vcflib::VariantAllele& edit : parsed_clump[variant][alt_index]) {
                                // For each VariantAllele used by the alt
                                #ifdef debug
                                cerr << "Apply " << edit.ref << " -> " << edit.alt << " @ " << edit.position << endl;
                                #endif

                                if (edit.alt != "") {
                                    // This is a visit to a node for the alt
                                    // We need a key to see if a node has been made for this edit already
                                    auto key = make_tuple(edit.position - 1 - chunk_offset, edit.ref, edit.alt);

                                    if (created_nodes.count(key) == 0) {
                                        // We don't have a run of nodes for this edit, so make one.
                                        vector<Node*> node_run = create_nodes(edit.alt);

                                        // Compute where the edit starts and ends in local chunk coordinates
                                        auto edit_start = edit.position - 1 - chunk_offset;
                                        auto edit_end = edit.position - 1 - chunk_offset + edit.ref.size() - 1;


                                        #ifdef debug
                                        cerr << "Created nodes running " << edit_start << " to " << edit_end << endl;
                                        #endif

                                        // Remember where the first one starts and the last one ends, for wiring up later.
                                        nodes_starting_at[edit_start].insert(node_run.front()->id());
                                        nodes_ending_at[edit_end].insert(node_run.back()->id());

                                        if (edit.ref == edit.alt) {
                                            // This edit is a no-op and so the node we just created is a reference run.
                                            // These can be necessary if insertions and deletions are part of the same record.
                                            // Remember the whole node run for inversion tracing
                                            #ifdef debug
                                            cerr << "Create ref run ending at " << edit_end << endl;
                                            #endif
                                            ref_runs_by_end[edit_end] = node_run;
                                        }

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

                                    if (alt_paths) {
                                        for (Node* node : created_nodes[key]) {
                                            // Add a visit to each node we created/found in
                                            // order to the path for this alt of this
                                            // variant.
                                            add_match(alt_path, node);
                                        }
                                    }

                                } else if (edit.ref != "") {
                                    // It's a deletion (and not a weird ""->"" edit).

                                    // Add an entry to the deletion arcs
                                    // What is the end position (last deleted)
                                    // Take the 0-based edit position, remove the chunk offset,
                                    // advance the ref bases, and then back up 1 base to the last deleted ref base.
                                    size_t arc_end = (edit.position - 1) - chunk_offset + edit.ref.size() - 1;
                                    // What is the before-the-beginning anchoring position (last non-deleted, may be -1)
                                    int64_t arc_start = (int64_t) edit.position - 1 - chunk_offset - 1;


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
                                    
                                    // No alt path mappings necessary for the
                                    // deletion
                                }

                            }
                        }

                    }
                }

                // Then after you finish all the alts, add the ref nodes that
                // haven't already been created, breaking wherever something can
                // come in or out.

                // We need a function to work that out

                /**
                 * Find the next breakpoint.
                 *
                 * Takes in the search position, like the position of the next
                 * un-made reference base. So searches from the left edge of
                 * the passed inclusive position.
                 *
                 * Finds the next required breakpoint within this clump, after
                 * the given position, given created nodes and deletions that
                 * already exist.
                 *
                 * Returns the inclusive position of the base to the left of
                 * this breakpoint, so the breakpoint is immediately to the
                 * right of the base at the returned position. This means that
                 * sometimes, such as if the next piece of the reference would
                 * be 1 bp long, this function will return the same value it
                 * was passed.
                 */
                auto next_breakpoint_after = [&](size_t position) -> size_t {
                    // If nothing else, we're going to break at the end of the last
                    // edit in the clump.
                    size_t to_return = last_edit_end;
                    
                    #ifdef debug
                    cerr << "Next breakpoint after " << position << " must be at or before " << last_edit_end << endl;
                    #endif

                    // See if any nodes are registered as starting after our
                    // position. They'll all start before the end of the clump, and
                    // we don't care if they start at our position since that
                    // breakpoint would be to the left and already happened.
                    auto next_start_iter = nodes_starting_at.upper_bound(position);

                    if(next_start_iter != nodes_starting_at.end()) {
                        // If we found something, walk back where the breakpoint
                        // needs to be so we break before that node starts.
                        to_return = min(to_return, next_start_iter->first - 1);
                        #ifdef debug
                        cerr << "Next node starts at " << next_start_iter->first - 1 << endl;
                        #endif
                    }

                    // See if any nodes are registered as ending at or after our
                    // position. We do care if they end at our position, since that
                    // means we need to break right here because that
                    // breakpoint would be to the right.
                    auto next_end_iter = nodes_ending_at.lower_bound(position);

                    if(next_end_iter != nodes_ending_at.end()) {
                        // If we found something, we need to break where that node
                        // ends.
                        to_return = min(to_return, next_end_iter->first );
                        #ifdef debug
                        cerr << "Next node ends at " << next_end_iter->first << endl;
                        #endif
                    }

                    // See if any deletions are registered as ending at or
                    // after here. If the deletion ends here, this is the last
                    // base deleted, and that creates a breeakpoint to our
                    // right.
                    auto deletion_end_iter = deletions_ending_at.lower_bound(position);

                    if(deletion_end_iter != deletions_ending_at.end()) {
                        // If we found something, walk back where the breakpoint
                        // needs to be so we break before the node after the
                        // deletion starts.
                        to_return = min(to_return, deletion_end_iter->first);
                        #ifdef debug
                        cerr << "Next deletion ends by deleting " << deletion_end_iter->first << endl;
                        #endif
                    }
                    
                    // See if any deletions are known to start at or after this
                    // base. We care about exact hits now, because deletions break
                    // to the right of the base they start at, since we are
                    // storing the base that the deletion arc attaches to the
                    // right side of.
                    auto deletion_start_iter = deletion_starts.lower_bound(position);
                    // We don't need to worry about -1s here. They won't be found
                    // with lower_bound on a size_t.

                    if(deletion_start_iter != deletion_starts.end()) {
                        // If we found something, walk back where the breakpoint
                        // needs to be so we break at the position the deletion
                        // needs to leave from.
                        to_return = min(to_return, (size_t)*deletion_start_iter);
                        #ifdef debug
                        cerr << "Next deletion starts after " << *deletion_start_iter << endl;
                        #endif
                    }

                    // Check to see if any inversions' last (largest
                    // coordinate) inverted bases are at or after this point
                    // Inversions break the reference twice, much like deletions.
                    // Since we store the final base that is inverted, and the
                    // inversion creates a breakpoint on the right side of that
                    // base, we care about exact hits.
                    auto inv_end_iter = inversions_ending.lower_bound(position);
                    if (inv_end_iter != inversions_ending.end()){
                        to_return = min(to_return, (size_t) inv_end_iter->first);
                        #ifdef debug
                        cerr << "Next inversion ends after inverting " << inv_end_iter->first << endl;
                        #endif
                    }

                    // Also look at inversions' first (smallest coordinate) bases.
                    // Inversions break just after the anchor the base they start at,
                    // so we care about exact hits and use lower_bound.
                    auto inv_start_iter = inversions_starting.lower_bound(position);
                    if (inv_start_iter != inversions_starting.end()){
                        to_return = min(to_return, (size_t) inv_start_iter->first);
                        #ifdef debug
                        cerr << "Next inversion starts after " << inv_start_iter->first << endl;
                        #endif
                    }
                    
                    #ifdef debug
                    cerr << "Selected " << to_return << " as breakpoint" << endl;
                    #endif

                    return to_return;

                };

                // Note that in some cases (i.e. pure inserts) we may not need any
                // ref nodes at all. Also note that in other cases (variants with
                // exterior matches) some ref nodes overlapping the variant may not
                // really belong on the ref path for the variant, because the alt
                // path for the variant starts/end further in.

                while (reference_cursor < last_edit_end + 1) {
                    // Until we hit the end

                    // Find where the next node run must end to attach to stuff
                    size_t next_end = next_breakpoint_after(reference_cursor);
                    
                    #ifdef debug
                    cerr << "Creating reference nodes from " << reference_cursor << " out to "
                        << next_end << "/" << reference_sequence.size() << endl;
                    #endif
                    
                    if (reference_cursor > reference_sequence.size()) {
                        #pragma omp critical (cerr)
                        cerr << "error:[vg::Constructor] On " << reference_path_name
                             << ", adding reference to last edit end, reference cursor is at " << reference_cursor
                             << " but reference is only " << reference_sequence.size() << " long" << endl;
                        exit(1);
                    }
                    if (next_end > reference_sequence.size()) {
                        #pragma omp critical (cerr)
                        cerr << "error:[vg::Constructor] On " << reference_path_name
                             << ", adding reference to last edit end, next end is at " << next_end
                             << " but reference is only " << reference_sequence.size() << " long" << endl;
                        exit(1);
                    }

                    // We need to have a reference node/run of nodes (which may have
                    // already been created by a reference match) between where the
                    // cursor is and where the next breakpoint has to be.
                    // This is the sequence it should have.
                    string run_sequence = reference_sequence.substr(reference_cursor, next_end - reference_cursor + 1);

                    // We need a key to see if a node (run) has been made for this sequece already
                    auto key = make_tuple(reference_cursor, run_sequence, run_sequence);
                    auto representative_nodes = created_nodes.find(key);
                    if (representative_nodes == created_nodes.end()) {
                        // We don't have a run of ref nodes up to the next break, so make one
                        vector<Node*> node_run = create_nodes(run_sequence);

                        // Remember where the first one starts and the last one ends, for wiring up later.
                        nodes_starting_at[reference_cursor].insert(node_run.front()->id());
                        nodes_ending_at[next_end].insert(node_run.back()->id());
                        
                        // Remember the whole node run for inversion tracing
                        #ifdef debug
                        cerr << "Create ref run ending at " << next_end << endl;
                        #endif
                        ref_runs_by_end[next_end] = node_run;
                        

#ifdef debug
                        cerr << "Created reference nodes running " << reference_cursor << " to " << next_end << endl;
#endif

                        // Save it in case any other alts also have this edit.
                        representative_nodes = created_nodes.insert(representative_nodes, {key, node_run});
                    } else {
#ifdef debug
                        cerr << "Reference nodes at  " << reference_cursor << " for constant " << run_sequence.size() << " bp sequence " << run_sequence << " already exist" << endl;
#endif
                    }

                    for (Node* node : representative_nodes->second) {
                        // Add a reference visit to each node we created/found
                        add_match(ref_path, node);
                    }

                    if (!representative_nodes->second.empty() && alt_paths) {
                        // Ref paths from other variants may need to visit these new nodes.
                        auto overlapping_intervals = variable_interval_tree.findOverlapping(reference_cursor, reference_cursor);
                        for (auto& interval : overlapping_intervals) {
                            if (interval.start <= reference_cursor && interval.stop >= reference_cursor && !skipped.count(interval.value)) {
                                // For each variant we might also be part of the ref allele of

                                // For unique variants that actually differ from reference,
                                // if this run of nodes starts within the variant's variable region...
                                // (We know if it starts in the variable region it has to
                                // end in the variant, because the variable region ends with
                                // a node break)

                                if (variant_ref_paths.count(interval.value) == 0) {
                                    // All unique variants ought to have a ref path created
                                    cerr << "error:[vg::Constructor] no ref path for " << *interval.value << endl;
                                    exit(1);
                                }
                                
                                for (Node* node : representative_nodes->second) {
                                    // Add a match along the variant's ref allele path
                                    add_match(variant_ref_paths[interval.value], node);
                                }
                            }
                        }
                    }

                    // Advance the reference cursor to after this run of reference nodes
                    reference_cursor = next_end + 1;
                    
                    if (reference_cursor > reference_sequence.size()) {
                        #pragma omp critical (cerr)
                        cerr << "error:[vg::Constructor] On " << reference_path_name
                             << ", after adding reference to last edit end, reference cursor is at " << reference_cursor
                             << " but reference is only " << reference_sequence.size() << " long" << endl;
                        exit(1);
                    }
                    
                    // Keep going until we have created reference nodes through to
                    // the end of the clump's interior edits.
                }

                // Now we have gotten through all the places where nodes start, before the end of the clump.
                
                for (auto& to_trace : inversion_trace_queue) {
                    // Now that all the ref nodes exist, create the path entries for inversion alt paths.
                    auto& alt_path = get<0>(to_trace);
                    auto& inv_start = get<1>(to_trace);
                    auto& inv_end = get<2>(to_trace);
                    
                    // We will walk this cursor back from the end of the
                    // inversion to the start, going backward through runs of
                    // reference nodes that end here.
                    // It will track the first base from right to left that we have yet to cover with our path.
                    // Our inversion end is inclusive and that base is inverted, so start past there.
                    int64_t inv_end_cursor = inv_end;
                    
                    
                    while (inv_end_cursor > inv_start) {
                        #ifdef debug
                        cerr << "Inversion cursor at " << inv_end_cursor << endl;
                        #endif
                    
                        // Get the next ref run on the right that the inversion has to visit
                        auto& trailing_run = ref_runs_by_end.at(inv_end_cursor);
                        
                        for (auto it = trailing_run.rbegin(); it != trailing_run.rend(); it++) {
                            // For each node in the run in reverse order
                            Node* node = *it;
                            
                            #ifdef debug
                            cerr << "Reverse node " << node->id() << endl;
                            #endif
                            
                            // Add a match to this node in its reverse orientation, since we are inverting.
                            add_match(alt_path, node, true);
                            
                            // Advance the cursor left after visiting the node
                            inv_end_cursor -= node->sequence().size();
                        }
                    }
                    
                    #ifdef debug
                    cerr << "Added inversion alt path from " << inv_end << " back to " << inv_start << " and arrived at "
                        << inv_end_cursor << endl;
                    #endif
                    
                    if (inv_end_cursor != inv_start) {
                        // Make sure we did it right
                        #pragma omp critical (cerr)
                        cerr << "error:[vg::Constructor] On " << reference_path_name << " near " << reference_cursor
                             << ", inversion end cursor " << inv_end_cursor << " did not reach inversion start " << inv_start << endl;
                        exit(1);
                    }
                
                }

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
                    
#ifdef debug
                    cerr << "Node " << right_node << " can start at " << kv.first << endl;
#endif

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

                    

                    // Now we do the deletions. We want to allow daisy-chaining
                    // deletions.

                    // We compose a set of deletion start points.
                    set<int64_t> possible_starts;

                    // We also keep a list of unexplored deletion end points to chain from.
                    list<int64_t> possible_ends;
                    possible_ends.push_back(kv.first - 1);

                    // And a set of explored ones
                    set<int64_t> explored_ends;

                    while (!possible_ends.empty()) {
                        // Find an unexplored place where we can find more daisy-
                        // chained deletions ending.
                        int64_t deletion_end = possible_ends.front();
                        possible_ends.pop_front();

#ifdef debug
                        cerr << deletions_ending_at[deletion_end].size() << " deletions end by deleting " << deletion_end << endl;
#endif
    
                        for (auto& deletion_start : deletions_ending_at[deletion_end]) {
                            // For every deletion start that can end there.

                            // Note that we can delete from there to our node along
                            // transitive deletions.
                            possible_starts.insert(deletion_start);

                            // We can daisy chain from deletions that end by
                            // deleting the anchor base that this deletion starts at.
                            int64_t possible_end = deletion_start;

                            if(chain_deletions && possible_end > 0 && !explored_ends.count(possible_end)) {
                                // Queue it up if not already queued. If we aren't
                                // chaining deletions, we'll only look at the starts
                                // accessible from the root of our searcj.
                                possible_ends.push_back(possible_end);
                                explored_ends.insert(possible_end);
                            }
                        }
                    }

                    for (auto& deletion_start : possible_starts) {
                        // For everywhere a deletion can start that comes to here

                        if (deletion_start == -1) {
                            // This node ought to be exposed on the left actually.
                            to_return.left_ends.insert(right_node);

                        } else {
                            // The deletion doesn't go all the way to the left edge
                            // but actually starts at a place where there are nodes.

                            for (auto& left_node : nodes_ending_at[deletion_start]) {
                                // For every node that the deletion could
                                // anchor from (because they end exactly where
                                // it starts with its anchor)

                                if (inserts.count(left_node)) {
                                    // Don't let an inserted node happen just before a deletion.
#ifdef debug
                                    cerr << "Skip insertion-deletion edge " << left_node << " -> " << right_node << endl;
#endif
                                    continue;
                                }

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
                    
                    #ifdef debug
                    for (auto& kv2 : inversions_starting) {
                        cerr << "Inversion can start after " << kv2.first << endl;
                    }
                    for (auto& kv2 : inversions_ending) {
                        cerr << "Inversion can end by inverting " << kv2.first << endl;
                    }
                    #endif

                    // Now do the inversions.
                    
                    // What do we hook up to the start of right_node, which starts at kv.first?
                    // For any inversions that end by inverting kv.first - 1, we hook up the starts of anything where the inversion started.
                    if (inversions_ending.count(kv.first - 1)) {
                        for (auto& inv_start : inversions_ending[kv.first - 1]) {
                            // For each inversion start position corresponding to inversions ending by inverting this base
                            
#ifdef debug
                            cerr << "Inversion ending by inverting " << kv.first - 1 << " is anchored at " << inv_start
                                << " after which " << nodes_starting_at[inv_start + 1].size() << " nodes start" << endl;
#endif
                            
                            for (auto& n : nodes_starting_at[inv_start + 1]) {
#ifdef debug
                                cerr << "Node " << n << " can start at " << (inv_start + 1) << " where inversion first inverts. "
                                 << "So link its start to our right_node's start." << endl;
#endif
                                
                                // For each node that starts at the inversion start position, link it up inverted.
                                auto* e = to_return.graph.add_edge();
                                e->set_from(n);
                                e->set_from_start(true);
                                e->set_to(right_node);
                                e->set_to_end(false);
                                
#ifdef debug
                                cerr << "Invert " << n << " to " << right_node << endl;
#endif
                            }
                        }
                    }
                
                
                }
                
                // Inversions continue with another loop over the left nodes
                for (auto& left_node : nodes_ending_at[kv.first - 1]) {
                
                    // What do we hook up to the end of left_node, which ends right before kv.first?
                    // For any inversions anchoring where this node ends, we hook up the ends of everything that is at where the inversion ends.
                    
                    if (inversions_starting.count(kv.first - 1)) {
                        for (auto& inv_end : inversions_starting[kv.first - 1]) {
                            // For each inversion end position corresponding to inversions starting by inverting this base
                            
#ifdef debug
                            cerr << "Inversion starting by inverting " << kv.first << " and anchored at " << (kv.first - 1)
                                << " can end at " << inv_end << " where " << nodes_ending_at[inv_end].size() << " nodes end" << endl;
#endif
                            
                            for (auto& n : nodes_ending_at[inv_end]) {
#ifdef debug
                                cerr << "Node " << n << " can end at " << inv_end << " where inversion does. "
                                 << "So link its end to " << left_node << "'s end at anchor point " << (kv.first - 1) << endl;
#endif
                            
                                // For each node that ends at that inversion end position, link it up inverted.
                                auto* e = to_return.graph.add_edge();
                                e->set_from(left_node);
                                e->set_from_start(false);
                                e->set_to(n);
                                e->set_to_end(true);
                                
#ifdef debug
                                cerr << "Invert " << left_node << " to " << n << endl;
#endif
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

        for(auto& deletion_start : deletions_ending_at[reference_sequence.size() - 1]) {
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
        
        // Filter out any empty variant Paths.
        // First stop pointing to them.
        max_rank.clear();
        // We have this many nonempty paths at the start of the collection at
        // the start of every loop.
        size_t nonempty_paths = 0;
        while (nonempty_paths < to_return.graph.path_size()) {
            if (to_return.graph.path(nonempty_paths).mapping_size() == 0) {
                // This is empty of mappings so swap it to the end and remove it.
                to_return.graph.mutable_path()->SwapElements(nonempty_paths, to_return.graph.path_size() - 1);
                to_return.graph.mutable_path()->RemoveLast();
                // Leave our cursor where it is; we have to check the element we swapped to here.
            } else {
                // This is nonempty so advance the cursor.
                nonempty_paths++;
            }
        }
        
        clear_crash_context();
        return to_return;
    }

    void Constructor::construct_graph(string vcf_contig, FastaReference& reference, VcfBuffer& variant_source,
        const vector<FastaReference*>& insertions, const function<void(Graph&)>& callback) {

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
        
#ifdef debug
        cerr << "building contig for chunk of reference " << reference_contig << " in interval " << leading_offset << " to " << reference_end << endl;
#endif

        // Set up a progress bar thhrough the chromosome
        create_progress("building graph for " + vcf_contig, reference_end - leading_offset);

        // Scan through variants until we find one that is on this contig and in this region.
        // If we're using an index, we ought to already be at the right place.
        variant_source.fill_buffer();
        while(variant_source.get() && (variant_source.get()->sequenceName != vcf_contig ||
                    variant_source.get()->zeroBasedPosition() < leading_offset ||
                    variant_source.get()->zeroBasedPosition() + variant_source.get()->ref.size() > reference_end)) {
            // This variant comes before or ends after our region

            // Discard variants that come out that are outside our region
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

        // And we need to do the same for ranks on the reference path? What's the
        // max rank used?
        size_t max_ref_rank = 0;

        // Whenever a chunk ends with a single node, with no edges to the end
        // or non-reference paths, we separate it out and buffer it here,
        // because we may need to glue it together with subsequent leading
        // nodes with no edges to the start or non-reference paths, to
        // eliminate spurious node breaks at chunk boundaries. 
        Node last_node_buffer;

        // Sometimes we need to emit single node reference chunks gluing things
        // together
        auto emit_reference_node = [&](Node& node) {

            if (node.id() == 0) {
                // Don't emit nonexistent nodes
                #pragma omp critical (cerr)
                cerr << "error:[vg::Constructor] On " << vcf_contig << " near " << chunk_start
                     << ", tried to produce a reference node without an ID" << endl;
                exit(1);
            }

            // Make a single node chunk for the node
            Graph chunk;
            *(chunk.add_node()) = node;

            // It needs a primary path mapping.
            Path* path = chunk.add_path();
            path->set_name(reference_contig);
            Mapping* mapping = path->add_mapping();
            mapping->mutable_position()->set_node_id(node.id());
            // With a rank
            mapping->set_rank(++max_ref_rank);
            // And an edit
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
            
            // If we have a single head or tail with no outer edges or
            // non-reference paths, we will fill in the ID here.
            nid_t head_id = 0;
            nid_t tail_id = 0;

            if (last_node_buffer.id() != 0 && chunk.left_ends.size() == 1) {
                // We have a single dangling node buffered that we can rewrite.
                // So we also want to know if we have a head to combine it with.
                // And it looks like we might.
                head_id = *chunk.left_ends.begin();
                #ifdef debug
                cerr << "Node " << head_id << " might be mergable with buffered node " << last_node_buffer.id() << endl;
                #endif
            }

            if (chunk.right_ends.size() == 1) {
                // We always need to know if we can buffer the tail.
                // Name this node a candidate tail.
                tail_id = *chunk.right_ends.begin();
                #ifdef debug
                cerr << "Node " << tail_id << " might be a tail we can buffer" << endl;
                #endif
            }

            for (auto& edge : chunk.graph.edge()) {            
                // Go through all edges and kick out head and tail candidates if they have any on the outside.
                if (head_id && ((edge.from() == head_id && edge.from_start()) || (edge.to() == head_id && !edge.to_end()))) {
                    // Edge connects to the start of the head candidate, so it fails.
                    #ifdef debug
                    cerr << "Node " << head_id << " has an edge to its left and so can't merge" << endl;
                    #endif
                    head_id = 0;
                }
                if (tail_id && ((edge.from() == tail_id && !edge.from_start()) || (edge.to() == tail_id && edge.to_end()))) {
                    // Edge connects to the end of the tail candidate, so it fails.
                    #ifdef debug
                    cerr << "Node " << tail_id << " has an edge to its right and so can't merge" << endl;
                    #endif
                    tail_id = 0;
                }
            }

            for (size_t i = 1; (head_id != 0 || tail_id != 0) && i < chunk.graph.path_size(); i++) {
                // Go through all paths other than the reference (which is 0)
                auto& path = chunk.graph.path(i);
                
                // Check the first and last steps on the path to see if they
                // touch our head/tail nodes. Other steps can't touch them
                // because of the edge restrictions we already checked.
                auto check_mapping = [&](size_t mapping_index) {
                    nid_t touched_node = path.mapping(mapping_index).position().node_id();
                    if (touched_node == head_id) {
                        #ifdef debug
                        cerr << "Node " << head_id << " is visited by path " << path.name() << " and so can't merge" << endl;
                        #endif
                        head_id = 0;
                    }
                    if (touched_node == tail_id) {
                        #ifdef debug
                        cerr << "Node " << tail_id << " is visited by path " << path.name() << " and so can't merge" << endl;
                        #endif
                        tail_id = 0;
                    }
                };

                // Sometimes the first and last step are the same step!
                if (path.mapping_size() > 0) {
                    // We have a first step
                    check_mapping(0);
                    if (path.mapping_size() > 1) {
                        // We have a distinct last step
                        check_mapping(path.mapping_size() - 1);
                    }
                }
            }


            if (last_node_buffer.id() != 0 && head_id != 0) {
                // We have a last node from the last chunk that we want to glom onto
                // this chunk, and we have a node to do it with.

                // We want to merge it with the single source node for this
                // chunk. But depending on the variant structure it may not be
                // the first node generated (because we generate variant alt
                // material first, and a variant may lead the chunk). So we do
                // a linear scan.
                
                // This seems slow, but actually shouldn't be: most of the time
                // we succeed on the first try, the whole process is linear in
                // graph size anyway, and we never have to scan through more
                // than a variant's worth of nodes.
                
                // We will fill this in
                Node* mutable_first_node = nullptr;
                for (size_t i = 0; i < chunk.graph.node_size(); i++) {
                    // Look at each node in turn
                    mutable_first_node = chunk.graph.mutable_node(i);
                    
                    if (mutable_first_node->id() == head_id) {
                        // We found the left end we want
                        break;
                    }
                }
                
                if (mutable_first_node == nullptr || mutable_first_node->id() != head_id) {
                    // Make sure we found it
                    #pragma omp critical (cerr)
                    cerr << "error:[vg::Constructor] On " << reference_contig
                         << ", could not find node " << head_id << endl;
                    exit(1);
                }

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
                if (path->name() != reference_contig) {
                    #pragma omp critical (cerr)
                    cerr << "error:[vg::Constructor] Expected path " << reference_contig
                         << " but found path " << path->name() << endl;
                    exit(1);
                }
                // Then the first mapping
                Mapping* mapping = path->mutable_mapping(0);
                if (mapping->position().node_id() != mutable_first_node->id()) {
                    #pragma omp critical (cerr)
                    cerr << "error:[vg::Constructor] On " << reference_contig
                         << ", expected node " << mutable_first_node->id()
                         << " but found node " << mapping->position().node_id() << endl;
                    exit(1);
                }
                if (mapping->edit_size() != 1) {
                    #pragma omp critical (cerr)
                    cerr << "error:[vg::Constructor] On " << reference_contig
                         << " at node " << mapping->position().node_id()
                         << ", expected 1 edit but found " << mapping->edit_size() << endl;
                    exit(1);
                }
                // Then the only edit
                Edit* edit = mapping->mutable_edit(0);
                // Correct its length
                edit->set_from_length(mutable_first_node->sequence().size());
                edit->set_to_length(mutable_first_node->sequence().size());
            } else if (last_node_buffer.id() != 0) {
                // There's no single leading node on this next chunk that we
                // are free to rewrite, but we still have a single trailing
                // node to emit.

                // Emit it
                emit_reference_node(last_node_buffer);
                // Clear it
                last_node_buffer = Node();
            }

            if (tail_id != 0) {
                // We need to pull out the last node in the chunk. Note that it may
                // also be the first node in the chunk...

                // We know it's the last node in the graph
                last_node_buffer = chunk.graph.node(chunk.graph.node_size() - 1);
                
                if (last_node_buffer.id() != tail_id) {
                    #pragma omp critical (cerr)
                    cerr << "error:[vg::Constructor] On " << reference_contig
                         << ", could not find right end for node " << last_node_buffer.id() << endl;
                    exit(1);
                }

                // Remove it
                chunk.graph.mutable_node()->RemoveLast();

                // Find the primary path
                Path* path = chunk.graph.mutable_path(0);
                if (path->name() != reference_contig) {
                    #pragma omp critical (cerr)
                    cerr << "error:[vg::Constructor] Expected path " << reference_contig
                         << " but found path " << path->name() << endl;
                    exit(1);
                }
                // Then drop last mapping, which has to be to this node
                if (path->mapping_size() == 0) {
                    #pragma omp critical (cerr)
                    cerr << "error:[vg::Constructor] On " << reference_contig
                         << ", found empty path" << endl;
                    exit(1);
                }
                if (path->mapping(path->mapping_size() - 1).position().node_id() != last_node_buffer.id()) {
                    #pragma omp critical (cerr)
                    cerr << "error:[vg::Constructor] On " << reference_contig
                         << ", expected last node" << last_node_buffer.id()
                         << " but found " << path->mapping(path->mapping_size() - 1).position().node_id() << endl;
                    exit(1);
                }
                path->mutable_mapping()->RemoveLast();

                // Update its ID separately, since it's no longer in the graph.
                last_node_buffer.set_id(last_node_buffer.id() + max_id);
                
                #ifdef debug
                cerr << "Buffered final node becomes: " << last_node_buffer.id() << endl;
                #endif
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

                    // Set the rank.
                    // TODO: we're just clobbering the ref path ranks that were generated in chunk construction.
                    mapping->set_rank(++max_ref_rank);
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

        bool do_external_insertions = false;
        FastaReference* insertion_fasta;


        if (insertions.size() == 1){
            // If we only get one fasta file for insertions, we'll
            // open it and take all insert sequences from there.
            do_external_insertions = true;
            insertion_fasta = insertions[0];
        }
        else if (insertions.size() > 1){
            // if we have more than one insertion fasta file, we can pull
            // sequences from the vcf:fasta pair (i.e. the same index in the vectors).
            do_external_insertions = true;
            cerr << "Passing multiple insertion files not implemented yet." << endl
            << "Please try combining all of your insertions fastas into one file." << endl;
            exit(1);
        }   
        else{
            // We didn't get any insertion fastas, so we will only handle
            // those with seqs in the vcf.

        }
        
        #ifdef debug 
        cerr << "Handling run of variants starting in region..." << endl;
        #endif

        while (variant_source.get() && variant_source.get()->sequenceName == vcf_contig &&
               variant_source.get()->zeroBasedPosition() >= leading_offset &&
               variant_source.get()->zeroBasedPosition() <= reference_end) {
               
            // For each variant that begins inside our region

            // Skip variants that don't fit in our range
            // (maybe there's one that does fit after, so we continue checking)
            if (variant_source.get()->zeroBasedPosition() + variant_source.get()->ref.size() > reference_end) {
                variant_source.handle_buffer();
                variant_source.fill_buffer();
                continue;
            }
                
            // While we have variants we want to include
            auto vvar = variant_source.get();

            // Fix the variant's canonical flag being uninitialized.
            vvar->canonical = false;

            // We need to decide if we want to use this variant. By default we will use all variants.
            bool variant_acceptable = true;
            
            if (vvar->alt.empty()) {
                // Variants with no alts are unimportant
                variant_acceptable = false;
            } else {
                for (string& alt : vvar->alt) {
                    // Variants with "." alts can't be processsed.
                    // TODO: handle remaining alts.
                    if (alt == ".") {
                        variant_acceptable = false;
                        if (vvar->alt.size() > 1) {
                            // Warn if there are more alts we will miss because of skipping
                            #pragma omp critical (cerr)
                            cerr << "warning:[vg::Constructor] Variant with '.' among multiple alts being skipped: "
                                << *vvar << endl;
                        }
                        break;
                    }
                }
            }
            
            if (variant_acceptable) {
                if (vvar->isSymbolicSV()) {
                    // We have a symbolic not-all-filled-in alt.
                    // We need to be processed as a symbolic SV
                    // It also might just have IUPAC codes in it and an SVTYPE and thus look symbolic.
                
                    if (this->do_svs) {
                        // We are actually going to try to handle this SV.
                        
                        if (vvar->alt.size() > 1) {
                            // vcflib will refuse to canonicalize multi-alt SVs.
                            #pragma omp critical (cerr)
                            {
                                if (!multiallelic_sv_warned) {
                                    cerr << "warning:[vg::Constructor] Multiallelic SVs cannot be canonicalized by vcflib; skipping variants like: " << *vvar << endl;
                                    multiallelic_sv_warned = true;
                                }
                            }
                            variant_acceptable = false;
                        }
                        
                        if (variant_acceptable) {
                            // Canonicalize the variant and see if that disqualifies it.
                            // This also takes care of setting the variant's alt sequences.
                            variant_acceptable = vvar->canonicalizable() && vvar->canonicalize(reference, insertions, true);
                            if (!variant_acceptable) {
                                #pragma omp critical (cerr)
                                {
                                    if (!uncanonicalizable_sv_warned) {
                                        cerr << "warning:[vg::Constructor] vcflib could not canonicalize some SVs to base-level sequence; skipping variants like: " << *vvar << endl;
                                        uncanonicalizable_sv_warned = true;
                                    }
                                }
                            }
                        }
         
                        if (variant_acceptable) {
                            // Worth checking for bounds problems.
                            // We have seen VCFs where the variant positions are on GRCh38 but the END INFO tags are on GRCh37.
                            // But for inserts the bounds will have the end right before the start, so we have to allow for those.
                            auto bounds = get_symbolic_bounds(*vvar);
                            if (bounds.second + 1 < bounds.first) {
                                #pragma omp critical (cerr)
                                cerr << "warning:[vg::Constructor] SV with end position " << bounds.second
                                    << " significantly before start " << bounds.first << " being skipped (check liftover?): "
                                    << *vvar << endl;
                                variant_acceptable = false;
                            }
                        }
                    } else {
                        // SV handling is off.
                        variant_acceptable = false;
                        
                        // Figure out exactly what to complain about.
                        for (string& alt : vvar->alt) {
                            // Validate each alt of the variant

                            if(!allATGCN(alt)) {
                                // This is our problem alt here.
                                // Either it's a symbolic alt or it is somehow lower case or something.
                                // It could be an IUPAC code, which we can't handle usually.
                                #pragma omp critical (cerr)
                                {
                                    bool warn = true;
                                    if (!alt.empty() && alt[0] == '<' && alt[alt.size()-1] == '>') {
                                        if (symbolic_allele_warnings.find(alt) != symbolic_allele_warnings.end()) {
                                            warn = false;
                                        } else {
                                            symbolic_allele_warnings.insert(alt);
                                        }
                                    }
                                    if (warn) {
                                        cerr << "warning:[vg::Constructor] Unsupported variant allele \""
                                            << alt << "\"; skipping variants like: " << *vvar <<" !" << endl;
                                    }
                                }
                                break;
                            }
                        }
                    }
                }
            }


            if (!variant_acceptable) {
                // Skip variants that have symbolic alleles or other nonsense we can't parse.
                variant_source.handle_buffer();
                variant_source.fill_buffer();
            } else if (!chunk_variants.empty() && chunk_end > vvar->zeroBasedPosition()) {
                // If the chunk is nonempty and this variant overlaps what's in there, put it in too and try the next.
                // TODO: this is a lot like the clumping code...

                // Add it in
                chunk_variants.push_back(*(vvar));
                // Expand out how big the chunk needs to be, so we can get other overlapping variants.
                chunk_end = max(chunk_end, chunk_variants.back().zeroBasedPosition() + chunk_variants.back().ref.size());

                // Try the next variant
                variant_source.handle_buffer();
                variant_source.fill_buffer();

            } else if(chunk_variants.size() < vars_per_chunk && variant_source.get()->zeroBasedPosition() < chunk_start + bases_per_chunk) {
                // Otherwise if this variant is close enough and the chunk isn't too big yet, put it in and try the next.

                // TODO: unify with above code?

                // Add it in
                chunk_variants.push_back(*(vvar));
                // Expand out how big the chunk needs to be, so we can get other overlapping variants.
                chunk_end = max(chunk_end, chunk_variants.back().zeroBasedPosition() + chunk_variants.back().ref.size());

                // Try the next variant
                variant_source.handle_buffer();
                variant_source.fill_buffer();

            } else {
                // This variant shouldn't go in this chunk.

                // Finish the chunk to a point before the next variant, before the
                // end of the reference, before the max chunk size, and after the
                // last variant the chunk contains.
                chunk_end = max(chunk_end,
                        min((size_t ) vvar->zeroBasedPosition(),
                            min((size_t) reference_end,
                                (size_t) (chunk_start + bases_per_chunk))));

                // Get the ref sequence we need
                auto chunk_ref = reference.getSubSequence(reference_contig, chunk_start, chunk_end - chunk_start);

                // Call the construction
                auto result = construct_chunk(chunk_ref, reference_contig, chunk_variants, chunk_start);

                // Wire up and emit the chunk graph
                wire_and_emit(result);

                // Say we've completed the chunk
                update_progress(chunk_end - leading_offset);

                // Set up a new chunk
                chunk_start = chunk_end;
                chunk_end = 0;
                chunk_variants.clear();

                // Loop again on the same variant.
            }
        }

        #ifdef debug
        cerr << "Variants in region depleted, which we know because we found a starting-after-region variant." << endl;
        #endif

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

            // Say we've completed the chunk
            update_progress(chunk_end - leading_offset);

            // Set up a new chunk
            chunk_start = chunk_end;
            chunk_end = 0;
            chunk_variants.clear();
        }

        // All the chunks have been wired and emitted.
        
        if (last_node_buffer.id() != 0) {
            // Now emit the very last node, if any
            emit_reference_node(last_node_buffer);
            // Update the max ID with that last node, so the next call starts at the next ID
            max_id = max(max_id, (id_t) last_node_buffer.id());
        }

        destroy_progress();

    }

    void Constructor::construct_graph(const vector<FastaReference*>& references,
        const vector<vcflib::VariantCallFile*>& variant_files, const vector<FastaReference*>& insertions,
        const function<void(Graph&)>& callback) {

        // Make a map from contig name to fasta reference containing it.
        map<string, FastaReference*> reference_for;
        for (size_t i = 0; i < references.size(); i++) {
            // For every FASTA reference, make sure it has an index
            auto* reference = references[i];
            if (!reference->index) {
                #pragma omp critical (cerr)
                cerr << "error:[vg::Constructor] Reference #" << i << " is missing its index" << endl;
                exit(1);
            }
            for (auto& kv : *(reference->index)) {
                // For every sequence name and index entry, point to this reference
                reference_for[kv.first] = reference;
#ifdef debug
                cerr << "Contig " << kv.first << " is in reference " << i << endl;
#endif
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

#ifdef debug
                cerr << "Make graph for " << vcf_name << " = " << fasta_name << endl;
#endif

                // Also the FASTA reference that has that sequence
                if (!reference_for.count(fasta_name)) {
                    cerr << "[vg::Constructor] Error: \"" << fasta_name << "\" not found in fasta file" <<endl;
                    exit(1);
                }
                FastaReference* reference = reference_for[fasta_name];

                // We'll set this to true if we actually find the VCF that contains
                // the variants for this sequence and successfully build the graph for it.
                bool built_region = false;
                
                for (auto& buffer : buffers) {
                    // For each VCF we are going to read
                    if(!buffer->has_tabix()) {
                        // Die if we don't have indexes for everyone.
                        // TODO: report errors to caller instead.
#pragma omp critical (cerr)
                        cerr << "[vg::Constructor] Error: all VCFs must be indexed when restricting to a region" << endl;
                        exit(1);
                    }

                    // We set this to true if this VCF contains this region.
                    bool in_this_vcf = false;

                    // Try seeking to the right contig/region
                    if (allowed_vcf_regions.count(vcf_name)) {
                        // Seek to just that region (0-based)
                        in_this_vcf = buffer->set_region(vcf_name, allowed_vcf_regions[vcf_name].first,
                                allowed_vcf_regions[vcf_name].second);
                    } else {
                        // Seek to just the whole contig
                        in_this_vcf = buffer->set_region(vcf_name);
                    }

                    if (in_this_vcf) {
                        // This VCF covers the region
                        
                        if (built_region) {
                            // The region has already been built; we are checking for conflicting VCFs and we found one.
                            // TODO: Use them all with some kind of on-the-fly merging version of the variant buffer.
#pragma omp critical (cerr)
                            cerr << "[vg::Constructor] Error: multiple VCFs cover selected region in " << vcf_name
                                << "; merge them before constructing the graph" << endl;
                            exit(1);
                        } else {
                            // This buffer is the one!
                            // Construct the graph for this contig with the FASTA and the VCF.
                            construct_graph(vcf_name, *reference, *buffer, insertions, callback);
                            
                            // Record that we built the region but check the
                            // other VCFs still so we can complain if the user
                            // gave us overlapping VCFs we can't use.
                            built_region = true;
                        }
                    }
                }

                if (!built_region) {
                    // None of the VCFs include variants on this sequence.
                    // Just build the graph for this sequence with no varaints.
                    VcfBuffer empty(nullptr);
                    construct_graph(vcf_name, *reference, empty, insertions, callback);
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
                    
                    if (constructed.count(vcf_contig)) {
                        // We already did this contig. The user must have
                        // passed us overlapping, unmerged VCFs which we can't
                        // support yet.
                        cerr << "[vg::Constructor] Error: multiple VCFs cover " << vcf_contig
                            << "; merge them before constructing the graph" << endl;
                        exit(1);
                    }

                    // Decide what FASTA contig that is and make sure we have it
                    string fasta_contig = vcf_to_fasta(vcf_contig);
                    if (!reference_for.count(fasta_contig)) {
                        cerr << "[vg::Constructor] Error: Reference contig \"" << vcf_contig << "\" in VCF not found in FASTA." << endl;
                        exit(1);
                    }
                    auto* reference = reference_for[fasta_contig];

                    // Construct on it with the appropriate FastaReference for that contig
                    construct_graph(vcf_contig, *reference, *buffer, insertions, callback);
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
                construct_graph(vcf_contig, *reference, empty, insertions, callback);
            }

            // Now we've constructed everything we can. We're done!


        }

    }

    void Constructor::construct_graph(const vector<string>& reference_filenames, const vector<string>& variant_filenames,
        const vector<string>& insertion_filenames, const function<void(Graph&)>& callback) {
        
        vector<unique_ptr<FastaReference>> references;
        for (auto& fasta_filename : reference_filenames) {
            // Open each FASTA file
            FastaReference* reference = new FastaReference();
            references.emplace_back(reference);
            reference->open(fasta_filename);
        }
        
        vector<unique_ptr<vcflib::VariantCallFile>> variant_files;
        for (auto& vcf_filename : variant_filenames) {
            // Make sure each VCF file exists. Otherwise Tabix++ may exit with a non-
            // helpful message.
            
            // We can't invoke stat woithout a place for it to write. But all we
            // really want is its return value.
            struct stat temp;
            if(stat(vcf_filename.c_str(), &temp)) {
                cerr << "error:[Constructor::construct_graph] file \"" << vcf_filename << "\" not found" << endl;
                exit(1);
            }
            vcflib::VariantCallFile* variant_file = new vcflib::VariantCallFile();
            variant_file->parseSamples = false; // Major speedup if there are many samples.
            variant_files.emplace_back(variant_file);
            // TODO: vcflib needs a non-const string for the filename for some reason. Fix that.
            string mutable_filename = vcf_filename;
            variant_file->open(mutable_filename);
            if (!variant_file->is_open()) {
                cerr << "error:[Constructor::construct_graph] could not open" << vcf_filename << endl;
                exit(1);
            }
        }
        
        vector<unique_ptr<FastaReference>> insertions;
        for (auto& insertion_filename : insertion_filenames){
            // Open up those insertion files
            FastaReference* insertion = new FastaReference();
            insertions.emplace_back(insertion);
            insertion->open(insertion_filename);
        }
        
        // Make vectors of just bare pointers
        vector<vcflib::VariantCallFile*> vcf_pointers;
        for(auto& vcf : variant_files) {
            vcf_pointers.push_back(vcf.get());
        }
        vector<FastaReference*> fasta_pointers;
        for(auto& fasta : references) {
            fasta_pointers.push_back(fasta.get());
        }
        vector<FastaReference*> ins_pointers;
        for (auto& ins : insertions){
            ins_pointers.push_back(ins.get());
        }
        
        // Construct the graph.
        construct_graph(fasta_pointers, vcf_pointers, ins_pointers, callback);
    }
    
    void Constructor::construct_graph(const vector<FastaReference*>& references,
        const vector<vcflib::VariantCallFile*>& variant_files, const vector<FastaReference*>& insertions,
        MutablePathMutableHandleGraph* destination) {
        
        vg::io::load_proto_to_graph(destination, [&](const function<void(Graph&)>& callback) {
            // Start a load of a stream of Protobuf Graphs, and when we get the
            // callback to handle them, construct into it.
            construct_graph(references, variant_files, insertions, callback);
        });
        
        // Now we did the construction and all the Graph chunks have been saved.
        // TODO: Refactor everything to not go through Graph chunks?
    }
    
    void Constructor::construct_graph(const vector<string>& reference_filenames, const vector<string>& variant_filenames,
        const vector<string>& insertion_filenames, MutablePathMutableHandleGraph* destination) {
        
        vg::io::load_proto_to_graph(destination, [&](const function<void(Graph&)>& callback) {
            // Start a load of a stream of Protobuf Graphs, and when we get the
            // callback to handle them, construct into it.
            construct_graph(reference_filenames, variant_filenames, insertion_filenames, callback);
        });
        
        // Now we did the construction and all the Graph chunks have been saved.
        // TODO: Refactor everything to not go through Graph chunks?
        
        // TODO: Deduplicate with the version that takes already-opened files somehow...
    }
}


