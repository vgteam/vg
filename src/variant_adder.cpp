#include "variant_adder.hpp"
#include "banded_global_aligner.hpp"
#include "mapper.hpp"
#include "algorithms/prune.hpp"

//#define debug

namespace vg {

using namespace std;
using namespace vg::io;

VariantAdder::VariantAdder(VG& graph) : graph(graph), sync([&](VG& g) -> VG& {
        // Dice nodes in the graph for GCSA indexing *before* constructing the synchronizer.
        handlealgs::chop(g, max_node_size);
        g.paths.compact_ranks();
        return g;
    }(this->graph)) {
    
    
    graph.paths.for_each_name([&](const string& name) {
        // Save the names of all the graph paths, so we don't need to lock the
        // graph to check them.
        path_names.insert(name);
    });
    
    // Show progress if the graph does.
    show_progress = graph.show_progress;
    
    // Configure the aligner to use a full length bonus
    aligner.full_length_bonus = 5;
}

const VG& VariantAdder::get_graph() const {
    return graph;
}

void VariantAdder::add_variants(vcflib::VariantCallFile* vcf) {
    
#ifdef debug
    // Count our current heads and tails
    vector<Node*> to_count;
    
    graph.head_nodes(to_count);
    size_t head_expected = to_count.size();
    to_count.clear();
    graph.tail_nodes(to_count);
    size_t tail_expected = to_count.size();
    to_count.clear();
    
    cerr << "Starting with heads: " << head_expected << " and tails: " << tail_expected << endl;
#endif
    
    // Make a buffer
    WindowedVcfBuffer buffer(vcf, variant_range);
    
    // Count how many variants we have done
    size_t variants_processed = 0;
    
    // Keep track of the previous contig name, so we know when to change our
    // progress bar.
    string prev_path_name;
    
    // We report when we skip contigs, but only once.
    set<string> skipped_contigs;
    
    while(buffer.next()) {
        // For each variant in its context of nonoverlapping variants
        vcflib::Variant* variant;
        vector<vcflib::Variant*> before;
        vector<vcflib::Variant*> after;
        tie(before, variant, after) = buffer.get_nonoverlapping();
        
        // Where is it?
        auto variant_path_name = vcf_to_fasta(variant->sequenceName);
        
        auto variant_path_offset = (variant->zeroBasedPosition()); // No longer made 0-based by variant buffer!
        
        if (!path_names.count(variant_path_name)) {
            // This variant isn't on a path we have.
            if (ignore_missing_contigs) {
                // That's OK. Just skip it.
                
                if (!skipped_contigs.count(variant_path_name)) {
                    // Warn first
                    
                    // Don't clobber an existing progress bar (which must be over since we must be on a new contig)
                    destroy_progress();
                    cerr << "warning:[vg::VariantAdder] skipping missing contig " << variant_path_name << endl;
                    skipped_contigs.insert(variant_path_name);
                }
                
                continue;
            } else {
                // Explode!
                throw runtime_error("Contig " + variant_path_name + " mentioned in VCF but not found in graph");
            }
        }
        
        if (variant->samples.empty()) {
            // Complain if the variant has no samples. If there are no samples
            // in the VCF, we can't generate any haplotypes to use to add the
            // variants.
            throw runtime_error("No samples in variant at " + variant_path_name +
                ":" + to_string(variant_path_offset) +"; can't make haplotypes");
        }
        
        // Grab the sequence of the path, which won't change
        const string& path_sequence = sync.get_path_sequence(variant_path_name);
    
        // Interlude: do the progress bar
        // TODO: not really thread safe
        if (variant_path_name != prev_path_name) {
            // Moved to a new contig
            prev_path_name = variant_path_name;
            destroy_progress();
            create_progress("contig " + variant_path_name, path_sequence.size());
        }
        update_progress(variant_path_offset);
        
        // Figure out what the actual bounds of this variant are. For big
        // deletions, the variant itself may be bigger than the window we're
        // using when looking for other local variants. For big insertions, we
        // might need to get a big subgraph to ensure we have all the existing
        // alts if they exist.
        
        // Make the list of all the local variants in one vector
        vector<vcflib::Variant*> local_variants = filter_local_variants(before, variant, after);
        
        // Get the unique haplotypes
        auto haplotypes = get_unique_haplotypes(local_variants, &buffer);
        
        // Track the total bp of haplotypes
        size_t total_haplotype_bases = 0;
        
        // Track the total graph size for the alignments
        size_t total_graph_bases = 0;
        
        // How many haplotypes actually pass any haplotype filtering?
        size_t used_haplotypes = 0;
            
#ifdef debug
        cerr << "Have " << haplotypes.size() << " haplotypes for variant "
            << variant->sequenceName << ":" << variant->position << endl;
#endif
        
        // Where does the group of nearby variants start?
        size_t group_start = local_variants.front()->zeroBasedPosition();
        // And where does it end (exclusive)? This is the latest ending point of any variant in the group...
        size_t group_end = local_variants.back()->zeroBasedPosition() + local_variants.back()->ref.size();
        
        // We need to make sure we also grab this much extra graph context,
        // since we count 2 radiuses + flank out from the ends of the group.
        size_t group_width = group_end - group_start;
        
        // Find the center and radius of the group of variants, so we know what graph part to grab.
        size_t overall_center;
        size_t overall_radius;
        tie(overall_center, overall_radius) = get_center_and_radius(local_variants);

        // Get the leading and trailing ref sequence on either side of this
        // group of variants (to pin the outside variants down).

        // On the left we want either flank_range bases, or all the bases before
        // the first base in the group.
        size_t left_context_length = min((int64_t) flank_range, (int64_t) group_start);
        // On the right we want either flank_range bases, or all the bases after
        // the last base in the group. We know nothing will overlap the end of
        // the last variant, because we grabbed nonoverlapping variants.
        size_t right_context_length = min(path_sequence.size() - group_end, (size_t) flank_range);
    
        // Turn those into desired substring bounds.
        // TODO: this is sort of just undoing some math we already did
        size_t left_context_start = group_start - left_context_length;
        size_t right_context_past_end = group_end + right_context_length;
            
#ifdef debug
            cerr << "Original context bounds: " << left_context_start << " - " << right_context_past_end << endl;
#endif

        // Find the reference sequence
        const string& ref_sequence = sync.get_path_sequence(vcf_to_fasta(variant->sequenceName));

        for (auto& haplotype : haplotypes) {
            // For each haplotype
            
            // TODO: since we lock repeatedly, neighboring variants will come in
            // in undefined order and our result is nondeterministic.
            
            // Only look at haplotypes that aren't pure reference.
            bool has_nonreference = false;
            for (auto& allele : haplotype) {
                if (allele != 0) {
                    has_nonreference = true;
                    break;
                }
            }
            if (!has_nonreference) {
                // Don't bother aligning all-ref haplotypes to the graph.
                // They're there already.
#ifdef debug
                cerr << "Skip all-reference haplotype." << endl;
#endif
                continue;
            }
            
#ifdef debug
            cerr << "Haplotype ";
            for (auto& allele_number : haplotype) {
                cerr << allele_number << " ";
            }
            cerr << endl;
#endif

            // This lets us know if we need to walk out more to find matchable sequence
            bool have_dangling_ends;
            do {
                // We need to be able to increase our bounds until we haven't
                // shifted an indel to the border of our context.
                
                // Round bounds to node start and endpoints.
                // This haplotype and all subsequent ones will be aligned with this wider context.
                sync.with_path_index(variant_path_name, [&](const PathIndex& index) {
                    tie(left_context_start, right_context_past_end) = index.round_outward(left_context_start,
                        right_context_past_end);
                });
                
#ifdef debug
                cerr << "New context bounds: " << left_context_start << " - " << right_context_past_end << endl;
#endif
                
                // Recalculate context lengths
                left_context_length = group_start - left_context_start;
                right_context_length = right_context_past_end - group_end;
                
                // Make sure we pull out out to the ends of the contexts
                overall_radius = max(overall_radius, max(overall_center - left_context_start,
                    right_context_past_end - overall_center));
                
                // Get actual context strings
                string left_context = path_sequence.substr(group_start - left_context_length, left_context_length);
                string right_context = path_sequence.substr(group_end, right_context_length);
                
                // Make the haplotype's combined string
                stringstream to_align;
                to_align << left_context << haplotype_to_string(haplotype, local_variants) << right_context;
                
#ifdef debug
                cerr << "Align " << to_align.str() << endl;
#endif

                // Make a request to lock the subgraph, leaving the nodes we rounded
                // to (or the child nodes they got broken into) as heads/tails.
                GraphSynchronizer::Lock lock(sync, variant_path_name, left_context_start, right_context_past_end);
                
#ifdef debug
                cerr << "Waiting for lock on " << variant_path_name << ":"
                    << left_context_start << "-" << right_context_past_end << endl;
#endif
                
                // Block until we get it
                lock_guard<GraphSynchronizer::Lock> guard(lock);
                
#ifdef debug
                cerr << "Got lock on " << variant_path_name << ":"
                    << left_context_start << "-" << right_context_past_end << endl;
#endif            
                    
#ifdef debug
                cerr << "Got " << lock.get_subgraph().length() << " bp in " << lock.get_subgraph().size() << " nodes" << endl;
#endif

#ifdef debug
                ofstream seq_dump("seq_dump.txt");
                seq_dump << to_align.str();
                seq_dump.close();

                sync.with_path_index(variant_path_name, [&](const PathIndex& index) {
                    // Make sure we actually have the endpoints we wanted
                    auto found_left = index.find_position(left_context_start);
                    auto found_right = index.find_position(right_context_past_end - 1);
                    assert(left_context_start == found_left->first);
                    assert(right_context_past_end == found_right->first + index.node_length(found_right));
                    
                    cerr << "Group runs " << group_start << "-" << group_end << endl;
                    cerr << "Context runs " << left_context_start << "-" << right_context_past_end << ": "
                        << right_context_past_end - left_context_start  << " bp" << endl;
                    cerr << "Sequence is " << to_align.str().size() << " bp" << endl;
                    cerr << "Leftmost node is " << found_left->second << endl;
                    cerr << "Leftmost Sequence: " << lock.get_subgraph().get_node(found_left->second.node)->sequence() << endl;
                    cerr << "Rightmost node is " << found_right->second << endl;
                    cerr << "Rightmost Sequence: " << lock.get_subgraph().get_node(found_right->second.node)->sequence() << endl;
                    cerr << "Left context: " << left_context << endl;
                    cerr << "Right context: " << right_context << endl;
                    
                    lock.get_subgraph().for_each_node([&](Node* node) {
                        // Look at nodes
                        if (index.by_id.count(node->id())) {
                            cerr << "Node " << node->id() << " at " << index.by_id.at(node->id()).first
                                << " orientation " << index.by_id.at(node->id()).second << endl;
                        } else {
                            cerr << "Node " << node->id() << " not on path" << endl;
                        }
                    });
                    
                    if (lock.get_subgraph().is_acyclic()) {
                        cerr << "Subgraph is acyclic" << endl;
                    } else {
                        cerr << "Subgraph is cyclic" << endl;
                    }
                });
#endif
                
                // Work out how far we would have to unroll the graph to account for
                // a giant deletion. We also want to account for alts that may
                // already be in the graph and need unrolling for a long insert.
                size_t max_span = max(right_context_past_end - left_context_start, to_align.str().size());
                
                // Do the alignment, dispatching cleverly on size
                Alignment aln = smart_align(lock.get_subgraph(), lock.get_endpoints(), to_align.str(), max_span);
                
#ifdef debug
                cerr << "Postprocessed: " << pb2json(aln) << endl;
#endif
                
                // Look at the ends of the alignment
                assert(aln.path().mapping_size() > 0);
                auto& last_mapping = aln.path().mapping(aln.path().mapping_size() - 1);
                assert(last_mapping.edit_size() > 0);
                auto& last_edit = last_mapping.edit(last_mapping.edit_size() - 1);
                auto& first_mapping = aln.path().mapping(0);
                assert(first_mapping.edit_size() > 0);
                auto& first_edit = first_mapping.edit(0);
                
                // Assume they aren't dangling
                have_dangling_ends = false;
                
                if (!edit_is_match(first_edit) && left_context_start > 0) {
                    // Actually the left end is dangling, so try looking left
                    have_dangling_ends = true;
                    left_context_start--;
#ifdef debug
                    cerr << "Left end dangled!" << endl;
#endif
                }
                
                if (!edit_is_match(last_edit) && right_context_past_end < ref_sequence.size()) {
                    // Actually the right end is dangling, so try looking right
                    have_dangling_ends = true;
                    right_context_past_end++;
#ifdef debug
                    cerr << "Right end dangled!" << endl;
#endif
                }
                
                if (!have_dangling_ends) {
                    
                    // Make this path's edits to the original graph. We don't need to do
                    // anything with the translations.
                    lock.apply_full_length_edit(aln.path(), max_node_size);
                    
                    // Count all the bases in the haplotype
                    total_haplotype_bases += to_align.str().size();
                    // Record the size of graph we're aligning to in bases
                    total_graph_bases += lock.get_subgraph().length();
                    // Record the haplotype as used
                    used_haplotypes++;
                    
                    
                } else {
#ifdef debug
                    cerr << "Expand context and retry" << endl;
#endif
                }
                // If we have dangling ends, we try again with our expanded context
                
            } while (have_dangling_ends);
        }
        
        if (print_updates && (variants_processed++ % 1000 == 0 || true)) {
            #pragma omp critical (cerr)
            cerr << "Variant " << variants_processed << ": " << haplotypes.size() << " haplotypes at "
                << variant->sequenceName << ":" << variant->position << ": "
                << (used_haplotypes ? (total_haplotype_bases / used_haplotypes) : 0) << " bp vs. "
                << (used_haplotypes ? (total_graph_bases / used_haplotypes) : 0) << " bp haplotypes vs. graphs average" << endl;
        }
        
        
#ifdef debug
        // Count our current heads and tails
        graph.head_nodes(to_count);
        size_t head_count = to_count.size();
        cerr << "Heads: ";
        for(auto head : to_count) {
            cerr << head->id() << " ";
        }
        cerr << endl;
        to_count.clear();
        graph.tail_nodes(to_count);
        cerr << "Tails: ";
        for(auto tail : to_count) {
            cerr << tail->id() << " ";
        }
        cerr << endl;
        size_t tail_count = to_count.size();
        to_count.clear();
        
        cerr << "Heads count: " << head_count << ", Tail count: " << tail_count << endl;

        if (head_count != head_expected || tail_count != tail_expected) {
            cerr << "Error! Count mismatch!" << endl;
            // Bail out but serialize the graph
            return;
        }
#endif
        
    }

    // Clean up after the last contig.
    destroy_progress();
    
}

void VariantAdder::align_ns(vg::VG& graph, Alignment& aln) {
    for (size_t i = 0; i < aln.path().mapping_size(); i++) {
        // For each mapping
        auto* mapping = aln.mutable_path()->mutable_mapping(i);
        // Start at its start position (and copy it)
        Position here = mapping->position();
        for (size_t j = 0; j < mapping->edit_size(); j++) {
            // For each edit
            auto* edit = mapping->mutable_edit(j);
            
#ifdef debug
            cerr << "Edit " << pb2json(*edit) << " at " << pb2json(here) << endl;
#endif
            
            if (edit_is_sub(*edit) && is_all_n(edit->sequence())) {
                // The edit is a substitution of N, so it might actually be an N/N match.
                
#ifdef debug
                cerr << "\tIs all Ns" << endl;
#endif
                
                // Is the source string also Ns?
                
                // Grab the node
                auto* node = graph.get_node(here.node_id());
                string original;
                if (here.is_reverse()) {
                    // Pull out the reverse complement of this bit counted from the end
                    original = reverse_complement(node->sequence().substr(
                        node->sequence().size() - here.offset() - edit->from_length(),
                        edit->from_length()));
                } else {
                    // Pull out this bit counted from the start
                    original = node->sequence().substr(here.offset(), edit->from_length());
                }
                
#ifdef debug
                cerr << "\tOriginal: " << original << endl;
#endif
                
                if (is_all_n(original)) {
                    // Yep it's an N for N replacement. Turn it into a match.
                    edit->set_sequence("");
                }
                    
            }
            // If the edit's sequence is a sub for Ns, and the sequence being
            // substituted is Ns, strip the sequence.
            
            
            // Update the position to consume the sequence used by this edit.
            here.set_offset(here.offset() + edit->from_length());
        }
    }
}

Alignment VariantAdder::smart_align(vg::VG& graph, pair<NodeSide, NodeSide> endpoints, const string& to_align, size_t max_span) {

    // We need this fro reverse compelmenting alignments
    auto node_length_function = [&](id_t id) {
        return graph.get_node(id)->sequence().size();
    };

    // Decide what kind of alignment we need to do. Whatever we pick,
    // we'll fill this in.
    Alignment aln;
    
#ifdef debug
    cerr << "Consider " << to_align.size() << " x " << graph.length() << " problem" << endl;
#endif
    
    if (to_align.size() <= whole_alignment_cutoff && graph.length() < whole_alignment_cutoff) {
        // If the graph and the string are short, do a normal banded global
        // aligner with permissive banding and the whole string length as
        // band padding. We can be inefficient but we won't bring down the
        // system.

#ifdef debug
        cerr << "\tUse full-scale " << to_align.size() << " x " << graph.length() << " alignment" << endl;
#endif
        
        // Do the alignment in both orientations
        
        // Align in the forward orientation using banded global aligner, unrolling for large deletions.
        aln = graph.align(to_align, &aligner, true, false, 0, false, false, true, 0, 0, max_span);
        // Align in the reverse orientation using banded global aligner, unrolling for large deletions.
        // TODO: figure out which way our reference path goes through our subgraph and do half the work.
        // Note that if we have reversing edges and a lot of unrolling, we might get the same alignment either way.
        Alignment aln2 = graph.align(reverse_complement(to_align), &aligner, true, false,
            0, false, false, true, 0, 0, max_span);
        
        // Note that the banded global aligner doesn't fill in identity.
        
#ifdef debug
        cerr << "\tScores: " << aln.score() << " fwd vs. " << aln2.score() << " rev" << endl;
#endif
            
        if (aln2.score() > aln.score()) {
            // The reverse alignment is better. But spit it back in the
            // forward orientation.
            aln = reverse_complement_alignment(aln2, node_length_function);
        }

#ifdef debug
        cerr << "\tSubgraph: " << pb2json(graph.graph) << endl;            
        cerr << "\tAlignment: " << pb2json(aln) << endl;
#endif
        
    } else {
        // Either the graph or the sequence to align is too big to just
        // throw in to the banded aligner with big bands.
        
        // First try the endpoint alignments and see if they look like the whole thing might be in the graph.
        
        
        // We need to figure out what bits we'll align
        string left_tail; 
        string right_tail;
        
        if (to_align.size() <= pinned_tail_size) {
            // Each tail will just be the whole string
            left_tail = to_align;
            right_tail = to_align;
        } else {
            // Cut off the tails
            left_tail = to_align.substr(0, pinned_tail_size);
            right_tail = to_align.substr(to_align.size() - pinned_tail_size);
        }
        
        // We don't want to try to align against truly massive graphs with
        // gssw because we can overflow. We also know our alignments need to
        // be near the ends of the extracted graph, so there's no point
        // aligning to the middle.
        
        // Extract one subgraph at each end of the big subgraph we're
        // aligning to. Since we know where we extracted the original
        // subgraph from, this is possible.
        VG left_subgraph;
        VG right_subgraph;
        left_subgraph.add_node(*graph.get_node(endpoints.first.node));
        right_subgraph.add_node(*graph.get_node(endpoints.second.node));
        graph.expand_context_by_length(left_subgraph, left_tail.size() * 2);
        graph.expand_context_by_length(right_subgraph, right_tail.size() * 2);
    
#ifdef debug
        cerr << "\tAttempt two smaller " << left_tail.size() << " x " << left_subgraph.length()
            << " and " << right_tail.size() << " x " << right_subgraph.length() << " alignments" << endl;
#endif
        
        // Do the two pinned tail alignments on the forward strand, pinning
        // opposite ends.
        Alignment aln_left = left_subgraph.align(left_tail, &aligner, true, false,
            0, true, true, false, 0, max_span);
        Alignment aln_right = right_subgraph.align(right_tail, &aligner, true, false,
            0, true, false, false, 0, max_span);
                
        if (aln_left.path().mapping_size() < 1 ||
            aln_left.path().mapping(0).position().node_id() != endpoints.first.node ||
            aln_left.path().mapping(0).edit_size() < 1 || 
            !edit_is_match(aln_left.path().mapping(0).edit(0))) {
        
            // The left alignment didn't start with a match to the correct
            // endpoint node. Try aligning it in reverse complement, and
            // pinning the other end.
            
#ifdef debug
            if (aln_left.path().mapping_size() < 1) {
                cerr << "\tLeft end alignment is null" << endl;
            }
            else {
                cerr << "\tLeft end initial mapping does not match expected end: " << pb2json(aln_left.path().mapping(0)) << endl;
            }
            cerr << "\tRealigning left end to reverse complement" << endl;
#endif
            
            // TODO: what if we have an exact palindrome over a reversing
            // edge, and we keep getting the same alignment arbitrarily no
            // matter how we flip the sequence to align?
            aln_left = reverse_complement_alignment(left_subgraph.align(reverse_complement(left_tail), &aligner,
                                                                        true, false, 0, true, false, false, 0, max_span), node_length_function);
                
        }
        
        // It's harder to do the same checks on the right side because we
        // can't just look at 0. Go find the rightmost mapping and edit, if
        // any.
        const Mapping* rightmost_mapping = nullptr;
        const Edit* rightmost_edit = nullptr;
        if (aln_right.path().mapping_size() > 0) {
            rightmost_mapping = &aln_right.path().mapping(aln_right.path().mapping_size() - 1);
        }
        if (rightmost_mapping != nullptr && rightmost_mapping->edit_size() > 0) {
            rightmost_edit = &rightmost_mapping->edit(rightmost_mapping->edit_size() - 1);
        }
        
        if (rightmost_mapping == nullptr ||
            rightmost_mapping->position().node_id() != endpoints.second.node ||
            rightmost_edit == nullptr || 
            !edit_is_match(*rightmost_edit)) {
        
            // The right alignment didn't end with a match to the correct
            // endpoint node. Try aligning it in reverse complement and
            // pinning the other end.
            
#ifdef debug
            if (rightmost_mapping == nullptr) {
                cerr << "\tRight end alignment is null" << endl;
            }
            else {
                cerr << "\tRight end final mapping does not match expected end: " << pb2json(*rightmost_mapping) << endl;
            }
            cerr << "\tRealigning right end to reverse complement" << endl;
#endif
            
            // TODO: what if we have an exact palindrome over a reversing
            // edge, and we keep getting the same alignment arbitrarily no
            // matter how we flip the sequence to align?
            aln_right = reverse_complement_alignment(right_subgraph.align(reverse_complement(right_tail), &aligner,
                                                                          true, false, 0, true, true, false, 0, max_span), node_length_function);
        }
        
        // Rescore the alignments as if N/N substitutions are matches so that we
        // can use the score as a proxy for how many matches there are.
        
        // Turn N/N subs into matches and score the alignments without end bonuses.
        align_ns(left_subgraph, aln_left);
        aln_left.set_score(aligner.score_contiguous_alignment(aln_left, true));
        align_ns(right_subgraph, aln_right);
        aln_right.set_score(aligner.score_contiguous_alignment(aln_right, true));
        
        
#ifdef debug
        cerr << "\t\tScores: " << aln_left.score() << "/" << left_tail.size() * aligner.match * min_score_factor
            << ", " << aln_right.score() << "/" << right_tail.size() * aligner.match * min_score_factor << endl;
#endif
        
        if (aln_left.score() > left_tail.size() * aligner.match * min_score_factor &&
            aln_right.score() > right_tail.size() * aligner.match * min_score_factor) {
        
            // Aligning the two tails suggests that the whole string might be in
            // the graph already.
            
            if (to_align.size() <= thin_alignment_cutoff && graph.length() < thin_alignment_cutoff) {
                // It's safe to try the tight banded alignment
        
                // We set this to true if we manage to find a valid alignment in the
                // narrow band.
                bool aligned_in_band;
                
                try {
                
#ifdef debug
                    cerr << "\tAttempt thin " << to_align.size() << " x " << graph.length() << " alignment" << endl;
#endif
                
                    // Throw it into the aligner with very restrictive banding to see if it's already basically present
                    aln = graph.align(to_align, &aligner,
                                      true, false, 0, false, false, true, large_alignment_band_padding, max_span);
                    Alignment aln2 = graph.align(reverse_complement(to_align), &aligner,
                                                 true, false, 0, false, false, true, large_alignment_band_padding, max_span);
                    if (aln2.score() > aln.score()) {
                        // The reverse alignment is better. But spit it back in the
                        // forward orientation.
                        aln = reverse_complement_alignment(aln2, node_length_function);
                    }
                    aligned_in_band = true;
                } catch(NoAlignmentInBandException ex) {
                    // If the aligner can't find any valid alignment in the restrictive
                    // band, we will need to knock together an alignment manually.
                    aligned_in_band = false;
#ifdef debug
                    cerr << "\t\tFailed." << endl;
#endif
                }

                // Rescore with Ns as matches again
                align_ns(left_subgraph, aln);
                aln.set_score(aligner.score_contiguous_alignment(aln, true));

#ifdef debug                
                if (aligned_in_band) {
                    cerr << "\tScore: " << aln.score() << "/" << (to_align.size() * aligner.match * min_score_factor) << endl;
                }
#endif
                
                if (aligned_in_band && aln.score() > to_align.size() * aligner.match * min_score_factor) {
                    // If we get a good score, use that alignment
#ifdef debug
                    cerr << "\tFound sufficiently good restricted banded alignment" << endl;
#endif
                    return aln;
                }
                
            } else if (to_align.size() < mapper_alignment_cutoff) {
                // It's safe to try the Mapper-based banded alignment
            
#ifdef debug
                cerr << "\tAttempt mapper-based " << to_align.size() << " x " << graph.length() << " alignment" << endl;
#endif
            
                // Otherwise, it's unsafe to try the tight banded alignment
                // (because our bands might get too big). Try a Mapper-based
                // fake-banded alignment, and return its alignment if it finds a
                // good one.
                
                // Generate an XG index
                xg::XG xg_index;
                xg_index.from_path_handle_graph(graph);

                // Generate a GCSA2 index
                gcsa::GCSA* gcsa_index = nullptr;
                gcsa::LCPArray* lcp_index = nullptr;
    
                if (edge_max) {
                    VG gcsa_graph = graph; // copy the graph
                    // remove complex components
                    algorithms::prune_complex_with_head_tail(gcsa_graph, kmer_size, edge_max);
                    if (subgraph_prune) {
                        algorithms::prune_short_subgraphs(gcsa_graph, subgraph_prune);
                    }
                        
                    // then index
#ifdef debug
                    cerr << "\tGCSA index size: " << gcsa_graph.length() << " bp" << endl;
#endif
                    build_gcsa_lcp(gcsa_graph, gcsa_index, lcp_index, kmer_size, doubling_steps);
                } else {
                    // if no complexity reduction is requested, just build the index
#ifdef debug
                    cerr << "\tGCSA index size: " << graph.length() << " bp" << endl;
#endif
                    build_gcsa_lcp(graph, gcsa_index, lcp_index, kmer_size, doubling_steps);
                }
                        
                // Make the Mapper
                Mapper mapper(&xg_index, gcsa_index, lcp_index);
                // Copy over alignment scores
                mapper.set_alignment_scores(aligner.match, aligner.mismatch, aligner.gap_open, aligner.gap_extension,
                                            aligner.full_length_bonus);
                
                // Map. Will invoke the banded aligner if the read is long, and
                // the normal index-based aligner otherwise.
                // Note: reverse complement is handled by the mapper.
                aln = mapper.align(to_align);
                
                // Clean up indexes
                delete lcp_index;
                delete gcsa_index;
                
                // Trace the path and complain somehow when it jumps.
                size_t discontinuities = 0;
                for (size_t i = 1; i < aln.path().mapping_size(); i++) {
                    // For every mapping after the first, look at where we came
                    // from
                    auto& last_position = aln.path().mapping(i-1).position();
                    // And how long we were
                    size_t consumed = mapping_from_length(aln.path().mapping(i-1));
                    // And where we are now
                    auto& next_position = aln.path().mapping(i).position();
                    
                    if (last_position.node_id() == next_position.node_id() &&
                        last_position.is_reverse() == next_position.is_reverse() &&
                        last_position.offset() + consumed == next_position.offset()) {
                        // The two mappings are on the same node with no discontinuity.
                        continue;
                    }
                    
                    if (last_position.offset() + consumed != graph.get_node(last_position.node_id())->sequence().size()) {
                        // We didn't use up all of the last node
                        discontinuities++;
                        continue;
                    }
                    
                    if (next_position.offset() != 0) {
                        // We aren't at the beginning of this node
                        discontinuities++;
                        continue;
                    }
                    
                    // If we get here we need to check for an edge
                    NodeSide last_side = NodeSide(last_position.node_id(), !last_position.is_reverse());
                    NodeSide next_side = NodeSide(next_position.node_id(), next_position.is_reverse());
                    
                    if (!graph.has_edge(last_side, next_side)) {
                        // There's no connecting edge between the nodes we hit
                        discontinuities++;
                    }
                }
                
                // Rescore with Ns as matches again
                align_ns(left_subgraph, aln);
                aln.set_score(aligner.score_contiguous_alignment(aln, true));
                
#ifdef debug
                cerr << "\tScore: " << aln.score() << "/" << (to_align.size() * aligner.match * min_score_factor)
                    << " with " << discontinuities << " breaks" << endl;
#endif
                
                if (aln.score() > to_align.size() * aligner.match * min_score_factor && discontinuities == 0) {
                    // This alignment looks good.
                    
                    return aln;
                }
                
                
            } else {
#ifdef debug
                cerr << "\tNo safe full alignment option available" << endl;
#endif
            }
        
        }
        
        // If we get here, we couldn't find a good banded alignment, or it looks
        // like the ends aren't present already, or we're just too big to try
        // banded aligning.
#ifdef debug
        cerr << "\tSplicing tail alignments" << endl;
#endif
        
        // Splice left and right tails together with any remaining sequence we didn't have
        
        // One problem we cah have is the fact that the alignments we are
        // splicing may overlap. This can happen e.g. if we see both of the TSD
        // copies from a mobile element insertion, and we align them both back
        // to the reference version.
        
        // It's going to be much more complicated to try and prevent this from
        // happening than it will be to just deal with the resulting cycles,
        // which can describe real homology anyway.
        
        // But we do have to think about overlap between the query sequences
        // themselves...
        
        // How much overlap is there between the two tails? May be negative.
        int64_t overlap = (int64_t) aln_left.sequence().size() +
            (int64_t) aln_right.sequence().size() - (int64_t) to_align.size();
        
        if (overlap >= 0) {
            // All of the string is accounted for in these two
            // alignments, and we can cut them and splice them.
            
            // Take half the overlap off each alignment and paste together
            aln = simplify(merge_alignments(strip_from_end(aln_left, overlap / 2),
                strip_from_start(aln_right, (overlap + 1) / 2)));
                
            // TODO: produce a better score!
            aln.set_score(aln_left.score() + aln_right.score());
            
#ifdef debug
            cerr << "\tSpliced overlapping end alignments" << endl;
#endif
            
        } else {
            // Not all of the string is accounted for in these two
            // alignments, so we will splice them together with any
            // remaining input sequence.
            
            string middle_sequence = to_align.substr(aln_left.sequence().size(), -overlap);
            
            // Make a big bogus alignment with an unplaced pure insert mapping.
            Alignment aln_middle;
            aln_middle.set_sequence(middle_sequence);
            auto* middle_mapping = aln_middle.mutable_path()->add_mapping();
            auto* middle_edit = middle_mapping->add_edit();
            middle_edit->set_sequence(middle_sequence);
            middle_edit->set_to_length(middle_sequence.size());
            
            // Paste everything together
            aln = simplify(merge_alignments(merge_alignments(aln_left, aln_middle), aln_right));
            
            // TODO: produce a better score!
            aln.set_score(aln_left.score() + aln_right.score());
            
#ifdef debug
            cerr << "\tSpliced disconnected end alignments" << endl;
#endif
            
        }
            
    }
    
    // TODO: check if we got alignments that didn't respect our specified
    // endpoints by one of the non-splicing-together alignment methods.
    
    return aln;

}

set<vector<int>> VariantAdder::get_unique_haplotypes(const vector<vcflib::Variant*>& variants, WindowedVcfBuffer* cache) const {
    set<vector<int>> haplotypes;
    
    if (variants.empty()) {
        // Nothing's there
        return haplotypes;
    }
    
    for (size_t sample_index = 0; sample_index < variants.front()->sampleNames.size(); sample_index++) {
        // For every sample
        auto& sample_name = variants.front()->sampleNames[sample_index];
        
        // Make its haplotype(s) on the region. We have a map from haplotype
        // number to actual vector. We'll tack stuff on the ends when they are
        // used, then throw out any that aren't full-length.
        map<size_t, vector<int>> sample_haplotypes;
        
        
        for (auto* variant : variants) {
            // Get the genotype for each sample
            const vector<int>* genotype;
            
            if (cache != nullptr) {
                // Use the cache provided by the buffer
                genotype = &cache->get_parsed_genotypes(variant).at(sample_index);
            } else {
                // Parse from the variant ourselves
                auto genotype_string = variant->getGenotype(sample_name);
            
                // Fake it being phased
                replace(genotype_string.begin(), genotype_string.end(), '/', '|');
                
                genotype = new vector<int>(vcflib::decomposePhasedGenotype(genotype_string));
            }
            
#ifdef debug
            cerr << "Genotype of " << sample_name << " at " << variant->position << ": ";
            for (auto& alt : *genotype) {
                cerr << alt << " ";
            }
            cerr << endl;
#endif
            
            for (size_t phase = 0; phase < genotype->size(); phase++) {
                // For each phase in the genotype
                
                // Get the allele number and ignore missing data
                int allele_index = (*genotype)[phase];
                if (allele_index == vcflib::NULL_ALLELE) {
                    allele_index = 0;
                }
                
                if (allele_index >= variant->alleles.size()) {
                    // This VCF has out-of-range alleles
                    cerr << "error:[vg::VariantAdder] variant " << variant->sequenceName << ":" << variant->position
                        << " has invalid allele index " << allele_index
                        << " but only " << variant->alt.size() << " alts" << endl;
                    exit(1);
                }
                
                if (skip_structural_duplications && 
                    variant->info.count("SVTYPE") &&
                    variant->info.at("SVTYPE").size() == 1 && 
                    (variant->info.at("SVTYPE").at(0) == "CNV" || variant->info.at("SVTYPE").at(0) == "DUP") &&
                    variant->alleles.at(allele_index).size() > variant->alleles.at(0).size()) {
                    
                    // This variant is a duplication (or a CNV) structural
                    // variant, and this is the duplicated version, and we want
                    // to skip those.
                    
                    // We only want to skip the duplication alts; deletion alts
                    // of CNVs are passed through.
                    
                    // Don't add to this haplotype. It will get filtered out for
                    // not being full length, and we won't consider the
                    // duplication.
                    continue;
                    
                }
                
                // Stick each allele number at the end of its appropriate phase
                sample_haplotypes[phase].push_back(allele_index);
            }
            
            if (cache == nullptr) {
                // We're responsible for this vector
                delete genotype;
            }
        }
        
        for (auto& kv : sample_haplotypes) {
            auto& haplotype = kv.second;
            // For every haplotype in this sample
            if (haplotype.size() != variants.size()) {
                // If it's not the full length, it means some variants don't
                // have it. Skip.
                continue;
            }
            
            // Otherwise, add it to the set of observed haplotypes
            haplotypes.insert(haplotype);
        }
    }
    
    // After processing all the samples, return the unique haplotypes
    return haplotypes;
    
}

string VariantAdder::haplotype_to_string(const vector<int>& haplotype, const vector<vcflib::Variant*>& variants) {
    // We'll fill this in with variants and separating sequences.
    stringstream result;
    
    // These lists need to be in 1 to 1 correspondence
    assert(haplotype.size() == variants.size());
    
    if (variants.empty()) {
        // No variants means no string representation.
        return "";
    }
    
    // Do the first variant
    result << variants.front()->alleles.at(haplotype.at(0));
    
    for (size_t i = 1; i < variants.size(); i++) {
        // For each subsequent variant
        auto* variant = variants.at(i);
        auto* last_variant = variants.at(i - 1);
        
        // Do the intervening sequence.
        // Where does that sequence start?
        size_t sep_start = last_variant->zeroBasedPosition() + last_variant->ref.size();
        // And how long does it run?
        size_t sep_length = variant->zeroBasedPosition() - sep_start;
        
        // Find the sequence to pull from
        auto& ref = sync.get_path_sequence(vcf_to_fasta(variant->sequenceName));
        
        // Pull out the separator sequence and tack it on.
        result << ref.substr(sep_start, sep_length);

        // Then put the appropriate allele of this variant.
        result << variant->alleles.at(haplotype.at(i));
        
        if (variant->alleles.at(0) != ref.substr(sep_start + sep_length, variant->alleles.at(0).size())) {
            // Complain if the variant reference doesn't match the real reference.
            // TODO: should avoid doing this for every single haplotype...
            throw runtime_error("Variant reference does not match actual reference at " +
                variant->sequenceName + ":" + to_string(variant->position));
        }
    }
    
    return result.str();
}

size_t VariantAdder::get_radius(const vcflib::Variant& variant) {
    // How long is the longest alt?
    size_t longest_alt_length = variant.ref.size();
    for (auto& alt : variant.alt) {
        // Take the length of the longest alt you find
        longest_alt_length = max(longest_alt_length, alt.size());
    }
    
    // Report half its length, and don't round down.
    return (longest_alt_length + 1) / 2;
}


size_t VariantAdder::get_center(const vcflib::Variant& variant) {
    // Where is the end of the variant in the reference?
    size_t path_last = variant.zeroBasedPosition() + variant.ref.size() - 1;
    
    // Where is the center of the variant in the reference?
    return (variant.zeroBasedPosition() + path_last) / 2;
}


pair<size_t, size_t> VariantAdder::get_center_and_radius(const vector<vcflib::Variant*>& variants) {

    // We keep track of the leftmost and rightmost positions we would need to
    // cover, which may be negative on the left.
    int64_t leftmost = numeric_limits<int64_t>::max();
    int64_t rightmost = 0;

    for (auto* v : variants) {
        // For every variant
        auto& variant = *v;
        
        // Work out its center (guaranteed positive)
        int64_t center = get_center(variant);
        // And its radius
        int64_t radius = get_radius(variant);
        
        // Expand the range of the block if needed
        leftmost = min(leftmost, center - radius);
        rightmost = max(rightmost, center + radius);
    }
    
    // Calculate the center between the two ends, and the radius needed to hit
    // both ends.
    size_t overall_center = (leftmost + rightmost) / 2;
    size_t overall_radius = (rightmost - leftmost + 1) / 2;
    
    return make_pair(overall_center, overall_radius);

}

vector<vcflib::Variant*> VariantAdder::filter_local_variants(const vector<vcflib::Variant*>& before,
    vcflib::Variant* variant, const vector<vcflib::Variant*>& after) const {

    // This is the filter we apply
    auto filter = [&](vcflib::Variant* v) {
        // Keep a variant if it isn't too big
        return get_radius(*v) <= max_context_radius;
    };

    // Make the list of all the local variants in one vector
    vector<vcflib::Variant*> local_variants;
    
    // Keep the nearby variants if they pass the test
    copy_if(before.begin(), before.end(), back_inserter(local_variants), filter);
    // And the main variant always
    local_variants.push_back(variant);
    copy_if(after.begin(), after.end(), back_inserter(local_variants), filter);

#ifdef debug
        cerr << "Local variants: ";
        for (auto* v : local_variants) {
            cerr << vcf_to_fasta(v->sequenceName) << ":" << v->position << " ";
        }
        cerr << endl;
#endif

    return local_variants;
}

}





















