#include "variant_adder.hpp"

namespace vg {

using namespace std;

VariantAdder::VariantAdder(VG& graph) : graph(graph), sync(graph) {
    graph.paths.for_each_name([&](const string& name) {
        // Save the names of all the graph paths, so we don't need to lock the
        // graph to check them.
        path_names.insert(name);
    });
    
    // Show progress if the graph does.
    show_progress = graph.show_progress;
}

void VariantAdder::add_variants(vcflib::VariantCallFile* vcf) {
    
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
        auto& variant_path_offset = variant->position; // Already made 0-based by the buffer
        
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
        vector<vcflib::Variant*> local_variants{before};
        local_variants.push_back(variant);
        copy(after.begin(), after.end(), back_inserter(local_variants));
        
#ifdef debug
        cerr << "Local variants: ";
        for (auto* v : local_variants) {
            cerr << vcf_to_fasta(v->sequenceName) << ":" << v->position << " ";
        }
        cerr << endl;
#endif
        
        // Where does the group of nearby variants start?
        size_t group_start = local_variants.front()->position;
        // And where does it end (exclusive)? This is the latest ending point of any variant in the group...
        size_t group_end = local_variants.back()->position + local_variants.back()->ref.size();
        
        // We need to make sure we also grab this much extra graph context,
        // since we count 2 radiuses + flank out from the ends of the group.
        size_t group_width = group_end - group_start;
        
        // Find the center and radius of the group of variants, so we know what graph part to grab.
        size_t overall_center;
        size_t overall_radius;
        tie(overall_center, overall_radius) = get_center_and_radius(local_variants);
        
        // Get the leading and trailing ref sequence on either side of this group of variants (to pin the outside variants down).

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
        
        // Round bounds to node start and endpoints.
        sync.with_path_index(variant_path_name, [&](const PathIndex& index) {
            tie(left_context_start, right_context_past_end) = index.round_outward(left_context_start, right_context_past_end);
        });
        
#ifdef debug
        cerr << "New context bounds: " << left_context_start << " - " << right_context_past_end << endl;
#endif
        
        // Recalculate context lengths
        left_context_length = group_start - left_context_start;
        right_context_length = right_context_past_end - group_end;
        
        // Make sure we pull out out to the ends of the contexts
        overall_radius = max(overall_radius, max(overall_center - left_context_start, right_context_past_end - overall_center));
        
        // Get actual context strings
        string left_context = path_sequence.substr(group_start - left_context_length, left_context_length);
        string right_context = path_sequence.substr(group_end, right_context_length);
        
        // Get the unique haplotypes
        auto haplotypes = get_unique_haplotypes(local_variants, &buffer);
        
        // Track the total bp of haplotypes
        size_t total_haplotype_bases = 0;
        
        // Track the total graph size for the alignments
        size_t total_graph_bases = 0;
        
#ifdef debug
        cerr << "Have " << haplotypes.size() << " haplotypes for variant "
            << variant->sequenceName << ":" << variant->position << endl;
#endif
        
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
            
            // Make its combined string
            stringstream to_align;
            to_align << left_context << haplotype_to_string(haplotype, local_variants) << right_context;
            
#ifdef debug
            cerr << "Align " << to_align.str() << endl;
#endif

            // Count all the bases
            total_haplotype_bases += to_align.str().size();
            
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
            
            // Calculate the expected min score, for a giant gap (i.e. this variant is an SV indel)
            int64_t min_length = min(right_context_past_end - left_context_start, to_align.str().size());
            int64_t max_length = max(right_context_past_end - left_context_start, to_align.str().size());
            // TODO: assumes a gap open of 6 and match and gap extend values of 1
            int64_t expected_score = min_length - 6 - (max_length - 1 - min_length);
                
            // Record the size of graph we're aligning to in bases
            total_graph_bases += lock.get_subgraph().length();
            
            // Work out how far we would have to unroll the graph to account for
            // a giant deletion. We also want to account for alts that may
            // alrteady be in the graph and need unrolling for a long insert.
            size_t max_span = max(right_context_past_end - left_context_start, to_align.str().size());
            
            // Do the alignment in both orientations
            Alignment aln;
            Alignment aln2;
            
            // Align in the forward orientation using banded global aligner, unrolling for large deletions.
            aln = lock.get_subgraph().align(to_align.str(), 0, false, false, 0, true, max_span);
            // Align in the reverse orientation using banded global aligner, unrolling for large deletions.
            // TODO: figure out which way our reference path goes through our subgraph and do half the work.
            // Note that if we have reversing edges and a lot of unrolling, we might get the same alignment either way.
            aln2 = lock.get_subgraph().align(reverse_complement(to_align.str()), 0, false, false, 0, true, max_span);
            
            // Note that the banded global aligner doesn't fill in identity.
            
#ifdef debug
            cerr << "Scores: " << aln.score() << " fwd vs. " << aln2.score() << " rev" << endl;
#endif
                
            if (aln2.score() > aln.score()) {
                // This is the better alignment
                swap(aln, aln2);
            }

#ifdef debug
            cerr << "Subgraph: " << pb2json(lock.get_subgraph().graph) << endl;            
            cerr << "Alignment: " << pb2json(aln) << endl;
#endif

            if (local_variants.size() == 1) {
                // We only have one variant here, so we ought to have at least the expected giant gap score
                assert(aln.score() >= expected_score);
            }
            
            // We shouldn't have dangling ends, really, but it's possible for
            // inserts that have copies already in the graph to end up producing
            // alignments just as good as the alignment we wanted that have
            // their gaps pushed to one end or the other, and we need to
            // tolerate them and make their insertions.
            
            // We know the aligner left-shifts the gaps for inserts, so make
            // sure that we at least *end* with a match.
            assert(aln.path().mapping_size() > 0);
            auto& last_mapping = aln.path().mapping(aln.path().mapping_size() - 1);
            assert(last_mapping.edit_size() > 0);
            auto& last_edit = last_mapping.edit(last_mapping.edit_size() - 1);
            assert(edit_is_match(last_edit));
            
            // Find the first edit and get its oriented node
            auto& first_mapping = aln.path().mapping(0);
            assert(first_mapping.edit_size() > 0);
            auto& first_edit = first_mapping.edit(0);
            
            // Construct the NodeSide on the left of the graph in the orientation the graph is aligned to.
            NodeSide left_of_alignment(first_mapping.position().node_id(), first_mapping.position().is_reverse());
            
            // Get all the NodeSides connected to it in the periphery of the
            // graph we extracted.
            set<NodeSide> connected = lock.get_peripheral_attachments(left_of_alignment);
            
            cerr << "Alignment starts at " << left_of_alignment << " which connects to ";
            for (auto& c : connected) {
                cerr << c << ", ";
            }
            cerr << endl;
            
            // Make this path's edits to the original graph. We don't need to do
            // anything with the translations. Handle insertions on the very
            // left by attaching them to whatever is attached to our leading
            // node.
            lock.apply_edit(aln.path(), connected);

        }
        
        if (variants_processed++ % 1000 == 0 || true) {
            cerr << "Variant " << variants_processed << ": " << haplotypes.size() << " haplotypes at "
                << variant->sequenceName << ":" << variant->position << ": "
                << (total_haplotype_bases / haplotypes.size()) << " bp vs. "
                << (total_graph_bases / haplotypes.size()) << " bp haplotypes vs. graphs average" << endl;
        }
        
    }

    // Clean up after the last contig.
    destroy_progress();
    
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
                    //cerr << "error:[vg::VariantAdder] variant " << variant->sequenceName << ":" << variant->position
                    //    << " has invalid allele index " << allele_index
                    //    << " but only " << variant->alt.size() << " alts" << endl;
                    // TODO: right now we skip them as a hack
                    allele_index = 0;
                    //exit(1);
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
        size_t sep_start = last_variant->position + last_variant->ref.size();
        // And how long does it run?
        size_t sep_length = variant->position - sep_start;
        
        // Pull out the separator sequence and tack it on.
        result << sync.get_path_sequence(vcf_to_fasta(variant->sequenceName)).substr(sep_start, sep_length);

        // Then put the appropriate allele of this variant.
        result << variant->alleles.at(haplotype.at(i));
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
    size_t path_last = variant.position + variant.ref.size() - 1;
    
    // Where is the center of the variant in the reference?
    return (variant.position + path_last) / 2;
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

}





















