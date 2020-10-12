#include "readfilter.hpp"
#include "IntervalTree.h"
#include "annotation.hpp"
#include "vg/io/alignment_emitter.hpp"
#include <vg/io/stream.hpp>

#include <fstream>
#include <sstream>

#include <htslib/khash.h>

namespace vg {

using namespace std;
using namespace vg::io;

bool ReadFilter::trim_ambiguous_ends(Alignment& alignment, int k) {
    assert(graph != nullptr);

    // Define a way to get node length, for flipping alignments
    function<int64_t(id_t)> get_node_length = [&](id_t node) {
        return graph->get_length(graph->get_handle(node));
    };

    // Because we need to flip the alignment, make sure it is new-style and
    // doesn't have any Mappings with no Edits.
    for(size_t i = 0; i < alignment.path().mapping_size(); i++) {
        if(alignment.path().mapping(i).edit_size() == 0) {
            // Complain!
            throw runtime_error("Found mapping with no edits in " + pb2json(alignment));
        }
    }

    // TODO: we're going to flip the alignment twice! This is a waste of time!
    // Define some kind of oriented view or something, or just two duplicated
    // trimming functions, so we can just trim once without flipping.

    // Trim the end
    bool end_changed = trim_ambiguous_end(alignment, k);
    // Flip and trim the start
    
    Alignment flipped = reverse_complement_alignment(alignment, get_node_length);
    
    if(trim_ambiguous_end(flipped, k)) {
        // The start needed trimming
        
        // Flip the trimmed flipped alignment back
        alignment = reverse_complement_alignment(flipped, get_node_length);
        // We definitely changed something
        return true;
    }
    
    // We maybe changed something
    return end_changed;
    
}

bool ReadFilter::trim_ambiguous_end(Alignment& alignment, int k) {
    // What mapping in the alignment is the leftmost one starting in the last k
    // bases? (Except when that would be the first mapping, we use the second.)
    // Start out with it set to the past-the-end value.
    size_t trim_start_mapping = alignment.path().mapping_size();
    
    // How many real non-softclip bases have we seen reading in from the end of
    // the read?
    size_t real_base_count = 0;
    // How many softclip bases have we seen in from the end of the read?
    size_t softclip_base_count = 0;
    for(size_t i = alignment.path().mapping_size() - 1; i != -1 && i != 0; i--) {
        // Scan in from the end of the read.
        
        auto* mapping = alignment.mutable_path()->mutable_mapping(i);
        
        // We should always have edits in our mappings.
        assert(mapping->edit_size() > 0);
        
        for(int j = mapping->edit_size() - 1; j != -1; j--) {
            // Visit every edit in the mapping
            auto& edit = mapping->edit(j);
            
            
            if(real_base_count == 0 && edit.from_length() == 0) {
                // This is a trailing insert. Put it as a softclip
                softclip_base_count += edit.to_length();
            } else {
                // This is some other kind of thing. Record it as real bases.
                real_base_count += edit.to_length();
            }
        }
        
        if(real_base_count <= k) {
            // This mapping starts fewer than k non-softclipped alignment
            // bases from the end of the read.
            trim_start_mapping = i;
        } else {
            // This mapping starts more than k in from the end. So the
            // previous one, if we had one, must be the right one.
            break;
        }
    }
    
    if(trim_start_mapping == alignment.path().mapping_size()) {
        // No mapping was found that starts within the last k non-softclipped
        // bases. So there's nothing to do.
        return false;
    }
    
    if(real_base_count == 0) {
        // We have an anchoring mapping, but all the mappings we could trim are
        // softclips, so there's no point. TODO: will we ever get softclips
        // placed as the only thing on a node?
        return false;
    }
    
    // Which is the last assumed-non-ambiguous mapping from which we can anchor
    // our search?
    size_t root_mapping = trim_start_mapping - 1;
    
    // What's the sequence, including that root node, that we are looking for?
    // We need the sequence of the nodes, rather than the read's sequence,
    // because you can still be ambiguous even if you have a SNP on top of the
    // ambiguous thing.
    
    // We need to ignore all the offsets and from_lengths, except for the from
    // length on the last node to let us know if we end early. It's sort of
    // nonsense to have offsets and non-full from_lengths on internal mappings,
    // and everything is easiest if we use the full length sequence of the root
    // node.
    stringstream target_sequence_stream;
    for(size_t i = root_mapping; i < alignment.path().mapping_size(); i++) {
        // Collect the appropriately oriented from sequence from each mapping
        auto& mapping = alignment.path().mapping(i);
        handle_t handle = graph->get_handle(mapping.position().node_id(),
                                            mapping.position().is_reverse());
        string sequence = graph->get_sequence(handle);
        
        if(i == root_mapping) {
            // Use the full length of the node and ignore any offset
            target_sequence_stream << sequence;
        } else {
            // Use the offset plus the total from_length of all the
            // edits (in case we're the last node and ending early). We made
            // sure all non-root nodes had edits earlier.
            
            size_t from_length = mapping.position().offset();
            for(size_t j = 0; j < mapping.edit_size(); j++) {
                from_length += mapping.edit(j).from_length();
            }
            
            // Put in the sequence that the mapping visits
            target_sequence_stream << sequence.substr(0, from_length);
        }
    }
    string target_sequence = target_sequence_stream.str();
        
#ifdef debug
    #pragma omp critical(cerr)
    cerr << "Need to look for " << target_sequence << " right of mapping " << root_mapping << endl;
#endif
    
    // We're not going to recurse hundreds of nodes deep, so we can use the real
    // stack and a real recursive function.
    
    // Do the DFS into the given node, after already having matched the given
    // number of bases of the target sequence. See if you can match any more
    // bases of the target sequence.
    
    // Return the total number of leaves in all subtrees that match the full
    // target sequence, and the depth in bases of the shallowest point at which
    // multiple subtrees with full lenght matches are unified.

    // We keep a maximum number of visited nodes here, just to prevent recursion
    // from going on forever in worst-case-type graphs
    size_t dfs_visit_count = 0;
    function<pair<size_t, size_t>(const handle_t&, size_t)> do_dfs =
        [&](const handle_t& handle, size_t matched) -> pair<size_t, size_t> {

        ++dfs_visit_count;
      
        // Grab the node sequence and match more of the target sequence.
        string node_sequence = graph->get_sequence(handle);
        
#ifdef debug
        #pragma omp critical(cerr)
        cerr << "Node " << graph->get_id(handle) <<  " " << (graph->get_is_reverse(handle) ? "rev" : "fwd") << ": "
            << node_sequence << " at offset " << matched << " in " << target_sequence << endl;
#endif
        
        // Now count up the new matches between this node and the target sequence.
        size_t new_matches;
        for(
            // Start with no matches
            new_matches = 0;
            // Keep going as long as we're inside both strings and haven't had a mismatch
            new_matches < node_sequence.size() && 
            matched + new_matches < target_sequence.size() && 
            node_sequence[new_matches] == target_sequence[matched + new_matches];
            // Count up all the matches we find
            new_matches++
        );

        if(matched + new_matches == target_sequence.size()) {
            // We found a tail end of a complete match of the target sequence
            // on this node.
            
#ifdef debug
            #pragma omp critical(cerr)
            cerr << "Node " << node_id << " is a matching leaf" << endl;
#endif
            
            // Return one match and unification at full length (i.e. nothing can
            // be discarded).
            return make_pair(1, target_sequence.size());
        }
        
        if(new_matches < node_sequence.size()) {
            // We didn't make it all the way through this node, nor did we
            // finish the target sequence; there's a mismatch between the node
            // and the target sequence.
            
#ifdef debug
            #pragma omp critical(cerr)
            cerr << "Node " << node_id << " has a mismatch" << endl;
#endif
            
            // If we mismatch, return 0 matches and unification at full length.
            return make_pair(0, target_sequence.size());
        }        
        
        // If we get through the whole node sequence without mismatching or
        // running out of target sequence, keep going.
        
#ifdef debug
        #pragma omp critical(cerr)
        cerr << "Node " << graph->get_id(handle) << " has " << new_matches << " internal new matches" << endl;
#endif
        
        // We're going to call all the children and collect the results, and
        // then aggregate them. It might be slightly faster to aggregate while
        // calling, but that might be less clear.
        vector<pair<size_t, size_t>> child_results;
        
        graph->follow_edges(handle, false, [&](const handle_t& next) {
            if (dfs_visit_count < defray_count) {
                child_results.push_back(do_dfs(next, matched + node_sequence.size()));
                return true;
            }
            else {
#ifdef debug
#pragma omp critical(cerr)
                cerr << "Aborting read filter DFS at node " << graph->get_id(next) << " after " << dfs_visit_count << " visited" << endl;
#endif
                return false;
            }
        });
        
        // Sum up the total leaf matches, which will be our leaf match count.
        size_t total_leaf_matches = 0;
        // If we don't find multiple children with leaf matches, report
        // unification at the min unification depth of any subtree (and there
        // will only be one that isn't at full length).
        size_t children_with_leaf_matches = 0;
        size_t unification_depth = target_sequence.size();
        
        for(auto& result : child_results) {
            total_leaf_matches += result.first;
            if(result.first > 0) {
                children_with_leaf_matches++;
            }
            unification_depth = min(unification_depth, result.second);
        }
        if(children_with_leaf_matches > 1) {
            // If multiple children have nonzero leaf match counts, report
            // unification at the end of this node.
            unification_depth = matched + node_sequence.size();
        }
        
        return make_pair(total_leaf_matches, unification_depth);
    };
    
    // Search from the root mapping's node looking right in its orientation in
    // the mapping
    const auto& root_pos = alignment.path().mapping(root_mapping).position();
    auto result = do_dfs(graph->get_handle(root_pos.node_id(), root_pos.is_reverse()), 0);
    
#ifdef debug
    #pragma omp critical(cerr)
    cerr << "Found " << result.first << " matching leaves with closest unification at " << result.second << endl;
#endif
    
    // We keep this much of the target sequence.
    size_t target_sequence_to_keep = result.second;
    
    if(target_sequence_to_keep == target_sequence.size()) {
        // Nothing to trim!
        return false;
    }
    
    // Figure out how many mappings we need to keep from the root in order to
    // get that much sequence; we know it is either full length or at a mapping
    // boundary. We handle the root special because it's always full length and
    // we have to cut after its end.
    size_t kept_sequence_accounted_for = graph->get_length(graph->get_handle(alignment.path().mapping(root_mapping).position().node_id()));
    size_t first_mapping_to_drop;
    for(first_mapping_to_drop = root_mapping + 1;
        first_mapping_to_drop < alignment.path().mapping_size();
        first_mapping_to_drop++) {
        // Consider starting dropping at each mapping after the root.
        if(kept_sequence_accounted_for == target_sequence_to_keep) {
            // OK this mapping really is the first one to drop.
            break;
        } else {
            // Keep going. Account for the sequence from this mapping.
            auto& mapping = alignment.path().mapping(first_mapping_to_drop);
            
            // We know it's not the root mapping, and it can't be the non-full-
            // length end mapping (because we would have kept the full length
            // target sequence and not had to cut). So assume full node is used.
            kept_sequence_accounted_for += graph->get_length(graph->get_handle(mapping.position().node_id()));
        }
    }
    
    // OK we know the first mapping to drop. We need to work out the to_size,
    // including all softclips, from there to the end, so we know how much to
    // trim off of the sequence and quality.
    size_t to_length_to_remove = 0;
    for(size_t i = first_mapping_to_drop; i < alignment.path().mapping_size(); i++) {
        // Go through all the mappings
        auto& mapping = alignment.path().mapping(i);
        for(size_t j = 0; j < mapping.edit_size(); j++) {
            // Add up the to_length of all the edits
            to_length_to_remove += mapping.edit(j).to_length();
        }
    }
    
#ifdef debug
    cerr << "Want to trim " << alignment.sequence().size() << " bp sequence and " << alignment.quality().size()
        << " quality values to remove " << to_length_to_remove << endl;
#endif
        
    // Make sure we have at least enough to trim.
    // Note that we allow the entire alignment to be trimmed away!
    assert(alignment.sequence().size() >= to_length_to_remove);
    assert(alignment.quality().empty() || alignment.quality().size() >= to_length_to_remove);
    // And that we made sence before trimming
    assert(alignment.quality().empty() || alignment.quality().size() == alignment.sequence().size());
   
    // Trim sequence
    alignment.set_sequence(alignment.sequence().substr(0, alignment.sequence().size() - to_length_to_remove));
    
    // Trim quality
    if(!alignment.quality().empty()) {
        alignment.set_quality(alignment.quality().substr(0, alignment.quality().size() - to_length_to_remove));
    }
    
    // Now we can discard the extra mappings
    size_t to_delete = alignment.path().mapping_size() - first_mapping_to_drop;
    alignment.mutable_path()->mutable_mapping()->DeleteSubrange(first_mapping_to_drop, to_delete);
    
    // Now the alignment is fixed!
    return true;
}

// quick and dirty filter to see if removing reads that can slip around
// and still map perfectly helps vg call.  returns true if at either
// end of read sequence, at least k bases are repetitive, checking repeats
// of up to size 2k
bool ReadFilter::has_repeat(Alignment& aln, int k) {
    if (k == 0) {
        return false;
    }
    const string& s = aln.sequence();
    for (int i = 1; i <= 2 * k; ++i) {
        int covered = 0;
        bool ffound = true;
        bool bfound = true;
        for (int j = 1; (ffound || bfound) && (j + 1) * i < s.length(); ++j) {
            ffound = ffound && s.substr(0, i) == s.substr(j * i, i);
            bfound = bfound && s.substr(s.length() - i, i) == s.substr(s.length() - i - j * i, i);
            if (ffound || bfound) {
                covered += i;
            }
        }
        if (covered >= k) {
            return true;
        }
    }
    return false;
}

bool ReadFilter::is_split(Alignment& alignment) {
    if(graph == nullptr) {
        // Can't tell if the read is split.
        throw runtime_error("HandleGraph (e.g. XG) required to check for split reads");
    }
    
    handle_t prev;
    for(size_t i = 0; i + 1 < alignment.path().mapping_size(); i++) {
        if (i == 0) {
            const auto& pos = alignment.path().mapping(i).position();
            prev = graph->get_handle(pos.node_id(), pos.is_reverse());
        }
        const auto& pos = alignment.path().mapping(i + 1).position();
        handle_t here = graph->get_handle(pos.node_id(), pos.is_reverse());
        
        // Can we find the same articulation of the edge as the alignment uses
        
        if(!graph->has_edge(prev, here)) {
            // We found a skip!
            if(verbose) {
                cerr << "Warning: read " << alignment.name() << " has an unknown edge "
                << graph->get_id(prev) << (graph->get_is_reverse(prev) ? "-" : "+") << " -> " << graph->get_id(here) << (graph->get_is_reverse(here) ? "-" : "+")
                    << ". Removing!" << endl;
            }
            return true;
        }
        
        prev = here;
    }
    
    // No wandering jumps between nodes found
    return false;
}


bool ReadFilter::sample_read(const Alignment& aln) {
    // Decide if the alignment is paired.
    // It is paired if fragment_next or fragment_prev point to something.
    bool is_paired = (!aln.fragment_prev().name().empty() || aln.fragment_prev().path().mapping_size() != 0 ||
        !aln.fragment_next().name().empty() || aln.fragment_next().path().mapping_size() != 0);

    // Compute the QNAME that samtools would use
    string qname;
    if (is_paired) {
        // Strip pair end identifiers like _1 or /2 that vg uses at the end of the name.    
        qname = regex_replace(aln.name(), regex("[/_][12]$"), "");
    } else {
        // Any _1 in the name is part of the actual read name.
        qname = aln.name();
    }
    
    // Now treat it as samtools would.
    // See https://github.com/samtools/samtools/blob/60138c42cf04c5c473dc151f3b9ca7530286fb1b/sam_view.c#L101-L104
    
    // Hash that with __ac_X31_hash_string from htslib and XOR against the seed mask
    auto masked_hash = __ac_X31_hash_string(qname.c_str()) ^ downsample_seed_mask;
    
    // Hash that again with __ac_Wang_hash from htslib, apparently to mix the bits.
    uint32_t mixed_hash = __ac_Wang_hash(masked_hash);
    
    // Take the low 24 bits and compute a double from 0 to 1
    const int32_t LOW_24_BITS = 0xffffff;
    double sample = ((double)(mixed_hash & LOW_24_BITS)) / (LOW_24_BITS + 1);
    
    // If the result is >= the portion to downsample to, discard the read.
    // Otherwise, keep it.
    return (sample < downsample_probability);
}

ReadFilter::Counts ReadFilter::filter_alignment(Alignment& aln) {
    Counts counts;
    
    double score = (double)aln.score();
    double denom = aln.sequence().length();
    // toggle substitution score
    if (sub_score == true) {
        // hack in ident to replace old counting logic.
        score = aln.identity() * aln.sequence().length();
        assert(score <= denom);
    } else if (rescore == true) {
        // We need to recalculate the score with the base aligner always
        const static Aligner unadjusted;
        GSSWAligner* aligner = (GSSWAligner*)&unadjusted;
            
        // Rescore and assign the score
        aln.set_score(aligner->score_contiguous_alignment(aln));
        // Also use the score
        score = aln.score();
    }

    // toggle absolute or fractional score
    if (frac_score == true) {
        if (denom > 0.) {
            score /= denom;
        }
        else {
            assert(score == 0.);
        }
    }
        
    ++counts.counts[Counts::FilterName::read];
    bool keep = true;
    // filter (current) alignment
    if (!name_prefixes.empty()) {
        // Make sure we match at least one name prefix
            
        bool found = false;
            
        // Do a binary search for the closest prefix and see if all of any prefix exists.
        // We assume the prefixes are sorted.
        size_t left_bound = 0;
        size_t left_match = 0;
        while (left_match < name_prefixes[left_bound].size() &&
               left_match < aln.name().size() &&
               name_prefixes[left_bound][left_match] == aln.name()[left_match]) {
            // Scan all the matches at the start
            left_match++;
        }
            
        size_t right_bound = name_prefixes.size() - 1;
        size_t right_match = 0;
        while (right_match < name_prefixes[right_bound].size() &&
               right_match < aln.name().size() &&
               name_prefixes[right_bound][right_match] == aln.name()[right_match]) {
            // Scan all the matches at the end
            right_match++;
        }
            
        if (left_match == name_prefixes[left_bound].size() || right_match == name_prefixes[right_bound].size()) {
            // We found a match already
            found = true;
        } else {
            while (left_bound + 1 < right_bound) {
                // Until we run out of unexamined prefixes, do binary search
                size_t center = (left_bound + right_bound) / 2;
                // No need to re-check any common prefix
                size_t center_match = min(left_match, right_match);
                    
                while (center_match < name_prefixes[center].size() &&
                       center_match < aln.name().size() &&
                       name_prefixes[center][center_match] == aln.name()[center_match]) {
                    // Scan all the matches here
                    center_match++;
                }
                    
                if (center_match == name_prefixes[center].size()) {
                    // We found a hit!
                    found = true;
                    break;
                }
                
                if (center_match == aln.name().size() ||
                    name_prefixes[center][center_match] > aln.name()[center_match]) {
                    // The match, if it exists, must be before us
                    right_bound = center;
                    right_match = center_match;
                }
                else {
                    // The match, if it exists, must be after us.
                    left_bound = center;
                    left_match = center_match;
                }
            }
        }
            
        if (!found) {
            // There are prefixes and we don't match any, so drop the read.
            ++counts.counts[Counts::FilterName::wrong_name];
            keep = false;
        }
    }
    if ((keep || verbose) && !excluded_refpos_contigs.empty() && aln.refpos_size() != 0) {
        // We have refpos exclusion filters and a refpos is set.
        // We need to bang every refpos anme against every filter.
            
        bool found_match = false;
        for (auto& expression : excluded_refpos_contigs) {
            for (auto& refpos : aln.refpos()) {
                if (regex_search(refpos.name(), expression)) {
                    // We don't want this read because of this match
                    found_match = true;
                    break;
                }
            }
            if (found_match) {
                break;
            }
        }
            
        if (found_match) {
            ++counts.counts[Counts::FilterName::wrong_refpos];
            keep = false;    
        }
    }
    if ((keep || verbose) && !excluded_features.empty()) {
        // Get all the feature tags on the read
        vector<string> features(get_annotation<vector<string>>(aln, "features"));
        
        for (auto& feature : features) {
            if (excluded_features.count(feature)) {
                // If the read has any banned features, fail it.
                ++counts.counts[Counts::FilterName::excluded_feature];
                keep = false;
                break;
            }
        }
    }
    if ((keep || verbose) && (!aln.is_secondary() && score < min_primary)) {
        ++counts.counts[Counts::FilterName::min_score];
        keep = false;
    }
    if ((keep || verbose) && (aln.is_secondary() && score < min_secondary)) {
        ++counts.counts[Counts::FilterName::min_sec_score];
        keep = false;
    }
    if (keep || verbose) {
        // compute overhang
        int overhang = 0;
        if (aln.path().mapping_size() > 0) {
            const auto& left_mapping = aln.path().mapping(0);
            if (left_mapping.edit_size() > 0) {
                overhang = left_mapping.edit(0).to_length() - left_mapping.edit(0).from_length();
            }
            const auto& right_mapping = aln.path().mapping(aln.path().mapping_size() - 1);
            if (right_mapping.edit_size() > 0) {
                const auto& edit = right_mapping.edit(right_mapping.edit_size() - 1);
                overhang = max(overhang, edit.to_length() - edit.from_length());
            }
        } else {
            overhang = aln.sequence().length();
        }
        if (overhang > max_overhang) {
            ++counts.counts[Counts::FilterName::max_overhang];
            keep = false;
        }
    }        
    if (keep || verbose) {
        // compute end matches.
        int end_matches = 0;
        // from the left
        for (int i = 0; i < aln.path().mapping_size() && end_matches < min_end_matches; ++i) {
            for (int j = 0; j < aln.path().mapping(i).edit_size() && end_matches < min_end_matches; ++j) {
                const Edit& edit = aln.path().mapping(i).edit(j);
                if (edit.from_length() == edit.to_length() && edit.sequence().empty()) {
                    end_matches += edit.to_length();
                } else {
                    i = aln.path().mapping_size();
                    break;
                }
            }
        }
        if (end_matches >= min_end_matches) {
            end_matches = 0;
            // from the right
            for (int i = aln.path().mapping_size() - 1; i >= 0 && end_matches < min_end_matches; --i) {
                for (int j = aln.path().mapping(i).edit_size() - 1; j >= 0 && end_matches < min_end_matches; --j) {
                    const Edit& edit = aln.path().mapping(i).edit(j);
                    if (edit.from_length() == edit.to_length() && edit.sequence().empty()) {
                        end_matches += edit.to_length();
                    } else {
                        i = -1;
                        break;
                    }
                }
            }
        }                        
        if (end_matches < min_end_matches) {
            ++counts.counts[Counts::FilterName::min_end_matches];
            keep = false;
        }
    }
    if ((keep || verbose) && aln.mapping_quality() < min_mapq) {
        ++counts.counts[Counts::FilterName::min_mapq];
        keep = false;
    }
    if ((keep || verbose) && min_base_quality > 0 && min_base_quality_fraction > 0.0) {
        int mq_count = 0;
        const string& base_qualities = aln.quality();
        for (int i = 0; i < base_qualities.length(); ++i) {
            if (short(base_qualities[i]) >= min_base_quality) {
                ++mq_count;
            }
        }
        if ((double)mq_count / (double)base_qualities.length() < min_base_quality_fraction) {
            ++counts.counts[Counts::FilterName::min_base_qual];
            keep = false;
        }
    }
    if ((keep || verbose) && drop_split && is_split(aln)) {
        ++counts.counts[Counts::FilterName::split];
        keep = false;
    }
    if ((keep || verbose) && has_repeat(aln, repeat_size)) {
        ++counts.counts[Counts::FilterName::repeat];
        keep = false;
    }
    if ((keep || verbose) && defray_length && trim_ambiguous_ends(aln, defray_length)) {
        ++counts.counts[Counts::FilterName::defray];
        // We keep these, because the alignments get modified.
        // Unless the *entire* read gets trimmed
        if (aln.sequence().size() == 0) {
            keep = false;
            ++counts.counts[Counts::FilterName::defray_all];
        }
    }
    if ((keep || verbose) && downsample_probability != 1.0 && !sample_read(aln)) {
        ++counts.counts[Counts::FilterName::random];
        keep = false;
    }
    if (!keep) {
        ++counts.counts[Counts::FilterName::filtered];
    }
    
    return counts;
}

int ReadFilter::filter(istream* alignment_stream) {

    if(defray_length > 0 && graph == nullptr) {
        cerr << "HandleGraph (e.g. XG) required for end de-fraying" << endl;
        return 1;
    }

    // Keep an AlignmentEmitter to multiplex output from multiple threads.
    unique_ptr<AlignmentEmitter> emitter;
    
    if (write_output) {
        emitter = get_non_hts_alignment_emitter("-", "GAM",  map<string, int64_t>(), get_thread_count());
    }

    // keep counts of what's filtered to report (in verbose mode)
    vector<Counts> counts_vec(threads);

    function<void(Alignment&)> lambda = [&](Alignment& aln) {
#ifdef debug
        cerr << "Encountered read named \"" << aln.name() << "\" with " << aln.sequence().size()
            << " bp sequence and " << aln.quality().size() << " quality values" << endl;
#endif    
        Counts aln_counts = filter_alignment(aln);
        counts_vec[omp_get_thread_num()] += aln_counts;
        if ((aln_counts.keep() != complement_filter) && write_output) {
            emitter->emit_single(std::move(aln));
        }
    };

    function<void(Alignment&, Alignment&)> pair_lambda = [&](Alignment& aln1, Alignment& aln2) {
        Counts aln_counts = filter_alignment(aln1);
        aln_counts += filter_alignment(aln2);
        if (filter_on_all) {
            // Unless both reads were filtered out (total filtered count == 2), keep the read.
            aln_counts.set_paired_all();
        } else {
            // Either read failing is sufficient to scuttle the pair.
            // So if we filter out one end for any reason, we filter out the other as well.
            aln_counts.set_paired_any();
        }
        counts_vec[omp_get_thread_num()] += aln_counts;
        if ((aln_counts.keep() != complement_filter) && write_output) {
            emitter->emit_pair(std::move(aln1), std::move(aln2));
        }
        
    };
    
    if (interleaved) {
        vg::io::for_each_interleaved_pair_parallel(*alignment_stream, pair_lambda);
    } else {
        vg::io::for_each_parallel(*alignment_stream, lambda);
    }

    if (verbose) {
        Counts& counts = counts_vec[0];
        for (int i = 1; i < counts_vec.size(); ++i) {
            counts += counts_vec[i];
        }
        cerr << counts;
    }    
    return 0;
}

ostream& operator<<(ostream& os, const ReadFilter::Counts& counts) {
    os << "Total Filtered:                " << counts.counts[ReadFilter::Counts::FilterName::filtered] << " / "
       << counts.counts[ReadFilter::Counts::FilterName::read] << endl
       << "Read Name Filter:              " << counts.counts[ReadFilter::Counts::FilterName::wrong_name] << endl
       << "refpos Contig Filter:          " << counts.counts[ReadFilter::Counts::FilterName::wrong_refpos] << endl
       << "Feature Filter:                " << counts.counts[ReadFilter::Counts::FilterName::excluded_feature] << endl
       << "Min Identity Filter:           " << counts.counts[ReadFilter::Counts::FilterName::min_score] << endl
       << "Min Secondary Identity Filter: " << counts.counts[ReadFilter::Counts::FilterName::min_sec_score] << endl
       << "Max Overhang Filter:           " << counts.counts[ReadFilter::Counts::FilterName::max_overhang] << endl
       << "Min End Match Filter:          " << counts.counts[ReadFilter::Counts::FilterName::min_end_matches] << endl
       << "Split Read Filter:             " << counts.counts[ReadFilter::Counts::FilterName::split] << endl
       << "Repeat Ends Filter:            " << counts.counts[ReadFilter::Counts::FilterName::repeat] << endl
       << "All Defrayed Filter:           " << counts.counts[ReadFilter::Counts::FilterName::defray_all] << endl
       << "Min Quality Filter:            " << counts.counts[ReadFilter::Counts::FilterName::min_mapq] << endl
       << "Min Base Quality Filter:       " << counts.counts[ReadFilter::Counts::FilterName::min_base_qual] << endl
       << "Random Filter:                 " << counts.counts[ReadFilter::Counts::FilterName::random] << endl
       << endl;
    return os;
}

}
