#include "flow_caller.hpp"
#include "algorithms/expand_context.hpp"
#include "annotation.hpp"
#include "gref.hpp"
#include "traversal_clusters.hpp"

//#define debug

namespace vg {

FlowCaller::FlowCaller(const PathPositionHandleGraph& graph,
                       SupportBasedSnarlCaller& snarl_caller,
                       SnarlManager& snarl_manager,
                       const string& sample_name,
                       TraversalFinder& traversal_finder,
                       const vector<string>& ref_paths,
                       const vector<size_t>& ref_path_offsets,
                       const vector<int>& ref_path_ploidies,
                       AlignmentEmitter* aln_emitter,
                       bool traversals_only,
                       bool gaf_output,
                       size_t trav_padding,
                       bool genotype_snarls,
                       const pair<size_t, size_t>& allele_length_range) :
    GraphCaller(snarl_caller, snarl_manager),
    VCFOutputCaller(sample_name),
    GAFOutputCaller(aln_emitter, sample_name, ref_paths, trav_padding),
    graph(graph),
    traversal_finder(traversal_finder),
    ref_paths(ref_paths),
    traversals_only(traversals_only),
    gaf_output(gaf_output),
    genotype_snarls(genotype_snarls),
    allele_length_range(allele_length_range),
    candidates(graph, traversal_finder, ref_path_set, allele_length_range)
{
    for (int i = 0; i < ref_paths.size(); ++i) {
        ref_offsets[ref_paths[i]] = i < ref_path_offsets.size() ? ref_path_offsets[i] : 0;
        ref_path_set.insert(ref_paths[i]);
        ref_ploidies[ref_paths[i]] = i < ref_path_ploidies.size() ? ref_path_ploidies[i] : 2;
    }

}
   
FlowCaller::FlowCaller(const PathPositionHandleGraph& graph,
                       SupportBasedSnarlCaller& snarl_caller,
                       SnarlManager& snarl_manager,
                       const string& sample_name,
                       TraversalFinder& traversal_finder,
                       const vector<string>& ref_paths,
                       const vector<size_t>& ref_path_offsets,
                       const vector<int>& ref_path_ploidies,
                       AlignmentEmitter* aln_emitter,
                       bool traversals_only,
                       bool gaf_output,
                       size_t trav_padding,
                       bool genotype_snarls,
                       const pair<size_t, size_t>& allele_length_range,
                       bool nested,
                       bool star_allele) :
    GraphCaller(snarl_caller, snarl_manager),
    VCFOutputCaller(sample_name),
    GAFOutputCaller(aln_emitter, sample_name, ref_paths, trav_padding),
    graph(graph),
    traversal_finder(traversal_finder),
    ref_paths(ref_paths),
    traversals_only(traversals_only),
    gaf_output(gaf_output),
    genotype_snarls(genotype_snarls),
    allele_length_range(allele_length_range),
    nested(nested),
    star_allele(star_allele),
    candidates(graph, traversal_finder, ref_path_set, allele_length_range)
{
    for (int i = 0; i < ref_paths.size(); ++i) {
        ref_offsets[ref_paths[i]] = i < ref_path_offsets.size() ? ref_path_offsets[i] : 0;
        ref_path_set.insert(ref_paths[i]);
        ref_ploidies[ref_paths[i]] = i < ref_path_ploidies.size() ? ref_path_ploidies[i] : 2;
    }
}

FlowCaller::~FlowCaller() {

}

bool FlowCaller::call_snarl(const Snarl& managed_snarl) {
    // Entry point: call with no parent context
    return call_snarl_internal(managed_snarl, "", make_pair(0, 0), nullptr);
}

bool FlowCaller::call_snarl_internal(const Snarl& managed_snarl,
                                      const string& parent_ref_path_name,
                                      pair<size_t, size_t> parent_ref_interval,
                                      const ChildTraversalSets* parent_child_trav_sets) {
    // The reference path, the site oriented along it, and its candidate traversals.
    CandidateFinder::Site site;
    if (!candidates.find(managed_snarl, parent_ref_path_name, parent_ref_interval,
                         parent_child_trav_sets, false, site)) {
        return false;
    }
    Snarl& snarl = site.snarl;
    const string& ref_path_name = site.ref_path_name;
    const tuple<int64_t, int64_t, bool, step_handle_t, step_handle_t>& ref_interval =
        site.ref_interval;
    const bool use_parent_interval = site.use_parent_interval;
    vector<SnarlTraversal>& travs = site.travs;
    int ref_trav_idx = site.ref_trav_idx;
    const auto& support_finder = dynamic_cast<SupportBasedSnarlCaller&>(snarl_caller).get_support_finder();

    bool ret_val = true;
    vector<int> trav_genotype;  // Declared outside block so we can pass to children
    // The contig's ploidy, or a per-region override's.
    int ploidy = ploidy_regions.ploidy_at(ref_path_name, get<0>(ref_interval),
                                          ref_offset_of(ref_offsets, ref_path_name),
                                          ref_ploidy_of(ref_ploidies, ref_path_name));

    if (traversals_only) {
        assert(gaf_output);
        pair<string, int64_t> pos_info = get_ref_position(graph, snarl, ref_path_name, ref_offset_of(ref_offsets, ref_path_name));
        emit_gaf_traversals(graph, print_snarl(snarl), travs, ref_trav_idx, pos_info.first, pos_info.second, &support_finder);
    } else if (parent_child_trav_sets != nullptr && !parent_child_trav_sets->empty()) {
        // Genotype using bounded search over traversal sets from parent
        // Each set contains traversals consistent with one parent allele
        ploidy = parent_child_trav_sets->size();

        // Track which set each traversal index belongs to (for phase consistency)
        // set_membership[i] = which parent allele set traversal i came from, or -1 if from finder
        vector<int> set_membership(travs.size(), -1);

        // Merge traversals from sets into travs, tracking membership

        for (int set_idx = 0; set_idx < ploidy; ++set_idx) {
            const TraversalSet& tset = (*parent_child_trav_sets)[set_idx];

            if (tset.empty()) {
                // Empty set means parent allele doesn't traverse this child (star allele)
                continue;
            }

            // Add traversals from this set to travs (avoiding duplicates)
            // Keep track of indices for this set
            for (const SnarlTraversal& trav : tset) {
                // Check if this traversal already exists in travs
                int match_idx = -1;
                for (int i = 0; i < travs.size() && match_idx < 0; ++i) {
                    if (travs[i] == trav) {
                        match_idx = i;
                    }
                }

                if (match_idx < 0) {
                    // New traversal - add it
                    match_idx = travs.size();
                    travs.push_back(trav);
                    set_membership.push_back(set_idx);
                } else if (set_membership[match_idx] < 0) {
                    // Traversal was from finder, now claim it for this set
                    set_membership[match_idx] = set_idx;
                }
                // Note: if already claimed by another set, that's fine (shared region)
            }
        }

        // Which parent haplotypes actually traverse this child? A parent allele with
        // an empty traversal set skips the child entirely, and gets a star or
        // missing allele rather than a genotype.
        vector<int> traversing_sets;
        for (int set_idx = 0; set_idx < ploidy; ++set_idx) {
            if (!(*parent_child_trav_sets)[set_idx].empty()) {
                traversing_sets.push_back(set_idx);
            }
        }

        unique_ptr<SnarlCaller::CallInfo> trav_call_info;
        int marker = star_allele ? STAR_ALLELE_MARKER : MISSING_ALLELE_MARKER;

        if (traversing_sets.empty()) {
            // No parent allele traverses this child at all.
            trav_genotype.assign(ploidy, marker);
        } else {
            // Genotype at the ploidy that passes through the site, not at the parent's ploidy, since
            // a site only one strand reaches is not diploid. genotype() returns a sorted multiset, so
            // the alleles are then placed on the strands that pass through.
            int effective_ploidy = (int)traversing_sets.size();
            vector<int> called_alleles;
            std::tie(called_alleles, trav_call_info) = snarl_caller.genotype(
                snarl, travs, ref_trav_idx, effective_ploidy, ref_path_name,
                make_pair(get<0>(ref_interval), get<1>(ref_interval)));

            // Scatter the called alleles back onto the traversing haplotypes,
            // leaving the others as star/missing.
            trav_genotype.assign(ploidy, marker);
            for (size_t j = 0; j < traversing_sets.size() && j < called_alleles.size(); ++j) {
                trav_genotype[traversing_sets[j]] = called_alleles[j];
            }
        }

        // Emit variant with selected genotype
        bool added = true;

        // Only emit VCF if snarl is on reference path
        if (use_parent_interval) {
            added = true;
        } else if (!gaf_output) {
            added = emit_variant(graph, snarl_caller, snarl, travs, trav_genotype, ref_trav_idx, trav_call_info, ref_path_name,
                                 ref_offset_of(ref_offsets, ref_path_name), genotype_snarls, ploidy);
        } else {
            pair<string, int64_t> pos_info = get_ref_position(graph, snarl, ref_path_name, ref_offset_of(ref_offsets, ref_path_name));
            emit_gaf_variant(graph, print_snarl(snarl), travs, trav_genotype, ref_trav_idx, pos_info.first, pos_info.second, &support_finder);
        }

        ret_val = trav_genotype.size() == ploidy && added;
    } else {
        // Top-level snarl or no parent context - genotype from scratch using support
        unique_ptr<SnarlCaller::CallInfo> trav_call_info;
        std::tie(trav_genotype, trav_call_info) = snarl_caller.genotype(snarl, travs, ref_trav_idx, ploidy, ref_path_name,
                                                                        make_pair(get<0>(ref_interval), get<1>(ref_interval)));

        assert(trav_genotype.empty() || trav_genotype.size() == ploidy);

        bool added = true;
        if (!gaf_output) {
            added = emit_variant(graph, snarl_caller, snarl, travs, trav_genotype, ref_trav_idx, trav_call_info, ref_path_name,
                                 ref_offset_of(ref_offsets, ref_path_name), genotype_snarls, ploidy);
        } else {
            pair<string, int64_t> pos_info = get_ref_position(graph, snarl, ref_path_name, ref_offset_of(ref_offsets, ref_path_name));
            emit_gaf_variant(graph, print_snarl(snarl), travs, trav_genotype, ref_trav_idx, pos_info.first, pos_info.second, &support_finder);
        }

        ret_val = trav_genotype.size() == ploidy && added;
    }

    // In nested mode, recursively call child snarls
    if (nested && !trav_genotype.empty()) {
        // Find the managed snarl pointer so we can get its children
        const Snarl* managed_ptr = snarl_manager.into_which_snarl(snarl.start().node_id(), snarl.start().backward());
        if (managed_ptr) {
            const vector<const Snarl*>& children = snarl_manager.children_of(managed_ptr);
            for (const Snarl* child : children) {
                if (child && !snarl_manager.is_trivial(child, graph)) {
                    // Build ChildTraversalSets: one set per parent allele
                    // Each set contains all traversals through child consistent with that parent allele
                    ChildTraversalSets child_trav_sets;
                    bool any_real_traversals = false;

                    for (int allele_idx : trav_genotype) {
                        if (allele_idx >= 0 && allele_idx < travs.size()) {
                            // Find all traversals through child consistent with this parent traversal
                            TraversalSet tset = candidates.find_child_traversal_set(travs[allele_idx], *child);
                            if (!tset.empty()) {
                                any_real_traversals = true;
                            }
                            child_trav_sets.push_back(std::move(tset));
                        } else {
                            // Star/missing allele - pass empty set
                            child_trav_sets.push_back(TraversalSet());
                        }
                    }

                    // If no genotyped alleles traverse the child, skip it
                    if (!any_real_traversals) {
                        continue;
                    }

                    // Recursively call child with traversal sets
                    call_snarl_internal(*child, ref_path_name,
                                        make_pair(get<0>(ref_interval), get<1>(ref_interval)),
                                        &child_trav_sets);
                }
            }
        }
    }

    return ret_val;
}

string FlowCaller::vcf_header(const PathHandleGraph& graph, const vector<string>& contigs,
                              const vector<size_t>& contig_length_overrides) const {
    return snarl_caller_vcf_header(graph, contigs, contig_length_overrides, snarl_caller);
}
}
