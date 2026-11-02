#include "flow_caller.hpp"
#include "algorithms/expand_context.hpp"
#include "annotation.hpp"
#include "gref.hpp"
#include "snarl_traversal_text.hpp"
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
    GAFOutputCaller(aln_emitter, sample_name, ref_paths, trav_padding),
    vcf_writer(sample_name),
    site_emitter(vcf_writer),
    graph(graph),
    traversal_finder(traversal_finder),
    ref_paths(ref_paths),
    traversals_only(traversals_only),
    gaf_output(gaf_output),
    genotype_snarls(genotype_snarls),
    allele_length_range(allele_length_range),
    candidates(graph, dynamic_cast<SiteTraversalFinder&>(traversal_finder), ref_path_set, allele_length_range)
{
    site_emitter.genotype_snarls = genotype_snarls;
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
    GAFOutputCaller(aln_emitter, sample_name, ref_paths, trav_padding),
    vcf_writer(sample_name),
    site_emitter(vcf_writer),
    graph(graph),
    traversal_finder(traversal_finder),
    ref_paths(ref_paths),
    traversals_only(traversals_only),
    gaf_output(gaf_output),
    genotype_snarls(genotype_snarls),
    allele_length_range(allele_length_range),
    nested(nested),
    star_allele(star_allele),
    candidates(graph, dynamic_cast<SiteTraversalFinder&>(traversal_finder), ref_path_set, allele_length_range)
{
    site_emitter.genotype_snarls = genotype_snarls;
    for (int i = 0; i < ref_paths.size(); ++i) {
        ref_offsets[ref_paths[i]] = i < ref_path_offsets.size() ? ref_path_offsets[i] : 0;
        ref_path_set.insert(ref_paths[i]);
        ref_ploidies[ref_paths[i]] = i < ref_path_ploidies.size() ? ref_path_ploidies[i] : 2;
    }
}

void FlowCaller::set_max_snarl_edges(size_t edges) {
    candidates.set_max_snarl_edges(edges);
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
                                      const vector<TraversalSet>* parent_child_trav_sets) {
    // The reference path, the site oriented along it, and its candidate traversals.
    CandidateFinder::Site site;
    if (!candidates.find(managed_snarl, parent_ref_path_name, parent_ref_interval,
                         parent_child_trav_sets, false, site)) {
        return false;
    }
    const string& ref_path_name = site.ref_path_name;
    const bool has_parent_sets = parent_child_trav_sets != nullptr && !parent_child_trav_sets->empty();

    // The ploidy of a region covering the site, or else of the reference path.
    int ploidy = ploidy_regions.ploidy_at(ref_path_name, get<0>(site.ref_interval),
                                          ref_offset_of(ref_offsets, ref_path_name),
                                          ref_ploidy_of(ref_ploidies, ref_path_name));
    if (ploidy == 0 && !has_parent_sets) {
        // The sample has no copy of this region, so the site gets no record and its children are
        // not called. A nested site takes its ploidy from the traversal sets of its parent.
        return true;
    }

    bool ret_val = true;
    vector<int> trav_genotype;
    if (traversals_only) {
        assert(gaf_output);
        const auto& support_finder = dynamic_cast<SupportBasedSnarlCaller&>(snarl_caller).get_support_finder();
        auto [contig, position] = get_ref_position(graph, visit_handle(graph, site.snarl.start()),
                                                   visit_handle(graph, site.snarl.end()),
                                                   ref_path_name,
                                                   ref_offset_of(ref_offsets, ref_path_name));
        emit_gaf_traversals(graph, snarl_name(site.snarl, vcf_writer.translation()), site.travs,
                            site.ref_trav_idx, contig, position, &support_finder);
    } else {
        unique_ptr<SnarlCaller::CallInfo> trav_call_info;
        bool added = true;
        if (has_parent_sets) {
            ploidy = (int)parent_child_trav_sets->size();
            trav_genotype = genotype_from_parent_sets(site, *parent_child_trav_sets, trav_call_info);
            // Only a site on a reference path of its own gets a record.
            if (!site.use_parent_interval) {
                added = emit_site(site, trav_genotype, trav_call_info, ploidy);
            }
        } else {
            // Top-level snarl or no parent context - genotype from scratch using support
            std::tie(trav_genotype, trav_call_info) = snarl_caller.genotype(site.snarl, site.travs, site.ref_trav_idx, ploidy, ref_path_name,
                                                                            make_pair(get<0>(site.ref_interval), get<1>(site.ref_interval)));
            assert(trav_genotype.empty() || trav_genotype.size() == ploidy);
            added = emit_site(site, trav_genotype, trav_call_info, ploidy);
        }
        ret_val = trav_genotype.size() == ploidy && added;
    }

    if (nested && !trav_genotype.empty()) {
        call_children(site, trav_genotype);
    }

    return ret_val;
}

/// Add to `travs` each traversal of the traversal sets that `travs` does not have yet.
static void add_traversal_sets(vector<SnarlTraversal>& travs, const vector<TraversalSet>& sets) {
    for (const TraversalSet& tset : sets) {
        for (const SnarlTraversal& trav : tset) {
            if (std::ranges::find(travs, trav) == travs.end()) {
                travs.push_back(trav);
            }
        }
    }
}

vector<int> FlowCaller::genotype_from_parent_sets(CandidateFinder::Site& site,
                                                  const vector<TraversalSet>& parent_sets,
                                                  unique_ptr<SnarlCaller::CallInfo>& call_info) {
    add_traversal_sets(site.travs, parent_sets);

    // The parent strands that visit this site. The others get a star or missing allele.
    const auto ploidy = (int)parent_sets.size();
    vector<int> traversing_sets;
    for (int set_idx = 0; set_idx < ploidy; ++set_idx) {
        if (!parent_sets[set_idx].empty()) {
            traversing_sets.push_back(set_idx);
        }
    }

    int marker = star_allele ? STAR_ALLELE_MARKER : MISSING_ALLELE_MARKER;
    vector<int> genotype;
    if (traversing_sets.empty()) {
        genotype.assign(ploidy, marker);
        return genotype;
    }
    // Genotype at the number of strands that visit the site, since a site that only one strand
    // visits is not diploid. genotype() returns a sorted multiset, so we then put the alleles on
    // the strands that visit the site, in order.
    auto effective_ploidy = (int)traversing_sets.size();
    vector<int> called_alleles;
    std::tie(called_alleles, call_info) = snarl_caller.genotype(
        site.snarl, site.travs, site.ref_trav_idx, effective_ploidy, site.ref_path_name,
        make_pair(get<0>(site.ref_interval), get<1>(site.ref_interval)));
    genotype.assign(ploidy, marker);
    for (size_t j = 0; j < traversing_sets.size() && j < called_alleles.size(); ++j) {
        genotype[traversing_sets[j]] = called_alleles[j];
    }
    return genotype;
}

bool FlowCaller::emit_site(const CandidateFinder::Site& site, const vector<int>& genotype,
                           const unique_ptr<SnarlCaller::CallInfo>& call_info, int ploidy) {
    const size_t ref_offset = ref_offset_of(ref_offsets, site.ref_path_name);
    if (!gaf_output) {
        CalledSite called{visit_handle(graph, site.snarl.start()),
                          visit_handle(graph, site.snarl.end()),
                          site.ref_path_name,
                          (int64_t)ref_offset,
                          genotype,
                          site.ref_trav_idx,
                          site.travs.size(),
                          ploidy};
        return site_emitter.emit(
            graph, called,
            [&](size_t allele) {
                return make_pair(trav_string(graph, site.travs[allele], vcf_writer.translation()),
                                 visits_of(site.travs[allele]));
            },
            [&](const vector<int>& allele_traversals, const vector<int>& vcf_genotype,
                vcflib::Variant& record) {
                snarl_caller.update_vcf_info(site.snarl,
                                             traversals_of(site.travs, allele_traversals),
                                             vcf_genotype, call_info,
                                             vcf_writer.sample_name(), record);
            });
    }
    const auto& support_finder = dynamic_cast<SupportBasedSnarlCaller&>(snarl_caller).get_support_finder();
    auto [contig, position] = get_ref_position(graph, visit_handle(graph, site.snarl.start()),
                                               visit_handle(graph, site.snarl.end()),
                                               site.ref_path_name, ref_offset);
    emit_gaf_variant(graph, snarl_name(site.snarl, vcf_writer.translation()), site.travs,
                     genotype, site.ref_trav_idx, contig, position, &support_finder);
    return true;
}

vector<TraversalSet> FlowCaller::child_traversal_sets(const CandidateFinder::Site& site,
                                                      const vector<int>& genotype,
                                                      const Snarl& child) const {
    vector<TraversalSet> child_trav_sets;
    for (int allele_idx : genotype) {
        if (allele_idx >= 0 && allele_idx < site.travs.size()) {
            child_trav_sets.push_back(candidates.find_child_traversal_set(site.travs[allele_idx], child));
        } else {
            child_trav_sets.emplace_back();
        }
    }
    return child_trav_sets;
}

void FlowCaller::call_children(const CandidateFinder::Site& site, const vector<int>& genotype) {
    const Snarl* managed_ptr = snarl_manager.into_which_snarl(site.snarl.start().node_id(),
                                                              site.snarl.start().backward());
    if (managed_ptr == nullptr) {
        return;
    }
    for (const Snarl* child : snarl_manager.children_of(managed_ptr)) {
        if (child == nullptr || snarl_manager.is_trivial(child, graph)) {
            continue;
        }
        vector<TraversalSet> child_trav_sets = child_traversal_sets(site, genotype, *child);
        // Skip a child that no called allele visits.
        if (std::ranges::any_of(child_trav_sets, [](const TraversalSet& tset) { return !tset.empty(); })) {
            call_snarl_internal(*child, site.ref_path_name,
                                make_pair(get<0>(site.ref_interval), get<1>(site.ref_interval)),
                                &child_trav_sets);
        }
    }
}

}
