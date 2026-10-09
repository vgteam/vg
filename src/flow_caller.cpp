#include <atomic>
#include <charconv>
#include <chrono>
#include <cstdio>
#include <limits>

#include <omp.h>

#include "flow_caller.hpp"
#include "symbolic_allele.hpp"
#include "read_likelihood_caller.hpp"
#include "algorithms/expand_context.hpp"
#include "annotation.hpp"
#include "gref.hpp"
#include "traversal_clusters.hpp"
#include "utility.hpp"

//#define debug

namespace vg {
static thread_local int g_descent_depth = 0;
/// The place in the nesting tree of the snarl the direct pass is genotyping on this thread. Descent
/// runs on the calling thread, so it is set just before a child is genotyped and restored after,
/// and no other thread sees it.
static thread_local NestingPlacement nested_context;

void FlowCaller::report_descent_instrumentation() const {
    size_t total = 0;
    for (int d = 0; d < 16; ++d) {
        total += descent_counters.depth_hist[d].load();
    }
    if (total == 0) {
        return;   // no symbolic descent in this run
    }
    cerr << "[vg call] descent depth:";
    for (int d = 1; d < 16; ++d) {
        size_t n = descent_counters.depth_hist[d].load();
        if (n > 0) {
            cerr << " " << d << "=" << n;
        }
    }
    cerr << " (" << total << " child calls)" << endl;
    if (descent_counters.child_multi_crossing.load() > 0) {
        cerr << "[vg call] descent: " << descent_counters.child_multi_crossing.load()
             << " children a called traversal enters more than once; visits after the first are"
             << " masked, so each contributes one copy and its first crossing's distance" << endl;
    }
    cerr << "[vg call] descent skipped: " << descent_counters.skipped_no_copy.load()
         << " children no called allele reaches, " << descent_counters.skipped_no_ref.load()
         << " with no reference path through them" << endl;
    if (descent_counters.off_reference.load() > 0 || descent_counters.no_ref_recorded.load() > 0) {
        cerr << "[vg call] off-reference nested: " << descent_counters.off_reference.load()
             << " chains the reference does not cross were descended into, "
             << descent_counters.no_ref_recorded.load() << " recorded into the linkage layer with no line;"
             << " copies 0/1/2 = " << descent_counters.no_ref_copies[0].load() << "/"
             << descent_counters.no_ref_copies[1].load() << "/" << descent_counters.no_ref_copies[2].load() << endl;
    }
}

// The anchor gqn column for a staged site.
//
// `gq_fraction` was computed in the direct pass for the reads' best genotype, so on a record whose
// genotype the linkage model changed, it describes the abandoned genotype. Such a record gets the
// signed margin of its chosen genotype instead, as the VCF's GQN does (see
// ReadLikelihoodSnarlCaller::rewrite_quality_for_chosen_genotype).
//
// Returns the direct pass's value when the model did not change the call; the recomputed signed margin
// when it did; and NaN, written as ".", when it did but the margin cannot be recomputed.
double FlowCaller::anchor_gqn_for(const StagedSite& rec,
                                  const vector<int>& chosen) const {
    const auto* info =
        dynamic_cast<const ReadLikelihoodSnarlCaller::ReadLikelihoodCallInfo*>(rec.call_info.get());
    // The direct pass's value, with its "no gap to normalise" value (-1) turned into NaN, written as ".",
    // so that it stays distinct from the signed range [-1, 1].
    const double direct_value = (info == nullptr || info->gq_fraction < 0.0)
        ? std::numeric_limits<double>::quiet_NaN()
        : info->gq_fraction;
    const double blank = std::numeric_limits<double>::quiet_NaN();
    if (!linker.enabled()) {
        return direct_value;
    }
    const auto& moved = linker.collector()->moved_quality();
    const auto found = moved.find(rec.record_key);
    if (found == moved.end()) {
        return direct_value;   // linkage left the call alone, so the direct pass's value still holds
    }
    // The model changed the call, so the direct pass's value describes the wrong genotype, and any
    // failure below gives NaN rather than falling back to it.
    if (info == nullptr || info->genotype_lls.empty()) {
        return blank;
    }
    // The divisor and share the VCF's GQN uses (see
    // ReadLikelihoodSnarlCaller::rewrite_quality_for_chosen_genotype), so that the two agree.
    const LinkageCollector::DirectQuality& direct = found->second.direct;
    if (!(direct.achievable_gap > 0.0)) {
        return blank;   // no scale, and no honest pre-linkage value to fall back on
    }
    const double achievable_phred = 10.0 * direct.achievable_gap / log(10.0);

    // The chosen genotype, not rec.genotype, which is the direct pass's call before the linkage model:
    // the reads prefer that call, so its margin would have the wrong sign.
    vector<int> called = chosen;
    sort(called.begin(), called.end());
    const auto mine = info->genotype_lls.find(called);
    if (mine == info->genotype_lls.end()) {
        return blank;
    }
    // Only genotypes over the written alleles, as the VCF's GL has: the reference traversal and
    // the ones the chosen genotype names. A traversal with no ALT could otherwise beat the call.
    set<int> emitted(called.begin(), called.end());
    if (rec.ref_trav_idx >= 0) {
        emitted.insert(rec.ref_trav_idx);
    }
    double best_other = -numeric_limits<double>::infinity();
    for (const auto& entry : info->genotype_lls) {
        if (entry.first == called) {
            continue;
        }
        bool all_emitted = true;
        for (int a : entry.first) {
            if (emitted.count(a) == 0) {
                all_emitted = false;
                break;
            }
        }
        if (all_emitted) {
            best_other = max(best_other, entry.second);
        }
    }
    if (!std::isfinite(best_other)) {
        return blank;
    }
    // Nats to phred, matching the VCF's GL, which is log10.
    const double margin_phred = 10.0 * (mine->second - best_other) / log(10.0);
    return min(1.0, max(-1.0, margin_phred / achievable_phred * direct.explained_share));
}

void FlowCaller::collect_anchors_for_record(const StagedSite& rec,
                                            const vector<int>& genotype) {
    if (!anchor_collector.is_enabled() || rec.call_info == nullptr) {
        return;
    }
    const auto* info =
        dynamic_cast<const ReadLikelihoodSnarlCaller::ReadLikelihoodCallInfo*>(rec.call_info.get());
    if (info == nullptr || info->anchor_evidence == nullptr) {
        // A genotype derived from a parent rather than scored here, or a run whose caller is not the
        // read-likelihood one. There are no per-read responsibilities to partition on.
        return;
    }
    anchor_collector.collect(*info->anchor_evidence, info->explained_share,
                             phase_table.phase_ordered_genotype(rec.record_key, genotype),
                             phase_table.haploid_slot(rec.record_key, genotype),
                             print_snarl(rec.snarl),
                             anchor_collector.wants_leaf_test() ? snarl_is_leaf(rec.snarl) : true,
                             anchor_gqn_for(rec, genotype), read_strands, rec.record_key);
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
    allele_length_range(allele_length_range)
{
    for (int i = 0; i < ref_paths.size(); ++i) {
        ref_offsets[ref_paths[i]] = i < ref_path_offsets.size() ? ref_path_offsets[i] : 0;
        ref_path_set.insert(ref_paths[i]);
        ref_ploidies[ref_paths[i]] = i < ref_path_ploidies.size() ? ref_path_ploidies[i] : 2;
    }
    install_record_steps();
    linker.set_site_reader(GenotypeLinker::SiteReader{
        .graph = &this->graph,
        .caller = &snarl_caller,
        .spell = [this](const SnarlTraversal& trav) { return trav_string(this->graph, trav); },
    });

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
    star_allele(star_allele)
{
    for (int i = 0; i < ref_paths.size(); ++i) {
        ref_offsets[ref_paths[i]] = i < ref_path_offsets.size() ? ref_path_offsets[i] : 0;
        ref_path_set.insert(ref_paths[i]);
        ref_ploidies[ref_paths[i]] = i < ref_path_ploidies.size() ? ref_path_ploidies[i] : 2;
    }
    install_record_steps();
    linker.set_site_reader(GenotypeLinker::SiteReader{
        .graph = &this->graph,
        .caller = &snarl_caller,
        .spell = [this](const SnarlTraversal& trav) { return trav_string(this->graph, trav); },
    });
}

FlowCaller::~FlowCaller() {

}

void FlowCaller::install_record_steps() {
    record_steps.phase = [this](const Snarl& site, const vector<int>& site_genotype,
                                const map<int, int>& trav_to_allele, string& gt) {
        return phase_record_genotype(site, site_genotype, trav_to_allele, gt);
    };
    // The read-likelihood genotyper writes GL in colexicographic order, and the support-based one
    // in i-major order.
    record_steps.gl_layout = [](const SnarlCaller::CallInfo* call_info) {
        return dynamic_cast<const ReadLikelihoodSnarlCaller::ReadLikelihoodCallInfo*>(call_info)
                   != nullptr
                   ? GLLayout::Colexicographic
                   : GLLayout::IMajor;
    };
    record_steps.write_blocks = [this](const PathPositionHandleGraph& graph, const Snarl& site,
                                       const vector<SnarlTraversal>& travs,
                                       const vector<int>& genotype, int ref_trav_idx,
                                       const SiteRecord& record, GLLayout gl_layout,
                                       bool genotype_snarls) {
        return block_records.write(graph, site, travs, genotype, ref_trav_idx, sample_name,
                                   translation, record, gl_layout, genotype_snarls,
                                   [this](vcflib::Variant& line, size_t block) {
                                       return add_variant(line, block);
                                   });
    };
    // The linkage model gets the site whether or not it has a line. A parent written as the
    // reference still has two alleles, which differ only inside its children, and the children
    // need them to know which strand carries the chain. In VCF allele numbering such a parent is
    // 0/0; only in traversal space is it heterozygous.
    record_steps.site_filed = [this](const Snarl& site, const map<int, int>& trav_to_allele,
                                     size_t traversal_count, bool has_line) {
        if (!linker.enabled()) {
            return;
        }
        // The site was recorded when it was genotyped. What remains is the traversal-to-VCF-allele
        // map, which depends on the alleles the record chose, and whether a line was written. A
        // site written as blocks gives an empty map, since each block numbers its own alleles.
        vector<int> trav_to_allele_vec(traversal_count, -1);
        for (const auto& kv : trav_to_allele) {
            if (kv.first >= 0 && (size_t)kv.first < trav_to_allele_vec.size()) {
                trav_to_allele_vec[kv.first] = kv.second;
            }
        }
        linker.collector()->set_allele_map(record_key_of(site), trav_to_allele_vec, has_line);
    };
}

void FlowCaller::call_top_level_snarls(const HandleGraph& graph, RecurseType recurse_type) {
    GraphCaller::call_top_level_snarls(graph, recurse_type);
    if (show_progress) {
        report_descent_instrumentation();
    }
}

bool FlowCaller::call_snarl(const Snarl& managed_snarl) {
    // Entry point: call with no parent context
    return call_snarl_internal(managed_snarl, "", make_pair(0, 0), nullptr);
}

TraversalSet FlowCaller::find_child_traversal_set(const SnarlTraversal& parent_trav,
                                                   const Snarl& child) const {
    TraversalSet result;

    // First, check if the parent traversal goes through this child snarl
    // by finding the child's start and end nodes in the parent
    nid_t child_start_id = child.start().node_id();
    nid_t child_end_id = child.end().node_id();
    bool found_start = false, found_end = false;

    for (int i = 0; i < parent_trav.visit_size(); ++i) {
        nid_t visit_id = parent_trav.visit(i).node_id();
        if (visit_id == child_start_id) found_start = true;
        if (visit_id == child_end_id) found_end = true;
    }

    // If parent doesn't traverse the child, return empty set (star allele case)
    if (!found_start || !found_end) {
        return result;
    }

    // Use the traversal finder to enumerate all traversals through the child
    FlowTraversalFinder* flow_finder = dynamic_cast<FlowTraversalFinder*>(&traversal_finder);
    if (flow_finder != nullptr) {
        auto weighted_travs = flow_finder->find_weighted_traversals(child, false);
        result = std::move(weighted_travs.first);
    } else {
        result = traversal_finder.find_traversals(child);
    }

    return result;
}

int64_t FlowCaller::base_offset_of_child(const SnarlTraversal& trav, const Snarl& child) const {
    const int entry = ChildPlacer::offset_of_child(trav, child);
    if (entry < 0) {
        return -1;
    }
    int64_t bases = 0;
    for (int i = 0; i < entry && i < trav.visit_size(); ++i) {
        if (trav.visit(i).has_snarl()) {
            continue;
        }
        bases += (int64_t)graph.get_length(graph.get_handle(trav.visit(i).node_id()));
    }
    return bases;
}

size_t FlowCaller::offset_along_genotype(const vector<SnarlTraversal>& travs,
                                         const vector<int>& genotype, const Snarl& child) const {
    for (int allele : genotype) {
        if (allele < 0 || allele >= (int)travs.size()) {
            continue;
        }
        const int64_t within = base_offset_of_child(travs[allele], child);
        if (within >= 0) {
            return (size_t)within;
        }
    }
    return 0;
}

int FlowCaller::child_ploidy(const vector<ChildPlacer::TraversalNodeIndex>& visits,
                             const vector<int>& genotype,
                             const Snarl& child, int cap) const {
    int copies = 0;
    bool capped = false;

    for (int allele : genotype) {
        if (allele < 0 || allele >= (int)visits.size()) {
            continue;   // star or missing: that haplotype contributes no copy here
        }
        int crossings = ChildPlacer::crossings_of_child(visits[allele], child);
        if (crossings > 1) {
            capped = true;
            crossings = 1;   // a cycle or tandem duplication; see the header comment
        }
        copies += crossings;
    }
    if (capped) {
        // Counted and reported once per run.
        ++descent_counters.child_multi_crossing;
    }
    return min(copies, cap);
}

void FlowCaller::set_stage_records(bool defer) {
    if (defer) {
        // Sized here rather than inside the parallel region that writes it.
        staged_sites.start(max((size_t)get_thread_count(), (size_t)omp_get_max_threads()));
    }
}

// The CallInfo is kept because update_vcf_info reads it when the record is rendered, to map the
// written alleles back to matrix columns, index GL and compute QUAL.
unique_ptr<StagedSite> FlowCaller::stage_render_record(
        const Snarl& snarl, const vector<int>& trav_genotype, int ref_trav_idx,
        unique_ptr<SnarlCaller::CallInfo>& call_info,
        const string& ref_path_name, int ref_offset, int ploidy) {
    if (!staged_sites.active()) {
        return nullptr;
    }
    unique_ptr<StagedSite> rec(new StagedSite());
    rec->snarl = snarl;
    rec->ref_path_name = ref_path_name;
    rec->ref_offset = ref_offset;
    rec->ref_trav_idx = ref_trav_idx;
    rec->genotype = trav_genotype;
    rec->ploidy = ploidy;
    rec->record_key = record_key_of(snarl);
    rec->level = 0;
    rec->call_info = std::move(call_info);
    // `travs` is not moved here: descent runs after the emit and reads `travs` to find which
    // children the called alleles reach. The caller completes the record after descent.
    return rec;
}

bool FlowCaller::snarl_is_leaf(const Snarl& snarl) const {
    // Through `manage`, not the address of this Snarl. `SnarlManager::record` casts a Snarl* to its
    // record, which is valid only for a Snarl the manager owns, and the Snarls here are copies.
    // `manage` throws for a snarl the manager does not own, as a nested chain reached by recursion
    // may be, so the call is guarded, and made only when --anchors-leaf-only needs the answer.
    try {
        const Snarl* managed = snarl_manager.manage(snarl);
        return managed != nullptr && snarl_manager.children_of(managed).empty();
    } catch (const std::runtime_error&) {
        // No answer, so treat it as a leaf rather than drop the site.
        return true;
    }
}

unordered_map<size_t, array<int, 3>> FlowCaller::chosen_snapshot() {
    // Each record's chosen pair and ploidy, keyed by record key.
    unordered_map<size_t, array<int, 3>> out;
    if (!linker.enabled()) {
        return out;
    }
    // Looked up on several threads, then filed in record order, so that the map is built exactly
    // as one loop over the records builds it.
    const vector<StagedSite*> records = staged_sites.in_order();
    vector<size_t> keys(records.size());
    for (size_t i = 0; i < records.size(); ++i) {
        keys[i] = records[i]->record_key;
    }
    vector<array<int, 3>> chosen;
    vector<char> found;
    linker.collector()->chosen_traversals_for(keys, chosen, found);
    for (size_t i = 0; i < keys.size(); ++i) {
        if (found[i]) {
            out[keys[i]] = chosen[i];
        }
    }
    return out;
}

size_t FlowCaller::snapshot_digest(const unordered_map<size_t, array<int, 3>>& snap) {
    // Independent of order, since the snapshot is a hash map: each record's contribution is
    // combined with a commutative mix.
    size_t acc = snap.size() * 1000003ULL;
    for (const auto& kv : snap) {
        size_t h = kv.first;
        h = h * 1000003ULL + (size_t)(kv.second[0] + 3);
        h = h * 1000003ULL + (size_t)(kv.second[1] + 3);
        h = h * 1000003ULL + (size_t)kv.second[2];
        acc ^= h + 0x9e3779b97f4a7c15ULL + (acc << 6) + (acc >> 2);
    }
    return acc;
}

size_t FlowCaller::chosen_changed(const unordered_map<size_t, array<int, 3>>& before,
                                  const unordered_map<size_t, array<int, 3>>& after) {
    size_t moved = 0;
    for (const auto& kv : after) {
        auto found = before.find(kv.first);
        if (found == before.end() || found->second != kv.second) {
            ++moved;
        }
    }
    // A record that had a chosen answer and now has none has changed too.
    for (const auto& kv : before) {
        if (after.count(kv.first) == 0) {
            ++moved;
        }
    }
    return moved;
}

void FlowCaller::apply_read_phasing() {
    if (!read_phasing || !linker.enabled() || phase_table.calls().empty()) {
        return;
    }
    // Reset, since re-genotyping calls this again on the new genotypes, and the report should
    // describe the phase the output carries.
    read_phasing_counters = ReadPhasingCounters();
    // Index the phasing by record key, the last one written winning.
    const std::unordered_map<size_t, size_t> phase_index = phase_table.index();
    const vector<LinkageCollector::PhaseCall>& calls = phase_table.calls();

    // Kept in `read_strands`, since re-genotyping uses these sites.
    vector<PhaseSite>& sites = read_strands.sites();
    sites.clear();
    // Each record's site is built from its own evidence alone, so the sites are built on several
    // threads, a block of records at a time into the block's own list, and then gathered in record
    // order. Phase sets are numbered in the order they are first seen, so a site's number is given
    // in that gathering pass, in record order.
    const vector<StagedSite*> records = staged_sites.in_order(true);
    const size_t block_records = 4096;
    const size_t n_blocks = (records.size() + block_records - 1) / block_records;
    vector<vector<PhaseSite>> block_sites(n_blocks);
    // For each site in a block, its PhaseCall's index in `calls`.
    vector<vector<size_t>> block_calls(n_blocks);
#pragma omp parallel for schedule(dynamic, 1)
    for (size_t b = 0; b < n_blocks; ++b) {
        const size_t end = min(records.size(), (b + 1) * block_records);
        for (size_t r = b * block_records; r < end; ++r) {
            const StagedSite& rec = *records[r];
            const auto found = phase_index.find(rec.record_key);
            if (found == phase_index.end()) {
                continue;
            }
            const LinkageCollector::PhaseCall& pc = calls[found->second];
            if (pc.ploidy != 2 || pc.trav_first < 0 || pc.trav_second < 0
                || pc.trav_first == pc.trav_second) {
                // Homozygous, haploid, or unplaced: no two strands to order.
                continue;
            }
            const auto* info = dynamic_cast<
                const ReadLikelihoodSnarlCaller::ReadLikelihoodCallInfo*>(rec.call_info.get());
            if (info == nullptr) {
                continue;
            }
            PhaseReadEvidence converted;
            const PhaseReadEvidence* pe = info->read_phasing_evidence(converted);
            if (pe == nullptr) {
                continue;
            }
            const size_t a0 = (size_t)pc.trav_first, a1 = (size_t)pc.trav_second;
            if (a0 >= pe->n_alleles || a1 >= pe->n_alleles) {
                // The PhaseCall names a traversal this site's matrix does not have. Skipped, since
                // reading another column would take the phase from the wrong allele.
                continue;
            }
            // Slot order is the PhaseCall's order, so slot 0 is strand 0, as for GT's first field and
            // the anchor file's slot column.
            PhaseSite site = reduce_to_pair(*pe, a0, a1);
            if (site.read_key.empty()) {
                continue;
            }
            site.record_key = rec.record_key;
            site.position = pc.position;
            block_sites[b].push_back(std::move(site));
            block_calls[b].push_back(found->second);
        }
    }
    size_t total_sites = 0;
    for (const vector<PhaseSite>& block : block_sites) {
        total_sites += block.size();
    }
    sites.reserve(total_sites);
    for (size_t b = 0; b < n_blocks; ++b) {
        for (size_t i = 0; i < block_sites[b].size(); ++i) {
            const LinkageCollector::PhaseCall& pc = calls[block_calls[b][i]];
            block_sites[b][i].phase_set = phase_set_id(pc.contig, pc.phase_set);
            sites.push_back(std::move(block_sites[b][i]));
        }
        vector<PhaseSite>().swap(block_sites[b]);
    }
    if (sites.empty()) {
        return;
    }

    read_strands.flips() = read_phase_flips(sites, read_phasing_params, read_phasing_counters);

    // Apply by swapping the chosen pair's order, and carry the swaps down the nesting tree. Nested
    // sites are reordered too. Under -A, block records spell the phase in their ALTs, so
    // reordering a nested site can change its GT's allele numbers. Every recorded chain is linked,
    // including one whose line an enclosing block's ALT spells (`reported_inline`): it still has
    // anchors, read from its strand, and its children's strands depend on its own. A dropped
    // chain is left out, since the sample does not carry it or anything inside it.
    vector<PhaseTable::NestedLink> links;
    staged_sites.for_each([&](const StagedSite& rec) {
        if (!rec.dropped && phase_index.count(rec.record_key) != 0) {
            links.push_back({rec.record_key, rec.parent_record_key, rec.level});
        }
    });
    read_phasing_counters.strands_rederived +=
        phase_table.swap_strands(read_strands.flips(), std::move(links));

    const ReadPhasingCounters& c = read_phasing_counters;
    cerr << "[vg call] read phasing: " << c.sites << " het sites, " << c.reliable
         << " reliable, " << c.chains << " blocks, " << c.breaks << " chain breaks ("
         << c.breaks_no_reads << " with no spanning read), " << c.hung
         << " sites hung off the chain (" << c.hung_no_reads << " with no read), " << c.flipped
         << " re-phased against the panel, " << c.strands_rederived
         << " nested strands carried with their parent" << endl;
    if (c.demoted_incoherent > 0) {
        cerr << "[vg call] read phasing: " << c.demoted_incoherent
             << " sites demoted for low phase coherence over " << c.coherence_rounds_run
             << " rounds"
             << (c.coherence_unconverged
                     ? " -- " + std::to_string(c.coherence_unconverged)
                           + " chains were STILL demoting at the round cap, so their chain is not"
                             " a coherent fixed point"
                     : " (every chain reached a coherent fixed point)")
             << endl;
    }
}

bool FlowCaller::apply_regenotyping() {
    const vector<LinkageCollector::PhaseCall>& calls = phase_table.calls();
    if (!regenotype || !linker.enabled() || calls.empty()) {
        return false;
    }
    const vector<PhaseSite>& phase_sites = read_strands.sites();
    const unordered_set<size_t>& phase_flips = read_strands.flips();
    // Reset the counters first, before `accumulate_lambda` fills the read counts, so that the report
    // describes this round. The calibration table and fitted temper are kept: they are set once, on
    // the first round.
    regenotype_counters = RegenotypeCounters();
    temper_fit.restore(regenotype_counters);

    // Lambda over every site read phasing covered, in one pass, into a table keyed by read.
    LambdaTable lambda;
    accumulate_lambda(phase_sites, phase_flips, lambda, regenotype_counters);

    // Each site's phase set, the last PhaseCall written winning, as in
    // `PhaseTable::freeze_for_render`. A read's strand is usable only at sites of the phase set it
    // was found in.
    unordered_map<size_t, size_t> site_phase_set;
    site_phase_set.reserve(calls.size() * 2);
    // And the allele the chain puts on strand 0 at each diploid site, against which the reads'
    // preferred order is reported.
    unordered_map<size_t, int> site_strand0;
    site_strand0.reserve(calls.size() * 2);
    for (const LinkageCollector::PhaseCall& pc : calls) {
        site_phase_set[pc.record_key] = phase_set_id(pc.contig, pc.phase_set);
        site_strand0[pc.record_key] = pc.ploidy == 2 ? pc.trav_first : -1;
    }

    // Which strand of its parent each nested ploidy-1 chain sits on. `nested_strand` was set in the
    // linkage pass and corrected when its parent's pair was swapped, so it is in the same frame as
    // Lambda for reads of the chain's phase set: strand 0 of that phase set.
    unordered_map<size_t, int> haploid_strand;
    if (regenotype_params.haploid_include) {
        for (const LinkageCollector::PhaseCall& pc : calls) {
            if (pc.ploidy == 1 && pc.nested_strand >= 0) {
                haploid_strand[pc.record_key] = pc.nested_strand == 0 ? 1 : -1;
            }
        }
    }

    double temper = regenotype_params.temper;
    double ceiling = regenotype_params.ceiling < 0.0 ? 1.0 : regenotype_params.ceiling;
    if (temper < 0.0) {
        // The temper is fitted once, on the first round, and kept. It describes how reliable the
        // reads' summed strand log-odds are, not which genotypes are called, and fitting it again
        // each round would feed each round's result into the next fit.
        if (regenotype_counters.fitted_temper > 0.0) {
            temper = regenotype_counters.fitted_temper;
            ceiling = regenotype_counters.fitted_ceiling;
        } else {
            fit_calibration(phase_sites, phase_flips, lambda, regenotype_params, temper, ceiling,
                            regenotype_counters);
        }
    } else {
        regenotype_counters.fitted_temper = temper;
        regenotype_counters.fitted_ceiling = ceiling;
    }
    temper_fit.keep(regenotype_counters);

    // Each site's own PhaseSite, so that its term can be subtracted from its reads' log-odds.
    unordered_map<size_t, const PhaseSite*> site_by_key;
    site_by_key.reserve(phase_sites.size() * 2);
    for (const PhaseSite& ps : phase_sites) {
        site_by_key[ps.record_key] = &ps;
    }

    ofstream ledger;
    const bool want_ledger = !regenotype_ledger.empty();
    if (want_ledger) {
        ledger.open(regenotype_ledger);
        if (!ledger) {
            cerr << "error [vg call]: cannot write --regeno-ledger " << regenotype_ledger << endl;
            exit(1);
        }
        ledger << "#regeno-ledger-version\t1" << endl;
        ledger << "#temper\t" << temper << endl;
        ledger << "#snarl\tcontig\tposition\tploidy\tcalled\tproposed\tdelta_ln\treads" << endl;
    }

    // Parallel over one flat list of records, strided across threads; `lambda`, `site_by_key` and
    // `phase_flips` are read only. Counters and ledger rows are kept per thread and merged
    // afterwards.
    const vector<StagedSite*> all_records = staged_sites.in_order();
    const size_t n_queues = max<size_t>(1, staged_sites.queue_count());
    // Only a read-likelihood caller makes the CallInfos corrected below.
    const auto* rl_caller = dynamic_cast<const ReadLikelihoodSnarlCaller*>(&snarl_caller);
    vector<RegenotypeCounters> thread_counters(n_queues);
    // Ledger rows are sorted before writing, since which thread handles a record depends on
    // scheduling.
    struct LedgerRow { string contig; size_t position; string snarl; string text; };
    vector<vector<LedgerRow>> thread_ledger(n_queues);
    vector<size_t> thread_moved(n_queues, 0);
#pragma omp parallel for schedule(dynamic, 1)
    for (size_t qi = 0; qi < n_queues; ++qi) {
        RegenotypeCounters& counters = thread_counters[qi];
        unordered_map<uint64_t, double> own;
        size_t moved = 0;
        for (size_t ri = qi; ri < all_records.size(); ri += n_queues) {
            StagedSite& rec = *all_records[ri];
            auto* info = dynamic_cast<ReadLikelihoodSnarlCaller::ReadLikelihoodCallInfo*>(
                rec.call_info.get());
            if (info == nullptr || rl_caller == nullptr) {
                continue;
            }
            // `converted` belongs to this iteration; nothing may point into it afterwards.
            PhaseReadEvidence converted;
            const PhaseReadEvidence* pe = info->read_phasing_evidence(converted);
            if (pe == nullptr) {
                continue;
            }
            // This site's own term, or none: a homozygote has no PhaseSite and contributed nothing to
            // Lambda, so there is nothing to subtract, and it can be corrected into a heterozygote.
            own.clear();
            auto found_site = site_by_key.find(rec.record_key);
            if (found_site != site_by_key.end()) {
                site_own_log_odds(*found_site->second, phase_flips.count(rec.record_key) != 0,
                                  own);
            }

            // The likelihoods before correction, copied only when something reads them.
            map<vector<int>, double> before;
            if (want_ledger) {
                before = info->genotype_lls;
            }
            // At --regeno-passes 1 the correction is computed and reported, and nothing is kept.
            // `genotype_lls` is what GL and QUAL are written from, so correcting it in place would
            // change them while the genotypes stood still.
            const bool keep = regenotype_passes >= 2;
            map<vector<int>, double> scratch;
            if (keep) {
                // Correct the direct pass's likelihoods every round, not the previous round's: the first
                // round saves them, and later rounds restore them before correcting. GQ is
                // restored with them, so that a GQ recomputed in an earlier round does not outlive
                // the correction it was computed from.
                if (info->uncorrected_lls == nullptr) {
                    info->uncorrected_lls.reset(
                        new map<vector<int>, double>(info->genotype_lls));
                } else {
                    info->genotype_lls = *info->uncorrected_lls;
                    rl_caller->recompute_gq(*info);
                }
            } else {
                scratch = info->genotype_lls;
            }
            map<vector<int>, double>& target = keep ? info->genotype_lls : scratch;
            const auto ps = site_phase_set.find(rec.record_key);
            const size_t phase_set = ps != site_phase_set.end() ? ps->second : NO_PHASE_SET;
            const auto s0 = site_strand0.find(rec.record_key);
            const int strand0_allele = s0 != site_strand0.end() ? s0->second : -1;
            const auto hap = haploid_strand.find(rec.record_key);
            const bool site_moved =
                hap != haploid_strand.end()
                    ? haploid_inclusion_correction(*pe, lambda, phase_set, own, temper, ceiling,
                                                   hap->second, regenotype_params, target,
                                                   counters)
                    : phase_aware_correction(*pe, lambda, phase_set, strand0_allele, own, temper,
                                             ceiling, regenotype_params, target, counters);
            if (keep && site_moved) {
                // The correction changed the best genotype, so GQ is recomputed from the corrected
                // likelihoods, as the direct pass computes it. GQI and GQN are not: GQN's achievable gap
                // assumes the site's own mixture weights, not per-read ones.
                rl_caller->recompute_gq(*info);
            }
            // Both ploidies, as in the direct pass, since the linkage pass can move a chain from
            // ploidy 1 to 2. At ploidy 1 the correction is zero, but it is applied the same way.
            if (keep && info->alt_ploidy_info != nullptr) {
                auto& alt = *info->alt_ploidy_info;
                if (alt.uncorrected_lls == nullptr) {
                    alt.uncorrected_lls.reset(new map<vector<int>, double>(alt.genotype_lls));
                } else {
                    alt.genotype_lls = *alt.uncorrected_lls;
                    rl_caller->recompute_gq(alt);
                }
                RegenotypeCounters ignored;
                if (phase_aware_correction(*pe, lambda, phase_set, strand0_allele, own, temper,
                                           ceiling, regenotype_params, alt.genotype_lls,
                                           ignored)) {
                    rl_caller->recompute_gq(alt);
                }
            }
            if (!site_moved) {
                continue;
            }
            ++moved;
            if (want_ledger) {
                auto best_of = [](const map<vector<int>, double>& gl) {
                    const vector<int>* b = nullptr;
                    double v = -numeric_limits<double>::infinity();
                    for (const auto& kv : gl) {
                        if (kv.second > v) { v = kv.second; b = &kv.first; }
                    }
                    return std::make_pair(b, v);
                };
                auto spell = [](const vector<int>* g) {
                    string out;
                    if (g == nullptr) {
                        return string(".");
                    }
                    for (size_t i = 0; i < g->size(); ++i) {
                        out += (i ? "/" : "") + std::to_string((*g)[i]);
                    }
                    return out;
                };
                const auto a = best_of(before);
                const auto b = best_of(target);
                const string snarl_id = print_snarl(rec.snarl);
                std::ostringstream row;
                row << snarl_id << "\t" << rec.ref_path_name << "\t"
                    << rec.ref_offset << "\t" << rec.ploidy << "\t" << spell(a.first) << "\t"
                    << spell(b.first) << "\t" << (b.second - a.second) << "\t"
                    << pe->num_reads();
                thread_ledger[qi].push_back(
                    LedgerRow{rec.ref_path_name, (size_t)rec.ref_offset, snarl_id, row.str()});
            }
        }
        thread_moved[qi] = moved;
    }
    size_t moved = 0;
    vector<LedgerRow> rows;
    for (size_t qi = 0; qi < n_queues; ++qi) {
        merge_counters(thread_counters[qi], regenotype_counters);
        moved += thread_moved[qi];
        std::move(thread_ledger[qi].begin(), thread_ledger[qi].end(), std::back_inserter(rows));
    }
    if (ledger.is_open()) {
        // The snarl ID breaks ties between records at the same position.
        std::sort(rows.begin(), rows.end(), [](const LedgerRow& x, const LedgerRow& y) {
            if (x.contig != y.contig) return x.contig < y.contig;
            if (x.position != y.position) return x.position < y.position;
            return x.snarl < y.snarl;
        });
        for (const LedgerRow& row : rows) {
            ledger << row.text << endl;
        }
        ledger.close();
    }

    const RegenotypeCounters& c = regenotype_counters;
    cerr << "[vg call] re-genotyping: temper " << temper << ", ceiling " << ceiling << ", "
         << c.reads_with_lambda
         << " reads carry a strand log-odds (" << c.reads_singleton
         << " span one site, so are inert; " << c.reads_multi_phase_set << " span two blocks), "
         << c.sites_corrected << " of " << c.sites_considered << " sites corrected, "
         << c.sites_would_move << " would move (" << c.moved_hom_to_het << " hom->het, "
         << c.moved_het_to_hom << " het->hom, " << c.moved_het_to_het << " het->het), "
         << c.order_reversed << " where the reads prefer the other order" << endl;
    if (regenotype_params.haploid_include) {
        cerr << "[vg call] re-genotyping: " << c.haploid_sites
             << " nested haploid chains weighted by whether the reads belong to their strand, "
             << c.haploid_would_move << " would move" << endl;
    }
    if (show_progress && !c.fit_count.empty()) {
        // Only under --progress: the calibration table is a diagnostic.
        cerr << "[vg call] re-genotyping calibration, |Lambda| / observed / predicted / n:";
        for (size_t i = 0; i < c.fit_count.size(); ++i) {
            cerr << "  " << c.fit_abs_lambda[i] << " " << c.fit_observed[i] << " "
                 << c.fit_predicted[i] << " " << c.fit_count[i];
        }
        cerr << endl;
    }
    return moved > 0;
}

void FlowCaller::rerun_linkage_pass() {
    if (!linker.enabled()) {
        return;
    }
    // Give the linkage model the corrected likelihoods, then run the linkage pass again in full, so that
    // every child is reassessed against its parent's new chosen pair, as on the first pass.
    linker.resync(staged_sites);
    run_linkage_pass();
}

void FlowCaller::phase_and_regenotype() {
    // Round 1 ends with read phasing; its linkage pass has already run.
    apply_read_phasing();
    // Each later round re-genotypes from the phase: the current phase gives every read its strand
    // log-odds, the correction rescores every site from the direct pass's likelihoods, the linkage
    // pass chooses the genotypes from the result and reassesses every nested child, and read
    // phasing runs again on the new genotypes. Rounds stop when the correction moves no site's
    // direct call, or when the chosen genotypes stop changing, return to an earlier round's, or
    // reach --regeno-passes rounds. With --regeno-passes 1 the correction is only computed and
    // reported.
    if (regenotype && regenotype_passes >= 2) {
        // Every state the rounds have reached, so that a cycle is recognised. The rounds can cycle:
        // dropping and reinstating a subtree is a discrete change, and the phase is a chain whose
        // links move with the genotypes, so no single quantity must increase.
        vector<size_t> seen_states;
        for (size_t round = 2; round <= regenotype_passes; ++round) {
            const auto before = chosen_snapshot();
            if (round == 2) {
                // Round 1's genotypes.
                seen_states.push_back(snapshot_digest(before));
            }
            const bool calls_moved = apply_regenotyping();
            // Chosen even when the correction moved no direct call: it has already changed every
            // site's likelihoods in place, and GL is written from them, so the genotypes are
            // chosen from them too.
            rerun_linkage_pass();
            apply_read_phasing();
            const auto after = chosen_snapshot();
            const size_t moved = chosen_changed(before, after);
            cerr << "[vg call] re-genotyping round " << round << ": " << moved
                 << " chosen genotypes moved" << endl;
            if (!calls_moved) {
                cerr << "[vg call] re-genotyping: the correction moved no site's direct call;"
                     << " stopping after round " << round << endl;
                break;
            }
            if (moved == 0) {
                cerr << "[vg call] re-genotyping: converged after " << round << " rounds" << endl;
                break;
            }
            const size_t digest = snapshot_digest(after);
            for (size_t i = 0; i < seen_states.size(); ++i) {
                if (seen_states[i] == digest) {
                    cerr << "[vg call] re-genotyping: LIMIT CYCLE of period "
                         << (seen_states.size() - i) << ", entered at round " << (i + 1)
                         << ". The iteration does not converge and no round of a cycle is more"
                         << " the answer than another; stopping here and reporting it rather than"
                         << " presenting round " << round << " as a fixed point" << endl;
                    goto regeno_done;
                }
            }
            seen_states.push_back(digest);
            if (round == regenotype_passes) {
                if (regenotype_passes == 2) {
                    // The default number of passes.
                    cerr << "[vg call] re-genotyping: one correction round applied; the iteration"
                         << " was not run further (--regeno-passes)" << endl;
                } else {
                    cerr << "[vg call] re-genotyping: NOT CONVERGED and no repeated state seen --"
                         << " still moving " << moved << " genotypes at the round cap of "
                         << regenotype_passes << ". A cycle longer than the rounds run"
                         << " cannot be detected, so raise the cap before concluding there is"
                         << " none" << endl;
                }
            }
        }
    regeno_done:;
    } else {
        // One round: compute and report the correction, and keep nothing.
        apply_regenotyping();
    }
}

void FlowCaller::render_retained_records() {
    // Each read's strand log-odds, for the anchors collected during the render and the hand-off.
    // Read phasing is done, so `read_strands` is final here.
    read_strands.build_lambda(
        phase_table.calls(),
        [&](const string& contig, size_t phase_set) { return phase_set_id(contig, phase_set); },
        regenotype_counters.fitted_temper, regenotype_counters.fitted_ceiling, regenotype_params);
    // The phase, before any record is built, so that each record is phased as it is rendered. Also
    // before the hand-off, which collects anchors for the records that get no line
    // (`reported_inline` and `no_reference`) and reads the frozen phase to order them. If read
    // phasing ran, the phase table already carries its swaps.
    phase_table.freeze_for_render(emit_phasing);
    // Every linkage pass is done, so the records move to the render, once, which also keeps their
    // anchors from being collected twice.
    hand_off_deferred_records();
    if (!staged_sites.active()) {
        return;
    }
    // `nested_context` describes the snarl a direct pass thread is recording, and only
    // `call_snarl_internal` reads it, which the render does not call. The loop still clears it and
    // restores it afterwards, so that it never runs under the context the thread's last swept snarl
    // left. The records are nested chains as well as top-level sites,
    // and each carries its own nesting in its `StagedSite`.
    const size_t n_threads = staged_sites.queue_count();
#pragma omp parallel for schedule(dynamic, 1)
    for (size_t t = 0; t < n_threads; ++t) {
        NestingPlacement saved_ctx = nested_context;
        nested_context = NestingPlacement();
        for (StagedSite& rec : staged_sites.queue(t)) {
            // The chosen pair, not the direct pass's. The ALT list, whether a line is written at
            // all, QUAL, and the arity of AD, GL and GQI are all built from the genotype passed in,
            // so they agree with the call.
            vector<int> genotype = linker.chosen_genotype(rec);
            // Before emit_variant, which passes the CallInfo on to update_vcf_info. The anchors are
            // collected in phase order, while `genotype` itself stays sorted, since emit_variant
            // builds the ALT list, AD, GL and QUAL from its order.
            collect_anchors_for_record(rec, genotype);
            emit_variant(graph, snarl_caller, rec.snarl, rec.travs, genotype, rec.ref_trav_idx,
                         rec.call_info, rec.ref_path_name, rec.ref_offset, genotype_snarls,
                         rec.ploidy);
        }
        nested_context = saved_ctx;
    }
    if (show_progress) {
        cerr << "[vg call] rendered " << staged_sites.queued_count()
             << " retained records after the direct pass" << endl;
    }
}

void FlowCaller::run_linkage_pass() {
    if (!staged_sites.active()) {
        return;
    }
    const GenotypeLinker::PassCounts counts = linker.link(staged_sites, phase_table, emit_phasing);
    vector<StagedSite>& pending = staged_sites.nested();
    size_t pass_inline_rederived = 0;

    // The exactly-once test, from the chosen genotypes the render builds each parent's blocks
    // from: `reported_inline` holds back a chain's line where an enclosing block's ALT spells it.
    // Without the linkage model every chosen genotype is the direct call the direct pass tested, so
    // there is nothing to redo.
    if (linker.enabled() && staged_sites.has_children()) {
        const unordered_map<size_t, StagedSite*> record_by_key = staged_sites.by_key();
        // Counted again from here, so that the report gives the chains held back now.
        block_records.restart_inline_count();
        // Parents before their children, so that a chain inherits its parent's final flag.
        staged_sites.for_each_parent_top_down(record_by_key, [&](const StagedSite& parent,
                                                                 const vector<size_t>& children) {
            if (parent.dropped) {
                return;   // its children were dropped with it
            }
            // The parts of the test that do not depend on the child, built once for this parent;
            // see BlockRecordWriter::ChainInlineContext.
            const BlockRecordWriter::ChainInlineContext ctx = block_records.chain_inline_context(
                parent.snarl, parent.travs, linker.chosen_genotype(parent), parent.ref_trav_idx);
            for (size_t ci : children) {
                StagedSite& child = pending[ci];
                if (child.dropped) {
                    continue;
                }
                const bool was = child.reported_inline;
                child.reported_inline = parent.reported_inline
                                        || block_records.chain_reported_inline(ctx, child.snarl);
                if (was != child.reported_inline) {
                    ++pass_inline_rederived;
                }
            }
        });
    }
    if (show_progress) {
        // The bytes kept for the staged sites, counted by walking the objects. They are walked in
        // parallel; the totals are sums, so they do not depend on how the walk is split.
        size_t retained_bytes = 0, retained_visits = 0, retained_gls = 0;
        auto measure = [](const StagedSite& rec, size_t& bytes, size_t& visits, size_t& gls) {
            bytes += sizeof(StagedSite) + rec.ref_path_name.capacity()
                     + rec.genotype.capacity() * sizeof(int)
                     + rec.panel_cache.capacity() * sizeof(int);
            bytes += rec.travs.capacity() * sizeof(SnarlTraversal);
            for (const SnarlTraversal& t : rec.travs) {
                visits += (size_t)t.visit_size();
                bytes += (size_t)t.visit_size() * sizeof(Visit);
            }
            const auto* rl = dynamic_cast<const ReadLikelihoodSnarlCaller::ReadLikelihoodCallInfo*>(
                rec.call_info.get());
            if (rl != nullptr) {
                for (const auto& kv : rl->genotype_lls) {
                    ++gls;
                    bytes += 48 + kv.first.capacity() * sizeof(int) + sizeof(double);
                }
                if (rl->anchor_evidence != nullptr) {
                    bytes += rl->anchor_evidence->bytes();
                }
                if (rl->phase_evidence != nullptr) {
                    bytes += rl->phase_evidence->bytes();
                }
                // The parts re-genotyping adds.
                auto gl_bytes = [](const map<vector<int>, double>& gl) {
                    size_t n = 0;
                    for (const auto& kv : gl) {
                        n += 48 + kv.first.capacity() * sizeof(int) + sizeof(double);
                    }
                    return n;
                };
                if (rl->uncorrected_lls != nullptr) {
                    bytes += gl_bytes(*rl->uncorrected_lls);
                }
                bytes += rl->scored_traversals.capacity() * sizeof(SnarlTraversal)
                         + rl->allele_support.capacity() * sizeof(double);
                if (rl->alt_ploidy_info != nullptr) {
                    // The alternate answer is kept too, with all its parts.
                    const auto& alt = *rl->alt_ploidy_info;
                    bytes += alt.scored_traversals.capacity() * sizeof(SnarlTraversal)
                             + alt.allele_support.capacity() * sizeof(double);
                    if (alt.uncorrected_lls != nullptr) {
                        bytes += gl_bytes(*alt.uncorrected_lls);
                    }
                }

                if (rl->alt_ploidy_info != nullptr) {
                    for (const auto& kv : rl->alt_ploidy_info->genotype_lls) {
                        ++gls;
                        bytes += 48 + kv.first.capacity() * sizeof(int) + sizeof(double);
                    }
                }
            }
        };
#pragma omp parallel reduction(+ : retained_bytes, retained_visits, retained_gls)
        {
            for (size_t q = 0; q < staged_sites.queue_count(); ++q) {
                const vector<StagedSite>& queue = staged_sites.queue(q);
#pragma omp for schedule(dynamic, 4096) nowait
                for (size_t r = 0; r < queue.size(); ++r) {
                    measure(queue[r], retained_bytes, retained_visits, retained_gls);
                }
            }
#pragma omp for schedule(dynamic, 4096) nowait
            for (size_t r = 0; r < pending.size(); ++r) {
                measure(pending[r], retained_bytes, retained_visits, retained_gls);
            }
        }
        // The read-phasing evidence. In the report below, the snarls the linkage pass will not revise
        // are the top-level ones and the children RecurseOnFail reaches without a ploidy override.
        for (const PhaseSite& ps : read_strands.sites()) {
            retained_bytes += sizeof(PhaseSite) + ps.read_key.capacity() * sizeof(uint64_t)
                              + ps.q0.capacity() * sizeof(float) + ps.c.capacity() * sizeof(float);
        }
        retained_bytes += read_strands.flips().size() * (sizeof(size_t) + 16);
        cerr << "[vg call] retained for rendering: " << staged_sites.queued_count()
             << " snarls the linkage pass will not revise, plus " << pending.size()
             << " nested chains; " << (retained_bytes / (1024.0 * 1024.0)) << " MB over "
             << retained_visits << " traversal visits and " << retained_gls
             << " genotype likelihoods" << endl;
        cerr << "[vg call] linkage pass exits: " << counts.no_crossing
             << " dropped because no parent candidate crosses them, " << counts.no_chosen
             << " whose parent's chosen pair could not be read, " << counts.unrenderable
             << " unrenderable so left unrevised, " << counts.ploidy_unscored
             << " stranded at a ploidy the direct pass never scored" << endl;
        if (pass_inline_rederived > 0) {
            cerr << "[vg call] linkage pass: " << pass_inline_rederived
                 << " children whose exactly-once suppression changed with their parent's"
                 << " chosen genotype" << endl;
        }
        cerr << "[vg call] single direct pass: " << pending.size() << " nested chains retained over "
             << (counts.levels + 1) << " levels; " << counts.revised << " revised, "
             << counts.gained << " reachable only under the chosen parent, " << counts.retracted
             << " retracted";
        if (counts.crossing_unknown > 0) {
            cerr << ", " << counts.crossing_unknown
                 << " with a crossing mask the direct pass could not compute";
        }
        if (counts.unspecifiable > 0) {
            cerr << ", " << counts.unspecifiable << " dropped from the layer because the site's "
                 << "compact allele space could not be built";
        }
        cerr << endl;
    }
}

void FlowCaller::hand_off_deferred_records() {
    if (!staged_sites.active()) {
        return;
    }
    const vector<StagedSite>& pending = staged_sites.nested();
    // A chain that gets no line still gets anchors: one whose variation an enclosing block's ALT
    // already spells, and one with no reference path, so no REF or POS, whose anchors are placed by
    // node ID. Each chain's anchors are its own and the anchor writer sorts every anchor before
    // writing, so they are collected on several threads, before the chains are handed over.
#pragma omp parallel for schedule(dynamic, 256)
    for (size_t i = 0; i < pending.size(); ++i) {
        const StagedSite& pr = pending[i];
        if (!pr.dropped && (pr.reported_inline || pr.no_reference)) {
            collect_anchors_for_record(pr, linker.chosen_genotype(pr));
        }
    }
    // Hand every surviving chain to the render, so that nested and top-level records are written
    // in one place from their chosen genotypes. A dropped chain is not handed over, since the
    // sample has no copy of it. A chain an enclosing block's ALT spells, or one with no reference
    // path, gets no line, since the render calls emit_variant for every record it is given; its
    // anchors were collected above.
    const StagedSiteTable::HandOff held = staged_sites.hand_off();
    if (held.inline_unrendered > 0 && show_progress) {
        cerr << "[vg call] block emission: " << held.inline_unrendered
             << " chains genotyped, recorded and phased, but left unrendered because an enclosing"
             << " block's ALT already spells them out" << endl;
    }
    if (held.no_ref_unrendered > 0 && show_progress) {
        cerr << "[vg call] off-reference nested: " << held.no_ref_unrendered
             << " chains chosen and left unrendered, having no reference position to write" << endl;
    }
}

bool FlowCaller::call_snarl_internal(const Snarl& managed_snarl,
                                      const string& parent_ref_path_name,
                                      pair<size_t, size_t> parent_ref_interval,
                                      const ChildTraversalSets* parent_child_trav_sets,
                                    int ploidy_override) {


    // todo: In order to experiment with merging consecutive snarls to make longer traversals,
    // I am experimenting with sending "fake" snarls through this code.  So make a local
    // copy to work on to do things like flip -- calling any snarl_manager code that
    // wants a pointer will crash.
    Snarl snarl = managed_snarl;

    // Staged in the nested branch below and completed after descent, which reads `travs`, since
    // this record then takes ownership of them.
    unique_ptr<StagedSite> pending_this;
    // The panel alleles `record_site` looked up for this snarl, if it recorded the site. The staged
    // record keeps them, so that re-genotyping does not look them up again; the record's
    // traversals are this snarl's `travs`, which do not change after the site is recorded.
    vector<int> site_panel;
    bool site_panel_set = false;
    // The same, for a snarl the linkage pass will not revise.
    unique_ptr<StagedSite> render_this;

#ifdef debug
    cerr << "call_snarl_internal on " << pb2json(snarl) << " with parent_ref_path=" << parent_ref_path_name
         << " parent_child_trav_sets=" << (parent_child_trav_sets ? "provided" : "null") << endl;
#endif

    if (snarl.start().node_id() == snarl.end().node_id() ||
        !graph.has_node(snarl.start().node_id()) || !graph.has_node(snarl.end().node_id())) {
        // can't call one-node or out-of graph snarls.
        return false;
    }

    // toggle average flow / flow width based on snarl length.  this is a bit inconsistent with
    // downstream which uses the longest traversal length, but it's a bit chicken and egg
    // todo: maybe use snarl length for everything?
    //
    // Only the flow traversal finder uses greedy_avg_flow, so the sum is computed only when there
    // is one.
    const auto& support_finder = dynamic_cast<SupportBasedSnarlCaller&>(snarl_caller).get_support_finder();
    FlowTraversalFinder* flow_trav_finder = dynamic_cast<FlowTraversalFinder*>(&traversal_finder);
    bool greedy_avg_flow = false;
    {
        auto snarl_contents = snarl_manager.deep_contents(&snarl, graph, false);
        if (snarl_contents.second.size() > max_snarl_edges) {
            // size cap needed as non-nested FlowCaller doesn't handle large snarls
            return false;
        }
        if (flow_trav_finder != nullptr) {
            size_t len_threshold = support_finder.get_average_traversal_support_switch_threshold();
            size_t length = 0;
            for (auto i = snarl_contents.first.begin();
                 i != snarl_contents.first.end() && length < len_threshold; ++i) {
                length += graph.get_length(graph.get_handle(*i));
            }
            greedy_avg_flow = length > len_threshold;
        }
    }
    
    handle_t start_handle = graph.get_handle(snarl.start().node_id(), snarl.start().backward());
    handle_t end_handle = graph.get_handle(snarl.end().node_id(), snarl.end().backward());

    // as we're writing to VCF, we need a reference path through the snarl.  we
    // look it up directly from the graph, and abort if we can't find one
    set<string> start_path_names;
    graph.for_each_step_on_handle(start_handle, [&](step_handle_t step_handle) {
            string name = graph.get_path_name(graph.get_path_handle_of_step(step_handle));
            if (!Paths::is_alt(name) && (ref_path_set.empty() || ref_path_set.count(name))) {
                start_path_names.insert(name);
            }
            return true;
        });
    
    set<string> end_path_names;
    if (!start_path_names.empty()) {
        graph.for_each_step_on_handle(end_handle, [&](step_handle_t step_handle) {
                string name = graph.get_path_name(graph.get_path_handle_of_step(step_handle));
                if (!Paths::is_alt(name) && (ref_path_set.empty() || ref_path_set.count(name))) {                
                    end_path_names.insert(name);
                }
                return true;
            });
    }
    
    // we do the full intersection (instead of more quickly finding the first common path)
    // so that we always take the lexicographically lowest path, rather than depending
    // on the order of iteration which could change between implementations / runs.
    vector<string> common_names;
    std::set_intersection(start_path_names.begin(), start_path_names.end(),
                          end_path_names.begin(), end_path_names.end(),
                          std::back_inserter(common_names));

    if (common_names.empty()) {
        // No reference path through snarl
        // If we have parent context, we can still process using parent's ref path
        // This test and the use_parent_interval test below must agree: otherwise get_ref_interval
        // would be called with the parent's reference path, which does not visit this snarl's
        // boundary nodes, and would assert.
        if ((parent_child_trav_sets == nullptr && !nested_context.no_reference)
            || parent_ref_path_name.empty()) {
#ifdef debug
            cerr << "  -> returning false: no common ref path and no parent context" << endl;
#endif
            return false;
        }
#ifdef debug
        cerr << "  -> using parent ref path: " << parent_ref_path_name << endl;
#endif
    }

    // Use parent's ref path if no direct path, otherwise prefer base reference over gref paths
    string ref_path_name;
    if (common_names.empty()) {
        ref_path_name = parent_ref_path_name;
    } else {
        // Prefer base reference paths over derived gref paths.  Test the whole gref
        // namespace, not just the fragment suffix: a gref copy of the reference sorts
        // before the path it was copied from (gref_x < x).
        // common_names is sorted, so we iterate to find first non-gref path
        ref_path_name = common_names.front();  // default to first (lexicographically smallest)
        for (const string& name : common_names) {
            if (!GrefCover::is_gref_derived(name)) {
                ref_path_name = name;
                break;
            }
        }
    }

    // find the reference traversal and coordinates using the path position graph interface
    tuple<int64_t, int64_t, bool, step_handle_t, step_handle_t> ref_interval;
    bool use_parent_interval = false;

    if (common_names.empty()) {
        // No direct reference path - use parent's interval and traversals directly
        ref_interval = make_tuple(parent_ref_interval.first, parent_ref_interval.second, false, step_handle_t(), step_handle_t());
        use_parent_interval = true;
    } else {
        ref_interval = get_ref_interval(graph, snarl, ref_path_name);
        if (get<0>(ref_interval) == -1) {
            // could not find reference path interval consistent with snarl due to orientation conflict
            return false;
        }
        if (get<2>(ref_interval) == true) {
            // calling code assumes snarl forward on reference
            flip_snarl(snarl);
            ref_interval = get_ref_interval(graph, snarl, ref_path_name);
        }
    }

    SnarlTraversal ref_trav;

    if (!use_parent_interval) {
        // Build reference traversal from path steps
        step_handle_t cur_step = get<3>(ref_interval);
        step_handle_t last_step = get<4>(ref_interval);
        if (get<2>(ref_interval)) {
            std::swap(cur_step, last_step);
        }
        bool start_backwards = snarl.start().backward() != graph.get_is_reverse(graph.get_handle_of_step(cur_step));

        while (true) {
            handle_t cur_handle = graph.get_handle_of_step(cur_step);
            Visit* visit = ref_trav.add_visit();
            visit->set_node_id(graph.get_id(cur_handle));
            visit->set_backward(start_backwards ? !graph.get_is_reverse(cur_handle) : graph.get_is_reverse(cur_handle));
            if (graph.get_id(cur_handle) == snarl.end().node_id()) {
                break;
            } else if (get<2>(ref_interval) == true) {
                if (!graph.has_previous_step(cur_step)) {
                    cerr << "Warning [vg call]: Unable, due to bug or corrupt path information, to trace reference path through snarl " << pb2json(managed_snarl) << endl;
                    return false;
                }
                cur_step = graph.get_previous_step(cur_step);
            } else {
                if (!graph.has_next_step(cur_step)) {
                    cerr << "Warning [vg call]: Unable, due to bug or corrupt path information, to trace reference path through snarl " << pb2json(managed_snarl) << endl;
                    return false;
                }
                cur_step = graph.get_next_step(cur_step);
            }
            // todo: we can compute flow at the same time
        }
        assert(ref_trav.visit(0) == snarl.start() && ref_trav.visit(ref_trav.visit_size() - 1) == snarl.end());
    }
    // If use_parent_interval, ref_trav stays empty - we'll use first parent traversal as pseudo-reference

    vector<SnarlTraversal> travs;
    if (flow_trav_finder != nullptr) {
        // find the max flow traversals using specialized interface that accepts avg heurstic toggle
        pair<vector<SnarlTraversal>, vector<double>> weighted_travs = flow_trav_finder->find_weighted_traversals(snarl, greedy_avg_flow);
        travs = std::move(weighted_travs.first);
    } else {
        // find the traversals using the generic interface
        travs = traversal_finder.find_traversals(snarl);
    }

    if (travs.empty()) {
        cerr << "Warning [vg call]: Unable, due to bug or corrupt graph, to search for any traversals through snarl " << pb2json(managed_snarl) << endl;
        return false;
    }
#ifdef debug
    cerr << "  found " << travs.size() << " traversals, use_parent_interval=" << use_parent_interval << endl;
#endif

    // optional traversal length clamp can, ex, avoid trying to resolve a giant snarl    
    if (allele_length_range.first > 0 || allele_length_range.second < numeric_limits<size_t>::max()) {
        size_t max_trav_len = 0;
        for (const SnarlTraversal & trav : travs) {
            size_t trav_len = 0;
            for (size_t i = 1; i < trav.visit_size() - 1; ++i) {
                trav_len += graph.get_length(graph.get_handle(trav.visit(i).node_id()));
            }
            max_trav_len = max(max_trav_len, trav_len);
            if (max_trav_len > allele_length_range.second) {
                return false;
            }
        }
        if (max_trav_len < allele_length_range.first) {
            return false;
        }
    }

    // find the reference traversal in the list of results from the traversal finder
    int ref_trav_idx = -1;

    if (use_parent_interval) {
        // No direct reference path - use first traversal from first non-empty set as pseudo-reference
        if (parent_child_trav_sets != nullptr) {
            for (const auto& tset : *parent_child_trav_sets) {
                if (!tset.empty()) {
                    const SnarlTraversal& first_trav = tset[0];
                    for (int i = 0; i < travs.size() && ref_trav_idx < 0; ++i) {
                        if (travs[i] == first_trav) {
                            ref_trav_idx = i;
                        }
                    }
                    if (ref_trav_idx < 0 && first_trav.visit_size() > 0) {
                        ref_trav_idx = travs.size();
                        travs.push_back(first_trav);
                    }
                    break;
                }
            }
        }
        // Left at -1 where no parent traversal set named a reference: travs[0] is the flow finder's
        // best-supported traversal, and making it REF would bias the genotyper's tie-break toward
        // it. The genotyper checks ref_trav_idx >= 0 before using it.
        if (ref_trav_idx < 0 && parent_child_trav_sets != nullptr) {
            ref_trav_idx = travs.empty() ? -1 : 0;
        }
    } else {
        for (int i = 0; i < travs.size() && ref_trav_idx < 0; ++i) {
            // todo: is there a way to speed this up?
            if (travs[i] == ref_trav) {
                ref_trav_idx = i;
            }
        }

        if (ref_trav_idx == -1) {
            ref_trav_idx = travs.size();
            // we didn't get the reference traversal from the finder, so we add it here
            travs.push_back(ref_trav);
        }
    }

    bool ret_val = true;
    vector<int> trav_genotype;  // Declared outside block so we can pass to children

    // A ploidy from the parent overrides the contig's or the region BED's: it is the number of the
    // parent's called alleles that reach this child. The region's is still the number of the
    // sample's haplotypes here, which the depth term needs.
    const int region_ploidy = ploidy_regions.ploidy_at(ref_path_name, get<0>(ref_interval),
                                                       ref_offset_of(ref_offsets, ref_path_name),
                                                       ref_ploidy_of(ref_ploidies, ref_path_name));
    int ploidy = ploidy_override >= 0 ? ploidy_override : region_ploidy;

    // What both the parent-traversal-set branch and the top-level branch do with their genotype.
    // `trav_call_info` differs between them, so it is a parameter. `snarl` is captured by reference;
    // `flip_snarl` may already have rewritten it above.
    auto stage_or_emit = [&](unique_ptr<SnarlCaller::CallInfo>& trav_call_info) -> bool {
        bool added;
        if (!gaf_output) {
            // Staged, not emitted: `render_retained_records` writes it after the direct pass.
            // `added` stands in for emit_variant's return value, which here only gates recursion;
            // a staged site counts as added.
            site_panel_set = linker.add(snarl, travs, trav_genotype, trav_call_info.get(),
                                        ref_trav_idx, ref_path_name,
                                        ref_offset_of(ref_offsets, ref_path_name),
                                        record_key_of(snarl), nested_context, false, 0,
                                        &site_panel);
            render_this = stage_render_record(snarl, trav_genotype, ref_trav_idx, trav_call_info,
                                              ref_path_name, ref_offset_of(ref_offsets, ref_path_name), ploidy);
            added = render_this != nullptr;
            if (!added) {
                added = emit_variant(graph, snarl_caller, snarl, travs, trav_genotype, ref_trav_idx,
                                     trav_call_info, ref_path_name, ref_offset_of(ref_offsets, ref_path_name),
                                     genotype_snarls, ploidy);
            }
        } else {
            added = true;
            pair<string, int64_t> pos_info = get_ref_position(graph, snarl, ref_path_name, ref_offset_of(ref_offsets, ref_path_name));
            emit_gaf_variant(graph, print_snarl(snarl), travs, trav_genotype, ref_trav_idx, pos_info.first, pos_info.second, &support_finder);
        }
        return added;
    };


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
            ReadLikelihoodSnarlCaller::set_region_ploidy(region_ploidy);
            std::tie(called_alleles, trav_call_info) = snarl_caller.genotype(
                snarl, travs, ref_trav_idx, effective_ploidy, ref_path_name,
                make_pair(get<0>(ref_interval), get<1>(ref_interval)));
            ReadLikelihoodSnarlCaller::set_region_ploidy(0);

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
        } else {
            added = stage_or_emit(trav_call_info);
        }

        ret_val = trav_genotype.size() == ploidy && added;
    } else if (ploidy_override >= 0) {
        // A nested chain, reached by descent, at the ploidy its parent implied. Only a nested chain
        // can have its ploidy revised at the linkage pass, so only it needs the other ploidy's answer.
        unique_ptr<SnarlCaller::CallInfo> trav_call_info;
        ReadLikelihoodSnarlCaller::set_want_alt_ploidy(true);
        ReadLikelihoodSnarlCaller::set_region_ploidy(region_ploidy);
        std::tie(trav_genotype, trav_call_info) = snarl_caller.genotype(
            snarl, travs, ref_trav_idx, ploidy, ref_path_name,
            make_pair(get<0>(ref_interval), get<1>(ref_interval)));
        ReadLikelihoodSnarlCaller::set_region_ploidy(0);
        ReadLikelihoodSnarlCaller::set_want_alt_ploidy(false);

        const bool retain_only = nested_context.retain_only;
        // Whether this snarl's own boundaries are on no reference path, checked from the graph for
        // each snarl.
        const bool no_ref_position = use_parent_interval;

        assert(trav_genotype.empty() || trav_genotype.size() == ploidy);
        bool added = true;
        if (no_ref_position) {
            // Genotyped and recorded, never written. Checked before retain_only, which does not
            // record. `added` is true, as for retain_only, since it gates descent into this chain's
            // children.
            site_panel_set = linker.add(
                snarl, travs, trav_genotype, trav_call_info.get(), ref_trav_idx, ref_path_name,
                ref_offset_of(ref_offsets, ref_path_name), record_key_of(snarl), nested_context,
                /*no_reference*/ true,
                // The parent's position, as `get_ref_position` gives it from the interval
                // `use_parent_interval` set, plus the chain's offset along its parent, as
                // `StagedSite::position_from_parent` has it.
                base_path_position(ref_path_name, get<0>(ref_interval)
                                                      + ref_offset_of(ref_offsets, ref_path_name))
                    + (int64_t)nested_context.parent_offset,
                &site_panel);
            ++descent_counters.no_ref_recorded;
            {
                int copies = 0;
                for (int a : trav_genotype) {
                    copies += (a >= 0);
                }
                descent_counters.no_ref_copies[copies < 3 ? copies : 2].fetch_add(1);
            }
            added = true;
        } else if (retain_only) {
            // No called parent allele reaches this chain, so nothing about it is written yet. It is
            // genotyped and kept, since the linkage model may still move the parent onto an allele
            // that reaches it.
            added = true;
        } else if (nested_context.reported_inline) {
            // An enclosing block's ALT already spells this chain, so it gets no line, but it is
            // genotyped and recorded, since its allele pair phases everything inside it. Checked
            // after retain_only, which does not record.
            site_panel_set = linker.add(snarl, travs, trav_genotype, trav_call_info.get(),
                                        ref_trav_idx, ref_path_name,
                                        ref_offset_of(ref_offsets, ref_path_name),
                                        record_key_of(snarl), nested_context, false, 0,
                                        &site_panel);
            added = true;
        } else if (!gaf_output) {
            // Recorded here rather than in emit_variant. A retained chain, on the path above, is
            // recorded only if the linkage pass later finds that the sample carries it.
            site_panel_set = linker.add(snarl, travs, trav_genotype, trav_call_info.get(),
                                        ref_trav_idx, ref_path_name,
                                        ref_offset_of(ref_offsets, ref_path_name),
                                        record_key_of(snarl), nested_context, false, 0,
                                        &site_panel);
            // Staged, not emitted, as at top level: the line is written after the linkage pass, from the
            // chosen genotype. `added` stands in for emit_variant's return value, which here only
            // gates recursion; a staged site counts as added.
            added = staged_sites.active();
            if (!added) {
                added = emit_variant(graph, snarl_caller, snarl, travs, trav_genotype, ref_trav_idx,
                                     trav_call_info, ref_path_name, ref_offset_of(ref_offsets, ref_path_name),
                                     genotype_snarls, ploidy);
            }
        } else {
            pair<string, int64_t> pos_info = get_ref_position(graph, snarl, ref_path_name,
                                                              ref_offset_of(ref_offsets, ref_path_name));
            emit_gaf_variant(graph, print_snarl(snarl), travs, trav_genotype, ref_trav_idx,
                             pos_info.first, pos_info.second, &support_finder);
        }

        // Stage the nested site without its traversals: descent below still reads `travs` to find
        // which children the called alleles reach, and they are moved in once descent is done.
        if (staged_sites.active()) {
            pending_this.reset(new StagedSite());
            pending_this->snarl = snarl;
            pending_this->ref_path_name = ref_path_name;
            pending_this->ref_offset = ref_offset_of(ref_offsets, ref_path_name);
            pending_this->ref_trav_idx = ref_trav_idx;
            pending_this->genotype = trav_genotype;
            pending_this->ploidy = ploidy;
            pending_this->record_key = record_key_of(snarl);
            pending_this->parent_record_key = nested_context.parent_record_key;
            pending_this->parent_crossing = nested_context.parent_crossing;
            pending_this->chain_key = nested_context.chain_key;
            pending_this->no_reference = no_ref_position;
            pending_this->reported_inline = nested_context.reported_inline;
            pending_this->position_from_parent =
                no_ref_position
                    ? base_path_position(ref_path_name,
                                         get<0>(ref_interval) + ref_offset_of(ref_offsets, ref_path_name))
                          + (int64_t)nested_context.parent_offset
                    : 0;
            pending_this->chain_offset = nested_context.parent_offset;
            pending_this->crossing_known = nested_context.crossing_known;
            pending_this->level = (uint8_t)min(nested_context.level, (size_t)255);
            pending_this->call_info = std::move(trav_call_info);
        }
        ret_val = trav_genotype.size() == ploidy && added;
    } else {
        // Top-level snarl or no parent context - genotype from scratch using support
        unique_ptr<SnarlCaller::CallInfo> trav_call_info;
        std::tie(trav_genotype, trav_call_info) = snarl_caller.genotype(snarl, travs, ref_trav_idx, ploidy, ref_path_name,
                                                                        make_pair(get<0>(ref_interval), get<1>(ref_interval)));

        assert(trav_genotype.empty() || trav_genotype.size() == ploidy);

        bool added = true;
        added = stage_or_emit(trav_call_info);

        ret_val = trav_genotype.size() == ploidy && added;
    }

    // Nested calling: descend into each child the called alleles reach, at the ploidy they reach
    // it with.
    //
    // Descent does not depend on whether a line was written: a parent written as the reference
    // still has children to call. Children are genotyped independently, with no parent traversal
    // sets. Only a successful call descends, since a failed snarl has no genotype to take a child's
    // ploidy from. RecurseOnFail calls the children of a failed top-level snarl as top-level
    // snarls, but nothing does so for a failed nested snarl: its children are not called.
    if (ret_val && symbolic_manager != nullptr && !trav_genotype.empty() &&
        parent_child_trav_sets == nullptr) {
        const Snarl* managed_ptr = snarl_manager.into_which_snarl(snarl.start().node_id(),
                                                                  snarl.start().backward());
        if (managed_ptr != nullptr) {

            // The child-independent parts of the exactly-once test, built once for this snarl.
            const BlockRecordWriter::ChainInlineContext inline_ctx =
                block_records.chain_inline_context(snarl, travs, trav_genotype, ref_trav_idx);
            // Also once for this snarl: see ChildPlacer::TraversalNodeIndex.
            vector<ChildPlacer::TraversalNodeIndex> trav_visits;
            trav_visits.reserve(travs.size());
            for (const SnarlTraversal& t : travs) {
                trav_visits.push_back(ChildPlacer::index_traversal_nodes(t));
            }
            for (const Snarl* child : snarl_manager.children_of(managed_ptr)) {
                if (child == nullptr || snarl_manager.is_trivial(child, graph)) {
                    continue;
                }
                // A chain that no reference path passes through has no REF or POS for its records,
                // so it is skipped unless off-reference descent is on.
                bool child_off_reference = false;
                if (ref_trav_idx >= 0 && ref_trav_idx < (int)travs.size()) {
                    vector<int> ref_only(1, ref_trav_idx);
                    if (child_ploidy(trav_visits, ref_only, *child, 1) == 0) {
                        // With off-reference descent, such a chain is genotyped and recorded but has
                        // no line.
                        if (!off_reference_nesting) {
                            ++descent_counters.skipped_no_ref;
                            continue;
                        }
                        child_off_reference = true;
                        ++descent_counters.off_reference;
                    }
                }
                // Inherited: everything under a chain the reference does not cross is also off it.
                if (nested_context.no_reference) {
                    child_off_reference = true;
                }

                // The exactly-once test: under block emission, a chain that every called strand
                // crosses only inside a difference block is already spelled by that block's ALT. It
                // holds back the chain's line, not its descent, so the chain is still genotyped,
                // recorded and phased. Inherited by chains inside it. Does nothing when block
                // emission is off, or for a snarl whose projection has no symbols.
                bool child_reported_inline =
                    nested_context.reported_inline
                    || block_records.chain_reported_inline(inline_ctx, *child);

                int copies = child_ploidy(trav_visits, trav_genotype, *child, ploidy);
                bool retain_only = nested_context.retain_only;
                if (copies <= 0) {
                    // No called allele reaches it yet. Visited anyway, while this window's reads are
                    // in memory, since the linkage model may move the parent onto an allele that
                    // does reach it. Nothing about it is written unless the linkage pass says so.
                    ++descent_counters.skipped_no_copy;
                    if (!staged_sites.active() || !linker.enabled()) {
                        // Without retention there is nothing to come back to. Without the linkage
                        // model nothing moves the parent after the direct pass, so the sample has no copy
                        // of this chain; the linkage pass, which has no chosen parent to read, would
                        // otherwise render it at the parent's ploidy.
                        continue;
                    }
                    retain_only = true;
                }

                // Saved and restored, since a child may descend further, and its own children must see
                // it as their parent.
                NestingPlacement saved = nested_context;
                nested_context.one_copy = (copies == 1);
                nested_context.parent_record_key = record_key_of(snarl);
                nested_context.retain_only = retain_only;
                nested_context.no_reference = child_off_reference;
                // Where this child starts along the first called allele that reaches it, added to
                // the offset of its parent. Only an off-reference chain uses it, but it is computed
                // for every chain, so that offsets add up down the tree.
                nested_context.parent_offset =
                    saved.parent_offset + offset_along_genotype(travs, trav_genotype, *child);
                nested_context.reported_inline = child_reported_inline;
                // The chain's identity, from its boundary nodes.
                {
                    const pair<nid_t, nid_t> cb = chain_bounds_of(child, snarl_manager);
                    nested_context.chain_key =
                        (size_t)((uint64_t)cb.first * 1000003ULL) ^ (size_t)(uint64_t)cb.second;
                }
                bool crossing_known = true;   // child_crossing_mask always sets it
                // The mask is over this snarl's own candidate traversals, which exist whether or not
                // a line was written.
                nested_context.parent_crossing =
                    ChildPlacer::child_crossing_mask(trav_visits, *child, &crossing_known);
                nested_context.crossing_known = crossing_known;
                nested_context.level = saved.level + 1;
                ++g_descent_depth;
                if (g_descent_depth < 16) {
                    ++descent_counters.depth_hist[g_descent_depth];
                }
                // `copies` is zero only for a chain no called parent allele reaches, which is still
                // genotyped; it then takes the parent's ploidy, the most copies a child can have.
                // The other ploidy's answer is computed as well, so the linkage pass can change it
                // later.
                call_snarl_internal(*child, ref_path_name,
                                    make_pair(get<0>(ref_interval), get<1>(ref_interval)),
                                    nullptr, copies >= 1 ? copies : ploidy);
                --g_descent_depth;
                nested_context = saved;
            }
        }
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
                            TraversalSet tset = find_child_traversal_set(travs[allele_idx], *child);
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

    // Descent above and the --top-down recursion, which builds each child's ChildTraversalSets
    // from `travs`, are done, so the staged site can take the traversals. At most one of these
    // is set.
    if (pending_this != nullptr) {
        pending_this->travs = std::move(travs);
        if (site_panel_set) {
            pending_this->panel_cache = std::move(site_panel);
            pending_this->panel_cached = true;
        }
        staged_sites.add_nested(std::move(*pending_this));
        pending_this.reset();
    } else if (render_this != nullptr) {
        render_this->travs = std::move(travs);
        if (site_panel_set) {
            render_this->panel_cache = std::move(site_panel);
            render_this->panel_cached = true;
        }
        staged_sites.add_top_level(std::move(*render_this));
        render_this.reset();
    }


    return ret_val;
}

string FlowCaller::vcf_header(const PathHandleGraph& graph, const vector<string>& contigs,
                              const vector<size_t>& contig_length_overrides) const {
    string header = VCFOutputCaller::vcf_header(graph, contigs, contig_length_overrides);
    header += "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n";
    snarl_caller.update_vcf_header(header);
    header += "##FILTER=<ID=PASS,Description=\"All filters passed\">\n";
    header += "##SAMPLE=<ID=" + sample_name + ">\n";
    header += "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + sample_name;
    assert(output_vcf.openForOutput(header));
    header += "\n";
    return header;
}
}

