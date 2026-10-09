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
    install_widgets();

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
    install_widgets();
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
    record_steps.gl_layout = [this](const SnarlCaller::CallInfo* call_info) {
        // Every call info a run with the read-likelihood genotyper writes is that genotyper's.
        return site_genotyper != nullptr && call_info != nullptr ? GLLayout::Colexicographic
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

void FlowCaller::install_widgets() {
    const SiteReader reader{
        .graph = &graph,
        .genotyper = site_genotyper.get(),
        .spell = [this](const SnarlTraversal& trav) { return trav_string(graph, trav); },
        .name = [this](const Snarl& site) { return print_snarl(site); },
    };
    linker.set_site_reader(reader);
    rescorer.set_site_reader(reader);
    record_renderer.configure(reader, [this](const Snarl& site) { return snarl_is_leaf(site); });
    child_placer.configure(&graph, &snarl_manager, &block_records, &descent_counters);
}

void FlowCaller::call_top_level_snarls(const HandleGraph& graph, RecurseType recurse_type) {
    GraphCaller::call_top_level_snarls(graph, recurse_type);
    if (show_progress) {
        report_descent_instrumentation();
    }
}

void FlowCaller::set_site_genotyper(ReadLikelihoodSnarlCaller& genotyper) {
    site_genotyper.reset(new SiteGenotyper(genotyper));
    install_widgets();
}

pair<vector<int>, unique_ptr<SnarlCaller::CallInfo>> FlowCaller::genotype_site(
    const Snarl& site, const vector<SnarlTraversal>& travs, int ref_trav_idx,
    const Ploidies& ploidies, const string& ref_path_name, pair<size_t, size_t> ref_range,
    SiteScore*& score) {
    score = nullptr;
    if (site_genotyper == nullptr) {
        // Another genotyper, which takes the ploidy alone.
        return snarl_caller.genotype(site, travs, ref_trav_idx, ploidies.ploidy, ref_path_name,
                                     ref_range);
    }
    auto called =
        site_genotyper->genotype(site, travs, ref_trav_idx, ploidies, ref_path_name, ref_range);
    score = called.second.get();
    return make_pair(std::move(called.first),
                     unique_ptr<SnarlCaller::CallInfo>(std::move(called.second)));
}

bool FlowCaller::call_snarl(const Snarl& managed_snarl) {
    // Entry point: call with no parent context
    return call_snarl_internal(managed_snarl, "", make_pair(0, 0), nullptr, -1,
                               NestingPlacement());
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
        unique_ptr<SnarlCaller::CallInfo>& call_info, SiteScore* score,
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
    rec->set_call(std::move(call_info), score);
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
    const auto phase_set_of = [&](const string& contig, size_t phase_set) {
        return phase_set_id(contig, phase_set);
    };
    // Read phasing and re-genotyping change the linkage model's phase and the likelihoods it
    // chooses from, so neither runs without it.
    const auto phase_reads = [&]() {
        if (linker.enabled()) {
            read_phaser.phase(staged_sites, phase_table, read_strands, phase_set_of);
        }
    };
    const auto rescore = [&]() {
        return linker.enabled()
               && rescorer.rescore(staged_sites, phase_table, read_strands, temper_fit,
                                   phase_set_of, show_progress);
    };
    // Round 1 ends with read phasing; its linkage pass has already run.
    phase_reads();
    // Each later round re-genotypes from the phase: the current phase gives every read its strand
    // log-odds, the correction rescores every site from the direct pass's likelihoods, the linkage
    // pass chooses the genotypes from the result and reassesses every nested child, and read
    // phasing runs again on the new genotypes. Rounds stop when the correction moves no site's
    // direct call, or when the chosen genotypes stop changing, return to an earlier round's, or
    // reach --regeno-passes rounds. With --regeno-passes 1 the correction is only computed and
    // reported.
    if (rescorer.enabled() && rescorer.passes() >= 2) {
        // Every state the rounds have reached, so that a cycle is recognised.
        RoundHistory history;
        for (size_t round = 2; round <= rescorer.passes(); ++round) {
            const RoundHistory::State before = RoundHistory::state(staged_sites, linker.collector());
            if (round == 2) {
                // Round 1's genotypes.
                history.remember(before);
            }
            const bool calls_moved = rescore();
            // Chosen even when the correction moved no direct call: it has already changed every
            // site's likelihoods in place, and GL is written from them, so the genotypes are
            // chosen from them too.
            rerun_linkage_pass();
            phase_reads();
            const RoundHistory::State after = RoundHistory::state(staged_sites, linker.collector());
            const size_t moved = RoundHistory::changed(before, after);
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
            const size_t entered = history.first_round_with(after);
            if (entered != 0) {
                cerr << "[vg call] re-genotyping: LIMIT CYCLE of period "
                     << (history.rounds() - entered + 1) << ", entered at round " << entered
                     << ". The iteration does not converge and no round of a cycle is more"
                     << " the answer than another; stopping here and reporting it rather than"
                     << " presenting round " << round << " as a fixed point" << endl;
                break;
            }
            history.remember(after);
            if (round == rescorer.passes()) {
                if (rescorer.passes() == 2) {
                    // The default number of passes.
                    cerr << "[vg call] re-genotyping: one correction round applied; the iteration"
                         << " was not run further (--regeno-passes)" << endl;
                } else {
                    cerr << "[vg call] re-genotyping: NOT CONVERGED and no repeated state seen --"
                         << " still moving " << moved << " genotypes at the round cap of "
                         << rescorer.passes() << ". A cycle longer than the rounds run"
                         << " cannot be detected, so raise the cap before concluding there is"
                         << " none" << endl;
                }
            }
        }
    } else {
        // One round: compute and report the correction, and keep nothing.
        rescore();
    }
}

void FlowCaller::render_retained_records() {
    // Each read's strand log-odds, for the anchors collected during the render and the hand-off.
    // Read phasing is done, so `read_strands` is final here.
    read_strands.build_lambda(
        phase_table.calls(),
        [&](const string& contig, size_t phase_set) { return phase_set_id(contig, phase_set); },
        rescorer.last_counters().fitted_temper, rescorer.last_counters().fitted_ceiling,
        rescorer.params());
    // The phase, before any record is built, so that each record is phased as it is rendered. Also
    // before the hand-off, which collects anchors for the records that get no line
    // (`reported_inline` and `no_reference`) and reads the frozen phase to order them. If read
    // phasing ran, the phase table already carries its swaps.
    phase_table.freeze_for_render(emit_phasing);
    record_renderer.render(
        staged_sites, phase_table, read_strands, linker,
        anchor_collector.is_enabled() ? &anchor_collector : nullptr,
        [&](const StagedSite& site, const vector<int>& genotype) {
            emit_variant(graph, snarl_caller, site.snarl, site.travs, genotype, site.ref_trav_idx,
                         site.call_info, site.ref_path_name, site.ref_offset, genotype_snarls,
                         site.ploidy);
        },
        show_progress);
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
            const SiteScore* rl = rec.score;
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

bool FlowCaller::call_snarl_internal(const Snarl& managed_snarl,
                                      const string& parent_ref_path_name,
                                      pair<size_t, size_t> parent_ref_interval,
                                      const ChildTraversalSets* parent_child_trav_sets,
                                      int ploidy_override, const NestingPlacement& placement) {


    // todo: In order to experiment with merging consecutive snarls to make longer traversals,
    // I am experimenting with sending "fake" snarls through this code.  So make a local
    // copy to work on to do things like flip -- calling any snarl_manager code that
    // wants a pointer will crash.
    Snarl snarl = managed_snarl;

    // Staged in the nested branch below and completed after descent, which reads `travs`, since
    // this record then takes ownership of them.
    unique_ptr<StagedSite> pending_this;
    // The panel alleles `GenotypeLinker::add` looked up for this snarl, if it filed the site. The staged
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
        if ((parent_child_trav_sets == nullptr && !placement.no_reference)
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
    auto stage_or_emit = [&](unique_ptr<SnarlCaller::CallInfo>& trav_call_info,
                             SiteScore* trav_score) -> bool {
        bool added;
        if (!gaf_output) {
            // Staged, not emitted: `render_retained_records` writes it after the direct pass.
            // `added` stands in for emit_variant's return value, which here only gates recursion;
            // a staged site counts as added.
            site_panel_set = linker.add(snarl, travs, trav_genotype, trav_score,
                                        ref_trav_idx, ref_path_name,
                                        ref_offset_of(ref_offsets, ref_path_name),
                                        record_key_of(snarl), placement, false, 0,
                                        &site_panel);
            render_this = stage_render_record(snarl, trav_genotype, ref_trav_idx, trav_call_info,
                                              trav_score,
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
        SiteScore* trav_score = nullptr;
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
            std::tie(called_alleles, trav_call_info) = genotype_site(
                snarl, travs, ref_trav_idx,
                Ploidies{.ploidy = effective_ploidy, .region_ploidy = region_ploidy}, ref_path_name,
                make_pair(get<0>(ref_interval), get<1>(ref_interval)), trav_score);

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
            added = stage_or_emit(trav_call_info, trav_score);
        }

        ret_val = trav_genotype.size() == ploidy && added;
    } else if (ploidy_override >= 0) {
        // A nested chain, reached by descent, at the ploidy its parent implied. Only a nested chain
        // can have its ploidy revised at the linkage pass, so only it needs the other ploidy's answer.
        unique_ptr<SnarlCaller::CallInfo> trav_call_info;
        SiteScore* trav_score = nullptr;
        std::tie(trav_genotype, trav_call_info) = genotype_site(
            snarl, travs, ref_trav_idx,
            Ploidies{.ploidy = ploidy, .region_ploidy = region_ploidy, .also_score_other = true},
            ref_path_name, make_pair(get<0>(ref_interval), get<1>(ref_interval)), trav_score);

        const bool retain_only = placement.retain_only;
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
                snarl, travs, trav_genotype, trav_score, ref_trav_idx, ref_path_name,
                ref_offset_of(ref_offsets, ref_path_name), record_key_of(snarl), placement,
                /*no_reference*/ true,
                // The parent's position, as `get_ref_position` gives it from the interval
                // `use_parent_interval` set, plus the chain's offset along its parent, as
                // `StagedSite::position_from_parent` has it.
                base_path_position(ref_path_name, get<0>(ref_interval)
                                                      + ref_offset_of(ref_offsets, ref_path_name))
                    + (int64_t)placement.parent_offset,
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
        } else if (placement.reported_inline) {
            // An enclosing block's ALT already spells this chain, so it gets no line, but it is
            // genotyped and recorded, since its allele pair phases everything inside it. Checked
            // after retain_only, which does not record.
            site_panel_set = linker.add(snarl, travs, trav_genotype, trav_score,
                                        ref_trav_idx, ref_path_name,
                                        ref_offset_of(ref_offsets, ref_path_name),
                                        record_key_of(snarl), placement, false, 0,
                                        &site_panel);
            added = true;
        } else if (!gaf_output) {
            // Recorded here rather than in emit_variant. A retained chain, on the path above, is
            // recorded only if the linkage pass later finds that the sample carries it.
            site_panel_set = linker.add(snarl, travs, trav_genotype, trav_score,
                                        ref_trav_idx, ref_path_name,
                                        ref_offset_of(ref_offsets, ref_path_name),
                                        record_key_of(snarl), placement, false, 0,
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
            pending_this->parent_record_key = placement.parent_record_key;
            pending_this->parent_crossing = placement.parent_crossing;
            pending_this->chain_key = placement.chain_key;
            pending_this->no_reference = no_ref_position;
            pending_this->reported_inline = placement.reported_inline;
            pending_this->position_from_parent =
                no_ref_position
                    ? base_path_position(ref_path_name,
                                         get<0>(ref_interval) + ref_offset_of(ref_offsets, ref_path_name))
                          + (int64_t)placement.parent_offset
                    : 0;
            pending_this->chain_offset = placement.parent_offset;
            pending_this->crossing_known = placement.crossing_known;
            pending_this->level = (uint8_t)min(placement.level, (size_t)255);
            pending_this->set_call(std::move(trav_call_info), trav_score);
        }
        ret_val = trav_genotype.size() == ploidy && added;
    } else {
        // Top-level snarl or no parent context - genotype from scratch using support
        unique_ptr<SnarlCaller::CallInfo> trav_call_info;
        SiteScore* trav_score = nullptr;
        std::tie(trav_genotype, trav_call_info) = genotype_site(
            snarl, travs, ref_trav_idx, Ploidies{.ploidy = ploidy}, ref_path_name,
            make_pair(get<0>(ref_interval), get<1>(ref_interval)), trav_score);

        assert(trav_genotype.empty() || trav_genotype.size() == ploidy);

        bool added = true;
        added = stage_or_emit(trav_call_info, trav_score);

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
        // A child no called allele reaches is kept only where the linkage pass can come back to
        // it: with staging and the linkage model.
        const vector<ChildPlacer::Placed> children = child_placer.place(
            snarl, record_key_of(snarl), travs, trav_genotype, ref_trav_idx, ploidy, placement,
            off_reference_nesting, staged_sites.active() && linker.enabled());
        for (const ChildPlacer::Placed& child : children) {
            if (child.placement.level < 16) {
                ++descent_counters.depth_hist[child.placement.level];
            }
            // The other ploidy's answer is computed as well, so the linkage pass can change it
            // later.
            call_snarl_internal(*child.snarl, ref_path_name,
                                make_pair(get<0>(ref_interval), get<1>(ref_interval)),
                                nullptr, child.ploidy, child.placement);
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
                                        &child_trav_sets, -1, placement);
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

