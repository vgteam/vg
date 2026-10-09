#include <atomic>
#include <limits>

#include <omp.h>

#include "multipass_caller.hpp"
#include "utility.hpp"

namespace vg {

MultiPassCaller::MultiPassCaller(const PathPositionHandleGraph& graph,
                                 ReadLikelihoodSnarlCaller& genotyper, SnarlManager& snarl_manager,
                                 const string& sample_name, TraversalFinder& traversal_finder,
                                 const vector<string>& ref_paths,
                                 const vector<size_t>& ref_path_offsets,
                                 const vector<int>& ref_path_ploidies, bool genotype_snarls,
                                 const pair<size_t, size_t>& allele_length_range, bool top_down,
                                 bool star_allele) :
    VCFOutputCaller(sample_name),
    graph(graph),
    snarl_caller(genotyper),
    snarl_manager(snarl_manager),
    ref_paths(ref_paths),
    genotype_snarls(genotype_snarls),
    top_down(top_down),
    star_allele(star_allele),
    candidates(graph, snarl_manager, traversal_finder, genotyper.get_support_finder(),
               ref_path_set, allele_length_range),
    site_genotyper(genotyper)
{
    for (int i = 0; i < ref_paths.size(); ++i) {
        ref_offsets[ref_paths[i]] = i < ref_path_offsets.size() ? ref_path_offsets[i] : 0;
        ref_path_set.insert(ref_paths[i]);
        ref_ploidies[ref_paths[i]] = i < ref_path_ploidies.size() ? ref_path_ploidies[i] : 2;
    }
    install_record_steps();
    install_widgets();
}

void MultiPassCaller::call(GraphCaller::RecurseType recurse_type,
                           const function<void()>& after_direct_pass) {
    // Sized here rather than inside the parallel region that writes it.
    staged_sites.start(max((size_t)get_thread_count(), (size_t)omp_get_max_threads()));
    // Configured here, once call_main has set nested calling and off-reference descent.
    tree_genotyper.configure(
        TreeGenotyper::Parts{
            .graph = &graph,
            .snarl_manager = &snarl_manager,
            .candidates = &candidates,
            .genotyper = &site_genotyper,
            .linker = &linker,
            .staged_sites = &staged_sites,
            .child_placer = &child_placer,
            .descent_counters = &descent_counters,
            .ploidy_regions = &ploidy_regions,
            .ref_offsets = &ref_offsets,
            .ref_ploidies = &ref_ploidies,
            .record_key_of = [this](const Snarl& site) { return record_key_of(site); },
        },
        TreeGenotyper::Options{
            .nested_calling = symbolic_manager != nullptr,
            .off_reference = off_reference_nesting,
            .top_down = top_down,
            .star_allele = star_allele,
        });

    // The direct pass.
    call_snarl_tree(graph, snarl_manager, recurse_type, snarl_batch_window, show_progress,
                    [&](const Snarl& site) { return tree_genotyper.genotype(site); });
    if (show_progress) {
        report_descent_instrumentation();
    }
    after_direct_pass();

    // Round 1's linkage pass, then read phasing and any further rounds, then the render.
    run_linkage_pass();
    phase_and_regenotype();
    render_retained_records();
    // Anchors are collected while records are built, so they are written afterwards.
    write_anchors();
}

string MultiPassCaller::vcf_header(const PathHandleGraph& graph, const vector<string>& contigs,
                                   const vector<size_t>& contig_length_overrides) const {
    return snarl_caller_vcf_header(graph, contigs, contig_length_overrides, snarl_caller);
}

void MultiPassCaller::report_descent_instrumentation() const {
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

void MultiPassCaller::install_record_steps() {
    record_steps.phase = [this](const Snarl& site, const vector<int>& site_genotype,
                                const map<int, int>& trav_to_allele, string& gt) {
        return phase_record_genotype(site, site_genotype, trav_to_allele, gt);
    };
    // The read-likelihood genotyper writes GL in colexicographic order, and the support-based one
    // in i-major order.
    record_steps.gl_layout = [this](const SnarlCaller::CallInfo* call_info) {
        // Every call info this caller writes is the read-likelihood genotyper's.
        return call_info != nullptr ? GLLayout::Colexicographic : GLLayout::IMajor;
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

void MultiPassCaller::install_widgets() {
    const SiteReader reader{
        .graph = &graph,
        .genotyper = &site_genotyper,
        .spell = [this](const SnarlTraversal& trav) { return trav_string(graph, trav); },
        .name = [this](const Snarl& site) { return print_snarl(site); },
    };
    linker.set_site_reader(reader);
    rescorer.set_site_reader(reader);
    record_renderer.configure(reader, [this](const Snarl& site) { return snarl_is_leaf(site); });
    child_placer.configure(&graph, &snarl_manager, &block_records, &descent_counters);
}

bool MultiPassCaller::snarl_is_leaf(const Snarl& snarl) const {
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

void MultiPassCaller::rerun_linkage_pass() {
    if (!linker.enabled()) {
        return;
    }
    // Give the linkage model the corrected likelihoods, then run the linkage pass again in full, so that
    // every child is reassessed against its parent's new chosen pair, as on the first pass.
    linker.resync(staged_sites);
    run_linkage_pass();
}

void MultiPassCaller::phase_and_regenotype() {
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

void MultiPassCaller::render_retained_records() {
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

void MultiPassCaller::run_linkage_pass() {
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

}
