#include <atomic>
#include <fstream>
#include <limits>

#include <omp.h>

#include "multipass_caller.hpp"
#include "utility.hpp"

namespace vg::multipass {

MultiPassCaller::MultiPassCaller(const PathPositionHandleGraph& graph,
                                 ReadLikelihoodSnarlCaller& genotyper,
                                 const SnarlDecomposition& decomposition,
                                 VCFOutputCaller& output, const string& sample_name,
                                 SiteTraversalFinder& traversal_finder,
                                 ReferencePathTable ref_paths, bool genotype_snarls,
                                 const pair<size_t, size_t>& allele_length_range, bool top_down,
                                 bool star_allele) :
    graph(graph),
    snarl_caller(genotyper),
    sites(decomposition, graph),
    walker(sites, graph),
    output(output),
    sample_name(sample_name),
    ref_paths(std::move(ref_paths)),
    genotype_snarls(genotype_snarls),
    top_down(top_down),
    star_allele(star_allele),
    splicer(graph, sites, phase_table),
    off_panel(graph, this->ref_paths),
    candidates(graph, traversal_finder, this->ref_paths, allele_length_range),
    record_renderer(phase_table, read_strands, linker, anchor_collector)
{
    install_header_steps();
    install_widgets();
}

void MultiPassCaller::call(GraphCaller::RecurseType recurse_type,
                           const function<void()>& after_direct_pass) {
    // Sized here rather than inside the parallel region that writes it.
    staged_sites.start(max((size_t)get_thread_count(), (size_t)omp_get_max_threads()));
    if (off_panel.is_enabled()) {
        off_panel.open(max((size_t)get_thread_count(), (size_t)omp_get_max_threads()));
    }
    // Configured here, once call_main has set nested calling and off-reference descent.
    tree_genotyper.configure(
        TreeGenotyper::Parts{
            .graph = &graph,
            .candidates = &candidates,
            .genotyper = &snarl_caller,
            .linker = &linker,
            .staged_sites = &staged_sites,
            .descent_policy = &descent_policy,
            .descent_counters = &descent_counters,
            .ploidy_regions = &ploidy_regions,
            .ref_paths = &ref_paths,
            .site_id = [this](const handle_t& start, const handle_t& end) {
                return output.print_snarl(&graph, start, end);
            },
            .edit_source = off_panel.source(),
            .read_source = read_source,
            .edit_outcomes = off_panel.outcomes(),
            .sites = &sites,
        },
        TreeGenotyper::Options{
            .nested_calling = nested_calling,
            .off_reference = off_reference_nesting,
            .top_down = top_down,
            .star_allele = star_allele,
            .edit_bar = off_panel.bar(),
        });

    // The direct pass.
    walker.walk(recurse_type, snarl_batch_window, show_progress,
                [&](const net_handle_t& site) { return tree_genotyper.genotype(site); });
    if (off_panel.is_enabled()) {
        // The edit sites no site holds, from the windows the sites fetched. A window no site
        // fetched holds no site, and is not fetched for its candidates alone.
        tree_genotyper.genotype_top_level_edits(off_panel.source()->counted_edit_candidates());
    }
    if (show_progress) {
        report_descent_instrumentation();
    }
    after_direct_pass();

    // Round 1's linkage pass, then read phasing and any further rounds, then the render.
    run_linkage_pass();
    phase_and_regenotype();
    render_retained_records();
    // Anchors are collected while records are built, so they are written afterwards. Does nothing
    // unless anchors are on.
    anchor_collector.write(sample_name);

    // The reports, and the mosaic, need to know which sites have a line, so they come last.
    finalise_linkage_outputs();
    off_panel.write(tree_genotyper.edit_nodes_after_pass());
    if (phase_declined.load() > 0 || quality_declined.load() > 0) {
        cerr << "[vg call] linkage: " << phase_declined.load()
             << " phases refused by the record they were rendered onto, and "
             << quality_declined.load() << " quality rewrites refused" << endl;
    }
    block_records.report();
}

unique_ptr<SiteTree> MultiPassCaller::site_tree() {
    return std::make_unique<StagedSiteTree>(
        staged_sites, sites, graph,
        [this](const handle_t& start, const handle_t& end) {
            return output.print_snarl(&graph, start, end);
        });
}

string MultiPassCaller::vcf_header(const PathHandleGraph& graph, const vector<string>& contigs,
                                   const vector<size_t>& contig_length_overrides) const {
    return output.genotyper_vcf_header(graph, contigs, contig_length_overrides,
                                       [this](string& header) { snarl_caller.update_vcf_header(header); });
}

void MultiPassCaller::report_descent_instrumentation() const {
    size_t total = 0;
    for (int d = 0; d < 16; ++d) {
        total += descent_counters.depth_hist[d].load();
    }
    if (total == 0) {
        // Nested descent did not run.
        return;
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

VCFOutputCaller::SiteRecordSteps MultiPassCaller::record_steps(
    const StagedSite& site, const AlleleSplicer::SplicedAlleles* spliced) {
    VCFOutputCaller::SiteRecordSteps steps;
    steps.id = &site.id;
    if (nested_calling) {
        steps.same_as_reference = [this, &site](int trav, int ref_trav_idx) {
            return is_symbolically_reference(site, trav, ref_trav_idx);
        };
        steps.count_site = [this, &site](const vector<int>&, int ref_trav_idx) {
            block_records.count_site(site.children, site.alleles.walks, ref_trav_idx);
        };
    }
    steps.phase = [this, &site](const vector<int>& site_genotype,
                                const map<int, int>& trav_to_allele, string& gt) {
        return phase_record_genotype(site.record_key, site_genotype, trav_to_allele, gt);
    };
    // The read-likelihood genotyper writes GL in colexicographic order, and the support-based one
    // in i-major order.
    steps.gl_layout = [](const SnarlCaller::CallInfo* call_info) {
        // Every call info this caller writes is the read-likelihood genotyper's.
        return call_info != nullptr ? GLLayout::Colexicographic : GLLayout::IMajor;
    };
    steps.write_blocks = [this, &site, spliced](const PathPositionHandleGraph& graph,
                                                const vector<int>& genotype, int ref_trav_idx,
                                                const SiteRecord& record, GLLayout gl_layout,
                                                bool genotype_snarls) {
        // The blocks spell the alleles the site record spells.
        const AlleleSequences site_steps =
            spliced != nullptr ? AlleleSequences() : site.alleles.all_step_sequences(graph);
        return block_records.write(graph, site.children,
                                   spliced != nullptr ? spliced->travs : site.alleles.walks,
                                   genotype,
                                   ref_trav_idx, sample_name, output.get_translation(), record,
                                   gl_layout, genotype_snarls,
                                   [this](vcflib::Variant& line, size_t block) {
                                       finish_moved_record(line);
                                       return output.add_variant(line, block);
                                   },
                                   spliced != nullptr ? &spliced->sequences : &site_steps);
    };
    steps.finish_record = [this](vcflib::Variant& record) {
        finish_moved_record(record);
    };
    // A written allele's scored allele is its traversal, since the site's walks are the ones it
    // was scored with, in order. Two alleles can share a walk when one has step sequences, so
    // they are told apart by index rather than by walk.
    steps.fill_info = [this, &site](const vector<int>& site_trav, const vector<int>& site_genotype,
                                    vcflib::Variant& variant) {
        snarl_caller.write_vcf_fields(site_trav, site_genotype, site.call_info.get(), sample_name,
                                      variant);
        if (sites.is_edit(site.net)) {
            variant.infoFlags["OP"] = true;
        }
    };
    // The linkage model gets the site whether or not it has a line. A parent written as the
    // reference still has two alleles, which differ only inside its children, and the children
    // need them to know which strand carries the chain. In VCF allele numbering such a parent is
    // 0/0; only in traversal space is it heterozygous.
    steps.site_filed = [this, &site](const map<int, int>& trav_to_allele,
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
        linker.collector()->set_allele_map(site.record_key, trav_to_allele_vec, has_line);
    };
    return steps;
}

void MultiPassCaller::finish_moved_record(vcflib::Variant& record) {
    if (!linker.enabled()) {
        return;
    }
    // The record already carries the chosen genotype, since it was built from it.
    const auto& quality = linker.collector()->moved_quality();
    if (quality.empty()) {
        return;
    }
    // Keyed by the site's record key: the hash of the record's ID, or of the site's ID for a block
    // record.
    auto found = quality.find(record_key_of(block_site_name(record.id)));
    if (found != quality.end()
        && !ReadLikelihoodSnarlCaller::rewrite_quality_for_chosen_genotype(
               record, sample_name, found->second, linkage_min_confidence)) {
        ++quality_declined;
    }
}

/// The header line of INFO/OP.
static const char* const OP_HEADER_LINE =
    "##INFO=<ID=OP,Number=0,Type=Flag,Description=\"Off-panel variant (--off-panel): "
    "the ALT is a change recurring in the reads that no walk of the graph spells, "
    "called at an edit site. Its ID starts with op\">";

/// The header line of INFO/SB.
static const char* const SB_HEADER_LINE =
    "##INFO=<ID=SB,Number=2,Type=Integer,Description=\"Index and count of this "
    "difference block within its snarl. A snarl is written as one record per "
    "difference block where the reference and the called haplotypes differ from each "
    "other in more than one place inside it, or where its own record would repeat a "
    "child snarl's, so the count can be 1. A block record's ID is the snarl's ID "
    "with _ and the index appended. DOUBLE COUNTING: the per-sample evidence is the "
    "SNARL's, repeated on every block, not apportioned between them -- AD, GL, GQ, "
    "GQI, GP and QUAL are identical across the set, because the genotype likelihood "
    "was computed over whole-snarl traversals and has no per-block decomposition. "
    "DP, DR and BL are per-site read counts and are site-level by definition. So any "
    "consumer that sums, averages or otherwise aggregates evidence across records "
    "must group by the snarl's ID first and count each snarl once. Records without "
    "SB are unaffected: they are the only record their snarl emitted.\">";

void MultiPassCaller::install_header_steps() {
    VCFOutputCaller::HeaderSteps steps;
    steps.format_header = [this]() {
        stringstream ss;
        if (emit_phasing) {
            // FORMAT/PS is the VCF phase set, which phasing tools read. It is unrelated to INFO/PS,
            // vg's parent-snarl field; the two are in different namespaces, so both are legal, and
            // their descriptions say which is which.
            ss << "##FORMAT=<ID=PS,Number=1,Type=Integer,Description=\"Phase set: the phase of a "
               << "genotype is comparable only with others carrying the same PS. One phase set per "
               << "chain, so blocks are chromosome-scale -- much longer than a read-based phaser "
               << "gives, because the phase comes from the haplotype panel rather than from reads "
               << "spanning consecutive sites. Not the INFO/PS emitted under -A, which is a parent "
               << "snarl pointer\">" << endl;
        }
        return ss.str();
    };
    steps.info_header = [this]() {
        stringstream ss;
        if (off_panel.is_enabled()) {
            ss << OP_HEADER_LINE << endl;
        }
        if (block_records.is_enabled()) {
            ss << SB_HEADER_LINE << endl;
        }
        return ss.str();
    };
    output.set_header_steps(std::move(steps));
}


void MultiPassCaller::install_widgets() {
    const SiteReader reader{
        .graph = &graph,
        .sites = &sites,
        .genotyper = &snarl_caller,
    };
    linker.set_site_reader(reader);
    rescorer.set_site_reader(reader);
    descent_policy.configure(&graph, &sites, &block_records, &descent_counters);
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

void MultiPassCaller::phase_reads() {
    // Read phasing changes the linkage model's phase, so it does not run without it.
    if (linker.enabled()) {
        read_phaser.phase(staged_sites, phase_table, read_strands,
                          [this](const string& contig, size_t phase_set) {
                              return phase_set_id(contig, phase_set);
                          });
    }
}

bool MultiPassCaller::rescore() {
    // Re-genotyping changes the likelihoods the linkage model chooses from, so it does not run
    // without it.
    return linker.enabled()
           && rescorer.rescore(staged_sites, phase_table, read_strands, temper_fit,
                               [this](const string& contig, size_t phase_set) {
                                   return phase_set_id(contig, phase_set);
                               },
                               show_progress);
}

bool MultiPassCaller::run_regenotyping_round(size_t round, RoundHistory& history) {
    const RoundHistory::State before = RoundHistory::state(staged_sites, linker.collector());
    if (round == 2) {
        // Remember round 1's genotypes.
        history.remember(before);
    }
    const bool calls_moved = rescore();
    // Choose the genotypes even when the correction moved no direct call. It has already changed
    // every site's likelihoods in place, and GL is written from them.
    rerun_linkage_pass();
    phase_reads();
    const RoundHistory::State after = RoundHistory::state(staged_sites, linker.collector());
    const size_t moved = RoundHistory::changed(before, after);
    cerr << "[vg call] re-genotyping round " << round << ": " << moved
         << " chosen genotypes moved" << endl;
    if (!calls_moved) {
        cerr << "[vg call] re-genotyping: the correction moved no site's direct call;"
             << " stopping after round " << round << endl;
        return false;
    }
    if (moved == 0) {
        cerr << "[vg call] re-genotyping: converged after " << round << " rounds" << endl;
        return false;
    }
    if (const size_t entered = history.first_round_with(after); entered != 0) {
        cerr << "[vg call] re-genotyping: LIMIT CYCLE of period "
             << (history.rounds() - entered + 1) << ", entered at round " << entered
             << ". The iteration does not converge and no round of a cycle is more"
             << " the answer than another; stopping here and reporting it rather than"
             << " presenting round " << round << " as a fixed point" << endl;
        return false;
    }
    history.remember(after);
    if (round == rescorer.passes()) {
        if (rescorer.passes() == 2) {
            // This is the default number of passes.
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
    return true;
}

void MultiPassCaller::phase_and_regenotype() {
    // Round 1 ends with read phasing. Its linkage pass has already run.
    phase_reads();
    if (rescorer.enabled() && rescorer.passes() >= 2) {
        // Remember every state the rounds reach, so that a cycle is recognised.
        RoundHistory history;
        for (size_t round = 2; round <= rescorer.passes(); ++round) {
            if (!run_regenotyping_round(round, history)) {
                break;
            }
        }
    } else {
        // With one pass, compute and report the correction, and keep nothing.
        rescore();
    }
}

void MultiPassCaller::render_retained_records() {
    // Each read's strand log-odds, for the anchors collected during the render and the hand-off.
    // Read phasing is done, so `read_strands` is final here.
    read_strands.build_lambda(
        phase_table.calls(),
        [&](const string& contig, size_t phase_set) { return phase_set_id(contig, phase_set); },
        rescorer.last_counters().calibration,
        rescorer.params());
    // The phase, before any record is built, so that each record is phased as it is rendered. Also
    // before the hand-off, which collects anchors for the records that get no line
    // (`reported_inline` and `no_reference`) and reads the frozen phase to order them. If read
    // phasing ran, the phase table already carries its swaps.
    phase_table.freeze_for_render(emit_phasing);
    record_renderer.render(
        staged_sites,
        [&](const StagedSite& site, const vector<int>& genotype) {
            // A parent spells its children's chosen alleles.
            AlleleSplicer::SplicedAlleles spliced;
            const bool is_spliced = splicer.splice(site, spliced);
            const vector<Traversal>& spelled = is_spliced ? spliced.travs : site.alleles.walks;
            const AlleleSequences site_steps =
                is_spliced ? AlleleSequences() : site.alleles.all_step_sequences(graph);
            const AlleleSequences& spelled_steps = is_spliced ? spliced.sequences : site_steps;
            output.emit_site(graph, sites.get_site_start(site.net), sites.get_site_end(site.net),
                             site.alleles.walks, genotype, site.alleles.ref_allele,
                             site.call_info.get(), site.reference.path_name,
                             site.reference.offset, genotype_snarls, site.ploidy,
                             record_steps(site, is_spliced ? &spliced : nullptr),
                             // Each allele is spelled as the site holds it: by its step sequences
                             // where it has them, or else by its walk.
                             [&](int trav, int) {
                                 return allele_sequence(graph, spelled[trav],
                                                        own_sequences(&spelled_steps, trav));
                             });
        },
        show_progress, [&]() { splicer.index(staged_sites, emit_phasing); });
    splicer.clear();
    staged_sites.release_unrendered();
}

size_t MultiPassCaller::rederive_reported_inline() {
    if (!linker.enabled() || !staged_sites.has_children()) {
        return 0;
    }
    vector<StagedSite>& pending = staged_sites.nested();
    size_t rederived = 0;
    const unordered_map<size_t, StagedSite*> record_by_key = staged_sites.by_key();
    // Count the chains held back from here, so that the report gives the chains held back now.
    block_records.restart_inline_count();
    // Visit parents before their children, so that a chain inherits its parent's final flag.
    staged_sites.for_each_parent_top_down(record_by_key, [&](const StagedSite& parent,
                                                             const vector<size_t>& children) {
        if (parent.dropped) {
            // Its children were dropped with it.
            return;
        }
        const BlockRecordWriter::ChainInlineContext ctx = block_records.chain_inline_context(
            graph, parent.children, parent.alleles.walks, linker.chosen_genotype(parent),
            parent.alleles.ref_allele);
        for (size_t ci : children) {
            StagedSite& child = pending[ci];
            if (child.dropped) {
                continue;
            }
            const bool was = child.placement.reported_inline;
            child.placement.reported_inline =
                parent.placement.reported_inline
                || block_records.chain_reported_inline(graph, ctx,
                                                       chain_holding(sites, graph, child.net));
            if (was != child.placement.reported_inline) {
                ++rederived;
            }
        }
    });
    return rederived;
}

/// Count the bytes of a map of genotype likelihoods, with an estimate of 48 bytes per map node.
static size_t gl_bytes(const map<vector<int>, double>& gl) {
    size_t n = 0;
    for (const auto& kv : gl) {
        n += 48 + kv.first.capacity() * sizeof(int) + sizeof(double);
    }
    return n;
}

/// Add the memory that the read-likelihood score of a staged site keeps to the totals.
static void measure_score(const SiteScore& rl, size_t& bytes, size_t& gls) {
    for (const auto& kv : rl.genotype_lls) {
        ++gls;
        bytes += 48 + kv.first.capacity() * sizeof(int) + sizeof(double);
    }
    if (rl.anchor_evidence != nullptr) {
        bytes += rl.anchor_evidence->bytes();
    }
    if (rl.phase_evidence != nullptr) {
        bytes += rl.phase_evidence->bytes();
    }
    if (rl.uncorrected_lls != nullptr) {
        bytes += gl_bytes(*rl.uncorrected_lls);
    }
    bytes += rl.scored_traversals.capacity() * sizeof(Traversal)
             + rl.allele_support.capacity() * sizeof(double);
    if (rl.alt_ploidy_info != nullptr) {
        const auto& alt = *rl.alt_ploidy_info;
        bytes += alt.scored_traversals.capacity() * sizeof(Traversal)
                 + alt.allele_support.capacity() * sizeof(double);
        if (alt.uncorrected_lls != nullptr) {
            bytes += gl_bytes(*alt.uncorrected_lls);
        }
        for (const auto& kv : alt.genotype_lls) {
            ++gls;
            bytes += 48 + kv.first.capacity() * sizeof(int) + sizeof(double);
        }
    }
}

/// Add the memory that a staged site keeps to the totals.
static void measure_staged(const StagedSite& rec, size_t& bytes, size_t& visits, size_t& gls) {
    bytes += sizeof(StagedSite) + rec.reference.path_name.capacity()
             + rec.genotype.capacity() * sizeof(int)
             + (rec.panel_cache ? rec.panel_cache->capacity() * sizeof(int) : 0);
    bytes += rec.alleles.walks.capacity() * sizeof(Traversal);
    for (const Traversal& t : rec.alleles.walks) {
        visits += t.size();
        bytes += t.capacity() * sizeof(handle_t);
    }
    bytes += rec.children.chains.capacity() * sizeof(ChildChain)
             + rec.children.entries.capacity() * sizeof(rec.children.entries[0]);
    if (rec.score != nullptr) {
        measure_score(*rec.score, bytes, gls);
    }
}

void MultiPassCaller::report_linkage_pass(const GenotypeLinker::PassCounts& counts,
                                          size_t inline_rederived) {
    const vector<StagedSite>& pending = staged_sites.nested();
    // Walk the objects in parallel. The totals are sums, so they do not depend on how the walk is
    // split.
    size_t retained_bytes = 0;
    size_t retained_visits = 0;
    size_t retained_gls = 0;
#pragma omp parallel reduction(+ : retained_bytes, retained_visits, retained_gls)
    {
        for (size_t q = 0; q < staged_sites.queue_count(); ++q) {
            const vector<StagedSite>& queue = staged_sites.queue(q);
#pragma omp for schedule(dynamic, 4096) nowait
            for (size_t r = 0; r < queue.size(); ++r) {
                measure_staged(queue[r], retained_bytes, retained_visits, retained_gls);
            }
        }
#pragma omp for schedule(dynamic, 4096) nowait
        for (size_t r = 0; r < pending.size(); ++r) {
            measure_staged(pending[r], retained_bytes, retained_visits, retained_gls);
        }
    }
    // Add the read-phasing evidence.
    for (const PhaseSite& ps : read_strands.sites()) {
        retained_bytes += sizeof(PhaseSite) + ps.read_key.capacity() * sizeof(uint64_t)
                          + ps.q0.capacity() * sizeof(float) + ps.c.capacity() * sizeof(float);
    }
    retained_bytes += read_strands.flips().size() * (sizeof(size_t) + 16);
    // The sites the linkage pass will not revise are the top-level ones, and the children that
    // RecurseOnFail reaches without a ploidy override.
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
    if (inline_rederived > 0) {
        cerr << "[vg call] linkage pass: " << inline_rederived
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

void MultiPassCaller::run_linkage_pass() {
    const GenotypeLinker::PassCounts counts = linker.link(staged_sites, phase_table, emit_phasing);
    const size_t inline_rederived = rederive_reported_inline();
    if (show_progress) {
        report_linkage_pass(counts, inline_rederived);
    }
}

void MultiPassCaller::set_linkage(LinkageCollector* collector, const gbwt::GBWT* gbwt,
                                  const vector<size_t>* sequence_to_haplotype) {
    this->panel_lookup = PanelLookup(&graph, gbwt, sequence_to_haplotype,
                                     collector != nullptr ? collector->panel_size() : 0);
    linker.configure(collector, &panel_lookup);
}

size_t MultiPassCaller::phase_set_id(const string& contig, size_t phase_set) {
    return phase_set_ids.emplace(make_pair(contig, phase_set), phase_set_ids.size()).first->second;
}

void MultiPassCaller::finalise_linkage_outputs() {
    // Built after every record has been rendered, since the mosaic needs to know which sites have
    // a line, which is not known while genotypes are being resolved.
    if (!linker.enabled()) {
        return;
    }
    // Read from the collector, since each PhaseCall's `emitted` was copied before any line was
    // written.
    const std::unordered_set<size_t> emitted_records = linker.collector()->emitted_records();
    size_t unexplained = 0;
    size_t order_arbitrary = 0;
    // Count the phased sites, separating those that became records from those that did not.
    size_t phased_unwritten = 0;
    for (const LinkageCollector::PhaseCall& pc : phase_table.calls()) {
        if (emitted_records.count(pc.record_key) == 0) {
            // Phased, since its children take their strand from it, but not a record, so it is kept
            // out of the mosaic and the record counts.
            ++phased_unwritten;
            continue;
        }
        // Count only the strands a site has. A haploid site has one strand and an empty slot, and
        // the empty slot can be either: a haploid contig fills the first slot, while a nested site
        // on its parent's second strand fills the second.
        const bool first_unexplained = pc.hap_first == LinkageModel::NO_HAPLOTYPE;
        const bool second_unexplained = pc.hap_second == LinkageModel::NO_HAPLOTYPE;
        if (pc.ploidy == 1 ? (first_unexplained && second_unexplained)
                           : (first_unexplained || second_unexplained)) {
            ++unexplained;
        }
        if (pc.order_arbitrary) {
            ++order_arbitrary;
        }
    }
    linker.report();
    if (emit_phasing) {
        // At sites where a strand names no panel haplotype, the phase across them rests on the
        // transitions alone.
        cerr << "[vg call] phasing: " << (phase_table.calls().size() - phased_unwritten)
             << " sites phased, " << unexplained
             << " with a strand the panel does not explain" << endl;
        if (phased_unwritten > 0) {
            // Sites that wrote no VCF line but are phased. A parent whose alleles differ only inside
            // its children is written as the reference and has no line, and its children still need
            // to know which of its strands carries the chain.
            cerr << "[vg call] phasing: " << phased_unwritten
                 << " collapsed sites phased with no line of their own, so their children can"
                 << " inherit a strand" << endl;
        }
        if (order_arbitrary > 0) {
            // Heterozygous sites where no panel haplotype on either strand carries either called
            // allele. The record is still phased and in the phase set, but its order came from
            // sorting the pair, so it is arbitrary.
            cerr << "[vg call] phasing: " << order_arbitrary
                 << " heterozygous sites carry an allele order the panel does not determine"
                 << endl;
        }
    }
    if (mosaic_writer.is_enabled()) {
        // Records only: the mosaic's segments are runs over sites of the call set, and it accounts
        // for exactly the written records.
        vector<LinkageCollector::PhaseCall> written;
        written.reserve(phase_table.calls().size());
        for (const LinkageCollector::PhaseCall& pc : phase_table.calls()) {
            if (emitted_records.count(pc.record_key) != 0) {
                written.push_back(pc);
            }
        }
        mosaic_writer.write(written, panel_lookup, sample_name);
    }
}

bool MultiPassCaller::is_symbolically_reference(const StagedSite& site, int trav_idx,
                                                int ref_trav_idx) const {
    // Only with nested calling.
    if (!nested_calling || ref_trav_idx < 0 || trav_idx < 0 ||
        ref_trav_idx >= (int)site.alleles.walks.size()
        || trav_idx >= (int)site.alleles.walks.size()) {
        return false;
    }
    // An allele that changes the sequences of its walk spells something the walk does not, so
    // sharing the reference's route does not make it the reference.
    if (!site.alleles.spelled_by_walk(trav_idx) || !site.alleles.spelled_by_walk(ref_trav_idx)) {
        return false;
    }
    return symbolically_equal(graph, site.alleles.walks[trav_idx], site.alleles.walks[ref_trav_idx],
                              site.children);
}

void MultiPassCaller::set_nested_calling(bool on) {
    nested_calling = on;
    block_records.set_nested(on);
}

int64_t MultiPassCaller::phase_record_genotype(size_t record_key,
                                               const vector<int>& site_genotype,
                                               const map<int, int>& trav_to_allele,
                                               string& gt) const {
    if (!emit_phasing || !phase_table.has_rendered()) {
        return -1;
    }
    const LinkageCollector::PhaseCall* found = phase_table.rendered(record_key);
    if (found == nullptr) {
        return -1;
    }
    const LinkageCollector::PhaseCall& phase = *found;
    if (phase.ploidy == 2 && phase.order_arbitrary
        && !read_strands.read_ordered().count(record_key)
        && !splicer.has_children(record_key)) {
        // Neither the panel nor the reads order the pair, so the record claims no phase. A site
        // with sites staged under it keeps its phase, since its one-copy children's `a|.` GTs name
        // its strands; so does a one-copy chain, whose GT shape says which strand it is on.
        return -1;
    }
    // `find`, since `operator[]` would insert a default 0 on a miss, and the map's size is not
    // a bound on traversal indices.
    const auto found_a = trav_to_allele.find(phase.trav_first);
    const auto found_b = trav_to_allele.find(phase.trav_second);
    const int a = (phase.trav_first >= 0 && found_a != trav_to_allele.end())
                      ? found_a->second : -1;
    const int b = (phase.trav_second >= 0 && found_b != trav_to_allele.end())
                      ? found_b->second : -1;
    // The phased genotype must be a permutation of the one this record carries, so that
    // phasing cannot change a genotype.
    bool same = false;
    if (phase.ploidy == 1 && site_genotype.size() == 1) {
        same = (a >= 0 && a == site_genotype[0]);
    } else if (phase.ploidy == 2 && site_genotype.size() == 2) {
        same = (a >= 0 && b >= 0)
               && ((a == site_genotype[0] && b == site_genotype[1])
                   || (a == site_genotype[1] && b == site_genotype[0]));
    }

    if (!same) {
        ++phase_declined;
        return -1;
    }
    if (phase.ploidy == 1 && phase.nested_strand >= 0) {
        // A nested ploidy-1 site is one strand of a diploid locus, since the parent's
        // other allele deletes the chain. Written as a phased pair with "." on the other
        // strand, which is how the VCF records which strand carries the allele.
        gt = nested_strand_genotype(a, phase.nested_strand);
    } else if (phase.ploidy == 1) {
        // A haploid locus: one allele and no order; PS labels its phase set. "a|a"
        // would claim a homozygous diploid call.
        gt = std::to_string(a);
    } else {
        gt = std::to_string(a) + "|" + std::to_string(b);
    }
    return (int64_t)phase.phase_set;
}

}
