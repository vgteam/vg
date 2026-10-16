#include <atomic>
#include <fstream>
#include <limits>

#include <omp.h>

#include "multipass_caller.hpp"
#include "utility.hpp"

namespace vg {
namespace multipass {

MultiPassCaller::MultiPassCaller(const PathPositionHandleGraph& graph,
                                 ReadLikelihoodSnarlCaller& genotyper,
                                 const SnarlDecomposition& decomposition,
                                 VCFOutputCaller& output, const string& sample_name,
                                 TraversalFinder& traversal_finder,
                                 const vector<string>& ref_paths,
                                 const vector<size_t>& ref_path_offsets,
                                 const vector<int>& ref_path_ploidies, bool genotype_snarls,
                                 const pair<size_t, size_t>& allele_length_range, bool top_down,
                                 bool star_allele) :
    graph(graph),
    snarl_caller(genotyper),
    decomposition(decomposition),
    walker(decomposition, graph),
    output(output),
    sample_name(sample_name),
    ref_paths(ref_paths),
    genotype_snarls(genotype_snarls),
    top_down(top_down),
    star_allele(star_allele),
    candidates(graph, traversal_finder, ref_path_set, allele_length_range),
    site_genotyper(genotyper)
{
    for (int i = 0; i < ref_paths.size(); ++i) {
        ref_offsets[ref_paths[i]] = i < ref_path_offsets.size() ? ref_path_offsets[i] : 0;
        ref_path_set.insert(ref_paths[i]);
        ref_ploidies[ref_paths[i]] = i < ref_path_ploidies.size() ? ref_path_ploidies[i] : 2;
    }
    install_header_steps();
    install_widgets();
}

void MultiPassCaller::call(GraphCaller::RecurseType recurse_type,
                           const function<void()>& after_direct_pass) {
    // Sized here rather than inside the parallel region that writes it.
    staged_sites.start(max((size_t)get_thread_count(), (size_t)omp_get_max_threads()));
    if (edit_source != nullptr) {
        edit_outcomes.reset(
            new EditOutcomes(max((size_t)get_thread_count(), (size_t)omp_get_max_threads())));
    }
    // Configured here, once call_main has set nested calling and off-reference descent.
    tree_genotyper.configure(
        TreeGenotyper::Parts{
            .graph = &graph,
            .candidates = &candidates,
            .genotyper = &site_genotyper,
            .linker = &linker,
            .staged_sites = &staged_sites,
            .child_placer = &child_placer,
            .descent_counters = &descent_counters,
            .ploidy_regions = &ploidy_regions,
            .ref_offsets = &ref_offsets,
            .ref_ploidies = &ref_ploidies,
            .site_id = [this](const SiteBounds& site) {
                return output.print_snarl(&graph, site.start, site.end);
            },
            .edit_source = edit_source,
            .read_source = read_source,
            .edit_outcomes = edit_outcomes.get(),
        },
        TreeGenotyper::Options{
            .nested_calling = nested_calling,
            .off_reference = off_reference_nesting,
            .top_down = top_down,
            .star_allele = star_allele,
            .edit_bar = edit_bar,
        });

    // The direct pass.
    walker.walk(recurse_type, snarl_batch_window, show_progress,
                [&](const SiteView& site) { return tree_genotyper.genotype(site); });
    if (edit_source != nullptr) {
        // The edit sites no site holds, from the windows the sites fetched. A window no site
        // fetched holds no site, and is not fetched for its candidates alone.
        tree_genotyper.genotype_top_level_edits(edit_source->counted_edit_candidates());
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
    write_edit_dump();
    if (phase_declined.load() > 0 || quality_declined.load() > 0) {
        cerr << "[vg call] linkage: " << phase_declined.load()
             << " phases refused by the record they were rendered onto, and "
             << quality_declined.load() << " quality rewrites refused" << endl;
    }
    block_records.report();
}

pair<string, int64_t> MultiPassCaller::reference_position(nid_t node, uint32_t offset) const {
    pair<string, int64_t> found;
    const handle_t handle = graph.get_handle(node);
    graph.for_each_step_on_handle(handle, [&](const step_handle_t& step) {
        const string path_name = graph.get_path_name(graph.get_path_handle_of_step(step));
        if (!ref_path_set.count(path_name)) {
            return true;
        }
        const size_t len = graph.get_length(handle);
        const bool reverse = graph.get_is_reverse(graph.get_handle_of_step(step));
        const int64_t along = (int64_t)graph.get_position_of_step(step)
                              + (int64_t)(reverse ? len - 1 - offset : offset)
                              + ref_offset_of(ref_offsets, path_name);
        string contig = PathMetadata::parse_locus_name(Paths::strip_subrange(path_name));
        found = make_pair(contig == PathMetadata::NO_LOCUS_NAME ? Paths::strip_subrange(path_name)
                                                                : contig,
                          base_path_position(path_name, along));
        return false;
    });
    return found;
}

void MultiPassCaller::write_edit_dump() const {
    if (edit_source == nullptr) {
        return;
    }
    {
        // The feature's one block of counters, written only when it is on.
        size_t top = 0, nested = 0, staged = 0, unplaced = 0;
        for (const EditOutcome& o : edit_outcomes->sorted()) {
            top += o.placed == "top";
            nested += o.placed == "nested";
            unplaced += o.placed != "top" && o.placed != "nested";
            staged += o.staged;
        }
        cerr << "[vg call] off-panel: " << edit_source->counted_edit_candidates().size()
             << " candidates counted; the call bar's genotyped as " << top << " top-level and "
             << nested << " nested sites, " << unplaced << " left unplaced; " << staged
             << " staged with their ALT; " << tree_genotyper.edit_nodes_after_pass()
             << " nodes' top-level sites left to after the direct pass" << endl;
        size_t indels = 0, indels_staged = 0;
        for (const EditCandidate& c : edit_source->counted_edit_candidates()) {
            indels += c.is_indel();
        }
        for (const EditOutcome& o : edit_outcomes->sorted()) {
            indels_staged += o.staged && o.kind != 'S';
        }
        cerr << "[vg call] off-panel: " << indels << " of the candidates are indels; "
             << indels_staged << " indels staged with their ALT" << endl;
    }
    if (edit_dump.empty()) {
        return;
    }
    ofstream out(edit_dump);
    if (!out) {
        cerr << "error [vg call]: could not open " << edit_dump << " for the off-panel dump" << endl;
        return;
    }
    out << "#contig\tpos\tnode\toffset\tref\talt\talt_frags\tref_frags\tother_frags"
           "\tdiscordant_frags\taf\tmean_alt_mapq\tplaced\tsite\tgt\tgq\tstaged\thomopolymer\tdr\n";
    const vector<EditOutcome> outcomes =
        edit_outcomes != nullptr ? edit_outcomes->sorted() : vector<EditOutcome>();
    auto outcome = outcomes.begin();
    for (const EditCandidate& c : edit_source->counted_edit_candidates()) {
        const pair<string, int64_t> at = reference_position(c.node, c.offset);
        out << (at.first.empty() ? "." : at.first) << "\t";
        if (at.first.empty()) {
            out << ".";
        } else {
            out << at.second;
        }
        out << "\t" << c.node << "\t" << c.offset << "\t";
        if (c.is_indel()) {
            // As VCF writes an indel: from its anchor base.
            out << c.ref << (c.kind == 'D' ? c.sequence : "") << "\t" << c.ref
                << (c.kind == 'I' ? c.sequence : "");
        } else {
            out << c.ref << "\t" << c.alt;
        }
        out << "\t" << c.alt_fragments << "\t" << c.ref_fragments << "\t" << c.other_fragments << "\t"
            << c.discordant_fragments << "\t" << c.fraction() << "\t" << c.mean_alt_mapq;
        // The candidate's outcome, if the call bar admitted it.
        while (outcome != outcomes.end()
               && tie(outcome->node, outcome->offset, outcome->kind, outcome->alt,
                      outcome->sequence)
                      < tie(c.node, c.offset, c.kind, c.alt, c.sequence)) {
            ++outcome;
        }
        if (outcome != outcomes.end() && outcome->node == c.node && outcome->offset == c.offset
            && outcome->kind == c.kind && outcome->alt == c.alt
            && outcome->sequence == c.sequence) {
            out << "\t" << outcome->placed << "\t" << (outcome->site.empty() ? "." : outcome->site)
                << "\t";
            for (size_t i = 0; i < outcome->genotype.size(); ++i) {
                out << (i ? "/" : "") << outcome->genotype[i];
            }
            if (outcome->genotype.empty()) {
                out << ".";
            }
            out << "\t" << outcome->gq << "\t" << (outcome->staged ? 1 : 0);
        } else {
            out << "\t" << (edit_bar.admits(c) ? "unplaced" : "below_bar") << "\t.\t.\t.\t0";
        }
        out << "\t" << c.homopolymer << "\t";
        if (outcome != outcomes.end() && outcome->node == c.node && outcome->offset == c.offset
            && outcome->kind == c.kind && outcome->alt == c.alt
            && outcome->sequence == c.sequence && !outcome->genotype.empty()) {
            out << outcome->depth_ratio;
        } else {
            out << ".";
        }
        out << "\n";
    }
}

unique_ptr<SiteTree> MultiPassCaller::site_tree() {
    return unique_ptr<SiteTree>(new StagedSiteTree(
        staged_sites, graph,
        [this](const SiteBounds& site) { return output.print_snarl(&graph, site.start, site.end); }));
}

string MultiPassCaller::vcf_header(const PathHandleGraph& graph, const vector<string>& contigs,
                                   const vector<size_t>& contig_length_overrides) const {
    return output.snarl_caller_vcf_header(graph, contigs, contig_length_overrides, snarl_caller);
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

VCFOutputCaller::SiteRecordSteps MultiPassCaller::record_steps(const StagedSite& site,
                                                              const SplicedAlleles* spliced) {
    // Each step reads the staged site rather than the Snarl and SnarlTraversals emit_variant
    // passes, which are the site's bounds and walks in another form.
    VCFOutputCaller::SiteRecordSteps steps;
    steps.id = &site.id;
    if (nested_calling) {
        steps.same_as_reference = [this, &site](const Snarl&, const vector<SnarlTraversal>&,
                                                int trav, int ref_trav_idx) {
            return is_symbolically_reference(site, trav, ref_trav_idx);
        };
        steps.count_site = [this, &site](const PathPositionHandleGraph&, const Snarl&,
                                         const vector<SnarlTraversal>&, const vector<int>&,
                                         int ref_trav_idx) {
            block_records.count_site(site.children, site.travs, ref_trav_idx);
        };
    }
    steps.phase = [this, &site](const Snarl&, const vector<int>& site_genotype,
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
                                                const Snarl&, const vector<SnarlTraversal>&,
                                                const vector<int>& genotype, int ref_trav_idx,
                                                const SiteRecord& record, GLLayout gl_layout,
                                                bool genotype_snarls) {
        // The blocks spell the alleles the site record spells.
        return block_records.write(graph, site.children,
                                   spliced != nullptr ? spliced->travs : site.travs, genotype,
                                   ref_trav_idx, sample_name, output.get_translation(), record,
                                   gl_layout, genotype_snarls,
                                   [this](vcflib::Variant& line, size_t block) {
                                       finish_moved_record(line);
                                       return output.add_variant(line, block);
                                   },
                                   spliced != nullptr ? &spliced->sequences : &site.sequences);
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
        if (site.off_panel) {
            variant.infoFlags["OP"] = true;
        }
    };
    // The linkage model gets the site whether or not it has a line. A parent written as the
    // reference still has two alleles, which differ only inside its children, and the children
    // need them to know which strand carries the chain. In VCF allele numbering such a parent is
    // 0/0; only in traversal space is it heterozygous.
    steps.site_filed = [this, &site](const Snarl&, const map<int, int>& trav_to_allele,
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
        if (edit_source != nullptr) {
            ss << "##INFO=<ID=OP,Number=0,Type=Flag,Description=\"Off-panel variant (--off-panel): "
               << "the ALT is a change recurring in the reads that no walk of the graph spells, "
               << "called at an edit site. Its ID starts with op\">" << endl;
        }
        if (block_records.is_enabled()) {
            ss << "##INFO=<ID=SB,Number=2,Type=Integer,Description=\"Index and count of this "
               << "difference block within its snarl. A snarl is written as one record per "
               << "difference block where the reference and the called haplotypes differ from each "
               << "other in more than one place inside it, or where its own record would repeat a "
               << "child snarl's, so the count can be 1. A block record's ID is the snarl's ID "
               << "with _ and the index appended. DOUBLE COUNTING: the per-sample evidence is the "
               << "SNARL's, repeated on every block, not apportioned between them -- AD, GL, GQ, "
               << "GQI, GP and QUAL are identical across the set, because the genotype likelihood "
               << "was computed over whole-snarl traversals and has no per-block decomposition. "
               << "DP, DR and BL are per-site read counts and are site-level by definition. So any "
               << "consumer that sums, averages or otherwise aggregates evidence across records "
               << "must group by the snarl's ID first and count each snarl once. Records without "
               << "SB are unaffected: they are the only record their snarl emitted.\">" << endl;
        }
        return ss.str();
    };
    output.set_header_steps(std::move(steps));
}

void MultiPassCaller::install_widgets() {
    const SiteReader reader{
        .graph = &graph,
        .genotyper = &site_genotyper,
        .spell = [this](const Traversal& walk, const StepSequences* own) {
            return allele_sequence(graph, walk, own);
        },
    };
    linker.set_site_reader(reader);
    rescorer.set_site_reader(reader);
    record_renderer.configure(reader);
    child_placer.configure(&graph, &decomposition, &block_records, &descent_counters);
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
    // reach `rescorer.passes()` rounds. With one pass the correction is only computed and reported.
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
            // A parent spells its children's chosen alleles.
            SplicedAlleles spliced;
            const bool is_spliced = splice(site, spliced);
            const vector<Traversal>& spelled = is_spliced ? spliced.travs : site.travs;
            const AlleleSequences& spelled_steps =
                is_spliced ? spliced.sequences : site.sequences;
            // The VCF writer and the genotyper's record fields still take a Snarl and
            // SnarlTraversals.
            output.emit_variant(graph, snarl_caller, snarl_of(graph, site.bounds),
                                snarl_traversals_of(graph, site.travs), genotype,
                                site.ref_trav_idx, site.call_info, site.ref_path_name,
                                site.ref_offset, genotype_snarls, site.ploidy,
                                record_steps(site, is_spliced ? &spliced : nullptr),
                                // Each allele is spelled as the site holds it: by its step
                                // sequences where it has them, or else by its walk.
                                [&](const vector<SnarlTraversal>&, const vector<int>&, int trav,
                                    int, int) {
                                    return allele_sequence(graph, spelled[trav],
                                                           own_sequences(&spelled_steps, trav));
                                });
        },
        show_progress, [&]() { index_splice_children(); });
    splice_children.clear();
    staged_sites.release_unrendered();
}

void MultiPassCaller::index_splice_children() {
    splice_children.clear();
    if (!emit_phasing || !phase_table.has_rendered()) {
        return;
    }
    auto file = [&](const StagedSite& site) {
        if (site.parent_record_key != 0 && !site.dropped) {
            splice_children[site.parent_record_key].push_back(&site);
        }
    };
    for (size_t q = 0; q < staged_sites.queue_count(); ++q) {
        for (const StagedSite& site : staged_sites.queue(q)) {
            file(site);
        }
    }
    for (const StagedSite& site : staged_sites.unrendered()) {
        file(site);
    }
}

bool MultiPassCaller::splice(const StagedSite& site, SplicedAlleles& out) const {
    if (splice_children.empty() || !splice_children.count(site.record_key)) {
        return false;
    }
    const LinkageCollector::PhaseCall* phase = phase_table.rendered(site.record_key);
    if (phase == nullptr) {
        return false;
    }
    bool changed = false;
    const int strands = phase->ploidy == 1 ? 1 : 2;
    int first_trav = -1;
    for (int strand = 0; strand < strands; ++strand) {
        const int trav = strand == 0 ? phase->trav_first : phase->trav_second;
        if (trav < 0 || trav >= (int)site.travs.size() || trav == site.ref_trav_idx
            || trav == first_trav) {
            // The reference is spelled as the reference, and an allele as on its first strand.
            continue;
        }
        first_trav = trav;
        Traversal walk;
        StepSequences steps;
        if (spliced_allele(site, trav, strand, strands == 1, walk, steps)) {
            if (!changed) {
                out.travs = site.travs;
                out.sequences.assign(site.travs.size(), StepSequences());
                for (size_t t = 0; t < site.travs.size(); ++t) {
                    if (const StepSequences* own = own_sequences(&site.sequences, t)) {
                        out.sequences[t] = *own;
                    }
                }
                changed = true;
            }
            out.travs[trav] = std::move(walk);
            out.sequences[trav] = std::move(steps);
        }
    }
    return changed;
}

bool MultiPassCaller::spliced_allele(const StagedSite& site, int trav, int strand, bool haploid,
                                     Traversal& walk, StepSequences& steps) const {
    auto kids = splice_children.find(site.record_key);
    if (kids == splice_children.end()) {
        return false;
    }
    const Traversal& parent_walk = site.travs[trav];
    // Where each node is first visited, which is where a child is entered.
    unordered_map<nid_t, int> first_visit;
    for (int i = (int)parent_walk.size() - 1; i >= 0; --i) {
        first_visit[graph.get_id(parent_walk[i])] = i;
    }
    // Each child's chosen allele on this strand, as the steps that replace the parent's visits
    // [begin, end).
    struct Replacement {
        int begin;
        int end;
        Traversal walk;
        StepSequences steps;
        bool inside_node;
        size_t key;
    };
    vector<Replacement> replacements;
    for (const StagedSite* child : kids->second) {
        const LinkageCollector::PhaseCall* phase = phase_table.rendered(child->record_key);
        if (phase == nullptr) {
            continue;
        }
        int child_trav;
        if (phase->ploidy == 1) {
            // One copy, on the parent's strand that carries it.
            if (!haploid && phase->nested_strand != strand) {
                continue;
            }
            child_trav = phase->trav_first;
        } else {
            child_trav = strand == 0 ? phase->trav_first : phase->trav_second;
        }
        if (child_trav < 0 || child_trav >= (int)child->travs.size()) {
            continue;
        }
        Replacement rep;
        rep.inside_node = child->bounds.inside_node;
        rep.key = child->record_key;
        if (!spliced_allele(*child, child_trav, strand, phase->ploidy == 1, rep.walk, rep.steps)) {
            rep.walk = child->travs[child_trav];
            if (const StepSequences* own = own_sequences(&child->sequences, child_trav)) {
                rep.steps = *own;
            } else {
                rep.steps.clear();
                for (const handle_t& h : rep.walk) {
                    rep.steps.push_back(graph.get_sequence(h));
                }
            }
        }
        if (rep.walk.empty()) {
            continue;
        }
        // Where the parent's walk crosses the child, and which way.
        auto at_start = first_visit.find(graph.get_id(rep.walk.front()));
        auto at_end = first_visit.find(graph.get_id(rep.walk.back()));
        if (at_start == first_visit.end() || at_end == first_visit.end()) {
            continue;
        }
        bool turned = false;
        if (rep.inside_node) {
            // Placed on the spliced walk below; here only which way the parent reads it.
            rep.begin = at_start->second;
            rep.end = rep.begin + 1;
            turned = graph.get_is_reverse(parent_walk[rep.begin])
                     != graph.get_is_reverse(rep.walk.front());
        } else if (at_start->second <= at_end->second
                   && parent_walk[at_start->second] == rep.walk.front()) {
            rep.begin = at_start->second;
            rep.end = at_end->second + 1;
        } else if (at_end->second <= at_start->second
                   && parent_walk[at_end->second] == graph.flip(rep.walk.back())) {
            rep.begin = at_end->second;
            rep.end = at_start->second + 1;
            turned = true;
        } else {
            continue;
        }
        if (turned) {
            // Read the child's allele the way the parent's walk reads it.
            std::reverse(rep.walk.begin(), rep.walk.end());
            std::reverse(rep.steps.begin(), rep.steps.end());
            for (size_t i = 0; i < rep.walk.size(); ++i) {
                rep.walk[i] = graph.flip(rep.walk[i]);
                rep.steps[i] = reverse_complement(rep.steps[i]);
            }
        }
        replacements.push_back(std::move(rep));
    }
    if (replacements.empty()) {
        return false;
    }
    // Routes first, in walk order, then the children inside a node on the spliced walk, so that
    // the result does not depend on the order the children were staged in.
    vector<Replacement> routes;
    vector<Replacement> edits;
    for (Replacement& rep : replacements) {
        (rep.inside_node ? edits : routes).push_back(std::move(rep));
    }
    std::sort(routes.begin(), routes.end(), [](const Replacement& a, const Replacement& b) {
        return std::tie(a.begin, a.end, a.key) < std::tie(b.begin, b.end, b.key);
    });
    std::sort(edits.begin(), edits.end(),
              [](const Replacement& a, const Replacement& b) { return a.key < b.key; });

    const StepSequences* own = own_sequences(&site.sequences, trav);
    walk.clear();
    steps.clear();
    bool changed = false;
    size_t r = 0;
    for (int i = 0; i < (int)parent_walk.size();) {
        if (r < routes.size() && routes[r].begin < i) {
            ++r;   // overlaps one already spliced
            continue;
        }
        if (r < routes.size() && routes[r].begin == i) {
            const Replacement& rep = routes[r];
            // The child's allele replaces the parent's route through it.
            bool same = rep.end - rep.begin == (int)rep.walk.size();
            for (int k = 0; same && k < (int)rep.walk.size(); ++k) {
                same = parent_walk[rep.begin + k] == rep.walk[k]
                       && rep.steps[k] == (own != nullptr ? own->at(rep.begin + k)
                                                          : graph.get_sequence(rep.walk[k]));
            }
            changed = changed || !same;
            walk.insert(walk.end(), rep.walk.begin(), rep.walk.end());
            steps.insert(steps.end(), rep.steps.begin(), rep.steps.end());
            i = rep.end;
            ++r;
            continue;
        }
        walk.push_back(parent_walk[i]);
        steps.push_back(own != nullptr ? own->at(i) : graph.get_sequence(parent_walk[i]));
        ++i;
    }
    // A child inside nodes changes the bases its allele differs from its nodes in, where the
    // spliced walk first visits its first node and then the rest of its walk, read the way the
    // spliced walk reads them. Each step's changes are made right to left, so that one that changes
    // the length does not move the offsets of those before it.
    struct StepEdit {
        size_t step;
        size_t at;
        size_t removed;
        string added;
        size_t key;
    };
    vector<StepEdit> step_edits;
    for (const Replacement& edit : edits) {
        // Its walk, and the same read the other way, as the spliced walk may read it.
        Traversal turned_walk;
        StepSequences turned_steps;
        for (size_t k = edit.walk.size(); k-- > 0;) {
            turned_walk.push_back(graph.flip(edit.walk[k]));
            turned_steps.push_back(reverse_complement(edit.steps[k]));
        }
        const Traversal* child_walk = nullptr;
        const StepSequences* child_steps = nullptr;
        size_t i = 0;
        for (int way = 0; way < 2 && child_walk == nullptr; ++way) {
            const Traversal& w = way == 0 ? edit.walk : turned_walk;
            const nid_t first = graph.get_id(w.front());
            i = 0;
            while (i < walk.size() && graph.get_id(walk[i]) != first) {
                ++i;
            }
            if (i + w.size() <= walk.size() && std::equal(w.begin(), w.end(), walk.begin() + i)) {
                child_walk = &w;
                child_steps = way == 0 ? &edit.steps : &turned_steps;
            }
        }
        if (child_walk == nullptr) {
            continue;
        }
        for (size_t k = 0; k < child_walk->size(); ++k) {
            const string node_seq = graph.get_sequence(walk[i + k]);
            const string& child_seq = (*child_steps)[k];
            if (child_seq == node_seq) {
                continue;
            }
            // The child's change on this node: its longest shared prefix and suffix aside.
            size_t prefix = 0;
            while (prefix < child_seq.size() && prefix < node_seq.size()
                   && child_seq[prefix] == node_seq[prefix]) {
                ++prefix;
            }
            size_t suffix = 0;
            while (suffix < child_seq.size() - prefix && suffix < node_seq.size() - prefix
                   && child_seq[child_seq.size() - 1 - suffix]
                          == node_seq[node_seq.size() - 1 - suffix]) {
                ++suffix;
            }
            step_edits.push_back(StepEdit{i + k, prefix, node_seq.size() - prefix - suffix,
                                          child_seq.substr(prefix,
                                                           child_seq.size() - prefix - suffix),
                                          edit.key});
        }
    }
    std::sort(step_edits.begin(), step_edits.end(), [](const StepEdit& a, const StepEdit& b) {
        // Right to left; at one offset, a replacement before an insertion in front of it.
        return std::tie(a.step, b.at, b.removed, a.key) < std::tie(b.step, a.at, a.removed, b.key);
    });
    for (size_t e = 0; e < step_edits.size(); ++e) {
        const StepEdit& edit = step_edits[e];
        // One already made to its right that it would overlap wins.
        if (e > 0 && step_edits[e - 1].step == edit.step
            && edit.at + edit.removed > step_edits[e - 1].at
            && !(edit.removed == 0 && step_edits[e - 1].removed == 0
                 && edit.at == step_edits[e - 1].at)) {
            continue;
        }
        string& spelled = steps[edit.step];
        if (edit.at + edit.removed > spelled.size()) {
            continue;
        }
        if (spelled.compare(edit.at, edit.removed, edit.added) != 0) {
            spelled.replace(edit.at, edit.removed, edit.added);
            changed = true;
        }
    }
    if (!changed) {
        walk.clear();
        steps.clear();
    }
    return changed;
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
                graph, parent.children, parent.travs, linker.chosen_genotype(parent),
                parent.ref_trav_idx);
            for (size_t ci : children) {
                StagedSite& child = pending[ci];
                if (child.dropped) {
                    continue;
                }
                const bool was = child.reported_inline;
                child.reported_inline = parent.reported_inline
                                        || (child.in_chain
                                            && block_records.chain_reported_inline(graph, ctx,
                                                                                   child.chain));
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
            bytes += rec.travs.capacity() * sizeof(Traversal);
            for (const Traversal& t : rec.travs) {
                visits += t.size();
                bytes += t.capacity() * sizeof(handle_t);
            }
            bytes += rec.children.chains.capacity() * sizeof(ChildChain)
                     + rec.children.entries.capacity() * sizeof(rec.children.entries[0]);
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
                bytes += rl->scored_traversals.capacity() * sizeof(Traversal)
                         + rl->allele_support.capacity() * sizeof(double);
                if (rl->alt_ploidy_info != nullptr) {
                    // The alternate answer is kept too, with all its parts.
                    const auto& alt = *rl->alt_ploidy_info;
                    bytes += alt.scored_traversals.capacity() * sizeof(Traversal)
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
        unexplained += (pc.ploidy == 1)
                       ? (pc.hap_first == LinkageModel::NO_HAPLOTYPE
                          && pc.hap_second == LinkageModel::NO_HAPLOTYPE)
                       : (pc.hap_first == LinkageModel::NO_HAPLOTYPE
                          || pc.hap_second == LinkageModel::NO_HAPLOTYPE);
        order_arbitrary += pc.order_arbitrary;
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
        ref_trav_idx >= (int)site.travs.size() || trav_idx >= (int)site.travs.size()) {
        return false;
    }
    // An allele with step sequences spells something its walk does not, so sharing the
    // reference's route does not make it the reference.
    if (own_sequences(&site.sequences, trav_idx) != nullptr
        || own_sequences(&site.sequences, ref_trav_idx) != nullptr) {
        return false;
    }
    return symbolically_equal(graph, site.travs[trav_idx], site.travs[ref_trav_idx],
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
        && !splice_children.count(record_key)) {
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
}
