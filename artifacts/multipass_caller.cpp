#include <atomic>
#include <cstdlib>
#include <fstream>
#include <limits>

#include <omp.h>

#include "multipass_caller.hpp"
#include "../src/utility.hpp"

namespace vg::multipass {

MultiPassCaller::MultiPassCaller(const PathPositionHandleGraph& graph,
                                 ReadLikelihoodSnarlCaller& genotyper,
                                 const SnarlDecomposition& decomposition,
                                 SortedVCFWriter& output, const string& sample_name,
                                 SiteTraversalFinder& traversal_finder,
                                 ReferencePathTable ref_paths, bool genotype_snarls,
                                 const pair<size_t, size_t>& allele_length_range) :
    graph(graph),
    snarl_caller(genotyper),
    sites(decomposition, graph),
    walker(sites, graph),
    output(output),
    site_emitter(*this, output),
    sample_name(sample_name),
    ref_paths(std::move(ref_paths)),
    genotype_snarls(genotype_snarls),
    splicer(graph, sites, phase_table),
    candidates(graph, traversal_finder, this->ref_paths, allele_length_range),
    record_renderer(phase_table, read_strands, linker, anchor_collector, fill_anchors)
{
    site_emitter.genotype_snarls = genotype_snarls;
    // Every call info this caller writes is the read-likelihood genotyper's, which writes GL in
    // colexicographic order.
    site_emitter.gl_layout = GLLayout::Colexicographic;
    install_widgets();
}

void MultiPassCaller::call(GraphCaller::RecurseType recurse_type,
                           const function<void()>& release_reads_now) {
    release_reads = release_reads_now;
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
                return site_name(graph, start, end, output.translation());
            },
            .edit_source = off_panel.source(),
            .read_source = read_source,
            .edit_outcomes = off_panel.outcomes(),
            .sites = &sites,
        },
        TreeGenotyper::Options{
            .nested_calling = nested_calling,
            .off_reference = off_reference_nesting,
            .group_candidates = group_candidates,
            .edit_bar = off_panel.bar(),
        });

    // Place the fill points before the direct pass, and have the read source pass each window it
    // fetches whole to the fill anchors, which place the reads crossing that window's points while
    // the window is in memory. The fill windows must be the source's windows, or
    // `resolve_window_reads` takes none of them; `snarl_batch_window` is the source's window (see
    // `set_snarl_batching`). The observer adds nothing to the source's tallies of read starts or to
    // its off-panel candidates, so it cannot change a genotype.
    const size_t fill_window = snarl_batch_window > 0 ? snarl_batch_window : FILL_READ_WINDOW;
    const bool fill_reads = fill_anchors.is_enabled() && read_source != nullptr;
    if (fill_reads) {
        fill_anchors.place_points(graph, fill_window);
        read_source->set_window_observer(
            [this](nid_t first, nid_t last, const SiteReadSource::WindowReads& reads) {
                fill_anchors.resolve_window_reads(first, last, reads);
            });
    }

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
    if (fill_reads) {
        // No fetch is running between the passes, so the observer can be removed. Then place the
        // reads at the fill points of the windows the direct pass never fetched whole, while the
        // read source is open. No pass after the direct pass reads reads, and this query adds
        // nothing to the source's tallies of read starts or to its off-panel candidates, so it
        // cannot change a genotype.
        read_source->set_window_observer(nullptr);
        fill_anchors.place(sites);
        fill_anchors.resolve_reads(*read_source);
    }
    release_reads();

    // Round 1's linkage pass, then read phasing and any further rounds, then the render.
    run_linkage_pass();
    phase_and_regenotype();
    render_retained_records();
    // Anchors are collected while records are built, so they are written afterwards. Does nothing
    // unless anchors are on. The render marked the called walks across the fill points, so the
    // fill anchors are made now.
    if (fill_anchors.is_enabled()) {
        fill_anchors.add_anchors(anchor_collector);
    }
    anchor_collector.write(sample_name);
    if (fill_anchors.is_enabled() && anchor_collector.is_enabled()) {
        fill_anchors.report(cerr);
    }

    // The reports, and the mosaic, need to know which sites have a line, so they come last.
    finalise_linkage_outputs();
    off_panel.report(tree_genotyper.edit_nodes_after_pass());
    if (phase_declined.load() > 0 || quality_declined.load() > 0) {
        cerr << "[vg call] linkage: " << phase_declined.load()
             << " phases refused by the record they were rendered onto, and "
             << quality_declined.load() << " quality rewrites refused" << endl;
    }
    block_records.report();
}

string MultiPassCaller::record_id_of_site(const net_handle_t& site) const {
    if (sites.is_edit(site)) {
        auto found = edit_record_ids.find(sites.canonical(site));
        return found != edit_record_ids.end() ? found->second : string();
    }
    return site_name(graph, sites.get_site_start(site), sites.get_site_end(site),
                     output.translation());
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

MultiPassCaller::StagedSiteEmitter::StagedSiteEmitter(MultiPassCaller& caller,
                                                      SortedVCFWriter& output) :
    SiteEmitter(output), caller(caller) {
    // Nothing to do
}

bool MultiPassCaller::StagedSiteEmitter::same_as_reference(const CalledSite& site, int trav) {
    return caller.nested_calling
           && caller.is_symbolically_reference(caller.rendered_site(site), trav,
                                               site.ref_trav_idx);
}

int64_t MultiPassCaller::StagedSiteEmitter::phase(const CalledSite& site,
                                                  const vector<int>& genotype,
                                                  const map<int, int>& trav_to_allele,
                                                  string& gt) {
    return caller.phase_record_genotype(caller.rendered_site(site).record_key, genotype,
                                        trav_to_allele, gt);
}

int MultiPassCaller::StagedSiteEmitter::write_blocks(const PathPositionHandleGraph& graph,
                                                     const CalledSite& site,
                                                     const SiteRecord& record) {
    const StagedSite& staged = caller.rendered_site(site);
    if (staged.flat) {
        // The flat records stand in for the site's record or blocks.
        caller.block_records.write_flat(
            graph, *staged.flat, output.sample_name(), output.translation(), record, gl_layout,
            [this](vcflib::Variant& line, size_t block) {
                caller.finish_moved_record(line);
                return output.add_variant(line, block);
            });
    }
    // A site whose lines are kept out is still cut into blocks, so that it is reported as it is
    // without flat records.
    const bool withheld = caller.withholds_lines(staged);
    // The blocks spell the alleles the site record spells: the spliced ones, if the site has
    // them, or else its own, by their step sequences where they have them.
    AlleleSequences own_steps;
    if (!staged.spliced && caller.block_records.is_enabled()) {
        own_steps = staged.alleles.all_step_sequences(graph);
    }
    return caller.block_records.write(
        graph, staged.children, staged.spliced ? staged.spliced->travs : staged.alleles.walks,
        site.genotype, site.ref_trav_idx, output.sample_name(), output.translation(), record,
        gl_layout, genotype_snarls,
        [this, withheld](vcflib::Variant& line, size_t block) {
            caller.finish_moved_record(line);
            return withheld || output.add_variant(line, block);
        },
        staged.spliced ? &staged.spliced->sequences : &own_steps, !staged.spliced);
}

bool MultiPassCaller::StagedSiteEmitter::withholds_record(const CalledSite& site) {
    return caller.withholds_lines(caller.rendered_site(site));
}

void MultiPassCaller::StagedSiteEmitter::finish_record(const CalledSite&, vcflib::Variant& record) {
    caller.finish_moved_record(record);
}

void MultiPassCaller::StagedSiteEmitter::on_site_done(const CalledSite& site,
                                                      const map<int, int>& trav_to_allele,
                                                      bool has_line) {
    // The linkage model gets the site whether or not it has a line. A parent written as the
    // reference still has two alleles, which differ only inside its children, and the children
    // need them to know which strand carries the chain. In VCF allele numbering such a parent is
    // 0/0; only in traversal space is it heterozygous.
    if (!caller.linker.enabled()) {
        return;
    }
    // The site was recorded when it was genotyped. What remains is the traversal-to-VCF-allele
    // map, which depends on the alleles the record chose, and whether a line was written. A
    // site written as blocks gives an empty map, since each block numbers its own alleles, and
    // the linkage model then gets an empty map too.
    vector trav_to_allele_vec(trav_to_allele.empty() ? 0 : site.traversal_count, -1);
    for (const auto& [trav, allele] : trav_to_allele) {
        if (trav >= 0 && (size_t)trav < trav_to_allele_vec.size()) {
            trav_to_allele_vec[trav] = allele;
        }
    }
    caller.linker.collector()->set_allele_map(caller.rendered_site(site).record_key,
                                              trav_to_allele_vec, has_line);
}

const StagedSite& MultiPassCaller::rendered_site(const CalledSite& site) const {
    return *rendered_sites.at(record_key_of(*site.id));
}

void MultiPassCaller::set_flatten_records(bool on) {
    flatten_records = on;
    block_records.set_flat_check(on && getenv("VG_CHECK_FLAT_RECORDS") != nullptr);
}

bool MultiPassCaller::can_flatten(const StagedSite& site) const {
    // The strands are those of the site's phase, which must carry the chosen genotype.
    const vector<int> genotype = linker.chosen_genotype(site);
    const LinkageCollector::PhaseCall* phase = phase_table.rendered(site.record_key);
    if (genotype.empty() || phase == nullptr) {
        block_records.count_flat_refusal(genotype.empty() ? FlatRefusal::NoGenotype
                                                          : FlatRefusal::NoPhase);
        return false;
    }
    vector<int> carried;
    if (phase->ploidy == genotype.size()) {
        for (int strand = 0; strand < (int)genotype.size(); ++strand) {
            const optional<int> trav = splicer.strand_traversal(site, strand);
            if (!trav) {
                break;
            }
            carried.push_back(*trav);
        }
    }
    vector<int> chosen = genotype;
    std::ranges::sort(chosen);
    std::ranges::sort(carried);
    if (carried != chosen) {
        block_records.count_flat_refusal(FlatRefusal::PhaseDiffers);
        return false;
    }
    // A record is placed by its offset from the first base of the site's reference allele, so
    // the resolved reference must be that allele.
    const optional<string> reference = splicer.reference_sequence(site);
    if (!reference || *reference != site.alleles.sequence(graph, site.alleles.ref_allele)) {
        block_records.count_flat_refusal(FlatRefusal::ReferenceDiffers);
        return false;
    }
    return true;
}

bool MultiPassCaller::orders_strands(const LinkageCollector::PhaseCall& phase) const {
    // Where neither the panel nor the reads order the pair, the record claims no phase. A site
    // with sites staged under it keeps its phase, since its one-copy children's `a|.` GTs name
    // its strands; so does a one-copy chain, whose GT shape says which strand it is on.
    return !(phase.ploidy == 2 && phase.order_arbitrary
             && !read_strands.read_ordered().count(phase.record_key)
             && !splicer.has_children(phase.record_key));
}

vector<MultiPassCaller::FlatRunPlan> MultiPassCaller::flat_runs(
    const vector<StagedSite*>& roots, const vector<char>& can) const {
    const size_t none = std::numeric_limits<size_t>::max();
    auto phase_of = [&](size_t i) { return phase_table.rendered(roots[i]->record_key); };
    auto reference_walk = [&](size_t i) -> const Traversal& {
        return roots[i]->alleles.walks[roots[i]->alleles.ref_allele];
    };
    auto bases_of = [&](const Traversal& walk) {
        size_t bases = 0;
        for (const handle_t& step : walk) {
            bases += graph.get_length(step);
        }
        return bases;
    };
    // A site can share a run where its record claims the phase its strands are spelled by.
    auto shares = [&](size_t i) {
        return can[i] && (phase_of(i)->trav_first == phase_of(i)->trav_second
                          || orders_strands(*phase_of(i)));
    };
    auto same_phase = [&](size_t a, size_t b) {
        const LinkageCollector::PhaseCall& first = *phase_of(a);
        const LinkageCollector::PhaseCall& second = *phase_of(b);
        return first.contig == second.contig && first.phase_set == second.phase_set
               && first.ploidy == second.ploidy;
    };
    // The reference path reads site `b` right after the nodes `nodes` that follow site `a`, so
    // that an offset along the run is an offset along the path.
    auto path_follows = [&](size_t a, size_t b, const vector<handle_t>& nodes) {
        const ReferencePlacement& first = roots[a]->reference;
        const ReferencePlacement& second = roots[b]->reference;
        if (first.no_reference || second.no_reference || first.path_name != second.path_name
            || first.offset != second.offset) {
            return false;
        }
        const auto here = get_ref_interval(graph, sites.get_site_start(roots[a]->net),
                                           sites.get_site_end(roots[a]->net), first.path_name);
        const auto there = get_ref_interval(graph, sites.get_site_start(roots[b]->net),
                                            sites.get_site_end(roots[b]->net), first.path_name);
        return get<0>(here) >= 0 && get<0>(there) >= 0 && !get<2>(here) && !get<2>(there)
               && get<0>(there) + (int64_t)graph.get_length(nodes.back())
                      == get<0>(here) + (int64_t)bases_of(reference_walk(a))
                             - (int64_t)graph.get_length(nodes.front())
                             + (int64_t)bases_of(nodes);
    };

    unordered_map<handle_t, size_t> by_start;
    unordered_map<nid_t, vector<size_t>> edits_on;
    for (size_t i = 0; i < roots.size(); ++i) {
        if (!sites.is_edit(roots[i]->net)) {
            by_start[sites.get_site_start(roots[i]->net)] = i;
            continue;
        }
        for (const handle_t& step : roots[i]->alleles.walks.front()) {
            edits_on[graph.get_id(step)].push_back(i);
        }
    }
    // The snarl after each snarl in its run, with the nodes and the edit sites between them.
    vector<size_t> next(roots.size(), none);
    vector<char> follows(roots.size(), 0);
    vector<vector<handle_t>> nodes_after(roots.size());
    vector<vector<size_t>> edits_after(roots.size());
    for (size_t i = 0; i < roots.size(); ++i) {
        if (sites.is_edit(roots[i]->net) || !shares(i)) {
            continue;
        }
        vector<handle_t> nodes;
        size_t bases = 0;
        size_t found = none;
        if (!walk_to_next_site(
                graph, sites.get_site_end(roots[i]->net), flatten_distance,
                [&](const handle_t& at) {
                    const auto starting = by_start.find(at);
                    if (starting == by_start.end() || starting->second == i) {
                        return false;
                    }
                    found = starting->second;
                    return true;
                },
                nodes, bases)
            || !shares(found) || !same_phase(i, found) || follows[found]
            || !path_follows(i, found, nodes)) {
            continue;
        }
        // Every edit site in the nodes between the two has to be in the run, with all its nodes
        // there, or its changes would lie where the run's records are.
        unordered_set<nid_t> between;
        for (const handle_t& node : nodes) {
            between.insert(graph.get_id(node));
        }
        vector<size_t> edits;
        bool all_in = true;
        for (const handle_t& node : nodes) {
            const auto on_node = edits_on.find(graph.get_id(node));
            for (size_t e = 0; on_node != edits_on.end() && e < on_node->second.size(); ++e) {
                const size_t edit = on_node->second[e];
                all_in = all_in && can[edit] && shares(edit) && same_phase(i, edit)
                         && std::ranges::all_of(roots[edit]->alleles.walks.front(),
                                                [&](const handle_t& step) {
                                                    return between.contains(graph.get_id(step));
                                                });
                if (std::ranges::find(edits, edit) == edits.end()) {
                    edits.push_back(edit);
                }
            }
        }
        if (!all_in) {
            continue;
        }
        next[i] = found;
        follows[found] = 1;
        nodes_after[i] = std::move(nodes);
        edits_after[i] = std::move(edits);
    }

    vector<FlatRunPlan> runs;
    vector<char> in_run(roots.size(), 0);
    for (size_t i = 0; i < roots.size(); ++i) {
        if (!can[i] || follows[i] || sites.is_edit(roots[i]->net)) {
            continue;
        }
        FlatRunPlan run;
        size_t bases = 0;
        for (size_t at = i, before = none; at != none; before = at, at = next[at]) {
            // The bases the site adds to the run's reference, after the node it starts at.
            const size_t own = bases_of(reference_walk(at));
            if (!run.sites.empty()) {
                const size_t more = bases_of(nodes_after[before]) + own
                                    - graph.get_length(nodes_after[before].front())
                                    - graph.get_length(nodes_after[before].back());
                if (bases + more > flat_run_max_bases) {
                    // The site would make the run too long, so it starts the next run.
                    runs.push_back(std::move(run));
                    run = FlatRunPlan();
                    bases = 0;
                } else {
                    run.junctions.push_back(nodes_after[before]);
                    for (const size_t edit : edits_after[before]) {
                        run.between.push_back(roots[edit]);
                        in_run[edit] = 1;
                    }
                    bases += more;
                }
            }
            if (run.sites.empty()) {
                bases = own;
            }
            run.sites.push_back(roots[at]);
        }
        runs.push_back(std::move(run));
    }
    for (size_t i = 0; i < roots.size(); ++i) {
        if (can[i] && !in_run[i] && sites.is_edit(roots[i]->net)) {
            runs.push_back(FlatRunPlan{{roots[i]}, {}, {}});
        }
    }
    return runs;
}

bool MultiPassCaller::align_flat_run(const FlatRunPlan& plan, FlatRefusal& refused) const {
    // The walk the run's reference spells, and the offset in it of each site's first base.
    Traversal along;
    vector<size_t> origins;
    size_t bases = 0;
    auto append = [&](const Traversal& walk, size_t from) {
        for (size_t k = from; k < walk.size(); ++k) {
            along.push_back(walk[k]);
            bases += graph.get_length(walk[k]);
        }
    };
    for (size_t m = 0; m < plan.sites.size(); ++m) {
        if (m > 0) {
            append(plan.junctions[m - 1], 1);
        }
        const CandidateAlleles& alleles = plan.sites[m]->alleles;
        // A site after the first starts at the last node so far.
        origins.push_back(m > 0 ? bases - graph.get_length(along.back()) : 0);
        append(alleles.walks[alleles.ref_allele], m > 0 ? 1 : 0);
    }
    string reference;
    for (const handle_t& step : along) {
        reference += graph.get_sequence(step);
    }
    vector<const StagedSite*> on_walk(plan.sites.begin(), plan.sites.end());
    on_walk.insert(on_walk.end(), plan.between.begin(), plan.between.end());
    const size_t ploidy = phase_table.rendered(plan.sites.front()->record_key)->ploidy;
    vector<string> strands;
    for (int strand = 0; strand < (int)ploidy; ++strand) {
        optional<string> sequence = splicer.sequence_along(along, on_walk, strand, ploidy == 1);
        if (!sequence) {
            refused = FlatRefusal::Unspelled;
            return false;
        }
        strands.push_back(std::move(*sequence));
    }
    // Guide each strand's alignment by the sites: the strand reaches the first and the last base
    // between a site's bound nodes where the reference does, moved by what the strand's alleles
    // at the sites before, and at that site, add to or take from the reference's length.
    vector<vector<pair<size_t, size_t>>> guides(ploidy);
    for (int strand = 0; strand < (int)ploidy; ++strand) {
        int64_t moved = 0;
        for (size_t m = 0; m < plan.sites.size(); ++m) {
            const StagedSite& site = *plan.sites[m];
            const Traversal& walk = site.alleles.walks[site.alleles.ref_allele];
            size_t length = 0;
            for (const handle_t& step : walk) {
                length += graph.get_length(step);
            }
            const size_t first = origins[m] + graph.get_length(walk.front());
            const size_t last = origins[m] + length - graph.get_length(walk.back());
            const optional<string> spelled = splicer.strand_sequence(site, strand);
            auto at = [&](size_t offset) {
                return make_pair(offset, (size_t)std::max<int64_t>(0, (int64_t)offset + moved));
            };
            guides[strand].push_back(at(first));
            moved += spelled ? (int64_t)spelled->size() - (int64_t)length : 0;
            guides[strand].push_back(at(last));
        }
    }
    const shared_ptr<FlatRun> run = block_records.align_flat(
        std::move(reference), std::move(strands), refused, guides);
    if (!run) {
        return false;
    }
    run->id = plan.sites.size() == 1
                  ? plan.sites.front()->id
                  : site_name(graph, sites.get_site_start(plan.sites.front()->net),
                              sites.get_site_end(plan.sites.back()->net), output.translation());
    run->check_waiting = plan.sites.size();
    // Whether the called genotype of a site of the run, or of a site under one, has an allele
    // built from the reads.
    auto note_consensus = [&](const StagedSite& at) {
        run->consensus_called =
            run->consensus_called
            || std::ranges::any_of(linker.chosen_genotype(at),
                                   [&](int allele) { return at.is_consensus_allele(allele); });
    };
    for (const StagedSite* site : on_walk) {
        note_consensus(*site);
        splicer.for_each_under(site->record_key, note_consensus);
    }
    // A site writes the records from its first base up to that of the next site.
    size_t cut = 0;
    for (size_t m = 0; m < plan.sites.size(); ++m) {
        FlatSite flat{run, origins[m], {}, cut, 0};
        while (cut < run->cuts.size()
               && (m + 1 == plan.sites.size()
                   || (size_t)run->cuts[cut].begin() < origins[m + 1])) {
            ++cut;
        }
        flat.cut_count = cut - flat.first_cut;
        for (int strand = 0; strand < (int)ploidy; ++strand) {
            flat.travs.push_back(*splicer.strand_traversal(*plan.sites[m], strand));
        }
        plan.sites[m]->flat = std::move(flat);
    }
    block_records.count_flat_run(*run, plan.sites.size());
    return true;
}

void MultiPassCaller::plan_flat_records() {
    flat_withheld.clear();
    if (!flatten_records) {
        return;
    }
    vector<StagedSite*> roots;
    for (size_t q = 0; q < staged_sites.queue_count(); ++q) {
        for (StagedSite& site : staged_sites.queue(q)) {
            if (site.placement.parent_record_key == 0) {
                roots.push_back(&site);
            }
        }
    }
    vector<char> can(roots.size(), 0);
#pragma omp parallel for schedule(dynamic, 64)
    for (size_t i = 0; i < roots.size(); ++i) {
        can[i] = can_flatten(*roots[i]);
    }
    const vector<FlatRunPlan> runs = flat_runs(roots, can);
    // The sites under the sites of each run that is written, and the edit sites between them.
    vector<vector<size_t>> under(runs.size());
    vector<vector<size_t>> between(runs.size());
#pragma omp parallel for schedule(dynamic, 1)
    for (size_t r = 0; r < runs.size(); ++r) {
        auto write = [&](const FlatRunPlan& plan) {
            FlatRefusal refused{};
            if (!align_flat_run(plan, refused)) {
                if (plan.sites.size() > 1) {
                    block_records.count_flat_run_refusal(refused);
                } else {
                    block_records.count_flat_refusal(refused);
                }
                return false;
            }
            for (const StagedSite* site : plan.sites) {
                splicer.for_each_under(site->record_key, [&](const StagedSite& nested) {
                    under[r].push_back(nested.record_key);
                });
            }
            for (const StagedSite* site : plan.between) {
                between[r].push_back(site->record_key);
            }
            return true;
        };
        if (!write(runs[r]) && runs[r].sites.size() > 1) {
            // Each site of the run, and each edit site between them, is a run of its own.
            for (StagedSite* site : runs[r].sites) {
                write(FlatRunPlan{{site}, {}, {}});
            }
            for (StagedSite* site : runs[r].between) {
                write(FlatRunPlan{{site}, {}, {}});
            }
        }
    }
    size_t nested = 0;
    size_t edits = 0;
    for (size_t r = 0; r < runs.size(); ++r) {
        for (const size_t key : under[r]) {
            // Only a site in a render queue has a line to leave out.
            if (rendered_sites.contains(key)) {
                nested += (size_t)flat_withheld.insert(key).second;
            }
        }
        for (const size_t key : between[r]) {
            edits += (size_t)flat_withheld.insert(key).second;
        }
    }
    block_records.count_flat_withheld(nested, edits);
}

RecordRenderer::WrittenAlleles MultiPassCaller::write_site(StagedSite& site,
                                                           const vector<int>& genotype) {
    if (nested_calling) {
        block_records.count_site(site.children, site.alleles.walks, site.alleles.ref_allele);
    }
    // A parent spells its children's chosen alleles.
    if (AlleleSplicer::SplicedAlleles spliced; splicer.splice(site, spliced)) {
        site.spliced = std::move(spliced);
    }
    // The alleles the record spells, and their step sequences.
    const vector<Traversal>& spelled = site.spliced ? site.spliced->travs : site.alleles.walks;
    AlleleSequences own_steps;
    if (!site.spliced) {
        own_steps = site.alleles.all_step_sequences(graph);
    }
    const AlleleSequences& steps = site.spliced ? site.spliced->sequences : own_steps;
    // Whether the called genotype has an allele built from the reads, here or, for a site
    // written as flat records, anywhere in its run.
    const bool consensus_called =
        read_consensus.enabled()
        && (std::ranges::any_of(genotype,
                                [&](int allele) { return site.is_consensus_allele(allele); })
            || (site.flat && site.flat->run->consensus_called));
    RecordRenderer::WrittenAlleles written;
    bool has_line = false;
    site_emitter.emit(
        graph,
        CalledSite{
            .start = sites.get_site_start(site.net),
            .end = sites.get_site_end(site.net),
            .ref_path_name = site.reference.path_name,
            .ref_offset = site.reference.offset,
            .genotype = genotype,
            .ref_trav_idx = site.alleles.ref_allele,
            .traversal_count = site.alleles.walks.size(),
            .ploidy = site.ploidy,
            .id = &site.id,
        },
        [&](size_t allele) {
            return make_pair(allele_sequence(graph, spelled[allele], own_sequences(&steps, allele)),
                             visits_of(graph, site.alleles.walks[allele]));
        },
        // A written allele's scored allele is its traversal, since the site's walks are the ones
        // it was scored with, in order. Two alleles can share a walk when one has step sequences,
        // so they are told apart by index rather than by walk.
        [&](const vector<int>& allele_traversals, const vector<int>& vcf_genotype,
            vcflib::Variant& record) {
            written = {allele_traversals, vcf_genotype};
            snarl_caller.write_vcf_fields(allele_traversals, vcf_genotype, site.call_info.get(),
                                          sample_name, record);
            if (sites.is_edit(site.net)) {
                record.infoFlags["OP"] = true;
            }
            if (consensus_called) {
                record.infoFlags["RC"] = true;
            }
        },
        &has_line);
    site.spliced.reset();
    site.flat.reset();
    return has_line ? written : RecordRenderer::WrittenAlleles();
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

/// The header line of INFO/RC.
static const char* const RC_HEADER_LINE =
    "##INFO=<ID=RC,Number=0,Type=Flag,Description=\"Read consensus (--read-consensus): the "
    "called genotype of the record's snarl has an allele that is the consensus of the reads "
    "spanning the snarl. Its AT is the walk the strand is called as, whose bases the allele "
    "changes, or the snarl's two bound nodes alone where the allele is given as its bases. The "
    "record's site can be a run of neighbouring snarls of one chain made one site, which is "
    "named by the run's two outer bound nodes. On a flat record (--flatten-records), the called "
    "genotype of a top-level snarl of the record's run, or of a snarl nested in one, has such "
    "an allele\">";

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

/// The header line of INFO/FL.
static const char* const FL_HEADER_LINE =
    "##INFO=<ID=FL,Number=0,Type=Flag,Description=\"Flat record (--flatten-records): a "
    "base-level difference between the strands and the reference over a run of neighbouring "
    "top-level sites of one chain, or over one top-level site, where each strand is spelled "
    "with the called alleles of the sites of the run, of the sites nested in them and of the "
    "edit sites between them, which have no records of their own. Each strand is aligned to "
    "the reference by itself. A run of mismatches is one record, and so is a gap, with the "
    "reference base before it. A difference both strands have is one record with GT 1|1, and "
    "any other is a record of its own strand with GT 1|0 or 0|1, where 0 says only that the "
    "other strand lacks this difference: two records can cover the same reference bases, one "
    "for each strand. The ID is that of the one site, or names the start bound node of the "
    "run's first site and the end bound node of its last, with _ and the record's index among "
    "the run's records (INFO/SB) appended. AT is . for each allele, since no walk of the "
    "graph spells them. QUAL and the FORMAT fields other than GT are those of the site of "
    "the run that the record starts in, or after, with AD and GL taken through the site "
    "allele each strand carries\">";

void MultiPassCaller::add_header_lines() {
    if (emit_phasing) {
        // FORMAT/PS is the VCF phase set, which phasing tools read. It is unrelated to INFO/PS,
        // vg's parent-snarl field; the two are in different namespaces, so both are legal, and
        // their descriptions say which is which.
        output.add_format_line(
            "##FORMAT=<ID=PS,Number=1,Type=Integer,Description=\"Phase set: the phase of a "
            "genotype is comparable only with others carrying the same PS. One phase set per "
            "chain, so blocks are chromosome-scale -- much longer than a read-based phaser "
            "gives, because the phase comes from the haplotype panel rather than from reads "
            "spanning consecutive sites. Not the INFO/PS emitted under -A, which is a parent "
            "snarl pointer\">");
    }
    if (off_panel.is_enabled()) {
        output.add_info_line(OP_HEADER_LINE);
    }
    if (read_consensus.enabled()) {
        output.add_info_line(RC_HEADER_LINE);
    }
    if (block_records.is_enabled() || flatten_records) {
        output.add_info_line(SB_HEADER_LINE);
    }
    if (flatten_records) {
        output.add_info_line(FL_HEADER_LINE);
    }
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

ReadConsensusStep::Parts MultiPassCaller::read_consensus_parts() {
    return ReadConsensusStep::Parts{
        .graph = &graph,
        .sites = &sites,
        .genotyper = &snarl_caller,
        .tree_genotyper = &tree_genotyper,
        .linker = &linker,
        .staged_sites = &staged_sites,
        .phase_table = &phase_table,
        .read_strands = &read_strands,
    };
}

void MultiPassCaller::add_read_consensus_alleles() {
    if (!read_consensus.enabled() || !linker.enabled()) {
        return;
    }
    // The step asks which strand each read is on, so the read strand table needs its sums, as
    // the render builds them. No round of re-genotyping has run, so the calibration is the
    // temper re-genotyping was given, or else one fitted here.
    const RegenotypeParams& regenotyping = rescorer.params();
    read_strands.build_lambda(
        phase_table.calls(),
        [&](const string& contig, size_t phase_set) { return phase_set_id(contig, phase_set); },
        StrandCalibration{std::max(regenotyping.temper, 0.0),
                          regenotyping.ceiling < 0.0 ? 1.0 : regenotyping.ceiling},
        regenotyping);
    const size_t changed = read_consensus.run(read_consensus_parts());
    // The step fetched reads, which no later pass does.
    release_reads();
    if (changed > 0) {
        // The sites that gained an allele have new calls, so choose the genotypes and phase the
        // reads again, as a re-genotyping round does after its correction.
        rerun_linkage_pass();
        phase_reads();
    }
}

void MultiPassCaller::phase_and_regenotype() {
    // Round 1 ends with read phasing. Its linkage pass has already run.
    phase_reads();
    add_read_consensus_alleles();
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
    if (read_consensus.enabled()) {
        read_consensus.count_called(read_consensus_parts());
        read_consensus.report(cerr);
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
        [&](StagedSite& site, const vector<int>& genotype) { return write_site(site, genotype); },
        show_progress,
        [&]() {
            splicer.index(staged_sites, emit_phasing);
            rendered_sites = staged_sites.by_key();
            plan_flat_records();
            if (mosaic_writer.is_exact()) {
                // An exact mosaic needs each strand as the records spell it at every top-level
                // site, and the sites are freed as they are rendered. The phase is final here.
                StrandAlleles& kept = mosaic_writer.strand_alleles();
                staged_sites.for_each([&](const StagedSite& site) {
                    if (site.placement.parent_record_key != 0) {
                        return;
                    }
                    const bool edit_site = sites.is_edit(site.net);
                    kept.unphased_parents += (size_t)(splicer.has_children(site.record_key)
                                                      && !phase_table.rendered(site.record_key));
                    for (int strand = 0; strand < 2; ++strand) {
                        Traversal walk;
                        StepSequences steps;
                        const optional<int> trav = splicer.strand_traversal(site, strand);
                        if (!trav) {
                            continue;
                        }
                        if (!edit_site) {
                            splicer.strand_walk(site, strand, walk, steps);
                            kept.add_site(graph, site.record_key, strand, walk, steps);
                        } else if (!site.alleles.spelled_by_walk(*trav)) {
                            kept.add_edit_site(graph, site.record_key, strand,
                                               site.alleles.walks[*trav],
                                               site.alleles.edits[*trav]);
                        }
                    }
                });
            }
        });
    rendered_sites.clear();
    flat_withheld.clear();
    staged_sites.for_each([&](StagedSite& site) {
        if (sites.is_edit(site.net)) {
            edit_record_ids.emplace(sites.canonical(site.net), site.id);
        }
    });
    splicer.clear();
    staged_sites.release_unrendered();
}

size_t MultiPassCaller::rederive_reported_inline() {
    if (!staged_sites.has_children()) {
        return 0;
    }
    vector<StagedSite>& pending = staged_sites.nested();
    size_t rederived = 0;
    const unordered_map<size_t, StagedSite*> record_by_key = staged_sites.by_key();
    // Count the chains held back from here, so that the report gives the chains held back now.
    block_records.restart_inline_count();
    // Whether the record that spells each chain held back spells one traversal for both strands.
    unordered_map<size_t, bool> spelled_once;
    // Visit parents before their children, so that a chain inherits its parent's final flag.
    staged_sites.for_each_parent_top_down(record_by_key, [&](const StagedSite& parent,
                                                             const vector<size_t>& children) {
        if (parent.dropped) {
            // Its children were dropped with it.
            return;
        }
        const vector<int> parent_genotype = linker.chosen_genotype(parent);
        const BlockRecordWriter::ChainInlineContext ctx = block_records.chain_inline_context(
            graph, parent.children, parent.alleles.walks, parent_genotype,
            parent.alleles.ref_allele);
        // A record spells a traversal that both strands carry once, as the first strand has it.
        // A chain whose strands differ is then not spelled for both, and keeps its line.
        const bool once = parent.placement.reported_inline
                              ? spelled_once[parent.record_key]
                              : parent_genotype.size() == 2
                                    && parent_genotype[0] == parent_genotype[1];
        for (size_t ci : children) {
            StagedSite& child = pending[ci];
            if (child.dropped) {
                continue;
            }
            const vector<int> genotype = linker.chosen_genotype(child);
            const bool differs = genotype.size() == 2 && genotype[0] != genotype[1];
            const bool was = child.placement.reported_inline;
            child.placement.reported_inline =
                !(once && differs)
                && (parent.placement.reported_inline
                    || block_records.chain_reported_inline(
                        graph, ctx, chain_holding(sites, graph, child.net)));
            spelled_once[child.record_key] = once;
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
        // An exact mosaic also covers the genotyped sites with no record, since a strand's
        // haplotype can take another allele than the called one there too.
        mosaic_writer.write(written, panel_lookup, sample_name, &phase_table.calls(), &graph);
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
    if (!orders_strands(phase)) {
        return -1;
    }
    if (!phased_gt(phase, site_genotype, trav_to_allele, gt)) {
        ++phase_declined;
        return -1;
    }
    return (int64_t)phase.phase_set;
}

}
