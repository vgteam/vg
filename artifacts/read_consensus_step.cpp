#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <fstream>

#include "read_consensus_step.hpp"
#include "site_values.hpp"

namespace vg::multipass {

using namespace std;

struct ReadConsensusStep::Finding {
    /// What the step found on one strand.
    struct Strand {
        /// The strand's spanning reads, and how many of them disagree.
        size_t reads = 0;
        size_t disagreeing = 0;
        bool candidate = false;
        bool built = false;
        bool unchanged = false;
        bool already_allele = false;
        bool covered = false;
        /// The bases the consensus puts between the bounds, where it is to become an allele.
        optional<string> allele;
    };
    Strand strands[2];
    /// Whether the span summary says a read can disagree.
    bool prefiltered = false;
    bool tested = false;
    bool over_base_cap = false;
    bool over_alignment_cap = false;
    /// Whether the site has two strands and all its reads were taken to be on strand 0.
    bool pooled = false;

    /// Return whether any strand has a consensus that is to become an allele.
    bool proposes() const { return strands[0].allele || strands[1].allele; }

    /// Return whether any strand is a candidate.
    bool candidate() const { return strands[0].candidate || strands[1].candidate; }
};

bool ReadConsensusStep::is_consensus_allele(const EditSnarlDecompositionOverlay& sites,
                                            const StagedSite& site, int allele) {
    return allele >= 0 && (size_t)allele < site.alleles.walks.size() && !sites.is_edit(site.net)
           && !site.alleles.spelled_by_walk((size_t)allele);
}

bool ReadConsensusStep::disagrees(const string& read, const string& called) const {
    const size_t longer = std::max(read.size(), called.size());
    if (longer - std::min(read.size(), called.size()) >= options.min_diff) {
        return true;
    }
    if (longer < options.min_diff) {
        // The two cannot differ by that many bases.
        return false;
    }
    const optional<int64_t> differing = edit_distance(read, called);
    return differing && (size_t)*differing >= options.min_diff
           && (double)*differing >= options.min_diff_share * (double)longer;
}

/// Return the bases `sequence` has between its first `before` and its last `after`, or nothing
/// if it has fewer bases than those.
static optional<string> between(const string& sequence, size_t before, size_t after) {
    if (sequence.size() < before + after) {
        return std::nullopt;
    }
    return sequence.substr(before, sequence.size() - before - after);
}

ReadConsensusStep::Finding ReadConsensusStep::examine(const Parts& parts,
                                                      const AlleleSplicer& splicer,
                                                      const StagedSite& site) const {
    Finding found;
    const HandleGraph& graph = *parts.graph;
    const AlleleReadLikelihoods::SpanSummary& spans = site.score->spans;
    const bool enough_span = spans.reads >= options.min_spanning;
    if (options.prefilter && !enough_span) {
        return found;
    }
    const handle_t start = parts.sites->get_site_start(site.net);
    const handle_t end = parts.sites->get_site_end(site.net);
    if (graph.get_id(start) == graph.get_id(end)) {
        return found;
    }
    const size_t before = graph.get_length(start);
    const size_t after = graph.get_length(end);

    // What each strand is called as between the bounds.
    optional<string> called[2];
    bool long_called = false;
    for (int strand = 0; strand < 2; ++strand) {
        if (const optional<string> spelled = splicer.strand_sequence(site, strand)) {
            called[strand] = between(*spelled, before, after);
        }
        long_called = long_called || (called[strand] && called[strand]->size() >= options.min_diff);
    }
    if (!called[0] && !called[1]) {
        return found;
    }
    found.prefiltered = enough_span && (long_called || spans.rank_bases >= options.min_diff);
    if (options.prefilter && !found.prefiltered) {
        return found;
    }
    for (const optional<string>& sequence : called) {
        if (sequence && sequence->size() > parts.max_bases) {
            found.over_base_cap = true;
            return found;
        }
    }

    const MaterializedSite view{start, before, end, 0, site.alleles};
    const vector<AlleleLikelihoodCalculator::SpanningRead> reads =
        parts.genotyper->calculator().spanning_reads(view);
    found.tested = true;

    // Put each read on a strand.
    const bool one_strand = !called[0] || !called[1];
    vector<int> strand_of(reads.size(), one_strand ? (called[0] ? 0 : 1) : -1);
    if (!one_strand) {
        unordered_map<uint64_t, double> own;
        const bool have_own = parts.read_strands->site_own_strand_log_odds(site.record_key, own);
        size_t told = 0;
        for (size_t r = 0; r < reads.size(); ++r) {
            const double log_odds = parts.read_strands->read_strand_log_odds(
                site.record_key, reads[r].name, have_own ? &own : nullptr);
            if (std::isfinite(log_odds) && log_odds != 0.0) {
                strand_of[r] = log_odds > 0.0 ? 0 : 1;
                ++told;
            }
        }
        if (2 * told < reads.size()) {
            found.pooled = true;
            std::ranges::fill(strand_of, 0);
        }
    }

    for (int strand = 0; strand < 2; ++strand) {
        if (!called[strand]) {
            continue;
        }
        Finding::Strand& on = found.strands[strand];
        vector<string> bases;
        for (size_t r = 0; r < reads.size(); ++r) {
            // A read with more bases than the scorer aligns does not span the site for it.
            if (strand_of[r] == strand && reads[r].length <= parts.max_bases) {
                bases.push_back(reads[r].bases);
                on.disagreeing += disagrees(reads[r].bases, *called[strand]) ? 1 : 0;
            }
        }
        on.reads = bases.size();
        on.candidate = on.reads >= options.min_spanning && on.disagreeing >= options.min_reads
                       && (double)on.disagreeing >= options.min_share * (double)on.reads;
        if (!on.candidate) {
            continue;
        }
        ReadConsensusParams building;
        building.max_reads = options.max_reads;
        optional<string> consensus = read_consensus(bases, *called[strand], building);
        const optional<int64_t> change =
            consensus ? edit_distance(*consensus, *called[strand]) : std::nullopt;
        if (!change) {
            continue;
        }
        on.built = true;
        if ((size_t)*change <= options.min_change) {
            on.unchanged = true;
        } else {
            on.allele = std::move(consensus);
        }
    }
    if (!found.proposes()) {
        return found;
    }

    // The scorer aligns the spanning reads to every allele, within its limits.
    size_t new_alleles = 0;
    for (size_t a = 0; a < site.alleles.walks.size(); ++a) {
        const optional<string> spelled = between(site.alleles.sequence(graph, a), before, after);
        if (!spelled || spelled->size() > parts.max_bases) {
            found.over_base_cap = true;
        }
        for (Finding::Strand& on : found.strands) {
            if (spelled && on.allele && *on.allele == *spelled) {
                on.allele.reset();
                on.already_allele = true;
            }
        }
    }
    for (const Finding::Strand& on : found.strands) {
        new_alleles += on.allele ? 1 : 0;
    }
    if (found.strands[0].allele && found.strands[1].allele
        && *found.strands[0].allele == *found.strands[1].allele) {
        --new_alleles;
    }
    if (new_alleles > 0 && !found.over_base_cap
        && reads.size() * (site.alleles.walks.size() + new_alleles) > parts.max_alignments) {
        found.over_alignment_cap = true;
    }
    if (found.over_base_cap || found.over_alignment_cap) {
        for (Finding::Strand& on : found.strands) {
            on.allele.reset();
        }
    }
    return found;
}

size_t ReadConsensusStep::run(const Parts& parts) {
    counted = ReadConsensusCounts();
    changed.clear();
    if (!options.enabled || parts.linker == nullptr || !parts.linker->enabled()) {
        return 0;
    }
    // The strands are those of the phase as it is now. The splicer reads the render phases.
    parts.phase_table->freeze_for_render(true);
    AlleleSplicer splicer(*parts.graph, *parts.sites, *parts.phase_table);
    splicer.index_before_hand_off(*parts.staged_sites);

    vector<StagedSite*> snarls;
    parts.staged_sites->for_each([&](StagedSite& site) {
        if (!site.dropped && site.score != nullptr && !site.genotype.empty()
            && !parts.sites->is_edit(site.net) && site.alleles.walks.size() > 0
            && parts.phase_table->rendered(site.record_key) != nullptr) {
            snarls.push_back(&site);
        }
    });
    counted.sites = snarls.size();

    if (const char* sites_name = getenv("VG_READ_CONSENSUS_SITES"); sites_name != nullptr) {
        // A table of every staged snarl on the reference, for laying other intervals over the
        // sites: where its bounds are, and how many bases its longest allele has between them.
        ofstream table(sites_name);
        table << "id\tcontig\tstart\tend\tlevel\tstart_node\tend_node\tmax_interior\tphased"
              << "\tdropped\n";
        parts.staged_sites->for_each([&](StagedSite& site) {
            const CandidateAlleles& alleles = site.alleles;
            if (parts.sites->is_edit(site.net) || site.reference.no_reference
                || alleles.ref_allele < 0 || (size_t)alleles.ref_allele >= alleles.walks.size()) {
                return;
            }
            const size_t bounds = parts.graph->get_length(parts.sites->get_site_start(site.net))
                                  + parts.graph->get_length(parts.sites->get_site_end(site.net));
            size_t longest = 0;
            for (size_t a = 0; a < alleles.walks.size(); ++a) {
                const size_t bases = alleles.sequence(*parts.graph, a).size();
                longest = std::max(longest, bases > bounds ? bases - bounds : 0);
            }
            const SiteLocus locus = site.reference.locus(*parts.graph, *parts.sites, site.net);
            table << site.id << '\t' << locus.contig << '\t' << locus.position << '\t'
                  << locus.position
                         + alleles.sequence(*parts.graph, (size_t)alleles.ref_allele).size()
                  << '\t' << site.placement.level << '\t'
                  << parts.graph->get_id(parts.sites->get_site_start(site.net)) << '\t'
                  << parts.graph->get_id(parts.sites->get_site_end(site.net)) << '\t' << longest
                  << '\t' << (parts.phase_table->rendered(site.record_key) != nullptr) << '\t'
                  << site.dropped << '\n';
        });
    }

    AlleleLikelihoodCalculator& calculator = parts.genotyper->calculator();
    calculator.set_counted_fetches(false);
    vector<Finding> findings(snarls.size());
#pragma omp parallel for schedule(dynamic, 64)
    for (size_t i = 0; i < snarls.size(); ++i) {
        findings[i] = examine(parts, splicer, *snarls[i]);
    }

    // Give each difference to the innermost site that holds it: decide the deepest sites first,
    // and leave a strand alone where a site under it gained an allele for the strand.
    vector<size_t> proposing;
    for (size_t i = 0; i < snarls.size(); ++i) {
        const Finding& found = findings[i];
        counted.prefiltered += found.prefiltered ? 1 : 0;
        counted.tested += found.tested ? 1 : 0;
        counted.over_base_cap += found.over_base_cap ? 1 : 0;
        counted.over_alignment_cap += found.over_alignment_cap ? 1 : 0;
        counted.candidate_sites += found.candidate() ? 1 : 0;
        counted.candidates_outside_prefilter += found.candidate() && !found.prefiltered ? 1 : 0;
        for (const Finding::Strand& on : found.strands) {
            counted.candidate_strands += on.candidate ? 1 : 0;
            counted.built += on.built ? 1 : 0;
            counted.no_consensus += on.candidate && !on.built ? 1 : 0;
            counted.unchanged += on.unchanged ? 1 : 0;
            counted.already_allele += on.already_allele ? 1 : 0;
        }
        if (found.proposes()) {
            proposing.push_back(i);
        }
    }
    std::ranges::stable_sort(proposing, [&](size_t a, size_t b) {
        return snarls[a]->placement.level > snarls[b]->placement.level;
    });
    // The strands each site gains an allele for, as a mask, by record key.
    unordered_map<size_t, int> gained;
    vector<size_t> gaining;
    for (const size_t i : proposing) {
        Finding& found = findings[i];
        int under = 0;
        splicer.for_each_under(snarls[i]->record_key, [&](const StagedSite& nested) {
            if (const auto mask = gained.find(nested.record_key); mask != gained.end()) {
                under |= mask->second;
            }
        });
        int mask = 0;
        for (int strand = 0; strand < 2; ++strand) {
            Finding::Strand& on = found.strands[strand];
            if (!on.allele) {
                continue;
            }
            // Pooled reads are those of both strands.
            if ((found.pooled ? under : under & (1 << strand)) != 0) {
                on.allele.reset();
                on.covered = true;
                ++counted.covered_by_nested;
            } else {
                mask |= found.pooled ? 3 : 1 << strand;
            }
        }
        if (mask != 0) {
            gained[snarls[i]->record_key] = mask;
            gaining.push_back(i);
        }
    }
    std::ranges::sort(gaining);

    // Add the alleles and genotype each site again. The sites differ, so they can be done at
    // once.
    vector<char> genotyped(gaining.size(), 0);
    vector<size_t> added(gaining.size(), 0);
#pragma omp parallel for schedule(dynamic, 1)
    for (size_t g = 0; g < gaining.size(); ++g) {
        StagedSite& site = *snarls[gaining[g]];
        const Finding& found = findings[gaining[g]];
        const CandidateAlleles kept = site.alleles;
        const size_t first_added = site.alleles.walks.size();
        for (int strand = 0; strand < 2; ++strand) {
            const optional<string>& allele = found.strands[strand].allele;
            if (allele && !(strand == 1 && found.strands[0].allele
                            && *found.strands[0].allele == *allele)) {
                added[g] += add_sequence_allele(*parts.graph, site.alleles, *allele) ? 1 : 0;
            }
        }
        if (added[g] > 0 && parts.tree_genotyper->regenotype(site, first_added)) {
            genotyped[g] = 1;
        } else {
            site.alleles = kept;
        }
    }
    calculator.set_counted_fetches(true);

    // File each site again, one at a time, in the order of the sites.
    LinkageCollector* model = parts.linker->collector();
    for (size_t g = 0; g < gaining.size(); ++g) {
        StagedSite& site = *snarls[gaining[g]];
        if (!genotyped[g]) {
            counted.regenotype_failed += added[g] > 0 ? 1 : 0;
            continue;
        }
        counted.added += added[g];
        ++counted.sites_with_added;
        changed.push_back(site.record_key);
        if (model->has_entry(site.record_key)) {
            // Retract before filing, because the model finds the first live entry for a key.
            model->retract(site.record_key);
            parts.linker->file(site, site.placement.level > 0 ? site.ploidy == 1
                                                              : site.placement.one_copy);
        }
    }

    if (const char* log_name = getenv("VG_READ_CONSENSUS_LOG"); log_name != nullptr) {
        // A table of every tested site, for comparing runs: where it is, and what was found on
        // each strand.
        ofstream log(log_name);
        log << "id\tcontig\tposition\tlevel\tprefiltered\tpooled\tcap\tstrand\treads\tdisagreeing"
            << "\tcandidate\toutcome\n";
        for (size_t i = 0; i < snarls.size(); ++i) {
            const Finding& found = findings[i];
            if (!found.tested) {
                continue;
            }
            const StagedSite& site = *snarls[i];
            const SiteLocus locus = site.reference.locus(*parts.graph, *parts.sites, site.net);
            const bool gains = std::ranges::find(changed, site.record_key) != changed.end();
            for (int strand = 0; strand < 2; ++strand) {
                const Finding::Strand& on = found.strands[strand];
                const char* outcome = on.allele ? (gains ? "added" : "not_genotyped")
                                      : on.covered ? "covered"
                                      : on.already_allele ? "already_allele"
                                      : on.unchanged ? "unchanged"
                                      : on.built ? "capped"
                                      : on.candidate ? "no_consensus" : "-";
                log << site.id << '\t' << locus.contig << '\t' << locus.position << '\t'
                    << site.placement.level << '\t' << found.prefiltered << '\t' << found.pooled
                    << '\t' << (found.over_base_cap ? "bases" : found.over_alignment_cap
                                ? "alignments" : "-")
                    << '\t' << strand << '\t' << on.reads << '\t' << on.disagreeing << '\t'
                    << on.candidate << '\t' << outcome << '\n';
            }
        }
    }

    // The phase is frozen again before the render, once it is final.
    parts.phase_table->freeze_for_render(false);
    return changed.size();
}

void ReadConsensusStep::count_called(const Parts& parts) {
    counted.sites_called = 0;
    if (changed.empty()) {
        return;
    }
    const unordered_map<size_t, StagedSite*> by_key = parts.staged_sites->by_key();
    for (const size_t key : changed) {
        const auto found = by_key.find(key);
        if (found == by_key.end() || found->second->dropped) {
            continue;
        }
        const StagedSite& site = *found->second;
        const vector<int> genotype = parts.linker->chosen_genotype(site);
        counted.sites_called += std::ranges::any_of(genotype, [&](int allele) {
            return is_consensus_allele(*parts.sites, site, allele);
        }) ? 1 : 0;
    }
}

void ReadConsensusStep::report(ostream& out) const {
    out << "[vg call] read consensus: " << counted.sites << " phased snarls, "
        << counted.prefiltered << " passed the prefilter, " << counted.tested << " tested, "
        << counted.candidate_strands << " candidate strands at " << counted.candidate_sites
        << " sites, " << counted.built << " consensuses built, " << counted.added
        << " alleles added at " << counted.sites_with_added << " sites, "
        << counted.sites_called << " sites called with a consensus allele; skipped: "
        << counted.over_base_cap << " sites over the base cap, " << counted.over_alignment_cap
        << " over the alignment cap, " << counted.covered_by_nested
        << " strands covered by a nested site, " << counted.no_consensus
        << " with no consensus, " << counted.unchanged << " consensuses within "
        << options.min_change << " bases of the call, " << counted.already_allele
        << " already an allele, " << counted.regenotype_failed << " sites not re-genotyped";
    if (!options.prefilter) {
        out << "; prefilter off: " << counted.candidates_outside_prefilter << " of "
            << counted.candidate_sites << " candidate sites fail the prefilter";
    }
    out << endl;
}

}
