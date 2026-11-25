#include <algorithm>
#include <cmath>
#include <limits>
#include <memory>
#include <numeric>
#include <unordered_set>

#include "read_consensus_step.hpp"
#include "site_values.hpp"

namespace vg::multipass {

using namespace std;

/// Return the bases `sequence` has between its first `before` and its last `after`, or nothing
/// if it has fewer bases than those.
static optional<string> between(const string& sequence, size_t before, size_t after) {
    if (sequence.size() < before + after) {
        return std::nullopt;
    }
    return sequence.substr(before, sequence.size() - before - after);
}

/// Add `consensus`, the bases a consensus has between the bounds of the site of `alleles`, as
/// `walk` with changes to its steps where `walk` is given, and as a sequence allele otherwise.
/// Count it in `entered` as threaded where it is the walk.
static optional<size_t> add_consensus(const HandleGraph& graph, CandidateAlleles& alleles,
                                      const optional<Traversal>& walk, const string& consensus,
                                      ReadConsensusCounts& entered) {
    if (walk) {
        if (optional<size_t> allele = add_walk_allele(graph, alleles, *walk, consensus,
                                                      ReadConsensusParams().aligner)) {
            ++entered.threaded;
            return allele;
        }
    }
    return add_sequence_allele(graph, alleles, consensus);
}

struct ReadConsensusStep::Run {
    /// The members, as indices of the phased snarls, in the order of their chain. A snarl on
    /// its own is a run of one member.
    vector<size_t> members;
    /// The nodes from each member's end bound to the next member's start bound, both included.
    /// It is the one node where the two share a bound.
    vector<vector<handle_t>> junctions;
};

struct ReadConsensusStep::RunFinding {
    /// What the step found on one strand.
    struct Strand {
        /// The strand's spanning reads.
        size_t reads = 0;
        /// Whether the consensus of the strand's reads differs from what the strand is called
        /// as.
        bool candidate = false;
        /// The bases the consensus puts between the run's two outer bounds, where the run is
        /// treated for it: it is to become an allele of the snarl, or of the one site that a
        /// run of several becomes.
        optional<string> consensus;
    };
    Strand strands[2];
    /// The bases each strand is called as between the two outer bounds of a run that is
    /// treated, where the strand crosses every member: at a run of one, where no allele of the
    /// snarl spells them.
    optional<string> called[2];
    /// The record keys of the edit sites inside the nodes between the members.
    vector<size_t> edits_between;
    /// Whether the run's reads were fetched and compared with the called sequences.
    bool fetched = false;
    /// Whether a called sequence or a consensus has more bases between the bounds than the
    /// scorer aligns of a read, or the reads and alleles need more alignments than it makes.
    bool over_cap = false;

    /// Return whether any strand has a consensus that the run is treated for.
    bool treats() const { return strands[0].consensus || strands[1].consensus; }
};

struct ReadConsensusStep::Alone {
    /// Whether the site has a called sequence and enough spanning reads to be in a run.
    bool measured = false;
    /// Whether the span summary says the site's reads can gain `min_gain` nats by an allele.
    bool can_disagree = false;
    /// The most bases a called sequence has between the bounds.
    size_t longest = 0;
    /// One read's say for each spanning read that no allele explains.
    double signal = 0.0;
    /// The longest gap that enough reads to gain `min_gain` nats have inside the node of the
    /// start bound or in the nodes before it, and inside that of the end bound or in the nodes
    /// after it (`SpanSummary::start_gap`, `end_gap`).
    size_t start_gap = 0;
    size_t end_gap = 0;
    /// Whether the site's reads were fetched and compared with the called sequences.
    bool tested = false;
    /// Whether the site is over a limit of the scorer (`RunFinding::over_cap`).
    bool over_cap = false;
    /// The number of the site's strands that are candidates.
    size_t candidates = 0;
    /// Whether a strand has a consensus that is to become an allele of the site.
    bool treated = false;
};

ReadConsensusStep::Alone ReadConsensusStep::measure(const Parts& parts,
                                                    const AlleleSplicer& splicer,
                                                    const StagedSite& site) const {
    Alone alone;
    const HandleGraph& graph = *parts.graph;
    const AlleleReadLikelihoods::SpanSummary& spans = site.score->spans;
    if (spans.reads < options.reads_needed(clamp)) {
        // Too few reads span the site to gain `min_gain` nats from any allele.
        return alone;
    }
    const handle_t start = parts.sites->get_site_start(site.net);
    const handle_t end = parts.sites->get_site_end(site.net);
    if (graph.get_id(start) == graph.get_id(end)) {
        return alone;
    }
    size_t shortest = std::numeric_limits<size_t>::max();
    for (int strand = 0; strand < 2; ++strand) {
        optional<string> called;
        if (const optional<string> spelled = splicer.strand_sequence(site, strand)) {
            called = between(*spelled, graph.get_length(start), graph.get_length(end));
        }
        if (called) {
            alone.longest = std::max(alone.longest, called->size());
            shortest = std::min(shortest, called->size());
        }
    }
    if (shortest == std::numeric_limits<size_t>::max()) {
        // No strand has a called sequence.
        return alone;
    }
    alone.measured = true;
    alone.start_gap = spans.start_gap;
    alone.end_gap = spans.end_gap;
    // Each read that no allele explains can gain a new allele one read's say.
    alone.signal = clamp * (double)spans.unexplained;
    alone.can_disagree = alone.signal >= options.min_gain;
    return alone;
}

void ReadConsensusStep::find_runs(const Parts& parts, const AlleleSplicer& splicer,
                                  const vector<StagedSite*>& snarls,
                                  const vector<Alone>& alone, vector<Run>& runs,
                                  vector<RunFinding>& run_findings) const {
    const HandleGraph& graph = *parts.graph;
    const size_t none = std::numeric_limits<size_t>::max();
    // A site can be in a run where enough reads span it and it has a called sequence of no
    // more bases than the scorer aligns of a read.
    auto can_join = [&](size_t i) {
        return alone[i].measured && alone[i].longest <= max_bases;
    };
    // Two sites can be in one run where one phase orders the strands of both.
    auto same_phase = [&](size_t a, size_t b) {
        const LinkageCollector::PhaseCall* first = parts.phase_table->rendered(
            snarls[a]->record_key);
        const LinkageCollector::PhaseCall* second = parts.phase_table->rendered(
            snarls[b]->record_key);
        return first->contig == second->contig && first->phase_set == second->phase_set
               && first->ploidy == second->ploidy
               && first->nested_strand == second->nested_strand;
    };
    unordered_map<handle_t, size_t> by_start;
    for (size_t i = 0; i < snarls.size(); ++i) {
        if (can_join(i)) {
            by_start[parts.sites->get_site_start(snarls[i]->net)] = i;
        }
    }
    // The signal of each site: the unexplained nats of its own span summary, or those of the
    // gap its reads have in a node between it and a neighbour, if that is more.
    vector<double> signal(snarls.size());
    const AlleleLikelihoodCalculator& calculator = parts.genotyper->calculator();
    // Return whether a gap of `bases` bases leaves a read unexplained: it costs the read
    // `min_gain` nats.
    auto unexplained_gap = [&](size_t bases) {
        const string gap(bases, 'A');
        const optional<double> spelled = calculator.sequence_nats(gap, gap);
        const optional<double> lacking = calculator.sequence_nats(gap, "");
        return bases > 0 && spelled && lacking && *spelled - *lacking >= options.min_gain;
    };
    // The nats of a gap that enough reads to gain `min_gain` have.
    const double gap_signal = clamp * (double)options.reads_needed(clamp);
    for (size_t i = 0; i < snarls.size(); ++i) {
        signal[i] = alone[i].signal;
    }
    // The site after each site in its chain, where the two are neighbours, and the nodes from
    // the one's end bound to the other's start bound.
    vector<size_t> next(snarls.size(), none);
    vector<char> follows(snarls.size(), 0);
    vector<vector<handle_t>> nodes_after(snarls.size());
    vector<size_t> bases_after(snarls.size(), 0);
    for (size_t i = 0; i < snarls.size(); ++i) {
        if (!can_join(i)) {
            continue;
        }
        vector<handle_t> nodes;
        size_t bases = 0;
        size_t found = none;
        const bool reached = walk_to_next_site(
            graph, parts.sites->get_site_end(snarls[i]->net), options.run_max_junction,
            [&](const handle_t& at) {
                const auto starting = by_start.find(at);
                if (starting == by_start.end() || starting->second == i) {
                    return false;
                }
                found = starting->second;
                return true;
            },
            nodes, bases);
        if (reached && same_phase(i, found) && !follows[found]) {
            next[i] = found;
            follows[found] = 1;
            nodes_after[i] = nodes;
            bases_after[i] = bases;
            // A gap that the reads have on the side of a bound that faces the other site is a
            // signal of both.
            if (unexplained_gap(std::max(alone[i].end_gap, alone[found].start_gap))) {
                signal[i] = std::max(signal[i], gap_signal);
                signal[found] = std::max(signal[found], gap_signal);
            }
        }
    }
    // A seed is a site whose reads have `min_gain` nats that nothing called explains.
    auto is_seed = [&](size_t i) { return signal[i] >= options.min_gain; };

    // A run begins and ends at a seed, and holds every site between two of its seeds. Two
    // seeds that follow each other are in one run where at most `run_max_junction` bases lie
    // between them, in the nodes of the junctions and in the sites that are no seeds, and the
    // run has no more bases between its outer bounds than the scorer aligns of a read. The
    // bases of a site are those of its longest called sequence. A run is tested where it has
    // two or more members.
    vector<Run> tested;
    for (size_t i = 0; i < snarls.size(); ++i) {
        if (!can_join(i) || follows[i]) {
            continue;
        }
        Run run;
        size_t bases = 0;
        // The sites after the run's last seed, which are no seeds, and the bases after that
        // seed: theirs, and those of the junctions.
        vector<size_t> skipped;
        size_t between = 0;
        auto close = [&]() {
            if (run.members.size() > 1) {
                tested.push_back(std::move(run));
            }
            run = Run();
            bases = 0;
            skipped.clear();
            between = 0;
        };
        size_t before = none;
        for (size_t at = i; at != none; before = at, at = next[at]) {
            if (!run.members.empty()) {
                between += bases_after[before];
            }
            if (!is_seed(at)) {
                if (!run.members.empty()) {
                    skipped.push_back(at);
                    between += alone[at].longest;
                }
                continue;
            }
            if (!run.members.empty()
                && (between > options.run_max_junction
                    || bases + between + alone[at].longest > max_bases)) {
                close();
            }
            if (!run.members.empty()) {
                skipped.push_back(at);
                for (const size_t member : skipped) {
                    run.junctions.push_back(nodes_after[run.members.back()]);
                    run.members.push_back(member);
                }
                bases += between;
                skipped.clear();
                between = 0;
            } else {
                run.members.push_back(at);
            }
            bases += alone[at].longest;
        }
        close();
    }

    // The runs share no site, so they are tested at once.
    vector<RunFinding> found_at(tested.size());
#pragma omp parallel for schedule(dynamic, 1)
    for (size_t r = 0; r < tested.size(); ++r) {
        found_at[r] = examine_run(parts, splicer, snarls, tested[r]);
    }
    for (size_t r = 0; r < tested.size(); ++r) {
        if (found_at[r].treats()) {
            runs.push_back(std::move(tested[r]));
            run_findings.push_back(std::move(found_at[r]));
        }
    }
}

ReadConsensusStep::RunFinding ReadConsensusStep::examine_run(const Parts& parts,
                                                             const AlleleSplicer& splicer,
                                                             const vector<StagedSite*>& snarls,
                                                             const Run& run) const {
    RunFinding found;
    const HandleGraph& graph = *parts.graph;
    const size_t count = run.members.size();

    // What each strand is called as across the run: its called sequence at each member, with
    // the junctions between them.
    string called[2];
    bool crosses[2] = {true, true};
    for (size_t m = 0; m < count; ++m) {
        const StagedSite& site = *snarls[run.members[m]];
        const size_t before = graph.get_length(parts.sites->get_site_start(site.net));
        const size_t after = graph.get_length(parts.sites->get_site_end(site.net));
        string junction;
        if (m > 0) {
            for (const handle_t& node : run.junctions[m - 1]) {
                junction += graph.get_sequence(node);
            }
        }
        for (int strand = 0; strand < 2; ++strand) {
            optional<string> interior;
            if (const optional<string> spelled = splicer.strand_sequence(site, strand)) {
                interior = between(*spelled, before, after);
            }
            if (!interior) {
                crosses[strand] = false;
            }
            if (crosses[strand]) {
                called[strand] += junction + *interior;
            }
        }
    }
    if (!crosses[0] && !crosses[1]) {
        return found;
    }
    for (int strand = 0; strand < 2; ++strand) {
        if (crosses[strand] && called[strand].size() > max_bases) {
            found.over_cap = true;
            return found;
        }
    }

    // The reads that span the run are those that span a snarl with the run's two outer bounds.
    // Its walks say which way the run reads each node: they are the members' walks, joined.
    CandidateAlleles across;
    size_t most_walks = 0;
    for (const size_t i : run.members) {
        most_walks = std::max(most_walks, snarls[i]->alleles.walks.size());
    }
    for (size_t w = 0; w < most_walks; ++w) {
        Traversal walk;
        for (size_t m = 0; m < count; ++m) {
            const vector<Traversal>& walks = snarls[run.members[m]]->alleles.walks;
            const Traversal& part = walks[std::min(w, walks.size() - 1)];
            if (m == 0) {
                walk = part;
                continue;
            }
            walk.insert(walk.end(), run.junctions[m - 1].begin() + 1, run.junctions[m - 1].end());
            walk.insert(walk.end(), part.begin() + 1, part.end());
        }
        across.walks.push_back(std::move(walk));
    }
    const handle_t start = parts.sites->get_site_start(snarls[run.members.front()]->net);
    const handle_t end = parts.sites->get_site_end(snarls[run.members.back()]->net);
    const MaterializedSite view{start, graph.get_length(start), end, 0, across};
    const vector<AlleleLikelihoodCalculator::SpanningRead> reads =
        parts.genotyper->calculator().spanning_reads(view);
    found.fetched = true;
    const AlleleLikelihoodCalculator& calculator = parts.genotyper->calculator();

    // Put each read on the strand that its strand log-odds at the members, summed, favour.
    const bool one_strand = !crosses[0] || !crosses[1];
    vector<int> strand_of(reads.size(), one_strand ? (crosses[0] ? 0 : 1) : -1);
    if (!one_strand) {
        vector<double> favour(reads.size(), 0.0);
        for (const size_t i : run.members) {
            const size_t key = snarls[i]->record_key;
            unordered_map<uint64_t, double> own;
            const bool have_own = parts.read_strands->site_own_strand_log_odds(key, own);
            for (size_t r = 0; r < reads.size(); ++r) {
                const double log_odds = parts.read_strands->read_strand_log_odds(
                    key, reads[r].name, have_own ? &own : nullptr);
                favour[r] += std::isfinite(log_odds) ? log_odds : 0.0;
            }
        }
        size_t told = 0;
        for (size_t r = 0; r < reads.size(); ++r) {
            if (favour[r] != 0.0) {
                strand_of[r] = favour[r] > 0.0 ? 0 : 1;
                ++told;
            }
        }
        if (2 * told < reads.size()) {
            std::ranges::fill(strand_of, 0);
        }
    }

    for (int strand = 0; strand < 2; ++strand) {
        if (!crosses[strand]) {
            continue;
        }
        RunFinding::Strand& on = found.strands[strand];
        // The strand's reads. Each gains one read's say at most.
        vector<string> bases;
        for (size_t r = 0; r < reads.size(); ++r) {
            if (strand_of[r] != strand || reads[r].length > max_bases) {
                continue;
            }
            bases.push_back(reads[r].bases);
        }
        on.reads = bases.size();
        if (clamp * (double)on.reads < options.min_gain) {
            continue;
        }
        optional<string> consensus = read_consensus(bases, called[strand]);
        if (!consensus || *consensus == called[strand]) {
            continue;
        }
        on.candidate = true;
        // The consensus becomes an allele where the reads gain `min_gain` nats by it.
        double gain = 0.0;
        for (size_t r = 0; r < bases.size(); ++r) {
            const optional<double> against = calculator.sequence_nats(bases[r], called[strand]);
            const optional<double> with = calculator.sequence_nats(bases[r], *consensus);
            if (with && against) {
                gain += std::clamp(*with - *against, -clamp, clamp);
            }
        }
        if (gain >= options.min_gain) {
            on.consensus = std::move(consensus);
        }
    }
    auto leave_alone = [&]() {
        for (RunFinding::Strand& on : found.strands) {
            on.consensus.reset();
        }
    };
    if (!found.treats()) {
        return found;
    }

    // The alleles the scorer aligns the spanning reads to, besides the consensuses.
    size_t alleles = 0;
    if (count == 1) {
        // A snarl on its own gains a consensus that none of its alleles spells. With it, it
        // gains what each strand is called as, where the sites nested in the snarl make that
        // a sequence none of its alleles spells, so that a strand can keep what it has.
        const StagedSite& site = *snarls[run.members.front()];
        const size_t walks = site.alleles.walks.size();
        for (int strand = 0; strand < 2; ++strand) {
            if (crosses[strand] && !(strand == 1 && crosses[0] && called[0] == called[1])) {
                found.called[strand] = called[strand];
            }
        }
        for (size_t a = 0; a < walks; ++a) {
            const optional<string> spelled = between(site.alleles.sequence(graph, a),
                                                     graph.get_length(start),
                                                     graph.get_length(end));
            for (int strand = 0; strand < 2; ++strand) {
                RunFinding::Strand& on = found.strands[strand];
                if (on.consensus && spelled && *on.consensus == *spelled) {
                    on.consensus.reset();
                }
                if (found.called[strand] && spelled && *found.called[strand] == *spelled) {
                    found.called[strand].reset();
                }
            }
        }
        alleles = walks + (found.called[0] ? 1 : 0) + (found.called[1] ? 1 : 0);
    } else if (can_merge(parts, snarls, run, found)) {
        // A run of several becomes one site, whose alleles are the reference across the run,
        // what each strand is called as, and each consensus.
        alleles = 3;
        for (int strand = 0; strand < 2; ++strand) {
            if (crosses[strand]) {
                found.called[strand] = called[strand];
            }
        }
    } else {
        leave_alone();
        return found;
    }
    for (RunFinding::Strand& on : found.strands) {
        if (on.consensus && on.consensus->size() > max_bases) {
            on.consensus.reset();
            found.over_cap = true;
        }
    }
    const RunFinding::Strand* strands = found.strands;
    const size_t added = (strands[0].consensus ? 1 : 0)
                         + (strands[1].consensus && strands[1].consensus != strands[0].consensus
                                ? 1 : 0);
    if (added > 0 && reads.size() * (alleles + added) > max_alignments) {
        found.over_cap = true;
        leave_alone();
    }
    return found;
}

bool ReadConsensusStep::can_merge(const Parts& parts, const vector<StagedSite*>& snarls,
                                  const Run& run, RunFinding& found) const {
    const StagedSite& first = *snarls[run.members.front()];
    for (const size_t i : run.members) {
        const StagedSite& member = *snarls[i];
        // The one site is placed on the reference as its first member is, and has a reference
        // allele where every member has one.
        if (member.reference.path_name != first.reference.path_name
            || member.reference.offset != first.reference.offset
            || member.reference.no_reference != first.reference.no_reference
            || member.placement.level != first.placement.level
            || member.placement.parent_record_key != first.placement.parent_record_key
            || member.ploidy != first.ploidy
            || (member.alleles.ref_allele < 0 && !member.reference.no_reference)) {
            return false;
        }
    }
    // An edit site inside a node between two members is no member, and what is called at it is
    // not in the strands' called sequences. The run is left as it is where one has a change
    // called.
    for (const vector<handle_t>& junction : run.junctions) {
        for (const handle_t& node : junction) {
            const auto on_node = edit_sites_on.find(parts.graph->get_id(node));
            if (on_node == edit_sites_on.end()) {
                continue;
            }
            for (const StagedSite* edit : on_node->second) {
                vector<int> carried = edit->genotype;
                if (const LinkageCollector::PhaseCall* phase =
                        parts.phase_table->rendered(edit->record_key)) {
                    carried = {phase->trav_first};
                    if (phase->ploidy != 1) {
                        carried.push_back(phase->trav_second);
                    }
                }
                if (std::ranges::any_of(carried, [&](int allele) {
                        return allele != edit->alleles.ref_allele;
                    })) {
                    return false;
                }
                if (std::ranges::find(found.edits_between, edit->record_key)
                    == found.edits_between.end()) {
                    found.edits_between.push_back(edit->record_key);
                }
            }
        }
    }
    return true;
}

optional<StagedSite> ReadConsensusStep::merge_run(const Parts& parts,
                                                  const AlleleSplicer& splicer,
                                                  const vector<StagedSite*>& snarls,
                                                  const Run& run,
                                                  const RunFinding& found,
                                                  ReadConsensusCounts& entered) const {
    const HandleGraph& graph = *parts.graph;
    const StagedSite& first = *snarls[run.members.front()];
    const StagedSite& last = *snarls[run.members.back()];
    const handle_t start = parts.sites->get_site_start(first.net);
    const handle_t end = parts.sites->get_site_end(last.net);
    const size_t before = graph.get_length(start);
    const size_t after = graph.get_length(end);

    // Return the walk across the run that takes the walk `pick` names at each member, or
    // nothing if it names none at one of them.
    auto joined = [&](const function<optional<int>(const StagedSite&)>& pick) {
        optional<Traversal> walk(std::in_place);
        for (size_t m = 0; m < run.members.size(); ++m) {
            const StagedSite& member = *snarls[run.members[m]];
            const optional<int> allele = pick(member);
            if (!allele || *allele < 0 || (size_t)*allele >= member.alleles.walks.size()) {
                walk.reset();
                return walk;
            }
            const Traversal& part = member.alleles.walks[*allele];
            if (m > 0) {
                walk->insert(walk->end(), run.junctions[m - 1].begin() + 1,
                             run.junctions[m - 1].end());
            }
            walk->insert(walk->end(), part.begin() + (m > 0 ? 1 : 0), part.end());
        }
        return walk;
    };
    // Return the bases the nodes of `walk` have between the run's bounds.
    auto interior_of = [&](const Traversal& walk) {
        string spelled;
        for (const handle_t& node : walk) {
            spelled += graph.get_sequence(node);
        }
        return between(spelled, before, after);
    };

    // The first allele is the reference where every member has one, and each member's first
    // walk otherwise.
    CandidateAlleles across;
    vector<string> interiors;
    bool all_reference = true;
    const optional<Traversal> base = joined([&](const StagedSite& member) {
        all_reference = all_reference && member.alleles.ref_allele >= 0;
        return optional<int>(std::max(member.alleles.ref_allele, 0));
    });
    const optional<string> base_bases = base ? interior_of(*base) : std::nullopt;
    if (!base_bases) {
        return std::nullopt;
    }
    across.walks.push_back(*base);
    across.ref_allele = all_reference ? 0 : -1;
    interiors.push_back(*base_bases);
    // Add the allele with `bases` between the bounds, unless there is one: as `walk` where
    // that spells them, and as a sequence allele otherwise.
    auto add = [&](const string& bases, const optional<Traversal>& walk) {
        if (std::ranges::find(interiors, bases) != interiors.end()) {
            return;
        }
        if (walk && interior_of(*walk) == bases) {
            across.walks.push_back(*walk);
        } else {
            add_sequence_allele(graph, across, bases);
        }
        interiors.push_back(bases);
    };
    // The walk each strand is called as across the run, where it spells what the strand is
    // called as.
    optional<Traversal> called_walk[2];
    for (int strand = 0; strand < 2; ++strand) {
        if (found.called[strand]) {
            called_walk[strand] = joined([&](const StagedSite& member) {
                return splicer.strand_traversal(member, strand);
            });
            add(*found.called[strand], called_walk[strand]);
            if (called_walk[strand] && interior_of(*called_walk[strand]) != found.called[strand]) {
                called_walk[strand].reset();
            }
        }
    }
    // The alleles so far spell what was called at the member sites, some of them as sequences.
    // Each consensus is the walk its strand is called as, with changes to its steps.
    const size_t first_consensus = across.walks.size();
    for (int strand = 0; strand < 2; ++strand) {
        const optional<string>& consensus = found.strands[strand].consensus;
        if (consensus && std::ranges::find(interiors, *consensus) == interiors.end()) {
            entered.unspelled += called_walk[strand] ? 0 : 1;
            add_consensus(graph, across, called_walk[strand], *consensus, entered);
            interiors.push_back(*consensus);
        }
    }
    if (across.walks.size() == first_consensus) {
        // Each consensus is what the other strand is called as, or the reference.
        return std::nullopt;
    }
    optional<StagedSite> site = parts.tree_genotyper->merged_site(first, last, std::move(across),
                                                                  first_consensus);
    if (site) {
        site->first_consensus_allele = first_consensus;
    }
    return site;
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
    // The edit sites inside each node, which a run is not made one site across where a change
    // is called at one.
    edit_sites_on.clear();
    parts.staged_sites->for_each([&](StagedSite& site) {
        if (!site.dropped && parts.sites->is_edit(site.net) && !site.alleles.walks.empty()) {
            for (const handle_t& node : site.alleles.walks.front()) {
                edit_sites_on[parts.graph->get_id(node)].push_back(&site);
            }
        }
    });

    AlleleLikelihoodCalculator& calculator = parts.genotyper->calculator();
    const AlleleLikelihoodCalculator::SequenceAlleleLimits limits =
        calculator.sequence_allele_limits();
    max_bases = limits.max_bases;
    max_alignments = limits.max_alignments;
    clamp = calculator.read_clamp();
    calculator.set_counted_fetches(false);
    // Test each site on its own, as a run of one, where its reads can gain by an allele.
    vector<Alone> alone(snarls.size());
    vector<unique_ptr<RunFinding>> found_alone(snarls.size());
#pragma omp parallel for schedule(dynamic, 64)
    for (size_t i = 0; i < snarls.size(); ++i) {
        alone[i] = measure(parts, splicer, *snarls[i]);
        if (!alone[i].can_disagree) {
            continue;
        }
        RunFinding found = examine_run(parts, splicer, snarls, Run{{i}, {}});
        alone[i].tested = found.fetched;
        alone[i].over_cap = found.over_cap;
        alone[i].candidates = (found.strands[0].candidate ? 1 : 0)
                              + (found.strands[1].candidate ? 1 : 0);
        alone[i].treated = found.treats();
        if (alone[i].treated) {
            found_alone[i] = std::make_unique<RunFinding>(std::move(found));
        }
    }

    // The runs that are treated: the runs of several sites with a consensus, and each
    // site with a consensus of its own that is in none of those.
    vector<Run> runs;
    vector<RunFinding> run_findings;
    find_runs(parts, splicer, snarls, alone, runs, run_findings);
    vector<char> in_run(snarls.size(), 0);
    for (const Run& run : runs) {
        for (const size_t i : run.members) {
            in_run[i] = 1;
        }
    }
    for (size_t i = 0; i < snarls.size(); ++i) {
        counted.tested += alone[i].tested ? 1 : 0;
        counted.over_cap += alone[i].over_cap ? 1 : 0;
        counted.candidate_strands += alone[i].candidates;
        if (found_alone[i] && !in_run[i]) {
            runs.push_back(Run{{i}, {}});
            run_findings.push_back(std::move(*found_alone[i]));
        }
    }
    found_alone.clear();

    // Every treated run is genotyped again, whatever is treated under it: the consensus allele
    // of a snarl crosses no site nested in the snarl, and the sites nested in a run that becomes
    // one site are removed with the run's members. The runs are taken in the order of their
    // first members.
    vector<size_t> tasks(runs.size());
    std::iota(tasks.begin(), tasks.end(), 0);
    std::ranges::sort(tasks, {}, [&](size_t r) { return runs[r].members.front(); });

    // Add the alleles and genotype each site again. The sites differ, so they can be done at
    // once.
    vector<char> genotyped(tasks.size(), 0);
    vector<size_t> added(tasks.size(), 0);
    // How each task entered its consensus alleles.
    vector<ReadConsensusCounts> entered(tasks.size());
    // The one site of each run that becomes one.
    vector<optional<StagedSite>> merged(tasks.size());
#pragma omp parallel for schedule(dynamic, 1)
    for (size_t t = 0; t < tasks.size(); ++t) {
        const Run& run = runs[tasks[t]];
        const RunFinding& found = run_findings[tasks[t]];
        if (run.members.size() > 1) {
            merged[t] = merge_run(parts, splicer, snarls, run, found, entered[t]);
            genotyped[t] = merged[t] ? 1 : 0;
            continue;
        }
        StagedSite& site = *snarls[run.members.front()];
        const CandidateAlleles kept = site.alleles;
        const size_t first_added = site.alleles.walks.size();
        for (const optional<string>& called : found.called) {
            if (called) {
                add_sequence_allele(*parts.graph, site.alleles, *called);
            }
        }
        const size_t first_consensus = site.alleles.walks.size();
        for (int strand = 0; strand < 2; ++strand) {
            const optional<string>& consensus = found.strands[strand].consensus;
            if (consensus && !(strand == 1 && consensus == found.strands[0].consensus)) {
                // The strand's called walk, where it spells what the strand is called as and
                // crosses no site nested in the snarl.
                optional<Traversal> walk(std::in_place);
                StepSequences steps;
                if (!site.is_leaf()) {
                    walk.reset();
                    ++entered[t].nested;
                } else if (!splicer.strand_walk(site, strand, *walk, steps) || !steps.empty()) {
                    walk.reset();
                    ++entered[t].unspelled;
                }
                added[t] += add_consensus(*parts.graph, site.alleles, walk, *consensus,
                                          entered[t]) ? 1 : 0;
            }
        }
        if (added[t] > 0 && parts.tree_genotyper->regenotype(site, first_added)) {
            genotyped[t] = 1;
            site.first_consensus_allele = first_consensus;
        } else {
            site.alleles = kept;
        }
    }
    calculator.set_counted_fetches(true);

    // File each site again, one at a time, in the order of the sites.
    LinkageCollector* model = parts.linker->collector();
    size_t filed = 0;
    // The record keys of the sites that a run made one site takes the place of.
    unordered_set<size_t> absorbed;
    for (size_t t = 0; t < tasks.size(); ++t) {
        const Run& run = runs[tasks[t]];
        counted.runs_treated += run.members.size() > 1 ? 1 : 0;
        StagedSite& site = *snarls[run.members.front()];
        if (!genotyped[t]) {
            continue;
        }
        ++filed;
        counted.threaded += entered[t].threaded;
        counted.unspelled += entered[t].unspelled;
        counted.nested += entered[t].nested;
        if (merged[t]) {
            // The one site takes the place of the run's first member. The other members, the
            // sites under any member and the edit sites between them are removed below.
            ++counted.runs_merged;
            for (const size_t i : run.members) {
                absorbed.insert(snarls[i]->record_key);
            }
            absorbed.insert(run_findings[tasks[t]].edits_between.begin(),
                            run_findings[tasks[t]].edits_between.end());
            model->retract(site.record_key);
            site = std::move(*merged[t]);
            changed.push_back(site.record_key);
            parts.linker->file(site, site.placement.level > 0 ? site.ploidy == 1
                                                              : site.placement.one_copy);
            continue;
        }
        if (added[t] > 0) {
            changed.push_back(site.record_key);
        }
        if (model->has_entry(site.record_key)) {
            // Retract before filing, because the model finds the first live entry for a key.
            model->retract(site.record_key);
            parts.linker->file(site, site.placement.level > 0 ? site.ploidy == 1
                                                              : site.placement.one_copy);
        }
    }

    // Remove what the merged sites take the place of, which can be a site that gained an allele
    // or was made of a run itself. This moves the staged sites, so nothing above is read after
    // it.
    edit_sites_on.clear();
    if (!absorbed.empty()) {
        unordered_set<size_t> removed;
        parts.staged_sites->remove_with_nested(absorbed, [&](const StagedSite& gone) {
            model->retract(gone.record_key);
            removed.insert(gone.record_key);
        });
        std::erase_if(changed, [&](size_t key) { return removed.contains(key); });
    }
    counted.sites_changed = changed.size();

    // The phase is frozen again before the render, once it is final.
    parts.phase_table->freeze_for_render(false);
    return filed;
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
            return site.is_consensus_allele(allele);
        }) ? 1 : 0;
    }
}

void ReadConsensusStep::report(ostream& out) const {
    out << "[vg call] read consensus: " << counted.sites << " phased snarls, " << counted.tested
        << " tested, " << counted.candidate_strands << " candidate strands, "
        << counted.runs_treated << " runs treated, " << counted.runs_merged
        << " of them made one site, " << counted.sites_changed
        << " sites gained an allele or were made of a run, " << counted.sites_called
        << " of them called with a consensus allele; consensus alleles: " << counted.threaded
        << " as a called walk with changes, " << counted.unspelled
        << " as bases where the called sequence is no walk, " << counted.nested
        << " as bases at a snarl with nested sites; left alone: " << counted.over_cap
        << " snarls over a limit of the scorer" << endl;
}

}
