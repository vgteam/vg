#include <algorithm>
#include <cmath>
#include <deque>
#include <limits>
#include <map>
#include <set>

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
    /// Whether the site has a called sequence and enough spanning reads to be in a run.
    bool measured = false;
    /// The most bases a called sequence has between the bounds.
    size_t longest = 0;
    /// The number of bases by which `min_reads` of the spanning reads have more between the
    /// bounds than the longest called sequence, or fewer than the shortest, whichever is more.
    size_t signal = 0;

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
    if (enough_span) {
        size_t shortest = std::numeric_limits<size_t>::max();
        for (const optional<string>& sequence : called) {
            if (sequence) {
                found.longest = std::max(found.longest, sequence->size());
                shortest = std::min(shortest, sequence->size());
            }
        }
        found.measured = true;
        found.signal = std::max(
            spans.rank_bases > found.longest ? spans.rank_bases - found.longest : 0,
            shortest > spans.low_bases ? shortest - spans.low_bases : 0);
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

struct ReadConsensusStep::Run {
    /// The members, as indices of the phased snarls, in the order of their chain.
    vector<size_t> members;
    /// The nodes from each member's end bound to the next member's start bound, both included.
    /// It is the one node where the two share a bound.
    vector<vector<handle_t>> junctions;
    /// The number of runs that were tested to grow this one, itself included.
    size_t tests = 0;
    /// Whether the run stopped growing on a side at the most bases the scorer aligns.
    bool capped = false;
};

namespace {

/// A sequence across a run, from after the first member's start bound to before the last
/// member's end bound, and where the bases between its members are in it.
struct RunSequence {
    string bases;
    /// The offset and the length in `bases` of each junction: the bases of the nodes from one
    /// member's end bound to the next member's start bound.
    vector<pair<size_t, size_t>> junctions;

    /// Return the bases each member has between its bounds.
    vector<string> pieces() const {
        vector<string> cut;
        size_t from = 0;
        for (const auto& [at, length] : junctions) {
            cut.push_back(bases.substr(from, at - from));
            from = at + length;
        }
        cut.push_back(bases.substr(from));
        return cut;
    }
};

/// Return where in the second string of the alignment `runs` each of `junctions`, stretches of
/// the first, lies, for each that the alignment matches whole, base for base.
vector<optional<size_t>> matched_junctions(const vector<AlignmentRun>& runs,
                                           const vector<pair<size_t, size_t>>& junctions) {
    vector<optional<size_t>> found(junctions.size());
    size_t i = 0;
    size_t j = 0;
    size_t next = 0;
    for (const AlignmentRun& run : runs) {
        if (run.op == AlignmentOp::Match) {
            // The junctions are in order, so those that start before this run are not in it.
            while (next < junctions.size() && junctions[next].first < i) {
                ++next;
            }
            for (size_t k = next; k < junctions.size() && junctions[k].first < i + run.length;
                 ++k) {
                if (junctions[k].first + junctions[k].second <= i + run.length) {
                    found[k] = j + (junctions[k].first - i);
                }
            }
        }
        i += run.op == AlignmentOp::Insertion ? 0 : run.length;
        j += run.op == AlignmentOp::Deletion ? 0 : run.length;
    }
    return found;
}

/**
 * Cut `consensus`, the consensus of the reads of a strand across a run, at the junctions of
 * `called`, the strand's called sequence across it.
 *
 * The two are aligned, and each junction of `called` must be matched whole, base for base, for
 * the bases on either side of it to belong to the members on either side. Where the best
 * alignment does not match a junction, the two are aligned again with the junctions pinned
 * (`align_global_pinned`), and that alignment is taken if it costs no more than the best one
 * plus `aligner.long_gap_open` for each junction the best one did not match: what opening one
 * more gap at each costs.
 *
 * @return The consensus with its junctions, or nothing if no such alignment is found.
 */
optional<RunSequence> cut_consensus(const RunSequence& called, const string& consensus,
                                    const StringAlignerParams& aligner) {
    const optional<StringAlignment> best = align_global(called.bases, consensus, aligner);
    if (!best) {
        return std::nullopt;
    }
    vector<optional<size_t>> at = matched_junctions(best->runs, called.junctions);
    auto unmatched = [](const optional<size_t>& column) { return !column.has_value(); };
    const auto missing = (size_t)std::ranges::count_if(at, unmatched);
    if (missing > 0) {
        vector<uint32_t> pinned(called.bases.size(), 0);
        for (size_t k = 0; k < called.junctions.size(); ++k) {
            std::fill_n(pinned.begin() + (ptrdiff_t)called.junctions[k].first,
                        called.junctions[k].second, (uint32_t)k + 1);
        }
        const optional<StringAlignment> kept = align_global_pinned(called.bases, consensus,
                                                                   pinned, aligner, best->band);
        if (!kept || kept->cost > best->cost + (int64_t)missing * aligner.long_gap_open) {
            return std::nullopt;
        }
        at = matched_junctions(kept->runs, called.junctions);
        if (std::ranges::any_of(at, unmatched)) {
            return std::nullopt;
        }
    }
    RunSequence cut;
    cut.bases = consensus;
    for (size_t k = 0; k < called.junctions.size(); ++k) {
        cut.junctions.emplace_back(*at[k], called.junctions[k].second);
    }
    return cut;
}

/// Return the bases `read` has between the bounds of each member of a run, given `runs`, an
/// alignment of `across` to `read`: its bases between the columns of the junctions of
/// `across`. The bases a read inserts beside a junction belong to the member on that side.
vector<string> cut_read(const RunSequence& across, const vector<AlignmentRun>& runs,
                        const string& read) {
    // The number of read bases before the column of each base of `across`, and through it.
    vector<size_t> before(across.bases.size());
    vector<size_t> through(across.bases.size());
    size_t i = 0;
    size_t j = 0;
    for (const AlignmentRun& run : runs) {
        if (run.op == AlignmentOp::Insertion) {
            j += run.length;
            continue;
        }
        for (size_t k = 0; k < run.length; ++k, ++i) {
            before[i] = j;
            j += run.op == AlignmentOp::Deletion ? 0 : 1;
            through[i] = j;
        }
    }
    vector<string> cut;
    size_t from = 0;
    for (const auto& [at, length] : across.junctions) {
        cut.push_back(read.substr(from, before[at] - from));
        from = through[at + length - 1];
    }
    cut.push_back(read.substr(from));
    return cut;
}

}

struct ReadConsensusStep::RunFinding {
    /// What the step found on one strand.
    struct Strand {
        /// The strand's spanning reads, and how many of them disagree.
        size_t reads = 0;
        size_t disagreeing = 0;
        bool candidate = false;
        bool built = false;
        bool unchanged = false;
        bool not_cuttable = false;
        /// The bases the consensus puts between the bounds of each member, where it was cut.
        vector<string> segments;
        /// Whether each of `segments` is no allele of its member.
        vector<char> is_new;
    };
    Strand strands[2];
    bool over_cap = false;
    /// Whether no strand has `min_spanning` reads that span the run.
    bool too_few = false;
    /// Whether the run has two strands and all its reads were taken to be on strand 0.
    bool pooled = false;
    /// The bases each read that spans the run has between the bounds of each member.
    vector<ReadInteriors> interiors;

    /// Return whether any strand has a consensus that was cut into the members.
    bool treats() const { return !strands[0].segments.empty() || !strands[1].segments.empty(); }
};

void ReadConsensusStep::grow_runs(const Parts& parts, const AlleleSplicer& splicer,
                                  const vector<StagedSite*>& snarls,
                                  const vector<Finding>& findings, vector<Run>& runs,
                                  vector<RunFinding>& run_findings) {
    const HandleGraph& graph = *parts.graph;
    const size_t none = std::numeric_limits<size_t>::max();
    // A site can be in a run where it is given no allele on its own.
    auto can_join = [&](size_t i) {
        const Finding& found = findings[i];
        return found.measured && !found.proposes() && found.longest <= parts.max_bases;
    };
    // A run grows from a site whose reads disagree with what is called: one with a candidate
    // strand, or one whose span summary shows a few bases of difference.
    auto is_seed = [&](size_t i) {
        return findings[i].candidate() || findings[i].signal >= options.run_site_diff;
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
    // The site after each site in its chain, where the two can be neighbours in a run, and the
    // nodes from the one's end bound to the other's start bound.
    vector<size_t> next(snarls.size(), none);
    vector<char> follows(snarls.size(), 0);
    vector<vector<handle_t>> nodes_after(snarls.size());
    vector<size_t> bases_after(snarls.size(), 0);
    for (size_t i = 0; i < snarls.size(); ++i) {
        if (!can_join(i)) {
            continue;
        }
        handle_t at = parts.sites->get_site_end(snarls[i]->net);
        vector<handle_t> nodes{at};
        size_t bases = graph.get_length(at);
        while (bases <= options.run_max_junction) {
            if (const auto found = by_start.find(at);
                found != by_start.end() && found->second != i) {
                if (same_phase(i, found->second) && !follows[found->second]) {
                    next[i] = found->second;
                    follows[found->second] = 1;
                    nodes_after[i] = nodes;
                    bases_after[i] = bases;
                }
                break;
            }
            // Go on to the next node where every walk of the chain takes it from here.
            handle_t only = at;
            size_t onward = 0;
            graph.follow_edges(at, false, [&](const handle_t& to) {
                only = to;
                return ++onward < 2;
            });
            size_t back = 0;
            if (onward == 1) {
                graph.follow_edges(only, true, [&](const handle_t&) { return ++back < 2; });
            }
            if (onward != 1 || back != 1) {
                break;
            }
            at = only;
            nodes.push_back(at);
            bases += graph.get_length(at);
        }
    }

    // The stretches of sites that follow each other, in which runs grow. A stretch is looked
    // at where it has a site with a candidate strand, or its sites' signals sum to `min_diff`.
    vector<vector<size_t>> stretches;
    for (size_t i = 0; i < snarls.size(); ++i) {
        if (!can_join(i) || follows[i] || next[i] == none) {
            continue;
        }
        vector<size_t> stretch;
        size_t signal = 0;
        bool candidate = false;
        bool seeded = false;
        for (size_t at = i; at != none; at = next[at]) {
            stretch.push_back(at);
            signal += findings[at].signal;
            candidate = candidate || findings[at].candidate();
            seeded = seeded || is_seed(at);
        }
        if (seeded && (candidate || signal >= options.min_diff || !options.prefilter)) {
            stretches.push_back(std::move(stretch));
        }
    }

    // Grow a run from each seed of a stretch, in the stretch's order: take the neighbour with
    // the greater signal, or on a tie the one on the side with fewer members, the later side if
    // those tie too, and test the run. Keep the first run that is treated, and stop where one
    // cannot be. The stretches share no site, so they are done at once.
    vector<vector<Run>> grown(stretches.size());
    vector<vector<RunFinding>> found_at(stretches.size());
    // The seeds of each stretch that the most bases the scorer aligns kept from growing at all.
    vector<size_t> capped_seeds(stretches.size(), 0);
#pragma omp parallel for schedule(dynamic, 1)
    for (size_t c = 0; c < stretches.size(); ++c) {
        const vector<size_t>& stretch = stretches[c];
        vector<char> used(stretch.size(), 0);
        // What was found at each run tested so far that was not treated, by its first and
        // last members, so that no run is tested or kept twice.
        map<pair<size_t, size_t>, RunFinding> tested;
        set<pair<size_t, size_t>> kept;
        for (size_t seed = 0; seed < stretch.size(); ++seed) {
            if (used[seed] || !is_seed(stretch[seed])) {
                continue;
            }
            size_t first = seed;
            size_t last = seed;
            size_t bases = findings[stretch[seed]].longest;
            Run run;
            RunFinding found;
            // Return the bases the neighbour on `side` would add, or `none` where it cannot
            // join.
            auto adds = [&](int side) {
                if (last - first + 1 >= options.run_max_members
                    || (side == 0 ? first == 0 || used[first - 1]
                                  : last + 1 == stretch.size() || used[last + 1])) {
                    return none;
                }
                const size_t more = side == 0 ? findings[stretch[first - 1]].longest
                                                    + bases_after[stretch[first - 1]]
                                              : findings[stretch[last + 1]].longest
                                                    + bases_after[stretch[last]];
                if (bases + more > parts.max_bases) {
                    run.capped = true;
                    return none;
                }
                return more;
            };
            // Return the run with the neighbour on `side` taken in, and what is found at it.
            auto with = [&](int side) {
                Run wider;
                wider.tests = run.tests;
                wider.capped = run.capped;
                const size_t from = first - (side == 0 ? 1 : 0);
                const size_t to = last + (side == 1 ? 1 : 0);
                wider.members.assign(stretch.begin() + (ptrdiff_t)from,
                                     stretch.begin() + (ptrdiff_t)to + 1);
                for (size_t m = from; m < to; ++m) {
                    wider.junctions.push_back(nodes_after[stretch[m]]);
                }
                if (const auto known = tested.find({from, to}); known != tested.end()) {
                    return make_pair(std::move(wider), known->second);
                }
                ++wider.tests;
                RunFinding at_wider = examine_run(parts, splicer, snarls, wider);
                if (!at_wider.treats()) {
                    tested.emplace(make_pair(from, to), at_wider);
                }
                return make_pair(std::move(wider), std::move(at_wider));
            };
            auto take = [&](int side, size_t more, pair<Run, RunFinding>&& wider) {
                bases += more;
                first -= side == 0 ? 1 : 0;
                last += side == 1 ? 1 : 0;
                run = std::move(wider.first);
                found = std::move(wider.second);
            };
            // Grow until a strand's consensus is cut into the members.
            while (!found.treats() && !found.over_cap && !found.too_few
                   && !found.strands[0].not_cuttable && !found.strands[1].not_cuttable) {
                const size_t more[2] = {adds(0), adds(1)};
                if (more[0] == none && more[1] == none) {
                    break;
                }
                int side = more[0] == none ? 1 : 0;
                if (more[0] != none && more[1] != none) {
                    const size_t before = findings[stretch[first - 1]].signal;
                    const size_t after = findings[stretch[last + 1]].signal;
                    side = after > before || (after == before && last - seed <= seed - first)
                               ? 1 : 0;
                }
                take(side, more[side], with(side));
            }
            // Then take in a neighbour while that gives the treated strands more reads that
            // disagree with what is called: reads whose difference the mapper put there.
            auto disagreeing = [](const RunFinding& at) {
                size_t reads = 0;
                for (const RunFinding::Strand& on : at.strands) {
                    reads += on.segments.empty() ? 0 : on.disagreeing;
                }
                return reads;
            };
            for (bool grew = found.treats(); grew;) {
                grew = false;
                for (const int side : {1, 0}) {
                    const size_t more = adds(side);
                    if (more == none) {
                        continue;
                    }
                    pair<Run, RunFinding> wider = with(side);
                    run.tests = wider.first.tests;
                    if (wider.second.treats()
                        && disagreeing(wider.second) > disagreeing(found)) {
                        take(side, more, std::move(wider));
                        grew = true;
                        break;
                    }
                }
            }
            if (found.treats()) {
                std::fill(used.begin() + (ptrdiff_t)first, used.begin() + (ptrdiff_t)last + 1, 1);
            }
            capped_seeds[c] += run.members.empty() && run.capped ? 1 : 0;
            if (!run.members.empty() && kept.insert({first, last}).second) {
                grown[c].push_back(std::move(run));
                found_at[c].push_back(std::move(found));
            }
        }
    }
    for (size_t c = 0; c < stretches.size(); ++c) {
        counted.runs_over_cap += capped_seeds[c];
        for (size_t r = 0; r < grown[c].size(); ++r) {
            ++counted.runs_grown;
            counted.runs_tested += grown[c][r].tests;
            counted.runs_over_cap += grown[c][r].capped ? 1 : 0;
            runs.push_back(std::move(grown[c][r]));
            run_findings.push_back(std::move(found_at[c][r]));
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

    // What each strand is called as across the run, and what each allele of each member has
    // between the member's bounds.
    RunSequence called[2];
    bool crosses[2] = {true, true};
    vector<vector<string>> allele_bases(count);
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
            if (!crosses[strand]) {
                continue;
            }
            if (m > 0) {
                called[strand].junctions.emplace_back(called[strand].bases.size(),
                                                      junction.size());
                called[strand].bases += junction;
            }
            called[strand].bases += *interior;
        }
        for (size_t a = 0; a < site.alleles.walks.size(); ++a) {
            optional<string> spelled = between(site.alleles.sequence(graph, a), before, after);
            if (!spelled || spelled->size() > parts.max_bases) {
                found.over_cap = true;
                return found;
            }
            allele_bases[m].push_back(std::move(*spelled));
        }
    }
    if (!crosses[0] && !crosses[1]) {
        return found;
    }
    for (int strand = 0; strand < 2; ++strand) {
        if (crosses[strand] && called[strand].bases.size() > parts.max_bases) {
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
            found.pooled = true;
            std::ranges::fill(strand_of, 0);
        }
    }

    ReadConsensusParams building;
    building.max_reads = options.max_reads;
    optional<RunSequence> consensus[2];
    for (int strand = 0; strand < 2; ++strand) {
        if (!crosses[strand]) {
            continue;
        }
        RunFinding::Strand& on = found.strands[strand];
        vector<string> bases;
        for (size_t r = 0; r < reads.size(); ++r) {
            if (strand_of[r] == strand && reads[r].length <= parts.max_bases) {
                bases.push_back(reads[r].bases);
                on.disagreeing += disagrees(reads[r].bases, called[strand].bases) ? 1 : 0;
            }
        }
        on.reads = bases.size();
        on.candidate = on.reads >= options.min_spanning && on.disagreeing >= options.min_reads
                       && (double)on.disagreeing >= options.min_share * (double)on.reads;
        if (!on.candidate) {
            continue;
        }
        const optional<string> built = read_consensus(bases, called[strand].bases, building);
        const optional<int64_t> change =
            built ? edit_distance(*built, called[strand].bases) : std::nullopt;
        if (!change) {
            continue;
        }
        on.built = true;
        if ((size_t)*change <= options.min_change) {
            on.unchanged = true;
            continue;
        }
        consensus[strand] = cut_consensus(called[strand], *built, building.aligner);
        if (!consensus[strand]) {
            on.not_cuttable = true;
            continue;
        }
        on.segments = consensus[strand]->pieces();
        for (size_t m = 0; m < count; ++m) {
            on.is_new.push_back(std::ranges::find(allele_bases[m], on.segments[m])
                                == allele_bases[m].end() ? 1 : 0);
        }
    }
    found.too_few = found.strands[0].reads < options.min_spanning
                    && found.strands[1].reads < options.min_spanning;
    if (!found.treats()) {
        return found;
    }

    // The scorer aligns the spanning reads to every allele of each member, within its limits.
    for (size_t m = 0; m < count; ++m) {
        size_t alleles = allele_bases[m].size();
        for (int strand = 0; strand < 2; ++strand) {
            const RunFinding::Strand& on = found.strands[strand];
            const bool repeats = strand == 1 && !found.strands[0].segments.empty()
                                 && !on.segments.empty()
                                 && found.strands[0].segments[m] == on.segments[m];
            alleles += !on.segments.empty() && on.is_new[m] && !repeats ? 1 : 0;
        }
        if (reads.size() * alleles > parts.max_alignments) {
            found.over_cap = true;
            for (RunFinding::Strand& on : found.strands) {
                on.segments.clear();
            }
            return found;
        }
    }

    // Place each read's bases in the members in one way: align them to each sequence known
    // across the run, and cut them at the junctions of the one they fit best, the earlier of
    // two they fit equally.
    vector<const RunSequence*> known;
    for (int strand = 0; strand < 2; ++strand) {
        if (crosses[strand]) {
            known.push_back(&called[strand]);
        }
    }
    for (int strand = 0; strand < 2; ++strand) {
        if (consensus[strand]) {
            known.push_back(&*consensus[strand]);
        }
    }
    found.interiors.resize(count);
    for (const AlleleLikelihoodCalculator::SpanningRead& read : reads) {
        if (read.length > parts.max_bases) {
            continue;
        }
        optional<StringAlignment> best;
        const RunSequence* fits = nullptr;
        for (const RunSequence* sequence : known) {
            if (fits != nullptr && sequence->bases == fits->bases) {
                continue;
            }
            optional<StringAlignment> alignment = align_global(sequence->bases, read.bases,
                                                               building.aligner);
            if (alignment && (!best || alignment->cost < best->cost)) {
                best = std::move(alignment);
                fits = sequence;
            }
        }
        if (fits == nullptr) {
            continue;
        }
        const vector<string> pieces = cut_read(*fits, best->runs, read.bases);
        for (size_t m = 0; m < count; ++m) {
            found.interiors[m].bases[read.key] = pieces[m];
        }
    }
    return found;
}

size_t ReadConsensusStep::run(const Parts& parts) {
    counted = ReadConsensusCounts();
    changed.clear();
    run_changed.clear();
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

    AlleleLikelihoodCalculator& calculator = parts.genotyper->calculator();
    calculator.set_counted_fetches(false);
    vector<Finding> findings(snarls.size());
#pragma omp parallel for schedule(dynamic, 64)
    for (size_t i = 0; i < snarls.size(); ++i) {
        findings[i] = examine(parts, splicer, *snarls[i]);
    }

    // Grow the runs around the sites whose reads disagree and that gain no allele on their
    // own.
    vector<Run> runs;
    vector<RunFinding> run_findings;
    if (options.runs) {
        grow_runs(parts, splicer, snarls, findings, runs, run_findings);
    }

    // Give each difference to the innermost site that holds it: decide the deepest sites first,
    // and leave a strand alone where a site under it gained an allele for the strand.
    struct Proposal {
        size_t level;
        bool is_run;
        size_t index;
    };
    vector<Proposal> proposing;
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
            proposing.push_back(Proposal{snarls[i]->placement.level, false, i});
        }
    }
    for (size_t r = 0; r < runs.size(); ++r) {
        const RunFinding& found = run_findings[r];
        counted.runs_over_cap += found.over_cap ? 1 : 0;
        for (const RunFinding::Strand& on : found.strands) {
            counted.run_candidate_strands += on.candidate ? 1 : 0;
            counted.run_built += on.built ? 1 : 0;
            counted.no_consensus += on.candidate && !on.built ? 1 : 0;
            counted.run_unchanged += on.unchanged ? 1 : 0;
            counted.runs_not_cuttable += on.not_cuttable ? 1 : 0;
        }
        if (found.treats()) {
            proposing.push_back(
                Proposal{snarls[runs[r].members.front()]->placement.level, true, r});
        }
    }
    std::ranges::stable_sort(proposing, [](const Proposal& a, const Proposal& b) {
        return a.level > b.level;
    });
    // The strands each site gains an allele for, as a mask, by record key. A member of a
    // treated run gains for the strands the run is treated on.
    unordered_map<size_t, int> gained;
    // The sites to genotype again: each with its index, and its run and place in it if it is a
    // member of one.
    struct Task {
        size_t site;
        size_t run;
        size_t member;
    };
    const size_t no_run = std::numeric_limits<size_t>::max();
    vector<Task> tasks;
    auto gained_under = [&](size_t i) {
        int under = 0;
        splicer.for_each_under(snarls[i]->record_key, [&](const StagedSite& nested) {
            if (const auto mask = gained.find(nested.record_key); mask != gained.end()) {
                under |= mask->second;
            }
        });
        return under;
    };
    for (const Proposal& proposal : proposing) {
        if (proposal.is_run) {
            const Run& run = runs[proposal.index];
            RunFinding& found = run_findings[proposal.index];
            int under = 0;
            for (const size_t i : run.members) {
                under |= gained_under(i);
            }
            int mask = 0;
            for (int strand = 0; strand < 2; ++strand) {
                RunFinding::Strand& on = found.strands[strand];
                if (on.segments.empty()) {
                    continue;
                }
                if ((found.pooled ? under : under & (1 << strand)) != 0) {
                    on.segments.clear();
                    ++counted.covered_by_nested;
                } else {
                    mask |= found.pooled ? 3 : 1 << strand;
                }
            }
            if (mask != 0) {
                ++counted.runs_treated;
                for (size_t m = 0; m < run.members.size(); ++m) {
                    gained[snarls[run.members[m]]->record_key] = mask;
                    tasks.push_back(Task{run.members[m], proposal.index, m});
                }
            }
            continue;
        }
        const size_t i = proposal.index;
        Finding& found = findings[i];
        const int under = gained_under(i);
        int mask = 0;
        for (int strand = 0; strand < 2; ++strand) {
            Finding::Strand& on = found.strands[strand];
            if (!on.allele) {
                continue;
            }
            // Pooled reads are those of both strands.
            if ((found.pooled ? under : under & (1 << strand)) != 0) {
                on.allele.reset();
                ++counted.covered_by_nested;
            } else {
                mask |= found.pooled ? 3 : 1 << strand;
            }
        }
        if (mask != 0) {
            gained[snarls[i]->record_key] = mask;
            tasks.push_back(Task{i, no_run, 0});
        }
    }
    std::ranges::sort(tasks, {}, &Task::site);

    // Add the alleles and genotype each site again. The sites differ, so they can be done at
    // once. A member of a run is genotyped from the bases its reads were given there.
    vector<char> genotyped(tasks.size(), 0);
    vector<size_t> added(tasks.size(), 0);
#pragma omp parallel for schedule(dynamic, 1)
    for (size_t t = 0; t < tasks.size(); ++t) {
        const Task& task = tasks[t];
        StagedSite& site = *snarls[task.site];
        const CandidateAlleles kept = site.alleles;
        const size_t first_added = site.alleles.walks.size();
        // The bases each strand's consensus has between the bounds, where they are a new
        // allele.
        const string* alleles[2] = {nullptr, nullptr};
        for (int strand = 0; strand < 2; ++strand) {
            if (task.run == no_run) {
                const optional<string>& allele = findings[task.site].strands[strand].allele;
                alleles[strand] = allele ? &*allele : nullptr;
            } else {
                const RunFinding::Strand& on = run_findings[task.run].strands[strand];
                alleles[strand] = !on.segments.empty() && on.is_new[task.member]
                                      ? &on.segments[task.member] : nullptr;
            }
        }
        for (int strand = 0; strand < 2; ++strand) {
            if (alleles[strand] != nullptr
                && !(strand == 1 && alleles[0] != nullptr && *alleles[0] == *alleles[1])) {
                added[t] += add_sequence_allele(*parts.graph, site.alleles, *alleles[strand])
                                ? 1 : 0;
            }
        }
        const bool again = task.run == no_run
                               ? added[t] > 0 && parts.tree_genotyper->regenotype(site, first_added)
                               : parts.tree_genotyper->regenotype(
                                     site, first_added,
                                     &run_findings[task.run].interiors[task.member]);
        if (again) {
            genotyped[t] = 1;
        } else {
            site.alleles = kept;
        }
    }
    calculator.set_counted_fetches(true);

    // File each site again, one at a time, in the order of the sites.
    LinkageCollector* model = parts.linker->collector();
    size_t filed = 0;
    for (size_t t = 0; t < tasks.size(); ++t) {
        StagedSite& site = *snarls[tasks[t].site];
        const bool in_run = tasks[t].run != no_run;
        if (!genotyped[t]) {
            counted.regenotype_failed += in_run || added[t] > 0 ? 1 : 0;
            continue;
        }
        ++filed;
        if (!in_run) {
            counted.added += added[t];
            ++counted.sites_with_added;
            changed.push_back(site.record_key);
        } else if (added[t] > 0) {
            ++counted.run_sites_gained;
            run_changed.push_back(site.record_key);
        }
        if (model->has_entry(site.record_key)) {
            // Retract before filing, because the model finds the first live entry for a key.
            model->retract(site.record_key);
            parts.linker->file(site, site.placement.level > 0 ? site.ploidy == 1
                                                              : site.placement.one_copy);
        }
    }

    // The phase is frozen again before the render, once it is final.
    parts.phase_table->freeze_for_render(false);
    return filed;
}

void ReadConsensusStep::count_called(const Parts& parts) {
    counted.sites_called = 0;
    counted.run_sites_called = 0;
    if (changed.empty() && run_changed.empty()) {
        return;
    }
    const unordered_map<size_t, StagedSite*> by_key = parts.staged_sites->by_key();
    auto called = [&](const vector<size_t>& keys) {
        size_t count = 0;
        for (const size_t key : keys) {
            const auto found = by_key.find(key);
            if (found == by_key.end() || found->second->dropped) {
                continue;
            }
            const StagedSite& site = *found->second;
            const vector<int> genotype = parts.linker->chosen_genotype(site);
            count += std::ranges::any_of(genotype, [&](int allele) {
                return is_consensus_allele(*parts.sites, site, allele);
            }) ? 1 : 0;
        }
        return count;
    };
    counted.sites_called = called(changed);
    counted.run_sites_called = called(run_changed);
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
        << " already an allele, " << counted.regenotype_failed << " sites not re-genotyped"
        << "; runs: " << counted.runs_grown << " grown, " << counted.runs_tested << " tested, "
        << counted.run_candidate_strands << " candidate strands, " << counted.run_built
        << " consensuses built, " << counted.run_unchanged << " within " << options.min_change
        << " bases of the call, " << counted.runs_not_cuttable << " not cuttable, "
        << counted.runs_over_cap << " over the cap, " << counted.runs_treated << " treated, "
        << counted.run_sites_gained << " member sites gained an allele, "
        << counted.run_sites_called << " member sites called with a consensus allele";
    if (!options.prefilter) {
        out << "; prefilter off: " << counted.candidates_outside_prefilter << " of "
            << counted.candidate_sites << " candidate sites fail the prefilter";
    }
    out << endl;
}

}
