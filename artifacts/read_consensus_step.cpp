#include <algorithm>
#include <cmath>
#include <limits>
#include <map>
#include <memory>
#include <numeric>
#include <unordered_set>

#include "read_consensus_step.hpp"
#include "site_values.hpp"

namespace vg::multipass {

using namespace std;

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

struct ReadConsensusStep::Run {
    /// The members, as indices of the phased snarls, in the order of their chain. A snarl on
    /// its own is a run of one member.
    vector<size_t> members;
    /// The nodes from each member's end bound to the next member's start bound, both included.
    /// It is the one node where the two share a bound.
    vector<vector<handle_t>> junctions;
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
    if (called.junctions.empty()) {
        // A run of one member has no junction to cut at.
        return RunSequence{consensus, {}};
    }
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
        bool not_cuttable = false;
        /// The bases the consensus puts between the bounds of each member, where it was cut.
        vector<string> segments;
        /// Whether each of `segments` is no allele of its member.
        vector<char> is_new;
        /// The bases the consensus puts between the run's two outer bounds, where the run is
        /// to become one site.
        optional<string> whole;

        /// Return whether the strand has a consensus that the run is treated for.
        bool treated() const { return !segments.empty() || whole.has_value(); }
    };
    Strand strands[2];
    /// The bases each strand is called as between the run's two outer bounds, where it crosses
    /// every member.
    optional<string> called[2];
    /// The record keys of the edit sites inside the nodes between the members.
    vector<size_t> edits_between;
    /// Whether the run's reads were fetched and compared with the called sequences.
    bool fetched = false;
    /// Whether a called sequence or an allele of a member has more bases between the bounds
    /// than the scorer aligns, or the reads and alleles need more alignments than it makes.
    bool over_cap = false;
    /// Whether no strand has `min_spanning` reads that span the run.
    bool too_few = false;
    /// Whether the run has two strands and all its reads were taken to be on strand 0.
    bool pooled = false;
    /// The bases each read that spans the run has between the bounds of each member.
    vector<ReadInteriors> interiors;

    /// Return whether any strand has a consensus that the run is treated for.
    bool treats() const { return strands[0].treated() || strands[1].treated(); }

    /// Return whether the run is to become one site.
    bool merges() const { return strands[0].whole.has_value() || strands[1].whole.has_value(); }
};

struct ReadConsensusStep::Alone {
    /// Whether the site has a called sequence and enough spanning reads to be in a run.
    bool measured = false;
    /// Whether the span summary says a read of the site can disagree with a called sequence.
    bool can_disagree = false;
    /// The most bases a called sequence has between the bounds.
    size_t longest = 0;
    /// The number of bases by which `min_reads` of the spanning reads have more between the
    /// bounds than the longest called sequence, or fewer than the shortest, whichever is more.
    size_t signal = 0;
    /// The longest gap that `min_reads` of the spanning reads have inside the node of the start
    /// bound, and inside that of the end bound.
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
    if (spans.reads < options.min_spanning) {
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
    alone.signal = std::max(
        spans.rank_bases > alone.longest ? spans.rank_bases - alone.longest : 0,
        shortest > spans.low_bases ? shortest - spans.low_bases : 0);
    // No read can disagree where neither a called sequence nor enough reads have `min_diff`
    // bases between the bounds.
    alone.can_disagree = alone.longest >= options.min_diff
                         || spans.rank_bases >= options.min_diff;
    return alone;
}

void ReadConsensusStep::grow_runs(const Parts& parts, const AlleleSplicer& splicer,
                                  const vector<StagedSite*>& snarls,
                                  const vector<Alone>& alone, vector<Run>& runs,
                                  vector<RunFinding>& run_findings) {
    const HandleGraph& graph = *parts.graph;
    const size_t none = std::numeric_limits<size_t>::max();
    // A site can be in a run where it is given no allele on its own.
    auto can_join = [&](size_t i) {
        return alone[i].measured && !alone[i].treated && alone[i].longest <= max_bases;
    };
    // The signal of each site: that of its own span summary, or the gap its reads have inside
    // a node between it and a neighbour, if that is `min_diff` bases and more.
    vector<size_t> signal(snarls.size());
    for (size_t i = 0; i < snarls.size(); ++i) {
        signal[i] = alone[i].signal;
    }
    // A run grows from a site whose reads disagree with what is called: one with a candidate
    // strand, or one whose signal is a few bases of difference.
    auto is_seed = [&](size_t i) {
        return alone[i].candidates > 0 || signal[i] >= options.run_site_diff;
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
        while (true) {
            if (const auto found = by_start.find(at);
                found != by_start.end() && found->second != i) {
                // The two are neighbours across few bases, or across any number where the
                // reads have a gap inside the node of a bound that faces the other.
                const size_t gap = std::max(alone[i].end_gap, alone[found->second].start_gap);
                const bool gapped = gap >= options.min_diff;
                if ((bases <= options.run_max_junction || gapped)
                    && same_phase(i, found->second) && !follows[found->second]) {
                    next[i] = found->second;
                    follows[found->second] = 1;
                    nodes_after[i] = nodes;
                    bases_after[i] = bases;
                    if (gapped) {
                        signal[i] = std::max(signal[i], gap);
                        signal[found->second] = std::max(signal[found->second], gap);
                    }
                }
                break;
            }
            if (bases > options.run_max_junction) {
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
        size_t signal_sum = 0;
        bool candidate = false;
        bool seeded = false;
        for (size_t at = i; at != none; at = next[at]) {
            stretch.push_back(at);
            signal_sum += signal[at];
            candidate = candidate || alone[at].candidates > 0;
            seeded = seeded || is_seed(at);
        }
        if (seeded && (candidate || signal_sum >= options.min_diff)) {
            stretches.push_back(std::move(stretch));
        }
    }

    // Grow a run from each seed of a stretch, in the stretch's order: take the neighbour with
    // the greater signal, or on a tie the one on the side with fewer members, the later side if
    // those tie too, and test the run. Keep the first run that is treated, and stop where one
    // cannot be. The stretches share no site, so they are done at once.
    vector<vector<Run>> grown(stretches.size());
    vector<vector<RunFinding>> found_at(stretches.size());
#pragma omp parallel for schedule(dynamic, 1)
    for (size_t c = 0; c < stretches.size(); ++c) {
        const vector<size_t>& stretch = stretches[c];
        vector<char> used(stretch.size(), 0);
        // What was found at each run tested so far that was not treated, by its first and
        // last members, so that no run is tested twice.
        map<pair<size_t, size_t>, RunFinding> tested;
        for (size_t seed = 0; seed < stretch.size(); ++seed) {
            if (used[seed] || !is_seed(stretch[seed])) {
                continue;
            }
            size_t first = seed;
            size_t last = seed;
            size_t bases = alone[stretch[seed]].longest;
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
                const size_t more = side == 0 ? alone[stretch[first - 1]].longest
                                                    + bases_after[stretch[first - 1]]
                                              : alone[stretch[last + 1]].longest
                                                    + bases_after[stretch[last]];
                return bases + more > max_bases ? none : more;
            };
            // Return the run with the neighbour on `side` taken in, and what is found at it.
            auto with = [&](int side) {
                Run wider;
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
                    const size_t before = signal[stretch[first - 1]];
                    const size_t after = signal[stretch[last + 1]];
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
                    reads += on.treated() ? on.disagreeing : 0;
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
                grown[c].push_back(std::move(run));
                found_at[c].push_back(std::move(found));
            }
        }
    }
    for (size_t c = 0; c < stretches.size(); ++c) {
        for (size_t r = 0; r < grown[c].size(); ++r) {
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
    // Whether an allele of a member has more bases between its bounds than the scorer aligns.
    bool long_allele = false;
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
            long_allele = long_allele || !spelled || spelled->size() > max_bases;
            allele_bases[m].push_back(spelled ? std::move(*spelled) : string());
        }
    }
    if (!crosses[0] && !crosses[1]) {
        found.over_cap = long_allele;
        return found;
    }
    for (int strand = 0; strand < 2; ++strand) {
        if (crosses[strand] && called[strand].bases.size() > max_bases) {
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

    const ReadConsensusParams building;
    optional<RunSequence> consensus[2];
    // Each strand's consensus across the run, where it is to become an allele.
    optional<string> whole[2];
    for (int strand = 0; strand < 2; ++strand) {
        if (!crosses[strand]) {
            continue;
        }
        RunFinding::Strand& on = found.strands[strand];
        vector<string> bases;
        for (size_t r = 0; r < reads.size(); ++r) {
            if (strand_of[r] == strand && reads[r].length <= max_bases) {
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
        if (!change || (size_t)*change <= options.min_change) {
            continue;
        }
        whole[strand] = *built;
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
    if (long_allele) {
        // The scorer cannot align the reads to every allele of each member.
        found.over_cap = true;
        for (RunFinding::Strand& on : found.strands) {
            on.segments.clear();
        }
        return found;
    }
    if (count == 1) {
        // A snarl on its own is treated only for a consensus that none of its alleles spells.
        // The members of a longer run have their reads placed again whatever they gain.
        for (RunFinding::Strand& on : found.strands) {
            if (!on.segments.empty() && !on.is_new[0]) {
                on.segments.clear();
            }
        }
    }
    if ((found.strands[0].not_cuttable || found.strands[1].not_cuttable)
        && can_merge(parts, snarls, run, found)) {
        // The run becomes one site, with each strand's consensus whole. The scorer aligns the
        // spanning reads to each allele of that site: the reference, what each strand is
        // called as, and each consensus.
        const size_t most_alleles = 5;
        for (int strand = 0; strand < 2; ++strand) {
            RunFinding::Strand& on = found.strands[strand];
            on.segments.clear();
            on.is_new.clear();
            if (crosses[strand]) {
                found.called[strand] = called[strand].bases;
            }
            if (whole[strand] && whole[strand]->size() <= max_bases
                && reads.size() * most_alleles <= max_alignments) {
                on.whole = std::move(whole[strand]);
            } else if (whole[strand]) {
                found.over_cap = true;
            }
        }
        return found;
    }
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
        if (reads.size() * alleles > max_alignments) {
            found.over_cap = true;
            for (RunFinding::Strand& on : found.strands) {
                on.segments.clear();
            }
            return found;
        }
    }

    if (count == 1) {
        // A snarl on its own is scored on the bases each read has between its bounds.
        return found;
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
        if (read.length > max_bases) {
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
                                                  const RunFinding& found) const {
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
    for (int strand = 0; strand < 2; ++strand) {
        if (found.called[strand]) {
            add(*found.called[strand], joined([&](const StagedSite& member) {
                    return splicer.strand_traversal(member, strand);
                }));
        }
    }
    // The alleles so far spell what was called at the member sites, some of them as sequences.
    const size_t first_consensus = across.walks.size();
    for (const RunFinding::Strand& on : found.strands) {
        if (on.whole) {
            add(*on.whole, std::nullopt);
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
    calculator.set_counted_fetches(false);
    // Test each site on its own, as a run of one, where a read of it can disagree.
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

    // The runs that are treated: the sites with a consensus of their own, and the runs grown
    // around the sites whose reads disagree and that gain no allele on their own.
    vector<Run> runs;
    vector<RunFinding> run_findings;
    for (size_t i = 0; i < snarls.size(); ++i) {
        counted.tested += alone[i].tested ? 1 : 0;
        counted.over_cap += alone[i].over_cap ? 1 : 0;
        counted.candidate_strands += alone[i].candidates;
        if (found_alone[i]) {
            runs.push_back(Run{{i}, {}});
            run_findings.push_back(std::move(*found_alone[i]));
        }
    }
    found_alone.clear();
    grow_runs(parts, splicer, snarls, alone, runs, run_findings);

    // Give each difference to the innermost site that holds it: decide the deepest runs first,
    // and leave a strand alone where a site under it gained an allele for the strand.
    vector<size_t> proposing(runs.size());
    std::iota(proposing.begin(), proposing.end(), 0);
    std::ranges::stable_sort(proposing, [&](size_t a, size_t b) {
        return snarls[runs[a].members.front()]->placement.level
               > snarls[runs[b].members.front()]->placement.level;
    });
    // The strands each site gains an allele for, as a mask, by record key. A member of a
    // treated run gains for the strands the run is treated on.
    unordered_map<size_t, int> gained;
    // The sites to genotype again: each with its index, its run and its place in the run.
    struct Task {
        size_t site;
        size_t run;
        size_t member;
    };
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
    for (const size_t r : proposing) {
        const Run& run = runs[r];
        RunFinding& found = run_findings[r];
        int under = 0;
        for (const size_t i : run.members) {
            under |= gained_under(i);
        }
        // One site takes the place of everything under its members, so a run that is to become
        // one is left alone where anything under it gained an allele. Pooled reads are those of
        // both strands.
        const bool merges = found.merges();
        int mask = 0;
        for (int strand = 0; strand < 2; ++strand) {
            RunFinding::Strand& on = found.strands[strand];
            if (!on.treated()) {
                continue;
            }
            if ((found.pooled || merges ? under : under & (1 << strand)) != 0) {
                on.segments.clear();
                on.whole.reset();
                ++counted.covered_by_nested;
            } else {
                mask |= found.pooled ? 3 : 1 << strand;
            }
        }
        if (mask == 0) {
            continue;
        }
        counted.runs_treated += run.members.size() > 1 ? 1 : 0;
        for (size_t m = 0; m < run.members.size(); ++m) {
            gained[snarls[run.members[m]]->record_key] = mask;
            if (!merges) {
                tasks.push_back(Task{run.members[m], r, m});
            }
        }
        if (merges) {
            // The one site is made in the place of the run's first member.
            tasks.push_back(Task{run.members.front(), r, 0});
        }
    }
    std::ranges::sort(tasks, {}, &Task::site);

    // Add the alleles and genotype each site again. The sites differ, so they can be done at
    // once.
    vector<char> genotyped(tasks.size(), 0);
    vector<size_t> added(tasks.size(), 0);
    // The one site of each run that becomes one.
    vector<optional<StagedSite>> merged(tasks.size());
#pragma omp parallel for schedule(dynamic, 1)
    for (size_t t = 0; t < tasks.size(); ++t) {
        const Task& task = tasks[t];
        const Run& run = runs[task.run];
        const RunFinding& found = run_findings[task.run];
        if (found.merges()) {
            merged[t] = merge_run(parts, splicer, snarls, run, found);
            genotyped[t] = merged[t] ? 1 : 0;
            continue;
        }
        StagedSite& site = *snarls[task.site];
        const CandidateAlleles kept = site.alleles;
        const size_t first_added = site.alleles.walks.size();
        // The bases each strand's consensus has between the bounds, where they are a new
        // allele.
        const string* alleles[2] = {nullptr, nullptr};
        for (int strand = 0; strand < 2; ++strand) {
            const RunFinding::Strand& on = found.strands[strand];
            alleles[strand] = !on.segments.empty() && on.is_new[task.member]
                                  ? &on.segments[task.member] : nullptr;
        }
        for (int strand = 0; strand < 2; ++strand) {
            if (alleles[strand] != nullptr
                && !(strand == 1 && alleles[0] != nullptr && *alleles[0] == *alleles[1])) {
                added[t] += add_sequence_allele(*parts.graph, site.alleles, *alleles[strand])
                                ? 1 : 0;
            }
        }
        // A member of a run of two or more is genotyped from the bases its reads were given
        // there, whether or not it gained an allele. A snarl on its own is genotyped again
        // where it gained one.
        const bool alone_in_run = run.members.size() == 1;
        const bool again = (!alone_in_run || added[t] > 0)
                           && parts.tree_genotyper->regenotype(
                                  site, first_added,
                                  alone_in_run ? nullptr : &found.interiors[task.member]);
        if (again) {
            genotyped[t] = 1;
            if (added[t] > 0) {
                site.first_consensus_allele = first_added;
            }
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
        StagedSite& site = *snarls[tasks[t].site];
        if (!genotyped[t]) {
            continue;
        }
        ++filed;
        if (merged[t]) {
            // The one site takes the place of the run's first member. The other members, the
            // sites under any member and the edit sites between them are removed below.
            const Run& run = runs[tasks[t].run];
            ++counted.runs_merged;
            for (const size_t i : run.members) {
                absorbed.insert(snarls[i]->record_key);
            }
            absorbed.insert(run_findings[tasks[t].run].edits_between.begin(),
                            run_findings[tasks[t].run].edits_between.end());
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

    // Remove what the merged sites take the place of. This moves the staged sites, so nothing
    // above is read after it.
    edit_sites_on.clear();
    if (!absorbed.empty()) {
        parts.staged_sites->remove_with_nested(
            absorbed, [&](const StagedSite& gone) { model->retract(gone.record_key); });
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
        << " of them called with a consensus allele; left alone: " << counted.over_cap
        << " snarls over a limit of the scorer, " << counted.covered_by_nested
        << " strands covered by a nested site" << endl;
}

}
