#include "read_phasing.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <numeric>
#include <optional>
#include <unordered_map>

namespace vg {

using std::log10;
using std::max;
using std::min;
using std::pair;
using std::unordered_map;

double PanelErrorModel::correlation(size_t distance) const {
    if (!fitted) {
        return 0.0;
    }
    const double kept = 1.0 - 2.0 * flip_rate;
    return kept * kept * std::exp(-2.0 * switch_rate * (double)distance);
}

double PanelErrorModel::vote(size_t distance) const {
    const double c = correlation(distance);
    return log10((1.0 + c) / (1.0 - c));
}

double read_confidence(float q0, float c) {
    const double win = max((double)q0, 1.0 - (double)q0) * (double)c;
    return -10.0 * log10(max(1.0 - win, 1e-12));
}

void merge_mates(PhaseSite& site) {
    const size_t n = site.read_key.size();
    vector<size_t> order(n);
    std::iota(order.begin(), order.end(), 0);
    // Sort by key, and within a key put the preferred row first. Every tie is broken on the rows'
    // values, so two rows that compare equal are identical and either may be kept.
    std::ranges::sort(order, [&](size_t x, size_t y) {
        if (site.read_key[x] != site.read_key[y]) {
            return site.read_key[x] < site.read_key[y];
        }
        const double cx = read_confidence(site.q0[x], site.c[x]);
        if (const double cy = read_confidence(site.q0[y], site.c[y]); cx != cy) {
            return cx > cy;
        }
        if (site.q0[x] != site.q0[y]) {
            return site.q0[x] > site.q0[y];
        }
        return site.c[x] > site.c[y];
    });
    vector<uint64_t> keys;
    vector<float> q0s;
    vector<float> cs;
    keys.reserve(n);
    q0s.reserve(n);
    cs.reserve(n);
    for (size_t r : order) {
        if (!keys.empty() && keys.back() == site.read_key[r]) {
            continue;
        }
        keys.push_back(site.read_key[r]);
        q0s.push_back(site.q0[r]);
        cs.push_back(site.c[r]);
    }
    site.read_key.swap(keys);
    site.q0.swap(q0s);
    site.c.swap(cs);
}

PhaseSite reduce_to_pair(const PhaseReadEvidence& pe, size_t a0, size_t a1) {
    PhaseSite site;
    // Slot 0 holds `a0` and slot 1 `a1`, so the weights follow the alleles into their slots.
    const vector<double> weight = allele_length_weights(
        pe.allele_length, pe.n_alleles, pe.mean_read_length, pe.length_weighted,
        vector<int>{(int)a0, (int)a1});
    for (size_t r = 0; r < pe.num_reads(); ++r) {
        const double e = pe.mismap[r];
        const double r0 = (1.0 - e) * weight[0] * (double)pe.rel_at(r, a0);
        const double r1 = (1.0 - e) * weight[1] * (double)pe.rel_at(r, a1);
        const double inside = r0 + r1;
        if (inside <= 0.0) {
            // The read fits neither chosen allele, so it says nothing about their order.
            continue;
        }
        site.read_key.push_back(pe.read_key[r]);
        site.q0.push_back((float)(r0 / inside));
        site.c.push_back((float)(inside / (inside + e)));
    }
    merge_mates(site);
    double score_sum = 0.0;
    for (size_t i = 0; i < site.read_key.size(); ++i) {
        score_sum += read_confidence(site.q0[i], site.c[i]);
    }
    if (!site.read_key.empty()) {
        site.reliability = score_sum / (double)site.read_key.size();
    }
    return site;
}

double phase_link(const PhaseSite& a, const PhaseSite& b, double cap) {
    // Both read lists are sorted by key, one row per key, so we merge them rather than look up
    // each read.
    double total = 0.0;
    size_t i = 0;
    size_t j = 0;
    while (i < a.read_key.size() && j < b.read_key.size()) {
        if (a.read_key[i] < b.read_key[j]) {
            ++i;
        } else if (b.read_key[j] < a.read_key[i]) {
            ++j;
        } else {
            const double qa = a.q0[i];
            const double qb = b.q0[j];
            // Probability that the read's two alleles lie on one strand in the current orders.
            // These sum to 1, so the ratio below needs no normalisation.
            const double same = qa * qb + (1.0 - qa) * (1.0 - qb);
            const double diff = 1.0 - same;
            // With probability `pr` the read reports the true relation, and otherwise a coin flip,
            // so that one mismapped read cannot dominate the sum.
            const double pr = (double)a.c[i] * (double)b.c[j];
            const double cis = pr * same + (1.0 - pr) * 0.5;
            if (const double trans = pr * diff + (1.0 - pr) * 0.5; cis > 0.0 && trans > 0.0) {
                total += log10(cis / trans);
            }
            ++i;
            ++j;
        }
    }
    if (cap > 0.0) {
        total = max(-cap, min(cap, total));
    }
    return total;
}

namespace {

/// Return the likelihood of a read's allele at a site if the read came from the strand of slot 0,
/// for `q` the probability that it carries the allele in slot 0 and `c` the probability that it
/// came from one of the two strands. For the strand of slot 1, pass 1 - q.
double strand_likelihood(double q, double c) {
    return c * q + (1.0 - c) * 0.5;
}

/// Add one phase set's counters into the totals for all phase sets.
///
/// `sites` is counted once for all phase sets, and `strands_rederived` belongs to the caller, so
/// neither is added here.
void add_phase_set_counters(const ReadPhasingCounters& from, ReadPhasingCounters& into) {
    into.chains += from.chains;
    into.reliable += from.reliable;
    into.breaks += from.breaks;
    into.breaks_no_reads += from.breaks_no_reads;
    into.hung += from.hung;
    into.hung_no_reads += from.hung_no_reads;
    into.flipped += from.flipped;
    into.demoted_incoherent += from.demoted_incoherent;
    into.coherence_unconverged += from.coherence_unconverged;
    into.coherence_rounds_run = max(into.coherence_rounds_run, from.coherence_rounds_run);
}

/// The fewest pairs of sites that `fit_panel_errors` fits rates to.
constexpr size_t PANEL_FIT_MIN_PAIRS = 200;
/// The number of bins, by distance, of the fit's table.
constexpr size_t PANEL_FIT_BINS = 8;

/// What phasing one phase set produces, gathered across the phase sets afterwards.
struct PhaseSetResult {
    ReadPhasingCounters counters;
    /// The record keys of the sites to flip.
    vector<size_t> flipped;
    /// The record keys of the sites the panel does not order whose reads ordered them.
    vector<size_t> read_ordered;
    /// For each site, its run (see `PhaseSite::run`).
    vector<size_t> run;
    /// For each run, its `PhaseSite::run_correlation`.
    vector<double> run_correlation;
    /// For each pair of phase-chain sites of one piece that the panel orders, with no such site
    /// between them, the distance between them and whether the reads order them the other way
    /// round from the panel.
    vector<pair<size_t, bool>> panel_pairs;
};

/// A tally of the votes on whether to flip one site.
struct FlipVote {
    /// Sum of the votes, in log10 units, positive for keeping the site's order and negative for
    /// flipping it.
    double total = 0.0;
    /// How many of the votes from reads were not 0: the link to one site, or the vote of one run's
    /// sites.
    size_t links = 0;
    /// The part of `total` from the reads' votes from each run's sites (see `PhaseSite::run`).
    vector<pair<size_t, double>> by_run;
    /// The part of `total` from the panel's vote.
    double panel = 0.0;
    /// The run of the phase-chain site whose panel order that vote keeps, or `NO_RUN` if the panel
    /// did not vote.
    size_t panel_run = NO_RUN;

    /// 1 to flip the site, 0 to keep its order. A tie keeps it.
    int flip() const { return total < 0.0 ? 1 : 0; }

    /// Add the vote `size`, positive for keeping the site's order, to the total and to the part
    /// of run `run`.
    void add(double size, size_t run) {
        total += size;
        if (run == NO_RUN) {
            return;
        }
        auto found = std::ranges::find(by_run, run, &pair<size_t, double>::first);
        if (found == by_run.end()) {
            by_run.emplace_back(run, 0.0);
            found = by_run.end() - 1;
        }
        found->second += size;
    }

    /// Add the panel's vote `size`, for keeping the order relative to a site of run `run`.
    void add_panel(double size, size_t run) {
        total += size;
        panel = size;
        panel_run = run;
    }

    /// Return the run the site's order is relative to: the run whose reads' vote agrees with the
    /// decision and is the largest, the lower run on a tie. If no run's vote agrees, the order came
    /// from the panel, and this is the panel's run if the panel agrees, or else `NO_RUN`.
    size_t deciding_run() const {
        const double sign = flip() ? -1.0 : 1.0;
        size_t best = NO_RUN;
        double best_size = 0.0;
        for (const auto& [run, size] : by_run) {
            const double agreeing = sign * size;
            if (agreeing > best_size || (agreeing == best_size && agreeing > 0.0 && run < best)) {
                best = run;
                best_size = agreeing;
            }
        }
        if (best == NO_RUN && sign * panel > 0.0) {
            best = panel_run;
        }
        return best;
    }
};

/// How a piece of the phase chain is related to the piece before it.
struct PieceJoin {
    /// Whether reads or the panel oriented it against the piece before.
    bool oriented = false;
    /// Whether it is in the run of the piece before: reads link the two, and it is oriented as
    /// their links alone say.
    bool joined = false;
    /// If it starts a run, the run's correlation with the run before (see
    /// `PhaseSite::run_correlation`).
    double correlation = 0.0;
};

/// One read's view of one phase-chain site, for measuring coherence.
struct ChainObservation {
    /// Index of the site in the phase chain.
    size_t chain_index;
    double q0;
    double c;
};

/**
 * A phaser that decides which sites of one phase set to flip against the panel's order.
 *
 * It works through the stages described in read_phasing.hpp. Positions passed to its methods are
 * offsets of sites in the phase set, which is sorted by position. Each site's flip is in `flip`,
 * and holds only once `decided` is set for it.
 */
class PhaseSetPhaser {
public:
    /// Prepare to phase the `site_count` sites starting at `first_site`, into `phase_set_result`.
    PhaseSetPhaser(const PhaseSite* first_site, size_t site_count,
                   const ReadPhasingParams& phasing_params, PhaseSetResult& phase_set_result);

    /// List in the result the pairs of sites that measure the panel against the reads.
    void measure_panel();

    /// Decide every site's flip and run, weighing the panel's order by `errors`.
    void phase(const PanelErrorModel& errors);

private:
    const PhaseSite* first;
    size_t n;
    const ReadPhasingParams& params;
    PhaseSetResult& result;
    /// How often the panel's order is wrong.
    PanelErrorModel panel;

    /// The number of nearest sites on each side that vote on a site off the phase chain.
    size_t per_side;

    /// For each site, 1 to flip it against the panel's order.
    vector<int> flip;
    /// For each site, whether `flip` has been decided.
    vector<char> decided;

    /// The phase chain: the reliable sites not removed for low coherence, in increasing order.
    vector<size_t> chain;
    /// Every other site, in increasing order.
    vector<size_t> unchained;
    /// For each pair of neighbouring phase-chain sites, the size of their link.
    vector<double> links;
    /// Where the phase chain's pieces start, as indexes in `chain`, then `chain.size()`.
    vector<size_t> bounds;
    /// For each piece, how it is related to the piece before it. The first piece's is unused.
    vector<PieceJoin> joins;
    /// The phase-chain sites the panel orders, in increasing order.
    vector<size_t> chain_panel_ordered;
    /// Breaks in the phase chain as the last round of the coherence stage left it.
    size_t breaks = 0;
    /// Of those, the breaks no read links across.
    size_t breaks_no_reads = 0;
    /// Rounds of the coherence stage that removed sites.
    size_t rounds = 0;

    const PhaseSite& site(size_t t) const { return first[t]; }

    /// Build the phase chain and decide its sites' flips, through every round of the coherence
    /// stage.
    void phase_chain();
    /// Measure the links of the phase chain, break it into pieces, and orient each piece's sites
    /// from its first.
    void link_pieces();
    /// Decide each break from the links across it and the panel's order, left to right, flipping
    /// the later piece if they say so.
    void relink_pieces();
    /// Return whether chain sites `before` and `after`, either side of the break before
    /// chain[b_start], are turned against each other relative to the sites either side of it.
    bool turned(size_t before, size_t after, size_t b_start) const;
    /// Return the sum of the links across break `b` between nearby sites of the pieces either side
    /// of it, each referred back to the sites either side of the break.
    double read_relink(size_t b) const;
    /// Return the nearest sites the panel orders either side of break `b`, the one before it at or
    /// after chain[stretch_start] and the one after it in the later piece, if there are both.
    std::optional<std::pair<size_t, size_t>> panel_pair(size_t b, size_t stretch_start) const;
    /// Flip each stretch of pieces that reads or the panel orient against each other as a whole,
    /// so that its first site the panel orders keeps the panel's order.
    void keep_panel_order();
    /// Number the runs of the phase chain as the last round left them, note how far the panel
    /// orients each against the one before, and list the sites that only the reads order.
    void number_runs();
    /// Remove the phase-chain sites whose coherence is too low, and reset every flip. Return
    /// whether any site was removed. `final_round` says this is the last round allowed.
    bool remove_incoherent(bool final_round);
    /// Count, for each phase-chain site, its reads that span another phase-chain site in `total`,
    /// and those that agree with the strand those other sites give them in `agree`.
    void count_agreeing_reads(vector<size_t>& agree, vector<size_t>& total) const;

    /// Hang each site off the phase chain from its nearest phase-chain sites, and return the sites
    /// the panel does not order that no read links to the chain, in increasing order.
    vector<size_t> hang_unchained();
    /// The nearest phase-chain site before or at `t` that the panel orders, or the one after it
    /// if it is nearer. On a tie, the one before.
    size_t nearest_panel_ordered(size_t t) const;

    /// Order by reads the sites in `unlinked`, which neither the panel nor the phase chain
    /// orders.
    void order_unlinked(const vector<size_t>& unlinked);
    /// Add to `vote` on flipping `t` the links to the nearest ordered sites on each side.
    ///
    /// The ordered sites are those in `known`, and those in `done`, which all lie before `t`. Both
    /// are sorted.
    void vote_from_ordered(size_t t, const vector<size_t>& known, const vector<size_t>& done,
                           FlipVote& vote) const;

    /// Add the link between sites `t` and `u` to `vote` on flipping `t`, given `u`'s flip.
    /// Return whether any read made the link.
    bool add_link_vote(size_t t, size_t u, FlipVote& vote) const;
    /// Add to `vote` on flipping `t` the reads' vote from the oriented sites `voters`, one part for
    /// each run among them.
    ///
    /// Each read that `t` shares with a run's voters votes once, however many of them it spans:
    /// its alleles there say which strand it came from, and its allele at `t` then favours `t`'s
    /// order or the flipped one. With one voter, this is the link between the two sites. Each
    /// run whose vote is not 0 counts in `FlipVote::links`.
    void add_read_votes(size_t t, const vector<size_t>& voters, FlipVote& vote) const;
    /// Multiply into `strand0` and `strand1`, for each read `t` shares with the oriented site `u`,
    /// the likelihood of its allele at `u` if it came from strand 0 and if it came from strand 1.
    /// Both are indexed by `t`'s rows.
    void add_strand_likelihoods(size_t t, size_t u, vector<double>& strand0,
                                vector<double>& strand1) const;
};

PhaseSetPhaser::PhaseSetPhaser(const PhaseSite* first_site, size_t site_count,
                               const ReadPhasingParams& phasing_params,
                               PhaseSetResult& phase_set_result)
    : first(first_site), n(site_count), params(phasing_params), result(phase_set_result),
      per_side(phasing_params.hang / 2 + 1), flip(site_count, 0), decided(site_count, 0) {
    result.run.assign(site_count, NO_RUN);
    // A site the panel does not order may join the phase chain too. Its reads order it against the
    // rest of the chain, and `keep_panel_order` takes the panel's order only from sites the panel
    // orders.
    for (size_t t = 0; t < n; ++t) {
        if (site(t).reliability >= params.reliability) {
            chain.push_back(t);
        } else {
            unchained.push_back(t);
        }
    }
}

void PhaseSetPhaser::measure_panel() {
    if (chain.empty()) {
        // No site is reliable enough to measure the panel against.
        return;
    }
    // Reads alone order the sites of a piece against each other.
    link_pieces();
    for (size_t b = 0; b + 1 < bounds.size(); ++b) {
        size_t previous = n;
        for (size_t m = bounds[b]; m < bounds[b + 1]; ++m) {
            const size_t t = chain[m];
            if (!site(t).panel_ordered) {
                continue;
            }
            if (previous != n) {
                result.panel_pairs.emplace_back(site(t).position - site(previous).position,
                                                flip[t] != flip[previous]);
            }
            previous = t;
        }
    }
}

void PhaseSetPhaser::phase(const PanelErrorModel& errors) {
    panel = errors;
    ++result.counters.chains;
    result.counters.reliable += chain.size();

    phase_chain();
    number_runs();
    result.counters.breaks += breaks;
    result.counters.breaks_no_reads += breaks_no_reads;
    result.counters.coherence_rounds_run = max(result.counters.coherence_rounds_run, rounds);

    vector<size_t> unlinked;
    if (chain.empty()) {
        for (size_t t : unchained) {
            if (!site(t).panel_ordered) {
                unlinked.push_back(t);
            }
        }
    } else {
        unlinked = hang_unchained();
    }
    if (!unlinked.empty()) {
        order_unlinked(unlinked);
    }

    for (size_t t = 0; t < n; ++t) {
        if (decided[t] && flip[t]) {
            result.flipped.push_back(site(t).record_key);
            ++result.counters.flipped;
        }
    }
}

void PhaseSetPhaser::phase_chain() {
    // One round with the coherence stage off, and otherwise one round more than the coherence
    // stage may remove sites in, so the last removal is followed by a fresh chain.
    const size_t max_pass = params.coherence_min > 0.0
                                ? max<size_t>(1, params.coherence_rounds) + 1
                                : 1;
    for (size_t pass = 0; pass < max_pass && !chain.empty(); ++pass) {
        link_pieces();
        relink_pieces();
        keep_panel_order();
        // The coherence stage is skipped in the last round, so that the chain it leaves is
        // linked again before anything hangs from it.
        if (pass + 1 >= max_pass || !remove_incoherent(pass + 2 == max_pass)) {
            break;
        }
    }
}

void PhaseSetPhaser::link_pieces() {
    links.assign(chain.empty() ? 0 : chain.size() - 1, 0.0);
    for (size_t m = 0; m + 1 < chain.size(); ++m) {
        links[m] = phase_link(site(chain[m]), site(chain[m + 1]), params.cap);
    }
    bounds.clear();
    bounds.push_back(0);
    for (size_t m = 0; m + 1 < chain.size(); ++m) {
        // A link no read makes is a break whatever the threshold, so that reads link the sites of
        // every piece.
        if (links[m] == 0.0 || std::fabs(links[m]) < params.break_threshold) {
            bounds.push_back(m + 1);
        }
    }
    bounds.push_back(chain.size());
    joins.assign(bounds.size() - 1, PieceJoin {});
    breaks = bounds.size() - 2;
    breaks_no_reads = 0;
    chain_panel_ordered.clear();
    std::ranges::copy_if(chain, std::back_inserter(chain_panel_ordered),
                         [&](size_t t) { return site(t).panel_ordered; });

    for (size_t b = 0; b + 1 < bounds.size(); ++b) {
        const size_t piece_start = bounds[b];
        const size_t piece_end = bounds[b + 1];
        flip[chain[piece_start]] = 0;
        decided[chain[piece_start]] = 1;
        for (size_t m = piece_start; m + 1 < piece_end; ++m) {
            flip[chain[m + 1]] = flip[chain[m]] ^ (links[m] < 0.0 ? 1 : 0);
            decided[chain[m + 1]] = 1;
        }
    }
}

bool PhaseSetPhaser::turned(size_t before, size_t after, size_t b_start) const {
    return (flip[before] ^ flip[chain[b_start - 1]]) ^ (flip[after] ^ flip[chain[b_start]]);
}

double PhaseSetPhaser::read_relink(size_t b) const {
    // Piece A is chain[bounds[b], b_start) and piece B chain[b_start, b_end).
    const size_t b_start = bounds[b + 1];
    const size_t b_end = bounds[b + 2];
    const size_t a_first = b_start > params.relink ? b_start - params.relink : 0;
    const size_t a_lo = max(bounds[b], a_first);
    const size_t b_hi = min(b_end, b_start + params.relink);
    double total = 0.0;
    for (size_t i = a_lo; i < b_start; ++i) {
        for (size_t j = b_start; j < b_hi; ++j) {
            const double dv = phase_link(site(chain[i]), site(chain[j]), params.cap);
            total += turned(chain[i], chain[j], b_start) ? -dv : dv;
        }
    }
    return total;
}

std::optional<std::pair<size_t, size_t>> PhaseSetPhaser::panel_pair(size_t b,
                                                                    size_t stretch_start) const {
    const auto after = std::ranges::lower_bound(chain_panel_ordered, chain[bounds[b + 1]]);
    if (after == chain_panel_ordered.begin() || after == chain_panel_ordered.end()
        || *(after - 1) < chain[stretch_start] || *after > chain[bounds[b + 2] - 1]) {
        return std::nullopt;
    }
    return std::pair {*(after - 1), *after};
}

void PhaseSetPhaser::relink_pieces() {
    // The question is the relative orientation of the last site before each break and the first
    // after it. The pieces from chain[stretch_start] to the break are oriented against each
    // other, so a site among them relates to the last one.
    size_t stretch_start = 0;
    for (size_t b = 0; b + 2 < bounds.size(); ++b) {
        const size_t b_start = bounds[b + 1];
        const double reads = read_relink(b);
        breaks_no_reads += reads == 0.0 ? 1 : 0;
        double total = reads;
        // The panel votes for the order it gives its nearest sites either side of the break.
        if (const auto pair = panel.fitted ? panel_pair(b, stretch_start) : std::nullopt) {
            const auto [before, after] = *pair;
            const size_t distance = site(after).position - site(before).position;
            total += panel.vote(distance) * (turned(before, after, b_start) ? -1.0 : 1.0);
        }
        if (total == 0.0) {
            // Neither reads nor the panel relate the later piece to the earlier, so
            // `keep_panel_order` orients it alone.
            stretch_start = b_start;
            continue;
        }
        PieceJoin& join = joins[b + 1];
        join.oriented = true;
        join.joined = reads != 0.0 && (reads < 0.0) == (total < 0.0);
        if (!join.joined) {
            // The decision is right with log10 odds |total|, which at a break no read links are the
            // panel's vote, and then this is the panel's correlation.
            join.correlation = 2.0 / (1.0 + std::pow(10.0, -std::fabs(total))) - 1.0;
        }
        // The later piece is flipped if its first site and the earlier piece's last are not
        // oriented against each other as the sum says.
        const bool apart = flip[chain[b_start - 1]] != flip[chain[b_start]];
        if ((total < 0.0) != apart) {
            for (size_t m = b_start; m < bounds[b + 2]; ++m) {
                flip[chain[m]] ^= 1;
            }
        }
    }
}

void PhaseSetPhaser::keep_panel_order() {
    // Reads, or the panel, orient each piece against the one before it, except where neither
    // says anything, and a stretch of pieces so oriented is turned as a whole. A stretch with no
    // site the panel orders keeps the order of its first site, which is arbitrary.
    for (size_t b = 0; b + 1 < bounds.size();) {
        size_t e = b + 1;
        while (e + 1 < bounds.size() && joins[e].oriented) {
            ++e;
        }
        // The stretch is chain[start, end).
        const size_t start = bounds[b];
        const size_t end = bounds[e];
        if (const auto framing = std::ranges::find_if(
                chain.begin() + start, chain.begin() + end,
                [this](size_t t) { return site(t).panel_ordered; });
            framing != chain.begin() + end && flip[*framing]) {
            for (size_t m = start; m < end; ++m) {
                flip[chain[m]] ^= 1;
            }
        }
        b = e;
    }
}

void PhaseSetPhaser::number_runs() {
    // Piece b starts a run unless it is joined to the piece before.
    result.run_correlation.clear();
    vector<size_t> run_sites;
    for (size_t b = 0; b + 1 < bounds.size(); ++b) {
        if (b == 0 || !joins[b].joined) {
            result.run_correlation.push_back(joins[b].correlation);
            run_sites.push_back(0);
        }
        for (size_t m = bounds[b]; m < bounds[b + 1]; ++m) {
            result.run[chain[m]] = result.run_correlation.size() - 1;
            ++run_sites.back();
        }
    }
    // A phase-chain site the panel does not order is ordered by reads if its run has another
    // site, since reads make every link inside a piece and every relink.
    for (size_t t : chain) {
        if (!site(t).panel_ordered && run_sites[result.run[t]] >= 2) {
            result.read_ordered.push_back(site(t).record_key);
        }
    }
}

void PhaseSetPhaser::count_agreeing_reads(vector<size_t>& agree, vector<size_t>& total) const {
    agree.assign(chain.size(), 0);
    total.assign(chain.size(), 0);
    unordered_map<uint64_t, vector<ChainObservation>> by_read;
    for (size_t m = 0; m < chain.size(); ++m) {
        const PhaseSite& st = site(chain[m]);
        for (size_t i = 0; i < st.read_key.size(); ++i) {
            by_read[st.read_key[i]].emplace_back(m, (double)st.q0[i], (double)st.c[i]);
        }
    }
    for (const auto& [read, obs] : by_read) {
        if (obs.size() < 2) {
            // No other site to hold this one out against.
            continue;
        }
        // Log10 likelihoods of the read coming from strand 0 and strand 1, from each site and in
        // total, in the current flips.
        double l0 = 0.0;
        double l1 = 0.0;
        vector<double> t0(obs.size());
        vector<double> t1(obs.size());
        for (size_t k = 0; k < obs.size(); ++k) {
            const double p = obs[k].c;
            const double q = flip[chain[obs[k].chain_index]] ? 1.0 - obs[k].q0 : obs[k].q0;
            t0[k] = log10(strand_likelihood(q, p));
            t1[k] = log10(strand_likelihood(1.0 - q, p));
            l0 += t0[k];
            l1 += t1[k];
        }
        for (size_t k = 0; k < obs.size(); ++k) {
            const size_t m = obs[k].chain_index;
            // The strand the read's other sites give it, leaving this site out.
            const bool hap1 = (l1 - t1[k]) > (l0 - t0[k]);
            const double q = flip[chain[m]] ? 1.0 - obs[k].q0 : obs[k].q0;
            const bool says1 = q < 0.5;
            ++total[m];
            agree[m] += says1 == hap1 ? 1 : 0;
        }
    }
}

bool PhaseSetPhaser::remove_incoherent(bool final_round) {
    if (params.coherence_min <= 0.0 || chain.size() < 3) {
        return false;
    }
    vector<size_t> agree;
    vector<size_t> total;
    count_agreeing_reads(agree, total);
    vector<size_t> keep;
    vector<size_t> drop;
    for (size_t m = 0; m < chain.size(); ++m) {
        // A site needs at least ten counted reads before low coherence can remove it.
        if (total[m] >= 10 && (double)agree[m] / (double)total[m] < params.coherence_min) {
            drop.push_back(chain[m]);
        } else {
            keep.push_back(chain[m]);
        }
    }
    if (drop.empty() || keep.size() < 2) {
        return false;
    }
    result.counters.demoted_incoherent += drop.size();
    ++rounds;
    if (final_round) {
        // Still removing sites in the last allowed round, so the chain has not converged.
        ++result.counters.coherence_unconverged;
    }
    unchained.insert(unchained.end(), drop.begin(), drop.end());
    std::ranges::sort(unchained);
    chain.swap(keep);
    std::ranges::fill(flip, 0);
    std::ranges::fill(decided, 0);
    return true;
}

bool PhaseSetPhaser::add_link_vote(size_t t, size_t u, FlipVote& vote) const {
    const double dv = phase_link(site(t), site(u), params.cap);
    if (dv == 0.0) {
        return false;
    }
    const int pred = flip[u] ^ (dv < 0.0 ? 1 : 0);
    vote.add(std::fabs(dv) * (pred == 0 ? 1.0 : -1.0), result.run[u]);
    ++vote.links;
    return true;
}

void PhaseSetPhaser::add_strand_likelihoods(size_t t, size_t u, vector<double>& strand0,
                                           vector<double>& strand1) const {
    const PhaseSite& target = site(t);
    const PhaseSite& voter = site(u);
    // Both read lists are sorted by key, one row per key, so we merge them.
    size_t i = 0;
    size_t j = 0;
    while (i < target.read_key.size() && j < voter.read_key.size()) {
        if (target.read_key[i] < voter.read_key[j]) {
            ++i;
        } else if (voter.read_key[j] < target.read_key[i]) {
            ++j;
        } else {
            const double c = voter.c[j];
            const double q = flip[u] ? 1.0 - voter.q0[j] : voter.q0[j];
            strand0[i] *= strand_likelihood(q, c);
            strand1[i] *= strand_likelihood(1.0 - q, c);
            ++i;
            ++j;
        }
    }
}

/// Return the vote of `target`'s reads, in log10 units, on keeping its order, from each read's
/// likelihoods `strand0` and `strand1` of having come from strand 0 and strand 1, limited to `cap`
/// if that is not 0.
double strand_vote(const PhaseSite& target, const vector<double>& strand0,
                   const vector<double>& strand1, double cap) {
    double total = 0.0;
    for (size_t i = 0; i < target.read_key.size(); ++i) {
        // The likelihoods of the read's allele at the site on the strand of slot 0 and on the
        // strand of slot 1.
        const double on0 = strand_likelihood(target.q0[i], target.c[i]);
        const double on1 = strand_likelihood(1.0 - target.q0[i], target.c[i]);
        // A read that shares no voting site has equal strand likelihoods and adds nothing.
        const double keep = strand0[i] * on0 + strand1[i] * on1;
        const double flipped = strand0[i] * on1 + strand1[i] * on0;
        if (keep > 0.0 && flipped > 0.0) {
            total += log10(keep / flipped);
        }
    }
    return cap > 0.0 ? max(-cap, min(cap, total)) : total;
}

void PhaseSetPhaser::add_read_votes(size_t t, const vector<size_t>& voters, FlipVote& vote) const {
    vector<size_t> runs;
    for (size_t u : voters) {
        if (std::ranges::find(runs, result.run[u]) == runs.end()) {
            runs.push_back(result.run[u]);
        }
    }
    // For each of `t`'s reads, the likelihoods of its alleles at one run's voters if it came from
    // strand 0 and if it came from strand 1.
    vector<double> strand0(site(t).read_key.size());
    vector<double> strand1(site(t).read_key.size());
    for (size_t run : runs) {
        std::ranges::fill(strand0, 1.0);
        std::ranges::fill(strand1, 1.0);
        for (size_t u : voters) {
            if (result.run[u] == run) {
                add_strand_likelihoods(t, u, strand0, strand1);
            }
        }
        if (const double total = strand_vote(site(t), strand0, strand1, params.cap); total != 0.0) {
            vote.add(total, run);
            ++vote.links;
        }
    }
}

size_t PhaseSetPhaser::nearest_panel_ordered(size_t t) const {
    const auto after = std::ranges::lower_bound(chain_panel_ordered, t);
    if (after == chain_panel_ordered.begin()) {
        return chain_panel_ordered.front();
    }
    if (after == chain_panel_ordered.end()) {
        return chain_panel_ordered.back();
    }
    const size_t before = *(after - 1);
    return t - before <= *after - t ? before : *after;
}

vector<size_t> PhaseSetPhaser::hang_unchained() {
    vector<size_t> unlinked;
    for (size_t t : unchained) {
        // The phase-chain sites before `t` are chain[0, split) and those after it
        // chain[split, end), since `chain` is sorted and does not hold `t`. The nearest on each
        // side vote, so that a site near an end of the chain still hears from both sides.
        const size_t split = std::ranges::lower_bound(chain, t) - chain.begin();
        const size_t lo = split > per_side ? split - per_side : 0;
        const vector<size_t> voters(chain.begin() + lo,
                                    chain.begin() + min(chain.size(), split + per_side));
        FlipVote vote;
        add_read_votes(t, voters, vote);
        const bool panel_ordered = site(t).panel_ordered;
        if (!chain_panel_ordered.empty() && panel_ordered) {
            // The panel votes for keeping this site's order relative to the nearest phase-chain
            // site the panel also orders, as strongly as its order of sites that far apart is
            // right.
            const size_t nearest = nearest_panel_ordered(t);
            const size_t distance = max(site(t).position, site(nearest).position)
                                    - min(site(t).position, site(nearest).position);
            const double weight = panel.fitted ? panel.vote(distance) : params.panel_weight;
            vote.add_panel(weight * (flip[nearest] == 0 ? 1.0 : -1.0), result.run[nearest]);
        }
        result.run[t] = vote.deciding_run();
        if (vote.links == 0 && !panel_ordered) {
            // `order_unlinked` counts it if no read orders it there either.
            unlinked.push_back(t);
        } else if (vote.links == 0) {
            ++result.counters.hung_no_reads;
        } else if (!panel_ordered) {
            result.read_ordered.push_back(site(t).record_key);
        }
        flip[t] = vote.flip();
        decided[t] = 1;
        ++result.counters.hung;
    }
    return unlinked;
}

void PhaseSetPhaser::order_unlinked(const vector<size_t>& unlinked) {
    // Such a site, as where no panel haplotype reaches, may still share reads with a site whose
    // order is known: a phase-chain site beyond the nearest few, a hung site, or one ordered
    // earlier here. We link each, in increasing order, to the nearest such sites on each side.
    // Sites that share no read with any of them are then ordered against each other, each
    // against the nearest before it, so that reads order each run of them among themselves.
    vector<char> is_unlinked(n, 0);
    for (size_t t : unlinked) {
        is_unlinked[t] = 1;
    }
    // Every other site is ordered, by the phase chain, by the panel, or by being hung.
    vector<size_t> known;
    known.reserve(n - unlinked.size());
    for (size_t t = 0; t < n; ++t) {
        if (!is_unlinked[t]) {
            known.push_back(t);
        }
    }
    vector<char> ordered(n, 0);

    // Order against the known sites. `done` holds the sites ordered here so far, which all lie
    // before `t`, so the nearest sites before `t` are merged from it and `known`.
    vector<size_t> done;
    vector<size_t> left;
    for (size_t t : unlinked) {
        FlipVote vote;
        vote_from_ordered(t, known, done, vote);
        if (vote.links == 0) {
            left.push_back(t);
            continue;
        }
        flip[t] = vote.flip();
        decided[t] = 1;
        ordered[t] = 1;
        result.run[t] = vote.deciding_run();
        done.push_back(t);
    }

    // Order the rest against each other. The first site of each run keeps its arbitrary order,
    // and a site with no read to any other site stays unordered.
    for (size_t i = 0; i < left.size(); ++i) {
        const size_t t = left[i];
        FlipVote vote;
        for (size_t k = 0; k < per_side && k < i; ++k) {
            const size_t u = left[i - 1 - k];
            if (add_link_vote(t, u, vote)) {
                ordered[u] = 1;
            }
        }
        flip[t] = vote.flip();
        decided[t] = 1;
        ordered[t] = ordered[t] || vote.links > 0;
        // Such sites are ordered only among themselves, so each set that reads join is a run, and
        // no panel order frames it.
        result.run[t] = vote.deciding_run();
        if (result.run[t] == NO_RUN) {
            result.run[t] = result.run_correlation.size();
            result.run_correlation.push_back(0.0);
        }
    }

    for (size_t t : unlinked) {
        if (ordered[t]) {
            result.read_ordered.push_back(site(t).record_key);
        } else if (!chain.empty()) {
            ++result.counters.hung_no_reads;
        }
    }
}

void PhaseSetPhaser::vote_from_ordered(size_t t, const vector<size_t>& known,
                                       const vector<size_t>& done, FlipVote& vote) const {
    // Merge the nearest sites before `t` from `known` and `done`, nearest first.
    const size_t split = std::ranges::lower_bound(known, t) - known.begin();
    size_t a = split;
    size_t b = done.size();
    size_t k = 0;
    while (k < per_side && (a > 0 || b > 0)) {
        if (b == 0 || (a > 0 && known[a - 1] > done[b - 1])) {
            --a;
            add_link_vote(t, known[a], vote);
        } else {
            --b;
            add_link_vote(t, done[b], vote);
        }
        ++k;
    }
    for (size_t j = split; j < known.size() && j - split < per_side; ++j) {
        add_link_vote(t, known[j], vote);
    }
}

/// Fit the panel's error rates to the pairs of `PhaseSetResult::panel_pairs`, into `counters`.
///
/// The rates are the maximum-likelihood ones under `PanelErrorModel`, on a grid: switch rates
/// from one per 10^9 bases to one per 100, and flip rates from 10^-4 to 10^-1, a floor that keeps
/// the panel's vote below 4 whatever the distance. With too few pairs, nothing is fitted. The
/// fit's table bins the pairs by distance, each bin holding an equal share.
void fit_panel_errors(vector<pair<size_t, bool>>& pairs, ReadPhasingCounters& counters) {
    counters.panel_pairs = pairs.size();
    counters.panel_pairs_reversed = std::ranges::count(pairs, true, &pair<size_t, bool>::second);
    if (pairs.size() < PANEL_FIT_MIN_PAIRS) {
        return;
    }
    double best = -std::numeric_limits<double>::infinity();
    vector<double> decay(pairs.size());
    for (int s = 0; s <= 140; ++s) {
        const double switch_rate = std::pow(10.0, -9.0 + s / 20.0);
        for (size_t i = 0; i < pairs.size(); ++i) {
            decay[i] = std::exp(-2.0 * switch_rate * (double)pairs[i].first);
        }
        for (int f = 0; f <= 30; ++f) {
            const double flip_rate = std::pow(10.0, -4.0 + f / 10.0);
            const double kept = (1.0 - 2.0 * flip_rate) * (1.0 - 2.0 * flip_rate);
            double total = 0.0;
            for (size_t i = 0; i < pairs.size(); ++i) {
                const double c = kept * decay[i];
                total += std::log((pairs[i].second ? 1.0 - c : 1.0 + c) / 2.0);
            }
            if (total > best) {
                best = total;
                counters.panel = PanelErrorModel{switch_rate, flip_rate, true};
            }
        }
    }

    std::ranges::sort(pairs);
    for (size_t bin = 0; bin < PANEL_FIT_BINS; ++bin) {
        const size_t begin = pairs.size() * bin / PANEL_FIT_BINS;
        const size_t end = pairs.size() * (bin + 1) / PANEL_FIT_BINS;
        double distance = 0.0;
        double agree = 0.0;
        double predicted = 0.0;
        for (size_t i = begin; i < end; ++i) {
            distance += (double)pairs[i].first;
            agree += pairs[i].second ? 0.0 : 1.0;
            predicted += (1.0 + counters.panel.correlation(pairs[i].first)) / 2.0;
        }
        const auto count = (double)(end - begin);
        counters.panel_fit.emplace_back(distance / count, agree / count, predicted / count,
                                        end - begin);
    }
}

/// Return the [begin, end) ranges of `sites` that hold one phase set each. `sites` must be sorted
/// by phase set.
vector<pair<size_t, size_t>> phase_set_ranges(const vector<PhaseSite>& sites) {
    vector<pair<size_t, size_t>> ranges;
    for (size_t begin = 0; begin < sites.size();) {
        size_t end = begin;
        while (end < sites.size() && sites[end].phase_set == sites[begin].phase_set) {
            ++end;
        }
        ranges.emplace_back(begin, end);
        begin = end;
    }
    return ranges;
}

}

unordered_set<size_t> read_phase_flips(vector<PhaseSite>& sites, const ReadPhasingParams& params,
                                       ReadPhasingCounters& counters,
                                       unordered_set<size_t>* read_ordered) {
    unordered_set<size_t> flips;
    if (sites.empty()) {
        return flips;
    }
    // `record_key` breaks ties between sites at the same position, so that the order does not
    // depend on the order in which the sites arrived.
    std::ranges::sort(sites, [](const PhaseSite& x, const PhaseSite& y) {
        if (x.phase_set != y.phase_set) {
            return x.phase_set < y.phase_set;
        }
        if (x.position != y.position) {
            return x.position < y.position;
        }
        return x.record_key < y.record_key;
    });
    // `phase_link` merges two sites' rows, so it needs them sorted, and coherence needs each read
    // once per site.
#pragma omp parallel for schedule(dynamic, 1024)
    for (PhaseSite& site : sites) {
        merge_mates(site);
    }
    counters.sites += sites.size();

    // Phase sets are independent, so we phase them in parallel, longest first. Each keeps its own
    // counters and lists, which we gather in phase-set order afterwards.
    const vector<pair<size_t, size_t>> phase_sets = phase_set_ranges(sites);
    vector<size_t> by_length(phase_sets.size());
    std::iota(by_length.begin(), by_length.end(), 0);
    std::ranges::stable_sort(by_length, [&](size_t x, size_t y) {
        return phase_sets[x].second - phase_sets[x].first > phase_sets[y].second - phase_sets[y].first;
    });
    // The panel's error rates are fitted from all phase sets before any is phased, since the
    // phasing weighs the panel's order by them.
    vector<PhaseSetResult> results(phase_sets.size());
#pragma omp parallel for schedule(dynamic, 1)
    for (size_t set_index : by_length) {
        const auto [begin, end] = phase_sets[set_index];
        PhaseSetPhaser(sites.data() + begin, end - begin, params, results[set_index])
            .measure_panel();
    }
    vector<pair<size_t, bool>> panel_pairs;
    for (const PhaseSetResult& result : results) {
        panel_pairs.insert(panel_pairs.end(), result.panel_pairs.begin(), result.panel_pairs.end());
    }
    fit_panel_errors(panel_pairs, counters);
#pragma omp parallel for schedule(dynamic, 1)
    for (size_t set_index : by_length) {
        const auto [begin, end] = phase_sets[set_index];
        PhaseSetPhaser(sites.data() + begin, end - begin, params, results[set_index])
            .phase(counters.panel);
    }

    for (size_t i = 0; i < phase_sets.size(); ++i) {
        const PhaseSetResult& result = results[i];
        for (size_t t = 0; t < result.run.size(); ++t) {
            PhaseSite& site = sites[phase_sets[i].first + t];
            site.run = result.run[t];
            site.run_correlation = site.run == NO_RUN ? 0.0 : result.run_correlation[site.run];
        }
        add_phase_set_counters(result.counters, counters);
        // Insert one at a time, as the set's iteration order can depend on how it grew.
        for (size_t key : result.flipped) {
            flips.insert(key);
        }
        if (read_ordered != nullptr) {
            read_ordered->insert(result.read_ordered.begin(), result.read_ordered.end());
        }
    }
    return flips;
}

vector<double> allele_length_weights(const vector<std::uint32_t>& allele_length, size_t n_alleles,
                                     float mean_read_length, bool length_weighted,
                                     const vector<int>& slot_allele) {
    const size_t n_slots = slot_allele.size();
    vector weight(n_slots, n_slots ? 1.0 / (double)n_slots : 0.0);
    if (length_weighted && mean_read_length > 0.0
        && allele_length.size() == n_alleles && n_slots > 1) {
        double total = 0.0;
        vector raw(n_slots, 0.0);
        for (size_t i = 0; i < n_slots; ++i) {
            raw[i] = (double)allele_length[slot_allele[i]]
                     + (double)mean_read_length - 1.0;
            total += raw[i];
        }
        if (total > 0.0) {
            for (size_t i = 0; i < n_slots; ++i) {
                weight[i] = raw[i] / total;
            }
        }
    }
    return weight;
}

}
