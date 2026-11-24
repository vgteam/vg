#include "read_phasing.hpp"

#include <algorithm>
#include <cmath>
#include <numeric>
#include <unordered_map>

namespace vg {

using std::log10;
using std::max;
using std::min;
using std::pair;
using std::unordered_map;

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

/// A tally of the links that vote on whether to flip one site.
struct FlipVote {
    /// Sum of the link sizes, in log10 units, positive for keeping the site's order and negative
    /// for flipping it.
    double total = 0.0;
    /// How many links had any read.
    size_t links = 0;

    /// 1 to flip the site, 0 to keep its order. A tie keeps it.
    int flip() const { return total < 0.0 ? 1 : 0; }
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
    /// Prepare to phase the `site_count` sites starting at `first_site`, adding what happens to
    /// `set_counters`, and the record keys of the sites only reads order to `set_read_ordered`.
    PhaseSetPhaser(const PhaseSite* first_site, size_t site_count,
                   const ReadPhasingParams& phasing_params, ReadPhasingCounters& set_counters,
                   vector<size_t>& set_read_ordered);

    /// Decide every site's flip, and append the record keys of the flipped sites to `flipped`.
    void phase(vector<size_t>& flipped);

private:
    const PhaseSite* first;
    size_t n;
    const ReadPhasingParams& params;
    ReadPhasingCounters& counters;
    /// Record keys of the sites the panel does not order whose reads ordered them.
    vector<size_t>& read_ordered;

    /// Nearest sites on each side that hanging a site links it to.
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
    /// For each piece, whether reads relinked it to the piece before it.
    vector<char> joined;
    /// The phase-chain sites the panel orders, in increasing order.
    vector<size_t> chain_panel_ordered;
    /// Record keys of the phase-chain sites the panel does not order whose reads ordered them.
    vector<size_t> chain_read_ordered;
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
    /// Decide each break from the links across it, left to right, flipping the later piece if
    /// they say so.
    void relink_pieces();
    /// Flip each run of pieces as a whole so that its first site the panel orders keeps the
    /// panel's order.
    void frame_runs();
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
};

PhaseSetPhaser::PhaseSetPhaser(const PhaseSite* first_site, size_t site_count,
                               const ReadPhasingParams& phasing_params,
                               ReadPhasingCounters& set_counters,
                               vector<size_t>& set_read_ordered)
    : first(first_site), n(site_count), params(phasing_params), counters(set_counters),
      read_ordered(set_read_ordered), per_side(phasing_params.hang / 2 + 1), flip(site_count, 0),
      decided(site_count, 0) {
    // A site the panel does not order may join the phase chain too. Its reads order it against
    // the rest of the chain, and `frame_runs` takes the panel's order only from sites the panel
    // orders.
    for (size_t t = 0; t < n; ++t) {
        if (site(t).reliability >= params.reliability) {
            chain.push_back(t);
        } else {
            unchained.push_back(t);
        }
    }
}

void PhaseSetPhaser::phase(vector<size_t>& flipped) {
    ++counters.chains;
    counters.reliable += chain.size();

    phase_chain();
    counters.breaks += breaks;
    counters.breaks_no_reads += breaks_no_reads;
    counters.coherence_rounds_run = max(counters.coherence_rounds_run, rounds);

    // A phase-chain site the panel does not order is ordered by reads if its run has another
    // site, since reads make every link inside a piece and every relink.
    read_ordered.insert(read_ordered.end(), chain_read_ordered.begin(), chain_read_ordered.end());

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
            flipped.push_back(site(t).record_key);
            ++counters.flipped;
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
        frame_runs();
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
    joined.assign(bounds.size() - 1, 0);
    breaks = bounds.size() - 2;
    breaks_no_reads = 0;

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

void PhaseSetPhaser::relink_pieces() {
    // Piece A is chain[a_start, b_start) and piece B chain[b_start, b_end).
    for (size_t b = 0; b + 2 < bounds.size(); ++b) {
        const size_t b_start = bounds[b + 1];
        const size_t b_end = bounds[b + 2];
        const size_t a_first = b_start > params.relink ? b_start - params.relink : 0;
        const size_t a_lo = max(bounds[b], a_first);
        const size_t b_hi = min(b_end, b_start + params.relink);
        // The question is the relative orientation of A's last site and B's first. Each link is
        // referred back to those two sites through the flips inside each piece.
        double total = 0.0;
        for (size_t i = a_lo; i < b_start; ++i) {
            for (size_t j = b_start; j < b_hi; ++j) {
                const double dv = phase_link(site(chain[i]), site(chain[j]), params.cap);
                if (dv == 0.0) {
                    continue;
                }
                const int turned = (flip[chain[i]] ^ flip[chain[b_start - 1]])
                                   ^ (flip[chain[j]] ^ flip[chain[b_start]]);
                total += turned ? -dv : dv;
            }
        }
        if (total == 0.0) {
            // `frame_runs` orients B from the panel instead.
            ++breaks_no_reads;
            continue;
        }
        joined[b + 1] = 1;
        const int x = total < 0.0 ? 1 : 0;
        if (x ^ flip[chain[b_start - 1]] ^ flip[chain[b_start]]) {
            for (size_t m = b_start; m < b_end; ++m) {
                flip[chain[m]] ^= 1;
            }
        }
    }
}

void PhaseSetPhaser::frame_runs() {
    // Reads fix the order of a run's sites relative to each other, and only the panel relates the
    // run to the runs either side. So we turn the run as a whole before the coherence stage
    // compares sites across runs. A run with no site the panel orders keeps the order of its first
    // site, which is arbitrary.
    chain_panel_ordered.clear();
    chain_read_ordered.clear();
    for (size_t b = 0; b + 1 < bounds.size();) {
        size_t e = b + 1;
        while (e + 1 < bounds.size() && joined[e]) {
            ++e;
        }
        // The run is chain[run_start, run_end).
        const size_t run_start = bounds[b];
        const size_t run_end = bounds[e];
        size_t framing = run_end;
        for (size_t m = run_start; m < run_end; ++m) {
            if (site(chain[m]).panel_ordered) {
                framing = min(framing, m);
                chain_panel_ordered.push_back(chain[m]);
            } else if (run_end - run_start >= 2) {
                chain_read_ordered.push_back(site(chain[m]).record_key);
            }
        }
        if (framing < run_end && flip[chain[framing]]) {
            for (size_t m = run_start; m < run_end; ++m) {
                flip[chain[m]] ^= 1;
            }
        }
        b = e;
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
            t0[k] = log10(p * q + (1.0 - p) * 0.5);
            t1[k] = log10(p * (1.0 - q) + (1.0 - p) * 0.5);
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
    counters.demoted_incoherent += drop.size();
    ++rounds;
    if (final_round) {
        // Still removing sites in the last allowed round, so the chain has not converged.
        ++counters.coherence_unconverged;
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
    vote.total += std::fabs(dv) * (pred == 0 ? 1.0 : -1.0);
    ++vote.links;
    return true;
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
        FlipVote vote;
        // The phase-chain sites before `t` are chain[0, split) and those after it
        // chain[split, end), since `chain` is sorted and does not hold `t`.
        const size_t split = std::ranges::lower_bound(chain, t) - chain.begin();
        // Link to the nearest phase-chain sites on both sides, nearest first, before `t` and then
        // after it. A site near an end of the chain thus still hears from both sides.
        for (size_t k = 0; k < per_side && k < split; ++k) {
            add_link_vote(t, chain[split - 1 - k], vote);
        }
        for (size_t k = split; k < chain.size() && k - split < per_side; ++k) {
            add_link_vote(t, chain[k], vote);
        }
        const bool panel_ordered = site(t).panel_ordered;
        if (params.panel_weight > 0.0 && !chain_panel_ordered.empty() && panel_ordered) {
            // The panel votes for keeping this site's order relative to the nearest phase-chain
            // site the panel also orders. The vote is weak, and decides only when reads say
            // nothing.
            vote.total += params.panel_weight * (flip[nearest_panel_ordered(t)] == 0 ? 1.0 : -1.0);
        }
        if (vote.links == 0 && !panel_ordered) {
            // `order_unlinked` counts it if no read orders it there either.
            unlinked.push_back(t);
        } else if (vote.links == 0) {
            ++counters.hung_no_reads;
        } else if (!panel_ordered) {
            read_ordered.push_back(site(t).record_key);
        }
        flip[t] = vote.flip();
        decided[t] = 1;
        ++counters.hung;
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
    }

    for (size_t t : unlinked) {
        if (ordered[t]) {
            read_ordered.push_back(site(t).record_key);
        } else if (!chain.empty()) {
            ++counters.hung_no_reads;
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
    vector<ReadPhasingCounters> set_counters(phase_sets.size());
    vector<vector<size_t>> set_flips(phase_sets.size());
    vector<vector<size_t>> set_read_ordered(phase_sets.size());

#pragma omp parallel for schedule(dynamic, 1)
    for (size_t set_index : by_length) {
        const auto [begin, end] = phase_sets[set_index];
        PhaseSetPhaser phaser(sites.data() + begin, end - begin, params, set_counters[set_index],
                              set_read_ordered[set_index]);
        phaser.phase(set_flips[set_index]);
    }

    for (size_t i = 0; i < phase_sets.size(); ++i) {
        add_phase_set_counters(set_counters[i], counters);
        // Insert one at a time, as the set's iteration order can depend on how it grew.
        for (size_t key : set_flips[i]) {
            flips.insert(key);
        }
        if (read_ordered != nullptr) {
            read_ordered->insert(set_read_ordered[i].begin(), set_read_ordered[i].end());
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
