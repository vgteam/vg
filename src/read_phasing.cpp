#include "read_phasing.hpp"

#include <algorithm>
#include <cmath>
#include <array>
#include <unordered_map>

namespace vg {

using std::log10;
using std::max;
using std::min;
using std::sort;
using std::unordered_map;

double phase_link(const PhaseSite& a, const PhaseSite& b, double cap) {
    // Both read lists are sorted by key, so this is a merge rather than a lookup per read.
    double total = 0.0;
    size_t i = 0, j = 0;
    while (i < a.read_key.size() && j < b.read_key.size()) {
        if (a.read_key[i] < b.read_key[j]) {
            ++i;
        } else if (b.read_key[j] < a.read_key[i]) {
            ++j;
        } else {
            const double qa = a.q0[i], qb = b.q0[j];
            // Probability the read's two alleles sit on ONE haplotype, against the two orderings.
            // These sum to 1, so the ratio below needs no normalisation.
            const double same = qa * qb + (1.0 - qa) * (1.0 - qb);
            const double diff = 1.0 - same;
            // With probability `pr` the read reports the true relation; otherwise it reports a coin
            // flip. That mixture is what keeps a single mismapped read from swamping the sum -- it
            // is the escape term, and without it the statistic is unusable.
            const double pr = (double)a.p[i] * (double)b.p[j];
            const double cis = pr * same + (1.0 - pr) * 0.5;
            const double trans = pr * diff + (1.0 - pr) * 0.5;
            if (cis > 0.0 && trans > 0.0) {
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

unordered_set<size_t> read_phase_flips(vector<PhaseSite>& sites, const ReadPhasingParams& params,
                                       ReadPhasingCounters& counters) {
    unordered_set<size_t> flips;
    if (sites.empty()) {
        return flips;
    }
    // Phase is only comparable inside a block, so the blocks are the unit of work.
    // `record_key` breaks the tie, and it has to: sites arrive in whatever order the parallel
    // render queues produced them, two sites can share a position, and an unstable sort over a
    // non-total order would then chain them differently from run to run. The phase would be
    // reproducible only by accident.
    sort(sites.begin(), sites.end(), [](const PhaseSite& x, const PhaseSite& y) {
        if (x.phase_set != y.phase_set) {
            return x.phase_set < y.phase_set;
        }
        if (x.position != y.position) {
            return x.position < y.position;
        }
        return x.record_key < y.record_key;
    });
    for (PhaseSite& s : sites) {
        // The merge in `phase_link` needs both sides ordered. Done once here rather than per link.
        vector<size_t> order(s.read_key.size());
        for (size_t i = 0; i < order.size(); ++i) {
            order[i] = i;
        }
        sort(order.begin(), order.end(),
             [&](size_t x, size_t y) { return s.read_key[x] < s.read_key[y]; });
        vector<uint64_t> k(order.size());
        vector<float> q(order.size()), pp(order.size());
        for (size_t i = 0; i < order.size(); ++i) {
            k[i] = s.read_key[order[i]];
            q[i] = s.q0[order[i]];
            pp[i] = s.p[order[i]];
        }
        s.read_key.swap(k);
        s.q0.swap(q);
        s.p.swap(pp);
    }
    counters.sites += sites.size();

    size_t begin = 0;
    while (begin < sites.size()) {
        size_t end = begin;
        while (end < sites.size() && sites[end].phase_set == sites[begin].phase_set) {
            ++end;
        }
        ++counters.chains;
        const size_t n = end - begin;

        // Which sites may carry a link at all.
        vector<size_t> rel;
        vector<size_t> unrel;
        for (size_t t = 0; t < n; ++t) {
            if (sites[begin + t].reliability >= params.reliability) {
                rel.push_back(t);
            } else {
                unrel.push_back(t);
            }
        }
        counters.reliable += rel.size();

        // `o[t] == 1` means "swap this site against the order the panel gave it". A chain's first
        // reliable site is pinned at 0, so with no read evidence nothing moves.
        vector<int> o(n, 0);
        vector<char> decided(n, 0);

        for (size_t pass = 0; pass < 2 && !rel.empty(); ++pass) {
            // --- stage 1: the chain of reliable sites, stepping over the rest ---
            vector<double> d(rel.size() > 0 ? rel.size() - 1 : 0, 0.0);
            for (size_t m = 0; m + 1 < rel.size(); ++m) {
                d[m] = phase_link(sites[begin + rel[m]], sites[begin + rel[m + 1]], params.cap);
            }
            vector<size_t> bounds;
            bounds.push_back(0);
            for (size_t m = 0; m + 1 < rel.size(); ++m) {
                if (std::fabs(d[m]) < params.break_threshold) {
                    bounds.push_back(m + 1);
                }
            }
            bounds.push_back(rel.size());
            counters.breaks += bounds.size() - 2;

            for (size_t b = 0; b + 1 < bounds.size(); ++b) {
                const size_t s0 = bounds[b], s1 = bounds[b + 1];
                o[rel[s0]] = 0;
                decided[rel[s0]] = 1;
                for (size_t m = s0; m + 1 < s1; ++m) {
                    o[rel[m + 1]] = o[rel[m]] ^ (d[m] < 0.0 ? 1 : 0);
                    decided[rel[m + 1]] = 1;
                }
            }

            // --- stage 2: relink consecutive reliable blocks ---
            for (size_t b = 0; b + 2 < bounds.size(); ++b) {
                const size_t ea = bounds[b + 1];              // one past A's last reliable site
                const size_t sb = bounds[b + 1];              // B's first reliable site
                const size_t eb = bounds[b + 2];
                const size_t a_first = ea > params.relink ? ea - params.relink : 0;
                const size_t a_lo = max(bounds[b], a_first);
                const size_t b_hi = min(eb, sb + params.relink);
                double total = 0.0;
                for (size_t i = a_lo; i < ea; ++i) {
                    for (size_t j = sb; j < b_hi; ++j) {
                        const double dv =
                            phase_link(sites[begin + rel[i]], sites[begin + rel[j]], params.cap);
                        if (dv == 0.0) {
                            continue;
                        }
                        // Referred back to the two blocks' own boundary sites, whose relative
                        // orientation is what is being decided; each block's internal parity is
                        // already settled so it cancels out of the question.
                        const int flip = (o[rel[i]] ^ o[rel[ea - 1]]) ^ (o[rel[j]] ^ o[rel[sb]]);
                        total += flip ? -dv : dv;
                    }
                }
                if (total == 0.0) {
                    ++counters.breaks_no_reads;
                    continue;                                  // the panel's frame stands
                }
                const int x = total < 0.0 ? 1 : 0;
                if (x ^ o[rel[ea - 1]] ^ o[rel[sb]]) {
                    for (size_t m = sb; m < eb; ++m) {
                        o[rel[m]] ^= 1;
                    }
                }
            }

            // --- changepoint pass: flip at junctions the whole-read evidence rejects ---
            //
            // Stages 1 and 2 decide a junction from the reads its two ADJACENT sites share. A read
            // spanning three or more sites also constrains the 1->3 relation, and nothing above
            // ever reads that. Here it is read.
            //
            // For read r assigned haplotype H, write a_s for log10 of its per-site term under H and
            // b_s for the term under the other haplotype. Flipping every orientation downstream of
            // junction j changes that read's likelihood by exactly the suffix sum of (b_s - a_s)
            // over its own sites at or after j -- and by nothing at all if it does not span j,
            // since a read wholly downstream is merely relabelled. So
            //
            //     S(j) = sum over reads spanning j of sum_{s >= j} (b_s - a_s)
            //
            // and S(j) > 0 means the flip raises sum_r max_H L(H, r). Greedily flipping at the
            // argmax junction therefore hill-climbs that objective and terminates.
            //
            // Accumulated with a difference array: read r contributes (b_i - a_i) to every junction
            // in (first site of r, i], which is one range update per site.
            if (params.changepoint_min > 0.0 && rel.size() >= 3) {
                // read key -> its (chain index, q0, p) over this chain's reliable sites
                unordered_map<uint64_t, vector<std::array<double, 3>>> by_read;
                by_read.reserve(rel.size() * 4);
                for (size_t m = 0; m < rel.size(); ++m) {
                    const PhaseSite& st = sites[begin + rel[m]];
                    for (size_t i = 0; i < st.read_key.size(); ++i) {
                        by_read[st.read_key[i]].push_back(
                            {(double)m, (double)st.q0[i], (double)st.p[i]});
                    }
                }
                for (size_t round = 0; round < params.changepoint_rounds; ++round) {
                    vector<double> diff(rel.size() + 2, 0.0);
                    for (auto& kv : by_read) {
                        auto& obs = kv.second;
                        if (obs.size() < 2) {
                            continue;   // spans no junction
                        }
                        // Orientation applied, then the read's own haplotype chosen.
                        double l0 = 0.0, l1 = 0.0;
                        for (const auto& o3 : obs) {
                            const size_t m = (size_t)o3[0];
                            const double p = o3[2];
                            const double q = o[rel[m]] ? 1.0 - o3[1] : o3[1];
                            l0 += log10(p * q + (1.0 - p) * 0.5);
                            l1 += log10(p * (1.0 - q) + (1.0 - p) * 0.5);
                        }
                        const bool hap1 = l1 > l0;
                        const size_t lo = (size_t)obs.front()[0];
                        for (size_t k = 1; k < obs.size(); ++k) {
                            const size_t m = (size_t)obs[k][0];
                            const double p = obs[k][2];
                            const double q = o[rel[m]] ? 1.0 - obs[k][1] : obs[k][1];
                            const double t_own = hap1 ? p * (1.0 - q) + (1.0 - p) * 0.5
                                                     : p * q + (1.0 - p) * 0.5;
                            const double t_alt = hap1 ? p * q + (1.0 - p) * 0.5
                                                     : p * (1.0 - q) + (1.0 - p) * 0.5;
                            const double v = log10(t_alt) - log10(t_own);
                            // junctions (lo, m], indexed by the site they precede
                            diff[lo + 1] += v;
                            diff[m + 1] -= v;
                        }
                    }
                    double run = 0.0, best = 0.0;
                    size_t best_j = 0;
                    for (size_t j = 1; j < rel.size(); ++j) {
                        run += diff[j];
                        if (run > best) {
                            best = run;
                            best_j = j;
                        }
                    }
                    if (best_j == 0 || best < params.changepoint_min) {
                        break;
                    }
                    for (size_t m = best_j; m < rel.size(); ++m) {
                        o[rel[m]] ^= 1;
                    }
                    ++counters.changepoints;
                    counters.changepoint_gain += best;
                    if (round + 1 == params.changepoint_rounds) {
                        ++counters.changepoint_capped;
                    }
                }
            }

            // --- coherence pass: demote sites whose reads disagree with their own neighbourhood ---
            //
            // `reliability` asks whether this site's reads separate its two alleles. It says nothing
            // about whether those reads belong where the rest of their evidence puts them, and
            // measured against chr20's switch positions it is nearly blind: worst-1% enrichment 1.1x,
            // and below the base rate at 5%. Coherence -- the share of a site's reads whose allele
            // here matches the haplotype their OTHER sites imply -- reaches 5.1x at 1% and 9.3x at
            // 0.1%. The two correlate at r = 0.437, so this is new information, not a restatement.
            //
            // Held out by construction: each read's haplotype is recomputed with THIS site's own
            // term removed, so a site never votes on itself.
            if (params.coherence_min > 0.0 && pass == 0 && rel.size() >= 3) {
                unordered_map<uint64_t, vector<std::array<double, 3>>> by_read;
                for (size_t m = 0; m < rel.size(); ++m) {
                    const PhaseSite& st = sites[begin + rel[m]];
                    for (size_t i = 0; i < st.read_key.size(); ++i) {
                        by_read[st.read_key[i]].push_back(
                            {(double)m, (double)st.q0[i], (double)st.p[i]});
                    }
                }
                vector<size_t> ok(rel.size(), 0), tot(rel.size(), 0);
                for (auto& kv : by_read) {
                    auto& obs = kv.second;
                    if (obs.size() < 2) {
                        continue;   // no other site to be held out against
                    }
                    double l0 = 0.0, l1 = 0.0;
                    vector<double> t0(obs.size()), t1(obs.size());
                    for (size_t k = 0; k < obs.size(); ++k) {
                        const size_t m = (size_t)obs[k][0];
                        const double p = obs[k][2];
                        const double q = o[rel[m]] ? 1.0 - obs[k][1] : obs[k][1];
                        t0[k] = log10(p * q + (1.0 - p) * 0.5);
                        t1[k] = log10(p * (1.0 - q) + (1.0 - p) * 0.5);
                        l0 += t0[k];
                        l1 += t1[k];
                    }
                    for (size_t k = 0; k < obs.size(); ++k) {
                        const size_t m = (size_t)obs[k][0];
                        const bool hap1 = (l1 - t1[k]) > (l0 - t0[k]);   // leave THIS site out
                        const double q = o[rel[m]] ? 1.0 - obs[k][1] : obs[k][1];
                        const bool says1 = q < 0.5;
                        ++tot[m];
                        if (says1 == hap1) {
                            ++ok[m];
                        }
                    }
                }
                vector<size_t> keep, drop;
                for (size_t m = 0; m < rel.size(); ++m) {
                    if (tot[m] >= 10
                        && (double)ok[m] / (double)tot[m] < params.coherence_min) {
                        drop.push_back(rel[m]);
                    } else {
                        keep.push_back(rel[m]);
                    }
                }
                if (!drop.empty() && keep.size() >= 2) {
                    counters.demoted_incoherent += drop.size();
                    for (size_t t : drop) {
                        unrel.push_back(t);
                    }
                    sort(unrel.begin(), unrel.end());
                    rel.swap(keep);
                    std::fill(o.begin(), o.end(), 0);
                    std::fill(decided.begin(), decided.end(), 0);
                    continue;                    // re-run stages 1 and 2 on the surviving sites
                }
            }
            break;                               // no demotion: one pass is all that is needed
        }

        if (!rel.empty()) {
            // --- stage 3: hang the unreliable sites off the settled chain ---
            for (size_t t : unrel) {
                double s = 0.0;
                size_t used = 0;
                // Nearest reliable neighbours on both sides. `hang` is split between them rather
                // than taken from one, so a site near a chain end is not decided one-sidedly.
                for (int dir = -1; dir <= 1; dir += 2) {
                    size_t taken = 0;
                    for (size_t idx = 0; idx < rel.size() && taken < params.hang / 2 + 1; ++idx) {
                        const size_t m = dir < 0 ? rel.size() - 1 - idx : idx;
                        const size_t u = rel[m];
                        if (dir < 0 ? !(u < t) : !(u > t)) {
                            continue;
                        }
                        const double dv = phase_link(sites[begin + t], sites[begin + u],
                                                     params.cap);
                        ++taken;
                        if (dv == 0.0) {
                            continue;
                        }
                        const int pred = o[u] ^ (dv < 0.0 ? 1 : 0);
                        s += std::fabs(dv) * (pred == 0 ? 1.0 : -1.0);
                        ++used;
                    }
                }
                if (params.panel_weight > 0.0 && !rel.empty()) {
                    // The panel says this site keeps the order it came in with, relative to its
                    // neighbour's. Weak, and only decisive when the reads say nothing.
                    size_t nearest = rel.front();
                    size_t best = nearest > t ? nearest - t : t - nearest;
                    for (size_t u : rel) {
                        const size_t dist = u > t ? u - t : t - u;
                        if (dist < best) {
                            best = dist;
                            nearest = u;
                        }
                    }
                    s += params.panel_weight * (o[nearest] == 0 ? 1.0 : -1.0);
                }
                if (used == 0) {
                    ++counters.hung_no_reads;
                }
                o[t] = s > 0.0 ? 0 : 1;
                decided[t] = 1;
                ++counters.hung;
            }
        }

        for (size_t t = 0; t < n; ++t) {
            if (decided[t] && o[t]) {
                flips.insert(sites[begin + t].record_key);
                ++counters.flipped;
            }
        }
        begin = end;
    }
    return flips;
}

}
