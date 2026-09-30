#include "read_phasing.hpp"

#include <algorithm>
#include <array>
#include <cmath>
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
            // Probability that the read's two alleles lie on one strand in the current orders.
            // These sum to 1, so the ratio below needs no normalisation.
            const double same = qa * qb + (1.0 - qa) * (1.0 - qb);
            const double diff = 1.0 - same;
            // With probability `pr` the read reports the true relation, and otherwise a coin flip,
            // so that one mismapped read cannot dominate the sum.
            const double pr = (double)a.c[i] * (double)b.c[j];
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
    // Phase is comparable only inside a phase set, so each phase set is worked on separately.
    // `record_key` breaks ties between sites at the same position, so that the order does not
    // depend on the order in which the sites arrived.
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
            pp[i] = s.c[order[i]];
        }
        s.read_key.swap(k);
        s.q0.swap(q);
        s.c.swap(pp);
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

        const size_t max_pass = params.coherence_min > 0.0
                                    ? max<size_t>(1, params.coherence_rounds) + 1
                                    : 1;
        for (size_t pass = 0; pass < max_pass && !rel.empty(); ++pass) {
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

            // --- stage 2: relink the unbroken pieces of the chain ---
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

            // --- stage 3: remove chain sites whose reads disagree with the rest of the chain ---
            //
            // A site's coherence is the share of its reads whose allele here matches the strand
            // their other chain sites imply. Each read's strand is computed with this site left
            // out, so a site does not vote on itself.
            if (params.coherence_min > 0.0 && pass + 1 < max_pass && rel.size() >= 3) {
                unordered_map<uint64_t, vector<std::array<double, 3>>> by_read;
                for (size_t m = 0; m < rel.size(); ++m) {
                    const PhaseSite& st = sites[begin + rel[m]];
                    for (size_t i = 0; i < st.read_key.size(); ++i) {
                        by_read[st.read_key[i]].push_back(
                            {(double)m, (double)st.q0[i], (double)st.c[i]});
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
                        const bool hap1 = (l1 - t1[k]) > (l0 - t0[k]);   // leave this site out
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
                    // A site needs at least ten counted reads before low coherence can remove it.
                    if (tot[m] >= 10 && (double)ok[m] / (double)tot[m] < params.coherence_min) {
                        drop.push_back(rel[m]);
                    } else {
                        keep.push_back(rel[m]);
                    }
                }
                if (!drop.empty() && keep.size() >= 2) {
                    counters.demoted_incoherent += drop.size();
                    ++counters.coherence_rounds_run;
                    if (pass + 2 == max_pass) {
                        // Still removing sites in the last allowed round, so the chain is
                        // reported as not converged.
                        ++counters.coherence_unconverged;
                    }
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
            // --- stage 4: hang the unreliable sites from the chain ---
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

vector<double> allele_length_weights(const vector<std::uint32_t>& allele_length, size_t n_alleles,
                                     float mean_read_length, bool length_weighted,
                                     const vector<int>& slot_allele) {
    // Expected share of the site's reads per slot, from the alleles' full lengths: an allele of
    // length L yields a read overlapping the site from L + R - 1 start positions. Flat when the
    // lengths are unavailable or under --flat-mixture.
    const size_t n_slots = slot_allele.size();
    vector<double> weight(n_slots, n_slots ? 1.0 / (double)n_slots : 0.0);
    if (length_weighted && mean_read_length > 0.0
        && allele_length.size() == n_alleles && n_slots > 1) {
        double total = 0.0;
        vector<double> raw(n_slots, 0.0);
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
