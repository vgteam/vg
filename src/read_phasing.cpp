#include "read_phasing.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <unordered_map>

namespace vg {

using std::log10;
using std::max;
using std::min;
using std::pair;
using std::sort;
using std::unordered_map;

double read_confidence(float q0, float c) {
    const double win = max((double)q0, 1.0 - (double)q0) * (double)c;
    return -10.0 * log10(max(1.0 - win, 1e-12));
}

void merge_mates(PhaseSite& site) {
    const size_t n = site.read_key.size();
    vector<size_t> order(n);
    for (size_t i = 0; i < n; ++i) {
        order[i] = i;
    }
    // By key, and within a key the preferred row first. Every tie is broken on the rows' values,
    // so two rows that compare equal are identical and either may be kept.
    sort(order.begin(), order.end(), [&](size_t x, size_t y) {
        if (site.read_key[x] != site.read_key[y]) {
            return site.read_key[x] < site.read_key[y];
        }
        const double cx = read_confidence(site.q0[x], site.c[x]);
        const double cy = read_confidence(site.q0[y], site.c[y]);
        if (cx != cy) {
            return cx > cy;
        }
        if (site.q0[x] != site.q0[y]) {
            return site.q0[x] > site.q0[y];
        }
        return site.c[x] > site.c[y];
    });
    vector<uint64_t> k;
    vector<float> q, pp;
    k.reserve(n);
    q.reserve(n);
    pp.reserve(n);
    for (size_t i = 0; i < n; ++i) {
        const size_t r = order[i];
        if (!k.empty() && k.back() == site.read_key[r]) {
            continue;
        }
        k.push_back(site.read_key[r]);
        q.push_back(site.q0[r]);
        pp.push_back(site.c[r]);
    }
    site.read_key.swap(k);
    site.q0.swap(q);
    site.c.swap(pp);
}

PhaseSite reduce_to_pair(const PhaseReadEvidence& pe, size_t a0, size_t a1) {
    PhaseSite site;
    // Slot 0 holds `a0` and slot 1 `a1`, so the weights follow the alleles into their slots.
    const vector<double> weight = allele_length_weights(
        pe.allele_length, pe.n_alleles, pe.mean_read_length, pe.length_weighted,
        vector<int>{(int)a0, (int)a1});
    for (size_t r = 0; r < pe.num_reads(); ++r) {
        const double e = (double)pe.mismap[r];
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
    // Both read lists are sorted by key, one row per key, so this is a merge rather than a lookup
    // per read.
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
    // The merge in `phase_link` needs both sides ordered, and coherence needs a read to appear
    // once per site. Done once here rather than per link, for each site on its own.
#pragma omp parallel for schedule(dynamic, 1024)
    for (size_t i = 0; i < sites.size(); ++i) {
        merge_mates(sites[i]);
    }
    counters.sites += sites.size();

    // Each phase set is a range of `sites` and is worked on independently of every other one, so
    // the phase sets are worked on in parallel, longest first. Each keeps its own counters and its
    // own list of the sites whose chosen pair is to be swapped. Afterwards the counters are added
    // up and the lists gathered in phase-set order, which gives what one loop over the phase sets,
    // in order, gives.
    vector<pair<size_t, size_t>> phase_sets;
    for (size_t begin = 0; begin < sites.size();) {
        size_t end = begin;
        while (end < sites.size() && sites[end].phase_set == sites[begin].phase_set) {
            ++end;
        }
        phase_sets.emplace_back(begin, end);
        begin = end;
    }
    vector<size_t> by_length(phase_sets.size());
    for (size_t i = 0; i < by_length.size(); ++i) {
        by_length[i] = i;
    }
    std::stable_sort(by_length.begin(), by_length.end(), [&](size_t x, size_t y) {
        return phase_sets[x].second - phase_sets[x].first > phase_sets[y].second - phase_sets[y].first;
    });
    vector<ReadPhasingCounters> set_counters(phase_sets.size());
    vector<vector<size_t>> set_flips(phase_sets.size());

#pragma omp parallel for schedule(dynamic, 1)
    for (size_t rank = 0; rank < by_length.size(); ++rank) {
        const size_t set_index = by_length[rank];
        const size_t begin = phase_sets[set_index].first;
        const size_t end = phase_sets[set_index].second;
        // This phase set's own counters and list (see above). The name `counters` hides the
        // function's parameter on purpose, so that the loop body below reads as it did when the
        // phase sets were worked on one at a time.
        ReadPhasingCounters& counters = set_counters[set_index];
        vector<size_t>& flipped_keys = set_flips[set_index];
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
        // Breaks in the chain as the last pass left it. A pass that demotes sites runs stages 1
        // and 2 again, so only the last pass's breaks are counted.
        size_t breaks = 0;
        size_t breaks_no_reads = 0;
        size_t rounds = 0;
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
            breaks = bounds.size() - 2;
            breaks_no_reads = 0;

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
                        // already chosen so it cancels out of the question.
                        const int flip = (o[rel[i]] ^ o[rel[ea - 1]]) ^ (o[rel[j]] ^ o[rel[sb]]);
                        total += flip ? -dv : dv;
                    }
                }
                if (total == 0.0) {
                    ++breaks_no_reads;
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
                    ++rounds;
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
        counters.breaks += breaks;
        counters.breaks_no_reads += breaks_no_reads;
        counters.coherence_rounds_run = max(counters.coherence_rounds_run, rounds);

        if (!rel.empty()) {
            // --- stage 4: hang the unreliable sites from the chain ---
            //
            // `rel` is in increasing order and holds no unreliable site, so the reliable sites
            // before `t` are rel[0, split) and those after it rel[split, end), found by binary
            // search rather than by walking `rel` for every unreliable site.
            const size_t per_side = params.hang / 2 + 1;
            for (size_t t : unrel) {
                double s = 0.0;
                size_t used = 0;
                const size_t split = std::lower_bound(rel.begin(), rel.end(), t) - rel.begin();
                auto link = [&](size_t u) {
                    const double dv = phase_link(sites[begin + t], sites[begin + u], params.cap);
                    if (dv == 0.0) {
                        return;
                    }
                    const int pred = o[u] ^ (dv < 0.0 ? 1 : 0);
                    s += std::fabs(dv) * (pred == 0 ? 1.0 : -1.0);
                    ++used;
                };
                // Nearest reliable neighbours on both sides, nearest first, before `t` and then
                // after it. `hang` is split between them rather than taken from one, so a site
                // near a chain end is not decided one-sidedly.
                for (size_t k = 0; k < per_side && k < split; ++k) {
                    link(rel[split - 1 - k]);
                }
                for (size_t k = split; k < rel.size() && k - split < per_side; ++k) {
                    link(rel[k]);
                }
                if (params.panel_weight > 0.0 && !rel.empty()) {
                    // The panel says this site keeps the order it came in with, relative to its
                    // neighbour's. Weak, and only decisive when the reads say nothing. The
                    // nearest reliable site is on one side of `t` or the other; on a tie, the
                    // one before it.
                    size_t nearest;
                    if (split == 0) {
                        nearest = rel.front();
                    } else if (split == rel.size()) {
                        nearest = rel.back();
                    } else {
                        nearest = t - rel[split - 1] <= rel[split] - t ? rel[split - 1]
                                                                       : rel[split];
                    }
                    s += params.panel_weight * (o[nearest] == 0 ? 1.0 : -1.0);
                }
                if (used == 0) {
                    ++counters.hung_no_reads;
                }
                // Flipped only on a vote against the panel's order; a tie keeps it.
                o[t] = s < 0.0 ? 1 : 0;
                decided[t] = 1;
                ++counters.hung;
            }
        }

        for (size_t t = 0; t < n; ++t) {
            if (decided[t] && o[t]) {
                flipped_keys.push_back(sites[begin + t].record_key);
                ++counters.flipped;
            }
        }
    }

    for (size_t i = 0; i < phase_sets.size(); ++i) {
        const ReadPhasingCounters& c = set_counters[i];
        counters.chains += c.chains;
        counters.reliable += c.reliable;
        counters.breaks += c.breaks;
        counters.breaks_no_reads += c.breaks_no_reads;
        counters.hung += c.hung;
        counters.hung_no_reads += c.hung_no_reads;
        counters.flipped += c.flipped;
        counters.demoted_incoherent += c.demoted_incoherent;
        counters.coherence_unconverged += c.coherence_unconverged;
        counters.coherence_rounds_run = max(counters.coherence_rounds_run, c.coherence_rounds_run);
        for (size_t key : set_flips[i]) {
            flips.insert(key);
        }
    }
    return flips;
}

vector<double> allele_length_weights(const vector<std::uint32_t>& allele_length, size_t n_alleles,
                                     float mean_read_length, bool length_weighted,
                                     const vector<int>& slot_allele) {
    // Expected share of the site's reads per slot, from the alleles' full lengths: an allele of
    // length L yields a read overlapping the site from L + R - 1 start positions. Flat when the
    // lengths are unavailable or `length_weighted` is false.
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
