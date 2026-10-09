#include "read_likelihood_caller.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <limits>
#include <iomanip>
#include <set>
#include <sstream>

#include "statistics.hpp"
#include "utility.hpp"

namespace vg {


using namespace std;

ReadLikelihoodSnarlCaller::ReadLikelihoodSnarlCaller(const PathHandleGraph& graph,
                                                     SnarlManager& snarl_manager,
                                                     TraversalSupportFinder& support_finder,
                                                     AlleleLikelihoodCalculator& likelihood_calculator)
    : SupportBasedSnarlCaller(graph, snarl_manager, support_finder),
      likelihood_calculator(likelihood_calculator) {
}

ReadLikelihoodSnarlCaller::~ReadLikelihoodSnarlCaller() {
}

void ReadLikelihoodSnarlCaller::set_likelihood_dump(ostream* dump_stream) {
    this->dump_stream = dump_stream;
}

void ReadLikelihoodSnarlCaller::set_share_discount(bool discount) {
    this->share_discount = discount;
}

void ReadLikelihoodSnarlCaller::set_depth_quality(double exponent, size_t min_length) {
    this->depth_quality = exponent;
    this->depth_quality_min_length = min_length;
}


double ReadLikelihoodSnarlCaller::explained_share(const ReadLikelihoodCallInfo& info,
                                                  const vector<int>& called) {
    // From the fractional allele_support rather than the rounded AD, and over the distinct
    // called alleles, so that a homozygote does not count its allele twice.
    double explained = 0.0;
    set<int> seen;
    for (int a : called) {
        if (a >= 0 && (size_t)a < info.allele_support.size() && seen.insert(a).second) {
            explained += info.allele_support[a];
        }
    }
    const size_t n = info.n_reads;
    // Clamped, since rounding in the fractional counts can take the sum just above the read
    // count, and a share above 1 would raise GQ.
    return n ? min(1.0, explained / (double)n) : 1.0;
}

double ReadLikelihoodSnarlCaller::discounted_gq(const ReadLikelihoodCallInfo& info,
                                                const vector<int>& called, double gap) const {
    double gq = gap;
    if (share_discount) {
        gq *= explained_share(info, called);
    }
    return gq * info.depth_discount;
}

double ReadLikelihoodSnarlCaller::gq_factor(const ReadLikelihoodCallInfo& info) const {
    return (share_discount ? info.explained_share : 1.0) * info.depth_discount;
}

void ReadLikelihoodSnarlCaller::recompute_gq(ReadLikelihoodCallInfo& info) const {
    const vector<int>* best = nullptr;
    double best_ll = -numeric_limits<double>::infinity();
    double second_ll = -numeric_limits<double>::infinity();
    for (const auto& kv : info.genotype_lls) {
        if (kv.second > best_ll) {
            second_ll = best_ll;
            best_ll = kv.second;
            best = &kv.first;
        } else if (kv.second > second_ll) {
            second_ll = kv.second;
        }
    }
    if (best != nullptr && std::isfinite(best_ll) && std::isfinite(second_ll)) {
        info.gq = discounted_gq(info, *best,
                                logprob_to_phred(second_ll) - logprob_to_phred(best_ll));
    }
}

void ReadLikelihoodSnarlCaller::set_min_confidence(double threshold) {
    this->min_confidence = threshold;
}

pair<vector<int>, unique_ptr<SnarlCaller::CallInfo>> ReadLikelihoodSnarlCaller::genotype(
    const Snarl&, const vector<SnarlTraversal>&, int, int, const string&, pair<size_t, size_t>) {
    throw std::logic_error("ReadLikelihoodSnarlCaller::genotype is not used; MultiPassCaller "
                           "genotypes through genotype_at");
}

pair<vector<int>, unique_ptr<ReadLikelihoodSnarlCaller::ReadLikelihoodCallInfo>>
ReadLikelihoodSnarlCaller::genotype_at(const SiteBounds& site,
                                       const vector<Traversal>& traversals, int ref_trav_idx,
                                       const Ploidies& ploidies,
                                       const vector<SiteBounds>& enclosing,
                                       const string& ref_path_name,
                                       pair<size_t, size_t> ref_range) {
    const int ploidy = ploidies.ploidy;
    ReadLikelihoodCallInfo* call_info = new ReadLikelihoodCallInfo();
    call_info->ploidy = ploidy;
    unique_ptr<ReadLikelihoodCallInfo> call_info_owner(call_info);

    if (traversals.empty() || ploidy < 1) {
        return make_pair(vector<int>(), std::move(call_info_owner));
    }

    // Build the reads x alleles matrix for this site.
    AlleleReadLikelihoods matrix = likelihood_calculator.compute(
        site, traversals, enclosing,
        ploidies.region_ploidy > 0 ? ploidies.region_ploidy : ploidy);

    // Per-allele read support and mean absolute fit. Neither enters the genotype likelihood;
    // both are written to the VCF, as AD and BL.
    call_info->allele_support.assign(matrix.num_alleles(), 0.0);
    double best_ln_total = 0.0;
    size_t best_ln_n = 0;
    for (size_t r = 0; r < matrix.num_reads(); ++r) {
        double best = 0.0;
        for (size_t a = 0; a < matrix.num_alleles(); ++a) {
            best = max(best, matrix.rel(r, a));
        }
        if (best > 0.0) {
            size_t winners = 0;
            for (size_t a = 0; a < matrix.num_alleles(); ++a) {
                if (matrix.rel(r, a) >= best) {
                    ++winners;
                }
            }
            for (size_t a = 0; a < matrix.num_alleles(); ++a) {
                if (matrix.rel(r, a) >= best) {
                    call_info->allele_support[a] += 1.0 / (double)winners;
                }
            }
        }
        double bl = matrix.best_ln_likelihood(r);
        if (std::isfinite(bl)) {
            best_ln_total += bl;
            ++best_ln_n;
        }
    }
    call_info->mean_best_ln = best_ln_n ? best_ln_total / (double)best_ln_n : 0.0;

    call_info->n_reads = matrix.num_reads();
    call_info->scored_traversals.assign(traversals.begin(), traversals.end());
    if (matrix.has_depth_context()) {
        call_info->depth_lengths.reserve(matrix.num_alleles());
        for (size_t a = 0; a < matrix.num_alleles(); ++a) {
            call_info->depth_lengths.push_back(matrix.traversal_length(a));
        }
        call_info->depth_rate = matrix.depth_rate_per_haplotype();
        call_info->depth_read_length = matrix.depth_read_length_used();
        call_info->depth_observed = matrix.observed_reads();
    }
    // The CallInfo outlives the matrix, so the per-read evidence moves into it: read phasing and
    // re-genotyping read it in the rounds, and the render builds the anchors from it.
    call_info->anchor_evidence = std::move(matrix.anchor_evidence);
    call_info->phase_evidence = std::move(matrix.phase_evidence);

    if (dump_stream != nullptr) {
        stringstream site_name;
        site_name << graph.get_id(site.start) << (graph.get_is_reverse(site.start) ? "-" : "+")
                  << "_" << graph.get_id(site.end) << (graph.get_is_reverse(site.end) ? "-" : "+");
#pragma omp critical (read_likelihood_dump)
        matrix.dump(*dump_stream, site_name.str());
    }

    if (matrix.num_reads() == 0) {
        // No read can tell the alleles apart here, so every genotype is equally likely,
        // and we make no call rather than calling the reference.
        return make_pair(vector<int>(), std::move(call_info_owner));
    }

    // Derive the call at ploidy p from the matrix, which does not depend on ploidy. It runs once
    // for the site's ploidy and, when `ploidies` asks for it, once for the other, so that
    // a nested site's record can later be built at whichever ploidy its parent's chosen genotype
    // gives it.
    auto derive = [&](int p, ReadLikelihoodCallInfo* info) -> vector<int> {
        // Score every genotype: for K alleles and ploidy 2 that is K(K+1)/2 genotypes, each
        // one pass over the reads.
            vector<pair<vector<int>, double>> scored = matrix.score_genotypes(p);
            if (scored.empty()) {
                return vector<int>();
            }

        double second_best_ll = -numeric_limits<double>::infinity();
        double total_ll = -numeric_limits<double>::infinity();

        // Break exact ties in favour of the all-reference genotype. When every read fits
        // every allele equally, all genotypes score the same, and taking whichever came
        // first would make a non-reference call on no evidence. traversals[0] is not
        // necessarily the reference; ref_trav_idx says which one is.
        size_t best_index = 0;
        if (ref_trav_idx >= 0 && (size_t)ref_trav_idx < traversals.size()) {
            vector<int> ref_genotype(p, ref_trav_idx);
            for (size_t i = 0; i < scored.size(); ++i) {
                if (scored[i].first == ref_genotype) {
                    best_index = i;
                    break;
                }
            }
        }
        double best_ll = scored[best_index].second;

        // Keep the runner-up's identity, not just its score: achievable_gap depends on how
        // the two genotypes differ.
        size_t second_best_index = scored.size();

        for (size_t i = 0; i < scored.size(); ++i) {
            double ll = scored[i].second;
            info->genotype_lls[scored[i].first] = ll;

            if (i != best_index) {
                if (ll > best_ll) {
                    second_best_ll = best_ll;
                    second_best_index = best_index;
                    best_ll = ll;
                    best_index = i;
                } else if (ll > second_best_ll) {
                    second_best_ll = ll;
                    second_best_index = i;
                }
            }

            total_ll = (i == 0) ? ll : add_log(total_ll, ll);
        }

        // GQ is the phred-scaled gap between the best and second best genotype, as in vg's
        // other callers.
        info->gq = 0;
        if (std::isfinite(best_ll) && std::isfinite(second_best_ll)) {
            info->gq = logprob_to_phred(second_best_ll) - logprob_to_phred(best_ll);
        }
        info->gq_undiscounted = info->gq;

        // GQN: the same gap as a fraction of the achievable gap. It is computed from the
        // unclamped log-likelihoods, since GQ is capped at 256 in the VCF. Both are in nats,
        // so the ratio needs no conversion. The explained-share discount is applied below.
        double achievable_gap = 0.0;
        if (second_best_index < scored.size() && std::isfinite(best_ll)
            && std::isfinite(second_best_ll)) {
            achievable_gap = matrix.achievable_gap(scored[best_index].first,
                                                   scored[second_best_index].first);
            info->achievable_gap = achievable_gap;
            if (achievable_gap > 0.0) {
                info->gq_fraction =
                    min(1.0, max(0.0, (best_ll - second_best_ll) / achievable_gap));
            }
        }

        {
            const vector<int>& called = scored[best_index].first;
            info->explained_share = explained_share(*info, called);

            // GQN takes the same discount as GQ, since the likelihood gap cannot see reads
            // that fit an uncalled allele. It is applied whatever --no-share-quality says,
            // since that option controls GQ only.
            if (info->gq_fraction >= 0.0) {
                info->gq_fraction *= info->explained_share;
            }

            info->depth_ratio = matrix.depth_ratio(called);

            // The depth discount (see set_depth_quality). A call's size is the largest
            // change in length of a called allele against the reference traversal, using the
            // lengths without boundary nodes that lambda uses, since the VCF alleles do not
            // exist yet.
            if (depth_quality > 0.0 && info->depth_ratio > 0.0 && ref_trav_idx >= 0) {
                size_t ref_len = matrix.traversal_length((size_t)ref_trav_idx);
                size_t change = 0;
                for (int a : called) {
                    if (a < 0) {
                        continue;
                    }
                    size_t len = matrix.traversal_length((size_t)a);
                    change = max(change, len > ref_len ? len - ref_len : ref_len - len);
                }
                if (change >= depth_quality_min_length) {
                    info->depth_discount = exp(-depth_quality * fabs(log(info->depth_ratio)));
                }
            }
            info->gq = discounted_gq(*info, called, info->gq_undiscounted);
        }

        // ln posterior under a uniform prior, which cancels, so the posterior is the
        // genotype's likelihood over the sum of all genotypes' likelihoods.
        info->posterior = std::isfinite(total_ll) ? best_ll - total_ll : 0.0;

            return scored[best_index].first;
    };

    vector<int> best_genotype = derive(ploidy, call_info);
    if (best_genotype.empty()) {
        return make_pair(vector<int>(), std::move(call_info_owner));
    }

    // The same site at the other ploidy, from the same matrix. The depth rate is per haplotype
    // of the region, so it does not change with the number of haplotypes crossing the site.
    if (ploidies.also_score_other && traversals.size() > 1) {
        int other = ploidy == 1 ? 2 : 1;
        auto alt = make_unique<ReadLikelihoodCallInfo>();
        // Copy the fields that depend only on the matrix, not on the ploidy.
        alt->n_reads = call_info->n_reads;
        alt->scored_traversals = call_info->scored_traversals;
        alt->allele_support = call_info->allele_support;
        // BL, also a property of the matrix.
        alt->mean_best_ln = call_info->mean_best_ln;
        alt->depth_lengths = call_info->depth_lengths;
        alt->depth_rate = call_info->depth_rate;
        alt->depth_read_length = call_info->depth_read_length;
        alt->depth_observed = call_info->depth_observed;
        alt->ploidy = other;
        vector<int> alt_best = derive(other, alt.get());
        if (!alt_best.empty()) {
            call_info->alt_ploidy_best = alt_best;
            call_info->alt_ploidy_info = std::move(alt);
        }
    }

    return make_pair(best_genotype, std::move(call_info_owner));
}

double ReadLikelihoodSnarlCaller::ReadLikelihoodCallInfo::depth_ratio_of(
        const vector<int>& scored_genotype) const {
    if (depth_lengths.empty()) {
        return -1.0;
    }
    // The arithmetic of AlleleReadLikelihoods::depth_ratio, so that the direct call gets back
    // exactly the `depth_ratio` it was given.
    double expected = AlleleReadLikelihoods::expected_reads_from(depth_lengths, depth_rate,
                                                                 depth_read_length,
                                                                 scored_genotype);
    return expected > 0.0 ? depth_observed / expected : -1.0;
}

const PhaseReadEvidence* ReadLikelihoodSnarlCaller::ReadLikelihoodCallInfo::read_phasing_evidence(
        PhaseReadEvidence& scratch) const {
    const PhaseReadEvidence* pe = phase_evidence.get();
    if (pe == nullptr && anchor_evidence != nullptr) {
        const AnchorSiteEvidence& ev = *anchor_evidence;
        scratch.n_alleles = ev.n_alleles;
        scratch.allele_length = ev.allele_length;
        scratch.mean_read_length = ev.mean_read_length;
        scratch.length_weighted = ev.length_weighted;
        scratch.rel = ev.rel;
        scratch.read_key.reserve(ev.reads.size());
        scratch.mismap.reserve(ev.reads.size());
        for (const AnchorRead& r : ev.reads) {
            // The same key the calculator gives a read when it fills `phase_evidence` itself.
            scratch.read_key.push_back((uint64_t)std::hash<string_view>{}(read_names().name(r.read)));
            scratch.mismap.push_back(r.mismap);
        }
        pe = &scratch;
    }
    if (pe != nullptr && (pe->n_alleles == 0 || pe->num_reads() == 0)) {
        return nullptr;
    }
    return pe;
}

void ReadLikelihoodSnarlCaller::update_vcf_info(const Snarl& snarl,
                                               const vector<SnarlTraversal>& traversals,
                                               const vector<int>& genotype,
                                               const unique_ptr<CallInfo>& call_info,
                                               const string& sample_name,
                                               vcflib::Variant& variant) {

    const ReadLikelihoodCallInfo* info =
        dynamic_cast<const ReadLikelihoodCallInfo*>(call_info.get());
    if (info == nullptr) {
        // A genotype derived from a parent site rather than scored here (nested
        // calling does this). There is no matrix to report.
        return;
    }

    // Number of informative reads at the site.
    variant.format.push_back("DP");
    variant.samples[sample_name]["DP"].push_back(std::to_string(info->n_reads));

    // Mean absolute fit. The model does not use it, but it says whether the reads fit any
    // allele, which GQ does not, so a filter can combine the two.
    variant.format.push_back("BL");
    {
        stringstream ss;
        ss << std::fixed << std::setprecision(2) << info->mean_best_ln;
        variant.samples[sample_name]["BL"].push_back(ss.str());
    }

    // Map each emitted VCF allele back to the matrix column it came from.
    //
    // emit_variant merged alleles with the same sequence and dropped uncalled ones, so
    // these indices are not the ones we genotyped. Traversals that match nothing, such as
    // the empty traversal of a star allele, stay unmapped.
    vector<int> site_to_scored(traversals.size(), -1);
    for (size_t s = 0; s < traversals.size(); ++s) {
        for (size_t k = 0; k < info->scored_traversals.size(); ++k) {
            if (same_walk(graph, traversals[s], info->scored_traversals[k])) {
                site_to_scored[s] = (int)k;
                break;
            }
        }
    }

    // Observed reads over what the written genotype predicts, written whether or not the depth
    // term is on. The written genotype is the direct call's unless the linkage model moved the
    // record. It is used where each of its alleles maps to a scored column and it has the
    // ploidy the site was genotyped at; otherwise DR is the direct call's.
    {
        double depth_ratio = info->depth_ratio;
        vector<int> scored_written;
        bool mapped = true;
        for (int allele : genotype) {
            if (allele < 0) {
                continue;
            }
            if ((size_t)allele >= site_to_scored.size() || site_to_scored[allele] < 0) {
                mapped = false;
                break;
            }
            scored_written.push_back(site_to_scored[allele]);
        }
        if (mapped && (int)scored_written.size() == info->ploidy) {
            // Sorted, as the direct call's genotype is, so that the sum runs in the same order.
            sort(scored_written.begin(), scored_written.end());
            double written = info->depth_ratio_of(scored_written);
            if (written >= 0.0) {
                depth_ratio = written;
            }
        }
        if (depth_ratio >= 0.0) {
            variant.format.push_back("DR");
            stringstream ss;
            ss << std::fixed << std::setprecision(3) << depth_ratio;
            variant.samples[sample_name]["DR"].push_back(ss.str());
        }
    }

    // AD over the emitted alleles, through the same mapping, rounded to integers. It need not
    // sum to DP: a read whose best allele was scored but not emitted has no entry. An allele
    // with no scored column, such as a star allele, gets 0, since AD needs one entry per
    // allele.
    {
        vector<long> ad(traversals.size(), 0);
        for (size_t s = 0; s < traversals.size(); ++s) {
            int k = site_to_scored[s];
            if (k >= 0 && (size_t)k < info->allele_support.size()) {
                ad[s] = lround(info->allele_support[k]);
            }
        }
        variant.format.push_back("AD");
        for (long v : ad) {
            variant.samples[sample_name]["AD"].push_back(std::to_string(v));
        }
    }

    // GL over the emitted alleles, in VCF genotype order: one entry per genotype of the
    // emitted alleles, looked up in the genotypes we scored.
    bool all_mapped = true;
    for (size_t s = 0; s < traversals.size(); ++s) {
        if (site_to_scored[s] < 0) {
            all_mapped = false;
            break;
        }
    }

    // A genotype with a star or missing allele was called at a lower ploidy than the record
    // shows: in nested calling only some of the parent's haplotypes pass through the child.
    // GL would then have the wrong length for GT's ploidy, so it is left out.
    bool genotype_has_marker =
        any_of(genotype.begin(), genotype.end(), [](int a) { return a < 0; });

    if (all_mapped && !genotype_has_marker && info->ploidy >= 1) {
        vector<string> gl_strings;
        bool complete = true;

        for (auto& site_genotype :
             AlleleReadLikelihoods::enumerate_genotypes(traversals.size(), info->ploidy)) {
            // Translate to matrix columns and sort, since genotype_lls is keyed
            // by the sorted multiset.
            vector<int> scored_genotype;
            scored_genotype.reserve(site_genotype.size());
            for (int site_allele : site_genotype) {
                scored_genotype.push_back(site_to_scored[site_allele]);
            }
            sort(scored_genotype.begin(), scored_genotype.end());

            auto found = info->genotype_lls.find(scored_genotype);
            if (found == info->genotype_lls.end()) {
                complete = false;
                break;
            }
            // VCF wants GL log10-scaled.
            gl_strings.push_back(std::to_string(ln_to_log10(found->second)));
        }

        if (complete && !gl_strings.empty()) {
            variant.format.push_back("GL");
            for (auto& gl : gl_strings) {
                variant.samples[sample_name]["GL"].push_back(gl);
            }
        }
    }

    variant.format.push_back("GQ");
    variant.samples[sample_name]["GQ"].push_back(
        std::to_string(min((int)256, max((int)0, (int)info->gq))));

    // GQ before both discounts, the explained share's (off under --no-share-quality) and
    // --depth-quality's, written whether or not either is on. It is the gap the site's own
    // likelihoods gave: when --regenotype recomputes GQ from corrected likelihoods, GQI is not
    // recomputed with it.
    variant.format.push_back("GQI");
    variant.samples[sample_name]["GQI"].push_back(
        std::to_string(min((int)256, max((int)0, (int)info->gq_undiscounted))));

    // GQN, written as a fraction rather than a phred score. "." where there was no gap to
    // normalise, which is different from 0.
    variant.format.push_back("GQN");
    {
        std::ostringstream gqn;
        if (info->gq_fraction < 0.0) {
            gqn << ".";
        } else {
            gqn << std::fixed << std::setprecision(3) << info->gq_fraction;
        }
        variant.samples[sample_name]["GQN"].push_back(gqn.str());
    }

    // GP is in natural log, unlike GL, which the VCF specification puts in log10. The header
    // says so.
    variant.format.push_back("GP");
    variant.samples[sample_name]["GP"].push_back(std::to_string(info->posterior));

    // QUAL as the phred-scaled probability that the site is not variant, taken
    // from the posterior of the all-reference genotype where we have it. Only the reference
    // allele's column is needed, so a record with an unscored allele, such as a star allele,
    // still gets one.
    variant.quality = 0;
    if (!genotype.empty()) {
        bool is_ref_call = all_of(genotype.begin(), genotype.end(), [](int a) { return a == 0; });
        double ref_posterior = 0;
        bool have_ref = false;
        if (!site_to_scored.empty() && site_to_scored[0] >= 0) {
            vector<int> ref_genotype(info->ploidy, site_to_scored[0]);
            sort(ref_genotype.begin(), ref_genotype.end());
            auto found = info->genotype_lls.find(ref_genotype);
            if (found != info->genotype_lls.end()) {
                // Renormalise against the full scored set.
                double total = -numeric_limits<double>::infinity();
                bool first = true;
                for (auto& entry : info->genotype_lls) {
                    total = first ? entry.second : add_log(total, entry.second);
                    first = false;
                }
                if (std::isfinite(total)) {
                    ref_posterior = found->second - total;
                    have_ref = true;
                }
            }
        }
        if (have_ref) {
            variant.quality = is_ref_call ? 0 : logprob_to_phred(ref_posterior);
        }
    }

    // Low-confidence records are marked, not dropped: the linkage model rewrites genotypes in
    // VCFOutputCaller::write_variants after this runs, and may fix them.
    variant.filter = "PASS";
    if (info->n_reads == 0) {
        variant.filter = "noreads";
    } else if (min_confidence > 0.0 && info->gq_fraction >= 0.0
               && info->gq_fraction < min_confidence) {
        // Only where GQN exists: a site with no gap to normalise has no GQN to compare.
        variant.filter = "lowconf";
    }
}

void ReadLikelihoodSnarlCaller::update_vcf_header(string& header) const {
    header += "##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Number of informative reads "
              "overlapping the site\">\n";
    header += "##FORMAT=<ID=AD,Number=R,Type=Integer,Description=\"Reads whose best-fitting "
        "allele is this one. **AD does not sum to DP**, for two reasons: a read fitting "
        "several alleles equally splits its vote between them, and more importantly only "
        "alleles that reached this record get a column, while the genotyper scored every "
        "allele the site offered. At a site where many alleles were enumerated and few "
        "emitted, most reads best-fit something absent here and the shortfall is large. "
        "That shortfall is itself informative: it is how much of the evidence the emitted "
        "alleles fail to explain. Not used by the genotype model, which gives each "
        "haplotype a share of the reads set by the length of sequence unique to its allele "
        "(an equal share under --flat-mixture), whatever the reads show\">\n";
    header += "##FORMAT=<ID=DR,Number=1,Type=Float,Description=\"Observed reads at this site "
              "divided by the number the called genotype predicts, from the rate of read starts per "
              "base in the site's rate window, the fixed 16,384 bp reference bucket holding the "
              "site and one bucket on each side, and the called alleles' traversal lengths. 1.0 "
              "means the read count is exactly what the call implies. Values well above 1 are "
              "collapsed repeats, where reads from several copies pile onto one; values near 0.5 "
              "are a genotype claiming twice the sequence actually covered, which is what a missed "
              "heterozygous deletion looks like. Both sides of the ratio count a read as 1 - e_r, "
              "the probability it came from this locus at all, so a site whose mapping quality "
              "matches its neighbourhood's is unaffected and only the difference shows; pass "
              "--depth-count-raw to count whole reads instead. Reported whether or not "
              "--depth-term is in use\">\n";
    header += "##FORMAT=<ID=BL,Number=1,Type=Float,Description=\"Mean over reads of the best "
        "raw alignment score any allele gave them. Measures whether reads fit anything at "
        "this site, where GQ measures only the gap between the top two genotypes, so the "
        "two are nearly independent. NOT normalised for site size or read overlap, so it "
        "is comparable between calls at similar sites rather than across a whole genome\">\n";
    header += "##FORMAT=<ID=GL,Number=G,Type=Float,Description=\"Genotype Likelihood, "
              "log10-scaled P(reads | genotype) from the read-level likelihood model. Useful for "
              "ranking; not a calibrated probability, and over-confident at high depth because "
              "reads are treated as independent\">\n";
    header += "##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"Genotype Quality, the "
              "phred-scaled gap between the best and second-best genotype likelihood, scaled "
              "by the fraction of reads the called genotype explains (sum(AD)/DP). The "
              "likelihood ratio alone cannot see reads that fit an allele outside the call, "
              "because those reads enter every genotype's likelihood and cancel; the scaling "
              "restores them. It also means GQ here is a quality score rather than a "
              "calibrated posterior. --no-share-quality turns this scaling off. With "
              "--depth-quality A in effect, records whose called alleles change length by at "
              "least 50 bp are also scaled by exp(-A * |ln DR|), so a call whose read count is "
              "implausible for the sequence it claims ranks lower. GQI is the value with neither "
              "scaling\">\n";
    header += "##FORMAT=<ID=GQI,Number=1,Type=Integer,Description=\"Genotype Quality from the "
              "likelihood ratio alone, with neither the explained-read-fraction scaling nor the "
              "--depth-quality scaling. Equals GQ under --no-share-quality, except where "
              "--depth-quality scales GQ\">\n";
    header += "##FORMAT=<ID=GQN,Number=1,Type=Float,Description=\"Normalised Genotype Quality: "
              "the likelihood-ratio gap as a fraction, in [0,1], of the gap this site could have "
              "produced had every read fitted the call perfectly. Unlike GQ it means the same "
              "thing at any depth and any ploidy, so one threshold works across a 5x diploid and "
              "a 15x haploid contig. GQ cannot: it scales with depth, and it scales with ploidy "
              "the other way, because at ploidy 1 the runner-up genotype is a different allele "
              "outright and every read discriminates fully, where a diploid heterozygote's "
              "runner-up differs on one strand only. On HG002 that makes hemizygous chrX calls "
              "run a median GQ of 247 where chr7 diploid homozygotes at the same depth run 46. "
              "Dividing GQ by depth does not fix this and measurably makes it worse, since it "
              "corrects the smaller axis and leaves the larger. '.' where there was no gap to "
              "normalise, which is not 0. NEGATIVE, and so in [-1,1] rather than [0,1], on a "
              "record the linkage layer moved: there GQN is the settled genotype's margin over "
              "the best alternative on the same scale, and a negative value says the panel "
              "overrode the reads. Those records are the caller's highest false-positive density "
              "-- 37.8% against 8.6% overall on 44x ONT chr20 -- so they are reported as a signed "
              "number rather than blanked, which is what left --min-confidence unable to see "
              "them\">\n";
    header += "##FORMAT=<ID=GP,Number=1,Type=Float,Description=\"Genotype Probability, the "
              "natural-log-scaled posterior of the called genotype under a uniform prior. "
              "Nats, not log10: unlike GL, exponentiate with e, so -2.303 means p=0.1\">\n";
    header += "##FILTER=<ID=lowconf,Description=\"GQN below --min-confidence: the call used less "
              "of the discrimination this site could offer than required. Unlike a GQ threshold "
              "this means the same thing at any depth and any ploidy -- requiring GQ >= 10 costs a "
              "5x diploid contig a third of its F1, where GQN >= 0.05 costs it 0.009 and gains a "
              "haploid contig 0.017. Marks rather than drops. Off unless --min-confidence is "
              "given\">\n";
    header += "##FILTER=<ID=noreads,Description=\"No informative read overlaps the site, so the "
              "read-level model has no evidence either way\">\n";
}


bool ReadLikelihoodSnarlCaller::rewrite_quality_for_chosen_genotype(
        string& vcf_line, const LinkageCollector::MovedQuality& moved, double lowconf_threshold) {
    // A record whose genotype the linkage model changed gets quality fields for its chosen
    // genotype, since the per-site quality fields describe the genotype the reads alone chose.
    vector<string> fields;
    split_delims_keep_empty(vcf_line, "\t", fields);
    if (fields.size() < 10) {
        return false;
    }
    vector<string> keys, values;
    split_delims_keep_empty(fields[8], ":", keys);
    split_delims_keep_empty(fields[9], ":", values);
    if (keys.size() != values.size()) {
        return false;
    }
    size_t gq_field = keys.size(), gqi_field = keys.size(), gqn_field = keys.size();
    size_t gl_field = keys.size(), gt_field = keys.size();
    for (size_t i = 0; i < keys.size(); ++i) {
        if (keys[i] == "GQ") {
            gq_field = i;
        } else if (keys[i] == "GQI") {
            gqi_field = i;
        } else if (keys[i] == "GQN") {
            gqn_field = i;
        } else if (keys[i] == "GL") {
            gl_field = i;
        } else if (keys[i] == "GT") {
            gt_field = i;
        }
    }

    // GQN's divisor, in phred units as GL's margin is below.
    const double achievable_phred = 10.0 * moved.direct.achievable_gap / log(10.0);
    if (gq_field != keys.size()) {
        // GQ becomes the phred-scaled complement of the posterior, multiplied by the factor the
        // direct call's GQ was, and capped at GQI. The posterior includes the panel's frequency
        // prior, so 1 - posterior can understate the uncertainty where the reads were weakest; the
        // cap keeps GQ within what the reads alone support, and makes GQ <= GQI hold on every
        // record rewritten here.
        const double posterior = moved.posterior;
        double q = posterior >= 1.0 ? 256.0 : -10.0 * log10(max(1.0 - posterior, 1e-26));
        q *= moved.direct.gq_factor;
        if (gqi_field != keys.size()) {
            try {
                q = min(q, stod(values[gqi_field]));
            } catch (const std::exception&) {
                // GQI absent or unparsable: keep the discounted posterior quality, uncapped.
            }
        }
        // Truncated, not rounded, as the per-site GQ is, so that equal qualities print the same.
        values[gq_field] = std::to_string((int)min(256.0, max(0.0, q)));
    }
    // GQN, recomputed for the chosen genotype: its likelihood margin over the best alternative,
    // as a fraction of the direct call's achievable gap, times its explained share, as the
    // per-site GQN is. It is negative when the linkage model moved the call against the reads,
    // which is why the sign is kept. It stays "." when GL is absent, the genotype cannot be read,
    // or the direct call had no achievable gap.
    bool gqn_known = false;
    double gqn_new = 0.0;
    if (gqn_field != keys.size() && gl_field != keys.size() && gt_field != keys.size()
        && achievable_phred > 0.0) {
        vector<double> gl;
        bool parsed = true;
        {
            size_t start = 0;
            while (parsed) {
                size_t comma = values[gl_field].find(',', start);
                string tok = values[gl_field].substr(
                    start, comma == string::npos ? string::npos : comma - start);
                try {
                    gl.push_back(stod(tok));
                } catch (const std::exception&) {
                    parsed = false;
                }
                if (comma == string::npos) {
                    break;
                }
                start = comma + 1;
            }
        }
        // The called genotype's index in the GL. A diploid record's GL is indexed j(j+1)/2 + i
        // for i <= j, and a haploid record's by allele; which applies is checked against the GL's
        // length (n against n(n+1)/2). A "." field is dropped rather than skipping the record,
        // since `1|.` is a nested chain on one strand of a diploid parent, with a real margin.
        vector<int> called;
        if (parsed) {
            size_t start = 0;
            while (true) {
                size_t sep = values[gt_field].find_first_of("/|", start);
                string tok = values[gt_field].substr(
                    start, sep == string::npos ? string::npos : sep - start);
                if (tok != "." && !tok.empty()) {
                    try {
                        called.push_back(std::stoi(tok));
                    } catch (const std::exception&) {
                        called.clear();
                        break;
                    }
                }
                if (sep == string::npos) {
                    break;
                }
                start = sep + 1;
            }
        }
        // How many alleles this record has, REF included, used only to recognise a haploid GL by
        // its length. The diploid case does not check the length: under --atomize-blocks a snarl's
        // block records share its GL while each has only its own ALTs, so the lengths need not
        // match.
        size_t n_alleles = 1;
        if (fields[4] != "." && !fields[4].empty()) {
            vector<string> alt_list;
            split_delims_keep_empty(fields[4], ",", alt_list);
            n_alleles += alt_list.size();
        }
        const bool diploid_gl = gl.size() == n_alleles * (n_alleles + 1) / 2;
        const bool haploid_gl = gl.size() == n_alleles;
        size_t idx = gl.size();
        if (parsed && called.size() == 2 && !gl.empty()) {
            int i = min(called[0], called[1]);
            int j = max(called[0], called[1]);
            idx = (size_t)(j * (j + 1) / 2 + i);
        } else if (parsed && called.size() == 1 && haploid_gl && !diploid_gl) {
            idx = (size_t)called[0];
        }
        if (idx < gl.size()) {
            double best_other = -std::numeric_limits<double>::infinity();
            for (size_t g = 0; g < gl.size(); ++g) {
                if (g != idx) {
                    best_other = max(best_other, gl[g]);
                }
            }
            if (best_other > -std::numeric_limits<double>::infinity()) {
                double margin_phred = 10.0 * (gl[idx] - best_other);
                gqn_new = min(1.0, max(-1.0, margin_phred / achievable_phred
                                                 * moved.direct.explained_share));
                gqn_known = true;
            }
        }
    }
    if (gqn_field != keys.size()) {
        if (gqn_known) {
            char buf[32];
            snprintf(buf, sizeof(buf), "%.3f", gqn_new);
            values[gqn_field] = buf;
        } else {
            values[gqn_field] = ".";
        }
    }
    // `lowconf` was set from the direct pass's GQN. Decide it again from the recomputed GQN where there
    // is one, and clear it where there is not.
    if (gqn_known && lowconf_threshold > 0.0) {
        fields[6] = gqn_new < lowconf_threshold ? "lowconf" : "PASS";
    } else if (fields[6] == "lowconf") {
        fields[6] = "PASS";
    }

    fields[9] = join_delim(values, ':');
    vcf_line = join_delim(fields, '\t');
    return true;
}

}
