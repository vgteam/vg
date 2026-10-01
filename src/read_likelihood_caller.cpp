#include "read_likelihood_caller.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <iomanip>
#include <set>
#include <sstream>

#include "statistics.hpp"

namespace vg {

thread_local bool ReadLikelihoodSnarlCaller::want_alt_ploidy = false;
thread_local int ReadLikelihoodSnarlCaller::region_ploidy = 0;

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

void ReadLikelihoodSnarlCaller::set_support_available(bool available) {
    this->support_available = available;
}

function<bool(const SnarlTraversal&, int)> ReadLikelihoodSnarlCaller::get_skip_allele_fn() const {
    if (support_available) {
        // A pack file is present, so use the inherited pruning.
        return SupportBasedSnarlCaller::get_skip_allele_fn();
    }
    // No real support to prune on. SnarlCaller::get_skip_allele_fn() would assert.
    return [](const SnarlTraversal&, int) { return false; };
}

bool ReadLikelihoodSnarlCaller::traversals_equal(const SnarlTraversal& a,
                                                 const SnarlTraversal& b) {
    // The protobuf operator== also compares visits to child snarls, which a node-by-node
    // comparison would miss. The traversals reaching this caller are node paths.
    return a == b;
}

pair<vector<int>, unique_ptr<SnarlCaller::CallInfo>> ReadLikelihoodSnarlCaller::genotype(
    const Snarl& snarl, const vector<SnarlTraversal>& traversals, int ref_trav_idx, int ploidy,
    const string& ref_path_name, pair<size_t, size_t> ref_range) {

    ReadLikelihoodCallInfo* call_info = new ReadLikelihoodCallInfo();
    call_info->ploidy = ploidy;
    unique_ptr<CallInfo> call_info_owner(call_info);

    if (traversals.empty() || ploidy < 1) {
        return make_pair(vector<int>(), std::move(call_info_owner));
    }

    // Build the reads x alleles matrix for this site.
    AlleleReadLikelihoods matrix = likelihood_calculator.compute(
        snarl, traversals, region_ploidy > 0 ? region_ploidy : ploidy);

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
    // Kept for later: the anchors are built when the record is rendered, from the settled
    // genotype.
    call_info->anchor_evidence = std::move(matrix.anchor_evidence);
    call_info->phase_evidence = std::move(matrix.phase_evidence);

    if (dump_stream != nullptr) {
        stringstream site_name;
        site_name << snarl.start().node_id() << (snarl.start().backward() ? "-" : "+") << "_"
                  << snarl.end().node_id() << (snarl.end().backward() ? "-" : "+");
#pragma omp critical (read_likelihood_dump)
        matrix.dump(*dump_stream, site_name.str());
    }

    if (matrix.num_reads() == 0) {
        // No read can tell the alleles apart here, so every genotype is equally likely,
        // and we make no call rather than calling the reference.
        return make_pair(vector<int>(), std::move(call_info_owner));
    }

    // Derive the call at ploidy p from the matrix, which does not depend on ploidy. It runs once
    // for the site's ploidy and, when set_want_alt_ploidy asked for it, once for the other, so
    // that a nested site's record can later be built at whichever ploidy its parent settles on.
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
    if (want_alt_ploidy && traversals.size() > 1) {
        int other = ploidy == 1 ? 2 : 1;
        auto alt = make_unique<ReadLikelihoodCallInfo>();
        // Copy the fields that depend only on the matrix, not on the ploidy.
        alt->n_reads = call_info->n_reads;
        alt->scored_traversals = call_info->scored_traversals;
        alt->allele_support = call_info->allele_support;
        // BL, also a property of the matrix.
        alt->mean_best_ln = call_info->mean_best_ln;
        alt->ploidy = other;
        vector<int> alt_best = derive(other, alt.get());
        if (!alt_best.empty()) {
            call_info->alt_ploidy_best = alt_best;
            call_info->alt_ploidy_info = std::move(alt);
        }
    }

    return make_pair(best_genotype, std::move(call_info_owner));
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

    // Observed reads over what the call predicts, written whether or not the depth term is on.
    if (info->depth_ratio >= 0.0) {
        variant.format.push_back("DR");
        stringstream ss;
        ss << std::fixed << std::setprecision(3) << info->depth_ratio;
        variant.samples[sample_name]["DR"].push_back(ss.str());
    }

    // Map each emitted VCF allele back to the matrix column it came from.
    //
    // emit_variant merged alleles with the same sequence and dropped uncalled ones, so
    // these indices are not the ones we genotyped. Traversals that match nothing, such as
    // the empty traversal of a star allele, stay unmapped.
    vector<int> site_to_scored(traversals.size(), -1);
    for (size_t s = 0; s < traversals.size(); ++s) {
        for (size_t k = 0; k < info->scored_traversals.size(); ++k) {
            if (traversals_equal(traversals[s], info->scored_traversals[k])) {
                site_to_scored[s] = (int)k;
                break;
            }
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
    // from the posterior of the all-reference genotype where we have it.
    variant.quality = 0;
    if (!genotype.empty()) {
        bool is_ref_call = all_of(genotype.begin(), genotype.end(), [](int a) { return a == 0; });
        double ref_posterior = 0;
        bool have_ref = false;
        if (all_mapped) {
            vector<int> ref_genotype(info->ploidy, site_to_scored.empty() ? 0 : site_to_scored[0]);
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
              "divided by the number the called genotype predicts, from a read rate measured over "
              "the read source's local fetch window and the called alleles' traversal lengths. 1.0 "
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
              "calibrated posterior. GQI is the unscaled value; --no-share-quality restores "
              "it as GQ. With --depth-quality A in effect, records whose called alleles change "
              "length by at least 50 bp are additionally scaled by exp(-A * |ln DR|), so a call "
              "whose read count is implausible for the sequence it claims ranks lower\">\n";
    header += "##FORMAT=<ID=GQI,Number=1,Type=Integer,Description=\"Genotype Quality from the "
              "likelihood ratio alone, with no explained-read-fraction scaling. Equals GQ "
              "when --no-share-quality is in effect\">\n";
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

}
