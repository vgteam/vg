#include "read_likelihood_caller.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <iomanip>
#include <set>
#include <sstream>

#include "statistics.hpp"

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


void ReadLikelihoodSnarlCaller::set_min_confidence(double threshold) {
    this->min_confidence = threshold;
}

void ReadLikelihoodSnarlCaller::set_support_available(bool available) {
    this->support_available = available;
}

function<bool(const SnarlTraversal&, int)> ReadLikelihoodSnarlCaller::get_skip_allele_fn() const {
    if (support_available) {
        // A pack file is present, so keep the established pruning behaviour.
        return SupportBasedSnarlCaller::get_skip_allele_fn();
    }
    // No real support to prune on. Note this must not fall through to
    // SnarlCaller::get_skip_allele_fn(), whose default asserts.
    return [](const SnarlTraversal&, int) { return false; };
}

bool ReadLikelihoodSnarlCaller::traversals_equal(const SnarlTraversal& a, const SnarlTraversal& b) {
    if (a.visit_size() != b.visit_size()) {
        return false;
    }
    for (int64_t i = 0; i < a.visit_size(); ++i) {
        if (a.visit(i).node_id() != b.visit(i).node_id() ||
            a.visit(i).backward() != b.visit(i).backward()) {
            return false;
        }
    }
    return true;
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
    AlleleReadLikelihoods matrix = likelihood_calculator.compute(snarl, traversals, ploidy);

    // Per-allele read support and mean absolute fit. Neither enters the genotype
    // likelihood -- the row normalisation divides the absolute fit out, and the mixture
    // uses only which alleles a read fits, not how the reads divide between them. Both
    // are reported so that a caller downstream can use evidence the model discards.
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

    call_info->n_informative = matrix.num_reads();
    call_info->n_unplaceable = matrix.num_unplaceable();
    call_info->scored_traversals.assign(traversals.begin(), traversals.end());

    if (dump_stream != nullptr) {
        stringstream site_name;
        site_name << snarl.start().node_id() << (snarl.start().backward() ? "-" : "+") << "_"
                  << snarl.end().node_id() << (snarl.end().backward() ? "-" : "+");
#pragma omp critical (read_likelihood_dump)
        matrix.dump(*dump_stream, site_name.str());
    }

    if (matrix.num_reads() == 0) {
        // No informative read overlaps this site, so every genotype is equally
        // likely. That is the correct answer for a depth-agnostic model, and it
        // is deliberately a no-call rather than a confident hom-ref: absence of
        // reads is absence of evidence here, not evidence of the reference.
        return make_pair(vector<int>(), std::move(call_info_owner));
    }

    // Enumerate every genotype exhaustively. For K alleles and ploidy 2 that is
    // K(K+1)/2 genotypes, each costing one pass over the reads, so there is no
    // need for the candidate pruning the Poisson caller does with top_k/top_m.
    vector<pair<vector<int>, double>> scored = matrix.score_genotypes(ploidy);
    if (scored.empty()) {
        return make_pair(vector<int>(), std::move(call_info_owner));
    }

    double second_best_ll = -numeric_limits<double>::infinity();
    double total_ll = -numeric_limits<double>::infinity();

    // Break exact ties in favour of the all-reference genotype.
    //
    // Reads can be present yet completely uninformative between alleles -- every
    // row flat, so every genotype scores identically. The likelihood genuinely has
    // no preference there, so something outside it has to choose, and taking
    // whichever traversal happens to sit at index 0 would emit a confident-looking
    // non-reference call on no evidence at all. Note traversals[0] is *not*
    // necessarily the reference; ref_trav_idx says which one is.
    //
    // This is a tie-break convention rather than a prior: it fires only on exact
    // equality, so it cannot move a site where the reads say anything at all.
    size_t best_index = 0;
    if (ref_trav_idx >= 0 && (size_t)ref_trav_idx < traversals.size()) {
        vector<int> ref_genotype(ploidy, ref_trav_idx);
        for (size_t i = 0; i < scored.size(); ++i) {
            if (scored[i].first == ref_genotype) {
                best_index = i;
                break;
            }
        }
    }
    double best_ll = scored[best_index].second;

    // The runner-up's identity, not just its score: achievable_gap needs to know which
    // genotype the gap is against, since what a site could have achieved depends on how
    // the two genotypes differ -- one strand or both.
    size_t second_best_index = scored.size();

    for (size_t i = 0; i < scored.size(); ++i) {
        double ll = scored[i].second;
        call_info->genotype_lls[scored[i].first] = ll;

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

    // GQ is the phred-scaled gap between the best and second best genotype,
    // matching the existing callers' convention.
    call_info->gq = 0;
    if (std::isfinite(best_ll) && std::isfinite(second_best_ll)) {
        call_info->gq = logprob_to_phred(second_best_ll) - logprob_to_phred(best_ll);
    }
    call_info->gq_undiscounted = call_info->gq;

    // The same gap as a fraction of what this site could have produced, which is what
    // makes it comparable across depth and ploidy. Computed here, from the raw
    // log-likelihoods, deliberately: GQ is clamped to 256 on the way into the VCF, and
    // that clamp censors 23% of haploid calls at full depth, so a consumer cannot
    // reconstruct this from the emitted fields.
    //
    // Both sides use natural-log likelihoods, so no phred conversion is needed -- the
    // ratio is scale-free and taking it in phred would give the same number by a longer
    // route. The explained-share discount is applied below, once it is known.
    double achievable_gap = 0.0;
    if (second_best_index < scored.size() && std::isfinite(best_ll)
        && std::isfinite(second_best_ll)) {
        achievable_gap = matrix.achievable_gap(scored[best_index].first,
                                               scored[second_best_index].first);
        if (achievable_gap > 0.0) {
            call_info->gq_fraction =
                min(1.0, max(0.0, (best_ll - second_best_ll) / achievable_gap));
        }
    }

    // How much of the pile-up the called genotype accounts for. Reads whose best-fitting
    // allele lies outside the call are counted in *every* genotype's likelihood and so
    // cancel out of the best-versus-second-best comparison; GQ cannot see them. See
    // set_share_discount for why this is worth correcting and what it costs.
    //
    // Computed from the fractional allele_support rather than the rounded AD, and over
    // the distinct called alleles, so a homozygote does not count its allele twice.
    {
        const vector<int>& called = scored[best_index].first;
        double explained = 0.0;
        set<int> seen;
        for (int a : called) {
            if (a >= 0 && (size_t)a < call_info->allele_support.size() && seen.insert(a).second) {
                explained += call_info->allele_support[a];
            }
        }
        size_t n = matrix.num_reads();
        // Clamped rather than trusted: ties split fractionally and floating-point
        // accumulation can land a hair above the read count, and a share above 1 would
        // *raise* GQ, which this is never allowed to do.
        call_info->explained_share = n ? min(1.0, explained / (double)n) : 1.0;

        // GQN carries the same discount as GQ, and is not merely the raw ratio. The
        // likelihood gap cannot see reads that fit an allele outside the call, because
        // those enter every genotype's likelihood and cancel; normalising a quantity
        // that is blind to them leaves it blind. Measured over the titration the
        // discount is worth most of the remaining calibration: 0.316 without it against
        // 0.260 with, where raw GQ scores 0.347. Applied whatever --no-share-quality
        // says, since that flag is about what GQ reports and GQN is a different field
        // whose whole purpose is to be comparable.
        if (call_info->gq_fraction >= 0.0) {
            call_info->gq_fraction *= call_info->explained_share;
        }

        if (share_discount) {
            call_info->gq *= call_info->explained_share;
        }

        call_info->depth_ratio = matrix.depth_ratio(called);

        // The depth-implausibility discount, gated on called-allele size. See
        // set_depth_quality for why it is gated and why it is not on by default.
        //
        // Size is measured as the largest length change against the reference
        // traversal, using the same interior lengths lambda uses -- so it is the change
        // in sequence a read could actually be recruited by, and it needs no reference
        // to the VCF alleles, which do not exist yet at this point.
        if (depth_quality > 0.0 && call_info->depth_ratio > 0.0 && ref_trav_idx >= 0) {
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
                call_info->gq *= exp(-depth_quality * fabs(log(call_info->depth_ratio)));
            }
        }
    }

    // Posterior under a uniform prior. Deliberately NOT the Poisson caller's
    // formula, which subtracts ln(number of candidates): under a uniform prior
    // that term cancels analytically and does not belong. It is nearly harmless
    // at a fixed candidate count, but this caller enumerates exhaustively, so
    // the count varies with the number of alleles and the term would make the
    // reported posterior vary with it too.
    call_info->posterior = std::isfinite(total_ll) ? best_ll - total_ll : 0.0;

    return make_pair(scored[best_index].first, std::move(call_info_owner));
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
    variant.samples[sample_name]["DP"].push_back(std::to_string(info->n_informative));

    // Mean absolute fit. Reported rather than used: the model divides it out, but it
    // separates true from false calls about as well as GQ does and is nearly
    // uncorrelated with it, so a downstream filter can combine them.
    variant.format.push_back("BL");
    {
        stringstream ss;
        ss << std::fixed << std::setprecision(2) << info->mean_best_ln;
        variant.samples[sample_name]["BL"].push_back(ss.str());
    }

    // Observed reads over what the call predicts. Diagnostic: the model uses it only
    // when --depth-term is armed, but it is always reported so the signal can be
    // measured before it is trusted.
    if (info->depth_ratio >= 0.0) {
        variant.format.push_back("DR");
        stringstream ss;
        ss << std::fixed << std::setprecision(3) << info->depth_ratio;
        variant.samples[sample_name]["DR"].push_back(ss.str());
    }

    // Map each emitted VCF allele back to the matrix column it came from.
    //
    // emit_variant deduplicated by allele string and dropped uncalled alleles, so
    // these indices are not the ones we genotyped. Structural comparison recovers
    // the mapping exactly, and is immune to however the allele strings were
    // flattened. Entries that match nothing -- notably the empty placeholder
    // traversal standing in for a star allele -- stay unmapped.
    vector<int> site_to_scored(traversals.size(), -1);
    for (size_t s = 0; s < traversals.size(); ++s) {
        for (size_t k = 0; k < info->scored_traversals.size(); ++k) {
            if (traversals_equal(traversals[s], info->scored_traversals[k])) {
                site_to_scored[s] = (int)k;
                break;
            }
        }
    }

    // AD over the emitted alleles, through the same remap. Rounded to integers for the
    // conventional Number=R Integer form. The total does not reconstruct DP, and at a
    // busy site falls a long way below it: reads whose best-fitting allele was scored
    // but never emitted have no column to land in, and a read that fits several alleles
    // equally splits its vote. The shortfall is the useful part -- it is the share of
    // reads the called genotype fails to explain -- so the header documents it rather
    // than papering over it. An allele that maps to no scored column (a star allele)
    // reports 0 rather than being omitted, since Number=R requires one entry per allele.
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

    // GL over the alleles actually emitted, in VCF genotype order. The spec
    // requires exactly one entry per genotype of the emitted allele set, which is
    // a subset of what we scored, so this is a lookup rather than a rescore --
    // recomputing would mean re-scoring every read.
    bool all_mapped = true;
    for (size_t s = 0; s < traversals.size(); ++s) {
        if (site_to_scored[s] < 0) {
            all_mapped = false;
            break;
        }
    }

    // A genotype carrying a star or missing allele was called at a lower ploidy
    // than the record reports: in nested mode only some parent haplotypes traverse
    // the child, so GT has an entry per parent haplotype while the likelihoods were
    // computed over the traversing ones only. Emitting GL here would give a vector
    // whose length disagrees with the ploidy GT implies, which is worse than
    // omitting it.
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
            gl_strings.push_back(std::to_string(found->second / 2.30258509299));
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

    // The likelihood-ratio quality before the explained-share discount. Emitted
    // unconditionally, including when the discount is off and the two are equal, so a
    // consumer never has to know which mode produced the file to know what GQ means.
    variant.format.push_back("GQI");
    variant.samples[sample_name]["GQI"].push_back(
        std::to_string(min((int)256, max((int)0, (int)info->gq_undiscounted))));

    // The depth- and ploidy-normalised quality. Emitted as a fraction rather than
    // rescaled into phred on purpose: it is a proportion of what this site could have
    // achieved, and dressing it as a phred score would invite exactly the cross-site
    // comparison of absolute quality that it exists to replace. "." where there was no
    // gap to normalise, which is not the same as 0.
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

    // Marked, never dropped, and this is a design decision rather than a convenience.
    // Records are buffered as text and the linkage layer rewrites their genotypes in
    // VCFOutputCaller::write_variants, *after* this runs -- so a record withheld here is
    // withheld before linkage ever sees it, and a low-confidence site is precisely the
    // kind linkage exists to fix. Dropping is one `bcftools view -f PASS` away for anyone
    // who wants it, and cannot be undone by anyone who does not.
    //
    variant.filter = "PASS";
    if (info->n_informative == 0) {
        variant.filter = "noreads";
    } else if (min_confidence > 0.0 && info->gq_fraction >= 0.0
               && info->gq_fraction < min_confidence) {
        // Only where GQN exists. A site with no gap to normalise reports '.', which is not
        // low confidence -- it is no measurement -- and must not be swept up by a threshold
        // as though it were zero.
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
        "alleles fail to explain. Not used by the genotype model, which assumes each "
        "haplotype contributed 1/ploidy of the reads whatever they show\">\n";
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
              "normalise, which is not 0\">\n";
    header += "##FORMAT=<ID=GP,Number=1,Type=Float,Description=\"Genotype Probability, the "
              "log-scaled posterior of the called genotype under a uniform prior\">\n";
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
