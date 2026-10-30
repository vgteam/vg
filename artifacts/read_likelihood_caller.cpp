#include "read_likelihood_caller.hpp"

#include <cstdio>
#include <algorithm>
#include <cmath>
#include <limits>
#include <optional>
#include <set>
#include <sstream>

#include "../src/statistics.hpp"
#include "../src/sorted_vcf_writer.hpp"
#include "../src/utility.hpp"
#include "../src/crash.hpp"

namespace vg {

/// Write a number with a fixed number of digits after the decimal point, as
/// printf's "%.*f" does.
static std::string fixed_decimal(double value, int digits) {
    char buffer[64];
    snprintf(buffer, sizeof(buffer), "%.*f", digits, value);
    return buffer;
}


using namespace std;

ReadLikelihoodSnarlCaller::ReadLikelihoodSnarlCaller(const PathHandleGraph& graph,
                                                     AlleleLikelihoodCalculator& likelihood_calculator)
    : graph(graph), likelihood_calculator(likelihood_calculator) {
}

ReadLikelihoodSnarlCaller::~ReadLikelihoodSnarlCaller() = default;

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
    gq *= info.depth_discount;
    // Both discounts are fractions, so GQ is never above the gap, and so never above GQI.
    crash_unless(gq <= gap);
    return gq;
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

/**
 * Fill the fields of `info` that come from the matrix alone, whatever the ploidy.
 *
 * Those are the read support of each allele (AD), the mean best fit (BL), the read count (DP), the
 * scored walks and the depth context. Neither AD nor BL enters the genotype likelihood.
 */
static void fill_matrix_fields(const AlleleReadLikelihoods& matrix,
                               const vector<Traversal>& traversals,
                               ReadLikelihoodSnarlCaller::ReadLikelihoodCallInfo& info) {
    info.allele_support.assign(matrix.num_alleles(), 0.0);
    double best_ln_total = 0.0;
    size_t best_ln_n = 0;
    for (size_t r = 0; r < matrix.num_reads(); ++r) {
        double best = 0.0;
        for (size_t a = 0; a < matrix.num_alleles(); ++a) {
            best = max(best, matrix.rel(r, a));
        }
        if (best > 0.0) {
            // Split the read's support among the alleles it fits best.
            size_t winners = 0;
            for (size_t a = 0; a < matrix.num_alleles(); ++a) {
                winners += matrix.rel(r, a) >= best ? 1 : 0;
            }
            for (size_t a = 0; a < matrix.num_alleles(); ++a) {
                if (matrix.rel(r, a) >= best) {
                    info.allele_support[a] += 1.0 / (double)winners;
                }
            }
        }
        if (double bl = matrix.best_ln_likelihood(r); std::isfinite(bl)) {
            best_ln_total += bl;
            ++best_ln_n;
        }
    }
    info.mean_best_ln = best_ln_n ? best_ln_total / (double)best_ln_n : 0.0;

    info.n_reads = matrix.num_reads();
    info.scored_traversals.assign(traversals.begin(), traversals.end());
    if (matrix.has_depth_context()) {
        info.depth = matrix.depth_context();
    }
}

/// Make a call info at ploidy `ploidy` that has the fields of `info` that come from the matrix
/// alone.
static unique_ptr<ReadLikelihoodSnarlCaller::ReadLikelihoodCallInfo> copy_matrix_fields(
    const ReadLikelihoodSnarlCaller::ReadLikelihoodCallInfo& info, int ploidy) {
    auto copy = make_unique<ReadLikelihoodSnarlCaller::ReadLikelihoodCallInfo>();
    copy->n_reads = info.n_reads;
    copy->scored_traversals = info.scored_traversals;
    copy->allele_support = info.allele_support;
    copy->mean_best_ln = info.mean_best_ln;
    copy->depth = info.depth;
    copy->ploidy = ploidy;
    return copy;
}

pair<vector<int>, unique_ptr<ReadLikelihoodSnarlCaller::ReadLikelihoodCallInfo>>
ReadLikelihoodSnarlCaller::genotype_at(const MaterializedSite& site, const Ploidies& ploidies,
                                       const vector<nid_t>& fallback_nodes) {
    const vector<Traversal>& traversals = site.alleles.walks;
    const int ref_trav_idx = site.alleles.ref_allele;
    const int ploidy = ploidies.ploidy;
    auto call_info = make_unique<ReadLikelihoodCallInfo>();
    call_info->ploidy = ploidy;

    if (traversals.empty() || ploidy < 1) {
        return make_pair(vector<int>(), std::move(call_info));
    }

    AlleleReadLikelihoods matrix = likelihood_calculator.compute(
        site, fallback_nodes, ploidies.region_ploidy > 0 ? ploidies.region_ploidy : ploidy);
    fill_matrix_fields(matrix, traversals, *call_info);
    // The call info outlives the matrix, so the per-read evidence moves into it. Read phasing and
    // re-genotyping read it in the rounds, and the render builds the anchors from it.
    call_info->anchor_evidence = matrix.take_anchor_evidence();
    call_info->phase_evidence = matrix.take_phase_evidence();

    if (dump_stream != nullptr) {
        stringstream site_name;
        site_name << graph.get_id(site.start) << (graph.get_is_reverse(site.start) ? "-" : "+")
                  << "_" << graph.get_id(site.end) << (graph.get_is_reverse(site.end) ? "-" : "+");
#pragma omp critical (read_likelihood_dump)
        matrix.dump(*dump_stream, site_name.str());
    }

    if (matrix.num_reads() == 0) {
        // No read can tell the alleles apart here, so every genotype is equally likely. Make no
        // call rather than calling the reference.
        return make_pair(vector<int>(), std::move(call_info));
    }

    vector<int> best_genotype =
        call_at_ploidy(matrix, ploidy, ref_trav_idx, traversals.size(), *call_info);
    if (best_genotype.empty()) {
        return make_pair(vector<int>(), std::move(call_info));
    }

    // Score the site at the other ploidy too, from the same matrix. The depth rate is per
    // haplotype of the region, so it does not change with the number of haplotypes crossing the
    // site.
    if (ploidies.also_score_other && traversals.size() > 1) {
        const int other = ploidy == 1 ? 2 : 1;
        unique_ptr<ReadLikelihoodCallInfo> alt = copy_matrix_fields(*call_info, other);
        vector<int> alt_best = call_at_ploidy(matrix, other, ref_trav_idx, traversals.size(), *alt);
        if (!alt_best.empty()) {
            call_info->alt_ploidy_best = alt_best;
            call_info->alt_ploidy_info = std::move(alt);
        }
    }

    return make_pair(best_genotype, std::move(call_info));
}

namespace {

/// The best and runner-up genotypes among the scored ones, and the total likelihood of all of
/// them.
struct GenotypeRanking {
    /// The index of the best genotype among the scored ones.
    size_t best_index = 0;
    /// The index of the runner-up, or the number of scored genotypes if there is none.
    size_t second_best_index = 0;
    /// The ln likelihood of the best genotype.
    double best_ll = 0;
    /// The ln likelihood of the runner-up, or -infinity if there is none.
    double second_best_ll = -numeric_limits<double>::infinity();
    /// The ln of the sum of the likelihoods of all the scored genotypes.
    double total_ll = -numeric_limits<double>::infinity();
};

/**
 * Rank the scored genotypes, starting from `start_best` as the best.
 *
 * A genotype displaces the best only with a strictly higher likelihood, so ties stay with
 * `start_best`.
 */
GenotypeRanking rank_genotypes(const vector<pair<vector<int>, double>>& scored,
                               size_t start_best) {
    GenotypeRanking ranking;
    ranking.best_index = start_best;
    ranking.best_ll = scored[start_best].second;
    ranking.second_best_index = scored.size();
    for (size_t i = 0; i < scored.size(); ++i) {
        const double ll = scored[i].second;
        if (i != ranking.best_index) {
            if (ll > ranking.best_ll) {
                ranking.second_best_ll = ranking.best_ll;
                ranking.second_best_index = ranking.best_index;
                ranking.best_ll = ll;
                ranking.best_index = i;
            } else if (ll > ranking.second_best_ll) {
                ranking.second_best_ll = ll;
                ranking.second_best_index = i;
            }
        }
        ranking.total_ll = (i == 0) ? ll : add_log(ranking.total_ll, ll);
    }
    return ranking;
}

/**
 * Find the all-reference genotype among the scored ones, or else return 0.
 *
 * Exact ties are broken in its favour. When every read fits every allele equally, all genotypes
 * score the same, and taking whichever came first would make a non-reference call on no
 * evidence.
 */
size_t reference_genotype_index(const vector<pair<vector<int>, double>>& scored, int ploidy,
                                int ref_trav_idx, size_t trav_count) {
    if (ref_trav_idx < 0 || (size_t)ref_trav_idx >= trav_count) {
        return 0;
    }
    const vector<int> ref_genotype(ploidy, ref_trav_idx);
    for (size_t i = 0; i < scored.size(); ++i) {
        if (scored[i].first == ref_genotype) {
            return i;
        }
    }
    return 0;
}

} // namespace

vector<int> ReadLikelihoodSnarlCaller::call_at_ploidy(const AlleleReadLikelihoods& matrix,
                                                      int ploidy, int ref_trav_idx,
                                                      size_t trav_count,
                                                      ReadLikelihoodCallInfo& info) const {
    // Score every genotype. For K alleles and ploidy 2 that is K(K+1)/2 genotypes, each one pass
    // over the reads.
    const vector<pair<vector<int>, double>> scored = matrix.score_genotypes(ploidy);
    if (scored.empty()) {
        return vector<int>();
    }
    for (const auto& [genotype, ll] : scored) {
        info.genotype_lls[genotype] = ll;
    }
    const GenotypeRanking ranking =
        rank_genotypes(scored, reference_genotype_index(scored, ploidy, ref_trav_idx, trav_count));
    const bool have_gap = std::isfinite(ranking.best_ll) && std::isfinite(ranking.second_best_ll);

    // GQ is the phred-scaled gap between the best and second best genotype, as in vg's other
    // callers.
    info.gq = 0;
    if (have_gap) {
        info.gq = logprob_to_phred(ranking.second_best_ll) - logprob_to_phred(ranking.best_ll);
    }
    info.gq_undiscounted = info.gq;

    // GQN is the same gap as a fraction of the achievable gap. Compute it from the unclamped
    // log-likelihoods, since GQ is capped at 256 in the VCF. Both are in nats.
    if (ranking.second_best_index < scored.size() && have_gap) {
        info.achievable_gap = matrix.achievable_gap(scored[ranking.best_index].first,
                                                    scored[ranking.second_best_index].first);
        if (info.achievable_gap > 0.0) {
            info.gq_fraction = min(
                1.0, max(0.0, (ranking.best_ll - ranking.second_best_ll) / info.achievable_gap));
        }
    }

    const vector<int>& called = scored[ranking.best_index].first;
    info.explained_share = explained_share(info, called);
    // GQN takes the same discount as GQ, since the likelihood gap cannot see reads that fit an
    // uncalled allele. `share_discount` controls GQ only.
    if (info.gq_fraction >= 0.0) {
        info.gq_fraction *= info.explained_share;
    }
    info.depth_ratio = matrix.depth_ratio(called);
    info.depth_discount = depth_discount_of(matrix, called, ref_trav_idx, info.depth_ratio);
    info.gq = discounted_gq(info, called, info.gq_undiscounted);

    // The uniform prior cancels, so the posterior is the genotype's likelihood over the sum of
    // the likelihoods of all genotypes.
    info.posterior = std::isfinite(ranking.total_ll) ? ranking.best_ll - ranking.total_ll : 0.0;
    return called;
}

double ReadLikelihoodSnarlCaller::depth_discount_of(const AlleleReadLikelihoods& matrix,
                                                    const vector<int>& called, int ref_trav_idx,
                                                    double depth_ratio) const {
    if (depth_quality <= 0.0 || depth_ratio <= 0.0 || ref_trav_idx < 0) {
        return 1.0;
    }
    // Use the lengths without boundary nodes that the expected read count uses, since the VCF
    // alleles do not exist yet.
    const size_t ref_len = matrix.traversal_length((size_t)ref_trav_idx);
    size_t change = 0;
    for (int a : called) {
        if (a < 0) {
            continue;
        }
        const size_t len = matrix.traversal_length((size_t)a);
        change = max(change, len > ref_len ? len - ref_len : ref_len - len);
    }
    if (change < depth_quality_min_length) {
        return 1.0;
    }
    return exp(-depth_quality * fabs(log(depth_ratio)));
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

namespace {

/**
 * Get DR for the written genotype.
 *
 * The written genotype is the direct call's unless the linkage model moved the record. Use it
 * where each of its alleles maps to a scored column and it has the ploidy the site was genotyped
 * at. Otherwise use the direct call's DR.
 *
 * @return DR, or a negative number if it is unknown.
 */
double written_depth_ratio(const ReadLikelihoodSnarlCaller::ReadLikelihoodCallInfo& info,
                           const vector<int>& site_to_scored, const vector<int>& genotype) {
    vector<int> scored_written;
    for (int allele : genotype) {
        if (allele < 0) {
            continue;
        }
        if ((size_t)allele >= site_to_scored.size() || site_to_scored[allele] < 0) {
            return info.depth_ratio;
        }
        scored_written.push_back(site_to_scored[allele]);
    }
    if ((int)scored_written.size() != info.ploidy) {
        return info.depth_ratio;
    }
    // Sort, as the direct call's genotype is, so that the sum runs in the same order.
    std::ranges::sort(scored_written);
    const double written = info.depth.ratio(scored_written);
    return written >= 0.0 ? written : info.depth_ratio;
}

/**
 * Get GL over the written alleles, in VCF genotype order.
 *
 * There is one entry per genotype of the written alleles, looked up in the scored genotypes.
 *
 * @return The GL entries, or none if any genotype was not scored.
 */
vector<string> written_gl(const ReadLikelihoodSnarlCaller::ReadLikelihoodCallInfo& info,
                          const vector<int>& site_to_scored) {
    vector<string> gl_strings;
    for (const auto& site_genotype :
         AlleleReadLikelihoods::enumerate_genotypes(site_to_scored.size(), info.ploidy)) {
        // Translate to matrix columns and sort, since genotype_lls is keyed by the sorted
        // multiset.
        vector<int> scored_genotype;
        scored_genotype.reserve(site_genotype.size());
        for (int site_allele : site_genotype) {
            scored_genotype.push_back(site_to_scored[site_allele]);
        }
        std::ranges::sort(scored_genotype);
        auto found = info.genotype_lls.find(scored_genotype);
        if (found == info.genotype_lls.end()) {
            return {};
        }
        // VCF wants GL log10-scaled.
        gl_strings.push_back(std::to_string(ln_to_log10(found->second)));
    }
    return gl_strings;
}

/**
 * Get QUAL, the phred-scaled probability that the site is not variant.
 *
 * Take it from the posterior of the all-reference genotype, renormalised against all the scored
 * genotypes. Only the column of the reference allele is needed, so a record with an unscored
 * allele, such as a star allele, still gets one.
 *
 * @return QUAL, or 0 for a reference call or where the posterior is unknown.
 */
double written_quality(const ReadLikelihoodSnarlCaller::ReadLikelihoodCallInfo& info,
                       const vector<int>& site_to_scored, const vector<int>& genotype) {
    if (genotype.empty() || site_to_scored.empty() || site_to_scored[0] < 0) {
        return 0;
    }
    vector<int> ref_genotype(info.ploidy, site_to_scored[0]);
    std::ranges::sort(ref_genotype);
    auto found = info.genotype_lls.find(ref_genotype);
    if (found == info.genotype_lls.end()) {
        return 0;
    }
    double total = -numeric_limits<double>::infinity();
    bool first = true;
    for (const auto& entry : info.genotype_lls) {
        total = first ? entry.second : add_log(total, entry.second);
        first = false;
    }
    if (!std::isfinite(total)) {
        return 0;
    }
    const bool is_ref_call = std::ranges::all_of(genotype, [](int a) { return a == 0; });
    return is_ref_call ? 0 : logprob_to_phred(found->second - total);
}

/// Add a FORMAT field with one value to the record.
void add_format(vcflib::Variant& variant, const string& sample_name, const string& key,
                const string& value) {
    variant.format.push_back(key);
    variant.samples[sample_name][key].push_back(value);
}

/// Format a phred quality as the VCF writes it, truncated and clamped to [0, 256].
string printed_quality(double quality) {
    return std::to_string(min((int)256, max((int)0, (int)quality)));
}

} // namespace

void ReadLikelihoodSnarlCaller::write_vcf_fields(const vector<int>& site_to_scored,
                                                 const vector<int>& genotype,
                                                 const SnarlCaller::CallInfo* call_info,
                                                 const string& sample_name,
                                                 vcflib::Variant& variant) const {
    const auto* info = dynamic_cast<const ReadLikelihoodCallInfo*>(call_info);
    if (info == nullptr) {
        // Another genotyper called the site, so there is no matrix to report.
        return;
    }

    add_format(variant, sample_name, "DP", std::to_string(info->n_reads));
    // BL says whether the reads fit any allele, which GQ does not, so a filter can combine the
    // two.
    add_format(variant, sample_name, "BL", fixed_decimal(info->mean_best_ln, 2));
    // Write DR whether or not the depth term is on.
    if (const double depth_ratio = written_depth_ratio(*info, site_to_scored, genotype);
        depth_ratio >= 0.0) {
        add_format(variant, sample_name, "DR", fixed_decimal(depth_ratio, 3));
    }

    // Write AD over the written alleles, rounded to integers. It need not sum to DP, since a read
    // whose best allele was scored but not written has no entry. An allele with no scored column,
    // such as a star allele, gets 0, since AD needs one entry per allele.
    variant.format.push_back("AD");
    for (int k : site_to_scored) {
        long v = 0;
        if (k >= 0 && (size_t)k < info->allele_support.size()) {
            v = lround(info->allele_support[k]);
        }
        variant.samples[sample_name]["AD"].push_back(std::to_string(v));
    }

    // Leave GL out if a written allele was not scored. Also leave it out for a genotype with a
    // star or missing allele, which was called at a lower ploidy than the record shows, so GL
    // would have the wrong length.
    const bool all_mapped = std::ranges::none_of(site_to_scored, [](int k) { return k < 0; });
    const bool genotype_has_marker = std::ranges::any_of(genotype, [](int a) { return a < 0; });
    if (all_mapped && !genotype_has_marker && info->ploidy >= 1) {
        if (vector<string> gl_strings = written_gl(*info, site_to_scored); !gl_strings.empty()) {
            variant.format.push_back("GL");
            for (const string& gl : gl_strings) {
                variant.samples[sample_name]["GL"].push_back(gl);
            }
        }
    }

    add_format(variant, sample_name, "GQ", printed_quality(info->gq));
    // GQI is the gap the site's own likelihoods gave, before either discount. Re-genotyping
    // recomputes GQ, but not GQI.
    add_format(variant, sample_name, "GQI", printed_quality(info->gq_undiscounted));
    // GQN is a fraction rather than a phred score, or "." where there was no gap to normalise.
    add_format(variant, sample_name, "GQN",
               info->gq_fraction < 0.0 ? string(".") : fixed_decimal(info->gq_fraction, 3));
    // GP is in natural log, unlike GL, which the VCF specification puts in log10. The header says
    // so.
    add_format(variant, sample_name, "GP", std::to_string(info->posterior));

    variant.quality = written_quality(*info, site_to_scored, genotype);

    // Mark low-confidence records rather than dropping them, since the linkage model may move
    // their genotypes after this runs.
    variant.filter = "PASS";
    if (info->n_reads == 0) {
        variant.filter = "noreads";
    } else if (min_confidence > 0.0 && info->gq_fraction >= 0.0
               && info->gq_fraction < min_confidence) {
        variant.filter = "lowconf";
    }
}

void ReadLikelihoodSnarlCaller::add_header_lines(SortedVCFWriter& writer) const {
    writer.add_format_line(
        "##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Number of informative reads "
        "overlapping the site\">");
    writer.add_format_line(
        "##FORMAT=<ID=AD,Number=R,Type=Integer,Description=\"Reads whose best-fitting allele is "
        "this one. **AD does not sum to DP**, for two reasons: a read fitting several alleles "
        "equally splits its vote between them, and more importantly only alleles that reached this "
        "record get a column, while the genotyper scored every allele the site offered. At a site "
        "where many alleles were enumerated and few emitted, most reads best-fit something absent "
        "here and the shortfall is large. That shortfall is itself informative: it is how much of "
        "the evidence the emitted alleles fail to explain. Not used by the genotype model, which "
        "gives each haplotype a share of the reads set by the length of sequence unique to its "
        "allele (an equal share under --flat-mixture), whatever the reads show\">");
    writer.add_format_line(
        "##FORMAT=<ID=DR,Number=1,Type=Float,Description=\"Observed reads at this site divided by "
        "the number the called genotype predicts, from the rate of read starts per base in the "
        "site's rate window, the fixed 16,384 bp reference bucket holding the site and one bucket "
        "on each side, and the called alleles' traversal lengths. 1.0 means the read count is "
        "exactly what the call implies. Values well above 1 are collapsed repeats, where reads "
        "from several copies pile onto one; values near 0.5 are a genotype claiming twice the "
        "sequence actually covered, which is what a missed heterozygous deletion looks like. Both "
        "sides of the ratio count a read as 1 - e_r, the probability it came from this locus at "
        "all, so a site whose mapping quality matches its neighbourhood's is unaffected and only "
        "the difference shows; pass --depth-count-raw to count whole reads instead. Reported "
        "whether or not --depth-term is in use\">");
    writer.add_format_line(
        "##FORMAT=<ID=BL,Number=1,Type=Float,Description=\"Mean over reads of the best raw "
        "alignment score any allele gave them. Measures whether reads fit anything at this site, "
        "where GQ measures only the gap between the top two genotypes, so the two are nearly "
        "independent. NOT normalised for site size or read overlap, so it is comparable between "
        "calls at similar sites rather than across a whole genome\">");
    writer.add_format_line(
        "##FORMAT=<ID=GL,Number=G,Type=Float,Description=\"Genotype Likelihood, log10-scaled "
        "P(reads | genotype) from the read-level likelihood model. Useful for ranking; not a "
        "calibrated probability, and over-confident at high depth because reads are treated as "
        "independent\">");
    writer.add_format_line(
        "##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"Genotype Quality, the phred-scaled "
        "gap between the best and second-best genotype likelihood, scaled by the fraction of reads "
        "the called genotype explains (sum(AD)/DP). The likelihood ratio alone cannot see reads "
        "that fit an allele outside the call, because those reads enter every genotype's "
        "likelihood and cancel; the scaling restores them. It also means GQ here is a quality "
        "score rather than a calibrated posterior. --no-share-quality turns this scaling off. With "
        "--depth-quality A in effect, records whose called alleles change length by at least 50 bp "
        "are also scaled by exp(-A * |ln DR|), so a call whose read count is implausible for the "
        "sequence it claims ranks lower. GQI is the value with neither scaling\">");
    writer.add_format_line(
        "##FORMAT=<ID=GQI,Number=1,Type=Integer,Description=\"Genotype Quality from the likelihood "
        "ratio alone, with neither the explained-read-fraction scaling nor the --depth-quality "
        "scaling. Equals GQ under --no-share-quality, except where --depth-quality scales GQ\">");
    writer.add_format_line(
        "##FORMAT=<ID=GQN,Number=1,Type=Float,Description=\"Normalised Genotype Quality: the "
        "likelihood-ratio gap as a fraction, in [0,1], of the gap this site could have produced "
        "had every read fitted the call perfectly. Unlike GQ it means the same thing at any depth "
        "and any ploidy, so one threshold works across a 5x diploid and a 15x haploid contig. GQ "
        "cannot: it scales with depth, and it scales with ploidy the other way, because at ploidy "
        "1 the runner-up genotype is a different allele outright and every read discriminates "
        "fully, where a diploid heterozygote's runner-up differs on one strand only. On HG002 that "
        "makes hemizygous chrX calls run a median GQ of 247 where chr7 diploid homozygotes at the "
        "same depth run 46. Dividing GQ by depth does not fix this and measurably makes it worse, "
        "since it corrects the smaller axis and leaves the larger. '.' where there was no gap to "
        "normalise, which is not 0. NEGATIVE, and so in [-1,1] rather than [0,1], on a record the "
        "linkage layer moved: there GQN is the settled genotype's margin over the best alternative "
        "on the same scale, and a negative value says the panel overrode the reads. Those records "
        "are the caller's highest false-positive density -- 37.8% against 8.6% overall on 44x ONT "
        "chr20 -- so they are reported as a signed number rather than blanked, which is what left "
        "--min-confidence unable to see them\">");
    writer.add_format_line(
        "##FORMAT=<ID=GP,Number=1,Type=Float,Description=\"Genotype Probability, the "
        "natural-log-scaled posterior of the called genotype under a uniform prior. Nats, not "
        "log10: unlike GL, exponentiate with e, so -2.303 means p=0.1\">");
    writer.add_filter_line(
        "##FILTER=<ID=lowconf,Description=\"GQN below --min-confidence: the call used less of the "
        "discrimination this site could offer than required. Unlike a GQ threshold this means the "
        "same thing at any depth and any ploidy -- requiring GQ >= 10 costs a 5x diploid contig a "
        "third of its F1, where GQN >= 0.05 costs it 0.009 and gains a haploid contig 0.017. Marks "
        "rather than drops. Off unless --min-confidence is given\">");
    writer.add_filter_line(
        "##FILTER=<ID=noreads,Description=\"No informative read overlaps the site, so the "
        "read-level model has no evidence either way\">");
}


namespace {

/// Print a field as the record's line prints it: a list joined with commas, or "." when empty.
string printed_field(const vector<string>& list) {
    string text;
    for (size_t i = 0; i < list.size(); ++i) {
        text += (i == 0 ? "" : ",") + list[i];
    }
    return list.empty() ? string(".") : text;
}

/**
 * Get GQ for a genotype the linkage model chose.
 *
 * It is the phred-scaled complement of the posterior, multiplied by the factor the direct call's
 * GQ was, and capped at GQI. The posterior includes the panel's frequency prior, so 1 - posterior
 * can understate the uncertainty where the reads were weakest. The cap keeps GQ within what the
 * reads alone support, and makes GQ <= GQI hold on every rewritten record.
 *
 * @param gqi The printed GQI, or null if the record has none.
 * @return GQ as the VCF prints it, truncated to a whole number and clamped to [0, 256].
 */
string chosen_gq(const LinkageCollector::MovedQuality& moved, const string* gqi) {
    const double posterior = moved.posterior;
    double q = posterior >= 1.0 ? 256.0 : -10.0 * log10(max(1.0 - posterior, 1e-26));
    q *= moved.direct.gq_factor;
    if (gqi != nullptr) {
        try {
            q = min(q, stod(*gqi));
        } catch (const std::exception&) {
            // GQI cannot be read, so leave GQ uncapped.
        }
    }
    // Truncate, as the per-site GQ is, so that equal qualities print the same.
    return std::to_string((int)min(256.0, max(0.0, q)));
}

/// Read the comma-separated numbers of a printed GL into `gl`.
///
/// @return False if any of them cannot be read.
bool parse_gl(const string& printed, vector<double>& gl) {
    size_t start = 0;
    bool parsed = true;
    while (parsed) {
        const size_t comma = printed.find(',', start);
        const string tok =
            printed.substr(start, comma == string::npos ? string::npos : comma - start);
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
    return parsed;
}

/// Read the alleles of a printed GT, leaving out "." alleles.
///
/// A "." allele is left out rather than failing the read, since `1|.` is a nested chain on one
/// strand of a diploid parent, with a real margin.
///
/// @return The alleles, or none if any cannot be read.
vector<int> parse_gt(const string& printed) {
    vector<int> called;
    size_t start = 0;
    while (true) {
        const size_t sep = printed.find_first_of("/|", start);
        const string tok = printed.substr(start, sep == string::npos ? string::npos : sep - start);
        if (tok != "." && !tok.empty()) {
            try {
                called.push_back(std::stoi(tok));
            } catch (const std::exception&) {
                return {};
            }
        }
        if (sep == string::npos) {
            return called;
        }
        start = sep + 1;
    }
}

/**
 * Find the index in GL of the genotype `called`.
 *
 * A diploid GL is indexed j(j+1)/2 + i for i <= j, and a haploid GL by allele. A haploid GL is
 * recognised by its length, n for n alleles against n(n+1)/2. A diploid GL's length is not
 * checked, since the block records of a site share its GL while each has only its own ALTs.
 *
 * @param printed_alt The ALT column, for counting the alleles.
 * @return The index, or `gl_size` if there is none.
 */
size_t gl_index(const vector<int>& called, size_t gl_size, const string& printed_alt) {
    size_t n_alleles = 1;
    if (printed_alt != "." && !printed_alt.empty()) {
        n_alleles += split_delims_keep_empty(printed_alt, ",").size();
    }
    const bool diploid_gl = gl_size == n_alleles * (n_alleles + 1) / 2;
    const bool haploid_gl = gl_size == n_alleles;
    if (called.size() == 2 && gl_size != 0) {
        const int i = min(called[0], called[1]);
        const int j = max(called[0], called[1]);
        return (size_t)(j * (j + 1) / 2 + i);
    }
    if (called.size() == 1 && haploid_gl && !diploid_gl) {
        return (size_t)called[0];
    }
    return gl_size;
}

/**
 * Compute GQN for a genotype the linkage model chose.
 *
 * It is the genotype's likelihood margin over the best alternative, as a fraction of the direct
 * call's achievable gap, times the explained share, as the per-site GQN is. It is negative when
 * the linkage model moved the call against the reads, so the sign is kept.
 *
 * @param achievable_phred The direct call's achievable gap, in phred units.
 * @param gqn Set to GQN, if it can be computed.
 * @return False if the GL or GT cannot be read, or the genotype has no entry in GL.
 */
bool chosen_gqn(const string& printed_gl, const string& printed_gt, const string& printed_alt,
                double achievable_phred, double explained_share, double& gqn) {
    vector<double> gl;
    if (!parse_gl(printed_gl, gl)) {
        return false;
    }
    const size_t idx = gl_index(parse_gt(printed_gt), gl.size(), printed_alt);
    if (idx >= gl.size()) {
        return false;
    }
    double best_other = -std::numeric_limits<double>::infinity();
    for (size_t g = 0; g < gl.size(); ++g) {
        if (g != idx) {
            best_other = max(best_other, gl[g]);
        }
    }
    if (best_other == -std::numeric_limits<double>::infinity()) {
        return false;
    }
    const double margin_phred = 10.0 * (gl[idx] - best_other);
    gqn = min(1.0, max(-1.0, margin_phred / achievable_phred * explained_share));
    return true;
}

} // namespace

bool ReadLikelihoodSnarlCaller::rewrite_quality_for_chosen_genotype(
        vcflib::Variant& variant, const string& sample,
        const LinkageCollector::MovedQuality& moved, double lowconf_threshold) {
    auto sample_found = variant.samples.find(sample);
    if (variant.format.empty() || sample_found == variant.samples.end()
        || sample_found->second.empty()) {
        return false;
    }
    map<string, vector<string>>& sample_fields = sample_found->second;
    // Read each field as the record's line prints it, or as nothing if FORMAT lacks its key. A key
    // in FORMAT with no values prints ".".
    auto printed = [&](const string& key) -> optional<string> {
        if (std::ranges::find(variant.format, key) == variant.format.end()) {
            return std::nullopt;
        }
        auto found = sample_fields.find(key);
        return printed_field(found == sample_fields.end() ? vector<string>() : found->second);
    };
    const optional<string> gq = printed("GQ");
    const optional<string> gqi = printed("GQI");
    const optional<string> gqn = printed("GQN");
    const optional<string> gl = printed("GL");
    const optional<string> gt = printed("GT");
    const string printed_alt = printed_field(variant.alt);
    string filter = variant.filter.empty() ? "." : variant.filter;

    if (gq) {
        sample_fields["GQ"] = {chosen_gq(moved, gqi ? &*gqi : nullptr)};
    }
    // GQN stays "." when it cannot be computed.
    bool gqn_known = false;
    double gqn_new = 0.0;
    // The divisor of GQN, in phred units as the GL margin is.
    const double achievable_phred = 10.0 * moved.direct.achievable_gap / log(10.0);
    if (gqn && gl && gt && achievable_phred > 0.0) {
        gqn_known = chosen_gqn(*gl, *gt, printed_alt, achievable_phred,
                               moved.direct.explained_share, gqn_new);
    }
    if (gqn) {
        sample_fields["GQN"] = {gqn_known ? fixed_decimal(gqn_new, 3) : string(".")};
    }
    // The direct pass set `lowconf` from its own GQN. Decide it again from the recomputed GQN
    // where there is one, and clear it where there is not.
    if (gqn_known && lowconf_threshold > 0.0) {
        filter = gqn_new < lowconf_threshold ? "lowconf" : "PASS";
    } else if (filter == "lowconf") {
        filter = "PASS";
    }
    variant.filter = filter;
    return true;
}

}
