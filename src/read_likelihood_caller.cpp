#include "read_likelihood_caller.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
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
    AlleleReadLikelihoods matrix = likelihood_calculator.compute(snarl, traversals);

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

    for (size_t i = 0; i < scored.size(); ++i) {
        double ll = scored[i].second;
        call_info->genotype_lls[scored[i].first] = ll;

        if (i != best_index) {
            if (ll > best_ll) {
                second_best_ll = best_ll;
                best_ll = ll;
                best_index = i;
            } else if (ll > second_best_ll) {
                second_best_ll = ll;
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

    variant.filter = "PASS";
    if (info->n_informative == 0) {
        variant.filter = "noreads";
    }
}

void ReadLikelihoodSnarlCaller::update_vcf_header(string& header) const {
    header += "##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Number of informative reads "
              "overlapping the site\">\n";
    header += "##FORMAT=<ID=GL,Number=G,Type=Float,Description=\"Genotype Likelihood, "
              "log10-scaled P(reads | genotype) from the read-level likelihood model. Useful for "
              "ranking; not a calibrated probability, and over-confident at high depth because "
              "reads are treated as independent\">\n";
    header += "##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"Genotype Quality, the "
              "phred-scaled gap between the best and second-best genotype likelihood\">\n";
    header += "##FORMAT=<ID=GP,Number=1,Type=Float,Description=\"Genotype Probability, the "
              "log-scaled posterior of the called genotype under a uniform prior\">\n";
    header += "##FILTER=<ID=noreads,Description=\"No informative read overlaps the site, so the "
              "read-level model has no evidence either way\">\n";
}

}
