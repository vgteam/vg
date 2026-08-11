#ifndef VG_LINKAGE_MODEL_HPP_INCLUDED
#define VG_LINKAGE_MODEL_HPP_INCLUDED

/** \file linkage_model.hpp
 * A Li-Stephens layer over per-site genotype likelihoods, so that consecutive calls are
 * judged against which combinations the haplotype panel actually carries.
 */

#include <cstddef>
#include <vector>

namespace vg {

using namespace std;

/**
 * Reweights per-site genotype likelihoods by linkage between adjacent sites.
 *
 * The caller genotypes each snarl independently, so the emitted call set is the concatenation
 * of per-site argmaxes -- a pair of haplotypes free to switch panel haplotype at every site, at
 * no cost. Measured on chr20 against a 33-haplotype panel, the typical adjacent called pair is
 * 1.8x more likely than independence predicts, and pairs that no single panel haplotype carries
 * are nine times commoner where the reads were undecided (2.8% at GQ < 10 against 0.3% at
 * GQ >= 40) than where they were confident. So the information is real, it is concentrated
 * exactly where the per-site likelihood is flat, and a prior acting there is not fighting the
 * reads.
 *
 * Structure follows PanGenie (Ebler et al. 2022): hidden states are pairs of panel haplotypes,
 * transitions are Li-Stephens, and the genotype is the argmax of the posterior summed over the
 * states implying it. The emission is **not** PanGenie's: it is this caller's own
 * `ln P(reads | G)`, which already carries per-read alignment likelihoods, a mismapping term, a
 * length-weighted mixture and a Poisson depth term. This is a transition model bolted onto a
 * richer emission, so the expected gain is correspondingly smaller than PanGenie's headline.
 *
 * Ordered rather than unordered pairs, deliberately. It doubles the state count to
 * `(H+1)^2`, which is nothing at panel sizes of tens, and it makes the transition factorise
 * exactly -- the forward step becomes O((H+1)^2) per site rather than O((H+1)^4). Genotypes are
 * symmetrised when the posterior is summed.
 *
 * Measured offline over the caller's own emitted likelihoods on chr20-34hap: small-variant
 * genotype F1 0.9546 -> 0.9575 and structural-variant F1 0.4655 -> 0.4697 at weight 2, while
 * changing 0.06% of genotypes with GQ >= 40. Small, and larger than nothing.
 */
class LinkageModel {
public:

    struct Params {
        /// Weight on the transition model. Zero disables it: every transition becomes uniform,
        /// the chain is memoryless, and the posterior collapses to the per-site likelihood, so
        /// the caller's behaviour is recovered exactly. One is the modelled rate; above one
        /// linkage tightens. Tempering the switch *probability* as `rho^weight` rather than
        /// scaling a log-transition keeps it a probability at every setting.
        double weight = 0.0;

        /// Probability that a haplotype-sampling block boundary changes which source assembly a
        /// panel haplotype continues from.
        ///
        /// **A per-graph property, not a constant, and zero by default.** In a haplotype-sampled
        /// GBZ the panel members are recombinations of real assemblies chosen in blocks: inside a
        /// block a haplotype is a contiguous piece of one assembly, so linkage there is genuine,
        /// while at a boundary the sampler continues the same haplotype only some of the time.
        /// A full pangenome has no such blocks and wants zero. Getting this wrong in the
        /// permissive direction over-penalises switching inside blocks, where linkage is real.
        double block_switch = 0.0;

        /// Block length in bp for the above. Only the product of `block_switch` and the number of
        /// boundaries crossed is identifiable from aggregate data, so these two are not
        /// separately meaningful unless boundary positions are known.
        double block_length = 10000.0;

        /// Distance scale of the biological term, in bp.
        double scale = 10000.0;

        /// Floor on the switch probability, so that a switch is never impossible.
        double rho_min = 1e-3;

        /// Mass on the wildcard haplotype, which may carry any allele at any site.
        ///
        /// Not optional in practice. A state implies a genotype, so without a wildcard a genotype
        /// no panel haplotype pair can spell is unreachable -- and the graph need not contain the
        /// sample being genotyped. Omitting it makes the model suppress novel alleles, which
        /// presents as a precision improvement while destroying recall.
        double escape = 1e-2;

        /// Weight on the allele-frequency prior the state space implies, from the number of
        /// haplotype pairs spelling each genotype. One keeps it, zero removes it.
        ///
        /// **Zero by default, and that is a correctness choice rather than caution.** When the
        /// panel was produced by haplotype sampling against the reads being genotyped, panel
        /// allele frequency is already conditioned on those reads, so using it as a prior and
        /// then multiplying by the read likelihood counts the same evidence twice. It is
        /// defensible over a full panel that was not selected using this sample.
        double freq_prior = 0.0;

        /// Sites of exact inference per window, and the margin discarded at each end.
        ///
        /// Exact inference over a whole chain would serialise a caller that is otherwise parallel
        /// over snarls, and chains here are chromosome arms. Linkage is spent by 10-30 kb, or
        /// tens of sites, so a window with a generous margin is near-exact; the margin is dropped
        /// so no posterior is read from a position that can see the artificial window edge.
        size_t window = 2000;
        size_t margin = 250;
    };

    /// One site's contribution. Deliberately free of graph and GBWT types: this layer is pure
    /// arithmetic over likelihoods and a panel matrix, which is what makes it testable.
    struct Site {
        /// Reference position, for the distance between sites.
        size_t position = 0;

        /// Number of alleles, so genotype indices can be decoded.
        size_t num_alleles = 0;

        /// ln P(reads | genotype), in VCF genotype order: index of (i,j) with i <= j is
        /// j*(j+1)/2 + i.
        vector<double> genotype_ln_likelihood;

        /// Allele carried by each panel haplotype, or -1 where the haplotype does not traverse
        /// this site. Absence is not the reference allele: a haplotype whose path ends here
        /// carries no allele, and treating that as reference would invent evidence.
        vector<int> haplotype_allele;
    };

    LinkageModel(const Params& params) : params(params) {}

    /// True when the transition model is armed. At zero weight the caller must be bit-for-bit
    /// unchanged, so it is worth asking rather than relying on the arithmetic to be neutral.
    bool active() const { return params.weight > 0.0; }

    /// Posterior over genotypes per site, in the same order as `genotype_ln_likelihood`.
    /// `sites` must be one chain in reference order. Returns an empty vector per site where no
    /// posterior could be formed.
    vector<vector<double>> posteriors(const vector<Site>& sites) const;

    /// VCF diploid genotype ordering: index of the genotype (i,j).
    static size_t genotype_index(size_t i, size_t j) {
        if (i > j) {
            size_t t = i; i = j; j = t;
        }
        return j * (j + 1) / 2 + i;
    }

    /// Per-strand switch probability between two sites `gap` bp apart, after weighting.
    double switch_probability(size_t gap) const;

private:

    /// Exact forward-backward over one window. `out` is filled for the whole window; the caller
    /// keeps only the interior.
    void window_posteriors(const vector<Site>& sites, size_t from, size_t to,
                           vector<vector<double>>& out) const;

    Params params;
};

}

#endif
