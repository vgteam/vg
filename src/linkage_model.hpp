#ifndef VG_LINKAGE_MODEL_HPP_INCLUDED
#define VG_LINKAGE_MODEL_HPP_INCLUDED

/** \file linkage_model.hpp
 * A Li-Stephens layer over per-site genotype likelihoods, so that consecutive calls are
 * judged against which combinations the haplotype panel actually carries.
 *
 * doc/read-likelihood-genotyping.md describes this layer as a user meets it, including how it
 * relates to the per-site objective it consumes. The comments here are the record of *why* each
 * default is what it is -- what was measured, on which panel, and what did not survive being
 * measured. That record belongs with the code that implements it and is not duplicated there.
 */

#include <cstddef>
#include <cstdint>
#include <mutex>
#include <string>
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

        /// Distance scale, in bp.
        ///
        /// There was once a second term beside this one -- a `block_switch` probability meant to
        /// charge extra for crossing a haplotype-sampling block boundary, on the theory that a
        /// sampled panel haplotype switches source assembly there. It is gone, for two reasons
        /// worth recording so that it is not reinvented.
        ///
        /// It was never a second parameter. Without boundary positions the crossing count has to
        /// be smeared as `gap / block_length`, and then
        ///
        ///     1 - rho' = (1-rho_min) exp(-g/scale) (1-block_switch)^(g/block_length)
        ///              = (1-rho_min) exp(-g/scale_eff),
        ///     1/scale_eff = 1/scale + -ln(1-block_switch)/block_length
        ///
        /// exactly. `block_switch = 0.57` over a 10 kb block was `scale = 5423` and nothing else,
        /// so a grid over the two axes was measuring one axis twice.
        ///
        /// The premise does not hold either, which is why it was not given real boundaries
        /// instead. Measured on the chr20 32-haplotype panel against its own subchain partition,
        /// with gap-matched strata and a permutation control, linkage across a boundary is weaker
        /// by 0.008 NMI at gaps under 5 kb (z = 1.1) -- and that range holds essentially every
        /// adjacent call pair. At a boundary the sampler switches to another assembly *in the
        /// same panel*, and human haplotypes agree at most sites, so switching source changes the
        /// allele only where the two assemblies differ. It largely preserves linkage rather than
        /// destroying it.
        ///
        /// The scale itself is nearly flat from 10 kb to 40 kb: about 0.001 of F1, 20 kb weakly
        /// best. It is left at 10 kb because that is what every measurement to date used.
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
        /// haplotype pairs spelling each genotype. Zero removes it, one keeps it as the state
        /// space presents it, and above one amplifies it beyond that -- the mass on a genotype
        /// is scaled by `multiplicity^freq_prior`, so this is an exponent, not a mixing weight,
        /// and nothing about it stops at one.
        ///
        /// **This is the dominant parameter of the model, and it was defaulted to zero on an
        /// argument that did not survive being measured.**
        ///
        /// The argument was that when the panel comes from haplotype sampling against the reads
        /// being genotyped, panel allele frequency is already conditioned on those reads, so
        /// using it as a prior counts the same evidence twice. That is true and it is not a
        /// reason to switch it off. Nothing connects the panel to the *truth set*: sampling reads
        /// k-mer counts, and the benchmark informed neither the graph nor the selection. Reusing
        /// one's own reads is what mapping and calling already do. What double counting can do is
        /// leave `GQ` overconfident while the genotype itself improves -- a claim about
        /// calibration, checkable against observed error rates per `GQ` bin, and separate from
        /// whether to turn this on.
        ///
        /// Measured on chr20 and chr6 against the 34-haplotype panels, crossed against the
        /// transition weight: it improves every variant class at every weight, monotonically,
        /// through 1 and well past it, peaking near 5-8. At the joint optimum (`weight` 2,
        /// `freq_prior` 5) it is most of the total gain -- small-variant genotype F1 +0.0099 on
        /// chr20 and +0.0074 on chr6 against no linkage at all, against +0.0047 and +0.0036 for
        /// the transition model alone. Beyond 8 it inverts: by 12 the prior overwhelms the reads,
        /// SNV F1 falls below the no-linkage baseline and structural-variant recall collapses.
        ///
        /// The effect is almost entirely small indels -- deletions most, insertions next, and
        /// SNVs flat to within 0.0002 across the whole axis. That is where the emission is flat
        /// and the panel has something to add; SNVs are already settled by the reads.
        ///
        /// One caveat that is not retired: every one of these numbers is from a 34-haplotype
        /// panel. Multiplicity is a far coarser statistic over three haplotypes than over
        /// thirty-three, and a large exponent over a count that barely varies is not the same
        /// operation. Measured on the two 4-haplotype graphs it is neither harmful nor useful: at
        /// `weight` 2 the difference between 0 and 5 is under 0.0004 on every class, in both
        /// directions. That is what the mechanism predicts -- with three haplotypes a genotype is
        /// spelled by at most a couple of pairs, so there is barely any multiplicity for an
        /// exponent to act on. Safe there, not helpful there.
        ///
        /// The struct default stays 0, so that a `LinkageModel` built directly is the plain HMM
        /// with no prior -- which is what the unit tests construct and compare against. `vg call`
        /// defaults its flag to 5.
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

        /// ln P(reads | genotype).
        ///
        /// At `ploidy` 2 this is in VCF genotype order: index of (i,j) with i <= j is
        /// j*(j+1)/2 + i. At `ploidy` 1 a genotype *is* an allele, so it is indexed by allele
        /// directly and has `num_alleles` entries.
        vector<double> genotype_ln_likelihood;

        /// 1 or 2. A whole chain shares one ploidy -- it is a property of the contig, not of the
        /// site -- but it is carried here because this struct is what the model is handed.
        ///
        /// Haploid chains are not a degenerate case to be skipped. chrX outside the
        /// pseudoautosomal regions and all of chrY are haploid in a male sample, and before this
        /// existed they were dropped from the linkage pass entirely: no transition model, and no
        /// mosaic, for about 5% of a genome.
        size_t ploidy = 2;

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

    /// One strand's assignment at one site: an index into the panel, or `WILDCARD`.
    struct Phase {
        size_t first = WILDCARD;
        size_t second = WILDCARD;
    };

    /// The wildcard haplotype's index, one past the panel. It carries any allele at any site, so
    /// a strand assigned to it is "explained by nothing in the panel" rather than by a haplotype.
    static constexpr size_t WILDCARD = (size_t)-1;

    /// Most probable path of haplotype pairs through the chain -- the *phasing*, as distinct from
    /// `posteriors()`, which decides each site on its own.
    ///
    /// The two answer different questions and do not agree in general. A sequence of per-site
    /// marginal argmaxes need not be spellable by any single pair of haplotypes: that is exactly
    /// the free-switching this whole layer exists to penalise, and penalising is not forbidding.
    /// So phasing needs max-product, not sum-product, and it needs its own traceback.
    ///
    /// `constraint[t]` is the genotype index the path must spell at site `t`, or `NO_CONSTRAINT`
    /// to leave the site free. Constraining every site to the called genotype is what makes the
    /// emitted phasing agree with the emitted VCF, which is the only reason to prefer this
    /// decoding over the unconstrained one -- unconstrained Viterbi maximises path probability
    /// and will happily contradict the calls.
    ///
    /// Feasible by construction wherever alleles come from the panel: if allele i is carried by
    /// some haplotype and j by another, the pair spelling (i,j) exists. The wildcard covers the
    /// rest, so a path is always returned.
    ///
    /// Diploid only, like the rest of this class: the states are ordered *pairs*.
    vector<Phase> phasing(const vector<Site>& sites, const vector<size_t>& constraint) const;

    /// Leaves a site's genotype unconstrained in `phasing()`.
    static constexpr size_t NO_CONSTRAINT = (size_t)-1;

    /// Posterior over *alleles* per site, for a haploid chain.
    ///
    /// The diploid model's state is a pair, so it cannot express a haploid chain at all: there is
    /// one strand, and a genotype is an allele rather than an unordered pair. That makes the
    /// haploid case structurally simpler, not a special case of the other -- `H+1` states rather
    /// than `(H+1)^2`, an ordinary Li-Stephens chain, and no symmetrisation.
    ///
    /// `sites[t].genotype_ln_likelihood` is read as one entry per allele here.
    vector<vector<double>> haploid_posteriors(const vector<Site>& sites) const;

    /// Most probable path of single haplotypes through a haploid chain, one per site.
    ///
    /// The haploid analogue of `phasing()`. There is no phase to infer on one strand, so what this
    /// produces is only the mosaic: which panel haplotype explains each stretch. That is still
    /// worth having -- it is the whole answer for chrY and for chrX outside the pseudoautosomal
    /// regions.
    vector<size_t> haploid_phasing(const vector<Site>& sites,
                                   const vector<size_t>& constraint) const;

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

    /// Max-product over one window, with an optional pinned state so consecutive windows join
    /// without inventing a switch at the seam. `out` is indexed from `from`.
    void window_phasing(const vector<Site>& sites, size_t from, size_t to,
                        const vector<size_t>& constraint,
                        size_t pin_index, const Phase& pin, vector<Phase>& out) const;

    /// Emission over single haplotypes for a haploid site: `e[a]` is the relative likelihood of
    /// the allele haplotype `a` carries, with the wildcard last.
    void haploid_emission(const Site& site, size_t n_hap, vector<double>& e,
                          vector<double>& per_allele) const;

    /// Forward-backward and max-product over one window of a haploid chain.
    void window_haploid_posteriors(const vector<Site>& sites, size_t from, size_t to,
                                   vector<vector<double>>& out) const;
    void window_haploid_phasing(const vector<Site>& sites, size_t from, size_t to,
                                const vector<size_t>& constraint,
                                size_t pin_index, size_t pin, vector<size_t>& out) const;

    Params params;
};

/**
 * Collects one compact record per site during calling, then re-decides genotypes once calling is
 * done.
 *
 * Two phases, because the alternative does not fit. Deferring emission until a chain completes
 * would mean holding a `ReadLikelihoodCallInfo` per site -- and that carries a
 * `vector<SnarlTraversal>` and a `map<vector<int>, double>`, which over a chromosome arm runs to
 * hundreds of megabytes and far worse at multi-allelic sites. Keeping only what the HMM consumes
 * costs about 80 bytes a site: roughly 8 MB for chr20, against ~8 GB peak for the run.
 *
 * The other reason is that per-window inference would not work here even though the caller
 * already partitions into windows and parallelises over them. Those windows are 256 node IDs
 * wide, which on chr20 is about eleven sites -- *shorter* than the 10-30 kb over which linkage
 * carries anything. Inference has to span more than one, and having threads genotype their
 * neighbours' sites to get a margin would duplicate the expensive half of the work.
 *
 * So: phase one genotypes exactly as before, in parallel, and records a compact site. Phase two
 * runs between calling and writing, groups by contig, sorts by reference position -- node-ID order
 * is close to reference order but not guaranteed, and transition probabilities are computed from
 * the gaps, so trusting it would silently use the wrong distances -- and reports the genotypes
 * that changed. Records are already buffered whole-genome as compressed strings until
 * `write_variants`, so patching between the two adds no new peak.
 *
 * Measured, chr20-34hap on 5 threads: 105251 sites, 7.87 MB retained -- so the 8 MB above was
 * right -- and 1.8 s in phase two against a total feature cost of 20.2 s, or 16% of the run.
 *
 * **The cost is not where this comment spends its words.** Phase two, the serial pass this design
 * is built around avoiding, is 9% of the overhead. The other 91% is phase one, and within it
 * almost all of that is `VCFOutputCaller::panel_alleles` asking the GBWT which haplotypes take
 * each allele -- three and a half times the cost of the per-read likelihood inner loop it exists
 * to annotate. Caching GBWT records per thread took 28% off; what remains is `locate` itself.
 * `record()`, the global mutex that looks like the obvious contention point, does not register in
 * a profile at all: one sample before the cache, three after.
 */
class LinkageCollector {
public:

    LinkageCollector(const LinkageModel::Params& params, size_t num_haplotypes)
        : params(params), model(params), n_haplotypes(num_haplotypes) {}

    /// A genotype the linkage pass wants changed, identified by the key the caller supplied.
    struct Change {
        size_t record_key = 0;
        /// Where to apply it. The output buffer already stores (contig, position) uncompressed as
        /// its sort key, so a change can be matched without keeping the record itself.
        string contig;
        size_t position = 0;
        /// The genotype the per-site model chose. Checked before patching: two records can share a
        /// position, and patching the wrong one would be silent.
        size_t called_i = 0;
        size_t called_j = 0;
        size_t allele_i = 0;
        size_t allele_j = 0;
        double posterior = 0.0;
        /// The site's explained-read share, carried through so the rewritten GQ can be
        /// discounted the same way the per-site GQ was. Without it a changed record
        /// reports an undiscounted quality while every other record reports a discounted
        /// one, and `GQ <= GQI` -- true everywhere else -- fails on about 5% of records.
        double explained_share = 1.0;
    };

    /// Record one site. Safe to call from several threads. `haplotype_allele` must have one entry
    /// per panel haplotype, -1 where the haplotype does not traverse the site; `called_i/j` is the
    /// genotype the per-site model chose, so a change can be detected without keeping the record.
    void record(const string& contig, size_t position, size_t num_alleles,
                const vector<double>& genotype_ln_likelihood,
                const vector<int>& haplotype_allele,
                size_t called_i, size_t called_j, size_t record_key,
                double explained_share, size_t ploidy = 2,
                int64_t start_node = 0, int64_t end_node = 0,
                bool nested = false, size_t parent_record_key = 0, size_t parent_slot = 0);

    /// One site's phasing: which strand carries which allele, and which panel haplotype explains
    /// each strand.
    ///
    /// The allele pair is ordered, and that order *is* the phase -- `allele_first` sits on the
    /// same strand as every other `allele_first` in the same phase set. The haplotype indices are
    /// what the mosaic output serialises; `LinkageModel::WILDCARD` means the panel does not
    /// explain that strand here.
    struct PhaseCall {
        size_t record_key = 0;
        string contig;
        size_t position = 0;
        size_t allele_first = 0;
        size_t allele_second = 0;
        size_t hap_first = LinkageModel::WILDCARD;
        size_t hap_second = LinkageModel::WILDCARD;
        /// 1 or 2. At 1 only the `_first` fields are meaningful: there is one strand.
        size_t ploidy = 2;
        /// The site's snarl boundary nodes. The mosaic output anchors on these rather than on
        /// reference positions, because a node ID is intrinsic to the graph while a position is
        /// a statement about one reference path.
        int64_t start_node = 0;
        int64_t end_node = 0;
        /// Identifies the phase block. Phase is only comparable within one.
        size_t phase_set = 0;
    };

    /// Run the model per contig and return only the genotypes that changed.
    ///
    /// With `phasing_out`, also returns a phasing of the **final** call set -- the genotypes after
    /// any change above has been applied, not before. Constraining to the final calls is what
    /// makes the phasing agree with the VCF; constraining to the pre-linkage calls would phase a
    /// genotype set that is never emitted.
    vector<Change> resolve(vector<PhaseCall>* phasing_out = nullptr) const;

    /// Retained bytes, for reporting. The point of the compact form is that this stays small, so
    /// it is worth being able to say what it actually is rather than trusting the estimate.
    size_t bytes() const;

    size_t num_sites() const { return entries.size(); }

private:

    struct Entry {
        uint32_t position = 0;
        uint32_t contig = 0;
        uint32_t gl_offset = 0;
        uint32_t hap_offset = 0;
        uint16_t num_alleles = 0;
        uint16_t called_i = 0;
        uint16_t called_j = 0;
        uint8_t ploidy = 2;
        /// True when this site's ploidy came from nested descent rather than from the contig or a
        /// --ploidy-bed region.
        ///
        /// The distinction decides whether a chain is cut here. A regional ploidy change -- chrX's
        /// pseudoautosomal boundary -- must cut: there is no haplotype correspondence across it. A
        /// nested ploidy-1 site must not: the surrounding diploid phase is continuous, and cutting
        /// there discards it. Treating the two alike fragmented chr20's autosomal phasing from 22
        /// blocks to 9,460, collapsing block N50 from 248 Mb to 1.08 Mb.
        bool nested = false;
        /// For a nested site, the parent record it hangs off and which of the parent's genotype
        /// slots crosses this child. A nested site has ploidy 1 precisely because exactly one parent
        /// allele crosses it, so the strand it belongs to is determined -- not a best fit -- once the
        /// parent has been phased.
        size_t parent_record_key = 0;
        uint8_t parent_slot = 0;
        float explained_share = 1.0f;
        size_t record_key = 0;
        /// Snarl boundary nodes, for the mosaic output's anchors. Costs 16 bytes a site, which
        /// `bytes()` reports rather than leaving to arithmetic.
        int64_t start_node = 0;
        int64_t end_node = 0;
    };

    LinkageModel::Params params;
    LinkageModel model;
    size_t n_haplotypes;

    /// Flat arenas rather than a vector per site: at 80 bytes a site the 48 bytes of overhead two
    /// empty vectors would add is most of the budget.
    vector<Entry> entries;
    vector<float> gl_arena;
    vector<int8_t> hap_arena;
    vector<string> contig_names;
    mutable std::mutex mutex;
};

}

#endif
