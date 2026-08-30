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

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <map>
#include <mutex>
#include <string>
#include <unordered_set>
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


        /// This site has no reference position at all, so `position` is not a coordinate and must
        /// not be differenced.
        ///
        /// Not inferred from `position == 0`. A position of zero is a legitimate coordinate at the
        /// head of a contig, and inferring absence from it is the same class of mistake as reading a
        /// crossing mask's 0 as "no allele crosses". Nothing sets this yet: a nested chain reaches
        /// the layer only where the reference crosses it, so every site is positioned today, and it
        /// exists for the population a covering reference admits.
        bool unpositioned = false;

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

        /// Fix this site's haplotype pair in `phasing()` rather than letting the path choose it.
        ///
        /// A generation-wise resolve decodes a later generation's phase against earlier
        /// generations that are already emitted, and an unpinned Viterbi is free to come back with
        /// a settled site's strands swapped -- which would write the new sites' phase in a frame
        /// the VCF does not use. This is the same mechanism the window seams already rely on,
        /// applied to every settled site instead of to one index per seam.
        ///
        /// `(size_t)-1` is `WILDCARD`, spelled out because the constant is declared further down.
        bool pinned = false;
        size_t pin_first = (size_t)-1;
        size_t pin_second = (size_t)-1;
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

    /// Posteriors over one segment of a chain, given the messages reaching its two ends. With the
    /// true alpha and beta a segment decode agrees with the whole-chain decode: that is the Markov
    /// property, and it holds for sum-product as much as for max-product.
    void segment_posteriors(const vector<Site>& sites, size_t from, size_t to,
                            const vector<double>* alpha_in, const vector<double>* beta_in,
                            vector<vector<double>>& out) const {
        // Sized here, not left to the caller: `window_posteriors` indexes `out[from + t]` with no
        // bounds check, and the first public caller that forgot segfaulted immediately. Grown, not
        // assigned, so several segments can be decoded into one buffer.
        if (out.size() < sites.size()) {
            out.resize(sites.size());
        }
        window_posteriors(sites, from, to, out, alpha_in, beta_in);
    }

private:

    /// Exact forward-backward over one window. `out` is filled for the whole window; the caller
    /// keeps only the interior.
    ///
    /// `alpha_in` / `beta_in` are the messages over haplotype PAIRS entering the window's two ends,
    /// m*m entries each in the same (a * m + b) layout as the emissions; uniform when null, which is
    /// what every existing caller wants.
    ///
    /// A message was once the only boundary object thought sufficient for cutting a chain, on the
    /// reasoning that a settled genotype leaves every pair spelling it live (hundreds of states at
    /// 34 haplotypes), a settled PAIR discards the mass on alternatives that are still plausible,
    /// and uniform discards everything. The middle of those three turned out to be wrong where it
    /// matters: conditioning a nested chain on a delta at its parent's settled pair, rather than on
    /// the parent's posterior over pairs, is byte-identical on chr6 and identical on every scored
    /// metric on chr20. So the harvest that produced those posteriors is gone, and `alpha_in` is now
    /// only ever a delta or null. The reasoning still stands for a WINDOW seam inside one chain,
    /// which is where it was measured and where `alpha_in` still carries a real distribution.
    void window_posteriors(const vector<Site>& sites, size_t from, size_t to,
                           vector<vector<double>>& out,
                           const vector<double>* alpha_in = nullptr,
                           const vector<double>* beta_in = nullptr) const;

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
/// One step of the forward/backward transition.
///
/// Declared here so the arithmetic can be unit-tested directly. `viterbi_step` takes the same pair
/// but stays in the .cpp's anonymous namespace with the top-2 helpers it depends on; its four
/// candidates are the four stay/jump combinations, so it is the same generalisation.
///
/// Takes a switch probability per
/// strand: the two haplotypes of a diploid sample recombine independently, so the distance each has
/// travelled since the previous site is its own. Passing the same value twice reproduces the single
/// scalar these used to take, up to floating-point re-association -- see the unit tests, which assert
/// agreement to 1e-12 rather than bit-identity for exactly that reason.
void transition_apply(const std::vector<double>& in, size_t m, double rho_a, double rho_b,
                      std::vector<double>& out);

class LinkageCollector {
public:

    LinkageCollector(const LinkageModel::Params& params, size_t num_haplotypes)
        : params(params), model(params), n_haplotypes(num_haplotypes) {}

    /// A genotype the linkage pass wants changed, identified by the key the caller supplied.

    /// Where a site sits in the tree, as one named object rather than ten trailing parameters.
    ///
    /// The parameters it replaces were all bool, int and size_t, and those convert to one another
    /// silently -- so a call site that dropped one still compiled, with every later argument
    /// shifted by one. That happened three times on this branch and the compiler caught none of
    /// them; designated initialisers make the same mistake a compile error, and a reader can tell
    /// what an argument means without counting commas.
    struct SiteContext {
        /// One copy of the chain, so its likelihoods are haploid.
        bool nested = false;
        size_t parent_record_key = 0;
        /// Bit t set iff the parent's candidate traversal t crosses this chain.
        uint64_t parent_crossing = 0;
        size_t generation = 0;
        /// Whether a VCF line was written for this site. Recorded either way.
        bool emitted = true;
        /// `position` is an anchor in the contig, not a coordinate for this snarl.
        bool unpositioned = false;
        /// The chain's own boundary pair, hashed. Its identity, and the group key.
        size_t chain_key = 0;
    };

    /// Record one genotyped site. Safe to call from several threads.
    ///
    /// Everything here is in the genotyper's own space -- candidate traversal indices -- not in VCF
    /// allele numbering. `haplotype_traversal` has one entry per panel haplotype, the candidate
    /// traversal it carries or -1 where it does not traverse the site. `called_trav_i/j` is the pair
    /// the per-site model chose, so a change can be detected without keeping the record.
    /// `traversal_to_allele` maps candidate traversals to the VCF alleles they were emitted as, for
    /// rendering only; pass it empty for a site that wrote no line, with `ctx.emitted` false.
    ///
    /// A site that emitted nothing is still recorded, because it still has an allele pair and that
    /// pair is what phases its children. Gating entry on "a line was written" is what left a
    /// symbolically-collapsed parent unphased and its children strandless.
    void record(const string& contig, size_t position,
                const map<vector<int>, double>& genotype_ln_likelihood,
                const vector<int>& haplotype_traversal,
                int called_trav_i, int called_trav_j,
                const vector<int>& traversal_to_allele,
                size_t record_key,
                double explained_share, size_t ploidy,
                int64_t start_node, int64_t end_node,
                const SiteContext& ctx);

    /// The compact allele space `record` builds for one site: the called pair plus every traversal
    /// some panel haplotype carries, deduplicated by traversal and sorted.
    ///
    /// Exposed for the unit tests, which is the only way to check the construction without a graph.
    /// `genotype_ln_likelihood` is keyed by *sorted* candidate-traversal vectors, as the genotyper
    /// produces it.
    static vector<int> compact_allele_space(const map<vector<int>, double>& genotype_ln_likelihood,
                                            const vector<int>& haplotype_traversal,
                                            int called_trav_i, int called_trav_j);

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
        /// The candidate traversal on each strand -- the same pair as `allele_*`, in the genotyper's
        /// numbering rather than the collector's compact one.
        ///
        /// This is the genome fact: a crossing mask is expressed in traversal terms, so a child's
        /// strand has to be derived against these and not against the compact indices, which agree
        /// only when every allele at the site is panel-carried. It is also what a per-haplotype path
        /// through the snarl is made of.
        int trav_first = -1;
        int trav_second = -1;
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
        /// True when nothing ordered the allele pair: the site is heterozygous and no panel
        /// haplotype on either strand spells either called allele, so `allele_first` is simply the
        /// smaller index and its pairing with `hap_first` is an accident.
        ///
        /// Which of the parent's two strands this site sits on, for a nested haploid site whose
        /// strand is determined; -1 for everything else -- a haploid contig, a nested site whose
        /// parent could not be found, and the incoherent ones that sit on both strands or neither.
        ///
        /// A nested site is haploid because the parent's other allele deletes the chain, so it is one
        /// strand of a diploid locus rather than a haploid locus. That is the difference between it
        /// and chrY, and the difference the VCF has to carry: a bare `GT` of `1` names no strand, so
        /// the strand existed only in the mosaic and no phasing tool could see it.
        int8_t nested_strand = -1;
        /// Matters to anything that reads the pair *as* the phase, which is what a phased `GT`
        /// claims: at these sites the emitted orientation is a placeholder and is indistinguishable
        /// from one the panel actually chose. It is also the leading explanation for why deriving a
        /// nested site's strand from the parent's phased pair measured worse than the traversal-order
        /// slot recorded at descent -- see the Stage 7 notes in the companion evaluation repository.
        bool order_arbitrary = false;

        /// How deep in the snarl tree this site sits: 0 for a top-level chain, 1 for a chain nested
        /// directly inside one, and so on.
        ///
        /// The mosaic carries it so that a recombination INSIDE a nested chain can be told from one
        /// between top-level sites. Both already appear as a segment boundary -- the haplotype
        /// changes either way -- but a consumer reading the file has no way to know which it is
        /// looking at, and the nested ones are the whole point of nested calling.
        uint8_t depth = 0;

        /// Whether a VCF line exists for this site. Sites that wrote none are still phased -- a
        /// parent whose alleles differ only inside its children has a real pair of haplotypes, and
        /// its children need to know which is which -- but they are not records, so nothing that
        /// counts or patches records may include them.
        bool emitted = true;
    };

    /// A nested site whose ploidy the final parent genotype does not support.
    ///
    // There was a `NestedIncoherence` here, with three kinds and three FILTERs, for the case where
    // the ploidy a nested child was called at disagreed with what its parent's final genotype
    // implies. It is gone, and not because the disagreement was decided to be acceptable: the child
    // now takes its ploidy and its strand from one reading of the parent's settled pair -- which of
    // that pair's traversals carries the chain -- so having one copy and sitting on that traversal's
    // strand are the same statement. There is nothing left for two derivations to disagree about.

    /// Run the model per contig and return only the genotypes that changed.
    ///
    /// With `phasing_out`, also returns a phasing of the **final** call set -- the genotypes after
    /// any change above has been applied, not before. Constraining to the final calls is what
    /// makes the phasing agree with the VCF; constraining to the pre-linkage calls would phase a
    /// genotype set that is never emitted.
    ///
    size_t resolve(vector<PhaseCall>* phasing_out = nullptr) {
        return resolve_generation(0, true, phasing_out);
    }

    /// Resolve one generation of sites, holding every earlier generation fixed.
    ///
    /// Sites of a later generation are not considered at all -- they do not exist yet, because the
    /// caller has not descended into them. Sites of an *earlier* generation are included in the
    /// chains but **clamped**: their emission becomes a delta at the genotype they were settled
    /// at, and their phase is pinned to the pair that was emitted for them. They therefore still
    /// carry transition context for this generation's sites -- which is the whole reason to include
    /// them, since a generation on its own is far too sparse for a 10 kb decay -- while being unable
    /// to move. That is what makes a parent's genotype final before its children are called, and it
    /// is the greedy step in `--nested-after-linkage`.
    ///
    /// Each site's `Change` and `PhaseCall` are produced exactly once, at its own generation.
    /// `phasing_out` accumulates across calls and must be passed back in each time: a nested site's
    /// strand is read off the parent's `PhaseCall`, and by this generation the parent belongs to an
    /// earlier one.
    ///
    /// `last` gates the reporting, the mosaic-facing sort of `phasing_out`, and the coherence
    /// counters, so a run of several generations reports once rather than once per generation.
    ///
    /// With every entry at generation 0 -- which is every run without `--nested-after-linkage` --
    /// this is exactly the single pass it replaces, since nothing is ever clamped or held back.
    /// Posterior and explained-read share for each site the model moved, by record key.
    ///
    /// The quality of a moved genotype is not derivable from the read likelihoods alone -- it is the
    /// phred complement of the HMM posterior, discounted by the explained share and capped at GQI --
    /// and the posterior exists only here. Exposed rather than pushed into a patch because the record
    /// is built after the decision, so whoever builds it can ask.
    const std::unordered_map<size_t, std::pair<double, double>>& moved_quality() const {
        return moved_quality_by_record;
    }

    /// Returns how many sites the model moved off their called genotype.
    size_t resolve_generation(size_t generation, bool last,
                                      vector<PhaseCall>* phasing_out = nullptr);

    /// Drop a site before its generation resolves: the settled parent genotype does not carry the
    /// chain at all, so the sample has no copy of it and no call belongs there.
    ///
    /// The entry is neutralised rather than erased -- the arenas are flat and offsets into them are
    /// held by every other entry -- so it is marked and then skipped by chain construction, phasing
    /// and every counter. Returns false for an unknown key.
    bool retract(size_t record_key);


    /// The pair of candidate traversals this site settled on, for a caller that builds the record
    /// after the decision instead of patching one written before it.
    ///
    /// Decoded out of the compact space here, because a compact index means nothing outside the
    /// collector -- that conflation is what put allele numbers past the end of an ALT list earlier on
    /// this branch. `final_i`/`final_j` are initialised to the called pair for every entry at the top
    /// of its generation's resolve, so an unmoved site returns what it came in with rather than
    /// nothing. Returns false for an unknown or retracted key.
    bool settled_traversals(size_t record_key, int* first, int* second, size_t* ploidy) const;

    /// Fill in the traversal-to-VCF-allele map for a site already recorded, and say whether a line
    /// exists for it.
    ///
    /// `record()` moved to the genotyping site, which is before any allele list exists: the emitted
    /// alleles are chosen while the record is built, and the map is a function of that choice. So the
    /// entry is created with no map and the emitter supplies one when it writes the line. Both halves
    /// are needed only by the patch path, and both go away with it.
    ///
    /// Returns false for an unknown key.
    bool set_allele_map(size_t record_key, const vector<int>& traversal_to_allele, bool emitted);

    /// The keys of every record that ended up with a VCF line.
    ///
    /// Read live rather than off a PhaseCall, because the phasing vector holds copies taken while
    /// the genotypes were being resolved and a record's line is built after that: the copy's
    /// `emitted` is a snapshot that is false for every site. Returned as a set in one pass over the
    /// entries rather than answered per key, which would be a scan per call over a quarter of a
    /// million entries.
    std::unordered_set<size_t> emitted_records() const;


    /// What a settled parent implies about one of its children: how many copies the sample carries,
    /// which of the parent's settled traversals carries it when that is one, and whether the
    /// question could be answered at all.
    ///
    /// DERIVED, never stored. `Entry::parent_trav` used to hold the carrying traversal and
    /// `Entry::ploidy` the copy count -- the same fact twice, written by a barrier pass with five
    /// early exits, so they could go stale and could go stale INCONSISTENTLY WITH EACH OTHER:
    /// entries with ploidy 2 and parent_trav >= 0 were measured, which the barrier's own comment
    /// says cannot happen. A value computed where it is consumed cannot do that.
    ///
    /// Landed first as a check (19,979 derivations on chr20, 52,800 over three contigs, zero
    /// disagreements), then consumed, then the stored field deleted -- in that order, because a
    /// check that does not exercise the substituted call site validates something adjacent to what
    /// is relied on.
    struct Relation {
        uint8_t copies = 0;
        int carrying_trav = -1;   ///< the traversal when copies == 1; -2 when both carry it
        bool known = false;       ///< false when the parent's settled pair could not be read
    };


    /// Whether an active (non-retracted) entry exists for this key. The barrier needs to tell
    /// "this site was never recorded" from "respecify refused a site that is already in the layer":
    /// in the second case the recorded entry describes a line the barrier has just replaced, so
    /// leaving it in place would patch the new line with the old line's allele numbering.
    bool has_entry(size_t record_key) const;

    /// How many sites belong to one generation, for reporting a per-generation pass.
    size_t num_sites_at(size_t generation) const;

    /// The highest generation any recorded site belongs to.
    size_t max_generation() const;

    /// Retained bytes, for reporting. The point of the compact form is that this stays small, so
    /// it is worth being able to say what it actually is rather than trusting the estimate.
    size_t bytes() const;

    size_t num_sites() const { return entries.size(); }

private:

    /// No entry. Chain terminator for `next_same_key` and the miss value of `live_index`.
    static constexpr uint32_t NO_ENTRY = (uint32_t)-1;

    struct Entry {
        uint32_t position = 0;
        /// Identity of this snarl's chain, from its boundary pair. The group key: sibling chains
        /// have no transition between them, so a chain needs only to be distinguishable from its
        /// siblings -- which the graph answers directly, with no alignment.
        size_t chain_key = 0;
        uint32_t contig = 0;
        uint32_t gl_offset = 0;
        uint32_t hap_offset = 0;
        /// Compact allele -> the candidate traversal it is, and -> the VCF allele it was emitted as.
        ///
        /// The site's allele space is a compact set of *distinct traversals*: the called pair plus
        /// every traversal some panel haplotype carries. That is the genotyper's own space, and it is
        /// the one the model has to work in, because symbolic collapsing maps distinct traversals
        /// onto one VCF allele -- so a parent whose haplotypes differ only inside its child chains is
        /// homozygous in VCF numbering and heterozygous in this one. Only the latter can phase its
        /// children.
        ///
        /// `trav_offset` is the genome fact: it answers "which traversal is this strand on", which is
        /// what a crossing mask tests and what a haplotype path needs. `allele_offset` is presentation
        /// only: -1 where a traversal reached the model but no VCF allele was emitted for it, which is
        /// possible now that the two spaces are not the same.
        uint32_t trav_offset = 0;
        uint32_t allele_offset = 0;
        uint16_t num_alleles = 0;
        uint16_t called_i = 0;
        uint16_t called_j = 0;
        uint8_t ploidy = 2;
        /// Which resolve pass settles this site: 0 for a top-level site and for anything descended
        /// inline behind a parent linkage cannot move, k for a child deferred behind k barriers.
        /// Every entry is 0 without --nested-after-linkage, so the default path resolves in one pass.
        uint8_t generation = 0;
        /// Dropped by `retract`: the settled parent does not carry this chain. Kept in place because
        /// the arenas are flat, and skipped everywhere a site is considered.
        bool retracted = false;
        /// The genotype this site was settled at, written back when its own generation resolves so
        /// that later generations can clamp it. Meaningless before that.
        uint16_t final_i = 0;
        uint16_t final_j = 0;
        /// Whether this site wrote a VCF line. A genotyped snarl enters the layer either way -- it
        /// has an allele pair, which is what phasing needs -- but only an emitted one has a record to
        /// patch, and only its `allele_offset` entries mean anything.
        bool emitted = true;
        /// This site has no reference position: `position` is not a coordinate and must not be
        /// differenced. Set only for a chain the reference does not cross, which reaches the layer
        /// only under --nested-no-ref.
        bool unpositioned = false;
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
        /// One bit per parent *candidate traversal*, set where that traversal crosses this child
        /// chain; 0 means descent could not say.
        ///
        /// Kept for the descent decision -- whether any traversal the parent could settle on reaches
        /// this chain -- and no longer consulted when the child is phased. It used to be re-tested
        /// there against the settled pair to check the child's ploidy, which was a second derivation
        /// of a fact the settled pair already carries, and the two could disagree.
        ///
        /// Placed next to the key rather than beside the narrow members: after an 8-byte member it
        /// needs no padding, where after a `uint8_t` it costs seven bytes of it -- 16 bytes a site
        /// instead of 8, which `bytes()` reported as 1.8 MB on chr20 rather than 0.9.
        uint64_t parent_crossing = 0;
        float explained_share = 1.0f;
        /// The next entry sharing this record key, or NO_ENTRY.
        ///
        /// Every by-key accessor used to scan the whole vector for the first live match. At
        /// 219 k entries, one `settled_traversals` and one `set_allele_map` per rendered
        /// record, and a global mutex held across the scan, that was quadratic and
        /// serialising at once: 31% of a chr20 run, with four threads in `__psynch_mutexwait`
        /// while the fifth walked 17 MB. The chain keeps first-match-in-insertion-order
        /// exactly, including duplicate keys, and fits in the padding before `record_key`.
        uint32_t next_same_key = NO_ENTRY;
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
    /// Per compact allele, the candidate traversal index it stands for. uint16 because a site's
    /// candidate list is the traversal finder's output, not the panel size.
    vector<uint16_t> trav_arena;
    /// Per compact allele, the VCF allele it was emitted as, or -1 for none.
    vector<int8_t> allele_arena;
    vector<string> contig_names;
    /// record key -> first and last entry carrying it, so a lookup is a hash probe and a
    /// walk of that key's chain. `last` is what makes appending O(1) rather than a walk.
    std::unordered_map<size_t, uint32_t> first_by_key;
    std::unordered_map<size_t, uint32_t> last_by_key;


    /// record key -> (posterior of the settled genotype, explained-read share).
    std::unordered_map<size_t, std::pair<double, double>> moved_quality_by_record;

    /// The first non-retracted entry with this key, or NO_ENTRY. Call with `mutex` held.
    uint32_t live_index(size_t record_key) const;

    /// Split a settled compact pair into the traversal and the VCF allele on each strand.
    void finish_phase_call(PhaseCall& pc, const Entry& e, size_t& fallbacks) const;

    /// The compact allele space a site is described in, with the genotype likelihoods folded into
    /// it. `record` and `respecify` built this identically -- 42 shared lines, 35 of them in two
    /// contiguous runs -- because a revised site has to be described exactly as a freshly recorded
    /// one would be. Now they cannot disagree about it.
    struct CompactSite {
        /// Compact allele -> candidate traversal, sorted by candidate index, so the numbering is a
        /// property of the site rather than of the order the genotypes arrived in.
        vector<int> space;
        /// Genotype likelihoods over the compact space, in `LinkageModel::genotype_index` layout.
        vector<float> gls;
        /// The called pair, compacted. Both are >= 0 whenever `ok`.
        int ci = -1, cj = -1;
        size_t site_ploidy = 0;
        bool ok = false;

        /// Candidate traversal -> compact allele, or -1. `space` is sorted, so this is a search
        /// rather than the `map<int, int>` both callers used to build and throw away -- one map
        /// and up to 127 node allocations per site, on every site of the contig.
        int compact_of(int trav) const {
            auto it = std::lower_bound(space.begin(), space.end(), trav);
            return (it == space.end() || *it != trav) ? -1 : (int)(it - space.begin());
        }
    };

    /// Build the compact space and fold the likelihoods into it. `ok` is false where the site
    /// cannot be described: no candidates, more than the 127 an int8 arena can name, or a called
    /// traversal that is not in its own space.
    CompactSite compact_site(const map<vector<int>, double>& genotype_ln_likelihood,
                             const vector<int>& haplotype_traversal,
                             int called_trav_i, int called_trav_j, size_t ploidy) const;

    mutable std::mutex mutex;
};

}

#endif
