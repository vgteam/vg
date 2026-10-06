#ifndef VG_LINKAGE_MODEL_HPP_INCLUDED
#define VG_LINKAGE_MODEL_HPP_INCLUDED

/** \file linkage_model.hpp
 * The linkage model: a Li-Stephens hidden Markov model over the haplotype panel, which
 * re-decides per-site genotypes using the combinations of alleles that panel haplotypes
 * carry at neighbouring sites. It also phases the calls.
 *
 * The model is described in doc/read-likelihood-linkage-model.md.
 */

#include <algorithm>
#include <atomic>
#include <cstddef>
#include <cstdint>
#include <map>
#include <mutex>
#include <string>
#include <unordered_set>
#include <vector>

namespace vg {

using namespace std;

/// Counters for the linkage pass, reported under --progress.
///
/// Below the top level, the model decodes one *group* at a time: the sites of one child chain at
/// one ploidy under one parent (see `LinkageCollector`). Several counters describe groups.
///
/// They are members of LinkageModel so that each run counts separately. They are atomic because
/// several threads count at once, and the model holds them as `mutable` so that its const
/// methods can count.
struct LinkageCounters {
    /// Pinned sites that `window_phasing` fixed to their pinned pair (`pin_applied`), and those it
    /// left free because that pair cannot spell the genotype the site is constrained to
    /// (`pin_declined`).
    std::atomic<size_t> pin_applied{0}, pin_declined{0};

    /// Groups whose parent had no PhaseCall to pin it to, and groups whose parent had one.
    std::atomic<size_t> group_parent_unpinned{0}, group_parent_pinned{0};

    /// Sites decoded alone rather than in a group: those with no parent key, those whose parent
    /// has no live entry, and the total of the two.
    std::atomic<size_t> grp_no_parent{0}, grp_no_entry{0}, grp_vetoed{0};

    /// Nested chains at ploidy 1 that were given a strand of their parent, and those whose
    /// parent is haploid, so that there is only one strand.
    std::atomic<size_t> nest_strand{0}, nest_one_hap{0};

    /// Nested chains at ploidy 1 under a diploid parent that were given no strand: those that
    /// both of the parent's chosen alleles cross, which happens where the linkage pass could not
    /// revise the chain's ploidy, and those whose parent's chosen pair could not be read.
    std::atomic<size_t> nest_both{0}, nest_unreadable{0};
};

class LinkageModel {
public:

    struct Params {
        /// Exponent on the switch probability, rho^weight (--linkage-weight). Zero turns the
        /// model off, and larger values make switches rarer. Only the switch probability is
        /// raised to the power, and the probability of staying is 1 minus the result, so the
        /// transitions still sum to 1.
        double weight = 2.0;

        /// Distance over which linkage decays, in bp (--linkage-scale).
        double scale = 10000.0;

        /// Floor on the switch probability, so that a switch is never impossible.
        double rho_min = 1e-3;

        /// Escape probability for each strand whose allele is unknown, because it copies the
        /// wildcard haplotype or a panel haplotype that does not pass through the site. The
        /// wildcard can carry any candidate allele, so a genotype that no panel pair spells
        /// can still be called.
        double escape = 1e-2;

        /// Exponent F on the allele-frequency prior that the states imply (--linkage-prior). The
        /// probability collected for a genotype that c ordered panel pairs spell is multiplied by
        /// c^(F-1). 1 keeps the prior as the states imply it, 0 removes it, and larger values
        /// strengthen it.
        double freq_prior = 5.0;

        /// Exponent used instead of `freq_prior` at a site whose alleles differ in the length of
        /// a homopolymer run (see `run_length_site`) (--hp-prior); 0 turns it off.
        ///
        /// Sequencing errors in a long run tend to recur in many reads at the same site, which
        /// the per-site likelihood counts as independent evidence, so its margin grows with
        /// depth. A larger exponent there keeps the panel's share of the decision.
        double hp_prior = 0.0;
        /// Shortest run, measured in the alleles' own sequence, to which `hp_prior` applies
        /// (--hp-prior-run).
        size_t hp_prior_run = 11;

        /// Sites per window of exact inference, and the sites discarded at each end.
        ///
        /// A long linkage chain is decoded as overlapping windows, each on its own, and only the
        /// posteriors away from a window's edges are kept. Linkage decays over a few sites, so a
        /// posterior with a margin of sites on each side is close to the whole chain's.
        size_t window = 2000;
        size_t margin = 250;
    };

    /// One site's input to the model. It holds no graph or GBWT types, so the model can be
    /// tested on numbers alone.
    struct Site {
        /// Where the site starts on the reference: the first base of its first boundary node.
        /// Distances between sites are differences of these.
        size_t position = 0;

        /// Number of alleles, so genotype indices can be decoded.
        size_t num_alleles = 0;

        /// The frequency exponent this site decodes with; negative means `Params::freq_prior`.
        double freq_prior = -1.0;

        /// ln P(reads | genotype).
        ///
        /// At `ploidy` 2 this is in VCF genotype order: index of (i,j) with i <= j is
        /// j*(j+1)/2 + i. At `ploidy` 1 a genotype *is* an allele, so it is indexed by allele
        /// directly and has `num_alleles` entries.
        vector<double> genotype_ln_likelihood;


        /// This site is in a chain that no reference path passes through, so `position` is its
        /// parent's reference start plus the chain's offset along the parent's allele. Two such
        /// sites of one chain are separated by the difference of their positions, and so are such
        /// a site and its positioned parent (`group_parent`). Between such a site and any other
        /// positioned site, no distance is known.
        bool unpositioned = false;

        /// This site is the parent its group is decoded under, held as the group's first site, so
        /// the next site is the parent's child and lies the child's offset along the parent's allele
        /// from it, whether or not the child is positioned.
        bool group_parent = false;

        /// 1 or 2. All sites of a linkage chain have the same ploidy.
        size_t ploidy = 2;

        /// Allele carried by each panel haplotype, or -1 where the haplotype does not pass
        /// through this site and so carries no allele here.
        vector<int> haplotype_allele;

        /// Fix this site's haplotype pair in `phasing()` to (pin_first, pin_second).
        ///
        /// A group is phased with its parent as its first site, and the parent's phase is already
        /// chosen. Pinning the parent keeps the path from swapping its strands, which would phase
        /// the group against the wrong strands. `(size_t)-1` is `WILDCARD`, which is declared
        /// below.
        bool pinned = false;
        size_t pin_first = (size_t)-1;
        size_t pin_second = (size_t)-1;
    };

    LinkageModel(const Params& params) : params(params) {}

    /// True when the model is on, that is, when its weight is positive. The caller checks this
    /// rather than running the model at weight 0.
    bool active() const { return params.weight > 0.0; }

    /// Posterior over genotypes per site, in the same order as `genotype_ln_likelihood`.
    /// `sites` must be one chain in reference order. Returns an empty vector per site where no
    /// posterior could be formed.
    ///
    /// `ploidy` selects the states: single panel haplotypes at 1, ordered pairs at 2. The two
    /// ploidies are decoded by separate functions behind this one entry point, so that an
    /// argument cannot be added to one and forgotten on the other.
    ///
    /// `alpha_in`, when given, replaces the uniform distribution over states at the chain's first
    /// site. The linkage pass builds it from the parent's chosen state. Where the chain's first
    /// site is the parent itself, it is a point mass there, which fixes that site's state; later
    /// sites can switch away from it. A ploidy-1 chain under a diploid parent does not hold the
    /// parent, so its message is the haplotype of the parent's strand that carries the chain,
    /// carried through one transition to the chain's first site.
    vector<vector<double>> posteriors(const vector<Site>& sites, size_t ploidy = 2,
                                      const vector<double>* alpha_in = nullptr) const;

    /// One strand's assignment at one site: an index into the panel, or `WILDCARD`.
    struct Phase {
        size_t first = WILDCARD;
        size_t second = WILDCARD;
    };

    /// The wildcard haplotype's index. It can carry any allele at any site, so a strand
    /// assigned to it is explained by no panel haplotype.
    static constexpr size_t WILDCARD = (size_t)-1;

    /// Whether the reference allele `ref` and another allele differ only in the length of one
    /// homopolymer run, by 1-49 copies of its base, in a run that is at least `min_run` long in
    /// the longer of the two or that reaches either end of that allele. Any pair of alleles counts
    /// when `ref` is out of range, as it is for a chain that no reference path passes through.
    ///
    /// A run that reaches an end of the allele continues into the neighbouring site, where the
    /// graph cut it, so its length is unknown and it counts as long.
    static bool run_length_site(const vector<string>& alleles, size_t min_run,
                                size_t ref = (size_t)-1);

    /// Most probable path of haplotype pairs through the chain, found by max-product (Viterbi)
    /// decoding. This is the phasing. `posteriors()` instead decides each site on its own, and
    /// its per-site answers need not form a path that one pair of haplotypes can spell.
    ///
    /// `constraint[t]` is the genotype index the path must spell at site `t`, or `NO_CONSTRAINT`
    /// to leave the site free. Constraining every site to its chosen genotype makes the
    /// phasing agree with the VCF. A constrained path always exists, because the wildcard can
    /// carry any allele.
    ///
    /// At `ploidy` 1 there is one strand, so the result gives only the panel haplotype it
    /// copies at each site, with `second` the wildcard. As for `posteriors()`, the two
    /// ploidies share this entry point.
    ///
    /// A group's parent fixes the group's starting state in one of two ways. A ploidy-2 group
    /// holds its parent as its first site, pinned to the parent's pair. A ploidy-1 group under a
    /// diploid parent cannot hold the parent, which has another ploidy, so `alpha_in` carries the
    /// parent's state instead, as for `posteriors()`. `alpha_in` is used at ploidy 1 only.
    vector<Phase> phasing(const vector<Site>& sites, const vector<size_t>& constraint,
                          size_t ploidy = 2, const vector<double>* alpha_in = nullptr) const;

    /// Leaves a site's genotype unconstrained in `phasing()`.
    static constexpr size_t NO_CONSTRAINT = (size_t)-1;

    /// VCF diploid genotype ordering: index of the genotype (i,j).
    static size_t genotype_index(size_t i, size_t j) {
        if (i > j) {
            size_t t = i; i = j; j = t;
        }
        return j * (j + 1) / 2 + i;
    }

    /// Per-strand switch probability between two sites `gap` bp apart, after weighting.
    double switch_probability(size_t gap) const;

    /// Posteriors over sites [from, to) decoded as a single window, starting from `alpha_in` and
    /// ending at `beta_in` (uniform where null). The linkage pass decodes a diploid group this
    /// way, from its parent's chosen state.
    void segment_posteriors(const vector<Site>& sites, size_t from, size_t to,
                            const vector<double>* alpha_in, const vector<double>* beta_in,
                            vector<vector<double>>& out) const {
        // `window_posteriors` writes `out[from + t]` unchecked, so size `out` here. It is grown,
        // not reassigned, so that several segments can be decoded into one buffer.
        if (out.size() < sites.size()) {
            out.resize(sites.size());
        }
        window_posteriors(sites, from, to, out, alpha_in, beta_in);
    }

    /// Counters; `mutable` so that const methods can count. See `LinkageCounters`.
    mutable LinkageCounters counters;

private:

    /// Exact forward-backward over one window. `out` is filled for the whole window; the caller
    /// keeps only the interior.
    ///
    /// `alpha_in` and `beta_in` are the messages over haplotype pairs entering the window's two
    /// ends, m*m entries each in the same (a * m + b) layout as the emissions; uniform when null.
    void window_posteriors(const vector<Site>& sites, size_t from, size_t to,
                           vector<vector<double>>& out,
                           const vector<double>* alpha_in = nullptr,
                           const vector<double>* beta_in = nullptr) const;

    /// Max-product over one window, with an optional pinned state so that consecutive windows
    /// join without a spurious switch between them. `out` is indexed from `from`.
    void window_phasing(const vector<Site>& sites, size_t from, size_t to,
                        const vector<size_t>& constraint,
                        size_t pin_index, const Phase& pin, vector<Phase>& out) const;


    /// Emission over single haplotypes for a haploid site: `e[a]` is the relative likelihood of
    /// the allele haplotype `a` carries, with the wildcard last.
    void haploid_emission(const Site& site, size_t n_hap, vector<double>& e,
                          vector<double>& per_allele) const;

    /// Forward-backward and max-product over one window of a ploidy-1 chain, whose states are
    /// single haplotypes.
    void window_haploid_posteriors(const vector<Site>& sites, size_t from, size_t to,
                                   vector<vector<double>>& out,
                                   const vector<double>* alpha_in = nullptr) const;
    void window_haploid_phasing(const vector<Site>& sites, size_t from, size_t to,
                                const vector<size_t>& constraint,
                                size_t pin_index, size_t pin, vector<size_t>& out,
                                const vector<double>* alpha_in = nullptr) const;

    Params params;
};

/// One step of the forward or backward transition over ordered pairs, with `m` states per strand.
///
/// The two strands switch independently, so each has its own switch probability, `rho_a` for
/// strand 0 and `rho_b` for strand 1. Declared here so that the unit tests can call it.
void transition_apply(const std::vector<double>& in, size_t m, double rho_a, double rho_b,
                      std::vector<double>& out);

/**
 * Keeps a compact entry for each genotyped site, and runs the linkage model over the entries.
 *
 * An entry holds only what the model needs: the genotype likelihoods over the site's *compact
 * allele space* (the called pair plus every allele some panel haplotype carries), and the allele
 * each panel haplotype carries. So the whole genome's entries stay small. Sites are recorded from
 * several threads in no fixed order.
 *
 * Sites are resolved one level at a time. Level 0's sites form the top-level linkage
 * chains: one contig, split where its ploidy changes, sorted by position, since the transitions
 * depend on distance. Each later site belongs to a *group*: the sites of one child chain at one
 * ploidy under one parent, decoded from the parent's chosen state. Resolving chooses each site's
 * genotype and phases it, as a `PhaseCall`.
 */
class LinkageCollector {
public:

    LinkageCollector(const LinkageModel::Params& params, size_t num_haplotypes)
        : params(params), model(params), n_haplotypes(num_haplotypes) {}

    /// How many panel haplotypes this collector was built for. The width of every per-haplotype
    /// row it hands out or takes in.
    size_t panel_size() const { return n_haplotypes; }

    /// The direct call's quality inputs, which `VCFOutputCaller::write_variants` needs to rewrite
    /// the quality fields of a record whose genotype the model changed. They come from the call the
    /// site was recorded with.
    struct DirectQuality {
        /// The fraction of the site's reads whose best allele is in the direct call.
        double explained_share = 1.0;
        /// The factor the direct call's GQ multiplied its likelihood gap by
        /// (`ReadLikelihoodSnarlCaller::gq_factor`).
        double gq_factor = 1.0;
        /// The divisor of the direct call's GQN, in nats; 0 where there was none.
        double achievable_gap = 0.0;
    };

    /// A moved site's chosen-genotype posterior, with its direct call's quality inputs.
    struct MovedQuality {
        double posterior = 0.0;
        DirectQuality direct;
    };

    /// Where a site sits in the snarl tree. Callers fill it with designated initialisers, so each
    /// field is named at the call site; as separate bool and integer arguments, a missing one
    /// would shift the rest and still compile.
    struct SiteContext {
        /// The site's chain has one copy: exactly one of the parent's alleles crosses it (the
        /// called alleles when the direct pass records it, the chosen ones when the linkage pass does).
        bool nested = false;
        size_t parent_record_key = 0;
        /// Bit t set iff the parent's candidate traversal t crosses this chain.
        uint64_t parent_crossing = 0;
        /// The site's depth of descent: its parent's level plus 1 for a chain reached by
        /// descent, and 0 for a site the direct pass calls directly, which includes the children of a
        /// snarl it could not genotype. The linkage pass chooses its genotype at its level.
        size_t level = 0;
        /// Whether a VCF line exists for this site. `vg call` records every site before its line
        /// is written, so it passes false and supplies the answer through `set_allele_map`.
        bool emitted = true;
        /// The site has no reference position, and `position` stands in for it, as described at
        /// `Site::unpositioned`.
        bool unpositioned = false;
        /// A hash of the chain's boundary nodes, which identifies the chain and groups its sites.
        size_t chain_key = 0;
        /// This site's frequency exponent (`Params::hp_prior` at a run-length site); negative
        /// means the model's `freq_prior`.
        double freq_prior = -1.0;
    };

    /// Record one genotyped site. Safe to call from several threads.
    ///
    /// Alleles are candidate traversal indices, not VCF allele numbers. `haplotype_traversal` has
    /// one entry per panel haplotype: the candidate traversal it carries, or -1 where it does not
    /// pass through the site. `called_trav_i/j` is the pair the per-site likelihood chose.
    /// `traversal_to_allele` maps candidate traversals to the VCF alleles they were written as;
    /// pass it empty, with `ctx.emitted` false, while no line exists.
    ///
    /// A site that gets no line is still recorded, because its allele pair phases its children.
    void record(const string& contig, size_t position,
                const map<vector<int>, double>& genotype_ln_likelihood,
                const vector<int>& haplotype_traversal,
                int called_trav_i, int called_trav_j,
                const vector<int>& traversal_to_allele,
                size_t record_key,
                const DirectQuality& direct, size_t ploidy,
                int64_t start_node, int64_t end_node,
                const SiteContext& ctx);

    /// The compact allele space `record` builds for one site: the called pair plus every traversal
    /// some panel haplotype carries, without duplicates and sorted. Public for the unit tests.
    /// `genotype_ln_likelihood` is keyed by sorted candidate-traversal vectors.
    static vector<int> compact_allele_space(const map<vector<int>, double>& genotype_ln_likelihood,
                                            const vector<int>& haplotype_traversal,
                                            int called_trav_i, int called_trav_j);

    /// One site's phasing: which strand carries which allele, and which panel haplotype each
    /// strand copies there.
    ///
    /// The pair is ordered, and the order is the phase: the `_first` allele is on the same strand
    /// as every other `_first` allele in the same phase set.
    struct PhaseCall {
        size_t record_key = 0;
        string contig;
        size_t position = 0;
        /// The VCF allele on each strand, or `LinkageModel::WILDCARD` where the site's
        /// traversal-to-allele map was not yet known when the site was phased (see
        /// `set_allele_map`). A record's GT is phased from `trav_first` and `trav_second`, which
        /// are always known.
        size_t allele_first = 0;
        size_t allele_second = 0;
        /// The candidate traversal on each strand. A record's GT is phased from these, and a child
        /// chain's strand is found from them, since the parent's crossing mask is indexed by
        /// candidate traversal.
        int trav_first = -1;
        int trav_second = -1;
        /// The panel haplotype each strand copies here, which the mosaic writes;
        /// `LinkageModel::WILDCARD` where no panel haplotype explains the strand.
        size_t hap_first = LinkageModel::WILDCARD;
        size_t hap_second = LinkageModel::WILDCARD;
        /// 1 or 2. At 1 only the `_first` fields are meaningful: there is one strand.
        size_t ploidy = 2;
        /// The site's boundary nodes. The mosaic locates sites by these, since a reference
        /// position depends on the reference path.
        int64_t start_node = 0;
        int64_t end_node = 0;
        /// The phase set. Phase is comparable only within one.
        size_t phase_set = 0;
        /// For a nested site at ploidy 1, which of the parent's two strands carries it, 0 or 1;
        /// -1 otherwise, including a nested site whose parent could not be found or whose strand
        /// is not determined. The VCF writes it as `a|.` or `.|a`.
        int8_t nested_strand = -1;
        /// True when nothing ordered the pair: the site is heterozygous and no panel haplotype on
        /// either strand carries either called allele, so the pair is in sorted order and the
        /// written phase is arbitrary.
        bool order_arbitrary = false;

        /// The site's level (see `SiteContext::level`). The mosaic writer uses it to
        /// find where a strand enters or leaves a nested chain, and to leave nested sites out
        /// under --no-mosaic-nested.
        uint8_t level = 0;
    };

    /// Resolve level 0 only, per contig, and return how many sites the model moved off their
    /// called genotype. For unit tests; `vg call` resolves each level with
    /// `resolve_level`.
    ///
    /// With `phasing_out`, also returns a phasing of the chosen genotypes, after the model's
    /// changes, so that the phasing agrees with the VCF.
    size_t resolve(vector<PhaseCall>* phasing_out = nullptr) {
        return resolve_level(0, true, phasing_out);
    }

    /// By record key, each live site whose genotype the model changed when its level was
    /// last resolved. `VCFOutputCaller::write_variants` rewrites GQ, GQN and FILTER on each such
    /// record's rendered line from these, since the posterior exists only here, and
    /// `FlowCaller::anchor_gqn_for` uses them for the anchors' gqn column.
    ///
    /// Resolving a site's level adds or removes its key, and retracting the site removes it,
    /// so after the linkage pass runs again the map describes that run alone.
    const std::unordered_map<size_t, MovedQuality>& moved_quality() const {
        return moved_quality_by_record;
    }

    /// Resolve one level of sites, holding every earlier level fixed.
    ///
    /// A group's parent, of the level before, is the one earlier site it holds, and only
    /// when their ploidies match. It is clamped: its emission becomes a point mass at its chosen
    /// genotype and its phase is pinned, so it starts the group from its chosen state and cannot
    /// change. Later levels are left out, since their ploidies are not known until this
    /// level is chosen.
    ///
    /// Each site of this level gets one `PhaseCall`, appended to `phasing_out`, which must be
    /// passed back in on every call of one linkage pass, since a nested site's strand is read from
    /// its parent's `PhaseCall`. `last` marks the pass's final level, whose call sorts
    /// `phasing_out` into reference order.
    ///
    /// Returns how many sites the model moved off their called genotype.
    size_t resolve_level(size_t level, bool last,
                                      vector<PhaseCall>* phasing_out = nullptr);

    /// Mark a site's live entry retracted, so that the model no longer sees it.
    ///
    /// The linkage pass retracts a chain before its level resolves when the parent's chosen
    /// genotype does not cross it, so the sample has no copy of it. It also retracts a chain's
    /// entry just before recording the chain again at a revised ploidy. A chain retracted on an
    /// earlier linkage pass that a later pass finds carried is recorded again. Lookups reach the
    /// first live entry for a key, so after retract and then `record` the new entry is the live
    /// one.
    ///
    /// The entry is marked rather than erased, since other entries hold offsets into the shared
    /// arrays, and everything that walks the entries skips it. Returns false for an unknown key.
    bool retract(size_t record_key);

    /// Replace a live entry's genotype likelihoods, and nothing else, as re-genotyping does.
    ///
    /// The compact allele space can change: it is the panel-carried traversals plus the called
    /// pair, so a new call on an allele no panel haplotype carries adds one allele. The new
    /// likelihoods and alleles are then appended to the arrays and the entry's offsets moved to
    /// them, since each entry's slice has a fixed width.
    ///
    /// Returns false only for a key with no live entry, or a space that cannot be compacted.
    bool rescore(size_t record_key, const map<vector<int>, double>& genotype_ln_likelihood,
                 const vector<int>& haplotype_traversal, int called_trav_i, int called_trav_j);


    /// The pair of candidate traversals the linkage model chose for this site.
    ///
    /// Translated out of the compact space here, since compact indices mean nothing outside the
    /// collector. A site the model did not change returns its called pair. Returns false for an
    /// unknown or retracted key.
    bool chosen_traversals(size_t record_key, int* first, int* second, size_t* ploidy) const;

    /// Fill in the traversal-to-VCF-allele map for a site already recorded, and say whether a line
    /// exists for it.
    ///
    /// Sites are recorded when they are genotyped, before their VCF alleles are chosen, so the
    /// writer supplies the map when it writes the line, after every site is phased. Returns false
    /// for an unknown key.
    bool set_allele_map(size_t record_key, const vector<int>& traversal_to_allele, bool emitted);

    /// The keys of every record that ended up with a VCF line, from `set_allele_map`.
    std::unordered_set<size_t> emitted_records() const;


    /// What a parent's chosen pair implies about one of its child chains: how many copies of the
    /// chain the sample carries, which of the parent's chosen traversals carries it when only one
    /// does, and whether the pair could be read at all. Computed where it is needed, from
    /// `relate_to_parent`, rather than stored.
    struct Relation {
        uint8_t copies = 0;
        int carrying_trav = -1;   ///< the traversal when copies == 1; -2 when both carry it
        bool known = false;       ///< false when the parent's chosen pair could not be read
    };

    /// Compute a Relation from the parent's crossing mask and chosen traversals. The linkage pass uses
    /// the copy count to set a child chain's ploidy, and `resolve_level` uses the carrying
    /// traversal to put a ploidy-1 group on one of its parent's strands
    /// (`PhaseCall::nested_strand`). Both call this, so they cannot disagree. `trav_b` is -1 for a
    /// haploid parent.
    static Relation relate_to_parent(uint64_t crossing, int trav_a, int trav_b);


    /// Whether a live (not retracted) entry exists for this key.
    bool has_entry(size_t record_key) const;

    /// Move the live entry for this key to another position on its contig, as when the place of a
    /// site with no reference position changes with its parent's chosen genotype. Returns false
    /// when the key has no live entry.
    bool set_position(size_t record_key, size_t position);

    /// How many sites belong to one level, for reporting a per-level pass.
    size_t num_sites_at(size_t level) const;

    /// The highest level any recorded site belongs to.
    size_t max_level() const;

    /// Bytes held by the collector, for reporting.
    size_t bytes() const;

    size_t num_sites() const { return entries.size(); }

    /// How many times `record()` filed a site under a key that already had a live entry, as a
    /// snarl recorded twice, or two snarls whose names hash alike, would be. `retract` and every
    /// lookup reach only the first live entry for a key, so the second is decoded but cannot be
    /// replaced. The count is reported with the linkage summary.
    size_t num_duplicate_live_keys() const { return duplicate_live_keys; }

    const LinkageModel::Params& model_params() const { return params; }

    /// Live entries that decode with a site-specific frequency exponent (`SiteContext::freq_prior`),
    /// for reporting.
    size_t num_site_prior_entries() const {
        size_t n = 0;
        for (const Entry& e : entries) {
            n += !e.retracted && e.freq_prior >= 0.0f;
        }
        return n;
    }

private:

    /// No entry. Chain terminator for `next_same_key` and the miss value of `live_index`.
    static constexpr uint32_t NO_ENTRY = (uint32_t)-1;

    struct Entry {
        uint32_t position = 0;
        /// `SiteContext::freq_prior`. Placed in the padding before `chain_key`, so it does not
        /// enlarge the entry.
        float freq_prior = -1.0f;
        /// `SiteContext::chain_key`: identifies the site's chain, so that sites of different
        /// chains, which have no transitions between them, are kept apart.
        size_t chain_key = 0;
        uint32_t contig = 0;
        uint32_t gl_offset = 0;
        uint32_t hap_offset = 0;
        /// Offsets into `trav_arena` and `allele_arena`, which map each compact allele to its
        /// candidate traversal and to the VCF allele it was written as.
        ///
        /// The compact allele space is a set of distinct traversals: the called pair plus every
        /// traversal some panel haplotype carries. The model works in this space rather than in VCF
        /// alleles, because symbolic collapsing can write two traversals as one VCF allele: a parent
        /// whose alleles differ only inside a child chain is homozygous in the VCF but heterozygous
        /// here, and only the second can phase the child. The VCF allele is -1 for a traversal that
        /// was not written as one.
        uint32_t trav_offset = 0;
        uint32_t allele_offset = 0;
        uint16_t num_alleles = 0;
        uint16_t called_i = 0;
        uint16_t called_j = 0;
        uint8_t ploidy = 2;
        /// The level at which the linkage pass chooses this site's genotype. Where nothing
        /// nests, every entry is 0 and one resolve chooses them all.
        uint8_t level = 0;
        /// Set by `retract`. The entry stays in place and is skipped.
        bool retracted = false;
        /// The chosen genotype, written when the site's level resolves, so that later
        /// levels can clamp it.
        uint16_t final_i = 0;
        uint16_t final_j = 0;
        /// Whether this site wrote a VCF line. Only then do its `allele_offset` entries mean
        /// anything.
        bool emitted = true;
        /// This site has no reference position, so `position` is not used for distances. See
        /// `SiteContext::unpositioned`.
        bool unpositioned = false;
        /// `SiteContext::nested`.
        bool nested = false;
        /// For a nested site, the record key of its parent site.
        size_t parent_record_key = 0;
        /// The crossing mask: one bit per parent candidate traversal, set where that traversal
        /// crosses this child chain; 0 when descent could not tell. Placed after an 8-byte member
        /// so that it needs no padding.
        uint64_t parent_crossing = 0;
        /// The `DirectQuality` the site was recorded with, as floats.
        float explained_share = 1.0f;
        float gq_factor = 1.0f;
        float achievable_gap = 0.0f;
        /// The next entry with the same record key, in insertion order, or NO_ENTRY. With
        /// `first_by_key`, it lets a lookup by key follow a short list rather than scan every
        /// entry. It fits in the padding before `record_key`.
        uint32_t next_same_key = NO_ENTRY;
        size_t record_key = 0;
        /// The site's boundary nodes, for the mosaic output.
        int64_t start_node = 0;
        int64_t end_node = 0;
    };

    LinkageModel::Params params;
    LinkageModel model;
    size_t n_haplotypes;

    /// The entries. Their variable-length data lives in shared arrays (the arenas below) rather
    /// than in vectors of their own, which would add more per entry than the entry itself.
    vector<Entry> entries;
    /// Entry indices by level, in the order the entries were appended, so that
    /// `resolve_level(k)` visits only level k's entries.
    vector<vector<uint32_t>> by_level;
    vector<float> gl_arena;
    vector<int8_t> hap_arena;
    /// Per compact allele, the candidate traversal index it stands for.
    vector<uint16_t> trav_arena;
    /// Per compact allele, the VCF allele it was emitted as, or -1 for none.
    vector<int8_t> allele_arena;
    vector<string> contig_names;
    /// Reverse of `contig_names`, so that `record()` does not scan it. A gRef cover can give a run
    /// thousands of contigs.
    unordered_map<string, uint32_t> contig_index;
    /// record key -> first and last entry carrying it, so a lookup is a hash probe and a
    /// walk of that key's chain. `last` is what makes appending O(1) rather than a walk.
    std::unordered_map<size_t, uint32_t> first_by_key;
    std::unordered_map<size_t, uint32_t> last_by_key;
    /// See `num_duplicate_live_keys`.
    size_t duplicate_live_keys = 0;


    /// See `moved_quality`.
    std::unordered_map<size_t, MovedQuality> moved_quality_by_record;

    /// The first non-retracted entry with this key, or NO_ENTRY. Call with `mutex` held.
    uint32_t live_index(size_t record_key) const;

    /// Split a chosen compact pair into the traversal and the VCF allele on each strand.
    void finish_phase_call(PhaseCall& pc, const Entry& e) const;

    /// The compact allele space of a site, with its genotype likelihoods translated into it. Built
    /// by `record` and `rescore` alike, so that a re-scored site is described as a newly recorded
    /// one would be.
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

        /// Candidate traversal -> compact allele, or -1, by binary search in `space`.
        int compact_of(int trav) const {
            auto it = std::lower_bound(space.begin(), space.end(), trav);
            return (it == space.end() || *it != trav) ? -1 : (int)(it - space.begin());
        }
    };

    /// Build the compact space and translate the likelihoods into it. `ok` is false where the site
    /// cannot be described: no candidates, more than the 127 an int8 arena can name, or a called
    /// traversal that is not in the space.
    CompactSite compact_site(const map<vector<int>, double>& genotype_ln_likelihood,
                             const vector<int>& haplotype_traversal,
                             int called_trav_i, int called_trav_j, size_t ploidy) const;

    mutable std::mutex mutex;
};

}

#endif
