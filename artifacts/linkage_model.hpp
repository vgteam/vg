#ifndef VG_LINKAGE_MODEL_HPP_INCLUDED
#define VG_LINKAGE_MODEL_HPP_INCLUDED

/** \file linkage_model.hpp
 * The linkage model: a Li-Stephens hidden Markov model over the haplotype panel, which
 * re-decides per-site genotypes using the combinations of alleles that panel haplotypes
 * carry at neighbouring sites. It also phases the calls.
 *
 * The model is described in artifacts/doc/read-likelihood-linkage-model.md.
 */

#include <algorithm>
#include <array>
#include <atomic>
#include <cstddef>
#include <cstdint>
#include <map>
#include <mutex>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

namespace vg {

using namespace std;

/// Counters for the linkage pass, reported with the progress messages.
///
/// Below the top level, the model decodes one *group* at a time: the sites of one child chain at
/// one ploidy under one parent (see `LinkageCollector`). Several counters describe groups.
///
/// Each LinkageModel holds its own counters, so that each run counts separately. They are atomic
/// because several threads count at once.
struct LinkageCounters {
    /// Pinned sites that `window_phasing` fixed to their pinned pair.
    std::atomic<size_t> pin_applied{0};
    /// Pinned sites that `window_phasing` left free, because their pinned pair cannot spell the
    /// genotype the site is constrained to.
    std::atomic<size_t> pin_declined{0};

    /// Groups whose parent had no PhaseCall to pin it to.
    std::atomic<size_t> group_parent_unpinned{0};
    /// Groups whose parent had a PhaseCall to pin it to.
    std::atomic<size_t> group_parent_pinned{0};

    /// Sites decoded alone because they have no parent key.
    std::atomic<size_t> grp_no_parent{0};
    /// Sites decoded alone because their parent has no live entry.
    std::atomic<size_t> grp_no_entry{0};
    /// All sites decoded alone rather than in a group.
    std::atomic<size_t> grp_vetoed{0};

    /// Nested chains at ploidy 1 that were given a strand of their diploid parent.
    std::atomic<size_t> nest_strand{0};
    /// Nested chains at ploidy 1 whose parent is haploid, so that there is only one strand.
    std::atomic<size_t> nest_one_hap{0};

    /// Nested chains at ploidy 1 under a diploid parent that both of the parent's chosen alleles
    /// cross, so that they get no strand. This happens where the linkage pass could not revise the
    /// chain's ploidy.
    std::atomic<size_t> nest_both{0};
    /// Nested chains at ploidy 1 under a diploid parent whose chosen pair could not be read, so
    /// that they get no strand.
    std::atomic<size_t> nest_unreadable{0};

    /// Phased strands whose panel haplotype carries an allele other than the strand's own, which
    /// the model explains as a mutation (`Params::mutation`).
    std::atomic<size_t> mutated_strands{0};
};

/// A Li-Stephens hidden Markov model over the haplotype panel, which decodes one chain of sites at
/// a time.
///
/// A state names the panel haplotype each strand copies: one haplotype at ploidy 1, an ordered
/// pair at ploidy 2. A strand switches haplotype between sites with a probability that grows with
/// their distance, and emits the reads through the genotype its haplotypes' alleles spell, allowing
/// for mutation. The model gives each site a posterior over genotypes, and the chain a phasing.
class LinkageModel {
public:

    /// The model's tuning parameters.
    struct Params {
        /// Exponent on the switch probability, rho^weight. Zero turns the model off, and larger
        /// values make switches rarer. Only the switch probability is raised to the power, and the
        /// probability of staying is 1 minus the result, so the transitions still sum to 1.
        double weight = 2.0;

        /// Distance over which linkage decays, in bp.
        double scale = 10000.0;

        /// Floor on the switch probability, so that a switch is never impossible.
        double rho_min = 1e-3;

        /// Escape probability for each strand whose allele is unknown, because it copies a panel
        /// haplotype that does not pass through the site. Such a strand carries each candidate
        /// allele with probability escape / (alleles).
        double escape = 1e-2;

        /// Li-Stephens mutation probability.
        ///
        /// A strand that copies a panel haplotype carries that haplotype's allele with
        /// probability 1 - mutation, and each of the site's other alleles with probability
        /// mutation / (alleles - 1). So an allele that no panel haplotype carries is explained by
        /// a mutation on the haplotype the strand copies, and the strand keeps copying it.
        ///
        /// Must be above 0 and at most 0.5. At 0 a genotype that no panel pair spells could not
        /// be called.
        double mutation = 3e-3;

        /// Exponent F on the allele-frequency prior that the states imply.
        ///
        /// The probability collected for a genotype that c ordered panel pairs spell is multiplied
        /// by c^(F-1). 1 keeps the prior as the states imply it, 0 removes it, and larger values
        /// strengthen it.
        double freq_prior = 5.0;

        /// Exponent used instead of `freq_prior` at a site whose alleles differ in the length of
        /// a homopolymer run (see `run_length_site`), or 0 for none.
        ///
        /// Sequencing errors in a long run tend to recur in many reads at the same site, which
        /// the per-site likelihood counts as independent evidence, so its margin grows with
        /// depth. A larger exponent there keeps the panel's share of the decision.
        double hp_prior = 0.0;
        /// Shortest run, measured in the alleles' own sequence, to which `hp_prior` applies.
        size_t hp_prior_run = 11;

        /// Sites per window of exact inference.
        ///
        /// A long linkage chain is decoded as overlapping windows, each on its own, and only the
        /// posteriors away from a window's edges are kept. Linkage decays over a few sites, so a
        /// posterior with a margin of sites on each side is close to the whole chain's.
        size_t window = 2000;
        /// Sites added to each end of a window, whose posteriors are discarded.
        size_t margin = 250;
    };

    /// One site's input to the model.
    ///
    /// It holds no graph or GBWT types, so the model can be tested on numbers alone.
    struct Site {
        /// Where the site starts on the reference, at the first base of its first boundary node.
        /// Distances between sites are differences of these.
        size_t position = 0;

        /// Number of alleles, so genotype indices can be decoded.
        size_t num_alleles = 0;

        /// The frequency exponent this site decodes with; negative means `Params::freq_prior`.
        double freq_prior = -1.0;

        /// ln P(reads | genotype).
        ///
        /// At `ploidy` 2 this is in VCF genotype order, by `genotype_index`. At `ploidy` 1 a
        /// genotype is an allele, so it is indexed by allele and has `num_alleles` entries.
        vector<double> genotype_ln_likelihood;

        /// Whether the site is in a chain that no reference path visits.
        ///
        /// Then `position` is its parent's reference start plus the chain's offset along the
        /// parent's allele. Two such sites of one chain are separated by the difference of their
        /// positions, and so are such a site and its positioned parent (`group_parent`). Between
        /// such a site and any other positioned site, no distance is known.
        bool unpositioned = false;

        /// Whether this site is the parent its group is decoded under, held as the group's first
        /// site.
        ///
        /// The next site is then the parent's child, and lies the child's offset along the
        /// parent's allele from it, whether or not the child is positioned.
        bool group_parent = false;

        /// The ploidy, 1 or 2. All sites of a linkage chain have the same ploidy.
        size_t ploidy = 2;

        /// Allele carried by each panel haplotype, or -1 where the haplotype does not pass
        /// through this site and so carries no allele here.
        vector<int> haplotype_allele;

        /// Whether to fix this site's haplotype pair in `phasing()` to (pin_first, pin_second).
        ///
        /// A group is phased with its parent as its first site, and the parent's phase is already
        /// chosen. Pinning the parent keeps the path from swapping its strands, which would phase
        /// the group against the wrong strands. A pin that names `NO_HAPLOTYPE` is not applied.
        bool pinned = false;
        /// The haplotype strand 0 is pinned to, or `NO_HAPLOTYPE` (declared below).
        size_t pin_first = (size_t)-1;
        /// The haplotype strand 1 is pinned to, or `NO_HAPLOTYPE`.
        size_t pin_second = (size_t)-1;
    };

    explicit LinkageModel(const Params& params) : params(params) {}

    /// Return whether the model is on, which it is when its weight is positive.
    ///
    /// The caller checks this rather than running the model at weight 0.
    bool active() const { return params.weight > 0.0; }

    /// Return the posterior over genotypes at each site of a chain.
    ///
    /// `sites` must be one chain in reference order. Each posterior is in the order of
    /// `genotype_ln_likelihood`, or empty where no posterior could be formed.
    ///
    /// `ploidy` selects the states: single panel haplotypes at 1, ordered pairs at 2. The two
    /// ploidies share this entry point, so that an argument cannot be added to one and forgotten
    /// on the other.
    ///
    /// `alpha_in`, when given, replaces the uniform distribution over states at the chain's first
    /// site. The linkage pass builds it from the parent's chosen state. Where the chain's first
    /// site is the parent itself, it is a point mass there, which fixes that site's state, though
    /// later sites can switch away from it. A ploidy-1 chain under a diploid parent does not hold
    /// the parent. Its message is the haplotype that the parent's strand carrying the chain
    /// copies, carried through one transition to the chain's first site.
    vector<vector<double>> posteriors(const vector<Site>& sites, size_t ploidy = 2,
                                      const vector<double>* alpha_in = nullptr) const;

    /// The panel haplotype each strand copies at one site.
    struct Phase {
        /// The haplotype strand 0 copies, or `NO_HAPLOTYPE`.
        size_t first = NO_HAPLOTYPE;
        /// The haplotype strand 1 copies, or `NO_HAPLOTYPE`.
        size_t second = NO_HAPLOTYPE;
    };

    /// The haplotype of a strand for which no panel haplotype is named, such as the second strand
    /// of a ploidy-1 site, or a nested strand that could not be placed.
    static constexpr size_t NO_HAPLOTYPE = (size_t)-1;

    /// Return whether the reference allele `ref` and another allele differ only in the length of
    /// one homopolymer run.
    ///
    /// They must differ by 1-49 copies of the run's base, in a run that is at least `min_run`
    /// long in the longer of the two or that reaches either end of that allele. Any pair of
    /// alleles counts when `ref` is out of range, as it is for a chain that no reference path
    /// visits.
    ///
    /// A run that reaches an end of the allele continues into the neighbouring site, where the
    /// graph cut it, so its length is unknown and it counts as long.
    static bool run_length_site(const vector<string>& alleles, size_t min_run,
                                size_t ref = (size_t)-1);

    /// Return the most probable path of haplotype pairs through the chain, which is its phasing.
    ///
    /// The path is found by max-product (Viterbi) decoding. `posteriors()` instead decides each
    /// site on its own, and its per-site answers need not form a path that one pair of haplotypes
    /// can spell.
    ///
    /// `constraint[t]` is the genotype index the path must spell at site `t`, or `NO_CONSTRAINT`
    /// to leave the site free. Constraining every site to its chosen genotype makes the
    /// phasing agree with the VCF. A constrained path always exists, because a mutation can
    /// carry any allele.
    ///
    /// At `ploidy` 1 there is one strand, so the result gives only the panel haplotype it
    /// copies at each site, with `second` `NO_HAPLOTYPE`. As for `posteriors()`, the two
    /// ploidies share this entry point.
    ///
    /// A group's parent fixes the group's starting state in one of two ways. A ploidy-2 group
    /// holds its parent as its first site, pinned to the parent's pair. A ploidy-1 group under a
    /// diploid parent cannot hold the parent, which has another ploidy, so `alpha_in` carries the
    /// parent's state instead, as for `posteriors()`. `alpha_in` is used at ploidy 1 only.
    ///
    /// @throws runtime_error If `ploidy` is 2 and the sites name more than 65536 panel haplotypes.
    vector<Phase> phasing(const vector<Site>& sites, const vector<size_t>& constraint,
                          size_t ploidy = 2, const vector<double>* alpha_in = nullptr) const;

    /// The constraint that leaves a site's genotype free in `phasing()`.
    static constexpr size_t NO_CONSTRAINT = (size_t)-1;

    /// Return the index of the diploid genotype (i, j) in VCF genotype order.
    static size_t genotype_index(size_t i, size_t j) {
        if (i > j) {
            size_t t = i; i = j; j = t;
        }
        return j * (j + 1) / 2 + i;
    }

    /// Return the per-strand switch probability between two sites `gap` bp apart, after
    /// weighting.
    double switch_probability(size_t gap) const;

    /// Fill `out` with the posteriors over diploid sites [from, to), decoded as a single window.
    ///
    /// The window starts from the message `alpha_in` and ends at `beta_in`, each uniform where
    /// null. The linkage pass decodes a diploid group this way, from its parent's chosen state.
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

    /// Get the counters of this model's runs, which its const methods also count in.
    LinkageCounters& counters() const { return counter_set; }

private:

    /// Fill `out` with posteriors over one window of diploid sites, by exact forward-backward.
    ///
    /// Site `from + t` goes to `out[from + t - out_base]`, so that a window can be decoded into a
    /// buffer of its own size.
    ///
    /// `alpha_in` and `beta_in` are the messages over haplotype pairs entering the window's two
    /// ends, with m*m entries each, at a * m + b. Each is uniform when null.
    void window_posteriors(const vector<Site>& sites, size_t from, size_t to,
                           vector<vector<double>>& out,
                           const vector<double>* alpha_in = nullptr,
                           const vector<double>* beta_in = nullptr,
                           size_t out_base = 0) const;

    /// Fill `out` with the most probable path of haplotype pairs through one window, indexed from
    /// `from`.
    ///
    /// The state at `pin_index` is pinned to `pin`, so that consecutive windows join without a
    /// spurious switch between them. `pin_index` is `(size_t)-1` for no pin.
    ///
    /// @throws runtime_error If the window's sites name more than 65536 panel haplotypes.
    void window_phasing(const vector<Site>& sites, size_t from, size_t to,
                        const vector<size_t>& constraint,
                        size_t pin_index, const Phase& pin, vector<Phase>& out) const;


    /// Fill `e` with the emission over single haplotypes at a haploid site, and `per_allele` with
    /// its allele likelihoods scaled so that the best is 1.
    ///
    /// `e[a]` is the relative likelihood of the reads given that the strand copies haplotype `a`.
    void haploid_emission(const Site& site, size_t n_hap, vector<double>& e,
                          vector<double>& per_allele) const;

    /// Fill `out` with posteriors over one window of a ploidy-1 chain, by forward-backward.
    ///
    /// The states are single haplotypes. Site `from + t` goes to `out[from + t - out_base]`.
    void window_haploid_posteriors(const vector<Site>& sites, size_t from, size_t to,
                                   vector<vector<double>>& out,
                                   const vector<double>* alpha_in = nullptr,
                                   size_t out_base = 0) const;

    /// Return the most probable path of single haplotypes through one window of a ploidy-1 chain.
    ///
    /// The result is indexed from `from`. The state at `pin_index` is pinned to haplotype `pin`,
    /// so that consecutive windows join without a spurious switch. `pin_index` is `(size_t)-1`
    /// for no pin.
    vector<size_t> window_haploid_phasing(const vector<Site>& sites, size_t from, size_t to,
                                          const vector<size_t>& constraint, size_t pin_index,
                                          size_t pin,
                                          const vector<double>* alpha_in = nullptr) const;

    Params params;

    /// The counts of what the model did, mutable so that the const decoding methods can count.
    mutable LinkageCounters counter_set;
};

/// Apply one step of the forward or backward transition to `in`, over ordered pairs with `m`
/// states per strand.
///
/// The two strands switch independently, so each has its own switch probability, `rho_a` for
/// strand 0 and `rho_b` for strand 1. Declared here so that the unit tests can call it.
void transition_apply(const std::vector<double>& in, size_t m, double rho_a, double rho_b,
                      std::vector<double>& out);

/**
 * A site translated into its compact allele space, with its genotype likelihoods.
 *
 * The compact allele space is the called pair plus every candidate traversal some panel haplotype
 * carries. `LinkageCollector::compact_site` builds one for each site it records or re-scores.
 */
struct CompactSite {
    /// Compact allele -> candidate traversal, sorted by candidate index, so the numbering is a
    /// property of the site rather than of the order the genotypes arrived in.
    vector<int> space;
    /// Genotype likelihoods over the compact space, in `LinkageModel::genotype_index` layout.
    vector<float> gls;
    /// The called pair, compacted. Both are >= 0 whenever `ok`.
    int ci = -1;
    int cj = -1;
    size_t site_ploidy = 0;
    bool ok = false;

    /// Candidate traversal -> compact allele, or -1, by binary search in `space`.
    int compact_of(int trav) const {
        auto it = std::ranges::lower_bound(space, trav);
        return (it == space.end() || *it != trav) ? -1 : (int)(it - space.begin());
    }
};

/**
 * The packed store of the linkage collector's site entries, found by record key or by level.
 *
 * Each site is one fixed-size `Entry`. Its variable-length data, which is its likelihoods, the
 * compact allele of each panel haplotype, and the traversal and VCF allele of each compact allele,
 * lives in shared arrays (the arenas), where the entry holds offsets to its slices. Separate
 * vectors per entry would cost more than the entry itself, and the whole genome's entries are held
 * at once.
 *
 * Entries are only appended. A retracted entry stays in place and is skipped by lookups, since
 * other entries hold offsets into the arenas. A slice is never moved or shrunk, so an entry's
 * offsets stay valid. A record key can have several entries. Its *live* entry is the first one,
 * in insertion order, that is not retracted.
 *
 * The table has no lock. Its owner serialises every write against every other access.
 */
class EntryTable {
public:

    /// No entry. Chain terminator for `Entry::next_same_key` and the miss value of `live_index`.
    static constexpr uint32_t NO_ENTRY = (uint32_t)-1;

    /// The compact record of one genotyped site.
    ///
    /// Fields are ordered to pack the record without padding, since the whole genome's entries are
    /// held at once.
    struct Entry {
        uint32_t position = 0;
        /// The site's frequency exponent from its `SiteContext`, or negative for the model's
        /// `freq_prior`. Placed in the padding before `chain_key`, so it does not enlarge the
        /// entry.
        float freq_prior = -1.0f;
        /// The ID of the start node of the site's chain, from its `SiteContext`. It keeps sites of
        /// different chains, which have no transitions between them, apart.
        size_t chain_key = 0;
        /// The contig, as numbered by `contig_id`.
        uint32_t contig = 0;
        uint32_t gl_offset = 0;
        uint32_t hap_offset = 0;
        /// Offsets into `trav_arena` and `allele_arena`, which map each compact allele to its
        /// candidate traversal and to the VCF allele it was written as.
        ///
        /// The model works in the compact allele space rather than in VCF alleles, because
        /// symbolic collapsing can write two traversals as one VCF allele. A parent whose alleles
        /// differ only inside a child chain is homozygous in the VCF but heterozygous here, and
        /// only the second can phase the child. The VCF allele is -1 for a traversal that was not
        /// written as one.
        uint32_t trav_offset = 0;
        uint32_t allele_offset = 0;
        uint16_t num_alleles = 0;
        /// The called pair, in compact alleles.
        uint16_t called_i = 0;
        uint16_t called_j = 0;
        uint8_t ploidy = 2;
        /// The level at which the linkage pass chooses this site's genotype. Where nothing
        /// nests, every entry is 0 and one resolve chooses them all.
        uint8_t level = 0;
        /// Set by `retract`. The entry stays in place and is skipped.
        bool retracted = false;
        /// The chosen pair, in compact alleles, written when the site's level resolves, so that
        /// later levels can clamp it.
        uint16_t final_i = 0;
        uint16_t final_j = 0;
        /// Whether this site wrote a VCF line. Only then do its `allele_offset` entries mean
        /// anything.
        bool emitted = true;
        /// Whether the site has no reference position, so that `position` is a stand-in. The
        /// difference between a stand-in and another position is a distance only within one
        /// chain, or to the parent of the chain.
        bool unpositioned = false;
        /// Whether the site's chain has one copy, so that exactly one of the parent's alleles
        /// crosses it.
        bool nested = false;
        /// For a nested site, the record key of its parent site.
        size_t parent_record_key = 0;
        /// The crossing mask: one bit per parent candidate traversal, set where that traversal
        /// crosses this child chain; 0 when descent could not tell. Placed after an 8-byte member
        /// so that it needs no padding.
        uint64_t parent_crossing = 0;
        /// The `LinkageCollector::DirectQuality` the site was recorded with, as floats.
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

        /// Return the genotype index of the called pair.
        size_t called_genotype() const;
        /// Return the genotype index of the chosen pair.
        size_t chosen_genotype() const;
    };

    /// Make an empty table whose entries each name an allele for `num_haplotypes` panel
    /// haplotypes.
    explicit EntryTable(size_t num_haplotypes) : n_haplotypes(num_haplotypes) {}

    /// Return how many panel haplotypes each entry names an allele for.
    size_t panel_size() const { return n_haplotypes; }

    /// Return the number of a contig, numbering it if it is new.
    uint32_t contig_id(const string& contig);

    /// Return the name of a contig numbered by `contig_id`.
    const string& contig_name(uint32_t contig) const { return contig_names[contig]; }

    /// Return how many contigs have been numbered. They are numbered from 0.
    size_t num_contigs() const { return contig_names.size(); }

    /// Append an entry for a site, with the site's data put in the arenas.
    ///
    /// Fill every field of `e` but the offsets and `num_alleles`, which this sets from `cs`.
    /// `haplotype_traversal` gives the candidate traversal each panel haplotype carries, or -1.
    /// `traversal_to_allele` maps candidate traversals to VCF alleles, and may be empty while no
    /// line exists.
    void add(Entry e, const CompactSite& cs, const vector<int>& haplotype_traversal,
             const vector<int>& traversal_to_allele);

    /// Replace the likelihoods and panel alleles of an entry with those of `cs`.
    ///
    /// Where the compact space is unchanged, overwrite the likelihoods in place. Otherwise append
    /// new slices and point the entry at them, since a slice of another width does not fit, and
    /// leave the VCF alleles unknown until `set_vcf_alleles`. Leave the called pair alone.
    ///
    /// @return False, changing nothing, where the space is unchanged but `cs` has the wrong number
    /// of likelihoods for the entry's ploidy.
    bool replace_site_data(Entry& e, const CompactSite& cs, const vector<int>& haplotype_traversal);

    /// Rewrite the VCF allele of each compact allele of an entry from `traversal_to_allele`.
    void set_vcf_alleles(const Entry& e, const vector<int>& traversal_to_allele);

    /// Return the index of the live entry for a key, or NO_ENTRY.
    uint32_t live_index(size_t record_key) const;

    /// Mark the live entry for a key retracted.
    ///
    /// The key's next entry, if there is one, becomes live. Return false for a key with no live
    /// entry.
    bool retract(size_t record_key);

    /// Return how many entries were added, including retracted ones.
    size_t size() const { return entries.size(); }
    bool empty() const { return entries.empty(); }

    Entry& operator[](size_t i) { return entries[i]; }
    const Entry& operator[](size_t i) const { return entries[i]; }

    /// Iterate over every entry, retracted ones included, in the order they were added.
    vector<Entry>::const_iterator begin() const { return entries.begin(); }
    vector<Entry>::const_iterator end() const { return entries.end(); }

    /// Return the indices of the entries of one level, retracted ones included, in the order they
    /// were added.
    const vector<uint32_t>& at_level(size_t level) const;

    /// Return the candidate traversal that compact allele `compact` of an entry stands for, or -1.
    int traversal_at(const Entry& e, size_t compact) const;

    /// Return the VCF allele that compact allele `compact` of an entry was written as, or -1.
    int vcf_allele_at(const Entry& e, size_t compact) const;

    /// Return the log-likelihood of genotype `g` of an entry, in `LinkageModel::genotype_index`
    /// layout.
    float likelihood_at(const Entry& e, size_t g) const { return gl_arena[e.gl_offset + g]; }

    /// Return the compact allele panel haplotype `h` carries at an entry, or -1 for none.
    int8_t haplotype_allele(const Entry& e, size_t h) const { return hap_arena[e.hap_offset + h]; }

    /// Return how many times `add` appended an entry under a key that already had a live entry.
    ///
    /// The second entry is not live, so it cannot be looked up or retracted, but it is still
    /// walked.
    size_t num_duplicate_live_keys() const { return duplicate_live_keys; }

    /// Return the bytes held by the entries, the arenas and the key index, for reporting.
    size_t bytes() const;

private:

    /// Append a site's likelihoods, panel alleles, traversals and VCF alleles to the arenas, and
    /// point the entry's offsets and allele count at them.
    void append_site_data(Entry& e, const CompactSite& cs, const vector<int>& haplotype_traversal,
                          const vector<int>& traversal_to_allele);

    size_t n_haplotypes;

    vector<Entry> entries;
    /// Entry indices by level, in the order the entries were appended, so that one level's
    /// entries can be visited without a scan.
    vector<vector<uint32_t>> by_level;
    /// Per entry, its genotype log-likelihoods over its compact allele space, in
    /// `LinkageModel::genotype_index` layout at its ploidy.
    vector<float> gl_arena;
    /// Per entry, the compact allele each panel haplotype carries, or -1 for none.
    vector<int8_t> hap_arena;
    /// Per compact allele, the candidate traversal index it stands for.
    vector<uint16_t> trav_arena;
    /// Per compact allele, the VCF allele it was emitted as, or -1 for none.
    vector<int8_t> allele_arena;
    /// Contig names, by contig number.
    vector<string> contig_names;
    /// Reverse of `contig_names`, so that numbering a contig does not scan it. A gRef cover can
    /// give a run thousands of contigs.
    unordered_map<string, uint32_t> contig_index;
    /// Record key -> the first entry carrying it, so a lookup is a hash probe and a walk of that
    /// key's chain.
    std::unordered_map<size_t, uint32_t> first_by_key;
    /// Record key -> the last entry carrying it, which makes appending O(1) rather than a walk.
    std::unordered_map<size_t, uint32_t> last_by_key;
    /// How many times `add` appended an entry under a key that already had a live entry.
    size_t duplicate_live_keys = 0;
};

/**
 * A site's place on the reference, as a contig and an offset along it.
 *
 * The linkage model orders the sites of a contig by position, measures the gaps between them as
 * differences of positions, and names each phase set (FORMAT/PS) by the position of its first
 * site.
 */
struct SiteLocus {
    /// The contig as the VCF names it.
    ///
    /// This is the locus of a PanSN path name, such as `chr20` for `CHM13#0#chr20`, or else the
    /// whole path name.
    string contig;
    /// A 0-based offset along the contig.
    ///
    /// For a site the reference path visits, it is where the first boundary node of the site
    /// starts on that path. It is not the POS of a record, which trimming the alleles can move.
    /// For a site that no reference path visits, it is a stand-in. The stand-in is the reference
    /// start of the parent plus the offset of the site's chain along the parent's allele. The
    /// difference between a stand-in and another position is a distance only within one chain, or
    /// to the parent of the chain.
    size_t position = 0;
};

/**
 * A ledger of compact entries for the genotyped sites, which runs the linkage model over them.
 *
 * An entry holds only what the model needs, so that the whole genome's entries stay small. That
 * is the genotype likelihoods over the site's *compact allele space*, and the allele each panel
 * haplotype carries. The compact allele space is the called pair plus every candidate traversal
 * some panel haplotype carries. Sites are recorded from several threads in no fixed order.
 *
 * Sites are resolved one level at a time. Level 0's sites form the top-level linkage chains, each
 * a run of one ploidy on one contig, sorted by position, since the transitions depend on distance.
 * Each later site belongs to a *group*: the sites of one child chain at one ploidy under one
 * parent, decoded from the parent's chosen state. Resolving chooses each site's genotype and
 * phases it, as a `PhaseCall`.
 */
class LinkageCollector {
public:

    /// Make an empty collector for a panel of `num_haplotypes` haplotypes, whose model decodes
    /// with `params`.
    LinkageCollector(const LinkageModel::Params& params, size_t num_haplotypes)
        : params(params), model(params), entries(num_haplotypes) {}

    /// Return how many panel haplotypes this collector was built for.
    ///
    /// This is the width of every per-haplotype row it hands out or takes in.
    size_t panel_size() const { return entries.panel_size(); }

    /// The quality inputs of the call a site was recorded with.
    ///
    /// `SortedVCFWriter::write_variants` needs them to rewrite the quality fields of a record
    /// whose genotype the model changed.
    struct DirectQuality {
        /// The fraction of the site's reads whose best allele is in the direct call.
        double explained_share = 1.0;
        /// The factor the direct call's GQ multiplied its likelihood gap by
        /// (`ReadLikelihoodSnarlCaller::gq_factor`).
        double gq_factor = 1.0;
        /// The divisor of the direct call's GQN, in nats; 0 where there was none.
        double achievable_gap = 0.0;
    };

    /// The per-site genotyper's verdict: the called pair, its ploidy, and the quality inputs it
    /// was written with.
    ///
    /// Alleles are candidate traversal indices.
    struct DirectCall {
        /// The first called allele.
        int called_i = -1;
        /// The second called allele, or `called_i` again at ploidy 1.
        int called_j = -1;
        /// The number of alleles in the call, 1 or 2.
        size_t ploidy = 2;
        /// The quality inputs of the call.
        DirectQuality quality;
    };

    /// The quality of a site whose genotype the model changed.
    struct MovedQuality {
        /// The posterior of the chosen genotype.
        double posterior = 0.0;
        /// The quality inputs of the call the site was recorded with.
        DirectQuality direct;
    };

    /// Where a site sits in the graph and in the snarl tree.
    ///
    /// Callers fill it with designated initialisers, so each field is named at the call site.
    struct SiteContext {
        /// The ID of the start boundary node of the site (see `PhaseCall::start_node`).
        int64_t start_node = 0;
        /// The ID of the end boundary node of the site.
        int64_t end_node = 0;
        /// Whether the site's chain has one copy, so that exactly one of the parent's alleles
        /// crosses it. Those are the called alleles when the direct pass records the site, and
        /// the chosen ones when the linkage pass does.
        bool nested = false;
        /// The record key of the site's parent, or 0 for none.
        size_t parent_record_key = 0;
        /// The crossing mask, with bit t set iff the parent's candidate traversal t crosses this
        /// chain.
        uint64_t parent_crossing = 0;
        /// The level at which the linkage pass chooses the site's genotype.
        ///
        /// It is the parent's level plus 1 for a chain reached by descent, and 0 for a site the
        /// direct pass calls directly, which includes the children of a snarl it could not
        /// genotype.
        size_t level = 0;
        /// Whether a VCF line exists for this site.
        ///
        /// `vg call` records every site before its line is written, so it passes false and
        /// supplies the answer through `set_allele_map`.
        bool emitted = true;
        /// Whether the site has no reference position, so that the `SiteLocus::position` it is
        /// recorded with is a stand-in. The stand-in is the reference start of the parent plus
        /// the offset of the site's chain along the parent's allele.
        bool unpositioned = false;
        /// The ID of the start node of the chain the site is in, which identifies the chain and
        /// groups its sites.
        size_t chain_key = 0;
        /// The site's frequency exponent, such as `Params::hp_prior` at a run-length site, or
        /// negative for the model's `freq_prior`.
        double freq_prior = -1.0;
    };

    /**
     * Record one genotyped site. Safe to call from several threads.
     *
     * A site that gets no line is still recorded, because its allele pair phases its children.
     * A site with no likelihoods or no first called allele is not recorded, and neither is one
     * whose compact allele space cannot be described.
     *
     * @param genotype_ln_likelihood The likelihood of each genotype of candidate traversals.
     * @param haplotype_traversal One entry per panel haplotype, holding the candidate traversal
     * it carries, or -1 where it does not visit the site.
     * @param traversal_to_allele Maps candidate traversals to the VCF alleles they were written
     * as. Pass it empty, with `ctx.emitted` false, while no line exists.
     */
    void record(const SiteLocus& locus,
                const map<vector<int>, double>& genotype_ln_likelihood,
                const vector<int>& haplotype_traversal, const DirectCall& call,
                const vector<int>& traversal_to_allele, size_t record_key,
                const SiteContext& ctx);

    /// Return the compact allele space `record` builds for one site.
    ///
    /// It holds the called pair plus every traversal some panel haplotype carries, without
    /// duplicates and sorted. Public for the unit tests.
    static vector<int> compact_allele_space(const vector<int>& haplotype_traversal,
                                            int called_trav_i, int called_trav_j);

    /// One site's phasing, which says which strand carries which allele, and which panel
    /// haplotype each strand copies there.
    ///
    /// The pair is ordered, and the order is the phase. The `_first` allele is on the same strand
    /// as every other `_first` allele in the same phase set.
    struct PhaseCall {
        size_t record_key = 0;
        string contig;
        size_t position = 0;
        /// The VCF allele on strand 0, or the largest `size_t` where the site's
        /// traversal-to-allele map was not yet known when the site was phased (see
        /// `set_allele_map`).
        size_t allele_first = 0;
        /// The VCF allele on strand 1, or the largest `size_t` where the site's
        /// traversal-to-allele map was not yet known when the site was phased.
        size_t allele_second = 0;
        /// The candidate traversal on strand 0, which is always known.
        ///
        /// A record's GT is phased from the traversals. A child chain's strand is found from them
        /// too, since the parent's crossing mask is indexed by candidate traversal.
        int trav_first = -1;
        /// The candidate traversal on strand 1.
        int trav_second = -1;
        /// The panel haplotype strand 0 copies here, or `LinkageModel::NO_HAPLOTYPE` where none
        /// is named. The mosaic writes it, whether or not it carries the strand's allele.
        size_t hap_first = LinkageModel::NO_HAPLOTYPE;
        /// The panel haplotype strand 1 copies here, or `LinkageModel::NO_HAPLOTYPE` where none
        /// is named. The mosaic writes it, whether or not it carries the strand's allele.
        size_t hap_second = LinkageModel::NO_HAPLOTYPE;
        /// The ploidy, 1 or 2. At 1 there is one strand, and only the `_first` fields mean
        /// anything.
        size_t ploidy = 2;
        /// The site's start boundary node. The mosaic locates sites by their boundary nodes,
        /// since a reference position depends on the reference path.
        int64_t start_node = 0;
        /// The site's end boundary node.
        int64_t end_node = 0;
        /// The phase set. Phase is comparable only within one.
        size_t phase_set = 0;
        /// For a nested site at ploidy 1, which of the parent's two strands carries it, 0 or 1.
        ///
        /// It is -1 otherwise, including for a nested site whose parent could not be found or
        /// whose strand is not determined. The VCF writes it as `a|.` or `.|a`.
        int8_t nested_strand = -1;
        /// Whether nothing ordered the pair, so that the written phase is arbitrary.
        ///
        /// This happens at a heterozygous site where no panel haplotype on either strand carries
        /// either called allele, and the pair is then in sorted order.
        bool order_arbitrary = false;

        /// The site's level (see `SiteContext::level`). The mosaic writer uses it to find where a
        /// strand enters or leaves a nested chain, and to leave nested sites out when
        /// `MosaicParams::keep_nested` is off.
        uint8_t level = 0;
    };

    /// Resolve level 0 only, and return how many sites the model moved off their called genotype.
    ///
    /// This is for unit tests. `vg call` resolves each level with `resolve_level`. With
    /// `phasing_out`, also return a phasing of the chosen genotypes, so that the phasing agrees
    /// with the VCF.
    size_t resolve(vector<PhaseCall>* phasing_out = nullptr) {
        return resolve_level(0, true, phasing_out);
    }

    /// Get the quality of each live site whose genotype the model changed when its level was last
    /// resolved, by record key.
    ///
    /// `SortedVCFWriter::write_variants` rewrites GQ, GQN and FILTER on each such record's
    /// rendered line from these, since the posterior exists only here. `RecordRenderer::anchor_gqn`
    /// uses them for the anchors' gqn column.
    ///
    /// Resolving a site's level adds or removes its key, and retracting the site removes it. So
    /// after the linkage pass runs again, the map describes that run alone.
    const std::unordered_map<size_t, MovedQuality>& moved_quality() const {
        return moved_quality_by_record;
    }

    /// Resolve one level of sites, holding every earlier level fixed.
    ///
    /// Resolve one level of sites, holding every earlier level fixed.
    ///
    /// A group holds its parent, of the level before, as its one earlier site, but only when
    /// their ploidies match. The parent is clamped. Its emission becomes a point mass at its
    /// chosen genotype and its phase is pinned, so it starts the group from its chosen state and
    /// cannot change. Later levels are left out, since their ploidies are not known until this
    /// level is chosen.
    ///
    /// Each site of this level gets one `PhaseCall`, appended to `phasing_out`. Pass the same
    /// `phasing_out` to every call of one linkage pass, since a nested site's strand is read from
    /// its parent's `PhaseCall`. `last` marks the pass's final level, whose call sorts
    /// `phasing_out` into reference order.
    ///
    /// @return How many sites the model moved off their called genotype.
    size_t resolve_level(size_t level, bool last, vector<PhaseCall>* phasing_out = nullptr);

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

    /// Replace a live entry's genotype likelihoods and called pair, as re-genotyping does.
    ///
    /// The compact allele space can change, since a new call on an allele no panel haplotype
    /// carries adds one allele. The new likelihoods and alleles are then appended to the arrays
    /// and the entry's offsets moved to them, since each entry's slice has a fixed width.
    ///
    /// Returns false only for a key with no live entry, or a space that cannot be compacted.
    bool rescore(size_t record_key, const map<vector<int>, double>& genotype_ln_likelihood,
                 const vector<int>& haplotype_traversal, int called_trav_i, int called_trav_j);


    /// Get the pair of candidate traversals the linkage model chose for this site, and its ploidy.
    ///
    /// A site the model did not change gives its called pair. Returns false for an unknown or
    /// retracted key.
    bool chosen_traversals(size_t record_key, int* first, int* second, size_t* ploidy) const;

    /// Run `chosen_traversals` for every key in `record_keys` at once, on several threads.
    ///
    /// `out[i]` gets the i-th key's {first, second, ploidy}, and `found[i]` is 1 where
    /// `chosen_traversals` would return true.
    void chosen_traversals_for(const vector<size_t>& record_keys, vector<std::array<int, 3>>& out,
                               vector<char>& found) const;

    /// Fill in the traversal-to-VCF-allele map for a site already recorded, and say whether a line
    /// exists for it.
    ///
    /// Sites are recorded when they are genotyped, before their VCF alleles are chosen, so the
    /// writer supplies the map when it writes the line, after every site is phased. Returns false
    /// for an unknown key.
    bool set_allele_map(size_t record_key, const vector<int>& traversal_to_allele, bool emitted);

    /// Return the keys of every record that ended up with a VCF line, from `set_allele_map`.
    std::unordered_set<size_t> emitted_records() const;


    /// What a parent's chosen pair implies about one of its child chains.
    ///
    /// It is computed where it is needed, by `relate_to_parent`, rather than stored.
    struct Relation {
        /// How many copies of the chain the sample carries.
        uint8_t copies = 0;
        /// The parent's chosen traversal that carries the chain when `copies` is 1, -2 when both
        /// carry it, and -1 otherwise.
        int carrying_trav = -1;
        /// Whether the parent's chosen pair could be read.
        bool known = false;
    };

    /// Compute a Relation from the parent's crossing mask and chosen traversals.
    ///
    /// The linkage pass uses the copy count to set a child chain's ploidy, and `resolve_level`
    /// uses the carrying traversal to put a ploidy-1 group on one of its parent's strands
    /// (`PhaseCall::nested_strand`). Both call this, so they cannot disagree. `trav_b` is -1 for a
    /// haploid parent.
    static Relation relate_to_parent(uint64_t crossing, int trav_a, int trav_b);


    /// Return whether a live entry exists for this key.
    bool has_entry(size_t record_key) const;

    /// Move the live entry for this key to another position on its contig.
    ///
    /// This is for a site with no reference position, whose place changes with its parent's
    /// chosen genotype. Returns false when the key has no live entry.
    bool set_position(size_t record_key, size_t position);

    /// Return how many live sites belong to one level, for reporting a per-level pass.
    size_t num_sites_at(size_t level) const;

    /// Return the highest level any recorded site belongs to.
    size_t max_level() const;

    /// Return the bytes held by the collector's entries, for reporting.
    size_t bytes() const { return entries.bytes(); }

    /// Return how many entries were recorded, including retracted ones.
    size_t num_sites() const { return entries.size(); }

    /// Return how many times `record()` added an entry under a key that already had a live entry.
    ///
    /// A snarl recorded twice, or two snarls whose names hash alike, would do that. `retract` and
    /// every lookup reach only the first live entry for a key, so the second is decoded but cannot
    /// be replaced.
    size_t num_duplicate_live_keys() const { return entries.num_duplicate_live_keys(); }

    /// Get the parameters of the model the collector runs.
    const LinkageModel::Params& model_params() const { return params; }

    /// Return how many phased strands so far the model explained by a mutation.
    size_t mutated_strands() const { return model.counters().mutated_strands.load(); }

    /// Return how many live entries decode with a site-specific frequency exponent
    /// (`SiteContext::freq_prior`), for reporting.
    size_t num_site_prior_entries() const {
        size_t n = 0;
        for (const Entry& e : entries) {
            if (!e.retracted && e.freq_prior >= 0.0f) {
                ++n;
            }
        }
        return n;
    }

private:

    /// The record of one site, as the table stores it.
    using Entry = EntryTable::Entry;

    LinkageModel::Params params;
    /// Decodes the linkage chains and groups.
    LinkageModel model;

    /// Holds every recorded site.
    ///
    /// `record`, `rescore`, `retract`, `set_allele_map`, `set_position`, `has_entry`,
    /// `emitted_records` and the `chosen_traversals` lookups hold `mutex` while they use it, so
    /// they can run on several threads at once. Every other method uses it without the lock, so it
    /// must not run while any other method does.
    EntryTable entries;

    /// The quality of each live site whose genotype the model changed when its level was last
    /// resolved, by record key.
    std::unordered_map<size_t, MovedQuality> moved_quality_by_record;

    /// What the parent's chosen pair implies about a child's chain, by `relate_to_parent`.
    Relation relation_to(const Entry& child, const Entry& parent) const;

    /// Set the traversal and the VCF allele on each strand of `pc`, from the compact alleles
    /// `c_first` and `c_second` chosen for its strands.
    ///
    /// Crossing masks are tested against the traversals, and the record writes the VCF alleles.
    void finish_phase_call(PhaseCall& pc, const Entry& e, size_t c_first, size_t c_second) const;

    /// Build the compact space and translate the likelihoods into it. `ok` is false where the site
    /// cannot be described: no candidates, more than the 127 an int8 arena can name, or a called
    /// traversal that is not in the space.
    CompactSite compact_site(const map<vector<int>, double>& genotype_ln_likelihood,
                             const vector<int>& haplotype_traversal,
                             int called_trav_i, int called_trav_j, size_t ploidy) const;

    //------------------------------------------------------------------------------
    // Resolving a level
    //
    // `resolve_level` cuts a level's sites into decode chains, decodes each chain with the model,
    // and applies the results to the entries and the phasing in chain order.

    /// The phase chosen for a site of an earlier level, by which a later level pins the site when
    /// it holds it, and places the site's ploidy-1 children on its strands.
    struct PinnedPhase {
        /// The panel haplotype each strand copies.
        size_t first;
        size_t second;
        /// The candidate traversal on each strand.
        int trav_first;
        int trav_second;
        size_t ploidy;
        /// For a nested site at ploidy 1, which of the parent's two strands carries it, 0 or 1, or
        /// -1 otherwise.
        int nested_strand;
        /// Whether nothing ordered the pair, so that the written phase is arbitrary.
        bool order_arbitrary;
        size_t phase_set;

        /// Return which of this parent's strands, 0 or 1, carries a nested chain at ploidy 1
        /// that the parent traversal `carrying` crosses, or -1 where no strand can be named.
        int strand_carrying(int carrying) const;
    };

    /// Where a nested ploidy-1 site sits on its parent's strands.
    struct NestedPlacement {
        /// The parent strand that carries the site, 0 or 1, or -1 where none is named.
        int strand = -1;
        /// Whether the parent's pair was in arbitrary order. The strand is then that arbitrary
        /// choice, one level down.
        bool order_arbitrary = false;
        /// Whether the site may name a panel haplotype. Under a haploid parent there is one
        /// strand and the child sits on it. Under a diploid parent, the site names none when no
        /// strand carries it alone.
        bool nameable = true;
    };

    /// A run of sites that the model decodes together: a top-level linkage chain, a nested chain
    /// grouped with its parent, or a site that could not be grouped.
    struct DecodeChain {
        /// Entry indices, in decoding order. A group holding its parent holds it first.
        vector<size_t> sites;
        /// The message over the model's states entering the first site, from the parent's
        /// chosen state, or empty for none.
        vector<double> entering;
        /// The phase set to give the chain's sites, or the largest `size_t` to name it by the
        /// chain's first position.
        size_t phase_set = (size_t)-1;
    };

    /// What decoding one chain gives, before it is applied to the entries.
    struct ChainDecode {
        /// The value of `best` for a site with no posterior.
        static constexpr size_t NO_POSTERIOR = (size_t)-1;

        /// The chain's sites of the level being resolved.
        size_t live_here = 0;
        /// The chain's sites of earlier levels that have a phase to pin them to.
        size_t pinned_here = 0;
        /// The genotype each site ends up with, whether or not linkage moved it. The phasing is
        /// constrained to these.
        vector<size_t> final_genotype;
        /// For each site of this level, the most probable genotype under its posterior, or
        /// `NO_POSTERIOR`.
        vector<size_t> best;
        /// The posterior probability of `best`.
        vector<double> best_posterior;
        /// The haplotypes the path chose, when phasing.
        vector<LinkageModel::Phase> phase;
        /// The allele the path's haplotype carries on each strand, or -1 where it names none.
        vector<int> allele_first;
        vector<int> allele_second;
        size_t phase_set = 0;
    };

    /// Whether entry `a` comes before entry `b` in position, then in record key.
    ///
    /// Sites arrive in whatever order the threads finished, and two records can share a position,
    /// so ordering on position alone would make the output depend on the thread count. The key is
    /// derived from the snarl ID, so it is a property of the site rather than of the run.
    bool before_in_position(size_t a, size_t b) const;

    /// The phase of each parent in `parent_keys`, from `phasing`. A key's last PhaseCall wins.
    static unordered_map<size_t, PinnedPhase> pinned_phases(const unordered_set<size_t>& parent_keys,
                                                            const vector<PhaseCall>& phasing);

    /// Cut the live sites of level 0 into linkage chains: maximal runs of one ploidy on one
    /// contig, in reference order.
    vector<DecodeChain> top_level_chains(const vector<uint32_t>& this_level);

    /// Group the live sites of a nested level by parent, chain, ploidy and carrying traversal,
    /// with each group conditioned on its parent's phase. A site that cannot be grouped becomes
    /// a chain of its own.
    ///
    /// Record each grouped ploidy-1 site in `placements`, and add the grouped sites and groups to
    /// `grouped_sites` and `grouped_groups`.
    vector<DecodeChain> parent_groups(const vector<uint32_t>& this_level,
                                      const unordered_set<size_t>& parent_keys,
                                      const unordered_map<size_t, PinnedPhase>& pinned_phase,
                                      unordered_map<size_t, NestedPlacement>& placements,
                                      size_t& grouped_sites, size_t& grouped_groups);

    /// Build the decode chain for one group of `members` under the entry `parent`, at
    /// `group_ploidy`. Sort `members` into position order.
    DecodeChain parent_group(size_t parent, size_t group_ploidy, vector<size_t>& members,
                             const unordered_map<size_t, PinnedPhase>& pinned_phase,
                             unordered_map<size_t, NestedPlacement>& placements);

    /// Return the message over single haplotypes entering a ploidy-1 group under the entry
    /// `parent`, whose phase is `pin`, or empty for none. Record the members' placements.
    vector<double> haploid_group_message(size_t parent, bool parent_in_group,
                                         const vector<size_t>& members, const PinnedPhase& pin,
                                         unordered_map<size_t, NestedPlacement>& placements);

    /// Build the model's input for the entries `indices`. Clamp sites of levels before `level` to
    /// their chosen genotypes, and pin those with a phase in `pinned_phase`.
    vector<LinkageModel::Site> model_sites(const vector<size_t>& indices, size_t level,
                                           const unordered_map<size_t, PinnedPhase>& pinned_phase)
        const;

    /// Fill the genotype fields of `d` from the posteriors of the sites `indices`.
    void choose_genotypes(const vector<size_t>& indices, size_t level,
                          const vector<vector<double>>& posteriors, ChainDecode& d) const;

    /// Decode one chain at `level`, and also phase it if `phase` is set. Writes nothing, so that
    /// chains can be decoded in parallel.
    ChainDecode decode_chain(const DecodeChain& chain, size_t level,
                             const unordered_map<size_t, PinnedPhase>& pinned_phase,
                             bool phase) const;

    /// Decode every chain at `level`, in parallel, and return the results in chain order.
    vector<ChainDecode> decode_chains(const vector<DecodeChain>& chains, size_t level,
                                      const unordered_map<size_t, PinnedPhase>& pinned_phase,
                                      bool phase) const;

    /// Write the genotypes a decode chose into the entries `indices` of `level`, and update
    /// `moved_quality_by_record`. Return how many sites moved off their called genotype.
    size_t apply_chosen_genotypes(const vector<size_t>& indices, const ChainDecode& d,
                                  size_t level);

    /// Build the PhaseCall of entry `e`, the site at index `t` of a decoded chain.
    PhaseCall phase_call_for(const Entry& e, const ChainDecode& d, size_t t,
                             const unordered_map<size_t, NestedPlacement>& placements) const;

    /// Append the PhaseCall of each site of `level` among the entries `indices` of a decoded
    /// chain to `phasing_out`.
    void append_phase_calls(const vector<size_t>& indices, const ChainDecode& d, size_t level,
                            const unordered_map<size_t, NestedPlacement>& placements,
                            vector<PhaseCall>& phasing_out) const;

    /// Return the record keys of the parents of the live sites among `this_level`.
    unordered_set<size_t> parent_keys_of(const vector<uint32_t>& this_level) const;

    /// Log what resolving `level` counted.
    void report_level(size_t level, size_t grouped_sites, size_t grouped_groups) const;

    /// Serialises the methods that may run on several threads at once. See `entries`.
    mutable std::mutex mutex;
};

}

#endif
