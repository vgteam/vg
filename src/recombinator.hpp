#ifndef VG_RECOMBINATOR_HPP_INCLUDED
#define VG_RECOMBINATOR_HPP_INCLUDED

/** \file recombinator.hpp
 * Tools for generating synthetic haplotypes as recombinations of existing
 * haplotypes.
 */

#include "gbwt_helper.hpp"
#include "gbwtgraph_helper.hpp"
#include "hash_map.hpp"
#include "recombinator_haplotypes.hpp"
#include "../artifacts/haplotype_minimizer_index.hpp"
#include "../artifacts/subchain_walk.hpp"
#include "../artifacts/recombinator_model.hpp"
#include "snarl_distance_index.hpp"

#include <functional>
#include <iostream>
#include <optional>

#include <gbwt/cached_gbwt.h>
#include <gbwtgraph/algorithms.h>

namespace vg {

//------------------------------------------------------------------------------

/**
 * A tool for transforming the haplotypes in a GBWT index into a `Haplotypes`
 * representation. Requires a GBZ graph, an r-index, a distance index, and a
 * minimizer index.
 */
class HaplotypePartitioner {
public:
    /// Approximate number of construction jobs to be created.
    constexpr static size_t APPROXIMATE_JOBS = 32;

    /// The amount of progress information that should be printed to stderr.
    typedef Haplotypes::Verbosity Verbosity;

    /// A GBWT sequence as (sequence identifier, offset in a node).
    typedef Haplotypes::sequence_type sequence_type;

    /// An encoded kmer.
    typedef Haplotypes::Subchain::kmer_type kmer_type;

    /// Minimizer index without payloads.
    using minimizer_index_type = HaplotypeMinimizerIndex;

    /**
     * A subchain is a substring of a top-level chain defined by at most two
     * boundary nodes.
     *
     * Normal subchains have two boundary nodes, which are assumed to be the
     * start node of a snarl and the end node of a possibly different snarl.
     * There are assumed to be haplotypes crossing the subchain. Prefixes and
     * suffixes lack one of the boundary nodes, while full haplotypes lack
     * both.
     *
     * When a top-level chain is partitioned into subchains, the boundary nodes
     * may either overlap or be connected by unary paths. If a snarl is not
     * connected, it may be presented as a suffix and a prefix.
     */
    struct Subchain {
        /// The type of this subchain.
        Haplotypes::Subchain::subchain_t type;

        /// Start node.
        handle_t start;

        /// End node.
        handle_t end;

        /// Shortest distance from the last base of `start` to the first base of `end`,
        /// if both are present.
        std::uint32_t length;

        /// Number of additional snarls included in the subchain to keep reversals
        /// within the subchain.
        std::uint32_t extra_snarls;

        /// Returns `true` if the subchain has a start node.
        bool has_start() const { return (this->type == Haplotypes::Subchain::normal || this->type == Haplotypes::Subchain::suffix); }

        /// Returns `true` if the subchain has an end node.
        bool has_end() const { return (this->type == Haplotypes::Subchain::normal || this->type == Haplotypes::Subchain::prefix); }
    };

    /// Creates a new `HaplotypePartitioner` using the given indexes.
    HaplotypePartitioner(const gbwtgraph::GBZ& gbz,
        const gbwt::FastLocate& r_index,
        const SnarlDistanceIndex& distance_index,
        const minimizer_index_type& minimizer_index,
        Verbosity verbosity);

    /// Parameters for `partition_haplotypes()`.
    struct Parameters {
        /// Target length for subchains (in bp).
        size_t subchain_length = Haplotypes::SUBCHAIN_LENGTH;

        /// Generate approximately this many jobs.
        size_t approximate_jobs = APPROXIMATE_JOBS;

        /// Avoid placing subchain boundaries in places where haplotypes would
        /// cross them multiple times.
        bool linear_structure = false;

        /// Represent the haplotypes with subchain kmers (see `Haplotypes`) instead of
        /// minimizers with a single occurrence in the graph.
        bool subchain_kmers = false;

        /// Print a description of the parameters.
        void print(std::ostream& out) const;
    };

    /**
     * Creates a `Haplotypes` representation of the haplotypes in the GBWT index.
     *
     * Top-level chains (weakly connected components in the graph) are assigned to
     * a number of jobs that can be later used as GBWT construction jobs. The
     * subchains of all chains are built as one set of tasks that OpenMP threads
     * run in parallel.
     *
     * Each top-level chain is partitioned into subchains that consist of one or
     * more snarls. Multiple snarls are combined into the same subchain if the
     * minimum distance over the subchain is at most the target length and there
     * are GBWT haplotypes that cross the subchain.
     *
     * With the right option, we keep extending the subchain if a haplotype would
     * cross the end in both directions. By doing this, we can avoid sequence loss
     * with haplotypes reversing their direction, while keeping kmers specific to
     * each subchain.
     *
     * If there are no snarls in a top-level chain, it is represented as a single
     * subchain without boundary nodes.
     *
     * Haplotypes crossing each subchain are represented using minimizers with a
     * single occurrence in the graph, or with subchain kmers (see `Parameters`).
     *
     * Throws `std::runtime_error` on error in single-threaded parts and exits
     * with `std::exit(EXIT_FAILURE)` in multi-threaded parts.
     */
    Haplotypes partition_haplotypes(const Parameters& parameters) const;

    const gbwtgraph::GBZ& gbz;
    gbwt::FragmentMap fragment_map;
    const gbwt::FastLocate& r_index;
    const SnarlDistanceIndex& distance_index;
    const minimizer_index_type& minimizer_index;

    Verbosity verbosity;

private:
    // Return the minimum distance from the last base of `from` to the first base of `to`.
    size_t get_distance(handle_t from, handle_t to) const;

    // Returns true if a haplotype visits the node in both orientations.
    bool contains_reversals(handle_t handle) const;

    // Partition the top-level chain into subchains.
    std::vector<Subchain> get_subchains(const gbwtgraph::TopLevelChain& chain, const Parameters& parameters) const;

    // Return (DA[i], i) for all GBWT sequences visiting a handle, sorted by sequence id
    // and the rank of the visit for the same sequence.
    std::vector<sequence_type> get_sequences(handle_t handle) const;

    // Get all GBWT sequences crossing the subchain.
    //
    // * If the subchain is a prefix (suffix), the sequences will be at the end
    //   (start) of the subchain.
    // * If the subchain is normal, the sequences will be at the start and
    //   correspond to minimal end-to-end visits to the subchain. A sequence
    //   that ends within the subchain may be selected if subsequent fragments
    //   of the same haplotype remain within the subchain and reach the end.
    std::vector<sequence_type> get_sequences(Subchain subchain) const;

    // The numbers of subchains with fragmented sequences, of fragmented sequences, and
    // of the fragments past the first in them.
    struct FragmentCounts {
        size_t subchains = 0;
        size_t sequences = 0;
        size_t additional = 0;

        void add(const FragmentCounts& another) {
            this->subchains += another.subchains;
            this->sequences += another.sequences;
            this->additional += another.additional;
        }
    };

    // Return the subchains of the top-level chain with the sequences crossing each. A
    // subchain that no haplotype crosses becomes a suffix and a prefix, and a chain
    // without snarls becomes a subchain of full haplotypes.
    std::vector<std::pair<Subchain, std::vector<sequence_type>>> chain_subchains(
        const gbwtgraph::TopLevelChain& chain, const std::string& contig_name, const Parameters& parameters
    ) const;

    // Build the output subchain from the sequences crossing the subchain. The kmers are
    // subchain kmers or minimizers with a single occurrence in the graph
    // (`Parameters::subchain_kmers`); a subchain of full haplotypes always uses the
    // latter. Returns the counts of its fragmented sequences.
    FragmentCounts build_subchain(
        const Subchain& subchain, const std::vector<sequence_type>& sequences, Haplotypes::Subchain& output,
        const Parameters& parameters
    ) const;
};

//------------------------------------------------------------------------------

/**
 * A class that creates synthetic haplotypes from a `Haplotypes` representation of
 * local haplotypes.
 */
class Recombinator {
public:
    /// Number of haplotypes to be generated.
    constexpr static size_t NUM_HAPLOTYPES = 4;

    /// A reasonable number of candidates for diploid sampling.
    constexpr static size_t NUM_CANDIDATES = 32;

    // TODO: Proper threshold?
    /// Badness threshold for subchains.
    constexpr static double BADNESS_THRESHOLD = 4.0;

    /// Expected kmer coverage. Use 0 to estimate from kmer counts.
    constexpr static size_t COVERAGE = 0;

    /// Multiplier to the score of a present kmer every time a haplotype with that
    /// kmer is selected.
    constexpr static double PRESENT_DISCOUNT = 0.9;

    /// Adjustment to the score of a heterozygous kmer every time a haplotype with
    /// (-) or without (+) that kmer is selected.
    constexpr static double HET_ADJUSTMENT = 0.05;

    /// Score for getting an absent kmer right/wrong. This should be less than 1, if
    /// we assume that having the right variants in the graph is more important than
    /// keeping wrong variants out.
    constexpr static double ABSENT_SCORE = 0.8;

    /// Threshold for homopolymer length in kmers that should be ignored when
    /// sampling with long reads.
    constexpr static size_t LONG_READ_HOMOPOLYMERS = 8;

    /// Minimum length (in bp) for a piece of a gref fragment to survive in the
    /// sampled graph. The same floor as `vg paths --min-gref-len`.
    constexpr static size_t MIN_GREF_LENGTH = 50;

    /// With the model pair, pair scores (log-likelihoods) are divided by this
    /// before they are combined with the transition probabilities.
    constexpr static double MODEL_TAU = 3.0;

    /// With the model pair, the rate per bp at which a member of the pair moves to
    /// another haplotype.
    constexpr static double MODEL_RHO = 3e-6;

    /// Number of haplotypes generated with the model preset.
    constexpr static size_t MODEL_NUM_HAPLOTYPES = 32;

    /// Continuity bonus with the model preset.
    constexpr static double MODEL_CONTINUITY_BONUS = 2.0;

    /// New-bp limit with the model preset.
    constexpr static size_t MODEL_MAX_NEW_BP = 1000;

    /// Escape from the new-bp limit with the model preset.
    constexpr static double MODEL_BP_ESCAPE = 300.0;

    /// N-base limit with the model preset.
    constexpr static size_t MODEL_MAX_N_BP = 0;

    /// The amount of progress information that should be printed to stderr.
    typedef Haplotypes::Verbosity Verbosity;

    /// A GBWT sequence as (sequence identifier, offset in a node).
    typedef Haplotypes::sequence_type sequence_type;

    /// Statistics on the generated haplotypes.
    struct Statistics {
        /// Number of top-level chains.
        size_t chains = 0;

        /// Number of subchains.
        size_t subchains = 0;

        /// Number of subchains exceeding the badness threshold.
        size_t bad_subchains = 0;

        /// Total number of fragments in the generated haplotypes.
        size_t fragments = 0;

        /// Number of top-level chains where full haplotypes were taken.
        /// These are not counted as fragments.
        size_t full_haplotypes = 0;

        /// Number of haplotypes generated.
        size_t haplotypes = 0;

        /// Number of additional haplotype fragments in bad subchains.
        size_t extra_fragments = 0;

        /// Number of times the same haplotype was extended from a subchain to the next subchain.
        size_t connections = 0;

        /// Number of reference paths included.
        size_t ref_paths = 0;

        /// Number of paths copied verbatim from excluded chains.
        size_t copied_paths = 0;

        /// Number of kmers selected.
        size_t kmers = 0;

        /// Total score for selected sequences.
        double score = 0.0;

        /// Combines the statistics into this object.
        void combine(const Statistics& another);

        /// Prints the statistics and returns the output stream.
        std::ostream& print(std::ostream& out) const;
    };

    /// Creates a new `Recombinator`.
    Recombinator(const gbwtgraph::GBZ& gbz, const Haplotypes& haplotypes, Verbosity verbosity);

    /// Parameters for `generate_haplotypes()`.
    struct Parameters {
        /// Number of haplotypes to be generated, or the number of candidates
        /// for diploid sampling.
        size_t num_haplotypes = NUM_HAPLOTYPES;

        /// Kmer coverage. Use 0 to estimate from kmer counts.
        size_t coverage = COVERAGE;

        /// Buffer size (in nodes) for GBWT construction.
        gbwt::size_type buffer_size = gbwt::DynamicGBWT::INSERT_BATCH_SIZE;

        /// Multiplicative factor for discounting the scores for present kmers after
        /// selecting a haplotype with that kmer.
        double present_discount = PRESENT_DISCOUNT;

        /// Additive term for adjusting the scores for heterozygous kmers after
        /// each haplotype to encourage even sampling of haplotypes with and without
        /// that kmer.
        double het_adjustment = HET_ADJUSTMENT;

        /// Score for absent kmers. This should be less than 1 if we assume that
        /// having the right variants in the graph is more important than keeping
        /// the wrong variants out.
        double absent_score = ABSENT_SCORE;

        /// Ignore kmers containing homopolymers of at least this length.
        size_t ignore_homopolymers = std::numeric_limits<size_t>::max();

        /// Score added to a candidate in a subchain if a generated haplotype took the
        /// same GBWT sequence in the previous subchain. This helps the generated
        /// haplotypes continue with the sequences they took.
        double continuity_bonus = 0.0;

        /// The new node sequence of a candidate is the total length of the nodes on its
        /// path through the subchain that no selected haplotype visits. After the first
        /// haplotype, a candidate is selected only if its new node sequence is at most
        /// this many bp, or if it passes `bp_escape`. When no candidate qualifies, the
        /// selected haplotypes are repeated.
        size_t max_new_bp = std::numeric_limits<size_t>::max();

        /// A candidate over `max_new_bp` still qualifies if it carries at least one new
        /// sample kmer (present or heterozygous, and carried by no selected haplotype)
        /// per this many bp of new node sequence. With 0, no candidate over the limit
        /// qualifies.
        double bp_escape = 0.0;

        /// A candidate is selected only if the nodes on its path through the subchain
        /// have at most this many N bases. If no candidate qualifies for the first
        /// haplotype, the best candidate is taken; if none qualifies later, the
        /// selected haplotypes are repeated.
        size_t max_n_bp = std::numeric_limits<size_t>::max();

        /// Returns true if there is a new-bp limit (`max_new_bp` is not the largest value).
        bool has_new_bp_limit() const { return this->max_new_bp != std::numeric_limits<size_t>::max(); }

        /// Returns true if there is an N-base limit (`max_n_bp` is not the largest value).
        bool has_n_base_limit() const { return this->max_n_bp != std::numeric_limits<size_t>::max(); }

        /// Match generated haplotypes with the haplotypes selected in a subchain by
        /// similarity. A generated haplotype first continues with the GBWT sequence it
        /// took in the previous subchain, if that sequence was selected. With this
        /// set, a generated haplotype whose GBWT sequence in the previous subchain
        /// also crosses this subchain then takes the remaining selection whose kmer
        /// presence in this subchain differs least from that sequence. The other
        /// generated haplotypes take the remaining selections in order.
        bool similar_threading = false;

        /// Select the first two haplotypes in each subchain with a model of the kmer
        /// counts and a hidden Markov model over pairs of haplotypes along each chain
        /// (see `recombinator_model.hpp`). A sequence with more than `max_n_bp` N
        /// bases is not a candidate for the pair, unless every sequence that can be
        /// selected has that many. The pair is exempt from `max_new_bp`. Chains with
        /// a non-standard scoring model do not use the model. It needs at least two
        /// haplotypes, and it cannot be combined with diploid sampling or haploid
        /// scoring.
        bool model_pair = false;

        /// Use the haploid scoring model. The most common kmer count is used as
        /// the coverage estimate. Kmers that would be classified as heterozygous
        /// are treated as homozygous.
        bool haploid_scoring = false;

        /// After selecting the initial `num_haplotypes` haplotypes, choose the
        /// highest-scoring pair out of them.
        bool diploid_sampling = false;

        /// Include named and reference paths. Gref fragments are clipped to the sampled
        /// graph, as copying them would keep the whole gref cover.
        bool include_reference = false;

        /// When using diploid sampling, include the remaining candidates as
        /// additional fragments in bad subchains.
        bool extra_fragments = false;

        /// Badness threshold for subchains when using diploid sampling.
        double badness_threshold = BADNESS_THRESHOLD;

        /// Minimum length (in bp) for a piece of a gref fragment to be kept.
        size_t min_gref_length = MIN_GREF_LENGTH;

        /// Samples whose haplotypes shouldn't be used, even if they score well.
        std::unordered_set<std::string> banned_samples;

        /// Kmer scoring model used for a chain.
        enum scoring_model_t {
            /// Standard model: absent, heterozygous, present, and frequent bands
            /// are scored according to the diploid/haploid coverage expectation.
            standard_scoring,
            /// High-coverage model. Kmers in the frequent category (count above
            /// the homozygous threshold) contribute the present score, while all
            /// other kmers contribute the absent score. This is intended for
            /// contigs such as chrM, where the true signal is in the frequent
            /// component and the haploid/diploid peaks are contamination (e.g.
            /// NuMTs or recurrent errors).
            high_coverage_scoring,
            /// Half-coverage model for heterogametic allosomes. The single true
            /// copy sits at ~cov/2, which the standard model labels heterozygous,
            /// so that band is rewarded as present. There is no real heterozygous
            /// component outside the PAR; the homozygous (~cov) band is paralog /
            /// contamination in the body and non-discriminative backbone in the
            /// PAR, so it is treated as uninformative (like the frequent band).
            half_coverage_scoring
        };

        /// Scoring model for the chain currently being processed. Set per chain
        /// from `high_coverage_chains` / `half_coverage_chains` before scoring.
        /// Diploid sampling is not used with the non-standard models.
        scoring_model_t scoring_model = standard_scoring;

        /// Top-level chains (by offset) to sample using the high-coverage model.
        std::unordered_set<size_t> high_coverage_chains;

        /// Number of haplotypes to generate for high-coverage chains.
        size_t high_coverage_num_haplotypes = NUM_HAPLOTYPES;

        /// Top-level chains (by offset) to sample using the half-coverage model.
        std::unordered_set<size_t> half_coverage_chains;

        /// Number of haplotypes to generate for half-coverage chains.
        size_t half_coverage_num_haplotypes = 2;

        /// Top-level chains (by offset) to copy through verbatim instead of
        /// personalizing. The reference/generic path and all haplotypes in the
        /// chain are preserved.
        std::unordered_set<size_t> excluded_chains;

        /// Contig names whose origin (first) fragment of every generated
        /// haplotype should be doubled, creating a self-loop that wraps the end
        /// of the sequence back onto its start. Intended for circular contigs
        /// such as chrM. Matched against the top-level chain contig name.
        std::unordered_set<std::string> wrap_contigs;

        /// Preset parameters for common use cases.
        enum preset_t : size_t {
            /// Default parameters.
            preset_default,
            /// Best practices for haploid sampling.
            preset_haploid,
            /// Best practices for diploid sampling.
            preset_diploid,
            /// Best practices for diploid sampling with HiFi reads.
            preset_hifi,
            /// Best practices for diploid sampling with ONT reads.
            preset_ont,
            /// The model pair and the selection limits, for an index built with
            /// subchain kmers.
            preset_model
        };

        /// A list of presets and their names.
        const static std::vector<std::pair<preset_t, std::string>> PRESETS;

        /// Returns the set of available presets as a string.
        static std::string preset_list(std::string_view separator);

        /// Parses a preset from its name.
        static std::optional<preset_t> parse_preset(std::string_view name);

        explicit Parameters(preset_t preset = preset_default);

        /// Print a description of the parameters.
        void print(std::ostream& out) const;
    };

    /**
     * Generates haplotypes based on the kmer counts in the given KFF file.
     *
     * Runs multiple GBWT construction jobs in parallel using OpenMP threads and
     * generates the specified number of haplotypes in each top-level chain
     * (component).
     *
     * Each generated haplotype has a single source haplotype in each subchain.
     * The source haplotype may consist of multiple fragments. Subchains are
     * by unary paths. Suffix / prefix subchains in the middle of a chain create
     * fragment breaks in every haplotype. If the chain starts without a prefix
     * (ends without a suffix), the haplotype chosen for the first (last)
     * subchain is used from the start (continued until the end).
     *
     * Throws `std::runtime_error` on error in single-threaded parts and exits
     * with `std::exit(EXIT_FAILURE)` in multi-threaded parts.
     */
    gbwt::GBWT generate_haplotypes(const std::string& kff_file, const Parameters& parameters) const;

    /// A local haplotype sequence within a single subchain.
    struct LocalHaplotype {
        /// Name of the haplotype.
        std::string name;

        /// Sequence in forward orientation.
        std::string sequence;

        /// (rank, score) in each round of haplotype selection this haplotype
        /// participates in.
        std::vector<std::pair<size_t, double>> scores;
    };

    /// What the haplotype selection in a subchain uses beyond the kmers.
    struct SelectionContext {
        /// For following a haplotype past the end of a fragment.
        const gbwt::FragmentMap& fragment_map;

        /// For each generated haplotype, the GBWT sequence it took in the previous
        /// subchain, or `gbwt::invalid_sequence()`.
        const std::vector<gbwt::size_type>* sources = nullptr;

        /// The pair chosen by the model, as sequence offsets in the subchain.
        const std::pair<size_t, size_t>* pair = nullptr;
    };

    const gbwtgraph::GBZ& gbz;
    const Haplotypes& haplotypes;
    gbwt::FragmentMap fragment_map;
    Verbosity verbosity;

    // A Haplotypes object contains a mapping from path ids to job ids.
    // This is a subset of the mapping for path handles / cached path offsets
    // corresponding to generic / reference paths in the current graph.
    // If the path is empty, the job id is haplotypes.jobs().
    std::vector<size_t> jobs_for_cached_paths;

private:
    // Return the pairs chosen by the count model (`Parameters::model_pair`) for each
    // chain, indexed by chain offset. The count model is fitted to the sample's counts
    // of the kmers in the chains sampled with the standard scoring model. Each of those
    // chains gets a pair of sequence offsets for each subchain, and the other chains get
    // no pairs. Throws `std::runtime_error` if the coverage is not positive.
    std::vector<std::vector<std::pair<size_t, size_t>>> model_pairs(
        const hash_map<Haplotypes::Subchain::kmer_type, size_t>& kmer_counts, double coverage, const Parameters& parameters
    ) const;

    // The haplotypes that `select_haplotypes()` selects in a subchain, and the kmers it
    // adds to `Statistics::kmers`.
    struct SubchainSelection {
        std::vector<std::pair<size_t, double>> haplotypes;
        size_t kmers = 0;
    };

    // Return the parameters for sampling the chain. High- and half-coverage chains use
    // their own scoring model, no diploid sampling, and their own haplotype count, which
    // are set in `chain_override` for them.
    const Parameters& chain_parameters(size_t chain_id, const Parameters& parameters, Parameters& chain_override) const;

    // Select the haplotypes in the subchains of the chains that are sampled, all
    // subchains in parallel, largest first. This is possible where the selection in a
    // subchain depends only on the subchain: without a continuity bonus. Returns the
    // selections by chain offset and subchain, and no selections for the other chains.
    std::vector<std::vector<SubchainSelection>> select_in_parallel(
        const hash_map<Haplotypes::Subchain::kmer_type, size_t>& kmer_counts, double coverage,
        const Parameters& parameters, const std::vector<std::vector<std::pair<size_t, size_t>>>& pairs
    ) const;

    // Generate haplotypes for the given chain. The model pairs are for the subchains
    // of the chain, or empty. The selections are those of `select_in_parallel()` for the
    // chain, or empty if the chain selects as it goes.
    Statistics generate_haplotypes(const Haplotypes::TopLevelChain& chain,
        const hash_map<Haplotypes::Subchain::kmer_type, size_t>& kmer_counts,
        gbwt::GBWTBuilder& builder, gbwtgraph::MetadataBuilder& metadata,
        const Parameters& parameters, double coverage,
        const std::vector<std::pair<size_t, size_t>>& pairs,
        const std::vector<SubchainSelection>& selections) const;

    // Copy the given chain through verbatim, preserving the reference/generic
    // path and all haplotypes without personalization.
    Statistics copy_chain(const Haplotypes::TopLevelChain& chain,
        gbwt::GBWTBuilder& builder, gbwtgraph::MetadataBuilder& metadata) const;
};

//------------------------------------------------------------------------------

/// Throw `std::runtime_error` if the parameters are invalid.
void validate_recombinator_parameters(const Recombinator::Parameters& parameters);

/**
 * Select haplotypes in a subchain. Return (sequence offset in the subchain, score)
 * for each selected haplotype. Haplotypes of banned samples and gref covers are never
 * selected.
 *
 * If the context has a pair, it becomes the first two haplotypes, and the function
 * throws `std::runtime_error` if a member of the pair cannot be selected. The rest are
 * selected greedily by kmer scores within the limits in the parameters. With diploid
 * sampling, the result is the best pair of the selected haplotypes, possibly followed by
 * extra fragments. Otherwise the selected haplotypes are repeated until there are
 * `parameters.num_haplotypes` of them.
 *
 * Update the statistics with the number of kmers used, if provided. Update the
 * local haplotypes with ranks and scores in each round of selection, if provided.
 */
std::vector<std::pair<size_t, double>> select_haplotypes(
    const gbwtgraph::GBZ& gbz,
    const Haplotypes::Subchain& subchain,
    const hash_map<Haplotypes::Subchain::kmer_type, size_t>& kmer_counts,
    double coverage, size_t k,
    Recombinator::Statistics* statistics,
    std::vector<Recombinator::LocalHaplotype>* local_haplotypes,
    const Recombinator::Parameters& parameters,
    const Recombinator::SelectionContext& context
);

/**
 * Match the unmatched generated haplotypes with the unused selections by similarity
 * (`Recombinator::Parameters::similar_threading`). A generated haplotype is
 * unmatched if `haplotype_to_selected` maps it to `sources.size()`, and selection
 * `i` is unused if `selected_in_use[i]` is 0.
 * A generated haplotype whose source (its GBWT sequence in the previous subchain) is
 * a sequence of the subchain can take an unused selection. Over all such pairs, the
 * pairs are taken in order of the number of kmers whose presence differs between the
 * source and the selected sequence, then by generated haplotype, and then by
 * sequence offset. A pair is taken if the generated haplotype is still unmatched and
 * an unused selection of that sequence remains; the selection with the lowest index
 * is used.
 */
void thread_by_similarity(
    const Haplotypes::Subchain& subchain,
    const std::vector<std::pair<size_t, double>>& selected_haplotypes,
    const std::vector<gbwt::size_type>& sources,
    std::vector<size_t>& haplotype_to_selected, sdsl::bit_vector& selected_in_use
);

/**
 * Return the candidates for the model pair in a subchain
 * (`Recombinator::Parameters::model_pair`), with pair scores from the table made
 * by `pair_model::count_table()`. The model uses the sequences that can be selected,
 * except those with more than `parameters.max_n_bp` N bases, unless every sequence
 * that can be selected has that many. If no sequence can be selected, the model uses
 * every sequence, because the hidden Markov model needs candidates in every subchain;
 * the N-base rule then applies to all of them. The candidates are the groups that
 * `pair_model::score_pairs()` keeps among the sequences the model uses. The haplotype
 * of a sequence is the sample and the phase of its path. The length of the subchain is
 * the median bases of the walks of the sequences, and the N bases of a sequence are those
 * of its walk (`Haplotypes::Subchain::median_walk_bases`). They are read from the subchain
 * if it has copies, and the sequences are walked otherwise.
 */
pair_model::Candidates model_candidates(
    const gbwtgraph::GBZ& gbz, const Haplotypes::Subchain& subchain,
    const hash_map<Haplotypes::Subchain::kmer_type, size_t>& kmer_counts, const pair_model::CountTable& table,
    const gbwt::FragmentMap& fragment_map, size_t k, const Recombinator::Parameters& parameters
);


//------------------------------------------------------------------------------

/// A maximal piece of a path whose nodes and edges all exist in a GBWT index.
struct PathPiece {
    /// Offset of the first node in the source path.
    size_t start = 0;

    /// Number of nodes in the piece.
    size_t nodes = 0;

    /// Offset (in bp) of the piece from the start of the source path.
    size_t bp_offset = 0;

    /// Length of the piece in bp.
    size_t bp_length = 0;
};

/**
 * Splits `path` into maximal pieces whose nodes and edges are all in `index`, so that
 * inserting a piece adds no nodes or edges. Offsets count every node of the path;
 * `node_length` gives the length of a node in bp.
 */
std::vector<PathPiece> clip_path_to_index(
    const gbwt::vector_type& path, const gbwt::DynamicGBWT& index,
    const std::function<size_t(gbwt::node_type)>& node_length
);

//------------------------------------------------------------------------------

} // namespace vg

#endif // VG_RECOMBINATOR_HPP_INCLUDED
