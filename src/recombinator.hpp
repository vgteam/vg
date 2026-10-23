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
#include "snarl_distance_index.hpp"

#include <functional>
#include <iostream>

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
    typedef gbwtgraph::MinimizerIndex<gbwtgraph::Key64> minimizer_index_type;

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

        /// Print a description of the parameters.
        void print(std::ostream& out) const;
    };

    /**
     * Creates a `Haplotypes` representation of the haplotypes in the GBWT index.
     *
     * Top-level chains (weakly connected components in the graph) are assigned to
     * a number of jobs that can be later used as GBWT construction jobs. Multiple
     * jobs are run in parallel using OpenMP threads.
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
     * single occurrence in the graph.
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

    // Return the sorted set of kmers that are minimizers in the sequence and have
    // a single occurrence in the graph.
    std::vector<kmer_type> unique_minimizers(
        gbwt::size_type sequence_id,
        gbwt::CachedGBWT& cached_gbwt,
        vg::hash_map<kmer_type, size_t>& kmer_cache
    ) const;

    // Returns the sorted set of kmers that are minimizers in the sequence over the
    // subchain and have a single occurrence in the graph. If the sequence does not
    // reach the end of the subchain, this will try to continue with the next fragment(s).
    //
    // Because we likely traverse the same nodes and see the same minimizers in
    // multiple haplotypes, we cache them at subchain level.
    //
    // Also reports the number of fragments that were used to generate the kmers.
    //
    // To avoid using kmers shared between all haplotypes in the subchain, and
    // potentially with neighboring subchains, this does not include kmers contained
    // entirely in the shared initial/final nodes.
    std::vector<kmer_type> unique_minimizers(
        sequence_type sequence, Subchain subchain, size_t& fragments,
        gbwt::CachedGBWT& cached_gbwt,
        vg::hash_map<kmer_type, size_t>& kmer_cache
    ) const;

    // Build subchains for a specific top-level chain.
    void build_subchains(const gbwtgraph::TopLevelChain& chain, Haplotypes::TopLevelChain& output, const Parameters& parameters) const;
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

    /// Minimum length (in bp) for a piece of a gref fragment to survive in the
    /// sampled graph. The same floor as `vg paths --min-gref-len`.
    constexpr static size_t MIN_GREF_LENGTH = 50;

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

        /// Use the haploid scoring model. The most common kmer count is used as
        /// the coverage estimate. Kmers that would be classified as heterozygous
        /// are treated as homozygous.
        bool haploid_scoring = false;

        /// After selecting the initial `num_haplotypes` haplotypes, choose the
        /// highest-scoring pair out of them.
        bool diploid_sampling = false;

        /// When using diploid sampling, include the remaining candidates as
        /// additional fragments in bad subchains.
        bool extra_fragments = false;

        /// Badness threshold for subchains when using diploid sampling.
        double badness_threshold = BADNESS_THRESHOLD;

        /// Include named and reference paths. Gref fragments are clipped to the sampled
        /// graph, as copying them would keep the whole gref cover.
        bool include_reference = false;

        /// Minimum length (in bp) for a piece of a gref fragment to be kept.
        size_t min_gref_length = MIN_GREF_LENGTH;

        /// Samples whose haplotypes shouldn't be used, even if they score well.
        unordered_set<std::string> banned_samples;

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

        // TODO: Should we use extra_fragments?
        /// Preset parameters for common use cases.
        enum preset_t {
            /// Default parameters.
            preset_default,
            /// Best practices for haploid sampling.
            preset_haploid,
            /// Best practices for diploid sampling.
            preset_diploid
        };

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
    // Generate haplotypes for the given chain.
    Statistics generate_haplotypes(const Haplotypes::TopLevelChain& chain,
        const hash_map<Haplotypes::Subchain::kmer_type, size_t>& kmer_counts,
        gbwt::GBWTBuilder& builder, gbwtgraph::MetadataBuilder& metadata,
        const Parameters& parameters, double coverage) const;

    // Copy the given chain through verbatim, preserving the reference/generic
    // path and all haplotypes without personalization.
    Statistics copy_chain(const Haplotypes::TopLevelChain& chain,
        gbwt::GBWTBuilder& builder, gbwtgraph::MetadataBuilder& metadata) const;
};

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
