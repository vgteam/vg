#ifndef VG_RECOMBINATOR_HPP_INCLUDED
#define VG_RECOMBINATOR_HPP_INCLUDED

/** \file recombinator.hpp
 * Tools for generating synthetic haplotypes as recombinations of existing
 * haplotypes.
 */

#include "gbwt_helper.hpp"
#include "gbwtgraph_helper.hpp"
#include "hash_map.hpp"
#include "snarl_distance_index.hpp"

#include <iostream>

#include <gbwtgraph/algorithms.h>

namespace vg {

//------------------------------------------------------------------------------

/**
 * A representation of the haplotypes in a graph.
 *
 * The graph is partitioned into top-level chains, which are further partitioned
 * into subchains. Each subchain contains a set of kmers and a collection of
 * sequences. Each sequence is defined by a bitvector marking the kmers that are
 * present.
 *
 * At the moment, the kmers are minimizers with a single occurrence in the graph.
 * The requirement is that each kmer is specific to a single subchain and does
 * not occur anywhere else in either orientation. (If no haplotype crosses a
 * snarl, that snarl is broken into a suffix and a prefix, and those subchains
 * may share kmers.)
 *
 * NOTE: This assumes that the top-level chains are linear, not cyclical.
 *
 * Versions:
 *
 * * Version 5: Every path in the graph is assigned to a construction job.
 *   This allows including reference paths that do not visit any snarls in the
 *   sampled graph. Not compatible with earlier versions.
 *
 * * Version 4: Subchains can have fragmented haplotypes instead of a single
 *   GBWT sequence always crossing from start to end. Compatible with version 3.
 *
 * * Version 3: Subchains use smaller integers when possible. Compatible with
 *   version 2.
 *
 * * Version 2: Top-level chains include a contig name. Compatible with version 1.
 *
 * * Version 1: Initial version.
 */
class Haplotypes {
public:
    /// The amount of progress information that should be printed to stderr.
    enum Verbosity : size_t {
        /// No progress information.
        verbosity_silent = 0,

        /// Basic information.
        verbosity_basic = 1,

        /// Basic information and detailed statistics.
        verbosity_detailed = 2,

        /// Basic information, detailed statistics, and debug information.
        verbosity_debug = 3,

        /// Hidden level; potentially tens of thousands of lines of debugging information.
        verbosity_extra_debug = 4
    };

    /// Header of the serialized file.
    struct Header {
        constexpr static std::uint32_t MAGIC_NUMBER = 0x4C504148; // "HAPL"
        constexpr static std::uint32_t VERSION = 5;
        constexpr static std::uint32_t MIN_VERSION = 5;
        constexpr static std::uint64_t DEFAULT_K = 29;

        /// A magic number that identifies the file.
        std::uint32_t magic_number = MAGIC_NUMBER;

        /// Version of the file.
        std::uint32_t version = VERSION;

        /// Number of top-level chains in the graph.
        std::uint64_t top_level_chains = 0;

        /// Number of GBWT construction jobs for the chains.
        std::uint64_t construction_jobs = 0;

        /// Total number of subchains in all chains.
        std::uint64_t total_subchains = 0;

        /// Total number of kmers in all subchains.
        std::uint64_t total_kmers = 0;

        /// Length of the kmers.
        std::uint64_t k = DEFAULT_K;
    };

    /// A GBWT sequence as (sequence identifier, offset in a node).
    typedef std::pair<gbwt::size_type, gbwt::size_type> sequence_type;

    /// A more space-efficient representation of `sequence_type`.
    typedef std::pair<std::uint32_t, std::uint32_t> compact_sequence_type;

    /// Representation of a subchain.
    struct Subchain {
        /// Subchain types.
        enum subchain_t : std::uint64_t {
            /// Normal subchain with two boundary nodes.
            normal = 0,

            /// A prefix with only an end node.
            prefix = 1,

            /// A suffix with only a start node.
            suffix = 2,

            /// A full haplotype with no boundary nodes.
            full_haplotype = 3
        };

        /// An encoded kmer.
        typedef gbwtgraph::Key64::value_type kmer_type;

        /// The type of this subchain.
        subchain_t type;

        /// Boundary nodes, or `gbwt::ENDMARKER` if not present.
        gbwt::node_type start, end;

        /// A vector of distinct kmers. For each kmer, list the kmer itself and the number
        /// of haplotypes it appears in.
        std::vector<kmer_type> kmers;

        /// Number of haplotypes each kmer appears in.
        sdsl::int_vector<0> kmer_counts;

        /// Sequences as (GBWT sequence id, offset in the relevant node).
        std::vector<compact_sequence_type> sequences;

        // TODO v6: Use an extra bit for each sequence to mark whether the presence for that sequence
        // is stored explicitly or relative to the last explicit sequence.
        // We need to cluster the sequences by similarity and store the clusters consecutively.
        // And then use sd_vector for the sequences with relative presence.
        // Decompress to a single bitvector when needed.
        /// A bit vector marking the presence of kmers in the sequences.
        /// Sequence `i` contains kmer `j` if and only if `kmers_present[i * kmers.size() + j] == 1`.
        sdsl::bit_vector kmers_present;

        /// Returns the start node as a GBWTGraph handle.
        handle_t start_handle() const { return gbwtgraph::GBWTGraph::node_to_handle(this->start); }

        /// Returns the end node as a GBWTGraph handle.
        handle_t end_handle() const { return gbwtgraph::GBWTGraph::node_to_handle(this->end); }

        /// Returns `true` if the subchain has a start node.
        bool has_start() const { return (this->type == normal || this->type == suffix); }

        /// Returns `true` if the subchain has an end node.
        bool has_end() const { return (this->type == normal || this->type == prefix); }

        /// Returns a string representation of the type and the boundary nodes.
        std::string to_string() const;

        /// Returns (sequence identifier, offset in a node) for the given sequence.
        sequence_type get_sequence(size_t i) const {
            return { this->sequences[i].first, this->sequences[i].second };
        }

        /// Returns the distance from the last base of `start` to the first base of
        /// `end` over the given sequence. Returns 0 if the subchain is not normal or
        /// if the sequence does not exist.
        size_t distance(const gbwtgraph::GBZ& gbz, size_t i) const;

        /// Returns an estimate of the badness of the subchain.
        /// The ideal value is 0.0, and higher values indicate worse subchains.
        /// The estimate is based on the following factors:
        /// * Length of the subchain.
        /// * Number of haplotypes relative to the expected number.
        /// * Information content of the kmers (disabled).
        double badness(const gbwtgraph::GBZ& gbz) const;

        /// Serializes the object to a stream in the Simple-SDS format.
        void simple_sds_serialize(std::ostream& out) const;

        /// Loads the object from a stream in the Simple-SDS format.
        void simple_sds_load(std::istream& in);

        /// Returns the size of the object in elements.
        size_t simple_sds_size() const;
    };

    /// Representation of a top-level chain.
    struct TopLevelChain {
        /// Offset in the child list of the root snarl.
        size_t offset;

        /// GBWT construction job for this chain.
        size_t job_id;

        /// Contig name corresponding to the chain.
        std::string contig_name;

        /// Subchains in the order they appear in.
        std::vector<Subchain> subchains;

        /// Serializes the object to a stream in the Simple-SDS format.
        void simple_sds_serialize(std::ostream& out) const;

        /// Loads the object from a stream in the Simple-SDS format.
        void simple_sds_load(std::istream& in);

        /// Returns the size of the object in elements.
        size_t simple_sds_size() const;
    };

    /// Returns the number of weakly connected components.
    size_t components() const { return this->header.top_level_chains; }

    /// Returns the number of GBWT construction jobs.
    size_t jobs() const { return this->header.construction_jobs; }

    /// Returns the length of the kmers.
    size_t k() const { return this->header.k; }

    /// Returns the number of kmers in the subchains.
    size_t kmers() const { return this->header.total_kmers; }

    Header header;

    // Job ids for each path in the GBWTGraph, or `jobs()` if the path is empty.
    std::vector<size_t> jobs_for_paths;

    std::vector<TopLevelChain> chains;

    /**
      * Returns a mapping from kmers to their counts in the given KFF file.
      * The counts include both the kmer and the reverse complement.
      *
      * Reads the KFF file using OpenMP threads. Exits with `std::exit()` if
      * the file cannot be opened and throws `std::runtime_error` if the kmer
      * counts cannot be used.
     */
    hash_map<Subchain::kmer_type, size_t> kmer_counts(const std::string& kff_file, Verbosity verbosity) const;

    /// Serializes the object to a stream in the Simple-SDS format.
    /// I/O errors can be detected by checking the stream state.
    void simple_sds_serialize(std::ostream& out) const;

    /// Serializes the object to a file in the Simple-SDS format.
    /// Prints an error message and exits the program on failure.
    void serialize_to(const std::string& filename) const;

    /// Loads the object from a stream in the Simple-SDS format.
    /// I/O errors can be detected by checking the stream state.
    /// Throws `sdsl::simple_sds::InvalidData` if data is unacceptable.
    void simple_sds_load(std::istream& in);

    /// Loads the object from a file in the Simple-SDS format.
    /// Prints an error message and exits the program on failure.
    void load_from(const std::string& filename);

    /// Returns the size of the object in elements.
    size_t simple_sds_size() const;

    /**
     * Assigns each reference and generic path in the graph to a GBWT construction job.
     *
     * For each path handle from 0 to gbz.named_paths() - 1, we assign the path to
     * the given construction job, or jobs() if the path is empty.
     */
    std::vector<size_t> assign_reference_paths(const gbwtgraph::GBZ& gbz, Verbosity verbosity) const;
};

//------------------------------------------------------------------------------

/**
 * A tool for transforming the haplotypes in a GBWT index into a `Haplotypes`
 * representation. Requires a GBZ graph, an r-index, a distance index, and a
 * minimizer index.
 */
class HaplotypePartitioner {
public:
    /// Target length of a subchain.
    constexpr static size_t SUBCHAIN_LENGTH = 10000;

    /// Approximate number of construction jobs to be created.
    constexpr static size_t APPROXIMATE_JOBS = 32;

    /// The amount of progress information that should be printed to stderr.
    typedef Haplotypes::Verbosity Verbosity;

    /// A GBWT sequence as (sequence identifier, offset in a node).
    typedef Haplotypes::sequence_type sequence_type;

    /// An encoded kmer.
    typedef Haplotypes::Subchain::kmer_type kmer_type;

    /// Minimizer index without payloads.
    typedef gbwtgraph::MinimizerIndex<gbwtgraph::Key64, gbwtgraph::Position> minimizer_index_type;

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
        size_t subchain_length = SUBCHAIN_LENGTH;

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
    std::vector<kmer_type> unique_minimizers(gbwt::size_type sequence_id) const;

    // Returns the sorted set of kmers that are minimizers in the sequence over the
    // subchain and have a single occurrence in the graph. If the sequence does not
    // reach the end of the subchain, this will try to continue with the next fragment(s).
    //
    // Also reports the number of fragments that were used to generate the kmers.
    //
    // To avoid using kmers shared between all haplotypes in the subchain, and
    // potentially with neighboring subchains, this does not include kmers contained
    // entirely in the shared initial/final nodes.
    std::vector<kmer_type> unique_minimizers(sequence_type sequence, Subchain subchain, size_t& fragments) const;

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

    /// Block size (in kmers) for reading KFF files.
    constexpr static size_t KFF_BLOCK_SIZE = 1000000;

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

        /// Include named and reference paths.
        bool include_reference = false;

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

    /// Kmer classification.
    enum kmer_presence { absent, heterozygous, present, frequent };

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
};

//------------------------------------------------------------------------------

} // namespace vg

#endif // VG_RECOMBINATOR_HPP_INCLUDED
