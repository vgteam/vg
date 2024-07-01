#ifndef VG_RECOMBINATOR_HPP_INCLUDED
#define VG_RECOMBINATOR_HPP_INCLUDED

/** \file 
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
        verbosity_debug = 3
    };

    /// Header of the serialized file.
    struct Header {
        constexpr static std::uint32_t MAGIC_NUMBER = 0x4C504148; // "HAPL"
        constexpr static std::uint32_t VERSION = 2;
        constexpr static std::uint32_t MIN_VERSION = 1;
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
        std::vector<std::pair<kmer_type, size_t>> kmers;

        // TODO: This could be smaller
        /// Sequences as (GBWT sequence id, offset in the relevant node).
        std::vector<sequence_type> sequences;

        // TODO: This needs to be compressed for larger datasets.
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

        /// Serializes the object to a stream in the simple-sds format.
        void simple_sds_serialize(std::ostream& out) const;

        /// Loads the object from a stream in the simple-sds format.
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

        /// Serializes the object to a stream in the simple-sds format.
        void simple_sds_serialize(std::ostream& out) const;

        /// Loads the object from a stream in the simple-sds format.
        void simple_sds_load(std::istream& in);

        /// Loads the old version without a contig name.
        void load_old(std::istream& in);

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

    // Job ids for each cached path in the GBWTGraph, or `jobs()` if the path is empty.
    std::vector<size_t> jobs_for_cached_paths;

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

    /// Serializes the object to a stream in the simple-sds format.
    void simple_sds_serialize(std::ostream& out) const;

    /// Loads the object from a stream in the simple-sds format.
    void simple_sds_load(std::istream& in);

    /// Returns the size of the object in elements.
    size_t simple_sds_size() const;
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

        /// Generate approximately this many  jobs.
        size_t approximate_jobs = APPROXIMATE_JOBS;
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
     * are GBWT haplotypes that cross the subchain. We also keep extending the
     * subchain if a haplotype would cross the end in both directions. By doing
     * this, we can avoid sequence loss with haplotypes reversing their direction,
     * while keeping kmers specific to each subchain.
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

    // Return (SA[i], i) for all GBWT sequences visiting a handle, sorted by sequence id
    // and the number of the visit.
    std::vector<sequence_type> get_sequence_visits(handle_t handle) const;

    // Return (DA[i], i) for all GBWT sequences visiting a handle, sorted by sequence id.
    std::vector<sequence_type> get_sequences(handle_t handle) const;

    // Get all GBWT sequences crossing the subchain. The sequences will be at
    // start for normal subchains and suffixes and at end for prefixes.
    std::vector<sequence_type> get_sequences(Subchain subchain) const;

    // Return the sorted set of kmers that are minimizers in the sequence and have
    // a single occurrence in the graph.
    std::vector<kmer_type> unique_minimizers(gbwt::size_type sequence_id) const;

    // Count the number of minimizers in the sequence over the subchain with a single
    // occurrence in the graph. Return the sorted set of kmers that are minimizers in
    // the sequence over the subchain and have a single occurrence in the graph.
    //
    // To avoid using kmers shared between all haplotypes in the subchain, and
    // potentially with neighboring subchains, this does not include kmers contained
    // entirely in the shared initial/final nodes.
    std::vector<kmer_type> unique_minimizers(sequence_type sequence, Subchain subchain) const;

    // Build subchains for a specific top-level chain.
    void build_subchains(const gbwtgraph::TopLevelChain& chain, Haplotypes::TopLevelChain& output, const Parameters& parameters) const;
};

//------------------------------------------------------------------------------

/*
  TODO:
    * Parameters should include the functionality for classifying kmers, scoring
      haplotypes, and samplign haplotypes.

  Models:
    * Diploid (current)
    * Haploid
*/

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

        /// Number of fragments.
        size_t fragments = 0;

        /// Number of top-level chains where full haplotypes were taken.
        size_t full_haplotypes = 0;

        /// Number of haplotypes generated.
        size_t haplotypes = 0;

        /// Number of times a haplotype was extended from a subchain to the next subchain.
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

        /// After selecting the initial `num_haplotypes` haplotypes, choose the
        /// highest-scoring pair out of them.
        bool diploid_sampling = false;

        /// Include named and reference paths.
        bool include_reference = false;
    };

    /**
     * Generates haplotypes based on the kmer counts in the given KFF file.
     *
     * Runs multiple GBWT construction jobs in parallel using OpenMP threads and
     * generates the specified number of haplotypes in each top-level chain
     * (component).
     *
     * Each generated haplotype has a single source haplotype in each subchain.
     * The subchains are connected by unary paths. Suffix / prefix subchains in
     * the middle of a chain create fragment breaks. If the chain starts without
     * a prefix (ends without a suffix), the haplotype chosen for the first (last)
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

    /**
     * Classifies the kmers used for describing the haplotypes according to
     * their frequency in the KFF file. Uses `A`, `H`, `P`, and `F` to represent
     * absent, heterozygous, present, and frequent kmers, respectively.
     *
     * Throws `std::runtime_error` on error.
     */
    std::vector<char> classify_kmers(const std::string& kff_file, const Parameters& parameters) const;

    /**
     * Extracts the local haplotypes in the given subchain. In addition to the
     * haplotype sequence, this also reports the name of the corresponding path
     * as well as (rank, score) for the haplotype in each round of haplotype
     * selection. The number of rounds is `parameters.num_haplotypes`, but if
     * the haplotype is selected earlier, it will not get further scores.
     *
     * Throws `std::runtime_error` on error.
     */
    std::vector<LocalHaplotype> extract_sequences(
        const std::string& kff_file, size_t chain_id, size_t subchain_id, const Parameters& parameters
    ) const;

    const gbwtgraph::GBZ& gbz;
    const Haplotypes& haplotypes;
    Verbosity verbosity;

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
