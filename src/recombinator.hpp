#ifndef VG_RECOMBINATOR_HPP_INCLUDED
#define VG_RECOMBINATOR_HPP_INCLUDED

/** \file 
 * Tools for generating synthetic haplotypes as recombinations of existing
 * haplotypes.
 */

#include "gbwt_helper.hpp"
#include "gbwtgraph_helper.hpp"
#include "snarl_distance_index.hpp"

#include <iostream>
#include <unordered_map>

#include <gbwtgraph/algorithms.h>

namespace vg {

//------------------------------------------------------------------------------

// FIXME tests
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
 * not occur anywhere else. (If no haplotype crosses a snarl, that snarl is
 * broken into a suffix and a prefix, and those subchains may share kmers.)
 */
class Haplotypes {
public:
    /// Header of the serialized file.
    struct Header {
        constexpr static std::uint32_t MAGIC_NUMBER = 0x4C504148; // "HAPL"
        constexpr static std::uint32_t VERSION = 1;
        constexpr static std::uint64_t DEFAULT_K = 29;

        /// A magic number that identifies the file.
        std::uint32_t magic_number = MAGIC_NUMBER;

        /// Version of the file.
        std::uint32_t version = VERSION;

        /// Number of top-level chains in the graph.
        std::uint64_t top_level_chains = 0;

        /// Number of GBWT construction jobs for the chains.
        std::uint64_t construction_jobs = 0;

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

        /// A vector of distinct kmers.
        std::vector<kmer_type> kmers;

        // TODO: This could be smaller
        /// Sequences as (GBWT sequence id, offset in the relevant node).
        std::vector<sequence_type> sequences;

        // TODO: This could be compressed
        /// Concatenated bitvectors for each sequence that mark the presence of each kmer
        /// in that sequence.
        sdsl::bit_vector kmers_present;

        /// Returns the start node as a GBWTGraph handle.
        handle_t start_handle() const { return gbwtgraph::GBWTGraph::node_to_handle(this->start); }

        /// Returns the end node as a GBWTGraph handle.
        handle_t end_handle() const { return gbwtgraph::GBWTGraph::node_to_handle(this->end); }

        /// Returns `true` if the subchain has a start node.
        bool has_start() const { return (this->type == normal || this->type == suffix); }

        /// Returns `true` if the subchain has an end node.
        bool has_end() const { return (this->type == normal || this->type == prefix); }

        /// Returns the type of the subchain as a string.
        std::string type_of() const;

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

        /// Subchains in the order they appear in.
        std::vector<Subchain> subchains;

        /// Serializes the object to a stream in the simple-sds format.
        void simple_sds_serialize(std::ostream& out) const;

        /// Loads the object from a stream in the simple-sds format.
        void simple_sds_load(std::istream& in);

        /// Returns the size of the object in elements.
        size_t simple_sds_size() const;
    };

    /// Kmer count for a kmer used for defining haplotypes in a subchain.
    struct KMerCount {
        /// Offset of the top-level chain in `Haplotypes::chains`.
        std::uint32_t chain;

        /// Offset of the subchain in `TopLevelChain::subchains`.
        std::uint32_t subchain;

        /// Offset of the kmer in `Subchain::kmers`.
        std::uint32_t kmer;

        /// Number of occurrences of the kmer.
        std::uint32_t count;
    };

    /// Returns the number of weakly connected components.
    size_t components() const { return this->header.top_level_chains; }

    /// Returns the number of GBWT construction jobs.
    size_t jobs() const { return this->header.construction_jobs; }

    /// Returns the length of the kmers.
    size_t k() const { return this->header.k; }

    Header header;
    std::vector<TopLevelChain> chains;

    /// Returns a mapping from canonical KFF kmers to their counts in the
    /// given file.
    ///
    /// A canonical KFF kmer is the minimum of a kmer and its reverse
    /// complement when the encoded kmers are interpreted as big-endian
    /// integers.
    ///
    /// Exits with `std::exit()` if the file cannot be opened and throws
    /// `std::runtime_error` if the kmer counts cannot be used.
    std::unordered_map<std::uint64_t, KMerCount> kmer_counts(const std::string& kff_file) const;

    /// Serializes the object to a stream in the simple-sds format.
    void simple_sds_serialize(std::ostream& out) const;

    /// Loads the object from a stream in the simple-sds format.
    void simple_sds_load(std::istream& in);

    /// Returns the size of the object in elements.
    size_t simple_sds_size() const;
};

//------------------------------------------------------------------------------

// FIXME tests
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

    /// A GBWT sequence as (sequence identifier, offset in a node).
    typedef Haplotypes::sequence_type sequence_type;

    /// An encoded kmer.
    typedef Haplotypes::Subchain::kmer_type kmer_type;

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

        /// Returns `true` if the subchain has a start node.
        bool has_start() const { return (this->type == Haplotypes::Subchain::normal || this->type == Haplotypes::Subchain::suffix); }

        /// Returns `true` if the subchain has an end node.
        bool has_end() const { return (this->type == Haplotypes::Subchain::normal || this->type == Haplotypes::Subchain::prefix); }
    };

    /// Creates a new `HaplotypePartitioner` using the given indexes.
    HaplotypePartitioner(const gbwtgraph::GBZ& gbz,
        const gbwt::FastLocate& r_index,
        const SnarlDistanceIndex& distance_index,
        const gbwtgraph::DefaultMinimizerIndex& minimizer_index,
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
     * are GBWT haplotypes that cross the subchain. If there are no snarls in a
     * top-level chain, it is represented as a single subchain without boundary
     * nodes.
     *
     * Haplotypes crossing each subchain are represented using minimizers with a
     * single occurrence in the graph.
     *
     * TODO: Are all chains in the right orientation?
     */
    Haplotypes partition_haplotypes(const Parameters& parameters) const;

    const gbwtgraph::GBZ& gbz;
    const gbwt::FastLocate& r_index;
    const SnarlDistanceIndex& distance_index;
    const gbwtgraph::DefaultMinimizerIndex& minimizer_index;

    Verbosity verbosity;

private:
    // Return the minimum distance from the last base of `from` to the first base of `to`.
    size_t get_distance(handle_t from, handle_t to) const;

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

// FIXME tests
/**
 * A class that creates synthetic haplotypes from a `Haplotypes` representation of
 * local haplotypes.
 */
class Recombinator {
public:
    /// Number of haplotypes to be generated.
    constexpr static size_t NUM_HAPLOTYPES = 16;

    /// A GBWT sequence as (sequence identifier, offset in a node).
    typedef Haplotypes::sequence_type sequence_type;

    /**
     * A haplotype beging generated as a GBWT path.
     *
     * GBWT metadata will be set as following:
     *
     * * Sample name is "recombination".
     * * Contig name is "chain_X", where X is the chain identifier.
     * * Haplotype identifier is set during construction.
     * * Fragment identifier is set as necessary.
     */
    struct Haplotype {
        /// Contig identifier in GBWT metadata.
        /// Offset of the top-level chain in the children of the root snarl.
        size_t chain;

        /// Haplotype identifier in GBWT metadata. 
        size_t id;

        /// Fragment identifier in GBWT metadata.
        /// If no original haplotype crosses a subchain, a new fragment will
        /// start after the subchain.
        size_t fragment;

        /// GBWT position at the end of the latest `extend()` call.
        /// `gbwt::invalid_edge()` otherwise.
        gbwt::edge_type position;

        /// The path being generated.
        gbwt::vector_type path;

        /**
         * Extends the haplotype over the given subchain by using the given
         * original haplotype.
         *
         * This assumes that the original haplotype crosses the subchain.
         *
         * If `extend()` has been called for this fragment, there must be a
         * unary path connecting the subchains, which will be used in the
         * generated haplotype.
         *
         * If `extend()` has not been called, the generated haplotype will
         * take the prefix of the original original haplotype until the start
         * of the subchain.
         */
        void extend(sequence_type sequence, const Haplotypes::Subchain& subchain, const Recombinator& recombinator, gbwt::GBWTBuilder& builder);

        /// Takes an existing haplotype from the GBWT index and inserts it into
        /// the builder. This is intended for fragments that do not contain
        /// subchains crossed by the original haplotypes. The call will fail if
        /// `extend()` has been called.
        void take(gbwt::size_type sequence_id, const Recombinator& recombinator, gbwt::GBWTBuilder& builder);

        /// Extends the original haplotype from the latest `extend()` call until
        /// the end, inserts it into the builder, and starts a new fragment.
        /// The call will fail if `extend()` has not been called for this
        /// fragment.
        void finish(const Recombinator& recombinator, gbwt::GBWTBuilder& builder);

    private:
        // Extends the haplotype over a unary path from a previous subchain.
        void connect(gbwt::node_type until, const gbwtgraph::GBWTGraph& graph);

        // Takes a prefix of a sequence.
        void prefix(gbwt::size_type sequence_id, gbwt::node_type until, const gbwt::GBWT& index);

        // Extends the haplotype from the previous subchain until the end.
        void suffix(const gbwt::GBWT& index);

        // Inserts the current fragment into the builder.
        void insert(gbwt::GBWTBuilder& builder);
    };

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

        /// Combines the statistics into this object.
        void combine(const Statistics& another);

        /// Prints the statistics and returns the output stream.
        std::ostream& print(std::ostream& out) const;
    };

    /// Creates a new `Recombinator`.
    Recombinator(const gbwtgraph::GBZ& gbz, HaplotypePartitioner::Verbosity verbosity);

    /// Parameters for `generate_haplotypes()`.
    struct Parameters {
        /// Number of haplotypes to be generated.
        size_t num_haplotypes = NUM_HAPLOTYPES;

        /// Buffer size (in nodes) for GBWT construction.
        gbwt::size_type buffer_size = gbwt::DynamicGBWT::INSERT_BATCH_SIZE;
    };

    // FIXME kmer input
    /**
     * Generates haplotypes based on the given `Haplotypes` representation.
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
     * TODO: Include reference paths?
     */
    gbwt::GBWT generate_haplotypes(const Haplotypes& haplotypes, const Parameters& parameters) const;

    const gbwtgraph::GBZ& gbz;
    HaplotypePartitioner::Verbosity verbosity;

private:
    // Generate haplotypes for the given chain.
    Statistics generate_haplotypes(const Haplotypes::TopLevelChain& chain, gbwt::GBWTBuilder& builder, const Parameters& parameters) const;
};

//------------------------------------------------------------------------------

} // namespace vg

#endif // VG_RECOMBINATOR_HPP_INCLUDED
