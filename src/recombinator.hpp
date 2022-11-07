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

#include <gbwtgraph/algorithms.h>

namespace vg {

//------------------------------------------------------------------------------

// FIXME tests, exceptions for error handling
/**
 * A class that generates synthetic haplotypes as recombinations of the haplotypes in
 * an existing GBWT index.
 *
 * Haplotypes will be generated separately for each weakly connected component in the
 * graph, which correspond to top-level chains in the snarl decomposition. Each
 * top-level chain is partitioned into subchains that start and end with a node. The
 * target length of a subchain is 10 kbp.
 *
 * Each generated haplotype has a single source haplotype in each subchain. The
 * subchains are connected by unary paths. If no haplotype crosses an entire subchain,
 * it creates a fragment break in the haplotype. In such cases, the source haplotype
 * for the previous subchain continues until its end, and the next fragment starts
 * from the beginning of the source haplotype for the next subchain.
 *
 * The generation of synthetic haplotypes can be parallelized using sets of top-level
 * chains as construction jobs.
 *
 * TODO: Parameterize number of haplotypes, subchain length, number of construction jobs.
 * TODO: Should the prefix and the suffix be from separate haplotypes?
 * TODO: Include reference paths?
 * TODO: Are all chains in the right orientation?
 */
class Recombinator {
public:
    /// Number of haplotypes to be generated.
    constexpr static size_t NUM_HAPLOTYPES = 16;

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
    typedef std::pair<gbwt::size_type, gbwt::size_type> sequence_type;

    /**
     * A subchain is a substring of a top-level chain defined by at one or two
     * boundary nodes.
     *
     * Normal subchains have two boundary nodes, which are assumed to be the
     * start node of a snarl and the end node of a possibly different snarl.
     * There are assumed to be haplotypes crossing the subchain. Prefixes and
     * suffixes lack one of the boundary nodes.
     *
     * When a top-level chain is partitioned into subchains, the boundary nodes
     * may either overlap or be connected by unary paths. If a snarl is not
     * connected, it may be presented as a suffix and a prefix.
     */
    struct Subchain {
        /// Type of a subchain.
        enum subchain_t { normal, prefix, suffix };

        /// The type of this subchain.
        subchain_t type;

        /// Start node.
        handle_t start;

        /// End node.
        handle_t end;

        /// Returns `true` if the subchain has a start node.
        bool has_start() const { return (this->type != prefix); }

        /// Returns `true` if the subchain has an end node.
        bool has_end() const { return (this->type != suffix); }
    };

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
        void extend(sequence_type sequence, Subchain subchain, const Recombinator& recombinator, gbwt::GBWTBuilder& builder);

        /// Takes an existing haplotype from the GBWT index and inserts it into
        /// the builder. This is intended for fragments that do not contain
        /// subchains crossed by the original haplotypes. The call will fail if
        /// `extend()` has been called.
        void take(gbwt::size_type sequence, const Recombinator& recombinator, gbwt::GBWTBuilder& builder);

        /// Extends the original haplotype from the latest `extend()` call until
        /// the end, inserts it into the builder, and starts a new fragment.
        /// The call will fail if `extend()` has not been called for this
        /// fragment.
        void finish(const Recombinator& recombinator, gbwt::GBWTBuilder& builder);

    private:
        // Extends the haplotype over a unary path from a previous subchain.
        void connect(handle_t until, const gbwtgraph::GBWTGraph& graph);

        // Takes a prefix of a sequence.
        void prefix(gbwt::size_type sequence, handle_t until, const gbwt::GBWT& index);

        // Extends the haplotype from the previous subchain until the end.
        void suffix(const gbwt::GBWT& index);

        // Inserts the current fragment into the builder.
        void insert(gbwt::GBWTBuilder& builder);
    };

    /// Statistics on the generated haplotypes.
    struct Statistics
    {
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

    /// FIXME: document
    Recombinator(const gbwtgraph::GBZ& gbz, const gbwt::FastLocate& r_index, const SnarlDistanceIndex& distance_index, Verbosity verbosity);

    // FIXME job()

    const gbwtgraph::GBZ&     gbz;
    const gbwt::FastLocate&   r_index;
    const SnarlDistanceIndex& distance_index;

    std::vector<std::vector<gbwtgraph::TopLevelChain>> chains_by_job;

    Verbosity verbosity;

    // FIXME should be private
    // Generate haplotypes for the given chain.
    Statistics generate_haplotypes(const gbwtgraph::TopLevelChain& chain, gbwt::GBWTBuilder& builder) const;

private:

    // Partition a chain into subchains.
    std::vector<Subchain> get_subchains(const gbwtgraph::TopLevelChain& chain) const;

    // Return the minimum distance over the candidate subchain.
    size_t get_distance(handle_t from, handle_t to) const;

    // Return all GBWT sequences visiting a handle.
    std::vector<sequence_type> get_sequences(handle_t handle) const;

    // Get all GBWT sequences crossing the subchain. The sequences will be at the
    // start for normal subchains and at the only boundary node for other subchains.
    std::vector<sequence_type> get_sequences(Subchain subchain) const;
};

//------------------------------------------------------------------------------

} // namespace vg

#endif // VG_RECOMBINATOR_HPP_INCLUDED
