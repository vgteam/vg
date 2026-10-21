#ifndef VG_RECOMBINATOR_HAPLOTYPES_HPP_INCLUDED
#define VG_RECOMBINATOR_HAPLOTYPES_HPP_INCLUDED

/** \file recombinator_haplotypes.hpp
 * Haplotype information files used in haplotype sampling.
 */

#include "gbwt_helper.hpp"
#include "gbwtgraph_helper.hpp"
#include "hash_map.hpp"

#include <functional>
#include <iostream>

namespace vg {

//------------------------------------------------------------------------------

// FIXME: Once we have determined that we sample the same haplotypes as with the
// old uncompressed kmer presence matrix, we could drop the mapping between
// the original order and the stored order and simply use the stored order
// in `Subchain`. (The existing order is already arbitrary.)
/**
 * A compressed binary matrix with haplotype sequences as rows and k-mers as columns.
 * Set bits indicate k-mers that are present in the corresponding sequence.
 *
 * The compression assumes that the matrix consists of clusters of similar sequences.
 * We choose one sequence from each cluster as a parent and store it explicitly.
 * For the other sequences (children), we only store the symmetric difference
 * with their parent sequence. The internal representation reorders the sequences so
 * that the parent of each cluster appears first, followed by its children.
 */
class KmerPresenceMatrix {
public:
    /// Kmer classification.
    enum kmer_classification : size_t { absent = 0, heterozygous = 1, present = 2, frequent = 3 };

    /// Default constructor that creates an empty k-mer presence matrix.
    KmerPresenceMatrix();

    // FIXME: Once we reorder the sequences and drop the mapping, the constructor
    // should take a `vector<sequence_type>&` argument and reorder it.
    /// Constructor from a given uncompressed matrix, where sequence `i` contains
    /// kmer `j` if and only if `matrix[i * num_kmers + j] == 1`.
    /// Throws `std::runtime_error` if the dimensions do not match.
    /// Moves the matrix into the object when `compress` is `false`.
    KmerPresenceMatrix(size_t num_sequences, size_t num_kmers, sdsl::bit_vector&& matrix, bool compress);

    KmerPresenceMatrix(const KmerPresenceMatrix& other) = default;
    KmerPresenceMatrix(KmerPresenceMatrix&& other) noexcept = default;
    KmerPresenceMatrix& operator=(const KmerPresenceMatrix& other) = default;
    KmerPresenceMatrix& operator=(KmerPresenceMatrix&& other) noexcept = default;

    bool operator==(const KmerPresenceMatrix& other) const;
    bool operator!=(const KmerPresenceMatrix& other) const {
        return !(*this == other);
    }

    /// Returns the size of the matrix (the number of cells in it).
    size_t size() const {
        return this->num_sequences * this->num_kmers;
    }

    /// Returns `true` if the matrix is empty.
    bool empty() const {
        return (this->size() == 0);
    }

    /// Returns the number of sequences (rows) in the matrix.
    size_t get_num_sequences() const {
        return this->num_sequences;
    }

    /// Returns the number of kmers (columns) in the matrix.
    size_t get_num_kmers() const {
        return this->num_kmers;
    }

    /// Returns the number of clusters, which is also the number of sequences
    /// stored explicitly as parents.
    size_t get_num_parents() const {
        if (this->num_sequences == 0) {
            return 0;
        }
        return this->parent_rank(this->num_sequences - 1) + 1;
    }

    /// Returns the number of present kmers in the given sequence.
    size_t get_num_present(size_t sequence) const;

    /// Returns the the total number of present kmers in the sequences.
    size_t get_total_present() const;

    /// Calls `callback(kmer_index, is_present)` for each kmer in the given sequence,
    /// ordered by their indices.
    void for_each_kmer(size_t sequence, const std::function<void(size_t, bool)>& callback) const;

    /// Scores every sequence in the subchain. The score of a sequence is the sum of
    /// `(present : 1 : -1) * kmer.second` over `kmers`. It is reported by calling
    /// `sequence_score(sequence, score)` in an arbitrary order.
    ///
    /// NOTE: Because we reuse the score of a parent when scoring its children, the
    /// scores are often different from what we would get by scoring each sequence
    /// independently. This is due to numerical instability with floating point arithmetic.
    void score_sequences(
        const std::vector<std::pair<kmer_classification, double>>& kmers,
        const std::function<void(size_t, double)>& sequence_score
    ) const;

    /// Scores a pair of sequences in the subchain. The score is the sum of
    /// `1 - abs(expected - observed)` over the kmers.
    std::int64_t diploid_score(
        size_t first, size_t second,
        const std::vector<std::pair<kmer_classification, double>>& kmers
    ) const;

    /// Serializes the object to a stream in the Simple-SDS format.
    void simple_sds_serialize(std::ostream& out) const;

    /// Loads the object from a stream in the Simple-SDS format.
    void simple_sds_load(std::istream& in);

    /// Returns the size of the object in elements.
    size_t simple_sds_size() const;
    
private:
    size_t num_sequences;
    size_t num_kmers;

    // FIXME: this can be replaced with just parent_rank
    // This is a concatenation of three vectors of length `num_sequences`.
    // 1. sequence_to_rank: sequence index to rank in stored order
    // 2. rank_to_sequence: rank in stored order to sequence index
    // 3. parent_rank: overall rank to parent rank
    sdsl::int_vector<> arrays;

    // Concatenated parent bitvectors for each cluster.
    sdsl::bit_vector parents;

    // Concatenated child bitvectors for each cluster.
    sdsl::sd_vector<> children;

    // Transforms the external sequence index to the internal rank in the stored order.
    size_t sequence_to_rank(size_t sequence) const {
        return this->arrays[sequence];
    }

    // Transforms the internal rank in the stored order to the external sequence index.
    size_t rank_to_sequence(size_t rank) const {
        return this->arrays[this->num_sequences + rank];
    }

    // Transforms the overall rank into the rank of the parent in `parents`.
    size_t parent_rank(size_t rank) const {
        return this->arrays[2 * this->num_sequences + rank];
    }

    // Transforms the overall rank into the rank in `children`, assuming that this is a child.
    size_t child_rank(size_t rank) const {
        return rank - this->parent_rank(rank) - 1;
    }

    // The first sequence in a cluster is the parent that is stored explicitly.
    bool is_parent(size_t rank) const {
        if(rank == 0) {
            return true;
        }
        return (this->parent_rank(rank) != this->parent_rank(rank - 1));
    }

    // An iterator over all kmers in a sequence.
    // If the sequence is a parent, `child_iter` is at the end and never matches.
    struct Iterator {
        const KmerPresenceMatrix& matrix;
        size_t parent_start;
        size_t child_start;
        sdsl::sd_vector<>::one_iterator child_iter;
        size_t kmer_id;
        bool is_present;

        Iterator(const KmerPresenceMatrix& matrix, size_t sequence);

        void operator++() {
            this->kmer_id++;
            this->update_is_present();
        }

        bool end() const {
            return (this->kmer_id >= this->matrix.num_kmers);
        }

    private:
        void update_is_present();
    };

    // An iterator over flipped bits in a child sequence.
    struct ChildIterator {
        const KmerPresenceMatrix& matrix;
        size_t start;
        sdsl::sd_vector<>::one_iterator iter;

        ChildIterator(const KmerPresenceMatrix& matrix, size_t rank);

        void operator++() {
            ++this->iter;
        }

        size_t kmer_id() const {
            return this->iter->second - this->start;
        }

        bool end() const {
            return (this->iter == this->matrix.children.one_end() || this->iter->second >= this->start + this->matrix.num_kmers);
        }
    };
};

//------------------------------------------------------------------------------

// FIXME: version 7: each top-level chain should know the default scoring model for that chain

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
 * * Version 7: Work in progress.
 *   Compressed kmer presence matrices, with faster haplotype scoring by
 *   taking advantage of the compression.
 *   Versions 5 and 6 can still be read.
 *
 * * Version 6: Tags for storing metadata, such as graph names. Version 5 can
 *   still be read.
 *
 * * Version 5: Every path in the graph is assigned to a construction job.
 *   This allows including reference paths that do not visit any snarls in the
 *   sampled graph. Earlier versions cannot be used anymore.
 *
 * * Version 4: Subchains can have fragmented haplotypes instead of a single
 *   GBWT sequence always crossing from start to end. Versions 1 to 3 can still
 *   be read.
 *
 * * Version 3: Subchains use smaller integers when possible. Versions 1 and 2
 *   can still be read.
 *
 * * Version 2: Top-level chains include a contig name. Versions 1 can still be
 *   read.
 *
 * * Version 1: Initial version.
 */
class Haplotypes {
public:
    /// Default target length of a subchain.
    constexpr static size_t SUBCHAIN_LENGTH = 10000;

    /// Block size (in kmers) for reading KFF files.
    constexpr static size_t KFF_BLOCK_SIZE = 1000000;

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
        constexpr static std::uint32_t VERSION = 7;
        constexpr static std::uint32_t COMPRESSED_VERSION = 7;
        constexpr static std::uint32_t VERSION_WITH_TAGS = 6;
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

    // FIXME: add a placeholder word for scoring model override
    /// Representation of a subchain.
    /// The subchain contains `num_sequences()` haplotype sequences and `num_kmers()` kmers.
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
        /// of sequences it appears in.
        std::vector<kmer_type> kmers;

        /// Number of sequences each kmer appears in.
        sdsl::int_vector<0> kmer_counts;

        /// Sequences as (GBWT sequence id, offset in the relevant node).
        std::vector<compact_sequence_type> sequences;

        /// A matrix indicating which kmers are present in which sequences.
        KmerPresenceMatrix kmers_present;

        /// Calls `callback(kmer_id, is_present)` for each kmer in the given sequence,
        /// ordered by their indices.
        void for_each_kmer(size_t sequence, const std::function<void(size_t, bool)>& callback) const {
            this->kmers_present.for_each_kmer(sequence, callback);
        }

        /// Scores every sequence in the subchain. The score of a sequence is the sum of
        /// `(present : 1 : -1) * kmer.second` over `kmers`. It is reported by calling
        /// `sequence_score(sequence, score)` in an arbitrary order.
        void score_sequences(
            const std::vector<std::pair<KmerPresenceMatrix::kmer_classification, double>>& kmers,
            const std::function<void(size_t, double)>& haplotype_score
        ) const {
            this->kmers_present.score_sequences(kmers, haplotype_score);
        }

        /// Scores a pair of sequences in the subchain. The score is the sum of
        /// `1 - abs(expected - observed)` over the kmers.
        std::int64_t diploid_score(
            size_t first, size_t second,
            const std::vector<std::pair<KmerPresenceMatrix::kmer_classification, double>>& kmers
        ) const {
            return this->kmers_present.diploid_score(first, second, kmers);
        }

        /// Returns the number of haplotype sequences in the subchain.
        size_t num_sequences() const { return this->sequences.size(); }

        /// Returns the number of kmers in the subchain.
        size_t num_kmers() const { return this->kmers.size(); }

        /// Returns the number of kmers present in the given sequence.
        size_t num_present(size_t i) const {
            return this->kmers_present.get_num_present(i);
        }

        /// Returns the size of the kmer presence matrix.
        size_t total_kmers() const { return this->kmers_present.size(); }

        /// Returns the total number of kmers present in the sequences.
        size_t total_present() const {
            return this->kmers_present.get_total_present();
        }

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
        /// As there is no header, the expected version must be provided.
        void simple_sds_load(std::istream& in, std::uint32_t version = Header::VERSION);

        /// Returns the size of the object in elements.
        size_t simple_sds_size() const;
    };

    // FIXME: add a placeholder word for default scoring model
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
        /// As there is no header, the expected version must be provided.
        void simple_sds_load(std::istream& in, std::uint32_t version = Header::VERSION);

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

    // Tags for storing metadata, such as graph names.
    gbwt::Tags tags;

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

    /**
     * Sets graph name (pggname) information based on the given GBZ graph.
     */
    void set_graph_name(const gbwtgraph::GBZ& gbz);

    /**
     * Returns the graph name (pggname) information stored in the tags.
     */
    gbwtgraph::GraphName graph_name() const { return gbwtgraph::GraphName(this->tags); }
};

//------------------------------------------------------------------------------

} // namespace vg

#endif // VG_RECOMBINATOR_HAPLOTYPES_HPP_INCLUDED
