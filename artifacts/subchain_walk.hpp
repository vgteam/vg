#ifndef VG_SUBCHAIN_WALK_HPP_INCLUDED
#define VG_SUBCHAIN_WALK_HPP_INCLUDED

/** \file
 * The walk of a haplotype through a subchain of the haplotype information, and the counts
 * of the kmers of the subchain in such walks.
 */

#include "../src/recombinator_haplotypes.hpp"

#include <functional>
#include <limits>
#include <string>
#include <tuple>
#include <vector>

namespace vg {

/**
 * The walk of a GBWT sequence through a subchain. It starts at node `start`, at the
 * offset of the sequence there, or at the start of the sequence without a start node
 * (`gbwt::ENDMARKER`). It ends at node `end`, or at the end of the sequence without an
 * end node. With both nodes, where a fragment of the haplotype ends, the walk continues
 * with the next fragment, because the generated haplotypes do the same. The part of the
 * walk in one fragment is a piece.
 */
class SubchainWalk {
public:
    /// An empty walk.
    SubchainWalk() = default;

    /// Walk the sequence through the subchain from `start` to `end`.
    SubchainWalk(
        const gbwt::CachedGBWT& index, const gbwt::FragmentMap& fragment_map,
        gbwt::node_type start, gbwt::node_type end, Haplotypes::sequence_type sequence
    );

    /// Call `visit(node, starts_piece)` for each node the walk of the sequence from `start`
    /// to `end` visits, in order, where `starts_piece` is true for the first node of each
    /// piece. Stop when `visit` returns false. Return false if `visit` stopped the walk.
    static bool for_each_visit(
        const gbwt::CachedGBWT& index, const gbwt::FragmentMap& fragment_map,
        gbwt::node_type start, gbwt::node_type end, Haplotypes::sequence_type sequence,
        const std::function<bool(gbwt::node_type, bool)>& visit
    );

    /// Return the ids of the nodes the walk visits, sorted and without duplicates.
    std::vector<nid_t> node_ids() const;

    /// Return the sequence of each piece. The first node of the walk contributes its
    /// last `start_max` bases and the end node its first `end_max` bases. A walk whose
    /// first node is the end node has no sequence, and neither has an empty walk.
    std::vector<std::string> sequence_pieces(const gbwtgraph::GBWTGraph& graph, size_t start_max, size_t end_max) const;

    /// Return the number of pieces in `sequence_pieces()`.
    size_t pieces() const;

    /// Two walks are equal if they visit the same nodes in the same pieces.
    bool operator==(const SubchainWalk& another) const = default;

    /// An order of the walks, for grouping equal ones.
    bool operator<(const SubchainWalk& another) const {
        return std::tie(this->visits, this->piece_starts, this->reaches_end) < std::tie(another.visits, another.piece_starts, another.reaches_end);
    }

private:
    // The nodes in the order the walk visits them, and the offset in `visits` where each
    // piece starts.
    std::vector<gbwt::node_type> visits;
    std::vector<std::uint32_t> piece_starts;

    // The walk ends at the end node.
    bool reaches_end = false;

    // Return true if `sequence_pieces()` returns any.
    bool has_sequence() const;
};

/**
 * A table from the kmers of a list, in both orientations, to their positions in the
 * list. Where a kmer is also the reverse complement of another kmer of the list, the
 * later position in the list wins. The table is an open-addressing hash table at most
 * half full, because it is looked up at every base of the sequences.
 */
class OrientedKmerTable {
public:
    /// The position returned for a kmer that is not in the table.
    constexpr static std::uint32_t NONE = std::numeric_limits<std::uint32_t>::max();

    /// Build the table for the kmers of length `k`.
    OrientedKmerTable(const std::vector<Haplotypes::Subchain::kmer_type>& kmers, size_t k);

    /// Return the position of the kmer in the list, or `NONE`.
    std::uint32_t find(Haplotypes::Subchain::kmer_type kmer) const {
        for (size_t i = this->slot(kmer); this->keys[i] != EMPTY; i = (i + 1) & this->mask) {
            if (this->keys[i] == kmer) {
                return this->positions[i];
            }
        }
        return NONE;
    }

private:
    // Kmers are at most 31 bases, so a key with every bit set is not a kmer.
    constexpr static Haplotypes::Subchain::kmer_type EMPTY = std::numeric_limits<Haplotypes::Subchain::kmer_type>::max();

    std::vector<Haplotypes::Subchain::kmer_type> keys;
    std::vector<std::uint32_t> positions;
    size_t mask = 0;
    size_t shift = 0;

    size_t slot(Haplotypes::Subchain::kmer_type kmer) const {
        return static_cast<size_t>((kmer * 0x9E3779B97F4A7C15ULL) >> this->shift) & this->mask;
    }
};

/// The number of occurrences of the kmers of an `OrientedKmerTable`, in either
/// orientation, in the pieces of one sequence, by the position of the kmer in the list
/// the table was built from. A kmer is counted where a piece has `k` valid characters
/// ending at it. Counting a new sequence replaces the counts of the previous one, and
/// takes time proportional to the length of the sequence, not to the number of kmers.
class OrientedKmerCounts {
public:
    /// Prepare to count the kmers of a table built from a list of `kmers` kmers.
    explicit OrientedKmerCounts(size_t kmers) : counts(kmers, 0) {}

    /// Count the kmers of the table in the pieces of a sequence.
    void count(const std::vector<std::string>& pieces, const OrientedKmerTable& positions, size_t k);

    /// Return (position, occurrences) for the kmers that occur in the sequence, in the
    /// order of their first occurrence.
    std::vector<std::pair<std::uint32_t, std::uint32_t>> occurrence_list() const;

    /// The number of occurrences of the kmer at the given position.
    std::uint32_t operator[](size_t position) const { return this->counts[position]; }

private:
    std::vector<std::uint32_t> counts;
    std::vector<std::uint32_t> occurring_positions;
};

} // namespace vg

#endif // VG_SUBCHAIN_WALK_HPP_INCLUDED
