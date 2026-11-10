#include "dedup_internal.hpp"

#include "parallel_sort.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <tuple>
#include <unordered_map>

#include <omp.h>

namespace vg {

using handlegraph::handle_t;
using handlegraph::HandleGraph;
using handlegraph::nid_t;

//------------------------------------------------------------------------------

void ContentChopParameters::validate() const {
    if (max_node_length == 0) {
        throw std::invalid_argument("ContentChopParameters: the maximum node length must be positive");
    }
    if (!std::isfinite(piece_weight) || piece_weight < 0.0) {
        throw std::invalid_argument("ContentChopParameters: the piece weight must be finite and non-negative");
    }
    if (kmer_length == 0 || kmer_length > 32) {
        throw std::invalid_argument("ContentChopParameters: the k-mer length must be between 1 and 32");
    }
}

//------------------------------------------------------------------------------

namespace {

/// Return the 2-bit code of a base, or 4 for a character that is not a base.
uint64_t base_code(char c) {
    switch (c) {
    case 'A': case 'a': return 0;
    case 'C': case 'c': return 1;
    case 'G': case 'g': return 2;
    case 'T': case 't': return 3;
    default: return 4;
    }
}

/// Return the hash of a packed k-mer. The MurmurHash3 finalizer is a bijection, so distinct
/// k-mers have distinct hashes.
uint64_t kmer_hash(uint64_t kmer) {
    kmer ^= kmer >> 33;
    kmer *= 0xFF51AFD7ED558CCDULL;
    kmer ^= kmer >> 33;
    kmer *= 0xC4CEB9FE1A85EC53ULL;
    kmer ^= kmer >> 33;
    return kmer;
}

}

std::vector<KmerAnchor> find_anchors(std::string_view sequence, size_t kmer_length, size_t radius) {
    // Collect the k-mers made of bases, with their hashes.
    struct HashedKmer {
        size_t offset;
        uint64_t kmer;
        uint64_t hash;
    };
    std::vector<HashedKmer> kmers;
    const uint64_t mask = (kmer_length >= 32 ? ~uint64_t(0) : (uint64_t(1) << (2 * kmer_length)) - 1);
    uint64_t packed = 0;
    size_t bases_in_a_row = 0;
    for (size_t i = 0; i < sequence.size(); i++) {
        uint64_t code = base_code(sequence[i]);
        if (code > 3) {
            bases_in_a_row = 0;
            continue;
        }
        packed = ((packed << 2) | code) & mask;
        bases_in_a_row++;
        if (bases_in_a_row >= kmer_length) {
            kmers.push_back({ i + 1 - kmer_length, packed, kmer_hash(packed) });
        }
    }

    // A k-mer is an anchor when the nearest k-mer with a hash at most its own is more than
    // `radius` bases away on both sides. We find the nearest such k-mer on each side with a stack
    // of k-mers whose hashes increase from the bottom. Before a k-mer goes on the stack, the
    // k-mers above it with larger hashes come off, and the k-mer left on top is the nearest one
    // with a hash at most its own.
    std::vector<bool> clear_on_left(kmers.size(), false);
    std::vector<size_t> stack;
    for (size_t i = 0; i < kmers.size(); i++) {
        while (!stack.empty() && kmers[stack.back()].hash > kmers[i].hash) {
            stack.pop_back();
        }
        clear_on_left[i] = (stack.empty() || kmers[i].offset - kmers[stack.back()].offset > radius);
        stack.push_back(i);
    }
    std::vector<bool> clear_on_right(kmers.size(), false);
    stack.clear();
    for (size_t i = kmers.size(); i-- > 0;) {
        while (!stack.empty() && kmers[stack.back()].hash > kmers[i].hash) {
            stack.pop_back();
        }
        clear_on_right[i] = (stack.empty() || kmers[stack.back()].offset - kmers[i].offset > radius);
        stack.push_back(i);
    }

    std::vector<KmerAnchor> anchors;
    for (size_t i = 0; i < kmers.size(); i++) {
        if (clear_on_left[i] && clear_on_right[i]) {
            anchors.push_back({ kmers[i].offset, kmers[i].kmer });
        }
    }
    return anchors;
}

//------------------------------------------------------------------------------

namespace {

/// The prime modulus of the fingerprints, 2^61 - 1.
constexpr uint64_t FINGERPRINT_PRIME = (uint64_t(1) << 61) - 1;

/// The base of the fingerprint polynomial.
constexpr uint64_t FINGERPRINT_BASE = 0x0123456789ABCDEFULL;

/// Return `x` modulo the prime, for `x < 2^64`. Since `2^61` is 1 modulo the prime, the bits
/// above the lowest 61 can be added to them.
uint64_t reduce_fingerprint(uint64_t x) {
    x = (x & FINGERPRINT_PRIME) + (x >> 61);
    return (x >= FINGERPRINT_PRIME ? x - FINGERPRINT_PRIME : x);
}

/// Return `a * b` modulo the prime, for `a` and `b` less than the prime.
uint64_t multiply_fingerprints(uint64_t a, uint64_t b) {
    __uint128_t product = static_cast<__uint128_t>(a) * b;
    uint64_t low = static_cast<uint64_t>(product) & FINGERPRINT_PRIME;
    uint64_t high = static_cast<uint64_t>(product >> 61);
    return reduce_fingerprint(low + high);
}

/// Return the fingerprint of the string `s + c`, given the fingerprint of `s`.
uint64_t extend_fingerprint(uint64_t fingerprint, char c) {
    return reduce_fingerprint(multiply_fingerprints(fingerprint, FINGERPRINT_BASE) + static_cast<unsigned char>(c));
}

}

uint64_t sequence_fingerprint(std::string_view sequence) {
    uint64_t result = 0;
    for (char c : sequence) {
        result = extend_fingerprint(result, c);
    }
    return result;
}

SubstringFingerprinter::SubstringFingerprinter(std::string_view sequence, size_t max_length) :
    prefix_fingerprints(sequence.size() + 1, 0), base_powers(std::min(max_length, sequence.size()) + 1, 1) {
    for (size_t i = 0; i < sequence.size(); i++) {
        prefix_fingerprints[i + 1] = extend_fingerprint(prefix_fingerprints[i], sequence[i]);
    }
    for (size_t i = 1; i < base_powers.size(); i++) {
        base_powers[i] = multiply_fingerprints(base_powers[i - 1], FINGERPRINT_BASE);
    }
}

uint64_t SubstringFingerprinter::fingerprint(size_t begin, size_t end) const {
    // The prefix of length `end` is the prefix of length `begin` times `B^(end - begin)`, plus
    // the substring.
    uint64_t shifted = multiply_fingerprints(prefix_fingerprints[begin], base_powers[end - begin]);
    return reduce_fingerprint(prefix_fingerprints[end] + FINGERPRINT_PRIME - shifted);
}

//------------------------------------------------------------------------------

CutCandidateSet::CutCandidateSet(std::string_view sequence, const ContentChopParameters& parameters) :
    max_length(parameters.max_node_length), fingerprinter(sequence, parameters.max_node_length) {
    // Add a candidate after the last one, with fillers before it where the gap would be longer
    // than a piece.
    auto add_candidate = [&](size_t offset) {
        if (offset == cut_offsets.back()) {
            // The candidate is already there: an anchor at offset 0, or the end of an empty
            // sequence.
            return;
        }
        while (offset - cut_offsets.back() > max_length) {
            cut_offsets.push_back(cut_offsets.back() + max_length);
        }
        cut_offsets.push_back(offset);
    };

    cut_offsets.push_back(0);
    for (const KmerAnchor& anchor : find_anchors(sequence, parameters.kmer_length, parameters.anchor_radius)) {
        add_candidate(anchor.offset);
    }
    add_candidate(sequence.size());
}

//------------------------------------------------------------------------------

SegmentCountTable::SegmentCountTable(std::vector<uint64_t> fingerprints) {
    parallel_sort(fingerprints.begin(), fingerprints.end(), std::less<uint64_t>());
    for (uint64_t fingerprint : fingerprints) {
        if (entries.empty() || entries.back().fingerprint != fingerprint) {
            entries.push_back({ fingerprint, 0 });
        }
        entries.back().count++;
    }
}

size_t SegmentCountTable::count(uint64_t fingerprint) const {
    auto iter = std::lower_bound(entries.begin(), entries.end(), fingerprint, [](const Entry& entry, uint64_t value) {
        return entry.fingerprint < value;
    });
    return (iter != entries.end() && iter->fingerprint == fingerprint ? iter->count : 0);
}

//------------------------------------------------------------------------------

std::vector<size_t> choose_cuts(const CutCandidateSet& candidates, const SegmentCountTable& counts,
                                double piece_weight) {
    const std::vector<size_t>& offsets = candidates.offsets();
    size_t max_length = candidates.max_segment_length();

    // The least cost of chopping the sequence up to each candidate, and the candidate where the
    // last piece of that chopping starts. Consecutive candidates are at most a segment apart, so
    // each candidate after the first can be reached.
    std::vector<double> best_cost(offsets.size(), std::numeric_limits<double>::infinity());
    std::vector<size_t> piece_start(offsets.size(), 0);
    best_cost[0] = 0.0;
    size_t first_in_reach = 0;
    for (size_t to = 1; to < offsets.size(); to++) {
        while (offsets[to] - offsets[first_in_reach] > max_length) {
            first_in_reach++;
        }
        // Trying the starts from the left, and replacing the best only with a strictly better
        // cost, prefers the longest last piece among equal costs.
        for (size_t from = first_in_reach; from < to; from++) {
            size_t occurrences = std::max<size_t>(counts.count(candidates.segment_fingerprint(from, to)), 1);
            double length = static_cast<double>(offsets[to] - offsets[from]);
            double cost = best_cost[from] + piece_weight + length / static_cast<double>(occurrences);
            if (cost < best_cost[to]) {
                best_cost[to] = cost;
                piece_start[to] = from;
            }
        }
    }

    std::vector<size_t> cuts;
    for (size_t at = offsets.size() - 1; at > 0 && piece_start[at] > 0; at = piece_start[at]) {
        cuts.push_back(offsets[piece_start[at]]);
    }
    std::reverse(cuts.begin(), cuts.end());
    return cuts;
}

//------------------------------------------------------------------------------

const PieceLocation* ChoppingPlan::find_piece(nid_t node, size_t node_offset) const {
    auto iter = std::lower_bound(pieces.begin(), pieces.end(), std::make_pair(node, node_offset),
                                 [](const PieceLocation& piece, const std::pair<nid_t, size_t>& key) {
        return std::make_pair(piece.node, piece.node_offset) < key;
    });
    if (iter == pieces.end() || iter->node != node || iter->node_offset != node_offset) {
        return nullptr;
    }
    return &*iter;
}

std::string_view ChoppingPlan::piece_sequence(const PieceLocation& piece) const {
    return std::string_view(strings[piece.string_id]).substr(piece.string_offset, piece.length);
}

size_t ChoppingPlan::stored_length() const {
    size_t total = 0;
    for (const std::string& s : strings) {
        total += s.size();
    }
    return total;
}

//------------------------------------------------------------------------------

namespace {

/// A place in the forward sequence of a node.
struct ForwardPosition {
    nid_t node;
    size_t offset;
};

/**
 * A PiecePacker is a packer that puts the sequences of the pieces of a graph into as little
 * stored text as it can. It decides, for each piece, whether its sequence is already stored, and
 * stores it otherwise.
 *
 * The stored text is the set of stored pieces. Stored pieces that are consecutive in a node form
 * one run of text, so a piece can be found across the boundary between two of them. The source
 * of each placed piece is a position where the stored text holds its sequence. A piece that is
 * stored is its own source.
 */
class PiecePacker {
public:
    /// Prepare to pack the pieces, given in order of node and offset. They must be non-empty, and
    /// they must cover each of their nodes with no gaps or overlaps.
    PiecePacker(const HandleGraph& graph, std::vector<PieceLocation> pieces, const ContentChopParameters& parameters);

    /// Place every piece, longest first, and fill the strings and the piece locations of the
    /// plan.
    void pack(ChoppingPlan& plan);

private:
    /// A piece and the progress of placing it.
    struct PackedPiece {
        nid_t node;
        size_t offset;
        size_t length;
        /// True if the sequence of the piece is part of the stored text.
        bool stored = false;
        /// The source of the piece, once it is placed.
        ForwardPosition source = { 0, 0 };
    };

    /// Find a source for the piece, and store the piece if its sequence is not stored yet.
    void place(size_t piece_id);

    /// Return true and set the source of the piece if an earlier piece has the same sequence.
    bool find_equal_piece(PackedPiece& piece, const std::string& sequence, uint64_t fingerprint) const;

    /// Return true and set the source of the piece if the stored text contains its sequence.
    bool find_in_stored_text(PackedPiece& piece, const std::string& sequence) const;

    /// Return true if stored pieces cover the interval [begin, end) of the node.
    bool is_stored(nid_t node, size_t begin, size_t end) const;

    /// Make the piece part of the stored text, and index the anchors in it.
    void store(PackedPiece& piece);

    /// Return the index of the piece of the node that contains the offset, or the number of
    /// pieces if there is none.
    size_t piece_containing(nid_t node, size_t offset) const;

    /// Return the sequence of the interval of a node.
    std::string node_sequence(nid_t node, size_t offset, size_t length) const;

    /// The most places with a matching anchor that we try for a piece. A k-mer of a satellite
    /// repeat can occur at very many places, and comparing the piece with all of them would take
    /// too long.
    static constexpr size_t MAX_TRIED_PLACES = 32;

    const HandleGraph& graph;
    size_t kmer_length;
    size_t anchor_radius;

    /// The pieces, in order of node and offset.
    std::vector<PackedPiece> pieces;

    /// The placed pieces with each fingerprint.
    std::unordered_map<uint64_t, std::vector<size_t>> pieces_by_fingerprint;

    /// The positions of the anchors in the stored text with each k-mer.
    std::unordered_map<uint64_t, std::vector<ForwardPosition>> stored_anchors;
};

PiecePacker::PiecePacker(const HandleGraph& graph, std::vector<PieceLocation> locations,
                         const ContentChopParameters& parameters) :
    graph(graph), kmer_length(parameters.kmer_length), anchor_radius(parameters.anchor_radius) {
    pieces.reserve(locations.size());
    for (const PieceLocation& location : locations) {
        pieces.push_back({ location.node, location.node_offset, location.length });
    }
}

void PiecePacker::pack(ChoppingPlan& plan) {
    // Long pieces go first. A short piece placed later can then be found inside them.
    std::vector<size_t> order(pieces.size());
    for (size_t i = 0; i < order.size(); i++) {
        order[i] = i;
    }
    std::sort(order.begin(), order.end(), [&](size_t a, size_t b) {
        return std::make_tuple(pieces[b].length, pieces[a].node, pieces[a].offset)
            < std::make_tuple(pieces[a].length, pieces[b].node, pieces[b].offset);
    });
    for (size_t piece_id : order) {
        place(piece_id);
    }

    // Each run of consecutive stored pieces in a node becomes a string. The runs are in order of
    // node and offset.
    struct StoredRun {
        nid_t node;
        size_t begin;
        size_t end;
    };
    std::vector<StoredRun> runs;
    for (const PackedPiece& piece : pieces) {
        if (!piece.stored) {
            continue;
        }
        if (!runs.empty() && runs.back().node == piece.node && runs.back().end == piece.offset) {
            runs.back().end += piece.length;
        } else {
            runs.push_back({ piece.node, piece.offset, piece.offset + piece.length });
        }
    }
    plan.strings.clear();
    plan.strings.reserve(runs.size());
    for (const StoredRun& run : runs) {
        plan.strings.push_back(node_sequence(run.node, run.begin, run.end - run.begin));
    }

    plan.pieces.clear();
    plan.pieces.reserve(pieces.size());
    for (const PackedPiece& piece : pieces) {
        // The run that contains the source is the last one starting at or before it.
        auto iter = std::upper_bound(runs.begin(), runs.end(), piece.source, [](const ForwardPosition& pos, const StoredRun& run) {
            return std::make_pair(pos.node, pos.offset) < std::make_pair(run.node, run.begin);
        });
        --iter;
        size_t string_id = iter - runs.begin();
        plan.pieces.push_back({ piece.node, piece.offset, string_id, piece.source.offset - iter->begin, piece.length });
    }
}

void PiecePacker::place(size_t piece_id) {
    PackedPiece& piece = pieces[piece_id];
    std::string sequence = node_sequence(piece.node, piece.offset, piece.length);
    uint64_t fingerprint = sequence_fingerprint(sequence);
    if (!find_equal_piece(piece, sequence, fingerprint) && !find_in_stored_text(piece, sequence)) {
        store(piece);
    }
    pieces_by_fingerprint[fingerprint].push_back(piece_id);
}

bool PiecePacker::find_equal_piece(PackedPiece& piece, const std::string& sequence, uint64_t fingerprint) const {
    auto iter = pieces_by_fingerprint.find(fingerprint);
    if (iter == pieces_by_fingerprint.end()) {
        return false;
    }
    for (size_t other_id : iter->second) {
        const PackedPiece& other = pieces[other_id];
        // Equal fingerprints almost always mean equal sequences, but we make sure.
        if (other.length == piece.length && node_sequence(other.node, other.offset, other.length) == sequence) {
            piece.source = other.source;
            return true;
        }
    }
    return false;
}

bool PiecePacker::find_in_stored_text(PackedPiece& piece, const std::string& sequence) const {
    // An anchor of the piece at least `anchor_radius` bases from both ends of the piece is an
    // anchor wherever the piece occurs. Any occurrence of the piece in the stored text then has
    // that anchor at the same offset, and the anchor is indexed. If some such anchor of the piece
    // is not in the stored text, neither is the piece.
    const std::vector<ForwardPosition>* rarest_places = nullptr;
    size_t rarest_offset = 0;
    for (const KmerAnchor& anchor : find_anchors(sequence, kmer_length, anchor_radius)) {
        if (anchor.offset < anchor_radius || anchor.offset + anchor_radius + kmer_length > sequence.size()) {
            continue;
        }
        auto iter = stored_anchors.find(anchor.kmer);
        if (iter == stored_anchors.end()) {
            return false;
        }
        if (rarest_places == nullptr || iter->second.size() < rarest_places->size()) {
            rarest_places = &iter->second;
            rarest_offset = anchor.offset;
        }
    }
    if (rarest_places == nullptr) {
        // The piece is too short or too repetitive to have anchors far enough from its ends.
        return false;
    }

    size_t places_to_try = std::min(rarest_places->size(), MAX_TRIED_PLACES);
    for (size_t i = 0; i < places_to_try; i++) {
        const ForwardPosition& place = (*rarest_places)[i];
        if (place.offset < rarest_offset) {
            continue;
        }
        size_t begin = place.offset - rarest_offset;
        if (!is_stored(place.node, begin, begin + sequence.size())) {
            continue;
        }
        if (node_sequence(place.node, begin, sequence.size()) == sequence) {
            piece.source = { place.node, begin };
            return true;
        }
    }
    return false;
}

bool PiecePacker::is_stored(nid_t node, size_t begin, size_t end) const {
    size_t i = piece_containing(node, begin);
    while (begin < end) {
        if (i >= pieces.size() || pieces[i].node != node || !pieces[i].stored) {
            return false;
        }
        begin = pieces[i].offset + pieces[i].length;
        i++;
    }
    return true;
}

void PiecePacker::store(PackedPiece& piece) {
    piece.stored = true;
    piece.source = { piece.node, piece.offset };

    // Whether a k-mer is an anchor of the node depends on the bases up to `anchor_radius` bases
    // before it and after it. We find the anchors in the piece and that much of the node on each
    // side.
    size_t node_length = graph.get_length(graph.get_handle(piece.node));
    size_t window_begin = (piece.offset > anchor_radius ? piece.offset - anchor_radius : 0);
    size_t window_end = std::min(piece.offset + piece.length + anchor_radius + kmer_length - 1, node_length);
    std::string window = node_sequence(piece.node, window_begin, window_end - window_begin);
    for (const KmerAnchor& anchor : find_anchors(window, kmer_length, anchor_radius)) {
        size_t offset = window_begin + anchor.offset;
        if (offset >= piece.offset && offset < piece.offset + piece.length) {
            stored_anchors[anchor.kmer].push_back({ piece.node, offset });
        }
    }
}

size_t PiecePacker::piece_containing(nid_t node, size_t offset) const {
    // The piece containing the offset is the last one starting at or before it.
    auto iter = std::upper_bound(pieces.begin(), pieces.end(), std::make_pair(node, offset),
                                 [](const std::pair<nid_t, size_t>& key, const PackedPiece& piece) {
        return key < std::make_pair(piece.node, piece.offset);
    });
    if (iter == pieces.begin()) {
        return pieces.size();
    }
    --iter;
    if (iter->node != node || offset >= iter->offset + iter->length) {
        return pieces.size();
    }
    return iter - pieces.begin();
}

std::string PiecePacker::node_sequence(nid_t node, size_t offset, size_t length) const {
    return graph.get_subsequence(graph.get_handle(node), offset, length);
}

}

//------------------------------------------------------------------------------

ChoppingPlan plan_content_chopping(const HandleGraph& graph, const ContentChopParameters& parameters) {
    parameters.validate();

    // The nodes to consider for chopping, in order of ID.
    std::vector<nid_t> long_nodes;
    graph.for_each_handle([&](const handle_t& handle) {
        if (graph.get_length(handle) > parameters.max_node_length) {
            long_nodes.push_back(graph.get_id(handle));
        }
    });
    std::sort(long_nodes.begin(), long_nodes.end());

    // Count the segments of the long nodes.
    std::vector<uint64_t> fingerprints;
    #pragma omp parallel
    {
        std::vector<uint64_t> thread_fingerprints;
        #pragma omp for schedule(dynamic, 1)
        for (size_t i = 0; i < long_nodes.size(); i++) {
            std::string sequence = graph.get_sequence(graph.get_handle(long_nodes[i]));
            CutCandidateSet candidates(sequence, parameters);
            candidates.for_each_segment([&](size_t from, size_t to) {
                thread_fingerprints.push_back(candidates.segment_fingerprint(from, to));
            });
        }
        #pragma omp critical
        {
            fingerprints.insert(fingerprints.end(), thread_fingerprints.begin(), thread_fingerprints.end());
        }
    }
    SegmentCountTable counts(std::move(fingerprints));

    // Choose the cuts in the long nodes. Keeping the candidate sets from counting would take
    // memory proportional to the total length of the long nodes, so we find them again.
    std::vector<std::vector<size_t>> cuts(long_nodes.size());
    #pragma omp parallel for schedule(dynamic, 1)
    for (size_t i = 0; i < long_nodes.size(); i++) {
        std::string sequence = graph.get_sequence(graph.get_handle(long_nodes[i]));
        CutCandidateSet candidates(sequence, parameters);
        cuts[i] = choose_cuts(candidates, counts, parameters.piece_weight);
    }

    ChoppingPlan plan;
    for (size_t i = 0; i < long_nodes.size(); i++) {
        if (!cuts[i].empty()) {
            plan.chop_points.emplace_back(long_nodes[i], cuts[i]);
        }
    }

    // List the pieces of the chopped nodes, and pack their sequences. The cuts are strictly
    // inside each node, so the pieces are non-empty.
    std::vector<PieceLocation> pieces;
    for (const auto& [node, node_cuts] : plan.chop_points) {
        size_t length = graph.get_length(graph.get_handle(node));
        size_t start = 0;
        for (size_t cut : node_cuts) {
            pieces.push_back({ node, start, 0, 0, cut - start });
            start = cut;
        }
        pieces.push_back({ node, start, 0, 0, length - start });
    }
    PiecePacker packer(graph, std::move(pieces), parameters);
    packer.pack(plan);

    return plan;
}

//------------------------------------------------------------------------------

}
