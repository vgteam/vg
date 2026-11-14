#include "dedup_internal.hpp"

#include "../artifacts/parallel_sort.hpp"

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
        throw std::invalid_argument("ContentChopParameters: maximum node length must be positive");
    }
    if (!std::isfinite(piece_weight) || piece_weight < 0.0) {
        throw std::invalid_argument("ContentChopParameters: piece weight must be finite and non-negative");
    }
    if (kmer_length == 0 || kmer_length > 32) {
        throw std::invalid_argument("ContentChopParameters: k-mer length must be between 1 and 32");
    }
}

//------------------------------------------------------------------------------

namespace {

std::uint64_t base_code(char c) {
    switch (c) {
    case 'A': case 'a': return 0;
    case 'C': case 'c': return 1;
    case 'G': case 'g': return 2;
    case 'T': case 't': return 3;
    default: return 4;
    }
}

// MurmurHash3 finalizer is a bijection, so distinct k-mers have distinct hashes.
std::uint64_t kmer_hash(std::uint64_t kmer) {
    kmer ^= kmer >> 33;
    kmer *= 0xFF51AFD7ED558CCDULL;
    kmer ^= kmer >> 33;
    kmer *= 0xC4CEB9FE1A85EC53ULL;
    kmer ^= kmer >> 33;
    return kmer;
}

}

std::vector<KmerAnchor> find_anchors(std::string_view sequence, size_t kmer_length, size_t radius) {
    // Collect valid forward strand k-mers.
    struct HashedKmer {
        size_t offset;
        std::uint64_t kmer;
        std::uint64_t hash;
    };
    std::vector<HashedKmer> kmers;
    const std::uint64_t mask = (kmer_length >= 32 ? ~std::uint64_t(0) : (std::uint64_t(1) << (2 * kmer_length)) - 1);
    std::uint64_t packed = 0;
    size_t bases_in_a_row = 0;
    for (size_t i = 0; i < sequence.size(); i++) {
        std::uint64_t code = base_code(sequence[i]);
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

    // Maintain a stack of k-mers with increasing hash values, up to the current k-mer.
    // Use it to select left-minimal k-mers (within the radius).
    std::vector<bool> minimal(kmers.size(), false);
    std::vector<size_t> stack;
    for (size_t i = 0; i < kmers.size(); i++) {
        while (!stack.empty() && kmers[stack.back()].hash > kmers[i].hash) {
            stack.pop_back();
        }
        minimal[i] = (stack.empty() || kmers[i].offset - kmers[stack.back()].offset > radius);
        stack.push_back(i);
    }

    // And now select right-minimal k-mers (within the radius).
    stack.clear();
    for (size_t i = kmers.size(); i-- > 0;) {
        while (!stack.empty() && kmers[stack.back()].hash > kmers[i].hash) {
            stack.pop_back();
        }
        minimal[i] = minimal[i] && (stack.empty() || kmers[stack.back()].offset - kmers[i].offset > radius);
        stack.push_back(i);
    }

    std::vector<KmerAnchor> anchors;
    for (size_t i = 0; i < kmers.size(); i++) {
        if (minimal[i]) {
            anchors.push_back({ kmers[i].offset, kmers[i].kmer });
        }
    }
    return anchors;
}

//------------------------------------------------------------------------------

namespace {

/// The prime modulus of the fingerprints, 2^61 - 1.
constexpr std::uint64_t FINGERPRINT_PRIME = (std::uint64_t(1) << 61) - 1;

/// The base of the fingerprint polynomial.
constexpr std::uint64_t FINGERPRINT_BASE = 0x0123456789ABCDEFULL;

/// Return `x` modulo the prime, for `x < 2^64`. Since `2^61` is 1 modulo the prime, the bits
/// above the lowest 61 can be added to them.
std::uint64_t reduce_fingerprint(std::uint64_t x) {
    x = (x & FINGERPRINT_PRIME) + (x >> 61);
    return (x >= FINGERPRINT_PRIME ? x - FINGERPRINT_PRIME : x);
}

/// Return `a * b` modulo the prime, for `a` and `b` less than the prime.
std::uint64_t multiply_fingerprints(std::uint64_t a, std::uint64_t b) {
    __uint128_t product = static_cast<__uint128_t>(a) * b;
    std::uint64_t low = static_cast<std::uint64_t>(product) & FINGERPRINT_PRIME;
    std::uint64_t high = static_cast<std::uint64_t>(product >> 61);
    return reduce_fingerprint(low + high);
}

/// Return the fingerprint of the string `s + c`, given the fingerprint of `s`.
std::uint64_t extend_fingerprint(std::uint64_t fingerprint, char c) {
    return reduce_fingerprint(multiply_fingerprints(fingerprint, FINGERPRINT_BASE) + static_cast<unsigned char>(c));
}

}

std::uint64_t sequence_fingerprint(std::string_view sequence) {
    std::uint64_t result = 0;
    for (char c : sequence) {
        result = extend_fingerprint(result, c);
    }
    return result;
}

SubstringFingerprinter::SubstringFingerprinter(std::string_view sequence, size_t max_length) :
    prefix_fingerprints(sequence.size() + 1, 0), base_powers(std::min(max_length, sequence.size()) + 1, 1) {
    for (size_t i = 0; i < sequence.size(); i++) {
        this->prefix_fingerprints[i + 1] = extend_fingerprint(this->prefix_fingerprints[i], sequence[i]);
    }
    for (size_t i = 1; i < this->base_powers.size(); i++) {
        this->base_powers[i] = multiply_fingerprints(this->base_powers[i - 1], FINGERPRINT_BASE);
    }
}

std::uint64_t SubstringFingerprinter::fingerprint(size_t begin, size_t end) const {
    // The prefix of length `end` is the prefix of length `begin` times `B^(end - begin)`, plus
    // the substring.
    std::uint64_t shifted = multiply_fingerprints(this->prefix_fingerprints[begin], this->base_powers[end - begin]);
    return reduce_fingerprint(this->prefix_fingerprints[end] + FINGERPRINT_PRIME - shifted);
}

//------------------------------------------------------------------------------

CutCandidateSet::CutCandidateSet(std::string_view sequence, const ContentChopParameters& parameters) :
    max_length(parameters.max_node_length), fingerprinter(sequence, parameters.max_node_length) {
    // Add the candidate and any necessary filler positions.
    auto add_candidate = [&](size_t offset) {
        if (offset == this->cut_offsets.back()) {
            return;
        }
        while (offset - this->cut_offsets.back() > this->max_length) {
            this->cut_offsets.push_back(this->cut_offsets.back() + this->max_length);
        }
        this->cut_offsets.push_back(offset);
    };

    this->cut_offsets.push_back(0);
    for (const KmerAnchor& anchor : find_anchors(sequence, parameters.kmer_length, parameters.anchor_radius)) {
        add_candidate(anchor.offset);
    }
    add_candidate(sequence.size());
}

//------------------------------------------------------------------------------

SegmentCountTable::SegmentCountTable(std::vector<std::uint64_t> fingerprints) {
    parallel_sort(fingerprints.begin(), fingerprints.end(), std::less<std::uint64_t>());
    for (std::uint64_t fingerprint : fingerprints) {
        if (this->entries.empty() || this->entries.back().fingerprint != fingerprint) {
            entries.push_back({ fingerprint, 0 });
        }
        this->entries.back().count++;
    }
}

size_t SegmentCountTable::count(std::uint64_t fingerprint) const {
    auto iter = std::lower_bound(this->entries.begin(), this->entries.end(), fingerprint, [](const Entry& entry, std::uint64_t value) {
        return entry.fingerprint < value;
    });
    return (iter != this->entries.end() && iter->fingerprint == fingerprint ? iter->count : 0);
}

//------------------------------------------------------------------------------

std::vector<size_t> choose_cuts(
    const CutCandidateSet& candidates,
    const SegmentCountTable& counts,
    double piece_weight
) {
    const std::vector<size_t>& offsets = candidates.offsets();
    size_t max_length = candidates.max_segment_length();

    // Dynamic programming over the offsets to find the least-cost chopping of the sequence.
    std::vector<double> best_cost(offsets.size(), std::numeric_limits<double>::infinity());
    std::vector<size_t> piece_start(offsets.size(), 0); // Previous offset in the best path to this offset.
    best_cost[0] = 0.0;
    size_t first_in_reach = 0;
    for (size_t to = 1; to < offsets.size(); to++) {
        while (offsets[to] - offsets[first_in_reach] > max_length) {
            // We are guaranteed to have a candidate within max_length.
            first_in_reach++;
        }
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
    auto iter = std::lower_bound(
        this->pieces.begin(), this->pieces.end(), std::make_pair(node, node_offset),
        [](const PieceLocation& piece, const std::pair<nid_t, size_t>& key
    ) {
        return std::make_pair(piece.node, piece.node_offset) < key;
    });
    if (iter == this->pieces.end() || iter->node != node || iter->node_offset != node_offset) {
        return nullptr;
    }
    return &*iter;
}

std::string_view ChoppingPlan::piece_sequence(const PieceLocation& piece) const {
    return std::string_view(this->strings[piece.string_id]).substr(piece.string_offset, piece.length);
}

size_t ChoppingPlan::stored_length() const {
    size_t total = 0;
    for (const std::string& s : this->strings) {
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

/// A structure for finding common superstrings of the sequences corresponding
/// to chopped pieces.
class PiecePacker {
public:
    /// Create a packer that reads node sequences from the graph.
    PiecePacker(const HandleGraph& graph, const ContentChopParameters& parameters);

    /// Packs the sequences for the pieces corresponding to the chop points in
    /// the plan. Stores the strings and the piece locations in the plan.
    ///
    /// @pre Chop points are sorted by node ID and offset in strictly
    /// increasing order. All chop points are strictly inside their nodes.
    void pack(ChoppingPlan& plan);

private:
    /// A piece with information on how the sequence is stored.
    struct PackedPiece {
        nid_t node;
        size_t offset;
        size_t length;
        /// True if we have chosen to store this piece directly.
        bool stored = false;
        /// Start position of a stored piece with matching sequence.
        ForwardPosition source = { 0, 0 };
    };

    /// Collects the pieces that the chop points in the plan cut the nodes into.
    /// Drops the results of any earlier packing.
    void collect_pieces(const ChoppingPlan& plan);

    /// Places a piece by:
    /// 1. Finding an already placed piece with identical sequence.
    /// 2. Finding the sequence as a substring of a run of stored pieces.
    /// 3. Storing the piece explicitly.
    void place(size_t piece_id);

    /// Returns true and sets the source of the piece if we have already placed
    /// a piece with identical sequence.
    bool find_equal_piece(PackedPiece& piece, const std::string& sequence, std::uint64_t fingerprint) const;

    // TODO: We currently cannot find occurrences in a run of stored pieces.
    /// Returns true and sets the source of the piece if we can find a run of
    /// stored pieces with this piece as a substring.
    bool find_in_stored_pieces(PackedPiece& piece, const std::string& sequence) const;

    /// Returns true if stored pieces cover interval [begin, end) of the node.
    bool is_stored(nid_t node, size_t begin, size_t end) const;

    /// Stores the piece explicitly.
    void store(PackedPiece& piece);

    /// Returns the index of the piece of the node that contains the offset,
    /// or the number of pieces if there is none.
    size_t piece_containing(nid_t node, size_t offset) const;

    /// Returns the sequence of an interval of a node.
    std::string node_sequence(nid_t node, size_t offset, size_t length) const;

    /// Maximum number of seed occurrences to extend.
    static constexpr size_t MAX_TRIED_PLACES = 32;

    const HandleGraph& graph;
    size_t kmer_length;
    size_t anchor_radius;

    /// The pieces, ordered by node and offset.
    std::vector<PackedPiece> pieces;

    /// Placed pieces with each fingerprint.
    std::unordered_map<std::uint64_t, std::vector<size_t>> pieces_by_fingerprint;

    /// Stored pieces indexed by anchors in them.
    std::unordered_map<std::uint64_t, std::vector<ForwardPosition>> stored_anchors;
};

PiecePacker::PiecePacker(const HandleGraph& graph, const ContentChopParameters& parameters) :
    graph(graph), kmer_length(parameters.kmer_length), anchor_radius(parameters.anchor_radius) {
}

void PiecePacker::pack(ChoppingPlan& plan) {
    collect_pieces(plan);

    // Place long pieces first and hope that we can find short pieces as their substrings.
    std::vector<size_t> order(this->pieces.size());
    for (size_t i = 0; i < order.size(); i++) {
        order[i] = i;
    }
    std::sort(order.begin(), order.end(), [&](size_t a, size_t b) {
        return std::make_tuple(this->pieces[b].length, this->pieces[a].node, this->pieces[a].offset)
            < std::make_tuple(this->pieces[a].length, this->pieces[b].node, this->pieces[b].offset);
    });
    for (size_t piece_id : order) {
        this->place(piece_id);
    }

    // Create a string for each run of consecutive pieces that we have chosen to store.
    // We process the pieces in sorted order by node and offset.
    struct StoredRun {
        nid_t node;
        size_t start_offset;
        size_t end_offset;
    };
    std::vector<StoredRun> runs;
    for (const PackedPiece& piece : this->pieces) {
        if (!piece.stored) {
            continue;
        }
        if (!runs.empty() && runs.back().node == piece.node && runs.back().end_offset == piece.offset) {
            runs.back().end_offset += piece.length;
        } else {
            runs.push_back({ piece.node, piece.offset, piece.offset + piece.length });
        }
    }
    plan.strings.clear();
    plan.strings.reserve(runs.size());
    for (const StoredRun& run : runs) {
        plan.strings.push_back(node_sequence(run.node, run.start_offset, run.end_offset - run.start_offset));
    }

    plan.pieces.clear();
    plan.pieces.reserve(this->pieces.size());
    for (const PackedPiece& piece : this->pieces) {
        // The run that contains the source is the last one starting at or before it.
        auto iter = std::upper_bound(runs.begin(), runs.end(), piece.source, [](const ForwardPosition& pos, const StoredRun& run) {
            return std::make_pair(pos.node, pos.offset) < std::make_pair(run.node, run.start_offset);
        });
        --iter;
        size_t string_id = iter - runs.begin();
        plan.pieces.push_back({ piece.node, piece.offset, string_id, piece.source.offset - iter->start_offset, piece.length });
    }
}

void PiecePacker::collect_pieces(const ChoppingPlan& plan) {
    this->pieces.clear();
    this->pieces_by_fingerprint.clear();
    this->stored_anchors.clear();
    for (const auto& [node, node_cuts] : plan.chop_points) {
        size_t length = this->graph.get_length(this->graph.get_handle(node));
        size_t start = 0;
        for (size_t cut : node_cuts) {
            this->pieces.push_back({ node, start, cut - start });
            start = cut;
        }
        this->pieces.push_back({ node, start, length - start });
    }
}

void PiecePacker::place(size_t piece_id) {
    PackedPiece& piece = this->pieces[piece_id];
    std::string sequence = this->node_sequence(piece.node, piece.offset, piece.length);
    std::uint64_t fingerprint = sequence_fingerprint(sequence);
    if (!this->find_equal_piece(piece, sequence, fingerprint) && !this->find_in_stored_pieces(piece, sequence)) {
        this->store(piece);
    }
    this->pieces_by_fingerprint[fingerprint].push_back(piece_id);
}

bool PiecePacker::find_equal_piece(PackedPiece& piece, const std::string& sequence, std::uint64_t fingerprint) const {
    auto iter = this->pieces_by_fingerprint.find(fingerprint);
    if (iter == this->pieces_by_fingerprint.end()) {
        return false;
    }
    for (size_t other_id : iter->second) {
        const PackedPiece& other = this->pieces[other_id];
        // Equal fingerprints almost always mean equal sequences, but we make sure.
        if (other.length == piece.length && this->node_sequence(other.node, other.offset, other.length) == sequence) {
            piece.source = other.source;
            return true;
        }
    }
    return false;
}

bool PiecePacker::find_in_stored_pieces(PackedPiece& piece, const std::string& sequence) const {
    // Find stored pieces by anchors in this piece and select the rarest anchor as a seed.
    const std::vector<ForwardPosition>* rarest_places = nullptr;
    size_t rarest_offset = 0;
    for (const KmerAnchor& anchor : find_anchors(sequence, kmer_length, anchor_radius)) {
        if (anchor.offset < anchor_radius || anchor.offset + anchor_radius + kmer_length > sequence.size()) {
            // If the entire window is not in the sequence, we could end up with a
            // different anchor in a full window.
            continue;
        }
        auto iter = this->stored_anchors.find(anchor.kmer);
        if (iter == this->stored_anchors.end()) {
            // All anchors must be present if this is a substring of a stored piece.
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

    // Try to extend the seed to a full match.
    size_t places_to_try = std::min(rarest_places->size(), MAX_TRIED_PLACES);
    for (size_t i = 0; i < places_to_try; i++) {
        const ForwardPosition& place = (*rarest_places)[i];
        if (place.offset < rarest_offset) {
            // The seed is too close to the start of other node.
            continue;
        }
        size_t begin = place.offset - rarest_offset;
        if (!is_stored(place.node, begin, begin + sequence.size())) {
            // Stored pieces do not cover the entire interval.
            continue;
        }
        if (this->node_sequence(place.node, begin, sequence.size()) == sequence) {
            piece.source = { place.node, begin };
            return true;
        }
    }
    return false;
}

bool PiecePacker::is_stored(nid_t node, size_t begin, size_t end) const {
    size_t i = this->piece_containing(node, begin);
    while (begin < end) {
        if (i >= this->pieces.size() || this->pieces[i].node != node || !this->pieces[i].stored) {
            return false;
        }
        begin = this->pieces[i].offset + this->pieces[i].length;
        i++;
    }
    return true;
}

void PiecePacker::store(PackedPiece& piece) {
    piece.stored = true;
    piece.source = { piece.node, piece.offset };

    // Index the piece by anchors of the underlying node such that the anchor starts
    // within the piece.
    size_t node_length = this->graph.get_length(this->graph.get_handle(piece.node));
    size_t window_begin = (piece.offset > anchor_radius ? piece.offset - anchor_radius : 0);
    size_t window_end = std::min(piece.offset + piece.length + anchor_radius + kmer_length - 1, node_length);
    std::string window = this->node_sequence(piece.node, window_begin, window_end - window_begin);
    for (const KmerAnchor& anchor : find_anchors(window, kmer_length, anchor_radius)) {
        size_t offset = window_begin + anchor.offset;
        if (offset >= piece.offset && offset < piece.offset + piece.length) {
            stored_anchors[anchor.kmer].push_back({ piece.node, offset });
        }
    }
}

size_t PiecePacker::piece_containing(nid_t node, size_t offset) const {
    // The piece containing the offset is the last one starting at or before it.
    auto iter = std::upper_bound(
        this->pieces.begin(), this->pieces.end(), std::make_pair(node, offset),
        [](const std::pair<nid_t, size_t>& key, const PackedPiece& piece) {
            return key < std::make_pair(piece.node, piece.offset);
        }
    );
    if (iter == this->pieces.begin()) {
        return this->pieces.size();
    }
    --iter;
    if (iter->node != node || offset >= iter->offset + iter->length) {
        return this->pieces.size();
    }
    return iter - this->pieces.begin();
}

std::string PiecePacker::node_sequence(nid_t node, size_t offset, size_t length) const {
    return this->graph.get_subsequence(this->graph.get_handle(node), offset, length);
}

}

//------------------------------------------------------------------------------

ChoppingPlan plan_content_chopping(const HandleGraph& graph, const ContentChopParameters& parameters) {
    parameters.validate();

    // Ignore nodes that are already short enough.
    std::vector<nid_t> long_nodes;
    graph.for_each_handle([&](const handle_t& handle) {
        if (graph.get_length(handle) > parameters.max_node_length) {
            long_nodes.push_back(graph.get_id(handle));
        }
    });
    std::sort(long_nodes.begin(), long_nodes.end());

    // Count the segments of the long nodes.
    std::vector<std::uint64_t> fingerprints;
    #pragma omp parallel
    {
        std::vector<std::uint64_t> thread_fingerprints;
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

    // Choose the chopping points.
    std::vector<std::vector<size_t>> cuts(long_nodes.size());
    #pragma omp parallel for schedule(dynamic, 1)
    for (size_t i = 0; i < long_nodes.size(); i++) {
        std::string sequence = graph.get_sequence(graph.get_handle(long_nodes[i]));
        // We could in principle reuse the candidate sets from the counting step,
        // but recomputing them is cheap enough.
        CutCandidateSet candidates(sequence, parameters);
        cuts[i] = choose_cuts(candidates, counts, parameters.piece_weight);
    }

    ChoppingPlan plan;
    for (size_t i = 0; i < long_nodes.size(); i++) {
        if (!cuts[i].empty()) {
            plan.chop_points.emplace_back(long_nodes[i], cuts[i]);
        }
    }

    // Build the stored strings and the piece locations.
    PiecePacker packer(graph, parameters);
    packer.pack(plan);

    return plan;
}

//------------------------------------------------------------------------------

}
