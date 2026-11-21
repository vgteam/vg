#include "subchain_walk.hpp"

#include <algorithm>

namespace vg {

SubchainWalk::SubchainWalk(
    const gbwt::CachedGBWT& index, const gbwt::FragmentMap& fragment_map,
    gbwt::node_type start, gbwt::node_type end, Haplotypes::sequence_type sequence
) {
    for_each_visit(index, fragment_map, start, end, sequence, [&](gbwt::node_type node, bool starts_piece) {
        if (starts_piece) {
            this->piece_starts.push_back(static_cast<std::uint32_t>(this->visits.size()));
        }
        this->visits.push_back(node);
        this->reaches_end = (node == end);
        return true;
    });
}

bool SubchainWalk::for_each_visit(
    const gbwt::CachedGBWT& index, const gbwt::FragmentMap& fragment_map,
    gbwt::node_type start, gbwt::node_type end, Haplotypes::sequence_type sequence,
    const std::function<bool(gbwt::node_type, bool)>& visit
) {
    bool single_fragment = (start == gbwt::ENDMARKER || end == gbwt::ENDMARKER);
    gbwt::size_type sequence_id = sequence.first;
    gbwt::edge_type pos = (start != gbwt::ENDMARKER ? gbwt::edge_type(start, sequence.second) : index.start(sequence_id));
    if (pos == gbwt::invalid_edge() || pos.first == gbwt::ENDMARKER) {
        return true;
    }
    bool starts_piece = true;
    while (true) {
        if (!visit(pos.first, starts_piece)) {
            return false;
        }
        if (pos.first == end) {
            return true;
        }
        starts_piece = false;
        pos = index.LF(pos);
        if (pos.first != gbwt::ENDMARKER) {
            continue;
        }
        if (single_fragment) {
            return true;
        }
        // Continue with the next fragment of the haplotype that is not empty.
        while (pos.first == gbwt::ENDMARKER) {
            sequence_id = fragment_map.oriented_next(sequence_id);
            if (sequence_id == gbwt::invalid_sequence()) {
                return true;
            }
            pos = index.start(sequence_id);
        }
        starts_piece = true;
    }
}

std::vector<nid_t> SubchainWalk::node_ids() const {
    std::vector<nid_t> result;
    result.reserve(this->visits.size());
    for (gbwt::node_type node : this->visits) {
        result.push_back(gbwt::Node::id(node));
    }
    gbwt::removeDuplicates(result, false);
    return result;
}

bool SubchainWalk::has_sequence() const {
    return !this->visits.empty() && !(this->reaches_end && this->visits.size() == 1);
}

std::vector<std::string> SubchainWalk::sequence_pieces(const gbwtgraph::GBWTGraph& graph, size_t start_max, size_t end_max) const {
    std::vector<std::string> result;
    if (!this->has_sequence()) {
        return result;
    }
    for (size_t piece = 0; piece < this->piece_starts.size(); piece++) {
        result.emplace_back();
        size_t limit = (piece + 1 < this->piece_starts.size() ? this->piece_starts[piece + 1] : this->visits.size());
        for (size_t i = this->piece_starts[piece]; i < limit; i++) {
            std::string_view view = graph.get_sequence_view(gbwtgraph::GBWTGraph::node_to_handle(this->visits[i]));
            if (i == 0) {
                size_t offset = (view.size() > start_max ? view.size() - start_max : 0);
                result.back().append(view.data() + offset, view.size() - offset);
            } else if (this->reaches_end && i + 1 == this->visits.size()) {
                result.back().append(view.data(), std::min(view.size(), end_max));
            } else {
                result.back().append(view.data(), view.size());
            }
        }
    }
    return result;
}

size_t SubchainWalk::pieces() const {
    return (this->has_sequence() ? this->piece_starts.size() : 0);
}

OrientedKmerTable::OrientedKmerTable(const std::vector<Haplotypes::Subchain::kmer_type>& kmers, size_t k) {
    size_t bits = std::max(size_t(sdsl::bits::length(std::max(4 * kmers.size(), size_t(1)) - 1)), size_t(1));
    this->keys.assign(size_t(1) << bits, EMPTY);
    this->positions.assign(this->keys.size(), NONE);
    this->mask = this->keys.size() - 1;
    this->shift = 64 - bits;
    auto insert = [this](Haplotypes::Subchain::kmer_type kmer, std::uint32_t position) {
        size_t i = this->slot(kmer);
        while (this->keys[i] != EMPTY && this->keys[i] != kmer) {
            i = (i + 1) & this->mask;
        }
        this->keys[i] = kmer;
        this->positions[i] = position;
    };
    for (size_t i = 0; i < kmers.size(); i++) {
        auto position = static_cast<std::uint32_t>(i);
        insert(kmers[i], position);
        insert(gbwtgraph::Key64(kmers[i]).reverse_complement(k).get_key(), position);
    }
}

void OrientedKmerCounts::count(const std::vector<std::string>& pieces, const OrientedKmerTable& positions, size_t k) {
    for (std::uint32_t position : this->occurring_positions) {
        this->counts[position] = 0;
    }
    this->occurring_positions.clear();
    for (const std::string& piece : pieces) {
        gbwtgraph::Key64 key;
        size_t valid_chars = 0;
        for (char c : piece) {
            key.forward(k, static_cast<unsigned char>(c), valid_chars);
            if (valid_chars < k) {
                continue;
            }
            if (std::uint32_t position = positions.find(key.get_key()); position != OrientedKmerTable::NONE) {
                if (this->counts[position]++ == 0) {
                    this->occurring_positions.push_back(position);
                }
            }
        }
    }
}

std::vector<std::pair<std::uint32_t, std::uint32_t>> OrientedKmerCounts::occurrence_list() const {
    std::vector<std::pair<std::uint32_t, std::uint32_t>> result;
    result.reserve(this->occurring_positions.size());
    for (std::uint32_t position : this->occurring_positions) {
        result.emplace_back(position, this->counts[position]);
    }
    return result;
}

} // namespace vg
