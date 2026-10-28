/**
 * \file unitigs.cpp
 * Haplotype unitigs, see unitigs.hpp.
 */

#include "unitigs.hpp"

#include <gbwtgraph/utils.h>

#include <algorithm>
#include <cctype>
#include <stdexcept>

#include <omp.h>

namespace vg {

using handlegraph::handle_t;
using handlegraph::nid_t;
using gbwtgraph::GBWTGraph;

//------------------------------------------------------------------------------

GrefIndex::GrefIndex(const gbwtgraph::GBZ& gbz, const std::string& name_prefix) : graph(gbz.graph) {
    const gbwt::GBWT& index = gbz.index;
    if (graph.get_node_count() > 0) {
        min_id = graph.min_node_id();
        positions.resize(graph.max_node_id() - min_id + 1);
    }
    if (!index.hasMetadata() || !index.metadata.hasPathNames()) {
        return;
    }

    gbwtgraph::sample_name_set reference_samples = gbwtgraph::parse_reference_samples_tag(index);
    for (gbwt::size_type path_id = 0; path_id < index.metadata.paths(); path_id++) {
        // A gref copy of a generic path has no sample, so match the whole name.
        handlegraph::PathSense sense = gbwtgraph::get_path_sense(index, path_id, reference_samples);
        std::string name = gbwtgraph::compose_path_name(index, path_id, sense);
        if (name.compare(0, name_prefix.size(), name_prefix) != 0) {
            continue;
        }
        uint64_t shift = 0;
        handlegraph::subrange_t subrange = gbwtgraph::get_path_subrange(index, path_id, sense);
        if (subrange != handlegraph::PathMetadata::NO_SUBRANGE) {
            // Report coordinates on the full path.
            shift = subrange.first;
            name = name.substr(0, name.rfind('['));
        }

        uint32_t rank = names.size();
        names.push_back(name);
        shifts.push_back(shift);
        step_starts.emplace_back();
        step_handles.emplace_back();
        uint64_t coord = 0;
        for (gbwt::node_type node : index.extract(gbwt::Path::encode(path_id, false))) {
            handle_t handle = GBWTGraph::node_to_handle(node);
            Position& position = positions[graph.get_id(handle) - min_id];
            if (position.path == NO_PATH) {
                position.path = rank;
                position.is_reverse = graph.get_is_reverse(handle);
                position.offset = coord;
            } else {
                // The gref should be node-disjoint. Keep the first visit.
                repeats++;
            }
            step_starts.back().push_back(coord);
            step_handles.back().push_back(handle);
            coord += graph.get_length(handle);
        }
    }
}

const GrefIndex::Position& GrefIndex::find(nid_t id) const {
    static const Position nowhere;
    if (id < min_id || size_t(id - min_id) >= positions.size()) {
        return nowhere;
    }
    return positions[id - min_id];
}

char GrefIndex::base(uint32_t path, uint64_t coord) const {
    const std::vector<uint64_t>& starts = step_starts[path];
    size_t step = std::upper_bound(starts.begin(), starts.end(), coord) - starts.begin() - 1;
    return graph.get_sequence_view(step_handles[path][step])[coord - starts[step]];
}

uint64_t GrefIndex::total_length() const {
    uint64_t total = 0;
    for (size_t path = 0; path < names.size(); path++) {
        if (!step_starts[path].empty()) {
            total += step_starts[path].back() + graph.get_length(step_handles[path].back());
        }
    }
    return total;
}

//------------------------------------------------------------------------------

namespace {

char complement(char base) {
    switch (std::toupper(base)) {
    case 'A': return 'T';
    case 'C': return 'G';
    case 'G': return 'C';
    case 'T': return 'A';
    default: return 'N';
    }
}

std::string reverse_complement(const std::string& sequence) {
    std::string result(sequence.rbegin(), sequence.rend());
    for (char& base : result) {
        base = complement(base);
    }
    return result;
}

/// A stretch of tig bases that lies contiguously along one gref path.
struct GrefRun {
    size_t begin, end;
    uint32_t path;
    /// Path coordinates decrease along the tig.
    bool is_reverse;
    /// Path coordinate of the base at begin.
    uint64_t coord;

    uint64_t coord_at(size_t position) const {
        return is_reverse ? coord - (position - begin) : coord + (position - begin);
    }
    uint64_t last() const { return coord_at(end - 1); }
};

/// A placement of part of a tig on one gref path.
struct Segment {
    uint32_t path;
    bool is_reverse;
    uint64_t low, high;
    /// Path coordinate the next base would have if there were no edit.
    int64_t next;
    std::vector<std::pair<size_t, char>> cigar;

    void add(size_t length, char op) {
        if (length == 0) {
            return;
        }
        if (!cigar.empty() && cigar.back().second == op) {
            cigar.back().first += length;
        } else {
            cigar.emplace_back(length, op);
        }
    }
};

struct Extension {
    handle_t handle;
    gbwt::BidirectionalState state;
};

/**
 * Finds the unitigs that start next to one node, in both orientations. Each
 * unitig is found from both of its ends, and reported from the end whose first
 * k-mer is smaller.
 */
class UnitigFinder {
public:
    UnitigFinder(const GBWTGraph& graph, size_t k, const GrefIndex* gref,
                 std::vector<Tig>& output, UnitigStatistics& statistics) :
        graph(graph), k(k), gref(gref), cache(graph.get_cache()), output(output), statistics(statistics) {
    }

    void visit(const handle_t& handle) {
        for (handle_t oriented : { handle, graph.flip(handle) }) {
            starting_at(oriented);
            ending_at(oriented);
        }
    }

private:
    const GBWTGraph& graph;
    size_t k;
    const GrefIndex* gref;
    gbwt::CachedGBWT cache;
    std::vector<Tig>& output;
    UnitigStatistics& statistics;

    /// Find the one-node extensions of the state. Returns true if some
    /// occurrences cannot be extended, because their haplotype ends there.
    bool extend(const gbwt::BidirectionalState& state, bool backward, std::vector<Extension>& extensions) {
        extensions.clear();
        gbwt::size_type extended = 0;
        graph.follow_paths(cache, state, backward, [&](const gbwt::BidirectionalState& next) -> bool {
            if (!next.empty()) {
                gbwt::node_type node = backward ? gbwt::Node::reverse(next.backward.node) : next.forward.node;
                extensions.push_back({ GBWTGraph::node_to_handle(node), next });
                extended += next.size();
            }
            return true;
        });
        return extended < state.size();
    }

    /// Does the k-mer starting at the first base of the walk start a unitig?
    /// The state is for the whole walk.
    bool starts_unitig(const std::vector<handle_t>& walk, size_t length, const gbwt::BidirectionalState& state,
                       std::vector<Extension>& extensions) {
        if (extend(state, true, extensions) || extensions.size() != 1) {
            return true;
        }
        // The predecessor is unique. It can only branch if it ends at the end of a
        // node, which happens when the k-mer ends at the first base of a node.
        if (length - graph.get_length(walk.back()) != k - 1) {
            return false;
        }
        std::vector<handle_t> predecessor { extensions.front().handle };
        predecessor.insert(predecessor.end(), walk.begin(), walk.end() - 1);
        return extend(graph.bd_find(cache, predecessor), false, extensions) || extensions.size() != 1;
    }

    /// Trace the unitigs starting with a k-mer at the first base of the node.
    void starting_at(const handle_t& handle) {
        gbwt::BidirectionalState state = graph.get_bd_state(cache, handle);
        if (state.empty()) {
            return;
        }
        struct Frame {
            std::vector<handle_t> walk;
            size_t length;
            gbwt::BidirectionalState state;
        };
        std::vector<Frame> stack;
        stack.push_back({ { handle }, graph.get_length(handle), state });
        std::vector<Extension> extensions;
        while (!stack.empty()) {
            Frame frame = std::move(stack.back());
            stack.pop_back();
            if (frame.length < k) {
                extend(frame.state, false, extensions);
                for (const Extension& extension : extensions) {
                    Frame next { frame.walk, frame.length + graph.get_length(extension.handle), extension.state };
                    next.walk.push_back(extension.handle);
                    stack.push_back(std::move(next));
                }
            } else if (starts_unitig(frame.walk, frame.length, frame.state, extensions)) {
                trace(std::move(frame.walk), 0);
            }
        }
    }

    /// Trace the unitigs starting with a k-mer that ends at the first base of
    /// the node and does not start at the first base of a node. Such a k-mer
    /// has a unique predecessor, which ends at the end of a node.
    void ending_at(const handle_t& handle) {
        gbwt::BidirectionalState state = graph.get_bd_state(cache, handle);
        if (state.empty()) {
            return;
        }
        struct Frame {
            std::vector<handle_t> walk;
            size_t before;
            gbwt::BidirectionalState state;
        };
        std::vector<Frame> stack;
        stack.push_back({ { handle }, 0, state });
        std::vector<Extension> extensions;
        while (!stack.empty()) {
            Frame frame = std::move(stack.back());
            stack.pop_back();
            if (frame.before < k - 1) {
                extend(frame.state, true, extensions);
                for (const Extension& extension : extensions) {
                    Frame next { { extension.handle }, frame.before + graph.get_length(extension.handle), extension.state };
                    next.walk.insert(next.walk.end(), frame.walk.begin(), frame.walk.end());
                    stack.push_back(std::move(next));
                }
                continue;
            }
            size_t offset = frame.before - (k - 1);
            if (offset == 0) {
                // Handled by starting_at().
                continue;
            }
            std::vector<handle_t> predecessor(frame.walk.begin(), frame.walk.end() - 1);
            if (extend(graph.bd_find(cache, predecessor), false, extensions) || extensions.size() != 1) {
                trace(std::move(frame.walk), offset);
            }
        }
    }

    /// Trace the unitig whose first k-mer starts at the given offset of the walk,
    /// and report it if this is its canonical end.
    void trace(std::vector<handle_t> walk, size_t first) {
        // starts[i] is where walk[i] starts in the walk; starts.back() is the length.
        std::vector<size_t> starts { 0 };
        for (const handle_t& handle : walk) {
            starts.push_back(starts.back() + graph.get_length(handle));
        }
        // Indexes of the nodes overlapping [from, to).
        auto overlapping = [&](size_t from, size_t to) {
            size_t i = std::upper_bound(starts.begin(), starts.end(), from) - starts.begin() - 1;
            size_t j = std::lower_bound(starts.begin(), starts.end(), to) - starts.begin();
            return std::make_pair(i, j);
        };
        auto state_of = [&](size_t kmer) {
            auto nodes = overlapping(kmer, kmer + k);
            std::vector<handle_t> pattern(walk.begin() + nodes.first, walk.begin() + nodes.second);
            return graph.bd_find(cache, pattern);
        };

        size_t support = state_of(first).size();
        size_t last = first;
        std::vector<Extension> extensions;
        while (true) {
            // Skip to the next k-mer where the unitig can end: one that ends at the
            // end of the walk, or one whose successor starts at the start of a node.
            size_t next = starts.back() - k;
            auto node_start = std::upper_bound(starts.begin(), starts.end() - 1, last);
            if (node_start != starts.end() - 1) {
                next = std::min(next, *node_start - 1);
            }
            last = next;
            if (last + k == starts.back()) {
                if (extend(state_of(last), false, extensions) || extensions.size() != 1) {
                    break;
                }
                walk.push_back(extensions.front().handle);
                starts.push_back(starts.back() + graph.get_length(walk.back()));
            }
            if (std::binary_search(starts.begin(), starts.end() - 1, last + 1)) {
                if (extend(state_of(last + 1), true, extensions) || extensions.size() != 1) {
                    break;
                }
            }
            last++;
        }

        // Report the unitig from the end whose first k-mer is smaller.
        auto head = overlapping(first, first + k);
        auto tail = overlapping(last, last + k);
        std::vector<uint64_t> forward_key { GBWTGraph::handle_to_node(walk[head.first]), first - starts[head.first] };
        for (size_t i = head.first + 1; i < head.second; i++) {
            forward_key.push_back(GBWTGraph::handle_to_node(walk[i]));
        }
        std::vector<uint64_t> reverse_key { gbwt::Node::reverse(GBWTGraph::handle_to_node(walk[tail.second - 1])),
                                            starts[tail.second] - (last + k) };
        for (size_t i = tail.second - 1; i-- > tail.first;) {
            reverse_key.push_back(gbwt::Node::reverse(GBWTGraph::handle_to_node(walk[i])));
        }
        if (reverse_key < forward_key) {
            return;
        }
        statistics.unitigs++;

        // Find where the unitig lies on the gref.
        size_t end = last + k;
        std::vector<GrefRun> runs;
        if (gref != nullptr) {
            auto nodes = overlapping(first, end);
            for (size_t i = nodes.first; i < nodes.second; i++) {
                const GrefIndex::Position& position = gref->find(graph.get_id(walk[i]));
                if (position.path == GrefIndex::NO_PATH) {
                    continue;
                }
                size_t from = std::max(starts[i], first), to = std::min(starts[i + 1], end);
                size_t node_length = starts[i + 1] - starts[i];
                bool is_reverse = (graph.get_is_reverse(walk[i]) != position.is_reverse);
                uint64_t coord = is_reverse ? position.offset + node_length - 1 - (from - starts[i])
                                            : position.offset + (from - starts[i]);
                if (!runs.empty() && runs.back().end == from && runs.back().path == position.path &&
                    runs.back().is_reverse == is_reverse && runs.back().coord_at(from) == coord) {
                    runs.back().end = to;
                } else {
                    runs.push_back({ from, to, position.path, is_reverse, coord });
                }
            }
        }

        // The k-mers within a run are in the gref sequence. Report the stretches between them.
        size_t piece_start = first;
        for (const GrefRun& run : runs) {
            if (run.end - run.begin < k) {
                continue;
            }
            if (run.begin > piece_start) {
                report(walk, starts, piece_start, run.begin - 1, support, runs);
            }
            piece_start = std::max(piece_start, run.end - k + 1);
        }
        if (piece_start <= last) {
            report(walk, starts, piece_start, last, support, runs);
        }
    }

    /// Report the tig consisting of the k-mers starting at [first, last].
    void report(const std::vector<handle_t>& walk, const std::vector<size_t>& starts,
                size_t first, size_t last, size_t support, const std::vector<GrefRun>& runs) {
        size_t end = last + k;
        size_t i = std::upper_bound(starts.begin(), starts.end(), first) - starts.begin() - 1;
        size_t j = std::lower_bound(starts.begin(), starts.end(), end) - starts.begin();

        Tig tig;
        tig.walk.assign(walk.begin() + i, walk.begin() + j);
        tig.offset = first - starts[i];
        tig.support = support;
        tig.sequence.reserve(end - first);
        for (size_t n = i; n < j; n++) {
            std::string_view view = graph.get_sequence_view(walk[n]);
            size_t from = std::max(starts[n], first) - starts[n];
            size_t to = std::min(starts[n + 1], end) - starts[n];
            tig.sequence.append(view.substr(from, to - from));
        }

        // Runs in tig coordinates.
        std::vector<GrefRun> local;
        for (const GrefRun& run : runs) {
            if (run.end <= first || run.begin >= end) {
                continue;
            }
            GrefRun clipped = run;
            clipped.begin = std::max(run.begin, first);
            clipped.end = std::min(run.end, end);
            clipped.coord = run.coord_at(clipped.begin);
            clipped.begin -= first;
            clipped.end -= first;
            local.push_back(clipped);
        }

        if (gref != nullptr) {
            // Use the orientation with the most bases on forward gref strands.
            size_t forward_bases = 0, reverse_bases = 0;
            for (const GrefRun& run : local) {
                (run.is_reverse ? reverse_bases : forward_bases) += run.end - run.begin;
            }
            if (reverse_bases > forward_bases) {
                flip(tig, local);
            }
            describe(tig, local);
        }

        statistics.tigs++;
        statistics.bases += tig.sequence.size();
        statistics.kmers += tig.sequence.size() - k + 1;
        statistics.longest = std::max(statistics.longest, tig.sequence.size());
        output.push_back(std::move(tig));
    }

    void flip(Tig& tig, std::vector<GrefRun>& runs) const {
        size_t walk_length = 0;
        for (const handle_t& handle : tig.walk) {
            walk_length += graph.get_length(handle);
        }
        size_t length = tig.sequence.size();
        std::reverse(tig.walk.begin(), tig.walk.end());
        for (handle_t& handle : tig.walk) {
            handle = graph.flip(handle);
        }
        tig.offset = walk_length - (tig.offset + length);
        tig.sequence = reverse_complement(tig.sequence);
        std::reverse(runs.begin(), runs.end());
        for (GrefRun& run : runs) {
            uint64_t coord = run.last();
            size_t begin = length - run.end;
            run.end = length - run.begin;
            run.begin = begin;
            run.coord = coord;
            run.is_reverse = !run.is_reverse;
        }
    }

    /// Describe where the tig lies on the gref, as segments with CIGARs.
    void describe(Tig& tig, const std::vector<GrefRun>& runs) {
        std::vector<Segment> segments;
        size_t position = 0;
        for (const GrefRun& run : runs) {
            size_t inserted = run.begin - position;
            bool continues = false;
            int64_t deleted = 0;
            if (!segments.empty() && segments.back().path == run.path && segments.back().is_reverse == run.is_reverse) {
                deleted = run.is_reverse ? segments.back().next - int64_t(run.coord) : int64_t(run.coord) - segments.back().next;
                continues = (deleted >= 0);
            }
            if (continues) {
                Segment& segment = segments.back();
                if (inserted == size_t(deleted)) {
                    // A substitution: compare the bases.
                    for (size_t t = 0; t < inserted; t++) {
                        int64_t coord = run.is_reverse ? segment.next - int64_t(t) : segment.next + int64_t(t);
                        char base = gref->base(run.path, coord);
                        if (run.is_reverse) {
                            base = complement(base);
                        }
                        segment.add(1, std::toupper(tig.sequence[position + t]) == std::toupper(base) ? '=' : 'X');
                    }
                } else {
                    segment.add(deleted, 'D');
                    segment.add(inserted, 'I');
                }
            } else {
                if (!segments.empty()) {
                    // Bases between placements on different paths go with the earlier one.
                    segments.back().add(inserted, 'I');
                }
                segments.push_back({ run.path, run.is_reverse, run.coord, run.coord, 0, {} });
                if (segments.size() == 1) {
                    segments.back().add(inserted, 'I');
                }
            }
            Segment& segment = segments.back();
            segment.add(run.end - run.begin, '=');
            segment.low = std::min({ segment.low, run.coord, run.last() });
            segment.high = std::max({ segment.high, run.coord, run.last() });
            segment.next = run.is_reverse ? int64_t(run.last()) - 1 : int64_t(run.last()) + 1;
            position = run.end;
        }
        size_t trailing = tig.sequence.size() - position;

        if (segments.empty()) {
            tig.gref = "*:" + std::to_string(trailing) + "I";
            statistics.unanchored++;
            return;
        }
        segments.back().add(trailing, 'I');
        tig.sort_path = segments.front().path;
        tig.sort_coord = segments.front().low;
        for (Segment& segment : segments) {
            if (segment.is_reverse) {
                // CIGARs are in path order.
                std::reverse(segment.cigar.begin(), segment.cigar.end());
            }
            if (!tig.gref.empty()) {
                tig.gref += ',';
            }
            uint64_t shift = gref->path_shift(segment.path) + 1;
            tig.gref += gref->path_name(segment.path) + ':' + std::to_string(segment.low + shift) + '-' +
                        std::to_string(segment.high + shift) + ':' + (segment.is_reverse ? '-' : '+') + ':';
            for (const auto& op : segment.cigar) {
                tig.gref += std::to_string(op.first) + op.second;
            }
        }
    }
};

}

//------------------------------------------------------------------------------

std::vector<Tig> haplotype_unitigs(const GBWTGraph& graph, size_t k, const GrefIndex* gref,
                                   UnitigStatistics* statistics) {
    if (k < 2) {
        throw std::invalid_argument("haplotype_unitigs(): k must be at least 2");
    }
    if (!graph.index->bidirectional()) {
        throw std::invalid_argument("haplotype_unitigs(): the GBWT must be bidirectional");
    }

    size_t threads = omp_get_max_threads();
    std::vector<std::vector<Tig>> results(threads);
    std::vector<UnitigStatistics> thread_statistics(threads);
    graph.for_each_handle([&](const handle_t& handle) {
        int thread = omp_get_thread_num();
        UnitigFinder finder(graph, k, gref, results[thread], thread_statistics[thread]);
        finder.visit(handle);
    }, true);

    std::vector<Tig> tigs;
    for (std::vector<Tig>& result : results) {
        std::move(result.begin(), result.end(), std::back_inserter(tigs));
        result = std::vector<Tig>();
    }
    // Sort along the gref, then by the walk, so that the order does not depend on threads.
    std::sort(tigs.begin(), tigs.end(), [](const Tig& a, const Tig& b) {
        if (a.sort_path != b.sort_path) {
            return a.sort_path < b.sort_path;
        }
        if (a.sort_coord != b.sort_coord) {
            return a.sort_coord < b.sort_coord;
        }
        for (size_t i = 0; i < std::min(a.walk.size(), b.walk.size()); i++) {
            gbwt::node_type x = GBWTGraph::handle_to_node(a.walk[i]), y = GBWTGraph::handle_to_node(b.walk[i]);
            if (x != y) {
                return x < y;
            }
        }
        if (a.walk.size() != b.walk.size()) {
            return a.walk.size() < b.walk.size();
        }
        if (a.offset != b.offset) {
            return a.offset < b.offset;
        }
        return a.sequence.size() < b.sequence.size();
    });

    if (statistics != nullptr) {
        *statistics = UnitigStatistics();
        for (const UnitigStatistics& thread : thread_statistics) {
            statistics->unitigs += thread.unitigs;
            statistics->tigs += thread.tigs;
            statistics->bases += thread.bases;
            statistics->kmers += thread.kmers;
            statistics->longest = std::max(statistics->longest, thread.longest);
            statistics->unanchored += thread.unanchored;
        }
    }
    return tigs;
}

}
