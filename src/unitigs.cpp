/**
 * \file unitigs.cpp
 * Haplotype unitigs, see unitigs.hpp.
 */

#include "unitigs.hpp"

#include <gbwtgraph/utils.h>

#include <algorithm>
#include <cctype>
#include <ostream>
#include <sstream>
#include <stdexcept>

#include <omp.h>

namespace vg {

using handlegraph::handle_t;
using handlegraph::nid_t;
using gbwtgraph::GBWTGraph;

//------------------------------------------------------------------------------

std::string WalkInterval::sequence(const GBWTGraph& graph) const {
    std::string result;
    result.reserve(length);
    size_t skip = offset;
    for (const handle_t& handle : walk) {
        std::string_view view = graph.get_sequence_view(handle);
        size_t take = std::min(view.size() - skip, length - result.size());
        result.append(view.substr(skip, take));
        skip = 0;
    }
    return result;
}

WalkInterval WalkInterval::reverse(const GBWTGraph& graph) const {
    WalkInterval result;
    result.walk.reserve(walk.size());
    for (auto it = walk.rbegin(); it != walk.rend(); ++it) {
        result.walk.push_back(graph.flip(*it));
    }
    result.offset = walk_length(graph) - (offset + length);
    result.length = length;
    return result;
}

WalkInterval WalkInterval::subinterval(const GBWTGraph& graph, size_t start, size_t end) const {
    WalkInterval result;
    result.length = end - start;
    // Walk coordinates of the first base and of the base after the last.
    size_t from = offset + start, to = offset + end;
    size_t node_start = 0;
    for (const handle_t& handle : walk) {
        size_t node_end = node_start + graph.get_length(handle);
        if (node_end > from) {
            if (result.walk.empty()) {
                result.offset = from - node_start;
            }
            result.walk.push_back(handle);
        }
        if (node_end >= to) {
            break;
        }
        node_start = node_end;
    }
    return result;
}

size_t WalkInterval::walk_length(const GBWTGraph& graph) const {
    size_t total = 0;
    for (const handle_t& handle : walk) {
        total += graph.get_length(handle);
    }
    return total;
}

void write_gaf_columns(std::ostream& out, const GBWTGraph& graph, const std::string& name,
                       const WalkInterval& interval) {
    out << name << '\t' << interval.length << "\t0\t" << interval.length << "\t+\t";
    for (const handle_t& handle : interval.walk) {
        out << (graph.get_is_reverse(handle) ? '<' : '>') << graph.get_id(handle);
    }
    out << '\t' << interval.walk_length(graph) << '\t' << interval.offset << '\t'
        << (interval.offset + interval.length) << '\t' << interval.length << '\t' << interval.length << "\t255";
}

WalkInterval parse_gaf_columns(const std::string& line, const GBWTGraph& graph, std::string& name) {
    auto fail = [&](const std::string& problem) {
        throw std::runtime_error(problem + " in GAF line: " + line);
    };
    std::vector<std::string> columns;
    std::stringstream stream(line);
    std::string column;
    while (columns.size() < 12 && std::getline(stream, column, '\t')) {
        columns.push_back(column);
    }
    if (columns.size() < 12) {
        fail("fewer than 12 columns");
    }
    auto number = [&](const std::string& text) -> size_t {
        size_t used = 0;
        size_t value = 0;
        try {
            value = std::stoull(text, &used);
        } catch (const std::exception&) {
            used = 0;
        }
        if (used == 0 || used != text.size()) {
            fail("bad number " + text);
        }
        return value;
    };

    name = columns[0];
    WalkInterval interval;
    const std::string& walk = columns[5];
    for (size_t i = 0; i < walk.size();) {
        if (walk[i] != '>' && walk[i] != '<') {
            fail("not a walk of node IDs");
        }
        size_t j = walk.find_first_of("<>", i + 1);
        if (j == std::string::npos) {
            j = walk.size();
        }
        nid_t id = number(walk.substr(i + 1, j - i - 1));
        if (!graph.has_node(id)) {
            fail("node " + std::to_string(id) + " not in the graph");
        }
        interval.walk.push_back(graph.get_handle(id, walk[i] == '<'));
        i = j;
    }
    size_t start = number(columns[7]), end = number(columns[8]);
    if (interval.walk.empty() || start >= end || end > interval.walk_length(graph)) {
        fail("bad walk or path coordinates");
    }
    interval.offset = start;
    interval.length = end - start;
    // Trim the walk to the nodes that hold the interval's bases.
    return interval.subinterval(graph, 0, interval.length);
}

void for_each_path_start(const gbwt::GBWT& index,
                         const std::function<void(gbwt::size_type, const std::string&, uint64_t)>& visit) {
    if (!index.hasMetadata() || !index.metadata.hasPathNames()) {
        return;
    }
    using handlegraph::PathMetadata;
    using handlegraph::PathSense;
    gbwtgraph::sample_name_set reference_samples = gbwtgraph::parse_reference_samples_tag(index);
    for (gbwt::size_type path_id = 0; path_id < index.metadata.paths(); path_id++) {
        PathSense sense = gbwtgraph::get_path_sense(index, path_id, reference_samples);
        std::string sample = gbwtgraph::get_path_sample_name(index, path_id, sense);
        std::string locus = gbwtgraph::get_path_locus_name(index, path_id, sense);
        size_t haplotype = gbwtgraph::get_path_haplotype(index, path_id, sense);
        uint64_t start = 0;
        if (sense == PathSense::HAPLOTYPE) {
            start = gbwtgraph::get_path_phase_block(index, path_id, sense);
            // Compose the name as a reference path name. That is the haplotype
            // path name without the phase block.
            sense = PathSense::REFERENCE;
        } else {
            handlegraph::subrange_t subrange = gbwtgraph::get_path_subrange(index, path_id, sense);
            if (subrange != PathMetadata::NO_SUBRANGE) {
                start = subrange.first;
            }
        }
        std::string name = PathMetadata::create_path_name(sense, sample, locus, haplotype,
                                                          PathMetadata::NO_PHASE_BLOCK, PathMetadata::NO_SUBRANGE);
        visit(path_id, name, start);
    }
}

//------------------------------------------------------------------------------

GrefIndex::GrefIndex(const gbwtgraph::GBZ& gbz, const std::string& name_prefix) : graph(gbz.graph) {
    const gbwt::GBWT& index = gbz.index;
    if (graph.get_node_count() > 0) {
        min_id = graph.min_node_id();
        positions.resize(graph.max_node_id() - min_id + 1);
    }

    for_each_path_start(index, [&](gbwt::size_type path_id, const std::string& name, uint64_t start) {
        if (name.compare(0, name_prefix.size(), name_prefix) != 0) {
            return;
        }
        uint32_t rank = paths.size();
        paths.push_back({ name, start, {} });
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
            paths.back().steps.push_back({ coord, handle });
            coord += graph.get_length(handle);
        }
    });

    // Find each path by the names vg paths gives a path with a subrange,
    // NAME[START] and NAME[START-END], and by NAME alone where it is unique.
    std::unordered_map<std::string, size_t> name_counts;
    for (const Path& path : paths) {
        name_counts[path.name]++;
    }
    for (uint32_t rank = 0; rank < paths.size(); rank++) {
        const Path& path = paths[rank];
        if (name_counts[path.name] == 1) {
            path_ranks[path.name] = rank;
        }
        if (path.shift > 0) {
            std::string start = std::to_string(path.shift);
            path_ranks[path.name + '[' + start + ']'] = rank;
            path_ranks[path.name + '[' + start + '-' + std::to_string(path.shift + path_length(rank)) + ']'] = rank;
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

uint64_t GrefIndex::path_length(uint32_t path) const {
    const std::vector<Step>& steps = paths[path].steps;
    return steps.empty() ? 0 : steps.back().start + graph.get_length(steps.back().handle);
}

uint32_t GrefIndex::find_path(const std::string& name) const {
    auto found = path_ranks.find(name);
    return found == path_ranks.end() ? NO_PATH : found->second;
}

namespace {

/// Return the index of the step that holds the given path coordinate.
template<typename Steps>
size_t step_at(const Steps& steps, uint64_t coord) {
    auto after = std::upper_bound(steps.begin(), steps.end(), coord, [](uint64_t c, const auto& step) {
        return c < step.start;
    });
    return after - steps.begin() - 1;
}

}

WalkInterval GrefIndex::interval(uint32_t path, uint64_t start, uint64_t end) const {
    const std::vector<Step>& steps = paths[path].steps;
    size_t first = step_at(steps, start), last = step_at(steps, end - 1);
    WalkInterval result;
    for (size_t step = first; step <= last; step++) {
        result.walk.push_back(steps[step].handle);
    }
    result.offset = start - steps[first].start;
    result.length = end - start;
    return result;
}

char GrefIndex::base(uint32_t path, uint64_t coord) const {
    const std::vector<Step>& steps = paths[path].steps;
    const Step& step = steps[step_at(steps, coord)];
    return graph.get_sequence_view(step.handle)[coord - step.start];
}

uint64_t GrefIndex::total_length() const {
    uint64_t total = 0;
    for (uint32_t path = 0; path < paths.size(); path++) {
        total += path_length(path);
    }
    return total;
}

std::vector<GrefIndex::Run> GrefIndex::runs(const WalkInterval& interval) const {
    std::vector<Run> result;
    size_t node_start = 0;
    for (const handle_t& handle : interval.walk) {
        size_t node_length = graph.get_length(handle);
        const Position& position = find(graph.get_id(handle));
        // Walk coordinates of the bases of the interval on this node.
        size_t from = std::max(node_start, interval.offset);
        size_t to = std::min(node_start + node_length, interval.offset + interval.length);
        if (position.path != NO_PATH && from < to) {
            bool is_reverse = (graph.get_is_reverse(handle) != position.is_reverse);
            uint64_t coord = is_reverse ? position.offset + node_length - 1 - (from - node_start)
                                        : position.offset + (from - node_start);
            from -= interval.offset;
            to -= interval.offset;
            if (!result.empty() && result.back().end == from && result.back().path == position.path &&
                result.back().is_reverse == is_reverse && result.back().coord_at(from) == coord) {
                result.back().end = to;
            } else {
                result.push_back({ from, to, position.path, is_reverse, coord });
            }
        }
        node_start += node_length;
    }
    return result;
}

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

/// A placement of part of an interval on one gref path.
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

}

GrefPlacement GrefIndex::place(const WalkInterval& interval) const {
    GrefPlacement placement;
    std::vector<Run> runs = this->runs(interval);
    for (const Run& run : runs) {
        (run.is_reverse ? placement.reverse_bases : placement.forward_bases) += run.end - run.begin;
    }
    placement.along_one_path = (runs.size() == 1 && runs.front().begin == 0 && runs.front().end == interval.length);
    // Only bases between two runs can be substitutions, which need the bases
    // to compare with the gref.
    std::string sequence;
    if (runs.size() > 1) {
        sequence = interval.sequence(graph);
    }

    std::vector<Segment> segments;
    size_t position = 0;
    for (const Run& run : runs) {
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
                    char gref_base = base(run.path, coord);
                    if (run.is_reverse) {
                        gref_base = complement(gref_base);
                    }
                    segment.add(1, std::toupper(sequence[position + t]) == std::toupper(gref_base) ? '=' : 'X');
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
    size_t trailing = interval.length - position;

    if (segments.empty()) {
        placement.text = "*:" + std::to_string(trailing) + "I";
        return placement;
    }
    segments.back().add(trailing, 'I');
    placement.path = segments.front().path;
    placement.coord = segments.front().low;
    for (Segment& segment : segments) {
        if (segment.is_reverse) {
            // CIGARs are in path order.
            std::reverse(segment.cigar.begin(), segment.cigar.end());
        }
        if (!placement.text.empty()) {
            placement.text += ',';
        }
        uint64_t shift = path_shift(segment.path) + 1;
        placement.text += path_name(segment.path) + ':' + std::to_string(segment.low + shift) + '-' +
                          std::to_string(segment.high + shift) + ':' + (segment.is_reverse ? '-' : '+') + ':';
        for (const auto& op : segment.cigar) {
            placement.text += std::to_string(op.first) + op.second;
        }
    }
    return placement;
}

//------------------------------------------------------------------------------

namespace {

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

        auto nodes = overlapping(first, last + k);
        WalkInterval unitig { std::vector<handle_t>(walk.begin() + nodes.first, walk.begin() + nodes.second),
                              first - starts[nodes.first], last + k - first };
        std::vector<GrefIndex::Run> runs;
        if (gref != nullptr) {
            runs = gref->runs(unitig);
        }

        // The k-mers within a run are in the gref sequence. Report the stretches
        // between them. Positions are in unitig coordinates from here on.
        size_t piece_start = 0, last_kmer = unitig.length - k;
        for (const GrefIndex::Run& run : runs) {
            if (run.end - run.begin < k) {
                continue;
            }
            if (run.begin > piece_start) {
                report(unitig, piece_start, run.begin - 1, support);
            }
            piece_start = std::max(piece_start, run.end - k + 1);
        }
        if (piece_start <= last_kmer) {
            report(unitig, piece_start, last_kmer, support);
        }
    }

    /// Report the tig made of the k-mers of the unitig that start at positions first to last.
    void report(const WalkInterval& unitig, size_t first, size_t last, size_t support) {
        Tig tig;
        tig.interval = unitig.subinterval(graph, first, last + k);
        tig.support = support;
        if (gref != nullptr) {
            // Use the orientation with the most bases on forward gref strands.
            tig.placement = gref->place(tig.interval);
            if (tig.placement.reverse_bases > tig.placement.forward_bases) {
                tig.interval = tig.interval.reverse(graph);
                tig.placement = gref->place(tig.interval);
            }
            if (tig.placement.path == GrefIndex::NO_PATH) {
                statistics.unanchored++;
            }
        }

        statistics.tigs++;
        statistics.bases += tig.interval.length;
        statistics.kmers += tig.interval.length - k + 1;
        statistics.longest = std::max(statistics.longest, tig.interval.length);
        output.push_back(std::move(tig));
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
        if (a.placement.path != b.placement.path) {
            return a.placement.path < b.placement.path;
        }
        if (a.placement.coord != b.placement.coord) {
            return a.placement.coord < b.placement.coord;
        }
        const WalkInterval& x = a.interval;
        const WalkInterval& y = b.interval;
        for (size_t i = 0; i < std::min(x.walk.size(), y.walk.size()); i++) {
            gbwt::node_type m = GBWTGraph::handle_to_node(x.walk[i]), n = GBWTGraph::handle_to_node(y.walk[i]);
            if (m != n) {
                return m < n;
            }
        }
        if (x.walk.size() != y.walk.size()) {
            return x.walk.size() < y.walk.size();
        }
        if (x.offset != y.offset) {
            return x.offset < y.offset;
        }
        return x.length < y.length;
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

std::vector<std::pair<std::string, size_t>> unsupported_gref_paths(const gbwtgraph::GBZ& gbz,
                                                                   const std::string& name_prefix) {
    std::vector<std::pair<std::string, size_t>> result;
    for_each_path_start(gbz.index, [&](gbwt::size_type path_id, const std::string& name, uint64_t) {
        if (name.compare(0, name_prefix.size(), name_prefix) != 0) {
            return;
        }
        std::vector<handle_t> walk;
        size_t length = 0;
        for (gbwt::node_type node : gbz.index.extract(gbwt::Path::encode(path_id, false))) {
            walk.push_back(GBWTGraph::node_to_handle(node));
            length += gbz.graph.get_length(walk.back());
        }
        // The gref path is one occurrence of its own walk. A second one means some
        // other path contains it, unless the walk is its own reverse, which is too
        // unlikely in a real gref to be worth a slower check.
        if (!walk.empty() && gbz.graph.bd_find(walk).size() < 2) {
            result.emplace_back(name, length);
        }
    });
    return result;
}

}
