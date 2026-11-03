/**
 * \file unitigs.cpp
 * Haplotype unitigs, see unitigs.hpp.
 */

#include "unitigs.hpp"

#include <algorithm>
#include <stdexcept>

#include <omp.h>

namespace vg {

using handlegraph::handle_t;
using gbwtgraph::GBWTGraph;

//------------------------------------------------------------------------------

namespace {

struct Extension {
    handle_t handle;
    gbwt::BidirectionalState state;
};

/**
 * A walk with the walk coordinate where each of its nodes starts.
 *
 * It finds the nodes that hold any stretch of the walk's bases.
 */
struct MeasuredWalk {
    std::vector<handle_t> walk;
    /// Where each node of the walk starts, followed by the length of the walk.
    std::vector<size_t> starts { 0 };

    MeasuredWalk(const GBWTGraph& graph, std::vector<handle_t> handles) : walk(std::move(handles)) {
        for (const handle_t& handle : walk) {
            starts.push_back(starts.back() + graph.get_length(handle));
        }
    }

    /// Return the length of the walk.
    size_t length() const {
        return starts.back();
    }
    /// Add a node to the end of the walk.
    void push_back(const GBWTGraph& graph, const handle_t& handle) {
        walk.push_back(handle);
        starts.push_back(starts.back() + graph.get_length(handle));
    }
    /// Return the indexes of the first node that holds bases [from, to) and of
    /// the node after the last one.
    std::pair<size_t, size_t> nodes(size_t from, size_t to) const {
        size_t i = std::upper_bound(starts.begin(), starts.end(), from) - starts.begin() - 1;
        size_t j = std::lower_bound(starts.begin(), starts.end(), to) - starts.begin();
        return std::make_pair(i, j);
    }
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
    bool extend(const gbwt::BidirectionalState& state, bool backward, std::vector<Extension>& extensions) const {
        extensions.clear();
        gbwt::size_type extended = 0;
        graph.follow_paths(cache, state, backward, [&](const gbwt::BidirectionalState& next) {
            if (!next.empty()) {
                gbwt::node_type node = backward ? gbwt::Node::reverse(next.backward.node) : next.forward.node;
                extensions.emplace_back(GBWTGraph::node_to_handle(node), next);
                extended += next.size();
            }
            return true;
        });
        return extended < state.size();
    }

    /// Does the k-mer starting at the first base of the walk start a unitig?
    /// The state is for the whole walk.
    bool starts_unitig(const std::vector<handle_t>& walk, size_t length, const gbwt::BidirectionalState& state,
                       std::vector<Extension>& extensions) const {
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

    /// Return the state of the k-mer that starts at the given base of the walk.
    gbwt::BidirectionalState state_of(const MeasuredWalk& walk, size_t kmer) const {
        auto [from, past] = walk.nodes(kmer, kmer + k);
        std::vector<handle_t> pattern(walk.walk.begin() + from, walk.walk.begin() + past);
        return graph.bd_find(cache, pattern);
    }

    /// Trace the unitig whose first k-mer starts at the given offset of the walk,
    /// and report it if this is its canonical end.
    void trace(std::vector<handle_t> handles, size_t first) {
        MeasuredWalk walk(graph, std::move(handles));
        size_t support = state_of(walk, first).size();
        size_t last = extend_to_end(walk, first);
        if (!is_canonical_end(walk, first, last)) {
            return;
        }
        statistics.unitigs++;
        auto [from, past] = walk.nodes(first, last + k);
        WalkInterval unitig { std::vector<handle_t>(walk.walk.begin() + from, walk.walk.begin() + past),
                              first - walk.starts[from], last + k - first };
        report_tigs(unitig, support);
    }

    /// Extend the walk to the end of the unitig whose first k-mer starts at
    /// base first, and return where its last k-mer starts.
    size_t extend_to_end(MeasuredWalk& walk, size_t first) const {
        size_t last = first;
        std::vector<Extension> extensions;
        while (true) {
            // Skip to the next k-mer where the unitig can end: one that ends at the
            // end of the walk, or one whose successor starts at the start of a node.
            size_t next = walk.length() - k;
            if (auto node_start = std::upper_bound(walk.starts.begin(), walk.starts.end() - 1, last);
                node_start != walk.starts.end() - 1) {
                next = std::min(next, *node_start - 1);
            }
            last = next;
            if (last + k == walk.length()) {
                if (extend(state_of(walk, last), false, extensions) || extensions.size() != 1) {
                    return last;
                }
                walk.push_back(graph, extensions.front().handle);
            }
            if (std::binary_search(walk.starts.begin(), walk.starts.end() - 1, last + 1) &&
                (extend(state_of(walk, last + 1), true, extensions) || extensions.size() != 1)) {
                return last;
            }
            last++;
        }
    }

    /// Return true if we report the unitig from this end, because its first
    /// k-mer is smaller than the reverse of its last. The first k-mer starts at
    /// base first of the walk, and the last at base last.
    bool is_canonical_end(const MeasuredWalk& walk, size_t first, size_t last) const {
        auto [head_from, head_past] = walk.nodes(first, first + k);
        auto [tail_from, tail_past] = walk.nodes(last, last + k);
        std::vector<uint64_t> forward_key { GBWTGraph::handle_to_node(walk.walk[head_from]),
                                            first - walk.starts[head_from] };
        for (size_t i = head_from + 1; i < head_past; i++) {
            forward_key.push_back(GBWTGraph::handle_to_node(walk.walk[i]));
        }
        std::vector<uint64_t> reverse_key { gbwt::Node::reverse(GBWTGraph::handle_to_node(walk.walk[tail_past - 1])),
                                            walk.starts[tail_past] - (last + k) };
        for (size_t i = tail_past - 1; i-- > tail_from;) {
            reverse_key.push_back(gbwt::Node::reverse(GBWTGraph::handle_to_node(walk.walk[i])));
        }
        return !(reverse_key < forward_key);
    }

    /// Report the tigs of the unitig. The k-mers within a run of at least k
    /// bases along the gref are in the gref sequence, so with a gref the tigs
    /// are the stretches between such runs.
    void report_tigs(const WalkInterval& unitig, size_t support) {
        std::vector<GrefIndex::Run> runs;
        if (gref != nullptr) {
            runs = gref->runs(unitig);
        }
        // Positions are in unitig coordinates.
        size_t piece_start = 0;
        size_t last_kmer = unitig.length - k;
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
        return a.interval.compare(b.interval) < 0;
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
