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

}
