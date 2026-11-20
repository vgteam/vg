/**
 * \file unitigs.cpp
 * Haplotype unitigs, see unitigs.hpp.
 */

#include "unitigs.hpp"

#include <algorithm>
#include <memory>
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

/// Return the sort key of a k-mer, which orders k-mers by where they start and
/// then by where they go: the GBWT node of the first node of the walk, the
/// offset of the first base in that node, then the GBWT nodes of the rest.
std::vector<uint64_t> kmer_key(const WalkInterval& kmer) {
    std::vector<uint64_t> key { GBWTGraph::handle_to_node(kmer.walk.front()), kmer.offset };
    for (size_t i = 1; i < kmer.walk.size(); i++) {
        key.push_back(GBWTGraph::handle_to_node(kmer.walk[i]));
    }
    return key;
}

/**
 * Finds the unitigs that start next to one node, in both orientations. Each
 * unitig is found from both of its ends, and reported from the end whose first
 * k-mer is smaller.
 */
class UnitigFinder {
public:
    UnitigFinder(const GBWTGraph& graph, size_t k, const GrefIndex* gref,
                 std::vector<Tig>& output, size_t& unitigs) :
        graph(graph), k(k), gref(gref), cache(graph.get_cache()), output(output), unitigs(unitigs) {
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
    /// Number of unitigs found, before the cutting that turns them into tigs.
    size_t& unitigs;

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

    /**
     * Enumerate the walks that start from the node and reach the needed
     * length, extending one node at a time: to the right of the node, or to
     * the left of it when backward is true.
     *
     * Call report(walk, length, state, extensions) on each walk that reaches
     * the needed length. The length is the length of the walk, and the state
     * is the GBWT state of the whole walk. The callback may use extensions as
     * scratch space.
     */
    template<typename Report>
    void for_each_walk(const handle_t& handle, bool backward, size_t needed, Report report) {
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
            if (frame.length >= needed) {
                report(frame.walk, frame.length, frame.state, extensions);
                continue;
            }
            extend(frame.state, backward, extensions);
            for (const Extension& extension : extensions) {
                Frame next { frame.walk, frame.length + graph.get_length(extension.handle), extension.state };
                next.walk.insert(backward ? next.walk.begin() : next.walk.end(), extension.handle);
                stack.push_back(std::move(next));
            }
        }
    }

    /// Trace the unitigs starting with a k-mer at the first base of the node.
    void starting_at(const handle_t& handle) {
        for_each_walk(handle, false, k, [&](std::vector<handle_t>& walk, size_t length,
                                            const gbwt::BidirectionalState& state,
                                            std::vector<Extension>& extensions) {
            if (starts_unitig(walk, length, state, extensions)) {
                trace(walk, 0);
            }
        });
    }

    /// Trace the unitigs starting with a k-mer that ends at the first base of
    /// the node and does not start at the first base of a node. Such a k-mer
    /// has a unique predecessor, which ends at the end of a node.
    void ending_at(const handle_t& handle) {
        // The k-mer ends at the first base of the node, so the walk needs k - 1
        // bases before the node.
        size_t node_length = graph.get_length(handle);
        for_each_walk(handle, true, node_length + k - 1, [&](std::vector<handle_t>& walk, size_t length,
                                                             const gbwt::BidirectionalState&,
                                                             std::vector<Extension>& extensions) {
            size_t offset = length - node_length - (k - 1);
            if (offset == 0) {
                // Handled by starting_at().
                return;
            }
            std::vector<handle_t> predecessor(walk.begin(), walk.end() - 1);
            if (extend(graph.bd_find(cache, predecessor), false, extensions) || extensions.size() != 1) {
                trace(walk, offset);
            }
        });
    }

    /// Return the state of the k-mer that starts at the given base of the walk.
    gbwt::BidirectionalState state_of(const PathSteps& walk, size_t kmer) const {
        return graph.bd_find(cache, walk.interval(kmer, kmer + k).walk);
    }

    /// Trace the unitig whose first k-mer starts at the given offset of the walk,
    /// and report it if this is its canonical end.
    void trace(const std::vector<handle_t>& handles, size_t first) {
        PathSteps walk(graph, handles);
        size_t support = state_of(walk, first).size();
        size_t last = extend_to_end(walk, first);
        if (!is_canonical_end(walk, first, last)) {
            return;
        }
        unitigs++;
        report_tigs(walk.interval(first, last + k), support);
    }

    /// Extend the walk to the end of the unitig whose first k-mer starts at
    /// base first, and return where its last k-mer starts.
    size_t extend_to_end(PathSteps& walk, size_t first) const {
        size_t last = first;
        std::vector<Extension> extensions;
        while (true) {
            // Skip to the next k-mer where the unitig can end: one that ends at the
            // end of the walk, or one whose successor starts at the start of a node.
            size_t after = walk.step_at(last) + 1;
            uint64_t next_node = after < walk.steps().size() ? walk.steps()[after].start : walk.length();
            last = std::min<size_t>(walk.length() - k, next_node - 1);
            if (last + k == walk.length()) {
                if (extend(state_of(walk, last), false, extensions) || extensions.size() != 1) {
                    return last;
                }
                walk.push_back(graph, extensions.front().handle);
            }
            // A successor that starts at the start of a node can have another
            // predecessor, which would end the unitig here.
            if (walk.steps()[walk.step_at(last + 1)].start == last + 1 &&
                (extend(state_of(walk, last + 1), true, extensions) || extensions.size() != 1)) {
                return last;
            }
            last++;
        }
    }

    /// Return true if we report the unitig from this end, because its first
    /// k-mer is not larger than the reverse of its last. The first k-mer starts
    /// at base first of the walk, and the last at base last.
    bool is_canonical_end(const PathSteps& walk, size_t first, size_t last) const {
        WalkInterval head = walk.interval(first, first + k);
        WalkInterval tail = walk.interval(last, last + k).reverse(graph);
        return !(kmer_key(tail) < kmer_key(head));
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
        }
        output.push_back(std::move(tig));
    }
};

}

//------------------------------------------------------------------------------

std::vector<Tig> haplotype_unitigs(const GBWTGraph& graph, size_t k, const GrefIndex* gref,
                                   size_t* unitig_count) {
    if (k < 2) {
        throw std::invalid_argument("haplotype_unitigs(): k must be at least 2");
    }
    if (!graph.index->bidirectional()) {
        throw std::invalid_argument("haplotype_unitigs(): the GBWT must be bidirectional");
    }

    size_t threads = omp_get_max_threads();
    std::vector<std::vector<Tig>> results(threads);
    std::vector<size_t> thread_unitigs(threads, 0);
    // One finder per thread, not one per node: a finder holds a gbwt::CachedGBWT,
    // whose record cache is worth keeping warm across the nodes a thread visits.
    std::vector<std::unique_ptr<UnitigFinder>> finders(threads);
    graph.for_each_handle([&](const handle_t& handle) {
        int thread = omp_get_thread_num();
        std::unique_ptr<UnitigFinder>& finder = finders[thread];
        if (!finder) {
            finder.reset(new UnitigFinder(graph, k, gref, results[thread], thread_unitigs[thread]));
        }
        finder->visit(handle);
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

    if (unitig_count != nullptr) {
        *unitig_count = 0;
        for (size_t found : thread_unitigs) {
            *unitig_count += found;
        }
    }
    return tigs;
}

}
