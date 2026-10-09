#include <algorithm>
#include <limits>

#include "child_placer.hpp"
#include "block_records.hpp"
#include "site_values.hpp"
#include "snarls.hpp"
#include "symbolic_allele.hpp"

namespace vg {

ChildPlacer::TraversalNodeIndex ChildPlacer::index_traversal_nodes(const HandleGraph& graph,
                                                                   const Traversal& walk) {
    TraversalNodeIndex visits;
    for (int i = 0; i < (int)walk.size(); ++i) {
        visits[graph.get_id(walk[i])].push_back(i);
    }
    return visits;
}

int ChildPlacer::crossings_of_child(const HandleGraph& graph, const TraversalNodeIndex& visits,
                                   const SiteBounds& child) {
    const nid_t start = graph.get_id(child.start);
    const nid_t end = graph.get_id(child.end);
    // Count crossings: an entry at one boundary followed by the other. Order matters: testing for
    // the two boundaries separately would count a walk that touches both on unrelated excursions,
    // as find_child_traversal_set does. Only visits to the two boundary nodes can change the
    // state, so we walk their positions merged in ascending order. When both boundaries are the
    // same node, one list plays both roles.
    static const vector<int> none;
    auto s = visits.find(start);
    auto e = visits.find(end);
    const vector<int>& sp = (s == visits.end()) ? none : s->second;
    const vector<int>& ep = (e == visits.end() || end == start) ? none : e->second;
    int crossings = 0;
    nid_t open = 0;
    size_t i = 0, j = 0;
    while (i < sp.size() || j < ep.size()) {
        nid_t node;
        if (j >= ep.size() || (i < sp.size() && sp[i] <= ep[j])) {
            node = start;
            ++i;
        } else {
            node = end;
            ++j;
        }
        if (open == 0 && (node == start || node == end)) {
            open = (node == start) ? end : start;
        } else if (open != 0 && node == open) {
            ++crossings;
            open = 0;
        }
    }
    return crossings;
}

int ChildPlacer::offset_of_child(const HandleGraph& graph, const Traversal& walk,
                                 const SiteBounds& child) {
    const nid_t start = graph.get_id(child.start);
    const nid_t end = graph.get_id(child.end);
    nid_t open = 0;
    int entry = -1;
    for (int i = 0; i < (int)walk.size(); ++i) {
        nid_t node = graph.get_id(walk[i]);
        if (open == 0 && (node == start || node == end)) {
            open = (node == start) ? end : start;
            entry = i;
        } else if (open != 0 && node == open) {
            return entry;   // first complete crossing, entry side
        }
    }
    return -1;
}

ChildPlacer::ChildOffsets::ChildOffsets(const HandleGraph& graph, const Traversal& walk) {
    bases_before.assign(walk.size() + 1, 0);
    for (int i = 0; i < (int)walk.size(); ++i) {
        visits_of[graph.get_id(walk[i])].push_back(i);
        bases_before[i + 1] = bases_before[i] + (int64_t)graph.get_length(walk[i]);
    }
}

int64_t ChildPlacer::ChildOffsets::base_offset(const HandleGraph& graph,
                                               const SiteBounds& child) const {
    // `offset_of_child`'s rule: the entry is the first visit to either boundary node, and it counts
    // only if the other boundary node is visited after it.
    const nid_t start = graph.get_id(child.start);
    const nid_t end = graph.get_id(child.end);
    auto start_visits = visits_of.find(start);
    auto end_visits = visits_of.find(end);
    const int none = numeric_limits<int>::max();
    const int first_start = start_visits == visits_of.end() ? none : start_visits->second.front();
    const int first_end = end_visits == visits_of.end() ? none : end_visits->second.front();
    const int entry = min(first_start, first_end);
    if (entry == none) {
        return -1;
    }
    const auto& closing = first_start <= first_end ? end_visits : start_visits;
    if (closing == visits_of.end()
        || std::upper_bound(closing->second.begin(), closing->second.end(), entry)
               == closing->second.end()) {
        return -1;
    }
    return bases_before[entry];
}

size_t ChildPlacer::offset_along_genotype(
    const HandleGraph& graph, const vector<Traversal>& travs, const vector<int>& genotype,
    const SiteBounds& child, unordered_map<const Traversal*, ChildOffsets>& offsets) {
    for (int allele : genotype) {
        if (allele < 0 || allele >= (int)travs.size()) {
            continue;
        }
        const Traversal* walk = &travs[allele];
        auto found = offsets.find(walk);
        if (found == offsets.end()) {
            found = offsets.emplace(walk, ChildOffsets(graph, *walk)).first;
        }
        const int64_t within = found->second.base_offset(graph, child);
        if (within >= 0) {
            return (size_t)within;
        }
    }
    return 0;
}

uint64_t ChildPlacer::child_crossing_mask(const HandleGraph& graph,
                                         const vector<TraversalNodeIndex>& visits,
                                         const SiteBounds& child, bool* known) {
    if (known != nullptr) {
        *known = true;
    }
    // One bit per candidate walk, not per VCF allele, since the linkage model chooses a site's
    // genotype as a pair of walks.
    if (visits.size() > 64) {
        // Unknown rather than none: the mask cannot index this site's candidates.
        if (known != nullptr) {
            *known = false;
        }
        return 0;
    }
    uint64_t mask = 0;
    for (size_t i = 0; i < visits.size(); ++i) {
        if (crossings_of_child(graph, visits[i], child) > 0) {
            mask |= (uint64_t)1 << i;
        }
    }
    return mask;
}

int64_t ChildPlacer::base_offset_of_child(const Traversal& walk, const SiteBounds& child) const {
    const int entry = ChildPlacer::offset_of_child(*graph, walk, child);
    if (entry < 0) {
        return -1;
    }
    int64_t bases = 0;
    for (int i = 0; i < entry && i < (int)walk.size(); ++i) {
        bases += (int64_t)graph->get_length(walk[i]);
    }
    return bases;
}

size_t ChildPlacer::offset_along_genotype(const vector<Traversal>& travs,
                                          const vector<int>& genotype,
                                          const SiteBounds& child) const {
    for (int allele : genotype) {
        if (allele < 0 || allele >= (int)travs.size()) {
            continue;
        }
        const int64_t within = base_offset_of_child(travs[allele], child);
        if (within >= 0) {
            return (size_t)within;
        }
    }
    return 0;
}

int ChildPlacer::child_ploidy(const vector<TraversalNodeIndex>& visits,
                              const vector<int>& genotype, const SiteBounds& child,
                              int cap) const {
    int copies = 0;
    bool capped = false;

    for (int allele : genotype) {
        if (allele < 0 || allele >= (int)visits.size()) {
            continue;   // star or missing: that haplotype contributes no copy here
        }
        int crossings = crossings_of_child(*graph, visits[allele], child);
        if (crossings > 1) {
            capped = true;
            crossings = 1;   // a cycle or tandem duplication; see the header comment
        }
        copies += crossings;
    }
    if (capped) {
        // Counted and reported once per run.
        ++counters->child_multi_crossing;
    }
    return min(copies, cap);
}

void ChildPlacer::configure(const HandleGraph* graph, const SnarlManager* manager,
                            const BlockRecordWriter* blocks, DescentCounters* counters) {
    this->graph = graph;
    this->manager = manager;
    this->blocks = blocks;
    this->counters = counters;
}

bool ChildPlacer::reported_inline(const BlockRecordWriter::ChainInlineContext& ctx,
                                  const Snarl& child) const {
    // The chain of the snarl a walk enters by the child's start.
    const Snarl* entered = manager->into_which_snarl(child.start().node_id(),
                                                     child.start().backward());
    return entered != nullptr
           && blocks->chain_reported_inline(*graph, ctx, chain_of_site(*manager, *graph, entered));
}

void ChildPlacer::place(const Snarl& site, const SiteChildren& children, size_t site_key,
                        const vector<Traversal>& travs,
                        const vector<int>& genotype, int ref_trav_idx, int ploidy,
                        const NestingPlacement& placement, bool off_reference,
                        bool keep_uncrossed, const function<void(const Placed& child)>& visit) const {
    const Snarl* managed_ptr = manager->into_which_snarl(site.start().node_id(),
                                                         site.start().backward());
    if (managed_ptr == nullptr) {
        return;
    }
    // The child-independent parts of the exactly-once test, built once for this snarl.
    const BlockRecordWriter::ChainInlineContext inline_ctx =
        blocks->chain_inline_context(*graph, children, travs, genotype, ref_trav_idx);
    // Also once for this snarl: see TraversalNodeIndex.
    vector<TraversalNodeIndex> trav_visits;
    trav_visits.reserve(travs.size());
    for (const Traversal& t : travs) {
        trav_visits.push_back(index_traversal_nodes(*graph, t));
    }
    for (const Snarl* child : manager->children_of(managed_ptr)) {
        if (child == nullptr || manager->is_trivial(child, *graph)) {
            continue;
        }
        const SiteBounds child_bounds = bounds_of(*graph, *child);
        // A chain that no reference path passes through has no REF or POS for its records,
        // so it is skipped unless off-reference descent is on.
        bool child_off_reference = false;
        if (ref_trav_idx >= 0 && ref_trav_idx < (int)travs.size()) {
            vector<int> ref_only(1, ref_trav_idx);
            if (child_ploidy(trav_visits, ref_only, child_bounds, 1) == 0) {
                // With off-reference descent, such a chain is genotyped and recorded but has
                // no line.
                if (!off_reference) {
                    ++counters->skipped_no_ref;
                    continue;
                }
                child_off_reference = true;
                ++counters->off_reference;
            }
        }
        // Inherited: everything under a chain the reference does not cross is also off it.
        if (placement.no_reference) {
            child_off_reference = true;
        }

        // The exactly-once test: under block emission, a chain that every called strand
        // crosses only inside a difference block is already spelled by that block's ALT. It
        // holds back the chain's line, not its descent, so the chain is still genotyped,
        // recorded and phased. Inherited by chains inside it. Does nothing when block
        // emission is off, or for a snarl whose projection has no symbols.
        bool child_reported_inline =
            placement.reported_inline || reported_inline(inline_ctx, *child);

        int copies = child_ploidy(trav_visits, genotype, child_bounds, ploidy);
        bool retain_only = placement.retain_only;
        if (copies <= 0) {
            // No called allele reaches it yet. Visited anyway, while this window's reads are
            // in memory, since the linkage model may move the parent onto an allele that
            // does reach it. Nothing about it is written unless the linkage pass says so.
            ++counters->skipped_no_copy;
            if (!keep_uncrossed) {
                // Without retention there is nothing to come back to. Without the linkage
                // model nothing moves the parent after the direct pass, so the sample has no copy
                // of this chain; the linkage pass, which has no chosen parent to read, would
                // otherwise render it at the parent's ploidy.
                continue;
            }
            retain_only = true;
        }

        Placed out;
        out.snarl = child;
        out.chain = chain_of_site(*manager, *graph, child);
        // `copies` is zero only for a chain no called parent allele reaches, which is still
        // genotyped; it then takes the parent's ploidy, the most copies a child can have.
        out.ploidy = copies >= 1 ? copies : ploidy;
        NestingPlacement& next = out.placement;
        next.one_copy = (copies == 1);
        next.parent_record_key = site_key;
        next.retain_only = retain_only;
        next.no_reference = child_off_reference;
        // Where this child starts along the first called allele that reaches it, added to
        // the offset of its parent. Only an off-reference chain uses it, but it is computed
        // for every chain, so that offsets add up down the tree.
        next.parent_offset =
            placement.parent_offset + offset_along_genotype(travs, genotype, child_bounds);
        next.reported_inline = child_reported_inline;
        // The chain's identity, from its boundary nodes.
        {
            const nid_t first = graph->get_id(out.chain.start);
            const nid_t second = graph->get_id(out.chain.end);
            next.chain_key = (size_t)((uint64_t)first * 1000003ULL) ^ (size_t)(uint64_t)second;
        }
        bool crossing_known = true;   // child_crossing_mask always sets it
        // The mask is over this snarl's own candidate traversals, which exist whether or not
        // a line was written.
        next.parent_crossing =
            child_crossing_mask(*graph, trav_visits, child_bounds, &crossing_known);
        next.crossing_known = crossing_known;
        next.level = placement.level + 1;
        visit(out);
    }
}

}
