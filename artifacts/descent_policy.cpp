#include <algorithm>
#include <limits>

#include "descent_policy.hpp"
#include "block_records.hpp"
#include "site_values.hpp"
#include "symbolic_allele.hpp"

namespace vg::multipass {

int DescentPolicy::child_ploidy(const vector<ChildCrossingFinder>& walks,
                                const vector<int>& genotype, const net_handle_t& child,
                                int cap) const {
    int copies = 0;
    bool capped = false;

    for (int allele : genotype) {
        if (allele < 0 || allele >= (int)walks.size()) {
            // A star or missing allele carries no copy.
            continue;
        }
        int crossings = (int)walks[allele].crossings(child).size();
        if (crossings > 1) {
            capped = true;
            crossings = 1;
        }
        copies += crossings;
    }
    if (capped) {
        ++counters->child_multi_crossing;
    }
    return min(copies, cap);
}

void DescentPolicy::configure(const HandleGraph* new_graph,
                              const EditSnarlDecompositionOverlay* new_decomposition,
                            const BlockRecordWriter* new_blocks, DescentCounters* new_counters) {
    graph = new_graph;
    decomposition = new_decomposition;
    blocks = new_blocks;
    counters = new_counters;
}

DescentPolicy::Nested DescentPolicy::nested_in(const net_handle_t& site,
                                             const handle_t& called_start,
                                             const handle_t& called_end) const {
    Nested nested;
    nested.children =
        site_children(*decomposition, *graph, site, called_start, called_end, &nested.sites);
    return nested;
}

NodePlacement DescentPolicy::placement_of(nid_t id) const {
    return placement_of_node(*decomposition, *graph, id);
}

bool DescentPolicy::reference_crosses(const Parent& parent,
                                      const vector<ChildCrossingFinder>& walks,
                                      const net_handle_t& child) const {
    if (parent.ref_trav_idx < 0 || parent.ref_trav_idx >= (int)parent.travs.size()) {
        // With no reference walk, there is nothing to test.
        return true;
    }
    const vector<int> ref_only(1, parent.ref_trav_idx);
    return child_ploidy(walks, ref_only, child, 1) != 0;
}

void DescentPolicy::place(const Parent& parent, bool off_reference, bool keep_uncrossed,
                        const function<void(const Placed& child)>& visit) const {
    if (!parent.nested.children.known) {
        return;
    }
    const vector<Traversal>& travs = parent.travs;
    const NestingPlacement& placement = parent.placement;
    // Build once for the parent the parts of the test for an inline-reported chain that do not
    // depend on the child.
    const BlockRecordWriter::ChainInlineContext inline_ctx = blocks->chain_inline_context(
        *graph, parent.nested.children, travs, parent.genotype, parent.ref_trav_idx);
    vector<ChildCrossingFinder> walks;
    walks.reserve(travs.size());
    for (const Traversal& t : travs) {
        walks.emplace_back(*decomposition, *graph, t);
    }
    // Make the crossing finder of each walk the first time a child's offset needs that walk.
    unordered_map<const Traversal*, ChildCrossingFinder> offset_finders;
    for (const net_handle_t& child : parent.nested.sites) {
        const ChildChain chain = chain_holding(*decomposition, *graph, child);
        // A chain the reference walk does not cross has no REF or POS for its records, so leave
        // it out unless off-reference descent is on. Everything under a chain off the reference
        // is off it too.
        bool child_off_reference = placement.no_reference;
        if (!reference_crosses(parent, walks, child)) {
            if (!off_reference) {
                ++counters->skipped_no_ref;
                continue;
            }
            child_off_reference = true;
            ++counters->off_reference;
        }

        // Under block emission, a difference block's ALT already spells a chain that every called
        // strand crosses only inside that block. That holds back the chain's line, but the chain
        // is still genotyped, filed and phased.
        const bool child_reported_inline =
            placement.reported_inline
            || blocks->chain_reported_inline(*graph, inline_ctx, chain);

        const int copies = child_ploidy(walks, parent.genotype, child, parent.ploidy);
        bool retain_only = placement.retain_only;
        if (copies <= 0) {
            ++counters->skipped_no_copy;
            if (!keep_uncrossed) {
                // Without the linkage model nothing moves the parent after the direct pass, so
                // the sample has no copy of this chain.
                continue;
            }
            // Genotype it anyway, while the reads of this window are in memory. Nothing about it
            // is written unless the linkage pass says so.
            retain_only = true;
        }

        Placed out;
        out.site = child;
        out.chain = chain;
        // A chain no called allele crosses takes the parent's ploidy, the most copies it can have.
        out.ploidy = copies >= 1 ? copies : parent.ploidy;
        NestingPlacement& next = out.placement;
        next.one_copy = (copies == 1);
        next.parent_record_key = parent.record_key;
        next.retain_only = retain_only;
        next.no_reference = child_off_reference;
        // Only an off-reference chain uses the offset, but compute it for every chain, so that
        // offsets add up down the tree.
        next.parent_offset =
            placement.parent_offset
            + entry_bases_along_genotype(*decomposition, *graph, travs, parent.genotype, child,
                                         offset_finders);
        next.reported_inline = child_reported_inline;
        next.chain_key = chain_key_of(*graph, out.chain);
        bool crossing_known = true;
        next.parent_crossing = crossing_mask(walks, child, &crossing_known);
        next.crossing_known = crossing_known;
        next.level = placement.level + 1;
        visit(out);
    }
}

}
