#include <algorithm>

#include "phase_table.hpp"

namespace vg {

unordered_map<size_t, size_t> PhaseTable::index() const {
    unordered_map<size_t, size_t> out;
    for (size_t i = 0; i < phase_calls.size(); ++i) {
        out[phase_calls[i].record_key] = i;
    }
    return out;
}

void PhaseTable::freeze_for_render(bool keep) {
    // Built from the calls the linkage pass accumulated, as read phasing left them. Sites with no
    // line are included; they are simply never looked up.
    render_phases.clear();
    if (!keep) {
        return;
    }
    render_phases.reserve(phase_calls.size() * 2);
    for (const PhaseCall& pc : phase_calls) {
        // Where a site has more than one PhaseCall, the last one written wins.
        render_phases[pc.record_key] = pc;
    }
}

const PhaseTable::PhaseCall* PhaseTable::rendered(size_t record_key) const {
    const auto found = render_phases.find(record_key);
    return found == render_phases.end() ? nullptr : &found->second;
}

vector<int> PhaseTable::phase_ordered_genotype(size_t record_key,
                                               const vector<int>& genotype) const {
    vector<int> ordered = genotype;
    if (ordered.size() != 2) {
        return ordered;
    }
    const auto found = render_phases.find(record_key);
    // Only on an exact reversal. A PhaseCall that is not a permutation of the chosen pair is left
    // alone, as `emit_variant` refuses to apply one. For a homozygote the swap changes nothing.
    if (found != render_phases.end() && found->second.ploidy == 2
        && found->second.trav_first == ordered[1]
        && found->second.trav_second == ordered[0]) {
        std::swap(ordered[0], ordered[1]);
    }
    return ordered;
}

int PhaseTable::haploid_slot(size_t record_key, const vector<int>& genotype) const {
    if (genotype.size() != 1) {
        return 0;
    }
    const auto found = render_phases.find(record_key);
    if (found == render_phases.end() || found->second.ploidy != 1
        || found->second.nested_strand < 0) {
        // No nested strand means a haploid locus, such as chrY or a region given ploidy 1, where
        // slot 1 means nothing.
        return 0;
    }
    // Only where the phase names the allele chosen for this site, as `phase_ordered_genotype` and
    // `emit_variant` require.
    if (found->second.trav_first != genotype[0]) {
        return 0;
    }
    return (int)found->second.nested_strand;
}

size_t PhaseTable::swap_strands(const unordered_set<size_t>& flips, vector<NestedLink> links) {
    // Each site's call, the last one winning, as in `freeze_for_render`.
    const unordered_map<size_t, size_t> phase_index = index();
    // The genotype is the same two traversals either way, so no call changes, only which strand
    // carries which allele.
    for (size_t key : flips) {
        const auto found = phase_index.find(key);
        if (found == phase_index.end()) {
            continue;
        }
        PhaseCall& pc = phase_calls[found->second];
        std::swap(pc.trav_first, pc.trav_second);
        std::swap(pc.allele_first, pc.allele_second);
        std::swap(pc.hap_first, pc.hap_second);
    }

    // A nested site's `nested_strand` was set from its parent's chosen pair when the linkage pass
    // resolved its level, so swapping the parent leaves it naming the other strand. Sites are
    // visited top-down by level, so a parent is done before its children, and each inverts
    // its strand where the meaning of its parent's strand 0 changed:
    //   * under a diploid parent, strand 0 is the parent's first allele, so it changed if the
    //     parent was swapped;
    //   * under a haploid parent, strand 0 is the grandparent's, so it changed if the parent's own
    //     `nested_strand` inverted.
    std::stable_sort(links.begin(), links.end(), [](const NestedLink& a, const NestedLink& b) {
        return a.level < b.level;
    });
    std::unordered_map<size_t, bool> frame_flipped;
    frame_flipped.reserve(links.size() * 2);
    size_t moved = 0;
    for (const NestedLink& link : links) {
        const auto index = phase_index.find(link.key);
        if (index == phase_index.end()) {
            continue;
        }
        PhaseCall& pc = phase_calls[index->second];
        bool parent_flipped = false;
        const auto at = frame_flipped.find(link.parent);
        if (at != frame_flipped.end()) {
            parent_flipped = at->second;
        }
        bool strand_moved = false;
        if (pc.nested_strand >= 0 && parent_flipped) {
            pc.nested_strand = pc.nested_strand == 0 ? 1 : 0;
            // The haplotype is held in the slot `nested_strand` names, and the other slot holds the
            // wildcard, which the mosaic reads as an empty strand.
            std::swap(pc.hap_first, pc.hap_second);
            strand_moved = true;
            ++moved;
        }
        frame_flipped[link.key] = pc.ploidy == 2 ? (flips.count(link.key) != 0) : strand_moved;
    }
    return moved;
}

}
