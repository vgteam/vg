#include <algorithm>

#include "flow_caller.hpp"

namespace vg {

void VCFOutputCaller::build_render_phases() {
    // Built from the phasing the linkage pass accumulated, as read phasing left it. Sites with no line
    // are included; they are simply never looked up.
    render_phases.clear();
    if (!emit_phasing) {
        return;
    }
    render_phases.reserve(linkage_phased.size() * 2);
    for (const LinkageCollector::PhaseCall& pc : linkage_phased) {
        // Where a site has more than one PhaseCall, the last one written wins.
        render_phases[pc.record_key] = pc;
    }
}

vector<int> VCFOutputCaller::phase_ordered_genotype(size_t record_key,
                                                    const vector<int>& genotype) const {
    vector<int> ordered = genotype;
    if (!emit_phasing || ordered.size() != 2) {
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

int VCFOutputCaller::phase_haploid_slot(size_t record_key, const vector<int>& genotype) const {
    if (!emit_phasing || genotype.size() != 1) {
        return 0;
    }
    const auto found = render_phases.find(record_key);
    if (found == render_phases.end() || found->second.ploidy != 1
        || found->second.nested_strand < 0) {
        // No nested strand means a haploid locus, such as chrY or a haploid --ploidy-bed region,
        // where slot 1 means nothing.
        return 0;
    }
    // Only where the phase names the allele chosen for this site, as `phase_ordered_genotype` and
    // `emit_variant` require.
    if (found->second.trav_first != genotype[0]) {
        return 0;
    }
    return (int)found->second.nested_strand;
}

size_t FlowCaller::cascade_nested_strands(vector<LinkageCollector::PhaseCall>& phased,
                                          const std::unordered_map<size_t, size_t>& phase_index,
                                          vector<NestedLink> links,
                                          const unordered_set<size_t>& flips) {
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
        LinkageCollector::PhaseCall& pc = phased[index->second];
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
