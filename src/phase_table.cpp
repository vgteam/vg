#include <algorithm>

#include "phase_table.hpp"

namespace vg::multipass {

unordered_map<size_t, size_t> PhaseTable::index() const {
    unordered_map<size_t, size_t> out;
    for (size_t i = 0; i < phase_calls.size(); ++i) {
        out[phase_calls[i].record_key] = i;
    }
    return out;
}

void PhaseTable::freeze_for_render(bool keep) {
    render_phases.clear();
    if (!keep) {
        return;
    }
    render_phases.reserve(phase_calls.size() * 2);
    for (const PhaseCall& pc : phase_calls) {
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
    // A phase naming any other traversals is ignored, as `emit_variant` ignores it.
    if (const auto found = render_phases.find(record_key);
        found != render_phases.end() && found->second.ploidy == 2
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
        // A haploid site with no nested strand is on a haploid locus, such as chrY, and has
        // just one strand.
        return 0;
    }
    // A phase of another allele is ignored, as `emit_variant` ignores it.
    if (found->second.trav_first != genotype[0]) {
        return 0;
    }
    return (int)found->second.nested_strand;
}

size_t PhaseTable::swap_strands(const unordered_set<size_t>& flips, vector<NestedLink> links) {
    const unordered_map<size_t, size_t> phase_index = index();
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

    // Move each nested site to its other strand where the meaning of its parent's strand 0
    // changed. Under a diploid parent, strand 0 is the parent's first allele, which changed if
    // the parent was swapped. Under a haploid parent, strand 0 is the grandparent's strand 0,
    // which changed if the parent moved. Visit sites by level so a parent is done before its
    // children.
    std::ranges::stable_sort(links, [](const NestedLink& a, const NestedLink& b) {
        return a.level < b.level;
    });
    // Whether each visited site's strand 0 changed meaning.
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
        if (const auto at = frame_flipped.find(link.parent); at != frame_flipped.end()) {
            parent_flipped = at->second;
        }
        bool strand_moved = false;
        if (pc.nested_strand >= 0 && parent_flipped) {
            pc.nested_strand = pc.nested_strand == 0 ? 1 : 0;
            // The haplotype is held in the slot `nested_strand` names, so it moves too.
            std::swap(pc.hap_first, pc.hap_second);
            strand_moved = true;
            ++moved;
        }
        frame_flipped[link.key] = pc.ploidy == 2 ? flips.contains(link.key) : strand_moved;
    }
    return moved;
}

}
