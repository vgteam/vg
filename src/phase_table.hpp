#ifndef VG_PHASE_TABLE_HPP_INCLUDED
#define VG_PHASE_TABLE_HPP_INCLUDED

#include <cstdint>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include "linkage_model.hpp"

namespace vg {

using namespace std;

/**
 * The phase of every phased site: which of the site's chosen alleles each strand carries.
 *
 * The linkage model fills it one level at a time, parents before their children, and read
 * phasing then swaps the strands at some sites. Before any record is rendered, the calls are
 * frozen into a copy keyed by record, which the render, the anchors and the record's GT read.
 */
class PhaseTable {
public:
    using PhaseCall = LinkageCollector::PhaseCall;

    /// A nested site's place in the nesting tree: its record key, its parent's, and its level.
    struct NestedLink {
        size_t key = 0;
        size_t parent = 0;
        uint8_t level = 0;
    };

    /// Every phase call, in the order the linkage model produced them. A site can have more
    /// than one; the last one is the site's phase.
    vector<PhaseCall>& calls() { return phase_calls; }
    const vector<PhaseCall>& calls() const { return phase_calls; }

    /// Each site's index in `calls()`, the last one winning where a site has more than one.
    unordered_map<size_t, size_t> index() const;

    /// Swap which strand carries which allele at each site in `flips`, and keep each nested
    /// site's strand pointing at the parent strand that carries it. A ploidy-1 nested site names
    /// one of its parent's two strands in `nested_strand` and holds its haplotype in the slot of
    /// that number, so where the parent's strands swapped, both move to the other strand; under a
    /// haploid parent, they move where the parent's own strand moved. `links` places each site in
    /// the nesting tree, in any order; a site left out stops the swap reaching its children.
    /// Returns how many nested strands moved.
    size_t swap_strands(const unordered_set<size_t>& flips, vector<NestedLink> links);

    /// Copy the calls into the table the render reads, keyed by record, the last call winning
    /// where a site has more than one. With `keep` false the render table is left empty, so
    /// nothing rendered is phased.
    void freeze_for_render(bool keep);

    /// Whether the render table has any call.
    bool has_rendered() const { return !render_phases.empty(); }

    /// The call the render table holds for a site, or null.
    const PhaseCall* rendered(size_t record_key) const;

    /// The chosen pair in phase order, for the anchors, which take each slot from the order of
    /// the pair they are given.
    ///
    /// `LinkageCollector::chosen_traversals` returns a sorted pair; the phase is in the render
    /// table. The pair is swapped when the record's PhaseCall names the same two traversals in
    /// the other order, and returned unchanged otherwise: no phasing, no PhaseCall, or a
    /// PhaseCall naming other traversals.
    vector<int> phase_ordered_genotype(size_t record_key, const vector<int>& genotype) const;

    /// Which strand a one-allele genotype sits on, for the anchors: 0 or 1.
    ///
    /// A nested chain at ploidy 1 is one strand of its parent, the one `nested_strand` names;
    /// `emit_variant` writes it as `a|.` or `.|a`, and this keeps the anchor's slot the same.
    /// Returns 0 when there is no phasing, no entry, a ploidy other than 1, no nested strand, or a
    /// phase that names a different allele from the chosen one.
    int haploid_slot(size_t record_key, const vector<int>& genotype) const;

private:
    vector<PhaseCall> phase_calls;
    /// `phase_calls` keyed by record, frozen by `freeze_for_render`. Keyed by record rather than
    /// by (contig, POS), since POS depends on which alleles the line carries.
    unordered_map<size_t, PhaseCall> render_phases;
};

}

#endif
