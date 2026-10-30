#ifndef VG_PHASE_TABLE_HPP_INCLUDED
#define VG_PHASE_TABLE_HPP_INCLUDED

#include <cstdint>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include "linkage_model.hpp"

namespace vg::multipass {

using namespace std;

/**
 * A table of the phase of every phased site, saying which chosen allele each strand carries.
 *
 * The linkage model fills it one level at a time, parents before their children. Read phasing
 * then swaps the strands at some sites. Before any record is rendered, `freeze_for_render()`
 * copies the final calls into a second table keyed by record, which record rendering reads.
 */
class PhaseTable {
public:
    using PhaseCall = LinkageCollector::PhaseCall;

    /// A nested site's place in the snarl tree.
    struct NestedLink {
        /// The site's record key.
        size_t key = 0;
        /// The record key of the site's parent site.
        size_t parent = 0;
        /// The level at which the linkage pass chooses the site's genotype.
        uint8_t level = 0;
    };

    /// Get every phase call, in the order the linkage model made them.
    vector<PhaseCall>& calls() { return phase_calls; }
    /// Get every phase call, in the order the linkage model made them.
    const vector<PhaseCall>& calls() const { return phase_calls; }

    /// Map each site's record key to the index in `calls()` of its phase.
    unordered_map<size_t, size_t> index() const;

    /// Swap which strand carries which allele at each site in `flips`.
    ///
    /// A ploidy-1 nested site names the parent strand that carries it in `nested_strand`. When
    /// the meaning of that parent strand changes, the nested site moves to the other strand, so
    /// that it stays with the same parent allele.
    ///
    /// @param links The nested sites, in any order. A site missing from `links` does not pass a
    /// swap on to its children.
    /// @return How many nested sites moved to the other strand.
    size_t swap_strands(const unordered_set<size_t>& flips, vector<NestedLink> links);

    /// Copy each site's phase into the table that record rendering reads.
    ///
    /// With `keep` false the table is left empty, and no record is written phased.
    void freeze_for_render(bool keep);

    /// Return whether the table for record rendering has any phase.
    bool has_rendered() const { return !render_phases.empty(); }

    /// Return the phase that record rendering uses for a site, or null if it has none.
    const PhaseCall* rendered(size_t record_key) const;

    /// Put a site's chosen pair of traversals in phase order.
    ///
    /// The pair is swapped if the site's phase names the same two traversals in the other order.
    /// Otherwise it comes back in the order given.
    vector<int> phase_ordered_genotype(size_t record_key, const vector<int>& genotype) const;

    /// Return which parent strand, 0 or 1, carries a site with a one-allele genotype.
    ///
    /// This is the phase's `nested_strand`, which the VCF writes as `a|.` or `.|a`. It is 0
    /// when the site has no such phase or its phase names a different allele.
    int haploid_slot(size_t record_key, const vector<int>& genotype) const;

private:
    /// The phase calls in the order they were made. Where a site has more than one, the last one
    /// is its phase.
    vector<PhaseCall> phase_calls;
    /// Each site's phase, keyed by record key, as `freeze_for_render()` last copied it.
    ///
    /// It is keyed by record and not by position, because a record's POS depends on which
    /// alleles it carries.
    unordered_map<size_t, PhaseCall> render_phases;
};

/**
 * Find the phased GT that `phase` gives a record.
 *
 * The record's genotype is `site_genotype`, in VCF allele numbers, and `trav_to_allele` gives the
 * VCF allele of each traversal in it. A one-strand site of a nested chain gets its allele on its
 * strand and "." on the other, and any other haploid site gets its allele alone.
 *
 * @return False, leaving `gt` alone, if the alleles of `phase` are not a reordering of
 * `site_genotype`, since phasing must not change a genotype.
 */
bool phased_gt(const LinkageCollector::PhaseCall& phase, const vector<int>& site_genotype,
               const map<int, int>& trav_to_allele, string& gt);

}

#endif
