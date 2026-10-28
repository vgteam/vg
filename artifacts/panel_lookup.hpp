#ifndef VG_PANEL_LOOKUP_HPP_INCLUDED
#define VG_PANEL_LOOKUP_HPP_INCLUDED

#include <vector>

#include <gbwt/cached_gbwt.h>

#include "../src/handle.hpp"
#include "materialized_site.hpp"
#include "site_values.hpp"
#include <vg/vg.pb.h>

namespace vg {

using namespace std;

/**
 * A directory of which allele of a site each panel haplotype carries.
 *
 * It answers by asking the GBWT which haplotypes follow each of the site's traversals. This is
 * the GBWT that haplotype enumeration draws the candidate alleles from.
 *
 * Each thread gets its own cache of decompressed GBWT records, because adjacent sites share most
 * of their records and `CachedGBWT` has no locking. The constructor builds the caches, so that a
 * parallel region never allocates one.
 */
class PanelLookup {
public:
    /// Make a lookup with no panel, whose `alleles` returns an empty row.
    PanelLookup() = default;

    /**
     * Make a lookup over the panel in `gbwt`, which holds walks in `graph`.
     *
     * GBWT sequence `s` belongs to panel haplotype `sequence_to_haplotype[s]`. No pointer is
     * owned. A null GBWT or haplotype map means no panel.
     *
     * @param panel_size The number of panel haplotypes, or 0 to use the length of
     * `sequence_to_haplotype`.
     */
    PanelLookup(const HandleGraph* graph, const gbwt::GBWT* gbwt,
                const vector<size_t>* sequence_to_haplotype, size_t panel_size);

    /**
     * Return which of the alleles `candidates` each panel haplotype carries.
     *
     * The row has one entry per panel haplotype, holding an index into the alleles, or -1 where
     * the haplotype does not traverse the site. A haplotype carries an allele when it follows the
     * allele's walk. An allele that changes the sequences of its walk is carried by no haplotype,
     * because no graph walk spells it. The row is empty when there is no panel.
     */
    vector<int> alleles(const CandidateAlleles& candidates) const;

    /// Get the GBWT, or null when there is no panel.
    const gbwt::GBWT* gbwt() const { return index; }

    /// Get the panel haplotype of each GBWT sequence, or null when there is no panel.
    const vector<size_t>* sequence_to_haplotype() const { return haplotype_of_sequence; }

private:
    /// Clear the calling thread's cache if the site that `travs` cross is far from where the
    /// cache was started, and start it again at the site.
    void recenter_cache(const vector<Traversal>& travs) const;

    /// Return the GBWT sequences that follow all of `walk`, which must not be empty.
    vector<gbwt::size_type> sequences_following(const Traversal& walk) const;

    const HandleGraph* graph = nullptr;
    const gbwt::GBWT* index = nullptr;
    const vector<size_t>* haplotype_of_sequence = nullptr;

    /// The number of panel haplotypes, which is the length of each row from `alleles`.
    ///
    /// Under a gRef cover the GBWT has more sequences than the panel has haplotypes.
    size_t panel_size = 0;

    /// The record cache of each thread, by OpenMP thread number.
    mutable vector<gbwt::CachedGBWT> cache;

    /// The node ID at which each thread's cache was started, or 0 if it was not started yet.
    ///
    /// `CachedGBWT` only grows, so a thread's cache is cleared and started again when its site
    /// is more than a fixed span of node IDs from here.
    mutable vector<nid_t> cache_origin;
};

}

#endif
