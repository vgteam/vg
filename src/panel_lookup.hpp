#ifndef VG_PANEL_LOOKUP_HPP_INCLUDED
#define VG_PANEL_LOOKUP_HPP_INCLUDED

#include <vector>

#include <gbwt/cached_gbwt.h>

#include "handle.hpp"
#include "site_values.hpp"
#include <vg/vg.pb.h>

namespace vg {

using namespace std;

/**
 * Which allele of a site each panel haplotype carries, found by asking the GBWT which haplotypes
 * take each of the site's traversals. It is the GBWT that haplotype enumeration draws the
 * candidate alleles from.
 *
 * Keeps one cache of decompressed GBWT records per thread, since adjacent sites share most of
 * their records, and per thread because `CachedGBWT` has no locking. The constructor builds the
 * caches, so that nothing allocates one in a parallel region.
 */
class PanelLookup {
public:
    /// No panel, so that `alleles` returns an empty row.
    PanelLookup() = default;

    /// The panel in `gbwt`, whose sequence `s` belongs to panel haplotype
    /// `sequence_to_haplotype[s]`, over walks in `graph`. `panel_size` is the number of
    /// haplotypes, or 0 to size each row by the length of `sequence_to_haplotype`. No pointer is
    /// owned, and a null GBWT or haplotype map means no panel.
    PanelLookup(const HandleGraph* graph, const gbwt::GBWT* gbwt,
                const vector<size_t>* sequence_to_haplotype, size_t panel_size);

    /// Which of the walks `travs` each panel haplotype follows, or -1 where it does not traverse
    /// the site. Empty when there is no panel.
    vector<int> alleles(const vector<Traversal>& travs) const;

    /// The GBWT, or null when there is no panel.
    const gbwt::GBWT* gbwt() const { return index; }

    /// The panel haplotype of each GBWT sequence, or null when there is no panel.
    const vector<size_t>* sequence_to_haplotype() const { return haplotype_of_sequence; }

private:
    const HandleGraph* graph = nullptr;
    const gbwt::GBWT* index = nullptr;
    const vector<size_t>* haplotype_of_sequence = nullptr;

    /// Panel size, so that `alleles` sizes its row by the haplotypes rather than by the GBWT's
    /// sequence count, which is larger under a gRef cover.
    size_t panel_size = 0;

    /// One cache per thread; see the class comment.
    mutable vector<gbwt::CachedGBWT> cache;

    /// The node ID at which each thread's cache was started. `CachedGBWT` only grows, so when a
    /// thread's site is more than a fixed span of node IDs from this, its cache is cleared and
    /// started again there.
    mutable vector<nid_t> cache_origin;
};

}

#endif
