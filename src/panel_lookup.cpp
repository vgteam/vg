#include <omp.h>

#include "panel_lookup.hpp"

namespace vg {

PanelLookup::PanelLookup(const HandleGraph* graph, const gbwt::GBWT* gbwt,
                         const vector<size_t>* sequence_to_haplotype, size_t panel_size)
    : graph(graph), index(gbwt), haplotype_of_sequence(sequence_to_haplotype),
      panel_size(panel_size) {
    if (gbwt != nullptr) {
        // One per thread, built here so the parallel region never allocates one.
        cache.reserve(omp_get_max_threads());
        for (int i = 0; i < omp_get_max_threads(); ++i) {
            cache.emplace_back(*gbwt);
        }
        cache_origin.assign(omp_get_max_threads(), 0);
    }
}

vector<int> PanelLookup::alleles(const vector<Traversal>& travs) const {
    vector<int> out;
    if (index == nullptr || haplotype_of_sequence == nullptr) {
        return out;
    }
    // -1 means the haplotype carries no allele here, which is different from carrying the
    // reference: a haplotype whose path ends inside the site has nothing to say. Sized by the
    // panel, since the row is indexed by haplotype.
    const size_t row = panel_size > 0 ? panel_size : haplotype_of_sequence->size();
    out.assign(row, -1);

    // The cache, not the index: same results, but records stay decompressed between sites.
    // Falls back to the index itself if the constructor was never given one to size the vector.
    int thread = omp_get_thread_num();
    const bool cached = (size_t)thread < cache.size();

    // CachedGBWT only grows, and with node-ID-ordered windows a thread does not come back to an
    // earlier window, so the cache is cleared when the site moves more than a fetch window past
    // where it was filled. Adjacent snarls still share records, and the cache stays to about one
    // window.
    if (cached && (size_t)thread < cache_origin.size() && !travs.empty()) {
        static const nid_t CACHE_ANCHOR_SPAN = 4096;
        const nid_t lead = travs[0].empty() ? 0 : graph->get_id(travs[0].front());
        if (lead != 0) {
            nid_t& anchor = cache_origin[thread];
            if (anchor == 0 || lead > anchor + CACHE_ANCHOR_SPAN
                || lead + CACHE_ANCHOR_SPAN < anchor) {
                cache[thread].clearCache();
                anchor = lead;
            }
        }
    }

    for (size_t a = 0; a < travs.size(); ++a) {
        const Traversal& walk = travs[a];
        if (walk.empty()) {
            continue;
        }
        gbwt::SearchState state;
        bool ok = true;
        for (size_t i = 0; i < walk.size(); ++i) {
            gbwt::node_type node = gbwt::Node::encode(graph->get_id(walk[i]),
                                                      graph->get_is_reverse(walk[i]));
            if (cached) {
                const gbwt::CachedGBWT& c = cache[thread];
                state = (i == 0) ? c.find(node) : c.extend(state, node);
            } else {
                state = (i == 0) ? index->find(node) : index->extend(state, node);
            }
            if (state.empty()) {
                ok = false;
                break;
            }
        }
        if (!ok || state.empty()) {
            continue;
        }
        vector<gbwt::size_type> seqs = cached ? cache[thread].locate(state)
                                              : index->locate(state);
        for (gbwt::size_type seq : seqs) {
            if (seq < haplotype_of_sequence->size()) {
                size_t hap = (*haplotype_of_sequence)[seq];
                if (hap < out.size()) {
                    // A haplotype stored as several fragments could reach one site twice, with
                    // two traversals; the last one written wins.
                    out[hap] = (int)a;
                }
            }
        }
    }
    return out;
}

}
