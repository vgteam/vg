#include <omp.h>

#include "panel_lookup.hpp"

namespace vg {

PanelLookup::PanelLookup(const HandleGraph* graph, const gbwt::GBWT* gbwt,
                         const vector<size_t>* sequence_to_haplotype, size_t panel_size)
    : graph(graph), index(gbwt), haplotype_of_sequence(sequence_to_haplotype),
      panel_size(panel_size) {
    if (gbwt != nullptr) {
        cache.reserve(omp_get_max_threads());
        for (int i = 0; i < omp_get_max_threads(); ++i) {
            cache.emplace_back(*gbwt);
        }
        cache_origin.assign(omp_get_max_threads(), 0);
    }
}

void PanelLookup::recenter_cache(const vector<Traversal>& travs) const {
    // Windows are in node ID order, so a thread does not come back to an earlier one. Clearing
    // the cache once the site is more than about a fetch window from its origin keeps the cache
    // to about one window, while adjacent snarls still share records.
    static const nid_t CACHE_ANCHOR_SPAN = 4096;
    const int thread = omp_get_thread_num();
    if ((size_t)thread >= cache.size() || (size_t)thread >= cache_origin.size() || travs.empty()) {
        return;
    }
    const nid_t lead = travs[0].empty() ? 0 : graph->get_id(travs[0].front());
    if (lead == 0) {
        return;
    }
    nid_t& anchor = cache_origin[thread];
    if (anchor == 0 || lead > anchor + CACHE_ANCHOR_SPAN || lead + CACHE_ANCHOR_SPAN < anchor) {
        cache[thread].clearCache();
        anchor = lead;
    }
}

vector<gbwt::size_type> PanelLookup::sequences_following(const Traversal& walk) const {
    // Search the thread's cache when it has one, which gives the same answers as the index but
    // keeps records decompressed between sites. A lookup built without a GBWT has no caches.
    const int thread = omp_get_thread_num();
    const bool cached = (size_t)thread < cache.size();
    gbwt::SearchState state;
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
            return {};
        }
    }
    return cached ? cache[thread].locate(state) : index->locate(state);
}

vector<int> PanelLookup::alleles(const CandidateAlleles& candidates) const {
    const vector<Traversal>& travs = candidates.walks;
    vector<int> out;
    if (index == nullptr || haplotype_of_sequence == nullptr) {
        return out;
    }
    // A haplotype whose path ends inside the site carries no allele there, which is different
    // from carrying the reference.
    const size_t row = panel_size > 0 ? panel_size : haplotype_of_sequence->size();
    out.assign(row, -1);

    recenter_cache(travs);

    for (size_t a = 0; a < travs.size(); ++a) {
        const Traversal& walk = travs[a];
        if (walk.empty() || !candidates.spelled_by_walk(a)) {
            continue;
        }
        for (gbwt::size_type seq : sequences_following(walk)) {
            if (seq >= haplotype_of_sequence->size()) {
                continue;
            }
            size_t hap = (*haplotype_of_sequence)[seq];
            if (hap < out.size()) {
                // A haplotype stored as several fragments can traverse one site twice, by two
                // traversals. The last one written wins.
                out[hap] = (int)a;
            }
        }
    }
    return out;
}

}
