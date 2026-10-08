#include <omp.h>

#include "vcf_output_caller.hpp"

namespace vg {

vector<int> VCFOutputCaller::panel_alleles(const HandleGraph& graph,
                                          const vector<SnarlTraversal>& travs) const {
    vector<int> out;
    if (linkage_gbwt == nullptr || linkage_sequence_to_haplotype == nullptr) {
        return out;
    }
    // -1 means the haplotype carries no allele here, which is different from carrying the
    // reference: a haplotype whose path ends inside the site has nothing to say. Sized by the
    // panel, since the row is indexed by haplotype.
    const size_t row = linkage_panel_size > 0 ? linkage_panel_size
                                              : linkage_sequence_to_haplotype->size();
    out.assign(row, -1);

    // The cache, not the index: same results, but records stay decompressed between sites.
    // Falls back to the index itself if set_linkage was never given one to size the vector.
    int thread = omp_get_thread_num();
    const bool cached = (size_t)thread < linkage_gbwt_cache.size();

    // CachedGBWT only grows, and with node-ID-ordered windows a thread does not come back to an
    // earlier window, so the cache is cleared when the site moves more than a fetch window past
    // where it was filled. Adjacent snarls still share records, and the cache stays to about one
    // window.
    if (cached && (size_t)thread < linkage_gbwt_cache_origin.size() && !travs.empty()) {
        static const nid_t CACHE_ANCHOR_SPAN = 4096;
        nid_t lead = 0;
        for (int64_t i = 0; i < travs[0].visit_size() && lead == 0; ++i) {
            lead = travs[0].visit(i).node_id();
        }
        if (lead != 0) {
            nid_t& anchor = linkage_gbwt_cache_origin[thread];
            if (anchor == 0 || lead > anchor + CACHE_ANCHOR_SPAN
                || lead + CACHE_ANCHOR_SPAN < anchor) {
                linkage_gbwt_cache[thread].clearCache();
                anchor = lead;
            }
        }
    }

    for (size_t a = 0; a < travs.size(); ++a) {
        const SnarlTraversal& trav = travs[a];
        if (trav.visit_size() < 1) {
            continue;
        }
        gbwt::SearchState state;
        bool ok = true;
        for (int64_t i = 0; i < trav.visit_size(); ++i) {
            const Visit& visit = trav.visit(i);
            if (visit.node_id() == 0) {
                // A visit to a child snarl rather than a node: the traversal is not expanded, so
                // it cannot be looked up in the GBWT.
                ok = false;
                break;
            }
            gbwt::node_type node = gbwt::Node::encode(visit.node_id(), visit.backward());
            if (cached) {
                const gbwt::CachedGBWT& c = linkage_gbwt_cache[thread];
                state = (i == 0) ? c.find(node) : c.extend(state, node);
            } else {
                state = (i == 0) ? linkage_gbwt->find(node) : linkage_gbwt->extend(state, node);
            }
            if (state.empty()) {
                ok = false;
                break;
            }
        }
        if (!ok || state.empty()) {
            continue;
        }
        vector<gbwt::size_type> seqs = cached ? linkage_gbwt_cache[thread].locate(state)
                                              : linkage_gbwt->locate(state);
        for (gbwt::size_type seq : seqs) {
            if (seq < linkage_sequence_to_haplotype->size()) {
                size_t hap = (*linkage_sequence_to_haplotype)[seq];
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
