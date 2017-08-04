#include "realigner.hpp"

namespace vg {


void Realigner::construct(void) {

    if (mapper) delete mapper;
    if (xgidx) delete xgidx;
    if (gcsaidx) delete gcsaidx;
    if (lcpidx) delete lcpidx;

    if (debug) cerr << "building xg index" << endl;
    xgidx = new xg::XG(graph->graph);
    if (debug) cerr << "building GCSA2 index" << endl;
    if (edge_max) {
        VG gcsa_graph = *graph; // copy the graph
        // remove complex components
        gcsa_graph.prune_complex_with_head_tail(idx_kmer_size, edge_max);
        if (subgraph_prune) gcsa_graph.prune_short_subgraphs(subgraph_prune);
        // then index
        gcsa_graph.build_gcsa_lcp(gcsaidx, lcpidx, idx_kmer_size, idx_path_only, false, doubling_steps);
    } else {
        // if no complexity reduction is requested, just build the index
        graph->build_gcsa_lcp(gcsaidx, lcpidx, idx_kmer_size, idx_path_only, false, doubling_steps);
    }
    mapper = new Mapper(xgidx, gcsaidx, lcpidx);
    /*
    // todo
    { // set mapper variables
        mapper->debug = debug_align;
        mapper->context_depth = context_depth;
        mapper->thread_extension = thread_extension;
        mapper->max_attempts = max_attempts;
        mapper->min_identity = min_identity;
        mapper->alignment_threads = alignment_threads;
        mapper->max_mem_length = max_mem_length;
        mapper->min_mem_length = min_mem_length;
        mapper->hit_max = hit_max;
        mapper->max_target_factor = max_target_factor;
        mapper->max_multimaps = max_multimaps;
        mapper->match = match;
        mapper->mismatch = mismatch;
        mapper->gap_open = gap_open;
        mapper->gap_extend = gap_extend;
    }
    */
}

Alignment Realigner::realign(const Alignment& aln) {
    assert(mapper != nullptr);
    double ident = identity(aln.path());
    int softclip = softclip_start(aln) + softclip_end(aln);
    if (ident < identity_trigger
        || softclip >= softclip_trigger) {
        auto alns = mapper[tid]->align_multi(aln);
        auto& raln = alns.front();
        double rident = identity(raln.path());
        int rsoftclip = softclip_start(raln) + softclip_end(raln);
        if (rident > ident || rsoftclip < softclip) {
            return raln;
        } else {
            return aln;
        }
    } else {
        return aln;
    }
}


}
