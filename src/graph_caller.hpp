#ifndef VG_GRAPH_CALLER_HPP_INCLUDED
#define VG_GRAPH_CALLER_HPP_INCLUDED

#include <atomic>
#include <iostream>
#include <algorithm>
#include <array>
#include <functional>
#include <cmath>
#include <limits>
#include <unordered_set>
#include <tuple>
#include <gbwt/cached_gbwt.h>
#include "handle.hpp"
#include "linkage_model.hpp"
#include "snarls.hpp"
#include "traversal_finder.hpp"
#include "anchor.hpp"
#include "read_phasing.hpp"
#include "regenotype.hpp"
#include "snarl_caller.hpp"
#include "symbolic_allele.hpp"
#include "region.hpp"
#include "zstdutil.hpp"
#include "vg/io/alignment_emitter.hpp"
#include "gref.hpp"
#include "vcf_genotype_likelihoods.hpp"

namespace vg {



using namespace std;

using vg::io::AlignmentEmitter;

/**
 * GraphCaller: Use the snarl decomposition to call snarls in a graph
 */
class GraphCaller {
public:

    enum RecurseType { RecurseOnFail, RecurseAlways, RecurseNever };
   
    GraphCaller(SnarlCaller& snarl_caller,
                SnarlManager& snarl_manager);

    virtual ~GraphCaller();

    /// Run call_snarl() on every top-level snarl in the manager, in parallel. Then call it on
    /// children, as `recurse_type` says: those of every snarl (RecurseAlways), those of snarls
    /// whose call returned false (RecurseOnFail), or none (RecurseNever), and so on down.
    ///
    /// By default each snarl is its own parallel job. After set_snarl_batching(w), the top-level
    /// snarls are grouped into batches instead: a batch holds the snarls whose lower boundary
    /// node ID falls in one window of w IDs, [k*w, (k+1)*w) for some k, each batch is one job,
    /// and a job calls its snarls in node-ID order. Children are still one job each, started in
    /// node-ID order.
    virtual void call_top_level_snarls(const HandleGraph& graph, RecurseType recurse_type = RecurseOnFail);

    /// For every chain, cut it up into pieces using max_edges and max_trivial to cap the size of
    /// each piece then make a fake snarl for each chain piece and call it.  If a fake snarl fails
    /// to call, It's child chains will be recursed on (if selected)_
    virtual void call_top_level_chains(const HandleGraph& graph,
                                       size_t max_edges,
                                       size_t max_trivial,
                                       RecurseType recurise_type = RecurseOnFail);

    /// Call a given snarl, and print the output to out_stream
    virtual bool call_snarl(const Snarl& snarl) = 0;

    /// toggle progress messages
    void set_show_progress(bool show_progress);

    /// Batch call_top_level_snarls' parallel jobs by windows of `window_size` node IDs (see
    /// call_top_level_snarls), so that snarls with nearby node IDs are called together, which
    /// suits a SnarlCaller that loads its reads a window of node IDs at a time. 0, the default,
    /// gives one job per snarl.
    void set_snarl_batching(size_t window_size);

protected:

    /// Break up a chain into bits that we want to call using size heuristics
    vector<Chain> break_chain(const HandleGraph& graph, const Chain& chain, size_t max_edges, size_t max_trivial);
    
protected:

    /// Our Genotyper
    SnarlCaller& snarl_caller;

    /// Our snarls
    SnarlManager& snarl_manager;

    /// See set_snarl_batching.
    size_t snarl_batch_window = 0;

    /// Toggle progress messages
    bool show_progress;
};

static void flip_snarl(Snarl& snarl) {
    Visit v = snarl.start();
    *snarl.mutable_start() = reverse(snarl.end());
    *snarl.mutable_end() = reverse(v);
}

/// The 1-based position on the base path of the base `along_path` bases into `ref_path_name`, a
/// path that may name a subrange of its base path.
static int64_t base_path_position(const string& ref_path_name, int64_t along_path) {
    subrange_t subrange;
    Paths::strip_subrange(ref_path_name, &subrange);
    const int64_t basepath_offset = subrange == PathMetadata::NO_SUBRANGE ? 0 : (int64_t)subrange.first;
    return along_path + 1 + basepath_offset;
}

/// Look a reference path up in one of the per-caller maps without inserting on a miss:
/// `operator[]` inserts, and these maps are read from worker threads.
static inline size_t ref_offset_of(const map<string, size_t>& offsets, const string& path) {
    auto it = offsets.find(path);
    return it != offsets.end() ? it->second : 0;
}
static inline int ref_ploidy_of(const map<string, int>& ploidies, const string& path) {
    auto it = ploidies.find(path);
    return it != ploidies.end() ? it->second : 0;
}












}

#endif
