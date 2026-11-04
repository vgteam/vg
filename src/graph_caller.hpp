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
#include "snarls.hpp"
#include "traversal_finder.hpp"
#include "snarl_caller.hpp"
#include "region.hpp"
#include "zstdutil.hpp"
#include "vg/io/alignment_emitter.hpp"
#include "gref.hpp"
#include "../artifacts/ploidy_regions.hpp"
#include "../artifacts/site_scheduler.hpp"

namespace vg {



using namespace std;

using vg::io::AlignmentEmitter;

/**
 * GraphCaller: Use the snarl decomposition to call snarls in a graph
 */
class GraphCaller {
public:

    using RecurseType = vg::RecurseType;
    static constexpr RecurseType RecurseOnFail = RecurseType::OnFail;
    static constexpr RecurseType RecurseAlways = RecurseType::Always;
    static constexpr RecurseType RecurseNever = RecurseType::Never;
   
    GraphCaller(SnarlCaller& snarl_caller,
                SnarlManager& snarl_manager);

    virtual ~GraphCaller();

    /**
     * Run call_snarl() on every top-level snarl in the manager, in parallel, and then on children
     * as `recurse_type` says.
     *
     * RecurseAlways calls the children of every snarl, RecurseOnFail the children of each snarl
     * whose call returned false, and RecurseNever no children. Calling continues down the snarl
     * tree the same way.
     *
     * Each snarl is its own parallel job, unless set_snarl_batching(w) was called. Then each job
     * calls, in node-ID order, the top-level snarls whose lower boundary node ID is in one window
     * [k*w, (k+1)*w). Children are always one job each, started in node-ID order.
     */
    virtual void call_top_level_snarls(const HandleGraph& graph, RecurseType recurse_type = RecurseOnFail);

    /// For every chain, cut it up into pieces using max_edges and max_trivial to cap the size of each piece
    /// then make a fake snarl for each chain piece and call it.  If a fake snarl fails to call,
    /// It's child chains will be recursed on (if selected)_
    virtual void call_top_level_chains(const HandleGraph& graph,
                                       size_t max_edges,
                                       size_t max_trivial,
                                       RecurseType recurise_type = RecurseOnFail);

    /// Call a given snarl, and print the output to out_stream
    virtual bool call_snarl(const Snarl& snarl) = 0;

    /// toggle progress messages
    void set_show_progress(bool show_progress);

    /// Make call_top_level_snarls call top-level snarls in batches, by windows of `window_size`
    /// node IDs.
    ///
    /// This suits a SnarlCaller that loads its reads one window of node IDs at a time. A size of
    /// 0 gives one job per snarl.
    void set_snarl_batching(size_t window_size);

    /// Set the ploidy of the regions of the reference that `regions` covers, in place of the
    /// ploidy of their paths.
    void set_ploidy_regions(PloidyRegions regions);

protected:

    /// Break up a chain into bits that we want to call using size heuristics
    vector<Chain> break_chain(const HandleGraph& graph, const Chain& chain, size_t max_edges, size_t max_trivial) const;
    
protected:

    /// Our Genotyper
    SnarlCaller& snarl_caller;

    /// Our snarls
    SnarlManager& snarl_manager;

    /// The width, in node IDs, of the windows that top-level snarls are batched by, or 0 for no
    /// batching.
    size_t snarl_batch_window = 0;

    /// Toggle progress messages
    bool show_progress;

    /// The ploidy of regions of the reference, which overrides the ploidy of their paths.
    PloidyRegions ploidy_regions;
};

static void flip_snarl(Snarl& snarl) {
    Visit v = snarl.start();
    *snarl.mutable_start() = reverse(snarl.end());
    *snarl.mutable_end() = reverse(v);
}

/// Return the offset of a reference path, or 0 if `offsets` has none for it.
///
/// Unlike `operator[]`, this does not insert, so worker threads can share the map.
static inline size_t ref_offset_of(const map<string, size_t>& offsets, const string& path) {
    auto it = offsets.find(path);
    return it != offsets.end() ? it->second : 0;
}
/// Return the ploidy of a reference path, or 0 if `ploidies` has none for it.
///
/// Unlike `operator[]`, this does not insert, so worker threads can share the map.
static inline int ref_ploidy_of(const map<string, int>& ploidies, const string& path) {
    auto it = ploidies.find(path);
    return it != ploidies.end() ? it->second : 0;
}












}

#endif
