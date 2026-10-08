#ifndef VG_BLOCK_RECORDS_HPP_INCLUDED
#define VG_BLOCK_RECORDS_HPP_INCLUDED

#include <atomic>
#include <vector>

#include "handle.hpp"
#include "snarls.hpp"

namespace vg {

using namespace std;

/// Counters for block emission. A member of each VCFOutputCaller, as MosaicCounters is.
struct AtomizeCounters {
    /// Sites that reached `tally_atomize`, so that the report can tell "nothing refused" from
    /// "never ran".
    std::atomic<size_t> sites{0};
    std::atomic<size_t> site_unresolvable{0};  // flip_snarl left projection with no symbols
    std::atomic<size_t> site_reversed{0};      // resolved only via the reversed pairing
    /// Sites written as blocks, and the lines they produced.
    std::atomic<size_t> split_sites{0}, split_lines{0};
    /// Chains whose own record is not written because a block's ALT already spells them.
    std::atomic<size_t> child_inlined{0};
    /// Why `emit_block_records` declined a site, by refusal point. Each means the site's single
    /// record is written instead.
    std::atomic<size_t> refuse[13] = {};
};

/// Counters for block emission: project the reference and each distinct called ALT traversal,
/// align them, and count. It changes no output.
void tally_atomize(const PathPositionHandleGraph& graph, const SnarlManager* mgr,
                   const Snarl& snarl, const vector<SnarlTraversal>& travs,
                   const vector<int>& genotype, int ref_trav_idx,
                   AtomizeCounters& atomize_counters);

}

#endif
