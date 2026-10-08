#ifndef VG_STAGED_SITE_HPP_INCLUDED
#define VG_STAGED_SITE_HPP_INCLUDED

#include <memory>
#include <string>
#include <vector>

#include "snarl_caller.hpp"
#include "snarls.hpp"

namespace vg {

using namespace std;

/// A staged site: one site's genotyping result from the direct pass, kept until the render
/// builds the site's records.
///
/// A nested chain's ploidy depends on its parent's chosen genotype, which is known only after
/// the direct pass, so the result is kept rather than computed again. The `CallInfo` has the
/// answer at both ploidies (see `alt_ploidy_info`), so the record can be rendered at whichever
/// ploidy the linkage pass gives the site. `snarl` is held by value because
/// `call_snarl_internal` may work on a flipped copy.
struct PendingRecord {
    Snarl snarl;
    string ref_path_name;
    int ref_offset = 0;
    vector<SnarlTraversal> travs;
    int ref_trav_idx = -1;
    /// The direct pass's genotype, before the linkage model, and its ploidy. A nested chain's
    /// ploidy here is the number of the parent's called alleles that cross it, or the parent's
    /// ploidy when none does; `call_info` also holds the answer at the other ploidy.
    vector<int> genotype;
    int ploidy = 2;
    unique_ptr<SnarlCaller::CallInfo> call_info;
    size_t record_key = 0;
    size_t parent_record_key = 0;
    /// See NestedContext::chain_key.
    size_t chain_key = 0;
    /// See NestedContext::reported_inline. Its line is held back, as for `no_reference`. The
    /// direct pass tests the parent's direct call; under the linkage model the linkage pass
    /// tests again with the parent's chosen genotype, which the parent's blocks are built from.
    bool reported_inline = false;
    /// This snarl has no reference path, so no line can be written for it, since REF and POS are
    /// undefined. It is still genotyped and recorded in the linkage model.
    bool no_reference = false;
    /// For a snarl with no reference path: its parent's reference start plus `chain_offset`,
    /// standing in for the position it lacks.
    int64_t position_from_parent = 0;
    /// `NestedContext::parent_offset` for this chain. The direct pass takes it from the parent's
    /// direct call, and the linkage pass computes it again from the parent's chosen genotype, so
    /// that an off-reference chain is placed along an allele of its parent's chosen genotype.
    size_t chain_offset = 0;
    /// See NestedContext::parent_crossing.
    uint64_t parent_crossing = 0;
    /// False when the parent has more than 64 candidate traversals, too many for
    /// `parent_crossing`. A 0 mask then means unknown, and the linkage pass leaves the chain at the
    /// ploidy the direct pass gave it. The linkage pass computes the mask again when it revises
    /// or first records the parent.
    bool crossing_known = true;
    /// The site's level.
    uint8_t level = 0;
    /// Set when the parent's chosen genotype, or an ancestor's, does not carry this chain, so
    /// the chain and its descendants do not exist in the sample and are not revised or
    /// written. Each linkage pass decides it again, so a chain dropped in one pass can come
    /// back in the next.
    bool dropped = false;
    /// `panel_lookup.alleles(travs)`, computed once: the traversals do not change after the
    /// direct pass, and each re-genotyping round would otherwise repeat the GBWT lookups.
    vector<int> panel_cache;
    bool panel_cached = false;

};

}

#endif
