#ifndef VG_SYMBOLIC_ALLELE_HPP_INCLUDED
#define VG_SYMBOLIC_ALLELE_HPP_INCLUDED

/**
 * \file symbolic_allele.hpp
 *
 * A traversal of a snarl, with each excursion through a nested chain replaced by a symbol for that
 * chain rather than the concrete path taken through it.
 *
 * A SnarlTraversal runs from snarl start to snarl end through every interior node, so nested
 * variation is baked into it: two traversals differing only *inside* a child chain are two distinct
 * traversals, and the caller emits them as two long, nearly identical alleles. Measured on HG002,
 * 55,222 of vg's 142,707 autosomal SNV false negatives sit inside a large allele vg itself emitted,
 * against a 0.6% rate among variants it calls correctly.
 *
 * Comparing symbolic forms instead answers a different question: not "is this the same sequence"
 * but "is this the same route at this level of the hierarchy". A traversal whose symbolic form
 * equals the reference traversal's differs from the reference only inside child chains, so it is
 * the reference allele *here* and its differences belong to those chains' own records. A traversal
 * that skips a chain, or crosses different ones, stays a genuine allele at this level -- so a real
 * deletion is still reported as a deletion.
 *
 * The whole design and the measurements behind it are in the companion evaluation repository, as
 * docs/nested-calling-design.md.
 */

#include <functional>
#include <utility>
#include <ostream>
#include <vector>

#include "handle.hpp"
#include "snarls.hpp"
#include <vg/vg.pb.h>

namespace vg {

using namespace std;

/**
 * One step of a symbolic allele: either a plain node, or a whole child chain collapsed to a symbol.
 *
 * A chain is identified by the boundary nodes of the chain itself, not of the child snarl the
 * traversal happened to enter first. Two traversals that enter a chain through the same boundary
 * and leave through the same boundary carry the same symbol however they cross it, which is the
 * entire point.
 */
struct SymbolicStep {
    /// Node id for a plain step; for a chain symbol, the chain's start node.
    nid_t id = 0;
    /// For a chain symbol, the chain's end node. 0 for a plain node step.
    nid_t end_id = 0;
    /// Orientation of the step as traversed.
    bool backward = false;

    bool is_chain() const { return end_id != 0; }

    bool operator==(const SymbolicStep& o) const {
        return id == o.id && end_id == o.end_id && backward == o.backward;
    }
    bool operator!=(const SymbolicStep& o) const { return !(*this == o); }
};

/// A traversal with child chains collapsed to symbols.
using SymbolicAllele = vector<SymbolicStep>;

/**
 * Project a traversal of `site` into symbolic form.
 *
 * Walks the traversal and, wherever a visit enters a child chain of `site`, emits one symbol for
 * that chain and resumes at the visit that leaves it. Visits that already carry a Snarl rather than
 * a node -- the protobuf supports both, and NestedFlowCaller emitted them -- are taken as symbols
 * directly.
 *
 * A chain entered but never left within this traversal (which a malformed or cyclic traversal can
 * produce) is emitted as a plain node step rather than swallowing the rest of the traversal: losing
 * the tail would silently make unrelated alleles compare equal, which is the one error mode this
 * must not have.
 */
/// Optionally reports, for each emitted step, the half-open range of `trav` visits it covers.
/// The ranges partition [0, visit_size) contiguously and in order, so a step range can be turned
/// straight into a sequence by concatenating those visits. Note a chain symbol's range is
/// [entry, exit): the exit boundary node belongs to whatever step comes next, because it is shared
/// between the chain and its successor.
SymbolicAllele symbolic_allele(const SnarlTraversal& trav, const Snarl& site,
                               const SnarlManager& snarl_manager,
                               vector<pair<int, int>>* out_visit_ranges = nullptr);

/// Whether `site` resolves to the snarl the manager knows it as, which is the precondition for
/// recognising any child chain at all. False means projection degenerates to the plain node list
/// with no symbols -- the case `flip_snarl` produces for a snarl whose reference path runs
/// backwards, where symbolic collapsing is therefore inert.
/// `out_reversed`, when given, reports whether the site resolved only through the REVERSED boundary
/// pairing -- i.e. whether this is one of the snarls `flip_snarl` reverses because the reference path
/// runs backwards through it.
///
/// Deliberately an out-parameter rather than a counter inside the resolver. The resolver is called
/// once per projection and projection runs per traversal, so counting there would count calls, not
/// sites, at several times the rate of the per-site counter it has to be compared against.
bool symbolic_site_resolvable(const Snarl& site, const SnarlManager& snarl_manager,
                              bool* out_reversed = nullptr);

/// The boundary node pair of the chain `child` belongs to, which is the identity a chain symbol
/// carries. Exposed so a caller holding a child snarl can find the symbol that child collapses to
/// without re-deriving the chain itself.
pair<nid_t, nid_t> chain_bounds_of(const Snarl* child, const SnarlManager& snarl_manager);

/// True if the two traversals are the same route through `site` at this level of the hierarchy.
bool symbolically_equal(const SnarlTraversal& a, const SnarlTraversal& b, const Snarl& site,
                        const SnarlManager& snarl_manager);

/// Whether any allele's traversal crosses a child chain, i.e. whether `site` is a non-leaf snarl
/// with something to symbolise. Leaf snarls are unaffected by any of this and on chr20 are 98.5%
/// of emitted records.
bool has_child_chain(const SnarlTraversal& trav, const Snarl& site,
                     const SnarlManager& snarl_manager);

/**
 * One difference between two symbolic alleles: a half-open step range on each side.
 *
 * Only differences are reported. The matched runs between them are implicit -- the gap between one
 * block's end and the next block's start is matched on both sides -- so an empty result means the
 * two alleles are the same route.
 *
 * Either range may be empty: an empty ref range is a pure insertion, an empty alt range a pure
 * deletion.
 */
struct DiffBlock {
    int ref_begin = 0;
    int ref_end = 0;
    int alt_begin = 0;
    int alt_end = 0;

    bool ref_empty() const { return ref_end == ref_begin; }
    bool alt_empty() const { return alt_end == alt_begin; }

    bool operator==(const DiffBlock& o) const {
        return ref_begin == o.ref_begin && ref_end == o.ref_end &&
               alt_begin == o.alt_begin && alt_end == o.alt_end;
    }
};

/**
 * Align two symbolic alleles and return the difference blocks between them, in reference order.
 *
 * The cost model is edit distance **with substitution at cost 1**, not the insert/delete-only model
 * a plain `diff` uses. That is a deliberate disambiguation rather than a different notion of
 * distance: under insert/delete-only, [a,b] against [b,b] has two minimal alignments of equal cost,
 * one giving a single replacement and one giving two replacements separated by a spurious match, and
 * nothing in "minimum edits" chooses between them. Substitution at 1 beats delete-plus-insert at 2,
 * so the single-block reading wins strictly. This encodes "prefer fewer, larger blocks", which is
 * the same preference the block aggregation already expresses.
 *
 * Ties remain, and they are broken **deterministically** by preferring, at each traceback step, the
 * diagonal (match, then substitute) over deletion over insertion. Determinism is the load-bearing
 * property: an unstable tie-break makes output depend on nothing the caller controls, and this
 * function's result decides how many records a snarl emits.
 *
 * The DP is O(|ref| x |alt|) in time and space. A traversal pair too large for that degrades to one
 * block spanning both alleles entirely -- which is exactly the whole-allele behaviour that predates
 * this function, so the caller stays correct rather than merely surviving. `out_degraded`, when
 * given, is set true in that case and only that case, so the population can be counted instead of
 * assumed to be empty.
 */
/// `out_alt_before_ref`, when given, is filled with |ref| + 1 entries: entry i is the number of alt
/// steps consumed strictly before reference step i, counting nothing inserted at boundary i. It is
/// what turns a reference step range into the alt step range aligned to it, which the diploid join
/// needs in order to express two haplotypes' alleles over one shared reference span.
vector<DiffBlock> symbolic_diff(const SymbolicAllele& ref, const SymbolicAllele& alt,
                                bool* out_degraded = nullptr,
                                vector<int>* out_alt_before_ref = nullptr);

/// For logging and tests.
ostream& operator<<(ostream& out, const SymbolicAllele& allele);
ostream& operator<<(ostream& out, const DiffBlock& block);

}

namespace std {
template<> struct hash<vg::SymbolicAllele> {
    size_t operator()(const vg::SymbolicAllele& a) const {
        size_t h = 1469598103934665603ULL;
        for (const auto& s : a) {
            for (uint64_t part : {(uint64_t)s.id, (uint64_t)s.end_id, (uint64_t)s.backward}) {
                h ^= part;
                h *= 1099511628211ULL;
            }
        }
        return h;
    }
};
}

#endif
