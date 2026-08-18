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
 * See doc/nested-calling-design.md for the whole design and the measurements behind it.
 */

#include <functional>
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
SymbolicAllele symbolic_allele(const SnarlTraversal& trav, const Snarl& site,
                               const SnarlManager& snarl_manager);

/// True if the two traversals are the same route through `site` at this level of the hierarchy.
bool symbolically_equal(const SnarlTraversal& a, const SnarlTraversal& b, const Snarl& site,
                        const SnarlManager& snarl_manager);

/// Whether any allele's traversal crosses a child chain, i.e. whether `site` is a non-leaf snarl
/// with something to symbolise. Leaf snarls are unaffected by any of this and on chr20 are 98.5%
/// of emitted records.
bool has_child_chain(const SnarlTraversal& trav, const Snarl& site,
                     const SnarlManager& snarl_manager);

/// For logging and tests.
ostream& operator<<(ostream& out, const SymbolicAllele& allele);

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
