#ifndef VG_SYMBOLIC_ALLELE_HPP_INCLUDED
#define VG_SYMBOLIC_ALLELE_HPP_INCLUDED

/**
 * \file symbolic_allele.hpp
 *
 * Symbolic alleles: a walk through a site with each pass through a child chain replaced by one
 * symbol for that chain, rather than the path taken through it.
 *
 * A walk runs through every interior node, so two walks that differ only inside a child chain are
 * different walks. Their symbolic forms are equal, since they take the same route at this level of
 * the snarl tree. A walk whose symbolic form equals the reference walk's is therefore the
 * reference allele at this site, and its differences are reported by the child chain's own
 * records. A walk that skips a chain, or passes through different ones, differs at this level, so
 * a deletion of a chain is still a deletion here.
 */

#include <functional>
#include <utility>
#include <ostream>
#include <vector>

#include "handle.hpp"
#include "site_values.hpp"

namespace vg {

using namespace std;

/**
 * One step of a symbolic allele: either a plain node, or a whole child chain as one symbol.
 *
 * A chain is identified by its own boundary nodes, not by those of the child snarl a traversal
 * enters first, so two traversals that enter and leave a chain through the same boundaries carry
 * the same symbol, however they cross it.
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
 * Project a walk through a site into symbolic form, given the site's `children`.
 *
 * Follows the walk and, wherever it enters a child site, emits one symbol for that site's chain
 * and resumes at the handle that leaves the chain.
 *
 * A chain entered but not left within the walk, as in a malformed or cyclic walk, is emitted as a
 * plain node step, so that the rest of the walk is not lost; losing it could make different
 * alleles compare equal.
 *
 * `out_visit_ranges`, when given, reports for each emitted step the half-open range of the walk's
 * handles it covers. The ranges partition [0, walk.size()) in order, so a step's sequence is the
 * concatenation of its handles'. A chain symbol's range is [entry, exit): the exit bound belongs
 * to the next step, since the chain shares it with its successor.
 */
SymbolicAllele symbolic_allele(const HandleGraph& graph, const Traversal& walk,
                               const SiteChildren& children,
                               vector<pair<int, int>>* out_visit_ranges = nullptr);

/// True if the two walks take the same route through the site at this level of the hierarchy.
bool symbolically_equal(const HandleGraph& graph, const Traversal& a, const Traversal& b,
                        const SiteChildren& children);

/**
 * One difference between two symbolic alleles: a half-open step range on each side.
 *
 * Only differences are reported; the steps between one block's end and the next block's start
 * match on both sides, so an empty result means the two alleles take the same route. Either range
 * may be empty: an empty ref range is an insertion, an empty alt range a deletion.
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
 * Align two symbolic alleles and return the blocks where they differ, in reference order.
 *
 * The alignment minimises edit distance with substitutions at cost 1, so that a replacement is
 * one block: with insertions and deletions only, [a,b] against [b,b] has a second minimal
 * alignment that gives two blocks around a spurious match. Remaining ties are broken the same way
 * every time, preferring a match, then a substitution, then a deletion, then an insertion, since
 * the result decides how many records a site writes.
 *
 * The dynamic program runs inside Ukkonen's band, in O((|ref| + |alt|) x D) time and
 * O(|ref| x D) space, where D is the edit distance, so its cost grows with the difference between
 * the two alleles rather than with their length.
 *
 * `out_alt_before_ref`, when given, is filled with |ref| + 1 entries: entry i is the number of alt
 * steps consumed before reference step i, not counting any inserted at boundary i. It turns a
 * reference step range into the alt step range aligned with it, as needed to write two
 * haplotypes' alleles over one reference span.
 */
vector<DiffBlock> symbolic_diff(const SymbolicAllele& ref, const SymbolicAllele& alt,
                                vector<int>* out_alt_before_ref = nullptr);

/// Print a symbolic allele or a difference block in a readable form.
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
