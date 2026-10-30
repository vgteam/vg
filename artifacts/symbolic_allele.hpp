#ifndef VG_SYMBOLIC_ALLELE_HPP_INCLUDED
#define VG_SYMBOLIC_ALLELE_HPP_INCLUDED

/**
 * \file symbolic_allele.hpp
 *
 * Symbolic alleles, which are walks through a site with each stretch inside a child chain replaced
 * by one symbol for that chain.
 *
 * Two walks that differ only inside a child chain have equal symbolic alleles, since they take the
 * same route at this level of the snarl tree. A walk whose symbolic allele equals the reference
 * walk's is the reference allele at this site, and the child chain's own records report its
 * differences. A walk that skips a chain, or crosses different ones, differs at this level, so a
 * deletion of a chain is still a deletion here.
 */

#include <functional>
#include <utility>
#include <ostream>
#include <vector>

#include "../src/handle.hpp"
#include "site_values.hpp"

namespace vg {

using namespace std;

/**
 * One step of a symbolic allele, which is a node or a whole child chain.
 *
 * A chain symbol names the chain by its own boundary nodes, so two walks that enter and leave a
 * chain through the same boundaries carry the same symbol, however they cross it.
 */
struct SymbolicStep {
    /// The node ID of a node step, or the start node of a chain symbol.
    nid_t id = 0;
    /// The end node of a chain symbol, or 0 for a node step.
    nid_t end_id = 0;
    /// Whether the walk reads the node or chain backward.
    bool backward = false;

    /// Return true if the step is a chain symbol.
    bool is_chain() const { return end_id != 0; }

    bool operator==(const SymbolicStep& o) const = default;
};

/// A traversal with child chains collapsed to symbols.
using SymbolicAllele = vector<SymbolicStep>;

/**
 * Make the symbolic allele of a walk through a site, given the site's `children`.
 *
 * Wherever the walk enters a child site, emit one symbol for that site's chain, and resume at the
 * handle that leaves the chain.
 *
 * If the walk enters a chain and never leaves it, emit the entering handle as a node step and go
 * on. Dropping the rest of the walk could make different alleles compare equal.
 *
 * @param out_visit_ranges If not null, filled with the half-open range of the walk's handles that
 * each emitted step covers. The ranges partition [0, walk.size()) in order, so the sequence of a
 * step is that of its handles. The range of a chain symbol is [entry, exit). The exit bound
 * belongs to the next step, since the chain shares it with what follows.
 */
SymbolicAllele symbolic_allele(const HandleGraph& graph, const Traversal& walk,
                               const multipass::SiteChildren& children,
                               vector<pair<int, int>>* out_visit_ranges = nullptr);

/// Return true if the two walks take the same route through the site at this level of the snarl
/// tree.
bool symbolically_equal(const HandleGraph& graph, const Traversal& a, const Traversal& b,
                        const multipass::SiteChildren& children);

/**
 * One difference between two symbolic alleles, as a half-open step range on each side.
 *
 * Either range may be empty. An empty ref range is an insertion, and an empty alt range a
 * deletion.
 */
struct DiffBlock {
    int ref_begin = 0;
    int ref_end = 0;
    int alt_begin = 0;
    int alt_end = 0;

    /// Return true if the block covers no reference steps.
    bool ref_empty() const { return ref_end == ref_begin; }
    /// Return true if the block covers no alt steps.
    bool alt_empty() const { return alt_end == alt_begin; }

    bool operator==(const DiffBlock& o) const = default;

    /// Print the block in a readable form.
    friend ostream& operator<<(ostream& out, const DiffBlock& block) {
        return out << "[ref " << block.ref_begin << "," << block.ref_end << " alt "
                   << block.alt_begin << "," << block.alt_end << "]";
    }
};

/**
 * Align two symbolic alleles and return the blocks where they differ, in reference order.
 *
 * The steps between blocks match on both sides, so no blocks means the two alleles take the same
 * route. The alignment minimises edit distance with substitutions at cost 1, so that a replacement
 * is one block. With insertions and deletions only, [a,b] against [b,b] would have a second
 * minimal alignment, giving two blocks around a spurious match. Ties are broken the same way every
 * time, preferring a match, then a substitution, then a deletion, then an insertion, since the
 * result decides how many records a site writes.
 *
 * The dynamic program runs inside Ukkonen's band, in O((|ref| + |alt|) x D) time and O(|ref| x D)
 * space, where D is the edit distance.
 *
 * @param out_alt_before_ref If not null, filled with |ref| + 1 entries. Entry i is the number of
 * alt steps consumed before reference step i, not counting any inserted at boundary i. It turns a
 * reference step range into the alt step range aligned with it.
 */
vector<DiffBlock> symbolic_diff(const SymbolicAllele& ref, const SymbolicAllele& alt,
                                vector<int>* out_alt_before_ref = nullptr);

/// Print a symbolic allele in a readable form.
ostream& operator<<(ostream& out, const SymbolicAllele& allele);

}

namespace std {
/// Hash a symbolic allele with FNV-1a over the fields of its steps.
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
