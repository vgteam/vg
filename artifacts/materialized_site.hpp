#ifndef VG_MATERIALIZED_SITE_HPP_INCLUDED
#define VG_MATERIALIZED_SITE_HPP_INCLUDED

/** \file
 * A site's alleles, and a view of a site with its alleles for code that works on one site at a
 * time.
 */

#include <vector>

#include "handle.hpp"
#include "site_values.hpp"
#include "traversal_finder.hpp"

namespace vg {

using namespace std;

/**
 * The alleles of a site: what the snarl decomposition does not know about it.
 *
 * The decomposition gives a site's bounds and its place in the tree. The alleles come from
 * elsewhere: a snarl's from the walks the candidate finder picks, and an edit site's from the
 * changes the reads make.
 */
struct CandidateAlleles {
    /// The walk of each allele, from the start bound of the site to its end bound.
    vector<Traversal> walks;
    /// The index of the reference allele in `walks`, or -1 for none.
    int ref_allele = -1;
    /// The step sequences of the alleles that have them, which an edit site's ALTs do.
    AlleleSequences sequences;

    /// Get the step sequences of allele `allele`, or null when its walk spells it.
    const StepSequences* own_sequences(size_t allele) const;
};

/**
 * A site with its alleles, written out for code that works on one site at a time.
 *
 * Code that walks the snarl tree holds a site as its net handle in the decomposition, and its
 * `CandidateAlleles`. Code that scores, genotypes or writes one site gets this instead, so that it
 * does not need the decomposition. A snarl and an edit site have the same form. The site is read in
 * one direction, which the walks of its alleles follow.
 *
 * At each end, the bounds are a node handle, oriented along the site, and the number of bases of
 * that node before the boundary. A snarl's bounds lie on the inward sides of its boundary nodes.
 * An edit site's bounds lie inside its first and last nodes.
 *
 * It refers to the alleles of its holder, and is valid while they are.
 */
struct MaterializedSite {
    /// The node holding the start boundary, oriented to read into the site.
    handle_t start;
    /// The number of bases of `start` before the start boundary.
    size_t start_offset;
    /// The node holding the end boundary, oriented to read out of the site.
    handle_t end;
    /// The number of bases of `end` before the end boundary.
    size_t end_offset;
    /// The alleles of the site.
    const CandidateAlleles& alleles;

    /// Return true if the bounds lie inside nodes, as an edit site's do, rather than on the sides
    /// of nodes, as a snarl's do.
    bool bounds_inside_nodes(const HandleGraph& graph) const;
};

}

#endif
