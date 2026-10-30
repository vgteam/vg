#ifndef VG_MATERIALIZED_SITE_HPP_INCLUDED
#define VG_MATERIALIZED_SITE_HPP_INCLUDED

/** \file
 * A site's alleles, and a view of a site with its alleles for code that works on one site at a
 * time.
 */

#include <string>
#include <vector>

#include "../src/handle.hpp"
#include "../src/traversal_finder.hpp"

namespace vg {

using namespace std;

/// The step sequences of an allele, in the order of its walk.
///
/// The step sequence of a visit of the walk holds the bases that the allele has at the visited
/// node, read the way the walk reads the node.
///
/// A walk spells an allele by reading the sequence of each node it visits. An edit found in the
/// reads, such as a substitution or a short indel, changes bases inside a node. An allele that
/// carries the edit visits the same nodes as an allele that does not, so its walk does not spell
/// it, and its step sequences do. A step sequence can be longer or shorter than the sequence of its
/// node, or empty, where the allele inserts or deletes bases.
///
/// The vector has one string for each visit of the walk, or is empty. An empty vector means that
/// the walk of the allele spells it.
using StepSequences = vector<string>;

/// The step sequences of each allele of a site, indexed by allele and then by visit of its walk.
///
/// An allele is numbered by the index of its walk among the walks of the site, not by its VCF
/// allele number. `sequences[a][i]` is the step sequence of visit `i` of the walk of allele `a`. An
/// allele whose entry is empty, or lies past the end of the vector, has no step sequences, and its
/// walk spells it. An empty vector means that the walks spell every allele.
///
/// An allele with step sequences is spelled and written as them, and reads are scored against
/// them. No panel haplotype is taken to carry it.
using AlleleSequences = vector<StepSequences>;

/// Get the step sequences of allele `allele` in `sequences`, or null if the allele has none or
/// `sequences` is null.
const StepSequences* own_sequences(const AlleleSequences* sequences, size_t allele);

/// Get the sequence an allele spells, from its step sequences `own` if not null, or else from its
/// walk.
string allele_sequence(const HandleGraph& graph, const Traversal& walk, const StepSequences* own);

/**
 * One change an allele makes to the sequence of one step of its walk.
 *
 * The sequence of a step is the sequence of its node, read the way the walk reads the node. The
 * change removes `removed` bases of that sequence from offset `at` on, and puts `added` in their
 * place. It inserts bases where `removed` is 0, and deletes them where `added` is empty.
 *
 * The offsets count bases of the node's own sequence, not of the sequence after other changes, so
 * a change means the same thing whatever other changes the allele makes.
 */
struct AlleleEdit {
    /// The index of the step in the walk.
    size_t step;
    /// The offset in the step's sequence where the change starts.
    size_t at;
    /// The number of bases the change removes.
    size_t removed;
    /// The bases the change puts in their place.
    string added;

    bool operator==(const AlleleEdit& other) const = default;
};

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
    /**
     * The changes each allele makes to the sequences of the steps of its walk, in any order.
     *
     * No two changes to one step may remove the same base. An insertion at the offset where
     * another change starts removing bases goes before the bases that change puts in.
     *
     * An edit site's ALTs differ from the graph inside nodes, so their walks do not spell them and
     * they have changes here. An allele whose entry is empty, or lies past the end of the vector,
     * is spelled by its walk. A snarl's alleles all are, and the vector is empty.
     */
    vector<vector<AlleleEdit>> edits;

    /// Return true if the walk of allele `allele` spells it, with no changes.
    bool spelled_by_walk(size_t allele) const;

    /// Get the sequence of each step of the walk of allele `allele`, with the allele's changes
    /// made, or an empty vector when its walk spells it.
    StepSequences step_sequences(const HandleGraph& graph, size_t allele) const;

    /// Get the step sequences of every allele, indexed like `walks`, or an empty vector when every
    /// walk spells its allele.
    AlleleSequences all_step_sequences(const HandleGraph& graph) const;

    /// Get the sequence allele `allele` spells.
    string sequence(const HandleGraph& graph, size_t allele) const;
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
