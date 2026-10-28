#ifndef VG_SITE_VALUES_HPP_INCLUDED
#define VG_SITE_VALUES_HPP_INCLUDED

/** \file
 * The types the read-likelihood caller copies a site's place in the decomposition into.
 *
 * The caller holds each site as its net handle in an `EditSnarlDecompositionOverlay`. It still
 * copies what a site's records need to know about the sites nested in it out of the overlay
 * (`SiteChildren`). An allele of a site is a `Traversal`, a walk of handles from the start bound of
 * the site to its end bound. A star allele, for a strand that does not cross the site, is an empty
 * walk.
 *
 * TODO: Hold a chain as its net handle too, so that `ChildChain` and `SiteBounds` are no longer
 * needed.
 *
 * Namespace vg holds how alleles are spelled and how edit sites are made. Namespace vg::multipass
 * holds the caller's own types, and the functions that read them from a decomposition.
 */

#include <utility>
#include <vector>

#include "decomposition_sites.hpp"
#include "edit_snarl_decomposition_overlay.hpp"
#include "handle.hpp"
#include "snarls.hpp"
#include "traversal_finder.hpp"
#include <vg/vg.pb.h>

namespace vg {
using namespace std;

/**
 * The sequence an allele has at each node its walk visits, for an allele the graph cannot spell.
 *
 * An allele of an edit site differs from the graph inside a node, such as by one changed base, so
 * no walk of the graph spells it. The allele still has a walk, which says which nodes it covers.
 * Element `i` is the allele's sequence where the walk visits node `walk[i]`, read the way the walk
 * reads that node: the node's sequence with the edit applied.
 *
 * An allele whose walk spells it has none, and this is empty.
 */
using StepSequences = vector<string>;

/**
 * The `StepSequences` of each allele of a site, indexed like the site's walks.
 *
 * For a site of the decomposition, every walk spells its allele, and this is empty. For an edit
 * site, the reference allele's entry is empty, and each other allele's entry holds its edited
 * sequences.
 *
 * Code that spells, scores or writes an allele uses its step sequences where it has them (see
 * `own_sequences`). No panel haplotype carries an allele that has them.
 */
using AlleleSequences = vector<StepSequences>;

/// Get the step sequences of allele `allele` in `sequences`, or null when its walk spells it or
/// `sequences` is null.
const StepSequences* own_sequences(const AlleleSequences* sequences, size_t allele);

/// Get the sequence an allele spells, from its step sequences `own` if not null, or else from its
/// walk.
string allele_sequence(const HandleGraph& graph, const Traversal& walk, const StepSequences* own);

/// Count the bases of the start node of a site inside nodes that come before the site's first base.
///
/// Count them the way `site.start` reads the node. Return 0 for any other site.
size_t bases_before_site(const HandleGraph& graph, const SiteBounds& site);

/// A site inside nodes, whose alternate alleles change the bases of the nodes.
///
/// All its alleles take the same walk. The bounds are the first and last handles of the walk. A
/// substitution lies inside one node. The walk of an indel runs from the node of its anchor base
/// to the node of the first base after its repeat tract. That way the site holds every placement
/// of the indel, and a base on each side.
struct EditSite {
    SiteBounds bounds;
    /// The walk of each allele, REF first.
    vector<Traversal> travs;
    /// The spelling of each allele. REF has none, because its walk spells it.
    AlleleSequences sequences;
};

/**
 * Make the edit site along `steps` whose ALTs `alts` each replace the `length` bases from `at` on.
 *
 * `at` counts bases along the steps as they read their nodes. Each ALT goes in the step holding
 * the base at `at`, or, for an insertion, the base before it. The replaced bases are cut from the
 * steps holding them.
 *
 * @param offset The first base of the site, along the forward strand of its first node.
 * @pre For an insertion (`length` 0), `at` is above 0.
 */
EditSite edit_site(const HandleGraph& graph, const Traversal& steps, uint32_t offset, size_t at,
                   size_t length, const vector<string>& alts);

/// Return true if `site` lies inside one node.
bool inside_one_node(const SiteBounds& site);

/**
 * Make the ID of an edit site, `op<node>.<offset>.<kind><alt>`.
 *
 * The offset and the ALT are along the forward strand of the node. A merged change has kind `M`,
 * and the number of bases it replaces in place of its ALT. The ID has no `_`, which
 * `block_site_name` cuts at, so it cannot be read as a pair of boundary visits.
 */
string edit_site_id(nid_t node, uint32_t offset, char kind, const string& alt);

namespace multipass {

/// A chain of sites nested in a site, as its bounds from `oriented_bounds`.
///
/// A chain of one site has that site's bounds. A symbolic allele names a chain by the node IDs of
/// these bounds.
using ChildChain = SiteBounds;

/**
 * A table of the ways into the child sites of a site, and the chains those child sites are in.
 *
 * It lets the records of a site replace the stretch of a walk inside a child chain with one
 * symbolic allele.
 */
struct SiteChildren {
    /// Whether the decomposition knows the site.
    ///
    /// For a site it does not know, the table is empty, so every walk through the site is written
    /// as its plain node list.
    bool known = false;
    /// Whether the decomposition knows the site only turned round, with its bounds swapped and
    /// flipped.
    ///
    /// The caller turns a site round to run forward along the reference path.
    bool reversed = false;
    /// The child chains.
    vector<ChildChain> chains;
    /// The entries, sorted by node and orientation.
    ///
    /// An entry is a node ID and orientation by which a walk enters a child site, with the index
    /// in `chains` of the chain that the child site is in.
    vector<pair<pair<nid_t, bool>, size_t>> entries;

    /// Get the chain whose child site a walk enters by reading node `id` in orientation
    /// `backward`, or null if that enters no child site.
    const ChildChain* entered_by(nid_t id, bool backward) const;
};

/// Get the bounds of the site or chain `net`, start to end as `decomposition` orients it,
/// whichever way `net` traverses it.
SiteBounds oriented_bounds(const SnarlDecomposition& decomposition, const HandleGraph& graph,
                           const net_handle_t& net);

/**
 * Get the `SiteChildren` of the site `site` in `decomposition`.
 *
 * Every two consecutive nodes of a child chain have a site between them, which the decomposition
 * may show as a trivial snarl or leave out. Both ways into that site are entries.
 *
 * @param as_called The bounds the caller has the site as: its oriented bounds, or those turned
 * round.
 * @param child_sites If not null, gets the child sites other than trivial snarls appended, in the
 * order of the decomposition.
 */
SiteChildren site_children(const SnarlDecomposition& decomposition, const HandleGraph& graph,
                           const net_handle_t& site, const SiteBounds& as_called,
                           vector<net_handle_t>* child_sites = nullptr);

/// Get the chain the site `site` is in.
ChildChain chain_of_site(const SnarlDecomposition& decomposition, const HandleGraph& graph,
                         const net_handle_t& site);

/// Get the chain that holds the site or edit `net` of `sites`.
///
/// An edit lies in the chain of the node that holds its start.
ChildChain chain_holding(const EditSnarlDecompositionOverlay& sites, const HandleGraph& graph,
                         const net_handle_t& net);

/// The place of a node in a decomposition: the sites that hold it and the chain it is in.
struct NodePlacement {
    /// The innermost site that holds the node, or the decomposition's root when no site does.
    net_handle_t site;
    /// The chain the node is in.
    ///
    /// For a node that bounds sites, this is the chain of those sites. For any other node, it is
    /// the trivial chain of the node.
    ChildChain chain;
};

/// Find the place of node `id` in `decomposition`.
///
/// A node that bounds a site lies in the chain of that site, and the parent of the chain holds it.
/// Any other node is a trivial chain inside the innermost site that holds it.
NodePlacement placement_of_node(const SnarlDecomposition& decomposition, const HandleGraph& graph,
                                nid_t id);

/// Make a key that tells `chain` from every other chain of the decomposition.
///
/// The key is the ID of the chain's first node. A trivial chain, of a node that bounds no site,
/// has a key like any other chain. A key is never 0, since node IDs start at 1.
size_t chain_key_of(const HandleGraph& graph, const ChildChain& chain);

}
}

#endif
