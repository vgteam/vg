#ifndef VG_SITE_VALUES_HPP_INCLUDED
#define VG_SITE_VALUES_HPP_INCLUDED

/** \file
 * The plain values the read-likelihood caller describes a site with, so that it needs no snarl
 * decomposition type once a site has been found.
 *
 * A site is its `SiteBounds`. An allele of a site is a `Traversal`: a walk of handles from the
 * site's start bound to its end bound. A star allele, for a strand that does not cross the site,
 * is an empty walk. The sites nested in a site sit in chains, and `SiteChildren` holds what the
 * site's records need to know about them.
 *
 * The conversions between these values and `Snarl` and `SnarlTraversal` are in namespace vg; the
 * caller's own types, and the functions that read them from a decomposition, are in
 * vg::multipass.
 */

#include <utility>
#include <vector>

#include "decomposition_sites.hpp"
#include "handle.hpp"
#include "snarls.hpp"
#include "traversal_finder.hpp"
#include <vg/vg.pb.h>

namespace vg {
using namespace std;

/// The bounds of `snarl` in `graph`.
SiteBounds bounds_of(const HandleGraph& graph, const Snarl& snarl);

/// `site` as a Snarl, for the code that still takes one.
Snarl snarl_of(const HandleGraph& graph, const SiteBounds& site);

/// A SnarlTraversal of node visits as a walk.
Traversal walk_of(const HandleGraph& graph, const SnarlTraversal& trav);

/// `walk_of` for each traversal.
vector<Traversal> walks_of(const HandleGraph& graph, const vector<SnarlTraversal>& travs);

/// A walk as a SnarlTraversal of node visits, for the code that still takes one.
SnarlTraversal snarl_traversal_of(const HandleGraph& graph, const Traversal& walk);

/// `snarl_traversal_of` for each walk.
vector<SnarlTraversal> snarl_traversals_of(const HandleGraph& graph,
                                           const vector<Traversal>& walks);

/// Whether a SnarlTraversal of node visits and a walk visit the same nodes in the same
/// orientations.
bool same_walk(const HandleGraph& graph, const SnarlTraversal& trav, const Traversal& walk);

/// The sequence of each step of an allele's walk, each read the way the walk reads its node, for
/// an allele that its walk does not spell: one that no graph walk spells, such as a walk with a
/// base changed inside one of its nodes. Empty for an allele that its walk spells.
using StepSequences = vector<string>;

/// The `StepSequences` of each of a site's alleles, by allele. Empty when every allele is spelled
/// by its walk, as every allele of a site the decomposition holds is. An allele with step
/// sequences is spelled, scored and written as them, and no panel haplotype carries it, whichever
/// walk it takes.
using AlleleSequences = vector<StepSequences>;

/// The step sequences of allele `allele` in `sequences`, or null when its walk spells it or
/// `sequences` is null.
const StepSequences* own_sequences(const AlleleSequences* sequences, size_t allele);

/// The sequence an allele spells: its step sequences `own` when not null, or else its walk's.
string allele_sequence(const HandleGraph& graph, const Traversal& walk, const StepSequences* own);

/// For a site inside a node, how many of the node's bases come before the site's first base,
/// read the way `site.start` reads the node. 0 for any other site.
size_t bases_before_site(const HandleGraph& graph, const SiteBounds& site);

/// A site inside nodes whose alternate allele changes their bases from `site.offset` on: its
/// bounds, its two alleles, REF then ALT, which both take the same walk, and the ALT's spelling,
/// as step sequences. A substitution lies inside one node. An indel's walk runs from its anchor
/// base's node to the node of the first base after its repeat tract, so the site holds every
/// placement of the indel and a base on each side; the bounds are the walk's first and last
/// handles.
struct EditSite {
    SiteBounds bounds;
    vector<Traversal> travs;
    AlleleSequences sequences;
};

/// The edit site that replaces the base at forward offset `offset` of the node `node` reads with
/// `alt`, given along the node's forward strand. The site is read the way `node` reads the node.
EditSite snv_edit_site(const HandleGraph& graph, handle_t node, uint32_t offset, char alt);

/// The edit site of an indel along `walk`, nodes read forward: `sequence` deleted, or inserted,
/// after the base at forward offset `offset` of the walk's first node, its anchor.
EditSite indel_edit_site(const HandleGraph& graph, const vector<nid_t>& walk, uint32_t offset,
                         bool deletion, const string& sequence);

/// Whether `site` lies inside one node.
bool inside_one_node(const SiteBounds& site);

/// The ID of an edit site: `op<node>.<offset>.<kind><alt>`, with the offset and the ALT along the
/// node's forward strand. It has no `_`, which `block_site_name` cuts at, and does not read as a
/// pair of boundary visits.
string edit_site_id(nid_t node, uint32_t offset, char kind, const string& alt);

namespace multipass {


/// A chain of sites nested in a site, by its bounds as `oriented_bounds` gives them. A chain of
/// one site has that site's bounds. A symbolic allele names a chain by the node IDs of these
/// bounds.
using ChildChain = SiteBounds;

/**
 * What a site's records need to know about the sites nested in it: for each way a walk can enter
 * a child site, the chain that child site belongs to. The symbolic layer reads it to replace a
 * walk's passage through a child chain with one symbol.
 */
struct SiteChildren {
    /// Whether the decomposition knows the site. When it does not, no child site is recognised,
    /// and every walk through the site projects to its plain node list.
    bool known = false;
    /// Whether the site is known only with its bounds swapped and reversed: the site was turned
    /// round to run forward along the reference path.
    bool reversed = false;
    /// The child chains.
    vector<ChildChain> chains;
    /// Each node and orientation by which a walk enters a child site, reading the child's start
    /// forward or its end backward, with the index in `chains` of that child site's chain. Sorted
    /// by node and orientation.
    vector<pair<pair<nid_t, bool>, size_t>> entries;

    /// The chain whose child site a walk enters by reading node `id` in orientation `backward`,
    /// or null when that enters no child site.
    const ChildChain* entered_by(nid_t id, bool backward) const;
};

/// A site as the read-likelihood caller visits it: the site in `decomposition`, its bounds as the
/// decomposition orients it, and the bounds of the sites enclosing it, innermost first.
struct SiteView {
    net_handle_t net;
    SiteBounds bounds;
    vector<SiteBounds> enclosing;
};

/// A site nested in another, as the decomposition gives it: the site, its bounds as the
/// decomposition orients it, and the bounds of the chain it is in.
struct ChildSite {
    net_handle_t net;
    SiteBounds bounds;
    ChildChain chain;
};

/// The bounds of the site or chain `net`, start to end as `decomposition` orients it, whichever
/// way `net` traverses it.
SiteBounds oriented_bounds(const SnarlDecomposition& decomposition, const HandleGraph& graph,
                           const net_handle_t& net);

/// The children of the site `site` in `decomposition`, read by a caller that has the site as
/// `as_called`: its oriented bounds, or those turned round. Every pair of consecutive nodes of a
/// child chain has a site between them, which the decomposition may not show, as SnarlManager's
/// adapter hides trivial snarls; both of that site's ways in are entries. If `child_sites` is
/// given, the child sites the decomposition shows are added to it, in its order.
SiteChildren site_children(const SnarlDecomposition& decomposition, const HandleGraph& graph,
                           const net_handle_t& site, const SiteBounds& as_called,
                           vector<ChildSite>* child_sites = nullptr);

/// The bounds of the sites enclosing `site` in `decomposition`, innermost first.
vector<SiteBounds> enclosing_sites(const SnarlDecomposition& decomposition,
                                   const HandleGraph& graph, const net_handle_t& site);

/// The chain the site `site` is in, by its oriented bounds.
ChildChain chain_of_site(const SnarlDecomposition& decomposition, const HandleGraph& graph,
                         const net_handle_t& site);

/// Where a node lies in a decomposition: the sites that hold it and the chain it is in.
struct NodePlacement {
    /// The innermost site that holds the node, or the decomposition's root when no site does.
    net_handle_t site;
    /// The bounds of `site` and of the sites enclosing it, innermost first, as the decomposition
    /// orients them, as `SiteView::enclosing` holds them for a site placed on the node. Empty when
    /// no site holds the node.
    vector<SiteBounds> enclosing;
    /// The chain the node is in, by its oriented bounds: the chain of the sites the node bounds,
    /// or, for a node that bounds no site, the node's own chain of one node.
    ChildChain chain;
};

/// Where node `id` lies in `decomposition`. A node that bounds a site lies in that site's chain,
/// between the site and its neighbour in the chain, and is held by the chain's parent; any other
/// node is a chain of its own inside the innermost site that holds it.
NodePlacement placement_of_node(const SnarlDecomposition& decomposition, const HandleGraph& graph,
                                nid_t id);

/// The key that tells `chain` from the other chains under the same parent, from the node IDs of
/// its bounds as `oriented_bounds` gives them. The linkage model groups a parent's child sites by
/// it. A node that bounds no site is a chain of one node, whose bounds are that node at
/// both ends, and it has a key like any other chain.
size_t chain_key_of(const HandleGraph& graph, const ChildChain& chain);

}
}

#endif
