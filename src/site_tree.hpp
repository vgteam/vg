#ifndef VG_SITE_TREE_HPP_INCLUDED
#define VG_SITE_TREE_HPP_INCLUDED

/** \file
 * The tree of nested sites that a VCF's nesting tags are computed from.
 *
 * A site is a region of the graph between two boundary nodes, and sites nest: a site can lie inside
 * another. The nesting tags of a record (LV, CH, PS, RC, RS and RD) describe where its site sits in
 * this tree relative to the sites of the other records.
 */

#include <functional>

#include "handle.hpp"
#include "snarls.hpp"

namespace vg {

using namespace std;

/// The two boundary visits of a site: the node it starts at and the node it ends at, each with
/// whether the site reads it backward.
struct SiteEnds {
    nid_t start_id;
    bool start_backward;
    nid_t end_id;
    bool end_backward;
};

/// The sites of a graph, nested in one another. A site is an opaque pointer, valid while the tree
/// exists.
class SiteTree {
public:
    using site_t = const void*;

    virtual ~SiteTree() = default;

    /// Call `visit` on every site. With `in_preorder`, the calls are on this thread and reach each
    /// site before the sites inside it; otherwise they are on several threads, in any order.
    virtual void for_each_site(const function<void(site_t)>& visit, bool in_preorder) const = 0;

    /// The site that `site` lies directly inside, or null for a site at the top of the tree.
    virtual site_t parent_of(site_t site) const = 0;

    /// The boundary visits of `site`.
    virtual SiteEnds ends_of(site_t site) const = 0;
};

/// The SiteTree of the snarls a SnarlManager holds. The manager must outlive it.
class SnarlManagerSiteTree : public SiteTree {
public:
    explicit SnarlManagerSiteTree(const SnarlManager& manager);

    void for_each_site(const function<void(site_t)>& visit, bool in_preorder) const override;
    site_t parent_of(site_t site) const override;
    SiteEnds ends_of(site_t site) const override;

private:
    const SnarlManager& manager;
};

}

#endif
