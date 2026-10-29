#ifndef VG_SITE_TREE_HPP_INCLUDED
#define VG_SITE_TREE_HPP_INCLUDED

/** \file
 * The tree of nested sites that a VCF's nesting tags are computed from.
 *
 * The nesting tags of a record (LV, CH, PS, RC, RS and RD) say where its site sits in this tree,
 * relative to the sites of the other records.
 */

#include <functional>

#include "../src/handle.hpp"
#include "../src/hash_map.hpp"
#include "../src/snarls.hpp"

namespace vg {

using namespace std;

/// The two boundary visits of a site, as node IDs and orientations.
struct SiteEnds {
    nid_t start_id;
    bool start_backward;
    nid_t end_id;
    bool end_backward;

    bool operator==(const SiteEnds& other) const = default;
};

/// The sites of a graph, nested in one another.
///
/// A site is an opaque pointer, valid while the tree exists.
class SiteTree {
public:
    using site_t = const void*;

    virtual ~SiteTree() = default;

    /// Call `visit` on every site.
    ///
    /// With `in_preorder`, call it on this thread, on each site before the sites inside it.
    /// Otherwise call it on several threads, in any order.
    virtual void for_each_site(const function<void(site_t)>& visit, bool in_preorder) const = 0;

    /// Get the site that `site` lies directly inside, or null for a site at the top of the tree.
    virtual site_t parent_of(site_t site) const = 0;

    /// Get the boundary visits of `site`.
    virtual SiteEnds ends_of(site_t site) const = 0;

    /// Get the ID that the records of `site` carry, or null if their ID is made from its boundary
    /// visits.
    ///
    /// Records are matched to a site with an ID of its own by that ID.
    virtual const string* id_of(site_t site) const { return nullptr; }
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

namespace std {

/// Hash a site's boundary visits.
template<>
struct hash<vg::SiteEnds> {
    size_t operator()(const vg::SiteEnds& ends) const {
        return hash<tuple<vg::nid_t, bool, vg::nid_t, bool>>()(
            make_tuple(ends.start_id, ends.start_backward, ends.end_id, ends.end_backward));
    }
};

}

#endif
