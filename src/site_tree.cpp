#include "site_tree.hpp"

namespace vg {

SnarlManagerSiteTree::SnarlManagerSiteTree(const SnarlManager& manager) : manager(manager) {
    // Nothing to do
}

void SnarlManagerSiteTree::for_each_site(const function<void(site_t)>& visit,
                                         bool in_preorder) const {
    auto visit_snarl = [&](const Snarl* snarl) {
        visit(snarl);
    };
    if (in_preorder) {
        manager.for_each_snarl_preorder(visit_snarl);
    } else {
        manager.for_each_snarl_unindexed_parallel(visit_snarl);
    }
}

SiteTree::site_t SnarlManagerSiteTree::parent_of(site_t site) const {
    return manager.parent_of(static_cast<const Snarl*>(site));
}

SiteEnds SnarlManagerSiteTree::ends_of(site_t site) const {
    const Snarl* snarl = static_cast<const Snarl*>(site);
    return SiteEnds{snarl->start().node_id(), snarl->start().backward(), snarl->end().node_id(),
                    snarl->end().backward()};
}

}
