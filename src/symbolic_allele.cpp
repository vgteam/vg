#include "symbolic_allele.hpp"

#include <unordered_map>

namespace vg {

/// Where each node id appears in the traversal, so "does this chain close later" is a lookup
/// rather than a rescan. Traversals of a large snarl run to thousands of visits and every visit
/// asks the question.
static unordered_map<nid_t, vector<int>> index_positions(const SnarlTraversal& trav) {
    unordered_map<nid_t, vector<int>> at;
    for (int i = 0; i < trav.visit_size(); ++i) {
        const Visit& v = trav.visit(i);
        if (!v.has_snarl()) {
            at[v.node_id()].push_back(i);
        }
    }
    return at;
}

/// The chain a child snarl belongs to, as its boundary node ids. Falls back to the snarl's own
/// boundaries when it is not in a chain the manager knows about, which is the trivial-chain case:
/// chains_of() presents unary snarls and trivial-chain snarls as their own chains anyway, so the
/// fallback agrees with the general case rather than being a special one.
static pair<nid_t, nid_t> chain_bounds(const Snarl* child, const SnarlManager& snarl_manager) {
    const Chain* chain = snarl_manager.chain_of(child);
    if (chain != nullptr && !chain->empty()) {
        return make_pair(get_start_of(*chain).node_id(), get_end_of(*chain).node_id());
    }
    return make_pair(child->start().node_id(), child->end().node_id());
}

SymbolicAllele symbolic_allele(const SnarlTraversal& trav, const Snarl& site,
                               const SnarlManager& snarl_manager) {
    SymbolicAllele out;
    unordered_map<nid_t, vector<int>> at = index_positions(trav);

    const nid_t site_start = site.start().node_id();
    const nid_t site_end = site.end().node_id();

    int i = 0;
    while (i < trav.visit_size()) {
        const Visit& v = trav.visit(i);

        // A visit already carrying a Snarl is a symbol as it stands.
        if (v.has_snarl()) {
            SymbolicStep step;
            step.id = v.snarl().start().node_id();
            step.end_id = v.snarl().end().node_id();
            step.backward = v.backward();
            out.push_back(step);
            ++i;
            continue;
        }

        const nid_t node = v.node_id();
        // The snarl entered by traversing into this node in this orientation, if any.
        const Snarl* child = snarl_manager.into_which_snarl(node, v.backward());
        bool symbolised = false;

        if (child != nullptr) {
            pair<nid_t, nid_t> bounds = chain_bounds(child, snarl_manager);
            // Never symbolise the site itself: a traversal of a snarl enters at its own start, and
            // collapsing that would reduce every allele to a single symbol and make them all equal.
            bool is_own = (bounds.first == site_start && bounds.second == site_end) ||
                          (bounds.first == site_end && bounds.second == site_start);
            if (!is_own) {
                // Leave the chain at whichever of its boundaries this traversal reaches next; a
                // chain can be crossed in either direction, so both are candidates.
                int exit = -1;
                for (nid_t boundary : {bounds.second, bounds.first}) {
                    if (boundary == node) {
                        continue;
                    }
                    auto found = at.find(boundary);
                    if (found == at.end()) {
                        continue;
                    }
                    for (int j : found->second) {
                        if (j > i && (exit < 0 || j < exit)) {
                            exit = j;
                        }
                    }
                }
                if (exit > i) {
                    SymbolicStep step;
                    step.id = bounds.first;
                    step.end_id = bounds.second;
                    // The chain is traversed backward when its recorded end is met before its start.
                    step.backward = (node == bounds.second);
                    out.push_back(step);
                    // Resume *at* the exit boundary, which belongs to both the chain and whatever
                    // follows it, so it is not consumed.
                    i = exit;
                    symbolised = true;
                }
                // A chain entered and never left within this traversal falls through deliberately
                // and is emitted as a plain node. Swallowing the remainder would drop real
                // differences and make unrelated alleles compare equal.
            }
        }

        if (!symbolised) {
            SymbolicStep step;
            step.id = node;
            step.backward = v.backward();
            out.push_back(step);
            ++i;
        }
    }
    return out;
}

bool symbolically_equal(const SnarlTraversal& a, const SnarlTraversal& b, const Snarl& site,
                        const SnarlManager& snarl_manager) {
    return symbolic_allele(a, site, snarl_manager) == symbolic_allele(b, site, snarl_manager);
}

bool has_child_chain(const SnarlTraversal& trav, const Snarl& site,
                     const SnarlManager& snarl_manager) {
    SymbolicAllele allele = symbolic_allele(trav, site, snarl_manager);
    for (const SymbolicStep& step : allele) {
        if (step.is_chain()) {
            return true;
        }
    }
    return false;
}

ostream& operator<<(ostream& out, const SymbolicAllele& allele) {
    for (const SymbolicStep& s : allele) {
        out << (s.backward ? '<' : '>');
        if (s.is_chain()) {
            out << "C" << s.id << "_" << s.end_id;
        } else {
            out << s.id;
        }
    }
    return out;
}

}
