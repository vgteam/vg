#include "symbolic_allele.hpp"

#include <algorithm>
#include <cstdint>
#include <limits>

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

/// The chain a child snarl belongs to, as its boundary node ids. A snarl that is in no chain the
/// manager knows is a chain of its own, so its own boundaries are used.
static pair<nid_t, nid_t> chain_bounds(const Snarl* child, const SnarlManager& snarl_manager) {
    const Chain* chain = snarl_manager.chain_of(child);
    if (chain != nullptr && !chain->empty()) {
        return make_pair(get_start_of(*chain).node_id(), get_end_of(*chain).node_id());
    }
    return make_pair(child->start().node_id(), child->end().node_id());
}

/// The snarl `site` refers to, as the manager knows it, or null when the two disagree. Shared by
/// projection and `symbolic_site_resolvable`, so that they agree.
static const Snarl* resolve_site(const Snarl& site, const SnarlManager& snarl_manager,
                                 bool* out_reversed = nullptr) {
    if (out_reversed != nullptr) {
        *out_reversed = false;
    }
    const Snarl* site_ptr = snarl_manager.into_which_snarl(site.start().node_id(),
                                                          site.start().backward());
    if (site_ptr == nullptr) {
        return nullptr;
    }
    // Check that the node led to this snarl rather than to another one reachable from it, by
    // comparing boundary nodes, in either order, since a snarl and its reversal are the same snarl.
    // `flip_snarl` (graph_caller.cpp) reverses a snarl whose reference path runs backwards, and the
    // caller then works on the reversed copy, whose start node is the original end node. The
    // boundary index maps a snarl's start and its reversed end to the same snarl, and everything
    // downstream is independent of the site's orientation.
    auto same_visit = [](const Visit& a, const Visit& b) {
        return a.node_id() == b.node_id() && a.backward() == b.backward();
    };
    // The forward pairing compares node ids only.
    const bool forward = site_ptr->start().node_id() == site.start().node_id() &&
                         site_ptr->end().node_id() == site.end().node_id();
    // The reversed pairing compares full visits, orientation included: the reversed snarl is the
    // canonical one with its boundaries reversed and swapped.
    const bool reversed = same_visit(site_ptr->start(), reverse(site.end())) &&
                          same_visit(site_ptr->end(), reverse(site.start()));
    if (!forward && !reversed) {
        return nullptr;
    }
    if (out_reversed != nullptr) {
        *out_reversed = reversed && !forward;
    }
    return site_ptr;
}

bool symbolic_site_resolvable(const Snarl& site, const SnarlManager& snarl_manager,
                              bool* out_reversed) {
    return resolve_site(site, snarl_manager, out_reversed) != nullptr;
}

pair<nid_t, nid_t> chain_bounds_of(const Snarl* child, const SnarlManager& snarl_manager) {
    return chain_bounds(child, snarl_manager);
}

SymbolicAllele symbolic_allele(const SnarlTraversal& trav, const Snarl& site,
                               const SnarlManager& snarl_manager,
                               vector<pair<int, int>>* out_visit_ranges) {
    SymbolicAllele out;
    if (out_visit_ranges != nullptr) {
        out_visit_ranges->clear();
    }
    unordered_map<nid_t, vector<int>> at = index_positions(trav);

    // The snarl we are projecting, as the manager knows it, so that a snarl the traversal enters can
    // be tested for being a child of this site. Null means no child is recognised and the
    // projection is the plain node list (see `symbolic_site_resolvable`).
    const Snarl* site_ptr = resolve_site(site, snarl_manager);

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
            if (out_visit_ranges != nullptr) {
                out_visit_ranges->emplace_back(i, i + 1);
            }
            ++i;
            continue;
        }

        const nid_t node = v.node_id();
        // The snarl entered by traversing into this node in this orientation, if any.
        const Snarl* child = snarl_manager.into_which_snarl(node, v.backward());
        bool symbolised = false;

        // Only a child of this site may become a symbol. Comparing the chain's boundaries with the
        // site's is not enough: a site that is itself in a longer chain would see that chain's
        // boundaries and collapse its own interior into one symbol, making all its alleles equal.
        bool is_child = child != nullptr && site_ptr != nullptr &&
                        snarl_manager.parent_of(child) == site_ptr;
        if (is_child) {
            pair<nid_t, nid_t> bounds = chain_bounds(child, snarl_manager);
            {
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
                    // index_positions fills each vector in increasing visit order, so the first
                    // entry past i is the nearest exit, found by binary search; a node can recur
                    // thousands of times in one traversal through a satellite repeat.
                    auto it = std::upper_bound(found->second.begin(), found->second.end(), i);
                    if (it != found->second.end() && (exit < 0 || *it < exit)) {
                        exit = *it;
                    }
                }
                if (exit > i) {
                    SymbolicStep step;
                    step.id = bounds.first;
                    step.end_id = bounds.second;
                    // The chain is traversed backward when its recorded end is met before its start.
                    step.backward = (node == bounds.second);
                    out.push_back(step);
                    if (out_visit_ranges != nullptr) {
                        out_visit_ranges->emplace_back(i, exit);
                    }
                    // Resume *at* the exit boundary, which belongs to both the chain and whatever
                    // follows it, so it is not consumed.
                    i = exit;
                    symbolised = true;
                }
                // A chain entered and not left within this traversal falls through and is emitted as
                // a plain node, so that the rest of the traversal is kept.
            }
        }

        if (!symbolised) {
            SymbolicStep step;
            step.id = node;
            step.backward = v.backward();
            out.push_back(step);
            if (out_visit_ranges != nullptr) {
                out_visit_ranges->emplace_back(i, i + 1);
            }
            ++i;
        }
    }
    return out;
}

bool symbolically_equal(const SnarlTraversal& a, const SnarlTraversal& b, const Snarl& site,
                        const SnarlManager& snarl_manager) {
    return symbolic_allele(a, site, snarl_manager) == symbolic_allele(b, site, snarl_manager);
}


vector<DiffBlock> symbolic_diff(const SymbolicAllele& ref, const SymbolicAllele& alt,
                                vector<int>* out_alt_before_ref) {

    const size_t m = ref.size();
    const size_t n = alt.size();

    // Filled for every exit path, so a caller never reads a stale or short vector.
    auto trivial_map = [&](size_t consumed_before_end) {
        if (out_alt_before_ref == nullptr) {
            return;
        }
        out_alt_before_ref->assign(m + 1, 0);
        (*out_alt_before_ref)[m] = (int)consumed_before_end;
    };

    if (m == 0 && n == 0) {
        trivial_map(0);
        return {};
    }
    if (m == 0 || n == 0) {
        // Wholly an insertion or wholly a deletion. No alignment to compute, and not a degradation.
        trivial_map(n);
        return {DiffBlock{0, (int)m, 0, (int)n}};
    }

    // Ukkonen's band. Every cell on an optimal path has |i - j| <= D, where D is the edit
    // distance, so a band of half-width k >= D holds the optimum. k starts at the length difference,
    // the least any alignment must spend, and doubles until the corner value is at most k, which
    // certifies it. Work is O((m + n) * D) rather than O(m * n).
    //
    // Cells outside the band read as INF. Their true values exceed D, so the traceback's tests
    // fail for them just as they would with the true values, and the traceback matches the full
    // matrix's.
    const uint32_t INF = std::numeric_limits<uint32_t>::max() / 4;
    const size_t max_k = std::max(m, n);
    size_t band_k = std::max<size_t>(1, m > n ? m - n : n - m);
    vector<uint32_t> cost;
    size_t stride = 0;

    // Read: INF outside the band. Write: only ever called in band.
    auto cell = [&](size_t i, size_t j) -> uint32_t {
        long long off = (long long)j - (long long)i + (long long)band_k;
        if (off < 0 || off >= (long long)stride) {
            return INF;
        }
        return cost[i * stride + (size_t)off];
    };
    auto put = [&](size_t i, size_t j, uint32_t v) {
        cost[i * stride + (size_t)((long long)j - (long long)i + (long long)band_k)] = v;
    };

    while (true) {
        stride = 2 * band_k + 1;
        cost.assign(stride * (m + 1), INF);
        // Row and column initialisations, clipped to the band.
        for (size_t i = 0; i <= m && i <= band_k; ++i) {
            put(i, 0, (uint32_t)i);
        }
        for (size_t j = 0; j <= n && j <= band_k; ++j) {
            put(0, j, (uint32_t)j);
        }
        for (size_t i = 1; i <= m; ++i) {
            const size_t jlo = (i > band_k) ? i - band_k : 1;
            const size_t jhi = std::min(n, i + band_k);
            for (size_t j = jlo == 0 ? 1 : jlo; j <= jhi; ++j) {
                uint32_t diag = cell(i - 1, j - 1);
                uint32_t del = cell(i - 1, j);
                uint32_t ins = cell(i, j - 1);
                if (diag < INF) {
                    diag += (ref[i - 1] == alt[j - 1] ? 0u : 1u);
                }
                if (del < INF) {
                    del += 1u;
                }
                if (ins < INF) {
                    ins += 1u;
                }
                put(i, j, std::min(diag, std::min(del, ins)));
            }
        }
        if (cell(m, n) <= (uint32_t)band_k || band_k >= max_k) {
            break;
        }
        band_k = std::min(max_k, band_k * 2);
    }
    auto at = [&](size_t i, size_t j) -> uint32_t { return cell(i, j); };

    // Traceback, with the header's tie-break order: diagonal first (match before substitution),
    // then deletion, then insertion. It walks backwards, so the ops are reversed afterwards.
    enum Op { OP_MATCH, OP_SUB, OP_DEL, OP_INS };
    vector<Op> ops;
    ops.reserve(m + n);
    {
        size_t i = m;
        size_t j = n;
        while (i > 0 || j > 0) {
            if (i > 0 && j > 0) {
                bool equal = ref[i - 1] == alt[j - 1];
                if (at(i, j) == at(i - 1, j - 1) + (equal ? 0u : 1u)) {
                    ops.push_back(equal ? OP_MATCH : OP_SUB);
                    --i;
                    --j;
                    continue;
                }
            }
            if (i > 0 && at(i, j) == at(i - 1, j) + 1u) {
                ops.push_back(OP_DEL);
                --i;
                continue;
            }
            // j > 0 necessarily: the row and column initialisations make the remaining move legal.
            ops.push_back(OP_INS);
            --j;
        }
    }
    std::reverse(ops.begin(), ops.end());

    if (out_alt_before_ref != nullptr) {
        // Entry i is the alt index on arrival at reference step i: after everything that consumed
        // reference steps below i, and before anything inserted at boundary i, so that an insertion
        // at a boundary belongs to the block that owns the boundary.
        out_alt_before_ref->assign(m + 1, 0);
        size_t ri = 0;
        size_t ai = 0;
        for (Op op : ops) {
            if (op == OP_INS) {
                ++ai;
            } else if (op == OP_DEL) {
                ++ri;
                (*out_alt_before_ref)[ri] = (int)ai;
            } else {
                ++ri;
                ++ai;
                (*out_alt_before_ref)[ri] = (int)ai;
            }
        }
    }

    // Each maximal run of non-match ops becomes one block. A substitution joins the run it is in,
    // so a mismatch inside a longer difference gives one record rather than three.
    vector<DiffBlock> out;
    size_t ri = 0;
    size_t ai = 0;
    size_t k = 0;
    while (k < ops.size()) {
        if (ops[k] == OP_MATCH) {
            ++ri;
            ++ai;
            ++k;
            continue;
        }
        DiffBlock block;
        block.ref_begin = (int)ri;
        block.alt_begin = (int)ai;
        while (k < ops.size() && ops[k] != OP_MATCH) {
            if (ops[k] == OP_SUB) {
                ++ri;
                ++ai;
            } else if (ops[k] == OP_DEL) {
                ++ri;
            } else {
                ++ai;
            }
            ++k;
        }
        block.ref_end = (int)ri;
        block.alt_end = (int)ai;
        out.push_back(block);
    }
    return out;
}

ostream& operator<<(ostream& out, const DiffBlock& block) {
    return out << "[ref " << block.ref_begin << "," << block.ref_end
               << " alt " << block.alt_begin << "," << block.alt_end << "]";
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
