#include "symbolic_allele.hpp"

#include <algorithm>
#include <cstdint>
#include <limits>

#include <unordered_map>

namespace vg {

/// Index the positions in the walk at which each node ID appears, in increasing order.
///
/// Finding where a walk leaves a chain is then a lookup rather than a scan of the walk.
static unordered_map<nid_t, vector<int>> index_positions(const HandleGraph& graph,
                                                          const Traversal& walk) {
    unordered_map<nid_t, vector<int>> at;
    for (int i = 0; i < (int)walk.size(); ++i) {
        at[graph.get_id(walk[i])].push_back(i);
    }
    return at;
}

/**
 * Find where a walk leaves a chain it enters at position `i`, by visiting node `node`.
 *
 * The walk leaves at its next visit, after `i`, to whichever of the chain's bounds other than
 * `node` comes first. A chain can be crossed in either direction, so both bounds are candidates.
 *
 * @return The position of the leaving visit, or -1 if the walk never leaves.
 */
static int chain_exit(const unordered_map<nid_t, vector<int>>& at,
                      const pair<nid_t, nid_t>& bounds, nid_t node, int i) {
    int exit = -1;
    for (nid_t boundary : {bounds.second, bounds.first}) {
        if (boundary == node) {
            continue;
        }
        auto found = at.find(boundary);
        if (found == at.end()) {
            continue;
        }
        // The first visit past i is the nearest exit. A node can recur thousands of times in one
        // walk through a satellite repeat, so find it by binary search.
        auto it = std::ranges::upper_bound(found->second, i);
        if (it != found->second.end() && (exit < 0 || *it < exit)) {
            exit = *it;
        }
    }
    return exit;
}

SymbolicAllele symbolic_allele(const HandleGraph& graph, const Traversal& walk,
                               const multipass::SiteChildren& children,
                               vector<pair<int, int>>* out_visit_ranges) {
    SymbolicAllele out;
    if (out_visit_ranges != nullptr) {
        out_visit_ranges->clear();
    }
    unordered_map<nid_t, vector<int>> at = index_positions(graph, walk);

    int i = 0;
    while (i < (int)walk.size()) {
        const nid_t node = graph.get_id(walk[i]);
        const bool backward = graph.get_is_reverse(walk[i]);
        // Find the chain of the child site this handle enters, if any. Only a child of this site
        // may become a symbol. Comparing a chain's bounds with the site's would not do, since a
        // site that is itself in a longer chain would collapse its own interior into one symbol,
        // making all its alleles equal.
        const multipass::ChildChain* chain = children.entered_by(node, backward);
        int exit = -1;
        pair<nid_t, nid_t> bounds;
        if (chain != nullptr) {
            bounds = {graph.get_id(chain->start), graph.get_id(chain->end)};
            exit = chain_exit(at, bounds, node, i);
        }

        if (exit > i) {
            SymbolicStep step;
            step.id = bounds.first;
            step.end_id = bounds.second;
            // The walk reads the chain backward when it enters at the chain's end.
            step.backward = (node == bounds.second);
            out.push_back(step);
            if (out_visit_ranges != nullptr) {
                out_visit_ranges->emplace_back(i, exit);
            }
            // Resume at the exit bound, which belongs to both the chain and whatever follows it.
            i = exit;
        } else {
            SymbolicStep step;
            step.id = node;
            step.backward = backward;
            out.push_back(step);
            if (out_visit_ranges != nullptr) {
                out_visit_ranges->emplace_back(i, i + 1);
            }
            ++i;
        }
    }
    return out;
}

bool symbolically_equal(const HandleGraph& graph, const Traversal& a, const Traversal& b,
                        const multipass::SiteChildren& children) {
    return symbolic_allele(graph, a, children) == symbolic_allele(graph, b, children);
}


namespace {

/**
 * A matrix of the edit distances between prefixes of two symbolic alleles, filled only inside
 * Ukkonen's band.
 *
 * Every cell on an optimal path has |i - j| <= D, where D is the edit distance, so a band of
 * half-width k >= D holds the optimum. Cells outside the band read as `INF`. Their true values
 * exceed D, so a traceback's tests fail for them as they would with the true values, and the
 * traceback matches that of the full matrix.
 */
class BandedEditMatrix {
public:
    /// A cost larger than any edit distance, which still has room to add to.
    static constexpr uint32_t INF = std::numeric_limits<uint32_t>::max() / 4;

    /// Fill the matrix for `ref` against `alt`, both nonempty.
    ///
    /// The band starts at the length difference, the least any alignment must spend, and doubles
    /// until the corner value is at most its half-width, which certifies the optimum.
    BandedEditMatrix(const SymbolicAllele& ref, const SymbolicAllele& alt) : ref(ref), alt(alt) {
        const size_t m = ref.size();
        const size_t n = alt.size();
        const size_t max_k = std::max(m, n);
        band_k = std::max<size_t>(1, m > n ? m - n : n - m);
        while (true) {
            fill();
            if (at(m, n) <= (uint32_t)band_k || band_k >= max_k) {
                break;
            }
            band_k = std::min(max_k, band_k * 2);
        }
    }

    /// Get the edit distance between the first `i` steps of ref and the first `j` of alt, or
    /// `INF` outside the band.
    uint32_t at(size_t i, size_t j) const {
        const long long offset = (long long)j - (long long)i + (long long)band_k;
        if (offset < 0 || offset >= (long long)stride) {
            return INF;
        }
        return cost[i * stride + (size_t)offset];
    }

private:
    /// Fill the cells inside the band of half-width `band_k`.
    void fill() {
        const size_t m = ref.size();
        const size_t n = alt.size();
        stride = 2 * band_k + 1;
        cost.assign(stride * (m + 1), INF);
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
                put(i, j, best_move(i, j));
            }
        }
    }

    /// Get the cost of the cheapest move into cell (i, j) from its filled neighbours.
    uint32_t best_move(size_t i, size_t j) const {
        uint32_t diag = at(i - 1, j - 1);
        uint32_t del = at(i - 1, j);
        uint32_t ins = at(i, j - 1);
        if (diag < INF) {
            diag += (ref[i - 1] == alt[j - 1] ? 0u : 1u);
        }
        if (del < INF) {
            del += 1u;
        }
        if (ins < INF) {
            ins += 1u;
        }
        return std::min(diag, std::min(del, ins));
    }

    /// Set cell (i, j), which must be inside the band.
    void put(size_t i, size_t j, uint32_t v) {
        cost[i * stride + (size_t)((long long)j - (long long)i + (long long)band_k)] = v;
    }

    const SymbolicAllele& ref;
    const SymbolicAllele& alt;
    /// The half-width of the band.
    size_t band_k = 0;
    /// The number of cells stored per row, 2 * `band_k` + 1.
    size_t stride = 0;
    /// The cells of each row inside the band, row by row.
    vector<uint32_t> cost;
};

/// One operation of an alignment of two symbolic alleles.
enum class EditOp { MATCH, SUB, DEL, INS };

/**
 * Trace an optimal alignment back through a filled matrix.
 *
 * Break ties diagonal first (a match before a substitution), then deletion, then insertion,
 * which is the order `symbolic_diff` promises.
 *
 * @return The operations of the alignment, from the start of both alleles.
 */
vector<EditOp> traceback(const BandedEditMatrix& matrix, const SymbolicAllele& ref,
                         const SymbolicAllele& alt) {
    vector<EditOp> ops;
    ops.reserve(ref.size() + alt.size());
    size_t i = ref.size();
    size_t j = alt.size();
    while (i > 0 || j > 0) {
        if (i > 0 && j > 0) {
            const bool equal = ref[i - 1] == alt[j - 1];
            if (matrix.at(i, j) == matrix.at(i - 1, j - 1) + (equal ? 0u : 1u)) {
                ops.push_back(equal ? EditOp::MATCH : EditOp::SUB);
                --i;
                --j;
                continue;
            }
        }
        if (i > 0 && matrix.at(i, j) == matrix.at(i - 1, j) + 1u) {
            ops.push_back(EditOp::DEL);
            --i;
            continue;
        }
        // Here j > 0, since the row and column initialisations make an insertion legal.
        ops.push_back(EditOp::INS);
        --j;
    }
    std::ranges::reverse(ops);
    return ops;
}

/**
 * Count, for each reference step, the alt steps an alignment consumes before it.
 *
 * Entry i is the alt index on arrival at reference step i. That is after everything that consumed
 * reference steps below i, and before anything inserted at boundary i, so that an insertion at a
 * boundary belongs to the block that owns the boundary.
 *
 * @return |ref| + 1 entries.
 */
vector<int> alt_before_ref(const vector<EditOp>& ops, size_t ref_size) {
    vector<int> out(ref_size + 1, 0);
    size_t ri = 0;
    size_t ai = 0;
    for (EditOp op : ops) {
        if (op == EditOp::INS) {
            ++ai;
        } else if (op == EditOp::DEL) {
            ++ri;
            out[ri] = (int)ai;
        } else {
            ++ri;
            ++ai;
            out[ri] = (int)ai;
        }
    }
    return out;
}

/**
 * Make a difference block of each maximal run of operations other than matches.
 *
 * A substitution joins the run it is in, so a mismatch inside a longer difference gives one
 * record rather than three.
 */
vector<DiffBlock> blocks_of(const vector<EditOp>& ops) {
    vector<DiffBlock> out;
    size_t ri = 0;
    size_t ai = 0;
    size_t k = 0;
    while (k < ops.size()) {
        if (ops[k] == EditOp::MATCH) {
            ++ri;
            ++ai;
            ++k;
            continue;
        }
        DiffBlock block;
        block.ref_begin = (int)ri;
        block.alt_begin = (int)ai;
        for (; k < ops.size() && ops[k] != EditOp::MATCH; ++k) {
            if (ops[k] != EditOp::INS) {
                ++ri;
            }
            if (ops[k] != EditOp::DEL) {
                ++ai;
            }
        }
        block.ref_end = (int)ri;
        block.alt_end = (int)ai;
        out.push_back(block);
    }
    return out;
}

} // namespace

vector<DiffBlock> symbolic_diff(const SymbolicAllele& ref, const SymbolicAllele& alt,
                                vector<int>* out_alt_before_ref) {
    const size_t m = ref.size();
    const size_t n = alt.size();
    if (m == 0 || n == 0) {
        // The alt is wholly an insertion or wholly a deletion, or the alleles are both empty.
        if (out_alt_before_ref != nullptr) {
            out_alt_before_ref->assign(m + 1, 0);
            (*out_alt_before_ref)[m] = (int)n;
        }
        if (m == 0 && n == 0) {
            return {};
        }
        return {DiffBlock{0, (int)m, 0, (int)n}};
    }

    const BandedEditMatrix matrix(ref, alt);
    const vector<EditOp> ops = traceback(matrix, ref, alt);
    if (out_alt_before_ref != nullptr) {
        *out_alt_before_ref = alt_before_ref(ops, m);
    }
    return blocks_of(ops);
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
