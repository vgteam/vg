#include <algorithm>
#include <limits>

#include "string_aligner.hpp"

extern "C" {
#include "wavefront/wavefront_align.h"
}

namespace vg {

namespace {

/// A cost no alignment has.
constexpr int32_t g_no_cost = 1 << 29;

/// Where the cost of a cell's match state came from, in the low bits of its trace byte.
enum TraceSource : uint8_t {
    FromDiagonal = 0,
    FromInsertion = 1,
    FromLongInsertion = 2,
    FromDeletion = 3,
    FromLongDeletion = 4,
    SourceMask = 7
};

/// The bits of a trace byte that say a gap state was reached by extending the same gap.
constexpr uint8_t g_insertion_extended = 8;
constexpr uint8_t g_long_insertion_extended = 16;
constexpr uint8_t g_deletion_extended = 32;
constexpr uint8_t g_long_deletion_extended = 64;

/// Add `length` bases of `op` to the end of `runs`, joining them to a last run of the same kind.
void add_run(vector<AlignmentRun>& runs, AlignmentOp op, size_t length) {
    if (length == 0) {
        return;
    }
    if (!runs.empty() && runs.back().op == op) {
        runs.back().length += length;
    } else {
        runs.push_back(AlignmentRun{op, length});
    }
}

/// The cost of a gap state reached by opening a gap from `from_match` or extending one from
/// `from_gap`, and whether extending is the cheaper.
struct GapStep {
    int32_t cost;
    bool extended;
};

GapStep gap_step(int32_t from_match, int32_t from_gap, int32_t open, int32_t extend) {
    const int32_t opened = from_match + open + extend;
    const int32_t extended = from_gap + extend;
    if (extended < opened) {
        return {std::min(extended, g_no_cost), true};
    }
    return {std::min(opened, g_no_cost), false};
}

/**
 * The costs of the cells of one row of an alignment matrix, by column.
 *
 * Entry `j + 1` is the cost of the cell in column `j`, and entry 0 that of column -1, which no
 * alignment reaches. Every cell outside the row's range holds `g_no_cost`, so a cell of the next
 * row can read the cells above it and above to its left without testing where they lie.
 */
struct CellRow {
    vector<int32_t> match;
    vector<int32_t> deletion;
    vector<int32_t> long_deletion;
    /// The columns the row last held, from `lo` to `hi`, or none where `lo > hi`.
    size_t lo = 1;
    size_t hi = 0;

    /// Make a row for a matrix of `m` + 1 columns, holding no cell.
    explicit CellRow(size_t m)
        : match(m + 2, g_no_cost), deletion(m + 2, g_no_cost), long_deletion(m + 2, g_no_cost) {}
};

/**
 * Fill the cells of row `i` from column `lo` to `hi` of the alignment of `a` to `b` into `here`,
 * from the row above, `above`, and write their trace bytes to `trace`, at `trace[j - trace_lo]`
 * for column `j`, if it is not null.
 *
 * A cell's match state comes from the cell above and to the left, then an insertion, a long
 * insertion, a deletion and a long deletion, each taken only where it is cheaper. The cells that
 * depend only on the row above are filled first, in a loop the compiler can vectorise, and the
 * insertions, which come from the cell to the left, then from left to right.
 */
void fill_row(const char* a, const char* b, size_t i, size_t lo, size_t hi,
              const StringAlignerParams& params, const CellRow& above, CellRow& here,
              uint8_t* trace, int64_t trace_lo) {
    // Clear what `here` held two rows ago.
    for (size_t j = here.lo; j <= here.hi; ++j) {
        here.match[j + 1] = g_no_cost;
        here.deletion[j + 1] = g_no_cost;
        here.long_deletion[j + 1] = g_no_cost;
    }
    here.lo = lo;
    here.hi = hi;
    const int32_t open = params.gap_open + params.gap_extend;
    const int32_t extend = params.gap_extend;
    const int32_t long_open = params.long_gap_open + params.long_gap_extend;
    const int32_t long_extend = params.long_gap_extend;
    const int32_t mismatch = params.mismatch;
    const char base = i > 0 ? a[i - 1] : 0;
    int32_t* match = here.match.data() + 1;
    int32_t* deletion = here.deletion.data() + 1;
    int32_t* long_deletion = here.long_deletion.data() + 1;
    const int32_t* match_up = above.match.data() + 1;
    const int32_t* deletion_up = above.deletion.data() + 1;
    const int32_t* long_deletion_up = above.long_deletion.data() + 1;
    // From the row above: the deletions and the diagonal, which `match` holds for now. Column 0
    // has no base of `b` before it.
    const size_t from = std::max<size_t>(lo, 1);
    if (lo == 0) {
        match[0] = i == 0 ? 0 : g_no_cost;
    }
    for (size_t j = from; j <= hi; ++j) {
        match[j] = std::min(g_no_cost, match_up[j - 1] + (base == b[j - 1] ? 0 : mismatch));
    }
    // The trace bytes start with whether each deletion extends one.
    uint8_t* traced = trace != nullptr ? trace + ((int64_t)lo - trace_lo) : nullptr;
    for (size_t j = lo; j <= hi; ++j) {
        const int32_t opened = match_up[j] + open;
        const int32_t extended = deletion_up[j] + extend;
        deletion[j] = std::min(g_no_cost, std::min(opened, extended));
        if (traced != nullptr) {
            traced[j - lo] = extended < opened ? g_deletion_extended : 0;
        }
    }
    if (params.long_gaps) {
        for (size_t j = lo; j <= hi; ++j) {
            const int32_t opened = match_up[j] + long_open;
            const int32_t extended = long_deletion_up[j] + long_extend;
            long_deletion[j] = std::min(g_no_cost, std::min(opened, extended));
            if (traced != nullptr) {
                traced[j - lo] |= extended < opened ? g_long_deletion_extended : 0;
            }
        }
    }
    // From the left: the insertions, the choice of source, and the trace.
    int32_t match_left = g_no_cost;
    int32_t insertion = g_no_cost;
    int32_t long_insertion = g_no_cost;
    for (size_t j = lo; j <= hi; ++j) {
        uint8_t bits = 0;
        const GapStep ins = gap_step(match_left, insertion, params.gap_open, params.gap_extend);
        insertion = ins.cost;
        bits |= ins.extended ? g_insertion_extended : 0;
        if (params.long_gaps) {
            const GapStep long_ins = gap_step(match_left, long_insertion, params.long_gap_open,
                                              params.long_gap_extend);
            long_insertion = long_ins.cost;
            bits |= long_ins.extended ? g_long_insertion_extended : 0;
        }
        int32_t best = match[j];
        uint8_t source = FromDiagonal;
        if (insertion < best) {
            best = insertion;
            source = FromInsertion;
        }
        if (long_insertion < best) {
            best = long_insertion;
            source = FromLongInsertion;
        }
        if (deletion[j] < best) {
            best = deletion[j];
            source = FromDeletion;
        }
        if (long_deletion[j] < best) {
            best = long_deletion[j];
            source = FromLongDeletion;
        }
        match[j] = best;
        match_left = best;
        if (traced != nullptr) {
            traced[j - lo] |= bits | source;
        }
    }
}

/**
 * Fill the band of diagonals `low` to `high` of the alignment of `a` to `b`.
 *
 * The diagonal of a cell is the number of bases of `b` before it minus the number of bases of
 * `a`. The band holds diagonal 0 and the diagonal of the last cell.
 *
 * @param trace If not null, set to one byte for each cell of each row of the band, which
 * `trace_back` reads.
 * @return The cost of the best alignment in the band.
 */
int32_t fill_band(const char* a, size_t n, const char* b, size_t m, int64_t low, int64_t high,
                  const StringAlignerParams& params, vector<uint8_t>* trace) {
    const size_t width = (size_t)(high - low + 1);
    if (trace != nullptr) {
        trace->assign((n + 1) * width, 0);
    }
    CellRow above(m);
    CellRow here(m);
    for (size_t i = 0; i <= n; ++i) {
        // The cells of the row that lie in the band and in the matrix.
        const int64_t row_lo = (int64_t)i + low;
        const auto lo = (size_t)std::max<int64_t>(0, row_lo);
        const auto hi = (size_t)std::min<int64_t>((int64_t)m, (int64_t)i + high);
        fill_row(a, b, i, lo, hi, params, above, here,
                 trace != nullptr ? trace->data() + i * width : nullptr, row_lo);
        std::swap(above, here);
    }
    return above.match[m + 1];
}

/// Trace the best alignment whose cells were traced in `trace` back from the last cell, where
/// `cell` gives the index in `trace` of cell (i, j), and return its runs, last first.
template<typename CellIndex>
vector<AlignmentRun> trace_cells(const vector<uint8_t>& trace, const char* a, size_t n,
                                 const char* b, size_t m, const CellIndex& cell) {
    vector<AlignmentRun> runs;
    // `state` is the state the alignment is in at cell (i, j).
    size_t i = n;
    size_t j = m;
    uint8_t state = FromDiagonal;
    while (i > 0 || j > 0) {
        const uint8_t traced = trace[cell(i, j)];
        if (state == FromDiagonal) {
            state = traced & SourceMask;
            if (state == FromDiagonal) {
                add_run(runs, a[i - 1] == b[j - 1] ? AlignmentOp::Match : AlignmentOp::Mismatch,
                        1);
                --i;
                --j;
            }
        } else if (state == FromInsertion || state == FromLongInsertion) {
            const uint8_t bit = state == FromInsertion ? g_insertion_extended
                                                       : g_long_insertion_extended;
            add_run(runs, AlignmentOp::Insertion, 1);
            --j;
            state = (traced & bit) ? state : (uint8_t)FromDiagonal;
        } else {
            const uint8_t bit = state == FromDeletion ? g_deletion_extended
                                                      : g_long_deletion_extended;
            add_run(runs, AlignmentOp::Deletion, 1);
            --i;
            state = (traced & bit) ? state : (uint8_t)FromDiagonal;
        }
    }
    return runs;
}

/// Trace the best alignment in the band whose cells `fill_band` traced in `trace` back from the
/// last cell, and return its runs, last first.
vector<AlignmentRun> trace_back(const vector<uint8_t>& trace, const char* a, size_t n,
                                const char* b, size_t m, int64_t low, int64_t high) {
    const size_t width = (size_t)(high - low + 1);
    return trace_cells(trace, a, n, b, m, [&](size_t i, size_t j) {
        return i * width + (size_t)((int64_t)j - (int64_t)i - low);
    });
}

/// The diagonals of a band.
struct Band {
    int64_t low;
    int64_t high;

    /// Return the number of cells the band has in a matrix of `n` + 1 rows, or the largest
    /// number if it does not fit one.
    size_t cells(size_t n) const {
        const auto width = (size_t)(high - low + 1);
        const size_t most = std::numeric_limits<size_t>::max();
        return width > most / (n + 1) ? most : width * (n + 1);
    }
};

/// Return the band of `reach` diagonals on each side of those between the two corners of the
/// matrix of `n` by `m` bases, cut to the diagonals the matrix has.
Band band_of(size_t n, size_t m, size_t reach) {
    const int64_t end = (int64_t)m - (int64_t)n;
    const auto side = (int64_t)std::min(reach, n + m);
    return Band{std::max(std::min<int64_t>(0, end) - side, -(int64_t)n),
                std::min(std::max<int64_t>(0, end) + side, (int64_t)m)};
}

/**
 * Return the least an alignment that leaves `band` can cost, or the largest number if the band
 * holds the whole matrix.
 *
 * An alignment that leaves the band above inserts enough bases to reach the first diagonal past
 * it, and deletes enough to come back to the last cell. One that leaves it below deletes and
 * then inserts. The gaps cost at least what one gap of each total length costs.
 */
int64_t cost_outside(const Band& band, size_t n, size_t m, const StringAlignerParams& params) {
    const int64_t end = (int64_t)m - (int64_t)n;
    int64_t outside = std::numeric_limits<int64_t>::max();
    if (band.high < (int64_t)m) {
        outside = std::min(outside, params.gap_cost((size_t)(band.high + 1))
                                        + params.gap_cost((size_t)(band.high + 1 - end)));
    }
    if (band.low > -(int64_t)n) {
        outside = std::min(outside, params.gap_cost((size_t)(1 - band.low))
                                        + params.gap_cost((size_t)(end - band.low + 1)));
    }
    return outside;
}

}

namespace {

/// The wavefront aligners of one thread: one for each set of costs it has aligned with.
struct WavefrontAligners {
    struct Entry {
        StringAlignerParams params;
        wavefront_aligner_t* aligner;
    };
    vector<Entry> entries;

    ~WavefrontAligners() {
        for (const Entry& entry : entries) {
            wavefront_aligner_delete(entry.aligner);
        }
    }

    /// Return the aligner with the costs of `params`.
    wavefront_aligner_t* with(const StringAlignerParams& params) {
        for (const Entry& entry : entries) {
            const StringAlignerParams& has = entry.params;
            if (has.mismatch == params.mismatch && has.gap_open == params.gap_open
                && has.gap_extend == params.gap_extend && has.long_gaps == params.long_gaps
                && (!has.long_gaps || (has.long_gap_open == params.long_gap_open
                                       && has.long_gap_extend == params.long_gap_extend))) {
                return entry.aligner;
            }
        }
        wavefront_aligner_attr_t attributes = wavefront_aligner_attr_default;
        attributes.alignment_scope = compute_alignment;
        attributes.alignment_form.span = alignment_end2end;
        attributes.heuristic.strategy = wf_heuristic_none;
        // The bidirectional aligner keeps memory in proportion to the cost of the alignment.
        attributes.memory_mode = wavefront_memory_ultralow;
        if (params.long_gaps) {
            attributes.distance_metric = gap_affine_2p;
            attributes.affine2p_penalties.match = 0;
            attributes.affine2p_penalties.mismatch = params.mismatch;
            attributes.affine2p_penalties.gap_opening1 = params.gap_open;
            attributes.affine2p_penalties.gap_extension1 = params.gap_extend;
            attributes.affine2p_penalties.gap_opening2 = params.long_gap_open;
            attributes.affine2p_penalties.gap_extension2 = params.long_gap_extend;
        } else {
            attributes.distance_metric = gap_affine;
            attributes.affine_penalties.match = 0;
            attributes.affine_penalties.mismatch = params.mismatch;
            attributes.affine_penalties.gap_opening = params.gap_open;
            attributes.affine_penalties.gap_extension = params.gap_extend;
        }
        entries.push_back({params, wavefront_aligner_new(&attributes)});
        return entries.back().aligner;
    }
};

/// Align all of `first` to all of `second` with the wavefront algorithm, which takes time and
/// memory in proportion to the cost of the alignment and not to the lengths of the strings.
optional<StringAlignment> align_wavefront(const string& first, const string& second,
                                          const StringAlignerParams& params) {
    StringAlignment alignment;
    if (first.empty() || second.empty()) {
        add_run(alignment.runs, AlignmentOp::Deletion, first.size());
        add_run(alignment.runs, AlignmentOp::Insertion, second.size());
        alignment.cost = params.gap_cost(first.size() + second.size());
        return alignment;
    }
    thread_local WavefrontAligners aligners;
    wavefront_aligner_t* aligner = aligners.with(params);
    // The aligner stops at a score. Its two directions can meet a little past the cost of the
    // alignment, by no more than the costs of opening gaps and of one step, so the score it
    // stops at is that much more than `max_cost`.
    const int64_t slack = 2 * ((int64_t)params.gap_open
                               + (params.long_gaps ? (int64_t)params.long_gap_open : 0)
                               + std::max<int64_t>(params.mismatch,
                                                   params.gap_open + params.gap_extend))
                          + 2;
    wavefront_aligner_set_max_alignment_steps(
        aligner, params.max_cost < 0
                     ? std::numeric_limits<int>::max()
                     : (int)std::min<int64_t>(std::numeric_limits<int>::max(),
                                              params.max_cost + slack));
    if (wavefront_align(aligner, first.data(), (int)first.size(), second.data(),
                        (int)second.size()) != WF_STATUS_ALG_COMPLETED) {
        return std::nullopt;
    }
    const cigar_t& cigar = *aligner->cigar;
    for (int i = cigar.begin_offset; i < cigar.end_offset; ++i) {
        const char op = cigar.operations[i];
        add_run(alignment.runs,
                op == 'M' ? AlignmentOp::Match
                          : op == 'X' ? AlignmentOp::Mismatch
                                      : op == 'I' ? AlignmentOp::Insertion
                                                  : AlignmentOp::Deletion,
                1);
    }
    for (const AlignmentRun& run : alignment.runs) {
        if (run.op == AlignmentOp::Mismatch) {
            alignment.cost += (int64_t)params.mismatch * (int64_t)run.length;
        } else if (run.op != AlignmentOp::Match) {
            alignment.cost += params.gap_cost(run.length);
        }
    }
    return alignment;
}

}

int64_t StringAlignerParams::gap_cost(size_t length) const {
    if (length == 0) {
        return 0;
    }
    const int64_t cost = gap_open + (int64_t)gap_extend * (int64_t)length;
    if (!long_gaps) {
        return cost;
    }
    return std::min(cost, long_gap_open + (int64_t)long_gap_extend * (int64_t)length);
}

optional<StringAlignment> align_global(const string& first, const string& second,
                                       const StringAlignerParams& params) {
    if (params.wavefront) {
        return align_wavefront(first, second, params);
    }
    // Some best alignment matches the bases the strings share at their starts, and then those
    // they share at their ends, so only what lies between needs a band.
    size_t prefix = 0;
    while (prefix < first.size() && prefix < second.size() && first[prefix] == second[prefix]) {
        ++prefix;
    }
    size_t suffix = 0;
    while (suffix < first.size() - prefix && suffix < second.size() - prefix
           && first[first.size() - 1 - suffix] == second[second.size() - 1 - suffix]) {
        ++suffix;
    }
    const size_t n = first.size() - prefix - suffix;
    const size_t m = second.size() - prefix - suffix;

    StringAlignment result;
    add_run(result.runs, AlignmentOp::Match, prefix);
    if (n == 0 || m == 0) {
        add_run(result.runs, n == 0 ? AlignmentOp::Insertion : AlignmentOp::Deletion,
                std::max(n, m));
        result.cost = params.gap_cost(std::max(n, m));
        add_run(result.runs, AlignmentOp::Match, suffix);
        return result;
    }

    const char* a = first.data() + prefix;
    const char* b = second.data() + prefix;
    vector<uint8_t> trace;
    // The least cost of any alignment, once a band wide enough to show it has been filled.
    optional<int64_t> least;
    int64_t narrower_cost = -1;
    for (size_t reach = std::max<size_t>(params.first_band, 1);; reach *= 2) {
        const Band band = band_of(n, m, reach);
        if (band.cells(n) > params.max_cells) {
            return std::nullopt;
        }
        const int64_t cost = fill_band(a, n, b, m, band.low, band.high, params, &trace);
        result.cells_filled += band.cells(n);
        result.cells_traced = std::max(result.cells_traced, band.cells(n));
        if (!least && cost > cost_outside(band, n, m, params) && cost == narrower_cost) {
            // Doubling the band found nothing cheaper, so this is likely the least cost. Fill,
            // without tracing, the narrowest band that an alignment of this cost cannot leave.
            // The best alignment in that band is the best of all.
            size_t enough = reach;
            size_t too_few = reach;
            while (cost_outside(band_of(n, m, enough), n, m, params) < cost) {
                too_few = enough;
                enough *= 2;
            }
            while (enough - too_few > 1) {
                const size_t middle = too_few + (enough - too_few) / 2;
                (cost_outside(band_of(n, m, middle), n, m, params) < cost ? too_few : enough) =
                    middle;
            }
            const Band proof = band_of(n, m, enough);
            if (proof.cells(n) > params.max_work) {
                return std::nullopt;
            }
            least = fill_band(a, n, b, m, proof.low, proof.high, params, nullptr);
            result.cells_filled += proof.cells(n);
            result.cells_proof += proof.cells(n);
        }
        if (cost <= cost_outside(band, n, m, params) || (least && cost == *least)) {
            const vector<AlignmentRun> middle = trace_back(trace, a, n, b, m, band.low,
                                                           band.high);
            for (size_t r = middle.size(); r-- > 0;) {
                add_run(result.runs, middle[r].op, middle[r].length);
            }
            add_run(result.runs, AlignmentOp::Match, suffix);
            result.cost = cost;
            result.band = reach;
            return result;
        }
        narrower_cost = cost;
    }
}

optional<StringAlignment> align_guided(const string& first, const string& second,
                                       const vector<pair<size_t, size_t>>& corners,
                                       size_t margin, const StringAlignerParams& params) {
    // As `align_global` does, match the bases the two share at their starts and then at their
    // ends, and align only what lies between.
    size_t prefix = 0;
    while (prefix < first.size() && prefix < second.size() && first[prefix] == second[prefix]) {
        ++prefix;
    }
    size_t suffix = 0;
    while (suffix < first.size() - prefix && suffix < second.size() - prefix
           && first[first.size() - 1 - suffix] == second[second.size() - 1 - suffix]) {
        ++suffix;
    }
    const size_t n = first.size() - prefix - suffix;
    const size_t m = second.size() - prefix - suffix;
    StringAlignment result;
    add_run(result.runs, AlignmentOp::Match, prefix);
    if (n == 0 || m == 0) {
        add_run(result.runs, n == 0 ? AlignmentOp::Insertion : AlignmentOp::Deletion,
                std::max(n, m));
        result.cost = params.gap_cost(std::max(n, m));
        add_run(result.runs, AlignmentOp::Match, suffix);
        return result;
    }
    // The corners in what lies between, from its first cell to its last, never going back in
    // either string.
    vector<pair<size_t, size_t>> at{{0, 0}};
    for (const auto& [i, j] : corners) {
        at.emplace_back(std::clamp(i > prefix ? i - prefix : 0, at.back().first, n),
                        std::clamp(j > prefix ? j - prefix : 0, at.back().second, m));
    }
    at.emplace_back(n, m);
    // Each row's cells are the columns that the bands of the corners hold in it. The band of two
    // corners that follow each other holds the rows from the one to the other, and the diagonals
    // from the one's to the other's, each widened by `margin`. The bands of the corners a row is
    // in meet, so its cells are one interval.
    vector<size_t> lo(n + 1, m);
    vector<size_t> hi(n + 1, 0);
    for (size_t e = 0; e + 1 < at.size(); ++e) {
        const int64_t from = (int64_t)at[e].second - (int64_t)at[e].first;
        const int64_t to = (int64_t)at[e + 1].second - (int64_t)at[e + 1].first;
        const int64_t low = std::min(from, to) - (int64_t)margin;
        const int64_t high = std::max(from, to) + (int64_t)margin;
        const size_t top = at[e].first > margin ? at[e].first - margin : 0;
        const size_t bottom = std::min(n, at[e + 1].first + margin);
        for (size_t i = top; i <= bottom; ++i) {
            lo[i] = std::min(lo[i], (size_t)std::clamp<int64_t>((int64_t)i + low, 0, (int64_t)m));
            hi[i] = std::max(hi[i], (size_t)std::clamp<int64_t>((int64_t)i + high, 0, (int64_t)m));
        }
    }
    // Where each row's cells start in the trace.
    vector<size_t> start(n + 2, 0);
    for (size_t i = 0; i <= n; ++i) {
        start[i + 1] = start[i] + (hi[i] - lo[i] + 1);
    }
    if (start[n + 1] > params.max_cells) {
        return std::nullopt;
    }
    const char* a = first.data() + prefix;
    const char* b = second.data() + prefix;
    vector<uint8_t> trace(start[n + 1]);
    CellRow above(m);
    CellRow here(m);
    for (size_t i = 0; i <= n; ++i) {
        fill_row(a, b, i, lo[i], hi[i], params, above, here, trace.data() + start[i],
                 (int64_t)lo[i]);
        std::swap(above, here);
    }
    result.cost = above.match[m + 1];
    const vector<AlignmentRun> runs = trace_cells(trace, a, n, b, m, [&](size_t i, size_t j) {
        return start[i] + (j - lo[i]);
    });
    for (size_t r = runs.size(); r-- > 0;) {
        add_run(result.runs, runs[r].op, runs[r].length);
    }
    add_run(result.runs, AlignmentOp::Match, suffix);
    result.cells_filled = start[n + 1];
    result.cells_traced = start[n + 1];
    return result;
}

void left_align_gaps(const string& first, const string& second, vector<AlignmentRun>& runs,
                     bool keep_first) {
    // One operation for each column of the alignment.
    vector<AlignmentOp> columns;
    for (const AlignmentRun& run : runs) {
        columns.insert(columns.end(), run.length, run.op);
    }
    const size_t floor = keep_first ? 1 : 0;
    // The numbers of bases of the two strings before the column.
    size_t i = 0;
    size_t j = 0;
    for (size_t c = 0; c < columns.size();) {
        const AlignmentOp op = columns[c];
        size_t end = c + 1;
        while (end < columns.size() && columns[end] == op) {
            ++end;
        }
        const size_t length = end - c;
        if (op == AlignmentOp::Match || op == AlignmentOp::Mismatch) {
            i += length;
            j += length;
            c = end;
            continue;
        }
        // Move the gap left while the base before it is matched and ends the gap as well.
        const string& gapped = op == AlignmentOp::Insertion ? second : first;
        size_t at = op == AlignmentOp::Insertion ? j : i;
        size_t start = c;
        while (start > floor && columns[start - 1] == AlignmentOp::Match
               && gapped[at - 1] == gapped[at + length - 1]) {
            columns[start - 1] = op;
            columns[start + length - 1] = AlignmentOp::Match;
            --start;
            --at;
        }
        // The matched bases the gap moved past now come after it, still before column `end`.
        (op == AlignmentOp::Insertion ? j : i) += length;
        c = end;
    }
    runs.clear();
    for (const AlignmentOp op : columns) {
        add_run(runs, op, 1);
    }
}

}
