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
    // The costs of the row above and of this row, by diagonal. Each has one cell past the band,
    // which no alignment reaches.
    vector<int32_t> match_above(width + 1, g_no_cost);
    vector<int32_t> match_here(width + 1, g_no_cost);
    vector<int32_t> deletion_above(width + 1, g_no_cost);
    vector<int32_t> deletion_here(width + 1, g_no_cost);
    vector<int32_t> long_deletion_above(width + 1, g_no_cost);
    vector<int32_t> long_deletion_here(width + 1, g_no_cost);

    for (size_t i = 0; i <= n; ++i) {
        std::ranges::fill(match_here, g_no_cost);
        std::ranges::fill(deletion_here, g_no_cost);
        std::ranges::fill(long_deletion_here, g_no_cost);
        // The cells of the row that lie in the band and in the matrix.
        const int64_t first = std::max<int64_t>(0, -(int64_t)i - low);
        const int64_t last = std::min<int64_t>((int64_t)width - 1, (int64_t)m - (int64_t)i - low);
        // The costs of the cell to the left.
        int32_t match_left = g_no_cost;
        int32_t insertion = g_no_cost;
        int32_t long_insertion = g_no_cost;
        uint8_t* trace_row = trace != nullptr ? trace->data() + i * width : nullptr;
        for (int64_t k = first; k <= last; ++k) {
            const size_t j = (size_t)((int64_t)i + low + k);
            uint8_t traced = 0;
            // An insertion comes from the cell to the left.
            const GapStep ins = gap_step(match_left, insertion, params.gap_open,
                                         params.gap_extend);
            insertion = ins.cost;
            traced |= ins.extended ? g_insertion_extended : 0;
            if (params.long_gaps) {
                const GapStep long_ins = gap_step(match_left, long_insertion,
                                                  params.long_gap_open, params.long_gap_extend);
                long_insertion = long_ins.cost;
                traced |= long_ins.extended ? g_long_insertion_extended : 0;
            }
            // A deletion comes from the cell above, which is one diagonal up.
            const GapStep del = gap_step(match_above[k + 1], deletion_above[k + 1],
                                         params.gap_open, params.gap_extend);
            deletion_here[k] = del.cost;
            traced |= del.extended ? g_deletion_extended : 0;
            if (params.long_gaps) {
                const GapStep long_del = gap_step(match_above[k + 1], long_deletion_above[k + 1],
                                                  params.long_gap_open, params.long_gap_extend);
                long_deletion_here[k] = long_del.cost;
                traced |= long_del.extended ? g_long_deletion_extended : 0;
            }
            // A match or mismatch comes from the cell above and to the left, on this diagonal.
            int32_t best = g_no_cost;
            if (i == 0 && j == 0) {
                best = 0;
            } else if (i > 0 && j > 0) {
                best = std::min(g_no_cost,
                                match_above[k] + (a[i - 1] == b[j - 1] ? 0 : params.mismatch));
            }
            uint8_t source = FromDiagonal;
            if (insertion < best) {
                best = insertion;
                source = FromInsertion;
            }
            if (long_insertion < best) {
                best = long_insertion;
                source = FromLongInsertion;
            }
            if (deletion_here[k] < best) {
                best = deletion_here[k];
                source = FromDeletion;
            }
            if (long_deletion_here[k] < best) {
                best = long_deletion_here[k];
                source = FromLongDeletion;
            }
            match_here[k] = best;
            match_left = best;
            if (trace_row != nullptr) {
                trace_row[k] = traced | source;
            }
        }
        match_above.swap(match_here);
        deletion_above.swap(deletion_here);
        long_deletion_above.swap(long_deletion_here);
    }
    return match_above[(size_t)((int64_t)m - (int64_t)n - low)];
}

/// Trace the best alignment in the band whose cells `fill_band` traced in `trace` back from the
/// last cell, and return its runs, last first.
vector<AlignmentRun> trace_back(const vector<uint8_t>& trace, const char* a, size_t n,
                                const char* b, size_t m, int64_t low, int64_t high) {
    const size_t width = (size_t)(high - low + 1);
    vector<AlignmentRun> runs;
    // `state` is the state the alignment is in at cell (i, j).
    size_t i = n;
    size_t j = m;
    uint8_t state = FromDiagonal;
    while (i > 0 || j > 0) {
        const uint8_t traced = trace[i * width + (size_t)((int64_t)j - (int64_t)i - low)];
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
