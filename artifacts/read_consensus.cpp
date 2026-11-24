/** \file
 * Implements the consensus of reads that span one window.
 */

#include "read_consensus.hpp"

#include <algorithm>
#include <array>
#include <cstdlib>
#include <limits>
#include <numeric>

namespace vg {

namespace {

/// The votes a base of the backbone can get: the four bases, any other character, a deletion.
constexpr size_t VOTES = 6;
constexpr size_t DELETION_VOTE = 5;

/// Return the vote that a read with `base` casts.
size_t vote_of(char base) {
    switch (base) {
    case 'A': return 0;
    case 'C': return 1;
    case 'G': return 2;
    case 'T': return 3;
    default: return 4;
    }
}

/// Return the base that vote `vote` stands for.
char base_of(size_t vote) {
    return "ACGTN"[vote];
}

/// Return the indexes of at most `limit` of `strings`, those whose lengths are nearest the
/// median length, in the order of the vector.
vector<size_t> nearest_median_length(const vector<const string*>& strings, size_t limit) {
    vector<size_t> order(strings.size());
    std::iota(order.begin(), order.end(), 0);
    if (strings.size() <= limit) {
        return order;
    }
    vector<size_t> lengths;
    lengths.reserve(strings.size());
    for (const string* s : strings) {
        lengths.push_back(s->size());
    }
    std::ranges::nth_element(lengths, lengths.begin() + (ptrdiff_t)(lengths.size() / 2));
    const auto median = (int64_t)lengths[lengths.size() / 2];
    std::ranges::stable_sort(order, [&](size_t a, size_t b) {
        return std::llabs((int64_t)strings[a]->size() - median)
               < std::llabs((int64_t)strings[b]->size() - median);
    });
    order.resize(limit);
    std::ranges::sort(order);
    return order;
}

/// The medoid of some strings, and the sum of its edit distances to the others.
struct Medoid {
    size_t index = 0;
    int64_t distance = 0;
};

/// Return the first of `strings` with the least summed edit distance to the others, or nothing
/// if an edit distance is not found within `max_cells`. A string that more than half of them
/// equal is the medoid, and is found without aligning.
optional<Medoid> medoid_of(const vector<const string*>& strings, size_t max_cells,
                           bool need_distance) {
    if (!need_distance) {
        for (size_t i = 0; i < strings.size(); ++i) {
            size_t same = 0;
            for (const string* other : strings) {
                same += *other == *strings[i];
            }
            if (same * 2 > strings.size()) {
                return Medoid{i, 0};
            }
        }
    }
    vector<int64_t> sums(strings.size(), 0);
    for (size_t i = 0; i < strings.size(); ++i) {
        for (size_t j = i + 1; j < strings.size(); ++j) {
            const optional<int64_t> distance = edit_distance(*strings[i], *strings[j], max_cells);
            if (!distance) {
                return std::nullopt;
            }
            sums[i] += *distance;
            sums[j] += *distance;
        }
    }
    const size_t best = (size_t)(std::ranges::min_element(sums) - sums.begin());
    return Medoid{best, sums[best]};
}

/// Return the result of one round: `backbone` with each base and each gap between bases set
/// to what most of `reads` have there. Nothing if an alignment is not found.
optional<string> vote(const vector<const string*>& reads, const string& backbone,
                      const ReadConsensusParams& params) {
    vector<std::array<uint32_t, VOTES>> columns(backbone.size(), std::array<uint32_t, VOTES>{});
    // The strings the reads insert before each base of the backbone, and after the last.
    vector<vector<string>> inserted(backbone.size() + 1);
    for (const string* read : reads) {
        optional<StringAlignment> alignment = align_global(backbone, *read, params.aligner);
        if (!alignment) {
            return std::nullopt;
        }
        left_align_gaps(backbone, *read, alignment->runs);
        size_t i = 0;
        size_t j = 0;
        for (const AlignmentRun& run : alignment->runs) {
            switch (run.op) {
            case AlignmentOp::Match:
            case AlignmentOp::Mismatch:
                for (size_t k = 0; k < run.length; ++k) {
                    ++columns[i++][vote_of((*read)[j++])];
                }
                break;
            case AlignmentOp::Deletion:
                for (size_t k = 0; k < run.length; ++k) {
                    ++columns[i++][DELETION_VOTE];
                }
                break;
            case AlignmentOp::Insertion:
                inserted[i].push_back(read->substr(j, run.length));
                j += run.length;
                break;
            }
        }
    }

    string result;
    result.reserve(backbone.size());
    for (size_t i = 0; i <= backbone.size(); ++i) {
        if (inserted[i].size() * 2 > reads.size()) {
            vector<const string*> all;
            for (const string& s : inserted[i]) {
                all.push_back(&s);
            }
            vector<const string*> compared;
            for (size_t index : nearest_median_length(all, params.max_medoid)) {
                compared.push_back(all[index]);
            }
            const optional<Medoid> medoid = medoid_of(compared, params.aligner.max_cells, false);
            if (!medoid) {
                return std::nullopt;
            }
            result += *compared[medoid->index];
        }
        if (i == backbone.size()) {
            break;
        }
        const size_t own = vote_of(backbone[i]);
        size_t best = own;
        for (size_t v = 0; v < VOTES; ++v) {
            if (columns[i][v] > columns[i][best]) {
                best = v;
            }
        }
        if (best == own) {
            result += backbone[i];
        } else if (best != DELETION_VOTE) {
            result += base_of(best);
        }
    }
    return result;
}

/// Return true if every base of `bases` is `base`.
bool all_are(const string& bases, size_t from, size_t length, char base) {
    for (size_t k = 0; k < length; ++k) {
        if (bases[from + k] != base) {
            return false;
        }
    }
    return true;
}

}

optional<int64_t> edit_distance(const string& first, const string& second, size_t max_cells) {
    StringAlignerParams unit;
    unit.mismatch = 1;
    unit.gap_open = 0;
    unit.gap_extend = 1;
    unit.long_gaps = false;
    unit.max_cells = max_cells;
    const optional<StringAlignment> alignment = align_global(first, second, unit);
    if (!alignment) {
        return std::nullopt;
    }
    return alignment->cost;
}

optional<string> keep_homopolymer_lengths(const string& sequence, const string& backbone,
                                          const StringAlignerParams& aligner, size_t min_run) {
    optional<StringAlignment> alignment = align_global(backbone, sequence, aligner);
    if (!alignment) {
        return std::nullopt;
    }
    vector<AlignmentRun>& runs = alignment->runs;
    left_align_gaps(backbone, sequence, runs);
    string result;
    result.reserve(sequence.size());
    size_t i = 0;
    size_t j = 0;
    for (size_t r = 0; r < runs.size(); ++r) {
        const AlignmentRun& run = runs[r];
        if (run.op == AlignmentOp::Match || run.op == AlignmentOp::Mismatch) {
            result.append(sequence, j, run.length);
            i += run.length;
            j += run.length;
            continue;
        }
        const bool insertion = run.op == AlignmentOp::Insertion;
        const char base = insertion ? sequence[j] : backbone[i];
        // The gap changes the length of a run if its bases are all one base and the matched
        // base before or after it is that base too.
        bool in_run = all_are(insertion ? sequence : backbone, insertion ? j : i, run.length, base);
        if (in_run) {
            const bool before = r > 0 && runs[r - 1].op == AlignmentOp::Match
                                && backbone[i - 1] == base;
            const size_t next = insertion ? i : i + run.length;
            const bool after = r + 1 < runs.size() && runs[r + 1].op == AlignmentOp::Match
                               && backbone[next] == base;
            in_run = before || after;
        }
        if (in_run) {
            // The run in the backbone that the gap lengthens or shortens.
            size_t first = i;
            while (first > 0 && backbone[first - 1] == base) {
                --first;
            }
            size_t past = i;
            while (past < backbone.size() && backbone[past] == base) {
                ++past;
            }
            in_run = past - first >= min_run;
        }
        if (insertion) {
            if (!in_run) {
                result.append(sequence, j, run.length);
            }
            j += run.length;
        } else {
            if (in_run) {
                result.append(backbone, i, run.length);
            }
            i += run.length;
        }
    }
    return result;
}

optional<string> read_consensus(const vector<string>& reads, const string& backbone,
                                const ReadConsensusParams& params) {
    if (reads.size() < std::max<size_t>(params.min_reads, 1)) {
        return std::nullopt;
    }
    vector<const string*> all;
    all.reserve(reads.size());
    for (const string& read : reads) {
        all.push_back(&read);
    }
    vector<const string*> used;
    for (size_t index : nearest_median_length(all, params.max_reads)) {
        used.push_back(all[index]);
    }

    string current = backbone;
    if (params.seed_from_reads) {
        // Start from the medoid read where the backbone is further from the compared reads
        // than it is.
        vector<const string*> compared;
        for (size_t index : nearest_median_length(used, params.max_medoid)) {
            compared.push_back(used[index]);
        }
        const optional<Medoid> medoid = medoid_of(compared, params.aligner.max_cells, true);
        if (!medoid) {
            return std::nullopt;
        }
        int64_t backbone_distance = 0;
        for (size_t k = 0; k < compared.size(); ++k) {
            if (k == medoid->index) {
                continue;
            }
            const optional<int64_t> distance = edit_distance(backbone, *compared[k],
                                                             params.aligner.max_cells);
            if (!distance) {
                return std::nullopt;
            }
            backbone_distance += *distance;
        }
        if (backbone_distance > medoid->distance) {
            current = *compared[medoid->index];
        }
    }

    for (size_t round = 0; round < params.max_rounds; ++round) {
        optional<string> next = vote(used, current, params);
        if (!next) {
            return std::nullopt;
        }
        if (*next == current) {
            break;
        }
        current = std::move(*next);
    }
    if (params.keep_homopolymer_lengths) {
        return keep_homopolymer_lengths(current, backbone, params.aligner,
                                        params.homopolymer_min_run);
    }
    return current;
}

}
