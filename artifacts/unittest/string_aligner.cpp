/// \file string_aligner.cpp
///
/// Unit tests for the global string aligner and for moving gaps left.

#include <random>

#include "../../src/unittest/catch.hpp"
#include "../string_aligner.hpp"

namespace vg {
namespace unittest {

using Op = AlignmentOp;

/// Return `length` random bases.
static string random_bases(size_t length, std::mt19937& generator) {
    string bases;
    for (size_t i = 0; i < length; ++i) {
        bases += "ACGT"[generator() % 4];
    }
    return bases;
}

/// Return the two strings that `runs` aligns, given the first, taking the bases of the second
/// from `second`, and require that the runs cover both and that each run is what it says.
static void require_spells(const string& first, const string& second,
                           const vector<AlignmentRun>& runs) {
    size_t i = 0;
    size_t j = 0;
    for (const AlignmentRun& run : runs) {
        REQUIRE(run.length > 0);
        for (size_t k = 0; k < run.length; ++k) {
            if (run.op == Op::Match || run.op == Op::Mismatch) {
                REQUIRE(i < first.size());
                REQUIRE(j < second.size());
                REQUIRE((first[i] == second[j]) == (run.op == Op::Match));
                ++i;
                ++j;
            } else if (run.op == Op::Insertion) {
                ++j;
            } else {
                ++i;
            }
        }
    }
    REQUIRE(i == first.size());
    REQUIRE(j == second.size());
}

/// Return the cost of the alignment `runs` under `params`.
static int64_t alignment_cost(const vector<AlignmentRun>& runs,
                              const StringAlignerParams& params) {
    int64_t cost = 0;
    for (const AlignmentRun& run : runs) {
        if (run.op == Op::Mismatch) {
            cost += (int64_t)params.mismatch * (int64_t)run.length;
        } else if (run.op != Op::Match) {
            cost += params.gap_cost(run.length);
        }
    }
    return cost;
}

TEST_CASE("The string aligner aligns identical and empty strings", "[string_aligner]") {
    const optional<StringAlignment> same = align_global("GATTACA", "GATTACA");
    REQUIRE(same);
    REQUIRE(same->cost == 0);
    REQUIRE(same->runs == vector<AlignmentRun>{{Op::Match, 7}});

    const optional<StringAlignment> nothing = align_global("", "");
    REQUIRE(nothing);
    REQUIRE(nothing->cost == 0);
    REQUIRE(nothing->runs.empty());

    const optional<StringAlignment> inserted = align_global("", "ACG");
    REQUIRE(inserted);
    REQUIRE(inserted->cost == 12);
    REQUIRE(inserted->runs == vector<AlignmentRun>{{Op::Insertion, 3}});

    const optional<StringAlignment> deleted = align_global("ACG", "");
    REQUIRE(deleted);
    REQUIRE(deleted->cost == 12);
    REQUIRE(deleted->runs == vector<AlignmentRun>{{Op::Deletion, 3}});
}

TEST_CASE("The string aligner finds a mismatch", "[string_aligner]") {
    const optional<StringAlignment> snv = align_global("GATTACAGATTACA", "GATTACCGATTACA");
    REQUIRE(snv);
    REQUIRE(snv->cost == 4);
    REQUIRE(snv->runs
            == vector<AlignmentRun>{{Op::Match, 6}, {Op::Mismatch, 1}, {Op::Match, 7}});
}

TEST_CASE("An insertion in a tandem repeat is placed at the start of the repeat",
          "[string_aligner]") {
    const string first = "GATTACACAGCAGCAGCAGCAGTTGCA";
    const string second = "GATTACACAGCAGCAGCAGCAGCAGTTGCA";
    optional<StringAlignment> aligned = align_global(first, second);
    REQUIRE(aligned);
    REQUIRE(aligned->cost == 12);
    left_align_gaps(first, second, aligned->runs);
    REQUIRE(aligned->runs
            == vector<AlignmentRun>{{Op::Match, 7}, {Op::Insertion, 3}, {Op::Match, 20}});
    require_spells(first, second, aligned->runs);

    // The same deletion, read the other way, goes to the same place.
    optional<StringAlignment> back = align_global(second, first);
    REQUIRE(back);
    left_align_gaps(second, first, back->runs);
    REQUIRE(back->runs
            == vector<AlignmentRun>{{Op::Match, 7}, {Op::Deletion, 3}, {Op::Match, 20}});

    // A gap that could move to the start of the strings can be made to keep one base before it.
    vector<AlignmentRun> runs{{Op::Match, 4}, {Op::Insertion, 1}};
    left_align_gaps("AAAA", "AAAAA", runs);
    REQUIRE(runs == vector<AlignmentRun>{{Op::Insertion, 1}, {Op::Match, 4}});
    runs = {{Op::Match, 4}, {Op::Insertion, 1}};
    left_align_gaps("AAAA", "AAAAA", runs, true);
    REQUIRE(runs == vector<AlignmentRun>{{Op::Match, 1}, {Op::Insertion, 1}, {Op::Match, 3}});
}

TEST_CASE("A long insertion is kept as one gap", "[string_aligner]") {
    std::mt19937 generator(7);
    const string left = random_bases(300, generator);
    const string right = random_bases(300, generator);
    const string inserted = random_bases(500, generator);
    const string first = left + right;
    const string second = left + inserted + right;
    optional<StringAlignment> aligned = align_global(first, second);
    REQUIRE(aligned);
    // The second piece of the gap cost is the cheaper for a gap this long.
    REQUIRE(aligned->cost == 24 + 500);
    left_align_gaps(first, second, aligned->runs);
    REQUIRE(aligned->runs.size() == 3);
    REQUIRE(aligned->runs[1].op == Op::Insertion);
    REQUIRE(aligned->runs[1].length == 500);
    require_spells(first, second, aligned->runs);
}

TEST_CASE("The string aligner widens its band until the alignment is the best",
          "[string_aligner]") {
    std::mt19937 generator(11);
    // The first string has 30 bases at its start that the second lacks, and the second has 30
    // at its end, so the best alignment runs 30 diagonals from the corners' diagonal.
    const string shared = random_bases(400, generator);
    const string first = random_bases(30, generator) + shared;
    const string second = shared + random_bases(30, generator);
    StringAlignerParams narrow;
    narrow.first_band = 2;
    const optional<StringAlignment> aligned = align_global(first, second, narrow);
    REQUIRE(aligned);
    REQUIRE(aligned->band > 30);
    REQUIRE(aligned->cost == 2 * (24 + 30));
    REQUIRE(alignment_cost(aligned->runs, narrow) == aligned->cost);
    require_spells(first, second, aligned->runs);

    // A band that holds the whole matrix finds the same cost.
    StringAlignerParams whole;
    whole.first_band = 1000;
    const optional<StringAlignment> full = align_global(first, second, whole);
    REQUIRE(full);
    REQUIRE(full->cost == aligned->cost);

    // With too few cells allowed, the aligner says so.
    StringAlignerParams capped;
    capped.first_band = 2;
    capped.max_cells = 431 * 20;
    REQUIRE(!align_global(first, second, capped));
}

TEST_CASE("The string aligner's banded alignments cost what a whole matrix finds",
          "[string_aligner]") {
    std::mt19937 generator(3);
    for (size_t trial = 0; trial < 40; ++trial) {
        // A string and a copy of it with mismatches, insertions and deletions.
        const string first = random_bases(200 + generator() % 200, generator);
        string second;
        for (size_t i = 0; i < first.size(); ++i) {
            const size_t roll = generator() % 60;
            if (roll == 0) {
                second += random_bases(1 + generator() % 40, generator);
            } else if (roll == 1) {
                i += generator() % 40;
                continue;
            }
            second += roll == 2 ? "ACGT"[generator() % 4] : first[i];
        }
        StringAlignerParams banded;
        banded.first_band = 1;
        StringAlignerParams whole;
        whole.first_band = 100000;
        optional<StringAlignment> narrow = align_global(first, second, banded);
        const optional<StringAlignment> full = align_global(first, second, whole);
        REQUIRE(narrow);
        REQUIRE(full);
        REQUIRE(narrow->cost == full->cost);
        REQUIRE(alignment_cost(narrow->runs, banded) == narrow->cost);
        require_spells(first, second, narrow->runs);
        left_align_gaps(first, second, narrow->runs, true);
        require_spells(first, second, narrow->runs);
        REQUIRE(alignment_cost(narrow->runs, banded) == narrow->cost);
    }
}

TEST_CASE("The string aligner aligns long strings with many differences", "[string_aligner]") {
    std::mt19937 generator(5);
    const string first = random_bases(25000, generator);
    string second;
    size_t changes = 0;
    for (size_t i = 0; i < first.size(); ++i) {
        const size_t roll = generator() % 40;
        if (roll == 0) {
            second += "ACGT"[generator() % 4];
            ++changes;
        } else if (roll == 1) {
            second += first[i];
            second += random_bases(1 + generator() % 3, generator);
            ++changes;
        } else if (roll != 2) {
            second += first[i];
        }
    }
    REQUIRE(changes > 1000);
    const optional<StringAlignment> aligned = align_global(first, second);
    REQUIRE(aligned);
    REQUIRE(alignment_cost(aligned->runs, StringAlignerParams()) == aligned->cost);
    require_spells(first, second, aligned->runs);
}

}

TEST_CASE("The wavefront aligner gives up an alignment that costs more than a bound",
          "[string_aligner]") {
    StringAlignerParams params;
    params.wavefront = true;
    const string first(400, 'A');
    const string second = string(200, 'A') + string(300, 'C') + string(200, 'A');
    const optional<StringAlignment> whole = align_global(first, second, params);
    REQUIRE(whole);
    REQUIRE(whole->cost == params.gap_cost(300));
    params.max_cost = whole->cost / 2;
    REQUIRE(!align_global(first, second, params));
    // The aligner aligns again after it gave up, and within a bound that the cost meets.
    params.max_cost = whole->cost;
    const optional<StringAlignment> bounded = align_global(first, second, params);
    REQUIRE(bounded);
    REQUIRE(bounded->cost == whole->cost);
    params.max_cost = -1;
    REQUIRE(align_global(first, first, params)->cost == 0);
}
}
