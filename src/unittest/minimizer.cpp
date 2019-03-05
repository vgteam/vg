/** \file
 *
 * Unit tests for MinimizerIndex, which is basically a hash table mapping kmers
 * to sets of pos_t.
 */

#include "../minimizer.hpp"
#include "../json2pb.h"
#include "../position.hpp"
#include "../utility.hpp"

#include "catch.hpp"

#include <map>
#include <set>
#include <vector>

namespace vg {
namespace unittest {

TEST_CASE("MinimizerIndex construction, assignment, and serialization", "[minimizer_index][indexing]") {

    SECTION("equality of empty indexes depends on parameter values") {
        MinimizerIndex default_index;
        MinimizerIndex default_copy(default_index);
        MinimizerIndex alt_index(15, 6);
        MinimizerIndex alt_copy(alt_index);
        REQUIRE(default_index == default_copy);
        REQUIRE(alt_index == alt_copy);
        REQUIRE(default_index != alt_index);
    }

    SECTION("contents make indexes different") {
        MinimizerIndex default_index;
        MinimizerIndex default_copy(default_index);

        // Different contents.
        default_index.insert(1, MinimizerIndex::make_pos_t(1, false, 3));
        REQUIRE(default_index != default_copy);

        // Same key, different value.
        default_copy.insert(1, MinimizerIndex::make_pos_t(2, false, 3));
        REQUIRE(default_index != default_copy);

        // Same contents.
        default_copy = default_index;
        REQUIRE(default_index == default_copy);
    }

    SECTION("index contents can be swapped") {
        MinimizerIndex first, second;
        first.insert(1, MinimizerIndex::make_pos_t(1, false, 3));
        second.insert(2, MinimizerIndex::make_pos_t(2, false, 3));

        MinimizerIndex first_copy(first), second_copy(second);
        first.swap(second);
        REQUIRE(first != first_copy);
        REQUIRE(first == second_copy);
        REQUIRE(second == first_copy);
        REQUIRE(second != second_copy);
    }

    SECTION("serialization preserves parameters and contents") {
        MinimizerIndex index(15, 6);
        index.insert(1, MinimizerIndex::make_pos_t(1, false, 3));
        index.insert(2, MinimizerIndex::make_pos_t(1, false, 3));
        index.insert(2, MinimizerIndex::make_pos_t(2, false, 3));

        std::string filename = temp_file::create("minimizer");
        std::ofstream out(filename, std::ios_base::binary);
        index.serialize(out);
        out.close();

        MinimizerIndex copy;
        std::ifstream in(filename, std::ios_base::binary);
        copy.load(in);
        in.close();
        temp_file::remove(filename);

        REQUIRE(index == copy);
    }
}

TEST_CASE("Minimizer extraction works correctly", "[minimizer_index][indexing]") {
    std::string str = "CGAATACAATACT";

    SECTION("minimizer is the leftmost occurrence") {
        MinimizerIndex index(5, 2);
        MinimizerIndex::minimizer_type correct(0 * 256 + 0 * 64 + 3 * 16 + 0 * 4 + 1 * 1, 2);
        MinimizerIndex::minimizer_type result = index.minimizer(str.begin(), str.end());
        REQUIRE(result == correct);
    }

    SECTION("all minimizers are found") {
        MinimizerIndex index(5, 2);
        std::vector<MinimizerIndex::minimizer_type> correct {
            MinimizerIndex::minimizer_type(1 * 256 + 2 * 64 + 0 * 16 + 0 * 4 + 3 * 1, 0),
            MinimizerIndex::minimizer_type(0 * 256 + 0 * 64 + 3 * 16 + 0 * 4 + 1 * 1, 2),
            MinimizerIndex::minimizer_type(0 * 256 + 3 * 64 + 0 * 16 + 1 * 4 + 0 * 1, 3),
            MinimizerIndex::minimizer_type(0 * 256 + 1 * 64 + 0 * 16 + 0 * 4 + 3 * 1, 5),
            MinimizerIndex::minimizer_type(0 * 256 + 0 * 64 + 3 * 16 + 0 * 4 + 1 * 1, 7)
        };
        std::vector<MinimizerIndex::minimizer_type> result = index.minimizers(str.begin(), str.end());
        REQUIRE(result == correct);
    }

    SECTION("minimizers depend on window length") {
        MinimizerIndex index(5, 3);
        std::vector<MinimizerIndex::minimizer_type> correct {
            MinimizerIndex::minimizer_type(0 * 256 + 0 * 64 + 3 * 16 + 0 * 4 + 1 * 1, 2),
            MinimizerIndex::minimizer_type(0 * 256 + 1 * 64 + 0 * 16 + 0 * 4 + 3 * 1, 5),
            MinimizerIndex::minimizer_type(0 * 256 + 0 * 64 + 3 * 16 + 0 * 4 + 1 * 1, 7)
        };
        std::vector<MinimizerIndex::minimizer_type> result = index.minimizers(str.begin(), str.end());
        REQUIRE(result == correct);
    }

    SECTION("minimizers cannot contain invalid characters") {
        std::string weird = "CGAATAxAATACT";
        MinimizerIndex index(5, 2);
        std::vector<MinimizerIndex::minimizer_type> correct {
            MinimizerIndex::minimizer_type(1 * 256 + 2 * 64 + 0 * 16 + 0 * 4 + 3 * 1, 0),
            MinimizerIndex::minimizer_type(2 * 256 + 0 * 64 + 0 * 16 + 3 * 4 + 0 * 1, 1),
            MinimizerIndex::minimizer_type(0 * 256 + 0 * 64 + 3 * 16 + 0 * 4 + 1 * 1, 7)
        };
        std::vector<MinimizerIndex::minimizer_type> result = index.minimizers(weird.begin(), weird.end());
        REQUIRE(result == correct);
    }
}

void check_minimizer_index(const MinimizerIndex& index, const std::map<size_t, std::set<pos_t>>& correct_values,
                           size_t keys, size_t values, size_t unique, size_t frequent) {
    REQUIRE(index.size() == keys);
    REQUIRE(index.values() == values);
    REQUIRE(index.unique_keys() == unique);
    REQUIRE(index.frequent_keys() == frequent);

    for (auto iter = correct_values.begin(); iter != correct_values.end(); ++iter) {
        std::vector<pos_t> result = index.find(iter->first);
        std::vector<pos_t> correct(iter->second.begin(), iter->second.end());
        bool correct_result_or_empty_value = (result == correct || (correct.empty() && result.size() == 1 && is_empty(result.front())));
        REQUIRE(correct_result_or_empty_value);
    }
}

TEST_CASE("Index contains the correct kmers", "[minimizer_index][indexing]") {
    constexpr size_t TOTAL_KEYS = 16;

    SECTION("unique keys") {
        MinimizerIndex index;
        size_t keys = 0, values = 0, unique = 0, frequent = 0;
        std::map<size_t, std::set<pos_t>> correct_values;

        for (size_t i = 1; i <= TOTAL_KEYS; i++) {
            pos_t pos = MinimizerIndex::make_pos_t(i, i & 1, i & MinimizerIndex::OFF_MASK);
            index.insert(i, pos);
            correct_values[i].insert(pos);
            keys++; values++; unique++;
        }
        check_minimizer_index(index, correct_values, keys, values, unique, frequent);
    }

    SECTION("missing keys") {
        MinimizerIndex index;
        for (size_t i = 1; i <= TOTAL_KEYS; i++) {
            index.insert(i, MinimizerIndex::make_pos_t(i, i & 1, i & MinimizerIndex::OFF_MASK));
        }
        for (size_t i = TOTAL_KEYS + 1; i <= 2 * TOTAL_KEYS; i++) {
            REQUIRE(index.find(i).empty());
        }
    }

    SECTION("empty keys and values") {
        MinimizerIndex index;

        index.insert(MinimizerIndex::NO_KEY, MinimizerIndex::make_pos_t(1, false, 0));
        REQUIRE(index.find(MinimizerIndex::NO_KEY).empty());

        index.insert(TOTAL_KEYS + 1, MinimizerIndex::decode(MinimizerIndex::NO_VALUE));
        REQUIRE(index.find(TOTAL_KEYS + 1).empty());
    }

    SECTION("multiple occurrences") {
        MinimizerIndex index;
        size_t keys = 0, values = 0, unique = 0, frequent = 0;
        std::map<size_t, std::set<pos_t>> correct_values;

        for (size_t i = 1; i <= TOTAL_KEYS; i++) {
            pos_t pos = MinimizerIndex::make_pos_t(i, i & 1, i & MinimizerIndex::OFF_MASK);
            index.insert(i, pos);
            correct_values[i].insert(pos);
            keys++; values++; unique++;
        }
        for (size_t i = 1; i <= TOTAL_KEYS; i += 2) {
            pos_t pos = MinimizerIndex::make_pos_t(i + 1, i & 1, (i + 1) & MinimizerIndex::OFF_MASK);
            index.insert(i, pos);
            correct_values[i].insert(pos);
            values++; unique--;
        }
        for (size_t i = 1; i <= TOTAL_KEYS; i += 4) {
            pos_t pos = MinimizerIndex::make_pos_t(i + 2, i & 1, (i + 2) & MinimizerIndex::OFF_MASK);
            index.insert(i, pos);
            correct_values[i].insert(pos);
            values++;
        }
        check_minimizer_index(index, correct_values, keys, values, unique, frequent);
    }

    SECTION("duplicate values") {
        MinimizerIndex index;
        size_t keys = 0, values = 0, unique = 0, frequent = 0;
        std::map<size_t, std::set<pos_t>> correct_values;

        for (size_t i = 1; i <= TOTAL_KEYS; i++) {
            pos_t pos = MinimizerIndex::make_pos_t(i, i & 1, i & MinimizerIndex::OFF_MASK);
            index.insert(i, pos);
            correct_values[i].insert(pos);
            keys++; values++; unique++;
        }
        for (size_t i = 1; i <= TOTAL_KEYS; i += 2) {
            pos_t pos = MinimizerIndex::make_pos_t(i + 1, i & 1, (i + 1) & MinimizerIndex::OFF_MASK);
            index.insert(i, pos);
            correct_values[i].insert(pos);
            values++; unique--;
        }
        for (size_t i = 1; i <= TOTAL_KEYS; i += 4) {
            pos_t pos = MinimizerIndex::make_pos_t(i + 1, i & 1, (i + 1) & MinimizerIndex::OFF_MASK);
            index.insert(i, pos);
        }
        check_minimizer_index(index, correct_values, keys, values, unique, frequent);
    }

    SECTION("rehashing") {
        MinimizerIndex index;
        size_t keys = 0, values = 0, unique = 0, frequent = 0;
        std::map<size_t, std::set<pos_t>> correct_values;
        size_t threshold = index.max_keys();

        for (size_t i = 1; i <= threshold; i++) {
            pos_t pos = MinimizerIndex::make_pos_t(i, i & 1, i & MinimizerIndex::OFF_MASK);
            index.insert(i, pos);
            correct_values[i].insert(pos);
            keys++; values++; unique++;
        }
        REQUIRE(index.max_keys() == threshold);

        {
            size_t i = threshold + 1;
            pos_t pos = MinimizerIndex::make_pos_t(i, i & 1, i & MinimizerIndex::OFF_MASK);
            index.insert(i, pos);
            correct_values[i].insert(pos);
            keys++; values++; unique++;
        }
        REQUIRE(index.max_keys() > threshold);

        check_minimizer_index(index, correct_values, keys, values, unique, frequent);
    }
}

void check_empty_values(const MinimizerIndex& index, size_t keys, size_t step) {
    for (size_t i = 1; i < keys; i += step) {
        std::vector<pos_t> result = index.find(i);
        bool empty_value = (result.size() == 1 && is_empty(result.front()));
        REQUIRE(empty_value);
    }
}

TEST_CASE("Setting max_values works", "[minimizer_index][indexing]") {
    constexpr size_t TOTAL_KEYS = 16;

    SECTION("max_values = 1") {
        MinimizerIndex index(MinimizerIndex::KMER_LENGTH, MinimizerIndex::WINDOW_LENGTH, 1);
        size_t keys = 0, values = 0, unique = 0, frequent = 0;
        std::map<size_t, std::set<pos_t>> correct_values;

        for (size_t i = 1; i <= TOTAL_KEYS; i++) {
            pos_t pos = MinimizerIndex::make_pos_t(i, i & 1, i & MinimizerIndex::OFF_MASK);
            index.insert(i, pos);
            correct_values[i].insert(pos);
            keys++; values++; unique++;
        }
        for (size_t i = 1; i <= TOTAL_KEYS; i += 2) {
            pos_t pos = MinimizerIndex::make_pos_t(i + 1, i & 1, (i + 1) & MinimizerIndex::OFF_MASK);
            index.insert(i, pos);
            correct_values[i].clear();
            values--; unique--; frequent++;
        }
        for (size_t i = 1; i <= TOTAL_KEYS; i += 4) {
            pos_t pos = MinimizerIndex::make_pos_t(i + 2, i & 1, (i + 2) & MinimizerIndex::OFF_MASK);
            index.insert(i, pos);
        }
        for (size_t i = 1; i <= TOTAL_KEYS; i += 8) {
            pos_t pos = MinimizerIndex::make_pos_t(i + 3, i & 1, (i + 3) & MinimizerIndex::OFF_MASK);
            index.insert(i, pos);
        }
        check_minimizer_index(index, correct_values, keys, values, unique, frequent);
        check_empty_values(index, TOTAL_KEYS, 2);
    }

    SECTION("max_values = 2") {
        MinimizerIndex index(MinimizerIndex::KMER_LENGTH, MinimizerIndex::WINDOW_LENGTH, 2);
        size_t keys = 0, values = 0, unique = 0, frequent = 0;
        std::map<size_t, std::set<pos_t>> correct_values;

        for (size_t i = 1; i <= TOTAL_KEYS; i++) {
            pos_t pos = MinimizerIndex::make_pos_t(i, i & 1, i & MinimizerIndex::OFF_MASK);
            index.insert(i, pos);
            correct_values[i].insert(pos);
            keys++; values++; unique++;
        }
        for (size_t i = 1; i <= TOTAL_KEYS; i += 2) {
            pos_t pos = MinimizerIndex::make_pos_t(i + 1, i & 1, (i + 1) & MinimizerIndex::OFF_MASK);
            index.insert(i, pos);
            correct_values[i].insert(pos);
            values++; unique--;
        }
        for (size_t i = 1; i <= TOTAL_KEYS; i += 4) {
            pos_t pos = MinimizerIndex::make_pos_t(i + 2, i & 1, (i + 2) & MinimizerIndex::OFF_MASK);
            index.insert(i, pos);
            correct_values[i].clear();
            values -= 2; frequent++;
        }
        for (size_t i = 1; i <= TOTAL_KEYS; i += 8) {
            pos_t pos = MinimizerIndex::make_pos_t(i + 3, i & 1, (i + 3) & MinimizerIndex::OFF_MASK);
            index.insert(i, pos);
        }
        check_minimizer_index(index, correct_values, keys, values, unique, frequent);
        check_empty_values(index, TOTAL_KEYS, 4);
    }

    SECTION("max_values = 3") {
        MinimizerIndex index(MinimizerIndex::KMER_LENGTH, MinimizerIndex::WINDOW_LENGTH, 3);
        size_t keys = 0, values = 0, unique = 0, frequent = 0;
        std::map<size_t, std::set<pos_t>> correct_values;

        for (size_t i = 1; i <= TOTAL_KEYS; i++) {
            pos_t pos = MinimizerIndex::make_pos_t(i, i & 1, i & MinimizerIndex::OFF_MASK);
            index.insert(i, pos);
            correct_values[i].insert(pos);
            keys++; values++; unique++;
        }
        for (size_t i = 1; i <= TOTAL_KEYS; i += 2) {
            pos_t pos = MinimizerIndex::make_pos_t(i + 1, i & 1, (i + 1) & MinimizerIndex::OFF_MASK);
            index.insert(i, pos);
            correct_values[i].insert(pos);
            values++; unique--;
        }
        for (size_t i = 1; i <= TOTAL_KEYS; i += 4) {
            pos_t pos = MinimizerIndex::make_pos_t(i + 2, i & 1, (i + 2) & MinimizerIndex::OFF_MASK);
            index.insert(i, pos);
            correct_values[i].insert(pos);
            values++;
        }
        for (size_t i = 1; i <= TOTAL_KEYS; i += 8) {
            pos_t pos = MinimizerIndex::make_pos_t(i + 3, i & 1, (i + 3) & MinimizerIndex::OFF_MASK);
            index.insert(i, pos);
            correct_values[i].clear();
            values -= 3; frequent++;
        }
        check_minimizer_index(index, correct_values, keys, values, unique, frequent);
        check_empty_values(index, TOTAL_KEYS, 8);
    }
}

}
}
