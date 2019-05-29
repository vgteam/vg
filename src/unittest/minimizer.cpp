/** \file
 *
 * Unit tests for MinimizerIndex, which is basically a hash table mapping kmers
 * to sets of pos_t.
 */

#include "../minimizer.hpp"
#include "../json2pb.h"
#include "../position.hpp"
#include "../utility.hpp"
#include "../wang_hash.hpp"

#include "catch.hpp"

#include <map>
#include <set>
#include <vector>

namespace vg {
namespace unittest {

namespace {

MinimizerIndex::minimizer_type get_minimizer(MinimizerIndex::key_type key, MinimizerIndex::offset_type offset = 0, bool orientation = false) {
    return { key, wang_hash_64(key), offset, orientation };
}

}

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
        default_index.insert(get_minimizer(1), make_pos_t(1, false, 3));
        REQUIRE(default_index != default_copy);

        // Same key, different value.
        default_copy.insert(get_minimizer(1), make_pos_t(2, false, 3));
        REQUIRE(default_index != default_copy);

        // Same contents.
        default_copy = default_index;
        REQUIRE(default_index == default_copy);
    }

    SECTION("index contents can be swapped") {
        MinimizerIndex first, second;
        first.insert(get_minimizer(1), make_pos_t(1, false, 3));
        second.insert(get_minimizer(2), make_pos_t(2, false, 3));

        MinimizerIndex first_copy(first), second_copy(second);
        first.swap(second);
        REQUIRE(first != first_copy);
        REQUIRE(first == second_copy);
        REQUIRE(second == first_copy);
        REQUIRE(second != second_copy);
    }

    SECTION("serialization preserves parameters and contents") {
        MinimizerIndex index(15, 6);
        index.insert(get_minimizer(1), make_pos_t(1, false, 3));
        index.insert(get_minimizer(2), make_pos_t(1, false, 3));
        index.insert(get_minimizer(2), make_pos_t(2, false, 3));

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

// FIXME orientation; same minimizers in both orientations
// wang_hash_64() order of 3-mers:
// AAT < TGT < TTG < TAT < ATA < TCG < ATT < ACA < GAA < ACT < TAC < CGA < CAA < GTA < TTC < AGT
TEST_CASE("Minimizer extraction works correctly", "[minimizer_index][indexing]") {
    std::string str = "CGAATACAATACT";
    std::string rev = reverse_complement(str);

    // Find the correct order of 3-mers.
    /*for (size_t i = 0; i < 64; i++) {
        std::string kmer;
        kmer += MinimizerIndex::PACK_TO_CHAR[(i / 16) & 3];
        kmer += MinimizerIndex::PACK_TO_CHAR[(i / 4) & 3];
        kmer += MinimizerIndex::PACK_TO_CHAR[(i) & 3];
        size_t hash = wang_hash_64(i);
        if (str.find(kmer) != std::string::npos || rev.find(kmer) != std::string::npos) {
            std::cerr << hash << " " << kmer << std::endl;
        }
    }*/

    SECTION("minimizer is the leftmost occurrence") {
        MinimizerIndex index(3, 2);
        MinimizerIndex::minimizer_type correct = get_minimizer(0 * 16 + 0 * 4 + 3 * 1, 2); // AAT
        MinimizerIndex::minimizer_type result = index.minimizer(str.begin(), str.end());
        REQUIRE(result == correct);
    }

    SECTION("all minimizers are found") {
        MinimizerIndex index(3, 2);
        std::vector<MinimizerIndex::minimizer_type> correct {
            get_minimizer(3 * 16 + 1 * 4 + 2 * 1, 2, true),   // TCG
            get_minimizer(0 * 16 + 0 * 4 + 3 * 1, 2, false),  // AAT
            get_minimizer(3 * 16 + 0 * 4 + 3 * 1, 5, true),   // TAT
            get_minimizer(3 * 16 + 2 * 4 + 3 * 1, 7, true),   // TGT
            get_minimizer(0 * 16 + 0 * 4 + 3 * 1, 7, false),  // AAT
            get_minimizer(3 * 16 + 0 * 4 + 3 * 1, 10, true),  // TAT
            get_minimizer(0 * 16 + 1 * 4 + 3 * 1, 10, false)  // ACT
        };
        std::vector<MinimizerIndex::minimizer_type> result = index.minimizers(str.begin(), str.end());
        REQUIRE(result == correct);
    }

    SECTION("minimizers depend on window length") {
        MinimizerIndex index(3, 3);
        std::vector<MinimizerIndex::minimizer_type> correct {
            get_minimizer(0 * 16 + 0 * 4 + 3 * 1, 2, false),  // AAT
            get_minimizer(3 * 16 + 2 * 4 + 3 * 1, 7, true),   // TGT
            get_minimizer(0 * 16 + 0 * 4 + 3 * 1, 7, false),  // AAT
            get_minimizer(3 * 16 + 0 * 4 + 3 * 1, 10, true)   // TAT
        };
        std::vector<MinimizerIndex::minimizer_type> result = index.minimizers(str.begin(), str.end());
        REQUIRE(result == correct);
    }

    SECTION("minimizers cannot contain invalid characters") {
        std::string weird = "CGAATAxAATACT";
        MinimizerIndex index(3, 2);
        std::vector<MinimizerIndex::minimizer_type> correct {
            get_minimizer(3 * 16 + 1 * 4 + 2 * 1, 2, true),   // TCG
            get_minimizer(0 * 16 + 0 * 4 + 3 * 1, 2, false),  // AAT
            get_minimizer(3 * 16 + 0 * 4 + 3 * 1, 5, true),   // TAT
            get_minimizer(0 * 16 + 0 * 4 + 3 * 1, 7, false),  // AAT
            get_minimizer(3 * 16 + 0 * 4 + 3 * 1, 10, true),  // TAT
            get_minimizer(0 * 16 + 1 * 4 + 3 * 1, 10, false)  // ACT
        };
        std::vector<MinimizerIndex::minimizer_type> result = index.minimizers(weird.begin(), weird.end());
        REQUIRE(result == correct);
    }

    SECTION("both orientations have the same minimizers") {
        MinimizerIndex index(3, 2);
        std::vector<MinimizerIndex::minimizer_type> forward_minimizers = index.minimizers(str.begin(), str.end());
        std::vector<MinimizerIndex::minimizer_type> reverse_minimizers = index.minimizers(rev.begin(), rev.end());
        REQUIRE(forward_minimizers.size() == reverse_minimizers.size());
        for (size_t i = 0; i < forward_minimizers.size(); i++) {
            MinimizerIndex::minimizer_type& f = forward_minimizers[i];
            MinimizerIndex::minimizer_type& r = reverse_minimizers[forward_minimizers.size() - 1 - i];
            REQUIRE(f.key == r.key);
            REQUIRE(f.offset == str.length() - 1 - r.offset);
            REQUIRE(f.is_reverse != r.is_reverse);
        }
    }
}

void check_minimizer_index(const MinimizerIndex& index, const std::map<size_t, std::set<pos_t>>& correct_values,
                           size_t keys, size_t values, size_t unique) {
    REQUIRE(index.size() == keys);
    REQUIRE(index.values() == values);
    REQUIRE(index.unique_keys() == unique);

    for (auto iter = correct_values.begin(); iter != correct_values.end(); ++iter) {
        std::vector<pos_t> result = index.find(get_minimizer(iter->first));
        std::vector<pos_t> correct(iter->second.begin(), iter->second.end());
        bool correct_result_or_empty_value = (result == correct || (correct.empty() && result.size() == 1 && is_empty(result.front())));
        REQUIRE(correct_result_or_empty_value);
    }
}

TEST_CASE("Index contains the correct kmers", "[minimizer_index][indexing]") {
    constexpr size_t TOTAL_KEYS = 16;

    SECTION("unique keys") {
        MinimizerIndex index;
        size_t keys = 0, values = 0, unique = 0;
        std::map<size_t, std::set<pos_t>> correct_values;

        for (size_t i = 1; i <= TOTAL_KEYS; i++) {
            pos_t pos = make_pos_t(i, i & 1, i & MinimizerIndex::OFF_MASK);
            index.insert(get_minimizer(i), pos);
            correct_values[i].insert(pos);
            keys++; values++; unique++;
        }
        check_minimizer_index(index, correct_values, keys, values, unique);
    }

    SECTION("missing keys") {
        MinimizerIndex index;
        for (size_t i = 1; i <= TOTAL_KEYS; i++) {
            index.insert(get_minimizer(i), make_pos_t(i, i & 1, i & MinimizerIndex::OFF_MASK));
        }
        for (size_t i = TOTAL_KEYS + 1; i <= 2 * TOTAL_KEYS; i++) {
            REQUIRE(index.find(get_minimizer(i)).empty());
        }
    }

    SECTION("empty keys and values") {
        MinimizerIndex index;

        index.insert(get_minimizer(MinimizerIndex::NO_KEY), make_pos_t(1, false, 0));
        REQUIRE(index.find(get_minimizer(MinimizerIndex::NO_KEY)).empty());

        index.insert(get_minimizer(TOTAL_KEYS + 1), MinimizerIndex::decode(MinimizerIndex::NO_VALUE));
        REQUIRE(index.find(get_minimizer(TOTAL_KEYS + 1)).empty());
    }

    SECTION("multiple occurrences") {
        MinimizerIndex index;
        size_t keys = 0, values = 0, unique = 0;
        std::map<size_t, std::set<pos_t>> correct_values;

        for (size_t i = 1; i <= TOTAL_KEYS; i++) {
            pos_t pos = make_pos_t(i, i & 1, i & MinimizerIndex::OFF_MASK);
            index.insert(get_minimizer(i), pos);
            correct_values[i].insert(pos);
            keys++; values++; unique++;
        }
        for (size_t i = 1; i <= TOTAL_KEYS; i += 2) {
            pos_t pos = make_pos_t(i + 1, i & 1, (i + 1) & MinimizerIndex::OFF_MASK);
            index.insert(get_minimizer(i), pos);
            correct_values[i].insert(pos);
            values++; unique--;
        }
        for (size_t i = 1; i <= TOTAL_KEYS; i += 4) {
            pos_t pos = make_pos_t(i + 2, i & 1, (i + 2) & MinimizerIndex::OFF_MASK);
            index.insert(get_minimizer(i), pos);
            correct_values[i].insert(pos);
            values++;
        }
        check_minimizer_index(index, correct_values, keys, values, unique);
    }

    SECTION("duplicate values") {
        MinimizerIndex index;
        size_t keys = 0, values = 0, unique = 0;
        std::map<size_t, std::set<pos_t>> correct_values;

        for (size_t i = 1; i <= TOTAL_KEYS; i++) {
            pos_t pos = make_pos_t(i, i & 1, i & MinimizerIndex::OFF_MASK);
            index.insert(get_minimizer(i), pos);
            correct_values[i].insert(pos);
            keys++; values++; unique++;
        }
        for (size_t i = 1; i <= TOTAL_KEYS; i += 2) {
            pos_t pos = make_pos_t(i + 1, i & 1, (i + 1) & MinimizerIndex::OFF_MASK);
            index.insert(get_minimizer(i), pos);
            correct_values[i].insert(pos);
            values++; unique--;
        }
        for (size_t i = 1; i <= TOTAL_KEYS; i += 4) {
            pos_t pos = make_pos_t(i + 1, i & 1, (i + 1) & MinimizerIndex::OFF_MASK);
            index.insert(get_minimizer(i), pos);
        }
        check_minimizer_index(index, correct_values, keys, values, unique);
    }

    SECTION("rehashing") {
        MinimizerIndex index;
        size_t keys = 0, values = 0, unique = 0;
        std::map<size_t, std::set<pos_t>> correct_values;
        size_t threshold = index.max_keys();

        for (size_t i = 1; i <= threshold; i++) {
            pos_t pos = make_pos_t(i, i & 1, i & MinimizerIndex::OFF_MASK);
            index.insert(get_minimizer(i), pos);
            correct_values[i].insert(pos);
            keys++; values++; unique++;
        }
        REQUIRE(index.max_keys() == threshold);

        {
            size_t i = threshold + 1;
            pos_t pos = make_pos_t(i, i & 1, i & MinimizerIndex::OFF_MASK);
            index.insert(get_minimizer(i), pos);
            correct_values[i].insert(pos);
            keys++; values++; unique++;
        }
        REQUIRE(index.max_keys() > threshold);

        check_minimizer_index(index, correct_values, keys, values, unique);
    }
}

}
}
