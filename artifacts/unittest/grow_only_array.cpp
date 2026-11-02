/// \file grow_only_array.cpp
///
/// Unit tests for GrowOnlyArray.

#include <atomic>
#include <string>
#include <vector>

#include <omp.h>

#include "../../src/unittest/catch.hpp"
#include "../grow_only_array.hpp"

namespace vg {
namespace unittest {

using namespace std;

namespace {

/// An element that counts how many of its kind are alive.
struct Counted {
    static inline atomic<int> alive {0};
    int value;

    explicit Counted(int value) noexcept : value(value) {
        ++alive;
    }
    ~Counted() {
        --alive;
    }
    Counted(const Counted&) = delete;
    Counted& operator=(const Counted&) = delete;
};

}

TEST_CASE("GrowOnlyArray hands out indices in order and reads them back across blocks",
          "[grow_only_array]") {
    // Blocks of 2, 4, 8, ... elements, so a few hundred elements span many blocks.
    GrowOnlyArray<string, 1> array;
    for (size_t i = 0; i < 300; ++i) {
        REQUIRE(array.emplace_back(to_string(i)) == i);
    }
    for (size_t i = 0; i < 300; ++i) {
        REQUIRE(array[i] == to_string(i));
    }
}

TEST_CASE("GrowOnlyArray elements keep their addresses as the array grows", "[grow_only_array]") {
    GrowOnlyArray<int, 1> array;
    vector<const int*> addresses;
    for (int i = 0; i < 100; ++i) {
        array.emplace_back(i);
        addresses.push_back(&array[i]);
    }
    for (int i = 100; i < 1000; ++i) {
        array.emplace_back(i);
    }
    for (int i = 0; i < 100; ++i) {
        REQUIRE(&array[i] == addresses[i]);
        REQUIRE(*addresses[i] == i);
    }
}

TEST_CASE("GrowOnlyArray takes elements from several threads at once", "[grow_only_array]") {
    const int per_thread = 5000;
    const int threads = 8;
    GrowOnlyArray<int, 2> array;
    // Each thread records the index of each value it added, and reads back its own as it goes.
    vector<vector<size_t>> indices(threads);
    atomic<bool> all_read {true};
#pragma omp parallel for num_threads(threads)
    for (int t = 0; t < threads; ++t) {
        for (int i = 0; i < per_thread; ++i) {
            const int value = t * per_thread + i;
            const size_t index = array.emplace_back(value);
            indices[t].push_back(index);
            if (array[index] != value) {
                all_read = false;
            }
        }
    }
    REQUIRE(all_read);
    // Every index was handed out once, and holds the value added with it.
    vector<bool> seen(threads * per_thread, false);
    for (int t = 0; t < threads; ++t) {
        for (int i = 0; i < per_thread; ++i) {
            const size_t index = indices[t][i];
            REQUIRE(index < seen.size());
            REQUIRE_FALSE(seen[index]);
            seen[index] = true;
            REQUIRE(array[index] == t * per_thread + i);
        }
    }
}

TEST_CASE("GrowOnlyArray destroys its elements with it", "[grow_only_array]") {
    Counted::alive = 0;
    {
        GrowOnlyArray<Counted, 1> array;
        for (int i = 0; i < 50; ++i) {
            array.emplace_back(i);
        }
        REQUIRE(Counted::alive == 50);
        REQUIRE(array[49].value == 49);
    }
    REQUIRE(Counted::alive == 0);
}

}
}
