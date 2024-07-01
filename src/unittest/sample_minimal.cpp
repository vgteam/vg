/// \file sample_minimal.cpp
///  
/// unit tests for minimizer (sub)sampling

#include "../algorithms/sample_minimal.hpp"
#include "catch.hpp"

#include <iostream>
#include <unordered_set>

namespace vg {
namespace unittest {

TEST_CASE("minimizer subsampling samples all tied minimizers", "[giraffe]") {
    // Say we have an element on every base of a sequence
    size_t sequence_length = 100;
    size_t element_length = 10;
    size_t element_count = sequence_length - element_length + 1;
    // This should work for any window size
    size_t window_size = 20;
    size_t window_count = sequence_length - window_size + 1;

    std::unordered_set<size_t> sampled_elements;

    algorithms::sample_minimal(element_count, element_length, window_size, sequence_length, [&](size_t i) {
        // Element i starts at offset i
        return i;
    }, [&](size_t a, size_t b) -> bool {
        // No element beats any other
        return false;
    }, [&](size_t sampled) {
        // Remember everything we sample
        sampled_elements.insert(sampled);
    });

    // If everything is tied, we should sample one element per window.
    REQUIRE(sampled_elements.size() == window_count);
}

TEST_CASE("minimizer subsampling samples both outer minimizers even if the first one is better", "[giraffe][subsampling]") {
    // Say we have an element on every base of a sequence
    size_t sequence_length = 100;
    size_t element_length = 10;
    std::vector<size_t> element_starts { 50, 55 };
    size_t element_count = element_starts.size();
    // Window should cover the whole clump of elements under test.
    size_t window_size = 20;

    std::unordered_set<size_t> sampled_elements;

    algorithms::sample_minimal(element_count, element_length, window_size, sequence_length, [&](size_t i) {
        return element_starts.at(i);
    }, [&](size_t a, size_t b) -> bool {
        // The first element beats all others
        return a == 0 && b != 0;
    }, [&](size_t sampled) {
        // Remember everything we sample
        sampled_elements.insert(sampled);
    });

    // We should sample both elements
    REQUIRE(sampled_elements.size() == 2);
    REQUIRE(sampled_elements.count(0));
    REQUIRE(sampled_elements.count(1));
}

TEST_CASE("minimizer subsampling samples both outer minimizers even if the second one is better", "[giraffe][subsampling]") {
    // Say we have an element on every base of a sequence
    size_t sequence_length = 100;
    size_t element_length = 10;
    std::vector<size_t> element_starts { 50, 55 };
    size_t element_count = element_starts.size();
    // Window should cover the whole clump of elements under test.
    size_t window_size = 20;

    std::unordered_set<size_t> sampled_elements;

    algorithms::sample_minimal(element_count, element_length, window_size, sequence_length, [&](size_t i) {
        return element_starts.at(i);
    }, [&](size_t a, size_t b) -> bool {
        // The second element beats all others
        return a == 1 && b != 1;
    }, [&](size_t sampled) {
        // Remember everything we sample
        sampled_elements.insert(sampled);
    });

    // We should sample both elements
    REQUIRE(sampled_elements.size() == 2);
    REQUIRE(sampled_elements.count(0));
    REQUIRE(sampled_elements.count(1));
}

TEST_CASE("minimizer subsampling samples only outer elements if a middle one is worst", "[giraffe][subsampling]") {
    // Say we have an element on every base of a sequence
    size_t sequence_length = 100;
    size_t element_length = 10;
    std::vector<size_t> element_starts { 50, 55, 58 };
    std::vector<size_t> element_goodness { 10, 0, 11 };
    size_t element_count = element_starts.size();
    // Window should cover the whole clump of elements under test.
    size_t window_size = 20;

    std::unordered_set<size_t> sampled_elements;

    algorithms::sample_minimal(element_count, element_length, window_size, sequence_length, [&](size_t i) {
        return element_starts.at(i);
    }, [&](size_t a, size_t b) -> bool {
        return element_goodness.at(a) > element_goodness.at(b);
    }, [&](size_t sampled) {
        // Remember everything we sample
        sampled_elements.insert(sampled);
    });

    // We should sample the outer elements
    REQUIRE(sampled_elements.size() == 2);
    REQUIRE(sampled_elements.count(0));
    REQUIRE(sampled_elements.count(2));
}

TEST_CASE("minimizer subsampling samples all 3 elements if the middle one is better than the first", "[giraffe][subsampling]") {
    // Say we have an element on every base of a sequence
    size_t sequence_length = 100;
    size_t element_length = 10;
    std::vector<size_t> element_starts { 50, 55, 58 };
    std::vector<size_t> element_goodness { 5, 10, 11 };
    size_t element_count = element_starts.size();
    // Window should cover the whole clump of elements under test.
    size_t window_size = 20;

    std::unordered_set<size_t> sampled_elements;

    algorithms::sample_minimal(element_count, element_length, window_size, sequence_length, [&](size_t i) {
        return element_starts.at(i);
    }, [&](size_t a, size_t b) {
        return element_goodness.at(a) > element_goodness.at(b);
    }, [&](size_t sampled) {
        // Remember everything we sample
        sampled_elements.insert(sampled);
    });

    // We should sample all the elements
    REQUIRE(sampled_elements.size() == 3);
    REQUIRE(sampled_elements.count(0));
    REQUIRE(sampled_elements.count(1));
    REQUIRE(sampled_elements.count(2));
}

TEST_CASE("minimizer subsampling samples all 3 elements if the middle one is better than the last", "[giraffe][subsampling]") {
    // Say we have an element on every base of a sequence
    size_t sequence_length = 100;
    size_t element_length = 10;
    std::vector<size_t> element_starts { 50, 55, 58 };
    std::vector<size_t> element_goodness { 11, 10, 5 };
    size_t element_count = element_starts.size();
    // Window should cover the whole clump of elements under test.
    size_t window_size = 20;

    std::unordered_set<size_t> sampled_elements;

    algorithms::sample_minimal(element_count, element_length, window_size, sequence_length, [&](size_t i) {
        return element_starts.at(i);
    }, [&](size_t a, size_t b) {
        return element_goodness.at(a) > element_goodness.at(b);
    }, [&](size_t sampled) {
        // Remember everything we sample
        sampled_elements.insert(sampled);
    });

    // We should sample all the elements
    REQUIRE(sampled_elements.size() == 3);
    REQUIRE(sampled_elements.count(0));
    REQUIRE(sampled_elements.count(1));
    REQUIRE(sampled_elements.count(2));
}


}
}
