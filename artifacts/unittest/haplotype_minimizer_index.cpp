/** \file
 *
 * Unit tests for haplotype_minimizer_index.cpp: the index must hold the keys and positions that
 * gbwtgraph's minimizer index holds for the same graph, whatever the number of threads.
 */

#include "../haplotype_minimizer_index.hpp"
#include "../../src/gbwt_helper.hpp"

#include "catch.hpp"

#include <gbwtgraph/index.h>

#include <omp.h>

#include <random>
#include <string>
#include <vector>

namespace vg::unittest {

namespace {

using key_type = HaplotypeMinimizerIndex::key_type;
using gbwtgraph_index_type = gbwtgraph::MinimizerIndex<key_type>;
using node_t = gbwt::vector_type::value_type;

// The number of bubbles in the test graph. The graph has more than 1024 nodes, so that the
// parallel window traversal hands nodes to several threads.
constexpr size_t BUBBLES = 1500;

// The number of haplotypes in the test graph.
constexpr size_t HAPLOTYPES = 6;

node_t forward(nid_t id) {
    return static_cast<node_t>(gbwt::Node::encode(id, false));
}

std::string random_sequence(std::mt19937_64& rng, size_t length) {
    static const std::string bases = "ACGT";
    std::string result;
    for (size_t i = 0; i < length; i++) {
        result.push_back(bases[rng() % bases.size()]);
    }
    return result;
}

size_t random_length(std::mt19937_64& rng, size_t low, size_t high) {
    return low + rng() % (high - low + 1);
}

// Add the nodes of a chain of bubbles to the graph. Shared node i has ID 4i + 1, and the alleles
// of bubble i have IDs 4i + 2 and 4i + 3, followed by shared node i + 1.
//
// A sequence repeated in every 20th shared node makes some long minimizers occur at several
// positions. Shared node 300 has N bases.
//
// Return, for each bubble, whether haplotypes may skip its alleles.
std::vector<bool> add_bubble_chain(gbwtgraph::NaiveGraph& graph, std::mt19937_64& rng) {
    std::string repeat = random_sequence(rng, 60);
    std::vector<bool> skippable(BUBBLES, false);
    for (size_t i = 0; i <= BUBBLES; i++) {
        std::string shared = random_sequence(rng, random_length(rng, 20, 120));
        if (i % 20 == 7) {
            shared += repeat;
        }
        if (i == 300) {
            shared += "NNNNN" + random_sequence(rng, 10);
        }
        graph.create_node(4 * i + 1, shared);
        if (i < BUBBLES) {
            graph.create_node(4 * i + 2, random_sequence(rng, random_length(rng, 1, 20)));
            graph.create_node(4 * i + 3, random_sequence(rng, random_length(rng, 1, 20)));
            skippable[i] = (rng() % 4 == 0);
        }
    }
    return skippable;
}

// Return the walk of a haplotype through the chain of bubbles, with a random allele in each.
//
// Haplotypes 1 and 3 visit the chain in reverse. Haplotype 4 covers only shared nodes 200 to
// 800. Haplotype 5 visits shared nodes 100 to 110 twice.
gbwt::vector_type haplotype_walk(size_t haplotype, const std::vector<bool>& skippable, std::mt19937_64& rng) {
    size_t first = (haplotype == 4 ? 200 : 0);
    size_t last = (haplotype == 4 ? 800 : BUBBLES);
    gbwt::vector_type walk;
    for (size_t i = first; i < last; i++) {
        walk.push_back(forward(4 * i + 1));
        if (haplotype == 5 && i == 110) {
            for (size_t j = 100; j <= 110; j++) {
                walk.push_back(forward(4 * j + 1));
            }
        }
        size_t allele = rng() % (skippable[i] ? 3 : 2);
        if (allele < 2) {
            walk.push_back(forward(4 * i + 2 + allele));
        }
    }
    walk.push_back(forward(4 * last + 1));
    if (haplotype == 1 || haplotype == 3) {
        gbwt::reversePath(walk);
    }
    return walk;
}

// Return a GBZ graph of a chain of bubbles with haplotypes that choose random alleles.
gbwtgraph::GBZ bubble_chain_graph() {
    std::mt19937_64 rng(42);
    gbwtgraph::NaiveGraph graph;
    std::vector<bool> skippable = add_bubble_chain(graph, rng);

    std::vector<gbwt::vector_type> walks;
    for (size_t haplotype = 0; haplotype < HAPLOTYPES; haplotype++) {
        walks.push_back(haplotype_walk(haplotype, skippable, rng));
    }
    gbwt::GBWT index = get_gbwt(walks);
    index.addMetadata();
    index.metadata.setSamples({ "s0", "s1", "s2", "s3", "s4", "s5" });
    index.metadata.setHaplotypes(walks.size());
    index.metadata.setContigs({ "chr" });
    for (size_t i = 0; i < walks.size(); i++) {
        index.metadata.addPath(i, 0, 0, 0);
    }
    return gbwtgraph::GBZ(index, graph);
}

// Build the index on the given number of threads.
HaplotypeMinimizerIndex build_on_threads(const gbwtgraph::GBZ& gbz, size_t k, size_t w, int threads) {
    int old_threads = omp_get_max_threads();
    omp_set_num_threads(threads);
    HaplotypeMinimizerIndex result(gbz, k, w);
    omp_set_num_threads(old_threads);
    return result;
}

// Build gbwtgraph's minimizer index with no payload.
gbwtgraph_index_type build_gbwtgraph_index(const gbwtgraph::GBZ& gbz, size_t k, size_t w) {
    gbwtgraph_index_type result(k, w, 0, false);
    gbwtgraph::index_haplotypes(gbz, result, [](const gbwtgraph::pos_t&) -> const gbwtgraph::KmerEncoding::code_type* {
        return nullptr;
    });
    return result;
}

// Return the positions in a list from find() of gbwtgraph's index, decoded.
std::vector<gbwtgraph::Position> positions_in(const gbwtgraph_index_type& index, gbwtgraph::KmerEncoding::multi_value_type values) {
    std::vector<gbwtgraph::Position> result;
    for (size_t i = 0; i < values.second; i++) {
        result.push_back(index.get_value(values, i).first);
    }
    return result;
}

// Return the positions of the key in the haplotype minimizer index, decoded.
std::vector<gbwtgraph::Position> positions_in(const HaplotypeMinimizerIndex& index, key_type key) {
    std::vector<gbwtgraph::Position> result;
    for (HaplotypeMinimizerIndex::code_type position : index.positions(key.get_key())) {
        result.push_back(gbwtgraph::Position(position));
    }
    return result;
}

// Require the index to hold the keys and positions of gbwtgraph's index.
void require_same_occurrences(const HaplotypeMinimizerIndex& index, const gbwtgraph_index_type& expected) {
    REQUIRE(index.k() == expected.k());
    REQUIRE(index.w() == expected.w());
    REQUIRE(index.size() == expected.size());
    REQUIRE(index.number_of_values() == expected.number_of_values());
    REQUIRE(index.unique_keys() == expected.unique_keys());

    // The key totals are equal, so checking every key of gbwtgraph's index also rules out extra keys here.
    size_t mismatches = 0;
    expected.for_each_minimizer([&index, &expected, &mismatches](key_type key, gbwtgraph::KmerEncoding::multi_value_type values) {
        if (index.count(key.get_key()) != values.second || positions_in(index, key) != positions_in(expected, values)) {
            mismatches++;
        }
        return true;
    });
    REQUIRE(mismatches == 0);
}

} // anonymous namespace

TEST_CASE("The haplotype minimizer index holds the occurrences of gbwtgraph's index", "[haplotypes][haplotype_minimizer_index]") {
    gbwtgraph::GBZ gbz = bubble_chain_graph();

    for (auto [k, w] : std::vector<std::pair<size_t, size_t>> { { 29, 11 }, { 7, 3 } }) {
        gbwtgraph_index_type expected = build_gbwtgraph_index(gbz, k, w);
        // Some keys must have several positions, or the comparison would not test how the positions
        // found in different windows and threads are merged.
        REQUIRE(expected.unique_keys() < expected.size());

        for (int threads : { 1, 4 }) {
            DYNAMIC_SECTION("k " << k << ", w " << w << ", " << threads << " threads") {
                HaplotypeMinimizerIndex index = build_on_threads(gbz, k, w, threads);
                require_same_occurrences(index, expected);
            }
        }
    }
}

TEST_CASE("The haplotype minimizer index reports no positions for an absent key", "[haplotypes][haplotype_minimizer_index]") {
    gbwtgraph::GBZ gbz = bubble_chain_graph();
    HaplotypeMinimizerIndex index(gbz, 29, 11);
    gbwtgraph_index_type expected = build_gbwtgraph_index(gbz, 29, 11);

    // A run of a single base is absent from random sequence.
    auto absent = expected.minimizers(std::string(39, 'A')).front();
    REQUIRE(expected.count(absent) == 0);
    REQUIRE(index.count(absent.key.get_key()) == 0);
    REQUIRE(index.positions(absent.key.get_key()).empty());
}

}
