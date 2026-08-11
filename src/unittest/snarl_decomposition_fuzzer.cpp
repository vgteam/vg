#include "catch.hpp"
#include "../handle.hpp"
#include <bdsg/hash_graph.hpp>

#include "support/snarl_decomposition_fuzzer.hpp"

#include <vector>
#include <unordered_set>

namespace vg {
namespace unittest {

using ET = DecompositionEventType;
using Event = DecompositionEvent;

TEST_CASE("ReplaySnarlFinder replays events faithfully", "[snarl_decomposition_fuzzer]") {
    // Build a small graph to get real handles
    bdsg::HashGraph graph;
    graph.create_handle("A", 1);
    graph.create_handle("C", 2);
    graph.create_handle("G", 3);
    graph.create_handle("T", 4);
    graph.create_handle("AA", 5);

    std::vector<Event> events = {
        {ET::BEGIN_CHAIN, 1, false},
          {ET::BEGIN_SNARL, 1, false},
            {ET::BEGIN_CHAIN, 2, true},
            {ET::END_CHAIN, 3, true},
          {ET::END_SNARL, 4, false},
          {ET::BEGIN_SNARL, 4, false},
          {ET::END_SNARL, 5, false},
        {ET::END_CHAIN, 5, false},
    };

    ReplaySnarlFinder finder(&graph, events);
    std::vector<Event> captured = capture_events(finder, graph);

    REQUIRE(captured == events);
}

TEST_CASE("SnarlDecompositionFuzzer passes through when nothing is flipped", "[snarl_decomposition_fuzzer]") {
    bdsg::HashGraph graph;
    graph.create_handle("A", 1);
    graph.create_handle("C", 2);
    graph.create_handle("G", 3);
    graph.create_handle("T", 4);
    graph.create_handle("AA", 5);

    std::vector<Event> events = {
        {ET::BEGIN_CHAIN, 1, false},
          {ET::BEGIN_SNARL, 1, false},
            {ET::BEGIN_CHAIN, 2, true},
            {ET::END_CHAIN, 3, true},
          {ET::END_SNARL, 4, false},
          {ET::BEGIN_SNARL, 4, false},
          {ET::END_SNARL, 5, false},
        {ET::END_CHAIN, 5, false},
    };

    ReplaySnarlFinder replay(&graph, events);

    // No chains to flip
    SnarlDecompositionFuzzer fuzzer(&graph, &replay, {});

    std::vector<Event> captured = capture_events(fuzzer, graph);

    REQUIRE(captured == events);
}

TEST_CASE("SnarlDecompositionFuzzer flips an outer chain", "[snarl_decomposition_fuzzer]") {
    // Graph:
    // Chain: 1fwd -> snarl(1fwd, 4fwd) -> snarl(4fwd, 5fwd) -> 5fwd
    // Inside first snarl: chain 2rev->3rev
    bdsg::HashGraph graph;
    graph.create_handle("A", 1);
    graph.create_handle("C", 2);
    graph.create_handle("G", 3);
    graph.create_handle("T", 4);
    graph.create_handle("AA", 5);
    
    std::vector<Event> events = {
        {ET::BEGIN_CHAIN, 1, false},
          {ET::BEGIN_SNARL, 1, false},
            {ET::BEGIN_CHAIN, 2, true},
            {ET::END_CHAIN, 3, true},
          {ET::END_SNARL, 4, false},
          {ET::BEGIN_SNARL, 4, false},
          {ET::END_SNARL, 5, false},
        {ET::END_CHAIN, 5, false},
    };

    ReplaySnarlFinder replay(&graph, events);

    SECTION("flip outer chain only") {
        // Flip the outer chain (1fwd -> 5fwd)
        std::unordered_set<nid_t> flips {1, 5};
        SnarlDecompositionFuzzer fuzzer(&graph, &replay, flips);

        std::vector<Event> captured = capture_events(fuzzer, graph);

        // Expected after flipping the outer chain:
        // Flipping a chain reverses everything inside it, including children.
        // The nested chain 2rev->3rev gets reversed to 3fwd->2fwd as
        // part of the parent flip.
        std::vector<Event> expected = {
            {ET::BEGIN_CHAIN, 5, true},
              {ET::BEGIN_SNARL, 5, true},
              {ET::END_SNARL, 4, true},
              {ET::BEGIN_SNARL, 4, true},
                {ET::BEGIN_CHAIN, 3, false},
                {ET::END_CHAIN, 2, false},
              {ET::END_SNARL, 1, true},
            {ET::END_CHAIN, 1, true},
        };

        REQUIRE(captured == expected);
    }

    SECTION("flip outer and nested chain") {
        // Flip outer chain (1fwd->5fwd) AND nested chain (2rev->3rev)
        std::unordered_set<nid_t> flips {1, 5, 2, 3};
        SnarlDecompositionFuzzer fuzzer(&graph, &replay, flips);

        std::vector<Event> captured = capture_events(fuzzer, graph);

        // Expected: outer chain flipped (reversing everything, including
        // the nested chain to 3fwd->2fwd), AND THEN the nested chain is
        // flipped again back to its original orientation 2rev->3rev.
        std::vector<Event> expected = {
            {ET::BEGIN_CHAIN, 5, true},
              {ET::BEGIN_SNARL, 5, true},
              {ET::END_SNARL, 4, true},
              {ET::BEGIN_SNARL, 4, true},
                {ET::BEGIN_CHAIN, 2, true},
                {ET::END_CHAIN, 3, true},
              {ET::END_SNARL, 1, true},
            {ET::END_CHAIN, 1, true},
        };

        REQUIRE(captured == expected);
    }

    SECTION("flip nested chain only") {
        // Flip only the nested chain (2rev->3rev), outer stays
        std::unordered_set<nid_t> flips {2, 3};
        SnarlDecompositionFuzzer fuzzer(&graph, &replay, flips);

        std::vector<Event> captured = capture_events(fuzzer, graph);

        // Outer chain not flipped, nested chain flipped
        std::vector<Event> expected = {
            {ET::BEGIN_CHAIN, 1, false},
              {ET::BEGIN_SNARL, 1, false},
                {ET::BEGIN_CHAIN, 3, false},
                {ET::END_CHAIN, 2, false},
              {ET::END_SNARL, 4, false},
              {ET::BEGIN_SNARL, 4, false},
              {ET::END_SNARL, 5, false},
            {ET::END_CHAIN, 5, false},
        };

        REQUIRE(captured == expected);
    }
}

TEST_CASE("SnarlDecompositionFuzzer handles empty chain", "[snarl_decomposition_fuzzer]") {
    bdsg::HashGraph graph;
    graph.create_handle("ACGT", 1);

    // An empty chain: begin and end with same handle, no snarls inside
    std::vector<Event> events = {
        {ET::BEGIN_CHAIN, 1, false},
        {ET::END_CHAIN, 1, false},
    };

    ReplaySnarlFinder replay(&graph, events);

    SECTION("flipping an empty chain") {
        std::unordered_set<nid_t> flips {1};
        SnarlDecompositionFuzzer fuzzer(&graph, &replay, flips);

        std::vector<Event> captured = capture_events(fuzzer, graph);

        std::vector<Event> expected = {
            {ET::BEGIN_CHAIN, 1, true},
            {ET::END_CHAIN, 1, true},
        };

        REQUIRE(captured == expected);
    }
}

TEST_CASE("SnarlDecompositionFuzzer handles multiple top-level chains", "[snarl_decomposition_fuzzer]") {
    bdsg::HashGraph graph;
    graph.create_handle("A", 1);
    graph.create_handle("C", 2);
    graph.create_handle("G", 3);
    graph.create_handle("T", 4);

    // Two top-level chains in the root snarl
    std::vector<Event> events = {
        // Chain 1: 1fwd -> snarl -> 2fwd
        {ET::BEGIN_CHAIN, 1, false},
          {ET::BEGIN_SNARL, 1, false},
          {ET::END_SNARL, 2, false},
        {ET::END_CHAIN, 2, false},
        // Chain 2: 3fwd -> snarl -> 4fwd
        {ET::BEGIN_CHAIN, 3, false},
          {ET::BEGIN_SNARL, 3, false},
          {ET::END_SNARL, 4, false},
        {ET::END_CHAIN, 4, false},
    };

    ReplaySnarlFinder replay(&graph, events);

    SECTION("flip only first chain") {
        std::unordered_set<nid_t> flips {1, 2};
        SnarlDecompositionFuzzer fuzzer(&graph, &replay, flips);

        std::vector<Event> captured = capture_events(fuzzer, graph);

        std::vector<Event> expected = {
            {ET::BEGIN_CHAIN, 2, true},
              {ET::BEGIN_SNARL, 2, true},
              {ET::END_SNARL, 1, true},
            {ET::END_CHAIN, 1, true},
            {ET::BEGIN_CHAIN, 3, false},
              {ET::BEGIN_SNARL, 3, false},
              {ET::END_SNARL, 4, false},
            {ET::END_CHAIN, 4, false},
        };

        REQUIRE(captured == expected);
    }

    SECTION("flip both chains") {
        std::unordered_set<nid_t> flips {1, 2, 3, 4};
        SnarlDecompositionFuzzer fuzzer(&graph, &replay, flips);

        std::vector<Event> captured = capture_events(fuzzer, graph);

        std::vector<Event> expected = {
            {ET::BEGIN_CHAIN, 2, true},
              {ET::BEGIN_SNARL, 2, true},
              {ET::END_SNARL, 1, true},
            {ET::END_CHAIN, 1, true},
            {ET::BEGIN_CHAIN, 4, true},
              {ET::BEGIN_SNARL, 4, true},
              {ET::END_SNARL, 3, true},
            {ET::END_CHAIN, 3, true},
        };

        REQUIRE(captured == expected);
    }
}

TEST_CASE("SnarlDecompositionFuzzer handles deeply nested chains", "[snarl_decomposition_fuzzer]") {
    bdsg::HashGraph graph;
    for (nid_t i = 1; i <= 8; i++) {
        graph.create_handle("A", i);
    }

    // Outer chain: 1->6
    //   Snarl(1,4)
    //     Inner chain: 2->3
    //       Snarl(2,3) [leaf snarl, no children]
    //   Snarl(4,6)
    //     Inner chain: 5->5 [empty/trivial]
    std::vector<Event> events = {
        {ET::BEGIN_CHAIN, 1, false},
          {ET::BEGIN_SNARL, 1, false},
            {ET::BEGIN_CHAIN, 2, false},
              {ET::BEGIN_SNARL, 2, false},
              {ET::END_SNARL, 3, false},
            {ET::END_CHAIN, 3, false},
          {ET::END_SNARL, 4, false},
          {ET::BEGIN_SNARL, 4, false},
            {ET::BEGIN_CHAIN, 5, false},
            {ET::END_CHAIN, 5, false},
          {ET::END_SNARL, 6, false},
        {ET::END_CHAIN, 6, false},
    };

    ReplaySnarlFinder replay(&graph, events);

    SECTION("flip outer chain only") {
        std::unordered_set<nid_t> flips {1, 6};
        SnarlDecompositionFuzzer fuzzer(&graph, &replay, flips);

        std::vector<Event> captured = capture_events(fuzzer, graph);

        // Inner chain and its snarls should flip too. 
        std::vector<Event> expected = {
            {ET::BEGIN_CHAIN, 6, true},
              {ET::BEGIN_SNARL, 6, true},
                {ET::BEGIN_CHAIN, 5, true},
                {ET::END_CHAIN, 5, true},
              {ET::END_SNARL, 4, true},
              {ET::BEGIN_SNARL, 4, true},
                {ET::BEGIN_CHAIN, 3, true},
                  {ET::BEGIN_SNARL, 3, true},
                  {ET::END_SNARL, 2, true},
                {ET::END_CHAIN, 2, true},
              {ET::END_SNARL, 1, true},
            {ET::END_CHAIN, 1, true},
        };

        REQUIRE(captured == expected);
    }

    SECTION("flip outer and inner chain") {
        std::unordered_set<nid_t> flips {1, 6, 2, 3};
        SnarlDecompositionFuzzer fuzzer(&graph, &replay, flips);

        std::vector<Event> captured = capture_events(fuzzer, graph);

        // Outer chain should flip but inner chain should flip back
        std::vector<Event> expected = {
            {ET::BEGIN_CHAIN, 6, true},
              {ET::BEGIN_SNARL, 6, true},
                {ET::BEGIN_CHAIN, 5, true},
                {ET::END_CHAIN, 5, true},
              {ET::END_SNARL, 4, true},
              {ET::BEGIN_SNARL, 4, true},
                {ET::BEGIN_CHAIN, 2, false},
                  {ET::BEGIN_SNARL, 2, false},
                  {ET::END_SNARL, 3, false},
                {ET::END_CHAIN, 3, false},
              {ET::END_SNARL, 1, true},
            {ET::END_CHAIN, 1, true},
        };

        REQUIRE(captured == expected);
    }
}

} // namespace unittest
} // namespace vg
