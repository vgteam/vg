#include "catch.hpp"
#include "../handle.hpp"
#include <bdsg/hash_graph.hpp>

#include "support/snarl_decomposition_fuzzer.hpp"

#include <vector>
#include <set>

namespace vg {
namespace unittest {

using namespace std;
using ET = ReplaySnarlFinder::EventType;
using Event = ReplaySnarlFinder::Event;

/// Capture all events emitted by a traverse_decomposition call into a vector.
static vector<pair<ET, handle_t>> capture_events(const HandleGraphSnarlFinder& finder) {
    vector<pair<ET, handle_t>> result;
    finder.traverse_decomposition(
        [&](handle_t h) { result.push_back({ET::BEGIN_CHAIN, h}); },
        [&](handle_t h) { result.push_back({ET::END_CHAIN, h}); },
        [&](handle_t h) { result.push_back({ET::BEGIN_SNARL, h}); },
        [&](handle_t h) { result.push_back({ET::END_SNARL, h}); }
    );
    return result;
}

TEST_CASE("ReplaySnarlFinder replays events faithfully", "[snarl_decomposition_fuzzer]") {
    // Build a small graph to get real handles
    bdsg::HashGraph graph;
    graph.create_handle("A", 10);
    graph.create_handle("C", 12);
    graph.create_handle("G", 15);
    graph.create_handle("T", 20);
    graph.create_handle("AA", 22);

    handle_t h10f = graph.get_handle(10, false);
    handle_t h12r = graph.get_handle(12, true);
    handle_t h15r = graph.get_handle(15, true);
    handle_t h20f = graph.get_handle(20, false);
    handle_t h22f = graph.get_handle(22, false);

    vector<Event> events = {
        {ET::BEGIN_CHAIN, h10f},
        {ET::BEGIN_SNARL, h10f},
        {ET::BEGIN_CHAIN, h12r},
        {ET::END_CHAIN, h15r},
        {ET::END_SNARL, h20f},
        {ET::BEGIN_SNARL, h20f},
        {ET::END_SNARL, h22f},
        {ET::END_CHAIN, h22f},
    };

    ReplaySnarlFinder finder(events);
    auto captured = capture_events(finder);

    REQUIRE(captured.size() == events.size());
    for (size_t i = 0; i < events.size(); i++) {
        REQUIRE(captured[i].first == events[i].type);
        REQUIRE(captured[i].second == events[i].handle);
    }
}

TEST_CASE("SnarlDecompositionFuzzer passes through when nothing is flipped", "[snarl_decomposition_fuzzer]") {
    bdsg::HashGraph graph;
    graph.create_handle("A", 10);
    graph.create_handle("C", 12);
    graph.create_handle("G", 15);
    graph.create_handle("T", 20);
    graph.create_handle("AA", 22);

    handle_t h10f = graph.get_handle(10, false);
    handle_t h12r = graph.get_handle(12, true);
    handle_t h15r = graph.get_handle(15, true);
    handle_t h20f = graph.get_handle(20, false);
    handle_t h22f = graph.get_handle(22, false);

    vector<Event> events = {
        {ET::BEGIN_CHAIN, h10f},
        {ET::BEGIN_SNARL, h10f},
        {ET::BEGIN_CHAIN, h12r},
        {ET::END_CHAIN, h15r},
        {ET::END_SNARL, h20f},
        {ET::BEGIN_SNARL, h20f},
        {ET::END_SNARL, h22f},
        {ET::END_CHAIN, h22f},
    };

    ReplaySnarlFinder replay(events);

    // No chains to flip
    set<pair<handle_t, handle_t>> no_flips;
    SnarlDecompositionFuzzer fuzzer(&graph, &replay, no_flips);

    auto captured = capture_events(fuzzer);

    REQUIRE(captured.size() == events.size());
    for (size_t i = 0; i < events.size(); i++) {
        REQUIRE(captured[i].first == events[i].type);
        REQUIRE(captured[i].second == events[i].handle);
    }
}

TEST_CASE("SnarlDecompositionFuzzer flips an outer chain", "[snarl_decomposition_fuzzer]") {
    // Graph:
    // Chain: 10fwd -> snarl(10fwd, 20fwd) -> snarl(20fwd, 22fwd) -> 22fwd
    // Inside first snarl: chain 12rev->15rev
    bdsg::HashGraph graph;
    graph.create_handle("A", 10);
    graph.create_handle("C", 12);
    graph.create_handle("G", 15);
    graph.create_handle("T", 20);
    graph.create_handle("AA", 22);

    handle_t h10f = graph.get_handle(10, false);
    handle_t h10r = graph.get_handle(10, true);
    handle_t h12r = graph.get_handle(12, true);
    handle_t h12f = graph.get_handle(12, false);
    handle_t h15r = graph.get_handle(15, true);
    handle_t h15f = graph.get_handle(15, false);
    handle_t h20f = graph.get_handle(20, false);
    handle_t h20r = graph.get_handle(20, true);
    handle_t h22f = graph.get_handle(22, false);
    handle_t h22r = graph.get_handle(22, true);

    vector<Event> events = {
        {ET::BEGIN_CHAIN, h10f},
          {ET::BEGIN_SNARL, h10f},
            {ET::BEGIN_CHAIN, h12r},
            {ET::END_CHAIN, h15r},
          {ET::END_SNARL, h20f},
          {ET::BEGIN_SNARL, h20f},
          {ET::END_SNARL, h22f},
        {ET::END_CHAIN, h22f},
    };

    ReplaySnarlFinder replay(events);

    SECTION("flip outer chain only") {
        // Flip the outer chain (10fwd -> 22fwd)
        set<pair<handle_t, handle_t>> flips = {{h10f, h22f}};
        SnarlDecompositionFuzzer fuzzer(&graph, &replay, flips);

        auto captured = capture_events(fuzzer);

        // Expected after flipping the outer chain:
        // begin_chain(22 rev)
        //   begin_snarl(22 rev)      -- was end_snarl(22 fwd), flipped
        //   end_snarl(20 rev)        -- was begin_snarl(20 fwd), flipped
        //   begin_snarl(20 rev)      -- was end_snarl(20 fwd), flipped
        //     begin_chain(12 rev)    -- NOT flipped (not in flip set)
        //     end_chain(15 rev)
        //   end_snarl(10 rev)        -- was begin_snarl(10 fwd), flipped
        // end_chain(10 rev)

        vector<pair<ET, handle_t>> expected = {
            {ET::BEGIN_CHAIN, h22r},
              {ET::BEGIN_SNARL, h22r},
              {ET::END_SNARL, h20r},
              {ET::BEGIN_SNARL, h20r},
                {ET::BEGIN_CHAIN, h12r},
                {ET::END_CHAIN, h15r},
              {ET::END_SNARL, h10r},
            {ET::END_CHAIN, h10r},
        };

        REQUIRE(captured.size() == expected.size());
        for (size_t i = 0; i < expected.size(); i++) {
            REQUIRE(captured[i].first == expected[i].first);
            REQUIRE(captured[i].second == expected[i].second);
        }
    }

    SECTION("flip outer and nested chain") {
        // Flip outer chain (10fwd->22fwd) AND nested chain (12rev->15rev)
        set<pair<handle_t, handle_t>> flips = {{h10f, h22f}, {h12r, h15r}};
        SnarlDecompositionFuzzer fuzzer(&graph, &replay, flips);

        auto captured = capture_events(fuzzer);

        // Expected: outer chain flipped, AND the nested chain is also flipped
        // The nested chain 12rev->15rev becomes: begin_chain(15fwd), end_chain(12fwd)
        vector<pair<ET, handle_t>> expected = {
            {ET::BEGIN_CHAIN, h22r},
              {ET::BEGIN_SNARL, h22r},
              {ET::END_SNARL, h20r},
              {ET::BEGIN_SNARL, h20r},
                {ET::BEGIN_CHAIN, h15f},
                {ET::END_CHAIN, h12f},
              {ET::END_SNARL, h10r},
            {ET::END_CHAIN, h10r},
        };

        REQUIRE(captured.size() == expected.size());
        for (size_t i = 0; i < expected.size(); i++) {
            REQUIRE(captured[i].first == expected[i].first);
            REQUIRE(captured[i].second == expected[i].second);
        }
    }

    SECTION("flip nested chain only") {
        // Flip only the nested chain (12rev->15rev), outer stays
        set<pair<handle_t, handle_t>> flips = {{h12r, h15r}};
        SnarlDecompositionFuzzer fuzzer(&graph, &replay, flips);

        auto captured = capture_events(fuzzer);

        // Outer chain not flipped, nested chain flipped
        vector<pair<ET, handle_t>> expected = {
            {ET::BEGIN_CHAIN, h10f},
              {ET::BEGIN_SNARL, h10f},
                {ET::BEGIN_CHAIN, h15f},
                {ET::END_CHAIN, h12f},
              {ET::END_SNARL, h20f},
              {ET::BEGIN_SNARL, h20f},
              {ET::END_SNARL, h22f},
            {ET::END_CHAIN, h22f},
        };

        REQUIRE(captured.size() == expected.size());
        for (size_t i = 0; i < expected.size(); i++) {
            REQUIRE(captured[i].first == expected[i].first);
            REQUIRE(captured[i].second == expected[i].second);
        }
    }
}

TEST_CASE("SnarlDecompositionFuzzer handles empty chain (no snarls)", "[snarl_decomposition_fuzzer]") {
    bdsg::HashGraph graph;
    graph.create_handle("ACGT", 5);

    handle_t h5f = graph.get_handle(5, false);
    handle_t h5r = graph.get_handle(5, true);

    // An empty chain: begin and end with same handle, no snarls inside
    vector<Event> events = {
        {ET::BEGIN_CHAIN, h5f},
        {ET::END_CHAIN, h5f},
    };

    ReplaySnarlFinder replay(events);

    SECTION("flipping an empty chain") {
        set<pair<handle_t, handle_t>> flips = {{h5f, h5f}};
        SnarlDecompositionFuzzer fuzzer(&graph, &replay, flips);

        auto captured = capture_events(fuzzer);

        // Flipped: begin_chain(flip(5fwd)) = begin_chain(5rev)
        //          end_chain(flip(5fwd)) = end_chain(5rev)
        vector<pair<ET, handle_t>> expected = {
            {ET::BEGIN_CHAIN, h5r},
            {ET::END_CHAIN, h5r},
        };

        REQUIRE(captured.size() == expected.size());
        for (size_t i = 0; i < expected.size(); i++) {
            REQUIRE(captured[i].first == expected[i].first);
            REQUIRE(captured[i].second == expected[i].second);
        }
    }
}

TEST_CASE("SnarlDecompositionFuzzer handles multiple top-level chains", "[snarl_decomposition_fuzzer]") {
    bdsg::HashGraph graph;
    graph.create_handle("A", 1);
    graph.create_handle("C", 2);
    graph.create_handle("G", 3);
    graph.create_handle("T", 4);

    handle_t h1f = graph.get_handle(1, false);
    handle_t h1r = graph.get_handle(1, true);
    handle_t h2f = graph.get_handle(2, false);
    handle_t h2r = graph.get_handle(2, true);
    handle_t h3f = graph.get_handle(3, false);
    handle_t h3r = graph.get_handle(3, true);
    handle_t h4f = graph.get_handle(4, false);
    handle_t h4r = graph.get_handle(4, true);

    // Two top-level chains in the root snarl
    vector<Event> events = {
        // Chain 1: 1fwd -> snarl -> 2fwd
        {ET::BEGIN_CHAIN, h1f},
          {ET::BEGIN_SNARL, h1f},
          {ET::END_SNARL, h2f},
        {ET::END_CHAIN, h2f},
        // Chain 2: 3fwd -> snarl -> 4fwd
        {ET::BEGIN_CHAIN, h3f},
          {ET::BEGIN_SNARL, h3f},
          {ET::END_SNARL, h4f},
        {ET::END_CHAIN, h4f},
    };

    ReplaySnarlFinder replay(events);

    SECTION("flip only first chain") {
        set<pair<handle_t, handle_t>> flips = {{h1f, h2f}};
        SnarlDecompositionFuzzer fuzzer(&graph, &replay, flips);

        auto captured = capture_events(fuzzer);

        vector<pair<ET, handle_t>> expected = {
            {ET::BEGIN_CHAIN, h2r},
              {ET::BEGIN_SNARL, h2r},
              {ET::END_SNARL, h1r},
            {ET::END_CHAIN, h1r},
            {ET::BEGIN_CHAIN, h3f},
              {ET::BEGIN_SNARL, h3f},
              {ET::END_SNARL, h4f},
            {ET::END_CHAIN, h4f},
        };

        REQUIRE(captured.size() == expected.size());
        for (size_t i = 0; i < expected.size(); i++) {
            REQUIRE(captured[i].first == expected[i].first);
            REQUIRE(captured[i].second == expected[i].second);
        }
    }

    SECTION("flip both chains") {
        set<pair<handle_t, handle_t>> flips = {{h1f, h2f}, {h3f, h4f}};
        SnarlDecompositionFuzzer fuzzer(&graph, &replay, flips);

        auto captured = capture_events(fuzzer);

        vector<pair<ET, handle_t>> expected = {
            {ET::BEGIN_CHAIN, h2r},
              {ET::BEGIN_SNARL, h2r},
              {ET::END_SNARL, h1r},
            {ET::END_CHAIN, h1r},
            {ET::BEGIN_CHAIN, h4r},
              {ET::BEGIN_SNARL, h4r},
              {ET::END_SNARL, h3r},
            {ET::END_CHAIN, h3r},
        };

        REQUIRE(captured.size() == expected.size());
        for (size_t i = 0; i < expected.size(); i++) {
            REQUIRE(captured[i].first == expected[i].first);
            REQUIRE(captured[i].second == expected[i].second);
        }
    }
}

TEST_CASE("SnarlDecompositionFuzzer handles deeply nested chains", "[snarl_decomposition_fuzzer]") {
    bdsg::HashGraph graph;
    for (nid_t i = 1; i <= 8; i++) {
        graph.create_handle("A", i);
    }

    handle_t h1f = graph.get_handle(1, false);
    handle_t h1r = graph.get_handle(1, true);
    handle_t h2f = graph.get_handle(2, false);
    handle_t h2r = graph.get_handle(2, true);
    handle_t h3f = graph.get_handle(3, false);
    handle_t h3r = graph.get_handle(3, true);
    handle_t h4f = graph.get_handle(4, false);
    handle_t h4r = graph.get_handle(4, true);
    handle_t h5f = graph.get_handle(5, false);
    handle_t h5r = graph.get_handle(5, true);
    handle_t h6f = graph.get_handle(6, false);
    handle_t h6r = graph.get_handle(6, true);

    // Outer chain: 1->6
    //   Snarl(1,4)
    //     Inner chain: 2->3
    //       Snarl(2,3) [leaf snarl, no children]
    //   Snarl(4,6)
    //     Inner chain: 5->5 [empty/trivial]
    vector<Event> events = {
        {ET::BEGIN_CHAIN, h1f},
          {ET::BEGIN_SNARL, h1f},
            {ET::BEGIN_CHAIN, h2f},
              {ET::BEGIN_SNARL, h2f},
              {ET::END_SNARL, h3f},
            {ET::END_CHAIN, h3f},
          {ET::END_SNARL, h4f},
          {ET::BEGIN_SNARL, h4f},
            {ET::BEGIN_CHAIN, h5f},
            {ET::END_CHAIN, h5f},
          {ET::END_SNARL, h6f},
        {ET::END_CHAIN, h6f},
    };

    ReplaySnarlFinder replay(events);

    SECTION("flip outer chain only") {
        set<pair<handle_t, handle_t>> flips = {{h1f, h6f}};
        SnarlDecompositionFuzzer fuzzer(&graph, &replay, flips);

        auto captured = capture_events(fuzzer);

        // Outer flipped: snarls reversed, but inner chains NOT flipped
        vector<pair<ET, handle_t>> expected = {
            {ET::BEGIN_CHAIN, h6r},
              {ET::BEGIN_SNARL, h6r},
                {ET::BEGIN_CHAIN, h5f},
                {ET::END_CHAIN, h5f},
              {ET::END_SNARL, h4r},
              {ET::BEGIN_SNARL, h4r},
                {ET::BEGIN_CHAIN, h2f},
                  {ET::BEGIN_SNARL, h2f},
                  {ET::END_SNARL, h3f},
                {ET::END_CHAIN, h3f},
              {ET::END_SNARL, h1r},
            {ET::END_CHAIN, h1r},
        };

        REQUIRE(captured.size() == expected.size());
        for (size_t i = 0; i < expected.size(); i++) {
            REQUIRE(captured[i].first == expected[i].first);
            REQUIRE(captured[i].second == expected[i].second);
        }
    }

    SECTION("flip outer and inner chain 2->3") {
        set<pair<handle_t, handle_t>> flips = {{h1f, h6f}, {h2f, h3f}};
        SnarlDecompositionFuzzer fuzzer(&graph, &replay, flips);

        auto captured = capture_events(fuzzer);

        // Outer flipped, inner chain 2->3 also flipped
        vector<pair<ET, handle_t>> expected = {
            {ET::BEGIN_CHAIN, h6r},
              {ET::BEGIN_SNARL, h6r},
                {ET::BEGIN_CHAIN, h5f},
                {ET::END_CHAIN, h5f},
              {ET::END_SNARL, h4r},
              {ET::BEGIN_SNARL, h4r},
                {ET::BEGIN_CHAIN, h3r},
                  {ET::BEGIN_SNARL, h3r},
                  {ET::END_SNARL, h2r},
                {ET::END_CHAIN, h2r},
              {ET::END_SNARL, h1r},
            {ET::END_CHAIN, h1r},
        };

        REQUIRE(captured.size() == expected.size());
        for (size_t i = 0; i < expected.size(); i++) {
            REQUIRE(captured[i].first == expected[i].first);
            REQUIRE(captured[i].second == expected[i].second);
        }
    }
}

} // namespace unittest
} // namespace vg
