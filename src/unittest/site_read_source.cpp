/// \file unittest/site_read_source.cpp
///
/// Unit tests for the windowing and caching shared by the on-demand read sources.
///
/// The two on-demand backends -- indexed GAM and GAF-Base -- differ only in how they
/// fetch a span of node IDs. Everything above that is shared, and everything above that
/// is where the subtle mistakes live: a query that straddles a window boundary, a read
/// returned by a window fetch that no site in the window actually wants, a cache that
/// answers with the wrong window's reads. Those are tested here against a fake backend
/// that records what was asked for, so they are tested without a GAM index, without a
/// GAF-Base, and without the gbz-base binary.
///
/// The point of the fake is that it can assert on *fetches*, which the real backends
/// cannot cheaply be asked about: the tests below pin how many times the backend was
/// hit and over what range, not just which reads came back.
///

#include <algorithm>
#include <map>
#include <mutex>
#include <stdexcept>
#include <string>
#include <tuple>
#include <vector>

#include "catch.hpp"
#include "site_read_source.hpp"

namespace vg {
namespace unittest {

using namespace std;

/// A read on the given nodes, in order, one base on each, with the given mapping quality.
static Alignment read_on(const string& name, const vector<nid_t>& nodes, int32_t mapq = 60) {
    Alignment aln;
    aln.set_name(name);
    aln.set_mapping_quality(mapq);
    aln.set_sequence(string(nodes.size(), 'A'));
    for (nid_t node : nodes) {
        auto* mapping = aln.mutable_path()->add_mapping();
        mapping->mutable_position()->set_node_id(node);
        auto* edit = mapping->add_edit();
        edit->set_from_length(1);
        edit->set_to_length(1);
    }
    return aln;
}

/// A read source over reads held in a vector, recording every span fetched.
///
/// Most reads are described by the single node they sit on, which is all the windowing logic
/// looks at; add() takes a read on several nodes. Sequence and quality are irrelevant except to
/// the read-start counts.
class FakeWindowedSource : public WindowedSiteReadSource {
public:
    FakeWindowedSource(const vector<pair<string, nid_t>>& reads, size_t window_size,
                       size_t cache_entries = 2)
        : WindowedSiteReadSource(SiteReadFilter(), window_size, cache_entries) {
        for (const auto& read : reads) {
            Alignment aln;
            aln.set_name(read.first);
            aln.set_mapping_quality(60);
            auto* mapping = aln.mutable_path()->add_mapping();
            mapping->mutable_position()->set_node_id(read.second);
            auto* edit = mapping->add_edit();
            edit->set_from_length(1);
            edit->set_to_length(1);
            held.push_back(aln);
        }
    }

    /// Hold another read, such as one on several nodes.
    void add(const Alignment& aln) {
        held.push_back(aln);
    }

    /// Make the next fetch throw, as a backend whose query fails would.
    void fail_next_fetch() {
        fail_next = true;
    }

    /// The (min, max) extent of every fetch_span call, in order. Recorded as an
    /// extent rather than the full range list because that is what the windowing
    /// assertions are about; get_fetch_ranges() has the detail.
    const vector<pair<nid_t, nid_t>>& get_fetches() const {
        return fetches;
    }

    /// Every range list fetch_span was called with, in order.
    const vector<vector<pair<nid_t, nid_t>>>& get_fetch_ranges() const {
        return fetch_ranges;
    }

protected:

    void fetch_span(const vector<pair<nid_t, nid_t>>& ranges,
                    const function<void(Alignment&)>& iteratee) const {
        nid_t min_id = ranges.front().first;
        nid_t max_id = ranges.front().second;
        for (const auto& range : ranges) {
            min_id = min(min_id, range.first);
            max_id = max(max_id, range.second);
        }
        {
            // Threads may fetch at once.
            lock_guard<std::mutex> guard(record_mutex);
            fetch_ranges.push_back(ranges);
            fetches.push_back(make_pair(min_id, max_id));
            if (fail_next) {
                fail_next = false;
                throw runtime_error("injected fetch failure");
            }
        }
        for (const Alignment& aln : held) {
            // A real backend returns a read touching the ranges anywhere along its path.
            bool in_range = false;
            for (const auto& mapping : aln.path().mapping()) {
                nid_t node = mapping.position().node_id();
                for (const auto& range : ranges) {
                    if (node >= range.first && node <= range.second) {
                        in_range = true;
                        break;
                    }
                }
            }
            if (in_range) {
                if (passes_filter(aln)) {
                    count_fetched();
                    // A copy, because `held` is this fake's permanent store and a
                    // caller is entitled to move from what fetch_span hands it. Real
                    // backends hand over a per-record scratch alignment instead. If
                    // this were `iteratee(aln)` the second fetch of a window would
                    // return empty reads, which is what the repeat-fetch cases check.
                    Alignment owned = aln;
                    iteratee(owned);
                }
            }
        }
    }

private:
    vector<Alignment> held;
    mutable vector<pair<nid_t, nid_t>> fetches;
    mutable vector<vector<pair<nid_t, nid_t>>> fetch_ranges;
    mutable bool fail_next = false;
    mutable std::mutex record_mutex;
};

/// Collect the names of the reads a query returns.
static vector<string> names_for(const SiteReadSource& source,
                                const vector<pair<nid_t, nid_t>>& ranges) {
    vector<string> names;
    source.for_each_alignment(ranges, [&](const Alignment& aln) {
        names.push_back(aln.name());
    });
    return names;
}

TEST_CASE("A fetch is quantised to a whole window, not the range asked for",
          "[site_read_source]") {
    // Window 100 means node 30's window is [0, 99], whatever narrow range was requested.
    FakeWindowedSource source({{"a", 30}}, 100);

    names_for(source, {{30, 31}});

    REQUIRE(source.get_fetches().size() == 1);
    REQUIRE(source.get_fetches()[0].first == 0);
    REQUIRE(source.get_fetches()[0].second == 99);
}

TEST_CASE("Reads in the window but outside the requested ranges are not returned",
          "[site_read_source]") {
    // This is the property that makes window fetching safe: the window is a
    // performance unit, not a change to what the query means. Without the narrowing,
    // a site would be handed reads from anywhere in its window.
    FakeWindowedSource source({{"wanted", 30}, {"same_window", 80}}, 100);

    vector<string> names = names_for(source, {{25, 35}});

    REQUIRE(names.size() == 1);
    REQUIRE(names[0] == "wanted");
}

TEST_CASE("A second query in the same window is served from the cache",
          "[site_read_source]") {
    FakeWindowedSource source({{"a", 30}, {"b", 40}}, 100);

    names_for(source, {{30, 31}});
    names_for(source, {{40, 41}});

    // One fetch, two queries: the second was answered from the first's reads.
    REQUIRE(source.get_fetches().size() == 1);
    REQUIRE(source.get_cache_hits() == 1);
    REQUIRE(source.get_cache_misses() == 1);
}

TEST_CASE("A cache hit still returns the right reads, not the whole window",
          "[site_read_source]") {
    // A cache that returned its whole window would pass the counting test above while
    // handing every site the wrong evidence, so check the reads too.
    FakeWindowedSource source({{"a", 30}, {"b", 40}}, 100);

    names_for(source, {{30, 31}});
    vector<string> names = names_for(source, {{40, 41}});

    REQUIRE(source.get_cache_hits() == 1);
    REQUIRE(names.size() == 1);
    REQUIRE(names[0] == "b");
}

TEST_CASE("Visiting windows in order fetches each exactly once", "[site_read_source]") {
    // The access pattern GraphCaller::set_node_id_ordering exists to produce. With two
    // cache slots and ascending visits, no window is ever fetched twice.
    FakeWindowedSource source({{"a", 10}, {"b", 110}, {"c", 210}, {"d", 310}}, 100);

    for (nid_t node : {10, 110, 210, 310}) {
        names_for(source, {{node, node}});
    }

    REQUIRE(source.get_fetches().size() == 4);
    REQUIRE(source.get_cache_misses() == 4);
    REQUIRE(source.get_cache_hits() == 0);
}

TEST_CASE("Revisiting a window beyond the cache's depth refetches it",
          "[site_read_source]") {
    // Two slots, three windows visited, then back to the first: it has been evicted.
    // Pinned so that the cost of an out-of-order visit is a known quantity rather than
    // a surprise, since that cost is the whole reason ordering was worth arranging.
    FakeWindowedSource source({{"a", 10}, {"b", 110}, {"c", 210}}, 100, 2);

    names_for(source, {{10, 10}});
    names_for(source, {{110, 110}});
    names_for(source, {{210, 210}});
    vector<string> names = names_for(source, {{10, 10}});

    REQUIRE(source.get_fetches().size() == 4);
    REQUIRE(names.size() == 1);
    REQUIRE(names[0] == "a");
}

TEST_CASE("A query straddling a window boundary fetches the exact span and does not cache",
          "[site_read_source]") {
    // Stitching windows together would mean de-duplicating reads that span the
    // boundary; fetching the span directly avoids that, since a backend emits each read
    // at most once per fetch.
    FakeWindowedSource source({{"a", 90}, {"b", 110}}, 100);

    vector<string> names = names_for(source, {{90, 110}});

    REQUIRE(source.get_fetches().size() == 1);
    REQUIRE(source.get_fetches()[0].first == 90);
    REQUIRE(source.get_fetches()[0].second == 110);
    REQUIRE(names.size() == 2);
    // Counted as a miss, and nothing was cached, so a repeat costs another fetch.
    REQUIRE(source.get_cache_hits() == 0);
    names_for(source, {{90, 110}});
    REQUIRE(source.get_fetches().size() == 2);
}

TEST_CASE("Several ranges are covered by one window fetch when they share a window",
          "[site_read_source]") {
    // A snarl's node set can arrive as several disjoint ranges. Their overall extent is
    // what picks the window, so ranges inside one window cost one fetch.
    FakeWindowedSource source({{"a", 10}, {"b", 50}, {"skipped", 30}}, 100);

    vector<string> names = names_for(source, {{5, 15}, {45, 55}});

    REQUIRE(source.get_fetches().size() == 1);
    REQUIRE(names.size() == 2);
    REQUIRE(names[0] == "a");
    REQUIRE(names[1] == "b");
}

TEST_CASE("An empty range list fetches nothing", "[site_read_source]") {
    FakeWindowedSource source({{"a", 10}}, 100);

    vector<string> names = names_for(source, {});

    REQUIRE(names.empty());
    REQUIRE(source.get_fetches().empty());
}

TEST_CASE("A window with no reads is still cached, so it is not fetched twice",
          "[site_read_source]") {
    // An empty result is a result. Refetching empty windows would make sparse regions
    // -- exactly where windowing is already a poor fit -- worse again.
    FakeWindowedSource source({{"far", 500}}, 100);

    names_for(source, {{10, 20}});
    names_for(source, {{30, 40}});

    REQUIRE(source.get_fetches().size() == 1);
    REQUIRE(source.get_cache_hits() == 1);
}

TEST_CASE("A read touching any node of a multi-node range is returned once",
          "[site_read_source]") {
    // A range spanning several nodes must yield each touching read exactly once, not once
    // per node it touches -- double-counting here would inflate every allele's read support.
    FakeWindowedSource source({{"a", 10}, {"b", 20}}, 100);

    REQUIRE(names_for(source, {{10, 20}}).size() == 2);
}

TEST_CASE("Paired mates sharing a read name are both kept", "[site_read_source]") {
    // Paired mates in GAF share a name, so de-duplication keyed on the name alone would drop one
    // mate of every pair. The de-duplication that must get this right is in GafBaseSiteReadSource's
    // split-query path, which needs the gbz-base binary; this checks the property it must keep: two
    // same-named reads at different positions are two reads.
    FakeWindowedSource source({{"pair", 10}, {"pair", 20}}, 100);

    vector<string> names = names_for(source, {{10, 20}});

    REQUIRE(names.size() == 2);
    REQUIRE(names[0] == "pair");
    REQUIRE(names[1] == "pair");
}

TEST_CASE("The MAPQ filter is applied by the base class, not left to each backend",
          "[site_read_source]") {
    // A filter applied in one backend and forgotten in the other would mean the two
    // paths genotyped different read sets, so it lives in one place and is checked here.
    SiteReadFilter filter;
    filter.min_mapq = 30;

    class LowMapqSource : public WindowedSiteReadSource {
    public:
        LowMapqSource(const SiteReadFilter& filter)
            : WindowedSiteReadSource(filter, 100, 2) {}
    protected:
        void fetch_span(const vector<pair<nid_t, nid_t>>& ranges,
                        const function<void(Alignment&)>& iteratee) const {
            for (int mapq : {10, 60}) {
                Alignment aln;
                aln.set_name("q" + to_string(mapq));
                aln.set_mapping_quality(mapq);
                auto* mapping = aln.mutable_path()->add_mapping();
                mapping->mutable_position()->set_node_id(10);
                if (passes_filter(aln)) {
                    count_fetched();
                    iteratee(aln);
                }
            }
        }
    };

    LowMapqSource source(filter);
    vector<string> names = names_for(source, {{10, 10}});

    REQUIRE(names.size() == 1);
    REQUIRE(names[0] == "q60");
    REQUIRE(source.get_filtered_count() == 1);
    REQUIRE(source.get_read_count() == 1);
}

/// The read starts a source reports on the nodes, summed by mapping quality, as
/// (reads, bases) per MAPQ. Summed because a source may group the starts any way it likes.
static map<int32_t, pair<size_t, size_t>> starts_on(const SiteReadSource& source,
                                                    const vector<nid_t>& nodes) {
    map<int32_t, pair<size_t, size_t>> by_mapq;
    source.for_each_read_start(nodes, [&](int32_t mapq, size_t reads, size_t bases) {
        by_mapq[mapq].first += reads;
        by_mapq[mapq].second += bases;
    });
    return by_mapq;
}

TEST_CASE("A windowed source counts read starts as the default does", "[site_read_source]") {
    // x begins in window 0 and runs into window 1; y and z begin in window 1; w begins
    // on a node not asked about. Each read must be counted once, on its first node, whichever
    // windows its path touches.
    // Node 105 also starts p at MAPQ 30 and q and r at MAPQ 60, and node 150 starts s, so
    // a node can hold several reads and several mapping qualities.
    vector<Alignment> reads{read_on("x", {95, 105}), read_on("y", {105, 106}),
                            read_on("z", {150}, 30), read_on("w", {30}),
                            read_on("p", {105}, 30), read_on("q", {105, 106, 107}),
                            read_on("r", {105}), read_on("s", {150, 151}, 30)};
    FakeWindowedSource windowed({}, 100);
    InMemorySiteReadSource in_memory;
    for (const Alignment& aln : reads) {
        windowed.add(aln);
        in_memory.add(aln);
    }

    vector<nid_t> nodes{95, 105, 150};
    auto expected = starts_on(in_memory, nodes);
    REQUIRE(expected.size() == 2);
    REQUIRE(expected[60] == make_pair<size_t, size_t>(4, 8));
    REQUIRE(expected[30] == make_pair<size_t, size_t>(3, 4));
    REQUIRE(starts_on(windowed, nodes) == expected);

    // The groups the windowed source reports on node 105 keep its two mapping qualities apart.
    vector<tuple<int32_t, size_t, size_t>> groups;
    windowed.for_each_read_start({105}, [&](int32_t mapq, size_t count, size_t bases) {
        groups.emplace_back(mapq, count, bases);
    });
    std::sort(groups.begin(), groups.end());
    REQUIRE(groups == vector<tuple<int32_t, size_t, size_t>>{{30, 1, 1}, {60, 3, 6}});
}

TEST_CASE("A read start is counted on the read's first node only", "[site_read_source]") {
    // x passes through node 105 but begins on node 95, so it is not a start on 105.
    FakeWindowedSource windowed({}, 100);
    windowed.add(read_on("x", {95, 105}));
    windowed.add(read_on("y", {105}));

    auto counts = starts_on(windowed, {105});

    REQUIRE(counts.size() == 1);
    REQUIRE(counts[60].first == 1);
}

TEST_CASE("Counting read starts fetches only windows never fetched whole", "[site_read_source]") {
    // Read starts are counted over spans wider than a site. Window 0 is already fetched for a
    // site, so counting fetches window 1 alone, whole, into the cache, where a later site in
    // window 1 finds it.
    FakeWindowedSource windowed({{"a", 30}, {"b", 95}, {"c", 105}}, 100);

    names_for(windowed, {{30, 30}});
    REQUIRE(windowed.get_fetches().size() == 1);

    auto counts = starts_on(windowed, {95, 105});
    REQUIRE(counts[60].first == 2);
    REQUIRE(windowed.get_fetches().size() == 2);
    REQUIRE(windowed.get_fetches()[1] == make_pair<nid_t, nid_t>(100, 199));

    // Counted again, nothing is fetched.
    starts_on(windowed, {95, 105});
    REQUIRE(windowed.get_fetches().size() == 2);

    vector<string> names = names_for(windowed, {{105, 105}});
    REQUIRE(names == vector<string>{"c"});
    REQUIRE(windowed.get_fetches().size() == 2);
}

TEST_CASE("A window's read starts outlive its eviction from the cache", "[site_read_source]") {
    // Two cache entries and three windows: window 0 is evicted, but its tally is kept, so
    // counting its starts fetches nothing.
    FakeWindowedSource windowed({{"a", 30}, {"b", 105}, {"c", 250}}, 100, 2);

    names_for(windowed, {{30, 30}});
    names_for(windowed, {{105, 105}});
    names_for(windowed, {{250, 250}});
    REQUIRE(windowed.get_fetches().size() == 3);

    auto counts = starts_on(windowed, {30});
    REQUIRE(counts[60].first == 1);
    REQUIRE(windowed.get_fetches().size() == 3);
}

TEST_CASE("A failed fetch leaves no query waiting for it", "[site_read_source]") {
    // A thread fetching a window claims it, and other threads wait for that fetch. If the fetch
    // fails, the claim must go, or the next query for the window would wait for ever.
    FakeWindowedSource windowed({{"a", 105}}, 100);

    windowed.fail_next_fetch();
    REQUIRE_THROWS(starts_on(windowed, {105}));
    REQUIRE(starts_on(windowed, {105})[60].first == 1);

    windowed.fail_next_fetch();
    REQUIRE_THROWS(names_for(windowed, {{5, 5}}));
    REQUIRE(names_for(windowed, {{5, 5}}).empty());
    REQUIRE(windowed.get_fetches().size() == 4);
}

TEST_CASE("Threads wanting a window at the same time fetch it once", "[site_read_source]") {
    // The cache is shared, so a thread that finds a window being fetched waits for that fetch
    // rather than fetching the window again.
    FakeWindowedSource windowed({{"a", 30}, {"b", 105}}, 100);

    vector<vector<string>> names(64);
    vector<size_t> starts(64);
#pragma omp parallel for num_threads(8)
    for (int i = 0; i < 64; ++i) {
        names[i] = names_for(windowed, {{30, 30}});
        starts[i] = starts_on(windowed, {30, 105})[60].first;
    }

    for (int i = 0; i < 64; ++i) {
        REQUIRE(names[i] == vector<string>{"a"});
        REQUIRE(starts[i] == 2);
    }
    REQUIRE(windowed.get_fetches().size() == 2);
}

}
}
