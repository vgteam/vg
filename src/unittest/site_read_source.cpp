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

#include <string>
#include <vector>

#include "catch.hpp"
#include "site_read_source.hpp"

namespace vg {
namespace unittest {

using namespace std;

/// A read source over reads held in a vector, recording every span fetched.
///
/// Reads are described by the single node they sit on, which is all the windowing logic
/// looks at. Sequence and quality are irrelevant here.
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

    /// Every (min, max) pair fetch_span was called with, in order.
    const vector<pair<nid_t, nid_t>>& get_fetches() const {
        return fetches;
    }

protected:

    void fetch_span(nid_t min_id, nid_t max_id,
                    const function<void(const Alignment&)>& iteratee) const {
        fetches.push_back(make_pair(min_id, max_id));
        for (const Alignment& aln : held) {
            nid_t node = aln.path().mapping(0).position().node_id();
            if (node >= min_id && node <= max_id) {
                if (passes_filter(aln)) {
                    count_fetched();
                    iteratee(aln);
                }
            }
        }
    }

private:
    vector<Alignment> held;
    mutable vector<pair<nid_t, nid_t>> fetches;
};

/// Collect the names of the reads a query returns.
static vector<string> names_for(const SiteReadSource& source,
                                const vector<pair<nid_t, nid_t>>& ranges) {
    vector<string> names;
    source.for_each_read(ranges, [&](const Alignment& aln) {
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
    // get_reads is the vector-returning convenience wrapper; check it agrees with the
    // callback form, since tests and diagnostics use it.
    FakeWindowedSource source({{"a", 10}, {"b", 20}}, 100);

    vector<Alignment> reads = source.get_reads({{10, 20}});

    REQUIRE(reads.size() == 2);
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
        void fetch_span(nid_t min_id, nid_t max_id,
                        const function<void(const Alignment&)>& iteratee) const {
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

}
}
