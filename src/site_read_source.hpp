#ifndef VG_SITE_READ_SOURCE_HPP_INCLUDED
#define VG_SITE_READ_SOURCE_HPP_INCLUDED

/** \file site_read_source.hpp
 *
 * Sources of read alignments by graph locality: given the node IDs of a site, deliver the reads
 * aligned to them. The read-likelihood genotyper uses these to score each read at each site.
 */

#include <atomic>
#include <chrono>
#include <condition_variable>
#include <cstdio>
#include <fstream>
#include <functional>
#include <memory>
#include <mutex>
#include <string>
#include <thread>
#include <unordered_map>
#include <vector>

#include <vg/vg.pb.h>

#include <htslib/hts.h>
#include <htslib/tbx.h>

#include "handle.hpp"
#include "off_panel.hpp"
#include "stream_index.hpp"

namespace vg {

using namespace std;

struct GafRecordScan;

/**
 * The rule for which reads count as evidence.
 *
 * Base qualities filter nothing, since the scoring weighs each base by its own quality.
 */
struct SiteReadFilter {
    /// Drop reads with mapping quality below this.
    int min_mapq = 0;
};

/**
 * A read handed to a site query, with an index of its mappings on the queried nodes.
 *
 * With the index, the work of scoring a read at a site grows with the read's overlap with the
 * site rather than with the read's length, which matters for long reads. A source that does not
 * index leaves `mappings` and `read_offsets` null, and consumers then walk the alignment.
 */
struct SiteRead {
    /// The alignment itself. Valid only for the duration of the callback.
    const Alignment* aln = nullptr;

    /// Ascending indices of the mappings whose node lies in the queried ranges, or null if this
    /// source does not index.
    const uint32_t* mappings = nullptr;
    /// How many indices `mappings` holds.
    size_t mapping_count = 0;

    /// How many read bases come before each mapping, by mapping index, or null if this source
    /// does not index.
    ///
    /// It has one more entry than the path has mappings, so the end of the last mapping can be
    /// read too.
    const uint32_t* read_offsets = nullptr;

    /// Whether the index is present.
    ///
    /// Sources supply both halves of the index or neither, so consumers have one condition to
    /// branch on.
    bool indexed() const { return mappings != nullptr && read_offsets != nullptr; }
};

/**
 * Random-access source of read alignments by graph locality.
 *
 * Implementations must be safe for concurrent read access: GraphCaller visits snarls on several
 * threads at once, so several threads ask for reads at different sites at the same time. Any
 * loading or index building must happen up front, before calling begins.
 *
 * Reads are handed to a callback rather than returned in a vector. Then a backend holding reads in
 * memory need not copy every read for every site, and a backend that decodes reads on demand can
 * hand out a transient reference. Sites are visited many times over a run, so this is the hot path.
 */
class SiteReadSource {
public:
    virtual ~SiteReadSource() = default;

    /// Visit every read with at least one mapping onto a node in any of the given inclusive node ID
    /// ranges.
    ///
    /// Each read is visited at most once, even if it touches several of the ranges. Must be safe
    /// to call concurrently.
    virtual void for_each_read(const vector<pair<nid_t, nid_t>>& ranges,
                               const function<void(const SiteRead&)>& iteratee) const = 0;

    /// Visit every read touching the ranges, as `for_each_read` does, without its index.
    ///
    /// This is named differently rather than overloaded, since a lambda converts to either
    /// `function` type.
    void for_each_alignment(const vector<pair<nid_t, nid_t>>& ranges,
                            const function<void(const Alignment&)>& iteratee) const {
        for_each_read(ranges, [&](const SiteRead& read) { iteratee(*read.aln); });
    }

    /// Visit the reads whose first mapping is onto one of the nodes, in groups.
    ///
    /// `nodes` must be sorted and free of duplicates. Each group shares a mapping quality, and
    /// comes as its number of reads and their total sequence length. Neither the grouping nor the
    /// order is fixed, so a caller can only add the groups up.
    ///
    /// The default fetches the reads touching the nodes and keeps those that begin on one.
    virtual void for_each_read_start(
        const vector<nid_t>& nodes,
        const function<void(int32_t mapq, size_t reads, size_t bases)>& iteratee) const;

    /// Return how many reads this source holds or can see, for logging.
    ///
    /// May be 0 if the backend cannot cheaply say.
    virtual size_t get_read_count() const = 0;

    /// Turn on off-panel detection, counting the reads' substitutions and indels on `graph`.
    ///
    /// Each read is counted once, while it is in memory. `graph` must outlive this. Call before
    /// any query.
    virtual void count_edits(const HandleGraph* graph, const EditCountParams& params) {
        // A source that cannot count edits finds no candidates.
    }

    /// Visit the off-panel candidates on `nodes`, in order.
    ///
    /// `nodes` must be sorted and free of duplicates. A node's candidates come from all the reads
    /// on it: reads not yet counted are fetched and counted first. Nothing is visited unless
    /// `count_edits` was called.
    virtual void for_each_edit_candidate(
        const vector<nid_t>& nodes, const function<void(const EditCandidate&)>& iteratee) const {
        // A source that cannot count edits has no candidates to visit.
    }

    /// Return whether a site query on node `a` also fetches the reads on node `b`.
    ///
    /// The default says yes, which is right where every read is in memory.
    virtual bool fetched_together(nid_t a, nid_t b) const { return true; }

    /// Return every off-panel candidate counted so far, sorted.
    virtual vector<EditCandidate> counted_edit_candidates() const { return {}; }

    /// Start holding reads on this thread for the sites of a subtree, until the matching
    /// `release_reads`.
    ///
    /// Holds nest, and are released innermost first. Within the innermost hold, the first fetch on
    /// this thread for a site too wide for one cached window keeps reads for the sites below it.
    /// Later queries on this thread are answered from the kept reads, with no fetch, when those
    /// include every read on the queried nodes. Read-start queries are answered from them the
    /// same way. With off-panel detection, the hold keeps at least the reads on the nodes of the
    /// candidates `bar` admits, since those nodes become edit sites that are queried later.
    ///
    /// The default holds nothing.
    virtual void hold_reads(const EditCallBar& bar) const {
        // A source that keeps every read in memory has nothing to hold.
    }

    /// End the innermost hold on this thread, dropping its reads.
    virtual void release_reads() const {
        // A source that holds nothing has nothing to release.
    }

};

/**
 * Reads held in memory, indexed by the node IDs they touch.
 *
 * Built in one pass over a GAM or GAF, with no index file. All the reads are held in memory, so a
 * whole-genome read set at ordinary depth will not fit. get_read_count() is logged so that the
 * size is visible.
 */
class InMemorySiteReadSource : public SiteReadSource {
public:

    /// Which reads are eligible.
    using Filter = SiteReadFilter;

    InMemorySiteReadSource() = default;

    /// Stream a GAM, retaining the reads that pass the filter.
    void load_gam(const string& filename, const Filter& filter = Filter());

    /// Stream a GAF, retaining the reads that pass the filter.
    ///
    /// The graph turns GAF records into Alignments.
    void load_gaf(const HandleGraph& graph, const string& filename, const Filter& filter = Filter());

    /// Retain a single read if it passes the filter, as loading would.
    ///
    /// This lets a source be assembled without a file.
    void add(const Alignment& aln, const Filter& filter = Filter());

    /// Visit the reads touching the ranges, without an index of their mappings.
    ///
    /// This source keeps every read for the whole run, and an index would add several bytes per
    /// mapping to that.
    void for_each_read(const vector<pair<nid_t, nid_t>>& ranges,
                       const function<void(const SiteRead&)>& iteratee) const override;

    size_t get_read_count() const override;

    /// Return how many reads the filter rejected, for logging.
    size_t get_filtered_count() const;

    /// Turn on off-panel detection.
    ///
    /// Every read is in memory, so all of them are counted at once, at the first query.
    void count_edits(const HandleGraph* graph, const EditCountParams& params) override;
    void for_each_edit_candidate(const vector<nid_t>& nodes,
                                 const function<void(const EditCandidate&)>& iteratee) const override;
    vector<EditCandidate> counted_edit_candidates() const override;

private:

    /// Retain a read if it passes the filter, indexing it by every node it touches.
    ///
    /// Not safe to call concurrently.
    void add_read(const Alignment& aln, const Filter& filter);

    /// Return the off-panel candidates of all the reads, counting them at the first call.
    const vector<EditCandidate>& all_edits() const;
    /// The graph the reads are counted on, or null when off-panel detection is off.
    const HandleGraph* edit_graph = nullptr;
    /// How the reads are counted.
    EditCountParams edit_params;
    /// Set once `edits` is filled.
    mutable std::once_flag edits_counted;
    /// The off-panel candidates of all the reads.
    mutable vector<EditCandidate> edits;

    /// The reads themselves, owned here and referenced by index below.
    vector<Alignment> reads;

    /// The indices into `reads` of the reads touching each node, by node ID.
    ///
    /// A read appears under every node it touches, so a read spanning n nodes costs n entries.
    unordered_map<nid_t, vector<size_t>> reads_by_node;

    /// Reads rejected by the filter.
    size_t filtered_count = 0;
};

/**
 * Base for on-demand sources, which fetch reads in fixed windows of consecutive node IDs and keep
 * the most recently used windows in a cache that all threads share.
 *
 * A backend query costs much more than one site's reads, because the backend over-fetches or
 * starts a process. So each window is fetched once and serves the sites inside it. That works when
 * sites are visited in node-ID order; see GraphCaller::set_snarl_batching. The cache is shared
 * because a thread also needs windows that other threads are working through, such as those
 * holding the reference nodes near its sites.
 *
 * A site whose nodes lie in more than one window is "wide". It gets its reads by a fetch of its
 * own, from its windows or from the backend directly.
 *
 * With off-panel detection, each window's reads are counted for candidates on the window's own
 * nodes, the first time a query other than a wide site's fetch uses the window. A wide site's fetch counts its nodes that no
 * window has counted. Either way each node is counted once, from every read on it.
 *
 * Subclasses supply fetch_span() and its per-thread resources. The windowing, the cache, and the
 * narrowing of a window's reads to the requested ranges live here, so that all backends interpret
 * a query the same way.
 */
class WindowedSiteReadSource : public SiteReadSource {
public:

    void for_each_read(const vector<pair<nid_t, nid_t>>& ranges,
                       const function<void(const SiteRead&)>& iteratee) const final;

    /// Visit the read starts on the nodes, from tallies of each window's read starts.
    ///
    /// A window is tallied the first time a query other than a wide site's fetch uses it, and the
    /// tallies are kept for the whole run. So a query fetches only windows that were never
    /// tallied.
    void for_each_read_start(
        const vector<nid_t>& nodes,
        const function<void(int32_t mapq, size_t reads, size_t bases)>& iteratee) const final;

    /// Return how many reads have been fetched from the backend so far, across all threads.
    ///
    /// An on-demand backend never knows the size of the read set.
    size_t get_read_count() const override;

    /// Turn on off-panel detection.
    ///
    /// The candidates are kept for the whole run.
    void count_edits(const HandleGraph* graph, const EditCountParams& params) final;
    void for_each_edit_candidate(const vector<nid_t>& nodes,
                                 const function<void(const EditCandidate&)>& iteratee) const final;
    /// Return whether the two nodes lie in one window.
    bool fetched_together(nid_t a, nid_t b) const final;
    vector<EditCandidate> counted_edit_candidates() const final;
    void hold_reads(const EditCallBar& bar) const final;
    void release_reads() const final;

    /**
     * A tally of the fetches and queries of a windowed source, for logging.
     *
     * The source keeps one with atomic counts, which every thread adds to, and hands out a
     * snapshot as a `FetchStats`.
     */
    template<typename Count>
    struct FetchTally {
        /// Reads fetched from the backend.
        Count fetched{0};
        /// Reads rejected by the filter.
        Count filtered{0};
        /// Site queries inside one window that found it cached.
        Count cache_hits{0};
        /// Site queries that fetched from the backend.
        Count cache_misses{0};
        /// Index entries examined for site queries.
        ///
        /// A read can appear under several of the nodes a site asks about, so the ratio of this to
        /// `delivered` says how much of the index a site reads again.
        Count scanned{0};
        /// Reads handed to site queries.
        Count delivered{0};
        /// Wide site queries fetched from the backend by their own ranges.
        Count straddles{0};
        /// Wide site queries answered from their windows.
        Count wide_from_windows{0};
        /// The total node-ID span, from lowest to highest, of the wide sites fetched by their own
        /// ranges.
        Count straddle_nodes{0};
        /// The total number of node IDs those sites asked for.
        Count straddle_wanted{0};
        /// Windows fetched whole.
        Count whole_fetches{0};
        /// Whole-window fetches of a window tallied before, which the cache had dropped.
        Count refetches{0};
        /// Whole-window fetches made only to count a window's off-panel candidates.
        Count edit_fetches{0};
        /// Site and read-start queries answered from held reads.
        Count held_queries{0};
        /// Reads in the windows fetched whole.
        Count whole_fetch_reads{0};
    };
    /// A snapshot of the tally of a source's fetches and queries.
    using FetchStats = FetchTally<size_t>;

    /// Return a snapshot of the tally of the fetches and queries so far, across all threads.
    FetchStats fetch_stats() const;

    /// Drop every cached window and return how many reads they held.
    ///
    /// Call this when calling is done. Otherwise the cache keeps its windows through the passes
    /// that follow calling, when the run uses the most memory. A later query fetches its window
    /// again.
    size_t drop_cached_windows();

    /// Sort the ranges and merge overlapping and adjacent ones.
    ///
    /// The result covers the same node IDs, and can be searched by `in_ranges` and `touches`.
    static vector<pair<nid_t, nid_t>> merge_ranges(vector<pair<nid_t, nid_t>> ranges);

    /// Return whether the node ID is in one of the ranges.
    ///
    /// @pre `ranges` comes from merge_ranges.
    static bool in_ranges(nid_t node_id, const vector<pair<nid_t, nid_t>>& ranges);

    /// Return whether the read has a mapping onto a node in the ranges.
    ///
    /// @pre `ranges` comes from merge_ranges.
    static bool touches(const Alignment& aln, const vector<pair<nid_t, nid_t>>& ranges);


protected:

    /// Make a source with windows of `window_nodes` node IDs, caching at most `max_cached_windows`
    /// windows in all.
    WindowedSiteReadSource(const SiteReadFilter& filter, size_t window_nodes,
                           size_t max_cached_windows);

    /// Fetch as `fetch_span` does, also giving each read a key, and return how many reads the
    /// fetch dropped for having the name and start of an earlier read.
    ///
    /// A read's key is its position in the backend's file. It is the same for a read whichever
    /// fetch returns it, so reads from several fetches can be put back in file order, once each.
    /// For a backend whose `has_record_keys` is false, the keys are all 0 and no dropped reads are
    /// reported.
    virtual size_t fetch_span_keyed(const vector<pair<nid_t, nid_t>>& ranges,
                                    const function<void(Alignment&, uint64_t)>& iteratee) const {
        fetch_span(ranges, [&](Alignment& aln) { iteratee(aln, 0); });
        return 0;
    }
    /// Return whether `fetch_span_keyed` gives reads real keys.
    virtual bool has_record_keys() const { return false; }

    /// Visit every read with a mapping onto a node in any of the inclusive ranges, having applied
    /// the filter.
    ///
    /// Each read must be visited at most once. Must be safe to call concurrently, so
    /// implementations own their per-thread resources. This takes ranges rather than one span,
    /// because a snarl's nodes can be spread thinly over a wide span of IDs.
    ///
    /// The iteratee takes a mutable reference, so the caller may move from the alignment. It is
    /// the backend's per-record scratch.
    virtual void fetch_span(const vector<pair<nid_t, nid_t>>& ranges,
                            const function<void(Alignment&)>& iteratee) const = 0;

    /// Apply the filter to a read, counting a rejection.
    ///
    /// Subclasses call this on each candidate read before handing it to fetch_span's iteratee.
    bool passes_filter(const Alignment& aln) const;

    /// Count a read as fetched.
    ///
    /// This is separate from passes_filter so a subclass can decide the order in which it filters
    /// and counts.
    void count_fetched() const;

private:

    /// One mapping of one read in a window, as the window's node index holds it.
    struct IndexEntry {
        nid_t node = 0;
        uint32_t read = 0;
        uint32_t mapping = 0;
        /// Order by node, then read, then mapping.
        auto operator<=>(const IndexEntry& other) const = default;
    };

    /// The reads of one fetched window, with an index of their mappings.
    ///
    /// It is not changed once fetched, so threads read it without a lock.
    struct CacheEntry {
        size_t window = 0;
        /// The reads touching the window's nodes, in the order the backend returned them.
        vector<Alignment> reads;
        /// Each read's key (see `fetch_span_keyed`), in step with `reads`.
        vector<uint64_t> keys;
        /// How many reads the window's fetch dropped for having the name and start of an earlier
        /// read.
        size_t dropped_repeats = 0;

        /// Every mapping in the window, sorted by node ID.
        ///
        /// A site finds the reads it wants, and which of their mappings, by binary search. Built
        /// once per fetch and used by every site in the window.
        vector<IndexEntry> node_index;

        /// The read offset before each mapping, for every read in the window end to end.
        ///
        /// Read r occupies `[offset_start[r], offset_start[r + 1])`. That is one entry longer than
        /// its path, so the end of the last mapping is readable too.
        vector<uint32_t> offsets;
        /// Where each read's run of `offsets` starts, and one more entry for the end of the last.
        vector<uint32_t> offset_start;
    };

    /**
     * A table of the fetched windows the source keeps, which drops the least recently used
     * beyond its capacity.
     *
     * A thread claims a window before fetching it, so that other threads wanting the window wait
     * for that fetch rather than fetching it again. A claimed window has a slot but no entry.
     * Call every method holding the source's `cache_mutex`.
     */
    class WindowCache {
    public:
        /// Make a table holding at most `capacity` fetched windows.
        explicit WindowCache(size_t capacity);

        /// Return the window's entry, marking it used, or claim the window and return null.
        ///
        /// While another thread has the window claimed, this waits, releasing `lock`.
        shared_ptr<const CacheEntry> find_or_claim(size_t window, unique_lock<std::mutex>& lock);

        /// Give up a claim on a window whose fetch failed, so that a thread waiting for the window
        /// fetches it itself.
        void release_claim(size_t window);

        /// Store the entry of a window this thread claimed, and drop the least recently used
        /// windows beyond the capacity.
        ///
        /// The dropped entries are moved into `dropped`, for the caller to free after releasing
        /// the lock. A thread still reading a dropped window keeps it alive through its own
        /// pointer.
        void store(size_t window, shared_ptr<const CacheEntry> entry,
                   vector<shared_ptr<const CacheEntry>>& dropped);

        /// Empty the table, returning the entries it held.
        vector<shared_ptr<const CacheEntry>> take_all();

    private:
        /// A window in the table. It has no entry while a thread is fetching it.
        struct Slot {
            shared_ptr<const CacheEntry> entry;
            /// When a query last used the window, as a value of `clock`.
            size_t last_use = 0;
        };

        /// The most fetched windows the table holds.
        size_t capacity;
        /// The windows in the table, by window index.
        unordered_map<size_t, Slot> slots;
        /// The number of uses so far, which orders them for eviction.
        size_t clock = 0;
        /// Signalled when a claimed window is stored or released.
        std::condition_variable filled;
    };

    /// Return the window's reads, from the cache or fetched now.
    ///
    /// Sets `was_fetched` if this call did the fetch. The entry stays valid while the caller holds
    /// it, even if the cache drops it meanwhile.
    ///
    /// Unless `counted` is false, this tallies the window's read starts and counts its off-panel
    /// candidates, if that has not been done. A wide site's fetch passes false, which leaves its
    /// windows as if it had never fetched them.
    shared_ptr<const CacheEntry> get_window(size_t window, bool& was_fetched,
                                            bool counted = true) const;

    /// Fetch a window this thread claimed, releasing the claim if the fetch throws.
    shared_ptr<CacheEntry> load_claimed_window(size_t window) const;

    /// Return whether the window's read starts are tallied and, with off-panel detection, its
    /// candidates counted.
    ///
    /// Call holding `cache_mutex`.
    bool summarized(size_t window) const;

    /// Hand the entry's reads that touch the ranges to the caller, in the order they were fetched.
    ///
    /// Reads are found through the entry's node index.
    void deliver(const CacheEntry& entry,
                 const vector<pair<nid_t, nid_t>>& ranges,
                 const function<void(const SiteRead&)>& iteratee) const;

    /// Hand a wide site's reads to the caller, from a fetch of its own.
    ///
    /// The reads come in file order, from the site's windows if there are few of them, or else
    /// from the backend by the site's ranges. This fills the innermost hold on this thread, if it
    /// is unfilled, and counts the site's nodes that no fetch has counted.
    void deliver_wide(const vector<pair<nid_t, nid_t>>& ranges, nid_t min_id, nid_t max_id,
                      const function<void(const SiteRead&)>& iteratee) const;

    /// Append this read's mappings to a window's index and its read-offset table.
    static void index_read(const Alignment& aln, uint32_t read_index, CacheEntry& entry);

    /// Return which window a node ID falls in.
    size_t window_of(nid_t id) const;

    /// Fetch a whole window from the backend and index it.
    CacheEntry load_window(size_t window) const;

    /// Find the windows holding the merged ranges' nodes, ascending, and return whether a query
    /// over them can be answered from them.
    ///
    /// That needs at most `WIDE_QUERY_MAX_WINDOWS` windows, and a backend that keys its reads.
    bool few_windows(const vector<pair<nid_t, nid_t>>& merged, vector<size_t>& windows) const;

    /// Visit the reads touching the merged ranges from the windows holding their nodes, as a fetch
    /// of the ranges from the backend would deliver them.
    ///
    /// The windows are fetched whole and cached like any other, with `counted` as for
    /// `get_window`. The reads come in file order, each once. A read with the name and start of an
    /// earlier one is dropped and counted in `dropped`.
    ///
    /// @return false, having visited nothing, if one of the windows' own fetches dropped a read.
    /// That window can then lack a read that the fetch by ranges delivers.
    bool for_each_window_read(const vector<pair<nid_t, nid_t>>& merged,
                              const vector<size_t>& windows, bool counted,
                              const function<void(const Alignment&)>& iteratee,
                              size_t& dropped) const;

    /// The reads beginning on one node with one mapping quality.
    struct StartTally {
        nid_t node = 0;
        int32_t mapq = 0;
        /// How many reads.
        size_t reads = 0;
        /// Their total sequence length.
        size_t bases = 0;
    };

    /// Tally the read starts of a fetched window, sorted by node, then mapping quality.
    ///
    /// A read is tallied under the window holding its first node. A window's fetch also returns
    /// reads that begin in other windows, and those are left for their own windows.
    vector<StartTally> tally_starts(const CacheEntry& entry) const;

    /// Return the window's read-start tallies, fetching the window if it was never tallied.
    ///
    /// The tallies are kept for the whole run and never changed, so the reference stays valid.
    const vector<StartTally>& window_starts(size_t window) const;

    /// Count the off-panel candidates on the window's own nodes, in its reads.
    ///
    /// Every read on one of the window's nodes is among its reads.
    vector<EditCandidate> count_window_edits(const CacheEntry& entry) const;

    /// Return the candidates of the window, fetching it if it was never counted.
    const vector<EditCandidate>& window_edits(size_t window) const;

    /// Answer the query from the innermost of this thread's filled holds that can answer it (see
    /// `hold_reads`), and return whether there was one.
    ///
    /// `wide` says whether the query is too wide for one window.
    bool deliver_held(const vector<pair<nid_t, nid_t>>& ranges, bool wide,
                      const function<void(const SiteRead&)>& iteratee) const;

    /// Return, as merged ranges, the nodes in the merged ranges that are on a candidate `bar`
    /// admits, among those counted so far.
    ///
    /// Call holding `cache_mutex`.
    vector<pair<nid_t, nid_t>> candidate_nodes(const vector<pair<nid_t, nid_t>>& merged,
                                               const EditCallBar& bar) const;

    /// Return, as merged ranges, the nodes in the merged ranges that no fetch has counted.
    ///
    /// Call holding `cache_mutex`.
    vector<pair<nid_t, nid_t>> uncounted(const vector<pair<nid_t, nid_t>>& merged) const;

    /// Record the candidates a wide site's fetch counted on the nodes `counted_nodes`.
    void record_span_edits(vector<pair<nid_t, nid_t>> counted_nodes,
                           const vector<EditCandidate>& counted) const;

    /// The rule for which fetched reads are kept.
    SiteReadFilter filter;

    /// How many node IDs each window holds.
    size_t window_size;

    /// The windows fetched and kept.
    mutable WindowCache cache;

    /// The tally of the fetches and queries, which `fetch_stats` reports.
    mutable FetchTally<atomic<size_t>> counters;

    /// Read-start tallies by window, added when a window is first tallied and kept for the whole
    /// run.
    mutable unordered_map<size_t, vector<StartTally>> starts;

    /// The graph the reads are counted on, or null when off-panel detection is off.
    const HandleGraph* edit_graph = nullptr;
    /// How the reads are counted.
    EditCountParams edit_params;
    /// The off-panel candidates by window, added when a window is first counted and kept for the
    /// whole run.
    mutable unordered_map<size_t, vector<EditCandidate>> edits;
    /// The nodes counted by wide sites' fetches, as merged ranges.
    mutable vector<pair<nid_t, nid_t>> span_counted;
    /// The candidates on `span_counted`, sorted.
    mutable vector<EditCandidate> span_edits;

    /// Guards `cache`, `starts`, `edits`, `span_counted` and `span_edits`.
    ///
    /// Fetches run without it.
    mutable std::mutex cache_mutex;
};

/**
 * A locker for each thread, holding the resources a backend gives each thread, such as an open
 * file that seeks.
 *
 * A thread gets its own `State` on its first call to `mine`, and keeps it for the life of the
 * lockers. Lockers are keyed by system thread, not OpenMP thread number. Threads in different
 * nested OpenMP teams can share a thread number, but each still needs a locker of its own.
 */
template<typename State>
class PerThreadStates {
public:

    /// Return the calling thread's state, made at the thread's first call.
    ///
    /// The reference stays valid until the lockers are destroyed, while other threads add theirs.
    State& mine() {
        std::scoped_lock guard(mutex);
        unique_ptr<State>& state = by_thread[std::this_thread::get_id()];
        if (!state) {
            state = make_unique<State>();
        }
        return *state;
    }

    /// Visit every thread's state.
    ///
    /// Call only when no thread is using its state, as when tearing down.
    template<typename Visit>
    void for_each(const Visit& visit) const {
        for (const auto& entry : by_thread) {
            visit(*entry.second);
        }
    }

private:

    /// Guards `by_thread`.
    std::mutex mutex;
    /// Each thread's state, by thread.
    ///
    /// The states are held by pointer, so adding a thread never moves one another thread is using.
    unordered_map<std::thread::id, unique_ptr<State>> by_thread;
};

/**
 * Reads fetched on demand from a sorted GAM and its `.gai` index (`vg gamsort -i`), so memory is
 * bounded by what one window needs.
 *
 * Two properties of StreamIndex shape the implementation:
 *
 * * One cursor per thread. Concurrent `find()` calls are safe, but a cursor seeks, so it cannot
 *   be shared. Cursors are created per thread on first use, as `vg chunk` does.
 * * It over-fetches. The index points to the first group of reads that may touch a node, not to
 *   the reads themselves, so a query reads groups in file order until their smallest node ID
 *   passes the query's. The windows of the base class limit the cost.
 */
class IndexedGamSiteReadSource : public WindowedSiteReadSource {
public:

    /// Make a source over a GAM sorted by `vg gamsort -i`, and its index.
    ///
    /// window_size is in node IDs, and cache_entries is how many windows the cache shared by all
    /// threads holds.
    IndexedGamSiteReadSource(const string& gam_filename, const string& index_filename,
                             const SiteReadFilter& filter = SiteReadFilter(),
                             size_t window_size = 256,
                             size_t cache_entries = 2);

protected:

    void fetch_span(const vector<pair<nid_t, nid_t>>& ranges,
                    const function<void(Alignment&)>& iteratee) const override;

private:

    /// One thread's open GAM and the index cursor reading it.
    struct ThreadState {
        unique_ptr<ifstream> stream;
        unique_ptr<GAMIndex::cursor_t> cursor;
    };

    /// Return this thread's state, opening its cursor on first use.
    ThreadState& thread_state() const;

    string gam_filename;
    unique_ptr<GAMIndex> index;

    /// Each thread's state.
    mutable PerThreadStates<ThreadState> threads;
};

/**
 * Reads fetched on demand from a GAF-Base database, by running `gbz-base query`.
 *
 * <https://github.com/jltsiren/gbz-base> stores alignments in SQLite and returns the reads that
 * overlap a set of nodes, without over-fetching. GAF-Base has no C API, and its file format may
 * change behind a version check, so we run its own binary rather than decode the format here.
 * This needs `gbz-base` at run time only when this source is used, and nothing at build time. The
 * binary writes GAF text, which is parsed, filtered, windowed and cached like the other sources.
 *
 * Starting a process is slow, so this fetches one window of node IDs per process.
 */
class GafBaseSiteReadSource : public WindowedSiteReadSource {
public:

    /// Make a source over a GAF-Base database.
    ///
    /// graph must outlive this, and must be the graph the alignments were made against. It
    /// supplies node lengths and sequences to turn GAF back into Alignments, and says which node
    /// IDs in a window exist.
    ///
    /// gbz_filename is a GBZ or a GBZ-Base (`gbz-base construct`). Prefer the latter: a plain GBZ
    /// is loaded in full on every query, while a GBZ-Base is random-access, and this issues many
    /// queries.
    ///
    /// window_size is in node IDs, and cache_entries is how many windows the cache shared by all
    /// threads holds.
    GafBaseSiteReadSource(const HandleGraph& graph,
                          const string& gaf_base_filename,
                          const string& gbz_filename,
                          const SiteReadFilter& filter = SiteReadFilter(),
                          size_t window_size = 256,
                          size_t cache_entries = 2,
                          const string& binary = "gbz-base");

    ~GafBaseSiteReadSource() override;

    /// Return how many gbz-base processes were started.
    ///
    /// This is the number that matters for run time, since each one costs milliseconds no matter
    /// how few reads it returns.
    size_t get_query_count() const;

    /// Return how many reads a fetch returned more than once and dropped.
    ///
    /// A split query returns a read that spans two of its pieces from both.
    size_t get_duplicate_count() const { return duplicates_dropped.load(); }

    /// Run one query up front, to check that the databases are readable and agree with the graph.
    ///
    /// A broken setup then fails immediately rather than on the first snarl in a worker thread.
    /// Throws with an actionable message on failure.
    void check_setup() const;

    /// Return whether gbz-base is handed the databases as SQLite immutable URIs (see `immutable`).
    bool immutable_databases() const { return immutable; }

    /// Where the time of the gbz-base queries went, summed over every query so far.
    ///
    /// It is collected with wait4() and two clock reads per query, so it costs nothing
    /// measurable.
    struct QueryTotals {
        size_t queries = 0;
        /// How long the calling threads waited for their children.
        double wait_s = 0;
        /// The CPU time the children spent in their own code.
        double child_user_s = 0;
        /// The CPU time the children spent in the kernel.
        double child_system_s = 0;
        /// The time spent parsing the GAF the children wrote.
        double parse_s = 0;
        /// The bytes of GAF the children wrote.
        size_t gaf_bytes = 0;
    };
    /// Return where the time of the queries went so far.
    QueryTotals get_query_totals() const;

private:

    /// One thread's GAF output files, one per query the thread can have running.
    ///
    /// They are reused across queries rather than created per query, since creating one is a
    /// filesystem round trip.
    struct ThreadState {
        vector<string> gaf_paths;
    };

    /// Return this thread's state.
    ThreadState& thread_state() const;

    /// Return this thread's output file for the given query slot, created on first use.
    const string& gaf_path(ThreadState& state, size_t slot) const;

    /// One `gbz-base query` child, started and not yet reaped.
    struct PendingQuery {
        pid_t pid = 0;
        /// The child's GAF output file.
        ///
        /// It is held by value, since starting a query in a later slot can grow the thread's
        /// vector of paths and move them.
        string gaf_path;
        /// The file the child's stderr goes to.
        string err_path;
        /// When the query started, for its wall time in `QueryTotals` and the query log.
        std::chrono::steady_clock::time_point started;
        /// When the query started, in seconds since the epoch, for the query log.
        double started_epoch = 0;
        /// How many nodes the query asks for.
        size_t n_nodes = 0;
        /// The first node the query asks for.
        nid_t first_node = 0;
        /// The last node the query asks for.
        nid_t last_node = 0;
        /// FNV-1a over the node IDs, so the same node set fetched twice shows up in the log.
        uint64_t nodes_hash = 0;
    };

    void fetch_span(const vector<pair<nid_t, nid_t>>& ranges,
                    const function<void(Alignment&)>& iteratee) const override;

    /// Start `gbz-base query` for these node IDs, writing to this thread's output file for `slot`.
    ///
    /// This does not wait for the child. The caller reaps it with reap_query, so a query split in
    /// pieces runs its pieces at the same time. Throws if the child cannot be started.
    PendingQuery spawn_query(ThreadState& state, size_t slot,
                             const vector<nid_t>& nodes) const;

    /// Wait for a started query and parse the GAF it wrote, applying the filter.
    ///
    /// Returns the number of records parsed. Throws if the child failed.
    size_t reap_query(const PendingQuery& pending,
                      const function<void(Alignment&)>& iteratee) const;

    /// Start a query and reap it. Throws on failure.
    size_t run_query(ThreadState& state, const vector<nid_t>& nodes,
                     const function<void(Alignment&)>& iteratee) const;

    /// Run a query, reporting an error and exiting instead of throwing.
    ///
    /// Calling happens inside an OpenMP parallel region, which an exception must not leave.
    void run_query_or_die(ThreadState& state, const vector<nid_t>& nodes,
                          const function<void(Alignment&)>& iteratee) const;

    const HandleGraph& graph;
    string gaf_base_filename;
    string gbz_filename;
    /// The gbz-base executable, as a path or a name to look up on PATH.
    string binary;

    /// Whether the databases are handed to gbz-base as SQLite immutable URIs.
    ///
    /// SQLite takes a POSIX advisory lock on the database file around every read transaction, and
    /// gbz-base looks up each node of a query in a transaction of its own. With one gbz-base per
    /// calling thread, all of those lock and unlock calls land on the same file. On a whole genome
    /// the kernel's bookkeeping for them took most of every query's time. vg only reads the
    /// databases, so by default it hands them over as immutable URIs, which take no locks.
    /// Setting VG_GAFBASE_LOCKING=1 in the environment passes the plain paths instead.
    bool immutable = true;
    /// The graph database argument gbz-base is given: the plain path, or `file:<path>?immutable=1`.
    string gbz_argument;
    /// The GAF-Base database argument gbz-base is given, in the same form as `gbz_argument`.
    string gaf_base_argument;
    /// Return `file:<absolute path>?immutable=1` for a SQLite database, or the path unchanged for
    /// any other file, such as a plain GBZ given as the graph.
    static string immutable_uri(const string& path);

    /// The parts of QueryTotals, in microseconds or bytes.
    mutable atomic<uint64_t> wait_us{0};
    mutable atomic<uint64_t> child_user_us{0};
    mutable atomic<uint64_t> child_system_us{0};
    mutable atomic<uint64_t> parse_us{0};
    mutable atomic<uint64_t> gaf_bytes_total{0};
    /// The query log, when VG_GAFBASE_QUERY_LOG names a file, or null.
    ///
    /// It has one line per query. The columns are named in its header line.
    FILE* query_log = nullptr;
    /// Guards `query_log`.
    mutable std::mutex query_log_mutex;

    /// Return the most node IDs one gbz-base command line can hold.
    ///
    /// That is a quarter of sysconf(_SC_ARG_MAX), leaving room for the environment and the fixed
    /// arguments, over the bytes each node costs ("-n" and its digits). It is computed once.
    static size_t argv_node_budget();

    /// The most node IDs per gbz-base process.
    ///
    /// A larger query is split into pieces run side by side.
    size_t max_query_nodes = argv_node_budget();

    /// Each thread's state.
    mutable PerThreadStates<ThreadState> threads;
    /// How many gbz-base processes were started.
    mutable atomic<size_t> queries{0};
    /// Reads a fetch returned more than once, identified by name and start, and dropped.
    mutable atomic<size_t> duplicates_dropped{0};
};

/**
 * Reads fetched on demand from a bgzipped GAF sorted by node ID, through its tabix index.
 *
 * The index comes from `tabix -p gaf`, which indexes each record by the smallest and largest node
 * ID on its path. The reads are read in vg's own threads. A query seeks to the compressed blocks
 * holding the records whose node interval overlaps it, so only those blocks are decompressed and
 * only those records parsed. A whole window costs one index lookup rather than a gbz-base process,
 * and a read is decoded only for the windows it touches. Reads come back in file order, as
 * GAF-Base returns them from the same sorted GAF.
 */
class TabixGafSiteReadSource : public WindowedSiteReadSource {
public:

    /// Make a source over a sorted, bgzipped GAF and its tabix index.
    ///
    /// graph must outlive this, and must be the graph the alignments were made against. It
    /// supplies the node sequences that turn GAF back into Alignments. index_filename is the
    /// `.tbi` or `.csi` beside gaf_filename. window_size is in node IDs, and cache_entries is how
    /// many windows the cache shared by all threads holds.
    TabixGafSiteReadSource(const HandleGraph& graph,
                           const string& gaf_filename,
                           const string& index_filename,
                           const SiteReadFilter& filter = SiteReadFilter(),
                           size_t window_size = 16384,
                           size_t cache_entries = 2);

    ~TabixGafSiteReadSource() override;

    /// Return how many index lookups were made.
    ///
    /// A fetch makes one per group of nearby node IDs it asks for.
    size_t get_query_count() const { return queries.load(); }

    /// Return how many records were dropped because their path visits none of a fetch's nodes.
    ///
    /// The index returns such a record when its node interval overlaps the fetch's ranges. Most
    /// of them are dropped before being parsed (see get_unparsed_count).
    size_t get_skipped_count() const { return skipped.load(); }

    /// Return how many reads a fetch returned more than once and dropped.
    ///
    /// Reads are identified by name and start.
    size_t get_duplicate_count() const { return duplicates_dropped.load(); }

    /// Return how many of the skipped records were dropped from their path text alone, before
    /// being parsed.
    size_t get_unparsed_count() const { return unparsed.load(); }

    /// Return the time all threads spent reading records through the index, seeking and
    /// decompressing.
    double get_read_seconds() const { return read_us.load() / 1e6; }
    /// Return the time all threads spent parsing records into Alignments.
    double get_parse_seconds() const { return parse_us.load() / 1e6; }
    /// Return the bytes of GAF read.
    size_t get_gaf_bytes() const { return gaf_bytes.load(); }

protected:

    void fetch_span(const vector<pair<nid_t, nid_t>>& ranges,
                    const function<void(Alignment&)>& iteratee) const override;

    /// Fetch as `fetch_span` does, keying each read by its record's virtual file offset after it
    /// was read.
    ///
    /// That offset is the same whichever lookup reads the record, and increases in file order.
    size_t fetch_span_keyed(const vector<pair<nid_t, nid_t>>& ranges,
                            const function<void(Alignment&, uint64_t)>& iteratee) const override;
    bool has_record_keys() const override { return true; }

private:

    /// One thread's handle on the GAF.
    ///
    /// A handle seeks, so it cannot be shared. The index can.
    struct ThreadState {
        htsFile* file = nullptr;
    };

    /// Return this thread's state, opening its handle on first use.
    ThreadState& thread_state() const;

    /// Read the records the index returns for one run of node IDs, handing each to `visit` with
    /// the file offset the reader reached after it.
    ///
    /// `scan` notes, for each record, whether it may touch the fetch's ranges. Adds the time spent
    /// reading to `read_time_us`. Throws if the index cannot be queried or the file cannot be read.
    void for_each_record_in_run(htsFile* file, const pair<nid_t, nid_t>& run, GafRecordScan& scan,
                                kstring_t& line, uint64_t& read_time_us,
                                const function<void(const kstring_t&, uint64_t)>& visit) const;

    const HandleGraph& graph;
    string gaf_filename;
    tbx_t* index = nullptr;

    /// Each thread's state.
    mutable PerThreadStates<ThreadState> threads;
    /// The counts and times the getters report.
    mutable atomic<size_t> queries{0};
    mutable atomic<size_t> skipped{0};
    mutable atomic<size_t> unparsed{0};
    mutable atomic<size_t> duplicates_dropped{0};
    mutable atomic<uint64_t> read_us{0};
    mutable atomic<uint64_t> parse_us{0};
    mutable atomic<size_t> gaf_bytes{0};
};

}

#endif
