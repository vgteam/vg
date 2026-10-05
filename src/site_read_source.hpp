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
#include <unordered_map>
#include <vector>

#include <vg/vg.pb.h>

#include "handle.hpp"
#include "stream_index.hpp"

namespace vg {

using namespace std;

/**
 * Which reads are used as evidence: those whose mapping quality reaches a minimum. Base qualities
 * filter nothing, since the scoring weighs each base by its own quality.
 */
struct SiteReadFilter {
    /// Drop reads with mapping quality below this (--read-min-mapq).
    int min_mapq = 0;
};

/**
 * A read handed to a site query, with an index of its mappings in the queried ranges, so that a
 * site's work grows with the read's overlap with the site rather than with the read's length,
 * which matters for long reads.
 *
 * `mappings` and `read_offsets` may be null for a source that does not index, and consumers then
 * walk the alignment instead.
 */
struct SiteRead {
    /// The alignment itself. Valid only for the duration of the callback.
    const Alignment* aln = nullptr;

    /// Ascending indices of the mappings whose node lies in the queried ranges, or
    /// null if this source does not index.
    const uint32_t* mappings = nullptr;
    size_t mapping_count = 0;

    /// How many read bases come before each mapping, indexed by mapping index and one
    /// longer than the path, or null if this source does not index.
    const uint32_t* read_offsets = nullptr;

    /// Whether the index is present. Both halves are supplied together or not at all,
    /// so consumers have one condition to branch on.
    bool indexed() const { return mappings != nullptr && read_offsets != nullptr; }
};

/**
 * Random-access source of read alignments by graph locality.
 *
 * Implementations must be safe for concurrent read access: GraphCaller visits
 * snarls on several threads at once, so several threads will be asking for
 * reads at different sites at the same time. Any loading or index building
 * must happen up front, before calling begins.
 *
 * Reads are handed to a callback rather than returned in a vector so that
 * backends holding reads in memory do not have to copy every read for every
 * site, while backends that decode reads on demand can hand out a transient
 * reference. Sites are visited many times over a run, so this is the hot path.
 */
class SiteReadSource {
public:
    virtual ~SiteReadSource() = default;

    /// Visit every read with at least one mapping onto a node in any of the
    /// given inclusive node ID ranges. Each read is visited at most once, even
    /// if it touches several of the ranges. Must be safe to call concurrently.
    virtual void for_each_read(const vector<pair<nid_t, nid_t>>& ranges,
                               const function<void(const SiteRead&)>& iteratee) const = 0;

    /// The same, for consumers with no use for the index. Named differently rather than
    /// overloaded, since a lambda converts to either `function` type.
    void for_each_alignment(const vector<pair<nid_t, nid_t>>& ranges,
                            const function<void(const Alignment&)>& iteratee) const {
        for_each_read(ranges, [&](const SiteRead& read) { iteratee(*read.aln); });
    }


    /// Visit the reads that begin on one of the nodes, meaning their first mapping is onto it.
    /// `nodes` must be sorted and free of duplicates. Reads come in groups sharing a mapping
    /// quality, as the group's size and its total sequence length. Neither the grouping nor the
    /// order is fixed, so a caller must only add the groups up.
    ///
    /// The default fetches the reads touching the nodes and keeps those that begin on one.
    /// WindowedSiteReadSource answers from per-window tallies instead.
    virtual void for_each_read_start(
        const vector<nid_t>& nodes,
        const function<void(int32_t mapq, size_t reads, size_t bases)>& iteratee) const;

    /// How many reads this source holds or can see, for logging. May be 0 if
    /// the backend cannot cheaply say.
    virtual size_t get_read_count() const = 0;

};

/**
 * Reads held in memory, indexed by the node IDs they touch.
 *
 * Built in one pass over a GAM or GAF, with no index file. All the reads are held in memory, so
 * a whole-genome read set at ordinary depth will not fit; get_read_count() is logged so that the
 * size is visible.
 */
class InMemorySiteReadSource : public SiteReadSource {
public:

    /// Which reads are eligible. Defined at namespace scope as SiteReadFilter;
    /// aliased here for readability at call sites.
    using Filter = SiteReadFilter;

    InMemorySiteReadSource() = default;

    /// Stream a GAM, retaining the reads that pass the filter.
    void load_gam(const string& filename, const Filter& filter = Filter());

    /// Stream a GAF. Needs the graph to turn GAF into Alignments.
    void load_gaf(const HandleGraph& graph, const string& filename, const Filter& filter = Filter());

    /// Retain a single read directly, applying the same filter as loading would, so that a
    /// source can be assembled without a file, as the unit tests do.
    void add(const Alignment& aln, const Filter& filter = Filter());

    /// Reads are delivered without an index: this source keeps every read for the whole run,
    /// and an index would add several bytes per mapping to that.
    void for_each_read(const vector<pair<nid_t, nid_t>>& ranges,
                       const function<void(const SiteRead&)>& iteratee) const;

    size_t get_read_count() const;

    /// How many reads the filter rejected, for logging.
    size_t get_filtered_count() const;

private:

    /// Retain a read if it passes the filter, indexing it by every node it
    /// touches. Not safe to call concurrently; loading is single-threaded.
    void add_read(const Alignment& aln, const Filter& filter);

    /// The reads themselves, owned here and referenced by index below.
    vector<Alignment> reads;

    /// Node ID to indices into `reads`. A read appears under every node it
    /// touches, so a read spanning n nodes costs n entries.
    unordered_map<nid_t, vector<size_t>> reads_by_node;

    /// Reads rejected by the filter.
    size_t filtered_count = 0;
};

/**
 * Base for on-demand sources: rounds each fetch out to fixed windows of consecutive node IDs and
 * keeps the most recently used windows in a cache that all threads share.
 *
 * A backend query costs much more than one site's reads, because the backend over-fetches or
 * starts a process, so each window is fetched once and serves the sites inside it. That works
 * when sites are visited in node-ID order; see GraphCaller::set_node_id_ordering. The cache is
 * shared because a thread also needs windows that other threads are working through, such as
 * those holding the reference nodes near its sites.
 *
 * Subclasses supply fetch_span() and its per-thread resources. The window arithmetic, the
 * cache, the handling of queries that cross a window boundary, and the narrowing of a window to
 * the requested ranges live here, so that all backends interpret a query the same way.
 */
class WindowedSiteReadSource : public SiteReadSource {
public:

    void for_each_read(const vector<pair<nid_t, nid_t>>& ranges,
                       const function<void(const SiteRead&)>& iteratee) const final;

    /// Answered from tallies of each window's read starts. A window is tallied the first time
    /// it is fetched whole, and the tallies are kept for the whole run, so a query fetches only
    /// windows that have never been fetched whole.
    void for_each_read_start(
        const vector<nid_t>& nodes,
        const function<void(int32_t mapq, size_t reads, size_t bases)>& iteratee) const final;

    /// Reads actually fetched from the backend so far, across all threads. Not the
    /// size of the read set, which an on-demand backend never knows.
    size_t get_read_count() const;

    size_t get_filtered_count() const;

    /// Queries served from the cache rather than the backend.
    size_t get_cache_hits() const;
    size_t get_cache_misses() const;

    /// Index entries examined while answering site queries. A read can appear under several of
    /// the nodes a site asks about, so the ratio of this to the reads delivered says how much of
    /// the index a site reads again. Counted per site query.
    size_t get_scanned_count() const;
    size_t get_delivered_count() const;

    /// Site queries that crossed a window boundary and so were fetched uncached, and
    /// the total node-ID span they covered.
    size_t get_straddle_count() const;
    size_t get_straddle_nodes() const;
    /// Node IDs those sites actually asked about, as opposed to the span they were
    /// collapsed to. The gap between the two is over-fetching.
    size_t get_straddle_wanted() const;

    /// Windows fetched whole, for sites or for read-start tallies, and the reads they held.
    size_t get_whole_fetches() const;
    size_t get_whole_fetch_reads() const;


protected:

    /// window_size is in node IDs; cache_entries is how many windows the cache holds in all.
    WindowedSiteReadSource(const SiteReadFilter& filter, size_t window_size,
                           size_t cache_entries);

    /// Visit every read with a mapping onto a node in any of the inclusive ranges,
    /// having applied the filter. Each read must be visited at most once. Must be safe
    /// to call concurrently: implementations own their per-thread resources.
    ///
    /// Ranges rather than one span, because a snarl's nodes can be spread thinly over a wide
    /// span of IDs.
    ///
    /// The iteratee takes a mutable reference: the alignment handed over is the
    /// backend's per-record scratch, and the caller may move from it.
    virtual void fetch_span(const vector<pair<nid_t, nid_t>>& ranges,
                            const function<void(Alignment&)>& iteratee) const = 0;

    /// Apply the filter, counting rejections. For subclasses to call on each
    /// candidate read before handing it to fetch_span's iteratee.
    bool passes_filter(const Alignment& aln) const;

    /// Count a read as fetched. Separate from passes_filter so a subclass can decide
    /// the order in which it filters and counts.
    void count_fetched() const;

    SiteReadFilter filter;

private:

    /// One mapping of one read in a window, as the window's index holds it.
    struct IndexEntry {
        nid_t node = 0;
        uint32_t read = 0;
        uint32_t mapping = 0;
        bool operator<(const IndexEntry& other) const {
            if (node != other.node) return node < other.node;
            if (read != other.read) return read < other.read;
            return mapping < other.mapping;
        }
    };

    /// One cached window fetch. Not changed once fetched, so threads read it without a lock.
    struct CacheEntry {
        size_t window = 0;
        vector<Alignment> reads;

        /// Every mapping in the window, sorted by node ID, so that a site finds the reads it
        /// wants, and which of their mappings, by binary search rather than by testing each read.
        /// Built once per fetch and used by every site in the window.
        vector<IndexEntry> node_index;

        /// Read offset before each mapping, for every read in the window end to end:
        /// read r occupies `[offset_start[r], offset_start[r + 1])`, which is one entry
        /// longer than its path so the last mapping's end is readable too.
        vector<uint32_t> offsets;
        vector<uint32_t> offset_start;
    };

    /// A window in the cache. It has no entry while a thread is fetching it.
    struct CacheSlot {
        shared_ptr<const CacheEntry> entry;
        /// When a query last used the window, for evicting the least recently used.
        size_t last_use = 0;
    };

    /// The window's reads, from the cache or fetched now. A thread wanting a window that another
    /// thread is fetching waits for that fetch rather than fetching it again. Sets `was_fetched` if
    /// this call did the fetch. The entry stays valid while the caller holds it, even if the
    /// cache drops it meanwhile.
    shared_ptr<const CacheEntry> get_window(size_t window, bool& was_fetched) const;

    /// Hand the entry's reads that touch the ranges to the caller, in the order they
    /// were fetched. Reads are found through the entry's node index, which lists a read
    /// under node n exactly when it has a mapping onto n.
    void deliver(const CacheEntry& entry,
                 const vector<pair<nid_t, nid_t>>& ranges,
                 const function<void(const SiteRead&)>& iteratee) const;

    /// Append this read's mappings to a window's index and its read-offset table.
    static void index_read(const Alignment& aln, uint32_t read_index, CacheEntry& entry);

    /// Which window a node ID falls in.
    size_t window_of(nid_t id) const;

    /// Does the read touch any node in the ranges?
    static bool touches(const Alignment& aln, const vector<pair<nid_t, nid_t>>& ranges);

    /// Fetch a whole window from the backend and index it.
    CacheEntry load_window(size_t window) const;

    /// Reads beginning on one node with one mapping quality: how many, and their total
    /// sequence length.
    struct StartTally {
        nid_t node = 0;
        int32_t mapq = 0;
        size_t reads = 0;
        size_t bases = 0;
    };

    /// Tally the read starts of a fetched window, sorted by node, then mapping quality. A read
    /// is tallied under the window holding its first node, though a window's fetch also returns
    /// reads that begin in other windows, so each read is tallied once.
    vector<StartTally> tally_starts(const CacheEntry& entry) const;

    /// The window's read-start tallies, fetching the window if it has never been fetched. The
    /// tallies are kept for the whole run and never changed, so the reference stays valid.
    const vector<StartTally>& window_starts(size_t window) const;

    size_t window_size;
    size_t cache_entries;

    /// Cached windows by window index, holding at most cache_entries fetched windows.
    /// unordered_map does not move its elements, so a reference to one stays valid while
    /// others are added.
    mutable unordered_map<size_t, CacheSlot> cache;
    /// Counts cache uses, to order them for eviction.
    mutable size_t cache_clock = 0;

    mutable atomic<size_t> fetched{0};
    mutable atomic<size_t> filtered{0};
    mutable atomic<size_t> cache_hits{0};
    mutable atomic<size_t> cache_misses{0};

    // Added once per site query rather than once per read, to keep atomic increments off the
    // per-read path.
    mutable atomic<size_t> scanned{0};
    mutable atomic<size_t> delivered{0};

    // Sites whose span crosses a window boundary, and the total ID span they asked for. These
    // bypass the cache.
    mutable atomic<size_t> straddles{0};
    mutable atomic<size_t> straddle_nodes{0};
    mutable atomic<size_t> straddle_wanted{0};
    mutable atomic<size_t> whole_fetches{0};
    mutable atomic<size_t> whole_fetch_reads{0};

    /// Read-start tallies by window, added when a window is first fetched and kept for the
    /// whole run.
    mutable unordered_map<size_t, vector<StartTally>> starts;

    /// Guards `cache`, `cache_clock` and `starts`. Fetches run without it.
    mutable std::mutex cache_mutex;
    /// Signalled when a fetch finishes or fails.
    mutable std::condition_variable cache_filled;
};

/**
 * Reads fetched on demand from a sorted GAM and its `.gai` index (`vg gamsort -i`), so memory
 * is bounded by what one window needs.
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

    /// gam_filename must be sorted (`vg gamsort -i`). window_size is in node IDs, and
    /// cache_entries is how many windows the cache shared by all threads holds.
    IndexedGamSiteReadSource(const string& gam_filename, const string& index_filename,
                             const SiteReadFilter& filter = SiteReadFilter(),
                             size_t window_size = 256,
                             size_t cache_entries = 2);

protected:

    void fetch_span(const vector<pair<nid_t, nid_t>>& ranges,
                    const function<void(Alignment&)>& iteratee) const;

private:

    /// Per-thread cursor. Mutable because fetch_span is logically const but seeks.
    struct ThreadState {
        unique_ptr<ifstream> stream;
        unique_ptr<GAMIndex::cursor_t> cursor;
    };

    ThreadState& thread_state() const;

    string gam_filename;
    unique_ptr<GAMIndex> index;

    mutable vector<ThreadState> threads;
};

/**
 * Reads fetched on demand from a GAF-Base database, by running `gbz-base query`.
 *
 * <https://github.com/jltsiren/gbz-base> stores alignments in SQLite and returns the reads that
 * overlap a set of nodes, without over-fetching. GAF-Base has no C API, and its file format may
 * change behind a version check, so we run its own binary rather than decode the format here.
 * This needs `gbz-base` at run time only when --gaf-base is given, and nothing at build time. The
 * binary writes GAF text, which is parsed, filtered, windowed and cached like the other sources.
 *
 * Starting a process is slow, so this fetches one window of node IDs per process.
 */
class GafBaseSiteReadSource : public WindowedSiteReadSource {
public:

    /// graph must outlive this, and must be the graph the alignments were made
    /// against: it supplies node lengths and sequences to turn GAF back into
    /// Alignments, and says which node IDs in a window actually exist.
    ///
    /// gbz_filename is a GBZ or a GBZ-Base (`gbz-base construct`). Prefer the
    /// latter: a plain GBZ is loaded in full on every query, while a GBZ-Base is
    /// random-access, and this issues many queries.
    ///
    /// window_size is in node IDs, and cache_entries is how many windows the cache shared by
    /// all threads holds.
    GafBaseSiteReadSource(const HandleGraph& graph,
                          const string& gaf_base_filename,
                          const string& gbz_filename,
                          const SiteReadFilter& filter = SiteReadFilter(),
                          size_t window_size = 256,
                          size_t cache_entries = 2,
                          const string& binary = "gbz-base");

    ~GafBaseSiteReadSource();

    /// Subprocesses spawned. The number that matters for run time, since each one
    /// costs milliseconds no matter how few reads it returns.
    size_t get_query_count() const;

    /// Run one query up front to check the databases are readable and agree with the
    /// graph, so a broken setup fails immediately rather than on the first snarl in
    /// a worker thread. Throws with an actionable message on failure.
    void check_setup() const;

    /// Whether gbz-base is handed the databases as SQLite immutable URIs (see `immutable`).
    bool immutable_databases() const { return immutable; }

    /// Where the time of the gbz-base queries went, summed over every query so far: how long the
    /// calling threads waited for their children, the CPU the children spent in their own code
    /// (user) and in the kernel (system), and the time spent parsing the GAF they wrote. Collected
    /// with wait4() and two clock reads per query, so it costs nothing measurable.
    struct QueryTotals {
        size_t queries = 0;
        double wait_s = 0;
        double child_user_s = 0;
        double child_system_s = 0;
        double parse_s = 0;
        size_t gaf_bytes = 0;
    };
    QueryTotals get_query_totals() const;

private:

    /// Per-thread GAF output files, one per query the thread can have in flight.
    /// Reused across queries rather than created per query, since creating one is a
    /// filesystem round trip.
    struct ThreadState {
        vector<string> gaf_paths;
    };

    ThreadState& thread_state() const;

    /// This thread's output file for the given in-flight slot, created on first use.
    const string& gaf_path(ThreadState& state, size_t slot) const;

    /// One `gbz-base query` child, spawned and not yet reaped.
    struct PendingQuery {
        pid_t pid = 0;
        /// By value, not a pointer into the thread's slot table: spawning a later slot
        /// can grow that vector and move what a pointer was aimed at.
        string gaf_path;
        string err_path;
        /// For QueryTotals and the per-query log.
        std::chrono::steady_clock::time_point started;
        double started_epoch = 0;
        size_t n_nodes = 0;
        nid_t first_node = 0;
        nid_t last_node = 0;
        /// FNV-1a over the node IDs, so the same node set fetched twice shows up in the log.
        uint64_t nodes_hash = 0;
    };

    void fetch_span(const vector<pair<nid_t, nid_t>>& ranges,
                    const function<void(Alignment&)>& iteratee) const;

    /// Start `gbz-base query` for these node IDs, writing to this thread's slot-th
    /// output file. Does not wait: the caller reaps it with reap_query, which is what
    /// lets a query that has to be split run its pieces at the same time rather than
    /// one after another. Throws if the child cannot be started.
    PendingQuery spawn_query(ThreadState& state, size_t slot,
                             const vector<nid_t>& nodes) const;

    /// Wait for a spawned query and parse the GAF it wrote, applying the filter.
    /// Returns the number of records parsed. Throws if the child failed.
    size_t reap_query(PendingQuery& pending,
                      const function<void(Alignment&)>& iteratee) const;

    /// Spawn and reap in one step. Throws on failure.
    size_t run_query(ThreadState& state, const vector<nid_t>& nodes,
                     const function<void(Alignment&)>& iteratee) const;

    /// run_query, but reporting and exiting instead of throwing. For use during
    /// calling, which happens inside an OpenMP parallel region.
    void run_query_or_die(ThreadState& state, const vector<nid_t>& nodes,
                          const function<void(Alignment&)>& iteratee) const;

    const HandleGraph& graph;
    string gaf_base_filename;
    string gbz_filename;
    string binary;

    /// SQLite takes a POSIX advisory lock on the database file around every read transaction, and
    /// gbz-base looks up each node of a query in a transaction of its own. With one gbz-base per
    /// calling thread, all of those lock and unlock calls land on the same file, and on a whole
    /// genome the kernel's bookkeeping for them took most of every query's time. vg only reads the
    /// databases, so by default it hands them to gbz-base as SQLite immutable URIs, which take no
    /// locks at all. Setting VG_GAFBASE_LOCKING=1 in the environment passes the plain paths instead.
    bool immutable = true;
    /// The database arguments gbz-base is given: the plain paths, or `file:<path>?immutable=1`.
    string gbz_argument;
    string gaf_base_argument;
    /// `file:<absolute path>?immutable=1` for a SQLite database, or the path unchanged for any other
    /// file, such as a plain GBZ given as the graph.
    static string immutable_uri(const string& path);

    /// QueryTotals, in microseconds.
    mutable atomic<uint64_t> wait_us{0};
    mutable atomic<uint64_t> child_user_us{0};
    mutable atomic<uint64_t> child_system_us{0};
    mutable atomic<uint64_t> parse_us{0};
    mutable atomic<uint64_t> gaf_bytes_total{0};
    /// One line per query, written when VG_GAFBASE_QUERY_LOG names a file. The columns are in its
    /// header line.
    FILE* query_log = nullptr;
    mutable std::mutex query_log_mutex;

    /// The most node IDs one gbz-base command line can hold: a quarter of sysconf(_SC_ARG_MAX),
    /// leaving room for the environment and the fixed arguments, over the bytes each node costs
    /// ("-n" and its digits). Computed once.
    static size_t argv_node_budget();

    /// Node IDs per gbz-base process. A larger query is split into pieces run side by side.
    size_t max_query_nodes = argv_node_budget();

    mutable vector<ThreadState> threads;
    mutable atomic<size_t> queries{0};
    /// Reads one query returned more than once, identified by name and start, and dropped. A
    /// split query returns a read that spans two of its pieces from both.
    mutable atomic<size_t> duplicates_dropped{0};

public:
    /// See duplicates_dropped.
    size_t get_duplicate_count() const { return duplicates_dropped.load(); }

private:
};

}

#endif
