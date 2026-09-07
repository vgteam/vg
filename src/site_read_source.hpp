#ifndef VG_SITE_READ_SOURCE_HPP_INCLUDED
#define VG_SITE_READ_SOURCE_HPP_INCLUDED

/** \file site_read_source.hpp
 *
 * Random-access sources of read alignments by graph locality, for callers that
 * need to reason about the individual reads overlapping a site.
 *
 * `vg call` historically consumed only a `vg pack` coverage index, which has no
 * read-level information at all. Read-level genotyping needs the actual reads
 * overlapping each snarl, which is what these classes provide.
 */

#include <atomic>
#include <fstream>
#include <functional>
#include <memory>
#include <string>
#include <unordered_map>
#include <vector>

#include <vg/vg.pb.h>

#include "handle.hpp"
#include "stream_index.hpp"

namespace vg {

using namespace std;

/**
 * Which reads are eligible to be used as evidence.
 *
 * These defaults mirror what `vg pack` does when building the support index,
 * because allele *enumeration* is pack-driven while *genotyping* is
 * read-source-driven. If the two disagree about which reads exist, the
 * traversals being scored were selected on evidence the genotyper cannot see,
 * and vice versa.
 *
 * Note `vg pack -Q` sets its minimum mapping quality and minimum base quality
 * from a single value; here they are separate, since base quality is used per
 * base by the scoring rather than as a read-level filter.
 */
struct SiteReadFilter {
    /// Drop reads with mapping quality below this. The `vg pack -Q` equivalent.
    int min_mapq = 0;
    /// Drop alignments flagged secondary. Nothing in vg's pack path does this,
    /// but a secondary alignment of a read already counted would break the
    /// per-read independence the genotype likelihood assumes.
    bool skip_secondary = true;
    /// Drop reads with no mappings at all.
    bool skip_unmapped = true;
};

/**
 * A read handed to a site query, with enough of an index to work on it without
 * re-walking the whole alignment.
 *
 * The index exists because a site's work is proportional to the read's *overlap* with
 * the site, while finding that overlap by scanning is proportional to the read's
 * *length*. On a 150 bp read the two are the same thing. On a 33 kb ONT read, which
 * carries some 3,000 mappings and is handed to the ~85 sites it spans, the scan is the
 * run time: four separate passes over every mapping -- selecting the site's steps,
 * slicing the site's stretch out for the reverse-complement, and resolving each of the
 * two anchor pins -- for a site that wants a handful of them.
 *
 * `mappings` and `read_offsets` are a hint, not a definition. A source that cannot
 * index cheaply leaves them null, and consumers must then fall back to walking the
 * alignment; a consumer must never treat a null index as "this read touches nothing".
 */
struct SiteRead {
    /// The alignment itself. Valid only for the duration of the callback.
    const Alignment* aln = nullptr;

    /// Ascending indices of the mappings whose node lies in the queried ranges, or
    /// null if this source does not index.
    const uint32_t* mappings = nullptr;
    size_t mapping_count = 0;

    /// Read offset -- summed `to_length` -- before each mapping, indexed by mapping
    /// index and one longer than the path, or null if this source does not index.
    const uint32_t* read_offsets = nullptr;

    /// Whether the index is present. Both halves are supplied together or not at all,
    /// so consumers have one condition to branch on.
    bool indexed() const { return mappings != nullptr && read_offsets != nullptr; }
};

/**
 * Random-access source of read alignments by graph locality.
 *
 * Implementations must be safe for concurrent read access: GraphCaller visits
 * snarls in parallel and in arbitrary order, so several threads will be asking
 * for reads at different sites at the same time. Any loading or index building
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

    /// The same, for consumers with no use for the index. Separately named rather
    /// than overloaded: a lambda converts to either `function` type, so an overload
    /// would make every call site ambiguous.
    void for_each_alignment(const vector<pair<nid_t, nid_t>>& ranges,
                            const function<void(const Alignment&)>& iteratee) const {
        for_each_read(ranges, [&](const SiteRead& read) { iteratee(*read.aln); });
    }


    /// How many reads this source holds or can see, for logging. May be 0 if
    /// the backend cannot cheaply say.
    virtual size_t get_read_count() const = 0;

    /// Width, in node IDs, of the locality this source fetches and caches around
    /// a request. 0 means "no such locality" -- an in-memory source answers each
    /// request exactly.
    ///
    /// Exposed so a caller can ask for a *neighbourhood* rather than a site and
    /// get it from the same cache entry the site already populated. That is what
    /// makes a local depth rate free: the reads are fetched either way, and a
    /// window is wide enough to be a meaningful denominator where a snarl is not.
};

/**
 * Reads held in memory, bucketed by the node IDs they touch.
 *
 * One streaming pass over a GAM or GAF, no index, no new dependency, and
 * correct by construction: there is no over-fetching to filter and no cursor
 * lifetime to get wrong. That makes it the right thing to develop and test
 * against, and it stays useful for regional and chunked workflows.
 *
 * The honest limit is memory: a whole-genome GAM at reasonable depth will not
 * fit. get_read_count() is logged so that limit is visible rather than
 * discovered as an out-of-memory kill.
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

    /// Retain a single read directly, applying the same filter as loading would.
    /// Lets callers assemble a source without going through a file, which is what
    /// makes the scoring unit-testable.
    void add(const Alignment& aln, const Filter& filter = Filter());

    /// Reads come through unindexed: this source keeps every read for the whole run,
    /// so a persistent per-read index would add several bytes per mapping to a peak
    /// that is already the reason to prefer an on-demand backend. Consumers fall back
    /// to walking the alignment, which is what they did before the index existed.
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
 * Base for on-demand backends: quantises fetches to fixed windows of node IDs and
 * caches the last few windows per thread.
 *
 * Every on-demand backend faces the same problem. A query costs far more than one
 * site's worth of reads -- because the backend over-fetches, or because it is a
 * process spawn -- so issuing one query per snarl rescans or respawns endlessly.
 * Quantising to windows fixes that, but only if the caller visits sites in node-ID
 * order, so that each window is asked for while it is still resident. That ordering
 * is not free and is arranged separately; see GraphCaller::set_node_id_ordering.
 *
 * Subclasses supply one primitive, fetch_span(), and manage whatever per-thread
 * resources it needs. Everything above it -- window arithmetic, the cache, the
 * boundary-straddling bypass, and narrowing a window back down to the ranges the
 * caller actually asked about -- lives here, so the two backends cannot drift apart
 * in how they interpret a query.
 */
class WindowedSiteReadSource : public SiteReadSource {
public:

    void for_each_read(const vector<pair<nid_t, nid_t>>& ranges,
                       const function<void(const SiteRead&)>& iteratee) const final;

    /// Reads actually fetched from the backend so far, across all threads. Not the
    /// size of the read set, which an on-demand backend never knows.
    size_t get_read_count() const;

    size_t get_filtered_count() const;

    /// Queries served from the cache rather than the backend. Low hit rates mean the
    /// caching assumption above does not hold for this workload, which is worth
    /// knowing rather than guessing.
    size_t get_cache_hits() const;
    size_t get_cache_misses() const;

    /// Index entries examined while answering site queries, against reads actually
    /// handed to the caller. A site's query walks the window's node index over the
    /// ranges it asked about, and one read can appear under several of those nodes,
    /// so the ratio between these is how much of the index a site re-reads -- close to
    /// 1 for short reads, higher for long ones that cross a site many times. Both are
    /// counted per site query, so a read in a window visited by many sites counts
    /// many times.
    size_t get_scanned_count() const;
    size_t get_delivered_count() const;

    /// Site queries that crossed a window boundary and so were fetched uncached, and
    /// the total node-ID span they covered.
    size_t get_straddle_count() const;
    size_t get_straddle_nodes() const;
    /// Node IDs those sites actually asked about, as opposed to the span they were
    /// collapsed to. The gap between the two is over-fetching.
    size_t get_straddle_wanted() const;


protected:

    WindowedSiteReadSource(const SiteReadFilter& filter, size_t window_size,
                           size_t cache_entries);

    /// Visit every read with a mapping onto a node in any of the inclusive ranges,
    /// having applied the filter. Each read must be visited at most once. Must be safe
    /// to call concurrently: implementations own their per-thread resources.
    ///
    /// Ranges rather than one span because a snarl's contents can be extremely sparse
    /// in ID space. On chr20, 215 sites spanning 13.2 M node IDs between them wanted
    /// only 133 k of those IDs -- a 99x over-fetch if the span is used instead. Both
    /// backends address ranges natively, so this costs them nothing.
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

    /// One cached window fetch.
    struct CacheEntry {
        size_t window = 0;
        bool valid = false;
        vector<Alignment> reads;

        /// Every mapping in the window, sorted by node ID, so a site finds the reads it
        /// wants -- and *which* of their mappings it wants -- by binary search instead of
        /// asking each read in turn whether it touches the site. The distinction is
        /// entirely one of read length: a 33 kb ONT read carries some 3,000 mappings, and
        /// testing it against a site walks half of them on average, so the scan that
        /// costs a short read a few hundred comparisons per window costs a long read
        /// millions. Built once per fetch and reused by every site inside the window.
        vector<IndexEntry> node_index;

        /// Read offset before each mapping, for every read in the window end to end:
        /// read r occupies `[offset_start[r], offset_start[r + 1])`, which is one entry
        /// longer than its path so the last mapping's end is readable too.
        vector<uint32_t> offsets;
        vector<uint32_t> offset_start;
    };

    /// Per-thread cache. Mutable because for_each_read is logically const but may
    /// populate the cache.
    struct CacheState {
        vector<CacheEntry> cache;
        size_t next_evict = 0;
        bool initialized = false;
    };

    CacheState& cache_state() const;

    /// Hand the entry's reads that touch the ranges to the caller, in the order they
    /// were fetched. Reads are found through the entry's node index, which is exact:
    /// a read is listed under node n precisely when it has a mapping onto n, so this
    /// selects the same set touches() would and needs no second adjudication.
    void deliver(const CacheEntry& entry,
                 const vector<pair<nid_t, nid_t>>& ranges,
                 const function<void(const SiteRead&)>& iteratee) const;

    /// Append this read's mappings to a window's index and its read-offset table.
    static void index_read(const Alignment& aln, uint32_t read_index, CacheEntry& entry);

    /// Which window a node ID falls in.
    size_t window_of(nid_t id) const;

    /// Does the read touch any node in the ranges?
    static bool touches(const Alignment& aln, const vector<pair<nid_t, nid_t>>& ranges);

    size_t window_size;
    size_t cache_entries;

    mutable vector<CacheState> caches;
    mutable atomic<size_t> fetched{0};
    mutable atomic<size_t> filtered{0};
    mutable atomic<size_t> cache_hits{0};
    mutable atomic<size_t> cache_misses{0};

    // Accumulated once per site query, not once per read: these count into the
    // hundreds of millions, and an atomic increment on that path would cost more than
    // the work it is measuring.
    mutable atomic<size_t> scanned{0};
    mutable atomic<size_t> delivered{0};

    // Sites whose span crosses a window boundary, and the total ID span they asked
    // for. These bypass the cache entirely, so if the span is large relative to the
    // window they are where the backend queries actually go.
    mutable atomic<size_t> straddles{0};
    mutable atomic<size_t> straddle_nodes{0};
    mutable atomic<size_t> straddle_wanted{0};
};

/**
 * Reads fetched on demand from a sorted GAM plus its `.gai` index.
 *
 * This is the backend that makes whole-genome work possible with no new dependency:
 * memory is bounded by what one window needs, rather than by the size of the read
 * set.
 *
 * Two properties of the index shape the implementation, and both come from
 * StreamIndex rather than from any choice here:
 *
 * * **One cursor per thread.** Concurrent `find()` calls are documented safe, but a
 *   cursor seeks, so it cannot be shared. Cursors are created lazily per thread, the
 *   pattern `vg chunk` uses.
 * * **It over-fetches.** The index can only give group start offsets, so a query
 *   scans groups and stops when a group's minimum node ID is too large. That is what
 *   the windowing in the base class is for.
 */
class IndexedGamSiteReadSource : public WindowedSiteReadSource {
public:

    /// gam_filename must be sorted (`vg gamsort -i`). window_size is in node IDs.
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
 * <https://github.com/jltsiren/gbz-base> stores alignments column-compressed in
 * SQLite and can return the reads overlapping a set of nodes. That is the access
 * pattern a caller wants, and unlike the GAM index it does not over-fetch.
 *
 * **This shells out to a binary rather than linking a library, and that is
 * deliberate.** GAF-Base has no C API, and its on-disk format is documented as able
 * to change without warning behind an exact-match version check. Reimplementing the
 * decoder in C++ would mean tracking a moving format; letting upstream's own binary
 * do the decoding costs us nothing when it changes. So this adds a *runtime*
 * dependency on `gbz-base` being on the PATH, and no build dependency at all --
 * nothing links, and users who do not pass --gaf-base never notice. If a C shim
 * appears upstream, it replaces run_query() and nothing else: both paths consume GAF
 * text, so the parsing, filtering, windowing, and caching are already shared.
 *
 * The cost of a spawn is why this derives from WindowedSiteReadSource. One process
 * per snarl would be hopeless; one per window of node IDs, visited in order, is a
 * handful of spawns for a whole contig.
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

    /// Node IDs per subprocess. A node list becomes argv, so it cannot grow without
    /// bound; queries larger than this are split, and their results de-duplicated.
    size_t max_query_nodes = 4096;

    mutable vector<ThreadState> threads;
    mutable atomic<size_t> queries{0};
};

}

#endif
