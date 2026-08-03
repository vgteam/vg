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
                               const function<void(const Alignment&)>& iteratee) const = 0;

    /// Convenience wrapper collecting for_each_read into a vector. Copies, so
    /// prefer for_each_read on the hot path; this is for tests and diagnostics.
    vector<Alignment> get_reads(const vector<pair<nid_t, nid_t>>& ranges) const;

    /// How many reads this source holds or can see, for logging. May be 0 if
    /// the backend cannot cheaply say.
    virtual size_t get_read_count() const = 0;
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

    void for_each_read(const vector<pair<nid_t, nid_t>>& ranges,
                       const function<void(const Alignment&)>& iteratee) const;

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
 * Reads fetched on demand from a sorted GAM plus its `.gai` index.
 *
 * This is the backend that makes whole-genome work possible: memory is bounded by
 * what one site needs plus a small cache, rather than by the size of the read set.
 *
 * Two properties of the index shape the implementation, and both come from
 * StreamIndex rather than from any choice here:
 *
 * * **One cursor per thread.** Concurrent `find()` calls are documented safe, but a
 *   cursor seeks, so it cannot be shared. Cursors are created lazily per thread, the
 *   pattern `vg chunk` uses.
 * * **It over-fetches.** The index can only give group start offsets, so a query
 *   scans groups and stops when a group's minimum node ID is too large. Issuing one
 *   query per snarl across millions of snarls therefore rescans the same GAM groups
 *   repeatedly. That is what the cache below is for: because GraphCaller visits a
 *   top-level snarl and then its descendants on the same thread, a child's node range
 *   is usually inside its parent's, so the parent's fetch can serve it.
 */
class IndexedGamSiteReadSource : public SiteReadSource {
public:

    /// gam_filename must be sorted (`vg gamsort -i`); index_filename defaults to
    /// gam_filename + ".gai".
    IndexedGamSiteReadSource(const string& gam_filename, const string& index_filename,
                             const SiteReadFilter& filter = SiteReadFilter(),
                             size_t cache_entries = 2);

    void for_each_read(const vector<pair<nid_t, nid_t>>& ranges,
                       const function<void(const Alignment&)>& iteratee) const;

    /// Reads actually fetched from the index so far, across all threads. Not the size
    /// of the read set, which this backend never knows.
    size_t get_read_count() const;

    size_t get_filtered_count() const;

    /// Queries served from the cache rather than the index. Low hit rates mean the
    /// caching assumption above does not hold for this workload, which is worth
    /// knowing rather than guessing.
    size_t get_cache_hits() const;
    size_t get_cache_misses() const;

private:

    /// One cached fetch: the span it covers and the reads it found.
    struct CacheEntry {
        nid_t min_id = 0;
        nid_t max_id = 0;
        bool valid = false;
        vector<Alignment> reads;
    };

    /// Per-thread cursor and cache. Mutable because for_each_read is logically const
    /// but must seek and may populate the cache.
    struct ThreadState {
        unique_ptr<ifstream> stream;
        unique_ptr<GAMIndex::cursor_t> cursor;
        vector<CacheEntry> cache;
        size_t next_evict = 0;
    };

    ThreadState& thread_state() const;

    /// True if every requested range lies inside the entry's span, so filtering the
    /// entry's reads gives exactly the right answer.
    static bool covers(const CacheEntry& entry, const vector<pair<nid_t, nid_t>>& ranges);

    /// Does the read touch any node in the ranges?
    static bool touches(const Alignment& aln, const vector<pair<nid_t, nid_t>>& ranges);

    string gam_filename;
    unique_ptr<GAMIndex> index;
    SiteReadFilter filter;
    size_t cache_entries;

    mutable vector<ThreadState> threads;
    mutable atomic<size_t> fetched{0};
    mutable atomic<size_t> filtered{0};
    mutable atomic<size_t> cache_hits{0};
    mutable atomic<size_t> cache_misses{0};
};

}

#endif
