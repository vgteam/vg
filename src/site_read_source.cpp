#include "site_read_source.hpp"

#include <algorithm>
#include <cerrno>
#include <cstring>
#include <sstream>
#include <unordered_set>

#include <fcntl.h>
#include <spawn.h>
#include <sys/wait.h>
#include <unistd.h>

/// The environment to hand a spawned gbz-base. Declared rather than included because
/// the header that provides it differs between platforms; the symbol does not.
extern char** environ;

#include <vg/io/alignment_io.hpp>
#include <vg/io/stream.hpp>

#include <omp.h>

#include "path.hpp"
#include "utility.hpp"

namespace vg {

using namespace std;


void SiteReadSource::for_each_read_start(
    const vector<nid_t>& nodes,
    const function<void(int32_t mapq, size_t reads, size_t bases)>& iteratee) const {

    // The nodes as coalesced ID ranges: the source visits each read once however many ranges it
    // touches.
    vector<pair<nid_t, nid_t>> ranges;
    for (nid_t id : nodes) {
        if (!ranges.empty() && ranges.back().second + 1 == id) {
            ranges.back().second = id;
        } else {
            ranges.emplace_back(id, id);
        }
    }
    if (ranges.empty()) {
        return;
    }
    for_each_alignment(ranges, [&](const Alignment& aln) {
        // The fetch also returns reads that only pass through the nodes.
        if (aln.path().mapping_size() == 0) {
            return;
        }
        nid_t start_node = aln.path().mapping(0).position().node_id();
        if (!std::binary_search(nodes.begin(), nodes.end(), start_node)) {
            return;
        }
        iteratee(aln.mapping_quality(), 1, aln.sequence().size());
    });
}

void InMemorySiteReadSource::add_read(const Alignment& aln, const Filter& filter) {
    // Secondary and unmapped alignments are always dropped. A secondary alignment repeats a read
    // that is already counted, which would break the independence of reads that the genotype
    // likelihood assumes.
    if (aln.is_secondary()) {
        ++filtered_count;
        return;
    }
    if (aln.mapping_quality() < filter.min_mapq) {
        ++filtered_count;
        return;
    }
    if (aln.path().mapping_size() == 0) {
        ++filtered_count;
        return;
    }

    size_t read_index = reads.size();
    reads.push_back(aln);

    // Index under every node the read touches. A read can visit the same node more than once (a
    // cycle, or a snarl traversed twice), so it is indexed once per node.
    nid_t prev_node = 0;
    bool have_prev = false;
    for (const auto& mapping : reads[read_index].path().mapping()) {
        nid_t node_id = mapping.position().node_id();
        if (have_prev && node_id == prev_node) {
            // Consecutive mappings on one node: already indexed.
            continue;
        }
        auto& bucket = reads_by_node[node_id];
        if (bucket.empty() || bucket.back() != read_index) {
            bucket.push_back(read_index);
        }
        prev_node = node_id;
        have_prev = true;
    }
}

void InMemorySiteReadSource::add(const Alignment& aln, const Filter& filter) {
    add_read(aln, filter);
}

void InMemorySiteReadSource::load_gam(const string& filename, const Filter& filter) {
    // Serial: filling a shared map in parallel would need a mutex, and this runs once at startup.
    get_input_file(filename, [&](istream& in) {
        vg::io::for_each<Alignment>(in, [&](Alignment& aln) {
            add_read(aln, filter);
        });
    });
}

void InMemorySiteReadSource::load_gaf(const HandleGraph& graph, const string& filename, const Filter& filter) {
    vg::io::gaf_unpaired_for_each(graph, filename, [&](Alignment& aln) {
        add_read(aln, filter);
    });
}

void InMemorySiteReadSource::for_each_read(const vector<pair<nid_t, nid_t>>& ranges,
                                          const function<void(const SiteRead&)>& iteratee) const {

    // Collect the matching read indices first, so a read touching several of
    // the ranges is only visited once.
    unordered_set<size_t> seen;

    for (const auto& range : ranges) {
        // Ranges are inclusive on both ends. Iterating IDs is right for the
        // small ranges a snarl produces; a range spanning much of the graph
        // would be better served by scanning the map, but no caller does that.
        for (nid_t node_id = range.first; node_id <= range.second; ++node_id) {
            auto it = reads_by_node.find(node_id);
            if (it == reads_by_node.end()) {
                continue;
            }
            for (size_t read_index : it->second) {
                if (seen.insert(read_index).second) {
                    SiteRead read;
                    read.aln = &reads[read_index];
                    iteratee(read);
                }
            }
        }
    }
}

size_t InMemorySiteReadSource::get_read_count() const {
    return reads.size();
}

size_t InMemorySiteReadSource::get_filtered_count() const {
    return filtered_count;
}

////////////////////////////////////////////////////////////////////////////////
// WindowedSiteReadSource
////////////////////////////////////////////////////////////////////////////////

WindowedSiteReadSource::WindowedSiteReadSource(const SiteReadFilter& filter,
                                               size_t window_size,
                                               size_t cache_entries)
    : filter(filter),
      window_size(max<size_t>(1, window_size)),
      cache_entries(max<size_t>(1, cache_entries)) {
}

size_t WindowedSiteReadSource::window_of(nid_t id) const {
    return (size_t)(id / (nid_t)window_size);
}

bool WindowedSiteReadSource::passes_filter(const Alignment& aln) const {
    if (aln.is_secondary()) {
        ++filtered;
        return false;
    }
    if (aln.mapping_quality() < filter.min_mapq) {
        ++filtered;
        return false;
    }
    if (aln.path().mapping_size() == 0) {
        ++filtered;
        return false;
    }
    return true;
}

void WindowedSiteReadSource::count_fetched() const {
    ++fetched;
}

bool WindowedSiteReadSource::touches(const Alignment& aln,
                                    const vector<pair<nid_t, nid_t>>& ranges) {
    for (const auto& mapping : aln.path().mapping()) {
        nid_t node_id = mapping.position().node_id();
        for (const auto& range : ranges) {
            if (node_id >= range.first && node_id <= range.second) {
                return true;
            }
        }
    }
    return false;
}

void WindowedSiteReadSource::for_each_read(
    const vector<pair<nid_t, nid_t>>& ranges,
    const function<void(const SiteRead&)>& iteratee) const {

    if (ranges.empty()) {
        return;
    }

    nid_t min_id = ranges.front().first;
    nid_t max_id = ranges.front().second;
    for (const auto& range : ranges) {
        min_id = min(min_id, range.first);
        max_id = max(max_id, range.second);
    }

    size_t first_window = window_of(min_id);
    size_t last_window = window_of(max_id);

    if (first_window != last_window) {
        // The site crosses a window boundary. Fetch its ranges directly rather than stitching
        // windows together: fetch_span emits each read once, and a region of unbounded size is not
        // worth caching. The ranges are passed as they are, rather than collapsed to one span,
        // since a snarl's nodes can be spread thinly over a wide span of IDs.
        ++cache_misses;
        ++straddles;
        straddle_nodes += (size_t)(max_id - min_id + 1);
        size_t wanted = 0;
        for (const auto& range : ranges) {
            wanted += (size_t)(range.second - range.first + 1);
        }
        straddle_wanted += wanted;
        size_t n_scanned = 0, n_delivered = 0;
        fetch_span(ranges, [&](Alignment& aln) {
            ++n_scanned;
            if (touches(aln, ranges)) {
                ++n_delivered;
                // Not indexed: this path serves the few sites too wide to cache, and an index for a
                // read seen once would cost as much as the walk it saves.
                SiteRead read;
                read.aln = &aln;
                iteratee(read);
            }
        });
        scanned += n_scanned;
        delivered += n_delivered;
        return;
    }

    // Usually served from the cache when sites are visited in node-ID order
    // (GraphCaller::set_node_id_ordering).
    bool was_fetched = false;
    shared_ptr<const CacheEntry> entry = get_window(first_window, was_fetched);
    if (was_fetched) {
        ++cache_misses;
    } else {
        ++cache_hits;
    }
    deliver(*entry, ranges, iteratee);
}

shared_ptr<const WindowedSiteReadSource::CacheEntry>
WindowedSiteReadSource::get_window(size_t window, bool& was_fetched) const {
    was_fetched = false;
    {
        unique_lock<std::mutex> lock(cache_mutex);
        while (true) {
            auto found = cache.find(window);
            if (found == cache.end()) {
                // No thread has the window: claim it, with no entry until it is fetched.
                cache.emplace(window, CacheSlot());
                break;
            }
            if (found->second.entry) {
                found->second.last_use = ++cache_clock;
                return found->second.entry;
            }
            // Another thread is fetching it.
            cache_filled.wait(lock);
        }
    }

    shared_ptr<CacheEntry> entry;
    bool tallied = false;
    try {
        entry = make_shared<CacheEntry>(load_window(window));
        {
            lock_guard<std::mutex> guard(cache_mutex);
            tallied = starts.count(window) > 0;
        }
    } catch (...) {
        // Release the claim, so that a thread waiting on this window fetches it itself rather
        // than waiting for ever.
        lock_guard<std::mutex> guard(cache_mutex);
        cache.erase(window);
        cache_filled.notify_all();
        throw;
    }
    // Only this thread can be fetching the window now, so no other thread can tally it first.
    vector<StartTally> tallies;
    if (!tallied) {
        tallies = tally_starts(*entry);
    }
    was_fetched = true;

    // Windows dropped from the cache are freed here, after the lock is released: freeing a
    // window's reads takes long enough that other threads would queue on the lock.
    vector<shared_ptr<const CacheEntry>> dropped;
    lock_guard<std::mutex> guard(cache_mutex);
    if (!tallied) {
        starts.emplace(window, std::move(tallies));
    }
    CacheSlot& slot = cache[window];
    slot.entry = entry;
    slot.last_use = ++cache_clock;

    // Drop the least recently used fetched windows beyond the cache's size. A window being
    // fetched has no entry and is not counted. A thread still reading a dropped window keeps it
    // alive through its own pointer.
    while (true) {
        size_t held = 0;
        auto oldest = cache.end();
        for (auto it = cache.begin(); it != cache.end(); ++it) {
            if (!it->second.entry) {
                continue;
            }
            ++held;
            if (it->first != window && (oldest == cache.end()
                                        || it->second.last_use < oldest->second.last_use)) {
                oldest = it;
            }
        }
        if (held <= cache_entries || oldest == cache.end()) {
            break;
        }
        dropped.push_back(std::move(oldest->second.entry));
        cache.erase(oldest);
    }
    cache_filled.notify_all();
    return entry;
}

WindowedSiteReadSource::CacheEntry WindowedSiteReadSource::load_window(size_t window) const {
    CacheEntry entry;
    entry.window = window;
    nid_t lo = (nid_t)(window * window_size);
    nid_t hi = lo + (nid_t)window_size - 1;
    fetch_span({{lo, hi}}, [&](Alignment& aln) {
        // Take ownership rather than copy the alignment. fetch_span hands out a mutable reference
        // so that this can move from it; the backends reuse one Alignment per record and clear it
        // before the next.
        index_read(aln, (uint32_t)entry.reads.size(), entry);
        entry.reads.push_back(std::move(aln));
    });
    entry.offset_start.push_back((uint32_t)entry.offsets.size());
    // Sorted once per fetch rather than per site query. There is one entry per mapping, so a read
    // that visits a node twice is listed twice.
    std::sort(entry.node_index.begin(), entry.node_index.end());
    ++whole_fetches;
    whole_fetch_reads += entry.reads.size();
    return entry;
}

vector<WindowedSiteReadSource::StartTally>
WindowedSiteReadSource::tally_starts(const CacheEntry& entry) const {
    vector<StartTally> tallies;
    for (const Alignment& aln : entry.reads) {
        // The fetch filtered the reads, and the filter drops reads with no mappings.
        nid_t node = aln.path().mapping(0).position().node_id();
        if (window_of(node) != entry.window) {
            // Tallied under its own window instead.
            continue;
        }
        StartTally tally;
        tally.node = node;
        tally.mapq = aln.mapping_quality();
        tally.reads = 1;
        tally.bases = aln.sequence().size();
        tallies.push_back(tally);
    }
    std::sort(tallies.begin(), tallies.end(), [](const StartTally& a, const StartTally& b) {
        return a.node != b.node ? a.node < b.node : a.mapq < b.mapq;
    });
    // Merge the reads sharing a node and a mapping quality.
    size_t kept = 0;
    for (size_t i = 0; i < tallies.size(); ++i) {
        if (kept > 0 && tallies[kept - 1].node == tallies[i].node
            && tallies[kept - 1].mapq == tallies[i].mapq) {
            tallies[kept - 1].reads += tallies[i].reads;
            tallies[kept - 1].bases += tallies[i].bases;
        } else {
            tallies[kept++] = tallies[i];
        }
    }
    tallies.resize(kept);
    return tallies;
}

const vector<WindowedSiteReadSource::StartTally>&
WindowedSiteReadSource::window_starts(size_t window) const {
    {
        lock_guard<std::mutex> guard(cache_mutex);
        auto found = starts.find(window);
        if (found != starts.end()) {
            return found->second;
        }
    }
    // Never fetched. A window's first fetch tallies it before the fetch is published, so once
    // get_window returns the tallies are in. The fetched window goes into the cache, where the
    // sites inside it will find it. Not a site query, so not counted as a cache hit or miss.
    bool was_fetched = false;
    get_window(window, was_fetched);
    lock_guard<std::mutex> guard(cache_mutex);
    return starts.at(window);
}

void WindowedSiteReadSource::for_each_read_start(
    const vector<nid_t>& nodes,
    const function<void(int32_t mapq, size_t reads, size_t bases)>& iteratee) const {

    size_t i = 0;
    while (i < nodes.size()) {
        // The nodes are sorted, so each window's nodes are consecutive.
        size_t window = window_of(nodes[i]);
        const vector<StartTally>& tallies = window_starts(window);
        for (; i < nodes.size() && window_of(nodes[i]) == window; ++i) {
            StartTally probe;
            probe.node = nodes[i];
            auto it = std::lower_bound(tallies.begin(), tallies.end(), probe,
                                       [](const StartTally& a, const StartTally& b) {
                                           return a.node < b.node;
                                       });
            for (; it != tallies.end() && it->node == nodes[i]; ++it) {
                iteratee(it->mapq, it->reads, it->bases);
            }
        }
    }
}

size_t WindowedSiteReadSource::get_whole_fetches() const {
    return whole_fetches.load();
}

size_t WindowedSiteReadSource::get_whole_fetch_reads() const {
    return whole_fetch_reads.load();
}

void WindowedSiteReadSource::index_read(const Alignment& aln, uint32_t read_index,
                                       CacheEntry& entry) {
    entry.offset_start.push_back((uint32_t)entry.offsets.size());
    uint32_t offset = 0;
    const Path& path = aln.path();
    for (int64_t i = 0; i < path.mapping_size(); ++i) {
        const Mapping& mapping = path.mapping(i);
        entry.node_index.push_back(IndexEntry{mapping.position().node_id(), read_index,
                                              (uint32_t)i});
        entry.offsets.push_back(offset);
        offset += (uint32_t)mapping_to_length(mapping);
    }
    // One past the end, so a consumer can read where the last mapping finishes without
    // a special case.
    entry.offsets.push_back(offset);
}

void WindowedSiteReadSource::deliver(const CacheEntry& entry,
                                     const vector<pair<nid_t, nid_t>>& ranges,
                                     const function<void(const SiteRead&)>& iteratee) const {
    // A local buffer rather than a thread-local one, since the iteratee is caller code that might
    // query the source again. A site names a handful of nodes, so the buffer is small.
    vector<pair<uint32_t, uint32_t>> hits;
    vector<uint32_t> mappings;

    for (const auto& range : ranges) {
        IndexEntry probe{range.first, 0, 0};
        auto it = std::lower_bound(entry.node_index.begin(), entry.node_index.end(), probe);
        for (; it != entry.node_index.end() && it->node <= range.second; ++it) {
            hits.emplace_back(it->read, it->mapping);
        }
    }

    // By read, then by mapping: reads are delivered in fetch order, and each read's mappings form
    // an ascending run the consumer can point into. There is one index entry per mapping and the
    // ranges do not overlap, so nothing repeats.
    std::sort(hits.begin(), hits.end());

    mappings.reserve(hits.size());
    for (const auto& hit : hits) {
        mappings.push_back(hit.second);
    }

    size_t n_delivered = 0;
    size_t i = 0;
    while (i < hits.size()) {
        size_t j = i;
        while (j < hits.size() && hits[j].first == hits[i].first) {
            ++j;
        }
        SiteRead read;
        read.aln = &entry.reads[hits[i].first];
        read.mappings = &mappings[i];
        read.mapping_count = j - i;
        read.read_offsets = &entry.offsets[entry.offset_start[hits[i].first]];
        iteratee(read);
        ++n_delivered;
        i = j;
    }

    // Added once per site query, to keep atomic increments off the per-entry path.
    scanned += hits.size();
    delivered += n_delivered;
}

size_t WindowedSiteReadSource::get_read_count() const {
    return fetched.load();
}

size_t WindowedSiteReadSource::get_filtered_count() const {
    return filtered.load();
}

size_t WindowedSiteReadSource::get_cache_hits() const {
    return cache_hits.load();
}

size_t WindowedSiteReadSource::get_cache_misses() const {
    return cache_misses.load();
}

size_t WindowedSiteReadSource::get_scanned_count() const {
    return scanned.load();
}

size_t WindowedSiteReadSource::get_delivered_count() const {
    return delivered.load();
}

size_t WindowedSiteReadSource::get_straddle_count() const {
    return straddles.load();
}

size_t WindowedSiteReadSource::get_straddle_nodes() const {
    return straddle_nodes.load();
}

size_t WindowedSiteReadSource::get_straddle_wanted() const {
    return straddle_wanted.load();
}

////////////////////////////////////////////////////////////////////////////////
// IndexedGamSiteReadSource
////////////////////////////////////////////////////////////////////////////////

IndexedGamSiteReadSource::IndexedGamSiteReadSource(const string& gam_filename,
                                                   const string& index_filename,
                                                   const SiteReadFilter& filter,
                                                   size_t window_size,
                                                   size_t cache_entries)
    : WindowedSiteReadSource(filter, window_size, cache_entries),
      gam_filename(gam_filename) {

    index.reset(new GAMIndex());
    get_input_file(index_filename, [&](istream& in) {
        index->load(in);
    });

    // One slot per thread, populated lazily: a cursor seeks, so it cannot be shared,
    // and opening one per thread up front would open files we may never use.
    threads.resize(max(1, get_thread_count()));
}

IndexedGamSiteReadSource::ThreadState& IndexedGamSiteReadSource::thread_state() const {
    int tid = omp_get_thread_num();
    if ((size_t)tid >= threads.size()) {
        // More threads than we sized for; grow rather than misbehave.
#pragma omp critical (indexed_gam_threads)
        if ((size_t)tid >= threads.size()) {
            threads.resize(tid + 1);
        }
    }
    ThreadState& state = threads[tid];
    if (!state.cursor) {
        state.stream.reset(new ifstream(gam_filename));
        if (!*state.stream) {
            throw runtime_error("could not open GAM for reading: " + gam_filename);
        }
        state.cursor.reset(new GAMIndex::cursor_t(*state.stream));
    }
    return state;
}

void IndexedGamSiteReadSource::fetch_span(const vector<pair<nid_t, nid_t>>& ranges,
                                          const function<void(Alignment&)>& iteratee) const {
    ThreadState& state = thread_state();
    // GAMIndex::find already takes a range list and de-duplicates across it, which is
    // exactly the contract fetch_span promises.
    vector<pair<id_t, id_t>> query;
    query.reserve(ranges.size());
    for (const auto& range : ranges) {
        query.emplace_back((id_t)range.first, (id_t)range.second);
    }
    index->find(*state.cursor, query, [&](const Alignment& aln) {
        if (!passes_filter(aln)) {
            return;
        }
        count_fetched();
        // GAMIndex hands out a const reference to its own buffer, so copy the read into a local
        // that the caller may move from.
        Alignment owned = aln;
        iteratee(owned);
    });
}

////////////////////////////////////////////////////////////////////////////////
// GafBaseSiteReadSource
////////////////////////////////////////////////////////////////////////////////

GafBaseSiteReadSource::GafBaseSiteReadSource(const HandleGraph& graph,
                                             const string& gaf_base_filename,
                                             const string& gbz_filename,
                                             const SiteReadFilter& filter,
                                             size_t window_size,
                                             size_t cache_entries,
                                             const string& binary)
    : WindowedSiteReadSource(filter, window_size, cache_entries),
      graph(graph),
      gaf_base_filename(gaf_base_filename),
      gbz_filename(gbz_filename),
      binary(binary) {

    threads.resize(max(1, get_thread_count()));
}

GafBaseSiteReadSource::~GafBaseSiteReadSource() {
    for (ThreadState& state : threads) {
        for (const string& path : state.gaf_paths) {
            temp_file::remove(path);
        }
    }
}

GafBaseSiteReadSource::ThreadState& GafBaseSiteReadSource::thread_state() const {
    int tid = omp_get_thread_num();
    if ((size_t)tid >= threads.size()) {
#pragma omp critical (gaf_base_threads)
        if ((size_t)tid >= threads.size()) {
            threads.resize(tid + 1);
        }
    }
    return threads[tid];
}

const string& GafBaseSiteReadSource::gaf_path(ThreadState& state, size_t slot) const {
    // One output file per in-flight query per thread, reused for every query that lands
    // in that slot. temp_file::create is mutex-guarded, so this is safe to race into.
    while (state.gaf_paths.size() <= slot) {
        state.gaf_paths.push_back(temp_file::create("vg-gafbase-reads-"));
    }
    return state.gaf_paths[slot];
}

size_t GafBaseSiteReadSource::argv_node_budget() {
    static const size_t budget = [](){
        long arg_max = sysconf(_SC_ARG_MAX);
        if (arg_max <= 0) {
            arg_max = 262144;   // conservative: POSIX only guarantees 4096
        }
        // "-n" and the id, each NUL-terminated. Node ids here run to ten digits.
        const size_t per_node = 3 + 11;
        const size_t from_limit = (size_t)arg_max / 4 / per_node;
        // At least 4096 and at most 65536 node IDs per child.
        return min<size_t>(max<size_t>(4096, from_limit), 65536);
    }();
    return budget;
}

GafBaseSiteReadSource::PendingQuery GafBaseSiteReadSource::spawn_query(
    ThreadState& state, size_t slot, const vector<nid_t>& nodes) const {

    // Build argv. --context 0 keeps the subgraph to the nodes we asked for, since the default
    // context would pull in reads no site here wants. --alignments overlapping returns each read
    // whole: the default, `clipped`, can cut one read into several pieces, which would put one
    // read in several rows of the likelihood matrix.
    vector<string> args{binary, "query", gbz_filename};
    args.reserve(args.size() + 2 * nodes.size() + 8);
    for (nid_t node : nodes) {
        args.push_back("-n");
        args.push_back(to_string(node));
    }
    args.push_back("--context");
    args.push_back("0");
    args.push_back("--gaf-base");
    args.push_back(gaf_base_filename);
    args.push_back("--gaf-output");
    args.push_back(gaf_path(state, slot));
    args.push_back("--alignments");
    args.push_back("overlapping");

    vector<const char*> argv;
    argv.reserve(args.size() + 1);
    for (const string& arg : args) {
        argv.push_back(arg.c_str());
    }
    argv.push_back(nullptr);

    // Capture stderr so a failure can say why, rather than just reporting a code.
    string err_path = gaf_path(state, slot) + ".err";

    // posix_spawn rather than fork and exec: fork() in a process whose threads are allocating
    // makes libc hold a lock around malloc, stalling the other threads. posix_spawn does not copy
    // the address space, so there is no such lock, and no stdio buffers need flushing.
    posix_spawn_file_actions_t actions;
    if (posix_spawn_file_actions_init(&actions) != 0) {
        throw runtime_error("posix_spawn_file_actions_init() failed: " + string(strerror(errno)));
    }
    // The subgraph GFA goes to stdout and we do not want it; only the separate
    // --gaf-output file interests us. stderr is captured so a failure can say why.
    posix_spawn_file_actions_addopen(&actions, STDOUT_FILENO, "/dev/null", O_WRONLY, 0);
    posix_spawn_file_actions_addopen(&actions, STDERR_FILENO, err_path.c_str(),
                                     O_WRONLY | O_CREAT | O_TRUNC, 0600);

    pid_t pid = 0;
    // posix_spawnp promises not to modify argv but cannot say so in C's type system;
    // see the same cast in index_registry.cpp's kmc call.
    int spawn_err = posix_spawnp(&pid, binary.c_str(), &actions, nullptr,
                                 (char* const*)&argv[0], environ);
    posix_spawn_file_actions_destroy(&actions);
    if (spawn_err != 0) {
        if (spawn_err == ENOENT) {
            throw runtime_error("could not execute '" + binary + "'. Install gbz-base "
                                "(https://github.com/jltsiren/gbz-base) and put it on your "
                                "PATH, or pass --gaf-base-binary with its location.");
        }
        throw runtime_error("posix_spawnp() failed for " + binary + ": " + strerror(spawn_err));
    }

    ++queries;
    PendingQuery pending;
    pending.pid = pid;
    pending.gaf_path = gaf_path(state, slot);
    pending.err_path = std::move(err_path);
    return pending;
}

size_t GafBaseSiteReadSource::reap_query(PendingQuery& pending,
                                         const function<void(Alignment&)>& iteratee) const {
    int child_stat = 0;
    while (waitpid(pending.pid, &child_stat, 0) == -1) {
        if (errno != EINTR) {
            throw runtime_error("waitpid() failed for " + binary + ": " + strerror(errno));
        }
    }
    const string& err_path = pending.err_path;

    int ret = WIFEXITED(child_stat) ? WEXITSTATUS(child_stat) : -1;
    if (ret != 0) {
        string message;
        {
            ifstream err_in(err_path);
            stringstream buffer;
            buffer << err_in.rdbuf();
            message = buffer.str();
        }
        unlink(err_path.c_str());
        // A missing binary is reported by posix_spawnp as ENOENT, so exit code 127 here comes from
        // gbz-base itself.
        throw runtime_error(binary + " query failed with exit code " + to_string(ret) +
                            (message.empty() ? "" : ": " + message));
    }
    unlink(err_path.c_str());

    // Parse the GAF text back into Alignments.
    size_t parsed = 0;
    vg::io::gaf_unpaired_for_each(graph, pending.gaf_path, [&](Alignment& aln) {
        ++parsed;
        if (!passes_filter(aln)) {
            return;
        }
        count_fetched();
        iteratee(aln);
    });
    return parsed;
}

size_t GafBaseSiteReadSource::run_query(ThreadState& state, const vector<nid_t>& nodes,
                                        const function<void(Alignment&)>& iteratee) const {
    if (nodes.empty()) {
        return 0;
    }
    PendingQuery pending = spawn_query(state, 0, nodes);
    return reap_query(pending, iteratee);
}

void GafBaseSiteReadSource::fetch_span(const vector<pair<nid_t, nid_t>>& ranges,
                                       const function<void(Alignment&)>& iteratee) const {
    // Ask only for node IDs that exist. The ranges are ranges of IDs, but ID space is
    // not dense, and gbz-base is entitled to complain about a node that is not there.
    vector<nid_t> nodes;
    for (const auto& range : ranges) {
        for (nid_t id = max<nid_t>(1, range.first); id <= range.second; ++id) {
            if (graph.has_node(id)) {
                nodes.push_back(id);
            }
        }
    }
    if (nodes.empty()) {
        return;
    }

    ThreadState& state = thread_state();

    // Drop reads returned more than once, whether or not the query was split: a repeated read
    // would be a second row for the same alignment in the likelihood matrix. Reads are keyed on
    // name and start position, not on name alone, since paired mates share a name.
    auto deduped = [&](const function<void(Alignment&)>& emit) {
        auto seen = make_shared<unordered_set<string>>();
        return [&emit, seen, this](Alignment& aln) {
            string key = aln.name();
            if (aln.path().mapping_size() > 0) {
                const Position& pos = aln.path().mapping(0).position();
                key += "\t" + to_string(pos.node_id()) + "\t" + to_string(pos.offset()) +
                       (pos.is_reverse() ? "-" : "+");
            }
            if (seen->insert(std::move(key)).second) {
                emit(aln);
            } else {
                ++duplicates_dropped;
            }
        };
    };

    if (nodes.size() <= max_query_nodes) {
        run_query_or_die(state, nodes, deduped(iteratee));
        return;
    }

    // Too many nodes for one command line, so split the query and run the pieces at the same
    // time; the work is mostly the children's, so waiting for each in turn would leave this thread
    // idle. Results are consumed in chunk order, so the order reads reach the caller does not
    // depend on which child finishes first. A read that spans a chunk boundary comes back from
    // both chunks, and is dropped the second time.
    try {
        vector<PendingQuery> pending;
        for (size_t start = 0, slot = 0; start < nodes.size();
             start += max_query_nodes, ++slot) {
            size_t end = min(start + max_query_nodes, nodes.size());
            vector<nid_t> chunk(nodes.begin() + start, nodes.begin() + end);
            pending.push_back(spawn_query(state, slot, chunk));
        }

        // Chunks overlap in reads rather than in nodes -- a read spanning a chunk boundary comes
        // back from both -- so one shared filter spans every chunk of this query.
        auto emit = deduped(iteratee);
        for (PendingQuery& query : pending) {
            reap_query(query, emit);
        }
    } catch (const std::exception& e) {
        // Calling happens inside an OpenMP parallel region, where an exception must not
        // propagate; see run_query_or_die.
        cerr << "error[vg::GafBaseSiteReadSource] " << e.what() << endl;
        exit(EXIT_FAILURE);
    }
}

void GafBaseSiteReadSource::run_query_or_die(ThreadState& state, const vector<nid_t>& nodes,
                                             const function<void(Alignment&)>& iteratee) const {
    // Calling runs inside an OpenMP parallel region, which an exception must not leave. A failed
    // query also means this site cannot be scored correctly, so report the error and exit.
    try {
        run_query(state, nodes, iteratee);
    } catch (const std::exception& e) {
        cerr << "error[vg::GafBaseSiteReadSource] " << e.what() << endl;
        exit(EXIT_FAILURE);
    }
}

void GafBaseSiteReadSource::check_setup() const {
    // Query the first node in the graph. A working setup returns cleanly; a missing
    // binary, unreadable database, or mismatched graph fails here rather than inside
    // a worker thread on the first snarl.
    nid_t probe = 0;
    graph.for_each_handle((function<bool(const handle_t&)>)[&](const handle_t& handle) -> bool {
        probe = graph.get_id(handle);
        return false;
    });
    if (probe == 0) {
        return;
    }

    ThreadState& state = thread_state();
    run_query(state, vector<nid_t>{probe}, [](const Alignment&) {});
}

size_t GafBaseSiteReadSource::get_query_count() const {
    return queries.load();
}

}
