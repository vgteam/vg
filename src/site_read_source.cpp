#include "site_read_source.hpp"

#include <algorithm>
#include <unordered_set>

#include <vg/io/alignment_io.hpp>
#include <vg/io/stream.hpp>

#include <omp.h>

#include "utility.hpp"

namespace vg {

using namespace std;

vector<Alignment> SiteReadSource::get_reads(const vector<pair<nid_t, nid_t>>& ranges) const {
    vector<Alignment> found;
    for_each_read(ranges, [&](const Alignment& aln) {
        found.push_back(aln);
    });
    return found;
}

void InMemorySiteReadSource::add_read(const Alignment& aln, const Filter& filter) {
    if (filter.skip_secondary && aln.is_secondary()) {
        ++filtered_count;
        return;
    }
    if (aln.mapping_quality() < filter.min_mapq) {
        ++filtered_count;
        return;
    }
    if (filter.skip_unmapped && aln.path().mapping_size() == 0) {
        ++filtered_count;
        return;
    }

    size_t read_index = reads.size();
    reads.push_back(aln);

    // Index under every node the read touches. A read can visit the same node
    // more than once (a cycle, or a snarl traversed twice), so only add the
    // index once per node or get_each_read would visit the read twice.
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
    // Deliberately serial rather than for_each_parallel: bucketing into a
    // shared map would need a mutex, and this runs once at startup where
    // simplicity is worth more than the wall clock.
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
                                          const function<void(const Alignment&)>& iteratee) const {

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
                    iteratee(reads[read_index]);
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
// IndexedGamSiteReadSource
////////////////////////////////////////////////////////////////////////////////

IndexedGamSiteReadSource::IndexedGamSiteReadSource(const string& gam_filename,
                                                   const string& index_filename,
                                                   const SiteReadFilter& filter,
                                                   size_t window_size,
                                                   size_t cache_entries)
    : gam_filename(gam_filename), filter(filter),
      window_size(max<size_t>(1, window_size)),
      cache_entries(max<size_t>(1, cache_entries)) {

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
        state.cache.resize(cache_entries);
    }
    return state;
}

size_t IndexedGamSiteReadSource::window_of(nid_t id) const {
    return (size_t)(id / (nid_t)window_size);
}

bool IndexedGamSiteReadSource::touches(const Alignment& aln,
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

void IndexedGamSiteReadSource::fetch_span(ThreadState& state, nid_t min_id, nid_t max_id,
                                         const function<void(const Alignment&)>& iteratee) const {
    vector<pair<id_t, id_t>> query{{(id_t)min_id, (id_t)max_id}};
    index->find(*state.cursor, query, [&](const Alignment& aln) {
        if (filter.skip_secondary && aln.is_secondary()) {
            ++filtered;
            return;
        }
        if (aln.mapping_quality() < filter.min_mapq) {
            ++filtered;
            return;
        }
        if (filter.skip_unmapped && aln.path().mapping_size() == 0) {
            ++filtered;
            return;
        }
        ++fetched;
        iteratee(aln);
    });
}

void IndexedGamSiteReadSource::fetch_window(ThreadState& state, size_t window,
                                           vector<Alignment>& out) const {
    nid_t lo = (nid_t)(window * window_size);
    nid_t hi = lo + (nid_t)window_size - 1;
    fetch_span(state, lo, hi, [&](const Alignment& aln) {
        out.push_back(aln);
    });
}

void IndexedGamSiteReadSource::for_each_read(
    const vector<pair<nid_t, nid_t>>& ranges,
    const function<void(const Alignment&)>& iteratee) const {

    if (ranges.empty()) {
        return;
    }

    ThreadState& state = thread_state();

    nid_t min_id = ranges.front().first;
    nid_t max_id = ranges.front().second;
    for (const auto& range : ranges) {
        min_id = min(min_id, range.first);
        max_id = max(max_id, range.second);
    }

    size_t first_window = window_of(min_id);
    size_t last_window = window_of(max_id);

    if (first_window != last_window) {
        // The site straddles a window boundary. Fetch the span directly rather than
        // stitching windows together: index->find already emits each read at most
        // once, so this sidesteps de-duplicating reads that span the boundary. A
        // site large enough to do this is rare, and caching a region of unbounded
        // size would defeat the point of windowing.
        ++cache_misses;
        fetch_span(state, min_id, max_id, [&](const Alignment& aln) {
            if (touches(aln, ranges)) {
                iteratee(aln);
            }
        });
        return;
    }

    // Serve from the cache if this window is already resident. With the caller
    // visiting sites in node-ID order (GraphCaller::set_node_id_ordering) this is the
    // common case, and each window is fetched exactly once.
    for (const CacheEntry& entry : state.cache) {
        if (entry.valid && entry.window == first_window) {
            ++cache_hits;
            for (const Alignment& aln : entry.reads) {
                if (touches(aln, ranges)) {
                    iteratee(aln);
                }
            }
            return;
        }
    }
    ++cache_misses;

    CacheEntry entry;
    entry.window = first_window;
    fetch_window(state, first_window, entry.reads);
    entry.valid = true;

    for (const Alignment& aln : entry.reads) {
        if (touches(aln, ranges)) {
            iteratee(aln);
        }
    }

    state.cache[state.next_evict] = std::move(entry);
    state.next_evict = (state.next_evict + 1) % state.cache.size();
}

size_t IndexedGamSiteReadSource::get_read_count() const {
    return fetched.load();
}

size_t IndexedGamSiteReadSource::get_filtered_count() const {
    return filtered.load();
}

size_t IndexedGamSiteReadSource::get_cache_hits() const {
    return cache_hits.load();
}

size_t IndexedGamSiteReadSource::get_cache_misses() const {
    return cache_misses.load();
}

}
