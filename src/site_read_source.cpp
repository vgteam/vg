#include "site_read_source.hpp"

#include <algorithm>
#include <cassert>
#include <cctype>
#include <cerrno>
#include <chrono>
#include <cstring>
#include <limits>
#include <sstream>
#include <unordered_set>

#include <cstdlib>
#include <fcntl.h>
#include <spawn.h>
#include <sys/resource.h>
#include <sys/stat.h>
#include <sys/time.h>
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

namespace {

/// One hold of `SiteReadSource::hold_reads`: the reads the fetch of a site too wide for one window
/// kept for the sites below it; the nodes a query on which it answers, as merged ranges, which with
/// off-panel detection are those of the candidates `bar` admits; the nodes the fetch covered, on
/// which a query too wide for one window is answered too, if the hold is exact; and the read starts
/// on every node the fetch covered, tallied as a window's are, by node, then mapping quality.
struct ReadHold {
    struct Start {
        nid_t node = 0;
        int32_t mapq = 0;
        size_t reads = 0;
        size_t bases = 0;
    };
    const void* source = nullptr;
    EditCallBar bar;
    bool filled = false;
    /// It keeps every read its fetch delivered, rather than only those on `nodes`: a fetch of the
    /// site's own ranges keeps them all, for the sites below, which would fetch them again, while
    /// one from cached windows leaves those sites to find the windows cached.
    bool all_reads = false;
    /// The fetch dropped no read for having the name and start of a read before it, so it kept
    /// every read on the nodes it covered, as a fetch of some of them would.
    bool exact = false;
    vector<pair<nid_t, nid_t>> nodes;
    vector<pair<nid_t, nid_t>> fetched;
    vector<Start> starts;

    /// The reads, parsed, in fetch order, with every mapping indexed by node and the read offset
    /// before each, as a cached window's are: a site below asks for a few of thousands of held
    /// reads, and is answered by binary search with its reads indexed, so that its work grows
    /// with its own size rather than with the held reads'.
    struct Visit {
        nid_t node = 0;
        uint32_t read = 0;
        uint32_t mapping = 0;
        bool operator<(const Visit& other) const {
            if (node != other.node) return node < other.node;
            if (read != other.read) return read < other.read;
            return mapping < other.mapping;
        }
    };
    vector<Alignment> alignments;
    vector<Visit> visits;
    vector<uint32_t> offsets;
    vector<uint32_t> offset_start;
};

/// This thread's holds, innermost last.
thread_local vector<ReadHold> read_holds;

/// Are all the nodes of the first merged ranges among those of the second?
bool ranges_within(const vector<pair<nid_t, nid_t>>& inner, const vector<pair<nid_t, nid_t>>& outer) {
    return std::all_of(inner.begin(), inner.end(), [&](const pair<nid_t, nid_t>& r) {
        auto it = std::upper_bound(outer.begin(), outer.end(),
                                   make_pair(r.first, numeric_limits<nid_t>::max()));
        return it != outer.begin() && std::prev(it)->first <= r.first
               && r.second <= std::prev(it)->second;
    });
}

/// What a fetch tells reads apart by, since paired mates share a name: the name and where the
/// read starts. A fetch passes on the first read with each, in file order.
string name_and_start(const Alignment& aln) {
    string key = aln.name();
    if (aln.path().mapping_size() > 0) {
        const Position& pos = aln.path().mapping(0).position();
        key += "\t" + to_string(pos.node_id()) + "\t" + to_string(pos.offset()) +
               (pos.is_reverse() ? "-" : "+");
    }
    return key;
}

/// For a search of sorted candidates by node.
bool candidate_node_below(const EditCandidate& c, nid_t node) {
    return c.node < node;
}

} // namespace


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

void InMemorySiteReadSource::count_edits(const HandleGraph* graph, const EditCountParams& params) {
    edit_graph = graph;
    edit_params = params;
}

const vector<EditCandidate>& InMemorySiteReadSource::all_edits() const {
    std::call_once(edits_counted, [&]() {
        if (edit_graph == nullptr) {
            return;
        }
        edits = vg::count_edits(
            *edit_graph,
            [&](const function<void(const Alignment&)>& visit) {
                for (const Alignment& aln : reads) {
                    visit(aln);
                }
            },
            [](nid_t) { return true; }, edit_params);
    });
    return edits;
}

void InMemorySiteReadSource::for_each_edit_candidate(
    const vector<nid_t>& nodes, const function<void(const EditCandidate&)>& iteratee) const {
    for_each_candidate_on(all_edits(), nodes, iteratee);
}

vector<EditCandidate> InMemorySiteReadSource::counted_edit_candidates() const {
    return all_edits();
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

vector<pair<nid_t, nid_t>> WindowedSiteReadSource::merge_ranges(vector<pair<nid_t, nid_t>> ranges) {
    std::sort(ranges.begin(), ranges.end());
    size_t kept = 0;
    for (size_t i = 0; i < ranges.size(); ++i) {
        if (kept > 0 && ranges[i].first <= ranges[kept - 1].second + 1) {
            ranges[kept - 1].second = max(ranges[kept - 1].second, ranges[i].second);
        } else {
            ranges[kept++] = ranges[i];
        }
    }
    ranges.resize(kept);
    return ranges;
}

bool WindowedSiteReadSource::in_ranges(nid_t node_id, const vector<pair<nid_t, nid_t>>& ranges) {
    // The last range starting at or before the ID is the only one that can hold it, since merged
    // ranges neither overlap nor touch.
    auto after = std::upper_bound(ranges.begin(), ranges.end(), node_id,
                                  [](nid_t id, const pair<nid_t, nid_t>& range) {
                                      return id < range.first;
                                  });
    return after != ranges.begin() && node_id <= (after - 1)->second;
}

bool WindowedSiteReadSource::touches(const Alignment& aln,
                                    const vector<pair<nid_t, nid_t>>& ranges) {
    for (const auto& mapping : aln.path().mapping()) {
        if (in_ranges(mapping.position().node_id(), ranges)) {
            return true;
        }
    }
    return false;
}

/// A query whose nodes lie in at most this many windows is answered from those windows, fetched
/// whole and cached, rather than fetched uncached. Neighbouring sites need the same windows, so
/// their reads are read once rather than once per site. A query over more windows is fetched by
/// its own node ranges, which reads less than its windows would.
static constexpr size_t WIDE_QUERY_MAX_WINDOWS = 16;

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

    const bool wide = first_window != last_window;
    if (!read_holds.empty() && deliver_held(ranges, wide, iteratee)) {
        return;
    }

    if (wide) {
        // The site crosses a window boundary, so it gets its reads by a fetch of its own. If its
        // nodes lie in a few windows, the fetch takes them from those windows, fetched whole and
        // cached as for a site inside one window, so that neighbouring sites read them once.
        // Otherwise its ranges are fetched directly, since a region of unbounded size is not worth
        // caching; they are passed as they are, rather than collapsed to one span, since a snarl's
        // nodes can be spread thinly over a wide span of IDs. Either way the site gets the same
        // reads in the same order.
        vector<pair<nid_t, nid_t>> merged = merge_ranges(ranges);
        vector<size_t> windows;
        const bool few = few_windows(merged, windows);
        // Off-panel detection counts the reads on the site's nodes that no fetch has counted,
        // which are kept for it once delivered.
        vector<pair<nid_t, nid_t>> to_count;
        // The innermost hold on this thread, if this is the first such fetch within it. Each read
        // keeps its place in the fetch, so that the held reads come in the order they came here.
        // With off-panel detection it answers any query on the nodes of the candidates its bar
        // admits, which alone become edit sites: those counted before, known now, and those this
        // fetch counts, known once it has. A fetch of the site's own ranges keeps every read, for
        // the sites below as wide, which would fetch them again; one from cached windows keeps only
        // the reads on those nodes, or without off-panel detection fills no hold, since the sites
        // below find the windows cached.
        ReadHold* hold = !read_holds.empty() && read_holds.back().source == this
                                 && !read_holds.back().filled
                             ? &read_holds.back()
                             : nullptr;
        ReadHold* filling = few && edit_graph == nullptr ? nullptr : hold;
        vector<pair<nid_t, nid_t>> known;
        if (edit_graph != nullptr) {
            lock_guard<std::mutex> guard(cache_mutex);
            to_count = uncounted(merged);
            if (hold != nullptr) {
                known = candidate_nodes(merged, hold->bar);
            }
        }
        bool keep_all = false;
        size_t n_scanned = 0, n_delivered = 0;
        vector<Alignment> kept;
        vector<size_t> kept_at;
        vector<pair<size_t, Alignment>> held;
        vector<ReadHold::Start> starts;
        // Each read the site gets, with `owned` the same read where it may be moved from.
        auto take = [&](const Alignment& aln, Alignment* owned) {
            ++n_delivered;
            if (filling != nullptr) {
                // The filter drops reads with no mappings.
                const nid_t first = aln.path().mapping(0).position().node_id();
                if (in_ranges(first, merged)) {
                    starts.push_back({first, aln.mapping_quality(), 1, aln.sequence().size()});
                }
            }
            // Not indexed: this path serves the few sites too wide for one window, and an index for
            // a read seen once would cost as much as the walk it saves.
            SiteRead read;
            read.aln = &aln;
            iteratee(read);
            const bool count = !to_count.empty() && touches(aln, to_count);
            const bool keep =
                filling != nullptr && (keep_all || (!known.empty() && touches(aln, known)));
            if (!count && !keep) {
                return;
            }
            // Moved where the fetch's own copy may be taken, copied from a cached window.
            Alignment taken = owned != nullptr ? Alignment(std::move(*owned)) : Alignment(aln);
            if (count) {
                kept_at.push_back(n_delivered);
                kept.push_back(std::move(taken));
            } else {
                held.emplace_back(n_delivered, std::move(taken));
            }
        };
        // The windows are fetched for this site alone, so they are not counted for off-panel
        // candidates, nor their read starts tallied, until another query uses them.
        size_t dropped = 0;
        const bool from_windows =
            few && for_each_window_read(merged, windows, false,
                                        [&](const Alignment& aln) {
                                            ++n_scanned;
                                            take(aln, nullptr);
                                        },
                                        dropped);
        if (from_windows) {
            ++wide_from_windows;
        } else {
            filling = hold;
            keep_all = true;
            ++cache_misses;
            ++straddles;
            straddle_nodes += (size_t)(max_id - min_id + 1);
            size_t wanted = 0;
            for (const auto& range : ranges) {
                wanted += (size_t)(range.second - range.first + 1);
            }
            straddle_wanted += wanted;
            dropped = fetch_span_keyed(ranges, [&](Alignment& aln, uint64_t) {
                ++n_scanned;
                if (touches(aln, merged)) {
                    take(aln, &aln);
                }
            });
        }
        scanned += n_scanned;
        delivered += n_delivered;
        if (!to_count.empty()) {
            vector<EditCandidate> counted = vg::count_edits(
                *edit_graph,
                [&](const function<void(const Alignment&)>& visit) {
                    for (const Alignment& aln : kept) {
                        visit(aln);
                    }
                },
                [&](nid_t node) { return in_ranges(node, to_count); }, edit_params);
            if (filling != nullptr) {
                // The nodes of the candidates just counted, on the site, join those known.
                for (const EditCandidate& c : counted) {
                    if (!filling->bar.admits(c)) {
                        continue;
                    }
                    known.emplace_back(c.node, c.node);
                    for (nid_t node : c.walk) {
                        known.emplace_back(node, node);
                    }
                }
                known.erase(std::remove_if(known.begin(), known.end(),
                                           [&](const pair<nid_t, nid_t>& r) {
                                               return r.first == r.second
                                                      && !in_ranges(r.first, merged);
                                           }),
                            known.end());
                known = merge_ranges(std::move(known));
                // Each list is in fetch order, so the two are merged rather than sorted.
                const size_t first_kept = held.size();
                for (size_t i = 0; i < kept.size(); ++i) {
                    if (keep_all || touches(kept[i], known)) {
                        held.emplace_back(kept_at[i], std::move(kept[i]));
                    }
                }
                std::inplace_merge(held.begin(), held.begin() + first_kept, held.end(),
                                   [](const pair<size_t, Alignment>& a,
                                      const pair<size_t, Alignment>& b) {
                                       return a.first < b.first;
                                   });
            }
            kept.clear();
            lock_guard<std::mutex> guard(cache_mutex);
            // Another thread may have counted some of the nodes meanwhile, alike.
            to_count.insert(to_count.end(), span_counted.begin(), span_counted.end());
            span_counted = merge_ranges(std::move(to_count));
            span_edits.insert(span_edits.end(), counted.begin(), counted.end());
            std::sort(span_edits.begin(), span_edits.end());
            span_edits.erase(std::unique(span_edits.begin(), span_edits.end(),
                                         [](const EditCandidate& a, const EditCandidate& b) {
                                             return !(a < b) && !(b < a);
                                         }),
                             span_edits.end());
        }
        if (filling != nullptr) {
            filling->filled = true;
            filling->all_reads = keep_all;
            // A fetch that dropped a read for repeating an earlier read's name and start can lack a
            // read that a fetch of fewer of its nodes, which the earlier read misses, delivers.
            filling->exact = has_record_keys() && dropped == 0;
            if (edit_graph != nullptr) {
                filling->nodes = std::move(known);
                filling->fetched = merged;
            } else if (filling->exact) {
                filling->fetched = merged;
            }
            std::sort(starts.begin(), starts.end(),
                      [](const ReadHold::Start& a, const ReadHold::Start& b) {
                          return a.node != b.node ? a.node < b.node : a.mapq < b.mapq;
                      });
            for (const ReadHold::Start& start : starts) {
                if (!filling->starts.empty() && filling->starts.back().node == start.node
                    && filling->starts.back().mapq == start.mapq) {
                    filling->starts.back().reads += start.reads;
                    filling->starts.back().bases += start.bases;
                } else {
                    filling->starts.push_back(start);
                }
            }
            filling->alignments.reserve(held.size());
            for (auto& entry : held) {
                const uint32_t r = (uint32_t)filling->alignments.size();
                filling->offset_start.push_back((uint32_t)filling->offsets.size());
                uint32_t offset = 0;
                const Path& path = entry.second.path();
                for (int64_t i = 0; i < path.mapping_size(); ++i) {
                    filling->visits.push_back(
                        ReadHold::Visit{path.mapping(i).position().node_id(), r, (uint32_t)i});
                    filling->offsets.push_back(offset);
                    offset += (uint32_t)mapping_to_length(path.mapping(i));
                }
                filling->offsets.push_back(offset);
                filling->alignments.push_back(std::move(entry.second));
            }
            filling->offset_start.push_back((uint32_t)filling->offsets.size());
            std::sort(filling->visits.begin(), filling->visits.end());
        }
        return;
    }

    // Usually served from the cache when sites are visited in node-ID order
    // (GraphCaller::set_snarl_batching).
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
WindowedSiteReadSource::get_window(size_t window, bool& was_fetched, bool counted) const {
    was_fetched = false;
    // A cached window that only sites too wide for one window have used, which leave it untallied
    // and uncounted: this use tallies and counts it, as fetching it would have.
    shared_ptr<const CacheEntry> found_entry;
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
                if (!counted || (starts.count(window) > 0
                                 && (edit_graph == nullptr || edits.count(window) > 0))) {
                    return found->second.entry;
                }
                found_entry = found->second.entry;
                break;
            }
            // Another thread is fetching it.
            cache_filled.wait(lock);
        }
    }

    shared_ptr<CacheEntry> entry;
    bool tallied = false;
    bool edits_counted = true;
    if (found_entry) {
        lock_guard<std::mutex> guard(cache_mutex);
        tallied = starts.count(window) > 0;
        edits_counted = edit_graph == nullptr || edits.count(window) > 0;
    } else {
        try {
            entry = make_shared<CacheEntry>(load_window(window));
            {
                lock_guard<std::mutex> guard(cache_mutex);
                tallied = starts.count(window) > 0;
                edits_counted = edit_graph == nullptr || edits.count(window) > 0;
            }
        } catch (...) {
            // Release the claim, so that a thread waiting on this window fetches it itself rather
            // than waiting for ever.
            lock_guard<std::mutex> guard(cache_mutex);
            cache.erase(window);
            cache_filled.notify_all();
            throw;
        }
        if (tallied) {
            ++refetches;
        }
    }
    const CacheEntry& reads = found_entry ? *found_entry : *entry;
    // Only this thread can be fetching the window now, so no other thread can tally it first, but
    // threads finding it cached untallied may tally it alike.
    vector<StartTally> tallies;
    const bool tally = counted && !tallied;
    if (tally) {
        tallies = tally_starts(reads);
    }
    // Likewise the off-panel candidates, on the window's own nodes only, since a read that reaches
    // into another window is fetched with that window too.
    vector<EditCandidate> window_candidates;
    const bool count = counted && !edits_counted;
    if (count) {
        window_candidates = vg::count_edits(
            *edit_graph,
            [&](const function<void(const Alignment&)>& visit) {
                for (const Alignment& aln : reads.reads) {
                    visit(aln);
                }
            },
            [&](nid_t node) { return window_of(node) == window; }, edit_params);
    }
    if (found_entry) {
        lock_guard<std::mutex> guard(cache_mutex);
        if (tally) {
            starts.emplace(window, std::move(tallies));
        }
        if (count) {
            edits.emplace(window, std::move(window_candidates));
        }
        return found_entry;
    }
    was_fetched = true;

    // Windows dropped from the cache are freed here, after the lock is released: freeing a
    // window's reads takes long enough that other threads would queue on the lock.
    vector<shared_ptr<const CacheEntry>> dropped;
    lock_guard<std::mutex> guard(cache_mutex);
    if (tally) {
        starts.emplace(window, std::move(tallies));
    }
    if (count) {
        edits.emplace(window, std::move(window_candidates));
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
    entry.dropped_repeats = fetch_span_keyed({{lo, hi}}, [&](Alignment& aln, uint64_t key) {
        // Take ownership rather than copy the alignment. fetch_span hands out a mutable reference
        // so that this can move from it; the backends reuse one Alignment per record and clear it
        // before the next.
        index_read(aln, (uint32_t)entry.reads.size(), entry);
        entry.reads.push_back(std::move(aln));
        entry.keys.push_back(key);
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

void WindowedSiteReadSource::count_edits(const HandleGraph* graph,
                                         const EditCountParams& params) {
    edit_graph = graph;
    edit_params = params;
}

const vector<EditCandidate>& WindowedSiteReadSource::window_edits(size_t window) const {
    {
        lock_guard<std::mutex> guard(cache_mutex);
        auto found = edits.find(window);
        if (found != edits.end()) {
            return found->second;
        }
    }
    // Never fetched, or fetched before counting was turned on. A window is counted before its
    // fetch is published, so once get_window returns it is in, unless the fetch came first.
    bool was_fetched = false;
    shared_ptr<const CacheEntry> entry = get_window(window, was_fetched);
    if (was_fetched) {
        ++edit_fetches;
    }
    {
        lock_guard<std::mutex> guard(cache_mutex);
        auto found = edits.find(window);
        if (found != edits.end()) {
            return found->second;
        }
    }
    vector<EditCandidate> counted = vg::count_edits(
        *edit_graph,
        [&](const function<void(const Alignment&)>& visit) {
            for (const Alignment& aln : entry->reads) {
                visit(aln);
            }
        },
        [&](nid_t node) { return window_of(node) == window; }, edit_params);
    lock_guard<std::mutex> guard(cache_mutex);
    return edits.emplace(window, std::move(counted)).first->second;
}

void WindowedSiteReadSource::hold_reads(const EditCallBar& bar) const {
    read_holds.emplace_back();
    read_holds.back().source = this;
    read_holds.back().bar = bar;
}

void WindowedSiteReadSource::release_reads() const {
    assert(!read_holds.empty() && read_holds.back().source == this);
    read_holds.pop_back();
}

bool WindowedSiteReadSource::deliver_held(const vector<pair<nid_t, nid_t>>& ranges, bool wide,
                                          const function<void(const SiteRead&)>& iteratee) const {
    if (std::none_of(read_holds.begin(), read_holds.end(),
                     [&](const ReadHold& hold) { return hold.source == this && hold.filled; })) {
        return false;
    }
    const vector<pair<nid_t, nid_t>> merged = merge_ranges(ranges);
    for (auto hold = read_holds.rbegin(); hold != read_holds.rend(); ++hold) {
        if (hold->source != this || !hold->filled) {
            continue;
        }
        // Every queried node lies on one of the hold's merged ranges, so every read on it was
        // kept. A query too wide for one window would get its reads by a fetch like the hold's,
        // and gets the same reads from the hold if the hold's fetch covered its nodes and kept
        // every read on them.
        if (!ranges_within(merged, hold->nodes)
            && !(wide && hold->all_reads && hold->exact && ranges_within(merged, hold->fetched))) {
            continue;
        }
        // As a cached window delivers: the held mappings on the queried nodes, by read, in fetch
        // order, each read with the index of its mappings in the ranges.
        vector<pair<uint32_t, uint32_t>> hits;
        for (const auto& range : merged) {
            ReadHold::Visit probe;
            probe.node = range.first;
            auto it = std::lower_bound(hold->visits.begin(), hold->visits.end(), probe);
            for (; it != hold->visits.end() && it->node <= range.second; ++it) {
                hits.emplace_back(it->read, it->mapping);
            }
        }
        std::sort(hits.begin(), hits.end());
        vector<uint32_t> mappings;
        mappings.reserve(hits.size());
        for (const auto& hit : hits) {
            mappings.push_back(hit.second);
        }
        // The hold's buffers, taken before handing out reads: the caller's code runs in between.
        const Alignment* alignments = hold->alignments.data();
        const uint32_t* offsets = hold->offsets.data();
        const uint32_t* offset_start = hold->offset_start.data();
        size_t n_delivered = 0;
        size_t i = 0;
        while (i < hits.size()) {
            size_t j = i;
            while (j < hits.size() && hits[j].first == hits[i].first) {
                ++j;
            }
            SiteRead read;
            read.aln = &alignments[hits[i].first];
            read.mappings = &mappings[i];
            read.mapping_count = j - i;
            read.read_offsets = &offsets[offset_start[hits[i].first]];
            iteratee(read);
            ++n_delivered;
            i = j;
        }
        ++held_queries;
        scanned += hits.size();
        delivered += n_delivered;
        return true;
    }
    return false;
}

vector<pair<nid_t, nid_t>> WindowedSiteReadSource::candidate_nodes(
    const vector<pair<nid_t, nid_t>>& merged, const EditCallBar& bar) const {
    vector<pair<nid_t, nid_t>> found;
    auto add = [&](const EditCandidate& c) {
        if (!bar.admits(c)) {
            return;
        }
        if (in_ranges(c.node, merged)) {
            found.emplace_back(c.node, c.node);
        }
        for (nid_t node : c.walk) {
            if (in_ranges(node, merged)) {
                found.emplace_back(node, node);
            }
        }
    };
    for (const auto& range : merged) {
        for (auto it = std::lower_bound(span_edits.begin(), span_edits.end(), range.first,
                                        candidate_node_below);
             it != span_edits.end() && it->node <= range.second; ++it) {
            add(*it);
        }
        for (size_t w = window_of(range.first); w <= window_of(range.second); ++w) {
            auto counted = edits.find(w);
            if (counted == edits.end()) {
                continue;
            }
            for (const EditCandidate& c : counted->second) {
                if (c.node >= range.first && c.node <= range.second) {
                    add(c);
                }
            }
        }
    }
    return merge_ranges(std::move(found));
}

vector<pair<nid_t, nid_t>>
WindowedSiteReadSource::uncounted(const vector<pair<nid_t, nid_t>>& merged) const {
    vector<pair<nid_t, nid_t>> left;
    for (const auto& range : merged) {
        for (nid_t id = range.first; id <= range.second; ++id) {
            if (edits.count(window_of(id)) || in_ranges(id, span_counted)) {
                continue;
            }
            if (!left.empty() && left.back().second + 1 == id) {
                left.back().second = id;
            } else {
                left.emplace_back(id, id);
            }
        }
    }
    return left;
}

void WindowedSiteReadSource::for_each_edit_candidate(
    const vector<nid_t>& nodes, const function<void(const EditCandidate&)>& iteratee) const {
    if (edit_graph == nullptr) {
        return;
    }
    size_t i = 0;
    while (i < nodes.size()) {
        // The nodes are sorted, so each window's nodes are consecutive.
        const size_t window = window_of(nodes[i]);
        size_t j = i;
        while (j < nodes.size() && window_of(nodes[j]) == window) {
            ++j;
        }
        const vector<nid_t> in_window(nodes.begin() + i, nodes.begin() + j);
        i = j;
        // From a wide site's fetch where that counted them and the window was never counted, so
        // that looking them up fetches nothing.
        vector<EditCandidate> spanned;
        bool from_span = false;
        {
            lock_guard<std::mutex> guard(cache_mutex);
            if (!edits.count(window)) {
                from_span = std::all_of(in_window.begin(), in_window.end(), [&](nid_t node) {
                    return in_ranges(node, span_counted);
                });
                if (from_span) {
                    // By node alone: a probe compared field by field would sort after the node's
                    // candidates at offset 0 of a lower kind or ALT, and skip them.
                    for (auto it = std::lower_bound(span_edits.begin(), span_edits.end(),
                                                    in_window.front(), candidate_node_below);
                         it != span_edits.end() && it->node <= in_window.back(); ++it) {
                        spanned.push_back(*it);
                    }
                }
            }
        }
        for_each_candidate_on(from_span ? spanned : window_edits(window), in_window, iteratee);
    }
}

bool WindowedSiteReadSource::fetched_together(nid_t a, nid_t b) const {
    return window_of(a) == window_of(b);
}

vector<EditCandidate> WindowedSiteReadSource::counted_edit_candidates() const {
    vector<EditCandidate> all;
    lock_guard<std::mutex> guard(cache_mutex);
    for (const auto& window : edits) {
        all.insert(all.end(), window.second.begin(), window.second.end());
    }
    // A node counted by a wide site's fetch and by its window's has the same candidates in both.
    all.insert(all.end(), span_edits.begin(), span_edits.end());
    std::sort(all.begin(), all.end());
    all.erase(std::unique(all.begin(), all.end(),
                          [](const EditCandidate& a, const EditCandidate& b) {
                              return !(a < b) && !(b < a);
                          }),
              all.end());
    return all;
}

void WindowedSiteReadSource::for_each_read_start(
    const vector<nid_t>& nodes,
    const function<void(int32_t mapq, size_t reads, size_t bases)>& iteratee) const {

    // From a hold on this thread whose fetch covered every node, which tallied their starts as
    // their windows would. Only with off-panel detection, whose holds always answered these: a
    // window's tallies can differ from a hold's where a fetch drops a read for repeating an earlier
    // read's name and start.
    for (auto hold = read_holds.rbegin(); edit_graph != nullptr && hold != read_holds.rend(); ++hold) {
        if (hold->source != this || !hold->filled
            || !std::all_of(nodes.begin(), nodes.end(),
                            [&](nid_t node) { return in_ranges(node, hold->fetched); })) {
            continue;
        }
        for (nid_t node : nodes) {
            auto it = std::lower_bound(hold->starts.begin(), hold->starts.end(), node,
                                       [](const ReadHold::Start& a, nid_t n) { return a.node < n; });
            for (; it != hold->starts.end() && it->node == node; ++it) {
                iteratee(it->mapq, it->reads, it->bases);
            }
        }
        ++held_queries;
        return;
    }

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

size_t WindowedSiteReadSource::drop_cached_windows() {
    // Moved out under the lock and freed after it, as get_window frees what it evicts.
    unordered_map<size_t, CacheSlot> dropped;
    {
        lock_guard<std::mutex> guard(cache_mutex);
        dropped.swap(cache);
    }
    size_t reads = 0;
    vector<shared_ptr<const CacheEntry>> entries;
    entries.reserve(dropped.size());
    for (auto& slot : dropped) {
        if (slot.second.entry) {
            reads += slot.second.entry->reads.size();
            entries.push_back(std::move(slot.second.entry));
        }
    }
    dropped.clear();
    // The windows are freed on several threads: on a whole genome they hold tens of gigabytes of
    // alignments, which are slow to free on one thread.
#pragma omp parallel for schedule(dynamic, 1)
    for (size_t i = 0; i < entries.size(); ++i) {
        entries[i].reset();
    }
    return reads;
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

bool WindowedSiteReadSource::few_windows(const vector<pair<nid_t, nid_t>>& merged,
                                         vector<size_t>& windows) const {
    // The windows' reads must be put back in file order, as a fetch of the ranges delivers them,
    // so the backend must key them.
    windows.clear();
    if (!has_record_keys()) {
        return false;
    }
    // Ascending, since the ranges are sorted and apart.
    for (const auto& range : merged) {
        size_t last = window_of(range.second);
        size_t w = window_of(range.first);
        if (!windows.empty() && w <= windows.back()) {
            w = windows.back() + 1;
        }
        for (; w <= last; ++w) {
            if (windows.size() == WIDE_QUERY_MAX_WINDOWS) {
                windows.clear();
                return false;
            }
            windows.push_back(w);
        }
    }
    return true;
}

bool WindowedSiteReadSource::for_each_window_read(
    const vector<pair<nid_t, nid_t>>& merged, const vector<size_t>& windows, bool counted,
    const function<void(const Alignment&)>& iteratee, size_t& dropped) const {

    // Held until the reads are delivered, even if the cache drops them meanwhile.
    vector<shared_ptr<const CacheEntry>> entries;
    entries.reserve(windows.size());
    for (size_t w : windows) {
        bool was_fetched = false;
        entries.push_back(get_window(w, was_fetched, counted));
        if (entries.back()->dropped_repeats > 0) {
            // A read this window's fetch dropped, for repeating an earlier read that touches the
            // window but not the ranges, is one that a fetch of the ranges delivers.
            return false;
        }
    }

    // Every read touching the ranges, from each window's node index. A read in several of the
    // windows is the same record, with the same key, in each.
    struct Hit {
        uint64_t key;
        uint32_t entry;
        uint32_t read;
    };
    vector<Hit> hits;
    vector<uint32_t> found;
    for (uint32_t e = 0; e < (uint32_t)entries.size(); ++e) {
        const CacheEntry& entry = *entries[e];
        const nid_t lo = (nid_t)(entry.window * window_size);
        const nid_t hi = lo + (nid_t)window_size - 1;
        found.clear();
        for (const auto& range : merged) {
            if (range.second < lo || range.first > hi) {
                continue;
            }
            IndexEntry probe{range.first, 0, 0};
            auto it = std::lower_bound(entry.node_index.begin(), entry.node_index.end(), probe);
            for (; it != entry.node_index.end() && it->node <= range.second; ++it) {
                found.push_back(it->read);
            }
        }
        std::sort(found.begin(), found.end());
        found.erase(std::unique(found.begin(), found.end()), found.end());
        for (uint32_t r : found) {
            hits.push_back(Hit{entry.keys[r], e, r});
        }
    }
    std::sort(hits.begin(), hits.end(), [](const Hit& a, const Hit& b) {
        return a.key != b.key ? a.key < b.key : a.entry < b.entry;
    });

    // As an uncached fetch delivers them: in file order, once each, and only the first of reads
    // sharing a name and a start.
    unordered_set<string> seen;
    dropped = 0;
    for (size_t i = 0; i < hits.size(); ++i) {
        if (i > 0 && hits[i].key == hits[i - 1].key) {
            continue;
        }
        const Alignment& aln = entries[hits[i].entry]->reads[hits[i].read];
        if (!seen.insert(name_and_start(aln)).second) {
            ++dropped;
            continue;
        }
        iteratee(aln);
    }
    return true;
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
      binary(binary),
      gbz_argument(gbz_filename),
      gaf_base_argument(gaf_base_filename) {

    threads.resize(max(1, get_thread_count()));

    const char* locking = getenv("VG_GAFBASE_LOCKING");
    immutable = !(locking != nullptr && string(locking) == "1");
    if (immutable) {
        gbz_argument = immutable_uri(gbz_filename);
        gaf_base_argument = immutable_uri(gaf_base_filename);
    }

    const char* log_path = getenv("VG_GAFBASE_QUERY_LOG");
    if (log_path != nullptr && *log_path != '\0') {
        query_log = fopen(log_path, "w");
        if (query_log == nullptr) {
            throw runtime_error("cannot write VG_GAFBASE_QUERY_LOG " + string(log_path) + ": " +
                                strerror(errno));
        }
        // Line-buffered, so the log is complete up to the last finished query however vg ends.
        setvbuf(query_log, nullptr, _IOLBF, 0);
        fprintf(query_log, "#start_epoch\twait_s\tuser_s\tsystem_s\tmax_rss_kb\texit\tn_nodes\t"
                           "first_node\tlast_node\tgaf_bytes\tnodes_hash\tmode\tparse_s\t"
                           "reads_parsed\tthread\n");
    }
}

GafBaseSiteReadSource::QueryTotals GafBaseSiteReadSource::get_query_totals() const {
    QueryTotals totals;
    totals.queries = queries.load();
    totals.wait_s = wait_us.load() / 1e6;
    totals.child_user_s = child_user_us.load() / 1e6;
    totals.child_system_s = child_system_us.load() / 1e6;
    totals.parse_s = parse_us.load() / 1e6;
    totals.gaf_bytes = gaf_bytes_total.load();
    return totals;
}

string GafBaseSiteReadSource::immutable_uri(const string& path) {
    // Only a SQLite database takes locks, and gbz-base passes a database path straight to SQLite,
    // which reads a `file:` name as a URI. Anything else, such as a plain GBZ given as the graph,
    // keeps its path: gbz-base recognises a GBZ by opening the file under the name it is given.
    {
        static const char SQLITE_HEADER[16] = {'S', 'Q', 'L', 'i', 't', 'e', ' ', 'f',
                                               'o', 'r', 'm', 'a', 't', ' ', '3', '\0'};
        char header[16] = {0};
        ifstream in(path, ios::binary);
        if (!in.read(header, sizeof(header)) || memcmp(header, SQLITE_HEADER, sizeof(header)) != 0) {
            return path;
        }
    }
    char* real = realpath(path.c_str(), nullptr);
    if (real == nullptr) {
        throw runtime_error("cannot resolve database path " + path + ": " + strerror(errno));
    }
    string absolute = real;
    free(real);
    if (absolute.find_first_of("?#%") != string::npos) {
        // These characters would have to be escaped in a URI. Rather than get that subtly wrong,
        // refuse, and say how to go back to plain paths.
        throw runtime_error("database path " + absolute + " contains ?, # or %; set "
                            "VG_GAFBASE_LOCKING=1 to give gbz-base the plain path");
    }
    return "file:" + absolute + "?immutable=1";
}

GafBaseSiteReadSource::~GafBaseSiteReadSource() {
    if (query_log != nullptr) {
        fclose(query_log);
    }
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

    // Build argv. A context of 0 keeps the subgraph to the nodes we asked for, since the default
    // context would pull in reads no site here wants. Asking for `overlapping` alignments returns
    // each read whole: the default, `clipped`, can cut one read into several pieces, which would
    // put one read in several rows of the likelihood matrix.
    vector<string> args{binary, "query", gbz_argument};
    args.reserve(args.size() + 2 * nodes.size() + 8);
    for (nid_t node : nodes) {
        args.push_back("-n");
        args.push_back(to_string(node));
    }
    args.push_back("--context");
    args.push_back("0");
    args.push_back("--gaf-base");
    args.push_back(gaf_base_argument);
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
    // The subgraph GFA goes to stdout and we do not want it; only the separate GAF output file
    // interests us. stderr is captured so a failure can say why.
    posix_spawn_file_actions_addopen(&actions, STDOUT_FILENO, "/dev/null", O_WRONLY, 0);
    posix_spawn_file_actions_addopen(&actions, STDERR_FILENO, err_path.c_str(),
                                     O_WRONLY | O_CREAT | O_TRUNC, 0600);

    PendingQuery pending;
    pending.started = std::chrono::steady_clock::now();
    {
        struct timeval now;
        gettimeofday(&now, nullptr);
        pending.started_epoch = now.tv_sec + now.tv_usec / 1e6;
    }

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
    pending.pid = pid;
    pending.n_nodes = nodes.size();
    pending.first_node = nodes.front();
    pending.last_node = nodes.back();
    uint64_t hash = 1469598103934665603ULL;
    for (nid_t node : nodes) {
        hash = (hash ^ (uint64_t)node) * 1099511628211ULL;
    }
    pending.nodes_hash = hash;
    pending.gaf_path = gaf_path(state, slot);
    pending.err_path = std::move(err_path);
    return pending;
}

size_t GafBaseSiteReadSource::reap_query(PendingQuery& pending,
                                         const function<void(Alignment&)>& iteratee) const {
    int child_stat = 0;
    // wait4 rather than waitpid: the same wait, and it also reports the child's CPU time.
    struct rusage usage;
    memset(&usage, 0, sizeof(usage));
    while (wait4(pending.pid, &child_stat, 0, &usage) == -1) {
        if (errno != EINTR) {
            throw runtime_error("wait4() failed for " + binary + ": " + strerror(errno));
        }
    }
    auto reaped = std::chrono::steady_clock::now();
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

    struct stat gaf_stat;
    size_t gaf_bytes = stat(pending.gaf_path.c_str(), &gaf_stat) == 0 ? (size_t)gaf_stat.st_size : 0;

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
    auto parsed_at = std::chrono::steady_clock::now();

    auto micros = [](auto d) {
        return (uint64_t)std::chrono::duration_cast<std::chrono::microseconds>(d).count();
    };
    auto tv_micros = [](const struct timeval& tv) {
        return (uint64_t)tv.tv_sec * 1000000 + (uint64_t)tv.tv_usec;
    };
    uint64_t waited = micros(reaped - pending.started);
    // The iteratee runs inside the parse, so this is parsing plus whatever the caller does with
    // each read.
    uint64_t parsing = micros(parsed_at - reaped);
    uint64_t user = tv_micros(usage.ru_utime);
    uint64_t system = tv_micros(usage.ru_stime);
    wait_us += waited;
    parse_us += parsing;
    child_user_us += user;
    child_system_us += system;
    gaf_bytes_total += gaf_bytes;
    if (query_log != nullptr) {
        std::lock_guard<std::mutex> lock(query_log_mutex);
        fprintf(query_log, "%.6f\t%.3f\t%.3f\t%.3f\t%ld\t%d\t%zu\t%lld\t%lld\t%zu\t%016llx\t%s\t%.3f\t%zu\t%d\n",
                pending.started_epoch, waited / 1e6, user / 1e6, system / 1e6, (long)usage.ru_maxrss, ret,
                pending.n_nodes, (long long)pending.first_node, (long long)pending.last_node, gaf_bytes,
                (unsigned long long)pending.nodes_hash, immutable ? "immutable" : "locking",
                parsing / 1e6, parsed, omp_get_thread_num());
    }
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

////////////////////////////////////////////////////////////////////////////////
// TabixGafSiteReadSource
////////////////////////////////////////////////////////////////////////////////

/// Node IDs closer than this are looked up in the index together, as gbz-base groups a query's
/// nodes. A site whose nodes are scattered over a wide span of IDs then costs a few lookups rather
/// than one per range, without reading every record in the whole span.
static constexpr nid_t TABIX_RUN_GAP = 1000;

/// What the index reader notes about each record it reads, for fetch_span.
struct GafRecordScan {
    /// The fetch's ranges, from merge_ranges.
    const vector<pair<nid_t, nid_t>>* merged = nullptr;
    /// Whether the last record read steps onto a node in the ranges. Also set for a path whose
    /// steps are not all node IDs, which is left to be parsed and tested properly.
    bool may_touch = true;
};

/// Read one record for the tabix iterator, in place of tabix's own reader (tbx_readrec), and give
/// its interval the same way: the smallest and largest number on its path, read as tbx_parse1 reads
/// them. While walking the path it also notes whether a step lands on a node of the fetch, so that
/// a path is read once rather than twice, and does without strtoll for the usual path of node IDs.
static int gaf_readrec(BGZF* fp, void* data, void* r, int* tid, hts_pos_t* beg, hts_pos_t* end) {
    kstring_t* line = (kstring_t*)r;
    GafRecordScan* scan = (GafRecordScan*)data;
    int ret = bgzf_getline(fp, '\n', line);
    if (ret < 0) {
        return ret;
    }
    const char* stop = line->s + line->l;
    const char* path = line->s;
    for (int column = 1; column < 6; ++column) {
        path = (const char*)memchr(path, '\t', stop - path);
        if (path == nullptr) {
            // tbx_readrec cannot place a record without a sixth column either.
            return -2;
        }
        ++path;
    }
    const char* path_end = (const char*)memchr(path, '\t', stop - path);
    if (path_end == nullptr) {
        path_end = stop;
    }

    int64_t lowest = -1, highest = -1;
    bool may_touch = false;
    bool plain = true;
    for (const char* p = path; p < path_end && plain; ) {
        // A step is an orientation followed by a node ID, with no leading zero, which strtoll's
        // base 0 would read as octal.
        if ((*p != '>' && *p != '<') || p + 1 >= path_end || !isdigit((unsigned char)p[1]) ||
            (p[1] == '0' && p + 2 < path_end && isdigit((unsigned char)p[2]))) {
            plain = false;
            break;
        }
        ++p;
        int64_t id = 0;
        for (; p < path_end && isdigit((unsigned char)*p); ++p) {
            id = id * 10 + (*p - '0');
        }
        if (lowest == -1) {
            lowest = highest = id;
        } else {
            lowest = min(lowest, id);
            highest = max(highest, id);
        }
        if (!may_touch && WindowedSiteReadSource::in_ranges((nid_t)id, *scan->merged)) {
            may_touch = true;
        }
    }
    if (!plain) {
        // Some step is not a plain node ID: take the interval exactly as tbx_parse1 does, and let
        // the record be parsed and tested properly.
        lowest = highest = -1;
        for (const char* p = path + 1; p < path_end; ) {
            char* after = nullptr;
            int64_t id = strtoll(p, &after, 0);
            if (lowest == -1) {
                lowest = highest = id;
            } else {
                lowest = min(lowest, id);
                highest = max(highest, id);
            }
            p = after + 1;
        }
        may_touch = true;
    }
    if (lowest < 0 || highest < 0) {
        return -2;
    }
    *tid = 0;
    *beg = lowest;
    *end = highest;
    scan->may_touch = may_touch;
    return ret;
}

TabixGafSiteReadSource::TabixGafSiteReadSource(const HandleGraph& graph,
                                               const string& gaf_filename,
                                               const string& index_filename,
                                               const SiteReadFilter& filter,
                                               size_t window_size,
                                               size_t cache_entries)
    : WindowedSiteReadSource(filter, window_size, cache_entries),
      graph(graph),
      gaf_filename(gaf_filename) {

    index = tbx_index_load2(gaf_filename.c_str(), index_filename.c_str());
    if (index == nullptr) {
        throw runtime_error("could not load the tabix index " + index_filename + " of " + gaf_filename);
    }
    if ((index->conf.preset & 0xffff) != TBX_GAF) {
        tbx_destroy(index);
        index = nullptr;
        throw runtime_error("the tabix index " + index_filename + " was not built for GAF: "
                            "index the GAF with 'tabix -p gaf'");
    }
    // An index addresses a bgzipped file; a gzipped or plain one cannot be read through it. Checked
    // here so that a wrong file fails now, rather than in a worker thread on the first site.
    htsFile* probe = hts_open(gaf_filename.c_str(), "r");
    bool bgzipped = probe != nullptr && hts_get_format(probe)->compression == bgzf;
    if (probe != nullptr) {
        hts_close(probe);
    }
    if (!bgzipped) {
        tbx_destroy(index);
        index = nullptr;
        throw runtime_error(gaf_filename + " cannot be read through an index: it must be the "
                            "sorted GAF compressed with bgzip");
    }

    // One slot per thread, opened lazily, as for an indexed GAM.
    threads.resize(max(1, get_thread_count()));
}

TabixGafSiteReadSource::~TabixGafSiteReadSource() {
    for (ThreadState& state : threads) {
        if (state.file != nullptr) {
            hts_close(state.file);
        }
    }
    if (index != nullptr) {
        tbx_destroy(index);
    }
}

TabixGafSiteReadSource::ThreadState& TabixGafSiteReadSource::thread_state() const {
    int tid = omp_get_thread_num();
    if ((size_t)tid >= threads.size()) {
#pragma omp critical (tabix_gaf_threads)
        if ((size_t)tid >= threads.size()) {
            threads.resize(tid + 1);
        }
    }
    ThreadState& state = threads[tid];
    if (state.file == nullptr) {
        state.file = hts_open(gaf_filename.c_str(), "r");
        if (state.file == nullptr) {
            throw runtime_error("could not open " + gaf_filename + " for reading");
        }
    }
    return state;
}

void TabixGafSiteReadSource::fetch_span(const vector<pair<nid_t, nid_t>>& ranges,
                                        const function<void(Alignment&)>& iteratee) const {
    fetch_span_keyed(ranges, [&](Alignment& aln, uint64_t) { iteratee(aln); });
}

size_t TabixGafSiteReadSource::fetch_span_keyed(
    const vector<pair<nid_t, nid_t>>& ranges,
    const function<void(Alignment&, uint64_t)>& iteratee) const {
    // Reads dropped for having the name and start of a read before them.
    size_t fetch_dropped = 0;
    // Fetches happen inside the OpenMP parallel region that visits the sites, which an exception
    // must not leave. A failed fetch also means a site cannot be scored, so report it and exit.
    try {
        ThreadState& state = thread_state();

        // The ranges merged, to test records against by binary search: a site can name tens of
        // thousands of ranges. Then grouped into runs of nearby node IDs (see TABIX_RUN_GAP), one
        // index lookup per run.
        vector<pair<nid_t, nid_t>> merged = merge_ranges(ranges);
        vector<pair<nid_t, nid_t>> runs(merged);
        size_t kept = 0;
        for (size_t i = 0; i < runs.size(); ++i) {
            if (kept > 0 && runs[i].first <= runs[kept - 1].second + TABIX_RUN_GAP) {
                runs[kept - 1].second = max(runs[kept - 1].second, runs[i].second);
            } else {
                runs[kept++] = runs[i];
            }
        }
        runs.resize(kept);

        // Reads are keyed on name and start rather than name alone, since paired mates share a
        // name, and are passed on once each, as GafBaseSiteReadSource does.
        unordered_set<string> seen;
        gafkluge::GafRecord record;
        Alignment aln;
        // Timings and counts for this fetch, added to the totals once at the end.
        uint64_t fetch_read_us = 0, fetch_parse_us = 0;
        size_t fetch_bytes = 0, fetch_skipped = 0, fetch_unparsed = 0;
        auto micros_since = [](std::chrono::steady_clock::time_point start) {
            return (uint64_t)std::chrono::duration_cast<std::chrono::microseconds>(
                std::chrono::steady_clock::now() - start).count();
        };
        // Takes the record as the reader's buffer holds it, and copies it into a string only to
        // parse it: most records the index returns are dropped unread.
        auto handle_record = [&](const char* text, size_t length, bool may_touch, uint64_t record_key) {
            if (gafkluge::is_gaf_header_line(text)) {
                return;
            }
            // The index says only that the record's node interval overlaps the query. Its path can
            // still step over every node the query names, and gbz-base does not return such a read.
            // Most such records are dropped from what the reader saw on their path, without being
            // parsed.
            if (!may_touch) {
                ++fetch_skipped;
                ++fetch_unparsed;
                return;
            }
            auto parse_start = std::chrono::steady_clock::now();
            gafkluge::parse_gaf_record(string(text, length), record);
            vg::io::gaf_to_alignment(graph, record, aln);
            fetch_parse_us += micros_since(parse_start);
            if (!touches(aln, merged)) {
                ++fetch_skipped;
                return;
            }
            if (!passes_filter(aln)) {
                return;
            }
            count_fetched();
            if (!seen.insert(name_and_start(aln)).second) {
                ++duplicates_dropped;
                ++fetch_dropped;
                return;
            }
            iteratee(aln, record_key);
        };

        // Records the index holds as [smallest node, largest node), and returns when that interval
        // overlaps the query's half-open interval. Widening the query by one ID at each end returns
        // every record whose smallest node is at most the run's last and whose largest node is at
        // least the run's first, including a record on a single node.
        //
        // With one run the records arrive in file order and are handled as they come. With several,
        // a record can come back for more than one run, and a later run can return a record that
        // lies earlier in the file, so they are gathered with the file offset the iterator reached
        // after each, which orders and identifies them, and handled in file order once all are in.
        bool gather = runs.size() > 1;
        struct Gathered {
            uint64_t offset;
            string text;
            bool may_touch;
        };
        vector<Gathered> gathered;
        kstring_t line = KS_INITIALIZE;
        GafRecordScan scan;
        scan.merged = &merged;
        for (const auto& run : runs) {
            hts_itr_t* itr = hts_itr_query(index->idx, 0, max<nid_t>(0, run.first - 1), run.second + 1,
                                           gaf_readrec);
            ++queries;
            if (itr == nullptr) {
                ks_free(&line);
                throw runtime_error("the tabix index of " + gaf_filename + " cannot be queried for "
                                    "node IDs " + to_string(run.first) + "-" + to_string(run.second));
            }
            int ret;
            while (true) {
                auto read_start = std::chrono::steady_clock::now();
                ret = hts_itr_next(hts_get_bgzfp(state.file), itr, &line, &scan);
                fetch_read_us += micros_since(read_start);
                if (ret < 0) {
                    break;
                }
                fetch_bytes += line.l;
                if (gather) {
                    // Only a record that may touch the site is kept whole; the rest are kept as their
                    // offset, to be counted once.
                    gathered.push_back(Gathered{itr->curr_off,
                                                scan.may_touch ? string(line.s, line.l) : string(),
                                                scan.may_touch});
                } else {
                    handle_record(line.s, line.l, scan.may_touch, (uint64_t)itr->curr_off);
                }
            }
            tbx_itr_destroy(itr);
            if (ret < -1) {
                ks_free(&line);
                throw runtime_error("error reading " + gaf_filename + " through its tabix index");
            }
        }
        ks_free(&line);

        if (gather) {
            std::sort(gathered.begin(), gathered.end(), [](const Gathered& a, const Gathered& b) {
                return a.offset < b.offset;
            });
            for (size_t i = 0; i < gathered.size(); ++i) {
                if (i > 0 && gathered[i].offset == gathered[i - 1].offset) {
                    continue;
                }
                handle_record(gathered[i].text.c_str(), gathered[i].text.size(), gathered[i].may_touch,
                              gathered[i].offset);
            }
        }
        read_us += fetch_read_us;
        parse_us += fetch_parse_us;
        gaf_bytes += fetch_bytes;
        skipped += fetch_skipped;
        unparsed += fetch_unparsed;
    } catch (const std::exception& e) {
        cerr << "error[vg::TabixGafSiteReadSource] " << e.what() << endl;
        exit(EXIT_FAILURE);
    }
    return fetch_dropped;
}

}
