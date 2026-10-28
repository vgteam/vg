#include "site_read_source.hpp"

#include <algorithm>
#include <array>
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
#include <sys/socket.h>
#include <sys/resource.h>
#include <sys/stat.h>
#include <sys/time.h>
#include <sys/wait.h>
#include <unistd.h>

#include <vg/io/alignment_io.hpp>
#include <vg/io/stream.hpp>

#include <omp.h>

#include "../src/path.hpp"
#include "../src/utility.hpp"

/// The environment to hand a spawned gaf-base.
///
/// It is declared rather than included because the header that provides it differs between
/// platforms.
extern char** environ;

namespace vg {

using namespace std;

namespace {

/// The most windows a wide site's nodes can lie in for the site to be answered from those windows.
///
/// Those windows are fetched whole and cached. Neighbouring sites need the same windows, so their
/// reads are read once rather than once per site. A site over more windows is fetched by its own
/// node ranges, which reads less than its windows would.
constexpr size_t WIDE_QUERY_MAX_WINDOWS = 16;

/**
 * One hold of `SiteReadSource::hold_reads`, which keeps reads on one thread for the sites of a
 * subtree.
 *
 * The first wide fetch on the thread within the hold fills it. A query is answered from the hold
 * when the hold kept every read on the queried nodes.
 */
struct ReadHold {
    /// The reads beginning on one node with one mapping quality.
    struct Start {
        nid_t node = 0;
        int32_t mapq = 0;
        /// How many reads.
        size_t reads = 0;
        /// Their total sequence length.
        size_t bases = 0;
    };

    /// One mapping of one held read, given by its node, the index of the read in `alignments` and
    /// the index of the mapping in the read's path.
    struct Visit {
        nid_t node = 0;
        uint32_t read = 0;
        uint32_t mapping = 0;
        /// Order by node, then read, then mapping.
        auto operator<=>(const Visit& other) const = default;
    };

    /// The source whose reads are held.
    const void* source = nullptr;
    /// The bar that decides which off-panel candidates become edit sites.
    EditCallBar bar;
    /// Whether a fetch has filled the hold.
    bool filled = false;
    /// Whether the hold kept every read its fetch delivered.
    ///
    /// A fetch from the backend by the site's ranges keeps every read, since the sites below would
    /// fetch them again. A fetch from cached windows keeps only the reads on `nodes`, since the
    /// sites below find the windows cached.
    bool all_reads = false;
    /// Whether the fetch kept every read on the nodes it covered, as a fetch of some of those nodes
    /// would deliver them.
    ///
    /// A fetch that drops a read for having the name and start of an earlier read can lack a read
    /// that a fetch of fewer nodes, which the earlier read misses, delivers.
    bool exact = false;
    /// The nodes on which the hold kept every read, as merged ranges.
    ///
    /// With off-panel detection, these are the site's nodes on candidates that `bar` admits.
    /// Without it, this is empty.
    vector<pair<nid_t, nid_t>> nodes;
    /// The nodes the fetch covered, as merged ranges, kept with off-panel detection or when the
    /// hold is exact.
    vector<pair<nid_t, nid_t>> fetched;
    /// The read starts on the fetched nodes, sorted by node, then mapping quality.
    vector<Start> starts;

    /// The held reads, in fetch order.
    vector<Alignment> alignments;
    /// Every mapping of the held reads, sorted.
    ///
    /// A site below asks for a few of thousands of held reads, and finds them by binary search, so
    /// that its work grows with its own size rather than with the hold's.
    vector<Visit> visits;
    /// The read offset before each mapping, and after the last, for every held read end to end.
    vector<uint32_t> offsets;
    /// Where each read's run of `offsets` starts, and one more entry for the end of the last.
    vector<uint32_t> offset_start;

    /// Add read starts of one read each, summing those that share a node and a mapping quality.
    void add_starts(vector<Start> new_starts) {
        std::ranges::sort(new_starts, [](const Start& a, const Start& b) {
            return a.node != b.node ? a.node < b.node : a.mapq < b.mapq;
        });
        for (const Start& start : new_starts) {
            if (!starts.empty() && starts.back().node == start.node
                && starts.back().mapq == start.mapq) {
                starts.back().reads += start.reads;
                starts.back().bases += start.bases;
            } else {
                starts.push_back(start);
            }
        }
    }

    /// Take and index the reads, which come in fetch order, each with its place in the fetch.
    void hold(vector<pair<size_t, Alignment>>& reads) {
        alignments.reserve(reads.size());
        for (auto& [place, aln] : reads) {
            const auto r = (uint32_t)alignments.size();
            offset_start.push_back((uint32_t)offsets.size());
            uint32_t offset = 0;
            const Path& path = aln.path();
            for (int64_t i = 0; i < path.mapping_size(); ++i) {
                visits.push_back(Visit{path.mapping(i).position().node_id(), r, (uint32_t)i});
                offsets.push_back(offset);
                offset += (uint32_t)mapping_to_length(path.mapping(i));
            }
            offsets.push_back(offset);
            alignments.push_back(std::move(aln));
        }
        offset_start.push_back((uint32_t)offsets.size());
        std::ranges::sort(visits);
    }
};

/// This thread's holds, innermost last.
thread_local vector<ReadHold> read_holds;

/**
 * A sieve for the reads of a wide site's fetch, which hands each read to the site and keeps back
 * the reads that off-panel counting or a hold needs.
 *
 * Kept reads keep their places in the fetch, so that a hold gets them in fetch order.
 */
class WideFetchSieve {
public:
    /// Make a sieve for a site.
    ///
    /// @param merged The site's nodes, as merged ranges.
    /// @param to_count The site's nodes that no fetch has counted, as merged ranges.
    /// @param known The site's nodes on candidates the hold's bar admits, among those counted so
    /// far, as merged ranges.
    /// @param filling The hold to fill, or null.
    /// @param iteratee The site's callback.
    WideFetchSieve(const vector<pair<nid_t, nid_t>>& merged,
                   const vector<pair<nid_t, nid_t>>& to_count, vector<pair<nid_t, nid_t>> known,
                   ReadHold* filling, const function<void(const SiteRead&)>& iteratee)
        : merged(merged), to_count(to_count), known(std::move(known)), filling(filling),
          iteratee(iteratee) {}

    /// Hand a read to the site, keeping it if counting or the hold needs it.
    ///
    /// `owned` is the same read where it may be moved from, or null.
    void take(const Alignment& aln, Alignment* owned) {
        ++n_delivered;
        if (filling != nullptr) {
            // The filter drops reads with no mappings.
            const nid_t first = aln.path().mapping(0).position().node_id();
            if (WindowedSiteReadSource::in_ranges(first, merged)) {
                starts.push_back({first, aln.mapping_quality(), 1, aln.sequence().size()});
            }
        }
        // The read goes without an index. A wide site is rare, and an index for a read the site
        // sees once would cost as much as the walk it saves.
        SiteRead read;
        read.aln = &aln;
        iteratee(read);
        const bool count = !to_count.empty() && WindowedSiteReadSource::touches(aln, to_count);
        const bool keep = filling != nullptr
                          && (keep_all
                              || (!known.empty() && WindowedSiteReadSource::touches(aln, known)));
        if (!count && !keep) {
            return;
        }
        Alignment taken = owned != nullptr ? Alignment(std::move(*owned)) : Alignment(aln);
        if (count) {
            kept_at.push_back(n_delivered);
            kept.push_back(std::move(taken));
        } else {
            held.emplace_back(n_delivered, std::move(taken));
        }
    }

    /// Fill `hold` from now on, with every read.
    void keep_all_for(ReadHold* hold) {
        filling = hold;
        keep_all = true;
    }

    /// Return how many reads were handed to the site.
    size_t delivered() const { return n_delivered; }

    /// Count the off-panel candidates on the uncounted nodes, in the reads kept for counting.
    ///
    /// The hold also gets the kept reads on the nodes of the candidates its bar admits, which
    /// become edit sites.
    vector<EditCandidate> count_edits(const HandleGraph& graph, const EditCountParams& params) {
        vector<EditCandidate> counted = vg::count_edits(
            graph,
            [&](const function<void(const Alignment&)>& visit) {
                for (const Alignment& aln : kept) {
                    visit(aln);
                }
            },
            [&](nid_t node) { return WindowedSiteReadSource::in_ranges(node, to_count); }, params);
        if (filling != nullptr) {
            learn_nodes(counted);
            hold_kept();
        }
        kept.clear();
        return counted;
    }

    /// Fill the hold, if there is one, with the reads kept for it.
    ///
    /// @param exact Whether the fetch kept every read on the nodes it covered.
    /// @param off_panel Whether off-panel detection is on.
    void fill_hold(bool exact, bool off_panel) {
        if (filling == nullptr) {
            return;
        }
        filling->filled = true;
        filling->all_reads = keep_all;
        filling->exact = exact;
        if (off_panel) {
            filling->nodes = std::move(known);
            filling->fetched = merged;
        } else if (exact) {
            filling->fetched = merged;
        }
        filling->add_starts(std::move(starts));
        filling->hold(held);
    }

private:
    /// Add the site's nodes on the newly counted candidates that the hold's bar admits to `known`.
    void learn_nodes(const vector<EditCandidate>& counted) {
        for (const EditCandidate& c : counted) {
            if (!filling->bar.admits(c)) {
                continue;
            }
            known.emplace_back(c.node, c.node);
            for (nid_t node : c.walk) {
                known.emplace_back(node, node);
            }
        }
        std::erase_if(known, [&](const pair<nid_t, nid_t>& r) {
            return r.first == r.second && !WindowedSiteReadSource::in_ranges(r.first, merged);
        });
        known = WindowedSiteReadSource::merge_ranges(std::move(known));
    }

    /// Move the reads kept for counting that the hold wants among those held for it.
    void hold_kept() {
        // Each list is in fetch order, so the two are merged rather than sorted.
        const size_t first_kept = held.size();
        for (size_t i = 0; i < kept.size(); ++i) {
            if (keep_all || WindowedSiteReadSource::touches(kept[i], known)) {
                held.emplace_back(kept_at[i], std::move(kept[i]));
            }
        }
        std::inplace_merge(held.begin(), held.begin() + first_kept, held.end(),
                           [](const pair<size_t, Alignment>& a, const pair<size_t, Alignment>& b) {
                               return a.first < b.first;
                           });
    }

    const vector<pair<nid_t, nid_t>>& merged;
    const vector<pair<nid_t, nid_t>>& to_count;
    /// The site's nodes on candidates the hold's bar admits, as merged ranges.
    vector<pair<nid_t, nid_t>> known;
    /// The hold being filled, or null.
    ReadHold* filling;
    /// Whether the hold gets every read.
    bool keep_all = false;
    const function<void(const SiteRead&)>& iteratee;
    /// How many reads were handed to the site.
    size_t n_delivered = 0;
    /// The reads kept for counting.
    vector<Alignment> kept;
    /// The places of the reads kept for counting, counting from 1, in step with `kept`.
    vector<size_t> kept_at;
    /// The reads kept for the hold alone, each with its place.
    vector<pair<size_t, Alignment>> held;
    /// The read starts on the site's nodes, one per read.
    vector<ReadHold::Start> starts;
};

/// Return whether all the nodes of the first merged ranges are among those of the second.
bool ranges_within(const vector<pair<nid_t, nid_t>>& inner, const vector<pair<nid_t, nid_t>>& outer) {
    return std::ranges::all_of(inner, [&](const pair<nid_t, nid_t>& r) {
        auto it = std::ranges::upper_bound(outer, make_pair(r.first, numeric_limits<nid_t>::max()));
        return it != outer.begin() && std::prev(it)->first <= r.first
               && r.second <= std::prev(it)->second;
    });
}

/// Return what a fetch tells reads apart by: the name, and where the read starts.
///
/// Paired mates share a name. A fetch passes on the first read with each key, in file order.
string name_and_start(const Alignment& aln) {
    string key = aln.name();
    if (aln.path().mapping_size() > 0) {
        const Position& pos = aln.path().mapping(0).position();
        key += "\t" + to_string(pos.node_id()) + "\t" + to_string(pos.offset()) +
               (pos.is_reverse() ? "-" : "+");
    }
    return key;
}

/// Add the nodes of the candidate that lie in the merged ranges to `found`, if `bar` admits it.
void add_admitted_nodes(const EditCandidate& c, const EditCallBar& bar,
                        const vector<pair<nid_t, nid_t>>& merged,
                        vector<pair<nid_t, nid_t>>& found) {
    if (!bar.admits(c)) {
        return;
    }
    if (WindowedSiteReadSource::in_ranges(c.node, merged)) {
        found.emplace_back(c.node, c.node);
    }
    for (nid_t node : c.walk) {
        if (WindowedSiteReadSource::in_ranges(node, merged)) {
            found.emplace_back(node, node);
        }
    }
}

/// Hand the caller each read among the hits, with the index of its hit mappings, and return how
/// many reads were handed over.
///
/// `hits` holds a read's index and a mapping's index for each hit mapping, sorted. A read's
/// alignment is in `alignments`, and its read offsets start at `offsets[offset_start[read]]`.
size_t deliver_hits(const vector<pair<uint32_t, uint32_t>>& hits, const Alignment* alignments,
                    const uint32_t* offsets, const uint32_t* offset_start,
                    const function<void(const SiteRead&)>& iteratee) {
    vector<uint32_t> mappings;
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
        read.aln = &alignments[hits[i].first];
        read.mappings = &mappings[i];
        read.mapping_count = j - i;
        read.read_offsets = &offsets[offset_start[hits[i].first]];
        iteratee(read);
        ++n_delivered;
        i = j;
    }
    return n_delivered;
}

} // namespace


void SiteReadSource::for_each_read_start(
    const vector<nid_t>& nodes,
    const function<void(int32_t mapq, size_t reads, size_t bases)>& iteratee) const {

    // Coalesce the nodes into ID ranges. The source visits each read once however many ranges it
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
        if (nid_t start_node = aln.path().mapping(0).position().node_id();
            !std::ranges::binary_search(nodes, start_node)) {
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

    // Index the read once under every node it touches, though it can visit a node more than once.
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
    // Load serially, since add_read fills a map with no lock.
    get_input_file(filename, [&](istream& in) {
        vg::io::for_each<Alignment>(in, [&](const Alignment& aln) {
            add_read(aln, filter);
        });
    });
}

void InMemorySiteReadSource::load_gaf(const HandleGraph& graph, const string& filename, const Filter& filter) {
    vg::io::gaf_unpaired_for_each(graph, filename, [&](const Alignment& aln) {
        add_read(aln, filter);
    });
}

void InMemorySiteReadSource::for_each_read(const vector<pair<nid_t, nid_t>>& ranges,
                                          const function<void(const SiteRead&)>& iteratee) const {

    // The reads visited so far, so that a read touching several of the ranges is visited once.
    unordered_set<size_t> seen;

    for (const auto& range : ranges) {
        // Look up each ID in the inclusive range. That suits the small ranges a snarl produces.
        // TODO: A range spanning much of the graph would be better served by scanning the map.
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

WindowedSiteReadSource::WindowCache::WindowCache(size_t capacity) : capacity(capacity) {
}

shared_ptr<const WindowedSiteReadSource::CacheEntry>
WindowedSiteReadSource::WindowCache::find_or_claim(size_t window, unique_lock<std::mutex>& lock) {
    auto found = slots.end();
    filled.wait(lock, [&]() {
        found = slots.find(window);
        return found == slots.end() || found->second.entry != nullptr;
    });
    if (found == slots.end()) {
        slots.try_emplace(window);
        return nullptr;
    }
    found->second.last_use = ++clock;
    return found->second.entry;
}

void WindowedSiteReadSource::WindowCache::release_claim(size_t window) {
    slots.erase(window);
    filled.notify_all();
}

void WindowedSiteReadSource::WindowCache::store(size_t window, shared_ptr<const CacheEntry> entry,
                                                vector<shared_ptr<const CacheEntry>>& dropped) {
    Slot& slot = slots[window];
    slot.entry = std::move(entry);
    slot.last_use = ++clock;

    // Drop the least recently used fetched windows other than this one, while there are too many.
    // A claimed window has no entry and does not count.
    while (true) {
        size_t held = 0;
        auto oldest = slots.end();
        for (auto it = slots.begin(); it != slots.end(); ++it) {
            if (!it->second.entry) {
                continue;
            }
            ++held;
            if (it->first != window && (oldest == slots.end()
                                        || it->second.last_use < oldest->second.last_use)) {
                oldest = it;
            }
        }
        if (held <= capacity || oldest == slots.end()) {
            break;
        }
        dropped.push_back(std::move(oldest->second.entry));
        slots.erase(oldest);
    }
    filled.notify_all();
}

vector<shared_ptr<const WindowedSiteReadSource::CacheEntry>>
WindowedSiteReadSource::WindowCache::take_all() {
    unordered_map<size_t, Slot> taken;
    taken.swap(slots);
    vector<shared_ptr<const CacheEntry>> entries;
    entries.reserve(taken.size());
    for (auto& [window, slot] : taken) {
        if (slot.entry) {
            entries.push_back(std::move(slot.entry));
        }
    }
    return entries;
}

WindowedSiteReadSource::WindowedSiteReadSource(const SiteReadFilter& filter,
                                               size_t window_nodes,
                                               size_t max_cached_windows)
    : filter(filter),
      window_size(max<size_t>(1, window_nodes)),
      cache(max<size_t>(1, max_cached_windows)) {
}

size_t WindowedSiteReadSource::window_of(nid_t id) const {
    return (size_t)(id / (nid_t)window_size);
}

bool WindowedSiteReadSource::passes_filter(const Alignment& aln) const {
    if (aln.is_secondary()) {
        ++counters.filtered;
        return false;
    }
    if (aln.mapping_quality() < filter.min_mapq) {
        ++counters.filtered;
        return false;
    }
    if (aln.path().mapping_size() == 0) {
        ++counters.filtered;
        return false;
    }
    return true;
}

void WindowedSiteReadSource::count_fetched() const {
    ++counters.fetched;
}

vector<pair<nid_t, nid_t>> WindowedSiteReadSource::merge_ranges(vector<pair<nid_t, nid_t>> ranges) {
    std::ranges::sort(ranges);
    // Compact in place. The ranges kept so far all lie before the one being read.
    size_t kept = 0;
    for (const auto& range : ranges) {
        if (kept > 0 && range.first <= ranges[kept - 1].second + 1) {
            ranges[kept - 1].second = max(ranges[kept - 1].second, range.second);
        } else {
            ranges[kept++] = range;
        }
    }
    ranges.resize(kept);
    return ranges;
}

bool WindowedSiteReadSource::in_ranges(nid_t node_id, const vector<pair<nid_t, nid_t>>& ranges) {
    // Merged ranges neither overlap nor touch, so only the last range starting at or before the ID
    // can hold it.
    auto after = std::ranges::upper_bound(ranges, node_id, std::less<>(),
                                          [](const pair<nid_t, nid_t>& range) {
                                              return range.first;
                                          });
    return after != ranges.begin() && node_id <= (after - 1)->second;
}

bool WindowedSiteReadSource::touches(const Alignment& aln,
                                    const vector<pair<nid_t, nid_t>>& ranges) {
    return std::ranges::any_of(aln.path().mapping(), [&](const Mapping& mapping) {
        return in_ranges(mapping.position().node_id(), ranges);
    });
}

void WindowedSiteReadSource::for_each_read(
    const vector<pair<nid_t, nid_t>>& ranges,
    const function<void(const SiteRead&)>& iteratee) const {

    if (ranges.empty()) {
        return;
    }

    nid_t min_id = ranges.front().first;
    nid_t max_id = ranges.front().second;
    for (const auto& [first, last] : ranges) {
        min_id = min(min_id, first);
        max_id = max(max_id, last);
    }

    const bool wide = window_of(min_id) != window_of(max_id);
    if (!read_holds.empty() && deliver_held(ranges, wide, iteratee)) {
        return;
    }
    if (wide) {
        deliver_wide(ranges, min_id, max_id, iteratee);
        return;
    }

    // A site inside one window usually finds it cached, when sites are visited in node-ID order
    // (GraphCaller::set_snarl_batching).
    bool was_fetched = false;
    shared_ptr<const CacheEntry> entry = get_window(window_of(min_id), was_fetched);
    if (was_fetched) {
        ++counters.cache_misses;
    } else {
        ++counters.cache_hits;
    }
    deliver(*entry, ranges, iteratee);
}

void WindowedSiteReadSource::deliver_wide(const vector<pair<nid_t, nid_t>>& ranges, nid_t min_id,
                                          nid_t max_id,
                                          const function<void(const SiteRead&)>& iteratee) const {
    // Either way the site gets the same reads in the same order. Reads from the site's windows
    // let neighbouring sites read them once. A site over many windows is fetched by its ranges,
    // since a region of unbounded size is not worth caching. The ranges go to the backend as they
    // are, rather than as one span, since a snarl's nodes can be spread thinly over a wide span of
    // IDs.
    vector<pair<nid_t, nid_t>> merged = merge_ranges(ranges);
    vector<size_t> windows;
    const bool few = few_windows(merged, windows);

    // Fill the innermost hold on this thread, if this is the first wide fetch within it.
    ReadHold* hold = nullptr;
    if (!read_holds.empty() && read_holds.back().source == this && !read_holds.back().filled) {
        hold = &read_holds.back();
    }
    vector<pair<nid_t, nid_t>> to_count;
    vector<pair<nid_t, nid_t>> known;
    if (edit_graph != nullptr) {
        std::scoped_lock guard(cache_mutex);
        to_count = uncounted(merged);
        if (hold != nullptr) {
            known = candidate_nodes(merged, hold->bar);
        }
    }
    // Reads from cached windows fill a hold only with off-panel detection. Without it, the sites
    // below find the windows cached.
    WideFetchSieve sieve(merged, to_count, std::move(known),
                         few && edit_graph == nullptr ? nullptr : hold, iteratee);

    size_t n_scanned = 0;
    size_t dropped = 0;
    bool from_windows = false;
    if (few) {
        // The windows are fetched for this site alone, so their read starts are left untallied and
        // their candidates uncounted until another query uses them.
        from_windows = for_each_window_read(merged, windows, false,
                                            [&](const Alignment& aln) {
                                                ++n_scanned;
                                                sieve.take(aln, nullptr);
                                            },
                                            dropped);
    }
    if (from_windows) {
        ++counters.wide_from_windows;
    } else {
        sieve.keep_all_for(hold);
        ++counters.cache_misses;
        ++counters.straddles;
        counters.straddle_nodes += (size_t)(max_id - min_id + 1);
        size_t wanted = 0;
        for (const auto& [first, last] : ranges) {
            wanted += (size_t)(last - first + 1);
        }
        counters.straddle_wanted += wanted;
        dropped = fetch_span_keyed(ranges, [&](Alignment& aln, uint64_t) {
            ++n_scanned;
            if (touches(aln, merged)) {
                sieve.take(aln, &aln);
            }
        });
    }
    counters.scanned += n_scanned;
    counters.delivered += sieve.delivered();
    if (!to_count.empty()) {
        vector<EditCandidate> counted = sieve.count_edits(*edit_graph, edit_params);
        record_span_edits(std::move(to_count), counted);
    }
    sieve.fill_hold(has_record_keys() && dropped == 0, edit_graph != nullptr);
}

void WindowedSiteReadSource::record_span_edits(vector<pair<nid_t, nid_t>> counted_nodes,
                                               const vector<EditCandidate>& counted) const {
    std::scoped_lock guard(cache_mutex);
    // Another thread may have counted some of the same nodes meanwhile. It found the same
    // candidates there, and the copies are merged away.
    counted_nodes.insert(counted_nodes.end(), span_counted.begin(), span_counted.end());
    span_counted = merge_ranges(std::move(counted_nodes));
    span_edits.insert(span_edits.end(), counted.begin(), counted.end());
    std::ranges::sort(span_edits);
    auto repeats = std::ranges::unique(span_edits);
    span_edits.erase(repeats.begin(), repeats.end());
}

shared_ptr<const WindowedSiteReadSource::CacheEntry>
WindowedSiteReadSource::get_window(size_t window, bool& was_fetched, bool counted) const {
    was_fetched = false;
    // A window cached by wide sites' fetches alone is cached unsummarized. This use summarizes it,
    // as fetching it would have.
    shared_ptr<const CacheEntry> found_entry;
    {
        unique_lock<std::mutex> lock(cache_mutex);
        found_entry = cache.find_or_claim(window, lock);
        if (found_entry && (!counted || summarized(window))) {
            return found_entry;
        }
    }

    shared_ptr<CacheEntry> entry;
    if (!found_entry) {
        entry = load_claimed_window(window);
    }
    bool tallied = false;
    bool edits_counted = true;
    {
        std::scoped_lock guard(cache_mutex);
        tallied = starts.contains(window);
        edits_counted = edit_graph == nullptr || edits.contains(window);
    }
    if (!found_entry && tallied) {
        ++counters.refetches;
    }

    // Only this thread can be fetching the window now. Other threads finding it cached may
    // summarize it at the same time, and find the same tallies and candidates.
    const CacheEntry& reads = found_entry ? *found_entry : *entry;
    vector<StartTally> tallies;
    const bool tally = counted && !tallied;
    if (tally) {
        tallies = tally_starts(reads);
    }
    vector<EditCandidate> window_candidates;
    const bool count = counted && !edits_counted;
    if (count) {
        window_candidates = count_window_edits(reads);
    }

    // Free the windows dropped from the cache after releasing the lock. Freeing a window's reads
    // takes long enough that other threads would queue on the lock.
    vector<shared_ptr<const CacheEntry>> dropped;
    std::scoped_lock guard(cache_mutex);
    if (tally) {
        starts.try_emplace(window, std::move(tallies));
    }
    if (count) {
        edits.try_emplace(window, std::move(window_candidates));
    }
    if (found_entry) {
        return found_entry;
    }
    was_fetched = true;
    cache.store(window, entry, dropped);
    return entry;
}

shared_ptr<WindowedSiteReadSource::CacheEntry>
WindowedSiteReadSource::load_claimed_window(size_t window) const {
    try {
        return make_shared<CacheEntry>(load_window(window));
    } catch (...) {
        std::scoped_lock guard(cache_mutex);
        cache.release_claim(window);
        throw;
    }
}

bool WindowedSiteReadSource::summarized(size_t window) const {
    return starts.contains(window) && (edit_graph == nullptr || edits.contains(window));
}

WindowedSiteReadSource::CacheEntry WindowedSiteReadSource::load_window(size_t window) const {
    CacheEntry entry;
    entry.window = window;
    auto lo = (nid_t)(window * window_size);
    nid_t hi = lo + (nid_t)window_size - 1;
    entry.dropped_repeats = fetch_span_keyed({{lo, hi}}, [&](Alignment& aln, uint64_t key) {
        // Take the alignment rather than copy it. The backends reuse one Alignment per record and
        // clear it before the next.
        index_read(aln, (uint32_t)entry.reads.size(), entry);
        entry.reads.push_back(std::move(aln));
        entry.keys.push_back(key);
    });
    entry.offset_start.push_back((uint32_t)entry.offsets.size());
    // Sort once per fetch rather than per site query.
    std::ranges::sort(entry.node_index);
    ++counters.whole_fetches;
    counters.whole_fetch_reads += entry.reads.size();
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
    std::ranges::sort(tallies, [](const StartTally& a, const StartTally& b) {
        return a.node != b.node ? a.node < b.node : a.mapq < b.mapq;
    });
    // Sum the reads sharing a node and a mapping quality, compacting in place.
    size_t kept = 0;
    for (const StartTally& tally : tallies) {
        if (kept > 0 && tallies[kept - 1].node == tally.node
            && tallies[kept - 1].mapq == tally.mapq) {
            tallies[kept - 1].reads += tally.reads;
            tallies[kept - 1].bases += tally.bases;
        } else {
            tallies[kept++] = tally;
        }
    }
    tallies.resize(kept);
    return tallies;
}

const vector<WindowedSiteReadSource::StartTally>&
WindowedSiteReadSource::window_starts(size_t window) const {
    {
        std::scoped_lock guard(cache_mutex);
        if (auto found = starts.find(window); found != starts.end()) {
            return found->second;
        }
    }
    // The window was never tallied. get_window tallies it before returning. This is not a site
    // query, so it is not counted as a cache hit or miss.
    bool was_fetched = false;
    get_window(window, was_fetched);
    std::scoped_lock guard(cache_mutex);
    return starts.at(window);
}

void WindowedSiteReadSource::count_edits(const HandleGraph* graph,
                                         const EditCountParams& params) {
    edit_graph = graph;
    edit_params = params;
}

vector<EditCandidate> WindowedSiteReadSource::count_window_edits(const CacheEntry& entry) const {
    return vg::count_edits(
        *edit_graph,
        [&](const function<void(const Alignment&)>& visit) {
            for (const Alignment& aln : entry.reads) {
                visit(aln);
            }
        },
        [&](nid_t node) { return window_of(node) == entry.window; }, edit_params);
}

const vector<EditCandidate>& WindowedSiteReadSource::window_edits(size_t window) const {
    {
        std::scoped_lock guard(cache_mutex);
        if (auto found = edits.find(window); found != edits.end()) {
            return found->second;
        }
    }
    // get_window counts the window before returning, unless it was cached before counting was
    // turned on.
    bool was_fetched = false;
    shared_ptr<const CacheEntry> entry = get_window(window, was_fetched);
    if (was_fetched) {
        ++counters.edit_fetches;
    }
    {
        std::scoped_lock guard(cache_mutex);
        if (auto found = edits.find(window); found != edits.end()) {
            return found->second;
        }
    }
    vector<EditCandidate> counted = count_window_edits(*entry);
    std::scoped_lock guard(cache_mutex);
    return edits.try_emplace(window, std::move(counted)).first->second;
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
    if (std::ranges::none_of(read_holds, [&](const ReadHold& hold) {
            return hold.source == this && hold.filled;
        })) {
        return false;
    }
    const vector<pair<nid_t, nid_t>> merged = merge_ranges(ranges);
    for (auto hold = read_holds.rbegin(); hold != read_holds.rend(); ++hold) {
        if (hold->source != this || !hold->filled) {
            continue;
        }
        // The hold has every read on the queried nodes if they all lie on its `nodes`. A wide
        // query would be fetched like the hold's fetch, so it also gets the same reads from a hold
        // that kept every read its fetch delivered, if that fetch covered the query's nodes.
        if (!ranges_within(merged, hold->nodes)
            && !(wide && hold->all_reads && hold->exact && ranges_within(merged, hold->fetched))) {
            continue;
        }
        // Collect the held mappings on the queried nodes, and sort them by read, then by mapping,
        // so that reads are delivered in fetch order and each read's mappings form an ascending
        // run.
        vector<pair<uint32_t, uint32_t>> hits;
        for (const auto& [first, last] : merged) {
            ReadHold::Visit probe;
            probe.node = first;
            auto it = std::ranges::lower_bound(hold->visits, probe);
            for (; it != hold->visits.end() && it->node <= last; ++it) {
                hits.emplace_back(it->read, it->mapping);
            }
        }
        std::ranges::sort(hits);
        // Take the hold's buffers before handing out reads. The caller's code can start a hold of
        // its own, which can move the ReadHold but not its buffers.
        const size_t n_delivered = deliver_hits(hits, hold->alignments.data(), hold->offsets.data(),
                                                hold->offset_start.data(), iteratee);
        ++counters.held_queries;
        counters.scanned += hits.size();
        counters.delivered += n_delivered;
        return true;
    }
    return false;
}

vector<pair<nid_t, nid_t>> WindowedSiteReadSource::candidate_nodes(
    const vector<pair<nid_t, nid_t>>& merged, const EditCallBar& bar) const {
    vector<pair<nid_t, nid_t>> found;
    for (const auto& [first, last] : merged) {
        for (auto it = std::ranges::lower_bound(span_edits, first, std::less<>(),
                                                &EditCandidate::node);
             it != span_edits.end() && it->node <= last; ++it) {
            add_admitted_nodes(*it, bar, merged, found);
        }
        for (size_t w = window_of(first); w <= window_of(last); ++w) {
            auto counted = edits.find(w);
            if (counted == edits.end()) {
                continue;
            }
            for (const EditCandidate& c : counted->second) {
                if (c.node >= first && c.node <= last) {
                    add_admitted_nodes(c, bar, merged, found);
                }
            }
        }
    }
    return merge_ranges(std::move(found));
}

vector<pair<nid_t, nid_t>>
WindowedSiteReadSource::uncounted(const vector<pair<nid_t, nid_t>>& merged) const {
    vector<pair<nid_t, nid_t>> left;
    for (const auto& [first, last] : merged) {
        for (nid_t id = first; id <= last; ++id) {
            if (edits.contains(window_of(id)) || in_ranges(id, span_counted)) {
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
        // Take the candidates from wide sites' fetches where those counted all the nodes and the
        // window was never counted, so that looking them up fetches nothing.
        vector<EditCandidate> spanned;
        bool from_span = false;
        {
            std::scoped_lock guard(cache_mutex);
            from_span = !edits.contains(window) && std::ranges::all_of(in_window, [&](nid_t node) {
                return in_ranges(node, span_counted);
            });
            if (from_span) {
                // Search by node alone. A probe compared field by field would sort after the
                // node's candidates at offset 0 of a lower kind or ALT, and skip them.
                for (auto it = std::ranges::lower_bound(span_edits, in_window.front(),
                                                        std::less<>(), &EditCandidate::node);
                     it != span_edits.end() && it->node <= in_window.back(); ++it) {
                    spanned.push_back(*it);
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
    std::scoped_lock guard(cache_mutex);
    for (const auto& [window, candidates] : edits) {
        all.insert(all.end(), candidates.begin(), candidates.end());
    }
    // A node counted by a wide site's fetch and by its window's has the same candidates in both.
    all.insert(all.end(), span_edits.begin(), span_edits.end());
    std::ranges::sort(all);
    auto repeats = std::ranges::unique(all);
    all.erase(repeats.begin(), repeats.end());
    return all;
}

void WindowedSiteReadSource::for_each_read_start(
    const vector<nid_t>& nodes,
    const function<void(int32_t mapq, size_t reads, size_t bases)>& iteratee) const {

    // Answer from a hold on this thread whose fetch covered every node, since it tallied their
    // starts as their windows would.
    //
    // A backend may keep only the first of the reads sharing a name and a start. They all start
    // on one node, so any fetch covering that node keeps one of them. The hold and the window
    // then count the same number of reads there. They can differ in which read they keep, and so
    // in its mapping quality and length, where the backend returns the reads in a different order
    // for different nodes. Neither choice is more right.
    for (auto hold = read_holds.rbegin(); hold != read_holds.rend(); ++hold) {
        if (hold->source != this || !hold->filled
            || !std::ranges::all_of(nodes, [&](nid_t node) { return in_ranges(node, hold->fetched); })) {
            continue;
        }
        for (nid_t node : nodes) {
            auto it = std::ranges::lower_bound(hold->starts, node, std::less<>(), &ReadHold::Start::node);
            for (; it != hold->starts.end() && it->node == node; ++it) {
                iteratee(it->mapq, it->reads, it->bases);
            }
        }
        ++counters.held_queries;
        return;
    }

    size_t i = 0;
    while (i < nodes.size()) {
        // The nodes are sorted, so each window's nodes are consecutive.
        size_t window = window_of(nodes[i]);
        const vector<StartTally>& tallies = window_starts(window);
        for (; i < nodes.size() && window_of(nodes[i]) == window; ++i) {
            auto it = std::ranges::lower_bound(tallies, nodes[i], std::less<>(), &StartTally::node);
            for (; it != tallies.end() && it->node == nodes[i]; ++it) {
                iteratee(it->mapq, it->reads, it->bases);
            }
        }
    }
}

size_t WindowedSiteReadSource::drop_cached_windows() {
    // Take the windows under the lock and free them after it. Freeing a window's reads takes long
    // enough that other threads would queue on the lock.
    vector<shared_ptr<const CacheEntry>> entries;
    {
        std::scoped_lock guard(cache_mutex);
        entries = cache.take_all();
    }
    size_t reads = 0;
    for (const auto& entry : entries) {
        reads += entry->reads.size();
    }
    // Free the windows on several threads. On a whole genome they hold tens of gigabytes of
    // alignments, which are slow to free on one thread.
#pragma omp parallel for schedule(dynamic, 1)
    for (size_t i = 0; i < entries.size(); ++i) {
        entries[i].reset();
    }
    return reads;
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
    // Add the end of the last mapping too.
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
    // The ranges are sorted and apart, so the windows come out ascending.
    for (const auto& [first, last] : merged) {
        size_t last_window = window_of(last);
        size_t w = window_of(first);
        if (!windows.empty() && w <= windows.back()) {
            w = windows.back() + 1;
        }
        for (; w <= last_window; ++w) {
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

    // Hold the windows until the reads are delivered, even if the cache drops them meanwhile.
    vector<shared_ptr<const CacheEntry>> entries;
    entries.reserve(windows.size());
    for (size_t w : windows) {
        bool was_fetched = false;
        entries.push_back(get_window(w, was_fetched, counted));
        if (entries.back()->dropped_repeats > 0) {
            // The window's fetch dropped a read for having the name and start of an earlier read.
            // If the earlier read touches the window but not the ranges, a fetch of the ranges
            // delivers the dropped read.
            return false;
        }
    }

    // Find every read touching the ranges, from each window's node index. A read in several of
    // the windows is the same record, with the same key, in each.
    struct Hit {
        uint64_t key;
        uint32_t entry;
        uint32_t read;
    };
    vector<Hit> hits;
    vector<uint32_t> found;
    for (uint32_t e = 0; e < (uint32_t)entries.size(); ++e) {
        const CacheEntry& entry = *entries[e];
        const auto lo = (nid_t)(entry.window * window_size);
        const nid_t hi = lo + (nid_t)window_size - 1;
        found.clear();
        for (const auto& [first, last] : merged) {
            if (last < lo || first > hi) {
                continue;
            }
            IndexEntry probe{first, 0, 0};
            auto it = std::ranges::lower_bound(entry.node_index, probe);
            for (; it != entry.node_index.end() && it->node <= last; ++it) {
                found.push_back(it->read);
            }
        }
        std::ranges::sort(found);
        auto repeats = std::ranges::unique(found);
        found.erase(repeats.begin(), repeats.end());
        for (uint32_t r : found) {
            hits.push_back(Hit{entry.keys[r], e, r});
        }
    }
    std::ranges::sort(hits, [](const Hit& a, const Hit& b) {
        return a.key != b.key ? a.key < b.key : a.entry < b.entry;
    });

    // Deliver them as a fetch of the ranges would: in file order, once each, and only the first of
    // reads sharing a name and a start.
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
    // Use local buffers, since the iteratee is caller code that might query the source again.
    vector<pair<uint32_t, uint32_t>> hits;
    // Merge the ranges, so that ranges that overlap find each index entry once.
    for (const auto& [first, last] : merge_ranges(ranges)) {
        IndexEntry probe{first, 0, 0};
        auto it = std::ranges::lower_bound(entry.node_index, probe);
        for (; it != entry.node_index.end() && it->node <= last; ++it) {
            hits.emplace_back(it->read, it->mapping);
        }
    }

    // Sort by read, then by mapping, so reads are delivered in fetch order and each read's
    // mappings form an ascending run. There is one index entry per mapping.
    std::ranges::sort(hits);

    const size_t n_delivered = deliver_hits(hits, entry.reads.data(), entry.offsets.data(),
                                            entry.offset_start.data(), iteratee);

    // Add once per site query, to keep atomic increments off the per-entry path.
    counters.scanned += hits.size();
    counters.delivered += n_delivered;
}

size_t WindowedSiteReadSource::get_read_count() const {
    return counters.fetched.load();
}

WindowedSiteReadSource::FetchStats WindowedSiteReadSource::fetch_stats() const {
    FetchStats stats;
    stats.fetched = counters.fetched.load();
    stats.filtered = counters.filtered.load();
    stats.cache_hits = counters.cache_hits.load();
    stats.cache_misses = counters.cache_misses.load();
    stats.scanned = counters.scanned.load();
    stats.delivered = counters.delivered.load();
    stats.straddles = counters.straddles.load();
    stats.wide_from_windows = counters.wide_from_windows.load();
    stats.straddle_nodes = counters.straddle_nodes.load();
    stats.straddle_wanted = counters.straddle_wanted.load();
    stats.whole_fetches = counters.whole_fetches.load();
    stats.refetches = counters.refetches.load();
    stats.edit_fetches = counters.edit_fetches.load();
    stats.held_queries = counters.held_queries.load();
    stats.whole_fetch_reads = counters.whole_fetch_reads.load();
    return stats;
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
}

IndexedGamSiteReadSource::ThreadState& IndexedGamSiteReadSource::thread_state() const {
    ThreadState& state = threads.mine();
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
    // GAMIndex::find takes a range list and visits each read once across it, as fetch_span must.
    vector<pair<id_t, id_t>> query;
    query.reserve(ranges.size());
    for (const auto& [first, last] : ranges) {
        query.emplace_back((id_t)first, (id_t)last);
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
// TabixGafSiteReadSource
////////////////////////////////////////////////////////////////////////////////

/// What the record reader notes about the last record it read, for the fetch reading it.
struct GafRecordScan {
    /// The fetch's ranges, from merge_ranges.
    const vector<pair<nid_t, nid_t>>* merged = nullptr;
    /// Whether the last record's path may visit a node in the ranges.
    ///
    /// It is false only when every step of the path is a node ID outside the ranges.
    bool may_touch = true;
};

namespace {

/// The most node IDs between two ranges that are looked up together, in the tabix index or in one
/// gaf-base query.
///
/// A site whose nodes are scattered over a wide span of IDs then costs a few lookups rather than
/// one per range, without reading every record in the whole span.
constexpr nid_t TABIX_RUN_GAP = 1000;

/// The smallest and largest node ID on a path, read step by step.
struct NodeInterval {
    /// The smallest ID, or -1 before any step is read.
    int64_t lowest = -1;
    /// The largest ID, or -1 before any step is read.
    int64_t highest = -1;

    /// Widen the interval to hold the ID.
    void add(int64_t id) {
        if (lowest == -1) {
            lowest = id;
            highest = id;
        } else {
            lowest = min(lowest, id);
            highest = max(highest, id);
        }
    }
};

/// Return the start of the sixth tab-separated column of the line, which holds the path, or null
/// if the line has fewer columns.
const char* find_path_column(const char* line, const char* stop) {
    const char* path = line;
    for (int column = 1; column < 6; ++column) {
        path = (const char*)memchr(path, '\t', stop - path);
        if (path == nullptr) {
            return nullptr;
        }
        ++path;
    }
    return path;
}

/// Read a path of plain steps into `interval`, setting `may_touch` if a step visits a node in
/// `merged`, and return false at the first step that is not plain.
///
/// A plain step is an orientation followed by a node ID with no leading zero.
bool scan_plain_path(const char* path, const char* path_end,
                     const vector<pair<nid_t, nid_t>>& merged, NodeInterval& interval,
                     bool& may_touch) {
    const char* p = path;
    while (p < path_end) {
        // A leading zero is not plain, since strtoll's base 0 would read the ID as octal.
        if ((*p != '>' && *p != '<') || p + 1 >= path_end || !isdigit((unsigned char)p[1]) ||
            (p[1] == '0' && p + 2 < path_end && isdigit((unsigned char)p[2]))) {
            return false;
        }
        ++p;
        int64_t id = 0;
        for (; p < path_end && isdigit((unsigned char)*p); ++p) {
            id = id * 10 + (*p - '0');
        }
        interval.add(id);
        if (!may_touch && WindowedSiteReadSource::in_ranges((nid_t)id, merged)) {
            may_touch = true;
        }
    }
    return true;
}

/// Return the interval of a path as tbx_parse1 reads it, from the number strtoll finds after each
/// step's first character.
NodeInterval tabix_path_interval(const char* path, const char* path_end) {
    NodeInterval interval;
    const char* p = path + 1;
    while (p < path_end) {
        char* after = nullptr;
        interval.add(strtoll(p, &after, 0));
        p = after + 1;
    }
    return interval;
}

/// Read one record for the tabix iterator, in place of tabix's own reader (tbx_readrec).
///
/// The record's interval is the smallest and largest number on its path, as tbx_parse1 reads them.
/// While reading the path, this also notes in the GafRecordScan `data` whether a step visits a node
/// of the fetch, so that a path is read once rather than twice.
int gaf_readrec(BGZF* fp, void* data, void* r, int* tid, hts_pos_t* beg, hts_pos_t* end) {
    auto* line = static_cast<kstring_t*>(r);
    auto* scan = static_cast<GafRecordScan*>(data);
    int ret = bgzf_getline(fp, '\n', line);
    if (ret < 0) {
        return ret;
    }
    const char* stop = line->s + line->l;
    const char* path = find_path_column(line->s, stop);
    if (path == nullptr) {
        // tbx_readrec cannot place a record without a sixth column either.
        return -2;
    }
    const auto* path_end = (const char*)memchr(path, '\t', stop - path);
    if (path_end == nullptr) {
        path_end = stop;
    }

    NodeInterval interval;
    bool may_touch = false;
    if (!scan_plain_path(path, path_end, *scan->merged, interval, may_touch)) {
        // Some step is not a plain node ID. Take the interval exactly as tbx_parse1 does, and let
        // the record be parsed and tested properly.
        interval = tabix_path_interval(path, path_end);
        may_touch = true;
    }
    if (interval.lowest < 0 || interval.highest < 0) {
        return -2;
    }
    *tid = 0;
    *beg = interval.lowest;
    *end = interval.highest;
    scan->may_touch = may_touch;
    return ret;
}

/// Return the microseconds since `start`.
uint64_t micros_since(std::chrono::steady_clock::time_point start) {
    return (uint64_t)std::chrono::duration_cast<std::chrono::microseconds>(
        std::chrono::steady_clock::now() - start).count();
}

/// Group merged ranges into runs, joining ranges less than TABIX_RUN_GAP node IDs apart.
vector<pair<nid_t, nid_t>> nearby_runs(const vector<pair<nid_t, nid_t>>& merged) {
    vector<pair<nid_t, nid_t>> runs;
    for (const auto& range : merged) {
        if (!runs.empty() && range.first <= runs.back().second + TABIX_RUN_GAP) {
            runs.back().second = max(runs.back().second, range.second);
        } else {
            runs.push_back(range);
        }
    }
    return runs;
}

/// An htslib line buffer, freed when it goes out of scope.
struct LineBuffer {
    kstring_t line = KS_INITIALIZE;
    LineBuffer() = default;
    LineBuffer(const LineBuffer&) = delete;
    LineBuffer& operator=(const LineBuffer&) = delete;
    ~LineBuffer() { ks_free(&line); }
};

/**
 * A parser that turns the GAF records of one fetch into Alignments, and drops the records whose
 * path visits none of the fetch's nodes.
 *
 * A tabix index says only that a record's node interval overlaps the fetch, and gaf-base is asked
 * for runs of nearby ranges as one range. Either way a record's path can step over every node the
 * fetch names. Most such records are dropped from what the reader saw on their path, without being
 * parsed.
 */
class GafFetchParser {
public:
    /// Make a parser for a fetch of the merged ranges, on the graph the reads were aligned to.
    ///
    /// The optional fields named in `dropped_tags` are left out of the Alignments. They are
    /// about where the backend keeps a record, not about the read.
    GafFetchParser(const HandleGraph& graph, const vector<pair<nid_t, nid_t>>& merged,
                   vector<string> dropped_tags = {})
        : graph(graph), merged(merged), dropped_tags(std::move(dropped_tags)) {}

    /// Parse a record, and return its alignment, or null if the record is a header or visits none
    /// of the fetch's nodes.
    ///
    /// `may_touch` is false where the reader found from the path text alone that the record visits
    /// none of them. The alignment is overwritten by the next call, and the caller may move from
    /// it.
    Alignment* parse(const char* text, size_t length, bool may_touch) {
        if (gafkluge::is_gaf_header_line(text)) {
            return nullptr;
        }
        if (!may_touch) {
            ++skipped;
            ++unparsed;
            return nullptr;
        }
        auto parse_start = std::chrono::steady_clock::now();
        gafkluge::parse_gaf_record(string(text, length), record);
        for (const string& tag : dropped_tags) {
            // An unknown tag would otherwise go into the Alignment's annotations.
            record.opt_fields.erase(tag);
        }
        vg::io::gaf_to_alignment(graph, record, aln);
        parse_us += micros_since(parse_start);
        if (!WindowedSiteReadSource::touches(aln, merged)) {
            ++skipped;
            return nullptr;
        }
        return &aln;
    }

    /// Return how many records were dropped for visiting none of the fetch's nodes.
    size_t skipped_count() const { return skipped; }
    /// Return how many of those were dropped without being parsed.
    size_t unparsed_count() const { return unparsed; }
    /// Return the time spent parsing, in microseconds.
    uint64_t parse_time_us() const { return parse_us; }

private:
    const HandleGraph& graph;
    const vector<pair<nid_t, nid_t>>& merged;
    /// The optional fields left out of the Alignments.
    vector<string> dropped_tags;
    /// The record being parsed.
    gafkluge::GafRecord record;
    /// The alignment of the last record parsed.
    Alignment aln;
    size_t skipped = 0;
    size_t unparsed = 0;
    uint64_t parse_us = 0;
};

} // namespace

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
    // An index addresses a bgzipped file, and a gzipped or plain one cannot be read through it.
    // Check now, so that a wrong file fails here rather than in a worker thread on the first site.
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
}

TabixGafSiteReadSource::~TabixGafSiteReadSource() {
    threads.for_each([](const ThreadState& state) {
        if (state.file != nullptr) {
            hts_close(state.file);
        }
    });
    if (index != nullptr) {
        tbx_destroy(index);
    }
}

TabixGafSiteReadSource::ThreadState& TabixGafSiteReadSource::thread_state() const {
    ThreadState& state = threads.mine();
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

void TabixGafSiteReadSource::for_each_record_in_run(
    htsFile* file, const pair<nid_t, nid_t>& run, GafRecordScan& scan, kstring_t& line,
    uint64_t& read_time_us, const function<void(const kstring_t&, uint64_t)>& visit) const {
    // The index holds each record as [smallest node, largest node), and returns a record when that
    // interval overlaps the query's half-open interval. Widening the run by one ID at each end
    // returns every record whose smallest node is at most the run's last and whose largest node is
    // at least the run's first, including a record on a single node.
    hts_itr_t* itr = hts_itr_query(index->idx, 0, max<nid_t>(0, run.first - 1), run.second + 1,
                                   gaf_readrec);
    ++queries;
    if (itr == nullptr) {
        throw runtime_error("the tabix index of " + gaf_filename + " cannot be queried for "
                            "node IDs " + to_string(run.first) + "-" + to_string(run.second));
    }
    int ret;
    while (true) {
        auto read_start = std::chrono::steady_clock::now();
        ret = hts_itr_next(hts_get_bgzfp(file), itr, &line, &scan);
        read_time_us += micros_since(read_start);
        if (ret < 0) {
            break;
        }
        visit(line, (uint64_t)itr->curr_off);
    }
    tbx_itr_destroy(itr);
    if (ret < -1) {
        throw runtime_error("error reading " + gaf_filename + " through its tabix index");
    }
}

size_t TabixGafSiteReadSource::fetch_span_keyed(
    const vector<pair<nid_t, nid_t>>& ranges,
    const function<void(Alignment&, uint64_t)>& iteratee) const {
    // Reads dropped for having the name and start of an earlier read.
    size_t fetch_dropped = 0;
    // A fetch happens inside the OpenMP parallel region that visits the sites, which an exception
    // must not leave. A failed fetch also means a site cannot be scored, so report it and exit.
    try {
        ThreadState& state = thread_state();

        // Merge the ranges, to test records against by binary search, since a site can name tens
        // of thousands of ranges.
        const vector<pair<nid_t, nid_t>> merged = merge_ranges(ranges);
        GafFetchParser parser(graph, merged);

        // Pass on only the first read with each name and start. GafBaseSiteReadSource does the
        // same.
        unordered_set<string> seen;
        auto handle_record = [&](const char* text, size_t length, bool may_touch, uint64_t record_key) {
            Alignment* aln = parser.parse(text, length, may_touch);
            if (aln == nullptr || !passes_filter(*aln)) {
                return;
            }
            count_fetched();
            if (!seen.insert(name_and_start(*aln)).second) {
                ++duplicates_dropped;
                ++fetch_dropped;
                return;
            }
            iteratee(*aln, record_key);
        };

        // Look up each run of nearby node IDs once. With one run the records arrive in file order
        // and are handled as they come. With several, a record can come back for more than one
        // run, and a later run can return a record that lies earlier in the file. So they are
        // gathered with the file offset the iterator reached after each, which orders and
        // identifies them, and handled in file order once all are in.
        const vector<pair<nid_t, nid_t>> runs = nearby_runs(merged);
        const bool gather = runs.size() > 1;
        struct Gathered {
            uint64_t offset;
            string text;
            bool may_touch;
        };
        vector<Gathered> gathered;
        LineBuffer buffer;
        GafRecordScan scan;
        scan.merged = &merged;
        uint64_t fetch_read_us = 0;
        size_t fetch_bytes = 0;
        for (const auto& run : runs) {
            for_each_record_in_run(state.file, run, scan, buffer.line, fetch_read_us,
                                   [&](const kstring_t& line, uint64_t offset) {
                fetch_bytes += line.l;
                if (gather) {
                    // Keep whole only a record that may touch the site. Keep the rest as their
                    // offset, to be counted once.
                    gathered.push_back(Gathered{offset,
                                                scan.may_touch ? string(line.s, line.l) : string(),
                                                scan.may_touch});
                } else {
                    handle_record(line.s, line.l, scan.may_touch, offset);
                }
            });
        }

        if (gather) {
            std::ranges::sort(gathered, std::less<>(), &Gathered::offset);
            for (size_t i = 0; i < gathered.size(); ++i) {
                if (i > 0 && gathered[i].offset == gathered[i - 1].offset) {
                    continue;
                }
                handle_record(gathered[i].text.c_str(), gathered[i].text.size(), gathered[i].may_touch,
                              gathered[i].offset);
            }
        }
        read_us += fetch_read_us;
        parse_us += parser.parse_time_us();
        gaf_bytes += fetch_bytes;
        skipped += parser.skipped_count();
        unparsed += parser.unparsed_count();
    } catch (const std::exception& e) {
        cerr << "error[vg::TabixGafSiteReadSource] " << e.what() << endl;
        exit(EXIT_FAILURE);
    }
    return fetch_dropped;
}

////////////////////////////////////////////////////////////////////////////////
// GafBaseSiteReadSource
////////////////////////////////////////////////////////////////////////////////

namespace {

/// The tag in which gaf-base batch gives the database row of each alignment's block.
constexpr char GAF_BASE_ROW_TAG[] = "ri";
/// The tag in which gaf-base batch gives each alignment's offset in its block.
constexpr char GAF_BASE_OFFSET_TAG[] = "ro";

/// Return the value of the integer tag `tag` on a GAF line, or -1 if the line has none.
///
/// gaf-base writes its position tags last, so the search runs from the end of the line. A field
/// holds no tab, so each tab starts a field.
int64_t find_int_tag(const char* line, size_t length, const char* tag) {
    for (size_t i = length; i-- > 0;) {
        if (line[i] != '\t' || i + 6 > length) {
            continue;
        }
        if (line[i + 1] == tag[0] && line[i + 2] == tag[1] && line[i + 3] == ':' && line[i + 4] == 'i'
            && line[i + 5] == ':') {
            int64_t value = 0;
            size_t j = i + 6;
            for (; j < length && j < i + 6 + 18 && isdigit((unsigned char)line[j]); ++j) {
                value = value * 10 + (line[j] - '0');
            }
            return j > i + 6 ? value : -1;
        }
    }
    return -1;
}

/// Return the key of a line from gaf-base batch, made from its position tags: the row of the
/// alignment's block in the high 32 bits and its offset in the block in the low 32 bits.
///
/// Keys order lines as the database holds them, which is the order of the GAF it was built from.
/// Throws if the line has no position.
uint64_t gaf_base_record_key(const char* line, size_t length) {
    const int64_t row = find_int_tag(line, length, GAF_BASE_ROW_TAG);
    const int64_t offset = find_int_tag(line, length, GAF_BASE_OFFSET_TAG);
    constexpr int64_t most = numeric_limits<uint32_t>::max();
    if (row < 0 || offset < 0 || row > most || offset > most) {
        throw runtime_error("gaf-base returned a line without its position in the database "
                            "(the ri and ro tags of --ids): " + string(line, min<size_t>(length, 60)));
    }
    return ((uint64_t)row << 32) | (uint64_t)offset;
}

/// Return whether the path of a GAF line may visit a node in the merged ranges, from its text.
///
/// It is false only when every step of the path is a plain node ID outside the ranges.
bool line_may_touch(const char* line, size_t length, const vector<pair<nid_t, nid_t>>& merged) {
    const char* stop = line + length;
    const char* path = find_path_column(line, stop);
    if (path == nullptr) {
        // Let the parser report the line.
        return true;
    }
    const auto* path_end = (const char*)memchr(path, '\t', stop - path);
    if (path_end == nullptr) {
        path_end = stop;
    }
    NodeInterval interval;
    bool may_touch = false;
    if (!scan_plain_path(path, path_end, merged, interval, may_touch)) {
        return true;
    }
    return may_touch;
}

/// Write all of `text` to a socket, and return false, with errno set, if that fails.
///
/// A socket whose other end has closed fails the write rather than raising SIGPIPE.
bool send_all(int fd, const string& text) {
#ifdef MSG_NOSIGNAL
    const int flags = MSG_NOSIGNAL;
#else
    // The socket was set up with SO_NOSIGPIPE instead.
    const int flags = 0;
#endif
    size_t sent = 0;
    while (sent < text.size()) {
        ssize_t n = send(fd, text.data() + sent, text.size() - sent, flags);
        if (n < 0) {
            if (errno == EINTR) {
                continue;
            }
            return false;
        }
        sent += (size_t)n;
    }
    return true;
}

/// Guards the making of the descriptors for a gaf-base process until the process has started.
///
/// They are made close-on-exec in a separate step, which is not atomic everywhere, so no other
/// process may start in between and inherit them. A process that inherited the end of another's
/// query socket would keep that process from ever seeing the end of its input.
std::mutex spawn_mutex;

/// Make a descriptor close in any process started from now on.
void close_on_exec(int fd) {
    fcntl(fd, F_SETFD, fcntl(fd, F_GETFD) | FD_CLOEXEC);
}

/// Count a fetch as in progress on its thread for as long as this exists (see
/// GafBaseSiteReadSource::ThreadState).
struct FetchDepth {
    size_t& depth;
    explicit FetchDepth(size_t& depth) : depth(depth) { ++depth; }
    ~FetchDepth() { --depth; }
};

} // namespace

GafBaseSiteReadSource::GafBaseSiteReadSource(const HandleGraph& graph,
                                             const string& gaf_base_filename,
                                             const SiteReadFilter& filter,
                                             size_t window_size,
                                             size_t cache_entries,
                                             const string& binary,
                                             const vector<string>& batch_arguments)
    : WindowedSiteReadSource(filter, window_size, cache_entries),
      graph(graph),
      gaf_base_filename(gaf_base_filename),
      binary(binary),
      batch_arguments(batch_arguments),
      database_argument(gaf_base_filename) {

    const char* locking = getenv("VG_GAFBASE_LOCKING");
    immutable = !(locking != nullptr && string(locking) == "1");
    if (immutable) {
        database_argument = immutable_uri(gaf_base_filename);
    }

    const char* progress = getenv("VG_GAFBASE_PROGRESS_DIR");
    if (progress != nullptr && *progress != '\0') {
        progress_dir = progress;
    }

    const char* log_path = getenv("VG_GAFBASE_QUERY_LOG");
    if (log_path != nullptr && *log_path != '\0') {
        query_log = fopen(log_path, "w");
        if (query_log == nullptr) {
            throw runtime_error("cannot write VG_GAFBASE_QUERY_LOG " + string(log_path) + ": " +
                                strerror(errno));
        }
        // Buffer by line, so the log is complete up to the last finished query however vg ends.
        setvbuf(query_log, nullptr, _IOLBF, 0);
        fprintf(query_log, "#start_epoch\tprocess\tfirst_node\tlast_node\tlines\tgaf_bytes\twait_s\t"
                           "total_s\n");
    }
}

GafBaseSiteReadSource::~GafBaseSiteReadSource() {
    stop_processes();
    if (query_log != nullptr) {
        fclose(query_log);
    }
}

string GafBaseSiteReadSource::immutable_uri(const string& path) {
    // Only a SQLite database takes locks. gaf-base passes the database argument straight to
    // SQLite, which reads a `file:` name as a URI. Any other file keeps its path, so gaf-base can
    // say what is wrong with it.
    {
        static constexpr array<char, 16> SQLITE_HEADER{'S', 'Q', 'L', 'i', 't', 'e', ' ', 'f',
                                                       'o', 'r', 'm', 'a', 't', ' ', '3', '\0'};
        array<char, 16> header{};
        ifstream in(path, ios::binary);
        if (!in.read(header.data(), header.size()) || header != SQLITE_HEADER) {
            return path;
        }
    }
    // realpath allocates its result with malloc, so it is released with free.
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
                            "VG_GAFBASE_LOCKING=1 to give gaf-base the plain path");
    }
    return "file:" + absolute + "?immutable=1";
}

unique_ptr<GafBaseSiteReadSource::Process> GafBaseSiteReadSource::start_process() const {
    auto process = make_unique<Process>();
    process->number = processes++;
    if (!progress_dir.empty()) {
        process->err_path = progress_dir + "/gaf-base." + to_string(process->number) + ".log";
        process->keep_err = true;
    } else {
        // temp_file::create holds a lock, so threads can create files at the same time.
        process->err_path = temp_file::create("vg-gafbase-err-");
    }

    // Ask for each line's position in the database, which keys the reads.
    vector<string> args{binary, "batch", "--ids"};
    if (!progress_dir.empty()) {
        args.emplace_back("--progress");
    }
    args.insert(args.end(), batch_arguments.begin(), batch_arguments.end());
    args.push_back(database_argument);
    vector<char*> argv;
    argv.reserve(args.size() + 1);
    for (string& arg : args) {
        argv.push_back(arg.data());
    }
    argv.push_back(nullptr);

    int query_sockets[2] = {-1, -1};
    int answer_pipe[2] = {-1, -1};
    pid_t pid = 0;
    int spawn_err = 0;
    {
        std::scoped_lock guard(spawn_mutex);
        if (socketpair(AF_UNIX, SOCK_STREAM, 0, query_sockets) != 0) {
            throw runtime_error("cannot make a socket for gaf-base: " + string(strerror(errno)));
        }
        if (pipe(answer_pipe) != 0) {
            int err = errno;
            close(query_sockets[0]);
            close(query_sockets[1]);
            throw runtime_error("cannot make a pipe for gaf-base: " + string(strerror(err)));
        }
        for (int fd : {query_sockets[0], query_sockets[1], answer_pipe[0], answer_pipe[1]}) {
            close_on_exec(fd);
        }

        // Use posix_spawn rather than fork and exec. fork() in a process whose threads are
        // allocating makes libc hold a lock around malloc, stalling the other threads. posix_spawn
        // does not copy the address space, so there is no such lock, and no stdio buffers need
        // flushing. dup2 gives the process its ends without the close-on-exec flag.
        posix_spawn_file_actions_t actions;
        if (posix_spawn_file_actions_init(&actions) != 0) {
            spawn_err = errno;
        } else {
            posix_spawn_file_actions_adddup2(&actions, query_sockets[1], STDIN_FILENO);
            posix_spawn_file_actions_adddup2(&actions, answer_pipe[1], STDOUT_FILENO);
            // Capture stderr so a failure can say why, rather than just reporting a code.
            posix_spawn_file_actions_addopen(&actions, STDERR_FILENO, process->err_path.c_str(),
                                             O_WRONLY | O_CREAT | O_TRUNC, 0600);
            spawn_err = posix_spawnp(&pid, binary.c_str(), &actions, nullptr, argv.data(), environ);
            posix_spawn_file_actions_destroy(&actions);
        }
    }
    // The process has its own copies of its ends.
    close(query_sockets[1]);
    close(answer_pipe[1]);
    if (spawn_err != 0) {
        close(query_sockets[0]);
        close(answer_pipe[0]);
        if (!process->keep_err) {
            temp_file::remove(process->err_path);
        }
        if (spawn_err == ENOENT) {
            throw runtime_error("could not execute '" + binary + "'. Install gbz-base "
                                "(https://github.com/jltsiren/gbz-base), which provides gaf-base, "
                                "and put gaf-base on your PATH, or pass --gaf-base-binary with its "
                                "location.");
        }
        throw runtime_error("posix_spawnp() failed for " + binary + ": " + strerror(spawn_err));
    }
#ifdef SO_NOSIGPIPE
    // Where send() has no MSG_NOSIGNAL, the socket itself is set not to raise SIGPIPE.
    int on = 1;
    setsockopt(query_sockets[0], SOL_SOCKET, SO_NOSIGPIPE, &on, sizeof(on));
#endif
#ifdef F_SETPIPE_SZ
    // A bigger pipe lets the process run further ahead of the parsing. Failing to grow it is fine.
    fcntl(answer_pipe[0], F_SETPIPE_SZ, 1 << 20);
#endif
    process->pid = pid;
    process->query_fd = query_sockets[0];
    process->answer_fd = answer_pipe[0];
    process->buffer.resize(1 << 20);
    return process;
}

string GafBaseSiteReadSource::process_errors(const Process& process) {
    string message;
    {
        ifstream err_in(process.err_path);
        stringstream buffer;
        buffer << err_in.rdbuf();
        message = buffer.str();
    }
    // Keep the end, which says what failed, and drop the trailing newline.
    constexpr size_t MOST = 2000;
    if (message.size() > MOST) {
        message = "..." + message.substr(message.size() - MOST);
    }
    while (!message.empty() && (message.back() == '\n' || message.back() == '\r')) {
        message.pop_back();
    }
    return message;
}

string GafBaseSiteReadSource::stop_process(Process& process) const {
    // The peak memory of gaf-base, read while the process is still there. wait4's ru_maxrss would
    // count vg's memory too: posix_spawn runs the new process in vg's address space until it
    // starts gaf-base, and the kernel keeps the larger peak. VmHWM belongs to gaf-base's own
    // address space. Where there is no /proc, it stays 0.
    size_t peak = 0;
    {
        ifstream status("/proc/" + to_string(process.pid) + "/status");
        for (string line; getline(status, line);) {
            if (line.rfind("VmHWM:", 0) == 0) {
                istringstream fields(line.substr(6));
                size_t kib = 0;
                fields >> kib;
                peak = kib * 1024;
                break;
            }
        }
    }

    // The end of its input ends gaf-base batch, as an empty query would.
    close(process.query_fd);
    process.query_fd = -1;
    // Read what it still writes, so it is not left blocked on a full pipe.
    while (true) {
        ssize_t got = read(process.answer_fd, process.buffer.data(), process.buffer.size());
        if (got > 0 || (got < 0 && errno == EINTR)) {
            continue;
        }
        break;
    }
    close(process.answer_fd);
    process.answer_fd = -1;

    int child_stat = 0;
    // wait4 rather than waitpid, since it also reports the process's CPU time.
    struct rusage usage;
    memset(&usage, 0, sizeof(usage));
    while (wait4(process.pid, &child_stat, 0, &usage) == -1) {
        if (errno != EINTR) {
            return "wait4() failed for " + binary + ": " + strerror(errno);
        }
    }
    auto tv_micros = [](const struct timeval& tv) {
        return (uint64_t)tv.tv_sec * 1000000 + (uint64_t)tv.tv_usec;
    };
    {
        std::scoped_lock guard(stopped_mutex);
        child_user_us += tv_micros(usage.ru_utime);
        child_system_us += tv_micros(usage.ru_stime);
        child_max_rss = max(child_max_rss, peak);
        child_rss_sum += peak;
    }

    string problem;
    if (!WIFEXITED(child_stat) || WEXITSTATUS(child_stat) != 0) {
        string errors = process_errors(process);
        problem = binary + " batch process " + to_string(process.number) + " ended with " +
                  (WIFEXITED(child_stat) ? "exit code " + to_string(WEXITSTATUS(child_stat))
                                         : "signal " + to_string(WTERMSIG(child_stat))) +
                  (errors.empty() ? "" : ": " + errors);
    }
    if (!process.keep_err) {
        temp_file::remove(process.err_path);
    }
    return problem;
}

void GafBaseSiteReadSource::stop_processes() {
    threads.for_each([&](ThreadState& state) {
        for (unique_ptr<Process>& process : state.processes) {
            string problem = stop_process(*process);
            if (!problem.empty()) {
                // The reads it returned were all used, so this is worth a word but not a failure.
                cerr << "warning[vg::GafBaseSiteReadSource] " << problem << endl;
            }
        }
        state.processes.clear();
    });
}

bool GafBaseSiteReadSource::read_line(Process& process, const char*& line, size_t& length,
                                      uint64_t& wait_time_us) const {
    while (true) {
        char* data = process.buffer.data();
        const auto* newline = (const char*)memchr(data + process.begin + process.scanned, '\n',
                                                  process.end - process.begin - process.scanned);
        if (newline != nullptr) {
            line = data + process.begin;
            length = (size_t)(newline - line);
            process.begin += length + 1;
            process.scanned = 0;
            return true;
        }
        process.scanned = process.end - process.begin;
        // Move the unfinished line to the front, and make room if it fills the buffer.
        if (process.begin > 0) {
            memmove(data, data + process.begin, process.end - process.begin);
            process.end -= process.begin;
            process.begin = 0;
        }
        if (process.end == process.buffer.size()) {
            process.buffer.resize(process.buffer.size() * 2);
            data = process.buffer.data();
        }
        auto wait_start = std::chrono::steady_clock::now();
        ssize_t got = read(process.answer_fd, data + process.end, process.buffer.size() - process.end);
        wait_time_us += micros_since(wait_start);
        if (got < 0) {
            if (errno == EINTR) {
                continue;
            }
            throw runtime_error("cannot read the answer of " + binary + ": " + strerror(errno));
        }
        if (got == 0) {
            return false;
        }
        process.end += (size_t)got;
    }
}

void GafBaseSiteReadSource::query(Process& process, nid_t first, nid_t last, uint64_t& wait_time_us,
                                  const function<void(const char*, size_t)>& visit) const {
    auto started = std::chrono::steady_clock::now();
    // gaf-base reads node IDs as unsigned, and IDs start from 1.
    first = max<nid_t>(0, first);
    last = max<nid_t>(0, last);
    const string request = to_string(first) + "\t" + to_string(last) + "\n";
    const string where = "node IDs " + to_string(first) + "-" + to_string(last);
    if (!send_all(process.query_fd, request)) {
        string errors = process_errors(process);
        throw runtime_error(binary + " batch stopped before a query for " + where + " (" +
                            strerror(errno) + ")" + (errors.empty() ? "" : ": " + errors));
    }
    ++queries;

    size_t lines = 0;
    size_t bytes = 0;
    uint64_t waited = 0;
    const char* line = nullptr;
    size_t length = 0;
    while (true) {
        if (!read_line(process, line, length, waited)) {
            string errors = process_errors(process);
            throw runtime_error(binary + " batch stopped in its answer for " + where +
                                (errors.empty() ? "" : ": " + errors));
        }
        if (length == 0) {
            // An empty line ends the answer.
            break;
        }
        ++lines;
        bytes += length + 1;
        visit(line, length);
    }
    wait_time_us += waited;

    if (query_log != nullptr) {
        double total_s = micros_since(started) / 1e6;
        struct timeval now;
        gettimeofday(&now, nullptr);
        double started_epoch = now.tv_sec + now.tv_usec / 1e6 - total_s;
        std::scoped_lock lock(query_log_mutex);
        fprintf(query_log, "%.6f\t%zu\t%lld\t%lld\t%zu\t%zu\t%.6f\t%.6f\n", started_epoch,
                process.number, (long long)first, (long long)last, lines, bytes, waited / 1e6, total_s);
    }
}

void GafBaseSiteReadSource::fetch_span(const vector<pair<nid_t, nid_t>>& ranges,
                                       const function<void(Alignment&)>& iteratee) const {
    fetch_span_keyed(ranges, [&](Alignment& aln, uint64_t) { iteratee(aln); });
}

size_t GafBaseSiteReadSource::fetch_span_keyed(
    const vector<pair<nid_t, nid_t>>& ranges,
    const function<void(Alignment&, uint64_t)>& iteratee) const {
    // Reads dropped for having the name and start of an earlier read.
    size_t fetch_dropped = 0;
    // A fetch happens inside the OpenMP parallel region that visits the sites, which an exception
    // must not leave. A failed fetch also means a site cannot be scored, so report it and exit.
    try {
        ThreadState& state = threads.mine();
        if (state.depth == state.processes.size()) {
            state.processes.push_back(start_process());
        }
        Process& process = *state.processes[state.depth];
        FetchDepth in_fetch(state.depth);

        // Merge the ranges, to test records against by binary search, since a site can name tens
        // of thousands of ranges.
        const vector<pair<nid_t, nid_t>> merged = merge_ranges(ranges);
        // gaf-base's position tags say where the database keeps a read, not what the read is.
        GafFetchParser parser(graph, merged, {GAF_BASE_ROW_TAG, GAF_BASE_OFFSET_TAG});

        // Pass on only the first read with each name and start. TabixGafSiteReadSource does the
        // same.
        unordered_set<string> seen;
        auto handle_record = [&](const char* text, size_t length, bool may_touch, uint64_t record_key) {
            Alignment* aln = parser.parse(text, length, may_touch);
            if (aln == nullptr || !passes_filter(*aln)) {
                return;
            }
            count_fetched();
            if (!seen.insert(name_and_start(*aln)).second) {
                ++duplicates_dropped;
                ++fetch_dropped;
                return;
            }
            iteratee(*aln, record_key);
        };

        // Ask for each run of nearby node IDs once, as the tabix source looks each up once. With one
        // run the reads arrive in file order and are handled as they come, while gaf-base decodes
        // the next ones. With several, a read can come back for more than one run, and a later run
        // can return a read that lies earlier in the file. So they are gathered with their keys,
        // and handled in file order once all are in.
        const vector<pair<nid_t, nid_t>> runs = nearby_runs(merged);
        const bool gather = runs.size() > 1;
        struct Gathered {
            uint64_t key;
            string text;
            bool may_touch;
        };
        vector<Gathered> gathered;
        uint64_t fetch_wait_us = 0;
        size_t fetch_bytes = 0;
        for (const auto& run : runs) {
            // gaf-base returns the reads with a node anywhere in the run, which takes in the IDs
            // between its ranges. A run of a single range needs no look at the paths.
            auto range = std::ranges::lower_bound(merged, run.first, std::less<>(),
                                                  &pair<nid_t, nid_t>::first);
            const bool one_range = range != merged.end() && *range == run;
            // The key of the last line of the answer, to check that the answer is in file order.
            uint64_t last_key = 0;
            bool any = false;
            query(process, run.first, run.second, fetch_wait_us, [&](const char* line, size_t length) {
                fetch_bytes += length + 1;
                const uint64_t key = gaf_base_record_key(line, length);
                const bool may_touch = one_range || line_may_touch(line, length, merged);
                if (gather) {
                    // Keep whole only a record that may touch the site. Keep the rest as their
                    // key, to be counted once.
                    gathered.push_back(Gathered{key, may_touch ? string(line, length) : string(),
                                                may_touch});
                } else {
                    // gaf-base answers in the order of its index of node intervals, which is the
                    // database's order only when it was built from a GAF sorted by that interval.
                    // Reads handled as they come must be in file order, or the calls would depend
                    // on how the windows fell.
                    if (any && key <= last_key) {
                        throw runtime_error("the GAF-Base " + gaf_base_filename + " returned reads "
                                            "out of order: build it from a GAF sorted by "
                                            "'vg gamsort -G' or 'gaf-base sort'");
                    }
                    last_key = key;
                    any = true;
                    handle_record(line, length, may_touch, key);
                }
            });
        }

        if (gather) {
            std::ranges::sort(gathered, std::less<>(), &Gathered::key);
            for (size_t i = 0; i < gathered.size(); ++i) {
                if (i > 0 && gathered[i].key == gathered[i - 1].key) {
                    continue;
                }
                handle_record(gathered[i].text.c_str(), gathered[i].text.size(), gathered[i].may_touch,
                              gathered[i].key);
            }
        }
        wait_us += fetch_wait_us;
        parse_us += parser.parse_time_us();
        gaf_bytes += fetch_bytes;
        skipped += parser.skipped_count();
    } catch (const std::exception& e) {
        cerr << "error[vg::GafBaseSiteReadSource] " << e.what() << endl;
        exit(EXIT_FAILURE);
    }
    return fetch_dropped;
}

void GafBaseSiteReadSource::check_setup() const {
    // Query the first node in the graph. A missing binary, a gaf-base whose batch cannot tag lines
    // with their positions, or an unreadable database fails here.
    nid_t probe = 0;
    graph.for_each_handle((function<bool(const handle_t&)>)[&](const handle_t& handle) {
        probe = graph.get_id(handle);
        return false;
    });
    if (probe == 0) {
        return;
    }

    ThreadState& state = threads.mine();
    if (state.processes.empty()) {
        state.processes.push_back(start_process());
    }
    Process& process = *state.processes.front();
    uint64_t waited = 0;
    try {
        query(process, probe, probe, waited, [](const char*, size_t) {
            // Only whether the query runs matters.
        });
    } catch (const std::exception& e) {
        string errors = process_errors(process);
        if (errors.find("unrecognized subcommand") != string::npos
            || errors.find("--ids") != string::npos) {
            throw runtime_error(binary + " is too old: vg needs `gaf-base batch --ids`, in gbz-base "
                                "from October 2026 (https://github.com/jltsiren/gbz-base). It said: " +
                                errors);
        }
        throw;
    }
}

GafBaseSiteReadSource::QueryTotals GafBaseSiteReadSource::get_query_totals() const {
    QueryTotals totals;
    totals.processes = processes.load();
    totals.queries = queries.load();
    totals.wait_s = wait_us.load() / 1e6;
    totals.parse_s = parse_us.load() / 1e6;
    totals.gaf_bytes = gaf_bytes.load();
    std::scoped_lock guard(stopped_mutex);
    totals.child_user_s = child_user_us / 1e6;
    totals.child_system_s = child_system_us / 1e6;
    totals.child_max_rss = child_max_rss;
    totals.child_rss_sum = child_rss_sum;
    return totals;
}

}
