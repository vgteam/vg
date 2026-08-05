#include "site_read_source.hpp"

#include <algorithm>
#include <cerrno>
#include <cstring>
#include <sstream>
#include <unordered_set>

#include <fcntl.h>
#include <sys/wait.h>
#include <unistd.h>

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
// WindowedSiteReadSource
////////////////////////////////////////////////////////////////////////////////

WindowedSiteReadSource::WindowedSiteReadSource(const SiteReadFilter& filter,
                                               size_t window_size,
                                               size_t cache_entries)
    : filter(filter),
      window_size(max<size_t>(1, window_size)),
      cache_entries(max<size_t>(1, cache_entries)) {

    caches.resize(max(1, get_thread_count()));
}

WindowedSiteReadSource::CacheState& WindowedSiteReadSource::cache_state() const {
    int tid = omp_get_thread_num();
    if ((size_t)tid >= caches.size()) {
        // More threads than we sized for; grow rather than misbehave.
#pragma omp critical (windowed_site_read_caches)
        if ((size_t)tid >= caches.size()) {
            caches.resize(tid + 1);
        }
    }
    CacheState& state = caches[tid];
    if (!state.initialized) {
        state.cache.resize(cache_entries);
        state.initialized = true;
    }
    return state;
}

size_t WindowedSiteReadSource::window_of(nid_t id) const {
    return (size_t)(id / (nid_t)window_size);
}

size_t WindowedSiteReadSource::get_cache_entries() const {
    return cache_entries;
}

bool WindowedSiteReadSource::passes_filter(const Alignment& aln) const {
    if (filter.skip_secondary && aln.is_secondary()) {
        ++filtered;
        return false;
    }
    if (aln.mapping_quality() < filter.min_mapq) {
        ++filtered;
        return false;
    }
    if (filter.skip_unmapped && aln.path().mapping_size() == 0) {
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
    const function<void(const Alignment&)>& iteratee) const {

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
        // The site straddles a window boundary. Fetch the span directly rather than
        // stitching windows together: fetch_span already emits each read at most
        // once, so this sidesteps de-duplicating reads that span the boundary. A
        // site large enough to do this is rare, and caching a region of unbounded
        // size would defeat the point of windowing.
        ++cache_misses;
        fetch_span(min_id, max_id, [&](const Alignment& aln) {
            if (touches(aln, ranges)) {
                iteratee(aln);
            }
        });
        return;
    }

    CacheState& state = cache_state();

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
    nid_t lo = (nid_t)(first_window * window_size);
    nid_t hi = lo + (nid_t)window_size - 1;
    fetch_span(lo, hi, [&](const Alignment& aln) {
        entry.reads.push_back(aln);
    });
    entry.valid = true;

    for (const Alignment& aln : entry.reads) {
        if (touches(aln, ranges)) {
            iteratee(aln);
        }
    }

    state.cache[state.next_evict] = std::move(entry);
    state.next_evict = (state.next_evict + 1) % state.cache.size();
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

void IndexedGamSiteReadSource::fetch_span(nid_t min_id, nid_t max_id,
                                          const function<void(const Alignment&)>& iteratee) const {
    ThreadState& state = thread_state();
    vector<pair<id_t, id_t>> query{{(id_t)min_id, (id_t)max_id}};
    index->find(*state.cursor, query, [&](const Alignment& aln) {
        if (!passes_filter(aln)) {
            return;
        }
        count_fetched();
        iteratee(aln);
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
        if (!state.gaf_path.empty()) {
            temp_file::remove(state.gaf_path);
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
    ThreadState& state = threads[tid];
    if (state.gaf_path.empty()) {
        // One output file per thread, reused for every query that thread makes.
        // temp_file::create is mutex-guarded, so this is safe to race into.
        state.gaf_path = temp_file::create("vg-gafbase-reads-");
    }
    return state;
}

size_t GafBaseSiteReadSource::run_query(ThreadState& state, const vector<nid_t>& nodes,
                                        const function<void(const Alignment&)>& iteratee) const {
    if (nodes.empty()) {
        return 0;
    }

    // Build argv. --context 0 matters: the default of 100bp would expand the
    // subgraph past the nodes we asked for and pull in reads no site here wants.
    // --alignments overlapping matters more: the CLI default of `clipped` can cut one
    // read into several fragments, which would put a single read in several rows of
    // the likelihood matrix and break the per-read independence the model assumes.
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
    args.push_back(state.gaf_path);
    args.push_back("--alignments");
    args.push_back("overlapping");

    vector<const char*> argv;
    argv.reserve(args.size() + 1);
    for (const string& arg : args) {
        argv.push_back(arg.c_str());
    }
    argv.push_back(nullptr);

    // Capture stderr so a failure can say why, rather than just reporting a code.
    string err_path = state.gaf_path + ".err";

    cout.flush();
    cerr.flush();
    fflush(nullptr);

    pid_t pid = fork();
    if (pid == -1) {
        throw runtime_error("fork() failed for " + binary + ": " + strerror(errno));
    } else if (pid == 0) {
        // Child. The subgraph GFA goes to stdout and we do not want it; only the
        // separate --gaf-output file interests us.
        int devnull = open("/dev/null", O_WRONLY);
        if (devnull != -1) {
            dup2(devnull, STDOUT_FILENO);
            close(devnull);
        }
        int errfd = open(err_path.c_str(), O_WRONLY | O_CREAT | O_TRUNC, 0600);
        if (errfd != -1) {
            dup2(errfd, STDERR_FILENO);
            close(errfd);
        }
        // execvp promises not to modify its arguments but cannot say so in C's type
        // system; see the same cast in index_registry.cpp's kmc call.
        execvp(binary.c_str(), (char**)&argv[0]);
        _exit(127);
    }

    int child_stat = 0;
    while (waitpid(pid, &child_stat, 0) == -1) {
        if (errno != EINTR) {
            throw runtime_error("waitpid() failed for " + binary + ": " + strerror(errno));
        }
    }
    ++queries;

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
        if (ret == 127) {
            throw runtime_error("could not execute '" + binary + "'. Install gbz-base "
                                "(https://github.com/jltsiren/gbz-base) and put it on your "
                                "PATH, or pass --gaf-base-binary with its location.");
        }
        throw runtime_error(binary + " query failed with exit code " + to_string(ret) +
                            (message.empty() ? "" : ": " + message));
    }
    unlink(err_path.c_str());

    // GAF text back to Alignments. This is why a C shim would change so little: it
    // would replace everything above and feed this same parse.
    size_t parsed = 0;
    vg::io::gaf_unpaired_for_each(graph, state.gaf_path, [&](Alignment& aln) {
        ++parsed;
        if (!passes_filter(aln)) {
            return;
        }
        count_fetched();
        iteratee(aln);
    });
    return parsed;
}

void GafBaseSiteReadSource::fetch_span(nid_t min_id, nid_t max_id,
                                       const function<void(const Alignment&)>& iteratee) const {
    // Ask only for node IDs that exist. The window is a range of IDs, but ID space is
    // not dense, and gbz-base is entitled to complain about a node that is not there.
    vector<nid_t> nodes;
    for (nid_t id = max<nid_t>(1, min_id); id <= max_id; ++id) {
        if (graph.has_node(id)) {
            nodes.push_back(id);
        }
    }
    if (nodes.empty()) {
        return;
    }

    ThreadState& state = thread_state();

    if (nodes.size() <= max_query_nodes) {
        run_query_or_die(state, nodes, iteratee);
        return;
    }

    // Too many nodes for one argv, so split. Chunks overlap in reads rather than in
    // nodes -- a read spanning a chunk boundary comes back from both -- so duplicates
    // have to be dropped.
    //
    // De-duplicate on name *and start position*, not name alone. Paired reads share a
    // name: in real Illumina GAF both mates carry the same identifier, so keying on the
    // name would silently discard one mate of every pair that reached this path,
    // halving the evidence at those sites. Two records sharing a name and a start
    // position are genuinely the same alignment returned twice.
    unordered_set<string> seen;
    for (size_t start = 0; start < nodes.size(); start += max_query_nodes) {
        size_t end = min(start + max_query_nodes, nodes.size());
        vector<nid_t> chunk(nodes.begin() + start, nodes.begin() + end);
        run_query_or_die(state, chunk, [&](const Alignment& aln) {
            string key = aln.name();
            if (aln.path().mapping_size() > 0) {
                const Position& pos = aln.path().mapping(0).position();
                key += "\t" + to_string(pos.node_id()) + "\t" + to_string(pos.offset()) +
                       (pos.is_reverse() ? "-" : "+");
            }
            if (seen.insert(std::move(key)).second) {
                iteratee(aln);
            }
        });
    }
}

void GafBaseSiteReadSource::run_query_or_die(ThreadState& state, const vector<nid_t>& nodes,
                                             const function<void(const Alignment&)>& iteratee) const {
    // Calling happens inside an OpenMP parallel region, and an exception must not
    // propagate out of one -- that is undefined behaviour, not a clean error. A failed
    // query also means we cannot score this site correctly, and carrying on with the
    // reads we did get would silently produce a wrong genotype. So report and stop.
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
