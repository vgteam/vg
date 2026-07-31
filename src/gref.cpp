#include "gref.hpp"
#include "path.hpp"
#include <cassert>
#include <sstream>
#include <algorithm>
#include <functional>
#include <queue>
#include <limits>
#include <iomanip>

//#define debug

namespace vg {

using namespace std;

const string GrefCover::gref_prefix = "gref_";
const string GrefCover::gref_suffix = "_alt";

bool GrefCover::is_gref_derived(const string& name) {
    return name.compare(0, gref_prefix.length(), gref_prefix) == 0;
}


string GrefCover::make_gref_name(const string& base_path_name, int64_t gref_index) {
    // New naming convention: {base}_{N}_alt
    return base_path_name + "_" + to_string(gref_index) + gref_suffix;
}

bool GrefCover::is_gref_name(const string& path_name) {
    // Check for pattern: _{digits}_alt at the end
    // Must end with "_alt"
    if (path_name.length() < 7) {  // minimum: "x_1_alt" (7 chars)
        return false;
    }
    if (path_name.substr(path_name.length() - 4) != gref_suffix) {
        return false;
    }
    // Find the underscore before the digits
    size_t alt_pos = path_name.length() - 4;  // position of "_alt"
    size_t underscore_pos = path_name.rfind('_', alt_pos - 1);
    if (underscore_pos == string::npos || underscore_pos == alt_pos - 1) {
        return false;  // no underscore found, or nothing between underscore and _alt
    }
    // Check that everything between underscore and _alt is digits
    for (size_t i = underscore_pos + 1; i < alt_pos; ++i) {
        if (!isdigit(path_name[i])) {
            return false;
        }
    }
    return true;
}

string GrefCover::parse_base_path(const string& gref_name) {
    if (!is_gref_name(gref_name)) {
        return gref_name;
    }
    // Find _{N}_alt and strip it
    size_t alt_pos = gref_name.length() - 4;  // position of "_alt"
    size_t underscore_pos = gref_name.rfind('_', alt_pos - 1);
    return gref_name.substr(0, underscore_pos);
}

int64_t GrefCover::parse_gref_index(const string& gref_name) {
    if (!is_gref_name(gref_name)) {
        return -1;
    }
    // Extract N from _{N}_alt
    size_t alt_pos = gref_name.length() - 4;  // position of "_alt"
    size_t underscore_pos = gref_name.rfind('_', alt_pos - 1);
    return stoll(gref_name.substr(underscore_pos + 1, alt_pos - underscore_pos - 1));
}

void GrefCover::set_verbose(bool verbose) {
    this->verbose = verbose;
}


void GrefCover::clear(MutablePathMutableHandleGraph* graph) {
    // Everything in the gref namespace goes, copied base contigs included: a
    // recomputed cover replaces all of it, and leaving stale copies behind would
    // let them outlive the reference they were copied from.
    vector<path_handle_t> gref_paths_to_remove;
    graph->for_each_path_handle([&](path_handle_t path_handle) {
        if (is_gref_derived(graph->get_path_name(path_handle))) {
            gref_paths_to_remove.push_back(path_handle);
        }
    });
    for (path_handle_t path_handle : gref_paths_to_remove) {
        graph->destroy_path(path_handle);
    }
}

void GrefCover::compute(const PathHandleGraph* graph,
                            const bdsg::SnarlDistanceIndex* distance_index,
                            const unordered_set<path_handle_t>& reference_paths,
                            int64_t minimum_length) {

    // start from scratch
    this->gref_intervals.clear();
    this->interval_snarl_bounds.clear();
    this->node_to_interval.clear();
    this->graph = graph;
    // Rank by name only (ignoring coverage) produces fewer, longer intervals
    // in practice: adjacent snarls tend to pick the same path, so same-path
    // merging succeeds more often during the fold.
    // start with the reference paths
    for (const path_handle_t& ref_path_handle : reference_paths) {
        // A reference contig whose own name ends in _{N}_alt puts its gref copy into the same
        // name space as the fragments: copy_base_paths_to_gref() would create
        // gref_CHM13#0#chr1_1_alt as a base copy while a fragment off CHM13#0#chr1 wants that
        // same name.  Refuse rather than let the two collide, which is what the per-base index
        // scans in apply() and write_gref_segments() used to paper over -- and they papered
        // over it inconsistently, because one ran before the base copies existed and one after,
        // so the segments table named a path holding entirely different sequence.
        // Test the copy name, not the path name: a PanSN path read from a GFA with no RS
        // header carries a phase block (CHM13#0#chr1_1_alt#0), which does not end in _alt even
        // though the copy made from it does.
        string ref_name = graph->get_path_name(ref_path_handle);
        string copy_name = make_gref_copy_name(ref_name);
        if (is_gref_name(copy_name)) {
            cerr << "[gref error]: reference path " << ref_name << " would be copied to "
                 << copy_name << ", which is a gref fragment name (_{N}_alt), so it would"
                 << " collide with the fragments hanging off "
                 << parse_base_path(copy_name) << ". Rename the contig, or select reference"
                 << " paths that are not in the gref namespace." << endl;
            exit(1);
        }
        this->gref_intervals.push_back(make_pair(graph->path_begin(ref_path_handle),
                                                    graph->path_end(ref_path_handle)));
        this->interval_snarl_bounds.push_back({0, 0});
        graph->for_each_step_in_path(ref_path_handle, [&](step_handle_t step_handle) {
            nid_t node_id = graph->get_id(graph->get_handle_of_step(step_handle));
            if (node_to_interval.count(node_id)) {
                cerr << "[gref error]: node " << node_id << " covered by two reference paths,"
                     << " including " << graph->get_path_name(ref_path_handle) << " and "
                     << graph->get_path_name(graph->get_path_handle_of_step(gref_intervals.at(node_to_interval.at(node_id)).first))
                     << ". Graph reference path support currently requires disjoint acyclic reference paths" << endl;
                exit(1);
            }
            node_to_interval[node_id] = gref_intervals.size() - 1;
        });
    }
    this->num_ref_intervals = this->gref_intervals.size();

#ifdef debug
#pragma omp critical(cerr)
    cerr << "[gref] Selected " << gref_intervals.size() << " rank=0 reference paths" << endl;
#endif

    // Reference nodes are now claimed, which is what makes the anchoring test meaningful, and
    // nothing has been covered yet, which is what makes withholding possible.
    this->excluded_nodes.clear();
    this->excluded_units = 0;
    this->excluded_bp = 0;
    if (distance_index != nullptr) {
        enforce_top_level_anchoring(*distance_index);
    }

    // The whole cover comes from fill_uncovered_nodes(): maximal runs of uncovered nodes
    // along each source path.  There is no snarl-traversal phase and no fold, because there
    // is nothing for them to add.
    //
    // A snarl traversal is bounded by its bubble, so it could never offer a candidate longer
    // than one snarl however far the haplotype actually agreed with the graph, and the cover
    // then had to reassemble the pieces afterwards -- which only worked when the same source
    // path happened to win every adjacent snarl.  That reassembly is what broke when ranking
    // changed in 8c42f5b54 (+36% fragments genome-wide).  Whole-path runs are maximal by
    // construction, so nothing needs reassembling.
    //
    // That also makes the candidate set complete rather than merely large: every admissible
    // fragment is contiguous on one source path, single-orientation, free of repeated node
    // ids and free of reference nodes, which is exactly what a maximal run is maximal under.
    // So every admissible fragment is a sub-run of exactly one enumerated run, and no other
    // candidate generator can beat it.

    // Build the cover.  minimum_length is deferred: a run that is too short on its own may
    // still be the longest thing available for its nodes, and filtering before selection
    // would free those nodes back into competition for no gain.
    fill_uncovered_nodes(1);

    // debug: verify all nodes are covered
    verify_cover(minimum_length);

    filter_short_intervals(minimum_length);

    // the cover is final here: nothing below adds, merges or extends an interval
    verify_disjoint();

    // With the cover final, record where each fragment sits in the top-level decomposition,
    // and how deeply it is nested inside other fragments.
    if (distance_index != nullptr) {
        assign_top_level_snarls(*distance_index);
        assign_nesting(*distance_index);
    }

    if (verbose) {
        int64_t final_alt = gref_intervals.size() - num_ref_intervals;
        cerr << "[gref] After length filter (min " << minimum_length << " bp): "
             << final_alt << " alt intervals" << endl;
    }
}

void GrefCover::fill_uncovered_nodes(int64_t minimum_length) {
    // Collect all uncovered nodes and the paths that pass through them
    unordered_set<nid_t> uncovered_nodes;
    map<string, path_handle_t> candidate_paths;  // sorted by name for deterministic ordering

    graph->for_each_handle([&](handle_t handle) {
        nid_t node_id = graph->get_id(handle);
        // Withheld nodes are never offered as candidates, so no run can enter one and no
        // fragment can be built from one.  See enforce_top_level_anchoring().
        if (!node_to_interval.count(node_id) && !this->excluded_nodes.count(node_id)) {
            uncovered_nodes.insert(node_id);
            graph->for_each_step_on_handle(handle, [&](step_handle_t step) {
                path_handle_t path_handle = graph->get_path_handle_of_step(step);
                string path_name = graph->get_path_name(path_handle);
                // Skip existing gref paths: covering a fragment with a copy of itself would
                // make the cover depend on whether one was already there.
                if (!is_gref_name(path_name)) {
                    candidate_paths[path_name] = path_handle;
                }
                return true;
            });
        }
    });

    if (uncovered_nodes.empty()) {
        return;
    }

    // Enumerate every maximal uncovered run on every candidate path, then select from all of
    // them at once, longest first.
    //
    // Walking the paths in name order and letting each claim greedily as it goes -- which is
    // what this used to do -- means the first path to reach a node wins it regardless of how
    // much either path could have covered in one piece.  With snarl-bounded candidates that
    // hardly mattered, because a run could not exceed its bubble anyway.  With whole-path
    // runs it decides the answer: a path that offers 40 kb in one piece has to beat one that
    // offers 200 bp of the same nodes, and only a global ordering can express that.
    struct FillRun {
        int64_t length;
        bool reverse;
        const string* name;
        step_handle_t start;
        step_handle_t end;   // one past the last step
        bool operator<(const FillRun& other) const {   // max-heap: "worse than"
            if (this->length != other.length) {
                return this->length < other.length;
            }
            // A run a path walks forwards needs no flipping to emit; prefer it, but only
            // between runs of equal length, never at the cost of covering less in one piece.
            if (this->reverse != other.reverse) {
                return this->reverse;
            }
            return *this->name > *other.name;   // flipped: low name wins
        }
    };

    // Split one path into maximal runs of currently-uncovered nodes, breaking at orientation
    // flips and at a repeated node id so no run can contain a cycle.
    auto enumerate_runs = [&](const string& path_name, path_handle_t path_handle,
                              vector<FillRun>& out) {
        bool in_run = false;
        step_handle_t run_start;
        int64_t run_length = 0;
        bool run_reverse = false;
        unordered_set<nid_t> run_nodes;
        step_handle_t run_end;

        auto close_run = [&]() {
            if (in_run) {
                if (run_length >= minimum_length) {
                    out.push_back({run_length, run_reverse, &path_name, run_start, run_end});
                }
                // Guarded: this is called for every already-covered step, which on a
                // haplotype means every reference node it walks.  clear() memsets the whole
                // bucket array and the table never shrinks, so an unguarded clear costs
                // O(largest run this path has produced) on each of those steps.  run_nodes
                // only ever gains entries on the path that also sets in_run, so in_run is
                // false exactly when the set is already empty.
                run_nodes.clear();
            }
            in_run = false;
            run_length = 0;
        };

        graph->for_each_step_in_path(path_handle, [&](step_handle_t step) {
            handle_t handle = graph->get_handle_of_step(step);
            nid_t node_id = graph->get_id(handle);
            bool is_reverse = graph->get_is_reverse(handle);
            if (!uncovered_nodes.count(node_id) ||
                (in_run && (run_nodes.count(node_id) || is_reverse != run_reverse))) {
                close_run();
            }
            if (uncovered_nodes.count(node_id)) {
                if (!in_run) {
                    in_run = true;
                    run_start = step;
                    run_reverse = is_reverse;
                }
                run_end = graph->get_next_step(step);
                run_length += graph->get_length(handle);
                run_nodes.insert(node_id);
            }
            return true;
        });
        close_run();
    };

    vector<FillRun> runs;
    for (const auto& name_path : candidate_paths) {
        enumerate_runs(name_path.first, name_path.second, runs);
    }
    std::make_heap(runs.begin(), runs.end());

    // Claim longest-first.  A run may have been partly taken since it was pushed, so verify
    // before claiming and re-push whatever is still free -- the same lazy-priority-queue
    // pattern compute_snarl() uses.
    while (!runs.empty()) {
        FillRun best = runs.front();
        std::pop_heap(runs.begin(), runs.end());
        runs.pop_back();

        // Re-split the run against what is still uncovered.
        bool chopped = false;
        bool in_piece = false;
        step_handle_t piece_start;
        step_handle_t piece_end;
        int64_t piece_length = 0;
        unordered_set<nid_t> piece_nodes;
        vector<FillRun> pieces;
        auto close_piece = [&]() {
            if (in_piece) {
                if (piece_length >= minimum_length) {
                    pieces.push_back({piece_length, best.reverse, best.name, piece_start, piece_end});
                }
                // Guarded for the same reason as close_run() above: called once per
                // already-claimed step, and an unguarded clear() is O(bucket count).
                piece_nodes.clear();
            }
            in_piece = false;
            piece_length = 0;
        };
        for (step_handle_t step = best.start; step != best.end; step = graph->get_next_step(step)) {
            handle_t handle = graph->get_handle_of_step(step);
            nid_t node_id = graph->get_id(handle);
            if (!uncovered_nodes.count(node_id) || piece_nodes.count(node_id)) {
                if (!uncovered_nodes.count(node_id)) {
                    chopped = true;
                }
                close_piece();
                if (!uncovered_nodes.count(node_id)) {
                    continue;
                }
            }
            if (!in_piece) {
                in_piece = true;
                piece_start = step;
            }
            piece_end = graph->get_next_step(step);
            piece_length += graph->get_length(handle);
            piece_nodes.insert(node_id);
        }
        close_piece();

        if (chopped) {
            // Something claimed part of it; the survivors go back and compete on their own
            // lengths rather than inheriting this run's priority.
            for (const FillRun& piece : pieces) {
                runs.push_back(piece);
                std::push_heap(runs.begin(), runs.end());
            }
            continue;
        }
        if (pieces.empty()) {
            continue;
        }

        for (const FillRun& piece : pieces) {
            add_interval(this->gref_intervals, this->node_to_interval,
                         make_pair(piece.start, piece.end),
                         &this->interval_snarl_bounds, {0, 0});
            for (step_handle_t step = piece.start; step != piece.end;
                 step = graph->get_next_step(step)) {
                uncovered_nodes.erase(graph->get_id(graph->get_handle_of_step(step)));
            }
        }
    }

#ifdef debug
#pragma omp critical(cerr)
    cerr << "[gref] fill_uncovered_nodes: " << uncovered_nodes.size() << " nodes still uncovered" << endl;
#endif
}


string GrefCover::make_gref_copy_name(const string& path_name) {
    PathSense sense;
    string sample, locus;
    size_t haplotype, phase_block;
    subrange_t subrange;
    PathMetadata::parse_path_name(path_name, sense, sample, locus, haplotype, phase_block, subrange);

    // The phase block has to go: make_gref_name() appends "_{N}_alt" to the base name
    // derived from this one, and a trailing phase block would turn SAMPLE#HAP#CONTIG#0
    // into SAMPLE#HAP#CONTIG#0_1_alt, which no longer parses as a path name at all (it
    // comes back GENERIC, with no sample and a '#' inside the locus).  The subrange
    // stays: it is what keeps subpaths of one contig distinct, and dropping it here
    // would collapse them onto one name.
    bool structured = sample != PathMetadata::NO_SAMPLE_NAME && locus != PathMetadata::NO_LOCUS_NAME;
    string copy_name = structured ?
        PathMetadata::create_path_name(PathSense::REFERENCE, sample, locus, haplotype,
                                       PathMetadata::NO_PHASE_BLOCK, subrange) :
        path_name;

    // Everything the cover writes lives in the gref namespace.  For a PanSN name the
    // prefix lands on the sample (GRCh38#0#chr1 -> gref_GRCh38#0#chr1), so the gref
    // paths form their own sample and both directions of the base <-> gref link stay
    // recoverable from the name alone.
    return gref_prefix + copy_name;
}

string GrefCover::make_gref_base_name(const string& path_name) {
    // A fragment name is this plus "_{N}_alt", which has to land on the locus, so the
    // subrange comes off here (and only here).
    subrange_t subrange;
    return Paths::strip_subrange(make_gref_copy_name(path_name), &subrange);
}

void GrefCover::apply(MutablePathMutableHandleGraph* mutable_graph) {
    if (this->graph != static_cast<const PathHandleGraph*>(mutable_graph)) {
        cerr << "[gref] error: apply() called with a different graph than compute()/load()" << endl;
        exit(1);
    }
#ifdef debug
    cerr << "applying gref cover with " << this->num_ref_intervals << " ref intervals "
         << " and " << this->gref_intervals.size() << " total intervals" << endl;
#endif

    // Copy the base reference paths into the gref namespace first, so the gref
    // reference is complete: the copies carry the reference sequence and the
    // fragments hang off them.
    {
        // Collect reference path handles from the reference intervals
        unordered_set<path_handle_t> reference_paths;
        for (int64_t i = 0; i < this->num_ref_intervals; ++i) {
            reference_paths.insert(graph->get_path_handle_of_step(gref_intervals[i].first));
        }
        copy_base_paths_to_gref(mutable_graph, reference_paths);
    }

    // Fragments are numbered from scratch.  clear() removes every gref path before the cover
    // is computed, and a reference path can no longer be named like a fragment, so there is
    // never a pre-existing index to continue from.
    base_path_gref_counter.clear();

    // write the gref paths
    int64_t written_intervals = 0;
    int64_t written_length = 0;
    int64_t skipped_intervals = 0;
    for (int64_t i = this->num_ref_intervals; i < this->gref_intervals.size(); ++i) {
        // Skip empty intervals (these can be created by defragment_intervals or merging)
        path_handle_t interval_path = graph->get_path_handle_of_step(gref_intervals[i].first);
        if (gref_intervals[i].first == graph->path_end(interval_path)) {
            skipped_intervals++;
            continue;
        }

        // Find the reference path this gref path extends from.  write_gref_segments() has
        // already predicted this name using the same call; see resolve_base_path_name().
        string base_path_name = this->resolve_base_path_name(i);

        // Check if this interval is all-reverse (needs to be flipped when writing)
        bool all_reverse = graph->get_is_reverse(graph->get_handle_of_step(gref_intervals[i].first));

        // Safety check: verify consistent orientation (should be guaranteed by upstream filtering)
        bool mixed = false;
        for (step_handle_t step_handle = gref_intervals[i].first; step_handle != gref_intervals[i].second;
             step_handle = graph->get_next_step(step_handle)) {
            if (graph->get_is_reverse(graph->get_handle_of_step(step_handle)) != all_reverse) {
                mixed = true;
                break;
            }
        }
        if (mixed) {
            // Mixed-orientation interval should not reach here; skip as safety net
            skipped_intervals++;
            continue;
        }

        // Get next available gref index for this base path
        int64_t gref_index = ++base_path_gref_counter[base_path_name];

        // Create the gref path name
        string gref_name = make_gref_name(base_path_name, gref_index);

        // Create the path as REFERENCE sense
        path_handle_t gref_handle = mutable_graph->create_path_handle(gref_name, false);

        int64_t interval_length = 0;
        if (!all_reverse) {
            // Forward interval: walk forward and append steps as-is
            for (step_handle_t step_handle = gref_intervals[i].first; step_handle != gref_intervals[i].second;
                 step_handle = mutable_graph->get_next_step(step_handle)) {
                mutable_graph->append_step(gref_handle, mutable_graph->get_handle_of_step(step_handle));
                interval_length += mutable_graph->get_length(mutable_graph->get_handle_of_step(step_handle));
            }
        } else {
            // All-reverse interval: collect handles, reverse order, flip each to forward
            vector<handle_t> handles;
            for (step_handle_t step_handle = gref_intervals[i].first; step_handle != gref_intervals[i].second;
                 step_handle = graph->get_next_step(step_handle)) {
                handles.push_back(mutable_graph->flip(mutable_graph->get_handle_of_step(step_handle)));
            }
            std::reverse(handles.begin(), handles.end());
            for (handle_t h : handles) {
                mutable_graph->append_step(gref_handle, h);
                interval_length += mutable_graph->get_length(h);
            }
        }
        written_intervals++;
        written_length += interval_length;
    }

#ifdef debug
    cerr << "[gref] apply: wrote " << written_intervals << " gref paths (" << written_length << " bp), skipped " << skipped_intervals << " empty intervals" << endl;
#endif
}












void GrefCover::add_interval(vector<pair<step_handle_t, step_handle_t>>& thread_gref_intervals,
                             unordered_map<nid_t, int64_t>& thread_node_to_interval,
                             const pair<step_handle_t, step_handle_t>& new_interval,
                             vector<pair<nid_t, nid_t>>* snarl_bounds_vec,
                             pair<nid_t, nid_t> snarl_bounds) {
    // Record the interval and give it every node it walks.  Nothing merges, extends or
    // decommissions: candidates are maximal runs of *uncovered* nodes, so a new interval
    // can never abut an existing one on the same path (the node between them would have to
    // be both covered and uncovered), and can never share a node with one.  Measured before
    // the merging code was removed: 708,097 add_interval calls on chrY produced 0 merges of
    // any kind, and 0 across all 60 nesting-fixture configurations.
    int64_t interval_idx = (int64_t)thread_gref_intervals.size();
    thread_gref_intervals.push_back(new_interval);
    if (snarl_bounds_vec != nullptr) {
        snarl_bounds_vec->push_back(snarl_bounds);
    }
    for (step_handle_t step = new_interval.first; step != new_interval.second;
         step = graph->get_next_step(step)) {
        thread_node_to_interval[graph->get_id(graph->get_handle_of_step(step))] = interval_idx;
    }
}


void GrefCover::filter_short_intervals(int64_t minimum_length) {
    if (minimum_length <= 1) {
        return;
    }
    // Build the surviving cover directly.  This used to decommission an interval by writing
    // path_end/path_front_end sentinels into it and then compact the sentinels away in a
    // second pass (defragment_intervals(), whose only caller this was).  The sentinel
    // protocol dated from when merging could retire an interval mid-computation; nothing
    // retires one now, so there is no state between the two passes worth having.
    vector<pair<step_handle_t, step_handle_t>> kept;
    vector<pair<nid_t, nid_t>> kept_bounds;
    kept.reserve(this->gref_intervals.size());
    kept_bounds.reserve(this->gref_intervals.size());

    // Reference intervals are never filtered: dropping one would shift a fragment into the
    // reference block and silently make every "idx < num_ref_intervals" test wrong.
    for (int64_t i = 0; i < this->num_ref_intervals; ++i) {
        kept.push_back(this->gref_intervals[i]);
        kept_bounds.push_back(this->interval_snarl_bounds[i]);
    }
    for (int64_t i = this->num_ref_intervals; i < (int64_t)this->gref_intervals.size(); ++i) {
        const pair<step_handle_t, step_handle_t>& interval = this->gref_intervals[i];
        int64_t length = 0;
        for (step_handle_t step = interval.first; step != interval.second;
             step = graph->get_next_step(step)) {
            length += graph->get_length(graph->get_handle_of_step(step));
        }
        if (length >= minimum_length) {
            kept.push_back(interval);
            kept_bounds.push_back(this->interval_snarl_bounds[i]);
        }
    }

    this->gref_intervals = std::move(kept);
    this->interval_snarl_bounds = std::move(kept_bounds);

    this->node_to_interval.clear();
    for (int64_t i = 0; i < (int64_t)this->gref_intervals.size(); ++i) {
        const pair<step_handle_t, step_handle_t>& interval = this->gref_intervals[i];
        for (step_handle_t step = interval.first; step != interval.second;
             step = graph->get_next_step(step)) {
            this->node_to_interval[graph->get_id(graph->get_handle_of_step(step))] = i;
        }
    }
}



string GrefCover::resolve_base_path_name(int64_t interval_index,
                                         path_handle_t* out_ref_path,
                                         nid_t* out_ref_node) const {
    // Trace back to the nearest reference node and name the fragment after the reference
    // path that node sits on.
    nid_t first_node = graph->get_id(graph->get_handle_of_step(this->gref_intervals[interval_index].first));
    vector<pair<int64_t, nid_t>> ref_nodes = this->get_reference_nodes(first_node, true);

    if (!ref_nodes.empty()) {
        nid_t ref_node_id = ref_nodes.at(0).second;
        int64_t ref_interval_idx = this->node_to_interval.at(ref_node_id);
        path_handle_t ref_path_handle =
            graph->get_path_handle_of_step(this->gref_intervals[ref_interval_idx].first);
        if (out_ref_path != nullptr) {
            *out_ref_path = ref_path_handle;
        }
        if (out_ref_node != nullptr) {
            *out_ref_node = ref_node_id;
        }
        return make_gref_base_name(graph->get_path_name(ref_path_handle));
    }

    // Nothing to trace back to: this fragment is in a component with no reference path in
    // it, which get_reference_nodes() notes is expected after clipping.  Name it after the
    // path it was taken from instead, and leave the out-parameters alone.
    path_handle_t source_path_handle =
        graph->get_path_handle_of_step(this->gref_intervals[interval_index].first);
    return make_gref_base_name(graph->get_path_name(source_path_handle));
}

vector<pair<int64_t, nid_t>> GrefCover::get_reference_nodes(nid_t node_id, bool first) const {

    // search back to reference in order to find the rank.
    unordered_set<nid_t> visited;
    // Min-heap by distance: explore nearest nodes first to find shortest path to reference
    priority_queue<pair<int64_t, nid_t>, vector<pair<int64_t, nid_t>>,
                   std::greater<pair<int64_t, nid_t>>> queue;
    queue.push(make_pair(0, node_id));

    nid_t current_id;
    int64_t distance = 0;

    // output reference intervals
    vector<pair<int64_t, nid_t>> output_reference_nodes;

    while (!queue.empty()) {
        std::tie(distance, current_id) = queue.top();
        queue.pop();

        if (!visited.count(current_id)) {

            visited.insert(current_id);

            if (this->node_to_interval.count(current_id)) {
                int64_t interval_idx = this->node_to_interval.at(current_id);

                const pair<step_handle_t, step_handle_t>& gref_interval = this->gref_intervals.at(interval_idx);

                // we've hit the reference, fish out its step and stop searching.
                if (interval_idx < this->num_ref_intervals) {
                    output_reference_nodes.push_back(make_pair(distance, current_id));
                    if (first) {
                        break;
                    }
                    continue;
                }

                // search out of the snarl -- any parent traversals will overlap here
                graph->follow_edges(graph->get_handle_of_step(gref_interval.first), true, [&](handle_t prev) {
                    queue.push(make_pair(distance + 1, graph->get_id(prev)));
                });
                // hack around gbwtgraph bug (feature?) that does not let you decrement path_end
                path_handle_t path_handle = graph->get_path_handle_of_step(gref_interval.first);
                step_handle_t last_step;
                if (gref_interval.second == graph->path_end(path_handle)) {
                    last_step = graph->path_back(path_handle);
                } else {
                    last_step = graph->get_previous_step(gref_interval.second);
                }
                graph->follow_edges(graph->get_handle_of_step(last_step), false, [&](handle_t next) {
                    queue.push(make_pair(distance + 1, graph->get_id(next)));
                });

            } else {
                // revert to graph search if node not in interval (distance doesn't increase -- we only count intervals)
                graph->follow_edges(graph->get_handle(current_id), false, [&](handle_t next) {
                    queue.push(make_pair(distance, graph->get_id(next)));
                });
                graph->follow_edges(graph->get_handle(current_id), true, [&](handle_t next) {
                    queue.push(make_pair(distance, graph->get_id(next)));
                });
            }

        }
    }

    // Note: output_reference_nodes may be empty if the node is in a disconnected
    // component that cannot trace back to any reference interval (e.g., after clipping)
    return output_reference_nodes;
}

void GrefCover::verify_disjoint() const {
    // One owner per node, established by walking the intervals rather than by trusting
    // node_to_interval -- see the header comment for why that map cannot detect this.
    unordered_map<nid_t, int64_t> owner;
    owner.reserve(this->node_to_interval.size());

    int64_t shared_nodes = 0;
    int64_t repeated_nodes = 0;
    // Report a bounded sample: on a whole chromosome an unbounded list would bury the count.
    const int64_t max_reported = 10;
    vector<string> reports;

    auto describe = [&](int64_t interval_idx) {
        const pair<step_handle_t, step_handle_t>& interval = this->gref_intervals.at(interval_idx);
        path_handle_t path = graph->get_path_handle_of_step(interval.first);
        string kind = interval_idx < this->num_ref_intervals ? "reference" : "fragment";
        return kind + " interval " + std::to_string(interval_idx) + " on " + graph->get_path_name(path);
    };

    for (int64_t i = 0; i < (int64_t)this->gref_intervals.size(); ++i) {
        const pair<step_handle_t, step_handle_t>& interval = this->gref_intervals[i];
        path_handle_t interval_path = graph->get_path_handle_of_step(interval.first);
        // Decommissioned intervals are empty and are skipped by apply(); skip them here too.
        if (interval.first == graph->path_end(interval_path)) {
            continue;
        }
        unordered_set<nid_t> seen_in_interval;
        for (step_handle_t step = interval.first; step != interval.second;
             step = graph->get_next_step(step)) {
            nid_t node_id = graph->get_id(graph->get_handle_of_step(step));

            // A3: no interval contains the same node twice.
            if (!seen_in_interval.insert(node_id).second) {
                ++repeated_nodes;
                if ((int64_t)reports.size() < max_reported) {
                    reports.push_back("node " + std::to_string(node_id) + " appears twice in "
                                      + describe(i));
                }
                continue;
            }

            auto found = owner.find(node_id);
            if (found == owner.end()) {
                owner.emplace(node_id, i);
            } else {
                ++shared_nodes;
                if ((int64_t)reports.size() < max_reported) {
                    reports.push_back("node " + std::to_string(node_id) + " is claimed by both "
                                      + describe(found->second) + " and " + describe(i));
                }
            }
        }
    }

    if (shared_nodes > 0 || repeated_nodes > 0) {
        cerr << "[gref error]: the cover is not node-disjoint: "
             << shared_nodes << " nodes claimed by more than one gref interval, "
             << repeated_nodes << " nodes repeated within one interval" << endl;
        for (const string& report : reports) {
            cerr << "[gref error]:   " << report << endl;
        }
        if (shared_nodes + repeated_nodes > (int64_t)reports.size()) {
            cerr << "[gref error]:   ... and "
                 << (shared_nodes + repeated_nodes - (int64_t)reports.size())
                 << " more" << endl;
        }
        exit(1);
    }

    if (verbose) {
        cerr << "[gref] Node-disjointness verified: " << owner.size()
             << " nodes, each claimed by exactly one of " << this->gref_intervals.size()
             << " intervals" << endl;
    }
}



void GrefCover::enforce_top_level_anchoring(const bdsg::SnarlDistanceIndex& distance_index) {
    using handlegraph::net_handle_t;
    this->excluded_nodes.clear();
    this->excluded_units = 0;
    this->excluded_bp = 0;

    auto interval_of = [&](nid_t node_id) {
        auto found = this->node_to_interval.find(node_id);
        return found == this->node_to_interval.end() ? (int64_t)-1 : found->second;
    };
    // Reference intervals are gref_intervals[0, num_ref_intervals), one per reference path,
    // so the interval index doubles as the identity of the reference contig.
    auto reference_contig_of = [&](nid_t node_id) {
        int64_t idx = interval_of(node_id);
        return (idx >= 0 && idx < this->num_ref_intervals) ? idx : (int64_t)-1;
    };

    // Withhold every node beneath net.  Reference nodes are already claimed and can never be
    // taken by a fragment, so excluding them would be meaningless; only alt sequence counts.
    std::function<void(const net_handle_t&)> withhold = [&](const net_handle_t& net) {
        if (distance_index.is_node(net)) {
            nid_t node_id = distance_index.node_id(net);
            if (reference_contig_of(node_id) < 0 && this->excluded_nodes.insert(node_id).second) {
                this->excluded_bp += graph->get_length(graph->get_handle(node_id));
            }
            return;
        }
        distance_index.for_each_child(net, [&](net_handle_t child) {
            withhold(child);
            return true;
        });
    };

    int64_t snarls_seen = 0;
    int64_t spine_nodes_seen = 0;
    int64_t non_reference_spine_nodes = 0;
    int64_t bridging_snarls = 0;
    vector<string> failures;
    const int64_t max_reported_failures = 5;
    net_handle_t root = distance_index.get_root();
    distance_index.for_each_child(root, [&](net_handle_t top_chain) {
        distance_index.for_each_child(top_chain, [&](net_handle_t child) {
            if (distance_index.is_snarl(child)) {
                ++snarls_seen;
                // The bounds are sentinels owned by the parent chain, so node_id() is valid
                // on them even though is_node() is not.
                nid_t start_id = distance_index.node_id(distance_index.get_bound(child, false, true));
                nid_t end_id = distance_index.node_id(distance_index.get_bound(child, true, false));
                int64_t start_contig = reference_contig_of(start_id);
                int64_t end_contig = reference_contig_of(end_id);
                if (start_contig < 0 || end_contig < 0 || start_contig != end_contig) {
                    ++this->excluded_units;
                    // Name what failed and how.  "Bridges two contigs" and "has a bound off
                    // the reference" are different problems with different causes, and a bare
                    // count cannot be acted on.
                    if ((int64_t)failures.size() < max_reported_failures) {
                        auto describe = [&](nid_t node_id, int64_t contig) {
                            if (contig < 0) {
                                return string("node ") + std::to_string(node_id) + " (off-reference)";
                            }
                            path_handle_t ref_path = graph->get_path_handle_of_step(
                                this->gref_intervals.at(contig).first);
                            return graph->get_path_name(ref_path);
                        };
                        failures.push_back(describe(start_id, start_contig) + " .. "
                                           + describe(end_id, end_contig));
                    }
                    if (start_contig >= 0 && end_contig >= 0) {
                        ++bridging_snarls;
                    }
                    withhold(child);
                }
            } else if (distance_index.is_node(child)) {
                ++spine_nodes_seen;
                // A bare chain node is linear spine sequence, not a site.  A non-reference one
                // does mean the spine is discontinuous here, but it costs nothing: a run of
                // spine nodes is covered by a single fragment, so there is no second fragment
                // to merge it with.  What has to be excluded is a snarl a fragment could cross
                // INTO, and that is handled above.  Measured: dangling_node, hap_extends_ref_-
                // start/end and unpathed_node each carry a non-reference spine node and
                // produce zero fragments spanning two units.  Withholding these would drop
                // sequence past a reference contig end for no gain.
                if (reference_contig_of(distance_index.node_id(child)) < 0) {
                    ++non_reference_spine_nodes;
                }
            } else {
                // Neither snarl nor node: not observed on any graph tested.  Withhold it
                // rather than assume it is safe.
                ++this->excluded_units;
                withhold(child);
            }
            return true;
        });
        return true;
    });

    if (this->excluded_units > 0) {
        cerr << "[gref] warning: " << this->excluded_units << " of " << snarls_seen
             << " top-level snarls do not have both bounds on one reference contig;"
             << " withholding " << this->excluded_nodes.size() << " nodes ("
             << this->excluded_bp << " bp) from the cover." << endl;
        cerr << "[gref] warning:   " << bridging_snarls << " bridge two reference contigs, "
             << (this->excluded_units - bridging_snarls)
             << " have a bound off the reference" << endl;
        for (const string& failure : failures) {
            cerr << "[gref] warning:   " << failure << endl;
        }
        if (this->excluded_units > (int64_t)failures.size()) {
            cerr << "[gref] warning:   ... and "
                 << (this->excluded_units - (int64_t)failures.size()) << " more" << endl;
        }
    } else if (verbose) {
        cerr << "[gref] Top-level anchoring: all " << snarls_seen
             << " snarls have both bounds on one reference contig ("
             << spine_nodes_seen << " spine nodes, " << non_reference_spine_nodes
             << " off-reference)" << endl;
    }
}

void GrefCover::assign_nesting(const bdsg::SnarlDistanceIndex& distance_index) {
    using handlegraph::net_handle_t;
    int64_t n = (int64_t)this->gref_intervals.size();
    this->interval_parent.assign(n, -1);
    this->interval_level.assign(n, 0);

    auto owner_of = [&](nid_t node_id) {
        auto found = this->node_to_interval.find(node_id);
        return found == this->node_to_interval.end() ? (int64_t)-1 : found->second;
    };

    // Each fragment's depth in the snarl tree: the shallowest depth any of its nodes sits at.
    // A fragment is a run of nodes and they need not all be at one depth, so taking the first
    // node's depth as the fragment's would be arbitrary -- and it was: doing that produced a
    // parent chain with a cycle on chr22, because a fragment could resolve to a parent that
    // was really nested inside it. Depth is what makes ancestry well founded, so it is
    // computed for the whole fragment and a parent must be strictly shallower.
    vector<int64_t> depth(n, std::numeric_limits<int64_t>::max());
    for (int64_t i = this->num_ref_intervals; i < n; ++i) {
        const pair<step_handle_t, step_handle_t>& interval = this->gref_intervals[i];
        path_handle_t interval_path = graph->get_path_handle_of_step(interval.first);
        if (interval.first == graph->path_end(interval_path)) {
            continue;
        }
        for (step_handle_t step = interval.first; step != interval.second;
             step = graph->get_next_step(step)) {
            nid_t node_id = graph->get_id(graph->get_handle_of_step(step));
            int64_t node_depth =
                (int64_t)distance_index.get_depth(distance_index.get_node_net_handle(node_id));
            depth[i] = min(depth[i], node_depth);
        }
    }

    // Pass 1: each fragment's parent.  Climb until a snarl whose bounds are owned by a gref
    // interval that is strictly shallower -- that interval's path is the one a caller reports
    // this fragment's enclosing site against.  Reference intervals are shallower than any
    // fragment by construction, so they always qualify.
    for (int64_t i = this->num_ref_intervals; i < n; ++i) {
        const pair<step_handle_t, step_handle_t>& interval = this->gref_intervals[i];
        path_handle_t interval_path = graph->get_path_handle_of_step(interval.first);
        if (interval.first == graph->path_end(interval_path)) {
            continue;
        }
        path_handle_t my_source = graph->get_path_handle_of_step(interval.first);
        auto qualifies = [&](int64_t owner) {
            if (owner < 0 || owner == i) {
                return false;
            }
            if (owner < this->num_ref_intervals) {
                return true;
            }
            // Two runs of the same haplotype are separated by sequence one of them does not
            // hold, so neither can contain the other.  Measured on chr22: 10 of 48 wrong
            // parents were the fragment's own source path, and 0 of 127 correct ones were.
            if (graph->get_path_handle_of_step(this->gref_intervals[owner].first) == my_source) {
                return false;
            }
            return depth[owner] < depth[i];
        };
        // A snarl counts as lying inside a fragment only when that fragment owns BOTH of its
        // boundary nodes.  Owning one is satisfied by a sibling allele that merely touches the
        // boundary, and a sibling is not an ancestor: on chr22 all 48 wrong parents shared a
        // top-level snarl with the fragment they were wrongly given.
        auto owner_of_snarl = [&](const net_handle_t& snarl) -> int64_t {
            nid_t s = distance_index.node_id(distance_index.get_bound(snarl, false, true));
            nid_t e = distance_index.node_id(distance_index.get_bound(snarl, true, false));
            int64_t so = owner_of(s), eo = owner_of(e);
            return (so == eo && qualifies(so)) ? so : (int64_t)-1;
        };

        // Start from the smallest net_handle containing the WHOLE fragment, not from its
        // first node.  A fragment is an allele of one site, and that site is the common
        // ancestor of all its nodes.  Climbing from the first node instead finds the first
        // ancestor bounded by any owned node, which can be a SIBLING allele that happens to
        // touch that boundary -- measured on chr22, that made chr22_1002_alt a child of
        // chr22_1756_alt when the VCF has both at CH=1 with the same PS, two alleles of one
        // site.
        nid_t first_node = graph->get_id(graph->get_handle_of_step(interval.first));
        vector<net_handle_t> lineage;
        unordered_map<size_t, int64_t> lineage_pos;
        for (net_handle_t up = distance_index.get_node_net_handle(first_node); ;
             up = distance_index.get_parent(up)) {
            lineage_pos[handlegraph::as_integer(up)] = (int64_t)lineage.size();
            lineage.push_back(up);
            if (distance_index.is_root(up)) {
                break;
            }
        }
        int64_t common = 0;
        for (step_handle_t step = interval.first; step != interval.second;
             step = graph->get_next_step(step)) {
            net_handle_t up = distance_index.get_node_net_handle(
                graph->get_id(graph->get_handle_of_step(step)));
            while (lineage_pos.find(handlegraph::as_integer(up)) == lineage_pos.end()) {
                if (distance_index.is_root(up)) {
                    break;
                }
                up = distance_index.get_parent(up);
            }
            auto found = lineage_pos.find(handlegraph::as_integer(up));
            if (found != lineage_pos.end()) {
                common = max(common, found->second);
            }
        }

        net_handle_t net = lineage[common];
        // The common ancestor may itself be the site.
        if (distance_index.is_snarl(net) && !distance_index.is_root(net)) {
            int64_t owner = owner_of_snarl(net);
            if (owner >= 0) {
                this->interval_parent[i] = owner;
                continue;
            }
        }
        while (!distance_index.is_root(net)) {
            net_handle_t up = distance_index.get_parent(net);
            if (distance_index.is_root(up)) {
                break;
            }
            if (distance_index.is_snarl(up)) {
                int64_t owner = owner_of_snarl(up);
                if (owner >= 0) {
                    this->interval_parent[i] = owner;
                    break;
                }
            }
            net = up;
        }
    }

    // Pass 2: levels, from the parent chain.  Memoised, and guarded against a cycle -- the
    // snarl tree is a tree so one cannot arise, but a wrong parent would otherwise hang here
    // rather than produce a visibly wrong number.
    const int64_t IN_PROGRESS = -1;
    std::function<int64_t(int64_t)> level_of = [&](int64_t i) -> int64_t {
        if (i < this->num_ref_intervals) {
            return 0;
        }
        if (this->interval_level[i] == IN_PROGRESS) {
            cerr << "[gref error]: nesting cycle at interval " << i
                 << "; the parent chain is not a tree" << endl;
            exit(1);
        }
        if (this->interval_level[i] != 0) {
            return this->interval_level[i];
        }
        this->interval_level[i] = IN_PROGRESS;
        int64_t parent = this->interval_parent[i];
        int64_t level = parent < 0 ? 1 : level_of(parent) + 1;
        this->interval_level[i] = level;
        return level;
    };
    for (int64_t i = this->num_ref_intervals; i < n; ++i) {
        level_of(i);
    }

    if (verbose) {
        map<int64_t, int64_t> histogram;
        for (int64_t i = this->num_ref_intervals; i < n; ++i) {
            ++histogram[this->interval_level[i]];
        }
        cerr << "[gref] Nesting levels:";
        for (const auto& level_count : histogram) {
            cerr << " " << level_count.first << ":" << level_count.second;
        }
        cerr << endl;
    }
}

void GrefCover::assign_top_level_snarls(const bdsg::SnarlDistanceIndex& distance_index) {
    // The unit of the top-level decomposition is a child of a top-level chain: either a
    // top-level snarl or a bare node of the chain (the spine).  Every node of the graph is
    // in exactly one unit.  A fragment can be processed independently of every other
    // fragment exactly when all of its nodes are in one unit -- that is the property this
    // measures, and the containing snarl it records is only meaningful when it holds.
    //
    // Spine nodes get one unit each rather than sharing a "spine" bucket, because two
    // distinct spine nodes are as separable as two snarls: whatever lies between them on
    // the chain is a top-level snarl a fragment would have had to cross to touch both.
    using handlegraph::net_handle_t;
    auto start_time = std::chrono::steady_clock::now();
    TopLevelSnarlStats stats;

    // unit >= 0 indexes snarl_bounds; unit <= -2 is a spine node.
    vector<pair<nid_t, nid_t>> snarl_bounds;
    int64_t next_spine_unit = -2;

    // Node -> unit, kept only for the nodes a fragment owns.  Reference nodes are
    // pre-assigned and verify_disjoint() forbids a fragment from holding one, so skipping
    // them loses nothing and holds the map to the size of the alt cover rather than the
    // graph (388k entries instead of 2.2M on chrY).
    unordered_map<nid_t, int64_t> node_unit;

    auto interval_of = [&](nid_t node_id) {
        auto found = this->node_to_interval.find(node_id);
        return found == this->node_to_interval.end() ? (int64_t)-1 : found->second;
    };
    auto is_fragment_node = [&](nid_t node_id) {
        return interval_of(node_id) >= this->num_ref_intervals;
    };

    // Record node -> unit, flagging any node two units both claim.  The decomposition is a
    // partition, so that cannot happen; if it did, every count below would be meaningless,
    // so it is counted rather than assumed.
    auto place = [&](nid_t node_id, int64_t unit) {
        auto placed = node_unit.emplace(node_id, unit);
        if (!placed.second && placed.first->second != unit) {
            ++stats.nodes_claimed_twice;
            placed.first->second = unit;
        }
    };

    // Claim every node beneath net for unit.  Top-down descent, not a per-node climb:
    // measured 10x cheaper, and it produces the node -> unit grouping directly.
    std::function<void(const net_handle_t&, int64_t)> claim =
        [&](const net_handle_t& net, int64_t unit) {
        if (distance_index.is_node(net)) {
            nid_t node_id = distance_index.node_id(net);
            if (is_fragment_node(node_id)) {
                place(node_id, unit);
            }
            return;
        }
        distance_index.for_each_child(net, [&](net_handle_t child) {
            claim(child, unit);
            return true;
        });
    };

    // One child of a top-level chain.
    auto handle_chain_child = [&](const net_handle_t& child) {
        if (distance_index.is_snarl(child)) {
            // The bounds are sentinels belonging to the parent chain, so they are spine
            // nodes: node_id() is valid on them even though is_node() is not.  That is what
            // makes "anchored" mean "both bounds are pre-assigned reference nodes", and
            // therefore uncrossable by any fragment.  Whether they are is checked and
            // reported by enforce_top_level_anchoring(), before the cover is built.
            nid_t start_id = distance_index.node_id(distance_index.get_bound(child, false, true));
            nid_t end_id = distance_index.node_id(distance_index.get_bound(child, true, false));
            int64_t unit = (int64_t)snarl_bounds.size();
            snarl_bounds.push_back({start_id, end_id});
            ++stats.top_level_snarls;
            claim(child, unit);
        } else if (distance_index.is_node(child)) {
            nid_t node_id = distance_index.node_id(child);
            ++stats.spine_nodes;
            if (is_fragment_node(node_id)) {
                place(node_id, next_spine_unit);
            }
            --next_spine_unit;
        } else {
            // Not observed on any graph tested.  Give it a unit of its own so the
            // per-fragment test stays total, and count it so it cannot pass unnoticed.
            ++stats.unexpected_chain_children;
            claim(child, next_spine_unit--);
        }
    };

    // Top-level chains are the children of the root, except that disconnected components
    // joined in the root hang one level deeper, under a root snarl.
    net_handle_t root = distance_index.get_root();
    distance_index.for_each_child(root, [&](net_handle_t child) {
        if (distance_index.is_root_snarl(child)) {
            distance_index.for_each_child(child, [&](net_handle_t grandchild) {
                if (distance_index.is_chain(grandchild)) {
                    ++stats.top_level_chains;
                    distance_index.for_each_child(grandchild, [&](net_handle_t c) {
                        handle_chain_child(c);
                        return true;
                    });
                } else {
                    handle_chain_child(grandchild);
                }
                return true;
            });
        } else if (distance_index.is_chain(child)) {
            ++stats.top_level_chains;
            distance_index.for_each_child(child, [&](net_handle_t c) {
                handle_chain_child(c);
                return true;
            });
        } else {
            handle_chain_child(child);
        }
        return true;
    });

    // Now test every fragment against the partition.
    const int64_t max_reported = 10;
    vector<string> reports;
    unordered_set<int64_t> units;
    unordered_set<int64_t> occupied_snarls;

    for (int64_t i = this->num_ref_intervals; i < (int64_t)this->gref_intervals.size(); ++i) {
        const pair<step_handle_t, step_handle_t>& interval = this->gref_intervals[i];
        this->interval_snarl_bounds[i] = {0, 0};
        if (interval.first == graph->path_end(graph->get_path_handle_of_step(interval.first))) {
            // decommissioned: apply() skips it, so it is not a fragment
            continue;
        }
        ++stats.fragments;

        units.clear();
        int64_t unmapped = 0;
        nid_t first_node = 0;
        nid_t last_node = 0;
        for (step_handle_t step = interval.first; step != interval.second;
             step = graph->get_next_step(step)) {
            nid_t node_id = graph->get_id(graph->get_handle_of_step(step));
            if (first_node == 0) {
                first_node = node_id;
            }
            last_node = node_id;
            auto found = node_unit.find(node_id);
            if (found == node_unit.end()) {
                ++unmapped;
            } else {
                units.insert(found->second);
            }
        }

        int64_t snarl_units = 0;
        int64_t only_snarl_unit = -1;
        for (int64_t unit : units) {
            if (unit >= 0) {
                ++snarl_units;
                only_snarl_unit = unit;
            }
        }
        stats.max_snarls_in_one_fragment = max(stats.max_snarls_in_one_fragment, snarl_units);
        stats.max_units_in_one_fragment = max(stats.max_units_in_one_fragment,
                                              (int64_t)units.size());
        if (unmapped > 0) {
            ++stats.fragments_with_unmapped_nodes;
        }
        if ((int64_t)units.size() > snarl_units) {
            ++stats.fragments_touching_spine;
        }
        if (snarl_units >= 2) {
            ++stats.fragments_spanning_snarls;
        }
        if ((int64_t)units.size() >= 2) {
            ++stats.fragments_spanning_units;
        }

        if (units.size() == 1 && unmapped == 0) {
            if (snarl_units == 1) {
                ++stats.fragments_in_one_snarl;
                occupied_snarls.insert(only_snarl_unit);
                // The only case where a containing top-level snarl exists.
                this->interval_snarl_bounds[i] = snarl_bounds[only_snarl_unit];
            } else {
                ++stats.fragments_on_one_spine_node;
            }
        } else if ((int64_t)reports.size() < max_reported) {
            reports.push_back("interval " + std::to_string(i) + " on "
                              + graph->get_path_name(graph->get_path_handle_of_step(interval.first))
                              + " nodes " + std::to_string(first_node) + ".." + std::to_string(last_node)
                              + " touches " + std::to_string(snarl_units) + " top-level snarl(s), "
                              + std::to_string((int64_t)units.size() - snarl_units) + " spine node(s), "
                              + std::to_string(unmapped) + " unplaced node(s)");
        }
    }

    stats.distinct_snarls_with_fragments = (int64_t)occupied_snarls.size();
    this->top_level_snarl_stats = stats;

    if (verbose) {
        double seconds = std::chrono::duration<double>(
            std::chrono::steady_clock::now() - start_time).count();
        cerr << "[gref] Mapped the cover onto the top-level decomposition in "
             << std::fixed << std::setprecision(2) << seconds << " s" << endl;
        cerr << "[gref] Top-level decomposition: " << stats.top_level_chains << " chains, "
             << stats.top_level_snarls << " top-level snarls, " << stats.spine_nodes
             << " spine nodes" << endl
             // Anchoring is reported by enforce_top_level_anchoring(), which tests the
             // condition that is actually enforced -- both bounds on ONE reference contig.
             // This used to print the weaker "has a non-reference bound" count beside it, so
             // chrOther-v2.1 showed 7 failures in one line and 2 in the next.
             << "[gref]   fragments: " << stats.fragments << " total, "
             << stats.fragments_in_one_snarl << " wholly inside one top-level snarl ("
             << stats.distinct_snarls_with_fragments << " distinct snarls), "
             << stats.fragments_on_one_spine_node << " on a single spine node" << endl
             << "[gref]   fragments spanning >1 top-level snarl: "
             << stats.fragments_spanning_snarls << " (max " << stats.max_snarls_in_one_fragment
             << " snarls in one fragment)" << endl
             << "[gref]   fragments spanning >1 unit: " << stats.fragments_spanning_units
             << " (max " << stats.max_units_in_one_fragment << " units); "
             << stats.fragments_touching_spine << " touch the spine; "
             << stats.fragments_with_unmapped_nodes << " have unplaced nodes" << endl;
        if (stats.unexpected_chain_children > 0) {
            cerr << "[gref]   warning: " << stats.unexpected_chain_children
                 << " top-level chain children were neither a snarl nor a node" << endl;
        }
        if (stats.nodes_claimed_twice > 0) {
            cerr << "[gref]   warning: " << stats.nodes_claimed_twice
                 << " nodes were placed in two different top-level units" << endl;
        }
        for (const string& report : reports) {
            cerr << "[gref]   " << report << endl;
        }
        if (stats.fragments - stats.fragments_in_one_snarl - stats.fragments_on_one_spine_node
            > (int64_t)reports.size()) {
            cerr << "[gref]   ... and "
                 << (stats.fragments - stats.fragments_in_one_snarl
                     - stats.fragments_on_one_spine_node - (int64_t)reports.size())
                 << " more not contained in a single unit" << endl;
        }
    }
}

void GrefCover::verify_cover(int64_t minimum_length) const {
    // Report sequence no fragment claims.  Only meaningful at minimum_length <= 1: above it,
    // short intervals are filtered after this runs and their nodes are expected to be
    // uncovered, so the count would be noise.
    //
    // There used to be a verbose summary here reporting node and interval totals.  It was
    // deleted rather than fixed: it ran before filter_short_intervals() and reported the
    // pre-filter interval count as though it were the result -- "169 ref + 45537 alt" on
    // chrOther-v2.1, against 92 actual fragments -- and its uncovered figure restated what
    // enforce_top_level_anchoring() had already printed three lines earlier.  The true count
    // is printed by compute() after the filter, and the useful checks live in
    // verify_disjoint() and assign_top_level_snarls().
    if (minimum_length > 1) {
        return;
    }
    int64_t uncovered_nodes = 0;
    int64_t uncovered_length = 0;
    graph->for_each_handle([&](handle_t handle) {
        if (!node_to_interval.count(graph->get_id(handle))) {
            ++uncovered_nodes;
            uncovered_length += graph->get_length(handle);
        }
    });
    if (uncovered_nodes > 0) {
        cerr << "[gref] warning: " << uncovered_nodes << " nodes ("
             << uncovered_length << " bp) not covered by gref paths" << endl;
    }
}

void GrefCover::copy_base_paths_to_gref(MutablePathMutableHandleGraph* mutable_graph,
                                            const unordered_set<path_handle_t>& reference_paths) {
    for (const path_handle_t& ref_path : reference_paths) {
        string ref_name = mutable_graph->get_path_name(ref_path);

        // Same naming rule the fragments use, minus the subrange strip, so a fragment
        // and the contig it hangs off always agree on their base name while subpaths of
        // one contig stay distinct.
        string new_name = make_gref_copy_name(ref_name);

        // Two reference paths mapping to one gref name would mean silently publishing
        // only one of them, so refuse rather than drop sequence on the floor.
        if (mutable_graph->has_path(new_name)) {
            cerr << "[gref] error: reference path " << ref_name << " maps to gref path "
                 << new_name << ", which already exists. Reference paths must have distinct"
                 << " gref names; run vg paths --compute-gref on a graph with no gref paths"
                 << " in it." << endl;
            exit(1);
        }

        // Create the new path with same sense as original
        path_handle_t new_path = mutable_graph->create_path_handle(new_name, false);

        // Copy all steps from original path to new path
        mutable_graph->for_each_step_in_path(ref_path, [&](step_handle_t step) {
            mutable_graph->append_step(new_path, mutable_graph->get_handle_of_step(step));
            return true;
        });

#ifdef debug
        cerr << "[gref] copy_base_paths_to_gref: copied " << ref_name << " -> " << new_name << endl;
#endif
    }
}

void GrefCover::write_gref_segments(ostream& os) {
    // Numbering starts from scratch, exactly as apply() does; see the note there.

    // Pre-compute reference node positions by walking all reference intervals once.
    // Maps node ID -> (ref_path_handle, offset of node start on ref path)
    unordered_map<nid_t, pair<path_handle_t, int64_t>> ref_node_positions;
    for (int64_t i = 0; i < this->num_ref_intervals; ++i) {
        const pair<step_handle_t, step_handle_t>& ref_interval = this->gref_intervals[i];
        path_handle_t ref_path_handle = graph->get_path_handle_of_step(ref_interval.first);
        int64_t offset = 0;
        for (step_handle_t step = ref_interval.first; step != ref_interval.second;
             step = graph->get_next_step(step)) {
            nid_t node_id = graph->get_id(graph->get_handle_of_step(step));
            ref_node_positions[node_id] = make_pair(ref_path_handle, offset);
            offset += graph->get_length(graph->get_handle_of_step(step));
        }
    }

    // Collect the specific steps whose offsets we need (interval boundaries only).
    // This avoids caching every step of every source path.
    unordered_set<step_handle_t> needed_steps;
    unordered_set<path_handle_t> source_paths_needed;
    for (int64_t i = this->num_ref_intervals; i < this->gref_intervals.size(); ++i) {
        const pair<step_handle_t, step_handle_t>& interval = this->gref_intervals[i];
        path_handle_t ph = graph->get_path_handle_of_step(interval.first);
        source_paths_needed.insert(ph);
        needed_steps.insert(interval.first);
        needed_steps.insert(interval.second);
    }
    // Walk each source path once, caching offsets only for needed steps.
    // Also cache path_end -> total path length for computing interval end offsets.
    unordered_map<step_handle_t, int64_t> step_offset_cache;
    // Where each source path touches the reference, in path order.  A fragment's enclosing
    // reference window is the span between the last such touch before it and the first
    // after it: the reference interval the fragment is an alternative to.  That is what the
    // snarl bounds used to approximate, derived from the path itself instead, so it is
    // available whether or not the cover was built from snarl traversals.
    unordered_map<path_handle_t, vector<pair<int64_t, nid_t>>> ref_anchors;
    for (const path_handle_t& ph : source_paths_needed) {
        int64_t offset = 0;
        vector<pair<int64_t, nid_t>>& anchors = ref_anchors[ph];
        for (step_handle_t step = graph->path_begin(ph);
             step != graph->path_end(ph);
             step = graph->get_next_step(step)) {
            if (needed_steps.count(step)) {
                step_offset_cache[step] = offset;
            }
            nid_t step_node = graph->get_id(graph->get_handle_of_step(step));
            if (ref_node_positions.count(step_node)) {
                anchors.emplace_back(offset, step_node);
            }
            offset += graph->get_length(graph->get_handle_of_step(step));
        }
        // Cache path_end sentinel with total path length
        step_handle_t path_end = graph->path_end(ph);
        if (needed_steps.count(path_end)) {
            step_offset_cache[path_end] = offset;
        }
    }

    // Cache resolved source display names per path handle to avoid
    // repeated parse_path_name/create_path_name per interval.
    // Stores (display_name, subrange_offset) per source path.
    unordered_map<path_handle_t, pair<string, int64_t>> source_display_cache;

    // Name every interval before emitting any, because a fragment's parent may sit at a
    // higher index than the fragment itself and the parent column needs its name.
    // The skip conditions must match the emit loop below exactly, or the counter -- and
    // therefore every later name -- drifts.
    vector<string> interval_name(this->gref_intervals.size());
    {
        unordered_map<string, int64_t> counter;
        for (int64_t i = this->num_ref_intervals; i < (int64_t)this->gref_intervals.size(); ++i) {
            const pair<step_handle_t, step_handle_t>& interval = this->gref_intervals[i];
            path_handle_t interval_path = graph->get_path_handle_of_step(interval.first);
            if (interval.first == graph->path_end(interval_path)) {
                continue;
            }
            bool all_reverse = graph->get_is_reverse(graph->get_handle_of_step(interval.first));
            bool mixed = false;
            for (step_handle_t step = interval.first; step != interval.second;
                 step = graph->get_next_step(step)) {
                if (graph->get_is_reverse(graph->get_handle_of_step(step)) != all_reverse) {
                    mixed = true;
                    break;
                }
            }
            if (mixed) {
                continue;
            }
            string base = this->resolve_base_path_name(i);
            interval_name[i] = make_gref_name(base, ++counter[base]);
        }
    }

    // The gref contig a fragment's coordinates are reported against: another fragment when it
    // is nested inside one, otherwise the gref copy of the reference path it hangs off.
    auto parent_name_of = [&](int64_t i) -> string {
        if (this->interval_parent.empty() || i >= (int64_t)this->interval_parent.size()) {
            return ".";
        }
        int64_t parent = this->interval_parent[i];
        if (parent < 0) {
            return ".";
        }
        if (parent < this->num_ref_intervals) {
            path_handle_t ref_path =
                graph->get_path_handle_of_step(this->gref_intervals[parent].first);
            return make_gref_copy_name(graph->get_path_name(ref_path));
        }
        return interval_name[parent].empty() ? "." : interval_name[parent];
    };

    // Write each gref interval
    for (int64_t i = this->num_ref_intervals; i < this->gref_intervals.size(); ++i) {
        const pair<step_handle_t, step_handle_t>& interval = this->gref_intervals[i];
        path_handle_t source_path_handle = graph->get_path_handle_of_step(interval.first);

        // Skip empty intervals
        if (interval.first == graph->path_end(source_path_handle)) {
            continue;
        }

        // Skip mixed-orientation intervals (must match apply() logic)
        bool all_reverse = graph->get_is_reverse(graph->get_handle_of_step(interval.first));
        bool mixed = false;
        for (step_handle_t step = interval.first; step != interval.second;
             step = graph->get_next_step(step)) {
            if (graph->get_is_reverse(graph->get_handle_of_step(step)) != all_reverse) {
                mixed = true;
                break;
            }
        }
        if (mixed) {
            continue;
        }

        // Look up pre-computed source path offsets
        int64_t source_start = step_offset_cache.at(interval.first);
        int64_t source_end = step_offset_cache.at(interval.second);

        // Name the fragment with the same call apply() will use, so the two agree by
        // construction rather than by two implementations happening to match.  Only the
        // reference *window* below varies with what information is available; the contig it
        // is measured on is always the one named in the gref path name.
        string ref_path_name = ".";
        int64_t ref_start = 0;
        int64_t ref_end = 0;
        path_handle_t ref_path_handle;
        nid_t traced_ref_node = 0;
        string base_path_name = this->resolve_base_path_name(i, &ref_path_handle, &traced_ref_node);

        if (traced_ref_node != 0) {
            ref_path_name = graph->get_path_name(ref_path_handle);

            // Narrowest honest window: the single reference node the trace landed on.  The
            // two branches below widen it when they can.
            ref_start = ref_node_positions[traced_ref_node].second;
            ref_end = ref_start + graph->get_length(graph->get_handle(traced_ref_node));

            // The window comes from this fragment's own source path: the last place it
            // touches this reference contig before the fragment and the first at or after
            // it.  Anchors on any other contig are skipped -- the window has to be measured
            // on the contig named in column 4.
            //
            // interval_snarl_bounds now holds the containing top-level snarl, and it is
            // deliberately NOT used here.  A top-level snarl is a child of the top-level
            // chain, not a bubble: on chrY the median one spans 185 kb of CHM13 and the
            // largest holds a quarter of the graph, because 77% of reference nodes live
            // *inside* top-level snarls rather than on the chain.  Measured over the 2,742
            // chrY fragments, taking the window from the snarl bounds instead of the anchors
            // moves the median window from 111 bp to 185,051 bp, the mean from 7.1 kb to
            // 425 kb, and the share of rows whose window is indistinguishable from another
            // row's from 41% to 83%.  It is a correct containment window and a useless one.
            const vector<pair<int64_t, nid_t>>& anchors = ref_anchors[source_path_handle];
            const pair<int64_t, nid_t>* left_anchor = nullptr;
            const pair<int64_t, nid_t>* right_anchor = nullptr;
            for (auto it = std::lower_bound(anchors.begin(), anchors.end(),
                                            make_pair(source_start, (nid_t)0));
                 it != anchors.begin(); ) {
                --it;
                if (ref_node_positions[it->second].first == ref_path_handle) {
                    left_anchor = &*it;
                    break;
                }
            }
            for (auto it = std::lower_bound(anchors.begin(), anchors.end(),
                                            make_pair(source_end, (nid_t)0));
                 it != anchors.end(); ++it) {
                if (ref_node_positions[it->second].first == ref_path_handle) {
                    right_anchor = &*it;
                    break;
                }
            }

            auto anchor_span = [&](const pair<int64_t, nid_t>* a) {
                const auto& pos = ref_node_positions[a->second];
                return make_pair(pos.second,
                                 pos.second + graph->get_length(graph->get_handle(a->second)));
            };
            if (left_anchor != nullptr || right_anchor != nullptr) {
                auto l = anchor_span(left_anchor != nullptr ? left_anchor : right_anchor);
                auto r = anchor_span(right_anchor != nullptr ? right_anchor : left_anchor);
                ref_start = min(l.first, r.first);
                ref_end = max(l.second, r.second);
            }
        }

        // Named in the pass above, which had to run first so parents could be looked up.
        const string& gref_name = interval_name[i];

        // Resolve source path name to full-path coordinates using cached result.
        // Parses path name once per source path to extract subrange offset and
        // strip the #0 phase block.
        auto cache_it = source_display_cache.find(source_path_handle);
        if (cache_it == source_display_cache.end()) {
            string source_path_name = graph->get_path_name(source_path_handle);
            PathSense sense;
            string sample, locus;
            size_t haplotype, phase_block;
            subrange_t subrange;
            PathMetadata::parse_path_name(source_path_name, sense, sample, locus, haplotype, phase_block, subrange);
            int64_t subrange_offset = (subrange != PathMetadata::NO_SUBRANGE) ? subrange.first : 0;
            if (phase_block == 0) {
                phase_block = PathMetadata::NO_PHASE_BLOCK;
                sense = PathSense::REFERENCE;
            }
            string display_name = PathMetadata::create_path_name(sense, sample, locus, haplotype, phase_block, PathMetadata::NO_SUBRANGE);
            cache_it = source_display_cache.emplace(source_path_handle, make_pair(std::move(display_name), subrange_offset)).first;
        }
        const string& display_source_name = cache_it->second.first;
        int64_t display_source_start = source_start + cache_it->second.second;
        int64_t display_source_end = source_end + cache_it->second.second;

        // The top-level snarl this fragment sits in, as its two boundary node ids.  Both are
        // reference nodes on one contig (enforce_top_level_anchoring), so the pair names the
        // site the fragment is an allele of, and every fragment inside that snarl carries the
        // same pair.  0 0 when no decomposition was supplied.
        //
        // Emitted in reference order.  The decomposition reports the bounds in the
        // orientation it traversed the chain, which can be either way round; ordering them by
        // where they sit on the reference makes the pair mean the same thing on every row.
        pair<nid_t, nid_t> top_snarl = this->interval_snarl_bounds[i];
        if (top_snarl.first != 0 && top_snarl.second != 0) {
            auto left = ref_node_positions.find(top_snarl.first);
            auto right = ref_node_positions.find(top_snarl.second);
            if (left != ref_node_positions.end() && right != ref_node_positions.end() &&
                left->second.first == right->second.first &&
                right->second.second < left->second.second) {
                std::swap(top_snarl.first, top_snarl.second);
            }
        }

        // Output the BED line
        os << display_source_name << "\t"
           << display_source_start << "\t"
           << display_source_end << "\t"
           << gref_name << "\t"
           << ref_path_name << "\t"
           << ref_start << "\t"
           << ref_end << "\t"
           << top_snarl.first << "\t"
           << top_snarl.second << "\t"
           << (this->interval_level.empty() ? 0 : this->interval_level[i]) << "\t"
           << parent_name_of(i) << "\n";
    }
}

}
