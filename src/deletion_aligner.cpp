/**
 * \file deletion_aligner.cpp
 *
 * Implements an aligner for global deletions
 *
 */

#include "deletion_aligner.hpp"

//#define debug_deletion_aligner

namespace vg {

DeletionAligner::DeletionAligner(int8_t gap_open, int8_t gap_extension)
    : gap_open(gap_open), gap_extension(gap_extension)
{
    
}

void DeletionAligner::align(Alignment& aln, const HandleGraph& graph) const {
    if (!aln.sequence().empty()) {
        cerr << "error: DeletionAligner can only be used for alignments of empty strings" << endl;
        exit(1);
    }
    auto traces = run_dp(graph, 1);
    trace_to_alignment(aln, traces.at(0), graph);
}

void DeletionAligner::align_multi(Alignment& aln, vector<Alignment>& alt_alns,
                                  const HandleGraph& graph, int32_t max_alt_alns) const {
    if (!aln.sequence().empty()) {
        cerr << "error: DeletionAligner can only be used for alignments of empty strings" << endl;
        exit(1);
    }
    auto traces = run_dp(graph, max_alt_alns);
    alt_alns.resize(traces.size());
    for (size_t i = 0; i < traces.size(); ++i) {
        alt_alns[i].set_sequence(aln.sequence());
        alt_alns[i].set_quality(aln.quality());
        trace_to_alignment(alt_alns[i], traces[i], graph);
    }
    *aln.mutable_path() = alt_alns.front().path();
    aln.set_score(alt_alns.front().score());
}

vector<vector<handle_t>> DeletionAligner::run_dp(const HandleGraph& graph,
                                                 int32_t max_tracebacks) const {
#ifdef debug_deletion_aligner
    cerr << "aligning deletions with " << max_tracebacks << " tracebacks" << endl;
#endif
    auto order = handlealgs::lazier_topological_order(&graph);
    
    if (order.empty() && max_tracebacks > 0) {
        // Turns out the graph is empty.
        // We need to produce one traceback of visiting nothing.
        return {{}};
    }
    
    unordered_map<handle_t, size_t> index_of;
    index_of.reserve(order.size());
    for (size_t i = 0; i < order.size(); ++i) {
        index_of[order[i]] = i;
    }
    vector<size_t> dists;
    vector<pair<size_t, size_t>> sinks;
    tie(dists, sinks) = min_dists(order, index_of, graph);
    return traceback(order, index_of, graph, dists, sinks, max_tracebacks);
}

pair<vector<size_t>, vector<pair<size_t, size_t>>> DeletionAligner::min_dists(
      const vector<handle_t>& order,
      const unordered_map<handle_t, size_t>& index_of,
      const HandleGraph& graph) const
{
#ifdef debug_deletion_aligner
    cerr << "finding min dists among " << order.size() << " handles" << endl;
#endif
    // use dynamic programming to compute the minimum distance from a source
    vector<size_t> dists(order.size(), numeric_limits<size_t>::max());
    vector<pair<size_t, size_t>> sinks;
    for (size_t i = 0; i < order.size(); ++i) {
        if (dists[i] == numeric_limits<size_t>::max()) {
            // nothing has replaced the starting value, must be a source
            dists[i] = 0;
#ifdef debug_deletion_aligner
            cerr << "Declare a source at " << graph.get_id(order[i]) << (graph.get_is_reverse(order[i]) ? '-' : '+') << endl;
#endif
        }
        size_t length_thru = dists[i] + graph.get_length(order[i]);
        bool is_sink = true;;
        graph.follow_edges(order[i], false, [&](const handle_t& next) {
            size_t j = index_of.at(next);
            dists[j] = min(dists[j], length_thru);
#ifdef debug_deletion_aligner
            cerr << "Edge from " << graph.get_id(order[i]) << (graph.get_is_reverse(order[i]) ? '-' : '+')
                << " to " << graph.get_id(next) << (graph.get_is_reverse(next) ? '-' : '+')
                << " assigns distance " << dists[j] << endl;
#endif
            is_sink = false;
        });
        if (is_sink) {
            // didn't find any edges forward
            sinks.emplace_back(i, length_thru);
#ifdef debug_deletion_aligner
            cerr << "Declare a sink at " << graph.get_id(order[i]) << (graph.get_is_reverse(order[i]) ? '-' : '+') << " with distance " << length_thru << endl;
#endif
        }
    }
#ifdef debug_deletion_aligner
    cerr << "min dist results:" << endl;
    for (size_t i = 0; i < order.size(); ++i) {
        cerr << "\t" << i << " (node " << graph.get_id(order[i]) << "): " << dists[i] << endl;
    }
    cerr << "sinks:" << endl;
    for (auto sink : sinks) {
        cerr << "\t" << sink.first << " " << sink.second << endl;
    }
#endif
    return make_pair(dists, sinks);
}

vector<vector<handle_t>> DeletionAligner::traceback(const vector<handle_t>& order,
                                                    const unordered_map<handle_t, size_t>& index_of,
                                                    const HandleGraph& graph,
                                                    const vector<size_t>& dists,
                                                    const vector<pair<size_t, size_t>>& sinks,
                                                    size_t max_tracebacks) const {

#ifdef debug_deletion_aligner
    cerr << "beginning multi-traceback" << endl;
#endif
    // records of (distance, deflections (from, to))
    structures::MinMaxHeap<pair<size_t, vector<pair<size_t, size_t>>>> heap;
    vector<vector<handle_t>> traces;
    
    // check if we want to take this deviation from the optimal traceback next time
    auto propose_deflection = [&](size_t from, size_t to, size_t dist,
                                  const vector<pair<size_t, size_t>>& curr_deflections) {
#ifdef debug_deletion_aligner
        cerr << "proposing deflection from " << from << " to " << to << " with dist " << dist << endl;
#endif
        if (heap.size() + traces.size() < max_tracebacks ||
            (!heap.empty() && heap.max().first > dist)) {
            // we've used all of the current deflections (so we can now propose more)
            // and we either haven't fully populated the heap, or this is better than
            // the worst
#ifdef debug_deletion_aligner
            cerr << "accepted deflection" << endl;
#endif
            vector<pair<size_t, size_t>> deflections = curr_deflections;
            deflections.emplace_back(from, to);
            heap.emplace(dist, std::move(deflections));
            if (heap.size() + traces.size() > max_tracebacks) {
#ifdef debug_deletion_aligner
                cerr << "ejecting traceback with dist " << heap.max().first << endl;
#endif
                heap.pop_max();
            }
        }
    };
    
    // get the next trace, either by taking a deflection or by doing traceback,
    // also propose deflections as needed
    auto get_next = [&](size_t at, size_t& deflxn, size_t curr_dist,
                        const vector<pair<size_t, size_t>>& curr_deflections) {
        if (deflxn < curr_deflections.size()
            && at == curr_deflections[deflxn].first) {
            return curr_deflections[deflxn++].second;
        }
        else {
            size_t next = numeric_limits<size_t>::max();
            size_t dist_here = dists[at];
            graph.follow_edges(order[at], true, [&](const handle_t& prev) {
                size_t idx = index_of.at(prev);
                size_t dist_thru = dists[idx] + graph.get_length(prev);
                if (next == numeric_limits<size_t>::max() && dist_thru == dist_here) {
                    next = idx;
                }
                else if (deflxn == curr_deflections.size()) {
                    propose_deflection(at, idx, curr_dist - dist_here + dist_thru,
                                       curr_deflections);
                }
                // keep looking if we haven't found the trace or we're looking
                // for deflections
                return (next == numeric_limits<size_t>::max() ||
                        deflxn == curr_deflections.size());
            });
            return next;
        }
    };
    
    // the first deflection is from a sentinel at the end to whichever
    // sink node we'll start from, propose these deflections to init the heap
    vector<pair<size_t, size_t>> deflections;
    for (auto& sink : sinks) {
        propose_deflection(order.size(), sink.first, sink.second,
                           deflections);
    }
    
    traces.reserve(max_tracebacks);
    while (!heap.empty()) {
        // init the next traceback
        traces.emplace_back();
        auto& trace = traces.back();
        size_t trace_dist;
        tie(trace_dist, deflections) = heap.min();
        heap.pop_min();
        
#ifdef debug_deletion_aligner
        cerr << "beginning trace of dist " << trace_dist << " with deflections:" << endl;
        for (auto d : deflections) {
            cerr << "\t" << d.first << " -> " << d.second << endl;
        }
#endif
        
        // start by taking deflection from the beginning sentinel
        size_t deflxn = 0;
        size_t tracer = get_next(order.size(), deflxn, trace_dist,
                                 deflections);
        while (tracer != numeric_limits<size_t>::max()) {
#ifdef debug_deletion_aligner
            cerr << "doing trace on " << tracer << endl;
#endif
            trace.push_back(order[tracer]);
            tracer = get_next(tracer, deflxn, trace_dist, deflections);
        }
    }
    return traces;
}

void DeletionAligner::trace_to_alignment(Alignment& aln, const vector<handle_t>& trace,
                                         const HandleGraph& graph) const {
    int64_t total_dist = 0;
    auto path = aln.mutable_path();
    // traces are constructed in reverse
    for (auto it = trace.rbegin(); it != trace.rend(); ++it) {
        handle_t handle = *it;
        auto mapping = path->add_mapping();
        auto position = mapping->mutable_position();
        position->set_node_id(graph.get_id(handle));
        position->set_is_reverse(graph.get_is_reverse(handle));
        auto edit = mapping->add_edit();
        edit->set_from_length(graph.get_length(handle));
        total_dist += edit->from_length();
    }
    // TODO: ideally scoring would live in the Aligner, but i guess it's okay
    aln.set_score(total_dist ? -gap_open - (total_dist - 1) * gap_extension : 0);
}

}
