#include <atomic>
#include <charconv>
#include <chrono>
#include <cstdio>
#include <limits>

#include <omp.h>

#include "graph_caller.hpp"
#include "symbolic_allele.hpp"
#include "read_likelihood_caller.hpp"
#include "algorithms/expand_context.hpp"
#include "annotation.hpp"
#include "gref.hpp"
#include "traversal_clusters.hpp"
#include "utility.hpp"

//#define debug

namespace vg {


GraphCaller::GraphCaller(SnarlCaller& snarl_caller,
                         SnarlManager& snarl_manager) :
    snarl_caller(snarl_caller), snarl_manager(snarl_manager), show_progress(false) {
}

GraphCaller::~GraphCaller() {
}

void GraphCaller::set_show_progress(bool show_progress) {
    this->show_progress = show_progress;
}

void GraphCaller::set_snarl_batching(size_t window_size) {
    snarl_batch_window = window_size;
}

/// Key a snarl by the lower of its two boundary node IDs, so that sorting on it puts snarls with
/// nearby node IDs next to each other.
static nid_t snarl_node_key(const Snarl* snarl) {
    return min(snarl->start().node_id(), snarl->end().node_id());
}

/// Sort snarls by `snarl_node_key`, so that snarls with nearby node IDs are called one after
/// another.
///
/// Each snarl's key is read once, beforehand: reading it at every comparison followed a pointer
/// into each snarl's scattered Snarl message, millions of times for a whole genome's top-level
/// snarls. The sort compares the same keys in the same sequence, so the result, including the
/// order of snarls with equal keys, is the one sorting the snarls themselves gives.
static void sort_snarls_by_node_id(vector<const Snarl*>& snarls) {
    vector<pair<nid_t, const Snarl*>> keyed(snarls.size());
#pragma omp parallel for schedule(static)
    for (size_t i = 0; i < snarls.size(); ++i) {
        keyed[i] = make_pair(snarl_node_key(snarls[i]), snarls[i]);
    }
    std::sort(keyed.begin(), keyed.end(),
              [](const pair<nid_t, const Snarl*>& a, const pair<nid_t, const Snarl*>& b) {
                  return a.first < b.first;
              });
    for (size_t i = 0; i < snarls.size(); ++i) {
        snarls[i] = keyed[i].second;
    }
}

void GraphCaller::call_top_level_snarls(const HandleGraph& graph, RecurseType recurse_type) {

    // Used to recurse on children of parents that can't be called
    size_t thread_count = get_thread_count();
    vector<vector<const Snarl*>> snarl_queue(thread_count);

    std::atomic<std::int64_t> top_snarl_count(0);
    std::atomic<std::int64_t> nested_snarl_count(0);
    bool top_level = true;

    // Run the snarl caller on a snarl, and queue up the children if it fails
    auto process_snarl = [&](const Snarl* snarl) {

        if (!snarl_manager.is_trivial(snarl, graph)) {

#ifdef debug
            cerr << "GraphCaller running call_snarl on " << pb2json(*snarl) << endl;
#endif

            bool was_called = call_snarl(*snarl);
            if (recurse_type == RecurseAlways || (!was_called && recurse_type == RecurseOnFail)) {
                const vector<const Snarl*>& children = snarl_manager.children_of(snarl);
                vector<const Snarl*>& thread_queue = snarl_queue[omp_get_thread_num()];
                thread_queue.insert(thread_queue.end(), children.begin(), children.end());
            }
            
            if (show_progress) {
                if (top_level) {
                    ++top_snarl_count;
                    if (top_snarl_count % 100000 == 0) {
#pragma omp critical (cerr)
                        cerr << "[vg call]: Processed " << top_snarl_count << " top-level snarls" << endl;
                    }
                } else {
                    ++nested_snarl_count;
                    if (nested_snarl_count % 100000 == 0) {
#pragma omp critical (cerr)                    
                        cerr << "[vg call]: Processed " << top_snarl_count << " nested snarls" << endl;
                    }
                }
            }
        }
    };

    // Start with the top level snarls
    if (snarl_batch_window > 0) {
        // Call in node-ID order, batched by window, so that snarls with nearby node IDs are
        // called together.
        vector<const Snarl*> roots;
        snarl_manager.for_each_top_level_snarl([&](const Snarl* snarl) {
            roots.push_back(snarl);
        });
        sort_snarls_by_node_id(roots);

        // Split the sorted snarls into runs that share a window, with one parallel job per run
        // rather than per snarl.
        vector<pair<size_t, size_t>> windows;
        size_t begin = 0;
        while (begin < roots.size()) {
            size_t window = (size_t)(snarl_node_key(roots[begin]) / (nid_t)snarl_batch_window);
            size_t end = begin + 1;
            while (end < roots.size() &&
                   (size_t)(snarl_node_key(roots[end]) / (nid_t)snarl_batch_window) == window) {
                ++end;
            }
            windows.emplace_back(begin, end);
            begin = end;
        }

#pragma omp parallel for schedule(dynamic, 1)
        for (int w = 0; w < (int)windows.size(); ++w) {
            for (size_t i = windows[w].first; i < windows[w].second; ++i) {
                process_snarl(roots[i]);
            }
        }
    } else {
        snarl_manager.for_each_top_level_snarl_parallel(process_snarl);
    }
    if (show_progress) cerr << "[vg call]: Finished processing " << top_snarl_count << " top-level snarls" << endl;

    top_level = false;

    // Then recurse on any children the snarl caller failed to handle
    while (!std::all_of(snarl_queue.begin(), snarl_queue.end(),
                        [](const vector<const Snarl*>& snarl_vec) {return snarl_vec.empty();})) {
        vector<const Snarl*> cur_queue;
        for (vector<const Snarl*>& thread_queue : snarl_queue) {
            cur_queue.reserve(cur_queue.size() + thread_queue.size());
            std::move(thread_queue.begin(), thread_queue.end(), std::back_inserter(cur_queue));
            thread_queue.clear();
        }

        if (snarl_batch_window > 0) {
            // Keep queued children in node-ID order too.
            sort_snarls_by_node_id(cur_queue);
        }

#pragma omp parallel for schedule(dynamic, 1)
        for (int i = 0; i < cur_queue.size(); ++i) {
            process_snarl(cur_queue[i]);
        }
    }
    if (show_progress && nested_snarl_count > 0) cerr << "[vg call]: Finished processing " << nested_snarl_count << " nested snarls" << endl;
  
}

void GraphCaller::call_top_level_chains(const HandleGraph& graph, size_t max_edges, size_t max_trivial, RecurseType recurse_type) {
    // Used to recurse on children of parents that can't be called
    size_t thread_count = get_thread_count();
    vector<vector<Chain>> chain_queue(thread_count);

    // Run the snarl caller on a chain. queue up the children if it fails
    auto process_chain = [&](const Chain* chain) {

#ifdef debug
        cerr << "calling top level chain ";
        for (const auto& i : *chain) {
            cerr << pb2json(*i.first) << "," << i.second << ",";
        }
        cerr << endl;
#endif
        // Break up the chain
        vector<Chain> chain_pieces = break_chain(graph, *chain, max_edges, max_trivial);

        for (Chain& chain_piece : chain_pieces) {
            // Make a fake snarl spanning the chain
            // It is important to remember that along with not actually being a snarl,
            // it's not managed by the snarl manager so functions looking into its nesting
            // structure will not work
            Snarl fake_snarl;
            *fake_snarl.mutable_start() = chain_piece.front().second == true ? reverse(chain_piece.front().first->end()) :
                chain_piece.front().first->start();
            *fake_snarl.mutable_end() = chain_piece.back().second == true ? reverse(chain_piece.back().first->start()) :
                chain_piece.back().first->end();

#ifdef debug
            cerr << "calling fake snarl " << pb2json(fake_snarl) << endl;
#endif
            
            bool was_called = call_snarl(fake_snarl);
            if (recurse_type == RecurseAlways || (!was_called && recurse_type == RecurseOnFail)) {
                vector<Chain>& thread_queue = chain_queue[omp_get_thread_num()];                
                for (pair<const Snarl*, bool> chain_link : chain_piece) {
                    const deque<Chain>& child_chains = snarl_manager.chains_of(chain_link.first);
                    thread_queue.insert(thread_queue.end(), child_chains.begin(), child_chains.end());
                }
            }
        }
    };

    // Start with the top level snarls
    snarl_manager.for_each_top_level_chain_parallel(process_chain);

    // Then recurse on any children the snarl caller failed to handle
    while (!std::all_of(chain_queue.begin(), chain_queue.end(),
                        [](const vector<Chain>& chain_vec) {return chain_vec.empty();})) {
        vector<Chain> cur_queue;
        for (vector<Chain>& thread_queue : chain_queue) {
            cur_queue.reserve(cur_queue.size() + thread_queue.size());
            std::move(thread_queue.begin(), thread_queue.end(), std::back_inserter(cur_queue));
            thread_queue.clear();
        }

#pragma omp parallel for schedule(dynamic, 1)
        for (int i = 0; i < cur_queue.size(); ++i) {
            process_chain(&cur_queue[i]);
        }
    
    }
}

vector<Chain> GraphCaller::break_chain(const HandleGraph& graph, const Chain& chain, size_t max_edges, size_t max_trivial) {
    
    vector<Chain> chain_frags;

    // keep track of the current fragment and add it to chain_frags as soon as it gets too big
    Chain frag;
    size_t frag_edge_count = 0;
    size_t frag_triv_count = 0;
    
    for (const pair<const Snarl*, bool>& link : chain) {
        // todo: we're getting the contents here as well as within the caller.
        auto contents = snarl_manager.deep_contents(link.first, graph, false);

        // todo: use annotation from snarl itself?
        bool trivial = contents.second.empty();

        if ((trivial && frag_triv_count > max_trivial) ||
            (contents.second.size() + frag_edge_count > max_edges)) {
            // adding anything more to the chain would make it too long, so we
            // add it to the output and clear the current fragment
            if (!frag.empty() && frag_triv_count < frag.size()) {
                chain_frags.push_back(frag);
            }
            frag.clear();
            frag_edge_count = 0;
            frag_triv_count = 0;
        }

        if (!trivial || (frag_triv_count < max_trivial)) {
            // we start a new fragment or add to an existing fragment
            frag.push_back(link);
            frag_edge_count += contents.second.size();
            if (trivial) {
                ++frag_triv_count;
            }
        }
    }

    // and the last one
    if (!frag.empty()) {
        chain_frags.push_back(frag);
    }

    return chain_frags;
}
}

