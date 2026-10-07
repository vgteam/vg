#include <atomic>
#include <charconv>
#include <chrono>
#include <cstdio>
#include <limits>

#include <omp.h>

#include "graph_caller.hpp"
#include "site_scheduler.hpp"
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

void GraphCaller::call_top_level_snarls(const HandleGraph& graph, RecurseType recurse_type) {

    std::atomic<std::int64_t> top_snarl_count(0);
    std::atomic<std::int64_t> nested_snarl_count(0);
    bool top_level = true;

    // Run the snarl caller on a snarl, and queue up the children if it fails
    auto process_snarl = [&](const Snarl* const& snarl, vector<const Snarl*>& thread_queue) {

        if (!snarl_manager.is_trivial(snarl, graph)) {

#ifdef debug
            cerr << "GraphCaller running call_snarl on " << pb2json(*snarl) << endl;
#endif

            bool was_called = call_snarl(*snarl);
            if (recurse_type == RecurseAlways || (!was_called && recurse_type == RecurseOnFail)) {
                const vector<const Snarl*>& children = snarl_manager.children_of(snarl);
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
    SiteScheduler<const Snarl*> scheduler(snarl_batch_window, snarl_node_key, process_snarl);

    // Start with the top level snarls
    scheduler.call_top_level(
        [&](const function<void(const Snarl* const&)>& lambda) {
            snarl_manager.for_each_top_level_snarl([&](const Snarl* snarl) { lambda(snarl); });
        },
        [&](const function<void(const Snarl* const&)>& lambda) {
            snarl_manager.for_each_top_level_snarl_parallel([&](const Snarl* snarl) { lambda(snarl); });
        });
    if (show_progress) cerr << "[vg call]: Finished processing " << top_snarl_count << " top-level snarls" << endl;

    top_level = false;

    // Then recurse on any children the snarl caller failed to handle
    scheduler.call_queued();
    if (show_progress && nested_snarl_count > 0) cerr << "[vg call]: Finished processing " << nested_snarl_count << " nested snarls" << endl;
  
}

void GraphCaller::call_top_level_chains(const HandleGraph& graph, size_t max_edges, size_t max_trivial, RecurseType recurse_type) {
    // Run the snarl caller on a chain. queue up the children if it fails
    auto process_chain = [&](const Chain& chain, vector<Chain>& thread_queue) {

#ifdef debug
        cerr << "calling top level chain ";
        for (const auto& i : chain) {
            cerr << pb2json(*i.first) << "," << i.second << ",";
        }
        cerr << endl;
#endif
        // Break up the chain
        vector<Chain> chain_pieces = break_chain(graph, chain, max_edges, max_trivial);

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
                for (pair<const Snarl*, bool> chain_link : chain_piece) {
                    const deque<Chain>& child_chains = snarl_manager.chains_of(chain_link.first);
                    thread_queue.insert(thread_queue.end(), child_chains.begin(), child_chains.end());
                }
            }
        }
    };

    // Chains are never batched.
    SiteScheduler<Chain> scheduler(0, nullptr, process_chain);

    // Start with the top level snarls
    scheduler.call_top_level(nullptr, [&](const function<void(const Chain&)>& lambda) {
        snarl_manager.for_each_top_level_chain_parallel([&](const Chain* chain) { lambda(*chain); });
    });

    // Then recurse on any children the snarl caller failed to handle
    scheduler.call_queued();
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

