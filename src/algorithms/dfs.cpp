#include "dfs.hpp"

namespace vg {
namespace algorithms {

using namespace std;

// depth first search across node traversals with interface to traversal tree via callback
void dfs(
    const HandleGraph& graph,
    const function<void(const handle_t&)>& handle_begin_fn,  // called when node orientation is first encountered
    const function<void(const handle_t&)>& handle_end_fn,    // called when node orientation goes out of scope
    const function<bool(void)>& break_fn,                    // called to check if we should stop the DFS
    const function<void(const edge_t&)>& edge_fn,            // called when an edge is encountered
    const function<void(const edge_t&)>& tree_fn,            // called when an edge forms part of the DFS spanning tree
    const function<void(const edge_t&)>& edge_curr_fn,       // called when we meet an edge in the current tree component
    const function<void(const edge_t&)>& edge_cross_fn,      // called when we meet an edge in an already-traversed tree component
    const vector<handle_t>& sources,                         // start only at these node traversals
    const unordered_set<handle_t>& sinks                     // when hitting a sink, don't keep walking
    ) {

    // to maintain search state
    enum SearchState { PRE = 0, CURR, POST };
    unordered_map<handle_t, SearchState> state; // implicitly constructed entries will be PRE.

    // to maintain stack frames
    struct Frame {
        handle_t handle;
        vector<edge_t>::iterator begin, end;
        Frame(handle_t h,
              vector<edge_t>::iterator b,
              vector<edge_t>::iterator e)
            : handle(h), begin(b), end(e) { }
    };

    // maintains edges while the node traversal's frame is on the stack
    unordered_map<handle_t, vector<edge_t> > edges;

    // do dfs from given root.  returns true if terminated via break condition, false otherwise
    function<bool(const handle_t&)> dfs_single_source = [&](const handle_t& root) {
                
        // to store the stack frames
        deque<Frame> todo;
        if (state[root] == SearchState::PRE) {
            state[root] = SearchState::CURR;
                
            // Collect all the edges attached to the outgoing side of the
            // traversal.
            auto& es = edges[root];
            // follow edges?
            graph.follow_edges(root, false, [&](const handle_t& next) {
                    es.push_back(graph.edge_handle(root, next));
                });
                
            todo.push_back(Frame(root, es.begin(), es.end()));
            // run our discovery-time callback
            handle_begin_fn(root);
            // and check if we should break
            if (break_fn()) {
                return true;
            }
        }
        // now begin the search rooted at this NodeTraversal
        while (!todo.empty()) {
            // get the frame
            auto& frame = todo.back();
            // and set up reference to it
            auto handle = frame.handle;
            auto edges_begin = frame.begin;
            auto edges_end = frame.end;
            todo.pop_back();
            // run through the edges to handle
            while (edges_begin != edges_end) {
                auto& edge = *edges_begin;
                // run the edge callback
                edge_fn(edge);
                    
                // what's the handle we'd get to following this edge
                const handle_t& target = edge.second;

                auto search_state = state[target];
                // if we've not seen it, follow it
                if (search_state == SearchState::PRE) {
                    tree_fn(edge);
                    // save the rest of the search for this handle on the stack
                    todo.push_back(Frame(handle, ++edges_begin, edges_end));
                    // switch our focus to the handle at the other end of the edge
                    handle = target;
                    // and store it on the stack
                    state[handle] = SearchState::CURR;
                    auto& es = edges[handle];

                    // only walk out of handle that are not the sink
                    if (sinks.empty() || sinks.count(handle) == false) {
                        // follow edges?
                        graph.follow_edges(handle, false, [&](const handle_t& next) {
                                es.push_back(graph.edge_handle(handle, next));
                            });
                    }
                    edges_begin = es.begin();
                    edges_end = es.end();
                    // run our discovery-time callback
                    handle_begin_fn(handle);
                } else if (search_state == SearchState::CURR) {
                    // if it's on the stack
                    edge_curr_fn(edge);
                    ++edges_begin;
                } else {
                    // it's already been handled, so in another part of the tree
                    edge_cross_fn(edge);
                    ++edges_begin;
                }
            }
            state[handle] = SearchState::POST;
            handle_end_fn(handle);
            edges.erase(handle); // clean up edge cache
        }

        return false;
    };

    if (sources.empty()) {
        // attempt the search rooted at all NodeTraversals
        graph.for_each_handle([&](const handle_t& handle_fwd) {
                handle_t handle_rev = graph.flip(handle_fwd);
                dfs_single_source(handle_fwd);
                dfs_single_source(handle_rev);
            });
    } else {
        for (auto source : sources) {
            dfs_single_source(source);
        }
    }
}

void dfs(const HandleGraph& graph,
                      const function<void(const handle_t&)>& handle_begin_fn,
                      const function<void(const handle_t&)>& handle_end_fn,
                      const vector<handle_t>& sources,
                      const unordered_set<handle_t>& sinks) {
    auto edge_noop = [](const edge_t& e) { };
    dfs(graph,
        handle_begin_fn,
        handle_end_fn,
        [](void) { return false; },
        edge_noop,
        edge_noop,
        edge_noop,
        edge_noop,
        sources,
        sinks);
}

void dfs(const HandleGraph& graph,
                      const function<void(const handle_t&)>& handle_begin_fn,
                      const function<void(const handle_t&)>& handle_end_fn,
                      const function<bool(void)>& break_fn) {
    auto edge_noop = [](const edge_t& e) { };
    vector<handle_t> empty_sources;
    unordered_set<handle_t> empty_sinks;
    dfs(graph,
        handle_begin_fn,
        handle_end_fn,
        break_fn,
        edge_noop,
        edge_noop,
        edge_noop,
        edge_noop,
        empty_sources,
        empty_sinks);
}

}
}
