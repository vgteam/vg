///
///  \file integrated_snarl_finder.cpp
///
///

#include "integrated_snarl_finder.hpp"

#include "algorithms/three_edge_connected_components.hpp"

#include <bdsg/overlays/overlay_helper.hpp>
#include <structures/union_find.hpp>

#include <utility>

namespace vg {

using namespace std;

class IntegratedSnarlFinder::MergedAdjacencyGraph {
protected:
    /// Hold onto the backing HandleGraph
    const HandleGraph* graph;
    
    /// Keep a vectorizable overlay over it to let us map between handles
    /// and union-find indices via handle ranking. The handles are all at index
    /// (rank - 1) * 2 + is_reverse.
    ///
    /// We rely on handles in the vectorizable overlay and handles in the
    /// backing graph being identical.
    ///
    /// TODO: make not-mutable when https://github.com/vgteam/libbdsg/issues/63
    /// is fixed.
    mutable bdsg::VectorizableOverlayHelper overlay_helper;
    
    /// Keep a union-find over the ranks of the merged oriented handles that
    /// make up each component. Runs with include_children=true so we can find
    /// all the members of each group.
    ///
    /// Needs to be mutable because union-find find operations do internal tree
    /// massaging and aren't const.
    /// TODO: this makes read operations not thread safe!
    mutable structures::UnionFind union_find;
    
    /// Get the rank corresponding to the given handle, in the union-find.
    /// Our ranks are 0-based.
    size_t uf_rank(handle_t into) const;
    
    /// Get the handle with the given rank in union-find space.
    /// Our ranks are 0-based.
    handle_t uf_handle(size_t rank) const;
    
public:
    /// Make a MergedAdjacencyGraph representing the graph of adjacency components of the given HandleGraph.
    MergedAdjacencyGraph(const HandleGraph* graph);
    
    /// Copy a MergedAdjacencyGraph by re-doing all the merges. Uses its own internal vectorization.
    MergedAdjacencyGraph(const MergedAdjacencyGraph& other);
    
    /// Given handles reading into two components, a and b, merge them into a single component.
    void merge(handle_t into_a, handle_t into_b);
    
    /// Find the handle heading the component that the given handle is in.
    handle_t find(handle_t into) const; 
    
    /// For each head, call the iteratee.
    void for_each_head(const function<void(handle_t)>& iteratee) const;
    
    /// For each item other than the head in the component headed by the given
    /// handle, calls the iteratee with the that other item. Does not call the
    /// iteratee for single-item components.
    void for_each_member(handle_t head, const function<void(handle_t)>& iteratee) const;
    
    /// For each item other than the head in each component, calls the iteratee
    /// with the head and the other item. Does not call the iteratee for
    /// single-item components.
    void for_each_membership(const function<void(handle_t, handle_t)>& iteratee) const;
    
    /// Return a vector of, for each connected component, the length in bases
    /// and the path edges of the longest simple cycle in it, if any.
    vector<pair<size_t, vector<handle_t>>> longest_cycles_in_connected_components() const;
    
    /// Return the path length (total edge length in bp) and edges for the longst path in each tree in a forest.
    vector<pair<size_t, vector<handle_t>>> longest_paths_in_forest() const;
};

size_t IntegratedSnarlFinder::MergedAdjacencyGraph::uf_rank(handle_t into) const {
    // We need to 0-base the backing rank, space it out, and make the low bit orientation
    return (overlay_helper.get()->id_to_rank(graph->get_id(into)) - 1) * 2 + (size_t) graph->get_is_reverse(into);
}

handle_t IntegratedSnarlFinder::MergedAdjacencyGraph::uf_handle(size_t rank) const {
    // We nefor here.ed to take the high bits and than make it 1-based, and get the orientation from the low bit
    return graph->get_handle(overlay_helper.get()->rank_to_id(rank / 2 + 1), rank % 2);
}

IntegratedSnarlFinder::MergedAdjacencyGraph::MergedAdjacencyGraph(const HandleGraph* graph) : graph(graph),
    overlay_helper(), union_find(graph->get_node_count() * 2, true) {
    
    // Make sure we have our vectorizable version of the graph.
    // TODO: remove const_cast when https://github.com/vgteam/libbdsg/issues/64 is fixed.
    overlay_helper.apply(const_cast<HandleGraph*>(graph));
    
    // TODO: we want the adjacency components that are just single edges
    // between two handles (i.e. trivial snarls) to be implicit, so we don't
    // have to do O(n) work for so much of the graph. But to do that we need a
    // union-find that lets us declare it over a potentially large space
    // without filling it all in.
    
    // So we do this the easy way and compute all the merges for all adjacency
    // components, including tiny/numerous ones, right now.
    
    // If we ever change this, we should also make MergedAdjacencyGraph
    // stackable, to save a copy when we want to do further merges but keep the
    // old state.
    
    graph->for_each_edge([&](const handlegraph::edge_t& e) {
        // Get the inward-facing version of the second handle
        auto into_b = graph->flip(e.second);
        
        // Merge to create initial adjacency components
        merge(e.first, into_b);
    });
}

IntegratedSnarlFinder::MergedAdjacencyGraph::MergedAdjacencyGraph(const MergedAdjacencyGraph& other) : MergedAdjacencyGraph(other.graph) {
    other.for_each_membership([&](handle_t head, handle_t member) {
        // For anything in a component, other than its head, do the merge with the head.
        merge(head, member);
    });
}

void IntegratedSnarlFinder::MergedAdjacencyGraph::merge(handle_t into_a, handle_t into_b) {
    // Get ranks and merge
    union_find.union_groups(uf_rank(into_a), uf_rank(into_b));
}

handle_t IntegratedSnarlFinder::MergedAdjacencyGraph::find(handle_t into) const {
    // Get rank, find head, and get handle
    return uf_handle(union_find.find_group(uf_rank(into)));
}

void IntegratedSnarlFinder::MergedAdjacencyGraph::for_each_head(const function<void(handle_t)>& iteratee) const {
    // TODO: is this better or worse than getting the vector of vectors of the whole union-find?
    // TODO: make iterating groups an actual capability that the union-find has, in O(groups).
    // This lets us do it in O(total items).
    
    // We track if we have seen a head yet. If we haven't, we emit it and mark it seen.
    vector<bool> seen_heads(union_find.size(), false);
    
    for (size_t i = 0; i < union_find.size(); i++) {
        // For each item in the union-find
        if (!seen_heads[i]) {
            // If we haven't emitted it, find the head of its group
            size_t head = union_find.find_group(i);
            if (!seen_heads[head]) {
                // If we haven't emitted that head either, say we have
                seen_heads[head] = true;
                // And emit its corresponding inward-facing handle
                iteratee(uf_handle(head));
            }
        }
    }
}

void IntegratedSnarlFinder::MergedAdjacencyGraph::for_each_member(handle_t head, const function<void(handle_t)>& iteratee) const {
    size_t head_rank = uf_rank(head);
    // Find the group the head is in
    vector<size_t> group = union_find.group(head_rank);
    for (auto& member_rank : group) {
        // And go through all the members
        if (member_rank != head_rank) {
            // We filter out the given head.
            // This function will happen to work for non-head inputs, leaving out that input, but we don't guarantee it!
            iteratee(uf_handle(member_rank));
        }
    }
}

void IntegratedSnarlFinder::MergedAdjacencyGraph::for_each_membership(const function<void(handle_t, handle_t)>& iteratee) const {
    // We do this weird iteration because it's vaguely efficient in the union-find we use.
    vector<vector<size_t>> uf_components = union_find.all_groups();
    
    for (auto& component : uf_components) {
        // For each component
        for (size_t i = 1; i < component.size(); i++) {
            // For everything other than the head, announce with the head.
            iteratee(uf_handle(component[0]), uf_handle(component[i]));
        }
    }
}

vector<pair<size_t, vector<handle_t>>> IntegratedSnarlFinder::MergedAdjacencyGraph::longest_cycles_in_connected_components() const {
    // Do a DFS over all connected components of the graph
    
    // We will fill this in
    vector<pair<size_t, vector<handle_t>>> to_return;
    
    // To let us measure cycles without walking them, we need to remember stack
    // frame index as our visit marker. No record = not visited, record =
    // visited.
    //
    // If you see something visited already, it must still be on the stack.
    // Otherwise it would have already visited you when it was on the stack.
    //
    // This is on heads, representing nodes.
    unordered_map<handle_t, size_t> visited_frame;
    
    // We need a stack.
    // Stack is actually in terms of inward edges followed.
    struct DFSFrame {
        handle_t here;
        vector<handle_t> todo;
        // Track the path length in bp from the root so we can easily measure cycles via difference.
        size_t path_length_bp = 0;
    };
    
    vector<DFSFrame> stack;
    
    for_each_head([&](handle_t component_root) {
        // For every node in the graph
        
        if (!visited_frame.count(component_root)) {
            // If it hasn't been searched yet, start a search of its connected component.
            stack.emplace_back();
            stack.back().here = component_root;
            
            // We'll put the result here, or get rid of it if we find no cycle.
            to_return.emplace_back();
            auto& longest_cycle = to_return.back();
            
            while (!stack.empty()) {
                // Until the DFS is done
                auto& frame = stack.back();
                // Find the node that following this edge got us to.
                auto frame_head = find(frame.here);
                auto frame_it = visited_frame.find(frame_head);
                if (frame_it == visited_frame.end()) {
                    // First visit to here.
                    
                    // Mark visited at this stack level
                    frame_it = visited_frame.emplace_hint(frame_it, frame_head, stack.size() - 1);
                    
                    if (stack.size() > 1) {
                        // Record path length, including the edge we took
                        frame.path_length_bp = stack[stack.size() - 2].path_length_bp + graph->get_length(frame.here);
                    }
                    
                    // Queue up edges
                    for_each_member(frame_head, [&](handle_t member) {
                        // Follow edge by flipping. But queue up the edge
                        // followed instead of the node reached (head), so we
                        // can emit the cycle later in terms of edges.
                        frame.todo.push_back(graph->flip(member));
                    });
                }
                
                if (!frame.todo.empty()) {
                    // Now do an edge
                    handle_t edge_into = frame.todo.back();
                    handle_t connected_head = find(edge_into);
                    frame.todo.pop_back();
                    
                    auto connected_it = visited_frame.find(connected_head);
                    
                    if (connected_it == visited_frame.end()) {
                        // Forward edge. Recurse.
                        // TODO: this immediately does a lookup in the hash table again.
                        stack.emplace_back();
                        stack.back().here = edge_into;
                    } else {
                        // Back edge
                        if (frame_it->second >= connected_it->second) {
                            // We have an edge to something that was visited at
                            // our stack level or earlier; it is still on the
                            // stack.
                        
                            // We only care about the back edges looking up, or
                            // being self loops. Only do them one way to avoid
                            // emitting cycles twice.
                            
                            // Get the stack frame we connected to (which may
                            // be the same as our frame, for self loops)
                            auto& connected_frame = stack[connected_it->second];
                            
                            // Measure the cycle: length down the path in the
                            // tree, plus length of the back edge.
                            size_t cycle_length = frame.path_length_bp - connected_frame.path_length_bp + graph->get_length(edge_into);
                            
                            if (cycle_length > longest_cycle.first || longest_cycle.second.empty()) {
                                // This is a new longest cycle.
                                // Record its length in bases.
                                longest_cycle.first = cycle_length;
                                // Allocate space for all its edges: the ones
                                // on the tree path between connected and here,
                                // plus one to close.
                                longest_cycle.second.resize(frame_it->second - connected_it->second + 1);
                                for (size_t i = 0; i + 1 < longest_cycle.second.size(); i++) {
                                    // Copy all the edges after connected into
                                    // all but the last slot in the cycle
                                    longest_cycle[i] = stack[connected_it->second + 1 + i].here;
                                }
                                // And at the end put the edge to close the
                                // cycle
                                longest_cycle.back() = edge_into;
                            }
                        }
                    }
                } else {
                    // Clean up
                    stack.pop_back();
                }
            }
            
            if (to_return.back().second.empty()) {
                // There's no cycle in this connected component.
                to_return.pop_back();
            }
        }
    });
    
    // Hand back all the longest cycles in connected components
    return to_return;
}

vector<pair<size_t, vector<handle_t>>> IntegratedSnarlFinder::MergedAdjacencyGraph::longest_paths_in_forest() const {
    
    // TODO: somehow unify DFS logic with cycle-finding DFS in a way that still allows us to inspect our stack in each case?
    
    // Going up the tree, we need to track the longest path from a leaf to the subtree root, and the longest path between leaves in the subtree.
    // These ought to overlap substantially, but either may be the real winner when we get to where we rooted the DFS.
    
    // When we find a longest path, we put its length and value in here.
    // We describe it as edges followed.
    vector<pair<size_t, vector<handle_t>>> to_return;
    
    // We only need the tree edges, so we only need visited flags on the heads representing nodes.
    unordered_set<handle_t> visited;
    
    // We need a stack.
    // Stack is actually in terms of inward edges followed.
    struct DFSFrame {
        handle_t here;
        // What edges still need to be followed
        vector<handle_t> todo;
        // When children finish they put entries here with path length and leaf-to-root edge path
        vector<pair<size_t, vector<handle_t>> child_results;
    };
    
    vector<DFSFrame> stack;
    
    for_each_head([&](handle_t head) {
        // For every node in the graph
        
        if (!visited.count(head)) {
            // If it hasn't been searched yet, start a search
            stack.emplace_back();
            stack.back().here = head;
            
            while (!stack.empty()) {
                // Until the DFS is done
                auto& frame = stack.back();
                // Find the node that following this edge got us to.
                auto frame_head = find(frame.here);
                
                if (!visited.count(frame_head)) {
                    // First visit to here.
                    
                    // Mark visited
                    visited.insert(frame_head);
                    
                    // Queue up edges
                    for_each_member(frame_head, [&](handle_t member) {
                        // Follow edge by flipping. But queue up the edge
                        // followed instead of the node reached (head), so we
                        // can emit the cycle later in terms of edges.
                        frame.todo.push_back(graph->flip(member));
                    });
                    
                    // Allocate locations for child returns
                    frame.child_results.resrve(frame.todo.size());
                }
                
                if (!frame.todo.empty()) {
                    // Now do an edge
                    handle_t edge_into = frame.todo.back();
                    handle_t connected_head = find(edge_into);
                    frame.todo.pop_back();
                    
                    if (!visited.count(connected_head)) {
                        // Forward edge. Recurse.
                        stack.emplace_back();
                        stack.back().here = edge_into;
                    }
                } else {
                    // No children left.
                    
                    if (stack.size() > 1) {
                        // We aren't the root.
                        
                        // Find our max subtree length result
                        size_t best_child = numeric_limits<size_t>::max();
                        for (size_t i = 0; i < frame.child_results.size(); i++) {
                            // Look at results from all children
                            if (best_child == numeric_limits<size_t>::max() || frame.child_results[best_child].first < frame.child_results[i].first) {
                                // Select the one with the biggest total length
                                best_child = i;
                            }
                        }
                    
                        // Record the distance along us into our parent.
                        auto& parent_frame = stack[stack.size() - 2];
                        parent_frame.emplace_back(0, vector<handle_t>());
                        auto& our_result = parent_frame.back();
                        
                        if (best_child != numeric_limits<size_t>::max()) {
                            // We have a child result to forward
                            our_result = std::move(frame.child_results[best_child]);
                        }
                        
                        // Count our length
                        our_result.first += graph->get_length(frame.here);
                        
                        // And add us to our path, going up
                        our_result.second.push_back(graph->flip(frame.here));
                    } else {
                        // When we finally get to the root, we do it a bit differently.
                        
                        // TODO: implement!
                        
                        // Find the best and second best 
                        
                        // Pair them up and combine them
                        
                        // Save it to return (or just yield it?)
                    }
                
                    // Clean up
                    stack.pop_back();
                }
            }
                
                
            }
            
        }
        
    });
    
    
}

    




////////////////////////////////////////////////////////////////////////////////////////////





IntegratedSnarlFinder::IntegratedSnarlFinder(const PathHandleGraph& graph) : graph(&graph) {
    // Nothing to do!
}

void IntegratedSnarlFinder::for_each_snarl_including_trivial(const function<void(handle_t, handle_t)>& iteratee) const {
    // Do the actual snarl finding work and then feed the iteratee our snarls.
    
    // We need a union-find over the adjacency components of the graph, in which we will build the cactus graph.
    MergedAdjacencyGraph cactus(graph);
    
    // It magically gets the adjacency components itself.
    
    // Now we need to do the 3 edge connected component merging, using Tsin's algorithm.
    // We don't really have a good dense rank space on the adjacency components, so we use the general version.
    // TODO: Somehow have a nice dense rank space?
    // We represent each adjacency component (node) by its heading handle.
    algorithms::three_edge_connected_component_merges<handle_t>([&](const function<void(handle_t)>& emit_node) {
        // Feed all the handles that head adjacency components into the algorithm
        cactus.for_each_head(emit_node);
    }, [&](handle_t node, const function<void(handle_t)>& emit_edge) {
        // When asked for edges
        
        // Track what we have announced
        unordered_set<handle_t> seen_connected_heads;
        
        // We know the handle we got is a head.
        
        // Try and follow it as an edge and get the headed component
        handle_t connected_head = cactus.find(graph->flip(node));
        
        // Announce that. No adjacency component ever has no edges.
        seen_connected_heads.insert(connected_head);
        emit_edge(connected_head);
        
        cactus.for_each_member(node, [&](handle_t other_member) {
            // For each other handle in the adjacency component that this handle is heading
            
            // Follow as an edge again, by flipping
            handle_t member_connected_head = cactus.find(graph->flip(node));
            
            // See if we announced the component we reached already
            auto found = seen_connected_heads.find(member_connected_head);
            if (found == seen_connected_heads.end()) {
                // If we haven't seen it, insert it as seen with a hint
                seen_connected_heads.insert(found, member_connected_head);
                // And announce it
                emit_edge(member_connected_head);
            }
        });
    }, [&](handle_t a, handle_t b) {
        // Now we got a merge to create the 3 edge connected components.
        // Tell the graph.
        cactus.merge(a, b);
    });
    
    // Now our 3-edge-connected components have been condensed, and we have a proper Cactus graph.
    
    // Then we need to copy the base Cactus graph so we can make the bridge forest
    MergedAdjacencyGraph forest(cactus);
    
    // We remember the longest cycle in each component, and its length
    vector<pair<size_t, vector<handle_t>>> longest_cycles;
    
    // Now we find the simple cycles in the cactus graph
    cactus.for_each_simple_cycle([&](const function<void(const function<void(handle_t)>&)>& for_each_step) {
        // When we get a simple cycle (in terms of handles representing edges followed)
        
        // Compute its length in terms of fixed sequence on the edges followed.
        size_t cycle_length = 0;
        
        // And make a vector of its edges in case it is the longest.
        vector<handle_t> cycle_edges;
        
        // We need to remember the previous handle as we go around the cycle.
        bool is_first = false;
        handle_t last_handle;
        
        for_each_step([&](handle_t edge) {
            if (!cycle_edges.empty()) {
                // Merge with where the previous edge went in the bridge forest
                forest.merge(cycle_edges.back(), edge);
            }
        
            // Record the edge length into the cycle length
            cycle_length += graph->get_length(edge);
            
            // Record the edge
            cycle_edges.push_back(edge);
        });
        
        if (cycle_length > longest_cycle_length || longest_cycle_edges.empty()) {
            // Adopt this as the new longest cycle
            longest_cycle_length = cycle_length;
            longest_cycle_edges = std::move(cycle_edges);
        }
    });
    
    // Now we find the longest path in each tree in the bridge forest, with its length in bases.
    vector<pair<size_t, vector<handle_t>> longest_paths;
  
    {
        // Pick an unvisited node to arbitrarily root
        
        //
    }
  
   
    // TODO: Then find the lengths of the bridge edge spanning trees. We want
    // to be able to know the lengths no matter where we root, so I guess for
    // each node we need to store the longest length in each direction.
    
    // TODO: For each bridge forest tree, find the longest bridge edge path and
    // the longest simple cycle, and root with whichever is longer.
    
    // TODO: Recursively root things down the bridge edge trees off what we've
    // already rooted.
    
    // TODO: When we encounter a snarl, feed it to the iteratee. 
}

SnarlManager IntegratedSnarlFinder::find_snarls() {
    // Start with an empty SnarlManager
    SnarlManager snarl_manager;
    
    for_each_snarl_including_trivial([&](handle_t start_inward, handle_t end_outward) {
        // For every snarl, including the trivial ones
        
        // Make a Protobuf version of it
        vg::Snarl proto_snarl;
        
        // Convert boundary handles to Visits
        proto_snarl.mutable_start()->set_node_id(graph->get_id(start_inward));
        proto_snarl.mutable_start()->set_backward(graph->get_is_reverse(start_inward));
        proto_snarl.mutable_end()->set_node_id(graph->get_id(end_outward));
        proto_snarl.mutable_end()->set_backward(graph->get_is_reverse(end_outward));
        
        // TODO: Determine snarl metadata somehow. Do a postorder traversal? Or make the manager somehow do it?
    
        // Add the Protobuf version of the snarl to the SNarl Manager
        snarl_manager.add_snarl(proto_snarl);
    });
    
    // Let the snarl manager compute all its indexes
    snarl_manager.finish();
    
    // Give it back
    return snarl_manager;
}

}
