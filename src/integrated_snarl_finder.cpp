///
///  \file integrated_snarl_finder.cpp
///
///

#include "integrated_snarl_finder.hpp"

#include "algorithms/three_edge_connected_components.hpp"
#include "subgraph_overlay.hpp"

#include <bdsg/overlays/overlay_helper.hpp>
#include <structures/union_find.hpp>

#include <array>
#include <iostream>

namespace vg {

//#define debug

using namespace std;

class IntegratedSnarlFinder::MergedAdjacencyGraph {
protected:
    /// Hold onto the backing RankedHandleGraph.
    /// Union find index is handle rank - 1 (to make it 0-based).
    const RankedHandleGraph* graph;
    
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
    /// Make a MergedAdjacencyGraph representing the graph of adjacency components of the given RankedHandleGraph.
    MergedAdjacencyGraph(const RankedHandleGraph* graph);
    
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
    void for_each_other_member(handle_t head, const function<void(handle_t)>& iteratee) const;
    
    /// For each item, including the head, in the component headed by the given
    /// handle, calls the iteratee with the that other item. Does not call the
    /// iteratee for single-item components.
    void for_each_member(handle_t head, const function<void(handle_t)>& iteratee) const;
    
    /// For each item other than the head in each component, calls the iteratee
    /// with the head and the other item. Does not call the iteratee for
    /// single-item components.
    void for_each_membership(const function<void(handle_t, handle_t)>& iteratee) const;
    
    /// In a graph where all 3-edge-connected components have had their nodes
    /// merged, find all the cycles. Cycles are guaranteed to overlap at at
    /// most one node, so no special handling of overlapping regions is done.
    ///
    /// Returns a list of cycle edge length (in bp) and an edge on each cycle
    /// (for the longest cycle in each connected component), and a map from
    /// each edge on the cycle to the next edge, going around each cycle in one
    /// direction (for all cycles).
    ///
    /// Ignores self loops.
    pair<vector<pair<size_t, handle_t>>, unordered_map<handle_t, handle_t>> cycles_in_cactus() const;
    
    /// Find a path of cycles connecting two components in a Cactus graph.
    /// Cycles are represented by the handle that brings that cyle into the component where it intersects the previous cycle.
    /// Because the graph is a Cactus graph, cycles are a tree and intersect at at most one node.
    /// Uses the given map of cycles, stored in one orientation only, to traverse cycles.
    vector<handle_t> find_cycle_path_in_cactus(const unordered_map<handle_t, handle_t>& next_along_cycle, handle_t start_cactus_head, handle_t end_cactus_head) const;
    
    /// Return the path length (total edge length in bp) and edges for the
    /// longest path in each tree in a forest. Ignores self loops on tree nodes.
    ///
    /// Also return the map from the head of each component to the edge into
    /// the child that is first along the longest path to a leaf. For
    /// components not themselves on the longest leaf-leaf path in their tree,
    /// these will always be dangling off/rooted by the longest leaf-leaf path
    /// or longest simple cycle merged away, whichever is longer.
    ///
    /// Needs access to the longest simple cycles that were merged out, if any.
    /// If a path in the forest doesn't match or beat the length of the cycle
    /// that lives in its tree, it is omitted.
    pair<vector<pair<size_t, vector<handle_t>>>, unordered_map<handle_t, handle_t>> longest_paths_in_forest(const vector<pair<size_t, handle_t>>& longest_simple_cycles) const;
    
    /// Describe the graph in dot format to the given stream;
    void to_dot(ostream& out) const;
};

void IntegratedSnarlFinder::MergedAdjacencyGraph::to_dot(ostream& out) const {
    out << "digraph G {" << endl;
    for_each_head([&](handle_t head) {
        // Have a node for every head
        out << "\tn" << graph->get_id(head) << (graph->get_is_reverse(head) ? "r" : "f") << "[shape=\"point\"];" << endl;
        for_each_member(head, [&](handle_t edge) {
            // For everything reading into here
            if (graph->get_is_reverse(edge)) {
                // If it is coming here backward (we are its start)
                handle_t flipped_head = find(graph->flip(edge));
                // Only draw edges in one direction. Label them with their nodes.
                out << "\tn" << graph->get_id(head) << (graph->get_is_reverse(head) ? "r" : "f")
                    << " -> n" << graph->get_id(flipped_head) << (graph->get_is_reverse(flipped_head) ? "r" : "f") 
                    << " [label=" << graph->get_id(edge) << "];" << endl;
            }
        });
    });
    out << "}" << endl;
}

size_t IntegratedSnarlFinder::MergedAdjacencyGraph::uf_rank(handle_t into) const {
    // We need to 0-base the backing rank
    return graph->handle_to_rank(into) - 1;
}

handle_t IntegratedSnarlFinder::MergedAdjacencyGraph::uf_handle(size_t rank) const {
    // We need to 1-base the rank and then get the handle.
    return graph->rank_to_handle(rank + 1);
}

IntegratedSnarlFinder::MergedAdjacencyGraph::MergedAdjacencyGraph(const RankedHandleGraph* graph) : graph(graph),
    union_find(graph->get_node_count() * 2, true) {
    
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

void IntegratedSnarlFinder::MergedAdjacencyGraph::for_each_other_member(handle_t head, const function<void(handle_t)>& iteratee) const {
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

void IntegratedSnarlFinder::MergedAdjacencyGraph::for_each_member(handle_t head, const function<void(handle_t)>& iteratee) const {
    size_t head_rank = uf_rank(head);
    // Find the group the head is in
    vector<size_t> group = union_find.group(head_rank);
    for (auto& member_rank : group) {
        // And go through all the members, including the head.
        // This function will happen to work for non-head inputs, leaving out that input, but we don't guarantee it!
        iteratee(uf_handle(member_rank));
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

pair<vector<pair<size_t, handle_t>>, unordered_map<handle_t, handle_t>> IntegratedSnarlFinder::MergedAdjacencyGraph::cycles_in_cactus() const {
    // Do a DFS over all connected components of the graph
    
    // We will fill this in
    pair<vector<pair<size_t, handle_t>>, unordered_map<handle_t, handle_t>> to_return;
    auto& longest_cycles = to_return.first;
    auto& next_edge = to_return.second;
    
    // When we see a back edge that isn't a self loop, we jump up the stack and
    // walk down it, writing the cycle relationships. We know the cycles can't
    // overlap (due to merged 3 edge connected components) so we also know we
    // won't ever re-walk the same part of the stack, so this is efficient.
    
    // To let us actually find the cycle paths when we see back edges, we need
    // to remember stack frame index as our visit marker. No record = not
    // visited, record = visited.
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
    };
    
    vector<DFSFrame> stack;
    
    for_each_head([&](handle_t component_root) {
        // For every node in the graph
        
        if (!visited_frame.count(component_root)) {
        
#ifdef debug
            cerr << "Root simple cycle search at " << graph->get_id(component_root) << (graph->get_is_reverse(component_root) ? "-" : "+") << endl; 
#endif
        
            // If it hasn't been searched yet, start a search of its connected component.
            stack.emplace_back();
            stack.back().here = component_root;
            
            // We'll put the longest cycle start edge result here, or get rid of it if we find no cycle.
            longest_cycles.emplace_back();
            auto& longest_cycle = longest_cycles.back();
            
            while (!stack.empty()) {
                // Until the DFS is done
                auto& frame = stack.back();
                // Find the node that following this edge got us to.
                auto frame_head = find(frame.here);
                
#ifdef debug
                cerr << "At stack frame " << stack.size() - 1 << " for edge " << graph->get_id(frame.here) << (graph->get_is_reverse(frame.here) ? "-" : "+") 
                    << " on component " << graph->get_id(frame_head) << (graph->get_is_reverse(frame_head) ? "-" : "+") << endl; 
#endif
                
                auto frame_it = visited_frame.find(frame_head);
                if (frame_it == visited_frame.end()) {
                    // First visit to here.
                    
#ifdef debug
                    cerr << "\tFirst visit" << endl; 
#endif
                    
                    // Mark visited at this stack level
                    frame_it = visited_frame.emplace_hint(frame_it, frame_head, stack.size() - 1);
                    
                    // Queue up edges
                    for_each_member(frame_head, [&](handle_t member) {
                        if (member != frame.here || stack.size() == 1) {
                            // If it's not just turning around and looking up
                            // the edge we took to get here, or if we're the
                            // top stack frame and we didn't come from anywhere
                            // anyway
                        
                            // Follow edge by flipping. But queue up the edge
                            // followed instead of the node reached (head), so we
                            // can emit the cycle later in terms of edges.
                            frame.todo.push_back(graph->flip(member));
                        
#ifdef debug
                            cerr << "\t\tNeed to follow " << graph->get_id(frame.todo.back()) << (graph->get_is_reverse(frame.todo.back()) ? "-" : "+") << endl; 
#endif
                        }
                    });
                }
                
                if (!frame.todo.empty()) {
                    // Now do an edge
                    handle_t edge_into = frame.todo.back();
                    handle_t connected_head = find(edge_into);
                    frame.todo.pop_back();
                    
#ifdef debug
                    cerr << "\tFollow " << graph->get_id(edge_into) << (graph->get_is_reverse(edge_into) ? "-" : "+")
                        << " to component " << graph->get_id(connected_head) << (graph->get_is_reverse(connected_head) ? "-" : "+") << endl; 
#endif
                    
                    auto connected_it = visited_frame.find(connected_head);
                    
                    if (connected_it == visited_frame.end()) {
                    
#ifdef debug
                        cerr << "\t\tNot yet visited. Recurse!" << endl; 
#endif
                    
                        // Forward edge. Recurse.
                        // TODO: this immediately does a lookup in the hash table again.
                        stack.emplace_back();
                        stack.back().here = edge_into;
                    } else {
                        // Back edge
                        if (frame_it->second > connected_it->second) {
                            // We have an edge to something that was visited above
                            // our stack level. It can't be a self loop, and it
                            // must close a unique cycle.
                            
#ifdef debug
                            cerr << "\tBack edge up stack to frame " << connected_it->second << endl; 
#endif
                        
#ifdef debug
                            cerr << "\t\tFound cycle:" << endl; 
#endif
                            
                            // Walk and measure the cycle. But don't count the
                            // frame we arrived at because its incoming edge
                            // isn't actually on the cycle.
                            size_t cycle_length_bp = graph->get_length(edge_into);
                            handle_t prev_edge = edge_into;
                            for (size_t i = connected_it->second + 1; i < stack.size(); i++) {
                                // For each edge along the cycle...
                                
#ifdef debug
                                cerr << "\t\t\t" << graph->get_id(stack[i].here) << (graph->get_is_reverse(stack[i].here) ? "-" : "+") << endl; 
#endif
                                
                                // Measure it
                                cycle_length_bp += graph->get_length(stack[i].here);
                                // Record the cycle membership
                                next_edge[prev_edge] = stack[i].here;
                                // Advance
                                prev_edge = stack[i].here;
                            }
                            // Close the cycle
                            next_edge[prev_edge] = edge_into;
                            
#ifdef debug
                            cerr << "\t\t\t" << graph->get_id(edge_into) << (graph->get_is_reverse(edge_into) ? "-" : "+") << endl; 
#endif
                            
#ifdef debug
                            cerr << "\t\tCycle length: " << cycle_length_bp << " bp" << endl; 
#endif
                            
                            if (cycle_length_bp > longest_cycle.first) {
                                // New longest cycle (or maybe only longest cycle).
                                
#ifdef debug
                                cerr << "\t\t\tNew longest cycle!" << endl; 
#endif
                                
                                // TODO: Assumes no cycles are 0-length
                                longest_cycle.first = cycle_length_bp;
                                longest_cycle.second = edge_into;
                            }
                        }
                    }
                } else {
                    // Now we're done with this stack frame.
                    
                    // Clean up
                    stack.pop_back();
                }
            }
            
            if (longest_cycle.first == 0) {
                // No (non-empty) nontrivial cycle found in this connected component.
                // Remove its spot.
                longest_cycles.pop_back();
            }
        }
    });
   
#ifdef debug
    cerr << "Cycle links:" << endl;
    for (auto& kv : next_edge) {
        cerr << "\t" << graph->get_id(kv.first) << (graph->get_is_reverse(kv.first) ? "-" : "+")
            << " -> " << graph->get_id(kv.second) << (graph->get_is_reverse(kv.second) ? "-" : "+") << endl;
    }
#endif
   
    return to_return;
}

vector<handle_t> IntegratedSnarlFinder::MergedAdjacencyGraph::find_cycle_path_in_cactus(const unordered_map<handle_t, handle_t>& next_along_cycle, handle_t start_head, handle_t end_head) const {
    // We fill this in with a path of cycles.
    // Each cycle is the edge on that cycle leading into the node that it shares with the previous cycle.
    vector<handle_t> cycle_path;
    
    // We just DFS through the cycle tree until we find one that touches the
    // other node. We represent each cycle by an edge on it into the node where
    // it overlaps the parent cycle, and store the current cycle and the other
    // cycles to do that share nodes.
    vector<tuple<handle_t, vector<handle_t>, bool>> cycle_stack;
    
    // We have a list of DFS roots we can stop early on.
    vector<handle_t> roots;
    for_each_member(start_head, [&](handle_t inbound) {
        if (next_along_cycle.count(inbound)) {
            // This edge is how a cycle comes into here. Consider this cycle.
            roots.push_back(inbound);
        }
    });
    
    for (auto& root : roots) {
        // Root at each root
        cycle_stack.emplace_back(root, vector<handle_t>(), false);
        while (!cycle_stack.empty()) {
            auto& cycle_frame = cycle_stack.back();
            if (!get<2>(cycle_frame)) {
                // First visit
                get<2>(cycle_frame) = true;
                
                // Need to fill in child cycles.
                for (auto it = next_along_cycle.find(get<0>(cycle_frame)); it->second != get<0>(cycle_frame); it = next_along_cycle.find(it->second)) {
                    // For each other edge around the cycle (in it->second) other than the one we started at
                    
                    handle_t node = find(it->second);
                    if (node == end_head) {
                        // This cycle intersects the destination. It is the last on the cycle path.
                        
                        // Copy the path on the stack over.
                        // Note that the first think on the path is in the
                        // start's component, but the last thing on the path
                        // isn't in the end's component.
                        cycle_path.reserve(cycle_stack.size());
                        for (auto& f : cycle_stack) {
                            cycle_path.push_back(get<0>(f));
                        }
                        
                        // Now the cycle path is done
                        return cycle_path;
                    }
                    
                    for_each_member(node, [&](handle_t inbound) {
                        // For each edge in the component it enters
                        if (inbound != it->second && next_along_cycle.count(inbound)) {
                            // This edge is a cycle coming into a node our current cycle touches.
                            get<1>(cycle_frame).push_back(inbound);
                        }
                    });
                }
            }
            if (!get<1>(cycle_frame).empty()) {
                // Need to recurse on a connected cycle.
                handle_t child = get<1>(cycle_frame).back();
                get<1>(cycle_frame).pop_back();
                cycle_stack.emplace_back(child, vector<handle_t>(), false);
            } else {
                // Need to clean up and return
                cycle_stack.pop_back();
            }
        }
    }
    
    // If we get here, we never found a path.
    // Complain! Something is wrong!
    throw runtime_error("Cound not find cycle path!");
}

pair<vector<pair<size_t, vector<handle_t>>>, unordered_map<handle_t, handle_t>> IntegratedSnarlFinder::MergedAdjacencyGraph::longest_paths_in_forest(
    const vector<pair<size_t, handle_t>>& longest_simple_cycles) const {
    
    // TODO: somehow unify DFS logic with cycle-finding DFS in a way that still
    // allows us to inspect our stack in each case?
    
    // Going up the tree, we need to track the longest path from a leaf to the
    // subtree root, and the longest path between leaves in the subtree. These
    // ought to overlap substantially, but either may be the real winner when
    // we get to where we rooted the DFS.
    
    // Set up the return value
    pair<vector<pair<size_t, vector<handle_t>>>, unordered_map<handle_t, handle_t>> to_return;
    
    // When we find a longest path in a connected component (tree), we put its
    // length and value in here. We describe it as edges followed.
    auto& longest_tree_paths = to_return.first;
    
    // We use this as part of our DFS scratch to record the first edge on the
    // deepest path to a leaf in a subtree. The actual length of that path is
    // stored in the main record for the head of the component the given edge
    // reaches. If we find a longest leaf-leaf path in the tree that beats the
    // simple cycle (if any), we rewrite this to be rooted somewhere along that
    // leaf-leaf path. Indexed by head.
    auto& deepest_child_edge = to_return.second;
    
    // The DFS also needs records, one per component, indexed by head.
    // Absence of a record = unvisited.
    struct DFSRecord {
        // Remember the edge to traverse to get back to the parent, so we can
        // find the path from the longest leaf-leaf path's converging node to
        // the DFS root if we need it.
        handle_t parent_edge;
        // How long is the deepest path to a leaf from here, plus the length of
        // the edge followed to here from the parent?
        // Filled in when we leave the stack, by looking at deepest_child_edge.
        size_t leaf_path_length = 0;
        // What edge goes to the second-deepest child, if we have one, to form
        // the longest leaf-leaf path converging here?
        handle_t second_deepest_child_edge;
        // And do we have such a second-deepest child?
        bool has_second_deepest_child = false;
        // And what head in the graph is the convergance point of the longest
        // leaf-leaf path in our subtree? If it points to us, and we don't have
        // a second deepest child, there is no leaf-leaf path in our subtree.
        //
        // Actually filled in when we finish a node. When children write to it
        // to max themselves in, they clobber it if it points back to us.
        handle_t longest_subtree_path_root;
        // We don't need to store this, because it's determined by the leaf
        // path lengths of the best and second best children of the longest
        // subtree path root, but to save a whole mess of transitive accesses
        // we track the longest subtree paht length here as well. Will be 0
        // when there is no subtree leaf-leaf path.
        size_t longest_subtree_path_length;
    };
    unordered_map<handle_t, DFSRecord> records;
    
    
    // We need a stack.
    // Stack is actually in terms of inward edges followed.
    struct DFSFrame {
        handle_t here;
        // What edges still need to be followed
        vector<handle_t> todo;
    };
    
    vector<DFSFrame> stack;
    
    // We have a function to try DFS from a root, if the root is unvisited.
    // If root_cycle_length is nonzero, we will not rewrite deepest_child_edge
    // to point towards the longest leaf-leaf path, if it isn't as long as the
    // cycle or longer.
    auto try_root = [&](handle_t traversal_root, size_t root_cycle_length) {
        if (!records.count(traversal_root)) {
            // If it hasn't been searched yet, start a search
            stack.emplace_back();
            stack.back().here = traversal_root;
            
#ifdef debug
            cerr << "Root bridge tree traversal at " << graph->get_id(traversal_root) << (graph->get_is_reverse(traversal_root) ? "-" : "+") << endl;
#endif
            
            while (!stack.empty()) {
                // Until the DFS is done
                auto& frame = stack.back();
                // Find the node that following this edge got us to.
                auto frame_head = find(frame.here);
                
#ifdef debug
                cerr << "At stack frame " << stack.size() - 1 << " for edge " << graph->get_id(frame.here) << (graph->get_is_reverse(frame.here) ? "-" : "+") 
                    << " into component with head " << graph->get_id(frame_head) << (graph->get_is_reverse(frame_head) ? "-" : "+") << endl;
#endif
                
                auto frame_it = records.find(frame_head);
                if (frame_it == records.end()) {
                    // First visit to here.
                    
#ifdef debug
                cerr << "\tFirst visit. Find edges." << endl;
#endif
                    
                    // Mark visited
                    frame_it = records.emplace_hint(frame_it, frame_head, DFSRecord());
                    // And fill it in with default references.
                    // Remember how to get back to the parent
                    frame_it->second.parent_edge = graph->flip(frame.here);
                    // Say there's no known leaf-leaf path converging anywhere under it yet.
                    frame_it->second.longest_subtree_path_root = frame_head;
                    
                    // Queue up edges
                    for_each_member(frame_head, [&](handle_t member) {
                        // Follow edge by flipping.
                        auto flipped = graph->flip(member);
                        
                        if (find(flipped) != frame_head) {
                            // Only accept non-self-loops.
                        
#ifdef debug
                            cerr << "\t\tNeed to follow " << graph->get_id(flipped) << (graph->get_is_reverse(flipped) ? "-" : "+") << endl;
#endif
                        
                            // Queue up the edge followed instead of the node
                            // reached (head), so we can emit the cycle later
                            // in terms of edges.
                            frame.todo.push_back(flipped);
                        }
                    });
                }
                
                auto& record = frame_it->second;
                
                if (!frame.todo.empty()) {
                    // Now do an edge
                    handle_t edge_into = frame.todo.back();
                    handle_t connected_head = find(edge_into);
                    frame.todo.pop_back();
                    
#ifdef debug
                    cerr << "\tFollowing " << graph->get_id(edge_into) << (graph->get_is_reverse(edge_into) ? "-" : "+") << endl;
#endif
                    
                    if (!records.count(connected_head)) {
                        // Forward edge. Recurse.
                        
#ifdef debug
                        cerr << "\t\tReaches unvisited " << graph->get_id(connected_head) << (graph->get_is_reverse(connected_head) ? "-" : "+") << "; Recurse!" << endl;
#endif
                        
                        stack.emplace_back();
                        stack.back().here = edge_into;
                    }
                } else {
                    // No children left.
                    
#ifdef debug
                    cerr << "\tDone with all children." << endl;
#endif
                    
                    // Did any of our children decalre themselves deepest?
                    // Or do we have no children.
                    auto deepest_child_edge_it = deepest_child_edge.find(frame_head);
                    
                    if (stack.size() > 1) {
                        // If we have a parent
                        auto& parent_frame = stack[stack.size() - 2];
                        auto parent_head = find(parent_frame.here);
                        auto& parent_record = records[parent_head];
                        
                        // The length of the path to a leaf will involve the edge from the parent to here.
                        record.leaf_path_length = graph->get_length(frame.here);
                        
#ifdef debug
                        cerr << "\t\tLength of path to deepest leaf is " << record.leaf_path_length << " bp" << endl;
#endif
                        
                        if (deepest_child_edge_it != deepest_child_edge.end()) {
                            // And if we have a child to go on with, we add the length of that path
                            record.leaf_path_length += records[find(deepest_child_edge_it->second)].leaf_path_length;
                            
#ifdef debug
                                cerr << "\t\t\tPlus length from here to leaf via "
                                    << graph->get_id(deepest_child_edge_it->second) << (graph->get_is_reverse(deepest_child_edge_it->second) ? "-" : "+")
                                    << " for " << record.leaf_path_length << " bp total" << endl;
#endif
                            
                        }
                        
                        // Fill in deepest_child_edge for the parent if not filled in already, or if we beat what's there.
                        // Also maintain parent's second_deepest_child_edge.
                        auto parent_deepest_child_it = deepest_child_edge.find(parent_head);
                        if (parent_deepest_child_it == deepest_child_edge.end()) {
                        
#ifdef debug
                            cerr << "\t\tWe are our parent's deepest child by default!" << endl;
#endif
                        
                            // Emplace in the map where we didn't find anything.
                            deepest_child_edge.emplace_hint(parent_deepest_child_it, parent_head, frame.here);
                        } else if(records[find(parent_deepest_child_it->second)].leaf_path_length < record.leaf_path_length) {
                            // We are longer than what's there now
                            
#ifdef debug
                            cerr << "\t\tWe are our parent's new deepest child!" << endl;
#endif
                            
                            // Demote what's there to second-best
                            parent_record.second_deepest_child_edge = parent_deepest_child_it->second;
                            parent_record.has_second_deepest_child = true;
                            
#ifdef debug
                            cerr << "\t\t\tWe demote "
                                << graph->get_id(parent_record.second_deepest_child_edge) << (graph->get_is_reverse(parent_record.second_deepest_child_edge) ? "-" : "+")
                                << " to second-deepest child" << endl;
#endif
                            
                            // Replace the value we found
                            parent_deepest_child_it->second = frame.here;
                        } else if (!parent_record.has_second_deepest_child) {
                            
#ifdef debug
                            cerr << "\t\tWe are our parent's second deepest child by default!" << endl;
#endif
                            
                            // There's no second-deepest recorded so we must be it.
                            parent_record.second_deepest_child_edge = frame.here;
                            parent_record.has_second_deepest_child = true;
                        } else if (records[find(parent_record.second_deepest_child_edge)].leaf_path_length < record.leaf_path_length) {
                            
#ifdef debug
                            cerr << "\t\tWe are our parent's new second deepest child!" << endl;
#endif
                            
                            // We are a new second deepest child.
                            parent_record.second_deepest_child_edge = frame.here;
                        }
                    }
                    
                    // The length of the longest leaf-leaf path converging at or under any child (if any) is in record.longest_subtree_path_length.
                    
                    if (record.has_second_deepest_child || stack.size() == 1) {
                        // If there's a second incoming leaf path there's a converging leaf-leaf path here.
                        // If we're the root and there *isn't* a second incoming leaf-leaf path, we are ourselves a leaf.
                        
                        // Grab the length of the longest leaf-leaf path converging exactly here.
                        // TODO: can we not look up the deepest child's record again?
                        size_t longest_here_path_length = 0;
                        if (deepest_child_edge_it != deepest_child_edge.end()) {
                            longest_here_path_length += records[find(deepest_child_edge_it->second)].leaf_path_length;
                        }
                        if (record.has_second_deepest_child) {
                            longest_here_path_length += records[find(record.second_deepest_child_edge)].leaf_path_length;
                        }
                        
#ifdef debug
                        cerr << "\t\tPaths converge here with total length " << longest_here_path_length << " bp" << endl;
#endif
                        
                        if (record.longest_subtree_path_root == frame_head || longest_here_path_length > record.longest_subtree_path_length) {
                            
#ifdef debug
                            cerr << "\t\t\tNew longest path in subtree!" << endl;
#endif
                            
                            // If there's no path from a child, or this path is
                            // longer, set record.longest_subtree_path_root
                            // (back) to frame_head to record that.
                            record.longest_subtree_path_root = frame_head;
                            // And save the length.
                            record.longest_subtree_path_length = longest_here_path_length;
                            
                            // Now we are the new root of the longest leaf-leaf path converging at or under us.
                        }
                    }
                    
                    if (stack.size() > 1 && record.longest_subtree_path_length > 0) {
                        // We have a leaf-leaf path converging at or under here, and we have a parent.
                        // TODO: we assume leaf-leaf paths are nonzero length here.
                        // TODO: save searching up the parent record again
                        auto& parent_frame = stack[stack.size() - 2];
                        auto parent_head = find(parent_frame.here);
                        auto& parent_record = records[parent_head];
                        
                        // Max our longest leaf-leaf path in against the paths contributed by previous children.
                        if (parent_record.longest_subtree_path_root == parent_head ||
                            parent_record.longest_subtree_path_length < record.longest_subtree_path_length) {
                            
#ifdef debug
                            cerr << "\t\tLongest path in our subtree converging at "
                                << graph->get_id(record.longest_subtree_path_root) << (graph->get_is_reverse(record.longest_subtree_path_root) ? "-" : "+")
                                << " is the new longest path in our parent's subtree." << endl;
#endif
                            
                            // No child has contributed their leaf-leaf path so far, or ours is better.
                            parent_record.longest_subtree_path_root = record.longest_subtree_path_root;
                            parent_record.longest_subtree_path_length = record.longest_subtree_path_length;
                        }
                    }
                   
                    if (stack.size() == 1) {
                        // When we get back to the root
                        
#ifdef debug
                        cerr << "\t\tWe were the root of the traversal." << endl;
#endif

                        if (record.longest_subtree_path_length >= root_cycle_length) {
                            // Either we didn't root at a cycle, or we found a longer leaf-leaf path that should be the decomposition root instead.
                            
#ifdef debug
                            cerr << "\t\t\tTree has leaf-leaf path that is as long as or longer than any cycle at root ("
                                << record.longest_subtree_path_length << "bp)." << endl;
#endif
                            
                            // We need to record the longest tree path.
                            longest_tree_paths.emplace_back();
                            longest_tree_paths.back().first = record.longest_subtree_path_length;
                            auto& path = longest_tree_paths.back().second;
                            
                            auto& path_root_frame = records[record.longest_subtree_path_root];
                            
                            if (path_root_frame.has_second_deepest_child) {
                                // This is an actual convergence point
                            
#ifdef debug
                                cerr << "\t\t\t\tConverges at real convergence point" << endl;
#endif
                                
                                // Collect the whole path down the second deepest child
                                path.push_back(path_root_frame.second_deepest_child_edge);
                                auto path_trace_it = deepest_child_edge.find(find(path.back()));
                                while (path_trace_it != deepest_child_edge.end()) {
                                    // Follow the deepest child relationships until they run out.
                                    path.push_back(path_trace_it->second);
                                    path_trace_it = deepest_child_edge.find(find(path.back()));
                                }
                                // Reverse what's there and flip all the edges
                                vector<handle_t> flipped;
                                flipped.reserve(path.size());
                                for (auto path_it = path.rbegin(); path_it != path.rend(); ++path_it) {
                                    flipped.push_back(graph->flip(*path_it));
                                }
                                path = std::move(flipped);
                            } else {
                                // There's no second-longest path; we statted at one of the most distant leaves.
#ifdef debug
                                cerr << "\t\t\t\tConverges at leaf" << endl;
#endif
                            }
                            
                            if (deepest_child_edge.count(record.longest_subtree_path_root)) {
                            
#ifdef debug
                                cerr << "\t\t\t\tNonempty path to distinct other leaf" << endl;
#endif
                            
                                // There's a nonempty path to another leaf,
                                // other than the furthest one (and we aren't
                                // just a point).
                                // Trace the actual longest path from root to leaf and add it on
                                path.push_back(deepest_child_edge[record.longest_subtree_path_root]);
                                auto path_trace_it = deepest_child_edge.find(find(path.back()));
                                while (path_trace_it != deepest_child_edge.end()) {
                                    // Follow the deepest child relationships until they run out.
                                    path.push_back(path_trace_it->second);
                                    path_trace_it = deepest_child_edge.find(find(path.back()));
                                }
                            }
                            
                            // OK now we have the longest leaf-leaf path saved.
                            
                            // We need to redo the path from the tree traversal
                            // root to the longest path convergence point, to
                            // fix up the subtree rooting information.
                            
                            // Go to the convergence point
                            handle_t cursor = record.longest_subtree_path_root;
                            
                            // This will be the path of edges to take from the convergence point (new root) to the traversal root (old root)
                            vector<handle_t> convergence_to_old_root;
                            while (cursor != frame_head) {
                                // Walk up the parent pointers to the traversal root and stack up the heads.
                                // We may get nothing if the root happened to already be on the longest leaf-leaf path.
                                auto& cursor_record = records[cursor];
                                convergence_to_old_root.push_back(cursor_record.parent_edge);
                                cursor = find(cursor_record.parent_edge);
                            }
                            
#ifdef debug
                            cerr << "\t\t\t\tRewrite along " << convergence_to_old_root.size() << " edges..." << endl;
#endif
                            
                            while (!convergence_to_old_root.empty()) {
                                // Then go down that stack
                                
                                // Define new child and parent
                                handle_t parent_child_edge = convergence_to_old_root.back();
                                handle_t child_head = find(parent_child_edge);
                                handle_t parent_head = find(graph->flip(parent_child_edge));
                                
                                // TODO: find a way to demote parent to child here on each iteration
                                auto& child_record = records[child_head];
                                auto& parent_record = records[parent_head];
                                
                                // If the deepest child of the child is actually the parent, disqualify it
                                deepest_child_edge_it = deepest_child_edge.find(child_head);
                                
                                if (deepest_child_edge_it != deepest_child_edge.end() && find(deepest_child_edge_it->second) == parent_head) {
                                    // The parent was the child's deepest child. Can't have that.
                                    if (child_record.has_second_deepest_child) {
                                        // Promote the second deepest child.
                                        deepest_child_edge_it->second = child_record.second_deepest_child_edge;
                                        child_record.has_second_deepest_child = false;
                                    } else {
                                        // No more deepest child
                                        deepest_child_edge.erase(deepest_child_edge_it);
                                        deepest_child_edge_it = deepest_child_edge.end();
                                    }
                                }
                                
                                // The child may not have had a parent before.
                                // So we need to fill in its longest leaf path
                                // length counting its new parent edge.
                                
                                // But we know all its children are done.
                                
                                // The length of the path to a leaf will involve the edge from the parent to the child
                                child_record.leaf_path_length = graph->get_length(parent_child_edge);
                                
                                if (deepest_child_edge_it != deepest_child_edge.end()) {
                                    // And if we have a child to go on with, we add the length of that path
                                    child_record.leaf_path_length += records[find(deepest_child_edge_it->second)].leaf_path_length;
                                }
                                
                                // Now we have to mix ourselves into the parent.
                                // We do it the same way as normal. Both the deepest and second-deepest child of the parent can't be the grandparent.
                                // So if they both beat us we can't be the real deepest child.
                                // If we beat the second deepest child, and the original deepest child gets disqualified for being the grandparent, we become the parent's deepest child.
                                // And if we beat both it doesn't matter whether either gets disqualified, because we win.
                                
                                // TODO: deduplicate code with the original DFS?
                                
                                // Fill in deepest_child_edge for the parent if not filled in already, or if we beat what's there.
                                // Also maintain parent's second_deepest_child_edge.
                                auto parent_deepest_child_it = deepest_child_edge.find(parent_head);
                                if (parent_deepest_child_it == deepest_child_edge.end()) {
                                    // Emplace in the map where we didn't find anything.
                                    deepest_child_edge.emplace_hint(parent_deepest_child_it, parent_head, parent_child_edge);
                                } else if(records[find(parent_deepest_child_it->second)].leaf_path_length < child_record.leaf_path_length) {
                                    // We are longer than what's there now
                                    
                                    // Demote what's there to second-best
                                    parent_record.second_deepest_child_edge = parent_deepest_child_it->second;
                                    parent_record.has_second_deepest_child = true;
                                    
                                    // Replace the value we found
                                    parent_deepest_child_it->second = parent_child_edge;
                                } else if (!parent_record.has_second_deepest_child) {
                                    // There's no second-deepest recorded so we must be it.
                                    parent_record.second_deepest_child_edge = parent_child_edge;
                                    parent_record.has_second_deepest_child = true;
                                } else if (records[find(parent_record.second_deepest_child_edge)].leaf_path_length < child_record.leaf_path_length) {
                                    // We are a new second deepest child.
                                    parent_record.second_deepest_child_edge = parent_child_edge;
                                }
                               
                                // Now the new child, if its path is deep enough, is the parent's new deepest or second deepest child edge.
                                // Go up a level, disqualify the grandparent, and see who wins.
                                // The new root has no parent itself, so all its edges are eligible and the one with the longest path wins.
                                convergence_to_old_root.pop_back();
                            }
                            
#ifdef debug
                            for (auto& item : path) {
                                cerr << "\t\t\t\tPath visits: "
                                    << graph->get_id(item) << (graph->get_is_reverse(item) ? "-" : "+")
                                    << " length " << graph->get_length(item) << endl;
                            }
#endif

                            if (path.empty()) {
                                // If the leaf-leaf path is empty, stick in a handle so we can actually find the single leaf in the bridge forest.
                                assert(longest_tree_paths.back().first == 0);
                                path.push_back(traversal_root);
                            } else {
                                // If anything is on the path, we shouldn't have 0 length.
                                assert(longest_tree_paths.back().first != 0);
                            }
                            
                        }
                    }
                    
                    // Now we're done with this stack frame.
                       
                    // Clean up
                    stack.pop_back();
                }
            }
        }
    };
    
    for (auto it = longest_simple_cycles.begin(); it != longest_simple_cycles.end(); ++it) {
        // Try it from the head of the component that each longest input simple
        // cycle got merged into. If we end up using that longest cycle to root
        // this component, we will have everything pointing the right way
        // already.
        try_root(find(it->second), it->first);
    }
    
    // And then try it on every head in general to mop up anything without a simple cycle in it
    for_each_head([&](handle_t head) {
        try_root(head, 0);
    });
    
    // The DFS records die with this function, but the rewritten deepest child
    // edges survive and let us root snarls having only their incoming ends.
    // And we have all the longest tree paths that beat their components
    // rooting cycles, if any.
    
#ifdef debug
    cerr << "Edges to deepest children in bridge forest:" << endl;
    for (auto& kv : deepest_child_edge) {
        cerr << "\t" << graph->get_id(kv.first) << (graph->get_is_reverse(kv.first) ? "-" : "+")
            << " -> " << graph->get_id(kv.second) << (graph->get_is_reverse(kv.second) ? "-" : "+") << endl;
    }
#endif
    
    return to_return;
}





////////////////////////////////////////////////////////////////////////////////////////////




IntegratedSnarlFinder::IntegratedSnarlFinder(const HandleGraph& graph) : HandleGraphSnarlFinder(&graph) {
    // Nothing to do!
}

void IntegratedSnarlFinder::traverse_decomposition(const function<void(handle_t)>& begin_chain, const function<void(handle_t)>& end_chain,
    const function<void(handle_t)>& begin_snarl, const function<void(handle_t)>& end_snarl) const {
    
    // Do the actual snarl finding work and then walk the bilayered tree.
    
#ifdef debug
    cerr << "Ranking graph handles." << endl;
#endif
    
    // First we need to ensure that our graph has dense handle ranks
    bdsg::RankedOverlayHelper overlay_helper;
    auto ranked_graph = overlay_helper.apply(graph);
    
#ifdef debug
    cerr << "Finding snarls." << endl;
#endif
    
    // We need a union-find over the adjacency components of the graph, in which we will build the cactus graph.
    MergedAdjacencyGraph cactus(ranked_graph);
    
#ifdef debug
    cerr << "Base adjacency components:" << endl;
    cactus.to_dot(cerr);
#endif
    
    // It magically gets the adjacency components itself.
    
#ifdef debug
    cerr << "Finding 3 edge connected components..." << endl;
#endif
    
    // Now we need to do the 3 edge connected component merging, using Tsin's algorithm.
    // We don't really have a good dense rank space on the adjacency components, so we use the general version.
    // TODO: Somehow have a nice dense rank space on components. Can we just use backing graph ranks and hope it's dense enough?
    // We represent each adjacency component (node) by its heading handle.
#ifdef debug
    size_t tecc_id = 0;
#endif
    // Buffer merges until the algorithm is done.
    vector<pair<handle_t, handle_t>> merge_list;
    algorithms::three_edge_connected_component_merges<handle_t>([&](const function<void(handle_t)>& emit_node) {
        // Feed all the handles that head adjacency components into the algorithm
        cactus.for_each_head([&](handle_t head) {
#ifdef debug
            cerr << "Three edge component node " << tecc_id << " is head " << graph->get_id(head) << (graph->get_is_reverse(head) ? "-" : "+") << endl;
            tecc_id++;
#endif
            emit_node(head);
        });
    }, [&](handle_t node, const function<void(handle_t)>& emit_edge) {
        // When asked for edges, don't deduplicate or filter. We want all multi-edges.
        cactus.for_each_member(node, [&](handle_t other_member) {
            // For each handle in the adjacency component that this handle is heading (including the head)
            
            // Follow as an edge again, by flipping
            handle_t member_connected_head = cactus.find(graph->flip(other_member));
            
            if (member_connected_head == node && graph->get_is_reverse(other_member)) {
                // For self loops, only follow them in one direction. Skip in the other.
                return;
            }
            
            // Announce it. Multi-edges are OK.
            emit_edge(member_connected_head);
        });
    }, [&](handle_t a, handle_t b) {
        // Now we got a merge to create the 3 edge connected components.
        // We can't actually do the merge now, because we can't let the merges
        // be visible to the algorithm while it is working. 
        merge_list.emplace_back(a, b);
    });
    
    // Now execute the merges, since the algorithm is done looking at the graph.
    for (auto& ab : merge_list) {
        cactus.merge(ab.first, ab.second);
    }
    merge_list.clear();
    
    // Now our 3-edge-connected components have been condensed, and we have a proper Cactus graph.
    
#ifdef debug
    cerr << "After 3ecc merging:" << endl;
    cactus.to_dot(cerr);
#endif
    
#ifdef debug
    cerr << "Creating bridge forest..." << endl;
#endif
    
    // Then we need to copy the base Cactus graph so we can make the bridge forest
    MergedAdjacencyGraph forest(cactus);
    
#ifdef debug
    cerr << "Finding simple cycles..." << endl;
#endif
    
    // Get cycle information: longest cycle in each connected component, and next edge along cycle for each edge (in one orientation)
    pair<vector<pair<size_t, handle_t>>, unordered_map<handle_t, handle_t>> cycles = cactus.cycles_in_cactus();
    auto& longest_cycles = cycles.first;
    auto& next_along_cycle = cycles.second;
    
    for (auto& kv : next_along_cycle) {
        // Merge along all cycles in the bridge forest
        forest.merge(kv.first, kv.second);
    }

#ifdef debug
    cerr << "Bridge forest:" << endl;
    forest.to_dot(cerr);
#endif
    
#ifdef debug
    cerr << "Finding bridge edge paths..." << endl;
#endif
    
    // Now we find the longest path in each tree in the bridge forest, with its
    // length in bases.
    //
    // We also find, for each bridge edge component head, the edge towards the
    // deepest bridge edge tree leaf, which lets us figure out how to root
    // dangly bits into chains.
    //
    // Make sure to root at the nodes corresponding to the collapsed longest
    // cycles, if the leaf-leaf paths don't win their components.
    //
    // For empty leaf-leaf paths, will emit a single node "path" with a length
    // of 0.
    pair<vector<pair<size_t, vector<handle_t>>>, unordered_map<handle_t, handle_t>> forest_paths = forest.longest_paths_in_forest(longest_cycles);
    auto& longest_paths = forest_paths.first;
    auto& towards_deepest_leaf = forest_paths.second;
    
#ifdef debug
    cerr << "Sorting candidate roots..." << endl;
#endif
    
    // Make sure we are looking at all the cycles and leaf-leaf paths in order.
    // We need basically priority queue between them. But we don't need insert so we jsut sort.
    std::sort(longest_cycles.begin(), longest_cycles.end());
    std::sort(longest_paths.begin(), longest_paths.end());
    
    // Now that we have computed the graphs we need, do the traversal of them.
    // This modifies the structures we have computed in-place, and also does
    // some extra merges in the cactus graph to make all chains cycles.
    traverse_computed_decomposition(cactus,
                                    forest,
                                    longest_paths,
                                    towards_deepest_leaf,
                                    longest_cycles,
                                    next_along_cycle,
                                    begin_chain,
                                    end_chain,
                                    begin_snarl,
                                    end_snarl);
}

/**
 * A set over the nodes in a handle graph.
 * All queries automatically ignore orientation.
 */
class HandleGraphNodeSet {
private:
    unordered_set<handle_t> visited;
    const HandleGraph* graph;
public:
    /**
     * Make a new set over the nodes of the given graph.
     */
    inline HandleGraphNodeSet(const HandleGraph* graph): graph(graph) {
        // Nothing to do
    }
    
    /**
     * Get the number of nodes in the set.
     */
    inline size_t size() const {
        return visited.size();
    }
    
    /**
     * Add a node to the set, given a handle to either orientation.
     */
    inline void insert(const handle_t& here) {
        visited.insert(graph->forward(here));
    }
    
    /**
     * Return whether a node is in the set, given a handle to either orientation.
     */
    inline bool count(const handle_t& here) const {
        return visited.count(graph->forward(here));
    }
};

void IntegratedSnarlFinder::traverse_computed_decomposition(MergedAdjacencyGraph& cactus,
    const MergedAdjacencyGraph& forest,
    vector<pair<size_t, vector<handle_t>>>& longest_paths,
    unordered_map<handle_t, handle_t>& towards_deepest_leaf,
    vector<pair<size_t, handle_t>>& longest_cycles,
    unordered_map<handle_t, handle_t>& next_along_cycle,
    const function<void(handle_t)>& begin_chain, const function<void(handle_t)>& end_chain,
    const function<void(handle_t)>& begin_snarl, const function<void(handle_t)>& end_snarl) const {
  
    // Now, keep a set of all the edges that have found a place in the decomposition.
    // Ignore handle orientation.
    // Because we don't want to mess up orientations, we only access the set through accessors.
    HandleGraphNodeSet visited(graph);
    
#ifdef debug
    cerr << "Traversing cactus graph..." << endl;
#endif

    // How many handle graph nodes need to be decomposed?
    size_t to_decompose = graph->get_node_count();
    while(visited.size() < to_decompose) {
        // While we haven't touched everything
        
#ifdef debug
        if (!longest_cycles.empty()) {
            cerr << "Longest cycle: " << longest_cycles.back().first << " bp" << endl;
        }
        
        if (!longest_paths.empty()) {
            cerr << "Longest path: " << longest_paths.back().first << " bp" << endl;
        }
#endif
        
        // We have a stack.
        struct SnarlChainFrame {
            // Set to true if this is a snarl being generated, and false if it is a chain.
            bool is_snarl = true;
            
            // Set to true if the children have already been enumerated.
            // If we get back to a frame, and this is true, and todo is empty, we are done with the frame.
            bool saw_children = false;
            
            // Into and out-of edges of this snarl or chain, within its parent.
            // Only set if we aren't the root frame on the stack.
            pair<handle_t, handle_t> bounds;
            
            // Edges denoting children to process.
            // If we are a snarl, an entry may be a bridge edge reading into us.
            // If so, we will transform it into a cycle.
            // If we are a snarl, an entry may be a cycle edge reading into us (with the next edge around the cycle reading out).
            // If so, we will recurse on the chain.
            // If we are a chain, an entry may be an edge reading into a child snarl.
            // If so, we will find the other side of the snarl and recurse on the snarl.
            vector<handle_t> todo;
        };
        vector<SnarlChainFrame> stack;
        
        if (longest_cycles.empty() || (!longest_paths.empty() && longest_cycles.back().first <= longest_paths.back().first)) {
            // There should be a path still
            assert(!longest_paths.empty());
            // It should not be empty. It should at least have a single bridge forest node to visit.
            assert(!longest_paths.back().second.empty());
            
            // We will root on a tip-tip path for its connected component, if
            // not already covered, because there isn't a longer cycle.
            
            if (!visited.count(longest_paths.back().second.front())) {
                // This connected component isn't already covered.
                
                handle_t first_edge = longest_paths.back().second.front();
                
                if (longest_paths.back().first == 0) {
                    // This is a 0-length path, but we want to root the decomposition here.
                    // This bridge tree has no nonempty cycles and no bridge edges. It's just all one adjacency component.
                    // All contents spill out into the root snarl as contained nodes.
                    
#ifdef debug
                    cerr << "Single node bridge tree with no real cycles for "
                        << graph->get_id(first_edge) << (graph->get_is_reverse(first_edge) ? "-" : "+") << endl;
                        
                    cerr << "\tSpilling contents into root snarl." << endl;
#endif
                    
                    cactus.for_each_member(cactus.find(first_edge), [&](handle_t inbound) {
                        // The contents are all self loops
                        assert(cactus.find(inbound) == cactus.find(graph->flip(inbound)));
                        if (!graph->get_is_reverse(inbound)) {
                            // We only want them forward so each becomes only one empty chain.
                        
#ifdef debug
                            cerr << "\t\tContain edge " << graph->get_id(inbound) << (graph->get_is_reverse(inbound) ? "-" : "+") << endl;
#endif
                        
                            begin_chain(inbound);
                            end_chain(inbound);
                            
                            visited.insert(inbound);
                        }
                    });
                } else {
                
                    // This is a real path between distinct bridge edge tree leaves
               
#ifdef debug
                    cerr << "Rooting component at tip-tip path starting with " << graph->get_id(first_edge) << (graph->get_is_reverse(first_edge) ? "-" : "+") << endl;
#endif
                    
                    for (size_t i = 1; i < longest_paths.back().second.size(); i++) {
                        // Rewrite the deepest bridge graph leaf path map to point from one end of the tip-tip path to the other
                        // TODO: bump this down into the bridge path finding function
                        
                        handle_t prev_path_edge = longest_paths.back().second[i - 1];
                        handle_t prev_head = forest.find(prev_path_edge);
                        handle_t next_path_edge = longest_paths.back().second[i];
                        
                        towards_deepest_leaf[prev_head] = next_path_edge;
                        
#ifdef debug
                        cerr << "\tEnforce leaf path goes " << graph->get_id(prev_path_edge) << (graph->get_is_reverse(prev_path_edge) ? "-" : "+")
                            << " with head " << graph->get_id(prev_head) << (graph->get_is_reverse(prev_head) ? "-" : "+")
                            << " to next edge " << graph->get_id(next_path_edge) << (graph->get_is_reverse(next_path_edge) ? "-" : "+") << endl;
#endif
                        
                    }
                
                    // Stack up a root/null snarl containing this bridge edge.
                    // Remember to queue it facing inward, toward the new new root at the start of the path.
                    stack.emplace_back();
                    stack.back().is_snarl = true;
                    stack.back().todo.push_back(graph->flip(first_edge));
                    
#ifdef debug
                    cerr << "\tPut cycles and self edges at tip into root snarl" << endl;
#endif
                    
                    // Find all the cycles and self edges that are also here and make sure to do them. Connectivity will be in the root snarl.
                    cactus.for_each_member(cactus.find(graph->flip(first_edge)), [&](handle_t inbound) {
                        if (inbound == graph->flip(first_edge)) {
                            // Skip the one bridge edge we started with
                            return;
                        }
                    
#ifdef debug
                        cerr << "\t\tLook at edge " << graph->get_id(inbound) << (graph->get_is_reverse(inbound) ? "-" : "+") << " on " << next_along_cycle.count(inbound) << " cycles" << endl;
#endif
                    
                        if (next_along_cycle.count(inbound)) {
                            // Put this cycle on the to do list also
                            
#ifdef debug
                            cerr << "\t\t\tLook at cycle edge " << graph->get_id(inbound) << (graph->get_is_reverse(inbound) ? "-" : "+") << endl;
#endif
                            
                            stack.back().todo.push_back(inbound);
                        } else if (cactus.find(inbound) == cactus.find(graph->flip(inbound)) && !graph->get_is_reverse(inbound)) {
                            // Self loop.
                            // We only want them forward so each becomes only one empty chain.
                            
#ifdef debug
                            cerr << "\t\t\tContain edge " << graph->get_id(inbound) << (graph->get_is_reverse(inbound) ? "-" : "+") << endl;
#endif
                        
                            begin_chain(inbound);
                            end_chain(inbound);
                            
                            visited.insert(inbound);
                        } 
                    });
                }
            }
            
            longest_paths.pop_back();
        } else {
            // We will root on a cycle for its component, if not already covered.
            
            if (!visited.count(longest_cycles.back().second)) {
                // This connected component hasn't been done yet.
                
#ifdef debug
                cerr << "Rooting component at cycle for " << graph->get_id(longest_cycles.back().second) << endl;
#endif
            
                // We have an edge on the longest cycle. But it may be reading into and out of nodes that also contains other cycles, bridge edges, and so on.
                // If we declare this longest cycle to be a chain, we need to make sure that both those nodes become snarls in a chain.
                // So we introduce a chain that starts and ends with this edge.
                // We can't quite articulate that as a todo list entry, so we forge two stack frames.
            
                // Stack up a root/null snarl containing this cycle as a chain.
                stack.emplace_back();
                stack.back().is_snarl = true;
                
                // Stack up a frame for doing the chain, with the cycle-closing edge as both ends.
                stack.emplace_back();
                stack.back().is_snarl = false;
                stack.back().bounds = make_pair(longest_cycles.back().second, longest_cycles.back().second);
                
                // We'll find all the self edges OK when we look in the first/last snarls on the chain.
            }
            
            longest_cycles.pop_back();
        }
        
        while (!stack.empty()) {
            auto& frame = stack.back();
            
#ifdef debug
            cerr << "At stack frame " << stack.size() - 1 << " for ";
            if (stack.size() == 1) {
                cerr << "root";
            } else {
                cerr << (frame.is_snarl ? "snarl" : "chain") << " " << graph->get_id(frame.bounds.first) << (graph->get_is_reverse(frame.bounds.first) ? "-" : "+")
                    << " to " << graph->get_id(frame.bounds.second) << (graph->get_is_reverse(frame.bounds.second) ? "-" : "+");
            }
            cerr << endl;
#endif
            
            if (stack.size() > 1 && !frame.saw_children) {
                // We need to queue up the children; this is the first time we are doing this frame.
                frame.saw_children = true;
                
#ifdef debug
                cerr << "\tAnnouncing entry..." << endl;
#endif
                
                // Announce entering this snarl or chain in the traversal
                (frame.is_snarl ? begin_snarl : begin_chain)(frame.bounds.first);
                
#ifdef debug
                cerr << "\tLooking for children..." << endl;
#endif
                
                if (frame.is_snarl) {
                    
                    // Visit the start and end of the snarl, for decomposition purposes.
                    visited.insert(frame.bounds.first);
                    visited.insert(frame.bounds.second);
                    // TODO: register as part of snarl in index
                    
                    // Make sure this isn't trying to be a unary snarl
                    assert(frame.bounds.first != frame.bounds.second);
                    
                    // For a snarl, we need to find all the bridge edges and all the incoming cycle edges
                    cactus.for_each_member(cactus.find(frame.bounds.first), [&](handle_t inbound) {
                        
                        if (inbound == frame.bounds.first || graph->flip(inbound) == frame.bounds.second) {
                            // This is our boundary; don't follow it as contents.
#ifdef debug
                            cerr << "\t\tStay inside snarl-bounding edge " << graph->get_id(inbound) << (graph->get_is_reverse(inbound) ? "-" : "+") << endl;
#endif
                        } else if (forest.find(graph->flip(inbound)) != forest.find(inbound)) {
                            // This is a bridge edge. The other side is a different component in the bridge graph.
                            
#ifdef debug
                            cerr << "\t\tLook at bridge edge " << graph->get_id(inbound) << (graph->get_is_reverse(inbound) ? "-" : "+") << endl;
#endif
                            
                            frame.todo.push_back(inbound);
                        } else if (next_along_cycle.count(inbound)) {
                            // This edge is the incoming edge for a cycle. Queue it up.
                            
#ifdef debug
                            cerr << "\t\tLook at cycle edge " << graph->get_id(inbound) << (graph->get_is_reverse(inbound) ? "-" : "+") << endl;
#endif
                            frame.todo.push_back(inbound);
                        } else if (cactus.find(graph->flip(inbound)) == cactus.find(inbound) && !graph->get_is_reverse(inbound)) {
                            // Count all self edges as empty chains, but only in one orientation.
                            
#ifdef debug
                            cerr << "\t\tContain edge " << graph->get_id(inbound) << (graph->get_is_reverse(inbound) ? "-" : "+") << endl;
#endif
                            
                            begin_chain(inbound);
                            end_chain(inbound);
                            
                            visited.insert(inbound);
                        }
                    });
                } else {
                    // For a chain, we need to queue up all the edges reading into child snarls, paired with the edges reading out of them.
                    // We know we're a cycle that can be followed.
                    handle_t here = frame.bounds.first;
                    unordered_set<handle_t> seen;
                    size_t region_start = frame.todo.size();
                    do {
                    
#ifdef debug
                        cerr << "\t\tLook at cycle edge " << graph->get_id(here) << (graph->get_is_reverse(here) ? "-" : "+") << endl;
#endif
                    
                        // We shouldn't loop around unless we hit the end of the chain.
                        assert(!seen.count(here));
                        seen.insert(here);
                    
                        // Queue up
                        frame.todo.push_back(here);
                        here = next_along_cycle.at(here);
                        // TODO: when processing entries, we're going to look them up in next_along_cycle again.
                        // Can we dispense with the todo list and create stack frames directly?
                        
                        // Keep going until we come to the end.
                        // We do this as a do-while because the start may be the end but we still want to go around the cycle.
                    } while (here != frame.bounds.second);
                    
                    // Now we have put all the snarls in the chain on the to
                    // do list. But we process the to do list from the end, so
                    // as is we're going to traverse them backward along the
                    // chain. We want to see them forward along the chain
                    // instead, so reverse this part of the vector.
                    // TODO: should we make the to do list a list? That would
                    // save a reverse but require a bunch of allocations and
                    // pointer follows.
                    std::reverse(frame.todo.begin() + region_start, frame.todo.end());
                }
                
            }
            
            if (!frame.todo.empty()) {
                // Until we run out of edges to work on
                handle_t task = frame.todo.back();
                frame.todo.pop_back();
                
                if (frame.is_snarl) {
                    // May have a bridge edge or a cycle edge, both inbound.
                    auto next_along_cycle_it = next_along_cycle.find(task);
                    if (next_along_cycle_it != next_along_cycle.end()) {
                        // To handle a cycle in the current snarl
                        
#ifdef debug
                        cerr << "\tHandle cycle edge " << graph->get_id(task) << (graph->get_is_reverse(task) ? "-" : "+") << endl;
#endif
                        
                        // We have the incoming edge, so find the outgoing edge along the same cycle
                        handle_t outgoing = next_along_cycle_it->second;
                        
#ifdef debug
                        cerr << "\t\tEnds chain starting at " << graph->get_id(outgoing) << (graph->get_is_reverse(outgoing) ? "-" : "+") << endl;
#endif

#ifdef debug
                        cerr << "\t\t\tRecurse on chain " << graph->get_id(outgoing) << (graph->get_is_reverse(outgoing) ? "-" : "+") << " to "
                            << graph->get_id(task) << (graph->get_is_reverse(task) ? "-" : "+") << endl;
#endif
                       
                        if (stack.size() > 1) {
                            // We have boundaries. Make sure we don't try and
                            // do a chain that starts or ends with our
                            // boundaries. That's impossible.
                            assert(frame.bounds.first != outgoing);
                            assert(frame.bounds.second != task);
                        }
                        
                        // Recurse on the chain bounded by those edges, as a child
                        stack.emplace_back();
                        stack.back().is_snarl = false;
                        stack.back().bounds = make_pair(outgoing, task);
                        
                    } else {
                        // To handle a bridge edge in the current snarl:
                        
#ifdef debug
                        cerr << "\tHandle bridge edge " << graph->get_id(task) << (graph->get_is_reverse(task) ? "-" : "+") << endl;
#endif
                        
                        // Flip it to look out
                        handle_t edge = graph->flip(task);
#ifdef debug
                        cerr << "\t\tWalk edge " << graph->get_id(edge) << (graph->get_is_reverse(edge) ? "-" : "+") << endl;
#endif
                        // Track the head in the Cactus graph for the bridge edges we walk.
                        handle_t cactus_head = cactus.find(edge);
                        // And track where its bridge forest component points to as towards the deepest leaf.
                        auto deepest_it = towards_deepest_leaf.find(forest.find(cactus_head));
                        while (deepest_it != towards_deepest_leaf.end()) {
                            // Follow its path down bridge graph heads, to the
                            // deepest bridge graph leaf head (which has no
                            // deeper child)
                            
                            // See what our next bridge edge comes out of in the Cactus graph
                            handle_t next_back_head = cactus.find(graph->flip(deepest_it->second));
                            
#ifdef debug
                            cerr << "\t\t\tHead: " << graph->get_id(cactus_head) << (graph->get_is_reverse(cactus_head) ? "-" : "+") << endl;
                            cerr << "\t\t\tNext edge back head: " << graph->get_id(next_back_head) << (graph->get_is_reverse(next_back_head) ? "-" : "+") << endl;
#endif
                            
                            if (cactus_head != next_back_head) {
                                // We skipped over a run of interlinked cycle in the bridge tree.
                                
                                // We need to find a path of cycles to complete the path in the bridge tree.
                                
                                // Each cycle needs to be cut into two pieces
                                // that can be alternatives in the snarl.
                                
#ifdef debug
                                cerr << "\t\t\tFind skipped cycle path" << endl;
#endif
                                
                                vector<handle_t> cycle_path = cactus.find_cycle_path_in_cactus(next_along_cycle, cactus_head, next_back_head);
                                
                                while (!cycle_path.empty()) {
                                    // Now pop stuff off the end of the path and
                                    // merge it with the component next_back_head
                                    // is in, making sure to pinch off the cycles
                                    // we cut as we do it.
                                    
                                    // Walk the cycle (again) to find where it hits the end component.
                                    // TODO: Save the first traversal we did!
                                    auto through_path_member = next_along_cycle.find(cycle_path.back());
                                    auto through_end = through_path_member;
                                    do {
                                        // Follow the cycle until we reach the edge going into the end component.
                                        through_end = next_along_cycle.find(through_end->second);
                                    } while (cactus.find(through_end->first) != cactus.find(next_back_head));
                                    
                                    // Now pinch the cycle
                                    
#ifdef debug
                                    cerr << "\t\t\tPinch cycle between " << graph->get_id(cycle_path.back()) << (graph->get_is_reverse(cycle_path.back()) ? "-" : "+")
                                        << " and " << graph->get_id(through_end->first) << (graph->get_is_reverse(through_end->first) ? "-" : "+") << endl;
#endif
                                    
                                    // Merge the two components where the bridge edges attach, to close the two new cycles.
                                    cactus.merge(cycle_path.back(), next_back_head);
                                    
#ifdef debug
                                    cerr << "\t\t\t\tExchange successors of " << graph->get_id(through_path_member->first) << (graph->get_is_reverse(through_path_member->first) ? "-" : "+")
                                        << " and " << graph->get_id(through_end->first) << (graph->get_is_reverse(through_end->first) ? "-" : "+") << endl;
#endif
                                    
                                    // Exchange their destinations to pinch the cycle in two.
                                    std::swap(through_path_member->second, through_end->second);
                                    
                                    if (through_path_member->first == through_path_member->second) {
                                        // Now a self loop cycle. Delete the cycle.
                                        
#ifdef debug
                                        cerr << "\t\t\t\t\tDelete self loop cycle " << graph->get_id(through_path_member->first) << (graph->get_is_reverse(through_path_member->first) ? "-" : "+") << endl;
#endif
                                        
                                        // Won't affect other iterators.
                                        next_along_cycle.erase(through_path_member);
                                    }
                                    
                                    if (through_end->first == through_end->second) {
                                        // Now a self loop cycle. Delete the cycle.
                                        
#ifdef debug
                                        cerr << "\t\t\t\t\tDelete self loop cycle " << graph->get_id(through_end->first) << (graph->get_is_reverse(through_end->first) ? "-" : "+") << endl;
#endif
                                        
                                        // Won't affect other iterators.
                                        next_along_cycle.erase(through_end);
                                    }
                                    
                                    // And pop it off and merge the end (which now includes it) with whatever came before it on the path.
                                    cycle_path.pop_back();
                                }
                            }
                            
                            // Record the new cycle we are making from this bridge path
                            next_along_cycle[edge] = deepest_it->second;
                            
                            // Advance along the bridge tree path.
                            edge = deepest_it->second;
#ifdef debug
                            cerr << "\t\tWalk edge " << graph->get_id(edge) << (graph->get_is_reverse(edge) ? "-" : "+") << endl;
#endif
                            cactus_head = cactus.find(edge);
                            deepest_it = towards_deepest_leaf.find(forest.find(cactus_head));
                        }
                        
                        // When you get to the end
                        
                        if (edge == graph->flip(task)) {
                            // It turns out there's only one edge here.
                            // It is going to become a contained self-loop, instead of a real cycle
                            
                            // Record we visited it.
                            visited.insert(edge);
                            
#ifdef debug
                            cerr << "\t\tContain new self-loop " << graph->get_id(edge) << (graph->get_is_reverse(edge) ? "-" : "+") << endl;
#endif

                            // Register as part of snarl in index
                            begin_chain(graph->forward(edge));
                            end_chain(graph->forward(edge));
                        } else {
                            // Close the cycle we are making out of the bridge
                            // forest path.
                            // The last edge crossed currently reads into the end
                            // component, but will read into us after the merge.
                            // The cycle comes in through there and leaves backward
                            // through the inbound bridge edge we started with.
                            next_along_cycle[edge] = graph->flip(task);
                            
#ifdef debug
                            cerr << "\t\tClose cycle between " << graph->get_id(edge) << (graph->get_is_reverse(edge) ? "-" : "+")
                                << " and " << graph->get_id(task) << (graph->get_is_reverse(task) ? "-" : "+") << endl;
#endif
                            
                        }
                        
                        // Merge the far end of the last bridge edge (which may have cycles on it) into the current snarl
                        
                        // First find all the new cycles this brings along.
                        // It can't bring any bridge edges.
                        // This will detect the cycle we just created.
                        cactus.for_each_member(cactus_head, [&](handle_t inbound) {
                            // TODO: deduplicate with snarl setup
                            if (next_along_cycle.count(inbound)) {
                            
#ifdef debug
                                cerr << "\t\tInherit cycle edge " << graph->get_id(inbound) << (graph->get_is_reverse(inbound) ? "-" : "+") << endl;
#endif
                            
                                // This edge is the incoming edge for a cycle. Queue it up.
                                frame.todo.push_back(inbound);
                            } else if (cactus.find(graph->flip(inbound)) == cactus.find(inbound) && !graph->get_is_reverse(inbound)) {
                            
#ifdef debug
                                cerr << "\t\tInherit contained edge " << graph->get_id(inbound) << (graph->get_is_reverse(inbound) ? "-" : "+") << endl;
#endif
                            
                                // Count all self edges as empty chains, but only from one side.
                                begin_chain(inbound);
                                end_chain(inbound);
                                
                                visited.insert(inbound);
                            }   
                        });
                        
                        // Then do the actual merge.
                        cactus.merge(edge, task);
                    
                        // Now we've queued up the cycle we just made out of
                        // the bridge edges, along with any cycles we picked up
                        // from the end of the bridge tree path.
                    }
                } else {
                
#ifdef debug
                    cerr << "\tHandle cycle edge " << graph->get_id(task) << (graph->get_is_reverse(task) ? "-" : "+") << endl;
#endif
                
                    // We're a chain, and WLOG a chain that represents a cycle.
                    // We have an edge.
                    // We need to find the other edge that defines the snarl, and recurse into the snarl.
                    handle_t out_edge = next_along_cycle.at(task);
                    
#ifdef debug
                    cerr << "\t\tRecurse on snarl " << graph->get_id(task) << (graph->get_is_reverse(task) ? "-" : "+") << " to "
                        << graph->get_id(out_edge) << (graph->get_is_reverse(out_edge) ? "-" : "+")<< endl;
#endif
                    
                    stack.emplace_back();
                    stack.back().is_snarl = true;
                    stack.back().bounds = make_pair(task, out_edge);
                }
            
            } else {
                // Now we have finished a stack frame!
                
                if (stack.size() > 1) {
                    // We have bounds
                    
#ifdef debug
                    cerr << "\tAnnouncing exit..." << endl;
#endif
                
                    // Announce leaving this snarl or chain in the traversal
                    (frame.is_snarl ? end_snarl : end_chain)(frame.bounds.second);
                
                }
                
#ifdef debug
                cerr << "\tReturn to parent frame" << endl;
#endif
                
                stack.pop_back();
            }
            
            
        }
    }
    
}

SnarlManager IntegratedSnarlFinder::find_snarls_parallel() {

    vector<unordered_set<id_t>> weak_components = handlealgs::weakly_connected_components(graph);
    vector<SnarlManager> snarl_managers(weak_components.size());

    #pragma omp parallel for schedule(dynamic, 1)
    for (size_t i = 0; i < weak_components.size(); ++i) {
        const HandleGraph* subgraph;
        if (weak_components.size() == 1) {
            subgraph = graph;
        } else {
            // turn the component into a graph
            subgraph = new SubgraphOverlay(graph, &weak_components[i]);
        }
        IntegratedSnarlFinder finder(*subgraph);
        // find the snarls without building the index
        snarl_managers[i] = finder.find_snarls_unindexed();
        if (weak_components.size() != 1) {
            // delete our component graph overlay
            delete subgraph;
        }
    }

    // merge the managers into the biggest one.
    size_t biggest_snarl_idx = 0;
    for (size_t i = 1; i < snarl_managers.size(); ++i) {
        if (snarl_managers[i].num_snarls() > snarl_managers[biggest_snarl_idx].num_snarls()) {
            biggest_snarl_idx = i;
        }
    }
    for (size_t i = 0; i < snarl_managers.size(); ++i) {
        if (i != biggest_snarl_idx) {
            snarl_managers[i].for_each_snarl_unindexed([&](const Snarl* snarl) {
                snarl_managers[biggest_snarl_idx].add_snarl(*snarl);
            });
        }
    }
    snarl_managers[biggest_snarl_idx].finish();
    return std::move(snarl_managers[biggest_snarl_idx]);
}

}
