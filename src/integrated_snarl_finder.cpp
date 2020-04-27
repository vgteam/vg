///
///  \file integrated_snarl_finder.cpp
///
///

#include "integrated_snarl_finder.hpp"

#include "algorithms/three_edge_connected_components.hpp"

#include <bdsg/overlays/overlay_helper.hpp>
#include <structures/union_find.hpp>

#include <array>

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
    /// If a path in the forest doesn't beat the length of the cycle that lives
    /// in its tree, it is omitted.
    pair<vector<pair<size_t, vector<handle_t>>>, unordered_map<handle_t, handle_t>> longest_paths_in_forest(const vector<pair<size_t, handle_t>>& longest_simple_cycles) const;
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
                auto frame_it = visited_frame.find(frame_head);
                if (frame_it == visited_frame.end()) {
                    // First visit to here.
                    
                    // Mark visited at this stack level
                    frame_it = visited_frame.emplace_hint(frame_it, frame_head, stack.size() - 1);
                    
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
                        if (frame_it->second > connected_it->second) {
                            // We have an edge to something that was visited above
                            // our stack level. It can't be a self loop, and it
                            // must close a unique cycle.
                        
                            // Get the stack frame we connected to (which may
                            // be the same as our frame, for self loops)
                            auto& connected_frame = stack[connected_it->second];
                            
                            // Walk and measure the cycle
                            size_t cycle_length_bp = graph->get_length(edge_into);
                            handle_t prev_edge = edge_into;
                            for (size_t i = connected_it->second; i < stack.size(); i++) {
                                // For each edge along the cycle...
                                // Measure it
                                cycle_length_bp += graph->get_length(stack[i].here);
                                // Record the cycle membership
                                next_edge[prev_edge] = stack[i].here;
                                // Advance
                                prev_edge = stack[i].here;
                            }
                            // Close the cycle
                            next_edge[prev_edge] = edge_into;
                            
                            if (cycle_length_bp > longest_cycle.first) {
                                // New longest cycle (or maybe only longest cycle).
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
   
    return to_return;
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
            
            while (!stack.empty()) {
                // Until the DFS is done
                auto& frame = stack.back();
                // Find the node that following this edge got us to.
                auto frame_head = find(frame.here);
                
                auto frame_it = records.find(frame_head);
                if (frame_it == records.end()) {
                    // First visit to here.
                    
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
                    
                    if (!records.count(connected_head)) {
                        // Forward edge. Recurse.
                        stack.emplace_back();
                        stack.back().here = edge_into;
                    }
                } else {
                    // No children left.
                    
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
                        
                        if (deepest_child_edge_it != deepest_child_edge.end()) {
                            // And if we have a child to go on with, we add the length of that path
                            record.leaf_path_length += records[find(deepest_child_edge_it->second)].leaf_path_length;
                        }
                        
                        // Fill in deepest_child_edge for the parent if not filled in already, or if we beat what's there.
                        // Also maintain parent's second_deepest_child_edge.
                        auto parent_deepest_child_it = deepest_child_edge.find(parent_head);
                        if (parent_deepest_child_it == deepest_child_edge.end()) {
                            // Emplace in the map where we didn't find anything.
                            deepest_child_edge.emplace_hint(parent_deepest_child_it, parent_head, frame.here);
                        } else if(records[find(parent_deepest_child_it->second)].leaf_path_length < record.leaf_path_length) {
                            // We are longer than what's there now
                            
                            // Demote what's there to second-best
                            parent_record.second_deepest_child_edge = parent_deepest_child_it->second;
                            parent_record.has_second_deepest_child = true;
                            
                            // Replace the value we found
                            parent_deepest_child_it->second = frame.here;
                        } else if (parent_record.has_second_deepest_child && records[find(parent_record.second_deepest_child_edge)].leaf_path_length < record.leaf_path_length) {
                            // We are a new second deepest child.
                            parent_record.second_deepest_child_edge = frame.here;
                        }
                    }
                    
                    // The length of the longest leaf-leaf path converging at or under any child (if any) is in record.longest_subtree_path_length.
                    
                    if (record.has_second_deepest_child) {
                        // If there's a second incoming leaf path there's a converging leaf-leaf path here.
                        // Grab the length of the longest leaf-leaf path converging exactly here.
                        // TODO: can we not look up he deepest child's record again?
                        size_t longest_here_path_length = records[find(deepest_child_edge_it->second)].leaf_path_length;
                        longest_here_path_length += records[find(record.second_deepest_child_edge)].leaf_path_length;
                        
                        if (record.longest_subtree_path_root == frame_head || longest_here_path_length > record.longest_subtree_path_length) {
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
                            
                            // No child has contributed their leaf-leaf path so far, or ours is better.
                            parent_record.longest_subtree_path_root = frame_head;
                            parent_record.longest_subtree_path_length = record.longest_subtree_path_length;
                        }
                    }
                   
                    if (stack.size() == 1) {
                        // When we get back to the root
                        
                        if (record.longest_subtree_path_length >= root_cycle_length) {
                            // Either we didn't root at a cycle, or we found a longer leaf-leaf path that should be the decomposition root instead.
                            
                            // We need to record the longest tree path.
                            longest_tree_paths.emplace_back();
                            longest_tree_paths.back().first = record.longest_subtree_path_length;
                            auto& path = longest_tree_paths.back().second;
                            
                            auto& path_root_frame = records[record.longest_subtree_path_root];
                            assert(path_root_frame.has_second_deepest_child);
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
                            // Now trace the actual longest path from root to leaf and add it on
                            path.push_back(deepest_child_edge[record.longest_subtree_path_root]);
                            path_trace_it = deepest_child_edge.find(find(path.back()));
                            while (path_trace_it != deepest_child_edge.end()) {
                                // Follow the deepest child relationships until they run out.
                                path.push_back(path_trace_it->second);
                                path_trace_it = deepest_child_edge.find(find(path.back()));
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
                                } else if (parent_record.has_second_deepest_child && records[find(parent_record.second_deepest_child_edge)].leaf_path_length < child_record.leaf_path_length) {
                                    // We are a new second deepest child.
                                    parent_record.second_deepest_child_edge = parent_child_edge;
                                }
                               
                                // Now the new child, if its path is deep enough, is the parent's new deepest or second deepest child edge.
                                // Go up a level, disqualify the grandparent, and see who wins.
                                // The new root has no parent itself, so all its edges are eligible and the one with the longest path wins.
                                convergence_to_old_root.pop_back();
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
    
    return to_return;
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
    
    // Get cycle information: longest cycle in each connected component, and next edge along cycle for each edge (in one orientation)
    pair<vector<pair<size_t, handle_t>>, unordered_map<handle_t, handle_t>> cycles = cactus.cycles_in_cactus();
    auto& longest_cycles = cycles.first;
    auto& next_along_cycle = cycles.second;
    
    for (auto& kv : next_along_cycle) {
        // Merge along all cycles in the bridge forest
        forest.merge(kv.first, kv.second);
    }
    
    // Now we find the longest path in each tree in the bridge forest, with its
    // length in bases.
    //
    // We also find, for each bridge edge component head, the edge towards the
    // deepest bridge edge tree leaf, which lets us figure out how to root
    // dangly bits into chains.
    //
    // Make sure to root at the nodes corresponding to the collapsed longest
    // cycles, if the leaf-leaf paths don't win their components.
    pair<vector<pair<size_t, vector<handle_t>>>, unordered_map<handle_t, handle_t>> forest_paths = forest.longest_paths_in_forest(longest_cycles);
    auto& longest_paths = forest_paths.first;
    auto& towards_deepest_leaf = forest_paths.second;
    
    // Make sure we are looking at all the cycles and leaf-leaf paths in order.
    // We need basically priority queue between them. But we don't need insert so we jsut sort.
    std::sort(longest_cycles.begin(), longest_cycles.end());
    std::sort(longest_paths.begin(), longest_paths.end());
  
    // Now, keep a set of all the edges that have found a place in the decomposition.
    unordered_set<handle_t> visited;
    
    while(visited.size() < graph->get_edge_count()) {
        // While we haven't touched everything
        
        // We have a stack
        
        // To handle a bridge edge in the current snarl:
            
            // Follow its path to the deepest bridge graph leaf
            
            // Along that path, merge inner ednpoints of successive bridge edges and pinch apart any cycles that shortens
            
            // When you get to the end, merge the far end of the last bridge edge into the current snarl
            
            // Write a cycle in the current snarl following all the bridge edges (and keep any cycles we got from the bridge graph tip node)
        
        // To handle a cycle in the current snarl:
        
            // Make a chain in the current snarl
        
            // For each pair of successive edges, recurse on that as a snarl in the chain
        
        if (longest_cycles.empty() || (!longest_paths.empty() && longest_cycles.back().first < longest_paths.back().first)) {
            assert(!longest_paths.empty());
            
            // We will root on a tip-tip path for its component, if not already covered.
            
            // Rewrite the deepest bridge graph leaf path map to point from one end of the tip-tip path to the other
                // TODO: bump this down into the bridge path finding function
            
            // Stack up a root/null snarl containing this bridge edge.
            
            longest_paths.pop_back();
        } else {
            // We will root on a cycle for its component, if not already covered.
            
            // Stack up a root/null snarl containing this cycle
            
            longest_cycles.pop_back();
        }
    }
    
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
