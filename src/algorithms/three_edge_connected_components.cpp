#include "three_edge_connected_components.hpp"

extern "C" {
#include "sonLib/sonLibList.h"
#include "sonLib/sonLibTuples.h"
#include "sonLib/3_Absorb3edge2x.h"
}

#include <structures/union_find.hpp>

#include <limits>
#include <cassert>
#include <iostream>
#include <sstream>

//#define debug

namespace vg {
namespace algorithms {

using namespace std;

void three_edge_connected_component_merges_dense(size_t node_count, size_t first_root, 
    const function<void(size_t, const function<void(size_t)>&)>& for_each_connected_node,
    const function<void(size_t, size_t)>& same_component) {
    
    // Independent implementation of Norouzi and Tsin (2014) "A simple 3-edge
    // connected component algorithm revisited", which can't really be
    // understood without Tsin (2007) "A Simple 3-Edge-Connected Component
    // Algorithm".
    
    // That algorithm assumes that all bridge edges are removed (i.e.
    // everything is at least 2-connected), but we hack it a bit to generalize
    // to graphs with bridge edges. It also assumes there are no self loops,
    // but this implementation detects and allows self loops.
    
    // The algorithm does a depth-first search through the graph, and is based
    // on this "absorb-eject" operation. You do it at a node, across ("on") an
    // edge. It (conceptually) steals all the edges from the node at the other
    // end of the edge, deletes the edge, and deletes the other node as well if
    // it has a degree greater than 2. (The original algorithm didn't have to
    // deal with degree 1; here we treat it about the same as degree 2 and
    // leave the node floating in its own 3 edge connected component, while
    // hiding the single edge from the real logic of the algorithm.)
    
    // Because of guarantees about the order in which we traverse the graph, we
    // don't actually have to *do* any of the absorb-eject graph topology
    // modifications. Instead, we just have to keep track of updates to nodes'
    // "effective degree" in what would be the modified graph, and allow
    // certain paths that we track during the algorithm to traverse the stolen
    // edges.
    
    // For each node, we keep track of a path. Because of guarantees we get
    // from the algorithm, we know that paths can safely share tails. So to
    // represent the tail of a path, we can just point to another node (or
    // nowhere if the path ends). The head of a path is tougher, because a
    // node's path can be empty. The node may not be on its own path. It is not
    // immediately clear from analyzing the algorithm whether a node can have a
    // nonempty path that it itself is not on, or be the tail of another node's
    // path without being on that path. To support those cases, we also give
    // each node a flag for whether it is on its own path.
    
    // TODO: should we template this on an integer size so we can fit more work
    // in less memory bandwidth when possible?
    using number_t = size_t;
    assert(node_count < numeric_limits<number_t>::max());
   
    /// This defines the data we track for each node in the graph
    struct TsinNode {
        /// When in the DFS were we first visited?
        number_t dfs_counter;
        /// When in the DFS were we last visited?
        /// Needed for finding replacement neighbors to implement path range
        /// absorption in part 1.3, when we're asked for a range to a neighbor
        /// that got eaten.
        number_t dfs_exit;
        /// What is our "low point" in the search. This is the earliest
        /// dfs_counter for a node that this node or any node in its DFS
        /// subtree has a back-edge to.
        number_t low_point;
        /// What is the effective degree of this node in the graph with all the
        /// absorb-eject modifications applied?
        number_t effective_degree = 0;
        /// What node has the continuation of this node's path? If equal to
        /// numeric_limits<number_t>::max(), the path ends after here.
        /// The node's path is the path from this node, into its DFS subtree,
        /// to (one of) the nodes in the subtree that has the back-edge that
        /// caused this node's low point to be so low. Basically a low point
        /// traceback.
        number_t path_tail;
        /// Is this node actually on its own path?
        /// Nodes can be removed from their paths if those nodes don't matter
        /// any more (i.e. got absorbed) but their paths still need to be tails
        /// for other paths.
        bool is_on_path;
        /// Has the node been visited yet? Must be 0. TODO: Move to its own
        /// vector to make zeroing them all free-ish with page table
        /// shenanigans.
        bool visited = false;
    };
    
    // We need to have all the nodes pre-allocated, so node references don't
    // invalidate when we follow edges.
    vector<TsinNode> nodes(node_count);
    
    // We need to say how to absorb-eject along a whole path.
    //
    // We let you specify the node to absorb into; if it isn't
    // numeric_limits<number_t>::max(), it is assumed to be the first node, and
    // actually on the path, and path_start (if itself on its path) is also
    // absorbed into it. This lets you absorb into a path with something
    // prepended, without constructing the path.
    //
    // Similarly, we let you specify a past end to stop before. If this isn't
    // numeric_limits<number_t>::max(), we stop and don't absorb the specified
    // node, if we reach it. This lets us implement absorbing a range of a
    // path, as called for in the algorithm.
    //
    // If you specify a past_end, and we never reach it, but also don't have
    // just a single-node, no-edge "null" path, then something has gone wrong
    // and we've violated known truths about the algorithm.
    auto absorb_all_along_path = [&](number_t into, number_t path_start, number_t path_past_end) {
        
#ifdef debug
        cerr << "(Absorbing all along path into " << into << " from " << path_start << " to before " << path_past_end << ")" << endl;
#endif
    
        // Set this to false as soon as we cross an edge
        bool path_null = true;
    
        number_t here = path_start;
        while (here != path_past_end) {
            // Until we hit the end of the path
            
#ifdef debug
            cerr << "(\tAt " << here  << ")" << endl;
#endif
            
            if (here == numeric_limits<number_t>::max()) {
                // We hit the end of the path and never saw path_past_end.
                
#ifdef debug
                cerr << "(\t\tReached path end and missed waypoint)" << endl;
#endif
                // Only allowed if the path was actually edge-free and no merges needed to happen.
                assert(path_null);
                
#ifdef debug
                cerr << "(\t\t\tBut path was empty of edges)" << endl;
#endif
                // Stop now.
                break;
            }
            
            // Find the node we are at
            auto& here_node = nodes[here];
            
            if (here_node.is_on_path) {
                // We're actually on the path.
                
#ifdef debug
                cerr << "(\t\tOn path)" << endl;
#endif

                if (into == numeric_limits<number_t>::max()) {
                    // We haven't found a first node to merge into yet; it is
                    // this one.
                    
#ifdef debug
                    cerr << "(\t\tUse as into)" << endl;
#endif

                    into = here;
                } else {
                    // We already have a first node to merge into, so merge.
                    
#ifdef debug
                    cerr << "(\t\tMerge with " << into << ")" << endl;
#endif

                    // We are doing a merge! We'd better actually find the
                    // ending range bound, or something is wrong with our
                    // implementation of the algorithm.
                    path_null = false;
                    
                    // Update the effective degrees as if we merged this node
                    // with the connected into node.
                    nodes[into].effective_degree = (nodes[into].effective_degree +
                                                    here_node.effective_degree - 2);

                    // Merge us into the same 3 edge connected component
                    same_component(into, here);
                }
            }
            
            // Advance to the tail of the path
            here = here_node.path_tail;
            
#ifdef debug
            cerr << "(\t\tNext: " << here << ")" << endl;
#endif

        }

#ifdef debug
        cerr << "(Done absorbing)" << endl;
#endif

    };
    
    // For debugging, we need to be able to dump a node's stored path
    auto path_to_string = [&](number_t node) {
        stringstream s;
        
        number_t here = node;
        bool first = true;
        while (here != numeric_limits<number_t>::max()) {
            if (nodes[here].is_on_path) {
                if (first && nodes[here].path_tail == numeric_limits<number_t>::max()) {
                    // Just a single node, no edge
                    s << "(just " << here << ")";
                    break;    
                }
                
                if (first) {
                    first = false;
                } else {
                    s << "-";
                }
                s << here;
            }
            here = nodes[here].path_tail;
        }
        
        return s.str();
    };
    
    // We need a DFS stack that we manage ourselves, to avoid stack-overflowing
    // as we e.g. walk along big cycles.
    struct DFSStackFrame {
        /// Track the node that this stack frame represents
        number_t current;
        /// Track all the neighbors left to visit.
        /// When we visit a neighbor we pop it off the back.
        vector<number_t> neighbors;
        /// When we look at the neighbors, we need to be able to tell the tree
        /// edge to the parent from further back edges to the parent. So we
        /// have a flag for whether we have seen the parent tree edge already,
        /// and the first neighbors entry that is our parent will get called
        /// the tree edge.
        bool saw_parent_tree_edge = false; 
        /// Track whether we made a recursive DFS call into the last neighbor
        /// or not. If we did, we need to do some work when we come out of it
        /// and return to this frame.
        bool recursing = false;
    };
    
    vector<DFSStackFrame> stack;
    
    // We need a way to produce unvisited nodes when we run out of nodes in a
    // connected component. This will always point to the next unvisited node
    // in order. If it points to node_count, all nodes are visited. When we
    // fisit this node, we have to scan ahead for the next unvisited node, in
    // number order.
    number_t next_unvisited = 0;
    
    // We also keep a global DFS counter, so we don't have to track parent
    // relationships when filling it in on the nodes.
    //
    // The paper starts it at 1, so we do too.
    number_t dfs_counter = 1;
    
    while (next_unvisited != node_count) {
        // We haven't visited everything yet.
        if (!nodes[first_root].visited) {
            // If possible start at the suggested root
            stack.emplace_back();
            stack.back().current = first_root;
        } else {
            // Stack up the next unvisited node.
            stack.emplace_back();
            stack.back().current = next_unvisited;
        }
        
#ifdef debug
        cerr << "Root a search at " << stack.back().current << endl;
#endif
        
        while (!stack.empty()) {
            // While there's still nodes on the DFS stack from the last component we broke into
            // Grab the stack frame.
            // Note that this reference will be invalidated if we add stuff to the stack!
            auto& frame = stack.back();
            // And the current node
            auto& node = nodes[frame.current];
            
            if (!node.visited) {
                // This is the first time we are in this stack frame. We need
                // to do the initial visit of the node and set up the frame
                // with the list of edges to do.
                node.visited = true;
                
#ifdef debug
                cerr << "First visit of node " << frame.current << endl;
#endif
                
                if (frame.current == next_unvisited) {
                    // We need to find the next unvisited node, if any, since
                    // we just visited what it used to be.
                    do {
                        next_unvisited++;
                    } while (next_unvisited != node_count && nodes[next_unvisited].visited);
                }
                
                node.dfs_counter = dfs_counter;
                dfs_counter++;
                node.low_point = node.dfs_counter;
                // Make sure the node's path is just itself
                node.path_tail = numeric_limits<number_t>::max();
                node.is_on_path = true;
                
#ifdef debug
                cerr << "\tDFS: " << node.dfs_counter
                    << " low point: " << node.low_point
                    << " degree: " << node.effective_degree
                    << " path: " << path_to_string(frame.current) << endl;
#endif
                
                // Stack up all the edges to follow.
                for_each_connected_node(frame.current, [&](size_t connected) {
                    frame.neighbors.push_back(connected);
                });
                
#ifdef debug
                cerr << "\tPut " << frame.neighbors.size() << " edges on to do list" << endl;
#endif
                
                // Now we're in a state where we can process edges.
                // So kick back to the work loop as if we just processed an edge.
                continue;
            } else {
                // We have (possibly 0) edges left to do for this node.
                if (!frame.neighbors.empty()) {
                
#ifdef debug
                    cerr << "Return to node " << frame.current << " with more edges to do" << endl;
                
                    cerr << "\tDFS: " << node.dfs_counter
                        << " low point: " << node.low_point
                        << " degree: " << node.effective_degree
                        << " path: " << path_to_string(frame.current) << endl;
#endif
                
                    // We have an edge to do!
                    // Look up the neighboring node.
                    number_t neighbor_number = frame.neighbors.back();
                    auto& neighbor = nodes[neighbor_number];
                    
                    if (!frame.recursing) {
                        // This is the first time we are thinking about this neighbor.
                        
#ifdef debug
                        cerr << "\tThink of edge to neighbor " << neighbor_number << " for the first time" << endl;
#endif
                    
                        // Increment degree of the node we're coming from
                        node.effective_degree++;
                        
#ifdef debug
                        cerr << "\t\tBump degree to " << node.effective_degree << endl;
#endif
                        
                        if (!neighbor.visited) {
                            // We need to recurse on this neighbor.
                            
#ifdef debug
                            cerr << "\t\tRecurse on unvisited neighbor" << endl;
#endif
                            
                            // So remember we are recursing.
                            frame.recursing = true;
                            // And set up the recursive frame.
                            stack.emplace_back();
                            stack.back().current = neighbor_number;
                            // Kick back to the work loop; we will see the
                            // unvisited node on top of the stack and do its
                            // visit and add its edges to its to do list.
                        } else {
                            // No need to recurse.This is either a back-edge or the back side of the tree edge to the parent.
                            
                            if (stack.size() > 1 && neighbor_number == stack[stack.size() - 2].current && !frame.saw_parent_tree_edge) {
                                // This is the edge we took to get here (tree edge)
#ifdef debug
                                cerr << "\t\tNeighbor is parent; this is the tree edge in." << endl;
#endif

                                // For tree edges, since they aren't either kind of back edge, neither 1.2 nor 1.3 fires.
                                // But the next edge to the parent will be a back edge.
                                frame.saw_parent_tree_edge = true;
                            } else if (neighbor.dfs_counter < node.dfs_counter) {
                                // The edge to the neighbor is an outgoing
                                // back-edge (i.e. the neighbor was visited
                                // first). Paper step 1.2.
                                
#ifdef debug
                                cerr << "\t\tNeighbor is upstream of us (outgoing back edge)." << endl;
#endif
                                
                                if (neighbor.dfs_counter < node.low_point) {
                                    // The neighbor is below our low point.
                                    
#ifdef debug
                                    cerr << "\t\t\tNeighbor has a lower low point ("
                                        << neighbor.dfs_counter << " < " << node.low_point << ")" << endl;
                                    
                                    cerr << "\t\t\t\tAbsorb along path to old low point source" << endl;
#endif

                                    // Absorb along our whole path.
                                    absorb_all_along_path(numeric_limits<number_t>::max(),
                                                          frame.current,
                                                          numeric_limits<number_t>::max());
                                    
                                    // Adopt the neighbor's DFS counter as our
                                    // new, lower low point.
                                    node.low_point = neighbor.dfs_counter;

#ifdef debug
                                    cerr << "\t\t\t\tNew lower low point " << node.low_point << endl;
#endif
                                    
                                    // Our path is now just us.
                                    node.is_on_path = true;
                                    node.path_tail = numeric_limits<number_t>::max();
                                    
#ifdef debug
                                    cerr << "\t\t\t\tNew path " << path_to_string(frame.current) << endl;
#endif

                                } else {
                                
#ifdef debug
                                    cerr << "\t\t\tWe have a sufficiently low low point" << endl;
#endif
                                
                                }
                            } else if (node.dfs_counter < neighbor.dfs_counter) {
                                // The edge to the neighbor is an incoming
                                // back-edge (i.e. we were visited first, but
                                // we recursed into something that got us to
                                // this neighbor already). Paper step 1.3.
                                
#ifdef debug
                                cerr << "\t\tWe are upstream of neighbor (incoming back edge)." << endl;
#endif
                                
                                // Drop our effective degree by 2 (I think
                                // we're closing a cycle or something?)
                                node.effective_degree -= 2;

#ifdef debug
                                cerr << "\t\t\tDrop degree to " << node.effective_degree << endl;
                                
                                cerr << "\t\t\tWant to absorb along path towards low point source through neighbor" << endl;
#endif
                                
                                // Now, the algorithm says to absorb
                                // "P_w[w..u]", a notation that it does not
                                // rigorously define. w is here, and u is the
                                // neighbor. The neighbor is not necessarily
                                // actually *on* our path at this point, not
                                // least of which because the neighbor may have
                                // already been eaten and merged into another
                                // node, which in theory adopted the back edge
                                // we are looking at. In practice we don't have
                                // the data structure to find that node. So
                                // here's the part where we have to do
                                // something clever to "allow certain paths
                                // that we track to traverse the stolen edges".
                                
                                // What we have to do is find the node that
                                // *is* along our path that either is or ate
                                // the neighbor. We don't track the union-find
                                // logic we would need to answer that question,
                                // but both 2007 algorithm implementations I've
                                // seen deal with this by tracking DFS counter
                                // intervals/subtree sizes, and deciding that
                                // the last thin on our path visited no later
                                // than the neighbor, and exited no earlier
                                // than the neighbor (i.e. the last ancestor of
                                // the neighbor on our path) should be our
                                // replacement neighbor.
                                
                                // This makes sense because if the neighbor
                                // merged into anything, it's an ancestor of
                                // the neighbor. So we go looking for it.
                                
                                // TODO: let absorb_all_along_path do this instead?
                                
                                // Start out with ourselves as the replacement neighbor ancestor.
                                number_t replacement_neighbor_number = frame.current;
                                // Consider the next candidate
                                number_t candidate = nodes[replacement_neighbor_number].path_tail;
                                while (candidate != numeric_limits<number_t>::max() &&
                                    nodes[candidate].dfs_counter <= neighbor.dfs_counter &&
                                    nodes[candidate].dfs_exit >= neighbor.dfs_exit) {
                                    
                                    // This candidate is a lower ancestor of the neighbor, so adopt it.
                                    replacement_neighbor_number = candidate;
                                    candidate = nodes[replacement_neighbor_number].path_tail;
                                }
                                
                                auto& replacement_neighbor = nodes[replacement_neighbor_number];
                                
#ifdef debug
                                cerr << "\t\t\tNeighbor currently belongs to node " << replacement_neighbor_number << endl;
                                
                                cerr << "\t\t\tAbsorb along path towards low point source through there" << endl;
#endif
                                
                                // Absorb along our path from ourselves to the
                                // replacement neighbor, inclusive.
                                // Ignores trivial paths.
                                absorb_all_along_path(numeric_limits<number_t>::max(),
                                                      frame.current,
                                                      replacement_neighbor.path_tail);
                                                      
                                // We also have to (or at least can) adopt the
                                // path of the replacement neighbor as our own
                                // path now. That's basically the rest of the
                                // path that we didn't merge.
                                // This isn't mentioned in the paper either,
                                // but I've seen the official implementation do
                                // it, and if we don't do it our path is going
                                // to go through a bunch of stuff we already
                                // merged, and waste time when we merge again.
                                
                                // If we ever merge us down our path again,
                                // continue with the part we didn't already
                                // eat.
                                node.path_tail = replacement_neighbor.path_tail;
                            } else {
                                // The other possibility is the neighbor is just
                                // us. Officially self loops aren't allowed, so
                                // we censor the edge.
                                
#ifdef debug
                                cerr << "\t\tWe are neighbor (self loop). Hide edge!" << endl;
#endif

                                node.effective_degree--;
                            }
                            
                            // Clean up the neighbor from the to do list; we
                            // finished it without recursing.
                            frame.neighbors.pop_back();
                            
                            // Kick back to the work loop to do the next
                            // neighbor, if any.
                        }
                    } else {
                        // We have returned from a recursive call on this neighbor.
                        
#ifdef debug
                        cerr << "\tReturned from recursion on neighbor " << neighbor_number << endl;
#endif

                        // Support bridge edges: detect if we are returning
                        // across a bridge edge and censor it. Norouzi and Tsin
                        // 2014 as written in the paper assumes no bridge
                        // edges, and what we're about to do relies on all
                        // neighbors connecting back somewhere.
                        if (neighbor.low_point == neighbor.dfs_counter) {
                            // It has no back-edges out of its own subtree, so it must be across a bridge.
#ifdef debug
                            cerr << "\t\tNeighbor is across a bridge edge! Hide edge!" << endl;
#endif
                            
                            // Hide the edge we just took from degree calculations.
                            neighbor.effective_degree--;
                            node.effective_degree--;
                            
                            // Don't do anything else with the edge
                        } else {
                            // Wasn't a bridge edge, so we care about more than just traversing that part of the graph.
                            
                            // Do steps 1.1.1 and 1.1.2 of the algorithm as described in the paper.
                            if (neighbor.effective_degree == 2) {
                                // This neighbor gets absorbed and possibly ejected.
                                
#ifdef debug
                                cerr << "\t\tNeighbor is on a stick" << endl;
                                
                                cerr << "\t\t\tEdge " << frame.current << "-" << neighbor_number << " should never be seen again" << endl;
#endif
                                
                                // Take it off of its own path.
                                neighbor.is_on_path = false;
                                
#ifdef debug
                                cerr << "\t\t\tNew neighbor path: " << path_to_string(neighbor_number) << endl;
#endif
                            }
                            
                            // Because we hid the bridge edges, degree 1 nodes should never happen
                            assert(neighbor.effective_degree != 1);
                            
                            if (node.low_point <= neighbor.low_point) {

#ifdef debug
                                cerr << "\t\tWe have a sufficiently low low point; neighbor comes back in in our subtree" << endl;
                                
                                cerr << "\t\t\tAbsorb us and then the neighbor's path to the end" << endl;
#endif
                                
                                // Absorb all along the path starting with here and
                                // continuing with this neighbor's path, to the
                                // end.
                                absorb_all_along_path(frame.current,
                                                     neighbor_number,
                                                     numeric_limits<number_t>::max()); 
                            } else {
#ifdef debug
                                cerr << "\t\tNeighbor has a lower low point ("
                                    << neighbor.low_point << " < " <<  node.low_point << "); comes back in outside our subtree" << endl;
#endif
                                
                                // Lower our low point to that of the neighbor
                                node.low_point = neighbor.low_point;
                                
#ifdef debug
                                cerr << "\t\t\tNew low point: " << node.low_point << endl;
                                
                                cerr << "\t\t\tAbsorb along path to old low point soure" << endl;
#endif

                                // Absorb all along our own path
                                absorb_all_along_path(numeric_limits<number_t>::max(),
                                                      frame.current,
                                                      numeric_limits<number_t>::max());
                                // Adjust our path to be us and then our neighbor's path
                                node.is_on_path = true;
                                node.path_tail = neighbor_number;
                                
#ifdef debug
                                cerr << "\t\t\tNew path " << path_to_string(frame.current) << endl;
#endif
                            }
                        }
                        
                        // Say we aren't coming back from a recursive call
                        // anymore.
                        frame.recursing = false;
                        
                        // Clean up the neighbor, 
                        frame.neighbors.pop_back();
                        
                        // Kick back to the work loop to do the next neighbor,
                        // if any.
                    }
                    
#ifdef debug
                    cerr << "\tDFS: " << node.dfs_counter
                        << " low point: " << node.low_point
                        << " degree: " << node.effective_degree
                        << " path: " << path_to_string(frame.current) << endl;
#endif
                    
                } else {
                    // All the neighbors left to do for this node are done.
                    
#ifdef debug
                    cerr << "\tNode is visited and no neighbors are on the to do list." << endl;
                    
                    cerr << "\tDFS: " << node.dfs_counter
                        << " low point: " << node.low_point
                        << " degree: " << node.effective_degree
                        << " path: " << path_to_string(frame.current) << endl;
#endif
                    
                    // This node is done.
                    
                    // Remember when we exited it
                    node.dfs_exit = dfs_counter;
                    
                    // Clean up the stack frame.
                    stack.pop_back();
                }
            }
        }
    }
    
    // When we run out of unvisited nodes and the stack is empty, we've
    // completed out search through all connected components of the graph.
}

void three_edge_connected_components_dense(size_t node_count, size_t first_root,
    const function<void(size_t, const function<void(size_t)>&)>& for_each_connected_node,
    const function<void(const function<void(const function<void(size_t)>&)>&)>& component_callback) {
    
    // Make a union-find over all the nodes
    structures::UnionFind uf(node_count, true);
    
    // Call Tsin's Algorithm
    three_edge_connected_component_merges_dense(node_count, first_root, for_each_connected_node, [&](size_t a, size_t b) {
        // When it says to do a merge, do it
        uf.union_groups(a, b);
    });
    
    for (auto& component : uf.all_groups()) {
        // Call the callback for each group
        component_callback([&](const function<void(size_t)>& emit_member) {
            // And whrn it asks for the members
            for (auto& member : component) {
                // Send them all
                emit_member(member);
            }
        });
    }
}

void three_edge_connected_components_dense_cactus(size_t node_count, 
    const function<void(size_t, const function<void(size_t)>&)>& for_each_connected_node,
    const function<void(const function<void(const function<void(size_t)>&)>&)>& component_callback) {
    
    // Use the known good pinchesAndCacti algorithm
    
    // Make the stList of all the vertices, where each vertex is an stList of single element stIntTuple items that point to the ranks of connected nodes.
    // When an item is removed, use the list destructor on it.
    stList* vertices = stList_construct3(0, (void(*)(void *)) stList_destruct);
    
    // TODO: No way to hint final size to the list, and we need the individual member lists to know their destructors for their elements.
    
#ifdef debug
    cerr << "Running Cactus 3ecc on " << node_count << " nodes" << endl;
#endif
    
    for (size_t rank = 0; rank < node_count; rank++) {
        while (rank >= stList_length(vertices)) {
            // Make sure we have an adjacency list allocated for the node
            // When an item in the node's adjacency list is destroyed, run the int tuple destructor.
            stList_append(vertices, stList_construct3(0, (void(*)(void *)) stIntTuple_destruct));
        }
        
        for_each_connected_node(rank, [&](size_t other_rank) {
#ifdef debug
            cerr << "Connect node " << rank << " to node " << other_rank << endl;
#endif
        
            // For each edge on the node, represent it as a 1-tuple in the node's list.
            stList_append((stList*) stList_get(vertices, rank), stIntTuple_construct1((int64_t) other_rank));
            // We don't have to do the back-edge now; we will do it when we visit the other node.
        });
    }
    

#ifdef debug
    for (size_t i = 0; i < stList_length(vertices); i++) {
        cerr << "Vertex " << i << " adjacent to:";
        stList* adjacencies = (stList*) stList_get(vertices, i);
        for (size_t j = 0; j < stList_length(adjacencies); j++) {
            stIntTuple* adj = (stIntTuple*) stList_get(adjacencies, j);
            cerr << " " << stIntTuple_get(adj, 0);
        }
        cerr << endl;
    }
#endif

    // Now we have the graph in the format Tsin's Algorithm wants, so run it.
    // The components come out as a list of lists, one for each component, with
    // the entries in each component's list being 1-element stIntTuples with
    // ranks in them.
    stList* components = computeThreeEdgeConnectedComponents(vertices);
    
#ifdef debug
    cerr << "Got back " << stList_length(components) << " components" << endl;
#endif
    
    for(size_t i = 0; i < stList_length(components); i++) {
        // For each component
        stList* component = (stList*) stList_get(components, i);
        // Announce the component
        component_callback([&](const function<void(size_t)>& visit_member) {
            // And when we get the function to feed the members to
            for (size_t j = 0; j < stList_length(component); j++) {
#ifdef debug
                cerr << "Component " << i << " contains node " << stIntTuple_get((stIntTuple*) stList_get(component, j), 0) << endl;
#endif
            
                // Call it with each member
                visit_member(stIntTuple_get((stIntTuple*) stList_get(component, j), 0));
            }
        });
    }

    // Clean up the component result
    stList_destruct(components);

    // Clean up the vertex data
    stList_destruct(vertices);
}

}
}
