#include "three_edge_connected_components.hpp"

// Grab the Tsin's Algorithm header out of pinchesAndCacti, which installs it into sonLib's include directory
extern "C" {
#include "sonLibList.h"
#include "sonLibTuples.h"
#include "3_Absorb3edge2x.h"
}

#include <structures/union_find.hpp>

#include <limits>
#include <cassert>
#include <iostream>

//#define debug

namespace vg {
namespace algorithms {

using namespace std;

void three_edge_connected_component_merges_dense(size_t node_count, 
    const function<void(size_t, const function<void(size_t)>&)>& for_each_connected_node,
    const function<void(size_t, size_t)>& same_component) {
    
    // Independent Implementation of Tsin and Norouzi (2014) "A simple 3-edge
    // connected component algorithm revisited", which can't reallyu be
    // understood without Tsin (2007) "A Simple 3-Edge-Connected Component
    // Algorithm".
    
    // That algorithm assumes that all bridge edges are removed (i.e.
    // everything is at least 2-connected), but we hack it a bit to generalize
    // to graphs with bridge edges.
    
    // The algorithm does a depth-first search through the graph, and is based
    // on this "absorb-eject" operation. You do it at a node, across ("on") an
    // edge. It (conceptually) steals all the edges from the node at the other
    // end of the edge, deletes the edge, and deletes the other node as well if
    // it has a degree greater than 2. (The original algorithm didn't have to
    // deal with degree 1; here we treat it the same as degree 2 and leave the
    // node floating in its own 3 edge connected component.)
    
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
    // nonempty path that it itself is not on. To support that case, we also
    // give each node a flag for whether it is on its own path.
    
    // In addition to the effective degrees of the nodes, we track a DFS
    // counter for each node saying when in the DFS it was visited, and a "low
    // point", called "lowpt" in the papers, and a visited flag.
    
    // TODO: should we template this on an integer size so we can fit more work
    // in less memory bandwidth when possible?
    using number_t = size_t;
    assert(node_count <= numeric_limits<number_t>::max());
   
    /// This defines the data we track for each node in the graph
    struct TsinNode {
        /// When in the DFS were we visited?
        number_t dfs_counter;
        /// What is our "low point" in the search. I don't quite understand
        /// what a low point is, but the paper explains how to calculate it.
        number_t low_point;
        /// What is the effective degree of this node in the graph with all the
        /// absorb-eject modifications applied?
        number_t effective_degree = 0;
        /// What node has the continuation of this node's path? If equal to
        /// numeric_limits<number_t>::max(), the path ends after here.
        number_t path_tail;
        /// Is this node actually on its own path?
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
        
        cerr << "Absorbing all along path into " << into << " from " << path_start << " to before " << path_past_end << endl;
    
        // Set this to false as soon as we cross an edge
        bool path_null = true;
    
        number_t here = path_start;
        while (here != path_past_end) {
            // Until we hit the end of the path
            
            cerr << "\tAt " << here << endl;
            
            if (here == numeric_limits<number_t>::max()) {
                // We hit the end of the path and never saw path_past_end.
                
                cerr << "\t\tReached path end and missed waypoint" << endl;
                
                // Only allowed if the path was actually edge-free and no merges needed to happen.
                assert(path_null);
                
                cerr << "\t\t\tBut path was empty of edges" << endl;
                
                // Stop now.
                break;
            }
            
            // Find the node we are at
            auto& here_node = nodes[here];
            
            if (here_node.is_on_path) {
                // We're actually on the path.
                
                cerr << "\t\tOn path" << endl;
                
                if (into == numeric_limits<number_t>::max()) {
                    // We haven't found a first node to merge into yet; it is
                    // this one.
                    
                    cerr << "\t\tUse as into" << endl;
                    
                    into = here;
                } else {
                    // We already have a first node to merge into, so merge.
                    
                    cerr << "\t\tMerge with " << into << endl;
                    
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
            
            cerr << "\t\tNext: " << here << endl;
        }
        cerr << "Done absorbing" << endl;
    };
    
    // We need a DFS stack that we manage ourselves, to avoid stack-overflowing
    // as we e.g. walk along big cycles.
    struct DFSStackFrame {
        /// Track the node that this stack frame represents
        number_t current;
        /// Track all the neighbors left to visit.
        /// When we visit a neighbor we pop it off the back.
        vector<number_t> neighbors;
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
    // The paper starts it at 1, but we start it at 0 for integer boundary
    // purposes.
    number_t dfs_counter = 0;
    
    while (next_unvisited != node_count) {
        // We haven't visited everything yet.
        // Stack up the next unvisited node.
        stack.emplace_back();
        stack.back().current = next_unvisited;
        
        cerr << "Root a search at " << next_unvisited << endl;
        
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
                
                cerr << "First visit of node " << frame.current << endl;
                
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
                
                cerr << "\tDFS: " << node.dfs_counter
                    << " low point: " << node.low_point
                    << " degree: " << node.effective_degree
                    << " path: " << node.is_on_path << " " << node.path_tail << endl;
                
                // Stack up all the edges to follow
                for_each_connected_node(frame.current, [&](size_t connected) {
                    frame.neighbors.push_back(connected);
                });
                
                cerr << "\tPut" << frame.neighbors.size() << " edges on to do list" << endl;
                
                // Now we're in a state where we can process edges.
                // So kick back to the work loop as if we just processed an edge.
                continue;
            } else {
                // We have (possibly 0) edges left to do for this node.
                if (!frame.neighbors.empty()) {
                
                    cerr << "Return to node " << frame.current << " with more edges to do" << endl;
                
                    // We have an edge to do!
                    // Look up the neighboring node.
                    number_t neighbor_number = frame.neighbors.back();
                    auto& neighbor = nodes[neighbor_number];
                    
                    if (!frame.recursing) {
                        // This is the first time we are thinking about this neighbor.
                        
                        cerr << "\tThink of neighbor " << neighbor_number << " for the first time" << endl;
                    
                        // Increment degree of the node we're coming from
                        node.effective_degree++;
                        
                        cerr << "\t\tBump degree to " << node.effective_degree << endl;
                        
                        if (!neighbor.visited) {
                            // We need to recurse on this neighbor.
                            
                            cerr << "\t\tRecurse on unvisited neighbor" << endl;
                            
                            // So remember we are recursing.
                            frame.recursing = true;
                            // And set up the recursive frame.
                            stack.emplace_back();
                            stack.back().current = neighbor_number;
                            // Kick back to the work loop; we will see the
                            // unvisited node on top of the stack and do its
                            // visit and add its edges to its to do list.
                            continue;
                        } else {
                            // No need to recurse. But this edge is a back-edge.
                            // Which way are we looking at it?
                            
                            cerr << "\t\tNeighbor has been visited. This is a back-edge." << endl;
                            
                            // Do steps 1.2 and 1.3 from the paper.
                            if (neighbor.dfs_counter < node.dfs_counter) {
                                // The edge to the neighbor is an outgoing
                                // back-edge (i.e. the neighbor was visited
                                // first)
                                
                                cerr << "\t\t\tNeighbor is upstream of us (outgoing back edge)." << endl;
                                
                                if (neighbor.dfs_counter < node.low_point) {
                                    // The neighbor is below our low point.
                                    
                                    cerr << "\t\t\t\tNeighbor has a lower low point" << endl;
                                    
                                    cerr << "\t\t\t\t\tAbsorb along our path to end" << endl;
                                    // Absorb along our whole path.
                                    absorb_all_along_path(numeric_limits<number_t>::max(),
                                                          frame.current,
                                                          numeric_limits<number_t>::max());
                                    
                                    // Adopt the neighbor's DFS counter as our
                                    // new, lower low point.
                                    node.low_point = neighbor.dfs_counter;
                                    cerr << "\t\t\t\t\tNew lower low point " << node.low_point << endl;
                                    
                                    // Our path is now just us.
                                    node.is_on_path = true;
                                    node.path_tail = numeric_limits<number_t>::max();
                                    cerr << "\t\t\t\t\tNew path " << node.is_on_path << " " << node.path_tail << endl;
                                } else {
                                    cerr << "\t\t\t\tWe have a sufficiently low low point" << endl;
                                }
                            } else if (node.dfs_counter < neighbor.dfs_counter) {
                                // The edge to the neighbor is an incoming
                                // back-edge (i.e. we were visited first, but
                                // we recursed into something that got us to
                                // this neighbor already).
                                
                                cerr << "\t\t\tWe are upstream of neighbor (incoming back edge)." << endl;
                                
                                // Drop our effective degree by 2 (I think
                                // we're closing a cycle or something?)
                                node.effective_degree -= 2;
                                cerr << "\t\t\t\tDrop degree to " << node.effective_degree << endl;
                                
                                
                                cerr << "\t\t\t\tAbsorb along our path through neighbor (if not single-node)" << endl;
                                // Absorb along our path from ourselves to the
                                // neighbor, inclusive. We can guarantee that
                                // either the path goes this way, or it's a
                                // single node (which the paper calls "null")
                                // path with no edges.
                                absorb_all_along_path(numeric_limits<number_t>::max(),
                                                      frame.current,
                                                      neighbor.path_tail);
                                
                            }
                            // The other possibility is the neighbor is just
                            // us. Then we don't do anything.
                            
                            // Clean up the neighbor from the to do list; we
                            // finished it without recursing.
                            frame.neighbors.pop_back();
                            
                            // Kick back to the work loop to do the next
                            // neighbor, if any.
                            continue;
                        }
                    } else {
                        // We have returned from a recursive call on this neighbor.
                        
                        cerr << "\tReturned from recursion on neighbor " << neighbor_number << endl;
                        
                        // Do steps 1.1.1 and 1.1.2 of the algorithm as described in the paper.
                        if (neighbor.effective_degree <= 2) {
                            // This neighbor gets absorbed and possibly ejected.
                            
                            cerr << "\t\tNeighbor is on a stick" << endl;
                            
                            cerr << "\t\t\tEdge " << frame.current << "-" << neighbor_number << " should never be seen again" << endl;
                            
                            // Take it off of its own path.
                            neighbor.is_on_path = false;
                            
                            cerr << "\t\t\tNew neighbor path: " << neighbor.is_on_path << " " << neighbor.path_tail << endl;
                            
                        }
                        if (node.low_point <= neighbor.low_point) {
                            cerr << "\t\tWe have a sufficiently low low point" << endl;
                        
                            
                            cerr << "\t\t\t\tAbsorb us and then the neighbor's path to the end" << endl;
                            // Absorb all along the path starting with here and
                            // continuing with this neighbor's path, to the
                            // end.
                            absorb_all_along_path(frame.current,
                                                 neighbor_number,
                                                 numeric_limits<number_t>::max()); 
                        } else {
                            cerr << "\t\tNeighbor has a lower low point" << endl;
                            
                            // Lower our low point to that of the neighbor
                            node.low_point = neighbor.low_point;
                            
                            cerr << "\t\t\tNew low point: " << node.low_point << endl;
                            
                            cerr << "\t\t\tAbsorb along our path to end" << endl;
                            // Absorb all along our own path
                            absorb_all_along_path(numeric_limits<number_t>::max(),
                                                  frame.current,
                                                  numeric_limits<number_t>::max());
                            // Adjust our path to be us and then our neighbor's path
                            node.is_on_path = true;
                            node.path_tail = neighbor_number;
                            cerr << "\t\t\tNew path " << node.is_on_path << " " << node.path_tail << endl;
                        }
                        
                        // Say we aren't coming back from a recursive call
                        // anymore.
                        frame.recursing = false;
                        
                        // Clean up the neighbor, 
                        frame.neighbors.pop_back();
                        
                        // Kick back to the work loop to do the next neighbor,
                        // if any.
                        continue;
                    }
                } else {
                    // All the neighbors left to do for this node are done.
                    
                    cerr << "\tNode is visited and no neighbors are on the to do list." << endl;
                    
                    // This node is done. Clean up the stack frame.
                    stack.pop_back();
                }
            }
        }
    }
    
    // When we run out of unvisited nodes and the stack is empty, we've
    // completed out search through all connected components of the graph.
}

void three_edge_connected_components_dense(size_t node_count, 
    const function<void(size_t, const function<void(size_t)>&)>& for_each_connected_node,
    const function<void(const function<void(const function<void(size_t)>&)>&)>& component_callback) {
    
    // Make a union-find over all the nodes
    structures::UnionFind uf(node_count, true);
    
    // Call Tsin's Algorithm
    three_edge_connected_component_merges_dense(node_count, for_each_connected_node, [&](size_t a, size_t b) {
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

}
}
