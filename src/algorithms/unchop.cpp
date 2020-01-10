/**
 * \file unchop.cpp
 *
 * Defines an algorithm to join adjacent handles.
 */

#include "unchop.hpp"

#include <unordered_set>
#include <list>
#include <set>
#include <iostream>
#include <sstream>

namespace vg {
namespace algorithms {

using namespace std;

/// Return true if nodes share all paths and the mappings they share in these paths
/// are adjacent, in the specified relative order and orientation.
static bool nodes_are_perfect_path_neighbors(PathHandleGraph* graph, handle_t left_handle, handle_t right_handle) {
    
#ifdef debug
    cerr << "Check if " << graph->get_id(left_handle) << (graph->get_is_reverse(left_handle) ? "-" : "+") << " and "
        << graph->get_id(right_handle) << (graph->get_is_reverse(right_handle) ? "-" : "+") << " are perfect path neighbors" << endl;
#endif
    
    // Set this false if we find an impermissible step
    bool ok = true;
    
    // Count the number of permissible steps on the next node we find
    size_t expected_next = 0;
    
    graph->for_each_step_on_handle(left_handle, [&](const step_handle_t& here) {
        // For each path step on the left
        
        // We need to work out if the path traverses this handle backward.
        bool step_is_to_reverse_of_handle = (graph->get_handle_of_step(here) != left_handle);
        
#ifdef debug
        cerr << "Consider visit of path " << graph->get_path_name(graph->get_path_handle_of_step(here))
            << " to " << (step_is_to_reverse_of_handle ? "reverse" : "forward") << " orientation of handle" << endl;
#endif
        
        if (!(step_is_to_reverse_of_handle ? graph->has_previous_step(here) : graph->has_next_step(here))) {
            // If there's no visit to the right of this handle, it can't be to the right next place
#ifdef debug
            cerr << "Path stops here so no subsequent handle is a perfect path neighbor" << endl;
#endif
            ok = false;
            return false;
        }
        
        // Walk along the path whatever direction is forward relative to our left handle.
        step_handle_t step_to_right = step_is_to_reverse_of_handle ? graph->get_previous_step(here) : graph->get_next_step(here);
        handle_t handle_to_right = graph->get_handle_of_step(step_to_right);
        if (step_is_to_reverse_of_handle) {
            // If we had to walk back along the reverse strand of the path, we have to flip where we ended up.
            handle_to_right = graph->flip(handle_to_right);
        }
        
        if (handle_to_right != right_handle) {
            // It goes to the wrong next place
            
#ifdef debug
            cerr << "Path goes to the wrong place ("
                << graph->get_id(handle_to_right) << (graph->get_is_reverse(handle_to_right) ? "-" : "+")
                << ") and so these nodes are not perfect path neighbors"  << endl;
#endif
            
            ok = false;
            return false;
        }
        
        // Otherwise, record a step that is allowed to exist on the next handle
        expected_next++;
        return true;
    });
    
    if (!ok) {
        // We found a bad step, or the path stopped.
        return false;
    }
    
    // Now count up the path steps on the right handle
    size_t observed_next = 0;
    graph->for_each_step_on_handle(right_handle, [&](const step_handle_t& ignored) {
#ifdef debug
        cerr << "Next node has path " << graph->get_path_name(graph->get_path_handle_of_step(ignored)) << endl;
#endif
        observed_next++;
    });
    
#ifdef debug
    if (observed_next != expected_next) {
        cerr << "Next node has " << observed_next << " path visits but should have " << expected_next << endl;
    }
#endif
    
    // If there are any steps on the right node that weren't accounted for on
    // the left node, fail. Otherwise, succeed.
    return observed_next == expected_next;

}

/// Get the set of components that could be merged into single nodes without
/// changing the path space of the graph. Emits oriented traversals of
/// nodes, in the order and orientation in which they are to be merged.
static set<list<handle_t>> simple_components(PathHandleGraph* graph, int min_size = 1) {

    // go around and establish groupings
    unordered_set<id_t> seen;
    set<list<handle_t>> components;
    graph->for_each_handle([&](const handle_t& n) {
            id_t n_id = graph->get_id(n);
    
            if (seen.count(n_id)) {
                // We already found this node in a previous component
                return;
            }
            
#ifdef debug
            cerr << "Component based on " << n_id << endl;
#endif
            
            seen.insert(n_id);
            // go left and right through each as far as we have only single edges connecting us
            // to nodes that have only single edges coming in or out
            // that go to other nodes
            list<handle_t> c;
            // go left
            {
                
                handle_t l = n;
                
                vector<handle_t> prev;
                graph->follow_edges(l, true, [&](const handle_t& h) {
                    prev.push_back(h);
                });
                
#ifdef debug
                cerr << "\tLeft: ";
                for (auto& x : prev) {
                    cerr << graph->get_id(x) << (graph->get_is_reverse(x) ? '-' : '+')
                        << "(" << graph->get_degree(x, false) << " edges right) ";
                    graph->follow_edges(x, false, [&](const handle_t& other) {
                        cerr << "(to " << graph->get_id(other) << (graph->get_is_reverse(other) ? '-' : '+') << ") ";
                    });
                }
                cerr << endl;
#endif
                while (prev.size() == 1
                       && graph->get_degree(prev.front(), false) == 1) {   
                       
                    // While there's only one node left of here, and one node right of that node...
                    auto last = l;
                    // Move over left to that node
                    l = prev.front();
                    // avoid merging if it breaks stored paths
                    if (!nodes_are_perfect_path_neighbors(graph, l, last)) {
#ifdef debug
                        cerr << "\tNot perfect neighbors!" << endl;
#endif
                        break;
                    }
                    // avoid merging if it's already in this or any other component (catch self loops)
                    if (seen.count(graph->get_id(l))) {
#ifdef debug
                        cerr << "\tAlready seen!" << endl;
#endif
                        break;
                    }
                    
                    prev.clear();
                    graph->follow_edges(l, true, [&](const handle_t& h) {
                        prev.push_back(h);
                    });
                    
#ifdef debug
                    cerr << "\tLeft: ";
                    for (auto& x : prev) {
                        cerr << graph->get_id(x) << (graph->get_is_reverse(x) ? '-' : '+')
                            << "(" << graph->get_degree(x, false) << " edges right) ";
                        graph->follow_edges(x, false, [&](const handle_t& other) {
                            cerr << "(to " << graph->get_id(other) << (graph->get_is_reverse(other) ? '-' : '+') << ") ";
                        });
                    }
                    cerr << endl;
#endif
                    c.push_front(l);
                    seen.insert(graph->get_id(l));
                }
            }
            // add the node (in the middle)
            c.push_back(n);
            // go right
            {
                handle_t r = n;
                vector<handle_t> next;
                graph->follow_edges(r, false, [&](const handle_t& h) {
                    next.push_back(h);
                });
                
#ifdef debug
                cerr << "\tRight: ";
                for (auto& x : next) {
                    cerr << graph->get_id(x) << (graph->get_is_reverse(x) ? '-' : '+')
                        << "(" << graph->get_degree(x, true) << " edges left) ";
                    graph->follow_edges(x, true, [&](const handle_t& other) {
                        cerr << "(to " << graph->get_id(other) << (graph->get_is_reverse(other) ? '-' : '+') << ") ";
                    });
                }
                cerr << endl;
#endif
                while (next.size() == 1
                       && graph->get_degree(next.front(), true) == 1) {   
                       
                    // While there's only one node right of here, and one node left of that node...
                    auto last = r;
                    // Move over right to that node
                    r = next.front();
                    // avoid merging if it breaks stored paths
                    if (!nodes_are_perfect_path_neighbors(graph, last, r)) {
#ifdef debug
                        cerr << "\tNot perfect neighbors!" << endl;
#endif
                        break;
                    }
                    // avoid merging if it's already in this or any other component (catch self loops)
                    if (seen.count(graph->get_id(r))) {
#ifdef debug
                        cerr << "\tAlready seen!" << endl;
#endif
                        break;
                    }
                    
                    next.clear();
                    graph->follow_edges(r, false, [&](const handle_t& h) {
                        next.push_back(h);
                    });
                    
#ifdef debug
                    cerr << "\tRight: ";
                    for (auto& x : next) {
                        cerr << graph->get_id(x) << (graph->get_is_reverse(x) ? '-' : '+')
                            << "(" << graph->get_degree(x, true) << " edges left) ";
                        graph->follow_edges(x, true, [&](const handle_t& other) {
                            cerr << "(to " << graph->get_id(other) << (graph->get_is_reverse(other) ? '-' : '+') << ") ";
                        });
                    }
                    cerr << endl;
#endif
                    c.push_back(r);
                    seen.insert(graph->get_id(r));
                }
            }
            if (c.size() >= min_size) {
                components.emplace(std::move(c));
            }
        });
#ifdef debug
    cerr << "components " << endl;
    for (auto& c : components) {
        for (auto x : c) {
            cerr << graph->get_id(x) << (graph->get_is_reverse(x) ? '-' : '+') << " ";
        }
        cerr << endl;
    }
#endif
    return components;
}

/// Concatenates the nodes into a new node with the same external linkage as
/// the provided component. All handles must be in left to right order and a
/// consistent orientation. All paths present must run all the way through the
/// run of nodes from start to end or end to start.
///
/// Returns the handle to the newly created node.
///
/// After calling this on a vg::VG, paths will be invalid until
/// Paths::compact_ranks() is called.
static handle_t concat_nodes(handlegraph::MutablePathDeletableHandleGraph* graph, const list<handle_t>& nodes) {

    // Make sure we have at least 2 nodes
    assert(!nodes.empty() && nodes.front() != nodes.back());
    
#ifdef debug
    cerr << "Paths before concat: " << endl;
   
    graph->for_each_path_handle([&](const path_handle_t p) {
        cerr << graph->get_path_name(p) << ": ";
        for (auto h : graph->scan_path(p)) {
            cerr << graph->get_id(h) << (graph->get_is_reverse(h) ? '-' : '+') << " ";
        }
        cerr << endl;
    });
   
#endif

    // We also require no edges enter or leave the run of nodes, but we can't check that now.
    
    // Make the new node
    handle_t new_node;
    {
        stringstream ss;
        for (auto& n : nodes) {
            ss << graph->get_sequence(n);
        }
        
        new_node = graph->create_handle(ss.str());
    }
    
#ifdef debug
    cerr << "Concatenating ";
    for (auto& n : nodes) {
        cerr << graph->get_id(n) << (graph->get_is_reverse(n) ? "-" : "+") << " ";
    }
    cerr << "into " << graph->get_id(new_node) << "+" << endl;
#endif
    
    // We should be able to rely on our handle graph to deduplicate edges, but see https://github.com/vgteam/libbdsg/issues/39
    // So we deduplicate ourselves.
    
    // Find all the neighbors. Make sure to translate edges to the other end of
    // the run, or self loops.
    unordered_set<handle_t> left_neighbors;
    graph->follow_edges(nodes.front(), true, [&](const handle_t& left_neighbor) {
        if (left_neighbor == nodes.back()) {
            // Loop back to the end
            left_neighbors.insert(new_node);
        } else if (left_neighbor == graph->flip(nodes.front())) {
            // Loop back to the front
            left_neighbors.insert(graph->flip(new_node));
        } else {
            // Normal edge
            left_neighbors.insert(left_neighbor);
        }
    });
    
    unordered_set<handle_t> right_neighbors;
    graph->follow_edges(nodes.back(), false, [&](const handle_t& right_neighbor) {
        if (right_neighbor == nodes.front()) {
            // Loop back to the front.
            // We will have seen it from the other side, so ignore it here.
        } else if (right_neighbor == graph->flip(nodes.back())) {
            // Loop back to the end
            right_neighbors.insert(graph->flip(new_node));
        } else {
            // Normal edge
            right_neighbors.insert(right_neighbor);
        }
    });
    
    // Make all the edges, now that we can't interfere with edge listing
    for (auto& n : left_neighbors) {
#ifdef debug
        cerr << "Creating edge " << graph->get_id(n) << (graph->get_is_reverse(n) ? "-" : "+") << " -> "
            <<  graph->get_id(new_node) << (graph->get_is_reverse(new_node) ? "-" : "+") << endl;
#endif
        graph->create_edge(n, new_node);
    }
    for (auto& n : right_neighbors) {
    
#ifdef debug
        cerr << "Creating edge " << graph->get_id(new_node) << (graph->get_is_reverse(new_node) ? "-" : "+") << " -> "
            <<  graph->get_id(n) << (graph->get_is_reverse(n) ? "-" : "+") << endl;
#endif
    
        graph->create_edge(new_node, n);
    }
    
    {
        // Collect the first and last visits along paths. TODO: this requires
        // walking each path all the way through.
        //
        // This contains the first and last handles in path orientation, and a flag
        // for if the path runs along the reverse strand of our run of nodes.
        vector<tuple<step_handle_t, step_handle_t, bool>> ranges_to_rewrite;
        
        graph->for_each_step_on_handle(nodes.front(), [&](const step_handle_t& front_step) {
            auto path = graph->get_path_handle_of_step(front_step);
#ifdef debug
            cerr << "Consider path " << graph->get_path_name(path) << endl;
#endif
        
            // If we don't get the same oriented node as where this step is
            // stepping, we must be stepping on the other orientation.
            bool runs_reverse = (graph->get_handle_of_step(front_step) != nodes.front());
       
            step_handle_t back_step = front_step;
            
            while(graph->get_handle_of_step(back_step) != (runs_reverse ? graph->flip(nodes.back()) : nodes.back())) {
                // Until we find the step on the path that visits the last node in our run.
                // Go along the path towards where our last node should be, in our forward orientation.
                back_step = runs_reverse ? graph->get_previous_step(back_step) : graph->get_next_step(back_step);
            }
            
            // Now we can record the range to rewrite
            // Make sure to put it into path-forward order
            if (runs_reverse) {
#ifdef debug
                cerr << "\tGoing to rewrite between " << graph->get_id(graph->get_handle_of_step(front_step)) << " and " << graph->get_id(graph->get_handle_of_step(back_step)) << " backward" << endl;
#endif
                ranges_to_rewrite.emplace_back(back_step, front_step, true);
            } else {
            
#ifdef debug
                cerr << "\tGoing to rewrite between " << graph->get_id(graph->get_handle_of_step(front_step)) << " and " << graph->get_id(graph->get_handle_of_step(back_step)) << endl;
#endif
                ranges_to_rewrite.emplace_back(front_step, back_step, false);
            }
        });
        
        for (auto& range : ranges_to_rewrite) {
            // Rewrite each range to visit the new node in the appropriate orientation instead of whatever it did before
            // Make sure to advance the end of the range because rewrite is end-exclusive (to allow insert).
            graph->rewrite_segment(get<0>(range), graph->get_next_step(get<1>(range)), {get<2>(range) ? graph->flip(new_node) : new_node});
        }
    }
    
    // Destroy all the old edges
    // We know they only exist to the left and right neighbors, and along the run
    for (auto& n : left_neighbors) {
        graph->destroy_edge(n, nodes.front());
    }
    for (auto& n : right_neighbors) {
        graph->destroy_edge(nodes.back(), n);
    }
    auto it = nodes.begin();
    auto next_it = it;
    ++next_it;
    while (next_it != nodes.end()) {
        graph->destroy_edge(*it, *next_it);
        it = next_it;
        ++next_it;
    }
    
    for (auto& n : nodes) {
        // Destroy all the old nodes
#ifdef debug
        cerr << "Destroying node " << graph->get_id(n) << endl;
#endif
        graph->destroy_handle(n);
    }

#ifdef debug
    cerr << "Paths after concat: " << endl;
   
    graph->for_each_path_handle([&](const path_handle_t p) {
        cerr << graph->get_path_name(p) << ": ";
        for (auto h : graph->scan_path(p)) {
            cerr << graph->get_id(h) << (graph->get_is_reverse(h) ? '-' : '+') << " ";
        }
        cerr << endl;
    });
   
#endif

    // Return the new handle we merged to.
    return new_node;
}

void unchop(handlegraph::MutablePathDeletableHandleGraph* graph) {
#ifdef debug
    cerr << "Running unchop" << endl;
#endif
    for (auto& comp : simple_components(graph, 2)) {
#ifdef debug
        cerr << "Unchop " << comp.size() << " nodes together" << endl;
#endif
        concat_nodes(graph, comp);
    }
}


}
}

