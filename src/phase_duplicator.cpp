#include "phase_duplicator.hpp"

namespace vg {

using namespace std;

PhaseDuplicator::PhaseDuplicator(const xg::XG& index) : index(index) {
    // Nothing to do!
}

pair<Graph, vector<Translation>> PhaseDuplicator::duplicate(const set<id_t>& subgraph, id_t& next_id) const {

    // Allocate space for the result
    pair<Graph, vector<Translation>> to_return;
    // Name the pair members
    Graph& duplicated = to_return.first;
    vector<Translation>& translations = to_return.second;
    
    // The duplicated graph is produced assuming the original nodes and their edges are all deleted.
    
    // Find all the traversals
    vector<pair<xg::XG::thread_t, int>> observed_threads = list_haplotypes(subgraph);
    
    for (auto& thread_and_count : observed_threads) {
        // For every unique haplotype in the subgraph
        
        // Grab it
        auto& thread = thread_and_count.first;
        
        // We want a translation
        translations.emplace_back();
        Translation& translation = translations.back();
        
        // As we scan along this thread, what was the last node we created?
        Node* last_node = nullptr;
        // Were we visiting that node in reverse?
        bool last_is_reverse = false;
        
        for (size_t i = 0; i < thread.size(); i++) {
            // For every visit in the thread...
            auto& visit = thread[i];
        
            // We'll duplicate the original node in its original orientation for
            // each visit.
            
            // Grab a new ID
            id_t new_id = ++next_id;
            
            // Make the node
            Node* node = duplicated.add_node();
            node->set_sequence(index.node_sequence(visit.node_id));
            node->set_id(new_id);
            
            if (last_node != nullptr) {
                // Add a connecting edge to the previous node
                Edge* edge = duplicated.add_edge();
                edge->set_from(last_node->id());
                edge->set_from_start(last_is_reverse);
                edge->set_to(node->id());
                edge->set_to_end(visit.is_reverse);
            }
            last_node = node;
            last_is_reverse = visit.is_reverse;
            
            
            if (i == 0 || i + 1 == thread.size()) {
                // We're at the start or end of a traversal, so we might be at a
                // border where we want to enter or leave the subgraph.
                
                // Get the edges on the outside ends of this node that cross borders
                vector<Edge> border_edges;
                
                if (i == 0) {
                    // If we're on the left edge of the thread, look at border-crossing edges on our left edge
                    auto off_left = find_border_edges(visit, true, subgraph);
                    copy(off_left.begin(), off_left.end(), back_inserter(border_edges));
                }
                
                if (i + 1 == thread.size()) {
                    // If we're on the right edge of the thread, look at border-crossing edges on our right edge
                    auto off_right = find_border_edges(visit, false, subgraph);
                    copy(off_right.begin(), off_right.end(), back_inserter(border_edges));
                }

                for (auto& border_edge : border_edges) {
                    // Make a new duplicate of each edge
                    Edge* new_edge = duplicated.add_edge();
                    *new_edge = border_edge;
                    
                    if (new_edge->from() == visit.node_id) {
                        // Update the from if necessary
                        new_edge->set_from(node->id());
                    }
                    
                    if (new_edge->to() == visit.node_id) {
                        // Update the to if necessary
                        new_edge->set_to(node->id());
                    }
                    
                }
            }
            
            // Add the pair of Mappings to the translation paths
            // We go "from" the original graph "to" the duplicated graph.
            Mapping* to_mapping = translation.mutable_to()->add_mapping();
            to_mapping->mutable_position()->set_node_id(new_id);
            to_mapping->mutable_position()->set_is_reverse(visit.is_reverse);
            Edit* to_edit = to_mapping->add_edit();
            to_edit->set_from_length(node->sequence().size());
            to_edit->set_to_length(node->sequence().size());
            
            Mapping* from_mapping = translation.mutable_from()->add_mapping();
            from_mapping->mutable_position()->set_node_id(visit.node_id);
            from_mapping->mutable_position()->set_is_reverse(visit.is_reverse);
            // Just copy the edit data because it's the same in both paths
            *(from_mapping->add_edit()) = *to_edit;
        }
    }
    
    // Ship out the result
    return to_return;
}

vector<pair<xg::XG::thread_t, int>> PhaseDuplicator::list_haplotypes(const set<id_t>& subgraph) const {
    // We need to look at all the nodes in the subgraph and extract all the
    // unique traversals not contained in another traversal.
    
    // What we need to do is look at all the threads coming in through border
    // nodes, and all the threads starting at interior nodes, and then orient
    // and deduplicate.
    
    // Get all the border traversals, and put them in a set for querying
    set<xg::XG::ThreadMapping> borders {find_borders(subgraph)};
    
    // This holds all the traversals observed (full and partial) in their canonical orientations, with their total counts
    map<xg::XG::thread_t, int> observed_traversals;
    
    for (auto& start : borders) {
        // For every incoming traversal of a node, get the haplotypes that go through it.
        vector<pair<xg::XG::thread_t, int>> threads_and_counts = list_haplotypes_through(start, subgraph);
        
        for (auto& found : threads_and_counts) {
            // Add each found traversal to the set in its canonical orientation
            observed_traversals[canonicalize(found.first)] += found.second;
        }
    }
    
    for (auto& start_node : subgraph) {
        // For every node we could start at
        for (auto& orientation : {false, true}) {
            // In every orientation
            
            // Consider starting there
            xg::XG::ThreadMapping start {start_node, orientation};
            
            if (borders.count(start)) {
                // We already handled all the threads coming through this
                // orientation, starting there or not.
                continue;
            }
            
            // Get the haplotypes that start with a traversal of this node
            vector<pair<xg::XG::thread_t, int>> threads_and_counts = list_haplotypes_from(start, subgraph);
            
            for (auto& found : threads_and_counts) {
                // Add each found traversal to the set in its canonical orientation.
                // TODO: should we remove prefixes or suffixes of existing traversals?
                observed_traversals[canonicalize(found.first)] += found.second;
            }
            
        }
    }
    
    // Vectorize the map
    // TODO: skip this copy somehow (maybe use sets throughout?)
    vector<pair<xg::XG::thread_t, int>> flat_traversals;
    for (auto& kv : observed_traversals) {
        // Blit all the pairs over
        flat_traversals.push_back(kv);
    }
    
    return flat_traversals;
}

vector<pair<xg::XG::thread_t, int>> PhaseDuplicator::list_haplotypes_through(xg::XG::ThreadMapping start_node,
    const set<id_t>& subgraph) const {
    
    // We want all haplotypes that start with a traversal through the given
    // node.
    
    // Start a search at the first node traversal
    xg::XG::ThreadSearchState start_state;
    index.extend_search(start_state, start_node);
    
    // Call the search starting from a state selecting everything on the side
    return list_haplotypes(start_node, start_state, subgraph);
}

vector<pair<xg::XG::thread_t, int>> PhaseDuplicator::list_haplotypes_from(xg::XG::ThreadMapping start_node,
    const set<id_t>& subgraph) const {
    
    // We want all haplotypes that start their threads actually at this side of this node.
    xg::XG::ThreadSearchState start_state = index.select_starting(start_node);
    
#ifdef debug
    cerr << "Starting search for threads from " << start_node.node_id << " "
        << start_node.is_reverse << " with range " << start_state.current_side 
        << ", " << start_state.range_start << ", " << start_state.range_end << endl;
#endif
    
    // Call the search starting with only the threads that begin here selected.
    return list_haplotypes(start_node, start_state, subgraph);
}
    

vector<pair<xg::XG::thread_t, int>> PhaseDuplicator::list_haplotypes(xg::XG::ThreadMapping start_node,
    xg::XG::ThreadSearchState start_state, const set<id_t>& subgraph) const {
    // Adapted from Yohei's "haplotype_extracter"
    
    // This holds extracted sub-haplotypes and the search states for looking right off their ends.
    vector<pair<xg::XG::thread_t, xg::XG::ThreadSearchState>> search_intermediates;
    // This holds the results we will return, as pairs of extracted haplotypes and haplotype counts
    vector<pair<xg::XG::thread_t, int>> search_results;
    
    if (start_state.is_empty()) {
        // Nothing visits the first node in the right direction
        return search_results;
    }
    // Otherwise, start with this as the first intermediate result
    search_intermediates.push_back(make_pair(xg::XG::thread_t{start_node}, start_state));
    
    while (search_intermediates.size() > 0) {
        // Until we've finished all our intermediates...
        // Pop one off
        pair<xg::XG::thread_t, xg::XG::ThreadSearchState> last(move(search_intermediates.back()));
        search_intermediates.pop_back();
        
        // Break out contents
        auto& last_thread = last.first;
        auto& last_node = last_thread.back();
        auto& last_state = last.second;
        
        // Remember how many we had besides this one, so we can see if we added
        // any descendants of this one, or if this one is a dead end.
        auto check_size = search_intermediates.size();
        
        // Grab the edges out of the last thing in the partial haplotype
        vector<Edge> edges = last_node.is_reverse ?
            index.edges_on_start(last_node.node_id) :
            index.edges_on_end(last_node.node_id);
            
        // Keep track of all the successful extensions
        vector<pair<xg::XG::ThreadMapping, xg::XG::ThreadSearchState>> extensions;
            
        for (int i = 0; i < edges.size(); i++) {
            // Follow each of them
            xg::XG::ThreadMapping next_node = traverse_edge(edges[i], last_node);
            
            if (!subgraph.count(next_node.node_id)) {
                // Skip edges that leave the subgraph
                continue;
            }
            
            // Try searching with where the edge goes
            xg::XG::ThreadSearchState new_state = last_state;
            index.extend_search(new_state, next_node);
            
            if (!new_state.is_empty()) {
                // We found something. Remember it
                extensions.push_back(make_pair(next_node, new_state));                
            }
        }
        
        while (extensions.size() > 1) {
            // Handle all but one extension by copying the thread
            xg::XG::thread_t new_thread(last_thread);
            // Stick the new ThreadMapping onto the end of the copy
            new_thread.push_back(extensions.back().first);
            // Move the longer thread and the new search state into the intermediates
            search_intermediates.emplace_back(move(new_thread), extensions.back().second);
            
            // We finished this extension
            extensions.pop_back();
        }
        if (extensions.size() == 1) {
            // Only one remaining successful extension, so we can just use the old thread
            last_thread.push_back(extensions.back().first);
            search_intermediates.emplace_back(move(last_thread), extensions.back().second);
        } else {
            // We couldn't find anywhere inside the subgraph that anyone
            // actually goes. So we have to end this traversal here. Send the
            // last thread out as the search result, together with its count.
            search_results.emplace_back(move(last_thread), last_state.count());
        }
    }
    return search_results;
}

set<xg::XG::ThreadMapping> PhaseDuplicator::find_borders(const set<id_t> subgraph) const {
    // We want all the traversals along which we can enter the subgraph.
    
    // We'll put them all in here
    set<xg::XG::ThreadMapping> borders;
    
    for (auto& node_id : subgraph) {
        // For every node
        for (auto is_reverse : {false, true}) {
            // In every direction
            
            // Get the *previous* edges that we could have read through to hit
            // this node in this orientation.
            vector<Edge> edges = is_reverse ?
                index.edges_on_end(node_id) :
                index.edges_on_start(node_id);
                
            // We set this to true if a predecessor edge escapes from the subgraph.
            bool is_border = false;
            for (auto& edge : edges) {
                // Look at each edge
                if (!subgraph.count(edge.from()) || !subgraph.count(edge.to())) {
                    // As soon as we find one that escapes the subgraph, we know we have a border here
                    is_border = true;
                    break;
                }
            }
            
            if (is_border) {
                // We can read into the subgraph and hit this traversal.
                xg::XG::ThreadMapping mapping {node_id, is_reverse};
                borders.insert(mapping);
            }
        }
    }
    
    return borders;
}

vector<Edge> PhaseDuplicator::find_border_edges(xg::XG::ThreadMapping mapping, bool on_start, const set<id_t>& subgraph) const {
        
    // We'll fill this in with the border edges
    vector<Edge> border_edges;
        
    // Get the edges on the start or end of the given ThreadMapping, as
    // appropriate.
    vector<Edge> edges = on_start != mapping.is_reverse ?
        index.edges_on_start(mapping.node_id) :
        index.edges_on_end(mapping.node_id);
        
    for (auto& edge : edges) {
        // Look at each edge
        if (!subgraph.count(edge.from()) || !subgraph.count(edge.to())) {
            // If either end escapes the subgraph, it's a border edge
            border_edges.push_back(edge);
        }
    }
    
    return border_edges;
    
}

xg::XG::ThreadMapping PhaseDuplicator::traverse_edge(const Edge& e, const xg::XG::ThreadMapping& prev) {
    // We're going to traverse this edge. We need to know if we take it forward
    // or backward. This is true if, reading through the prev node in the given
    // orientation, we could traverse the edge forwards. We may also be able to
    // traverse it backwards, but in that case it's a single-side self loop, so
    // direction doesn't really matter.
    bool traverse_forward = (prev.node_id == e.from() && prev.is_reverse == e.from_start());
    
    // We'll fill this in with where we land
    xg::XG::ThreadMapping next;
    
    if (traverse_forward) {
        // Land at the end of the edge
        next.node_id = e.to();
        next.is_reverse = e.to_end();
    } else {
        // Land at the start of the edge
        next.node_id = e.from();
        next.is_reverse = !e.from_start();
    }
    
    return next;
}

xg::XG::thread_t PhaseDuplicator::canonicalize(const xg::XG::thread_t& thread) {
    xg::XG::thread_t reversed;
    
    for (auto i = thread.rbegin(); i != thread.rend(); ++i) {
        // Reverse every mapping in the original thread
        xg::XG::ThreadMapping reversed_mapping = *i;
        reversed_mapping.is_reverse = !reversed_mapping.is_reverse;
        
        // And put it in the new one
        reversed.push_back(reversed_mapping);
    }
    
    if (reversed < thread) {
        // The reversed version compares smaller and is thus canonical.
        return reversed;
    } else {
        // The original compares smaller and is thus canonical.
        return thread;
    }
}

} 
