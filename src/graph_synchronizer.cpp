#include "graph_synchronizer.hpp"
#include "algorithms/extract_connecting_graph.hpp"

#include <iterator>

namespace vg {

using namespace std;

GraphSynchronizer::GraphSynchronizer(VG& graph) : graph(graph) {
    // Because in general paths can overlap each other, and because we can't
    // build a path index after a path has been modified (since we don't keep
    // the ranks up to date internally), we need to build all the indexes up
    // front, even if we're just working on a single path.
    graph.for_each_path_handle([&](const path_handle_t& path) {
        string name = graph.get_path_name(path);
        if (!Paths::is_alt(name)) {
            // We only care about reference paths.
            get_path_index(name);
        }
    });
}

void GraphSynchronizer::with_path_index(const string& path_name, const function<void(const PathIndex&)>& to_run) {
    // Get a reader lock on the graph
    std::lock_guard<std::mutex> guard(whole_graph_lock);
    to_run(get_path_index(path_name));
}

const string& GraphSynchronizer::get_path_sequence(const string& path_name) {
    // Lock the whole graph
    std::lock_guard<std::mutex> guard(whole_graph_lock);
    
    // Get (and possibly generate from the graph) the index, and return its
    // sequence string (which won't change)
    return get_path_index(path_name).sequence;
}

    
// We need a function to grab the index for a path
PathIndex& GraphSynchronizer::get_path_index(const string& path_name) {
    
    // We don't work on alt paths; there could be too many to pre-index.
    assert(!Paths::is_alt(path_name));
    
    if (!indexes.count(path_name)) {
        // Not already made. Generate it.
        indexes.emplace(piecewise_construct,
            forward_as_tuple(path_name), // Make the key
            forward_as_tuple(graph, path_name, true)); // Make the PathIndex
    }
    return indexes.at(path_name);
}

void GraphSynchronizer::update_path_indexes(const vector<Translation>& translations) {
    for (auto& kv : indexes) {
        // We need to touch every index (IN PLACE!)
        
        // Feed each index all the translations, which it will parse into node-
        // partitioning translations and then apply.
        kv.second.apply_translations(translations);
    }
}

GraphSynchronizer::Lock::Lock(GraphSynchronizer& synchronizer,
    const string& path_name, size_t path_offset, size_t context_bases, bool reflect) : 
    synchronizer(synchronizer), path_name(path_name), path_offset(path_offset), 
    context_bases(context_bases), reflect(reflect) {
    
    // Nothing to do. We've saved all the details on the request.
}

GraphSynchronizer::Lock::Lock(GraphSynchronizer& synchronizer, const string& path_name, size_t start, size_t past_end) :
    synchronizer(synchronizer), path_name(path_name), start(start), past_end(past_end) {

    // Nothing to do. We've saved all the details on the request.
}

void GraphSynchronizer::Lock::lock() {
    // Now we have to block until a lock is obtained.
    
    if (!locked_nodes.empty()) {
        // We already have a lock
        return;
    }
    
    // What we do is, we lock the graph and wait on the condition variable, with
    // the check code being that we find the subgraph and immediate neighbors
    // and verify none of its nodes are locked
    
    // Lock the whole graph
    std::unique_lock<std::mutex> lk(synchronizer.whole_graph_lock);
    synchronizer.wait_for_region.wait(lk, [&]{
        // Now we have exclusive use of the graph and indexes, and we need to
        // see if anyone else is using any nodes we need.
        
        
        
        // Extract the context around that node
        VG context;
        
        if (start != 0 || past_end != 0) {
            // We want to extract a range
            
            // Find the outer ends of this range
            NodeSide start_left = synchronizer.get_path_index(path_name).at_position(start);
            NodeSide end_right = synchronizer.get_path_index(path_name).at_position(past_end == 0 ? 0 : past_end - 1).flip();
            
            // Fill in the endpoints pair
            endpoints = make_pair(start_left, end_right);
            
#ifdef debug
            cerr << "Endpoints: " << start_left << ", " << end_right << endl;
            
            // Trace the path in the index to say what should be found
            cerr << "Path: " << endl;
            auto it = synchronizer.get_path_index(path_name).find_position(start);
            while(it != synchronizer.get_path_index(path_name).end() && it->second.flip() != end_right) {
                cerr << "\tVisit " << it->second;
                auto it2 = it;
                ++it2;
                if (it2 != synchronizer.get_path_index(path_name).end()) {
                    // Make sure we have an edge from this node to the next on the path
                    assert(synchronizer.graph.has_edge(it->second.flip(), it2->second));
                    cerr << " " << pb2json(*synchronizer.graph.get_edge(it->second.flip(), it2->second));
                }
                ++it;
                cerr << endl;
            }
#endif

            // Make them into pos_ts that point left to right, the way Jordan thinks.
            pos_t left_pos = make_pos_t(start_left.node, start_left.is_end, 0);
            pos_t right_pos = make_pos_t(end_right.node, !end_right.is_end,
                synchronizer.graph.get_node(end_right.node)->sequence().size());
            
            // Since these are already at node ends, we don't need to worry about node cuts.
            
            // Extract paths out to the length we need to connect the ends, or a bit further.
            // TODO: be sure to extract really big indels somehow...
            auto translator = algorithms::extract_connecting_graph(&synchronizer.graph,
                &context,
                (past_end - start) * 2,
                left_pos,
                right_pos,
                false); // We don't care about being strictly less than the specified length
                
#ifdef debug
            cerr << "Extracted " << context.graph.node_size() << " nodes and " << context.graph.edge_size() << " edges between " << path_name << ":" << start << "-" << past_end << endl;
#endif
                
            // Any ID mismatch is going to mess things up, since we need
            // operations on the new graph to make sense in the original graph.
            // We need all the entries in this translation map to be no-ops, so
            // we translate all IDs back to their original (which is possible
            // because we chose extraction paramters that never duplicate nodes).
            for (auto& kv : translator) {
                if (kv.first != kv.second) {
                    context.swap_node_id(kv.first, kv.second);
                }
            }
            
        } else {
            // We want to extract a radius
            
            // Find the center node, at the position we want to lock out from
            NodeSide center = synchronizer.get_path_index(path_name).at_position(path_offset);
            
            synchronizer.graph.nonoverlapping_node_context_without_paths(synchronizer.graph.get_node(center.node), context);
            synchronizer.graph.expand_context_by_length(context, context_bases, false, reflect);
        }
        
        // Also remember all the nodes connected to but not in the context,
        // which also need to be locked.
        periphery.clear();
        peripheral_attachments.clear();
        
        // We set this to false if a node we want is taken
        bool nodes_available = true;
        
        context.for_each_node([&](Node* node) {
            // For every node in the graph
            
            if (synchronizer.locked_nodes.count(node->id())) {
                // Someone else already has this node. So our condition is false
                // and we need to wait.
                nodes_available = false;
                return;
            }
            
            for (auto* edge : synchronizer.graph.edges_from(node)) {
                if (!context.has_node(edge->to())) {
                    // This is connected but not in the actual context graph. So it's on the periphery.
                    
                    if (synchronizer.locked_nodes.count(edge->to())) {
                        // Someone else already has this node. So our condition
                        // is false and we need to wait.
                        nodes_available = false;
                        return;
                    }
                    
                    // The destination of the edge is in the periphery
                    periphery.insert(edge->to());
                    // And you get to it from this side of this graph node.
                    peripheral_attachments[NodeSide(edge->from(), !edge->from_start())].insert(
                        NodeSide(edge->to(), edge->to_end()));
                }
            }
            for (auto* edge : synchronizer.graph.edges_to(node)) {
                if (!context.has_node(edge->from())) {
                    // This is connected but not in the actual context graph. So it's on the periphery.
                    
                    if (synchronizer.locked_nodes.count(edge->from())) {
                        // Someone else already has this node. So our condition
                        // is false and we need to wait.
                        nodes_available = false;
                        return;
                    }
                    
                    // The source of the edge is in the periphery
                    periphery.insert(edge->from());
                    // And you get to it from this side of this graph node.
                    peripheral_attachments[NodeSide(edge->to(), edge->to_end())].insert(
                        NodeSide(edge->from(), !edge->from_start()));
                }
            }
        });
        
        if (nodes_available) {
            // We can have the nodes we need. Remember what they are.
            subgraph = std::move(context);
        }
        
        // Return whether we were successful
        return nodes_available;
    });

    // Once we get here, we have a lock on the whole graph, and nobody else has
    // claimed our nodes. Our graph and periphery have been filled in, so we
    // just have to record our nodes as locked.
    
    for(id_t id : periphery) {
        // Mark the periphery
        synchronizer.locked_nodes.insert(id);
        locked_nodes.insert(id);
    }
    
    subgraph.for_each_node([&](Node* node) {
        // Mark the actual graph
        synchronizer.locked_nodes.insert(node->id());
        locked_nodes.insert(node->id());
    });
    
    // We should have actually grabbed something.
    if (locked_nodes.empty()) {
        cerr << "error:[vg::GraphSynchronizer] No nodes locked for " << path_name << ":" << start << "-" << past_end << endl;
        throw runtime_error("No nodes locked!");
    }
    
    // Now we know nobody else can touch those nodes and we can safely release
    // our lock on the main graph by letting it leave scope.
}

void GraphSynchronizer::Lock::unlock() {
    // Get the main graph lock
    std::unique_lock<std::mutex> lk(synchronizer.whole_graph_lock);
    
    // Release all the nodes
    for (id_t locked : locked_nodes) {
        synchronizer.locked_nodes.erase(locked);
    }
    
    // Clear our locked nodes
    locked_nodes.clear();
    
    // Notify anyone waiting, so they can all check to see if now they can go.
    lk.unlock();
    synchronizer.wait_for_region.notify_all();
}

VG& GraphSynchronizer::Lock::get_subgraph() {
    if (locked_nodes.empty()) {
        // Make sure we're actually locked
        throw runtime_error("No nodes are locked! Can't get graph!");
    }
    
    return subgraph;
}

pair<NodeSide, NodeSide> GraphSynchronizer::Lock::get_endpoints() const {
    if (locked_nodes.empty()) {
        // Make sure we're actually locked
        throw runtime_error("No nodes are locked! Can't get endpoints!");
    }
    
    if (start == 0 && past_end == 0) {
        // Make sure we used the endpoint-based constructor
        throw runtime_error("Graph was not locked with endpoints. Can't return andpoints!");
    }
    
    return endpoints;
}


set<NodeSide> GraphSynchronizer::Lock::get_peripheral_attachments(NodeSide graph_side) {
    if (peripheral_attachments.count(graph_side)) {
        // This side actually touches the periphery
        return peripheral_attachments.at(graph_side);
    } else {
        // This side doesn't touch anything outside the graph.
        return set<NodeSide>();
    }
}

vector<Translation> GraphSynchronizer::Lock::apply_edit(const Path& path, size_t max_node_size) {
    set<NodeSide> dangling;
    return apply_edit(path, dangling, max_node_size);
}

vector<Translation> GraphSynchronizer::Lock::apply_edit(const Path& path, set<NodeSide>& dangling, size_t max_node_size) {
    // Make sure we have exclusive ownership of the graph itself since we're
    // going to be modifying its data structures.
    std::lock_guard<std::mutex> guard(synchronizer.whole_graph_lock);
    
    for (size_t i = 0; i < path.mapping_size(); i++) {
        // Check each Mapping to make sure it's on a locked node
        auto node_id = path.mapping(i).position().node_id();
        if (!locked_nodes.count(node_id)) {
            throw runtime_error("Cannot edit unlocked node " + to_string(node_id));
        }
    }
    
    // Make all the edits, passing along the dangling node set.
    vector<Translation> translations = synchronizer.graph.edit_fast(path, dangling, max_node_size);
    
    // Lock all the nodes that result from the translations. They're guaranteed
    // to either be nodes we already have or novel nodes with fresh IDs.
    for (auto& translation : translations) {
        // For every translation's to path
        auto& new_path = translation.to();
        
        for (size_t i = 0; i < new_path.mapping_size(); i++) {
            // For every mapping to a node on that path
            auto node_id = new_path.mapping(i).position().node_id();
            
            if (!locked_nodes.count(i)) {
                // If it's not already locked, lock it.
                locked_nodes.insert(i);
                synchronizer.locked_nodes.insert(i);
            }
        }
    }
    
    // Apply the edits to the path indexes
    synchronizer.update_path_indexes(translations);
    
    // Spit out the translations to the caller. Maybe they can use them on their subgraph or something?
    return translations;
}

vector<Translation> GraphSynchronizer::Lock::apply_full_length_edit(const Path& path, size_t max_node_size) {
    // Find the left and right outer nodesides of the subgraph
    auto ends = get_endpoints();
    
    // Find everything attached to the left
    auto dangling = get_peripheral_attachments(ends.first);
    
    // Apply the edit, attaching its left end to the stuff attached to the left
    // end of the graph. Get back in the dangling set where the right end of the
    // edit's material is.
    auto translations = apply_edit(path, dangling, max_node_size);
    
    // Get the places that the right end of the graph attaches to
    auto right_periphery = get_peripheral_attachments(ends.second);
    
    // Get ownership of the graph because we're making edges
    std::lock_guard<std::mutex> guard(synchronizer.whole_graph_lock);
    
    for (const NodeSide& dangled : dangling) {
        // For every dangling NodeSide
        for (const NodeSide& attached : right_periphery) {
            // Attach it to each NodeSide the right end of the graph is attached to
            synchronizer.graph.create_edge(dangled, attached);
        }
    }
    
    return translations;
}

}





















