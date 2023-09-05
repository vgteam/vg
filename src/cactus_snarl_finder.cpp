///
///  \file cactus_snarl_finder.cpp
///
///

//#define debug

#include "subgraph_overlay.hpp"
#include "handle.hpp"

#include "cactus_snarl_finder.hpp"

namespace vg {

using namespace std;

CactusSnarlFinder::CactusSnarlFinder(const PathHandleGraph& graph, const string& hint_path) :
    graph(&graph) {
    if (!hint_path.empty()) {
        hint_paths.insert(hint_path);
        // TODO: actually use it
    }
}

SnarlManager CactusSnarlFinder::find_snarls_impl(bool known_single_component, bool finish_index) {
    
    if (graph->get_node_count() <= 1) {
        // No snarls here!
        return SnarlManager();
    }
    // convert to cactus
    pair<stCactusGraph*, stList*> cac_pair = handle_graph_to_cactus(*graph, hint_paths, known_single_component);
    stCactusGraph* cactus_graph = cac_pair.first;
    stList* telomeres = cac_pair.second;

    // get the snarl decomposition as a C struct
    stSnarlDecomposition *snarls = stCactusGraph_getSnarlDecomposition(cactus_graph, telomeres);
    
    // Get a non-owning pointer to the list of chains (which are themselves lists of snarls).
    stList* cactus_chains_list = snarls->topLevelChains;
    
    // And one to the list of top-level unary snarls
    stList* cactus_unary_snarls_list = snarls->topLevelUnarySnarls;
    
    
    // We'll fill this with all the snarls
    SnarlManager snarl_manager;
    
    // Fill the manager with all of the snarls, recursively.
    recursively_emit_snarls(Visit(), Visit(), Visit(), Visit(), cactus_chains_list, cactus_unary_snarls_list, snarl_manager);
    
    // Free the decomposition
    stSnarlDecomposition_destruct(snarls);
    
    // Free the telomeres
    stList_destruct(telomeres);

    // free the cactus graph
    stCactusGraph_destruct(cactus_graph);

    if (finish_index) {
        // Finish the SnarlManager
        snarl_manager.finish();
    }
    
    // Return the completed SnarlManager
    return snarl_manager;
    
}

SnarlManager CactusSnarlFinder::find_snarls() {
    return find_snarls_impl(false, true);
}

SnarlManager CactusSnarlFinder::find_snarls_parallel() {

    vector<unordered_set<id_t>> weak_components = handlealgs::weakly_connected_components(graph);
    vector<SnarlManager> snarl_managers(weak_components.size());

#pragma omp parallel for schedule(dynamic, 1)
    for (size_t i = 0; i < weak_components.size(); ++i) {
        const PathHandleGraph* subgraph;
        if (weak_components.size() == 1) {
            subgraph = graph;
        } else {
            // turn the component into a graph
            subgraph = new PathSubgraphOverlay(graph, &weak_components[i]);
        }
        string hint_path = !hint_paths.empty() ? *hint_paths.begin() : "";
        CactusSnarlFinder finder(*subgraph, hint_path);
        // find the snarls, telling the finder that the graph is a single component
        // and that we don't want to finish the snarl index
        snarl_managers[i] = finder.find_snarls_impl(true, false);
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


const Snarl* CactusSnarlFinder::recursively_emit_snarls(const Visit& start, const Visit& end,
                                                        const Visit& parent_start, const Visit& parent_end,
                                                        stList* chains_list, stList* unary_snarls_list, SnarlManager& destination) {
        
#ifdef debug    
    cerr << "Explore snarl " << start << " -> " << end << endl;
#endif
           
    // This is the snarl we are filling in to add to the SnarlManger, or an
    // empty snarl if we're a fake root snarl.
    Snarl snarl;
        
    if (start.node_id() != 0 && end.node_id() != 0) {
        // This is a real snarl
                
        // Set up the start and end
        *snarl.mutable_start() = start;
        *snarl.mutable_end() = end;
        
        if (parent_start.node_id() != 0 && parent_end.node_id() != 0) {
            // We have a parent that isn't the fake root, so fill in its ends
            *snarl.mutable_parent()->mutable_start() = parent_start;
            *snarl.mutable_parent()->mutable_end() = parent_end;
        }
    } 
    
    // This will hold the pointer to the copy of the snarl in the SnarlManager,
    // or null if the snarl is a fake root and we don't add it.
    const Snarl* managed = nullptr;
    
    // Before we can pass our snarl to the snarl manager, we need to look at all
    // its children so we can get connectivity info.
    
    // We have a vector of the snarls made for the child snarls in each ordinary
    // chain, plus trivial chains for the unary snarls.
    vector<Chain> child_chains;
    
#ifdef debug
    cerr << "Look at " << stList_length(chains_list) << " child chains" << endl;
#endif
    
    int chain_offset = 0;
    for (int64_t i = 0; i < stList_length(chains_list); i++) {
        // For each child chain
        stList* cactus_chain = (stList*)stList_get(chains_list, i);
            
        // Make a new chain.
        // We aren't going to pass it on to the snarl manager, because chains need to be recomputed for consistency.
        // But we need it for computing the internal snarl connectivity.
        child_chains.emplace_back();
        auto& chain = child_chains.back();
        
#ifdef debug
        cerr << "Chain " << i << " has " << stList_length(cactus_chain) << " child snarls" << endl;
#endif
        
        for (int64_t j = 0; j < stList_length(cactus_chain); j++) {
            // for each child snarl in the chain
            stSnarl* child_snarl = (stSnarl*)stList_get(cactus_chain, j);

            // scrape the vg coordinate information out of the cactus ends where we stuck
            // it during cactus construction
            CactusSide* cac_child_side1 = (CactusSide*)stCactusEdgeEnd_getObject(child_snarl->edgeEnd1);
            CactusSide* cac_child_side2 = (CactusSide*)stCactusEdgeEnd_getObject(child_snarl->edgeEnd2);
            
            // Convert from CactusSide (the interior endpoint of each node) to Visit (inward at start, outward at end)
            Visit child_start;
            child_start.set_node_id(cac_child_side1->node);
            // Start is backward if the interior is not an end
            child_start.set_backward(!cac_child_side1->is_end);
            Visit child_end;
            child_end.set_node_id(cac_child_side2->node);
            // End is backward if the interior is an end
            child_end.set_backward(cac_child_side2->is_end);
                
            // Recursively create a snarl for the child
            const Snarl* converted_child = recursively_emit_snarls(child_start, child_end, start, end,
                                                             child_snarl->chains, child_snarl->unarySnarls, destination);
            // Work out if it should be backward in the chain
            bool backward_in_chain = false;
            if (!chain.empty()) {
                 bool last_backward_in_chain = chain.back().second;
                 auto dangling_id = last_backward_in_chain ? chain.back().first->end().node_id() : chain.back().first->start().node_id();
                 // We are backward if our end is shared with the previous snarl in the chain.
                 backward_in_chain = converted_child->end().node_id() == dangling_id;
            }
            
            // And then add it to this chain.
            chain.emplace_back(converted_child, backward_in_chain);
        }
    }
    
#ifdef debug
    cerr << "Look at " << stList_length(unary_snarls_list) << " child unary snarls" << endl;
#endif
    
    for (int64_t i = 0; i < stList_length(unary_snarls_list); i++) {
        // for each child unary snarl
        stSnarl* child_snarl = (stSnarl*)stList_get(unary_snarls_list, i);

        // TODO: deduplicate this code

        // scrape the vg coordinate information out of the cactus ends where we stuck
        // it during cactus construction
        CactusSide* cac_child_side1 = (CactusSide*)stCactusEdgeEnd_getObject(child_snarl->edgeEnd1);
        CactusSide* cac_child_side2 = (CactusSide*)stCactusEdgeEnd_getObject(child_snarl->edgeEnd2);
        
        // Convert from CactusSide (the interior endpoint of each node) to Visit (inward at start, outward at end)
        Visit child_start;
        child_start.set_node_id(cac_child_side1->node);
        // Start is backward if the interior is not an end
        child_start.set_backward(!cac_child_side1->is_end);
        Visit child_end;
        child_end.set_node_id(cac_child_side2->node);
        // End is backward if the interior is an end
        child_end.set_backward(cac_child_side2->is_end);
        
        // Make a trivial chain
        child_chains.emplace_back();
        auto& chain = child_chains.back();
        
        // Recursively create a snarl for the child, and then add it to the trivial chain as forward
        chain.emplace_back(recursively_emit_snarls(child_start, child_end, start, end,
                                                   child_snarl->chains, child_snarl->unarySnarls, destination), false);
    }

    if (snarl.start().node_id() != 0 || snarl.end().node_id() != 0) {
        // This snarl is real, we care about type and connectivity.

        // First determine connectivity
        {

            // Make a net graph for the snarl that uses internal connectivity
            NetGraph connectivity_net_graph(start, end, child_chains, graph, true);
            
            // Evaluate connectivity
            // A snarl is minimal, so we know out start and end will be normal nodes.
            handle_t start_handle = connectivity_net_graph.get_handle(start.node_id(), start.backward());
            handle_t end_handle = connectivity_net_graph.get_handle(end.node_id(), end.backward());
            
            // Start out by assuming we aren't connected
            bool connected_start_start = false;
            bool connected_end_end = false;
            bool connected_start_end = false;
            
            // We do a couple of direcred walk searches to test connectivity.
            list<handle_t> queue{start_handle};
            unordered_set<handle_t> queued{start_handle};
            auto handle_edge = [&](const handle_t& other) {
#ifdef debug
                cerr << "\tCan reach " << connectivity_net_graph.get_id(other)
                << " " << connectivity_net_graph.get_is_reverse(other) << endl;
#endif
                
                // Whenever we see a new node orientation, queue it.
                if (!queued.count(other)) {
                    queue.push_back(other);
                    queued.insert(other);
                }
            };
            
#ifdef debug
            cerr << "Looking for start-start turnarounds and through connections from "
                 << connectivity_net_graph.get_id(start_handle) << " " <<
                connectivity_net_graph.get_is_reverse(start_handle) << endl;
#endif
            
            while (!queue.empty()) {
                handle_t here = queue.front();
                queue.pop_front();
                
                if (here == end_handle) {
                    // Start can reach the end
                    connected_start_end = true;
                }
                
                if (here == connectivity_net_graph.flip(start_handle)) {
                    // Start can reach itself the other way around
                    connected_start_start = true;
                }
                
                if (connected_start_end && connected_start_start) {
                    // No more searching needed
                    break;
                }
                
                // Look at everything reachable on a proper rightward directed walk.
                connectivity_net_graph.follow_edges(here, false, handle_edge);
            }
            
            auto end_inward = connectivity_net_graph.flip(end_handle);
            
#ifdef debug
            cerr << "Looking for end-end turnarounds from " << connectivity_net_graph.get_id(end_inward)
                 << " " << connectivity_net_graph.get_is_reverse(end_inward) << endl;
#endif
            
            // Reset and search the other way from the end to see if it can find itself.
            queue = {end_inward};
            queued = {end_inward};
            while (!queue.empty()) {
                handle_t here = queue.front();
                queue.pop_front();
                
#ifdef debug
                cerr << "Got to " << connectivity_net_graph.get_id(here) << " "
                     << connectivity_net_graph.get_is_reverse(here) << endl;
#endif
                
                if (here == end_handle) {
                    // End can reach itself the other way around
                    connected_end_end = true;
                    break;
                }
                
                // Look at everything reachable on a proper rightward directed walk.
                connectivity_net_graph.follow_edges(here, false, handle_edge);
            }
            
            // Save the connectivity info. TODO: should the connectivity flags be
            // calculated based on just the net graph, or based on actual connectivity
            // within child snarls.
            snarl.set_start_self_reachable(connected_start_start);
            snarl.set_end_self_reachable(connected_end_end);
            snarl.set_start_end_reachable(connected_start_end);

#ifdef debug
            cerr << "Connectivity: " << connected_start_start << " " << connected_end_end << " " << connected_start_end << endl;
#endif
            
        
        }
        
        {
            // Determine cyclicity/acyclicity
        
            // Make a net graph that just pretends child snarls/chains are ordinary nodes
            NetGraph flat_net_graph(start, end, child_chains, graph);
            
            // This definitely should be calculated based on the internal-connectivity-ignoring net graph.
            snarl.set_directed_acyclic_net_graph(handlealgs::is_directed_acyclic(&flat_net_graph));
        }

        // Now we need to work out if the snarl can be a unary snarl or an ultrabubble or what.
        if (start.node_id() == end.node_id()) {
            // Snarl has the same start and end (or no start or end, in which case we don't care).
            snarl.set_type(UNARY);
#ifdef debug
            cerr << "Snarl is UNARY" << endl;
#endif
        } else if (!snarl.start_end_reachable()) {
            // Can't be an ultrabubble if we're not connected through.
            snarl.set_type(UNCLASSIFIED);
#ifdef debug
            cerr << "Snarl is UNCLASSIFIED because it doesn't connect through" << endl;
#endif
        } else if (snarl.start_self_reachable() || snarl.end_self_reachable()) {
            // Can't be an ultrabubble if we have these cycles
            snarl.set_type(UNCLASSIFIED);
            
#ifdef debug
            cerr << "Snarl is UNCLASSIFIED because it allows turning around, creating a directed cycle" << endl;
#endif

        } else {
            // See if we have all ultrabubble children
            bool all_ultrabubble_children = true;
            for (auto& chain : child_chains) {
                for (auto& child : chain) {
                    if (child.first->type() != ULTRABUBBLE) {
                        all_ultrabubble_children = false;
                        break;
                    }
                }
                if (!all_ultrabubble_children) {
                    break;
                }
            }
            
            // Note that ultrabubbles *can* loop back on their start or end.
            
            if (!all_ultrabubble_children) {
                // If we have non-ultrabubble children, we can't be an ultrabubble.
                snarl.set_type(UNCLASSIFIED);
#ifdef debug
                cerr << "Snarl is UNCLASSIFIED because it has non-ultrabubble children" << endl;
#endif
            } else if (!snarl.directed_acyclic_net_graph()) {
                // If all our children are ultrabubbles but we ourselves are cyclic, we can't be an ultrabubble
                snarl.set_type(UNCLASSIFIED);
                
#ifdef debug
                cerr << "Snarl is UNCLASSIFIED because it is not directed-acyclic" << endl;
#endif
            } else {
                // We have only ultrabubble children and are acyclic.
                // We're an ultrabubble.
                snarl.set_type(ULTRABUBBLE);
#ifdef debug
                cerr << "Snarl is an ULTRABUBBLE" << endl;
#endif
            }
        }
        
        // Now we know enough about the snarl to actually put it in the SnarlManager
        managed = destination.add_snarl(snarl);
        
    }
    
    // Return a pointer to the managed snarl.
    return managed;
}

}
