//
//  snarls.cpp
//
//

//#define debug

#include "snarls.hpp"
#include "json2pb.h"
#include "algorithms/topological_sort.hpp"
#include "algorithms/is_acyclic.hpp"

namespace vg {

CactusSnarlFinder::CactusSnarlFinder(VG& graph) :
    graph(graph) {
    // Make sure the graph is sorted.
    algorithms::topological_sort(&graph);
}

CactusSnarlFinder::CactusSnarlFinder(VG& graph, const string& hint_path) :
    CactusSnarlFinder(graph) {
    
    // Save the hint path
    hint_paths.insert(hint_path);
    
    // TODO: actually use it
}

SnarlManager CactusSnarlFinder::find_snarls() {
    
    if (graph.size() == 0) {
        // No snarls here!
        return SnarlManager();
    }
    // convert to cactus
    pair<stCactusGraph*, stList*> cac_pair = handle_graph_to_cactus(graph, hint_paths);
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
    
    // Finish the SnarlManager
    snarl_manager.finish();
    
    // Return the completed SnarlManager
    return snarl_manager;
    
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
            NetGraph connectivity_net_graph(start, end, child_chains, &graph, true);
            
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
            NetGraph flat_net_graph(start, end, child_chains, &graph);
            
            // This definitely should be calculated based on the internal-connectivity-ignoring net graph.
            snarl.set_directed_acyclic_net_graph(algorithms::is_directed_acyclic(&flat_net_graph));
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

bool start_backward(const Chain& chain) {
    // The start snarl is backward if it is marked backward.
    return !chain.empty() && chain.front().second;
}
    
bool end_backward(const Chain& chain) {
    // The end snarl is backward if it is marked backward.
    return !chain.empty() && chain.back().second;
}

Visit get_start_of(const Chain& chain) {
    // Get a bounding visit and return it.
    return chain.front().second ? reverse(chain.front().first->end()) : chain.front().first->start();
}
    
Visit get_end_of(const Chain& chain) {
    // Get a bounding visit and return it.
    return chain.back().second ? reverse(chain.back().first->start()) : chain.back().first->end();
}
    
bool ChainIterator::operator==(const ChainIterator& other) const {
    return (tie(go_left, pos, chain_start, chain_end, is_rend, complement) ==
            tie(other.go_left, other.pos, other.chain_start, other.chain_end, other.is_rend, other.complement));
}
    
bool ChainIterator::operator!=(const ChainIterator& other) const {
    auto unequal = !(*this == other);
    return unequal;
}
    
ChainIterator& ChainIterator::operator++() {
    if (go_left) {
        // Walk left
            
        if (pos == chain_start) {
            if (is_rend) {
                throw runtime_error("Walked off start!");
            }
                
            // We're already at the start, so next is just going to become rend
            is_rend = true;
        } else {
            // There's actually something to the left of us
            --pos;
        }
    } else {
        // Walk right
            
        if (pos == chain_end) {
            throw runtime_error("Walked off end!");
        }
            
        ++pos;
    }
        
    return *this;
}
    
pair<const Snarl*, bool> ChainIterator::operator*() const {
    return make_pair(pos->first, pos->second != complement);
}
    
const pair<const Snarl*, bool>* ChainIterator::operator->() const {
    // Make the pair we need
    scratch = *(*this);
    // Return a pointer to it.
    return &scratch;
}
    
ChainIterator chain_begin(const Chain& chain) {
    ChainIterator to_return{
        false, // Don't go left
        chain.begin(), // Be at the start of the chain
        chain.begin(), // Here's the chain's start
        chain.end(), // And its end
        false, // This is not a reverse end
        false // Do not complement snarl orientations
    };
        
    return to_return;
}
    
ChainIterator chain_end(const Chain& chain) {
    ChainIterator to_return{
        false, // Don't go left
        chain.end(), // Be at the end of the chain
        chain.begin(), // Here's the chain's start
        chain.end(), // And its end
        false, // This is not a reverse end
        false // Do not complement snarl orientations
    };
        
    return to_return;
}
    
ChainIterator chain_rbegin(const Chain& chain) {
    if (chain.empty()) {
        // If it's empty we should be the rend past-the-end reverse iterator
        return chain_rend(chain);
    }
        
    // Otherwise there's at least one element so point to the last.
    ChainIterator to_return{
        true, // Go left
        --chain.end(), // Be at the last real thing in the chain
        chain.begin(), // Here's the chain's start
        chain.end(), // And its end
        false, // This is not a reverse end
        false // Do not complement snarl orientations
    };
        
    return to_return;
}
    
ChainIterator chain_rend(const Chain& chain) {
    ChainIterator to_return{
        true, // Go left
        chain.begin(), // Be at the start of the chain
        chain.begin(), // Here's the chain's start
        chain.end(), // And its end
        true, // This is a reverse end
        false // Do not complement snarl orientations
    };
        
    return to_return;
}
    
ChainIterator chain_rcbegin(const Chain& chain) {
    ChainIterator to_return = chain_rbegin(chain);
    to_return.complement = true;
    return to_return;
}
    
ChainIterator chain_rcend(const Chain& chain) {
    ChainIterator to_return = chain_rend(chain);
    to_return.complement = true;
    return to_return;
}
    
ChainIterator chain_begin_from(const Chain& chain, const Snarl* start_snarl, bool snarl_orientation) {
    assert(!chain.empty());
    if (start_snarl == chain.front().first && snarl_orientation == start_backward(chain)) {
        // We are at the left end of the chain, in the correct orientation, so go forward
        return chain_begin(chain);
    } else if (start_snarl == chain.back().first) {
        // We are at the right end of the chain, so go reverse complement
        return chain_rcbegin(chain);
    } else {
        throw runtime_error("Tried to view a chain from a snarl not at either end!");
    }
    return ChainIterator();
}
    
ChainIterator chain_end_from(const Chain& chain, const Snarl* start_snarl, bool snarl_orientation) {
    assert(!chain.empty());
    if (start_snarl == chain.front().first && snarl_orientation == start_backward(chain)) {
        // We are at the left end of the chain, so go forward
        return chain_end(chain);
    } else if (start_snarl == chain.back().first) {
        // We are at the right end of the chain, so go reverse complement
        return chain_rcend(chain);
    } else {
        throw runtime_error("Tried to view a chain from a snarl not at either end!");
    }
    return ChainIterator();
} 

// TODO: this is duplicative with the other constructor, but protobuf won't let me make
// a deserialization iterator to match its signature because its internal file streams
// disallow copy constructors
SnarlManager::SnarlManager(istream& in) {
    // Add snarls to master list
    for (stream::ProtobufIterator<Snarl> iter(in); iter.has_next(); iter.get_next()) {
        // Add each snarl
        add_snarl(*iter);
    }
    // Record the tree structure and build the other indexes
    finish();
}
    
const vector<const Snarl*>& SnarlManager::children_of(const Snarl* snarl) const {
    if (snarl == nullptr) {
        // Looking for top level snarls
        return roots;
    }
    return record(snarl)->children;
}
    
const Snarl* SnarlManager::parent_of(const Snarl* snarl) const {
    return record(snarl)->parent;
}
    
const Snarl* SnarlManager::snarl_sharing_start(const Snarl* here) const {
    // Look out the start and see what we come to
    const Snarl* next = into_which_snarl(here->start().node_id(), !here->start().backward());
        
    // Return it, unless it's us, in which case we're a unary snarl that should go nowhere.
    return next == here ? nullptr : next;
        
}
    
const Snarl* SnarlManager::snarl_sharing_end(const Snarl* here) const {
    // Look out the end and see what we come to
    const Snarl* next = into_which_snarl(here->end().node_id(), here->end().backward());
        
    // Return it, unless it's us, in which case we're a unary snarl that should go nowhere.
    return next == here ? nullptr : next;
}
    
const Chain* SnarlManager::chain_of(const Snarl* snarl) const {
    return record(snarl)->parent_chain;
}

bool SnarlManager::chain_orientation_of(const Snarl* snarl) const { 
    const Chain* chain = chain_of(snarl);
    if (chain != nullptr) {
        // Go get the orientation flag.
        return chain->at(record(snarl)->parent_chain_index).second;
    }
    return false;
}

bool SnarlManager::in_nontrivial_chain(const Snarl* here) const {
    return chain_of(here)->size() > 1;
}
    
Visit SnarlManager::next_snarl(const Visit& here) const {
    // Must be a snarl visit
    assert(here.node_id() == 0);
    const Snarl* here_snarl = manage(here.snarl());
    assert(here_snarl != nullptr);
        
    // What visit are we going to return?
    Visit to_return;
        
    // And what snarl are we visiting next?
    const Snarl* next = here.backward() ? snarl_sharing_start(here_snarl) : snarl_sharing_end(here_snarl);

    if (next == nullptr) {
        // Nothing next
        return to_return;
    }
        
    // Fill in the start and end of the next snarl;
    transfer_boundary_info(*next, *to_return.mutable_snarl());
        
    if (here.backward()) {
        // We came out our start. So the next thing is also backward as long as its end matches our start. 
        to_return.set_backward(next->end().node_id() == here_snarl->start().node_id());
    } else {
        // We came out our end. So the next thing is backward if its start doesn't match our end.
        to_return.set_backward(next->start().node_id() != here_snarl->end().node_id());
    }
    
    return to_return;
}
    
Visit SnarlManager::prev_snarl(const Visit& here) const {
    return reverse(next_snarl(reverse(here)));
}
    
const deque<Chain>& SnarlManager::chains_of(const Snarl* snarl) const {
    if (snarl == nullptr) {
        // We want the root chains
        return root_chains;
    }
    
    // Otherwise, go look up the child chains of this snarl.
    return record(snarl)->child_chains;
}
    
NetGraph SnarlManager::net_graph_of(const Snarl* snarl, const HandleGraph* graph, bool use_internal_connectivity) const {
    // Just get the chains and forward them on to the NetGraph.
    // TODO: The NetGraph ends up computing its own indexes.
    return NetGraph(snarl->start(), snarl->end(), chains_of(snarl), graph, use_internal_connectivity);
}
    
bool SnarlManager::is_leaf(const Snarl* snarl) const {
    return record(snarl)->children.size() == 0;
}
    
bool SnarlManager::is_root(const Snarl* snarl) const {
    return parent_of(snarl) == nullptr;
}

bool SnarlManager::is_trivial(const Snarl* snarl, VG& graph) const {
    // If it's an ultrabubble with no children and no contained nodes, it is a trivial snarl.
    return snarl->type() == ULTRABUBBLE &&
        is_leaf(snarl)
        && shallow_contents(snarl, graph, false).first.size() == 0;
}

bool SnarlManager::all_children_trivial(const Snarl* snarl, VG& graph) const {
    for (auto& child : children_of(snarl)) {
        if (!is_trivial(child, graph)) {
            return false;
        }
    }
    return true;
}
    
const vector<const Snarl*>& SnarlManager::top_level_snarls() const {
    return roots;
}
    
void SnarlManager::for_each_top_level_snarl_parallel(const function<void(const Snarl*)>& lambda) const {
#pragma omp parallel for
    for (int i = 0; i < roots.size(); i++) {
        lambda(roots[i]);
    }
}
    
void SnarlManager::for_each_top_level_snarl(const function<void(const Snarl*)>& lambda) const {
    for (const Snarl* snarl : roots) {
        lambda(snarl);
    }
}
    
void SnarlManager::for_each_snarl_preorder(const function<void(const Snarl*)>& lambda) const {
    // We define a recursive function to apply the lambda in a preorder traversal of the snarl tree.
    std::function<void(const Snarl*)> process = [&](const Snarl* parent) {
        // Do the parent
        lambda(parent);
        for (auto child : children_of(parent)) {
            // Then do each child
            process(child);
        }
    };
        
    // Then we run that on the root of each tree of snarls.
    for_each_top_level_snarl(process);
}
    
void SnarlManager::for_each_snarl_parallel(const function<void(const Snarl*)>& lambda) const {
    // We define a recursive function to apply the lambda in a preorder traversal of the snarl tree.
    std::function<void(const Snarl*)> process = [&](const Snarl* parent) {
        // Do the parent
        lambda(parent);
            
        auto& children = children_of(parent);
            
#pragma omp parallel for
        for (size_t i = 0; i < children.size(); i++) {
            // Then do each child in parallel
            process(children[i]);
        }
    };
        
    for_each_top_level_snarl_parallel(process);
}

void SnarlManager::for_each_chain(const function<void(const Chain*)>& lambda) const {

    // We define a function to run a bunch of chains in serial
    auto do_chain_list = [&](const deque<Chain>& chains) {
        for (size_t i = 0; i < chains.size(); i++) {
            lambda(&chains[i]);
        }
    };
    
    // Do our top-level chains
    do_chain_list(root_chains);
    
    for_each_snarl_preorder([&](const Snarl* snarl) {
        // Then in preorder order through all the snarls, do the child chains.
        do_chain_list(chains_of(snarl));
    });
}

void SnarlManager::for_each_chain_parallel(const function<void(const Chain*)>& lambda) const {

    // We define a function to run a bunch of chains in parallel
    auto do_chain_list = [&](const deque<Chain>& chains) {
#pragma omp parallel for
        for (size_t i = 0; i < chains.size(); i++) {
            lambda(&chains[i]);
        }
    };
    
    // Do our top-level chains in parallel.
    do_chain_list(root_chains);
    
    for_each_snarl_parallel([&](const Snarl* snarl) {
        // Then in parallel through all the snarls, do the child chains in parallel.
        do_chain_list(chains_of(snarl));
    });
}

    
void SnarlManager::flip(const Snarl* snarl) {
        
    // Get a non-const pointer to the SnarlRecord, which we own.
    // Allowed because we ourselves aren't const.
    SnarlRecord* to_flip = (SnarlRecord*) record(snarl);
    // swap and reverse the start and end Visits
    int64_t start_id = to_flip->start().node_id();
    bool start_orientation = to_flip->start().backward();
        
    to_flip->mutable_start()->set_node_id(to_flip->end().node_id());
    to_flip->mutable_start()->set_backward(!to_flip->end().backward());
        
    to_flip->mutable_end()->set_node_id(start_id);
    to_flip->mutable_end()->set_backward(!start_orientation);
    
    if (to_flip->parent_chain != nullptr) {
        // Work out where we keep the orientation of this snarl in its parent chain
        bool& to_invert = (*to_flip->parent_chain)[to_flip->parent_chain_index].second;
        // And flip it
        to_invert = !to_invert;
    }
        
    // Note: snarl_into index is invariant to flipping.
    // All the other indexes live in the SnarlRecords and don't need to change.
}

void SnarlManager::flip(const Chain* chain) {

    if (chain->empty()) {
        // Empty chains are already flipped
        return;
    }

    // Get ahold of a non-const version of the chain, without casting.
    Chain* mutable_chain = record(chain_begin(*chain)->first)->parent_chain;

    // Bust open the chain abstraction and flip it.
    // First reverse the order
    std::reverse(mutable_chain->begin(), mutable_chain->end());
    for (auto& chain_entry : *mutable_chain) {
        // Flip all the orientation flags so now we know the snarls are the other way relative to their chain.
        chain_entry.second = !chain_entry.second;
        
        // Get a mutable snarl record
        SnarlRecord* mutable_snarl = (SnarlRecord*) record(chain_entry.first);
        
        // Flip around its index in its chain so it can find its record again.
        mutable_snarl->parent_chain_index = chain->size() - mutable_snarl->parent_chain_index - 1;
    }

}
    
const Snarl* SnarlManager::add_snarl(const Snarl& new_snarl) {
    // Don't let anyone add snarls if we are already finished.
    assert(!finished);

    // Allocate a default SnarlRecord
    snarls.emplace_back();
    
    SnarlRecord* new_record = &snarls.back();
    
    // Hackily copy the snarl in
    *new_record = new_snarl;
    
    // TODO: Should this be a non-default SnarlRecord constructor?

#ifdef debug
    cerr << "Adding snarl " << new_snarl.start().node_id() << " " << new_snarl.start().backward() << " -> "
         << new_snarl.end().node_id() << " " << new_snarl.end().backward() << endl;
#endif
        
    // We will set the parent and children and snarl_into and chain info when we finish().

    return unrecord(new_record);
}

void SnarlManager::finish() {
    // We can only do this once
    assert(!finished);
    
    // Mark ourselves finished so nobody can add more snarls.
    // Do it before indexing so we can use out partly built indexes without hitting asserts.
    finished = true;

    // Build all the indexes from the snarls we were given
    build_indexes();
    
    // Clean up the snarl and chain orientations so everything is predictably and intuitively oriented
    regularize();
}

const Snarl* SnarlManager::into_which_snarl(int64_t id, bool reverse) const {
    assert(finished);
    return snarl_into.count(make_pair(id, reverse)) ? snarl_into.at(make_pair(id, reverse)) : nullptr;
}
    
const Snarl* SnarlManager::into_which_snarl(const Visit& visit) const {
    return visit.has_snarl() ? manage(visit.snarl()) : into_which_snarl(visit.node_id(), visit.backward());
}
    
unordered_map<pair<int64_t, bool>, const Snarl*> SnarlManager::snarl_boundary_index() const {
    unordered_map<pair<int64_t, bool>, const Snarl*> index;
    for (const Snarl& snarl : snarls) {
        index[make_pair(snarl.start().node_id(), snarl.start().backward())] = &snarl;
        index[make_pair(snarl.end().node_id(), !snarl.end().backward())] = &snarl;
    }
    return index;
}
    
unordered_map<pair<int64_t, bool>, const Snarl*> SnarlManager::snarl_end_index() const {
    unordered_map<pair<int64_t, bool>, const Snarl*> index;
    for (const Snarl& snarl : snarls) {
        index[make_pair(snarl.end().node_id(), !snarl.end().backward())] = &snarl;
    }
    return index;
}
    
unordered_map<pair<int64_t, bool>, const Snarl*> SnarlManager::snarl_start_index() const {
    unordered_map<pair<int64_t, bool>, const Snarl*> index;
    for (const Snarl& snarl : snarls) {
        index[make_pair(snarl.start().node_id(), snarl.start().backward())] = &snarl;
    }
    return index;
}
    
void SnarlManager::build_indexes() {
#ifdef debug
    cerr << "Building SnarlManager index of " << snarls.size() << " snarls" << endl;
#endif

    // Reserve space for the snarl_into index, so we hopefully don't need to rehash or move anything.
    snarl_into.reserve(snarls.size() * 2);

    for (SnarlRecord& rec : snarls) {
        Snarl& snarl = *unrecord(&rec);
    
        // Build the snarl_into index first so we can manage() to resolve populated-snarl cross-references to parents later.
        snarl_into[make_pair(snarl.start().node_id(), snarl.start().backward())] = &snarl;
        snarl_into[make_pair(snarl.end().node_id(), !snarl.end().backward())] = &snarl;
#ifdef debug
        cerr << snarl.start().node_id() << " " << snarl.start().backward() << " reads into " << pb2json(snarl) << endl;
        cerr << snarl.end().node_id() << " " << !snarl.end().backward() << " reads into " << pb2json(snarl) << endl;
#endif
    }
    
        
    for (SnarlRecord& rec : snarls) {
        Snarl& snarl = *unrecord(&rec);
            
#ifdef debug
        cerr << pb2json(snarl) << endl;
#endif
            
        // is this a top-level snarl?
        if (snarl.has_parent()) {
            // add this snarl to the parent-to-children index
#ifdef debug
            cerr << "\tSnarl is a child" << endl;
#endif
            
            // Record it as a child of its parent
            SnarlRecord* parent = (SnarlRecord*) record(manage(snarl.parent()));
            parent->children.push_back(&snarl);
            
            // And that its parent is its parent
            record(&snarl)->parent = parent;
        }
        else {
            // record top level status
#ifdef debug
            cerr << "\tSnarl is top-level" << endl;
#endif
            roots.push_back(&snarl);
            
            record(&snarl)->parent = nullptr;
        }
    }
        
    // Compute the chains using the into and out-of indexes.
    
    // Compute the chains for the root level snarls
    root_chains = compute_chains(roots);
        
    // Build the back index from root snarl to containing chain
    for (Chain& chain : root_chains) {
        for (size_t i = 0; i < chain.size(); i++) {
            auto& oriented_snarl = chain[i];
            
            // Get the mutable record for each child snarl in the chain
            SnarlRecord* child_record = (SnarlRecord*) record(oriented_snarl.first);
            
            // Give it a pointer to the chain.
            child_record->parent_chain = &chain;
            // And where it is in it
            child_record->parent_chain_index = i;
        }
    }
    
    for (SnarlRecord& rec : snarls) {
        if (rec.children.empty()) {
            // Only look at snarls with children.
            continue;
        }
        
        // Compute the chains among the children
        rec.child_chains = compute_chains(rec.children);
        
        // Build the back index from child snarl to containing chain
        for (Chain& chain : rec.child_chains) {
            for (size_t i = 0; i < chain.size(); i++) {
                auto& oriented_snarl = chain[i];
                
                // Get the mutable record for each child snarl in the chain
                SnarlRecord* child_record = (SnarlRecord*) record(oriented_snarl.first);
                
                // Give it a pointer to the chain.
                child_record->parent_chain = &chain;
                // And where it is in it
                child_record->parent_chain_index = i;
            }
        }
    }
}
    
deque<Chain> SnarlManager::compute_chains(const vector<const Snarl*>& input_snarls) {
    // We populate this
    deque<Chain> to_return;
        
    // We track the snarls we have seen in chain traversals so we only have to see each chain once.
    unordered_set<const Snarl*> seen;
        
    for (const Snarl* snarl : input_snarls) {
        // For every snarl in this snarl (or, if snarl is null, every top level snarl)
            
        if (seen.count(snarl)) {
            // Already in a chain
            continue;
        }
            
        // Make a new chain for this child, with it in the forward direction in the chain.
        list<pair<const Snarl*, bool>> chain{{snarl, false}};
            
        // Mark it as seen
        seen.insert(snarl);
            
        // Make a visit to the child in forward orientation
        Visit here;
        transfer_boundary_info(*snarl, *here.mutable_snarl());
        // The default is already not-backward, but we set it anyway
        here.set_backward(false);
        
        for (Visit walk_left = prev_snarl(here);
             walk_left.has_snarl() && !seen.count(manage(walk_left.snarl()));
             walk_left = prev_snarl(walk_left)) {
            
            // For everything in the chain left from here, until we hit the
            // end or come back to the start
             
            // Add it to the chain in the orientation we find it
            chain.emplace_front(manage(walk_left.snarl()), walk_left.backward());
            // Mark it as seen
            seen.insert(chain.front().first);
        }
            
        for (Visit walk_right = next_snarl(here);
             walk_right.has_snarl() && !seen.count(manage(walk_right.snarl()));
             walk_right = next_snarl(walk_right)) {
                
            // For everything in the chain right from here, until we hit the
            // end or come back to the start
            
            // Add it to the chain in the orientation we find it
            chain.emplace_back(manage(walk_right.snarl()), walk_right.backward());
            // Mark it as seen
            seen.insert(chain.back().first);
        }
            
        // Copy from the list into a vector
        to_return.emplace_back(chain.begin(), chain.end());
    }
        
    return to_return;
}

void SnarlManager::regularize() {
    // Algorithm:
    // For each chain
    // Flip any snarls that are backward in the chain
    // If now the majority of the snarls end lower than they start, flip all the snarls and invert the chain.
    
#ifdef debug
    cerr << "Regularizing snarls and chains" << endl;
#endif
    
    for_each_chain_parallel([&](const Chain* chain) {
        // For every chain
        
        // Make a list of snarls to flip
        vector<const Snarl*> backward;
        // And a list of snarls to not flip
        vector<const Snarl*> forward;
        
        // Count the snarls that go low to high, as they should
        size_t correctly_oriented = 0;
        
        auto chain_start = chain_begin(*chain);
        auto chain_stop = chain_end(*chain);
        for (auto it = chain_start; it != chain_stop; ++it) {
            // For each snarl in the chain
            if (it->second) {
                // If it is backward, remember to flip it
                backward.push_back(it->first);
                
#ifdef debug
                cerr << "Snarl " << it->first->start() << " -> " << it->first->end() << " is backward in chain " << chain << endl;
#endif
                
                if (it->first->end().node_id() <= it->first->start().node_id()) {
                    // Count it as correctly oriented if it will be
#ifdef debug
                    cerr << "\tWill be graph-ascending when brought in line with chain" << endl;
#endif
                    correctly_oriented++;
                }
            } else {
                // If it is forward, remember that
                forward.push_back(it->first);
#ifdef debug
                cerr << "Snarl " << it->first->start() << " -> " << it->first->end() << " is forward in chain " << chain << endl;
#endif
                
                if (it->first->start().node_id() <= it->first->end().node_id()) {
                    // Count it as correctly oriented if it is
#ifdef debug
                    cerr << "\tIs graph-ascending already" << endl;
#endif
                    correctly_oriented++;
                }
            }
        }
        
#ifdef debug
        cerr << "Found " << correctly_oriented << "/" << chain->size() << " snarls of chain in graph-ascending orientation" << endl;
#endif
        
        if (correctly_oriented * 2 < chain->size()) {
            // Fewer than half the snarls are pointed the right way when they
            // go with the chain. (Don't divide chain size because then a chain
            // size of 1 requires 0 correctly oriented sanrls.)
            
#ifdef debug
            cerr << "Chain is in the worse orientation overall. Flipping!" << endl;
#endif
            
            // Really we want to invert the entire chain around the snarls, and
            // then only flip the formerly-chain-forward snarls.
            flip(chain);
            
            // Now set up to flip the other set of snarls
            backward.swap(forward);
            
        }
        
        for (auto& to_flip : backward) {
            // Flip all the snarls we found to flip to agree with the chain,
            // while not looping over the chain.
            flip(to_flip);
            
#ifdef debug
            cerr << "Flipped snarl to produce " << to_flip->start() << " " << to_flip->end() << endl;
#endif
        }
    });
    
}
    
pair<unordered_set<Node*>, unordered_set<Edge*> > SnarlManager::shallow_contents(const Snarl* snarl, VG& graph,
                                                                                 bool include_boundary_nodes) const {
        
    pair<unordered_set<Node*>, unordered_set<Edge*> > to_return;
        
    unordered_set<Node*> already_stacked;
        
    // initialize stack for DFS traversal of site
    vector<Node*> stack;
        
    Node* start_node = graph.get_node(snarl->start().node_id());
    Node* end_node = graph.get_node(snarl->end().node_id());
        
    // mark the boundary nodes as already stacked so that paths will terminate on them
    already_stacked.insert(start_node);
    already_stacked.insert(end_node);
        
    // add boundary nodes as directed
    if (include_boundary_nodes) {
        to_return.first.insert(start_node);
        to_return.first.insert(end_node);
    }
        
    vector<Edge*> edges_of_node;
        
    // stack up the nodes one edge inside the snarl from the start
    graph.edges_of_node(start_node, edges_of_node);
    for (Edge* edge : edges_of_node) {
            
        // does the edge point into the snarl?
        if (edge->from() == snarl->start().node_id() && edge->from_start() == snarl->start().backward()) {

            Node* node = graph.get_node(edge->to());
                
            if (!already_stacked.count(node)) {
                stack.push_back(node);
                already_stacked.insert(node);
            }

            to_return.second.insert(edge);
        }
        else if (edge->to() == snarl->start().node_id() && edge->to_end() != snarl->start().backward()) {

            Node* node = graph.get_node(edge->from());
                
            if (!already_stacked.count(node)) {
                stack.push_back(node);
                already_stacked.insert(node);
            }
                
            to_return.second.insert(edge);
        }
    }
    edges_of_node.clear();
        
    // stack up the nodes one edge inside the snarl from the end
    graph.edges_of_node(end_node, edges_of_node);
    for (Edge* edge : edges_of_node) {
        // does the edge point into the snarl?
        if (edge->from() == snarl->end().node_id() && edge->from_start() != snarl->end().backward()) {
                
            Node* node = graph.get_node(edge->to());
                
            if (!already_stacked.count(node)) {
                stack.push_back(node);
                already_stacked.insert(node);
            }
                
            to_return.second.insert(edge);
        }
        else if (edge->to() == snarl->end().node_id() && edge->to_end() == snarl->end().backward()) {
                
            Node* node = graph.get_node(edge->from());
                
            if (!already_stacked.count(node)) {
                stack.push_back(node);
                already_stacked.insert(node);
            }
                
            to_return.second.insert(edge);
        }
    }
    edges_of_node.clear();
        
    // traverse the snarl with DFS, skipping over any child snarls
    // do not pay attention to valid walks since we also want to discover any tips
    while (stack.size()) {
            
        // pop the top node off the stack
        Node* node = stack.back();
        stack.pop_back();
            
        // record that this node is in the snarl
        to_return.first.insert(node);
            
        const Snarl* forward_snarl = into_which_snarl(node->id(), false);
        const Snarl* backward_snarl = into_which_snarl(node->id(), true);
        if (forward_snarl) {
            // this node points into a snarl
                
            // What's on the other side of the snarl?
            id_t other_id = forward_snarl->start().node_id() == node->id() ? forward_snarl->end().node_id() : forward_snarl->start().node_id();
                
            // stack up the node on the opposite side of the snarl
            // rather than traversing it
            Node* opposite_node = graph.get_node(other_id);
            if (!already_stacked.count(opposite_node)) {
                stack.push_back(opposite_node);
                already_stacked.insert(opposite_node);
            }
        }
            
        if (backward_snarl) {
            // the reverse of this node points into a snarl
                
            // What's on the other side of the snarl?
            id_t other_id = backward_snarl->end().node_id() == node->id() ? backward_snarl->start().node_id(): backward_snarl->end().node_id();
                
            // stack up the node on the opposite side of the snarl
            // rather than traversing it
            Node* opposite_node = graph.get_node(other_id);
            if (!already_stacked.count(opposite_node)) {
                stack.push_back(opposite_node);
                already_stacked.insert(opposite_node);
            }
        }
            
        graph.edges_of_node(node, edges_of_node);
            
        for (Edge* edge : edges_of_node) {
            // which end of the edge is the current node?
            if (edge->from() == node->id()) {
                // does this edge point forward or backward?
                if ((edge->from_start() && !backward_snarl) ||
                    (!edge->from_start() && !forward_snarl)) {
                        
                    to_return.second.insert(edge);
                    Node* next_node = graph.get_node(edge->to());
                        
                    if (!already_stacked.count(next_node)) {
                            
                        stack.push_back(next_node);
                        already_stacked.insert(next_node);
                    }
                }
            }
            else {
                // does this edge point forward or backward?
                if ((edge->to_end() && !forward_snarl) ||
                    (!edge->to_end() && !backward_snarl)) {
                        
                    to_return.second.insert(edge);
                    Node* next_node = graph.get_node(edge->from());
                        
                    if (!already_stacked.count(next_node)) {
                            
                        stack.push_back(next_node);
                        already_stacked.insert(next_node);
                    }
                }
            }
        }
            
        edges_of_node.clear();
    }
        
    return to_return;
}
    
pair<unordered_set<Node*>, unordered_set<Edge*> > SnarlManager::deep_contents(const Snarl* snarl, VG& graph,
                                                                              bool include_boundary_nodes) const {
        
    pair<unordered_set<Node*>, unordered_set<Edge*> > to_return;
        
    unordered_set<Node*> already_stacked;
        
    // initialize stack for DFS traversal of site
    vector<Node*> stack;
        
    Node* start_node = graph.get_node(snarl->start().node_id());
    Node* end_node = graph.get_node(snarl->end().node_id());
        
    // mark the boundary nodes as already stacked so that paths will terminate on them
    already_stacked.insert(start_node);
    already_stacked.insert(end_node);
        
    // add boundary nodes as directed
    if (include_boundary_nodes) {
        to_return.first.insert(start_node);
        to_return.first.insert(end_node);
    }
        
    vector<Edge*> edges_of_node;
        
    // stack up the nodes one edge inside the snarl from the start
    graph.edges_of_node(start_node, edges_of_node);
    for (Edge* edge : edges_of_node) {
        // does the edge point into the snarl?
        if (edge->from() == snarl->start().node_id() && edge->from_start() == snarl->start().backward()) {
                
            Node* node = graph.get_node(edge->to());
                
            if (!already_stacked.count(node)) {
                stack.push_back(node);
                already_stacked.insert(node);
            }
                
            to_return.second.insert(edge);
        }
        else if (edge->to() == snarl->start().node_id() && edge->to_end() != snarl->start().backward()) {
                
            Node* node = graph.get_node(edge->from());
                
            if (!already_stacked.count(node)) {
                stack.push_back(node);
                already_stacked.insert(node);
            }
                
            to_return.second.insert(edge);
        }
    }
    edges_of_node.clear();
        
    // stack up the nodes one edge inside the snarl from the end
    graph.edges_of_node(end_node, edges_of_node);
    for (Edge* edge : edges_of_node) {
        // does the edge point into the snarl?
        if (edge->from() == snarl->end().node_id() && edge->from_start() != snarl->end().backward()) {
                
            Node* node = graph.get_node(edge->to());
                
            if (!already_stacked.count(node)) {
                stack.push_back(node);
                already_stacked.insert(node);
            }
                
            to_return.second.insert(edge);
        }
        else if (edge->to() == snarl->end().node_id() && edge->to_end() == snarl->end().backward()) {
                
            Node* node = graph.get_node(edge->from());
                
            if (!already_stacked.count(node)) {
                stack.push_back(node);
                already_stacked.insert(node);
            }
                
            to_return.second.insert(edge);
        }
    }
    edges_of_node.clear();
        
    // traverse the snarl with DFS, skipping over any child snarls
    // do not pay attention to valid walks since we also want to discover any tips
    while (stack.size()) {
            
        // pop the top node off the stack
        Node* node = stack.back();
        stack.pop_back();
            
        // record that this node is in the snarl
        to_return.first.insert(node);
            
        graph.edges_of_node(node, edges_of_node);
            
        for (Edge* edge : edges_of_node) {
            to_return.second.insert(edge);
            // get the other end of the edge
            Node* next_node = edge->from() == node->id() ? graph.get_node(edge->to()) :
                graph.get_node(edge->from());
            if (!already_stacked.count(next_node)) {
                stack.push_back(next_node);
                already_stacked.insert(next_node);
            }
        }
            
        edges_of_node.clear();
    }
        
    return to_return;
}
    
const Snarl* SnarlManager::manage(const Snarl& not_owned) const {
    // TODO: keep the Snarls in some kind of sorted order to make lookup
    // efficient. We could also have a map<Snarl, Snarl*> but that would be
    // a tremendous waste of space.
    
    // Work out how to read into the snarl
    pair<int64_t, bool> reading_in = make_pair(not_owned.start().node_id(), not_owned.start().backward());
    
    // Get the cannonical pointer to the snarl that we are reading into with the start, inward visit.
    auto it = snarl_into.find(reading_in);
        
    if (it == snarl_into.end()) {
        // It's not there. Someone is trying to manage a snarl we don't
        // really own. Complain.
        throw runtime_error("Unable to find snarl " +  pb2json(not_owned) + " in SnarlManager");
    }
        
    // Return the official copy of that snarl
    return it->second;
}
    
vector<Visit> SnarlManager::visits_right(const Visit& visit, VG& graph, const Snarl* in_snarl) const {
        
#ifdef debug
    cerr << "Look right from " << visit << endl;
#endif
        
    // We'll populate this
    vector<Visit> to_return;
        
    // Find the right side of the visit we're on
    NodeSide right_side = to_right_side(visit);
        
    if (visit.node_id() == 0) {
        // We're leaving a child snarl, so we are going to need to check if
        // another child snarl shares this boundary node in the direction
        // we're going.
            
        const Snarl* child = into_which_snarl(right_side.node, !right_side.is_end);
        if (child != nullptr && child != in_snarl
            && into_which_snarl(right_side.node, right_side.is_end) != in_snarl) {
            // We leave the one child and immediately enter another!
                
            // Make a visit to it
            Visit child_visit;
            transfer_boundary_info(*child, *child_visit.mutable_snarl());
                
            if (right_side.node == child->end().node_id()) {
                // We came in its end
                child_visit.set_backward(true);
            } else {
                // We should have come in its start
                assert(right_side.node == child->start().node_id());
            }
                
            // Bail right now, so we don't try to explore inside this child snarl.
            to_return.push_back(child_visit);
            return to_return;
                
        }
            
    }
        
    for (auto attached : graph.sides_of(right_side)) {
        // For every NodeSide attached to the right side of this visit
            
#ifdef debug
        cerr << "\tFind NodeSide " << attached << endl;
#endif
            
        const Snarl* child = into_which_snarl(attached.node, attached.is_end);
        if (child != nullptr && child != in_snarl
            && into_which_snarl(attached.node, !attached.is_end) != in_snarl) {
            // We're reading into a child
                
#ifdef debug
            cerr << "\t\tGoes to Snarl " << *child << endl;
#endif
                
            if (attached.node == child->start().node_id()) {
                // We're reading into the start of the child
                    
                // Make a visit to the child snarl
                Visit child_visit;
                transfer_boundary_info(*child, *child_visit.mutable_snarl());
                    
#ifdef debug
                cerr << "\t\tProduces Visit " << child_visit << endl;
#endif
                    
                // Put it in in the forward orientation
                to_return.push_back(child_visit);
            } else if (attached.node == child->end().node_id()) {
                // We're reading into the end of the child
                    
                // Make a visit to the child snarl
                Visit child_visit;
                transfer_boundary_info(*child, *child_visit.mutable_snarl());
                child_visit.set_backward(true);
                    
#ifdef debug
                cerr << "\t\tProduces Visit " << child_visit << endl;
#endif
                    
                // Put it in in the reverse orientation
                to_return.push_back(child_visit);
            } else {
                // Should never happen
                throw runtime_error("Read into child " + pb2json(*child) + " with non-matching traversal");
            }
        } else {
            // We just go into a normal node
            to_return.emplace_back();
            Visit& next_visit = to_return.back();
            next_visit.set_node_id(attached.node);
            next_visit.set_backward(attached.is_end);
            
#ifdef debug
            cerr << "\t\tProduces Visit " << to_return.back() << endl;
#endif
                
        }
    }
        
    return to_return;
        
}
    
vector<Visit> SnarlManager::visits_left(const Visit& visit, VG& graph, const Snarl* in_snarl) const {
        
    // Get everything right of the reversed visit
    vector<Visit> to_return = visits_right(reverse(visit), graph, in_snarl);
        
    // Un-reverse them so they are in the correct orientation to be seen
    // left of here.
    for (auto& v : to_return) {
        v = reverse(v);
    }
        
    return to_return;
        
}
    
NetGraph::NetGraph(const Visit& start, const Visit& end, const HandleGraph* graph, bool use_internal_connectivity) :
    graph(graph),
    start(graph->get_handle(start.node_id(), start.backward())),
    end(graph->get_handle(end.node_id(), end.backward())),
    use_internal_connectivity(use_internal_connectivity) {
    // Nothing to do!
}
    
NetGraph::NetGraph(const Visit& start, const Visit& end,
                   const vector<vector<pair<Snarl, bool>>>& child_chains,
                   const vector<Snarl>& child_unary_snarls,
                   const HandleGraph* graph,
                   bool use_internal_connectivity) :
                   NetGraph(start, end, graph, use_internal_connectivity) {
        
    for (auto& unary : child_unary_snarls) {
        add_unary_child(&unary);
    }
        
    for (auto& chain : child_chains) {
        Chain converted_chain;
        for (auto& item : chain) {
            // Convert from actual snarls to pointers
            converted_chain.emplace_back(&item.first, item.second);
        }
        add_chain_child(converted_chain);
    }
   
}
    
void NetGraph::add_unary_child(const Snarl* unary) {
    // For each unary snarl, make its bounding handle
    handle_t snarl_bound = graph->get_handle(unary->start().node_id(), unary->start().backward());
        
    // Get its ID
    id_t snarl_id = unary->start().node_id();
        
    // Make sure it is properly specified to be unary (in and out the same node in opposite directions)
    assert(unary->end().node_id() == snarl_id);
    assert(unary->end().backward() == !unary->start().backward());
        
    // Save it as a unary snarl
    unary_boundaries.insert(snarl_bound);
        
    if (use_internal_connectivity) {
        // Save its connectivity
        connectivity[snarl_id] = make_tuple(unary->start_self_reachable(), unary->end_self_reachable(),
                                            unary->start_end_reachable());
    } else {
        // Use the connectivity of an ordinary node that has a different
        // other side. Don't set connected_start_end because, for a real
        // unary snarl, the end and the start are the same, so that
        // means you can turn aroiund.
        connectivity[snarl_id] = make_tuple(false, false, false);
    }
}
    
void NetGraph::add_chain_child(const Chain& chain) {
    // For every chain, get its bounding handles in the base graph
    handle_t chain_start_handle = graph->get_handle(get_start_of(chain));
    handle_t chain_end_handle = graph->get_handle(get_end_of(chain));
        
    // Save the links that let us cross the chain.
    chain_ends_by_start[chain_start_handle] = chain_end_handle;
    chain_end_rewrites[graph->flip(chain_end_handle)] = graph->flip(chain_start_handle);
        
    if (use_internal_connectivity) {
        
        // Determine child snarl connectivity.
        bool connected_left_left = false;
        bool connected_right_right = false;
        bool connected_left_right = true;
            
        for (auto it = chain_begin(chain); it != chain_end(chain); ++it) {
            // Go through the oriented child snarls from left to right
            const Snarl* child = it->first;
            bool backward = it->second;
                
            // Unpack the child's connectivity
            bool start_self_reachable = child->start_self_reachable();
            bool end_self_reachable = child->end_self_reachable();
            bool start_end_reachable = child->start_end_reachable();
                
            if (backward) {
                // Look at the connectivity in reverse
                std::swap(start_self_reachable, end_self_reachable);
            }
                
            if (start_self_reachable) {
                // We found a turnaround from the left
                connected_left_left = true;
            }
                
            if (!start_end_reachable) {
                // There's an impediment to getting through.
                connected_left_right = false;
                // Don't keep looking for turnarounds
                break;
            }
        }
            
        for (auto it = chain_rbegin(chain); it != chain_rend(chain); ++it) {
            // Go through the oriented child snarls from left to right
            const Snarl* child = it->first;
            bool backward = it->second;
                
            // Unpack the child's connectivity
            bool start_self_reachable = child->start_self_reachable();
            bool end_self_reachable = child->end_self_reachable();
            bool start_end_reachable = child->start_end_reachable();
                
            if (backward) {
                // Look at the connectivity in reverse
                std::swap(start_self_reachable, end_self_reachable);
            }
                
            if (end_self_reachable) {
                // We found a turnaround from the right
                connected_right_right = true;
                break;
            }
                
            if (!start_end_reachable) {
                // Don't keep looking for turnarounds
                break;
            }
        }
            
        // Save the connectivity
        connectivity[graph->get_id(chain_start_handle)] = tie(connected_left_left,
                                                              connected_right_right, connected_left_right);
    } else {
        // Act like a normal connected-through node.
        connectivity[graph->get_id(chain_start_handle)] = make_tuple(false, false, true);
    }
}
    
handle_t NetGraph::get_handle(const id_t& node_id, bool is_reverse) const {
    // We never let anyone see any node IDs that aren't assigned to child snarls/chains or content nodes.
    return graph->get_handle(node_id, is_reverse);
}
    
id_t NetGraph::get_id(const handle_t& handle) const {
    // We just use the handle/ID mapping of the backing graph
    return graph->get_id(handle);
}
    
bool NetGraph::get_is_reverse(const handle_t& handle) const {
    // We just use the handle/orientation mapping of the backing graph
    return graph->get_is_reverse(handle);
}
    
handle_t NetGraph::flip(const handle_t& handle) const {
    // We just use the flip implementation of the backing graph
    return graph->flip(handle);
}
    
size_t NetGraph::get_length(const handle_t& handle) const {
    // TODO: We don't really want to support this operation; maybe lengths
    // and sequences should be factored out into another interface.
    throw runtime_error("Cannot expose sequence lengths via NetGraph");
}
    
string NetGraph::get_sequence(const handle_t& handle) const {
    // TODO: We don't really want to support this operation; maybe lengths
    // and sequences should be factored out into another interface.
    throw runtime_error("Cannot expose sequences via NetGraph");
}
    
bool NetGraph::follow_edges(const handle_t& handle, bool go_left, const function<bool(const handle_t&)>& iteratee) const {
    // Now we do the real work.
        
#ifdef debug
    cerr << "Look for edges on " << graph->get_id(handle) << " " << graph->get_is_reverse(handle)
         << " going " << (go_left ? "left" : "right") << endl;
#endif
        
    // We also need to deduplicate edges. Maybe the start and end of a chain
    // connect to the same next node, and we could read out both traversing
    // the chain.
    unordered_set<handle_t> seen;
            
    // This deduplicates and emits the edges, and also handles rewriting
    // visits to the ends of chains as visits to the start, which we use to
    // represent the whole chain.
    auto handle_edge = [&](const handle_t& other) -> bool {
            
        handle_t real_handle = other;
        if (chain_end_rewrites.count(other)) {
            // We're reading into the end of a chain.
            // Warp to the start.
            real_handle = chain_end_rewrites.at(other);
        } else if (chain_end_rewrites.count(graph->flip(other))) {
            // We're backing into the end of a chain.
            // Warp to the start.
            real_handle = graph->flip(chain_end_rewrites.at(graph->flip(other)));
        }
            
#ifdef debug
        cerr << "Found edge " << (go_left ? "from " : "to ") << graph->get_id(other) << " " << graph->get_is_reverse(other) << endl;
#endif
            
        if (!seen.count(real_handle)) {
            seen.insert(real_handle);
#ifdef debug
            cerr << "Report as " << graph->get_id(real_handle) << " " << graph->get_is_reverse(real_handle) << endl;
#endif
                
            return iteratee(real_handle);
        } else {
#ifdef debug
            cerr << "Edge has been seen" << endl;
#endif
            return true;
        }
    };
        
    // This does the same as handle_edge, but flips the real handle
    auto flip_and_handle_edge = [&](const handle_t& other) -> bool {
            
        handle_t real_handle = other;
        if (chain_end_rewrites.count(other)) {
            // We're reading into the end of a chain.
            // Warp to the start.
            real_handle = chain_end_rewrites.at(other);
        } else if (chain_end_rewrites.count(graph->flip(other))) {
            // We're backing into the end of a chain.
            // Warp to the start.
            real_handle = graph->flip(chain_end_rewrites.at(graph->flip(other)));
        }
            
        real_handle = graph->flip(real_handle);
            
#ifdef debug
        cerr << "Found edge " << (go_left ? "from " : "to ") << graph->get_id(other) << " " << graph->get_is_reverse(other) << endl;
#endif
            
        if (!seen.count(real_handle)) {
            seen.insert(real_handle);
#ifdef debug
            cerr << "Report as " << graph->get_id(real_handle) << " " << graph->get_is_reverse(real_handle) << endl;
#endif
                
            return iteratee(real_handle);
        } else {
#ifdef debug
            cerr << "Edge has been seen" << endl;
#endif
            return true;
        }
    };
        
    // Each way of doing things needs to support going either left or right
        
    if (end != start &&
      ((handle == end && !go_left) || (handle == graph->flip(end) && go_left) ||
        (handle == graph->flip(start) && !go_left) || (handle == start && go_left))) {
        // If we're looking outside of the snarl this is the net graph for, don't admit to having any edges.
        //If start and end are the same, all edges are within the net graph
            
#ifdef debug
        cerr << "We are at the bound of the graph so don't say anything" << endl;
#endif
        return true;
    }
        
    if (chain_ends_by_start.count(handle) || chain_ends_by_start.count(graph->flip(handle))) {
        // If we have an associated chain end for this start, we have to use chain connectivity to decide what to do.
            
#ifdef debug
        cerr << "We are a chain node" << endl;
#endif
            
        bool connected_start_start;
        bool connected_end_end;
        bool connected_start_end;
        tie(connected_start_start, connected_end_end, connected_start_end) = connectivity.at(graph->get_id(handle));
            
#ifdef debug
        cerr << "Connectivity: " << connected_start_start << " " << connected_end_end << " " << connected_start_end << endl;
#endif
            
        if (chain_ends_by_start.count(handle)) {
            // We visit the chain in its forward orientation
                
#ifdef debug
            cerr << "We are visiting the chain forward" << endl;
#endif
                
            if (go_left) {
                // We want predecessors.
                // So we care about end-end connectivity (how could we have left our end?)
                    
#ifdef debug
                cerr << "We are going left from a forward chain" << endl;
#endif
                    
                if (connected_end_end) {
                    
#ifdef debug
                    cerr << "We can reverse and go back out the end" << endl;
#endif
                    
                    // Anything after us but in its reverse orientation could be our predecessor
                    // But a thing after us may be a chain, in which case we need to find its head before flipping.
                    if (!graph->follow_edges(chain_ends_by_start.at(handle), false, flip_and_handle_edge)) {
                        // Iteratee is done
                        return false;
                    }
                }
                    
                if (connected_start_end) {
                    
#ifdef debug
                    cerr << "We can continue through and go out the start" << endl;
#endif
                    
                    // Look left out of the start of the chain (which is the handle we really are on)
                    if (!graph->follow_edges(handle, true, handle_edge)) {
                        // Iteratee is done
                        return false;
                    }
                }
                    
            } else {
                // We want our successors
                    
#ifdef debug
                cerr << "We are going right from a forward chain" << endl;
#endif
                    
                if (connected_start_start) {
                    
#ifdef debug
                    cerr << "We can reverse and go back out the start" << endl;
#endif
                    
                    // Anything before us but in its reverse orientation could be our successor
                    // But a thing before us may be a chain, in which case we need to find its head before flipping.
                    if (!graph->follow_edges(handle, true, flip_and_handle_edge)) {
                        // Iteratee is done
                        return false;
                    }
                }
                    
                if (connected_start_end) {
                    
#ifdef debug
                    cerr << "We can continue through and go out the end" << endl;
#endif
                    
                    // Look right out of the end of the chain (which is the handle we really are on)
                    if (!graph->follow_edges(chain_ends_by_start.at(handle), false, handle_edge)) {
                        // Iteratee is done
                        return false;
                    }
                }
                    
            }
                
        } else {
            // We visit the chain in its reverse orientation.
            // Just flip the cases of above and reverse all the emitted orientations.
                
#ifdef debug
            cerr << "We are visiting the chain in reverse" << endl;
#endif

            if (go_left) {
                // We want predecessors of the reverse version (successors, but flipped)
                    
#ifdef debug
                cerr << "We are going left from a reverse chain" << endl;
#endif
                    
                if (connected_start_start) {
                    
#ifdef debug
                    cerr << "We can reverse and go back out the start" << endl;
#endif
                    
                    if (!graph->follow_edges(handle, false, flip_and_handle_edge)) {
                        // Iteratee is done
                        return false;
                    }
                }
                    
                if (connected_start_end) {
                    
#ifdef debug
                    cerr << "We can continue through and go out the end" << endl;
#endif
                    
                    if (!graph->follow_edges(chain_ends_by_start.at(graph->flip(handle)), false, flip_and_handle_edge)) {
                        // Iteratee is done
                        return false;
                    }
                }
                    
            } else {
                // We want successors of the reverse version (predecessors, but flipped)
                    
#ifdef debug
                cerr << "We are going right from a reverse chain" << endl;
#endif
                    
                if (connected_end_end) {
                    
#ifdef debug
                    cerr << "We can reverse and go back out the end" << endl;
#endif
                    
                    if (!graph->follow_edges(chain_ends_by_start.at(graph->flip(handle)), false, handle_edge)) {
                        // Iteratee is done
                        return false;
                    }
                }
                    
                if (connected_start_end) {
                    
#ifdef debug
                    cerr << "We can continue through and go out the start" << endl;
#endif
                    
                    if (!graph->follow_edges(handle, false, handle_edge)) {
                        // Iteratee is done
                        return false;
                    }
                }
                    
            }
                
        }
            
        return true;
    }
        
    if (unary_boundaries.count(handle) || unary_boundaries.count(graph->flip(handle))) {
        // We are dealign with a node representing a unary child snarl.
            
#ifdef debug
        cerr << "We are looking at a unary snarl" << endl;
#endif
            
        // We have to use chain connectivity to decide what to do.
        bool connected_start_start;
        bool connected_end_end;
        bool connected_start_end;
        tie(connected_start_start, connected_end_end, connected_start_end) = connectivity.at(graph->get_id(handle));
            
        if (unary_boundaries.count(handle)) {
            // We point into a unary snarl
            if (go_left) {
                // We want the predecessors
                    
                if (!use_internal_connectivity) {
                    // We treat this as a normal node
                    // Get the real predecessors
                    if (!graph->follow_edges(handle, true, handle_edge)) {
                        // Iteratee is done
                        return false;
                    }
                }
                // Otherwise just think about what we can traverse to (i.e. nothing)
                    
                // Can't read a forward oriented unary boundary as a
                // predecessor, so no need to support read through.
                    
            } else {
                // We want the successors
                    
                // No real successors
                    
                if (connected_start_start || connected_end_end || connected_start_end) {
                    // Successors also include our predecessors but backward
                    if (!graph->follow_edges(handle, true, flip_and_handle_edge)) {
                        // Iteratee is done
                        return false;
                    }
                        
                }
            }
        } else {
            // We point out of a unary snarl.
            // Reverse of above. Sort of.
            if (go_left) {
                if (connected_start_start || connected_end_end || connected_start_end) {
                    if (!graph->follow_edges(handle, false, flip_and_handle_edge)) {
                        // Iteratee is done
                        return false;
                    }
                        
                }
                    
            } else {
                if (!use_internal_connectivity) {

                    if (!graph->follow_edges(handle, false, handle_edge)) {
                        // Iteratee is done
                        return false;
                    }
                }
                
            }
            
        }
            
        return true;
    }
        
#ifdef debug
    cerr << "We are an ordinary node" << endl;
#endif
        
    // Otherwise, this is an ordinary snarl content node
    return graph->follow_edges(handle, go_left, handle_edge);
}
    
void NetGraph::for_each_handle(const function<bool(const handle_t&)>& iteratee, bool parallel) const {
    // Find all the handles by a traversal.
        
    // We have to do the traversal on the underlying backing graph, because
    // the traversal functions we implemented on the graph we present will
    // maybe use internal child snarl connectivity, which can mean parts of
    // the graph are in this snarl but not actually reachable.
        
    // We let both the starts and ends of child snarls into the queue.
    // But we'll only reveal the starts to our iteratee.
    list<handle_t> queue;
    unordered_set<id_t> queued;
        
    // Start at both the start and the end of the snarl.
    queue.push_back(start);
    queued.insert(graph->get_id(start));
    queue.push_back(end);
    queued.insert(graph->get_id(end));
        
    while (!queue.empty()) {
        handle_t here = queue.front();
        queue.pop_front();
            
        if (unary_boundaries.count(graph->flip(here))) {
            // This is a backward unary child snarl head, so we need to look at it the other way around.
            here = graph->flip(here);
        } else if (chain_ends_by_start.count(graph->flip(here))) {
            // This is a backward child chain head, so we need to look at it the other way around.
            here = graph->flip(here);
        } else if (chain_end_rewrites.count(graph->flip(here))) {
            // This is a backward child chain tail, so we need to look at it the other way around.
            here = graph->flip(here);
        }
            
        if (!chain_end_rewrites.count(here)) {
            // This is not a chain end, so it's either a real contained node or a chain head.
            // We can emit it.
                
            if (graph->get_is_reverse(here)) {
                if (!iteratee(graph->flip(here))) {
                    // Run the iteratee on the forward version, and stop if it wants to stop
                    break;
                }
            } else {
                if (!iteratee(here)) {
                    // Run the iteratee, and stop if it wants to stop.
                    break;
                }
            }
                
        } 
            
        // We define a function to queue up nodes we could visit next
        auto handle_edge = [&](const handle_t& other) {
            // Whenever we see a new node, add it to the queue
            if (!queued.count(graph->get_id(other))) {
                queue.push_back(other);
                queued.insert(graph->get_id(other));
            }
        };
            
        // We already have flipped any backward heads or tails frontward. So
        // we don't need to check if the backward version of us is in
        // anything.
            
        if ( ((start != end && here != end && here != graph->flip(start)) ||
               start == end)
            && !unary_boundaries.count(here) &&
            !chain_ends_by_start.count(here)  && !chain_end_rewrites.count(here)) {
                
            // We have normal graph to our right and not the exterior of this snarl or the interior of a child.
            graph->follow_edges(here, false, handle_edge);
        }
            
        if ((start != end && here != start && here != graph->flip(end)) ||
             start == end) {
            // We have normal graph to our left.
            graph->follow_edges(here, true, handle_edge);
        }
            
        if (chain_end_rewrites.count(here)) {
            // We need to look right off the reverse head of this child snarl.
            graph->follow_edges(chain_end_rewrites.at(here), false, handle_edge);
        }
            
        if (chain_ends_by_start.count(here)) {
            // We need to look right off the (reverse) tail of this child snarl.
            graph->follow_edges(chain_ends_by_start.at(here), false, handle_edge);
        }
    }
}
    
size_t NetGraph::node_size() const {
    // TODO: this is inefficient!
    size_t size = 0;
    for_each_handle([&](const handle_t& ignored) {
            size++;
        });
    return size;
}

id_t NetGraph::min_node_id() const {
    // TODO: this is inefficient!
    id_t winner = numeric_limits<id_t>::max();
    for_each_handle([&](const handle_t& handle) {
            winner = min(winner, this->get_id(handle));
        });
    return winner;
}

id_t NetGraph::max_node_id() const {
    // TODO: this is inefficient!
    id_t winner = numeric_limits<id_t>::min();
    for_each_handle([&](const handle_t& handle) {
            winner = max(winner, this->get_id(handle));
        });
    return winner;
}
    
const handle_t& NetGraph::get_start() const {
    return start;
}
    
const handle_t& NetGraph::get_end() const {
    return end;
}
    
bool NetGraph::is_child(const handle_t& handle) const {
    // It's a child if we're going forward or backward through a chain, or into a unary snarl.
    
#ifdef debug
    cerr << "Is " << graph->get_id(handle) << " " << graph->get_is_reverse(handle) << " a child?" << endl;
    
    for (auto& kv : chain_ends_by_start) {
        cerr << "\t" << graph->get_id(kv.first) << " " << graph->get_is_reverse(kv.first) << " is a child." << endl;
    }
    
    for (auto& kv : chain_ends_by_start) {
        auto flipped = graph->flip(kv.first);
        cerr << "\t" << graph->get_id(flipped) << " " << graph->get_is_reverse(flipped) << " is a child." << endl;
    }
    
    for (auto& boundary : unary_boundaries) {
        cerr << "\t" << graph->get_id(boundary) << " " << graph->get_is_reverse(boundary) << " is a child." << endl;
    }
#endif
    
    return chain_ends_by_start.count(handle) ||
        chain_ends_by_start.count(flip(handle)) ||
        unary_boundaries.count(handle);
}

handle_t NetGraph::get_inward_backing_handle(const handle_t& child_handle) const {
    if (chain_ends_by_start.count(child_handle)) {
        // Reading into a chain, so just return this
        return child_handle;
    } else if (chain_ends_by_start.count(flip(child_handle))) {
        // Reading out of a chain, so get the outward end of the chain and flip it
        return graph->flip(chain_ends_by_start.at(flip(child_handle)));
    } else if (unary_boundaries.count(child_handle)) {
        // Reading into a unary snarl.
        // Always already facing inward.
        // TODO: what if we're reading out of a chain *and* into a unary snarl?
        return child_handle;
    } else {
        throw runtime_error("Cannot get backing handle for a handle that is not a handle to a child's node in the net graph");
    }
}

handle_t NetGraph::get_handle_from_inward_backing_handle(const handle_t& backing_handle) const {
    if (chain_ends_by_start.count(backing_handle)) {
        // If it's a recorded chain start it gets passed through
        return backing_handle;
    } else if (chain_end_rewrites.count(backing_handle)) {
        // If it's a known chain end, we produce the start in reverse orientation, which we stored.
        return chain_end_rewrites.at(backing_handle);
    } else if (unary_boundaries.count(backing_handle)) {
        // Unary snarl handles are passed through too.
        return backing_handle;
    } else {
        throw runtime_error("Cannot assign backing handle to a child chain or unary snarl");
    }
}
    
bool operator==(const Visit& a, const Visit& b) {
    // IDs and orientations have to match, and nobody has a snarl or the
    // snarls match.
    return a.node_id() == b.node_id() &&
        a.backward() == b.backward() &&
        ((!a.has_snarl() && !b.has_snarl()) ||
         a.snarl() == b.snarl());
}
    
bool operator!=(const Visit& a, const Visit& b) {
    return !(a == b);
}
    
bool operator<(const Visit& a, const Visit& b) {
    if (!a.has_snarl() && !b.has_snarl()) {
        // Compare everything but the snarl
        return make_tuple(a.node_id(), a.backward()) < make_tuple(b.node_id(), b.backward());
    } else {
        // Compare including the snarl
        return make_tuple(a.node_id(), a.snarl(), a.backward()) < make_tuple(b.node_id(), b.snarl(), b.backward());
    }        
}
    
ostream& operator<<(ostream& out, const Visit& visit) {
    if (!visit.has_snarl()) {
        // Use the node ID
        out << visit.node_id();
    } else {
        // Use the snarl
        out << visit.snarl();
    }
    out << " " << (visit.backward() ? "rev" : "fwd");
    return out;
}
    
bool operator==(const SnarlTraversal& a, const SnarlTraversal& b) {
    if (a.visit_size() != b.visit_size()) {
        return false;
    }
    for (size_t i = 0; i < a.visit_size(); i++) {
        if (a.visit(i) != b.visit(i)) {
            return false;
        }
    }
    // Otherwise everything we can think of matches
    return true;
}
    
bool operator!=(const SnarlTraversal& a, const SnarlTraversal& b) {
    return !(a == b);
}
    
bool operator<(const SnarlTraversal& a, const SnarlTraversal& b) {
    for (size_t i = 0; i < b.visit_size(); i++) {
        if (i >= a.visit_size()) {
            // A has run out and B is still going
            return true;
        }
            
        if (a.visit(i) < b.visit(i)) {
            return true;
        } else if (b.visit(i) < a.visit(i)) {
            return false;
        }
    }
        
    // If we get here either they're equal or A has more visits than B
    return false;
}
    
bool operator==(const Snarl& a, const Snarl& b) {
    if (a.type() != b.type()) {
        return false;
    }
    if (a.start() != b.start()) {
        return false;
    }
    if (a.end() != b.end()) {
        return false;
    }
    if (a.has_parent() || b.has_parent()) {
        // Someone has a parent so we must compare them.
        return a.parent() == b.parent();
    }
    return true;
}
    
bool operator!=(const Snarl& a, const Snarl& b) {
    return !(a == b);
}
    
bool operator<(const Snarl& a, const Snarl& b) {
    if (!a.has_parent() && !b.has_parent()) {
        // Compare without parent
        return make_tuple(a.type(), a.start(), a.end()) < make_tuple(b.type(), b.start(), b.end());
    } else {
        // Compare with parent
        return make_tuple(a.type(), a.start(), a.end(), a.parent()) < make_tuple(b.type(), b.start(), b.end(), b.parent());
    }
}
    
ostream& operator<<(ostream& out, const Snarl& snarl) {
    return out << snarl.start() << "-" << snarl.end();
}
}






