///
///  \file snarls.cpp
///
///

// #define debug

#include <vg/io/protobuf_emitter.hpp>

#include "snarls.hpp"
#include "vg/io/json2pb.h"
#include "subgraph_overlay.hpp"

namespace vg {

SnarlManager SnarlFinder::find_snarls_parallel() {
    // By default, just use a single thread, unless this finder has a parallel
    // overriding implementation.
    return find_snarls();
}

HandleGraphSnarlFinder::HandleGraphSnarlFinder(const HandleGraph* graph) : graph(graph) {
    // Nothing to do!
}

SnarlManager HandleGraphSnarlFinder::find_snarls_unindexed() {
    // Start with an empty SnarlManager
    SnarlManager snarl_manager;
    
    // We need a stack with the information we need to translate the traversal
    // into vg::Snarl and vg::Chain objects, so we can compute connectivity and
    // snarl classification as we go up.
    struct TranslationFrame {
        // This will hold the unmanaged scratch snarl we pass to the manager.
        Snarl snarl;
        // This will hold all the child snarls that need their parent information filled in before they can become managed.
        // They are sorted by chain.
        vector<vector<Snarl>> child_chains;
        // For creating the current chain for this frame, we need to know where the chain claimed to start.
        // If the start = the end and the chain is inside a snarl, it's just a trivial chain (single node) and we drop it.
        handle_t current_chain_start;
    };
    
    // Stack that lets us connect snarls to their parents.
    // Holds each snarl and the child snarls we have finished for it so far.
    vector<TranslationFrame> stack;
   
    traverse_decomposition([&](handle_t chain_start) {
        // We got the start of a (possibly empty) chain.
        if (!stack.empty()) {
            // We're in a snarl, so we're a chain that we need for snarl connectivity/classification.
            stack.back().current_chain_start = chain_start;
            
            // Allocate a place to store the snarls in the chain.
            stack.back().child_chains.emplace_back();
        }
    }, [&](handle_t chain_end) {
        // We got the end of a (possibly empty) chain.
        if (!stack.empty() && stack.back().current_chain_start == chain_end) {
            // We're an empty chain in an actual snarl.
            // Get rid of our empty chain vector that got no snarls in it
            assert(stack.back().child_chains.back().empty());
            stack.back().child_chains.pop_back();
        }
    }, [&](handle_t snarl_start) {
        // Stack up a snarl
        stack.emplace_back();
        // And fill in its start
        auto& snarl = stack.back().snarl;
        snarl.mutable_start()->set_node_id(graph->get_id(snarl_start));
        snarl.mutable_start()->set_backward(graph->get_is_reverse(snarl_start));
    }, [&](handle_t snarl_end) {
        // Fill in its end
        auto& snarl = stack.back().snarl;
        snarl.mutable_end()->set_node_id(graph->get_id(snarl_end));
        snarl.mutable_end()->set_backward(graph->get_is_reverse(snarl_end));
        
        // We need to manage all our children and put them in Chain objects that net graphs can understand.
        vector<Chain> managed_child_chains;
        
        for (auto& child_chain : stack.back().child_chains) {
            // For every child chain
            
            // Make a translated version
            managed_child_chains.emplace_back();
            for (auto& child : child_chain) {
                // For each child snarl, fill us in as the parent (before we have connectivity info filled in)
                *child.mutable_parent() = snarl;
                // And report it to the manager with the cross-reference to us filled in.
                const Snarl* managed_child = snarl_manager.add_snarl(child);
                // And save it in the child chain.
                // We know it must be forward in the chain.
                managed_child_chains.back().emplace_back(managed_child, false);
            }
        }
        
        // This snarl is real, we care about type and connectivity.
        // All its children are done.

        /////
        // Determine connectivity
        /////
        
        // Make a net graph for the snarl that uses internal connectivity
        NetGraph connectivity_net_graph(snarl.start(), snarl.end(), managed_child_chains, graph, true);
        
        // Evaluate connectivity
        // A snarl is minimal, so we know out start and end will be normal nodes.
        handle_t start_handle = connectivity_net_graph.get_handle(snarl.start().node_id(), snarl.start().backward());
        handle_t end_handle = connectivity_net_graph.get_handle(snarl.end().node_id(), snarl.end().backward());
        
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

        /////
        // Determine tip presence
        /////
        
        // Make a net graph that just pretends child snarls/chains are ordinary nodes
        NetGraph flat_net_graph(snarl.start(), snarl.end(), managed_child_chains, graph);
        
        // Having internal tips in the net graph disqualifies a snarl from being an ultrabubble
        auto tips = handlealgs::find_tips(&flat_net_graph);

#ifdef debug
        cerr << "Tips: " << endl;
        for (auto& tip : tips) {
            cerr << "\t" << flat_net_graph.get_id(tip) << (flat_net_graph.get_is_reverse(tip) ? '-' : '+') << endl;
        }
#endif

        // We should have at least the bounding nodes.
        assert(tips.size() >= 2);
        bool has_internal_tips = (tips.size() > 2); 
        
        /////
        // Determine cyclicity/acyclicity
        /////
    
        // This definitely should be calculated based on the internal-connectivity-ignoring net graph.
        snarl.set_directed_acyclic_net_graph(handlealgs::is_directed_acyclic(&flat_net_graph));

        /////
        // Determine classification
        /////

        // Now we need to work out if the snarl can be a unary snarl or an ultrabubble or what.
        if (snarl.start().node_id() == snarl.end().node_id()) {
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
            for (auto& chain : managed_child_chains) {
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
            
            if (!all_ultrabubble_children) {
                // If we have non-ultrabubble children, we can't be an ultrabubble.
                snarl.set_type(UNCLASSIFIED);
#ifdef debug
                cerr << "Snarl is UNCLASSIFIED because it has non-ultrabubble children" << endl;
#endif
            } else if (has_internal_tips) {
                // If we have internal tips, we can't be an ultrabubble
                snarl.set_type(UNCLASSIFIED);
                
#ifdef debug
                cerr << "Snarl is UNCLASSIFIED because it contains internal tips" << endl;
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
        
        // Now we know all about our snarl, but we don't know about our parent.
        
        if (stack.size() > 1) {
            // We have a parent. Join it as a child, at the end of the current chain
            assert(!stack[stack.size() - 2].child_chains.empty());
            stack[stack.size() - 2].child_chains.back().emplace_back(std::move(snarl));
        } else {
            // Just manage ourselves now, because our parent can't manage us.
            snarl_manager.add_snarl(snarl);
        }
        
        // Leave the stack
        stack.pop_back();
    });
    
    // Give it back
    return snarl_manager;
}

SnarlManager HandleGraphSnarlFinder::find_snarls() {
    // Find all the snarls
    auto snarl_manager(find_snarls_unindexed());
    
    // Index them
    snarl_manager.finish();
    
    // Return the finished SnarlManager
    return snarl_manager;
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

SnarlManager::SnarlManager(istream& in) : SnarlManager([&in](const function<void(Snarl&)>& consume_snarl) -> void {
    // Find all the snarls in the input stream and use each of them in the callback-based constructor
    for (vg::io::ProtobufIterator<Snarl> iter(in); iter.has_current(); iter.advance()) {
        consume_snarl(*iter);
    }
}) {
    // Nothing to do!
}

SnarlManager::SnarlManager(const function<void(const function<void(Snarl&)>&)>& for_each_snarl) {
    for_each_snarl([&](Snarl& snarl) {
        // Add each snarl to us
        add_snarl(snarl);
    });
    // Record the tree structure and build the other indexes
    finish();
}

void SnarlManager::serialize(ostream& out) const {
    
    vg::io::ProtobufEmitter<Snarl> emitter(out);
    list<const Snarl*> stack;

    for (const Snarl* root : top_level_snarls()) {
        stack.push_back(root);
        
        while (!stack.empty()) {
            // Grab a snarl from the stack
            const Snarl* snarl = stack.back();
            stack.pop_back();
            
            // Write out the snarl
            emitter.write_copy(*root);

            for (const Snarl* child_snarl : children_of(snarl)) {
                // Stack up its children
                stack.push_back(child_snarl);
            }
        }
    }
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

size_t SnarlManager::chain_rank_of(const Snarl* snarl) const { 
    const Chain* chain = chain_of(snarl);
    if (chain != nullptr) {
        // The index is a perfectly good rank.
        return record(snarl)->parent_chain_index;
    }
    // If you're in a single-snarl chain you are at index 0.
    return 0;
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

bool SnarlManager::is_trivial(const Snarl* snarl, const HandleGraph& graph) const {
    // If it's an ultrabubble with no children and no contained nodes, it is a trivial snarl.
    return snarl->type() == ULTRABUBBLE &&
        is_leaf(snarl)
        && shallow_contents(snarl, graph, false).first.size() == 0;
}

bool SnarlManager::all_children_trivial(const Snarl* snarl, const HandleGraph& graph) const {
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
    #pragma omp parallel
    {
        #pragma omp single
        {
            for (int i = 0; i < roots.size(); i++) {
                #pragma omp task firstprivate(i)
                {
                    lambda(roots[i]);
                }
            }
        }
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

void SnarlManager::for_each_top_level_chain(const function<void(const Chain*)>& lambda) const {
    for (const Chain& chain : root_chains) {
        lambda(&chain);
    }    
}

void SnarlManager::for_each_top_level_chain_parallel(const function<void(const Chain*)>& lambda) const {
#pragma omp parallel for schedule(dynamic, 1)
    for (size_t i = 0; i < root_chains.size(); ++i) {
        lambda(&root_chains[i]);
    }
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

void SnarlManager::for_each_snarl_unindexed(const function<void(const Snarl*)>& lambda) const {
    for (const SnarlRecord& snarl_record : snarls) {
        lambda(unrecord(&snarl_record));
    }
}

const Snarl* SnarlManager::discrete_uniform_sample(minstd_rand0& random_engine)const{
    // have to set the seed to the random engine in the unit tests , pass the random engine 

    int number_of_snarls = num_snarls();
#ifdef debug
    cerr << "number_of_snarls "<< number_of_snarls <<endl;
    for (int i =0; i< snarls.size(); i++){
        const Snarl* snarl  = unrecord(&snarls[i]);
        cerr << snarl->start().node_id() << " -> "<<snarl->end().node_id() <<endl;
    }
#endif

    // if we have no snarls we return a flag 
    if(number_of_snarls ==0){
        return nullptr;
    }
    
    // we choose a snarl from the master list of snarls in the graph at random uniformly
    // unif[a,b],  deque starts at index 0 so upperbound is size-1
    uniform_int_distribution<int> distribution(0, number_of_snarls-1);  
    int random_num = distribution(random_engine);
#ifdef debug
    cerr << "modifying snarl num " << random_num << endl;  
    if(unrecord(&snarls[random_num]) == nullptr){
        cerr << "unrecorded snarl is null" <<endl;
    }else{
       const Snarl* snarl  = unrecord(&snarls[random_num]);
       cerr << snarl->start() << endl;
       cerr << snarl->end() <<endl;
    }
#endif

    return unrecord(&snarls[random_num]);

} 

int SnarlManager::num_snarls()const{
    // get size of snarls in the master list deque<SnarlRecord> snarls 
    int num_snarls = this->snarls.size();
    return num_snarls;

}

    
void SnarlManager::flip(const Snarl* snarl) {
        
    // Get a non-const pointer to the SnarlRecord, which we own.
    // Allowed because we ourselves aren't const.
    SnarlRecord* to_flip = (SnarlRecord*) record(snarl);
    // Get the Snarl of it
    Snarl* to_flip_snarl = unrecord(to_flip);
    // swap and reverse the start and end Visits
    int64_t start_id = to_flip_snarl->start().node_id();
    bool start_orientation = to_flip_snarl->start().backward();
        
    to_flip_snarl->mutable_start()->set_node_id(to_flip_snarl->end().node_id());
    to_flip_snarl->mutable_start()->set_backward(!to_flip_snarl->end().backward());
        
    to_flip_snarl->mutable_end()->set_node_id(start_id);
    to_flip_snarl->mutable_end()->set_backward(!start_orientation);
    
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

    // Allocate a default SnarlRecord
    snarls.emplace_back();
    
    SnarlRecord* new_record = &snarls.back();
    
    // Hackily copy the snarl in
    *new_record = new_snarl;

    // Initialized snarl number for each record as deque is being filled
    new_record->snarl_number = (size_t)snarls.size()-1;
    
    // TODO: Should this be a non-default SnarlRecord constructor?

#ifdef debug
    cerr << "Adding snarl " << new_snarl.start().node_id() << " " << new_snarl.start().backward() << " -> "
         << new_snarl.end().node_id() << " " << new_snarl.end().backward() << endl;
#endif
        
    // We will set the parent and children and snarl_into and chain info when we finish().

    return unrecord(new_record);
}

void SnarlManager::finish() {
    // Build all the indexes from the snarls we were given
    build_indexes();
    
    // Clean up the snarl and chain orientations so everything is predictably and intuitively oriented
    regularize();

}

const Snarl* SnarlManager::into_which_snarl(int64_t id, bool reverse) const {
    return snarl_into.count(make_pair(id, reverse)) ? snarl_into.at(make_pair(id, reverse)) : nullptr;
}
    
const Snarl* SnarlManager::into_which_snarl(const Visit& visit) const {
    return visit.has_snarl() ? manage(visit.snarl()) : into_which_snarl(visit.node_id(), visit.backward());
}
    
unordered_map<pair<int64_t, bool>, const Snarl*> SnarlManager::snarl_boundary_index() const {
    unordered_map<pair<int64_t, bool>, const Snarl*> index;
    for (const SnarlRecord& snarl_record : snarls) {
        const Snarl& snarl = *unrecord(&snarl_record);
        index[make_pair(snarl.start().node_id(), snarl.start().backward())] = &snarl;
        index[make_pair(snarl.end().node_id(), !snarl.end().backward())] = &snarl;
    }
    return index;
}
    
unordered_map<pair<int64_t, bool>, const Snarl*> SnarlManager::snarl_end_index() const {
    unordered_map<pair<int64_t, bool>, const Snarl*> index;
    for (const SnarlRecord& snarl_record : snarls) {
        const Snarl& snarl = *unrecord(&snarl_record);
        index[make_pair(snarl.end().node_id(), !snarl.end().backward())] = &snarl;
    }
    return index;
}
    
unordered_map<pair<int64_t, bool>, const Snarl*> SnarlManager::snarl_start_index() const {
    unordered_map<pair<int64_t, bool>, const Snarl*> index;
    for (const SnarlRecord& snarl_record : snarls) {
        const Snarl& snarl = *unrecord(&snarl_record);
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
            rec.parent = unrecord(parent);
        }
        else {
            // record top level status
#ifdef debug
            cerr << "\tSnarl is top-level" << endl;
#endif
            roots.push_back(&snarl);
            
            rec.parent = nullptr;
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
    
pair<unordered_set<id_t>, unordered_set<edge_t> > SnarlManager::shallow_contents(const Snarl* snarl, const HandleGraph& graph,
                                                                                 bool include_boundary_nodes) const {
    
    pair<unordered_set<id_t>, unordered_set<edge_t> > to_return;
        
    unordered_set<id_t> already_stacked;
        
    // initialize stack for DFS traversal of site
    vector<handle_t> stack;
        
    handle_t start_node = graph.get_handle(snarl->start().node_id());
    handle_t end_node = graph.get_handle(snarl->end().node_id());
        
    // mark the boundary nodes as already stacked so that paths will terminate on them
    already_stacked.insert(graph.get_id(start_node));
    already_stacked.insert(graph.get_id(end_node));
        
    // add boundary nodes as directed
    if (include_boundary_nodes) {
        to_return.first.insert(graph.get_id(start_node));
        to_return.first.insert(graph.get_id(end_node));
    }

    // stack up the nodes one edge inside the snarl from the start
    graph.follow_edges(start_node, snarl->start().backward(), [&](const handle_t& node) {            

            if (!already_stacked.count(graph.get_id(node))) {
                stack.push_back(node);
                already_stacked.insert(graph.get_id(node));
            }
            if (snarl->start().backward()) {
                to_return.second.insert(graph.edge_handle(node, start_node));
            } else {
                to_return.second.insert(graph.edge_handle(start_node, node));
            }
        });
      
    // stack up the nodes one edge inside the snarl from the end
    graph.follow_edges(end_node, !snarl->end().backward(), [&](const handle_t& node) {
            
            if (!already_stacked.count(graph.get_id(node))) {
                stack.push_back(node);
                already_stacked.insert(graph.get_id(node));
            }
            if (snarl->end().backward()) {
                to_return.second.insert(graph.edge_handle(end_node, node));
            } else {
                to_return.second.insert(graph.edge_handle(node, end_node));
            }
        });
        
    // traverse the snarl with DFS, skipping over any child snarls
    // do not pay attention to valid walks since we also want to discover any tips
    while (stack.size()) {
            
        // pop the top node off the stack
        handle_t node = stack.back();
        stack.pop_back();
            
        // record that this node is in the snarl
        to_return.first.insert(graph.get_id(node));
            
        const Snarl* forward_snarl = into_which_snarl(graph.get_id(node), false);
        const Snarl* backward_snarl = into_which_snarl(graph.get_id(node), true);
        if (forward_snarl) {
            // this node points into a snarl
                
            // What's on the other side of the snarl?
            id_t other_id = forward_snarl->start().node_id() == graph.get_id(node) ? forward_snarl->end().node_id() : forward_snarl->start().node_id();
                
            // stack up the node on the opposite side of the snarl
            // rather than traversing it
            handle_t opposite_node = graph.get_handle(other_id);
            if (!already_stacked.count(other_id)) {
                stack.push_back(opposite_node);
                already_stacked.insert(other_id);
            }
        }
            
        if (backward_snarl) {
            // the reverse of this node points into a snarl
                
            // What's on the other side of the snarl?
            id_t other_id = backward_snarl->end().node_id() == graph.get_id(node) ? backward_snarl->start().node_id(): backward_snarl->end().node_id();
                
            // stack up the node on the opposite side of the snarl
            // rather than traversing it
            handle_t opposite_node = graph.get_handle(other_id);
            if (!already_stacked.count(other_id)) {
                stack.push_back(opposite_node);
                already_stacked.insert(other_id);
            }
        }
            
        graph.follow_edges(node, false, [&](const handle_t& next_node) {
                edge_t edge = graph.edge_handle(node, next_node);
                // does this edge point forward or backward?
                if ((graph.get_is_reverse(node) && !backward_snarl) ||
                    (!graph.get_is_reverse(node) && !forward_snarl)) {

                        to_return.second.insert(edge);
                        
                        if (!already_stacked.count(graph.get_id(next_node))) {
                            
                            stack.push_back(next_node);
                            already_stacked.insert(graph.get_id(next_node));
                        }
                }
            });
        
        graph.follow_edges(node, true, [&](const handle_t& prev_node) {
                edge_t edge = graph.edge_handle(prev_node, node);
                // does this edge point forward or backward?
                if ((graph.get_is_reverse(node) && !forward_snarl) ||
                    (!graph.get_is_reverse(node) && !backward_snarl)) {

                    to_return.second.insert(edge);
                        
                    if (!already_stacked.count(graph.get_id(prev_node))) {
                            
                        stack.push_back(prev_node);
                        already_stacked.insert(graph.get_id(prev_node));
                    }
                }
            });
    }
        
    return to_return;
}
    
pair<unordered_set<id_t>, unordered_set<edge_t> > SnarlManager::deep_contents(const Snarl* snarl, const HandleGraph& graph,
                                                                              bool include_boundary_nodes) const {
        
    pair<unordered_set<id_t>, unordered_set<edge_t> > to_return;
        
    unordered_set<id_t> already_stacked;
        
    // initialize stack for DFS traversal of site
    vector<handle_t> stack;

    handle_t start_node = graph.get_handle(snarl->start().node_id());
    handle_t end_node = graph.get_handle(snarl->end().node_id());
        
    // mark the boundary nodes as already stacked so that paths will terminate on them
    already_stacked.insert(graph.get_id(start_node));
    already_stacked.insert(graph.get_id(end_node));
        
    // add boundary nodes as directed
    if (include_boundary_nodes) {
        to_return.first.insert(graph.get_id(start_node));
        to_return.first.insert(graph.get_id(end_node));
    }

    // stack up the nodes one edge inside the snarl from the start
    graph.follow_edges(start_node, snarl->start().backward(), [&](const handle_t& node) {            

            if (!already_stacked.count(graph.get_id(node))) {
                stack.push_back(node);
                already_stacked.insert(graph.get_id(node));
            }
            if (snarl->start().backward()) {
                to_return.second.insert(graph.edge_handle(node, start_node));
            } else {
                to_return.second.insert(graph.edge_handle(start_node, node));
            }
        });
      
    // stack up the nodes one edge inside the snarl from the end
    graph.follow_edges(end_node, !snarl->end().backward(), [&](const handle_t& node) {
            
            if (!already_stacked.count(graph.get_id(node))) {
                stack.push_back(node);
                already_stacked.insert(graph.get_id(node));
            }
            if (snarl->end().backward()) {
                to_return.second.insert(graph.edge_handle(end_node, node));
            } else {
                to_return.second.insert(graph.edge_handle(node, end_node));
            }
        });
        
    // traverse the snarl with DFS, skipping over any child snarls
    // do not pay attention to valid walks since we also want to discover any tips
    while (stack.size()) {
            
        // pop the top node off the stack
        handle_t node = stack.back();
        stack.pop_back();
            
        // record that this node is in the snarl
        to_return.first.insert(graph.get_id(node));

        graph.follow_edges(node, false, [&] (const handle_t& next_node) {
            edge_t edge = graph.edge_handle(node, next_node);
            to_return.second.insert(edge);
            if (!already_stacked.count(graph.get_id(next_node))) {
                stack.push_back(next_node);
                already_stacked.insert(graph.get_id(next_node));
            }
            });
        
        graph.follow_edges(node, true, [&] (const handle_t& prev_node) {
            edge_t edge = graph.edge_handle(prev_node, node);
            to_return.second.insert(edge);
            if (!already_stacked.count(graph.get_id(prev_node))) {
                stack.push_back(prev_node);
                already_stacked.insert(graph.get_id(prev_node));
            }
            });
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
    
vector<Visit> SnarlManager::visits_right(const Visit& visit, const HandleGraph& graph, const Snarl* in_snarl) const {
        
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

    graph.follow_edges(graph.get_handle(right_side.node), !right_side.is_end, [&](const handle_t& next_handle) {
        // For every NodeSide attached to the right side of this visit
        NodeSide attached(graph.get_id(next_handle), right_side.is_end ? graph.get_is_reverse(next_handle) : !graph.get_is_reverse(next_handle));
            
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
    });
        
    return to_return;
        
}
    
vector<Visit> SnarlManager::visits_left(const Visit& visit, const HandleGraph& graph, const Snarl* in_snarl) const {
        
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
    
#ifdef debug
    cerr << "Creating net graph of " << graph->get_id(this->start) << (graph->get_is_reverse(this->start) ? "-" : "+")
        << "->" << graph->get_id(this->end) << (graph->get_is_reverse(this->end) ? "-" : "+") << endl;
#endif
    
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
    
#ifdef debug
    cerr << "\tAdd unary child snarl on " << graph->get_id(snarl_bound) << (graph->get_is_reverse(snarl_bound) ? "-" : "+") << endl;
#endif
        
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
    auto start_visit = get_start_of(chain);
    handle_t chain_start_handle = graph->get_handle(start_visit.node_id(), start_visit.backward());
    auto end_visit = get_end_of(chain);
    handle_t chain_end_handle = graph->get_handle(end_visit.node_id(), end_visit.backward());
        
    // Save the links that let us cross the chain.
    chain_ends_by_start[chain_start_handle] = chain_end_handle;
    chain_end_rewrites[graph->flip(chain_end_handle)] = graph->flip(chain_start_handle);
    
#ifdef debug
    cerr << "\tAdd child chain " << graph->get_id(chain_start_handle) << (graph->get_is_reverse(chain_start_handle) ? "-" : "+")
        << " -> " << graph->get_id(chain_end_handle) << (graph->get_is_reverse(chain_end_handle) ? "-" : "+") << endl;
#endif
        
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

bool NetGraph::has_node(id_t node_id) const {
    return graph->has_node(node_id);
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
    
bool NetGraph::follow_edges_impl(const handle_t& handle, bool go_left, const function<bool(const handle_t&)>& iteratee) const {
    // Now we do the real work.
        
#ifdef debug
    cerr << "Look for edges in net graph of " << graph->get_id(start) << (graph->get_is_reverse(start) ? "-" : "+")
        << "->" << graph->get_id(end) << (graph->get_is_reverse(end) ? "-" : "+") << " on " << graph->get_id(handle) << (graph->get_is_reverse(handle) ? "-" : "+")
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
            
#ifdef debug
            cerr << "\tRead into chain end; warp to start" << endl;
#endif
            
            // Warp to the start.
            real_handle = chain_end_rewrites.at(other);
        } else if (chain_end_rewrites.count(graph->flip(other))) {
            // We're backing into the end of a chain.
            
#ifdef debug
            cerr << "\tBack into chain end; warp to start" << endl;
#endif
            
            // Warp to the start.
            real_handle = graph->flip(chain_end_rewrites.at(graph->flip(other)));
        }
            
#ifdef debug
        cerr << "\tFound edge " << (go_left ? "from " : "to ") << graph->get_id(other) << (graph->get_is_reverse(other) ? "-" : "+") << endl;
#endif
            
        if (!seen.count(real_handle)) {
            seen.insert(real_handle);
#ifdef debug
            cerr << "\t\tReport as " << graph->get_id(real_handle) << (graph->get_is_reverse(real_handle) ? "-" : "+") << endl;
#endif
                
            return iteratee(real_handle);
        } else {
#ifdef debug
            cerr << "\t\tEdge has been seen" << endl;
#endif
            return true;
        }
    };
        
    // This does the same as handle_edge, but flips the real handle
    auto flip_and_handle_edge = [&](const handle_t& other) -> bool {
            
        handle_t real_handle = other;
        if (chain_end_rewrites.count(other)) {
            // We're reading into the end of a chain.
#ifdef debug
        cerr << "\tRead into chain end; warp to start" << endl;
#endif
            // Warp to the start.
            real_handle = chain_end_rewrites.at(other);
        } else if (chain_end_rewrites.count(graph->flip(other))) {
            // We're backing into the end of a chain.
            
#ifdef debug
            cerr << "\tBack into chain end; warp to start" << endl;
#endif
            
            // Warp to the start.
            real_handle = graph->flip(chain_end_rewrites.at(graph->flip(other)));
        }
            
        real_handle = graph->flip(real_handle);
            
#ifdef debug
        cerr << "\tFound edge " << (go_left ? "from " : "to ") << graph->get_id(other) << (graph->get_is_reverse(other) ? "-" : "+") << endl;
#endif
            
        if (!seen.count(real_handle)) {
            seen.insert(real_handle);
#ifdef debug
            cerr << "\t\tReport as " << graph->get_id(real_handle) << (graph->get_is_reverse(real_handle) ? "-" : "+") << endl;
#endif
                
            return iteratee(real_handle);
        } else {
#ifdef debug
            cerr << "\t\tEdge has been seen" << endl;
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
        cerr << "\tWe are at the bound of the graph so don't say anything" << endl;
#endif
        return true;
    }
        
    if (chain_ends_by_start.count(handle) || chain_ends_by_start.count(graph->flip(handle))) {
        // If we have an associated chain end for this start, we have to use chain connectivity to decide what to do.
            
#ifdef debug
        cerr << "\tWe are a chain node" << endl;
#endif
            
        bool connected_start_start;
        bool connected_end_end;
        bool connected_start_end;
        tie(connected_start_start, connected_end_end, connected_start_end) = connectivity.at(graph->get_id(handle));
            
#ifdef debug
        cerr << "\t\tConnectivity: " << connected_start_start << " " << connected_end_end << " " << connected_start_end << endl;
#endif
            
        if (chain_ends_by_start.count(handle)) {
            // We visit the chain in its forward orientation
                
#ifdef debug
            cerr << "\t\tWe are visiting the chain forward" << endl;
#endif
                
            if (go_left) {
                // We want predecessors.
                // So we care about end-end connectivity (how could we have left our end?)
                    
#ifdef debug
                cerr << "\t\t\tWe are going left from a forward chain" << endl;
#endif
                    
                if (connected_end_end) {
                    
#ifdef debug
                    cerr << "\t\t\t\tWe can reverse and go back out the end" << endl;
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
                    cerr << "\t\t\t\tWe can continue through and go out the start" << endl;
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
                cerr << "\t\t\tWe are going right from a forward chain" << endl;
#endif
                    
                if (connected_start_start) {
                    
#ifdef debug
                    cerr << "\t\t\t\tWe can reverse and go back out the start" << endl;
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
                    cerr << "\t\t\t\tWe can continue through and go out the end" << endl;
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
            cerr << "\t\tWe are visiting the chain in reverse" << endl;
#endif

            if (go_left) {
                // We want predecessors of the reverse version (successors, but flipped)
                    
#ifdef debug
                cerr << "\t\t\tWe are going left from a reverse chain" << endl;
#endif
                    
                if (connected_start_start) {
                    
#ifdef debug
                    cerr << "\t\t\t\tWe can reverse and go back out the start" << endl;
#endif
                    
                    if (!graph->follow_edges(handle, false, flip_and_handle_edge)) {
                        // Iteratee is done
                        return false;
                    }
                }
                    
                if (connected_start_end) {
                    
#ifdef debug
                    cerr << "\t\t\t\tWe can continue through and go out the end" << endl;
#endif
                    
                    if (!graph->follow_edges(chain_ends_by_start.at(graph->flip(handle)), false, flip_and_handle_edge)) {
                        // Iteratee is done
                        return false;
                    }
                }
                    
            } else {
                // We want successors of the reverse version (predecessors, but flipped)
                    
#ifdef debug
                cerr << "\t\t\tWe are going right from a reverse chain" << endl;
#endif
                    
                if (connected_end_end) {
                    
#ifdef debug
                    cerr << "\t\t\t\tWe can reverse and go back out the end" << endl;
#endif
                    
                    if (!graph->follow_edges(chain_ends_by_start.at(graph->flip(handle)), false, handle_edge)) {
                        // Iteratee is done
                        return false;
                    }
                }
                    
                if (connected_start_end) {
                    
#ifdef debug
                    cerr << "\t\t\t\tWe can continue through and go out the start" << endl;
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
        cerr << "\tWe are looking at a unary snarl" << endl;
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
    cerr << "\tWe are an ordinary node" << endl;
#endif
        
    // Otherwise, this is an ordinary snarl content node
    return graph->follow_edges(handle, go_left, handle_edge);
}
    
bool NetGraph::for_each_handle_impl(const function<bool(const handle_t&)>& iteratee, bool parallel) const {
    // Find all the handles by a traversal.
        
#ifdef debug
    cerr << "Look for contents of net graph of " << graph->get_id(start) << (graph->get_is_reverse(start) ? "-" : "+")
        << "->" << graph->get_id(end) << (graph->get_is_reverse(end) ? "-" : "+") << endl;
#endif
        
    // We have to do the traversal on the underlying backing graph, because
    // the traversal functions we implemented on the graph we present will
    // maybe use internal child snarl connectivity, which can mean parts of
    // the graph are in this snarl but not actually reachable.
        
    // We let both the starts and ends of child snarls into the queue.
    // But we'll only reveal the starts to our iteratee.
    list<handle_t> queue;
    unordered_set<id_t> queued;
        
    // We define a function to queue up nodes we could visit next
    auto see_node = [&](const handle_t& other) {
        // Whenever we see a new node, add it to the queue
        auto found = queued.find(graph->get_id(other));
        if (found == queued.end()) {

#ifdef debug
            cerr << "\t\t\tFound new contained node " << graph->get_id(other) << (graph->get_is_reverse(other) ? "-" : "+") << endl;
#endif
        
            queue.push_back(other);
            queued.emplace_hint(found, graph->get_id(other));
        }
    };
        
    // Start at both the start and the end of the snarl.
    see_node(start);
    see_node(end);
        
    while (!queue.empty()) {
        handle_t here = queue.front();
        queue.pop_front();
            
#ifdef debug
        cerr << "\tVisit node " << graph->get_id(here) << (graph->get_is_reverse(here) ? "-" : "+") << endl;
#endif
            
        if (unary_boundaries.count(graph->flip(here))) {
            // This is a backward unary child snarl head, so we need to look at it the other way around.
            here = graph->flip(here);
            
#ifdef debug
            cerr << "\t\tReverse to match unary child boundary" << endl;
#endif
            
        } else if (chain_ends_by_start.count(graph->flip(here))) {
            // This is a backward child chain head, so we need to look at it the other way around.
            here = graph->flip(here);
            
#ifdef debug
            cerr << "\t\tReverse to match child chain head" << endl;
#endif
            
        } else if (chain_end_rewrites.count(graph->flip(here))) {
            // This is a backward child chain tail, so we need to look at it the other way around.
            here = graph->flip(here);
            
#ifdef debug
            cerr << "\t\tReverse to match child chain tail" << endl;
#endif
            
        }
            
        if (!chain_end_rewrites.count(here)) {
            // This is not a chain end, so it's either a real contained node or a chain head.
            // We can emit it.
            
#ifdef debug
            cerr << "\t\tVisit forward version" << endl;
#endif
                
            if (graph->get_is_reverse(here)) {
                if (!iteratee(graph->flip(here))) {
                    // Run the iteratee on the forward version, and stop if it wants to stop
                    return false;
                }
            } else {
                if (!iteratee(here)) {
                    // Run the iteratee, and stop if it wants to stop.
                    return false;
                }
            }
                
        } else {
#ifdef debug
            cerr << "\t\tSkip chain end but see start at " << graph->get_id(chain_end_rewrites.at(here)) << (graph->get_is_reverse(chain_end_rewrites.at(here)) ? "-" : "+") << endl;
#endif
            // If we reach a chain end, make sure to eventually visit the chain start.
            // There might not be any other edges to it.
            see_node(chain_end_rewrites.at(here));
        }
            
        // We already have flipped any backward heads or tails frontward. So
        // we don't need to check if the backward version of us is in
        // anything.
            
        if ( ((start != end && here != end && here != graph->flip(start)) ||
               start == end)
            && !unary_boundaries.count(here) &&
            !chain_ends_by_start.count(here)  && !chain_end_rewrites.count(here)) {
            
#ifdef debug
            cerr << "\t\tRight side faces into net graph" << endl;
#endif
                
            // We have normal graph to our right and not the exterior of this snarl or the interior of a child.
            graph->follow_edges(here, false, see_node);
        }
            
        if ((start != end && here != start && here != graph->flip(end)) ||
             start == end) {
             
#ifdef debug
            cerr << "\t\tLeft side faces into net graph" << endl;
#endif
             
            // We have normal graph to our left.
            graph->follow_edges(here, true, see_node);
        }
            
        if (chain_end_rewrites.count(here)) {
        
#ifdef debug
            cerr << "\t\tWe are chain end; look right off reverse start at " << graph->get_id(chain_end_rewrites.at(here)) << (graph->get_is_reverse(chain_end_rewrites.at(here)) ? "-" : "+") << endl;
#endif
        
            // We need to look right off the reverse head of this child snarl.
            graph->follow_edges(chain_end_rewrites.at(here), false, see_node);
        }
            
        if (chain_ends_by_start.count(here)) {
        
#ifdef debug
            cerr << "\t\tWe are chain start; look right off end at " << graph->get_id(chain_ends_by_start.at(here)) << (graph->get_is_reverse(chain_ends_by_start.at(here)) ? "-" : "+") << endl;
#endif
        
            // We need to look right off the (reverse) tail of this child snarl.
            graph->follow_edges(chain_ends_by_start.at(here), false, see_node);
        }
    }
    
    return true;
}
    
size_t NetGraph::get_node_count() const {
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

edge_t to_edge(const HandleGraph& graph, const Visit& v1, const Visit& v2) {

    id_t prev_id;
    bool prev_back;
    if (v1.node_id() != 0) {
        prev_id = v1.node_id();
        prev_back = v1.backward();
    } else {
        const Snarl& prev_snarl = v1.snarl();
        if (v1.backward()) {
            prev_id = prev_snarl.start().node_id();
            prev_back = !prev_snarl.start().backward();
        } else {
            prev_id = prev_snarl.end().node_id();
            prev_back = prev_snarl.end().backward();
        }
    }
    
    id_t cur_id;
    bool cur_back;                
    if (v2.node_id() != 0) {
        cur_id = v2.node_id();
        cur_back = v2.backward();
    } else {
        const Snarl& cur_snarl = v2.snarl();
        if (v2.backward()) {
            cur_id = cur_snarl.end().node_id();
            cur_back = !cur_snarl.end().backward();
        } else {
            cur_id = cur_snarl.start().node_id();
            cur_back = cur_snarl.start().backward();
        }
    }

    return graph.edge_handle(graph.get_handle(prev_id, prev_back),
                             graph.get_handle(cur_id, cur_back));

}

bool operator==(const Visit& a, const Visit& b) {
    // IDs and orientations have to match, and nobody has a snarl or the
    // snarls match.
    return a.node_id() == b.node_id() &&
        a.name() == b.name() &&
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
        return make_tuple(a.node_id(), a.backward(), a.name()) < make_tuple(b.node_id(), b.backward(), b.name());
    } else {
        // Compare including the snarl
        return make_tuple(a.node_id(), a.snarl(), a.backward(), a.name()) < make_tuple(b.node_id(), b.snarl(), b.backward(), b.name());
    }        
}
    
ostream& operator<<(ostream& out, const Visit& visit) {
    if (!visit.has_snarl()) {
        if (visit.name().empty()) {
            // Use the node ID
            out << visit.node_id();
        } else {
            // Use the name
            out << visit.name();
        }
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






