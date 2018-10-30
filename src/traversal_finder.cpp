#include "traversal_finder.hpp"
#include "genotypekit.hpp"
#include "algorithms/topological_sort.hpp"
#include "algorithms/is_acyclic.hpp"
#include "cactus.hpp"

//#define debug

namespace vg {

using namespace std;

PathBasedTraversalFinder::PathBasedTraversalFinder(vg::VG& g, SnarlManager& sm) : graph(g), snarlmanager(sm){
}

vector<SnarlTraversal> PathBasedTraversalFinder::find_traversals(const Snarl& site){
    // Goal: enumerate traversals in the snarl supported by paths in the graph
    // that may not cover the ends of the snarl.
    // Label the Snarl's name as tthe hash of the variant and the SnarlTraversal's name
    // as the name of the alt_path (i.e. "_alt_[a-z0-9]*_[0-9]*")s
    vector<SnarlTraversal> ret;

    // If the snarl is not an ultrabubble, just return an empty set of traversals.
    if (site.type() != ULTRABUBBLE){
        return ret;
    }
   
    // Get the Snarl's nodes
    unordered_set<int64_t> snarl_node_ids;
    pair<unordered_set<Node*>, unordered_set<Edge*> > contents = snarlmanager.shallow_contents(&site, graph, true);


    // Get the variant paths at the snarl nodes.
    set<string> var_path_names;
    regex front ("(_alt_)(.*)");
    regex alt_str ("(_alt_)");
    regex back ("(_[0-9]*)");
    auto& gpaths = graph.paths._paths;
    set<string> gpath_names;
    for (auto x : gpaths){
        gpath_names.insert(x.first);
    }
    map<string, set<string> > basename_to_pathnames;
    map<string, bool> path_processed;


    // Collect our paths which cross our snarl's nodes.
    for (auto node : contents.first){
        //cerr << "Processing node " << id << endl;
        set<string> p_of_n = graph.paths.of_node(node->id());

        for (auto pn : p_of_n){
            if (!std::regex_match(pn, front)){
                // don't include reference paths
                continue;
            }
            string variant_hash = std::regex_replace(pn, alt_str, "");
            variant_hash = std::regex_replace(variant_hash, back, "");
            //cerr << variant_hash << endl;
            regex varbase(variant_hash);

            path_processed[pn] = false;
            basename_to_pathnames[variant_hash].insert(pn);
            var_path_names.insert(pn);
            for (auto g : gpath_names){
                if (std::regex_search(g, varbase)){
                    basename_to_pathnames[variant_hash].insert(g);
                    path_processed[g] = false;
                    var_path_names.insert(g);
                }
            }
        }
    }
    // for (auto p : var_path_names){
    //     cerr << p << endl;
    // }
    // exit(1);
    for (auto cpath : var_path_names){

        //cerr << "Working on path " << cpath << endl;
        if (!std::regex_match(cpath, front) || path_processed[cpath]){
            cerr << "Path already processed " << cpath << endl;
            continue;
        }

        if (std::regex_match(cpath, front)){
            // We found an alt path at this location
            // We need to check if it has any paths in the graph with no paths.
            string variant_hash = std::regex_replace(cpath, alt_str, "");
            variant_hash = std::regex_replace(variant_hash, back, "");
            set<string> allele_path_names = basename_to_pathnames[variant_hash];
            for (auto a : allele_path_names){
                //cerr << "Processing path " << a << endl;
                // for each allele, generate a traversal
                SnarlTraversal fresh_trav;
                fresh_trav.set_name(a);
                
                // Add the start node to the traversal
                *fresh_trav.add_visit() = site.start();
                // Fill in our traversal
                auto& ms = gpaths[a];
                for (auto m : ms){
                    int64_t n_id = m.node_id();
                    bool backward = m.is_reverse();
                    Visit* v = fresh_trav.add_visit();
                    v->set_node_id(n_id);
                    v->set_backward(backward);
                }
                // Add the end node to the traversal
                *fresh_trav.add_visit() = site.end();
                ret.push_back(fresh_trav);
                //cerr << "Finished path: " << a << endl;
                path_processed[a] = true;
            }

        }
    }
        

    // Check to make sure we got all our paths
    for (auto p : path_processed){
        if (!path_processed[p.first]){
            cerr << "VARIANT PATH MISSED: " << p.first << endl;
            exit(1617);
        }
    }

    return ret;
}
   
ExhaustiveTraversalFinder::ExhaustiveTraversalFinder(VG& graph, SnarlManager& snarl_manager,
                                                     bool include_reversing_traversals) :
    graph(graph), snarl_manager(snarl_manager),
    include_reversing_traversals(include_reversing_traversals) {
    // nothing more to do
}
    
ExhaustiveTraversalFinder::~ExhaustiveTraversalFinder() {
    // no heap objects
}

void ExhaustiveTraversalFinder::stack_up_valid_walks(NodeTraversal walk_head, vector<NodeTraversal>& stack) {
    
    id_t head_id = walk_head.node->id();

    // get all edges involving this node so we can filter them down to valid walks
    for (Edge* edge : graph.edges_of(walk_head.node)) {
        Node* next_node = nullptr;
        bool next_backward;
        bool from_start;
        // determine id and orientation of our nodes given that they can
        // be either from or to in the edge
        if (head_id == edge->from()) {
            next_node = graph.get_node(edge->to());
            next_backward = edge->to_end();
            from_start = edge->from_start();
        } else {
            next_node = graph.get_node(edge->from());
            next_backward = !edge->from_start();
            from_start = !edge->to_end();
        }
        // are we walking the same direction relative to head_id?
        if (walk_head.backward == from_start &&
            // derived classes can use this to filter search
            visit_next_node(next_node, edge)) {
            // add the next traversal in the walk to the stack            
            stack.push_back(NodeTraversal(next_node, next_backward));
        }
    }
}

void ExhaustiveTraversalFinder::add_traversals(vector<SnarlTraversal>& traversals,
                                               NodeTraversal traversal_start,
                                               set<NodeTraversal>& stop_at,
                                               set<NodeTraversal>& yield_at) {
    // keeps track of the walk of the DFS traversal
    list<Visit> path;
    
    // these mark the start of the edges out of the node that is on the head of the path
    // they can be used to see how many nodes we need to peel off the path when we're
    // backtracking
    NodeTraversal stack_sentinel(nullptr);
    
    // initialize stack for DFS traversal of site
    vector<NodeTraversal> stack{traversal_start};
    
    while (stack.size()) {
        
        NodeTraversal node_traversal = stack.back();
        stack.pop_back();
        
        // we have traversed all of edges out of the head of the path, so we can pop it off
        if (node_traversal == stack_sentinel) {
            path.pop_back();
            continue;
        }
        
        // have we finished a traversal through the site?
        if (stop_at.count(node_traversal)) {
            if (yield_at.count(node_traversal)) {
                // yield path as a snarl traversal
                traversals.emplace_back();
                
                // record the traversal in the return value
                for (auto iter = path.begin(); iter != path.end(); iter++) {
                    *traversals.back().add_visit() = *iter;
                }
                // add the final visit
                *traversals.back().add_visit() = to_visit(node_traversal);
            }
            
            // don't proceed to add more onto the DFS stack
            continue;
        }
        
        // mark the beginning of this node's edges forward in the stack
        stack.push_back(stack_sentinel);
        
        // make a visit through the node traversal and add it to the path
        path.emplace_back();
        path.back().set_node_id(node_traversal.node->id());
        path.back().set_backward(node_traversal.backward);
        
        // does this traversal point into a child snarl?
        const Snarl* into_snarl = snarl_manager.into_which_snarl(node_traversal.node->id(),
                                                                 node_traversal.backward);
                                                                 
#ifdef debug
        cerr << "Traversal " << node_traversal.node->id() << " " << node_traversal.backward << " enters";
        if (into_snarl != nullptr) {
            cerr << " " << pb2json(*into_snarl) << endl;
        } else {
            cerr << " NULL" << endl; 
        }
#endif
                                                                 
        if (into_snarl && !(node_traversal == traversal_start)) {
            // add a visit for the child snarl
            path.emplace_back();
            *path.back().mutable_snarl()->mutable_start() = into_snarl->start();
            *path.back().mutable_snarl()->mutable_end() = into_snarl->end();
            
            // mark the beginning of this child snarls edges forward in the stack
            stack.push_back(stack_sentinel);
            
            // which side of the snarl does the traversal point into?
            if (into_snarl->start().node_id() == node_traversal.node->id()
                && into_snarl->start().backward() == node_traversal.backward) {
                // Into the start
#ifdef debug
                cerr << "Entered child through its start" << endl;
#endif
                if (into_snarl->start_end_reachable()) {
                    // skip to the other side and proceed in the orientation that the end node takes.
                    stack.push_back(to_node_traversal(into_snarl->end(), graph));
#ifdef debug
                    cerr << "Stack up " << stack.back().node->id() << " " << stack.back().backward << endl;
#endif
                }
                
                // if the same side is also reachable, add it to the stack too
                if (into_snarl->start_self_reachable()) {
                    // Make sure to flip it around so we come out of the snarl instead of going in again,
                    stack.push_back(to_rev_node_traversal(into_snarl->start(), graph).reverse());
#ifdef debug
                    cerr << "Stack up " << stack.back().node->id() << " " << stack.back().backward << endl;
#endif
                }
                
            }
            else {
                // Into the end
#ifdef debug
                cerr << "Entered child through its end" << endl;
#endif
                if (into_snarl->start_end_reachable()) {
                    // skip to the other side and proceed in the orientation
                    // *opposite* what the start node takes (i.e. out of the
                    // snarl)
                    stack.push_back(to_node_traversal(into_snarl->start(), graph).reverse());
#ifdef debug
                    cerr << "Stack up " << stack.back().node->id() << " " << stack.back().backward << endl;
#endif
                }
                
                // if the same side is also reachable, add it to the stack too
                if (into_snarl->end_self_reachable()) {
                    stack.push_back(to_rev_node_traversal(into_snarl->end(), graph));
#ifdef debug
                    cerr << "Stack up " << stack.back().node->id() << " " << stack.back().backward << endl;
#endif
                }
            }
        }
        else {
            // add all of the node traversals we can reach through valid walks to stack
            stack_up_valid_walks(node_traversal, stack);
        }
    }
}
    
vector<SnarlTraversal> ExhaustiveTraversalFinder::find_traversals(const Snarl& site) {

    vector<SnarlTraversal> to_return;
    
    NodeTraversal site_end = to_node_traversal(site.end(), graph);
    NodeTraversal site_start = to_node_traversal(site.start(), graph);
    NodeTraversal site_rev_start = NodeTraversal(site_start.node, !site_start.backward);
    
    // stop searching when the traversal is leaving the site
    set<NodeTraversal> stop_at;
    stop_at.insert(site_end);
    stop_at.insert(site_rev_start);
    
    // choose which side(s) can be the end of the traversal
    set<NodeTraversal> yield_at;
    yield_at.insert(site_end);
    if (include_reversing_traversals) {
        yield_at.insert(site_rev_start);
    }
    
    // search forward from the start and add any traversals that leave the indicated boundaries
    add_traversals(to_return, site_start, stop_at, yield_at);

    if (site.end_self_reachable() && include_reversing_traversals) {
        // if the end is reachable from itself, also look for traversals that both enter and
        // leave through the end
        yield_at.erase(site_rev_start);
        add_traversals(to_return, NodeTraversal(site_end.node, !site_end.backward),
                       stop_at, yield_at);
    }
    
    return to_return;
}

SupportRestrictedTraversalFinder::SupportRestrictedTraversalFinder(AugmentedGraph& augmented_graph,
                                                                   SnarlManager& snarl_manager,
                                                                   int min_node_support,
                                                                   int min_edge_support,
                                                                   bool include_reversing_traversals) :
    ExhaustiveTraversalFinder(augmented_graph.graph,
                              snarl_manager,
                              include_reversing_traversals),
    aug(augmented_graph),
    min_node_support(min_node_support),
    min_edge_support(min_edge_support) {
}

SupportRestrictedTraversalFinder::~SupportRestrictedTraversalFinder() {}

bool SupportRestrictedTraversalFinder::visit_next_node(const Node* node, const Edge* edge) {
    return aug.get_alignments(node->id()).size() >= min_node_support &&
        aug.get_alignments(NodeSide::pair_from_edge(*edge)).size() >= min_edge_support;
}


PathRestrictedTraversalFinder::PathRestrictedTraversalFinder(VG& graph,
                                                             SnarlManager& snarl_manager,
                                                             map<string, const Alignment*>& reads_by_name,
                                                             int min_recurrence, int max_path_search_steps) :
    graph(graph),
    snarl_manager(snarl_manager),
    reads_by_name(reads_by_name),
    min_recurrence(min_recurrence),
    max_path_search_steps(max_path_search_steps) {
    // nothing else to do
}

PathRestrictedTraversalFinder::~PathRestrictedTraversalFinder() {
    // no heap variables
}

static bool mapping_enters_side(const Mapping& mapping, const handle_t& side, const HandleGraph* graph) {

#ifdef debug
#pragma omp critical (cerr)
     cerr << "Does mapping " << pb2json(mapping) << " enter " << graph->get_id(side) << " " << graph->get_is_reverse(side) << endl;
#endif
    
    bool enters = mapping.position().node_id() == graph->get_id(side) &&
        mapping.position().offset() == 0;
       
#ifdef debug
#pragma omp critical (cerr)
    cerr << enters << endl;
#endif
    
    return enters;
}

static bool mapping_exits_side(const Mapping& mapping, const handle_t& side, const HandleGraph* graph) {
    
#ifdef debug
#pragma omp critical (cerr)
    cerr << "Does mapping " << pb2json(mapping) << " exit " << graph->get_id(side) << " " << graph->get_is_reverse(side) << endl;
#endif
    
    bool exits = mapping.position().node_id() == graph->get_id(side) &&
        mapping.position().offset() + mapping_from_length(mapping) == graph->get_length(side) &&
        mapping.edit(mapping.edit_size() - 1).from_length() ==
        mapping.edit(mapping.edit_size() - 1).to_length();
        
#ifdef debug
#pragma omp critical (cerr)
    cerr << exits << endl;
#endif
    
    return exits;
}

// replaces get_paths_through_site from genotyper
vector<SnarlTraversal> PathRestrictedTraversalFinder::find_traversals(const Snarl& site) {
    // We're going to emit traversals supported by any paths in the graph.
    
    // Put all our subpaths in here to deduplicate them by sequence they spell
    // out. And to count occurrences. Note that the occurrence count will be
    // boosted to min_recurrence if a non-read path in the graph supports a
    // certain traversal string, so we don't end up dropping unsupported ref
    // alleles.
    map<string, pair<SnarlTraversal, int>> results;

#ifdef debug
#pragma omp critical (cerr)
    cerr << "Looking for paths between " << site.start() << " and " << site.end() << endl;
#endif

    if(graph.paths.has_node_mapping(graph.get_node(site.start().node_id())) &&
        graph.paths.has_node_mapping(graph.get_node(site.end().node_id()))) {
        // If we have some paths that visit both ends (in some orientation)

        // Get all the mappings to the end node, by path name
        auto endmappings_by_name = graph.paths.get_node_mapping_by_path_name(graph.get_node(site.end().node_id()));

        for(auto name_and_mappings : graph.paths.get_node_mapping_by_path_name(graph.get_node(site.start().node_id()))) {
            // Go through the paths that visit the start node

            // Grab their names
            auto& name = name_and_mappings.first;

            if(!endmappings_by_name.count(name_and_mappings.first)) {
                //cerr << "no endmappings match" << endl;
                // No path by this name has any mappings to the end node. Skip
                // it early.
                continue;
            }

            for(auto* mapping : name_and_mappings.second) {
                // Start at each mapping in the appropriate orientation

#ifdef debug
#pragma omp critical (cerr)
                cerr << "Trying mapping of read/path " << name_and_mappings.first << " to " << mapping->node_id() << (mapping->is_reverse() ? "-" : "+") << endl;
#endif

                // How many times have we gone to the next mapping looking for a
                // mapping to the end node in the right orientation?
                size_t traversal_count = 0;

                // Do we want to go left (true) or right (false) from this
                // mapping? If start is a forward traversal and we found a
                // forward mapping, we go right. If either is backward we go
                // left, and if both are backward we go right again.
                bool traversal_direction = mapping->is_reverse() != site.start().backward();

#ifdef debug
#pragma omp critical (cerr)
                cerr << "Traversal direction: " << traversal_direction << endl;
#endif

                // Now work out if we are entering the snarl or not
                if (traversal_direction) {
                
                    // We are going left in the read but right in the snarl, so
                    // we want to enter the snarl's start node
                    bool enter_start = mapping_enters_side(mapping->to_mapping(), graph.get_handle(site.start()), &graph);

#ifdef debug
#pragma omp critical (cerr)
                    cerr << "Enter start: " << enter_start << endl;
#endif
                    
                    if (!enter_start) {
                        // We only want reads that enter the snarl
                        continue;
                    }
                    
                } else {
                    // We are going right, so we want to exit the snarl's start
                    // node
                    bool exit_start = mapping_exits_side(mapping->to_mapping(), graph.get_handle(site.start()), &graph);
                    
#ifdef debug
#pragma omp critical (cerr)
                    cerr << "Exit start: " << exit_start << endl;
#endif
                    
                    if (!exit_start) {
                        // We are only interested in reads that exit the snarl
                        continue;
                    }
                }

                // What orientation would we want to find the end node in? If
                // we're traveling backward, we expect to find it in the
                // opposnarl direction to the one we were given.
                bool expected_end_orientation = site.end().backward() != traversal_direction;

                // We're going to fill in this SnarlTraverasl with visits.
                SnarlTraversal path_traversed;

                while(mapping != nullptr && traversal_count < max_path_search_steps) {
                    // Traverse along until we hit the end traversal or take too
                    // many steps

#ifdef debug
#pragma omp critical (cerr)
                    cerr << "\tTraversing " << mapping->node_id() << (mapping->is_reverse() ? "-" : "+") << endl;
#endif

                    // Say we visit this node along the path, in this orientation
                    *path_traversed.add_visit() = to_visit(!traversal_direction ? mapping->to_mapping() :
                                                           reverse_complement_mapping(mapping->to_mapping(),[this](id_t node_id) {
                                                                   return this->graph.get_node(node_id)->sequence().length();
                                                               }), true);

                    if(mapping->node_id() == site.end().node_id() && mapping->is_reverse() == expected_end_orientation) {
                        // Does our mapping actually cross through the ending side?
                        // It has to either enter the end node, or exit the end
                        // node, depending on which way in the read we read. And
                        // if it doesn't we try again.
                        if (!traversal_direction &&
                            !mapping_enters_side(mapping->to_mapping(), graph.get_handle(site.end()), &graph) ||
                            traversal_direction && 
                            !mapping_exits_side(mapping->to_mapping(), graph.get_handle(site.end()), &graph)) {
                            break;
                        }

                        // get the string for the sequence
                        string allele_seq;
                        for (size_t i = 0; i < path_traversed.visit_size(); i++) {
                            auto& path_visit = path_traversed.visit(i);
                            Node* map_node = graph.get_node(path_visit.node_id());
                            allele_seq += path_visit.backward() ? reverse_complement(map_node->sequence()) : map_node->sequence();
                        }
                        
                        // We have stumbled upon the end node in the orientation we wanted it in.
                        if(results.count(allele_seq)) {
                            // It is already there! Increment the observation count.
#ifdef debug
#pragma omp critical (cerr)
                            cerr << "\tFinished; got known sequence " << allele_seq << endl;
#endif
                            
                            if(reads_by_name.count(name)) {
                                // We are a read. Just increment count
                                results[allele_seq].second++;
                            } else {
                                // We are a named path (like "ref")
                                if(results[allele_seq].second < min_recurrence) {
                                    // Ensure that this allele doesn't get
                                    // eliminated, since ref or some other named
                                    // path supports it.
                                    results[allele_seq].second = min_recurrence;
                                } else {
                                    results[allele_seq].second++;
                                }
                            }
                        } else {
                            // Add it in. Give it a count of 1 if we are a read,
                            // and a count of min_recurrence (so it doesn't get
                            // filtered later) if we are a named non-read path
                            // (like "ref").
                            results[allele_seq] = make_pair(path_traversed, reads_by_name.count(name) ? 1 : min_recurrence);
#ifdef debug
#pragma omp critical (cerr)
                            cerr << "\tFinished; got novel sequence " << allele_seq << endl;
#endif
                        }

                        // Then try the next embedded path
                        break;
                    }

                    // Otherwise just move to the right (or left)
                    if(traversal_direction) {
                        // We're going backwards
                        mapping = graph.paths.traverse_left(mapping);
                    } else {
                        // We're going forwards
                        mapping = graph.paths.traverse_right(mapping);
                    }
                    // Tick the counter so we don't go really far on long paths.
                    traversal_count++;

                }


            }
        }

    }

    // Now collect the unique results
    vector<SnarlTraversal> to_return;

    for(auto& result : results) {
        // Break out each result
        const string& seq = result.first;
        auto& traversals = result.second.first;
        auto& count = result.second.second;

        if(count < min_recurrence) {
            // We don't have enough initial hits for this sequence to justify
            // trying to re-align the rest of the reads. Skip it. Note that the
            // reference path (and other named paths) will stuff in at least
            // min_recurrence to make sure we don't throw out their alleles.
            continue;
        }

        // Send out each list of traversals
        to_return.emplace_back(std::move(traversals));
    }

    return to_return;
}

ReadRestrictedTraversalFinder::ReadRestrictedTraversalFinder(AugmentedGraph& augmented_graph,
                                                             SnarlManager& snarl_manager,
                                                             int min_recurrence, int max_path_search_steps) :
    aug(augmented_graph),
    snarl_manager(snarl_manager),
    min_recurrence(min_recurrence),
    max_path_search_steps(max_path_search_steps) {
    // nothing else to do
}

ReadRestrictedTraversalFinder::~ReadRestrictedTraversalFinder() {
    // no heap variables
}
    
// replaces get_paths_through_site from genotyper
vector<SnarlTraversal> ReadRestrictedTraversalFinder::find_traversals(const Snarl& site) {
    
    // We're going to emit traversals supported by reads in the AugmentedGraph.
    
    // Put all our subpaths in here to deduplicate them by sequence they spell
    // out. And to count occurrences. Note that the occurrence count will be
    // boosted to min_recurrence if a non-read path in the graph supports a
    // certain traversal string, so we don't end up dropping unsupported ref
    // alleles.
    map<string, pair<SnarlTraversal, int>> results;

#ifdef debug
#pragma omp critical (cerr)
    cerr << "Looking for reads between " << site.start() << " and " << site.end() << endl;
#endif

    // Get all the alignments that touch each end of the snarl.
    auto start_alignments = aug.get_alignments(site.start().node_id());
    auto end_alignments = aug.get_alignments(site.end().node_id());
    
#ifdef debug
#pragma omp critical (cerr)
    cerr << "Found " << start_alignments.size() << " reads on start and " << end_alignments.size() << " reads on end" << endl;
#endif
    
    // Turn one into a set
    unordered_set<const Alignment*> end_set{end_alignments.begin(), end_alignments.end()};

    for(const Alignment* alignment : start_alignments) {
        // Go through the alignments that visit the start node
        if(!end_set.count(alignment)) {
            // Skip alignments that don't also visit the end node
            continue;
        }
        
#ifdef debug
#pragma omp critical (cerr)
        cerr << "Alignment " << alignment->name() << " touches both start and end" << endl;
#endif

        // This alignment visits both sides of the snarl so it probably goes through it.

        // We will clip out the traversing part of the read, if any, into lists of visits.
        list<Visit> building;
        
        // We need to support multiple traversals of a snarl. So when we finish a traversal it goes in here.
        list<list<Visit>> completed;
        
        // We set this to true when we go inside the snarl, and false again when we leave.
        // Note that
        bool in_snarl = false;
        
        // We set this if we enter the snarl through the end and not through the start.
        bool in_end = false;

        for(const Mapping& mapping : alignment->path().mapping()) {
            // For each mapping in the read's path in order
            
            if (!in_snarl) {
                // If we aren't yet in the snarl, see if we can enter.
                if (mapping.position().node_id() == site.start().node_id() &&
                    mapping.position().is_reverse() == site.start().backward()) {
                    
                    // Enter through the start
                    in_snarl = true;
                    in_end = false;
                } else if (mapping.position().node_id() == site.end().node_id() &&
                    mapping.position().is_reverse() == !site.end().backward()) {
                    
                    // Enter through the end
                    in_snarl = true;
                    in_end = true;
                }
            }
            
            if (in_snarl) {
                // If we are now in the snarl, add this mapping to the list as a
                // visit. Make sure to convert mappings with offsets/edits to
                // visits, because we may not have cut the graph at alignment
                // ends.
                building.push_back(to_visit(mapping, true));
                
                // Now see if we should leave the snarl
                if ((mapping.position().node_id() == site.start().node_id() &&
                    mapping.position().is_reverse() == !site.start().backward()) ||
                    (mapping.position().node_id() == site.end().node_id() &&
                    mapping.position().is_reverse() == site.end().backward())) {
                    // We are leaving through the start or end
                    in_snarl = false;
                    
                    if (in_end) {
                        // Flip the traversal to the snarl's orientation
                        building.reverse();
                        for (auto& visit : building) {
                            visit = reverse(visit);
                        }
                    }
                    
                    // Put the traversal on the list of completed traversals
                    completed.emplace_back(std::move(building));
                    building.clear();
                }
            }
        }
        
        // Now turn all the traversals of the snarl that we found in this read into real SnarlTraversals
        for (auto& visits : completed) {
            // Convert each list of visits to a SnarlTraversal
            SnarlTraversal converted;
            for (auto& visit : visits) {
                // Transfer over all the visits
                *converted.add_visit() = visit;
            }
            
            // Work out what sequence we have found
            auto sequence = traversal_to_string(aug.graph, converted);
            
            if (results.count(sequence)) {
                // This is a known traversal so up the count
                results[sequence].second++;
            } else {
                // This is a new thing so add it with one occurrence
                results[sequence] = make_pair(converted, 1);
            }
        }        
    }

    // Now collect the unique results
    vector<SnarlTraversal> to_return;

    for(auto& result : results) {
        // Break out each result
        const string& seq = result.first;
        auto& traversals = result.second.first;
        auto& count = result.second.second;

        if(count < min_recurrence) {
            // We don't have enough initial hits for this sequence to justify
            // trying to re-align the rest of the reads. Skip it.
            continue;
        }

        // Send out each list of traversals
        to_return.emplace_back(std::move(traversals));
    }
    
#ifdef debug
#pragma omp critical (cerr)
    cerr << "Found " << to_return.size() << " traversals based on reads" << endl;
#endif

    return to_return;
}

TrivialTraversalFinder::TrivialTraversalFinder(VG& graph) : graph(graph) {
    // Nothing to do!
}

vector<SnarlTraversal> TrivialTraversalFinder::find_traversals(const Snarl& site) {
    assert(site.type() == ULTRABUBBLE);
    
    // We'll fill this in and send it back
    vector<SnarlTraversal> to_return;
    
    // We don't want to be duplicating partial paths, so we store for each
    // NodeTraversal we can reach the previous NodeTraversal we can reach it
    // from.
    map<NodeTraversal, NodeTraversal> previous;
    
    list<NodeTraversal> stack{to_node_traversal(site.start(), graph)};
    
    while (!stack.empty()) { 
        // While there's still stuff on the stack
        
        // Grab the first thing
        NodeTraversal here = stack.front();
        stack.pop_front();
        
        if (here.node->id() == site.end().node_id()) {
            // Trace back a path
            list<NodeTraversal> path;
            
            while (true) {
                // Until we get to the start of the site
                
                // Put this traversal on the front of the path
                path.push_front(here);
                
                if (here.node->id() == site.start().node_id()) {
                    // Stop when we've reached the start of the site
                    break;
                }
                
                // Trace back
                here = previous.at(here);
            }
            
            // Initialize a SnarlTraversal in the return value
            to_return.emplace_back();
            
            // Translate the path into the traversal
            for (NodeTraversal node_traversal : path) {
                *(to_return.back().add_visit()) = to_visit(node_traversal);
            }
            
            // Stop early after having found one path
            break;
        } else {
            // We haven't reached the end of the site
            
            for (NodeTraversal next : graph.nodes_next(here)) {
                // Look at all the places we can go from this node
                if (previous.count(next)) {
                    // We already know how to get there.
                    continue;
                }
                
                // Remember how we got there
                previous[next] = here;
                // Explore it, depth first
                stack.push_front(next);
            }
        }
    }
    
    // When we get here, either we found a path, or there isn't one.
    return to_return;
}


RepresentativeTraversalFinder::RepresentativeTraversalFinder(AugmentedGraph& augmented,
                                                             SnarlManager& snarl_manager, size_t max_depth, size_t max_width, size_t max_bubble_paths,
                                                             function<PathIndex*(const Snarl&)> get_index) : augmented(augmented), snarl_manager(snarl_manager),
                                                                                                             max_depth(max_depth), max_width(max_width), max_bubble_paths(max_bubble_paths), get_index(get_index) {
    
    // Nothing to do!

}

Path RepresentativeTraversalFinder::find_backbone(const Snarl& site) {
    
    // TODO: this cheats and uses certain things that happen to be true about
    // the TrivialTraversalFinder in order to work.

    // Find a traversal, ignoring the fact that child sites ought to own their
    // nodes.
    TrivialTraversalFinder finder(augmented.graph);
    auto traversals = finder.find_traversals(site);
    assert(!traversals.empty());
    auto& traversal = traversals.front();
    
    // Convert it into a path that includes the boundary nodes
    Path to_return;
    for (size_t i = 0; i < traversal.visit_size(); i++) {
        *to_return.add_mapping() = to_mapping(traversal.visit(i), augmented.graph);
    }
    
    return to_return;
    
}

vector<SnarlTraversal> RepresentativeTraversalFinder::find_traversals(const Snarl& site) {
    
    // TODO: we can only do ultrabubbles right now. Other snarls may not have
    // traversals through from end to end.
    assert(site.type() == ULTRABUBBLE);
    
    const Snarl* managed_site = snarl_manager.manage(site);
    
    // We may need to make a new index for a backbone for this site, if it's not
    // on the primary path.
    unique_ptr<PathIndex> backbone_index;
    
    // See what the index for the appropriate primary path, if any, is. If we
    // get something non-null the site must be threaded on it.
    PathIndex* primary_path_index = get_index(site);
    
    if (primary_path_index == nullptr ||
        !primary_path_index->by_id.count(site.start().node_id()) ||
        !primary_path_index->by_id.count(site.end().node_id())) {
        // This site is not strung along the primary path, so we will need to
        // generate a backbone traversal of it to structure our search for
        // representative traversals (because they always want to come back to
        // the backbone as soon as possible).
        
        // TODO: we don't handle children correctly (we just glom them into
        // ourselves).
        Path backbone = find_backbone(site);
        
        // Index the backbone (but don't bother with the sequence)
        backbone_index = unique_ptr<PathIndex>(new PathIndex(backbone));
    }
    
    // Determine what path will be the path we use to scaffold the traversals:
    // the primary path index by default, or the backbone index if we needed one.
    PathIndex& index = (backbone_index.get() != nullptr ? *backbone_index : *primary_path_index);
    
    // Get the site's nodes and edges, including our outer boundary nodes, not used inside children.
    // TODO: can we not include the child boundaries? Would that make things easier?
    pair<unordered_set<Node*>, unordered_set<Edge*>> contents = snarl_manager.shallow_contents(&site, augmented.graph, true);
    
    // Copy its node set
    unordered_set<Node*> nodes_left(contents.first);

    // Trace the ref path through the site.
    vector<Visit> ref_path_for_site;
    
    // First figure where the site starts and ends in the selected path
    size_t site_start = index.by_id.at(site.start().node_id()).first;
    size_t site_end = index.by_id.at(site.end().node_id()).first;
    
#ifdef debug
    cerr << "Site starts with " << to_node_traversal(site.start(), augmented.graph)
         << " at " << site_start
         << " and ends with " << to_node_traversal(site.end(), augmented.graph)
         << " at " << site_end << endl;
        
    for (auto* node : nodes_left) {
        cerr << "\tContains node " << node->id() << endl;
    }
#endif

    // The primary path may go through the site backward. So get the primary min and max coords
    size_t primary_min = min(site_start, site_end);
    size_t primary_max = max(site_start, site_end);

    // Then walk nodes from min coordinate to max coordinate. This holds the
    // start coordinate of the current node.
    int64_t ref_node_start = primary_min;
    while(ref_node_start <= primary_max) {
    
        // Find the reference node starting here or later.
        auto found = index.by_start.lower_bound(ref_node_start);
        if(found == index.by_start.end()) {
            throw runtime_error("No backbone node found when tracing through site!");
        }
#ifdef debug
        cerr << "Ref node: " << found->second << " at " << ref_node_start << "/" << primary_max << endl;
#endif
        if((*found).first > primary_max) {
            // The next reference node we can find is out of the space
            // being replaced. We're done.
            if (verbose) {
                cerr << "Stopping for out-of-bounds node" << endl;
            }
            break;
        }
        
        // Get the corresponding Visit
        Visit found_visit = found->second.to_visit();
        
        // What node did we hit?
        Node* visited_node = augmented.graph.get_node(found_visit.node_id());
        
        const Snarl* child = snarl_manager.into_which_snarl(found_visit);
        if (child != nullptr && child != managed_site &&
            snarl_manager.into_which_snarl(reverse(found_visit)) != managed_site &&
            !(eat_trivial_children && snarl_manager.is_trivial(child, augmented.graph))) {
            // If the node in this orientation enters a child, and it's not a
            // trivial child we are taking care of ourselves
        
            // Visit the child
            Visit child_visit;
            *child_visit.mutable_snarl()->mutable_start() = child->start();
            *child_visit.mutable_snarl()->mutable_end() = child->end();
            if (found_visit == child->start()) {
                // We enter the child on its left
                child_visit.set_backward(false);
            } else {
                assert(found_visit == reverse(child->end()));
                // We enter the child on its right
                child_visit.set_backward(true);
            }
            ref_path_for_site.push_back(child_visit);
        
            // And skip to its other end.
            // TODO: the path is not allowed to end inside the snarl.
            Node* here = visited_node;
            do {
#ifdef debug
                cerr << "at node " << pb2json(*here) << endl;
#endif
                // Advance
                ref_node_start = found->first + here->sequence().size();
                // And look at what we get
                found = index.by_start.lower_bound(ref_node_start);
                assert(found != index.by_start.end());
                // And grab out the node
                found_visit = found->second.to_visit();
                here = augmented.graph.get_node(found_visit.node_id());
                // Until we find something in this parent again that isn't the
                // closing visit of a child snarl. We'll look at what we find
                // next.
            } while (!contents.first.count(here));
            
            if (snarl_manager.into_which_snarl(reverse(found_visit)) != nullptr) {
                // We hit the end node of the child snarl.
                
                if (snarl_manager.into_which_snarl(found_visit) == nullptr) {
                    // We don't have another child snarl immediately. Look at the node after this one.
                    ref_node_start = found->first + here->sequence().size();
                    found = index.by_start.lower_bound(ref_node_start);
                    assert(found != index.by_start.end());
                    found_visit = found->second.to_visit();
                    here = augmented.graph.get_node(found_visit.node_id());
                } else {
                    // It's also the start node of another child snarl, so loop
                    // on it again. Do nothing here.
#ifdef debug
                    cerr << "Back-to-back child snarls!" << endl;
#endif
                }
            }
            
            // Make sure we actually found something meeting those criteria.
            // TODO: the path is not allowed to end inside the snarl.
            assert(contents.first.count(here));
        } else {
            // Otherwise, visit this node
            
            if (nodes_left.count(visited_node)) {
                // If the node is one we still expect to see, drop it from the
                // set of nodes in the site we haven't visited.
                nodes_left.erase(visited_node);
            }
            
            // Add the traversal to the ref path through the site
            ref_path_for_site.push_back(found_visit);
            
            
            // Next iteration look where this node ends.
            ref_node_start = found->first + visited_node->sequence().size();
        }
        
#ifdef debug
        cerr << "Added visit: " << pb2json(ref_path_for_site.back()) << endl; 
#endif   
    }
    
    // We leave the ref path in backbone-relative forward orientation, because
    // all our bubbles we find will also be in backbone-relative forward
    // orientation.
    
    for(auto node : nodes_left) {
        // Make sure none of the nodes in the site that we didn't visit
        // while tracing along the ref path are on the ref path.
        
        if (snarl_manager.into_which_snarl(node->id(), true) || snarl_manager.into_which_snarl(node->id(), false)) {
            // Skip child boundary nodes.
            continue;
        }
        
        if(index.by_id.count(node->id())) {
            cerr << "error[RepresentativeTraversalFinder]: Node " << node->id() << " is on backbone path at "
                 << index.by_id.at(node->id()).first << " but not traced in site "
                 << to_node_traversal(site.start(), augmented.graph) << " to " 
                 << to_node_traversal(site.end(), augmented.graph) << " that contains it." << endl;
            cerr << "error[RepresentativeTraversalFinder]: This can happen when the path you are calling "
                 << "against traverses the same part of your graph twice." << endl;
            throw runtime_error("Extra ref node found");
        }
    }
    
    // We need to know all the full-length traversals we're going to consider.
    // XREF states will have to be calculated later, over the whole traversal.
    set<vector<Visit>> site_traversal_set;
    
    // We have this function to extend a partial traversal into a full
    // traversal and add it as a path. The path must already be rooted on
    // the reference in the correct order and orientation.
    auto extend_into_allele = [&](vector<Visit> path) {
        // Splice the ref path through the site and the bubble's path
        // through the site together.
        vector<Visit> extended_path;
#ifdef debug
        cerr << "Input path: " << endl;
        for(auto& visit : path) {
            if(visit.node_id() != 0 && index.by_id.count(visit.node_id())) {
                cerr << "\tPath member " << visit << " lives on backbone at "
                     << index.by_id.at(visit.node_id()).first << endl;
            } else {
                cerr << "\tPath member " << visit << " does not live on backbone" << endl;
            }
        }
#endif
        
        for(auto& visit : path) {
            if (visit.node_id() != 0) {
                // Make sure the site actually has the nodes we're visiting.
                assert(contents.first.count(augmented.graph.get_node(visit.node_id())));
            }
            // Child snarls will have ownership of their end nodes, so they won't be part of our contents.
        }
        
        size_t ref_path_index = 0;
        size_t ref_path_end = ref_path_for_site.size() - 1;
        size_t bubble_path_index = 0;
        size_t bubble_path_end = path.size() - 1;
        
        // Function that returns the node visit on the on one end of a visit, even
        // if the visit is of a child snarl
        auto frontier_visit = [&](const Visit& visit, bool left_side) {
            if (visit.node_id() != 0) {
                return visit;
            }
            else if (visit.backward() && !left_side) {
                return reverse(visit.snarl().start());
            }
            else if (!visit.backward() && !left_side) {
                return visit.snarl().end();
            }
            else if (visit.backward()) {
                return reverse(visit.snarl().end());
            }
            else {
                return visit.snarl().start();
            }
        };
        
        while(frontier_visit(ref_path_for_site[ref_path_index], false) != frontier_visit(path[bubble_path_index], true) &&
              !(path[bubble_path_index].node_id() == 0 &&
                frontier_visit(ref_path_for_site[ref_path_index], false) == frontier_visit(path[bubble_path_index], false))) {
            // Collect NodeTraversals from the ref path until we hit the one
            // at which the bubble path starts.
#ifdef debug
            cerr << "Before path: " << pb2json(ref_path_for_site[ref_path_index]) << endl;
#endif
            extended_path.push_back(ref_path_for_site[ref_path_index++]);
        }
        
        if (ref_path_for_site[ref_path_index].node_id() == 0) {
            // The last Visit we traversed from the ref was a Snarl, so it already
            // includes the first node of the path as one of its boundaries. We need
            // to add the ref visit and exclude the bubble visit unless it is also of
            // a child Snarl and that Snarl is different from the ref path
            
#ifdef debug
            cerr << "Adding final ref child visit " << pb2json(ref_path_for_site[ref_path_index]) << endl;
#endif
            extended_path.push_back(ref_path_for_site[ref_path_index]);
            
            if (path.front().node_id() != 0 || (path.front().snarl().start() == extended_path.back().snarl().start()
                                                && path.front().snarl().end() == extended_path.back().snarl().end())) {
#ifdef debug
                cerr << "Skipping bubble visit " << pb2json(path[bubble_path_index]) << endl;
#endif
                bubble_path_index++;
            }
        }
        
        while(bubble_path_index < path.size()) {
            // Then take the whole bubble path
#ifdef debug
            cerr << "In path: " << pb2json(path[bubble_path_index]) << endl;
#endif
            extended_path.push_back(path[bubble_path_index++]);
        }
        while(ref_path_index < ref_path_for_site.size()) {
            // Then skip ahead to the matching point in the ref path, which may
            // be either a full visit match, or a match to the exit node of a
            // visited child snarl.
            
            // Check each reference visit
            auto& ref_visit = ref_path_for_site.at(ref_path_index);
            
#ifdef debug
            cerr << "At ref: " << pb2json(ref_visit) << " with frontier visit " << pb2json(frontier_visit(ref_visit, true)) << " looking for " << pb2json(frontier_visit(path.back(), false)) << endl;
#endif
            
            if (frontier_visit(ref_visit, true) == frontier_visit(path.back(), false) ||
                (path.back().node_id() == 0 && frontier_visit(ref_visit, false) == frontier_visit(path.back(), false))) {
                // We found the exit node. Put the after-the-bubble visits
                // continuing from here.
                break;
            }
            
            // Otherwise this ref visit isn't the right one to match up with our
            // bubble's traversal.
#ifdef debug
            cerr << "Skip ref: " << pb2json(ref_path_for_site[ref_path_index]) << endl;
            cerr << "\tWant: " << pb2json(path.back()) << endl;
#endif
            ref_path_index++;
        }
        
        if(ref_path_index == ref_path_for_site.size()) {
            // We ran out of ref path before finding the place to leave the alt.
            // This must be a backtracking loop or something; start over from the beginning.
            cerr << "ran out of ref" << endl;
            ref_path_index = 0;
            
            while(ref_path_index < ref_path_for_site.size() &&
                  frontier_visit(ref_path_for_site.at(ref_path_index), true) != frontier_visit(path.back(), false)) {
                // Then skip ahead to the matching point in the ref path
                ref_path_index++;
            }
            
            if(ref_path_index == ref_path_for_site.size()) {
                // Still couldn't find it!
                stringstream err;
                err << "Couldn't find " << path.back() << " in backbone path of site "
                << site.start()
                << " to " << site.end() << endl;
                throw runtime_error(err.str());
            }
        }
        
        if (ref_path_for_site[ref_path_index].node_id() == 0) {
            // The next Visit on the ref path is to a Snarl, so the final Visit we added
            // from the bubble will be redundant with its boundary nodes unless that Visit
            // was also to a Snarl
            if (extended_path.back().node_id() != 0 ||
                (ref_path_for_site[ref_path_index].snarl().start() == extended_path.back().snarl().start()
                 && ref_path_for_site[ref_path_index].snarl().end() == extended_path.back().snarl().end())) {
#ifdef debug
                cerr << "Removing bubble visit " << pb2json(extended_path.back()) << endl;
#endif
                extended_path.pop_back();
            }
#ifdef debug
            cerr << "Adding adjacent ref child visit" << pb2json(ref_path_for_site[ref_path_index]) << endl;
#endif
            extended_path.push_back(ref_path_for_site[ref_path_index]);
        }
        // Skip the matching NodeTraversal
        ref_path_index++;
        while(ref_path_index < ref_path_for_site.size()) {
            // Then take the entier rest of the ref path
#ifdef debug
            cerr << "After path: " << pb2json(ref_path_for_site[ref_path_index]) << endl;
#endif
            extended_path.push_back(ref_path_for_site[ref_path_index++]);
        }
        
#ifdef debug
        cerr << "Output path: " << endl;
        for (auto& v : extended_path) {
            cerr << "\t" << pb2json(v) << endl;
        }
#endif
        
        // Now add it to the set
        site_traversal_set.insert(extended_path);
        
    };

#ifdef debug
    cerr << "Explore " << contents.first.size() << " nodes" << endl;
#endif

    for (Node* node : contents.first) {
        // Find the bubble for each node
        
        if (snarl_manager.into_which_snarl(node->id(), true) || snarl_manager.into_which_snarl(node->id(), false)) {
            // Don't start from nodes that are child boundaries
            continue;
        }
        
        if (augmented.has_supports() && total(augmented.get_support(node)) == 0) {
            // Don't bother with unsupported nodes
            continue;
        }
        
        if (index.by_id.count(node->id())) {
            // Don't try to pathfind to the backbone for backbone nodes.
            continue;
        }
        
#ifdef debug
        cerr << "Base path on " << node->id() << endl;
#endif
        
        // Find bubbles that backend into the backbone path
        pair<Support, vector<Visit>> sup_path = find_bubble(node, nullptr, nullptr, index, site);

        vector<Visit>& path = sup_path.second;
        
        if(path.empty()) {
            // We couldn't find a path back to the primary path. Discard
            // this material.
            if (verbose) {
                cerr << "Warning: No path found for node " << node->id() << endl;
            }
            // TODO: record the node's bases as lost.
            
            // TODO: what if it's already in another bubble/the node is deleted?
            continue;
        }
        
        // Extend it out into an allele
        extend_into_allele(path);
        
    }
    

#ifdef debug
    cerr << "Explore " << contents.second.size() << " edges" << endl;
#endif

    for(Edge* edge : contents.second) {
        // Go through all the edges
        
        if(augmented.has_supports() && total(augmented.get_support(edge)) == 0) {
            // Don't bother with unsupported edges
#ifdef debug
            cerr << "Skip unsupported edge " << edge->from() << " -> " << edge->to() << endl;
#endif
            continue;
        }
        
        if(!index.by_id.count(edge->from()) || !index.by_id.count(edge->to())) {
            // Edge doesn't touch backbone at both ends. Don't use it
            // because for some reason it makes performance worse
            // overall.
#ifdef debug
            cerr << "Skip off-backbone edge " << edge->from() << " -> " << edge->to() << endl;
#endif
            continue;
        }
        
#ifdef debug
        cerr << "Base path on " << edge->from() << " -> " << edge->to() << endl;
#endif
        
        // Find a path based around this edge
        pair<Support, vector<Visit>> sup_path = find_bubble(nullptr, edge, nullptr, index, site);
        vector<Visit>& path = sup_path.second;
        
#ifdef debug
        cerr << "Edge " << edge->from() << " to " << edge->to() << " yields:" << endl;
        for(auto& visit : path) {
            cerr << "\t" << visit << endl;
        }
#endif
        
        if(path.empty()) {
            // We couldn't find a path back to the primary path. Discard
            // this material.
            if (verbose) {
                cerr << "Warning: No path found for edge " << edge->from() << "," << edge->to() << endl;
            }
            // TODO: bases lost
            continue;
        }
        
        // Extend it out into an allele
        extend_into_allele(path);
    }
    
    auto children = snarl_manager.children_of(&site);

#ifdef debug
    cerr << "Explore " << children.size() << " children" << endl;
#endif

    for (const Snarl* child : children) {
        // Go through all the child snarls

        if (eat_trivial_children && snarl_manager.is_trivial(child, augmented.graph)) {
            // Skip trivial children
            continue;
        }
        
#ifdef debug
        cerr << "Base path on " << *child << endl;
#endif
        
        // Find a path based around this child snarl
        pair<Support, vector<Visit>> sup_path = find_bubble(nullptr, nullptr, child, index, site);
        vector<Visit>& path = sup_path.second;
        
        if(path.empty()) {
            // We couldn't find a path back to the primary path. Discard
            // this material.
            if (verbose) {
                cerr << "Warning: No path found for child snarl " << *child << endl;
            }
            // TODO: bases lost
            continue;
        }
        
        // Extend it out into an allele
        extend_into_allele(path);
    }
    
    
    // Now convert to SnarlTraversals
    vector<SnarlTraversal> unique_traversals;
    
    // Have a function to convert a vector of NodeTraversals including the snarl
    // ends into a SnarlTraversal
    auto emit_traversal = [&](vector<Visit> visits) {
        // Make it into this SnarlTraversal
        SnarlTraversal trav;
        
#ifdef debug
        cerr << "Unique traversal's visits:" << endl;
        for(auto& visit : visits) {
            cerr << "\t" << visit << endl;
        }
#endif
        
        // Convert the path to Visits.
        // TODO: think about nested sites?
        
        if (site_start > site_end) {
            // The primary path runs backward, so we need to emit all our lists
            // of NodeTraversals backward so they come out in the snarl's
            // orientation and not the primary path's.
            
            for(size_t i = 0; i < visits.size(); i++) {
                // Record a Visit for each Visit but the first and last,
                // but backward and in reverse order.
                *trav.add_visit() = reverse(visits[visits.size() - i - 1]);
            }
            
        } else {
            // The primary path and the snarl use the same orientation
            
            for(size_t i = 0; i < visits.size(); i++) {
                // Make a Visit for each NodeTraversal but the first and last
                *trav.add_visit() = visits[i];
            }
            
        }
        
#ifdef debug
        cerr << "Unique traversal: " << pb2json(trav) << endl;
#endif

        // Save the SnarlTraversal
        unique_traversals.push_back(trav);
    };
    
    
    // Do the ref path first
    emit_traversal(ref_path_for_site);
    for(auto& visits : site_traversal_set) {
        // Look at each vector of Visits
        if (visits != ref_path_for_site) {
            // And do everything other than the ref path
            emit_traversal(visits);            
        }
        
    }
    
    return unique_traversals;
}

pair<Support, vector<Visit>> RepresentativeTraversalFinder::find_bubble(Node* node, Edge* edge,
                                                                        const Snarl* snarl, PathIndex& index, const Snarl& site) {

    // What are we going to find our left and right path halves based on?
    Visit left_visit;
    Visit right_visit;
    
    const Snarl* managed_site = snarl_manager.manage(site);
    
    if (edge != nullptr) {
        // Be edge-based
        
        // Find the nodes at the ends of the edges. Look at them traversed in the
        // edge's local orientation.
        left_visit = to_visit(edge->from(), edge->from_start());
        right_visit = to_visit(edge->to(), edge->to_end());
        
        const Snarl* right_child = snarl_manager.into_which_snarl(right_visit);
        const Snarl* left_child = snarl_manager.into_which_snarl(left_visit);
        
        if (right_child != nullptr && right_child != managed_site
            && snarl_manager.into_which_snarl(reverse(right_visit)) != managed_site) {
            // We're reading into a child snarl on the right.
#ifdef debug
            cerr << "Child to right of edge " << pb2json(*right_child) << endl;
#endif
            
            // Make a visit.
            Visit right_child_visit = to_visit(*right_child);
            
            if (to_left_side(right_visit) != to_left_side(right_child_visit)) {
                // We really go through it the other way around
                right_child_visit = reverse(right_child_visit);
            }
            assert(to_left_side(right_visit) == to_left_side(right_child_visit));
            // Use the snarl visit instead of the visit to the entering node.
            right_visit = right_child_visit;
        }
        
        if (left_child != nullptr && left_child != managed_site
            && snarl_manager.into_which_snarl(reverse(left_visit)) != managed_site) {
            // We're reading out of a child snarl on the left.
#ifdef debug
            cerr << "Child to left of edge " << pb2json(*left_child) << endl;
#endif
            
            // Make a visit.
            Visit left_child_visit = to_visit(*left_child);
            
            if (to_right_side(left_visit) != to_right_side(left_child_visit)) {
                // We really go through it the other way around
                left_child_visit = reverse(left_child_visit);
            }
            assert(to_right_side(left_visit) == to_right_side(left_child_visit));
            // Use the snarl visit instead of the visit to the exiting node.
            left_visit = left_child_visit;
        }
#ifdef debug
        cerr << "Edge becomes " << left_visit << " -> " << right_visit << endl;
#endif
        
    } else if (node != nullptr) {
        // Be node-based. TODO: we trust the caller not to feed us nodes that
        // are part of/boundaries of child snarls.
        left_visit = right_visit = to_visit(node->id(), false);
    } else {
        // Be snarl-based
        assert(snarl != nullptr);
        left_visit = right_visit = to_visit(*snarl);
    }
    
#ifdef debug
    cerr << "Starting from: " << left_visit << ", " << right_visit << endl;
#endif

    // Find paths on both sides, with nodes or snarls on the primary path at the
    // outsides and this visit in the middle. Returns path lengths and paths in
    // pairs in a set.
    auto leftPaths = bfs_left(left_visit, index, false, managed_site);
    auto rightPaths = bfs_right(right_visit, index, false, managed_site);
    
    // Find a combination of two paths which gets us to the reference in a
    // consistent orientation (meaning that when you look at the ending nodes'
    // Mappings in the reference path, the ones with minimal ranks have the same
    // orientations) and which doesn't use the same nodes on both sides.
    // Track support of up to max_bubble_paths combinations, and return the
    // highest
    pair<Support, vector<Visit> > bestBubblePath;
    int bubbleCount = 0;
    
    // We need to look in different combinations of lists.
    auto testCombinations = [&](const list<list<Visit>>& leftList,
                                const list<list<Visit>>& rightList) {
        
#ifdef debug        
        cerr << "Combine " << leftList.size() << " left sides and "
        << rightList.size() << " right sides" << endl;
#endif
        
        // We know the left list starts and the right list ends with an actual
        // node visit, if only to the snarl's start or end.

        for(auto leftPath : leftList) {
            // Figure out the relative orientation for the leftmost node.
#ifdef debug        
            cerr << "Left path: " << endl;
            for(auto visit : leftPath ) {
                cerr << "\t" << visit << endl;
            }
#endif    
            
            // Find what node side actually represents the end
            auto leftSide = to_left_side(leftPath.front());
            
            // Split out the node's orientation
            bool leftOrientation = leftSide.is_end;
            
            // Get where it falls in the reference as a position, orientation pair.
            auto leftRefPos = index.by_id.at(leftSide.node);
            
            // We have a backward orientation relative to the reference path if we
            // were traversing the anchoring node backwards, xor if it is backwards
            // in the reference path.
            bool leftRelativeOrientation = leftOrientation != leftRefPos.second;
            
            // Make a set of all the nodes in the left path
            set<int64_t> leftPathNodes;
            // And one of all the snarls (with bounding visits set)
            set<Snarl> leftPathSnarls;
            for(auto visit : leftPath) {
                if (visit.node_id() != 0) {
                    // It's a node visit
                    leftPathNodes.insert(visit.node_id());
                } else {
                    // It's a snarl visit
                    leftPathSnarls.insert(visit.snarl());
                }
            }

            // Get the minimum support in the left path
            Support minLeftSupport = min_support_in_path(leftPath);
            
            for(auto rightPath : rightList) {
                // Figure out the relative orientation for the rightmost node.
#ifdef debug            
                cerr << "Right path: " << endl;
                for(auto visit : rightPath ) {
                    cerr << "\t" << visit << endl;
                }
#endif            
                
                // Find what node side actually represents the end
                // Remember it's at the end of this path.
                auto rightSide = to_right_side(rightPath.back());
                
                // Split out the node's orientation
                bool rightOrientation = !rightSide.is_end;
                
                // Get where it falls in the reference as a position, orientation pair.
                auto rightRefPos = index.by_id.at(rightSide.node);
                
                // We have a backward orientation relative to the reference path if we
                // were traversing the anchoring node backwards, xor if it is backwards
                // in the reference path.
                bool rightRelativeOrientation = rightOrientation != rightRefPos.second;

                // Get the minimum support in the right path
                Support minRightSupport = min_support_in_path(rightPath);
                
                if(leftRelativeOrientation == rightRelativeOrientation &&
                   ((!leftRelativeOrientation && leftRefPos.first < rightRefPos.first) ||
                    (leftRelativeOrientation && leftRefPos.first > rightRefPos.first))) {
                    // We found a pair of paths that get us to and from the
                    // reference without turning around, and that don't go back to
                    // the reference before they leave.

                    // Get the minimum support of combined left and right paths
                    Support minFullSupport = support_min(minLeftSupport, minRightSupport);
                    
                    // Start with the left path
                    vector<Visit> fullPath{leftPath.begin(), leftPath.end()};
                    
                    // We need to detect overlap with the left path
                    bool overlap = false;
                    
                    // If we're starting from an edge, we keep the first visit
                    // on the right path. If we're starting from a node or
                    // snarl, we need to discard it because it's just another
                    // copy of our visit we're starting with.
                    for(auto it = (edge != nullptr ? rightPath.begin() : ++(rightPath.begin())); it != rightPath.end(); ++it) {
                        // For all but the first node on the right path, add that in
                        fullPath.push_back(*it);
                        
                        if (it->node_id() != 0) {
                            // This right-side visit hits a node
                            if(leftPathNodes.count(it->node_id())) {
                                // We already visited this node on the left side. Try
                                // the next right path instead.
                                overlap = true;
                            }
                        } else {
                            // This right-side visit hits a snarl
                            if(leftPathSnarls.count(it->snarl())) {
                                // We already visited this snarl on the left side. Try
                                // the next right path instead.
                                overlap = true;
                            }
                        }
                    }
                    
                    if(overlap) {
                        // Can't combine this right with this left, as they
                        // share nodes or child snarls and we can't handle the
                        // copy number implications. Try the next right. TODO:
                        // handle the copy number implications.
                        // TODO: This shouldn't happen in ultrabubbles.
                        continue;
                    }
                    
                    if(leftRelativeOrientation) {
                        // Turns out our anchored path is backwards.
                        
                        // Reorder everything the other way
                        reverse(fullPath.begin(), fullPath.end());
                        
                        for(auto& visit : fullPath) {
                            // Flip each Visit
                            visit = reverse(visit);
                        }
                    }

                    
#ifdef debug        
                    cerr << "Merged path:" << endl;
                    for(auto visit : fullPath) {
                        cerr << "\t" << visit << endl;
                    }                    
#endif
                    // Update our best path by seeing if we've found one with
                    // higher min support. Make sure to replace the empty path
                    // even if we only find a traversal with 0 support (because
                    // maybe we have no support data at all).
                    if (total(minFullSupport) > total(bestBubblePath.first) ||
                        (total(minFullSupport) == total(bestBubblePath.first) && bestBubblePath.second.empty())) {
                        bestBubblePath.first = minFullSupport;
                        bestBubblePath.second = fullPath;
                    }

                    // keep things from getting out of hand
                    if (++bubbleCount >= max_bubble_paths) {
                        return bestBubblePath;
                    }
                }
            }
        }
        
        // Return the best path along with its min support
        // (could be empty)
        return bestBubblePath;
        
    };
    
    // Convert sets to lists, which requires a copy again...
    // TODO: Can we just completely remove the length calculation?
    list<list<Visit>> leftConverted;
    for(auto lengthAndPath : leftPaths) {
        leftConverted.emplace_back(move(lengthAndPath.second));
    }
    list<list<Visit>> rightConverted;
    for(auto lengthAndPath : rightPaths) {
        rightConverted.emplace_back(move(lengthAndPath.second));
    }
    
    // Look for a valid combination, or return an empty path if one iesn't
    // found.
    return testCombinations(leftConverted, rightConverted);
    
}

Support RepresentativeTraversalFinder::min_support_in_path(const list<Visit>& path) {
    
    if (path.empty()) {
        // No support if we visit nothing!
        return Support();
    }
    
    // Get an iterator to the current visit on the path
    auto cur = path.begin();
    // And to the next visit on the path
    auto next = path.begin();
    ++next;
    
    // Have we found anything with a support yet?
    bool supportFound = false;
    // If we have, this holds the min support we have found.
    Support minSupport;
    
    if (cur->node_id() != 0) {
        // We're at a node visit, so we have a support to start with
        minSupport = augmented.get_support(augmented.graph.get_node(cur->node_id()));
        supportFound = true;
    }
    
    for (; next != path.end(); ++cur, ++next) {
        // For each visit and its next visit
    
        if (next->node_id() != 0) {
            // The next visit is to a node, so get its support
            Support nextSupport = augmented.get_support(augmented.graph.get_node(next->node_id()));
            
            if (supportFound) {
                // Min it against existing support
                minSupport = support_min(minSupport, nextSupport);
            } else {
                // Take as the found support
                minSupport = nextSupport;
                supportFound = true;
            }
        }
        
        // TODO: Support for child snarls!
    
        // check the edge support
        Edge* edge = augmented.graph.get_edge(to_right_side(*cur), to_left_side(*next));
        
        if (edge != nullptr) {
            // The edge exists (because we aren't back-to-back child snarls)
            Support edgeSupport = augmented.get_support(edge);
            
            if (supportFound) {
                // Min it against existing support
                minSupport = support_min(minSupport, edgeSupport);
            } else {
                // Take as the found support
                minSupport = edgeSupport;
                supportFound = true;
            }
        }
    }

    // This may be 0 if we hit no nodes or edges, but I guess that's OK...
    return minSupport;
}

set<pair<size_t, list<Visit>>> RepresentativeTraversalFinder::bfs_left(Visit visit,
                                                                       PathIndex& index, bool stopIfVisited, const Snarl* in_snarl) {

    // Holds partial paths we want to return, with their lengths in bp.
    set<pair<size_t, list<Visit>>> toReturn;
    
    // Do a BFS
    
    // This holds the paths to get to NodeTraversals to visit (all of which will
    // end with the node we're starting with).
    list<list<Visit>> toExtend;
    
    // This keeps a set of all the oriented nodes we already got to and don't
    // need to queue again.
    set<Visit> alreadyQueued;
    
    // Start at this node at depth 0
    toExtend.emplace_back(list<Visit> {visit});
    // Mark this traversal as already queued
    alreadyQueued.insert(visit);
    
    // How many ticks have we spent searching?
    size_t searchTicks = 0;
    
#ifdef debug
    cerr << "Start BFS" << endl;
#endif

    // Track how many options we have because size may be O(n).
    size_t stillToExtend = toExtend.size();
    
    while (!toExtend.empty()) {
        // Keep going until we've visited every node up to our max search depth.
        
        searchTicks++;
        

#ifdef debug
        // Report on how much searching we are doing.
        cerr << "Search tick " << searchTicks << ", " << stillToExtend << " options." << endl;
#endif

        
        // Dequeue a path to extend.
        // Make sure to move out of the list to avoid a useless copy.
        list<Visit> path(move(toExtend.front()));
        toExtend.pop_front();
        stillToExtend--;
        
        // We can't just throw out longer paths, because shorter paths may need
        // to visit a node twice (in opposite orientations) and thus might get
        // rejected later. Or they might overlap with paths on the other side.
        
        // Look up and see if the front node on the path is on our reference
        // path
        if (path.front().node_id() != 0 && index.by_id.count(path.front().node_id())) {
            // This visit is to a node, which is on the reference path.
            
#ifdef debug
            cerr << "Reached anchoring node " << path.front().node_id() << endl;
            cerr << "Emit path of length " << path.size() << endl;
#endif
            
            // Say we got to the right place
            toReturn.emplace(bp_length(path), move(path));
            
            // Don't bother looking for extensions, we already got there.
        } else if (path.front().node_id() == 0 && !path.front().backward() &&
                   index.by_id.count(path.front().snarl().start().node_id())) {
            // This visit is to a snarl, which is on the reference path on its
            // left end.
            
#ifdef debug
            cerr << "Reached start of anchoring snarl " << path.front().snarl() << endl;
#endif
            
            // Say we got to the right place
            toReturn.emplace(bp_length(path), move(path));
            
            // Don't bother looking for extensions, we already got there.
        } else if (path.front().node_id() == 0 && path.front().backward() &&
                   index.by_id.count(path.front().snarl().end().node_id())) {
            // This visit is to a snarl in reverse, which is on the reference
            // path on its right end.
            
#ifdef debug
            cerr << "Reached end of anchoring snarl " << path.front().snarl() << endl;
#endif
            
            // Say we got to the right place
            toReturn.emplace(bp_length(path), move(path));
            
            // Don't bother looking for extensions, we already got there.
        } else if (path.size() <= max_depth) {
            // We haven't hit the reference path yet, but we also haven't hit
            // the max depth. Extend with all the possible extensions.
            
            // Look left, possibly entering child snarls
            vector<Visit> prevVisits = snarl_manager.visits_left(path.front(), augmented.graph, in_snarl);
            
#ifdef debug
            cerr << "Consider " << prevVisits.size() << " prev visits" << endl;
            for (Visit vis : prevVisits) {
                cerr << "\t" << vis << endl;
            }
#endif
            
            for (auto prevVisit : prevVisits) {
                // For each node we can get to
                
                if (prevVisit.node_id() != 0) {
                    // This is a visit to a node
                    
                    // Make sure the edge is real, since it can't be a back-to-
                    // back site
                    Edge* edge = augmented.graph.get_edge(to_right_side(prevVisit), to_left_side(path.front()));
                    assert(edge != NULL);
                
                    // Fetch the actual node
                    Node* prevNode = augmented.graph.get_node(prevVisit.node_id());
                    
                    if (augmented.has_supports() && 
                        (total(augmented.get_support(prevNode)) == 0 || total(augmented.get_support(edge)) == 0)) {
                        // We have no support at all for visiting this node by this
                        // edge (but we do have some read support data)
                        
#ifdef debug
                        cerr << "Reject " << prevNode->id() << " with no support" << endl;
#endif
                        
                        continue;
                    }
                } else {
                    // This is a visit to a child snarl
                    
                    // Look at the node we would leave the child snarl on
                    // That node can't be shared with a snarl we are already at.
                    Node* prevNode = augmented.graph.get_node(to_left_side(prevVisit).node);
                    
                    if (augmented.has_supports() && total(augmented.get_support(prevNode)) == 0) {
                        // We have no support at all for visiting the far node of this snarl
                        
#ifdef debug
                        cerr << "Reject " << prevVisit.snarl() << " with no support on far node" << endl;
#endif
                        
                        continue;
                    }
                    
                    // TODO: when snarls are not back-to-back, check the connecting edges
                }
                
                
                if (stopIfVisited && alreadyQueued.count(prevVisit)) {
                    // We already have a way to get here.
                    
#ifdef debug
                    cerr << "Reject " << prevVisit << " which is already seen" << endl;
#endif
                    
                    continue;
                }
                
                if (stillToExtend >= max_width) {
                    // Don't consider this extension because we already have too
                    // many. But since every time we go around the outer loop
                    // here we pop something off the queue to consider extending
                    // it, we'll always have at least one slot left to fill with
                    // an extension, so we'll at least keep exploring one path.

                
#ifdef debug
                    cerr << "Reject " << prevVisit << " because queue is full" << endl;
#endif

                    
                    continue;
                }
                
#ifdef debug
                cerr << "Accept " << prevVisit << endl;
#endif
            
                // Make a new path extended left with the node
                list<Visit> extended(path);
                extended.push_front(prevVisit);
                toExtend.emplace_back(move(extended));
                stillToExtend++;
                
                // Remember we found a way to this node, so we don't try and
                // visit it other ways.
                alreadyQueued.insert(prevVisit);
            }
        } else if (path.size() >= max_depth) {
#ifdef debug
            cerr << "Path has reached max depth! Aborting!" << endl;
#endif
        } else {
            // We should have handled all the possibilities.
            assert(false);
        }
        
    }
    
    return toReturn;
}

set<pair<size_t, list<Visit>>> RepresentativeTraversalFinder::bfs_right(Visit visit, PathIndex& index, bool stopIfVisited,
                                                                        const Snarl* in_snarl) {

    // Look left from the backward version of the visit.
    auto toConvert = bfs_left(reverse(visit), index, stopIfVisited, in_snarl);
    
    // Since we can't modify set records in place, we need to do a copy
    set<pair<size_t, list<Visit>>> toReturn;
    
    for(auto lengthAndPath : toConvert) {
        // Flip every path to run the other way
        lengthAndPath.second.reverse();
        for(auto& v : lengthAndPath.second) {
            // And invert the orientation of every visit in the path in place.
            v = reverse(v);
        }
        // Stick it in the new set
        toReturn.emplace(move(lengthAndPath));
    }
    
    return toReturn;
}

size_t RepresentativeTraversalFinder::bp_length(const list<Visit>& path) {
    size_t length = 0;
    for(auto& visit : path) {
        // Sum up length of each node's sequence
        if (visit.node_id() != 0) {
            length += augmented.graph.get_node(visit.node_id())->sequence().size();
        }
        // TODO: handle nested sites
    }
    return length;
}


}
