#include "traversal_finder.hpp"
#include "genotypekit.hpp"
#include "algorithms/k_widest_paths.hpp"
#include "cactus.hpp"
#include "gbwt_helper.hpp"
#include "haplotype_extracter.hpp"
//#define debug

namespace vg {

using namespace std;

PathBasedTraversalFinder::PathBasedTraversalFinder(const PathHandleGraph& g, SnarlManager& sm) : graph(g), snarlmanager(sm){
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
    pair<unordered_set<id_t>, unordered_set<edge_t> > contents = snarlmanager.shallow_contents(&site, graph, true);


    // Get the variant paths at the snarl nodes.
    set<string> var_path_names;
    regex front ("(_alt_)(.*)");
    regex alt_str ("(_alt_)");
    regex back ("(_[0-9]*)");
    set<string> gpath_names;
    graph.for_each_path_handle([&](const path_handle_t& path_handle) {
        gpath_names.insert(graph.get_path_name(path_handle));
    });
    
    map<string, set<string> > basename_to_pathnames;
    map<string, bool> path_processed;


    // Collect our paths which cross our snarl's nodes.
    for (id_t node_id : contents.first){
        handle_t node = graph.get_handle(node_id);
        //cerr << "Processing node " << id << endl;
        set<string> p_of_n;
        graph.for_each_step_on_handle(node, [&](const step_handle_t& step) {
            p_of_n.insert(graph.get_path_name(graph.get_path_handle_of_step(step)));
        });

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
                for (auto h : graph.scan_path(graph.get_path_handle(a))) {
                    int64_t n_id = graph.get_id(h);
                    bool backward = graph.get_is_reverse(h);
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
   
ExhaustiveTraversalFinder::ExhaustiveTraversalFinder(const HandleGraph& graph, SnarlManager& snarl_manager,
                                                     bool include_reversing_traversals) :
    graph(graph), snarl_manager(snarl_manager),
    include_reversing_traversals(include_reversing_traversals) {
    // nothing more to do
}
    
ExhaustiveTraversalFinder::~ExhaustiveTraversalFinder() {
    // no heap objects
}

void ExhaustiveTraversalFinder::stack_up_valid_walks(handle_t walk_head, vector<Visit>& stack) {
    
    id_t head_id = graph.get_id(walk_head);

    // get all edges involving this node so we can filter them down to valid walks
    graph.follow_edges(walk_head, false, [&](const handle_t next_node) {
            if (visit_next_node(next_node)) {
                Visit next_visit;
                next_visit.set_node_id(graph.get_id(next_node));
                next_visit.set_backward(graph.get_is_reverse(next_node));
                stack.push_back(next_visit);
            }
        });
}

void ExhaustiveTraversalFinder::add_traversals(vector<SnarlTraversal>& traversals,
                                               handle_t traversal_start,
                                               unordered_set<handle_t>& stop_at,
                                               unordered_set<handle_t>& yield_at) {
    // keeps track of the walk of the DFS traversal
    list<Visit> path;
    
    // these mark the start of the edges out of the node that is on the head of the path
    // they can be used to see how many nodes we need to peel off the path when we're
    // backtracking
    Visit stack_sentinel;
    
    // initialize stack for DFS traversal of site
    vector<Visit> stack{to_visit(graph, traversal_start)};
    
    while (stack.size()) {
        
        Visit node_visit = stack.back();
        stack.pop_back();
        
        // we have traversed all of edges out of the head of the path, so we can pop it off
        if (node_visit == stack_sentinel) {
            path.pop_back();
            continue;
        }
        
        // have we finished a traversal through the site?
        handle_t node_handle = graph.get_handle(node_visit.node_id(), node_visit.backward());
        if (stop_at.count(node_handle)) {
            if (yield_at.count(node_handle)) {
                // yield path as a snarl traversal
                traversals.emplace_back();
                
                // record the traversal in the return value
                for (auto iter = path.begin(); iter != path.end(); iter++) {
                    *traversals.back().add_visit() = *iter;
                }
                // add the final visit
                *traversals.back().add_visit() = node_visit;
            }
            
            // don't proceed to add more onto the DFS stack
            continue;
        }
        
        // mark the beginning of this node's edges forward in the stack
        stack.push_back(stack_sentinel);
        
        // make a visit through the node traversal and add it to the path
        path.emplace_back();
        path.back().set_node_id(node_visit.node_id());
        path.back().set_backward(node_visit.backward());
        
        // does this traversal point into a child snarl?
        const Snarl* into_snarl = snarl_manager.into_which_snarl(node_visit.node_id(),
                                                                 node_visit.backward());
                                                                 
#ifdef debug
        cerr << "Traversal " << node_visit.node_id() << " " << node_visit.backward() << " enters";
        if (into_snarl != nullptr) {
            cerr << " " << pb2json(*into_snarl) << endl;
        } else {
            cerr << " NULL" << endl; 
        }
#endif
                                                                 
        if (into_snarl && !(node_handle == traversal_start)) {
            // add a visit for the child snarl
            path.emplace_back();
            *path.back().mutable_snarl()->mutable_start() = into_snarl->start();
            *path.back().mutable_snarl()->mutable_end() = into_snarl->end();
            
            // mark the beginning of this child snarls edges forward in the stack
            stack.push_back(stack_sentinel);
            
            // which side of the snarl does the traversal point into?
            if (into_snarl->start().node_id() == node_visit.node_id()
                && into_snarl->start().backward() == node_visit.backward()) {
                // Into the start
#ifdef debug
                cerr << "Entered child through its start" << endl;
#endif
                if (into_snarl->start_end_reachable()) {
                    // skip to the other side and proceed in the orientation that the end node takes.
                    stack.push_back(into_snarl->end());
#ifdef debug
                    cerr << "Stack up " << stack.back().node_id() << " " << stack.back().backward() << endl;
#endif
                }
                
                // if the same side is also reachable, add it to the stack too
                if (into_snarl->start_self_reachable()) {
                    // Make sure to flip it around so we come out of the snarl instead of going in again,
                    Visit rev_visit = into_snarl->start();
                    rev_visit.set_backward(!rev_visit.backward());
                    stack.push_back(rev_visit);
#ifdef debug
                    cerr << "Stack up " << stack.back().node_id() << " " << stack.back().backward() << endl;
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
                    Visit rev_visit = into_snarl->start();
                    rev_visit.set_backward(!rev_visit.backward());
                    stack.push_back(rev_visit);
#ifdef debug
                    cerr << "Stack up " << stack.back().node_id() << " " << stack.back().backward() << endl;
#endif
                }
                
                // if the same side is also reachable, add it to the stack too
                if (into_snarl->end_self_reachable()) {
                    Visit rev_visit = into_snarl->end();
                    rev_visit.set_backward(!rev_visit.backward());
                    stack.push_back(rev_visit);
#ifdef debug
                    cerr << "Stack up " << stack.back().node_id() << " " << stack.back().backward() << endl;
#endif
                }
            }
        }
        else {
            // add all of the node traversals we can reach through valid walks to stack
            stack_up_valid_walks(node_handle, stack);
        }
    }
}
    
vector<SnarlTraversal> ExhaustiveTraversalFinder::find_traversals(const Snarl& site) {

    vector<SnarlTraversal> to_return;

    handle_t site_end = graph.get_handle(site.end().node_id(), site.end().backward());
    handle_t site_start = graph.get_handle(site.start().node_id(), site.start().backward());
    handle_t site_rev_start = graph.get_handle(site.start().node_id(), !site.start().backward());
    
    // stop searching when the traversal is leaving the site
    unordered_set<handle_t> stop_at;
    stop_at.insert(site_end);
    stop_at.insert(site_rev_start);
    
    // choose which side(s) can be the end of the traversal
    unordered_set<handle_t> yield_at;
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
        add_traversals(to_return, graph.get_handle(graph.get_id(site_end), !graph.get_is_reverse(site_end)),
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
                                                             int min_recurrence, int max_path_search_steps,
                                                             bool allow_duplicates) :
    graph(graph),
    snarl_manager(snarl_manager),
    reads_by_name(reads_by_name),
    min_recurrence(min_recurrence),    
    max_path_search_steps(max_path_search_steps),
    allow_duplicates(allow_duplicates) {
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
pair<vector<SnarlTraversal>, vector<string>> PathRestrictedTraversalFinder::find_named_traversals(const Snarl& site) {
    // We're going to emit traversals supported by any paths in the graph.
    
    // Put all our subpaths in here to deduplicate them by sequence they spell
    // out. And to count occurrences. Note that the occurrence count will be
    // boosted to min_recurrence if a non-read path in the graph supports a
    // certain traversal string, so we don't end up dropping unsupported ref
    // alleles.
    map<string, tuple<SnarlTraversal, int, string>> results;

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
                    bool enter_start = mapping_enters_side(mapping->to_mapping(),
                        graph.get_handle(site.start().node_id(), site.start().backward()), &graph);

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
                    bool exit_start = mapping_exits_side(mapping->to_mapping(),
                        graph.get_handle(site.start().node_id(), site.start().backward()), &graph);
                    
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
                            !mapping_enters_side(mapping->to_mapping(),
                                graph.get_handle(site.end().node_id(), site.end().backward()), &graph) ||
                            traversal_direction && 
                            !mapping_exits_side(mapping->to_mapping(),
                                graph.get_handle(site.end().node_id(), site.end().backward()), &graph)) {
                            break;
                        }

                        // get the string for the sequence
                        string allele_seq;
                        for (size_t i = 0; i < path_traversed.visit_size(); i++) {
                            auto& path_visit = path_traversed.visit(i);
                            Node* map_node = graph.get_node(path_visit.node_id());
                            allele_seq += path_visit.backward() ? reverse_complement(map_node->sequence()) : map_node->sequence();
                        }

                        // hack for allow_duplicates toggle
                        if (!reads_by_name.count(name) && allow_duplicates) {
                            allele_seq = name;
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
                                get<1>(results[allele_seq])++;
                            } else {
                                // We are a named path (like "ref")
                                if(get<1>(results[allele_seq]) < min_recurrence) {
                                    // Ensure that this allele doesn't get
                                    // eliminated, since ref or some other named
                                    // path supports it.
                                    get<1>(results[allele_seq]) = min_recurrence;
                                } else {
                                    get<1>(results[allele_seq])++;
                                }
                            }
                        } else {
                            // Add it in. Give it a count of 1 if we are a read,
                            // and a count of min_recurrence (so it doesn't get
                            // filtered later) if we are a named non-read path
                            // (like "ref").
                            int trav_occ = reads_by_name.count(name) ? 1 : min_recurrence;
                            results[allele_seq] = std::tie(path_traversed, trav_occ, name);
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
#ifdef debug
#pragma omp critical (cerr)
                        cerr << "traversing left to mapping " << *mapping << endl;
#endif
                    } else {
                        // We're going forwards
                        mapping = graph.paths.traverse_right(mapping);
#ifdef debug
#pragma omp critical (cerr)
                        cerr << "traversing right to mapping " << *mapping << endl;
#endif
                    }
                    // Tick the counter so we don't go really far on long paths.
                    traversal_count++;

                }


            }
        }

    }

    // Now collect the unique results
    pair<vector<SnarlTraversal>, vector<string>> to_return;

    for(auto& result : results) {
        // Break out each result
        const string& seq = result.first;
        auto& traversals = get<0>(result.second);
        auto& count = get<1>(result.second);
        auto& name = get<2>(result.second);

        if(count < min_recurrence) {
            // We don't have enough initial hits for this sequence to justify
            // trying to re-align the rest of the reads. Skip it. Note that the
            // reference path (and other named paths) will stuff in at least
            // min_recurrence to make sure we don't throw out their alleles.
            continue;
        }

        // Send out each list of traversals
        to_return.first.emplace_back(std::move(traversals));
        to_return.second.push_back(name);
    }

    return to_return;
}

vector<SnarlTraversal> PathRestrictedTraversalFinder::find_traversals(const Snarl& site) {
    return find_named_traversals(site).first;
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

PathTraversalFinder::PathTraversalFinder(const PathHandleGraph& graph, SnarlManager& snarl_manager,
                                         const vector<string>& path_names)  :
    graph(graph), snarl_manager(snarl_manager) {
    for (const string& path_name : path_names) {
        assert(graph.has_path(path_name));
        paths.insert(graph.get_path_handle(path_name));
    }
}

vector<SnarlTraversal> PathTraversalFinder::find_traversals(const Snarl& site) {
    return find_path_traversals(site).first;
}

pair<vector<SnarlTraversal>, vector<pair<step_handle_t, step_handle_t> > > PathTraversalFinder::find_path_traversals(const Snarl& site) {

    handle_t start_handle = graph.get_handle(site.start().node_id(), site.start().backward());
    handle_t end_handle = graph.get_handle(site.end().node_id(), site.end().backward());
    
    vector<step_handle_t> start_steps = graph.steps_of_handle(start_handle);
    vector<step_handle_t> end_steps = graph.steps_of_handle(end_handle);

    pair<unordered_set<id_t>, unordered_set<edge_t> > snarl_contents = snarl_manager.deep_contents(&site, graph, true);
    
    // use this to skip paths that don't reach the end node
    unordered_set<path_handle_t> end_path_handles;
    for (const step_handle_t& step : end_steps) {
        end_path_handles.insert(graph.get_path_handle_of_step(step));
    }

#ifdef debug
    cerr << "Finding traversals of " << pb2json(site) << " using PathTraversalFinder" << endl
         << " - there are " << start_steps.size() << " start_steps, " << end_steps.size() << " end_steps"
         << " and " << end_path_handles.size() << " end_path_handles" << endl;
#endif

    vector<SnarlTraversal> out_travs;
    vector<pair<step_handle_t, step_handle_t> > out_steps;

    for (const step_handle_t& start_step : start_steps) {
        path_handle_t start_path_handle = graph.get_path_handle_of_step(start_step);
        // only crawl paths that have a chance of reaching the end
        if ((paths.empty() || paths.count(start_path_handle)) && end_path_handles.count(start_path_handle)) {

            handle_t end_check = end_handle;

#ifdef debug
            cerr << " - considering path " << graph.get_path_name(start_path_handle) << endl;
#endif
            // try to make a traversal by walking forward
            SnarlTraversal trav;
            bool can_continue = true;
            step_handle_t step = start_step;
            while (can_continue) {
                handle_t handle = graph.get_handle_of_step(step);
                Visit* start_visit = trav.add_visit();
                start_visit->set_node_id(graph.get_id(handle));
                start_visit->set_backward(graph.get_is_reverse(handle));

                can_continue = false;
                if (graph.has_next_step(step) && handle != end_handle) {
                    step_handle_t next_step = graph.get_next_step(step);
                    handle_t next_handle = graph.get_handle_of_step(next_step);
                    if (snarl_contents.first.count(graph.get_id(next_handle)) &&
                        snarl_contents.second.count(graph.edge_handle(handle, next_handle))) {
                        step = next_step;
                        can_continue = true;
                    } 
                }
            }

            if (graph.get_handle_of_step(step) != end_check) {
#ifdef debug
                cerr << "     - failed to find forward traversal of path " << graph.get_path_name(start_path_handle) << endl;
#endif
                // try to make a traversal by walking backward
                end_check = graph.flip(end_handle);
                
                trav.Clear();
                can_continue = true;
                step = start_step;
                while (can_continue) {
                    handle_t handle = graph.flip(graph.get_handle_of_step(step));
                    
                    Visit* start_visit = trav.add_visit();
                    start_visit->set_node_id(graph.get_id(handle));
                    start_visit->set_backward(graph.get_is_reverse(handle));

                    can_continue = false;
                    if (graph.has_previous_step(step) && handle != end_handle) {
                        step_handle_t prev_step = graph.get_previous_step(step);
                        handle_t prev_handle = graph.flip(graph.get_handle_of_step(prev_step));

                        if (snarl_contents.first.count(graph.get_id(prev_handle)) &&
                            snarl_contents.second.count(graph.edge_handle(handle, prev_handle))) {
                            step = prev_step;
                            can_continue = true;
                        } 
                    }
                }
            }
            if (graph.get_handle_of_step(step) == end_check) {
                out_travs.push_back(trav);
                out_steps.push_back(make_pair(start_step, step));
            } 
        }
    }
    
    return make_pair(out_travs, out_steps);
}

TrivialTraversalFinder::TrivialTraversalFinder(const HandleGraph& graph) : graph(graph) {
    // Nothing to do!
}

vector<SnarlTraversal> TrivialTraversalFinder::find_traversals(const Snarl& site) {
    assert(site.start_end_reachable());
    assert(site.directed_acyclic_net_graph());
    
    // We'll fill this in and send it back
    vector<SnarlTraversal> to_return;
    
    // We don't want to be duplicating partial paths, so we store for each
    // NodeTraversal we can reach the previous NodeTraversal we can reach it
    // from.
    unordered_map<handle_t, handle_t> previous;
    
    list<handle_t> stack{graph.get_handle(site.start().node_id(), site.start().backward())};
    
    while (!stack.empty()) { 
        // While there's still stuff on the stack
        
        // Grab the first thing
        handle_t here = stack.front();
        stack.pop_front();
        
        if (graph.get_id(here) == site.end().node_id()) {
            // Trace back a path
            list<handle_t> path;
            
            while (true) {
                // Until we get to the start of the site
                
                // Put this traversal on the front of the path
                path.push_front(here);
                
                if (graph.get_id(here) == site.start().node_id()) {
                    // Stop when we've reached the start of the site
                    break;
                }
                
                // Trace back
                here = previous.at(here);
            }
            
            // Initialize a SnarlTraversal in the return value
            to_return.emplace_back();
            
            // Translate the path into the traversal
            for (handle_t node_traversal : path) {
                *(to_return.back().add_visit()) = to_visit(graph, node_traversal);
            }
            
            // Stop early after having found one path
            break;
        } else {
            // We haven't reached the end of the site

            graph.follow_edges(here, false, [&] (const handle_t& next) {
                    // Look at all the places we can go from this node
                    if (!previous.count(next)) {
                        // Remember how we got there
                        previous[next] = here;
                        // Explore it, depth first
                        stack.push_front(next);
                    }
                });
        }
    }
    
    // When we get here, either we found a path, or there isn't one.
    return to_return;
}

RepresentativeTraversalFinder::RepresentativeTraversalFinder(const PathHandleGraph& graph,
                                                             SnarlManager& snarl_manager,
                                                             size_t max_depth,
                                                             size_t max_width,
                                                             size_t max_bubble_paths,
                                                             size_t min_node_support,
                                                             size_t min_edge_support,
                                                             function<PathIndex*(const Snarl&)> get_index,
                                                             function<Support(id_t)> get_node_support,
                                                             function<Support(edge_t)> get_edge_support) :
  graph(graph), snarl_manager(snarl_manager), max_depth(max_depth), max_width(max_width),
  max_bubble_paths(max_bubble_paths), min_node_support(min_node_support), min_edge_support(min_edge_support),
  get_index(get_index), get_node_support(get_node_support), get_edge_support(get_edge_support) {
    has_supports = this->get_node_support != nullptr && this->get_edge_support != nullptr;
    
    // Nothing to do!

}

Path RepresentativeTraversalFinder::find_backbone(const Snarl& site) {
    
    // TODO: this cheats and uses certain things that happen to be true about
    // the TrivialTraversalFinder in order to work.

    // Find a traversal, ignoring the fact that child sites ought to own their
    // nodes.
    TrivialTraversalFinder finder(graph);
    auto traversals = finder.find_traversals(site);
    assert(!traversals.empty());
    auto& traversal = traversals.front();
    
    // Convert it into a path that includes the boundary nodes
    Path to_return;
    for (size_t i = 0; i < traversal.visit_size(); i++) {
        *to_return.add_mapping() = to_mapping(traversal.visit(i), graph);
    }
    
    return to_return;
    
}

vector<SnarlTraversal> RepresentativeTraversalFinder::find_traversals(const Snarl& site) {
    
    // We can only do snarls with start-to-end traversals.
    assert(site.start_end_reachable());
    // And that aren't themselves directed-cyclic
    // TODO: We don't ignore children, but this check does!
    assert(site.directed_acyclic_net_graph());
    
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
    pair<unordered_set<id_t>, unordered_set<edge_t>> contents = snarl_manager.shallow_contents(&site, graph, true);
    
    // Copy its node set
    unordered_set<id_t> nodes_left;
    for (id_t node_id : contents.first) {
        nodes_left.insert(node_id);
    }

    // Trace the ref path through the site.
    vector<Visit> ref_path_for_site;
    
    // First figure where the site starts and ends in the selected path
    size_t site_start = index.by_id.at(site.start().node_id()).first;
    size_t site_end = index.by_id.at(site.end().node_id()).first;
    
#ifdef debug
    cerr << "Site starts with " << pb2json(site.start())
         << " at " << site_start
         << " and ends with " << pb2json(site.end())
         << " at " << site_end << endl;
        
    for (id_t node : nodes_left) {
        cerr << "\tContains node " << node << endl;
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
        id_t visited_node = found_visit.node_id();
        
        const Snarl* child = snarl_manager.into_which_snarl(found_visit);
        if (child != nullptr && child != managed_site &&
            snarl_manager.into_which_snarl(reverse(found_visit)) != managed_site &&
            !(eat_trivial_children && snarl_manager.is_trivial(child, graph))) {
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
            id_t here = visited_node;
            do {
#ifdef debug
                cerr << "at node " << here << endl;
#endif
                // Advance
                ref_node_start = found->first + graph.get_length(graph.get_handle(here));
                // And look at what we get
                found = index.by_start.lower_bound(ref_node_start);
                assert(found != index.by_start.end());
                // And grab out the node
                found_visit = found->second.to_visit();
                here = found_visit.node_id();
                // Until we find something in this parent again that isn't the
                // closing visit of a child snarl. We'll look at what we find
                // next.
            } while (!contents.first.count(here));
            
            if (snarl_manager.into_which_snarl(reverse(found_visit)) != nullptr) {
                // We hit the end node of the child snarl.
                
                if (snarl_manager.into_which_snarl(found_visit) == nullptr) {
                    // We don't have another child snarl immediately. Look at the node after this one.
                    ref_node_start = found->first + graph.get_length(graph.get_handle(here));
                    found = index.by_start.lower_bound(ref_node_start);
                    assert(found != index.by_start.end());
                    found_visit = found->second.to_visit();
                    here = found_visit.node_id();
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
            ref_node_start = found->first + graph.get_length(graph.get_handle(visited_node));
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
        
        if (snarl_manager.into_which_snarl(node, true) || snarl_manager.into_which_snarl(node, false)) {
            // Skip child boundary nodes.
            continue;
        }
        
        if(index.by_id.count(node)) {
            cerr << "error[RepresentativeTraversalFinder]: Node " << node << " is on backbone path at "
                 << index.by_id.at(node).first << " but not traced in site "
                 << pb2json(site) << endl;
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
            if(visit.node_id() != 0) {
                auto found = index.find_in_orientation(visit.node_id(), visit.backward());
                if (found != index.end()) {
                    cerr << "\tPath member " << visit << " lives on backbone at "
                         << found->first << endl;
                } else {
                    cerr << "\tPath member " << visit << " does not live on backbone" << endl;
                }
            } else {
                cerr << "\tPath member " << visit << " is to a child snarl" << endl;
                
                auto found_start = index.find_in_orientation(visit.snarl().start().node_id(),
                    visit.snarl().start().backward() != visit.backward());
                    
                if (found_start != index.end()) {
                    cerr << "\t\tStart lives on backbone at "
                         << found_start->first << endl;
                } else {
                    cerr << "\t\tStart does not live on backbone" << endl;
                }
                
                auto found_end = index.find_in_orientation(visit.snarl().end().node_id(),
                    visit.snarl().end().backward() != visit.backward());
                
                if (found_end != index.end()) {
                    cerr << "\t\tEnd lives on backbone at "
                         << found_end->first << endl;
                } else {
                    cerr << "\t\tEnd does not live on backbone" << endl;
                }
                    
                
            }
        }
#endif
        
        for(auto& visit : path) {
            if (visit.node_id() != 0) {
                // Make sure the site actually has the nodes we're visiting.
                if (!contents.first.count(visit.node_id())) {
                    cerr << "error[RepresentativeTraversalFinder::find_traversals]: Node "
                        << visit.node_id() << " not in snarl " << pb2json(site) << " contents:" << endl;
                    
                    for (auto& node_id : contents.first) {
                        cerr << "\t" << node_id << "," << endl;
                    }
                    
                    cerr << "children:" << endl;
                    
                    for (auto& snarl_ptr : snarl_manager.children_of(&site)) {
                        cerr << pb2json(*snarl_ptr) << endl;
                    }
                    
                    cerr << "Input path: " << endl;
                    for(auto& visit : path) {
                        if(visit.node_id() != 0) {
                            auto found = index.find_in_orientation(visit.node_id(), visit.backward());
                            if (found != index.end()) {
                                cerr << "\tPath member " << visit << " lives on backbone at "
                                     << found->first << endl;
                            } else {
                                cerr << "\tPath member " << visit << " does not live on backbone" << endl;
                            }
                        } else {
                            cerr << "\tPath member " << visit << " is to a child snarl" << endl;
                            
                            auto found_start = index.find_in_orientation(visit.snarl().start().node_id(),
                                visit.snarl().start().backward() != visit.backward());
                                
                            if (found_start != index.end()) {
                                cerr << "\t\tStart lives on backbone at "
                                     << found_start->first << endl;
                            } else {
                                cerr << "\t\tStart does not live on backbone" << endl;
                            }
                            
                            auto found_end = index.find_in_orientation(visit.snarl().end().node_id(),
                                visit.snarl().end().backward() != visit.backward());
                            
                            if (found_end != index.end()) {
                                cerr << "\t\tEnd lives on backbone at "
                                     << found_end->first << endl;
                            } else {
                                cerr << "\t\tEnd does not live on backbone" << endl;
                            }
                                
                            
                        }
                    }
                
                    assert(false);
                }
            }
            // Child snarl end nodes will still appear in our contents.
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
        
#ifdef debug
        cerr << "Ref path length: " << ref_path_for_site.size() << " visits" << endl;
        cerr << "Path to be anchored: " << path.size() << " visits" << endl;
        cerr << "Looking for " << frontier_visit(path.at(bubble_path_index), true)
            << " or " << frontier_visit(path.at(bubble_path_index), false) << " exiting an anchoring snarl" << endl;

        cerr << "Check pos " << ref_path_index << " on ref path and " << bubble_path_index << " on path to be anchored" << endl;
#endif
        
        
        
        while(frontier_visit(ref_path_for_site.at(ref_path_index), false) != frontier_visit(path.at(bubble_path_index), true) &&
              !(path.at(bubble_path_index).node_id() == 0 &&
                frontier_visit(ref_path_for_site.at(ref_path_index), false) == frontier_visit(path.at(bubble_path_index), false))) {
            // The right visit of where we are on the ref path isn't the left
            // visit of where we want to start, nor is it the end of a snarl
            // and the right visit of where we want to start.
            
            // Collect NodeTraversals from the ref path until we hit the one
            // at which the bubble path starts.
#ifdef debug
            cerr << "Before path: " << pb2json(ref_path_for_site.at(ref_path_index)) << endl;
#endif
            extended_path.push_back(ref_path_for_site.at(ref_path_index++));
            
#ifdef debug
            cerr << "Check pos " << ref_path_index << " on ref path and " << bubble_path_index << " on path to be anchored" << endl;
#endif
            if (ref_path_index >= ref_path_for_site.size()) {
                // We hit the end of the reference path. If the path we are
                // trying to anchor actually starts and ends along the
                // reference in the right orientation, this should never
                // happen.
                throw runtime_error("Ran out of reference path when looking for start of path to be anchored");
            }
        }
        
        if (ref_path_for_site.at(ref_path_index).node_id() == 0) {
            // The last Visit we traversed from the ref was a Snarl, so it already
            // includes the first node of the path as one of its boundaries. We need
            // to add the ref visit and exclude the bubble visit unless it is also of
            // a child Snarl and that Snarl is different from the ref path
            
#ifdef debug
            cerr << "Adding final ref child visit " << pb2json(ref_path_for_site.at(ref_path_index)) << endl;
#endif
            extended_path.push_back(ref_path_for_site.at(ref_path_index));
            
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
            cerr << "Skip ref: " << pb2json(ref_path_for_site.at(ref_path_index)) << endl;
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
        
        if (ref_path_for_site.at(ref_path_index).node_id() == 0) {
            // The next Visit on the ref path is to a Snarl, so the final Visit we added
            // from the bubble will be redundant with its boundary nodes unless that Visit
            // was also to a Snarl
            if (extended_path.back().node_id() != 0 ||
                (ref_path_for_site.at(ref_path_index).snarl().start() == extended_path.back().snarl().start()
                 && ref_path_for_site.at(ref_path_index).snarl().end() == extended_path.back().snarl().end())) {
#ifdef debug
                cerr << "Removing bubble visit " << pb2json(extended_path.back()) << endl;
#endif
                extended_path.pop_back();
            }
#ifdef debug
            cerr << "Adding adjacent ref child visit" << pb2json(ref_path_for_site.at(ref_path_index)) << endl;
#endif
            extended_path.push_back(ref_path_for_site.at(ref_path_index));
        }
        // Skip the matching NodeTraversal
        ref_path_index++;
        while(ref_path_index < ref_path_for_site.size()) {
            // Then take the entier rest of the ref path
#ifdef debug
            cerr << "After path: " << pb2json(ref_path_for_site.at(ref_path_index)) << endl;
#endif
            extended_path.push_back(ref_path_for_site.at(ref_path_index++));
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

    for (id_t node_id : contents.first) {
        // Find the bubble for each node
        
        if (snarl_manager.into_which_snarl(node_id, true) || snarl_manager.into_which_snarl(node_id, false)) {
            // Don't start from nodes that are child boundaries
            continue;
        }
        
        if (has_supports && total(get_node_support(node_id)) < min_node_support) {
            // Don't bother with unsupported nodes
            continue;
        }
        
        if (index.by_id.count(node_id)) {
            // Don't try to pathfind to the backbone for backbone nodes.
            continue;
        }
        
#ifdef debug
        cerr << "Base path on " << node_id << endl;
#endif
        
        // Find bubbles that backend into the backbone path
        pair<Support, vector<Visit>> sup_path = find_bubble(node_id, nullptr, nullptr, index, site);

        vector<Visit>& path = sup_path.second;
        
        if(path.empty()) {
            // We couldn't find a path back to the primary path. Discard
            // this material.
            if (verbose) {
                cerr << "Warning: No path found for node " << node_id << endl;
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

    for(const edge_t edge : contents.second) {
        // Go through all the edges
        
        if(has_supports && total(get_edge_support(edge)) < min_edge_support) {
            // Don't bother with unsupported edges
#ifdef debug
            cerr << "Skip unsupported edge " << graph.get_id(edge.first) << ":" << graph.get_is_reverse(edge.first)
                 << " -> " << graph.get_id(edge.second) << ":" << graph.get_is_reverse(edge.second) << endl;
#endif
            continue;
        }
        
        if(!index.by_id.count(graph.get_id(edge.first)) || !index.by_id.count(graph.get_id(edge.second))) {
            // Edge doesn't touch backbone at both ends. Don't use it
            // because for some reason it makes performance worse
            // overall.
#ifdef debug
            cerr << "Skip off-backbone edge " << graph.get_id(edge.first) << ":" << graph.get_is_reverse(edge.first)
                 << " -> " << graph.get_id(edge.second) << ":" << graph.get_is_reverse(edge.second) << endl;
#endif
            continue;
        }
        
#ifdef debug
        cerr << "Base path on " << graph.get_id(edge.first) << ":" << graph.get_is_reverse(edge.first)
             << " -> " << graph.get_id(edge.second) << ":" << graph.get_is_reverse(edge.second) << endl;
#endif
        
        // Find a path based around this edge
        pair<Support, vector<Visit>> sup_path = find_bubble(0, &edge, nullptr, index, site);
        vector<Visit>& path = sup_path.second;
        
#ifdef debug
        cerr << "Edge " << graph.get_id(edge.first) << ":" << graph.get_is_reverse(edge.first)
             << " -> " << graph.get_id(edge.second) << ":" << graph.get_is_reverse(edge.second)
            << " yields:" << endl;
        for(auto& visit : path) {
            cerr << "\t" << visit << endl;
        }
#endif
        
        if(path.empty()) {
            // We couldn't find a path back to the primary path. Discard
            // this material.
            if (verbose) {
                cerr << "Warning: No path found for edge " << graph.get_id(edge.first) << ":" << graph.get_is_reverse(edge.first)
                     << " -> " << graph.get_id(edge.second) << ":" << graph.get_is_reverse(edge.second) << endl;
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

        if (eat_trivial_children && snarl_manager.is_trivial(child, graph)) {
            // Skip trivial children
            continue;
        }
        
#ifdef debug
        cerr << "Base path on " << *child << endl;
#endif
        
        // Find a path based around this child snarl
        pair<Support, vector<Visit>> sup_path = find_bubble(0, nullptr, child, index, site);
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

pair<Support, vector<Visit>> RepresentativeTraversalFinder::find_bubble(id_t node, const edge_t* edge,
                                                                        const Snarl* snarl, PathIndex& index, const Snarl& site) {

    // What are we going to find our left and right path halves based on?
    Visit left_visit;
    Visit right_visit;
    
    const Snarl* managed_site = snarl_manager.manage(site);
    
    if (edge != nullptr) {
        // Be edge-based
        
        // Find the nodes at the ends of the edges. Look at them traversed in the
        // edge's local orientation.
        left_visit = to_visit(graph, edge->first);
        right_visit = to_visit(graph, edge->second);
        
        // Find any child snarls looking out form the edge
        const Snarl* right_child = snarl_manager.into_which_snarl(right_visit);
        const Snarl* left_child = snarl_manager.into_which_snarl(reverse(left_visit));
        
        if (right_child != nullptr && right_child != managed_site
            && snarl_manager.into_which_snarl(reverse(right_visit)) != managed_site) {
            // We're reading into a child snarl on the right.
            // And we're not reading out of ourselves.
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
            && snarl_manager.into_which_snarl(left_visit) != managed_site) {
            // We're reading out of a child snarl on the left.
            // And we're not reading into ourselves.
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
        
    } else if (node != 0) {
        // Be node-based. TODO: we trust the caller not to feed us nodes that
        // are part of/boundaries of child snarls.
        left_visit = right_visit = to_visit(node, false);
    } else {
        // Be snarl-based
        assert(snarl != nullptr);
        left_visit = right_visit = to_visit(*snarl);
    }
    
#ifdef debug
    cerr << "Starting from: " << left_visit << ", " << right_visit << endl;
#endif

    // Find paths on both sides, with nodes or snarls on the primary path at
    // the outsides and this visit in the middle. Returns path lengths,
    // orientatioins in which the reference path was encountered, and paths in
    // tuples in a set.
    // Make sure to keep looking for the other orientation of the ref path
    // after we find a first one, to handle inversions of up to a certain size.
    auto leftPaths = bfs_left(left_visit, index, false, managed_site, other_orientation_timeout);
    auto rightPaths = bfs_right(right_visit, index, false, managed_site, other_orientation_timeout);
    
    // Sort out the paths not just by whether they are left or right from here,
    // but also by whether they hit the reference path in forward or reverse
    // ref-path-relative orientation.
    // TODO: give ImmutableList a .back() so we can avoid converting to real lists here.
    list<list<Visit>> left_forward;
    list<list<Visit>> right_forward;
    list<list<Visit>> left_reverse;
    list<list<Visit>> right_reverse;
    
    for (auto& annotatedPath : leftPaths) {
        // Break up the paths on the left by orientation
        auto& ref_reverse = get<1>(annotatedPath);
        auto& path = get<2>(annotatedPath);
        // TODO: ImmutableList iterators don't actually satisfy
        // https://en.cppreference.com/w/cpp/named_req/Iterator because they
        // lack the tyypedefs for std::iterator_traits. So we can't use them to
        // construct lists. So we have to build the lists manually.
        list<Visit> converted;
        for (auto& item : path) {
            converted.push_back(item);
        }
        (ref_reverse ? left_reverse : left_forward).emplace_back(move(converted));
    }
    
    for (auto& annotatedPath : rightPaths) {
        // Break up the paths on the right by orientation
        auto& ref_reverse = get<1>(annotatedPath);
        auto& path = get<2>(annotatedPath);
        // TODO: ImmutableList iterators don't actually satisfy
        // https://en.cppreference.com/w/cpp/named_req/Iterator because they
        // lack the tyypedefs for std::iterator_traits. So we can't use them to
        // construct lists. So we have to build the lists manually.
        list<Visit> converted;
        for (auto& item : path) {
            converted.push_back(item);
        }
        (ref_reverse ? right_reverse : right_forward).emplace_back(move(converted));
    }
    
    // We need to look in different combinations of lists.
    auto testCombinations = [&](const list<list<Visit>>& leftList,
                                const list<list<Visit>>& rightList) -> pair<Support, vector<Visit>> {
                                
                                
        // Find a combination of two paths which gets us to the reference and
        // which doesn't use the same nodes on both sides. Track support of up
        // to max_bubble_paths combinations, and return the highest. Always
        // returns the combined path in a valid reference-relative-forward
        // orientation.
        
        // Because we do our own identification of the anchoring reverence
        // occurrences, we may produce a reference-relative-forward path from
        // what was supposed to be a reference-relative-backward pair of
        // partial paths.
        
        // TODO: Fix that by making the BFS code pass along the particular
        // anchoring occurrences it finds
        
        pair<Support, vector<Visit> > bestBubblePath;
        int bubbleCount = 0;
        
#ifdef debug        
        cerr << "Combine " << leftList.size() << " left sides and "
        << rightList.size() << " right sides" << endl;
#endif
        
        // We know the left list starts and the right list ends with an actual
        // node visit, if only to the snarl's start or end.

        for(auto leftPath : leftList) {
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
            
            // TODO: We're using the first occurrence, because that's what
            // we'll encounter as we scan along the reference path and what
            // we'll use to try and build the final traversal. This may NOT be
            // the occurence that got this partial path into the collection for
            // this particular relative orientation. So we still need to check
            // on orientation consistency.
            
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
                        
#ifdef debug
                        cerr << "Anchored to ref path backward! Reverse combination!" << endl;
#endif
                        
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
    
    // Find the best valid combination, if any, in each orientation
#ifdef debug
    cerr << "Combine forward paths" << endl;
#endif
    pair<Support, vector<Visit> > best_forward = testCombinations(left_forward, right_forward);
#ifdef debug
    cerr << "Combine reverse paths" << endl;
#endif
    pair<Support, vector<Visit> > best_reverse = testCombinations(left_reverse, right_reverse);
    
#ifdef debug
    cerr << "Best forward path:" << endl;
    for (auto& visit : best_forward.second) {
        cerr << "\t" << visit << endl;
    }
    cerr << "Best reverse path (in forward orientation):" << endl;
    for (auto& visit : best_reverse.second) {
        cerr << "\t" << visit << endl;
    }
#endif
    
    if (total(best_forward.first) > total(best_reverse.first) || best_reverse.second.empty()) {
        // The forward orientation wins
        return best_forward;
    } else {
        // The reverse orientation wins.
        // testCombinations already made it be reference-forward.
        return best_reverse;
    }
}

Support RepresentativeTraversalFinder::min_support_in_path(const list<Visit>& path) {
    
    if (path.empty() || !has_supports) {
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
        minSupport = get_node_support(cur->node_id());
        if (cur->backward()) {
            // Put the support in the path forward direction
            minSupport = flip(minSupport);
        }
        supportFound = true;
    }
    
    for (; next != path.end(); ++cur, ++next) {
        // For each visit and its next visit
    
        if (next->node_id() != 0) {
            // The next visit is to a node, so get its support
            Support nextSupport = get_node_support(next->node_id());
            
            if (next->backward()) {
                // This occurs backward on the path, so flip its support
                nextSupport = flip(nextSupport);
            }
            
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
        NodeSide from_side = to_right_side(*cur);
        NodeSide to_side = to_left_side(*next);
        edge_t edge = graph.edge_handle(graph.get_handle(from_side.node, !from_side.is_end),
                                        graph.get_handle(to_side.node, to_side.is_end));
                          
        if (graph.has_edge(edge.first, edge.second)) {
            // The edge exists (because we aren't back-to-back child snarls)
            Support edgeSupport = get_edge_support(edge);
            
            if (cur->node_id() > next->node_id() || (cur->node_id() == next->node_id() && cur->backward())) {
                // We are taking the edge backward, so flip its support
                edgeSupport = flip(edgeSupport);
            }
            
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

set<tuple<size_t, bool, structures::ImmutableList<Visit>>>
RepresentativeTraversalFinder::bfs_left(Visit visit,
                                        PathIndex& index, bool stop_if_visited, const Snarl* in_snarl,
                                        size_t both_orientations_distance) {

    // Holds partial paths we want to return, with their lengths in bp and
    // target-path-relative orientations.
    set<tuple<size_t, bool, structures::ImmutableList<Visit>>> toReturn;
    
    // Do a BFS
    
    // Define a stack frame fro the BFS to track an outstanding path being explored.
    // Stores the path, the path length in nodes, the countdown to reach the
    // target path in the other orientation (or 0 if we have not yet reached
    // the target path), and a flag for if the target path was reached in
    // reverse orientation when we started the countdown.  
    using frame_t = tuple<structures::ImmutableList<Visit>, size_t, size_t, bool>;
    
    // This holds the paths to get to NodeTraversals to visit (all of which will
    // end with the node we're starting with).
    list<frame_t> toExtend;
    
    // This keeps a set of all the oriented nodes we already got to and don't
    // need to queue again.
    set<Visit> alreadyQueued;
    
    // Start at this node, with no visits and no countdown running.
    toExtend.emplace_back(visit, 0, 0, false);
    // Mark this traversal as already queued
    alreadyQueued.insert(visit);
    
    // How many ticks have we spent searching?
    size_t searchTicks = 0;
    
#ifdef debug
    cerr << "Start BFS left from " << visit << endl;
    
    if (in_snarl != nullptr) {
        cerr << "Stay inside " << pb2json(*in_snarl) << endl;
    }
    
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

        
        // Dequeue a frame to extend.
        // Make sure to move out of the list to avoid a useless copy.
        frame_t frame(move(toExtend.front()));
        toExtend.pop_front();
        stillToExtend--;
        
        // Unpack it
        auto& path = get<0>(frame);
        auto& path_length = get<1>(frame);
        auto& countdown = get<2>(frame);
        auto& first_found_reverse = get<3>(frame);
        
        // Determine an effective node ID and orientation. Either the one we
        // actually have, or the one on the other side of the snarl we're
        // visiting if we have no node.
        size_t node_id;
        bool is_reverse;
        if (path.front().node_id() != 0) {
            // We are visiting a node so just unpack it
            node_id = path.front().node_id();
            is_reverse = path.front().backward();
        } else {
            // We are visiting a snarl
            
            if (path.front().backward()) {
                // We are using the snarl in reverse. Since we are searching
                // left, we want what's on the left of the reversed snarl,
                // which is its end node.
                node_id = path.front().snarl().end().node_id();
                // Since we're using the snarl in reverse, we invert that end node's orientation
                is_reverse = !path.front().snarl().end().backward();
            } else {
                // We're visiting a snarl forward, so we use its start node in the orientation the snarl does.
                node_id = path.front().snarl().start().node_id();
                is_reverse = path.front().snarl().start().backward();
            }
        }
       
        // Determine if we connect to the forward orientation and/or reverse orientation of the target path
        pair<bool, bool> orientations = index.get_contained_orientations(node_id);
        if (is_reverse) {
            // We actually hit the path in the opposite orientation.
            // So if the node we hit is on the path in reverse, we hit the path forward, because we used the node in reverse.
            std::swap(orientations.first, orientations.second);
        }
       
        // Now that we know where we are, work out where we want to be
       
        // Determine if we want forward orientation hits
        bool want_forward = countdown == 0 || first_found_reverse;
        // And if we want reverse orientation hits
        bool want_reverse = countdown == 0 || !first_found_reverse;
            
        // This flag will determine if we want to extend from here
        bool extend = false;
            
        // We want to lazily compute the path length in bases
        size_t length = 0;
            
        if (want_forward) {
            if (orientations.first) {
                // Process a wanted hit in forward orientation
                
                length = bp_length(path);
                
#ifdef debug
                cerr << "Reached anchoring node " << node_id << " on target path forward" << endl;
                cerr << "Emit path of length " << path_length << " visits and " << length << " bp to forward" << endl;
#endif

                toReturn.emplace(length, false, path);
                
                if (want_reverse) {
                    // We found forward first and now need to start the countdown
                    first_found_reverse = false;
                    countdown = both_orientations_distance;
                }
                
            } else {
                // We still want forward.
                if (want_reverse) {
                    // We also want reverse; No countdown is active.
                    extend = true;
                } else {
                    // We are doing a countdown because we already found reverse
                    if (countdown > 1) {
                        // There is still time, so extend.
                        countdown--;
                        extend = true;
                    }
                }
            }
        }
        
        if (want_reverse) {
            if(orientations.second) {
                // Process a wanted hit in reverse orientation
                
                if (length == 0) {
                    length = bp_length(path);
                }
            
#ifdef debug
                cerr << "Reached anchoring node " << node_id << " on target path forward" << endl;
                cerr << "Emit path of length " << path_length << " visits and " << length << " bp to reverse" << endl;
#endif

                toReturn.emplace(length, true, path);
                
                if (want_forward) {
                    // We found reverse first and now need to start the countdown
                    first_found_reverse = true;
                    countdown = both_orientations_distance;
                }
            } else {
                // We still want reverse.
                if (want_forward) {
                    // We also want forward; No countdown is active.
                    extend = true;
                } else {
                    // We are doing a countdown because we already found forward
                    if (countdown > 1) {
                        // There is still time, so extend.
                        countdown--;
                        extend = true;
                    }
                }
            }
        }
        
        
        if (path_length >= max_depth) {
#ifdef debug
            cerr << "Path has reached max depth! Aborting!" << endl;
#endif
        } else if (!extend) {
            // We chose not to extend.
#ifdef debug
            cerr << "Choosing not to extend" << endl;
#endif
        } else if (in_snarl != nullptr &&
            ((node_id == in_snarl->start().node_id() && is_reverse == in_snarl->start().backward()) ||
            (node_id == in_snarl->end().node_id() && is_reverse != in_snarl->end().backward()))) {
            // We hit a boundary node of the snarl we are working on, and are
            // headed out of the snarl (i.e. we're at the start or end in the
            // into-snarl orientation).
#ifdef debug
            cerr << "Path has reached containing snarl boundary! Aborting!" << endl;
#endif
        } else {
            // We haven't hit the reference path yet in all orientations, but
            // we also haven't hit the max depth or the snarl bounds. Extend
            // with all the possible extensions.
            
            // Look left, possibly entering child snarls
            vector<Visit> prevVisits = snarl_manager.visits_left(path.front(), graph, in_snarl);
            
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

                    NodeSide from_side = to_right_side(prevVisit);
                    NodeSide to_side = to_left_side(path.front());
                    edge_t edge = graph.edge_handle(graph.get_handle(from_side.node, !from_side.is_end),
                                        graph.get_handle(to_side.node, to_side.is_end));

                    // Make sure the edge is real, since it can't be a back-to-
                    // back site
                    assert(graph.has_edge(edge.first, edge.second));
                
                    // Fetch the actual node
                    id_t prevNode = prevVisit.node_id();
                    
                    if (has_supports &&
                        (total(get_node_support(prevNode)) < min_node_support ||
                         total(get_edge_support(edge)) < min_edge_support)) {
                        // We have no support at all for visiting this node by this
                        // edge (but we do have some read support data)
                        
#ifdef debug
                        cerr << "Reject " << prevNode << " with no support" << endl;
#endif
                        
                        continue;
                    }
                } else {
                    // This is a visit to a child snarl
                    
                    // Look at the node we would leave the child snarl on
                    // That node can't be shared with a snarl we are already at.
                    id_t prevNode = to_left_side(prevVisit).node;
                    
                    if (has_supports && total(get_node_support(prevNode)) < min_node_support) {
                        // We have no support at all for visiting the far node of this snarl
                        
#ifdef debug
                        cerr << "Reject " << prevVisit.snarl() << " with no support on far node" << endl;
#endif
                        
                        continue;
                    }
                    
                    // TODO: when snarls are not back-to-back, check the connecting edges
                }
                
                
                if (stop_if_visited && alreadyQueued.count(prevVisit)) {
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
                toExtend.emplace_back(path.push_front(prevVisit), path_length + 1, countdown, first_found_reverse);
                stillToExtend++;
                
                // Remember we found a way to this node, so we don't try and
                // visit it other ways.
                alreadyQueued.insert(prevVisit);
            }
        } 
    }
    
    return toReturn;
}

set<tuple<size_t, bool, structures::ImmutableList<Visit>>>
RepresentativeTraversalFinder::bfs_right(Visit visit, PathIndex& index, bool stop_if_visited,
                                         const Snarl* in_snarl, size_t both_orientations_distance) {

    // Look left from the backward version of the visit.
    auto toConvert = bfs_left(reverse(visit), index, stop_if_visited, in_snarl, both_orientations_distance);
    
    // Since we can't modify set records in place, we need to do a copy
    set<tuple<size_t, bool, structures::ImmutableList<Visit>>> toReturn;
    
    for(auto lengthOrientationAndPath : toConvert) {
        // Unpack
        auto& length = get<0>(lengthOrientationAndPath);
        auto& orientation = get<1>(lengthOrientationAndPath);
        auto& path = get<2>(lengthOrientationAndPath);
        
        // Flip every path to run the other way
        // TODO: this duplicates previously shared nodes...
        structures::ImmutableList<Visit> reverse_path;
        for (auto& item : path) {
            // While we're at it, reverse each visit
            reverse_path = reverse_path.push_front(reverse(item));
        }
        
        // Stick it in the new set.
        // Also flip the orientation flag, since if we encountered the forward
        // version of a path searchilg left, we should hit the reverse version
        // searching right, and visa versa.
        toReturn.emplace(length, !orientation, reverse_path);
    }
    
    return toReturn;
}

size_t RepresentativeTraversalFinder::bp_length(const structures::ImmutableList<Visit>& path) {
    size_t length = 0;
    for(auto& visit : path) {
        // Sum up length of each node's sequence
        if (visit.node_id() != 0) {
            length += graph.get_length(graph.get_handle(visit.node_id()));
        }
        // TODO: handle nested sites
    }
    return length;
}


VCFTraversalFinder::VCFTraversalFinder(const PathHandleGraph& graph, SnarlManager& snarl_manager,
                                       vcflib::VariantCallFile& vcf,
                                       const vector<string>& ref_path_names,
                                       FastaReference* ref_fasta,
                                       FastaReference* ins_fasta,
                                       function<bool(const SnarlTraversal&, int)> skip_alt,
                                       size_t max_traversal_cutoff) :
    graph(graph),
    snarl_manager(snarl_manager),
    skip_alt(skip_alt),
    max_traversal_cutoff(max_traversal_cutoff),
    path_finder(graph, snarl_manager, ref_path_names) {

    create_variant_index(vcf, ref_fasta, ins_fasta);
}

VCFTraversalFinder::~VCFTraversalFinder() {
    delete_variant_index();
}

void VCFTraversalFinder::create_variant_index(vcflib::VariantCallFile& vcf, FastaReference* ref_fasta,
                                              FastaReference* ins_fasta) {

    vcflib::Variant var;
#ifdef debug
    cerr << "indexing vcf using alt-path information from graph" << endl;
#endif
    vector<FastaReference*> insertion_fastas;
    if (ins_fasta != nullptr) {
        insertion_fastas.push_back(ins_fasta);
    }
    
    while (vcf.getNextVariant(var)) {
        bool path_found = false;
        path_handle_t path_handle;

        // we need to run this in order for symbolic alleles to get the same hashes as in construct
        if (var.isSymbolicSV()) {
            if (ref_fasta == nullptr) {
                cerr << "[VCFTraversalFinder] Warning: Unable to canonicalize symbolic variant because no reference fasta"
                     << " was given:\n" << var << endl;
                continue;
            }
            bool could_canonicalize = var.canonicalize(*ref_fasta, insertion_fastas, true);
            if (!could_canonicalize) {
                cerr << "[VCFTraversalFinder] Warning: Failed to canonicalize symbolic variant:\n" << var << endl;
                continue;
            }
        }

        // scan paths in the graph for any alt path that could have come from this variant
        // then add any node id from the path to our index
        // we add the first id we find under the assumption that alt paths are entirely contained within sites
        for (int allele = 0; !path_found && allele < var.alleles.size(); ++allele) {
            string alt_path_name = "_alt_" + make_variant_id(var) + "_" + to_string(allele);
            if (graph.has_path(alt_path_name)) {
                path_handle_t path_handle = graph.get_path_handle(alt_path_name);
                if (!graph.is_empty(path_handle)) {
                    path_found = true;
                    step_handle_t step_handle = graph.path_begin(path_handle);
                    handle_t handle = graph.get_handle_of_step(step_handle);
                    id_t node_id = graph.get_id(handle);
                    // copy our variant just this once, and add its new pointer to our map
                    if (node_to_variant.count(node_id)) {
                        node_to_variant[node_id].push_back(new vcflib::Variant(var));
                    } else {
                        node_to_variant[node_id] = list<vcflib::Variant*>({new vcflib::Variant(var)});
                    }
                }
            }
        }
        if (!path_found) {
            cerr << "[VCFTraversalFinder] Warning: No alt path (prefix="
                 << ("_alt_" + make_variant_id(var) + "_") << ") found in graph for variant.  It will be ignored:\n"
                 << var << endl;
        }
    }
#ifdef debug
    cerr << "Indexed " << node_to_variant.size() << " nodes" << endl;
#endif
}

void VCFTraversalFinder::delete_variant_index() {
    for (auto nv : node_to_variant) {
        for (auto var : nv.second) {
            delete var;
        }
    }
    node_to_variant.clear();
}

vector<vcflib::Variant*> VCFTraversalFinder::get_variants_in_site(const Snarl& site) {
    vector<vcflib::Variant*> site_variants;

    pair<unordered_set<id_t>, unordered_set<edge_t> > contents = snarl_manager.deep_contents(&site, graph, false);

    for (auto node_id : contents.first) {
        auto map_it = node_to_variant.find(node_id);
        if (map_it != node_to_variant.end()) {
            for (auto var : map_it->second) {
                site_variants.push_back(var);
            }
        }
    }
            
    return site_variants;
}

pair<vector<pair<SnarlTraversal, vector<int>>>, vector<vcflib::Variant*>>
VCFTraversalFinder::find_allele_traversals(Snarl site) {

    vector<pair<SnarlTraversal, vector<int>>> output_traversals;

    // This traversal finder is pretty simple-minded.  It's expecting forward-oriented variation relative
    // to the reference.  We flip our snarl to canonicalize it if possible.
    if (site.start().backward() && site.end().backward()) {
        Visit start = site.start();
        *site.mutable_start() = site.end();
        *site.mutable_end() = start;
        site.mutable_start()->set_backward(false);
        site.mutable_end()->set_backward(false);
    }
    
    vector<vcflib::Variant*> site_variants = get_variants_in_site(site);

    if (site_variants.empty()) {
        return make_pair(output_traversals, site_variants);
    }

    pair<vector<SnarlTraversal>, vector<pair<step_handle_t, step_handle_t> > > trav_steps = path_finder.find_path_traversals(site);

    if (trav_steps.first.empty()) {
        return make_pair(output_traversals, site_variants);
    }

    // we can certainly relax this if needed
    assert(trav_steps.first.size() == 1);

    step_handle_t start_step = trav_steps.second[0].first;
    step_handle_t end_step = trav_steps.second[0].second;
    path_handle_t ref_path = graph.get_path_handle_of_step(start_step);

#ifdef debug
    cerr << "Computing alt traversals for site " << pb2json(site) << " that spans the following "
         << site_variants.size() << " variants:\n";
    for (auto site_var : site_variants) {
      cerr << " ->" << *site_var << endl;
    }
#endif

    // fill in the alt traversals
    brute_force_alt_traversals(site, site_variants, ref_path, start_step, end_step, output_traversals);

    return make_pair(output_traversals, site_variants);
}

vector<SnarlTraversal> VCFTraversalFinder::find_traversals(const Snarl& site) {
    pair<vector<pair<SnarlTraversal, vector<int>>>, vector<vcflib::Variant*>> allele_travs = find_allele_traversals(site);
    vector<SnarlTraversal> traversals;
    traversals.reserve(allele_travs.first.size());
    for (auto& trav : allele_travs.first) {
        traversals.push_back(trav.first);
    }
    return traversals;
}

void VCFTraversalFinder::brute_force_alt_traversals(
    const Snarl& site,
    const vector<vcflib::Variant*>& site_variants,
    path_handle_t ref_path,
    step_handle_t start_step,
    step_handle_t end_step,
    vector<pair<SnarlTraversal, vector<int> > >& output_traversals) {

    // the haplotype we're going to look for a traversal for
    // it's in terms of alt_alleles below (and not the VCF), so needs
    // to be converted back
    vector<int> haplotype(site_variants.size(), 0);

    // use our skip_alt() method (if defined) to prune the search space
    vector<vector<int>> alt_alleles = get_pruned_alt_alleles(site, site_variants, ref_path);
    assert(alt_alleles.size() == haplotype.size());
    
    // if we failed to prune enough, we print a message here:
    // todo: we can move to a ranking (eg by support), where instead of filtering, we just
    // take the K most supported traversals.  this would avoid ever skipping a site
    if (!check_max_trav_cutoff(alt_alleles)) {
        cerr << "[VCFTraversalFinder] Warning: Site " << pb2json(site) << " with " << site_variants.size()
             << " variants contains too many traversals (>" << max_traversal_cutoff 
             << ") to enumerate so it will be skipped:\n";
        for (auto site_var : site_variants) {
            cerr << "   " << *site_var << endl;
        }
        output_traversals.clear();
        return;
    }

    // increment the haplotype.  we can use this to loop over every possible haplotype
    auto next_haplotype = [&] () -> bool {
        // do this by "adding" 1 to our haplotype. each digit is in base-|alleles|
        for (int i = alt_alleles.size() - 1; i >= 0; --i) {
            if (haplotype[i] < alt_alleles[i].size() - 1) {
                // add to column
                ++haplotype[i];
                return true;
            } else if (i > 0) {
                // carry 1 to left
                haplotype[i] = 0;
            }
        }
        return false;
    };

    do {
        // convert back to vcf allele offsets
        // todo: can we change the enumeration to avoid this?
        vector<int> vcf_haplotype(haplotype.size());
        for (int i = 0; i < site_variants.size(); ++i) {
            vcf_haplotype[i] = (alt_alleles[i][haplotype[i]]);
            assert(skip_alt != nullptr || vcf_haplotype[i] == haplotype[i]);
        }
        // make sure we don't forget the reference.  I'm sure there's a more elegant way to
        // do this, but it's fussy in light of pruning logic
        if (output_traversals.empty() &&
            !std::all_of(vcf_haplotype.begin(), vcf_haplotype.end(), [](int i) {return i == 0;})) {
            vector<int> ref_haplotype(vcf_haplotype.size(), 0);
            pair<SnarlTraversal, bool> alt_traversal = get_alt_traversal(
                site, site_variants, ref_path, start_step, end_step, ref_haplotype);
            assert(alt_traversal.second == true);
            output_traversals.push_back(make_pair(alt_traversal.first, ref_haplotype));
        }
        
        pair<SnarlTraversal, bool> alt_traversal = get_alt_traversal(site, site_variants, ref_path,
                                                                     start_step, end_step, vcf_haplotype);
#ifdef debug
        cerr << "bf haplotype <";
        for (auto allele : vcf_haplotype) {
            cerr << allele << ",";
        }
        cerr << "> gives " << (alt_traversal.second ? "valid" : "invalid") <<  " trav: "
             << pb2json(alt_traversal.first) << endl;
#endif
        if (alt_traversal.second) {
            output_traversals.push_back(make_pair(alt_traversal.first, vcf_haplotype));
        }
    } while (next_haplotype());   
}

pair <SnarlTraversal, bool> VCFTraversalFinder::get_alt_traversal(const Snarl& site,
                                                                  const vector<vcflib::Variant*>& site_variants,
                                                                  path_handle_t ref_path,
                                                                  step_handle_t start_step,
                                                                  step_handle_t end_step,
                                                                  const vector<int>& haplotype) {

    // Find the alt paths that we must cover if we traverse this haplotype
    pair<unordered_set<handle_t>, unordered_set<pair<handle_t, handle_t>>> alt_contents =
        get_haplotype_alt_contents(site_variants, haplotype, ref_path);
    unordered_set<handle_t>& alt_nodes = alt_contents.first;
    unordered_set<pair<handle_t, handle_t>>& alt_edges = alt_contents.second;

    // the edges of our reference path.  we must follow these in our traversal
    // unless we're going to an alt
    unordered_set<pair<handle_t, handle_t> > ref_edges;
    // nodes of the reference path.  we use these to see if we can go from
    // an alt path back to the reference
    unordered_set<handle_t> ref_nodes;
    for (auto step = start_step; step != end_step; step = graph.get_next_step(step)) {
        auto next = graph.get_next_step(step);

        // todo: assuming forward ref path
        ref_edges.insert(graph.edge_handle(graph.get_handle_of_step(step),
                                           graph.get_handle_of_step(next)));

        ref_nodes.insert(graph.get_handle_of_step(step));
    }
    ref_nodes.insert(graph.get_handle_of_step(end_step));

#ifdef debug
    cerr << "  alt nodes: ";
    for (auto alt_node : alt_nodes) {
        cerr << graph.get_id(alt_node) << ":" << graph.get_is_reverse(alt_node) << ",";
    }
    cerr << endl << "  alt edges: ";
    for (auto alt_edge : alt_edges) {
        cerr << graph.get_id(alt_edge.first) << ":" << graph.get_is_reverse(alt_edge.first) << "-"
             << graph.get_id(alt_edge.second) << ":" << graph.get_is_reverse(alt_edge.second) << ",";
    }
    cerr << endl << "  ref edges: ";
    for (auto ref_edge : ref_edges) {
        cerr << graph.get_id(ref_edge.first) << ":" << graph.get_is_reverse(ref_edge.first) << "-"
             << graph.get_id(ref_edge.second) << ":" << graph.get_is_reverse(ref_edge.second) << ",";
    }            
#endif
    
    // we walk by always following reference edges unless we can step into an alt
    // there are some simplifying assumptions about alt paths at play here, like
    // how there are unique hooks between them and the reference.
    bool in_alt_path = false;
    auto walk_forward = [&] (Visit& visit) {
        handle_t handle = graph.get_handle(visit.node_id(), visit.backward());
        // take an alt edge if we find one
        bool found_edge = !graph.follow_edges(handle, false, [&] (const handle_t& next) {
                auto edge = graph.edge_handle(handle, next);
                bool ret = true;
                if (alt_edges.count(edge)) {
                    ret = false; // stop, we found deletion edge
                } else if (alt_nodes.count(next)) {
                    in_alt_path = true;
                    ret = false; // stop, we found an edge to an alt path node
                }
                if (ret == false) {
                    // make sure we never cross this node/edge again in our traversal
                    alt_edges.erase(edge);
                    alt_nodes.erase(graph.get_handle(graph.get_id(next), false));
                    alt_nodes.erase(graph.get_handle(graph.get_id(next), true));
                    visit.set_node_id(graph.get_id(next));
                    visit.set_backward(graph.get_is_reverse(next));
                }
                return ret;
            });
        if (!found_edge) {
            // no alt edge found, take a reference edge
            found_edge = !graph.follow_edges(handle, false, [&] (const handle_t& next) {
                auto edge = graph.edge_handle(handle, next);
                bool ret = true;
                if (ref_edges.count(edge) && ref_nodes.count(next)) {
                    ret = false; // stop, we found a reference edge
                } else if (in_alt_path && ref_nodes.count(next)) {
                    in_alt_path = false;
                    ret = false; // stop, we found a reference node after our alt path
                }
                if (ret == false) {
                    // make sure we never cross this node/edge again in our traversal
                    ref_edges.erase(edge);
                    ref_nodes.erase(graph.get_handle(graph.get_id(next), false));
                    ref_nodes.erase(graph.get_handle(graph.get_id(next), true));
                    visit.set_node_id(graph.get_id(next));
                    visit.set_backward(graph.get_is_reverse(next));
                }
                return ret;
            });
        }
        return found_edge;
    };

    
    Visit visit;
    SnarlTraversal traversal;
    
    // start at the start
    // todo: should make sure this works if our snarl is backward on reference
    visit.set_node_id(graph.get_id(graph.get_handle_of_step(start_step)));
    visit.set_backward(graph.get_is_reverse(graph.get_handle_of_step(start_step)));
    ref_nodes.erase(graph.get_handle(visit.node_id(), false));
    ref_nodes.erase(graph.get_handle(visit.node_id(), true));
    

    if (include_endpoints) {
        *traversal.add_visit() = visit;
    }

#ifdef debug
    cerr << "  start walk: " << pb2json(visit) << endl;
#endif
    
    // walk our traversal
    bool found_end;
    while (walk_forward(visit)) {
#ifdef debug
        cerr << "  visit: " << pb2json(visit) << endl;
#endif
        if (visit.node_id() != graph.get_id(graph.get_handle_of_step(end_step))) {
            *traversal.add_visit() = visit;
        } else {
            found_end = true;
            break;
        }
    }

    if (include_endpoints) {
        Visit* visit = traversal.add_visit();
        handle_t end_handle = graph.get_handle_of_step(end_step);
        visit->set_node_id(graph.get_id(end_handle));
        visit->set_backward(graph.get_is_reverse(end_handle));
    }

    // sanity check: we compare the output to something gotten directly from the
    // path index when doing the reference haplotype.
    if (all_of(haplotype.begin(), haplotype.end(), [] (int i) {return i == 0;})) { 
        SnarlTraversal ref_trav;
        step_handle_t step = graph.get_next_step(start_step);
        if (include_endpoints) {
            Visit* visit = ref_trav.add_visit();
            visit->set_node_id(graph.get_id(graph.get_handle_of_step(start_step)));
            visit->set_backward(graph.get_is_reverse(graph.get_handle_of_step(start_step)));
        }
        for (; step != end_step; step = graph.get_next_step(step)) {
            Visit* visit = ref_trav.add_visit();
            visit->set_node_id(graph.get_id(graph.get_handle_of_step(step)));
            // todo: do we get an orientation out of the path index?  
            visit->set_backward(graph.get_is_reverse(graph.get_handle_of_step(step)));
        }
        if (include_endpoints) {
            Visit* visit = ref_trav.add_visit();
            visit->set_node_id(graph.get_id(graph.get_handle_of_step(end_step)));
            visit->set_backward(graph.get_is_reverse(graph.get_handle_of_step(end_step)));
        }
        assert(found_end && ref_trav == traversal);
    }
    
    return make_pair(traversal, found_end && alt_nodes.empty() && alt_edges.empty());    
}

pair<unordered_set<handle_t>, unordered_set<pair<handle_t, handle_t> > >
VCFTraversalFinder::get_haplotype_alt_contents(
    const vector<vcflib::Variant*>& site_variants,
    const vector<int>& haplotype,
    path_handle_t ref_path) {

    assert(haplotype.size() == site_variants.size());

    unordered_set<handle_t> alt_nodes;
    unordered_set<pair<handle_t, handle_t> > alt_deletion_edges;

    for (size_t allele = 0; allele < haplotype.size(); ++allele) {
        // ignore reference alleles
        if (haplotype[allele] == 0) {
            continue;
        }
        vcflib::Variant* var = site_variants[allele];

        // get the alt path information out of the graph
        pair<SnarlTraversal, vector<edge_t>> alt_path_info = get_alt_path(var, haplotype[allele], ref_path);
        if (alt_path_info.first.visit_size() == 0) {
            // skip deletion alt path where we can't find the deletion edge in the graph
            continue;
        }
        SnarlTraversal& alt_traversal = alt_path_info.first;
        bool is_deletion = !alt_path_info.second.empty();

        if (!is_deletion) {
            // fill in the nodes from the path
            for (size_t i = 0; i < alt_traversal.visit_size(); ++i) {
                alt_nodes.insert(graph.get_handle(alt_traversal.visit(i).node_id(),
                                                  alt_traversal.visit(i).backward()));
            }
        }  else {
            // add the deletion edges from the path
            for (auto deletion_edge : alt_path_info.second) {
                alt_deletion_edges.insert(deletion_edge);
            }
        }
    }

    return make_pair(alt_nodes, alt_deletion_edges);
}

pair<SnarlTraversal, vector<edge_t>> VCFTraversalFinder::get_alt_path(vcflib::Variant* var, int allele,
                                                            path_handle_t ref_path) {

    SnarlTraversal alt_path;
    vector<edge_t> deletion_edges;
    
    string alt_path_name = "_alt_" + make_variant_id(*var) + "_" + to_string(allele);
    if (graph.has_path(alt_path_name) && !graph.is_empty(graph.get_path_handle(alt_path_name))) {
        // if there's an alt path, then we're dealing with a snp or insertion.
        // we take the edges from the path, as well as those back to the reference
        for (handle_t handle : graph.scan_path(graph.get_path_handle(alt_path_name))) {
            // fill in the nodes from the path
            Visit* visit = alt_path.add_visit();
            visit->set_node_id(graph.get_id(handle));
            visit->set_backward(graph.get_is_reverse(handle));
        }
    }  else {
        // there's no alt path, it must be a deletion (if our input allele != 0)
        // in this case we use the reference allele path, and try to find an edge that spans
        // it.  this will be our alt edge
        // todo: put an alt path name maker into utility.hpp
        bool is_deletion = allele != 0;
        alt_path_name = "_alt_" + make_variant_id(*var) + "_0";
        // allele 0 can be empty for an insertion.  we don't complain if it's not in the graph
        assert(allele == 0 || graph.has_path(alt_path_name));

        if (graph.has_path(alt_path_name)) {
            path_handle_t path_handle = graph.get_path_handle(alt_path_name);            
            if (!graph.is_empty(path_handle)) {
                // find where this path begins and ends in the reference path index
                auto first_step_found = step_in_path(graph.get_handle_of_step(graph.path_begin(path_handle)), ref_path);
                assert(first_step_found.second);
                step_handle_t first_step = first_step_found.first;
                auto last_step_found = step_in_path(graph.get_handle_of_step(graph.path_back(path_handle)), ref_path);
                assert(last_step_found.second);
                step_handle_t last_step = last_step_found.first;
                
                // todo: logic needed here if want to support non-forward reference paths.
                first_step = graph.get_previous_step(first_step);
                last_step = graph.get_next_step(last_step);
                if (allele == 0) {
                    handle_t left = graph.get_handle_of_step(first_step);
                    handle_t right = graph.get_handle_of_step(last_step);
                    
                    Visit* visit = alt_path.add_visit();
                    visit->set_node_id(graph.get_id(left));
                    visit->set_backward(graph.get_is_reverse(left));
                    visit = alt_path.add_visit();
                    visit->set_node_id(graph.get_id(right));
                    visit->set_backward(graph.get_is_reverse(right));
                } else {
                    // alt paths don't always line up with deletion edges, so we hunt for
                    // our deletion edge using the path_index here.
                    pair<SnarlTraversal, vector<edge_t>> scan_deletion = scan_for_deletion(var, allele, ref_path,
                                                                                           first_step, last_step);
                    if (scan_deletion.first.visit_size() == 0) {
                        cerr << "[VCFTraversalFinder] Warning: Could not find deletion edge that matches allele "
                             << allele << " of\n" << *var << "\naround alt path" << alt_path_name << ":";
                    }
                    alt_path = std::move(scan_deletion.first);
                    deletion_edges = std::move(scan_deletion.second);
                }
            } else {
                assert(allele == 0);
            }
        }
    }

    return make_pair(alt_path, deletion_edges);
}

pair<SnarlTraversal, vector<edge_t>> VCFTraversalFinder::scan_for_deletion(vcflib::Variant* var, int allele, path_handle_t ref_path,
                                                                           step_handle_t first_path_step, step_handle_t last_path_step) {
    assert(allele > 0);

    // if our path matches an edge, we don't need to do anything
    edge_t spanning_edge = graph.edge_handle(graph.get_handle_of_step(first_path_step),
                                             graph.get_handle_of_step(last_path_step));
    if (graph.has_edge(spanning_edge)) {
        SnarlTraversal traversal;
        Visit* visit = traversal.add_visit();
        visit->set_node_id(graph.get_id(graph.get_handle_of_step(first_path_step)));
        visit->set_backward(graph.get_is_reverse(graph.get_handle_of_step(first_path_step)));
        visit = traversal.add_visit();
        visit->set_node_id(graph.get_id(graph.get_handle_of_step(last_path_step)));
        visit->set_backward(graph.get_is_reverse(graph.get_handle_of_step(last_path_step)));
        return make_pair(traversal, vector<edge_t>(1, spanning_edge));
    }

    // we're doing everything via length comparison, so keep track of the length we're
    // looking for (vcf_deletion_length) and the length we have (path_deletion_length)
    int vcf_deletion_length = var->alleles[0].length() - var->alleles[allele].length();

    // make our search window by scanning out the ends of our path
    step_handle_t first_window_step = first_path_step;
    for (int i = 0; i < max_deletion_scan_nodes && graph.has_previous_step(first_window_step); ++i) {
        first_window_step = graph.get_previous_step(first_window_step);
    }
    step_handle_t last_window_step = last_path_step;
    for (int i = 0; i < max_deletion_scan_nodes && graph.has_next_step(last_window_step); ++i) {
        last_window_step = graph.get_next_step(last_window_step);
    }
    
    // index our reference offsets (assuming forward reference path with no cycles)
    // needing this logic doesn't happen very often, otherwise would consider
    // requiring path position interface
    unordered_map<nid_t, int> ref_offsets;
    unordered_map<nid_t, step_handle_t> node_to_step;
    int offset = 0;
    int first_offset = 0;
    int last_offset = 0;
    step_handle_t cur_step = first_window_step;
    while (true) { // want to iterate last-inclusive
        handle_t cur_handle = graph.get_handle_of_step(cur_step);
        assert(graph.get_is_reverse(cur_handle) == false);
        ref_offsets[graph.get_id(cur_handle)] = offset;
        node_to_step[graph.get_id(cur_handle)] = cur_step;
        if (cur_step == first_path_step) {
            first_offset = offset + graph.get_length(cur_handle);
        }
        if (cur_step == last_path_step) {
            last_offset = offset;
        }
        offset += graph.get_length(cur_handle);
        if (cur_step == last_window_step) {
            break;
        } else {
            cur_step = graph.get_next_step(cur_step);
        }
    }

    // find our deletions, and index them by how close they match our given alt path
    // the delta is the min length of the deletion's two endpoints to the paths enpoints. 
    multimap<int, edge_t> delta_to_deletion; // index deleions by their distance from ref path ends
    unordered_map<edge_t, size_t> deletion_to_length; // lookup deletion sizes
    for (cur_step = first_window_step; cur_step != last_window_step; cur_step = graph.get_next_step(cur_step)) {
        handle_t cur_handle = graph.get_handle_of_step(cur_step);
        handle_t next_handle = graph.get_handle_of_step(graph.get_next_step(cur_step));        
        graph.follow_edges(cur_handle, false, [&] (const handle_t& edge_next_handle) {
                if (!graph.get_is_reverse(edge_next_handle) && // ignore inversions
		    graph.get_id(next_handle) != graph.get_id(edge_next_handle) &&
                    ref_offsets.count(graph.get_id(edge_next_handle))) {
                    // we are in a deletion that's contained in the window
                    int deletion_start_offset = ref_offsets[graph.get_id(cur_handle)] + graph.get_length(cur_handle);
                    int deletion_end_offset = ref_offsets[graph.get_id(edge_next_handle)];
                    int delta = std::min(std::abs(deletion_start_offset - first_offset), std::abs(deletion_end_offset - last_offset));
                    delta_to_deletion.insert(make_pair(delta, make_pair(cur_handle, edge_next_handle)));
                    deletion_to_length[make_pair(cur_handle, edge_next_handle)] = deletion_end_offset - deletion_start_offset;
                }
            });
    }

    // our goal is to find a traversal that threads the deletions we find.  in order to do this, our deltions
    // can't overlap
    function<bool(edge_t, const vector<edge_t>&)> doesnt_intersect = [&](edge_t edge, const vector<edge_t>& edge_set) {
        int edge_start = ref_offsets[graph.get_id(edge.first)] + graph.get_length(edge.first);
        int edge_end = ref_offsets[graph.get_id(edge.second)];
        // because of previous assumptins, are edges shoudl always be forward alogn the path.
        // note they are not made wtih graph.edge_handle, and are just oriented in the scan order
        assert(edge_start <= edge_end);
        for (edge_t other_edge : edge_set) {
            int other_start = ref_offsets[graph.get_id(other_edge.first)] + graph.get_length(other_edge.first);
            int other_end = ref_offsets[graph.get_id(other_edge.second)];
            if ((other_start >= edge_start && other_start < edge_end) || (other_end > edge_start && other_end <= edge_end) ||
                (edge_start >= other_start && edge_start < other_end) || (edge_end > other_start && edge_end >= other_end)) {
                return false;
            }
        }
        return true;
    };
    
    // greedily try to find some deletions that add up to the desired length, and are close to spanning
    // the alt path
    int best_delta = numeric_limits<int>::max();
    vector<edge_t> best_set;
    for (auto delta_edge : delta_to_deletion) {
        // can do better than quadratic here, but the sizes should be small enough not to matter
        vector<edge_t> candidate_set = {delta_edge.second};
        size_t total_size = deletion_to_length[delta_edge.second];
        int total_delta = delta_edge.first;
        for (auto delta_edge2 : delta_to_deletion) {
            if (delta_edge2 != delta_edge) {
                if (total_size + deletion_to_length[delta_edge2.second] <= vcf_deletion_length &&
                    // make that cubic...
                    doesnt_intersect(delta_edge2.second, candidate_set)) {
                        total_size += deletion_to_length[delta_edge2.second];
                        total_delta += delta_edge2.first;
                        candidate_set.push_back(delta_edge2.second);
                    }
            }
            if (total_size == vcf_deletion_length) {
                break;
            }
        }
        if (total_delta < best_delta) {
            best_set = candidate_set;
            best_delta = total_delta;
        }
    }

    // sort the edges along the path
    std::sort(best_set.begin(), best_set.end(), [&](edge_t e1, edge_t e2) {
            return ref_offsets[graph.get_id(e1.first)] < ref_offsets[graph.get_id(e2.first)]; });

    SnarlTraversal traversal;
    Visit* visit;
    // fill out the traversal
    for (int i = 0; i < best_set.size(); ++i) {
        // add a visit for each edge endpoint
        if (i == 0 || best_set[i].first != best_set[i-1].second) {
            visit = traversal.add_visit();
            visit->set_node_id(graph.get_id(best_set[i].first));
            visit->set_backward(graph.get_is_reverse(best_set[i].first));
        }
        visit = traversal.add_visit();
        visit->set_node_id(graph.get_id(best_set[i].second));
        visit->set_backward(graph.get_is_reverse(best_set[i].second));
        // the fill in the reference path to the next edge
        if (i < best_set.size() - 1) {
            step_handle_t next_step = node_to_step[graph.get_id(best_set[i + 1].first)];
            step_handle_t cur_step = node_to_step[graph.get_id(best_set[i].second)];
            if (cur_step != next_step) {
                for (cur_step = graph.get_next_step(cur_step); cur_step != next_step; cur_step = graph.get_next_step(cur_step)) {
                    visit = traversal.add_visit();
                    visit->set_node_id(graph.get_id(graph.get_handle_of_step(cur_step)));
                    visit->set_backward(graph.get_is_reverse(graph.get_handle_of_step(cur_step)));
                }
            }
        }
    }

    return make_pair(traversal, best_set);
}


vector<vector<int>> VCFTraversalFinder::get_pruned_alt_alleles(
    const Snarl& site,
    const vector<vcflib::Variant*>& site_variants,
    path_handle_t ref_path) {

    vector<vector<int> > alt_alleles(site_variants.size());

    for (int var_i = 0; var_i < site_variants.size(); ++var_i) {
        for (int allele = 0; allele < site_variants[var_i]->alleles.size(); ++allele) {
            alt_alleles[var_i].push_back(allele);
        }
    }

    // only invoke pruning if we exceed our cutoff.  fairly rare on most graphs
    for (int prune_it = 0; prune_it < max_prune_iterations && !check_max_trav_cutoff(alt_alleles); ++prune_it) {
        for (auto& alleles : alt_alleles) {
            alleles.clear();
        }
        
        for (int var_i = 0; var_i < site_variants.size(); ++var_i) {
            for (int allele = 0; allele < site_variants[var_i]->alleles.size(); ++allele) {
                if (skip_alt == nullptr ||
                    skip_alt(get_alt_path(site_variants[var_i], allele, ref_path).first, prune_it) == false) {
                    alt_alleles[var_i].push_back(allele);
                }
#ifdef debug
                else {
                    cerr << "Pruning allele " << allele << " from variant " << site_variants[var_i]->id << endl;
                }
#endif
            }
            // always leave at least one path through the site, even if that means
            // going through a reference allele that fails the skip_alt check.
            if (alt_alleles[var_i].empty()) {
                alt_alleles[var_i].push_back(0);
            }
        }
    }

    return alt_alleles;
}

bool VCFTraversalFinder::check_max_trav_cutoff(const vector<vector<int> >& alleles) {
    if (alleles.empty()) {
        return true;
    }    
    size_t count = 1;

    for (int i = 0; i < alleles.size(); ++i) {
        count *= alleles[i].size();
        if (count > max_traversal_cutoff) {
            return false;
        }
    }

    return true;
}

pair<step_handle_t, bool> VCFTraversalFinder::step_in_path(handle_t handle, path_handle_t path_handle) const {
    vector<step_handle_t> steps = graph.steps_of_handle(handle);
    // must be a cyclic!
    for (auto step : steps) {
        if (graph.get_path_handle_of_step(step) == path_handle) {
            return make_pair(step, true);
        }
    }
    return make_pair(step_handle_t(), false);
}


FlowTraversalFinder::FlowTraversalFinder(const HandleGraph& graph, SnarlManager& snarl_manager,
                                         size_t K,
                                         function<double(handle_t)> node_weight_callback,
                                         function<double(edge_t)> edge_weight_callback) :
    graph(graph),
    snarl_manager(snarl_manager),
    K(K),
    node_weight_callback(node_weight_callback),
    edge_weight_callback(edge_weight_callback)  {
    
}

void FlowTraversalFinder::setK(size_t k) {
    K = k;
}

vector<SnarlTraversal> FlowTraversalFinder::find_traversals(const Snarl& site) {
    return find_weighted_traversals(site).first;
}

pair<vector<SnarlTraversal>, vector<double>> FlowTraversalFinder::find_weighted_traversals(const Snarl& site, bool greedy_avg,
                                                                                           const HandleGraph* overlay) {

    // option to use the overlay graph for the search
    const HandleGraph* use_graph = overlay != nullptr ? overlay : & graph;
    
    handle_t start_handle = use_graph->get_handle(site.start().node_id(), site.start().backward());
    handle_t end_handle = use_graph->get_handle(site.end().node_id(), site.end().backward());
    
    vector<pair<double, vector<handle_t>>> widest_paths = algorithms::yens_k_widest_paths(use_graph, start_handle, end_handle, K,
                                                                                          node_weight_callback,
                                                                                          edge_weight_callback,
                                                                                          greedy_avg);

    vector<SnarlTraversal> travs;
    travs.reserve(widest_paths.size());
    vector<double> weights;
    weights.reserve(widest_paths.size());

    for (const auto& wp : widest_paths) {
        weights.push_back(wp.first);
        travs.emplace_back();
        for (const auto& h : wp.second) {
            Visit* visit = travs.back().add_visit();
            visit->set_node_id(use_graph->get_id(h));
            visit->set_backward(use_graph->get_is_reverse(h));
        }
    }

    return make_pair(travs, weights);
}        

GBWTTraversalFinder::GBWTTraversalFinder(const HandleGraph& graph, const gbwt::GBWT& gbwt) : 
    graph(graph),
    gbwt(gbwt) {

}
    
GBWTTraversalFinder::~GBWTTraversalFinder() {

}

pair<vector<SnarlTraversal>, vector<vector<gbwt::size_type>>>
GBWTTraversalFinder::find_gbwt_traversals(const Snarl& site, bool return_paths) {

    // follow all gbwt threads from start to end
    vector<pair<vector<gbwt::node_type>, gbwt::SearchState> > forward_traversals = list_haplotypes(
        graph,
        gbwt,                                                                                                   
        graph.get_handle(site.start().node_id(), site.start().backward()),
        [&] (const vector<gbwt::node_type>& new_thread) {
            return gbwt::Node::id(new_thread.back()) == site.end().node_id() &&
            gbwt::Node::is_reverse(new_thread.back()) == site.end().backward();
        });

    // follow all gbwt threads from end to start
    vector<pair<vector<gbwt::node_type>, gbwt::SearchState> > backward_traversals;
    if (!gbwt.bidirectional()) {
        backward_traversals = list_haplotypes(
            graph,
            gbwt,
            graph.get_handle(site.end().node_id(), !site.end().backward()),
            [&] (const vector<gbwt::node_type>& new_thread) {
                return gbwt::Node::id(new_thread.back()) == site.start().node_id() &&
                gbwt::Node::is_reverse(new_thread.back()) == !site.start().backward();
            });
    }

    // store them all as snarltraversals
    vector<SnarlTraversal> traversals;
    vector<vector<gbwt::size_type>> gbwt_paths;
    traversals.reserve(forward_traversals.size() + backward_traversals.size());

    // copy the forward traversals from gbwt vectors to snarl traversals
    for (int i = 0; i < forward_traversals.size(); ++i) {
        traversals.emplace_back();
        for (auto j = forward_traversals[i].first.begin(); j != forward_traversals[i].first.end(); ++j) {
            Visit* visit = traversals.back().add_visit();
            *visit = to_visit(gbwt::Node::id(*j), gbwt::Node::is_reverse(*j));
        }
        if (return_paths) {
            gbwt_paths.push_back(gbwt.locate(forward_traversals[i].second));
        }
    }

    if (!backward_traversals.empty()) {

        // want to check we don't have the same element twice
        std::sort(forward_traversals.begin(), forward_traversals.end(),
                  [&](const pair<vector<gbwt::node_type>, gbwt::SearchState>& t1,
                      const pair<vector<gbwt::node_type>, gbwt::SearchState>& t2) {
                      return t1.first < t2.first; });
        
        // copy and reverse the backward traversals into the snarl traversals
        for (int i = 0; i < backward_traversals.size(); ++i) {

            vector<gbwt::size_type> gbwt_path;
            if (return_paths) {
                gbwt_path = gbwt.locate(backward_traversals[i].second);
            }
            
            // orient along the snarl
            std::reverse(backward_traversals[i].first.begin(), backward_traversals[i].first.end());
            for (auto& gnode : backward_traversals[i].first) {
                gnode = gbwt::Node::encode(gbwt::Node::id(gnode), !gbwt::Node::is_reverse(gnode));
            }

            // search in the forward traversals
            auto si = std::lower_bound(forward_traversals.begin(), forward_traversals.end(), backward_traversals[i],
                                       [&](const pair<vector<gbwt::node_type>, gbwt::SearchState>& t1,
                                           const pair<vector<gbwt::node_type>, gbwt::SearchState>& t2) {
                                           return t1.first < t2.first; });
            if (si != forward_traversals.end() && si->first == backward_traversals[i].first) {
                // we found and exact forward match, just add in the paths
                if (return_paths) {
                    size_t idx = si - forward_traversals.begin();
                    gbwt_paths[idx].insert(gbwt_paths[idx].end(), gbwt_path.begin(), gbwt_path.end());
                }
            } else {
                // insert if not duplicate of existing forward traversal
                traversals.emplace_back();
                for (auto j = backward_traversals[i].first.begin(); j != backward_traversals[i].first.end(); ++j) {
                    Visit* visit = traversals.back().add_visit();
                    *visit = to_visit(gbwt::Node::id(*j), gbwt::Node::is_reverse(*j));
                }
                if (return_paths) {
                    gbwt_paths.push_back(gbwt.locate(backward_traversals[i].second));
                }
            }
        }
    }
    return make_pair(traversals, gbwt_paths);
}

vector<SnarlTraversal> GBWTTraversalFinder::find_traversals(const Snarl& site) {
    return find_gbwt_traversals(site, false).first;
}

pair<vector<SnarlTraversal>, vector<gbwt::size_type>> GBWTTraversalFinder::find_path_traversals(const Snarl& site) {
    // get the unique traversals
    pair<vector<SnarlTraversal>, vector<vector<gbwt::size_type>>> gbwt_traversals = find_gbwt_traversals(site, true);

    // expand them out to one per path (this is to be consistent with PathTraversalFinder as used in deconstruct)
    pair<vector<SnarlTraversal>, vector<gbwt::size_type>> path_traversals;
    for (size_t i = 0; i < gbwt_traversals.first.size(); ++i) {
        SnarlTraversal& trav = gbwt_traversals.first[i];
        vector<gbwt::size_type>& paths = gbwt_traversals.second[i];
        for (size_t j = 0; j < paths.size(); ++j) {
            path_traversals.first.push_back(trav);
            path_traversals.second.push_back(paths[j]);
        }
    }
    
    return path_traversals;
}

}



