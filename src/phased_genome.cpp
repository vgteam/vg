//
//  phased_genome.cpp
//  
//
//  Created by Jordan Eizenga on 10/14/16.
//
//

#include "phased_genome.hpp"


namespace vg; {
    
    HaplotypeNode::HaplotypeNode(NodeTraversal node_traversal, HaplotypeNode* next, HaplotypeNode* prev) :
                                 node_traversal(node_traversal), next(next), prev(prev)
    {
        // nothing to do
    }
    
    HaplotypeNode::~HaplotypeNode() {
        // nothing to do
    }
    
    Haplotype::Haplotype(NodeTraversal node_traversal) {
        // construct seed node
        left_telomere_node = new HaplotypeNode(node_traversal, nullptr, nullptr);
        right_telomere_node = left_telomere_node;
    }
    
    Haplotype::Haplotype(list<NodeTraversal>& node_path) {
        if (node_path.empty()) {
            cerr << "error:[Haplotype] cannot construct a haplotype with an empty list of nodes" << endl;
            assert(0);
        }
        
        // construct seed node
        auto iter = node_path.begin();
        left_telomere_node = left_telomere_node = new HaplotypeNode(*iter, nullptr, nullptr);
        right_telomere_node = left_telomere_node;
        iter++;
        
        // add subsequent nodes
        for (; iter != node_path.end(); iter++) {
            append_right(*iter);
        }
    }
    
    Haplotype::~Haplotype() {
        
        // destruct nodes
        HaplotypeNode* haplo_node = left_telomere_node;
        while (haplo_node) {
            HaplotypeNode* next = haplo_node->next;
            delete haplo_node;
            haplo_node = next;
        }
    }
    
    HaplotypeNode* Haplotype::append_left(NodeTraversal node_traversal) {
        left_telomere_node = new HaplotypeNode(node_traversal, left_telomere_node, nullptr);
        left_telomere_node->next->prev = left_telomere_node
        return left_telomere_node;
    }
    
    HaplotypeNode* Haplotype::append_right(NodeTraversal node_traversal) {
        right_telomere_node = new HaplotypeNode(node_traversal, nullptr, right_telomere_node);
        right_telomere_node->prev->next = right_telomere_node;
        return right_telomere_node;
    }
    
    PhasedGenome::PhasedGenome(const VG& graph) : graph(graph) {
        // nothing to do before
    }
    
    template <typename NodeTraversalIterator>
    void PhasedGenome::add_haplotype(NodeTraversalIterator first, NodeTraversalIterator last) {
        
        if (first == last) {
            cerr << "error:[PhasedGenome] cannot insert an empty haplotype into genome"
            assert(0);
        }
        
        // initialize haplotype with one node
        haplotypes.emplace_back(*first);
        Haplotype& haplotype = haplotypes.back();
        first++;
        
        // add each subsequence node
        for (; first != last; first++) {
            haplotype.append_right(*first);
        }
    }
    
    void PhasedGenome::build_site_indices(vector<NestedSite>& top_level_sites) {
        
        // construct the start and end of site indices
        for (NestedSite& site : top_level_sites) {
            build_site_indices_internal(site);
        }
        
        // label the sites on each haplotype
        for (Haplotype& haplotype : haplotypes) {
            
            // keep track of where we enter and leave sites
            unordered_map<NestedSite, HaplotypeNode*, HashSite> site_sides;
            
            // iterate along the node path of the haplotype
            HaplotypeNode* haplo_node = haplotype.left_telomere_node;
            while (haplo_node != nullptr) {
                int64_t node_id = haplo_node->node_traversal.node->id();
                
                // mark this instance of the node in the node location index
                node_locations[node_id].push_back(haplo_node);
                
                // are we at the start or end of a site?
                NestedSite* site = nullptr;
                if (site_starts.count(node_id)) {
                    site = site_starts[node_id];
                }
                else if (site_ends.count(node_id)) {
                    site = site_ends[node_id];
                }
                
                if (site) {
                    // are we leaving or entering the site?
                    if (site_sides.count(*site)) {
                        // leaving: put the site in the index in the orientation of haplotype travesal
                        HaplotypeNode* other_side_node = site_sides[*site];
                        site_sides.erase(*site);
                        sites[*site] = make_pair(other_side_node, haplo_node);
                    }
                    else {
                        // entering: mark the node in the haplotype path where we entered
                        site_sides[*subsite] = haplo_node;
                    }
                }
                
                // advance to next node in haplotype path
                haplo_node = haplo_node->next;
            }
        }
    }
    
    void PhasedGenome::build_site_indices_internal(NestedSite& site) {
        
        // record start and end of site
        site_starts[site.start->id()] == &site;
        site_ends[site.end->id()] == &site;
        
        // recurse through child sites
        for (NestedSite& subsite : site.children) {
            build_site_indices_internal(subsite);
        }
        
    }
    
    void PhasedGenome::insert_left(NodeTraversal node_traversal, HaplotypeNode* haplo_node) {
        
        HaplotypeNode* new_node = new HaplotypeNode(node_traversal, haplo_node, haplo_node->prev);
        haplo_node->prev->next = new_node;
        haplo_node->prev = new_node;
        
        node_locations[node->id()].push_back(new_node);
    }
    
    void PhasedGenome::insert_right(NodeTraversal node_traversal, HaplotypeNode* haplo_node) {
        
        HaplotypeNode* new_node = new HaplotypeNode(node_traversal, haplo_node->next, haplo_node);
        haplo_node->next->prev = new_node;
        haplo_node->next = new_node;
        
        node_locations[node->id()].push_back(new_node);
    }
    
    void PhasedGenome::swap_label(NestedSite& site, Haplotype& haplotype_1, Haplotype& haplotype_2) {
        
        // update index for this site
        auto temp = haplotype_1.sites[site];
        haplotype_1.sites[site] = haplotype_2.sites[site];
        haplotype_2.sites[site] = temp;
        
        // update index for child sites
        for (NestedSite& child_site : site.children) {
            // not all child sites will be in the haplotype at a given time
            if (haplotype_1.sites.count(NestedSite) && haplotype_2.sites.count(NestedSite)) {
                swap_label(child_site, haplotype_1, haplotype_2);
            }
        }
    }
    
    void PhasedGenome::swap_alleles(NestedSite& site, int haplotype_1, int haplotype_2) {
        Haplotype& haplo_1 = haplotypes[haplotype_1];
        Haplotype& haplo_2 = haplotypes[haplotype_2];
        
        pair<HaplotypeNode*, HaplotypeNode*> haplo_nodes_1 = haplo_1.sites[site];
        pair<HaplotypeNode*, HaplotypeNode*> haplo_nodes_2 = haplo_2.sites[site];
        
        // rewire the paths to switch
        HaplotypeNode* temp = haplo_nodes_1.first->prev;
        haplo_nodes_1.first->prev = haplo_nodes_2.first->prev;
        haplo_nodes_2.first->prev = temp;
        
        temp = haplo_nodes_1.second->next;
        haplo_nodes_1.first->next = haplo_nodes_2.second->next;
        haplo_nodes_2.second->next = temp;
        
        // update index for this site
        haplo_1.sites[site] = haplo_nodes_2;
        haplo_2.sites[site] = haplo_nodes_1;
        
        // update index for child sites
        for (NestedSite& child_site : site.children) {
            // not all child sites will be in the haplotype at a given time
            if (haplo_1.sites.count(NestedSite) && haplo_2.sites.count(NestedSite)) {
                swap_label(child_site, haplo_1, haplo_2);
            }
        }
    }
    
    
    // TODO: it seems like there should be a better way to do this than completely erasing
    // and then rewriting the site (especially at long sites)
    void PhasedGenome::set_allele(NestedSite& site, list<NodeTraversal>& allele, int which_haplotype) {
        
        Haplotype& haplotype = haplotypes[which_haplotype];
        
        // can only set the allele of a site that already is in the haplotype
        assert(haplotype.sites.count(site));
        
        pair<HaplotypeNode*, HaplotypeNode*> haplo_site = haplotype.sites[site];
        
        // remove the current site
        HaplotypeNode* haplo_node = haplo_site.first->next;
        while (haplo_node != haplo_site.second) {
            // don't need to worry about erasing the same site twice (once on start and
            // once on end) since unordered_map.erase is defined in both cases
            if (site_starts.count(haplo_node.node->id())) {
                NestedSite* subsite = site_starts[haplo_node.node->id()];
                haplotype.sites.erase(*subsite);
            }
            else if (site_ends.count(haplo_node.node->id())) {
                NestedSite* subsite = site_ends[haplo_node.node->id()];
                haplotype.sites.erase(*subsite);
            }
            
            // cut out each node
            remove(haplo_node);
        }
        
        // gives the haplotype node where we entered a site and thereby also
        // indicates whether we have entered the site yet
        unordered_map<NestedSite, HaplotypeNode*, HashSite> subsite_side;
        
        haplo_node = haplo_site.first;
        // is site in forward or reverse direction on haplotype?
        if (haplo_site.first->node_traversal.node == site.start.node) {
            for (auto iter = allele.begin(); iter != allele.end(); iter++) {
                
                // create the next node in the allele
                insert_right(*iter, haplo_node);
                // move the pointer onto the new node
                haplo_node = haplo_node->next;
                int64_t node_id = haplo_node.node_traversal.node->id();
                
                // does a site start or end here?
                if (site_starts.count(node_id)) {
                    NestedSite* subsite = site_starts[node_id];
                    // are we entering or leaving this site?
                    if (subsite_side.count(*subsite)) {
                        // add this site into the haplotype's site index
                        HaplotypeNode* end_haplo_node = subsite_side[*subsite];
                        haplotype.sites[*subsite] = make_pair(end_haplo_node, haplo_node);
                        // note: the site sides should always be entered in the order that
                        // they occur in the haplotype
                    }
                    else {
                        // we are entering a site, mark the location of the entrance
                        subsite_side[*subsite] = haplo_node;
                    }
                }
                else if (site_ends.count(node_id)) {
                    NestedSite* subsite = site_ends[node_id];
                    // are we entering or leaving this site?
                    if (subsite_side.count(*subsite)) {
                        // add this site into the haplotype's site index
                        HaplotypeNode* start_haplo_node = subsite_side[*subsite];
                        haplotype.sites[*subsite] = make_pair(start_haplo_node, haplo_node);
                    }
                    else {
                        // we are entering a site, mark the location of the entrance
                        subsite_side[*subsite] = haplo_node;
                    }
                }
                
            }
        }
        else {
            // TODO: reimplement this in a template so I don't need to have it repeated
            // for the two iterator types, or with advancing functions I can switch around
            for (auto iter = allele.rbegin(); iter != allele.rend(); iter++) {
                // create the next node in the allele
                insert_right(*iter, haplo_node);
                // move the pointer onto the new node
                haplo_node = haplo_node->next;
                int64_t node_id = haplo_node.node_traversal.node->id();
                
                // does a site start or end here?
                if (site_starts.count(node_id)) {
                    NestedSite* subsite = site_starts[node_id];
                    // are we entering or leaving this site?
                    if (subsite_side.count(*subsite)) {
                        // add this site into the haplotype's site index
                        HaplotypeNode* end_haplo_node = subsite_side[*subsite];
                        haplotype.sites[*subsite] = make_pair(end_haplo_node, haplo_node);
                        // note: the site sides should always be entered in the order that
                        // they occur in the haplotype
                    }
                    else {
                        // we are entering a site, mark the location of the entrance
                        subsite_side[*subsite] = haplo_node;
                    }
                }
                else if (site_ends.count(node_id)) {
                    NestedSite* subsite = site_ends[node_id];
                    // are we entering or leaving this site?
                    if (subsite_side.count(*subsite)) {
                        // add this site into the haplotype's site index
                        HaplotypeNode* start_haplo_node = subsite_side[*subsite];
                        haplotype.sites[*subsite] = make_pair(start_haplo_node, haplo_node);
                    }
                    else {
                        // we are entering a site, mark the location of the entrance
                        subsite_side[*subsite] = haplo_node;
                    }
                }
            }
        }
    }
    
    int32_t PhasedGenome::optimal_score_on_genome(const MultipathAlignment& multipath_aln) {
        
        int32_t optimal_score = 0;
        
        // find the places in the path where the alignment might start
        unordered_set<HaplotypeNode*> candidate_start_positions;
        for (int i = 0; i < multipath_aln.start_size(); i++) {
            // a starting subpath in the multipath alignment
            const Subpath& start_subpath = multipath_aln.subpath(multipath_aln.start(i));
            const Position& start_pos = start_subpath.path().mapping(0).position();
            
            // add each location the start nodes occur in the path to the candidate starts
            if (node_locations.count(start_pos.node_id())) {
                const vector<HaplotypeNode*>& subpath_matches = node_locations[start_pos.node_id()]
                candidate_start_positions.insert(subpath_matches.begin(), subpath_matches.end());
            }
        }
        
        // check alignments starting at each node in the path that has a source subpath starting on it
        for (HaplotypeNode* path_start_node : candidate_start_positions) {
            
            // match up forward and backward traversal on the path to forward and backward traversal through
            // the multipath alignment
            auto move_forward = [path](HaplotypeNode*& path_node) { path_node = path_node->next; };
            auto move_backward = [path](HaplotypeNode*& path_node) { path_node = path_node->prev; };
            
            if (start_pos.is_reverse() != path_node->node_traversal.backward) {
                swap(move_forward, move_backward);
            }
            
            // initialize 2 dynamic programming structures:
            // pointer to place in path corresponding to the beginning of a subpath
            vector<HaplotypeNode*> subpath_nodes = vector<HaplotypeNode*>(multipath_aln.subpath_size(), nullptr);
            // score of the best preceding path before this subpath
            vector<int32_t> subpath_prefix_score = vector<int32_t>(multipath_aln.subpath_size(), 0);
            
            // set DP base case with the subpaths that start at this path node
            for (int i = 0; i < multipath_aln.start_size(); i++) {
                Subpath* start_subpath = multipath_aln.subpath(multipath_aln.start(i));
                const Position& start_pos = first_path_position(start_subpath.path());
                
                if (start_pos.node_id() == (*path_start_node).node->id()) {
                    subpath_nodes[multipath_aln.start(i)] = path_start_node;
                }
            }
            
            for (int i = 0; i < multipath_aln.subpath_size(); i++) {
                HaplotypeNode* subpath_node = subpath_nodes[i];
                
                // this subpath may be unreachable from subpaths consistent with the path
                if (subpath_node == nullptr) {
                    continue;
                }
                
                const Subpath& subpath = multipath_aln.subpath(i);
                
                // iterate through mappings in this subpath (assumes one mapping per node)
                bool subpath_follows_path = true;
                for (int j = 0; j < subpath.path().mapping_size(); j++, move_forward(subpath_node)) {
                    // check if mapping corresponds to the next node in the path
                    if (subpath.path().mapping(j).position().node_id() != subpath_node->node_traversal.node->id()) {
                        subpath_follows_path = false;
                        break;
                    }
                }
                
                // update subsequent subpaths or record completed alignment if this subpath didn't fall off the path
                if (subpath_follows_path) {
                    int32_t extended_prefix_score = subpath_prefix_score[j] + subpath.score();
                    if (subpath.next_size() == 0) {
                        // reached a sink subpath (thereby completing an alignment), check for optimality
                        if (extended_prefix_score > optimal_score) {
                            optimal_score = extended_prefix_score;
                        }
                    }
                    else {
                        // edge case: check if subpath_node was improperly incremented from a mapping that ended in the
                        // middle of a node
                        Position end_pos = last_path_position(subpath.path());
                        if (end_pos.offset() != graph.get_node(end_pos.node_id())->sequence().length()) {
                            move_backward(subpath_node);
                        }
                        // TODO: this could be a problem if the next node is the end of a chromosome (will seg fault
                        // because can't get the previous node from nullptr)
                        
                        // mark which node the next subpath starts at
                        for (int j = 0; j < subpath.next_size(); j++) {
                            if (subpath_prefix_score[subpath.next(j)] < extended_prefix_score) {
                                subpath_prefix_score[subpath.next(j)] = extended_prefix_score;
                                subpath_nodes[subpath.next(j)] = subpath_node;
                            }
                        }
                    }
                }
            }
        }
        
        return optimal_score;
    }
}




