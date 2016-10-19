//
//  phased_genome.cpp
//

#include "phased_genome.hpp"

using namespace std;

namespace vg {
    
    PhasedGenome::HaplotypeNode::HaplotypeNode(NodeTraversal node_traversal, HaplotypeNode* next, HaplotypeNode* prev) :
                                 node_traversal(node_traversal), next(next), prev(prev) {
        // nothing to do
    }
    
    PhasedGenome::HaplotypeNode::~HaplotypeNode() {
        // nothing to do
    }
    
    PhasedGenome::Haplotype::Haplotype(NodeTraversal node_traversal) {
        // construct seed node
        left_telomere_node = new HaplotypeNode(node_traversal, nullptr, nullptr);
        right_telomere_node = left_telomere_node;
    }
    
    template <typename NodeTraversalIterator>
    PhasedGenome::Haplotype::Haplotype(NodeTraversalIterator first, NodeTraversalIterator last) {
        if (first == last) {
            cerr << "error:[PhasedGenome] cannot construct haplotype with 0 nodes" << endl;
            assert(0);
        }
        
        // construct seed node
        left_telomere_node = new HaplotypeNode(*first, nullptr, nullptr);
        right_telomere_node = left_telomere_node;
        first++;
        
        // add each subsequent node
        for (; first != last; first++) {
            append_right(*first);
        }
    }
    
    PhasedGenome::Haplotype::~Haplotype() {
        
        // destruct nodes
        HaplotypeNode* haplo_node = left_telomere_node;
        while (haplo_node) {
            HaplotypeNode* next = haplo_node->next;
            delete haplo_node;
            haplo_node = next;
        }
    }
    
    inline PhasedGenome::HaplotypeNode* PhasedGenome::Haplotype::append_left(NodeTraversal node_traversal) {
        left_telomere_node = new HaplotypeNode(node_traversal, left_telomere_node, nullptr);
        left_telomere_node->next->prev = left_telomere_node;
        return left_telomere_node;
    }
    
    inline PhasedGenome::HaplotypeNode* PhasedGenome::Haplotype::append_right(NodeTraversal node_traversal) {
        right_telomere_node = new HaplotypeNode(node_traversal, nullptr, right_telomere_node);
        right_telomere_node->prev->next = right_telomere_node;
        return right_telomere_node;
    }
    
    PhasedGenome::PhasedGenome(VG& graph) : graph(graph) {
        // nothing to do
    }
    
    template <typename NodeTraversalIterator>
    void PhasedGenome::add_haplotype(NodeTraversalIterator first, NodeTraversalIterator last) {
        
        haplotypes.emplace_back(first, last);
        
    }
    
    void PhasedGenome::build_indices(vector<NestedSite>& top_level_sites) {
        
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
                
                // are we at the start of a site?
                NestedSite* site = nullptr;
                if (site_starts.count(node_id)) {
                    site = site_starts[node_id];
                    // are we leaving or entering the site?
                    if (site_sides.count(*site)) {
                        // leaving: put the site in the index in the orientation of haplotype travesal
                        HaplotypeNode* other_side_node = site_sides[*site];
                        site_sides.erase(*site);
                        haplotype.sites[*site] = make_pair(other_side_node, haplo_node);
                    }
                    else {
                        // entering: mark the node in the haplotype path where we entered
                        site_sides[*site] = haplo_node;
                    }
                }
                
                // are we at the start of a site?
                if (site_ends.count(node_id)) {
                    site = site_ends[node_id];
                    // are we leaving or entering the site?
                    if (site_sides.count(*site)) {
                        // leaving: put the site in the index in the orientation of haplotype travesal
                        HaplotypeNode* other_side_node = site_sides[*site];
                        site_sides.erase(*site);
                        haplotype.sites[*site] = make_pair(other_side_node, haplo_node);
                    }
                    else {
                        // entering: mark the node in the haplotype path where we entered
                        site_sides[*site] = haplo_node;
                    }
                }
                
                // advance to next node in haplotype path
                haplo_node = haplo_node->next;
            }
        }
    }
    
    void PhasedGenome::build_site_indices_internal(NestedSite& site) {
        
        // record start and end of site
        site_starts[site.start.node->id()] == &site;
        site_ends[site.end.node->id()] == &site;
        
        // recurse through child sites
        for (NestedSite& subsite : site.children) {
            build_site_indices_internal(subsite);
        }
        
    }
    
    inline void PhasedGenome::insert_left(NodeTraversal node_traversal, HaplotypeNode* haplo_node) {
        
        HaplotypeNode* new_node = new HaplotypeNode(node_traversal, haplo_node, haplo_node->prev);
        haplo_node->prev->next = new_node;
        haplo_node->prev = new_node;
        
        node_locations[node_traversal.node->id()].push_back(new_node);
    }
    
    inline void PhasedGenome::insert_right(NodeTraversal node_traversal, HaplotypeNode* haplo_node) {
        
        HaplotypeNode* new_node = new HaplotypeNode(node_traversal, haplo_node->next, haplo_node);
        haplo_node->next->prev = new_node;
        haplo_node->next = new_node;
        
        node_locations[node_traversal.node->id()].push_back(new_node);
    }
    
    inline void PhasedGenome::remove(HaplotypeNode* haplo_node) {
        // note: updating site indices is handled outside this function
        
        haplo_node->next->prev = haplo_node->prev;
        haplo_node->prev->next = haplo_node->next;
        
        // remove the node from the node locations index
        int64_t node_id = haplo_node->node_traversal.node->id();
        list<HaplotypeNode*>& node_occurrences = node_locations[node_id];
        for (auto iter = node_occurrences.begin(); iter != node_occurrences.end(); iter++) {
            if (*iter == haplo_node) {
                node_occurrences.erase(iter);
            }
        }
        
        delete haplo_node;
    }
    
    void PhasedGenome::swap_label(NestedSite& site, Haplotype& haplotype_1, Haplotype& haplotype_2) {
        
        // update index for this site
        if (haplotype_1.sites.count(site)) {
            if (haplotype_2.sites.count(site)) {
                // site is in both haplotypes, swap the labels
                swap(haplotype_1.sites[site], haplotype_2.sites[site]);
            }
            else {
                // site is only in haplotype 1, transfer the label
                haplotype_2.sites[site] = haplotype_1.sites[site];
                haplotype_1.sites.erase(site);
            }
        }
        else if (haplotype_2.sites.count(site)) {
            // site is only in haplotype 2, transfer the label
            haplotype_1.sites[site] = haplotype_2.sites[site];
            haplotype_2.sites.erase(site);
        }
        else {
            // site is in neither haplotype
            return;
        }
        
        // update index for child sites
        for (NestedSite& child_site : site.children) {
            swap_label(child_site, haplotype_1, haplotype_2);
        }
    }
    
    void PhasedGenome::swap_alleles(NestedSite& site, int haplotype_1, int haplotype_2) {
        Haplotype& haplo_1 = haplotypes[haplotype_1];
        Haplotype& haplo_2 = haplotypes[haplotype_2];
        
        pair<HaplotypeNode*, HaplotypeNode*> haplo_nodes_1 = haplo_1.sites[site];
        pair<HaplotypeNode*, HaplotypeNode*> haplo_nodes_2 = haplo_2.sites[site];
        
        // rewire the paths to switch out the alleles
        haplo_nodes_1.first->prev->next = haplo_nodes_2.first;
        haplo_nodes_2.first->prev->next = haplo_nodes_1.first;
        haplo_nodes_1.second->next->prev = haplo_nodes_2.second;
        haplo_nodes_2.second->next->prev = haplo_nodes_1.second;
        swap(haplo_nodes_1.first->prev, haplo_nodes_2.first->prev);
        swap(haplo_nodes_1.second->next, haplo_nodes_2.second->next);
        
        // update index for this site
        haplo_1.sites[site] = haplo_nodes_2;
        haplo_2.sites[site] = haplo_nodes_1;
        
        // update index for child sites
        for (NestedSite& child_site : site.children) {
            swap_label(child_site, haplo_1, haplo_2);
        }
    }
    
    // TODO: it seems like there should be a better way to do this than completely erasing
    // and then rewriting the site (especially at long sites)
    template <typename NodeTraversalIterator>
    void PhasedGenome::set_allele(NestedSite& site, NodeTraversalIterator first, NodeTraversalIterator last,
                                  int which_haplotype) {
        
        Haplotype& haplotype = haplotypes[which_haplotype];
        
        // can only set the allele of a site that already is in the haplotype
        assert(haplotype.sites.count(site));
        
        pair<HaplotypeNode*, HaplotypeNode*> haplo_site = haplotype.sites[site];
        
        // remove the current site
        HaplotypeNode* haplo_node = haplo_site.first->next;
        while (haplo_node != haplo_site.second) {
            // don't need to worry about erasing the same site twice (once on start and
            // once on end) since unordered_map.erase is defined in both cases
            int64_t node_id = haplo_node->node_traversal.node->id();
            if (site_starts.count(node_id)) {
                NestedSite* subsite = site_starts[node_id];
                haplotype.sites.erase(*subsite);
            }
            else if (site_ends.count(node_id)) {
                NestedSite* subsite = site_ends[node_id];
                haplotype.sites.erase(*subsite);
            }
            // cut out each node and iterate to next node
            HaplotypeNode* next_haplo_node = haplo_node->next;
            remove(haplo_node);
            haplo_node = next_haplo_node;
        }
        
        // is site in forward or reverse direction on haplotype?
        bool forward = (haplo_site.first->node_traversal.node == site.start.node);
        
        // orient traversal through the allele to traversal along the haplotype
        haplo_node = forward ? haplo_site.first : haplo_site.second;
        auto oriented_insert = forward ? [this](NodeTraversal node_traversal, HaplotypeNode* haplo_node)
                                               {
                                                   this->insert_right(node_traversal, haplo_node);
                                                   return haplo_node->next;
                                               }
                                       : [this](NodeTraversal node_traversal, HaplotypeNode* haplo_node)
                                               {
                                                   this->insert_left(node_traversal, haplo_node);
                                                   return haplo_node->prev;
                                               };
        auto oriented_site = forward ? [](HaplotypeNode* haplo_node_1, HaplotypeNode* haplo_node_2)
                                         {
                                             return make_pair(haplo_node_1, haplo_node_2);
                                         }
                                     : [](HaplotypeNode* haplo_node_1, HaplotypeNode* haplo_node_2)
                                         {
                                             return make_pair(haplo_node_2, haplo_node_1);
                                         };
        
        // keeps track of the haplotype node where we entered a site and thereby also
        // indicates whether we have entered the site yet
        unordered_map<NestedSite, HaplotypeNode*, HashSite> subsite_side;
        
        // start inserting nodes from the first position in the allele
        for (; first != last; first++) {
            
            // create the next node in the allele and move the pointer onto the new node
            haplo_node = oriented_insert(*first, haplo_node);
            
            int64_t node_id = haplo_node->node_traversal.node->id();
            
            // does a site start here?
            if (site_starts.count(node_id)) {
                NestedSite* subsite = site_starts[node_id];
                // are we entering or leaving this site?
                if (subsite_side.count(*subsite)) {
                    // add this site into the haplotype's site index
                    HaplotypeNode* other_side_node = subsite_side[*subsite];
                    haplotype.sites[*subsite] = oriented_site(other_side_node, haplo_node);
                    // note: the site sides should always be entered in the order that
                    // they occur in the haplotype
                }
                else {
                    // we are entering a site, mark the location of the entrance
                    subsite_side[*subsite] = haplo_node;
                }
            }
            
            // does a site end here?
            if (site_ends.count(node_id)) {
                NestedSite* subsite = site_ends[node_id];
                // are we entering or leaving this site?
                if (subsite_side.count(*subsite)) {
                    // add this site into the haplotype's site index
                    HaplotypeNode* other_side_node = subsite_side[*subsite];
                    haplotype.sites[*subsite] = oriented_site(other_side_node, haplo_node);
                    // note: the site sides should always be entered in the order that
                    // they occur in the haplotype
                }
                else {
                    // we are entering a site, mark the location of the entrance
                    subsite_side[*subsite] = haplo_node;
                }
            }
        }
    }
    
    int32_t PhasedGenome::optimal_score_on_genome(const MultipathAlignment& multipath_aln) {
        
        // hash function for a traversal along a haplotype
        struct HashHaplotypeNodeTraversal {
            size_t operator()(const pair<HaplotypeNode*, bool>& haplo_node_traversal) const {
                return (size_t) 1099511628211ull * ((uintptr_t) haplo_node_traversal.first + haplo_node_traversal.second * 16777619ull) + 14695981039346656037ull;
            }
        };
        
        // iteration functions to facilitate iterating on forward/reverse strands
        auto move_right = [](HaplotypeNode*& path_node) { path_node = path_node->next; };
        auto move_left = [](HaplotypeNode*& path_node) { path_node = path_node->prev; };
        
        int32_t optimal_score = 0;
        
        // find the places in the path where the alignment might start
        unordered_map< pair<HaplotypeNode*, bool>, vector<int>, HashHaplotypeNodeTraversal> candidate_start_positions;
        for (int i = 0; i < multipath_aln.start_size(); i++) {
            // a starting subpath in the multipath alignment
            const Subpath& start_subpath = multipath_aln.subpath(multipath_aln.start(i));
            const Position& start_pos = start_subpath.path().mapping(0).position();
            
            // add each location the start nodes occur in the path to the candidate starts
            for ( HaplotypeNode* haplo_node : node_locations[start_pos.node_id()] ) {
                
                // mark the start locations orientation relative to the start node
                if (haplo_node->node_traversal.backward == start_pos.is_reverse()) {
                    candidate_start_positions[make_pair(haplo_node, true)].push_back(i);
                }
                else {
                    candidate_start_positions[make_pair(haplo_node, false)].push_back(i);
                }
            }
        }
        
        // check alignments starting at each node in the path that has a source subpath starting on it
        for (pair< pair<HaplotypeNode*, bool>, vector<int> > path_starts : candidate_start_positions) {
            
            HaplotypeNode* path_start_node = path_starts.first.first;
            bool oriented_forward = path_starts.first.second;
            const vector<int>& aln_starts = path_starts.second;
            
            // match up forward and backward traversal on the path to forward and backward traversal through
            // the multipath alignment
            auto move_forward = oriented_forward ? move_right : move_left;
            auto move_backward = oriented_forward ? move_left : move_right;
            
            // initialize dynamic programming structures:
            // pointer to place in haplotype path corresponding to the beginning of a subpath
            vector<HaplotypeNode*> subpath_nodes = vector<HaplotypeNode*>(multipath_aln.subpath_size(), nullptr);
            // score of the best preceding path before this subpath
            vector<int32_t> subpath_prefix_score = vector<int32_t>(multipath_aln.subpath_size(), 0);
            
            // set DP base case with the subpaths that start at this path node
            for (int i : aln_starts) {
                subpath_nodes[multipath_aln.start(i)] = path_start_node;
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
                    // check if mapping corresponds to the next node in the path in the correct orientation
                    const Position& position = subpath.path().mapping(j).position();
                    if (position.node_id() != subpath_node->node_traversal.node->id()
                        || ((position.is_reverse() == subpath_node->node_traversal.backward) != oriented_forward)) {
                        subpath_follows_path = false;
                        break;
                    }
                }
                
                // if subpath followed haplotype path, extend to subsequent subpaths or record completed alignment
                if (subpath_follows_path) {
                    int32_t extended_prefix_score = subpath_prefix_score[i] + subpath.score();
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




