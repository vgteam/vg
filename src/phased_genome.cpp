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
    
    PhasedGenome::Haplotype::~Haplotype() {
        // destruct nodes
        HaplotypeNode* haplo_node = this->left_telomere_node;
        while (haplo_node) {
            HaplotypeNode* next = haplo_node->next;
            delete haplo_node;
            haplo_node = next;
        }
    }
    
    PhasedGenome::PhasedGenome(VG& graph) : graph(graph) {
        // nothing to do
    }
    
    PhasedGenome::~PhasedGenome() {
        
        for (Haplotype* haplotype : haplotypes) {
            delete haplotype;
        }
        
    }
    
    void PhasedGenome::build_indices(vector<NestedSite>& top_level_sites) {
        
#ifdef debug_phased_genome
        cerr << "[PhasedGenome::build_indices]: building node id to site index" << endl;
#endif
        // construct the start and end of site indices
        for (NestedSite& site : top_level_sites) {
            build_site_indices_internal(&site);
        }
        
#ifdef debug_phased_genome
        cerr << "[PhasedGenome::build_indices]: building site to haplotype node index" << endl;
#endif
        
        // label the sites on each haplotype
        for (Haplotype* haplotype : haplotypes) {
#ifdef debug_phased_genome
            cerr << "[PhasedGenome::build_indices]: traversing haplotype at " << haplotype << endl;
#endif
            // keep track of where we enter and leave sites
            unordered_map<NestedSite, HaplotypeNode*, HashSite> site_start_sides;
            unordered_map<NestedSite, HaplotypeNode*, HashSite> site_end_sides;
            
            // iterate along the node path of the haplotype
            HaplotypeNode* haplo_node = haplotype->left_telomere_node;
            while (haplo_node != nullptr) {
                int64_t node_id = haplo_node->node_traversal.node->id();
                
                // mark this instance of the node in the node location index
                node_locations[node_id].push_back(haplo_node);
                
                
                // are we at the start of a site?
                NestedSite* site = nullptr;
                if (site_starts.count(node_id)) {
                    site = site_starts[node_id];

                    // are we leaving or entering the site?
                    if (site_start_sides.count(*site)) {
#ifdef debug_phased_genome
                        cerr << "[PhasedGenome::build_indices]: leaving at start of site " << site->start.node->id() << "->" << site->end.node->id() << endl;
#endif
                        // leaving: put the site in the index in the orientation of haplotype travesal
                        HaplotypeNode* other_side_node = site_start_sides[*site];
                        site_start_sides.erase(*site);
                        haplotype->sites[*site] = make_pair(other_side_node, haplo_node);
                    }
                    else {
#ifdef debug_phased_genome
                        cerr << "[PhasedGenome::build_indices]: entering at start of site " << site->start.node->id() << "->" << site->end.node->id() << endl;
#endif
                        // entering: mark the node in the haplotype path where we entered
                        site_end_sides[*site] = haplo_node;
                    }
                }
                // are we at the start of a site?
                if (site_ends.count(node_id)) {
                    site = site_ends[node_id];
                    // are we leaving or entering the site?
                    if (site_end_sides.count(*site)) {
#ifdef debug_phased_genome
                        cerr << "[PhasedGenome::build_indices]: leaving at end of site " << site->start.node->id() << "->" << site->end.node->id() << endl;
#endif
                        // leaving: put the site in the index in the orientation of haplotype travesal
                        HaplotypeNode* other_side_node = site_end_sides[*site];
                        site_end_sides.erase(*site);
                        haplotype->sites[*site] = make_pair(other_side_node, haplo_node);
                    }
                    else {
#ifdef debug_phased_genome
                        cerr << "[PhasedGenome::build_indices]: entering at end of site " << site->start.node->id() << "->" << site->end.node->id() << endl;
#endif
                        // entering: mark the node in the haplotype path where we entered
                        site_start_sides[*site] = haplo_node;
                    }
                }
                
                // advance to next node in haplotype path
                haplo_node = haplo_node->next;
            }
        }
    }
    
    void PhasedGenome::build_site_indices_internal(NestedSite* site) {
        
#ifdef debug_phased_genome
        cerr << "[PhasedGenome::build_indices]: recording start and end of site " << site->start.node->id() << "->" << site->end.node->id() << " at " << site << endl;
#endif
        
        // record start and end of site
        site_starts[site->start.node->id()] = site;
        site_ends[site->end.node->id()] = site;
        
        // recurse through child sites
        for (NestedSite& subsite : site->children) {
            build_site_indices_internal(&subsite);
        }
    }
    
    void PhasedGenome::swap_label(NestedSite& site, Haplotype& haplotype_1, Haplotype& haplotype_2) {
        
        // update index for this site
        if (haplotype_1.sites.count(site)) {
            if (haplotype_2.sites.count(site)) {
#ifdef debug_phased_genome
                cerr << "[PhasedGenome::swap_label]: site " << site.start.node->id() << "->" << site.end.node->id() << " is present on both alleles, swapping the label" << endl;
#endif
                // site is in both haplotypes, swap the labels
                swap(haplotype_1.sites[site], haplotype_2.sites[site]);
            }
            else {
#ifdef debug_phased_genome
                cerr << "[PhasedGenome::swap_label]: site " << site.start.node->id() << "->" << site.end.node->id() << " is present on only haplotype 1, transferring the label" << endl;
#endif
                // site is only in haplotype 1, transfer the label
                haplotype_2.sites[site] = haplotype_1.sites[site];
                haplotype_1.sites.erase(site);
            }
        }
        else if (haplotype_2.sites.count(site)) {
#ifdef debug_phased_genome
            cerr << "[PhasedGenome::swap_label]: site " << site.start.node->id() << "->" << site.end.node->id() << " is present on only haplotype 2, transferring the label" << endl;
#endif
            // site is only in haplotype 2, transfer the label
            haplotype_1.sites[site] = haplotype_2.sites[site];
            haplotype_2.sites.erase(site);
        }
        else {
#ifdef debug_phased_genome
            cerr << "[PhasedGenome::swap_label]: site " << site.start.node->id() << "->" << site.end.node->id() << " is present on neither allele, ending recursion" << endl;
#endif
            // site is in neither haplotype
            return;
        }
        
        // update index for child sites
        for (NestedSite& child_site : site.children) {
            swap_label(child_site, haplotype_1, haplotype_2);
        }
    }
    
    size_t PhasedGenome::num_haplotypes() {
        return haplotypes.size();
    }
    
    PhasedGenome::iterator PhasedGenome::begin(int which_haplotype) {
        return iterator(1, which_haplotype, haplotypes[which_haplotype]->left_telomere_node);
    }
    
    PhasedGenome::iterator PhasedGenome::end(int which_haplotype) {
        return iterator(0, which_haplotype, nullptr);
    }
    
    void PhasedGenome::swap_alleles(NestedSite& site, int haplotype_1, int haplotype_2) {
        Haplotype& haplo_1 = *haplotypes[haplotype_1];
        Haplotype& haplo_2 = *haplotypes[haplotype_2];
        
        pair<HaplotypeNode*, HaplotypeNode*> haplo_nodes_1 = haplo_1.sites[site];
        pair<HaplotypeNode*, HaplotypeNode*> haplo_nodes_2 = haplo_2.sites[site];
        
#ifdef debug_phased_genome
        cerr << "[PhasedGenome::swap_alleles]: swapping allele at site " << site.start.node->id() << "->" << site.end.node->id() << " between chromosomes " << haplotype_1 << " and " << haplotype_2 << " with haplotype nodes " << haplo_nodes_1.first << "->" << haplo_nodes_1.second << " and " << haplo_nodes_2.first << "->" << haplo_nodes_2.second << endl;
#endif
        
        bool is_deletion_1 = (haplo_nodes_1.first->next == haplo_nodes_1.second);
        bool is_deletion_2 = (haplo_nodes_2.first->next == haplo_nodes_2.second);
        
        if (is_deletion_1 && !is_deletion_2) {
            haplo_nodes_2.first->next->prev = haplo_nodes_1.first;
            haplo_nodes_2.second->prev->next = haplo_nodes_1.second;
            
            haplo_nodes_1.first->next = haplo_nodes_2.first->next;
            haplo_nodes_1.second->prev = haplo_nodes_2.second->prev;
            
            haplo_nodes_2.first->next = haplo_nodes_2.second;
            haplo_nodes_2.second->prev = haplo_nodes_2.first;
        }
        else if (is_deletion_2 && !is_deletion_1) {
            haplo_nodes_1.first->next->prev = haplo_nodes_2.first;
            haplo_nodes_1.second->prev->next = haplo_nodes_2.second;
            
            haplo_nodes_2.first->next = haplo_nodes_1.first->next;
            haplo_nodes_2.second->prev = haplo_nodes_1.second->prev;
            
            haplo_nodes_1.first->next = haplo_nodes_1.second;
            haplo_nodes_1.second->prev = haplo_nodes_1.first;
        }
        else if (!is_deletion_1 && !is_deletion_2) {
            haplo_nodes_2.first->next->prev = haplo_nodes_1.first;
            haplo_nodes_2.second->prev->next = haplo_nodes_1.second;
            
            haplo_nodes_1.first->next->prev = haplo_nodes_2.first;
            haplo_nodes_1.second->prev->next = haplo_nodes_2.second;
            
            swap(haplo_nodes_1.first->next, haplo_nodes_2.first->next);
            swap(haplo_nodes_1.second->prev, haplo_nodes_2.second->prev);
        }
        // else two deletions and nothing will change
        
#ifdef debug_phased_genome
        cerr << "[PhasedGenome::swap_alleles]: swapping labels on nested sites" << endl;
#endif
        
        // update index for child sites
        for (NestedSite& child_site : site.children) {
            swap_label(child_site, haplo_1, haplo_2);
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
            
#ifdef debug_phased_genome
            cerr << "[PhasedGenome::optimal_score_on_genome]: looking for candidate start positions for subpath " << multipath_aln.start(i) << " on node " << start_pos.node_id() << endl;
#endif
            
            // add each location the start nodes occur in the path to the candidate starts
            for ( HaplotypeNode* haplo_node : node_locations[start_pos.node_id()] ) {
#ifdef debug_phased_genome
                cerr << "[PhasedGenome::optimal_score_on_genome]: marking candidate start position at " << haplo_node->node_traversal.node->id() << " on haplotype node at " << haplo_node << endl;
#endif
                // mark the start locations orientation relative to the start node
                candidate_start_positions[make_pair(haplo_node, haplo_node->node_traversal.backward == start_pos.is_reverse())].push_back(i);
            }
        }
        
        // check alignments starting at each node in the path that has a source subpath starting on it
        for (pair< pair<HaplotypeNode*, bool>, vector<int> > path_starts : candidate_start_positions) {
            
#ifdef debug_phased_genome
            cerr << "[PhasedGenome::optimal_score_on_genome]: checking for an alignment at candidate start position on node " << path_starts.first.first->node_traversal.node->id() << " on " << (path_starts.first.second ? "forward" : "reverse" ) << " strand of haplotype" << endl;
#endif
            
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
#ifdef debug_phased_genome
                    cerr << "[PhasedGenome::optimal_score_on_genome]: subpath " << i << " is unreachable through consistent paths" << endl;
#endif
                    continue;
                }
                
#ifdef debug_phased_genome
                cerr << "[PhasedGenome::optimal_score_on_genome]: checking subpath " << i << " for consistent paths" << endl;
#endif
                
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
#ifdef debug_phased_genome
                        cerr << "[PhasedGenome::optimal_score_on_genome]: subpath " << i << " is inconsistent with haplotype" << endl;
#endif
                    }
                }
                
                // if subpath followed haplotype path, extend to subsequent subpaths or record completed alignment
                if (subpath_follows_path) {
#ifdef debug_phased_genome
                    cerr << "[PhasedGenome::optimal_score_on_genome]: subpath " << i << " is consistent with haplotype" << endl;
#endif
                    int32_t extended_prefix_score = subpath_prefix_score[i] + subpath.score();
                    if (subpath.next_size() == 0) {
#ifdef debug_phased_genome
                        cerr << "[PhasedGenome::optimal_score_on_genome]: sink path, considering score of " << extended_prefix_score << endl;
#endif
                        // reached a sink subpath (thereby completing an alignment), check for optimality
                        if (extended_prefix_score > optimal_score) {
                            optimal_score = extended_prefix_score;
                        }
                    }
                    else {
                        
#ifdef debug_phased_genome
                        cerr << "[PhasedGenome::optimal_score_on_genome]: non sink path, extending score of " << extended_prefix_score << endl;
#endif
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
    
    
    PhasedGenome::iterator::iterator() : rank(0), haplotype_number(-1), haplo_node(nullptr) {
    
    }
    
    PhasedGenome::iterator::iterator(size_t rank, int haplotype_number, HaplotypeNode* haplo_node) :
                                    rank(rank), haplotype_number(haplotype_number), haplo_node(haplo_node) {
        
    }
    
    PhasedGenome::iterator::iterator(const iterator& other) : rank(other.rank),
                                                              haplotype_number(other.haplotype_number),
                                                              haplo_node(other.haplo_node) {
        
    }
    
    PhasedGenome::iterator::~iterator() {
    
    }
    
}




