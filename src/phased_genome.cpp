//
//  phased_genome.cpp
//

#include "phased_genome.hpp"

//#define debug_phased_genome

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
    
    PhasedGenome::PhasedGenome(const SnarlManager& snarl_manager) : snarl_manager(&snarl_manager) {
        // nothing to do
    }
    
    PhasedGenome::~PhasedGenome() {
        
        for (Haplotype* haplotype : haplotypes) {
            delete haplotype;
        }
        
    }
    
    PhasedGenome::PhasedGenome(PhasedGenome& rhs){
        
        *this = rhs;
    }
    PhasedGenome& PhasedGenome::operator = (PhasedGenome& phased_genome){

        snarl_manager = phased_genome.snarl_manager;

        for (Haplotype* haplotype : haplotypes) {
            delete haplotype;
        }
        node_locations.clear();
        site_starts.clear();
        site_ends.clear();
        haplotypes.clear();
        
        for(int i = 0; i < phased_genome.haplotypes.size(); i++ ){
            // build haplotypes    
            Haplotype* new_haplo = new Haplotype(phased_genome.begin(i), phased_genome.end(i)); 
            
            haplotypes.push_back(new_haplo);
        }
                     
        // build indices on the new object 
        build_indices();
        
        return phased_genome;
    } 
    
    void PhasedGenome::build_indices() {
        
#ifdef debug_phased_genome
        cerr << "[PhasedGenome::build_indices]: building node id to site index" << endl;
#endif
        // construct the start and end of site indices
        for (const Snarl* snarl : snarl_manager->top_level_snarls()) {
            build_site_indices_internal(snarl);
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
            unordered_map<const Snarl*, HaplotypeNode*> site_start_sides;
            unordered_map<const Snarl*, HaplotypeNode*> site_end_sides;
            
            // iterate along the node path of the haplotype
            HaplotypeNode* haplo_node = haplotype->left_telomere_node;
            while (haplo_node != nullptr) {
                int64_t node_id = haplo_node->node_traversal.node->id();
#ifdef debug_phased_genome
                cerr << "[PhasedGenome::build_indices]: recording an instance of node " << node_id << " in a haplotype node at " << haplo_node << endl;
#endif
                // mark this instance of the node in the node location index
                node_locations[node_id].push_back(haplo_node);
                
                
                // are we at the start of a site?
                if (site_starts.count(node_id)) {
                    const Snarl* site = site_starts[node_id];

                    // are we leaving or entering the site?
                    if (site_start_sides.count(site)) {
#ifdef debug_phased_genome
                        cerr << "[PhasedGenome::build_indices]: leaving at start of site " << site->start().node_id() << "->" << site->end().node_id() << endl;
#endif
                        // leaving: put the site in the index in the orientation of haplotype travesal
                        HaplotypeNode* other_side_node = site_start_sides[site];
                        site_start_sides.erase(site);
                        haplotype->sites[site] = make_pair(other_side_node, haplo_node);
                    }
                    else {
#ifdef debug_phased_genome
                        cerr << "[PhasedGenome::build_indices]: entering at start of site " << site->start().node_id() << "->" << site->end().node_id() << endl;
#endif
                        // entering: mark the node in the haplotype path where we entered
                        site_end_sides[site] = haplo_node;
                    }
                }
                // are we at the start of a site?
                if (site_ends.count(node_id)) {
                    const Snarl* site = site_ends[node_id];
                    // are we leaving or entering the site?
                    if (site_end_sides.count(site)) {
#ifdef debug_phased_genome
                        cerr << "[PhasedGenome::build_indices]: leaving at end of site " << site->start().node_id() << "->" << site->end().node_id() << endl;
#endif
                        // leaving: put the site in the index in the orientation of haplotype travesal
                        HaplotypeNode* other_side_node = site_end_sides[site];
                        site_end_sides.erase(site);
                        haplotype->sites[site] = make_pair(other_side_node, haplo_node);
                    }
                    else {
#ifdef debug_phased_genome
                        cerr << "[PhasedGenome::build_indices]: entering at end of site " << site->start().node_id() << "->" << site->end().node_id() << endl;
#endif
                        // entering: mark the node in the haplotype path where we entered
                        site_start_sides[site] = haplo_node;
                    }
                }
                
                // advance to next node in haplotype path
                haplo_node = haplo_node->next;
            }
        }
    }
    
    void PhasedGenome::build_site_indices_internal(const Snarl* snarl) {
        
#ifdef debug_phased_genome
        cerr << "[PhasedGenome::build_indices]: recording start and end of site " << snarl->start().node_id() << "->" << snarl->end().node_id() << " at " << snarl << endl;
#endif
        
        // record start and end of site
        site_starts[snarl->start().node_id()] = snarl;
        site_ends[snarl->end().node_id()] = snarl;
        
        // recurse through child sites
        for (const Snarl* subsnarl : snarl_manager->children_of(snarl)) {
            build_site_indices_internal(subsnarl);
        }
    }
    
    void PhasedGenome::swap_label(const Snarl& site, Haplotype& haplotype_1, Haplotype& haplotype_2) {
        
        // update index for this site
        if (haplotype_1.sites.count(&site)) {
            if (haplotype_2.sites.count(&site)) {
#ifdef debug_phased_genome
                cerr << "[PhasedGenome::swap_label]: site " << site.start().node_id() << "->" << site.end().node_id() << " is present on both alleles, swapping the label" << endl;
#endif
                // site is in both haplotypes, swap the labels
                swap(haplotype_1.sites[&site], haplotype_2.sites[&site]);
            }
            else {
#ifdef debug_phased_genome
                cerr << "[PhasedGenome::swap_label]: site " << site.start().node_id() << "->" << site.end().node_id() << " is present on only haplotype 1, transferring the label" << endl;
#endif
                // site is only in haplotype 1, transfer the label
                haplotype_2.sites[&site] = haplotype_1.sites[&site];
                haplotype_1.sites.erase(&site);
            }
        }
        else if (haplotype_2.sites.count(&site)) {
#ifdef debug_phased_genome
            cerr << "[PhasedGenome::swap_label]: site " << site.start().node_id() << "->" << site.end().node_id() << " is present on only haplotype 2, transferring the label" << endl;
#endif
            // site is only in haplotype 2, transfer the label
            haplotype_1.sites[&site] = haplotype_2.sites[&site];
            haplotype_2.sites.erase(&site);
        }
        else {
#ifdef debug_phased_genome
            cerr << "[PhasedGenome::swap_label]: site " << site.start().node_id() << "->" << site.end().node_id() << " is present on neither allele, ending recursion" << endl;
#endif
            // site is in neither haplotype
            return;
        }
        
        // update index for child sites
        for (const Snarl* child_site : snarl_manager->children_of(&site)) {
            swap_label(*child_site, haplotype_1, haplotype_2);
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

    vector<id_t> PhasedGenome::get_haplotypes_with_snarl(const Snarl* snarl_to_find){

        // a vector that will hold the haplotype IDs of haplotypes found to traverse through the snarl
        vector<id_t> matched_haplotype_ids;

        // interate through the vector of haplotype pointers and do a lookup for snarl_to_find
        // if found then we add it to the list of matched haplotypes 
        unordered_map<const Snarl*, pair<HaplotypeNode*, HaplotypeNode*> >::iterator it;
        id_t id = 0;
        for (Haplotype* haplotype : haplotypes){
            bool found = haplotype->sites.count(snarl_to_find);

            if(found){
                // add the ID to the haplotype to the vector 
                matched_haplotype_ids.push_back(id);
            }
            id++;
        }
        
        return matched_haplotype_ids;

    }

    void PhasedGenome::print_phased_genome(){
         // output number of haplotypes contained in phased genome 
        size_t haplo_num = num_haplotypes();
        //cerr << "The haplotype num is: " << haplo_num << endl;

        // iterate through the genome and all its haplotypes
        for(int i = 0; i < haplo_num; i++){
            cerr << "Haplotype ID: " << i <<endl;
            // iterate through each node in the haplotype
            for(auto iter = begin(i); iter != end(i); iter++ ){    
                cerr << "node " << (*iter).node->id() << ": " <<(*iter).node->sequence() <<endl;
            }
            cerr<<endl;
            cerr<<endl;
        }     


    }
    vector<NodeTraversal> PhasedGenome::get_allele(const Snarl& site, int which_haplotype) {
        
        Haplotype& haplotype = *haplotypes[which_haplotype];
        
        // can only get the allele of a site that already is in the haplotype
        assert(haplotype.sites.count(&site));
        
        // get the allele
        pair<HaplotypeNode*, HaplotypeNode*> haplo_site = haplotype.sites[&site];
        vector<NodeTraversal> allele;
        for (HaplotypeNode* haplo_node = haplo_site.first->next; haplo_node != haplo_site.second;
             haplo_node = haplo_node->next) {
            allele.push_back(haplo_node->node_traversal);
        }
        
        // is site in the reverse direction on haplotype?
        if (haplo_site.first->node_traversal.node->id() != site.start().node_id()) {
            // reverse the allele
            reverse(allele.begin(), allele.end());
            // swap the orientation
            for (NodeTraversal& node_traversal : allele) {
                node_traversal.backward = !node_traversal.backward;
            }
        }
        return allele;
    }
    
    void PhasedGenome::swap_alleles(const Snarl& site, int haplotype_1, int haplotype_2) {
        Haplotype& haplo_1 = *haplotypes[haplotype_1];
        Haplotype& haplo_2 = *haplotypes[haplotype_2];
        
        pair<HaplotypeNode*, HaplotypeNode*> haplo_nodes_1 = haplo_1.sites[&site];
        pair<HaplotypeNode*, HaplotypeNode*> haplo_nodes_2 = haplo_2.sites[&site];
        
#ifdef debug_phased_genome
        cerr << "[PhasedGenome::swap_alleles]: swapping allele at site " << site.start().node_id() << "->" << site.end().node_id() << " between chromosomes " << haplotype_1 << " and " << haplotype_2 << " with haplotype nodes " << haplo_nodes_1.first << "->" << haplo_nodes_1.second << " and " << haplo_nodes_2.first << "->" << haplo_nodes_2.second << endl;
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
            
            std::swap(haplo_nodes_1.first->next, haplo_nodes_2.first->next);
            std::swap(haplo_nodes_1.second->prev, haplo_nodes_2.second->prev);
        }
        // else two deletions and nothing will change
        
#ifdef debug_phased_genome
        cerr << "[PhasedGenome::swap_alleles]: swapping labels on nested sites" << endl;
#endif
        
        // update index for child sites
        for (const Snarl* child_site : snarl_manager->children_of(&site)) {
            swap_label(*child_site, haplo_1, haplo_2);
        }
    }

    double PhasedGenome::read_log_likelihood(const multipath_alignment_t& multipath_aln, double log_base) {
#ifdef debug_phased_genome
        cerr << "[PhasedGenome::read_log_likelihood] computing read likelihood with log base " << log_base << endl;
        cerr << "read:" << endl;
        cerr << debug_string(multipath_aln) << endl;
        
#endif
        
        if (multipath_aln.mapping_quality() == 0) {
            // this is the answer we'll produce anyway and handling it as an edge case
            // avoids numerical problems at the end
            return 0.0;
        }
        
        // an accumulator that we will sum the log-likelihood of each alignment into
        double log_likelihood = numeric_limits<double>::lowest();
        
        // iteration functions to facilitate iterating on forward/reverse strands
        auto move_right = [](HaplotypeNode*& path_node) { path_node = path_node->next; };
        auto move_left = [](HaplotypeNode*& path_node) { path_node = path_node->prev; };
        
        /*
         * STEP 1: find which subpaths are represented on the current haplotypes and which
         * subpaths are adjacent to each order (expressed as "links")
         */
        
        // records of: (subpath coming from, which match of that subpath, node the next should be on, should it be matching the orientation?)
        vector<vector<tuple<size_t, size_t, HaplotypeNode*, bool>>> possible_forward_links(multipath_aln.subpath_size());
        
        // for each subpath
        //     for each copy of it in the phased genome
        //         the (subpath index, which copy)'s that this copy links back to
        vector<vector<vector<pair<size_t, size_t>>>> backward_links(multipath_aln.subpath_size());
        
        for (size_t i = 0; i < multipath_aln.subpath_size(); ++i) {
#ifdef debug_phased_genome
            cerr << "[PhasedGenome::read_log_likelihood] looking for matches of subpath " << i << endl;
            
#endif
            
            const subpath_t& subpath = multipath_aln.subpath(i);
            const path_t& path = subpath.path();

            // check all of the locations of this subpath among the haplotypes
            for (HaplotypeNode* starting_haplo_node : node_locations[path.mapping(0).position().node_id()]) {

                // are we traversing forward or backward along the haplotype?
                bool matches_orientation = (starting_haplo_node->node_traversal.backward == path.mapping(0).position().is_reverse());
                auto move_forward = matches_orientation ? move_right : move_left;
       
#ifdef debug_phased_genome
                cerr << "[PhasedGenome::read_log_likelihood] found starting haplo node that has " << (matches_orientation ? "matching" : "reverse") << " orientation at " << starting_haplo_node << " with trav " << starting_haplo_node->node_traversal << endl;
#endif
                
                // determine whether the rest of the subpath matches
                bool full_match = true;
                HaplotypeNode* haplo_node = starting_haplo_node;
                for (size_t j = 1; j < path.mapping_size(); ++j) {
                    
                    // advance to the haplo node that we would find next
                    move_forward(haplo_node);
                    
#ifdef debug_phased_genome
                    cerr << "[PhasedGenome::read_log_likelihood] moving forward to haplo node at " << haplo_node << " with trav " << haplo_node->node_traversal << endl;
#endif
                    
                    const position_t& pos = path.mapping(j).position();
                    if (haplo_node->node_traversal.node->id() != pos.node_id() ||
                        (haplo_node->node_traversal.backward == pos.is_reverse()) != matches_orientation) {
                        // the subpath doesn't match the path of the haplotype here
                        full_match = false;
                        break;
                    }
                }
                
                if (!full_match) {
                    // scores aren't necessarily dynamic programmable except on full matches, so
                    // we'll leave this alone
                    // TODO: allow for alignments to partial subpaths?
                    continue;
                }
                
#ifdef debug_phased_genome
                cerr << "[PhasedGenome::read_log_likelihood] found a full match" << endl;
#endif
                
                const path_mapping_t& final_mapping = path.mapping(path.mapping_size() - 1);
                size_t final_offset = mapping_from_length(final_mapping) + final_mapping.position().offset();
                if (final_offset == haplo_node->node_traversal.node->sequence().size()) {
                    // the last mapping hits the end of its node, so we expect to find
                    // the next subpath on the following node
                    move_forward(haplo_node);
#ifdef debug_phased_genome
                    cerr << "[PhasedGenome::read_log_likelihood] any subsequent match is expected to be on the next haplo node at " << haplo_node << endl;
#endif
                }
                
                // we found a full match of this subpath on the phased genome, add a match
                size_t match_num = backward_links[i].size();
                backward_links[i].emplace_back();
                vector<pair<size_t, size_t>>& links = backward_links[i].back();
                
                // let's check if we could have extended along any of the forward links we
                // previously identified
                for (auto& possible_link : possible_forward_links[i]) {
                    // TODO: make this a multimap from (node,orientation) instead?
                    if (get<2>(possible_link) == starting_haplo_node &&
                        get<3>(possible_link) == matches_orientation) {
                        // we started on haplo node and orientation that we would have expected
                        // if we followed this forward link
                        links.emplace_back(get<0>(possible_link), get<1>(possible_link));
#ifdef debug_phased_genome
                        cerr << "[PhasedGenome::read_log_likelihood] confirmed a link from subpath " << get<0>(possible_link) << ", match number " << get<1>(possible_link) << endl;
#endif
                    }
                }
                
                // add a candidate forward link from here to each subsequent subpath
                for (size_t j : subpath.next()) {
#ifdef debug_phased_genome
                    cerr << "[PhasedGenome::read_log_likelihood] adding a possible link to " << j << endl;
#endif
                    possible_forward_links[j].emplace_back(i, match_num, haplo_node, matches_orientation);
                }
            }
        }
        
        /*
         * STEP 2: follow the links we discovered in step 1 to compute the
         * scores of the alignments consistent with the current haplotypes
         * and add them to the log likelihood
         */
        
        // to keep track of which partial alignments we've already accounted for
        vector<vector<bool>> traversed(multipath_aln.subpath_size());
        for (size_t i = 0; i < traversed.size(); ++i) {
            traversed[i].resize(backward_links[i].size(), false);
        }
        
        // iterate backwards over the backward facing links for each subpath
        for (int64_t i = backward_links.size() - 1; i >= 0; --i) {
            auto& links = backward_links[i];
            // iterate over the links for each match we found for this subpath
            for (size_t j = 0; j < links.size(); ++j) {
                if (traversed[i][j]) {
                    // we already scored the alignments corresponding to these
                    // links
                    continue;
                }
#ifdef debug_phased_genome
                cerr << "[PhasedGenome::read_log_likelihood] checking backward links from " << i << ", match num " << j << endl;
#endif
                
                // use DFS to generate all longest-possible alignments along
                // the backward links
                
                // records of (subpath idx, copy of subpath, index of next link to take)
                vector<tuple<size_t, size_t, size_t>> stack;
                stack.emplace_back(i, j, 0);
                while (!stack.empty()) {
                    auto& record = stack.back();
#ifdef debug_phased_genome
                    cerr << "[PhasedGenome::read_log_likelihood] unstacking " << get<0>(record) << ", " << get<1>(record) << ", " << get<2>(record) << endl;
#endif
                    
                    auto& links = backward_links[get<0>(record)][get<1>(record)];
                    if (get<2>(record) < links.size()) {
                        auto next = links[get<2>(record)++];
                        stack.emplace_back(next.first, next.second, 0);
#ifdef debug_phased_genome
                        cerr << "[PhasedGenome::read_log_likelihood] following link to " << next.first << ", " << next.second << endl;
#endif
                    }
                    else {
                        // we've finished traversing this
                        traversed[get<0>(record)][get<1>(record)] = true;
#ifdef debug_phased_genome
                        cerr << "[PhasedGenome::read_log_likelihood] finished traversing " << get<0>(record) << ", " << get<1>(record) << endl;
#endif
                        if (links.empty()) {
                            // this is the final subpath in the alignment, so the stack now
                            // represents a valid alignment. we can compute the highest scoring
                            // segment of the alignment with dynamic programming
#ifdef debug_phased_genome
                            cerr << "[PhasedGenome::read_log_likelihood] completed an alignment, current alignment is " << endl;
                            for (auto& stack_record : stack) {
                                cerr << "\t" << get<0>(stack_record) << ", " << get<1>(stack_record) << ", " << get<2>(stack_record) << endl;
                            }
#endif
                            
                            int32_t current_score = multipath_aln.subpath(get<0>(stack[0])).score();
                            int32_t max_score = current_score;
                            for (size_t k = 1; k < stack.size(); ++k) {
                                int32_t subpath_score = multipath_aln.subpath(get<0>(stack[k])).score();
                                current_score = max(current_score + subpath_score, subpath_score);
                                max_score = max(max_score, current_score);
                            }
                            
                            // TODO: do i need to check if the same alignment sub sequence gets
                            // double-counted here?
                            log_likelihood = add_log(log_likelihood, max_score * log_base);
#ifdef debug_phased_genome
                            cerr << "[PhasedGenome::read_log_likelihood] added max score " << max_score << " to update likelihood to " << log_likelihood << endl;
#endif
                        }
                        stack.pop_back();
                    }
                }
            }
        }
        
        // adjust by the mapping quality (likelihood = (1-err)*l + err)
        double mapping_err_log_prob = phred_to_logprob(multipath_aln.mapping_quality());
        return add_log(subtract_log(0.0, mapping_err_log_prob) + log_likelihood,
                       mapping_err_log_prob);
    }
    
    int32_t PhasedGenome::optimal_score_on_genome(const multipath_alignment_t& multipath_aln, VG& graph) {
        
        
        // must have identified start subpaths before computing optimal score   
        assert(multipath_aln.start_size() > 0);
        
        // iteration functions to facilitate iterating on forward/reverse strands
        auto move_right = [](HaplotypeNode*& path_node) { path_node = path_node->next; };
        auto move_left = [](HaplotypeNode*& path_node) { path_node = path_node->prev; };
        
        int32_t optimal_score = 0;
        
        // find the places in the path where the alignment might start
        unordered_map< pair<HaplotypeNode*, bool>, vector<int>> candidate_start_positions;
        for (int i = 0; i < multipath_aln.start_size(); i++) {
            // a starting subpath in the multipath alignment
            const subpath_t& start_subpath = multipath_aln.subpath(multipath_aln.start(i));
            const position_t& start_pos = start_subpath.path().mapping(0).position();
            
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
                
                const subpath_t& subpath = multipath_aln.subpath(i);
                
                // iterate through mappings in this subpath (assumes one mapping per node)
                bool subpath_follows_path = true;
                for (int j = 0; j < subpath.path().mapping_size(); j++) {
                    // check if mapping corresponds to the next node in the path in the correct orientation
                    const path_mapping_t& mapping = subpath.path().mapping(j);
                    const position_t& position = mapping.position();
                    if (position.node_id() != subpath_node->node_traversal.node->id()
                        || ((position.is_reverse() == subpath_node->node_traversal.backward) != oriented_forward)) {
                        subpath_follows_path = false;
                        
#ifdef debug_phased_genome
                        cerr << "[PhasedGenome::optimal_score_on_genome]: subpath " << i << " is inconsistent with haplotype" << endl;
                        
#endif              
                        break;
                    }
                    
                    if (position.offset() + mapping_from_length(mapping) == subpath_node->node_traversal.node->sequence().size()) {
                        move_forward(subpath_node);
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




