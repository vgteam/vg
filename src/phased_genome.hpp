//
//  phased_genome.hpp
//  
// Defines a object that represents a collection of haplotypes as walks through a variation graph.
// The object is designed to be dynamic and support fast editing operations, such as changing the
// allele at a site and swapping alleles between chromosomes.
//

#ifndef phased_genome_hpp
#define phased_genome_hpp

#include <stdio.h>
#include <cassert>
#include <list>
#include "vg.pb.h"
#include "vg.hpp"
#include "nodetraversal.hpp"
#include "genotypekit.hpp"

//#define debug_phased_genome

using namespace std;

namespace vg {
    
    // a collection of haplotypes that represent all of the chromosomes of a genome (including
    // phasing) as walks through a variation graph, also maintains indices for sites and nodes
    class PhasedGenome {
        
        // TODO: Do I also need to keep track of the ploidy of sites in case the path includes them
        // multiple times?
        // I think for now it will suffice to keep them in an arbitrary order as long as they are accessible by it
        // pay attention to logic for relabeling though, that might get hairy
        
    public:
        
        // unidirectional iterator along a haplotype
        class iterator;
        
        /*
         *  CONSTRUCTION METHODS
         */
        
        // constructor
        PhasedGenome(VG& graph);
        ~PhasedGenome();
        
        // build a haplotype in place from an iterator that returns NodeTraversal objects from its
        // dereference operator (allows construction without instantiating the haplotype elsewhere)
        // returns the numerical id of the new haplotype
        //
        // note: the haplotype must have at least one node
        template <typename NodeTraversalIterator>
        int add_haplotype(NodeTraversalIterator first, NodeTraversalIterator last);
        
        // TODO: make an interface for constructing a set of homologous haplotypes at the same time
        // from a single iterator so that they don't need to be instantiated anywhere else
        
        // construct the site ends, node locations, and haplotype site location indices
        // intended to be called one time after building haplotypes, after which the indices
        // are maintained automatically during edit operations
        void build_indices(vector<NestedSite>& top_level_sites);
        
        /*
         *  ITERATION AND UTILITY METHODS
         */
        
        size_t num_haplotypes();
        
        // unidirectional iterators beginnning at the left telomere and moving to the right
        iterator begin(int which_haplotype);
        iterator end(int which_haplotype);
        
        /*
         *  HAPLOTYPE EDITING METHODS
         */
        
        // swap the allele from one haplotype to the other, maintaining all indices
        // if a nested site is being swapped, only the topmost site should be given
        void swap_alleles(NestedSite& site, int haplotype_1, int haplotype_2);
        
        // set the allele at a site with an iterator that yields its node sequence (not including the
        // source or sink nodes of the site, e.g. an empty iterator for a deletion)
        // allele should be provided in the order indicated by the NestedSite (i.e. from start to end)
        //
        // note: function does not check that the allele path takes only edges that are actually
        // included in the graph, so client must ensure this itself
        template <typename NodeTraversalIterator>
        void set_allele(NestedSite& site, NodeTraversalIterator first, NodeTraversalIterator last,
                        int which_haplotype);
        
        // returns the score of the highest scoring alignment contained in the multipath alignment
        // that is restricted to the phased genome's paths through the variation graph
        //
        // note: assumes that MultipathAlignment has 'start' field filled in
        int32_t optimal_score_on_genome(const MultipathAlignment& multipath_aln);
        
        // TODO: make a local subalignment optimal score function (main obstacle is scoring partial subpaths)
        
    private:
        
        struct HaplotypeNode;
        class Haplotype;
        
        // graph though which the genome is threaded
        VG& graph;
        
        // all haplotypes in the genome (generally 2 per chromosome)
        vector<Haplotype*> haplotypes;
        
        // an index of where nodes from the graph occur in the graph
        unordered_map<int64_t, list<HaplotypeNode*> > node_locations;
        
        // sufficient for these purposes to maintain only node ids instead of node sides
        // since the path must go through the site either before or after entering here
        unordered_map<int64_t, NestedSite*> site_starts;
        unordered_map<int64_t, NestedSite*> site_ends;
        
        // helper function
        void build_site_indices_internal(NestedSite* site);
        
        // no safety checks that the adjacent nodes aren't null (should only be used in interior of haplotype)
        // these operations maintain the node location indices but no others
        inline void insert_left(NodeTraversal node_traversal, HaplotypeNode* haplo_node);
        inline void insert_right(NodeTraversal node_traversal, HaplotypeNode* haplo_node);
        inline void remove(HaplotypeNode* haplo_node);
        
        // swap a site location in indices after swapping an allele
        void swap_label(NestedSite& site, Haplotype& haplotype_1, Haplotype& haplotype_2);
        
    };
    
    /*
     *  INTERNAL OBJECTS (NO PUBLIC-FACING API BELOW HERE)
     */
    
    // a node in walk through the graph taken by a haplotype
    struct PhasedGenome::HaplotypeNode {
    public:
        
        NodeTraversal node_traversal;
        HaplotypeNode* next;
        HaplotypeNode* prev;
        
        HaplotypeNode(NodeTraversal node_traversal, HaplotypeNode* next, HaplotypeNode* prev);
        ~HaplotypeNode();
        
    };
    
    // hash function for node traversals
    struct HashNodeTraversal {
        size_t operator()(const NodeTraversal& node_traversal) const {
            return (size_t) 1099511628211ull * ((uintptr_t) node_traversal.node + node_traversal.backward * 16777619ull) + 14695981039346656037ull;
        }
    };
    
    // hash function for sites
    struct HashSite {
        size_t operator()(const NestedSite& site) const {
            HashNodeTraversal hsh;
            return hsh(site.start) ^ (hsh(site.end) << 1);
        }
    };
    
    // a specialized linked list that tracks a walk through the variation graph and maintains
    // an index of sites
    class PhasedGenome::Haplotype {
        
        PhasedGenome::HaplotypeNode* left_telomere_node;
        PhasedGenome::HaplotypeNode* right_telomere_node;
        
        // assumes pairs are given such that .first->prev and .second->next are outside the bubble
        unordered_map<NestedSite, pair<HaplotypeNode*, HaplotypeNode*>, HashSite> sites;
        
    public:
        // construct a haplotype with a single node
        Haplotype(NodeTraversal node_traversal);
        
        // construct a haplotype with a list of nodes
        template <typename NodeTraversalIterator>
        Haplotype(NodeTraversalIterator first, NodeTraversalIterator last);
        
        ~Haplotype();
        
        // intended for use by the PhasedGenome to build haploytpes initially, returns new node
        inline HaplotypeNode* append_left(NodeTraversal node_traversal);
        inline HaplotypeNode* append_right(NodeTraversal node_traversal);
        
        friend class PhasedGenome;
        friend class HaplotypeNode;
        
    };
    
    // iterator to obtain the NodeTraversals of a haplotype
    class PhasedGenome::iterator {
    private:
        
        size_t rank;
        int haplotype_number;
        HaplotypeNode* haplo_node;
        
        iterator(size_t rank, int haplotype_number, HaplotypeNode* haplo_node);
        
    public:
        
        iterator();
        iterator(const iterator& other);
        ~iterator();
        
        // UNIDIRECTIONAL ITERATOR INTERFACE
        
        inline iterator& operator=(const iterator& other) {
            rank = other.rank;
            haplotype_number = other.haplotype_number;
            haplo_node = other.haplo_node;
            return *this;
        }
        
        inline bool operator==(const iterator& other) const {
            return rank == other.rank && haplotype_number == other.haplotype_number;
        }
        
        inline bool operator!=(const iterator& other) const {
            return rank != other.rank || haplotype_number != other.haplotype_number;
        }
        
        inline iterator operator++() {
            haplo_node = haplo_node->next;
            rank = (haplo_node == nullptr) ? 0 : rank + 1;
            return *this;
        }
        
        inline iterator operator++( int ) {
            iterator temp = *this;
            haplo_node = haplo_node->next;
            rank = (haplo_node == nullptr) ? 0 : rank + 1;
            return temp;
        }
        
        inline NodeTraversal operator*(){
            return haplo_node->node_traversal;
        }
        
        inline int which_haplotype() {
            return haplotype_number;
        }
        
        friend class PhasedGenome;
        friend class Haplotype;
    };
    
    /*
     *  TEMPLATE FUNCTIONS
     */
    
    template <typename NodeTraversalIterator>
    int PhasedGenome::add_haplotype(NodeTraversalIterator first, NodeTraversalIterator last) {
        
#ifdef debug_phased_genome
        cerr << "[PhasedGenome::add_haplotype]: adding haplotype number " << haplotypes.size() << endl;
#endif
        
        Haplotype* haplotype = new Haplotype(first, last);
        haplotypes.push_back(haplotype);
        
        return haplotypes.size() - 1;
    }
    
    // TODO: it seems like there should be a better way to do this than completely erasing
    // and then rewriting the site (especially at long sites)
    template <typename NodeTraversalIterator>
    void PhasedGenome::set_allele(NestedSite& site, NodeTraversalIterator first, NodeTraversalIterator last,
                                  int which_haplotype) {
        
        Haplotype& haplotype = *haplotypes[which_haplotype];
        
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
        
        // keeps track of the haplotype node where we entered a site and thereby also
        // indicates whether we have entered the site yet
        unordered_map<NestedSite, HaplotypeNode*, HashSite> subsite_start_side;
        unordered_map<NestedSite, HaplotypeNode*, HashSite> subsite_end_side;
        
        // start inserting nodes from the first position in the allele
        for (; first != last; first++) {
            
            // create the next node in the allele and move the pointer onto the new node
            if (forward) {
                insert_right(*first, haplo_node);
                haplo_node = haplo_node->next;
            }
            else {
                insert_left(*first, haplo_node);
                haplo_node = haplo_node->prev;
            }
            
            int64_t node_id = haplo_node->node_traversal.node->id();
            
            // does a site start here?
            if (site_starts.count(node_id)) {
                NestedSite* subsite = site_starts[node_id];
                // are we entering or leaving this site?
                if (subsite_end_side.count(*subsite)) {
                    // add this site into the haplotype's site index
                    HaplotypeNode* other_side_node = subsite_end_side[*subsite];
                    // the site sides should always be entered in the order that
                    // they occur in the haplotype
                    haplotype.sites[*subsite] = forward ? make_pair(other_side_node, haplo_node)
                                                        : make_pair(haplo_node, other_side_node);
                }
                else {
                    // we are entering a site, mark the location of the entrance
                    subsite_start_side[*subsite] = haplo_node;
                }
            }
            
            // does a site end here?
            if (site_ends.count(node_id)) {
                NestedSite* subsite = site_ends[node_id];
                // are we entering or leaving this site?
                if (subsite_start_side.count(*subsite)) {
                    // add this site into the haplotype's site index
                    HaplotypeNode* other_side_node = subsite_start_side[*subsite];
                    // the site sides should always be entered in the order that
                    // they occur in the haplotype
                    haplotype.sites[*subsite] = forward ? make_pair(other_side_node, haplo_node)
                                                        : make_pair(haplo_node, other_side_node);
                }
                else {
                    // we are entering a site, mark the location of the entrance
                    subsite_end_side[*subsite] = haplo_node;
                }
            }
        }
    }
    
    template <typename NodeTraversalIterator>
    PhasedGenome::Haplotype::Haplotype(NodeTraversalIterator first, NodeTraversalIterator last) {
        
#ifdef debug_phased_genome
        cerr << "[Haplotype::Haplotype]: constructing Haplotype at " << this << endl;
#endif
        
        if (first == last) {
            cerr << "error:[PhasedGenome] cannot construct haplotype with 0 nodes" << endl;
            assert(0);
        }
        
        // construct seed node
        left_telomere_node = new HaplotypeNode(*first, nullptr, nullptr);
        right_telomere_node = left_telomere_node;
        first++;
        
#ifdef debug_phased_genome
        cerr << "[Haplotype::Haplotype]: initialized haplotype with node " << right_telomere_node->node_traversal.node->id() << " with sequence " << right_telomere_node->node_traversal.node->sequence() << " in " << (right_telomere_node->node_traversal.backward ? "reverse" : "forward") << " orientation at memory location " << right_telomere_node << endl;
#endif
        
        // add each subsequent node
        for (; first != last; first++) {
            append_right(*first);
        }
    }
    
    /*
     *   INLINE FUNCTIONS
     */
    
    inline PhasedGenome::HaplotypeNode* PhasedGenome::Haplotype::append_left(NodeTraversal node_traversal) {
        left_telomere_node = new HaplotypeNode(node_traversal, left_telomere_node, nullptr);
        left_telomere_node->next->prev = left_telomere_node;
        return left_telomere_node;
    }
    
    inline PhasedGenome::HaplotypeNode* PhasedGenome::Haplotype::append_right(NodeTraversal node_traversal) {
        right_telomere_node = new HaplotypeNode(node_traversal, nullptr, right_telomere_node);
        right_telomere_node->prev->next = right_telomere_node;
        
#ifdef debug_phased_genome
        cerr << "[Haplotype::append_right]: appended to right side node " << node_traversal.node->id() << " with sequence " << node_traversal.node->sequence() << " in " << (node_traversal.backward ? "reverse" : "forward") << " orientation at memory location " << right_telomere_node << endl;
#endif
        
        return right_telomere_node;
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
}



#endif /* phased_genome_hpp */
