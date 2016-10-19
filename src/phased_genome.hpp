//
//  phased_genome.hpp
//  
// defines a object that represents a collection of haplotypes as walks through a variation graph
// the object is designed to be dynamic and support fast editing operations, such as changing the
// allele at a site and swapping alleles between chromosomes
//

#ifndef phased_genome_hpp
#define phased_genome_hpp

#include <stdio.h>
#include <cassert>
#include <list>
#include "vg.pb.h"
#include "vg.hpp"
#include "genotypekit.hpp"

using namespace std;

namespace vg {
    
    // a collection of haplotypes that represent all of the chromosomes of a genome (including
    // phasing) as walks through a variation graph, also maintains indices for sites and nodes
    class PhasedGenome {
        
    public:
        
        // constructor
        PhasedGenome(VG& graph);
        ~PhasedGenome();
        
        // TODO: Do I also need to keep track of the ploidy of sites in case the path includes them
        // multiple times?
        // I think for now it will suffice to keep them in an arbitrary order as long as they are accessible by it
        // pay attention to logic for relabeling though, that might get hairy
        
        
        // build a haplotype in place from an iterator that returns NodeTraversal objects from its
        // dereference operator
        //
        // note: the haplotype must have at least one node
        template <typename NodeTraversalIterator>
        void add_haplotype(NodeTraversalIterator first, NodeTraversalIterator last);
        
        // construct the site ends, node locations, and haplotype site location indices
        // intended to be called one time after building haplotypes, after which the indices
        // are maintained automatically during edit operations
        void build_indices(vector<NestedSite>& top_level_sites);
        
        // swap the allele from one haplotype to the other, maintaining all indices
        // if a nested site is being swapped, only the topmost site should be given
        void swap_alleles(NestedSite& site, int haplotype_1, int haplotype_2);
        
        // set the allele at a site with an iterator that yields its (not including the source
        // or sink nodes of the site)
        // allele should be provided in the order indicated by the NestedSite (i.e. from start to end)
        //
        // note: function does not check that the allele path takes only edges that are actually
        // included in the graph, so applications must ensure this themselves
        template <typename NodeTraversalIterator>
        void set_allele(NestedSite& site, NodeTraversalIterator first, NodeTraversalIterator last,
                        int which_haplotype);
        
        // returns the score of the highest scoring alignment contained in the multipath alignment
        // that is restricted to the phased genome's paths through the variation graph
        //
        // note: assumes that MultipathAlignment has 'start' field filled in
        int32_t optimal_score_on_genome(const MultipathAlignment& multipath_aln);
        
    private:
        
        struct HaplotypeNode;
        class Haplotype;
        
        // graph though which the genome is threaded
        VG& graph;
        
        // all haplotypes in the genome (generally 2 per chromosome)
        vector<Haplotype> haplotypes;
        
        // an index of where nodes from the graph occur in the graph
        unordered_map<int64_t, list<HaplotypeNode*> > node_locations;
        
        // sufficient for these purposes to maintain only node ids instead of node sides
        // since the path must go through the site either before or after entering here
        unordered_map<int64_t, NestedSite*> site_starts;
        unordered_map<int64_t, NestedSite*> site_ends;
        
        // helper function
        void build_site_indices_internal(NestedSite& site);
        
        // no safety checks that the adjacent nodes aren't null (should only be used in interior of haplotype)
        // these operations maintain the node location indices but no others
        inline void insert_left(NodeTraversal node_traversal, HaplotypeNode* haplo_node);
        inline void insert_right(NodeTraversal node_traversal, HaplotypeNode* haplo_node);
        inline void remove(HaplotypeNode* haplo_node);
        
        void swap_label(NestedSite& site, Haplotype& haplotype_1, Haplotype& haplotype_2);
        
    };
    
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
        
    public:
        
        PhasedGenome::HaplotypeNode* left_telomere_node;
        PhasedGenome::HaplotypeNode* right_telomere_node;
        
        // assumes pairs are given such that .first->prev and .second->next are outside the bubble
        unordered_map<NestedSite, pair<HaplotypeNode*, HaplotypeNode*>, HashSite> sites;
        
        // construct a haplotype with a single node
        Haplotype(NodeTraversal node_traversal);
        
        // construct a haplotype with a list of nodes
        template <typename NodeTraversalIterator>
        Haplotype(NodeTraversalIterator first, NodeTraversalIterator last);
        
        ~Haplotype();
        
        // intended for use by the PhasedGenome to build haploytpes initially, returns new node
        inline HaplotypeNode* append_left(NodeTraversal node_traversal);
        inline HaplotypeNode* append_right(NodeTraversal node_traversal);
        
    };
}



#endif /* phased_genome_hpp */
