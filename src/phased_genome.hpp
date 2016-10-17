//
//  phased_genome.hpp
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
    
    // TODO: only keep track of paths here, leave probability calculations and sampling
    // to other modules
    
    struct HaplotypeNode {
    public:
        NodeTraversal node_traversal;
        HaplotypeNode* next;
        HaplotypeNode* prev;
        
        HaplotypeNode(NodeTraversal node_traversal, HaplotypeNode* next, HaplotypeNode* prev);
        ~HaplotypeNode();
        
        
    }
    
    // forward declaration
    struct HashSite;
    
    class Haplotype {
        
    public:
        
        HaplotypeNode* left_telomere_node;
        HaplotypeNode* right_telomere_node;
        
        // assumes pairs are given such that .first->prev and .second->next are outside the bubble
        unordered_map<NestedSite, pair<HaplotypeNode*, HaplotypeNode*>, HashSite> sites;
        
        Haplotype(NodeTraversal node_traversal);
        Haplotype(list<NodeTraversal>& node_path);
        ~Haplotype();
        
        // intended for use by the PhasedGenome to built haploytpes initially
        HaplotypeNode* append_left(NodeTraversal node_traversal);
        HaplotypeNode* append_right(NodeTraversal node_traversal);
        
    };
    
    class PhasedGenome {
    private:
        
        vector<Haplotype> haplotypes;
        
        unordered_map<int64_t, vector<HaplotypeNode*> > node_locations;
        
        const VG& graph;
        
        // sufficient for these purposes to maintain only node ids instead of node sides
        // since the path must go through the site either before or after entering here
        unordered_map<int64_t, NestedSite*> site_starts;
        unordered_map<int64_t, NestedSite*> site_ends;
        
        // construct the site ends, node locations, and haplotype site location indices
        // intended to be called one time after building haplotypes, after which the indices
        // are maintained automatically
        void build_site_indices(vector<NestedSite>& top_level_sites);
        // helper function for the above
        void build_site_indices_internal(NestedSite& site);
        
        // no safety checks that the adjacent nodes aren't null (should only be used in interior of haplotype)
        void insert_left(NodeTraversal node_traversal, HaplotypeNode* haplo_node);
        void insert_right(NodeTraversal node_traversal, HaplotypeNode* haplo_node);
        void remove(HaplotypeNode* haplo_node);
        
        void swap_label(NestedSite& site, Haplotype& haplotype_1, Haplotype& haplotype_2);
        
    public:
        
        PhasedGenome(const VG& graph);
        ~PhasedGenome();
        
        // TODO: Do I also need to keep track of the ploidy of sites in case the path includes them
        // multiple times?
        
        // build a haplotype in place from an iterator that returns NodeTraversal objects from its
        // dereference operator
        // note: the haplotype must have at least one node
        template <typename NodeTraversalIterator>
        void add_haplotype(NodeTraversalIterator first, NodeTraversalIterator last);
        
        // swap the allele from one haplotype to the other, maintaining all indices
        // if a nested site is being swapped, only the topmost site should be given
        void swap_alleles(NestedSite& site, int haplotype_1, int haplotype_2);
        
        // set the allele at a site (allele should not include the source or sink nodes of the site)
        void set_allele(NestedSite& site, list<NodeTraversal>& allele, int which_haplotype);
        
        // returns the score of the highest scoring alignment contained in the multipath alignment
        // that is restricted to the phased genome's paths through the variation graph
        //
        // note: assumes that MultipathAlignment has 'start' field filled in
        //
        //  Args:
        //    multipath_aln     multipath alignment to score
        //
        int32_t optimal_score_on_genome(const MultipathAlignment& multipath_aln);
        
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
    
}



#endif /* phased_genome_hpp */
