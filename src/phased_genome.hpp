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
#include <algorithm>
#include <vg/vg.pb.h>
#include "vg.hpp"
#include "nodetraversal.hpp"
#include "genotypekit.hpp"
#include "hash_map.hpp"
#include "snarls.hpp"
#include "multipath_alignment.hpp"
#include "statistics.hpp"


using namespace std;

namespace vg {
    
    /** 
     * A collection of haplotypes that represent all of the chromosomes of a genome (including
     * phasing) as walks through a variation graph. Designed for fast editing at a site level,
     * so it maintains indices of sites for that purpose.
     *
     */
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
        
        /// Constructor
        PhasedGenome(const SnarlManager& snarl_manager);
        ~PhasedGenome();
        /// overloaded constructor
        PhasedGenome(PhasedGenome& phased_genome);
        /// move assignment ctor
        PhasedGenome(PhasedGenome&& other) = delete;
        /// move assignment operator 
        PhasedGenome& operator =(PhasedGenome&& phased_genome) = delete;

        /// copy assignment operator
        PhasedGenome& operator =(PhasedGenome& phased_genome);
        
        /// Build a haplotype in place from an iterator that returns NodeTraversal objects from its
        /// dereference operator (allows construction without instantiating the haplotype elsewhere)
        /// returns the numerical id of the new haplotype
        ///
        /// note: the haplotype must have at least one node
        template <typename NodeTraversalIterator>
        int add_haplotype(NodeTraversalIterator first, NodeTraversalIterator last);
        
        // TODO: make an interface for constructing a set of homologous haplotypes at the same time
        // from a single iterator so that they don't need to be instantiated anywhere else
        
        /// Construct the site ends, node locations, and haplotype site location indices. This method
        /// is intended to be called one time after building haplotypes. After this, the are maintained
        /// automatically during edit operations.
        void build_indices();
        
        /*
         *  ITERATION AND UTILITY METHODS
         */
        
        size_t num_haplotypes();
        
        /// Unidirectional iterator starting at the left telomere and moving to the right. Requires
        /// a haplotype id as input
        iterator begin(int which_haplotype);
        /// Iterator representing the past-the-last position of the given haplotype, with the last
        /// position being the right telomere node.
        iterator end(int which_haplotype);

        /// Check which haplotypes a snarl is found in 
        // Returns a list of haplotype IDs 
        vector<id_t> get_haplotypes_with_snarl(const Snarl* snarl_to_find);

        /// Prints out the haplotypes, and node values
        void print_phased_genome();
        
        /*
         *  HAPLOTYPE EDITING METHODS
         */
        
        /// Swap the allele from between two haplotypes, maintaining all indices. If a nested site is
        /// being swapped, this method should be called only once for the top-most site. Child sites
        /// are swapped along with the top-most site automatically.
        void swap_alleles(const Snarl& site, int haplotype_1, int haplotype_2);
        
        /// Set the allele at a site with an iterator that yields its node sequence. The allele should be
        /// provided in the order indicated by the Snarl (i.e. from start to end) and it should not
        /// include the boundary nodes of the Snarl.
        ///
        /// Note: This method does not check that the allele path takes only edges that are actually
        /// included in the graph, so client must ensure this itself.
        template <typename NodeTraversalIterator>
        void set_allele(const Snarl& site, NodeTraversalIterator first, NodeTraversalIterator last,
                        int which_haplotype);
        
        /// Returns a vector of node traversals representing the allele at the indicated site on
        /// the indicated haplotype. The allele is returned in the orientation of the Snarl (start
        /// to end), not necessarily the orientation on the haplotype. The start and end nodes
        /// of the Snarl are not included in the allele.
        vector<NodeTraversal> get_allele(const Snarl& site, int which_haplotype);
        
        /// Returns the score of the highest scoring alignment contained in the multipath alignment
        /// that is restricted to the phased genome's paths through the variation graph.
        ///
        /// Note: assumes that multipath_alignment_t has 'start' field filled in
        int32_t optimal_score_on_genome(const multipath_alignment_t& multipath_aln, VG& graph);
        
        // TODO: make a local subalignment optimal score function (main obstacle is scoring partial subpaths)
        
        /// Returns the sum of the log-likelihoods of all of the alignments expressed in a multipath
        /// alignment, given a
        double read_log_likelihood(const multipath_alignment_t& multipath_aln, double log_base);
        
    private:
        
        struct HaplotypeNode;
        class Haplotype;
        
        const SnarlManager* snarl_manager;
        
        /// All haplotypes in the genome (generally 2 per chromosome)
        vector<Haplotype*> haplotypes;
        
        /// Index of where nodes from the graph occur in the phased genome
        unordered_map<int64_t, list<HaplotypeNode*> > node_locations;
        
        /// Index of which nodes are starts of Snarls
        unordered_map<int64_t, const Snarl*> site_starts;
        /// Index of which nodes are ends of Snarls
        unordered_map<int64_t, const Snarl*> site_ends;
        // note: sufficient for these purposes to maintain only node ids instead of node sides
        // since the path must go through the site either before or after entering here
        
        // Helper function
        void build_site_indices_internal(const Snarl* snarl);
        
        // Editing methods:
        // note: no safety checks that the adjacent nodes aren't null (should only be used in interior
        // of haplotype) these operations maintain the node location indices but no others
        
        /// Insert a node traversal to the left of this haplotype node and update indices.
        inline void insert_left(NodeTraversal node_traversal, HaplotypeNode* haplo_node);
        /// Insert a node traversal to the right of this haplotype node and update indices.
        inline void insert_right(NodeTraversal node_traversal, HaplotypeNode* haplo_node);
        /// Remove this haplotype node from its haplotype and update indices.
        inline void remove(HaplotypeNode* haplo_node);
        
        /// Update a subsite's location in indices after swapping its parent allele
        void swap_label(const Snarl& site, Haplotype& haplotype_1, Haplotype& haplotype_2);
        
    };
    
    /*
     *  INTERNAL OBJECTS (NO PUBLIC-FACING API BELOW HERE)
     */
    
    /**
     * A node in walk through the graph taken by a haplotype.
     *
     */
    struct PhasedGenome::HaplotypeNode {
    public:
        
        /// Node and strand
        NodeTraversal node_traversal;
        /// Next node in walk
        HaplotypeNode* next;
        /// Previous node in walk
        HaplotypeNode* prev;
        
        /// Constructor
        HaplotypeNode(NodeTraversal node_traversal, HaplotypeNode* next, HaplotypeNode* prev);
        /// Destructor
        ~HaplotypeNode();
        
    };
    
//    /// Hash function for node traversals for use in indices
//    struct HashNodeTraversal {
//        size_t operator()(const NodeTraversal& node_traversal) const {
//            return hash(make_pair(node_traversal.node, node_traversal.backward));
//        }
//    };
//    
//    /// Hash function for sites for use in indices
//    struct HashSite {
//        size_t operator()(const NestedSite& site) const {
//            return hash(make_pair(make_pair(site.start.node, site.start.backward),
//                                  make_pair(site.end.node, site.end.backward)));
//        }
//    };
    
    /**
     * Specialized linked list that tracks a walk through the variation graph and maintains
     * an index of sites.
     *
     */
    class PhasedGenome::Haplotype {
        
    private:
        /// Leftmost node in walk
        PhasedGenome::HaplotypeNode* left_telomere_node;
        /// Rightmost node in walk
        PhasedGenome::HaplotypeNode* right_telomere_node;
        
        /// Index of the location in the haplotype of nested sites. Locations of sites are stored
        /// as the nodes on haplotype that correspond to the start and end node of the site. The pair
        /// of haplotype nodes is stored in left-to-right order along the haplotype (i.e. .first->prev
        /// and .second->next are outside the bubble).
        unordered_map<const Snarl*, pair<HaplotypeNode*, HaplotypeNode*> > sites;
        
    public:
        /// Construct a haplotype with a single node
        Haplotype(NodeTraversal node_traversal);
        
        /// Construct a haplotype with an iterator that yields NodeTraversals
        template <typename NodeTraversalIterator>
        Haplotype(NodeTraversalIterator first, NodeTraversalIterator last);
        
        ~Haplotype();
        
        /// Add a haplotype node for this node traversal to left end of haplotye and return new node.
        /// Does not maintain indices; intended for use by the PhasedGenome for initial haplotype build.
        inline HaplotypeNode* append_left(NodeTraversal node_traversal);
        /// Add a haplotype node for this node traversal to right end of haplotye and return new node.
        /// Does not maintain indices; intended for use by the PhasedGenome for initial haplotype build.
        inline HaplotypeNode* append_right(NodeTraversal node_traversal);
        
        friend class PhasedGenome;
        friend class HaplotypeNode;
        
    };
    
    /**
     * Unidirectional iterator to obtain the NodeTraversals of a haplotype. Can be come invalid
     * if the PhasedGenome is edited while iterating.
     *
     */
    class PhasedGenome::iterator {
    private:
        
        /// Ordinal position along the haplotype (to distinguish the same node with multiple copies)
        size_t rank;
        /// The ID of the haplotype
        int haplotype_number;
        /// The position along the haplotype
        HaplotypeNode* haplo_node;
        
        iterator(size_t rank, int haplotype_number, HaplotypeNode* haplo_node);
        
    public:

        using value_type = NodeTraversal;   

        /// Default constructor
        iterator();
        /// Copy constructor
        iterator(const iterator& other);
        /// Destructor
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
    void PhasedGenome::set_allele(const Snarl& site, NodeTraversalIterator first, NodeTraversalIterator last,
                                  int which_haplotype) {
#ifdef debug_phased_genome
        cerr << "[PhasedGenome::set_allele]: setting allele on haplotype " << which_haplotype << endl;
#endif
        Haplotype& haplotype = *haplotypes[which_haplotype];
        
        // can only set the allele of a site that already is in the haplotype
        assert(haplotype.sites.count(&site));
        
        pair<HaplotypeNode*, HaplotypeNode*> haplo_site = haplotype.sites[&site];
        
#ifdef debug_phased_genome
        cerr << "[PhasedGenome::set_allele]: deleting allele at site " << haplo_site.first->node_traversal.node->id() << "->" << haplo_site.second->node_traversal.node->id() << endl;
#endif
        
        // remove the current site
        HaplotypeNode* haplo_node = haplo_site.first->next;
        while (haplo_node != haplo_site.second) {
            // don't need to worry about erasing the same site twice (once on start and
            // once on end) since unordered_map.erase is defined in both cases
            int64_t node_id = haplo_node->node_traversal.node->id();
            if (site_starts.count(node_id)) {
#ifdef debug_phased_genome
                cerr << "[PhasedGenome::set_allele]: deleting nested site starting at node " << node_id << " from index" << endl;
#endif
                const Snarl* subsite = site_starts[node_id];
                haplotype.sites.erase(subsite);
            }
            else if (site_ends.count(node_id)) {
#ifdef debug_phased_genome
                cerr << "[PhasedGenome::set_allele]: deleting nested site ending at node " << node_id << " from index" << endl;
#endif
                const Snarl* subsite = site_ends[node_id];
                haplotype.sites.erase(subsite);
            }
#ifdef debug_phased_genome
            cerr << "[PhasedGenome::set_allele]: deleting haplotype node " << node_id << endl;
#endif
            // cut out each node and iterate to next node
            HaplotypeNode* next_haplo_node = haplo_node->next;
            remove(haplo_node);
            haplo_node = next_haplo_node;
        }
        
        // is site in forward or reverse direction on haplotype?
        bool forward = (haplo_site.first->node_traversal.node->id() == site.start().node_id());
        
        // orient traversal through the allele to traversal along the haplotype
        haplo_node = forward ? haplo_site.first : haplo_site.second;
        
        // keeps track of the haplotype node where we entered a site and thereby also
        // indicates whether we have entered the site yet
        unordered_map<const Snarl*, HaplotypeNode*> subsite_start_side;
        unordered_map<const Snarl*, HaplotypeNode*> subsite_end_side;
        
        // start inserting nodes from the first position in the allele
        for (; first != last; first++) {
            
            // create the next node in the allele and move the pointer onto the new node
            if (forward) {
#ifdef debug_phased_genome
                cerr << "[PhasedGenome::set_allele]: inserting in forward direction " << (*first).node->id() << ((*first).backward ? "-" : "+") << endl;
#endif
                insert_right(*first, haplo_node);
                haplo_node = haplo_node->next;
            }
            else {
#ifdef debug_phased_genome
                cerr << "[PhasedGenome::set_allele]: inserting in reverse direction " << (*first).node->id() << ((*first).backward ? "-" : "+") << endl;
#endif
                insert_left(*first, haplo_node);
                haplo_node = haplo_node->prev;
            }
            
            int64_t node_id = haplo_node->node_traversal.node->id();
            
            // does a site start here?
            if (site_starts.count(node_id)) {
                const Snarl* subsite = site_starts[node_id];
                // are we entering or leaving this site?
                if (subsite_end_side.count(subsite)) {
#ifdef debug_phased_genome
                    cerr << "[PhasedGenome::set_allele]: detected leaving a site in its reverse orientation" << endl;
#endif
                    // add this site into the haplotype's site index
                    HaplotypeNode* other_side_node = subsite_end_side[subsite];
                    // the site sides should always be entered in the order that
                    // they occur in the haplotype
                    haplotype.sites[subsite] = forward ? make_pair(other_side_node, haplo_node)
                                                        : make_pair(haplo_node, other_side_node);
                }
                else {
#ifdef debug_phased_genome
                    cerr << "[PhasedGenome::set_allele]: detected entering a site in its forward orientation" << endl;
#endif
                    // we are entering a site, mark the location of the entrance
                    subsite_start_side[subsite] = haplo_node;
                }
            }
            
            // does a site end here?
            if (site_ends.count(node_id)) {
                const Snarl* subsite = site_ends[node_id];
                // are we entering or leaving this site?
                if (subsite_start_side.count(subsite)) {
#ifdef debug_phased_genome
                    cerr << "[PhasedGenome::set_allele]: detected leaving a site in its forward orientation" << endl;
#endif
                    // add this site into the haplotype's site index
                    HaplotypeNode* other_side_node = subsite_start_side[subsite];
                    // the site sides should always be entered in the order that
                    // they occur in the haplotype
                    haplotype.sites[subsite] = forward ? make_pair(other_side_node, haplo_node)
                                                        : make_pair(haplo_node, other_side_node);
                }
                else {
#ifdef debug_phased_genome
                    cerr << "[PhasedGenome::set_allele]: detected entering a site in its reverse orientation" << endl;
#endif
                    // we are entering a site, mark the location of the entrance
                    subsite_end_side[subsite] = haplo_node;
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
        
#ifdef debug_phased_genome
        cerr << "[PhasedGenome::insert_left]: recording inserted node " << node_traversal.node->id() << " in haplotype node at " << haplo_node << endl;
#endif
        
        node_locations[node_traversal.node->id()].push_back(new_node);
    }
    
    inline void PhasedGenome::insert_right(NodeTraversal node_traversal, HaplotypeNode* haplo_node) {
        
        HaplotypeNode* new_node = new HaplotypeNode(node_traversal, haplo_node->next, haplo_node);
        haplo_node->next->prev = new_node;
        haplo_node->next = new_node;
        
#ifdef debug_phased_genome
        cerr << "[PhasedGenome::insert_right]: recording inserted node " << node_traversal.node->id() << " in haplotype node at " << haplo_node << endl;
#endif
        node_locations[node_traversal.node->id()].push_back(new_node);
    }
    
    inline void PhasedGenome::remove(HaplotypeNode* haplo_node) {
        // note: updating site indices is handled outside this function
        
        haplo_node->next->prev = haplo_node->prev;
        haplo_node->prev->next = haplo_node->next;
        
        // remove the node from the node locations index
        int64_t node_id = haplo_node->node_traversal.node->id();
        
        list<HaplotypeNode*>& node_occurrences = node_locations[node_id];
#ifdef debug_phased_genome
        cerr << "[PhasedGenome::remove]: node " << node_id << " occurs in " << node_occurrences.size() << " places, must search through each to find which to delete from indices" << endl;
#endif
        for (auto iter = node_occurrences.begin(); iter != node_occurrences.end(); iter++) {
            if (*iter == haplo_node) {
                node_occurrences.erase(iter);
                break;
            }
        }
        
        delete haplo_node;
    }
}

namespace std{
    template<>
    struct iterator_traits<vg::PhasedGenome::iterator>{
        using value_type = vg::NodeTraversal;   
        using iterator_category = forward_iterator_tag;
    };
}

    


#endif /* phased_genome_hpp */
