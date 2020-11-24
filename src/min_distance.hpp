#ifndef VG_MIN_DISTANCE_HPP_INCLUDED
#define VG_MIN_DISTANCE_HPP_INCLUDED

#include "snarls.hpp"

using namespace sdsl;
namespace vg { 

/**
 * The distance index. Stores minimum distances among nodes in each 
 * netgraph and chain.
 */

class MinimumDistanceIndex {

    public:

    /** 
     * A structure representing a node/snarl/chain in a snarl tree
     */
    //TODO: This should probably be packed into a size_t
    //TODO: This could also be a net_handle_t, since I think a net_handle_t would just need this
    //      plus the orientation of traversal
    struct snarl_tree_node_t {

        //This is a node id if it represents a node, component number if it is a root,
        //otherwise the index of the snarl or chain (index into snarl/chain_records)
        //Component number is the index of a component into root_snarl
        size_t id;

        //Which end of the node has edges (within the parent)?
        //Boundary nodes of snarls and tips within a snarl will have left/right_only
        enum connectivity {LEFT_ONLY, RIGHT_ONLY, BOTH, NONE};


        //What kind of node is this? 
        //Root is only for the parent field of a snarl or chain record
        enum node_type {NODE, SNARL, CHAIN, ROOT};

    }
 
    /**
     * Structure for accessing indexes and interpreting a vector of size_t as a node, snarl, 
     * or chain record
     */ 
    class DistanceRecord {
        //TODO: Or maybe the net_handle_t should be a pared down version of this
        //Or since the SnarlDecomposition already needs a bunch of functions taking net_handle_t's,
        //maybe this shouldn't exist and each of these functions should be separate and parse the 
        //values from the snarl_tree_node_t individually
        //
        //TODO: might want to be able to make this for a trivial snarl even though we don't have an actual record for it

        public:


        //Accessors for properties of the snarl tree
        virtual size_t get_node_count();
        virtual int64_t get_min_length();
        virtual int64_t get_max_length();
        virtual snarl_tree_node_t get_parent();
        virtual distance_record_t get_parent_record();
        virtual size_t get_snarl_tree_depth();
        //TODO: Some sort of iterator for children maybe

        //Distance calculations
        virtual int64_t distance_between_nodes (id_t node1, bool rev1, id_t node2, bool rev2);
        virtual int64_t distance_to_start (id_t node, bool rev);
        virtual int64_t distance_to_end (id_t node, bool rev);

        private:

        snarl_tree_node_t self;  //What node/snarl/chain this record represents
                             //This also knows where the record starts, and can find where it ends
                             //based on the number of nodes
                             //TODO: Be careful not to access past the end of the record
        }
    }
    class SnarlRecord : DistanceRecord {
    }
    class ChainRecord : DistanceRecord {
        //TODO: This also needs to be able to access each of the loop distances for clustering
    }



    protected:

    size_t minimum_node_id;
    size_t max_depth;

    /**
     * All root-level structures. Can contain nodes, snarls, and chains.
     * Connectivity will always be "none" (Or I could just take it out)
     */
    vector<snarl_tree_node_t> root_snarl;


    /**
     * This stores each of the snarls as a snarl record, which will be vector<size_t>'s
     *
     * Each snarl record contains:
     * (
     *  # of nodes, 
     *  parent (a snarl_tree_node_t minus connectivity, node_type=root if it has no parent),
     *  depth in snarl tree (starting at 0 as root and including chains and snarls),
     *  maximum length of snarl (start to end, including boundary nodes), 
     *  bool is simple snarl (If it is simple, don't store distances), //TODO: Need to be careful about orientation of child nodes
     *  bool is unary snarl, TODO: Do we still have unary snarls? Or snarl where the start and end are the same node in different directions
     *  [list of child nodes as snarl_tree_node_t's]
     *  [distances]
     *  )
     *
     *  The list of child nodes will be unordered, except that the start and end boundary nodes
     *  will be first and last
     *  The orientation of the boundary nodes can be found from the list of snarl_tree_node_t's,
     *  since they won't allow edges going out of the snarl, and all other nodes we'll assume to be
     *  oriented left (or start node according to child) to right (or end node)
     *  TODO: I think this will be fine
     *
     *  The distances will be ordered based on the order of the child nodes
     *  Distances are offset by 1, 0 is used for no path between nodes
     *
     *
     *  TODO: In the current distance index, each node in the snarl also needs to know its rank (to
     *  know where to find the correct distance in the distance vector), but  here I"m only 
     *  storing a list of nodes, then I plan on walking through the list to get the  ranks. It might
     *  be a bit slow but we're already relying on the snarls being small so I don't think it will
     *  be that bad
     */
    vector<size_t> snarl_records;



    /**
     * This stores all of the chains as chain records
     *
     * Each chain record is a vector<size_t> representing:
     *
     * (
     *  # nodes,
     *  parent (a snarl_tree_node_t minus connectivity, node_type=root if it has no parent),
     *  depth in snarl tree,
     *  maximum length of chain,
     *  is looping chain (true if the start and end node are the same),
     *  ((node id, prefix sum value, loop left, loop right) bool has snarl )xN
     * )
     *
     * The distances are stored as a list of nodes (id and maybe distances), interspersed
     * with booleans indicating whether there is a snarl between two nodes
     *
     * TODO: Could store the index of the snarl instead of a bool, if it's a bool we can look up the
     * snarl from either node id - it would be the primary snarl of the node that came first in the 
     * chain and the secondary snarl of the node that came second
     *
     */
    vector<size_t> chain_records;



    /**
     * Assigns each node to a snarl/secondary snarl or chain or to the root if it's not connected to
     * anything. This looks up one level in the snarl tree
     *
     * A node can have a secondary snarl if it is involved in more than one snarl as a boundary of
     * two snarls or if it represents a child of a snarl 
     * TODO: We might not need it to know the grandparent snarl, just the parent one
     *
     * Each node assignment record contains:
     *
     * ( 
     *  bool in snarl,            //True if this node is in a snarl, false if it is a root
     *  bool in trivial snarl,    //True if the snarl it is in is trivial
     *  node length,
     *  primary snarl index         OR   chain index     OR  component number
     *  secondary snarl index       OR   rank in chain   OR  -
     *  
     * )
     *
     * Stores primary and secondary snarl indexes if (in_snarl && !in_trivial_snarl), 
     *        chain index and rank if (in_snarl && in_trivial_snarl), 
     *        component number if (!in_snarl)
     *
     * Snarl/chain index is the index of the start of the snarl/chain record in snarl/chain_records
     * The length of the entire record can be found from the first value (# nodes), and the record
     * will be a slice of the whole vector
     *  
     * Nodes are indexed as node_id - minimum_node_id
     *
     */
    vector<size_t> node_records;


};

}

#endif
