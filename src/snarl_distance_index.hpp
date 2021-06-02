#ifndef VG_SNARL_DISTANCE_HPP_INCLUDED
#define VG_SNARL_DISTANCE_HPP_INCLUDED

//#define debug_indexing

#include <handlegraph/snarl_decomposition.hpp>
#include <structures/union_find.hpp>
#include "position.hpp"
#include "snarls.hpp"


using namespace sdsl;
using namespace handlegraph;
namespace vg { 

/**
 * The distance index. Stores minimum distances among nodes in each 
 * netgraph and chain.
 * Also used to store the snarl tree
 */


class SnarlDistanceIndex : public SnarlDecomposition {

public:
    ~SnarlDistanceIndex();
    SnarlDistanceIndex();

    //Fill in the index
    SnarlDistanceIndex(const HandleGraph* graph, const HandleGraphSnarlFinder* snarl_finder, size_t size_limit = 500);

private:
    
    /**
     *
     * This stores all records for the root, nodes, chains/snarls, and snarls' children
     * 
     * It's really made up of five types of vectors: 
     * 
     * - The (single) root vector has the format:
     *   [root tag, # connected components, # nodes, min_node_id, [pointer to node/snarl/chain record] x N]]
     *   The root vector stores the root of every connected component, which can be a 
     *   node, snarl, or chain
     * 
     * - The (single) node vector stores a record for each node and has the format:
     *   [[node tag, pointer to parent, node length, rank in parent,
     *       component #]] x N TODO: Could add extra stuff here, I'm not sure where to define it
     *   Rank in a chain is actually the offset in the chain, so it points to that node in the chain
     * 
     * - A chain record for each chain, which is interspersed with snarl records:
     *   [chain tag, #nodes, pointer to parent, min length, max length, rank in parent, start, end, 
     *       [node id, prefix sum, fd loop, rev loop, snarl record size, (snarl record, snarl record size)] x N] 
     *          (plus an extra node id, prefix sum, fd loop, rev loop at the end for the last node), 0, 0
     *    The two 0's at the end indicates that it's the end of the chain
     *    snarl_record_size is the number of things in the snarl record (so 0 if there is no snarl there)
     *    start/end include the orientations
     * 
     * - A snarl record for each snarl, which are stuck in chains
     *   [snarl tag, # nodes, pointer to parent, min length, max length, rank in parent, 
     *      pointer to children (in child vector), start, end, distance vector]
     *   Trivial snarls will still be contained within a chain
     *   Single nodes will technically be contained in a chain, but not stored as a chain
     *   The rank of the start node will always be 0 and the end node rank will be 1
     *   start/end are the start and end nodes, include the orientations
     *   For the first and last nodes, we only care about the node sides pointing in
     *   Each node side in the actual distance matrix will be 2*(rank-1) for the left side, and 
     *   2rank+1 for the right side, and 0 for the start, 1 for the end, where we only keep the 
     *   inner node side of the start and end
     *   Node count is the number of nodes, not including boundary nodes
     * 
     * - The (single) child vector, listing children in snarls
     *   [child vector tag, (pointer to records) x N
     *   Each snarl will have a pointer into here, and will also know how many children it has
     * 
     * 
     *   For each of the "rank_in_parent" fields, the last bit is 1 if it is reversed in the parent
     *
     */
    
    vector<size_t> snarl_tree_records;

    /* If this is 0, then don't store distances
     * Otherwise, for snarls with more children than snarl_size_limit, only store the distances
     * that include boundary nodes (OVERSIZED_SNARL)
     */
    size_t snarl_size_limit = 500;


    
     //TODO: Get rid of this
     public:
        void print_self() {
            cerr << "Entire index:" << endl;
            for (size_t i = 0 ; i < snarl_tree_records.size() ; i++) {
                cerr << "(" << i << "):" << snarl_tree_records[i] << "   " ;
            }
            cerr << endl;
        }
    private:

    /*
     *
     * The "tags" for defining what kind of record we're looking at will be a record_t and a 
     * bit vector indicating connectivity. The bit vector will be the last 6 bits of the tag
     * 
     * Each bit represents one type of connectivity:
     * start-start, start-end, start-tip, end-end, end-tip, tip-tip
     * std::bitset<6> connectivity;
     * 
     * The remainder of the tag will be the record_t of the record
     * NODE, SNARL, and CHAIN indicate that they don't store distances.
     * SIMPLE_SNARL is a snarl with all children connecting only to the boundary nodes
     * OVERSIZED_SNARL only stores distances to the boundaries
     * ROOT_SNARL represents a connected component of the root. It has no start or end node so 
     *    its children actually belong to the root
     * MULTICOMPONENT_CHAIN can represent a chain with snarls that are not start-end connected
     *     The chain is split up into components between these snarls, each node is tagged with
     *     which component it belongs to
     */
    //TODO: Make simple snarls work
    //TODO: Maybe also add a tag for trivial snarls/chains
    //TODO: Unary snarls? Looping chains?
    enum record_t {ROOT=1, 
                   NODE, DISTANCED_NODE, 
                   SNARL, DISTANCED_SNARL, SIMPLE_SNARL, OVERSIZED_SNARL, ROOT_SNARL, DISTANCED_ROOT_SNARL,
                   CHAIN, DISTANCED_CHAIN, MULTICOMPONENT_CHAIN,
                   CHILDREN};

    enum connectivity_t { START_START=1, START_END, START_TIP, 
                            END_START, END_END, END_TIP, 
                            TIP_START, TIP_END, TIP_TIP};
    //Type of a net_handle_t. This is to allow a node record to be seen as a chain from the 
    //perspective of a handle
    enum net_handle_record_t {ROOT_HANDLE=0, NODE_HANDLE, SNARL_HANDLE, CHAIN_HANDLE, SENTINEL_HANDLE};

private:
    /*Give each of the enum types a name for debugging */
    vector<string> record_t_as_string = {"ROOT", "NODE", "DISTANCED_NODE", 
                     "SNARL", "DISTANCED_SNARL", "SIMPLE_SNARL", "OVERSIZED_SNARL", 
                     "ROOT_SNARL", "DISTANCED_ROOT_SNARL",
                     "CHAIN", "DISTANCED_CHAIN", "MULTICOMPONENT_CHAIN",
                     "CHILDREN"};
    vector<string> connectivity_t_as_string = { "START_START", "START_END", "START_TIP", 
                            "END_START", "END_END", "END_TIP", 
                            "TIP_START", "TIP_END", "TIP_TIP"};
    vector<string> net_handle_record_t_string = {"ROOT_HANDLE", "NODE_HANDLE", "SNARL_HANDLE", 
                                                "CHAIN_HANDLE", "SENTINEL_HANDLE"};

public:


////////////////// SnarlDecomposition methods

    /**
     * Get a net handle referring to a tip-to-tip traversal of the contents of the root snarl.
     * TODO: Special handling for circular things in the root snarl? Circular traversal type?
     */
    net_handle_t get_root() const ;
    
    /**
     * Return true if the given handle refers to (a traversal of) the root
     * snarl, and false otherwise.
     */
    bool is_root(const net_handle_t& net) const;
    
    /**
     * Returns true if the given net handle refers to (a traversal of) a snarl.
     */
    bool is_snarl(const net_handle_t& net) const;
    /**
     * Returns true if the given net handle refers to (a traversal of) a chain.
     */
    bool is_chain(const net_handle_t& net) const;
    /**
     * Returns true if the given net handle refers to (a traversal of) a trivial chain that represents a single node.
     */
    bool is_trivial_chain(const net_handle_t& net) const;
    /**
     * Returns true if the given net handle refers to (a traversal of) a single node, and thus has a corresponding handle_t.
     */
    bool is_node(const net_handle_t& net) const;
    /**
     * Return true if the given net handle is a snarl bound sentinel (in either
     * inward or outward orientation), and false otherwise.
     */
    bool is_sentinel(const net_handle_t& net) const;
    
    /**
     * Turn a handle to an oriented node into a net handle for a start-to-end or end-to-start traversal of the node, as appropriate.
     */
    net_handle_t get_net(const handle_t& handle, const handlegraph::HandleGraph* graph) const;
    
    /**
     * For a net handle to a traversal of a single node, get the handle for that node in the orientation it is traversed.
     * May not be called for other net handles.
     */
    handle_t get_handle(const net_handle_t& net, const handlegraph::HandleGraph* graph) const;
    
    /**
     * Get the parent snarl of a chain, or the parent chain of a snarl or node.
     * If the child is start-to-end or end-to-start, and the parent is a chain,
     * the chain comes out facing the same way, accounting for the relative
     * orientation of the child snarl or node in the chain. Otherwise,
     * everything is produced as start-to-end, even if that is not actually a
     * realizable traversal of a snarl or chain. May not be called on the root
     * snarl.
     *
     * Also works on snarl boundary sentinels.
     */
    net_handle_t get_parent(const net_handle_t& child) const;
    
    
    // We have sentinel net_handle_t values for the start/end of each snarl, so
    // that we can tell which last edge a traversal of the contents of a snarl
    // takes when we represent it as a list of net_handle_t items. We also use
    // these to query what's attached to the snarl start/end when traversing,
    // and to see self loops immediately inside the snarl. These may actually
    // just be the handles for the nodes at the start'end of the snarl with
    // special flags set.
    //
    // For chains, we use the net handles to the appropriate first/last nodes
    // in the appropriate orientation.
    
    /**
     * Get the bounding handle for the snarl or chain referenced by the given
     * net handle, getting the start or end facing in or out as appropriate.
     *
     * For snarls, returns the bounding sentinel net handles. For chains,
     * returns net handles for traversals of the bounding nodes of the chain.
     *
     * Ignores traversal type.
     *
     * May not be called on traversals of individual nodes.
     */
    net_handle_t get_bound(const net_handle_t& snarl, bool get_end, bool face_in) const;
    
    /**
     * Return a net handle to the same snarl/chain/node in the opposite orientation.
     * No effect on tip-to-tip, start-to-start, or end-to-end net handles. Flips all the others.
     */
    net_handle_t flip(const net_handle_t& net) const;
    
    /**
     * Get a canonical traversal handle from any net handle. All handles to the
     * same net graph element have the same canonical traversal. That canonical
     * traversal must be realizable, and might not always be start-to-end or
     * even consistently be the same kind of traversal for different snarls,
     * chains, or nodes. Mostly useful to normalize for equality comparisons.
     */
    net_handle_t canonical(const net_handle_t& net) const;

    /**
     * Return the kind of location at which the given traversal starts.
     */
    endpoint_t starts_at(const net_handle_t& traversal) const;
    
    /**
     * Return the kind of location at which the given traversal ends.
     */
    endpoint_t ends_at(const net_handle_t& traversal) const;

protected:
    /**
     * Internal implementation for for_each_child.
     */
    bool for_each_child_impl(const net_handle_t& traversal, const std::function<bool(const net_handle_t&)>& iteratee) const;

    /**
     * Internal implementation for for_each_traversal.
     */
    bool for_each_traversal_impl(const net_handle_t& item, const std::function<bool(const net_handle_t&)>& iteratee) const;

    /**
     * Internal implementation for follow_net_edges.
     */
    bool follow_net_edges_impl(const net_handle_t& here, const handlegraph::HandleGraph* graph, bool go_left, const std::function<bool(const net_handle_t&)>& iteratee) const;

public:


        /**
     * Get a net handle for traversals of a snarl or chain that contains
     * the given oriented bounding node traversals or sentinels. Given two
     * sentinels for a snarl, produces a net handle to a start-to-end,
     * end-to-end, end-to-start, or start-to-start traversal of that snarl.
     * Given handles to traversals of the bounding nodes of a chain, similarly
     * produces a net handle to a traversal of the chain.
     *
     * For a chain, either or both handles can also be a snarl containing tips,
     * for a tip-to-start, tip-to-end, start-to-tip, end-to-tip, or tip-to-tip
     * traversal. Similarly, for a snarl, either or both handles can be a chain
     * in the snarl that contains internal tips, or that has no edges on the
     * appropriate end.
     *
     * May only be called if a path actually exists between the given start
     * and end.
     */
    net_handle_t get_parent_traversal(const net_handle_t& traversal_start, const net_handle_t& traversal_end) const;


///////////////////////////// More public functions for distance calculations using net_handle_t's 


public:
    /**
     * Find the distance between the two child node sides in the parent, facing each other, not 
     * including the lengths of the nodes
     * This only takes into account the endpoint of the net_handle_t traversal, it does not care if the traversal
     * was possible. Doesn't allow you to find the distance from a traversal ending/starting in a tip
     * requires that the children are children of the parent
     * For chains, the distance to/from a snarl is really the distance from the outer node of the snarl
     * Returns std::numeric_limits<int64_t>::max() if there is not path between them in the parent 
     * or if they are not children of the parent
     */
    int64_t distance_in_parent(const net_handle_t& parent, const net_handle_t& child1, const net_handle_t& child2) const;


    /** 
     * For two net handles, get a net handle lowest common ancestor
     * If the lowest common ancestor is the root, then the two handles may be in 
     * different connected components. In this case, return false
     */
    pair<net_handle_t, bool> lowest_common_ancestor(const net_handle_t& net1, const net_handle_t& net2) const;


    /**
     * Return the length of the net, which must represent a node (or sentinel of a snarl)
     */
    int64_t node_length(const net_handle_t& net) const ;


    /**
     * This is also the length of a net, but it can also be a snarl or chain. 
     * The length of a chain includes the boundary nodes, a snarl does not
     * A looping chain only includes the start/end node once
     */
    int64_t minimum_length(const net_handle_t& net) const;

    /**
     * What is the node id of the node represented by this net handle
     * net must be a node or a sentinel
     */
    int64_t node_id(const net_handle_t& net) const ;

    bool has_connectivity(const net_handle_t& net, endpoint_t start, endpoint_t end) const ;
    bool has_external_connectivity(const net_handle_t& net, endpoint_t start, endpoint_t end) const ; 


    /**
     * Get the minimum distance between two positions in the graph
     * If unoriented_distance is true, then ignore the orientations of the positions
     * Otherwise, distance is calculated from the first position going forward to the second position going forward
     * The distance includes one of the positions; the distance from one position to itself is 1
     * Returns std::numeric_limits<int64_t>::max() if there is no path between the two positions
     */
    //TODO: The positions can't be const?
    int64_t minimum_distance(pos_t pos1, pos_t pos2, bool unoriented_distance = false) const ;



////////////////////////////// How to interpret net_handle_ts
//TODO: Does this depend on endianness???
//TODO: Should this also know what kind of node it's pointing to?
//Last 3 bits are the net_handle_record_t, next four are the connectivity_t, last are the offset into snarl_tree_records
//
//For sentinel net_handle_t's, which represent the boundaries of snarls (just inside of the boundary nodes),
//The record points to the snarl containing them, and the connectivity indicates which bound 
//we're looking at (START_END for start in, START_START for start out, etc)

private:

    const static size_t get_record_offset (const handlegraph::net_handle_t& net_handle) {
        return as_integer(net_handle) >> 7;
    }
    const static connectivity_t get_connectivity (const handlegraph::net_handle_t& net_handle){
        size_t connectivity_as_int = (as_integer(net_handle)>>3) & 15; //Get last 4 bits
        assert (connectivity_as_int <= 9);
        return static_cast<connectivity_t>(connectivity_as_int);
    }
    const static net_handle_record_t get_handle_type (const handlegraph::net_handle_t& net_handle) {
        size_t connectivity_as_int = as_integer(net_handle) & 7; //Get last 3 bits
        assert (connectivity_as_int <= 4);
        return static_cast<net_handle_record_t>(connectivity_as_int);
    }

    const static handlegraph::net_handle_t get_net_handle(size_t pointer, connectivity_t connectivity, net_handle_record_t type) {
        net_handle_t handle =  as_net_handle(((((pointer << 4) | connectivity)<<3) | type)); 
        return handle;
    
    }
    handlegraph::net_handle_t get_net_handle(size_t pointer, connectivity_t connectivity) const  {
        net_handle_record_t type = SnarlTreeRecord(pointer, &snarl_tree_records).get_record_handle_type(); 
        return get_net_handle(pointer, connectivity, type); 
    
    }


private:


/////////////////////////////// My methods for interpreting the snarl tree records
//// Uses a SnarlTreeRecord as the main class for defining and interpreting the records
    
    //TODO: Replace this once everything is set
    const static size_t ROOT_RECORD_SIZE = 4;
    const static size_t COMPONENT_COUNT_OFFSET = 1;
    const static size_t NODE_COUNT_OFFSET = 2;
    const static size_t MIN_NODE_ID_OFFSET = 3;
    //The offset of each value in snarl_tree_records, offset from the start of the record
    const static size_t NODE_RECORD_SIZE = 5;
    const static size_t NODE_PARENT_OFFSET = 1;
    const static size_t NODE_LENGTH_OFFSET = 2;
    const static size_t NODE_RANK_OFFSET = 3;
    const static size_t NODE_COMPONENT_OFFSET = 4;
    
    const static size_t SNARL_RECORD_SIZE = 9;
    const static size_t SNARL_NODE_COUNT_OFFSET = 1;
    const static size_t SNARL_PARENT_OFFSET = 2;
    const static size_t SNARL_MIN_LENGTH_OFFSET = 3;
    const static size_t SNARL_MAX_LENGTH_OFFSET = 4;
    const static size_t SNARL_RANK_OFFSET = 5;
    const static size_t SNARL_CHILD_RECORD_OFFSET = 6;
    //TODO: This could also be found from the list of the snarl's children, but probably better here, even if it's duplicative
    const static size_t SNARL_START_NODE_OFFSET = 7;
    const static size_t SNARL_END_NODE_OFFSET = 8;
    
    const static size_t CHAIN_RECORD_SIZE = 8;
    const static size_t CHAIN_NODE_COUNT_OFFSET = 1;
    const static size_t CHAIN_PARENT_OFFSET = 2;
    const static size_t CHAIN_MIN_LENGTH_OFFSET = 3;
    const static size_t CHAIN_MAX_LENGTH_OFFSET = 4;
    const static size_t CHAIN_RANK_OFFSET = 5;
    const static size_t CHAIN_START_NODE_OFFSET = 6;
    const static size_t CHAIN_END_NODE_OFFSET = 7;

    const static size_t CHAIN_NODE_RECORD_SIZE = 5; //# things for a node (not including snarl record)
    const static size_t CHAIN_NODE_PREFIX_SUM_OFFSET = 1;
    const static size_t CHAIN_NODE_FORWARD_LOOP_OFFSET = 2;
    const static size_t CHAIN_NODE_REVERSE_LOOP_OFFSET = 3;
    const static size_t CHAIN_NODE_SNARL_SIZE_OFFSET = 4;
    //This is only for multicomponent chains (in which case CHAIN_NODE_RECORD_SIZE should be 6)
    const static size_t CHAIN_NODE_COMPONENT_OFFSET = 4;
    //If there is a snarl, the snarl size is repeated after the snarl record


    //Get the offset into snarl_tree_records for a node record
    size_t get_offset_from_node_id (const id_t& id) const {
        size_t node_records_offset = snarl_tree_records.at(COMPONENT_COUNT_OFFSET) + ROOT_RECORD_SIZE; 
        size_t offset = (id-snarl_tree_records.at(MIN_NODE_ID_OFFSET)) * NODE_RECORD_SIZE;
        return node_records_offset + offset; 
    }
    //And its inverse, get the id from the offset of the node record
    //TODO: Do I want to add the min_node_id?
    id_t get_node_id_from_offset(size_t offset) const {
        size_t min_node_id = snarl_tree_records.at(MIN_NODE_ID_OFFSET);
        size_t node_records_offset = snarl_tree_records.at(COMPONENT_COUNT_OFFSET) + ROOT_RECORD_SIZE; 
        return ((offset-node_records_offset) / NODE_RECORD_SIZE) + min_node_id;
    }

    //Define a struct for interpreting each type of snarl tree node record (For node, snarl, chain)
    //These will be used with net_handle_t's, which store a pointer into snarl_tree_records 
    //as an offset. The last 4 bits of the net handle will be the connectivity of the handle,
    //as a connectivity_t
    //TODO: This might be overkill but I want it to be easy to change what gets stored in the index
    //
    struct SnarlTreeRecord {


        //The offset of the start of this record in snarl_tree_records
        size_t record_offset;
        const vector<size_t>* records;

        //Constructors assuming that this record already exists
        SnarlTreeRecord(){};
        SnarlTreeRecord (size_t pointer, const vector<size_t>* tree_records){
            record_offset = pointer;
            records = tree_records;

            record_t type = get_record_type();
            assert(type >= 1 && type <= 13 );
        }
        SnarlTreeRecord (const net_handle_t& net, const vector<size_t>* tree_records){
            record_offset = get_record_offset(net);
            records = tree_records;
            record_t type = get_record_type();
            assert(type >= 1 && type <= 13 );
        }

        virtual size_t get_offset() {
            return record_offset;
        }

        //What type of snarl tree node is this?
        //This will be the first value of any record
        virtual record_t get_record_type() const {
            return static_cast<record_t>(records->at(record_offset) >> 9);
        }

        //The name is a bit misleading, it is the handle type that the record thinks it is, 
        //not necessarily the record type of the net_handle_t that was used to produce it
        virtual net_handle_record_t get_record_handle_type() const {
            record_t type= get_record_type();
            //TODO: I"m not sure if a root snarl should be a root or a snarl
            if (type == ROOT || type == ROOT_SNARL || type == DISTANCED_ROOT_SNARL) {
                return ROOT_HANDLE;
            } else if (type == NODE || type == DISTANCED_NODE) {
                return NODE_HANDLE;
            } else if (type == SNARL || type == DISTANCED_SNARL || type ==  SIMPLE_SNARL ||type ==  OVERSIZED_SNARL){
                return SNARL_HANDLE;
            } else if (type == CHAIN || type == DISTANCED_CHAIN || type == MULTICOMPONENT_CHAIN) {
                return CHAIN_HANDLE;
            } else {
                throw runtime_error("error: trying to get the handle type of a list of children");
            }
        }

        //Get the internal connectivity of the structure
        virtual bool is_start_start_connected() const {return records->at(record_offset) & 32;}
        virtual bool is_start_end_connected() const {return records->at(record_offset) & 16;}
        virtual bool is_start_tip_connected() const {return records->at(record_offset) & 8;}
        virtual bool is_end_end_connected() const {return records->at(record_offset) & 4;}
        virtual bool is_end_tip_connected() const {return records->at(record_offset) & 2;}
        virtual bool is_tip_tip_connected() const {return records->at(record_offset) & 1;}

        //And the external connectivity. This is only relevant for root-level structures
        //since it would otherwise be captured by the containing snarl
        virtual bool is_externally_start_end_connected() const {return records->at(record_offset) & 64;}
        virtual bool is_externally_start_start_connected() const {return records->at(record_offset) & 128;}
        virtual bool is_externally_end_end_connected() const {return records->at(record_offset) & 256;}

        virtual bool has_connectivity(connectivity_t connectivity) const {
            if (connectivity == START_START) {
                return is_start_start_connected();
            } else if (connectivity == START_END || connectivity == END_START) {
                return is_start_end_connected();
            } else if (connectivity == START_TIP || connectivity == TIP_START) {
                return is_start_tip_connected();
            } else if (connectivity == END_END) {
                return is_end_end_connected();
            } else if (connectivity == END_TIP || connectivity == TIP_END) {
                return is_end_tip_connected(); 
            } else if (connectivity == TIP_TIP) {
                return is_tip_tip_connected();
            } else {
                throw runtime_error("error: Invalid connectivity");
            }
        }

        virtual bool has_connectivity(endpoint_t start, endpoint_t end){
            return has_connectivity(endpoints_to_connectivity(start, end));
        }

        virtual bool has_external_connectivity(endpoint_t start, endpoint_t end) {
            if (start == START && end == START) {
                return is_externally_start_start_connected();
            } else if (start == END && end == END) {
                return is_externally_end_end_connected();
            } else if ( (start == START && end == END ) || (start == END && end == START)) {
                return is_externally_start_end_connected();
            } else {
                return false;
            }
        }

        //Get and set a pointer to this node's parent, including its orientation
        //TODO: I don't think it matters if a chain is reversed or not, also it might not matter if a snarl is
        virtual size_t get_parent_record_offset() const {
            record_t type = get_record_type();
            if (type == ROOT) {
                return 0;
            } else if (type == NODE || type == DISTANCED_NODE) {
                return (records->at(record_offset + NODE_PARENT_OFFSET));
            } else if (type == SNARL || type == DISTANCED_SNARL || type == OVERSIZED_SNARL  
                    || type == ROOT_SNARL || type == DISTANCED_ROOT_SNARL)  {
                return (records->at(record_offset + SNARL_PARENT_OFFSET));
            } else if (type == CHAIN || type == DISTANCED_CHAIN || type == MULTICOMPONENT_CHAIN)  {
                return (records->at(record_offset + CHAIN_PARENT_OFFSET));
            } else {
                throw runtime_error("error: trying to access a snarl tree node of the wrong type");
            }
        };

        //Get and set the minimum length (distance from start to end, including boundaries for 
        // chains but not snarls, just node length for nodes)
        virtual size_t get_min_length() const {
            record_t type = get_record_type();
            if (type == DISTANCED_NODE) {
                return records->at(record_offset + NODE_LENGTH_OFFSET);
            } else if (type == DISTANCED_SNARL || type == OVERSIZED_SNARL)  {
                return records->at(record_offset + SNARL_MIN_LENGTH_OFFSET);
            } else if (type == DISTANCED_CHAIN || type == MULTICOMPONENT_CHAIN)  {
                return records->at(record_offset + CHAIN_MIN_LENGTH_OFFSET);
            } else if (type == NODE || type == SNARL || type == CHAIN) {
                throw runtime_error("error: trying to access get distance in a distanceless index");
            } else if (type == ROOT || type == ROOT_SNARL || type == DISTANCED_ROOT_SNARL) {
                throw runtime_error("error: trying to find the length of the root");
            } else {
                throw runtime_error("error: trying to access a snarl tree node of the wrong type");
            }
        };
        
        //Get and set this node's maximum length
        //This isn't actually a maximum, it's the maximum among minimum distance paths 
        //through each node in the snarl/chain
        virtual size_t get_max_length() const {
            record_t type = get_record_type();
            if (type == DISTANCED_NODE) {
                return records->at(record_offset + NODE_LENGTH_OFFSET);
            } else if (type == DISTANCED_SNARL || type == OVERSIZED_SNARL)  {
                return records->at(record_offset + SNARL_MAX_LENGTH_OFFSET);
            } else if (type == DISTANCED_CHAIN || type == MULTICOMPONENT_CHAIN)  {
                return records->at(record_offset + CHAIN_MAX_LENGTH_OFFSET);
            } else if (type == NODE || type == SNARL || type == CHAIN) {
                throw runtime_error("error: trying to access get distance in a distanceless index");
            }  else if (type == ROOT || type == ROOT_SNARL || type == DISTANCED_ROOT_SNARL) {
                throw runtime_error("error: trying to find the length of the root");
            } else {
                throw runtime_error("error: trying to access a snarl tree node of the wrong type");
            }
        };

        //Get and set this structure's rank in its parent
        //For children of snarls, this means the actual rank
        //For children of chains, it points to the node in the chain
        virtual size_t get_rank_in_parent() const {
            record_t type = get_record_type();
            if (type == NODE || type == DISTANCED_NODE) {
                return records->at(record_offset + NODE_RANK_OFFSET) >> 1;
            } else if (type == SNARL || type == DISTANCED_SNARL || type == OVERSIZED_SNARL
                    || type == ROOT_SNARL || type == DISTANCED_ROOT_SNARL)  {
                return records->at(record_offset + SNARL_RANK_OFFSET) >> 1;
            } else if (type == CHAIN || type == DISTANCED_CHAIN || type == MULTICOMPONENT_CHAIN)  {
                return records->at(record_offset + CHAIN_RANK_OFFSET) >> 1;
            } else {
                throw runtime_error("error: trying to access a snarl tree node of the wrong type");
            }
        };

        //Is this node reversed in its parent
        //TODO: I think a node is the only thing that can be reversed in its parent, and only if it is in a chain
        virtual bool get_is_reversed_in_parent() const {
            record_t type = get_record_type();
            if (type == NODE || type == DISTANCED_NODE) {
                return records->at(record_offset + NODE_RANK_OFFSET) & 1;
            } else if (type == SNARL || type == DISTANCED_SNARL || type == OVERSIZED_SNARL 
                    || type == ROOT_SNARL || type == DISTANCED_ROOT_SNARL)  {
                return records->at(record_offset + SNARL_RANK_OFFSET) & 1;
            } else if (type == CHAIN || type == DISTANCED_CHAIN || type == MULTICOMPONENT_CHAIN)  {
                return records->at(record_offset + CHAIN_RANK_OFFSET) & 1;
            } else {
                throw runtime_error("error: trying to access a snarl tree node of the wrong type");
            }
        };

        //Get the node id of the start/end of this structure (start node of a snarl/chain)
        //TODO: Also need to add min node id
        virtual id_t get_start_id() const {
            record_t type = get_record_type();
            if (type == ROOT) {
                //TODO: Also not totally sure what this should do
                throw runtime_error("error: trying to get the start node of the root");
            } else if (type == NODE || type == DISTANCED_NODE) {
                //TODO: I'm not sure if I want to allow this
                cerr << "warning: Looking for the start of a node" << endl;
                return get_id_from_offset(record_offset);
            } else if (type == SNARL || type == DISTANCED_SNARL || type == OVERSIZED_SNARL)  {
                return (records->at(record_offset + SNARL_START_NODE_OFFSET)) >> 1;
            } else if (type == CHAIN || type == DISTANCED_CHAIN || type == MULTICOMPONENT_CHAIN)  {
                return (records->at(record_offset + CHAIN_START_NODE_OFFSET)) >> 1;
            } else if (type == ROOT_SNARL || type == DISTANCED_ROOT_SNARL) {
                throw runtime_error("error: trying to find the start node of a root snarl");
            } else {
                throw runtime_error("error: trying to access a snarl tree node of the wrong type");
            }
        }
        //True if the node is traversed backwards to enter the structure
        virtual bool get_start_orientation() const {
            record_t type = get_record_type();
            if (type == ROOT) {
                //TODO: Also not totally sure what this should do
                throw runtime_error("error: trying to get the start node of the root");
            } else if (type == NODE || type == DISTANCED_NODE) {
                //TODO: I'm not sure if I want to allow this
                cerr << "warning: Looking for the start of a node" << endl;
                return false;
            } else if (type == SNARL || type == DISTANCED_SNARL || type == OVERSIZED_SNARL)  {
                return (records->at(record_offset + SNARL_START_NODE_OFFSET) & 1);
            } else if (type == CHAIN || type == DISTANCED_CHAIN || type == MULTICOMPONENT_CHAIN)  {
                return (records->at(record_offset + CHAIN_START_NODE_OFFSET)) & 1;
            } else if (type == ROOT_SNARL || type == DISTANCED_ROOT_SNARL) {
                throw runtime_error("error: trying to find the start node of a root snarl");
            }else {
                throw runtime_error("error: trying to access a snarl tree node of the wrong type");
            }
        }
        virtual id_t get_end_id() const {
            record_t type = get_record_type();
            if (type == ROOT) {
                //TODO: Also not totally sure what this should do
                throw runtime_error("error: trying to get the end node of the root");
            } else if (type == NODE || type == DISTANCED_NODE) {
                //TODO: I'm not sure if I want to allow this
                cerr << "warning: Looking for the end of a node" << endl;
                //TODO: Put this in its own function? Also double check for off by ones
                //Offset of the start of the node vector
                return get_id_from_offset(record_offset);
            } else if (type == SNARL || type == DISTANCED_SNARL || type == OVERSIZED_SNARL)  {
                return (records->at(record_offset + SNARL_END_NODE_OFFSET)) >> 1;
            } else if (type == CHAIN || type == DISTANCED_CHAIN || type == MULTICOMPONENT_CHAIN)  {
                return (records->at(record_offset + CHAIN_END_NODE_OFFSET)) >> 1;
            } else if (type == ROOT_SNARL || type == DISTANCED_ROOT_SNARL) {
                throw runtime_error("error: trying to find the end node of a root snarl");
            } else {
                throw runtime_error("error: trying to access a snarl tree node of the wrong type");
            }
        }
        //Return true if the end node is traversed backwards to leave the snarl
        virtual id_t get_end_orientation() const {
            record_t type = get_record_type();
            if (type == ROOT) {
                //TODO: Also not totally sure what this should do
                throw runtime_error("error: trying to get the end node of the root");
            } else if (type == NODE || type == DISTANCED_NODE) {
                //TODO: I'm not sure if I want to allow this
                cerr << "warning: Looking for the end of a node" << endl;
                //TODO: Put this in its own function? Also double check for off by ones
                //Offset of the start of the node vector
                return get_id_from_offset(record_offset);
            } else if (type == SNARL || type == DISTANCED_SNARL || type == OVERSIZED_SNARL)  {
                return (records->at(record_offset + SNARL_END_NODE_OFFSET)) & 1;
            } else if (type == CHAIN || type == DISTANCED_CHAIN || type == MULTICOMPONENT_CHAIN)  {
                return (records->at(record_offset + CHAIN_END_NODE_OFFSET)) & 1;
            } else if (type == ROOT_SNARL || type == DISTANCED_ROOT_SNARL) {
                throw runtime_error("error: trying to find the end node of a root snarl");
            } else {
                throw runtime_error("error: trying to access a snarl tree node of the wrong type");
            }
        }

        //TODO: These are redeclared so that I don't need to pass the SnarlTreeRecord the actual distance index
        //Get the offset into snarl_tree_records for a node record
        virtual size_t get_offset_from_id (const id_t id) const {
            size_t node_records_offset = records->at(COMPONENT_COUNT_OFFSET) + ROOT_RECORD_SIZE; 
            size_t offset = (id-records->at(MIN_NODE_ID_OFFSET)) * NODE_RECORD_SIZE;
            return node_records_offset + offset; 

        }
        //And its inverse, get the id from the offset of the node record
        //TODO: Do I want to add the min_node_id?
        virtual id_t get_id_from_offset(size_t offset) const {
        size_t min_node_id = records->at(MIN_NODE_ID_OFFSET);
        size_t node_records_offset = records->at(COMPONENT_COUNT_OFFSET) + ROOT_RECORD_SIZE; 
            return ((offset-node_records_offset) / NODE_RECORD_SIZE) + min_node_id;
        }
    };

    //Record interpreter that has a non-const reference to snarl_tree_records, so it can
    //also add things
    struct SnarlTreeRecordConstructor {

        //TODO: This might be bad but this has the same members as SnarlTreeRecord
        //and does basically the same thing but doesn't inherit from it
        size_t record_offset;
        vector<size_t>* records;

        //Constructors assuming that this record already exists
        SnarlTreeRecordConstructor() {};
        SnarlTreeRecordConstructor (size_t pointer, vector<size_t>* tree_records){
            record_offset = pointer;
            records = tree_records;

            record_t type = get_record_type();
            assert(type == ROOT || type == NODE || type == DISTANCED_NODE || type == SNARL || 
                    type == DISTANCED_SNARL || type == OVERSIZED_SNARL || 
                    type == ROOT_SNARL || type == DISTANCED_ROOT_SNARL || type == CHAIN || 
                    type == DISTANCED_CHAIN || type == MULTICOMPONENT_CHAIN);
        }
        SnarlTreeRecordConstructor (const net_handle_t& net, vector<size_t>* tree_records){
            record_offset = get_record_offset(net);
            records = tree_records;
            record_t type = get_record_type();
            assert(type == ROOT || type == NODE || type == DISTANCED_NODE || type == SNARL || 
                    type == DISTANCED_SNARL || type == OVERSIZED_SNARL || 
                    type == ROOT_SNARL || type == DISTANCED_ROOT_SNARL || type == CHAIN || 
                    type == DISTANCED_CHAIN || type == MULTICOMPONENT_CHAIN);
        }
        //What type of snarl tree node is this?
        //This will be the first value of any record
        //TODO: This copies from SnarlTreeRecord
        virtual record_t get_record_type() const {
            return static_cast<record_t>(records->at(record_offset) >> 9);
        }
        virtual void set_start_start_connected() {
#ifdef debug_indexing
            cerr << record_offset << " set start_start_connected" << endl;
#endif
            records->at(record_offset) = records->at(record_offset) | 32;
        }
        virtual void set_start_end_connected() {
#ifdef debug_indexing
            cerr << record_offset << " set start_end_connected" << endl;
#endif
            records->at(record_offset) = records->at(record_offset) | 16;
        }
        virtual void set_start_tip_connected() {
#ifdef debug_indexing
            cerr << record_offset << " set start_tip_connected" << endl;
#endif
            records->at(record_offset) = records->at(record_offset) | 8;
        }
        virtual void set_end_end_connected() {
#ifdef debug_indexing
            cerr << record_offset << " set end_end_connected" << endl;
#endif
            records->at(record_offset) = records->at(record_offset) | 4;
        }
        virtual void set_end_tip_connected() {
#ifdef debug_indexing
            cerr << record_offset << " set end_tip_connected" << endl;
#endif
            records->at(record_offset) = records->at(record_offset) | 2;
        }
        virtual void set_tip_tip_connected() {
#ifdef debug_indexing
            cerr << record_offset << " set tpi_tip_connected" << endl;
#endif
            records->at(record_offset) = records->at(record_offset) | 1;
        }
        virtual void set_externally_start_end_connected() {
#ifdef debug_indexing
            cerr << record_offset << " set externally start_end connected" << endl;
#endif
            records->at(record_offset) = records->at(record_offset) | 64;
        }
        virtual void set_externally_start_start_connected() const {
#ifdef debug_indexing
            cerr << record_offset << " set externally start_start connected" << endl;
#endif
            records->at(record_offset) = records->at(record_offset) | 128;
        }
        virtual void set_externally_end_end_connected() const {
#ifdef debug_indexing
            cerr << record_offset << " set externally end_end connected" << endl;
#endif
            records->at(record_offset) = records->at(record_offset) | 256;
        }
        virtual void set_record_type(record_t type) {
#ifdef debug_indexing
            cerr << record_offset << " set record type to be " << type << endl;
#endif
            assert(records->at(record_offset) == 0);
            records->at(record_offset) = ((static_cast<size_t>(type) << 9) | (records->at(record_offset) & 511));
        }

        virtual void set_min_length(size_t length) {
            record_t type = get_record_type();
            size_t offset;
            if (type == DISTANCED_NODE) {
                offset = record_offset + NODE_LENGTH_OFFSET;
            } else if (type == DISTANCED_SNARL || type == OVERSIZED_SNARL)  {
                offset = record_offset + SNARL_MIN_LENGTH_OFFSET;
            } else if (type == DISTANCED_CHAIN || type == MULTICOMPONENT_CHAIN)  {
                offset = record_offset + CHAIN_MIN_LENGTH_OFFSET;
            } else if (type == NODE || type == SNARL || type == CHAIN ) {
                throw runtime_error("error: trying to access get distance in a distanceless index");
            } else if (type == ROOT_SNARL || type == DISTANCED_ROOT_SNARL) {
                throw runtime_error("error: set the length of a root snarl");
            } else {
                throw runtime_error("error: trying to access a snarl tree node of the wrong type");
            }
#ifdef debug_indexing
            cerr << offset << " set min length to be " << length << endl;
            assert(records->at(offset) == 0);
#endif
            records->at(offset) = length;
        };
        virtual void set_max_length(size_t length) {
            record_t type = get_record_type();
            size_t offset;
            if (type == DISTANCED_NODE) {
                offset = record_offset + NODE_LENGTH_OFFSET;
            } else if (type == DISTANCED_SNARL || type == OVERSIZED_SNARL)  {
                offset = record_offset + SNARL_MAX_LENGTH_OFFSET;
            } else if (type == DISTANCED_CHAIN || type == MULTICOMPONENT_CHAIN)  {
                offset = record_offset + CHAIN_MAX_LENGTH_OFFSET;
            } else if (type == DISTANCED_NODE || type == SNARL || type == CHAIN) {
                throw runtime_error("error: trying to access get distance in a distanceless index");
            }  else if (type == ROOT_SNARL || type == DISTANCED_ROOT_SNARL) {
                throw runtime_error("error: set the length of a root snarl");
            } else {
                throw runtime_error("error: trying to access a snarl tree node of the wrong type");
            }
#ifdef debug_indexing
            cerr << offset << " set max length to be " << length << endl;
            assert(records->at(offset) == 0);
#endif
            records->at(offset) = length;
        };
        virtual void set_rank_in_parent(size_t rank) {
            record_t type = get_record_type();
            size_t offset;
            if (type == NODE || type == DISTANCED_NODE) {
                offset = record_offset + NODE_RANK_OFFSET;
            } else if (type == SNARL || type == DISTANCED_SNARL || type == OVERSIZED_SNARL
                    || type == ROOT_SNARL || type == DISTANCED_ROOT_SNARL)  {
                offset = record_offset + SNARL_RANK_OFFSET;
            } else if (type == CHAIN || type == DISTANCED_CHAIN || type == MULTICOMPONENT_CHAIN)  {
                offset = record_offset + CHAIN_RANK_OFFSET;
            } else {
                throw runtime_error("error: trying to access a snarl tree node of the wrong type");
            }
#ifdef debug_indexing
            cerr << offset << " set rank in parent to be " << rank << endl;
            assert(records->at(offset) >> 1 == 0);
#endif
            bool rev = records->at(offset) & 1;
            records->at(offset) = (rank << 1) | rev; 
        };
        virtual void set_is_reversed_in_parent(bool rev) {
            record_t type = get_record_type();
            size_t offset;
            if (type == NODE || type == DISTANCED_NODE) {
                offset = record_offset + NODE_RANK_OFFSET;
            } else if (type == SNARL || type == DISTANCED_SNARL || type == OVERSIZED_SNARL 
                    || type == ROOT_SNARL || type == DISTANCED_ROOT_SNARL)  {
                offset = record_offset + SNARL_RANK_OFFSET;
            } else if (type == CHAIN || type == DISTANCED_CHAIN || type == MULTICOMPONENT_CHAIN)  {
                offset = record_offset + CHAIN_RANK_OFFSET;
            } else {
                throw runtime_error("error: trying to access a snarl tree node of the wrong type");
            }
#ifdef debug_indexing
            cerr << offset << " set rev in parent to be " << rev << endl;
#endif
            records->at(offset) =  ((records->at(offset)>>1)<<1) | rev; 
        };
        virtual void set_parent_record_offset(size_t pointer){
            record_t type = get_record_type();
            size_t offset; 
            if (type == NODE || type == DISTANCED_NODE) {
                offset = record_offset + NODE_PARENT_OFFSET;
            } else if (type == SNARL || type == DISTANCED_SNARL || type == OVERSIZED_SNARL
                    || type == ROOT_SNARL || type == DISTANCED_ROOT_SNARL)  {
#ifdef debug_indexing
                if (type == ROOT_SNARL || type == DISTANCED_ROOT_SNARL) {
                    assert(pointer == 0);
                }
#endif
                offset = record_offset + SNARL_PARENT_OFFSET;
            } else if (type == CHAIN || type == DISTANCED_CHAIN || type == MULTICOMPONENT_CHAIN)  {
                offset = record_offset + CHAIN_PARENT_OFFSET;
            } else {
                throw runtime_error("error: trying to access a snarl tree node of the wrong type");
            }
#ifdef debug_indexing
            cerr << offset << " set parent offset to be " << pointer << endl;
            assert(records->at(offset) == 0);
#endif
            records->at(offset) = pointer;
        };
        //Rev is true if the node is traversed backwards to enter the snarl
        virtual void set_start_node(id_t id, bool rev) {
            record_t type = get_record_type();
            size_t offset;
            if (type == ROOT || type == NODE || type == DISTANCED_NODE) {
                throw runtime_error("error: trying to set the node id of a node or root");
            } else if (type == SNARL || type == DISTANCED_SNARL || type == OVERSIZED_SNARL)  {
                offset = record_offset + SNARL_START_NODE_OFFSET;
            } else if (type == CHAIN || type == DISTANCED_CHAIN || type == MULTICOMPONENT_CHAIN)  {
                offset = record_offset + CHAIN_START_NODE_OFFSET;
            }  else if (type == ROOT_SNARL || type == DISTANCED_ROOT_SNARL) {
                throw runtime_error("error: set the start node of a root snarl");
            } else {
                throw runtime_error("error: trying to access a snarl tree node of the wrong type");
            }
#ifdef debug_indexing
            cerr << offset << " set start node to be " << id << " facing " << (rev ? "rev" : "fd") << endl;
            assert(records->at(offset) == 0);
#endif
            records->at(offset) = (id << 1) | rev;
        }
        //Rev is true if the node is traversed backwards to leave the snarl
        virtual void set_end_node(id_t id, bool rev) const {
            record_t type = get_record_type();
            size_t offset;
            if (type == ROOT || type == NODE || type == DISTANCED_NODE) {
                throw runtime_error("error: trying to set the node id of a node or root");
            } else if (type == SNARL || type == DISTANCED_SNARL || type == OVERSIZED_SNARL)  {
                offset = record_offset + SNARL_END_NODE_OFFSET;
            } else if (type == CHAIN || type == DISTANCED_CHAIN || type == MULTICOMPONENT_CHAIN)  {
                offset = record_offset + CHAIN_END_NODE_OFFSET;
            }  else if (type == ROOT_SNARL || type == DISTANCED_ROOT_SNARL) {
                throw runtime_error("error: set the end node of a root snarl");
            } else {
                throw runtime_error("error: trying to access a snarl tree node of the wrong type");
            }
#ifdef debug_indexing
            cerr << offset << " set end node to be " << id << " facing " << (rev ? "rev" : "fd") << endl;
            assert(records->at(offset) == 0);
#endif
            records->at(offset) = (id << 1) | rev;
        }
    };

    struct RootRecord : SnarlTreeRecord {

        RootRecord (){};
        RootRecord (size_t pointer, const vector<size_t>* tree_records){
            record_offset = pointer;
            records = tree_records;
            assert(get_record_type() == ROOT);
        }
        RootRecord (net_handle_t net, const vector<size_t>* tree_records){
            record_offset = get_record_offset(net);
            records = tree_records;
            assert(get_record_type() == ROOT);
        }
        virtual size_t get_connected_component_count() const {
            return records->at(record_offset+COMPONENT_COUNT_OFFSET);
        }
        virtual size_t get_node_count() const {
            return records->at(record_offset+NODE_COUNT_OFFSET);
        }
        virtual size_t get_min_node_id() const {
            return records->at(record_offset+MIN_NODE_ID_OFFSET);
        }

        virtual SnarlTreeRecord get_component_record(size_t component_number) const {
            //TODO: Maybe should be +1 if the component numbers start at 1?
            return SnarlTreeRecord(records->at(record_offset+2+component_number), records);
        }
        virtual bool for_each_child(const std::function<bool(const handlegraph::net_handle_t&)>& iteratee) const {
            size_t connected_component_count = get_connected_component_count();
            for (size_t i = 0 ; i < connected_component_count ; i++) {

                size_t child_offset = records->at(record_offset + ROOT_RECORD_SIZE + i);
                net_handle_record_t type = SnarlTreeRecord(child_offset, records).get_record_handle_type();
                record_t record_type = SnarlTreeRecord(child_offset, records).get_record_type();


                if (record_type == ROOT_SNARL || record_type == DISTANCED_ROOT_SNARL) {
                    //This is a bunch of root components that are connected, so go through each 
                    SnarlRecord snarl_record(child_offset, records);
                    if (! snarl_record.for_each_child(iteratee)) {
                        return false;
                    }
                } else { 
                    //Otherwise, it is a separate connected component
                    net_handle_t child_handle =  get_net_handle(child_offset, START_END, type);
                    if (!iteratee(child_handle)) {
                        return false;
                    }
                }
            }
            return true;
        }

    };

    struct RootRecordConstructor : RootRecord, SnarlTreeRecordConstructor {

        //Constructor meant for creating a new record, at the end of snarl_tree_records
        //TODO: The way I wrote this pointer should be 0
        RootRecordConstructor (size_t pointer, size_t connected_component_count, size_t node_count, 
                    id_t min_node_id, vector<size_t>* records){
            RootRecord::record_offset = pointer;
            SnarlTreeRecordConstructor::record_offset = pointer;
            SnarlTreeRecordConstructor::records = records;
            RootRecord::records = records;
            //Allocate memory for the root vector and for all of the nodes
#ifdef debug_indexing
            cerr << " Resizing array to add root: length " << SnarlTreeRecordConstructor::records->size() << " -> " << SnarlTreeRecordConstructor::records->size() + ROOT_RECORD_SIZE + connected_component_count + (NODE_RECORD_SIZE * node_count) << endl;
#endif
            SnarlTreeRecordConstructor::records->resize(SnarlTreeRecordConstructor::records->size() + ROOT_RECORD_SIZE + connected_component_count + (NODE_RECORD_SIZE * node_count) , 0);
            set_record_type(ROOT);
            set_min_node_id(min_node_id);
            set_node_count(node_count);
            set_connected_component_count(connected_component_count);
        }
        virtual void set_connected_component_count(size_t connected_component_count) {
#ifdef debug_indexing
            cerr << RootRecord::record_offset+COMPONENT_COUNT_OFFSET << " set connected component to be " << connected_component_count << endl;
            assert(SnarlTreeRecordConstructor::records->at(RootRecord::record_offset+COMPONENT_COUNT_OFFSET) == 0);
#endif
            SnarlTreeRecordConstructor::records->at(RootRecord::record_offset+COMPONENT_COUNT_OFFSET)=connected_component_count;
        }
        virtual void set_node_count(size_t node_count) {
#ifdef debug_indexing
            cerr << RootRecord::record_offset+NODE_COUNT_OFFSET << " set node count to be " << node_count << endl;
            assert(RootRecord::records->at(RootRecord::record_offset+NODE_COUNT_OFFSET) == 0);
#endif
            SnarlTreeRecordConstructor::records->at(RootRecord::record_offset+NODE_COUNT_OFFSET)=node_count;
        }
        virtual void set_min_node_id(id_t node_id) {
#ifdef debug_indexing
            cerr << RootRecord::record_offset+MIN_NODE_ID_OFFSET << " set min node id to be " << node_id << endl;
            assert(RootRecord::records->at(RootRecord::record_offset+MIN_NODE_ID_OFFSET) == 0);
#endif
            SnarlTreeRecordConstructor::records->at(RootRecord::record_offset+MIN_NODE_ID_OFFSET)=node_id;
        }
        virtual void add_component(size_t index, size_t offset) {
#ifdef debug_indexing
            cerr << RootRecord::record_offset+ROOT_RECORD_SIZE+index << " set new component " << offset << endl;
            assert(SnarlTreeRecordConstructor::records->at(RootRecord::record_offset+ROOT_RECORD_SIZE+index) == 0);
            assert(index < get_connected_component_count());
#endif
            SnarlTreeRecordConstructor::records->at(RootRecord::record_offset+ROOT_RECORD_SIZE+index)
                = offset;

        }
    };
    struct NodeRecord : SnarlTreeRecord {


        NodeRecord() {};
        NodeRecord (size_t pointer, const vector<size_t>* tree_records){
            record_offset = pointer;
            records = tree_records;

            assert(get_record_type() == NODE || get_record_type() == DISTANCED_NODE);
        }
        NodeRecord (net_handle_t net, const vector<size_t>* tree_records){
            records = tree_records;
            record_offset = get_record_offset(net);

            assert(get_handle_type(net) == NODE_HANDLE);
            assert(get_record_type() == NODE || get_record_type() == DISTANCED_NODE);
            assert(get_connectivity(net) == START_END || get_connectivity(net) == END_START);
        }

        virtual id_t get_node_id() const {
            return get_id_from_offset(record_offset);
        }

        //TODO: This one is a bit redundant but fine I think
        virtual int64_t get_node_length() const {
            return records->at(record_offset + NODE_LENGTH_OFFSET);
        }

        virtual size_t get_root_component() const {
            return records->at(record_offset + NODE_COMPONENT_OFFSET);
        }

        virtual bool in_chain() const {
            return SnarlTreeRecord(get_parent_record_offset(), records).get_record_handle_type() == CHAIN_HANDLE;
        }

    };

    struct NodeRecordConstructor : NodeRecord , SnarlTreeRecordConstructor {


        //Constructor meant for creating a new record, at the end of snarl_tree_records
        //The memory for all nodes has already been allocated by the root
        NodeRecordConstructor (size_t pointer, record_t type, vector<size_t>* records){
            SnarlTreeRecordConstructor::records = records;
            NodeRecord::records = records;
            //TODO: Only one get_offset_from_id?
            NodeRecord::record_offset =  pointer;
            SnarlTreeRecordConstructor::record_offset =  pointer;
            set_record_type(type);
            set_start_end_connected();
        }
        virtual void set_node_length(size_t length) {
#ifdef debug_indexing
            cerr <<NodeRecord::record_offset + NODE_LENGTH_OFFSET << " set node length " << length << endl;
            assert(SnarlTreeRecordConstructor::records->at(NodeRecord::record_offset + NODE_LENGTH_OFFSET) == 0);
#endif
            SnarlTreeRecordConstructor::records->at(NodeRecord::record_offset + NODE_LENGTH_OFFSET) = length;
        }
        virtual void set_root_component(size_t component) {
#ifdef debug_indexing
            cerr <<NodeRecord::record_offset + NODE_COMPONENT_OFFSET << " set root component " << component << endl;
            assert(SnarlTreeRecordConstructor::records->at(NodeRecord::record_offset + NODE_COMPONENT_OFFSET) == 0);
#endif
            SnarlTreeRecordConstructor::records->at(NodeRecord::record_offset + NODE_COMPONENT_OFFSET) = component;
        }
    };

    struct SnarlRecord : SnarlTreeRecord {


        SnarlRecord() {};
        SnarlRecord (size_t pointer, const vector<size_t>* tree_records){
            record_offset = pointer;
            records = tree_records;
            record_t type = get_record_type();
            assert(type == SNARL || type == DISTANCED_SNARL || type == OVERSIZED_SNARL || type == ROOT_SNARL
                || type == ROOT_SNARL || type == DISTANCED_ROOT_SNARL);
        }

        SnarlRecord (net_handle_t net, const vector<size_t>* tree_records){
            record_offset = get_record_offset(net);
            records = tree_records;
            net_handle_record_t type = get_handle_type(net);
            assert(type == SNARL_HANDLE || type == SENTINEL_HANDLE || type == ROOT_HANDLE);
        }

        //How big is the entire snarl record?
        static size_t distance_vector_size(record_t type, size_t node_count) {
            if (type == SNARL || type == ROOT_SNARL){
                //For a normal snarl, its just the record size and the pointers to children
                return 0; 
            } else if (type == DISTANCED_SNARL) {
                //For a normal min distance snarl, record size and the pointers to children, and
                //matrix of distances
                size_t node_side_count = node_count * 2 + 2;
                return  (((node_side_count+1)*node_side_count) / 2);
            } else if (type ==  OVERSIZED_SNARL){
                //For a large min_distance snarl, record the side, pointers to children, and just 
                //the min distances from each node side to the two boundary nodes
                size_t node_side_count = node_count * 2 + 2;
                return  (node_side_count * 2);
            } else if (type == DISTANCED_ROOT_SNARL) {
                size_t node_side_count = node_count * 2;
                return  (((node_side_count+1)*node_side_count) / 2);
            } else {
                throw runtime_error ("error: this is not a snarl");
            }
        }
        static size_t record_size (record_t type, size_t node_count) {
            return SNARL_RECORD_SIZE + distance_vector_size(type, node_count) + node_count;
        }
        virtual size_t record_size() { 
            record_t type = get_record_type();
            return record_size(type, get_node_count()); 
        }


        //Get the index into the distance vector for the calculating distance between the given node sides
        static int64_t get_distance_vector_offset(size_t rank1, bool right_side1, size_t rank2, 
                bool right_side2, size_t node_count, record_t type) {

            //how many node sides in this snarl
            size_t node_side_count = type == DISTANCED_ROOT_SNARL ? node_count * 2 : node_count * 2 + 2; 

            //make sure we're looking at the correct node side
            //If this is the start or end node, then we don't adjust based on orientation
            //because we only care about the inner sides. If this is a root snarl, then
            //there is no start or end node and the ranks 0 and 1 are not special
             if (type == DISTANCED_ROOT_SNARL) { 
                rank1 = rank1 * 2;
                if (right_side1) {
                    rank1 += 1;
                }
             } else if (rank1 != 0 && rank1 != 1) {
                rank1 = (rank1-1) * 2;
                if (right_side1) {
                    rank1 += 1;
                }
            }
            if (type == DISTANCED_ROOT_SNARL) {
                rank2 = rank2 * 2;
                if (right_side2) {
                    rank2 += 1;
                }
            } else if (rank2 != 0 && rank2 != 1) {
                rank2 = (rank2-1) * 2;
                if (right_side2) {
                    rank2 += 1;
                }
            }

            //reverse order of ranks if necessary
            if (rank1 > rank2) {
                size_t tmp = rank1;
                rank1 = rank2;
                rank2 = tmp;
            }
            
            if (type == SNARL || type == ROOT_SNARL) {
                throw runtime_error("error: trying to access distance in a distanceless snarl tree");
            } else if (type == DISTANCED_SNARL || type == DISTANCED_ROOT_SNARL) {
                //normal distance index
                size_t k = node_side_count-rank1;
                return (((node_side_count+1) * node_side_count)/2) - (((k+1)*k) / 2) + rank2 - rank1;
            } else if (type ==  OVERSIZED_SNARL) {
                //abbreviated distance index storing only the distance from each node side to the
                //start and end
                if (rank1 == 0) {
                    return rank2;
                } else if (rank1 == 1) {
                    return node_count + 2 + rank2;
                } else {
                    throw runtime_error("error: trying to access distance in an oversized snarl");
                }
            } else {
                throw runtime_error("error: trying to distance from something that isn't a snarl");
            }
        } 
        virtual int64_t get_distance_vector_offset(size_t rank1, bool right_side1, 
                size_t rank2, bool right_side2) const {
            return get_distance_vector_offset(rank1, right_side1, rank2, right_side2, 
                get_node_count(), get_record_type()); 
        }

        //TODO: I want to also add a function to do this given a node id or net_handle_t instead of rank
        //Get and set the distances between two node sides in the graph
        //Ranks identify which node, sides indicate node side: false for left, true for right
        virtual int64_t get_distance(size_t rank1, bool right_side1, size_t rank2, bool right_side2) const {

            //offset of the start of the distance vectors in snarl_tree_records
            size_t distance_vector_start = record_offset + SNARL_RECORD_SIZE + get_node_count();
            //Offset of this particular distance in the distance vector
            size_t distance_vector_offset = get_distance_vector_offset(rank1, right_side1, rank2, right_side2);

            size_t val = records->at(distance_vector_start + distance_vector_offset);
            return  val == 0 ? std::numeric_limits<int64_t>::max() : val-1;

        }

        virtual size_t get_node_count() const {
            return records->at(record_offset + SNARL_NODE_COUNT_OFFSET);
        }

        virtual size_t get_child_record_pointer() const {
            return records->at(record_offset+SNARL_CHILD_RECORD_OFFSET) ;
        }

        virtual bool for_each_child(const std::function<bool(const net_handle_t&)>& iteratee) const {
            size_t child_count = get_node_count();
            size_t child_record_offset = get_child_record_pointer();
            for (size_t i = 0 ; i < child_count ; i++) {
                size_t child_offset =  records->at(child_record_offset + i);
                net_handle_record_t type = SnarlTreeRecord(child_offset, records).get_record_handle_type(); 
                assert(type == NODE_HANDLE || type == CHAIN_HANDLE);
                net_handle_t child_handle =  get_net_handle (child_offset, START_END, CHAIN_HANDLE);
                bool result = iteratee(child_handle); 
                if (result == false) {
                    return false;
                }
            }
            return true;
        }

    };

    struct SnarlRecordConstructor : SnarlRecord , SnarlTreeRecordConstructor {


        SnarlRecordConstructor();

        SnarlRecordConstructor (size_t node_count, vector<size_t>* records, record_t type){
            //Constructor for making a new record, including allocating memory.
            //Assumes that this is the latest record being made, so pointer will be the end of
            //the array and we need to allocate extra memory past it
            //TODO:I'm not sure yet how memory will actually be allocated
            SnarlTreeRecordConstructor::records = records;
            SnarlRecord::records = records;

            SnarlRecord::record_offset = records->size();
            SnarlTreeRecordConstructor::record_offset = records->size();

            size_t extra_size = record_size(type, node_count);
#ifdef debug_indexing
            cerr << " Resizing array to add snarl: length " << SnarlTreeRecordConstructor::records->size() << " -> "  << SnarlTreeRecordConstructor::records->size() + extra_size << endl;
#endif
            SnarlTreeRecordConstructor::records->resize(SnarlTreeRecordConstructor::records->size() + extra_size, 0);
            set_node_count(node_count);
            set_record_type(type);
        }
        SnarlRecordConstructor(vector<size_t>* records, size_t pointer) {
            //Make a constructor for a snarl record that has already been allocated.
            //For adding children to an existing snarl record
            SnarlRecord::record_offset = pointer;
            SnarlTreeRecordConstructor::record_offset = pointer;
            SnarlTreeRecordConstructor::records = records;
            SnarlRecord::records = records;
        }
        void set_distance(size_t rank1, bool right_side1, size_t rank2, bool right_side2, int64_t distance) {
            //offset of the start of the distance vectors in snarl_tree_records
            size_t distance_vector_start = SnarlRecord::record_offset + SNARL_RECORD_SIZE + get_node_count();
            //Offset of this particular distance in the distance vector
            size_t distance_vector_offset = get_distance_vector_offset(rank1, right_side1, rank2, right_side2);

            size_t val = distance == std::numeric_limits<int64_t>::max() ? 0 : distance+1;
#ifdef debug_indexing
            cerr <<  distance_vector_start + distance_vector_offset << " set distance_value " << val << endl;
            //assert(SnarlTreeRecordConstructor::records->at(distance_vector_start + distance_vector_offset) == 0 || 
            //        SnarlTreeRecordConstructor::records->at(distance_vector_start + distance_vector_offset) == val);
#endif

            SnarlTreeRecordConstructor::records->at(distance_vector_start + distance_vector_offset) = val;
        }

        //Node count is the number of nodes in the snarl, not including boundary nodes
        void set_node_count(size_t node_count) {
#ifdef debug_indexing
            cerr << SnarlRecord::record_offset + SNARL_NODE_COUNT_OFFSET << " set snarl node count " << node_count << endl;
            assert(node_count > 0);//TODO: Don't bother making a record for trivial snarls
            assert(SnarlTreeRecordConstructor::records->at(SnarlRecord::record_offset + SNARL_NODE_COUNT_OFFSET) == 0);
#endif
            SnarlTreeRecordConstructor::records->at(SnarlRecord::record_offset + SNARL_NODE_COUNT_OFFSET) = node_count;
        }


        void set_child_record_pointer(size_t pointer) {
#ifdef debug_indexing
            cerr << SnarlRecord::record_offset+SNARL_CHILD_RECORD_OFFSET << " set snarl child record pointer " << pointer << endl;
            assert(SnarlTreeRecordConstructor::SnarlTreeRecordConstructor::records->at(SnarlRecord::record_offset+SNARL_CHILD_RECORD_OFFSET) == 0);
#endif
            SnarlTreeRecordConstructor::records->at(SnarlRecord::record_offset+SNARL_CHILD_RECORD_OFFSET) = pointer;
        }
        //Add a reference to a child of this snarl. Assumes that the index is completed up
        //to here
        void add_child(size_t pointer){
#ifdef debug_indexing
            cerr << SnarlTreeRecordConstructor::records->size() << " Adding child pointer to the end of the array " << endl;
#endif
            SnarlTreeRecordConstructor::records->emplace_back(pointer);
        }
    };

    struct ChainRecord : SnarlTreeRecord {

        ChainRecord() {};
        ChainRecord (size_t pointer, const vector<size_t>* tree_records){
            record_offset = pointer;
            records = tree_records;
            net_handle_record_t record_type= get_record_handle_type();
            if (record_type == NODE_HANDLE) {
                net_handle_record_t parent_type = SnarlTreeRecord(
                    NodeRecord(pointer, records).get_parent_record_offset(), records
                ).get_record_handle_type();
                assert(parent_type == ROOT_HANDLE || parent_type == SNARL_HANDLE);
            } else {
                assert(get_record_handle_type() == CHAIN_HANDLE);
            }
        }
        ChainRecord (net_handle_t net, const vector<size_t>* tree_records){
            record_offset = get_record_offset(net);
            records = tree_records;

            net_handle_record_t record_type = get_record_handle_type();
            if (record_type == NODE_HANDLE) {
                net_handle_record_t parent_type = SnarlTreeRecord(
                    NodeRecord(record_offset, records).get_parent_record_offset(), records
                ).get_record_handle_type();
                assert(parent_type == ROOT_HANDLE || parent_type == SNARL_HANDLE);
            } else {
                assert(get_record_handle_type() == CHAIN_HANDLE);
            }
        }

        virtual size_t get_node_count() const {
            if (get_record_handle_type() == NODE_HANDLE) {
                return 1;
            } else {
                return records->at(record_offset + CHAIN_NODE_COUNT_OFFSET);
            }
        }

        //Get the prefix sum value for this node (boundary node of a snarl in the chain)
        //pointer is a pointer into snarl_tree_records, to the beginning of the record for this node
        //So it'll point to the node id of the node we're looking at
        virtual int64_t get_prefix_sum_value(size_t pointer) const {
            if (get_record_handle_type() == NODE_HANDLE) {
                throw runtime_error("error: Trying to get chain distances from a node");
            }
            size_t val = records->at(pointer+CHAIN_NODE_PREFIX_SUM_OFFSET);
            return val == 0 ? std::numeric_limits<int64_t>::max() : val-1; 
        }
        virtual int64_t get_forward_loop_value(size_t pointer) const {
            if (get_record_handle_type() == NODE_HANDLE) {
                throw runtime_error("error: Trying to get chain distances from a node");
            }
            size_t val = records->at(pointer+CHAIN_NODE_FORWARD_LOOP_OFFSET);
            return val == 0 ? std::numeric_limits<int64_t>::max() : val-1; 
        }
        virtual int64_t get_reverse_loop_value(size_t pointer) const {
            if (get_record_handle_type() == NODE_HANDLE) {
                throw runtime_error("error: Trying to get chain distances from a node");
            }
            size_t val =  records->at(pointer+CHAIN_NODE_REVERSE_LOOP_OFFSET); 
            return val == 0 ? std::numeric_limits<int64_t>::max() : val-1; 
        }
        virtual size_t get_chain_component(size_t pointer) const {
            if (get_record_handle_type() == NODE_HANDLE) {
                throw runtime_error("error: Trying to get chain distances from a node");
            }
            if (get_record_type() != MULTICOMPONENT_CHAIN){
                throw runtime_error("error: Trying to get the component of a single component chain");
            }
            return records->at(pointer+CHAIN_NODE_COMPONENT_OFFSET); 
        }

        //Get the distance between the given node sides (relative to the orientation of the chain)
        //Nodes represent a tuple of <pointer, right_side, and length of the node>
        //This is the distance between the node sides, leaving the first and entering the second,
        //not including node lengths
        //TODO: I don't think we're allowing looping chains so I'm going to ignore them for now
        //TODO: Double check finding the distance for the same node
        virtual int64_t get_distance(tuple<size_t, bool, int64_t> node1, 
                             tuple<size_t, bool, int64_t> node2) const {
            if (get_record_handle_type() == NODE_HANDLE) {
                throw runtime_error("error: Trying to get chain distances from a node");
            } else if (get_record_type() == MULTICOMPONENT_CHAIN) {
                if (records->at(std::get<0>(node1) + CHAIN_NODE_COMPONENT_OFFSET) != 
                    records->at(std::get<0>(node2) + CHAIN_NODE_COMPONENT_OFFSET)) {
                    return std::numeric_limits<int64_t>::max();
                }
            }

            if (std::get<0>(node1) > std::get<0>(node2)) {
                //If the first node comes after the second in the chain, reverse them
                tuple<size_t, bool, size_t> tmp = node1;
                node1 = node2;
                node2 = tmp;

            }

            int64_t distance;
            bool is_looping_chain = get_start_id() == get_end_id(); 

            if (!std::get<1>(node1) && std::get<1>(node2)) {
                //Right of 1 and left of 2, so a simple forward traversal of the chain
                if (std::get<0>(node1) == std::get<0>(node2)) {
                    //If these are the same node, then the path would need to go around the node
                    distance = sum({get_forward_loop_value(std::get<0>(node1)), 
                                get_reverse_loop_value(std::get<0>(node1)),
                                std::get<2>(node1)});
                } else {
                    distance = minus(get_prefix_sum_value(std::get<0>(node2)) - get_prefix_sum_value(std::get<0>(node1)),
                         std::get<2>(node1));
                }
            } else if (!std::get<1>(node1) && !std::get<1>(node2)) {
                //Right side of 1 and right side of 2
                if (std::get<0>(node1) == std::get<0>(node2)) {
                    distance = get_forward_loop_value(std::get<0>(node2));

                } else {
                    distance = minus( sum({get_prefix_sum_value(std::get<0>(node2)) - get_prefix_sum_value(std::get<0>(node1)) , 
                                       std::get<2>(node2), 
                                       get_forward_loop_value(std::get<0>(node2))}), 
                                 std::get<2>(node1));
                }
            } else if (std::get<1>(node1) && std::get<1>(node2)) {
                //Left side of 1 and left side of 2
                if (std::get<0>(node1) == std::get<0>(node2)) {
                    distance = get_reverse_loop_value(std::get<0>(node1));

                } else {
                    distance = sum({get_prefix_sum_value(std::get<0>(node2)) - get_prefix_sum_value(std::get<0>(node1)),
                                    get_reverse_loop_value(std::get<0>(node1))});
                }
                if (is_looping_chain) {
                    //Check distance for taking loop in chain: from the first node left to the start, around the 
                    //chain loop, then the reverse loop of the second node
                    //This assumes that the length of the chain only includes the start/end node's length once,
                    //which it does but might change
                    distance = std::min(distance,
                                    sum({get_prefix_sum_value(std::get<0>(node1)), 
                                         minus(get_min_length(), get_prefix_sum_value(std::get<0>(node2))),
                                         get_reverse_loop_value(std::get<0>(node1))}));
                    cerr << " Looping chain distance: " << sum({get_prefix_sum_value(std::get<0>(node1)), 
                                         minus(get_min_length(), get_prefix_sum_value(std::get<0>(node2))),
                                         get_reverse_loop_value(std::get<0>(node1))}) << endl;
                }
            } else {
                assert(std::get<1>(node1) && !std::get<1>(node2));
                //Left side of 1 and right side of 2
                distance = sum({get_prefix_sum_value(std::get<0>(node2)) - get_prefix_sum_value(std::get<0>(node1)), 
                                get_reverse_loop_value(std::get<0>(node1)),
                                get_forward_loop_value(std::get<0>(node1)),
                                std::get<2>(node1)});

                if (is_looping_chain) {
                    //Check the distance going backwards around the chain
                    distance = std::min(distance,
                                    sum({get_prefix_sum_value(std::get<0>(node1)), 
                                         minus(minus(get_min_length(), get_prefix_sum_value(std::get<0>(node2))),
                                          std::get<2>(node2))}));
                        cerr << " Looping chain distance: " << sum({get_prefix_sum_value(std::get<0>(node1)), 
                                         minus(minus(get_min_length(), get_prefix_sum_value(std::get<0>(node2))),
                                          std::get<2>(node2))}) << endl;
                }
            }
            return distance;
        }

        ////////////////////////// methods for navigating the snarl tree from this chain

        //Get the offset into snarl_tree_records of the first node in the chain
        virtual size_t get_first_node_offset() const {
            if (get_record_handle_type() == NODE_HANDLE) {
                throw runtime_error("error: Trying to get chain traversal from a node");
            }
            return record_offset + CHAIN_RECORD_SIZE;
        }

        //Given a pointer to a child node, return the next child. This includes nodes and snarls.
        //The node pointer points to the node id, snarl pointer points to the start of the snarl record
        //bool is true if it is a snarl, false if it is a node
        //go_left is true if we are traversing the chain right to left
        //returns 0 if this was the last (/first) node in the chain
        virtual pair<size_t, bool> get_next_child(const pair<size_t, bool> pointer, bool go_left) const {
            //If this is a multicomponent chain, then the size of each node record in the chain is bigger
            size_t extra_node_record_size = get_record_type() == MULTICOMPONENT_CHAIN ? 1 : 0;
            if (get_record_handle_type() == NODE_HANDLE) {
                throw runtime_error("error: Trying to get chain traversal from a node");
            }
            if (pointer.second) {
                //This is a snarl
                if (go_left) {
                    return make_pair(pointer.first - CHAIN_NODE_RECORD_SIZE - extra_node_record_size, false);
                } else {
                    size_t snarl_record_length = SnarlRecord(pointer.first, records).record_size();
                    return make_pair(pointer.first + snarl_record_length + 1, false);
                }
            } else {
                //This is a node
                if (go_left) {
                    //Looking left in the chain
                    if (pointer.first == record_offset + CHAIN_RECORD_SIZE) {
                        //If this is the first node in the chain
                        return make_pair(0, false);
                    }
                    size_t snarl_record_size = records->at(pointer.first-1);
                    if (snarl_record_size == 0) {
                        //Just another node to the left
                        return make_pair(pointer.first-CHAIN_NODE_RECORD_SIZE - extra_node_record_size, false);
                    } else {
                        //There is a snarl to the left of this node
                        return make_pair(pointer.first - snarl_record_size - 1, true);
                    }
                } else {
                    //Looking right in the chain
                    if (records->at(pointer.first+CHAIN_NODE_SNARL_SIZE_OFFSET + extra_node_record_size) == 0 &&
                        records->at(pointer.first+CHAIN_NODE_RECORD_SIZE+ extra_node_record_size) == 0) {
                        //If this is the last node in the chain
                        return make_pair(0, false);
                    }
                    size_t snarl_record_size = records->at(pointer.first+CHAIN_NODE_SNARL_SIZE_OFFSET+extra_node_record_size);
                    return make_pair(pointer.first+CHAIN_NODE_RECORD_SIZE+extra_node_record_size, snarl_record_size != 0);
                }
            }
        }
        //The same thing but take and return a net_handle_t
        //return the same net_handle_t if this is the end of the chain
        virtual net_handle_t get_next_child(const net_handle_t& net_handle, bool go_left) const {
            //get the next child in the chain. net_handle must point to a snarl or node in the chain 
            net_handle_record_t handle_type = get_handle_type(net_handle);
            net_handle_record_t record_type = get_record_handle_type();
            if (record_type == NODE_HANDLE) {
                //If this record is a node pretending to be a chain, then there is no next child
                assert(handle_type == NODE_HANDLE);
                assert(get_record_offset(net_handle) == record_offset);
                return net_handle;
            }

            //Get the current pointer, pointing at the net_handle in the chain
            if (handle_type == NODE_HANDLE) {
                //If this net handle is a node, then it's rank in parent points to it in the chain
                NodeRecord node_record(get_record_offset(net_handle), records);
                assert(node_record.get_parent_record_offset() == record_offset);
                pair<size_t, bool> next_pointer = get_next_child(
                    make_pair(node_record.get_rank_in_parent(), true), go_left);
                if (next_pointer.first == 0 ){
                    return net_handle;
                }
                bool next_is_reversed_in_parent = NodeRecord(
                        get_offset_from_id(records->at(next_pointer.first)), records
                    ).get_is_reversed_in_parent();
                
                return get_net_handle(get_offset_from_id(records->at(next_pointer.first)),
                                  go_left == next_is_reversed_in_parent ? START_END : END_START,
                                  next_pointer.first ? SNARL_HANDLE : NODE_HANDLE);
            } else {
                //Otherwise, it is a snarl and we can use the snarl's offset, since it exists in 
                //the chain
                assert(handle_type == SNARL_HANDLE) ;
                pair<size_t, bool> next_pointer = get_next_child(
                    make_pair(get_record_offset(net_handle), true), go_left);
                if (next_pointer.first == 0 ){
                    return net_handle;
                }
                return get_net_handle(next_pointer.first,
                                  go_left ? END_START : START_END,
                                  next_pointer.first ? SNARL_HANDLE : NODE_HANDLE);
            }
        }

        virtual bool for_each_child(const std::function<bool(const net_handle_t&)>& iteratee) const {

            if (get_record_handle_type() == NODE_HANDLE) {
                //If this is a node pretending to be a chain, just do it for the node
                return iteratee(get_net_handle(record_offset, START_END, NODE_HANDLE));
            }


            //If this is a node, then the offset of the node in the chain, false
            //If it is a snarl, then the offset of the snarl record, true
            pair<size_t, bool> current_child (get_first_node_offset(), false);
            bool is_first = true;

            while (current_child.first != 0) {
                net_handle_t child_handle = current_child.second 
                    ? get_net_handle (current_child.first, START_END, SNARL_HANDLE) 
                    : get_net_handle (get_offset_from_id(records->at(current_child.first)), START_END, NODE_HANDLE);
                if (!is_first && current_child == make_pair(get_first_node_offset(), false)){
                    //Don't look at the first node a second time
                    return true;
                }

                bool result = iteratee(child_handle); 
                if (result == false) {
                    return false;
                }
                current_child = get_next_child(current_child, false); 
                is_first = false;
            }
            return true;
        }
    };

    struct ChainRecordConstructor : ChainRecord , SnarlTreeRecordConstructor {
        //Constructor for a chain record
        //Since the size of the vector will be unknown (since we don't know how big the snarls are),
        //Add nodes and snarls as they come up. Assumes that the memory has already been reserved but
        //not allocated yet

        //TODO: I don't think I even need node count
        ChainRecordConstructor() {}
        ChainRecordConstructor (size_t pointer, record_t type, size_t node_count, vector<size_t>* records){
            assert(type == CHAIN || 
                   type == DISTANCED_CHAIN ||
                   type == MULTICOMPONENT_CHAIN);
            ChainRecord::record_offset = pointer;
            SnarlTreeRecordConstructor::record_offset = pointer;
            SnarlTreeRecordConstructor::records = records;
            ChainRecord::records = records;
#ifdef debug_indexing
            cerr << " Resizing array to add chain: length " << SnarlTreeRecordConstructor::records->size() << " -> "  << SnarlTreeRecordConstructor::records->size() + CHAIN_RECORD_SIZE << endl;
#endif
            SnarlTreeRecordConstructor::records->resize(SnarlTreeRecordConstructor::records->size() + CHAIN_RECORD_SIZE, 0);
            set_node_count(node_count);
            set_record_type(type);
        }


        /* Functions to add children to a chain. Assumes that the chain is well formed up to here
         * These will always be called in order going forward in the chain, starting and ending 
         * with add_node, with either add_snarl or add_trivial_snarl in between them. 
         * Nodes will be [node id, prefix sum, fd loop, rev loop], and snarls will be flanked by
         * the record size. A trivial snarl will only be a record size of 0. 
         * Two 0's where the snarl size and next node id should be indicates the end of the chain
         */

        void add_node(id_t id, int64_t prefix_sum, int64_t forward_loop, int64_t reverse_loop) {
            assert(ChainRecord::get_record_type() != MULTICOMPONENT_CHAIN);
#ifdef debug_indexing
            cerr << SnarlTreeRecordConstructor::records->size() << " - " << SnarlTreeRecordConstructor::records->size() + 3 << " Adding chain's child node " << id << " to the end of the array (values: " << prefix_sum << "," << forward_loop << "," << reverse_loop << ")" << endl;
#endif
            SnarlTreeRecordConstructor::records->emplace_back(id);
            SnarlTreeRecordConstructor::records->emplace_back(prefix_sum==std::numeric_limits<int64_t>::max() ? 0 : prefix_sum+1);
            SnarlTreeRecordConstructor::records->emplace_back(forward_loop==std::numeric_limits<int64_t>::max() ? 0 : forward_loop+1);
            SnarlTreeRecordConstructor::records->emplace_back(reverse_loop==std::numeric_limits<int64_t>::max() ? 0 : reverse_loop+1);
        }
        void add_node(id_t id, int64_t prefix_sum, int64_t forward_loop, int64_t reverse_loop, size_t component) {
            assert(ChainRecord::get_record_type() == MULTICOMPONENT_CHAIN);
#ifdef debug_indexing
            cerr << SnarlTreeRecordConstructor::records->size() << " - " << SnarlTreeRecordConstructor::records->size() + 4 << " Adding chain's child node the end of the array " << endl;
#endif
            SnarlTreeRecordConstructor::records->emplace_back(id);
            SnarlTreeRecordConstructor::records->emplace_back(prefix_sum==std::numeric_limits<int64_t>::max() ? 0 : prefix_sum+1);
            SnarlTreeRecordConstructor::records->emplace_back(forward_loop==std::numeric_limits<int64_t>::max() ? 0 : forward_loop+1);
            SnarlTreeRecordConstructor::records->emplace_back(reverse_loop==std::numeric_limits<int64_t>::max() ? 0 : reverse_loop+1);
            SnarlTreeRecordConstructor::records->emplace_back(component);
        }
        void set_node_count(size_t node_count) {
#ifdef debug_indexing
            cerr << ChainRecord::record_offset + CHAIN_NODE_COUNT_OFFSET << " set chain node count " << node_count << endl;
            assert(SnarlTreeRecordConstructor::records->at(ChainRecord::record_offset + CHAIN_NODE_COUNT_OFFSET) == 0);
#endif
            SnarlTreeRecordConstructor::records->at(ChainRecord::record_offset + CHAIN_NODE_COUNT_OFFSET) = node_count;
        }
        void add_trivial_snarl() {
#ifdef debug_indexing
            cerr << SnarlTreeRecordConstructor::records->size() << "  Adding chain's trivial snarl to the end of the array " << endl;
#endif
            SnarlTreeRecordConstructor::records->emplace_back(0);
        }
        //Add a snarl to the end of the chain and return a SnarlRecordConstructor pointing to it
        SnarlRecordConstructor add_snarl(size_t snarl_size, record_t type) {
            size_t snarl_record_size = SnarlRecord::record_size(type, snarl_size);
#ifdef debug_indexing
            cerr << SnarlTreeRecordConstructor::records->size() << " Adding child snarl length to the end of the array " << endl;
#endif
            SnarlTreeRecordConstructor::records->emplace_back(snarl_record_size);
            SnarlRecordConstructor snarl_record(snarl_size, SnarlTreeRecordConstructor::records, type);
#ifdef debug_indexing
            cerr << SnarlTreeRecordConstructor::records->size() << " Adding child snarl length to the end of the array " << endl;
#endif
            SnarlTreeRecordConstructor::records->emplace_back(snarl_record_size);
            return snarl_record;
        }

        void finish_chain(){
#ifdef debug_indexing
            cerr << SnarlTreeRecordConstructor::records->size()  << " - " <<  SnarlTreeRecordConstructor::records->size()+1 << " Adding the last two chain 0's to the end of the array " << endl;
#endif
            SnarlTreeRecordConstructor::records->emplace_back(0);
            SnarlTreeRecordConstructor::records->emplace_back(0);
        }
    };


private:
    ////////////////////// More methods for dealing with net_handle_ts
    SnarlTreeRecord get_snarl_tree_record(const handlegraph::net_handle_t& net_handle) const {
        return SnarlTreeRecord(get_record_offset(net_handle), &snarl_tree_records);
    }
    SnarlTreeRecord get_node_record(const handlegraph::net_handle_t& net_handle) const {
        return NodeRecord(get_record_offset(net_handle), &snarl_tree_records); 
    }
    SnarlTreeRecord get_snarl_record(const handlegraph::net_handle_t& net_handle) const {
        return SnarlRecord(get_record_offset(net_handle), &snarl_tree_records); 
    }
    SnarlTreeRecord get_chain_record(const handlegraph::net_handle_t& net_handle) const {
        return ChainRecord(get_record_offset(net_handle), &snarl_tree_records); 
    }

    const static connectivity_t endpoints_to_connectivity(endpoint_t start, endpoint_t end) {
        if (start == START && end == START) {
            return START_START;
        } else if (start == START && end == END) {
            return START_END;
        } else if (start == START && end == TIP) {
            return START_TIP;
        } else if (start == END && end == START) {
            return END_START;
        } else if (start == END && end == END) {
            return END_END;
        } else if (start == END && end == TIP) {
            return END_TIP;
        } else if (start == TIP && end == START) {
            return TIP_START;
        } else if (start == TIP && end == END) {
            return TIP_END;
        } else if (start == TIP && end == TIP) {
            return TIP_TIP;
        } else {
            throw runtime_error("error: invalid endpoints");
        }
    }
    const static endpoint_t get_start_endpoint(connectivity_t connectivity) {
        endpoint_t start_endpoint;
        if (connectivity == START_START || connectivity == START_END || connectivity == START_TIP){
            start_endpoint = START;
        } else if (connectivity == END_START || connectivity == END_END || connectivity == END_TIP){
            start_endpoint = END;
        } else if (connectivity == TIP_START || connectivity == TIP_END || connectivity == TIP_TIP){
            start_endpoint = TIP;
        } else {
            throw runtime_error("error: invalid connectivity");
        }
        return start_endpoint;
    }
    const static endpoint_t get_start_endpoint(net_handle_t net) {
        return get_start_endpoint(get_connectivity(net));
    }
    const static endpoint_t get_end_endpoint(connectivity_t connectivity) {
        endpoint_t end_endpoint;
        if (connectivity == START_START || connectivity == END_START || connectivity == TIP_START){
            end_endpoint = START;
        } else if (connectivity == START_END || connectivity == END_END || connectivity == TIP_END){
            end_endpoint = END;
        } else if (connectivity == START_TIP || connectivity == END_TIP || connectivity == TIP_TIP){
            end_endpoint = TIP;
        } else {
            throw runtime_error("error: invalid connectivity");
        }
        return end_endpoint;
    }
    const static endpoint_t get_end_endpoint(const net_handle_t& net) {
        return get_end_endpoint(get_connectivity(net));
    }
    const static pair<endpoint_t, endpoint_t> connectivity_to_endpoints(const connectivity_t& connectivity) {
        return make_pair(get_start_endpoint(connectivity), get_end_endpoint(connectivity));
    }
    public:
    string net_handle_as_string(const net_handle_t& net) const {
        net_handle_record_t type = get_handle_type(net);
        SnarlTreeRecord record (net, &snarl_tree_records);
        net_handle_record_t record_type = record.get_record_handle_type();
        string result;
        if (type == ROOT_HANDLE) {
            return "root";
        } else if (type == NODE_HANDLE) {
            return  "node " + std::to_string( get_node_id_from_offset(get_record_offset(net))) + (ends_at(net) == START ? "rev" : "fd");
        } else if (type == SNARL_HANDLE) {
            if (record.get_record_type() == ROOT) {
                return "root snarl";
            }
            result += "snarl ";        
        } else if (type == CHAIN_HANDLE && record_type == NODE_HANDLE) {
            return  "node " + std::to_string( get_node_id_from_offset(get_record_offset(net))) 
                   + std::to_string( get_node_id_from_offset(get_record_offset(net))) + (ends_at(net) == START ? "rev" : "fd")
                   + " pretending to be a chain";
        } else if (type == CHAIN_HANDLE) {
            result += "chain ";
        } else if (type == SENTINEL_HANDLE) {
            result += "sentinel of snarl ";
        } else { 
            throw runtime_error("error: Unknown net_handle_t type");
        }
        result += ( std::to_string(record.get_start_id()) 
                + (record.get_start_orientation() ? "rev" : "fd") 
                + "->"
                + std::to_string(record.get_end_id())   
                + (record.get_end_orientation() ? "rev" : "fd"));
        result += "traversing ";
        result += (starts_at(net) == START ? "start" : (starts_at(net) == END ? "end" : "tip"));
        result += "->";
        result += (ends_at(net) == START ? "start" : (ends_at(net) == END ? "end" : "tip"));
        return result;
    }


protected:

    /*
     * A structure to store everything in the distance index, but not in just one vector to make it easier to construct
     * This can also be used to combine distance indexes
     */
    enum temp_record_t {TEMP_CHAIN=0, TEMP_SNARL, TEMP_NODE, TEMP_ROOT};

    class TemporaryDistanceIndex{
    public:
        TemporaryDistanceIndex();
        ~TemporaryDistanceIndex();
        TemporaryDistanceIndex(const HandleGraph* graph, const HandleGraphSnarlFinder* snarl_finder, size_t size_limit);

        //Get a string of the start and end of a structure
        string structure_start_end_as_string(pair<temp_record_t, size_t> index) const;
    
    protected:
        id_t min_node_id;
        id_t max_node_id;
        size_t root_structure_count=0; //How many things are in the root
        size_t index_size;//TODO: This will change depending on how the actual index is represented

        //This will actually store each individual record separately, and each 
        //will have real pointers to their parents/children (as opposed to offsets)
        struct TemporaryRecord {
        };
        struct TemporaryChainRecord : TemporaryRecord {
            id_t start_node_id;
            bool start_node_rev;
            id_t end_node_id;
            bool end_node_rev;
            int64_t end_node_length;
            //Type of the parent and offset into the appropriate vector
            //(TEMP_ROOT, 0) if this is a root level chain
            pair<temp_record_t, size_t> parent;
            int64_t min_length;//Including boundary nodes
            int64_t max_length;
            vector<pair<temp_record_t, size_t>> children; //All children, both nodes and snarls, in order
            //Distances for the chain, one entry per node
            vector<int64_t> prefix_sum;
            vector<int64_t> forward_loops;
            vector<int64_t> backward_loops;
            vector<size_t> chain_components;//Which component does each node belong to, usually all 0s
            size_t rank_in_parent;
            bool reversed_in_parent;
            bool is_trivial;
            bool is_tip = false;
            //What is the index of this record in root_snarl_components
            size_t root_snarl_index = std::numeric_limits<size_t>::max();
        };
        struct TemporarySnarlRecord : TemporaryRecord{
            id_t start_node_id;
            bool start_node_rev;
            int64_t start_node_length;
            id_t end_node_id;
            bool end_node_rev;
            int64_t end_node_length;
            size_t node_count;
            int64_t min_length; //Not including boundary nodes
            int64_t max_length;
            pair<temp_record_t, size_t> parent;
            vector<pair<temp_record_t, size_t>> children; //All children, nodes and chains, in arbitrary order
            unordered_set<size_t> tippy_child_ranks; //The ranks of children that are tips
            unordered_map<pair<pair<size_t, bool>, pair<size_t, bool>>, int64_t> distances;
            size_t rank_in_parent;
            bool reversed_in_parent;
            //The minimum distances to go into and out of the snarl from the start or end nodes
            int64_t loop_start;
            int64_t loop_end;
            bool is_trivial;
            bool is_tip = false;
            bool is_root_snarl = false;
        };
        struct TemporaryNodeRecord : TemporaryRecord{
            TemporaryNodeRecord() :
            node_id(0), parent(make_pair(TEMP_ROOT, 0)), node_length(0), 
            rank_in_parent(0), reversed_in_parent(false){
            }
            id_t node_id;
            pair<temp_record_t, size_t> parent;
            size_t node_length;
            size_t rank_in_parent;
            bool reversed_in_parent;
            bool is_tip = false;
        };


        vector<pair<temp_record_t, size_t>> components;
        vector<pair<temp_record_t, size_t>> root_snarl_components;
        vector<TemporaryChainRecord> temp_chain_records;
        vector<TemporarySnarlRecord> temp_snarl_records;
        vector<TemporaryNodeRecord> temp_node_records;
        friend class SnarlDistanceIndex;

        //Fill in the temporary snarl record with distances
        void populate_snarl_index(TemporaryDistanceIndex::TemporarySnarlRecord& temp_snarl_record, 
            pair<temp_record_t, size_t> snarl_index, size_t size_limit, const HandleGraph* graph) ;
    };

    //Given an arbitrary number of temporary indexes, produce the final one
    //Each temporary index must be a separate connected component
    vector<size_t> get_snarl_tree_records(const vector<const TemporaryDistanceIndex*>& temporary_indexes, const HandleGraph* graph);
    friend class TemporaryDistanceIndex;

    private:
    ///Add integers, returning max() if any of them are max()
    static int64_t sum(const vector<int64_t> vals) {
        int64_t sum = 0;
        for (const int64_t& x : vals) {
            if (x ==  std::numeric_limits<int64_t>::max()) {
                return std::numeric_limits<int64_t>::max();
            } else {
                sum += x;
            }
        }
        return sum;
    }
    static int64_t minus(int64_t x, int64_t y) {
        if (x == std::numeric_limits<int64_t>::max()) {
            return numeric_limits<int64_t>::max();
        } else if (y == std::numeric_limits<int64_t>::max()) {
            return -numeric_limits<int64_t>::max();
        } else {
            return x - y;
        }
    }
};


}

#endif
