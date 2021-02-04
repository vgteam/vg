#ifndef VG_MIN_DISTANCE_HPP_INCLUDED
#define VG_MIN_DISTANCE_HPP_INCLUDED

#include <handlegraph/snarl_decomposition.hpp>
#include "snarls.hpp"

using namespace sdsl;
using namespace handlegraph;
namespace vg { 

/**
 * The distance index. Stores minimum distances among nodes in each 
 * netgraph and chain.
 */

//TODO: Check about looping chains and unary snarls

class SnarlDistanceIndex : public SnarlDecomposition {

private:
    
    //This stores all records for the root, nodes, chains/snarls, and snarls' children
    //
    //It's really made up of five types of vectors: 
    //
    //- The (single) root vector has the format:
    //  [root tag, # connected components, [pointer to node/snarl/chain record] x N]
    //  The root vector stores the root of every connected component, which can be a 
    //  node, snarl, or chain
    //
    //- The (single) node vector stores a record for each node and has the format:
    //  [# nodes, min_node_id, [node tag, pointer to parent, node length, rank in parent,
    //      component #] x N TODO: Could add extra stuff here, I'm not sure where to define it
    //
    //- A chain record for each chain, which is interspersed with snarl records:
    //  [chain tag, #nodes, pointer to parent, min length, max length, rank in parent, 
    //      [node id, prefix sum, fd loop, rev loop, has snarl, (snarl record), # nodes in snarl,
    //       has snarl] x N] (plus an extra node id, prefix sum, fd loop, rev loop at the end for 
    //       the last node)
    //
    //- A snarl record for each snarl, which are stuck in chains
    //  [snarl tag, # nodes, pointer to parent, min length, max length, rank in parent, 
    //   pointer to children (in child vector), distance vector]
    //   Trivial snarls will still be contained within a chain
    //   Single nodes will technically be contained in a chain, but not stored as a chain
    //   THe first node in a snarl (rank 0 and first in child vector) will always be the start node,
    //   and the last node will always be the end node
    //   TODO: Might also add start and end node ids and maybe orientations too
    //   For the first and last nodes, we only care about the node sides pointing in
    //   Each node side in the actual distance matrix will be 2*rank-1 for the left side, and 
    //   2*rank for the right side, and 0 for the start, 2*(num_nodes-1)-1 for the end
    //
    //- The (single) child vector, listing children in snarls
    //  [child vector tag, (pointer to records) x N
    //  Each snarl will have a pointer into here, and will also know how many children it has

    //TODO: I'm not sure this should be static, at least for construction?
    static vector<size_t> snarl_tree_records;

    //The "tags" for defining what kind of record we're looking at will be a RecordType and a 
    //bit vector indicating connectivity. The bit vector will be the last 6 bits of the tag
    //
    //Each bit represents one type of connectivity:
    //start-start, start-end, start-tip, end-end, end-tip, tip-tip
    //std::bitset<6> connectivity;
    //
    //The remainder of the tag will be the RecordType of the record
    //NODE, SNARL, and CHAIN indicate that they don't store distances.
    //SIMPLE_SNARL is a snarl with all children connecting only to the boundary nodes
    //OVERSIZED_SNARL only stores distances to the boundaries
    //TODO: Make simple snarls work
    //TODO: Maybe also add a tag for trivial snarls/chains, or for defining connectivity
    //TODO: Unary snarls? Looping chains?
    enum RecordType {ROOT=1, 
                     NODE, DISTANCED_NODE, 
                     SNARL, DISTANCED_SNARL, SIMPLE_SNARL, OVERSIZED_SNARL, 
                     CHAIN, DISTANCED_CHAIN, 
                     CHILDREN};

    enum ConnectivityType { START_START=1, START_END, START_TIP, 
                            END_START, END_END, END_TIP, 
                            TIP_START, TIP_END, TIP_TIP};

public:

private:



    //Define a struct for interpreting each type of snarl tree node record (For node, snarl, chain)
    //TODO: This might be overkill but I want it to be easy to change what gets stored in the index
    //
    struct snarl_tree_record_t {


        //The offset of the start of this record in snarl_tree_records
        size_t record_offset;

        snarl_tree_record_t();
        snarl_tree_record_t (size_t pointer){
            record_offset = pointer;
            RecordType type = static_cast<RecordType>(snarl_tree_records[record_offset]>>6);
            assert(type == ROOT || type == NODE ||  type == SNARL || type == DISTANCED_SNARL || 
                   type == OVERSIZED_SNARL || type == CHAIN || type == DISTANCED_CHAIN );
        }


        //What type of snarl tree node is this?
        //This will be the first value of any record
        RecordType get_record_type() {
            return static_cast<RecordType>(snarl_tree_records[record_offset] >> 6);
        }

        bool is_start_start_connected() {return snarl_tree_records[record_offset] & 32;}
        bool is_start_end_connected() {return snarl_tree_records[record_offset] & 16;}
        bool is_start_tip_connected() {return snarl_tree_records[record_offset] & 8;}
        bool is_end_end_connected() {return snarl_tree_records[record_offset] & 4;}
        bool is_end_tip_connected() {return snarl_tree_records[record_offset] & 2;}
        bool is_tip_connected() {return snarl_tree_records[record_offset] & 1;}
        void set_start_start_connected() {snarl_tree_records[record_offset] = snarl_tree_records[record_offset] | 32;}
        void set_start_end_connected() {snarl_tree_records[record_offset] = snarl_tree_records[record_offset] | 16;}
        void set_start_tip_connected() {snarl_tree_records[record_offset] = snarl_tree_records[record_offset] | 8;}
        void set_end_end_connected() {snarl_tree_records[record_offset] = snarl_tree_records[record_offset] | 4;}
        void set_end_tip_connected() {snarl_tree_records[record_offset] = snarl_tree_records[record_offset] | 2;}
        void set_tip_connected() {snarl_tree_records[record_offset] = snarl_tree_records[record_offset] | 1;}

        //How many things get stored in this record_t?
        //This is used to allocate memory for the node's record
        size_t record_size;

        //Get and set a pointer to this chain's parent
        size_t get_parent_record_pointer() {
            RecordType type = static_cast<RecordType>(snarl_tree_records[record_offset] >> 6);
            if (type == ROOT ) {
                return 0;
            } else if (type == NODE) {
                return snarl_tree_records[record_offset + 1];
            } else if (type == SNARL || type == DISTANCED_SNARL || type == OVERSIZED_SNARL)  {

                return snarl_tree_records[record_offset + 2];
            } else if (type == CHAIN || type == DISTANCED_CHAIN)  {

                return snarl_tree_records[record_offset + 2];
            } else {
                throw runtime_error("error: trying to access a snarl tree node of the wrong type");
            }
        };
        void set_parent_record_pointer(size_t pointer){
            RecordType type = static_cast<RecordType>(snarl_tree_records[record_offset] >> 6);
            if (type == NODE) {
                snarl_tree_records[record_offset + 1] = pointer;
            } else if (type == SNARL || type == DISTANCED_SNARL || type == OVERSIZED_SNARL)  {

                snarl_tree_records[record_offset + 2] = pointer;
            } else if (type == CHAIN || type == DISTANCED_CHAIN)  {

                snarl_tree_records[record_offset + 2] = pointer;
            } else {
                throw runtime_error("error: trying to access a snarl tree node of the wrong type");
            }
        };

        //Get and set the minimum length (distance from start to end, including boundaries for 
        //snarls and chains, just node length for nodes)
        size_t get_min_length() {
            RecordType type = static_cast<RecordType>(snarl_tree_records[record_offset] >> 6);
            if (type == NODE) {
                return snarl_tree_records[record_offset + 2];
            } else if (type == SNARL || type == DISTANCED_SNARL || type == OVERSIZED_SNARL)  {
                return snarl_tree_records[record_offset + 3];
            } else if (type == CHAIN || type == DISTANCED_CHAIN)  {
                return snarl_tree_records[record_offset + 3];
            } else {
                throw runtime_error("error: trying to access a snarl tree node of the wrong type");
            }
        };
        void set_min_length(size_t length) {
            RecordType type = static_cast<RecordType>(snarl_tree_records[record_offset] >> 6);
            if (type == NODE) {
                snarl_tree_records[record_offset + 2] = length;
            } else if (type == SNARL || type == DISTANCED_SNARL || type == OVERSIZED_SNARL)  {
                snarl_tree_records[record_offset + 3] = length;
            } else if (type == CHAIN || type == DISTANCED_CHAIN)  {
                snarl_tree_records[record_offset + 3] = length;
            } else {
                throw runtime_error("error: trying to access a snarl tree node of the wrong type");
            }
        };

        //Get and set this snarl's maximum length
        //This isn't actually a maximum, it's the maximum among minimum distance paths 
        //through each node in the snarl
        size_t get_max_length() {
            RecordType type = static_cast<RecordType>(snarl_tree_records[record_offset] >> 6);
            if (type == NODE) {
                return snarl_tree_records[record_offset + 2];
            } else if (type == SNARL || type == DISTANCED_SNARL || type == OVERSIZED_SNARL)  {
                return snarl_tree_records[record_offset + 4];
            } else if (type == CHAIN || type == DISTANCED_CHAIN)  {
                return snarl_tree_records[record_offset + 4];
            } else {
                throw runtime_error("error: trying to access a snarl tree node of the wrong type");
            }
        };
        void set_max_length(size_t length) {
            RecordType type = static_cast<RecordType>(snarl_tree_records[record_offset] >> 6);
            if (type == NODE) {
                snarl_tree_records[record_offset + 2] = length;
            } else if (type == SNARL || type == DISTANCED_SNARL || type == OVERSIZED_SNARL)  {
                snarl_tree_records[record_offset + 4] = length;
            } else if (type == CHAIN || type == DISTANCED_CHAIN)  {
                snarl_tree_records[record_offset + 4] = length;
            } else {
                throw runtime_error("error: trying to access a snarl tree node of the wrong type");
            }
        };

        //Get and set this snarl's rank in its parent
        size_t get_rank_in_parent() {
            RecordType type = static_cast<RecordType>(snarl_tree_records[record_offset] >> 6);
            if (type == NODE) {
                return snarl_tree_records[record_offset + 3];
            } else if (type == SNARL || type == DISTANCED_SNARL || type == OVERSIZED_SNARL)  {
                return snarl_tree_records[record_offset + 5];
            } else if (type == CHAIN || type == DISTANCED_CHAIN)  {
                return snarl_tree_records[record_offset + 5];
            } else {
                throw runtime_error("error: trying to access a snarl tree node of the wrong type");
            }
        };
        void set_rank_in_parent(size_t rank) {
            RecordType type = static_cast<RecordType>(snarl_tree_records[record_offset] >> 6);
            if (type == NODE) {
                snarl_tree_records[record_offset + 3] = rank;
            } else if (type == SNARL || type == DISTANCED_SNARL || type == OVERSIZED_SNARL)  {
                snarl_tree_records[record_offset + 5] = rank;
            } else if (type == CHAIN || type == DISTANCED_CHAIN)  {
                snarl_tree_records[record_offset + 5] = rank;
            } else {
                throw runtime_error("error: trying to access a snarl tree node of the wrong type");
            }
        };
    };

    struct node_record_t : snarl_tree_record_t {

        size_t record_size = 5; //TODO: These numbers include distances

        node_record_t (size_t pointer, bool new_record=false) {
            record_offset = pointer;
            if (new_record) {
                snarl_tree_records.resize(snarl_tree_records.size() + record_size, 0);
                snarl_tree_records[record_offset] = NODE << 6;
            }
            assert(static_cast<RecordType>(snarl_tree_records[record_offset]>>6) == NODE);
        }

        //TODO: This one is a bit redundant but fine I think
        size_t get_node_length() {
            return snarl_tree_records[record_offset + 2];
        }
        void set_node_length(size_t length) {
            snarl_tree_records[record_offset + 2] = length;
        }

        size_t get_root_component() {
            return snarl_tree_records[record_offset + 4];
        }
        void set_root_component(size_t component) {
            snarl_tree_records[record_offset + 4] = component;
        }

    };

    struct snarl_record_t : snarl_tree_record_t {

        size_t record_size = 7;

        snarl_record_t (size_t pointer){
            record_offset = pointer;
            RecordType type = static_cast<RecordType>(snarl_tree_records[record_offset]>>6);
            assert(type == SNARL || 
                   type == DISTANCED_SNARL ||
                   type == OVERSIZED_SNARL);
        }

        snarl_record_t (size_t pointer, RecordType type, size_t node_count){
            //Constructor for making a new record, including allocating memory.
            //Assumes that this is the latest record being made, so pointer will be the end of
            //the array and we need to allocate extra memory past it
            //TODO:I'm not sure yet how memory will actually be allocated
            record_offset = pointer;
            size_t extra_size = snarl_record_size(type, node_count);
            snarl_tree_records.resize(snarl_tree_records.size() + extra_size, 0);
            set_node_count(node_count);
            snarl_tree_records[record_offset] = type << 6;
        }

        //How big is the entire snarl record?
        size_t snarl_record_size (RecordType type, size_t node_count) {
            if (type == SNARL){
                //For a normal snarl, its just the record size and the pointers to children
                return record_size + node_count; 
            } else if (type == DISTANCED_SNARL) {
                //For a normal min distance snarl, record size and the pointers to children, and
                //matrix of distances
                size_t node_side_count = node_count * 2 - 2;
                return record_size + node_count + (((node_side_count+1)*node_side_count) / 2);
            } else if (type ==  OVERSIZED_SNARL){
                //For a large min_distance snarl, record the side, pointers to children, and just 
                //the min distances from each node side to the two boundary nodes
                size_t node_side_count = node_count * 2 - 2;
                return record_size + node_count + (node_side_count * 2);
            }
        }
        size_t snarl_record_size() { 
            RecordType type = static_cast<RecordType>(snarl_tree_records[record_offset]>>6);
            return snarl_record_size(type, get_node_count()); 
        }


        //Get the index into the distance vector for the calculating distance between the given node sides
        int64_t get_distance_vector_offset(size_t rank1, bool right_side1, size_t rank2, bool right_side2) {

            //how many node sides in this snarl
            size_t node_side_count = (get_node_count()-1) * 2; 

            //make sure we're looking at the correct node side
            if (rank1 == get_node_count()) {
                rank1 = node_side_count-1;
            } else if (rank1 != 0) {
                rank1 = rank1 * 2 - 1;
                if (right_side1) {
                    rank1 += 1;
                }
            }
            if (rank2 == get_node_count()) {
                rank2 = node_side_count-1;
            } else if (rank2 != 0) {
                rank2 = rank2 * 2 - 1;
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
            
            RecordType type = static_cast<RecordType>(snarl_tree_records[record_offset]>>6);
            if (type == SNARL) {
                throw runtime_error("error: trying to access distance in a distanceless snarl tree");
            } else if (type == DISTANCED_SNARL) {
                //normal distance index
                size_t k = node_side_count-rank1;
                return (((node_side_count+1) * node_side_count)/2) - (((k+1)*k) / 2) + rank2 - rank1;
            } else if (type ==  OVERSIZED_SNARL) {
                //abbreviated distance index storing only the distance from each node side to the
                //start and end
                if (rank1 == 0) {
                    return rank2;
                } else if (rank2 == node_side_count-1) {
                    return get_node_count() + rank1;
                } else {
                    throw runtime_error("error: trying to access distance in an oversized snarl");
                }
            }
        }

        //TODO: I want to also add a function to do this given a node id or net_handle_t instead of rank
        //Get and set the distances between two node sides in the graph
        //Ranks identify which node, sides indicate node side: false for left, true for right
        int64_t get_distance(size_t rank1, bool right_side1, size_t rank2, bool right_side2) {

            //offset of the start of the distance vectors in snarl_tree_records
            size_t distance_vector_start = record_offset + get_node_count();
            //Offset of this particular distance in the distance vector
            size_t distance_vector_offset = get_distance_vector_offset(rank1, right_side1, rank2, right_side2);

            return snarl_tree_records[distance_vector_start + distance_vector_offset] + 1;

        }
        void set_distance(size_t rank1, bool right_side1, size_t rank2, bool right_side2, int64_t distance) {
            //offset of the start of the distance vectors in snarl_tree_records
            size_t distance_vector_start = record_offset + get_node_count();
            //Offset of this particular distance in the distance vector
            size_t distance_vector_offset = get_distance_vector_offset(rank1, right_side1, rank2, right_side2);

            snarl_tree_records[distance_vector_start + distance_vector_offset] = distance - 1;
        }

        size_t get_node_count() {
            return snarl_tree_records[record_offset + 1];
        }
        void set_node_count(size_t parent_pointer) {
            snarl_tree_records[record_offset + 1] = parent_pointer;
        }

        void set_child_record_pointer(size_t pointer) {
            snarl_tree_records[record_offset+6] = pointer;
        }

    };

    struct chain_record_t : snarl_tree_record_t {

        size_t record_size = 6;

        chain_record_t (size_t pointer){
            record_offset = pointer;
            RecordType type = static_cast<RecordType>(snarl_tree_records[record_offset]>>6);
            assert(type == CHAIN || 
                   type == DISTANCED_CHAIN);
        }

        size_t get_node_count() {
            return snarl_tree_records[record_offset + 1];
        }
        void set_node_count(size_t parent_pointer) {
            snarl_tree_records[record_offset + 1] = parent_pointer;
        }

        //Get the prefix sum value for this node (boundary node of a snarl in the chain)
        //pointer is a pointer into snarl_tree_records, to the beginning of the record for this node
        //So it'll point to the node id of the node we're looking at
        int64_t get_prefix_sum_value(size_t pointer) {
            return snarl_tree_records[pointer+1]-1; 
        }
        void set_prefix_sum_value(size_t pointer, int64_t value) {
            snarl_tree_records[pointer+1] = value + 1; 
        }
        int64_t get_forward_loop_value(size_t pointer) {
            return snarl_tree_records[pointer+2]-1; 
        }
        void set_forward_loop_value(size_t pointer, int64_t value) {
            snarl_tree_records[pointer+2] = value + 1; 
        }
        int64_t get_reverse_loop_value(size_t pointer) {
            return snarl_tree_records[pointer+3]-1; 
        }
        void set_reverse_loop_value(size_t pointer, int64_t value) {
            snarl_tree_records[pointer+3] = value + 1; 
        }

        //Get the distance between the given node sides (relative to the orientation of the chain)
        //Nodes represent a tuple of <pointer, rank, right_side, and length of the node>
        //This is the distance between the node sides, leaving the first and entering the second,
        //not including node lengths
        //TODO: I don't think we're allowing looping chains so I'm going to ignore them for now
        int64_t get_distance(tuple<size_t, size_t, bool, size_t> node1, 
                             tuple<size_t, size_t, bool, size_t> node2) {

            if (std::get<1>(node1) > std::get<1>(node2)) {
                //If the first node comes after the second in the chain, reverse them
                tuple<size_t,size_t, bool, size_t> tmp = node1;
                node1 = node2;
                node2 = tmp;
            }

            if (std::get<2(node1) && !std::get<2(node2)) {
                //Right of 1 and left of 2, so a simple forward traversal of the chain
                return get_prefix_sum_value(std::get<0>(node2) 
                     - get_prefix_sum_value(std::get<0>(node1)  
                     - std::get<3>(node1);
            } else if (std::get<2(node1) && std::get<2(node2)) {
                //Right side of 1 and right side of 2
                return get_prefix_sum_value(std::get<0>(node2) 
                     - get_prefix_sum_value(std::get<0>(node1)  
                     - std::get<3>(node1) + std::get<3>(node2) 
                     + get_forward_loop_value(std::get<0>(node2);
            } else if (!std::get<2(node1) && !std::get<2(node2)) {
                //Left side of 1 and left side of 2
                return get_prefix_sum_value(std::get<0>(node2) 
                     - get_prefix_sum_value(std::get<0>(node1)  
                     + get_reverse_loop_value(std::get<0>(node1);
            } else if (!std::get<2(node1) && std::get<2(node2)) {
                //Left side of 1 and right side of 2
                return get_prefix_sum_value(std::get<0>(node2) 
                     - get_prefix_sum_value(std::get<0>(node1)  
                     + get_reverse_loop_value(std::get<0>(node1)
                     + get_forward_loop_value(std::get<0>(node1) 
                     + std::get<3>(node1);
            }
        }
    };


    struct trivial_chain_record_t{
        //Struct for a chain record of a trivial chain, that's actually just a node
    };

    //Define what a net_handle_t is:
    //last 4 bits will be the connectivity of the handle (exactly one bit will be set), 
    //the rest will be a pointer into snarl_tree_records
    //
    handlegraph::net_handle_t get_net_handle(size_t pointer, ConnectivityType connectivity) {
        return (pointer << 4) & (size_t)connectivity;
    }

    ConnectivityType get_connectivity (handlegraph::net_handle_t net_handle) {
        size_t connectivity_as_int = net_handle & 15; //Get last 4 bits
        assert (connectivity_as_int <= 9);
        return static_cast<ConnectivityType>(connectivity_as_int);
    }

    snarl_tree_record_t get_snarl_tree_record(handlegraph::net_handle_t net_handle) {
        return snarl_tree_record_t((size_t)(net_handle >> 4)); 
    }
    snarl_tree_record_t get_node_record(handlegraph::net_handle_t net_handle) {
        return node_record_t((size_t)(net_handle >> 4)); 
    }
    snarl_tree_record_t get_snarl_record(handlegraph::net_handle_t net_handle) {
        return snarl_record_t((size_t)(net_handle >> 4)); 
    }
    snarl_tree_record_t get_chain_record(handlegraph::net_handle_t net_handle) {
        return chain_record_t((size_t)(net_handle >> 4)); 
    }
};

}

#endif
