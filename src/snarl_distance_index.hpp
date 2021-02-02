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

class MinimumDistanceIndex : public SnarlDecomposition {

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
    //
    //- The (single) child vector, listing children in snarls
    //  [child vector tag, (pointer to records) x N
    //  Each snarl will have a pointer into here, and will also know how many children it has

    int_vector<> snarl_tree_records;


public:

private:

    //A tag for defining what kind of record we're looking at
    //TODO: Maybe also add a tag for trivial snarls/chains, or for defining connectivity
    //TODO: Also need to put a bool or something somewhere that says whether we store distances
    enum RecordType {ROOT, NODE, SNARL, OVERSIZED_SNARL, CHAIN, CHILDREN}

    //Define a struct for interpreting each type of snarl tree node record (For node, snarl, chain)
    //TODO: This might be overkill but I want it to be easy to change what gets stored in the index
    //
    struct snarl_tree_record_t {


        //The offset of the start of this record in snarl_tree_records
        size_t record_offset;

        //What type of snarl tree node is this?
        //This will be the first value of any record
        RecordType type() {
            return snarl_tree_records[record_offset];
        }

        //How many things get stored in this record_t?
        //This is used to allocate memory for the node's record
        size_t record_size;

        //Get and set a pointer to this chain's parent
        size_t get_parent_record_pointer();
        void set_parent_record_pointer(size_t pointer);

        //Get and set the minimum length (distance from start to end, including boundaries for 
        //snarls and chains, just node length for nodes)
        size_t get_min_length();
        void set_min_length(size_t length);

        //Get and set this snarl's maximum length
        //This isn't actually a maximum, it's the maximum among minimum distance paths 
        //through each node in the snarl
        size_t get_max_length();
        void set_max_length(size_t length);

        //Get and set this snarl's rank in its parent
        size_t get_rank_in_parent();
        void set_rank_in_parent(size_t rank);
    };

    struct node_record_t : snarl_tree_record_t {

        size_t record_size = 4;

        node_record_t (size_t pointer){
            record_offset = pointer;
            assert(snarl_tree_records[record_offset] == NODE);
        }

        size_t get_parent_record_pointer() {
            return snarl_tree_records[record_offset + 1];
        }
        void set_parent_record_pointer(size_t parent_pointer) {
            snarl_tree_records[record_offset + 1] = parent_pointer;
        }

        size_t get_node_length() {
            return snarl_tree_records[record_offset + 2];
        }
        void set_node_length(size_t length) {
            snarl_tree_records[record_offset + 2] = length;
        }
        //TODO: Adding these to be consistent but unecessary?
        size_t get_min_length() {
            return snarl_tree_records[record_offset + 2];
        }
        void set_min_length(size_t length) {
            snarl_tree_records[record_offset + 2] = length;
        }
        size_t get_max_length() {
            return snarl_tree_records[record_offset + 2];
        }
        void set_max_length(size_t length) {
            snarl_tree_records[record_offset + 2] = length;
        }

        size_t get_rank_in_parent() {
            return snarl_tree_records[record_offset + 3];
        }
        void set_rank_in_parent(size_t rank) {
            snarl_tree_records[record_offset + 3] = rank;
        }

        size_t get_root_component() {
            return snarl_tree_records[record_offset + 4];
        }
        void set_root_component(size_t component) {
            snarl_tree_records[record_offset + 4] = component;
        }

    };

    struct snarl_record_t : snarl_tree_record_t {

        size_t record_size = 6;

        snarl_record_t (size_t pointer){
            record_offset = pointer;
            assert(snarl_tree_records[record_offset] == SNARL ||
                   snarl_tree_records[record_offset] == OVERSIZED_SNARL);
        }

        size_t get_node_count() {
            return snarl_tree_records[record_offset + 1];
        }
        void set_node_count(size_t parent_pointer) {
            snarl_tree_records[record_offset + 1] = parent_pointer;
        }

        size_t get_parent_record_pointer() {
            return snarl_tree_records[record_offset + 2];
        }
        void set_parent_record_pointer(size_t pointer) {
            snarl_tree_records[record_offset + 2] = pointer;
        }

        size_t get_min_length() {
            return snarl_tree_records[record_offset + 3];
        }
        void set_min_length(size_t length) {
            snarl_tree_records[record_offset + 3] = length;
        }

        size_t get_max_length() {
            return snarl_tree_records[record_offset + 4];
        }
        void set_max_length(size_t length) {
            snarl_tree_records[record_offset + 4] = length;
        }

        size_t get_rank_in_parent() {
            return snarl_tree_records[record_offset + 5];
        }
        void set_rank_in_parent(size_t rank) {
            snarl_tree_records[record_offset + 5] = rank;
        }

        child_record_t get_child_records(){
            return child_record_t(snarl_tree_records[record_offset+6], get_node_count());
        }
        void set_child_record_pointer(size_t pointer) {
            snarl_tree_records[record_offset+6] = pointer;
        }

    };

    struct chain_record_t : snarl_tree_record_t {

        size_t record_size = 5;

        chain_record_t (size_t pointer){
            record_offset = pointer;
            assert(snarl_tree_records[record_offset] == CHAIN);
        }

        size_t get_node_count() {
            return snarl_tree_records[record_offset + 1];
        }
        void set_node_count(size_t parent_pointer) {
            snarl_tree_records[record_offset + 1] = parent_pointer;
        }

        size_t get_parent_record_pointer() {
            return snarl_tree_records[record_offset + 2];
        }
        void set_parent_record_pointer(size_t pointer) {
            snarl_tree_records[record_offset + 2] = pointer;
        }

        size_t get_min_length() {
            return snarl_tree_records[record_offset + 3];
        }
        void set_min_length(size_t length) {
            snarl_tree_records[record_offset + 3] = length;
        }

        size_t get_max_length() {
            return snarl_tree_records[record_offset + 4];
        }
        void set_max_length(size_t length) {
            snarl_tree_records[record_offset + 4] = length;
        }

        size_t get_rank_in_parent() {
            return snarl_tree_records[record_offset + 5];
        }
        void set_rank_in_parent(size_t rank) {
            snarl_tree_records[record_offset + 5] = rank;
        }

    };

    struct child_record_t {
        size_t start_offset;
        size_t end_offset;

        child_record_t (size_t start, size_t count) {
            start_offset = start;
            end_offset = start + count;
        }
    };

    struct trivial_chain_record_t{
        //Struct for a chain record of a trivial chain, that's actually just a node
    };

};

}

#endif
