#include "zip_code.hpp"

//#define DEBUG_ZIPCODE

namespace vg{
using namespace std;

void ZipCode::fill_in_zipcode (const SnarlDistanceIndex& distance_index, const pos_t& pos, bool fill_in_decoder) {

    std::vector<net_handle_t> ancestors;
    net_handle_t current_handle = distance_index.get_node_net_handle(id(pos));

    //Put all ancestors of the node in a vector, starting from the node, and not including the root
    while (!distance_index.is_root(current_handle)) {
        ancestors.emplace_back(distance_index.start_end_traversal_of(current_handle));
        current_handle = distance_index.get_parent(current_handle);
    }


    //Now add the root-level snarl or chain
    if (distance_index.is_root_snarl(current_handle)) {
        //First thing is a snarl, so add the snarl's connected component number
        zipcode.add_value(0);
#ifdef DEBUG_ZIPCODE
        cerr << "Adding code for top-level snarl " << distance_index.net_handle_as_string(current_handle) << endl;
#endif
        zipcode.add_value(distance_index.get_connected_component_number(current_handle));
    } else {
        //First thing is a chain so add its connected component number and remove the chain from the stack
        zipcode.add_value(1);

        //If the root-level structure is actually a chain, then save the connected component number and take out
        //the chain from the stack
        //If the root-level structure is a trivial chain, then just store the node (as a chain, which will have the 
        //connected-component number as the rank in the snarl anyways)
        zipcode.add_value(distance_index.get_connected_component_number(ancestors.back()));
        if (ancestors.size() == 2 && distance_index.is_trivial_chain(ancestors.back())) {
#ifdef DEBUG_ZIPCODE
           cerr << "Adding code for top-level trivial chain" << endl;
#endif
            zipcode.add_value(distance_index.minimum_length(ancestors.back())+1);
            size_t connectivity = 0;
            if ( distance_index.is_externally_start_end_connected(ancestors.back())) {
                connectivity = connectivity | 1;
            }
            if ( distance_index.is_externally_start_start_connected(ancestors.back())) {
                connectivity = connectivity | 2;
            }
            if ( distance_index.is_externally_end_end_connected(ancestors.back())) {
                connectivity = connectivity | 4;
            }
 
            zipcode.add_value(connectivity);
            if (fill_in_decoder) {
                fill_in_full_decoder();
            }
            return;
        } else {
#ifdef DEBUG_ZIPCODE
            cerr << "Adding code for top-level chain" << endl;
#endif

            size_t component = distance_index.get_chain_component(distance_index.get_bound(ancestors.back(), true, false), true);
            component = component == std::numeric_limits<size_t>::max() ? 0 : component*2;
            bool is_looping_chain = distance_index.is_looping_chain(ancestors.back());
            if (is_looping_chain) {
                component += 1;
            }
            zipcode.add_value(component);

            size_t connectivity = 0;
            if (is_looping_chain) {
                //For a looping chain, the "connectivity" is the length of the last component
                size_t length = distance_index.chain_minimum_length(ancestors.back());
                zipcode.add_value(length == std::numeric_limits<size_t>::max() ? 0 : length+1);
            } else {
                //For a non-looping chain, it is actually the connectivity
                if ( distance_index.is_externally_start_end_connected(ancestors.back())) {
                    connectivity = connectivity | 1;
                }
                if ( distance_index.is_externally_start_start_connected(ancestors.back())) {
                    connectivity = connectivity | 2;
                }
                if ( distance_index.is_externally_end_end_connected(ancestors.back())) {
                    connectivity = connectivity | 4;
                }
                zipcode.add_value(connectivity);
            }
        }
        ancestors.pop_back();
    }

    //Go through the ancestors top (root) down and add them to the zip code
    //ancestors has everything but the root-level snarl/chain
    for (int i = ancestors.size()-1 ; i >= 0 ; i--) {
        net_handle_t current_ancestor = ancestors[i];
#ifdef DEBUG_ZIPCODE
        cerr << "Adding code for " << distance_index.net_handle_as_string(current_ancestor) << endl;
#endif
        if (distance_index.is_node(current_ancestor)) {
            node_code_t node_code = get_node_code(current_ancestor, distance_index);
            zipcode.add_value(node_code.get_raw_prefix_sum_or_identifier());
            zipcode.add_value(node_code.get_raw_length());
            zipcode.add_value(node_code.get_raw_is_reversed());
            zipcode.add_value(node_code.get_raw_chain_component());

        } else if (distance_index.is_chain(current_ancestor)) {
            chain_code_t chain_code = get_chain_code(current_ancestor, distance_index);
            zipcode.add_value(chain_code.get_raw_snarl_rank_or_identifier());
            zipcode.add_value(chain_code.get_raw_length());
            zipcode.add_value(chain_code.get_raw_last_component());
#ifdef DEBUG_ZIPCODE
                assert(to_add.size() == ZipCode::CHAIN_SIZE);
#endif
            if (distance_index.is_trivial_chain(current_ancestor)) {
                if (fill_in_decoder) {
                    fill_in_full_decoder();
                }
                return;
            }
        } else if (distance_index.is_regular_snarl(current_ancestor)) {
            snarl_code_t snarl_code = get_regular_snarl_code(current_ancestor, ancestors[i-1], distance_index); 
            zipcode.add_value(snarl_code.get_raw_code_type());
            zipcode.add_value(snarl_code.get_raw_prefix_sum_or_identifier());
            zipcode.add_value(snarl_code.get_raw_length());
            zipcode.add_value(snarl_code.get_raw_child_count());
            zipcode.add_value(snarl_code.get_raw_chain_component());
            zipcode.add_value(snarl_code.get_raw_is_reversed());
        } else {
#ifdef DEBUG_ZIPCODE
            assert(distance_index.is_snarl(current_ancestor));
#endif
            snarl_code_t snarl_code = get_irregular_snarl_code(current_ancestor, ancestors[i-1], distance_index); 
            zipcode.add_value(snarl_code.get_raw_code_type());
            zipcode.add_value(snarl_code.get_raw_prefix_sum_or_identifier());
            zipcode.add_value(snarl_code.get_raw_length());
            zipcode.add_value(snarl_code.get_raw_child_count());
            zipcode.add_value(snarl_code.get_raw_chain_component());
            zipcode.add_value(snarl_code.get_raw_record_offset());
            zipcode.add_value(snarl_code.get_raw_distance_start_left());
            zipcode.add_value(snarl_code.get_raw_distance_end_left());
            zipcode.add_value(snarl_code.get_raw_distance_start_right());
            zipcode.add_value(snarl_code.get_raw_distance_end_right());
        }
    }
    if (fill_in_decoder) {
        fill_in_full_decoder();
    }
}

std::vector<size_t> ZipCode::to_vector() const {
    return zipcode.to_vector();
}

void ZipCode::from_vector(const std::vector<size_t>& values) {
    zipcode.from_vector(values);
}


void ZipCode::fill_in_full_decoder() {
    if (byte_count() == 0 || finished_decoding) {
        //If the zipcode is empty
        return;
    }
    decoder.reserve(byte_count() / 4);
    bool done=false;
    while (!done) {
        done = fill_in_next_decoder();
    }
    finished_decoding = true;
}

bool ZipCode::fill_in_next_decoder() {
#ifdef DEBUG_ZIPCODE
    cerr << "Decode one more thing in the zipcode. Currently decoded " << decoder_length() << " things" << endl;
#endif
    if (finished_decoding) {
        return true;
    }
    
    //The zipcode may be partially or fully filled in already, so first
    //check to see how much has been filled in
    size_t zip_length = decoder_length();

    //Does the most recent thing in the zip_index point to a chain/node?
    bool previous_is_chain;

    size_t zip_index=0;
    size_t zip_value;

    if (zip_length == 0) {
        //If there is nothing in the decoder yet, then the first thing will start at 0
        for (size_t i = 0 ; i <= ZipCode::ROOT_IS_CHAIN_OFFSET ; i++) {
            std::tie(zip_value, zip_index) = zipcode.get_value_and_next_index(zip_index);
        }

        //Is the root a chain/node?
        previous_is_chain = zip_value;
        decoder.emplace_back(previous_is_chain, 0);

#ifdef DEBUG_ZIPCODE
cerr << "\tadding the root, which is a " << (previous_is_chain ? "chain or node" : "snarl") << endl;
#endif
        //There might be something else but we're done for now
        return false;
    } else if (zip_length == 1) {
        //If there is one thing in the zipcode
        previous_is_chain = decoder.back().is_chain;

        //If the top-level structure is a chain, it might actually be a node, in which case
        //the only other thing that got stored is the length
        if (previous_is_chain) {
            //Get to the end of the root chain
            assert(ZipCode::ROOT_CHAIN_SIZE==ZipCode::ROOT_NODE_SIZE);//This is true for now but all this will change if it isn't

            for (size_t i = 0 ; i < ZipCode::ROOT_NODE_SIZE ; i++) {
                std::tie(zip_value, zip_index) = zipcode.get_value_and_next_index(zip_index);
            }
            if (zip_index == std::numeric_limits<size_t>::max()) {
                //If the zip code ends here (after the length), then this was a node and we're done
#ifdef DEBUG_ZIPCODE
cerr << "\tThe last thing was a root-level node, so nothing else" << endl;
#endif
                finished_decoding = true;
                return true;
            } else {
                //Otherwise, check if this is a node or a snarl. If it is a node, then there are three things remaining
                size_t start_index = zip_index;

                //If it's a node, then there are three remaining things in the index 
                //If it were a snarl, then there are more than three things
                for (size_t i = 0 ; i < ZipCode::NODE_SIZE ; i++) {
                    std::tie(zip_value, zip_index) = zipcode.get_value_and_next_index(zip_index);
                }


                //Return the start of this thing, and true if it was a node
                decoder.emplace_back(zip_index == std::numeric_limits<size_t>::max(), start_index);
#ifdef DEBUG_ZIPCODE
                cerr << "\tAdding a " << (zip_index == std::numeric_limits<size_t>::max() ? "node" : "snarl") << endl;
#endif
                //If this was a node, then we're done so return true. Otherwise, it was a snarl to return false
                return zip_index == std::numeric_limits<size_t>::max();
            }
        } else {
            //Otherwise, the top-level thing is a snarl and the next thing is a chain 
            for (size_t i = 0 ; i < ZipCode::ROOT_SNARL_SIZE ; i++) {
                std::tie(zip_value, zip_index) = zipcode.get_value_and_next_index(zip_index);
            }
            decoder.emplace_back(!previous_is_chain, zip_index);
            return false;
        }
    } else {
        //If there was already stuff in the decoder, then figure out where the last thing
        //is and set values
        previous_is_chain = decoder.back().is_chain;
        zip_index = decoder.back().offset;
#ifdef DEBUG_ZIPCODE
        cerr << "Last thing was a " << (previous_is_chain ? "chain or node" : "snarl") << " starting at " << zip_index << endl;
#endif

        //get to the end of the current thing, add the next thing to the decoder and return

        if (previous_is_chain) {
            //If the current zip_index points to a chain, then either it points to a node, or to 
            //a chain that is followed by a node or snarl
            //The node is the shorter of the two, so if the zipcode ends after the node, then it was
            //a node and otherwise, it was an actual chain

            //This must be true in order for this to work
            assert(std::min(ZipCode::CHAIN_SIZE + ZipCode::REGULAR_SNARL_SIZE,
                            ZipCode::CHAIN_SIZE + ZipCode::IRREGULAR_SNARL_SIZE) > ZipCode::NODE_SIZE);

            //Get to the end of the "node". If it is the end of the zipcode, then it was a node
            //Otherwise, it was a snarl
            //The node could actually be a chain in a snarl, in which case the zipcode ends after the 
            //chain
            size_t check_zip_index = zip_index;
            for (size_t i = 0 ; i < std::min(ZipCode::CHAIN_SIZE, ZipCode::NODE_SIZE) ; i++) {
                check_zip_index = zipcode.get_value_and_next_index(check_zip_index).second;
            }
            //If the zipcode ends after a chain
            if (check_zip_index == std::numeric_limits<size_t>::max()) {
#ifdef DEBUG_ZIPCODE
                cerr << "\tThe last thing was a chain pretending to be a node so we're done" << endl;
#endif
                finished_decoding = true;
                return true;
            }
            //Now check if it was actually a real node
            for (size_t i = 0 ; i < std::max(ZipCode::NODE_SIZE, ZipCode::CHAIN_SIZE)
                                    - std::min(ZipCode::NODE_SIZE, ZipCode::CHAIN_SIZE); i++) {
                check_zip_index = zipcode.get_value_and_next_index(check_zip_index).second;
            }

            //This might be a node that is a child of the chain, in which case there is one
            //more thing in the zip code

            if (check_zip_index == std::numeric_limits<size_t>::max()) {
                //If the zip code ends here, then this was a node and we're done
                //This should never really happen since it would have returned true when
                //adding the node, but I'll leave in just in case someone calls this when they
                //shouldn't have
#ifdef DEBUG_ZIPCODE
                cerr << "\tThe last thing was a node so we're done" << endl;
#endif
                finished_decoding = true;
                return true;
            } else {
                //Otherwise, the last thing was a chain
                //Get to the end of the chain
                for (size_t i = 0 ; i < ZipCode::CHAIN_SIZE  ; i++) {
                    std::tie(zip_value, zip_index) = zipcode.get_value_and_next_index(zip_index);
                }

                //zip_index is now the start of the current thing that we want to add - the thing after the chain

                //The current thing can be either a snarl or a node. If it is a node, then the zipcode
                //ends after the node. If it is a snarl, then the shortest the remaining zipcocde can be
                //is the size of a snarl and a chain
                //This must be true in order for this to work
                assert(std::min(ZipCode::CHAIN_SIZE + ZipCode::REGULAR_SNARL_SIZE,
                                ZipCode::CHAIN_SIZE + ZipCode::IRREGULAR_SNARL_SIZE) > ZipCode::NODE_SIZE);

                //Check if the current thing is a node
                check_zip_index = zip_index;
                for (size_t i = 0 ; i < ZipCode::NODE_SIZE ; i++) {
                    check_zip_index = zipcode.get_value_and_next_index(check_zip_index).second;
                }

                //Return the start of this thing, and true if it was a node
                decoder.emplace_back(check_zip_index == std::numeric_limits<size_t>::max(), zip_index);
#ifdef DEBUG_ZIPCODE
                cerr << "\tAdd a " << (check_zip_index == std::numeric_limits<size_t>::max() ? "node" : "snarl") << endl;
#endif
                //If this was a node, then we're done so return true. Otherwise, it was a snarl to return false
                return check_zip_index == std::numeric_limits<size_t>::max();
            }
        } else {
            //If !previous_is_chain, then the current zip_index points to a snarl

            //The regular/irregular snarl tag
            for (size_t i = 0 ; i <= ZipCode::SNARL_IS_REGULAR_OFFSET  ; i++) {
                std::tie(zip_value, zip_index) = zipcode.get_value_and_next_index(zip_index);
            }

            if (zip_value == 1) {
#ifdef DEBUG_ZIPCODE
                cerr << "\tAdd a node child of a regular snarl" << endl;
#endif
                //Regular snarl, so 2 remaining things in the code
                for (size_t i = 0 ; i < ZipCode::REGULAR_SNARL_SIZE - ZipCode::SNARL_IS_REGULAR_OFFSET - 1; i++) {
                    std::tie(zip_value, zip_index) = zipcode.get_value_and_next_index(zip_index);
                }
                decoder.emplace_back(!previous_is_chain, zip_index);
                return false;
            } else {
#ifdef DEBUG_ZIPCODE
                cerr << "\tAdd the child of " << (decoder.size() == 2 ? "a top-level " : "an" ) << " irregular snarl" << endl;
#endif
                //If the decoder has two things in it (top-level chain and the current snarl), then this
                //is a top-level irregular snarl. Otherwise a normal irregular snarl
                size_t code_size = ZipCode::IRREGULAR_SNARL_SIZE;
                for (size_t i = 0 ; i < code_size - ZipCode::SNARL_IS_REGULAR_OFFSET - 1; i++) {
                    std::tie(zip_value, zip_index) = zipcode.get_value_and_next_index(zip_index);
                }
                decoder.emplace_back(!previous_is_chain, zip_index);
                return false;
            }
        }
    }    
}

size_t ZipCode::max_depth() const {
    return decoder_length()-1;

}

ZipCode::code_type_t ZipCode::get_code_type(const size_t& depth) const {

    //Now get the code type
    //A snarl is always a snarl. A chain could actually be a node
    if (depth == 0) {
        //If it is a root snarl/chain
        if (decoder[0].is_chain) {
            //If it says it's a chain, then it might be a chain or a node

            //If there is still only one thing in the decoder, then it's a node
            if (decoder_length() == 1) {
                return ZipCode::ROOT_NODE;
            } else {
                return ZipCode::ROOT_CHAIN;
            }
        } else {
            return ZipCode::ROOT_SNARL;
        }
    } else {
        if (decoder[depth].is_chain) {
            //is_chain so could be a chain or a node
            if (decoder[depth-1].is_chain) {
                //If the thing before this was also a chain, then it is a node
                return ZipCode::NODE;
            } else {
                //Otherwise it's a chain
                return ZipCode::CHAIN;
            }
        } else {
            //Definitely a snarl
            size_t zip_value;
            size_t zip_index = decoder[depth].offset;
            for (size_t i = 0 ; i <= ZipCode::SNARL_IS_REGULAR_OFFSET; i++) {
                std::tie(zip_value, zip_index) = zipcode.get_value_and_next_index(zip_index);
            }
            if (zip_value == 0) {
                return ZipCode::IRREGULAR_SNARL;
            } else if (zip_value == 1) {
                return ZipCode::REGULAR_SNARL;
            } else {
                return  ZipCode::CYCLIC_SNARL;
            }
        }
    }
}

size_t ZipCode::get_length(const size_t& depth, bool get_chain_component_length) const {

    if (depth == 0) {
        //If this is the root chain/snarl/node

        if (decoder_length() == 1) {
            //If the length is 1, then it's a node
            size_t zip_value;
            size_t zip_index = decoder[depth].offset;
            for (size_t i = 0 ; i <= ZipCode::ROOT_NODE_LENGTH_OFFSET ; i++) {
                std::tie(zip_value, zip_index) = zipcode.get_value_and_next_index(zip_index);
            }
            return zip_value == 0 ? std::numeric_limits<size_t>::max() : zip_value-1;

        } else {
            //Otherwise, if it is a looping chain then we stored the "chain component length"
            if (get_chain_component_length) {
                size_t zip_value;
                size_t zip_index = decoder[depth].offset;
                for (size_t i = 0 ; i <= ZipCode::ROOT_CHAIN_COMPONENT_COUNT_OFFSET ; i++) {
                    std::tie(zip_value, zip_index) = zipcode.get_value_and_next_index(zip_index);
                }
                if (zip_value % 2) {
                    //If it was a looping chain
                    for (size_t i = ZipCode::CHAIN_COMPONENT_COUNT_OFFSET+1 ; i <= ZipCode::ROOT_NODE_OR_CHAIN_CONNECTIVITY_OR_LENGTH_OFFSET; i++) {
                        std::tie(zip_value, zip_index) = zipcode.get_value_and_next_index(zip_index);
                    }
                    return zip_value == 0 ? std::numeric_limits<size_t>::max() : zip_value - 1;
                }
            }
            
            //Otherwise, we didn't store the length
            throw std::runtime_error("zipcodes don't store lengths of top-level chains or snarls. Do your zipcode, minimizer, and graph files match?");
        }
    } else if (decoder[depth].is_chain) {
        //If this is a chain/node

        //If this is a chain or a node, then the length will be the second thing
        size_t zip_value;
        size_t zip_index = decoder[depth].offset;

        for (size_t i = 0 ; i <= ZipCode::CHAIN_LENGTH_OFFSET ; i++) {
            std::tie(zip_value, zip_index) = zipcode.get_value_and_next_index(zip_index);
        }
        
        size_t len = zip_value == 0 ? std::numeric_limits<size_t>::max() : zip_value-1;
        if (get_chain_component_length || (depth != 0 && decoder[depth-1].is_chain)) {
            //If this is a node or we want the component length that got saved, return the actual saved value
            return len;
        } else {
            //If we want the length of the last component of the chain, check if it is a multicopmonent chain
            std::tie(zip_value, zip_index) = zipcode.get_value_and_next_index(zip_index);
            if (zip_value != 0) {
                //If this is a multicomponent (or looping chain, which also must be a multicomponent chain)
                return std::numeric_limits<size_t>::max();
            } else {
                return len;
            }
        }

    } else {
        //If this is a snarl

        size_t zip_value;
        size_t zip_index = decoder[depth].offset;

        for (size_t i = 0 ; i <= ZipCode::SNARL_LENGTH_OFFSET ; i++) {
            std::tie(zip_value, zip_index) = zipcode.get_value_and_next_index(zip_index);
        }

        return zip_value == 0 ? std::numeric_limits<size_t>::max() : zip_value-1;
    }
}

size_t ZipCode::get_rank_in_snarl(const size_t& depth) const {


    if (depth == 0) {
        //If this is the root chain/snarl/node
        throw std::runtime_error("zipcodes don't store ranks of top-level chains or snarls. Do your zipcode, minimizer, and graph files match?");

    } else if (decoder[depth].is_chain) {
        //If this is a chain/node

        if (decoder[depth-1].is_chain) {
            throw std::runtime_error("zipcodes trying to find the rank in snarl of a node in a chain. Do your zipcode, minimizer, and graph files match?");
        }

        size_t zip_value;
        size_t zip_index = decoder[depth].offset;
        for (size_t i = 0 ; i <= ZipCode::CHAIN_RANK_IN_SNARL_OFFSET ; i++) {
            std::tie(zip_value, zip_index) = zipcode.get_value_and_next_index(zip_index);
        }
        return zip_value;
    } else {
        //If this is a snarl
        throw std::runtime_error("zipcodes don't store snarl ranks for snarls. Do your zipcode, minimizer, and graph files match?");
    }
}

size_t ZipCode::get_snarl_child_count(const size_t& depth, const SnarlDistanceIndex* distance_index) const {


    if (depth == 0) {
        //TODO: This could be actually saved in the zipcode but I'll have to go to the distance index anyway
        assert(distance_index != nullptr);
        size_t child_count = 0;
        distance_index->for_each_child(get_net_handle(depth, distance_index), [&] (const net_handle_t& child) {
            child_count++;
        });
        return child_count;

    } else if (!decoder[depth].is_chain) {
        //If this is a snarl

        size_t zip_value;
        size_t zip_index = decoder[depth].offset;
        for (size_t i = 0 ; i <= ZipCode::SNARL_CHILD_COUNT_OFFSET ; i++) {
            std::tie(zip_value, zip_index) = zipcode.get_value_and_next_index(zip_index);
        }
        return zip_value;
    } else {
        //If this is not a snarl
        throw std::runtime_error("trying to get the snarl child count of a non-snarl zipcode. Do your zipcode, minimizer, and graph files match?");
    }
}

size_t ZipCode::get_offset_in_chain(const size_t& depth, const SnarlDistanceIndex* distance_index) const {


    if (depth == 0) {
        //If this is the root chain/snarl/node
        throw std::runtime_error("zipcodes don't have chain offsets for roots. Do your zipcode, minimizer, and graph files match?");

    } else if (decoder[depth].is_chain) {
        //If this is a chain/node

        if (!decoder[depth-1].is_chain) {
            throw std::runtime_error("zipcodes trying to find the offset in child of a snarl. Do your zipcode, minimizer, and graph files match?");
        }
        size_t zip_value;
        size_t zip_index = decoder[depth].offset;
        for (size_t i = 0 ; i <= ZipCode::NODE_OFFSET_OFFSET ; i++) {
            std::tie(zip_value, zip_index) = zipcode.get_value_and_next_index(zip_index);
        }

        return zip_value == 0 ? std::numeric_limits<size_t>::max() : zip_value-1;
    } else {
        //If this is a snarl

        size_t zip_value;
        size_t zip_index = decoder[depth].offset;
        for (size_t i = 0 ; i <= ZipCode::SNARL_OFFSET_IN_CHAIN_OFFSET ; i++) {
                std::tie(zip_value, zip_index) = zipcode.get_value_and_next_index(zip_index);
        }

        return zip_value == 0 ? std::numeric_limits<size_t>::max() : zip_value-1;
    }
}
size_t ZipCode::get_chain_component(const size_t& depth) const {


    if (depth == 0) {
        //If this is the root chain/snarl/node
        throw std::runtime_error("zipcodes don't have chain offsets for roots. Do your zipcode, minimizer, and graph files match?");

    } else if (decoder[depth].is_chain) {
        //If this is a chain/node

        if (!decoder[depth-1].is_chain) {
            throw std::runtime_error("zipcodes trying to find the offset in child of a snarl. Do your zipcode, minimizer, and graph files match?");
        }
        size_t zip_value;
        size_t zip_index = decoder[depth].offset;
        for (size_t i = 0 ; i <= ZipCode::NODE_CHAIN_COMPONENT_OFFSET ; i++) {
            std::tie(zip_value, zip_index) = zipcode.get_value_and_next_index(zip_index);
        }

        return zip_value;
    } else {
        //If this is a snarl

        size_t zip_value;
        size_t zip_index = decoder[depth].offset;
        for (size_t i = 0 ; i <= ZipCode::SNARL_CHAIN_COMPONENT_OFFSET ; i++) {
                std::tie(zip_value, zip_index) = zipcode.get_value_and_next_index(zip_index);
        }

        return zip_value;
    }
}

size_t ZipCode::get_last_chain_component(const size_t& depth, bool get_end) const {

    if (!decoder[depth].is_chain) {
        throw std::runtime_error("zipcodes trying to find the last chain component a snarl. Do your zipcode, minimizer, and graph files match?");
    }
    size_t zip_value;
    size_t zip_index = decoder[depth].offset;
    for (size_t i = 0 ; i <= ZipCode::CHAIN_COMPONENT_COUNT_OFFSET ; i++) {
        std::tie(zip_value, zip_index) = zipcode.get_value_and_next_index(zip_index);
    }
    if (zip_value % 2) {
        if (!get_end) {
            return 0;
        } else {
            zip_value -= 1;
        }
    }
    
    return zip_value / 2;
}

bool ZipCode::get_is_looping_chain(const size_t& depth) const {

    if (!decoder[depth].is_chain) {
        throw std::runtime_error("zipcodes trying to find the last chain component a snarl. Do your zipcode, minimizer, and graph files match?");
    }
    size_t zip_value;
    size_t zip_index = decoder[depth].offset;
    for (size_t i = 0 ; i <= ZipCode::CHAIN_COMPONENT_COUNT_OFFSET ; i++) {
        std::tie(zip_value, zip_index) = zipcode.get_value_and_next_index(zip_index);
    }
    return zip_value % 2;
}
bool ZipCode::get_is_reversed_in_parent(const size_t& depth) const {


    if (depth == 0) {
        //If this is the root chain/snarl/node
        return false;

    } else if (decoder[depth].is_chain) {
        //If this is a chain/node

        if (decoder[depth-1].is_chain) {
            //If the parent is a chain, then this is a node and we need to check its orientation

            size_t zip_value;
            size_t zip_index = decoder[depth].offset;
            for (size_t i = 0 ; i <= ZipCode::NODE_IS_REVERSED_OFFSET ; i++) {
                std::tie(zip_value, zip_index) = zipcode.get_value_and_next_index(zip_index);
            }
            return zip_value;
        } else {
            //If the parent is a snarl, then this might be a chain in a regular snarl
            size_t zip_value;
            size_t zip_index = decoder[depth-1].offset;
            //zip_value is true if the parent is a regular snarl
            for (size_t i = 0 ; i <= ZipCode::SNARL_IS_REGULAR_OFFSET ; i++) {
                std::tie(zip_value, zip_index) = zipcode.get_value_and_next_index(zip_index);
            }
            if (zip_value == 1) {
                //The parent is a regular snarl, which stores is_reversed for the child
                
                for (size_t i = 0 ; i <= ZipCode::REGULAR_SNARL_IS_REVERSED_OFFSET -
                                         ZipCode::SNARL_IS_REGULAR_OFFSET - 1 ; i++) {
                    std::tie(zip_value, zip_index) = zipcode.get_value_and_next_index(zip_index);
                }
                return zip_value;
            } else {
                //The parent is an irregular snarl, so it isn't reversed
                return false;
            }
        }
    } else {
        //If this is a snarl
        return false;
    }
}

net_handle_t ZipCode::get_net_handle(const size_t& depth, const SnarlDistanceIndex* distance_index) const {
    //get_net_handle_slow does the same thing so if this gets changed need to change that too


    if (depth == 0) {
        //If this is the root chain/snarl/node

        size_t zip_value, zip_index = 0;
        for (size_t i = 0 ; i <= ZipCode::ROOT_IDENTIFIER_OFFSET ; i++) {
            std::tie(zip_value, zip_index) = zipcode.get_value_and_next_index(zip_index);
        }
        return distance_index->get_handle_from_connected_component(zip_value);

    } else if (decoder[depth].is_chain) {
        //If this is a chain/node

        throw std::runtime_error("zipcodes trying to get a handle of a chain or node. Do your zipcode, minimizer, and graph files match?");
    } else {
        //If this is a snarl

        size_t zip_value; 
        size_t zip_index = decoder[depth].offset;
        //zip_value is is_regular_snarl
        for (size_t i = 0 ; i <= ZipCode::SNARL_IS_REGULAR_OFFSET ; i++) {
            std::tie(zip_value, zip_index) = zipcode.get_value_and_next_index(zip_index);
        }
        if (zip_value == 1) {
            //If this is a regular snarl

            throw std::runtime_error("zipcodes trying to get a handle of a regular snarl. Do your zipcode, minimizer, and graph files match?");
        } else {
            //Irregular snarl

            //zip_value is distance index offset
            for (size_t i = 0 ; i <= ZipCode::IRREGULAR_SNARL_RECORD_OFFSET -
                                     ZipCode::SNARL_IS_REGULAR_OFFSET-1 ; i++) {
                std::tie(zip_value, zip_index) = zipcode.get_value_and_next_index(zip_index);
            }
            net_handle_t snarl_handle = distance_index->get_net_handle_from_values(zip_value, SnarlDistanceIndex::START_END, SnarlDistanceIndex::SNARL_HANDLE);
            return snarl_handle;
        }
    }
}

net_handle_t ZipCode::get_net_handle_slow(nid_t id, const size_t& depth, const SnarlDistanceIndex* distance_index, const net_handle_t* child_handle) const {
    //This is just copying get_net_handle except adding a slower version for the things we don't remember

    if (depth == 0) {
        //If this is the root chain/snarl/node

        size_t zip_value, zip_index = 0;
        for (size_t i = 0 ; i <= ZipCode::ROOT_IDENTIFIER_OFFSET ; i++) {
            std::tie(zip_value, zip_index) = zipcode.get_value_and_next_index(zip_index);
        }
        return distance_index->get_handle_from_connected_component(zip_value);

    } else if (decoder[depth].is_chain) {
        //If this is a chain/node
        if (child_handle != nullptr) {
            return distance_index->get_parent(*child_handle);
        }

        net_handle_t n = distance_index->get_node_net_handle(id);
        for (size_t d = max_depth() ; d > depth ; d--) {
            n = distance_index->get_parent(n);
            if (distance_index->is_trivial_chain(n)){
                n = distance_index->get_parent(n);
            }
        }
        return n;
    } else {
        //If this is a snarl

        size_t zip_value; 
        size_t zip_index = decoder[depth].offset;
        //zip_value is is_regular_snarl
        for (size_t i = 0 ; i <= ZipCode::SNARL_IS_REGULAR_OFFSET ; i++) {
            std::tie(zip_value, zip_index) = zipcode.get_value_and_next_index(zip_index);
        }
        if (zip_value == 1) {
            //If this is a regular snarl
            if (child_handle != nullptr) {
                return distance_index->get_parent(*child_handle);
            }

            net_handle_t n = distance_index->get_node_net_handle(id);
            for (size_t d = max_depth() ; d > depth ; d--) {
                n = distance_index->get_parent(n);
                if (distance_index->is_trivial_chain(n)){
                    n = distance_index->get_parent(n);
                }
            }
            return n;
        } else {
            //Irregular snarl

            //zip_value is distance index offset
            for (size_t i = 0 ; i <= ZipCode::IRREGULAR_SNARL_RECORD_OFFSET -
                                     ZipCode::SNARL_IS_REGULAR_OFFSET-1 ; i++) {
                std::tie(zip_value, zip_index) = zipcode.get_value_and_next_index(zip_index);
            }
            net_handle_t snarl_handle = distance_index->get_net_handle_from_values(zip_value, SnarlDistanceIndex::START_END, SnarlDistanceIndex::SNARL_HANDLE);
            return snarl_handle;
        }
    }
}


size_t ZipCode::get_distance_index_address(const size_t& depth) const {


    if (depth == 0) {
        //If this is the root chain/snarl/node

        size_t zip_value, zip_index = 0;
        for (size_t i = 0 ; i <= ZipCode::ROOT_IDENTIFIER_OFFSET ; i++) {
            std::tie(zip_value, zip_index) = zipcode.get_value_and_next_index(zip_index);
        }
        return zip_value;

    } else if (decoder[depth].is_chain) {
        //If this is a chain/node

        throw std::runtime_error("zipcodes trying to get a handle of a chain or node. Do your zipcode, minimizer, and graph files match?");
    } else {
        //If this is a snarl

        size_t zip_value;
        size_t zip_index = decoder[depth].offset;
        //zip_value is is_regular_snarl
        for (size_t i = 0 ; i <= ZipCode::SNARL_IS_REGULAR_OFFSET ; i++) {
            std::tie(zip_value, zip_index) = zipcode.get_value_and_next_index(zip_index);
        }
        if (zip_value == 1) {
            //If this is a regular snarl

            throw std::runtime_error("zipcodes trying to get a handle of a regular snarl. Do your zipcode, minimizer, and graph files match?");
        } else {
            //Irregular snarl

            //zip_value is distance index offset
            for (size_t i = 0 ; i <= ZipCode::IRREGULAR_SNARL_RECORD_OFFSET -
                                     ZipCode::SNARL_IS_REGULAR_OFFSET-1 ; i++) {
                std::tie(zip_value, zip_index) = zipcode.get_value_and_next_index(zip_index);
            }
            return zip_value;
        }
    }
}
size_t ZipCode::get_distance_to_snarl_bound(const size_t& depth, bool snarl_start, bool left_side) const {

#ifdef DEBUG_ZIPCODE
    assert(depth > 0);
    assert((get_code_type(depth-1) == ZipCode::IRREGULAR_SNARL || get_code_type(depth-1) == ZipCode::REGULAR_SNARL || get_code_type(depth-1) == ZipCode::CYCLIC_SNARL)); 
#endif
     size_t zip_value;
     size_t zip_index = decoder[depth-1].offset;
     //zip_value is 1 if the parent is a regular snarl
     for (size_t i = 0 ; i <= ZipCode::SNARL_IS_REGULAR_OFFSET ; i++) {
         std::tie(zip_value, zip_index) = zipcode.get_value_and_next_index(zip_index);
     }
     if (zip_value == 1) {
         //The parent is a regular snarl, which stores is_reversed for the child
         for (size_t i = 0 ; i <= ZipCode::REGULAR_SNARL_IS_REVERSED_OFFSET -
                                  ZipCode::SNARL_IS_REGULAR_OFFSET - 1 ; i++) {
             std::tie(zip_value, zip_index) = zipcode.get_value_and_next_index(zip_index);
         }
         //Zip value is true if the child is reversed

         if ((snarl_start && left_side) || (!snarl_start && !left_side)) {
             return zip_value ? std::numeric_limits<size_t>::max() : 0;
         } else {
             assert((snarl_start && !left_side) || (!snarl_start && left_side));
             return zip_value ? 0 : std::numeric_limits<size_t>::max();
         }
     } else {
        //If the parent is an irregular snarl (or cyclic, which is the same), get the saved value
        size_t distance_offset;
        if (snarl_start && left_side) {
            distance_offset = ZipCode::IRREGULAR_SNARL_DISTANCE_LEFT_START_OFFSET;
        } else if (snarl_start && !left_side) {
            distance_offset = ZipCode::IRREGULAR_SNARL_DISTANCE_RIGHT_START_OFFSET;
        } else if (!snarl_start && left_side) {
            distance_offset = ZipCode::IRREGULAR_SNARL_DISTANCE_LEFT_END_OFFSET;
        } else {
            distance_offset = ZipCode::IRREGULAR_SNARL_DISTANCE_RIGHT_END_OFFSET;
        }
        for (size_t i = 0 ; i <= distance_offset - ZipCode::SNARL_IS_REGULAR_OFFSET -1 ; i++) {
            std::tie(zip_value, zip_index) = zipcode.get_value_and_next_index(zip_index);
        }
        return zip_value == 0 ? std::numeric_limits<size_t>::max() : zip_value - 1;
     }
}

bool ZipCode::is_externally_start_end_connected (const size_t& depth) const {
    assert(depth == 0);
    assert(decoder[0].is_chain);
    assert(CHAIN_COMPONENT_COUNT_OFFSET < ROOT_NODE_OR_CHAIN_CONNECTIVITY_OR_LENGTH_OFFSET);
    size_t zip_value;
    size_t zip_index = decoder[depth].offset;
    for (size_t i = 0 ; i <= ZipCode::CHAIN_COMPONENT_COUNT_OFFSET ; i++) {
        std::tie(zip_value, zip_index) = zipcode.get_value_and_next_index(zip_index);
    }
    if (zip_value % 2) {
        //If it is a looping chain, then it is technically start-end connected
        return true;
    }
    for (size_t i = ZipCode::CHAIN_COMPONENT_COUNT_OFFSET+1 ; i <= ZipCode::ROOT_NODE_OR_CHAIN_CONNECTIVITY_OR_LENGTH_OFFSET; i++) {
        std::tie(zip_value, zip_index) = zipcode.get_value_and_next_index(zip_index);
    }
    return (zip_value & 1) != 0;
}
bool ZipCode::is_externally_start_start_connected (const size_t& depth) const {
    assert(depth == 0);
    assert(decoder[0].is_chain);
    assert(CHAIN_COMPONENT_COUNT_OFFSET < ROOT_NODE_OR_CHAIN_CONNECTIVITY_OR_LENGTH_OFFSET);
    size_t zip_value;
    size_t zip_index = decoder[depth].offset;
    for (size_t i = 0 ; i <= ZipCode::CHAIN_COMPONENT_COUNT_OFFSET ; i++) {
        std::tie(zip_value, zip_index) = zipcode.get_value_and_next_index(zip_index);
    }
    if (zip_value % 2) {
        //If it is a looping chain, then it can't be start-start connected
        return false;
    }
    for (size_t i = ZipCode::CHAIN_COMPONENT_COUNT_OFFSET+1 ; i <= ZipCode::ROOT_NODE_OR_CHAIN_CONNECTIVITY_OR_LENGTH_OFFSET; i++) {
        std::tie(zip_value, zip_index) = zipcode.get_value_and_next_index(zip_index);
    }
    return (zip_value & 2) != 0;
}
bool ZipCode::is_externally_end_end_connected (const size_t& depth) const {
    assert(depth == 0);
    assert(decoder[0].is_chain);
    assert(CHAIN_COMPONENT_COUNT_OFFSET < ROOT_NODE_OR_CHAIN_CONNECTIVITY_OR_LENGTH_OFFSET);
    size_t zip_value;
    size_t zip_index = decoder[depth].offset;
    for (size_t i = 0 ; i <= ZipCode::CHAIN_COMPONENT_COUNT_OFFSET ; i++) {
        std::tie(zip_value, zip_index) = zipcode.get_value_and_next_index(zip_index);
    }
    if (zip_value % 2) {
        //If it is a looping chain, then it can't be end-end connected
        return false;
    }
    for (size_t i = ZipCode::CHAIN_COMPONENT_COUNT_OFFSET+1 ; i <= ZipCode::ROOT_NODE_OR_CHAIN_CONNECTIVITY_OR_LENGTH_OFFSET; i++) {
        std::tie(zip_value, zip_index) = zipcode.get_value_and_next_index(zip_index);
    }
    return (zip_value & 4) != 0;
}

const bool ZipCode::is_equal(const ZipCode& zip1, const ZipCode& zip2,
                                        const size_t& depth) {

    if (zip1.max_depth() < depth && zip2.max_depth() < depth ) {
        return false;
    }

    //First, check if the code types are the same
    ZipCode::code_type_t type1 = zip1.get_code_type(depth);
    ZipCode::code_type_t type2 = zip2.get_code_type(depth);
    if (type1 != type2) {
        return false;
    }

    if (type1 == ZipCode::ROOT_NODE || type1 == ZipCode::ROOT_CHAIN || type1 == ZipCode::ROOT_SNARL || type1 == ZipCode::IRREGULAR_SNARL || type1 == ZipCode::CYCLIC_SNARL ) {
        //If the codes are for root-structures or irregular/cyclic snarls, just check if the 
        //connected component numbers are the same
        return zip1.get_distance_index_address(depth) == zip2.get_distance_index_address(depth);
    } else {
        //Check the parent type. If the parent is a snarl, then check rank. If it's a chain,
        //then check the prefix sum
        if (zip1.get_code_type(depth-1) == ZipCode::REGULAR_SNARL ||
            zip1.get_code_type(depth-1) == ZipCode::IRREGULAR_SNARL ||
            zip1.get_code_type(depth-1) == ZipCode::CYCLIC_SNARL ||
            zip1.get_code_type(depth-1) == ZipCode::ROOT_SNARL) {
            //If the parent is a snarl, then check the rank
            return zip1.get_rank_in_snarl(depth) == zip2.get_rank_in_snarl(depth);
        } else {
            //Otherwise, check the offset in the chain
            //Since the type is the same, this is sufficient
            return zip1.get_offset_in_chain(depth) == zip2.get_offset_in_chain(depth);
        }
    }
}

void ZipCode::dump(std::ostream& out) const {
    std::vector<size_t> numbers = to_vector();
    // Print out the numbers in a way that is easy to copy-paste as a vector literal.
    out << "<zipcode {";
    for (size_t i = 0; i < numbers.size(); i++) {
        out << numbers[i];
        if (i + 1 < numbers.size()) {
            out << ", ";
        }
    }
    out << "}>";
}

std::ostream& operator<<(std::ostream& out, const ZipCode& zip) {
    return out << "<zipcode {" << zip.get_identifier(zip.max_depth())<< "}>";
}


ZipCode::node_code_t ZipCode::get_node_code(const net_handle_t& node, const SnarlDistanceIndex& distance_index) {
#ifdef DEBUG_ZIPCODE
    assert(!distance_index.is_trivial_chain(node));
    assert((distance_index.is_chain(distance_index.get_parent(node)) || distance_index.is_root(distance_index.get_parent(node))));
#endif
    //Node code is: offset in chain, length, is reversed
    node_code_t node_code;
    //Assume this node is in a regular chain
    node_code.set_prefix_sum_or_identifier(distance_index.get_prefix_sum_value(node)); 

    node_code.set_length(distance_index.minimum_length(node));

    node_code.set_is_reversed(distance_index.is_reversed_in_parent(node));
    node_code.set_chain_component(distance_index.get_chain_component(node));

    return node_code;

}
ZipCode::chain_code_t ZipCode::get_chain_code(const net_handle_t& chain, const SnarlDistanceIndex& distance_index) {
    //Chain code is: rank in snarl, length
    chain_code_t chain_code;
    chain_code.set_snarl_rank_or_identifier(distance_index.get_rank_in_parent(chain));

    bool is_trivial = distance_index.is_trivial_chain(chain) ;

    chain_code.set_length(is_trivial ? distance_index.minimum_length(chain) : distance_index.chain_minimum_length(chain));

    bool is_looping_chain(is_trivial ? false : distance_index.is_looping_chain(chain));
    size_t component = is_trivial
                       ? 0 
                       : distance_index.get_chain_component(distance_index.get_bound(chain, true, false), true);
    chain_code.set_last_component(component, is_looping_chain);

    return chain_code;

}
ZipCode::snarl_code_t ZipCode::get_regular_snarl_code(const net_handle_t& snarl, const net_handle_t& snarl_child, const SnarlDistanceIndex& distance_index) {
    //Regular snarl code is 1, offset in chain, length, is reversed
    snarl_code_t snarl_code;

    //Tag to say that it's a regular snarl
    snarl_code.set_code_type(1);

    //The number of children
    size_t child_count = 0;
    distance_index.for_each_child(snarl, [&] (const net_handle_t& child) {
        child_count++;
    });
    snarl_code.set_child_count(child_count);

    //Chain prefix sum value for the start of the snarl, which is the prefix sum of the start node + length of the start node
    net_handle_t start_node = distance_index.get_node_from_sentinel(distance_index.get_bound(snarl, false, false));

    snarl_code.set_prefix_sum_or_identifier(SnarlDistanceIndex::sum(distance_index.get_prefix_sum_value(start_node), distance_index.minimum_length(start_node)));

    snarl_code.set_chain_component(distance_index.get_chain_component(start_node));

    //Length of the snarl
    snarl_code.set_length(distance_index.minimum_length(snarl));

    //Is the child of the snarl reversed in the snarl
#ifdef DEBUG_ZIPCODE
    assert(distance_index.is_chain(snarl_child));
#endif
    snarl_code.set_is_reversed((distance_index.distance_in_parent(snarl, 
                                                        distance_index.get_bound(snarl, false, true),
                                                        distance_index.flip(distance_index.canonical(snarl_child))) != 0));

    return snarl_code;

}
ZipCode::snarl_code_t ZipCode::get_irregular_snarl_code(const net_handle_t& snarl, const net_handle_t& snarl_child, 
                                                 const SnarlDistanceIndex& distance_index) {
    snarl_code_t snarl_code;

    //Tag to say that it's an irregular snarl
    snarl_code.set_code_type(distance_index.is_dag(snarl) ? 0 : 2);

    //The number of children
    size_t child_count = 0;
    distance_index.for_each_child(snarl, [&] (const net_handle_t& child) {
        child_count++;
    });
    snarl_code.set_child_count(child_count);

    //Chain prefix sum value for the start of the snarl, which is the prefix sum of the start node + length of the start node
    net_handle_t start_node = distance_index.get_node_from_sentinel(distance_index.get_bound(snarl, false, false));

    snarl_code.set_prefix_sum_or_identifier(SnarlDistanceIndex::sum(distance_index.get_prefix_sum_value(start_node), distance_index.minimum_length(start_node)));


    snarl_code.set_chain_component(distance_index.get_chain_component(start_node) );

    //Length of the snarl
    snarl_code.set_length(distance_index.minimum_length(snarl));


    //Record offset to look up distances in the index later
    snarl_code.set_record_offset(distance_index.get_record_offset(snarl));

    snarl_code.set_distance_start_left(distance_index.distance_to_parent_bound(snarl, true, distance_index.flip(snarl_child)));
    snarl_code.set_distance_end_left(distance_index.distance_to_parent_bound(snarl, false, distance_index.flip(snarl_child)));
    snarl_code.set_distance_start_right(distance_index.distance_to_parent_bound(snarl, true, snarl_child));
    snarl_code.set_distance_end_right(distance_index.distance_to_parent_bound(snarl, false, snarl_child));


    return snarl_code;
}

ZipCode::node_code_t ZipCode::unpack_node_code(size_t zipcode_level) const {
    node_code_t node_code;
    if (zipcode_level == 0) {
        throw std::runtime_error("error: Unpacking a root node. Use a chain instead");
    } else { 

        size_t zip_index = decoder[zipcode_level].offset;
        size_t zip_value;
        //Prefix sum
        std::tie(zip_value, zip_index) = zipcode.get_value_and_next_index(zip_index);
        node_code.set_raw_prefix_sum_or_identifier(zip_value);
        //Length
        std::tie(zip_value, zip_index) = zipcode.get_value_and_next_index(zip_index);
        node_code.set_raw_length(zip_value);

        //Is reversed
        std::tie(zip_value, zip_index) = zipcode.get_value_and_next_index(zip_index);
        node_code.set_raw_is_reversed(zip_value);
        //Chain component
        std::tie(zip_value, zip_index) = zipcode.get_value_and_next_index(zip_index);
        node_code.set_raw_chain_component(zip_value);
    }
    return node_code;

}

ZipCode::chain_code_t ZipCode::unpack_chain_code(size_t zipcode_level) const {
    chain_code_t chain_code;
    size_t zip_index = decoder[zipcode_level].offset;
    size_t zip_value;
    if (zipcode_level == 0 && decoder.size() == 1) {
        //Root node
        //is_chain
        std::tie(zip_value, zip_index) = zipcode.get_value_and_next_index(zip_index);
        //Identifier
        std::tie(zip_value, zip_index) = zipcode.get_value_and_next_index(zip_index);
        chain_code.set_raw_snarl_rank_or_identifier(zip_value);

        //Node length
        std::tie(zip_value, zip_index) = zipcode.get_value_and_next_index(zip_index);
        chain_code.set_raw_length(zip_value);

        //Connectivity
        std::tie(zip_value, zip_index) = zipcode.get_value_and_next_index(zip_index);
        chain_code.set_raw_connectivity (zip_value);

        //No component
        chain_code.set_last_component(0, false);

    } else if (zipcode_level == 0) {
        //Root chain
        //is_chain
        std::tie(zip_value, zip_index) = zipcode.get_value_and_next_index(zip_index);
        //Identifier
        std::tie(zip_value, zip_index) = zipcode.get_value_and_next_index(zip_index);
        chain_code.set_raw_snarl_rank_or_identifier(zip_value);

        //Component count
        std::tie(zip_value, zip_index) = zipcode.get_value_and_next_index(zip_index);
        chain_code.set_raw_last_component(zip_value);

        //Connectivity
        std::tie(zip_value, zip_index) = zipcode.get_value_and_next_index(zip_index);
        chain_code.set_raw_connectivity (zip_value);

        //No Node length
        chain_code.set_length(std::numeric_limits<size_t>::max());
    } else {
        //Nested chain
        //Rank in snarl
        std::tie(zip_value, zip_index) = zipcode.get_value_and_next_index(zip_index);
        chain_code.set_raw_snarl_rank_or_identifier(zip_value);

        //Node length
        std::tie(zip_value, zip_index) = zipcode.get_value_and_next_index(zip_index);
        chain_code.set_raw_length(zip_value);


        //Component count
        std::tie(zip_value, zip_index) = zipcode.get_value_and_next_index(zip_index);
        chain_code.set_raw_last_component(zip_value);

        //No connectivity
        chain_code.set_connectivity (0);

    } 

    return chain_code;
}

ZipCode::snarl_code_t ZipCode::unpack_snarl_code(size_t zipcode_level) const {
    snarl_code_t snarl_code;
    size_t zip_index = decoder[zipcode_level].offset;
    size_t zip_value;
    if (zipcode_level == 0) {
        //Root snarl
        //is_chain
        std::tie(zip_value, zip_index) = zipcode.get_value_and_next_index(zip_index);
        //Identifier
        std::tie(zip_value, zip_index) = zipcode.get_value_and_next_index(zip_index);
        snarl_code.set_raw_prefix_sum_or_identifier(zip_value);

        //Nothing else gets stored so set everything else to inf
        snarl_code.set_length(std::numeric_limits<size_t>::max());
        snarl_code.set_distance_start_left(std::numeric_limits<size_t>::max());
        snarl_code.set_distance_start_right(std::numeric_limits<size_t>::max());
        snarl_code.set_distance_end_left(std::numeric_limits<size_t>::max());
        snarl_code.set_distance_end_right(std::numeric_limits<size_t>::max());
        snarl_code.set_record_offset(std::numeric_limits<size_t>::max());
        snarl_code.set_child_count(std::numeric_limits<size_t>::max());
        snarl_code.set_chain_component(std::numeric_limits<size_t>::max());
        snarl_code.set_code_type(std::numeric_limits<size_t>::max());

    } else {
        //Nested snarl

        //Snarl is regular
        std::tie(zip_value, zip_index) = zipcode.get_value_and_next_index(zip_index);
        snarl_code.set_raw_code_type(zip_value);

        //Offset in chain
        std::tie(zip_value, zip_index) = zipcode.get_value_and_next_index(zip_index);
        snarl_code.set_raw_prefix_sum_or_identifier(zip_value);

        //Length
        std::tie(zip_value, zip_index) = zipcode.get_value_and_next_index(zip_index);
        snarl_code.set_raw_length(zip_value);

        //Child count 
        std::tie(zip_value, zip_index) = zipcode.get_value_and_next_index(zip_index);
        snarl_code.set_raw_child_count(zip_value);

        //Chain component
        std::tie(zip_value, zip_index) = zipcode.get_value_and_next_index(zip_index);
        snarl_code.set_raw_chain_component(zip_value);

        if (snarl_code.get_code_type() == 1) {
            //Regular snarl

            //Is-reversed
            std::tie(zip_value, zip_index) = zipcode.get_value_and_next_index(zip_index);
            snarl_code.set_raw_is_reversed(zip_value);
        } else {
            //Irregular/cyclic snarl

            //Record offset
            std::tie(zip_value, zip_index) = zipcode.get_value_and_next_index(zip_index);
            snarl_code.set_raw_record_offset(zip_value);

            //distance left start
            std::tie(zip_value, zip_index) = zipcode.get_value_and_next_index(zip_index);
            snarl_code.set_raw_distance_start_left(zip_value);

            //distance left end
            std::tie(zip_value, zip_index) = zipcode.get_value_and_next_index(zip_index);
            snarl_code.set_raw_distance_end_left(zip_value);

            //distance right start
            std::tie(zip_value, zip_index) = zipcode.get_value_and_next_index(zip_index);
            snarl_code.set_raw_distance_start_right(zip_value);

            //Distance right end
            std::tie(zip_value, zip_index) = zipcode.get_value_and_next_index(zip_index);
            snarl_code.set_raw_distance_end_right(zip_value);
        }

    }
    return snarl_code;
}



size_t ZipCode::minimum_distance_between(ZipCode& zip1, const pos_t& pos1,   
    ZipCode& zip2, const pos_t& pos2, const SnarlDistanceIndex& distance_index,
    size_t distance_limit, bool undirected_distance, const HandleGraph* graph){


#ifdef DEBUG_ZIPCODE
//Make sure that the zip codes actually correspond to the positions
    ZipCode check_zip1;
    check_zip1.fill_in_zipcode(distance_index, pos1);
    assert(zip1 == check_zip1);

    ZipCode check_zip2;
    check_zip2.fill_in_zipcode(distance_index, pos2);
    assert(zip2 == check_zip2);

    cerr << endl << "Minimum distance between " << pos1 << " and " << pos2 << " using zipcodes" << endl;
    cerr << "Ancestors for " << pos1 << endl;
    net_handle_t net1 = distance_index.get_node_net_handle(id(pos1));
    while ( !distance_index.is_root(net1)){
        cerr << "\t" << distance_index.net_handle_as_string(net1) << endl;
        net1 = distance_index.get_parent(net1);
    }
    cerr << "\t" << distance_index.net_handle_as_string(net1) << endl;
    cerr << "Ancestors for " << pos2 << endl;
    net_handle_t net2 = distance_index.get_node_net_handle(id(pos2));
    while ( !distance_index.is_root(net2)){
        cerr << "\t" << distance_index.net_handle_as_string(net2) << endl;
        net2 = distance_index.get_parent(net2);
    }
    cerr << "\t" << distance_index.net_handle_as_string(net2) << endl;
#endif

    //Helper function to update the distances to the ends of the parent
    //distance_start and distance_end get updated
    auto update_distances_to_ends_of_parent = [&] (ZipCode& zip, const size_t& child_depth, 
                                            size_t& distance_to_start, size_t& distance_to_end) {
#ifdef DEBUG_ZIPCODE
        cerr << "Update distance to ends of parent at depth " << child_depth << endl;
#endif
        //The distances from the start/end of current child to the start/end(left/right) of the parent
        size_t distance_start_left = std::numeric_limits<size_t>::max();
        size_t distance_start_right = std::numeric_limits<size_t>::max();
        size_t distance_end_left = std::numeric_limits<size_t>::max();
        size_t distance_end_right = std::numeric_limits<size_t>::max();

        code_type_t parent_type = zip.get_code_type(child_depth-1);

        if (parent_type == IRREGULAR_SNARL || parent_type == CYCLIC_SNARL) {
            //If the parent is an irregular snarl
            net_handle_t parent_handle = zip.get_net_handle(child_depth-1, &distance_index);
            size_t child_rank = zip.get_rank_in_snarl(child_depth);
            distance_start_left = distance_index.distance_in_snarl(parent_handle, 
                                    child_rank, false, 0, false, graph);
            distance_start_right = distance_index.distance_in_snarl(parent_handle, 
                                    child_rank, false, 1, false, graph);
            distance_end_right = distance_index.distance_in_snarl(parent_handle, 
                                    child_rank, true, 1, false, graph);
            distance_end_left = distance_index.distance_in_snarl(parent_handle, 
                                    child_rank, true, 0, false, graph);
#ifdef DEBUG_ZIPCODE
            cerr << "Distances to parent irregular snarl: " << distance_start_left << " " << distance_start_right << " " << distance_end_left << " " << distance_end_right << endl;
#endif
        } else if (parent_type == REGULAR_SNARL) {
            //If it's a regular snarl, then the distances to the ends are either 0 or inf
            //For a regular snarl, the snarl stores if the child was reversed, rather than the child
            if (zip.get_is_reversed_in_parent(child_depth)) {
                distance_start_left = std::numeric_limits<size_t>::max();
                distance_start_right = 0;
                distance_end_right = std::numeric_limits<size_t>::max();
                distance_end_left = 0;
            } else {
                distance_start_left = 0;
                distance_start_right = std::numeric_limits<size_t>::max();
                distance_end_right = 0;
                distance_end_left = std::numeric_limits<size_t>::max();
            }
#ifdef DEBUG_ZIPCODE
            cerr << "Distances to parent regular snarl: " << distance_start_left << " " << distance_start_right << " " << distance_end_left << " " << distance_end_right << endl;
#endif
        } else if (parent_type == CHAIN) {
            if (zip.get_code_type(child_depth) == NODE && 
                zip.get_is_reversed_in_parent(child_depth)){ 
                //If this is reversed in the chain

                distance_start_left = std::numeric_limits<size_t>::max();
                distance_end_right = std::numeric_limits<size_t>::max();
                //Prefix sum of the child
                distance_end_left = zip.get_offset_in_chain(child_depth, &distance_index);
                //Length of the chain - prefix sum of the child - length of the child
                distance_start_right = SnarlDistanceIndex::minus(SnarlDistanceIndex::minus(
                        zip.get_length(child_depth-1, &distance_index), 
                        zip.get_offset_in_chain(child_depth, &distance_index)), 
                        zip.get_length(child_depth, &distance_index));
            } else {
                //If it is a node that isn't reversed in the chain, or it's a snarl which is never reversed
                distance_end_left = std::numeric_limits<size_t>::max();
                distance_start_right = std::numeric_limits<size_t>::max();
                //Prefix sum of the child
                distance_start_left = zip.get_offset_in_chain(child_depth, &distance_index);
                //Length of the chain - prefix sum of the child - length of the child
                distance_end_right = SnarlDistanceIndex::minus(SnarlDistanceIndex::minus(
                        zip.get_length(child_depth-1, &distance_index), 
                        zip.get_offset_in_chain(child_depth, &distance_index)), 
                        zip.get_length(child_depth, &distance_index));
            }
#ifdef DEBUG_ZIPCODE
            cerr << "Distances to parent chain: " << distance_start_left << " " << distance_start_right << " " << distance_end_left << " " << distance_end_right << endl;
#endif
        }


        size_t new_distance_to_start = std::min(SnarlDistanceIndex::sum(distance_start_left, distance_to_start),
                                      SnarlDistanceIndex::sum(distance_end_left, distance_to_end));
        size_t new_distance_to_end = std::min(SnarlDistanceIndex::sum(distance_start_right, distance_to_start),
                                      SnarlDistanceIndex::sum(distance_end_right, distance_to_end));
        distance_to_start = new_distance_to_start;
        distance_to_end = new_distance_to_end;


    };

    if (zip1.get_distance_index_address(0) != zip2.get_distance_index_address(0)) {
#ifdef DEBUG_ZIPCODE
        cerr << "Zip codes are on different connected components" << endl;
#endif
        return std::numeric_limits<size_t>::max();
    }

    //The two positions are in the same connected component so now fill in the rest
    //of the decoder and try to find the distance
    zip1.fill_in_full_decoder();
    zip2.fill_in_full_decoder();

    //Now find the lowest common ancestor of the two zipcodes
    size_t lowest_common_ancestor_depth = 0;
    bool still_equal = true;
    while (still_equal) {

        if (lowest_common_ancestor_depth == zip1.decoder_length()-1 ||
            lowest_common_ancestor_depth == zip2.decoder_length()-1 ||
            !ZipCode::is_equal(zip1, zip2, lowest_common_ancestor_depth+1)) {
            //If we've hit the end of either decoder or if they are no longer equal,
            //Then break the loop and keep the current lowest_common_ancestor_depth
            still_equal = false;
        } else {
            //Otherwise increment lowest_common_ancestor_depth and keep going
            lowest_common_ancestor_depth ++;
        }
    }

#ifdef DEBUG_ZIPCODE
    vector<net_handle_t> ancestors;
    net_handle_t ancestor = distance_index.get_node_net_handle(id(pos1));
    while (!distance_index.is_root(ancestor)) {
        ancestors.push_back(ancestor);
        ancestor = distance_index.get_parent(ancestor);
    }
    ancestors.push_back(ancestor);
    cerr << "The lowest common ancestor is the " << lowest_common_ancestor_depth << "th thing from the root" << endl;
    cerr << "That should be " << distance_index.net_handle_as_string(ancestors[ancestors.size() - lowest_common_ancestor_depth - 1]) << endl; 
#endif


    if (distance_limit != std::numeric_limits<size_t>::max() &&
        lowest_common_ancestor_depth < zip1.decoder_length()-1){
        //If we're aborting when the distance is definitely too far,
        code_type_t ancestor_type = zip1.get_code_type(lowest_common_ancestor_depth);
        if  (ancestor_type == CHAIN || ancestor_type == ROOT_CHAIN) {
            //If the current ancestor is a chain, then check the distance
            size_t prefix_sum1 = zip1.get_offset_in_chain(lowest_common_ancestor_depth+1, &distance_index);
            size_t prefix_sum2 = zip2.get_offset_in_chain(lowest_common_ancestor_depth+1, &distance_index);
            size_t distance_in_chain; 
            if (prefix_sum1 < prefix_sum2) {
                //zip1 comes before zip2
                distance_in_chain = SnarlDistanceIndex::minus(
                    prefix_sum2, 
                    SnarlDistanceIndex::sum(prefix_sum1, 
                                            zip1.get_length(lowest_common_ancestor_depth+1, &distance_index)));
            } else {
                //zip2 comes before zip1
                distance_in_chain = SnarlDistanceIndex::minus(
                    prefix_sum1, 
                    SnarlDistanceIndex::sum(prefix_sum2, 
                                            zip2.get_length(lowest_common_ancestor_depth+1, &distance_index)));
            }
            if (distance_in_chain > distance_limit) {
                return std::numeric_limits<size_t>::max();
            }
        }
    }

    //Start from the nodes
    size_t distance_to_start1 = is_rev(pos1) 
        ? zip1.get_length(zip1.decoder_length()-1, &distance_index) - offset(pos1) 
        : offset(pos1) + 1;
    size_t distance_to_end1 = is_rev(pos1) ? offset(pos1) + 1 
         : zip1.get_length(zip1.decoder_length()-1, &distance_index) - offset(pos1);
    size_t distance_to_start2 = is_rev(pos2) 
         ? zip2.get_length(zip2.decoder_length()-1, &distance_index) - offset(pos2) 
         : offset(pos2) + 1;
    size_t distance_to_end2 = is_rev(pos2) ? offset(pos2) + 1 
         : zip2.get_length(zip2.decoder_length()-1, &distance_index) - offset(pos2);

    if (!undirected_distance) {
        //These are directed distances so set backwards distances to inf
        if (is_rev(pos1)) {
            distance_to_end1 = std::numeric_limits<size_t>::max();
        } else {
            distance_to_start1 = std::numeric_limits<size_t>::max();
        }
        if (is_rev(pos2)) {
            distance_to_start2 = std::numeric_limits<size_t>::max();
        } else {
            distance_to_end2 = std::numeric_limits<size_t>::max();
        }

    }
#ifdef DEBUG_ZIPCODE
cerr << "Distances in nodes: " << distance_to_start1 << " " << distance_to_end1 << " " << distance_to_start2 << " " << distance_to_end2 << endl;
cerr << "Finding distances to ancestors of first position" << endl;
#endif


    //Now walk up the snarl tree from each position to one level below the lowest common ancestor
    for (int i = zip1.decoder_length()-2 ; i > 0 && i > lowest_common_ancestor_depth ; i--) {
        //the parent snarl tree node is at index i
        //The distances are currently to the ends of the current node
        //Find the distances to the ends of the parent
        update_distances_to_ends_of_parent(zip1, i+1, distance_to_start1, distance_to_end1);
    }
#ifdef DEBUG_ZIPCODE
cerr << "Finding distances to ancestors of second position" << endl;
#endif
    //The same thing for the second position
    for (int i = zip2.decoder_length()-2 ; i > 0 && i > lowest_common_ancestor_depth ; i--) {
        //the parent snarl tree node is at index i
        //The distances are currently to the ends of the current node
        //Find the distances to the ends of the parent

        update_distances_to_ends_of_parent(zip2, i+1, distance_to_start2, distance_to_end2);
    }


    //Distances are now the distances to the ends of a child of the common ancestor

#ifdef DEBUG_ZIPCODE
    cerr << "Distances in children of common ancestor: " << distance_to_start1 << " " << distance_to_end1 << " " << distance_to_start2 << " " << distance_to_end2 << endl;
    //Check that the current nodes are actually children of the lca
    assert(ZipCode::is_equal(zip1, zip2, lowest_common_ancestor_depth));
#endif

    //Find the distance between them in the lowest common ancestor

    size_t distance_between = std::numeric_limits<size_t>::max();

    //Walk up the snarl tree from the lca and find the distance between the common ancestor
    for (int depth = lowest_common_ancestor_depth ; depth >= 0 ; depth--) {
        //Depth is the depth of a common ancestor. Current distances are to the ends of
        //a child of the common ancestor, at depth depth+1
#ifdef DEBUG_ZIPCODE
        cerr << "At " << depth << "st/th ancestor" << endl;
        cerr << "\tdistances are " << distance_to_start1 << " " << distance_to_end1 << " " << distance_to_start2 << " " << distance_to_end2 << endl;
#endif
        if (depth == zip1.decoder_length()-1) {
            //If the lca is a node that both positions are on

#ifdef DEBUG_ZIPCODE
            //If the lca is a node, then both the zipcode nodes should be the same node
            assert(ZipCode::is_equal(zip1, zip2, depth));
            assert(depth == zip2.decoder_length()-1);
            cerr << "\tAncestor should be a node" << endl;
#endif
            size_t d1 = SnarlDistanceIndex::sum(distance_to_end1, distance_to_start2);
            size_t d2 = SnarlDistanceIndex::sum(distance_to_end2, distance_to_start1);
            size_t node_length = zip1.get_length(depth, &distance_index);
            if (d1 > node_length) {
                distance_between = std::min(distance_between,
                                            SnarlDistanceIndex::minus(SnarlDistanceIndex::minus(d1, node_length),1));
            } 
            if (d2 > node_length) {
                distance_between = std::min(distance_between,
                                            SnarlDistanceIndex::minus(SnarlDistanceIndex::minus(d2, node_length),1));
            }
        } else if ( zip1.decoder[depth].is_chain) {
#ifdef DEBUG_ZIPCODE
            cerr << "\tancestor should be a chain" << endl;
#endif
            //If this ancestor is a chain

            //If the children are reversed in the chain, then flip their distances
            bool rev1 = (zip1.get_code_type(depth+1) == NODE && 
                zip1.get_is_reversed_in_parent(depth+1));
            size_t dist_start1 = rev1 ? distance_to_end1 : distance_to_start1;
            size_t dist_end1 = rev1 ? distance_to_start1 : distance_to_end1;

            bool rev2 = zip2.get_code_type(depth+1) == NODE && 
                zip2.get_is_reversed_in_parent(depth+1);
            size_t dist_start2 = rev2 ? distance_to_end2 : distance_to_start2;
            size_t dist_end2 = rev2 ? distance_to_start2 : distance_to_end2;

            //If they are the same child, then there is no path between them in the chain because we don't allow loops
            //So first check that they aren't the same
            if (!(ZipCode::is_equal(zip1, zip2, depth+1) 
                )){//TODO: I think this is unnecessary || (zip1.get_code_type(depth+1) == NODE && id(pos1) == id(pos2)))) 
                size_t prefix_sum1 = zip1.get_offset_in_chain(depth+1, &distance_index);
                size_t prefix_sum2 = zip2.get_offset_in_chain(depth+1, &distance_index);
                code_type_t code_type1 = zip1.get_code_type(depth+1);
                code_type_t code_type2 = zip2.get_code_type(depth+1);

                if (prefix_sum1 < prefix_sum2 ||
                    (prefix_sum1 == prefix_sum2 &&
                     (code_type1 == REGULAR_SNARL || code_type1 == IRREGULAR_SNARL || code_type1 == CYCLIC_SNARL)
                     && code_type2 == NODE)) {
                    //First child comes first in the chain
                    
                    if (code_type1 == REGULAR_SNARL || code_type1 == IRREGULAR_SNARL || code_type1 == CYCLIC_SNARL) {
                        //If the first thing is a snarl, then we need to take into account the length of the snarl
                        //(prefix sum 2 + distance left 2) - (prefix sum 1 + length 1) + distance right 1

#ifdef DEBUG_ZIPCODE
                        cerr << "First child comes first in the chain and it is a snarl" << endl;
                        cerr << "Find distances from : " << prefix_sum2 << " " << dist_start2 << " " << prefix_sum1 << " " << zip1.get_length(depth+1, &distance_index)  << " " << dist_end1 << endl;
#endif
                        if (dist_start2 != std::numeric_limits<size_t>::max()
                            && dist_end1 != std::numeric_limits<size_t>::max()) {
                            distance_between = std::min(distance_between,
                                                        SnarlDistanceIndex::minus(SnarlDistanceIndex::sum(
                                                            SnarlDistanceIndex::minus(
                                                                SnarlDistanceIndex::sum(prefix_sum2, 
                                                                                        dist_start2), 
                                                                SnarlDistanceIndex::sum(prefix_sum1,
                                                                                        zip1.get_length(depth+1, &distance_index))),
                                                             dist_end1),1));
                        }
                    } else {
                        //Otherwise, all that matters is the prefix sums
                        //(Prefix sum 2  + distance left 2) - (prefix sum1+ length 1) + distance right 1
#ifdef DEBUG_ZIPCODE
                        cerr << "First child comes first in the chain and it isn't a snarl" << endl;
                        cerr << "Find distances from : " << prefix_sum2 << " " << dist_start2 << " " << prefix_sum1 << " " << dist_end1 << " " << zip1.get_length(depth+1, &distance_index) << endl;
#endif
                        if (dist_start2 != std::numeric_limits<size_t>::max()
                            && dist_end1 != std::numeric_limits<size_t>::max()) {
                            distance_between = std::min(distance_between,
                                                        SnarlDistanceIndex::minus(
                                                        SnarlDistanceIndex::sum(
                                                        SnarlDistanceIndex::minus(
                                                            SnarlDistanceIndex::sum(prefix_sum2, 
                                                                                    dist_start2),
                                                            SnarlDistanceIndex::sum(prefix_sum1,
                                                                                    zip1.get_length(depth+1, &distance_index))), 

                                                            dist_end1),1) );
                        }
                    }
                } else {
                    //Second child comes first in the chain, or they are the same (doesn't matter)
                    if (code_type2 == REGULAR_SNARL || code_type2 == IRREGULAR_SNARL || code_type2 == CYCLIC_SNARL) {
                        //If the first thing is a snarl, then we need to take into account the length of the snarl
                        //(prefix sum 1 + distance left 1) - (prefix sum 2 + length 2) + distance right 2
#ifdef DEBUG_ZIPCODE
                        cerr << "Second child comes first in the chain and it is a snarl" << endl;
                        cerr << "Find distances from : " << prefix_sum1 << " " << dist_start1 << " " << prefix_sum2 << " " << zip2.get_length(depth+1, &distance_index)  << " " << dist_end2 << endl;
#endif
                        if (dist_start1 != std::numeric_limits<size_t>::max() 
                             && dist_end2 != std::numeric_limits<size_t>::max() ){
                            distance_between = std::min(distance_between,
                                                        SnarlDistanceIndex::minus(SnarlDistanceIndex::sum(
                                                            SnarlDistanceIndex::minus(
                                                                SnarlDistanceIndex::sum(prefix_sum1, 
                                                                                        dist_start1), 
                                                                SnarlDistanceIndex::sum(prefix_sum2,
                                                                                        zip2.get_length(depth+1, &distance_index))),
                                                             dist_end2), 1));
                        }
                    } else {
                        //Otherwise, all that matters is the prefix sums
                        //(Prefix sum 1  + distance left 1) - (prefix sum2 + length 2) + distance right 2
#ifdef DEBUG_ZIPCODE
                        cerr << "Second child comes first in the chain and it isn't a snarl" << endl;
                        cerr << "Find distances from : " << prefix_sum1 << " " << dist_start1 << " " << prefix_sum2 << " " << dist_end2 << endl;
#endif
                        if (dist_start1 != std::numeric_limits<size_t>::max() 
                             && dist_end2 != std::numeric_limits<size_t>::max() ){
                            distance_between = std::min(distance_between,
                                                        SnarlDistanceIndex::minus(
                                                        SnarlDistanceIndex::sum(
                                                        SnarlDistanceIndex::minus(
                                                            SnarlDistanceIndex::sum(prefix_sum1, 
                                                                                    dist_start1),
                                                            SnarlDistanceIndex::sum(prefix_sum2,
                                                                                    zip2.get_length(depth+1, &distance_index))), 

                                                            dist_end2),1) );
                        }
                    }
                }
            }
            //Update distances from the ends of the children (at depth+1) to parent (depth)
            update_distances_to_ends_of_parent(zip1, depth+1, distance_to_start1, distance_to_end1);
            update_distances_to_ends_of_parent(zip2, depth+1, distance_to_start2, distance_to_end2);
        } else {

#ifdef DEBUG_ZIPCODE
            cerr << "\tancestor is a snarl" << endl;
#endif
            //If the ancestor is a snarl
            
            //If the parent is a regular snarl, then there is no path between them so
            //just update the distances to the ends of the parent 
            if (zip1.get_code_type(depth) != REGULAR_SNARL) {
                //Parent may be an irregular snarl or a root snarl (which is also irregular)
                net_handle_t parent_handle = zip1.get_net_handle(depth, &distance_index);
                size_t rank1 = zip1.get_rank_in_snarl(depth+1);
                size_t rank2 = zip2.get_rank_in_snarl(depth+1);
#ifdef DEBUG_ZIPCODE
                cerr << "irregular snarl so find distances in the distance index: " << distance_index.net_handle_as_string(parent_handle) << endl;
                cerr << "\t at offset " << distance_index.get_record_offset(parent_handle) << endl;
                cerr << "ranks: " << rank1 << " and " << rank2 << endl;
#endif

                size_t distance_start_start = distance_index.distance_in_snarl(parent_handle, 
                                    rank1, false, rank2, false, graph);
                size_t distance_start_end = distance_index.distance_in_snarl(parent_handle, 
                                    rank1, false, rank2, true, graph);
                size_t distance_end_start = distance_index.distance_in_snarl(parent_handle, 
                                    rank1, true, rank2, false, graph);
                size_t distance_end_end = distance_index.distance_in_snarl(parent_handle, 
                                    rank1, true, rank2, true, graph);
                size_t distance_between_snarl = std::min( SnarlDistanceIndex::sum(SnarlDistanceIndex::sum(
                                                distance_to_start1, distance_to_start2), distance_start_start),
                                   std::min( SnarlDistanceIndex::sum(SnarlDistanceIndex::sum(
                                                distance_to_start1, distance_to_end2), distance_start_end),
                                   std::min( SnarlDistanceIndex::sum(SnarlDistanceIndex::sum(
                                                distance_to_end1, distance_to_start2), distance_end_start),
                                             SnarlDistanceIndex::sum(SnarlDistanceIndex::sum(
                                                distance_to_end1, distance_to_end2), distance_end_end))));

                distance_between = std::min(distance_between,
                                            SnarlDistanceIndex::minus(distance_between_snarl, 1));
            }
#ifdef DEBUG_ZIPCODE
            else {
                cerr << "\tAncestor is a regular snarl so there is no path between the children" << endl;
            }
#endif
            //Update distances from the ends of the children (at depth+1) to parent (depth)
            update_distances_to_ends_of_parent(zip1, depth+1, distance_to_start1, distance_to_end1);
            update_distances_to_ends_of_parent(zip2, depth+1, distance_to_start2, distance_to_end2);
        }
#ifdef DEBUG_ZIPCODE
        cerr << "distance in ancestor: " << distance_between << endl;
#endif
    }

    return distance_between;
}

bool ZipCode::is_farther_than(const ZipCode& zip1, const ZipCode& zip2, const size_t& limit){
#ifdef DEBUG_ZIPCODE
    cerr << "Checking if two zip codes are farther than " << limit << endl;
#endif

    size_t zip_index1 = 0; size_t zip_index2 = 0;
    size_t zip_value1 = std::numeric_limits<size_t>::max();
    size_t zip_value2 = std::numeric_limits<size_t>::max();

    //If the two positions aren't on the same connected component, then we're done
    for (size_t i = 0 ; i <= ROOT_IS_CHAIN_OFFSET ; i++) { 
        std::tie(zip_value1, zip_index1) = zip1.zipcode.get_value_and_next_index(zip_index1);
        std::tie(zip_value2, zip_index2) = zip2.zipcode.get_value_and_next_index(zip_index2);
    }
    if (zip_value1 != zip_value2) {
#ifdef DEBUG_ZIPCODE
        cerr << "Zip codes are on different connected components" << endl;
#endif
        return true;
    }

    bool is_top_level_chain = zip_value1;
    for (size_t i = 0 ; i <= ROOT_IDENTIFIER_OFFSET - ROOT_IS_CHAIN_OFFSET - 1; i++) { 
        std::tie(zip_value1, zip_index1) = zip1.zipcode.get_value_and_next_index(zip_index1);
        std::tie(zip_value2, zip_index2) = zip2.zipcode.get_value_and_next_index(zip_index2);
    }
    if (zip_value1 != zip_value2) {
#ifdef DEBUG_ZIPCODE
        cerr << "Zip codes are on different connected components" << endl;
#endif
        return true;
    }

    if (!is_top_level_chain) {
        //If the top-level thing is a snarl, then check if the zips are in the same chain. 
        //If they are, then proceed from the shared chain

        //The next thing will be the identifier for the chain
        for (size_t i = 0 ; i <= CHAIN_RANK_IN_SNARL_OFFSET; i++) { 
            std::tie(zip_value1, zip_index1) = zip1.zipcode.get_value_and_next_index(zip_index1);
            std::tie(zip_value2, zip_index2) = zip2.zipcode.get_value_and_next_index(zip_index2);
        }
        if (zip_value1 != zip_value2) {
            //We can't tell
            return false;
        }
        //Next is the length of the chain
        for (size_t i = 0 ; i <= CHAIN_LENGTH_OFFSET - CHAIN_RANK_IN_SNARL_OFFSET - 1; i++) { 
            std::tie(zip_value1, zip_index1) = zip1.zipcode.get_value_and_next_index(zip_index1);
            std::tie(zip_value2, zip_index2) = zip2.zipcode.get_value_and_next_index(zip_index2);
        }
        if (zip_value1 < limit) {
            return true;
        }

        //The zips now point to the children of the shared chain, so we can proceed as if the top-level
        //structure was a chain

    } else {
        //If it is a chain, get two more things to get to the end of the chain
        for (size_t i = 0 ; i < 2 ; ++i) {
            std::tie(zip_value1, zip_index1) = zip1.zipcode.get_value_and_next_index(zip_index1);
            std::tie(zip_value2, zip_index2) = zip2.zipcode.get_value_and_next_index(zip_index2);
        }
    }

    //Both zips now point to a thing in a shared chain
    //Get the minimum possible distance between the structures on the chain
    //For a lower bound, this assumes that the positions are as close as they can be on the structure in the chain
    size_t prefix_sum1, prefix_sum2, length1, length2, component1, component2;

    //The next thing could either be a snarl or a node. If it is a node, 
    vector<size_t> next_values;
    for (size_t i = 0 ; i < NODE_SIZE ; i++ ) {
#ifdef DEBUG_ZIPCODE
        assert(zip_index1 != std::numeric_limits<size_t>::max());
#endif
        std::tie(zip_value1, zip_index1) = zip1.zipcode.get_value_and_next_index(zip_index1);
        next_values.emplace_back(zip_value1);
    }
    if (zip_index1 == std::numeric_limits<size_t>::max()) {
#ifdef DEBUG_ZIPCODE
        cerr << "zip1 is a node in a chain" << endl;
#endif
        //If the last thing was a node
        prefix_sum1 = next_values[0];
        length1 = next_values[1];
        component1 = next_values[2];
        prefix_sum1 = prefix_sum1 == 0 ? std::numeric_limits<size_t>::max() : prefix_sum1-1;
        length1 = length1 == 0 ? std::numeric_limits<size_t>::max() : length1-1;
    } else {
#ifdef DEBUG_ZIPCODE
        cerr << "zip1 is in a snarl in a chain" << endl;
#endif
        //If the last thing was a snarl
        if (next_values[0]) {
            //If the next thing was a regular snarl
            prefix_sum1 = next_values[1];
            length1 = next_values[2];
            std::tie(zip_value1, zip_index1) = zip1.zipcode.get_value_and_next_index(zip_index1);
            component1 = zip_value1;
            prefix_sum1 = prefix_sum1 == 0 ? std::numeric_limits<size_t>::max() : prefix_sum1-1;
            length1 = length1 == 0 ? std::numeric_limits<size_t>::max() : length1-1;
        } else {
            //If the next thing was an irregular snarl
            //TODO: If it's an irregular snarl, then we don't actually store the relevant values so we can't tell. Could look it up in the distance index or store it
            return false;
        }
    }
    
    //Do the same for the other zip
    next_values.clear();
    for (size_t i = 0 ; i < NODE_SIZE ; i++ ) {
#ifdef DEBUG_ZIPCODE
        assert(zip_index2 != std::numeric_limits<size_t>::max());
#endif
        std::tie(zip_value2, zip_index2) = zip2.zipcode.get_value_and_next_index(zip_index2);
        next_values.emplace_back(zip_value2);
    }
    if (zip_index2 == std::numeric_limits<size_t>::max()) {
#ifdef DEBUG_ZIPCODE
        cerr << "zip2 is a node in a chain" << endl;
#endif
        //If the last thing was a node
        prefix_sum2 = next_values[0];
        length2 = next_values[1];
        component2 = next_values[2];
        prefix_sum2 = prefix_sum2 == 0 ? std::numeric_limits<size_t>::max() : prefix_sum2-1;
        length2 = length2 == 0 ? std::numeric_limits<size_t>::max() : length2-1;
    } else {
#ifdef DEBUG_ZIPCODE
        cerr << "zip2 is in a snarl in a chain" << endl;
#endif
        //If the last thing was a snarl
        if (next_values[0]) {
            //If the next thing was a regular snarl
            prefix_sum2 = next_values[1];
            length2 = next_values[2];
            std::tie(zip_value2, zip_index2) = zip2.zipcode.get_value_and_next_index(zip_index2);
            component2 = zip_value2;
            prefix_sum2 = prefix_sum2 == 0 ? std::numeric_limits<size_t>::max() : prefix_sum2-1;
            length2 = length2 == 0 ? std::numeric_limits<size_t>::max() : length2-1;
        } else {
            //If the next thing was an irregular snarl
            //TODO: If it's an irregular snarl, then we don't actually store the relevant values so we can't tell. Could look it up in the distance index or store it
            return false;
        }
    }
#ifdef DEBUG_ZIPCODE
    cerr << "Finding distance in chain between " << prefix_sum1 << " " << length1 << " and " << prefix_sum2 << " and " << length2 << endl;
#endif

    if (component1 != component2 ||
        prefix_sum1 == std::numeric_limits<size_t>::max() ||
        prefix_sum2 == std::numeric_limits<size_t>::max() ||
        length1 == std::numeric_limits<size_t>::max() ||
        length2 == std::numeric_limits<size_t>::max()) {
        //If anything is infinite, then we can't tell
        return false;
    }


    if (prefix_sum1 < prefix_sum2) {
        //If 1 comes first

        if (prefix_sum1 + length1 > prefix_sum2) {
            //They might be close
            return false;
        } else {
            //Return true if the distance between is greater than the limit
            return prefix_sum2 - (prefix_sum1 + length1) > limit; 
        }
    } else {
        //If 2 comes first

        if (prefix_sum2 + length2 > prefix_sum1) {
            //They might be close
            return false;
        } else {
            //Return true if the distance between is greater than the limit
            return prefix_sum1 - (prefix_sum2 + length2) > limit; 
        }
    }
}

gbwtgraph::Payload ZipCode::get_payload_from_zip() const {
    varint_vector_t decoder_vector;
    //The zipcode decoder's is_chain will always alternate is_chain between levels, except for the very end,
    // which may have two is_chains in a row for a trivial chain. So we can store the whole series in two bits.
    //For the decoder, we never need to know the byte count, since the value in the decoder is never 0


    //TODO: This is assuming the decoder is filled in already
    bool is_root_chain = decoder[0].is_chain;
    bool is_trivial_chain = decoder.size() > 1 && decoder[decoder.size()-1].is_chain && decoder[decoder.size()-2].is_chain;
    size_t is_chain_value = 0;
    if (is_root_chain) {
        is_chain_value |= 1;
    }
    if (is_trivial_chain) {
        is_chain_value |= 1<<1;
    }
    decoder_vector.add_value(is_chain_value);
    //The first offset is always 0 so ignore it
    for (const ZipCode::decoder_t& d : decoder) {
        if (d.offset != 0) {
            decoder_vector.add_value(d.offset);
        }
    }

    //First byte is for the byte_count
    if (byte_count() + decoder_vector.byte_count() > 15) {
        //If there aren't enough bits to represent the zip code
        return MIPayload::NO_CODE;
    }

    //Encode it as the byte count of the zipcode, the zipcode, and the decoder
    
    //Index and value as we walk through the zip code
    size_t index = 0;
    size_t value;
    
    //The values that get returned
    code_type encoded1 = 0;
    code_type encoded2 = 0;

    encoded1 |= byte_count();

    size_t encoded_bytes = 1;

    for (size_t i = 0 ; i < zipcode.data.size() ; i++ ) {
       size_t byte = static_cast<size_t> (zipcode.data[i]); 
       if ( encoded_bytes < 8 ) {
            //Add to first code
            encoded1 |= (byte << (encoded_bytes*8));

        } else {
            //Add to second code
            encoded2 |= (byte << ((encoded_bytes-8)*8));
        }
        encoded_bytes++;
    
    }
    for (size_t i = 0 ; i < decoder_vector.data.size() ; i++) {
       size_t byte = static_cast<size_t> (decoder_vector.data[i]); 
       if ( encoded_bytes < 8 ) {
            //Add to first code
            encoded1 |= (byte << (encoded_bytes*8));

        } else {
            //Add to second code
            encoded2 |= (byte << ((encoded_bytes-8)*8));
        }
        encoded_bytes++;
    }
    assert(encoded_bytes <= 16);
    return {encoded1, encoded2};

}

void ZipCode::fill_in_zipcode_from_payload(const gbwtgraph::Payload& payload) {
    assert(payload != MIPayload::NO_CODE);
    zipcode.data.reserve(16);

    size_t decoded_bytes = 0;

    //get one byte at a time from the payload and add it to the zip code
    size_t bit_mask = (1 << 8) - 1;
    size_t byte_count = payload.first & bit_mask;
    decoded_bytes++;
    for (size_t i = 0 ; i < byte_count ; i++) {
        if (decoded_bytes < 8) {
            zipcode.add_one_byte((payload.first >> (decoded_bytes*8)) & bit_mask);
        } else {
            zipcode.add_one_byte((payload.second >> ((decoded_bytes-8)*8)) & bit_mask);
        }
        decoded_bytes++;
    }

    //Find the booleans specifying the is_chain values
    uint8_t is_chain_val = 0;
    if (decoded_bytes < 8) {
        is_chain_val = (payload.first >> (decoded_bytes*8)) & bit_mask;
    } else {
        is_chain_val = (payload.second >> ((decoded_bytes-8)*8)) & bit_mask;
    }
    decoded_bytes++;
    bool is_chain = is_chain_val & 1;
    bool is_trivial_chain = is_chain_val & (1<<1);

    //Get the decoder offsets
    varint_vector_t decoder_vector;
    decoder_vector.data.reserve(16-decoded_bytes);
    for (size_t i = decoded_bytes ; i <16 ; i++) {
        uint8_t saved_byte;
        if (decoded_bytes < 8) {
            saved_byte = (payload.first >> (decoded_bytes*8)) & bit_mask;
        } else {
            saved_byte = (payload.second >> ((decoded_bytes-8)*8)) & bit_mask;
        }
        if (saved_byte != 0) {
            decoder_vector.add_one_byte(saved_byte);
        }
        
        decoded_bytes++;
    }
    //Now go through the varint vector up and add anything that isn't 0
    size_t varint_value= 1;
    size_t varint_index = 0;
    decoder.reserve(decoder_vector.byte_count());
    decoder.emplace_back(is_chain, 0);
    is_chain = !is_chain;
    if (decoder_vector.byte_count() != 0) {
        while (varint_index != std::numeric_limits<size_t>::max() && varint_value != 0) {
            std::tie(varint_value, varint_index) = decoder_vector.get_value_and_next_index(varint_index);
            
            decoder.emplace_back(is_chain, varint_value);

            is_chain = !is_chain;
        }
    }
    if (is_trivial_chain) {
        assert(!decoder.back().is_chain);
        decoder.back().is_chain = true;
    }

}

std::ostream& operator<<(std::ostream& out, const ZipCode::code_type_t& type) {
    if (type == ZipCode::NODE) {
        return out << "NODE";
    } else if (type == ZipCode::CHAIN) {
        return out << "CHAIN";
    } else if (type == ZipCode::REGULAR_SNARL) {
        return out << "REGULAR_SNARL";
    } else if (type == ZipCode::IRREGULAR_SNARL) {
        return out << "IRREGULAR_SNARL";
    } else if (type == ZipCode::CYCLIC_SNARL) {
        return out << "CYCLIC_SNARL";
    } else if (type == ZipCode::ROOT_SNARL) {
        return out << "ROOT_SNARL";
    } else if (type == ZipCode::ROOT_CHAIN) {
        return out << "ROOT_CHAIN";
    } else if (type == ZipCode::ROOT_NODE) {
        return out << "ROOT_NODE";
    } else if (type == ZipCode::EMPTY) {
        return out << "EMPTY";
    } else {
        throw std::runtime_error("error: Trying to print an invalid code_type_t");
    }
}


void ZipCodeCollection::serialize(std::ostream& out) const {
    //The zipcode vector will be serialized as a bunch of varint_vector_ts
    //The first varint_vector_t will have one value, which will be the length of the
    //zipcode that follows it

    //First serialize the header, which is the magic number and version
    uint32_t magic = magic_number;
    uint32_t vers = version;
    out.write(reinterpret_cast<const char*>(&magic), sizeof(magic));
    out.write(reinterpret_cast<const char*>(&vers), sizeof(vers));


    for (const ZipCode& zip : zipcodes) {
    
        //How many bytes are going to be saved for the zipcode? 
        size_t byte_count = zip.byte_count();

        varint_vector_t size_vector;
        size_vector.add_value(byte_count);
        //Write the number of bytes about to be saved 
        for (const uint8_t& byte : size_vector.data) {
            out << char(byte);
        } 

        //Write the zipcode 
#ifdef DEBUG_ZIPCODE
        size_t zip_byte_count = 0;
#endif
        for (const uint8_t& byte : zip.zipcode.data ) {
#ifdef DEBUG_ZIPCODE
            zip_byte_count++;
#endif
            out << char(byte);
        }
#ifdef DEBUG_ZIPCODE
        assert(byte_count == zip_byte_count); 
#endif

        //Also save the decoder
        varint_vector_t decoder_vector;
        for (const ZipCode::decoder_t& d : zip.decoder) {
            decoder_vector.add_value(d.is_chain);
            decoder_vector.add_value(d.offset);
        }

        //Write the number of bytes for the zipcode
        varint_vector_t decoder_byte_count;
        decoder_byte_count.add_value(decoder_vector.byte_count());
        for (const uint8_t& byte : decoder_byte_count.data) {
            out << char(byte);
        }


        //Write the decoder
        for (const uint8_t& byte : decoder_vector.data ) {
            out << char(byte);
        }
    }

}
void ZipCodeCollection::deserialize(std::istream& in) {

    //Check the magic number and version
    uint32_t saved_magic_number, saved_version;
    in.read(reinterpret_cast<char*>(&saved_magic_number), sizeof(saved_magic_number));
    if (saved_magic_number != magic_number) {
        throw std::runtime_error("error: Loading the wrong type of file when looking for zipcodes");
    }

    in.read(reinterpret_cast<char*>(&saved_version), sizeof(saved_version));
    if (saved_version != version) {
        throw std::runtime_error("error: Loading the wrong zipcode version");
    }

    while (in.peek() != EOF) {

        //First, get the number of bytes used by the zipcode
        //This will be a varint_vector_t with one value, which is the number of bytes in the zipcode
        //Each byte in the varint_vector_t starts with 0 if it is the last bit in the 
        //number, and 1 if the next byte is included
        varint_vector_t byte_count_vector;
        while (in.peek() & (1<<7)) {
            //If the first bit in the byte is 1, then add it, stop once the first bit is 0 
            char c;
            in.get(c);
            byte_count_vector.add_one_byte((uint8_t)c);
        }
        assert(! (in.peek() & (1<<7))); 
        //The next byte has a 0 as its first bit, so add it
        char c;
        in.get(c);
        byte_count_vector.add_one_byte((uint8_t)c);

        //The first (and only) value in the vector is the length of the zipcode
        size_t zipcode_byte_count = byte_count_vector.get_value_and_next_index(0).first;

#ifdef DEBUG_ZIPCODE
        cerr << "Get zipcode of " << zipcode_byte_count << " bytes" << endl;
        //assert(zipcode_byte_count >= 15);
        assert(byte_count_vector.get_value_and_next_index(0).second == std::numeric_limits<size_t>::max());
#endif

        char line [zipcode_byte_count];

        in.read(line, zipcode_byte_count);

        ZipCode zip;
        for (const char& character : line) {
            zip.zipcode.add_one_byte(uint8_t(character));
        }


        //Now get the decoder

        varint_vector_t decoder_byte_count_vector;
        while (in.peek() & (1<<7)) {
            //If the first bit in the byte is 1, then add it, stop once the first bit is 0 
            char ch;
            in.get(ch);
            decoder_byte_count_vector.add_one_byte((uint8_t)ch);
        }
        assert(! (in.peek() & (1<<7))); 
        //The next byte has a 0 as its first bit, so add it
        char ch;
        in.get(ch);
        decoder_byte_count_vector.add_one_byte((uint8_t)ch);

        //The first (and only) value in the vector is the length of the zipcode
        size_t decoder_byte_count = decoder_byte_count_vector.get_value_and_next_index(0).first;

#ifdef DEBUG_ZIPCODE
        cerr << "Get decoder of " << decoder_byte_count << " bytes" << endl;
        //assert(decoder_byte_count >= 15);
        assert(decoder_byte_count_vector.get_value_and_next_index(0).second == std::numeric_limits<size_t>::max());
#endif

        char line1 [decoder_byte_count];

        in.read(line1, decoder_byte_count);

        varint_vector_t decoder_vector;
        for (const char& character : line1) {
            decoder_vector.add_one_byte(uint8_t(character));
        }

        if (decoder_vector.byte_count() != 0) {
            size_t index = 0;
            while (index != std::numeric_limits<size_t>::max()) {
                size_t is_chain, offset;
                std::tie(is_chain, index) = decoder_vector.get_value_and_next_index(index);
                std::tie(offset, index) = decoder_vector.get_value_and_next_index(index);
                zip.decoder.emplace_back(is_chain != 0, offset);
            }
        }
        zip.finished_decoding=true;


        zipcodes.emplace_back(std::move(zip));
    }

}
MIPayload ZipCode::get_payload_from_zipcode(nid_t id, const SnarlDistanceIndex& distance_index, hash_map<size_t, net_handle_t>* component_to_net_handle) const {
    MIPayload payload;

    if (decoder_length() == 1) {
        //If the root-level structure is a node
        payload.parent_is_root = true;
        payload.parent_is_chain = true;

        //Walk through the zipcode to get values
        size_t zip_value;
        size_t zip_index = decoder[0].offset;
        //Root is chain
        std::tie(zip_value, zip_index) = zipcode.get_value_and_next_index(zip_index);
        //root_identifier
        std::tie(zip_value, zip_index) = zipcode.get_value_and_next_index(zip_index);
        if (component_to_net_handle!= nullptr && component_to_net_handle->count(zip_value)) {
            payload.node_handle = component_to_net_handle->at(zip_value);
        } else {
            payload.node_handle = distance_index.get_net_handle_from_values(distance_index.get_record_offset(distance_index.get_handle_from_connected_component(zip_value)),
                                                                          SnarlDistanceIndex::START_END,
                                                                          SnarlDistanceIndex::CHAIN_HANDLE);
            if (component_to_net_handle!= nullptr) {
                component_to_net_handle->emplace(zip_value, payload.node_handle);
            }
        }

        //Root node length
        std::tie(zip_value, zip_index) = zipcode.get_value_and_next_index(zip_index);

        payload.node_length =  zip_value == std::numeric_limits<size_t>::max() ? 0 : zip_value-1;
        payload.is_trivial_chain = true;
        payload.is_reversed = false;
        payload.parent_handle = distance_index.get_root();
        payload.parent_type = ZipCode::ROOT_NODE;
        payload.parent_record_offset = 0;

    } else if (decoder[max_depth() - 1].is_chain) {
        //If the parent is a chain
        payload.node_handle = distance_index.get_node_net_handle(id);
        payload.parent_is_chain = true;
        payload.parent_is_root = false;

        //Walk through the zipcode to get values
        size_t zip_value;
        size_t zip_index = decoder[max_depth()-1].offset;
        //is_chain/rank in snarl
        std::tie(zip_value, zip_index) = zipcode.get_value_and_next_index(zip_index);

        //root_identifier for root, chain length for anything else
        std::tie(zip_value, zip_index) = zipcode.get_value_and_next_index(zip_index);

        if (decoder_length() == 2) {
            //If the node is a child of the root chain
            if (component_to_net_handle!= nullptr && component_to_net_handle->count(zip_value)) {
                payload.parent_handle = component_to_net_handle->at(zip_value);
            } else {
                payload.parent_handle = distance_index.start_end_traversal_of(distance_index.get_handle_from_connected_component(zip_value));
                if (component_to_net_handle!= nullptr) {
                    component_to_net_handle->emplace(zip_value, payload.parent_handle);
                }
            }
            payload.parent_type = ZipCode::ROOT_CHAIN;
            payload.parent_is_root = true;
            std::tie(zip_value, zip_index) = zipcode.get_value_and_next_index(zip_index);
        } else {
            payload.parent_handle = distance_index.start_end_traversal_of(distance_index.get_parent(payload.node_handle));
            payload.parent_type = ZipCode::CHAIN;
        }
        payload.parent_record_offset = distance_index.get_record_offset(payload.parent_handle);

        //chain component count
        std::tie(zip_value, zip_index) = zipcode.get_value_and_next_index(zip_index);

        //Node prefix sum
        std::tie(zip_value, zip_index) = zipcode.get_value_and_next_index(zip_index);
        payload.prefix_sum = zip_value == std::numeric_limits<size_t>::max() ? 0 : zip_value-1;
        //Node length
        std::tie(zip_value, zip_index) = zipcode.get_value_and_next_index(zip_index);
        payload.node_length = zip_value == std::numeric_limits<size_t>::max() ? 0 : zip_value-1;
        //is_reversed
        std::tie(zip_value, zip_index) = zipcode.get_value_and_next_index(zip_index);
        //TODO: For top-level chains we got this from the distance index
        payload.is_reversed = zip_value;

        std::tie(zip_value, zip_index) = zipcode.get_value_and_next_index(zip_index);
        payload.chain_component = zip_value;



    } else {
        //If the node is a child of a snarl
        
        payload.node_handle = distance_index.get_node_net_handle(id);
        payload.parent_handle = distance_index.get_net_handle_from_values(distance_index.get_record_offset(payload.node_handle),
                                                         SnarlDistanceIndex::START_END,
                                                         SnarlDistanceIndex::CHAIN_HANDLE,
                                                         distance_index.get_node_record_offset(payload.node_handle));
        payload.parent_is_chain = false;
        payload.parent_is_root = decoder_length() == 2;
        payload.is_trivial_chain = true;


        size_t zip_value;
        size_t zip_index;
        if (payload.parent_is_root) {
            //is_chain
            zip_index = decoder[0].offset;
            std::tie(zip_value, zip_index) = zipcode.get_value_and_next_index(zip_index);
            //Identifier for root snarl
            std::tie(zip_value, zip_index) = zipcode.get_value_and_next_index(zip_index);
            payload.node_handle = payload.parent_handle;
            if (component_to_net_handle!= nullptr && component_to_net_handle->count(zip_value)) {
                payload.parent_handle = component_to_net_handle->at(zip_value);
                payload.parent_record_offset = distance_index.get_record_offset(payload.parent_handle);
            } else {
                    payload.parent_record_offset = distance_index.get_record_offset(distance_index.get_handle_from_connected_component(zip_value));
                    payload.parent_handle = distance_index.get_net_handle_from_values(payload.parent_record_offset,
                                                                                      SnarlDistanceIndex::START_END,
                                                                                      SnarlDistanceIndex::ROOT_HANDLE);
                if (component_to_net_handle!= nullptr) {
                    component_to_net_handle->emplace(zip_value, payload.parent_handle);
                }
            }
            payload.parent_type = ZipCode::ROOT_SNARL;
        } else {
            zip_index = decoder[max_depth()-1].offset;
            //is_regular
            std::tie(zip_value, zip_index) = zipcode.get_value_and_next_index(zip_index);
            //If this is a non-root snarl, get as much as we can from it 
            payload.parent_type = ZipCode::EMPTY;
            if (zip_value == 0) {
                payload.parent_type = ZipCode::IRREGULAR_SNARL;
            } else if (zip_value == 1) {
                payload.parent_type = ZipCode::REGULAR_SNARL;
            } else {
                payload.parent_type = ZipCode::CYCLIC_SNARL;
            }

            //Snarl prefix sum
            std::tie(zip_value, zip_index) = zipcode.get_value_and_next_index(zip_index);

            payload.prefix_sum = 0; //TODO: SHould use this zip_value == std::numeric_limits<size_t>::max() ? 0 : zip_value-1;

            //Snarl length
            std::tie(zip_value, zip_index) = zipcode.get_value_and_next_index(zip_index);
            //Snarl child_count
            std::tie(zip_value, zip_index) = zipcode.get_value_and_next_index(zip_index);
            //Chain component of the snarl
            std::tie(zip_value, zip_index) = zipcode.get_value_and_next_index(zip_index);
            //TODO: SHould use this somehow
            payload.chain_component = 0;
            //is_reversed for regular snarl and record offset for irregular/cyclic snarl
            std::tie(zip_value, zip_index) = zipcode.get_value_and_next_index(zip_index);

            if (payload.parent_type == ZipCode::REGULAR_SNARL) {
                //Snarl is reversed
                net_handle_t grandparent_handle = distance_index.get_parent(payload.parent_handle);
                //Simple and regular snarls are different for clustering
                if (distance_index.is_simple_snarl(grandparent_handle)) {
                    payload.is_reversed = zip_value;
                    payload.parent_is_chain=true;
                    payload.parent_record_offset = distance_index.get_record_offset(distance_index.get_parent(grandparent_handle));
                } else {
                    payload.is_reversed = false;
                    payload.parent_record_offset = distance_index.get_record_offset(grandparent_handle);
                }

            } else {
                payload.is_reversed = false;
                payload.parent_record_offset = zip_value;
            }

        }
        //We should be at the node/trivial chain now
        zip_index = decoder[max_depth()].offset;
        //Chain rank in snarl
        std::tie(zip_value, zip_index) = zipcode.get_value_and_next_index(zip_index);
        //Chain length
        std::tie(zip_value, zip_index) = zipcode.get_value_and_next_index(zip_index);
        payload.node_length = zip_value == std::numeric_limits<size_t>::max() ? 0 : zip_value-1;

        //Get the rest as default values

    }
    payload.parent_depth = 0;
    for (size_t d = 0 ; d <= max_depth() ; d++) {
        auto type = get_code_type(d);
        if (type == ZipCode::CHAIN || type == ZipCode::ROOT_CHAIN || type == ZipCode::ROOT_NODE) {
            payload.parent_depth++;
        }
    }



    return payload;
}

net_identifier_t ZipCode::get_identifier(size_t depth) const {
    if (depth == std::numeric_limits<size_t>::max()) {
        //This is equivalent to distance_index.get_root()
        return "ROOT";
    }
    string result = "";
    for (size_t d = 0 ; d < depth ; d++) {
        result += (decoder[d].is_chain ? "1" : "0");
        if (d == 0) {
            //Root structure
            size_t zip_value;
            size_t zip_index = decoder[d].offset;
            for (size_t i = 0 ; i <= ZipCode::ROOT_IDENTIFIER_OFFSET; i++) {
                std::tie(zip_value, zip_index) = zipcode.get_value_and_next_index(zip_index);
                result += std::to_string(zip_value);
            }
        } else if (decoder[d].is_chain) {
            //is_chain so could be a chain or a node
            if (decoder[d-1].is_chain) {
                //If the thing before this was also a chain, then it is a node
                size_t zip_value;
                size_t zip_index = decoder[d].offset;
                for (size_t i = 0 ; i <= ZipCode::NODE_OFFSET_OFFSET; i++) {
                    std::tie(zip_value, zip_index) = zipcode.get_value_and_next_index(zip_index);
                    result += std::to_string(zip_value);
                }
            } else {
                //Otherwise it's a chain
                size_t zip_value;
                size_t zip_index = decoder[d].offset;
                for (size_t i = 0 ; i <= ZipCode::CHAIN_RANK_IN_SNARL_OFFSET; i++) {
                    std::tie(zip_value, zip_index) = zipcode.get_value_and_next_index(zip_index);
                    result += std::to_string(zip_value);
                }
            }
        } else {
            //Definitely a snarl
            size_t zip_value;
            size_t zip_index = decoder[d].offset;
            for (size_t i = 0 ; i <= ZipCode::SNARL_OFFSET_IN_CHAIN_OFFSET; i++) {
                std::tie(zip_value, zip_index) = zipcode.get_value_and_next_index(zip_index);
                result += std::to_string(zip_value);
            }
        }
        if (d < std::min(depth, max_depth())) {
            result += ".";
        }
        
    }
    if (depth > max_depth()) {
        //If this was node that's in a trivial chain
        result += ".n";
    }

    return result;
}

const net_identifier_t ZipCode::get_parent_identifier(const net_identifier_t& child) {
    if (child == "ROOT") {
        throw std::runtime_error("error: trying to get the parent of the root net_identifier_t. Do your zipcode, minimizer, and graph files match?");
    }
    for (int i = child.size()-1 ; i >= 0 ; i--) {
        if (child[i] == '.') {
            return (net_identifier_t) string(child, 0, i);
        }
    }
    //If we didn't find a '.', then the parent is just the root
    return "ROOT";
}



}
