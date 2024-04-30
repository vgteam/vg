#include "zip_code.hpp"

//#define DEBUG_ZIPCODE

namespace vg{
using namespace std;

void ZipCode::fill_in_zipcode (const SnarlDistanceIndex& distance_index, const pos_t& pos) {

    std::vector<net_handle_t> ancestors;
    net_handle_t current_handle = distance_index.get_node_net_handle(id(pos));

    //Put all ancestors of the node in a vector, starting from the node, and not including the root
    while (!distance_index.is_root(current_handle)) {
        ancestors.emplace_back(current_handle);
        current_handle = distance_index.get_parent(current_handle);
    }


    //Now add the root-level snarl or chain
    if (distance_index.is_root_snarl(current_handle)) {
        //FIrst thing is a snarl, so add the snarl's connected component number
        zipcode.add_value(0);
#ifdef DEBUG_ZIPCODE
        cerr << "Adding code for top-level snarl " << distance_index.net_handle_as_string(current_handle) << endl;
#endif
        zipcode.add_value(distance_index.get_connected_component_number(current_handle));
    } else {
        //FIrst thing is a chain so add its connected component number and remove the chain from the stack
        zipcode.add_value(1);

        //If the root-level structure is actually a chain, then save the connected component number and take out
        //the chain from the stack
        //If the root-level structure is a trivial chain, then just store the node (as a chain, which will have the 
        //connected-component number as the rank in the snarl anyways)
        if (!distance_index.is_trivial_chain(ancestors.back())) {
#ifdef DEBUG_ZIPCODE
        cerr << "Adding code for top-level chain" << endl;
#endif
            zipcode.add_value(distance_index.get_connected_component_number(ancestors.back()));
            ancestors.pop_back();
        }
    }

    //Go through the ancestors top (root) down and add them to the zip code
    //ancestors has everything but the root-level snarl/chain
    for (int i = ancestors.size()-1 ; i >= 0 ; i--) {
        net_handle_t current_ancestor = ancestors[i];
#ifdef DEBUG_ZIPCODE
        cerr << "Adding code for " << distance_index.net_handle_as_string(current_ancestor) << endl;
#endif
        if (distance_index.is_node(current_ancestor)) {
            vector<size_t> to_add = get_node_code(current_ancestor, distance_index);
            for (auto& x : to_add) {
                zipcode.add_value(x);
            }
#ifdef DEBUG_ZIPCODE
                assert(to_add.size() == ZipCode::NODE_SIZE);
#endif
        } else if (distance_index.is_chain(current_ancestor)) {
            vector<size_t> to_add = get_chain_code(current_ancestor, distance_index);
            for (auto& x : to_add) {
                zipcode.add_value(x);
            }
#ifdef DEBUG_ZIPCODE
                assert(to_add.size() == ZipCode::CHAIN_SIZE);
#endif
            if (distance_index.is_trivial_chain(current_ancestor)) {
                return;
            }
        } else if (distance_index.is_regular_snarl(current_ancestor)) {
            vector<size_t> to_add = get_regular_snarl_code(current_ancestor, ancestors[i-1], distance_index); 
            for (auto& x : to_add) {
                zipcode.add_value(x);
            }
#ifdef DEBUG_ZIPCODE
                assert(to_add.size() == ZipCode::REGULAR_SNARL_SIZE);
#endif
        } else {
#ifdef DEBUG_ZIPCODE
            assert(distance_index.is_snarl(current_ancestor));
#endif
            vector<size_t> to_add = get_irregular_snarl_code(current_ancestor, ancestors[i-1], distance_index);
#ifdef DEBUG_ZIPCODE
            assert(to_add.size() == ZipCode::IRREGULAR_SNARL_SIZE);
#endif
            for (auto& x : to_add) {
                zipcode.add_value(x);
            }
        }
    }
}

std::vector<size_t> ZipCode::to_vector() const {
    return zipcode.to_vector();
}

void ZipCode::from_vector(const std::vector<size_t>& values) {
    zipcode.from_vector(values);
}

ZipCodeDecoder::ZipCodeDecoder(const ZipCode* zipcode) :
    zipcode(zipcode), decoder(0) {
    fill_in_full_decoder();
}

void ZipCodeDecoder::fill_in_full_decoder() {
    if (zipcode->byte_count() == 0) {
        //If the zipcode is empty
        return;
    }
    bool done=false;
    while (!done) {
        done = fill_in_next_decoder();
    }
}

bool ZipCodeDecoder::fill_in_next_decoder() {
#ifdef DEBUG_ZIPCODE
    cerr << "Decode one more thing in the zipcode. Currently decoded " << decoder_length() << " things" << endl;
#endif
    
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
            std::tie(zip_value, zip_index) = zipcode->zipcode.get_value_and_next_index(zip_index);
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

        //Get the first value, which is 1 if the top-level structure is a chain
        for (size_t i = 0 ; i <= ZipCode::ROOT_IS_CHAIN_OFFSET ; i++) {
            std::tie(previous_is_chain, zip_index) = zipcode->zipcode.get_value_and_next_index(0);
        }
        //The next thing is the connected-component number
        for (size_t i = 0 ; i <= ZipCode::ROOT_IDENTIFIER_OFFSET - ZipCode::ROOT_IS_CHAIN_OFFSET -1; i++) {
            std::tie(zip_value, zip_index) = zipcode->zipcode.get_value_and_next_index(zip_index);
        }

        //If the top-level structure is a chain, it might actually be a node, in which case
        //the only other thing that got stored is the length
        if (previous_is_chain) {
            if (zipcode->zipcode.get_value_and_next_index(zip_index).second == std::numeric_limits<size_t>::max()) {
                //If the zip code ends here, then this was a node and we're done
#ifdef DEBUG_ZIPCODE
cerr << "\tThe last thing was a root-level node, so nothing else" << endl;
#endif
                return true;
            } else {
                //Otherwise, check if this is a node or a snarl. If it is a node, then there are three things remaining
                size_t start_index = zip_index;

                //If it's a node, then there are three remaining things in the index 
                //If it were a snarl, then there are more than three things
                for (size_t i = 0 ; i < ZipCode::NODE_SIZE ; i++) {
                    std::tie(zip_value, zip_index) = zipcode->zipcode.get_value_and_next_index(zip_index);
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
            decoder.emplace_back(!previous_is_chain, zip_index);
            return false;
        }
    } else {
        //If there was already stuff in the decoder, then figure out where the last thing
        //is and set values
        previous_is_chain = decoder.back().first;
        zip_index = decoder.back().second;
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
                check_zip_index = zipcode->zipcode.get_value_and_next_index(check_zip_index).second;
            }
            //If the zipcode ends after a chain
            if (check_zip_index == std::numeric_limits<size_t>::max()) {
#ifdef DEBUG_ZIPCODE
                cerr << "\tThe last thing was a chain pretending to be a node so we're done" << endl;
#endif
                return true;
            }
            //Now check if it was actually a real node
            for (size_t i = 0 ; i < std::max(ZipCode::NODE_SIZE, ZipCode::CHAIN_SIZE)
                                    - std::min(ZipCode::NODE_SIZE, ZipCode::CHAIN_SIZE); i++) {
                check_zip_index = zipcode->zipcode.get_value_and_next_index(check_zip_index).second;
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
                return true;
            } else {
                //Otherwise, the last thing was a chain
                //Get to the end of the chain
                for (size_t i = 0 ; i < ZipCode::CHAIN_SIZE  ; i++) {
                    std::tie(zip_value, zip_index) = zipcode->zipcode.get_value_and_next_index(zip_index);
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
                    check_zip_index = zipcode->zipcode.get_value_and_next_index(check_zip_index).second;
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
                std::tie(zip_value, zip_index) = zipcode->zipcode.get_value_and_next_index(zip_index);
            }

            if (zip_value == 1) {
#ifdef DEBUG_ZIPCODE
                cerr << "\tAdd a node child of a regular snarl" << endl;
#endif
                //Regular snarl, so 2 remaining things in the code
                for (size_t i = 0 ; i < ZipCode::REGULAR_SNARL_SIZE - ZipCode::SNARL_IS_REGULAR_OFFSET - 1; i++) {
                    std::tie(zip_value, zip_index) = zipcode->zipcode.get_value_and_next_index(zip_index);
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
                    std::tie(zip_value, zip_index) = zipcode->zipcode.get_value_and_next_index(zip_index);
                }
                decoder.emplace_back(!previous_is_chain, zip_index);
                return false;
            }
        }
    }    
}

size_t ZipCodeDecoder::max_depth() {
    return decoder_length()-1;

}

ZipCode::code_type_t ZipCodeDecoder::get_code_type(const size_t& depth) {

    //Now get the code type
    //A snarl is always a snarl. A chain could actually be a node
    if (depth == 0) {
        //If it is a root snarl/chain
        if (decoder[0].first) {
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
        if (decoder[depth].first) {
            //is_chain so could be a chain or a node
            if (decoder[depth-1].first) {
                //If the thing before this was also a chain, then it is a node
                return ZipCode::NODE;
            } else {
                //Otherwise it's a chain
                return ZipCode::CHAIN;
            }
        } else {
            //Definitely a snarl
            size_t zip_value;
            size_t zip_index = decoder[depth].second;
            for (size_t i = 0 ; i <= ZipCode::SNARL_IS_REGULAR_OFFSET; i++) {
                std::tie(zip_value, zip_index) = zipcode->zipcode.get_value_and_next_index(zip_index);
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

size_t ZipCodeDecoder::get_length(const size_t& depth, const SnarlDistanceIndex* distance_index) {

    if (depth == 0) {
        //If this is the root chain/snarl/node

        if (decoder_length() == 1) {
            //If the length is 1, then it's a node
            size_t zip_value;
            size_t zip_index = decoder[depth].second;
            for (size_t i = 0 ; i <= ZipCode::ROOT_NODE_LENGTH_OFFSET ; i++) {
                std::tie(zip_value, zip_index) = zipcode->zipcode.get_value_and_next_index(zip_index);
            }
            return zip_value == 0 ? std::numeric_limits<size_t>::max() : zip_value-1;

        } else {
            //Otherwise, we didn't store the length
            throw std::runtime_error("zipcodes don't store lengths of top-level chains or snarls");
        }
    } else if (decoder[depth].first) {
        //If this is a chain/node

        //If this is a chain or a node, then the length will be the second thing
        size_t zip_value;
        size_t zip_index = decoder[depth].second;

        for (size_t i = 0 ; i <= ZipCode::CHAIN_LENGTH_OFFSET ; i++) {
            std::tie(zip_value, zip_index) = zipcode->zipcode.get_value_and_next_index(zip_index);
        }
        return zip_value == 0 ? std::numeric_limits<size_t>::max() : zip_value-1;
    } else {
        //If this is a snarl

        size_t zip_value;
        size_t zip_index = decoder[depth].second;

        for (size_t i = 0 ; i <= ZipCode::SNARL_LENGTH_OFFSET ; i++) {
            std::tie(zip_value, zip_index) = zipcode->zipcode.get_value_and_next_index(zip_index);
        }

        return zip_value == 0 ? std::numeric_limits<size_t>::max() : zip_value-1;
    }
}

size_t ZipCodeDecoder::get_rank_in_snarl(const size_t& depth) {


    if (depth == 0) {
        //If this is the root chain/snarl/node
        throw std::runtime_error("zipcodes don't store ranks of top-level chains or snarls");

    } else if (decoder[depth].first) {
        //If this is a chain/node

        if (decoder[depth-1].first) {
            throw std::runtime_error("zipcodes trying to find the rank in snarl of a node in a chain");
        }

        size_t zip_value;
        size_t zip_index = decoder[depth].second;
        for (size_t i = 0 ; i <= ZipCode::CHAIN_RANK_IN_SNARL_OFFSET ; i++) {
            std::tie(zip_value, zip_index) = zipcode->zipcode.get_value_and_next_index(zip_index);
        }
        return zip_value;
    } else {
        //If this is a snarl
        throw std::runtime_error("zipcodes don't store snarl ranks for snarls");
    }
}

size_t ZipCodeDecoder::get_snarl_child_count(const size_t& depth, const SnarlDistanceIndex* distance_index) {


    if (depth == 0) {
        //TODO: This could be actually saved in the zipcode but I'll have to go to the distance index anyway
        assert(distance_index != nullptr);
        size_t child_count = 0;
        distance_index->for_each_child(get_net_handle(depth, distance_index), [&] (const net_handle_t& child) {
            child_count++;
        });
        return child_count;

    } else if (!decoder[depth].first) {
        //If this is a snarl

        size_t zip_value;
        size_t zip_index = decoder[depth].second;
        for (size_t i = 0 ; i <= ZipCode::SNARL_CHILD_COUNT_OFFSET ; i++) {
            std::tie(zip_value, zip_index) = zipcode->zipcode.get_value_and_next_index(zip_index);
        }
        return zip_value;
    } else {
        //If this is not a snarl
        throw std::runtime_error("trying to get the snarl child count of a non-snarl zipcode");
    }
}

size_t ZipCodeDecoder::get_offset_in_chain(const size_t& depth, const SnarlDistanceIndex* distance_index) {


    if (depth == 0) {
        //If this is the root chain/snarl/node
        throw std::runtime_error("zipcodes don't have chain offsets for roots");

    } else if (decoder[depth].first) {
        //If this is a chain/node

        if (!decoder[depth-1].first) {
            throw std::runtime_error("zipcodes trying to find the offset in child of a snarl");
        }
        size_t zip_value;
        size_t zip_index = decoder[depth].second;
        for (size_t i = 0 ; i <= ZipCode::CHAIN_RANK_IN_SNARL_OFFSET ; i++) {
            std::tie(zip_value, zip_index) = zipcode->zipcode.get_value_and_next_index(zip_index);
        }

        return zip_value == std::numeric_limits<size_t>::max() ? 0 : zip_value-1;
    } else {
        //If this is a snarl

        size_t zip_value;
        size_t zip_index = decoder[depth].second;
        for (size_t i = 0 ; i <= ZipCode::SNARL_OFFSET_IN_CHAIN_OFFSET ; i++) {
                std::tie(zip_value, zip_index) = zipcode->zipcode.get_value_and_next_index(zip_index);
        }

        return zip_value == 0 ? std::numeric_limits<size_t>::max() : zip_value-1;
    }
}
bool ZipCodeDecoder::get_is_reversed_in_parent(const size_t& depth) {


    if (depth == 0) {
        //If this is the root chain/snarl/node
        return false;

    } else if (decoder[depth].first) {
        //If this is a chain/node

        if (decoder[depth-1].first) {
            //If the parent is a chain, then this is a node and we need to check its orientation

            size_t zip_value;
            size_t zip_index = decoder[depth].second;
            for (size_t i = 0 ; i <= ZipCode::NODE_IS_REVERSED_OFFSET ; i++) {
                std::tie(zip_value, zip_index) = zipcode->zipcode.get_value_and_next_index(zip_index);
            }
            return zip_value;
        } else {
            //If the parent is a snarl, then this might be a chain in a regular snarl
            size_t zip_value;
            size_t zip_index = decoder[depth-1].second;
            //zip_value is true if the parent is a regular snarl
            for (size_t i = 0 ; i <= ZipCode::SNARL_IS_REGULAR_OFFSET ; i++) {
                std::tie(zip_value, zip_index) = zipcode->zipcode.get_value_and_next_index(zip_index);
            }
            if (zip_value == 1) {
                //The parent is a regular snarl, which stores is_reversed for the child
                
                for (size_t i = 0 ; i <= ZipCode::REGULAR_SNARL_IS_REVERSED_OFFSET -
                                         ZipCode::SNARL_IS_REGULAR_OFFSET - 1 ; i++) {
                    std::tie(zip_value, zip_index) = zipcode->zipcode.get_value_and_next_index(zip_index);
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

net_handle_t ZipCodeDecoder::get_net_handle(const size_t& depth, const SnarlDistanceIndex* distance_index) {


    if (depth == 0) {
        //If this is the root chain/snarl/node

        size_t zip_value, zip_index = 0;
        for (size_t i = 0 ; i <= ZipCode::ROOT_IDENTIFIER_OFFSET ; i++) {
            std::tie(zip_value, zip_index) = zipcode->zipcode.get_value_and_next_index(zip_index);
        }
        return distance_index->get_handle_from_connected_component(zip_value);

    } else if (decoder[depth].first) {
        //If this is a chain/node

        throw std::runtime_error("zipcodes trying to get a handle of a chain or node");
    } else {
        //If this is a snarl

        size_t zip_value; 
        size_t zip_index = decoder[depth].second;
        //zip_value is is_regular_snarl
        for (size_t i = 0 ; i <= ZipCode::SNARL_IS_REGULAR_OFFSET ; i++) {
            std::tie(zip_value, zip_index) = zipcode->zipcode.get_value_and_next_index(zip_index);
        }
        if (zip_value == 1) {
            //If this is a regular snarl

            throw std::runtime_error("zipcodes trying to get a handle of a regular snarl");
        } else {
            //Irregular snarl

            //zip_value is distance index offset
            for (size_t i = 0 ; i <= ZipCode::IRREGULAR_SNARL_RECORD_OFFSET -
                                     ZipCode::SNARL_IS_REGULAR_OFFSET-1 ; i++) {
                std::tie(zip_value, zip_index) = zipcode->zipcode.get_value_and_next_index(zip_index);
            }
            net_handle_t snarl_handle = distance_index->get_net_handle_from_values(zip_value, SnarlDistanceIndex::START_END, SnarlDistanceIndex::SNARL_HANDLE);
            return snarl_handle;
        }
    }
}

size_t ZipCodeDecoder::get_distance_index_address(const size_t& depth) {


    if (depth == 0) {
        //If this is the root chain/snarl/node

        size_t zip_value, zip_index = 0;
        for (size_t i = 0 ; i <= ZipCode::ROOT_IDENTIFIER_OFFSET ; i++) {
            std::tie(zip_value, zip_index) = zipcode->zipcode.get_value_and_next_index(zip_index);
        }
        return zip_value;

    } else if (decoder[depth].first) {
        //If this is a chain/node

        throw std::runtime_error("zipcodes trying to get a handle of a chain or node");
    } else {
        //If this is a snarl

        size_t zip_value;
        size_t zip_index = decoder[depth].second;
        //zip_value is is_regular_snarl
        for (size_t i = 0 ; i <= ZipCode::SNARL_IS_REGULAR_OFFSET ; i++) {
            std::tie(zip_value, zip_index) = zipcode->zipcode.get_value_and_next_index(zip_index);
        }
        if (zip_value == 1) {
            //If this is a regular snarl

            throw std::runtime_error("zipcodes trying to get a handle of a regular ansl");
        } else {
            //Irregular snarl

            //zip_value is distance index offset
            for (size_t i = 0 ; i <= ZipCode::IRREGULAR_SNARL_RECORD_OFFSET -
                                     ZipCode::SNARL_IS_REGULAR_OFFSET-1 ; i++) {
                std::tie(zip_value, zip_index) = zipcode->zipcode.get_value_and_next_index(zip_index);
            }
            return zip_value;
        }
    }
}
size_t ZipCodeDecoder::get_distance_to_snarl_bound(const size_t& depth, bool snarl_start, bool left_side) {

#ifdef DEBUG_ZIPCODE
    assert(depth > 0);
    assert((get_code_type(depth-1) == ZipCode::IRREGULAR_SNARL || get_code_type(depth-1) == ZipCode::REGULAR_SNARL || get_code_type(depth-1) == ZipCode::CYCLIC_SNARL)); 
#endif
     size_t zip_value;
     size_t zip_index = decoder[depth-1].second;
     //zip_value is 1 if the parent is a regular snarl
     for (size_t i = 0 ; i <= ZipCode::SNARL_IS_REGULAR_OFFSET ; i++) {
         std::tie(zip_value, zip_index) = zipcode->zipcode.get_value_and_next_index(zip_index);
     }
     if (zip_value == 1) {
         //The parent is a regular snarl, which stores is_reversed for the child
         for (size_t i = 0 ; i <= ZipCode::REGULAR_SNARL_IS_REVERSED_OFFSET -
                                  ZipCode::SNARL_IS_REGULAR_OFFSET - 1 ; i++) {
             std::tie(zip_value, zip_index) = zipcode->zipcode.get_value_and_next_index(zip_index);
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
            std::tie(zip_value, zip_index) = zipcode->zipcode.get_value_and_next_index(zip_index);
        }
        return zip_value == 0 ? std::numeric_limits<size_t>::max() : zip_value - 1;
     }
}

const bool ZipCodeDecoder::is_equal(ZipCodeDecoder& decoder1, ZipCodeDecoder& decoder2,
                                        const size_t& depth) {

    if (decoder1.max_depth() < depth && decoder2.max_depth() < depth ) {
        return false;
    }

    //First, check if the code types are the same
    ZipCode::code_type_t type1 = decoder1.get_code_type(depth);
    ZipCode::code_type_t type2 = decoder2.get_code_type(depth);
    if (type1 != type2) {
        return false;
    }

    if (type1 == ZipCode::ROOT_NODE || type1 == ZipCode::ROOT_CHAIN || type1 == ZipCode::ROOT_SNARL || type1 == ZipCode::IRREGULAR_SNARL || type1 == ZipCode::CYCLIC_SNARL ) {
        //If the codes are for root-structures or irregular/cyclic snarls, just check if the 
        //connected component numbers are the same
        return decoder1.get_distance_index_address(depth) == decoder2.get_distance_index_address(depth);
    } else {
        //Check the parent type. If the parent is a snarl, then check rank. If it's a chain,
        //then check the prefix sum
        if (decoder1.get_code_type(depth-1) == ZipCode::REGULAR_SNARL ||
            decoder1.get_code_type(depth-1) == ZipCode::IRREGULAR_SNARL ||
            decoder1.get_code_type(depth-1) == ZipCode::CYCLIC_SNARL ||
            decoder1.get_code_type(depth-1) == ZipCode::ROOT_SNARL) {
            //If the parent is a snarl, then check the rank
            return decoder1.get_rank_in_snarl(depth) == decoder2.get_rank_in_snarl(depth);
        } else {
            //Otherwise, check the offset in the chain
            //Since the type is the same, this is sufficient
            return decoder1.get_offset_in_chain(depth) == decoder2.get_offset_in_chain(depth);
        }
    }
}

void ZipCodeDecoder::dump(std::ostream& out) const {
    if (!zipcode) {
        // We're decoding nothing
        out << *this;
    } else {
        std::vector<size_t> numbers = zipcode->to_vector();
        // Print out the numbers in a way that is easy to copy-paste as a vector literal.
        out << "<decoder for {";
        for (size_t i = 0; i < numbers.size(); i++) {
            out << numbers[i];
            if (i + 1 < numbers.size()) {
                out << ", ";
            }
        }
        out << "}>";
    }
}

std::ostream& operator<<(std::ostream& out, const ZipCodeDecoder& decoder) {
    return out << "<decoder for " << decoder.zipcode << ">";
}


vector<size_t> ZipCode::get_node_code(const net_handle_t& node, const SnarlDistanceIndex& distance_index) {
#ifdef DEBUG_ZIPCODE
    assert(!distance_index.is_trivial_chain(node));
    assert((distance_index.is_chain(distance_index.get_parent(node)) || distance_index.is_root(distance_index.get_parent(node))));
#endif
    //Node code is: offset in chain, length, is reversed
    vector<size_t> node_code;
    //Assume this node is in a regular chain
    size_t prefix_sum = distance_index.get_prefix_sum_value(node); 
    node_code.emplace_back(prefix_sum == std::numeric_limits<size_t>::max() ? 0 : prefix_sum+1);
    node_code.emplace_back(distance_index.minimum_length(node)+1);
    node_code.emplace_back(distance_index.is_reversed_in_parent(node));
    return node_code;

}
vector<size_t> ZipCode::get_chain_code(const net_handle_t& chain, const SnarlDistanceIndex& distance_index) {
    //Chain code is: rank in snarl, length
    vector<size_t> chain_code;
    chain_code.emplace_back(distance_index.get_rank_in_parent(chain));
    size_t len = distance_index.minimum_length(chain);
    chain_code.emplace_back(len == std::numeric_limits<size_t>::max() ? 0 : len+1);
    return chain_code;

}
vector<size_t> ZipCode::get_regular_snarl_code(const net_handle_t& snarl, const net_handle_t& snarl_child, const SnarlDistanceIndex& distance_index) {
    //Regular snarl code is 1, offset in chain, length, is reversed
    vector<size_t> snarl_code (REGULAR_SNARL_SIZE);

    //Tag to say that it's a regular snarl
    snarl_code[SNARL_IS_REGULAR_OFFSET] = 1;

    //The number of children
    size_t child_count = 0;
    distance_index.for_each_child(snarl, [&] (const net_handle_t& child) {
        child_count++;
    });
    snarl_code[SNARL_CHILD_COUNT_OFFSET] = child_count;

    //Chain prefix sum value for the start of the snarl, which is the prefix sum of the start node + length of the start node
    net_handle_t start_node = distance_index.get_node_from_sentinel(distance_index.get_bound(snarl, false, false));
    size_t prefix_sum = SnarlDistanceIndex::sum(distance_index.get_prefix_sum_value(start_node), distance_index.minimum_length(start_node));
    snarl_code[SNARL_OFFSET_IN_CHAIN_OFFSET] = (prefix_sum == std::numeric_limits<size_t>::max() ? 0 : prefix_sum+1);

    //Length of the snarl
    size_t len = distance_index.minimum_length(snarl);
    snarl_code[SNARL_LENGTH_OFFSET] = (len == std::numeric_limits<size_t>::max() ? 0 : len+1);

    //Is the child of the snarl reversed in the snarl
#ifdef DEBUG_ZIPCODE
    assert(distance_index.is_chain(snarl_child));
#endif
    snarl_code[REGULAR_SNARL_IS_REVERSED_OFFSET] = (distance_index.distance_in_parent(snarl, 
                                                        distance_index.get_bound(snarl, false, true),
                                                        distance_index.flip(distance_index.canonical(snarl_child))) != 0);

    return snarl_code;

}
vector<size_t> ZipCode::get_irregular_snarl_code(const net_handle_t& snarl, const net_handle_t& snarl_child, 
                                                 const SnarlDistanceIndex& distance_index) {
    vector<size_t> snarl_code (IRREGULAR_SNARL_SIZE);

    //Tag to say that it's an irregular snarl
    snarl_code[SNARL_IS_REGULAR_OFFSET] = distance_index.is_dag(snarl) ? 0 : 2;

    //The number of children
    size_t child_count = 0;
    distance_index.for_each_child(snarl, [&] (const net_handle_t& child) {
        child_count++;
    });
    snarl_code[SNARL_CHILD_COUNT_OFFSET] = child_count;

    //Chain prefix sum value for the start of the snarl, which is the prefix sum of the start node + length of the start node
    net_handle_t start_node = distance_index.get_node_from_sentinel(distance_index.get_bound(snarl, false, false));
    size_t prefix_sum = SnarlDistanceIndex::sum(distance_index.get_prefix_sum_value(start_node), distance_index.minimum_length(start_node));
    snarl_code[SNARL_OFFSET_IN_CHAIN_OFFSET] = (prefix_sum == std::numeric_limits<size_t>::max() ? 0 : prefix_sum+1);

    //Length of the snarl
    size_t len = distance_index.minimum_length(snarl);
    snarl_code[SNARL_LENGTH_OFFSET] = (len == std::numeric_limits<size_t>::max() ? 0 : len+1);


    //Record offset to look up distances in the index later
    snarl_code[IRREGULAR_SNARL_RECORD_OFFSET] = (distance_index.get_record_offset(snarl));

    snarl_code[IRREGULAR_SNARL_DISTANCE_LEFT_START_OFFSET] = distance_index.distance_to_parent_bound(snarl, true, distance_index.flip(snarl_child));
    snarl_code[IRREGULAR_SNARL_DISTANCE_LEFT_END_OFFSET] = distance_index.distance_to_parent_bound(snarl, false, distance_index.flip(snarl_child));
    snarl_code[IRREGULAR_SNARL_DISTANCE_RIGHT_START_OFFSET] = distance_index.distance_to_parent_bound(snarl, true, snarl_child);
    snarl_code[IRREGULAR_SNARL_DISTANCE_RIGHT_END_OFFSET] = distance_index.distance_to_parent_bound(snarl, false, snarl_child);

    //Add 1 to values to store inf properly
    snarl_code[IRREGULAR_SNARL_DISTANCE_LEFT_START_OFFSET] = 
        snarl_code[IRREGULAR_SNARL_DISTANCE_LEFT_START_OFFSET] == std::numeric_limits<size_t>::max() 
        ? 0 
        : snarl_code[IRREGULAR_SNARL_DISTANCE_LEFT_START_OFFSET] + 1;
    snarl_code[IRREGULAR_SNARL_DISTANCE_RIGHT_START_OFFSET] = 
        snarl_code[IRREGULAR_SNARL_DISTANCE_RIGHT_START_OFFSET] == std::numeric_limits<size_t>::max() 
        ? 0 
        : snarl_code[IRREGULAR_SNARL_DISTANCE_RIGHT_START_OFFSET] + 1;
    snarl_code[IRREGULAR_SNARL_DISTANCE_LEFT_END_OFFSET] = 
        snarl_code[IRREGULAR_SNARL_DISTANCE_LEFT_END_OFFSET] == std::numeric_limits<size_t>::max() 
        ? 0 
        : snarl_code[IRREGULAR_SNARL_DISTANCE_LEFT_END_OFFSET] + 1;
    snarl_code[IRREGULAR_SNARL_DISTANCE_RIGHT_END_OFFSET] = 
        snarl_code[IRREGULAR_SNARL_DISTANCE_RIGHT_END_OFFSET] == std::numeric_limits<size_t>::max() 
        ? 0 
        : snarl_code[IRREGULAR_SNARL_DISTANCE_RIGHT_END_OFFSET] + 1;

    return snarl_code;

}

size_t ZipCode::minimum_distance_between(ZipCodeDecoder& zip1_decoder, const pos_t& pos1,   
    ZipCodeDecoder& zip2_decoder, const pos_t& pos2, const SnarlDistanceIndex& distance_index,
    size_t distance_limit, bool undirected_distance, const HandleGraph* graph){


#ifdef DEBUG_ZIPCODE
//Make sure that the zip codes actually correspond to the positions
    ZipCode check_zip1;
    check_zip1.fill_in_zipcode(distance_index, pos1);
    assert(*zip1_decoder.zipcode == check_zip1);

    ZipCode check_zip2;
    check_zip2.fill_in_zipcode(distance_index, pos2);
    assert(*zip2_decoder.zipcode == check_zip2);

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
    auto update_distances_to_ends_of_parent = [&] (ZipCodeDecoder& decoder, const size_t& child_depth, 
                                            size_t& distance_to_start, size_t& distance_to_end) {
#ifdef DEBUG_ZIPCODE
        cerr << "Update distance to ends of parent at depth " << child_depth << endl;
#endif
        //The distances from the start/end of current child to the start/end(left/right) of the parent
        size_t distance_start_left = std::numeric_limits<size_t>::max();
        size_t distance_start_right = std::numeric_limits<size_t>::max();
        size_t distance_end_left = std::numeric_limits<size_t>::max();
        size_t distance_end_right = std::numeric_limits<size_t>::max();

        code_type_t parent_type = decoder.get_code_type(child_depth-1);

        if (parent_type == IRREGULAR_SNARL || parent_type == CYCLIC_SNARL) {
            //If the parent is an irregular snarl
            net_handle_t parent_handle = decoder.get_net_handle(child_depth-1, &distance_index);
            size_t child_rank = decoder.get_rank_in_snarl(child_depth);
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
            //If its a regular snarl, then the distances to the ends are either 0 or inf
            //For a regular snarl, the snarl stores if the child was reversed, rather than the child
            if (decoder.get_is_reversed_in_parent(child_depth)) {
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
            if (decoder.get_code_type(child_depth) == NODE && 
                decoder.get_is_reversed_in_parent(child_depth)){ 
                //If this is reversed in the chain

                distance_start_left = std::numeric_limits<size_t>::max();
                distance_end_right = std::numeric_limits<size_t>::max();
                //Prefix sum of the child
                distance_end_left = decoder.get_offset_in_chain(child_depth, &distance_index);
                //Length of the chain - prefix sum of the child - length of the child
                distance_start_right = SnarlDistanceIndex::minus(SnarlDistanceIndex::minus(
                        decoder.get_length(child_depth-1, &distance_index), 
                        decoder.get_offset_in_chain(child_depth, &distance_index)), 
                        decoder.get_length(child_depth, &distance_index));
            } else {
                //If it is a node that isn't reversed in the chain, or it's a snarl which is never reversed
                distance_end_left = std::numeric_limits<size_t>::max();
                distance_start_right = std::numeric_limits<size_t>::max();
                //Prefix sum of the child
                distance_start_left = decoder.get_offset_in_chain(child_depth, &distance_index);
                //Length of the chain - prefix sum of the child - length of the child
                distance_end_right = SnarlDistanceIndex::minus(SnarlDistanceIndex::minus(
                        decoder.get_length(child_depth-1, &distance_index), 
                        decoder.get_offset_in_chain(child_depth, &distance_index)), 
                        decoder.get_length(child_depth, &distance_index));
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

    if (zip1_decoder.get_distance_index_address(0) != zip2_decoder.get_distance_index_address(0)) {
#ifdef DEBUG_ZIPCODE
        cerr << "Zip codes are on different connected components" << endl;
#endif
        return std::numeric_limits<size_t>::max();
    }

    //The two positions are in the same connected component so now fill in the rest
    //of the decoder and try to find the distance
    zip1_decoder.fill_in_full_decoder();
    zip2_decoder.fill_in_full_decoder();

    //Now find the lowest common ancestor of the two zipcodes
    size_t lowest_common_ancestor_depth = 0;
    bool still_equal = true;
    while (still_equal) {

        if (lowest_common_ancestor_depth == zip1_decoder.decoder_length()-1 ||
            lowest_common_ancestor_depth == zip2_decoder.decoder_length()-1 ||
            !ZipCodeDecoder::is_equal(zip1_decoder, zip2_decoder, 
                                         lowest_common_ancestor_depth+1)) {
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
        lowest_common_ancestor_depth < zip1_decoder.decoder_length()-1){
        //If we're aborting when the distance is definitely too far,
        code_type_t ancestor_type = zip1_decoder.get_code_type(lowest_common_ancestor_depth);
        if  (ancestor_type == CHAIN || ancestor_type == ROOT_CHAIN) {
            //If the current ancestor is a chain, then check the distance
            size_t prefix_sum1 = zip1_decoder.get_offset_in_chain(lowest_common_ancestor_depth+1, &distance_index);
            size_t prefix_sum2 = zip2_decoder.get_offset_in_chain(lowest_common_ancestor_depth+1, &distance_index);
            size_t distance_in_chain; 
            if (prefix_sum1 < prefix_sum2) {
                //zip1 comes before zip2
                distance_in_chain = SnarlDistanceIndex::minus(
                    prefix_sum2, 
                    SnarlDistanceIndex::sum(prefix_sum1, 
                                            zip1_decoder.get_length(lowest_common_ancestor_depth+1, &distance_index)));
            } else {
                //zip2 comes before zip1
                distance_in_chain = SnarlDistanceIndex::minus(
                    prefix_sum1, 
                    SnarlDistanceIndex::sum(prefix_sum2, 
                                            zip2_decoder.get_length(lowest_common_ancestor_depth+1, &distance_index)));
            }
            if (distance_in_chain > distance_limit) {
                return std::numeric_limits<size_t>::max();
            }
        }
    }

    //Start from the nodes
    size_t distance_to_start1 = is_rev(pos1) 
        ? zip1_decoder.get_length(zip1_decoder.decoder_length()-1, &distance_index) - offset(pos1) 
        : offset(pos1) + 1;
    size_t distance_to_end1 = is_rev(pos1) ? offset(pos1) + 1 
         : zip1_decoder.get_length(zip1_decoder.decoder_length()-1, &distance_index) - offset(pos1);
    size_t distance_to_start2 = is_rev(pos2) 
         ? zip2_decoder.get_length(zip2_decoder.decoder_length()-1, &distance_index) - offset(pos2) 
         : offset(pos2) + 1;
    size_t distance_to_end2 = is_rev(pos2) ? offset(pos2) + 1 
         : zip2_decoder.get_length(zip2_decoder.decoder_length()-1, &distance_index) - offset(pos2);

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
    for (int i = zip1_decoder.decoder_length()-2 ; i > 0 && i > lowest_common_ancestor_depth ; i--) {
        //the parent snarl tree node is at index i
        //The distances are currently to the ends of the current node
        //FInd the distances to the ends of the parent
        update_distances_to_ends_of_parent(zip1_decoder, i+1, distance_to_start1, distance_to_end1);
    }
#ifdef DEBUG_ZIPCODE
cerr << "Finding distances to ancestors of second position" << endl;
#endif
    //The same thing for the second position
    for (int i = zip2_decoder.decoder_length()-2 ; i > 0 && i > lowest_common_ancestor_depth ; i--) {
        //the parent snarl tree node is at index i
        //The distances are currently to the ends of the current node
        //FInd the distances to the ends of the parent

        update_distances_to_ends_of_parent(zip2_decoder, i+1, distance_to_start2, distance_to_end2);
    }


    //Distances are now the distances to the ends of a child of the common ancestor

#ifdef DEBUG_ZIPCODE
    cerr << "Distances in children of common ancestor: " << distance_to_start1 << " " << distance_to_end1 << " " << distance_to_start2 << " " << distance_to_end2 << endl;
    //Check that the current nodes are actually children of the lca
    assert(ZipCodeDecoder::is_equal(zip1_decoder, zip2_decoder, lowest_common_ancestor_depth));
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
        if (depth == zip1_decoder.decoder_length()-1) {
            //If the lca is a node that both positions are on

#ifdef DEBUG_ZIPCODE
            //If the lca is a node, then both the zipcode nodes should be the same node
            assert(ZipCodeDecoder::is_equal(zip1_decoder, zip2_decoder, depth));
            assert(depth == zip2_decoder.decoder_length()-1);
            cerr << "\tAncestor should be a node" << endl;
#endif
            size_t d1 = SnarlDistanceIndex::sum(distance_to_end1, distance_to_start2);
            size_t d2 = SnarlDistanceIndex::sum(distance_to_end2, distance_to_start1);
            size_t node_length = zip1_decoder.get_length(depth, &distance_index);
            if (d1 > node_length) {
                distance_between = std::min(distance_between,
                                            SnarlDistanceIndex::minus(SnarlDistanceIndex::minus(d1, node_length),1));
            } 
            if (d2 > node_length) {
                distance_between = std::min(distance_between,
                                            SnarlDistanceIndex::minus(SnarlDistanceIndex::minus(d2, node_length),1));
            }
        } else if ( zip1_decoder.decoder[depth].first) {
#ifdef DEBUG_ZIPCODE
            cerr << "\tancestor should be a chain" << endl;
#endif
            //If this ancestor is a chain

            //If the children are reversed in the chain, then flip their distances
            bool rev1 = (zip1_decoder.get_code_type(depth+1) == NODE && 
                zip1_decoder.get_is_reversed_in_parent(depth+1));
            size_t dist_start1 = rev1 ? distance_to_end1 : distance_to_start1;
            size_t dist_end1 = rev1 ? distance_to_start1 : distance_to_end1;

            bool rev2 = zip2_decoder.get_code_type(depth+1) == NODE && 
                zip2_decoder.get_is_reversed_in_parent(depth+1);
            size_t dist_start2 = rev2 ? distance_to_end2 : distance_to_start2;
            size_t dist_end2 = rev2 ? distance_to_start2 : distance_to_end2;

            //If they are the same child, then there is no path between them in the chain because we don't allow loops
            //So first check that they aren't the same
            if (!(ZipCodeDecoder::is_equal(zip1_decoder, zip2_decoder, depth+1) 
                )){//TODO: I think this is unnecessary || (zip1_decoder.get_code_type(depth+1) == NODE && id(pos1) == id(pos2)))) 
                size_t prefix_sum1 = zip1_decoder.get_offset_in_chain(depth+1, &distance_index);
                size_t prefix_sum2 = zip2_decoder.get_offset_in_chain(depth+1, &distance_index);
                code_type_t code_type1 = zip1_decoder.get_code_type(depth+1);
                code_type_t code_type2 = zip2_decoder.get_code_type(depth+1);

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
                        cerr << "Find distances from : " << prefix_sum2 << " " << dist_start2 << " " << prefix_sum1 << " " << zip1_decoder.get_length(depth+1, &distance_index)  << " " << dist_end1 << endl;
#endif
                        if (dist_start2 != std::numeric_limits<size_t>::max()
                            && dist_end1 != std::numeric_limits<size_t>::max()) {
                            distance_between = std::min(distance_between,
                                                        SnarlDistanceIndex::minus(SnarlDistanceIndex::sum(
                                                            SnarlDistanceIndex::minus(
                                                                SnarlDistanceIndex::sum(prefix_sum2, 
                                                                                        dist_start2), 
                                                                SnarlDistanceIndex::sum(prefix_sum1,
                                                                                        zip1_decoder.get_length(depth+1, &distance_index))),
                                                             dist_end1),1));
                        }
                    } else {
                        //Otherwise, all that matters is the prefix sums
                        //(Prefix sum 2  + distance left 2) - (prefix sum1+ length 1) + distance right 1
#ifdef DEBUG_ZIPCODE
                        cerr << "First child comes first in the chain and it isn't a snarl" << endl;
                        cerr << "Find distances from : " << prefix_sum2 << " " << dist_start2 << " " << prefix_sum1 << " " << dist_end1 << " " << zip1_decoder.get_length(depth+1, &distance_index) << endl;
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
                                                                                    zip1_decoder.get_length(depth+1, &distance_index))), 

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
                        cerr << "Find distances from : " << prefix_sum1 << " " << dist_start1 << " " << prefix_sum2 << " " << zip2_decoder.get_length(depth+1, &distance_index)  << " " << dist_end2 << endl;
#endif
                        if (dist_start1 != std::numeric_limits<size_t>::max() 
                             && dist_end2 != std::numeric_limits<size_t>::max() ){
                            distance_between = std::min(distance_between,
                                                        SnarlDistanceIndex::minus(SnarlDistanceIndex::sum(
                                                            SnarlDistanceIndex::minus(
                                                                SnarlDistanceIndex::sum(prefix_sum1, 
                                                                                        dist_start1), 
                                                                SnarlDistanceIndex::sum(prefix_sum2,
                                                                                        zip2_decoder.get_length(depth+1, &distance_index))),
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
                                                                                    zip2_decoder.get_length(depth+1, &distance_index))), 

                                                            dist_end2),1) );
                        }
                    }
                }
            }
            //Update distances from the ends of the children (at depth+1) to parent (depth)
            update_distances_to_ends_of_parent(zip1_decoder, depth+1, distance_to_start1, distance_to_end1);
            update_distances_to_ends_of_parent(zip2_decoder, depth+1, distance_to_start2, distance_to_end2);
        } else {

#ifdef DEBUG_ZIPCODE
            cerr << "\tancestor is a snarl" << endl;
#endif
            //If the ancestor is a snarl
            
            //If the parent is a regular snarl, then there is no path between them so
            //just update the distances to the ends of the parent 
            if (zip1_decoder.get_code_type(depth) != REGULAR_SNARL) {
                //Parent may be an irregular snarl or a root snarl (which is also irregular)
                net_handle_t parent_handle = zip1_decoder.get_net_handle(depth, &distance_index);
                size_t rank1 = zip1_decoder.get_rank_in_snarl(depth+1);
                size_t rank2 = zip2_decoder.get_rank_in_snarl(depth+1);
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
            update_distances_to_ends_of_parent(zip1_decoder, depth+1, distance_to_start1, distance_to_end1);
            update_distances_to_ends_of_parent(zip2_decoder, depth+1, distance_to_start2, distance_to_end2);
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

    }

    //Both zips now point to a thing in a shared chain
    //Get the minimum possible distance between the structures on the chain
    //For a lower bound, this assumes that the positions are as close as they can be on the structure in the chain
    size_t prefix_sum1, prefix_sum2, length1, length2;

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

    if (prefix_sum1 == std::numeric_limits<size_t>::max() ||
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
    if (byte_count() > 15) {
        //If there aren't enough bits to represent the zip code
        return MIPayload::NO_CODE;
    }
    
    //Index and value as we walk through the zip code
    size_t index = 0;
    size_t value;
    
    //The values that get returned
    code_type encoded1 = 0;
    code_type encoded2 = 0;

    encoded1 |= byte_count();

    for (size_t i = 0 ; i < zipcode.data.size() ; i++ ) {
       size_t byte = static_cast<size_t> (zipcode.data[i]); 
       if ( i < 7 ) {
            //Add to first code
            encoded1 |= (byte << ((i+1)*8));

        } else {
            //Add to second code
            encoded2 |= (byte << ((i-7)*8));
        }
    
    }
    return {encoded1, encoded2};

}

void ZipCode::fill_in_zipcode_from_payload(const gbwtgraph::Payload& payload) {
    assert(payload != MIPayload::NO_CODE);

    //get one byte at a time from the payload and add it to the zip code
    size_t bit_mask = (1 << 8) - 1;
    size_t byte_count = payload.first & bit_mask;
    for (size_t i = 1 ; i <= byte_count ; i++) {
        if (i < 8) {
            zipcode.add_one_byte((payload.first >> (i*8)) & bit_mask);
        } else {
            zipcode.add_one_byte((payload.second >> ((i-8)*8)) & bit_mask);
        }

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
        zipcodes.emplace_back(std::move(zip));
    }

}

size_t MIPayload::record_offset(const ZipCode& code, const SnarlDistanceIndex& distance_index, const nid_t& id ) {

    //TODO: This is pointless but I'll keep it until I fix everything
    net_handle_t node_handle = distance_index.get_node_net_handle(id);
    return distance_index.get_record_offset(node_handle);
}

size_t MIPayload::parent_record_offset(const ZipCode& zip, const SnarlDistanceIndex& distance_index, const nid_t& id) {

 
    ZipCodeDecoder decoder (&zip);

    bool root_is_chain = decoder.get_code_type(0) != ZipCode::ROOT_SNARL;

    if (decoder.decoder_length() == 1) {
        //If the root-level structure is a node

        return  0;
    } else if (decoder.decoder_length() == 2 && root_is_chain) {
        //If this is a node in the top-level chain
#ifdef DEBUG_ZIPCODE
        assert(distance_index.get_record_offset(decoder.get_net_handle(0, &distance_index)) ==
              distance_index.get_record_offset(distance_index.get_parent(distance_index.get_node_net_handle(id))));
#endif

        return distance_index.get_record_offset(decoder.get_net_handle(0, &distance_index));
    
    } else if (decoder.decoder_length() == 2 && !root_is_chain) {
        //If the node is the child of the root snarl
#ifdef DEBUG_ZIPCODE
        assert(distance_index.get_record_offset(decoder.get_net_handle(0, &distance_index)) ==
              distance_index.get_record_offset(distance_index.get_parent(distance_index.get_parent(distance_index.get_node_net_handle(id)))));
#endif
        
        return  distance_index.get_record_offset(decoder.get_net_handle(0, &distance_index));

    } else {
        //Otherwise, check the last thing in the zipcode to get the node values
        size_t node_depth = decoder.decoder_length()-1;

        ZipCode::code_type_t parent_type = decoder.get_code_type(node_depth-1);
        if (parent_type == ZipCode::IRREGULAR_SNARL || parent_type == ZipCode::CYCLIC_SNARL) {
            //If the parent is an irregular snarl
            return decoder.get_distance_index_address(node_depth-1);

        } else  {
            //TODO: I'm not sure about what to do about this, I don't like doing it here
            net_handle_t node_handle = distance_index.get_node_net_handle(id);
            net_handle_t parent = distance_index.get_parent(node_handle);
            if (distance_index.is_trivial_chain(parent)) {
                net_handle_t grandparent = distance_index.get_parent(parent);
                if (distance_index.is_simple_snarl(grandparent)) {
                    return distance_index.get_record_offset(distance_index.get_parent(grandparent));

                } else {
                    return distance_index.get_record_offset(grandparent);
                }
            } else {
                return distance_index.get_record_offset(parent);
            }
        }
    }
}

size_t MIPayload::node_record_offset(const ZipCode& zip, const SnarlDistanceIndex& distance_index, const nid_t& id) {

    //TODO: This is pointless but I'll keep it until I fix everything
    net_handle_t node_handle = distance_index.get_node_net_handle(id);
    return distance_index.get_node_record_offset(node_handle);
}

size_t MIPayload::node_length(const ZipCode& zip) {
    ZipCodeDecoder decoder (&zip);

    if (decoder.decoder_length() == 1) {
        //If the root-level structure is a node

        return decoder.get_length(0);

    } else if (decoder.decoder_length() == 2) {
        //If this is a node in the top-level chain

        return decoder.get_length(1);

   } else {
        //Otherwise, check the last thing in the zipcode to get the node values
        size_t node_depth = decoder.decoder_length()-1;

        return decoder.get_length(node_depth);
    }
}

bool MIPayload::is_reversed(const ZipCode& zip, const SnarlDistanceIndex& distance_index, const nid_t& id) {
 
    ZipCodeDecoder decoder (&zip);

    bool root_is_chain = decoder.get_code_type(0) != ZipCode::ROOT_SNARL;
    if (decoder.decoder_length() == 1) {
        //If the root-level structure is a node

        return false;

    } else if (decoder.decoder_length() == 2 && root_is_chain) {
        //If this is a node in the top-level chain

        return decoder.get_is_reversed_in_parent(1);
    
    } else if (decoder.decoder_length() == 2 && !root_is_chain) {
        //If the node is the child of the root snarl
        
        return false;
    } else {
        //Otherwise, check the last thing in the zipcode to get the node values
        size_t node_depth = decoder.decoder_length()-1;

        ZipCode:: code_type_t parent_type = decoder.get_code_type(node_depth-1);
        if (parent_type == ZipCode::IRREGULAR_SNARL || parent_type == ZipCode::CYCLIC_SNARL) {
            //If the parent is an irregular snarl
            return false;

        } else if (parent_type == ZipCode::REGULAR_SNARL) {
            //If the parent is a regular snarl

            //Because I'm storing "regular" and not "simple", need to check this
            if (distance_index.is_simple_snarl(distance_index.get_parent(distance_index.get_parent(distance_index.get_node_net_handle(id))))) {
                return decoder.get_is_reversed_in_parent(node_depth);
            } else {
                return false;
            }
        } else {
            //If the parent is a chain
            //If this was a node in a chain
            return decoder.get_is_reversed_in_parent(node_depth);
        }
    }
}

bool MIPayload::is_trivial_chain(const ZipCode& zip) {
 
    ZipCodeDecoder decoder (&zip);

    bool root_is_chain = decoder.get_code_type(0) != ZipCode::ROOT_SNARL;
    if (decoder.decoder_length() == 1) {
        //If the root-level structure is a node

        return true;
    } else if (decoder.decoder_length() == 2 && root_is_chain) {
        //If this is a node in the top-level chain

        return false;
    
    } else if (decoder.decoder_length() == 2 && !root_is_chain) {
        //If the node is the child of the root snarl
        
        return true;

    } else {
        //Otherwise, check the last thing in the zipcode to get the node values
        size_t node_depth = decoder.decoder_length()-1;

        ZipCode::code_type_t parent_type = decoder.get_code_type(node_depth-1);
        if (parent_type == ZipCode::IRREGULAR_SNARL || parent_type == ZipCode::CYCLIC_SNARL) {
            //If the parent is an irregular snarl
            return true;

        } else if (parent_type == ZipCode::REGULAR_SNARL) {
            //If the parent is a regular snarl
            return true;

        } else {
            //If the parent is a chain
            //If this was a node in a chain
            return false;
        }
    }
}

bool MIPayload::parent_is_chain(const ZipCode& zip, const SnarlDistanceIndex& distance_index, const nid_t& id) {
 
    ZipCodeDecoder decoder (&zip);

    bool root_is_chain = decoder.get_code_type(0) != ZipCode::ROOT_SNARL;
    if (decoder.decoder_length() == 1) {
        //If the root-level structure is a node

        return true;

    } else if (decoder.decoder_length() == 2 && root_is_chain) {
        //If this is a node in the top-level chain

        return true;
    
    } else if (decoder.decoder_length() == 2 && !root_is_chain) {
        //If the node is the child of the root snarl
        
        return false;

    } else {
        //Otherwise, check the last thing in the zipcode to get the node values
        size_t node_depth = decoder.decoder_length()-1;

        ZipCode::code_type_t parent_type = decoder.get_code_type(node_depth-1);
        if (parent_type == ZipCode::IRREGULAR_SNARL || parent_type == ZipCode::CYCLIC_SNARL) {
            //If the parent is an irregular snarl

            return false;

        } else if (parent_type == ZipCode::REGULAR_SNARL) {

            net_handle_t node_handle = distance_index.get_node_net_handle(id);
            net_handle_t parent = distance_index.get_parent(node_handle);
            if (distance_index.is_trivial_chain(parent)) {
                if (distance_index.is_simple_snarl(distance_index.get_parent(parent))) {
                    return true;
                } else {
                    return false;
                }
            } else {
                return true;
            }

        } else {
            //If the parent is a chain
            //If this was a node in a chain
            return true;

        }
    }
}


bool MIPayload::parent_is_root(const ZipCode& zip) {
 
    ZipCodeDecoder decoder (&zip);

    bool root_is_chain = decoder.get_code_type(0) != ZipCode::ROOT_SNARL;
    if (decoder.decoder_length() == 1) {
        //If the root-level structure is a node

        return true;

    } else if (decoder.decoder_length() == 2 && root_is_chain) {
        //If this is a node in the top-level chain

        return false;
    
    } else if (decoder.decoder_length() == 2 && !root_is_chain) {
        //If the node is the child of the root snarl
        
        return true;

    } else {

        return false;
    }
}


size_t MIPayload::prefix_sum(const ZipCode& zip, const SnarlDistanceIndex& distance_index, const nid_t& id) {
 
    ZipCodeDecoder decoder (&zip);

    bool root_is_chain = decoder.get_code_type(0) != ZipCode::ROOT_SNARL;
    if (decoder.decoder_length() == 1) {
        //If the root-level structure is a node
        return  0;

    } else if (decoder.decoder_length() == 2 && root_is_chain) {
        //If this is a node in the top-level chain

        return decoder.get_offset_in_chain(1);

    } else if (decoder.decoder_length() == 2 && !root_is_chain) {
        //If the node is the child of the root snarl
        return 0;
    } else {
        //Otherwise, check the last thing in the zipcode to get the node values
        size_t node_depth = decoder.decoder_length()-1;

        ZipCode::code_type_t parent_type = decoder.get_code_type(node_depth-1);
        if (parent_type == ZipCode::IRREGULAR_SNARL || parent_type == ZipCode::CYCLIC_SNARL) {
            return 0;
        } else if (parent_type == ZipCode::REGULAR_SNARL) {
            //If the parent is a snarl
            //Because I'm storing "regular" and not "simple", need to check this
            if (distance_index.is_simple_snarl(distance_index.get_parent(distance_index.get_parent(distance_index.get_node_net_handle(id))))) {
                return decoder.get_offset_in_chain(node_depth-1);
            } else {
                return 0;
            }
        } else {
            //If the parent is a chain
            //If this was a node in a chain
            return decoder.get_offset_in_chain(node_depth);
        }
    }
}

size_t MIPayload::chain_component(const ZipCode& zip, const SnarlDistanceIndex& distance_index, const nid_t& id) {
 
    ZipCodeDecoder decoder (&zip);

    bool root_is_chain = decoder.get_code_type(0) != ZipCode::ROOT_SNARL;

    if (decoder.decoder_length() == 1) {
        //If the root-level structure is a node

        return 0;

    } else if (decoder.decoder_length() == 2 && root_is_chain) {
        //If this is a node in the top-level chain

        net_handle_t net_handle = distance_index.get_node_net_handle(id);
        net_handle_t parent = distance_index.get_parent(net_handle);
        return distance_index.is_multicomponent_chain(parent) 
                ? distance_index.get_chain_component(net_handle)
                : 0;
    
    } else if (decoder.decoder_length() == 2 && !root_is_chain) {
        //If the node is the child of the root snarl
        
        return 0;
    } else {
        //Otherwise, check the last thing in the zipcode to get the node values
        size_t node_depth = decoder.decoder_length()-1;

        net_handle_t net_handle = distance_index.get_node_net_handle(id);
        net_handle_t parent = distance_index.get_parent(net_handle);
        return distance_index.is_multicomponent_chain(parent) 
                ? distance_index.get_chain_component(net_handle)
                : 0;
    }
}


}
