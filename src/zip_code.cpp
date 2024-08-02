#include "zip_code.hpp"

//#define DEBUG_ZIPCODE

namespace vg{
using namespace std;

void ZipCode::fill_in_zipcode (const SnarlDistanceIndex& distance_index, const pos_t& pos, bool fill_in_decoder) {

    std::vector<net_handle_t> ancestors;
    net_handle_t current_handle = distance_index.get_node_net_handle(id(pos));

    //Put all ancestors of the node in a vector, starting from the node, and not including the root
    while (!distance_index.is_root(current_handle)) {
        ancestors.emplace_back(current_handle);
        current_handle = distance_index.get_parent(current_handle);
    }

    //Make a temporary zipcode that will turn into the real one
    vector<size_t> temp_zipcode;
    temp_zipcode.reserve(ancestors.size() * 4);
    //Remember the maximum value we see to set the bitwidth when we make the real zipcode
    size_t max_value = 0;


    //Now add the root-level snarl or chain
    if (distance_index.is_root_snarl(current_handle)) {
        //FIrst thing is a snarl, so add the snarl's connected component number
        temp_zipcode.emplace_back(0);
#ifdef DEBUG_ZIPCODE
        cerr << "Adding code for top-level snarl " << distance_index.net_handle_as_string(current_handle) << endl;
#endif
        temp_zipcode.emplace_back(distance_index.get_connected_component_number(current_handle));
        max_value = std::max(max_value, temp_zipcode.back());
    } else {
#ifdef DEBUG_ZIPCODE
        cerr << "Adding code for top-level chain " << distance_index.net_handle_as_string(current_handle) << endl;
#endif
        //FIrst thing is a chain so add its connected component number and remove the chain from the stack
        temp_zipcode.emplace_back(1);
        max_value = std::max(max_value, temp_zipcode.back());

        //If the root-level structure is actually a chain, then save the connected component number and take out
        //the chain from the stack
        //If the root-level structure is a trivial chain, then just store the node (as a chain, which will have the 
        //connected-component number as the rank in the snarl anyways)
        temp_zipcode.emplace_back(distance_index.get_connected_component_number(ancestors.back()));
        max_value = std::max(max_value, temp_zipcode.back());
        if (ancestors.size() == 2 && distance_index.is_trivial_chain(ancestors.back())) {
#ifdef DEBUG_ZIPCODE
           cerr << "Adding code for top-level trivial chain" << endl;
#endif
            temp_zipcode.emplace_back(distance_index.minimum_length(ancestors.back())+1);
            max_value = std::max(max_value, temp_zipcode.back());
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
 
            temp_zipcode.emplace_back(connectivity);
            max_value = std::max(max_value, temp_zipcode.back());
            zipcode.from_vector(temp_zipcode, max_value);
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
            if (distance_index.is_looping_chain(ancestors.back())) {
                component += 1;
            }
            temp_zipcode.emplace_back(component);
            max_value = std::max(max_value, temp_zipcode.back());

        }
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
 
        temp_zipcode.emplace_back(connectivity);
        max_value = std::max(max_value, temp_zipcode.back());
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
            get_node_code(current_ancestor, distance_index, temp_zipcode, max_value);
        } else if (distance_index.is_chain(current_ancestor)) {
            get_chain_code(current_ancestor, distance_index, temp_zipcode, max_value);

            if (distance_index.is_trivial_chain(current_ancestor)) {
                zipcode.from_vector(temp_zipcode, max_value);
                if (fill_in_decoder) {
                    fill_in_full_decoder();
                }
                return;
            }
        } else if (distance_index.is_regular_snarl(current_ancestor)) {
            get_regular_snarl_code(current_ancestor, ancestors[i-1], distance_index, temp_zipcode, max_value);
        } else {
#ifdef DEBUG_ZIPCODE
            assert(distance_index.is_snarl(current_ancestor));
#endif
            get_irregular_snarl_code(current_ancestor, ancestors[i-1], distance_index, temp_zipcode, max_value);
        }
    }
    zipcode.from_vector(temp_zipcode, max_value);

    if (fill_in_decoder) {
        fill_in_full_decoder();
    }
}

void ZipCode::from_vector(const std::vector<size_t>& values, size_t max_value) {
    zipcode.from_vector(values, max_value);
}


void ZipCode::fill_in_full_decoder() {
    if (zipcode.size() == 0 || finished_decoding) {
        //If the zipcode is empty
        return;
    }
    decoder.reserve(zipcode.size() / 4);
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


    if (zip_length == 0) {
        //If there is nothing in the decoder yet, then the first thing will start at 0

        //Is the root a chain/node?
        decoder.emplace_back(zipcode.at(ROOT_IS_CHAIN_OFFSET), 0);

#ifdef DEBUG_ZIPCODE
cerr << "\tadding the root, which is a " << (decoder.back().first ? "chain or node" : "snarl") << endl;
#endif
        if (zipcode.size() == ROOT_NODE_SIZE) {
            //If this was a root node, then we're done
            finished_decoding = true;
            return true;
        } else {
            //There might be something else but we're done for now
            return false;
        }
    } else {
        //This is not a root
        bool previous_is_chain = decoder.back().first;
        size_t previous_start = decoder.back().second;

        if (previous_is_chain) {
            //If the last thing was chain, then either the chain was the last thing in the zipcode
            // (if it was the child of a snarl) or the next thing is either a node or snarl

            assert(std::min(ZipCode::CHAIN_SIZE + ZipCode::REGULAR_SNARL_SIZE,
                            ZipCode::CHAIN_SIZE + ZipCode::IRREGULAR_SNARL_SIZE) > ZipCode::NODE_SIZE);

            size_t this_size = zip_length == 1 ? ROOT_CHAIN_SIZE : CHAIN_SIZE;
            if (zipcode.size() == previous_start + this_size) {
                //If the zipcode ends here
#ifdef DEBUG_ZIPCODE
            cerr << "The last thing was a trivial chain so we're done" << endl;
#endif
                finished_decoding = true;
                return true;
            } else if (zipcode.size() == previous_start + this_size + NODE_SIZE) {
                //If the zipcode ends after the node, add the node and we're done
#ifdef DEBUG_ZIPCODE
            cerr << "Adding a node and we're done" << endl;
#endif
                decoder.emplace_back(true, previous_start + this_size);
                finished_decoding = true;
                return true;
            } else {
                //Otherwise, this is a snarl and we're not done
#ifdef DEBUG_ZIPCODE
            cerr << "Adding a snarl starting at " << (previous_start + this_size) << endl;
#endif
                decoder.emplace_back(false, previous_start + this_size);
                return false;
            }
        } else {
            //Otherwise, the last thing was a snarl
            size_t next_start = previous_start;

            //The regular/irregular snarl tag
            if (zip_length == 1) {
                //IF this was a root snarl
                next_start += ROOT_SNARL_SIZE;
            } else if (zipcode.at(previous_start + SNARL_IS_REGULAR_OFFSET) == 1) {
                //If this was a regular snarl
                next_start += REGULAR_SNARL_SIZE;
            } else {
                //Technically it could be irregular or cyclic but it doesn't matter because the codes are the same
                next_start += IRREGULAR_SNARL_SIZE;            
            }
            decoder.emplace_back(true, next_start);
            return false;
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
            size_t code_type_int = zipcode.at(decoder[depth].second + ZipCode::SNARL_IS_REGULAR_OFFSET);
            if (code_type_int == 0) {
                return IRREGULAR_SNARL;
            } else if (code_type_int == 1) {
                return REGULAR_SNARL;
            } else {
                return CYCLIC_SNARL;
            }
        }
    }
}

size_t ZipCode::get_length(const size_t& depth, const SnarlDistanceIndex* distance_index) const {

    if (depth == 0) {
        //If this is the root chain/snarl/node

        if (decoder_length() == 1) {
            //If the length is 1, then it's a node
            size_t zip_value = zipcode.at(decoder[depth].second + ROOT_NODE_LENGTH_OFFSET);
            return zip_value == 0 ? std::numeric_limits<size_t>::max() : zip_value-1;

        } else {
            //Otherwise, we didn't store the length
            throw std::runtime_error("zipcodes don't store lengths of top-level chains or snarls");
        }
    } else if (decoder[depth].first) {
        //If this is a chain/node

        //If this is a chain or a node, then the length will be the second thing
        assert(CHAIN_LENGTH_OFFSET == NODE_LENGTH_OFFSET);
        size_t zip_value = zipcode.at(decoder[depth].second + CHAIN_LENGTH_OFFSET);
        return zip_value == 0 ? std::numeric_limits<size_t>::max() : zip_value-1;
    } else {
        //If this is a snarl

        size_t zip_value = zipcode.at(decoder[depth].second + SNARL_LENGTH_OFFSET);
        return zip_value == 0 ? std::numeric_limits<size_t>::max() : zip_value-1;
    }
}

size_t ZipCode::get_rank_in_snarl(const size_t& depth) const {


    if (depth == 0) {
        //If this is the root chain/snarl/node
        throw std::runtime_error("zipcodes don't store ranks of top-level chains or snarls");

    } else if (decoder[depth].first) {
        //If this is a chain/node

        if (decoder[depth-1].first) {
            throw std::runtime_error("zipcodes trying to find the rank in snarl of a node in a chain");
        }

        return zipcode.at(decoder[depth].second + CHAIN_RANK_IN_SNARL_OFFSET);
    } else {
        //If this is a snarl
        throw std::runtime_error("zipcodes don't store snarl ranks for snarls");
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

    } else if (!decoder[depth].first) {
        //If this is a snarl

        return zipcode.at(decoder[depth].second + SNARL_CHILD_COUNT_OFFSET);
    } else {
        //If this is not a snarl
        throw std::runtime_error("trying to get the snarl child count of a non-snarl zipcode");
    }
}

size_t ZipCode::get_offset_in_chain(const size_t& depth, const SnarlDistanceIndex* distance_index) const {


    if (depth == 0) {
        //If this is the root chain/snarl/node
        throw std::runtime_error("zipcodes don't have chain offsets for roots");

    } else if (decoder[depth].first) {
        //If this is a chain/node

        if (!decoder[depth-1].first) {
            throw std::runtime_error("zipcodes trying to find the offset in child of a snarl");
        }
        size_t zip_value = zipcode.at(decoder[depth].second + NODE_OFFSET_OFFSET);

        return zip_value == 0 ? std::numeric_limits<size_t>::max() : zip_value-1;
    } else {
        //If this is a snarl

        size_t zip_value = zipcode.at(decoder[depth].second + SNARL_OFFSET_IN_CHAIN_OFFSET);

        return zip_value == 0 ? std::numeric_limits<size_t>::max() : zip_value-1;
    }
}
size_t ZipCode::get_chain_component(const size_t& depth) const {


    if (depth == 0) {
        //If this is the root chain/snarl/node
        throw std::runtime_error("zipcodes don't have chain offsets for roots");

    } else if (decoder[depth].first) {
        //If this is a chain/node

        if (!decoder[depth-1].first) {
            throw std::runtime_error("zipcodes trying to find the offset in child of a snarl");
        }
        return zipcode.at(decoder[depth].second + NODE_CHAIN_COMPONENT_OFFSET);
    } else {
        //If this is a snarl

        return zipcode.at(decoder[depth].second + SNARL_CHAIN_COMPONENT_OFFSET);
    }
}

size_t ZipCode::get_last_chain_component(const size_t& depth, bool get_end) const {

    if (!decoder[depth].first) {
        throw std::runtime_error("zipcodes trying to find the last chain component a snarl");
    }
    size_t zip_value = zipcode.at(decoder[depth].second + CHAIN_COMPONENT_COUNT_OFFSET);
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

    if (!decoder[depth].first) {
        throw std::runtime_error("zipcodes trying to find the last chain component a snarl");
    }
    return zipcode.at(decoder[depth].second + CHAIN_COMPONENT_COUNT_OFFSET) % 2;
}
bool ZipCode::get_is_reversed_in_parent(const size_t& depth) const {


    if (depth == 0) {
        //If this is the root chain/snarl/node
        return false;

    } else if (decoder[depth].first) {
        //If this is a chain/node

        if (decoder[depth-1].first) {
            //If the parent is a chain, then this is a node and we need to check its orientation

            return zipcode.at(decoder[depth].second + NODE_IS_REVERSED_OFFSET);
        } else {
            //If the parent is a snarl, then this might be a chain in a regular snarl

            size_t snarl_type = zipcode.at(decoder[depth-1].second + SNARL_IS_REGULAR_OFFSET);
            if (snarl_type == 1) {
                //The parent is a regular snarl, which stores is_reversed for the child
                
                return zipcode.at(decoder[depth-1].second + REGULAR_SNARL_IS_REVERSED_OFFSET);
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

        return distance_index->get_handle_from_connected_component(zipcode.at(ROOT_IDENTIFIER_OFFSET));

    } else if (decoder[depth].first) {
        //If this is a chain/node

        throw std::runtime_error("zipcodes trying to get a handle of a chain or node");
    } else {
        //If this is a snarl

        size_t snarl_type = zipcode.at(decoder[depth].second + SNARL_IS_REGULAR_OFFSET);
        if (snarl_type == 1) {
            //If this is a regular snarl

            throw std::runtime_error("zipcodes trying to get a handle of a regular snarl");
        } else {
            //Irregular snarl

            size_t zip_value = zipcode.at(decoder[depth].second + IRREGULAR_SNARL_RECORD_OFFSET);
            net_handle_t snarl_handle = distance_index->get_net_handle_from_values(zip_value, 
                                                                                   SnarlDistanceIndex::START_END, 
                                                                                   SnarlDistanceIndex::SNARL_HANDLE);
            return snarl_handle;
        }
    }
}

net_handle_t ZipCode::get_net_handle_slow(nid_t id, const size_t& depth, const SnarlDistanceIndex* distance_index) const {
    //This is just copying get_net_handle except adding a slower version for the things we don't remember

    if (depth == 0) {
        //If this is the root chain/snarl/node

        return distance_index->get_handle_from_connected_component(zipcode.at(ROOT_IDENTIFIER_OFFSET));

    } else if (decoder[depth].first) {
        //If this is a chain/node

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

        size_t snarl_type = zipcode.at(decoder[depth].second + SNARL_IS_REGULAR_OFFSET);
        if (snarl_type == 1) {
            //If this is a regular snarl

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

            size_t zip_value = zipcode.at(decoder[depth].second + IRREGULAR_SNARL_RECORD_OFFSET);
            net_handle_t snarl_handle = distance_index->get_net_handle_from_values(zip_value, 
                                                                                   SnarlDistanceIndex::START_END, 
                                                                                   SnarlDistanceIndex::SNARL_HANDLE);
            return snarl_handle;
        }
    }
}


size_t ZipCode::get_distance_index_address(const size_t& depth) const {


    if (depth == 0) {
        //If this is the root chain/snarl/node

        return zipcode.at(ROOT_IDENTIFIER_OFFSET);

    } else if (decoder[depth].first) {
        //If this is a chain/node

        throw std::runtime_error("zipcodes trying to get a handle of a chain or node");
    } else {
        //If this is a snarl

        size_t snarl_type = zipcode.at(decoder[depth].second + SNARL_IS_REGULAR_OFFSET);
        if (snarl_type == 1) {
            //If this is a regular snarl

            throw std::runtime_error("zipcodes trying to get a handle of a regular ansl");
        } else {
            //Irregular snarl

            return zipcode.at(decoder[depth].second + IRREGULAR_SNARL_RECORD_OFFSET);
        }
    }
}
size_t ZipCode::get_distance_to_snarl_bound(const size_t& depth, bool snarl_start, bool left_side) const {

#ifdef DEBUG_ZIPCODE
    assert(depth > 0);
    assert((get_code_type(depth-1) == ZipCode::IRREGULAR_SNARL || get_code_type(depth-1) == ZipCode::REGULAR_SNARL || get_code_type(depth-1) == ZipCode::CYCLIC_SNARL)); 
#endif
     size_t snarl_type = zipcode.at(decoder[depth-1].second + SNARL_IS_REGULAR_OFFSET);
     if (snarl_type == 1) {
         //The parent is a regular snarl, which stores is_reversed for the child

         size_t zip_value = zipcode.at(decoder[depth-1].second + REGULAR_SNARL_IS_REVERSED_OFFSET);
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
        size_t zip_value = zipcode.at(decoder[depth-1].second + distance_offset);
        return zip_value == 0 ? std::numeric_limits<size_t>::max() : zip_value - 1;
     }
}

bool ZipCode::is_externally_start_end_connected (const size_t& depth) const {
    assert(depth == 0);
    assert(decoder[0].first);
    size_t zip_value = zipcode.at(decoder[depth].second + ROOT_NODE_OR_CHAIN_CONNECTIVITY_OFFSET);
    return (zip_value & 1) != 0;
}
bool ZipCode::is_externally_start_start_connected (const size_t& depth) const {
    assert(depth == 0);
    assert(decoder[0].first);
    size_t zip_value = zipcode.at(decoder[depth].second + ROOT_NODE_OR_CHAIN_CONNECTIVITY_OFFSET);
    return (zip_value & 2) != 0;
}
bool ZipCode::is_externally_end_end_connected (const size_t& depth) const {
    assert(depth == 0);
    assert(decoder[0].first);
    size_t zip_value = zipcode.at(decoder[depth].second + ROOT_NODE_OR_CHAIN_CONNECTIVITY_OFFSET);
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
    // Print out the numbers in a way that is easy to copy-paste as a vector literal.
    out << "<zipcode {";
    for (size_t i = 0; i < zipcode.size(); i++) {
        out << zipcode.at(i);
        if (i + 1 < zipcode.size()) {
            out << ", ";
        }
    }
    out << "}>";
}

std::ostream& operator<<(std::ostream& out, const ZipCode& zip) {
    return out << "<zipcode {" << zip.get_identifier(zip.max_depth())<< "}>";
}


void ZipCode::get_node_code(const net_handle_t& node, const SnarlDistanceIndex& distance_index,
                            vector<size_t>& temp_zipcode, size_t& max_value) {
#ifdef DEBUG_ZIPCODE
    assert(!distance_index.is_trivial_chain(node));
    assert((distance_index.is_chain(distance_index.get_parent(node)) || distance_index.is_root(distance_index.get_parent(node))));
#endif
    size_t start_i = temp_zipcode.size();
    temp_zipcode.resize(start_i + NODE_SIZE);
    //Node code is: offset in chain, length, is reversed, chain component

    //Assume this node is in a regular chain
    size_t prefix_sum = distance_index.get_prefix_sum_value(node); 
    temp_zipcode[start_i + NODE_OFFSET_OFFSET] = prefix_sum == std::numeric_limits<size_t>::max() ? 0 : prefix_sum+1;
    max_value = std::max(max_value, temp_zipcode[start_i + NODE_OFFSET_OFFSET]);

    temp_zipcode[start_i + NODE_LENGTH_OFFSET] = distance_index.minimum_length(node)+1;
    max_value = std::max(max_value, temp_zipcode[start_i + NODE_LENGTH_OFFSET]);

    temp_zipcode[start_i + NODE_IS_REVERSED_OFFSET] = distance_index.is_reversed_in_parent(node);
    max_value = std::max(max_value, temp_zipcode[start_i + NODE_IS_REVERSED_OFFSET]);

    size_t component = distance_index.get_chain_component(node);
    temp_zipcode[start_i + NODE_CHAIN_COMPONENT_OFFSET] = component == std::numeric_limits<size_t>::max() ? 0 : component;
    max_value = std::max(max_value, temp_zipcode[start_i + NODE_CHAIN_COMPONENT_OFFSET]);

    return;

}
void ZipCode::get_chain_code(const net_handle_t& chain, const SnarlDistanceIndex& distance_index,
                            vector<size_t>& temp_zipcode, size_t& max_value) {
    //Chain code is: rank in snarl, length

    size_t start_i = temp_zipcode.size();
    temp_zipcode.resize(start_i + CHAIN_SIZE);

    //Rank in snarl
    temp_zipcode[start_i + CHAIN_RANK_IN_SNARL_OFFSET] = distance_index.get_rank_in_parent(chain);
    max_value = std::max(max_value, temp_zipcode[start_i + CHAIN_RANK_IN_SNARL_OFFSET]);

    //Length
    size_t len = distance_index.minimum_length(chain);
    temp_zipcode[start_i + CHAIN_LENGTH_OFFSET] = len == std::numeric_limits<size_t>::max() ? 0 : len+1;
    max_value = std::max(max_value, temp_zipcode[start_i + CHAIN_LENGTH_OFFSET]);

    //Component count and if it loops
    bool is_trivial = distance_index.is_trivial_chain(chain) ;
    size_t component = is_trivial
                       ? 0 
                       : distance_index.get_chain_component(distance_index.get_bound(chain, true, false), true);
    component = component == std::numeric_limits<size_t>::max() ? 0 : component*2;
    if (!is_trivial && distance_index.is_looping_chain(chain)) {
        component += 1;
    }
    temp_zipcode[start_i + CHAIN_COMPONENT_COUNT_OFFSET] = component;
    max_value = std::max(max_value, component);

    return;

}
void ZipCode::get_regular_snarl_code(const net_handle_t& snarl, const net_handle_t& snarl_child, 
                            const SnarlDistanceIndex& distance_index,
                            vector<size_t>& temp_zipcode, size_t& max_value) {

    size_t start_i = temp_zipcode.size();
    temp_zipcode.resize(start_i + REGULAR_SNARL_SIZE);


    //Tag to say that it's a regular snarl
    temp_zipcode[start_i + SNARL_IS_REGULAR_OFFSET] = 1;

    //The number of children
    size_t child_count = 0;
    distance_index.for_each_child(snarl, [&] (const net_handle_t& child) {
        child_count++;
    });
    temp_zipcode[start_i + SNARL_CHILD_COUNT_OFFSET] = child_count;
    max_value = std::max(max_value, child_count);

    //Chain prefix sum value for the start of the snarl, which is the prefix sum of the start node + length of the start node
    net_handle_t start_node = distance_index.get_node_from_sentinel(distance_index.get_bound(snarl, false, false));
    size_t prefix_sum = SnarlDistanceIndex::sum(distance_index.get_prefix_sum_value(start_node), distance_index.minimum_length(start_node));
    temp_zipcode[start_i + SNARL_OFFSET_IN_CHAIN_OFFSET] = (prefix_sum == std::numeric_limits<size_t>::max() ? 0 : prefix_sum+1);
    max_value = std::max(max_value, temp_zipcode[start_i + SNARL_OFFSET_IN_CHAIN_OFFSET]);

    size_t component = distance_index.get_chain_component(start_node);
    temp_zipcode[start_i + SNARL_CHAIN_COMPONENT_OFFSET] = component == std::numeric_limits<size_t>::max() ? 0 : component;
    max_value = std::max(max_value, temp_zipcode[start_i + SNARL_CHAIN_COMPONENT_OFFSET]);

    //Length of the snarl
    size_t len = distance_index.minimum_length(snarl);
    temp_zipcode[start_i + SNARL_LENGTH_OFFSET] = (len == std::numeric_limits<size_t>::max() ? 0 : len+1);
    max_value = std::max(max_value, temp_zipcode[start_i + SNARL_LENGTH_OFFSET]);

    //Is the child of the snarl reversed in the snarl
#ifdef DEBUG_ZIPCODE
    assert(distance_index.is_chain(snarl_child));
#endif
    temp_zipcode[start_i + REGULAR_SNARL_IS_REVERSED_OFFSET] = (distance_index.distance_in_parent(snarl, 
                                                        distance_index.get_bound(snarl, false, true),
                                                        distance_index.flip(distance_index.canonical(snarl_child))) != 0);
    max_value = std::max(max_value, temp_zipcode[start_i + REGULAR_SNARL_IS_REVERSED_OFFSET]);

    return;

}
void ZipCode::get_irregular_snarl_code(const net_handle_t& snarl, const net_handle_t& snarl_child, 
                                                 const SnarlDistanceIndex& distance_index,
                            vector<size_t>& temp_zipcode, size_t& max_value) {

    size_t start_i = temp_zipcode.size();
    temp_zipcode.resize(start_i + IRREGULAR_SNARL_SIZE);

    //Tag to say that it's an irregular snarl
    temp_zipcode[start_i + SNARL_IS_REGULAR_OFFSET] = distance_index.is_dag(snarl) ? 0 : 2;
    max_value = std::max(max_value, temp_zipcode[start_i + SNARL_IS_REGULAR_OFFSET]);

    //The number of children
    size_t child_count = 0;
    distance_index.for_each_child(snarl, [&] (const net_handle_t& child) {
        child_count++;
    });
    temp_zipcode[start_i + SNARL_CHILD_COUNT_OFFSET] = child_count;
    max_value = std::max(max_value, child_count);

    //Chain prefix sum value for the start of the snarl, which is the prefix sum of the start node + length of the start node
    net_handle_t start_node = distance_index.get_node_from_sentinel(distance_index.get_bound(snarl, false, false));
    size_t prefix_sum = SnarlDistanceIndex::sum(distance_index.get_prefix_sum_value(start_node), distance_index.minimum_length(start_node));
    temp_zipcode[start_i + SNARL_OFFSET_IN_CHAIN_OFFSET] = (prefix_sum == std::numeric_limits<size_t>::max() ? 0 : prefix_sum+1);
    max_value = std::max(max_value, temp_zipcode[start_i + SNARL_OFFSET_IN_CHAIN_OFFSET]);

    size_t component = distance_index.get_chain_component(start_node);
    temp_zipcode[start_i + SNARL_CHAIN_COMPONENT_OFFSET] = component == std::numeric_limits<size_t>::max() ? 0 : component;
    max_value = std::max(max_value, temp_zipcode[start_i + SNARL_CHAIN_COMPONENT_OFFSET]);

    //Length of the snarl
    size_t len = distance_index.minimum_length(snarl);
    temp_zipcode[start_i + SNARL_LENGTH_OFFSET] = (len == std::numeric_limits<size_t>::max() ? 0 : len+1);
    max_value = std::max(max_value, temp_zipcode[start_i + SNARL_LENGTH_OFFSET]);


    //Record offset to look up distances in the index later
    temp_zipcode[start_i + IRREGULAR_SNARL_RECORD_OFFSET] = (distance_index.get_record_offset(snarl));
    max_value = std::max(max_value, temp_zipcode[start_i + IRREGULAR_SNARL_RECORD_OFFSET]);

    temp_zipcode[start_i + IRREGULAR_SNARL_DISTANCE_LEFT_START_OFFSET] = distance_index.distance_to_parent_bound(snarl, true, distance_index.flip(snarl_child));
    temp_zipcode[start_i + IRREGULAR_SNARL_DISTANCE_LEFT_END_OFFSET] = distance_index.distance_to_parent_bound(snarl, false, distance_index.flip(snarl_child));
    temp_zipcode[start_i + IRREGULAR_SNARL_DISTANCE_RIGHT_START_OFFSET] = distance_index.distance_to_parent_bound(snarl, true, snarl_child);
    temp_zipcode[start_i + IRREGULAR_SNARL_DISTANCE_RIGHT_END_OFFSET] = distance_index.distance_to_parent_bound(snarl, false, snarl_child);


    //Add 1 to values to store inf properly
    temp_zipcode[start_i + IRREGULAR_SNARL_DISTANCE_LEFT_START_OFFSET] = 
        temp_zipcode[start_i + IRREGULAR_SNARL_DISTANCE_LEFT_START_OFFSET] == std::numeric_limits<size_t>::max() 
        ? 0 
        : temp_zipcode[start_i + IRREGULAR_SNARL_DISTANCE_LEFT_START_OFFSET] + 1;
    temp_zipcode[start_i + IRREGULAR_SNARL_DISTANCE_RIGHT_START_OFFSET] = 
        temp_zipcode[start_i + IRREGULAR_SNARL_DISTANCE_RIGHT_START_OFFSET] == std::numeric_limits<size_t>::max() 
        ? 0 
        : temp_zipcode[start_i + IRREGULAR_SNARL_DISTANCE_RIGHT_START_OFFSET] + 1;
    temp_zipcode[start_i + IRREGULAR_SNARL_DISTANCE_LEFT_END_OFFSET] = 
        temp_zipcode[start_i + IRREGULAR_SNARL_DISTANCE_LEFT_END_OFFSET] == std::numeric_limits<size_t>::max() 
        ? 0 
        : temp_zipcode[start_i + IRREGULAR_SNARL_DISTANCE_LEFT_END_OFFSET] + 1;
    temp_zipcode[start_i + IRREGULAR_SNARL_DISTANCE_RIGHT_END_OFFSET] = 
        temp_zipcode[start_i + IRREGULAR_SNARL_DISTANCE_RIGHT_END_OFFSET] == std::numeric_limits<size_t>::max() 
        ? 0 
        : temp_zipcode[start_i + IRREGULAR_SNARL_DISTANCE_RIGHT_END_OFFSET] + 1;

    max_value = std::max(max_value, temp_zipcode[start_i + IRREGULAR_SNARL_DISTANCE_LEFT_START_OFFSET]);
    max_value = std::max(max_value, temp_zipcode[start_i + IRREGULAR_SNARL_DISTANCE_RIGHT_START_OFFSET]);
    max_value = std::max(max_value, temp_zipcode[start_i + IRREGULAR_SNARL_DISTANCE_LEFT_END_OFFSET]);
    max_value = std::max(max_value, temp_zipcode[start_i + IRREGULAR_SNARL_DISTANCE_RIGHT_END_OFFSET]);

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
            //If its a regular snarl, then the distances to the ends are either 0 or inf
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
        //FInd the distances to the ends of the parent
        update_distances_to_ends_of_parent(zip1, i+1, distance_to_start1, distance_to_end1);
    }
#ifdef DEBUG_ZIPCODE
cerr << "Finding distances to ancestors of second position" << endl;
#endif
    //The same thing for the second position
    for (int i = zip2.decoder_length()-2 ; i > 0 && i > lowest_common_ancestor_depth ; i--) {
        //the parent snarl tree node is at index i
        //The distances are currently to the ends of the current node
        //FInd the distances to the ends of the parent

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
        } else if ( zip1.decoder[depth].first) {
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

    if (zip1.decoder[0].first != zip2.decoder[0].first) {
#ifdef DEBUG_ZIPCODE
        cerr << "Zip codes are on different connected components" << endl;
#endif
        return true;
    }

    if (zip1.get_distance_index_address(0) != zip2.get_distance_index_address(0)) {
#ifdef DEBUG_ZIPCODE
        cerr << "Zip codes are on different connected components" << endl;
#endif
        return true;
    }

    //The depth of a chain that both zips are on 
    size_t shared_depth = 0;

    if (!zip1.decoder[0].first) {
        //If the top-level thing is a snarl, then check if the zips are in the same chain. 
        //If they are, then proceed from the shared chain

        if (zip1.get_rank_in_snarl(1) != zip2.get_rank_in_snarl(1)) {
            //We can't tell
            return false;
        }
        //Next check the length of the chain
        if (zip1.get_length(1) < limit) {
            return true;
        }
        //The two zipcodes are on the same chain at depth 1
        shared_depth = 1;

        //The zips now point to the children of the shared chain, so we can proceed as if the top-level
        //structure was a chain

    }

    //Both zips now point to a thing in a shared chain
    //Get the minimum possible distance between the structures on the chain
    //For a lower bound, this assumes that the positions are as close as they can be on the structure in the chain
    size_t prefix_sum1 = zip1.get_offset_in_chain(shared_depth+1);
    size_t prefix_sum2 = zip2.get_offset_in_chain(shared_depth+1);
    size_t length1 = zip1.get_length(shared_depth+1); 
    size_t length2 = zip2.get_length(shared_depth+1); 
    size_t component1 = zip1.get_chain_component(shared_depth+1); 
    size_t component2 = zip2.get_chain_component(shared_depth+1); 

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
#ifdef DEBUG_ZIPCODE
    cerr << "Encode integers: ";
    for (size_t i = 0 ; i < zipcode.size() ; i++) {
        cerr << zipcode.at(i) << " ";
    }
    cerr << endl;
#endif
    if (bit_count() > 112) {
        //If there aren't enough bits to represent the zip code
        return MIPayload::NO_CODE;
    }
    //The values that get returned
    code_type encoded1 = 0;
    code_type encoded2 = 0;

    //The first (leftmost of first int) 8 bits is the width
    encoded1 |= zipcode.get_bit_width();

    //Left shift by 8 to make space for the next thing we're adding
    encoded1 <<= 8;
    //The second 8 bits is the number of items in the vector (not the number of bits)
    encoded1 |= zipcode.size(); 
    encoded1 <<= 1;

#ifdef DEBUG_ZIPCODE
cerr << "Encode the bit width "<< ((size_t) zipcode.get_bit_width()) << " and size " << zipcode.size() << endl;
cerr << "\t";
#endif
    

    //16 bits are set, so 112 left
    //Now add each bit one by one and left shift to make space for the next one
    for (size_t i = 0 ; i < 112 ; i++ ) {
        if ( i < 48 ) {
            //Add to first code, just one bit to the end
            if (i < zipcode.get_bit_count() && zipcode.bit_at(i)) {
                encoded1 |= 1;
#ifdef DEBUG_ZIPCODE
                cerr << "1";
#endif
            }
#ifdef DEBUG_ZIPCODE
            else {
                cerr << "0";
            }
#endif
            //Left shift by one after everything except the last bit
            if (i != 47) {
                encoded1 <<= 1;
            }
        } else {
            //Add to second code
            if (i < zipcode.get_bit_count() && zipcode.bit_at(i)) {
                encoded2 |= 1;
#ifdef DEBUG_ZIPCODE
                cerr << "1";
#endif
            }
#ifdef DEBUG_ZIPCODE
            else {
                cerr << "0";
            }
#endif
            if ( i != 111) {
                encoded2 <<= 1;
            }
        }
    
    }
#ifdef DEBUG_ZIPCODE
    cerr << endl;
    cerr << "Actual ints being stored: " << encoded1 << " and " << encoded2 << ": ";
    for (int i = 63 ; i >= 0 ; --i) {
        if (((size_t) 1 << i) & encoded1) {
            cerr << "1";
        } else {
            cerr << "0";
        }
    }
    for (int i = 63 ; i >= 0 ; --i) {
        if (((size_t) 1 << i) & encoded2) {
            cerr << "1";
        } else {
            cerr << "0";
        }
    }
    cerr << endl;
#endif
    return {encoded1, encoded2};

}

void ZipCode::fill_in_zipcode_from_payload(const gbwtgraph::Payload& payload) {
    assert(payload != MIPayload::NO_CODE);

    //First 8 bits of first int is the width
    size_t width = payload.first >> 56;
    zipcode.set_bit_width((uint8_t)width);

    //Second 8 bits is the item count
    size_t item_count = (payload.first >> 48) & ((1 << 8)-1);

    //bit count is the product of the two
    size_t bit_count = (size_t)width * (size_t)item_count;
    zipcode.set_bitvector_length(bit_count);

#ifdef DEBUG_ZIPCODE
    cerr << "Get zipcode from payload " << payload.first << " and " << payload.second<< " with  width: " << width << " item count " << item_count << " meaning " << bit_count << " bits" << endl;
    cerr << "\t";
#endif


    //Mask for checking the relevant bit
    //Start by checking the 17th bit from the left
    //Right shift by one for each bit we look at
    uint64_t mask1 = (uint64_t)1 << 47;
    uint64_t mask2 = (uint64_t)1 << 63;
    //get one bit at a time from the payload and add it to the zip code
    for (size_t i = 0 ; i < bit_count ; i++) {
        if (i < 48) {
            if ((payload.first & mask1) != 0) {
                zipcode.set_bit_at(i);
#ifdef DEBUG_ZIPCODE
                cerr << "1";
#endif
            }
#ifdef DEBUG_ZIPCODE
            else {
                cerr << "0";
            }
#endif
            mask1 >>= 1;
        } else {
            if ((payload.second & mask2) != 0) {
                zipcode.set_bit_at(i);
#ifdef DEBUG_ZIPCODE
                cerr << "1";
#endif
            }
#ifdef DEBUG_ZIPCODE
            else {
                cerr << "0";
            }
#endif
            mask2 >>= 1;
        }
    }
#ifdef DEBUG_ZIPCODE
    cerr << endl;
    cerr << "Found encoded integers: ";
    for (size_t i = 0 ; i < zipcode.size() ; i++) {
        cerr << zipcode.at(i) << " ";
    }
    cerr << endl;
#endif
    return;
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
    //The zipcode vector will be serialized as a bunch of min_width_int_vector_ts
    //The first min_width_int_vector_t will have one value, which will be the length of the
    //zipcode that follows it
#ifdef DEBUG_ZIPCODE
    cerr << "Serialize zipcode collection" << endl;
#endif

    //First serialize the header, which is the magic number and version
    uint32_t magic = magic_number;
    uint32_t vers = version;
    out.write(reinterpret_cast<const char*>(&magic), sizeof(magic));
    out.write(reinterpret_cast<const char*>(&vers), sizeof(vers));


    for (const ZipCode& zip : zipcodes) {

        //Write the width
        uint8_t width = (uint8_t) zip.zipcode.get_bit_width();
        out.write(reinterpret_cast<const char*>(&width), sizeof(width));
    
        //How many values are in the vector. Used with width to get the bit count
        size_t item_count = zip.zipcode.size();

        out.write(reinterpret_cast<const char*>(&item_count), sizeof(item_count));


        //Write the zipcode 
#ifdef DEBUG_ZIPCODE
        cerr << "Write width " << (size_t) width << " and item count " << item_count << " and zipcode: " << endl;
        cerr << "\t";
        for (size_t i = 0 ; i < zip.zipcode.size() ; i++) {
            cerr << zip.zipcode.at(i) << " ";
        }
        cerr << endl << "\t";
        size_t zip_byte_count = 0;
#endif
        size_t bit_count = zip.zipcode.get_bit_count();
        for (size_t i = 0 ; i < bit_count ; i += 8) {
#ifdef DEBUG_ZIPCODE
            zip_byte_count++;
#endif
            uint8_t result = 0;
            for (size_t j = 0 ; j < 8 ; j++) {
                result <<= 1;
                if (i+j < bit_count && zip.zipcode.bit_at(i+j)) {
#ifdef DEBUG_ZIPCODE
                    cerr << "1";
#endif
                    result |= 1;
                }
#ifdef DEBUG_ZIPCODE
                else {
                    cerr << "0";
                }
#endif
            }
            out << char(result);
        }
#ifdef DEBUG_ZIPCODE
        cerr << endl;
        assert(zip_byte_count == ceil((float)bit_count / 8)); 
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

        //First, get the bitwidth of the vector
        uint8_t width;
        in.read(reinterpret_cast<char*>(&width), sizeof(width));

        //Next, get the number of items in the zipcode
        size_t item_count;
        in.read(reinterpret_cast<char*>(&item_count), sizeof(item_count)); 

        size_t bit_count = (size_t)width * item_count;

        //How many bytes were used to store all the bits in the zipcode bit vector
        size_t byte_count = (size_t) std::ceil((float)bit_count / 8);



#ifdef DEBUG_ZIPCODE
        cerr << "Get zipcode with width " << (size_t) width << " and item count " << item_count << endl << "\t";
#endif

        char line [byte_count];

        in.read(line, byte_count);

        ZipCode zip;
        zip.zipcode.set_bit_width(width);
        zip.zipcode.set_bitvector_length(bit_count);
        size_t added_bits = 0;
        for (const char& character : line) {
            for (int i = 7 ; i >= 0 ; i--) {
                if (added_bits < bit_count) {
                    if  (((uint8_t)character & ((uint8_t)1 << i)) != 0)  {
                        zip.zipcode.set_bit_at(added_bits);
#ifdef DEBUG_ZIPCODE
                        cerr << "1";
#endif
                    }
#ifdef DEBUG_ZIPCODE
                    else {
                        cerr << "0";
                    }
#endif
                    added_bits++;
                }
            }
        }
#ifdef DEBUG_ZIPCODE
        cerr << endl <<"\t";
        for (size_t i = 0 ; i < zip.zipcode.size() ; i++) {
            cerr << zip.zipcode.at(i) << " ";
        }
        cerr << endl;
#endif

        zipcodes.emplace_back(std::move(zip));
    }

}
MIPayload ZipCode::get_payload_from_zipcode(nid_t id, const SnarlDistanceIndex& distance_index) const {
    MIPayload payload;

    if (decoder_length() == 1) {
        //If the root-level structure is a node
        payload.parent_is_root = true;
        payload.parent_is_chain = true;

        payload.node_handle = distance_index.get_net_handle_from_values(
                distance_index.get_record_offset(distance_index.get_handle_from_connected_component(get_distance_index_address(0))),
                SnarlDistanceIndex::START_END,
                SnarlDistanceIndex::CHAIN_HANDLE);

        payload.node_length = get_length(0);
        payload.is_trivial_chain = true;
        payload.is_reversed = false;
        payload.parent_handle = distance_index.get_root();
        payload.parent_type = ZipCode::ROOT_NODE;
        payload.parent_record_offset = 0;

    } else if (decoder[max_depth() - 1].first) {
        //If the parent is a chain
        payload.node_handle = distance_index.get_node_net_handle(id);
        payload.parent_is_chain = true;
        payload.parent_is_root = false;

        size_t parent_depth = max_depth() - 1;

        if (decoder_length() == 2) {
            //If the node is a child of the root chain
            payload.parent_handle = distance_index.start_end_traversal_of(
                                        distance_index.get_handle_from_connected_component(get_distance_index_address(0)));
            payload.parent_type = ZipCode::ROOT_CHAIN;
            payload.parent_is_root = true;
        } else {
            payload.parent_handle = distance_index.start_end_traversal_of(distance_index.get_parent(payload.node_handle));
            payload.parent_type = ZipCode::CHAIN;
        }
        payload.parent_record_offset = distance_index.get_record_offset(payload.parent_handle);

        payload.prefix_sum = get_offset_in_chain(parent_depth+1);

        //Node length
        payload.node_length = get_length(parent_depth+1);
        //is_reversed
        //TODO: For top-level chains we got this from the distance index
        payload.is_reversed = get_is_reversed_in_parent(parent_depth+1);

        payload.chain_component = get_chain_component(parent_depth+1);



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


        if (payload.parent_is_root) {
            //is_chain
            payload.node_handle = payload.parent_handle;
            payload.parent_record_offset = distance_index.get_record_offset(
                                                distance_index.get_handle_from_connected_component(
                                                    get_distance_index_address(0)));
            payload.parent_handle = distance_index.get_net_handle_from_values(payload.parent_record_offset,
                                            SnarlDistanceIndex::START_END,
                                            SnarlDistanceIndex::ROOT_HANDLE);
            payload.parent_type = ZipCode::ROOT_SNARL;
        } else {
            size_t parent_depth = max_depth() - 1;
            payload.parent_type = get_code_type(parent_depth);

            payload.prefix_sum = 0;

            payload.chain_component = 0;

            if (payload.parent_type == ZipCode::REGULAR_SNARL) {
                //Snarl is reversed
                net_handle_t grandparent_handle = distance_index.get_parent(payload.parent_handle);
                //Simple and regular snarls are different for clustering
                if (distance_index.is_simple_snarl(grandparent_handle)) {
                    payload.is_reversed = get_is_reversed_in_parent(parent_depth+1);
                    payload.parent_is_chain=true;
                    payload.parent_record_offset = distance_index.get_record_offset(distance_index.get_parent(grandparent_handle));
                } else {
                    payload.is_reversed = false;
                    payload.parent_record_offset = distance_index.get_record_offset(grandparent_handle);
                }

            } else {
                payload.is_reversed = false;
                payload.parent_record_offset = get_distance_index_address(parent_depth);
            }

        }
        payload.node_length = get_length(max_depth());

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
        result += (decoder[d].first ? "1" : "0");
        if (d == 0) {
            //Root structure
            result += std::to_string(get_distance_index_address(0));
        } else if (decoder[d].first) {
            //is_chain so could be a chain or a node
            if (decoder[d-1].first) {
                //If the thing before this was also a chain, then it is a node
                result += std::to_string(get_offset_in_chain(d));
            } else {
                //Otherwise it's a chain
                result += std::to_string(get_rank_in_snarl(d));
            }
        } else {
            //Definitely a snarl
            result += std::to_string(get_offset_in_chain(d));
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
        throw std::runtime_error("error: trying to get the parent of the root net_identifier_t");
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
