#ifndef VG_STREAM_VPKG_HPP_INCLUDED
#define VG_STREAM_VPKG_HPP_INCLUDED

/**
 * \file vpkg.hpp: frontend load/save interface for multi-type type-tagged files
 */


#include "registry.hpp"
#include "message_iterator.hpp"

#include <iostream>
#include <tuple>
#include <vector>
#include <memory>

namespace vg {

namespace stream {

using namespace std;

/**
 * Interface for reading/writing type-tagged files.
 *
 * Allows you to load a HandleGraph from a file with the implementation being
 * auto-selected based on what is in the file.
 *
 * Allows you to load multiple indexes or objects from the same file, as they
 * are found.
 *
 * Allows saving multiple indexes or objects to a file with one call.
 */
class VPKG {
public:
    /**
     * Allocate and load one or more objects from the given stream on a single pass.
     * Retuens a tuple of pointers to loaded objects, or null if they could not be found.
     */
    template<typename... Wanted>
    static tuple<unique_ptr<Wanted>...> load_all(istream& in);
    
private:

    /**
     * If the one item of type One at index i in the destination tuple can be
     * filled from the given MessageIterator, and is empty, fill it.
     * 
     * Returns false if it can't be filled, or is already filled, or the
     * iterator is over.
     */
    template<typename One, typename... TupleTypes>
    static bool load_into_one(MessageIterator& it, size_t i, tuple<unique_ptr<TupleTypes>...>& dest);
};

/////////////
// Template implementations
/////////////

template<typename... Wanted>
auto VPKG::load_all(istream& in) -> tuple<unique_ptr<Wanted>...> {
    // Make an iterator
    MessageIterator it(in);
    
    // Make a destination tuple
    tuple<unique_ptr<Wanted>...> to_return;
    
    bool keep_going = false;

    do {

        // We exploit initializer list evaluation order to be able to tell
        // individual calls resulting from a ... variadic template argument
        // expansion what number they are, so they can index into a tuple.
        // See https://stackoverflow.com/a/21194071
        
        size_t tuple_index = 0;
        
        // Call the load function for each type, and get the statuses
        vector<bool> load_statuses = {load_into_one<Wanted, Wanted...>(it, tuple_index++, to_return)...};
        
        for (bool status : load_statuses) {
            // OR together all the statuses so we know if we need to continue for anything.
            keep_going |= status;
        }
        
    } while (keep_going);
    
    // Now all the unique_ptrs that can be filled in are filled in
    return to_return;
}

template<typename One, typename... TupleTypes>
auto load_into_one(MessageIterator& it, size_t i, tuple<unique_ptr<TupleTypes>...>& dest) -> bool {
    
    // Find the pointer to load
    unique_ptr<One>& ptr = get<i>(dest);
    
    if (ptr.get() != nullptr) {
        // If it's already loaded, we're done
        return false;
    }
    
    if (!it.has_next()) {
        // If there's nothing to look at, we're done
        return false;
    }
    
    // Grab and cache the tag
    string tag_to_load = (*it).first;
    
    // Look for a loader in the registry based on the tag.
    auto* loader = Registry::find_loader<One>(tag_to_load);
    
    if (loader == nullptr) {
        // We can't load from this. Try again later.
        return true;
    }
    
    // Otherwise we can load, so do it.
    ptr = unique_ptr<One>((*loader)([&](const message_consumer_function_t& handle_message) {
        while (it.has_next() && (*it).first == tag_to_load) {
            // Feed in messages from the file until we run out or the tag changes
            handle_message((*it).second);
            ++it;
        }
    }));
    
    // Now there's nothing left to load
    return false;
}

}

}

#endif
