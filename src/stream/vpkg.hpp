#ifndef VG_STREAM_VPKG_HPP_INCLUDED
#define VG_STREAM_VPKG_HPP_INCLUDED

/**
 * \file vpkg.hpp: frontend load/save interface for multi-type type-tagged files
 */


#include "registry.hpp"
#include "message_iterator.hpp"

#include <iostream>
#include <tuple>
#include <list>
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
     * If any type in Wanted can be loaded from the message the given iterator
     * is on, load it into the corresponding empty unique_ptr<T> in the given
     * list.
     *
     * Returns true if somethign was loaded, or if something is still wanted
     * but the current tag can't load anything.
     *
     * Returns false when all the unique_ptrs have been fileld, or the iterator
     * hits EOF.
     */
    template<typename... Wanted>
    static bool load_into(MessageIterator& it, list<void*>::iterator unique_ptrs_begin, list<void*>::iterator unique_ptrs_end);
    
    template<typename One>
    static bool load_into_one(MessageIterator& it, list<void*>::iterator unique_ptrs_begin, list<void*>::iterator unique_ptrs_end);

    /**
     * We need a way to scan over the tuple entries, but std::tuple can't do
     * first/rest. So we have this tool to allocate us a bunch of empty
     * unique_ptrs of various types and return a list of pointers to them.
     */
    template<typename... Wanted>
    static list<void*> allocate_unique_ptrs();
    
    /**
     * And then we have this to turn those void*s back into a tuple of
     * unique_ptrs, by moving and de-allocating the ones in the list.
     */
    template<typename... Wanted>
    static tuple<unique_ptr<Wanted>...> deallocate_unique_ptrs(list<void*>::iterator unique_ptrs_begin,
        list<void*>::iterator unique_ptrs_end);
};

/////////////
// Template implementations
/////////////

template<typename... Wanted>
auto VPKG::load_all(istream& in) -> tuple<unique_ptr<Wanted>...> {

    // Allocate all the unique pointers
    list<void*> unique_ptr_ptrs = allocate_unique_ptrs<Wanted...>();
    
    // Make an iterator
    MessageIterator it(in);
    
    // Repeatedly recurse over all our types.
    // When we find one that we can read from the current tag, read it and make the appropriate unique_ptr own it
    // Repeat until we run out of file or all our unique_ptrs are filled.
    while(load_into<Wanted...>(it, unique_ptr_ptrs.begin(), unique_ptr_ptrs.end())) {
        // Nothing to do!
    }
    
    // Go from unique pointers floating to unique pointers in a tuple
    return deallocate_unique_ptrs<Wanted...>(unique_ptr_ptrs.begin(), unique_ptr_ptrs.end());
}

// Load one thing into a unique_ptr in a void pointer list
template<typename One>
auto VPKG::load_into_one(MessageIterator& it, list<void*>::iterator unique_ptrs_begin,
    list<void*>::iterator unique_ptrs_end) -> bool {
    
    // Find the pointer to load
    unique_ptr<One>& ptr = *((unique_ptr<One>*)*unique_ptrs_begin);
    
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

template<typename... Wanted>
auto VPKG::load_into(MessageIterator& it, list<void*>::iterator unique_ptrs_begin,
    list<void*>::iterator unique_ptrs_end) -> bool {
   
    // TODO: Work out a way to get through the list without relying on argument evaluation order.
    
    return false;
}



// Allocate a possibly empty list of unique_ptrs.
template<typename... Wanted>
auto VPKG::allocate_unique_ptrs() -> list<void*> {
    // Try and use magic espression looping in an initializer list.
    return list<void*>{((void*) new unique_ptr<Wanted>())...};
}

// Move any and all uniqur_ptrs from the void* list to a tuple
template<typename... Wanted>
auto VPKG::deallocate_unique_ptrs(list<void*>::iterator unique_ptrs_begin,
    list<void*>::iterator unique_ptrs_end) -> tuple<unique_ptr<Wanted>...> {
    
    // TODO: Work out a way to get through the list without relying on argument evaluation order.
    
    return tuple<unique_ptr<Wanted>...>();
}

}

}

#endif
