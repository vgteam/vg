#ifndef VG_STREAM_VPKG_HPP_INCLUDED
#define VG_STREAM_VPKG_HPP_INCLUDED

/**
 * \file vpkg.hpp: frontend load/save interface for multi-type type-tagged files
 */


#include "registry.hpp"
#include "message_iterator.hpp"
#include "message_emitter.hpp"

#include <iostream>
#include <tuple>
#include <vector>
#include <deque>

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
    static tuple<unique_ptr<Wanted>...> load_all(istream& in) {
         // Make an iterator
        MessageIterator it(in);
        
        // Create a collection of null void*s that will hold the allocated objects we want to load when we think we can load them.
        deque<void*> to_fill { (void*)(Wanted*)nullptr... };
        
        bool keep_going = false;

        do {

            // We exploit initializer list evaluation order to be able to tell
            // individual calls resulting from a ... variadic template argument
            // expansion what number they are, so they can index into a data
            // structure. See https://stackoverflow.com/a/21194071
            
            size_t index = 0;
            
            // Call the load function for each type, and get the statuses
            vector<bool> load_statuses { load_into_one<Wanted>(it, index++, to_fill)... };
            
            for (bool status : load_statuses) {
                // OR together all the statuses so we know if we need to continue for anything.
                keep_going |= status;
            }
            
        } while (keep_going);
        
        // Now all the unique_ptrs that can be filled in are filled in.
        // Convert to a tuple and return.
        return to_tuple<Wanted...>(to_fill);
    }
    
    /**
     * Save an object to the given stream, using the appropriate saver.
     */
    template<typename Have>
    static void save(const Have& have, ostream& out) {
        // Look for a saver in the registry
        auto* tag_and_saver = Registry::find_saver<Have>();
        
        // We shouldn't ever be saving something we don't know how to save.
        assert(tag_and_saver != nullptr);
        
        // Make an emitter to emit tagged messages
        MessageEmitter emitter(out);
        
        // Start the save
        tag_and_saver->second((const void*)&have, [&](const string& message) {
            // For each message that we have to output during the save, output it via the emitter with the selected tag.
            // TODO: We copy the data string.
            emitter.write_copy(tag_and_saver->first, message);
        });
    }
    
private:

    /**
     * Given a collection of void pointers, give ownership of the objects they point to, if any, to unique_ptrs in a tuple.
     */
    template<typename... TupleTypes>
    static tuple<unique_ptr<TupleTypes>...> to_tuple(deque<void*> items) {
        // Use initializer list expansion to repeatedly pop the first thing off the collection and type it correctly.
        tuple<unique_ptr<TupleTypes>...> to_return { extract_first<TupleTypes>(items)... };
        return to_return;
    }
    
    /**
     * Pop off the first item in the given collection and wrap it in a typed unique_ptr.
     */
    template<typename Pointed>
    static unique_ptr<Pointed> extract_first(deque<void*>& pointers) {
        // Grab off the first thing
        void* got = pointers.front();
        pointers.pop_front();
        // Wrap it in a properly typed unique_ptr;
        return unique_ptr<Pointed>((Pointed*) got);
    }

    /**
     * If the null slot at index i in the given collection of void*s can be
     * filled with an object of type One from the given MessageIterator, fill
     * it.
     * 
     * Returns false if it can't be filled, or is already filled, or the
     * iterator is over.
     */
    template<typename One>
    static bool load_into_one(MessageIterator& it, size_t i, deque<void*>& dest) {
        // Find the slot to load into
        void*& slot = dest[i];
        
        if (slot != nullptr) {
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
        slot = (*loader)([&](const message_consumer_function_t& handle_message) {
            while (it.has_next() && (*it).first == tag_to_load) {
                // Feed in messages from the file until we run out or the tag changes
                handle_message((*it).second);
                ++it;
            }
        });
        
        // Now there's nothing left to load
        return false;
    }
};

}

}

#endif
