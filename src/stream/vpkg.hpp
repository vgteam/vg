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
     * Returns a tuple of pointers to loaded objects, or null if they could not be found.
     * Only works on VPKG-formatted files.
     */
    template<typename... Wanted>
    static tuple<unique_ptr<Wanted>...> load_all(istream& in) {
         // Make an iterator
        MessageIterator it(in);
        
        // Create a collection of null void*s that will hold the allocated objects we want to load when we think we can load them.
        deque<void*> to_fill { (void*)(Wanted*)nullptr... };
        
        // We set this to false if we should stop calling the loaders because we hit EOF or they all loaded.
        bool keep_going = true;

        do {

            // We exploit initializer list evaluation order to be able to tell
            // individual calls resulting from a ... variadic template argument
            // expansion what number they are, so they can index into a data
            // structure. See https://stackoverflow.com/a/21194071
            size_t index = 0;
            
            // Call the load function for each type, and get the statuses
            vector<LoadProgress> load_statuses { load_into_one<Wanted>(it, index++, to_fill)... };
            
            // We keep this true if we should skip the current tag because all
            // the loaders looked at and rejected it with LOAD_SEARCHING or are
            // LOAD_FINISHED. So we set it false if anybody is LOAD_SUCCESS.
            bool skip_tag = true;
            
            // We keep this false if all the loaders are LOAD_FINISHED.
            keep_going = false;
            
            for (auto status : load_statuses) {
                if (status != LOAD_FINISHED) {
                    // Somebody wants to keep going
                    keep_going = true;
                }
                
                if (status == LOAD_SUCCESS) {
                    // We made progress, so don't skip the tag.
                    skip_tag = false;
                }
            }
            
            if (keep_going && skip_tag) {
                // Nobody made progress, so skip all messages tagged with the current tag
                string to_skip = (*it).first;
                while (it.has_next() && (*it).first == to_skip) {
                    ++it;
                }
            }
            
            // Loop until everybody is quiescent.
        } while (keep_going);
        
        // Now all the unique_ptrs that can be filled in are filled in.
        // Convert to a tuple and return.
        return to_tuple<Wanted...>(to_fill);
    }
    
    /**
     * Load an object of the given type form a stream.
     * The stream may be VPKG with the appropriate tag, or a bare non-VPKG stream understood by the loader.
     * Tagged messages that can't be used to load the thing we are looking for are skipped.
     */
    template<typename Wanted>
    static unique_ptr<Wanted> load_one(istream& in) {
        // Check if the thing we want can be loaded from a bare stream
        auto* bare_loader = Registry::find_bare_loader<Wanted>();
        
        if (bare_loader != nullptr) {
        
            // If so, sniff out whether this is actually a BGZF-compressed
            // type-tagged message file.
            
            if (!BlockedGzipInputStream::SmellsLikeGzip(in)) {
                // If it is not BGZF-compressed, try loading it directly with the loader.
                return unique_ptr<Wanted>((Wanted*)(*bare_loader)(in));
            }
            
            // TODO: We assume that if it starts with the GZIP magic number
            // (0x1F 0x8B) it is GZIP and therefore BGZF type-tagged message
            // data. For some of our old formats that didn't include their own
            // leading magic numbers (GCSA, LCP), this assumption may break
            // down for some files!
        
        }
        
        // If it is BGZF-compressed, or we don't have a loader from a bare stream, then make the MessageIterator.
        MessageIterator it(in);
        
        while (it.has_next()) {
            // Scan through kinds of tagged messages
            string current_tag = (*it).first;
            
            // See if we have one that has a registered loader for this type.
            auto* loader = Registry::find_loader<Wanted>(current_tag);
            
            if (loader == nullptr) {
                // Skip all these messages with this tag
                while (it.has_next() && (*it).first == current_tag) {
                    ++it;
                }
            } else {
                // Load with it and return a unique_ptr for the result.
                return unique_ptr<Wanted>((Wanted*)(*loader)([&](const message_consumer_function_t& handle_message) {
                    while (it.has_next() && (*it).first == current_tag) {
                        // Feed in messages from the file until we run out or the tag changes
                        handle_message((*it).second);
                        ++it;
                    }
                }));
            }
        }
        
        // If we get here, nothing with an appropriate tag could be found, and it wasn;t a bare loadable file.
        return unique_ptr<Wanted>(nullptr);
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
     * Return type to represent whether the loader is making progress. Lets us
     * know when we have loaded items vs. when nobody has anything to load so
     * we can skip unwanted sections.
     */
    enum LoadProgress {
        LOAD_SUCCESS,
        LOAD_FINISHED,
        LOAD_SEARCHING
    };

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
     * Returns LOAD_SUCCESS if it finds and loads something, LOAD_FINISHED if
     * it is filled in or hits EOF, and LOAD_SEARCHING if it is still unfilled
     * but can't be loaded from the current tag.
     */
    template<typename One>
    static LoadProgress load_into_one(MessageIterator& it, size_t i, deque<void*>& dest) {
        // Find the slot to load into
        void*& slot = dest[i];
        
        if (slot != nullptr) {
            // If it's already loaded, we're done
            return LOAD_FINISHED;
        }
        
        if (!it.has_next()) {
            // If there's nothing to look at, we're done
            return LOAD_FINISHED;
        }
        
        // Grab and cache the tag
        string tag_to_load = (*it).first;
        
        // Look for a loader in the registry based on the tag.
        auto* loader = Registry::find_loader<One>(tag_to_load);
        
        if (loader == nullptr) {
            // We can't load from this. Try again later.
            return LOAD_SEARCHING;
        }
        
        // Otherwise we can load, so do it.
        slot = (*loader)([&](const message_consumer_function_t& handle_message) {
            while (it.has_next() && (*it).first == tag_to_load) {
                // Feed in messages from the file until we run out or the tag changes
                handle_message((*it).second);
                ++it;
            }
        });
        
        // Say we loaded something
        return LOAD_SUCCESS;
    }
};

}

}

#endif
