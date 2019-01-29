#ifndef VG_STREAM_REGISTRY_HPP_INCLUDED
#define VG_STREAM_REGISTRY_HPP_INCLUDED

/**
 * \file registry.hpp
 * Handles bookkeeping for the various data type tags in stream files.
 *
 * TO ADD A PROTOBUF TYPE:
 * - Add a register_protobuf<Type>("TAG");
 * - Call in Registry::register_everything() in registry.cpp
 * TO ADD A NON-PROTOBUF LOADER/SAVER:
 * - Create a register_loader_saver_whatever.cpp/.hpp defining register_loader_savert_whatever().
 * - Make it call Registry::register_loader_saver<Type>(tag, load_function, save_function)
 * - Call register_loader_saver_whatever() in Registry::register_everything() in registry.cpp
 * TO LOAD/SAVE SOMETHING:
 */

#include <string>
#include <unordered_map>
#include <typeinfo>
#include <typeindex>
#include <type_traits>
#include <functional>
#include <iostream>

#include "../handle.hpp"

namespace vg {

namespace stream {

using namespace std;

// We *want* to use magic registration like we do in the Subcommand system,
// where you just drop a cpp file in the project and it registers load/save
// functions for HandleGraph implementations. But that won't work by default,
// because we're building this stuff into libvg.a, and .o files in the library
// that aren't actually referenced won't get pulled into the final executable.
// You would have to use --whole-archive (gcc) or -force_load (clang) to pull
// them all in. See
// <https://www.bfilipek.com/2018/02/static-vars-static-lib.html> and
// <https://stackoverflow.com/questions/805555/ld-linker-question-the-whole-archive-option/842770#842770>.


// So all the registration actually happens in registry.cpp. It includes
// headers for all the things we can load, and calls out manually to register
// them. Each header defines a register_whatever() function.

// We define some functional-programming-style types to build the type of a
// loader function.
//
// A loader function takes a callback that loops over incoming messages. It
// then returns a pointer to the loaded object.
//
// The callback loops over incoming messages by itself calling a callback.
//
// This lets you get per-load state (as locals in the outer function) without
// defining a real loader interface and a bunch of implementation child classes
// that are only used in one place.

/// This is the type of a function that can be fed a series of messages
using message_consumer_function_t = function<void(const string&)>;
/// This is the type of a function that can be given a message consumer to feed messages to.
using message_sender_function_t = function<void(const message_consumer_function_t&)>;
/// This is the type of a function that can allocate and load an object of unspecified type from a message source.
using load_function_t = function<void*(const message_sender_function_t&)>;

/// This is the type of a function that can serialize an object of unspecified type to a message consumer.
using save_function_t = function<void(const void*, const message_consumer_function_t&)>;

/// This is the type of a function that can load an object of unspecified type from a bare input stream.
using bare_load_function_t = function<void*(istream&)>;

/// This is the type of a function that can save an object of unspecified type to a bare output stream.
using bare_save_function_t = function<void(const void*, ostream&)>;

/**
 * We also have an adapter that takes a function from an istream& to a void*
 * object, and runs that in a thread to adapt it to the message consuming shape
 * of interface. It captures the wrapped function by value.
 */
load_function_t wrap_bare_loader(bare_load_function_t istream_loader);

/**
 * And an adapter that takes a function of void* and ostream&, and adapts that
 * to a message consumer destination.
 */
save_function_t wrap_bare_saver(function<void(const void*, ostream&)> ostream_saver);

/**
 * A registry mapping between tag strings for serialized message groups and the
 * Protobuf types or serialization/deserialization handlers to use for them.
 * All registration anmd lookup is done through static methods.
 * Static methods are safe to call from other static initialization code.
 *
 * We handle two kinds of registration right now:
 * - Registration of Protobuf types so we can tag them appropriately on serialization and
 *   detect/skip mismatches on deserialization
 * - Registration of Loader implementations for chunked binary files, so we can
 *   select the appropriate implementation if we are asked to read a HandleGraph
 *   or an LCP or whatever from a stream file.
 */
class Registry {
public:

    /////////
    // Registration functions
    /////////
    
    /**
     * Register everything. Main entry point.
     * Returns true on success.
     */
    static bool register_everything();

    /**
     * Register a Protobuf Message type to be associated with the given string
     * tag. By default, Protobuf types use the name string from the Protobuf
     * library as a tag. But this can and should be overridden with something
     * short that will never change.
     */
    template<typename Message>
    static void register_protobuf(const string& tag);
    
    /**
     * Register a loading function and a saving function with the given tag
     * for the given object type.
     */
    template<typename Handled>
    static void register_loader_saver(const string& tag, load_function_t loader, save_function_t saver);
    
    /**
     * Register a loading function and a saving function with the given
     * collection of tags for the given object type. The first tag in the
     * collection will be used for saving. If "" appears in the list of tags,
     * the loader can be deployed on untagged message groups (for backward compatibility).
     */
    template<typename Handled>
    static void register_loader_saver(const vector<string>& tags, load_function_t loader, save_function_t saver);
    
    /**
     * Register a loading function and a saving function with the given tag for
     * the given object type. The functions operate on bare streams; conversion
     * to type-tagged messages of chunks of stream data is performed
     * automatically. The load function will also be registered to load from
     * non-type-tagged-message-format files, for backward compatibility.
     */
    template<typename Handled>
    static void register_bare_loader_saver(const string& tag, bare_load_function_t loader, bare_save_function_t saver);
    
    /////////
    // Lookup functions
    /////////
    
    /**
     * Determine if the given tag string loaded from a file is a valid/possible
     * tag, or if it should be interpreted as message data from a pre-tags VG
     * file instead. Only tag values literally registered with the registry are
     * valid.
     * NOT thread-safe to run simultaneously with tag registrations.
     */
    static bool is_valid_tag(const string& tag);
    
    /**
     * Get the correct tag to use when serializing Protobuf messages of the
     * given type.
     */
    template<typename Message>
    static const string& get_protobuf_tag();

    /**
     * Check to see if the given tag is expected when deserializing Protobuf
     * messages of the given tag.
     */
    template<typename Message>
    static bool check_protobuf_tag(const string& tag);

    /**
     * Look up the appropriate loader function to use to load an object of the
     * given type from data with the given tag. If there is one registered,
     * return a pointer to it. The caller has to call it and cast the result to
     * the right type. If there isn't, returns nullptr.
     */
    template<typename Want>
    static const load_function_t* find_loader(const string& tag);
    
    /**
     * Look up the appropriate loader function to use to load an object of the
     * given type from a bare stream. If there is one registered, return a
     * pointer to it. The caller has to call it and cast the result to the
     * right type. If there isn't, returns nullptr.
     */
    template<typename Want>
    static const bare_load_function_t* find_bare_loader();
    
    /**
     * Look up the appropriate saver function to use to save an object of the
     * given type. If there is one registered, return a pointer to a pair of
     * the tag to use and the function. The caller has to call it and cast the
     * result to the right type. If there isn't, returns nullptr.
     */
    template<typename Have>
    static const pair<string, save_function_t>* find_saver();
    
private:
    
    /**
     * Holds the actual singleton registry tables.
     */
    struct Tables {
        /// Maps from tag string to Protobuf type type_index that it indicates, if any.
        unordered_map<string, type_index> tag_to_protobuf;
        /// Maps from Protobuf type type_index back to string tag.
        unordered_map<type_index, string> protobuf_to_tag;
        
        /// Maps from tag to a map from type_index we want to load to a loading
        /// function that can load it from data with that tag.
        unordered_map<string, unordered_map<type_index, load_function_t>> tag_to_loader;
        /// Maps from type to a single tag and save function pair to use when outputting that type.
        unordered_map<type_index, pair<string, save_function_t>> type_to_saver;
        
        /// Maps from type_index we want to load from a old,
        /// non-tagged-message-format file to a "bare" loader taht can load the
        /// desired thing from an istream.
        unordered_map<type_index, bare_load_function_t> type_to_bare_loader;
    };
    
    /**
     * Get or create the registry tables in which things are registerd.
     */
    static Tables& get_tables();
    
    /**
     * Register a load function for a tag. The empty tag means it can run on
     * untagged message groups.
     */
    template<typename Handled>
    static void register_loader(const string& tag, load_function_t loader);
    
    /**
     * Register a load function for loading a type from non-type-tagged-message
     * "bare" streams.
     */
    template<typename Handled>
    static void register_bare_loader(bare_load_function_t loader);
    
    /**
     * Register a save function to save a type with a given tag. The empty tag
     * is not permitted.
     */
    template<typename Handled>
    static void register_saver(const string& tag, save_function_t saver);
};

/////////////
// Template implementations
/////////////

template<typename Message>
void Registry::register_protobuf(const string& tag) {
    // Get our state
    Tables& tables = get_tables();

    // Register in both directions
    tables.tag_to_protobuf.emplace(tag, type_index(typeid(Message)));
    tables.protobuf_to_tag.emplace(type_index(typeid(Message)), tag);
    
#ifdef debug
    cerr << "Registered " << Message::descriptor()->full_name() << " as " << tag << endl;
#endif
}

template<typename Handled>
void Registry::register_loader(const string& tag, load_function_t loader) {
    // Get our state
    Tables& tables = get_tables();
    
    // Save the loading function to load the given type
    tables.tag_to_loader[tag][type_index(typeid(Handled))] = loader;
    
    // Special checks: if the class we got was a HandleGraph implementation, register for the base HandleGraph types.
    // TODO: Don't copy the loader unnecessarily.
    
    if (is_base_of<HandleGraph, Handled>::value) {
        // It implements HandleGraph, so it can get us a HandleGraph
        tables.tag_to_loader[tag][type_index(typeid(HandleGraph))] = loader;
    }
    
    if (is_base_of<PathHandleGraph, Handled>::value) {
        // It implements PathHandleGraph, so it can get us a PathHandleGraph
        tables.tag_to_loader[tag][type_index(typeid(PathHandleGraph))] = loader;
    }
}

template<typename Handled>
void Registry::register_bare_loader(bare_load_function_t loader) {
    // Get our state
    Tables& tables = get_tables();
    
    // Save the function in the table
    tables.type_to_bare_loader.emplace(type_index(typeid(Handled)), loader);
}

template<typename Handled>
void Registry::register_saver(const string& tag, save_function_t saver) {
    // Prohibit the empty tag here.
    assert(!tag.empty());
    
    // Get our state
    Tables& tables = get_tables();
    
    // Save the saving function to save the given type
    tables.type_to_saver.emplace(type_index(typeid(Handled)), make_pair(tag, saver));
}

template<typename Handled>
void Registry::register_loader_saver(const string& tag, load_function_t loader, save_function_t saver) {
    // Dispatch to the vector implementation
    register_loader_saver<Handled>(vector<string>{tag}, loader, saver);
}

template<typename Handled>
void Registry::register_loader_saver(const vector<string>& tags, load_function_t loader, save_function_t saver) {
    // There must be tags
    assert(!tags.empty());

    // The first must be a real tag we can save with
    assert(!tags.front().empty());
    
    // The first tag gets the loader and saver
    register_loader<Handled>(tags.front(), loader);
    register_saver<Handled>(tags.front(), saver);
    
    for (size_t i = 1; i < tags.size(); i++) {
        // Other tags just get loaders
        register_loader<Handled>(tags.front(), loader);
    }
}

template<typename Handled>
void Registry::register_bare_loader_saver(const string& tag, bare_load_function_t loader, bare_save_function_t saver) {

    // Register the type-tagged wrapped functions
    register_loader_saver<Handled>(tag, wrap_bare_loader(loader), wrap_bare_saver(saver));
    
    // Register the bare stream loader
    register_bare_loader<Handled>(loader);

}

template<typename Want>
const load_function_t* Registry::find_loader(const string& tag) {
    // Get our state
    Tables& tables = get_tables();
    
    auto found_tag = tables.tag_to_loader.find(tag);
    if (found_tag != tables.tag_to_loader.end()) {
        // We can load this tag to something.
        // Grab the map from type_index to actual loader
        auto& loaders = found_tag->second;
        
        auto found_loader = loaders.find(type_index(typeid(Want)));
        if (found_loader != loaders.end()) {
            // We can load this tag to the requested type.
            // Return a pointer to the function that does it.
            return &found_loader->second;
        }
    }
    
    // We can't find the right function. Return null.
    return nullptr;
}

template<typename Want>
const bare_load_function_t* Registry::find_bare_loader() {
    // Get our state
    Tables& tables = get_tables();
    
    // Look for a loader for this type from bare streams
    auto found = tables.type_to_bare_loader.find(type_index(typeid(Want)));
    
    if (found != tables.type_to_bare_loader.end()) {
        // We found one. Return a pointer to it.
        return &found->second;
    }
    
    // We don't have a loader to load this from a bare file.
    return nullptr;
}
    

template<typename Have>
const pair<string, save_function_t>* Registry::find_saver() {
    // Get our state
    Tables& tables = get_tables();
   
    // Look for a saver for this templated type.
    auto found = tables.type_to_saver.find(type_index(typeid(Have)));
    if (found != tables.type_to_saver.end()) {
        // We only have one. Return a pointer to the pair of the tag to apply and the function to call.
        return &found->second;
    }
    
    // Otherwise we didn't find anything to use to save this.
    return nullptr;
}

template<typename Message>
const string& Registry::get_protobuf_tag() {
    // Get our state
    Tables& tables = get_tables();

    // See if we have a tag defined
    auto found = tables.protobuf_to_tag.find(type_index(typeid(Message)));
    
    if (found != tables.protobuf_to_tag.end()) {
        // There is a custom tag registered.
#ifdef debug
        cerr << "Tag found for type " << Message::descriptor()->full_name() << ": " << found->second << endl;
#endif
        return found->second;
    } else {
        // Use the default name from the Protobuf library as a tag
#ifdef debug
        cerr << "Tag not found for " << Message::descriptor()->full_name() << endl;
#endif
        return Message::descriptor()->full_name();
    }
}

template<typename Message>
bool Registry::check_protobuf_tag(const string& tag) {
    if (tag.empty()) {
        // For reading old tagless files, "" is always a valid tag for Protobuf data.
        return true;
    }

    // Get our state
    Tables& tables = get_tables();
    
    // See if we have a protobuf defined for this tag as an override.
    auto found = tables.tag_to_protobuf.find(tag);
    
    if (found != tables.tag_to_protobuf.end()) {
        // We do have a Protobuf specifically assigned to this tag
        
        // Return true iff it is the same Protobuf type we are checking
        return type_index(typeid(Message)) == found->second;
    } else {
        // Return if the tag is the Protobuf type's fully qualified name
        return tag == Message::descriptor()->full_name();
    }
}

}

}

#endif
