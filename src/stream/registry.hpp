#ifndef VG_STREAM_REGISTRY_HPP_INCLUDED
#define VG_STREAM_REGISTRY_HPP_INCLUDED

/**
 * \file registry.hpp
 * Handles bookkeeping for the various data type tags in stream files.
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

// We define some functional-programming-style types to build the type of a loader function.
// A loader function takes a callback that loops over incoming messages. It then returns a pointer to the loaded object.
// The callback loops over incoming messages by itself calling a callback.
// This lets you get per-load state (as locals in the outer function) without defining a real loader interface and a bunch of implementation child classes that are only used in one place.

using message_consumer_function_t = function<void(const string&)>;
using message_sender_function_t = function<void(const message_consumer_function_t&)>;
using load_function_t = function<void*(const message_sender_function_t&)>;


/**
 * We also have an adapter that takes a function from an istream& to a void*
 * object, and runs that in a thread to adapt it to the message consuming shape
 * of interface. It captures the wrapped function by value.
 */
load_function_t wrap_stream_loader(function<void*(istream&)> istream_loader);

/**
 * Register a function to be used to load the given type of object from the
 * given tag. Like the Subcommand class, you should statically construct one of
 * these, and it will register itself when linked into the binary. Unless your
 * C++ compiler doesn't actually bother making the static objects, which is
 * technically allowed.
 */
template<typename Loaded>
class Loader {
public:
    /**
     * Register the given function to load objects of our template type from
     * message streams with the given tag. The registration will be permanent,
     * even if the Loader object is destroyed. The load function will be copied
     * and so should not reference anything that isn't static/global.
     */
    Loader(const string& tag, load_function_t load_function);
};

/**
 * Register a tag for a particular Protobuf type.
 * Like the Loader, ought to be constructed as a static global at file scope.
 */
template<typename Message>
class ProtobufTag {
public:
    /**
     * Register the given tag as associated with our Protobuf Message type.
     */
    ProtobufTag(const string& tag);
};

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
     * Register a Protobuf Message type to be associated with the given string tag.
     * By default, Protobuf types use the name string from the Protobuf library as a tag.
     * But this can and should be overridden with something short that will never change.
     */
    template<typename Message>
    static void register_protobuf(const string& tag);
    
    /**
     * Register a loader function to load the given type from messages with the given tag.
     * Automatically handles required cleverness for registering for base classes if warranted.
     * This is the backend for the Loader class, which you should user instead of calling this directly.
     */
    template<typename Loaded>
    static void register_loader(const string& tag, load_function_t loader);
    
    /////////
    // Lookup functions
    /////////
    
    /**
     * Get the correct tag to use when serializing Protobuf messages of the given type.
     */
    template<typename Message>
    static const string& get_protobuf_tag();

    /**
     * Check to see if the given tag is expected when deserializing Protobuf messages of the given tag.
     */
    template<typename Message>
    static bool check_protobuf_tag(const string& tag);

    /**
     * Look up the appropriate loader function to use to load an object of the given type from data with the given tag.
     * If there is one registered, return a pointer to it. The caller has to call it and cast the result to the right type.
     * If there isn't, returns nullptr.
     */
    template<typename Wanted>
    static const load_function_t* find_loader(const string& tag);
    
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
    };
    
    /**
     * Get or create the registry tables in which things are registerd.
     */
    static Tables& get_tables();
    
    
};

/////////////
// Template implementations
/////////////

template<typename Loaded>
Loader<Loaded>::Loader(const string& tag, load_function_t load_function) {
    // Pass the load function off to the registry, which will keep a copy.
    Registry::register_loader<Loaded>(tag, load_function);
}

template<typename Message>
ProtobufTag<Message>::ProtobufTag(const string& tag) {
    // Pass the Protobuf type tag off to the registry
    Registry::register_protobuf<Message>(tag);
}

template<typename Message>
void Registry::register_protobuf(const string& tag) {
    // Get our state
    Tables& tables = get_tables();

    // Register in both directions
    tables.tag_to_protobuf.emplace(tag, type_index(typeid(Message)));
    tables.protobuf_to_tag.emplace(type_index(typeid(Message)), tag);
    
}

template<typename Loaded>
void Registry::register_loader(const string& tag, load_function_t loader) {
    // Get our state
    Tables& tables = get_tables();
    
    // Save the loading function to load the given type
    tables.tag_to_loader[tag][type_index(typeid(Loaded))] = loader;
    
    // Special checks: if the class we got was a HandleGraph implementation, register for the base HandleGraph types.
    // TODO: Don't copy the loader unnecessarily.
    
    if (is_base_of<HandleGraph, Loaded>::value) {
        // It implements HandleGraph, so it can get us a HandleGraph
        tables.tag_to_loader[tag][type_index(typeid(HandleGraph))] = loader;
    }
    
    if (is_base_of<PathHandleGraph, Loaded>::value) {
        // It implements PathHandleGraph, so it can get us a PathHandleGraph
        tables.tag_to_loader[tag][type_index(typeid(PathHandleGraph))] = loader;
    }
    
}

template<typename Message>
const string& Registry::get_protobuf_tag() {
    // Get our state
    Tables& tables = get_tables();

    // See if we have a tag defined
    auto found = tables.protobuf_to_tag.find(type_index(typeid(Message)));
    
    if (found != tables.protobuf_to_tag.end()) {
        // There is a custom tag registered.
        return found->second;
    } else {
        // Use the default name from the Protobuf library as a tag
        return Message::descriptor()->full_name();
    }
}

template<typename Message>
bool Registry::check_protobuf_tag(const string& tag) {
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
