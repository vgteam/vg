#ifndef VG_OPTION_H
#define VG_OPTION_H

/**
 * \file option.hpp: Command-line options that can be attached to configurable thing-doer calsses (variant callers, mappers, etc.)
 *
 * To use: make your configurable object inherit from Configurable, and make all your command-line-option fields into Option<whatever type>.
 */

#include <string>
#include <vector>
#include <set>
#include <map>
#include <functional>
#include <iostream>
#include <sstream>
#include <getopt.h>
 
namespace vg {

using namespace std;

/**
 * All of the option templates inherit from this base class, which the command-
 * line parser uses to feed them strings.
 */
class OptionInterface {
public:
    /**
     * Get the long oiption text without --, like "foos-to-bar".
     */
    virtual const string& get_long_option() const = 0;
    
    /**
     * Gets a list of short option characters that the option wants, in priority
     * order. If none of these are available, the option will be automatically
     * assigned some other free character.
     */
    virtual const string& get_short_options() const = 0;

    /**
     * Get the description, like "number of foos to bar per frobnitz".
     */
    virtual const string& get_description() const = 0;
    
    /**
     * Get the default value as a string. May be generated on the fly.
     */
    virtual string get_default_value() const = 0;
    
    /**
     * Returns true if the option takes an argument, and false otherwise.
     */
    virtual bool has_argument() const = 0;
    
    /**
     * Called for no-argument options when the parser encounters them.
     */
    virtual void parse() = 0;
    
    /**
     * Called for argument-having options when the parser encounters them. The
     * passed reference is only valid during the function call, so the option
     * should make a copy.
     */
    virtual void parse(const string& arg) = 0;
    
    /**
     * Everyone needs a virtual destructor!
     */
    virtual ~OptionInterface() = default;
};

/**
 * Represents a thing-doing class that has options. You construct the class,
 * configure its options, and then set it going.
 */
class Configurable {
public:
    /**
     * Each option will call this on construction and register with its owner.
     */
    virtual void register_option(OptionInterface* option);
    
    /**
     * Get all the options for this class. These pointers are only valid unless
     * or until the underlying Configurable object moves or is assigned to.
     */
    virtual vector<OptionInterface*> get_options();
    
private:
    /// We really should be using something like CRTP and pointer-to-members on
    /// the actual class that's configurable, but for now we'll just exploit
    /// int<->pointer conversion and store the offset from us in memory to the
    /// location of the option in memory. If the option is correctly a mamber of
    /// a class derived from this one, it will work correctly because all
    /// instances share the same memory layout for members.
    vector<ptrdiff_t> option_offsets;
};

/**
 * Actual implementation class that connects a bunch of Configurable things to
 * getopt.
 *
 * TODO: right now every option must have a long option and a char short option.
 * In theory we can use any int as an ID for a long option. We should support
 * long options with untypable or non-char option IDs.
 */
class ConfigurableParser {
public:

    /**
     * Constructor that sets up our magic option assigning code. Makes a parser
     * that will parse the given base Getopt null-terminated long and short
     * options with Getopt and call the given function with each.
     *
     * The base option handler, if specified, should error on unrecognized
     * options.
     */
    ConfigurableParser(const char* base_short_options = nullptr,
        const struct option* base_long_options = nullptr,
        function<void(int)> handle_base_option = [](int c) { throw runtime_error("Invalid option: " + string(1, (char) c)); });
        
    /**
     * Register a Configurable thing and allocate characters for all its
     * options. The Configurable is not allowed to move after registration.
     */
    void register_configurable(Configurable* configurable);
    
    /**
     * Print option descriptions for registered Configurables to the given
     * stream.
     */
    void print_help(ostream& out) const;
    
    /**
     * Parse the options with getopt. Uses the existing value of the global
     * optind, and updates the global optind, just as getopt would.
     */
    void parse(int argc, char** argv);

private:
    /// Holds all the long option structs for options we have
    vector<struct option> long_options;
    /// Holds all the available short option characters that have not been
    /// assigned. Stored as ints so we can check ints against it safely.
    set<int> available_short_options;
    /// Holds the short options string (with colons) for the short options we
    /// have.
    string short_options;
    
    /// Keep a map from assigned character to actual option.
    map<int, OptionInterface*> options_by_code;
    
    /// And a reverse map for when we look up assigned codes for printing.
    /// TODO: No forgery!
    map<OptionInterface*, int> codes_by_option;
    
    /// And a vector of Configurables in the order registered so we can
    /// interrogate them and group their options when printing help.
    vector<Configurable*> configurables;
    
    /// Handle an option not assigned to an Option object.
    function<void(int)> handle_base_option;
};

/**
 * This class holds static methods explaining how to parse a type.
 */
template <typename Value>
class OptionValueParser {
public:
    // In general, try using << and >>.
    
    /**
     * Return true if we need an argument and false otherwise.
     */
    static bool has_argument() {
        return true;
    }
    
    /**
     * Parse from no argument, but a default value.
     */
    static void parse_default(const Value& default_value, Value& value) {
        throw runtime_error("Argument required for option!");
    }
    
    /**
     * Parse from an argument.
     */
    static void parse(const string& arg, Value& value) {
        stringstream s(arg);
        s >> value;
    }
    
    /**
     * Stringify a default value.
     */
    static string unparse(const Value& value) {
        stringstream s;
        s << value;
        return s.str();
    }
};

// For bools we override some things in the basic OptionValueParser.
// If these weren't inline we'd have to put the code in the .cpp file.

/**
 * Bool options don't need arguments.
 */
template<>
inline bool OptionValueParser<bool>::has_argument() {
    return false;
}

/**
 * When someone gives a bool option they mean to invert its default value.
 */
template<>
inline void OptionValueParser<bool>::parse_default(const bool& default_value, bool& value) {
    value = !default_value;
}

/**
 * If someone gives an option to a bool, explode.
 */
template<>
inline void OptionValueParser<bool>::parse(const string& arg, bool& value) {
    throw runtime_error("Boolean options should not take arguments.");
}

/**
 * For vector options, we recurse.
 */
template <typename Item>
class OptionValueParser<vector<Item>> {
public:
    
    /**
     * Return true if we need an argument and false otherwise.
     */
    static bool has_argument() {
        return true;
    }
    
    /**
     * Parse from no argument, but a default value.
     */
    static void parse_default(const vector<Item>& default_value, vector<Item>& value) {
        throw runtime_error("Argument required for option!");
    }
    
    /**
     * Parse from an argument.
     */
    static void parse(const string& arg, vector<Item>& value) {
        Item entry;
        OptionValueParser<Item>::parse(arg, entry);
        value.push_back(entry);
    }
    
    /**
     * Stringify a default value.
     */
    static string unparse(const vector<Item>& value) {
        stringstream s;
        // Enclose the vector in braces
        s << "{";
        for (size_t i = 0; i < value.size(); i++) {
            // Put every item
            s << OptionValueParser<Item>::unparse(value[i]);
            if (i + 1 < value.size()) {
                // Add comma separators
                s << ",";
            }
        }
        s << "}";
        return s.str();
    }
};

/**
 * The correct user entry point is one of the Option<> specializations, not this
 * class.
 *
 * Represents a configurable parameter for a class that we might want to expose
 * on the command line.
 *
 * We need to allow these to be defined with default values in the class's
 * header, but be overridable by simple assignment of a value of the right type.
 * But we also need to be able to interrogate an instance of the class to get
 * all its options that need to be filled in and their help text, and ideally to
 * choose short options to assign to them that don't conflict.
 *
 * So our approach is to have the class keep track of all its option objects,
 * and to have the options register themselves at construction time with a
 * passed this pointer (since this is in scope in the class body in the header).
 *
 * Option instances MUST exist as MEMBERS of classes that inherit from
 * Configurable! They can't live in vectors or anything, or the magic code that
 * tracks them as the enclosing object moves won't work!
 *
 * We'll have an assignment operator from the wrapped type to make setting
 * options manually easy.
 *
 * We also have magical conversion to the wrapped type.
 *
 * We wrap all the type-specific parsing into a parser template. You could
 * specify your own if you want custom parsing logic for like a map or
 * something.
 *
 */
template<typename Value, typename Parser = OptionValueParser<Value>>
class BaseOption : public OptionInterface {

public:

    /**
     * No default constructor.
     */
    BaseOption() = delete;
    
    /**
     * Default destructor.
     */
    virtual ~BaseOption() = default;

    /**
     * Make a new BaseOption that lives in a class, with the given name, short
     * option characters, and default value.
     */
    BaseOption(Configurable* owner, const string& long_opt, const string& short_opts, const Value& default_value,
        const string& description) : long_opt(long_opt), short_opts(short_opts), description(description),
        default_value(default_value), value(default_value) {
        
        // Register with the owner
        owner->register_option(this);
    }
    
    /**
     * Assignment from an Option of the correct type.
     */
    BaseOption<Value, Parser>& operator=(const BaseOption<Value, Parser>& other) = default;
    
    /**
     * Assignment from an unwrapped value.
     */
    BaseOption<Value, Parser>& operator=(const Value& other) {
        // Just use the value type's assignment operator.
        value = other;
        return *this;
    }
    
    /**
     * Conversion to the wrapped type.
     */
    operator Value&() {
        return value;
    }
    
    // Implementation for the option interface follows
    
    /**
     * Get the long oiption text without --, like "foos-to-bar".
     */
    virtual const string& get_long_option() const {
        return long_opt;
    }
    
    /**
     * Gets a list of short option characters that the option wants, in priority
     * order. If none of these are available, the option will be automatically
     * assigned some other free character.
     */
    virtual const string& get_short_options() const {
        return short_opts;
    }

    /**
     * Get the description, like "number of foos to bar per frobnitz".
     */
    virtual const string& get_description() const {
        return description;
    }
    
    /**
     * Get the default value as a string.
     */
    virtual string get_default_value() const {
        return Parser::unparse(default_value);
    };
    
    /**
     * Returns true if the option takes an argument, and false otherwise.
     */
    virtual bool has_argument() const {
        return Parser::has_argument();
    }
    
    /**
     * Called for no-argument options when the parser encounters them.
     */
    virtual void parse() {
        Parser::parse_default(default_value, value);
    }
    
    /**
     * Called for argument-having options when the parser encounters them. The
     * passed reference is only valid during the function call, so the option
     * should make a copy.
     */
    virtual void parse(const string& arg) {
        Parser::parse(arg, value);
    }
    
protected:
    /// What is the options's long option name
    string long_opt;
    /// What is the option's short option name
    string short_opts;
    /// How is this option described to the user?
    string description;
    /// We keep a value of our chosen type. It just has to be copyable and
    /// assignable.
    Value value;
    /// We also keep the default value around, in case we somehow get assigned
    /// and then get asked to report our metadata.
    /// TODO: do this as a string instead?
    Value default_value;
    

};

/**
 * Represents an option for a type with no extra methods.
 */
template<typename Value, typename Parser = OptionValueParser<Value>>
class Option : public BaseOption<Value, Parser> {
public:
    // Inherit constructor and stuff
    using BaseOption<Value, Parser>::BaseOption;
    
    Option() = delete;
    virtual ~Option() = default;
    
};

/**
 * Specialize and add vector methods.
 * TODO: magically autodetect if a container type has things and expose them.
 * TODO: switch to operator* and operator->.
 */
template<typename Item, typename Parser>
class Option<vector<Item>, Parser> : public BaseOption<vector<Item>, Parser> {
public:
    using Value = vector<Item>;
    
    // Inherit constructor and stuff
    using BaseOption<Value, Parser>::BaseOption;
    
    Option() = delete;
    virtual ~Option() = default;
    
    size_t size() const {
        return this->value.size();
    }
    
    bool empty() const {
        return this->value.empty();
    }
    
    typename Value::value_type& at(size_t i) {
        return this->value.at(i);
    }
    
    const typename Value::value_type& at(size_t i) const {
        return this->value.at(i);
    }
    
    typename Value::iterator begin() {
        return this->value.begin();
    }
    
    typename Value::iterator end() {
        return this->value.end();
    }
    
    typename Value::const_iterator begin() const {
        return this->value.begin();
    }
    
    typename Value::const_iterator end() const {
        return this->value.end();
    }
    
};



}

#endif
