#ifndef VG_SUBCOMMAND_OPTIONS_HPP_INCLUDED
#define VG_SUBCOMMAND_OPTIONS_HPP_INCLUDED

/**
 *\file
 * options.hpp: option parser system.
 *
 * Make a BaseOptionGroup, and use add_group<Receiver>(heading) to add
 * subgroups that apply parsed options to instances of a given class.
 *
 * Set up each option in the group with add_option(), add_range() (for an
 * option that can be cycled through a range of values for a grid search), or
 * add_flag(). Every option always has a logn option name; short option
 * character is optional and comes after it. Options take a pointer-to-member
 * into the group's type (where the value will be ultimately written) and a
 * default value, a help stiring, and an optional "validator function" which
 * gets to check the parsed value and should raise std::domain_error with a
 * complaint if the value isn't acceptable.
 *
 * Get the help with get_help() on the root group, and use the free
 * print_table() function to print it to a stream with headings and automatic
 * indentation.
 *
 * To parse options, use make_long_options() and make_short_options() to adjust
 * some (possibly already partly populated) getopt_long() inputs,
 * null-terminate the long options array, and call getopt_long() as normal.
 * Show every option number and optarg value to the root group's parse(), and
 * then continue with your own parsing if it returns false because it does not
 * consume that option.
 *
 * To apply presets before parsing, make a Preset and fill it in with
 * add_entry<T>(option, value). Then call apply() on the preset with the root
 * group.
 *
 * To read and write option values manually, use the get_option_value() and
 * set_option_value() methods on the root option group.
 *
 * To log option values, there is a print_options() method on the root option
 * group.
 *
 * EXAMPLE
 *
 * struct ThingDoer {
 *     static constexpr int default_count = 5;
 *     int count = default_count;
 * };
 *
 * vg::subcommand::BaseOptionGroup get_parser() {
 *     vg::subcommand::BaseOptionGroup parser;
 *     auto thing_doer_opts = parser.add_group<ThingDoer>("thing doer configuration");
 *     thing_doer_opts.add_option(
 *         "count", 'c',
 *         &ThingDoer::count,
 *         ThingDoer::default_count,
 *         "number of things to do",
 *         vg::subcommand::int_is_nonnegative
 *     );
 *     return parser;
 * }
 *
 * int main(int argc, char** argv) {
 *     auto parser = get_parser();
 *     std::vector<struct option> long_options;
 *     parser.make_long_options(long_options);
 *     long_options.push_back({0, 0, 0, 0});
 *     std::string short_options;
 *     parser.make_short_optins(short_options);
 *     int c;
 *     while (true) {
 *         int option_index = 0;
 *         int option_id = getopt_long (argc, argv, short_options.c_str(),
 *                                      &long_options[0], &option_index);
 *         if (option_id == -1) break;
 *         if (parser.parse(option_id, optarg)) continue;
 *         switch (option_id) {
 *         default:
 *             vg::subcommand::print_table(parser.get_help(), std::cerr);
 *             return 1;
 *         }
 *     }
 *     ThingDoer thing_doer;
 *     parser.apply(thing_doer);
 *     std::cout << "Doing " << thing_doer.count << " things" << std::endl;
 *     return 0;
 * }
 */

#include "../utility.hpp"

#include <map>
#include <functional>
#include <string>
#include <iostream>

#include <getopt.h>

namespace vg {
namespace subcommand {

/**
 * Interface for things that form a "chain" that can be "ticked".
 *
 * Each link in the chain works like a digit place in a number, and ticking increments the number.
 * This lets us do gird search over a bunch of values of different types without a bunch of nexted loops.
 */
struct TickChainLink {
    /// This will be called when we want to reset_chain what we are chained onto.
    std::function<void(void)> reset_chain_parent = []() {
    };
    /// This will be called when we need to tick_chain our parent
    std::function<bool(void)> tick_chain_parent = []() {
        return false;
    };
    
    /// Reset the chain to its initial values.
    virtual void reset_chain();
    
    /// Tick the chain. Return true if there's still a value for the chain, and
    /// false if the chain is out of values.
    virtual bool tick_chain();
    
    /// Add a thing to the chain after us.
    /// Return that thing.
    virtual TickChainLink& chain(TickChainLink& next);
    
    /// Get a function that runs another function for each combination of
    /// values for this Range and all Ranges it has been chained onto.
    virtual std::function<void(const std::function<void(void)>&)> get_iterator();
};

}
}

namespace vg {
// TODO: If Range isn't in the vg namespace, then vg::parse<>'s specialization
// for it doesn't actually get treated as a specialization, and we don't
// instantiate the template when we need to, and we get linker errors because
// the linker can't find the instantiation of the template. Someday we should
// work out why and fix this.

/**
 * Tickable link that represents a single value or a range of values.
 * Range rusn from start to <=end, going up by step.
 * You can set the range to s aingle value or to a full range, and when you read it you see the current value.
 */
template<typename Number>
struct Range : public subcommand::TickChainLink {
    
    // Expose the thing we are a range of
    using type = Number;

    /// Represents the start of the range
    Number start = 0;
    /// Represents the inclusive end of the range
    Number end = 0;
    /// Represents the step to move by each tick
    Number step = 1;
    
    /// Represents the current value the range is at
    Number here = 0;
    /// Determines if we are running or not (i.e. is here valid)
    bool running = false;
    
    /// Default constructor
    Range() {
        // Nothing to do!
    }
    
    /// Construct from a single value
    Range(const Number& val): start(val), end(val) {
        // Nothing to do!
    }
    
    /// Copy, preserving destination links
    Range(const Range& other): start(other.start), end(other.end), step(other.step) {
        // Nothing to do
    }
    
    /// Move, preserving destination links
    Range(Range&& other): start(other.start), end(other.end), step(other.step) {
        // Nothing to do
    }
    
    /// Copy assignment, preserving destination links
    Range& operator=(const Range& other) {
        start = other.start;
        end = other.end;
        step = other.step;
        return *this;
    }
    
    /// Move assignment, preserving destination links
    Range& operator=(Range&& other) {
        start = other.start;
        end = other.end;
        step = other.step;
        return *this;
    }
    
    /// Check the range for usefulness
    inline bool is_valid() {
        if (start != end && step == 0) {
            // We'll never make it
            cerr << "Invalid range (no movement): " << start << " to " << end << " step " << step << endl;
            return false;
        }
        
        if (start > end && step > 0) {
            // We're going the wrong way
            cerr << "Invalid range (need to go down): " << start << " to " << end << " step " << step << endl;
            return false;
        }
        
        if (start < end && step < 0) {
            // We're going the other wrong way
            cerr << "Invalid range (need to go up): " << start << " to " << end << " step " << step << endl;
            return false;
        }
        
        return true;
    }
    
    /// Convert to Number with the current value
    operator Number() const {
        if (running) {
            return here;
        } else {
            return start;
        }
    }
    
    /// Start at our start value
    void reset() {
        here = start;
        running = true;
    }
    
    /// Start us and all the things we are chained onto at their start values
    void reset_chain() {
        reset();
        reset_chain_parent();
    }
    
    /// Increment our value.
    /// Returns true if the new value needs processing, and false if we have left or would leave the range.
    bool tick() {
        if (here == end) {
            // We are at the end
            return false;
        }
        
        here += step;
        if ((step > 0 && here > end) || (step < 0 && here < end)) {
            // We have passed the end (for things like double)
            return false;
        }
        
        return true;
    }
    
    /// Increment our value.
    /// If it overflows, tick_chain whatever we are chained onto, and reset and succeed if that succeeds.
    bool tick_chain() {
        if (tick()) {
            // We could change
            return true;
        } else {
            // We couldn't change.
            if (tick_chain_parent()) {
                // We have a parent we could advance.
                reset();
                return true;
            } else {
                // Our parent couldn't advance either.
                return false;
            }
        }
    }
};

}

namespace vg {

// Define a way to test if a type is an instantiation of a template on a type
// See https://stackoverflow.com/a/25803794

// In general, things aren't instantiations of things
template <typename Subject, template<typename...> class Predicate>
struct is_instantiation_of : std::false_type {
};

// Except things that are instantiations of things with some arguments
template <template<typename... > class Predicate, class... PredicateArgs>
struct is_instantiation_of<Predicate<PredicateArgs...>, Predicate> : std::true_type {
};

/// Parse a range as start[:end[:step]]
template<typename Result>
inline bool parse(const string& arg, typename enable_if<is_instantiation_of<Result, Range>::value, Result>::type& dest) {

    auto colon1 = arg.find(':');
    
    if (colon1 == string::npos) {
        // No colons here. Parse one number.
        if (!parse<typename Result::type>(arg, dest.start)) {
            return false;
        }
        dest.end = dest.start;
        dest.step = 0;
        return dest.is_valid();
    } else if (colon1 == arg.size()) {
        // Can't end in a colon
        return false;
    } else {
        // Look for another colon
        auto colon2 = arg.find(':', colon1 + 1);
        if (colon2 == string::npos) {
            // Just a range of two things
            if (!parse<typename Result::type>(arg.substr(0, colon1), dest.start)) {
                return false;
            }
            if (!parse<typename Result::type>(arg.substr(colon1 + 1), dest.end)) {
                return false;
            }
            dest.step = 1;
            return dest.is_valid();
        } else if (colon2 == arg.size()) {
            // Can't end in a colon
            return false;
        } else {
            // We have 3 numbers
            if (!parse<typename Result::type>(arg.substr(0, colon1), dest.start)) {
                return false;
            }
            if (!parse<typename Result::type>(arg.substr(colon1 + 1, colon2 - colon1 - 1), dest.end)) {
                return false;
            }
            if (!parse<typename Result::type>(arg.substr(colon2 + 1), dest.step)) {
                return false;
            }
            
            return dest.is_valid();
        }
    }
}

}

namespace vg {
namespace subcommand {

/// Get a new unique option ID.
int get_option_id();

/**
 * Get a string "metavar" placeholder for a command line option, appropriate to its type.
 */
template<typename T>
const char* get_metavar();

template<>
const char* get_metavar<size_t>();

template<>
const char* get_metavar<int>();

template<>
const char* get_metavar<bool>();

template<>
const char* get_metavar<double>();

template<>
const char* get_metavar<std::string>();

/**
 * Represents an option being set to a value. Base interface.
 */
struct BaseValuation {
    /// Make a new BaseValuation for the given option
    BaseValuation(const std::string& option);
    virtual ~BaseValuation() = default;
    
    /// Long option to give a value to
    std::string option;
};

/**
 * Represents an option being set to a value. Actually has the value.
 */
template<typename T>
struct Valuation : public BaseValuation {
    /// Make a preset entry that sets the given long option to the given value.
    Valuation(const std::string& option, const T& value) : BaseValuation(option), value(value) {
        // Nothing to do
    }
    
    virtual ~Valuation() = default;

    /// Value for the option
    T value;
};

/// Function type used to validate arguments. Throw std::domain_error if not allowed, explaining why.
template<typename T>
using ValidatorFunction = std::function<void(const T&)>;

/// Validate that a double is positive, or throw std::domain_error
extern const ValidatorFunction<double> double_is_positive;

/// Validate that a double is not negative, or throw std::domain_error
extern const ValidatorFunction<double> double_is_nonnegative;

/// Validate that a size_t is not zero, or throw std::domain_error
extern const ValidatorFunction<size_t> size_t_is_nonzero;

/// Validate that an int is not negative, or throw std::domain_error;
extern const ValidatorFunction<int> int_is_nonnegative;

/**
 * Interface for a command-line argument that goes into a field on an object of
 * the given type.
 */
template<typename Receiver>
struct BaseArgSpec : public TickChainLink {
    /// Make an option with a long and short option name
    BaseArgSpec(const std::string& option, char short_option, const std::string& help) : option(option), help(help), short_option(short_option), option_id(short_option != '\0' ? short_option : get_option_id()) {
        // Nothing to do
    }
    /// Make an option with a long option name only
    BaseArgSpec(const std::string& option, const std::string& help) : BaseArgSpec(option, '\0', help) {
        // Nothing to do
    }
    virtual ~BaseArgSpec() = default;
    
    /// Parse the argument's value from the command line.
    /// Throws std::domain_error if validation fails.
    virtual void parse(const char* optarg) = 0;
    /// Apply a preset item, or fail if it doesn't match.
    /// The preset value will sit under any parsed value but above the default.
    virtual void preset(const BaseValuation& entry) = 0;
    /// Apply a valuation, or fail if it doesn't match.
    /// The value will replace any parsed value!
    /// Validation will not be run!
    virtual void set(const BaseValuation& entry) = 0;
    /// Put our current effective value into the given BaseValuation, which
    /// must be for the right option and have the right type.
    virtual void query(BaseValuation& entry) const = 0;
    /// Apply the value to the right field of the given object.
    virtual void apply(Receiver& receiver) const = 0;
     /// Print value to the given stream after the given separator.
    virtual void print_value(ostream& out, const char* sep = "") const = 0;
    /// Print value metavar placeholder to the given stream after the given separator.
    virtual void print_metavar(ostream& out, const char* sep = "") const = 0;
    /// Print default value to the given stream, if appropriate.
    virtual void print_default(ostream& out) const = 0;
    /// Print option and value to the given stream, without newlines, between the given separators.
    /// If slug is set, use short option if available and don't include spaces.
    virtual void print(ostream& out, const char* sep = "", const char* after = "", bool slug = false) const {
        out << sep;
        if (slug && short_option != '\0') {
            out << "-" << short_option;
        } else {
            out << "--" << option;
        }
        this->print_value(out, slug ? "" : " ");
        out << after;
    }
    /// Get the getopt structure for this option. Option must outlive it and not move.
    virtual struct option get_option_struct() const = 0;
    
    /// Name of the option (long opt)
    std::string option;
    /// Help for the option
    std::string help;
    /// Character of the option (short opt), or 0
    char short_option;
    /// Int value to represent the option
    int option_id;
};

/**
 * Interface for a command-line argument that corresponds to a value of a given type.
 * Storage method is left to be implemented by inheritor.
 */
template<typename T, typename Receiver>
struct ArgSpec : public BaseArgSpec<Receiver> {
    /// Make an option with a long and short option name
    ArgSpec(const std::string& option, char short_option, T Receiver::*dest, const T& default_value, const std::string& help,  const ValidatorFunction<T>& validator) : BaseArgSpec<Receiver>(option, short_option, help), dest(dest), default_value(default_value), validator(validator) {
        // Nothing to do!
    }
    /// Make an option with a long option name only
    ArgSpec(const std::string& option, T Receiver::*dest, const T& default_value, const std::string& help, const ValidatorFunction<T>& validator) : ArgSpec(option, '\0', dest, default_value, help, validator) {
        // Nothing to do!
    }
    
    virtual ~ArgSpec() = default;
    
    /// Allow setting our stored value
    virtual void set_value(const T& value) = 0;
    /// And getting our current effective value
    virtual T get_value() const = 0;
    /// Return true if a value has been set from parsing or a preset.
    virtual bool was_set() const = 0;
    
    virtual void preset(const BaseValuation& entry) {
        // Needs to be a preset for the right option
        assert(entry.option == this->option);
        const Valuation<T>* as_typed = dynamic_cast<const Valuation<T>*>(&entry);
        if (as_typed) {
            if (!this->was_set()) {
                // Apply the preset value, if nothing is set yet.
                this->set_value(as_typed->value);
            }
        } else {
            throw std::runtime_error("Could not cast valuation");
        }
    }
    
    virtual void set(const BaseValuation& entry) {
        // Needs to be for the right option
        assert(entry.option == this->option);
        const Valuation<T>* as_typed = dynamic_cast<const Valuation<T>*>(&entry);
        if (as_typed) {
            // Apply the value
            this->set_value(as_typed->value);
        } else {
            throw std::runtime_error("Could not cast valuation");
        }
    }
    
    virtual void query(BaseValuation& entry) const {
        // Needs to be a valuation for the right option
        assert(entry.option == this->option);
        Valuation<T>* as_typed = dynamic_cast<Valuation<T>*>(&entry);
        if (as_typed) {
            // Put our value in there.
            as_typed->value = this->get_value();
        } else {
            throw std::runtime_error("Could not cast valuation");
        }
    }
    
    /// Field in the receiving type we set.
    T Receiver::*dest;
    /// Original default value.
    T default_value;
    /// Function to check value with
    ValidatorFunction<T> validator;
};

/**
 * Definition structure for normal value-having options. Lets you specify
 * storage type for the actual value.
 */
template<typename T, typename Receiver, typename Holder=T>
struct ValueArgSpec : public ArgSpec<T, Receiver> {
    /// Make an option with a long and short option name
    ValueArgSpec(const std::string& option, char short_option, T Receiver::*dest, const T& default_value, const std::string& help, const ValidatorFunction<T>& validator) : ArgSpec<T, Receiver>(option, short_option, dest, default_value, help, validator), value(default_value)  {
        // Nothing to do
    }
    /// Make an option with a long option name only
    ValueArgSpec(const std::string& option, T Receiver::*dest, const T& default_value, const std::string& help, const ValidatorFunction<T>& validator) : ValueArgSpec(option, '\0', dest, default_value, help, validator) {
        // Nothing to do
    }
    virtual ~ValueArgSpec() = default;
    
    virtual void set_value(const T& replacement) {
        // We assume the holder supports assignment.
        this->value = replacement;
        // Remember we got a value applied. Presets shouldn't clobber it.
        this->value_set = true;
    }
    
    virtual T get_value() const {
        return this->value;
    }
    
    virtual bool was_set() const {
        return this->value_set;
    }
    
    virtual void parse(const char* optarg) {
        try {
            if (!optarg) {
                // Protect against nulls
                throw std::domain_error("requires a value");
            }
       
            this->value = vg::parse<Holder>(optarg);
            this->validator(this->value);
            this->value_set = true;
        } catch (std::domain_error& e) {
            cerr << "error: option ";
                if (this->short_option) {
                    cerr << "-" << this->short_option << "/";
                }
                cerr << "--" << this->option << " ";
                cerr << e.what() << endl;
                exit(1);
        }
    }
    
    virtual void apply(Receiver& receiver) const {
        receiver.*(this->dest) = value;
    }
    virtual void print_metavar(ostream& out, const char* sep = "") const {
        out << sep << get_metavar<T>();
    }
    virtual void print_value(ostream& out, const char* sep = "") const {
        out << sep << value;
    }
    virtual void print_default(ostream& out) const {
        out << " [" << this->default_value << "]";
    }
    virtual struct option get_option_struct() const {
        return {this->option.c_str(), required_argument, 0, this->option_id};
    }
    
    
    Holder value;
    bool value_set = false;
};

/**
 * Definition structure for value-having options that can run through a range.
 */
template<typename T, typename Receiver>
struct RangeArgSpec : public ValueArgSpec<T, Receiver, Range<T>> {
    using Holder = Range<T>;
    
    using ValueArgSpec<T, Receiver, Range<T>>::ValueArgSpec;
    virtual ~RangeArgSpec() = default;
    
    virtual TickChainLink& chain(TickChainLink& next) {
        // Wire our value range into the chain.
        TickChainLink::chain(this->value);
        this->value.chain(next);
        return next;
    }
};

/**
 * Definition structure for flag options that flip a default value.
 */
template<typename Receiver>
struct FlagArgSpec : public ValueArgSpec<bool, Receiver> {
    using T = bool;
    using Holder = T;
    
    using ValueArgSpec<bool, Receiver>::ValueArgSpec;
    virtual ~FlagArgSpec() = default;
    
    virtual void parse(const char* optarg) {
        // When parsing, flip stored default.
        this->set_value(!this->default_value);
    }
    virtual void print_metavar(ostream& out, const char* sep = "") const {
        // Don't do anything
    }
    virtual void print_value(ostream& out, const char* sep = "") const {
        // Don't do anything
    }
    virtual void print_default(ostream& out) const {
        // Don't do anything
    }
    virtual void print(ostream& out, const char* sep = "", const char* after = "", bool slug = false) const {
        // Override print to just print the flag when used
        if (this->value != this->default_value) {
            out << sep;
            if (slug && this->short_option != '\0') {
                out << "-" << this->short_option;
            } else {
                out << "--" << this->option;
            }
            out << after;
        }
    }
    virtual struct option get_option_struct() const {
        return {this->option.c_str(), no_argument, 0, this->option_id};
    }
};

/**
 * Represents a set of command-line options.
 */
struct BaseOptionGroup : public TickChainLink {

    virtual ~BaseOptionGroup() = default;
    
    /// Parse the given option ID, with the given value if needed.
    /// Return true if we matched the ID, and false otherwise.
    virtual bool parse(int option_id, const char* optarg) = 0;
    
    /// Apply a preset value to its option. Returns true if it was found, and
    /// false otherwies.
    virtual bool preset(const BaseValuation& entry) = 0;
    
    /// Apply a value to its option. Returns true if it was found, and false
    /// otherwies.
    virtual bool set(const BaseValuation& entry) = 0;
    
    /// Fill in entry with the value of the correspondign option, if we have
    /// that option. If so, return true.
    virtual bool query(BaseValuation& entry) const = 0;
    
    /// Print all options set.
    /// By default, prints one option per line.
    /// If slug is set, prints short options, all on one line.
    virtual void print_options(ostream& out, bool slug = false) const = 0;
    
    /// Get help, in the form of pairs of options and descriptions.
    /// Headings are descriptions without options.
    virtual std::vector<std::pair<std::string, std::string>> get_help() const = 0;
    
    /// Add options to non-null-terminated input for getopt_long
    virtual void make_long_options(std::vector<struct option>& dest) const = 0;
    
    /// Add options to string input for getopt_long
    virtual void make_short_options(std::string& dest) const = 0;
    
    /// Allow the user to query an option value by name.
    /// Would be simpler if we could override template methods but we can't.
    template<typename T>
    T get_option_value(const std::string& option) const {
        Valuation<T> question(option, T());
        bool found = this->query(question);
        if (!found) {
            throw std::runtime_error("Undefined option: " + option);
        }
        return question.value;
    }
    
    /// Allow the user to manually set an option value
    template<typename T>
    void set_option_value(const std::string& option, const T& value) {
        Valuation<T> setter(option, value);
        bool found = this->set(setter);
        if (!found) {
            throw std::runtime_error("Undefined option: " + option);
        }
    }
};

/**
 * Represents a set of command-line options that can be applied to an object.
 * Internal values can be ranges that can be ticked.
 * Comes with a heading.
 */
template<typename Receiver>
struct OptionGroup : public BaseOptionGroup {

    virtual ~OptionGroup() = default;

    /// Make an option group woith the given heading.
    OptionGroup(const std::string& heading) : heading(heading) {
        // Nothing to do!
    }
   
    /// Chain through all options 
    virtual TickChainLink& chain(TickChainLink& next) {
        if (args.empty()) {
            // Just chain through
            return TickChainLink::chain(next);
        } else {
            // Chain us to first arg, and last arg to next.
            TickChainLink::chain(*args.front());
            args.back()->chain(next);
            return next;
        }
    }
    
    // We need to take default_value by value, and not by reference, because we
    // often want to pass stuff that is constexpr and trying to use a reference
    // will make us try to link against it.
    // TODO: C++17 fixes this, so fix when we use that.
    
    /// Add a new option that goes to the given field, with the given default.
    template<typename T, typename Spec = ValueArgSpec<T, Receiver>>
    void add_option(const std::string& name, char short_option, T Receiver::*dest, T default_value, const std::string& help, const ValidatorFunction<T>& validator = [](const T& ignored) {}) {
        args.emplace_back(new Spec(name, short_option, dest, default_value, help, validator));
        if (args.size() > 1) {
            // Chain onto previous option
            args[args.size() - 2]->chain(*args[args.size() - 1]);
        }
        // Index it by option ID
        id_to_index.emplace(args[args.size() - 1]->option_id, args.size() - 1);
        // And option name
        option_to_index.emplace(args[args.size() - 1]->option, args.size() - 1);
    }
    
    /// Add a new option that goes to the given field, with the given default.
    template<typename T, typename Spec = ValueArgSpec<T, Receiver>>
    void add_option(const std::string& name, T Receiver::*dest, T default_value, const std::string& help, const ValidatorFunction<T>& validator = [](const T& ignored) {}) {
        add_option<T, Spec>(name, '\0', dest, default_value, help, validator);
    }
    
    /// Add a new option that handles range values
    template<typename T>
    void add_range(const std::string& name, char short_option, T Receiver::*dest, T default_value, const std::string& help, const ValidatorFunction<T>& validator = [](const T& ignored) {}) {
        add_option<T, RangeArgSpec<T, Receiver>>(name, short_option, dest, default_value, help, validator);
    }
    /// Add a new option that handles range values
    template<typename T>
    void add_range(const std::string& name, T Receiver::*dest, T default_value, const std::string& help, const ValidatorFunction<T>& validator = [](const T& ignored) {}) {
        add_range<T>(name, '\0', dest, default_value, help, validator);
    }
    
    /// Add a new option that is a boolean flag
    void add_flag(const std::string& name, char short_option, bool Receiver::*dest, bool default_value, const std::string& help, const ValidatorFunction<bool>& validator = [](const bool& ignored) {}) {
        add_option<bool, FlagArgSpec<Receiver>>(name, short_option, dest, default_value, help, validator);
    }
    /// Add a new option that is a boolean flag
    void add_flag(const std::string& name, bool Receiver::*dest, bool default_value, const std::string& help, const ValidatorFunction<bool>& validator = [](const bool& ignored) {}) {
        add_flag(name, '\0', dest, default_value, help, validator);
    }
    
    /// Parse the given option ID, with the given value if needed.
    /// Return true if we matched the ID, and false otherwise.
    virtual bool parse(int option_id, const char* optarg) {
        auto found = id_to_index.find(option_id);
        if (found != id_to_index.end()) {
            // We have this option, so parse.
            args.at(found->second)->parse(optarg);
            return true;
        } else {
            // We don't have this option, maybe someone else does.
            return false;
        }
    }
    
    virtual bool preset(const BaseValuation& entry) {
        auto found = option_to_index.find(entry.option);
        if (found != option_to_index.end()) {
            // We have this option, so assign the preset.
            args.at(found->second)->preset(entry);
            return true;
        } else {
            // We don't have this option, maybe someone else does.
            return false;
        }
    }
    
    virtual bool set(const BaseValuation& entry) {
        auto found = option_to_index.find(entry.option);
        if (found != option_to_index.end()) {
            // We have this option, so assign the preset.
            args.at(found->second)->set(entry);
            return true;
        } else {
            // We don't have this option, maybe someone else does.
            return false;
        }
    }
    
    virtual bool query(BaseValuation& entry) const {
        auto found = option_to_index.find(entry.option);
        if (found != option_to_index.end()) {
            // We have this option, so get its value.
            args.at(found->second)->query(entry);
            return true;
        } else {
            // We don't have this option, maybe someone else does.
            return false;
        }
    }
    
    /// Print all options set, one per line
    virtual void print_options(ostream& out, bool slug = false) const {
        if (slug) {
            for (auto& arg : args) {
                // Print unseparated short options
                arg->print(out, "", "", true);
            }
        } else {
            for (auto& arg : args) {
                // Print long options, one per line
                arg->print(out, "", "\n");
            }
        }
    }
    
    /// Apply all flags to the receiver
    void apply(Receiver& receiver) const {
        for (auto& arg : args) {
            arg->apply(receiver);
        }
    }
    
    /// Get help, in the form of pairs of options and descriptions.
    /// Headings are descriptions without options.
    virtual std::vector<std::pair<std::string, std::string>> get_help() const {
        std::vector<std::pair<std::string, std::string>> to_return;
        to_return.reserve(args.size() + 1);
        
        // Put the heading
        to_return.emplace_back("", heading + ":");
        
        for (auto& arg : args) {
            // Show the option
            std::stringstream opt_stream;
            if (arg->short_option) {
                opt_stream << "-" << arg->short_option << ", ";
            }
            opt_stream << "--" << arg->option;
            arg->print_metavar(opt_stream, " ");
            
            // Show the help text with default value
            std::stringstream help_stream;
            help_stream << arg->help;
            arg->print_default(help_stream);
            
            to_return.emplace_back(opt_stream.str(), help_stream.str());
        }
        
        return to_return;
    }
    
    /// Add to non-null-terminated input for getopt_long
    virtual void make_long_options(std::vector<struct option>& dest) const {
        dest.reserve(dest.size() + args.size());
        for (auto& arg : args) {
            // Collect them from all the options
            dest.emplace_back(arg->get_option_struct());
        }
    }
    
    /// Add to string input for getopt_long
    virtual void make_short_options(std::string& dest) const {
        for (auto& arg : args) {
            struct option long_spec = arg->get_option_struct();
            if (long_spec.val < std::numeric_limits<char>::max()) {
                // This has a short option. Encode the short option string.
                dest.push_back(long_spec.val);
                switch (long_spec.has_arg) {
                case optional_argument:
                    dest.push_back(':');
                    // Fall-through
                case required_argument:
                    dest.push_back(':');
                    // Fall-through
                case no_argument:
                    break;
                }
            }
        }
    }
    
    /// Heading we will appear under in the help.
    std::string heading;
    /// Holds the argument definitions and parsing destinations
    std::vector<std::unique_ptr<BaseArgSpec<Receiver>>> args;
    /// Map from option ID to option index
    std::unordered_map<int, size_t> id_to_index;
    /// Map from long option to option index, to allow applying presets.
    std::unordered_map<std::string, size_t> option_to_index;
};

/**
 * Represents a group of groups of options.
 *
 * Also doubles as the main parser type; you can make one of these and populate
 * it with subgroups and options, and then use get_help() and print_table() to
 * do help, and make_long_options(), make_short_options(), and parse() to parse
 * with getopt_long(), and then you can apply() the options to objects they
 * eventually belong in.
 */
struct GroupedOptionGroup : public BaseOptionGroup {
    
    // We can't copy because we contain unique_ptr values
    GroupedOptionGroup() = default;
    GroupedOptionGroup(const GroupedOptionGroup& other) = delete;
    GroupedOptionGroup& operator=(GroupedOptionGroup& other) = delete;
    GroupedOptionGroup(GroupedOptionGroup&& other) = default;
    GroupedOptionGroup& operator=(GroupedOptionGroup&& other) = default;
    virtual ~GroupedOptionGroup() = default;

    /// Create a new child group with a new heading, which we can add options
    /// to.
    template<typename Receiver>
    OptionGroup<Receiver>& add_group(const std::string& heading) {
        OptionGroup<Receiver>* new_group = new OptionGroup<Receiver>(heading);
        subgroups.emplace_back(new_group);
        if (subgroups.size() > 1) {
            // Chain the groups
            subgroups[subgroups.size() - 2]->chain(*subgroups[subgroups.size() - 1]);
        }
        return *new_group;
    }
    
    /// Apply all options that go on an object of this type to the given object.
    template<typename Receiver>
    void apply(Receiver& receiver) {
        for (auto& group : subgroups) {
            OptionGroup<Receiver>* as_relevant_leaf = dynamic_cast<OptionGroup<Receiver>*>(group.get());
            if (as_relevant_leaf) {
                // This is a group that cares about this type
                as_relevant_leaf->apply(receiver);
            }
            GroupedOptionGroup* as_internal_node = dynamic_cast<GroupedOptionGroup*>(group.get());
            if (as_internal_node) {
                // This is a group that has child groups
                as_internal_node->apply(receiver);
            }
        }
    }
    
    /// Chain through all subgroups 
    virtual TickChainLink& chain(TickChainLink& next);
    
    virtual bool parse(int option_id, const char* optarg); 
    
    virtual bool preset(const BaseValuation& entry);
    
    virtual bool set(const BaseValuation& entry);
    
    virtual bool query(BaseValuation& entry) const;
    
    virtual void print_options(ostream& out, bool slug = false) const;
    
    virtual std::vector<std::pair<std::string, std::string>> get_help() const;
    
    virtual void make_long_options(std::vector<struct option>& dest) const;
    
    virtual void make_short_options(std::string& dest) const;
    
    /// Holds all the child groups of options
    std::vector<std::unique_ptr<BaseOptionGroup>> subgroups;
};

/**
 * Represents a named preset of command-line option default overrides. Options
 * are organized by long option name.
 *
 * Make one of these, and use add_entry() to add values to it, and then use
 * apply(root_option_group) to apply it.
 */
struct Preset {
    /// As part of this preset, set the given option to the given value.
    template<typename T>
    Preset& add_entry(const std::string& option, const T& value) {
        Valuation<T>* entry = new Valuation<T>(option, value);
        entries.emplace_back(entry);
        return *this;
    }
    
    /// Apply stored presets to the given parser
    void apply(BaseOptionGroup& parser) const {
        for (auto& entry : entries) {
            // Apply the entry
            bool applied = parser.preset(*entry);
            // Make sure it worked
            assert(applied);
        }
    }
    
    std::vector<std::unique_ptr<BaseValuation>> entries;
};

/**
 * Print a table of rows, with each column starting at the same character on the line.
 *
 * Prints the help from get_help() on an option parsing group in a nice way.
 */
void print_table(const std::vector<std::pair<std::string, std::string>>& rows, ostream& out);


}
}

#endif
