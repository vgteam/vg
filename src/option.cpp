#include <typeinfo>
#include <cxxabi.h>
#include <cassert>

#include "option.hpp"

/// \file option.cpp: Definitions for our thing-doer-class-attached option
/// system

namespace vg {

using namespace std;

void Configurable::register_option(OptionInterface* option) {
    // Turn into an offset from this
    ptrdiff_t offset = (ptrdiff_t) option - (ptrdiff_t) this;
    option_offsets.push_back(offset);
}

vector<OptionInterface*> Configurable::get_options() {
    vector<OptionInterface*> reconstructed;
    
    for (ptrdiff_t offset : option_offsets) {
        // Convert each offset back to a pointer
        OptionInterface* option = (OptionInterface*)((ptrdiff_t) this + offset);
        reconstructed.push_back(option);
    }
    
    return reconstructed;
}

string Configurable::get_name() {
    // We'll fill this in
    string to_return;

    // Grab the mangled name
    // TODO: On some compilers this isn;t actually mangled.
    auto type_name = typeid(*this).name();
    
    int demangle_error = 0;
    char* demangled = abi::__cxa_demangle(type_name, 0, 0, &demangle_error);
    
    if (demangle_error != 0) {
        // Demangling failed
        if (demangled != nullptr) {
            free(demangled);
        }
        // Use just the original mangled name
        to_return = type_name;
    } else {
        // Demangling succeeded
        assert(demangled != nullptr);
        // Copy name to a string
        string demangled_str(demangled);
        free(demangled);
        // Use the demangled name
        to_return = demangled_str;
    }
    
    if (to_return.size() >= 4 && to_return.substr(0, 4) == "vg::") {
        // Drop any leading vg
        to_return = to_return.substr(4);
    }
    
    return to_return;
}

ConfigurableParser::ConfigurableParser(const char* base_short_options, const struct option* base_long_options,
    function<void(int)> handle_base_option): handle_base_option(handle_base_option) {
    
    // Initialize the set of unused option characters.
    for (char i = 'a'; i <= 'z'; i++) {
        available_short_options.insert(i);
    }
    for (char i = 'A'; i <= 'Z'; i++) {
        available_short_options.insert(i);
    }
    for (char i = '0'; i <= '9'; i++) {
        available_short_options.insert(i);
    }
    
    if (base_short_options != nullptr) {
        // Use these base short options
        short_options = string(base_short_options);
        
        for (auto c : short_options) {
            // Look at each base short option
            if (c == ':') {
                // Ignore the colons marking things as needing arguments
                continue;
            }
            
            if (available_short_options.count(c)) {
                // We would want to assign this, but it's used already.
                available_short_options.erase(c);
            }
        }
    }
    
    if (base_long_options != nullptr) {
        // Use these base long options
        const struct option* long_option = base_long_options;
        while (long_option->name != nullptr || long_option->has_arg != 0 ||
            long_option->flag != nullptr || long_option->val != 0) {
            // This long option isn't the null terminator, so keep it
            long_options.push_back(*long_option);
            
            if (long_option->flag == nullptr && available_short_options.count(long_option->val)) {
                // This long option uses a short option character we might want
                // to assign to something later. Don't re-assign it.
                available_short_options.erase(long_option->val);
            }
            
            // Check the next long option
            long_option++;
        }
    }
}

void ConfigurableParser::register_configurable(Configurable* configurable) {

    // Remember that we have options for this thing.
    configurables.push_back(configurable);

    for (OptionInterface* option : configurable->get_options()) {
        // For each option...
        
        if (long_options_used.count(option->get_long_option())) {
            // TODO: we should try and reassign long options or something if
            // there are collisions.
            throw runtime_error("Duplicate long option: --" +
                option->get_long_option());
        }
        long_options_used.insert(option->get_long_option());
        
        // Try and assign it a short option character.
        char assigned = 0;
        for (char wanted : option->get_short_options()) {
            // First try the ones it wants
            if (available_short_options.count(wanted)) {
                assigned = wanted;
                break;
            }
        }
        if (assigned == 0 && !available_short_options.empty()) {
            // Then fall back to any available
            assigned = (char) *available_short_options.begin();
        }
        if (assigned == 0) {
            // Die because we ran out of room
            throw runtime_error("Unable to assign short option to " + option->get_long_option());
        }
        // Mark the character as used
        available_short_options.erase(assigned);
        // Remember the option for the code
        options_by_code[assigned] = option;
        // And the code for the option
        codes_by_option[option] = assigned;
        
        // Add it to the short options config string
        short_options.push_back(assigned);
        if (option->has_argument()) {
            short_options.push_back(':');
        }
        
        // Add it to the long options struct list
        long_options.push_back({
            option->get_long_option().c_str(),
            option->has_argument() ? required_argument : no_argument,
            nullptr,
            assigned
        });
    }
}

void ConfigurableParser::print_help(ostream& out) const {
    // Print help for all our registered options.
    
    for (Configurable* component : configurables) {
        // For each thing that has options
        
        // Give it a header
        out << component->get_name() << " options:" << endl;
        
        // Determine the padding we need to use to get the option descriptions
        // to line up.
        size_t padding_length = 0;
        for (OptionInterface* option : component->get_options()) {
            // Count the option text
            size_t option_text_length = option->get_long_option().size();
            if (option->has_argument()) {
                // And any argument placeholder
                option_text_length += 4;
            }
            
            // We need to pad to the longest option's length
            padding_length = max(padding_length, option_text_length);
        }
        // Add a couple spaces for visual separation, and so the longest option
        // doesn't run into its description.
        padding_length += 2;
        
        for (OptionInterface* option : component->get_options()) {
            // For each option, find its code.
            // TODO: we assume each option actually gets assigned a code.
            auto& code = codes_by_option.at(option);
        
            // First indent
            out << "    "; 
            
            // Then the short option
            out << "-" << (char) code;
            
            // Then the separator
            out << ", ";
            
            // Then the long option
            out << "--" << option->get_long_option();
            
            if (option->has_argument()) {
                // Then the argument placeholder if needed
                out << " ARG";
            }
            
            // Calculate how long this option's text bit is
            size_t option_text_length = option->get_long_option().size();
            if (option->has_argument()) {
                // And any argument placeholder
                option_text_length += 4;
            }
            
            // Then some spaces to pad to the required padding length
            for (size_t i = 0; i < (padding_length - option_text_length); i++) {
                out << " ";
            }
            
            // Then the description
            out << option->get_description();
            
            // Finish up with the default value, if the option has an argument.
            // Otherwise it must just be a flag; you can't specify an argument.
            if (option->has_argument()) {
                out << " [" << option->get_default_value() << "]";
            }
            
            out << endl;
        }
        
    }
    
    
}

void ConfigurableParser::parse(int argc, char** argv) {
    
    // Stick the null terminator on the long options
    long_options.push_back({0, 0, 0, 0});
    
    while (true) {
        // For each option we get, save its index
        int option_index = 0;
        // And its code
        int c = getopt_long (argc, argv, short_options.c_str(), &long_options[0], &option_index);

        // A code of -1 indicates there are no more options.
        if (c == -1) {
            break;
        }

        // Otherwise, see if the code is associated with an option object.
        if (options_by_code.count(c)) {
            // If so, dispatch to that option
            OptionInterface* option = options_by_code.at(c);
            
            if (option->has_argument()) {
                // Parse with the option's argument
                option->parse(optarg);
            } else {
                // Parse without an argument (i.e. set a bool)
                option->parse();
            }
        } else {
            // This option is to be handled manually, or is unrecognized
            handle_base_option(c);
        }
    }

    // Pop the null terminator again
    long_options.pop_back();
}








    
}
