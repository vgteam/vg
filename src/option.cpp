#include "option.hpp"

/// \file option.cpp: Definitions for our thing-doer-class-attached option
/// system

namespace vg {

using namespace std;

void Configurable::register_option(BaseOption* option) {
    // Turn into an offset from this
    ptrdiff_t offset = (ptrdiff_t) option - (ptrdiff_t) this;
    option_offsets.push_back(offset);
}

vector<BaseOption*> Configurable::get_options() {
    vector<BaseOption*> reconstructed;
    
    for (ptrdiff_t offset : option_offsets) {
        // Convert each offset back to a pointer
        BaseOption* option = (BaseOption*)((ptrdiff_t) this + offset);
        reconstructed.push_back(option);
    }
    
    return reconstructed;
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

    for (BaseOption* option : configurable->get_options()) {
        // For each option...
        
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
        out << "configurable component options:" << endl;
        
        for (BaseOption* option : component->get_options()) {
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
            
            // Then some spaces.
            // TODO: how many?
            for (size_t i = 0; i < 4; i++) {
                out << " ";
            }
            
            // Then the description
            out << option->get_description();
            
            // Finish up with the default value
            out << " [" << option->get_default_value() << "]" << endl;
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
            BaseOption* option = options_by_code.at(c);
            
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
