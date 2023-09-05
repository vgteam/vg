// subcommand.cpp: subcommand registry system implementation

#include "subcommand.hpp"

#include <algorithm>
#include <utility>
#include <vector>
#include <limits>

namespace vg {
namespace subcommand {

std::ostream& operator<<(std::ostream& out, const CommandCategory& category) {
    switch(category) {
    case PIPELINE:
        out << "main mapping and calling pipeline";
        break;
    case TOOLKIT:
        out << "useful graph tools";
        break;
    case WIDGET:
        out << "specialized graph tools";
        break;
    case DEVELOPMENT:
        out << "developer commands";
        break;
    case DEPRECATED:
        // we don't show these
        break;
    }
    
    return out;
}

Subcommand::Subcommand(std::string name, std::string description,
    CommandCategory category, int priority,
    std::function<int(int, char**)> main_function) : name(name),
    category(category), priority(priority), description(description),
    main_function(main_function) {
    
    // Add this subcommand to the registry
    Subcommand::get_registry()[name] = this;
}

Subcommand::Subcommand(std::string name, std::string description,
    CommandCategory category, 
    std::function<int(int, char**)> main_function) : Subcommand(name,
    description, category, std::numeric_limits<int>::max(), main_function) {
    
    // Nothing to do!
}

Subcommand::Subcommand(std::string name, std::string description,
    std::function<int(int, char**)> main_function) : Subcommand(name, description, WIDGET, main_function) {
    // Nothing to do!
}

const std::string& Subcommand::get_name() const {
    return name;
}

const std::string& Subcommand::get_description() const {
    return description;
}

const CommandCategory& Subcommand::get_category() const {
    return category;
}

const int& Subcommand::get_priority() const {
    return priority;
}

const int Subcommand::operator()(int argc, char** argv) const {
    return main_function(argc, argv);
}

const Subcommand* Subcommand::get(int argc, char** argv) {
    if(argc < 2) {
        // We don't have a subcommand name
        return nullptr;
    }
    
    if(Subcommand::get_registry().count(argv[1])) {
        // We have a matching subcommand pointer, so return it.
        return Subcommand::get_registry()[argv[1]];
    } else {
        // No matching subcommand was found
        return nullptr;
    }
}

void Subcommand::for_each(const std::function<void(const Subcommand&)>& lambda) {
    for(auto& kv : Subcommand::get_registry()) {
        // For every subcommand, call the callback
        lambda(*kv.second);
    }
}

void Subcommand::for_each(CommandCategory category, const std::function<void(const Subcommand&)>& lambda) {
    if (category == PIPELINE) {
        // Pipeline commands get a special priority order
        
        // We will store them with their priorities and sort them.
        // Easier than writing a custom comparator.
        std::vector<std::pair<int, const Subcommand*>> by_priority;
        
        for_each([&](const Subcommand& command) {
            // Loop over all the subcommands
            if (command.category == category) {
                // And add the ones we care about by priority
                by_priority.push_back(std::make_pair(command.priority, &command));
            }
        });
        
        std::sort(by_priority.begin(), by_priority.end());
        
        for (auto& kv : by_priority) {
            // Now in order of decreasing priority
            // Run the lambda
            lambda(*kv.second);
        }
        
    } else {
        // All other categories just list in alphabetical order
        for_each([&](const Subcommand& command) {
            // Loop over all the subcommands
            if (command.category == category) {
                // And subset to the ones we want
                lambda(command);
            }
        });
    }
}

std::map<std::string, Subcommand*>& Subcommand::get_registry() {
    // We keep a static local, which gets initialized when we get called.
    static std::map<std::string, Subcommand*> registry;
    
    // Return a reference to it
    return registry;
}

}
}
