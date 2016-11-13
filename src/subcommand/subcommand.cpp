// subcommand.cpp: subcommand registry system implementation

#include "subcommand.hpp"

namespace vg {
namespace subcommand {

Subcommand::Subcommand(std::string name, std::string description,
    std::function<int(int, char**)> main_function) : name(name),
    description(description), main_function(main_function) {
    
    // Add this subcommand to the registry
    Subcommand::get_registry()[name] = this;
}

const std::string& Subcommand::get_name() const {
    return name;
}

const std::string& Subcommand::get_description() const {
    return description;
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

std::map<std::string, Subcommand*>& Subcommand::get_registry() {
    // We keep a static local, which gets initialized when we get called.
    static std::map<std::string, Subcommand*> registry;
    
    // Return a reference to it
    return registry;
}

}
}
