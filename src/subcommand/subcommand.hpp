#ifndef VG_SUBCOMMAND_SUBCOMMAND_H
#define VG_SUBCOMMAND_SUBCOMMAND_H

/**
 * subcommand.hpp: defines a system for registering subcommands of the vg
 * command (vg construct, vg view, etc.) at compile time. Replaces the system of
 * defining two functions and a giant run of if statements in main.cpp.
 *
 * main.cpp does *not* need to include any subcommand headers!
 *
 * Subcommands are created as static global objects in their own compilation
 * units, which have to be explicitly linked into the binary (they won't be
 * pulled out of a library if nothing references their symbols).
 *
 * Subcommands are responsible for printing their own help; we can do "vg help"
 * and print all the subcommands that exist (via a help subcommand), but we
 * can't do "vg help subcommand" (because the help subcommand doesn't know how
 * to get help info on the others).
 *
 * Subcommands get passed all of argv, so they have to skip past their names
 * when parsing arguments.
 *
 * To make a subcommand, do something like this in a cpp file in this
 * "subcommand" directory:
 * 
 * #include "subcommand.hpp"
 * using namespace vg::subcommand;
 * 
 * int main_frobnicate(int argc, char** argv) {
 *     return 0;
 * }
 *
 * static Subcommand vg_frobnicate("frobnicate", "frobnicate nodes and edges",
 *     main_frobnicate);
 */
 
#include <map>
#include <functional>
#include <string>

namespace vg {
namespace subcommand {

/**
 * Represents a subcommand with a name, a description, and some functions.
 * Registers itself on construction in a static registry, and provides static
 * functions for enumerating through that registry.
 */
class Subcommand {

public:
    /**
     * Make and register a subcommand with the given name and description, which
     * calls the given main function when invoked.
     */
    Subcommand(std::string name, std::string description,
        std::function<int(int, char**)> main_function);

    /**
     * Get the name of a subcommand.
     */
    const std::string& get_name() const;
    
    /**
     * Get the description of a subcommand.
     */
    const std::string& get_description() const;
    
    /**
     * Run the main function of a subcommand. Return the return code.
     */
    const int operator()(int argc, char** argv) const;
    
    /**
     * Get the appropriate subcommand to handle the given arguments, or nullptr
     * if no matching subcommand is found.
     */
    static const Subcommand* get(int argc, char** argv);
    
    /**
     * Call the given lambda with each known subcommand, in order.
     */
    static void for_each(const std::function<void(const Subcommand&)>& lambda);


private:
    /**
     * Since we can't rely on a static member field being constructed before any
     * static code that creates actual subcommands gets run, we rely on keeping
     * the registry in a static variable inside a static method, so it gets
     * constructed on first use. Note that at shutdown some of the poinbters in
     * the registry may be to already-destructed static objects.
     */
    static std::map<std::string, Subcommand*>& get_registry();
    
    // These hold the actual fields defining the subcommand
    std::string name;
    std::string description;
    std::function<int(int, char**)> main_function;
    
    /**
     * Get the main function of a subcommand.
     */
    const std::function<int(int, char**)>& get_main() const;
};

}
}

#endif
