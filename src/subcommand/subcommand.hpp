#ifndef VG_SUBCOMMAND_SUBCOMMAND_HPP_INCLUDED
#define VG_SUBCOMMAND_SUBCOMMAND_HPP_INCLUDED

/** \file
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
 * can't do "vg help subcommand" and have that be equivalent to "vg subcommand
 * --help" (because the help subcommand doesn't know how to get help info on the
 * others).
 *
 * We have a subcommand importance/category system, so we can tell people about
 * just the main pipeline and keep the subcommands they don't want out of their
 * brains and off their screen.
 *
 * Subcommands get passed all of argv, so they have to skip past their names
 * when parsing arguments.
 *
 * To make a subcommand, do something like this in a cpp file in this
 * "subcommand" directory:
 * 
 *     #include "subcommand.hpp"
 *     using namespace vg::subcommand;
 * 
 *     int main_frobnicate(int argc, char** argv) {
 *         return 0;
 *     }
 *
 *     static Subcommand vg_frobnicate("frobnicate", "frobnicate nodes and edges",
 *         main_frobnicate);
 * 
 */
 
#include <map>
#include <functional>
#include <string>
#include <iostream>

namespace vg {
namespace subcommand {

/**
 * Defines what kind of command each subcommand is.
 */
enum CommandCategory {
    /// Some commands are part of the main build-graph, align, call variants pipeline
    PIPELINE, 
    /// Some subcommands are important parts of the toolkit/swiss army knife for working with graphs and data
    TOOLKIT,
    /// Some commands are less important but potentially useful widgets that let you do a thing you might need
    WIDGET,
    /// Some commands are useful really only for developers
    DEVELOPMENT,
    /// Some commands we're trying to move away from
    DEPRECATED
};

/// Define a way to print the titles of the different categories
std::ostream& operator<<(std::ostream& out, const CommandCategory& category);

/**
 * Represents a subcommand with a name, a description, and some functions.
 * Registers itself on construction in a static registry, and provides static
 * functions for enumerating through that registry.
 */
class Subcommand {

public:
    
    /**
     * Make and register a subcommand with the given name and description, in
     * the given category, with the given priority (lower is better), which
     * calls the given main function when invoked.
     */
    Subcommand(std::string name, std::string description,
        CommandCategory category, int priority,
        std::function<int(int, char**)> main_function);
    
    /**
     * Make and register a subcommand with the given name and description, in
     * the given category, with worst priority, which calls the given main
     * function when invoked.
     */
    Subcommand(std::string name, std::string description,
        CommandCategory category,
        std::function<int(int, char**)> main_function);
    
    /**
     * Make and register a subcommand with the given name and description, in
     * the WIDGET category, with worst priority, which calls the given main
     * function when invoked.
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
     * Get the category of a subcommand, which determines who might want to use
     * it and why.
     */
    const CommandCategory& get_category() const;
    
    /**
     * Get the priority level of a subcommand (lower is more important).
     */
    const int& get_priority() const;
    
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
    
    /**
     * Call the given lambda with each known subcommand in the given category,
     * in order.
     */
    static void for_each(CommandCategory category, const std::function<void(const Subcommand&)>& lambda);


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
    CommandCategory category;
    int priority;
    std::function<int(int, char**)> main_function;
    
    /**
     * Get the main function of a subcommand.
     */
    const std::function<int(int, char**)>& get_main() const;
};

}
}

#endif
