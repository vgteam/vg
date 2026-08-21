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
 * and print all the subcommands that exist (via a help subcommand), and we
 * can do "vg help subcommand" by calling a saved help function.
 *
 * We have a subcommand importance/category system, so we can tell people about
 * just the main pipeline and keep the subcommands they don't want out of their
 * brains and off their screen.
 *
 * Subcommands get passed all of argv, so they have to skip past their names
 * when parsing arguments.
 *
 * To make a subcommand, do something like this in a *_main.cpp file in this
 * "subcommand" directory:
 * 
 *     #include "subcommand.hpp"
 *     using namespace vg::subcommand;
 * 
 *     void help_frobnicate(char** argv) {
 *         cerr << "usage: " << argv[0] << " frobnicate" << endl
 *              << "Foo the bar" << endl
 *              << "  -b, --bar          the bar to foo" << endl
 *              << "  -h, --help         print this help message to stderr and exit" << endl
 *              << endl;
 *     }
 * 
 *     int main_frobnicate(int argc, char** argv) {
 *         return 0;
 *     }
 *
 *     static Subcommand vg_frobnicate("frobnicate", "frobnicate nodes and edges",
 *                                     help_frobnicate, main_frobnicate);
 * 
 * All src/subcommand/{subcommand}_main.cpp files must pass the checks
 * (formatting etc.) in scripts/lint.py as part of an automated test.
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

/**
 * The sub-lists on the manpage to organize subcommands
 */
enum ManpageSection {
    // Graph construction and indexing
    SET_UP_GRAPH,
    // Read mapping
    MAP_READS,
    // Downstream analyses
    DOWNSTREAM,
    // Working with read alignments
    MANIPULATE_ALN,
    // Graph and read statistics
    GET_STATS,
    // Manipulate a graph
    MANIPULATE_GRAPH,
    // Conversion between formats
    CONVERT_FORMAT,
    // Subgraph extraction
    EXTRACT_GRAPH,
    // Extremely specific analyses
    RARE_NEEDS,
    // Developer tools
    DEV_TOOLS
};

/// Represents a list item in the manpage
struct manpage_item {
    /// Sub-list to put this item in
    ManpageSection section;
    /// Description of subcommand usage
    std::string blurb;
    /// A relevant wiki page (optional)
    std::string wiki_link;
};

const static std::vector<std::pair<ManpageSection, std::string>> MANPAGE_CATEGORY_HEADERS{
    {SET_UP_GRAPH, "Graph construction and indexing"},
    {MAP_READS, "Read mapping"},
    {DOWNSTREAM, "Downstream analyses"},
    {MANIPULATE_ALN, "Work with read alignments"},
    {GET_STATS, "Graph and read statistics"},
    {MANIPULATE_GRAPH, "Manipulate a graph"},
    {CONVERT_FORMAT, "Convert between formats"},
    {EXTRACT_GRAPH, "Subgraph extraction"},
    {RARE_NEEDS, "Extremely specific analyses"},
    {DEV_TOOLS, "Developer tools"}
};

const static std::map<std::string, std::string> REMOVED_CMD_MESSAGES{
    {"explode", std::string("Please use \"vg chunk -C source.vg -b part_dir/component\" "
                            "for the same* functionality as \"vg explode source.vg part_dir\"\n"
                            "* (unlike explode, the output directory must already exist when running chunk)")
    },
    {"msga", std::string("vg msga was an early prototype for constructing genome graphs "
                         "from multiple sequence alignments, but VG team members have developed "
                         "improved graph construction algorithms in Cactus and PGGB, "
                         "and several other tools have been developed by other groups.")}
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
     * the given category, with the given priority (lower is better), 
     * with the given manpage entries,
     * with the given default help function (used for vg help)
     * which calls the given main function when invoked.
     */
    Subcommand(std::string name, std::string description,
        CommandCategory category, int priority,
        std::vector<manpage_item> manpage_entries,
        std::function<void(char**)> help_function,
        std::function<int(int, char**)> main_function);
    
    /**
     * Make and register a subcommand with the given name and description, in
     * the given category, with the given manpage entries,
     * with the given help function, with worst priority,
     * which calls the given main function when invoked.
     */
    Subcommand(std::string name, std::string description,
        CommandCategory category,
        std::vector<manpage_item> manpage_entries,
        std::function<void(char**)> help_function,
        std::function<int(int, char**)> main_function);

    /**
     * Make and register a subcommand with the given name and description,
     * with the given help function, in the DEPRECATED category, with worst priority,
     * which calls the given main function when invoked.
     */
    Subcommand(std::string name, std::string description,
        std::function<void(char**)> help_function,
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
     * Get the manpage elements for this subcommand
     */
    const std::vector<manpage_item>& get_manpage_entries() const;
    
    /**
     * Get the priority level of a subcommand (lower is more important).
     */
    const int& get_priority() const;

    /**
     * Run the subcommand's help function (print to stderr)
     */
    const void run_help(char** argv) const;
    
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
     * constructed on first use. Note that at shutdown some of the pointers in
     * the registry may be to already-destructed static objects.
     */
    static std::map<std::string, Subcommand*>& get_registry();
    
    // These hold the actual fields defining the subcommand
    std::string name;
    std::string description;
    CommandCategory category;
    int priority;
    std::function<void(char**)> help_function;
    std::function<int(int, char**)> main_function;
    /// Things to put in the manpage (vg help --man)
    /// Stored as (section, blurb, wiki link)
    /// Not shown for DEPRECATED subcommands
    std::vector<manpage_item> manpage_entries;
};

}
}

#endif
