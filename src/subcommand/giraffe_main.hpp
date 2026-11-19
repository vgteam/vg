#ifndef VG_SUBCOMMAND_GIRAFFE_MAIN_HPP_INCLUDED
#define VG_SUBCOMMAND_GIRAFFE_MAIN_HPP_INCLUDED

/** \file
 * Lets another subcommand run the Giraffe mapper with its own additions, as `vg safari` does.
 */

#include <functional>
#include <ostream>
#include <string>

namespace vg {

class MinimizerMapper;

namespace subcommand {

/// Additions that a subcommand makes to Giraffe. The defaults give plain `vg giraffe`.
struct GiraffeVariant {
    /// The subcommand name shown in usage lines.
    std::string name = "giraffe";
    /// If set, prints the subcommand's own help before Giraffe's.
    std::function<void(std::ostream&)> print_help;
    /// If set, changes the mapper once Giraffe has made it and before any read is mapped. It is
    /// told whether the reads are paired.
    std::function<void(MinimizerMapper&, bool paired)> configure_mapper;
};

/// Runs Giraffe on the command line, as `vg giraffe` does, with the variant's additions.
int run_giraffe(int argc, char** argv, const GiraffeVariant& variant);

}
}

#endif
