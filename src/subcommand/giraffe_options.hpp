#ifndef VG_SUBCOMMAND_GIRAFFE_OPTIONS_HPP_INCLUDED
#define VG_SUBCOMMAND_GIRAFFE_OPTIONS_HPP_INCLUDED

/**
 * \file giraffe_options.hpp
 * vg giraffe's mapping parameters and presets, for other subcommands that run
 * a MinimizerMapper.
 */

#include <optional>
#include <string>
#include <utility>
#include <vector>
#include <memory>

#include <gbwtgraph/minimizer.h>
#include "options.hpp"

namespace vg {

class MinimizerMapper;

namespace subcommand {

/// Get the names of vg giraffe's parameter presets.
std::vector<std::string> giraffe_preset_names();


/// mapping parameters
struct GiraffeParameters {
    // name of the giraffe parameter preset
    std::string preset;

    /// options set on top of the presets. Each one of the options is a giraffe long option name without the dashes
    /// example for one vector element: {"max-multimaps", "20"}. A flag has no value.
    /// Also the options replace preset values if duplicated.
    std::vector<std::pair<std::string, std::string>> options;

    /// remcombination aware mode on or off.
    /// It is separate because vg giraffe does not keep --rec-mode in the option parser.
    std::optional<bool> rec_mode;
};


/// make vg giraffe's option parser and fill it. This helps when we cannot use the argv for all the parameters
/// like what vg giraffe does right now.
/// Throws std::invalid_argument for unknown preset or option or bad value
std::unique_ptr<GroupedOptionGroup> make_giraffe_parser(const GiraffeParameters& params);

/// Return true if the preset maps long reads.
bool giraffe_preset_maps_long_reads(const std::string& preset);

/**
 * Decide whether to map recombination-aware.
 *
 * Use rec_mode if it is set. Otherwise map recombination-aware for long reads
 * when the minimizer payload stores paths.
 *
 * @throws std::invalid_argument if recombination-aware mapping is asked for
 * with short reads, or the minimizer payload does not suit the chosen mode.
 */
bool giraffe_rec_mode(const gbwtgraph::DefaultMinimizerIndex& minimizer_index,
                      bool map_long_reads, std::optional<bool> rec_mode);

/// Write the parser's values into the mapper, check that the range options
/// are sensible, and set the alignment scores and recombination-aware mode.
/// use_rec_mode is not a parser option. It comes from giraffe_rec_mode(),
/// which needs the loaded minimizer index.
void set_up_giraffe_mapper(GroupedOptionGroup& parser, bool use_rec_mode, MinimizerMapper& mapper);

}
}

#endif
