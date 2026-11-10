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

#include <gbwtgraph/minimizer.h>
#include "options.hpp"

namespace vg {

class MinimizerMapper;

namespace subcommand {

/// Write the parser's values into the mapper, check that the range options
/// are sensible, and set the alignment scores and recombination-aware mode.
/// use_rec_mode is not a parser option. It comes from giraffe_rec_mode(),
/// which needs the loaded minimizer index.
void set_up_giraffe_mapper(GroupedOptionGroup& parser, bool use_rec_mode, MinimizerMapper& mapper);

/// Get the names of vg giraffe's parameter presets.
std::vector<std::string> giraffe_preset_names();

/**
 * Check a preset and its overrides, so bad parameters fail before any index
 * is loaded.
 *
 * @throws std::invalid_argument for an unknown preset or option, or a value
 * that does not parse.
 */
void check_giraffe_parameters(const std::string& preset,
                              const std::vector<std::pair<std::string, std::string>>& overrides);

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

/**
 * Set the mapper's parameters and alignment scores as vg giraffe does for a
 * preset and overrides.
 *
 * Each override is a vg giraffe long option name, without dashes, and its
 * value. An override replaces the preset's value.
 *
 * @throws std::invalid_argument for an unknown preset or option, or a value
 * that does not parse.
 */
void configure_giraffe_mapper(MinimizerMapper& mapper, const std::string& preset,
                              const std::vector<std::pair<std::string, std::string>>& overrides = {});

}
}

#endif
