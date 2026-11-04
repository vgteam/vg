#ifndef VG_SUBCOMMAND_GIRAFFE_OPTIONS_HPP_INCLUDED
#define VG_SUBCOMMAND_GIRAFFE_OPTIONS_HPP_INCLUDED

/**
 * \file giraffe_options.hpp
 * vg giraffe's mapping parameters and presets, for other subcommands that drive
 * a MinimizerMapper without copying them.
 */

#include <optional>
#include <string>
#include <utility>
#include <vector>

#include <gbwtgraph/minimizer.h>

namespace vg {

class MinimizerMapper;

namespace subcommand {

/// Names of vg giraffe's parameter presets ("default", "hifi", "r10", ...).
std::vector<std::string> giraffe_preset_names();

/// Throw std::invalid_argument if the preset or any override is unknown, so a
/// caller can reject bad parameters before loading indexes.
void check_giraffe_parameters(const std::string& preset,
                              const std::vector<std::pair<std::string, std::string>>& overrides);

/// True for presets that map long reads (hifi, r10, and blat, which builds on hifi).
bool giraffe_preset_maps_long_reads(const std::string& preset);

/// Whether to map recombination-aware, as vg giraffe decides it: rec_mode if
/// given, else on for long reads when the index stores paths. Throws
/// std::invalid_argument if the minimizer index cannot be used that way.
bool giraffe_rec_mode(const gbwtgraph::DefaultMinimizerIndex& minimizer_index,
                      bool map_long_reads, std::optional<bool> rec_mode);

/// Configure mapper as `vg giraffe --parameter-preset <preset>` would, then
/// apply overrides, given as (vg giraffe long option name, value) pairs.
/// Throws std::invalid_argument for an unknown preset or option.
void configure_giraffe_mapper(MinimizerMapper& mapper, const std::string& preset,
                              const std::vector<std::pair<std::string, std::string>>& overrides = {});

}
}

#endif
