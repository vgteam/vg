/** \file safari_main.cpp
 *
 * Defines the "vg safari" subcommand, which maps ancient DNA reads with Giraffe plus SAFARI's
 * damage-aware seeding. It takes every Giraffe option, and these of its own:
 *
 *   -q, --rymer-name FILE           RYmer index from vg rymer (required)
 *       --deam-5p FILE              5' deamination profile (required)
 *       --deam-3p FILE              3' deamination profile (required)
 *   -j, --posterior-threshold X     accept a graph k-mer above this posterior [0.5]
 *   -V, --spurious-prior X          prior probability of a spurious RYmer match [0.5]
 *
 * It takes these out of the command line, loads what they name, and hands the rest to Giraffe,
 * which calls back to attach the seed rescue to the mapper. See artifacts/rymer_seed_rescue.hpp.
 */

#include "subcommand.hpp"
#include "giraffe_main.hpp"

#include <fstream>
#include <iostream>
#include <memory>
#include <optional>
#include <stdexcept>
#include <string>
#include <vector>

#include <vg/io/vpkg.hpp>
#include <gbwtgraph/minimizer.h>

#include "../log.hpp"
#include "../minimizer_mapper.hpp"
#include "../utility.hpp"
#include "../../artifacts/deamination_profile.hpp"
#include "../../artifacts/rymer_index.hpp"
#include "../../artifacts/rymer_seed_rescue.hpp"

using namespace vg;
using namespace vg::subcommand;

namespace {

/// SAFARI keeps minimizers with many more hits than Giraffe does, since RYmer seeds come from the
/// graph k-mers and do not share the read's errors.
const std::string SAFARI_DEFAULT_HARD_HIT_CAP = "50000";

constexpr double SAFARI_DEFAULT_POSTERIOR_THRESHOLD = 0.5;
constexpr double SAFARI_DEFAULT_SPURIOUS_PRIOR = 0.5;

void print_safari_help(std::ostream& out) {
    out << "usage:" << std::endl
        << "  vg safari -Z graph.gbz -d graph.dist -m graph.min -q graph.ry \\" << std::endl
        << "      --deam-5p sample_5p.prof --deam-3p sample_3p.prof -f reads.fq [options] > output.gam" << std::endl
        << std::endl
        << "Maps ancient DNA reads with Giraffe, adding seeds that deamination damage would hide:" << std::endl
        << "a read k-mer is matched to graph k-mers with the same purine/pyrimidine (RYmer)" << std::endl
        << "encoding, and a match is kept when a damage model makes it likely enough." << std::endl
        << "Only single-end reads get these extra seeds; merge read pairs first." << std::endl
        << std::endl
        << "SAFARI options:" << std::endl
        << "  -q, --rymer-name FILE         use this RYmer index, from vg rymer (required)" << std::endl
        << "      --deam-5p FILE            5' deamination profile, a .prof table (required)" << std::endl
        << "      --deam-3p FILE            3' deamination profile, a .prof table (required)" << std::endl
        << "  -j, --posterior-threshold X   keep an RYmer match whose posterior probability of" << std::endl
        << "                                being damage is above X [" << SAFARI_DEFAULT_POSTERIOR_THRESHOLD << "]" << std::endl
        << "  -V, --spurious-prior X        prior probability that an RYmer match is spurious [" << SAFARI_DEFAULT_SPURIOUS_PRIOR << "]" << std::endl
        << std::endl
        << "vg safari sets Giraffe's --hard-hit-cap to " << SAFARI_DEFAULT_HARD_HIT_CAP
        << " unless it, or a parameter preset, is given." << std::endl
        << std::endl;
    print_safari_citation(out);
    out << std::endl
        << "All vg giraffe options below are also accepted." << std::endl
        << std::endl;
}

/// The SAFARI settings found on the command line.
struct SafariOptions {
    std::string rymer_name;
    std::string five_prime_profile;
    std::string three_prime_profile;
    double posterior_threshold = SAFARI_DEFAULT_POSTERIOR_THRESHOLD;
    double spurious_prior = SAFARI_DEFAULT_SPURIOUS_PRIOR;
    bool wants_help = false;
    bool sets_hard_hit_cap = false;
    bool sets_preset = false;
    bool shows_progress = false;
};

/// If `args[i]` is option `short_name` (none if '\0') or one of `long_names`, given as "-x VALUE",
/// "-xVALUE", "--name VALUE" or "--name=VALUE", returns its value and advances `i` past it.
/// Otherwise returns nothing and leaves `i` alone. Exits with an error if the value is missing.
std::optional<std::string> take_value(const std::vector<std::string>& args, size_t& i, char short_name,
                                      const std::vector<std::string>& long_names, const Logger& logger) {
    const std::string& arg = args[i];
    for (const std::string& name : long_names) {
        std::string flag = "--" + name;
        if (arg == flag) {
            if (i + 1 >= args.size()) {
                logger.error() << "option " << flag << " needs a value" << std::endl;
            }
            i += 2;
            return args[i - 1];
        }
        if (arg.compare(0, flag.size() + 1, flag + "=") == 0) {
            i += 1;
            return arg.substr(flag.size() + 1);
        }
    }
    if (short_name != '\0' && arg.size() >= 2 && arg[0] == '-' && arg[1] == short_name) {
        if (arg.size() > 2) {
            i += 1;
            return arg.substr(2);
        }
        if (i + 1 >= args.size()) {
            logger.error() << "option -" << short_name << " needs a value" << std::endl;
        }
        i += 2;
        return args[i - 1];
    }
    return std::nullopt;
}

bool is_option(const std::string& arg, char short_name, const std::string& long_name) {
    return (arg.size() >= 2 && arg[0] == '-' && arg[1] == short_name)
        || arg == "--" + long_name
        || arg.compare(0, long_name.size() + 3, "--" + long_name + "=") == 0;
}

double parse_probability(const std::string& value, const std::string& option, const Logger& logger) {
    // parse() exits with an error if the value is not a number.
    double result = parse<double>(value);
    if (result < 0.0 || result > 1.0) {
        logger.error() << "option " << option << " must be from 0 to 1" << std::endl;
    }
    return result;
}

/// Takes the SAFARI options out of `args`, which holds the command line after "vg safari", and
/// returns them. Leaves everything else, in order, for Giraffe.
SafariOptions take_safari_options(std::vector<std::string>& args, const Logger& logger) {
    SafariOptions options;
    std::vector<std::string> rest;
    size_t i = 0;
    while (i < args.size()) {
        const std::string& arg = args[i];
        if (arg == "--") {
            // Everything after this is a positional argument.
            rest.insert(rest.end(), args.begin() + i, args.end());
            break;
        }
        std::optional<std::string> value;
        if ((value = take_value(args, i, 'q', {"rymer-name"}, logger))) {
            options.rymer_name = *value;
        } else if ((value = take_value(args, i, '\0', {"deam-5p"}, logger))) {
            options.five_prime_profile = *value;
        } else if ((value = take_value(args, i, '\0', {"deam-3p"}, logger))) {
            options.three_prime_profile = *value;
        } else if ((value = take_value(args, i, 'j', {"posterior-threshold", "posterior-odds-threshold"}, logger))) {
            options.posterior_threshold = parse_probability(*value, "--posterior-threshold", logger);
        } else if ((value = take_value(args, i, 'V', {"spurious-prior", "spurious-alignment-prior"}, logger))) {
            options.spurious_prior = parse_probability(*value, "--spurious-prior", logger);
        } else {
            options.wants_help |= (arg == "-h" || arg == "-?" || arg == "--help");
            options.sets_hard_hit_cap |= is_option(arg, 'C', "hard-hit-cap");
            options.sets_preset |= is_option(arg, 'b', "parameter-preset");
            options.shows_progress |= (arg == "-p" || arg == "--progress");
            rest.push_back(arg);
            i++;
        }
    }
    args = std::move(rest);
    return options;
}

void require_readable(const std::string& filename, const std::string& what, const Logger& logger) {
    if (!std::ifstream(filename).is_open()) {
        logger.error() << "cannot open " << what << " " << filename << std::endl;
    }
}

}

int main_safari(int argc, char** argv) {
    Logger logger("vg safari");

    GiraffeVariant variant;
    variant.name = "safari";
    variant.print_help = print_safari_help;

    // With no arguments, Giraffe prints the short help, with ours first.
    if (argc <= 2) {
        return run_giraffe(argc, argv, variant);
    }

    std::vector<std::string> args(argv + 2, argv + argc);
    SafariOptions options = take_safari_options(args, logger);

    // Loaded once here and shared by every mapping thread. They stay alive until Giraffe, and
    // the mapper it makes, are done.
    std::unique_ptr<gbwtgraph::DefaultMinimizerIndex> rymer_index;
    std::unique_ptr<DeaminationProfile> profile;
    std::unique_ptr<RymerSeedRescue> rescue;

    if (!options.wants_help) {
        if (options.rymer_name.empty()) {
            logger.error() << "an RYmer index (-q) is required; build one with vg rymer" << std::endl;
        }
        if (options.five_prime_profile.empty() || options.three_prime_profile.empty()) {
            logger.error() << "deamination profiles for both ends (--deam-5p and --deam-3p) are required" << std::endl;
        }
        require_readable(options.rymer_name, "RYmer index", logger);
        require_readable(options.five_prime_profile, "5' deamination profile", logger);
        require_readable(options.three_prime_profile, "3' deamination profile", logger);

        if (options.shows_progress) {
            logger.info() << "Loading RYmer index from " << options.rymer_name << std::endl;
        }
        rymer_index = vg::io::VPKG::load_one<gbwtgraph::DefaultMinimizerIndex>(options.rymer_name);
        if (!is_rymer_index(*rymer_index)) {
            logger.error() << options.rymer_name << " is not an RYmer index; build one with vg rymer" << std::endl;
        }
        if (rymer_index->uses_syncmers()) {
            logger.error() << options.rymer_name << " uses syncmers, which SAFARI does not support" << std::endl;
        }

        try {
            profile = std::make_unique<DeaminationProfile>(options.five_prime_profile, options.three_prime_profile);
        } catch (const std::runtime_error& e) {
            logger.error() << e.what() << std::endl;
        }

        DamageModelParameters model;
        model.spurious_prior = options.spurious_prior;
        rescue = std::make_unique<RymerSeedRescue>(*rymer_index, *profile, options.posterior_threshold, model);

        if (!options.sets_hard_hit_cap && !options.sets_preset) {
            args.insert(args.begin(), {"--hard-hit-cap", SAFARI_DEFAULT_HARD_HIT_CAP});
        }
    }

    variant.configure_mapper = [&rescue, &logger](MinimizerMapper& mapper, bool paired) {
        if (!rescue) {
            return;
        }
        if (paired) {
            logger.warn() << "RYmer seeds are only added to single-end reads; "
                          << "these paired reads are mapped as by vg giraffe" << std::endl;
        }
        const RymerSeedRescue* seed_rescue = rescue.get();
        mapper.stand_in_sequences = [seed_rescue](const Alignment& aln) {
            return seed_rescue->stand_in_sequences(aln.sequence());
        };
    };

    // Rebuild the command line for Giraffe. getopt may reorder the pointers, so keep the
    // strings they point into alive and separate.
    std::vector<std::string> giraffe_args;
    giraffe_args.emplace_back(argv[0]);
    giraffe_args.emplace_back(argv[1]);
    giraffe_args.insert(giraffe_args.end(), args.begin(), args.end());
    std::vector<char*> giraffe_argv;
    for (std::string& arg : giraffe_args) {
        giraffe_argv.push_back(arg.data());
    }
    giraffe_argv.push_back(nullptr);

    return run_giraffe(static_cast<int>(giraffe_argv.size() - 1), giraffe_argv.data(), variant);
}

// Register subcommand
static Subcommand vg_safari("safari", "map ancient DNA reads with damage-aware seeding (SAFARI)",
                            PIPELINE, 7, main_safari);
