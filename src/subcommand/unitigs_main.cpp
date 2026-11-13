/**
 * \file unitigs_main.cpp: cover the haplotype k-mers of a GBZ with graph walks
 */

#include <getopt.h>
#include <omp.h>

#include <algorithm>
#include <fstream>
#include <iostream>
#include <span>

#include "subcommand.hpp"
#include "../gbwtgraph_helper.hpp"
#include "../gref.hpp"
#include "../unitigs.hpp"
#include "../utility.hpp"

using namespace std;
using namespace vg;
using namespace vg::subcommand;

void help_unitigs(char** argv) {
    cerr << "usage: " << argv[0] << " unitigs [options] graph.gbz > tigs.fa" << endl
         << "Cover every k-mer spelled by the haplotypes in a GBZ with graph walks." << endl
         << "Scanning the FASTA output sees each canonical haplotype k-mer once." << endl
         << endl
         << "options:" << endl
         << "  -k, --kmer-length N   k-mer length, at least 2 (required)" << endl
         << "  -r, --gref            leave out k-mers that lie along a gref path (made by" << endl
         << "                        vg paths -u), so the output complements the gref" << endl
         << "  -g, --gaf FILE        write the graph walk of each tig to FILE as GAF, with" << endl
         << "                        the graph's name and this command in its header" << endl
         << "  -n, --name-prefix S   name the tigs S1, S2, ... [u]" << endl
         << "  -t, --threads N       number of threads to use [all available]" << endl
         << "  -p, --progress        show progress and statistics" << endl
         << "  -h, --help            print this help message to stderr and exit" << endl;
}

/// The options and arguments of vg unitigs.
struct UnitigsConfig {
    string graph_file;
    size_t k = 0;
    bool use_gref = false;
    string gaf_file;
    string name_prefix = "u";
    bool progress = false;

    /**
     * Read the command line.
     *
     * Print the help and exit on -h or an unknown option. Report an error on a
     * k below 2, and on a missing or extra argument.
     */
    UnitigsConfig(int argc, char** argv, const Logger& logger);
};

UnitigsConfig::UnitigsConfig(int argc, char** argv, const Logger& logger) {
    optind = 2; // force optind past command positional argument
    while (true) {
        static struct option long_options[] = {
            {"kmer-length", required_argument, nullptr, 'k'},
            {"gref", no_argument, nullptr, 'r'},
            {"gaf", required_argument, nullptr, 'g'},
            {"name-prefix", required_argument, nullptr, 'n'},
            {"threads", required_argument, nullptr, 't'},
            {"progress", no_argument, nullptr, 'p'},
            {"help", no_argument, nullptr, 'h'},
            {0, 0, 0, 0}
        };
        int option_index = 0;
        int c = getopt_long(argc, argv, "k:rg:n:t:ph?", long_options, &option_index);
        if (c == -1) {
            break;
        }
        switch (c) {
        case 'k':
            k = parse<size_t>(optarg);
            break;
        case 'r':
            use_gref = true;
            break;
        case 'g':
            gaf_file = ensure_writable(logger, optarg);
            break;
        case 'n':
            name_prefix = optarg;
            break;
        case 't':
            set_thread_count(logger, optarg);
            break;
        case 'p':
            progress = true;
            break;
        case 'h':
        case '?':
        default:
            help_unitigs(argv);
            std::exit(EXIT_FAILURE);
        }
    }
    if (optind + 1 != argc) {
        logger.error() << "expected exactly one GBZ file" << endl;
    }
    if (k < 2) {
        logger.error() << "-k/--kmer-length must be at least 2" << endl;
    }
    graph_file = require_exists(logger, argv[optind]);
}

/**
 * Warn if some gref paths are not contained in any other path.
 *
 * Every path in the GBWT counts as a haplotype, gref paths included. That is
 * harmless while each gref path is still contained in the path it was cut
 * from.
 */
static void warn_unsupported_gref(const gbwtgraph::GBZ& gbz, const string& graph_file, const Logger& logger) {
    if (auto unsupported = unsupported_gref_paths(gbz, GrefCover::gref_prefix); !unsupported.empty()) {
        size_t unsupported_length = 0;
        for (const auto& [name, length] : unsupported) {
            unsupported_length += length;
        }
        logger.warn() << unsupported.size() << " gref paths (" << unsupported_length << " bp, e.g. "
                      << unsupported.front().first << ") are not contained in any other path, so the"
                      << " output includes k-mers that no haplotype in " << graph_file << " has."
                      << " Was the GBZ subsampled after vg paths -u?" << endl;
    }
}

/// Index the gref paths of the GBZ, and report an error if it has none.
static unique_ptr<GrefIndex> load_gref(const gbwtgraph::GBZ& gbz, const string& graph_file, const Logger& logger,
                                       bool progress) {
    auto gref = make_unique<GrefIndex>(gbz, GrefCover::gref_prefix);
    if (gref->path_count() == 0) {
        logger.error() << "no gref paths (named " << GrefCover::gref_prefix << "*) in "
                       << graph_file << "; make them with vg paths -u" << endl;
    }
    if (progress) {
        cerr << "[vg unitigs] gref: " << gref->path_count() << " paths, "
             << gref->total_length() << " bp" << endl;
    }
    if (gref->repeated_visits() > 0) {
        logger.warn() << gref->repeated_visits() << " node visits on the gref were to nodes already on it;"
                      << " k-mers along them are reported as tigs" << endl;
    }
    return gref;
}

/**
 * Write the header of the GAF of the tigs.
 *
 * The graph's name lets vg lift check that the tigs belong to its graph. The
 * command says which cover the tig names refer to.
 */
static void write_gaf_header(ostream& gaf, const gbwtgraph::GBZ& gbz, span<char* const> command) {
    for (const string& line : gbz.graph_name().gaf_header_lines()) {
        gaf << line << '\n';
    }
    gaf << "@CO\t";
    for (size_t i = 0; i < command.size(); i++) {
        string argument = command[i];
        replace(argument.begin(), argument.end(), '\t', ' ');
        gaf << (i > 0 ? " " : "") << argument;
    }
    gaf << '\n';
}

int main_unitigs(int argc, char** argv) {
    Logger logger("vg unitigs");
    if (argc == 2) {
        help_unitigs(argv);
        return 1;
    }
    UnitigsConfig config(argc, argv, logger);

    gbwtgraph::GBZ gbz;
    load_gbz(gbz, config.graph_file, config.progress);
    warn_unsupported_gref(gbz, config.graph_file, logger);
    unique_ptr<GrefIndex> gref;
    if (config.use_gref) {
        gref = load_gref(gbz, config.graph_file, logger, config.progress);
    }

    double start = gbwt::readTimer();
    UnitigStatistics statistics;
    vector<Tig> tigs = haplotype_unitigs(gbz.graph, config.k, gref.get(), &statistics);
    if (config.progress) {
        cerr << "[vg unitigs] " << statistics.unitigs << " unitigs, " << statistics.tigs << " tigs, "
             << statistics.bases << " bp, " << statistics.kmers << " " << config.k << "-mers, longest "
             << statistics.longest << " bp";
        if (gref) {
            cerr << ", " << statistics.unanchored << " unanchored";
        }
        cerr << " in " << (gbwt::readTimer() - start) << " s" << endl;
    }

    ofstream gaf;
    if (!config.gaf_file.empty()) {
        gaf.open(config.gaf_file);
        if (!gaf) {
            logger.error() << "could not open " << config.gaf_file << " for writing" << endl;
        }
        write_gaf_header(gaf, gbz, span<char* const>(argv, static_cast<size_t>(argc)));
    }
    for (size_t i = 0; i < tigs.size(); i++) {
        const Tig& tig = tigs[i];
        string name = config.name_prefix + to_string(i + 1);
        const string& placement = tig.placement.text;
        cout << '>' << name;
        if (!placement.empty()) {
            cout << ' ' << placement;
        }
        cout << '\n' << tig.interval.sequence(gbz.graph) << '\n';

        if (gaf.is_open()) {
            write_gaf_columns(gaf, gbz.graph, name, tig.interval);
            gaf << "\tsp:i:" << tig.support;
            if (!placement.empty()) {
                gaf << "\tgr:Z:" << placement;
            }
            gaf << '\n';
        }
    }
    cout.flush();

    return 0;
}

// Register subcommand
static Subcommand vg_unitigs("unitigs", "cover haplotype k-mers with walks for linear scanning", WIDGET, main_unitigs);
