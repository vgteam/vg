/** \file combine_main.cpp
 *
 * Defines the "vg combine" subcommand
 */


#include <unistd.h>
#include <getopt.h>

#include <iostream>
#include <string>
#include <vector>

#include <vg/io/vpkg.hpp>

#include "subcommand.hpp"

#include "../combine.hpp"
#include "../handle.hpp"
#include "../io/save_handle_graph.hpp"

using namespace std;
using namespace vg;
using namespace vg::subcommand;

void help_combine(char** argv) {
    cerr << "usage: " << argv[0] << " combine [options] <graph1.vg> [graph2.vg ...] >merged.vg" << endl
         << "Combines one or more graphs into a single file." << endl
         << "Input graphs must not share any path name." << endl
         << endl
         << "Options:" << endl
         << "  -s, --seam POLICY     how to resolve overlap between adjacent chunks:" << endl
         << "                          renumber  assume no overlap [default]" << endl
         << "                          shared    overlap by shared nodes" << endl
         << "                          trim      overlap in REFERENCE coordinates" << endl
         << "  -u, --fuse            requires --seam trim. Fuse the boundary nodes into" << endl
         << "                        one node instead of adding an edge." << endl
         << "  -P, --phase-block-is-offset" << endl
         << "                        treat a path's phase block as a start coordinate, so" << endl
         << "                        fragments differing only by phase block are merged." << endl
         << "  -h, --help            print this help message to stderr and exit" << endl;
}

/// Map a --seam argument to its policy, or exit naming the legal values.
static SeamPolicy parse_seam_policy(const string& value, const Logger& logger) {
    if (value == "renumber") {
        return SeamPolicy::RENUMBER;
    }
    if (value == "shared") {
        return SeamPolicy::SHARED_IDS;
    }
    if (value == "trim") {
        return SeamPolicy::COORD_TRIM;
    }
    logger.error() << "unrecognized --seam policy \"" << value
                   << "\"; expected one of: renumber, shared, trim." << endl;
    return SeamPolicy::RENUMBER; // not reached; logger.error() exits
}

int main_combine(int argc, char** argv) {
    Logger logger("vg combine");

    if (argc == 2) {
        help_combine(argv);
        return 1;
    }

    SeamPolicy seam = SeamPolicy::RENUMBER;
    bool fuse = false;
    bool phase_block_is_offset = false;

    int c;
    optind = 2; // force optind past command positional argument

    // Parse command line options.
    while (true) {
        static struct option long_options[] =
        {
            {"help", no_argument, 0, 'h'},
            {"seam", required_argument, 0, 's'},
            {"fuse", no_argument, 0, 'u'},
            {"phase-block-is-offset", no_argument, 0, 'P'},
            {0, 0, 0, 0}
        };

        int option_index = 0;
        c = getopt_long (argc, argv, "h?s:uP",
                         long_options, &option_index);

        // Detect the end of the options.
        if (c == -1)
            break;

        switch (c)
        {
        case 's':
            seam = parse_seam_policy(optarg, logger);
            break;
        case 'u':
            fuse = true;
            break;
        case 'P':
            phase_block_is_offset = true;
            break;
        case 'h':
        case '?':
            help_combine(argv);
            exit(1);
            break;

        default:
            abort ();
        }
    }

    // --fuse requires --seam trim.
    if (fuse && seam != SeamPolicy::COORD_TRIM) {
        logger.error() << "--fuse/-u requires --seam trim." << endl;
    }

    vector<GraphCombiner::Input> inputs;
    while (optind < argc) {
        GraphCombiner::Input in;
        in.name = get_input_file_name(optind, argc, argv);
        in.graph = vg::io::VPKG::load_one<MutablePathDeletableHandleGraph>(in.name);
        inputs.push_back(std::move(in));
    }

    GraphCombiner combiner(seam, fuse, phase_block_is_offset);
    unique_ptr<MutablePathDeletableHandleGraph> combined =
        combiner.combine(std::move(inputs));
    vg::io::save_handle_graph(combined.get(), cout);
    return 0;
}

// Register subcommand
static Subcommand vg_combine("combine", "merge multiple graph files together", main_combine);
