/**
 * \file unitigs_main.cpp: cover the haplotype k-mers of a GBZ with graph walks
 */

#include <getopt.h>
#include <omp.h>

#include <fstream>
#include <iostream>

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
         << "Scanning the FASTA output sees each haplotype k-mer exactly once." << endl
         << endl
         << "options:" << endl
         << "  -k, --kmer-length N   k-mer length, at least 2 (required)" << endl
         << "  -r, --gref            leave out k-mers that lie along a gref path (made by" << endl
         << "                        vg paths -u), so the output complements the gref" << endl
         << "  -g, --gaf FILE        write the graph walk of each tig to FILE as GAF" << endl
         << "  -t, --threads N       number of threads to use [all available]" << endl
         << "  -p, --progress        show progress and statistics" << endl
         << "  -h, --help            print this help message to stderr and exit" << endl;
}

int main_unitigs(int argc, char** argv) {
    Logger logger("vg unitigs");
    if (argc == 2) {
        help_unitigs(argv);
        return 1;
    }

    size_t k = 0;
    bool use_gref = false;
    string gaf_file;
    bool progress = false;

    int c;
    optind = 2; // force optind past command positional argument
    while (true) {
        static struct option long_options[] = {
            {"kmer-length", required_argument, 0, 'k'},
            {"gref", no_argument, 0, 'r'},
            {"gaf", required_argument, 0, 'g'},
            {"threads", required_argument, 0, 't'},
            {"progress", no_argument, 0, 'p'},
            {"help", no_argument, 0, 'h'},
            {0, 0, 0, 0}
        };
        int option_index = 0;
        c = getopt_long(argc, argv, "k:rg:t:ph?", long_options, &option_index);
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
            gaf_file = optarg;
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
            return 1;
        }
    }
    if (optind + 1 != argc) {
        logger.error() << "expected exactly one GBZ file" << endl;
    }
    if (k < 2) {
        logger.error() << "-k/--kmer-length must be at least 2" << endl;
    }
    string graph_file = require_exists(logger, argv[optind]);

    gbwtgraph::GBZ gbz;
    load_gbz(gbz, graph_file, progress);

    // Every path in the GBWT counts as a haplotype, gref paths included. That is harmless
    // while each gref path is still contained in the path it was cut from.
    auto unsupported = unsupported_gref_paths(gbz, GrefCover::gref_prefix);
    if (!unsupported.empty()) {
        size_t unsupported_length = 0;
        for (const auto& path : unsupported) {
            unsupported_length += path.second;
        }
        logger.warn() << unsupported.size() << " gref paths (" << unsupported_length << " bp, e.g. "
                      << unsupported.front().first << ") are not contained in any other path, so the"
                      << " output includes k-mers that no haplotype in " << graph_file << " has."
                      << " Was the GBZ subsampled after vg paths -u?" << endl;
    }

    unique_ptr<GrefIndex> gref;
    if (use_gref) {
        gref = make_unique<GrefIndex>(gbz, GrefCover::gref_prefix);
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
    }

    double start = gbwt::readTimer();
    UnitigStatistics statistics;
    vector<Tig> tigs = haplotype_unitigs(gbz.graph, k, gref.get(), &statistics);
    if (progress) {
        cerr << "[vg unitigs] " << statistics.unitigs << " unitigs, " << statistics.tigs << " tigs, "
             << statistics.bases << " bp, " << statistics.kmers << " " << k << "-mers, longest "
             << statistics.longest << " bp";
        if (gref) {
            cerr << ", " << statistics.unanchored << " unanchored";
        }
        cerr << " in " << (gbwt::readTimer() - start) << " s" << endl;
    }

    ofstream gaf;
    if (!gaf_file.empty()) {
        gaf.open(gaf_file);
        if (!gaf) {
            logger.error() << "could not open " << gaf_file << " for writing" << endl;
        }
    }
    for (size_t i = 0; i < tigs.size(); i++) {
        const Tig& tig = tigs[i];
        string name = "u" + to_string(i + 1);
        cout << '>' << name;
        if (!tig.gref.empty()) {
            cout << ' ' << tig.gref;
        }
        cout << '\n' << tig_sequence(gbz.graph, tig) << '\n';

        if (gaf.is_open()) {
            size_t path_length = 0;
            gaf << name << '\t' << tig.length << "\t0\t" << tig.length << "\t+\t";
            for (const handle_t& handle : tig.walk) {
                gaf << (gbz.graph.get_is_reverse(handle) ? '<' : '>') << gbz.graph.get_id(handle);
                path_length += gbz.graph.get_length(handle);
            }
            gaf << '\t' << path_length << '\t' << tig.offset << '\t' << (tig.offset + tig.length)
                << '\t' << tig.length << '\t' << tig.length << "\t255"
                << "\tsp:i:" << tig.support;
            if (!tig.gref.empty()) {
                gaf << "\tgr:Z:" << tig.gref;
            }
            gaf << '\n';
        }
    }
    cout.flush();

    return 0;
}

// Register subcommand
static Subcommand vg_unitigs("unitigs", "cover haplotype k-mers with walks for linear scanning", WIDGET, main_unitigs);
