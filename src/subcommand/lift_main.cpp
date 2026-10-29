/**
 * \file lift_main.cpp: lift features found on a cover to the graph, the gref and the assemblies
 */

#include <getopt.h>
#include <omp.h>

#include <fstream>
#include <iostream>
#include <unordered_set>

#include "subcommand.hpp"
#include "../cover_lift.hpp"
#include "../gbwtgraph_helper.hpp"
#include "../gref.hpp"
#include "../utility.hpp"

using namespace std;
using namespace vg;
using namespace vg::subcommand;

void help_lift(char** argv) {
    cerr << "usage: " << argv[0] << " lift [options] graph.gbz features.bed > features.gaf" << endl
         << "Lift features found by scanning a cover (the gref paths of the GBZ and the" << endl
         << "tigs of vg unitigs) to the graph and the gref, and to each path that visits" << endl
         << "them. A feature found on several cover sequences is reported once." << endl
         << endl
         << "Each feature gets an ID (f1, f2, ...) and a GAF record with that query name," << endl
         << "its graph walk, and the tags sp:i: (visits by paths other than gref paths)," << endl
         << "gr:Z: (gref placement, see vg unitigs) and an:Z: (BED name, if any)." << endl
         << endl
         << "options:" << endl
         << "  -u, --tigs FILE         GAF of the tigs, from vg unitigs -g" << endl
         << "  -a, --assemblies FILE   write each place where a path visits a feature to FILE" << endl
         << "                          as BED on the path's full sequence (for haplotypes," << endl
         << "                          the assembly contig): the feature's name, score and" << endl
         << "                          strand (., 0 and . if it has none), any further BED" << endl
         << "                          columns, and then its ID" << endl
         << "  -Q, --paths-by PREFIX   only write places on paths whose name starts with" << endl
         << "                          PREFIX (may repeat)" << endl
         << "  -t, --threads N         number of threads to use [all available]" << endl
         << "  -p, --progress          show progress" << endl
         << "  -h, --help              print this help message to stderr and exit" << endl;
}

int main_lift(int argc, char** argv) {
    Logger logger("vg lift");
    if (argc == 2) {
        help_lift(argv);
        return 1;
    }

    string tigs_file;
    string assemblies_file;
    vector<string> path_prefixes;
    bool progress = false;

    int c;
    optind = 2; // force optind past command positional argument
    while (true) {
        static struct option long_options[] = {
            {"tigs", required_argument, 0, 'u'},
            {"assemblies", required_argument, 0, 'a'},
            {"paths-by", required_argument, 0, 'Q'},
            {"threads", required_argument, 0, 't'},
            {"progress", no_argument, 0, 'p'},
            {"help", no_argument, 0, 'h'},
            {0, 0, 0, 0}
        };
        int option_index = 0;
        c = getopt_long(argc, argv, "u:a:Q:t:ph?", long_options, &option_index);
        if (c == -1) {
            break;
        }
        switch (c) {
        case 'u':
            tigs_file = require_exists(logger, optarg);
            break;
        case 'a':
            assemblies_file = optarg;
            break;
        case 'Q':
            path_prefixes.push_back(optarg);
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
            help_lift(argv);
            return 1;
        }
    }
    if (optind + 2 != argc) {
        logger.error() << "expected a GBZ file and a BED file" << endl;
    }
    string graph_file = require_exists(logger, argv[optind]);
    string features_file = argv[optind + 1];

    gbwtgraph::GBZ gbz;
    load_gbz(gbz, graph_file, progress);
    const gbwtgraph::GBWTGraph& graph = gbz.graph;

    GrefIndex gref(gbz, GrefCover::gref_prefix);
    if (progress) {
        cerr << "[vg lift] gref: " << gref.path_count() << " paths, " << gref.total_length() << " bp" << endl;
    }
    Cover cover(graph, gref.path_count() > 0 ? &gref : nullptr);

    vector<CoverFeature> cover_features;
    get_input_file(features_file, [&](istream& in) {
        try {
            cover_features = read_bed_features(in);
        } catch (const runtime_error& e) {
            logger.error() << e.what() << endl;
        }
    });
    if (progress) {
        cerr << "[vg lift] " << cover_features.size() << " features" << endl;
    }

    // Load only the tigs that features are on.
    unordered_set<string> wanted;
    for (const CoverFeature& feature : cover_features) {
        if (!cover.is_gref_path(feature.sequence)) {
            wanted.insert(feature.sequence);
        }
    }
    if (!wanted.empty() && !tigs_file.empty()) {
        get_input_file(tigs_file, [&](istream& in) {
            try {
                cover.add_tigs(in, [&](const string& name) { return wanted.count(name) > 0; });
            } catch (const runtime_error& e) {
                logger.error() << e.what() << endl;
            }
        });
    }
    for (const CoverFeature& feature : cover_features) {
        optional<size_t> length = cover.sequence_length(feature.sequence);
        if (!length) {
            logger.error() << "feature on " << feature.sequence << ", which is neither a gref path in "
                           << graph_file << " nor a tig" << (tigs_file.empty() ? " (use -u)" : " in " + tigs_file) << endl;
        }
        if (feature.end > *length) {
            logger.error() << "feature " << feature.sequence << ":" << feature.start << "-" << feature.end
                           << " ends past the end of the sequence (" << *length << " bp)" << endl;
        }
    }

    vector<GraphFeature> features(cover_features.size());
    #pragma omp parallel for schedule(dynamic, 1024)
    for (size_t i = 0; i < cover_features.size(); i++) {
        features[i] = cover.lift(cover_features[i]);
    }
    cover_features = vector<CoverFeature>();
    sort_and_merge_features(features);
    if (progress) {
        cerr << "[vg lift] " << features.size() << " distinct features in the graph" << endl;
    }

    auto feature_id = [](size_t i) {
        return "f" + to_string(i + 1);
    };
    for (size_t i = 0; i < features.size(); i++) {
        const GraphFeature& feature = features[i];
        write_gaf_columns(cout, graph, feature_id(i), feature.interval);
        cout << "\tsp:i:" << feature.support;
        if (!feature.placement.text.empty()) {
            cout << "\tgr:Z:" << feature.placement.text;
        }
        if (!feature.fields.empty()) {
            cout << "\tan:Z:" << feature.fields.front();
        }
        cout << '\n';
    }
    cout.flush();

    if (!assemblies_file.empty()) {
        ofstream out(assemblies_file);
        if (!out) {
            logger.error() << "could not open " << assemblies_file << " for writing" << endl;
        }
        auto include = [&](const string& name) {
            if (name.compare(0, GrefCover::gref_prefix.size(), GrefCover::gref_prefix) == 0) {
                return false;
            }
            if (path_prefixes.empty()) {
                return true;
            }
            for (const string& prefix : path_prefixes) {
                if (name.compare(0, prefix.size(), prefix) == 0) {
                    return true;
                }
            }
            return false;
        };
        size_t written = 0;
        for_each_path_occurrences(gbz, features, include, [&](const string& name, const vector<PathOccurrence>& occurrences) {
            for (const PathOccurrence& occurrence : occurrences) {
                const GraphFeature& feature = features[occurrence.feature];
                const vector<string>& fields = feature.fields;
                out << name << '\t' << occurrence.start << '\t' << (occurrence.start + feature.interval.length)
                    << '\t' << (fields.size() > 0 ? fields[0] : ".")
                    << '\t' << (fields.size() > 1 ? fields[1] : "0") << '\t';
                if (feature.is_stranded()) {
                    out << (occurrence.is_reverse ? '-' : '+');
                } else {
                    out << '.';
                }
                for (size_t j = BED_STRAND_FIELD + 1; j < fields.size(); j++) {
                    out << '\t' << fields[j];
                }
                out << '\t' << feature_id(occurrence.feature) << '\n';
            }
            written += occurrences.size();
        });
        if (progress) {
            cerr << "[vg lift] " << written << " places on paths" << endl;
        }
    }

    return 0;
}

// Register subcommand
static Subcommand vg_lift("lift", "lift features found on a unitig cover to the graph, gref and assemblies", WIDGET, main_lift);
