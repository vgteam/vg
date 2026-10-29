/**
 * \file lift_main.cpp: lift features on a cover or on paths to the graph, the gref and the paths
 */

#include <getopt.h>
#include <omp.h>

#include <fstream>
#include <iostream>
#include <unordered_set>

#include "subcommand.hpp"
#include "../feature_lift.hpp"
#include "../gbwtgraph_helper.hpp"
#include "../gref.hpp"
#include "../utility.hpp"

using namespace std;
using namespace vg;
using namespace vg::subcommand;

void help_lift(char** argv) {
    cerr << "usage: " << argv[0] << " lift [options] graph.gbz features.bed > features.gaf" << endl
         << "Lift features to the graph and the gref, and to each path that visits them." << endl
         << "A feature can be on a cover sequence (a gref path of the GBZ or a tig of" << endl
         << "vg unitigs) or on any path in the GBZ, named without a subrange or phase" << endl
         << "block and in coordinates of the full path (for haplotypes, the assembly" << endl
         << "contig). A feature on several sequences at the same graph place is reported" << endl
         << "once." << endl
         << endl
         << "Each feature gets an ID (f1, f2, ...) and a GAF record with that query name," << endl
         << "its graph walk, and the tags sp:i: (visits of the whole walk by paths other" << endl
         << "than gref paths), gr:Z: (gref placement, see vg unitigs) and an:Z: (BED name)." << endl
         << endl
         << "options:" << endl
         << "  -u, --tigs FILE         GAF of the tigs, from vg unitigs -g" << endl
         << "  -a, --assemblies FILE   write each place where a path visits a feature to FILE" << endl
         << "                          as BED on the full path: the feature's name, score" << endl
         << "                          and strand (., 0 and . if it has none), any further" << endl
         << "                          BED columns, and then its ID" << endl
         << "  -m, --min-match F       write a place wherever a path visits at least F of the" << endl
         << "                          feature's bases in order (liftover), and end each row" << endl
         << "                          with that fraction. Without -m, a place must hold the" << endl
         << "                          feature's whole walk. A gap on the path longer than" << endl
         << "                          the feature splits a place." << endl
         << "  -b, --blocks            write places as BED12, with a block for each stretch" << endl
         << "                          that holds feature bases" << endl
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
    double min_match = 0;
    bool blocks = false;
    bool progress = false;

    int c;
    optind = 2; // force optind past command positional argument
    while (true) {
        static struct option long_options[] = {
            {"tigs", required_argument, 0, 'u'},
            {"assemblies", required_argument, 0, 'a'},
            {"min-match", required_argument, 0, 'm'},
            {"blocks", no_argument, 0, 'b'},
            {"paths-by", required_argument, 0, 'Q'},
            {"threads", required_argument, 0, 't'},
            {"progress", no_argument, 0, 'p'},
            {"help", no_argument, 0, 'h'},
            {0, 0, 0, 0}
        };
        int option_index = 0;
        c = getopt_long(argc, argv, "u:a:m:bQ:t:ph?", long_options, &option_index);
        if (c == -1) {
            break;
        }
        switch (c) {
        case 'u':
            tigs_file = require_exists(logger, optarg);
            break;
        case 'a':
            assemblies_file = ensure_writable(logger, optarg);
            break;
        case 'm':
            min_match = parse<double>(optarg);
            if (min_match <= 0 || min_match > 1) {
                logger.error() << "-m/--min-match must be in (0, 1]" << endl;
            }
            break;
        case 'b':
            blocks = true;
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
    const GrefIndex* gref_pointer = gref.path_count() > 0 ? &gref : nullptr;
    if (progress) {
        cerr << "[vg lift] gref: " << gref.path_count() << " paths, " << gref.total_length() << " bp" << endl;
    }
    Cover cover(graph, gref_pointer);

    vector<BedFeature> bed_features;
    get_input_file(features_file, [&](istream& in) {
        try {
            bed_features = read_bed_features(in);
        } catch (const runtime_error& e) {
            logger.error() << e.what() << endl;
        }
    });
    if (progress) {
        cerr << "[vg lift] " << bed_features.size() << " features" << endl;
    }
    for (const BedFeature& feature : bed_features) {
        if (feature.fields.size() >= 9) {
            logger.warn() << "BED12 blocks are not lifted; columns from the 7th on are copied as they are" << endl;
            break;
        }
    }

    // Find the sequence of each feature: a gref path, a tig, or a path in the
    // GBZ. Load only the tigs and paths that features are on.
    unordered_set<string> wanted;
    for (const BedFeature& feature : bed_features) {
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
    PathDirectory paths(gbz, [&](const string& name) {
        return wanted.count(name) > 0 && !cover.sequence_length(name);
    });
    for (const BedFeature& feature : bed_features) {
        optional<size_t> length = cover.sequence_length(feature.sequence);
        if (length && feature.end > *length) {
            logger.error() << "feature " << feature.sequence << ":" << feature.start << "-" << feature.end
                           << " ends past the end of the sequence (" << *length << " bp)" << endl;
        }
        if (!length && !paths.has_path(feature.sequence)) {
            logger.error() << "feature on " << feature.sequence << ", which is not a gref path or a path in "
                           << graph_file << ", or a tig" << (tigs_file.empty() ? " (use -u)" : " in " + tigs_file) << endl;
        }
    }

    // A feature on a path lifts only if one fragment of the path holds it.
    vector<GraphFeature> features(bed_features.size());
    vector<bool> lifted(bed_features.size(), true);
    #pragma omp parallel for schedule(dynamic, 1024)
    for (size_t i = 0; i < bed_features.size(); i++) {
        BedFeature& feature = bed_features[i];
        optional<WalkInterval> interval;
        if (cover.sequence_length(feature.sequence)) {
            interval = cover.interval(feature.sequence, feature.start, feature.end);
        } else {
            interval = paths.interval(feature.sequence, feature.start, feature.end);
        }
        if (interval) {
            features[i] = make_graph_feature(graph, gref_pointer, std::move(*interval), std::move(feature.fields));
        } else {
            lifted[i] = false;
        }
    }
    size_t outside = 0;
    for (size_t i = 0; i < features.size(); i++) {
        if (!lifted[i]) {
            outside++;
        } else if (outside > 0) {
            features[i - outside] = std::move(features[i]);
        }
    }
    features.resize(features.size() - outside);
    if (outside > 0) {
        logger.warn() << outside << " features on paths are not within one fragment of the path in "
                      << graph_file << " and were skipped" << endl;
    }
    bed_features = vector<BedFeature>();
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

    if (assemblies_file.empty()) {
        return 0;
    }
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
    auto write_places = [&](const string& name, const vector<PathOccurrence>& places) {
        for (const PathOccurrence& place : places) {
            const GraphFeature& feature = features[place.feature];
            const vector<string>& fields = feature.fields;
            out << name << '\t' << place.start << '\t' << place.end
                << '\t' << (fields.size() > 0 ? fields[0] : ".")
                << '\t' << (fields.size() > 1 ? fields[1] : "0") << '\t';
            if (feature.is_stranded()) {
                out << (place.is_reverse ? '-' : '+');
            } else {
                out << '.';
            }
            if (blocks) {
                out << '\t' << place.start << '\t' << place.end << "\t0\t" << place.blocks.size() << '\t';
                for (const auto& block : place.blocks) {
                    out << (block.second - block.first) << ',';
                }
                out << '\t';
                for (const auto& block : place.blocks) {
                    out << (block.first - place.start) << ',';
                }
            }
            for (size_t j = BED_STRAND_FIELD + 1; j < fields.size(); j++) {
                out << '\t' << fields[j];
            }
            out << '\t' << feature_id(place.feature);
            if (min_match > 0) {
                out << '\t' << double(place.matched) / feature.interval.length;
            }
            out << '\n';
        }
        written += places.size();
    };
    if (min_match > 0) {
        for_each_path_liftover(gbz, features, min_match, include, write_places);
    } else {
        for_each_path_occurrences(gbz, features, include, write_places);
    }
    if (progress) {
        cerr << "[vg lift] " << written << " places on paths" << endl;
    }

    return 0;
}

// Register subcommand
static Subcommand vg_lift("lift", "lift features on a cover or on paths to the graph, gref and assemblies", WIDGET, main_lift);
