/** \file hub_eval_main.cpp
 *
 * Defines the "vg hub-eval" subcommand, which compares hub labeling methods on
 * the oversized snarls of a graph without building a distance index.
 */

#include <getopt.h>

#include <iostream>
#include <string>
#include <unordered_map>

#include "subcommand.hpp"

#include "../utility.hpp"
#include "../integrated_snarl_finder.hpp"
#include "../snarl_distance_index.hpp"
#include "../gbwtgraph_helper.hpp"
#include <vg/io/vpkg.hpp>
#include <bdsg/hub_label_eval.hpp>

using namespace std;
using namespace vg;
using namespace vg::subcommand;

void help_hub_eval(char** argv) {
    std::cerr << "usage: " << argv[0] << " hub-eval [options] <graph.gbz> >report.tsv" << endl
         << "Build hub labels for each oversized snarl with each labeling method," << endl
         << "check them against Dijkstra, and report sizes and times. No index is written." << endl
         << "options:" << endl
         << "  -l, --snarl-limit N         snarls with more than N nodes are oversized [50000]" << endl
         << "  -m, --methods LIST          comma-separated labeling methods [all]" << endl
         << "  -M, --min-vertices N        skip snarls whose overlay has fewer than N vertices [0]" << endl
         << "  -F, --full-verify-limit N   check all pairs in overlays up to N vertices [2000]" << endl
         << "  -V, --verify-sources N      else check N random sources against all targets [50]" << endl
         << "  -q, --query-pairs N         random queries to time per snarl [100000]" << endl
         << "  -P, --path-prefix NAME      upweight tips of reference paths with this prefix," << endl
         << "                              as vg index does, so snarls come out the same" << endl
         << "  -p, --progress              show progress" << endl
         << "  -h, --help                  print this help message to stderr and exit" << endl
         << "methods:";
    for (const auto& method : default_labeling_methods()) {
        std::cerr << " " << method.name;
    }
    std::cerr << endl;
}

int main_hub_eval(int argc, char** argv) {
    Logger logger("vg hub-eval");
    bool show_progress = false;
    size_t snarl_limit = 50000;
    size_t min_vertices = 0;
    string method_names;
    string ref_prefix;
    EvalOptions eval_options;

    int c;
    optind = 2; // force optind past command positional argument
    while (true) {
        static struct option long_options[] =
            {
                {"snarl-limit",       required_argument, 0, 'l'},
                {"methods",           required_argument, 0, 'm'},
                {"min-vertices",      required_argument, 0, 'M'},
                {"full-verify-limit", required_argument, 0, 'F'},
                {"verify-sources",    required_argument, 0, 'V'},
                {"query-pairs",       required_argument, 0, 'q'},
                {"path-prefix",       required_argument, 0, 'P'},
                {"progress",          no_argument,       0, 'p'},
                {"help",              no_argument,       0, 'h'},
                {0, 0, 0, 0}
            };

        int option_index = 0;
        c = getopt_long(argc, argv, "l:m:M:F:V:q:P:ph?",
                        long_options, &option_index);

        if (c == -1)
            break;

        switch (c) {
        case 'l':
            snarl_limit = parse<size_t>(optarg);
            break;
        case 'm':
            method_names = optarg;
            break;
        case 'M':
            min_vertices = parse<size_t>(optarg);
            break;
        case 'F':
            eval_options.full_verify_limit = parse<size_t>(optarg);
            break;
        case 'V':
            eval_options.verify_sources = parse<size_t>(optarg);
            break;
        case 'q':
            eval_options.query_pairs = parse<size_t>(optarg);
            break;
        case 'P':
            ref_prefix = optarg;
            break;
        case 'p':
            show_progress = true;
            break;
        case 'h':
        case '?':
            help_hub_eval(argv);
            exit(1);
            break;
        default:
            abort();
        }
    }

    if (optind != argc - 1) {
        help_hub_eval(argv);
        exit(1);
    }
    string graph_name = require_exists(logger, argv[optind]);

    vector<LabelingMethod> methods;
    try {
        methods = method_names.empty() ? default_labeling_methods() : labeling_methods_by_name(method_names);
    } catch (const std::invalid_argument& e) {
        logger.error() << e.what() << endl;
    }

    // May be GBZ or a HandleGraph.
    auto options = vg::io::VPKG::try_load_first<gbwtgraph::GBZ, handlegraph::HandleGraph>(graph_name);
    const HandleGraph* graph = nullptr;
    if (get<0>(options)) {
        graph = &(get<0>(options)->graph);
    } else if (get<1>(options)) {
        graph = get<1>(options).get();
    } else {
        logger.error() << "input is not a graph or GBZ" << endl;
    }

    // Upweight the tips of reference paths, as vg index does.
    constexpr size_t EXTRA_WEIGHT = 10000000000;
    std::unordered_map<nid_t, size_t> extra_node_weight;
    if (!ref_prefix.empty()) {
        const PathHandleGraph* path_graph = dynamic_cast<const PathHandleGraph*>(graph);
        if (path_graph == nullptr) {
            logger.error() << "-P cannot be used because graph format does not support paths" << endl;
        }
        path_graph->for_each_path_of_sense({PathSense::REFERENCE, PathSense::GENERIC}, [&](const path_handle_t& path_handle) {
            string path_name = path_graph->get_path_name(path_handle);
            if (path_name.compare(0, ref_prefix.size(), ref_prefix) == 0 && !path_graph->is_empty(path_handle)) {
                extra_node_weight[path_graph->get_id(path_graph->get_handle_of_step(path_graph->path_begin(path_handle)))] += EXTRA_WEIGHT;
                extra_node_weight[path_graph->get_id(path_graph->get_handle_of_step(path_graph->path_back(path_handle)))] += EXTRA_WEIGHT;
            }
        });
    }

    IntegratedSnarlFinder snarl_finder(*graph, extra_node_weight);

    bool header = true;
    size_t evaluated = 0;
    size_t skipped = 0;
    size_t mismatches = 0;
    OversizedSnarlCallback evaluate_snarl = [&](const CHOverlay& ov, const SnarlDistanceIndex::TemporaryDistanceIndex::TemporarySnarlRecord& record) {
        if (num_vertices(ov) < min_vertices) {
            skipped++;
            return;
        }
        string label = record.is_root_snarl ? "root" :
            to_string(record.start_node_id) + (record.start_node_rev ? "-" : "+") + "_" +
            to_string(record.end_node_id) + (record.end_node_rev ? "-" : "+");
        if (show_progress) {
            logger.info() << "Evaluating snarl " << label << " with " << num_vertices(ov) << " overlay vertices" << endl;
        }
        auto results = evaluate(ov, methods, eval_options);
        write_tsv(cout, label, results, header);
        cout.flush();
        header = false;
        evaluated++;
        for (const auto& result : results) {
            mismatches += result.mismatches;
        }
    };

    // This still builds the default labels for each oversized snarl, which the
    // temporary index needs, but never makes or writes the final index.
    make_temporary_distance_index(graph, &snarl_finder, snarl_limit, false, HubOrdering::CH, evaluate_snarl);

    if (show_progress) {
        logger.info() << "Evaluated " << evaluated << " oversized snarls, skipped " << skipped << endl;
    }
    if (mismatches > 0) {
        logger.error() << mismatches << " hub label distances did not match Dijkstra" << endl;
    }
    return 0;
}

// Register subcommand
static Subcommand vg_hub_eval("hub-eval", "compare hub labeling methods on oversized snarls", DEVELOPMENT, main_hub_eval);
