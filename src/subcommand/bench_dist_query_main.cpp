/** \file bench_dist_query_main.cpp
 *
 * Defines the "vg bench-dist-query" subcommand, which benchmarks distance query speed across multiple indexes.
 */

#include <omp.h>
#include <unistd.h>
#include <getopt.h>

#include <iostream>
#include <fstream>
#include <filesystem>
#include <sstream>
#include <chrono>

#include "subcommand.hpp"

#include "../benchmark.hpp"
#include "../version.hpp"

#include "../snarl_distance_index.hpp"
#include "../integrated_snarl_finder.hpp"
#include "../algorithms/gfa_to_handle.hpp"
#include <vg/io/vpkg.hpp>
#include <bdsg/hash_graph.hpp>
#include "../gbwtgraph_helper.hpp"



using namespace std;
using namespace vg;
using namespace vg::subcommand;

void help_bench_dist_query(char** argv) {
    std::cerr << "usage: " << argv[0] << " bench-dist-query -g <graph.gbz> -d <index1.dist> [-d <index2.dist> ...] [options] >report.tsv" << endl
         << "options:" << endl
         << "  -g, --graph FILE         path to input GBZ graph file" << endl
         << "  -d, --dist FILE          path to distance index file (repeatable)" << endl
         << "  -q, --numQueries N       number of queries to run (default: 10000)" << endl
         << "  -s, --save-queries FILE  save generated queries to FILE for reproducibility" << endl
         << "  -Q, --load-queries FILE  load queries from FILE instead of generating new ones" << endl
         << "  -p, --progress           show progress" << endl
         << "  -h, --help               print this help message to stderr and exit" << endl;
}


int main_bench_dist_query(int argc, char** argv) {
    Logger logger("vg bench-dist-query");
    bool show_progress = false;

    string graph_path = "";
    vector<string> dist_paths;
    int num_queries = 10000;
    string save_queries_path = "";
    string load_queries_path = "";

    int c;
    optind = 2; // force optind past command positional argument
    while (true) {
        static struct option long_options[] =
            {
                {"progress",      no_argument,       0, 'p'},
                {"help",          no_argument,       0, 'h'},
                {"graph",         required_argument, 0, 'g'},
                {"dist",          required_argument, 0, 'd'},
                {"numQueries",    required_argument, 0, 'q'},
                {"save-queries",  required_argument, 0, 's'},
                {"load-queries",  required_argument, 0, 'Q'},
                {0, 0, 0, 0}
            };

        int option_index = 0;
        c = getopt_long(argc, argv, "g:d:q:s:Q:ph?",
                        long_options, &option_index);

        if (c == -1)
            break;

        switch (c) {
        case 'g':
            graph_path = require_exists(logger, optarg);
            break;
        case 'd':
            dist_paths.push_back(require_exists(logger, optarg));
            break;
        case 'q':
            {
                num_queries = stoi(optarg);
            }
            break;
        case 's':
            save_queries_path = ensure_writable(logger, optarg);
            break;
        case 'Q':
            load_queries_path = require_exists(logger, optarg);
            break;
        case 'p':
            show_progress = true;
            break;
        case 'h':
        case '?':
            help_bench_dist_query(argv);
            exit(1);
            break;
        default:
            abort();
        }
    }

    if (graph_path.empty()) {
        logger.error() << "a GBZ graph file is required (-g)" << endl;
    }

    if (dist_paths.empty()) {
        logger.error() << "at least one distance index file is required (-d)" << endl;
    }

    // Load GBZ graph
    if (show_progress) {
        logger.info() << "Loading GBZ graph from " << graph_path << "..." << endl;
    }
    gbwtgraph::GBZ gbz;
    load_gbz(gbz, graph_path, show_progress);
    const HandleGraph& graph = gbz.graph;
    logger.info() << "Loaded graph with " << graph.get_node_count() << " nodes" << endl;

    // Collect all node IDs
    vector<nid_t> all_node_ids;
    graph.for_each_handle([&](handle_t h) {
        all_node_ids.push_back(graph.get_id(h));
    });

    using QueryEntry = pair<tuple<nid_t, bool, size_t>, tuple<nid_t, bool, size_t>>;
    vector<QueryEntry> queries;

    if (!load_queries_path.empty()) {
        if (show_progress) {
            logger.info() << "Loading queries from " << load_queries_path << "..." << endl;
        }
        ifstream qf(load_queries_path);
        string line;
        while (getline(qf, line)) {
            if (line.empty()) continue;
            istringstream iss(line);
            nid_t id1, id2; int rev1, rev2; size_t off1, off2;
            iss >> id1 >> rev1 >> off1 >> id2 >> rev2 >> off2;
            QueryEntry q;
            q.first  = make_tuple(id1, (bool)rev1, off1);
            q.second = make_tuple(id2, (bool)rev2, off2);
            queries.push_back(q);
        }
        logger.info() << "Loaded " << queries.size() << " queries from " << load_queries_path << endl;
    } else {
        if (show_progress) {
            logger.info() << "Generating " << num_queries << " queries..." << endl;
        }
        queries.resize(num_queries);
        for (auto& query : queries) {
            nid_t node1 = all_node_ids[rand() % all_node_ids.size()];
            nid_t node2 = all_node_ids[rand() % all_node_ids.size()];
            size_t len1 = graph.get_length(graph.get_handle(node1));
            size_t len2 = graph.get_length(graph.get_handle(node2));
            query.first  = make_tuple(node1, rand() % 2 == 1, len1 > 0 ? rand() % len1 : 0);
            query.second = make_tuple(node2, rand() % 2 == 1, len2 > 0 ? rand() % len2 : 0);
        }
        logger.info() << "Generated " << queries.size() << " queries" << endl;
    }

    if (!save_queries_path.empty()) {
        ofstream qf(save_queries_path);
        for (auto& query : queries) {
            auto& [id1, rev1, off1] = query.first;
            auto& [id2, rev2, off2] = query.second;
            qf << id1 << "\t" << (int)rev1 << "\t" << off1 << "\t"
               << id2 << "\t" << (int)rev2 << "\t" << off2 << "\n";
        }
        logger.info() << "Saved " << queries.size() << " queries to " << save_queries_path << endl;
    }

    // Output header
    cout << "dist_index\tavg_query_us" << endl;

    // Benchmark each distance index
    for (const auto& dist_path : dist_paths) {
        if (show_progress) {
            logger.info() << "Loading distance index from " << dist_path << "..." << endl;
        }
        SnarlDistanceIndex distance_index;
        distance_index.deserialize(dist_path);
        logger.info() << "Loaded distance index from " << dist_path << endl;

        // Pull the whole index into the OS page cache so timings reflect
        // unavoidable query cost, not avoidable first-touch I/O
        distance_index.preload(true);

        // Time all queries
        auto start = chrono::high_resolution_clock::now();
        for (auto& query : queries) {
            auto& [node1, node1_rev, node1_offset] = query.first;
            auto& [node2, node2_rev, node2_offset] = query.second;
            distance_index.minimum_distance(
                node1, node1_rev, node1_offset,
                node2, node2_rev, node2_offset,
                false,
                nullptr
            );
        }
        auto end = chrono::high_resolution_clock::now();

        double total_us = chrono::duration<double, micro>(end - start).count();
        double avg_us = total_us / queries.size();

        filesystem::path dist_fs_path(dist_path);
        cout << dist_fs_path.filename().string() << "\t" << avg_us << endl;
        logger.info() << dist_path << ": avg query time = " << avg_us << " us" << endl;
    }

    return 0;
}

// Register subcommand
static Subcommand vg_bench_dist_query("bench-dist-query", "benchmark distance query speed across multiple indexes", DEVELOPMENT, main_bench_dist_query);
