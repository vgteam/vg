/** \file benchmark_main.cpp
 *
 * Defines the "vg benchmark" subcommand, which runs and reports on microbenchmarks.
 */

#include <omp.h>
#include <unistd.h>
#include <getopt.h>

#include <iostream>

#include "subcommand.hpp"

#include "../benchmark.hpp"
#include "../version.hpp"

#include "../vg.hpp"
#include "../xg.hpp"
#include "../algorithms/extract_connecting_graph.hpp"



using namespace std;
using namespace vg;
using namespace vg::subcommand;

void help_benchmark(char** argv) {
    cerr << "usage: " << argv[0] << " benchmark [options] >report.tsv" << endl
         << "options:" << endl
         << "    -p, --progress         show progress" << endl;
}

int main_benchmark(int argc, char** argv) {

    bool show_progress = false;
    
    int c;
    optind = 2; // force optind past command positional argument
    while (true) {
        static struct option long_options[] =
            {
                {"progress",  no_argument, 0, 'p'},
                {"help", no_argument, 0, 'h'},
                {0, 0, 0, 0}
            };

        int option_index = 0;
        c = getopt_long (argc, argv, "ph?",
                         long_options, &option_index);

        /* Detect the end of the options. */
        if (c == -1)
            break;

        switch (c)
        {

        case 'p':
            show_progress = true;
            break;
            
        case 'h':
        case '?':
            /* getopt_long already printed an error message. */
            help_benchmark(argv);
            exit(1);
            break;

        default:
            abort ();

        }
    }
    
    if (optind != argc) {
        // Extra arguments found
        help_benchmark(argv);
        exit(1);
    }
    
    // Do all benchmarking on one thread
    omp_set_num_threads(1);
    
    // Turn on nested parallelism, so we can parallelize over VCFs and over alignment bands
    omp_set_nested(1);
    
    // Generate a test graph
    VG vg;
    for (size_t i = 1; i < 101; i++) {
        // It will have 100 nodes
        vg.create_node("ACGTACGT", i);
    }
    size_t bits = 1;
    for (size_t i = 1; i < 101; i++) {
        for (size_t j = 1; j < 101; j++) {
            if ((bits ^ (i + (j << 3))) % 50 == 0) {
                // Make some arbitrary edges
                vg.create_edge(i, j, false, false);
            }
            // Shifts and xors make good PRNGs right?
            bits = bits ^ (bits << 13) ^ j;            
        }
    }
    
    // And a test XG of it
    xg::XG xg_index(vg.graph);
    
    vector<BenchmarkResult> results;
    
    // Do the control against itself
    results.push_back(run_benchmark("control", 1000, benchmark_control));
    
    results.push_back(run_benchmark("VG::get_node", 1000, [&]() {
        for (size_t rep = 0; rep < 100; rep++) {
            for (size_t i = 1; i < 101; i++) {
                vg.get_node(i);
            }
        }
    }));
    
    
    results.push_back(run_benchmark("algorithms::extract_connecting_graph on xg", 1000, [&]() {
        pos_t pos_1 = make_pos_t(55, false, 0);
        pos_t pos_2 = make_pos_t(32, false, 0);
        
        int64_t max_len = 5;
        
        Graph g;
        
        auto trans = algorithms::extract_connecting_graph(&xg_index, g, max_len, pos_1, pos_2, false, false, true, true, true);
    
    }));
    
    results.push_back(run_benchmark("algorithms::extract_connecting_graph on vg", 1000, [&]() {
        pos_t pos_1 = make_pos_t(55, false, 0);
        pos_t pos_2 = make_pos_t(32, false, 0);
        
        int64_t max_len = 5;
        
        Graph g;
        
        auto trans = algorithms::extract_connecting_graph(vg, g, max_len, pos_1, pos_2, false, false, true, true, true);
    
    }));
    
                

    cout << "# Benchmark results for vg " << VG_VERSION_STRING << endl;
    cout << "# runs\ttest(us)\tstddev(us)\tcontrol(us)\tstddev(us)\tscore\terr\tname" << endl;
    for (auto& result : results) {
        cout << result << endl;
    }

    return 0;
}

// Register subcommand
static Subcommand vg_benchmark("benchmark", "run and report on performance benchmarks", main_benchmark);

