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
#include "../indexed_vg.hpp"
#include "../algorithms/extract_connecting_graph.hpp"
#include "../algorithms/topological_sort.hpp"
#include "../algorithms/weakly_connected_components.hpp"



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
    
    // Which experiments should we run?
    bool sort_and_order_experiment = false;
    bool get_sequence_experiment = true;
    
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
    
    
    
    vector<BenchmarkResult> results;
    
    // Generate a test graph
    VG vg_mut;
    for (size_t i = 1; i < 101; i++) {
        // It will have 100 nodes in ID order
        vg_mut.create_node("ACGTACGT", i);
    }
    size_t bits = 1;
    for (size_t i = 1; i < 101; i++) {
        for (size_t j = 1; j < 101; j++) {
            if ((bits ^ (i + (j << 3))) % 50 == 0) {
                // Make some arbitrary edges
                vg_mut.create_edge(i, j, false, false);
            }
            // Shifts and xors make good PRNGs right?
            bits = bits ^ (bits << 13) ^ j;            
        }
    }
    
    // Save a constant copy to reset it with
    const VG vg(vg_mut);
    
    // Get an XG for it
    const xg::XG xg_index(vg_mut.graph);
    
    // And an IndexedVG
    auto filename = temp_file::create();
    vg_mut.serialize_to_file(filename, 10);
    const IndexedVG indexed_vg(filename);
    
    if (sort_and_order_experiment) {
    
        results.push_back(run_benchmark("vg::algorithms topological_order", 1000, [&]() {
            vector<handle_t> order = algorithms::topological_order(&vg);
            assert(order.size() == vg.node_size());
        }));
        
        results.push_back(run_benchmark("vg::algorithms topological_sort", 1000, [&]() {
            vg_mut = vg;
        }, [&]() {
            algorithms::topological_sort(&vg_mut);
        }));
        
        results.push_back(run_benchmark("vg::algorithms orient_nodes_forward", 1000, [&]() {
            vg_mut = vg;
        }, [&]() {
            algorithms::orient_nodes_forward(&vg_mut);
        }));
        
        results.push_back(run_benchmark("vg::algorithms weakly_connected_components", 1000, [&]() {
            auto components = algorithms::weakly_connected_components(&vg);
            assert(components.size() == 1);
            assert(components.front().size() == vg.node_size());
        }));
        
    }
    
    if (get_sequence_experiment) {
    
        results.push_back(run_benchmark("VG::get_sequence", 1000, [&]() {
            for (size_t i = 1; i < 101; i++) {
                handle_t handle = vg.get_handle(i);
                string sequence = vg.get_sequence(handle);
            }
        }));
        
        results.push_back(run_benchmark("XG::get_sequence", 1000, [&]() {
            for (size_t i = 1; i < 101; i++) {
                handle_t handle = xg_index.get_handle(i);
                string sequence = xg_index.get_sequence(handle);
            }
        }));
        
        results.push_back(run_benchmark("IndexedVG::get_sequence", 1000, [&]() {
            for (size_t i = 1; i < 101; i++) {
                handle_t handle = indexed_vg.get_handle(i);
                string sequence = indexed_vg.get_sequence(handle);
            }
        }));
        
    }
    
    // Do the control against itself
    results.push_back(run_benchmark("control", 1000, benchmark_control));

    cout << "# Benchmark results for vg " << Version::get_short() << endl;
    cout << "# runs\ttest(us)\tstddev(us)\tcontrol(us)\tstddev(us)\tscore\terr\tname" << endl;
    for (auto& result : results) {
        cout << result << endl;
    }
    
    return 0;
}

// Register subcommand
static Subcommand vg_benchmark("benchmark", "run and report on performance benchmarks", DEVELOPMENT, main_benchmark);

