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

#include "../unittest/test_aligner.hpp"
#include "../vg.hpp"



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
    
    vg::unittest::TestAligner aligner_source;
    Aligner* aligner = (Aligner*) aligner_source.get_regular_aligner();
    
    vg::VG graph;

    vg::Node* n0 = graph.create_node("AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA");
    vg::Node* n1 = graph.create_node("CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC");
    vg::Node* n2 = graph.create_node("AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA");
    vg::Node* n3 = graph.create_node("TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT");
    vg::Node* n4 = graph.create_node("AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA");
   
    
    graph.create_edge(n0, n1);
    graph.create_edge(n0, n3);
    graph.create_edge(n1, n2);
    graph.create_edge(n3, n4);

    vg::Node* last = n4;
    for (size_t i = 0; i < 100; i++) {
        vg::Node* next = graph.create_node("GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG");
        graph.create_edge(last, next);
        last = next;
    }

    string read = string("AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAACCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA");
    Alignment aln;
    aln.set_sequence(read);
    
    vector<BenchmarkResult> results;
        
    results.push_back(run_benchmark("map against forking graph", 1000, [&]() {
        aligner->align_pinned(aln, graph, true, true, false); 
    }));

    results.push_back(run_benchmark("map against forking graph with node drop", 1000, [&]() {
        aligner->align_pinned(aln, graph, true, true, true); 
    }));
        
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

