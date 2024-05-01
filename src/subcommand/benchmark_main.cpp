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

#include <vg/io/json2pb.h>

#include <bdsg/hash_graph.hpp>
#include <vg/io/vpkg.hpp>



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
    
    vg::unittest::TestAligner aligner_source;
    const Aligner* aligner = aligner_source.get_regular_aligner();
    
    // Read the whole graph
    std::unique_ptr<HandleGraph> graph = vg::io::VPKG::load_one<HandleGraph>("test/alignment/pinned.vg");\
    assert(graph);

    // Read the whole read text.
    // See <https://stackoverflow.com/a/2912614>
    std::ifstream read_text_file("test/alignment/pinned.txt");
    std::string read_text((std::istreambuf_iterator<char>(read_text_file)), (std::istreambuf_iterator<char>()));
    while(!read_text.empty() && read_text.back() == '\n') {
        read_text.pop_back();
    }
    assert(!read_text.empty());
    
    vector<BenchmarkResult> results;

    Alignment aln;
    aln.set_sequence(read_text);

    /*results.push_back(run_benchmark("align to graph with node drop, 1k gap", 10, [&]() {
        aligner->align_pinned(aln, *graph, false, true, true, 1000);
    }));*/

    results.push_back(run_benchmark("align to graph with node drop, 9437 gap", 1, [&]() {
        aligner->align_pinned(aln, *graph, false, true, true, 9437);
    }));

    results.push_back(run_benchmark("align to graph with node drop, 9437 gap, again", 1, [&]() {
        aligner->align_pinned(aln, *graph, false, true, true, 9437);
    }));

    results.push_back(run_benchmark("align to graph with node drop, 9437 gap, repeatedly", 10, [&]() {
        aligner->align_pinned(aln, *graph, false, true, true, 9437);
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

