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

#include "../gbwt_extender.hpp"
#include "../gbwt_helper.hpp"



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
    
    size_t node_count = 10;
    size_t node_length = 32;
    
    // Prepare a GBWT of one long path
    std::vector<gbwt::vector_type> paths;
    paths.emplace_back();
    for (size_t i = 0; i < node_count; i++) {
        paths.back().push_back(gbwt::Node::encode(i + 1, false));
    }
    gbwt::GBWT index = get_gbwt(paths);
    
    // Turn it into a GBWTGraph.
    // Make a SequenceSource we will consult later for getting sequence.
    gbwtgraph::SequenceSource source;
    uint32_t bits = 0xcafebebe;
    auto step_rng = [&bits]() {
        // Try out <https://stackoverflow.com/a/69142783>
        bits = (bits * 73 + 1375) % 477218579;
    };
    for (size_t i = 0; i < node_count; i++) {
        std::stringstream ss;
        for (size_t j = 0; j < node_length; j++) {
            // Pick a deterministic character
            ss << "ACGT"[bits & 0x3];
            step_rng();
        }
        source.add_node(i + 1, ss.str());
    }
    // And then make the graph
    gbwtgraph::GBWTGraph graph(index, source);
    
    // Decide what we are going to align
    pos_t from_pos = make_pos_t(1, false, 3);
    pos_t to_pos = make_pos_t(node_count, false, 11);
    
    // Synthesize a sequence
    std::stringstream seq_stream;
    seq_stream << source.get_sequence(get_id(from_pos)).substr(get_offset(from_pos) + 1);
    for (nid_t i = get_id(from_pos) + 1; i < get_id(to_pos); i++) {
        std::string seq = source.get_sequence(i);
        // Add some errors
        if (bits & 0x1) {
            int offset = bits % seq.size();
            step_rng();
            char replacement = "ACGT"[bits & 0x3];
            step_rng();
            if (bits & 0x1) {
                seq[offset] = replacement;
            } else {
                step_rng();
                if (bits & 0x1) {
                    seq.insert(offset, 1, replacement);
                } else {
                    seq.erase(offset);
                }
            }
        }
        step_rng();
        // And keep the sequence
        seq_stream << seq;
    }
    seq_stream << source.get_sequence(get_id(to_pos)).substr(0, get_offset(to_pos)); 
    
    std::string to_connect = seq_stream.str();
    
    // Make the Aligner and Extender
    Aligner aligner;
    WFAExtender extender(graph, aligner);
    
    results.push_back(run_benchmark("connect() on long sequence", 1000, [&]() {
        // Do the alignment
        WFAAlignment aligned = extender.connect(to_connect, from_pos, to_pos);
        // Make sure it succeeded
        assert(aligned);
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

