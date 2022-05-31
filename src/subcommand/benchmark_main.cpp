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
#include "../integrated_snarl_finder.hpp"
#include "../min_distance.hpp"

#include <bdsg/internal/mapped_structs.hpp>



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
    
    char filename[] = "tmpXXXXXX";
    int tmpfd = mkstemp(filename);
    assert(tmpfd != -1);
    
    {
        // Make all the vectors to benchmark
        sdsl::int_vector<> sdsl_vec_holder;
        auto sdsl_vec = &sdsl_vec_holder;
        bdsg::CompatIntVector<> compat_vec_holder;
        auto compat_vec = &compat_vec_holder;
        bdsg::yomo::UniqueMappedPointer<bdsg::MappedIntVector> mapped_vec;
        mapped_vec.construct("");
        mapped_vec.save(tmpfd);
        
        size_t vector_size = 1024 * 100;
        size_t access_iters = 1000;
        
        sdsl_vec->width(24);
        sdsl_vec->resize(vector_size);
        for (size_t i = 0; i < sdsl_vec->size(); i++) {
            (*sdsl_vec)[i] = i % 1<<20;
        }
        compat_vec->width(24);
        compat_vec->resize(vector_size);
        for (size_t i = 0; i < compat_vec->size(); i++) {
            (*compat_vec)[i] = i % 1<<20;
        }
        mapped_vec->width(24);
        mapped_vec->resize(vector_size);
        for (size_t i = 0; i < sdsl_vec->size(); i++) {
            (*mapped_vec)[i] = i % 1<<20;
        }
        
        size_t bits = 1;
        
        results.push_back(run_benchmark("SDSL random access", 1000, [&]() {
            for (size_t i = 1; i < access_iters; i++) {
                size_t pos = bits % vector_size;
                size_t read = (*sdsl_vec)[pos];
                if (read != pos % 1<<20) {
                    throw std::runtime_error("Incorrect read!");
                }
                bits = bits ^ (bits << 13) ^ i;
            }
        }));
        
        results.push_back(run_benchmark("SDSL sequential access", 1000, [&]() {
            for (size_t i = 0; i < std::min(vector_size, access_iters); i++) {
                size_t read = (*sdsl_vec)[i];
                if (read != i % 1<<20) {
                    throw std::runtime_error("Incorrect read!");
                }
                bits = bits ^ (bits << 13) ^ i;
            }
        }));
        
        results.push_back(run_benchmark("Compat random access", 1000, [&]() {
            for (size_t i = 1; i < access_iters; i++) {
                size_t pos = bits % vector_size;
                size_t read = (*compat_vec)[pos];
                if (read != pos % 1<<20) {
                    throw std::runtime_error("Incorrect read!");
                }
                bits = bits ^ (bits << 13) ^ i;
            }
        }));
        
        results.push_back(run_benchmark("Compat sequential access", 1000, [&]() {
            for (size_t i = 0; i < std::min(vector_size, access_iters); i++) {
                size_t read = (*compat_vec)[i];
                if (read != i % 1<<20) {
                    throw std::runtime_error("Incorrect read!");
                }
                bits = bits ^ (bits << 13) ^ i;
            }
        }));
        
        results.push_back(run_benchmark("Mapped random access", 1000, [&]() {
            for (size_t i = 1; i < access_iters; i++) {
                size_t pos = bits % vector_size;
                size_t read = (*mapped_vec)[pos];
                if (read != pos % 1<<20) {
                    throw std::runtime_error("Incorrect read!");
                }
                bits = bits ^ (bits << 13) ^ i;
            }
        }));
        
        results.push_back(run_benchmark("Mapped sequential access", 1000, [&]() {
            for (size_t i = 0; i < std::min(vector_size, access_iters); i++) {
                size_t read = (*mapped_vec)[i];
                if (read != i % 1<<20) {
                    throw std::runtime_error("Incorrect read!");
                }
                bits = bits ^ (bits << 13) ^ i;
            }
        }));
        
        // Reloading the mapped vector will get it all into the same chain link
        // and make sure we don't need to consult global tables on pointer
        // access.
        mapped_vec.reset();
        mapped_vec.load(tmpfd, "");
        
        results.push_back(run_benchmark("Mapped reloaded random access", 1000, [&]() {
            for (size_t i = 1; i < access_iters; i++) {
                size_t pos = bits % vector_size;
                size_t read = (*mapped_vec)[pos];
                if (read != pos % 1<<20) {
                    throw std::runtime_error("Incorrect read!");
                }
                bits = bits ^ (bits << 13) ^ i;
            }
        }));
        
        results.push_back(run_benchmark("Mapped reloaded sequential access", 1000, [&]() {
            for (size_t i = 0; i < std::min(vector_size, access_iters); i++) {
                size_t read = (*mapped_vec)[i];
                if (read != i % 1<<20) {
                    throw std::runtime_error("Incorrect read!");
                }
            }
        }));
    }
    
    close(tmpfd);
    unlink(filename);
    
    
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

