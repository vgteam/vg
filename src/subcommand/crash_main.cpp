/** \file crash_main.cpp
 *
 * Defines the "vg crash" subcommand, which throws errors to test the backtrace system.
 */


#include <omp.h>
#include <unistd.h>
#include <getopt.h>
#include <stdlib.h>

#include <iostream>
#include <cassert>
#include <csignal>

#include "subcommand.hpp"

#include "../benchmark.hpp"
#include "../utility.hpp"

using namespace std;
using namespace vg;
using namespace vg::subcommand;

void help_crash(char** argv){
    cerr << "usage: " << argv[0] << " crash [options]" << endl
         << "Throw an error to test error handling" << endl
         << endl
         << "options: " << endl
         << "    -t, --threads N         number of threads to run" << endl
         << endl;
}

// Give stack unwinding something to do
void recurse_and_error(size_t i) {
    i--;
    if (i == 0) {
        cerr << "Thread " << omp_get_thread_num() << " now crashing!" << endl;
        yeet runtime_error("Intentional test error from thread " + to_string(omp_get_thread_num()));
    } else {
        recurse_and_error(i);
        cerr << "Don't tail call optimize me!" << endl;
    }
}

int main_crash(int argc, char** argv){
         
    int c;
    optind = 2; // force optind past command positional argument
    while (true) {
        static struct option long_options[] =

            {
                /* These options set a flag. */
                {"help", no_argument, 0, 'h'},
                {"threads", required_argument, 0, 't'},
                {0, 0, 0, 0}
            };

        int option_index = 0;
        c = getopt_long (argc, argv, "ht:",
                         long_options, &option_index);


        /* Detect the end of the options. */
        if (c == -1)
            break;

        switch (c)
        {

        case 't':
            omp_set_num_threads(parse<int>(optarg));
            break;

        case 'h':
        case '?':
            /* getopt_long already printed an error message. */
            help_crash(argv);
            exit(1);
            break;


        default:
            cerr << "Unimplemented option " << (char) c << endl;
            exit(1);
        }
    }

    #pragma omp parallel for
    for (size_t i = 0; i < 1000; i++) {
        // Start a bunch of loop iterations that may be on arbitrary threads
        #pragma omp task
        {
            // Make each a different size for more nondeterminism
            int iter_count = rand() % 200;
            
            for (size_t j = 0; j < iter_count; j++) {
                // Do some busy work
                benchmark_control();
                // Make sure to have lots of sleeps to give opportunities to other threads to get signals
                usleep(1);
            }
            if (i == 432) {
                // In one thread, throw an error
                recurse_and_error(10); 
            }
            for (size_t j = 0; j < iter_count; j++) {
                // Do more busy work
                benchmark_control();
                usleep(1);
            }
        }
    }
    cerr << "Waiting for tasks" << endl;
    #pragma omp taskwait

    
    
    return 0;
}

// Register subcommand
static Subcommand vg_crash("crash", "throw an error", DEVELOPMENT, main_crash);

