// SPDX-FileCopyrightText: 2014 Erik Garrison
//
// SPDX-License-Identifier: MIT

#include "benchmark.hpp"
#include <vector>
#include <iostream>
#include <cassert>
#include <numeric>
#include <cmath>
#include <iomanip>

/**
 * \file benchmark.hpp: implementations of benchmarking functions
 */
 
namespace vg {
using namespace std;

double BenchmarkResult::score() const {
    // We comnpute a score in points by comparing the experimental and control runtimes.
    // Higher is better.
    
    // How many tests can we run per control run?
    double tests_per_control = (double)control_mean.count() / (double)test_mean.count();
    return tests_per_control * 1000;
}

double BenchmarkResult::score_error() const {
    // Do error propagation for the score calculation
    
    // Set up some abstract variables according to the notation at
    // <https://en.wikipedia.org/wiki/Propagation_of_uncertainty#Example_formulas>
    double A = control_mean.count();
    double stddev_A = control_stddev.count();
    double B = test_mean.count();
    double stddev_B = test_stddev.count();
    double f = score() / 1000;
    
    // Do the error propagation
    // Assume the variables are uncorrelated (covariance = 0)
    // Not because that's really true but because we don't track the covariance
    double err = abs(f) * sqrt(pow(stddev_A / A, 2) + pow(stddev_B / B, 2));
    
    // Scale up because of scaling factor on the score.
    return err * 1000;
}

ostream& operator<<(ostream& out, const BenchmarkResult& result) {
    // Dump it as a partial TSV line
    
    // We want to report times in fractional us
    using frac_secs = chrono::duration<double, std::micro>;
    
    // Save stream settings
    auto initial_precision = out.precision();
    auto initial_flags = out.flags();
    
    // Set up formatting
    cout << setprecision(2) << scientific;
    
    out << result.runs;
    out << "\t";
    out << chrono::duration_cast<frac_secs>(result.test_mean).count();
    out << "\t";
    out << chrono::duration_cast<frac_secs>(result.test_stddev).count();
    out << "\t";
    out << chrono::duration_cast<frac_secs>(result.control_mean).count();
    out << "\t";
    out << chrono::duration_cast<frac_secs>(result.control_stddev).count();
    out << "\t";
    
    // Scores get different formatting
    cout << fixed;
    
    out << result.score();
    out << "\t";
    out << result.score_error();
    out << "\t";
    out << result.name;
    
    out.precision(initial_precision);
    out.flags(initial_flags);
    
    return out;
}

void benchmark_control() {
    // We need to do something that takes time.
    
    vector<size_t> things;
    
    for (size_t i = 0; i < 100; i++) {
        things.push_back(i ^ (i/5));
    }
    
    size_t max_thing = 0;
    
    for (auto& thing : things) {
        for (auto& other_thing : things) {
            max_thing = max(max_thing, max(thing, other_thing));
            other_thing = other_thing ^ (thing << 5);
        }
    }
    
    size_t total = 0;
    for (auto& thing : things) {
        total += thing;
    }
    
    // These are the results of that arbitrary math
    assert(max_thing == 18444166782598024656U);
    assert(total == 17868247911721767448U);
    
}

BenchmarkResult run_benchmark(const string& name, size_t iterations, const function<void(void)>& under_test) {
    return run_benchmark(name, iterations, []() {}, under_test);
}

BenchmarkResult run_benchmark(const string& name, size_t iterations, const function<void(void)>& setup,
    const function<void(void)>& under_test) {

    // We'll fill this in with the results of the benchmark run
    BenchmarkResult to_return;
    to_return.runs = iterations;
    to_return.name = name;
    
    // Where do we put our test runtime samples?
    // They need to be normal integral types so we can feasibly square them.
    vector<benchtime::rep> test_samples;
    // And the samples for the control?
    vector<benchtime::rep> control_samples;
    
    // We know how big they will be
    test_samples.reserve(iterations);
    control_samples.reserve(iterations);
    
    for (size_t i = 0; i < iterations; i++) {
        // For each iteration
        
        // Run the setup function
        setup();
        
        // Run the function under test
        auto test_start = chrono::high_resolution_clock::now();
        under_test();
        auto test_stop = chrono::high_resolution_clock::now();
        
        // And run the control
        auto control_start = chrono::high_resolution_clock::now();
        benchmark_control();
        auto control_stop = chrono::high_resolution_clock::now();
        
        // Make sure time went forward and nobody changed the clock (as far as we know)
        assert(test_stop > test_start);
        assert(control_stop > control_start);
        
        // Compute the runtimes and save them
        test_samples.push_back(chrono::duration_cast<benchtime>(test_stop - test_start).count());
        control_samples.push_back(chrono::duration_cast<benchtime>(control_stop - control_start).count());
    }
    
    // Calculate the moments with magic numeric algorithms
    // See <https://stackoverflow.com/a/7616783>
    
    // Total up the ticks for the test
    benchtime::rep test_total = accumulate(test_samples.begin(), test_samples.end(),
        benchtime::zero().count());
    // And make a duration for the mean
    to_return.test_mean = benchtime(test_total / iterations);
    
    // Then total up the squares of the tick counts
    benchtime::rep test_square_total = inner_product(test_samples.begin(), test_samples.end(),
        test_samples.begin(), benchtime::zero().count());
    // Calculate the standard deviation in ticks, and represent it as a duration
    to_return.test_stddev = benchtime((benchtime::rep) sqrt(test_square_total / iterations -
        to_return.test_mean.count() * to_return.test_mean.count()));
    
    // Similarly for the control
    benchtime::rep control_total = accumulate(control_samples.begin(), control_samples.end(),
        benchtime::zero().count());
    to_return.control_mean = benchtime(control_total / iterations);
    
    benchtime::rep control_square_total = inner_product(control_samples.begin(), control_samples.end(), 
        control_samples.begin(), benchtime::zero().count());
    to_return.control_stddev = benchtime((benchtime::rep) sqrt(control_square_total / iterations -
        to_return.control_mean.count() * to_return.control_mean.count()));
    
    return to_return;
    
}

}













