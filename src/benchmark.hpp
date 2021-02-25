// SPDX-FileCopyrightText: 2014 Erik Garrison
//
// SPDX-License-Identifier: MIT

#ifndef VG_BENCHMARK_HPP_INCLUDED
#define VG_BENCHMARK_HPP_INCLUDED

#include <chrono>
#include <functional>
#include <iostream>
#include <string>

/** 
 * \file benchmark.hpp
 * Contains utilities for running (micro)benchmarks of vg code.
 */
 
namespace vg {

using namespace std;

/// We define a duration type for expressing benchmark times in.
using benchtime = chrono::nanoseconds;

/**
 * Represents the results of a benchmark run. Tracks the mean and standard
 * deviation of a number of runs of a function under test, interleaved with runs
 * of a standard control function.
 */ 
struct BenchmarkResult {
    /// How many runs were done
    size_t runs;
    /// What was the mean runtime of each test run
    benchtime test_mean;
    /// What was the standard deviation of test run times
    benchtime test_stddev;
    /// What was the mean runtime of each control run
    benchtime control_mean;
    /// What was the standard deviation of control run times
    benchtime control_stddev;
    /// What was the name of the test being run
    string name;
    /// How many control-standardized "points" do we score?
    double score() const;
    /// What is the uncertainty on the score?
    double score_error() const;
};

/**
 * Benchmark results can be output to streams
 */
ostream& operator<<(ostream& out, const BenchmarkResult& result);

/**
 * The benchmark control function, designed to take some amount of time that might vary with CPU load.
 */
void benchmark_control();

/**
 * Run the given function the given number of times, interleaved with runs of
 * the control function, and return a BenchmarkResult describing its
 * performance.
 */
BenchmarkResult run_benchmark(const string& name, size_t iterations, const function<void(void)>& under_test);

/**
 * Run a benchmark with a setup function.
 */
BenchmarkResult run_benchmark(const string& name, size_t iterations, const function<void(void)>&  setup, const function<void(void)>& under_test);


}

#endif
