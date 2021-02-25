// SPDX-FileCopyrightText: 2014 Erik Garrison
//
// SPDX-License-Identifier: MIT

/** \file test_main.cpp
 *
 * Defines the "vg test" subcommand, which runs unit tests.
 */


#include <omp.h>
#include <unistd.h>
#include <getopt.h>

#include <iostream>

#include "subcommand.hpp"

#include "../unittest/driver.hpp"

using namespace std;
using namespace vg;
using namespace vg::subcommand;

// No help_test is necessary because the unit testing library takes care of
// complaining about missing options.

int main_test(int argc, char** argv){
    // Forward arguments along to the main unit test driver
    return vg::unittest::run_unit_tests(argc, argv);
}

// Register subcommand
static Subcommand vg_test("test", "run unit tests", DEVELOPMENT, main_test);

