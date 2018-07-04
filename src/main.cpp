// Needed for crash.hpp to work because it uses newer types
#define _POSIX_C_SOURCE 200809L

#include <iostream>
#include <fstream>
#include <ctime>
#include <cstdio>
#include <csignal>
#include <getopt.h>
#include <sys/stat.h>

#include "google/protobuf/stubs/common.h"
#include "version.hpp"
#include "utility.hpp"
#include "crash.hpp"

// New subcommand system provides all the subcommands that used to live here
#include "subcommand/subcommand.hpp"

using namespace std;
using namespace google::protobuf;
using namespace vg;


void vg_help(char** argv) {
    cerr << "vg: variation graph tool, version " << Version::get_short() << endl
         << endl
         << "usage: " << argv[0] << " <command> [options]" << endl
         << endl
         << vg::subcommand::PIPELINE << ":" << endl;
         
     vg::subcommand::Subcommand::for_each(vg::subcommand::PIPELINE, [](const vg::subcommand::Subcommand& command) {
        // Announce every subcommand we have
        
        // Pad all the names so the descriptions line up
        string name = command.get_name();
        name.resize(14, ' ');
        cerr << "  -- " << name << command.get_description() << endl;
     });
     
     cerr << endl << "For more commands, type `vg help`." << endl;
 }

int main(int argc, char *argv[])
{

    // Set up stack trace support from crash.hpp
    enable_crash_handling();

    // set a higher value for tcmalloc warnings
    setenv("TCMALLOC_LARGE_ALLOC_REPORT_THRESHOLD", "1000000000000000", 1);

    if (argc == 1) {
        vg_help(argv);
        return 1;
    }
    
    auto* subcommand = vg::subcommand::Subcommand::get(argc, argv);
    if (subcommand != nullptr) {
        // We found a matching subcommand, so run it
        return (*subcommand)(argc, argv);
    } else {
        // No subcommand found
        string command = argv[1];
        cerr << "error:[vg] command " << command << " not found" << endl;
        vg_help(argv);
        return 1;
    }

}
