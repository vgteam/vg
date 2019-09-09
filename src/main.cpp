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
#include "preflight.hpp"
#include "io/register_libvg_io.hpp"

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

// We make sure to compile main for the lowest common denominator architecture.
// This works on GCC and Clang. But we have to declare main and then define it.
// This *doesn't* work on Mac with GNU GCC and Apple libc++, so we exclude that combination.
#if (!defined(__GNUC__) || !defined(_LIBCPP_VERSION) || !defined(__APPLE__))
int main(int argc, char *argv[]) __attribute__((__target__("arch=x86-64")));
#endif

int main(int argc, char *argv[]) {

    // Make sure the system meets system requirements (i.e. has all the instructions we need)
    preflight_check();

    // Set up stack trace support from crash.hpp
    enable_crash_handling();
    
    // Tell the IO library about libvg types.
    // TODO: Make a more generic libvg startup function?
    if (!vg::io::register_libvg_io()) {
        cerr << "error[vg]: Could not register libvg types with libvgio" << endl;
        return 1;
    }
    
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
