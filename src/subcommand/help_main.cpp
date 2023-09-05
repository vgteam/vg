/** \file help_main.cpp
 *
 * Defines the "vg help" subcommand, which describes subcommands.
 */


#include <omp.h>
#include <unistd.h>
#include <getopt.h>
#include <stdlib.h>

#include <iostream>

#include "subcommand.hpp"
#include "../version.hpp"

using namespace std;
using namespace vg;
using namespace vg::subcommand;

void help_help(char** argv){
    cerr << "usage: " << argv[0] << " help" << endl
         << endl;
}

int main_help(int argc, char** argv){

    cerr << "vg: variation graph tool, version " << Version::get_short() << endl
         << endl
         << "usage: " << argv[0] << " <command> [options]" << endl
         << endl;
         
     for (auto category : {PIPELINE, TOOLKIT, WIDGET, DEVELOPMENT}) {

         cerr << category << ":" << endl;
         
         vg::subcommand::Subcommand::for_each(category, [](const vg::subcommand::Subcommand& command) {
            // Announce every subcommand we have
            
            // Pad all the names so the descriptions line up
            string name = command.get_name();
            name.resize(14, ' ');
            cerr << "  -- " << name << command.get_description() << endl;
         });
         
         cerr << endl;
     }
     
     cerr << "For technical support, please visit: https://www.biostars.org/tag/vg/" << endl << endl;
     
    
    return 0;
}

// Register subcommand
static Subcommand vg_help("help", "show all subcommands", PIPELINE, main_help);

