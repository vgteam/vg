#include <iostream>
#include <fstream>
#include <ctime>
#include <cstdio>
#include <getopt.h>
#include <sys/stat.h>
#include "gcsa/gcsa.h"
#include "gcsa/algorithms.h"
#include "json2pb.h"
#include "vg.hpp"
#include "vg.pb.h"
#include "vg_set.hpp"
#include "index.hpp"
#include "mapper.hpp"
#include "Variant.h"
#include "Fasta.h"
#include "stream.hpp"
#include "alignment.hpp"
#include "convert.hpp"
#include "caller.hpp"
#include "deconstructor.hpp"
#include "filter.hpp"
#include "google/protobuf/stubs/common.h"
#include "progress_bar.hpp"
#include "version.hpp"
#include "genotyper.hpp"
#include "bubbles.hpp"
#include "readfilter.hpp"
#include "distributions.hpp"
// New subcommand system provides all the subcommands that used to live here
#include "subcommand/subcommand.hpp"

using namespace std;
using namespace google::protobuf;
using namespace vg;

    void help_sv(char** argv){
        cerr << "usage: " << argv[0] << " sv [options] <aln.gam>" << endl
            << "options: " << endl
            << " -g --graph <graph>.vg " << endl
            << " -m --mask <vcf>.vcf" << endl
            << endl;
    }


void vg_help(char** argv) {
    cerr << "vg: variation graph tool, version " << VG_VERSION_STRING << endl
         << endl
         << "usage: " << argv[0] << " <command> [options]" << endl
         << endl
         << "commands:" << endl;
         
     vg::subcommand::Subcommand::for_each([](const vg::subcommand::Subcommand& command) {
        // Announce every subcommand we have
        
        // Pad all the names so the descriptions line up
        string name = command.get_name();
        name.resize(14, ' ');
        cerr << "  -- " << name << command.get_description() << endl;
     });
         
     // Also announce all the old-style hardcoded commands
         cerr << "  -- genotype      compute genotypes from aligned reads" << endl;
 }

int main(int argc, char *argv[])
{

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
    }
    
    // Otherwise, fall abck on the old chain of if statements.

    //omp_set_dynamic(1); // use dynamic scheduling

    string command = argv[1];
    cerr << "error:[vg] command " << command << " not found" << endl;
    vg_help(argv);
    return 1;

}
