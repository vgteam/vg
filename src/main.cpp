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
#include "unittest/driver.hpp"
// New subcommand system provides all the subcommands that used to live here
#include "subcommand/subcommand.hpp"
#include "flow_sort.hpp"


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


void help_sort(char** argv){
    cerr << "usage: " << argv[0] << " sort [options] -i <input_file> -r <reference_name> > sorted.vg " << endl
         << "options: " << endl
         << "           -g, --gfa              input in GFA format" << endl
         << "           -i, --in               input file" << endl
         << "           -r, --ref              reference name" << endl
         << "           -w, --without-grooming no grooming mode" << endl
         << "           -f, --fast             sort using Eades algorithm, otherwise max-flow sorting is used" << endl   
         << endl;
}

int main_sort(int argc, char *argv[]) {

    //default input format is vg
    bool gfa_input = false;
    string file_name = "";
    string reference_name = "";
    bool without_grooming = false;
    bool use_fast_algorithm = false;
    int c;
    while (true) {
        static struct option long_options[] =
            {
                {"gfa", no_argument, 0, 'g'},
                {"in", required_argument, 0, 'i'},
                {"ref", required_argument, 0, 'r'},
                {"without-grooming", no_argument, 0, 'w'},
                {"fast", no_argument, 0, 'f'},
                {0, 0, 0, 0}
            };

        int option_index = 0;
        c = getopt_long (argc, argv, "i:r:gwf",
                         long_options, &option_index);

        /* Detect the end of the options. */
        if (c == -1)
            break;

        switch (c)
        {
        case 'g':
            gfa_input = true;
            break;
        case 'r':
            reference_name = optarg;
            break;
        case 'i':
            file_name = optarg;
            break;
        case 'w':
            without_grooming = true;
            break;
        case 'f':
            use_fast_algorithm = true;
            break;
        case 'h':
        case '?':
            /* getopt_long already printed an error message. */
            help_sort(argv);
            exit(1);
            break;
        default:
            abort ();
        }
    }
  
    if (reference_name.empty() || file_name.empty()) {
        help_sort(argv);
        exit(1);
    }
    
    ifstream in;
    std::unique_ptr<VG> graph;
    {
        in.open(file_name.c_str());        
        if (gfa_input) {
            graph.reset(new VG());
            graph->from_gfa(in);
        } else {
            graph.reset(new VG(in));
        }
    }
    FlowSort flow_sort(*graph.get());
    if (use_fast_algorithm) {
        flow_sort.fast_linear_sort(reference_name, !without_grooming);
    } else {
        flow_sort.max_flow_sort(reference_name);
    }
    
    graph->serialize_to_ostream(std::cout);
    in.close();
    return 0;
}

// No help_test is necessary because the unit testing library takes care of
// complaining about missing options.

int main_test(int argc, char** argv){
    // Forward arguments along to the main unit test driver
    return vg::unittest::run_unit_tests(argc, argv);
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
         cerr << "  -- genotype      compute genotypes from aligned reads" << endl
         << "  -- sort          sort variant graph using max flow algorithm or Eades fast heuristic algorithm" << endl
         << "  -- test          run unit tests" << endl;
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
    if (command == "test") {
        return main_test(argc, argv);
    } else if (command == "sort") {
        return main_sort(argc, argv);
    } else {
        cerr << "error:[vg] command " << command << " not found" << endl;
        vg_help(argv);
        return 1;
    }

    return 0;

}
