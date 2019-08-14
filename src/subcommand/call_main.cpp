/** \file call_main.cpp
 *
 * Defines the "vg call" subcommand, which calls variation from an augmented graph
 */

#include <omp.h>
#include <unistd.h>
#include <getopt.h>

#include <list>
#include <fstream>

#include "subcommand.hpp"

#include "../option.hpp"

#include "../vg.hpp"
#include "../xg.hpp"
#include "../support_caller.hpp"
#include <vg/io/stream.hpp>
#include <vg/io/vpkg.hpp>



using namespace std;
using namespace vg;
using namespace vg::subcommand;

void help_call(char** argv, ConfigurableParser& parser) {
    cerr << "usage: " << argv[0] << " call [options] <augmented-graph.vg> > output.vcf" << endl
         << "Output variant calls in VCF or Loci format given a graph and pileup" << endl
         << endl
         << "genotyper options:" << endl
        
         << "general options:" << endl
         << "    -z, --translation FILE      input translation table" << endl
         << "    -b, --base-graph FILE       base graph.  currently needed for XREF tag" << endl
         << "    -g, --snarls FILE           snarls of input graph (to save re-computing)" << endl
         << "    -h, --help                  print this help message" << endl
         << "    -p, --progress              show progress" << endl
         << "    -v, --verbose               print information and warnings about vcf generation" << endl
         << "    -t, --threads N             number of threads to use" << endl;
     
     // Then report more options
     parser.print_help(cerr);
}

int main_call(int argc, char** argv) {

    string translation_file_name;

    // todo: get rid of this.  only used to check if deletion edge is novel.  must be some
    // way to get that out of the translations
    string base_graph_file_name;

    string snarl_file_name;
    
    // This manages conversion from an augmented graph to a VCF, and makes the
    // actual calls.
    SupportCaller support_caller;

    bool show_progress = false;
    int thread_count = 0;

    static const struct option long_options[] = {
        {"base-graph", required_argument, 0, 'b'},
        {"translation", required_argument, 0, 'z'},
        {"snarls", required_argument, 0, 'g'},
        {"progress", no_argument, 0, 'p'},
        {"verbose", no_argument, 0, 'v'},
        {"threads", required_argument, 0, 't'},
        {"help", no_argument, 0, 'h'},
        {0, 0, 0, 0}
    };
    static const char* short_options = "z:b:pvt:hg:";
    optind = 2; // force optind past command positional arguments

    // This is our command-line parser
    ConfigurableParser parser(short_options, long_options, [&](int c) {
        // Parse all the options we have defined here.
        switch (c)
        {
        case 'z':
            translation_file_name = optarg;
            break;
        case 'b':
            base_graph_file_name = optarg;
            break;
        case 'g':
            snarl_file_name = optarg;
            break;
        case 'p':
            show_progress = true;
            break;
        case 'v':
            support_caller.verbose = true;
            break;
        case 't':
            thread_count = parse<int>(optarg);
            break;
        case 'h':
        case '?':
            /* getopt_long already printed an error message. */
            help_call(argv, parser);
            exit(1);
            break;
        default:
          abort ();
        }
    });
    // Register the support_caller converter for configuring with its options.
    parser.register_configurable(&support_caller);

    if (argc <= 3) {
        help_call(argv, parser);
        return 1;
    }
    
    // Parse the command line options, updating optind.
    parser.parse(argc, argv);

    if (thread_count != 0) {
        // Use a non-default number of threads
        omp_set_num_threads(thread_count);
    }
    thread_count = get_thread_count();

    // Parse the arguments
    if (optind >= argc) {
        help_call(argv, parser);
        return 1;
    }
    string graph_file_name = get_input_file_name(optind, argc, argv);

    if (string(support_caller.support_file_name).empty() ==
        string(support_caller.pack_file_name).empty()) {
        cerr << "[vg call]: Support file must be specified with either -s (or -P)" << endl;
        return 1;
    }
    

    if (string(support_caller.pack_file_name).empty() && translation_file_name.empty()) {
        cerr << "[vg call]: Translation file must be specified with -Z" << endl;
        return 1;
    }

    // read the vcf and accompanying fasta files if specified
    if (!((string)(support_caller.recall_vcf_filename)).empty()) {
        if (support_caller.variant_offset != 0) {
            cerr << "error: [vg call] vcf offset (-o) not supported in recall mode (-f)" << endl;
            return 1;
        }
        support_caller.variant_file.parseSamples = false;
        support_caller.variant_file.open(support_caller.recall_vcf_filename);
        if (!support_caller.variant_file.is_open()) {
            cerr << "error: [vg call] could not open " << (string)support_caller.recall_vcf_filename << endl;
            return 1;
        }

        // load up the fasta
        if (!((string)support_caller.recall_ref_fasta_filename).empty()) {
            support_caller.ref_fasta = unique_ptr<FastaReference>(new FastaReference);
            support_caller.ref_fasta->open(support_caller.recall_ref_fasta_filename);
        }
        if (!((string)support_caller.recall_ins_fasta_filename).empty()) {
            support_caller.ins_fasta = unique_ptr<FastaReference>(new FastaReference);
            support_caller.ins_fasta->open(support_caller.recall_ins_fasta_filename);
        }
    }
    
    // read the graph
    if (show_progress) {
        cerr << "Reading input graph" << endl;
    }
    VG* graph;
    get_input_file(graph_file_name, [&](istream& in) {
        graph = new VG(in);
    });

    // and the base graph
    VG* base_graph = NULL;
    if (!base_graph_file_name.empty()) {
        get_input_file(base_graph_file_name, [&](istream& in) {
                base_graph = new VG(in);
            });
    }

    SupportAugmentedGraph augmented_graph;
    // Move our input graph into the augmented graph
    swap(augmented_graph.graph, *graph); 

    augmented_graph.base_graph = base_graph;
    
    delete graph;

    // Load the supports
    if (!string(support_caller.support_file_name).empty()) {
        ifstream support_file(support_caller.support_file_name);
        if (!support_file) {
            cerr << "[vg call]: Unable to load supports file: "
                 << string(support_caller.support_file_name) << endl;
            return 1;
        }
        augmented_graph.load_supports(support_file);
    } else {
        assert(!string(support_caller.pack_file_name).empty());
        if (string(support_caller.xg_file_name).empty()) {
            cerr << "[vg call]: pack support (-P) requires xg index (-x)" << endl;
            return 1;
        }
        unique_ptr<PathPositionHandleGraph> xgidx = vg::io::VPKG::load_one<PathPositionHandleGraph>(support_caller.xg_file_name);
        augmented_graph.load_pack_as_supports(support_caller.pack_file_name, xgidx.get());
        // make sure we're ignoring quality, as it's not read from the pack
        bool& usc = support_caller.use_support_count;
        usc = true;
    }        

    // Load the translations
    if (!translation_file_name.empty()) {
        ifstream translation_file(translation_file_name.c_str());
        if (!translation_file) {
            cerr << "[vg call]: Unable to load translations file: " << translation_file_name << endl;
            return 1;
        }
        augmented_graph.load_translations(translation_file);
    }

    // Load or compute the snarls
    unique_ptr<SnarlManager> snarl_manager;    
    if (!snarl_file_name.empty()) {
        ifstream snarl_file(snarl_file_name.c_str());
        if (!snarl_file) {
            cerr << "[vg call]: Unable to load snarls file: " << snarl_file_name << endl;
            return 1;
        }
        snarl_manager = vg::io::VPKG::load_one<SnarlManager>(snarl_file);
    } else {
        CactusSnarlFinder finder(augmented_graph.graph);
        snarl_manager = unique_ptr<SnarlManager>(new SnarlManager(std::move(finder.find_snarls())));
    }
    
    
    if (show_progress) {
        cerr << "Calling variants with support caller" << endl;
    }

    // project the augmented graph to a reference path
    // in order to create a VCF of calls.
    support_caller.call(augmented_graph, *snarl_manager.get(), {});

    delete base_graph;
    return 0;
}

// Register subcommand
static Subcommand vg_call("call", "call variants on an augmented graph", PIPELINE, 5, main_call);

