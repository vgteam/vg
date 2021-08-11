// mod.cpp: define the "normalize" subcommand, which realigns snarls to produce more
// efficient representations of snarls.

#include <getopt.h>
#include <omp.h>
#include <unistd.h>

#include <gbwtgraph/gbwtgraph.h>

#include "subcommand.hpp"

//todo: should be able to remove '../../include/...' and replace with e.g. <bdsg/hash...> 
#include "../../include/bdsg/hash_graph.hpp"
#include "../../include/vg/io/vpkg.hpp"
// #include "../algorithms/0_draft_haplotype_realignment.hpp"
#include "../algorithms/0_draft_snarl_normalization_evaluation.cpp"
#include "../algorithms/0_oo_normalize_snarls.hpp"
#include "../gbwt_helper.hpp"

#include "../io/save_handle_graph.hpp"

#include "../algorithms/0_snarl_analyzer.hpp"

#include <chrono> // for high_resolution_clock

using namespace std;
using namespace vg;
using namespace vg::subcommand;

void help_normalize(char **argv) {
    cerr
        << "usage: " << argv[0] << " normalize [options] <graph.vg> >[mod.vg]" << endl
        << "Modifies snarls, outputs modified on stdout." << endl
        << endl
        << "options:" << endl
        << "    -g, --gbwt       gbwt corresponding to hashgraph." << endl
        << "    -s, --snarls       snarls file corresponding to hashgraph." << endl
        << "    -m, --max_alignment_size       limits the number of threads that will "
           "be aligned in any snarl. If exceeded, program skips snarl. Default is none "
           "threads. If you don't want to skip any snarls based on thread count, enter 0."
        << endl;
}

int main_normalize(int argc, char **argv) {

    if (argc == 2) {
        help_normalize(argv);
        return 1;
    }

    int max_alignment_size = INT_MAX; // default cutoff used to be 200 threads in a snarl.
    string gbwt;
    string snarls;
    string normalize_type = "all";
    int source = NULL; //todo: do something other than NULL to avoid the compiler warnings. 
    int sink = NULL;
    bool paths_right_to_left = false;
    bool evaluate = false;
    string snarl_sizes;
    bool snarl_sizes_skip_source_sink = false;
    int max_handle_size = INT_MAX;
    bool handles_in_snarl = false;


    int c;
    optind = 2; // force optind past command positional argument
    while (true) {
        static struct option long_options[] =

            {{"help", no_argument, 0, 'h'},
             {"gbwt", required_argument, 0, 'g'},
             {"snarls", required_argument, 0, 's'},
             {"normalize_type", required_argument, 0, 'n'},
             {"source", required_argument, 0, 'a'},
             {"sink", required_argument, 0, 'b'},
             {"paths_right_to_left", no_argument, 0, 'p'},
             {"max_alignment_size", required_argument, 0, 'm'},
             {"evaluate", no_argument, 0, 'e'},
             {"snarl_sizes", required_argument, 0, 'i'},
             {"snarl_sizes_skip_source_sink", no_argument, 0, 'k'},
             {"max_handle_size", required_argument, 0, 'h'}, // currently, default is INT_MAX. This is for compatibility with changes in graph size measures (0_snarl_analyzer). Eventually should change to handle size standard. 
             {"handles_in_snarl", no_argument, 0, 'x'}, // used in conjunction with arguments source and sink. Will print all the node ids in between source and sink, inclusive.
             {0, 0, 0, 0}};

        int option_index = 0;
        c = getopt_long(argc, argv, "hg:s:n:a:b:pm:ei:kh:x", long_options, &option_index);

        // Detect the end of the options.
        if (c == -1)
            break;

        switch (c) {

        case 'g':
            gbwt = optarg;
            break;

        case 's':
            snarls = optarg;
            break;
            
        case 'n':
            // can be "all" (default), "one", or "none" (latter for if you just want to run evaluate) //todo: make better!
            normalize_type = optarg;
            break;
            
        case 'a':
            source = parse<int>(optarg);
            break;

        case 'b':
            sink = parse<int>(optarg);
            break;
            
        case 'p':
            paths_right_to_left = true;
            break;
            
        case 'm':
            max_alignment_size = parse<int>(optarg);
            // if max_alignment_size is 0, then that signifies that it should actually be
            // infinite, i.e. that we should not exclude any snarls.
            if (max_alignment_size == 0) {
                max_alignment_size = INT_MAX;
            }
            break;

        case 'e':
            evaluate = true;
            break;
        
        case 'i':
            snarl_sizes = optarg;
            normalize_type = "none";
            break;

        case 'k':
            snarl_sizes_skip_source_sink = true;
            break;
            
        case 'h':
            max_handle_size = parse<int>(optarg);
            break;
            
        case 'x':
            handles_in_snarl = true;
            normalize_type = "none";
            break;
            
        default:
            cerr << "error:[vg normalize] abort" << endl;
            abort();
        }
    }

    //getting graph of any type, except non-mutable graphs (e.g., xg)
    unique_ptr<MutablePathDeletableHandleGraph> graph;
    get_input_file(optind, argc, argv, [&](istream &in) {
        graph = vg::io::VPKG::load_one<MutablePathDeletableHandleGraph>(in);
    });

    if (normalize_type!="all" && normalize_type!="one" && normalize_type != "none") 
    {
        cerr << "please enter a valid normalize_type: all, one, or none." << endl;
    }

    if (normalize_type == "all" || normalize_type == "one") {
        cerr << "running normalize!" << endl;

        /// Build the gbwt:
        ifstream gbwt_stream;
        gbwt_stream.open(gbwt);

        // Load the GBWT from its container
        unique_ptr<gbwt::GBWT> gbwt;
        gbwt = vg::io::VPKG::load_one<gbwt::GBWT>(gbwt_stream);
        gbwtgraph::GBWTGraph haploGraph = gbwtgraph::GBWTGraph(*gbwt, *graph);

        std::ifstream snarl_stream;
        string snarl_file = snarls;
        snarl_stream.open(snarl_file);

        if (!snarl_stream) {
            cerr << "error:[vg normalize] Cannot open Snarls file " << snarl_file << endl;
            exit(1);
        }
        // Record start time
        auto start = chrono::high_resolution_clock::now();

        algorithms::SnarlNormalizer normalizer =
            algorithms::SnarlNormalizer(*graph, haploGraph, max_alignment_size, max_handle_size);

        if (normalize_type == "all")
        {
            normalizer.normalize_top_level_snarls(snarl_stream);
        }
        else if (normalize_type == "one")
        {
            if (source == NULL && sink == NULL)
            {
                cerr << "ERROR: please provide a source and sink for the snarl you want to normalize." << endl;
                return 0;
            }
            vector<int> error_record = normalizer.normalize_snarl(source, sink, paths_right_to_left);
            if (!(error_record[0] || error_record[1] ||
                    error_record[2] || error_record[3] ||
                    error_record[6]))
            {
                cerr << "snarl starting at " << source << " and ending at " << sink << " normalized." << endl;
                cerr << "amount of sequence in normalized snarl before normalization: "
                    << error_record[4] << endl;
                cerr << "amount of sequence in normalized snarl after normalization: "
                    << error_record[5] << endl;            }
            else
            {
                //todo: make it so it only prints the relevant message:
                cerr << "snarl skipped because...\nthey exceeded the size limit ("
                    << error_record[0] << " snarls),\n"
                    << "had haplotypes starting/ending in the middle of the snarl ("
                    << error_record[1] << "),\n"
                    << "the snarl was cyclic (" 
                    << error_record[3] << " snarls),\n"
                    << " there were handles not connected by the gbwt info ("
                    << error_record[2] << " snarls),\n" 
                    << "the snarl was cyclic (" << error_record[3] << " snarls),\n"
                    << "or the snarl was trivial - composed of only one or two nodes ("
                    << error_record[6] << " snarls)."
                    << endl;
            }

        }
        // // run test code on all snarls in graph. (non obj-oriented code)
        // disambiguate_top_level_snarls(*graph, haploGraph, snarl_stream,
        // max_alignment_size);

        // Record end time
        auto finish = std::chrono::high_resolution_clock::now();
        chrono::duration<double> elapsed = finish - start;
        cerr << "Elapsed time: " << elapsed.count() << " s\n";
    }

    if (evaluate) {
        std::ifstream snarl_stream;
        string snarl_file = snarls;
        snarl_stream.open(snarl_file);
        cerr << "about to evaluate normalized snarls" << endl;
        vg::evaluate_normalized_snarls(snarl_stream);
    }

    // snarl_sizes identifies the size of every snarl, outputs in a document specified with format "source\tsink\tsize\n"
    if (snarl_sizes.size() != 0) 
    { //todo: add "evaluate" functions to snarl_analyzer; rename snarl_sizes doc as snarl_analyzer.hpp,etc
        std::ifstream snarl_stream;
        snarl_stream.open(snarls);
        if (!snarl_stream) {
            cerr << "error:[vg normalize] Cannot open Snarls file " << snarls << endl;
            exit(1);
        }

        algorithms::SnarlAnalyzer sizes = algorithms::SnarlAnalyzer(*graph, snarl_stream, snarl_sizes_skip_source_sink);

        sizes.output_snarl_sizes(snarl_sizes);

    }
    
    if (handles_in_snarl)
    {
        if (source == NULL && sink == NULL)
        {
            cerr << "error:[vg normalize] please enter a values for source and sink to define the snarl." << endl;
        }
        else
        {
            algorithms::print_handles_in_snarl(*graph, source, sink);
        }
    }

    if (normalize_type!="none") {
        // Save the modified graph
        // vg::io::save_handle_graph(graph.get(), std::cout);
        
        //todo: maybe rewrite to mimic mod_main.
        // vg::io::VPKG::save(*dynamic_cast<bdsg::HashGraph *>(graph.get()), cout);

        // graph->serialize(std::cout);
    }


    // delete graph;

    return 0;
}

// Register subcommand
static Subcommand vg_normalize("normalize",
                               "edit snarls to reduce information duplication", TOOLKIT,
                               main_normalize);