#include <iostream>
#include <vector>
#include <getopt.h>
#include <functional>
#include <regex>
#include "subcommand.hpp"
#include "srpe.hpp"
#include "stream.hpp"
#include "index.hpp"
#include "position.hpp"
#include "vg.pb.h"
#include "path.hpp"
#include "genotypekit.hpp"
#include "genotyper.hpp"
#include "path_index.hpp"
#include "vg.hpp"
#include <math.h>
#include "srpe.hpp"
#include "filter.hpp"
#include "utility.hpp"
#include "Variant.h"
#include "translator.hpp"
#include "Fasta.h"
#include "IntervalTree.h"

using namespace std;
using namespace vg;
using namespace vg::subcommand;


void help_srpe(char** argv){
    cerr << "Usage: " << argv[0] << " srpe [options] <data.gam> <graph.vg>" << endl
    << "Options: " << endl 
    << "  -p / --ref-path" << endl
    << "  -x / --xg" << endl 
    << "  -g / --gcsa" << endl 
    << endl;
}



int main_srpe(int argc, char** argv){
    string alignment_file = "";
    string gam_index_name = "";
    string graph_name = "";
    string xg_name = "";
    string gcsa_name = "";
    string lcp_name = "";

    string spec_vcf = "";
    string ref_fasta = "";
    string ins_fasta = "";

    string augmented_graph_name = "";
    bool augment_paths = true;

    string ref_path = "";

    int max_iter = 2;
    int max_frag_len = 10000;
    int min_soft_clip = 20;
    bool remap = false;

    bool do_all = false;

    vector<string> search_types;
    search_types.push_back("DEL");

    int threads = 1;

    if (argc <= 2) {
        help_srpe(argv);
        return 1;
    }

    int c;
    optind = 2; // force optind past command positional argument
    while (true) {
        static struct option long_options[] =
        {
            {"max-iter", required_argument, 0, 'm'},
            {"xg-index", required_argument, 0, 'x'},
            {"augmented", required_argument, 0, 'a'},
            {"help", no_argument, 0, 'h'},
            {"gcsa-index", required_argument, 0, 'g'},
            {"specific", required_argument, 0, 'S'},
            {"recall", no_argument, 0, 'R'},
            {"insertions", required_argument, 0, 'I'},
            {"reference", required_argument, 0, 'r'},
            {"threads", required_argument, 0, 't'},
            {"ref-path", required_argument, 0, 'p'},
            {"remap", no_argument, 0, 'z'},
            {0, 0, 0, 0}

        };
        int option_index = 0;
        c = getopt_long (argc, argv, "hzx:g:m:S:RI:r:t:a:wp:",
                long_options, &option_index);

        // Detect the end of the options.
        if (c == -1)
            break;

        switch (c)
        {
            case 'a':
                augmented_graph_name = optarg;
                break;
            case 'z':
                remap = true;
                break;
            case 'm':
                max_iter = atoi(optarg);
                break;

            case 't':
                threads = atoi(optarg);
                break;

            case 'R':
                do_all = true;
                break;
            case 'x':
                xg_name = optarg;
                break;
            case 'g':
                gcsa_name = optarg;
                break;
            case 'S':
                spec_vcf = optarg;
                break;
            case 'r':
                ref_fasta = optarg;
                break;
            case 'I':
                ins_fasta = optarg;
                break;
            case 'p':
                ref_path = optarg;
                break;
            case 'h':
            case '?':
            default:
                help_srpe(argv);
                abort();
        }

    }

    omp_set_num_threads(threads);


    SRPE srpe;

    


    alignment_file = argv[optind];
    //gam_index_name = argv[++optind];
    graph_name = argv[++optind];

    xg::XG* xg_ind = new xg::XG();
    Index gamind;

    vg::VG* graph;

    if (!xg_name.empty()){
        ifstream in(xg_name);
        xg_ind->load(in);
        srpe.ff.set_my_xg_idx(xg_ind);
    }
    // Set GCSA indexes
    if (!gcsa_name.empty()){
            ifstream gcsa_stream(gcsa_name);
            srpe.ff.gcsa_ind = new gcsa::GCSA();
            srpe.ff.gcsa_ind->load(gcsa_stream);
            string lcp_name = gcsa_name + ".lcp";
            ifstream lcp_stream(lcp_name);
            srpe.ff.lcp_ind = new gcsa::LCPArray();
            srpe.ff.lcp_ind->load(lcp_stream);
    }
    if (!xg_name.empty()){
        ifstream xgstream(xg_name);
        xg_ind->load(xgstream);
        srpe.ff.set_my_xg_idx(xg_ind);
    }
    srpe.ff.init_mapper();
    // else{

    // }

    

    return 0;
}

static Subcommand vg_srpe ("srpe", "graph-external SV detection", main_srpe);

