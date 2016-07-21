#include <iostream>
#include <fstream>
#include <ctime>
#include <cstdio>
#include <getopt.h>
#include <sys/stat.h>
#include "gcsa.h"
// From gcsa2
#include "files.h"
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
#include "pileup.hpp"
#include "caller.hpp"
#include "deconstructor.hpp"
#include "vectorizer.hpp"
#include "sampler.hpp"
#include "filter.hpp"
#include "google/protobuf/stubs/common.h"
#include "progress_bar.hpp"
#include "vg_git_version.hpp"
#include "IntervalTree.h"
#include "genotyper.hpp"
#include "bubbles.hpp"
#include "translator.hpp"
#include "homogenizer.hpp"

using namespace std;
using namespace vg;

void help_homogenize(char** argv){
    cerr << "Usage: " << argv[0] << " homogenize [options] <graph.vg>" << endl
        << "Options:" << endl
        << " -x, --xg   an xg index for generating mappings." << endl
        << " -g, --gcsa a gcsa2 index for generating mappings" << endl
        << " -k, --kmer a kmer size to use for in-memory GCSA/XG construction, if no indicies are provided." << endl
        << endl;
}


int main_homogenize(int argc, char** argv){
    string gcsa_name = "";
    string lcp_name = "";
    string xg_name = "";

    int kmer_size = 11;

    bool debug = false;
    bool in_mem_path_only = false;
    bool build_in_memory = false;


    if (argc <= 2) {
        help_homogenize(argv);
        return 1;
    }

    int c;
    optind = 2; // force optind past command positional argument
    while (true) {
        static struct option long_options[] =
        {
            {"xg", required_argument, 0, 'x'},
            {"gcsa", required_argument, 0, 'g'},
            {"kmer", required_argument, 0, 'k'},
            {0, 0, 0, 0}

        };
        int option_index = 0;
        c = getopt_long (argc, argv, "g:x:k:",
                long_options, &option_index);

        // Detect the end of the options.
        if (c == -1)
            break;

        switch (c)
        {
            case 'x':
                xg_name = optarg;
                break;
            case 'g':
                gcsa_name = optarg;
                lcp_name = gcsa_name + ".lcp";
                break;
            case 'k':
                kmer_size = atoi(optarg);
                break;
            case 'h':
            case '?':
            default:
                help_homogenize(argv);
                abort();
        }

    }

    VG* graph;
    string file_name = argv[optind];
    if (file_name == "-") {
        graph = new VG(std::cin);
    } else {
        ifstream in;
        in.open(file_name.c_str());
        graph = new VG(in);
    }

    gcsa::GCSA* gcsa_index = nullptr;
    gcsa::LCPArray* gcsa_lcp = nullptr;
    xg::XG* xg_index = nullptr;


    // TODO filter by whether a read is on the ref path
    // i.e.:
    // 1. Copy the graph
    VG* o_graph = new VG(*graph);
    cerr << "GRAPH COPIED" << endl;
    // 2. Cache the reference path(s)
    // (Not sure how to grab these, so for now just grab the first path in the graph)
    map<string, list<Mapping> > cached_paths;
    set<string> kept_paths;
    set<string> removed_paths;
    int i = 0;
    for (auto x : o_graph->paths._paths){
        cached_paths[x.first] = x.second;
        if (i == 0){
            kept_paths.insert(x.first);
        }
        else{
            removed_paths.insert(x.first);
        }
        i++;
    }
    //
    Paths p = graph->paths;
    //
    // 3. Remove all paths in the graph, except the reference
    (o_graph->paths).remove_paths(removed_paths);
    //
    //
    //graph = o_graph;
    delete graph;

    if (xg_name.empty() && gcsa_name.empty()){
        build_in_memory = true;
    }

    if (build_in_memory) {
        xg_index = new xg::XG(o_graph->graph);
        assert(kmer_size);
        int doubling_steps = 2;
        o_graph->build_gcsa_lcp(gcsa_index, gcsa_lcp, kmer_size, in_mem_path_only, false, 2);
        //delete graph;
    } else {
        // We try opening the file, and then see if it worked
        ifstream xg_stream(xg_name);

        if(xg_stream) {
            // We have an xg index!
            if(debug) {
                cerr << "Loading xg index " << xg_name << "..." << endl;
            }
            xg_index = new xg::XG(xg_stream);
        }

        ifstream gcsa_stream(gcsa_name);
        if(gcsa_stream) {
            // We have a GCSA index too!
            if(debug) {
                cerr << "Loading GCSA2 index " << gcsa_name << "..." << endl;
            }
            gcsa_index = new gcsa::GCSA();
            gcsa_index->load(gcsa_stream);
        }

        string lcp_name = gcsa_name + ".lcp";
        ifstream lcp_stream(lcp_name);
        if (lcp_stream) {
            if(debug) {
                cerr << "Loading LCP index " << gcsa_name << "..." << endl;
            }
            gcsa_lcp = new gcsa::LCPArray();
            gcsa_lcp->load(lcp_stream);
        }
    }

    Homogenizer hh;

    hh.homogenize(o_graph, xg_index, gcsa_index, gcsa_lcp);

    /* stream out graph */

    o_graph->serialize_to_ostream(std::cout);
    delete o_graph;

}
