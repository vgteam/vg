/**
 * \file cluster_main.cpp: experimental snarl clustering test harness
 */

#include <omp.h>
#include <unistd.h>
#include <getopt.h>

#include "subcommand.hpp"

#include "../seed_clusterer.hpp"
#include "../stream/vpkg.hpp"

#include <iostream>

using namespace std;
using namespace vg;
using namespace vg::subcommand;

void help_cluster(char** argv) {
    cerr
    << "usage: " << argv[0] << " cluster [options] input.gam > output.gam" << endl
    << "Find and cluster mapping seeds." << endl
    << endl
    << "basic options:" << endl
    << "  -x, --xg-name FILE            use this xg index (required)" << endl
    << "  -g, --gcsa-name FILE          use this GCSA2/LCP index pair (required; both FILE and FILE.lcp)" << endl
    << "  -s, --snarls FILE             cluster using these snarls (required)" << endl
    << "  -d, --dist-name FILE          cluster using this distance index (required)" << endl
    << "computational parameters:" << endl
    << "  -t, --threads INT             number of compute threads to use" << endl;
}

int main_cluster(int argc, char** argv) {

    if (argc == 2) {
        help_cluster(argv);
        return 1;
    }

    // initialize parameters with their default options
    string xg_name;
    string gcsa_name;
    string snarls_name;
    string distance_name;
    
    int c;
    optind = 2; // force optind past command positional argument
    while (true) {
        static struct option long_options[] =
        {
            {"help", no_argument, 0, 'h'},
            {"xg-name", required_argument, 0, 'x'},
            {"gcsa-name", required_argument, 0, 'g'},
            {"snarls", required_argument, 0, 's'},
            {"dist-name", required_argument, 0, 'd'},
            {"threads", required_argument, 0, 't'},
            {0, 0, 0, 0}
        };

        int option_index = 0;
        c = getopt_long (argc, argv, "hx:g:s:d:t:",
                         long_options, &option_index);


        // Detect the end of the options.
        if (c == -1)
            break;

        switch (c)
        {
            case 'x':
                xg_name = optarg;
                if (xg_name.empty()) {
                    cerr << "error:[vg cluster] Must provide XG file with -x." << endl;
                    exit(1);
                }
                break;
                
            case 'g':
                gcsa_name = optarg;
                if (gcsa_name.empty()) {
                    cerr << "error:[vg cluster] Must provide GCSA file with -g." << endl;
                    exit(1);
                }
                break;
                
            case 's':
                snarls_name = optarg;
                if (snarls_name.empty()) {
                    cerr << "error:[vg cluster] Must provide snarl file with -s." << endl;
                    exit(1);
                }
                break;
                
            case 'd':
                distance_name = optarg;
                if (snarls_name.empty()) {
                    cerr << "error:[vg cluster] Must provide distance index file with -d." << endl;
                    exit(1);
                }
                break;
                
            case 't':
            {
                int num_threads = parse<int>(optarg);
                if (num_threads <= 0) {
                    cerr << "error:[vg cluster] Thread count (-t) set to " << num_threads << ", must set to a positive integer." << endl;
                    exit(1);
                }
                omp_set_num_threads(num_threads);
            }
                break;
                
            case 'h':
            case '?':
            default:
                help_cluster(argv);
                exit(1);
                break;
        }
    }
    
    
    if (xg_name.empty()) {
        cerr << "error:[vg cluster] Finding clusters requires an XG index, must provide XG file (-x)" << endl;
        exit(1);
    }
    
    if (gcsa_name.empty()) {
        cerr << "error:[vg cluster] Finding clusters requires a GCSA2 index, must provide GCSA2 file (-g)" << endl;
        exit(1);
    }
    
    if (snarls_name.empty()) {
        cerr << "error:[vg cluster] Finding clusters requires snarls, must provide snarls file (-s)" << endl;
        exit(1);
    }
    
    if (distance_name.empty()) {
        cerr << "error:[vg cluster] Finding clusters requires a distance index, must provide distance index file (-d)" << endl;
        exit(1);
    }
    
    // create in-memory objects
    unique_ptr<xg::XG> xg_index = stream::VPKG::load_one<xg::XG>(xg_name);
    unique_ptr<gcsa::GCSA> gcsa_index = stream::VPKG::load_one<gcsa::GCSA>(gcsa_name);
    unique_ptr<gcsa::LCPArray> lcp_array = stream::VPKG::load_one<gcsa::LCPArray>(gcsa_name + ".lcp");
    unique_ptr<SnarlManager> snarl_manager = stream::VPKG::load_one<SnarlManager>(snarls_name);
    
    // Now load the distance index.
    unique_ptr<DistanceIndex> distance_index;
    {
        ifstream distance_stream(distance_name);
        if (!distance_stream) {
            cerr << "error:[vg cluster] Could not open " << distance_name << endl;
            exit(1);
        }
        // TODO: let this go through VPKG as soon as we get it loadable without passing a graph
        distance_index = make_unique<DistanceIndex>(xg_index.get(), snarl_manager.get(), distance_stream);
    }
   
    return 0;
}

// Register subcommand
static Subcommand vg_cluster("cluster", "find and cluster mapping seeds", DEVELOPMENT, main_cluster);


