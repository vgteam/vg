/** \file crash_main.cpp
 *
 * Defines the "vg crash" subcommand, which throws errors to test the backtrace system.
 */


#include <omp.h>
#include <unistd.h>
#include <getopt.h>
#include <stdlib.h>

#include <iostream>
#include <cassert>
#include <csignal>
#include <random>
#include <time.h>

#include "subcommand.hpp"
#include "distance.hpp"


#include "../benchmark.hpp"
#include "../utility.hpp"

using namespace std;
using namespace vg;
using namespace vg::subcommand;

void help_distcompare(char** argv){
    cerr << "usage: " << argv[0] << " distcompare [options]" << endl
         << "Output info about distance calculations" << endl
         << endl
         << "options: " << endl
         << "    -x, --xg-name, FILE       use this xg index (required)" << endl
         << "    -v, --vg-name, FILE       use this vg index (required)" << endl
         << endl;
}

int main_distcompare(int argc, char** argv){
         
    string xg_name;
    string vg_name;

    int c;
    optind = 2; // force optind past command positional argument
    while (true) {
        static struct option long_options[] =

            {
                /* These options set a flag. */
                {"help", no_argument, 0, 'h'},
                {"xg-name", required_argument, 0, 'x'},
                {"vg-name", required_argument, 0, 'v'},
                {0, 0, 0, 0}
            };

        int option_index = 0;
        c = getopt_long (argc, argv, "hx:v:",
                         long_options, &option_index);


        /* Detect the end of the options. */
        if (c == -1)
            break;

        switch (c)
        {



        case 'v':
            vg_name = optarg;
            if (vg_name.empty()) {
                cerr << "error: [vg distcompare] Must provide VG file with -v." << endl;
                exit(1);
            }
            break;

        case 'x':
            xg_name = optarg;
            if (xg_name.empty()) {
                cerr << "error: [vg distcompare] Must provide XG file with -x." << endl;
                exit(1);
            }
            break;



        case 'h':
        case '?':


        default:
            help_distcompare(argv);
            exit(1);
            break;
        }
    }

    //Check that all required arguments are given


    if (xg_name.empty()) {
        cerr << "error: [vg distcompare] Must provide XG file with -x." << endl;
        exit(1);
    }


    if (vg_name.empty()) {
        cerr << "error: [vg distcompare] Must provide VG file with -v." << endl;
        exit(1);
    }
    
    //Get file streams for xg and vg


    ifstream xg_stream(xg_name);
    if (!xg_stream) {
        cerr << "error:[vg distcompare] Cannot open XG file " << xg_name << endl;
        exit(1);
    }


    ifstream vg_stream(vg_name);
    if (!vg_stream) {
        cerr << "error:[vg distcompare] Cannot open VG file " << vg_name << endl;
        exit(1);
    }


    //Get vg and xg objects
    VG vg(vg_stream);

    xg::XG xg_index(xg_stream);

    CactusSnarlFinder bubble_finder(vg);
    SnarlManager snarl_manager = bubble_finder.find_snarls();

    clock_t t = clock();
    DistanceIndex di (&vg, &snarl_manager, 100); 
    t = clock() - t;



    int64_t size = di.sizeOf();
    
    cout << "Time to create distance index: " << t << endl;
    cout << "Space taken up by distance index: " << size << endl;

    cout << endl;
 

    random_device seed_source;
    default_random_engine generator(seed_source());

    size_t count = 0;
    while (count < 10000) {
        //Find distances between random positions in the graph
        //Outputs: my distance /t old distance /t time for my calculation /t 
        //           time for old calculation /t
        
    
        //Generate random positions by choosing snarls then nodes, then position

        size_t maxPos = xg_index.seq_length;


        size_t offset1 = uniform_int_distribution<int>(1, maxPos)(generator);
        size_t offset2 = uniform_int_distribution<int>(1, maxPos)(generator);

        vg::id_t nodeID1 = xg_index.node_at_seq_pos(offset1);
        vg::id_t nodeID2 = xg_index.node_at_seq_pos(offset2);
 
        pos_t pos1 = make_pos_t(nodeID1, false, 0);
        pos_t pos2 = make_pos_t(nodeID2, false, 0);

       
          
        count ++;

        clock_t start2 = clock();
        int64_t oldDist = abs(xg_index.closest_shared_path_oriented_distance(
                    nodeID1, 0, false, nodeID2, 0, false));
        clock_t end2 = clock();
        clock_t t2 = end2 - start2;


        clock_t start3 = clock();
//        const Snarl* snarl1 = di.snarlOf(get_id(pos1));
//        const Snarl* snarl2 = di.snarlOf(get_id(pos2));
        di.distance(pos1, pos2);
        clock_t end3 = clock();
        clock_t t3 = end3 - start3;


        clock_t start1 = clock();
        int64_t myDist = di.distance(pos1, pos2);
        clock_t end1 = clock();
        clock_t t1 = end1 - start1;
        
        cout  << myDist << "\t" << oldDist << "\t" << t1 << "\t" << t2 << "\t" << t3;
        cout << endl;
    }
    


    return 0;
}

// Register subcommand
static Subcommand vg_crash("distcompare", "compare distance calculations", DEVELOPMENT, main_distcompare);

