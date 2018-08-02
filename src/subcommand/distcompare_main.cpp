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

/*
        case 'x':
            xg_name = optarg;
            if (xg_name.empty()) {
                cerr << "error: [vg distcompare] Must provide XG file with -x." << endl;
                exit(1);
            }
            break;

*/

        case 'v':
            vg_name = optarg;
            if (vg_name.empty()) {
                cerr << "error: [vg distcompare] Must provide VG file with -v." << endl;
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

/*
    if (xg_name.empty()) {
        cerr << "error: [vg distcompare] Must provide XG file with -x." << endl;
        exit(1);
    }
*/

    if (vg_name.empty()) {
        cerr << "error: [vg distcompare] Must provide VG file with -v." << endl;
        exit(1);
    }
    
    //Get file streams for xg and vg

/*
    ifstream xg_stream(xg_name);
    if (!xg_stream) {
        cerr << "error:[vg distcompare] Cannot open XG file " << xg_name << endl;
        exit(1);
    }
*/

    ifstream vg_stream(vg_name);
    if (!vg_stream) {
        cerr << "error:[vg distcompare] Cannot open VG file " << vg_name << endl;
        exit(1);
    }


    //Get vg and xg objects
    VG vg(vg_stream);

//TODO: shouldn't need to add paths

Visit visit = to_visit(1, false);
NodeTraversal n1 = to_node_traversal(visit, vg); 
list<NodeTraversal> traversal;
traversal.push_back(n1);

handle_t currHandle;

bool flag = true;
auto addVisit = [&] (const handle_t& h)-> bool {
    
    Visit v = to_visit(vg.get_id(h), vg.get_is_reverse(h));
    NodeTraversal n = to_node_traversal(v, vg);
    traversal.push_back(n1);
    currHandle = h;
    flag = true;
    return true;
};
Path p = path_from_node_traversals(traversal);
vg.include(p);
while (flag) {

    flag = false;
    vg.follow_edges(currHandle, false, addVisit);
    
}

    xg::XG xg_index(vg.graph);//xg_stream);

    CactusSnarlFinder bubble_finder(vg);
    SnarlManager snarl_manager = bubble_finder.find_snarls();

    clock_t t = clock();
    DistanceIndex di (&vg, &snarl_manager, 100); 
    t = clock() - t;
//TODO: record time to do this, and memory 
     cout << "Time to create distance index: " << t << endl;
     cout << "Space taken up by distance index: " << endl;
     cout << endl;


    vector<const Snarl*> allSnarls;
    auto addSnarl = [&] (const Snarl* s) {
        allSnarls.push_back(s);
    };
    snarl_manager.for_each_snarl_preorder(addSnarl);
   
    random_device seed_source;
    uniform_int_distribution<int> randomSnarlIndex(0, allSnarls.size() - 1);
    default_random_engine generator(seed_source());

    size_t count = 0;
    while (count < 1000) {
        //Find distances between random positions in the graph
        /*Outputs: my distance /t old distance /t time for my calculation /t 
                   time for old calculation /t
        */
    
        //Generate random positions by choosing snarls then nodes, then position
        const Snarl* snarl1 = allSnarls[randomSnarlIndex(generator)];
        const Snarl* snarl2 = allSnarls[randomSnarlIndex(generator)];
   
        pair<unordered_set<Node*>, unordered_set<Edge*>> contents1 = 
              snarl_manager.shallow_contents(snarl1, vg, true);
        pair<unordered_set<Node*>, unordered_set<Edge*>> contents2 = 
              snarl_manager.shallow_contents(snarl2, vg, true);
   
        vector<Node*> nodes1 (contents1.first.begin(), contents1.first.end());
        vector<Node*> nodes2 (contents2.first.begin(), contents2.first.end());

        uniform_int_distribution<int> randNodeIndex1(0, nodes1.size()-1);
        uniform_int_distribution<int> randNodeIndex2(0, nodes2.size()-1);
   
        Node* node1 = nodes1[randNodeIndex1(generator)];
        Node* node2 = nodes2[randNodeIndex2(generator)];
  
        vg::id_t nodeID1 = node1->id();
        vg::id_t nodeID2 = node2->id();


        vg::off_t offset1 = uniform_int_distribution<int>(0, node1->sequence().size()-1)(generator);
        vg::off_t offset2 = uniform_int_distribution<int>(0, node2->sequence().size()-1)(generator);
 
        bool rev1 = uniform_int_distribution<int>(0,1)(generator) == 0;
        bool rev2 = uniform_int_distribution<int>(0,1)(generator) == 0;
        pos_t pos1 = make_pos_t(nodeID1, rev1, offset1);
        pos_t pos2 = make_pos_t(nodeID2, rev2, offset2);

       
        //Check that positions chosen are in the snarl, not child snarl
        if (!(nodeID1 != snarl1->start().node_id() &&
            (snarl_manager.into_which_snarl(nodeID1, false) != NULL ||
             snarl_manager.into_which_snarl(nodeID1, true) != NULL)) &&
            ! (nodeID2 != snarl2->start().node_id() &&
               (snarl_manager.into_which_snarl(nodeID2, false) != NULL ||
                snarl_manager.into_which_snarl(nodeID2, true) != NULL))) {
            //Nodes aren't on child snarls
          
            count ++;

            clock_t t = clock();
            int64_t myDist = di.distance(snarl1, snarl2, pos1, pos2);
            clock_t t1 = clock() - t;
            cout << myDist << "\t";

            t = clock();
            int64_t oldDist = xg_index.closest_shared_path_oriented_distance(
                    nodeID1, offset1, rev1, nodeID2, offset2, rev2, false, 10000); 
            clock_t t2 = clock() - t;

            cout << oldDist << "\t" << t1 << "\t" << t2;
            cout << endl;
        }
    
    }
    return 0;
}

// Register subcommand
static Subcommand vg_crash("distcompare", "compare distance calculations", DEVELOPMENT, main_distcompare);

