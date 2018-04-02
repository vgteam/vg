/** \file circularize_main.cpp
 *
 * Defines the "vg circularize" subcommand
 */


#include <omp.h>
#include <unistd.h>
#include <getopt.h>

#include <iostream>

#include "subcommand.hpp"

#include "../vg.hpp"

using namespace std;
using namespace vg;
using namespace vg::subcommand;

void help_circularize(char** argv){
    cerr << "usage: " << argv[0] << " circularize [options] <graph.vg> > [circularized.vg]" << endl
        << "Makes specific paths or nodes in a graph circular." << endl
        << endl
        << "options:" << endl
        << "    -p  --path  <PATHNAME>  circularize the path by connecting its head/tail node." << endl
        << "    -P, --pathfile <PATHSFILE> circularize all paths in the provided file." << endl
        << "    -a, --head  <node_id>   circularize a head and tail node (must provide a tail)." << endl
        << "    -z, --tail  <tail_id>   circularize a head and tail node (must provide a head)." << endl
        << "    -d  --describe          list all the paths in the graph."   << endl
        << endl;
    exit(1);
}

int main_circularize(int argc, char** argv){
    if (argc == 2){
        help_circularize(argv);
        exit(1);
    }

    string path = "";
    string pathfile = "";
    bool describe = false;
    vg::id_t head = -1;
    vg::id_t tail = -1;


    int c;
    optind = 2;
    while (true){
        static struct option long_options[] =
        {
            {"path", required_argument, 0, 'p'},
            {"pathfile", required_argument, 0, 'P'},
            {"head", required_argument, 0, 'a'},
            {"tail", required_argument, 0, 'z'},
            {"describe", required_argument, 0, 'd'},
            {0,0,0,0}
        };


    int option_index = 0;
    c = getopt_long (argc, argv, "hdp:P:a:z:",
            long_options, &option_index);
    if (c == -1){
        break;
    }

        switch(c){
            case 'a':
                head = atoi(optarg);
                break;
            case 'z':
                tail = atoi(optarg);
                break;
            case 'p':
                path = optarg;
                break;
            case 'P':
                pathfile = optarg;
                break;
            case 'd':
                describe = true;
                break;
            case 'h':
            case '?':
                help_circularize(argv);
                exit(1);
                break;

            default:
                abort();
        }
    }

    vector<string> paths_to_circularize;
    if (!((head * tail) > 0)){
        cerr << "Both a head and tail node must be provided" << endl;
        help_circularize(argv);
        exit(1);
    }
    if  (pathfile != ""){
        string line;
        ifstream pfi;
        pfi.open(pathfile);
        if (!pfi.good()){
            cerr << "There is an error with the input file." << endl;
            help_circularize(argv);
        }
        while (getline(pfi, line)){
            paths_to_circularize.push_back(line);
        }
        pfi.close();

    }
    else if (path != ""){
        paths_to_circularize.push_back(path);
    }

    VG* graph;
    get_input_file(optind, argc, argv, [&](istream& in) {
        graph = new VG(in);
    });

    // Check if paths are in graph:
    for (string p : paths_to_circularize){
        bool paths_in_graph = true;
        if (!graph->paths.has_path(p)){
            cerr << "ERROR: PATH NOT IN GRAPH - " << p << endl;
            paths_in_graph = false;
        }

        if (!paths_in_graph){
            exit(1);
        }

    }

    if (describe){
       for (auto& p : graph->paths._paths){
            cout << p.first << endl;
       }
       exit(0);
    }

    if (head > 0 && tail > head){
        graph->circularize(head, tail);
    }
    else{
        graph->circularize(paths_to_circularize);
    }

    graph->serialize_to_ostream(std::cout);
    delete graph;

    return 0;
}

// Register subcommand
static Subcommand vg_circularize("circularize", "circularize a path within a graph", main_circularize);

