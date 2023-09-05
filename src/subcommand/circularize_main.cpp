/** \file circularize_main.cpp
 *
 * Defines the "vg circularize" subcommand
 */


#include <omp.h>
#include <unistd.h>
#include <getopt.h>

#include <iostream>

#include "subcommand.hpp"

#include "../utility.hpp"
#include "../handle.hpp"
#include "../vg.hpp"
#include <vg/io/stream.hpp>
#include <vg/io/vpkg.hpp>

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
                head = parse<int>(optarg);
                break;
            case 'z':
                tail = parse<int>(optarg);
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

    // TODO: if we settle on a uniform serialzation method that covers the VG class, the code is ready to be switched
    VG* graph;
    get_input_file(optind, argc, argv, [&](istream& in) {
        graph = new VG(in);
    });

    // Check if paths are in graph:
    for (const string& p : paths_to_circularize){
        if (!graph->has_path(p)){
            cerr << "ERROR: PATH NOT IN GRAPH - " << p << endl;
            exit(1);
        }
    }

    if (describe){
        graph->for_each_path_handle([&](const path_handle_t& path_handle) {
            cout << graph->get_path_name(path_handle) << endl;
        });
       exit(0);
    }

    if (head > 0 && tail > head){
        graph->create_edge(graph->get_handle(tail), graph->get_handle(head));
    }
    else{
        for (const auto& path_name : paths_to_circularize) {
            path_handle_t path = graph->get_path_handle(path_name);
            if (graph->get_step_count(path) > 0) {
                graph->create_edge(graph->get_handle_of_step(graph->path_back(path)),
                                   graph->get_handle_of_step(graph->path_begin(path)));
            }
            graph->set_circularity(path, true);
        }
    }
    
    graph->serialize_to_ostream(cout);
//    SerializableHandleGraph* to_serialize = dynamic_cast<SerializableHandleGraph*>(&(*graph));
//    if (!to_serialize) {
//        cerr << "error: graph format is not serializable!" << endl;
//        return 1;
//    }
//    to_serialize->serialize(std::cout);
    
    return 0;
}

// Register subcommand
static Subcommand vg_circularize("circularize", "circularize a path within a graph", main_circularize);

