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

void help_circularize(char** argv) {
    cerr << "usage: " << argv[0] << " circularize [options] <graph.vg> > [circularized.vg]" << endl
         << "Make specific paths or nodes in a graph circular by connecting head/tail." << endl
         << endl
         << "options:" << endl
         << "  -p, --path NAME         circularize the path [may repeat]" << endl
         << "  -P, --pathfile FILE     circularize all paths in the provided file" << endl
         << "  -a, --head ID           circularize a head and tail node (must provide a tail)" << endl
         << "  -z, --tail ID           circularize a head and tail node (must provide a head)" << endl
         << "  -h, --help              print this help message to stderr and exit" << endl;
    exit(1);
}

int main_circularize(int argc, char** argv) {
    Logger logger("vg circularize");
    if (argc == 2){
        help_circularize(argv);
        exit(1);
    }

    vector<string> paths_to_circularize;
    string pathfile = "";
    const vg::id_t DEFAULT_ID = std::numeric_limits<nid_t>::max();
    vg::id_t head = DEFAULT_ID;
    vg::id_t tail = DEFAULT_ID;


    int c;
    optind = 2;
    while (true) {
        static struct option long_options[] =
        {
            {"help", no_argument, 0, 'h'},
            {"path", required_argument, 0, 'p'},
            {"pathfile", required_argument, 0, 'P'},
            {"head", required_argument, 0, 'a'},
            {"tail", required_argument, 0, 'z'},
            {"describe", no_argument, 0, 'd'},
            {0,0,0,0}
        };


        int option_index = 0;
        c = getopt_long (argc, argv, "h?dp:P:a:z:",
                        long_options, &option_index);
        if (c == -1) {
            break;
        }

        switch(c) {
            case 'a':
                head = parse<int>(optarg);
                break;
            case 'z':
                tail = parse<int>(optarg);
                break;
            case 'p':
                paths_to_circularize.emplace_back(optarg);
                break;
            case 'P':
                pathfile = require_exists(logger, optarg);
                break;
            case 'd':
                logger.error() << "vg circularize --describe has been removed."
                               << " Use vg paths --list" << std::endl;
                break;
            case 'h':
            case '?':
                help_circularize(argv);
                exit(1);

            default:
                abort();
        }
    }

    if ((head == DEFAULT_ID) != (tail == DEFAULT_ID)) {
        help_circularize(argv);
        logger.error() << "Both a head and tail node must be provided" << endl;
    } else if (tail < head) {
        logger.error() << "Tail " << tail << " is smaller than head " << head << endl;
    }

    if (pathfile != "") {
        string line;
        ifstream pfi;
        pfi.open(pathfile);
        if (!pfi.good()) {
            help_circularize(argv);
            logger.error() << "There is an error with the input file." << endl;
        }
        while (getline(pfi, line)) {
            paths_to_circularize.push_back(line);
        }
        pfi.close();
    }

    // TODO: if we settle on a uniform serialzation method that covers the VG class, the code is ready to be switched
    VG* graph;
    get_input_file(optind, argc, argv, [&](istream& in) {
        graph = new VG(in);
    });

    if (head != DEFAULT_ID) {
        graph->create_edge(graph->get_handle(tail), graph->get_handle(head));
    }

    for (const auto& path_name : paths_to_circularize) {
        if (!graph->has_path(path_name)) {
            logger.error() << "Path not in graph \"" << path_name << "\"" << endl;
        }

        path_handle_t path = graph->get_path_handle(path_name);
        if (graph->get_step_count(path) > 0) {
            graph->create_edge(graph->get_handle_of_step(graph->path_back(path)),
                                graph->get_handle_of_step(graph->path_begin(path)));
        }
        graph->set_circularity(path, true);
    }
    
    graph->serialize_to_ostream(cout);
//    SerializableHandleGraph* to_serialize = dynamic_cast<SerializableHandleGraph*>(&(*graph));
//    if (!to_serialize) {
//        logger.error() << "graph format is not serializable!" << endl;
//    }
//    to_serialize->serialize(std::cout);
    
    return 0;
}

// Register subcommand
static Subcommand vg_circularize("circularize", "circularize a path within a graph", main_circularize);

