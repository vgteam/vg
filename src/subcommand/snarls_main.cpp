// snarl_main.cpp: define the "vg snarls" subcommand, which outputs snarls and bubbles

#include <omp.h>
#include <unistd.h>
#include <getopt.h>

#include <list>
#include <fstream>
#include <regex>

#include "subcommand.hpp"

#include "../vg.hpp"
#include "vg.pb.h"
#include "../traversal_finder.hpp"


using namespace std;
using namespace vg;
using namespace vg::subcommand;

void help_snarl(char** argv) {
    cerr << "usage: " << argv[0] << " snarls [options] graph.vg > snarls.pb" << endl
         << "       By default, a list of protobuf Snarls is written" << endl
         << "options:" << endl
         << "    -p, --pathnames        output variant paths as SnarlTraversals to STDOUT" << endl
         << "    -r, --traversals FILE  output SnarlTraversals for ultrabubbles." << endl
         << "    -l, --leaf-only        restrict traversals to leaf ultrabubbles." << endl
         << "    -o, --top-level        restrict traversals to top level ultrabubbles" << endl
         << "    -m, --max-nodes N      only compute traversals for snarls with <= N nodes [10]" << endl
         << "    -t, --include-trivial  report snarls that consist of a single edge" << endl
         << "    -s, --sort-snarls      return snarls in sorted order by node ID (for topologically ordered graphs)" << endl;
}

int main_snarl(int argc, char** argv) {

    if (argc == 2) {
        help_snarl(argv);
        return 1;
    }

    static const int buffer_size = 100;
    
    string traversal_file;
    bool leaf_only = false;
    bool top_level_only = false;
    int max_nodes = 10;
    bool filter_trivial_snarls = true;
    bool sort_snarls = false;
    bool fill_path_names = false;

    int c;
    optind = 2; // force optind past command positional argument
    while (true) {
        static struct option long_options[] =
            {
                {"traversals", required_argument, 0, 'r'},
		        {"pathnames", no_argument, 0, 'p'},
                {"leaf-only", no_argument, 0, 'l'},
                {"top-level", no_argument, 0, 'o'},
                {"max-nodes", required_argument, 0, 'm'},
                {"include-trivial", no_argument, 0, 't'},
                {"sort-snarls", no_argument, 0, 's'},
                {0, 0, 0, 0}
            };

        int option_index = 0;

        c = getopt_long (argc, argv, "sr:ltopm:h?",
                         long_options, &option_index);

        /* Detect the end of the options. */
        if (c == -1)
            break;

        switch (c)
        {
            
        case 'r':
            traversal_file = optarg;
            break;

        case 'l':
            leaf_only = true;
            break;

        case 'o':
            top_level_only = true;
            break;

        case 'm':
            max_nodes = atoi(optarg);
            break;
            
        case 't':
            filter_trivial_snarls = false;
            break;
            
        case 's':
            sort_snarls = true;
            break;
        case 'p':
            fill_path_names = true;
            break;
            
        case 'h':
        case '?':
            /* getopt_long already printed an error message. */
            help_snarl(argv);
            exit(1);
            break;

        default:
            abort ();

        }
    }

    // Prepare traversal output stream
    ofstream trav_stream;
    if (!traversal_file.empty()) {
        trav_stream.open(traversal_file);
        if (!trav_stream) {
            cerr << "error:[vg snarl]: Could not open \"" << traversal_file
                 << "\" for writing" << endl;
            return 1;
        }
    }

    // Read the graph
    VG* graph;
    get_input_file(optind, argc, argv, [&](istream& in) {
            graph = new VG(in);
        });
    
    if (graph == nullptr) {
        cerr << "error:[vg snarl]: Could not load graph" << endl;
        exit(1);
    }

    // The only implemented snarl finder:
    SnarlFinder* snarl_finder = new CactusSnarlFinder(*graph);
    
    // Load up all the snarls
    SnarlManager snarl_manager = snarl_finder->find_snarls();
    vector<const Snarl*> snarl_roots = snarl_manager.top_level_snarls();
    if (fill_path_names){
        TraversalFinder* trav_finder = new PathBasedTraversalFinder(*graph, snarl_manager);
        for (const Snarl* snarl : snarl_roots ){
            if (filter_trivial_snarls) {
                auto contents = snarl_manager.shallow_contents(snarl, *graph, false);
                if (contents.first.empty()) {
                    // Nothing but the boundary nodes in this snarl
                    continue;
                }
            }
            vector<SnarlTraversal> travs =  trav_finder->find_traversals(*snarl);
            stream::write_buffered(cout, travs, 0);
        }

        delete trav_finder;
        delete snarl_finder;
        delete graph;

        exit(0);
    }


    TraversalFinder* trav_finder = new ExhaustiveTraversalFinder(*graph, snarl_manager);
    
    // Sort the top level Snarls
    if (sort_snarls) {
        // Ensure that all snarls are stored in sorted order
        list<const Snarl*> snarl_stack;
        for (const Snarl* root : snarl_roots) {
            snarl_stack.push_back(root);
            while (!snarl_stack.empty()) {
                const Snarl* snarl = snarl_stack.back();
                snarl_stack.pop_back();
                if (snarl->start().node_id() > snarl->end().node_id()) {
                    snarl_manager.flip(snarl);
                }
                for (const Snarl* child_snarl : snarl_manager.children_of(snarl)) {
                    snarl_stack.push_back(child_snarl);
                }
            }
        }
        
        // Sort the snarls by node ID
        std::sort(snarl_roots.begin(), snarl_roots.end(), [](const Snarl* snarl_1, const Snarl* snarl_2) {
            return snarl_1->start().node_id() < snarl_2->end().node_id();
        });
    }

  

    // Protobuf output buffers
    vector<Snarl> snarl_buffer;
    vector<SnarlTraversal> traversal_buffer;
    
    list<const Snarl*> stack;

    for (const Snarl* root : snarl_roots) {
        
        stack.push_back(root);
        
        while (!stack.empty()) {
            const Snarl* snarl = stack.back();
            stack.pop_back();
            
            if (filter_trivial_snarls) {
                auto contents = snarl_manager.shallow_contents(snarl, *graph, false);
                if (contents.first.empty()) {
                    // Nothing but the boundary nodes in this snarl
                    continue;
                }
            }
            
            // Write our snarl tree
            snarl_buffer.push_back(*snarl);
            stream::write_buffered(cout, snarl_buffer, buffer_size);
            
            // Optionally write our traversals
            if (!traversal_file.empty() && snarl->type() == ULTRABUBBLE &&
                (!leaf_only || snarl_manager.is_leaf(snarl)) &&
                (!top_level_only || snarl_manager.is_root(snarl)) &&
                (snarl_manager.deep_contents(snarl, *graph, true).first.size() < max_nodes)) {
                
#ifdef debug
                cerr << "Look for traversals of " << pb2json(*snarl) << endl;
#endif
                vector<SnarlTraversal> travs = trav_finder->find_traversals(*snarl);
#ifdef debug        
                cerr << "Found " << travs.size() << endl;
#endif
                
                traversal_buffer.insert(traversal_buffer.end(), travs.begin(), travs.end());
                stream::write_buffered(trav_stream, traversal_buffer, buffer_size);
            }
            
            // Sort the child snarls by node ID?
            if (sort_snarls) {
                vector<const Snarl*> children = snarl_manager.children_of(snarl);
                std::sort(children.begin(), children.end(), [](const Snarl* snarl_1, const Snarl* snarl_2) {
                    return snarl_1->start().node_id() < snarl_2->end().node_id();
                });
                
                for (const Snarl* child_snarl : children) {
                    stack.push_back(child_snarl);
                }
            }
            else {
                for (const Snarl* child_snarl : snarl_manager.children_of(snarl)) {
                    stack.push_back(child_snarl);
                }
            }
        }
        
        
    }
    // flush
    stream::write_buffered(cout, snarl_buffer, 0);
    if (!traversal_file.empty()) {
        stream::write_buffered(trav_stream, traversal_buffer, 0);
    }
    
    delete snarl_finder;
    delete trav_finder;
    delete graph;

    return 0;
}

// Register subcommand
static Subcommand vg_snarl("snarls", "compute snarls and their traversals", main_snarl);

