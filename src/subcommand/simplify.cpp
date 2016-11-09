// simplify.cpp: define the "vg simplify" subcommand, which removes small variation

#include <omp.h>
#include <unistd.h>
#include <getopt.h>

#include <list>

#include "subcommand.hpp"

#include "../vg.hpp"
// This provides the CactusSiteFinder
#include "../genotypekit.hpp"


using namespace std;
using namespace vg;
using namespace vg::subcommand;

void help_simplify(char** argv) {
    cerr << "usage: " << argv[0] << " simplify [options] <old.vg >new.vg" << endl
         << "options:" << endl
         << "    -m, --min-size N      remove bubbles with fewer than N bases involved (default: 10)" << endl
         << "    -p, --progress        show progress" << endl
         << "    -t, --threads N       use N threads to construct graph (defaults to numCPUs)" << endl;
}

int main_simplify(int argc, char** argv) {

    if (argc == 2) {
        help_simplify(argv);
        return 1;
    }


    size_t min_size = 10;
    bool show_progress = false;

    int c;
    while (true) {
        static struct option long_options[] =
            {
                {"min-size", required_argument, 0, 'm'},
                {"progress",  no_argument, 0, 'p'},
                {"threads", required_argument, 0, 't'},
                {"help", no_argument, 0, 'h'},
                {0, 0, 0, 0}
            };

        int option_index = 0;
        c = getopt_long (argc, argv, "m:pt:h?",
                         long_options, &option_index);

        /* Detect the end of the options. */
        if (c == -1)
            break;

        switch (c)
        {

        case 'm':
            min_size = atoi(optarg);
            break;

        case 'p':
            show_progress = true;
            break;

        case 't':
            omp_set_num_threads(atoi(optarg));
            break;

        case 'h':
        case '?':
            /* getopt_long already printed an error message. */
            help_simplify(argv);
            exit(1);
            break;

        default:
            abort ();

        }
    }
    
    // TODO: move all this to a simplifier object
    
    // Load the graph
    VG* graph;
    string file_name = argv[optind];
    if (file_name == "-") {
        graph = new VG(std::cin);
    } else {
        ifstream in;
        in.open(file_name.c_str());
        graph = new VG(in);
    }
    
    // We need this to get the bubble tree
    CactusSiteFinder site_finder(*graph, "");
    
    // We need this to find traversals for sites.
    TrivialTraversalFinder traversal_finder(*graph);
    
    // Make a list of leaf sites
    list<NestedSite> leaves;
    
    site_finder.for_each_site_parallel([&](NestedSite root) {
        // For every tree of sites
        
        // We keep a queue of sites to process.
        list<NestedSite*> to_check{&root};
        
        while (!to_check.empty()) {
            // Until we have seen all the sites, check one
            NestedSite& check = *(to_check.front());
            to_check.pop_front();
            
            if (check.children.empty()) {
                // This is a leaf. Copy it into the leaves queue.
                #pragma omp critical (leaves)
                leaves.push_back(check);
            } else {
                for (auto& child : check.children) {
                    // Check all the children, depth-first.
                    to_check.push_front(&child);
                }
            }
        }
        
    });
    
    // Now we have a list of all the leaf sites.
    
    for (auto& leaf : leaves) {
        // Look at all the leaves
        
        // For each leaf, calculate its total size.
        size_t total_size = 0;
        for (auto* node : leaf.nodes) {
            // For each node
            if (node == leaf.start.node || node == leaf.end.node) {
                // That isn't a start or end
                continue;
            }
            // Include it in the size figure
            total_size += node->sequence().size();
        }
        
        if (total_size >= min_size) {
            // This site is too big to remove
            continue;
        }
        
        // Otherwise we want to simplify this site away
        
        // Identify the replacement traversal for the bubble
        auto traversals = traversal_finder.find_traversals(leaf);
        
        if (traversals.empty()) {
            // We couldn't find any paths through the site.
            continue;
        }
        
        // Get the collection of visits in the traversal we want to keep
        auto& visits = traversals.front().visits;
        
        // Now we have to rewrite paths that visit nodes/edges not on this
        // traversal, or in a different order, or whatever. To be safe we'll
        // just rewrite all paths.
        
        // Find all the paths that visit nodes in this region
        
        // For each path
            // Find all the places it visits the start node of the site.
            // For each, determine what orientation we're going to scan in
            // Tracing along forward/backward from each as appropriate, see if the end of the site is found in the expected orientation (or if the path ends first).
            // If we found the end, remove all the mappings encountered.
            // Then insert mappings for the official traversal we picked, in the appropriate orientation.
            
        // Now delete all the nodes that aren't on the blessed traversal.
        
        // Now delete all edges that aren't connecting adjacent nodes on the blessed traversal.
    }
    
    // Reset the ranks in the graph, since we rewrote paths
    graph->paths.clear_mapping_ranks();
    
    // Serialize the graph

    // NB: If you worry about "still reachable but possibly lost" warnings in valgrind,
    // this would free all the memory used by protobuf:
    //ShutdownProtobufLibrary();

    return 0;
}

// Register subcommand
static Subcommand vg_construct("simplify", "graph simplification", main_simplify);

