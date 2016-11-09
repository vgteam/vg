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
    optind = 2; // force optind past command positional argument
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
        if (show_progress) {
            cerr << "Reading standard input..." << endl;
        }
        graph = new VG(std::cin);
    } else {
        if (show_progress) {
            cerr << "Reading " << file_name << "..." << endl;
        }
        ifstream in;
        in.open(file_name.c_str());
        graph = new VG(in);
    }
    
    // Set the progress flag
    graph->show_progress = show_progress;
    
    // We need this to get the bubble tree
    CactusSiteFinder site_finder(*graph, "");
    
    // We need this to find traversals for sites.
    TrivialTraversalFinder traversal_finder(*graph);
    
    // Make a list of leaf sites
    list<NestedSite> leaves;
    
    if (show_progress) {
        cerr << "Scanning " << graph->node_count() << " nodes for sites..." << endl;
    }
        
    
    site_finder.for_each_site_parallel([&](NestedSite root) {
        // For every tree of sites
        
        // We keep a queue of sites to process.
        list<NestedSite*> to_check{&root};
        
#ifdef debug
        #pragma omp critical (cerr)
        cerr << "Found root site with " << root.children.size() << " children" << endl;
#endif
        
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
    
    if (show_progress) {
        cerr << "Found " << leaves.size() << " leaves" << endl;
    }
    
    // Now we have a list of all the leaf sites.
    graph->create_progress("simplify leaves", leaves.size());
    
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
        
#ifdef debug
        cerr << "Found " << total_size << " bp leaf" << endl;
        for (auto* node : leaf.nodes) {
            cerr << "\t" << node->id() << ": " << node->sequence() << endl;
        }
#endif
        
        if (total_size == 0) {
            // This site is just the start and end nodes, so it doesn't make
            // sense to try and remove it.
            continue;
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
        
        // Find all the paths that traverse this region.
        
        // We start at the start node. Copy out all the mapping pointers on that
        // node, so we can go through them while tampering with them.
        map<string, set<Mapping*> > mappings_by_path = graph->paths.get_node_mapping(leaf.start.node);
        
        for (auto& kv : mappings_by_path) {
            // For each path that hits the start node
            
            // Unpack the name
            auto& path_name = kv.first;
            
            for (Mapping* start_mapping : kv.second) {
                // For each visit to the start node
                
                // Determine what orientation we're going to scan in
                bool backward = start_mapping->position().is_reverse();
                
                // We're going to fill this list with the mappings we need to
                // remove and replace in this path for this traversal. Initially
                // runs from start of site to end of site, but later gets
                // flipped into path-local orientation.
                list<Mapping*> existing_mappings;
                
                // Tracing along forward/backward from each as appropriate, see
                // if the end of the site is found in the expected orientation
                // (or if the path ends first).
                bool found_end = false;
                Mapping* here = start_mapping;
                while (here) {
                    // Until we hit the start/end of the path or the mapping we want
                    
                    // Remember the mapping so we can remove it later.
                    existing_mappings.push_back(here);
                    if (here->position().node_id() == leaf.end.node->id() &&
                        here->position().is_reverse() == (leaf.end.backward != backward)) {
                        // We have encountered the end of the site in the
                        // orientation we expect, given the orientation we saw
                        // for the start.
                        
                        found_end = true;
                        
                        // Stop scanning!
                        break;
                    }
                    
                    // Scan left along ther path if we found the site start backwards, and right if we found it forwards.
                    here = backward ? graph->paths.traverse_left(here) : graph->paths.traverse_right(here);
                }
                
                if (!found_end) {
                    // This path only partly traverses the site.
                                 
                    for(auto* mapping : existing_mappings) {
                        // Trim the path out of the site
                        graph->paths.remove_mapping(mapping);
                    }
                          
                    // Maybe the next time the path visits the site it will go
                    // all the way through.
                    continue;
                }
                
                // If we found the end, remove all the mappings encountered, in
                // order so that the last one removed is the last one along the
                // path.
                if (backward) {
                    // Make sure the last mapping in the list is the last
                    // mapping to occur along the path.
                    existing_mappings.reverse();
                }
                
                // Where will we insert the new site traversal into the path?
                list<Mapping>::iterator insert_position;
                for (auto* mapping : existing_mappings) {
                    // Remove each mapping from left to right along the path,
                    // saving the position after the mapping we just removed. At
                    // the end we'll have the position after the whole site
                    // traversal.
                    
#ifdef debug
                    cerr << path_name << ": Drop mapping " << pb2json(*mapping) << endl;
#endif
                    
                    insert_position = graph->paths.remove_mapping(mapping);
                }
                
                // Then insert mappings for the official traversal we picked,
                // from right to left in the path's local orientation.
                for (auto i = visits.rbegin(); i != visits.rend(); ++i) {
                    // For each visit in the official traversal for this site, right to left
                    
                    // Make a Mapping to represent it
                    Mapping new_mapping;
                    new_mapping.mutable_position()->set_node_id(i->node->id());
                    new_mapping.mutable_position()->set_is_reverse(i->backward);
                    
                    // Add an edit
                    Edit* edit = new_mapping.add_edit();
                    edit->set_from_length(i->node->sequence().size());
                    edit->set_to_length(i->node->sequence().size());
                    
#ifdef debug
                    cerr << path_name << ": Add mapping " << pb2json(new_mapping) << endl;
#endif
                    
                    // Insert the mapping in the path, moving right to left
                    insert_position = graph->paths.insert_mapping(insert_position, path_name, new_mapping);
                }
                
            }
        }
           
        // Now delete all edges that aren't connecting adjacent nodes on the
        // blessed traversal (before we delete their nodes).
        set<Edge*> blessed_edges;
        for (auto i = visits.begin(); i != --visits.end(); ++i) {
            // For each node and the next node (which won't be the end)
            auto next = i;
            next++;
            
            assert(next != visits.end());
            
            // Find the edge between them
            NodeTraversal here(i->node, i->backward);
            NodeTraversal next_traversal(next->node, next->backward);
            Edge* edge = graph->get_edge(here, next_traversal);
            assert(edge != nullptr);
            
            // Remember we need it
            blessed_edges.insert(edge);
        }
        
        for (auto* edge : leaf.edges) {
            if (!blessed_edges.count(edge)) {
                // Get rid of all the edges not needed for the one true traversal
                graph->destroy_edge(edge);
            }
        }
       
           
        // Now delete all the nodes that aren't on the blessed traversal.
        
        // What nodes are on it?
        set<Node*> blessed_nodes;
        for (auto& visit : visits) {
            blessed_nodes.insert(visit.node);
        }
        
        for (auto* node : leaf.nodes) {
            // For every node in the site
            if (!blessed_nodes.count(node)) {
                // If we don't need it for the chosen path, destroy it
                graph->destroy_node(node);
            }
        }
        
        // OK we finished a leaf
        graph->increment_progress();
    }
    
    graph->destroy_progress();
    
    // Reset the ranks in the graph, since we rewrote paths
    graph->paths.clear_mapping_ranks();
    
    // Serialize the graph
    graph->serialize_to_ostream(std::cout);
    
    delete graph;

    // NB: If you worry about "still reachable but possibly lost" warnings in valgrind,
    // this would free all the memory used by protobuf:
    //ShutdownProtobufLibrary();

    return 0;
}

// Register subcommand
static Subcommand vg_construct("simplify", "graph simplification", main_simplify);

