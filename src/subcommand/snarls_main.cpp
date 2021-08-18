// snarl_main.cpp: define the "vg snarls" subcommand, which outputs snarls and bubbles

#include <omp.h>
#include <unistd.h>
#include <getopt.h>

#include <list>
#include <fstream>
#include <regex>

#include "subcommand.hpp"

#include "../vg.hpp"
#include <vg/vg.pb.h>
#include "../traversal_finder.hpp"
#include "../cactus_snarl_finder.hpp"
#include "../integrated_snarl_finder.hpp"
#include "../gbwtgraph_helper.hpp"
#include <vg/io/stream.hpp>
#include <vg/io/vpkg.hpp>

//#define debug

using namespace std;
using namespace vg;
using namespace vg::subcommand;

void help_snarl(char** argv) {
    cerr << "usage: " << argv[0] << " snarls [options] graph > snarls.pb" << endl
         << "       By default, a list of protobuf Snarls is written" << endl
         << "options:" << endl
         << "    -A, --algorithm NAME   compute snarls using 'cactus' or 'integrated' algorithms (default: integrated)" << endl
         << "    -p, --pathnames        output variant paths as SnarlTraversals to STDOUT" << endl
         << "    -r, --traversals FILE  output SnarlTraversals for ultrabubbles." << endl
         << "    -e, --path-traversals  only consider traversals that correspond to paths in the graph. (-m ignored)" << endl
         << "    -l, --leaf-only        restrict traversals to leaf ultrabubbles." << endl
         << "    -o, --top-level        restrict traversals to top level ultrabubbles" << endl
         << "    -a, --any-snarl-type   compute traversals for any snarl type (not limiting to ultrabubbles)" << endl
         << "    -m, --max-nodes N      only compute traversals for snarls with <= N nodes (with degree > 1) [10]" << endl
         << "    -T, --include-trivial  report snarls that consist of a single edge" << endl
         << "    -s, --sort-snarls      return snarls in sorted order by node ID (for topologically ordered graphs)" << endl
         << "    -v, --vcf FILE         use vcf-based instead of exhaustive traversal finder with -r" << endl
         << "    -f  --fasta FILE       reference in FASTA format (required for SVs by -v)" << endl
         << "    -i  --ins-fasta FILE   insertion sequences in FASTA format (required for SVs by -v)" << endl
         << "    -t, --threads N        number of threads to use [all available]" << endl;
}

int main_snarl(int argc, char** argv) {

    if (argc == 2) {
        help_snarl(argv);
        return 1;
    }

    static const int buffer_size = 100;
    
    string algorithm = "integrated";
    string traversal_file;
    bool leaf_only = false;
    bool top_level_only = false;
    bool ultrabubble_only = true;
    int max_nodes = 10;
    bool filter_trivial_snarls = true;
    bool sort_snarls = false;
    bool fill_path_names = false;
    string vcf_filename;
    string ref_fasta_filename;
    string ins_fasta_filename;
    bool path_traversals = false;
        
    int c;
    optind = 2; // force optind past command positional argument
    while (true) {
        static struct option long_options[] =
            {
                {"algorithm", required_argument, 0, 'A'},
                {"traversals", required_argument, 0, 'r'},
                {"pathnames", no_argument, 0, 'p'},
                {"leaf-only", no_argument, 0, 'l'},
                {"top-level", no_argument, 0, 'o'},
                {"any-snarl-type", no_argument, 0, 'a'},
                {"max-nodes", required_argument, 0, 'm'},
                {"include-trivial", no_argument, 0, 'T'},
                {"sort-snarls", no_argument, 0, 's'},
                {"vcf", required_argument, 0, 'v'},
                {"fasta", required_argument, 0, 'f'},
                {"ins-fasta", required_argument, 0, 'i'},
                {"path-traversals", no_argument, 0, 'e'},                
                {"threads", required_argument, 0, 't'},
                {0, 0, 0, 0}
            };

        int option_index = 0;

        c = getopt_long (argc, argv, "A:sr:laTopm:v:f:i:eh?t:",
                         long_options, &option_index);

        /* Detect the end of the options. */
        if (c == -1)
            break;

        switch (c)
        {
        
        case 'A':
            algorithm = optarg;
            break;

        case 'r':
            traversal_file = optarg;
            break;

        case 'l':
            leaf_only = true;
            break;

        case 'o':
            top_level_only = true;
            break;

        case 'a':
            ultrabubble_only = false;
            break;

        case 'm':
            max_nodes = parse<int>(optarg);
            break;
            
        case 'T':
            filter_trivial_snarls = false;
            break;
            
        case 's':
            sort_snarls = true;
            break;
        case 'p':
            fill_path_names = true;
            break;
        case 'v':
            vcf_filename = optarg;
            break;
        case 'f':
            ref_fasta_filename = optarg;
            break;
        case 'i':
            ins_fasta_filename = optarg;
            break;
        case 'e':
            path_traversals = true;
            break;

        case 't':
        {
            int num_threads = parse<int>(optarg);
            if (num_threads <= 0) {
                cerr << "error:[vg snarls] Thread count (-t) set to " << num_threads << ", must set to a positive integer." << endl;
                exit(1);
            }
            omp_set_num_threads(num_threads);
            break;
        }
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
            cerr << "error: [vg snarls] Could not open \"" << traversal_file
                 << "\" for writing" << endl;
            return 1;
        }
    }

    // Read the graph into one of these two typed pointers.
    unique_ptr<PathHandleGraph> path_handle_graph;
    unique_ptr<gbwtgraph::GBZ> gbz;
    // And set this to point to the actual handle graph to use.
    // TODO: Should GBZ be a handle graph proxy?
    HandleGraph* handle_graph = nullptr;
    
    string graph_filename = get_input_file_name(optind, argc, argv);
    auto input = vg::io::VPKG::try_load_first<gbwtgraph::GBZ, PathHandleGraph>(graph_filename);
    if (get<0>(input)) {
        // We encountered a GBZ
        gbz = std::move(get<0>(input));
        handle_graph = &(gbz->graph);
        
        // Make sure we can still do what we wanted to do.
        // GBZ input does not work with all options.
        if (algorithm == "cactus") {
            cerr << "error: [vg snarls] the cactus algorithm does not work with GBZ graphs" << endl;
            return 1;
        }
        if (fill_path_names || path_traversals || !vcf_filename.empty()) {
            cerr << "error: [vg snarls] GBZ graphs do not have embedded paths required by -p, -r, and -v" << endl;
            return 1;
        }
    } else if (get<1>(input)) {
        // We encountered a PathHandleGraph
        path_handle_graph = std::move(get<1>(input));
        handle_graph = path_handle_graph.get();
    } else {
        // We didn't get anything.
        cerr << "error: [vg snarls] input graph is not in GBZ or PathHandleGraph format" << endl;
        return 1;
    }
    
    // Pick a SnalrFinder
    unique_ptr<SnarlFinder> snarl_finder;
    
    if (algorithm == "cactus") {
        snarl_finder.reset(new CactusSnarlFinder(*path_handle_graph));
    } else if (algorithm == "integrated") {
        snarl_finder.reset(new IntegratedSnarlFinder(*handle_graph));
    } else {
        cerr << "error:[vg snarls]: Algorithm must be 'cactus' or 'integrated', not '" << algorithm << "'" << endl;
        return 1;
    }
    if (!vcf_filename.empty() && path_traversals) {
        cerr << "error:[vg snarls]: -v cannot be used with -e" << endl;
        return 1;
    }
    if (path_traversals && traversal_file.empty()) {
        cerr << "error:[vg snarls]: -e requires -r" << endl;
        return 1;
    }
    if (!vcf_filename.empty() && traversal_file.empty()) {
        cerr << "error:[vg snarls]: -v requires -r" << endl;
        return 1;
    }

    unique_ptr<TraversalFinder> trav_finder;
    vcflib::VariantCallFile variant_file;
    unique_ptr<FastaReference> ref_fasta;
    unique_ptr<FastaReference> ins_fasta;
    
    if (!vcf_filename.empty()) {
        variant_file.parseSamples = false;
        variant_file.open(vcf_filename);
        if (!variant_file.is_open()) {
            cerr << "error: [vg snarls] could not open " << vcf_filename << endl;
            return 1;
        }

        // load up the fasta
        if (!ref_fasta_filename.empty()) {
            ref_fasta = unique_ptr<FastaReference>(new FastaReference);
            ref_fasta->open(ref_fasta_filename);
        }
        if (!ins_fasta_filename.empty()) {
            ins_fasta = unique_ptr<FastaReference>(new FastaReference);
            ins_fasta->open(ins_fasta_filename);
        }
    }
    
    // Load up all the snarls
    SnarlManager snarl_manager = snarl_finder->find_snarls_parallel();
    
    // Get all snarls in top level chains, in chain order
    vector<const Snarl*> snarl_roots;
    
    snarl_manager.for_each_top_level_chain([&](const Chain* chain) {
        // For each top level chain of 1 or more snarls
        for(auto here = chain_begin(*chain); here != chain_end(*chain); ++here) {
            // For each snarl in the chain in order
            // Remember to visit it (ignoring chain-relative orientation)
            snarl_roots.push_back(here->first);
        }
    });
    
    if (fill_path_names){
        // This finder needs a vg::VG
        trav_finder.reset(new PathBasedTraversalFinder(*path_handle_graph, snarl_manager));
        for (const Snarl* snarl : snarl_roots ){
            if (filter_trivial_snarls) {
                auto contents = snarl_manager.shallow_contents(snarl, *handle_graph, false);
                if (contents.first.empty()) {
                    // Nothing but the boundary nodes in this snarl
                    continue;
                }
            }
            vector<SnarlTraversal> travs =  trav_finder->find_traversals(*snarl);
            vg::io::write_buffered(cout, travs, 0);
        }

        exit(0);
    }

    if (path_traversals) {
        // Limit traversals to embedded paths
        trav_finder.reset(new PathTraversalFinder(*path_handle_graph, snarl_manager));
    } else if (vcf_filename.empty()) {
        // This finder works with any backing graph
        trav_finder.reset(new ExhaustiveTraversalFinder(*handle_graph, snarl_manager));
    } else {
        // This should effectively be the same as above, and is included in this tool
        // for testing purposes.  The VCFTraversalFinder differs from Exhaustive in that
        // it's easier to limit traversals using read support, and it takes care of
        // mapping back to the VCF via the alt paths.
        vector<string> ref_paths;
        path_handle_graph->for_each_path_handle([&](path_handle_t path_handle) {
            const string& name = path_handle_graph->get_path_name(path_handle);
            if (!Paths::is_alt(name)) {
              ref_paths.push_back(name);
            }
        });
        trav_finder.reset(new VCFTraversalFinder(*path_handle_graph, snarl_manager, variant_file, ref_paths,
                                                 ref_fasta.get(), ins_fasta.get()));
    }
    
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
  

    // Now we have to output stuff.
    // TODO: remove extra features and just use SnarlManager::serialize()

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
                auto contents = snarl_manager.shallow_contents(snarl, *handle_graph, false);
                if (contents.first.empty()) {
                    // Nothing but the boundary nodes in this snarl
                    continue;
                }
            }
            
            // Write our snarl tree
            snarl_buffer.push_back(*snarl);
            vg::io::write_buffered(cout, snarl_buffer, buffer_size);

            auto check_max_nodes = [&handle_graph, &max_nodes](const unordered_set<vg::id_t>& nodeset)  {
                int node_count = 0;
                for (auto node_id : nodeset) {
                    handle_t node = handle_graph->get_handle(node_id);
                    if (handle_graph->get_degree(node, false) > 1 || handle_graph->get_degree(node, true) > 1) {
                        ++node_count;
                        if (node_count > max_nodes) {
                            return false;
                        }
                    }
                }
                return true;
            };

            // Optionally write our traversals
            if (!traversal_file.empty() &&
                (!ultrabubble_only || snarl->type() == ULTRABUBBLE) &&
                (!leaf_only || snarl_manager.is_leaf(snarl)) &&
                (!top_level_only || snarl_manager.is_root(snarl)) &&
                (path_traversals || check_max_nodes(snarl_manager.deep_contents(snarl, *handle_graph, true).first))) { 
                
#ifdef debug
                cerr << "Look for traversals of " << pb2json(*snarl) << endl;
#endif
                vector<SnarlTraversal> travs = trav_finder->find_traversals(*snarl);
#ifdef debug        
                cerr << "Found " << travs.size() << endl;
#endif
                
                traversal_buffer.insert(traversal_buffer.end(), travs.begin(), travs.end());
                vg::io::write_buffered(trav_stream, traversal_buffer, buffer_size);
            }
            
            if (sort_snarls) {
                // Sort the child snarls by node ID
                vector<const Snarl*> children = snarl_manager.children_of(snarl);
                std::sort(children.begin(), children.end(), [](const Snarl* snarl_1, const Snarl* snarl_2) {
                    return snarl_1->start().node_id() < snarl_2->end().node_id();
                });
                
                for (const Snarl* child_snarl : children) {
                    stack.push_back(child_snarl);
                }
            }
            else {
                // Visit the child chains in contiguous blocks
                for (auto chain : snarl_manager.chains_of(snarl)) {
                    // For every child chain
                    for (auto here = chain_rbegin(chain); here != chain_rend(chain); ++here) {
                        // Stack up its child snarls in reverse order, so we visit them in forward order
                        stack.push_back(here->first);
                    }
                }
            }
        }
        
    }
    // flush
    vg::io::write_buffered(cout, snarl_buffer, 0);
    if (!traversal_file.empty()) {
        vg::io::write_buffered(trav_stream, traversal_buffer, 0);
    }
    
    return 0;
}

// Register subcommand
static Subcommand vg_snarl("snarls", "compute snarls and their traversals", TOOLKIT, main_snarl);

