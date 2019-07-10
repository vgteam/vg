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
         << "    -p, --pathnames        output variant paths as SnarlTraversals to STDOUT" << endl
         << "    -r, --traversals FILE  output SnarlTraversals for ultrabubbles." << endl
         << "    -l, --leaf-only        restrict traversals to leaf ultrabubbles." << endl
         << "    -o, --top-level        restrict traversals to top level ultrabubbles" << endl
         << "    -a, --any-snarl-type   compute traversals for any snarl type (not limiting to ultrabubbles)" << endl
         << "    -m, --max-nodes N      only compute traversals for snarls with <= N nodes (with degree > 1) [10]" << endl
         << "    -t, --include-trivial  report snarls that consist of a single edge" << endl
         << "    -s, --sort-snarls      return snarls in sorted order by node ID (for topologically ordered graphs)" << endl
         << "    -v, --vcf FILE         use vcf-based instead of exhaustive traversal finder with -r" << endl
         << "    -f  --fasta FILE       reference in FASTA format (required for SVs by -v)" << endl
         << "    -i  --ins-fasta FILE   insertion sequences in FASTA format (required for SVs by -v)" << endl;
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
    bool ultrabubble_only = true;
    int max_nodes = 10;
    bool filter_trivial_snarls = true;
    bool sort_snarls = false;
    bool fill_path_names = false;
    string vcf_filename;
    string ref_fasta_filename;
    string ins_fasta_filename;

    int c;
    optind = 2; // force optind past command positional argument
    while (true) {
        static struct option long_options[] =
            {
                {"traversals", required_argument, 0, 'r'},
                {"pathnames", no_argument, 0, 'p'},
                {"leaf-only", no_argument, 0, 'l'},
                {"top-level", no_argument, 0, 'o'},
                {"any-snarl-type", no_argument, 0, 'a'},
                {"max-nodes", required_argument, 0, 'm'},
                {"include-trivial", no_argument, 0, 't'},
                {"sort-snarls", no_argument, 0, 's'},
                {"vcf", required_argument, 0, 'v'},
                {"fasta", required_argument, 0, 'f'},
                {"ins-fasta", required_argument, 0, 'i'},
                {0, 0, 0, 0}
            };

        int option_index = 0;

        c = getopt_long (argc, argv, "sr:latopm:v:f:i:h?",
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

        case 'a':
            ultrabubble_only = false;
            break;

        case 'm':
            max_nodes = parse<int>(optarg);
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
        case 'v':
            vcf_filename = optarg;
            break;
        case 'f':
            ref_fasta_filename = optarg;
            break;
        case 'i':
            ins_fasta_filename = optarg;
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
    VG* vg_graph = nullptr;
    unique_ptr<PathHandleGraph> graph;
    get_input_file(optind, argc, argv, [&](istream& in) {
            graph = vg::io::VPKG::load_one<PathHandleGraph>(in);
        });

    // this is hopefully temporary, pending traversal finder support
    vg_graph = dynamic_cast<VG*>(graph.get());
    
    // The only implemented snarl finder:
    SnarlFinder* snarl_finder = new CactusSnarlFinder(*graph);

    TraversalFinder* trav_finder = nullptr;
    vcflib::VariantCallFile variant_file;
    unique_ptr<FastaReference> ref_fasta;
    unique_ptr<FastaReference> ins_fasta;

    // vcftraversal finder is dependent on this relic of the support_caller (for now)
    map<string, PathIndex*> path_index;
    
    if (!vcf_filename.empty()) {
        if (vg_graph == nullptr) {
            cerr << "error: [vg snarls] -v requires .vg graph input" << endl;
            return 1;
        }
        variant_file.parseSamples = false;
        variant_file.open(vcf_filename);
        if (!variant_file.is_open()) {
            cerr << "error: [vg snarls] could not open " << vcf_filename << endl;
            return 1;
        }

        // load every reference path into the index
        vg_graph->paths.for_each_name([&] (const string& path_name) {
                if (!Paths::is_alt(path_name)) {
                    path_index[path_name] = new PathIndex(*vg_graph, path_name);
                }
            });

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
    auto delete_path_index = [&] () {
        for (auto name_index : path_index) {
            delete name_index.second;
        }
    };
    auto get_path_index = [&] (const Snarl& site) -> PathIndex* {
        for(auto name_index : path_index) { 
            if (name_index.second->by_id.count(site.start().node_id()) &&
                name_index.second->by_id.count(site.end().node_id())) {
                return name_index.second;
            }
        }
        return nullptr;
    };
    
    // Load up all the snarls
    SnarlManager snarl_manager = snarl_finder->find_snarls();
    vector<const Snarl*> snarl_roots = snarl_manager.top_level_snarls();
    if (fill_path_names){
        if (vg_graph == nullptr) {
            cerr << "error: [vg snarls] -p requires .vg graph input" << endl;
            return 1;
        }
        trav_finder = new PathBasedTraversalFinder(*vg_graph, snarl_manager);
        for (const Snarl* snarl : snarl_roots ){
            if (filter_trivial_snarls) {
                auto contents = snarl_manager.shallow_contents(snarl, *graph, false);
                if (contents.first.empty()) {
                    // Nothing but the boundary nodes in this snarl
                    continue;
                }
            }
            vector<SnarlTraversal> travs =  trav_finder->find_traversals(*snarl);
            vg::io::write_buffered(cout, travs, 0);
        }

        delete trav_finder;
        delete snarl_finder;
        delete_path_index();

        exit(0);
    }

    if (vcf_filename.empty()) {
        trav_finder = new ExhaustiveTraversalFinder(*graph, snarl_manager);
    } else {
        // This shuold effectively be the same as above, and is included in this tool
        // for testing purposes.  The VCFTraversalFinder differs from Exhaustive in that
        // it's easier to limit traversals using read support, and it takes care of
        // mapping back to the VCF via the alt paths. 
        trav_finder = new VCFTraversalFinder(*vg_graph, snarl_manager, variant_file, get_path_index,
                                             ref_fasta.get(), ins_fasta.get());
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
            vg::io::write_buffered(cout, snarl_buffer, buffer_size);

            auto check_max_nodes = [&graph, &max_nodes](const unordered_set<vg::id_t>& nodeset)  {
                int node_count = 0;
                for (auto node_id : nodeset) {
                    handle_t node = graph->get_handle(node_id);
                    if (graph->get_degree(node, false) > 1 || graph->get_degree(node, true) > 1) {
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
                (check_max_nodes(snarl_manager.deep_contents(snarl, *graph, true).first))) { 
                
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
    vg::io::write_buffered(cout, snarl_buffer, 0);
    if (!traversal_file.empty()) {
        vg::io::write_buffered(trav_stream, traversal_buffer, 0);
    }
    
    delete snarl_finder;
    delete trav_finder;
    delete_path_index();

    return 0;
}

// Register subcommand
static Subcommand vg_snarl("snarls", "compute snarls and their traversals", main_snarl);

