// mod.cpp: define the "vg mod" subcommand, which modifies vg graphs

#include <omp.h>
#include <unistd.h>
#include <getopt.h>

#include <string>
#include <vector>
#include <regex>

#include "subcommand.hpp"

#include "../vg.hpp"
#include "../cactus.hpp"
#include <vg/io/stream.hpp>
#include <vg/io/vpkg.hpp>
#include <handlegraph/mutable_path_deletable_handle_graph.hpp>
#include "../handle.hpp"
#include "../utility.hpp"
#include "../algorithms/simplify_siblings.hpp"
#include "../algorithms/normalize.hpp"
#include "../algorithms/prune.hpp"
#include "../io/save_handle_graph.hpp"

using namespace std;
using namespace vg;
using namespace vg::subcommand;

void help_mod(char** argv) {
    cerr << "usage: " << argv[0] << " mod [options] <graph.vg> >[mod.vg]" << endl
         << "Modifies graph, outputs modified on stdout." << endl
         << endl
         << "options:" << endl
         << "    -P, --label-paths       don't edit with -i alignments, just use them for labeling the graph" << endl
         << "    -c, --compact-ids       should we sort and compact the id space? (default false)" << endl
         << "    -b, --break-cycles      use an approximate topological sort to break cycles in the graph" << endl
         << "    -n, --normalize         normalize the graph so that edges are always non-redundant" << endl
         << "                            (nodes have unique starting and ending bases relative to neighbors," << endl
         << "                            and edges that do not introduce new paths are removed and neighboring" << endl
         << "                            nodes are merged)" << endl
         << "    -U, --until-normal N    iterate normalization until convergence, or at most N times" << endl
         << "    -z, --nomerge-pre STR   do not let normalize (-n, -U) zip up any pair of nodes that both belong to path with prefix STR" << endl
         << "    -E, --unreverse-edges   flip doubly-reversing edges so that they are represented on the" << endl
         << "                            forward strand of the graph" << endl
         << "    -s, --simplify          remove redundancy from the graph that will not change its path space" << endl
         << "    -d, --dagify-step N     copy strongly connected components of the graph N times, forwarding" << endl
         << "                            edges from old to new copies to convert the graph into a DAG" << endl
         << "    -w, --dagify-to N       copy strongly connected components of the graph forwarding" << endl
         << "                            edges from old to new copies to convert the graph into a DAG" << endl
         << "                            until the shortest path through each SCC is N bases long" << endl
         << "    -L, --dagify-len-max N  stop a dagification step if the unrolling component has this much sequence" << endl
         << "    -f, --unfold N          represent inversions accessible up to N from the forward" << endl
         << "                            component of the graph" << endl
         << "    -O, --orient-forward    orient the nodes in the graph forward" << endl
         << "    -N, --remove-non-path   keep only nodes and edges which are part of paths" << endl
         << "    -A, --remove-path       keep only nodes and edges which are not part of any path" << endl
         << "    -k, --keep-path NAME    keep only nodes and edges in the path" << endl
         << "    -R, --remove-null       removes nodes that have no sequence, forwarding their edges" << endl
         << "    -g, --subgraph ID       gets the subgraph rooted at node ID, multiple allowed" << endl
         << "    -x, --context N         steps the subgraph out by N steps (default: 1)" << endl
         << "    -p, --prune-complex     remove nodes that are reached by paths of --length which" << endl
         << "                            cross more than --edge-max edges" << endl
         << "    -S, --prune-subgraphs   remove subgraphs which are shorter than --length" << endl
         << "    -l, --length N          for pruning complex regions and short subgraphs" << endl
         << "    -X, --chop N            chop nodes in the graph so they are not more than N bp long" << endl
         << "    -u, --unchop            where two nodes are only connected to each other and by one edge" << endl
         << "                            replace the pair with a single node that is the concatenation of their labels" << endl
         << "    -e, --edge-max N        only consider paths which make edge choices at <= this many points" << endl
         << "    -M, --max-degree N      unlink nodes that have edge degree greater than N" << endl
         << "    -m, --markers           join all head and tails nodes to marker nodes" << endl
         << "                            ('###' starts and '$$$' ends) of --length, for debugging" << endl
         << "    -y, --destroy-node ID   remove node with given id" << endl
         << "    -a, --cactus            convert to cactus graph representation" << endl
         << "    -v, --sample-vcf FILE   for a graph with allele paths, compute the sample graph from the given VCF" << endl
         << "    -G, --sample-graph FILE subset an augmented graph to a sample graph using a Locus file" << endl
         << "    -t, --threads N         for tasks that can be done in parallel, use this many threads" << endl;
}

int main_mod(int argc, char** argv) {

    if (argc == 2) {
        help_mod(argv);
        return 1;
    }

    string path_name;
    bool label_paths = false;
    bool compact_ids = false;
    bool prune_complex = false;
    int path_length = 0;
    int edge_max = 0;
    int chop_to = 0;
    bool add_start_and_end_markers = false;
    bool prune_subgraphs = false;
    bool simplify_graph = false;
    bool unchop = false;
    bool normalize_graph = false;
    bool remove_non_path = false;
    bool remove_path = false;
    bool compact_ranks = false;
    vector<nid_t> root_nodes;
    int32_t context_steps;
    bool remove_null = false;
    uint32_t unfold_to = 0;
    bool break_cycles = false;
    uint32_t dagify_steps = 0;
    uint32_t dagify_to = 0;
    uint32_t dagify_component_length_max = 0;
    bool orient_forward = false;
    nid_t destroy_node_id = 0;
    int until_normal_iter = 0;
    bool flip_doubly_reversed_edges = false;
    bool cactus = false;
    string vcf_filename;
    string loci_filename;
    int max_degree = 0;
    string nomerge_prefix;

    int c;
    optind = 2; // force optind past command positional argument
    while (true) {
        static struct option long_options[] =

        {
            {"help", no_argument, 0, 'h'},
            {"include-aln", required_argument, 0, 'i'},
            {"include-loci", required_argument, 0, 'q'},
            {"include-gt", required_argument, 0, 'Q'},
            {"compact-ids", no_argument, 0, 'c'},
            {"compact-ranks", no_argument, 0, 'C'},
            {"keep-path", required_argument, 0, 'k'},
            {"remove-orphans", no_argument, 0, 'o'},
            {"prune-complex", no_argument, 0, 'p'},
            {"prune-subgraphs", no_argument, 0, 'S'},
            {"length", required_argument, 0, 'l'},
            {"edge-max", required_argument, 0, 'e'},
            {"chop", required_argument, 0, 'X'},
            {"markers", no_argument, 0, 'm'},
            {"threads", no_argument, 0, 't'},
            {"label-paths", no_argument, 0, 'P'},
            {"simplify", no_argument, 0, 's'},
            {"unchop", no_argument, 0, 'u'},
            {"normalize", no_argument, 0, 'n'},
            {"until-normal", required_argument, 0, 'U'},
            {"nomerge-pre", required_argument, 0, 'z'},
            {"remove-non-path", no_argument, 0, 'N'},
            {"remove-path", no_argument, 0, 'A'},
            {"orient-forward", no_argument, 0, 'O'},
            {"unfold", required_argument, 0, 'f'},
            {"subgraph", required_argument, 0, 'g'},
            {"context", required_argument, 0, 'x'},
            {"remove-null", no_argument, 0, 'R'},
            {"dagify-steps", required_argument, 0, 'd'},
            {"dagify-to", required_argument, 0, 'w'},
            {"dagify-len-max", required_argument, 0, 'L'},
            {"break-cycles", no_argument, 0, 'b'},
            {"destroy-node", required_argument, 0, 'y'},
            {"translation", required_argument, 0, 'Z'},
            {"unreverse-edges", required_argument, 0, 'E'},
            {"cactus", no_argument, 0, 'a'},
            {"sample-vcf", required_argument, 0, 'v'},
            {"sample-graph", required_argument, 0, 'G'},
            {"max-degree", required_argument, 0, 'M'},
            {"drop-paths", no_argument, 0, 'D'},
            {"retain-path", required_argument, 0, 'r'},
            {"retain-complement", no_argument, 0, 'I'},
            {0, 0, 0, 0}
        };

        int option_index = 0;
        c = getopt_long (argc, argv, "hk:oi:q:Q:cpl:e:mt:SX:KPsunz:NAf:Cg:x:RTU:Bbd:Ow:L:y:Z:Eav:G:M:Dr:I",
                long_options, &option_index);


        // Detect the end of the options.
        if (c == -1)
            break;

        switch (c)
        {

        case 'i':
            cerr << "[vg mod] error: vg mod -i is deprecated.  please switch to vg augment" << endl;
            exit(1);

        case 'q':
            cerr << "[vg mod] error: vg mod -q is deprecated.  please switch to vg augment -l" << endl;
            exit(1);

        case 'Q':
            cerr << "[vg mod] error: vg mod -Q is deprecated.  please switch to vg augment -L" << endl;
            exit(1);
            break;

        case 'Z':
            cerr << "[vg mod] error: vg mod -Z is deprecated.  please switch to vg augment -Z" << endl;
            exit(1);
            break;

        case 'D':
            cerr << "[vg mod] error: vg mod -D is deprecated.  please switch to vg paths -d" << endl;
            exit(1);
            break;

        case 'r':
            cerr << "[vg mod] error: vg mod -r is deprecated.  please switch to vg paths -r" << endl;
            exit(1);
            break;

        case 'I':
            cerr << "[vg mod] error: vg mod -I is deprecated.  please switch to vg paths -d" << endl;
            exit(1);
            break;

        case 'c':
            compact_ids = true;
            break;

        case 'k':
            path_name = optarg;
            break;

        case 'o':
            cerr << "warning[vg mod]: -o is deprecated. Dangling edges are now automatically removed." << endl;
            break;

        case 'p':
            prune_complex = true;
            break;

        case 'S':
            prune_subgraphs = true;
            break;

        case 'l':
            path_length = parse<int>(optarg);
            break;

        case 'X':
            chop_to = parse<int>(optarg);
            break;

        case 'u':
            unchop = true;
            break;

        case 'E':
            flip_doubly_reversed_edges = true;
            break;

        case 'e':
            edge_max = parse<int>(optarg);
            break;

        case 'm':
            add_start_and_end_markers = true;
            break;

        case 't':
            omp_set_num_threads(parse<int>(optarg));
            break;

        case 'f':
            unfold_to = parse<int>(optarg);
            break;

        case 'O':
            orient_forward = true;
            break;

        case 'P':
            cerr << "[vg mod] warning: vg mod -P is deprecated and will soon be removed.  please switch to vg augment -B" << endl;
            label_paths = true;
            break;

        case 's':
            simplify_graph = true;
            break;

        case 'n':
            normalize_graph = true;
            break;

        case 'z':
            nomerge_prefix = optarg;
            break;

        case 'N':
            remove_non_path = true;
            break;
            
        case 'A':
            remove_path = true;
            break;

        case 'U':
            until_normal_iter = parse<int>(optarg);
            break;

        case 'd':
            dagify_steps = parse<int>(optarg);
            break;

        case 'w':
            dagify_to = parse<int>(optarg);
            break;

        case 'L':
            dagify_component_length_max = parse<int>(optarg);
            break;

        case 'b':
            break_cycles = true;
            break;

        case 'g':
            root_nodes.push_back(parse<int>(optarg));
            break;

        case 'x':
            context_steps = parse<int>(optarg);
            break;

        case 'R':
            remove_null = true;
            break;

        case 'y':
            destroy_node_id = parse<int>(optarg);
            break;

        case 'a':
            cactus = true;
            break;

        case 'v':
            vcf_filename = optarg;
            break;
            
        case 'G':
            loci_filename = optarg;
            break;

        case 'M':
            max_degree = parse<int>(optarg);
            break;

        case 'h':
        case '?':
            help_mod(argv);
            exit(1);
            break;

        default:
            abort ();
        }
    }

    unique_ptr<handlegraph::MutablePathDeletableHandleGraph> graph;
    string graph_filename = get_input_file_name(optind, argc, argv);
    graph = vg::io::VPKG::load_one<handlegraph::MutablePathDeletableHandleGraph>(graph_filename);

    if (!vcf_filename.empty()) {
        // We need to throw out the parts of the graph that are on alt paths,
        // but not on alt paths for alts used by the first sample in the VCF.

        // This is called with the entire path name string to detect alt
        // paths.
        const function<bool(const string&)>& is_alt = Paths::is_alt;

        // This holds the VCF file we read the variants from. It needs to be the
        // same one used to construct the graph.
        vcflib::VariantCallFile variant_file;
        variant_file.open(vcf_filename);
        if (!variant_file.is_open()) {
            cerr << "error:[vg mod] could not open" << vcf_filename << endl;
            return 1;
        }

        // Now go through and prune down the varaints.

        // How many phases are there?
        size_t num_samples = variant_file.sampleNames.size();
        // TODO: we can only handle single-sample VCFs
        assert(num_samples == 1);

        // This will hold the IDs of all nodes visited by alt paths that aren't used.
        set<vg::id_t> alt_path_ids;

        graph->for_each_path_handle([&](const path_handle_t& p) {
            auto name = graph->get_path_name(p);
            // For every path name in the graph

            if(is_alt(name)) {
                // If it's an alt path, walk it
                
                graph->for_each_step_in_path(p, [&](const step_handle_t& s) {
                    // Mark all nodes that are part of it as on alt paths
                    alt_path_ids.insert(graph->get_id(graph->get_handle_of_step(s)));
                });
            }
        });

        // We also have a function to handle each variant as it comes in.
        auto handle_variant = [&](vcflib::Variant& variant) {
            // So we have a variant

            if(variant.alleles.size() < 2) {
                // Skip non-variable variants.
                return;
            }

            // Grab its id, or make one by hashing stuff if it doesn't
            // have an ID.
            string var_name = make_variant_id(variant);
            
            // For now always work on sample 0. TODO: let the user specify a
            // name and find it.
            int sample_number = 0;

            // What sample is it?
            string& sample_name = variant_file.sampleNames[sample_number];

            // Parse it out and see if it's phased.
            string genotype = variant.getGenotype(sample_name);

            // Tokenize into allele numbers
            // The token iterator can't hold the regex
            regex allele_separator("[|/]");
            for (sregex_token_iterator it(genotype.begin(), genotype.end(), allele_separator, -1);
                it != sregex_token_iterator(); ++it) {
                // For every token separated by / or |
                int allele_number;
                if(it->str() == ".") {
                    // Unknown; pretend it's ref for the purposes of making a
                    // sample graph.
                    allele_number = 0;
                } else {
                    // Parse the allele number
                    allele_number = stoi(it->str());
                }

                // Make the name for its alt path
                string alt_path_name = "_alt_" + var_name + "_" + to_string(allele_number);
                if (graph->has_path(alt_path_name)) {
                    // This alt path is existent and may be nonempty.
                    graph->for_each_step_in_path(graph->get_path_handle(alt_path_name), [&](const step_handle_t& s) {
                        // Un-mark all nodes that are on this alt path, since it is used by the sample.
                        alt_path_ids.erase(graph->get_id(graph->get_handle_of_step(s)));
                    });
                }
            }

        };


        // Allocate a place to store actual variants
        vcflib::Variant var(variant_file);

        while (variant_file.is_open() && variant_file.getNextVariant(var)) {
            // this ... maybe we should remove it as for when we have calls against N
            bool isDNA = allATGC(var.ref);
            for (vector<string>::iterator a = var.alt.begin(); a != var.alt.end(); ++a) {
                if (!allATGC(*a)) isDNA = false;
            }
            // only work with DNA sequences
            if (!isDNA) {
                continue;
            }

            // Handle the variant
            handle_variant(var);
        }


        for(auto& node_id : alt_path_ids) {
            // And delete all the nodes that were used by alt paths that weren't
            // in the genotype of the first sample.
            
            // TODO: keep handles instead maybe?
            handle_t node = graph->get_handle(node_id);
            
            // We need to estroy all paths that touch this node.
            // But we can't do it while iterating (that's asking a lot of the handle graph implementation).
            // So we collect path handles and then destroy them.
            // They need to be a set so we can deduplicate multiple visits of a path.
            unordered_set<path_handle_t> paths_to_remove;
            graph->for_each_step_on_handle(node, [&](const step_handle_t& s) {
                // For every path that touches the node we're destroying,
                // destroy the path. We can't leave it because it won't be the
                // same path without this node.
                paths_to_remove.emplace(graph->get_path_handle_of_step(s));
#ifdef debug
                cerr << "Node " << node_id << " was on path " << graph->get_path_name(graph->get_path_handle_of_step(s)) << endl;
#endif
            });

            for(auto& path : paths_to_remove) {
                graph->destroy_path(path);
            }

            // Actually get rid of the node once its paths are gone.
            graph->destroy_handle(node);
        }

    }
    
    if (!loci_filename.empty()) {
        // Open the file
        ifstream loci_file(loci_filename);
        assert(loci_file.is_open());
    
        // What nodes and edges are called as present by the loci?
        set<handle_t> called_nodes;
        set<edge_t> called_edges;
    
        function<void(Locus&)> lambda = [&](Locus& locus) {
            // For each locus
            
            if (locus.genotype_size() == 0) {
                // No call made here. Just remove all the nodes/edges. TODO:
                // should we keep them all if we don't know if they're there or
                // not? Or should the caller call ref with some low confidence?
                return;
            }
            
            const Genotype& gt = locus.genotype(0);
            
            for (size_t j = 0; j < gt.allele_size(); j++) {
                // For every allele called as present
                int allele_number = gt.allele(j);
                const Path& allele = locus.allele(allele_number);
                
                for (size_t i = 0; i < allele.mapping_size(); i++) {
                    // For every Mapping in the allele
                    const Mapping& m = allele.mapping(i);
                    
                    // Remember to keep this node
                    called_nodes.insert(graph->get_handle(m.position().node_id()));
                    
                    if (i + 1 < allele.mapping_size()) {
                        // Look at the next mapping, which exists
                        const Mapping& m2 = allele.mapping(i + 1);
                        
                        // Find the edge from the last Mapping's node to this one and mark it as used
                        called_edges.insert(graph->edge_handle(graph->get_handle(m.position().node_id(), m.position().is_reverse()),
                            graph->get_handle(m2.position().node_id(), m2.position().is_reverse())));
                    }
                }
            }
        };
        vg::io::for_each(loci_file, lambda);
        
        // Collect all the unused nodes and edges (so we don't try to delete
        // while iterating...)
        set<handle_t> unused_nodes;
        set<edge_t> unused_edges;
        
        graph->for_each_handle([&](const handle_t& n) {
            if (!called_nodes.count(n)) {
                unused_nodes.insert(n);
            }
        });
        
        graph->for_each_edge([&](const edge_t& e) {
            if (!called_edges.count(e)) {
                unused_edges.insert(e);
            }
        });
        
        // Destroy all the extra edges (in case they use extra nodes)
        for (auto& e : unused_edges) {
            graph->destroy_edge(e);
        }
        
        for (auto& n : unused_nodes) {
            graph->destroy_handle(n);
        }
    }
    
    // Some stuff below here needs a vg graph.
    VG* vg_graph = dynamic_cast<vg::VG*>(graph.get());
    
    // Call this to populate the vg_graph if it isn't populated.
    auto ensure_vg = [&]() -> vg::VG* {
        if (vg_graph == nullptr) {
            // Copy instead.
            vg_graph = new vg::VG();
            handlealgs::copy_path_handle_graph(graph.get(), vg_graph);
            // Give the unique_ptr ownership and delete the graph we loaded.
            graph.reset(vg_graph);
            // Make sure the paths are all synced up
            vg_graph->paths.to_graph(vg_graph->graph);
        }
        return vg_graph;
    };
    
    if (!path_name.empty()) {
        // TODO: turn into an algorithm or reimplement
        ensure_vg();
        vg_graph->keep_path(path_name);
    }

    if (unchop) {
        handlealgs::unchop(*graph);
    }

    if (simplify_graph) {
        // Run at up to twice to try and get both ends of nodes.
        // This could be a loop until everything that can simplify does.
        algorithms::simplify_siblings(graph.get()) && algorithms::simplify_siblings(graph.get());
    }

    // check if a handle is contained within a path whose name has nomerge_prefix
    function<bool(const handle_t&)> check_prefix = [&nomerge_prefix, &graph](const handle_t& handle) {
        bool has_prefix = false;
        graph->for_each_step_on_handle(handle, [&nomerge_prefix, &graph, &has_prefix](const step_handle_t& step_handle) {
                string path_name = graph->get_path_name(graph->get_path_handle_of_step(step_handle));
                if (path_name.compare(0, nomerge_prefix.length(), nomerge_prefix) == 0) {
                    has_prefix = true;
                }
                return !has_prefix;
            });
        return has_prefix;
    };
    function<bool(const handle_t&, const handle_t&)> can_merge = nullptr;
    if (!nomerge_prefix.empty()) {        
        can_merge = [&nomerge_prefix, &graph, &check_prefix](const handle_t& h1, const handle_t& h2) {
            return !check_prefix(h1) || !check_prefix(h2);
        };
    }
    
    if (normalize_graph) {
        algorithms::normalize(graph.get(), 1, false, can_merge);
    }

    if (until_normal_iter) {
        // TODO: This doesn't work with vg::VG due to its paths needing re-syncing
        assert(vg_graph == nullptr);
        algorithms::normalize(graph.get(), until_normal_iter, true, can_merge);
    }

    if (remove_non_path) {
        // TODO: turn into an algorithm
        ensure_vg()->remove_non_path();
    }
    
    if (remove_path) {
        // TODO: turn into an algorithm
        ensure_vg()->remove_path();
    }

    if (orient_forward) {
        handlealgs::apply_orientations(graph.get(), handlealgs::topological_order(graph.get()));
    }

    if (flip_doubly_reversed_edges) {
        // TODO: turn into an algorithm
        ensure_vg()->flip_doubly_reversed_edges();
    }

    if (dagify_steps) {
        unordered_map<nid_t, pair<nid_t, bool> > node_translation;
        // TODO: turn into an algorithm
        ensure_vg();
        *vg_graph = vg_graph->dagify(dagify_steps, node_translation, 0, dagify_component_length_max);
    }

    if (dagify_to) {
        unordered_map<nid_t, pair<nid_t, bool> > node_translation;
        // use the walk as our maximum number of steps; it's the worst case
        // TODO: turn into an algorithm
        ensure_vg();
        *vg_graph = vg_graph->dagify(dagify_to, node_translation, dagify_to, dagify_component_length_max);
    }

    if (unfold_to) {
        unordered_map<nid_t, pair<nid_t, bool> > node_translation;
        // TODO: turn into an algorithm
        ensure_vg();
        *vg_graph = vg_graph->unfold(unfold_to, node_translation);
    }

    if (remove_null) {
        // TODO: turn into an algorithm
        ensure_vg()->remove_null_nodes_forwarding_edges();
    }

    if (break_cycles) {
        // TODO: turn into an algorithm
        ensure_vg()->break_cycles();
    }

    // to subset the graph
    if (!root_nodes.empty()) {
        VG g;
        // TODO: turn into an algorithm
        ensure_vg();
        for (auto root : root_nodes) {
            vg_graph->nonoverlapping_node_context_without_paths(vg_graph->get_node(root), g);
            vg_graph->expand_context(g, max(context_steps, 1));
            g.remove_orphan_edges();
        }
        *vg_graph = g;
    }

    // and optionally compact ids
    if (compact_ids) {
        // Sort and compact IDs.
        // TODO: This differs from vg ids! Make an alforithm.
        graph->apply_ordering(handlealgs::topological_order(graph.get()), true);
    }

    if (prune_complex) {
        if (!(path_length > 0 && edge_max > 0)) {
            cerr << "[vg mod]: when pruning complex regions you must specify a --length and --edge-max" << endl;
            return 1;
        }
        algorithms::prune_complex_with_head_tail(*graph, path_length, edge_max);
    }

    if (max_degree) {
        algorithms::remove_high_degree_nodes(*graph, max_degree);
    }

    if (prune_subgraphs) {
        algorithms::prune_short_subgraphs(*graph, path_length);
    }

    if (chop_to) {
        MutablePathDeletableHandleGraph* chop_graph = graph.get();
        if (vg_graph != nullptr) {
            chop_graph = vg_graph;
        }
        
        handlealgs::chop(*chop_graph, chop_to);
        
        if (chop_graph == vg_graph) {
            vg_graph->paths.compact_ranks();
        }
    }

    if (add_start_and_end_markers) {
        if (!(path_length > 0)) {
            cerr << "[vg mod]: when adding start and end markers you must provide a --length" << endl;
            return 1;
        }
        // TODO: replace this with the SourceSinkOverlay, accounting somehow for its immutability.
        Node* head_node = NULL;
        Node* tail_node = NULL;
        vg::id_t head_id = 0, tail_id = 0;
        ensure_vg()->add_start_end_markers(path_length, '#', '$', head_node, tail_node, head_id, tail_id);
    }

    if (destroy_node_id > 0) {
        graph->destroy_handle(graph->get_handle(destroy_node_id));
    }

    if (cactus) {
        // TODO: turn into an algorithm
        ensure_vg();
        // ensure we're sorted
        vg_graph->sort();
        *vg_graph = cactusify(*vg_graph);
        // no paths survive, make sure they are erased
        vg_graph->paths = Paths();
    }

    // Save the modified graph
    vg::io::save_handle_graph(graph.get(), std::cout);

    return 0;
}

// Register subcommand
static Subcommand vg_mod("mod", "filter, transform, and edit the graph", TOOLKIT, main_mod);
