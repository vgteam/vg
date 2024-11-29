#include "subcommand.hpp"
#include "../crash.hpp"
#include "../utility.hpp"
#include "../mapper.hpp"
#include <vg/io/stream.hpp>
#include <vg/io/vpkg.hpp>
#include <vg/io/protobuf_emitter.hpp>
#include "../io/save_handle_graph.hpp"
#include <gbwt/gbwt.h>
#include <gcsa/support.h>
#include "../region.hpp"
#include "../stream_index.hpp"
#include "../algorithms/subgraph.hpp"
#include "../algorithms/sorted_id_ranges.hpp"
#include "../algorithms/approx_path_distance.hpp"
#include "../algorithms/extract_connecting_graph.hpp"
#include "../algorithms/walk.hpp"
#include <bdsg/overlays/overlay_helper.hpp>
#include <htslib/hts.h>

#include <unistd.h>
#include <getopt.h>

using namespace vg;
using namespace vg::subcommand;
using namespace vg::io;

void help_find(char** argv) {
    cerr << "usage: " << argv[0] << " find [options] >sub.vg" << endl
         << "options:" << endl
         << "graph features:" << endl
         << "    -x, --xg-name FILE     use this xg index or graph (instead of rocksdb db)" << endl
         << "    -n, --node ID          find node(s), return 1-hop context as graph" << endl
         << "    -N, --node-list FILE   a white space or line delimited list of nodes to collect" << endl
         << "        --mapping FILE     also include nodes that map to the selected node ids" << endl
         << "    -e, --edges-end ID     return edges on end of node with ID" << endl
         << "    -s, --edges-start ID   return edges on start of node with ID" << endl
         << "    -c, --context STEPS    expand the context of the subgraph this many steps" << endl
         << "    -L, --use-length       treat STEPS in -c or M in -r as a length in bases" << endl
         << "    -P, --position-in PATH find the position of the node (specified by -n) in the given path" << endl
         << "    -I, --list-paths       write out the path names in the index" << endl
         << "    -r, --node-range N:M   get nodes from N to M" << endl
         << "    -G, --gam GAM          accumulate the graph touched by the alignments in the GAM" << endl
         << "    --connecting-start POS find the graph connecting from POS (node ID, + or -, node offset) to --connecting-end" << endl
         << "    --connecting-end POS   find the graph connecting to POS (node ID, + or -, node offset) from --connecting-start" << endl
         << "    --connecting-range INT traverse up to INT bases when going from --connecting-start to --connecting-end (default: 100)" << endl
         << "subgraphs by path range:" << endl
         << "    -p, --path TARGET      find the node(s) in the specified path range(s) TARGET=path[:pos1[-pos2]]" << endl
         << "    -R, --path-bed FILE    read our targets from the given BED FILE" << endl
         << "    -E, --path-dag         with -p or -R, gets any node in the partial order from pos1 to pos2, assumes id sorted DAG" << endl
         << "    -W, --save-to PREFIX   instead of writing target subgraphs to stdout," << endl
         << "                           write one per given target to a separate file named PREFIX[path]:[start]-[end].vg" << endl
         << "    -K, --subgraph-k K     instead of graphs, write kmers from the subgraphs" << endl
         << "    -H, --gbwt FILE        when enumerating kmers from subgraphs, determine their frequencies in this GBWT haplotype index" << endl
         << "alignments:" << endl
         << "    -l, --sorted-gam FILE  use this sorted, indexed GAM file" << endl
         << "    -F, --sorted-gaf FILE  use this sorted, indexed GAF file" << endl
         << "    -o, --alns-on N:M      write alignments which align to any of the nodes between N and M (inclusive)" << endl
         << "    -A, --to-graph VG      get alignments to the provided subgraph" << endl
         << "sequences:" << endl
         << "    -g, --gcsa FILE        use this GCSA2 index of the sequence space of the graph (required for sequence queries)" << endl
         << "    -S, --sequence STR     search for sequence STR using" << endl
         << "    -M, --mems STR         describe the super-maximal exact matches of the STR (gcsa2) in JSON" << endl
         << "    -B, --reseed-length N  find non-super-maximal MEMs inside SMEMs of length at least N" << endl
         << "    -f, --fast-reseed      use fast SMEM reseeding algorithm" << endl
         << "    -Y, --max-mem N        the maximum length of the MEM (default: GCSA2 order)" << endl
         << "    -Z, --min-mem N        the minimum length of the MEM (default: 1)" << endl
         << "    -D, --distance         return distance on path between pair of nodes (-n). if -P not used, best path chosen heurstically" << endl
         << "    -Q, --paths-named S    return all paths whose names are prefixed with S (multiple allowed)" << endl;

}

int main_find(int argc, char** argv) {

    if (argc == 2) {
        help_find(argv);
        return 1;
    }

    string sequence;
    vector<string> kmers;
    vector<vg::id_t> node_ids;
    string node_list_file, node_mapping_file;
    int context_size=0;
    bool use_length = false;
    bool kmer_table = false;
    vector<string> targets_str;
    vector<Region> targets;
    string path_name;
    bool position_in = false;
    string range;
    string gcsa_in;
    string xg_name;
    bool get_mems = false;
    int mem_reseed_length = 0;
    bool use_fast_reseed = true;
    string sorted_gam_name;
    string sorted_gaf_name;
    bool get_mappings = false;
    string aln_on_id_range;
    vg::id_t start_id = 0;
    vg::id_t end_id = 0;
    bool pairwise_distance = false;
    string gam_file;
    pos_t connecting_start = make_pos_t(0, false, 0);
    pos_t connecting_end = make_pos_t(0, false, 0);
    size_t connecting_range = 100;
    int max_mem_length = 0;
    int min_mem_length = 1;
    string to_graph_file;
    bool extract_threads = false;
    vector<string> extract_thread_patterns;
    bool extract_paths = false;
    vector<string> extract_path_patterns;
    bool list_path_names = false;
    bool path_dag = false;
    string bed_targets_file;
    string save_to_prefix;
    int subgraph_k = 0;
    string gbwt_name;

    constexpr int OPT_MAPPING = 1000;
    constexpr int OPT_CONNECTING_START = 1001;
    constexpr int OPT_CONNECTING_END = 1002;
    constexpr int OPT_CONNECTING_RANGE = 1003;

    int c;
    optind = 2; // force optind past command positional argument
    while (true) {
        static struct option long_options[] =
            {
                //{"verbose", no_argument,       &verbose_flag, 1},
                {"xg-name", required_argument, 0, 'x'},
                {"gcsa", required_argument, 0, 'g'},
                {"node", required_argument, 0, 'n'},
                {"node-list", required_argument, 0, 'N'},
                {"mapping", required_argument, 0, OPT_MAPPING},
                {"edges-end", required_argument, 0, 'e'},
                {"edges-start", required_argument, 0, 's'},
                {"sequence", required_argument, 0, 'S'},
                {"mems", required_argument, 0, 'M'},
                {"reseed-length", required_argument, 0, 'B'},
                {"fast-reseed", no_argument, 0, 'f'},
                {"context", required_argument, 0, 'c'},
                {"use-length", no_argument, 0, 'L'},
                {"path", required_argument, 0, 'p'},
                {"path-bed", required_argument, 0, 'R'},
                {"path-dag", no_argument, 0, 'E'},
                {"save-to", required_argument, 0, 'W'},
                {"position-in", required_argument, 0, 'P'},
                {"node-range", required_argument, 0, 'r'},
                {"sorted-gam", required_argument, 0, 'l'},
                {"sorted-gaf", required_argument, 0, 'F'},
                {"mappings", no_argument, 0, 'm'},
                {"alns-on", required_argument, 0, 'o'},
                {"distance", no_argument, 0, 'D'},
                {"gam", required_argument, 0, 'G'},
                {"connecting-start", required_argument, 0, OPT_CONNECTING_START},
                {"connecting-end", required_argument, 0, OPT_CONNECTING_END},
                {"connecting-range", required_argument, 0, OPT_CONNECTING_RANGE},
                {"to-graph", required_argument, 0, 'A'},
                {"max-mem", required_argument, 0, 'Y'},
                {"min-mem", required_argument, 0, 'Z'},
                {"paths-named", required_argument, 0, 'Q'},
                {"list-paths", no_argument, 0, 'I'},
                {"subgraph-k", required_argument, 0, 'K'},
                {"gbwt", required_argument, 0, 'H'},
                {0, 0, 0, 0}
            };

        int option_index = 0;
        c = getopt_long (argc, argv, "x:n:e:s:o:hc:LS:p:P:r:l:F:mg:M:B:fDG:N:A:Y:Z:IQ:ER:W:K:H:",
                         long_options, &option_index);

        // Detect the end of the options.
        if (c == -1)
            break;

        switch (c)
        {

        case 'x':
            xg_name = optarg;
            break;

        case 'g':
            gcsa_in = optarg;
            break;

        case 'S':
            sequence = optarg;
            break;

        case 'M':
            sequence = optarg;
            get_mems = true;
            break;
            
        case 'B':
            mem_reseed_length = parse<int>(optarg);
            break;
            
        case 'f':
            use_fast_reseed = true;
            break;

        case 'Y':
            max_mem_length = parse<int>(optarg);
            break;
            
        case 'Z':
            min_mem_length = parse<int>(optarg);
            break;

        case 'p':
            targets_str.push_back(optarg);
            break;

        case 'R':
            bed_targets_file = optarg;
            break;

        case 'E':
            path_dag = true;
            break;

        case 'P':
            path_name = optarg;
            position_in = true;
            break;

        case 'W':
            save_to_prefix = optarg;
            break;

        case 'c':
            context_size = parse<int>(optarg);
            break;

        case 'L':
            use_length = true;
            break;

        case 'n':
            node_ids.push_back(parse<int>(optarg));
            break;

        case 'N':
            node_list_file = optarg;
            break;

        case OPT_MAPPING:
            node_mapping_file = optarg;
            break;

        case 'e':
            end_id = parse<int>(optarg);
            break;

        case 's':
            start_id = parse<int>(optarg);
            break;

        case 'r':
            range = optarg;
            break;
        
        case 'l':
            sorted_gam_name = optarg;
            break;

        case 'F':
            sorted_gaf_name = optarg;
            break;

        case 'I':
            list_path_names = true;
            break;

        case 'm':
            get_mappings = true;
            break;

        case 'o':
            aln_on_id_range = optarg;
            break;

        case 'D':
            pairwise_distance = true;
            break;

        case 'Q':
            extract_paths = true;
            extract_path_patterns.push_back(optarg);
            break;

        case 'G':
            gam_file = optarg;
            break;
            
        case OPT_CONNECTING_START:
            connecting_start = parse<pos_t>(optarg);
            break;
            
        case OPT_CONNECTING_END:
            connecting_end = parse<pos_t>(optarg);
            break;
        
        case OPT_CONNECTING_RANGE:
            connecting_range = parse<size_t>(optarg);
            break;

        case 'A':
            to_graph_file = optarg;
            break;

        case 'K':
            subgraph_k = atoi(optarg);
            break;

        case 'H':
            gbwt_name = optarg;
            break;

        case 'h':
        case '?':
            help_find(argv);
            exit(1);
            break;

        default:
            abort ();
        }
    }
    if (optind < argc) {
        cerr << "[vg find] find does not accept positional arguments" << endl;
        return 1;
    }

    if (gcsa_in.empty() && xg_name.empty() && sorted_gam_name.empty() && sorted_gaf_name.empty()) {
        cerr << "[vg find] find requires -g, -x, -l, or -F to know where to find its database" << endl;
        return 1;
    }

    if (context_size > 0 && use_length == true && xg_name.empty()) {
        cerr << "[vg find] error, -L not supported without -x" << endl;
        exit(1);
    }
    
    if (xg_name.empty() && mem_reseed_length) {
        cerr << "error:[vg find] SMEM reseeding requires an XG index. Provide XG index with -x." << endl;
        exit(1);
    }
    
    if ((id(connecting_start) == 0) != (id(connecting_end) == 0)) {
        cerr << "error:[vg find] --connecting-start and --connecting-end must be specified together." << endl;
        exit(1);
    }
    
    // process input node list
    if (!node_list_file.empty()) {
        ifstream nli;
        nli.open(node_list_file);
        if (!nli.good()){
            cerr << "[vg find] error, unable to open the node list input file." << endl;
            exit(1);
        }
        string line;
        while (getline(nli, line)){
            for (auto& idstr : split_delims(line, " \t")) {
                node_ids.push_back(parse<nid_t>(idstr.c_str()));
            }
        }
        nli.close();
    }

    // Add the duplicate nodes that map to the original node ids according to the
    // provided node mapping.
    if (!node_mapping_file.empty() && !node_ids.empty()) {
        gcsa::NodeMapping mapping;
        sdsl::load_from_file(mapping, node_mapping_file);
        std::unordered_set<nid_t> original_ids(node_ids.begin(), node_ids.end());
        for (gcsa::size_type id = mapping.begin(); id < mapping.end(); id++) {
            if (original_ids.find(mapping(id)) != original_ids.end()) {
                node_ids.push_back(id);
            }
        }
    }

    PathPositionHandleGraph* xindex = nullptr;
    unique_ptr<PathHandleGraph> path_handle_graph;
    bdsg::PathPositionOverlayHelper overlay_helper;
    bool input_gfa = false;
    if (!xg_name.empty()) {
        path_handle_graph = vg::io::VPKG::load_one<PathHandleGraph>(xg_name);
        input_gfa = dynamic_cast<GFAHandleGraph*>(path_handle_graph.get()) != nullptr;
        xindex = overlay_helper.apply(path_handle_graph.get());

        // Remove node ids that do not exist in the graph.
        std::vector<nid_t> final_ids;
        for (nid_t id : node_ids) {
            if (xindex->has_node(id)) {
                final_ids.push_back(id);
            } else {
                std::cerr << "warning: [vg find] no node with id " << id << " in the graph" << std::endl;
            }
        }
        node_ids = final_ids;
    }
    function<unique_ptr<MutablePathDeletableHandleGraph>()> get_output_graph = [&]() {
        if (input_gfa) {
            return unique_ptr<MutablePathDeletableHandleGraph>(new GFAHandleGraph());
        }
        // todo: move away from VG here
        return unique_ptr<MutablePathDeletableHandleGraph>(new VG());
    };

    unique_ptr<gbwt::GBWT> gbwt_index;
    if (!gbwt_name.empty()) {
        // We are tracing haplotypes, and we want to use the GBWT instead of the old gPBWT.
        // Load the GBWT from its container
        gbwt_index = vg::io::VPKG::load_one<gbwt::GBWT>(gbwt_name.c_str());
        if (gbwt_index.get() == nullptr) {
            // Complain if we couldn't.
            cerr << "error:[vg find] unable to load gbwt index file" << endl;
            return 1;
        }
    }
    
    unique_ptr<GAMIndex> gam_index;
    unique_ptr<vg::io::ProtobufIterator<Alignment>> gam_cursor;
    if (!sorted_gam_name.empty()) {
        // Load the GAM index
        gam_index = unique_ptr<GAMIndex>(new GAMIndex());
        get_input_file(sorted_gam_name + ".gai", [&](istream& in) {
            // We get it form the appropriate .gai, which must exist
            gam_index->load(in); 
        });
    }

    // load GAF index
    tbx_t *gaf_tbx = NULL;
    htsFile *gaf_fp = NULL;
    if (!sorted_gaf_name.empty()){
        gaf_tbx = tbx_index_load3(sorted_gaf_name.c_str(), NULL, 0);
        if ( !gaf_tbx ){
            cerr << "Could not load .tbi/.csi index of " << sorted_gaf_name << endl;
            exit(1);
        }
        int nseq;
        gaf_fp = hts_open(sorted_gaf_name.c_str(),"r");
        if ( !gaf_fp ) {
            cerr << "Could not open " << sorted_gaf_name << endl;
            exit(1);
        }
    }
    
    if (!aln_on_id_range.empty()) {
        // Parse the range
        vector<string> parts = split_delims(aln_on_id_range, ":");
        if (parts.size() == 1) {
            convert(parts.front(), start_id);
            end_id = start_id;
        } else {
            convert(parts.front(), start_id);
            convert(parts.back(), end_id);
        }
        if (gam_index.get() != nullptr) {
            // Find in sorted GAM
            
            get_input_file(sorted_gam_name, [&](istream& in) {
                // Make a cursor for input
                // TODO: Refactor so we can put everything with the GAM index inside one get_input_file call to deduplicate code
                vg::io::ProtobufIterator<Alignment> cursor(in);
                
                // Find the alignments and dump them to cout
                gam_index->find(cursor, start_id, end_id, vg::io::emit_to<Alignment>(cout));
            });
        } else if (!sorted_gaf_name.empty()) {
            // read GAF slice in region 'reg'
            string reg = "{node}:" + convert(start_id) + "-" + convert(end_id);
            hts_itr_t *itr = tbx_itr_querys(gaf_tbx, reg.c_str());
            kstring_t str = {0,0,0};
            if ( itr ) {
                while (tbx_itr_next(gaf_fp, gaf_tbx, itr, &str) >= 0) {
                    puts(str.s);
                }
                tbx_itr_destroy(itr);
            }
        } else {
            cerr << "error [vg find]: Cannot find alignments on range without a sorted GAM or GAF" << endl;
            exit(1);
        }
    }
    
    if (!to_graph_file.empty()) {
        // Find alignments touching a graph
        
        // Load up the graph
        auto graph = vg::io::VPKG::load_one<PathHandleGraph>(to_graph_file);

        if (gam_index.get() != nullptr | !sorted_gaf_name.empty()) {
            // Get the ID ranges from the graph
            auto ranges = vg::algorithms::sorted_id_ranges(graph.get());
            // Throw out the graph
            graph.reset();
            
            if (gam_index.get() != nullptr) {
                // Find in sorted GAM
                get_input_file(sorted_gam_name, [&](istream& in) {
                    // Make a cursor for input
                    // TODO: Refactor so we can put everything with the GAM index inside one get_input_file call to deduplicate code
                    vg::io::ProtobufIterator<Alignment> cursor(in);
                
                    // Find the alignments and send them to cout
                    gam_index->find(cursor, ranges, vg::io::emit_to<Alignment>(cout)); 
                });
              
            } else if (!sorted_gaf_name.empty()) {
                // Find in sorted GAF
                // loop over ranges and print GAF records
                for (auto range : ranges) {
                    string reg = "{node}:" + convert(range.first) + "-" + convert(range.second);
                    hts_itr_t *itr = tbx_itr_querys(gaf_tbx, reg.c_str());
                    kstring_t str = {0,0,0};
                    if ( itr ) {
                        while (tbx_itr_next(gaf_fp, gaf_tbx, itr, &str) >= 0) {
                            puts(str.s);
                        }
                        tbx_itr_destroy(itr);
                    }
                }
            }
        } else {
            cerr << "error [vg find]: Cannot find alignments on graph without a sorted GAM" << endl;
            exit(1);
        }
    }

    if (!xg_name.empty()) {
        if (!node_ids.empty() && path_name.empty() && !pairwise_distance) {
            auto output_graph = get_output_graph();
            auto& graph = *output_graph;
            for (auto node_id : node_ids) {
                graph.create_handle(xindex->get_sequence(xindex->get_handle(node_id)), node_id);
            }
            if (context_size > 0) {
                if (use_length) {
                    vg::algorithms::expand_subgraph_by_length(*xindex, graph, context_size);
                } else {
                    vg::algorithms::expand_subgraph_by_steps(*xindex, graph, context_size);
                }
            } else {
                vg::algorithms::add_connecting_edges_to_subgraph(*xindex, graph);
            }
            vg::algorithms::add_subpaths_to_subgraph(*xindex, graph);

            VG* vg_graph = dynamic_cast<VG*>(&graph);
            if (vg_graph) {
                vg_graph->remove_orphan_edges();
            
                // Order the mappings by rank. TODO: how do we handle breaks between
                // different sections of a path with a single name?
                vg_graph->paths.sort_by_mapping_rank();
            }
            
            // return it
            vg::io::save_handle_graph(&graph, cout);
            // TODO: We're serializing graphs all with their own redundant EOF markers if we use multiple functions simultaneously.
        } else if (end_id != 0) {
            xindex->follow_edges(xindex->get_handle(end_id), false, [&](handle_t next) {
                    edge_t e = xindex->edge_handle(xindex->get_handle(end_id, false), next);
                    cout << (xindex->get_is_reverse(e.first) ? -1 : 1) * xindex->get_id(e.first) << "\t"
                         << (xindex->get_is_reverse(e.second) ? -1 : 1) * xindex->get_id(e.second) << endl;
                });
        } else if (start_id != 0) {
            xindex->follow_edges(xindex->get_handle(start_id), true, [&](handle_t next) {
                    edge_t e = xindex->edge_handle(xindex->get_handle(start_id, true), next);
                    cout << (xindex->get_is_reverse(e.first) ? -1 : 1) * xindex->get_id(e.first) << "\t"
                         << (xindex->get_is_reverse(e.second) ? -1 : 1) * xindex->get_id(e.second) << endl;
                });
        }
        if (!node_ids.empty() && !path_name.empty() && !pairwise_distance && position_in) {
            // Go get the positions of these nodes in this path
            if (xindex->has_path(path_name) == false) {
                // This path doesn't exist, and we'll get a segfault or worse if
                // we go look for positions in it.
                cerr << "[vg find] error, path \"" << path_name << "\" not found in index" << endl;
                exit(1);
            }
            
            // Note: this isn't at all consistent with -P option with rocksdb, which couts a range
            // and then mapping, but need this info right now for scripts/chunked_call
            path_handle_t path_handle = xindex->get_path_handle(path_name);
            for (auto node_id : node_ids) {
                cout << node_id;
                assert(position_in);
                xindex->for_each_step_on_handle(xindex->get_handle(node_id), [&](step_handle_t step_handle) {
                        if (xindex->get_path_handle_of_step(step_handle) == path_handle) {
                            cout << "\t" << xindex->get_position_of_step(step_handle);
                        }
                    });
                cout << endl;
            }
        }
        if (pairwise_distance) {
            if (node_ids.size() != 2) {
                cerr << "[vg find] error, exactly 2 nodes (-n) required with -D" << endl;
                exit(1);
            }
            cout << vg::algorithms::min_approx_path_distance(dynamic_cast<PathPositionHandleGraph*>(&*xindex), make_pos_t(node_ids[0], false, 0), make_pos_t(node_ids[1], false, 0), 1000) << endl;
            return 0;
        }
        if (list_path_names) {
            xindex->for_each_path_handle([&](path_handle_t path_handle) {
                    cout << xindex->get_path_name(path_handle) << endl;
                });
        }
        // handle targets from BED
        if (!bed_targets_file.empty()) {
            parse_bed_regions(bed_targets_file, targets);
        }
        // those given on the command line
        for (auto& target : targets_str) {
            Region region;
            parse_region(target, region);
            targets.push_back(region);
        }
        if (!targets.empty()) {
            auto output_graph = get_output_graph();
            auto& graph = *output_graph;
            auto prep_graph = [&](void) {
                if (context_size > 0) {
                    if (use_length) {
                        vg::algorithms::expand_subgraph_by_length(*xindex, graph, context_size);
                    } else {
                        vg::algorithms::expand_subgraph_by_steps(*xindex, graph, context_size);
                    }
                } else {
                    vg::algorithms::add_connecting_edges_to_subgraph(*xindex, graph);
                }
                vg::algorithms::add_subpaths_to_subgraph(*xindex, graph);
                VG* vg_graph = dynamic_cast<VG*>(&graph);
                if (vg_graph) {
                    vg_graph->remove_orphan_edges();
                    
                    // Order the mappings by rank. TODO: how do we handle breaks between
                    // different sections of a path with a single name?
                    vg_graph->paths.sort_by_mapping_rank();
                }
            };
            for (auto& target : targets) {
                // Grab each target region
                if(xindex->has_path(target.seq) == false) { 
                    // Passing a nonexistent path to get_path_range produces Undefined Behavior
                    cerr << "[vg find] error, path " << target.seq << " not found in index" << endl;
                    exit(1);
                }
                path_handle_t path_handle = xindex->get_path_handle(target.seq);
                // no coordinates given, we do whole thing (0,-1)
                if (target.start < 0 && target.end < 0) {
                    target.start = 0;
                }
                vg::algorithms::extract_path_range(*xindex, path_handle, target.start, target.end, graph);
                if (path_dag) {
                    // find the start and end node of this
                    // and fill things in
                    nid_t id_start = std::numeric_limits<nid_t>::max();
                    nid_t id_end = 1;
                    graph.for_each_handle([&](handle_t handle) {
                            nid_t id = graph.get_id(handle);
                            id_start = std::min(id_start, id);
                            id_end = std::max(id_end, id);
                        });
                    vg::algorithms::extract_id_range(*xindex, id_start, id_end, graph);
                }
                if (!save_to_prefix.empty()) {
                    prep_graph();
                    // write to our save_to file
                    stringstream s;
                    s << save_to_prefix << target.seq;
                    if (target.end >= 0) s << ":" << target.start << ":" << target.end;
                    s << ".vg";
                    ofstream out(s.str().c_str());
                    vg::io::save_handle_graph(&graph, out);
                    out.close();
                    // reset our graph so it has no nodes or paths anymore
                    graph.clear();
                    crash_unless(graph.get_path_count() == 0);
                }
                if (subgraph_k) {
                    prep_graph(); // don't forget to prep the graph, or the kmer set will be wrong[
                    // enumerate the kmers, calculating including their start positions relative to the reference
                    // and write to stdout?
                    bool use_gbwt = false;
                    if (!gbwt_name.empty()) {
                        use_gbwt = true;
                    }
                    vg::algorithms::for_each_walk(
                        graph, subgraph_k, 0,
                        [&](const vg::algorithms::walk_t& walk) {
                            // get the reference-relative position
                            string start_str, end_str;
                            for (auto& p : vg::algorithms::nearest_offsets_in_paths(xindex, walk.begin, subgraph_k*2)) {
                                const uint64_t& start_p = p.second.front().first;
                                const bool& start_rev = p.second.front().second;
                                if (p.first == path_handle && (!start_rev && start_p >= target.start || start_rev && start_p <= target.end)) {
                                    start_str = target.seq + ":" + std::to_string(start_p) + (p.second.front().second ? "-" : "+");
                                }
                            }
                            for (auto& p : vg::algorithms::nearest_offsets_in_paths(xindex, walk.end, subgraph_k*2)) {
                                const uint64_t& end_p = p.second.front().first;
                                const bool& end_rev = p.second.front().second;
                                if (p.first == path_handle && (!end_rev && end_p <= target.end || end_rev && end_p >= target.start)) {
                                    end_str = target.seq + ":" + std::to_string(end_p) + (p.second.front().second ? "-" : "+");
                                }
                            }
                            if (!start_str.empty() && !end_str.empty()) {
                                stringstream ss;
                                ss << target.seq << ":" << target.start << "-" << target.end << "\t"
                                   << walk.seq << "\t" << start_str << "\t" << end_str << "\t";
                                uint64_t on_path = 0;
                                for (auto& h : walk.path) {
                                    xindex->for_each_step_on_handle(xindex->get_handle(graph.get_id(h), graph.get_is_reverse(h)),
                                                                    [&](const step_handle_t& step) {
                                                                        if (xindex->get_path_handle_of_step(step) == path_handle) {
                                                                            ++on_path;
                                                                        }
                                                                    });
                                }
                                // get haplotype frequency
                                if (use_gbwt) {
                                    ss << walk_haplotype_frequency(graph, *gbwt_index, walk) << "\t";
                                } else {
                                    ss << 0 << "\t";
                                }
                                if (on_path == walk.path.size()) {
                                    ss << "ref" << "\t";
                                } else {
                                    ss << "non.ref" << "\t";
                                }
                                for (auto& h : walk.path) {
                                    ss << graph.get_id(h) << (graph.get_is_reverse(h)?"-":"+") << ",";
                                }
                                if (use_gbwt) {
                                    ss << "\t";
                                    for (auto& name : walk_haplotype_names(graph, *gbwt_index, walk)) {
                                        ss << name << ",";
                                    }
                                }
                                // write our record
#pragma omp critical (cout)
                                cout << ss.str() << std::endl;
                            }
                        });
                }
            }
            if (save_to_prefix.empty() && !subgraph_k) {
                prep_graph();
                vg::io::save_handle_graph(&graph, cout);
            }
        }
        if (!range.empty()) {
            auto output_graph = get_output_graph();
            auto& graph = *output_graph;
            nid_t id_start=0, id_end=0;
            vector<string> parts = split_delims(range, ":");
            if (parts.size() == 1) {
                cerr << "[vg find] error, format of range must be \"N:M\" where start id is N and end id is M, got " << range << endl;
                exit(1);
            }
            convert(parts.front(), id_start);
            convert(parts.back(), id_end);
            if (!use_length) {
                vg::algorithms::extract_id_range(*xindex, id_start, id_end, graph);
            } else {
                // treat id_end as length instead.
                size_t length = 0;
                nid_t found_id_end = id_start;
                for (nid_t cur_id = id_start; length < id_end; ++cur_id) {
                    if (!xindex->has_node(cur_id)) {
                        break;
                    }
                    length += xindex->get_length(xindex->get_handle(cur_id));
                    found_id_end = cur_id;
                }
                vg::algorithms::extract_id_range(*xindex, id_start, found_id_end, graph);
            }
            if (context_size > 0) {
                if (use_length) {
                    vg::algorithms::expand_subgraph_by_length(*xindex, graph, context_size);
                } else {
                    vg::algorithms::expand_subgraph_by_steps(*xindex, graph, context_size);
                }
            } else {
                vg::algorithms::add_connecting_edges_to_subgraph(*xindex, graph);
            }
            vg::algorithms::add_subpaths_to_subgraph(*xindex, graph);

            VG* vg_graph = dynamic_cast<VG*>(&graph);
            if (vg_graph) {
                vg_graph->remove_orphan_edges();
            }
            vg::io::save_handle_graph(&graph, cout);
        }
        if (extract_paths) {
            for (auto& pattern : extract_path_patterns) {
            
                // We want to write uncompressed protobuf Graph objects containing our paths.
                vg::io::ProtobufEmitter<Graph> out(cout, false);
            
                xindex->for_each_path_handle([&](path_handle_t path_handle) {
                        string path_name = xindex->get_path_name(path_handle);
                        if (pattern.length() <= path_name.length() && path_name.compare(0, pattern.length(), pattern) == 0) {
                            // We need a Graph for serialization purposes.
                            Graph g;
                            *g.add_path() = path_from_path_handle(*xindex, path_handle);
                            out.write(std::move(g));
                        }
                    });
            }
        }
        if (!gam_file.empty()) {
            set<vg::id_t> nodes;
            function<void(Alignment&)> lambda = [&nodes](Alignment& aln) {
                // accumulate nodes matched by the path
                auto& path = aln.path();
                if (path.mapping_size() == 1 && !path.mapping(0).has_position() &&
                    path.mapping(0).edit_size() == 1 && edit_is_insertion(path.mapping(0).edit(0))) {
                    // This read is a (presumably full length) insert to no
                    // position. The aligner (used to?) generate these for some
                    // unmapped reads. We should skip it.
                    return;
                }
                for (int i = 0; i < path.mapping_size(); ++i) {
                    nodes.insert(path.mapping(i).position().node_id());
                }
            };
            if (gam_file == "-") {
                vg::io::for_each(std::cin, lambda);
            } else {
                ifstream in;
                in.open(gam_file.c_str());
                if(!in.is_open()) {
                    cerr << "[vg find] error: could not open alignments file " << gam_file << endl;
                    exit(1);
                }
                vg::io::for_each(in, lambda);
            }
            // now we have the nodes to get
            auto output_graph = get_output_graph();
            auto& graph = *output_graph;
            for (auto& node : nodes) {
                handle_t node_handle = xindex->get_handle(node);
                graph.create_handle(xindex->get_sequence(node_handle), xindex->get_id(node_handle));
            }
            vg::algorithms::expand_subgraph_by_steps(*xindex, graph, max(1, context_size)); // get connected edges
            vg::algorithms::add_connecting_edges_to_subgraph(*xindex, graph);
            vg::io::save_handle_graph(&graph, cout);
        }
        if (id(connecting_start) != 0) {
            // Extract a connecting graph
            auto output_graph = get_output_graph();
            auto& graph = *output_graph;
            vg::algorithms::extract_connecting_graph(xindex, &graph, connecting_range, connecting_start, connecting_end);
            vg::io::save_handle_graph(&graph, cout);
        }
    }

    // todo cleanup if/else logic to allow only one function

    if (!sequence.empty()) {
        if (gcsa_in.empty()) {
            cerr << "error:[vg find] need GCSA index to query sequences" << endl;
            return 1;
        }
        
        // Configure GCSA2 verbosity so it doesn't spit out loads of extra info
        gcsa::Verbosity::set(gcsa::Verbosity::SILENT);
        
        // Open it
        auto gcsa_index = vg::io::VPKG::load_one<gcsa::GCSA>(gcsa_in);
        // default LCP is the gcsa base name +.lcp
        auto lcp_index = vg::io::VPKG::load_one<gcsa::LCPArray>(gcsa_in + ".lcp");
        
        //range_type find(const char* pattern, size_type length) const;
        //void locate(size_type path, std::vector<node_type>& results, bool append = false, bool sort = true) const;
        //locate(i, results);
        if (!get_mems) {
            auto paths = gcsa_index->find(sequence.c_str(), sequence.length());
            //cerr << paths.first << " - " << paths.second << endl;
            for (gcsa::size_type i = paths.first; i <= paths.second; ++i) {
                std::vector<gcsa::node_type> ids;
                gcsa_index->locate(i, ids);
                for (auto id : ids) {
                    cout << gcsa::Node::decode(id) << endl;
                }
            }
        } else {
            // for mems we need to load up the gcsa and lcp structures into the mapper
            Mapper mapper(xindex, gcsa_index.get(), lcp_index.get());
            mapper.fast_reseed = use_fast_reseed;
            // get the mems
            double lcp_avg, fraction_filtered;
            auto mems = mapper.find_mems_deep(sequence.begin(), sequence.end(), lcp_avg, fraction_filtered, max_mem_length, min_mem_length, mem_reseed_length);
            
            // dump them to stdout
            cout << mems_to_json(mems) << endl;
            
        }
    }
    
    return 0;

}

static Subcommand vg_msga("find", "use an index to find nodes, edges, kmers, paths, or positions", DEVELOPMENT, main_find);
