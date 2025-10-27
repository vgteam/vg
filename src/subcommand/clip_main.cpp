#include "subcommand.hpp"
#include "../vg.hpp"
#include "../utility.hpp"
#include "../integrated_snarl_finder.hpp"
#include "../io/save_handle_graph.hpp"
#include <vg/io/stream.hpp>
#include <vg/io/vpkg.hpp>
#include <vg/io/alignment_emitter.hpp>
#include "../clip.hpp"
#include <bdsg/overlays/overlay_helper.hpp>

#include <unistd.h>
#include <getopt.h>

using namespace vg;
using namespace vg::subcommand;
using namespace vg::io;

void help_clip(char** argv) {
    cerr << "usage: " << argv[0] << " [options] <graph>" << endl
         << "Chop out variation within path intervals of a vg graph" << endl
         << endl
         << "input options: " << endl
         << "  -b, --bed FILE             BED regions for path intervals of graph to target" << endl
         << "  -r, --snarls FILE          snarls from vg snarls (recomputed if unspecified" << endl
         << "                             and -n, -e, -N, -E, -a, -l, -L, or -D used)" << endl
         << "depth clipping options: " << endl
         << "  -d, --depth N              clip out nodes and edges with path depth below N" << endl
         << "                             use snarl selection/BED options to target regions" << endl
         << "stub clipping options:" << endl
         << "  -s, --stubs                clip out all stubs, i.e. nodes with degree-0 sides" << endl
         << "                             that aren't on reference" << endl
         << "  -S, --stubbify-paths       clip out all edges neccessary so that selected" << endl
         << "                             reference paths have exactly two stubs" << endl
         << "snarl selection and clipping options:" << endl
         << "  -n, --min-nodes N          clip out snarls with >N nodes" << endl
         << "  -e, --min-edges N          clip out snarls with >N edges" << endl
         << "  -i, --min-bases N          clip out snarls with >N bases" << endl     
         << "  -N, --min-nodes-shallow N  clip out snarls with >N nodes," << endl
         << "                             ignoring nested snarls" << endl
         << "  -E, --min-edges-shallow N  clip out snarls with >N edges," << endl
         << "                             ignoring nested snarls" << endl
         << "  -a, --min-avg-degree N     clip out snarls with average degree > N" << endl
         << "  -A, --min-reflen N         ignore snarls with reference traversal span < N bp" << endl     
         << "  -l, --max-reflen-prop F    ignore snarls with reference traversal span > F" << endl
         << "                             (0<=F<=1) of the whole reference path" << endl
         << "  -L, --max-reflen N         ignore snarls with reference traversal span > N bp" << endl
         << "  -g, --net-edges            only clip net-edges inside snarls" << endl
         << "  -G, --top-net-edges        only clip net-edges inside top-level snarls" << endl
         << "big deletion edge clipping options:" << endl
         << "  -D, --max-deletion-edge N  clip out all edges whose endpoints have" << endl
         << "                             distance > N on a reference path" << endl
         << "  -c, --context N            search up to at most N steps from reference paths" << endl
         << "                             for candidate deletion edges [1]" << endl
         << "general options: " << endl
         << "  -P, --path-prefix STR      do not clip out alleles on paths beginning with STR" << endl
         << "                             (specify such references with -P or -b; may repeat)" << endl
         << "  -m, --min-fragment-len N   don't write novel path fragment if < N bp long" << endl
         << "  -B, --output-bed           write BED-style file of intervals," << endl
         << "                             instead of clipped grap. columns 4-9 are:" << endl
         << "                             snarl node-count edge-count shallow-node-count" << endl
         << "                             shallow-edge-count avg-degree" << endl
         << "  -t, --threads N            number of threads to use [all available]" << endl
         << "  -v, --verbose              print some logging messages" << endl
         << "  -h, --help                 print this help message to stderr and exit" << endl;
}    

int main_clip(int argc, char** argv) {
    Logger logger("vg clip");

    string bed_path;
    string snarls_path;
    vector<string> ref_prefixes;
    int64_t min_depth = -1;
    int64_t min_fragment_len = 0;
    bool verbose = false;
    bool depth_clipping = false;
    bool stub_clipping = false;
    bool stubbify_reference = false;

    size_t min_nodes = 0;
    size_t min_edges = 0;
    size_t min_bases = 0;
    size_t min_nodes_shallow = 0;
    size_t min_edges_shallow = 0;
    double min_avg_degree = 0.;
    size_t min_reflen = 0;
    double max_reflen_prop = 1.;    
    size_t max_reflen = numeric_limits<size_t>::max();
    bool only_net_edges = false;
    bool only_top_net_edges = false;
    bool out_bed = false;
    bool snarl_option = false;
    
    int64_t max_deletion = -1;
    int64_t context_steps = -1;

    if (argc == 2) {
        help_clip(argv);
        return 1;
    }
    
    int c;
    optind = 2; // force optind past command positional argument
    while (true) {
        static struct option long_options[] =
        {
            {"help", no_argument, 0, 'h'},
            {"bed", required_argument, 0, 'b'},
            {"depth", required_argument, 0, 'd'},
            {"stubs", no_argument, 0, 's'},
            {"stubbify-paths", no_argument, 0, 'S'},
            {"min-nodes", required_argument, 0, 'n'},
            {"min-edges", required_argument, 0, 'e'},
            {"min-bases", required_argument, 0, 'i'},
            {"min-nodes-shallow", required_argument, 0, 'N'},
            {"min-edges-shallow", required_argument, 0, 'E'},
            {"min-avg-degree", required_argument, 0, 'a'},
            {"max-reflen-prop", required_argument, 0, 'l'},
            {"max-reflen", required_argument, 0, 'L'},
            {"min-reflen", required_argument, 0, 'A'},
            {"net-edges", no_argument, 0, 'g'},
            {"top-net-edges", no_argument, 0, 'G'},            
            {"max-deletion-edge", required_argument, 0, 'D'},
            {"context", required_argument, 0, 'c'},
            {"path-prefix", required_argument, 0, 'P'},
            {"snarls", required_argument, 0, 'r'},
            {"min-fragment-len", required_argument, 0, 'm'},
            {"output-bed", no_argument, 0, 'B'},
            {"threads", required_argument, 0, 't'},
            {"verbose", no_argument, 0, 'v'},
            {0, 0, 0, 0}

        };
        int option_index = 0;
        c = getopt_long (argc, argv, "h?b:d:sSn:e:i:N:E:a:l:L:A:gGD:c:P:r:m:Bt:v",
                         long_options, &option_index);

        // Detect the end of the options.
        if (c == -1)
            break;

        switch (c)
        {

        case '?':
        case 'h':
            help_clip(argv);
            return 1;
        case 'b':
            bed_path = require_exists(logger, optarg);
            break;
        case 'd':
            min_depth = parse<size_t>(optarg);
            break;
        case 's':
            stub_clipping = true;
            break;
        case 'S':
            stubbify_reference = true;
            break;            
        case 'n':
            min_nodes = parse<size_t>(optarg);
            snarl_option = true;
            break;
        case 'e':
            min_edges = parse<size_t>(optarg);
            snarl_option = true;
            break;
        case 'i':
            min_bases = parse<size_t>(optarg);
            snarl_option = true;
            break;            
        case 'N':
            min_nodes_shallow = parse<size_t>(optarg);
            snarl_option = true;
            break;
        case 'E':
            min_edges_shallow = parse<size_t>(optarg);
            snarl_option = true;
            break;            
        case 'a':
            min_avg_degree = parse<double>(optarg);
            snarl_option = true;
            break;
        case 'l':
            max_reflen_prop = parse<double>(optarg);
            snarl_option = true;
            break;
        case 'L':
            max_reflen = parse<size_t>(optarg);
            snarl_option = true;
            break;
        case 'A':
            min_reflen = parse<size_t>(optarg);
            snarl_option = true;
            break;            
        case 'g':
            only_net_edges = true;
            break;
        case 'G':
            only_top_net_edges = true;
            break;                        
        case 'D':
            max_deletion = parse<size_t>(optarg);
            break;
        case 'c':
            context_steps = parse<size_t>(optarg);
            break;            
        case 'P':
            ref_prefixes.push_back(optarg);
            break;
        case 'r':
            snarls_path = require_exists(logger, optarg);
            break;
        case 'm':
            min_fragment_len = parse<int>(optarg);
            break;
        case 'B':
            out_bed = true;
            break;
        case 'v':
            verbose = true;
            break;
        case 't':
            set_thread_count(logger, optarg);
            break;       
        default:
            abort();
        }
    }

    if (bed_path.empty() == ref_prefixes.empty()) {
        logger.error() << "Reference intervals must be specified with one of -b or -P" << endl;
    }

    if ((max_deletion >= 0 || stub_clipping || stubbify_reference) && (snarl_option || out_bed)) {
        logger.error() << "BED output (-B) and snarl complexity options (-n, -e, -N, -E, -a, -l, -L, -A) "
                       << "cannot be used with -D, -s or -S" << endl;
    }

    // ditto about combining
    if ((stub_clipping || stubbify_reference) && (min_depth >= 0 || max_deletion >= 0)) {
        logger.error() << "-s and -S cannot (yet?) be used with -d or -D" << endl;
    }
    
    if (context_steps >= 0 && max_deletion < 0) {
        logger.error() << "-c can only be used with -D" << endl;
    }

    if (stubbify_reference && ref_prefixes.empty()) {
        logger.error() << "-S can only be used with -P" << endl;
    }

    if (only_net_edges && only_top_net_edges) {
        logger.error() << "-g and -G cannot be used together: choose one" << endl;
    }

    // default to same
    if (max_deletion > 0 && context_steps < 0) {
        context_steps = max_deletion;
    }

    // load the graph
    string graph_path = get_input_file_name(optind, argc, argv);
    unique_ptr<MutablePathMutableHandleGraph> graph = vg::io::VPKG::load_one<MutablePathMutableHandleGraph>(graph_path);

    // optional overlay only needed with bed regions input
    bdsg::PathPositionOverlayHelper overlay_helper;
    PathPositionHandleGraph* pp_graph = nullptr;

    unique_ptr<SnarlManager> snarl_manager;
    vector<Region> bed_regions;

    // need the path positions unless we're doing depth, deletion or stub clipping without regions
    bool need_pp = !(bed_path.empty() && (min_depth >= 0 || max_deletion >= 0 || stub_clipping));

    // need snarls if input regions are provided, or doing snarl based clipping
    bool need_snarls = snarl_option || !bed_path.empty();

    // TodO: FIX!!  shouldn't need pp without BED coordinates
    need_pp = need_pp || need_snarls;

    if (need_pp) {
        pp_graph = overlay_helper.apply(graph.get());
        if (verbose) {
            logger.info() << "Computed path position overlay of input graph" << endl;
        }
    }

    if (need_snarls) {
        // Load or compute the snarls which are required for targetting BED regions
        if (!snarls_path.empty()) {
            ifstream snarl_file(snarls_path.c_str());
            snarl_manager = vg::io::VPKG::load_one<SnarlManager>(snarl_file);
            if (verbose) {
                logger.info() << "Loaded " << snarl_manager->num_snarls() << " snarls" << endl;
            }
        } else {
            IntegratedSnarlFinder finder(*graph);
            snarl_manager = unique_ptr<SnarlManager>(new SnarlManager(std::move(finder.find_snarls_parallel())));
            if (verbose) {
                logger.info() << "Computed " << snarl_manager->num_snarls() << " snarls" << endl;
            }
        }
        
        // load the BED file
        if (!bed_path.empty()) {
            parse_bed_regions(bed_path, bed_regions);
            if (verbose) {
                logger.info() << "Loaded " << bed_regions.size() << " BED regions" << endl;
            }
            // contig names left in this set are *not* in the graph
            unordered_set<string> contig_set;
            for (const Region& region : bed_regions) {
                contig_set.insert(region.seq);
            }
            graph->for_each_path_handle([&] (path_handle_t path_handle) {
                    string base_name = Paths::strip_subrange(graph->get_path_name(path_handle));
                    if (contig_set.count(base_name)) {
                        // todo: should take into account coordinate comp
                        contig_set.erase(base_name);
                    }
                });
            vector<Region> bed_regions_in_graph;
            for (const Region& region : bed_regions) {
                if (!contig_set.count(region.seq)) {
                    bed_regions_in_graph.push_back(region);
                }
            }
            if (bed_regions_in_graph.size() != bed_regions.size()) {
                if (verbose) {
                    logger.info() << "Dropped " << (bed_regions.size() - bed_regions_in_graph.size()) 
                                  << " BED regions whose sequence names "
                                  << "do not correspond to paths in the graph" << endl;
                }
                if (bed_regions_in_graph.empty()) {
                    logger.warn() << "No BED region found that lies on path in graph "
                                  << "(use vg paths -Lv to list paths that are in the graph)" << endl;
                }
            }
            swap(bed_regions, bed_regions_in_graph);
        } else {
            assert(need_pp);
            assert(!ref_prefixes.empty());
            // load the BED regions from the reference path prefix
            pp_graph->for_each_path_handle([&](path_handle_t path_handle) {
                    string path_name = pp_graph->get_path_name(path_handle);
                    subrange_t subrange;
                    path_name = Paths::strip_subrange(path_name, &subrange);
                    int64_t offset = subrange == PathMetadata::NO_SUBRANGE ? 0 : subrange.first;
                    for (const string& ref_prefix : ref_prefixes) {
                        if (path_name.compare(0, ref_prefix.length(), ref_prefix) == 0) {
                            Region region = {path_name, offset, offset + (int64_t)pp_graph->get_path_length(path_handle) - 1};
                            bed_regions.push_back(region);
                            break;
                        }
                    }
                });
            if (verbose) {
                logger.info() << "Inferred " << bed_regions.size() 
                              << " BED regions from paths in the graph" << endl;
            }
        }
    }        
    
    if (min_depth >= 0 && max_deletion < 0) {
        // run the depth clipping       
        if (bed_regions.empty()) {            
            // do the whole graph
            clip_low_depth_nodes_and_edges(graph.get(), min_depth, ref_prefixes, min_fragment_len, verbose);
        } else {
            // do the contained snarls
            clip_contained_low_depth_nodes_and_edges(graph.get(), pp_graph, bed_regions, *snarl_manager, false, 
                                                     min_depth, min_fragment_len, min_nodes, min_edges, min_bases,
                                                     min_nodes_shallow, min_edges_shallow, min_avg_degree, min_reflen,
                                                     max_reflen_prop, max_reflen, only_net_edges || only_top_net_edges,
                                                     only_top_net_edges, out_bed, verbose);
        }
        
    } else if (max_deletion >= 0) {
        // run the deletion edge clipping on the whole graph
        clip_deletion_edges(graph.get(), max_deletion, context_steps, min_depth, ref_prefixes, min_fragment_len, verbose);
    } else if (stub_clipping || stubbify_reference) {
        // run the stub clipping
        if (bed_path.empty()) {            
            // do the whole graph
            if (stubbify_reference) {
                // important that this is done first, as it can actually create non-reference stubs that'd need removal below
                stubbify_ref_paths(graph.get(), ref_prefixes, min_fragment_len, verbose);
            }
            if (stub_clipping) {
                clip_stubs(graph.get(), ref_prefixes, min_fragment_len, verbose);
            }
        } else {
            assert(stub_clipping && !stubbify_reference);
            // do the contained snarls
            clip_contained_stubs(graph.get(), pp_graph, bed_regions, *snarl_manager, false, min_fragment_len, verbose);
        }        
    }else {
        // run the alt-allele clipping
        clip_contained_snarls(graph.get(), pp_graph, bed_regions, *snarl_manager, false, min_fragment_len,
                              min_nodes, min_edges, min_bases, min_nodes_shallow, min_edges_shallow,
                              min_avg_degree, min_reflen, max_reflen_prop, max_reflen,
                              only_net_edges || only_top_net_edges, only_top_net_edges, out_bed, verbose);
    }

    // write the graph
    if (!out_bed) {
        vg::io::save_handle_graph(graph.get(), std::cout);
    }
    
    return 0;
}


// Register subcommand
static Subcommand vg_clip("clip", "remove BED regions (other other nodes from their snarls) from a graph", main_clip);
