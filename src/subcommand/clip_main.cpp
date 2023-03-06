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
       << "    -b, --bed FILE            BED regions corresponding to path intervals of the graph to target" << endl
       << "    -r, --snarls FILE         Snarls from vg snarls (recomputed if not given unless -d and -P used)." << endl
       << "depth clipping options: " << endl
       << "    -d, --depth N             Clip out nodes and edges with path depth below N" << endl
       << "stub clipping options:" << endl
       << "    -s, --stubs               Clip out all stubs (nodes with degree-0 sides that aren't on reference)" << endl
       << "snarl complexity clipping options: [default mode]" << endl
       << "    -n, --max-nodes N         Only clip out snarls with > N nodes" << endl
       << "    -e, --max-edges N         Only clip out snarls with > N edges" << endl
       << "    -N  --max-nodes-shallow N Only clip out snarls with > N nodes not including nested snarls" << endl
       << "    -E  --max-edges-shallow N Only clip out snarls with > N edges not including nested snarls" << endl
       << "    -a, --max-avg-degree N    Only clip out snarls with average degree > N" << endl
       << "    -l, --max-reflen-prop F   Ignore snarls whose reference traversal spans more than F (0<=F<=1) of the whole reference path" << endl
       << "    -L, --max-reflen N        Ignore snarls whose reference traversal spans more than N bp" << endl
       << "big deletion edge clipping options:" << endl
       << "    -D, --max-deletion-edge N Clip out all edges whose endpoints have distance > N on a reference path" << endl
       << "    -c, --context N           Search up to at most N steps from reference paths for candidate deletion edges [1]" << endl
       << "general options: " << endl
       << "    -P, --path-prefix STRING  Do not clip out alleles on paths beginning with given prefix (such references must be specified either with -P or -b). Multiple allowed" << endl
       << "    -m, --min-fragment-len N  Don't write novel path fragment if it is less than N bp long" << endl
       << "    -B, --output-bed          Write BED-style file of affected intervals instead of clipped graph. " << endl
       << "                              Columns 4-9 are: snarl node-count edge-count shallow-node-count shallow-edge-count avg-degree" << endl
       << "    -t, --threads N           number of threads to use [default: all available]" << endl
       << "    -v, --verbose             Print some logging messages" << endl
       << endl;
}    

int main_clip(int argc, char** argv) {

    string bed_path;
    string snarls_path;
    vector<string> ref_prefixes;
    int64_t min_depth = -1;
    int64_t min_fragment_len = 0;
    bool verbose = false;
    bool depth_clipping = false;
    bool stub_clipping = false;

    size_t max_nodes = 0;
    size_t max_edges = 0;
    size_t max_nodes_shallow = 0;
    size_t max_edges_shallow = 0;
    double max_avg_degree = 0.;
    double max_reflen_prop = numeric_limits<double>::max();
    size_t max_reflen = numeric_limits<size_t>::max();
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
            {"max-nodes", required_argument, 0, 'n'},
            {"max-edges", required_argument, 0, 'e'},
            {"max-nodes-shallow", required_argument, 0, 'N'},
            {"max-edges-shallow", required_argument, 0, 'E'},
            {"max-avg-degree", required_argument, 0, 'a'},
            {"max-reflen-prop", required_argument, 0, 'l'},
            {"max-reflen", required_argument, 0, 'L'},
            {"max-deletion", required_argument, 0, 'D'},
            {"context", required_argument, 0, 'c'},
            {"path-prefix", required_argument, 0, 'P'},
            {"snarls", required_argument, 0, 'r'},
            {"min-fragment-len", required_argument, 0, 'm'},
            {"output-bed", no_argument, 0, 'B'},
            {"threads", required_argument, 0, 't'},
            {"verbose", required_argument, 0, 'v'},
            {0, 0, 0, 0}

        };
        int option_index = 0;
        c = getopt_long (argc, argv, "hb:d:sn:e:N:E:a:l:L:D:c:P:r:m:Bt:v",
                long_options, &option_index);

        // Detect the end of the options.
        if (c == -1)
            break;

        switch (c)
        {

        case '?':
        case 'h':
            help_clip(argv);
            return 0;
        case 'b':
            bed_path = optarg;
            break;
        case 'd':
            min_depth = parse<size_t>(optarg);
            break;
        case 's':
            stub_clipping = true;
            break;
        case 'n':
            max_nodes = parse<size_t>(optarg);
            snarl_option = true;
            break;
        case 'e':
            max_edges = parse<size_t>(optarg);
            snarl_option = true;
            break;
        case 'N':
            max_nodes_shallow = parse<size_t>(optarg);
            snarl_option = true;
            break;
        case 'E':
            max_edges_shallow = parse<size_t>(optarg);
            snarl_option = true;
            break;            
        case 'a':
            max_avg_degree = parse<double>(optarg);
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
            snarls_path = optarg;
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
        {
            int num_threads = parse<int>(optarg);
            if (num_threads <= 0) {
                cerr << "error:[vg clip] Thread count (-t) set to " << num_threads << ", must set to a positive integer." << endl;
                exit(1);
            }
            omp_set_num_threads(num_threads);
            break;
        }            
        default:
            abort();
        }
    }

    if (bed_path.empty() == ref_prefixes.empty()) {
        cerr << "error:[vg-clip] Reference intervals must be specified with one of -b or -P" << endl;
        return 1;
    }

    if ((min_depth >= 0 || max_deletion >= 0 || stub_clipping) && (snarl_option || out_bed)) {
        cerr << "error:[vg-clip] bed output (-B) and snarl complexity options (-n, -e, -N, -E, -a, -l, -L) cannot be used with -d, -D or -s" << endl;
        return 1;
    }

    // to do: I think it could be a good idea to combine these options
    if (min_depth >= 0 && max_deletion >= 0) {
        cerr << "error:[vg-clip] -d cannot (yet?) be used with -D" << endl;
        return 1;
    }

    // ditto about combining
    if (stub_clipping && (min_depth >= 0 || max_deletion >= 0)) {
        cerr << "error:[vg-clip] -s cannot (yet?) be used with -d or -D" << endl;
        return 1;
    }
    
    if (context_steps >= 0 && max_deletion < 0) {
        cerr << "error:[vg-clip] -c can only be used with -D" << endl;
        return 1;
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
    bool need_snarls = !bed_path.empty() || (min_depth < 0 && max_deletion < 0 && !stub_clipping);

    if (need_pp) {
        pp_graph = overlay_helper.apply(graph.get());
        if (verbose) {
            cerr << "[vg clip]: Computed path position overlay of input graph" << endl;
        }
    }

    if (need_snarls) {
        // Load or compute the snarls which are required for targetting bed regions
        if (!snarls_path.empty()) {
            ifstream snarl_file(snarls_path.c_str());
            if (!snarl_file) {
                cerr << "Error [vg clip]: Unable to load snarls file: " << snarls_path << endl;
                return 1;
            }
            snarl_manager = vg::io::VPKG::load_one<SnarlManager>(snarl_file);
            if (verbose) {
                cerr << "[vg clip]: Loaded " << snarl_manager->num_snarls() << " snarls" << endl;
            }
        } else {
            IntegratedSnarlFinder finder(*graph);
            snarl_manager = unique_ptr<SnarlManager>(new SnarlManager(std::move(finder.find_snarls_parallel())));
            if (verbose) {
                cerr << "[vg clip]: Computed " << snarl_manager->num_snarls() << " snarls" << endl;
            }
        }
        
        // load the bed file
        if (!bed_path.empty()) {
            parse_bed_regions(bed_path, bed_regions);
            if (verbose) {
                cerr << "[vg clip]: Loaded " << bed_regions.size() << " BED regions" << endl;
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
                    cerr << "[vg clip]: Dropped " << (bed_regions.size() - bed_regions_in_graph.size()) << " BED regions whose sequence names do not correspond to paths in the graph" << endl;
                }
                if (bed_regions_in_graph.empty()) {
                    cerr << "warning:[vg-clip] No BED region found that lies on path in graph (use vg paths -Lv to list paths that are in the graph)" << endl;
                }
            }
            swap(bed_regions, bed_regions_in_graph);
        } else {
            assert(need_pp);
            assert(!ref_prefixes.empty());
            // load the bed regions from the reference path prefix
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
                cerr << "[vg clip]: Inferred " << bed_regions.size() << " BED regions from paths in the graph" << endl;
            }
        }
    }        

    if (min_depth >= 0) {
        // run the depth clipping       
        if (bed_path.empty()) {            
            // do the whole graph
            clip_low_depth_nodes_and_edges(graph.get(), min_depth, ref_prefixes, min_fragment_len, verbose);
        } else {
            // do the contained snarls
            clip_contained_low_depth_nodes_and_edges(graph.get(), pp_graph, bed_regions, *snarl_manager, false, min_depth, min_fragment_len, verbose);
        }
        
    } else if (max_deletion >= 0) {
        // run the deletion edge clipping on the whole graph
        clip_deletion_edges(graph.get(), max_deletion, context_steps, ref_prefixes, min_fragment_len, verbose);
    } else if (stub_clipping) {
        // run the stub clipping
        if (bed_path.empty()) {            
            // do the whole graph
            clip_stubs(graph.get(), ref_prefixes, min_fragment_len, verbose);
        } else {
            // do the contained snarls
            clip_contained_stubs(graph.get(), pp_graph, bed_regions, *snarl_manager, false, min_fragment_len, verbose);
        }        
    }else {
        // run the alt-allele clipping
        clip_contained_snarls(graph.get(), pp_graph, bed_regions, *snarl_manager, false, min_fragment_len,
                              max_nodes, max_edges, max_nodes_shallow, max_edges_shallow, max_avg_degree, max_reflen_prop, max_reflen, out_bed, verbose);
    }

    // write the graph
    if (!out_bed) {
        vg::io::save_handle_graph(graph.get(), std::cout);
    }
    
    return 0;
}


// Register subcommand
static Subcommand vg_clip("clip", "remove BED regions (other other nodes from their snarls) from a graph", main_clip);
