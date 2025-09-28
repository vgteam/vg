// index.cpp: define the "vg index" subcommand, which makes XG, GCSA2, and distance indexes

#include <omp.h>
#include <unistd.h>
#include <getopt.h>

#include <random>
#include <string>
#include <vector>

#include "subcommand.hpp"

#include "../vg.hpp"
#include "xg.hpp"
#include <vg/io/stream.hpp>
#include <vg/io/vpkg.hpp>
#include "../io/save_handle_graph.hpp"
#include "../stream_index.hpp"
#include "../vg_set.hpp"
#include "../utility.hpp"
#include "../region.hpp"
#include "../integrated_snarl_finder.hpp"
#include "../snarl_distance_index.hpp"
#include "../source_sink_overlay.hpp"
#include "../gbwtgraph_helper.hpp"
#include "../gcsa_helper.hpp"

#include <gcsa/algorithms.h>
#include <bdsg/overlays/packed_subgraph_overlay.hpp>
#include <handlegraph/algorithms/weakly_connected_components.hpp>

using namespace std;
using namespace vg;
using namespace vg::subcommand;

const string context = "[vg index]";

void help_index(char** argv) {
    cerr << "usage: " << argv[0] << " index [options] <graph1.vg> [graph2.vg ...]" << endl
         << "Creates an index on the specified graph or graphs. All graphs indexed must " << endl
         << "already be in a joint ID space." << endl
         << "general options:" << endl
         << "  -h, --help                print this help message to stderr and exit" << endl
         << "  -b, --temp-dir DIR        use DIR for temporary files" << endl
         << "  -t, --threads N           number of threads to use" << endl
         << "  -p, --progress            show progress" << endl
         << "xg options:" << endl
         << "  -x, --xg-name FILE        use this file to store a succinct, queryable version" << endl
         << "                            of graph(s), or read for GCSA or distance indexing" << endl
         << "  -L, --xg-alts             include alt paths in XG" << endl
         << "gcsa options:" << endl
         << "  -g, --gcsa-out FILE       output GCSA2 (FILE) & LCP (FILE.lcp) indexes" << endl
       //<< "  -i, --dbg-in FILE         use kmers from FILE instead of input VG (may repeat)" << endl
         << "  -f, --mapping FILE        use this node mapping in GCSA2 construction" << endl
         << "  -k, --kmer-size N         index kmers of size N in the graph [" << gcsa::Key::MAX_LENGTH << "]" << endl
         << "  -X, --doubling-steps N    use N doubling steps for GCSA2 construction "
                                     << "[" << gcsa::ConstructionParameters::DOUBLING_STEPS << "]" << endl
         << "  -Z, --size-limit N        limit temp disk space usage to N GB "
                                     << "[" << gcsa::ConstructionParameters::SIZE_LIMIT << "]" << endl
         << "  -V, --verify-index        validate the GCSA2 index using the input kmers" << endl
         << "                            (important for testing)" << endl
         << "gam indexing options:" << endl
         << "  -l, --index-sorted-gam    input is sorted .gam format alignments," << endl
         << "                            store a GAI index of the sorted GAM in INPUT.gam.gai" << endl
         << "vg in-place indexing options:" << endl
         << "      --index-sorted-vg     input is ID-sorted .vg format graph chunks" << endl
         << "                            store a VGI index of the sorted vg in INPUT.vg.vgi" << endl
         << "snarl distance index options" << endl
         << "  -j, --dist-name FILE      use this file to store a snarl-based distance index" << endl
         << "      --snarl-limit N       don't store distances for snarls > N nodes [10000]" << endl
         << "                            if 0 then don't store distances, only the snarl tree" << endl
         << "      --no-nested-distance  only store distances along the top-level chain" << endl
         << "  -w, --upweight-node N     upweight the node with ID N to push it to be part" << endl
         << "                            of a top-level chain (may repeat)" << endl;
}

int main_index(int argc, char** argv) {

    if (argc == 2) {
        help_index(argv);
        return 1;
    }

    constexpr int OPT_BUILD_VGI_INDEX = 1000;
    constexpr int OPT_RENAME_VARIANTS = 1001;
    constexpr int OPT_DISTANCE_SNARL_LIMIT = 1002;
    constexpr int OPT_DISTANCE_NESTING = 1003;

    // Which indexes to build.
    bool build_xg = false, build_gcsa = false, build_dist = false;

    // Files we should read.
    string vcf_name, mapping_name;
    vector<string> dbg_names;

    // Files we should write.
    string xg_name, gcsa_name, dist_name;

    // General
    bool show_progress = false;

    // GCSA
    gcsa::size_type kmer_size = gcsa::Key::MAX_LENGTH;
    gcsa::ConstructionParameters params;
    bool verify_gcsa = false;
    
    // Gam index (GAI)
    bool build_gai_index = false;
    
    // VG in-place index (VGI)
    bool build_vgi_index = false;

    // Include alt paths in xg
    bool xg_alts = false;

    //Distance index
    size_t snarl_limit = 50000;
    bool only_top_level_chain_distances = false;
    std::unordered_map<nid_t, size_t> extra_node_weight;
    // We will put this amount of extra weight on upweighted nodes. It should
    // be longer than the maximum plausible distracting path or spurious bridge
    // edge cycle, but small enough that several of it fit in a size_t.
    // TODO: Expose to command line.
    constexpr size_t EXTRA_WEIGHT = 10000000000;

    int c;
    optind = 2; // force optind past command positional argument
    while (true) {
        static struct option long_options[] =
        {
            // General
            {"temp-dir", required_argument, 0, 'b'},
            {"threads", required_argument, 0, 't'},
            {"progress",  no_argument, 0, 'p'},
            {"help",  no_argument, 0, 'h'},

            // XG
            {"xg-name", required_argument, 0, 'x'},
            {"xg-alts", no_argument, 0, 'L'},

            // GBWT. These have been removed and will return an error.
            {"vcf-phasing", required_argument, 0, 'v'},
            {"ignore-missing", no_argument, 0, 'W'},
            {"store-threads", no_argument, 0, 'T'},
            {"store-gam", required_argument, 0, 'M'},
            {"store-gaf", required_argument, 0, 'F'},
            {"gbwt-name", required_argument, 0, 'G'},
            {"actual-phasing", no_argument, 0, 'z'},
            {"force-phasing", no_argument, 0, 'P'},
            {"discard-overlaps", no_argument, 0, 'o'},
            {"batch-size", required_argument, 0, 'B'},
            {"buffer-size", required_argument, 0, 'u'},
            {"id-interval", required_argument, 0, 'n'},
            {"range", required_argument, 0, 'R'},
            {"rename", required_argument, 0, 'r'},
            {"rename-variants", no_argument, 0, OPT_RENAME_VARIANTS},
            {"region", required_argument, 0, 'I'},
            {"exclude", required_argument, 0, 'E'},

            // GCSA
            {"gcsa-out", required_argument, 0, 'g'},
            {"dbg-in", required_argument, 0, 'i'},
            {"mapping", required_argument, 0, 'f'},
            {"kmer-size", required_argument, 0, 'k'},
            {"doubling-steps", required_argument, 0, 'X'},
            {"size-limit", required_argument, 0, 'Z'},
            {"verify-index", no_argument, 0, 'V'},
            
            // GAM index (GAI)
            {"index-sorted-gam", no_argument, 0, 'l'},
            
            // VG in-place index (VGI)
            {"index-sorted-vg", no_argument, 0, OPT_BUILD_VGI_INDEX},

            //Snarl distance index
            {"snarl-limit", required_argument, 0, OPT_DISTANCE_SNARL_LIMIT},
            {"dist-name", required_argument, 0, 'j'},
            {"no-nested-distance", no_argument, 0, OPT_DISTANCE_NESTING},
            {"upweight-node", required_argument, 0, 'w'},
            {0, 0, 0, 0}
        };

        int option_index = 0;
        c = getopt_long (argc, argv, "b:t:px:Lv:WTM:F:G:zPoB:u:n:R:r:I:E:g:i:f:k:X:Z:Vlj:w:h?",
                         long_options, &option_index);

        // Detect the end of the options.
        if (c == -1)
            break;

        switch (c)
        {
        // General
        case 'b':
            temp_file::set_dir(optarg);
            break;
        case 't':
            set_thread_count(context, optarg);
            break;
        case 'p':
            show_progress = true;
            break;

        // XG
        case 'x':
            build_xg = true;
            // This may be an input *or* output
            xg_name = optarg;
            break;
        case 'L':
            xg_alts = true;
            break;

        // GBWT. The options remain, but they are no longer supported.
        case 'v': // Fall through
        case 'W': // Fall through
        case 'T': // Fall through
        case 'M': // Fall through
        case 'F': // Fall through
        case 'G': // Fall through
        case 'z': // Fall through
        case 'P': // Fall through
        case 'o': // Fall through
        case 'B': // Fall through
        case 'u': // Fall through
        case 'n': // Fall through
        case 'R': // Fall through
        case 'r': // Fall through
        case OPT_RENAME_VARIANTS: // Fall through
        case 'I': // Fall through
        case 'E':
            fatal_error(context) << "GBWT construction options have been removed; use vg gbwt instead" << endl;
            break;

        // GCSA
        case 'g':
            build_gcsa = true;
            gcsa_name = ensure_writable(context, optarg);
            // We also write to gcsa_name + ".lcp"
            ensure_writable(context, gcsa_name + ".lcp");
            break;
        case 'i':
            warning(context) << "-i option is deprecated" << endl;
            dbg_names.push_back(optarg);
            break;
        case 'f':
            mapping_name = require_exists(context, optarg);
            break;
        case 'k':
            kmer_size = std::max(parse<size_t>(optarg), 1ul);
            break;
        case 'X':
            params.setSteps(parse<size_t>(optarg));
            break;
        case 'Z':
            params.setLimit(parse<size_t>(optarg));
            break;
        case 'V':
            verify_gcsa = true;
            break;
            
        // Gam index (GAI)
        case 'l':
            build_gai_index = true;
            break;
            
        // VGI index
        case OPT_BUILD_VGI_INDEX:
            build_vgi_index = true;
            break;

        //Snarl distance index
        case 'j':
            build_dist = true;
            dist_name = ensure_writable(context, optarg);
            break;
        case OPT_DISTANCE_SNARL_LIMIT:
            snarl_limit = parse<int>(optarg);
            break;
        case OPT_DISTANCE_NESTING:
            only_top_level_chain_distances = true;
            break;
        case 'w':
            // We use += so you can repeat a node and make it even more
            // heavier.
            extra_node_weight[parse<nid_t>(optarg)] += EXTRA_WEIGHT;
            break;

        case 'h':
        case '?':
            help_index(argv);
            exit(1);
            break;
        default:
            abort ();
        }
    }

    vector<string> file_names;
    while (optind < argc) {
        string file_name = get_input_file_name(optind, argc, argv);
        file_names.push_back(file_name);
    }


    if (xg_name.empty() && gcsa_name.empty() && !build_gai_index && !build_vgi_index && dist_name.empty()) {
        fatal_error(context) << "index type not specified" << endl;
    }

    if (file_names.size() <= 0 && dbg_names.empty()){
        //fatal_error(context) << "No graph provided for indexing. "
        //                     << "Please provide a .vg file or GCSA2-format deBruijn graph to index." << endl;
    }
    
    if (file_names.size() != 1 && build_gai_index) {
        fatal_error(context) << "can only index exactly one sorted GAM file at a time" << endl;
    }
    
    if (file_names.size() != 1 && build_vgi_index) {
        fatal_error(context) << "can only index exactly one sorted VG file at a time" << endl;
    }
    
    if (file_names.size() > 1 && build_dist) {
        // Allow zero filenames for the index-from-xg mode
        fatal_error(context) << "can only create one distance index at a time" << endl;
    }
    
    if (build_gcsa && kmer_size > gcsa::Key::MAX_LENGTH) {
        fatal_error(context) << "GCSA2 cannot index with kmer size greater than "
                             << gcsa::Key::MAX_LENGTH << endl;
    }

    if (!build_dist && !extra_node_weight.empty()) {
        fatal_error(context) << "cannot up-weight nodes for snarl finding if not building distance index" << endl;
    }
    
    if (build_xg && build_gcsa && file_names.empty()) {
        // Really we want to build a GCSA by *reading* an XG
        build_xg = false;
        // We'll continue in the build_gcsa section
        warning(context) << "providing input XG with option -x is deprecated" << endl;
    }
    if (build_dist && file_names.empty()) {
        //If we want to build the distance index from the xg
        build_xg = false;
        warning(context) << "providing input XG with option -x is deprecated" << endl;
    }


    // Build XG. Include alt paths in the XG if requested with -L.
    if (build_xg) {
        ensure_writable(context, xg_name);
        if (file_names.empty()) {
            // VGset or something segfaults when we feed it no graphs.
            fatal_error(context) << "at least one graph is required to build an XG index" << endl;
        }
        if (show_progress) {
            basic_log(context) << "Building XG index" << endl;
        }
        xg::XG xg_index;
        VGset graphs(file_names);
        graphs.to_xg(xg_index, (xg_alts ? [](const string&) {return false;} : Paths::is_alt), nullptr);
        if (show_progress) {
            basic_log(context) << "Saving XG index to " << xg_name << endl;
        }
        // Save the XG.
        vg::io::save_handle_graph(&xg_index, xg_name);
    }

    // Build GCSA
    if (build_gcsa) {

        // Configure GCSA2 verbosity so it doesn't spit out loads of extra info
        if (!show_progress) {
            gcsa::Verbosity::set(gcsa::Verbosity::SILENT);
        }

        double start = gcsa::readTimer();

        // Generate temporary kmer files
        bool delete_kmer_files = false;
        if (dbg_names.empty()) {
            if (show_progress) {
                basic_log(context) << "Generating kmer files..." << endl;
            }
            
            if (!file_names.empty()) {
                // Get the kmers from a VGset.
                VGset graphs(file_names);
                size_t kmer_bytes = params.getLimitBytes();
                dbg_names = graphs.write_gcsa_kmers_binary(kmer_size, kmer_bytes);
                params.reduceLimit(kmer_bytes);
                delete_kmer_files = true;
            } else if (!xg_name.empty()) {
                // Get the kmers from an XG or other single graph
                
                // Load the graph
                require_exists(context, xg_name);
                auto single_graph = vg::io::VPKG::load_one<HandleGraph>(xg_name);
                
                auto make_kmers_for_component = [&](const HandleGraph* g) {
                    // Make an overlay on it to add source and sink nodes
                    // TODO: Don't use this directly; unify this code with VGset's code.
                    SourceSinkOverlay overlay(g, kmer_size);
                    
                    // Get the size limit
                    size_t kmer_bytes = params.getLimitBytes();
                    
                    // Write the kmer temp file
                    dbg_names.push_back(write_gcsa_kmers_to_tmpfile(overlay, kmer_size, kmer_bytes,
                        overlay.get_id(overlay.get_source_handle()),
                        overlay.get_id(overlay.get_sink_handle())));
                        
                    // Feed back into the size limit
                    params.reduceLimit(kmer_bytes);
                    delete_kmer_files = true;
                };
                
                if (show_progress) {
                    basic_log(context) << "Finding connected components..." << endl;
                }
                
                // Get all the components in the graph, which we can process separately to save memory.
                std::vector<std::unordered_set<nid_t>> components = \
                    handlealgs::weakly_connected_components(single_graph.get());
                
                if (components.size() == 1) {
                    // Only one component
                    if (show_progress) {
                        basic_log(context) << "Processing single component graph..." << endl;
                    }
                    make_kmers_for_component(single_graph.get());
                } else {
                    for (size_t i = 0; i < components.size(); i++) {
                        // Run separately on each component.
                        // Don't run in parallel or size limit tracking won't work.
                
                        if (show_progress) {
                            basic_log(context) << "Selecting component " 
                                               << i << "/" << components.size() << "..." << endl;
                        }
                        
                        bdsg::PackedSubgraphOverlay component_graph(single_graph.get());
                        for (auto& id : components[i]) {
                            // Add each node to the subgraph.
                            // TODO: use a handle-returning component
                            // finder so we don't need to get_handle here.
                            component_graph.add_node(single_graph->get_handle(id, false));
                        }
                        
                        if (show_progress) {
                            cerr << "Processing component " << i << "/" << components.size() << "..." << endl;
                        }
                        
                        make_kmers_for_component(&component_graph);
                    }
                }
            } else {
                fatal_error(context) << "cannot generate GCSA index without either a VG or an XG" << endl;
            }
        }

        // Build the index
        if (show_progress) {
            basic_log(context) << "Building the GCSA2 index..." << endl;
        }
        gcsa::InputGraph input_graph(dbg_names, true, params, gcsa::Alphabet(), mapping_name);
        gcsa::GCSA gcsa_index(input_graph, params);
        gcsa::LCPArray lcp_array(input_graph, params);
        if (show_progress) {
            double seconds = gcsa::readTimer() - start;
            basic_log(context) << "GCSA2 index built in " << seconds << " seconds, "
                               << gcsa::inGigabytes(gcsa::memoryUsage()) << " GB" << endl;
            basic_log(context) << "I/O volume:" << gcsa::inGigabytes(gcsa::readVolume()) << " GB read, "
                               << gcsa::inGigabytes(gcsa::writeVolume()) << " GB write" << endl;
        }

        // Save the indexes
        save_gcsa(gcsa_index, gcsa_name, show_progress);
        save_lcp(lcp_array, gcsa_name + ".lcp", show_progress);

        // Verify the index
        if (verify_gcsa) {
            if (show_progress) {
                basic_log(context) << "Verifying the index..." << endl;
            }
            if (!gcsa::verifyIndex(gcsa_index, &lcp_array, input_graph)) {
                warning(context) << "GCSA2 index verification failed" << endl;
            }
        }

        // Delete the temporary kmer files
        if (delete_kmer_files) {
            for (auto& filename : dbg_names) {
                temp_file::remove(filename);
            }
        }
    }
    
    if (build_gai_index) {
        // Index a sorted GAM file.
        
        get_input_file(file_names.at(0), [&](istream& in) {
            // Grab the input GAM stream and wrap it in a cursor
            vg::io::ProtobufIterator<Alignment> cursor(in);
            
            // Index the file
            StreamIndex<Alignment> index;
            index.index(cursor);
 
            // Save the GAM index in the appropriate place.
            // TODO: Do we really like this enforced naming convention just beacuse samtools does it?
            ofstream index_out(file_names.at(0) + ".gai");
            if (!index_out.good()) {
                fatal_error(context) << "could not open " << file_names.at(0) << ".gai for writing" << endl;
            }
            index.save(index_out);
        });
    }
    
    if (build_vgi_index) {
        // Index an ID-sorted VG file.
        get_input_file(file_names.at(0), [&](istream& in) {
            // Grab the input VG stream and wrap it in a cursor
            vg::io::ProtobufIterator<Graph> cursor(in);
            
            // Index the file
            StreamIndex<Graph> index;
            index.index(cursor);
 
            // Save the index in the appropriate place.
            // TODO: Do we really like this enforced naming convention just beacuse samtools does it?
            string index_name = file_names.at(0) + ".vgi";
            ensure_writable(context, index_name);
            ofstream index_out(index_name);
            index.save(index_out);
        });
        
    }

    //Build a snarl-based minimum distance index
    if (build_dist) {
        if (file_names.empty() && xg_name.empty()) {
            fatal_error(context) << "one graph is required to build a distance index" << endl;
        } else if (file_names.size() > 1 || (file_names.size() == 1 && !xg_name.empty())) {
            fatal_error(context) << "only one graph at a time can be used to build a distance index" << endl;
        } else if (dist_name.empty()) {
            fatal_error(context) << "distance index requires an output file" << endl;
            
        } else  {
            //Get graph and build dist index

            if (file_names.empty() && !xg_name.empty()) {
                // We were given a -x specifically to read as XG
                
                auto xg = vg::io::VPKG::load_one<xg::XG>(xg_name);

                IntegratedSnarlFinder snarl_finder(*xg.get(), extra_node_weight);
                // Create the SnarlDistanceIndex
                SnarlDistanceIndex distance_index;

                //Fill it in
                fill_in_distance_index(&distance_index, xg.get(), &snarl_finder, snarl_limit, only_top_level_chain_distances, false);
                // Save it
                distance_index.serialize(dist_name);
            } else {
                // May be GBZ or a HandleGraph.
                auto options = vg::io::VPKG::try_load_first<gbwtgraph::GBZ, handlegraph::HandleGraph>(file_names.at(0));
                
                if (get<0>(options)) {
                    // We have a GBZ graph
                    auto& gbz = get<0>(options);
                    
                    // Create the SnarlDistanceIndex
                    IntegratedSnarlFinder snarl_finder(gbz->graph, extra_node_weight);

                    //Make a distance index and fill it in
                    SnarlDistanceIndex distance_index;
                    fill_in_distance_index(&distance_index, &(gbz->graph), &snarl_finder, snarl_limit, only_top_level_chain_distances, false);
                    // Save it
                    distance_index.serialize(dist_name);
                } else if (get<1>(options)) {
                    // We were given a graph generically
                    auto& graph = get<1>(options);
                    
                    // Create the SnarlDistanceIndex
                    IntegratedSnarlFinder snarl_finder(*graph.get(), extra_node_weight);

                    //Make a distance index and fill it in
                    SnarlDistanceIndex distance_index;
                    fill_in_distance_index(&distance_index, graph.get(), &snarl_finder, snarl_limit, only_top_level_chain_distances, false);
                    // Save it
                    distance_index.serialize(dist_name);
                } else {
                    fatal_error(context) << "input is not a graph or GBZ" << endl;
                }
            }
        }

    }
    if (show_progress) {
        basic_log(context) << "Memory usage: " << gcsa::inGigabytes(gcsa::memoryUsage()) << " GB" << endl;
    }
    return 0;
}

// Register subcommand
static Subcommand vg_construct("index", "index graphs or alignments for random access or mapping", PIPELINE, 4, main_index);
