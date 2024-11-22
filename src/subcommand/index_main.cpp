// index.cpp: define the "vg index" subcommand, which makes xg, GCSA2, and distance indexes

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

void help_index(char** argv) {
    cerr << "usage: " << argv[0] << " index [options] <graph1.vg> [graph2.vg ...]" << endl
         << "Creates an index on the specified graph or graphs. All graphs indexed must " << endl
         << "already be in a joint ID space." << endl
         << "general options:" << endl
         << "    -b, --temp-dir DIR     use DIR for temporary files" << endl
         << "    -t, --threads N        number of threads to use" << endl
         << "    -p, --progress         show progress" << endl
         << "xg options:" << endl
         << "    -x, --xg-name FILE     use this file to store a succinct, queryable version of the graph(s), or read for GCSA or distance indexing" << endl
         << "    -L, --xg-alts          include alt paths in xg" << endl
         << "gcsa options:" << endl
         << "    -g, --gcsa-out FILE    output a GCSA2 index to the given file" << endl
         //<< "    -i, --dbg-in FILE      use kmers from FILE instead of input VG (may repeat)" << endl
         << "    -f, --mapping FILE     use this node mapping in GCSA2 construction" << endl
         << "    -k, --kmer-size N      index kmers of size N in the graph (default " << gcsa::Key::MAX_LENGTH << ")" << endl
         << "    -X, --doubling-steps N use this number of doubling steps for GCSA2 construction (default " << gcsa::ConstructionParameters::DOUBLING_STEPS << ")" << endl
         << "    -Z, --size-limit N     limit temporary disk space usage to N gigabytes (default " << gcsa::ConstructionParameters::SIZE_LIMIT << ")" << endl
         << "    -V, --verify-index     validate the GCSA2 index using the input kmers (important for testing)" << endl
         << "gam indexing options:" << endl
         << "    -l, --index-sorted-gam input is sorted .gam format alignments, store a GAI index of the sorted GAM in INPUT.gam.gai" << endl
         << "vg in-place indexing options:" << endl
         << "    --index-sorted-vg      input is ID-sorted .vg format graph chunks, store a VGI index of the sorted vg in INPUT.vg.vgi" << endl
         << "snarl distance index options" << endl
         << "    -j  --dist-name FILE   use this file to store a snarl-based distance index" << endl
         << "        --snarl-limit N    don't store snarl distances for snarls with more than N nodes (default 10000)" << endl
         << "                           if N is 0 then don't store distances, only the snarl tree" << endl;
}

int main_index(int argc, char** argv) {

    if (argc == 2) {
        help_index(argv);
        return 1;
    }

    #define OPT_BUILD_VGI_INDEX  1000
    #define OPT_RENAME_VARIANTS  1001
    #define OPT_DISTANCE_SNARL_LIMIT 1002

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

    int c;
    optind = 2; // force optind past command positional argument
    while (true) {
        static struct option long_options[] =
        {
            // General
            {"temp-dir", required_argument, 0, 'b'},
            {"threads", required_argument, 0, 't'},
            {"progress",  no_argument, 0, 'p'},

            // XG
            {"xg-name", required_argument, 0, 'x'},
            {"thread-db", required_argument, 0, 'F'},
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
            {0, 0, 0, 0}
        };

        int option_index = 0;
        c = getopt_long (argc, argv, "b:t:px:Lv:WTM:F:G:zPoB:u:n:R:r:I:E:g:i:f:k:X:Z:Vlj:h",
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
            omp_set_num_threads(parse<int>(optarg));
            break;
        case 'p':
            show_progress = true;
            break;

        // XG
        case 'x':
            build_xg = true;
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
            std::cerr << "error: [vg index] GBWT construction options have been removed; use vg gbwt instead" << std::endl;
            std::exit(EXIT_FAILURE);
            break;

        // GCSA
        case 'g':
            build_gcsa = true;
            gcsa_name = optarg;
            break;
        case 'i':
            cerr << "warning: -i option is deprecated" << endl;
            dbg_names.push_back(optarg);
            break;
        case 'f':
            mapping_name = optarg;
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
            dist_name = optarg;
            break;
        case OPT_DISTANCE_SNARL_LIMIT:
            snarl_limit = parse<int>(optarg);
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
        cerr << "error: [vg index] index type not specified" << endl;
        return 1;
    }

    if (file_names.size() <= 0 && dbg_names.empty()){
        //cerr << "No graph provided for indexing. Please provide a .vg file or GCSA2-format deBruijn graph to index." << endl;
        //return 1;
    }
    
    if (file_names.size() != 1 && build_gai_index) {
        cerr << "error: [vg index] can only index exactly one sorted GAM file at a time" << endl;
        return 1;
    }
    
    if (file_names.size() != 1 && build_vgi_index) {
        cerr << "error: [vg index] can only index exactly one sorted VG file at a time" << endl;
        return 1;
    }
    
    if (file_names.size() > 1 && build_dist) {
        // Allow zero filenames for the index-from-xg mode
        cerr << "error: [vg index] can only create one distance index at a time" << endl;
        return 1;
    }
    
    if (build_gcsa && kmer_size > gcsa::Key::MAX_LENGTH) {
        cerr << "error: [vg index] GCSA2 cannot index with kmer size greater than " << gcsa::Key::MAX_LENGTH << endl;
        return 1;
    }
    
    if (build_xg && build_gcsa && file_names.empty()) {
        // Really we want to build a GCSA by *reading* and XG
        build_xg = false;
        // We'll continue in the build_gcsa section
        std::cerr << "warning: [vg index] providing input XG with option -x is deprecated" << std::endl;
    }
    if (build_dist && file_names.empty()) {
        //If we want to build the distance index from the xg
        build_xg = false;
        std::cerr << "warning: [vg index] providing input XG with option -x is deprecated" << std::endl;
    }


    // Build XG. Include alt paths in the XG if requested with -L.
    if (build_xg) {
        if (file_names.empty()) {
            // VGset or something segfaults when we feed it no graphs.
            cerr << "error: [vg index] at least one graph is required to build an xg index" << endl;
            return 1;
        }
        if (show_progress) {
            cerr << "Building XG index" << endl;
        }
        xg::XG xg_index;
        VGset graphs(file_names);
        graphs.to_xg(xg_index, (xg_alts ? [](const string&) {return false;} : Paths::is_alt), nullptr);
        if (show_progress) {
            cerr << "Saving XG index to " << xg_name << endl;
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
                cerr << "Generating kmer files..." << endl;
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
                    cerr << "Finding connected components..." << endl;
                }
                
                // Get all the components in the graph, which we can process separately to save memory.
                std::vector<std::unordered_set<nid_t>> components = handlealgs::weakly_connected_components(single_graph.get());
                
                if (components.size() == 1) {
                    // Only one component
                    if (show_progress) {
                        cerr << "Processing single component graph..." << endl;
                    }
                    make_kmers_for_component(single_graph.get());
                } else {
                    for (size_t i = 0; i < components.size(); i++) {
                        // Run separately on each component.
                        // Don't run in parallel or size limit tracking won't work.
                
                        if (show_progress) {
                            cerr << "Selecting component " << i << "/" << components.size() << "..." << endl;
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
                cerr << "error: [vg index] cannot generate GCSA index without either a vg or an xg" << endl;
                exit(1);
            }
        }

        // Build the index
        if (show_progress) {
            cerr << "Building the GCSA2 index..." << endl;
        }
        gcsa::InputGraph input_graph(dbg_names, true, params, gcsa::Alphabet(), mapping_name);
        gcsa::GCSA gcsa_index(input_graph, params);
        gcsa::LCPArray lcp_array(input_graph, params);
        if (show_progress) {
            double seconds = gcsa::readTimer() - start;
            cerr << "GCSA2 index built in " << seconds << " seconds, "
                 << gcsa::inGigabytes(gcsa::memoryUsage()) << " GB" << endl;
            cerr << "I/O volume: " << gcsa::inGigabytes(gcsa::readVolume()) << " GB read, "
                 << gcsa::inGigabytes(gcsa::writeVolume()) << " GB write" << endl;
        }

        // Save the indexes
        save_gcsa(gcsa_index, gcsa_name, show_progress);
        save_lcp(lcp_array, gcsa_name + ".lcp", show_progress);

        // Verify the index
        if (verify_gcsa) {
            if (show_progress) {
                cerr << "Verifying the index..." << endl;
            }
            if (!gcsa::verifyIndex(gcsa_index, &lcp_array, input_graph)) {
                cerr << "warning: [vg index] GCSA2 index verification failed" << endl;
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
                cerr << "error: [vg index] could not open " << file_names.at(0) << ".gai" << endl;
                exit(1);
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
            ofstream index_out(file_names.at(0) + ".vgi");
            if (!index_out.good()) {
                cerr << "error: [vg index] could not open " << file_names.at(0) << ".vgi" << endl;
                exit(1);
            }
            index.save(index_out);
        });
        
    }

    //Build a snarl-based minimum distance index
    if (build_dist) {
        if (file_names.empty() && xg_name.empty()) {
            cerr << "error: [vg index] one graph is required to build a distance index" << endl;
            return 1;
        } else if (file_names.size() > 1 || (file_names.size() == 1 && !xg_name.empty())) {
            cerr << "error: [vg index] only one graph at a time can be used to build a distance index" << endl;
        } else if (dist_name.empty()) {
            cerr << "error: [vg index] distance index requires an output file" << endl;
            return 1;
            
        } else  {
            //Get graph and build dist index

            if (file_names.empty() && !xg_name.empty()) {
                // We were given a -x specifically to read as XG
                
                auto xg = vg::io::VPKG::load_one<xg::XG>(xg_name);

                IntegratedSnarlFinder snarl_finder(*xg.get());
                // Create the SnarlDistanceIndex
                SnarlDistanceIndex distance_index;

                //Fill it in
                fill_in_distance_index(&distance_index, xg.get(), &snarl_finder, snarl_limit, false);
                // Save it
                distance_index.serialize(dist_name);
            } else {
                // May be GBZ or a HandleGraph.
                auto options = vg::io::VPKG::try_load_first<gbwtgraph::GBZ, handlegraph::HandleGraph>(file_names.at(0));
                
                if (get<0>(options)) {
                    // We have a GBZ graph
                    auto& gbz = get<0>(options);
                    
                    // Create the SnarlDistanceIndex
                    IntegratedSnarlFinder snarl_finder(gbz->graph);

                    //Make a distance index and fill it in
                    SnarlDistanceIndex distance_index;
                    fill_in_distance_index(&distance_index, &(gbz->graph), &snarl_finder, snarl_limit);
                    // Save it
                    distance_index.serialize(dist_name);
                } else if (get<1>(options)) {
                    // We were given a graph generically
                    auto& graph = get<1>(options);
                    
                    // Create the SnarlDistanceIndex
                    IntegratedSnarlFinder snarl_finder(*graph.get());

                    //Make a distance index and fill it in
                    SnarlDistanceIndex distance_index;
                    fill_in_distance_index(&distance_index, graph.get(), &snarl_finder, snarl_limit);
                    // Save it
                    distance_index.serialize(dist_name);
                } else {
                    cerr << "error: [vg index] input is not a graph or GBZ" << endl;
                    return 1;
                }
            }
        }

    }
    if (show_progress) {
        cerr << "Memory usage: " << gcsa::inGigabytes(gcsa::memoryUsage()) << " GB" << endl;
    }
    return 0;
}

// Register subcommand
static Subcommand vg_construct("index", "index graphs or alignments for random access or mapping", PIPELINE, 4, main_index);
