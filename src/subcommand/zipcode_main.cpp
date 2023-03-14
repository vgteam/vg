/**
 * \file zipcode.cpp: experimental zipcode test harness
 */

#include <omp.h>
#include <unistd.h>
#include <getopt.h>
#include <iostream>
#include <cassert>
#include <vector>
#include <unordered_set>
#include <chrono>

#include "subcommand.hpp"

#include "../zip_code.hpp"
#include "../mapper.hpp"
#include "../annotation.hpp"
#include <vg/io/vpkg.hpp>
#include <vg/io/stream.hpp>
#include <vg/io/protobuf_emitter.hpp>


#include <gbwtgraph/minimizer.h>
#include <bdsg/overlays/overlay_helper.hpp>

//#define USE_CALLGRIND

#ifdef USE_CALLGRIND
#include <valgrind/callgrind.h>
#endif

using namespace std;
using namespace vg;
using namespace vg::subcommand;

void help_zipcode(char** argv) {
    cerr
    << "usage: " << argv[0] << " test zipcodes on minimizers from reads [options] input.gam > output.gam" << endl
    << endl
    << "basic options:" << endl
    << "  -x, --xg-name FILE            use this xg index or graph (required)" << endl
    << "  -m, --minimizer-name FILE     use this minimizer index" << endl
    << "  -d, --dist-name FILE          use this distance index (required)" << endl
    << "  -c, --hit-cap INT             ignore minimizers with more than this many locations [10]" << endl
    << "computational parameters:" << endl
    << "  -t, --threads INT             number of compute threads to use" << endl;
}

int main_zipcode(int argc, char** argv) {

    if (argc == 2) {
        help_zipcode(argv);
        return 1;
    }

    // initialize parameters with their default options
    string xg_name;
    string gcsa_name;
    string minimizer_name;
    string distance_name;
    size_t hit_cap = 10;
    
    int c;
    optind = 2; // force optind past command positional argument
    while (true) {
        static struct option long_options[] =
        {
            {"help", no_argument, 0, 'h'},
            {"xg-name", required_argument, 0, 'x'},
            {"gcsa-name", required_argument, 0, 'g'},
            {"minimizer-name", required_argument, 0, 'm'},
            {"dist-name", required_argument, 0, 'd'},
            {"hit-cap", required_argument, 0, 'c'},
            {"threads", required_argument, 0, 't'},
            {0, 0, 0, 0}
        };

        int option_index = 0;
        c = getopt_long (argc, argv, "hx:g:m:d:c:t:",
                         long_options, &option_index);


        // Detect the end of the options.
        if (c == -1)
            break;

        switch (c)
        {
            case 'x':
                xg_name = optarg;
                if (xg_name.empty()) {
                    cerr << "error:[vg zipcode] Must provide XG file with -x." << endl;
                    exit(1);
                }
                break;
                
            case 'g':
                gcsa_name = optarg;
                if (gcsa_name.empty()) {
                    cerr << "error:[vg zipcode] Must provide GCSA file with -g." << endl;
                    exit(1);
                }
                break;
            
            case 'm':
                minimizer_name = optarg;
                if (minimizer_name.empty()) {
                    cerr << "error:[vg zipcode] Must provide minimizer file with -m." << endl;
                    exit(1);
                }
                break;
                
            case 'd':
                distance_name = optarg;
                if (distance_name.empty()) {
                    cerr << "error:[vg zipcode] Must provide distance index file with -d." << endl;
                    exit(1);
                }
                break;
            
            case 'c':
                hit_cap = parse<size_t>(optarg);
                break;
                
            case 't':
            {
                int num_threads = parse<int>(optarg);
                if (num_threads <= 0) {
                    cerr << "error:[vg zipcode] Thread count (-t) set to " << num_threads << ", must set to a positive integer." << endl;
                    exit(1);
                }
                omp_set_num_threads(num_threads);
            }
                break;
                
            case 'h':
            case '?':
            default:
                help_zipcode(argv);
                exit(1);
                break;
        }
    }
    
    
    if (xg_name.empty()) {
        cerr << "error:[vg zipcode] Finding zipcodes requires an XG index, must provide XG file (-x)" << endl;
        exit(1);
    }
    
    if (gcsa_name.empty() && minimizer_name.empty()) {
        cerr << "error:[vg zipcode] Finding zipcodes requires a GCSA2 index or minimizer index (-g, -m)" << endl;
        exit(1);
    }
    
    
    if (distance_name.empty()) {
        cerr << "error:[vg zipcode] Finding zipcodes requires a distance index, must provide distance index file (-d)" << endl;
        exit(1);
    }
    
    // create in-memory objects
    unique_ptr<PathHandleGraph> path_handle_graph = vg::io::VPKG::load_one<PathHandleGraph>(xg_name);
    bdsg::PathPositionOverlayHelper overlay_helper;
    PathPositionHandleGraph* xg_index = overlay_helper.apply(path_handle_graph.get());
    unique_ptr<gcsa::GCSA> gcsa_index;
    unique_ptr<gcsa::LCPArray> lcp_index;
    if (!gcsa_name.empty()) {
        gcsa_index = vg::io::VPKG::load_one<gcsa::GCSA>(gcsa_name);
        lcp_index = vg::io::VPKG::load_one<gcsa::LCPArray>(gcsa_name + ".lcp");
    }
    unique_ptr<gbwtgraph::DefaultMinimizerIndex> minimizer_index;
    if (!minimizer_name.empty()) {
        minimizer_index = vg::io::VPKG::load_one<gbwtgraph::DefaultMinimizerIndex>(minimizer_name);
    }
    unique_ptr<SnarlDistanceIndex> distance_index = vg::io::VPKG::load_one<SnarlDistanceIndex>(distance_name);
    distance_index->preload(true);
    
    // Make a Mapper to look up MEM seeds
    unique_ptr<Mapper> mapper;
    if (gcsa_index) {
        // We will find MEMs using a Mapper
        mapper = make_unique<Mapper>(xg_index, gcsa_index.get(), lcp_index.get());
    }
    // Otherwise we will find minimizers using the minimizer_index
    
    get_input_file(optind, argc, argv, [&](istream& in) {
        // Open up the input GAM
        
        // Make the output emitter
        vg::io::ProtobufEmitter<Alignment> emitter(cout);
        
#ifdef USE_CALLGRIND
        // We want to profile the zipcodes and the code around it.
        CALLGRIND_START_INSTRUMENTATION;
#endif
        
        vg::io::for_each_parallel<Alignment>(in, [&](Alignment& aln) {
            // For each input alignment
            
            // We will find all the seed hits
            vector<pos_t> seeds;
            
            // If working with MEMs, this will hold all the MEMs
            vector<MaximalExactMatch> mems;
            // If working with minimizers, this will hold all the minimizers in the query
            vector<gbwtgraph::DefaultMinimizerIndex::minimizer_type> minimizers;
            // And either way this will map from seed to MEM or minimizer that generated it
            vector<size_t> seed_to_source;
            
            if (mapper) {
                // Find MEMs
                double lcp_avg, fraction_filtered;
                mems = mapper->find_mems_deep(aln.sequence().begin(), aln.sequence().end(), lcp_avg, fraction_filtered);
                
                // Convert to position seeds
                for (size_t i = 0; i < mems.size(); i++) {
                    auto& mem = mems[i];
                    for (gcsa::node_type n : mem.nodes) {
                        // Convert from GCSA node_type packing to a pos_t
                        seeds.push_back(make_pos_t(n));
                        // And remember which MEM the seed came from.
                        seed_to_source.push_back(i);
                    }
                }
            } else {
                // Find minimizers
                assert(minimizer_index);
                
                // Find minimizers in the query
                minimizers = minimizer_index->minimizers(aln.sequence());
                
                for (size_t i = 0; i < minimizers.size(); i++) {
                    // For each minimizer
                    if (hit_cap != 0 && minimizer_index->count(minimizers[i]) <= hit_cap) {
                        // The minimizer is infrequent enough to be informative
                        
                        // Locate it in the graph. We do not have to reverse the hits for a
                        // reverse minimizers, as the zipcodes only cares about node ids.
                        for (auto& hit : minimizer_index->find(minimizers[i])) {
                            // For each position, remember it and what minimizer it came from
                            seeds.push_back(hit.first);
                            seed_to_source.push_back(i);
                        }
                    }
                }
                
            }
            vector<double> elapsed_seconds_zip;
            vector<double> elapsed_seconds_index;
            vector<double> depths;
            vector<double> has_irregular_snarl;
            size_t count = 0;
            for (pos_t pos1 : seeds) {
                for (pos_t pos2 : seeds) {
                    count++;

                    //Get zip codes
                    ZipCode zip1;
                    zip1.fill_in_zipcode(*distance_index, pos1);
                    ZipCode zip2;
                    zip2.fill_in_zipcode(*distance_index, pos2);

                    //Time finding distance with the zip codes
                    std::chrono::time_point<std::chrono::system_clock> start = std::chrono::system_clock::now();
                    size_t zip_distance = ZipCode::minimum_distance_between(zip1, pos1, zip2, pos2, *distance_index);
                    std::chrono::time_point<std::chrono::system_clock> end = std::chrono::system_clock::now();
                    std::chrono::duration<double> elapsed_seconds = end-start;
                    elapsed_seconds_zip.emplace_back(elapsed_seconds.count());


                    //Time finding the distance with the index
                    start = std::chrono::system_clock::now();
                    size_t index_distance = minimum_distance(*distance_index, pos1, pos2);
                    end = std::chrono::system_clock::now();
                    elapsed_seconds = end-start;

                    elapsed_seconds_index.emplace_back(elapsed_seconds.count());
                    net_handle_t net1 = distance_index->get_node_net_handle(id(pos1));
                    net_handle_t net2 = distance_index->get_node_net_handle(id(pos2));
                    size_t depth = std::max(distance_index->get_depth(net1), 
                                            distance_index->get_depth(net2));
                    depths.emplace_back(depth);

                    bool is_irregular = false;
                    while(!distance_index->is_root(net1)){
                        if (distance_index->is_snarl(net1) && !distance_index->is_regular_snarl(net1)) {
                            is_irregular = true;
                        }
                        net1 = distance_index->get_parent(net1);
                    }
                    while(!distance_index->is_root(net2)){
                        if (distance_index->is_snarl(net2) && !distance_index->is_regular_snarl(net2)) {
                            is_irregular = true;
                        }
                        net2 = distance_index->get_parent(net2);
                    }
                    has_irregular_snarl.emplace_back(is_irregular);
                }
            }
                       
            // Tag the alignment times
            set_annotation(aln, "seconds_zip", elapsed_seconds_zip);
            set_annotation(aln, "seconds_index", elapsed_seconds_index);
            set_annotation(aln, "depths", depths);
            set_annotation(aln, "irregular", has_irregular_snarl);
            
            
            // TODO: parallelize this
            #pragma omp critical (cout)
            emitter.write(std::move(aln));
        });
    });
    
    return 0;
}

// Register subcommand
static Subcommand vg_zipcode("zipcode", "find distances between seeds using zipcodes", DEVELOPMENT, main_zipcode);


