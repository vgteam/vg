/**
 * \file cluster_main.cpp: experimental snarl clustering test harness
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

#include "../snarl_seed_clusterer.hpp"
#include "../mapper.hpp"
#include "../annotation.hpp"
#include "../xg.hpp"
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

void help_cluster(char** argv) {
    cerr
    << "usage: " << argv[0] << " cluster [options] input.gam > output.gam" << endl
    << "Find and cluster mapping seeds." << endl
    << endl
    << "basic options:" << endl
    << "  -x, --xg-name FILE            use this xg index or graph (required)" << endl
    << "  -g, --gcsa-name FILE          use this GCSA2/LCP index pair (both FILE and FILE.lcp)" << endl
    << "  -m, --minimizer-name FILE     use this minimizer index" << endl
    << "  -d, --dist-name FILE          cluster using this distance index (required)" << endl
    << "  -c, --hit-cap INT             ignore minimizers with more than this many locations [10]" << endl
    << "computational parameters:" << endl
    << "  -t, --threads INT             number of compute threads to use" << endl;
}

int main_cluster(int argc, char** argv) {

    if (argc == 2) {
        help_cluster(argv);
        return 1;
    }

    // initialize parameters with their default options
    string xg_name;
    string gcsa_name;
    string minimizer_name;
    string distance_name;
    // How close should two hits be to be in the same cluster?
    size_t distance_limit = 1000;
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
                    cerr << "error:[vg cluster] Must provide XG file with -x." << endl;
                    exit(1);
                }
                break;
                
            case 'g':
                gcsa_name = optarg;
                if (gcsa_name.empty()) {
                    cerr << "error:[vg cluster] Must provide GCSA file with -g." << endl;
                    exit(1);
                }
                break;
            
            case 'm':
                minimizer_name = optarg;
                if (minimizer_name.empty()) {
                    cerr << "error:[vg cluster] Must provide minimizer file with -m." << endl;
                    exit(1);
                }
                break;
                
            case 'd':
                distance_name = optarg;
                if (distance_name.empty()) {
                    cerr << "error:[vg cluster] Must provide distance index file with -d." << endl;
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
                    cerr << "error:[vg cluster] Thread count (-t) set to " << num_threads << ", must set to a positive integer." << endl;
                    exit(1);
                }
                omp_set_num_threads(num_threads);
            }
                break;
                
            case 'h':
            case '?':
            default:
                help_cluster(argv);
                exit(1);
                break;
        }
    }
    
    
    if (xg_name.empty()) {
        cerr << "error:[vg cluster] Finding clusters requires an XG index, must provide XG file (-x)" << endl;
        exit(1);
    }
    
    if (gcsa_name.empty() && minimizer_name.empty()) {
        cerr << "error:[vg cluster] Finding clusters requires a GCSA2 index or minimizer index (-g, -m)" << endl;
        exit(1);
    }
    
    
    if (distance_name.empty()) {
        cerr << "error:[vg cluster] Finding clusters requires a distance index, must provide distance index file (-d)" << endl;
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
    
    // Make the clusterer
    NewSnarlSeedClusterer clusterer(*distance_index);
    
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
        // We want to profile the clustering and the code around it.
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
                        // The minimizer is infrequent enough to be informative, so feed it into clustering
                        
                        // Locate it in the graph. We do not have to reverse the hits for a
                        // reverse minimizers, as the clusterer only cares about node ids.
                        for (auto& hit : minimizer_index->find(minimizers[i])) {
                            // For each position, remember it and what minimizer it came from
                            seeds.push_back(hit.first);
                            seed_to_source.push_back(i);
                        }
                    }
                }
                
            }
            vector<NewSnarlSeedClusterer::Seed> seed_clusters;
            for (pos_t pos : seeds) {
                seed_clusters.emplace_back();
                seed_clusters.back().pos = pos;
            }

            
            // Cluster the seeds. Get sets of input seed indexes that go together.
            // Make sure to time it.
            std::chrono::time_point<std::chrono::system_clock> start = std::chrono::system_clock::now();
            vector<NewSnarlSeedClusterer::Cluster> clusters = clusterer.cluster_seeds(seed_clusters, distance_limit);
            std::chrono::time_point<std::chrono::system_clock> end = std::chrono::system_clock::now();
            std::chrono::duration<double> elapsed_seconds = end-start;
            
            // Compute the covered portion of the read represented by each cluster
            vector<double> read_coverage_by_cluster;
            for (auto& cluster : clusters) {
                // We set bits in here to true when query anchors cover them
                vector<bool> covered(aln.sequence().size());
                // We use this to convert iterators to indexes
                auto start = aln.sequence().begin();
                
                for (auto hit_index : cluster.seeds) {
                    // For each hit in the cluster, work out what anchor sequence it is from.
                    size_t source_index = seed_to_source.at(hit_index);
                    
                    if (mapper) {
                        // Using MEMs
                        for (size_t i = (mems[source_index].begin - start); i < (mems[source_index].end - start); i++) {
                            // Set all the bits in read space for that MEM
                            covered[i] = true;
                        }
                    } else {
                        // Using minimizers
                        // The offset of a reverse minimizer is the endpoint of the kmer
                        size_t start_offset = minimizers[source_index].offset;
                        if (minimizers[source_index].is_reverse) {
                            start_offset = start_offset + 1 - minimizer_index->k();
                        }
                        for (size_t i = start_offset; i < start_offset + minimizer_index->k(); i++) {
                            // Set all the bits in read space for that minimizer.
                            // Each minimizr is a length-k exact match starting at a position
                            covered[i] = true;
                        }
                    }
                }
                
                // Count up the covered positions
                size_t covered_count = 0;
                for (auto bit : covered) {
                    covered_count += bit;
                }
                
                // Turn that into a fraction
                read_coverage_by_cluster.push_back(covered_count / (double) covered.size());
            }
            
            // Make a vector of cluster indexes to sort
            vector<size_t> cluster_indexes_in_order;
            for (size_t i = 0; i < clusters.size(); i++) {
                cluster_indexes_in_order.push_back(i);
            }
        
            // Put the most covering cluster's index first
            std::sort(cluster_indexes_in_order.begin(), cluster_indexes_in_order.end(), [&](const size_t& a, const size_t& b) -> bool {
                // Return true if a must come before b, and false otherwise
                return read_coverage_by_cluster.at(a) > read_coverage_by_cluster.at(b);
            });
            
            // Find the seeds in the clusters tied for best.
            vector<pos_t> best;
            if (!clusters.empty()) {
                // How much does the best cluster cover
                double best_coverage = read_coverage_by_cluster.at(cluster_indexes_in_order.front());
                for (size_t i = 0; i < cluster_indexes_in_order.size() &&
                    read_coverage_by_cluster.at(cluster_indexes_in_order[i]) >= best_coverage; i++) {
                    
                    // For each cluster covering that much or more of the read
                    for (auto seed_index : clusters.at(cluster_indexes_in_order[i]).seeds) {
                        // For each seed in those clusters
                        
                        // Mark that seed as being part of the best cluster(s)
                        best.push_back(seeds.at(seed_index));
                    }
                    
                }
                
            }
            
            // Decide if they are in the right place for the original alignment or not
            unordered_set<vg::id_t> true_nodes;
            for (auto& mapping : aln.path().mapping()) {
                true_nodes.insert(mapping.position().node_id());
            }
            // We are in the right place if we share any nodes
            bool have_overlap = false;
            for (auto& pos : best) {
                if (true_nodes.count(get_id(pos))) {
                    // The cluster had a position on a node that the real alignment had.
                    have_overlap = true;
                }
            }
            
            // We also want to know if we overlap any non-filtered hit
            bool have_hit_overlap = false;
            for (auto& pos : seeds) {
                if (true_nodes.count(get_id(pos))) {
                    // The hit set had a position on a node that the real alignment had.
                    have_hit_overlap = true;
                }
            }
            
            // And we need a vector of cluster sizes
            vector<double> cluster_sizes;
            cluster_sizes.reserve(clusters.size());
            for (auto& cluster : clusters) {
                cluster_sizes.push_back((double)cluster.seeds.size());
            }
            
            // Tag the alignment with cluster accuracy
            set_annotation(aln, "best_cluster_overlap", have_overlap);
            // And with any-hit overlap
            set_annotation(aln, "any_seed_overlap", have_hit_overlap);
            // And with cluster time
            set_annotation(aln, "cluster_seconds", elapsed_seconds.count());
            // And with hit count clustered
            set_annotation(aln, "seed_count", (double)seeds.size());
            // And with cluster count returned
            set_annotation(aln, "cluster_count", (double)clusters.size());
            // And with size of each cluster
            set_annotation(aln, "cluster_sizes", cluster_sizes);
            // And with the coverage of the read in the best cluster
            set_annotation(aln, "best_cluster_coverage", clusters.empty() ? 0.0 :
                read_coverage_by_cluster.at(cluster_indexes_in_order.front()));
            
            
            // TODO: parallelize this
            #pragma omp critical (cout)
            emitter.write(std::move(aln));
        });
    });
    
    return 0;
}

// Register subcommand
static Subcommand vg_cluster("cluster", "find and cluster mapping seeds", DEVELOPMENT, main_cluster);


