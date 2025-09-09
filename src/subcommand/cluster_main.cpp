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
#include "../zip_code_tree.hpp"
#include "../mapper.hpp"
#include "../annotation.hpp"
#include "../xg.hpp"
#include "../minimizer_mapper.hpp"
#include "../index_registry.hpp"
#include <vg/io/vpkg.hpp>
#include <vg/io/stream.hpp>
#include <vg/io/protobuf_emitter.hpp>

#include <gbwtgraph/minimizer.h>
#include <gbwtgraph/gbwtgraph.h>
#include <gbwtgraph/gbz.h>
#include <bdsg/overlays/overlay_helper.hpp>

//#define USE_CALLGRIND

#ifdef USE_CALLGRIND
#include <valgrind/callgrind.h>
#endif

using namespace std;
using namespace vg;
using namespace vg::subcommand;

const string context = "[vg cluster]";

void help_cluster(char** argv) {
    cerr << "usage: " << argv[0] << " cluster [options] input.gam > output.gam" << endl
         << "Find and cluster mapping seeds." << endl
         << endl
         << "basic options:" << endl
         << "  -h, --help                    print this help message to stderr and exit" << endl
         << "  -x, --xg-name FILE            use this xg index or graph (required)" << endl
         << "  -f, --gbz-format              input graph is GBZ format (has graph & GBWT)" << endl
         << "  -g, --gcsa-name FILE          use FILE & FILE.lcp GCSA2/LCP index pair" << endl
         << "  -G, --gbwt-name FILE          use this GBWT" << endl
         << "  -B, --gbwtgraph-name FILE     use this GBWTGraph" << endl
         << "  -m, --minimizer-name FILE     use this minimizer index" << endl
         << "  -l, --long-reads              minimizer index is for long reads" << endl
         << "  -d, --dist-name FILE          cluster using this distance index (required)" << endl
         << "  -p, --prefix PREFIX           prefix to find indexes at" << endl
         << "  -c, --hit-cap INT             use all minimizers with at most INT hits [10]" << endl
         << "  -C, --hard-hit-cap INT        ignore minimizers with more than INT hits [500]" << endl
         << "  -a, --hits-above INT          output minimizers with at least INT hits [inf]" << endl
         << "  -S, --sequences-only          only output -s sequences of minimizers," << endl
         << "                                no clustering/ziptree (overrides -Z)" << endl
         << "  -F, --score-fraction FLOAT    select minimizers between -c/-C until" << endl
         << "                                score is FLOAT of total [0.9]" << endl
         << "  -U, --max-min INT             use at most INT minimizers, 0 for no limit [500]" << endl
         << "  -b, --num-bp-per-min INT      use maximum of number minimizers calculated by" << endl
         << "                                READ_LENGTH / INT and --max-min [1000]" << endl
         << "  -D, --downsample-min INT      downsample minimizers with windows of length" << endl
         << "                                read length/INT, 0 for no downsampling [0]" << endl
         << "  -z, --zip-codes FILE          file containing ip codes not in the minimizers" << endl
         << "  -Z, --zip-tree                create a zipcode tree instead of clustering" << endl
         << "computational parameters:" << endl 
         << "  -t, --threads INT             number of compute threads to use" << endl;
}

int main_cluster(int argc, char** argv) {

    if (argc == 2) {
        help_cluster(argv);
        return 1;
    }

    // initialize parameters with their default options
    string graph_name = "";
    bool gbz_format = false;
    string gcsa_name;
    string zipcode_name;
    string minimizer_name = "";
    bool long_reads = false;
    string distance_name;
    // How close should two hits be to be in the same cluster?
    size_t distance_limit = 1000;
    size_t hit_cap = 10;
    size_t hard_hit_cap = 500;
    size_t hits_above_threshold = std::numeric_limits<size_t>::max();
    float score_fraction = 0.9;
    size_t max_min = 500;
    size_t num_bp_per_min = 1000;
    size_t downsample_min = 0;
    bool output_sequences_only = false;
    bool make_zip_tree = false;
 
    //Get an index registry to keep track of all the indexes
    IndexRegistry registry = VGIndexes::get_vg_index_registry();

   
    int c;
    optind = 2; // force optind past command positional argument
    while (true) {
        static struct option long_options[] =
        {
            {"help", no_argument, 0, 'h'},
            {"xg-name", required_argument, 0, 'x'},
            {"gbz-format", no_argument, 0, 'f'},
            {"gcsa-name", required_argument, 0, 'g'},
            {"gbwt-name", required_argument, 0, 'G'},
            {"gbwtgraph-name", required_argument, 0, 'B'},
            {"minimizer-name", required_argument, 0, 'm'},
            {"long-reads", no_argument, 0, 'l'},
            {"dist-name", required_argument, 0, 'd'},
            {"prefix", required_argument, 0, 'p'},
            {"hit-cap", required_argument, 0, 'c'},
            {"hard-hit-cap", required_argument, 0, 'C'},
            {"hits-above", required_argument, 0, 'a'},
            {"sequences-only", no_argument, 0, 'S'},
            {"score-fraction", required_argument, 0, 'F'},
            {"max-min", required_argument, 0, 'U'},
            {"num-bp-per-min", required_argument, 0, 'b'},
            {"downsample-min", required_argument, 0, 'D'},
            {"zip-codes", required_argument, 0, 'z'},
            {"zip-tree", no_argument, 0, 'Z'},
            {"threads", required_argument, 0, 't'},
            {0, 0, 0, 0}
        };

        int option_index = 0;
        c = getopt_long (argc, argv, "h?x:fg:G:B:m:ld:p:c:C:a:SF:U:b:D:z:Zt:",
                         long_options, &option_index);


        // Detect the end of the options.
        if (c == -1)
            break;

        switch (c)
        {
            case 'x':
                // Remember the string for MEMs
                graph_name = require_exists(context, optarg);
                break;

            case 'f':
                gbz_format = true;
                break;
                
            case 'g':
                registry.provide("GCSA", optarg);
                break;
            
            case 'G':
                registry.provide("Giraffe GBWT", optarg);
                break;
                

            case 'B':
                registry.provide("GBWTGraph", optarg);
                
                // But if we have a GBWTGraph we probably want to use *its* name as the base name.
                // Whichever is specified last will win, unless we also have a FASTA input name.
                registry.set_prefix(split_ext(optarg).first);
                
                break;

            
            case 'm':
                minimizer_name = require_exists(context, optarg);
                break;
                
            case 'l':
                long_reads = true;
                break;
                
            case 'd':
                distance_name = require_exists(context, optarg);
                registry.provide("Giraffe Distance Index", optarg);
                break;

            case 'p':
                if (!optarg || !*optarg) {
                    fatal_error(context) << "Must provide prefix with -p." << endl;
                }
                registry.set_prefix(optarg);
                break;
            
            case 'c':
                hit_cap = parse<size_t>(optarg);
                break;

            case 'C':
                hard_hit_cap = parse<size_t>(optarg);
                break;

            case 'a':
                hits_above_threshold = parse<size_t>(optarg);
                break;
            
            case 'S':
                output_sequences_only = true;
                break;

            case 'F':
                score_fraction = parse<double>(optarg);
                break;

            case 'U':
                max_min = parse<size_t>(optarg);
                break;

            case 'b':
                num_bp_per_min = parse<size_t>(optarg);
                break;

            case 'D':
                downsample_min = parse<size_t>(optarg);
                break;
            
            case 'z':
                zipcode_name = optarg;
                break;

            case 'Z':
                make_zip_tree = true;
                break;
                
            case 't':
                omp_set_num_threads(parse_thread_count(context, optarg));
                break;
                
            case 'h':
            case '?':
            default:
                help_cluster(argv);
                exit(1);
                break;
        }
    }
    
    if (graph_name.empty()) {
        fatal_error(context) << "Must provide a graph file with -x." << endl;
    }

    if (hits_above_threshold != std::numeric_limits<size_t>::max()) {
        if (hard_hit_cap < hits_above_threshold) {
            fatal_error(context) << "Hard hit cap (-C) must be greater than or equal "
                                 << "to hits-above threshold (-a)." << endl;
        }
    } else {
        if (output_sequences_only) {
            fatal_error(context) << "Cannot use -S (sequences only) without a hits-above threshold (-a)." << endl;
        }
    }

    // Give the file to the index registry for clustering minimizers
    if (gbz_format) {
        registry.provide("GBZ", graph_name);
    } else {
        registry.provide("XG", graph_name);
    }

    std::string minimizer_index_type = long_reads ? "Long Read Minimizers" : "Short Read Minimizers";

    if (!minimizer_name.empty()) {
        registry.provide(minimizer_index_type, minimizer_name);
    }

    // We define a child class to expose protected stuff
    // This is copied from the minimizer mapper unit tests
    class TestMinimizerMapper : public MinimizerMapper {
    public:
        TestMinimizerMapper(
            gbwtgraph::GBWTGraph gbwt_graph,
            gbwtgraph::DefaultMinimizerIndex minimizer_index,
            SnarlDistanceIndex* distance_index,
            PathPositionHandleGraph* handle_graph)
            : MinimizerMapper(gbwt_graph, minimizer_index, distance_index, nullptr, handle_graph){};
        using MinimizerMapper::MinimizerMapper;
        using MinimizerMapper::Minimizer;
        using MinimizerMapper::find_minimizers;
        using MinimizerMapper::sort_minimizers_by_score;
        using MinimizerMapper::find_seeds;
        using MinimizerMapper::hit_cap;
        using MinimizerMapper::hard_hit_cap;
        using MinimizerMapper::minimizer_score_fraction;
        using MinimizerMapper::max_unique_min;
        using MinimizerMapper::num_bp_per_min;
        using MinimizerMapper::minimizer_downsampling_window_count;
        using MinimizerMapper::track_provenance;

    };
    
    // create in-memory objects for mems
    unique_ptr<PathHandleGraph> path_handle_graph;
    bdsg::PathPositionOverlayHelper overlay_helper;
    PathPositionHandleGraph* xg_index;
    unique_ptr<gcsa::GCSA> gcsa_index;
    unique_ptr<gcsa::LCPArray> lcp_index;

    // The IndexRegistry doesn't try to infer index files based on the
    // basename, so do that here. We can have multiple extension options that
    // we try in order of priority.
    unordered_map<string, vector<string>> indexes_and_extensions = {
        {"Giraffe GBZ", {"giraffe.gbz", "gbz"}},
        {"Giraffe Distance Index", {"dist"}},
        {minimizer_index_type, {long_reads ? "longread.withzip.min" : "shortread.withzip.min"}},
        {long_reads ? "Long Read Zipcodes" : "Short Read Zipcodes", 
            {long_reads ? "longread.zipcodes" : "shortread.zipcodes"}}
    };
    //Get minimizer indexes
    for (auto& completed : registry.completed_indexes()) {
        // Drop anything we already got from the list
        indexes_and_extensions.erase(completed);
    }
    for (auto& index_and_extensions : indexes_and_extensions) {
        // For each index type
        for (auto& extension : index_and_extensions.second) {
            // For each extension in priority order
            string inferred_filename = registry.get_prefix() + "." + extension;
            if (ifstream(inferred_filename).is_open()) {
                // A file with the appropriate name exists and we can read it
                registry.provide(index_and_extensions.first, inferred_filename);
                // Report it because this may not be desired behavior
                cerr << "Guessing that " << inferred_filename << " is " << index_and_extensions.first << endl;
                // Skip other extension options for the index
                break;
            }
        }
    }
    // create in-memory objects

    // Don't try and use all the memory.
    // TODO: add memory options like autoindex?
    registry.set_target_memory_usage(IndexRegistry::get_system_memory() / 2);

    auto index_targets = long_reads 
                         ? VGIndexes::get_default_long_giraffe_indexes() 
                         : VGIndexes::get_default_short_giraffe_indexes();

    //Make sure we have all necessary indexes
    try {
        registry.make_indexes(index_targets);
    }
    catch (InsufficientInputException ex) {
        fatal_error(context) << "Input is not sufficient to create indexes:\n" << ex.what() << endl;
    }

    //Get the minimizer index
    auto minimizer_index = vg::io::VPKG::load_one<gbwtgraph::DefaultMinimizerIndex>(registry.require(minimizer_index_type).at(0));

    //Get the zipcodes
    ZipCodeCollection oversized_zipcodes;
    if (!zipcode_name.empty()) {

        ifstream zip_in (zipcode_name);
        oversized_zipcodes.deserialize(zip_in);
        zip_in.close();
    }

    // Grab the GBZ
    auto gbz = vg::io::VPKG::load_one<gbwtgraph::GBZ>(registry.require("Giraffe GBZ").at(0));

    //Get the distance index
    auto distance_index = vg::io::VPKG::load_one<SnarlDistanceIndex>(registry.require("Giraffe Distance Index").at(0));



    // Make the clusterer
    SnarlDistanceIndexClusterer clusterer(*distance_index);

    
    // Make a Mapper to look up MEM seeds
    unique_ptr<Mapper> mapper;
    if (gcsa_index) {
        // We will find MEMs using a Mapper
        mapper = make_unique<Mapper>(xg_index, gcsa_index.get(), lcp_index.get());
        if (output_sequences_only) {
            fatal_error(context) << "Cannot output minimizers (-S) with a GCSA index (-g)." << endl;
        }
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
            vector<pos_t> positions;


            //Make a vector of seeds for using minimizer to cluster
            vector<SnarlDistanceIndexClusterer::Seed> seeds;
            
            // If working with MEMs, this will hold all the MEMs
            vector<MaximalExactMatch> mems;
            // If working with minimizers, this will hold all the minimizers in the query
            vector<MinimizerMapper::Minimizer> minimizers_in_read;
            // And either way this will map from seed to MEM or minimizer that generated it
            vector<size_t> seed_to_source;
            VectorView<MinimizerMapper::Minimizer> minimizers;
            std::vector<size_t> minimizer_score_order;
            
            if (mapper) {
                // Find MEMs
                double lcp_avg, fraction_filtered;
                mems = mapper->find_mems_deep(aln.sequence().begin(), aln.sequence().end(), lcp_avg, fraction_filtered);
                
                // Convert to position seeds
                for (size_t i = 0; i < mems.size(); i++) {
                    auto& mem = mems[i];
                    for (gcsa::node_type n : mem.nodes) {
                        // Convert from GCSA node_type packing to a pos_t
                        positions.push_back(make_pos_t(n));
                        // And remember which MEM the seed came from.
                        seed_to_source.push_back(i);
                    }
                }
            } else {
                // Find minimizers
                assert(minimizer_index);

                //Use a MinimizerMapper to find the minimizers, using the provided parameters
                //This will have an empty gbwtgraph::GBWTGraph, so it shouldn't be used
                //for anything except finding minimizers
                TestMinimizerMapper minimizer_mapper(gbz->graph, *minimizer_index, &(*distance_index), &oversized_zipcodes, nullptr);

                //Set the parameters
                minimizer_mapper.hit_cap = hit_cap;
                minimizer_mapper.hard_hit_cap = hard_hit_cap;
                minimizer_mapper.minimizer_score_fraction = score_fraction;
                minimizer_mapper.max_unique_min = max_min;
                minimizer_mapper.num_bp_per_min = num_bp_per_min;
                minimizer_mapper.minimizer_downsampling_window_count = downsample_min;
                minimizer_mapper.track_provenance = true;
                Funnel funnel;
                funnel.start(aln.name());

                //Find the minimizers and then the seeds using the minimizer mapper
                minimizers_in_read = minimizer_mapper.find_minimizers(aln.sequence(), funnel);
                if (hits_above_threshold < hard_hit_cap) {
                    // minimizer : (hits in read, hits in graph)
                    unordered_map<std::string, std::pair<size_t, size_t>> minimizer_hit_counts;
                    std::string cur_minimizer;
                    for (auto& m : minimizers_in_read) {
                        cur_minimizer = m.value.key.decode(m.length);
                        if (minimizer_hit_counts.find(cur_minimizer) == minimizer_hit_counts.end()) {
                            // If we haven't seen this minimizer before, initialize its count
                            minimizer_hit_counts[cur_minimizer] = make_pair(1, m.hits);
                        }
                        minimizer_hit_counts[cur_minimizer].first++;
                    }

                    std::stringstream high_hit_minimizers;
                    for (auto& m : minimizer_hit_counts) {
                        // For each minimizer, if it has more hits than the threshold, add it to the annotation
                        if (m.second.second * m.second.first > hits_above_threshold) {
                            high_hit_minimizers << m.first << "=" << m.second.first << "x" 
                                                << m.second.second << ",";
                        }
                    }
                    set_annotation(aln, "high_hit_minimizers", high_hit_minimizers.str());
                }

                if (output_sequences_only) {
                    // Skip finding seeds, we only care about minimizers
                    #pragma omp critical (cout)
                    emitter.write(std::move(aln));
                    
                    funnel.stop();
                } else {
                    // Indexes of minimizers, sorted into score order, best score first
                    LazyRNG rng([&]() {
                        return aln.sequence();
                    });
                    minimizer_score_order = minimizer_mapper.sort_minimizers_by_score(minimizers_in_read, rng);

                    // Minimizers sorted by best score first
                    minimizers = {minimizers_in_read, minimizer_score_order};
                    
                    // Find the seeds and mark the minimizers that were located.
                    seeds = minimizer_mapper.find_seeds(minimizers_in_read, minimizers, aln, funnel);

                    //Fill in seeds_to_source using the funnel
                    vector<vector<size_t>> seed_to_source_vector = funnel.map_stage_results_to_previous_stage("seed");

                    //This was a vector of vectors, but each seed came from just one minimizer, so flatten the vector
                    for (auto& v : seed_to_source_vector) {
                        assert(v.size() == 1);
                        seed_to_source.emplace_back(v.front());
                    }
                    assert(seed_to_source.size() == seeds.size());
                    funnel.stop();
                }
            }

            
                
            if (!output_sequences_only) {
                if (make_zip_tree) {
                    //Time making the zipcode tree

                    ZipCodeForest zip_forest;

                    std::chrono::time_point<std::chrono::system_clock> start = std::chrono::system_clock::now();
                    zip_forest.fill_in_forest(seeds, minimizers, *distance_index, std::numeric_limits<size_t>::max());
                    std::chrono::time_point<std::chrono::system_clock> end = std::chrono::system_clock::now();
                    std::chrono::duration<double> elapsed_seconds = end-start;

                    std::pair<size_t, size_t> dag_non_dag_count (0, 0);
                    for (const auto& zip_tree : zip_forest.trees) {
                        pair<size_t, size_t> tree_count = zip_tree.dag_and_non_dag_snarl_count(seeds, *distance_index);
                        dag_non_dag_count.first += tree_count.first;
                        dag_non_dag_count.second += tree_count.second;
                    }

                    std::string cyclic_snarl_sizes_string;
                    for (const auto& zip_tree : zip_forest.trees) {
                        for (const auto& size : zip_tree.cyclic_snarl_sizes(seeds, *distance_index)) {
                            cyclic_snarl_sizes_string += std::to_string(size) + ",";
                        }
                    }

                    // And with hit count clustered
                    set_annotation(aln, "seed_count", (double)seeds.size());

                    // Annotate with the time spent making the zip tree
                    set_annotation(aln, "zip_tree_construction_seconds", elapsed_seconds.count());

                    //The number of snarls that are dags
                    set_annotation(aln, "zip_tree_dag_count", dag_non_dag_count.first);

                    //The number of snarls that aren't dags
                    set_annotation(aln, "zip_tree_non_dag_count", dag_non_dag_count.second);

                    //CSV of the sizes of the cyclic snarls
                    set_annotation(aln, "zip_tree_cyclic_snarl_sizes", cyclic_snarl_sizes_string);

                    // TODO: parallelize this
                    #pragma omp critical (cout)
                    emitter.write(std::move(aln));

                } else {
                    // Cluster the seeds. Get sets of input seed indexes that go together.
                    // Make sure to time it.
                    std::chrono::time_point<std::chrono::system_clock> start = std::chrono::system_clock::now();
                    vector<SnarlDistanceIndexClusterer::Cluster> clusters = clusterer.cluster_seeds(seeds,  distance_limit);
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
                                size_t start_offset = minimizers_in_read[source_index].forward_offset();
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
                    std::sort(cluster_indexes_in_order.begin(), cluster_indexes_in_order.end(), 
                    [&](const size_t& a, const size_t& b) -> bool {
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
                                best.push_back(positions.at(seed_index));
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
                    for (auto& pos : positions) {
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
                    set_annotation(aln, "seed_count", (double)positions.size());
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
                }
            }
        });
    });
    
    return 0;
}

// Register subcommand
static Subcommand vg_cluster("cluster", "find and cluster mapping seeds", DEVELOPMENT, main_cluster);


