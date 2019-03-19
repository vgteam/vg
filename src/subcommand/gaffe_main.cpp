/**
 * \file gaffe_main.cpp: GFA (Graph Alignment Format) Fast Emitter: a new mapper that will be *extremely* fast once we actually write it
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

#include "../seed_clusterer.hpp"
#include "../mapper.hpp"
#include "../annotation.hpp"
#include "../minimizer.hpp"
#include "../stream/vpkg.hpp"
#include "../stream/stream.hpp"
#include "../alignment_emitter.hpp"
#include "../gapless_extender.hpp"

using namespace std;
using namespace vg;
using namespace vg::subcommand;

void help_gaffe(char** argv) {
    cerr
    << "usage: " << argv[0] << " gaffe [options] > output.gam" << endl
    << "Map unpaired reads using minimizers and gapless extension." << endl
    << endl
    << "basic options:" << endl
    << "  -x, --xg-name FILE            use this xg index (required)" << endl
    << "  -H, --gbwt-name FILE          use this GBWT index (required)" << endl
    << "  -m, --minimizer-name FILE     use this minimizer index (required)" << endl
    << "  -s, --snarls FILE             cluster using these snarls (required)" << endl
    << "  -d, --dist-name FILE          cluster using this distance index (required)" << endl
    << "  -c, --hit-cap INT             ignore minimizers with more than this many locations [10]" << endl
    << "input options:" << endl
    << "  -G, --gam-in FILE             read and realign GAM-format reads from FILE (may repeat)" << endl
    << "  -f, --fastq-in FILE           read and align FASTQ-format reads from FILE (may repeat)" << endl
    << "output options:" << endl
    << "  -M, --max-multimaps INT       produce up to INT alignments for each read [1]"
    << "  -N, --sample NAME             add this sample name" << endl
    << "  -R, --read-group NAME         add this read group" << endl
    << "computational parameters:" << endl
    << "  -t, --threads INT             number of compute threads to use" << endl;
}

int main_gaffe(int argc, char** argv) {

    if (argc == 2) {
        help_gaffe(argv);
        return 1;
    }

    // initialize parameters with their default options
    string xg_name;
    string gbwt_name;
    string minimizer_name;
    string snarls_name;
    string distance_name;
    // How close should two hits be to be in the same cluster?
    size_t distance_limit = 1000;
    size_t hit_cap = 10;
    // What GAMs should we realign?
    vector<string> gam_filenames;
    // What FASTQs should we align.
    // Note: multiple FASTQs are not interpreted as paired.
    vector<string> fastq_filenames;
    // How many mappigns per read can we emit?
    size_t max_multimaps = 1;
    // What sample name if any should we apply?
    string sample_name;
    // What read group if any should we apply?
    string read_group;
    
    int c;
    optind = 2; // force optind past command positional argument
    while (true) {
        static struct option long_options[] =
        {
            {"help", no_argument, 0, 'h'},
            {"xg-name", required_argument, 0, 'x'},
            {"gbwt-name", required_argument, 0, 'H'},
            {"minimizer-name", required_argument, 0, 'm'},
            {"snarls", required_argument, 0, 's'},
            {"dist-name", required_argument, 0, 'd'},
            {"hit-cap", required_argument, 0, 'c'},
            {"gam-in", required_argument, 0, 'G'},
            {"fastq-in", required_argument, 0, 'f'},
            {"max-multimaps", required_argument, 0, 'M'},
            {"sample", required_argument, 0, 'N'},
            {"read-group", required_argument, 0, 'R'},
            {"threads", required_argument, 0, 't'},
            {0, 0, 0, 0}
        };

        int option_index = 0;
        c = getopt_long (argc, argv, "hx:H:m:s:d:c:G:f:M:t:",
                         long_options, &option_index);


        // Detect the end of the options.
        if (c == -1)
            break;

        switch (c)
        {
            case 'x':
                xg_name = optarg;
                if (xg_name.empty()) {
                    cerr << "error:[vg gaffe] Must provide XG file with -x." << endl;
                    exit(1);
                }
                break;
                
            case 'H':
                gbwt_name = optarg;
                if (gbwt_name.empty()) {
                    cerr << "error:[vg gaffe] Must provide GBWT file with -H." << endl;
                    exit(1);
                }
                break;
                
            case 'm':
                minimizer_name = optarg;
                if (minimizer_name.empty()) {
                    cerr << "error:[vg gaffe] Must provide minimizer file with -m." << endl;
                    exit(1);
                }
                break;
                
            case 's':
                snarls_name = optarg;
                if (snarls_name.empty()) {
                    cerr << "error:[vg gaffe] Must provide snarl file with -s." << endl;
                    exit(1);
                }
                break;
                
            case 'd':
                distance_name = optarg;
                if (distance_name.empty()) {
                    cerr << "error:[vg gaffe] Must provide distance index file with -d." << endl;
                    exit(1);
                }
                break;
            
            case 'c':
                hit_cap = parse<size_t>(optarg);
                break;
                
            case 'G':
                gam_filenames.push_back(optarg);
                break;
            
            case 'f':
                fastq_filenames.push_back(optarg);
                break;
                
            case 'M':
                max_multimaps = parse<size_t>(optarg);
                break;
            
            case 'N':
                sample_name = optarg;
                break;
                
            case 'R':
                read_group = optarg;
                break;
                
            case 't':
            {
                int num_threads = parse<int>(optarg);
                if (num_threads <= 0) {
                    cerr << "error:[vg gaffe] Thread count (-t) set to " << num_threads << ", must set to a positive integer." << endl;
                    exit(1);
                }
                omp_set_num_threads(num_threads);
            }
                break;
                
            case 'h':
            case '?':
            default:
                help_gaffe(argv);
                exit(1);
                break;
        }
    }
    
    
    if (xg_name.empty()) {
        cerr << "error:[vg gaffe] Mapping requires an XG index (-x)" << endl;
        exit(1);
    }
    
    if (gbwt_name.empty()) {
        cerr << "error:[vg gaffe] Mapping requires a GBWT index (-H)" << endl;
        exit(1);
    }
    
    if (minimizer_name.empty()) {
        cerr << "error:[vg gaffe] Mapping requires a minimizer index (-m)" << endl;
        exit(1);
    }
    
    if (snarls_name.empty()) {
        cerr << "error:[vg gaffe] Mapping requires snarls (-s)" << endl;
        exit(1);
    }
    
    if (distance_name.empty()) {
        cerr << "error:[vg gaffe] Mapping requires a distance index (-d)" << endl;
        exit(1);
    }
    
    // create in-memory objects
    unique_ptr<xg::XG> xg_index = stream::VPKG::load_one<xg::XG>(xg_name);
    unique_ptr<gbwt::GBWT> gbwt_index = stream::VPKG::load_one<gbwt::GBWT>(gbwt_name);
    unique_ptr<MinimizerIndex> minimizer_index = stream::VPKG::load_one<MinimizerIndex>(minimizer_name);
    unique_ptr<SnarlManager> snarl_manager = stream::VPKG::load_one<SnarlManager>(snarls_name);
    unique_ptr<DistanceIndex> distance_index = stream::VPKG::load_one<DistanceIndex>(distance_name);
    
    // Connect the DistanceIndex to the other things it needs to work.
    distance_index->setGraph(xg_index.get());
    distance_index->setSnarlManager(snarl_manager.get());
    
    // Make a GBWTGraph over the GBWT and the XG
    GBWTGraph gbwt_graph(*gbwt_index, *xg_index);
    
    // Make a gapless extender to extend seed hits in haplotype space.
    GaplessExtender extender(gbwt_graph);
    
    // Make the clusterer
    SnarlSeedClusterer clusterer;
    
    // Set up output to an emitter that will handle serialization
    unique_ptr<AlignmentEmitter> alignment_emitter = get_alignment_emitter("-", "GAM", {});
    
    // Define how to align and output a read, in a thread.
    auto map_read = [&](Alignment& aln) {
        // For each input alignment
            
        // We will find all the seed hits
        vector<pos_t> seeds;
        
        // This will hold all the minimizers in the query
        vector<MinimizerIndex::minimizer_type> minimizers;
        // And either way this will map from seed to minimizer that generated it
        vector<size_t> seed_to_source;
        
        // Find minimizers in the query
        minimizers = minimizer_index->minimizers(aln.sequence());
        
        for (size_t i = 0; i < minimizers.size(); i++) {
            // For each minimizer
            if (hit_cap != 0 && minimizer_index->count(minimizers[i].first) <= hit_cap) {
                // The minimizer is infrequent enough to be informative, so feed it into clustering
                
                // Locate it in the graph
                for (auto& hit : minimizer_index->find(minimizers[i].first)) {
                    // For each position, remember it and what minimizer it came from
                    seeds.push_back(hit);
                    seed_to_source.push_back(i);
                }
            }
        }
            
        // Cluster the seeds. Get sets of input seed indexes that go together.
        vector<hash_set<size_t>> clusters = clusterer.cluster_seeds(seeds, distance_limit, *snarl_manager, *distance_index);
        
        // Compute the covered portion of the read represented by each cluster.
        // TODO: Put this and sorting into the clusterer to deduplicate with vg cluster.
        vector<double> read_coverage_by_cluster;
        for (auto& cluster : clusters) {
            // We set bits in here to true when query anchors cover them
            vector<bool> covered(aln.sequence().size());
            // We use this to convert iterators to indexes
            auto start = aln.sequence().begin();
            
            for (auto& hit_index : cluster) {
                // For each hit in the cluster, work out what anchor sequence it is from.
                size_t source_index = seed_to_source.at(hit_index);
                
                for (size_t i = minimizers[source_index].second; i < minimizers[source_index].second + minimizer_index->k(); i++) {
                    // Set all the bits in read space for that minimizer.
                    // Each minimizr is a length-k exact match starting at a position
                    covered[i] = true;
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
        
        // We will fill this with the output alignments (primary and secondaries) in score order.
        vector<Alignment> aligned;
        aligned.reserve(cluster_indexes_in_order.size());
        
        // Annotate the original read with metadata before copying
        if (!sample_name.empty()) {
            aln.set_sample_name(sample_name);
        }
        if (!read_group.empty()) {
            aln.set_read_group(read_group);
        }
        
        // How many multimaps will we produce?
        size_t actual_multimaps = min(cluster_indexes_in_order.size(), max_multimaps);
        
        for (size_t i = 0; i < actual_multimaps; i++) {
            // For each cluster
            hash_set<size_t>& cluster = clusters[i];
            
            // Pack the seeds into (read position, graph position) pairs.
            vector<pair<size_t, pos_t>> seed_matchings;
            seed_matchings.reserve(cluster.size());
            for (auto& seed_index : cluster) {
                // For each seed in the cluster, generate its matching pair
                seed_matchings.emplace_back(minimizers[seed_to_source[seed_index]].second, seeds[seed_index]);
            }
            
            // Extend seed hits in the cluster into a real alignment path and mismatch count.
            std::pair<Path, size_t> extended = extender.extend_seeds(seed_matchings, aln.sequence());
            auto& path = extended.first;
            auto& mismatch_count = extended.second;
            
            // Produce an output Alignment
            aligned.emplace_back(aln);
            Alignment& out = aligned.back();
            
            if (path.mapping_size() != 0) {
                // We have a mapping
                
                // Compute a score based on the sequence length and mismatch count.
                // Alignments will only contain matches and mismatches.
                int alignment_score = default_match * (aln.sequence().size() - mismatch_count) - default_mismatch * extended.second;
                // Compute identity from mismatch count.
                double identity = aln.sequence().size() == 0 ? 0.0 : (aln.sequence().size() - mismatch_count) / (double) aln.sequence().size();
                
                // Fill in the extension info
                *out.mutable_path() = path;
                out.set_score(alignment_score);
                out.set_identity(identity);
                
            } else {
                // Read was not able to be mapped; mismatch_count is not true.
                // Make an un-aligned alignment.
                
                out.clear_path();
                out.set_score(0);
                out.set_identity(0);
            }
            
            // Set other fields
            out.clear_refpos();
        }
        
        // Sort again by actual score instead of cluster coverage
        std::sort(aligned.begin(), aligned.end(), [](const Alignment& a, const Alignment& b) -> bool {
            // Return true if a must come before b (i.e. it has a larger score)
            return a.score() > b.score();
        });
        
        for (size_t i = 0; i < aligned.size(); i++) {
            // For each output alignment in score order
            auto& out = aligned[i];
            
            // Assign primary and secondary status
            out.set_is_secondary(i > 0);
            out.set_mapping_quality(0);
        }
        
        // Ship out all the aligned alignments
        alignment_emitter->emit_mapped_single(std::move(aligned));
    };
    
    for (auto& gam_name : gam_filenames) {
        // For every GAM file to remap
        get_input_file(gam_name, [&](istream& in) {
            // Open it and map all the reads in parallel.
            stream::for_each_parallel<Alignment>(in, map_read);
        });
    }
    
    for (auto& fastq_name : fastq_filenames) {
        // For every FASTQ file to map, map all its reads in parallel.
        fastq_unpaired_for_each_parallel(fastq_name, map_read);
    }
        
    return 0;
}

// Register subcommand
static Subcommand vg_gaffe("gaffe", "Graph Alignment Format Fast Emitter", DEVELOPMENT, main_gaffe);


