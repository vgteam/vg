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
#include <mutex>

#include "subcommand.hpp"

#include "../seed_clusterer.hpp"
#include "../mapper.hpp"
#include "../annotation.hpp"
#include "../minimizer.hpp"
#include <vg/io/vpkg.hpp>
#include <vg/io/stream.hpp>
#include "../alignment_emitter.hpp"
#include "../gapless_extender.hpp"
#include "../minimizer_mapper.hpp"

//#define USE_CALLGRIND

#ifdef USE_CALLGRIND
#include <valgrind/callgrind.h>
#endif

using namespace std;
using namespace vg;
using namespace vg::subcommand;

void help_gaffe(char** argv) {
    cerr
    << "usage: " << argv[0] << " gaffe [options] > output.gam" << endl
    << "Map unpaired reads using minimizers and gapless extension." << endl
    << endl
    << "basic options:" << endl
    << "  -x, --xg-name FILE            use this xg index (required if -g not specified)" << endl
    << "  -g, --graph-name FILE         use this GBWTGraph (required if -x not specified)" << endl
    << "  -H, --gbwt-name FILE          use this GBWT index (required)" << endl
    << "  -m, --minimizer-name FILE     use this minimizer index (required)" << endl
    << "  -d, --dist-name FILE          cluster using this distance index (required)" << endl
    << "  -p, --progress                show progress" << endl
    << "input options:" << endl
    << "  -G, --gam-in FILE             read and realign GAM-format reads from FILE (may repeat)" << endl
    << "  -f, --fastq-in FILE           read and align FASTQ-format reads from FILE (may repeat)" << endl
    << "output options:" << endl
    << "  -M, --max-multimaps INT       produce up to INT alignments for each read [1]" << endl
    << "  -N, --sample NAME             add this sample name" << endl
    << "  -R, --read-group NAME         add this read group" << endl
    << "  -n, --discard                 discard all output alignments (for profiling)" << endl
    << "computational parameters:" << endl
    << "  -c, --hit-cap INT             use all minimizers with at most INT hits [10]" << endl
    << "  -C, --hard-hit-cap INT        ignore all minimizers with more than INT hits [300]" << endl
    << "  -F, --score-fraction FLOAT    select minimizers between hit caps until score is FLOAT of total [0.6]" << endl
    << "  -e, --max-extensions INT      extend up to INT clusters [48]" << endl
    << "  -a, --max-alignments INT      align up to INT extensions [8]" << endl
    << "  -s, --cluster-score INT       only extend clusters if they are within cluster-score of the best score" << endl
    << "  -u, --cluster-coverage FLOAT  only extend clusters if they are within cluster-coverage of the best read coverage" << endl
    << "  -v, --extension-score INT     only align extensions if their score is within extension-score of the best score [1]" << endl
    << "  -w, --extension-set INT       only align extension sets if their score is within extension-set of the best score" << endl
    << "  -O, --no-dp                   disable all gapped alignment" << endl
    << "  --track-provenance            track how internal intermediate alignment candidates were arrived at" << endl
    << "  --track-correctness           track if internal intermediate alignment candidates are correct (implies --track-provenance)" << endl
    << "  -t, --threads INT             number of compute threads to use" << endl;
}

int main_gaffe(int argc, char** argv) {

    std::chrono::time_point<std::chrono::system_clock> launch = std::chrono::system_clock::now();

    if (argc == 2) {
        help_gaffe(argv);
        return 1;
    }

    #define OPT_TRACK_PROVENANCE 1001
    #define OPT_TRACK_CORRECTNESS 1002

    // initialize parameters with their default options
    string xg_name;
    string graph_name;
    string gbwt_name;
    string minimizer_name;
    string distance_name;
    // How close should two hits be to be in the same cluster?
    size_t distance_limit = 1000;
    size_t hit_cap = 10, hard_hit_cap = 300;
    double minimizer_score_fraction = 0.6;
    bool progress = false;
    // Should we try chaining or just give up if we can't find a full length gapless alignment?
    bool do_dp = true;
    // What GAMs should we realign?
    vector<string> gam_filenames;
    // What FASTQs should we align.
    // Note: multiple FASTQs are not interpreted as paired.
    vector<string> fastq_filenames;
    // How many mappings per read can we emit?
    size_t max_multimaps = 1;
    // How many clusters should we extend?
    size_t max_extensions = 48;
    // How many extended clusters should we align, max?
    size_t max_alignments = 8;
    //Throw away cluster with scores that are this amount below the best
    double cluster_score = 0;
    //Throw away clusters with coverage this amount below the best 
    double cluster_coverage = 0;
    //Throw away extension sets with scores that are this amount below the best
    double extension_set = 0;
    //Throw away extensions with scores that are this amount below the best
    int extension_score = 1;
    // What sample name if any should we apply?
    string sample_name;
    // What read group if any should we apply?
    string read_group;
    // Should we throw out our alignments instead of outputting them?
    bool discard_alignments = false;
    // Should we track candidate provenance?
    bool track_provenance = false;
    // Should we track candidate correctness?
    bool track_correctness = false;
    
    vector<size_t> threads_to_run;
    
    int c;
    optind = 2; // force optind past command positional argument
    while (true) {
        static struct option long_options[] =
        {
            {"help", no_argument, 0, 'h'},
            {"xg-name", required_argument, 0, 'x'},
            {"graph-name", required_argument, 0, 'g'},
            {"gbwt-name", required_argument, 0, 'H'},
            {"minimizer-name", required_argument, 0, 'm'},
            {"dist-name", required_argument, 0, 'd'},
            {"progress", no_argument, 0, 'p'},
            {"gam-in", required_argument, 0, 'G'},
            {"fastq-in", required_argument, 0, 'f'},
            {"max-multimaps", required_argument, 0, 'M'},
            {"sample", required_argument, 0, 'N'},
            {"read-group", required_argument, 0, 'R'},
            {"discard", no_argument, 0, 'n'},
            {"hit-cap", required_argument, 0, 'c'},
            {"hard-hit-cap", required_argument, 0, 'C'},
            {"max-extensions", required_argument, 0, 'e'},
            {"max-alignments", required_argument, 0, 'a'},
            {"cluster-score", required_argument, 0, 's'},
            {"cluster-coverage", required_argument, 0, 'u'},
            {"extension-score", required_argument, 0, 'v'},
            {"extension-set", required_argument, 0, 'w'},
            {"score-fraction", required_argument, 0, 'F'},
            {"no-dp", no_argument, 0, 'O'},
            {"track-provenance", no_argument, 0, OPT_TRACK_PROVENANCE},
            {"track-correctness", no_argument, 0, OPT_TRACK_CORRECTNESS},
            {"threads", required_argument, 0, 't'},
            {0, 0, 0, 0}
        };

        int option_index = 0;
        c = getopt_long (argc, argv, "hx:g:H:m:s:d:pG:f:M:N:R:nc:C:F:e:a:s:u:v:w:Ot:",
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

            case 'g':
                graph_name = optarg;
                if (graph_name.empty()) {
                    cerr << "error:[vg gaffe] Must provide GBGTGraph file with -g." << endl;
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
                
            case 'd':
                distance_name = optarg;
                if (distance_name.empty()) {
                    cerr << "error:[vg gaffe] Must provide distance index file with -d." << endl;
                    exit(1);
                }
                break;

            case 'p':
                progress = true;
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
                
            case 'n':
                discard_alignments = true;
                break;

            case 'c':
                {
                    size_t cap = parse<size_t>(optarg);
                    if (cap <= 0) {
                        cerr << "error: [vg gaffe] Hit cap (" << cap << ") must be a positive integer" << endl;
                        exit(1);
                    }
                    hit_cap = cap;
                }
                break;

            case 'C':
                {
                    size_t cap = parse<size_t>(optarg);
                    if (cap <= 0) {
                        cerr << "error: [vg gaffe] Hard hit cap (" << cap << ") must be a positive integer" << endl;
                        exit(1);
                    }
                    hard_hit_cap = cap;
                }
                break;

            case 'F':
                minimizer_score_fraction = parse<double>(optarg);
                break;

            case 'e':
                {
                    size_t extensions = parse<size_t>(optarg);
                    if (extensions <= 0) {
                        cerr << "error: [vg gaffe] Number of extensions (" << extensions << ") must be a positive integer" << endl;
                        exit(1);
                    }
                    max_extensions = extensions;
                }
                break;

            case 'a':
                {
                    size_t alignments = parse<size_t>(optarg);
                    if (alignments <= 0) {
                        cerr << "error: [vg gaffe] Number of alignments (" << alignments << ") must be a positive integer" << endl;
                        exit(1);
                    }
                    max_alignments = alignments;
                }
                break;

            case 's':
                {
                    double score = parse<double>(optarg);
                    if (score < 0) {
                        cerr << "error: [vg gaffe] Cluster score threshold (" << score << ") must be positive" << endl;
                        exit(1);
                    }
                    cluster_score = score;
                }
                break;

            case 'u':
                {
                    double score = parse<double>(optarg);
                    if (score < 0) {
                        cerr << "error: [vg gaffe] Cluster coverage threshold (" << score << ") must be positive" << endl;
                        exit(1);
                    }
                    cluster_coverage = score;
                }
                break;
            case 'v':
                {
                    double score = parse<double>(optarg);
                    if (score < 0) {
                        cerr << "error: [vg gaffe] Extension score threshold (" << score << ") must be positive" << endl;
                        exit(1);
                    }
                    extension_score = score;
                }
                break;
            case 'w':
                {
                    int score = parse<int>(optarg);
                    if (score < 0) {
                        cerr << "error: [vg gaffe] Extension set score threshold (" << score << ") must be positive" << endl;
                        exit(1);
                    }
                    extension_set = score;
                }
                break;
                
            case 'O':
                do_dp = false;
                break;
                
            case OPT_TRACK_PROVENANCE:
                track_provenance = true;
                break;
            
            case OPT_TRACK_CORRECTNESS:
                track_provenance = true;
                track_correctness = true;
                break;
                
            case 't':
            {
                int num_threads = parse<int>(optarg);
                if (num_threads <= 0) {
                    cerr << "error:[vg gaffe] Thread count (-t) set to " << num_threads << ", must set to a positive integer." << endl;
                    exit(1);
                }
                threads_to_run.push_back(num_threads);
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
    
    
    if (xg_name.empty() && graph_name.empty()) {
        cerr << "error:[vg gaffe] Mapping requires an XG index (-x) or a GBWTGraph (-g)" << endl;
        exit(1);
    }
    
    if (track_correctness && xg_name.empty()) {
        cerr << "error:[vg gaffe] Tracking correctness requires and XG index (-x)" << endl;
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
    
    if (distance_name.empty()) {
        cerr << "error:[vg gaffe] Mapping requires a distance index (-d)" << endl;
        exit(1);
    }
    
    // create in-memory objects
    if (progress && !xg_name.empty()) {
        cerr << "Loading XG index " << xg_name << endl;
    }
    unique_ptr<PathPositionHandleGraph> xg_index = (xg_name.empty() ? nullptr : vg::io::VPKG::load_one<PathPositionHandleGraph>(xg_name));

    if (progress) {
        cerr << "Loading GBWT index " << gbwt_name << endl;
    }
    unique_ptr<gbwt::GBWT> gbwt_index = vg::io::VPKG::load_one<gbwt::GBWT>(gbwt_name);

    if (progress) {
        cerr << "Loading minimizer index " << minimizer_name << endl;
    }
    unique_ptr<MinimizerIndex> minimizer_index = vg::io::VPKG::load_one<MinimizerIndex>(minimizer_name);

    if (progress) {
        cerr << "Loading distance index " << distance_name << endl;
    }
    MinimumDistanceIndex distance_index;
    ifstream dist_in (distance_name);
    distance_index.load(dist_in);
    //unique_ptr<MinimumDistanceIndex> distance_index = vg::io::VPKG::load_one<MinimumDistanceIndex>(distance_name);
    
    // Build or load the GBWTGraph.
    unique_ptr<GBWTGraph> gbwt_graph = nullptr;
    if (graph_name.empty()) {
        if (progress) {
            cerr << "Building GBWTGraph" << endl;
        }
        gbwt_graph.reset(new GBWTGraph(*gbwt_index, *xg_index));
    } else {
        if (progress) {
            cerr << "Loading GBWTGraph " << graph_name << endl;
        }
        gbwt_graph = vg::io::VPKG::load_one<GBWTGraph>(graph_name);
        gbwt_graph->set_gbwt(*gbwt_index);
    }

    // Set up the mapper
    if (progress) {
        cerr << "Initializing MinimizerMapper" << endl;
    }
    MinimizerMapper minimizer_mapper(*gbwt_graph, *minimizer_index, distance_index, xg_index.get());


    if (progress) {
        cerr << "--hit-cap " << hit_cap << endl;
    }
    minimizer_mapper.hit_cap = hit_cap;

    if (progress) {
        cerr << "--hard-hit-cap " << hard_hit_cap << endl;
    }
    minimizer_mapper.hard_hit_cap = hard_hit_cap;

    if (progress) {
        cerr << "--score-fraction " << minimizer_score_fraction << endl;
    }
    minimizer_mapper.minimizer_score_fraction = minimizer_score_fraction;

    if (progress) {
        cerr << "--max-extensions " << max_extensions << endl;
    }
    minimizer_mapper.max_extensions = max_extensions;

    if (progress) {
        cerr << "--max-alignments " << max_alignments << endl;
    }
    minimizer_mapper.max_alignments = max_alignments;

    if (progress) {
        cerr << "--cluster-score " << cluster_score << endl;
    }
    minimizer_mapper.cluster_score_threshold = cluster_score;

    if (progress) {
        cerr << "--cluster-coverage " << cluster_coverage << endl;
    }
    minimizer_mapper.cluster_coverage_threshold = cluster_coverage;

    if (progress) {
        cerr << "--extension-score " << extension_score << endl;
    }
    minimizer_mapper.extension_score_threshold = extension_score;

    if (progress) {
        cerr << "--extension-set " << extension_set << endl;
    }
    minimizer_mapper.extension_set_score_threshold = extension_set;

    if (progress && !do_dp) {
        cerr << "--no-dp " << endl;
    }
    minimizer_mapper.do_dp = do_dp;

    if (progress) {
        cerr << "--max-multimaps " << max_multimaps << endl;
    }
    minimizer_mapper.max_multimaps = max_multimaps;

    if (progress) {
        cerr << "--distance-limit " << distance_limit << endl;
    }
    minimizer_mapper.distance_limit = distance_limit;
    
    if (progress && track_provenance) {
        cerr << "--track-provenance " << endl;
    }
    minimizer_mapper.track_provenance = track_provenance;
    
    if (progress && track_correctness) {
        cerr << "--track-correctness " << endl;
    }
    minimizer_mapper.track_correctness = track_correctness;

    minimizer_mapper.sample_name = sample_name;
    minimizer_mapper.read_group = read_group;

    std::chrono::time_point<std::chrono::system_clock> init = std::chrono::system_clock::now();
    std::chrono::duration<double> init_seconds = init - launch;
    if (progress) {
        cerr << "Loading and initialization: " << init_seconds.count() << " seconds" << endl;
    }


    if (threads_to_run.empty()) {
        // Use 0 to represent the default.
        threads_to_run.push_back(0);
    }
    
    for (size_t thread_count_current : threads_to_run) {
    
        if (thread_count_current != 0) {
            omp_set_num_threads(thread_count_current);
        }
        
        // Work out the number of threads we will have
        size_t thread_count = omp_get_max_threads();

        // Set up counters per-thread for total reads mapped
        vector<size_t> reads_mapped_by_thread(thread_count, 0);
        
        // Have a place to log start time
        std::chrono::time_point<std::chrono::system_clock> start;
        
        {
            // Set up output to an emitter that will handle serialization
            // Discard alignments if asked to. Otherwise spit them out in GAM format.
            unique_ptr<AlignmentEmitter> alignment_emitter = discard_alignments ?
                make_unique<NullAlignmentEmitter>() :
                get_alignment_emitter("-", "GAM", {}, thread_count);

#ifdef USE_CALLGRIND
            // We want to profile the alignment, not the loading.
            CALLGRIND_START_INSTRUMENTATION;
#endif

            // Start timing overall mapping time now that indexes are loaded.
            start = std::chrono::system_clock::now();
            
            // Define how to align and output a read, in a thread.
            auto map_read = [&](Alignment& aln) {
                // Map the read with the MinimizerMapper.
                minimizer_mapper.map(aln, *alignment_emitter);
                // Record that we mapped a read.
                reads_mapped_by_thread.at(omp_get_thread_num())++;
            };
                
            for (auto& gam_name : gam_filenames) {
                // For every GAM file to remap
                get_input_file(gam_name, [&](istream& in) {
                    // Open it and map all the reads in parallel.
                    vg::io::for_each_parallel<Alignment>(in, map_read);
                });
            }
            
            for (auto& fastq_name : fastq_filenames) {
                // For every FASTQ file to map, map all its reads in parallel.
                fastq_unpaired_for_each_parallel(fastq_name, map_read);
            }
        
        } // Make sure alignment emitter is destroyed and all alignments are on disk.
        
        // Now mapping is done
        std::chrono::time_point<std::chrono::system_clock> end = std::chrono::system_clock::now();
        std::chrono::duration<double> elapsed_seconds = end-start;
        
        // How many reads did we map?
        size_t total_reads_mapped = 0;
        for (auto& reads_mapped : reads_mapped_by_thread) {
            total_reads_mapped += reads_mapped;
        }
        
        // Produce a report
        if (progress) {
            cerr << "Mapped " << total_reads_mapped << " reads across "
                << thread_count << " threads in "
                << elapsed_seconds.count() << " seconds." << endl;
            
            cerr << "Mapping speed: " << ((total_reads_mapped / elapsed_seconds.count()) / thread_count)
                << " reads per second per thread" << endl;

            cerr << "Memory footprint: " << gbwt::inGigabytes(gbwt::memoryUsage()) << " GB" << endl;
        }
        
    }
        
    return 0;
}

// Register subcommand
static Subcommand vg_gaffe("gaffe", "Graph Alignment Format Fast Emitter", DEVELOPMENT, main_gaffe);


