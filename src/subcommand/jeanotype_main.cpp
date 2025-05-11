/** \file jeanotype_main.cpp
 *
 * Defines the "vg jeanotype" subcommand, which genotypes variation from a pangenome
 */

#include <getopt.h>
#include "subcommand.hpp"
#include "../utility.hpp"
#include "../graph_caller.hpp"
#include "../gbzgraph.hpp"
#include <vg/io/vpkg.hpp>
#include "../integrated_snarl_finder.hpp"
#include <bdsg/overlays/overlay_helper.hpp>

using namespace std;
using namespace vg;
using namespace vg::subcommand;

/// Helper to ensure that a PathHandleGraph has the VectorizableHandleGraph and PathPositionHandleGraph interfaces.
typedef bdsg::PairOverlayHelper<PathPositionHandleGraph, bdsg::PackedReferencePathOverlay, PathHandleGraph,
                                VectorizableHandleGraph, bdsg::PathPositionVectorizableOverlay, PathPositionHandleGraph> ReferencePathVectorizableOverlayHelper;

void help_jeanotype(char** argv) {
  cerr << "usage: " << argv[0] << " jeanotype [options] <graph> > output.vcf" << endl
       << "Genotype variants" << endl
       << endl
       << "    -g, --gaf FILE          read alignments from this indexed bgzipped GAF file" << endl
       << "    -b, --bed-name FILE     BED defining regions to call together (short tandem repeats)" << endl
       << "    -r, --snarls FILE       Snarls (from vg snarls) to avoid recomputing." << endl
       << "    -d, --dist-name FILE    cluster using this distance index" << endl
       << "    -l, --block-size N      call blocks of snarls covering ~N bp. Default is 1000000" << endl
       << "    -a, --genotype-snarls   genotype every snarl, including reference calls (use to compare multiple samples)" << endl
       << "    -c, --min-length N      genotype only snarls with at least one traversal of length >= N" << endl
       << "    -C, --max-length N      genotype only snarls where all traversals have length <= N" << endl
       << "    -s, --sample NAME       Sample name [default=SAMPLE]" << endl
       << "    -p, --ref-path NAME     Reference path to call on (multipile allowed. defaults to all paths)" << endl
       << "    -S, --ref-sample NAME   Call on all paths with given sample name (cannot be used with -p)" << endl
       << "        --progress          Show progress" << endl
       << "    -t, --threads N         number of threads to use" << endl;
}    


int main_jeanotype(int argc, char** argv) {

    string sample_name = "SAMPLE";
    string gaf_filename;
    string dist_filename;
    string bed_filename;
    string snarl_filename;
    vector<string> ref_paths;
    string ref_sample;
    vector<size_t> ref_path_offsets;
    vector<size_t> ref_path_lengths;
    bool genotype_snarls = false;
    size_t min_allele_len = 0;
    size_t max_allele_len = numeric_limits<size_t>::max();
    int ploidy = 2;
    int block_size = 1000000;
    std::vector<std::pair<std::regex, size_t>> ploidy_rules;
    int num_threads = 1;
    bool show_progress = false;
    const int OPT_PROGRESS = 1000;

    int c;
    optind = 2; // force optind past command positional argument
    while (true) {

        static const struct option long_options[] = {
            {"gaf", required_argument, 0, 'g'},
            {"genotype-snarls", no_argument, 0, 'a'},
            {"dist-name", required_argument, 0, 'd'},
            {"min-length", required_argument, 0, 'c'},
            {"max-length", required_argument, 0, 'C'},
            {"bed-name", required_argument, 0, 'b'},
            {"sample", required_argument, 0, 's'},            
            {"snarls", required_argument, 0, 'r'},
            {"block-size", required_argument, 0, 'l'},
            {"ref-path", required_argument, 0, 'p'},
            {"ref-sample", required_argument, 0, 'S'},            
            {"threads", required_argument, 0, 't'},
            {"progress", no_argument, 0, OPT_PROGRESS },
            {"help", no_argument, 0, 'h'},
            {0, 0, 0, 0}
        };

        int option_index = 0;

        c = getopt_long (argc, argv, "ag:d:b:s:r:c:l:p:C:S:t:h",
                         long_options, &option_index);

        // Detect the end of the options.
        if (c == -1)
            break;

        switch (c)
        {
        case 'g':
            if (!optarg || !*optarg) {
                cerr << "error: Must provide input GAF file with -a." << endl;
                exit(1);
            }
            gaf_filename = optarg;
            break;
        case 'd':
            if (!optarg || !*optarg) {
                cerr << "error:Must provide distance index file with -d." << endl;
                exit(1);
            }
            dist_filename = optarg;
            break;
        case 'a':
            genotype_snarls = true;
            break;
        case 'b':
            bed_filename = optarg;
            break;
        case 'l':
            block_size = parse<int>(optarg);;
            break;
        case 'c':
            min_allele_len = parse<size_t>(optarg);
            break;
        case 'C':
            max_allele_len = parse<size_t>(optarg);
            break;
        case 's':
            sample_name = optarg;
            break;
        case 'r':
            snarl_filename = optarg;
            break;
        case 'p':
            ref_paths.push_back(optarg);
            break;
        case 'S':
            ref_sample = optarg;
            break;            
        case 't':
        {
            num_threads = parse<int>(optarg);
            if (num_threads <= 0) {
                cerr << "error:[vg jeanotype] Thread count (-t) set to " << num_threads << ", must set to a positive integer." << endl;
                exit(1);
            }
            break;
        }
        case OPT_PROGRESS:
            show_progress = true;
            break;
        case 'h':
        case '?':
            /* getopt_long already printed an error message. */
            help_jeanotype(argv);
            exit(1);
            break;
        default:
            abort ();
        }
    }

    if (argc <= 2) {
        help_jeanotype(argv);
        return 1;
    }

    if (!ref_paths.empty() && !ref_sample.empty()) {
        cerr << "error [vg jeanotype]: -S cannot be used with -p" << endl;
        return 1;
    }

    
    if (show_progress) cerr << "Using " << num_threads << " threads." << endl;
    omp_set_num_threads(num_threads);
    
    string graph_filename = get_input_file_name(optind, argc, argv);
    
    const size_t max_chain_edges = 1000; 
    const size_t max_chain_trivial_travs = 5;

    if (show_progress) cerr << "Graph: " << graph_filename << endl;
    if (show_progress) cerr << "GAF: " << gaf_filename << endl;
    
    // load graph
    unique_ptr<GBZGraph> gbz_graph;
    gbwt::GBWT* gbwt_index = nullptr;
    PathHandleGraph* graph = nullptr;
    if (show_progress) cerr << "Loading GBZ" << endl;
    auto input = vg::io::VPKG::try_load_first<GBZGraph, PathHandleGraph>(graph_filename);
    if (get<0>(input)) {        
        gbz_graph = std::move(get<0>(input));
        graph = gbz_graph.get();
        gbwt_index = &gbz_graph->gbz.index;
        if (show_progress) cerr << "Loaded GBZ" << endl;
    } else {
        cerr << "Error: only GBZ graphs accepted." << endl;
        return 1;            
    }

    ReferencePathVectorizableOverlayHelper ppv_overlay_helper;
    graph = dynamic_cast<PathHandleGraph*>(ppv_overlay_helper.apply(graph));

    // in order to add subpath support, we let all ref_paths be subpaths and then convert coordinates
    // on vcf export.  the exception is writing the header where we need base paths. we keep
    // track of them the best we can here (just for writing the ##contigs)
    unordered_map<string, size_t> basepath_length_map;

    // call doesn't always require path positions .. .don't change that now
    function<size_t(path_handle_t)> compute_path_length = [&] (path_handle_t path_handle) {
        PathPositionHandleGraph* pp_graph = dynamic_cast<PathPositionHandleGraph*>(graph);
        if (pp_graph) {
            return pp_graph->get_path_length(path_handle);
        } else {
            size_t len = 0;
            graph->for_each_step_in_path(path_handle, [&] (step_handle_t step) {
                    len += graph->get_length(graph->get_handle_of_step(step));
                });
            return len;
        }
    };

    if (ref_paths.empty()) {
        set<string> ref_sample_names;
        graph->for_each_path_handle([&](path_handle_t path_handle) {
            const string& name = graph->get_path_name(path_handle);
            PathSense path_sense = PathMetadata::parse_sense(name);
            if (!Paths::is_alt(name) && path_sense != PathSense::HAPLOTYPE) {
                string sample_name = PathMetadata::parse_sample_name(name);
                if (ref_sample.empty() || sample_name == ref_sample) {                        
                    ref_paths.push_back(name);
                    // keep track of length best we can using maximum coordinate in event of subpaths
                    subrange_t subrange;
                    string base_name = Paths::strip_subrange(name, &subrange);
                    size_t offset = subrange == PathMetadata::NO_SUBRANGE ? 0 : subrange.first;
                    size_t& cur_len = basepath_length_map[base_name];
                    cur_len = max(cur_len, compute_path_length(path_handle) + offset);
                    if (sample_name != PathMetadata::NO_SAMPLE_NAME) {
                        ref_sample_names.insert(sample_name);
                    }
                }
            }
        });
        if (ref_sample_names.size() > 1 && ref_sample.empty()) {
            cerr << "error [vg jeanotype]: Multiple reference samples detected: [";
            size_t count = 0;
            for (const string& n : ref_sample_names) {                
                cerr << n;
                if (++count >= std::min(ref_sample_names.size(), (size_t)5)) {
                    if (ref_sample_names.size() > 5) {
                        cerr << ", ...";
                    }
                    break;
                } else {
                    cerr << ", ";
                }
            }
            cerr << "]. Please use -S to specify a single reference sample or use -p to specify reference paths";
            return 1;
        }                
    } else {
        // if paths are given, we convert them to subpaths so that ref paths list corresponds
        // to path names in graph.  subpath handling will only be handled when writing the vcf
        // (this way, we add subpath support without changing anything in between)
        vector<string> ref_subpaths;
        unordered_map<string, bool> ref_path_set;
        for (const string& ref_path : ref_paths) {
            ref_path_set[ref_path] = false;
        }
        graph->for_each_path_handle([&](path_handle_t path_handle) {
                const string& name = graph->get_path_name(path_handle);
                subrange_t subrange;
                string base_name = Paths::strip_subrange(name, &subrange);
                size_t offset = subrange == PathMetadata::NO_SUBRANGE ? 0 : subrange.first;
                if (ref_path_set.count(base_name)) {
                    ref_subpaths.push_back(name);
                    // keep track of length best we can
                    if (ref_path_lengths.empty()) {
                        size_t& cur_len = basepath_length_map[base_name];
                        cur_len = max(cur_len, compute_path_length(path_handle) + offset);
                    }
                    ref_path_set[base_name] = true;
                }
            });

        // if we have reference lengths, great!  this will be the only way to get a correct header in the presence of supbpaths
        if (!ref_path_lengths.empty()) {
            assert(ref_path_lengths.size() == ref_paths.size());
            for (size_t i = 0; i < ref_paths.size(); ++i) {
                basepath_length_map[ref_paths[i]] = ref_path_lengths[i];
            }
        }

        // Check our paths
        for (const auto& ref_path_used : ref_path_set) {
            if (!ref_path_used.second) {
                cerr << "error [vg jeanotype]: Reference path \"" << ref_path_used.first << "\" not found in graph" << endl;
                return 1;
            }
        }
        
        swap(ref_paths, ref_subpaths);
    }
        
    // make sure we have some ref paths
    if (ref_paths.empty()) {
        if (!ref_sample.empty()) {
            cerr << "error [vg jeanotype]: No paths with selected reference sample \"" << ref_sample << "\" found. "
                 << "Try using vg paths -M to see which samples are in your graph" << endl;
            return 1;
        }
        cerr << "error [vg jeanotype]: No reference paths found" << endl;
        return 1;
    }

    // build table of ploidys
    vector<int> ref_path_ploidies;
    for (const string& ref_path : ref_paths) {
        int path_ploidy = ploidy;
        for (auto& rule : ploidy_rules) {
            if (std::regex_match(ref_path, rule.first)) {
                path_ploidy = rule.second;
                break;
            }
        }
        ref_path_ploidies.push_back(path_ploidy);
    }

    // Load or compute the snarls
    unique_ptr<SnarlManager> snarl_manager;    
    // if (!snarl_filename.empty()) {
    //     ifstream snarl_file(snarl_filename.c_str());
    //     if (!snarl_file) {
    //         cerr << "Error [vg jeanotype]: Unable to load snarls file: " << snarl_filename << endl;
    //         return 1;
    //     }
    //     if (show_progress) cerr << "Loading snarls from " << snarl_filename << endl;
    //     snarl_manager = vg::io::VPKG::load_one<SnarlManager>(snarl_file);
    //     if (show_progress) cerr << "Loaded snarls" << endl;
    // } else {
    //     if (show_progress) cerr << "Computing snarls" << endl;
    //     IntegratedSnarlFinder finder(*graph);
    //     if (show_progress) cerr << "Computed snarls" << endl;
    //     snarl_manager = unique_ptr<SnarlManager>(new SnarlManager(std::move(finder.find_snarls_parallel())));
    // }
    
    unique_ptr<TraversalFinder> traversal_finder;
    GBWTTraversalFinder* gbwt_traversal_finder = new GBWTTraversalFinder(*graph, *gbwt_index);
    traversal_finder = unique_ptr<TraversalFinder>(gbwt_traversal_finder);

    // load read index
    if (show_progress) cerr << "Loading read index" << endl;
    GAFindex* reads = new GAFindex();
    reads->load(gaf_filename);   
    if (show_progress) cerr << "Loaded read index" << endl;
    
    auto read_caller = new ReadBasedSnarlCaller(*graph);
    unique_ptr<SnarlCaller> snarl_caller = unique_ptr<SnarlCaller>(read_caller);

    SnarlDistanceIndex* distance_index = new SnarlDistanceIndex();
    distance_index->deserialize(dist_filename);
    
    HapCaller* hap_caller = new HapCaller(*dynamic_cast<PathPositionHandleGraph*>(graph),
                                          *dynamic_cast<ReadBasedSnarlCaller*>(snarl_caller.get()),
                                          *snarl_manager, *distance_index, *gbwt_traversal_finder,
                                          *reads,
                                          sample_name,
                                          ref_paths, ref_path_offsets,
                                          ref_path_ploidies,
                                          genotype_snarls,
                                          make_pair(min_allele_len, max_allele_len));

    string header;
    // Init The VCF       
    VCFOutputCaller* vcf_caller = dynamic_cast<VCFOutputCaller*>(hap_caller);
    assert(vcf_caller != nullptr);
    // Make sure the basepath information we inferred above goes directy to the VCF header
    // (and that it does *not* try to read it from the graph paths)
    vector<string> header_ref_paths;
    vector<size_t> header_ref_lengths;
    bool need_overrides = dynamic_cast<VCFGenotyper*>(hap_caller) == nullptr;
    for (const auto& path_len : basepath_length_map) {
        header_ref_paths.push_back(path_len.first);
        if (need_overrides) {
            header_ref_lengths.push_back(path_len.second);
        }
    }
    header = vcf_caller->vcf_header(*graph, header_ref_paths, header_ref_lengths);
    
    hap_caller->set_show_progress(show_progress);

    // add regions if provided
    vector<Region> regions;
    if (!bed_filename.empty()) {
        parse_bed_regions(bed_filename, regions);
        hap_caller->update_regions(regions);
    }
    
    if (show_progress) cerr << "Calling snarls" << endl;
    // hap_caller->call_top_level_snarls(*graph, GraphCaller::RecurseAlways);
    // hap_caller->call_top_level_chains(*graph, 200, 200, GraphCaller::RecurseAlways);
    hap_caller->call_top_level_snarl_block(block_size);
    if (show_progress) cerr << "Called snarls" << endl;

    // use from rpvg in the snarl_caller later
    // AlignmentPathFinder<vg::Alignment> align_path_finder(paths_index, library_type, score_not_qual, use_allelic_mapq, pre_frag_length_dist.maxLength(), max_partial_offset, est_missing_noise_prob, max_score_diff, min_best_score_filter);
    // unaligned_read_count = findAlignmentPaths<vg::Alignment>(alignments_istream, align_paths_buffer_queue, align_path_finder, num_threads);


    // PathEstimator * path_estimator;
    // path_estimator = new PathGroupPosteriorEstimator(ploidy, use_hap_gibbs, prob_precision);
    // path_estimator->estimate(&(path_cluster_estimates->back().second), read_path_cluster_probs, &mt_rng);

    // mt19937 mt_rng = mt19937(rng_seed + i);
    // path_estimator->estimate(&(path_cluster_estimates->back().second), read_path_cluster_probs, &mt_rng);

    // Output VCF
    cout << header << flush;
    if (show_progress) cerr << "Writing VCF Variants" << endl;
    vcf_caller->write_variants(cout, snarl_manager.get());
    if (show_progress) cerr << "VCF complete" << endl;        
    
    return 0;
}

// Register subcommand
static Subcommand vg_jeanotype("jeanotype", "genotype variants", PIPELINE, 10, main_jeanotype);


