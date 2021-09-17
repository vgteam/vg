/** \file call_main.cpp
 *
 * Defines the "vg call" subcommand, which calls variation from an augmented graph
 */

#include <omp.h>
#include <unistd.h>
#include <getopt.h>
#include <regex>
#include <list>
#include <fstream>

#include "subcommand.hpp"
#include "../path.hpp"
#include "../graph_caller.hpp"
#include "../integrated_snarl_finder.hpp"
#include "../xg.hpp"
#include <vg/io/stream.hpp>
#include <vg/io/vpkg.hpp>
#include <bdsg/overlays/overlay_helper.hpp>

using namespace std;
using namespace vg;
using namespace vg::subcommand;

void help_call(char** argv) {
  cerr << "usage: " << argv[0] << " call [options] <graph> > output.vcf" << endl
       << "Call variants or genotype known variants" << endl
       << endl
       << "support calling options:" << endl
       << "    -k, --pack FILE          Supports created from vg pack for given input graph" << endl
       << "    -m, --min-support M,N    Minimum allele support (M) and minimum site support (N) for call [default = 2,4]" << endl
       << "    -e, --baseline-error X,Y Baseline error rates for Poisson model for small (X) and large (Y) variants [default= 0.005,0.01]" << endl
       << "    -I, --insertion-bias     Error rate multiplier for insertions for small (X) and large (Y) variants to compensate for breakpoint uncertainty [defaul=1,10]" << endl
       << "    -B, --bias-mode          Use old ratio-based genotyping algorithm as opposed to porbablistic model" << endl
       << "    -b, --het-bias M,N       Homozygous alt/ref allele must have >= M/N times more support than the next best allele [default = 6,6]" << endl
       << "GAF options:" << endl
       << "    -G, --gaf               Output GAF genotypes instead of VCF" << endl
       << "    -T, --traversals        Output all candidate traversals in GAF without doing any genotyping" << endl
       << "    -M, --trav-padding N    Extend each flank of traversals (from -T) with reference path by N bases if possible" << endl
       << "general options:" << endl
       << "    -v, --vcf FILE          VCF file to genotype (must have been used to construct input graph with -a)" << endl
       << "    -a, --genotype-snarls   Genotype every snarl, including reference calls (use to compare multiple samples)" << endl
       << "    -f, --ref-fasta FILE    Reference fasta (required if VCF contains symbolic deletions or inversions)" << endl
       << "    -i, --ins-fasta FILE    Insertions fasta (required if VCF contains symbolic insertions)" << endl
       << "    -s, --sample NAME       Sample name [default=SAMPLE]" << endl
       << "    -r, --snarls FILE       Snarls (from vg snarls) to avoid recomputing." << endl
       << "    -g, --gbwt FILE         Only call genotypes that are present in given GBWT index." << endl
       << "    -p, --ref-path NAME     Reference path to call on (multipile allowed.  defaults to all paths)" << endl
       << "    -o, --ref-offset N      Offset in reference path (multiple allowed, 1 per path)" << endl
       << "    -l, --ref-length N      Override length of reference in the contig field of output VCF" << endl
       << "    -d, --ploidy N          Ploidy of sample.  Only 1 and 2 supported. (default: 2)" << endl
       << "    -R, --ploidy-regex RULES    use the given comma-separated list of colon-delimited REGEX:PLOIDY rules to assign" << endl
       << "                                ploidies to contigs not visited by the selected samples, or to all contigs simulated" << endl
       << "                                from if no samples are used. Unmatched contigs get ploidy 2 (or that from -d)." << endl
       << "    -n, --nested            Activate nested calling mode (experimental)" << endl
       << "    -t, --threads N         number of threads to use" << endl;
}    

int main_call(int argc, char** argv) {

    string pack_filename;
    string vcf_filename;
    string sample_name = "SAMPLE";
    string snarl_filename;
    string gbwt_filename;
    string ref_fasta_filename;
    string ins_fasta_filename;
    vector<string> ref_paths;
    vector<size_t> ref_path_offsets;
    vector<size_t> ref_path_lengths;
    string min_support_string;
    string baseline_error_string;
    string insertion_bias_string;
    string bias_string;
    // require at least some support for all breakpoint edges
    // inceases sv precision, but at some recall cost.
    // think this is worth leaving on by default and not adding an option (famouse last words)
    bool expect_bp_edges = true;
    bool ratio_caller = false;
    bool legacy = false;
    int ploidy = 2;
    // copied over from vg sim
    std::vector<std::pair<std::regex, size_t>> ploidy_rules;

    bool traversals_only = false;
    bool gaf_output = false;
    size_t trav_padding = 0;
    bool genotype_snarls = false;
    bool nested = false;

    // constants
    const size_t avg_trav_threshold = 50;
    const size_t avg_node_threshold = 50;
    const size_t min_depth_bin_width = 50;
    const size_t max_depth_bin_width = 50000000;
    const double depth_scale_fac = 1.5;
    const size_t max_yens_traversals = traversals_only ? 100 : 50;
    // used to merge up snarls from chains when generating traversals
    const size_t max_chain_edges = 1000; 
    const size_t max_chain_trivial_travs = 5;
    // used to decide if site is insertion in order to apply bias
    // (is insertion if alt allele len > ref_allele len * insertion_threshold)
    const double insertion_threshold = 5;
    
    int c;
    optind = 2; // force optind past command positional argument
    while (true) {

        static const struct option long_options[] = {
            {"pack", required_argument, 0, 'k'},
            {"bias-mode", no_argument, 0, 'B'},
            {"baseline-error", required_argument, 0, 'e'},
            {"insertion-bias", required_argument, 0, 'I'},
            {"het-bias", required_argument, 0, 'b'},
            {"min-support", required_argument, 0, 'm'},
            {"vcf", required_argument, 0, 'v'},
            {"genotype-snarls", no_argument, 0, 'a'},
            {"ref-fasta", required_argument, 0, 'f'},
            {"ins-fasta", required_argument, 0, 'i'},
            {"sample", required_argument, 0, 's'},            
            {"snarls", required_argument, 0, 'r'},
            {"gbwt", required_argument, 0, 'g'},
            {"ref-path", required_argument, 0, 'p'},
            {"ref-offset", required_argument, 0, 'o'},
            {"ref-length", required_argument, 0, 'l'},
            {"ploidy", required_argument, 0, 'd'},
            {"ploidy-regex", required_argument, 0, 'R'},
            {"gaf", no_argument, 0, 'G'},
            {"traversals", no_argument, 0, 'T'},
            {"min-trav-len", required_argument, 0, 'M'},
            {"legacy", no_argument, 0, 'L'},
            {"nested", no_argument, 0, 'n'},
            {"threads", required_argument, 0, 't'},
            {"help", no_argument, 0, 'h'},
            {0, 0, 0, 0}
        };

        int option_index = 0;

        c = getopt_long (argc, argv, "k:Be:I:b:m:v:af:i:s:r:g:p:o:l:d:R:GTLM:nt:h",
                         long_options, &option_index);

        // Detect the end of the options.
        if (c == -1)
            break;

        switch (c)
        {
        case 'k':
            pack_filename = optarg;
            break;
        case 'B':
            ratio_caller = true;
            break;
        case 'b':
            bias_string = optarg;
            break;
        case 'm':
            min_support_string = optarg;
            break;
        case 'I':
            insertion_bias_string = optarg;
            break;            
        case 'e':
            baseline_error_string = optarg;
            break;            
        case 'v':
            vcf_filename = optarg;
            break;
        case 'a':
            genotype_snarls = true;
            break;
        case 'f':
            ref_fasta_filename = optarg;
            break;
        case 'i':
            ins_fasta_filename = optarg;
            break;
        case 's':
            sample_name = optarg;
            break;
        case 'r':
            snarl_filename = optarg;
            break;
        case 'g':
            gbwt_filename = optarg;
            break;
        case 'p':
            ref_paths.push_back(optarg);
            break;
        case 'o':
            ref_path_offsets.push_back(parse<int>(optarg));
            break;
        case 'l':
            ref_path_lengths.push_back(parse<int>(optarg));
            break;
        case 'd':
            ploidy = parse<int>(optarg);
            break;
        case 'R':
            for (auto& rule : split_delims(optarg, ",")) {
                // For each comma-separated rule
                auto parts = split_delims(rule, ":");
                if (parts.size() != 2) {
                    cerr << "error: ploidy rules must be REGEX:PLOIDY" << endl;
                    exit(1);
                }
                try {
                    // Parse the regex
                    std::regex match(parts[0]);
                    size_t weight = parse<size_t>(parts[1]);
                    // Save the rule
                    ploidy_rules.emplace_back(match, weight);
                } catch (const std::regex_error& e) {
                    // This is not a good regex
                    cerr << "error: unacceptable regular expression \"" << parts[0] << "\": " << e.what() << endl;
                    exit(1);
                }
            }
            break;            
        case 'G':
            gaf_output = true;
            break;
        case 'T':
            traversals_only = true;
            gaf_output = true;
            break;
        case 'M':
            trav_padding = parse<size_t>(optarg);
            break;
        case 'L':
            legacy = true;
            break;
        case 'n':
            nested =true;
            break;
        case 't':
        {
            int num_threads = parse<int>(optarg);
            if (num_threads <= 0) {
                cerr << "error:[vg call] Thread count (-t) set to " << num_threads << ", must set to a positive integer." << endl;
                exit(1);
            }
            omp_set_num_threads(num_threads);
            break;
        }
        case 'h':
        case '?':
            /* getopt_long already printed an error message. */
            help_call(argv);
            exit(1);
            break;
        default:
            abort ();
        }
    }

    if (argc <= 2) {
        help_call(argv);
        return 1;
    }

    // parse the supports (stick together to keep number of options down)
    vector<string> support_toks = split_delims(min_support_string, ",");
    double min_allele_support = -1;
    double min_site_support = -1;
    if (support_toks.size() >= 1) {
        min_allele_support = parse<double>(support_toks[0]);
        min_site_support = min_allele_support;
    }
    if (support_toks.size() == 2) {
        min_site_support = parse<double>(support_toks[1]);
    } else if (support_toks.size() > 2) {
        cerr << "error [vg call]: -m option expects at most two comma separated numbers M,N" << endl;
        return 1;
    }
    // parse the biases
    vector<string> bias_toks = split_delims(bias_string, ",");
    double het_bias = -1;
    double ref_het_bias = -1;
    if (bias_toks.size() >= 1) {
        het_bias = parse<double>(bias_toks[0]);
        ref_het_bias = het_bias;
    }
    if (bias_toks.size() == 2) {
        ref_het_bias = parse<double>(bias_toks[1]);
    } else if (bias_toks.size() > 2) {
        cerr << "error [vg call]: -b option expects at most two comma separated numbers M,N" << endl;
        return 1;
    }
    // parse the baseline errors (defaults are in snarl_caller.hpp)
    vector<string> error_toks = split_delims(baseline_error_string, ",");
    double baseline_error_large = -1;
    double baseline_error_small = -1;
    if (error_toks.size() == 2) {
        baseline_error_small = parse<double>(error_toks[0]);
        baseline_error_large = parse<double>(error_toks[1]);
        if (baseline_error_small > baseline_error_large) {
            cerr << "warning [vg call]: with baseline error -e X,Y option, small variant error (X) normally less than large (Y)" << endl;
        }
    } else if (error_toks.size() != 0) {
        cerr << "error [vg call]: -e option expects exactly two comma-separated numbers X,Y" << endl;
        return 1;
    }

    // parse the insertion bias (defaults are in snarl_caller.hpp)
    error_toks = split_delims(insertion_bias_string, ",");
    double insertion_bias_large = -1;
    double insertion_bias_small = -1;
    if (error_toks.size() == 2) {
        insertion_bias_small = parse<double>(error_toks[0]);
        insertion_bias_large = parse<double>(error_toks[1]);
    } else if (error_toks.size() != 0) {
        cerr << "error [vg call]: -I option expects exactly two comma-separated numbers X,Y" << endl;
        return 1;
    }
    
    if (trav_padding > 0 && traversals_only == false) {
        cerr << "error [vg call]: -M option can only be used in conjunction with -T" << endl;
        return 1;
    }

    if (!vcf_filename.empty() && genotype_snarls) {
        cerr << "error [vg call]: -v and -a options cannot be used together" << endl;
        return 1;
    }
    
    // Read the graph
    unique_ptr<PathHandleGraph> path_handle_graph;
    string graph_filename = get_input_file_name(optind, argc, argv);
    path_handle_graph = vg::io::VPKG::load_one<PathHandleGraph>(graph_filename);
    PathHandleGraph* graph = path_handle_graph.get();

    // Apply overlays as necessary
    bool need_path_positions = vcf_filename.empty();
    bool need_vectorizable = !pack_filename.empty();
    bdsg::PathPositionOverlayHelper pp_overlay_helper;
    bdsg::PathPositionVectorizableOverlayHelper ppv_overlay_helper;
    bdsg::PathVectorizableOverlayHelper pv_overlay_helper;
    if (need_path_positions && need_vectorizable) {
        graph = dynamic_cast<PathHandleGraph*>(ppv_overlay_helper.apply(graph));
    } else if (need_path_positions && !need_vectorizable) {
        graph = dynamic_cast<PathHandleGraph*>(pp_overlay_helper.apply(graph));
    } else if (!need_path_positions && need_vectorizable) {
        graph = dynamic_cast<PathHandleGraph*>(pv_overlay_helper.apply(graph));
    }
    
    // Check our paths
    for (const string& ref_path : ref_paths) {
        if (!graph->has_path(ref_path)) {
            cerr << "error [vg call]: Reference path \"" << ref_path << "\" not found in graph" << endl;
            return 1;
        }
    }
    // Check our offsets
    if (ref_path_offsets.size() != 0 && ref_path_offsets.size() != ref_paths.size()) {
        cerr << "error [vg call]: when using -o, the same number paths must be given with -p" << endl;
        return 1;
    }
    if (!ref_path_offsets.empty() && !vcf_filename.empty()) {
        cerr << "error [vg call]: -o cannot be used with -v" << endl;
        return 1;
    }
    // Check our ref lengths
    if (ref_path_lengths.size() != 0 && ref_path_lengths.size() != ref_paths.size()) {
        cerr << "error [vg call]: when using -l, the same number paths must be given with -p" << endl;
        return 1;
    }
    // Check bias option
    if (!bias_string.empty() && !ratio_caller) {
        cerr << "error [vg call]: -b can only be used with -B" << endl;
        return 1;
    }
    // Check ploidy option
    if (ploidy < 1 || ploidy > 2) {
        cerr << "error [vg call]: ploidy (-d) must be either 1 or 2" << endl;
        return 1;
    }
    if (ratio_caller == true && ploidy != 2) {
        cerr << "error [vg call]: ploidy (-d) must be 2 when using ratio caller (-B)" << endl;
        return 1;
    }
    if (legacy == true && ploidy != 2) {
        cerr << "error [vg call]: ploidy (-d) must be 2 when using legacy caller (-L)" << endl;
        return 1;
    }
    if (!vcf_filename.empty() && !gbwt_filename.empty()) {
        cerr << "error [vg call]: gbwt (-g) cannot be used when genotyping VCF (-v)" << endl;
        return 1;
    }
    if (legacy == true && !gbwt_filename.empty()) {
        cerr << "error [vg call]: gbwt (-g) cannot be used with legacy caller (-L)" << endl;
        return 1;
    }

    // No paths specified: use them all
    if (ref_paths.empty()) {
        graph->for_each_path_handle([&](path_handle_t path_handle) {
                const string& name = graph->get_path_name(path_handle);
                if (!Paths::is_alt(name)) {
                    ref_paths.push_back(name);
                }
            });
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
    if (!snarl_filename.empty()) {
        ifstream snarl_file(snarl_filename.c_str());
        if (!snarl_file) {
            cerr << "Error [vg call]: Unable to load snarls file: " << snarl_filename << endl;
            return 1;
        }
        snarl_manager = vg::io::VPKG::load_one<SnarlManager>(snarl_file);
    } else {
        IntegratedSnarlFinder finder(*graph);
        snarl_manager = unique_ptr<SnarlManager>(new SnarlManager(std::move(finder.find_snarls_parallel())));
    }
    
    // Make a Packed Support Caller
    unique_ptr<SnarlCaller> snarl_caller;
    algorithms::BinnedDepthIndex depth_index;

    unique_ptr<Packer> packer;
    unique_ptr<TraversalSupportFinder> support_finder;
    if (!pack_filename.empty()) {        
        // Load our packed supports (they must have come from vg pack on graph)
        packer = unique_ptr<Packer>(new Packer(graph));
        packer->load_from_file(pack_filename);
        if (nested) {
            // Make a nested packed traversal support finder (using cached veresion important for poisson caller)
            support_finder.reset(new NestedCachedPackedTraversalSupportFinder(*packer, *snarl_manager));
        } else {
            // Make a packed traversal support finder (using cached veresion important for poisson caller)
            support_finder.reset(new CachedPackedTraversalSupportFinder(*packer, *snarl_manager));
        }
                
        // need to use average support when genotyping as small differences in between sample and graph
        // will lead to spots with 0-support, espeically in and around SVs. 
        support_finder->set_support_switch_threshold(avg_trav_threshold, avg_node_threshold);

        // upweight breakpoint edges even when taking average support otherwise
        support_finder->set_min_bp_edge_override(expect_bp_edges);

        // todo: toggle between min / average (or thresholds) via command line
        
        SupportBasedSnarlCaller* packed_caller = nullptr;

        if (ratio_caller == false) {
            // Make a depth index
            depth_index = algorithms::binned_packed_depth_index(*packer, ref_paths, min_depth_bin_width, max_depth_bin_width,
                                                                depth_scale_fac, 0, true, true);
            // Make a new-stype probablistic caller
            auto poisson_caller = new PoissonSupportSnarlCaller(*graph, *snarl_manager, *support_finder, depth_index,
                                                                //todo: qualities need to be used better in conjunction with
                                                                //expected depth.
                                                                //packer->has_qualities());
                                                                false);

            // Pass the errors through
            poisson_caller->set_baseline_error(baseline_error_small, baseline_error_large);
            poisson_caller->set_insertion_bias(insertion_threshold, insertion_bias_small, insertion_bias_large);
                
            packed_caller = poisson_caller;
        } else {
            // Make an old-style ratio support caller
            auto ratio_caller = new RatioSupportSnarlCaller(*graph, *snarl_manager, *support_finder);
            if (het_bias >= 0) {
                ratio_caller->set_het_bias(het_bias, ref_het_bias);
            }
            packed_caller = ratio_caller;
        }
        if (min_allele_support >= 0) {
            packed_caller->set_min_supports(min_allele_support, min_allele_support, min_site_support);
        }
        
        snarl_caller = unique_ptr<SnarlCaller>(packed_caller);
    }

    if (!snarl_caller) {
        cerr << "error [vg call]: pack file (-k) is required" << endl;
        return 1;
    }

    unique_ptr<AlignmentEmitter> alignment_emitter;
    if (gaf_output) {
      alignment_emitter = vg::io::get_non_hts_alignment_emitter("-", "GAF", {}, get_thread_count(), graph);
    }

    unique_ptr<GraphCaller> graph_caller;
    unique_ptr<TraversalFinder> traversal_finder;
    unique_ptr<gbwt::GBWT> gbwt_index;

    vcflib::VariantCallFile variant_file;
    unique_ptr<FastaReference> ref_fasta;
    unique_ptr<FastaReference> ins_fasta;
    if (!vcf_filename.empty()) {
        // Genotype the VCF
        variant_file.parseSamples = false;
        variant_file.open(vcf_filename);
        if (!variant_file.is_open()) {
            cerr << "error: [vg call] could not open " << vcf_filename << endl;
            return 1;
        }

        // load up the fasta
        if (!ref_fasta_filename.empty()) {
            ref_fasta = unique_ptr<FastaReference>(new FastaReference);
            ref_fasta->open(ref_fasta_filename);
        }
        if (!ins_fasta_filename.empty()) {
            ins_fasta = unique_ptr<FastaReference>(new FastaReference);
            ins_fasta->open(ins_fasta_filename);
        }
        
        VCFGenotyper* vcf_genotyper = new VCFGenotyper(*graph, *snarl_caller,
                                                       *snarl_manager, variant_file,
                                                       sample_name, ref_paths, ref_path_ploidies,
                                                       ref_fasta.get(),
                                                       ins_fasta.get(),
                                                       alignment_emitter.get(),
                                                       traversals_only,
                                                       gaf_output,
                                                       trav_padding);
        graph_caller = unique_ptr<GraphCaller>(vcf_genotyper);
    } else if (legacy) {
        // de-novo caller (port of the old vg call code, which requires a support based caller)
        LegacyCaller* legacy_caller = new LegacyCaller(*dynamic_cast<PathPositionHandleGraph*>(graph),
                                                       *dynamic_cast<SupportBasedSnarlCaller*>(snarl_caller.get()),
                                                       *snarl_manager,
                                                       sample_name, ref_paths, ref_path_offsets, ref_path_ploidies);
        graph_caller = unique_ptr<GraphCaller>(legacy_caller);
    } else {
        // flow caller can take any kind of traversal finder.  two are supported for now:
        
        if (!gbwt_filename.empty()) {
            // GBWT traversals
            gbwt_index = vg::io::VPKG::load_one<gbwt::GBWT>(gbwt_filename);
            if (gbwt_index.get() == nullptr) {
                cerr << "error:[vg call] unable to load gbwt index file: " << gbwt_filename << endl;
                return 1;
            }
            GBWTTraversalFinder* gbwt_traversal_finder = new GBWTTraversalFinder(*graph, *gbwt_index.get());
            traversal_finder = unique_ptr<TraversalFinder>(gbwt_traversal_finder);
        } else {
            // Flow traversals (Yen's algorithm)
            
            // todo: do we ever want to toggle in min-support?
            function<double(handle_t)> node_support = [&] (handle_t h) {
                return support_finder->support_val(support_finder->get_avg_node_support(graph->get_id(h)));
            };
            
            function<double(edge_t)> edge_support = [&] (edge_t e) {
                return support_finder->support_val(support_finder->get_edge_support(e));
            };

            // create the flow traversal finder
            FlowTraversalFinder* flow_traversal_finder = new FlowTraversalFinder(*graph, *snarl_manager, max_yens_traversals,
                                                                                 node_support, edge_support);
            traversal_finder = unique_ptr<TraversalFinder>(flow_traversal_finder);
        }

        if (nested) {
            graph_caller.reset(new NestedFlowCaller(*dynamic_cast<PathPositionHandleGraph*>(graph),
                                                    *dynamic_cast<SupportBasedSnarlCaller*>(snarl_caller.get()),
                                                    *snarl_manager,
                                                    sample_name, *traversal_finder, ref_paths, ref_path_offsets,
                                                    ref_path_ploidies,
                                                    alignment_emitter.get(),
                                                    traversals_only,
                                                    gaf_output,
                                                    trav_padding,
                                                    genotype_snarls));
        } else {
            graph_caller.reset(new FlowCaller(*dynamic_cast<PathPositionHandleGraph*>(graph),
                                              *dynamic_cast<SupportBasedSnarlCaller*>(snarl_caller.get()),
                                              *snarl_manager,
                                              sample_name, *traversal_finder, ref_paths, ref_path_offsets,
                                              ref_path_ploidies,
                                              alignment_emitter.get(),
                                              traversals_only,
                                              gaf_output,
                                              trav_padding,
                                              genotype_snarls));            
        }
    }

    string header;
    if (!gaf_output) {
        // Init The VCF       
        VCFOutputCaller* vcf_caller = dynamic_cast<VCFOutputCaller*>(graph_caller.get());
        assert(vcf_caller != nullptr);
        header = vcf_caller->vcf_header(*graph, ref_paths, ref_path_lengths);
    }

    // Call the graph
    if (!traversals_only) {

        // Call each snarl
        // (todo: try chains in normal mode)
        graph_caller->call_top_level_snarls(*graph, ploidy);
    } else {
        // Attempt to call chains instead of snarls so that the output traversals are longer
        // Todo: this could probably help in some cases when making VCFs too
        graph_caller->call_top_level_chains(*graph, ploidy, max_chain_edges,  max_chain_trivial_travs);
    }

    if (!gaf_output) {
        // Output VCF
        VCFOutputCaller* vcf_caller = dynamic_cast<VCFOutputCaller*>(graph_caller.get());
        assert(vcf_caller != nullptr);
        cout << header << flush;
        vcf_caller->write_variants(cout);
    }
    
    return 0;
}

// Register subcommand
static Subcommand vg_call("call", "call or genotype VCF variants", PIPELINE, 10, main_call);

