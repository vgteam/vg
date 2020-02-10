/** \file call_main.cpp
 *
 * Defines the "vg call" subcommand, which calls variation from an augmented graph
 */

#include <omp.h>
#include <unistd.h>
#include <getopt.h>

#include <list>
#include <fstream>

#include "subcommand.hpp"
#include "../path.hpp"
#include "../graph_caller.hpp"
#include "../xg.hpp"
#include <vg/io/stream.hpp>
#include <vg/io/vpkg.hpp>
#include <bdsg/overlay_helper.hpp>

using namespace std;
using namespace vg;
using namespace vg::subcommand;

void help_call(char** argv) {
  cerr << "usage: " << argv[0] << " call [options] <graph> > output.vcf" << endl
       << "Call variants or genotype known variants" << endl
       << endl
       << "support calling options:" << endl
       << "    -k, --pack FILE         Supports created from vg pack for given input graph" << endl
       << "    -m, --min-support M,N   Minimum allele support (M) and minimum site support (N) for call [default = 1,4]" << endl
       << "    -B, --bias-mode         Use old ratio-based genotyping algorithm as opposed to porbablistic model" << endl
       << "    -b, --het-bias M,N      Homozygous alt/ref allele must have >= M/N times more support than the next best allele [default = 6,6]" << endl
       << "general options:" << endl
       << "    -v, --vcf FILE          VCF file to genotype (must have been used to construct input graph with -a)" << endl
       << "    -f, --ref-fasta FILE    Reference fasta (required if VCF contains symbolic deletions or inversions)" << endl
       << "    -i, --ins-fasta FILE    Insertions fasta (required if VCF contains symbolic insertions)" << endl
       << "    -s, --sample NAME       Sample name [default=SAMPLE]" << endl
       << "    -r, --snarls FILE       Snarls (from vg snarls) to avoid recomputing." << endl
       << "    -p, --ref-path NAME     Reference path to call on (multipile allowed.  defaults to all paths)" << endl
       << "    -o, --ref-offset N      Offset in reference path (multiple allowed, 1 per path)" << endl
       << "    -l, --ref-length N      Override length of reference in the contig field of output VCF" << endl
       << "    -t, --threads N         number of threads to use" << endl;
}    

int main_call(int argc, char** argv) {

    string pack_filename;
    string vcf_filename;
    string sample_name = "SAMPLE";
    string snarl_filename;
    string ref_fasta_filename;
    string ins_fasta_filename;
    vector<string> ref_paths;
    vector<size_t> ref_path_offsets;
    vector<size_t> ref_path_lengths;
    string min_support_string;
    string bias_string;
    bool ratio_caller = false;
    bool legacy = false;
    
    int c;
    optind = 2; // force optind past command positional argument
    while (true) {

        static const struct option long_options[] = {
            {"pack", required_argument, 0, 'k'},
            {"bias-mode", no_argument, 0, 'B'},
            {"het-bias", required_argument, 0, 'b'},
            {"min-support", required_argument, 0, 'm'},
            {"vcf", required_argument, 0, 'v'},
            {"ref-fasta", required_argument, 0, 'f'},
            {"ins-fasta", required_argument, 0, 'i'},
            {"sample", required_argument, 0, 's'},            
            {"snarls", required_argument, 0, 'r'},
            {"ref-path", required_argument, 0, 'p'},
            {"ref-offset", required_argument, 0, 'o'},
            {"ref-length", required_argument, 0, 'l'},
            {"legacy", no_argument, 0, 'L'},
            {"threads", required_argument, 0, 't'},
            {"help", no_argument, 0, 'h'},
            {0, 0, 0, 0}
        };

        int option_index = 0;

        c = getopt_long (argc, argv, "k:Bb:m:v:f:i:s:r:p:o:l:Lt:h",
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
        case 'v':
            vcf_filename = optarg;
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
        case 'p':
            ref_paths.push_back(optarg);
            break;
        case 'o':
            ref_path_offsets.push_back(parse<int>(optarg));
            break;
        case 'l':
            ref_path_lengths.push_back(parse<int>(optarg));
            break;
        case 'L':
            legacy = true;
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
    
    // Read the graph
    unique_ptr<PathHandleGraph> path_handle_graph;
    get_input_file(optind, argc, argv, [&](istream& in) {
            path_handle_graph = vg::io::VPKG::load_one<PathHandleGraph>(in);
        });
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

    // No paths specified: use them all
    if (ref_paths.empty()) {
        graph->for_each_path_handle([&](path_handle_t path_handle) {
                const string& name = graph->get_path_name(path_handle);
                if (!Paths::is_alt(name)) {
                    ref_paths.push_back(name);
                }
            });
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
        CactusSnarlFinder finder(*graph);
        snarl_manager = unique_ptr<SnarlManager>(new SnarlManager(std::move(finder.find_snarls_parallel())));
    }
    
    unique_ptr<GraphCaller> graph_caller;
    unique_ptr<SnarlCaller> snarl_caller;
    algorithms::BinnedDepthIndex depth_index;

    // Make a Packed Support Caller
    unique_ptr<Packer> packer;
    unique_ptr<TraversalSupportFinder> support_finder;
    if (!pack_filename.empty()) {        
        // Load our packed supports (they must have come from vg pack on graph)
        packer = unique_ptr<Packer>(new Packer(graph));
        packer->load_from_file(pack_filename);
        // Make a packed traversal support finder (using cached veresion important for poisson caller)
        PackedTraversalSupportFinder* packed_support_finder = new CachedPackedTraversalSupportFinder(*packer, *snarl_manager);
        support_finder = unique_ptr<TraversalSupportFinder>(packed_support_finder);
        
        SupportBasedSnarlCaller* packed_caller = nullptr;

        if (ratio_caller == false) {
            // Make a depth index
            depth_index = algorithms::binned_packed_depth_index(*packer, ref_paths, 50, 0, true, true);
            // Make a new-stype probablistic caller
            auto poisson_caller = new PoissonSupportSnarlCaller(*graph, *snarl_manager, *packed_support_finder, depth_index);
            packed_caller = poisson_caller;
        } else {
            // Make an old-style ratio support caller
            auto ratio_caller = new RatioSupportSnarlCaller(*graph, *snarl_manager, *packed_support_finder);
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
                                                       sample_name, ref_paths,
                                                       ref_fasta.get(),
                                                       ins_fasta.get());
        graph_caller = unique_ptr<GraphCaller>(vcf_genotyper);
    } else if (legacy) {
        // de-novo caller (port of the old vg call code, which requires a support based caller)
        LegacyCaller* legacy_caller = new LegacyCaller(*dynamic_cast<PathPositionHandleGraph*>(graph),
                                                       *dynamic_cast<SupportBasedSnarlCaller*>(snarl_caller.get()),
                                                       *snarl_manager,
                                                       sample_name, ref_paths, ref_path_offsets);
        graph_caller = unique_ptr<GraphCaller>(legacy_caller);
    } else {
        FlowCaller* flow_caller = new FlowCaller(*dynamic_cast<PathPositionHandleGraph*>(graph),
                                                 *dynamic_cast<SupportBasedSnarlCaller*>(snarl_caller.get()),
                                                 *snarl_manager,
                                                 sample_name, 100, ref_paths, ref_path_offsets);
        graph_caller = unique_ptr<GraphCaller>(flow_caller);
    }

    // Call the graph
    graph_caller->call_top_level_snarls();

    // VCF output is our only supported output
    VCFOutputCaller* vcf_caller = dynamic_cast<VCFOutputCaller*>(graph_caller.get());
    assert(vcf_caller != nullptr);
    cout << vcf_caller->vcf_header(*graph, ref_paths, ref_path_lengths) << flush;
    vcf_caller->write_variants(cout);
        
    return 0;
}

// Register subcommand
static Subcommand vg_call("call", "call or genotype VCF variants", PIPELINE, 7, main_call);

