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

using namespace std;
using namespace vg;
using namespace vg::subcommand;

void help_call2(char** argv) {
  cerr << "usage: " << argv[0] << " call [options] <graph> > output.vcf" << endl
       << "Call variants or genotype known variants" << endl
       << endl
       << "support calling options:" << endl
       << "    -k, --pack FILE         Supports created from vg pack (input graph must be xg)" << endl
       << "general options:" << endl
       << "    -v, --vcf FILE          VCF file to genotype (must have been used to construct input graph with -a)" << endl
       << "    -f, --ref-fasta FILE    Reference fasta (required if VCF contains symbolic deletions or inversions)" << endl
       << "    -i, --ins-fasta FILE    Insertions fasta (required if VCF contains symbolic insertions)" << endl
       << "    -s, --sample NAME       Sample name [default=SAMPLE]" << endl
       << "    -r, --snarls FILE       Snarls (from vg snarls) to avoid recomputing." << endl
       << "    -p, --ref-path NAME     Reference path to call on (multipile allowed.  defaults to all paths)" << endl
       << "    -t, --threads N         number of threads to use" << endl;
}    

int main_call2(int argc, char** argv) {

    string pack_filename;
    string vcf_filename;
    string sample_name = "SAMPLE";
    string snarl_filename;
    string ref_fasta_filename;
    string ins_fasta_filename;
    vector<string> ref_paths;

    int c;
    optind = 2; // force optind past command positional argument
    while (true) {

        static const struct option long_options[] = {
            {"pack", required_argument, 0, 'k'},
            {"vcf", required_argument, 0, 'v'},
            {"ref-fasta", required_argument, 0, 'f'},
            {"ins-fasta", required_argument, 0, 'i'},
            {"sample", required_argument, 0, 's'},            
            {"snarls", required_argument, 0, 'r'},
            {"ref-path", required_argument, 0, 'p'},
            {"threads", required_argument, 0, 't'},
            {"help", no_argument, 0, 'h'},
            {0, 0, 0, 0}
        };

        int option_index = 0;

        c = getopt_long (argc, argv, "k:v:f:i:s:r:p:t:h",
                         long_options, &option_index);

        // Detect the end of the options.
        if (c == -1)
            break;

        switch (c)
        {
        case 'k':
            pack_filename = optarg;
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
            help_call2(argv);
            exit(1);
            break;
        default:
            abort ();
        }
    }

    if (argc <= 2) {
        help_call2(argv);
        return 1;
    }

    // Read the graph
    unique_ptr<PathHandleGraph> graph;
    get_input_file(optind, argc, argv, [&](istream& in) {
            graph = vg::io::VPKG::load_one<PathHandleGraph>(in);
        });

    // Check our paths
    for (const string& ref_path : ref_paths) {
        if (!graph->has_path(ref_path)) {
            cerr << "error [vg call]: Reference path \"" << ref_path << "\" not found in graph" << endl;
            return 1;
        }
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
        snarl_manager = unique_ptr<SnarlManager>(new SnarlManager(std::move(finder.find_snarls())));
    }
    
    unique_ptr<GraphCaller> graph_caller;
    unique_ptr<SnarlCaller> snarl_caller;

    // Make a Packed Support Caller
    unique_ptr<Packer> packer;
    if (!pack_filename.empty()) {
        if (dynamic_cast<xg::XG*>(graph.get()) == nullptr) {
            cerr << "error [vg call]: Packed support option requires input graph in XG format" << endl;
            return 1;
        }
        // Load our packed supports (they must have come from vg pack on graph)
        packer = unique_ptr<Packer>(new Packer(dynamic_cast<xg::XG*>(graph.get())));
        packer->load_from_file(pack_filename);
        PackedSupportSnarlCaller* packed_caller = new PackedSupportSnarlCaller(*packer, *snarl_manager);
        snarl_caller = unique_ptr<SnarlCaller>(packed_caller);
    }

    if (!snarl_caller) {
        cerr << "error [vg call]: pack file (-p) is required" << endl;
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
    } else {
        // de-novo caller (port of the old vg call code, which requires a support based caller)
        LegacyCaller* legacy_caller = new LegacyCaller(*dynamic_cast<PathPositionHandleGraph*>(graph.get()),
                                                       *dynamic_cast<SupportBasedSnarlCaller*>(snarl_caller.get()),
                                                        *snarl_manager,
                                                        sample_name, ref_paths);
        graph_caller = unique_ptr<GraphCaller>(legacy_caller);
    }

    // Call the graph
    cout << graph_caller->vcf_header(*graph, ref_paths) << flush;
    graph_caller->call_top_level_snarls();
        
    return 0;
}

// Register subcommand
static Subcommand vg_call2("call2", "call2 variants", PIPELINE, 5, main_call2);

