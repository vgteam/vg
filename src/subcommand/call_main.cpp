/** \file call_main.cpp
 *
 * Defines the "vg call" subcommand, which calls variation from a graph and a pileup.
 */

#include <omp.h>
#include <unistd.h>
#include <getopt.h>

#include <list>
#include <fstream>

#include "subcommand.hpp"

#include "../vg.hpp"
#include "../caller.hpp"



using namespace std;
using namespace vg;
using namespace vg::subcommand;

void help_call(char** argv) {
    cerr << "usage: " << argv[0] << " call [options] <graph.vg> <pileup.vgpu> > output.vcf" << endl
         << "Output variant calls in VCF or Loci format given a graph and pileup" << endl
         << endl
         << "options:" << endl
         << "    -V, --no-vcf               output variants in binary Loci format instead of text VCF format" << endl
         << "    -d, --min-depth INT        minimum depth of pileup [" << Caller::Default_min_depth <<"]" << endl
         << "    -e, --max-depth INT        maximum depth of pileup [" << Caller::Default_max_depth <<"]" << endl
         << "    -s, --min-support INT      minimum number of reads required to support snp [" << Caller::Default_min_support <<"]" << endl
         << "    -f, --min-frac FLOAT       minimum percentage of reads required to support snp[" << Caller::Default_min_frac <<"]" << endl
         << "    -q, --default-read-qual N  phred quality score to use if none found in the pileup ["
         << (int)Caller::Default_default_quality << "]" << endl
         << "    -b, --max-strand-bias FLOAT limit to absolute difference between 0.5 and proportion of supporting reads on reverse strand. [" << Caller::Default_max_strand_bias << "]" << endl
         << "    -a, --link-alts            add all possible edges between adjacent alts" << endl
         << "    -A, --aug-graph FILE       write out the agumented graph in vg format" << endl
         << "    -r, --ref PATH             use the given path name as the reference path" << endl
         << "    -c, --contig NAME          use the given name as the VCF contig name" << endl
         << "    -S, --sample NAME          name the sample in the VCF with the given name [SAMPLE]" << endl
         << "    -o, --offset INT           offset variant positions by this amount in VCF [0]" << endl
         << "    -l, --length INT           override total sequence length in VCF" << endl
         << "    -U, --subgraph             expect a subgraph and ignore extra pileup entries outside it" << endl
         << "    -P, --pileup               write pileup under VCF lines (for debugging, output not valid VCF)" << endl
         << "    -D, --depth INT            maximum depth for path search [default 10 nodes]" << endl
         << "    -F, --min-cov-frac FLOAT   min fraction of average coverage at which to call [0.0]" << endl
         << "    -H, --max-het-bias FLOAT   max imbalance factor between alts to call heterozygous [3]" << endl
         << "    -R, --max-ref-bias FLOAT   max imbalance factor between ref and alts to call heterozygous ref [4]" << endl
         << "    -M, --bias-mult FLOAT      multiplier for bias limits for indels as opposed to substitutions [1]" << endl
         << "    -n, --min-count INT        min total supporting read count to call a variant [1]" << endl
         << "    -B, --bin-size  INT        bin size used for counting coverage [250]" << endl
         << "    -C, --exp-coverage INT     specify expected coverage (instead of computing on reference)" << endl
         << "    -O, --no-overlap           don't emit new variants that overlap old ones" << endl
         << "    -u, --use-avg-support      use average instead of minimum support" << endl
         << "    -E, --min_mad              min. minimum allele depth required to PASS filter [5]" << endl
         << "    -h, --help                 print this help message" << endl
         << "    -p, --progress             show progress" << endl
         << "    -v, --verbose              print information and warnings about vcf generation" << endl
         << "    -t, --threads N            number of threads to use" << endl;
}

int main_call(int argc, char** argv) {

    if (argc <= 3) {
        help_call(argv);
        return 1;
    }

    double het_prior = Caller::Default_het_prior;
    int min_depth = Caller::Default_min_depth;
    int max_depth = Caller::Default_max_depth;
    int min_support = Caller::Default_min_support;
    double min_frac = Caller::Default_min_frac;
    int default_read_qual = Caller::Default_default_quality;
    double max_strand_bias = Caller::Default_max_strand_bias;
    string aug_file;
    bool bridge_alts = false;
    
    
    
    // Should we expect a subgraph and ignore pileups for missing nodes/edges?
    bool expectSubgraph = false;
    
    // Should we annotate the VCF with pileup info?
    bool pileupAnnotate = false;

    // This manages conversion from an augmented graph to a VCF, and makes the
    // actual calls.
    Call2Vcf call2vcf;

    bool show_progress = false;
    int thread_count = 1;

    int c;
    optind = 2; // force optind past command positional arguments
    while (true) {
        static struct option long_options[] =
            {
                {"no-vcf", no_argument, 0, 'V'},
                {"min-depth", required_argument, 0, 'd'},
                {"max-depth", required_argument, 0, 'e'},
                {"min-support", required_argument, 0, 's'},
                {"min-frac", required_argument, 0, 'f'},
                {"default-read-qual", required_argument, 0, 'q'},
                {"max-strand-bias", required_argument, 0, 'b'},
                {"aug-graph", required_argument, 0, 'A'},
                {"link-alts", no_argument, 0, 'a'},
                {"progress", no_argument, 0, 'p'},
                {"verbose", no_argument, 0, 'v'},
                {"threads", required_argument, 0, 't'},
                {"ref", required_argument, 0, 'r'},
                {"contig", required_argument, 0, 'c'},
                {"sample", required_argument, 0, 'S'},
                {"offset", required_argument, 0, 'o'},
                {"depth", required_argument, 0, 'D'},
                {"length", required_argument, 0, 'l'},
                {"subgraph", no_argument, 0, 'U'},
                {"pileup", no_argument, 0, 'P'},
                {"min-cov-frac", required_argument, 0, 'F'},
                {"max-het-bias", required_argument, 0, 'H'},
                {"max-ref-bias", required_argument, 0, 'R'},
                {"bias-mult", required_argument, 0, 'M'},
                {"min-count", required_argument, 0, 'n'},
                {"bin-size", required_argument, 0, 'B'},
                {"avg-coverage", required_argument, 0, 'C'},
                {"no-overlap", no_argument, 0, 'O'},
                {"use-avg-support", no_argument, 0, 'u'},
                {"min-mad", required_argument, 0, 'E'},
                {"help", no_argument, 0, 'h'},
                {0, 0, 0, 0}
            };

        int option_index = 0;
        c = getopt_long (argc, argv, "Vd:e:s:f:q:b:A:apvt:r:c:S:o:D:l:UPF:H:R:M:n:B:C:OuE:h",
                         long_options, &option_index);

        /* Detect the end of the options. */
        if (c == -1)
            break;

        switch (c)
        {
        case 'V':
            call2vcf.convert_to_vcf = false;
            break;
        case 'd':
            min_depth = atoi(optarg);
            break;
        case 'e':
            max_depth = atoi(optarg);
            break;
        case 's':
            min_support = atoi(optarg);
            break;
        case 'f':
            min_frac = atof(optarg);
            break;
        case 'q':
            default_read_qual = atoi(optarg);
            break;
        case 'b':
            max_strand_bias = atof(optarg);
            break;
        case 'A':
            aug_file = optarg;
            break;
        case 'a':
            bridge_alts = true;
            break;
        // old glenn2vcf opts start here
        case 'r':
            // Set the reference path name
            call2vcf.refPathName = optarg;
            break;
        case 'c':
            // Set the contig name
            call2vcf.contigName = optarg;
            break;
        case 'S':
            // Set the sample name
            call2vcf.sampleName = optarg;
            break;
        case 'o':
            // Offset variants
            call2vcf.variantOffset = std::stoll(optarg);
            break;
        case 'D':
            // Limit max depth for pathing to primary path
            call2vcf.maxDepth = std::stoll(optarg);
            break;
        case 'l':
            // Set a length override
            call2vcf.lengthOverride = std::stoll(optarg);
            break;
        case 'U':
            expectSubgraph = true;
            break;
        case 'P':
            pileupAnnotate = true;
            break;
        case 'F':
            // Set min fraction of average coverage for a call
            call2vcf.minFractionForCall = std::stod(optarg);
            break;
        case 'H':
            // Set max factor between reads on one alt and reads on the other
            // alt for calling a het.
            call2vcf.maxHetBias = std::stod(optarg);
            break;
        case 'R':
            // Set max factor between reads on ref and reads on the other
            // alt for calling a homo ref.
            call2vcf.maxRefHetBias = std::stod(optarg);
            break;
        case 'M':
            // Set multiplier for bias limits for indels
            call2vcf.indelBiasMultiple = std::stod(optarg);
            break;
        case 'n':
            // How many reads need to touch an allele before we are willing to
            // call it?
            call2vcf.minTotalSupportForCall = std::stoll(optarg);
            break;
        case 'B':
            // Set the reference bin size
            call2vcf.refBinSize = std::stoll(optarg);
            break;
        case 'C':
            // Override expected coverage
            call2vcf.expCoverage = std::stoll(optarg);
            break;
        case 'O':
            // Suppress variants that overlap others
            call2vcf.suppress_overlaps = true;
            break;
        case 'u':
            // Average (isntead of min) support
            call2vcf.useAverageSupport = true;
            break;
        case 'E':
            // Minimum min-allele-depth required to give Filter column a PASS
            call2vcf.min_mad_for_filter = std::stoi(optarg);
            break;
        case 'p':
            show_progress = true;
            break;
        case 'v':
            call2vcf.verbose = true;
            break;
        case 't':
            thread_count = atoi(optarg);
            break;
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
    omp_set_num_threads(thread_count);
    thread_count = get_thread_count();

    // Parse the arguments
    if (optind >= argc) {
        help_call(argv);
        return 1;
    }
    string graph_file_name = get_input_file_name(optind, argc, argv);
    if (optind >= argc) {
        help_call(argv);
        return 1;
    }
    string pileup_file_name = get_input_file_name(optind, argc, argv);
    
    if (pileup_file_name == "-" && graph_file_name == "-") {
        cerr << "error: graph and pileup can't both be from stdin." << endl;
        exit(1);
    }
    
    // read the graph
    if (show_progress) {
        cerr << "Reading input graph" << endl;
    }
    VG* graph;
    get_input_file(graph_file_name, [&](istream& in) {
        graph = new VG(in);
    });

    if (show_progress) {
        cerr << "Computing augmented graph" << endl;
    }
    Caller caller(graph,
                  het_prior, min_depth, max_depth, min_support,
                  min_frac, Caller::Default_min_log_likelihood,
                  default_read_qual, max_strand_bias, bridge_alts);

    // setup pileup stream
    get_input_file(pileup_file_name, [&](istream& pileup_stream) {
        // compute the augmented graph
        function<void(Pileup&)> lambda = [&](Pileup& pileup) {
            for (int i = 0; i < pileup.node_pileups_size(); ++i) {
                if (!graph->has_node(pileup.node_pileups(i).node_id())) {
                    // This pileup doesn't belong in this graph
                    if(!expectSubgraph) {
                        throw runtime_error("Found pileup for nonexistent node " + to_string(pileup.node_pileups(i).node_id()));
                    }
                    // If that's expected, just skip it
                    continue;
                }
                // Send approved pileups to the caller
                caller.call_node_pileup(pileup.node_pileups(i));
            }
            for (int i = 0; i < pileup.edge_pileups_size(); ++i) {
                if (!graph->has_edge(pileup.edge_pileups(i).edge())) {
                    // This pileup doesn't belong in this graph
                    if(!expectSubgraph) {
                        throw runtime_error("Found pileup for nonexistent edge " + pb2json(pileup.edge_pileups(i).edge()));
                    }
                    // If that's expected, just skip it
                    continue;
                }
                // Send approved pileups to the caller
                caller.call_edge_pileup(pileup.edge_pileups(i));
            }
        };
        stream::for_each(pileup_stream, lambda);
    });
    
    // map the edges from original graph
    if (show_progress) {
        cerr << "Mapping edges into augmented graph" << endl;
    }
    caller.update_augmented_graph();

    // map the paths from the original graph
    if (show_progress) {
        cerr << "Mapping paths into augmented graph" << endl;
    }
    caller.map_paths();

    if (!aug_file.empty()) {
        // write the augmented graph
        if (show_progress) {
            cerr << "Writing augmented graph" << endl;
        }
        ofstream aug_stream(aug_file.c_str());
        caller.write_augmented_graph(aug_stream, false);
    }
    
    if (show_progress) {
        cerr << "Calling variants" << endl;
    }

    // project the augmented graph to a reference path
    // in order to create a VCF of calls.
    call2vcf.call(caller._augmented_graph,
        pileupAnnotate ? pileup_file_name : string());

    return 0;
}

// Register subcommand
static Subcommand vg_call("call", "call variants on a graph from a pileup", main_call);

