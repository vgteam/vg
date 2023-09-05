
#include <getopt.h>
#include "subcommand.hpp"
#include <vg/io/stream.hpp>
#include "../genotyper.hpp"
#include "../genotypekit.hpp"
#include "../variant_recall.hpp"
/**
* GAM sort main
*/

using namespace std;
using namespace vg;
using namespace vg::subcommand;
void help_genotype(char** argv) {
    cerr << "usage: " << argv[0] << " genotype [options] <graph.vg> alignments.gam > <calls.vcf>" << endl
         << "Compute genotypes from a graph and a collection of reads" << endl
         << endl
         << "options:" << endl
         << "    -j, --json              output in JSON" << endl
         << "    -v, --vcf               output in VCF" << endl
         << "    -V, --recall-vcf VCF    recall variants in a specific VCF file." << endl
         << "    -F, --fasta  FASTA" << endl
         << "    -I, --insertions INS" << endl
         << "    -r, --ref PATH          use the given path name as the reference path" << endl
         << "    -c, --contig NAME       use the given name as the VCF contig name" << endl
         << "    -s, --sample NAME       name the sample in the VCF with the given name" << endl
         << "    -o, --offset INT        offset variant positions by this amount" << endl
         << "    -l, --length INT        override total sequence length" << endl
         << "    -a, --augmented FILE    dump augmented graph to FILE" << endl
         << "    -Q, --ignore_mapq       do not use mapping qualities" << endl
         << "    -A, --no_indel_realign  disable indel realignment" << endl
         << "    -d, --het_prior_denom   denominator for prior probability of heterozygousness" << endl
         << "    -P, --min_per_strand    min unique reads per strand for a called allele to accept a call" << endl
         << "    -E, --no_embed          don't embed gam edits into graph" << endl
         << "    -T, --traversal         traversal finder to use {reads, exhaustive, representative, adaptive} (adaptive)" << endl
         << "    -p, --progress          show progress" << endl
         << "    -t, --threads N         number of threads to use" << endl;
}

int main_genotype(int argc, char** argv) {

    if (argc <= 2) {
        help_genotype(argv);
        return 1;
    }
    // Should we output genotypes in JSON (true) or Protobuf (false)?
    bool output_json = false;
    // Should we output VCF instead of protobuf?
    bool output_vcf = false;
    // Should we show progress with a progress bar?
    bool show_progress = false;
    // How many threads should we use?
    int thread_count = 0;

    // What reference path should we use
    string ref_path_name;
    // What sample name should we use for output
    string sample_name;
    // What contig name override do we want?
    string contig_name;
    // What offset should we add to coordinates
    int64_t variant_offset = 0;
    // What length override should we use
    int64_t length_override = 0;
    // Should we embed gam edits (for debugging as we move to further decouple augmentation and calling)
    bool embed_gam_edits = true;
    // Which traversal finder should we use
    string traversal_finder = "adaptive";

    // Should we we just do a quick variant recall,
    // based on this VCF and GAM, then exit?
    string recall_vcf;
    string gam_file;
    string fasta;
    string insertions_file;

    // Should we use mapping qualities?
    bool use_mapq = true;
    // Should we do indel realignment?
    bool realign_indels = true;

    // Should we dump the augmented graph to a file?
    string augmented_file_name;
    
    // What should the heterozygous genotype prior be? (1/this)
    double het_prior_denominator = 10.0;
    // At least how many reads must be unique support for a called allele per strand for a call?
    size_t min_unique_per_strand = 2;

    bool just_call = false;
    int c;
    optind = 2; // force optind past command positional arguments
    while (true) {
        static struct option long_options[] =
            {
                {"json", no_argument, 0, 'j'},
                {"vcf", no_argument, 0, 'v'},
                {"ref", required_argument, 0, 'r'},
                {"contig", required_argument, 0, 'c'},
                {"sample", required_argument, 0, 's'},
                {"offset", required_argument, 0, 'o'},
                {"length", required_argument, 0, 'l'},
                {"augmented", required_argument, 0, 'a'},
                {"ignore_mapq", no_argument, 0, 'Q'},
                {"no_indel_realign", no_argument, 0, 'A'},
                {"het_prior_denom", required_argument, 0, 'd'},
                {"min_per_strand", required_argument, 0, 'P'},
                {"progress", no_argument, 0, 'p'},
                {"threads", required_argument, 0, 't'},
                {"recall-vcf", required_argument, 0, 'V'},
                {"fasta", required_argument, 0, 'F'},
                {"insertions", required_argument, 0, 'I'},
                {"call", no_argument, 0, 'z'},
                {"no_embed", no_argument, 0, 'E'},
                {"traversal", required_argument, 0, 'T'},
                {0, 0, 0, 0}
            };

        int option_index = 0;
        c = getopt_long (argc, argv, "hjvr:c:s:o:l:a:QAd:P:pt:V:I:F:zET:",
                         long_options, &option_index);

        /* Detect the end of the options. */
        if (c == -1)
            break;

        switch (c)
        {
        case 'j':
            output_json = true;
            break;
        case 'v':
            output_vcf = true;
            break;
        case 'r':
            // Set the reference path name
            ref_path_name = optarg;
            break;
        case 'c':
            // Set the contig name for output
            contig_name = optarg;
            break;
        case 's':
            // Set the sample name
            sample_name = optarg;
            break;
        case 'o':
            // Offset variants
            variant_offset = parse<int64_t>(optarg);
            break;
        case 'l':
            // Set a length override
            length_override = parse<int64_t>(optarg);
            break;
        case 'a':
            // Dump augmented graph
            augmented_file_name = optarg;
            break;
        case 'Q':
            // Ignore mapping qualities
            use_mapq = false;
            break;
        case 'z':
            just_call = true;
            break;
        case 'A':
            // Don't do indel realignment
            realign_indels = false;
            break;
        case 'd':
            // Set heterozygous genotype prior denominator
            het_prior_denominator = parse<double>(optarg);
            break;
        case 'P':
            // Set min consistent reads per strand required to keep an allele
            min_unique_per_strand = parse<int64_t>(optarg);
            break;
        case 'p':
            show_progress = true;
            break;
        case 't':
            thread_count = parse<int>(optarg);
            break;
        case 'V':
            recall_vcf = optarg;
            break;
        case 'I':
            insertions_file = optarg;
            break;
        case 'F':
            fasta = optarg;
            break;
        case 'E':
            embed_gam_edits = false;
            break;
        case 'T':
            traversal_finder = optarg;
            break;
        case 'h':
        case '?':
            /* getopt_long already printed an error message. */
            help_genotype(argv);
            exit(1);
            break;
        default:
          abort ();
        }
    }

    if(thread_count > 0) {
        omp_set_num_threads(thread_count);
    }

    // read the graph
    if (optind >= argc) {
        help_genotype(argv);
        return 1;
    }

    if (show_progress) {
        cerr << "Reading input graph..." << endl;
    }
    VG* graph;
    get_input_file(optind, argc, argv, [&](istream& in) {
        graph = new VG(in);
    });
    
    // get GAM
    if (optind < argc){
        gam_file = get_input_file_name(optind, argc, argv);
        
    } else {
        cerr << "[vg genotype] GAM file must be specified as positional argument" << endl;
        return 1;
    }

    if (just_call){
        string gamfi(gam_file);
        string rstr(ref_path_name);
        genotype_svs(graph, gamfi, rstr);
        exit(0);
    }

    // Build the set of all the node IDs to operate on
    vector<vg::id_t> graph_ids;
    graph->for_each_node([&](Node* node) {
        // Put all the ids in the set
        graph_ids.push_back(node->id());
    });

    if (!(recall_vcf.empty() || fasta.empty())){
        vcflib::VariantCallFile* vars = new vcflib::VariantCallFile();
        vars->open(recall_vcf);
        FastaReference* lin_ref = new FastaReference();
        lin_ref->open(fasta);

        vector<FastaReference*> insertions;
        if (!insertions_file.empty()){
            FastaReference* ins = new FastaReference();
            insertions.emplace_back(ins);
            ins->open(insertions_file);
        }
        variant_recall(graph, vars, lin_ref, insertions, gam_file);
        return 0;

    }

    // Load all the reads matching the graph into memory
    vector<Alignment> alignments;

    if(show_progress) {
        cerr << "Loading reads..." << endl;
    }

    function<bool(const Alignment&)> alignment_contained = [&graph](const Alignment& alignment) {
        for(size_t i = 0; i < alignment.path().mapping_size(); i++) {
            if(!graph->has_node(alignment.path().mapping(i).position().node_id())) {
                return false;
            }
        }
        return alignment.path().mapping_size() > 0;
    };

    // load in all reads (activated by passing GAM directly with -G).
    // This is used by, ex., toil-vg, which has already used the gam index
    // to extract relevant reads
    ifstream gam_reads(gam_file.c_str());
    if (!gam_reads) {
        cerr << "[vg genotype] Error opening gam: " << gam_file << endl;
        return 1;
    }
    vg::io::for_each<Alignment>(gam_reads, [&alignments, &alignment_contained](Alignment& alignment) {
        if (alignment_contained(alignment)) {
            alignments.push_back(alignment);
        }
    });
    
    if(show_progress) {
        cerr << "Loaded " << alignments.size() << " alignments" << endl;
    }
    
    // Make a Genotyper to do the genotyping
    Genotyper genotyper;
    // Configure it
    genotyper.use_mapq = use_mapq;
    genotyper.realign_indels = realign_indels;
    assert(het_prior_denominator > 0);
    genotyper.het_prior_logprob = prob_to_logprob(1.0/het_prior_denominator);
    genotyper.min_unique_per_strand = min_unique_per_strand;
    if (traversal_finder == "reads") {
        genotyper.traversal_alg = Genotyper::TraversalAlg::Reads;
    } else if (traversal_finder == "exhaustive") {
        genotyper.traversal_alg = Genotyper::TraversalAlg::Exhaustive;
    } else if (traversal_finder == "representative") {
        genotyper.traversal_alg = Genotyper::TraversalAlg::Representative;
    } else if (traversal_finder == "adaptive") {
      genotyper.traversal_alg = Genotyper::TraversalAlg::Adaptive;
    } else {
        cerr << "Invalid value for traversal finder: " << traversal_finder
             << ".  Must be in {reads, representative, exhaustive, adaptive}" << endl;
        return 1;
    }
    genotyper.show_progress = show_progress;

    // Guess the reference path if not given
    if(ref_path_name.empty()) {
        // Guess the ref path name
        if(graph->paths.size() == 1) {
            // Autodetect the reference path name as the name of the only path
            ref_path_name = (*graph->paths._paths.begin()).first;
        } else {
            ref_path_name = "ref";
        }
    }

    // Augment the graph with all the reads
    AugmentedGraph augmented_graph;

    // Move our input graph into the augmented graph
    swap(augmented_graph.graph, *graph); 

    // Do the actual augmentation using vg edit. If augmentation was already
    // done, just embeds the reads. Reads will be taken by the AugmentedGraph
    // and stored in it.
    augmented_graph.augment_from_alignment_edits(alignments, true, !embed_gam_edits);
    
    // TODO: move arguments below up into configuration
    genotyper.run(augmented_graph,
                  cout,
                  ref_path_name,
                  contig_name,
                  sample_name,
                  augmented_file_name,
                  output_vcf,
                  output_json,
                  length_override,
                  variant_offset);

    delete graph;

    return 0;
}


static Subcommand vg_genotype("genotype", "Genotype (or type) graphs, GAMS, and VCFs.", main_genotype);
