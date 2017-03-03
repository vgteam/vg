#include <iostream>
#include <vector>
#include <getopt.h>
#include <functional>
#include <regex>
//#include "intervaltree.hpp"
#include "subcommand.hpp"
#include "stream.hpp"
#include "index.hpp"
#include "position.hpp"
#include "vg.pb.h"
#include "path.hpp"
//#include "genotyper.hpp"
#include "genotypekit.hpp"
#include "path_index.hpp"
#include "vg.hpp"
#include "srpe.hpp"
#include "filter.hpp"
#include "utility.hpp"
#include "Variant.h"
#include "translator.hpp"
#include "Fasta.h"

using namespace std;
using namespace vg;
using namespace vg::subcommand;

#define DEBUG

void help_srpe(char** argv){
    cerr << "Usage: " << argv[0] << " srpe [options] <data.gam> <graph.vg>" << endl
        << "Options: " << endl
        << "Variant recall:" << endl
        << "   -S / --specific    <VCF>       look up variants in <VCF> in the graph and report only those." << endl
        << "   -R / --recall                  recall (i.e. type) all variant paths stored in the graph." << endl
        << "   -r / --reference   <REF>       reference genome to pull structural variants from." << endl
        << "   -I / --insertions  <INS>       fasta file containing insertion sequences." << endl
        << "Smart genotyping:" << endl
        << "   -a / --augmented   <AUG>       write the intermediate augmented graph to <AUG>." << endl
        << "   -p / --ref-path   <PATHNAME>   find variants relative to <PATHNAME>" << endl
        << endl;
    //<< "-S / --SV-TYPE comma separated list of SV types to detect (default: all)." << endl



}



int main_srpe(int argc, char** argv){
    string gam_name = "";
    string gam_index_name = "";
    string graph_name = "";
    string xg_name = "";
    string gcsa_name = "";
    string lcp_name = "";

    string spec_vcf = "";
    string ref_fasta = "";
    string ins_fasta = "";

    string augmented_graph_name = "";
    bool augment_paths = true;

    string ref_path = "";

    int max_iter = 2;
    int max_frag_len = 10000;
    int min_soft_clip = 12;

    bool do_all = false;

    vector<string> search_types;
    search_types.push_back("DEL");

    int threads = 1;

    if (argc <= 2) {
        help_srpe(argv);
        return 1;
    }

    int c;
    optind = 2; // force optind past command positional argument
    while (true) {
        static struct option long_options[] =
        {
            {"max-iter", required_argument, 0, 'm'},
            {"xg-index", required_argument, 0, 'x'},
            {"augmented", required_argument, 0, 'a'},
            {"help", no_argument, 0, 'h'},
            {"gcsa-index", required_argument, 0, 'g'},
            {"specific", required_argument, 0, 'S'},
            {"recall", no_argument, 0, 'R'},
            {"insertions", required_argument, 0, 'I'},
            {"reference", required_argument, 0, 'r'},
            {"threads", required_argument, 0, 't'},
            {"ref-path", required_argument, 0, 'p'},
            {0, 0, 0, 0}

        };
        int option_index = 0;
        c = getopt_long (argc, argv, "hx:g:m:S:RI:r:t:a:wp:",
                long_options, &option_index);

        // Detect the end of the options.
        if (c == -1)
            break;

        switch (c)
        {
            case 'a':
                augmented_graph_name = optarg;
                break;
            case 'm':
                max_iter = atoi(optarg);
                break;

            case 't':
                threads = atoi(optarg);
                break;

            case 'R':
                do_all = true;
                break;
            case 'x':
                xg_name = optarg;
                break;
            case 'g':
                gcsa_name = optarg;
                break;
            case 'S':
                spec_vcf = optarg;
                break;
            case 'r':
                ref_fasta = optarg;
                break;
            case 'I':
                ins_fasta = optarg;
                break;
            case 'p':
                ref_path = optarg;
                break;
            case 'h':
            case '?':
            default:
                help_srpe(argv);
                abort();
        }

    }

    omp_set_num_threads(threads);


    //SRPE srpe;


    gam_name = argv[optind];
    //gam_index_name = argv[++optind];
    graph_name = argv[++optind];

    xg::XG* xg_ind = new xg::XG();
    Index gamind;

    vg::VG* graph;

    // if (!xg_name.empty()){
    //     ifstream in(xg_name);
    //     xg_ind->load(in);
    // }
    // if (!gam_index_name.empty()){
    //     gamind.open_read_only(gam_index_name);
    // }
    // else{

    // }

    if (!graph_name.empty()){
        ifstream in(graph_name);
        graph = new VG(in, false);
    }

    // Open a variant call file,
    // hash each variant to an hash ID
    // have in if in the loop below.
    map<string, Locus> name_to_loc;
    // Makes a pathindex, which allows us to query length and push out a VCF with a position
    map<string, PathIndex*> pindexes;
    regex is_alt ("_alt_.*");

    vector<FastaReference*> insertions;
    if (!ins_fasta.empty()){
        FastaReference* ins = new FastaReference();
        insertions.emplace_back(ins);
        ins->open(ins_fasta);

    }

    if (!spec_vcf.empty() && ref_fasta.empty()){
        cerr << "Error: option -S requires a fasta reference using the -r <reference> flag" << endl;
    }
    else if (!spec_vcf.empty()){

        FastaReference* linear_ref = new FastaReference();
        linear_ref->open(ref_fasta);

        for (auto r_path : (graph->paths)._paths){
            if (!regex_match(r_path.first, is_alt)){
                pindexes[r_path.first] = new PathIndex(*graph, r_path.first, true);
            }
        }

        vcflib::VariantCallFile* variant_file = new vcflib::VariantCallFile();
        variant_file->open(spec_vcf);

        string descrip = "";
        descrip = "##INFO=<ID=AD,Number=R,Type=Integer,Description=\"Allele depth for each allele.\"\\>";
        variant_file->addHeaderLine(descrip);

        cout << variant_file->header << endl;

        unordered_map<string, list<Mapping> > graphpaths( (graph->paths)._paths.begin(), (graph->paths)._paths.end() );
        map<string, vcflib::Variant> hash_to_var;
        unordered_map<string, vector<int64_t> > varname_to_nodeid;
        unordered_map<int64_t, int32_t> node_to_depth;
        vector<int64_t> variant_nodes;
        variant_nodes.reserve(10000);

        // Hash a variant from the VCF
        vcflib::Variant var;
        while (variant_file->getNextVariant(var)){
            var.position -= 1;
            var.canonicalize_sv(*linear_ref, insertions, -1);
            string var_id = make_variant_id(var);
            hash_to_var[ var_id ] = var;
            for (int alt_ind = 0; alt_ind <= var.alt.size(); alt_ind++){
                string alt_id = "_alt_" + var_id + "_" + std::to_string(alt_ind);
                list<Mapping> x_path = graphpaths[ alt_id ];
                for (Mapping x_m : x_path){
                    variant_nodes.push_back(x_m.position().node_id());
                    varname_to_nodeid[ alt_id ].push_back(x_m.position().node_id());
                }
            }

        }
        std::function<void(const Alignment& a)> incr = [&](const Alignment& a){
            for (int i = 0; i < a.path().mapping_size(); i++){
                node_to_depth[ a.path().mapping(i).position().node_id() ] += 1;
            }
        };

        gamind.for_alignment_to_nodes(variant_nodes, incr);

        cerr << node_to_depth.size () << " reads in count map" << endl;

        for (auto it : hash_to_var){
            for (int i = 0; i <= it.second.alt.size(); i++){
                int32_t sum_reads = 0;
                string alt_id = "_alt_" + it.first + "_" + std::to_string(i);
                for (int i = 0; i < varname_to_nodeid[ alt_id ].size(); i++){
                    sum_reads += node_to_depth[varname_to_nodeid[ alt_id ][i]];

                }
                #pragma omp critical
                it.second.info["AD"].push_back(std::to_string(sum_reads));
            }
            cout << it.second << endl;
        }

    }
    else if (do_all){
        vector<Support> supports;

        for (auto r_path : (graph->paths)._paths){
            if (!regex_match(r_path.first, is_alt)){
                pindexes[r_path.first] = new PathIndex(*graph, r_path.first, true);
            }
        }
        for (auto x_path : (graph->paths)._paths){
            cerr << x_path.first << endl;
            int32_t support = 0;
            if (regex_match(x_path.first, is_alt)){
                vector<Alignment> alns;
                vector<int64_t> var_node_ids;
                for (Mapping x_m : x_path.second){
                    var_node_ids.push_back(x_m.position().node_id()); 
                }

                std::function<void(const Alignment&)> incr = [&](const Alignment& a){
                    ++support;
                };
                gamind.for_alignment_to_nodes(var_node_ids, incr);
                cout << support << " reads support " << x_path.first << endl;
            }
        }
    }
    else{

        // First, slurp in our discordant, split-read,
        // one-end-anchored, etc reads.
        // basename = * 
        // *.gam.split
        // *.gam.oea
        // *.gam.unmapped
        // *.gam.discordant
        
        ifstream all_reads;
        all_reads.open(gam_name);
        // Reads with no mismatches.
        vector<Alignment> perfects;
        // Reads with anchored mismatches or indels.
        vector<Alignment> simple_mismatches;
        // All the other reads
        vector<Alignment> complex_reads;

        Filter ff;
        vector<Path> simple_paths;

        // TODO expose at CLI
        int min_depth = 5;

        std::function<void(Alignment&)> get_simples_and_perfects = [&](Alignment& aln){
            if (ff.perfect_filter(aln)){
                #pragma omp critical
                perfects.push_back(aln);
            }
            else if (ff.simple_filter(aln)){
                #pragma omp critical
                {
                    simple_mismatches.push_back(aln);
                    Path p = trim_hanging_ends(aln.path());
                    p.set_name(aln.name());
                    simple_paths.push_back(p);
                    
                }
            }
            else{
                #pragma omp critical
                {
                complex_reads.push_back(aln);
                Path p = trim_hanging_ends(aln.path());
                p.set_name(aln.name());
                simple_paths.push_back(p);
                }
                
            }
        };

        stream::for_each_parallel(all_reads, get_simples_and_perfects);

        if (true){
            cerr << perfects.size() << " perfect reads." << endl;
            cerr << simple_mismatches.size() << " simple reads." << endl; 
            cerr << "Addding " << simple_paths.size() << " paths to graph." << endl;
            cerr << complex_reads.size() << " more complex reads." << endl;
        }

        // Grab our perfect and simple reads
        // incorporate our simple reads.
        // optionally, incorporate all reads and their paths.
        Translator tt;
        vector<Translation> translations = graph->edit(simple_paths);
        tt.load(translations);
        graph->paths.rebuild_mapping_aux();
        cerr << "Loaded " << translations.size() << " translations." << endl;

        if (augment_paths){

        }


        
        
        // Add both simple and perfect reads to the depth map.
        DepthMap dm(graph->size());
        for (auto x : complex_reads){
            list<Mapping>& mappings = graph->paths.get_path(x.name());
            for (auto m : mappings){
                dm.fill(m);
            }
        }
        cerr << "Filled depth map." << endl;
        


        
        // Save memory (maybe???) by disposing of our perfects and simples.
        // TODO is clear enough?
        //simple_mismatches.clear();
        //perfects.clear();

        // Rip out any nodes supported by an insufficient number of reads.
        std::function<void(Node*)> remove_low_depth_nodes = [&](Node* n){
            if (dm.get_depth(n->id()) < min_depth){
                #pragma omp critical
                {
                    cerr << "Destroying node: " << n->id() << endl;
                    graph->destroy_node(n);
                }
            }
        };
        graph->for_each_node(remove_low_depth_nodes);
        //graph->remove_orphan_edges();
        cerr << "Low-depth nodes/edges removed from graph" << endl;

        if (!augmented_graph_name.empty()){
            ofstream augstream;
            augstream.open(augmented_graph_name);
            graph->serialize_to_ostream(augstream);
        }
        // You could call SNPs right here. We should also check known variants at this stage.
        // Call snarlmanager, tossing out trivial snarls
        // Report them as either loci or VCF records
        SnarlFinder* snf = new CactusUltrabubbleFinder(*graph, ref_path, true);

        SnarlManager snarl_manager = snf->find_snarls();
        vector<const Snarl*> snarl_roots = snarl_manager.top_level_snarls();
        TraversalFinder* trav_finder;
        if (augment_paths){

        }
        else{
            trav_finder = new ExhaustiveTraversalFinder(*graph, snarl_manager);

        }

        vector<vector<Locus> > locus_buffer;
        int thread_count = omp_get_num_threads();
        locus_buffer.resize(thread_count);

        vcflib::VariantCallFile* vcf = nullptr;
        vector<PathIndex*> pindexes;


        bool justSNPs = false;
        if (justSNPs){

            for (auto x : snarl_roots){
                vector<SnarlTraversal> site_traversals = trav_finder->find_traversals(*x);
                for (auto t : site_traversals){
                    //cout << t.start() << " ";
                    for (auto zed : t.visits()){
                        //cout << zed.node_id() << ", ";
                    }
                    //cout <<  endl;
                }
            }
        }

        // Read in all our nasty discordant/split/clipped reads and
        // filter the nice ones from the nasty ones.

        double expected_insert_size = 500.0;
        double expected_insert_sd = 100.0;

        std::unordered_map<string, Alignment> mates_by_name;
        vector<Alignment> complex_singles;
        vector<Alignment> complex_pairs;

        // for every read in our nasty set, try to normalize it. If we succeed,
        // grab its mate (if appropriate) and generate an IMPRECISE variant based on the two.
        // Place that IMPRECISE candidate somewhere within its predicted range in the graph.

        // If desired, remap our unmapped and nasty reads. This will help us to augment our depth map.


        // call variants 
    
    }


    return 0;
}

static Subcommand vg_srpe ("srpe", "graph-external SV detection", main_srpe);

