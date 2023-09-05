/**
 * \file mcmc_main.cpp: GFA (Graph Alignment Format) Fast Emitter: a new mapper that will be *extremely* fast once we actually write it
 */

#include <omp.h>
#include <unistd.h>
#include <getopt.h>
#include <iostream>
#include <typeinfo>
#include <cassert>
#include <vector>
#include "subcommand.hpp"
#include <vg/io/vpkg.hpp>
#include "../mcmc_genotyper.hpp"
#include "../vg.hpp"
#include "../multipath_alignment.hpp"
#include "../mcmc_caller.hpp"
#include "../graph_caller.hpp"
#include <vg/io/stream.hpp>
#include <bdsg/overlays/overlay_helper.hpp>
#include <fstream> 
#include <iostream>

using namespace std;
using namespace vg;
using namespace vg::subcommand;

void help_mcmc(char** argv) {
    cerr
    << "usage: " << argv[0] << " mcmc [options] multipath_alns.mgam graph.vg sites.snarls > graph_with_paths.vg" << endl
    << "Finds haplotypes based on reads using MCMC methods" << endl
    << endl
    << "basic options:" << endl
    << "  -i, --iteration-number INT        tells us the number of iterations to run mcmc_genotyper with" <<endl
    << "  -r, --seed INT                    the seed we will use for the random number generator " << endl
    << "  -s, --sample NAME                 sample name [default=SAMPLE]" << endl
    << "  -p  --ref-path NAME               reference path to call on (multipile allowed.  defaults to all paths)"<< endl
    << "  -o, --ref-offset N                offset in reference path (multiple allowed, 1 per path)" << endl
    << "  -l, --ref-length N                override length of reference in the contig field of output VCF" << endl
    << "  -v, --vcf-out FILE                write VCF output to this file" << endl
    << "  -b, --burn-in INT                 number of iterations to run original sample proposal only" <<endl
    << "  -g, --gamma-freq INT              the frequency (every n iterations) for which to re-make the gamma set (starts after burn-in)" <<endl;
}

int main_mcmc(int argc, char** argv) {

    vector<string> ref_paths;
    vector<size_t> ref_path_offsets;
    vector<size_t> ref_path_lengths;
    
    string vcf_out;
    int burn_in;
    int gamma_freq;

    if (argc < 7) {
        help_mcmc(argv);
        return 1;
    }

    // initialize parameters with their default options
    int n_iterations = 1000;
    int seed = std::chrono::system_clock::now().time_since_epoch().count();
    string sample_name = "SAMPLE";
    

    int c;
    optind = 2; // force optind past command positional argument
    while (true) {
        static struct option long_options[] =
        {
            {"help", no_argument, 0, 'h'},
            {"iteration-number", required_argument, 0, 'i'},
            {"seed", required_argument, 0, 'r'},
            {"sample", required_argument, 0, 's'}, 
            {"ref-path", required_argument, 0, 'p'},
            {"ref-offset", required_argument, 0, 'o'},
            {"ref-length", required_argument, 0, 'l'}, 
            {"vcf-out", required_argument, 0, 'v'},
            {"burn-in", required_argument, 0, 'b'},
            {"gamma-freq", required_argument, 0, 'g'},
            {0, 0, 0, 0}
        };

        int option_index = 0;
        c = getopt_long (argc, argv, "hi:s:p:o:l:r:v:b:g:",
                         long_options, &option_index);


        // Detect the end of the options.
        if (c == -1)
            break;

        switch (c)
        {
            case 'i':
                n_iterations = parse<int>(optarg);
                break;
            case 'r':
                seed = parse<int>(optarg);
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
            case 's':
                sample_name = optarg;
                break;  
            case 'v':
                vcf_out = optarg;
                break;
            case 'b':
                burn_in = parse<int>(optarg);
                break;
            case 'g':
                gamma_freq = parse<int>(optarg);
                break;
            case 'h':
            case '?':
            default:
                help_mcmc(argv);
                exit(1);
                break;
        }
    }

    string multipath_file =  get_input_file_name(optind, argc, argv);
    string graph_file =  get_input_file_name(optind, argc, argv);
    string snarls_file =  get_input_file_name(optind, argc, argv);
   
    unique_ptr<SnarlManager> snarls = (vg::io::VPKG::load_one<SnarlManager>(snarls_file));

    // // create a PathHandleGraph 
    unique_ptr<PathHandleGraph> path_hgraph;
    bdsg::PathPositionOverlayHelper overlay_helper;
    path_hgraph = vg::io::VPKG::load_one<PathHandleGraph>(graph_file);

    // Some stuff below here needs a vg graph.
    VG* vg_graph = dynamic_cast<vg::VG*>(path_hgraph.get());
    
    // Call this to populate the vg_graph if it isn't populated.
    auto ensure_vg = [&]() -> vg::VG* {
        if (vg_graph == nullptr) {
            // Copy instead.
            vg_graph = new vg::VG();
            handlealgs::copy_path_handle_graph(path_hgraph.get(), vg_graph);
            // Give the unique_ptr ownership and delete the graph we loaded.
            path_hgraph.reset(vg_graph);
            // Make sure the paths are all synced up
            vg_graph->paths.to_graph(vg_graph->graph);
        }
        return vg_graph;
    };
    
    //convert to VG graph if needed
    ensure_vg();

    if(vg_graph == nullptr || vg_graph == 0){
        cerr << "Graph is NULL" <<endl;
        exit(1);
    }
    PathPositionHandleGraph* graph = nullptr;
    graph = overlay_helper.apply(vg_graph);
    
     
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
    // Check our ref lengths
    if (ref_path_lengths.size() != 0 && ref_path_lengths.size() != ref_paths.size()) {
        cerr << "error [vg call]: when using -l, the same number paths must be given with -p" << endl;
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
    
    // Check if VCF output file is specified 
    ofstream vcf_file_out;
    if(!vcf_out.empty()){
        vcf_file_out.open(vcf_out, ios::out);
    }
    
      /*
    *########################################################################################
    *                      GENOTYPING
    *########################################################################################
    **/

    vector<multipath_alignment_t> reads;
    get_input_file(multipath_file, [&] (istream& open_file){
        io::ProtobufIterator<MultipathAlignment> iter (open_file);
        while(iter.has_current()){
            reads.emplace_back();
            from_proto_multipath_alignment(*iter, reads.back());
            // vg::view_multipath_alignment_as_dot(cerr,*iter,true);
            ++iter;
        }
    });
    double log_base = gssw_dna_recover_log_base(1,4,.5,1e-12);
    // invoke run genotyper 
    MCMCGenotyper mcmc_genotyper(*snarls, *vg_graph, n_iterations, seed, burn_in, gamma_freq);
    unique_ptr<PhasedGenome> genome = mcmc_genotyper.run_genotype(reads, log_base);
    
    // genome->print_phased_genome();

    /*
    *########################################################################################
    *                      VCF OUTPUT
    *########################################################################################
    **/
    
    // Create MCMC_Caller object
    MCMCCaller mcmc_caller(graph, *genome, *snarls, sample_name, ref_paths, ref_path_offsets, ref_path_lengths, cout);

    // Write header to ofstream  
    vcf_file_out << mcmc_caller.vcf_header(*graph, ref_paths, ref_path_lengths);
    
    //current implimentation is writing vcf record after each variant processed
    mcmc_caller.call_top_level_snarls();

    // mcmc_caller.write_variants(cerr);
    mcmc_caller.write_variants(vcf_file_out);
    
    //close the vcf file
    vcf_file_out.close();
   
    // will output a graph w/ embedded paths
    vg_graph->serialize_to_ostream(std::cout);

    return 0;
}

// Register subcommand
static Subcommand vg_mcmc("mcmc", "Finds haplotypes based on reads using MCMC methods", DEVELOPMENT, main_mcmc);


