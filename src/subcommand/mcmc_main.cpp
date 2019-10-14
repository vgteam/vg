/**
 * \file mcmc_main.cpp: GFA (Graph Alignment Format) Fast Emitter: a new mapper that will be *extremely* fast once we actually write it
 */

#include <omp.h>
#include <unistd.h>
#include <getopt.h>
#include <iostream>
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
    << "  -s, --seed INT                    the seed we will use for the random number generator " << endl
    << "  -p  --ref-path NAME               reference path to call on (multipile allowed.  defaults to all paths)"<< endl
    << "  -o, --ref-offset N                offset in reference path (multiple allowed, 1 per path)" << endl
    << "  -l, --ref-length N                override length of reference in the contig field of output VCF" << endl;
}

int main_mcmc(int argc, char** argv) {

    // holds the ref-path name 
    vector<string> ref_paths;
    vector<size_t> ref_path_offsets;
    vector<size_t> ref_path_lengths;

    if (argc < 5) {
        help_mcmc(argv);
        return 1;
    }

    // initialize parameters with their default options
    int n_iterations = 1000;
    int seed = 1;
    
    int c;
    optind = 2; // force optind past command positional argument
    while (true) {
        static struct option long_options[] =
        {
            {"help", no_argument, 0, 'h'},
            {"iteration-number", required_argument, 0, 'i'},
            {"seed", required_argument, 0, 's'},
            {"ref-path", required_argument, 0, 'p'},
            {"ref-offset", required_argument, 0, 'o'},
            {"ref-length", required_argument, 0, 'l'}, 
            {0, 0, 0, 0}
        };

        int option_index = 0;
        c = getopt_long (argc, argv, "hi:s:p:",
                         long_options, &option_index);


        // Detect the end of the options.
        if (c == -1)
            break;

        switch (c)
        {
            case 'i':
                n_iterations = parse<int>(optarg);
                break;

            case 's':
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

    unique_ptr<VG> graph = (vg::io::VPKG::load_one<VG>(graph_file));
    unique_ptr<SnarlManager> snarls = (vg::io::VPKG::load_one<SnarlManager>(snarls_file));

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
    

    vector<MultipathAlignment> reads;
    get_input_file(multipath_file, [&] (istream& open_file){
        io::ProtobufIterator<MultipathAlignment> iter (open_file);
        while(iter.has_current()){
            reads.push_back(*iter);
            // vg::view_multipath_alignment_as_dot(cerr,*iter,true);
            ++iter;
        }
    });
    double log_base = gssw_dna_recover_log_base(1,4,.5,1e-12);
    // invoke run genotyper 
    MCMCGenotyper mcmc_genotyper(*snarls, *graph, n_iterations, seed);
    unique_ptr<PhasedGenome> genome = mcmc_genotyper.run_genotype(reads, log_base );
    
    
    for (int i = 0; i < 2; i++){
        
        // create two paths tracing each haplotype
        path_handle_t path_handle = graph->create_path_handle("H" + to_string(i));
        
        for(auto iter = genome->begin(i); iter!= genome->end(i); iter++){
        
        graph->append_step(path_handle, graph->get_handle((*iter).node->id()));
        
    }
    }

    unique_ptr<MCMCCaller> mcmc_caller; //TODO: initialize this to something.....
    VCFOutputCaller* vcf_caller = dynamic_cast<VCFOutputCaller*>(mcmc_caller.get());
    assert(vcf_caller != nullptr);
    cout << vcf_caller->vcf_header(*graph, ref_paths, ref_path_lengths) << flush;
    vcf_caller->write_variants(cout);
        
    return 0;


    
    // will output a graph w/ embedded paths
    graph->serialize_to_ostream(std::cout);
   
    return 0;
}

// Register subcommand
static Subcommand vg_mcmc("mcmc", "Finds haplotypes based on reads using MCMC methods", DEVELOPMENT, main_mcmc);


