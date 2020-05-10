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
    << "  -s, --seed INT                    the seed we will use for the random number generator " << endl;
}

int main_mcmc(int argc, char** argv) {


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
            {0, 0, 0, 0}
        };

        int option_index = 0;
        c = getopt_long (argc, argv, "hi:s:",
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
    MCMCGenotyper mcmc_genotyper(*snarls, *graph, n_iterations, seed);
    unique_ptr<PhasedGenome> genome = mcmc_genotyper.run_genotype(reads, log_base );
    
    
    for (int i = 0; i < 2; i++){
        
        // create two paths tracing each haplotype
        path_handle_t path_handle = graph->create_path_handle("H" + to_string(i));
        
        for(auto iter = genome->begin(i); iter!= genome->end(i); iter++){
        
        graph->append_step(path_handle, graph->get_handle((*iter).node->id()));
        
    }
    }
    
    
    
    // will output a graph w/ embedded paths
    graph->serialize_to_ostream(std::cout);
   
    return 0;
}

// Register subcommand
static Subcommand vg_mcmc("mcmc", "Finds haplotypes based on reads using MCMC methods", DEVELOPMENT, main_mcmc);


