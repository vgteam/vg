#pragma once // TODO: remove this, to avoid warnings + maybe bad coding practice?
#include <getopt.h>
#include <omp.h>
#include <unistd.h>

#include <gbwtgraph/gbwtgraph.h>

#include "subcommand.hpp"

#include "../../include/bdsg/hash_graph.hpp"
#include "../../include/vg/io/vpkg.hpp"
// #include "../algorithms/0_draft_haplotype_realignment.hpp"
#include "../algorithms/0_oo_normalize_snarls.hpp"
#include "../gbwt_helper.hpp"

#include "../../include/vg/io/vpkg.hpp"

using namespace std;
using namespace vg;
using namespace vg::subcommand;

void test(const HandleGraph &graph){
    handle_t handle = graph.get_handle(2);
    cerr << "seq of handle 2 " << graph.get_sequence(handle);
}

void help_vpkg_wrap(char **argv) {
    cerr
        << "usage: " << argv[0] << " vpkg_wrap [options] <graph.vg> >[vpkg_wrapped.hg]" << endl
        << "Wraps given vg graph into vpkg format, saves as hg." << endl
        << endl
        << "options:" << endl;
        // << "    -v, --input_vg       unwrapped vg." << endl;
}


int main_vpkg_wrap(int argc, char **argv) {

    if (argc == 2) {
        help_vpkg_wrap(argv);
        return 1;
    }

    string input_vg;

    int c;
    optind = 2; // force optind past command positional argument
    while (true) {
        static struct option long_options[] =

            {{"help", no_argument, 0, 'h'},
             {"input_vg", required_argument, 0, 'v'},
             {0, 0, 0, 0}};

        int option_index = 0;
        c = getopt_long(argc, argv, "hv:", long_options, &option_index);

        // Detect the end of the options.
        if (c == -1)
            break;

        switch (c) {

        case 'v':
            input_vg = optarg;
            break;

        default:
            abort();
        }
    }

    // if (input_vg.size() > 0) {
        // //get a HashGraph graph:
        // cerr << "getting a HashGraph graph" << endl;
        // unique_ptr<bdsg::HashGraph> graph;
        // get_input_file(optind, argc, argv, [&](istream &in) {
            // graph = vg::io::VPKG::load_one<bdsg::HashGraph>(in);
        // });

        // cerr << "got a HashGraph graph" << endl;
        // PathHandleGraph *hand_graph = dynamic_cast<PathHandleGraph*>(graph.get());
        // cerr << "pointer location " << hand_graph << endl;
        // // test(*hand_graph);
        //get a VG graph:
        cerr << "getting a vg graph" << endl;
        unique_ptr<VG> graph;
        get_input_file(optind, argc, argv, [&](istream &in) {
            graph = vg::io::VPKG::load_one<VG>(in);
        });
        cerr << "got a vg graph" << endl;
        // PathHandleGraph *hand_graph = dynamic_cast<PathHandleGraph*>(graph.get());

        //write graph
        // vg::io::VPKG::save(*dynamic_cast<bdsg::HashGraph *>(graph.get()), cout);
        vg::io::VPKG::save(*dynamic_cast<VG *>(graph.get()), cout);
        cerr << "saved the vg graph, now wrapped." << endl;
    // }
    return 0;
}

// Register subcommand
static Subcommand vg_vpkg_wrap("vpkg_wrap",
                               "Wraps given vg graph into vpkg format, saves as hg.", TOOLKIT,
                               main_vpkg_wrap);