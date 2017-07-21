/** \file sim_main.cpp
 *
 * Defines the "vg sim" subcommand, which generates potential reads from a graph.
 */


#include <omp.h>
#include <unistd.h>
#include <getopt.h>

#include <list>
#include <fstream>

#include "subcommand.hpp"

#include "../vg.hpp"
#include "../mapper.hpp"
#include "../sampler.hpp"

using namespace std;
using namespace vg;
using namespace vg::subcommand;

void help_sim(char** argv) {
    cerr << "usage: " << argv[0] << " sim [options]" << endl
         << "Samples sequences from the xg-indexed graph." << endl
         << endl
         << "options:" << endl
         << "    -x, --xg-name FILE    use the xg index in FILE" << endl
         << "    -l, --read-length N   write reads of length N" << endl
         << "    -n, --num-reads N     simulate N reads" << endl
         << "    -s, --random-seed N   use this specific seed for the PRNG" << endl
         << "    -e, --base-error N    base substitution error rate (default 0.0)" << endl
         << "    -i, --indel-error N   indel error rate (default 0.0)" << endl
         << "    -f, --forward-only    don't simulate from the reverse strand" << endl
         << "    -p, --frag-len N      make paired end reads with given fragment length N" << endl
         << "    -v, --frag-std-dev N  use this standard deviation for fragment length estimation" << endl
         << "    -N, --allow-Ns        allow reads to be sampled from the graph with Ns in them" << endl
         << "    -a, --align-out       generate true alignments on stdout rather than reads" << endl
         << "    -J, --json-out        write alignments in json" << endl;
}

int main_sim(int argc, char** argv) {

    if (argc == 2) {
        help_sim(argv);
        return 1;
    }

    int read_length = 100;
    int num_reads = 1;
    int seed_val = time(NULL);
    double base_error = 0;
    double indel_error = 0;
    bool forward_only = false;
    bool align_out = false;
    bool json_out = false;
    int fragment_length = 0;
    double fragment_std_dev = 0;
    bool reads_may_contain_Ns = false;
    string xg_name;

    int c;
    optind = 2; // force optind past command positional argument
    while (true) {
        static struct option long_options[] =
        {
            {"help", no_argument, 0, 'h'},
            {"xg-name", required_argument, 0, 'x'},
            {"read-length", required_argument, 0, 'l'},
            {"num-reads", required_argument, 0, 'n'},
            {"random-seed", required_argument, 0, 's'},
            {"forward-only", no_argument, 0, 'f'},
            {"align-out", no_argument, 0, 'a'},
            {"json-out", no_argument, 0, 'J'},
            {"allow-Ns", no_argument, 0, 'N'},
            {"base-error", required_argument, 0, 'e'},
            {"indel-error", required_argument, 0, 'i'},
            {"frag-len", required_argument, 0, 'p'},
            {"frag-std-dev", required_argument, 0, 'v'},
            {0, 0, 0, 0}
        };

        int option_index = 0;
        c = getopt_long (argc, argv, "hl:n:s:e:i:fax:Jp:v:N",
                long_options, &option_index);

        // Detect the end of the options.
        if (c == -1)
            break;

        switch (c)
        {

        case 'x':
            xg_name = optarg;
            break;

        case 'l':
            read_length = atoi(optarg);
            break;

        case 'n':
            num_reads = atoi(optarg);
            break;

        case 's':
            seed_val = atoi(optarg);
            break;

        case 'e':
            base_error = atof(optarg);
            break;

        case 'i':
            indel_error = atof(optarg);
            break;

        case 'f':
            forward_only = true;
            break;

        case 'a':
            align_out = true;
            break;

        case 'J':
            json_out = true;
            align_out = true;
            break;

        case 'N':
            reads_may_contain_Ns = true;
            break;

        case 'p':
            fragment_length = atoi(optarg);
            break;

        case 'v':
            fragment_std_dev = atof(optarg);
            break;

        case 'h':
        case '?':
            help_sim(argv);
            exit(1);
            break;

        default:
            abort ();
        }
    }

    if (xg_name.empty()) {
        cerr << "[vg sim] error: we need an xg index to sample reads from" << endl;
        return 1;
    }

    mt19937 rng;
    rng.seed(seed_val);

    xg::XG* xgidx = nullptr;
    ifstream xg_stream(xg_name);
    if(xg_stream) {
        xgidx = new xg::XG(xg_stream);
    }
    if (!xg_stream || xgidx == nullptr) {
        cerr << "[vg sim] error: could not open xg index" << endl;
        return 1;
    }

    Sampler sampler(xgidx, seed_val, forward_only, reads_may_contain_Ns);
    size_t max_iter = 1000;
    int nonce = 1;
    for (int i = 0; i < num_reads; ++i) {
        if (fragment_length) {
            auto alns = sampler.alignment_pair(read_length, fragment_length, fragment_std_dev, base_error, indel_error);
            size_t iter = 0;
            while (iter++ < max_iter) {
                if (alns.front().sequence().size() < read_length
                    || alns.back().sequence().size() < read_length) {
                    alns = sampler.alignment_pair(read_length, fragment_length, fragment_std_dev, base_error, indel_error);
                }
            }
            // write the alignment or its string
            if (align_out) {
                // write it out as requested
                if (json_out) {
                    cout << pb2json(alns.front()) << endl;
                    cout << pb2json(alns.back()) << endl;
                } else {
                    function<Alignment(uint64_t)> lambda = [&alns](uint64_t n) { return alns[n]; };
                    stream::write(cout, 2, lambda);
                }
            } else {
                cout << alns.front().sequence() << "\t" << alns.back().sequence() << endl;
            }
        } else {
            auto aln = sampler.alignment_with_error(read_length, base_error, indel_error);
            size_t iter = 0;
            while (iter++ < max_iter) {
                if (aln.sequence().size() < read_length) {
                    auto aln_prime = sampler.alignment_with_error(read_length, base_error, indel_error);
                    if (aln_prime.sequence().size() > aln.sequence().size()) {
                        aln = aln_prime;
                    }
                }
            }
            // write the alignment or its string
            if (align_out) {
                // write it out as requested
                if (json_out) {
                    cout << pb2json(aln) << endl;
                } else {
                    function<Alignment(uint64_t)> lambda = [&aln](uint64_t n) { return aln; };
                    stream::write(cout, 1, lambda);
                }
            } else {
                cout << aln.sequence() << endl;
            }
        }
    }

    return 0;
}

// Register subcommand
static Subcommand vg_sim("sim", "simulate reads from a graph", main_sim);

