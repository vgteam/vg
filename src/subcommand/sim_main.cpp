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
         << "    -x, --xg-name FILE          use the xg index in FILE" << endl
         << "    -F, --fastq FILE            superpose errors matching the error profile of NGS reads in FILE (ignores -l,-f)" << endl
         << "    -P, --path PATH             simulate from the given names path (multiple allowed)" << endl
         << "    -l, --read-length N         write reads of length N" << endl
         << "    -n, --num-reads N           simulate N reads or read pairs" << endl
         << "    -s, --random-seed N         use this specific seed for the PRNG" << endl
         << "    -e, --sub-rate FLOAT        base substitution rate (default 0.0)" << endl
         << "    -i, --indel-rate FLOAT      indel rate (default 0.0)" << endl
         << "    -d, --indel-err-prop FLOAT  proportion of trained errors from -F that are indels (default 0.0)" << endl
         << "    -S, --scale-err FLOAT       scale trained error probabilities from -F by this much (default 1.0)" << endl
         << "    -f, --forward-only          don't simulate from the reverse strand" << endl
         << "    -p, --frag-len N            make paired end reads with given fragment length N" << endl
         << "    -v, --frag-std-dev FLOAT    use this standard deviation for fragment length estimation" << endl
         << "    -N, --allow-Ns              allow reads to be sampled from the graph with Ns in them" << endl
         << "    -a, --align-out             generate true alignments on stdout rather than reads" << endl
         << "    -J, --json-out              write alignments in json" << endl;
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
    bool strip_bonuses = false;
    double indel_prop = 0.0;
    double error_scale_factor = 1.0;
    string fastq_name;
    // What path should we sample from? Empty string = the whole graph.
    vector<string> path_names;

    int c;
    optind = 2; // force optind past command positional argument
    while (true) {
        static struct option long_options[] =
        {
            {"help", no_argument, 0, 'h'},
            {"xg-name", required_argument, 0, 'x'},
            {"fastq", required_argument, 0, 'F'},
            {"path", required_argument, 0, 'P'},
            {"read-length", required_argument, 0, 'l'},
            {"num-reads", required_argument, 0, 'n'},
            {"random-seed", required_argument, 0, 's'},
            {"forward-only", no_argument, 0, 'f'},
            {"align-out", no_argument, 0, 'a'},
            {"json-out", no_argument, 0, 'J'},
            {"allow-Ns", no_argument, 0, 'N'},
            {"base-rate", required_argument, 0, 'e'},
            {"indel-rate", required_argument, 0, 'i'},
            {"indel-err-prop", required_argument, 0, 'd'},
            {"scale-err", required_argument, 0, 'S'},
            {"frag-len", required_argument, 0, 'p'},
            {"frag-std-dev", required_argument, 0, 'v'},
            {0, 0, 0, 0}
        };

        int option_index = 0;
        c = getopt_long (argc, argv, "hl:n:s:e:i:fax:Jp:v:Nd:F:P:S:",
                long_options, &option_index);

        // Detect the end of the options.
        if (c == -1)
            break;

        switch (c)
        {

        case 'x':
            xg_name = optarg;
            break;
            
        case 'F':
            fastq_name = optarg;
            break;
            
        case 'P':
            path_names.push_back(optarg);
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
            
        case 'd':
            indel_prop = atof(optarg);
            break;
            
        case 'S':
            error_scale_factor = atof(optarg);
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

    for (auto& path_name : path_names) {
        if (xgidx->path_rank(path_name) == 0) {
            cerr << "[vg sim] error: path \""<< path_name << "\" not found in index" << endl;
            return 1;
        }
    }

    // Make a sample to sample reads with
    Sampler sampler(xgidx, seed_val, forward_only, reads_may_contain_Ns);
    
    // Make a Mapper to score reads, with the default parameters
    Mapper rescorer(xgidx, nullptr, nullptr);
    // We define a function to score a generated alignment under the mapper
    auto rescore = [&] (Alignment& aln) {
        // Score using exact distance.
        aln.set_score(rescorer.score_alignment(aln, false));
    };
    
    if (fastq_name.empty()) {
        // Use the fixed error rate sampler
        
        // Make a sample to sample reads with
        Sampler sampler(xgidx, seed_val, forward_only, reads_may_contain_Ns, path_names);
        
        // Make a Mapper to score reads, with the default parameters
        Mapper rescorer(xgidx, nullptr, nullptr);
        // Override the "default" full length bonus, just like every other subcommand that uses a mapper ends up doing.
        // TODO: is it safe to change the default?
        rescorer.set_alignment_scores(default_match, default_mismatch, default_gap_open, default_gap_extension, 5);
        // Include the full length bonuses if requested.
        rescorer.strip_bonuses = strip_bonuses;
        // We define a function to score a generated alignment under the mapper
        auto rescore = [&] (Alignment& aln) {
            // Score using exact distance.
            aln.set_score(rescorer.score_alignment(aln, false));
        };
        
        size_t max_iter = 1000;
        int nonce = 1;
        for (int i = 0; i < num_reads; ++i) {
            // For each read we are going to generate
            
            if (fragment_length) {
                // fragment_lenght is nonzero so make it two paired reads
                auto alns = sampler.alignment_pair(read_length, fragment_length, fragment_std_dev, base_error, indel_error);
                
                size_t iter = 0;
                while (iter++ < max_iter) {
                    // For up to max_iter iterations
                    if (alns.front().sequence().size() < read_length
                        || alns.back().sequence().size() < read_length) {
                        // If our read was too short, try again
                        alns = sampler.alignment_pair(read_length, fragment_length, fragment_std_dev, base_error, indel_error);
                    }
                }
                
                // write the alignment or its string
                if (align_out) {
                    // write it out as requested
                    
                    // We will need scores
                    rescore(alns.front());
                    rescore(alns.back());
                    
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
                // Do single-end reads
                auto aln = sampler.alignment_with_error(read_length, base_error, indel_error);
                
                size_t iter = 0;
                while (iter++ < max_iter) {
                    // For up to max_iter iterations
                    if (aln.sequence().size() < read_length) {
                        // If our read is too short, try again
                        auto aln_prime = sampler.alignment_with_error(read_length, base_error, indel_error);
                        if (aln_prime.sequence().size() > aln.sequence().size()) {
                            // But only keep the new try if it is longer
                            aln = aln_prime;
                        }
                    }
                }
                
                // write the alignment or its string
                if (align_out) {
                    // write it out as requested
                    
                    // We will need scores
                    rescore(aln);
                    
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
        
    }
    else {
        // Use the trained error rate
        
        Aligner aligner(default_match, default_mismatch, default_gap_open, default_gap_extension, 5);
        
        NGSSimulator sampler(*xgidx, fastq_name, path_names, base_error, indel_error, indel_prop,
                             fragment_length ? fragment_length : std::numeric_limits<double>::max(), // suppresses warnings about fragment length
                             fragment_std_dev ? fragment_std_dev : 0.000001, // eliminates errors from having 0 as stddev without substantial difference
                             error_scale_factor, !reads_may_contain_Ns, seed_val);
        
        if (fragment_length) {
            for (size_t i = 0; i < num_reads; i++) {
                pair<Alignment, Alignment> read_pair = sampler.sample_read_pair();
                read_pair.first.set_score(aligner.score_ungapped_alignment(read_pair.first, strip_bonuses));
                read_pair.second.set_score(aligner.score_ungapped_alignment(read_pair.second, strip_bonuses));
                
                if (align_out) {
                    if (json_out) {
                        cout << pb2json(read_pair.first) << endl;
                        cout << pb2json(read_pair.second) << endl;
                    }
                    else {
                        function<Alignment(uint64_t)> lambda = [&read_pair](uint64_t n) {
                            return n % 2 ? read_pair.first : read_pair.second;
                        };
                        stream::write(cout, 2, lambda);
                    }
                }
                else {
                    cout << read_pair.first.sequence() << "\t" << read_pair.second.sequence() << endl;
                }
            }
        }
        else {
            for (size_t i = 0; i < num_reads; i++) {
                Alignment read = sampler.sample_read();
                read.set_score(aligner.score_ungapped_alignment(read, strip_bonuses));
                
                if (align_out) {
                    if (json_out) {
                        cout << pb2json(read) << endl;
                    }
                    else {
                        function<Alignment(uint64_t)> lambda = [&read](uint64_t n) {
                            return read;
                        };
                        stream::write(cout, 1, lambda);
                    }
                }
                else {
                    cout << read.sequence() << endl;
                }
            }
        }
    }
    

    return 0;
}

// Register subcommand
static Subcommand vg_sim("sim", "simulate reads from a graph", main_sim);

