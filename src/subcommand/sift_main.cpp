#include <iostream>
#include <fstream>
#include <ctime>
#include <cstdio>
#include <getopt.h>
#include <functional>
#include <omp.h>
#include "subcommand.hpp"
#include "gcsa.h"
#include "stream.hpp"
// From gcsa2
#include "files.h"
#include "json2pb.h"
#include "vg.hpp"
#include "vg.pb.h"
#include "filter.hpp"
#include "alignment.hpp"

using namespace std;
using namespace vg;
using namespace vg::subcommand;

//TODO ideal behavior is to filter READ PAIRS
//when a mate fails one of the individual read filters.
//
//Max path len should be TODO s > avg_p_len + 3*sigma
//and there should be min path len, essentially discordant path lengths

/*
 *
 
bool is_significant(double zed, double cutoff){
    if (zed > cutoff || -1.0 * zed > cutoff){
        return true;
    }
    return false;
}
double zed_score(double mu, double sigma, double val){
    return (val - mu) / sigma;
}

double reservoir_avg(int sz, vector<double> vals){

}

*/

void help_sift(char** argv){
    cerr << "Usage: " << argv[0] << " sift [options] <alignments.gam>" << endl
        << "Sift through a GAM and select / remove reads with particular properties." << endl
        << "General Options: " << endl
        << "    -t / --threads" << endl
        << "    -v / --inverse" << endl
        << "    -x / --xg" << endl
        << "    -P / --paired" << endl
        << "Paired-end options:" << endl
        << "    -I / --insert-size <INSRTSZ>" << endl
        << "    -w / --insert-size-sigma <SIGMA>" << endl
        << "    -O / --one-end-anchored" << endl
        << "    -C / --interchromosomal" << endl
        << "    -D / --discordant-orientation" << endl
        << "    -F / --force-paired" << endl
        << "    -M / --max-path-length <MXPTHLN>" << endl
        << "Single-end / individual read options:" << endl
        << "    -c / --softclip <CLIPLEN>" << endl
        << "    -s / --split-read <SPLITLEN>" << endl
        << "    -q / --quality <QUAL> " << endl
        << "    -d / --depth <DEPTH> " << endl
        << "    -p / --percent-identity <PCTID>" << endl
        << "    -a / --average" << endl
        << "    -w / --window <WINDOWLEN>" << endl
        << "    -r / --reversing" << endl
        << endl;
}


int main_sift(int argc, char** argv){

    if (argc <= 2){
        help_sift(argv);
        return 1;
    }

    string alignment_file = "";

    bool inverse = false;
    bool is_paired = false;
    bool forced = false;

    bool do_orientation = false;
    bool do_oea = false;
    int threads = 1;

    int insert_size = 300;
    int insert_size_sigma = 50;
    bool do_insrt_sz = false;

    bool do_interchromosomal = false;

    int soft_clip_max = -1;
    int max_path_length = 0;
    bool do_max_path = false;

    int split_read_limit = -1;
    double depth = -1;
    int window = 0;
    bool use_avg = false;
    double pct_id = 0;
    double quality = 0.0;




    int c;
    optind = 2; // force optind past command positional argument
    while (true) {
        static struct option long_options[] =
        {
            {"help", no_argument, 0, 'h'},
            {"one-end-anchored", no_argument, 0, 'O'},
            {"interchromosomal", no_argument, 0, 'C'},
            {"discordant-orientation", no_argument, 0, 'D'},
            {"paired", no_argument, 0, 'P'},
            {"insert-size", required_argument, 0, 'I'},
            {"insert-size-sigma", required_argument, 0, 'w'},
            {"soft-clip", required_argument, 0, 'c'},
            {"max-path-length", required_argument, 0, 'M'},
            {0, 0, 0, 0}

        };
        int option_index = 0;
        c = getopt_long (argc, argv, "htvM:I:w:ODPFc:",
                long_options, &option_index);

        // Detect the end of the options.
        if (c == -1)
            break;

        switch (c)
        {
            case '?':
            case 'h':
                help_sift(argv);
                return 1;
            case 't':
                threads = atoi(optarg);
                break;
            case 'v':
                inverse = true;
                break;
            case 's':
                soft_clip_max = atoi(optarg);
                break;
            case 'F':
                is_paired = true;
                forced = true;
                break;
            case 'O':
                do_oea = true;
                break;
            case 'I':
                insert_size = atoi(optarg);
                do_insrt_sz = true;
                break;
            case 'w':
                insert_size_sigma = atoi(optarg);
                do_insrt_sz = true;
                break;
            case 'P':
                is_paired = true;
                break;
            case 'D':
                do_orientation = true;
                break;
            case 'M':
                do_max_path = true;
                max_path_length = atoi(optarg);
                break;
            case 'c':
                soft_clip_max = atoi(optarg);
                break;
            default:
                help_sift(argv);
                exit(1);
                break;

        }

    }

    if (optind >= argc){
        cerr << "Error: no alignment file given" << endl;
        exit(1);
    }
    alignment_file = argv[optind];

    cerr << "Filtering  " << alignment_file << endl;


    omp_set_num_threads(threads);

    Filter ff;
    ff.set_inverse(inverse);
    ff.max_path_length = max_path_length;

    if (soft_clip_max >= 0){
        ff.set_soft_clip_limit(soft_clip_max);
    }

    vector<Alignment> buffer;
    static const int buffer_size = 1000; // we let this be off by 1
    function<Alignment&(uint64_t)> write_buffer = [&buffer](uint64_t i) -> Alignment& {
        return buffer[i];
    };



    std::function<void(Alignment&, Alignment&)> pair_filters = [&](Alignment& alns_first, Alignment& alns_second){

        /**  if (aln.sequence().size() > 0){
          buffer.push_back(aln);
          }
          };*/

        if (do_orientation){
            std::pair<Alignment, Alignment> filt = ff.orientation_filter(alns_first, alns_second);
            alns_first = filt.first;
            alns_second = filt.second;
            if (!alns_first.name().empty()){
                buffer.push_back(alns_first);
                buffer.push_back(alns_second);
            }
        }
        if (do_oea){
            std::pair<Alignment, Alignment> filt = ff.one_end_anchored_filter(alns_first, alns_second);
            alns_first = filt.first;
            alns_second = filt.second;
            if (!alns_first.name().empty()){
                buffer.push_back(alns_first);
                buffer.push_back(alns_second);
            }
        }
        if (do_insrt_sz){

        }
        if (do_interchromosomal){

        }
        if (do_max_path ){
            std::pair<Alignment, Alignment> filt = ff.path_length_filter(alns_first, alns_second);
            alns_first = filt.first;
            alns_second = filt.second;
            if (!alns_first.name().empty()){
                buffer.push_back(alns_first);
                buffer.push_back(alns_second);
            }
        }
    //return std::make_pair(alns.first, alns.second);
};

std::function<void(Alignment&, Alignment&)> se_pair_filters = [&](Alignment& alns_first, Alignment& alns_second){
    

};

std::function<void(Alignment&)> single_filters = [&](Alignment& aln){
    if (soft_clip_max >=0){
        aln = ff.soft_clip_filter(aln);    
    }
    buffer.push_back(aln);
};



if (alignment_file == "-"){
}
else{
    ifstream in;
    in.open(alignment_file);
    if (in.good()){

        if (is_paired && forced){
            cerr << "Processing..." << endl;
            stream::for_each_interleaved_pair_parallel(in, pair_filters);
        }
        else if (is_paired){
            stream::for_each_interleaved_pair_parallel(in, pair_filters);
        }
        else{
            stream::for_each_parallel(in, single_filters);

        }
    }
    else{
        cerr << "Could not open " << alignment_file << endl;
        help_sift(argv);
    }
}

if (buffer.size() > 0) {
    stream::write(cout, buffer.size(), write_buffer);
    buffer.clear();
}



return 0;
}

static Subcommand vg_sift("sift", "GAM filter / scrubber", main_sift);


