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
        << "    -p / --paired" << endl
        << "    -R / --remap" << endl
        << "    -o / --output <PREFIX>" << endl
        << "Paired-end options:" << endl
        << "    -I / --insert-size <INSRTSZ>" << endl
        << "    -W / --insert-size-sigma <SIGMA>" << endl
        << "    -O / --one-end-anchored" << endl
        << "    -C / --interchromosomal" << endl
        << "    -D / --discordant-orientation" << endl
        << "Single-end / individual read options:" << endl
        << "    -c / --softclip <MAXCLIPLEN>" << endl
        << "    -s / --split-read <SPLITLEN>" << endl
        << "    -q / --quality <QUAL> " << endl
        << "    -d / --depth <DEPTH> " << endl
        << "    -i / --percent-identity <PCTID>" << endl
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
    int threads = 1;

    bool inverse = false;
    bool is_paired = false;
    bool remap = false;


    bool do_all = true;

    bool do_orientation = false;
    bool do_oea = false;
    bool do_insert_size = false;
    bool do_softclip = false;
    bool do_reversing = false;
    bool do_interchromosomal = false;
    bool do_split_read = false;
    bool do_unmapped = true;
    bool do_depth = false;
    bool do_percent_id = false;
    bool do_quality = false;


    float insert_size = 300;
    float insert_size_sigma = 50;


    int softclip_max = -1;
    int max_path_length = 0;

    int split_read_limit = -1;
    double depth = -1;


    int window = 0;
    bool use_avg = false;
    double pct_id = 0;
    double quality = 0.0;

    Filter ff;


    int c;
    optind = 2; // force optind past command positional argument
    while (true) {
        static struct option long_options[] =
        {
            {"help", no_argument, 0, 'h'},
            {0, 0, 0, 0}

        };
        int option_index = 0;
        c = getopt_long (argc, argv, "ht:vx:g:pRo:I:W:OCDc:s:q:d:i:aw:r",
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

            case 'u':
                do_unmapped = true;
                do_all = false;
                break;
            case 'c':
                do_softclip = true;
                softclip_max = atoi(optarg);
                ff.set_soft_clip_limit(softclip_max);
                break;
            case 's':
                do_split_read = true;
                break;
            case 'q':
                do_quality = true;
                quality = atof(optarg);
                break;
            case 'd':
                do_depth = true;
                depth = atof(optarg);
                break;
            case 'R':
                remap = true;
                break;
            case 'r':
                do_reversing = true;
                break;
            case 'p':
                is_paired = true;
                break;
            case 'C':
                do_interchromosomal = true;
                break;
            case 'I':
                insert_size = atof(optarg);
                break;
            case 'W':
                insert_size_sigma = atof(optarg);
                break;
            case 'O':
                do_oea = true;
                break;
            case 'D':
                do_orientation = true;
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
    ff.set_inverse(inverse);




    vector<Alignment> buffer;
    
    

    function<bool(Alignment&)> single_diff = [&](Alignment& aln){
        bool anchored = false;
        for (int i = 0; i < aln.path().mapping_size(); ++i){
            // check if anchored on both sides
            // check for small mismatches (to/from len diff or sequence)
        }
        return true;
    };

    vector<Alignment> orphaned_selected;
    vector<Alignment> discordant_selected;
    vector<Alignment> split_selected;
    vector<Alignment> reversing_selected;
    vector<Alignment> clipped_selected;
    vector<Alignment> one_end_anchored;
    vector<Alignment> quality_selected;
    vector<Alignment> insert_selected;
    vector<Alignment> depth_selected;
    vector<Alignment> single_bp_diff;
    vector<Alignment> unmapped_selected;

    string unmapped_fn = alignment_file + ".unmapped";
    string discordant_fn = alignment_file + ".discordant";
    string split_fn = alignment_file + ".split";
    string reversing_fn = alignment_file +  ".reversing";
    string oea_fn = alignment_file + ".oea";
    string clipped_fn = alignment_file + ".clipped";
    string insert_fn = alignment_file + ".insert_size";
    string quality_fn = alignment_file + ".quality";
    string depth_fn = alignment_file + ".depth";

    ofstream unmapped_stream;
    ofstream discordant_stream;
    ofstream split_stream;
    ofstream reversing_stream;
    ofstream oea_stream;
    ofstream clipped_stream;
    ofstream insert_stream;
    ofstream quality_stream;
    ofstream depth_stream;

    std::function<void()> write_all = [&](){

    };




    std::function<void(Alignment&, Alignment&)> pair_filters = [&](Alignment& alns_first, Alignment& alns_second){
        pair<Alignment, Alignment> ret;
        if (do_orientation){
            
            ret = ff.pair_orientation_filter(alns_first, alns_second);
            #pragma omp critical (discordant_selected)
            {
                discordant_selected.push_back(alns_first);
                discordant_selected.push_back(alns_second);
            }

        }
        if (do_oea){
            ret = ff.one_end_anchored_filter(alns_first, alns_second);
            #pragma omp critical (oea_selected)
            {
                one_end_anchored.push_back(alns_first);
                one_end_anchored.push_back(alns_second);
            }

        }
        if (do_insert_size){

            #pragma omp critial (insert_selected)
            {
                insert_selected.push_back(alns_first);
                insert_selected.push_back(alns_second);
            }


        }
        if (do_split_read){
            // first check softclip with a default size of 15
            // then get the position of the clip start
            // then remap the clipped/unclipped portions if required.
            // Then binary-search map if needed.

            if (remap){
                // ff.split_remap(alns_first, alns_second);

            }

            
            #pragma omp critical (split_selected)
            {
                split_selected.push_back(alns_first);
                split_selected.push_back(alns_second);
            }

        }
        if (do_reversing){
            Alignment x = ff.reversing_filter(alns_first);
            Alignment y = ff.reversing_filter(alns_second);
            if (x.name() != "" || y.name() != ""){
            #pragma omp critical (reversing_selected)
            {
                reversing_selected.push_back(alns_first);
                reversing_selected.push_back(alns_second);
            }
            }
       


        }
        if (do_softclip){

            Alignment x = ff.soft_clip_filter(alns_first);
            Alignment y = ff.soft_clip_filter(alns_second);

            #pragma omp critical (clipped_selected)
            {
                clipped_selected.push_back(alns_first);
                clipped_selected.push_back(alns_second);
            }

        }
        if (do_quality){

            #pragma omp critical (quality_selected)
            {
                quality_selected.push_back(alns_first);
                quality_selected.push_back(alns_second);
            }

        }
        if (do_depth){

            #pragma omp critical (depth_selected)
            {
                depth_selected.push_back(alns_first);
                depth_selected.push_back(alns_second);
            }
            
        }
        if (do_unmapped){
            if (alns_first.path().mapping_size() == 0 ||
                alns_second.path().mapping_size() == 0){

                    #pragma omp critical (unmapped_selected)
                    {
                        unmapped_selected.push_back(alns_first);
                        unmapped_selected.push_back(alns_second);
                    }
                }
        }


    };

    std::function<void(Alignment&)> single_filters = [&](Alignment& aln){
        if (do_split_read){
            split_selected.push_back(aln);
        }
        if (do_reversing){
            reversing_selected.push_back(aln);
        }
        if (do_softclip){
           if (!ff.soft_clip_filter(aln).name().empty()){
                clipped_selected.push_back(aln);
           }
           stream::write_buffered(clipped_stream, clipped_selected, 1000);

        }
        if (do_quality){
            quality_selected.push_back(aln);
        }
        if (do_depth){
            depth_selected.push_back(aln);
        }

    };



if (alignment_file == "-"){
    stream::for_each_interleaved_pair_parallel(cin, pair_filters);
}
else{
    ifstream in;
    in.open(alignment_file);
    if (in.good()){

        if (is_paired){
            cerr << "Processing..." << endl;
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


    buffer.clear();




return 0;
}

static Subcommand vg_sift("sift", "GAM filter", main_sift);


