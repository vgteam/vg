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

void help_sift(char** argv){
    cerr << "Usage: " << argv[0] << " sift [options] <alignments.gam>" << endl
        << "Sift through a GAM and select / remove reads with particular properties." << endl
        << "General Options: " << endl
        << "    -t / --threads  <MTHRDS>    number of OMP threads (not all algorithms are parallelized)." << endl
        << "    -v / --inverse      return the inverse of a query (like grep -v)"   << endl
        << "    -x / --xg  <MYXG>   An XG index (for realignment of split reads)" << endl
        << "    -p / --paired       Input reads are paired-end" << endl
        << "    -R / --remap        Remap (locally) any soft-clipped, split, or discordant read pairs." << endl
        << "    -o / --output <PREFIX>" << endl
        << "Paired-end options:" << endl
        << "    -I / --insert-size <INSRTSZ>        Insert size mean. Flag reads where ((I - insrtsz) / W) > 2.95" << endl
        << "    -W / --insert-size-sigma <SIGMA>    Standard deviation of insert size." << endl
        << "    -O / --one-end-anchored             Flag reads where one read of the pair is mapped and the other is unmapped." << endl
        << "    -C / --interchromosomal             Flag reads mapping to two distinct Paths" << endl
        << "    -D / --discordant-orientation       Flag reads that do not have the expected --> <-- orientation." << endl
        << "Single-end / individual read options:" << endl
        << "    -c / --softclip <MAXCLIPLEN>        Flag reads with softclipped sections longer than MAXCLIPLEN" << endl
        << "    -s / --split-read <SPLITLEN>        Flag reads with softclips that map independently of the anchored portion (requires -x, -g)." << endl
        << "    -q / --quality <QUAL>               Flag reads with a single base quality below <QUAL>" << endl
        << "    -d / --depth <DEPTH>                Flag reads which have a low-depth Pos+Edit combo." << endl
        << "    -i / --percent-identity <PCTID>     Flag reads with percent identity to their primary path beliw <PCTID>" << endl
        //<< "    -a / --average" << endl
        //<< "    -w / --window <WINDOWLEN>" << endl
        << "    -r / --reversing                    Flag reads with sections that reverse within the read, as with small inversions ( --->|<-|--> )" << endl
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
    bool do_depth = false;
    bool do_percent_id = false;
    bool do_quality = false;
    bool do_unmapped = false;


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
                do_all = false;
                break;
            case 'q':
                do_quality = true;
                do_all = false;
                quality = atof(optarg);
                break;
            case 'd':
                do_depth = true;
                do_all = false;
                depth = atof(optarg);
                break;
            case 'R':
                remap = true;
                do_all = false;
                break;
            case 'r':
                do_reversing = true;
                do_all = false;
                break;
            case 'p':
                is_paired = true;
                break;
            case 'C':
                do_interchromosomal = true;
                do_all = false;
                break;
            case 'I':
                insert_size = atof(optarg);
                do_all = false;
                break;
            case 'W':
                insert_size_sigma = atof(optarg);
                do_all = false;
                break;
            case 'O':
                do_oea = true;
                do_all = false;
                break;
            case 'D':
                do_orientation = true;
                do_all = false;
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

    if (softclip_max >= 0){
        ff.set_soft_clip_limit(softclip_max);
    }


    if (do_all){
    do_orientation = true;
    do_oea = true;
    do_insert_size = true;
    do_softclip = true;
    do_reversing = true;
    do_interchromosomal = true;
    do_split_read = true;
    do_depth = true;
    do_percent_id = true;
    do_quality = true;
    do_unmapped = true;

    }

    vector<Alignment> buffer;
    
    

    function<bool(Alignment&)> single_diff = [&](Alignment& aln){
        bool anchored = false;
        for (int i = 0; i < aln.path().mapping_size(); ++i){
            // check if anchored on both sides
            // check for small mismatches (to/from len diff or sequence)
        }
        return true;
    };

    //vector<Alignment> orphaned_selected;
    vector<Alignment> just_fine;
    vector<Alignment> perfect;
    vector<Alignment> simple_mismatch;
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
    string oea_fn = alignment_file + ".one_end_anchored";
    string clipped_fn = alignment_file + ".softclipped";
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

    unmapped_stream.open(unmapped_fn);

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

        stream::write_buffered(unmapped_stream, unmapped_selected, 100);


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
    stream::write_buffered(unmapped_stream, unmapped_selected, 0);
    stream::write_buffered(oea_stream, one_end_anchored, 0);
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

static Subcommand vg_sift("sift", "Filter Alignments by various metrics related to variant calling.", main_sift);


