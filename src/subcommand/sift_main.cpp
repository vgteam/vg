#include <iostream>
#include <fstream>
#include <ctime>
#include <cstdio>
#include <getopt.h>
#include <functional>
#include <omp.h>
#include "subcommand.hpp"
#include <gcsa/gcsa.h>
#include <vg/io/stream.hpp>
#include "vg/io/json2pb.h"
#include "../vg.hpp"
#include <vg/vg.pb.h>
#include "../filter.hpp"
#include "../alignment.hpp"

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
        //<< "    -v / --inverse      return the inverse of a query (like grep -v)"   << endl
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
        //<< "    -q / --quality <QUAL>               Flag reads with a single base quality below <QUAL>" << endl
        //<< "    -d / --depth <DEPTH>                Flag reads which have a low-depth Pos+Edit combo." << endl
        //<< "    -i / --percent-identity <PCTID>     Flag reads with percent identity to their primary path beliw <PCTID>" << endl
        //<< "    -a / --average" << endl
        //<< "    -w / --window <WINDOWLEN>" << endl
        << "    -r / --reversing                    Flag reads with sections that reverse within the read, as with small inversions ( --->|<-|-->  or ---->||<-- )" << endl
        << "Helpful helpers: " << endl
        << "    -1 / --calc-insert                  Calculate and print the insert size mean / sd every 1000 reads." << endl
        << endl;
}


int main_sift(int argc, char** argv){

    if (argc <= 2){
        help_sift(argv);
        return 1;
    }

    string alignment_file = "";
    string graph_name = "";
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

    bool just_calc_insert = false;


    float insert_size = 300;
    float insert_size_sigma = 50;


    int softclip_max = 15;
    int max_path_length = 200;

    int split_read_limit = 15;
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
        c = getopt_long (argc, argv, "hut:vgG:pRo:I:W:OCDc:s:q:d:i:aw:r1",
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
                threads = parse<int>(optarg);
                break;
            case 'G':
                graph_name = optarg;
                break;
            case 'u':
                do_unmapped = true;
                do_all = false;
                break;
            case 'c':
                do_softclip = true;
                do_all = false;
                softclip_max = parse<int>(optarg);
                ff.set_soft_clip_limit(softclip_max);
                break;
            case 's':
                do_split_read = true;
                split_read_limit = parse<int>(optarg);
                if (softclip_max < 0){
                    softclip_max = 15;
                    ff.set_soft_clip_limit(15);
                }
                do_all = false;
                break;
            case 'q':
                do_quality = true;
                do_all = false;
                quality = parse<double>(optarg);
                break;
            case 'd':
                do_depth = true;
                do_all = false;
                depth = parse<double>(optarg);
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
                insert_size = parse<double>(optarg);
                ff.insert_mean = insert_size;
                do_all = false;
                do_insert_size = true;
                break;
            case 'W':
                insert_size_sigma = parse<double>(optarg);
                ff.insert_sd = insert_size_sigma;
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
            case '1':
                just_calc_insert = true;
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

    vg::VG* graph;
    if (!graph_name.empty()){
        ifstream gstream(graph_name);
        ff.my_vg = new vg::VG(gstream);
    }

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
    
    

    //vector<Alignment> orphaned_selected;
    vector<Alignment> clean;
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
    string clean_fn = alignment_file + ".clean";
    string perfect_fn = alignment_file + ".perfect";

    ofstream unmapped_stream;
    ofstream discordant_stream;
    ofstream split_stream;
    ofstream reversing_stream;
    ofstream oea_stream;
    ofstream clipped_stream;
    ofstream insert_stream;
    ofstream quality_stream;
    ofstream depth_stream;
    ofstream clean_stream;
    ofstream perfect_stream;

    if (do_reversing){
        reversing_stream.open(reversing_fn);
    }

    if (do_unmapped){
        unmapped_stream.open(unmapped_fn);
    }

    if (do_orientation){
        discordant_stream.open(discordant_fn);
    }

    if (do_oea){
        oea_stream.open(oea_fn);
    }

    if (do_split_read){
        split_stream.open(split_fn);
    }
    if (do_insert_size){
        insert_stream.open(insert_fn);
    }
    if (do_softclip){
        clipped_stream.open(clipped_fn);
    }

    std::function<bool(Alignment, Alignment)> normalish = [&](Alignment a, Alignment b){
        bool a_forward = true;
        bool b_reverse = false;
        for (int i = 0; i < a.path().mapping_size(); i++){
            if (a.path().mapping(i).position().is_reverse()){
                a_forward = false;
            }
        }
        for (int i = 0; i < b.path().mapping_size(); i++){
            if (b.path().mapping(i).position().is_reverse()){
                b_reverse = true;
            }
        }
        return (a_forward != b_reverse);
    };


    std::function<void(Alignment&, Alignment&)> pair_filters = [&](Alignment& alns_first, Alignment& alns_second){
        bool ret;
        bool flagged = false;

        if (do_unmapped && !flagged){
            if (ff.unmapped_filter(alns_first) && ff.unmapped_filter(alns_second)){

                    #pragma omp critical (unmapped_selected)
                    {
                        flagged = true;
                        alns_first.set_read_mapped(false);
                        alns_first.set_mate_unmapped(true);
                        alns_second.set_read_mapped(false);
                        alns_second.set_mate_unmapped(false);
                        unmapped_selected.push_back(alns_first);
                        unmapped_selected.push_back(alns_second);
                        
                    }
                }
        }

        if (do_orientation && !flagged){
            
            ret = ff.pair_orientation_filter(alns_first, alns_second);
            if (ret){
                #pragma omp critical (discordant_selected)
            {
                flagged = true;
                discordant_selected.push_back(alns_first);
                discordant_selected.push_back(alns_second);
            }
            }
            

        }
        if (do_oea && !flagged){
            ret = ff.one_end_anchored_filter(alns_first, alns_second);
            if (ret){
                #pragma omp critical (oea_selected)
                {
                    one_end_anchored.push_back(alns_first);
                    one_end_anchored.push_back(alns_second);
                }
            }
            

        }
        if (do_insert_size && !flagged){
            
            ret = ff.insert_size_filter(alns_first, alns_second);
            if (ret){
                #pragma omp critical (insert_selected)
                {
                    insert_selected.push_back(alns_first);
                    insert_selected.push_back(alns_second);
                }
            }
            


        }
        if (do_split_read && !flagged){

            if (ff.split_read_filter(alns_first)){
                #pragma omp critical (split_selected)
                {
                    split_selected.push_back(alns_first);
                    split_selected.push_back(alns_second);
                    flagged = true;

                }
                
            }
            if (ff.split_read_filter(alns_second)){
                #pragma omp critical (split_selected)
                {
                    split_selected.push_back(alns_first);
                    split_selected.push_back(alns_second);
                    flagged = true;

                }
            }
            

        }
        if (do_reversing && !flagged){
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
        if (do_softclip && !flagged){

            bool x = ff.soft_clip_filter(alns_first);
            bool y = ff.soft_clip_filter(alns_second);
            if (x){
                #pragma omp critical (clipped_selected)
                {
                    clipped_selected.push_back(alns_first);
                    flagged = true;
                }
            } 
            if (y){
                #pragma omp critical (clipped_selected)
                {
                    flagged = true;
                    clipped_selected.push_back(alns_second);
                }
            } 

        }
        if (do_quality && !flagged){

            #pragma omp critical (quality_selected)
            {
                quality_selected.push_back(alns_first);
                quality_selected.push_back(alns_second);
            }

        }
        if (do_depth && !flagged){

            #pragma omp critical (depth_selected)
            {
                depth_selected.push_back(alns_first);
                depth_selected.push_back(alns_second);
            }
            
        }
        if (!flagged){
            // Check if read is perfect

            if (ff.perfect_filter(alns_first) && ff.perfect_filter(alns_second)){

            }
            else{
                // otherwise, it's pretty clean, so place it in the clean pile.
                #pragma omp critical (clean)
                {
                    clean.push_back(alns_first);
                    clean.push_back(alns_second);
                }
            }
            
        }
        

        vg::io::write_buffered(unmapped_stream, unmapped_selected, 100);
        vg::io::write_buffered(discordant_stream, discordant_selected, 100);
        vg::io::write_buffered(oea_stream, one_end_anchored, 100);
        vg::io::write_buffered(insert_stream, insert_selected, 100);
        vg::io::write_buffered(split_stream, split_selected, 100);
        vg::io::write_buffered(clipped_stream, clipped_selected, 100);
        vg::io::write_buffered(clean_stream, clean, 100);
        vg::io::write_buffered(reversing_stream, reversing_selected, 100);
        vg::io::write_buffered(perfect_stream, perfect, 100);
    };

    std::function<void(Alignment&)> single_filters = [&](Alignment& aln){
        if (do_split_read){
            split_selected.push_back(aln);
        }
        if (do_reversing){
            reversing_selected.push_back(aln);
        }
        if (do_softclip){
           if (ff.soft_clip_filter(aln)){
                clipped_selected.push_back(aln);
           }
           //vg::io::write_buffered(clipped_stream, clipped_selected, 1000);

        }
        if (do_quality){
            quality_selected.push_back(aln);
        }
        if (do_depth){
            depth_selected.push_back(aln);
        }

    };

    // double curr_insert = 0.0;
    // double insert_tot = 0.0;
    // int insert_count = 0;
    // double curr_sigma = 0.0;
    // std::function<void(Alignment&, Alignment&)> calc_insert = [&](Alignment& a, Alignment& b){
    //     if (a.fragment_size() > 1){
    //         insert_count += 1;
    //         insert_tot += (double) a.fragment(0).length();   
    //     }
    //     else if (b.fragment_size() > 1){

    //     }
    // };





if (alignment_file == "-"){
    vg::io::for_each_interleaved_pair_parallel(cin, pair_filters);
}
else{
    ifstream in;
    in.open(alignment_file);
    if (in.good()){

        // if (just_calc_insert){
        //     vg::io::for_each_interleaved_pair_parallel(in, calc_insert);
        //     exit(0);
        // }

        if (is_paired){
            cerr << "Processing..." << endl;
            vg::io::for_each_interleaved_pair_parallel(in, pair_filters);
        }
        else{
            vg::io::for_each_parallel(in, single_filters);

        }
    }
    else{
        cerr << "Could not open " << alignment_file << endl;
        help_sift(argv);
    }
}
    vg::io::write_buffered(unmapped_stream, unmapped_selected, 0);
    vg::io::write_buffered(discordant_stream, discordant_selected, 0);
    vg::io::write_buffered(oea_stream, one_end_anchored, 0);
    vg::io::write_buffered(insert_stream, insert_selected, 0);
    vg::io::write_buffered(split_stream, split_selected, 0);
    vg::io::write_buffered(clipped_stream, clipped_selected, 0);
    vg::io::write_buffered(clean_stream, clean, 0);
    vg::io::write_buffered(reversing_stream, reversing_selected, 0);
    vg::io::write_buffered(perfect_stream, perfect, 0);

    buffer.clear();




    return 0;
}

static Subcommand vg_sift("sift", "Filter Alignments by various metrics related to variant calling.", DEPRECATED, main_sift);


