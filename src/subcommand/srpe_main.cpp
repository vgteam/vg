#include <iostream>
#include <vector>
#include <getopt.h>
#include <functional>
#include <regex>
#include "subcommand.hpp"
#include "srpe.hpp"
#include "stream.hpp"
#include "index.hpp"
#include "position.hpp"
#include "vg.pb.h"
#include "path.hpp"
#include "genotypekit.hpp"
#include "genotyper.hpp"
#include "path_index.hpp"
#include "vg.hpp"
#include <math.h>
#include "srpe.hpp"
#include "filter.hpp"
#include "utility.hpp"
#include "Variant.h"
#include "translator.hpp"
#include "Fasta.h"
#include "IntervalTree.h"

using namespace std;
using namespace vg;
using namespace vg::subcommand;


void help_srpe(char** argv){
    cerr << "Usage: " << argv[0] << " srpe [options] <data.gam> <graph.vg>" << endl
    << "Options: " << endl 
    << "  -p / --ref-path" << endl
    << "  -x / --xg" << endl 
    << "  -g / --gcsa" << endl 
    << endl;
}



int main_srpe(int argc, char** argv){
    string alignment_file = "";
    string gam_index_name = "";
    string graph_name = "";
    string xg_name = "";
    string gcsa_name = "";
    string lcp_name = "";

    string spec_vcf = "";
    string ref_fasta = "";
    string ins_fasta = "";

    string augmented_graph_name = "";
    bool augment_paths = true;

    string ref_path = "";

    int max_iter = 2;
    int max_frag_len = 10000;
    int min_soft_clip = 12;

    bool do_all = false;

    vector<string> search_types;
    search_types.push_back("DEL");

    int threads = 1;

    if (argc <= 2) {
        help_srpe(argv);
        return 1;
    }

    int c;
    optind = 2; // force optind past command positional argument
    while (true) {
        static struct option long_options[] =
        {
            {"max-iter", required_argument, 0, 'm'},
            {"xg-index", required_argument, 0, 'x'},
            {"augmented", required_argument, 0, 'a'},
            {"help", no_argument, 0, 'h'},
            {"gcsa-index", required_argument, 0, 'g'},
            {"specific", required_argument, 0, 'S'},
            {"recall", no_argument, 0, 'R'},
            {"insertions", required_argument, 0, 'I'},
            {"reference", required_argument, 0, 'r'},
            {"threads", required_argument, 0, 't'},
            {"ref-path", required_argument, 0, 'p'},
            {0, 0, 0, 0}

        };
        int option_index = 0;
        c = getopt_long (argc, argv, "hx:g:m:S:RI:r:t:a:wp:",
                long_options, &option_index);

        // Detect the end of the options.
        if (c == -1)
            break;

        switch (c)
        {
            case 'a':
                augmented_graph_name = optarg;
                break;
            case 'm':
                max_iter = atoi(optarg);
                break;

            case 't':
                threads = atoi(optarg);
                break;

            case 'R':
                do_all = true;
                break;
            case 'x':
                xg_name = optarg;
                break;
            case 'g':
                gcsa_name = optarg;
                break;
            case 'S':
                spec_vcf = optarg;
                break;
            case 'r':
                ref_fasta = optarg;
                break;
            case 'I':
                ins_fasta = optarg;
                break;
            case 'p':
                ref_path = optarg;
                break;
            case 'h':
            case '?':
            default:
                help_srpe(argv);
                abort();
        }

    }

    omp_set_num_threads(threads);


    SRPE srpe;

    


    alignment_file = argv[optind];
    //gam_index_name = argv[++optind];
    graph_name = argv[++optind];

    xg::XG* xg_ind = new xg::XG();
    Index gamind;

    vg::VG* graph;

    if (!xg_name.empty()){
        ifstream in(xg_name);
        xg_ind->load(in);
        srpe.ff.set_my_xg_idx(xg_ind);
    }
    // Set GCSA indexes
    if (!gcsa_name.empty()){
            ifstream gcsa_stream(gcsa_name);
            srpe.ff.gcsa_ind = new gcsa::GCSA();
            srpe.ff.gcsa_ind->load(gcsa_stream);
            string lcp_name = gcsa_name + ".lcp";
            ifstream lcp_stream(lcp_name);
            srpe.ff.lcp_ind = new gcsa::LCPArray();
            srpe.ff.lcp_ind->load(lcp_stream);
    }
    if (!xg_name.empty()){
        ifstream xgstream(xg_name);
        xg_ind->load(xgstream);
        srpe.ff.set_my_xg_idx(xg_ind);
    }
    srpe.ff.init_mapper();
    // else{

    // }

    if (!graph_name.empty()){
        ifstream in(graph_name);
        graph = new VG(in, false);
    }

    // Open a variant call file,
    // hash each variant to an hash ID
    // have in if in the loop below.
    map<string, Locus> name_to_loc;
    // Makes a pathindex, which allows us to query length and push out a VCF with a position
    map<string, PathIndex*> pindexes;
    Genotyper gg;
    regex is_alt ("_alt_.*");

    vector<FastaReference*> insertions;
    if (!ins_fasta.empty()){
        FastaReference* ins = new FastaReference();
        insertions.emplace_back(ins);
        ins->open(ins_fasta);

    }

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

    ifstream unmapped_stream;
    ifstream discordant_stream;
    ifstream split_stream;
    ifstream reversing_stream;
    ifstream oea_stream;
    ifstream clipped_stream;
    ifstream insert_stream;
    ifstream quality_stream;
    ifstream depth_stream;
    vector<Interval<int> > inters;

    double insert_size = 0.0;
    double insert_var = 0.0;


    
    
    std::function<vector<BREAKPOINT> (vector<BREAKPOINT>)> merge_breakpoints = [&](vector<BREAKPOINT> bps){
        vector<BREAKPOINT> ret;
        BREAKPOINT sent;
        sent.start = -100;
        ret.push_back(sent);
        for (int i = 0; i < bps.size(); i++){
            BREAKPOINT a = bps[i];
            bool merged = false;
            for (int j = 0; j < ret.size(); j++){
                if (ret[j].overlap(a, 20)){
                    ret[j].other_supports += 1;
                    merged = true;
                }
            }
            if (!merged){
                ret.push_back(a);
            }
        }
        return ret;
    };

    std::function<vector<INS_INTERVAL>(vector<INS_INTERVAL>)> merge_intervals = [&](vector<INS_INTERVAL> ins){
        vector<INS_INTERVAL> ret;
        // sentinel interval
        INS_INTERVAL sent;
        sent.start = -1;
        sent.end = -1;
        for (int i = 0; i < ins.size(); i++){
            bool merged = false;
            for (int j = 0; j < ret.size(); j++){
                if (ret[j].overlap(ins[i])){
                    ret[j].other_supports+=1;
                    merged = true;
                }
            }
            if (!merged){
                ret.push_back(ins[i]);
            }

        }
        return ret;
    };

    // Set up path index
    srpe.ff.set_my_vg(graph);
    srpe.ff.soft_clip_limit = 20;
    srpe.ff.fill_node_to_position(ref_path);
    

    vector<INS_INTERVAL> ins;
    vector<BREAKPOINT> bps;

    std::function<void(vector<INS_INTERVAL>&, map<int64_t, vector<INS_INTERVAL> >&)> merge = [&](vector<INS_INTERVAL>& ins, map<int64_t, vector<INS_INTERVAL> >& start_to_interval){

    };
    
    std::function<void(Alignment&, Alignment&)> calc_insert_size = [&](Alignment& a, Alignment& b){
        insert_size = 1000.0;
        insert_var = 300.0;
    };

    std::function<void(Alignment&, Alignment&)> insert_sz_func = [&](Alignment& a, Alignment& b){
        if (a.fragment_size() > 0 && a.mapping_quality() > 0
            && b.mapping_quality() >  0){
            int frag_diff = abs(a.fragment(0).length() - floor(insert_size));
            // Get mappings of node from graph

            // Get position (using Mapping* and path index)
            INS_INTERVAL i;
            i.start = srpe.ff.node_to_position[ a.path().mapping(0).position().node_id()] ;
            i.end = srpe.ff.node_to_position[ b.path().mapping(0).position().node_id()] ;
            i.fragl_supports +=1;
            // Set start to final position of first mate
            // Set end position to start + (frag_diff)
            // Length to frag_diff
            // Increment supports
            ins.push_back(i);
        }
    };

    std::vector<Alignment> splits;
    std::vector<INS_INTERVAL> split_intervals;
    std::function<void(Alignment&)> clip_read_func = [&](Alignment& a){
        // remap in clipped reads
        string a_clip = srpe.ff.get_clipped_seq(a);
        //cerr << a.DebugString() << endl;
        //cerr << a_clip << endl;
        if (a_clip.length() > 25){
            //cerr << "Remapping " << a_clip << endl;
            BREAKPOINT b;
            b.start = srpe.ff.get_clipped_position(a);
            b.split_supports += 1;
            bps.push_back(b);
            //vector<Alignment> a_remap = srpe.ff.remap(a_clip);
            //if (a_remap.size() > 0 && a_remap[0].score() > 0){
                //cerr << (a_remap[0].score() > 0 ? "MAPPED!" : "unmapped :(") << endl;
                // INS_INTERVAL i;
                // i.start = srpe.ff.get_clipped_position(a);
                // i.end = srpe.ff.node_to_position[ a_remap[0].path().mapping(0).position().node_id() ];
                // i.split_supports +=1;
                // ins.push_back(i);
                // cerr << i.to_string() << "\t" << a_remap.size() << endl;
                
            //}
        }
        
        // Set putative breakpoint evidence based on those (i.e. make interval for each split read)
        // Merge overlapping intervals.
        // Search for additional support in split reads.

    };



   
    // Once merged, we can output putative breakpoints
    // using the combined information.

    // Get splits that map near the breakpoint


    // INSERTIONS
    // Open relevant GAMfiles
    insert_stream.open(discordant_fn);
    clipped_stream.open(clipped_fn);
    if (!insert_stream.good()){
        cerr << "Must provide a .discordant file.";
    }
    else{
        // Collect all long insert reads
        stream::for_each_interleaved_pair_parallel(insert_stream, insert_sz_func);
        stream::for_each_parallel(clipped_stream, clip_read_func);
    }

    std::function<void(Alignment&, Alignment&)> oea_func = [&](Alignment& a, Alignment& b){
        BREAKPOINT bp;
        if (srpe.ff.soft_clip_filter(a).sequence() != ""){
            bp.start = srpe.ff.get_clipped_position(a);
        }
        else if (srpe.ff.soft_clip_filter(b).sequence() != ""){
            bp.start = srpe.ff.get_clipped_position(b);
        }
        else{
            if (a.read_mapped()){
                // take the right end of the read
                bp.start = srpe.ff.node_to_position[a.path().mapping( a.path().mapping_size() - 1).position().node_id() ];
            }
            else if (b.read_mapped()){
                // take the left end of the read
                bp.start = srpe.ff.node_to_position[a.path().mapping( b.path().mapping_size() - 1).position().node_id() ];
            }
        }
        bps.push_back(bp);
    };

    oea_stream.open(oea_fn);
    if (!oea_stream.good()){
        cerr << "No one-end-anchored file found." << endl
        << "Please provide one using [ vg sift -p -O x.gam ]" << endl;
        exit(1);
    }
    else{
        // Get all OEA reads
        // and place a breakpoint <insertsize> bp from the mapped read.
        // check if it's clipped
        stream::for_each_interleaved_pair_parallel(oea_stream, oea_func);
    }

    clipped_stream.open(clipped_fn);
    if(!clipped_stream.good()){
    }
    else{
        // collect some clipped reads
        // Split-read map their clipped portions
        // check if the seq mapped; if not, move on
        // if it did:
        //   find the interval where we see the break
        //   add some SR evidence
        //   refine its start / end / CI
    }

    

     // Merge insertion intervals and supports
    vector<INS_INTERVAL> merged = merge_intervals(ins);
    vector<BREAKPOINT> merged_bps = merge_breakpoints(bps);

    cerr << "Found " << merged.size() << " candidate intervals." << endl;
    for (auto x : merged){
        cerr << x.to_string() << endl;
    }
    cerr << "Found " << bps.size() << " candidate breakpoints." << endl;
    int count = 0;
    for (auto x : merged_bps){
        if (++count > 100){
            break;
        }
        cerr << x.to_string() << endl;
    }
     
    return 0;
}

static Subcommand vg_srpe ("srpe", "graph-external SV detection", main_srpe);

