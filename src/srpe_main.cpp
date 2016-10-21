#include <iostream>
#include <vector>
#include <getopt.h>
#include <functional>
#include "stream.hpp"
#include "mapper.hpp"
#include "index.hpp"
#include "position.hpp"
#include "vg.pb.h"
#include "vg.hpp"
#include "srpe.hpp"

using namespace std;
using namespace vg;

void help_srpe(char** argv){
    cerr << "Usage: " << argv[0] << " srpe [options] <data.gam> <data.gam.index> <graph.vg>" << endl
        << "Options: " << endl
        << endl;

    /**
     * --sv-type
     * --threads
     * --max-iter
     * --realign
     */

}



int main_srpe(int argc, char** argv){
    string gam_name = "";
    string gam_index_name = "";
    string graph_name = "";
    string xg_name = "";
    string gcsa_name = "";
    string lcp_name = "";

    int max_iter = 2;
    int max_frag_len = 10000;
    int min_soft_clip = 12;

    /*
     * int sc_cutoff
     * int max_dist between reads
     *
     */
    vector<string> search_types;
    search_types.push_back("DEL");

    if (argc <= 2) {
        help_srpe(argv);
        return 1;
    }

    int c;
    optind = 2; // force optind past command positional argument
    while (true) {
        static struct option long_options[] =
        {
            {"xg-index", required_argument, 0, 'x'},
            {"help", no_argument, 0, 'h'},
            {"gcsa-index", required_argument, 0, 'g'},
            {0, 0, 0, 0}

        };
        int option_index = 0;
        c = getopt_long (argc, argv, "hx:g:",
                long_options, &option_index);

        // Detect the end of the options.
        if (c == -1)
            break;

        switch (c)
        {
            case 'x':
                xg_name = optarg;
                break;
            case 'g':
                gcsa_name = optarg;
                break;
            case 'h':
            case '?':
            default:
                help_srpe(argv);
                abort();
        }

    }


    SRPE srpe;



    gam_name = argv[optind];
    gam_index_name = argv[++optind];
    graph_name = argv[++optind];

    xg::XG* xg_ind = new xg::XG();
    Index gamind;

    vg::VG* graph;

    if (!xg_name.empty()){
        ifstream in(xg_name);
        xg_ind->load(in);
    }
    if (!gam_index_name.empty()){
        gamind.open_read_only(gam_index_name);
    }
    else{

    }

    if (!graph_name.empty()){
        ifstream in(graph_name);
        graph = new VG(in, false);
    }
    else{

    }
    /***
     * First we need to find discordant Aligments
     * so that we can generate signatures for each SVTYPE.
     * We assume that our GAM contains these.
     */
    vector<pair<Alignment, Alignment> > discords;
    vector<pair<Alignment, Alignment> > fraggles;
    vector<pair<Alignment, Alignment> > clippies;
    std::function<void(Alignment&, Alignment&)> lambda = [&](Alignment& aln_one, Alignment& aln_two){

        /* For each alignment in the gam:
         * Find the node the alignment maps to
         * find its mate read, if it has one
         * Check if alignment/mate fail any filters
         * if they do:
         *  find all the alignments that map to that node
         *  Look for matching signatures on other reads at the node
         *  Check the depth at the Locus, and update it as well
         *
         */
#pragma omp critical
        {
        srpe.ff.set_inverse(false);
        pair<Alignment, Alignment> intermeds = srpe.ff.deletion_filter(aln_one, aln_two);
        if (intermeds.first.name() != ""){
            discords.push_back(intermeds);
        }

        srpe.ff.set_soft_clip_limit(min_soft_clip);
        intermeds = srpe.ff.soft_clip_filter(aln_one, aln_two);
        if (intermeds.first.name() != ""){
            clippies.push_back(intermeds);
        }

        //TODO set path_length limit
        srpe.ff.max_path_length = 500;
        intermeds = srpe.ff.path_length_filter(aln_one, aln_two);
        if (intermeds.first.name() != ""){
            fraggles.push_back(intermeds);
        }

        srpe.ff.set_inverse(false);
    }
    };


    if (!gam_name.empty()){
        ifstream gamfi(gam_name);
        gam_paired_interleaved_for_each_parallel(gamfi, lambda);
    }
    else{
        cerr << "NO GAM PROVIDED" << endl;
    }


        gcsa::GCSA* gcsa_ind;
        gcsa_ind = new gcsa::GCSA();
        ifstream ifstr_gcsa (gcsa_name);
        gcsa_ind->load(ifstr_gcsa);

        lcp_name = gcsa_name + ".lcp";
        gcsa::LCPArray * lcp_ind;
        lcp_ind = new gcsa::LCPArray();
        ifstream ifstr_lcp (lcp_name);
        lcp_ind->load(ifstr_lcp);


        
        // Create a mapper
        Mapper* mapper;
        mapper = new Mapper(xg_ind, gcsa_ind, lcp_ind);
        


    for (int i = 0; i < clippies.size(); i++){
        Alignment a = clippies[i].first;


        Path add_me;

        int64_t clipped_id = 0;
        string clip = "";
        string unclip = "";
        if (a.path().mapping_size() > 0){
            Path path = a.path();
            Edit left_edit = path.mapping(0).edit(0);
            Edit right_edit = path.mapping(path.mapping_size() - 1).edit(path.mapping(path.mapping_size() - 1).edit_size() - 1);
            int left_overhang = left_edit.to_length() - left_edit.from_length();
            int right_overhang = right_edit.to_length() - right_edit.from_length();
            clipped_id = (left_overhang >= right_overhang) ? path.mapping(0).position().node_id() : path.mapping(path.mapping_size() - 1).position().node_id();
            clip = (left_overhang >= right_overhang) ? a.sequence().substr(a.sequence().length() - left_overhang, left_overhang) : a.sequence().substr(a.sequence().length() - right_overhang, right_overhang);
            unclip = (left_overhang >= right_overhang) ? a.sequence().substr(0, a.sequence().length() - left_overhang) : a.sequence().substr(0, a.sequence().length() - right_overhang);
            if (left_overhang >= right_overhang){
                //mapper->align(unclip);            
            }
            else{
                //Alignment unclipped_aln = mapper->align(unclip);
                //add_me = unclipped_aln.path();
            }
        }

        int est_frag = 0;
        std::function<void(const Alignment&)> frag_len_fil = [&](const Alignment& x){
            Alignment ret = srpe.ff.path_length_filter( (Alignment&) x);
            if (ret.name() != ""){
                for (int fi = 0; fi < ret.fragment_size(); fi++)
                est_frag = ret.fragment(i).length();
            }
        };

        // extract sot-clip sequences


        // add 500bp to fragment len as a margin of error.
        Alignment clip_aln = mapper->align(clip);

        Alignment unclipped_aln = mapper->align(unclip);

        bool forward = clip_aln.path().mapping(0).position().node_id() > unclipped_aln.path().mapping(0).position().node_id() ? true : false;

        if (forward){
        for (int ti = 0; ti < unclipped_aln.path().mapping_size(); ti++){
            Mapping mi = unclipped_aln.path().mapping(ti);
            Mapping* mm_temp = add_me.add_mapping();
            *mm_temp = mi;
        }
        // Create edge if we find a SMEM

        //Mapping* mm = add_me.add_mapping();
        if (true){
        //TODO: chance to dynamically find path name
        xg::size_t p_rank = xg_ind->id_to_rank( xg_ind->path_rank("HPV16") );
        //TODO should use most recently matched ID instead
        //*mm->mutable_position() = make_position(clipped_id, false, 0);
        int64_t next_id = clipped_id;
        while (next_id != 0 && next_id != clip_aln.path().mapping(0).position().node_id()){
            if (next_id != clipped_id){
                Mapping* m_next = add_me.add_mapping();
                *m_next->mutable_position() = make_position(next_id, false, 0);
                Edit* ee = m_next->add_edit();
                ee->set_to_length(0);
                ee->set_from_length( xg_ind->node_length(next_id) );
            }
            ++next_id;
            next_id = xg_ind->next_path_node_by_id(p_rank, next_id);
        }
        }
        for (int m_i = 0; m_i < clip_aln.path().mapping_size(); m_i++){
            Mapping* mm = add_me.add_mapping();
            Mapping clip_mapping = clip_aln.path().mapping(m_i);
            *mm = clip_mapping;
        }


        cerr << a.name() <<"\t" <<  a.sequence() << "\t" << a.score() << "\t" << clip_aln.DebugString() << "\t" << unclipped_aln.DebugString() << endl;
        //exit(1);
        // get subgraph, index, remap
        vg::VG subg;
        xg_ind->neighborhood(clipped_id, 10000, subg.graph, false);
        vector<Path> edit_mes;
        edit_mes.push_back(add_me);
        cerr << add_me.DebugString() << endl;
        subg.edit(edit_mes);

        gcsa::GCSA* sub_gcsa;
        gcsa::LCPArray* sub_lcp;
        xg::XG* sub_xg = new xg::XG(subg.graph);
        int doubling_steps = 2;
        subg.build_gcsa_lcp(sub_gcsa, sub_lcp, 11, true, false, 2);
        Mapper* sub_mapper = new Mapper(sub_xg, sub_gcsa, sub_lcp);
        Alignment realn = sub_mapper->align(a.sequence());
        cerr << realn.score() << endl;


        subg.serialize_to_ostream(cout);
   

        edit_mes.clear();
        }

        // create a Path from the last non-clipped base in X and the first non-clipped base in X'
        // where X ---->     X' <-----
        // Include the match portion of the reads because why not.




        //gamind.for_alignment_in_range(p.node_id() - 1, p.node_id() + 1, frag_len_fil);
        // gam_index.for_alignment_in_range();  // get all those alignments within the clipped region and 1 adj node each side,
        // and check if they have discordant frag_lens
        // if they do, get those fragment lens, determine the insert size,
        // and tuck an edge there. Remap that soft-clipped  alignments; expect their score to go up (hopefully).
        //
        // vector<MaximalExactMatch> find_smems(const string& seq, int max_length);
        //
        // 
    }



    /**
     * Next we will convert those alignments to Locus records.
     */

    }
