#include <iostream>
#include <vector>
#include <getopt.h>
#include "stream.hpp"
#include "index.hpp"
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
    string gca_lcp_name = "";

    int max_iter = 0;
    int max_frag_len = -1;
    int max_softclip = -1;

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

            {0, 0, 0, 0}

        };
        int option_index = 0;
        c = getopt_long (argc, argv, "hx:",
                long_options, &option_index);

        // Detect the end of the options.
        if (c == -1)
            break;

        switch (c)
        {
            case 'x':
                xg_name = optarg;
                break;
            case 'h':
            case '?':
            default:
                help_srpe(argv);
                abort();
        }

    }


    SRPE srpe;
    vector<pair<Alignment, Alignment> > discords;
    std::function<void(Alignment&, Alignment&)> lambda = [&srpe, &discords](Alignment& aln_one, Alignment& aln_two){
        
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
        srpe.ff.set_inverse(true);
        pair<Alignment, Alignment> intermeds = srpe.ff.deletion_filter(aln_one, aln_two);
        if (aln_one.name() != ""){
            discords.push_back(intermeds);
        }
        srpe.ff.set_inverse(false);

    };

    
    /*
     * now that we've been through every read, check for erroneous depth signals
     * read the gam
     * the gam index
     * and the graph
     */


    gam_name = argv[optind];
    gam_index_name = argv[++optind];
    graph_name = argv[++optind];

    xg::XG xg_ind;
    Index gamind;
    
    if (!xg_name.empty()){
        ifstream in(xg_name);
        xg_ind.load(in);
    }
    if (!gam_index_name.empty()){
        gamind.open_read_only(gam_index_name);
    }
    else{

    }

    if (!graph_name.empty()){

    }
    else{

    }
    /***
     * First we need to find discordant Aligments
     * so that we can generate signatures for each SVTYPE.
     * We assume that our GAM contains these.
     */


    vector<Locus> sigs;
    std::function<void(Alignment&, Alignment&)> aln_to_sig_del = [&xg_ind, &sigs](Alignment& aln_one, Alignment& aln_two){
            // create a Path from the last non-clipped base in X and the first non-clipped base in X'
            // where X ---->     X' <-----
            // Include the match portion of the reads because why not.

#pragma omp critical
        {
            Locus ll;
            ll.set_name( aln_one.name() + "_L" );
            Position p_left;
            Position p_right;
            Path l_p;

            Path p_one = aln_one.path();
            Path p_two = aln_two.path();
            bool found_one = false;
            bool found_two = false;
            // Add first matches
            for (int i =  p_one.mapping_size() - 1; i >= 0; --i){
                Mapping m_curr = p_one.mapping(i);
                if (m_curr.edit_size() > 0){
                for (int j = m_curr.edit_size() - 1; j >= 0; --j){
                    Edit e = m_curr.edit(j);
                    if (e.from_length() != e.to_length()){
                        continue;
                    }
                    else{
                        p_left = m_curr.position();
                        Mapping* m_new = l_p.add_mapping();
                        *m_new->mutable_position() = p_left;
                        Edit* e_new = m_new->add_edit();
                        e_new->set_from_length(1);
                        e_new->set_to_length(1);
                        found_one = true;
                        break;
                    }
                }
                }
            }
                        //Add the second match portion.
                        //
            Mapping* m_new_right;
            for (int i = 0; i < p_two.mapping_size(); i++){
                Mapping m = p_two.mapping(i);
                if (m.edit_size() > 0){
                for (int j = m.edit_size() - 1; j >= 0; --j){
                    Edit e = m.edit(j);
                    if (e.from_length() == e.to_length()){
                        continue;
                    }
                    else{
                        p_right = m.position();
                        found_two = true;
                        //m_new_right = l_p.add_mapping();
                        //*m_new_right->mutable_position() = p_right;
                        //Edit* e_new = m_new_right->add_edit();
                        //e_new->set_from_length(1);
                        //e_new->set_to_length(1);
                        break;
                    }   
                }
                }
            }
            
            // Add a large edge via:
            // 1. A from_length the size of the deletion.
            // Get length of farthest clip -> nearest clip
            if (found_one && found_two){
                map<string, vector<size_t> > path_to_node_pos_left = xg_ind.node_positions_in_paths(p_left.node_id());
                map<string, vector<size_t> > path_to_node_pos_right = xg_ind.node_positions_in_paths(p_right.node_id());
                Path p;
                for (auto x_path : path_to_node_pos_left){
                    p = xg_ind.path(x_path.first);
                }
            // Add an edit with from_length of this length
            // and to_length = 0
            xg::size_t p_rank = xg_ind.path_rank(p.name());
            cerr << p.name() << " has been found." << endl;
            if ((p_left.node_id() < p_right.node_id()) && !p_left.is_reverse()){
                int64_t next_id = p_left.node_id();
                while (next_id != p_right.node_id()){
                    if (next_id == p_left.node_id()){

                    }
                    else{
                    Mapping mm = *l_p.add_mapping();
                    Position pp = *mm.mutable_position();
                    pp.set_node_id(next_id);
                    Edit* ee = mm.add_edit();
                    ee->set_to_length(0);
                    ee->set_from_length( xg_ind.node_length(next_id) );

                    cerr << next_id << " -> ";
                    }
                    next_id = xg_ind.prev_path_node_by_id(p_rank, next_id);

                }
            }
            cerr << endl;

            *ll.add_allele() = l_p;
           //err << ll.DebugString(); 
            sigs.push_back(ll);

            }
            else {
                cerr << "only one read contains a mapping." << endl;
            }
        }
    };


    /**
     * Next we will convert those alignments to Locus records.
     */


    if (!gam_name.empty()){
        ifstream gamfi(gam_name);
        gam_paired_interleaved_for_each_parallel(gamfi, aln_to_sig_del);
    }
    else{
        cerr << "NO GAM PROVIDED" << endl;
    }



    /**
     * We'll need to homogenize nearby records.
     */
    
    // TODO CHECK SIGS; if two within 10bp of start and 10 bp of end, homogenize them.


    /**
     * At this point we can emit preliminary calls,
     * or we can:
     *  1. Modify the graph by calling vg::edit on our newly created loci
     *  2. Remap reads around these modifications
     *  3. Re-call our loci
     *  4. Add in new variation / remove variation with low scores.
     */ 

    for (auto x : sigs){
    }

}
