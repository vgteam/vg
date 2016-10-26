#include <iostream>
#include <vector>
#include <getopt.h>
#include <functional>
#include "stream.hpp"
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
        //<< "-S / --SV-TYPE comma separated list of SV types to detect (default: all)." << endl
        


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
    int max_soft_clip = -1;

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
            {"max-iter", required_argument, 0, 'm'},
            {"xg-index", required_argument, 0, 'x'},
            {"help", no_argument, 0, 'h'},

            {0, 0, 0, 0}

        };
        int option_index = 0;
        c = getopt_long (argc, argv, "m:x:",
                long_options, &option_index);

        // Detect the end of the options.
        if (c == -1)
            break;

        switch (c)
        {
            case 'm':
                max_iter = atoi(optarg);
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

    vg::VG* graph;

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
    std::function<void(Alignment&, Alignment&)> lambda = [&srpe, &discords, &clippies, &fraggles, &max_soft_clip](Alignment& aln_one, Alignment& aln_two){

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
        if (intermeds.first.name() != ""){
            discords.push_back(intermeds);
        }

        srpe.ff.set_soft_clip_limit(max_soft_clip);
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

        
    };


    gam_name = argv[optind];
    gam_index_name = argv[++optind];
    graph_name = argv[++optind];

    // Read in GAM and filter for Sigs
    
    // Convert sigs to edges and nodes

    // Find newly incorporated nodes / edges within X nodes or Y bp of each other
    // and remap reads to refine them.
    // Refinement: for N variant edges/nodes within X bp/nodes of one another, calculate the
    // sig quality ~ (reads mapped, mapping qual of reads, soft clipping, fragment length consistency,
    // )
    // Take the X highest-scoring edges / nodes and remap reads to just these.
    // Repeat this process until convergence, or after max-iter iterations

    // Can now emit refined variant calls.

    for (int i = 0; i < clippies.size(); i++){
        Alignment a = clippies[i].first;
        int64_t clipped_id = 0;
        if (a.path().mapping_size() > 0){
            Path path = a.path();
            Edit left_edit = path.mapping(0).edit(0);
            Edit right_edit = path.mapping(path.mapping_size() - 1).edit(path.mapping(path.mapping_size() - 1).edit_size() - 1);
            int left_overhang = left_edit.to_length() - left_edit.from_length();
            int right_overhang = right_edit.to_length() - right_edit.from_length();
            clipped_id = (left_overhang >= right_overhang) ? path.mapping(0).position().node_id() : path.mapping(path.mapping_size() - 1).position().node_id();
        }

        int est_frag = 0;
        std::function<void(const Alignment&)> frag_len_fil = [&](const Alignment& x){
            Alignment ret = srpe.ff.path_length_filter( (Alignment&) x);
            if (ret.name() != ""){
                for (int fi = 0; fi < ret.fragment_size(); fi++)
                est_frag = ret.fragment(i).length();
            }
        };

        // Create a mapper
        //
        // find SMEMs
        //
        // Create edge
        //
        // get subgraph, index, remap


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




    vector<Locus> sigs;
    std::function<void(Alignment&, Alignment&)> aln_to_sig_del = [&xg_ind, &sigs](Alignment& aln_one, Alignment& aln_two){
        // create a Path from the last non-clipped base in X and the first non-clipped base in X'
        // where X ---->     X' <-----
        // Include the match portion of the reads because why not.

        //#pragma omp critical
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
            //Mapping* m_new_right;
            Mapping m;
            for (int i = 0; i < p_two.mapping_size(); i++){
                m = p_two.mapping(i);
                if (m.edit_size() > 0){
                    for (int j = m.edit_size() - 1; j >= 0; --j){
                        Edit e = m.edit(j);
                        if (e.from_length() != e.to_length()){
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
                //cerr << p_left.node_id() << " <- " << p_right.node_id() << endl;
                int64_t low_pos = (p_left.node_id() < p_right.node_id()) ? p_left.node_id() : p_right.node_id();
                int64_t high_pos = (p_left.node_id() > p_right.node_id()) ? p_left.node_id() : p_right.node_id();
                //           map<string, vector<size_t> > path_to_node_pos_left = xg_ind.node_positions_in_paths(p_left.node_id());
                //           map<string, vector<size_t> > path_to_node_pos_right = xg_ind.node_positions_in_paths(p_right.node_id());
                //           Path p;
                //           for (auto x_path : path_to_node_pos_left){
                //               p = xg_ind.path(x_path.first);
                //           }
                // Add an edit with from_length of this length
                // and to_length = 0
                vector<size_t> p_of_n = xg_ind.paths_of_node(low_pos);
                int64_t p_id = xg_ind.rank_to_id(p_of_n[0]);
                xg::size_t p_rank = xg_ind.id_to_rank(p_id);
                //cerr << xg_ind.path_name(p_rank) << " has been found." << endl;

                int64_t next_id = low_pos;
                while ((next_id != high_pos) && (next_id != 0)){
                    if (next_id != low_pos){
                        Mapping* mm = l_p.add_mapping();
                        Position pp;
                        pp.set_node_id(next_id);
                        *mm->mutable_position() = make_position(next_id, false, 0);
                        Edit* ee = mm->add_edit();
                        ee->set_to_length(0);
                        ee->set_from_length( xg_ind.node_length(next_id) );

                    }

                    //cerr << next_id << " -> ";
                    ++next_id;
                    next_id = xg_ind.prev_path_node_by_id(p_rank, next_id);

                }

                Mapping* far_right = l_p.add_mapping();
                *far_right = m;
                //cerr << l_p.DebugString() << endl;
                //cerr << endl;

                *ll.add_allele() = l_p;
                //err << ll.DebugString(); 
#pragma omp critical
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
        gam_paired_interleaved_for_each_parallel(gamfi, lambda);
    }
    else{
        cerr << "NO GAM PROVIDED" << endl;
    }


    std::function<double(vector<pair<Alignment, Alignment> >)> score_paired_alignments = [](vector<pair<Alignment, Alignment> > alns){
        double MAX_MAPQ = 684;
        double MAX_SCORE = MAX_MAPQ * alns.size();
        double score = 0.0;

        double sum_mapq = 0.0;
        double sum_softlcip = 0.0;
        double cost_inverted_readpair = 2.0;
        double cost_one_end_anchored = 4.0;

        int ploidy = 2;
        double ploidy_score = 0.0 ;

        for (auto a : alns){
            sum_mapq += a.first.mapping_quality();
            sum_mapq += a.second.mapping_quality();
            if (a.first.path().mapping_size() == 0){

            }
            else{

            }
            if (a.second.path().mapping_size() == 0){

            }
            else{

            }
        }

        return (score / MAX_SCORE);
    };

    std::function<double(vector<Alignment>)> score_single_alignments = [](vector<Alignment> alns){
        double MAX_MAPQ = 684;
        double MAX_SCORE = MAX_MAPQ * alns.size();
        double score = 0.0;

        double sum_mapq = 0.0;
        double sum_softlcip = 0.0;

        int ploidy = 2;
        double ploidy_score = 0.0 ;

        for (auto a : alns){
            sum_mapq += a.mapping_quality();
            sum_mapq += a.mapping_quality();
            if (a.path().mapping_size() == 0){

            }
            else{

            }
        }

        return (score / MAX_SCORE);

    };



    vector<Path> paths_to_add;
    std::function<void(vector<Locus>)> homogenize_loci = [&paths_to_add](vector<Locus> loci){

    };

    std::function<void(int64_t, double&)> get_local_depth_at_node = [](int64_t node_id, double& avg_d){

    };

    //std::function<void(int64_t, int64_t, double&) frac_homozygous = [](int64_t start, int64_t end, double& score{
    //
    //};

    std::function<void(int64_t, int64_t, double&)> depth_in_range = [](int64_t front, int64_t back, double& avg_d){};

    double best_score = 0;
    for (auto ll : sigs){
        // Add that sig
        paths_to_add.push_back(ll.allele(0));
        //graph->edit(paths_to_add);
        ///remap_subgraph();
        // Score that sig
        //score = score_alns_paired();
        // Is it the best sig?
        // if (score > best_score){
        //  // Pop the old sig
        //
        //  // Push the new sig
        //  // best_score = score
        // }
        // else if (score == best_score){
        //   I guess include both?
        // }
        // else{
        //  continue;
        // } 
    }

    graph->edit(paths_to_add);

    graph->serialize_to_ostream(cout);
    delete graph;


    /**
     * We'll need to homogenize nearby records.
     * Homogenization process:
     *  1. for all variation w/in <WINDOW> bp:
     *      - Get subgraph
     *      - Remove all other variation in this window.
     *      - Map reads over nodes in thei region (and unmapped reads!)
     *      - Score it
     *  Take the variant(s) with the highest score; no idea how to decide this.
     *  Maybe just take the top N within some limit of each other.
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

    //for (auto x : sigs){
    //cerr << x.name() << endl;
    //}

}
