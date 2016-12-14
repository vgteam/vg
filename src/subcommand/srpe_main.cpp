#include <iostream>
#include <vector>
#include <getopt.h>
#include <functional>
#include <regex>
#include "subcommand.hpp"
#include "stream.hpp"
#include "mapper.hpp"
#include "index.hpp"
#include "position.hpp"
#include "vg.pb.h"
#include "vg.hpp"
#include "srpe.hpp"
#include "filter.hpp"

using namespace std;
using namespace vg;
using namespace vg::subcommand;

void help_srpe(char** argv){
    cerr << "Usage: " << argv[0] << " srpe [options] <data.gam> <data.gam.index> <graph.vg>" << endl
        << "Options: " << endl
        
        << "   -S/ --specific <VCF>    look up variants in <VCF> in the graph and report only those." << endl
       << endl;
        //<< "-S / --SV-TYPE comma separated list of SV types to detect (default: all)." << endl
        


}



int main_srpe(int argc, char** argv){
    string gam_name = "";
    string gam_index_name = "";
    string graph_name = "";
    string xg_name = "";
    string gcsa_name = "";
    string lcp_name = "";

    string spec_vcf = "";

    int max_iter = 2;
    int max_frag_len = 10000;
    int min_soft_clip = 12;

    bool do_all = false;

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
            {"gcsa-index", required_argument, 0, 'g'},
            {"specific", required_argument, 0, 'S'},
            {"all", no_argument, 0, 'A'},
            {0, 0, 0, 0}

        };
        int option_index = 0;
        c = getopt_long (argc, argv, "hx:g:m:S:A",
                long_options, &option_index);

        // Detect the end of the options.
        if (c == -1)
            break;

        switch (c)
        {
            case 'm':
                max_iter = atoi(optarg);
                break;

            case 'A':
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
                optind++;
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

    auto var_t_from_path = [](){

    };

    auto get_nodes_from_path = [](list<Mapping> x){
        
    };

    regex is_alt ("_alt_.*");

    if (do_all){
        vector<Support> supports;
        for (auto x_path : (graph->paths)._paths){
            cerr << x_path.first << endl;
            if (regex_match(x_path.first, is_alt)){
                vector<Alignment> alns;
                vector<int64_t> var_node_ids;
                for (Mapping x_m : x_path.second){
                    var_node_ids.push_back(x_m.position().node_id()); 
                }
               
                std::function<void(const Alignment&)> incr = [&](const Alignment& a){
                    cerr << a.name() << endl;
                };
                gamind.for_alignment_to_nodes(var_node_ids, incr);
            }

        }
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
        srpe.ff.set_inverse(false);
        pair<Alignment, Alignment> intermeds = srpe.ff.deletion_filter(aln_one, aln_two);
        if (intermeds.first.name() != ""){
#pragma omp critical
            discords.push_back(intermeds);
        }

        srpe.ff.set_soft_clip_limit(min_soft_clip);
        intermeds = srpe.ff.soft_clip_filter(aln_one, aln_two);
        if (intermeds.first.name() != ""){
#pragma omp critical
            clippies.push_back(intermeds);
        }

        //TODO set path_length limit
        srpe.ff.max_path_length = 500;
        intermeds = srpe.ff.path_length_filter(aln_one, aln_two);
        if (intermeds.first.name() != ""){
#pragma omp critical
            fraggles.push_back(intermeds);
        }

        srpe.ff.set_inverse(false);

        
    };



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
    };

    auto se_score_func = [&](vector<Alignment> locals, Mapper* mp){
        double score = 0.0;
        for (auto xx : locals){
            score += (mp->align(xx.sequence())).score();
        }
        return score;
    };

    auto score_func = [&](vector<pair<Alignment, Alignment> > discordos, Mapper* mp){
        double score = 0.0;
        for (auto xx : discordos){
            Alignment tmp = mp->align(xx.first.sequence());
#pragma omp atomic update
            score += tmp.score();
            tmp = mp->align(xx.second.sequence());
#pragma omp atomic update
            score += tmp.score();
        }
        return score;
    };


    if (!gam_name.empty()){
        ifstream gamfi(gam_name);
        stream::for_each_interleaved_pair_parallel(gamfi, lambda);
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


    double max_aln_score = 0.0;
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
        if (clip_aln.path().mapping_size() == 0 || unclipped_aln.path().mapping_size() == 0){
            continue;
        }
        bool forward = clip_aln.path().mapping(0).position().node_id() > unclipped_aln.path().mapping(0).position().node_id() ? true : false;

        if (forward){

            for (int ti = 0; ti < unclipped_aln.path().mapping_size(); ti++){
                Mapping mi = unclipped_aln.path().mapping(ti);
                Mapping* mm_temp = add_me.add_mapping();
                *mm_temp = mi;
            }
            // Create edge if we find a SMEM

            //Mapping* mm = add_me.add_mapping();
            string onpath;
            if (true){
                //TODO: chance to dynamically find path name
                vector<xg::size_t> paths_of_clippy = xg_ind->paths_of_node(clipped_id);
                int iter = 0;
                int64_t next_path_node = clipped_id;
                while (iter < 5 && paths_of_clippy.empty()){
                    ++next_path_node;
                    paths_of_clippy = xg_ind->paths_of_node(next_path_node);
                }
                if (paths_of_clippy.size() > 1){
                    cerr << "Too many paths - SRPE only works with one path at a time. Exiting." << endl;
                    exit(1);
                }
                else if (paths_of_clippy.empty()){
                    cerr << "No path found within five nodes. Exiting." << endl;
                    exit(1);
                }
                xg::size_t p_rank = xg_ind->id_to_rank( paths_of_clippy[0]);
                onpath = paths_of_clippy[0];
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
            


            cerr << a.name() <<"\t" <<  a.sequence() << "\t" << a.score() << endl;
            //exit(1);
            // get subgraph, index, remap
            vg::VG* subg = new vg::VG();
            xg_ind->neighborhood(clipped_id == 1 ? clipped_id : clipped_id - 1, est_frag + 1000, subg->graph, false);
            //void expand_context(Graph& g, size_t dist, bool add_paths = true, bool use_steps = true,
            //                        bool expand_forward = true, bool expand_backward = true,
            //                                                int64_t until_node = 0) const;
            xg_ind->expand_context(subg->graph, 1, true, true, true, false, 0);
            // get_connected_nodes(Graph& g) const;
            //xg_ind->get_connected_nodes(subg.graph);

            subg->remove_orphan_edges();
            subg->rebuild_indexes();
            //subg->rebuild_mapping_aux();
            vector<Path> edit_mes;
            edit_mes.push_back(add_me);
            //cerr << add_me.DebugString() << endl;
            subg->edit(edit_mes);
            subg->compact_ids();

            subg->serialize_to_ostream(cout);

            gcsa::GCSA* sub_gcsa;
            gcsa::LCPArray* sub_lcp;
            xg::XG* sub_xg = new xg::XG();
            sub_xg->from_graph(subg->graph, false, false, false, false);
            int doubling_steps = 2;
            subg->build_gcsa_lcp(sub_gcsa, sub_lcp, 11, true, false, 2);
            Mapper* sub_mapper = new Mapper(sub_xg, sub_gcsa, sub_lcp);
            // our_mapper->align_paired_multi(aln1, aln2, queued_resolve_later, kmer_size, kmer_stride, max_mem_length, band_width, pair_window);
            bool qrl = false;
            vector<Alignment> score_me;
            vector<int64_t> p_node_ids;
            score_me.reserve(1000);
            xg::size_t clip_start = xg_ind->node_start(clipped_id);
            int64_t stop_id = xg_ind->node_at_path_position("HPV16",2000 + clip_start + est_frag + 1000);
            auto addr = [&score_me](const Alignment& aln){
                score_me.push_back(aln);
            };
            /*for (int64_t jj = clipped_id; jj < 35; jj++){
              vector<Alignment> tmp;
              gamind.for_alignments(jj, tmp);
              cerr << tmp.size() << endl;
              score_me.insert(score_me.end(), tmp.begin(), tmp.end());
              }*/
            vector<int64_t> ns;
            for (int i = 0; i < 35; i++){
                ns.push_back(i);
            }
            gamind.for_alignment_to_nodes(ns, addr);
            //gamind.get_alignments(clipped_id, xg_ind->node_at_path_position("HPV16", clip_start + est_frag + 1000), score_me);
            ////void for_alignment_in_range(int64_t id1, int64_t id2, std::function<void(const Alignment&)> lambda);
            /** for (int n_i = 0; n_i < 40; n_i++){
              p_node_ids.push_back((int64_t) n_i);
              vector<Alignment> tmp;
              gamind.get_alignments(n_i, tmp);
              cerr << tmp.size() << endl;
              score_me.insert(score_me.end(), tmp.begin(), tmp.end());
              }*/
            //double realn_score = score_func(clippies, sub_mapper);
            cerr << "Score_me size: " << score_me.size() << endl;
            double realn_score = se_score_func(score_me, sub_mapper);
            if (max_aln_score < realn_score){
                max_aln_score = realn_score;
            }
            cerr << max_aln_score << endl;;
            Alignment realn = sub_mapper->align(a.sequence());
            //cerr << realn.score() << endl;


            subg->serialize_to_ostream(cout);


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
    }
    return 0;
}

static Subcommand vg_srpe ("srpe", "graph-external SV detection", main_srpe);

