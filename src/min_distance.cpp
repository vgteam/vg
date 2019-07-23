//#define debugIndex
//#define debugDistance

#include "min_distance.hpp"

using namespace std;
namespace vg {

/*TODO: Remove old distance index from vg index 
 * Also change how nodes are stored in chain - in case of loops/unary snarls -might not actually need this
 * Make snarls/chains represented by the node id in netgraph
 */
MinimumDistanceIndex::MinimumDistanceIndex(const HandleGraph* graph, 
                             const SnarlManager* snarl_manager,
                             int64_t cap) {
    /*Constructor for the distance index given a VG and snarl manager
    */
    
    #ifdef debugIndex
    assert(graph != nullptr);
    assert(snarl_manager!= nullptr);
    #endif

    min_node_id = graph->min_node_id();
    max_node_id = graph->max_node_id();

    primary_snarl_assignments.resize(max_node_id - min_node_id + 1);
    primary_snarl_ranks.resize(max_node_id - min_node_id + 1);
    sdsl::util::set_to_value(primary_snarl_assignments, 0);
    sdsl::util::set_to_value(primary_snarl_ranks, 0);

    secondary_snarl_assignments.resize(max_node_id - min_node_id + 1);
    secondary_snarl_ranks.resize(max_node_id - min_node_id + 1);
    sdsl::util::set_to_value(secondary_snarl_assignments, 0);
    sdsl::util::set_to_value(secondary_snarl_ranks, 0);
    has_secondary_snarl_bv.resize(max_node_id - min_node_id + 1);
    sdsl::util::set_to_value(has_secondary_snarl_bv, 0);

    chain_assignments.resize(max_node_id - min_node_id + 1);
    chain_ranks.resize(max_node_id - min_node_id + 1);
    sdsl::util::set_to_value(chain_assignments, 0);
    sdsl::util::set_to_value(chain_ranks, 0);
    has_chain_bv.resize(max_node_id - min_node_id + 1);
    sdsl::util::set_to_value(has_chain_bv, 0);

    tree_depth = 0;


    #ifdef debugIndex
        cerr << endl << "Creating distance index"<< endl;
    #endif

    //Calculate minimum distance index
    const vector<const Snarl*> top_snarls = snarl_manager->top_level_snarls();

    unordered_set<const Snarl*> seen_snarls;
    for (const Snarl* snarl : top_snarls) {
       //Make an index for each disconnected snarl/chain
 
       if (seen_snarls.count(snarl) == 0){
          if (snarl_manager->in_nontrivial_chain(snarl)){
              const Chain* chain = snarl_manager->chain_of(snarl);
              calculateMinIndex(graph, snarl_manager, chain, 
                     0, false, false, 0);
              for (auto s : *chain) {
                  seen_snarls.insert(s.first);
              }
           } else {
               Chain curr_chain;
               curr_chain.emplace_back(snarl, false);
               calculateMinIndex(graph, snarl_manager, &curr_chain,
                            0, false, true, 0);
               seen_snarls.insert(snarl);
           }
       }
    }

    #ifdef debugIndex
    //Every node should be assigned to a snarl
    for (SnarlIndex& si : snarl_indexes) {
        assert (tree_depth >=
    }
    auto check_assignments = [&](const handle_t& h)-> bool {
        id_t id = graph->get_id(h); 
        assert( primary_snarl_assignments[id - min_node_id] != 0);
        assert( primary_snarl_ranks[id - min_node_id] != 0);

        if( secondary_snarl_assignments[id - min_node_id] != 0) {
            assert( secondary_snarl_ranks[id - min_node_id] != 0);
        }

        if( chain_assignments[id - min_node_id] != 0) {
            assert( chain_ranks[id - min_node_id] != 0);
        }

        handle_t handle = graph->get_handle(id, false);
        if ( snarl_indexes[
                   primary_snarl_assignments[id - min_node_id]-1].nodeLength(
                   primary_snarl_ranks[id-min_node_id]-1) 
               !=  graph->get_length(handle)
              ){
            cerr << id << ": predicted " << snarl_indexes[
                   primary_snarl_assignments[id - min_node_id]-1].nodeLength(
                   primary_snarl_ranks[id-min_node_id]-1)  << " actual " << graph->get_length(handle);
        }
        assert( snarl_indexes[
                   primary_snarl_assignments[id - min_node_id]-1].nodeLength(
                   primary_snarl_ranks[id-min_node_id]-1) 
               ==  graph->get_length(handle)
              );
        return true;
               
    };
    graph->for_each_handle(check_assignments);
    #endif


    util::assign(has_secondary_snarl, 
                 rank_support_v<1>(&has_secondary_snarl_bv));
    util::assign(has_chain, rank_support_v<1>(&has_chain_bv));

    //Remove empty entries of chain/secondary snarl assignments and ranks
    int_vector<> filtered_secondary_assignments 
                  (has_secondary_snarl.rank(max_node_id-min_node_id)+1, 0);

    size_t i = 0;
    for (auto x : secondary_snarl_assignments) {
        if (x != 0) {
            filtered_secondary_assignments[i] = x;
            i++;
        }
    }
    secondary_snarl_assignments = move(filtered_secondary_assignments);

    int_vector<> filtered_secondary_ranks 
                  (has_secondary_snarl.rank(max_node_id-min_node_id)+1, 0); 

    i = 0;
    for (auto x : secondary_snarl_ranks) {
        if (x != 0) {
            filtered_secondary_ranks[i] = x;
            i++;
        }
    }
    secondary_snarl_ranks = move(filtered_secondary_ranks);

    int_vector<> filtered_chain_assignments 
                    (has_chain.rank(max_node_id-min_node_id)+1, 0);
    i = 0;
    for (auto x : chain_assignments) {
        if (x != 0) {
            filtered_chain_assignments[i] = x;
            i++;
        }
    }
    chain_assignments = move(filtered_chain_assignments);

    int_vector<> filtered_chain_ranks 
                    (has_chain.rank(max_node_id-min_node_id)+1, 0);
    i = 0;
    for (auto x : chain_ranks) {
        if (x != 0) {
            filtered_chain_ranks[i] = x;
            i++;
        }
    }
    chain_ranks = move(filtered_chain_ranks);

    util::bit_compress(primary_snarl_assignments);
    util::bit_compress(primary_snarl_ranks);
    util::bit_compress(secondary_snarl_assignments);
    util::bit_compress(secondary_snarl_ranks);
    util::bit_compress(chain_assignments);
    util::bit_compress(chain_ranks);
    util::bit_compress(has_chain_bv);
    util::bit_compress(has_secondary_snarl_bv);


    if (cap > 0) {
        include_maximum = true;
        calculateMaxIndex(graph, cap);
    } else {
        include_maximum = false;
    }

    #ifdef debugIndex
    if (include_maximum) {
        //Every node should have a min and max dist to source 
        auto check_max = [&](const handle_t& h)-> bool {
            id_t id = graph->get_id(h); 
                cerr << id << endl;
            assert(max_distances[id-min_node_id] > 0);
            assert(min_distances[id-min_node_id] > 0);
            assert (max_distances[id-min_node_id] >= 
                    min_distances[id-min_node_id]);
            return true;
        };
        graph->for_each_handle(check_max);
    }
    #endif
};

MinimumDistanceIndex::MinimumDistanceIndex () {
    // Nothing to do
}

MinimumDistanceIndex::MinimumDistanceIndex (istream& in) : MinimumDistanceIndex() {
    // Load the index
    load(in);
}

  
void MinimumDistanceIndex::load(istream& in){
    //Load serialized index from an istream
    size_t num_snarls;
    sdsl::read_member(num_snarls, in);
    snarl_indexes.reserve(num_snarls);
    for (size_t i = 0 ; i < num_snarls ; i++) {
        snarl_indexes.emplace_back(); 
        snarl_indexes.back().load(in);
    }
    primary_snarl_assignments.load(in); 
    primary_snarl_ranks.load(in);
    secondary_snarl_assignments.load(in);
    secondary_snarl_ranks.load(in);
    has_secondary_snarl_bv.load(in);
    has_secondary_snarl.load(in);
    util::assign(has_secondary_snarl, 
                 rank_support_v<1>(&has_secondary_snarl_bv));

    //Load serialized chains
    size_t num_chains;
    sdsl::read_member(num_chains, in);
    chain_indexes.reserve(num_chains);
    for (size_t i = 0 ; i < num_chains ; i++) {
        chain_indexes.emplace_back();
        chain_indexes.back().load(in);
    }
    chain_assignments.load(in);
    chain_ranks.load(in);
    has_chain_bv.load(in);
    has_chain.load(in);
    util::assign(has_chain, 
                 rank_support_v<1>(&has_chain_bv));

    sdsl::read_member(min_node_id, in );
    sdsl::read_member(max_node_id, in );

    sdsl::read_member(tree_depth, in );

    sdsl::read_member(include_maximum, in );

    if (include_maximum) {
        min_distances.load(in);
        max_distances.load(in);
    }



};

void MinimumDistanceIndex::serialize(ostream& out) const {

    //Serialize snarls
    sdsl::write_member(snarl_indexes.size(), out);
    
    for (auto& snarl_index: snarl_indexes) {
        snarl_index.serialize(out);
    }
    primary_snarl_assignments.serialize(out);
    primary_snarl_ranks.serialize(out);
    secondary_snarl_assignments.serialize(out);
    secondary_snarl_ranks.serialize(out);
    has_secondary_snarl_bv.serialize(out);
    has_secondary_snarl.serialize(out);

    //Serialize chains 
    sdsl::write_member(chain_indexes.size(), out);
    
    for (auto& chain_index: chain_indexes) {
        chain_index.serialize(out);
    }
    chain_assignments.serialize(out);
    chain_ranks.serialize(out);
    has_chain_bv.serialize(out);
    has_chain.serialize(out);

    sdsl::write_member(min_node_id, out);
    sdsl::write_member(max_node_id, out);

    sdsl::write_member(tree_depth, out);

    sdsl::write_member(include_maximum, out);
    if (include_maximum) {
        min_distances.serialize(out);
        max_distances.serialize(out);
    }

};

/////////////////////////    MINIMUM INDEX    ///////////////////////////////


int64_t MinimumDistanceIndex::calculateMinIndex(const HandleGraph* graph,
                                    const SnarlManager* snarl_manager,
                                    const Chain* chain, size_t parent_id,
                                    bool rev_in_parent, bool trivial_chain, 
                                    size_t depth) {
    /*Populate the MinimumDistanceIndex
     * Compute the ChainIndex for this chain and recursively calculate the 
     * SnarlIndexes for all snarls within the chain
     * parentId is the id of this chain's parent snarl where the snarl is the 
     * node's primary snarl
     * trivialChain is true if the chain is really just a single snarl
    */

    #ifdef debugIndex
        cerr << "starting ";
        if (trivial_chain) { cerr << "snarl at ";}
        else {cerr << "chain at ";}
        cerr << get_start_of(*chain) << endl;
    #endif
    tree_depth = std::max(depth, tree_depth);

    auto cmp = [] (pair<pair<id_t, bool>,int64_t> x,
                                           pair<pair<id_t, bool>,int64_t> y) {
        //Comparison function for the priority of a pair of handle, distance
        return (x.second > y.second);
    };

 
    if (!trivial_chain) {
        //If this is a chain, initialize a new ChainIndex object

        //Get the start of the chain
        auto first_visit = get_start_of(*chain);
        chain_indexes.emplace_back(parent_id, first_visit.node_id(),
                                   rev_in_parent,
                                   first_visit.node_id() 
                                            == get_end_of(*chain).node_id(),
                                   chain->size());

        chain_assignments[first_visit.node_id()-min_node_id] = 
                                                       chain_indexes.size();
        chain_ranks[first_visit.node_id()-min_node_id] = 1;
        has_chain_bv[first_visit.node_id()-min_node_id] = 1; 

        handle_t first_node = graph->get_handle(first_visit.node_id(), 
                                           first_visit.backward());
        chain_indexes.back().prefix_sum[0] = graph->get_length(first_node) + 1;
    }
    size_t curr_chain_assignment = chain_indexes.size() - 1;
    size_t curr_chain_rank = 0;

    ChainIterator c_end = chain_end(*chain);
    for (ChainIterator c = chain_begin(*chain); c != c_end; ++c) {
        //for each snarl in the chain 

        const Snarl* snarl = c->first;
        bool snarl_rev_in_chain = c->second;

        id_t snarl_start_id = snarl->start().node_id();
        bool snarl_start_rev = snarl->start().backward(); //into snarl
        id_t snarl_end_id = snarl->end().node_id();
        bool snarl_end_rev = snarl->end().backward();   //pointing out
        //Id of boundary node that occurs second in the chain
        id_t second_id = snarl_rev_in_chain ? snarl_start_id : snarl_end_id;
            #ifdef debugIndex
                cerr << "At snarl " << snarl_start_id << " with rank "
                    << curr_chain_rank << " in chain. Snarl starts: " << snarl->start() << " and ends at " << snarl->end() << endl;
            #endif

        if (!trivial_chain && chain_assignments[second_id - min_node_id] 
                               == 0){
            //Store the index of the start of the snarl only if it hasn't
            //already been seen (if the chain loops)
            chain_assignments[second_id-min_node_id] = curr_chain_assignment+1;
            chain_ranks[second_id - min_node_id] = curr_chain_rank + 2;
            has_chain_bv[snarl_end_id - min_node_id] = 1;
           
        } 

        NetGraph ng = NetGraph(snarl->start(), snarl->end(), 
                                snarl_manager->chains_of(snarl), graph);

        //Get all the nodes in the snarl
        //TODO: Make this a vector. Make sure net graph for each handle only traverses each node once in unary snarls
        //
        hash_set<pair<id_t, bool>> all_nodes;

        size_t snarl_assignment = snarl_indexes.size();
        auto add_node = [&](const handle_t& h)-> bool {
            id_t id = ng.get_id(h); 
            if (id != snarl_start_id && id != snarl_end_id) {
                const Snarl* temp_snarl = snarl_manager->into_which_snarl(
                                          id, false);
                const Snarl* curr_snarl = temp_snarl == NULL ? 
                      snarl_manager->into_which_snarl(id, true) :
                      temp_snarl; 
                if (curr_snarl != NULL) {
                    //If this node represents a snarl or chain, then this snarl
                    //is a secondary snarl
                    has_secondary_snarl_bv[id-min_node_id] = 1;
                    secondary_snarl_assignments[id - min_node_id] 
                                                          = snarl_assignment+1;
                    secondary_snarl_ranks[id - min_node_id]= all_nodes.size()+1;
                } else {
                    //Otherwise this is the node's primary snarl
                    primary_snarl_assignments[id-min_node_id] = 
                                                            snarl_assignment+1;
                    primary_snarl_ranks[id - min_node_id] = all_nodes.size()+1;
                }
                all_nodes.emplace(id, false);
                all_nodes.emplace(id, true);
            }
            return true;
                   
        };

        //Put all visits in the snarl into a vector, ensuring that the
        //inward start and end visits are at the beginning and end of the list
        all_nodes.emplace(snarl_start_id, snarl_start_rev);
        all_nodes.emplace(snarl_start_id, !snarl_start_rev);
        ng.for_each_handle(add_node);
        all_nodes.emplace(snarl_end_id, snarl_end_rev);
        all_nodes.emplace(snarl_end_id, !snarl_end_rev);


        id_t start_in_chain = snarl_rev_in_chain ? snarl_end_id : snarl_start_id; 
        id_t end_in_chain = snarl_rev_in_chain ? snarl_start_id : snarl_end_id; 
        //Assign the second boundary node (relative to the chain) to this snarl
        //This will replace the first node in a chain if the chain loops
        primary_snarl_assignments[start_in_chain-min_node_id] = 
                                                           snarl_assignment + 1;
        if (start_in_chain == snarl_start_id) {
            primary_snarl_ranks[start_in_chain-min_node_id] = 
                                                        snarl_start_rev ? 2 : 1;
        } else {
            primary_snarl_ranks[start_in_chain-min_node_id] = snarl_end_rev ? 
                                    all_nodes.size() : all_nodes.size() - 1;
        }
        if (primary_snarl_assignments[end_in_chain-min_node_id] == 0 ){
            //If the 2nd boundary node doesn't already have a primary snarl,
            //then assign it to this snarl
            primary_snarl_assignments[end_in_chain-min_node_id] = 
                                                             snarl_assignment+1;
            primary_snarl_ranks[end_in_chain-min_node_id] = 
                 end_in_chain == snarl_end_id ?  
                 (snarl_end_rev ? all_nodes.size() : all_nodes.size() - 1) :
                 (snarl_start_rev ? 2 : 1);
        } 
        if (!trivial_chain &&
             secondary_snarl_assignments[end_in_chain-min_node_id] == 0){
            //Otherwise, assign the first boundary node a secondary snarl
            secondary_snarl_assignments[end_in_chain-min_node_id] = 
                                                            snarl_assignment+1;
            secondary_snarl_ranks[end_in_chain-min_node_id] = 
                 end_in_chain == snarl_end_id ? 
                 (snarl_end_rev ? all_nodes.size()  : all_nodes.size() - 1) :
                 (snarl_start_rev ? 1 : 0);
            has_secondary_snarl_bv[end_in_chain-min_node_id] = 1;
        }

        //Make the snarl index
        if (trivial_chain) {
            //The parent is the parent snarl
            snarl_indexes.emplace_back(parent_id, rev_in_parent, 
                           snarl_start_id, snarl_start_id == snarl_end_id, 
                           depth, all_nodes.size()/2, false);
        } else {
            //The parent is the chain
            snarl_indexes.emplace_back(get_start_of(*chain).node_id(), 
                               snarl_rev_in_chain, start_in_chain, 
                               snarl_start_id == snarl_end_id,
                               depth, all_nodes.size()/2, true);
        }


        for (pair<id_t, bool> start_id : all_nodes){

            //Use each node in the snarl as start of djikstra search

            //Index of the start node in the current snarl
            size_t start_node_rank = 
                   primary_snarl_assignments[start_id.first - min_node_id]-1
                                                      == snarl_assignment
                 ? primary_snarl_ranks[start_id.first-min_node_id] - 1
                 : secondary_snarl_ranks[start_id.first-min_node_id]-1;


            if (start_id.second) {
                start_node_rank = start_node_rank % 2 == 0 ? start_node_rank + 1
                                                          : start_node_rank - 1;
            }
            handle_t start_handle = 
                           graph->get_handle(start_id.first, start_id.second);
            //Priority queue of reachable nodes (pair of node id and direction)
            priority_queue<  pair<pair<id_t, bool>, int64_t>,  
                         vector<pair<pair<id_t, bool>, int64_t>>, 
                                     decltype(cmp)> reachable(cmp);
            reachable.push(make_pair(start_id, 0));

            #ifdef debugIndex
                cerr << "  Start Node: " << start_id.first << "," 
                    << start_id.second << endl;
                assert( primary_snarl_assignments[start_id.first - min_node_id]-1 == snarl_assignment ||  secondary_snarl_assignments[start_id.first - min_node_id]-1 == snarl_assignment);
            #endif
            bool first_loop = true;
            unordered_set<pair<id_t, bool>> seen_nodes;

            while (reachable.size() > 0) {
                pair<pair<id_t, bool>, int64_t> next = reachable.top();
                reachable.pop();
                pair<id_t, bool> curr_id = next.first;
                handle_t curr_handle = ng.get_handle(curr_id.first, 
                                                        curr_id.second);
                int64_t curr_dist = next.second;
                if ( seen_nodes.count(curr_id) == 0) {
                    //If node has not already been found:

                    //Record distance from start to current node 
                    if (!first_loop) {
                        size_t curr_node_rank = 
                          primary_snarl_assignments[curr_id.first-min_node_id]-1
                                                       == snarl_assignment
                          ? primary_snarl_ranks[curr_id.first-min_node_id] -1
                          : secondary_snarl_ranks[curr_id.first-min_node_id]-1;

                        if (curr_id.second) {
                            curr_node_rank = curr_node_rank % 2 == 0 
                                           ? curr_node_rank + 1
                                           : curr_node_rank - 1;
                        }

                        snarl_indexes[snarl_assignment].insertDistance(
                                    start_node_rank, curr_node_rank, curr_dist);
                        seen_nodes.insert(curr_id);

                    }

                    
                    int64_t node_len; //length of the current node
                       
                    int64_t loop_dist = -1;
                         //Dist to enter curr node then exit at same side 

                    //Get the snarl that the node represents, if any
                    const Snarl* temp_snarl = snarl_manager->into_which_snarl(
                                              curr_id.first, curr_id.second);
                    const Snarl* curr_snarl = temp_snarl == NULL ? 
                          snarl_manager->into_which_snarl(curr_id.first, 
                                                          !curr_id.second) :
                          temp_snarl; 

                    if (curr_id.first != snarl_start_id &&
                           curr_id.first != snarl_end_id && curr_snarl != NULL) {
                        //If current node is a child snarl/chain


                        if (snarl_manager->in_nontrivial_chain(curr_snarl)) {
                           //The node is a chain

                            const Chain* curr_chain= snarl_manager->chain_of(
                                                                    curr_snarl);
                            size_t chain_start = get_start_of(*curr_chain).node_id();

                            if (chain_assignments[chain_start-min_node_id]!= 0){
                                //Length of chain has already been found
                                ChainIndex& chain_dists = chain_indexes[
                                  chain_assignments[chain_start-min_node_id]-1];
                               
                                //Get the length of the node (chain)
                                node_len = chain_dists.chainLength();
       
                                //Get loop dist- enter and exit chain at same side
                                if (get_start_of(*curr_chain).backward() 
                                                            == curr_id.second) {
                                    //If traversing chain forward in snarl

                                    loop_dist = chain_dists.loop_fd[0] - 1;

                                    if (loop_dist != -1) {
                                        auto visit = get_start_of(*curr_chain);
                                        handle_t temp_handle= graph->get_handle(
                                                             visit.node_id(),
                                                             visit.backward());

                                       loop_dist = loop_dist + 
                                               graph->get_length(temp_handle) ;
                                     }

                                } else {
                                    loop_dist = chain_dists.loop_rev[
                                             chain_dists.loop_rev.size()-1] - 1;

                                    if (loop_dist != -1) { 

                                        auto end_visit= get_end_of(*curr_chain);
                                        handle_t temp_handle= graph->get_handle(
                                                   end_visit.node_id(),
                                                   end_visit.backward());
                                        loop_dist = loop_dist + 
                                           graph->get_length(temp_handle);
                                     }
                                }

                            } else {//haven't recursed on this chain yet
                                #ifdef debugIndex
                                    cerr << " recurse" << endl;
                                #endif
                                bool rev_in_snarl = curr_id.first ==
                                            get_start_of(*curr_chain).node_id()
                                          ? get_start_of(*curr_chain).backward()
                                          : !get_end_of(*curr_chain).backward();
                                node_len = calculateMinIndex(graph, 
                                             snarl_manager, curr_chain, 
                                             start_in_chain, rev_in_snarl,
                                             false, depth + 1);

                                ChainIndex& curr_chain_dists = chain_indexes[
                                  chain_assignments[chain_start-min_node_id]-1];
                                if (get_start_of( *curr_chain).backward()
                                                 == curr_id.second) {
                                    //If traversing snarl forward in chain

                                    loop_dist = curr_chain_dists.loop_fd[0] - 1;

                                    if (loop_dist != -1) {
                                        auto visit = get_start_of(*curr_chain);
                                        handle_t temp_handle= graph->get_handle(
                                                             visit.node_id(),
                                                             visit.backward());

                                       loop_dist = loop_dist + 
                                               graph->get_length(temp_handle) ;
                                    }
                                } else {

                                    loop_dist = curr_chain_dists.loop_rev[
                                        curr_chain_dists.loop_rev.size()-1] - 1;

                                    if (loop_dist != -1) {
                                        auto end_visit = get_end_of(*curr_chain);
                                        handle_t temp_handle = graph->get_handle(
                                                     end_visit.node_id(),
                                                     end_visit.backward());
                                        loop_dist = loop_dist +
                                           graph->get_length(temp_handle);
                                     }
                                } 
                            }
                        } else {//Snarl

                            id_t snarl_id = curr_snarl->start().node_id();
                            bool snarl_rev = curr_snarl->start().backward();
                            id_t end_id = curr_snarl->end().node_id();
                            bool end_rev = curr_snarl->end().backward();
  

                            if (primary_snarl_assignments[snarl_id-min_node_id]
                                != 0) {
                                //Already found
                                SnarlIndex& snarl_dists = snarl_indexes[
                                                  primary_snarl_assignments[
                                                       snarl_id-min_node_id]-1];
                                node_len = snarl_dists.snarlLength();

                                //Find the distance to enter and exit snarl
                                //at the same side
                                if (curr_id.second == snarl_rev) { 
                                    //If traversing snarl forward
                                    loop_dist = snarl_dists.snarlDistance(0, 1);

                                     if (loop_dist != -1) { 
                                         handle_t temp_handle =
                                              graph->get_handle(
                                                curr_snarl->start().node_id(),
                                                curr_snarl->start().backward());
                                         loop_dist = loop_dist
                                             + 2*graph->get_length(temp_handle);
                                     }
                                } else {
                                    size_t end_in = snarl_dists.is_unary_snarl ?
                                              0 : snarl_dists.num_nodes * 2 - 1;
                                    size_t end_out= snarl_dists.is_unary_snarl ?
                                              1 : snarl_dists.num_nodes * 2 - 2;
                                    loop_dist = snarl_dists.snarlDistance(
                                             end_in, end_out);

                                     if (loop_dist != -1) {
                                         handle_t temp_handle =
                                              graph->get_handle(
                                                  curr_snarl->end().node_id(),
                                                  curr_snarl->end().backward());
                                         loop_dist = loop_dist + 
                                                2*graph->get_length(temp_handle);
                                     }
                                }
                            } else {//Haven't recursed on snarl yet
                                #ifdef debugIndex
                                    cerr << " recurse" << endl;
                                #endif
                                
                                //Create chain to recurse on and recurse
                                Chain curr_chain;

                                curr_chain.emplace_back(curr_snarl, false);
                                bool rev_in_snarl = curr_id.first == snarl_id 
                                              ? snarl_rev 
                                              : !end_rev;
                                calculateMinIndex(graph, snarl_manager,
                                                 &curr_chain, start_in_chain,
                                                 rev_in_snarl, true, depth + 1);

                                SnarlIndex& curr_snarl_dists = snarl_indexes[
                                          primary_snarl_assignments[
                                                snarl_id-min_node_id]-1];

                                node_len = curr_snarl_dists.snarlLength(); 

                                //Find the distance to enter and exit snarl
                                //at the same side
                                if (curr_id.second == snarl_rev) {

                                    loop_dist = curr_snarl_dists.snarlDistance(0, 1);

                                    handle_t temp_handle = 
                                          graph->get_handle(
                                             curr_snarl->start().node_id(),
                                             curr_snarl->start().backward());
                                     if (loop_dist != -1) { 
                                         loop_dist = loop_dist 
                                             + 2*graph->get_length(curr_handle);
                                     }

                                 } else {

                                    size_t end_in = 
                                         curr_snarl_dists.is_unary_snarl ?
                                         0 : curr_snarl_dists.num_nodes * 2 - 1;
                                    size_t end_out = 
                                             curr_snarl_dists.is_unary_snarl ?
                                         1 : curr_snarl_dists.num_nodes * 2 - 2;
                                     loop_dist = curr_snarl_dists.snarlDistance(
                                               end_in, end_out);

                                     if (loop_dist != -1) { 
                                         handle_t temp_handle = 
                                               graph->get_handle(
                                                  curr_snarl->end().node_id(),
                                                  curr_snarl->end().backward());
                                         loop_dist = loop_dist + 
                                               2*graph->get_length(temp_handle);
                                      }
                                 }
                            }
                                        
                        }
                    } else { //Node is just a node
                        node_len = graph->get_length(curr_handle);
                    }
 
                    if (curr_id == start_id) {
                        snarl_indexes[snarl_assignment].distances[
                                            start_node_rank/2]  = node_len + 1; 
                    }
       

                    if (loop_dist != -1 && !first_loop) {
                        /*If there is a path within the current node that loops 
                          to enter the node and exit it at the same side - add
                          reachable nodes from current node in reverse 
                          Do not add this distance if the current node is the 
                          starting node */

                        handle_t rev_handle = ng.get_handle(
                                          ng.get_id(curr_handle), 
                                          !ng.get_is_reverse(curr_handle)); 
                            

                        auto add_rev_handle = [&](const handle_t& h)-> bool {
                            pair<id_t, bool> node = make_pair(
                                        ng.get_id(h), ng.get_is_reverse(h));
                            reachable.push(make_pair(node, 
                                                     curr_dist + loop_dist));
 

                             return true;
                        };

                        ng.follow_edges(rev_handle, false, add_rev_handle);
                    }

                    //Add reachable nodes to priority queue
                    auto add_handle = [&](const handle_t& h)-> bool {
                         pair<id_t, bool> node = make_pair(
                                    ng.get_id(h), ng.get_is_reverse(h));
                       if (node_len != -1) {
                       reachable.push(make_pair(node, curr_dist + node_len));
                       }
                      
                         #ifdef debugIndex
                             cerr << node.first << " " << node.second << ", ";
                         #endif
                         return true;
                    };
                    //Add reachable nodes to priority queue for unary snarl that doesn't loop - 0 distance
                    auto add_handle0 = [&](const handle_t& h)-> bool {
                         pair<id_t, bool> node = make_pair(
                                    ng.get_id(h), ng.get_is_reverse(h));
                       reachable.push(make_pair(node, 0));
                       
                         #ifdef debugIndex
                             cerr << node.first << " " << node.second << ", ";
                         #endif
                         return true;
                    };


#ifdef debugIndex

     cerr << "    From start node " << start_id.first << " " << start_id.second 
        << " in snarl " << snarl_indexes[snarl_assignment].id_in_parent
        << " at " << ng.get_id(curr_handle) << " " << ng.get_is_reverse(curr_handle) << endl; 
     cerr << "        Adding next nodes:  ";
#endif
                    if ((node_len == -1 && first_loop) || curr_id == start_id) {
                        //If the nodeLen is -1 then node is a unary snarl that 
                        //doesn't have a path from start to end. If this is the
                        //start of the distance calculation then add subsequent
                        //nodes assuming that the node length was 0
                        //Or if this is the starting node

                        ng.follow_edges(curr_handle, false, add_handle0);

                    } else  {

                        ng.follow_edges(curr_handle, false, add_handle);
                    }  
                        

                    //Add edges between the boundary nodes that are not in 
                    //the net graph
                    int64_t next_dist = curr_id == start_id ? 0 
                                                           : curr_dist+node_len;

                    if ((curr_id.first == snarl_start_id &&
                        curr_id.second != snarl_start_rev) ||
                         ( curr_id.first == snarl_end_id &&
                                 curr_id.second == snarl_end_rev )  ) {
                           
                        //If currently leaving the snarl
                        auto add_handle_end = [&](const handle_t& h)-> bool {
                            pair<id_t, bool> node = make_pair(
                                    ng.get_id(h), ng.get_is_reverse(h));
                             if ( node.first == snarl_start_id || 
                                  node.first == snarl_end_id ) {
                               reachable.push(make_pair(node, next_dist));
                            }
                            return true;
                        };
                        graph->follow_edges(curr_handle, false, add_handle_end);

                    }                      
#ifdef debugIndex
     cerr << "    prev dist: " << curr_dist << "+ new dist " << node_len << endl;
#endif
                } 
                first_loop = false;
            }//End while loop
        }//End for loop over starting node/directions in a snarl
#ifdef debugIndex
    cerr << "End snarl " << snarl_indexes[snarl_assignment].id_in_parent << endl;
#endif

        if (!trivial_chain) {
            // Add to prefix sum the distance to the beginning and end of the 
            // last node in the current snarl
            SnarlIndex& sd = snarl_indexes[snarl_assignment];

            size_t num_nodes = snarl_indexes[snarl_assignment].num_nodes;
            int64_t dist = snarl_rev_in_chain
                ? sd.snarlDistance(num_nodes * 2 - 1, 1) 
                                      + sd.nodeLength(num_nodes * 2 - 1)
                : sd.snarlDistance( 0, num_nodes * 2 - 2) + sd.nodeLength(0);
            #ifdef debugIndex
                cerr << "Prefix sum before snarl: " 
                << chain_indexes[curr_chain_assignment].prefix_sum[
                                                curr_chain_rank + 1] << endl;
            #endif
            chain_indexes[curr_chain_assignment].prefix_sum[curr_chain_rank+1] =
                               curr_chain_rank == 0 ? dist + 1:
                                chain_indexes[curr_chain_assignment].prefix_sum[
                                                        curr_chain_rank]+dist;


            //Add the reverse loop distance
            if ( curr_chain_rank == 0) {
               //If this is the first snarl in the chain, get the loop distance
               //of the first node
               int64_t first_rev_dist;
               if (snarl_rev_in_chain){ 
                    first_rev_dist = sd.snarlDistance( sd.num_nodes * 2 - 2,
                                                       sd.num_nodes * 2 - 1);
                    first_rev_dist = first_rev_dist == -1 ? -1 :
                          first_rev_dist + sd.nodeLength(sd.num_nodes * 2 - 2);
                } else {
                    first_rev_dist = sd.snarlDistance( 1, 0);
                    first_rev_dist = first_rev_dist == -1 ? -1 :
                                     first_rev_dist + sd.nodeLength(0);
                }
                chain_indexes[curr_chain_assignment].loop_rev[0] = first_rev_dist + 1;
            }

            int64_t rev_loop_dist;
            if ( snarl_rev_in_chain ) {
    
                rev_loop_dist = sd.snarlDistance(0, 1);
                rev_loop_dist = rev_loop_dist == -1 ? -1 : 
                                     rev_loop_dist + sd.nodeLength(0);
            } else {
                rev_loop_dist = sd.snarlDistance(sd.num_nodes * 2 - 1,
                                                 sd.num_nodes * 2 - 2);
                rev_loop_dist = rev_loop_dist == -1 ? -1 : 
                            rev_loop_dist + sd.nodeLength(sd.num_nodes * 2 - 2);
            }
     
    
            //Loop distance of the previous node
            int64_t last_loop = chain_indexes[curr_chain_assignment].loop_rev[
                                            curr_chain_rank] - 1;
            if (last_loop == -1) {
    
                chain_indexes[curr_chain_assignment].loop_rev[curr_chain_rank+1]
                            = rev_loop_dist + 1;
    
            } else {
    
                //Push the minimum of the loop distance of the current snarl and
                //the loop distance of the previous snarl + dist to and from loop 
                int64_t dist_to_end = sd.snarlDistance(0, sd.num_nodes * 2 - 2);
                dist_to_end = dist_to_end == -1 ? -1 
                            : dist_to_end + dist_to_end + sd.nodeLength(0) 
                                    + sd.nodeLength(sd.num_nodes * 2 - 1);


                int64_t loop_distance = minPos({rev_loop_dist, 
                                                last_loop + dist_to_end});
               chain_indexes[curr_chain_assignment].loop_rev[curr_chain_rank+1] = loop_distance + 1;
            }
        }
        
        //Bit compress distance matrix of snarl index
        util::bit_compress(snarl_indexes[snarl_assignment].distances);

        curr_chain_rank ++;
    }//End for loop over snarls in chain

    if (!trivial_chain){
        //Get the distances for loops in the chain
        ChainIndex& cd = chain_indexes[curr_chain_assignment];

        //Add the length of the last node to chain prefix sum
        auto last_visit = get_end_of(*chain);
        handle_t last_node = graph->get_handle(last_visit.node_id(), false);
        cd.prefix_sum[cd.prefix_sum.size() - 1] = 
               cd.prefix_sum[cd.prefix_sum.size() - 2] + graph->get_length(last_node);
    
        if (get_start_of(*chain).node_id() == get_end_of(*chain).node_id()) {
            //If the chain loops, then the reverse loop distances might include
            //looping through the chain
            size_t curr_chain_rank = 0;
            for (ChainIterator c = chain_begin(*chain); c != c_end; ++c) {
                //Loop through the chain forward 
                const Snarl* snarl = c->first; 
                bool snarl_rev_in_chain = c->second;
    
                //Snarl is the primary snarl of the first node in the chain
                auto& sd = snarl_rev_in_chain ? 
                  snarl_indexes[getPrimaryAssignment(snarl->start().node_id())] :
                  snarl_indexes[getPrimaryAssignment(snarl->end().node_id())];
                int64_t new_loop;
                if (curr_chain_rank == 0) {
                    new_loop = cd.loop_rev[cd.loop_rev.size() - 1] - 1;
                } else {
                    int64_t prev_loop = cd.loop_rev[curr_chain_rank - 1] - 1;
                    int64_t dist_to_end = sd.snarlDistance(0, sd.num_nodes * 2 - 2);
                    dist_to_end = dist_to_end == -1 ? -1 
                                : dist_to_end + dist_to_end + sd.nodeLength(0) 
                                        + sd.nodeLength(sd.num_nodes * 2 - 1);
    
                    new_loop = prev_loop  == -1 ? -1 : prev_loop + dist_to_end;
                }

                if (cd.loop_rev[curr_chain_rank] == 0 ||
                     new_loop < cd.loop_rev[curr_chain_rank]) {
                    cd.loop_rev[curr_chain_rank] = new_loop + 1;
                } else {
                    //If this isn't an improvement, subsequent entries will not
                    //be improved either
                    break;
                }
                curr_chain_rank ++;

            }
        }
        //Add forward loop distances 
       
        //Check if there is an edge traversing last node in chain fd -> rev 
    
        curr_chain_rank = chain->size();
        ChainIterator chain_start_r = chain_rend(*chain);
        for (ChainIterator c = chain_rbegin(*chain); c != chain_start_r; ++c) {
            //Loop through the chain in reverse
            const Snarl* snarl = c->first; 
            bool snarl_rev_in_chain = c->second;
            id_t snarl_start_id = snarl->start().node_id();
            id_t snarl_end_id = snarl->end().node_id();

            //Snarl is the primary snarl of the first node in the chain
            auto& sd = snarl_rev_in_chain ? 
               snarl_indexes[primary_snarl_assignments[snarl_end_id-min_node_id]-1] :
               snarl_indexes[primary_snarl_assignments[snarl_start_id-min_node_id]-1];
            NetGraph ng (snarl->start(), snarl->end(),
                          snarl_manager->chains_of(snarl), graph);
    
                                          
    
            if (c == chain_rbegin(*chain)) {
                //If this is the last snarl in the chain, push loop for last node
    
                int64_t loop_dist_last; 
                if (snarl_rev_in_chain) {
           
                    loop_dist_last = sd.snarlDistance( 1, 0 );
                    loop_dist_last = loop_dist_last == -1 ? -1 : 
                                     loop_dist_last + sd.nodeLength(0);
                } else {
    
                    loop_dist_last = sd.snarlDistance(sd.num_nodes * 2 - 2,
                                                      sd.num_nodes * 2 - 1);
                    loop_dist_last = loop_dist_last == -1 ? -1 : 
                           loop_dist_last + sd.nodeLength(sd.num_nodes * 2 - 2);
                }
    
                if (get_start_of(*chain).node_id() 
                       == get_end_of(*chain).node_id()) {
                    //If the chain loops, might need distance from first snarl

                    ChainIterator chain_start = chain_begin(*chain);
                    const Snarl* first_snarl = chain_start->first;
                    bool first_snarl_rev = chain_start->second;
                        
                    SnarlIndex& sd_first = snarl_rev_in_chain ? 
                           snarl_indexes[primary_snarl_assignments[snarl_start_id-min_node_id]-1] :
                           snarl_indexes[primary_snarl_assignments[snarl_end_id-min_node_id]-1];
                    if (first_snarl_rev) {
                        int64_t new_dist = 
                              sd_first.snarlDistance(sd_first.num_nodes * 2 - 1,
                                                    sd_first.num_nodes * 2 - 2);
                        new_dist = new_dist == -1 ? -1 : 
                           new_dist + sd_first.nodeLength(sd_first.num_nodes * 2 - 2);
                        loop_dist_last = minPos({loop_dist_last, new_dist });
                    } else {
                        int64_t new_dist = sd_first.snarlDistance(0, 1);
                        new_dist = new_dist == -1 ? -1 : 
                                              new_dist + sd_first.nodeLength(0);
                        loop_dist_last = minPos({loop_dist_last, new_dist });
   
                    }
                  
                }
                cd.loop_fd[curr_chain_rank] = loop_dist_last + 1;
            }
    
            int64_t fd_loop_dist;
    
    
            if (snarl_rev_in_chain) {
                //If the snarl is reversed in the chain
                fd_loop_dist = sd.snarlDistance(sd.num_nodes * 2 - 1, 
                                                sd.num_nodes * 2 - 2);
                fd_loop_dist = fd_loop_dist == -1 ? -1 :
                           fd_loop_dist + sd.nodeLength(sd.num_nodes*2-1);
            } else {
                fd_loop_dist = sd.snarlDistance(0, 1);
                fd_loop_dist = fd_loop_dist == -1 ? -1 :
                                               fd_loop_dist + sd.nodeLength(0);
            }
    
            int64_t last_loop = cd.loop_fd[curr_chain_rank] - 1;
            curr_chain_rank--;
    
            if (last_loop == -1) {
    
                cd.loop_fd[curr_chain_rank] = fd_loop_dist + 1;
    
            } else {
            //push dist to end of snarl + loop dist + dist to start of snarl 
    
                int64_t dist_end_start = 
                                  sd.snarlDistance( sd.num_nodes * 2 - 1, 1);
                dist_end_start = dist_end_start == -1 ? -1 :
                        dist_end_start + sd.nodeLength(sd.num_nodes * 2 - 1);


                int64_t dist_start_end = 
                               sd.snarlDistance(0, sd.num_nodes * 2 - 2);
                dist_start_end = dist_start_end == -1 ? -1 :
                        dist_start_end + sd.nodeLength(0);
                int64_t loop_distance = minPos({fd_loop_dist, 
                                  last_loop + dist_end_start + dist_start_end});
                cd.loop_fd[curr_chain_rank] = loop_distance + 1;
            }           
          
        }
        util::bit_compress(cd.prefix_sum);
        util::bit_compress(cd.loop_fd);
        util::bit_compress(cd.loop_rev);
    }
 
    //return length of entire chain
    
    return trivial_chain ? 
       snarl_indexes[primary_snarl_assignments[
                  get_start_of(*chain).node_id()-min_node_id]-1].snarlLength() :
       chain_indexes[curr_chain_assignment].prefix_sum[chain_indexes[
                              curr_chain_assignment].prefix_sum.size() - 1] - 1;
};



//////////////////    Distance Calculations
//

int64_t MinimumDistanceIndex::maxDistance(pos_t pos1, pos_t pos2) {
    if (!include_maximum) {
        return -1;
    } else {
        id_t id1 = get_id(pos1);
        id_t id2 = get_id(pos2);
        int64_t len1 = snarl_indexes[getPrimaryAssignment(id1)].nodeLength(
                                                         getPrimaryRank(id1));
        int64_t len2 = snarl_indexes[getPrimaryAssignment(id2)].nodeLength(
                                                         getPrimaryRank(id2));

        len1 = std::max((int64_t)(get_offset(pos1)+1), (int64_t)(len1-get_offset(pos1)));
        len2 = std::max((int64_t)(get_offset(pos2)+1), (int64_t)(len2-get_offset(pos2)));

        id1 -= min_node_id;
        id2 -= min_node_id;

        int64_t max_dist = std::max(
                         (int) max_distances[id1] - (int)min_distances[id2], 
                         (int) max_distances[id2] - (int)min_distances[id1]);

        return len1 + len2 + max_dist;

    }
}

int64_t MinimumDistanceIndex::minDistance(pos_t pos1, pos_t pos2) {
    /*Minimum distance between positions not including the position itself*/
    
    int64_t shortest_distance = -1; 

    id_t node_id1 = get_id(pos1);
    bool node_rev1 = is_rev(pos1);
    id_t node_id2 = get_id(pos2); 
    bool node_rev2 = is_rev(pos2);

    if (node_id1 == node_id2 && node_rev1 == node_rev2 ) {
        //if positions are on the same node and strand
        int64_t offset1 = get_offset(pos1);
        int64_t offset2 = get_offset(pos2);

        if (offset1 <= offset2) {
            shortest_distance = offset2-offset1+1; //+1 to be consistent
        }

    }
    
    //Index into snarl_indexes/chain_indexes of the common ancestor snarl/chain
    //true if common ancestor is a chain
    pair<id_t, bool> common_ancestor ( 0, false);


#ifdef debugDistance
    cerr << endl << "Start distance calculation from " << pos1 << "->" <<
         pos2 << endl;

    cerr << "Shortes distance within same node: " << shortest_distance<<  endl;

    cerr << "Find common ancestor" << endl;
#endif


    //// Find common ancestor of the two snarls
    unordered_set<pair<id_t, bool>> ancestors1;

    //set of all ancestor snarls and chains (true) of node1
    pair<id_t, bool> ancestor1 ( snarl_indexes[getPrimaryAssignment(
                                         node_id1)].id_in_parent,
                                 false);

#ifdef debugDistance
    cerr << "Ancestors of 1: ";
#endif


    while (ancestor1.first != 0) {
#ifdef debugDistance
        if (ancestor1.second) {cerr << "chain ";}
        else {cerr << "snarl ";}
        cerr << ancestor1.first << " ";
#endif
        if (ancestor1.second) {
            //If ancestor1 is a chain
            size_t chain_assignment = getChainAssignment(ancestor1.first);
            ancestors1.emplace(chain_assignment, true );
            ancestor1 = make_pair(chain_indexes[chain_assignment].parent_id,
                                   false);
        } else {
            size_t snarl_assignment = getPrimaryAssignment(ancestor1.first);
            ancestors1.emplace(snarl_assignment, false );
            SnarlIndex& si = snarl_indexes[snarl_assignment];
            ancestor1 = make_pair(si.parent_id, si.in_chain);
        }
    }

#ifdef debugDistance
      cerr << endl << "ancestors of 2: ";
#endif


    pair<id_t, bool> ancestor2 ( snarl_indexes[getPrimaryAssignment(
                                      node_id2)].id_in_parent,
                                  false);
    while (ancestor2.first != 0) {

#ifdef debugDistance
        if (ancestor2.second) {cerr << "chain ";}
        else {cerr << "snarl ";}
        cerr << ancestor2.first << " ";
#endif

        if (ancestor2.second) {
            //If this ancestor is a chain
            size_t chain_assignment = getChainAssignment(ancestor2.first);
            if (ancestors1.count(make_pair(chain_assignment, true)) != 0) {
                common_ancestor = ancestor2;
                break;
            }
            ancestor2 = make_pair(chain_indexes[chain_assignment].parent_id,
                                  false);
        } else { 
            size_t snarl_assignment = getPrimaryAssignment(ancestor2.first);
            if (ancestors1.count(make_pair(snarl_assignment, false)) != 0) {
                common_ancestor = ancestor2;
                break;
            }
            SnarlIndex& si = snarl_indexes[snarl_assignment];
            ancestor2 = make_pair(si.parent_id, si.in_chain);
        }
    }
    if (common_ancestor.first == 0) {
        //If the two positions don't share a common ancestor
        return -1;
    }

#ifdef debugDistance 
    cerr << endl;
    if (common_ancestor.second) {
        cerr << "common ancestor chain ";
    } else { 
        cerr << "common ancestor snarl ";
    }
    cerr << common_ancestor.first << endl;
    
#endif

    pair<size_t, bool> ancestor ( common_ancestor.second 
                             ? getChainAssignment(common_ancestor.first)
                             : getPrimaryAssignment(common_ancestor.first),
                       common_ancestor.second);


    //Find distances from pos1 and pos2 to ends of child snarls of ancestor
    int64_t distL1; int64_t distR1; pair<id_t, bool> snarl_tree_node1;
    tie (distL1, distR1, snarl_tree_node1) = 
                                    distToCommonAncestor(ancestor, pos1, false);

     
    int64_t distL2; int64_t distR2; pair<id_t, bool> snarl_tree_node2;
    tie (distL2, distR2, snarl_tree_node2) = 
                                     distToCommonAncestor(ancestor, pos2, true);

    pair<id_t, bool> parent = common_ancestor;
    bool lowest_ancestor = true;
    while (parent.first != 0) {
        //snarl_tree_nodes 1 and 2 are children of parent snarl or chain 
        if (parent.second) {
            //If the parent is a chain and both snarl_tree_nodes are snarls

            ChainIndex& chain_index = chain_indexes[
                                           getChainAssignment(parent.first)]; 

            //If the two nodes are snarls in the common ancestor chain
            //find the distance between them in the chain
            SnarlIndex& snarl_index1 = snarl_indexes[getPrimaryAssignment(
                                 snarl_tree_node1.first)];
            size_t start_rank1 = getChainRank(snarl_index1.id_in_parent);
            size_t end_rank1 = start_rank1 + 1; 
            if (snarl_index1.rev_in_parent) {
                int64_t temp = distL1;
                distL1 = distR1;
                distR1 = temp;
            }
            int64_t start_len1 = snarl_index1.nodeLength(
                                  snarl_index1.rev_in_parent ?
                                            snarl_index1.num_nodes * 2 - 1 : 0);
            int64_t end_len1 = snarl_index1.nodeLength(snarl_index1.rev_in_parent ?
                                            0 : snarl_index1.num_nodes * 2 - 1);
            size_t start_rank2, end_rank2;
            int64_t start_len2, end_len2;

            if (lowest_ancestor) {
                //If this is the lowest common ancestor, then there are two
                //separate nodes that need to be found
               
               SnarlIndex& snarl_index2 = snarl_indexes[
                             getPrimaryAssignment( snarl_tree_node2.first)];
               start_rank2 = getChainRank(snarl_index2.id_in_parent);
               end_rank2 = start_rank2+1; 
               if (snarl_index2.rev_in_parent) {
                   int64_t temp = distL2;
                   distL2 = distR2;
                   distR2 = temp;
               }
               start_len2 = snarl_index2.nodeLength(snarl_index2.rev_in_parent ?
                                            snarl_index2.num_nodes * 2 - 1 : 0);
               end_len2 = snarl_index2.nodeLength(snarl_index2.rev_in_parent ?
                                               0 : snarl_index2.num_nodes * 2 - 1);
            } else {
                //Otherwise just copy the first one
                start_rank2 = start_rank1;
                end_rank2 = end_rank1;
                start_len2 = start_len1;
                end_len2 = end_len1;
            }

            //Distance from left of s1 (reverse), left of s2 (forward)
            int64_t d1 = chain_index.chainDistance(
                   make_pair(start_rank1, true), 
                   make_pair(start_rank2, false), start_len1, start_len2);
            d1 = (distL1 == -1 || distL2 == -1 || d1 == -1) ? -1 : 
                                    distL1 + distL2 + d1 - start_len1;

            //Distance from left of s1 (reverse) to right of s2 (reverse)
            int64_t d2;
            if (start_rank1 == end_rank2) {
                //If snarls share a node, then the distances up to this point
                //both include the length of the shared node
                d2 = (distL1 == -1 || distR2 == -1) ? -1 : 
                                       distL1 + distR2 - start_len1; 
            } else {
                d2 = chain_index.chainDistance(
                   make_pair(start_rank1, true), 
                   make_pair(end_rank2, true), start_len1, end_len2);
                d2 = (distL1 == -1 || distR2 == -1 || d2 == -1) ? -1 : 
                                       distL1 + distR2 + d2 - start_len1;
            }

            //Distance from right of s1 (fd) to left of s2 (fd)
            int64_t d3;
            if (end_rank1 == start_rank2) {
                d3 = (distR1 == -1 || distL2 == -1) ? -1 : 
                                       distR1 + distL2 - end_len1; 
            } else {
                d3 = chain_index.chainDistance(
                   make_pair(end_rank1, false), 
                   make_pair(start_rank2, false), end_len1, start_len2);

                d3 = (distR1 == -1 || distL2 == -1 || d3 == -1) ? -1 : 
                                       distR1 + distL2 + d3 - end_len1; 
            }

            //Distance from right of s1 (fd) to right of s2 (rev)
            int64_t d4 =  chain_index.chainDistance(
                   make_pair(end_rank1, false), 
                  make_pair(end_rank2, true), end_len1, end_len2);
            d4 = (distR1 == -1 || distR2 == -1 || d4 == -1) ? -1 : 
                                       distR1 + distR2 + d4 - end_len1;
            
                       
            shortest_distance = minPos({d1, d2, d3, d4, shortest_distance});


            //Extend distances of both positions to the ends of the chain

            //Find the rank of the end node and current node in the chain
            size_t chain_end_rank = chain_index.prefix_sum.size() - 2;
                                    
            //Get the lengths of start, end, and current node
            int64_t chain_start_len = snarl_indexes[getPrimaryAssignment(
                   chain_index.id_in_parent)].nodeLength(
                       getPrimaryRank(chain_index.id_in_parent));

            int64_t chain_end_len = chain_index.prefix_sum[chain_index.prefix_sum.size()-1] 
                                - chain_index.prefix_sum[chain_index.prefix_sum.size()-2];

            int64_t dsl = chain_index.chainDistance(make_pair(0, false), 
                              make_pair(start_rank1, false),
                                          chain_start_len, start_len1);
            int64_t dsr = chain_index.chainDistance(make_pair(0, false), 
                                make_pair(end_rank1, true),
                                          chain_start_len, end_len1);
            int64_t der = chain_index.chainDistance(
                                make_pair(chain_end_rank, true), 
                                make_pair(end_rank1, true),
                                           chain_end_len, end_len1);
            int64_t del = chain_index.chainDistance(make_pair(chain_end_rank, true),
                              make_pair(start_rank1, false),
                                          chain_end_len, start_len1);

            int64_t dsl1 = dsl == -1 || distL1 == -1? -1 : distL1 + dsl; 
            int64_t dsr1 =  dsr == -1 || distR1 == -1? -1 : distR1 + dsr; 
            int64_t der1 = der == -1 || distR1 == -1? -1 : distR1 + der; 
            int64_t del1 = del == -1 || distL1 == -1? -1 : distL1 + del; 
 
            distL1 = minPos({dsr1, dsl1});
            distR1 = minPos({der1, del1}); 

            if (lowest_ancestor) {
                //If the two snarl tree nodes are different, need to find 
                //distances to ends for the second one

                dsl = chain_index.chainDistance(make_pair(0, false), 
                                 make_pair(start_rank2, false),
                                  chain_start_len, start_len2);
                dsr = chain_index.chainDistance(make_pair(0, false),
                                  make_pair(end_rank2, true),
                                    chain_start_len, end_len2);
                der = chain_index.chainDistance(make_pair(chain_end_rank, true), 
                                 make_pair(end_rank2, true),
                                 chain_end_len, end_len2);
                del = chain_index.chainDistance(make_pair(chain_end_rank, true),
                                make_pair(start_rank2, false),
                                chain_end_len, start_len2);

            }

            int64_t dsl2 = dsl == -1 || distL2 == -1? -1 : distL2 + dsl; 
            int64_t dsr2 = dsr == -1 || distR2 == -1? -1 : distR2 + dsr; 
            int64_t der2 = der == -1 || distR2 == -1? -1 : distR2 + der; 
            int64_t del2 = del == -1 || distL2 == -1? -1 : distL2 + del; 

            distL2 = minPos({dsr2, dsl2});
            distR2 = minPos({der2, del2});

            snarl_tree_node1 = make_pair(chain_index.id_in_parent,
                                         chain_index.rev_in_parent);


#        ifdef debugDistance
            cerr << "At ancestor chain " << chain_indexes[
                     getChainAssignment(parent.first)].id_in_parent << endl;
            cerr << "Ranks: " << start_rank1 << " " << start_rank2 << endl;
            
            cerr << "  Distances within ancestor: " << d1 << ", " << d2
                                                << ", " << d3 << ", " << d4 << endl;
            cerr << "  Shortest dist: " << shortest_distance
                << endl;
            cerr << "  Distances to ends of ancestor: " << distL1 << " " << distR1
                 << " " << distL2 << " " << distR2 
                << endl;
#        endif
            parent = make_pair(chain_index.parent_id, false);
        } else {
            //The two nodes are in the parent snarl

            size_t parent_snarl_index = getPrimaryAssignment(parent.first);
            SnarlIndex& snarl_index = snarl_indexes[parent_snarl_index];

            size_t rank1 = getPrimaryAssignment(snarl_tree_node1.first) 
                             == parent_snarl_index
               ? getPrimaryRank(snarl_tree_node1.first) :
                  getSecondaryRank(snarl_tree_node1.first);
            size_t rev_rank1;
            if (snarl_tree_node1.second) {
                //If this node is reversed
                rev_rank1 = rank1;
                rank1 = rev_rank1 % 2 == 0 ? rev_rank1 + 1 : rev_rank1 - 1;
            } else {
                rev_rank1 = rank1 % 2 == 0 ? rank1 + 1 : rank1 - 1;
            }

            size_t rev_rank2 ;
            size_t rank2;
            if (lowest_ancestor) {
                rank2 = getPrimaryAssignment(snarl_tree_node2.first )
                                     == parent_snarl_index
                     ? getPrimaryRank(snarl_tree_node2.first): 
                      getSecondaryRank(snarl_tree_node2.first);
                if (snarl_tree_node2.second) {
                    //If this node is reversed
                    rev_rank2 = rank2;
                    rank2 = rev_rank2 % 2 == 0 ?  rev_rank2 + 1 : rev_rank2 - 1;
                } else {
                    rev_rank2 = rank2 % 2 == 0 ?  rank2 + 1 : rank2 - 1;
                }
            } else {
                rank2 = rank1;
                rev_rank2 = rev_rank1;
            }


            int64_t d1 = snarl_index.snarlDistance(rank1, rank2);
            d1 = (distR1 == -1 || distL2 == -1 || d1 == -1) ? -1 : 
                                                          distR1 + distL2 + d1; 

            int64_t d2 = snarl_index.snarlDistance(rank1, rev_rank2);

            d2 = (distR1 == -1 || distR2 == -1 || d2 == -1) ? -1 :
                                                           distR1 + distR2 + d2;
            int64_t d3 = snarl_index.snarlDistance(rev_rank1, rank2);
            d3 = (distL1 == -1 || distL2 == -1 || d3 == -1) ? -1 : 
                                                      distL1 + distL2 + d3; 
            int64_t d4 = snarl_index.snarlDistance(rev_rank1, rev_rank2);
            d4 = (distL1 == -1 || distR2 == -1 || d4 == -1) ? -1 : 
                                                          distL1 + distR2 + d4; 

            shortest_distance =  minPos({d1, d2, d3, d4, shortest_distance});

            
            
            //Find distances to the ends of the parent snarl
            tie (distL1, distR1) = snarl_index.distToEnds(rank1, distL1, distR1);

            tie (distL2, distR2) = snarl_index.distToEnds(rank2, distL2, distR2);
             

            snarl_tree_node1 = make_pair(snarl_index.id_in_parent,
                                         snarl_index.rev_in_parent);

#        ifdef debugDistance
            cerr << "At ancestor snarl " << snarl_indexes[getPrimaryAssignment( parent.first)].id_in_parent << endl;
            cerr << "  Distances within ancestor: " << d1 << ", " << d2
                                                << ", " << d3 << ", " << d4 << endl;
            cerr << "  Shortest dist: " << shortest_distance
                << endl;
            cerr << "  Distances to ends of ancestor: " << distL1 << " " << distR1
                 << " " << distL2 << " " << distR2 
                << endl;
#        endif
            parent = make_pair(snarl_index.parent_id, snarl_index.in_chain);

        }
        lowest_ancestor = false;
    }

    shortest_distance = shortest_distance == -1 ? -1 : shortest_distance - 1;
    return shortest_distance;

};

tuple<int64_t, int64_t, pair<id_t, bool>> MinimumDistanceIndex::distToCommonAncestor(
          pair<size_t, bool> common_ancestor, pos_t& pos, bool rev){

    /* Find the distance from pos to either end of a snarl tree node 
     * (snarl or chain) in common_ancestor snarl or chain
     * common_ancestor is the index of the start node of the snarl/chain
     * that is the common ancestor and true if the common ancestor is a chain
       Return the two distances and the snarl/chain whose parent is the is
       common_ancestor or commonAncestor if the position is on a node (not a
       snarl) in common_ancestor
    */

    int64_t distL; //Dist from pos1 to boundaries of curr snarl 
    int64_t distR; //To start and end of snarls, not necessarily left/right
    id_t node_id = get_id(pos); 
    
    size_t snarl_assignment = getPrimaryAssignment(node_id);
    size_t snarl_rank = getPrimaryRank(node_id);
    SnarlIndex& snarl_index = snarl_indexes[snarl_assignment];

    int64_t offset = get_offset(pos);
    #ifdef debugDistance
    cerr << "Finding the distance to the common ancestor " << endl;
    cerr << "  Dist to the "; 
    rev ? cerr << "end " : cerr << "start ";
    cerr << "of node " << get_id(pos) << " offset " << offset <<" reversed " 
         << is_rev(pos) << " in snarl " << snarl_index.id_in_parent << ": " ;
    
    #endif
    if (is_rev(pos)) {
        distR = offset+1;
        distL = snarl_index.nodeLength(snarl_rank)-offset;
    } else {
        distL = offset+1;
        distR = snarl_index.nodeLength(snarl_rank)-offset;
    }
    //Make sure that the distance calculation will be directed - must traverse
    //positions in the correct direction
    if (rev == is_rev(pos)) {
        distL = -1;
    } else {
        distR = -1;
    }
#ifdef debugDistance
cerr << distL << " " << distR << endl;
#endif

    bool node_rev = false;
    int64_t start_len = 0;
    int64_t end_len = 0;
    pair<size_t, bool> parent (snarl_assignment, false);
    bool at_node = true;
    while (parent != common_ancestor) {
        //Traverse the snarl tree up to the common ancestor and accumulate 
        //distances to the ends of each ancestor snarl/chain of the position
        if (parent.second) {
            //If the parent is a chain
            ChainIndex& chain_index = chain_indexes[parent.first];

            size_t start_rank = getChainRank(node_id);
            size_t end_rank = start_rank + 1; 
            if (node_rev) {
                //If the snarl is traversed backwards in the chain, then the
                //left and right distances must be switched
                int64_t temp = distL;
                distL = distR;
                distR = temp;
            }


            //Find the rank of the end node and current node in the chain
            size_t chain_end_rank = chain_index.prefix_sum.size() - 2;
                                    
            //Get the lengths of start, end, and current node

            int64_t chain_start_len = chain_index.prefix_sum[0];

            int64_t chain_end_len = 
                     chain_index.prefix_sum[chain_index.prefix_sum.size()-1] 
                    - chain_index.prefix_sum[chain_index.prefix_sum.size()-2];

            int64_t dsl = chain_index.chainDistance(make_pair(0, false), 
                                          make_pair(start_rank, false),
                                          chain_start_len, start_len);
            int64_t dsr = chain_index.chainDistance(make_pair(0, false), 
                                          make_pair(end_rank, true),
                                          chain_start_len, end_len);
            int64_t der =chain_index.chainDistance(
                                           make_pair(chain_end_rank, true),
                                           make_pair(end_rank, true),
                                           chain_end_len, end_len);
            int64_t del =chain_index.chainDistance(
                                          make_pair(chain_end_rank, true),
                                          make_pair(start_rank, false),
                                          chain_end_len, start_len);

            dsl = dsl == -1 || distL == -1? -1 : distL + dsl;
            dsr = dsr == -1 || distR == -1? -1 : distR + dsr;
            der = der == -1 || distR == -1? -1 : distR + der;
            del = del == -1 || distL == -1? -1 : distL + del;
 
            distL = minPos({dsr, dsl});
            distR = minPos({der, del});

            node_id = chain_index.id_in_parent;
            node_rev = chain_index.rev_in_parent;
            parent = make_pair(getPrimaryAssignment(chain_index.parent_id),
                               false);

    #ifdef debugDistance
        cerr << "  At ancestor chain " << chain_index.id_in_parent ;
        cerr << ": " << distL << " " << distR << endl;
    #endif

        } else {
            //If the parent is a snarl
            
            SnarlIndex& snarl_index = snarl_indexes[parent.first];
            size_t node_rank = at_node 
                              ? getPrimaryRank(node_id) 
                              : getSecondaryRank(node_id);
            if (node_rev) {
                //Get the rank of this node in reverse
                node_rank = node_rank % 2 == 1 ? node_rank - 1 : node_rank + 1;
            }
            tie(distL, distR) = snarl_index.distToEnds(node_rank, distL, distR);
            start_len = snarl_index.nodeLength(0);
            end_len = snarl_index.nodeLength(snarl_index.num_nodes*2 - 1);

            node_id = snarl_index.id_in_parent;
            node_rev = snarl_index.rev_in_parent;
            if (snarl_index.in_chain) {
                parent = make_pair(getChainAssignment(snarl_index.parent_id), 
                                    true);
            } else {
                parent = make_pair(getPrimaryAssignment(snarl_index.parent_id), false);
            }

    #ifdef debugDistance
        cerr << "  At ancestor snarl " << snarl_index.id_in_parent ;
        cerr << ": " << distL << " " << distR << endl;
    #endif
        }
        at_node = false;
    }

    return make_tuple (distL, distR, make_pair(node_id, node_rev));
};

int64_t MinimumDistanceIndex::minPos (vector<int64_t> vals) {
    /*return the minimum value in vals that is not -1, returns -1 if all
     values are -1 */
    return accumulate(vals.begin(), vals.end(), -1, 
          [](int x, int y) {if (x==-1) {return y;} 
                            else if (y == -1) {return x;}
                            else {return min(x, y);}} 
          ); 
   
};



size_t MinimumDistanceIndex::getPrimaryAssignment(id_t i) {
    return primary_snarl_assignments[i-min_node_id]- 1;
}
size_t MinimumDistanceIndex::getPrimaryRank(id_t i) {
    return primary_snarl_ranks[i-min_node_id] - 1;
}
size_t MinimumDistanceIndex::getChainAssignment(id_t i) {
    return chain_assignments[has_chain.rank(i-min_node_id)] - 1;
}
size_t MinimumDistanceIndex::getChainRank(id_t i) {
    return chain_ranks[has_chain.rank(i-min_node_id)] - 1;
}
size_t MinimumDistanceIndex::getSecondaryAssignment(id_t i) {
    return secondary_snarl_assignments[has_secondary_snarl.rank(i-min_node_id)]
            - 1;
}
size_t MinimumDistanceIndex::getSecondaryRank(id_t i) {
    return secondary_snarl_ranks[has_secondary_snarl.rank(i-min_node_id)] - 1;
}

void MinimumDistanceIndex::printSelf() {
    //TODO: DOn't actually know the node ids when we're printing things out
    cerr << "node id \t primary snarl \t rank \t secondary snarl \t rank \t chain \t rank" << endl;
    for (size_t i = 0 ; i < primary_snarl_assignments.size() ; i ++ ) {
        if (primary_snarl_assignments[i] != 0){
            cerr << i + min_node_id << "\t";

             cerr << snarl_indexes[primary_snarl_assignments[i]-1].id_in_parent 
                  << "\t" << primary_snarl_ranks[i]-1 << "\t";

            if (has_secondary_snarl_bv[i] == 0) {
                cerr << "/\t/\t";
            } else {
                cerr << snarl_indexes[secondary_snarl_assignments[
                                     has_secondary_snarl.rank(i)]-1].id_in_parent 
                     << "\t" << secondary_snarl_ranks[has_secondary_snarl.rank(i)]-1 << "\t";
            }
            if (has_chain_bv[i] == 0) {
                cerr << "/\t/\t";
            } else {
                cerr << chain_indexes[chain_assignments[has_chain.rank(i)]-1].id_in_parent 
                     << "\t" << chain_ranks[has_chain.rank(i)]-1 << "\t";
            }
            cerr << endl;
        }
    }

    cerr << endl << "Snarls: " << endl;
    for (auto snarls : snarl_indexes) {
        snarls.printSelf();
    }
    cerr << endl << "Chains:" << endl;
    for (auto chains : chain_indexes) {
        chains.printSelf();
    }
}


MinimumDistanceIndex::SnarlIndex::SnarlIndex(
                       id_t parent_id, bool rev_in_parent, id_t id_in_parent, 
                       bool is_unary_snarl, size_t depth, size_t num_nodes, 
                       bool in_chain) :
                       parent_id(parent_id), id_in_parent(id_in_parent), 
                       rev_in_parent(rev_in_parent), 
                       is_unary_snarl(is_unary_snarl), in_chain(in_chain), 
                       depth(depth), num_nodes(num_nodes) {
    /*Constructor for SnarlIndex object that stores distances between
        nodes in a snarl */
    size_t size = num_nodes * 2;
    util::assign(distances, int_vector<>((((size+1)*size)/2) + (size/2), 0));
}

MinimumDistanceIndex::SnarlIndex::SnarlIndex()  {
}
void MinimumDistanceIndex::SnarlIndex::load(istream& in){
    /*Load contents of SnarlIndex from serialization */
    
    distances.load(in);

    sdsl::read_member(in_chain, in);
    sdsl::read_member(parent_id, in);
    sdsl::read_member(rev_in_parent, in);
    sdsl::read_member(id_in_parent, in);
    sdsl::read_member(num_nodes, in);
    sdsl::read_member(depth, in);
    sdsl::read_member(is_unary_snarl, in);
}

void MinimumDistanceIndex::SnarlIndex::serialize(ostream& out) const {
    /* Serialize object to out stream
      Vector contains a header of four ints: #nodes, start node, end node, parent
                  a vector representing visitToIndex [node1, node2, ...] where                          the nodes are ordered by the index they map to
                  a vector representing distances*/

    distances.serialize(out);

    sdsl::write_member(in_chain, out);
    sdsl::write_member(parent_id, out);
    sdsl::write_member(rev_in_parent, out);
    sdsl::write_member(id_in_parent, out);
    sdsl::write_member(num_nodes, out);
    sdsl::write_member(depth, out);
    sdsl::write_member(is_unary_snarl, out);

}


size_t MinimumDistanceIndex::SnarlIndex::index(size_t start, size_t end) {
    /*Get the index of dist from start to end in a snarl distance matrix
      given the node ids + direction */
      
    //The second node must be reversed so that the distance matrix is
    // symmetrical. Since the distance from n1 fd to n2 fd is the same as
    // n2 rev to n1 rev, only one of these is stored
    //
    // Ranks of nodes in the snarl are stored in dist_index as 1fd, 1rev, 2fd...
    // but in the matrix, they are 1fd, 2fd, ... 1rev, 2rev... so that half of 
    // the matrix can be discarded.
    bool rev1 = start % 2 == 1;
    bool rev2 = end % 2 == 1;
    size_t i1 = rev1 ? start/2 + num_nodes : start/2;
    size_t i2 = rev2 ? end/2 : end/2 + num_nodes;
    if (i1 > i2) {
        //Reverse order of nodes
        size_t tmp = i1;
        i1 = i2;
        i2 = tmp;
    }
    
    size_t length = num_nodes * 2;
    size_t k = length - i1;
    return ( ((length + 1) * length ) / 2 ) - ( ((k + 1) * k ) / 2 ) + i2-i1 +
             (length/2);
}

void MinimumDistanceIndex::SnarlIndex::insertDistance(size_t start, 
                                          size_t end, int64_t dist) {
    //Assign distance between start and end
    size_t i = index(start, end);

    distances[i] = dist + 1;
}
   
int64_t MinimumDistanceIndex::SnarlIndex::snarlDistance(size_t start, size_t end) {
    /*Distance between beginnings of two nodes start and end in snarl
     * given their rank
    */
    size_t i = index(start, end);
    return int64_t(distances[i])-1;
}
int64_t MinimumDistanceIndex::SnarlIndex::nodeLength(size_t i){

    return distances[i/2] - 1;
}

int64_t MinimumDistanceIndex::SnarlIndex::snarlLength() {
    //Return the length of the snarl- dist from beginning of start to end of end
    int64_t dist = snarlDistance(0, num_nodes * 2 - 2);
    
     //length of snarl
    if (dist == -1) {
        return -1;
    } else {
        if (num_nodes == 1) {
            return distances[0]-1;
        } else {
            int64_t node_len = nodeLength(num_nodes * 2 - 1) + nodeLength(0);
            return dist + node_len; 
        }
    }
 
}

pair<int64_t, int64_t> MinimumDistanceIndex::SnarlIndex::distToEnds(size_t rank, 
                                                 int64_t distL, int64_t distR) {
    /* Given the distances to either end of a node, find the distances to 
       either end of the snarl
       Rev is true if the node is reversed in the snarl
    */
    int64_t start_len = nodeLength(0);
    int64_t end_len = is_unary_snarl ? start_len 
                                     : nodeLength(num_nodes * 2 - 1);
    size_t rev_rank = rank % 2 == 0 ? rank + 1 : rank - 1;
    size_t end_rank = is_unary_snarl ? 0 : num_nodes * 2 - 1;
    
    int64_t dsl = snarlDistance(0, rank); 

    int64_t dsr = snarlDistance(0, rev_rank);

    int64_t der = snarlDistance(end_rank, rev_rank);

    int64_t del = snarlDistance(end_rank, rank);

    dsl = dsl == -1 ? -1 : dsl + start_len;
    dsr = dsr == -1 ? -1 : dsr + start_len;
    der = der == -1 ? -1 : der + end_len;
    del = del == -1 ? -1 : del + end_len;

    //If the current node is already the start or end position of the snarl
    //then there may be no path between them in the index but the distance is 0
    if (rank == 0) {
        dsl = 0;
    } else if (rev_rank == 0) {
        dsr = 0;
    }
    if (rank == end_rank) {
        del = 0;
    } else if (rev_rank == end_rank) {
        der = 0;
    }
 
    dsl = dsl == -1 || distL == -1? -1 : distL + dsl;
    dsr = dsr == -1 || distR == -1? -1 : distR + dsr;
    der = der == -1 || distR == -1? -1 : distR + der;
    del = del == -1 || distL == -1? -1 : distL + del;

    int64_t dist_start = minPos({dsr, dsl});

    int64_t dist_end = minPos({der, del});

    return make_pair(dist_start, dist_end);
}

void MinimumDistanceIndex::SnarlIndex::printSelf() {
    //Print the nodes contained in SnarlDistance
    cerr << endl;
    if (is_unary_snarl) {
        cerr << "Unary snarl starting at " << id_in_parent;
    } else {
        cerr << "Snarl starting at " << id_in_parent;
    }

    if (in_chain) {
        cerr << endl << "Parent chain: " << parent_id;
    } else {
        cerr << endl << "Parent snarl: " << parent_id;
    }
    cerr << endl << "Length of snarl : " << snarlLength() << endl;

    cerr << "Node lengths; " << endl;
    for (size_t n = 0 ; n < num_nodes; n++) {
        cerr << distances[n]-1 << "\t";
    }
    cerr << endl;
    cerr << "Distances:" << endl;
    cerr << "\t";
    for (size_t n = 0 ; n < num_nodes * 2; n++) {
        cerr << n << "\t";
    }
    cerr << endl;
    for (size_t n1 = 0 ; n1 < num_nodes* 2 ; n1++) {
        cerr << n1 << "\t";
        for (size_t n2 = 0 ; n2 < num_nodes * 2 ; n2++) {
            cerr << snarlDistance(n1, n2) << "\t"; 
        }
        cerr << endl;
    }
    cerr << endl; 
}

//ChainDistance methods

MinimumDistanceIndex::ChainIndex::ChainIndex( 
                       size_t parent_id, size_t id_in_parent, 
                       bool rev_in_parent, bool loops, size_t length):
                parent_id(parent_id), id_in_parent(id_in_parent), 
                is_looping_chain(loops), rev_in_parent(rev_in_parent) {
    

    util::assign(prefix_sum, int_vector<>(length+2, 0));
    util::assign(loop_fd, int_vector<>(length+1, 0));
    util::assign(loop_rev, int_vector<>(length+1, 0));

}
MinimumDistanceIndex::ChainIndex::ChainIndex()  {
}
void MinimumDistanceIndex::ChainIndex::load(istream& in){
    //Populate object from serialization 
    //
    prefix_sum.load(in);
    loop_fd.load(in);
    loop_rev.load(in);

    sdsl::read_member(parent_id, in);
    sdsl::read_member(rev_in_parent, in);
    sdsl::read_member(id_in_parent, in);
    sdsl::read_member(is_looping_chain, in);
}

void MinimumDistanceIndex::ChainIndex::serialize(ostream& out) const {
    /* Serialize the chain index to a file
     * Store startID + endID + parent + 
     * prefix_sum + loop_fd + loop_rev + 
     * snarlToIndex as an int_vector of ids in order of traversal
     */
    prefix_sum.serialize(out);
    loop_fd.serialize(out);
    loop_rev.serialize(out);

    sdsl::write_member(parent_id, out);
    sdsl::write_member(rev_in_parent, out);
    sdsl::write_member(id_in_parent, out);
    sdsl::write_member(is_looping_chain, out);
   

}
int64_t MinimumDistanceIndex::ChainIndex::loopDistance(
         pair<size_t, bool> start, pair<size_t, bool> end, 
         int64_t start_len, int64_t end_len) {

    if (start.first == 0 ) {
        return chainDistance(make_pair(prefix_sum.size() - 2, start.second), 
                             end, start_len,
                             end_len, true);
    } else if (end.first == 0) {
        return chainDistance(start,make_pair(prefix_sum.size() - 2, end.second),
                     start_len, end_len, true);

    } else if (start.first < end.first && start.second) {

        return chainDistance(start, make_pair(0, start.second),
                             start_len, prefix_sum[0]-1, true)
              + 
                chainDistance(make_pair(prefix_sum.size() - 2, start.second),
                             end, prefix_sum[prefix_sum.size() - 1]-1, 
                             end_len, true);
    } else if (start.first > end.first && !start.second) {
        return chainDistance(start, 
                          make_pair(prefix_sum.size() - 2, start.second),
                          start_len, prefix_sum[prefix_sum.size() - 1]-1, true) 
              + 
                chainDistance(make_pair(0, start.second),
                             end, prefix_sum[0]-1, end_len, true);
    } else {
        return -1;
    }

}
int64_t MinimumDistanceIndex::ChainIndex::chainDistance(
         pair<size_t, bool> start, pair<size_t, bool> end, 
         int64_t start_len, int64_t end_len, bool check_loop) {

    /*
     * Return the distance between the given node sides, except node side is
     * specified relative to the reading orientation of the chain that the
     * nodes are in. 
     */

    int64_t start_sum = start.first == 0 ? 0 : prefix_sum[start.first] - 1;
    int64_t end_sum = end.first == 0 ? 0 : prefix_sum[end.first] - 1;
    int64_t loop_dist = -1;
    if (is_looping_chain && !check_loop) {
        loop_dist = loopDistance(start, end, start_len, end_len);
    }
    if (!start.second && !end.second) {
        //If start and end are facing forward relative to the start of the chain
        if (start.first <= end.first) {
            return minPos({loop_dist, end_sum - start_sum});
        } else {
            int64_t rev1 = loop_fd[start.first] - 1;
            int64_t rev2 = loop_rev[end.first] - 1;

            int64_t chain_dist = (start_sum + start_len) - (end_sum + end_len); 
            return minPos({loop_dist, 
                   (rev1 == -1 || rev2 == -1) ? -1 : chain_dist + rev1 + rev2});
        }

    } else if (start.second && end.second ){
        //If start and end are both reversed relative to the start of the chain
        if (start.first >= end.first) {
            return minPos({loop_dist, 
                            (start_sum + start_len) - (end_sum + end_len)});
            
        } else {
            int64_t rev1 = loop_rev[start.first] - 1;
            int64_t rev2 = loop_fd[end.first] - 1;
            int64_t chain_dist = end_sum - start_sum; 
            return minPos({loop_dist, 
                    (rev1 == -1 || rev2 == -1) ? -1 : chain_dist+ rev1 + rev2});
        }
    } else if (!start.second && end.second) {
        //Start is forward, end is reversed
        if (start.first <= end.first) {
            int64_t rev = loop_fd[end.first] - 1;
            int64_t chain_dist = end_sum - start_sum;
            return minPos({loop_dist, rev == -1 ? -1 : rev + chain_dist });
        } else {
            int64_t rev = loop_fd[start.first] - 1;
            int64_t chain_dist = (start_sum+start_len) - (end_sum+end_len);
            return minPos({loop_dist, rev == -1 ? -1 : rev + chain_dist});
        }
        
    } else {
        //start is reverse, end is forward
        if (start.first <= end.first) {
            int64_t rev = loop_rev[start.first] - 1;
            int64_t chain_dist = end_sum - start_sum;
            return minPos({loop_dist, rev == -1 ? -1 : rev + chain_dist});

            
        } else {
            int64_t rev = loop_rev[end.first] - 1; 
            int64_t chain_dist = (start_sum+start_len) - (end_sum+end_len);
            return minPos({loop_dist, rev == -1 ? -1 : rev + chain_dist});
        }
    }
}


int64_t MinimumDistanceIndex::ChainIndex::chainLength() {

    //Get the length of a chain including length of last node
    //TODO: if there is a unary snarl then this should be -1
    return prefix_sum[prefix_sum.size()-1] - 1 ;
}

void MinimumDistanceIndex::ChainIndex::printSelf() {
    //Print the contenst of ChainDistance
   
    if (is_looping_chain) {
        cerr << "Looping chain starting at " << id_in_parent << endl;
    } else {
        cerr << "Chain starting at " << id_in_parent << endl;
    }
    cerr << "Parent snarl " << parent_id << endl;
    
    cerr << "Distances:" << endl;
    cerr << endl;
    for (auto n : prefix_sum) {
        cerr << (int64_t)n - 1 << " ";
    }
    cerr << endl; 
    cerr << "Loop Forward:" << endl;
    cerr << endl;
    for (auto n : loop_fd) {
        cerr << (int64_t)n - 1 << " ";
    }
    cerr << endl; 
    cerr << "Loop Reverse:" << endl;
    cerr << endl;
    for (auto n : loop_rev) {
        cerr << (int64_t)n - 1 << " ";
    }
    cerr << endl;
}


///////////////////  Maximum Distance Index 

void MinimumDistanceIndex::calculateMaxIndex(const HandleGraph* graph, int64_t cap) {

    min_distances.resize(max_node_id - min_node_id + 1);
    max_distances.resize(max_node_id - min_node_id + 1);

    sdsl::util::set_to_value(min_distances, 0);
    sdsl::util::set_to_value(max_distances, 0);


    unordered_map<id_t, pair<id_t, bool>> split_to_id;
    
    //TODO: Unecessary???
    //Make the graph single stranded - each node is traversed in the forward
    //direction

    /* 
    //Make the graph single stranded - each node is traversed in the forward
    //direction
    StrandSplitGraph split_graph (graph);

    unordered_map<id_t, id_t> acyclic_to_id;
    bool is_acyclic = algorithms::is_directed_acyclic(&split_graph);

    //If the graph has directed cycles, dagify it
    VG dagified; 
    HandleGraph* dag;
    if (!is_acyclic) {
#ifdef debugIndex
cerr << "Making graph acyclic" << endl;
#endif
        acyclic_to_id = algorithms::dagify(&split_graph, &dagified, cap);
        dag = &dagified;
    } else {
        dag = &split_graph;
    }

     */
    sglib::HashGraph split;
    split_to_id = algorithms::split_strands(graph, &split);
    graph = &split;


    unordered_map<id_t, id_t> acyclic_to_id;
    bool is_acyclic = algorithms::is_directed_acyclic(&split);

    //If the graph has directed cycles, dagify it
    HandleGraph* dag;
    sglib::HashGraph dagified; 
    if (!is_acyclic) {
#ifdef debugIndex
cerr << "Making graph acyclic" << endl;
#endif
        acyclic_to_id = algorithms::dagify(&split, &dagified, cap);
        dag = &dagified;
    } else {
        dag = &split;
    }

#ifdef debugIndex
    assert(algorithms::is_single_stranded(dag));
    assert(algorithms::is_directed_acyclic(dag));
#endif

    //In DAGified graph, go through graph in topological order and get the 
    //minimum and maximum distances to tips 

    vector<handle_t> order = algorithms::topological_order(dag);

    for (handle_t curr_handle : order) {

        //Get the correct id in the original graph
        id_t curr_id = dag->get_id(curr_handle);
        if (!is_acyclic) {
            curr_id = acyclic_to_id[curr_id];
        }
        curr_id = split_to_id[curr_id].first;

        //Get the previous values for this node
        int64_t curr_min = min_distances[curr_id-min_node_id];
        curr_min = curr_min == 0 ? std::numeric_limits<int64_t>::max() 
                                 : curr_min; 
                              
        int64_t curr_max = max_distances[curr_id-min_node_id];

        //True if this has no incoming nodes 
        bool is_source = true;

        auto check_prev = [&](const handle_t& h)-> bool {
            is_source = false;
            id_t prev_id = dag->get_id(h);
            if (!is_acyclic) {
                prev_id = acyclic_to_id[prev_id];
            }
            prev_id = split_to_id[prev_id].first;
            
            int64_t node_len = dag->get_length(h);
            int64_t prev_min = min_distances[prev_id-min_node_id]+node_len;
            int64_t prev_max = max_distances[prev_id-min_node_id]+node_len;
            curr_min = std::min(curr_min, prev_min);
            curr_max = std::max(curr_max, prev_max);


             return true;
        };
        dag->follow_edges(curr_handle, true, check_prev);
        if (is_source) {
            //If this is a source node, distance is 1 (increment everything by
            //1 so that 0 represents an empty node
            curr_min = 1; 
            curr_max = 1;
        }
        min_distances[curr_id-min_node_id] = curr_min; 
        max_distances[curr_id-min_node_id] = curr_max;
        
    }

}



void MinimumDistanceIndex::printSnarlStats() {
    //Print out stats bout the snarl
   
    cout << "Node count of snarls: " << endl;
    for (auto snarls : snarl_indexes) {
       cout << snarls.num_nodes << "\t"; 
    }
    cout << endl << "Snarl count of chains: " << endl;
    for (auto chains : chain_indexes) {
        cout << chains.loop_fd.size() << "\t";
    }
    cout << endl;
}

}
