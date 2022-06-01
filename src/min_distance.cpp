//#define debugIndex
//#define debugDistance
//#define debugSubgraph

#include "min_distance.hpp"

using namespace std;
namespace vg {

//#define debugIndex

/*TODO: Remove old distance index from vg index 
 * Also change how nodes are stored in chain - in case of loops/unary snarls -might not actually need this
 * Make snarls/chains represented by the node id in netgraph
 */
MinimumDistanceIndex::MinimumDistanceIndex(const HandleGraph* graph, 
                             const SnarlManager* snarl_manager,
                             int64_t cap) {
    /*Constructor for the distance index given a handle graph and snarl manager
    */
    
    #ifdef debugIndex
    assert(graph != nullptr);
    assert(snarl_manager!= nullptr);
    #endif

    min_node_id = graph->min_node_id();
    max_node_id = graph->max_node_id();

    node_to_component.resize(max_node_id - min_node_id + 1);
    sdsl::util::set_to_value(node_to_component, 0);

    component_to_chain_index.resize(24);
    sdsl::util::set_to_value(component_to_chain_index, 0);

    component_to_chain_length.resize(24);
    sdsl::util::set_to_value(component_to_chain_length, 0);

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
    size_t curr_component = 1;//Assign each connected component a unique identifier
    for (const Snarl* snarl : top_snarls) {
       //Make an index for each disconnected snarl/chain
        
        if (seen_snarls.count(snarl) == 0){

            //Assign this connected component to the chain index we're about to make
            if (component_to_chain_index.size() < curr_component) {
                component_to_chain_index.resize(curr_component);
            }
            component_to_chain_index[curr_component-1] = chain_indexes.size();
            int64_t chain_length;

            //Calculate the index for this connected component
            if (snarl_manager->in_nontrivial_chain(snarl)){
                //If this is an actual chain
                const Chain* chain = snarl_manager->chain_of(snarl);

                chain_length = calculate_min_index(graph, snarl_manager, chain, 
                       0, false, false, 0, curr_component);
                for (auto s : *chain) {
                    //TODO: Do this earlier instead of another sweep through the chain
                    seen_snarls.insert(s.first);
                }
            } else {
                //If this is a trivial chain, pretend that its a chain
                Chain curr_chain;
                curr_chain.emplace_back(snarl, false);
                chain_length = calculate_min_index(graph, snarl_manager, &curr_chain,
                             0, false, true, 0, curr_component);
                seen_snarls.insert(snarl);
            }

            //Give this connected component a length
            if (component_to_chain_length.size() < curr_component) {
                component_to_chain_length.resize(curr_component+1);
            }
            component_to_chain_length[curr_component-1] = chain_length;
            curr_component++;

        }
    }

    auto add_single_nodes = [&](const handle_t& h)-> bool {
        id_t id = graph->get_id(h); 
        if (primary_snarl_assignments[id - min_node_id] == 0) {
            //If this node hasn't already been added to the distance index, make a fake snarl for it
            handle_t handle = graph->get_handle(id, false);
            int64_t node_len = graph->get_length(handle);

            //TODO: Do components properly

            //Make a new connected component for this one node
            size_t component_num = component_to_chain_index.size()+1;
            component_to_chain_index.resize(curr_component);
            component_to_chain_length.resize(curr_component);
            //Assign it to a chain that doesn't exist?
            component_to_chain_index[curr_component-1] = chain_indexes.size();
            component_to_chain_length[curr_component-1] = node_len;
            //Also assign this node to a connected component
            node_to_component[id - min_node_id] = component_num;
            

            //Make a snarl index for it
            size_t snarl_assignment = snarl_indexes.size();
            primary_snarl_assignments[id-min_node_id] = snarl_assignment+1;
            primary_snarl_ranks[id - min_node_id] = 1;

            snarl_indexes.emplace_back(0, false, id, id, false, 0, 1, false);
            snarl_indexes.back().distances[0]  = node_len + 1; 
        }
        return true;
               
    };
    graph->for_each_handle(add_single_nodes);

    #ifdef debugIndex
    //Every node should be assigned to a snarl
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
        if ( node_length(id)!=  graph->get_length(handle)){
            cerr << id << ": predicted " <<
                   node_length(id)  << " actual " << graph->get_length(handle);
        }
        assert( node_length(id)  ==  graph->get_length(handle) );
        return true;
               
    };
    graph->for_each_handle(check_assignments);
    #endif


    util::assign(has_secondary_snarl, rank_support_v<1>(&has_secondary_snarl_bv));
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
    util::bit_compress(node_to_component);
    util::bit_compress(component_to_chain_index);
    util::bit_compress(component_to_chain_length);


    if (cap > 0) {
        include_maximum = true;
        calculate_max_index(graph, cap);
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

 
    //Check the file's header to make sure it's the correct version
    if (!in) {
        throw runtime_error("Could not load distance index");
    } else {
        //Check that the header is correct
        size_t char_index = 0;
        //TODO: We're only checking up to the last two so if the header changes this needs to change too
        while (in.peek() != EOF && char_index < file_header.size()-2) {
            char next = (char) in.get();
            if ( next != file_header[char_index]) {
                throw runtime_error ("Distance index file is outdated");
            }
            char_index ++;
        }
        if (in.peek() == '.') {
            char next = (char) in.get();
            if (next != '.' || (char)in.get() != '2') {
                throw runtime_error ("Distance index file is outdated");
            }
            include_component = true;
        } else {
            cerr << "warning: Loading an out-of-date distance index" << endl;
            include_component = false;
        }
        char_index+=2;
        if (char_index < file_header.size()) {
            throw runtime_error ("Distance index file is outdated");
        }
    }

    size_t num_snarls;
    sdsl::read_member(num_snarls, in);
    snarl_indexes.reserve(num_snarls);
    for (size_t i = 0 ; i < num_snarls ; i++) {
        snarl_indexes.emplace_back(); 
        snarl_indexes.back().load(in, include_component);
    }
    primary_snarl_assignments.load(in); 
    primary_snarl_ranks.load(in);
    secondary_snarl_assignments.load(in);
    secondary_snarl_ranks.load(in);
    has_secondary_snarl_bv.load(in);
    has_secondary_snarl.load(in);
    util::assign(has_secondary_snarl, 
                 rank_support_v<1>(&has_secondary_snarl_bv));

    if (include_component) {
        node_to_component.load(in);
        component_to_chain_index.load(in);
        component_to_chain_length.load(in);
    }
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

    //Write the header to the serialized file
    out << file_header;
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
    node_to_component.serialize(out);
    component_to_chain_index.serialize(out);
    component_to_chain_length.serialize(out);

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


//TODO: Add seen nodes here instead of going through them earlier
int64_t MinimumDistanceIndex::calculate_min_index(const HandleGraph* graph,
                                    const SnarlManager* snarl_manager,
                                    const Chain* chain, size_t parent_id,
                                    bool rev_in_parent, bool trivial_chain, 
                                    size_t depth, size_t component_num) {
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


 
    if (!trivial_chain) {
        //If this is a chain, initialize a new ChainIndex object

        //Get the start of the chain
        auto first_visit = get_start_of(*chain);
        chain_indexes.emplace_back(parent_id, first_visit.node_id(), get_end_of(*chain).node_id(),
                                   rev_in_parent,first_visit.node_id()  == get_end_of(*chain).node_id(),
                                   chain->size());

        chain_assignments[first_visit.node_id()-min_node_id] = chain_indexes.size();
        chain_ranks[first_visit.node_id()-min_node_id] = 1;
        has_chain_bv[first_visit.node_id()-min_node_id] = 1; 

        handle_t first_node = graph->get_handle(first_visit.node_id(), first_visit.backward());
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

        if (!trivial_chain && chain_assignments[second_id - min_node_id] == 0){
            //Store the index of the start of the snarl only if it hasn't
            //already been seen (if the chain loops)
            chain_assignments[second_id-min_node_id] = curr_chain_assignment+1;
            chain_ranks[second_id - min_node_id] = curr_chain_rank + 2;
            has_chain_bv[snarl_end_id - min_node_id] = 1;
           
        } 

        NetGraph ng = NetGraph(snarl->start(), snarl->end(), snarl_manager->chains_of(snarl), graph);

        //Get all the nodes in the snarl
        hash_set<pair<id_t, bool>> all_nodes;

        size_t snarl_assignment = snarl_indexes.size();
        auto add_node = [&](const handle_t& h)-> bool {
            id_t id = ng.get_id(h); 
            if (id != snarl_start_id && id != snarl_end_id) {
                const Snarl* temp_snarl = snarl_manager->into_which_snarl(id, false);
                const Snarl* curr_snarl = temp_snarl == NULL ? snarl_manager->into_which_snarl(id, true) : temp_snarl; 
                if (curr_snarl != NULL) {
                    //If this node represents a snarl or chain, then this snarl
                    //is a secondary snarl
                    has_secondary_snarl_bv[id-min_node_id] = 1;
                    secondary_snarl_assignments[id - min_node_id] = snarl_assignment+1;
                    secondary_snarl_ranks[id - min_node_id] = all_nodes.size()+1;
                } else {
                    //Otherwise this is the node's primary snarl
                    primary_snarl_assignments[id-min_node_id] = snarl_assignment+1;
                    primary_snarl_ranks[id - min_node_id] = all_nodes.size()+1;
                    //Also assign this node to a connected component
                    node_to_component[id - min_node_id] = component_num;

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
        primary_snarl_assignments[start_in_chain-min_node_id] = snarl_assignment + 1;
        node_to_component[start_in_chain - min_node_id] = component_num;
        if (start_in_chain == snarl_start_id) {
            primary_snarl_ranks[start_in_chain-min_node_id] =  snarl_start_rev ? 2 : 1;
        } else {
            primary_snarl_ranks[start_in_chain-min_node_id] = snarl_end_rev ? 
                                    all_nodes.size() : all_nodes.size() - 1;
        }
        if (primary_snarl_assignments[end_in_chain-min_node_id] == 0 ){
            //If the 2nd boundary node doesn't already have a primary snarl,
            //then assign it to this snarl
            primary_snarl_assignments[end_in_chain-min_node_id] =  snarl_assignment+1;
            node_to_component[end_in_chain-min_node_id] = component_num;
            primary_snarl_ranks[end_in_chain-min_node_id] =  end_in_chain == snarl_end_id ?  
                 (snarl_end_rev ? all_nodes.size() : all_nodes.size() - 1) :
                 (snarl_start_rev ? 2 : 1);
        } 
        if (!trivial_chain &&
             secondary_snarl_assignments[end_in_chain-min_node_id] == 0){
            //Otherwise, assign the first boundary node a secondary snarl
            secondary_snarl_assignments[end_in_chain-min_node_id] = snarl_assignment+1;
            secondary_snarl_ranks[end_in_chain-min_node_id] = end_in_chain == snarl_end_id ? 
                 (snarl_end_rev ? all_nodes.size()  : all_nodes.size() - 1) :
                 (snarl_start_rev ? 1 : 0);
            has_secondary_snarl_bv[end_in_chain-min_node_id] = 1;
        }

        //Make the snarl index
        if (trivial_chain) {
            //The parent is the parent snarl
            snarl_indexes.emplace_back(parent_id, rev_in_parent, 
                           snarl_start_id, snarl_end_id, snarl_start_id == snarl_end_id, 
                           depth, all_nodes.size()/2, false);
        } else {
            //The parent is the chain
            snarl_indexes.emplace_back(get_start_of(*chain).node_id(), 
                               snarl_rev_in_chain, start_in_chain, end_in_chain, 
                               snarl_start_id == snarl_end_id,
                               depth, all_nodes.size()/2, true);
        }
        populate_snarl_index(graph, snarl_manager, ng, snarl, snarl_rev_in_chain, snarl_assignment, all_nodes, depth, component_num);
#ifdef debugIndex

    snarl_indexes[snarl_assignment].print_self();
    cerr << snarl_indexes[snarl_assignment].max_width << " " << snarl_indexes[snarl_assignment].snarl_length() << endl; 
    assert(snarl_indexes[snarl_assignment].max_width >= snarl_indexes[snarl_assignment].snarl_length());
    cerr << "End snarl " << snarl_indexes[snarl_assignment].id_in_parent << endl;
#endif

        if (!trivial_chain) {
            // Add to prefix sum the distance to the beginning and end of the 
            // last node in the current snarl
            const SnarlIndex& sd = snarl_indexes[snarl_assignment];

            size_t num_nodes = snarl_indexes[snarl_assignment].num_nodes;
            int64_t dist = snarl_rev_in_chain
                ? sd.snarl_distance(num_nodes * 2 - 1, 1) + sd.node_length(num_nodes * 2 - 1)
                : sd.snarl_distance( 0, num_nodes * 2 - 2) + sd.node_length(0);
            #ifdef debugIndex
                cerr << "Prefix sum before snarl: " 
                << chain_indexes[curr_chain_assignment].prefix_sum[curr_chain_rank + 1] << endl;
            #endif
            chain_indexes[curr_chain_assignment].prefix_sum[curr_chain_rank+1] =
                               curr_chain_rank == 0 ? dist + 1 : chain_indexes[curr_chain_assignment].prefix_sum[curr_chain_rank]+dist;


            //Add the reverse loop distance
            if ( curr_chain_rank == 0) {
               //If this is the first snarl in the chain, get the loop distance
               //of the first node
               int64_t first_rev_dist;
               if (snarl_rev_in_chain){ 
                    first_rev_dist = sd.snarl_distance( sd.num_nodes * 2 - 2, sd.num_nodes * 2 - 1);
                    first_rev_dist = first_rev_dist == -1 ? -1 : first_rev_dist + sd.node_length(sd.num_nodes * 2 - 2);
                } else {
                    first_rev_dist = sd.snarl_distance( 1, 0);
                    first_rev_dist = first_rev_dist == -1 ? -1 : first_rev_dist + sd.node_length(0);
                }
                chain_indexes[curr_chain_assignment].loop_rev[0] = first_rev_dist + 1;
            }

            int64_t rev_loop_dist;
            if ( snarl_rev_in_chain ) {
    
                rev_loop_dist = sd.snarl_distance(0, 1);
                rev_loop_dist = rev_loop_dist == -1 ? -1 : rev_loop_dist + sd.node_length(0);
            } else {
                rev_loop_dist = sd.snarl_distance(sd.num_nodes * 2 - 1, sd.num_nodes * 2 - 2);
                rev_loop_dist = rev_loop_dist == -1 ? -1 : rev_loop_dist + sd.node_length(sd.num_nodes * 2 - 2);
            }
     
    
            //Loop distance of the previous node
            int64_t last_loop = chain_indexes[curr_chain_assignment].loop_rev[curr_chain_rank] - 1;

            if (last_loop == -1) {
                chain_indexes[curr_chain_assignment].loop_rev[curr_chain_rank+1] = rev_loop_dist + 1;
            } else {
    
                //Push the minimum of the loop distance of the current snarl and
                //the loop distance of the previous snarl + dist to and from loop 
                int64_t dist_to_end = sd.snarl_distance(0, sd.num_nodes * 2 - 2);
                dist_to_end = dist_to_end == -1 ? -1  
                              : dist_to_end + dist_to_end + sd.node_length(0) + sd.node_length(sd.num_nodes * 2 - 1);


                int64_t loop_distance = min_pos(rev_loop_dist, last_loop + dist_to_end);
               chain_indexes[curr_chain_assignment].loop_rev[curr_chain_rank+1] = loop_distance + 1;
            }
            if ( c == chain_begin(*chain)) {
                //If this is the first snarl, include the length of the start node
                chain_indexes[curr_chain_assignment].max_width += sd.max_width; 
            } else {
                //Otherwise don't
                int64_t start_len = snarl_rev_in_chain
                    ? sd.node_length(num_nodes * 2 - 1) : sd.node_length(0);
                chain_indexes[curr_chain_assignment].max_width += sd.max_width - start_len; 
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
        cd.prefix_sum[cd.prefix_sum.size() - 1] =  cd.prefix_sum[cd.prefix_sum.size() - 2] + graph->get_length(last_node);
    
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
                  snarl_indexes[get_primary_assignment(snarl->start().node_id())] :
                  snarl_indexes[get_primary_assignment(snarl->end().node_id())];
                int64_t new_loop;
                if (curr_chain_rank == 0) {
                    new_loop = cd.loop_rev[cd.loop_rev.size() - 1] - 1;
                } else {
                    int64_t prev_loop = cd.loop_rev[curr_chain_rank - 1] - 1;
                    int64_t dist_to_end = sd.snarl_distance(0, sd.num_nodes * 2 - 2);
                    dist_to_end = dist_to_end == -1 ? -1 : dist_to_end + dist_to_end + sd.node_length(0)  + sd.node_length(sd.num_nodes * 2 - 1);
    
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
            NetGraph ng (snarl->start(), snarl->end(), snarl_manager->chains_of(snarl), graph);
    
                                          
    
            if (c == chain_rbegin(*chain)) {
                //If this is the last snarl in the chain, push loop for last node
    
                int64_t loop_dist_last; 
                if (snarl_rev_in_chain) {
           
                    loop_dist_last = sd.snarl_distance( 1, 0 );
                    loop_dist_last = loop_dist_last == -1 ? -1 : loop_dist_last + sd.node_length(0);
                } else {
    
                    loop_dist_last = sd.snarl_distance(sd.num_nodes * 2 - 2, sd.num_nodes * 2 - 1);
                    loop_dist_last = loop_dist_last == -1 ? -1 : 
                           loop_dist_last + sd.node_length(sd.num_nodes * 2 - 2);
                }
    
                if (get_start_of(*chain).node_id()  == get_end_of(*chain).node_id()) {
                    //If the chain loops, might need distance from first snarl

                    ChainIterator chain_start = chain_begin(*chain);
                    const Snarl* first_snarl = chain_start->first;
                    bool first_snarl_rev = chain_start->second;
                        
                    const SnarlIndex& sd_first = snarl_rev_in_chain ? 
                           snarl_indexes[primary_snarl_assignments[snarl_start_id-min_node_id]-1] :
                           snarl_indexes[primary_snarl_assignments[snarl_end_id-min_node_id]-1];
                    if (first_snarl_rev) {
                        int64_t new_dist = sd_first.snarl_distance(sd_first.num_nodes * 2 - 1, sd_first.num_nodes * 2 - 2);
                        new_dist = new_dist == -1 ? -1 :  new_dist + sd_first.node_length(sd_first.num_nodes * 2 - 2);
                        loop_dist_last = min_pos(loop_dist_last, new_dist);
                    } else {
                        int64_t new_dist = sd_first.snarl_distance(0, 1);
                        new_dist = new_dist == -1 ? -1 : new_dist + sd_first.node_length(0);
                        loop_dist_last = min_pos(loop_dist_last, new_dist);
   
                    }
                  
                }
                cd.loop_fd[curr_chain_rank] = loop_dist_last + 1;
            }
    
            int64_t fd_loop_dist;
    
    
            if (snarl_rev_in_chain) {
                //If the snarl is reversed in the chain
                fd_loop_dist = sd.snarl_distance(sd.num_nodes * 2 - 1,  sd.num_nodes * 2 - 2);
                fd_loop_dist = fd_loop_dist == -1 ? -1 : fd_loop_dist + sd.node_length(sd.num_nodes*2-1);
            } else {
                fd_loop_dist = sd.snarl_distance(0, 1);
                fd_loop_dist = fd_loop_dist == -1 ? -1 : fd_loop_dist + sd.node_length(0);
            }
    
            int64_t last_loop = cd.loop_fd[curr_chain_rank] - 1;
            curr_chain_rank--;
    
            if (last_loop == -1) {
    
                cd.loop_fd[curr_chain_rank] = fd_loop_dist + 1;
    
            } else {
            //push dist to end of snarl + loop dist + dist to start of snarl 
    
                int64_t dist_end_start = sd.snarl_distance( sd.num_nodes * 2 - 1, 1);
                dist_end_start = dist_end_start == -1 ? -1 : dist_end_start + sd.node_length(sd.num_nodes * 2 - 1);


                int64_t dist_start_end =  sd.snarl_distance(0, sd.num_nodes * 2 - 2);
                dist_start_end = dist_start_end == -1 ? -1 : dist_start_end + sd.node_length(0);
                int64_t loop_distance = min_pos(fd_loop_dist, last_loop + dist_end_start + dist_start_end);
                cd.loop_fd[curr_chain_rank] = loop_distance + 1;
            }           
          
        }
        util::bit_compress(cd.prefix_sum);
        util::bit_compress(cd.loop_fd);
        util::bit_compress(cd.loop_rev);
    }
 
    //return length of entire chain
    
    return trivial_chain ? 
       snarl_indexes[primary_snarl_assignments[ get_start_of(*chain).node_id()-min_node_id]-1].snarl_length() :
       chain_indexes[curr_chain_assignment].prefix_sum[chain_indexes[curr_chain_assignment].prefix_sum.size() - 1] - 1;
};

void MinimumDistanceIndex::populate_snarl_index(const HandleGraph* graph, const SnarlManager* snarl_manager, const NetGraph& ng,
                                              const Snarl* snarl, bool snarl_rev_in_chain, size_t snarl_assignment, 
                                              hash_set<pair<id_t, bool>>& all_nodes, size_t depth, size_t component_num) {
    //Fill in the given snarl index's distances by doing a dijkstra search starting from each node in the snarl


    auto cmp = [] (pair<pair<id_t, bool>,int64_t> x,
                                           pair<pair<id_t, bool>,int64_t> y) {
        //Comparison function for the priority of a pair of handle, distance
        return (x.second > y.second);
    };

   id_t snarl_start_id = snarl->start().node_id();
   bool snarl_start_rev = snarl->start().backward(); //into snarl
   id_t snarl_end_id = snarl->end().node_id();
   bool snarl_end_rev = snarl->end().backward();   //pointing out
   id_t start_in_chain = snarl_rev_in_chain ? snarl_end_id : snarl_start_id; 

   //Id of boundary node that occurs second in the chain
   id_t second_id = snarl_rev_in_chain ? snarl_start_id : snarl_end_id;
   for (pair<id_t, bool> start_id : all_nodes){

       //Use each node in the snarl as start of dijkstra search

       //Index of the start node in the current snarl
       size_t start_node_rank = primary_snarl_assignments[start_id.first - min_node_id] - 1 == snarl_assignment
            ? primary_snarl_ranks[start_id.first-min_node_id] - 1
            : secondary_snarl_ranks[start_id.first-min_node_id] - 1;

       if (start_id.second) {
           start_node_rank = start_node_rank % 2 == 0 ? start_node_rank + 1 : start_node_rank - 1;
       }
       handle_t start_handle = graph->get_handle(start_id.first, start_id.second);
       //Priority queue of reachable nodes (pair of node id and direction)
       priority_queue< pair<pair<id_t, bool>, int64_t>, vector<pair<pair<id_t, bool>, int64_t>>, decltype(cmp)> reachable(cmp);
       reachable.emplace(start_id, 0);

#ifdef debugIndex
           cerr << "  Start Node: " << start_id.first << ","  << start_id.second << endl;
           assert( primary_snarl_assignments[start_id.first - min_node_id]-1 == snarl_assignment 
                  || secondary_snarl_assignments[start_id.first - min_node_id]-1 == snarl_assignment);
#endif
        bool first_loop = true;
        unordered_set<pair<id_t, bool>> seen_nodes;

         
        while (reachable.size() > 0) {
            pair<pair<id_t, bool>, int64_t> next = reachable.top();
            reachable.pop();
            pair<id_t, bool> curr_id = next.first;
            handle_t curr_handle = ng.get_handle(curr_id.first, curr_id.second);
            int64_t curr_dist = next.second;

            if (start_id.first != snarl_start_id && start_id.first != snarl_end_id &&
                curr_id.first != snarl_start_id && curr_id.first != snarl_end_id &&
                !first_loop ) {
                //If we started from an interior node and ended at another interior node
                snarl_indexes[snarl_assignment].is_simple_snarl = false;
            }


            if ( seen_nodes.count(curr_id) == 0) {
                //If node has not already been found:

                //Record distance from start to current node 
                if (!first_loop) {
                    size_t curr_node_rank =  primary_snarl_assignments[curr_id.first-min_node_id]-1 == snarl_assignment
                      ? primary_snarl_ranks[curr_id.first-min_node_id] -1 : secondary_snarl_ranks[curr_id.first-min_node_id]-1;

                    if (curr_id.second) {
                        curr_node_rank = curr_node_rank % 2 == 0 ? curr_node_rank + 1 : curr_node_rank - 1;
                    }

                    snarl_indexes[snarl_assignment].insert_distance(start_node_rank, curr_node_rank, curr_dist);
                    seen_nodes.insert(curr_id);

                }

                
                int64_t node_len; //length of the current node
                   
                int64_t loop_dist = -1;
                     //Dist to enter curr node then exit at same side 

                //Get the snarl that the node represents, if any
                const Snarl* temp_snarl = snarl_manager->into_which_snarl(curr_id.first, curr_id.second);
                const Snarl* curr_snarl = temp_snarl == NULL ? 
                      snarl_manager->into_which_snarl(curr_id.first, !curr_id.second) : temp_snarl; 

                if (curr_id.first != snarl_start_id && curr_id.first != snarl_end_id && curr_snarl != NULL) {
                    //If current node is a child snarl/chain

                    snarl_indexes[snarl_assignment].is_simple_snarl = false;


                    if (snarl_manager->in_nontrivial_chain(curr_snarl)) {
                       //The node is a chain

                        const Chain* curr_chain= snarl_manager->chain_of(curr_snarl);
                        size_t chain_start = get_start_of(*curr_chain).node_id();
                        const ChainIndex* chain_dists;
                        if (chain_assignments[chain_start-min_node_id]!= 0){
                            //Length of chain has already been found
                            chain_dists = &chain_indexes[chain_assignments[chain_start-min_node_id]-1];
                           

                        } else {//haven't recursed on this chain yet
                            #ifdef debugIndex
                                cerr << " recurse" << endl;
                            #endif
                            bool rev_in_snarl = curr_id.first == get_start_of(*curr_chain).node_id()
                                      ? get_start_of(*curr_chain).backward() : !get_end_of(*curr_chain).backward();
                            node_len = calculate_min_index(graph,  snarl_manager, curr_chain, 
                                         start_in_chain, rev_in_snarl, false, depth + 1, component_num);

                             chain_dists = &chain_indexes[chain_assignments[chain_start-min_node_id]-1]; 
                        }
                        //Get the length of the node (chain)
                        node_len = chain_dists->chain_length();
                        if (get_start_of( *curr_chain).backward() == curr_id.second) {
                            //If traversing snarl forward in chain

                            loop_dist = chain_dists->loop_fd[0] - 1;

                            if (loop_dist != -1) {
                                auto visit = get_start_of(*curr_chain);
                                handle_t temp_handle= graph->get_handle(visit.node_id(), visit.backward());

                                loop_dist = loop_dist + graph->get_length(temp_handle) ;
                            }
                        } else {

                            loop_dist = chain_dists->loop_rev[chain_dists->loop_rev.size()-1] - 1;

                            if (loop_dist != -1) {
                                auto end_visit = get_end_of(*curr_chain);
                                handle_t temp_handle = graph->get_handle(end_visit.node_id(), end_visit.backward());
                                loop_dist = loop_dist + graph->get_length(temp_handle);
                             }
                        }
                    } else {//Snarl

                        id_t snarl_id = curr_snarl->start().node_id();
                        bool snarl_rev = curr_snarl->start().backward();
                        id_t end_id = curr_snarl->end().node_id();
                        bool end_rev = curr_snarl->end().backward();
  

                        const SnarlIndex* snarl_dists;
                        if (primary_snarl_assignments[snarl_id-min_node_id] != 0) {
                            //Already found
                            snarl_dists = &snarl_indexes[primary_snarl_assignments[snarl_id-min_node_id]-1];
                        } else {//Haven't recursed on snarl yet
                            #ifdef debugIndex
                                cerr << " recurse" << endl;
                            #endif
                            
                            //Create chain to recurse on and recurse
                            Chain curr_chain;

                            curr_chain.emplace_back(curr_snarl, false);
                            bool rev_in_snarl = curr_id.first == snarl_id ? snarl_rev : !end_rev;
                            calculate_min_index(graph, snarl_manager, &curr_chain, start_in_chain,
                                             rev_in_snarl, true, depth + 1, component_num);

                            snarl_dists = &snarl_indexes[primary_snarl_assignments[snarl_id-min_node_id]-1];
                        }
                        node_len = snarl_dists->snarl_length();

                        //Find the distance to enter and exit snarl
                        //at the same side
                        if (curr_id.second == snarl_rev) { 
                            //If traversing snarl forward
                            loop_dist = snarl_dists->snarl_distance(0, 1);

                             if (loop_dist != -1) { 
                                 handle_t temp_handle = graph->get_handle(curr_snarl->start().node_id(), curr_snarl->start().backward());
                                 loop_dist = loop_dist + 2*graph->get_length(temp_handle);
                             }
                        } else {
                            size_t end_in = snarl_dists->is_unary_snarl ? 0 : snarl_dists->num_nodes * 2 - 1;
                            size_t end_out= snarl_dists->is_unary_snarl ? 1 : snarl_dists->num_nodes * 2 - 2;
                            loop_dist = snarl_dists->snarl_distance(end_in, end_out);

                             if (loop_dist != -1) {
                                 handle_t temp_handle = graph->get_handle( curr_snarl->end().node_id(), curr_snarl->end().backward());
                                 loop_dist = loop_dist + 2*graph->get_length(temp_handle);
                             }
                        }
                                    
                    }
                } else { //Node is just a node
                    node_len = graph->get_length(curr_handle);
                }
 
                if (curr_id == start_id) {
                    snarl_indexes[snarl_assignment].distances[start_node_rank/2]  = node_len + 1; 
                }
   

                if (loop_dist != -1 && !first_loop) {
                    /*If there is a path within the current node that loops 
                      to enter the node and exit it at the same side - add
                      reachable nodes from current node in reverse 
                      Do not add this distance if the current node is the 
                      starting node */

                    handle_t rev_handle = ng.get_handle(ng.get_id(curr_handle), !ng.get_is_reverse(curr_handle)); 
                        

                    auto add_rev_handle = [&](const handle_t& h)-> bool {
                        pair<id_t, bool> node = make_pair(ng.get_id(h), ng.get_is_reverse(h));
                        reachable.emplace(node, curr_dist + loop_dist);
 
                         return true;
                    };

                    ng.follow_edges(rev_handle, false, add_rev_handle);
                }

                
                //Add reachable nodes to priority queue
                auto add_handle = [&](const handle_t& h)-> bool {

                    pair<id_t, bool> node = make_pair(ng.get_id(h), ng.get_is_reverse(h));
                    if (node_len != -1) {
                        reachable.emplace(node, curr_dist + node_len);
                    }
                  
                    
                    #ifdef debugIndex
                        cerr << node.first << " " << node.second << ", ";
                    #endif
                    return true;
                };
                //Add reachable nodes to priority queue for unary snarl that doesn't loop - 0 distance
                auto add_unary_handle = [&](const handle_t& h)-> bool {
                    pair<id_t, bool> node = make_pair(ng.get_id(h), ng.get_is_reverse(h));
                    reachable.emplace(node, 0);
                   
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

                    ng.follow_edges(curr_handle, false, add_unary_handle);

                } else  {
                    ng.follow_edges(curr_handle, false, add_handle);
                }  
                    

                //Add edges between the boundary nodes that are not in 
                //the net graph
                int64_t next_dist = curr_id == start_id ? 0 : curr_dist+node_len;

                if ((curr_id.first == snarl_start_id && curr_id.second != snarl_start_rev) ||
                     ( curr_id.first == snarl_end_id && curr_id.second == snarl_end_rev )  ) {
                       
                    //If currently leaving the snarl
                    graph->follow_edges(curr_handle, false, [&](const handle_t& h)-> bool {
                        pair<id_t, bool> node = make_pair(ng.get_id(h), ng.get_is_reverse(h));
                        if ( node.first == snarl_start_id || node.first == snarl_end_id ) {
                           reachable.emplace(node, next_dist);
                        }
                        return true;
                    });

                }                      
#ifdef debugIndex
     cerr << "    prev dist: " << curr_dist << "+ new dist " << node_len << endl;
#endif
            }
            first_loop = false;
        }//End while loop
        if ( snarl_indexes[snarl_assignment].snarl_distance(0, start_node_rank) != -1 && 
            snarl_indexes[snarl_assignment].snarl_distance(snarl_indexes[snarl_assignment].num_nodes * 2 - 1, start_node_rank) != -1) {
            //If there is an edge to both start and end
            snarl_indexes[snarl_assignment].is_simple_snarl = false;
        }
    }//End for loop over starting node/directions in a snarl


    //Find the maximum width of the snarl
    //Max width is the maximum among all minimum distance paths - for each node, the sum of the minimum distance to each end
    SnarlIndex& curr_snarl_index = snarl_indexes[snarl_assignment];

    ng.for_each_handle([&](const handle_t& h)-> bool {

        int64_t node_len;
        size_t node_rank;
        id_t curr_id = ng.get_id(h);
        bool curr_rev = ng.get_is_reverse(h);

        const Snarl* temp_snarl = snarl_manager->into_which_snarl(curr_id, curr_rev);
        const Snarl* curr_snarl = temp_snarl == NULL ? snarl_manager->into_which_snarl(curr_id, !curr_rev) : temp_snarl; 
        if (curr_id != snarl_start_id && curr_id != snarl_end_id && curr_snarl != NULL) {
             //This is a snarl or chain
             if (snarl_manager->in_nontrivial_chain(curr_snarl)) {
                 //The node is a chain
                 size_t chain_start = get_start_of(*snarl_manager->chain_of(curr_snarl)).node_id();
                 const ChainIndex& chain_dists = chain_indexes[chain_assignments[chain_start-min_node_id]-1];
                 node_rank = secondary_snarl_ranks[chain_dists.id_in_parent-min_node_id] - 1;
                 node_len = chain_dists.max_width;

             } else {
                 //This node is a snarl
                 id_t snarl_id = curr_snarl->start().node_id();
                 const SnarlIndex& snarl_dists = snarl_indexes[primary_snarl_assignments[snarl_id-min_node_id]-1];
                 node_rank = secondary_snarl_ranks[snarl_dists.id_in_parent-min_node_id] - 1;
                 node_len = snarl_dists.max_width;
             }
         } else {
             //This is a node
             node_rank = get_primary_rank(curr_id);
             node_len = graph->get_length(h);
         }
         pair<int64_t, int64_t> dists_to_end = curr_snarl_index.dist_to_ends(node_rank, 0, 0);
         int64_t min_max_dist = dists_to_end.first == -1 || dists_to_end.second == -1 ? -1 : dists_to_end.first + dists_to_end.second + node_len;
        
        //Keep track of the maximum of all minimum distances
        curr_snarl_index.max_width = std::max( snarl_indexes[snarl_assignment].max_width,min_max_dist);
        return true;
    });


}

//////////////////    Distance Calculations
//

int64_t MinimumDistanceIndex::max_distance(pos_t pos1, pos_t pos2) const {
    if (!include_maximum) {
        return -1;
    } else {
        id_t id1 = get_id(pos1);
        id_t id2 = get_id(pos2);
        int64_t len1 = node_length(id1);
        int64_t len2 = node_length(id2);

        len1 = std::max((int64_t)(get_offset(pos1)+1), (int64_t)(len1-get_offset(pos1)));
        len2 = std::max((int64_t)(get_offset(pos2)+1), (int64_t)(len2-get_offset(pos2)));

        id1 -= min_node_id;
        id2 -= min_node_id;

        int64_t max_dist = std::max((int) max_distances[id1] - (int)min_distances[id2], 
                                    (int) max_distances[id2] - (int)min_distances[id1]);

        return len1 + len2 + max_dist;

    }
}

tuple<id_t, bool, bool> MinimumDistanceIndex::into_which_snarl(id_t node_id, bool reverse) const {
    size_t primary_assignment = get_primary_assignment(node_id);
    const SnarlIndex& primary_snarl_index = snarl_indexes[primary_assignment];

    //Rank always returns the rank of the forward orientation of the node
    //The first (0) and last (num_nodes*2+1) rank will be pointing in, so:
    //If this is the start node of snarl, 
    // rank is 0 if it is traversed forward to enter the snarl, 1 if it is traversed backward
    //If it is the end node of the primary snarl, rank is even if it is traversed backward to enter the snarl
    pair<id_t, bool> primary_start (primary_snarl_index.id_in_parent, 
                     get_primary_assignment(primary_snarl_index.id_in_parent) == primary_assignment ? 
                                    get_primary_rank(primary_snarl_index.id_in_parent) == 1 :
                                    get_secondary_rank(primary_snarl_index.id_in_parent) == 1 );
    pair<id_t, bool> primary_end (primary_snarl_index.end_id, 
                     get_primary_assignment(primary_snarl_index.end_id) == primary_assignment ? 
                                  get_primary_rank(primary_snarl_index.end_id)%2 == 0 : 
                                  get_secondary_rank(primary_snarl_index.end_id)%2 == 0);
    if (primary_snarl_index.is_unary_snarl) {
        primary_end = primary_start;
    }

    if ((node_id == primary_start.first && reverse == primary_start.second) || 
        (node_id == primary_end.first && reverse == primary_end.second )) {
        //If this is then end node and it points in
        return make_tuple(primary_start.first, primary_start.second, primary_snarl_index.is_trivial_snarl());
    }

    if (has_secondary_snarl_bv[node_id-min_node_id]){ 
        size_t secondary_assignment = get_secondary_assignment(node_id);
        const SnarlIndex& secondary_snarl_index = snarl_indexes[secondary_assignment];
        pair<id_t, bool> secondary_start (secondary_snarl_index.id_in_parent, 
                         get_primary_assignment(secondary_snarl_index.id_in_parent) == secondary_assignment ?
                         get_primary_rank(secondary_snarl_index.id_in_parent) == 1 :
                         get_secondary_rank(secondary_snarl_index.id_in_parent) == 1);
        pair<id_t, bool> secondary_end (secondary_snarl_index.end_id, 
                         get_primary_assignment(secondary_snarl_index.id_in_parent) == secondary_assignment ?
                         get_primary_rank(secondary_snarl_index.end_id)%2 == 0 :
                         get_secondary_rank(secondary_snarl_index.end_id)%2 == 0);
        if (secondary_snarl_index.is_unary_snarl) {
            secondary_end = secondary_start;
        }
                            
        if ( (node_id == secondary_start.first && reverse == secondary_start.second) ||
             (node_id == secondary_end.first && reverse == secondary_end.second)) {
            //If this is the inward facing start or end node of the secondary snarl
            return make_tuple(secondary_start.first, secondary_start.second, secondary_snarl_index.is_trivial_snarl());
        }
    }
    //This does not point into a snarl
    return make_tuple(0, false, false);
}

int64_t MinimumDistanceIndex::node_length(id_t id) const {
    return snarl_indexes[get_primary_assignment(id)].node_length(get_primary_rank(id));
}
int64_t MinimumDistanceIndex::min_distance(pos_t pos1, pos_t pos2) const {
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
    pair<id_t, bool> ancestor1 ( snarl_indexes[get_primary_assignment(node_id1)].id_in_parent, false);

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
            size_t chain_assignment = get_chain_assignment(ancestor1.first);
            ancestors1.emplace(chain_assignment, true );
            ancestor1 = make_pair(chain_indexes[chain_assignment].parent_id, false);

        } else {
            size_t snarl_assignment = get_primary_assignment(ancestor1.first);
            ancestors1.emplace(snarl_assignment, false );
            const SnarlIndex& si = snarl_indexes[snarl_assignment];
            ancestor1 = make_pair(si.parent_id, si.in_chain);
        }
    }

#ifdef debugDistance
      cerr << endl << "ancestors of 2: ";
#endif


    pair<id_t, bool> ancestor2 ( snarl_indexes[get_primary_assignment(node_id2)].id_in_parent,false);
    while (ancestor2.first != 0) {

#ifdef debugDistance
        if (ancestor2.second) {cerr << "chain ";}
        else {cerr << "snarl ";}
        cerr << ancestor2.first << " ";
#endif

        if (ancestor2.second) {
            //If this ancestor is a chain
            size_t chain_assignment = get_chain_assignment(ancestor2.first);
            if (ancestors1.count(make_pair(chain_assignment, true)) != 0) {
                common_ancestor = ancestor2;
                break;
            }
            ancestor2 = make_pair(chain_indexes[chain_assignment].parent_id, false);
        } else { 
            size_t snarl_assignment = get_primary_assignment(ancestor2.first);
            if (ancestors1.count(make_pair(snarl_assignment, false)) != 0) {
                common_ancestor = ancestor2;
                break;
            }
            const SnarlIndex& si = snarl_indexes[snarl_assignment];
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

    pair<size_t, bool> ancestor ( common_ancestor.second  ? get_chain_assignment(common_ancestor.first)
                                                         : get_primary_assignment(common_ancestor.first),
                                  common_ancestor.second);


    //Find distances from pos1 and pos2 to ends of child snarls of ancestor
    int64_t distL1; int64_t distR1; pair<id_t, bool> snarl_tree_node1;
    tie (distL1, distR1, snarl_tree_node1) = dist_to_common_ancestor(ancestor, pos1, false);

     
    int64_t distL2; int64_t distR2; pair<id_t, bool> snarl_tree_node2;
    tie (distL2, distR2, snarl_tree_node2) = dist_to_common_ancestor(ancestor, pos2, true);

    pair<id_t, bool> parent = common_ancestor;
    bool lowest_ancestor = true;
    while (parent.first != 0) {
        //snarl_tree_nodes 1 and 2 are children of parent snarl or chain 
        if (parent.second) {
            //If the parent is a chain and both snarl_tree_nodes are snarls

            const ChainIndex& chain_index = chain_indexes[get_chain_assignment(parent.first)]; 

            //If the two nodes are snarls in the common ancestor chain
            //find the distance between them in the chain
            const SnarlIndex& snarl_index1 = snarl_indexes[get_primary_assignment( snarl_tree_node1.first)];
            size_t start_rank1 = get_chain_rank(snarl_index1.id_in_parent);
            size_t end_rank1 = start_rank1 + 1; 
            if (snarl_index1.rev_in_parent) {
                int64_t temp = distL1;
                distL1 = distR1;
                distR1 = temp;
            }
            int64_t start_len1 = snarl_index1.node_length(snarl_index1.rev_in_parent ? snarl_index1.num_nodes * 2 - 1 : 0);
            int64_t end_len1 = snarl_index1.node_length(snarl_index1.rev_in_parent ? 0 : snarl_index1.num_nodes * 2 - 1);
            size_t start_rank2, end_rank2;
            int64_t start_len2, end_len2;

            if (lowest_ancestor) {
                //If this is the lowest common ancestor, then there are two
                //separate nodes that need to be found
               
               const SnarlIndex& snarl_index2 = snarl_indexes[get_primary_assignment( snarl_tree_node2.first)];
               start_rank2 = get_chain_rank(snarl_index2.id_in_parent);
               end_rank2 = start_rank2+1; 
               if (snarl_index2.rev_in_parent) {
                   int64_t temp = distL2;
                   distL2 = distR2;
                   distR2 = temp;
               }
               start_len2 = snarl_index2.node_length(snarl_index2.rev_in_parent ? snarl_index2.num_nodes * 2 - 1 : 0);
               end_len2 = snarl_index2.node_length(snarl_index2.rev_in_parent ? 0 : snarl_index2.num_nodes * 2 - 1);
            } else {
                //Otherwise just copy the first one
                start_rank2 = start_rank1;
                end_rank2 = end_rank1;
                start_len2 = start_len1;
                end_len2 = end_len1;
            }

            //Distance from left of s1 (reverse), left of s2 (forward)
            int64_t d1 = chain_index.chain_distance(make_pair(start_rank1, true), make_pair(start_rank2, false), start_len1, start_len2);
            d1 = (distL1 == -1 || distL2 == -1 || d1 == -1) ? -1 : distL1 + distL2 + d1 - start_len1;

            //Distance from left of s1 (reverse) to right of s2 (reverse)
            int64_t d2;
            if (start_rank1 == end_rank2) {
                //If snarls share a node, then the distances up to this point
                //both include the length of the shared node
                d2 = (distL1 == -1 || distR2 == -1) ? -1 :  distL1 + distR2 - start_len1; 
            } else {
                d2 = chain_index.chain_distance(make_pair(start_rank1, true),  make_pair(end_rank2, true), start_len1, end_len2);
                d2 = (distL1 == -1 || distR2 == -1 || d2 == -1) ? -1 :  distL1 + distR2 + d2 - start_len1;
            }

            //Distance from right of s1 (fd) to left of s2 (fd)
            int64_t d3;
            if (end_rank1 == start_rank2) {
                d3 = (distR1 == -1 || distL2 == -1) ? -1 :  distR1 + distL2 - end_len1; 
            } else {
                d3 = chain_index.chain_distance(make_pair(end_rank1, false), make_pair(start_rank2, false), end_len1, start_len2);

                d3 = (distR1 == -1 || distL2 == -1 || d3 == -1) ? -1 : distR1 + distL2 + d3 - end_len1; 
            }

            //Distance from right of s1 (fd) to right of s2 (rev)
            int64_t d4 =  chain_index.chain_distance(make_pair(end_rank1, false),  make_pair(end_rank2, true), end_len1, end_len2);
            d4 = (distR1 == -1 || distR2 == -1 || d4 == -1) ? -1 : distR1 + distR2 + d4 - end_len1;
            
                       
            shortest_distance = min_pos({d1, d2, d3, d4, shortest_distance});


            //Extend distances of both positions to the ends of the chain

            //Find the rank of the end node and current node in the chain
            size_t chain_end_rank = chain_index.prefix_sum.size() - 2;
                                    
            //Get the lengths of start, end, and current node
            int64_t chain_start_len = node_length(chain_index.id_in_parent);

            int64_t chain_end_len = chain_index.prefix_sum[chain_index.prefix_sum.size()-1] 
                                - chain_index.prefix_sum[chain_index.prefix_sum.size()-2];

            int64_t dsl = chain_index.chain_distance(make_pair(0, false),  make_pair(start_rank1, false), chain_start_len, start_len1);
            int64_t dsr = chain_index.chain_distance(make_pair(0, false), make_pair(end_rank1, true), chain_start_len, end_len1);
            int64_t der = chain_index.chain_distance( make_pair(chain_end_rank, true), make_pair(end_rank1, true), chain_end_len, end_len1);
            int64_t del = chain_index.chain_distance(make_pair(chain_end_rank, true), make_pair(start_rank1, false), chain_end_len, start_len1);

            int64_t dsl1 = dsl == -1 || distL1 == -1? -1 : distL1 + dsl; 
            int64_t dsr1 =  dsr == -1 || distR1 == -1? -1 : distR1 + dsr; 
            int64_t der1 = der == -1 || distR1 == -1? -1 : distR1 + der; 
            int64_t del1 = del == -1 || distL1 == -1? -1 : distL1 + del; 
 
            distL1 = min_pos(dsr1, dsl1);
            distR1 = min_pos(der1, del1);

            if (lowest_ancestor) {
                //If the two snarl tree nodes are different, need to find 
                //distances to ends for the second one

                dsl = chain_index.chain_distance(make_pair(0, false),  make_pair(start_rank2, false), chain_start_len, start_len2);
                dsr = chain_index.chain_distance(make_pair(0, false), make_pair(end_rank2, true), chain_start_len, end_len2);
                der = chain_index.chain_distance(make_pair(chain_end_rank, true),  make_pair(end_rank2, true), chain_end_len, end_len2);
                del = chain_index.chain_distance(make_pair(chain_end_rank, true), make_pair(start_rank2, false), chain_end_len, start_len2);

            }

            int64_t dsl2 = dsl == -1 || distL2 == -1? -1 : distL2 + dsl; 
            int64_t dsr2 = dsr == -1 || distR2 == -1? -1 : distR2 + dsr; 
            int64_t der2 = der == -1 || distR2 == -1? -1 : distR2 + der; 
            int64_t del2 = del == -1 || distL2 == -1? -1 : distL2 + del; 

            distL2 = min_pos(dsr2, dsl2);
            distR2 = min_pos(der2, del2);

            snarl_tree_node1 = make_pair(chain_index.id_in_parent,chain_index.rev_in_parent);


#        ifdef debugDistance
            cerr << "At ancestor chain " << chain_indexes[
                     get_chain_assignment(parent.first)].id_in_parent << endl;
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

            size_t parent_snarl_index = get_primary_assignment(parent.first);
            const SnarlIndex& snarl_index = snarl_indexes[parent_snarl_index];

            size_t rank1 = get_primary_assignment(snarl_tree_node1.first)  == parent_snarl_index
                           ? get_primary_rank(snarl_tree_node1.first) : get_secondary_rank(snarl_tree_node1.first);
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
                rank2 = get_primary_assignment(snarl_tree_node2.first ) == parent_snarl_index
                         ? get_primary_rank(snarl_tree_node2.first) : get_secondary_rank(snarl_tree_node2.first);
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


            int64_t d1 = snarl_index.snarl_distance(rank1, rank2);
            d1 = (distR1 == -1 || distL2 == -1 || d1 == -1) ? -1 : distR1 + distL2 + d1; 

            int64_t d2 = snarl_index.snarl_distance(rank1, rev_rank2);

            d2 = (distR1 == -1 || distR2 == -1 || d2 == -1) ? -1 : distR1 + distR2 + d2;
            int64_t d3 = snarl_index.snarl_distance(rev_rank1, rank2);
            d3 = (distL1 == -1 || distL2 == -1 || d3 == -1) ? -1 :  distL1 + distL2 + d3; 
            int64_t d4 = snarl_index.snarl_distance(rev_rank1, rev_rank2);
            d4 = (distL1 == -1 || distR2 == -1 || d4 == -1) ? -1 :  distL1 + distR2 + d4; 

            shortest_distance =  min_pos({d1, d2, d3, d4, shortest_distance});

            
            
            //Find distances to the ends of the parent snarl
            tie (distL1, distR1) = snarl_index.dist_to_ends(rank1, distL1, distR1);

            tie (distL2, distR2) = snarl_index.dist_to_ends(rank2, distL2, distR2);
             

            snarl_tree_node1 = make_pair(snarl_index.id_in_parent, snarl_index.rev_in_parent);

#        ifdef debugDistance
            cerr << "At ancestor snarl " << snarl_indexes[get_primary_assignment( parent.first)].id_in_parent << endl;
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

tuple<int64_t, int64_t, pair<id_t, bool>> MinimumDistanceIndex::dist_to_common_ancestor(
          pair<size_t, bool> common_ancestor, pos_t& pos, bool rev) const {

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
    
    size_t snarl_assignment = get_primary_assignment(node_id);
    size_t snarl_rank = get_primary_rank(node_id);
    const SnarlIndex& snarl_index = snarl_indexes[snarl_assignment];

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
        distL = snarl_index.node_length(snarl_rank)-offset;
    } else {
        distL = offset+1;
        distR = snarl_index.node_length(snarl_rank)-offset;
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
            const ChainIndex& chain_index = chain_indexes[parent.first];

            size_t start_rank = get_chain_rank(node_id);
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

            int64_t chain_end_len = chain_index.prefix_sum[chain_index.prefix_sum.size()-1] 
                                        - chain_index.prefix_sum[chain_index.prefix_sum.size()-2];

            int64_t dsl = chain_index.chain_distance(make_pair(0, false), make_pair(start_rank, false),chain_start_len, start_len);
            int64_t dsr = chain_index.chain_distance(make_pair(0, false), make_pair(end_rank, true), chain_start_len, end_len);
            int64_t der =chain_index.chain_distance(make_pair(chain_end_rank, true),make_pair(end_rank, true), chain_end_len, end_len);
            int64_t del =chain_index.chain_distance(make_pair(chain_end_rank, true), make_pair(start_rank, false), chain_end_len, start_len);

            dsl = dsl == -1 || distL == -1? -1 : distL + dsl;
            dsr = dsr == -1 || distR == -1? -1 : distR + dsr;
            der = der == -1 || distR == -1? -1 : distR + der;
            del = del == -1 || distL == -1? -1 : distL + del;
 
            distL = min_pos(dsr, dsl);
            distR = min_pos(der, del);

            node_id = chain_index.id_in_parent;
            node_rev = chain_index.rev_in_parent;
            parent = make_pair(get_primary_assignment(chain_index.parent_id), false);

    #ifdef debugDistance
        cerr << "  At ancestor chain " << chain_index.id_in_parent ;
        cerr << ": " << distL << " " << distR << endl;
    #endif

        } else {
            //If the parent is a snarl
            
            const SnarlIndex& snarl_index = snarl_indexes[parent.first];
            size_t node_rank = at_node ? get_primary_rank(node_id) : get_secondary_rank(node_id);
            if (node_rev) {
                //Get the rank of this node in reverse
                node_rank = node_rank % 2 == 1 ? node_rank - 1 : node_rank + 1;
            }
            tie(distL, distR) = snarl_index.dist_to_ends(node_rank, distL, distR);
            start_len = snarl_index.node_length(0);
            end_len = snarl_index.node_length(snarl_index.num_nodes*2 - 1);

            node_id = snarl_index.id_in_parent;
            node_rev = snarl_index.rev_in_parent;
            if (snarl_index.in_chain) {
                parent = make_pair(get_chain_assignment(snarl_index.parent_id),  true);
            } else {
                parent = make_pair(get_primary_assignment(snarl_index.parent_id), false);
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

int64_t MinimumDistanceIndex::min_pos (vector<int64_t> vals) {
    /*return the minimum value in vals that is not -1, returns -1 if all
     values are -1 */
    return accumulate(vals.begin(), vals.end(), -1, 
          [](int64_t x, int64_t y) {
              return min_pos(x, y);
          });
   
};

void MinimumDistanceIndex::print_self() const {
    cerr << "node id \t primary snarl \t rank \t secondary snarl \t rank \t chain \t rank" << endl;
    for (size_t i = 0 ; i < primary_snarl_assignments.size() ; i ++ ) {
        if (primary_snarl_assignments[i] != 0){
            cerr << i + min_node_id << "\t";

             cerr << snarl_indexes[primary_snarl_assignments[i]-1].id_in_parent  << "\t" << primary_snarl_ranks[i]-1 << "\t";

            if (has_secondary_snarl_bv[i] == 0) {
                cerr << "/\t/\t";
            } else {
                cerr << snarl_indexes[secondary_snarl_assignments[has_secondary_snarl.rank(i)]-1].id_in_parent 
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
        snarls.print_self();
    }
    cerr << endl << "Chains:" << endl;
    for (auto chains : chain_indexes) {
        chains.print_self();
    }
}


MinimumDistanceIndex::SnarlIndex::SnarlIndex(
                       id_t parent_id, bool rev_in_parent, id_t id_in_parent, id_t end_id, 
                       bool is_unary_snarl, size_t depth, size_t num_nodes, 
                       bool in_chain) :
                       parent_id(parent_id), id_in_parent(id_in_parent), end_id(end_id), 
                       rev_in_parent(rev_in_parent), 
                       is_unary_snarl(is_unary_snarl), in_chain(in_chain), 
                       depth(depth), num_nodes(num_nodes), max_width(-1) {
    /*Constructor for SnarlIndex object that stores distances between
        nodes in a snarl */
    size_t size = num_nodes * 2;
    is_simple_snarl = true;
    util::assign(distances, int_vector<>((((size+1)*size)/2) + (size/2), 0));
}

MinimumDistanceIndex::SnarlIndex::SnarlIndex()  {
}
void MinimumDistanceIndex::SnarlIndex::load(istream& in, bool include_component){
    /*Load contents of SnarlIndex from serialization */
    
    distances.load(in);

    sdsl::read_member(in_chain, in);
    sdsl::read_member(parent_id, in);
    sdsl::read_member(rev_in_parent, in);
    sdsl::read_member(id_in_parent, in);
    sdsl::read_member(end_id, in);
    sdsl::read_member(num_nodes, in);
    sdsl::read_member(depth, in);
    sdsl::read_member(is_unary_snarl, in);
    if (include_component) {
        sdsl::read_member(is_simple_snarl, in);
    }
    sdsl::read_member(max_width, in);
}

void MinimumDistanceIndex::SnarlIndex::serialize(ostream& out) const {
    /* Serialize object to out stream
      Vector contains a header of four ints: #nodes, start node, end node, parent
                  a vector representing visitToIndex [node1, node2, ...] where
                  the nodes are ordered by the index they map to
                  a vector representing distances*/

    distances.serialize(out);

    sdsl::write_member(in_chain, out);
    sdsl::write_member(parent_id, out);
    sdsl::write_member(rev_in_parent, out);
    sdsl::write_member(id_in_parent, out);
    sdsl::write_member(end_id, out);
    sdsl::write_member(num_nodes, out);
    sdsl::write_member(depth, out);
    sdsl::write_member(is_unary_snarl, out);
    sdsl::write_member(is_simple_snarl, out);
    sdsl::write_member(max_width, out);

}


size_t MinimumDistanceIndex::SnarlIndex::index(size_t start, size_t end) const {
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
    return ( ((length + 1) * length ) / 2 ) - ( ((k + 1) * k ) / 2 ) + i2-i1 + (length/2);
}

void MinimumDistanceIndex::SnarlIndex::insert_distance(size_t start, 
                                          size_t end, int64_t dist) {
    //Assign distance between start and end
    size_t i = index(start, end);

    distances[i] = dist + 1;
}

int64_t MinimumDistanceIndex::SnarlIndex::snarl_length() const {
    //Return the length of the snarl- dist from beginning of start to end of end
    int64_t dist = is_unary_snarl ? snarl_distance(0, 1) : snarl_distance(0, num_nodes * 2 - 2);
    
     //length of snarl
    if (dist == -1) {
        return -1;
    } else {
        if (num_nodes == 1) {
            return distances[0]-1;
        } else {
            int64_t node_len = is_unary_snarl ? node_length(0) * 2 : node_length(num_nodes * 2 - 1) + node_length(0);
            return dist + node_len; 
        }
    }
 
}

pair<int64_t, int64_t> MinimumDistanceIndex::SnarlIndex::dist_to_ends(size_t rank, 
                                                 int64_t distL, int64_t distR) const {
    /* Given the distances to either end of a node, find the distances to 
       either end of the snarl
       Rev is true if the node is reversed in the snarl
    */
    int64_t start_len = node_length(0);
    int64_t end_len = is_unary_snarl ? start_len : node_length(num_nodes * 2 - 1);
    size_t rev_rank = rank % 2 == 0 ? rank + 1 : rank - 1;
    size_t end_rank = is_unary_snarl ? 0 : num_nodes * 2 - 1;
    
    int64_t dsl = snarl_distance(0, rank); 

    int64_t dsr = snarl_distance(0, rev_rank);

    int64_t der = snarl_distance(end_rank, rev_rank);

    int64_t del = snarl_distance(end_rank, rank);

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
 
    dsl = dsl == -1 || distL == -1 ? -1 : distL + dsl;
    dsr = dsr == -1 || distR == -1 ? -1 : distR + dsr;
    der = der == -1 || distR == -1 ? -1 : distR + der;
    del = del == -1 || distL == -1 ? -1 : distL + del;

    int64_t dist_start = min_pos(dsr, dsl);

    int64_t dist_end = min_pos(der, del);

    return make_pair(dist_start, dist_end);
}

bool MinimumDistanceIndex::SnarlIndex::is_trivial_snarl() const {
    return num_nodes == 2 && !is_unary_snarl && snarl_distance(0,1) == -1 && snarl_distance(4,3) == -1; 
}

json_t*  MinimumDistanceIndex::SnarlIndex::snarl_to_json() {
    json_t* out_json = json_object();
    json_object_set_new(out_json, "type", json_string(is_unary_snarl ? "unary snarl" : "snarl"));

    json_t* start_node = json_object();
    json_object_set_new(start_node, "node_id", json_integer(id_in_parent));
    json_t* end_node = json_object();
    json_object_set_new(end_node, "node_id", json_integer(end_id));
    json_t* parent = json_object();
    if (parent_id == 0 && depth == 0 && !in_chain ) {
        json_object_set_new(parent, "type", json_string("root"));
    } else {
        json_object_set_new(parent, "type", json_string(in_chain ? "chain" : "snarl"));
        json_object_set_new(parent, "node_id", json_integer(parent_id));
    }

    json_object_set_new(out_json, "start", start_node);
    json_object_set_new(out_json, "end", end_node);
    json_object_set_new(out_json, "parent", parent);

    json_object_set_new(out_json, "node_count", json_integer(num_nodes));
    json_object_set_new(out_json, "depth", json_integer(depth));
    json_object_set_new(out_json, "minimum_length", json_integer(snarl_length()));
    json_object_set_new(out_json, "maximum_length", json_integer(max_width));


    return out_json;
}

void MinimumDistanceIndex::SnarlIndex::print_self() const {
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
    cerr << endl << "Length of snarl : " << snarl_length() << endl;

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
            cerr << snarl_distance(n1, n2) << "\t"; 
        }
        cerr << endl;
    }
    cerr << endl; 
}

//ChainDistance methods

MinimumDistanceIndex::ChainIndex::ChainIndex( 
                       id_t parent_id, id_t id_in_parent, id_t end_id, 
                       bool rev_in_parent, bool loops, size_t length):
                parent_id(parent_id), id_in_parent(id_in_parent), end_id(end_id), 
                is_looping_chain(loops), rev_in_parent(rev_in_parent), max_width(0) {
    

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
    sdsl::read_member(end_id, in);
    sdsl::read_member(is_looping_chain, in);
    sdsl::read_member(max_width, in);
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
    sdsl::write_member(end_id, out);
    sdsl::write_member(is_looping_chain, out);
    sdsl::write_member(max_width, out);
   

}
int64_t MinimumDistanceIndex::ChainIndex::loop_distance(
         pair<size_t, bool> start, pair<size_t, bool> end, 
         int64_t start_len, int64_t end_len) const {

    if (start.first == 0 ) {
        return chain_distance(make_pair(prefix_sum.size() - 2, start.second), end, start_len, end_len, true);
    } else if (end.first == 0) {
        return chain_distance(start,make_pair(prefix_sum.size() - 2, end.second), start_len, end_len, true);

    } else if (start.first < end.first && start.second) {

        return chain_distance(start, make_pair(0, start.second),start_len, prefix_sum[0]-1, true) 
               + chain_distance(make_pair(prefix_sum.size() - 2, start.second),end, prefix_sum[prefix_sum.size() - 1]-1, end_len, true);
    } else if (start.first > end.first && !start.second) {
        return chain_distance(start,  make_pair(prefix_sum.size() - 2, start.second), start_len, prefix_sum[prefix_sum.size() - 1]-1, true) 
              + chain_distance(make_pair(0, start.second),end, prefix_sum[0]-1, end_len, true);
    } else {
        return -1;
    }

}
int64_t MinimumDistanceIndex::ChainIndex::chain_distance(
         pair<size_t, bool> start, pair<size_t, bool> end, 
         int64_t start_len, int64_t end_len, bool check_loop) const {

    /*
     * Return the distance between the given node sides, except node side is
     * specified relative to the reading orientation of the chain that the
     * nodes are in. 
     */

    int64_t start_sum = start.first == 0 ? 0 : prefix_sum[start.first] - 1;
    int64_t end_sum = end.first == 0 ? 0 : prefix_sum[end.first] - 1;
    int64_t loop_dist = -1;
    if (is_looping_chain && !check_loop) {
        loop_dist = loop_distance(start, end, start_len, end_len);
    }
    if (!start.second && !end.second) {
        //If start and end are facing forward relative to the start of the chain
        if (start.first <= end.first) {
            return min_pos(loop_dist, end_sum - start_sum);
        } else {
            int64_t rev1 = loop_fd[start.first] - 1;
            int64_t rev2 = loop_rev[end.first] - 1;

            int64_t chain_dist = (start_sum + start_len) - (end_sum + end_len); 
            return min_pos(loop_dist, (rev1 == -1 || rev2 == -1) ? -1 : chain_dist + rev1 + rev2);
        }

    } else if (start.second && end.second ){
        //If start and end are both reversed relative to the start of the chain
        if (start.first >= end.first) {
            return min_pos(loop_dist, (start_sum + start_len) - (end_sum + end_len));
            
        } else {
            int64_t rev1 = loop_rev[start.first] - 1;
            int64_t rev2 = loop_fd[end.first] - 1;
            int64_t chain_dist = end_sum - start_sum; 
            return min_pos(loop_dist, (rev1 == -1 || rev2 == -1) ? -1 : chain_dist+ rev1 + rev2);
        }
    } else if (!start.second && end.second) {
        //Start is forward, end is reversed
        if (start.first <= end.first) {
            int64_t rev = loop_fd[end.first] - 1;
            int64_t chain_dist = end_sum - start_sum;
            return min_pos(loop_dist, rev == -1 ? -1 : rev + chain_dist);
        } else {
            int64_t rev = loop_fd[start.first] - 1;
            int64_t chain_dist = (start_sum+start_len) - (end_sum+end_len);
            return min_pos(loop_dist, rev == -1 ? -1 : rev + chain_dist);
        }
        
    } else {
        //start is reverse, end is forward
        if (start.first <= end.first) {
            int64_t rev = loop_rev[start.first] - 1;
            int64_t chain_dist = end_sum - start_sum;
            return min_pos(loop_dist, rev == -1 ? -1 : rev + chain_dist);

            
        } else {
            int64_t rev = loop_rev[end.first] - 1; 
            int64_t chain_dist = (start_sum+start_len) - (end_sum+end_len);
            return min_pos(loop_dist, rev == -1 ? -1 : rev + chain_dist);
        }
    }
}


json_t*  MinimumDistanceIndex::ChainIndex::chain_to_json() {
    json_t* out_json = json_object();
    json_object_set_new(out_json, "type", json_string("chain"));

    json_t* start_node = json_object();
    json_object_set_new(start_node, "node_id", json_integer(id_in_parent));
    json_t* end_node = json_object();
    json_object_set_new(end_node, "node_id", json_integer(end_id));
    json_t* parent = json_object();
    if (parent_id == 0) {
        json_object_set_new(parent, "type", json_string("root"));
    } else {
        json_object_set_new(parent, "type", json_string("snarl"));
        json_object_set_new(parent, "node_id", json_integer(parent_id));
    }

    json_object_set_new(out_json, "start", start_node);
    json_object_set_new(out_json, "end", end_node);
    json_object_set_new(out_json, "parent", parent);

    json_object_set_new(out_json, "snarl_count", json_integer(prefix_sum.size() - 1));
    json_object_set_new(out_json, "minimum_length", json_integer(prefix_sum[prefix_sum.size() - 1]));
    json_object_set_new(out_json, "maximum_length", json_integer(max_width));

    return out_json;
}
void MinimumDistanceIndex::ChainIndex::print_self() const {
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

void MinimumDistanceIndex::calculate_max_index(const HandleGraph* graph, int64_t cap) {

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
    bool is_acyclic = handlealgs::is_directed_acyclic(&split_graph);

    //If the graph has directed cycles, dagify it
    VG dagified; 
    HandleGraph* dag;
    if (!is_acyclic) {
#ifdef debugIndex
cerr << "Making graph acyclic" << endl;
#endif
        acyclic_to_id = handlealgs::dagify(&split_graph, &dagified, cap);
        dag = &dagified;
    } else {
        dag = &split_graph;
    }

     */
    bdsg::HashGraph split;
    {
        // TODO: hack to avoid recoding this with the new return type
        auto split_to_id_tmp = handlealgs::split_strands(graph, &split);
        
        split_to_id.reserve(split_to_id_tmp.size());
        for (const auto& trans : split_to_id_tmp) {
            split_to_id[split.get_id(trans.first)] = make_pair(graph->get_id(trans.second),
                                                               graph->get_is_reverse(trans.second));
        }
    }
    graph = &split;


    unordered_map<id_t, id_t> acyclic_to_id;
    bool is_acyclic = handlealgs::is_directed_acyclic(&split);

    //If the graph has directed cycles, dagify it
    HandleGraph* dag;
    bdsg::HashGraph dagified; 
    if (!is_acyclic) {
#ifdef debugIndex
cerr << "Making graph acyclic" << endl;
#endif
        acyclic_to_id = handlealgs::dagify(&split, &dagified, cap);
        dag = &dagified;
    } else {
        dag = &split;
    }

#ifdef debugIndex
    assert(handlealgs::is_single_stranded(dag));
    assert(handlealgs::is_directed_acyclic(dag));
#endif

    //In DAGified graph, go through graph in topological order and get the 
    //minimum and maximum distances to tips 

    vector<handle_t> order = handlealgs::topological_order(dag);

    for (handle_t curr_handle : order) {

        //Get the correct id in the original graph
        id_t curr_id = dag->get_id(curr_handle);
        if (!is_acyclic) {
            curr_id = acyclic_to_id[curr_id];
        }
        curr_id = split_to_id[curr_id].first;

        //Get the previous values for this node
        int64_t curr_min = min_distances[curr_id-min_node_id];
        curr_min = curr_min == 0 ? std::numeric_limits<int64_t>::max() : curr_min; 
                              
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


void MinimumDistanceIndex::write_snarls_to_json() {
    for (auto& snarl : snarl_indexes) {

        cout << json_dumps(snarl.snarl_to_json(), JSON_ENCODE_ANY) << endl;
    }
    for (auto& chain : chain_indexes) {
        json_t* chain_json = chain.chain_to_json();
        size_t depth = snarl_indexes[get_primary_assignment(chain.id_in_parent)].depth; 
        json_object_set_new(chain_json, "depth", json_integer(depth));
        cout << json_dumps(chain_json, JSON_ENCODE_ANY) << endl;
    }

}

void MinimumDistanceIndex::print_snarl_stats() {
    //Print out stats bout the snarl
    //

    cout << "#Distance index has " << snarl_indexes.size() << " snarls and " << chain_indexes.size() << " chains" << endl;
    cout << "#start_id\tend_id\tsnarl_size\tsnarl_depth" << endl;
   
    for (auto snarls : snarl_indexes) {
       cout << snarls.id_in_parent << "\t" << snarls.end_id << "\t" << snarls.num_nodes << "\t" << snarls.depth << endl; 
    }
}

void MinimumDistanceIndex::add_nodes_in_range(const HandleGraph* super_graph, int64_t min_distance, int64_t max_distance, 
                                      std::unordered_set<id_t>& sub_graph, vector<tuple<handle_t, int64_t>>& start_nodes, 
                                      hash_set<pair<id_t, bool>>& seen_nodes) {
    //Starting from a given handle in the super_graph, traverse the graph and add all nodes within the distance range to sub_graph
    
#ifdef debugSubgraph
    cerr << "Starting search from nodes ";
    for (auto& start_handle : start_nodes) {
        cerr << super_graph->get_id(std::get<0>(start_handle)) << " " << super_graph->get_is_reverse(std::get<0>(start_handle))
             << " with distance " << std::get<1>(start_handle) << endl;
    }
#endif
    auto cmp =  [] (const tuple<handle_t, int64_t> a, const tuple<handle_t, int64_t> b ) {
            return std::get<1>(a) > std::get<1>(b);
        };
    priority_queue< tuple<handle_t, int64_t>, vector<tuple<handle_t, int64_t>>, decltype(cmp)> next_handles (cmp);
    for (auto& start_handle : start_nodes) {
        next_handles.emplace(start_handle);
    }
    bool first_node = true;

    while (next_handles.size() > 0) {
        //Traverse the graph, adding nodes if they are within the range
        handle_t curr_handle; int64_t curr_distance;
        std::tie(curr_handle, curr_distance) = next_handles.top();
        next_handles.pop();
        if (seen_nodes.count(make_pair(super_graph->get_id(curr_handle), super_graph->get_is_reverse(curr_handle))) == 0) {
            seen_nodes.emplace(super_graph->get_id(curr_handle), super_graph->get_is_reverse(curr_handle));

            int64_t node_len = super_graph->get_length(curr_handle);
            int64_t curr_distance_end = curr_distance + node_len - 1;
            if ((curr_distance >= min_distance && curr_distance <= max_distance) ||
                 (curr_distance_end >= min_distance && curr_distance_end <= max_distance) ||
                 (curr_distance <= min_distance && curr_distance_end >= max_distance)) {
#ifdef debugSubgraph
                cerr << "\tadding node " << super_graph->get_id(curr_handle) << " " << super_graph->get_is_reverse(curr_handle) << " with distance "
                     << curr_distance << " and node length " << node_len << endl;
#endif
                sub_graph.insert(super_graph->get_id(curr_handle));
               
            }
#ifdef debugSubgraph
            else {
                cerr << "\tdisregarding node " << super_graph->get_id(curr_handle) << " " << super_graph->get_is_reverse(curr_handle) 
                     << " with distance " << curr_distance << " and node length " << node_len << endl;
            }
#endif
            
            curr_distance += node_len;

            //If the end of this node is still within the range, add the next nodes that are within 
            if (curr_distance-1 <= max_distance ) {
                super_graph->follow_edges(curr_handle, false, [&](const handle_t& next) {
                    id_t next_id = super_graph->get_id(next);
                    if (seen_nodes.count(make_pair(next_id, super_graph->get_is_reverse(next))) == 0) {
                        next_handles.emplace(next, curr_distance); 
                    }
                    return true;
                });
            }
            first_node = false;
        }

    }
}

void MinimumDistanceIndex::subgraph_in_range(const Path& path, const HandleGraph* super_graph, int64_t min_distance, 
                                           int64_t max_distance, std::unordered_set<id_t>& sub_graph, bool look_forward){

    //Get the subgraph of all nodes for which the minimum distance to any position in the node is within the distance range
    //Algorithm proceeds in two phases: First, traverse up the snarl tree and get the distance to the ends of each snarl
    //or chain. When the distance from either end of the current structure plus the distance to any other node in the parent
    //is greater than the lower end of the distance range, start a traversal of the graph from the ends of the current node
    //and look for all nodes within the distance range (add_nodes_in_range)

    pos_t start_pos;
    int64_t node_len;
    if (look_forward ){
        start_pos = initial_position(path);
        node_len = super_graph->get_length(super_graph->get_handle(get_id(start_pos)));
    } else {
        start_pos = final_position(path);
        node_len = super_graph->get_length(super_graph->get_handle(get_id(start_pos)));
        start_pos = reverse_base_pos(start_pos, node_len);
    }


    //The distance from the position to the ends of the current node(/snarl/chain)
    int64_t dist_to_curr_start = is_rev(start_pos) ? node_len - get_offset(start_pos) + 1 : -1 ;
    int64_t dist_to_curr_end = is_rev(start_pos) ? -1 : node_len - get_offset(start_pos) + 1; 

    //Graph node of the start and end of the current node(/snarl/chain) pointing out
    pair<id_t, bool> start_node (get_id(start_pos), true); 
    pair<id_t, bool> end_node (get_id(start_pos), false); 
    int64_t start_len = node_len; int64_t end_len = node_len;

    //The rank of the current node in its parent 
    size_t curr_rank = get_primary_rank(get_id(start_pos));
    bool curr_rev = false;
    bool passed_root = false;//This becomes true when we need to stop traversing up the tree
    //The parent of the current node. bool is true if it is a chain, size_t is the index to snarl/chain_indexes
    pair<bool, size_t> parent_structure (false, get_primary_assignment(get_id(start_pos)));

#ifdef debugSubgraph
    cerr << endl << endl << endl << "Finding subgraph within range " << min_distance << "-" << max_distance << " starting from " << start_pos << endl;
    cerr << "\tstart with distances on node " << dist_to_curr_start << " " << dist_to_curr_end << endl;
#endif

    vector<tuple<handle_t, int64_t>> search_start_nodes;
    hash_set<pair<id_t, bool>> seen_nodes;//Nodes that are too close and should be avoided when doing search
    //TODO: Need this to pass unit tests, since the distance from a node to itself is 0, but might not actually want it
    seen_nodes.emplace(start_node.first, is_rev(start_pos));

    //Walk up the snarl tree until the distance to the end of the containing structure is at least the min_distance
    while (!passed_root) {
        //TODO: I think this should just walk all the way up the snarl tree 
        if (parent_structure.first) {
            //If the parent is a chain
            
            const ChainIndex& chain_index = chain_indexes[parent_structure.second];

            //Find the rank of the end node and current node in the chain
            size_t chain_end_rank = chain_index.prefix_sum.size() - 2;
                                    
            //Get the lengths of start, end, and current node

            int64_t chain_start_len = chain_index.prefix_sum[0]-1;

            int64_t chain_end_len =  chain_index.prefix_sum[chain_index.prefix_sum.size()-1] 
                                   - chain_index.prefix_sum[chain_index.prefix_sum.size()-2];

            int64_t dsl = chain_index.chain_distance(make_pair(0, false), make_pair(curr_rank, false), chain_start_len, start_len);
            int64_t dsr = chain_index.chain_distance(make_pair(0, false),  make_pair(curr_rank + 1, true), chain_start_len, end_len);
            int64_t der = chain_index.chain_distance(make_pair(chain_end_rank, true), make_pair(curr_rank + 1, true), chain_end_len, end_len);
            int64_t del = chain_index.chain_distance(make_pair(chain_end_rank, true),make_pair(curr_rank, false),chain_end_len, start_len);

            dsl = dsl == -1 || dist_to_curr_start == -1 ? -1 : dist_to_curr_start + dsl;
            dsr = dsr == -1 || dist_to_curr_end == -1 ? -1 : dist_to_curr_end + dsr;
            der = der == -1 || dist_to_curr_end == -1 ? -1 : dist_to_curr_end + der;
            del = del == -1 || dist_to_curr_start == -1 ? -1 : dist_to_curr_start + del;
 
            int64_t dist_start = min_pos(dsr, dsl);
            int64_t dist_end = min_pos(der, del);

            if ((chain_index.max_width < min_distance) && parent_structure.second != 0) {
                //If we haven't reached the min distance yet, update distances and nodes
                dist_to_curr_start = dist_start;
                dist_to_curr_end   =  dist_end;
#ifdef debugSubgraph
                cerr << "At node " << start_node.first << "(" << start_len << ")" << "->" << end_node.first << "(" << end_len << ") "
                     << " in chain " << chain_index.id_in_parent << ", reached ends with distances " << dist_to_curr_start << " and " 
                     << dist_to_curr_end << ", moving up snarl tree..." << endl;
#endif
                if (chain_index.parent_id == 0) {
                    passed_root = true;
                } else {
                    parent_structure = make_pair(false, get_primary_assignment(chain_index.parent_id));
                    curr_rank = get_secondary_rank(chain_index.id_in_parent);
                    start_node = make_pair(chain_index.id_in_parent, !get_primary_rank(chain_index.id_in_parent) % 2);
                    end_node = make_pair(chain_index.end_id, get_primary_rank(chain_index.end_id) % 2); 
                    start_len = chain_start_len;
                    end_len = chain_end_len;

                }

            } else {
                //The distance to at least one end of the chain passes the lower boundary of the distance range, so we 
                //will start a search from inside the chain
#ifdef debugSubgraph
                cerr << "At node " << start_node.first << "(" << start_len << ")" << "->" << end_node.first << "(" << end_len << ") "
                     << " in chain " << chain_index.id_in_parent  << ", reached ends with distances " 
                     << dist_start << " and " << dist_end << ", within distance range" << endl;
#endif

                //Update left and right distances with the loop distances in the chain
                int64_t right_loop = chain_index.chain_distance(make_pair(curr_rank+1, false), make_pair(curr_rank, true), end_len, start_len);
                int64_t left_loop = chain_index.chain_distance(make_pair(curr_rank, true), make_pair(curr_rank+1, false), start_len, end_len); 
                int64_t tmp_dist = dist_to_curr_end;
                dist_to_curr_end = left_loop == -1 || dist_to_curr_start == -1 ? dist_to_curr_end : min(dist_to_curr_end, dist_to_curr_start + left_loop + end_len);
                dist_to_curr_start = right_loop == -1 || tmp_dist == -1 ? dist_to_curr_start : min(dist_to_curr_start, tmp_dist + right_loop+start_len);

#ifdef debugSubgraph
                cerr << "  With distances to ends of snarl " << dist_to_curr_start << " and " << dist_to_curr_end << endl;
#endif

                vector<tuple<handle_t, int64_t>> search_start_nodes;

                    
#ifdef debugSubgraph
                cerr << "\tsearching left in chain starting from rank " << curr_rank << endl;
#endif
                bool got_start = dist_to_curr_start == -1;
                pair<id_t, bool> last_end = end_node;
                for (size_t i = curr_rank ; i > 0 && !got_start; i --) {
                    //Start from i = the start node of the original snarl
                    //The current snarl is defined by (index i-1 -> i) and go left
                    SnarlIndex& snarl_index = snarl_indexes[get_secondary_assignment(start_node.first)];
                    id_t next_start = snarl_index.id_in_parent; //Start id of this snarl (i)
                    int64_t next_len = snarl_index.node_length(get_primary_rank(next_start)); //len of i
                    int64_t prev_len = snarl_index.node_length(get_secondary_rank(start_node.first));//len of i-1

                    //The maximum minimum distance from the start position to any node in the current snarl
                    int64_t max_dist_left = dist_to_curr_start == -1 ? -1 : dist_to_curr_start + snarl_index.max_width - prev_len;
                    //Distance to go into the snarl and back out from i
                    int64_t loop_dist = snarl_index.rev_in_parent 
                                        ? snarl_index.snarl_distance(0, 1) 
                                        : snarl_index.snarl_distance(snarl_index.num_nodes*2 - 1, snarl_index.num_nodes*2-2);
#ifdef debugSubgraph
                    cerr << "\t\tat rank " << i << " at snarl between nodes " << start_node.first << " and " << next_start 
                         << " the max dist to node " << next_start << " is " << max_dist_left  << endl;
#endif

                    if ((max_dist_left != -1 && max_dist_left >= min_distance) || 
                        (loop_dist != -1)){
                        //If the start of this snarl is within the distance range, start search from end of the snarl (start_node.first)


#ifdef debugSubgraph
                        cerr << "\t Add node to start search: " << start_node.first << " " << start_node.second 
                             << " distance to the start(?) of the node is " << dist_to_curr_start << endl;
#endif
                        search_start_nodes.emplace_back(super_graph->get_handle(start_node.first, start_node.second),
                                                        dist_to_curr_start-prev_len);
                        seen_nodes.erase(start_node);
                        got_start = true;
                    }
                    dist_to_curr_start = dist_to_curr_start == -1 ? -1 : snarl_index.snarl_length() + dist_to_curr_start - prev_len;
                    last_end = make_pair(start_node.first, !start_node.second);
                    start_node = make_pair(next_start, !(get_primary_rank(next_start) % 2));
                }
                if (!got_start){ 
                    search_start_nodes.emplace_back(super_graph->get_handle(chain_index.id_in_parent, 
                                !get_primary_rank(chain_index.id_in_parent) % 2), dist_start-chain_start_len+1);
                    seen_nodes.erase(make_pair(chain_index.id_in_parent, !get_primary_rank(chain_index.id_in_parent) % 2));
                }
#ifdef debugSubgraph
                cerr << "\tsearching right in chain:  " << endl;
#endif
                bool got_end = dist_to_curr_end == -1;
                for (size_t i = curr_rank + 2 ; i < chain_index.prefix_sum.size() - 1 && !got_end ; i++) {
                    //Start with i being the end of the snarl after the original snarl
                    //Looking at snarl with indices i-1 to i
                    SnarlIndex& snarl_index = snarl_indexes[get_primary_assignment(end_node.first)];
                    id_t next_end = snarl_index.end_id; // i+1, end_node is i
                    int64_t last_len = snarl_index.node_length(get_primary_rank(end_node.first));

                    int64_t max_dist_right = dist_to_curr_end == -1 ? -1 : dist_to_curr_end + snarl_index.max_width - last_len;
                    //Distance to go into the snarl and back out from i-1
                    int64_t loop_dist = snarl_index.rev_in_parent 
                                        ? snarl_index.snarl_distance(snarl_index.num_nodes*2 - 1, snarl_index.num_nodes*2-2)
                                        : snarl_index.snarl_distance(0, 1) ;
#ifdef debugSubgraph
                    cerr << "\t\tat rank " << i << " snarl between " << end_node.first << " and " << next_end << " the max dist to node " << next_end << " is " << max_dist_right << " and the loop dist is " << loop_dist << endl;
#endif

                    if ((max_dist_right != -1 &&  max_dist_right >= min_distance) || 
                        (loop_dist != -1 )) {

#ifdef debugSubgraph
                        cerr << "\t Add node to start search:  " << end_node.first << " " << end_node.second 
                             <<  " dist is " << dist_to_curr_end << endl;
#endif
                        search_start_nodes.emplace_back(super_graph->get_handle(end_node.first, end_node.second), dist_to_curr_end-last_len);
                        seen_nodes.erase(end_node);

                        got_end = true;
                    }
                    dist_to_curr_end = dist_to_curr_end == -1 ? -1 : dist_to_curr_end + snarl_index.snarl_length() - last_len;
                    end_node = make_pair(next_end, get_secondary_rank(next_end) % 2);
                }
                if (!got_end) {
                    search_start_nodes.emplace_back(super_graph->get_handle(chain_index.end_id, get_primary_rank(chain_index.end_id) % 2), 
                                                    dist_end-chain_end_len);
                    seen_nodes.erase(make_pair(chain_index.end_id, get_primary_rank(chain_index.end_id) % 2));
                }

                add_nodes_in_range(super_graph, min_distance, max_distance, sub_graph, search_start_nodes, seen_nodes);
                return;
            }
            start_len = chain_start_len;
            end_len = chain_end_len;
        } else {
            //If the parent is a snarl
            
            const SnarlIndex& snarl_index = snarl_indexes[parent_structure.second];
            pair<int64_t, int64_t> new_dists = snarl_index.dist_to_ends(curr_rank, dist_to_curr_start, dist_to_curr_end);

            if ((new_dists.first != -1 && new_dists.first > min_distance) || (new_dists.second != -1 && new_dists.second > min_distance)
                || ( snarl_index.max_width > min_distance) 
                || (snarl_index.is_unary_snarl && !snarl_index.in_chain && snarl_indexes[get_primary_assignment(snarl_index.parent_id)].is_unary_snarl)) {
                //If this goes past the minimum distance or the width of the snarl plus distance we've found already is in the range
                //Or this is a nested unary snarl
#ifdef debugSubgraph
                cerr << "At node " << start_node.first << "(" << start_len << ")" << "->" << end_node.first << "(" << end_len << ") "
                     << " in snarl between " << snarl_index.id_in_parent << " and " << snarl_index.end_id 
                     << ", reached ends with distances " << new_dists.first << " and " << new_dists.second  << " within range " << endl;
#endif
                vector<tuple<handle_t, int64_t>> search_start_nodes;
                
                if (dist_to_curr_start != -1 ) {
                    search_start_nodes.emplace_back(super_graph->get_handle(start_node.first, start_node.second), 
                                                    dist_to_curr_start-start_len);
                    seen_nodes.erase(start_node);
#ifdef debugSubgraph
                    cerr << "\t Add snarl start node to start search: " << start_node.first << " " << start_node.second <<  " dist is " << dist_to_curr_start << " and node len is " << start_len << endl;
#endif
                }
                if (dist_to_curr_end != -1 ) {
                    search_start_nodes.emplace_back(super_graph->get_handle(end_node.first, end_node.second),
                                                    dist_to_curr_end-end_len);
                    seen_nodes.erase(end_node);
#ifdef debugSubgraph
                    cerr << "\tAdd snarl end node to start search: " << end_node.first << " " << end_node.second <<  " dist is " << dist_to_curr_end 
                         << " and node len is " << end_len << endl;
#endif
                }
                add_nodes_in_range(super_graph, min_distance, max_distance, sub_graph, search_start_nodes, seen_nodes);
                return;
            } else {
                //Update current distances and nodes
#ifdef debugSubgraph
                cerr << "At node " << start_node.first << "(" << start_len << ")" << "->" << end_node.first << "(" << end_len << ") "
                     << " in " << (snarl_index.is_unary_snarl ? "unary" : "" ) << " snarl " 
                     << snarl_index.id_in_parent <<  ", reached ends with distances " << new_dists.first << " and " 
                     << new_dists.second << " from distances " << dist_to_curr_start << " and " << dist_to_curr_end << ", moving up snarl tree.. " << endl;
#endif


                //We can reach the ends of this snarl without hitting the distance range, so make sure we don't look at the nodes
                //immediately after this node again
                //TODO: This prevents us from reaching nodes that are too close from their minimum distance path from 
                //the original position but we could still reach them from other paths
                pair<id_t, bool> snarl_start (snarl_index.id_in_parent, get_primary_rank(snarl_index.id_in_parent)%2==0);
                pair<id_t, bool> snarl_end (snarl_index.end_id, get_primary_rank(snarl_index.end_id)%2==1);
                if (dist_to_curr_start != -1 && start_node != snarl_start && start_node != snarl_end){
                    super_graph->follow_edges(super_graph->get_handle(start_node.first, start_node.second), false, [&] (const handle_t next) {
                        seen_nodes.emplace(super_graph->get_id(next), super_graph->get_is_reverse(next));
#ifdef debugSubgraph
                        cerr << "Adding seen node " << super_graph->get_id(next) << " " << super_graph->get_is_reverse(next) << endl;
#endif
                        return true;
                    });
                }
                if (dist_to_curr_end != -1 && end_node != snarl_start && end_node != snarl_end){
                    super_graph->follow_edges(super_graph->get_handle(end_node.first, end_node.second), false, [&] (const handle_t next) {
                        seen_nodes.emplace(super_graph->get_id(next), super_graph->get_is_reverse(next));
#ifdef debugSubgraph
                        cerr << "Adding seen node " << super_graph->get_id(next) << " " << super_graph->get_is_reverse(next) << endl;
#endif
                        return true;
                    });
                }
                dist_to_curr_start = snarl_index.rev_in_parent ? new_dists.second : new_dists.first;
                dist_to_curr_end = snarl_index.rev_in_parent ? new_dists.first : new_dists.second;

                start_node = std::move(snarl_start);
                end_node = std::move(snarl_end) ; 
                start_len = snarl_index.node_length(0);
                end_len = snarl_index.node_length(snarl_index.num_nodes*2 - 1);

                if (snarl_index.parent_id == 0) {
                    passed_root = true;
                } else {
                    if (snarl_index.is_unary_snarl) {
                        end_node = start_node;
                        end_len = start_len;
                        dist_to_curr_start = min_pos(dist_to_curr_start, dist_to_curr_end);
                        dist_to_curr_end = dist_to_curr_start;
                    }
                    if (snarl_index.in_chain) {
                        //If the parent is a chain
                        parent_structure = make_pair(true, get_chain_assignment(snarl_index.parent_id));
                        curr_rank = get_chain_rank(snarl_index.id_in_parent);

                    } else {
                        parent_structure = make_pair(false, get_primary_assignment(snarl_index.parent_id));
                        curr_rank = get_secondary_rank(snarl_index.id_in_parent);
                        if (snarl_index.rev_in_parent) {

                            auto tmp = start_node;
                            start_node = end_node;
                            end_node = tmp;

                            int64_t tmp_len = start_len;
                            start_len = end_len;
                            end_len = tmp_len;
                        }
                    }
                }
            }
        }
    }
    return;
}

tuple<bool, size_t, size_t, bool, size_t, size_t, size_t, size_t, bool> MinimumDistanceIndex::get_minimizer_distances (pos_t pos) {
    if (node_to_component.size() == 0) {
        throw runtime_error("error: distance index is out-of-date");
    }
    id_t id = get_id(pos);
    SnarlIndex& snarl_index = snarl_indexes[get_primary_assignment(id)];
    size_t snarl_rank =  get_primary_rank(id);

    bool is_boundary_node = snarl_rank == 0 || snarl_rank == 1 || 
                            snarl_rank == snarl_index.num_nodes*2-1 || snarl_rank == snarl_index.num_nodes*2-2;
    size_t component = node_to_component[id - min_node_id]; 


    if (component != 0 && snarl_index.depth == 0 && snarl_index.in_chain && is_boundary_node && !chain_indexes[get_chain_assignment(snarl_index.id_in_parent)].is_looping_chain ) {
        //If this node is a boundary node of a top-level chain
        int64_t node_offset = get_offset(pos);
        bool node_is_rev_in_snarl = snarl_rank% 2;
        node_is_rev_in_snarl = is_rev(pos) ? !node_is_rev_in_snarl : node_is_rev_in_snarl;
        bool node_is_rev_in_chain = node_is_rev_in_snarl ? !snarl_index.rev_in_parent : snarl_index.rev_in_parent;
        if (node_is_rev_in_chain){
            node_offset = snarl_index.node_length(snarl_rank) - node_offset;
        } else {
            node_offset += 1;
        }
  
        size_t length = component_to_chain_length[component-1];
        size_t chain_rank = get_chain_rank(id); 
        size_t offset = chain_rank == 0 ? 0 : chain_indexes[component_to_chain_index[component-1]].prefix_sum[chain_rank] - 1;

        return tuple<bool, size_t, size_t, bool, size_t, size_t, size_t, size_t, bool>(true, component, offset + node_offset,
            false, MinimumDistanceIndex::MIPayload::NO_VALUE, MinimumDistanceIndex::MIPayload::NO_VALUE, MinimumDistanceIndex::MIPayload::NO_VALUE, MinimumDistanceIndex::MIPayload::NO_VALUE, false);

    } else if (component != 0 && snarl_index.depth == 0 && snarl_index.in_chain && snarl_index.is_simple_snarl && !chain_indexes[get_chain_assignment(snarl_index.id_in_parent)].is_looping_chain ) {
        //This node is on a top-level simple snarl

        size_t chain_rank = get_chain_rank(snarl_index.id_in_parent);
        size_t start_len = snarl_index.rev_in_parent ? snarl_index.node_length(snarl_index.num_nodes*2-1) : snarl_index.node_length(0);
        size_t end_len = snarl_index.rev_in_parent ? snarl_index.node_length(0) : snarl_index.node_length(snarl_index.num_nodes*2-1);
        size_t node_len = snarl_index.node_length(snarl_rank);
        bool rev_in_snarl = snarl_index.snarl_distance(0, snarl_rank) == -1;
        bool is_rev = rev_in_snarl ? !snarl_index.rev_in_parent : snarl_index.rev_in_parent;
        return tuple<bool, size_t, size_t, bool, size_t, size_t, size_t, size_t, bool>(false, MinimumDistanceIndex::MIPayload::NO_VALUE, MinimumDistanceIndex::MIPayload::NO_VALUE,
                true, chain_rank, start_len, end_len, node_len, is_rev);

    } else {
        //If this is a nested position
        return tuple<bool, size_t, size_t, bool, size_t, size_t, size_t, size_t, bool>(false, MinimumDistanceIndex::MIPayload::NO_VALUE, MinimumDistanceIndex::MIPayload::NO_VALUE,
                false, MinimumDistanceIndex::MIPayload::NO_VALUE, MinimumDistanceIndex::MIPayload::NO_VALUE, MinimumDistanceIndex::MIPayload::NO_VALUE, MinimumDistanceIndex::MIPayload::NO_VALUE, false);
    }
}

int64_t MinimumDistanceIndex::top_level_chain_length(id_t node_id) {
    if (node_to_component.size() == 0) {
        throw runtime_error("error: distance index is out-of-date");
    }
    size_t component = node_to_component[node_id-min_node_id];
    return component == 0 ? -1 : component_to_chain_length[component-1];
}
size_t MinimumDistanceIndex::get_connected_component(id_t node_id) {
    if (node_to_component.size() == 0) {
        throw runtime_error("error: distance index is out-of-date");
    }
    return node_to_component[node_id-min_node_id];
}

constexpr MinimumDistanceIndex::MIPayload::code_type MinimumDistanceIndex::MIPayload::NO_CODE;
constexpr size_t MinimumDistanceIndex::MIPayload::NO_VALUE;

}
