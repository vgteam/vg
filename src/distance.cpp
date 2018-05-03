#include "distance.hpp"
using namespace std;
namespace vg {

/*
SnarlDistances::SnarlDistances(unordered_set<Node*>* allNodes) {
    //Assign all nodes+direction in snarl to an index
    unordered_map<pair<id_t, bool>, size_t> visit_to_index_obj;
    visit_to_index = &visit_to_index_obj;
    size_t snarlIndex = 0;
    for (unordered_set<Node*>::iterator n = allNodes->begin(); 
            n != allNodes->end(); ++n) {

        const Node* node = *n;
        id_t nID = node->id();

        visit_to_index->insert(make_pair(make_pair(nID, true), snarlIndex++));
        visit_to_index->insert(make_pair(make_pair(nID, false), snarlIndex++));
    }

    int size = visit_to_index->size();
    //Initialize all distances to -1
    vector<int64_t> distance_obj (size*size, -1);
    distances = &distance_obj;//((size+1) * size / 2,-1);

}
*/
SnarlDistances::SnarlDistances(unordered_map<pair<id_t, bool>, size_t>* vti, 
                   vector<int64_t>* dists) {
    visit_to_index = vti;
    distances = dists;
    assert(distances->size() == visit_to_index->size() ^ 2);
    
}
size_t SnarlDistances::index(pair<id_t, bool> start, pair<id_t, bool> end) {
    /*Get the index of dist from start to end in a snarl distance matrix
      given the node ids + direction */
    size_t length = visit_to_index->size();
    size_t i1 = visit_to_index->at(start);//TODO Sometimes doesn't work?????
    size_t i2 = visit_to_index->at(end);
    //return (length + 1) * (length / 2) + (i1+1)* (i1/2) + i2 - 1;
    /*TODO Currently store entire distance matrix since it's simpler
       dist from n1 -> n2 = -n2 -> -n1 - len(n1) + len(n2) */
    return i1 * length + i2;
}
   
int64_t SnarlDistances::snarlDistance(pair<id_t, bool> start,
                                                      pair<id_t, bool> end) {
    /*Distance between beginnings of two nodes n1 and n2 in snarl
    */
    size_t i = index(start, end);
    return distances->at(i); 
}

void SnarlDistances::insertDistance(pair<id_t, bool> start, 
                                           pair<id_t, bool> end, int64_t dist) {
    //Assign distance between start and end
    size_t i = index(start, end);
    distances->at(i) = dist;
}

int64_t SnarlDistances::snarlDistance_short(VG* graph, pair<id_t, bool> start,
                                                      pair<id_t, bool> end) {
    /*Distance between end of node n1 and beginning of node n2 in snarl
    */
    size_t i = index(start, end);
    int64_t totalDist = distances->at(i); 
    if (totalDist == -1) {
        return -1;

    } else {
       /* If there is a path, return the distance between two nodes - length of
         first node */
       handle_t handle = graph->get_handle(start.first, start.second);
       vector<const handle_t*> nextNodes;//vector of next nodes
       auto addHandle = [&](const handle_t& h)-> bool {
           nextNodes.push_back(&h);
           return true;
       };
       graph->follow_edges(handle, false, addHandle);
       handle_t nextHandle = *nextNodes[0];
       size_t j = index(start, make_pair(graph->get_id(nextHandle), 
                                         graph->get_is_reverse(nextHandle)));

       int64_t firstDist = distances->at(j);//Length of first node 
       return totalDist - firstDist;
    }
}

ChainDistances::ChainDistances(unordered_map<id_t, size_t>* s, 
                                                          vector<int64_t>* p) {
    
    snarl_to_index = s;
    prefix_sum = p;
}

int64_t ChainDistances::chainDistance(pair<id_t, bool> start, 
                                                      pair<id_t, bool> end) {
    //Distance between start of start and start of end in a chain
    size_t i1 = snarl_to_index->at(start.first);
    size_t i2 = snarl_to_index->at(end.first); 
    if ((start.second && end.second) || (!start.second && !end.second)) {
        //If start and end are facing the same direction
        int64_t dNoRev = abs(prefix_sum->at(i1) - prefix_sum->at(i2)); 
        return dNoRev;
    } else {
        //TODO Currently assumes no reversing edges in the chain
        return -1;
    }
}

int64_t ChainDistances::chainDistance_short(VG* graph, pair<id_t, bool> start,
                                                      pair<id_t, bool> end) {
    //Distance between end of start node to beginning of end node in chain
    
    size_t i1 = snarl_to_index->at(start.first);
    size_t i2 = snarl_to_index->at(end.first); 
    int64_t d1 = chainDistance(start, end);
    int64_t d2 = chainDistance(end, start);
    if (d1 == -1 && d2 == -1) {
        return -1;
    } else if (d1 == -1) {
        return d2 - graph->get_node(i1)->sequence().size();
    } else if (d2 == -1) {
        return d1 - graph->get_node(i2)->sequence().size();
    } else {
        return min(d2 - graph->get_node(i1)->sequence().size(), 
                   d1 - graph->get_node(i2)->sequence().size());
    }
}

int64_t ChainDistances::chainLength() {
    //Get the length of a chain including length of last node
    return prefix_sum->back();
}


int64_t calculateIndex(DistanceIndex* di, VG* graph, SnarlManager* sm, 
         const Chain* chain) {
    /*Add chain to distance index and recursively add all distances in 
      snarls contained within chain
    */
    auto cmp = [] (pair<const handle_t*,int64_t> x,
                                           pair<const handle_t*,int64_t> y) {
        //Comparison function for the priority of a pair of handle, distance
        return (x.second < y.second); 
    };

    vector<int64_t> chain_prefix_sum (1,0); //initialize to [0]

    unordered_map<id_t, size_t> chain_snarl_to_index;
    for (int i = 0; i < chain->size(); i++) { 
        //for each snarl in the chain TODO ChainIterator?
        const Snarl* snarl = chain->at(i);
        chain_snarl_to_index.insert( make_pair<id_t, size_t>
                     (snarl->start().node_id(), chain_snarl_to_index.size()));
        //TODO use visit as key to identify a node+direction? or node id? 
        //Get all nodes in snarl includes both boundary nodes of child snarls
        pair<unordered_set<Node*>, unordered_set<Edge*>> contents = 
                        sm->shallow_contents(snarl, *graph, true);
        unordered_set<Node*> allNodes = contents.first;
        /*TODO better way of getting all the nodes in the snarl? Currently
          don't include end boundary node of a child snarl in distance 
          matrix because I don't think I need them?
        */ 

         // Create SnarlDistance object with all nodes in this snarl
        unordered_map <pair<id_t, bool>, size_t> vti;
        size_t snarlIndex = 0;  
        for (unordered_set<Node*>::iterator n = allNodes.begin(); 
            n != allNodes.end(); ++n) {

            const Node* node = *n;
            id_t nID = node->id();

            vti.insert(make_pair(make_pair(nID, true), vti.size()));
            vti.insert(make_pair(make_pair(nID, false), vti.size()));
        }
        int size = vti.size();
        //Initialize all distances to -1
        vector<int64_t> distance_obj (size*size, -1);
        SnarlDistances sd = SnarlDistances(&vti, &distance_obj);
        //add this snarl distance to distance index
        di->sd.insert(make_pair<id_t, SnarlDistances*>
                                                (snarl->start().node_id(),&sd));

        //Create a NetGraph for current snarl
        NetGraph ng = 
            NetGraph(snarl->start(), snarl->end(), sm->chains_of(snarl), graph);


        // For each node in snarl calculate distance to every reachable node
        for (unordered_set<Node*>::iterator n = allNodes.begin(); 
                                                     n != allNodes.end(); ++n){
        bool bools [2] = {true, false};
        for (bool rev : bools) {
        //TODO Better way to loop over all nodes/direction??

            Node* startNode = *n;
            pair<id_t, bool> startID (startNode->id(), rev); 
            const handle_t handle = ng.get_handle(startNode->id(), rev);

            //Priority queue of reachable nodes (handles)
            priority_queue<   pair<const handle_t*, int64_t>,  
                         vector<pair<const handle_t*, int64_t>>, 
                                     decltype(cmp)> reachable(cmp);
            pair<const handle_t*, int64_t> init (&handle, 0);
            reachable.push(init);

            while (reachable.size() > 0) {
                pair<const handle_t*, int64_t> next = reachable.top();
                reachable.pop();
                handle_t currHandle = *next.first;
                int64_t currDist = next.second;
                pair<id_t, bool> currID (ng.get_id(currHandle),
                                            ng.get_is_reverse(currHandle));
            
                if ( sd.snarlDistance(startID, currID) == -1) {//??????
                    //If node has not already been found:
                    //Save distance from start to current node 
                    sd.insertDistance(startID, currID, currDist);
                   
                    //Get the length of the current node
                    int64_t nodeLen;

                    const Snarl* currSnarl = sm->into_which_snarl(
                                                  currID.first, currID.second);
                    //TODO:If is_child, then this shouldn't be null?????
                    // assert(currSnarl != NULL);
                    if (currID.first != snarl->start().node_id() &&
                                  currSnarl != NULL && ng.is_child(currHandle)){
                    //If a child snarl/chain begins at the current node

                        if (sm->in_nontrivial_chain(currSnarl)) {//Chain
                            /*TODO assuming start and end are consistent 
                              within the graph, not relative to the current
                              orientation
                            */
                            const Chain* currChain = sm->chain_of(currSnarl);
                            unordered_map<id_t, ChainDistances*>::iterator 
                                 chainDists = di->cd.find( get_start_of(
                                     *currChain).node_id());
                            if (chainDists != di->cd.end()) {
                                //Length of chain has already been found
                                nodeLen =chainDists->second->chainLength();
                                //last element should be length of chain

                            } else {//haven't recursed on this chain yet
                                nodeLen = calculateIndex( 
                                                      di, graph, sm, currChain);
                            }
                        } else {//Snarl
                            unordered_map<id_t, SnarlDistances*>::iterator 
                                 snarlDists = di->sd.find(currID.first);
                            if (snarlDists != di->sd.end()) {//Already found
                                nodeLen = snarlDists->second->snarlDistance(
                                  make_pair(currSnarl->start().node_id(), true),
                                  make_pair(currSnarl->end().node_id(), true))
                                +
                                  graph->get_node(currSnarl->end().node_id()
                                                       )->sequence().size();

                            } else {//Haven't recursed on snarl yet
                                Chain currChain;
                                currChain.push_back(currSnarl);
                                nodeLen = calculateIndex( 
                                                 di, graph, sm, &currChain);
                            }
                        }
                    } else { //Node is just a node
                        Node* n = graph->get_node(currID.first);
                        nodeLen = n->sequence().size();
                        //TODO NetGraph.get_length doesn't work?
                    }
                    if (ng.get_id(currHandle) != snarl->end().node_id() &&
                        !(ng.get_id(currHandle) == snarl->start().node_id()
                         && currID.second == true)) {
                    //Make sure that the next nodes are within the same snarl???
                    //If curr node is not end, add next nodes to priority queue
                        vector<const handle_t*> nextNodes;//vector of next nodes
                        auto addHandle = [&](const handle_t& h)-> bool {
                             nextNodes.push_back(&h);
                             return true;
                        };
                        ng.follow_edges(currHandle, false, addHandle);
                        for (int i = 0; i < nextNodes.size(); i++) {
                            const handle_t* h = nextNodes[i]; 
                            pair<const handle_t*, int64_t> p = 
                                            make_pair (h, currDist + nodeLen);
                            reachable.push(p);
                        }
                    }
                } 
            }//End while loop
        }}//End for loop over starting node/directions in a snarl
        /* add length of snarl (start of start node to start of end node) to
           the chain prefix sum vector */

        int64_t dist = sd.snarlDistance(
                  make_pair (snarl->start().node_id(), true), 
                  make_pair (snarl->end().node_id(), true) );
        chain_prefix_sum.push_back(chain_prefix_sum.back() + dist);

    }//End for loop over snarls in chain
    //Add the length of the last node in chain to get length of entire chain
    Visit lastVisit = get_end_of(*chain);
    //TODO Visit points <-, snarl as well??
    Node* lastNode = graph->get_node(lastVisit.node_id());
    
    chain_prefix_sum.push_back(chain_prefix_sum.back() + 
                                                  lastNode->sequence().size());

    if (chain_prefix_sum.size() > 2) { //If chain and not just one snarl
        ChainDistances cd = ChainDistances(&chain_snarl_to_index, 
                                                             &chain_prefix_sum);
        pair <id_t, ChainDistances*> p = 
                                 make_pair(get_start_of(*chain).node_id(), &cd);
        di->cd.insert(p);
    }
    return chain_prefix_sum.back();//return length of entire chain
};

DistanceIndex makeDistanceIndex(VG* graph, SnarlManager* sm) {
    /*Wrapper for creating the distance index given a VG and snarl manager
    */
    DistanceIndex distances;
    const vector<const Snarl*> topSnarls = sm->top_level_snarls();
    calculateIndex(&distances, graph, sm, &topSnarls);
    return distances;
};


//////////////////    Calculate distances

int64_t minPos (vector<int64_t> vals) {
    /*return the minimum value in vals that is not -1, returns -1 if all
     values are -1 */
    return accumulate(vals.begin(), vals.end(), -1, 
          [](int x, int y) {if (x==-1) {return y;} 
                            else if (y == -1) {return x;}
                            else {return min(x, y);}} 
          ); 
   
};

pair<int64_t, int64_t> minEndDists (int64_t dsl, int64_t dsr, int64_t del, 
      int64_t der, int64_t distL, int64_t distR) {
    /*Given the distances from a starting point to the boundary nodes of a 
     ancestor snarl and the distances from the bounary nodes of the ancestor 
     snarl to the boundary nodes of its parent snarl, find the minimum
     distances to end nodes not including -1 distances
     */
            dsl = dsl == -1 || distL == -1? -1 : distL + dsl; 
            dsr =  dsr == -1 || distR == -1? -1 : distR + dsr; 
            der = der == -1 || distR == -1? -1 : distR + der; 
            del = del == -1 || distL == -1? -1 : distL + del; 
            int64_t left[] = {dsl, dsr};
            int64_t right[] = {del, der}; 
            vector<int64_t> leftV (left, left+sizeof(left)/sizeof(int64_t));
            vector<int64_t> rightV (right, right+sizeof(right)/sizeof(int64_t));
            distL = minPos(leftV);
            distR = minPos(rightV);
            return make_pair(distL, distR);
};


pair<pair<int64_t, int64_t>, pair<const Snarl*, id_t>> distToCommonAncestor(
           VG* graph, DistanceIndex* di, SnarlManager* sm, const Snarl* snarl, 
           const Snarl* commonAncestor, pos_t pos, id_t node){

    /* Find the distance from pos to both boundary nodes of commonAncestor.
       Return the two distances and the node_id and Snarl of the ancestor snarl
       whose parent is the commonAncestor
    */
    int64_t distL; 
    int64_t distR;//Dist from pos1 to boundaries of curr snarl 
    if (is_rev(pos)) {
        distL = graph->get_node(get_id(pos))->sequence().size() -
                           get_offset(pos);
        distR = get_offset(pos);
    } else {
        distR = graph->get_node(get_id(pos))->sequence().size() -
                                                           get_offset(pos);
        distL = get_offset(pos);
    }
 
    id_t startID = snarl->start().node_id(); 
    id_t endID = snarl->end().node_id();
    unordered_map<id_t, SnarlDistances*>::iterator 
                         snarlDists = di->sd.find(snarl->start().node_id());

    int64_t dsl = snarlDists->second->snarlDistance( 
                             make_pair(startID, false), make_pair(node, false)); 
    int64_t dsr = snarlDists->second->snarlDistance( 
                          make_pair(startID, false), make_pair(node, true));

    int64_t der = snarlDists->second->snarlDistance(  
                        make_pair(endID, true), make_pair(node, true));

    int64_t del = snarlDists->second->snarlDistance( 
                        make_pair(endID, true), make_pair(node, false));

    pair<int64_t, int64_t> endDists = minEndDists(dsl, dsr, del,
                     der, distL, distR);    
    distL = endDists.first;
    distR = endDists.second;
    while (sm->parent_of(snarl) != commonAncestor) {
        int64_t dsl; int64_t dsr; int64_t der; int64_t del;
        if (sm->in_nontrivial_chain(snarl)) {
            const Chain* chain = sm->chain_of(snarl);
            id_t chainStartID = get_start_of(*chain).node_id();
            id_t chainEndID = get_end_of(*chain).node_id();
            unordered_map<id_t, ChainDistances*>::iterator chainDists = 
                          di->cd.find( chainStartID);
                
            dsl = chainDists->second->chainDistance(make_pair(
                     chainStartID, false), make_pair(node, false));
            der = chainDists->second->chainDistance(make_pair(
                        chainStartID, false), make_pair(node, true));
            der = chainDists->second->chainDistance(make_pair(
                         chainEndID, false), make_pair(node, true));
            del = chainDists->second->chainDistance(make_pair(
                          chainEndID, false), make_pair(node, false));
            pair<int64_t, int64_t> endDists = minEndDists(dsl, dsr, del,
                           der, distL, distR);    
            distL = endDists.first;   
            distR = endDists.second;
            node = chainStartID;
        }
        snarl = sm->parent_of(snarl);
        id_t startNodeID = snarl->start().node_id();
        id_t endNodeID = snarl->end().node_id();
        unordered_map<id_t, SnarlDistances*>::iterator snarlDists = 
                                      di->sd.find(snarl->start().node_id());
        dsl = snarlDists->second->snarlDistance(
                      make_pair(startNodeID, false), make_pair(node, false));
        dsr = snarlDists->second->snarlDistance(
                         make_pair(startNodeID, false), make_pair(node, true));
        der = snarlDists->second->snarlDistance(
                          make_pair(endNodeID, true), make_pair(node, true));
        del = snarlDists->second->snarlDistance(
                          make_pair(endNodeID, true), make_pair(node, false));
            
       pair<int64_t, int64_t> endDists = minEndDists(dsl, dsr, del,	
                                                            der, distL, distR);
          
       distL = endDists.first;
       distR = endDists.second;
       node = startNodeID;
     }
     return make_pair(make_pair(distL, distR), make_pair(snarl, node));
};

int64_t distance(VG* graph, SnarlManager* sm, DistanceIndex* di, 
                                               pos_t pos1, pos_t pos2) {
    if (get_id(pos1) == get_id(pos2)) { //if positions are on the same node
        int64_t offset1;
        if (is_rev(pos1)) {
            offset1 = graph->get_node(get_id(pos1))->sequence().size() -
                      get_offset(pos1);//Len of node - offset 
        } else {
            offset1 = get_offset(pos1);
        }
        int64_t offset2;
        if (is_rev(pos2)) {
            offset2 = graph->get_node(get_id(pos2))->sequence().size() - 
                      get_offset(pos2);
        } else {
            offset2 = get_offset(pos2);
        }
        return abs(offset1-offset2);

    } else { //Positions are on different nodes
        id_t node1 = get_id(pos1);
        id_t node2 = get_id(pos2); 
        const Snarl* snarl1 = sm->into_which_snarl(node1, is_rev(pos1));
        const Snarl* snarl2 = sm->into_which_snarl(node2, is_rev(pos2));
        const Snarl* commonAncestor; 

        //// Find common ancestor
        unordered_set<const Snarl*> ancestors1;
                      //set of all ancestor snarls of node1
        const Snarl* ancestor1 = sm->into_which_snarl(node1, is_rev(pos1));
        while (ancestor1 != NULL) {
            ancestors1.emplace(ancestor1);
            ancestor1 = sm->parent_of(ancestor1);
        }

        const Snarl* ancestor2 = sm->into_which_snarl(node2, is_rev(pos2));
        while (ancestor2 != NULL) {

            if (ancestors1.count(ancestor2) > 0) { 
                commonAncestor = ancestor2;
                break;
            }
        } 

        //Find distances from pos1 and pos2 to ends of child snarls of ancestor
        pair<pair<int64_t, int64_t>, pair<const Snarl*, id_t>> p1 = 
        distToCommonAncestor(graph, di, sm, snarl1, commonAncestor, pos1,node1);
        pair<int64_t, int64_t> d1 = p1.first; 
        pair<const Snarl*, id_t> d2 = p1.second;
        int64_t distL1 = d1.first; int64_t distR1 = d1.second;
        snarl1 = d2.first; node1 = d2.second;

        pair<pair<int64_t, int64_t>, pair<const Snarl*, id_t>> p2 = 
        distToCommonAncestor(graph, di, sm, snarl2, commonAncestor, pos2,node2);
        pair<int64_t, int64_t> d3 = p2.first; 
        pair<const Snarl*, id_t> d4 = p2.second;
        int64_t distL2 = d3.first; int64_t distR2 = d3.second;
        snarl2 = d4.first; node2 = d4.second;
        //Snarl1 and snarl2 are children of common ancestor

        int64_t chainDist = -1; 
        //Find shortest distance between boundary nodes of snarls containing pos
        if (sm->in_nontrivial_chain(snarl1) && sm->in_nontrivial_chain(snarl2)
               && sm->chain_of(snarl1) == sm->chain_of(snarl2)) {
            const Chain* chain = sm->chain_of(snarl1);
            id_t chainStartID = get_start_of(*chain).node_id();
            unordered_map<id_t, ChainDistances*>::iterator chainDists = 
                  di->cd.find( chainStartID); 
            int64_t d1 = chainDists->second->chainDistance_short(graph,
                   make_pair(node1, false), make_pair(node2, false));
            d1 = (distL1 == -1 || distL2 == -1 || d1 == -1) ? -1 : 
                                    distL1 + distL2 + d1; 
            int64_t d2 = chainDists->second->chainDistance_short(graph,
                   make_pair(node1, false), make_pair(node2, true));
            d2 = (distL1 == -1 || distR2 == -1 || d2 == -1) ? -1 : 
                                       distL1 + distR2 + d2; 
            int64_t d3 = chainDists->second->chainDistance_short(graph,
                   make_pair(node1, true), make_pair(node2, false));
            d3 = (distR1 == -1 || distL2 == -1 || d3 == -1) ? -1 : 
                                       distR1 + distL2 + d3; 
            int64_t d4 =  chainDists->second->chainDistance_short(graph,
                   make_pair(node1, true), make_pair(node2, true));
            d4 = (distR1 == -1 || distR2 == -1 || d4 == -1) ? -1 : 
                                       distR1 + distR2 + d4;
           int64_t vals[] {d1, d2, d3, d4}; 
           vector<int64_t> valsv (vals, vals+sizeof(vals)/sizeof(int64_t));
           chainDist = minPos(valsv);

        } else { //if snarls are not in the same chain

            //Get dist from pos1 to ends of chain 
            if (sm->in_nontrivial_chain(snarl1)) {
                const Chain* chain = sm->chain_of(snarl1);
                id_t chainStartID = 
                                get_start_of(*chain).node_id();
                id_t chainEndID = 
                                get_end_of(*chain).node_id();
                unordered_map<id_t, ChainDistances*>::iterator chainDists = 
                  di->cd.find( chainStartID);

                int64_t dsl = chainDists->second->chainDistance(make_pair(
                               chainStartID, false), make_pair(node1, false));
                int64_t dsr = chainDists->second->chainDistance(make_pair(
                               chainStartID, false), make_pair(node1, true));
                int64_t der = chainDists->second->chainDistance(make_pair(
                               chainEndID, false), make_pair(node1, true));
                int64_t del = chainDists->second->chainDistance(make_pair(
                               chainEndID, false), make_pair(node1, false));

                pair<int64_t, int64_t> endDists = minEndDists(dsl, dsr, del,
                                                      der, distL1, distR1);    
                distL1 = endDists.first;
                distR1 = endDists.second;
                node1 = chainStartID;
            }
            //Get dist from pos2 to ends of its chain 
            if (sm->in_nontrivial_chain(snarl2)) {
                const Chain* chain = sm->chain_of(snarl2);
                id_t chainStartID = 
                                get_start_of(*chain).node_id();
                id_t chainEndID = 
                                get_end_of(*chain).node_id();
                unordered_map<id_t, ChainDistances*>::iterator chainDists = 
                  di->cd.find( chainStartID);

                int64_t dsl =  chainDists->second->chainDistance(make_pair(
                               chainStartID, false), make_pair(node2, false));
                int64_t dsr = chainDists->second->chainDistance(make_pair(
                               chainStartID, false), make_pair(node2, true));
                int64_t der = chainDists->second->chainDistance(make_pair(
                               chainEndID, false), make_pair(node2, true));
                int64_t del = chainDists->second->chainDistance(make_pair(
                               chainEndID, false), make_pair(node2, false));
                pair<int64_t, int64_t> endDists = minEndDists(dsl, dsr, del,
                                                      der, distL2, distR2);    
                distL2 = endDists.first;
                distR2 = endDists.second;

                node2 = chainStartID;
            }
            //Get distance between end nodes in snarl
            unordered_map<id_t, SnarlDistances*>::iterator 
                 snarlDists = di->sd.find(commonAncestor->start().node_id());
            int64_t d1 = snarlDists->second->snarlDistance_short(graph,
                     make_pair(node1, false), make_pair(node2, false));
            d1 = (distL1 == -1 || distL2 == -1 || d1 == -1) ? -1 : 
                                                   distL1 + distL2 + d1; 
            int64_t d2 = snarlDists->second->snarlDistance_short(graph,
                     make_pair(node1, false), make_pair(node2, true));

            d2 = (distL1 == -1 || distR1 == -1 || d2 == -1) ? -1 :
                                                   distL1 + distR2 + d2;
            int64_t d3 = snarlDists->second->snarlDistance_short(graph,
                     make_pair(node1, true), make_pair(node2, false));
            d3 = (distR1 == -1 || distL2 == -1 || d3 == -1) ? -1 : 
                                                  distR1 + distL2 + d3; 
            int64_t d4 = snarlDists->second->snarlDistance_short(graph,
                     make_pair(node1, true), make_pair(node2, true));
            d4 = (distR1 == -1 || distR2 == -1 || d4 == -1) ? -1 : 
                                                  distR1 + distR2 + d4; 
            vector<int64_t> vals {d1, d2, d3, d4, chainDist}; 
            return minPos(vals);//TODO -1
            
        }
    }
};
}
