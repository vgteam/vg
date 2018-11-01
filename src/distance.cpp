//#define indexTraverse
//#define printDistances

#include "distance.hpp"

using namespace std;
namespace vg {

/*TODO: 
 * Also change how nodes are stored in chain - in case of loops/unary snarls -might not actually need this
 * Make snarls/chains represented by the node id in netgraph
 */
DistanceIndex::DistanceIndex(HandleGraph* vg, SnarlManager* snarlManager, uint64_t cap){
    /*Constructor for the distance index given a VG and snarl manager
      cap is the largest distance that the maximum distance estimation will be
      accurate to 
    */

    if (snarlManager->top_level_snarls().size() == 0 && maxNodeID > 1) {
 
        throw runtime_error("Snarl manager is empty");       
    }
    minNodeID = -1;
    maxNodeID = -1;


    graph = vg;
    sm = snarlManager;

    #ifdef indexTraverse
        cerr << endl << "Creating distance index"<< endl;
    #endif

    //Calculate minimum distance index
    const vector<const Snarl*> topSnarls = sm->top_level_snarls();

    unordered_set<const Snarl*> seenSnarls;
    for (const Snarl* snarl : topSnarls) {
       //Make an index for each disconnected snarl/chain
 
       if (seenSnarls.count(snarl) == 0){
          if (sm->in_nontrivial_chain(snarl)){
              const Chain* chain = sm->chain_of(snarl);
              calculateMinIndex(chain);
              for (auto s : *chain) {
                  seenSnarls.insert(s.first);
              }
           } else {
               Chain currChain;
               currChain.emplace_back(snarl, false);
               calculateMinIndex(&currChain);
               seenSnarls.insert(snarl);
           }
           
       }
        
    }
    nodeToSnarl = calculateNodeToSnarl(snarlManager);
    //TODO: Cap should be given
    maxIndex = MaxDistanceIndex (this, topSnarls, cap);
  
};


DistanceIndex::DistanceIndex(HandleGraph* vg, SnarlManager* snarlManager, istream& in) {

    /*Constructor for the distance index given a VG, snarl manager, and a vector
      of ints from serialization
    */
   
    graph = vg;
    sm = snarlManager;
    load(in);

    for ( auto x : snarlDistances ) {
    //Check that vg and snarl manager match the distance index
  
        pair<id_t, bool> node = x.first;
/* TODO: Make test that uses handle graph
        if (!graph->has_node(node.first)) {

            throw runtime_error("Distance index does not match vg");

        } else */if (sm->into_which_snarl(node.first, node.second) == NULL) {

            throw runtime_error("Distance index does not match snarl manager");
        }

    }

}
void DistanceIndex::load(istream& in){
    //Load serialized index from an istream
    auto toInt = [] (uint64_t uval) {
        /*convert unsigned representation of signed int back to int64_t*/
        int64_t val = uval / 2;
        if (uval % 2 == 1) {val = -val;}
        return val;
    };

    int_vector<> d1;
    int_vector<> d2;
    int_vector<> d3;
    int_vector<> d4;
  
    d1.load(in);
    d2.load(in);
    d3.load(in);
    d4.load(in);
    nodeToSnarl.load(in);
    maxIndex.nodeToComponent.load(in);
    maxIndex.minDistances.load(in);
    maxIndex.maxDistances.load(in);

    vector<int64_t> snarlNodes(d1.size()-4, 0); 
    vector<int64_t> snarlVector(d2.size(), 0);
    vector<int64_t> chainNodes(d3.size(), 0);
    vector<int64_t> chainVector(d4.size(), 0);

    
 
    minNodeID = d1[0];
    maxNodeID = d1[1];
    maxIndex.numCycles = d1[2];
    maxIndex.cap = d1[3];
    

    for (size_t i = 4; i < d1.size(); i++) {
        uint64_t uval = d1[i];
        int64_t val = toInt(uval);
        snarlNodes[i-4] = val;
    }

    for (size_t i = 0; i < d2.size(); i++) {
        snarlVector[i] = toInt(d2[i]);
    }
    for (size_t i = 0; i < d3.size(); i++) {
        chainNodes[i] = toInt(d3[i]);
    }
    for (size_t i = 0; i < d4.size(); i++) {
        chainVector[i] = toInt(d4[i]);
    }

  
    //Construct snarl index
    size_t snarlI = 0;//Index into snarlVector
    for (size_t i = 0; i < snarlNodes.size()/2; i++) {
        int64_t snarlInt = snarlNodes[2*i];
  
        pair<id_t, bool> node = snarlInt < 0 ? 
                                     make_pair( (id_t) abs(snarlInt), true) :
                                     make_pair( (id_t) abs(snarlInt), false);
        size_t nextIndex = snarlI + snarlNodes[2*i+1];


        vector<int64_t> snarlv;//get subvector
        snarlv.resize(nextIndex - snarlI);
        for (size_t j = 0; j < nextIndex-snarlI; j++) {
            snarlv[j] = snarlVector[snarlI + j]; 
        }
         

        //Create SnarlIndex object from vector and assign
        snarlDistances.insert(make_pair(node, SnarlIndex (this, snarlv)));

        snarlI = nextIndex; 
      
    }
    
    
    size_t chainI = 0; //Index into chainVector
    //Construct chain index
    for (size_t i = 0; i < chainNodes.size()/2; i++) {
        id_t chainID = (id_t) chainNodes[2*i];
        
        size_t nextIndex = chainI + chainNodes[2*i + 1];
   
        vector<int64_t> chainv;
        chainv.resize(nextIndex - chainI);
        for ( size_t j = 0; j < nextIndex - chainI; j++) {

            chainv[j] = chainVector[chainI + j];
        }

        //Create chaindistances object and assign in index
        chainDistances.insert(make_pair(chainID, ChainIndex(this, chainv)));

        chainI = nextIndex;

    }
    maxIndex.distIndex = this;
};

void DistanceIndex::serialize(ostream& out) {

    auto toUint = [](int64_t val) {
        /* convert signed integer into unsigned representation where last bit 
           represents sign*/ 
        uint64_t uval= abs(val) * 2;
        if (val < 0) { uval += 1; }
        return uval;
    };

    vector<int64_t> d1; //Snarl nodes as a vector [node, length, node, ...]
    vector<int64_t> d2; //Snarl distances as a vector
    vector<int64_t> d3; //Chain nodes as a vector
    vector<int64_t> d4; //Chain distances as a vector

    size_t snarlNodesI = 0;
    size_t snarlVectorI = 0;
    //Serialize snarls
    d1.resize(2*snarlDistances.size());

    for (pair<pair<id_t, bool>, SnarlIndex> snarlPair : snarlDistances) {
        int64_t nodeInt = snarlPair.first.second ? 
            - (int64_t) snarlPair.first.first : (int64_t) snarlPair.first.first;
        vector<int64_t> currVector = snarlPair.second.toVector();
        
        d1[snarlNodesI++] = nodeInt;
        d1[snarlNodesI++] = (int64_t) currVector.size();
    
        d2.resize(d2.size() + currVector.size());
        //Concatenate new vector onto whole snarl vector
        for (auto x : currVector) {
            d2[snarlVectorI++] = x;
        }
    }

    size_t chainNodesI = 0;
    size_t chainVectorI = 0;
    //Serialize chains
    for (pair<id_t, ChainIndex> chainP: chainDistances) {
        vector<int64_t> currVector = chainP.second.toVector();
        
        d3.resize(d3.size() + 2 );
        d3[chainNodesI++] = (int64_t) chainP.first;
        d3[chainNodesI++] = (int64_t) currVector.size();
        
        d4.resize(d4.size() + currVector.size()); 
        for (auto x : currVector) {
            d4[chainVectorI ++ ] = x;
        } 
    }
    
    //Convert vectors of ints to int_vector TODO: Start with int_vector
    int_vector<> snarlNodes;
    int_vector<> snarlVector;
    int_vector<> chainNodes;
    int_vector<> chainVector;

    util::assign(snarlNodes, int_vector<>(d1.size() + 4));
    util::assign(snarlVector, int_vector<>(d2.size()));
    util::assign(chainNodes, int_vector<>(d3.size()));
    util::assign(chainVector, int_vector<>(d4.size()));

    snarlNodes[0] = minNodeID;
    snarlNodes[1] = maxNodeID;
    snarlNodes[2] = maxIndex.numCycles;
    snarlNodes[3] = maxIndex.cap;
    //Copy vector<int64_t> into int_vector
    for (size_t i = 0; i < d1.size(); i++) {
        snarlNodes[i+4] = toUint(d1[i]);
    }

    for (size_t i = 0; i < d2.size(); i++) {
        snarlVector[i] = toUint(d2[i]);
    }
    for (size_t i = 0; i < d3.size(); i++) {
        chainNodes[i] = toUint(d3[i]);
    }
    for (size_t i = 0; i < d4.size(); i++) {
        chainVector[i] = toUint(d4[i]);
    }
        
    util::bit_compress(snarlNodes);
    util::bit_compress(snarlVector);
    util::bit_compress(chainNodes);
    util::bit_compress(chainVector);

    //Serialize
    
    snarlNodes.serialize(out, NULL, "snarl_nodes");
    snarlVector.serialize(out, NULL, "snarl_vector");
    chainNodes.serialize(out, NULL, "chain_nodes");
    chainVector.serialize(out, NULL, "chain_vector");
    nodeToSnarl.serialize(out, NULL, "node_to_snarl");
    maxIndex.nodeToComponent.serialize(out, NULL, "node_to_comp");
    maxIndex.minDistances.serialize(out, NULL, "minDists");
    maxIndex.maxDistances.serialize(out, NULL, "maxDists");
}


int_vector<> DistanceIndex::calculateNodeToSnarl(SnarlManager* sm){

    auto toUint = [](int64_t val) {
        /* convert signed integer into unsigned representation where last bit 
           represents sign*/ 
        uint64_t uval= abs(val) * 2;
        if (val < 0) { uval += 1; }
        return uval;
    };

    int_vector<> result(maxNodeID - minNodeID + 1, 0);

    const vector<const Snarl*> topSnarls = sm->top_level_snarls();
    vector<const Snarl*> allSnarls(topSnarls.begin(), topSnarls.end());

    while (allSnarls.size() > 0) {
        const Snarl* snarl = allSnarls.back();
        allSnarls.pop_back();

        int64_t currSnarlID = snarl->start().backward() ? 
                         -snarl->start().node_id() : snarl->start().node_id();

        NetGraph ng = NetGraph(snarl->start(), snarl->end(), sm->chains_of(snarl), graph);

        //Get all the nodes in the snarl
        
        vector<id_t> allNodes;
   
        auto addNode = [&](const handle_t& h)-> bool {
            allNodes.push_back(ng.get_id(h));
            return true;
                   
        };
        ng.for_each_handle(addNode);

        for ( id_t nodeID : allNodes) {
            
            const Snarl* tempSnarl = sm->into_which_snarl(nodeID, true); 
            const Snarl* nextSnarl = tempSnarl == NULL ? 
                                sm->into_which_snarl(nodeID, false) : tempSnarl;

            if (nodeID != snarl->start().node_id() && 
                        nodeID != snarl->end().node_id() && nextSnarl != NULL) {
                //If this node is a child snarl
                if (sm->in_nontrivial_chain(nextSnarl)) {
                    const Chain* chain = sm->chain_of(nextSnarl);
                    for (auto s : *chain) {
                        allSnarls.push_back(s.first);
                    }
                } else {
                    allSnarls.push_back(nextSnarl);
                }

            } else {

                //If this node is just a node
                result[nodeID - minNodeID] = toUint(currSnarlID); 

            }
        }
    }
    util::bit_compress(result);
    return result;

}



/////////////////////////    MINIMUM INDEX    ///////////////////////////////



int64_t DistanceIndex::calculateMinIndex(const Chain* chain) {
    /*Populate the DistanceIndex
      Add chain to distance index and recursively add all distances in 
      snarls contained within chain
    */
    auto cmp = [] (pair<pair<id_t, bool>,int64_t> x,
                                           pair<pair<id_t, bool>,int64_t> y) {
        //Comparison function for the priority of a pair of handle, distance
        return (x.second > y.second);
    };

    //Vectors for chain distance index
    //initialize chain prefix sum to [0, len of first node in chain]
    vector<int64_t> chainPrefixSum (1,0); 
    vector<int64_t> chainLoopFd; 
    vector<int64_t> chainLoopRev; 
    handle_t firstNode = graph->get_handle(get_start_of(*chain));
    
    
    hash_map<id_t, size_t> snarlToIndex;
    snarlToIndex[get_start_of(*chain).node_id()] = 0;
    #ifdef indexTraverse
        cerr << "Node " << get_start_of(*chain).node_id() << " represents snarl at index "
            << snarlToIndex[get_start_of(*chain).node_id()] << endl;
    #endif

    ChainIterator chainEnd = chain_end(*chain);
    for (ChainIterator c = chain_begin(*chain); c != chainEnd; ++c) {
        //for each snarl in the chain 


        const Snarl* snarl = c->first;
        bool snarlRevInChain = c->second;

        id_t snarlStartID = snarl->start().node_id();
        bool snarlStartRev = snarl->start().backward(); //into snarl
        id_t snarlEndID = snarl->end().node_id();
        bool snarlEndRev = snarl->end().backward();   //pointing out


/*TODO: Make a test that uses handle graph
        if (!( graph->has_node(snarlStartID) && graph->has_node(snarlEndID))) {
            //Make sure that vg contains the boundary nodes of this snarl
            throw runtime_error("Snarl manager does not match vg");
        }
*/

        if (snarlToIndex.find(snarlEndID) == snarlToIndex.end()){
            //Store the index of the start of the snarl only if it hasn't
            //already been seen (if the chain loops)
            size_t nextIndex = snarlToIndex.size();
            snarlToIndex[snarlEndID] = nextIndex;
            
            #ifdef indexTraverse
                cerr << "Node " << snarlEndID << " represents snarl at index "
                    << snarlToIndex[snarlEndID] << endl;
            #endif
        } else {
            #ifdef indexTraverse
                cerr << "Node " << snarlEndID << " already represents a snarl, at index " << snarlToIndex[snarlEndID] << endl;
            #endif
        }


        NetGraph ng = 
            NetGraph(snarl->start(), snarl->end(), sm->chains_of(snarl), graph);



        //Get all the nodes in the snarl
        
        unordered_set<pair<id_t, bool>> allNodes;
   
        auto addNode = [&](const handle_t& h)-> bool {
            allNodes.insert(make_pair(ng.get_id(h), ng.get_is_reverse(h)));
            allNodes.insert(make_pair(ng.get_id(h), !ng.get_is_reverse(h)));
            return true;
                   
        };
        ng.for_each_handle(addNode);

        //Create snarl distance obj for current snarl and add to distance index
        snarlDistances.insert(make_pair( make_pair(snarlStartID,snarlStartRev),
            SnarlIndex(this,allNodes, make_pair(snarlStartID,snarlStartRev),
                 make_pair(snarlEndID,  snarlEndRev))));

        SnarlIndex& sd = snarlDistances.at(make_pair(snarlStartID, snarlStartRev));

        #ifdef indexTraverse
            cerr << "Snarl at " << snarl->start() << " -> " << snarl->end() << endl;
            cerr << "    Contains nodes : ";
            {
                unordered_set<id_t> reported;
                for (pair<id_t, bool> node: allNodes) {
                    if (!reported.count(node.first)) {
                        cerr << node.first << " ";
                        reported.insert(node.first);
                    }
                }
            }
            cerr << endl;
        #endif

        for (pair<id_t, bool> startID : allNodes){
            //Use each node in the snarl as start of djikstra search


            minNodeID = minNodeID == -1 ? startID.first : 
                                   min(minNodeID, startID.first);
            maxNodeID = max(maxNodeID, startID.first);

            handle_t startHandle = 
                           graph->get_handle(startID.first, startID.second);
            //Priority queue of reachable nodes (pair of node id and direction)
            priority_queue<  pair<pair<id_t, bool>, int64_t>,  
                         vector<pair<pair<id_t, bool>, int64_t>>, 
                                     decltype(cmp)> reachable(cmp);
            reachable.push(make_pair(startID, 0));

            #ifdef indexTraverse
                cerr << "  Start Node: " << startID.first << "," 
                    << startID.second << endl;
            #endif
            bool firstLoop = true;
            unordered_set<pair<id_t, bool>> seenNodes;

            while (reachable.size() > 0) {
                pair<pair<id_t, bool>, int64_t> next = reachable.top();
                reachable.pop();
                pair<id_t, bool> currID= next.first;
                handle_t currHandle = graph->get_handle(currID.first, 
                                                        currID.second);
                int64_t currDist = next.second;
                if ( seenNodes.count(currID) == 0) {
                    //If node has not already been found:

                    //Save distance from start to current node 
                    if (!firstLoop) {

                        sd.insertDistance(startID, currID, currDist);
                        seenNodes.insert(currID);

                    }
               

                    
                    int64_t nodeLen; //length of the current node
                       
                    int64_t loopDist = -1;
                         //Dist to enter curr node then exit at same side 

                    //Get the snarl that the node represents, if any
                    const Snarl* tempSnarl = sm->into_which_snarl(
                                              currID.first, currID.second);
                    const Snarl* currSnarl = tempSnarl == NULL ? 
                        sm->into_which_snarl(currID.first, !currID.second) :
                        tempSnarl; 

                    if (currID.first != snarlStartID &&
                            currID.first != snarlEndID &&
                            currSnarl != NULL) {
                        //If current node is a child snarl/chain


                        if (sm->in_nontrivial_chain(currSnarl)) {
                           //The node is a chain

                            const Chain* currChain= sm->chain_of(currSnarl);
                            auto chainDists = chainDistances.find( get_start_of(
                                     *currChain).node_id());

                            if (chainDists != chainDistances.end()) {
                                //Length of chain has already been found
                               
                                //Get the length of the node (chain)
                                nodeLen = chainDists->second.chainLength();
       
                                //Get loop dist- enter and exit chain at same side
                                if (get_start_of(
                                     *currChain).backward() == currID.second) {
                                    //If traversing snarl forward in chain

                                    loopDist = chainDists->second.loopFd[0] - 1;

                                    if (loopDist != -1) {
                                       loopDist = loopDist  + graph->get_length(
                                                                  currHandle) ;
                                     }

                                } else {
                                    loopDist = chainDists->second.loopRev[
                                      chainDists->second.loopRev.size()-1] - 1;
                                      handle_t tempHandle = graph->get_handle(
                                                        get_end_of(*currChain));

                                    if (loopDist != -1) { 
                                        loopDist = loopDist + 
                                           graph->get_length(tempHandle);
                                     }
                                }

                            } else {//haven't recursed on this chain yet
                                #ifdef indexTraverse
                                    cerr << " recurse" << endl;
                                #endif
                                nodeLen = calculateMinIndex(currChain);

                                ChainIndex& currChainDists = 
                                         chainDistances.at( get_start_of(
                                                    *currChain).node_id());
                                if (get_start_of(
                                     *currChain).backward() == currID.second) {
                                    //If traversing snarl forward in chain

                                    loopDist = currChainDists.loopFd[0] - 1;

                                    if (loopDist != -1) {
                                        loopDist = loopDist + graph->get_length(
                                                                   currHandle);
                                    }
                                } else {

                                    loopDist = currChainDists.loopRev[
                                           currChainDists.loopRev.size()-1] - 1;

                                    if (loopDist != -1) { 
                                        handle_t tempHandle = graph->get_handle(
                                                     get_end_of(*currChain));
                                        loopDist = loopDist +
                                           graph->get_length(tempHandle);
                                     }
                                } 
                            }
                        } else {//Snarl

                            id_t snarlID = currSnarl->start().node_id();
                            bool snarlRev = currSnarl->start().backward();
                            id_t endID = currSnarl->end().node_id();
                            bool endRev = currSnarl->end().backward();
  
                            auto snarlDists = snarlDistances.find(make_pair(
                                        snarlID, snarlRev));

                            if (snarlDists != snarlDistances.end()) {//Already found
                                nodeLen = snarlDists->second.snarlLength(
                                                                    graph, &ng);

                                //Find the distance to enter and exit snarl
                                //at the same side
                                if (currID.second == snarlRev) { 
                                    //If traversing snarl forward
                                    loopDist = snarlDists->second.snarlDistance(
                                          graph, &ng, 
                                           make_pair(snarlID, snarlRev), 
                                           make_pair(snarlID, !snarlRev));

                                     if (loopDist != -1) { 
                                         loopDist = loopDist+ graph->get_length(
                                                                    currHandle);
                                     }
                                } else {
                                    loopDist = snarlDists->second.snarlDistance(
                                           graph, &ng, 
                                           make_pair(endID, !endRev), 
                                           make_pair(endID, endRev));

                                     if (loopDist != -1) { 
                                         handle_t tempHandle = 
                                            graph->get_handle(currSnarl->end());
                                         loopDist = loopDist + 
                                            graph->get_length(tempHandle); 
                                     }
                                }
                            } else {//Haven't recursed on snarl yet
                                #ifdef indexTraverse
                                    cerr << " recurse" << endl;
                                #endif
                                
                                //Create chain to recurse on and recurse
                                Chain currChain;

                                currChain.emplace_back(currSnarl, false);
                                calculateMinIndex(&currChain);

                                SnarlIndex& currSnarlDists = 
                                     snarlDistances.at(make_pair(snarlID,snarlRev));

                                nodeLen = currSnarlDists.snarlLength(graph, &ng); 

                                //Find the distance to enter and exit snarl
                                //at the same side
                                if (currID.second == snarlRev) {

                                    loopDist = currSnarlDists.snarlDistance(
                                        graph, &ng, 
                                       make_pair(snarlID, snarlRev), 
                                       make_pair(snarlID, !snarlRev));

                                     if (loopDist != -1) { 
                                         loopDist = loopDist +graph->get_length(
                                                                 currHandle);
                                     }

                                 } else {

                                     loopDist = currSnarlDists.snarlDistance(
                                         graph, &ng, 
                                        make_pair(endID, !endRev), 
                                        make_pair(endID, endRev));

                                     if (loopDist != -1) { 
                                         handle_t tempHandle = graph->get_handle
                                                    (currSnarl->end());
                                         loopDist = loopDist+ 
                                                  graph->get_length(tempHandle);
                                      }
                                 }
                            }
                                        
                        }
                    } else { //Node is just a node
                        nodeLen = graph->get_length(currHandle); 
                    }
       

                    const handle_t currHandle = ng.get_handle(currID.first, 
                                                               currID.second);


                    if (loopDist != -1 && !firstLoop) {
                        /*If there is a path within the current node that loops 
                          to enter the node and exit it at the same side - add
                          reachable nodes from current node in reverse 
                          Do not add this distance if the current node is the 
                          starting node */

                        handle_t revHandle = ng.get_handle(
                                          ng.get_id(currHandle), 
                                          !ng.get_is_reverse(currHandle)); 
                            

                        auto addRevHandle = [&](const handle_t& h)-> bool {
                            pair<id_t, bool> node = make_pair(
                                        ng.get_id(h), ng.get_is_reverse(h));
                            reachable.push(make_pair(node, 
                                                     currDist + loopDist));
 

                             return true;
                        };

                        ng.follow_edges(revHandle, false, addRevHandle);
                    }

                    //Add reachable nodes to priority queue
                    auto addHandle = [&](const handle_t& h)-> bool {
                         pair<id_t, bool> node = make_pair(
                                    ng.get_id(h), ng.get_is_reverse(h));
                       if (nodeLen != -1) {
                       reachable.push(make_pair(node, currDist + nodeLen));
                       }
                      
                         #ifdef indexTraverse
                             cerr << node.first << " " << node.second << ", ";
                         #endif
                         return true;
                    };
                    //Add reachable nodes to priority queue for unary snarl that doesn't loop - 0 distance
                    auto addHandle0 = [&](const handle_t& h)-> bool {
                         pair<id_t, bool> node = make_pair(
                                    ng.get_id(h), ng.get_is_reverse(h));
                       reachable.push(make_pair(node, 0));
                       
                         #ifdef indexTraverse
                             cerr << node.first << " " << node.second << ", ";
                         #endif
                         return true;
                    };


                    if ( (nodeLen == -1 && firstLoop) || currID == startID) {
                        //If the nodeLen is -1 then node is a unary snarl that doesn't have a path from start to end. If this is the start of the distance calculation then add subsequent nodes assuming that the node length was 0
                        //Or if this is the starting node

                    #ifdef indexTraverse

                         cerr << "    From start node " << startID.first << " " << startID.second 
                            << " in snarl " << snarl->start() << " -> " << snarl->end()
                            << " at " << ng.get_id(currHandle) << " " << ng.get_is_reverse(currHandle) << endl; 
                         cerr << "        Adding next nodes:  ";
                    #endif
                        ng.follow_edges(currHandle, false, addHandle0);

                    } else  {

                    #ifdef indexTraverse
                         cerr << "    From start node " << startID.first << " " << startID.second<< " in snarl " << snarl->start().node_id() << " at " << ng.get_id(currHandle) << " " << ng.get_is_reverse(currHandle) << endl; 
                         cerr << "        Adding next nodes:  ";
                    #endif
                        ng.follow_edges(currHandle, false, addHandle);
                    }  

                    if (nodeLen != -1) {
                        ng.follow_edges(currHandle, false, addHandle);
                    } else if (firstLoop) {
                        //If the nodeLen is -1 then node is a unary snarl that
                        //doesn't have a path from start to end. If this is the
                        //start of the distance calculation then add subsequent
                        //nodes assuming that the node length was 0
                        ng.follow_edges(currHandle, false, addHandle0);
                    } 
                        

                    //Add edges between the boundary nodes that are not in 
                    //the net graph
                    int64_t nextDist = currID == startID ? 0 : currDist+nodeLen;

                    if ((currID.first == snarlStartID &&
                        currID.second != snarlStartRev) ||
                         ( currID.first == snarlEndID &&
                                 currID.second == snarlEndRev )  ) {
                           
                        //If currently leaving start of snarl
                        auto addHandleEnd = [&](const handle_t& h)-> bool {
                            pair<id_t, bool> node = make_pair(
                                    ng.get_id(h), ng.get_is_reverse(h));
                             if ( node.first == snarlStartID || 
                                  node.first == snarlEndID ) {
                               reachable.push(make_pair(node, nextDist));
                            }
                            return true;
                        };
                        graph->follow_edges(currHandle, false, addHandleEnd);

                    }                      
                    #ifdef indexTraverse
                         cerr << "    prev dist: " << currDist << "+ new dist " << nodeLen << endl;
                    #endif
                } 
                firstLoop = false;
            }//End while loop
        }//End for loop over starting node/directions in a snarl

        #ifdef indexTraverse
            cerr << "End snarl " << snarl->start() << " -> " << snarl->end() << endl;
        #endif

        /*Add to prefix sum the distance to the beginning and end of the last
            node in the current snarl
        */

        int64_t dist;
        if (snarlRevInChain) {
            //If traversing snarl backwards 
            dist = sd.snarlDistance( graph, &ng, 
                 make_pair (snarlEndID, !snarlEndRev) , 
                 make_pair (snarlStartID, !snarlStartRev));

            chainPrefixSum.push_back(chainPrefixSum.back()+ dist);

            #ifdef indexTraverse
                cerr << "Prefix sum before snarl reverse start: " << chainPrefixSum.back() << endl;
            #endif
        
        } else { 
            dist = sd.snarlDistance( graph, &ng,
                  make_pair (snarlStartID, snarlStartRev),
                  make_pair (snarlEndID, snarlEndRev) );

            chainPrefixSum.push_back(chainPrefixSum.back()+dist);
            #ifdef indexTraverse
                cerr << "Prefix sum before snarl end: " << chainPrefixSum.back() << endl;
            #endif
        }
        
        //Bit compress distance matrix of snarl index
        util::bit_compress(sd.distances);

    }//End for loop over snarls in chain


    //Add reverse loop distances


    for (ChainIterator c = chain_begin(*chain); c != chainEnd; ++c) {
        //Loop through the chain in reverse
        const Snarl* snarl = c->first; 
        bool snarlRevInChain = c->second;
        id_t snarlStartID = snarl->start().node_id();
        bool snarlStartRev = snarl->start().backward();
        id_t snarlEndID = snarl->end().node_id();
        bool snarlEndRev = snarl->end().backward();
        auto& sd =  snarlDistances.at(make_pair(snarlStartID, snarlStartRev)); 
        NetGraph ng (snarl->start(), snarl->end(),sm->chains_of(snarl), graph);
        //Add reverse loop distances- from start node rev to start node fd



        if ( chainLoopRev.size() == 0) {
           int64_t firstRevDist;
           if (snarlRevInChain){ 
                //If this is the first snarl in the chain 
                firstRevDist = sd.snarlDistance( graph, &ng,
                          make_pair(snarlEndID, snarlEndRev),
                          make_pair(snarlEndID, !snarlEndRev));
            } else {
                firstRevDist = sd.snarlDistance( graph, &ng,
                   make_pair (snarlStartID, !snarlStartRev),
                   make_pair (snarlStartID, snarlStartRev));
            }


            if (get_start_of(*chain).node_id() == 
                                                 get_end_of(*chain).node_id()) {
                //If the chain loops, might need distance from last snarl
                ChainIterator chainEndR = chain_rbegin(*chain);
                const Snarl* lastSnarl = chainEndR->first;
                bool lastRev = chainEndR->second;
                    
                id_t lastStartID = lastSnarl->start().node_id();
                bool lastStartRev = lastSnarl->start().backward();
                id_t lastEndID = lastSnarl->end().node_id();
                bool lastEndRev = lastSnarl->end().backward();
                auto& sdLast =  snarlDistances.at(make_pair(
                                              lastStartID, lastStartRev)); 

                if (lastRev) {
                    firstRevDist = minPos({firstRevDist, 
                             sdLast.snarlDistance(graph, &ng,
                                  make_pair(lastStartID, lastStartRev), 
                                   make_pair(lastStartID, !lastStartRev)) });
                  
              
                } else { 
                    firstRevDist = minPos({firstRevDist, 
                             sdLast.snarlDistance(graph, &ng,
                                  make_pair(lastEndID, !lastEndRev), 
                                   make_pair(lastEndID, lastEndRev)) });

                }
            }
            chainLoopRev.push_back(firstRevDist);
        }
        int64_t revLoopDist;

        if ( snarlRevInChain ) {

            revLoopDist = sd.snarlDistance(graph, &ng,
              make_pair (snarlStartID, snarlStartRev),
             make_pair (snarlStartID, !snarlStartRev));
        } else {
            revLoopDist = sd.snarlDistance( graph, &ng, 
              make_pair (snarlEndID, !snarlEndRev),
             make_pair (snarlEndID, snarlEndRev));
        }
 
 

        int64_t lastLoop = chainLoopRev.back();

        if (lastLoop == -1) {

            chainLoopRev.push_back(revLoopDist);

        } else {

            //Push the minimum of the loop distance of the current snarl and
            //the loop distance of the previous snarl + dist to and from loop 
            int64_t loopDistance = minPos({revLoopDist, lastLoop +  
                 sd.snarlDistance(graph, &ng,
                make_pair (snarlEndID, !snarlEndRev),
               make_pair (snarlStartID, !snarlStartRev))
            + 
                 sd.snarlDistance(graph, &ng,
                make_pair (snarlStartID, snarlStartRev),
                 make_pair (snarlEndID, snarlEndRev))});
            chainLoopRev.push_back(loopDistance);
        }
    }
    //Add forward loop distances 
   
    //Check if there is an edge traversing last node in chain fd -> rev 

    ChainIterator chainStartR = chain_rend(*chain);
    for (ChainIterator c = chain_rbegin(*chain); c != chainStartR; ++c) {
        //Loop through the chain in reverse
        const Snarl* snarl = c->first; 
        bool snarlRevInChain = c->second;
        id_t snarlStartID = snarl->start().node_id();
        bool snarlStartRev = snarl->start().backward();
        id_t snarlEndID = snarl->end().node_id();
        bool snarlEndRev = snarl->end().backward();
        auto& sd =  snarlDistances.at(make_pair(snarlStartID, snarlStartRev)); 
        NetGraph ng (snarl->start(), snarl->end(),sm->chains_of(snarl), graph);

                                      

        if (c == chain_rbegin(*chain)) {
            //If this is the last snarl in the chain, push loop for last node

            int64_t loopDistLast; 
            if (snarlRevInChain) {
       
                loopDistLast = sd.snarlDistance( graph, &ng,
                         make_pair(snarlStartID, !snarlStartRev), 
                         make_pair(snarlStartID, snarlStartRev));
            } else {

                loopDistLast = sd.snarlDistance(graph, &ng,
                         make_pair(snarlEndID, snarlEndRev), 
                         make_pair(snarlEndID, !snarlEndRev));
            }

            if (get_start_of(*chain).node_id() == get_end_of(*chain).node_id()) {
                //If the chain loops, might need distance from first snarl
                ChainIterator chainStart = chain_begin(*chain);
                const Snarl* firstSnarl = chainStart->first;
                bool firstSnarlRev = chainStart->second;
                    
                id_t firstStartID = firstSnarl->start().node_id();
                bool firstStartRev = firstSnarl->start().backward();
                id_t firstEndID = firstSnarl->end().node_id();
                bool firstEndRev = firstSnarl->end().backward();
                auto& sdFirst =  snarlDistances.at(make_pair(
                               firstStartID, firstStartRev)); 
                if (firstSnarlRev) {
                    loopDistLast = minPos({loopDistLast, 
                             sdFirst.snarlDistance(graph, &ng,
                                  make_pair(firstEndID, !firstEndRev), 
                                   make_pair(firstEndID, firstEndRev)) });
                } else {
                    loopDistLast = minPos({loopDistLast, 
                             sdFirst.snarlDistance(graph, &ng,
                                   make_pair(firstStartID, firstStartRev), 
                                   make_pair(firstStartID, !firstStartRev)) });

                }
              
            }
            chainLoopFd.push_back(loopDistLast);
        }

        int64_t fdLoopDist;


        if (snarlRevInChain) {
            //If the snarl is reversed in the chain
            fdLoopDist = sd.snarlDistance(graph, &ng,
                  make_pair (snarlEndID, !snarlEndRev),
                  make_pair (snarlEndID, snarlEndRev));
        } else {
            fdLoopDist = sd.snarlDistance(graph, &ng,
                  make_pair (snarlStartID, snarlStartRev),
                  make_pair (snarlStartID, !snarlStartRev));
        }

        int64_t lastLoop = chainLoopFd.back();

        if (lastLoop == -1) {

            chainLoopFd.push_back(fdLoopDist);

        } else {
        //push dist to end of snarl + loop dist + dist to start of snarl 

            int64_t loopDistance = minPos({fdLoopDist, lastLoop + 
                    sd.snarlDistance(graph, &ng, 
                                     make_pair(snarlEndID, !snarlEndRev),
                                     make_pair(snarlStartID, !snarlStartRev)) + 
                    sd.snarlDistance(graph, &ng,
                                     make_pair(snarlStartID, snarlStartRev), 
                                     make_pair(snarlEndID, snarlEndRev))});
            chainLoopFd.push_back(loopDistance); 
        }           
      
    }
    reverse(chainLoopFd.begin(), chainLoopFd.end()); 
 
    id_t s = get_start_of(*chain).node_id();
    id_t e = get_end_of(*chain).node_id();
    if (chainPrefixSum.size() > 2) { //If chain and not just one snarl
        chainDistances.insert(make_pair(s, 
                 ChainIndex(this, s, e, snarlToIndex, chainPrefixSum, 
                                                   chainLoopFd, chainLoopRev)));
    }
    //return length of entire chain
    return chainPrefixSum.back() + graph->get_length(graph->get_handle(e,
                                                                false));};



//////////////////    Distance Calculations

int64_t DistanceIndex::maxDistance(pos_t pos1, pos_t pos2) {
    //Get the upper bound of the distance between two positions

/* TODO: Make test that uses handle graph
    if (!(graph->has_node(get_id(pos1)) && graph->has_node(get_id(pos2)))) {
        throw runtime_error("Node not in graph");       
    }

    int64_t minDist = minDistance(pos1, pos2);

    if (minDist == -1) { 
        return -1;
    } else if (minDist >= maxIndex.cap) {
        return minDist;
    }
*/
     
    return maxIndex.maxDistance(pos1, pos2);

}

int64_t DistanceIndex::minDistance(pos_t pos1, pos_t pos2) {
    /*Minimum distance between positions not including the position itself*/
    const Snarl* snarl1 = snarlOf(get_id(pos1));
    const Snarl* snarl2 = snarlOf(get_id(pos2)); 
    return minDistance(snarl1, snarl2, pos1,pos2);
}
int64_t DistanceIndex::minDistance(const Snarl* snarl1, const Snarl* snarl2, 
                                   pos_t pos1, pos_t pos2) {
    /*Find the shortest distance between two positions
      pos1 and pos2 must be on nodes contained in snarl1/snarl2 */
    
    int64_t shortestDistance = -1; 

    id_t nodeID1 = get_id(pos1);
    bool nodeRev1 = is_rev(pos1);
    id_t nodeID2 = get_id(pos2); 
    bool nodeRev2 = is_rev(pos2);

    if (nodeID1 == nodeID2 && nodeRev1 == nodeRev2 ) {
        //if positions are on the same node and strand
        int64_t offset1 = get_offset(pos1);
        int64_t offset2 = get_offset(pos2);

        if (offset1 <= offset2) {
            shortestDistance = offset2-offset1+1; //+1 to be consistent
        }

    }



    const Snarl* commonAncestor = NULL; 


#ifdef printDistances
    cerr << endl << "Start distance calculation from " << nodeID1 << "->" <<
         nodeID2 << endl;

    cerr << "Shortes distance within same node: " << shortestDistance<<  endl;

    cerr << "Find common ancestor" << endl;
#endif


    //// Find common ancestor of the two snarls
    unordered_set<pair<id_t, bool>> ancestors1;//set of all ancestor snarls of node1
    const Snarl* ancestor1 = snarl1;

#ifdef printDistances
    cerr << "Ancestors of 1: ";
#endif


    while (ancestor1 != NULL) {
#ifdef printDistances
        cerr << ancestor1->start().node_id() << " ";
#endif
        ancestors1.emplace(make_pair(ancestor1->start().node_id(),
                                     ancestor1->start().backward()));
        ancestor1 = sm->parent_of(ancestor1);
    }


#ifdef printDistances
      cerr << endl << "ancestors of 2: ";
#endif


    const Snarl* ancestor2 = snarl2;
    while (ancestor2 != NULL) {


#ifdef printDistances
         cerr << ancestor2->start().node_id() << " ";
#endif


        if (ancestors1.count(make_pair(ancestor2->start().node_id(),
                                       ancestor2->start().backward())) > 0) { 
            commonAncestor = ancestor2;
            break;
        }
        ancestor2 = sm->parent_of(ancestor2); 
    }

#ifdef printDistances 
    cerr << endl;
    if (commonAncestor == NULL) {
        cerr << "common ancestor found: NULL" << endl;
    } else { 
        cerr << "common ancestor found: " << 
                           commonAncestor->start().node_id()<< endl;
    }

    cerr << "  Snarl1: " << snarl1->start().node_id() << " Snarl2: "
                                           << snarl2->start().node_id() << endl;
#endif


    //Find distances from pos1 and pos2 to ends of child snarls of ancestor
    pair<pair<int64_t, int64_t>, const Snarl*> p1 = 
                             distToCommonAncestor(snarl1, commonAncestor, pos1, false);
    pair<int64_t, int64_t> temp1 = p1.first; 
    snarl1 = p1.second;

    nodeRev1 = false;
    if (snarl1 != commonAncestor) {
        nodeID1 = snarl1->start().node_id();
        nodeRev1 = snarl1->start().backward();
    }
    int64_t distL1 = temp1.first; int64_t distR1 = temp1.second;
    
    pair<pair<int64_t, int64_t>, const Snarl*> p2 = 
                             distToCommonAncestor(snarl2, commonAncestor, pos2, true);
    nodeRev2 = false;
    pair<int64_t, int64_t> temp3 = p2.first; 
    snarl2 = p2.second;
    if (snarl2 != commonAncestor) {
        nodeID2 = snarl2->start().node_id();
        nodeRev2 = snarl2->start().backward();
    }
    int64_t distL2 = temp3.first; int64_t distR2 = temp3.second;
    

    id_t endID1 = snarl1->end().node_id();
    bool endRev1 = snarl1->end().backward();
    id_t endID2 = snarl2->end().node_id();
    bool endRev2 = snarl2->end().backward();

    //Snarl1 and snarl2 are children of common ancestor or common ancestor

#ifdef printDistances
    cerr << "Distances to snarl in common ancestor: " << distL1 << ", " <<
         distR1 << "   " << distL2 << ", " << distR2 << endl;
#endif
    int64_t chainDist = -1; 

    //Find shortest distance between boundary nodes of snarls containing pos
    // within the common ancestor snarl

    if (snarl1 != commonAncestor && snarl2 != commonAncestor && 
          sm->in_nontrivial_chain(snarl1) && sm->in_nontrivial_chain(snarl2)
           && sm->chain_of(snarl1) == sm->chain_of(snarl2)) {

        //If positions are in the same chain within common ancestor

        const Chain* chain = sm->chain_of(snarl1);
        id_t chainStartID = get_start_of(*chain).node_id();


        ChainIndex& chainDists = chainDistances.at( chainStartID); 

        //Distance from left of s1 (reverse), left of s2 (forward)
        int64_t d1 = chainDists.chainDistanceShort(graph,
               make_pair(nodeID1, !nodeRev1), 
               make_pair(nodeID2, nodeRev2), snarl1, snarl2);
        d1 = (distL1 == -1 || distL2 == -1 || d1 == -1) ? -1 : 
                                distL1 + distL2 + d1; 

        //Distance from left of s1 (reverse) to right of s2 (reverse)
        int64_t d2 = chainDists.chainDistanceShort(graph,
               make_pair(nodeID1, !nodeRev1), 
               make_pair(endID2, !endRev2), snarl1, snarl2);
        if (nodeID1 == endID2) {
            //If snarls share a node, chainDistanceShort returns length of 
            //shared node
            d2 = (distL1 == -1 || distR2 == -1 || d2 == -1) ? -1 : 
                                   distL1 + distR2 - d2; 
        } else {
            d2 = (distL1 == -1 || distR2 == -1 || d2 == -1) ? -1 : 
                                   distL1 + distR2 + d2; 
        }

        //Distance from right of s1 (fd) to left of s2 (fd)
        int64_t d3 = chainDists.chainDistanceShort(graph,
                  make_pair(endID1, endRev1), 
                  make_pair(nodeID2, nodeRev2), snarl1, snarl2);
        if (endID1 == nodeID2) {
            d3 = (distR1 == -1 || distL2 == -1 || d3 == -1) ? -1 : 
                                   distR1 + distL2 - d3; 
        } else {
            d3 = (distR1 == -1 || distL2 == -1 || d3 == -1) ? -1 : 
                                   distR1 + distL2 + d3; 
        }

        //Distance from right of s1 (fd) to right of s2 (rev)
        int64_t d4 =  chainDists.chainDistanceShort(graph,
              make_pair(endID1, endRev1),
              make_pair(endID2, !endRev2), snarl1, snarl2);
        d4 = (distR1 == -1 || distR2 == -1 || d4 == -1) ? -1 : 
                                   distR1 + distR2 + d4;
        
                   
        chainDist = minPos({d1, d2, d3, d4});

#ifdef printDistances
        cerr << "    Possible distances within chain: " << d1 << " " << d2
             << " " << d3 << " " << d4 << endl;
        cerr << "Chain distance in common ancestor: " << chainDist << endl;
        
#endif

    }
    if (commonAncestor == NULL) {
        chainDist = chainDist == -1 ? -1 : chainDist - 1;
        shortestDistance = shortestDistance == -1 ? -1 : shortestDistance - 1;
        return minPos({chainDist, shortestDistance});
    }  

    //Get dist from pos1 to ends of its chain 
    if (snarl1 != commonAncestor && sm->in_nontrivial_chain(snarl1)) {
        const Chain* chain = sm->chain_of(snarl1);



        Visit chainStart = get_start_of(*chain);
        Visit chainEnd = get_end_of(*chain);



        pair<id_t, bool> chainStartIn (chainStart.node_id(), 
                                       chainStart.backward());
        pair<id_t, bool> chainEndIn (chainEnd.node_id(), !chainEnd.backward());

        const Snarl* startSnarl = sm->into_which_snarl(chainStart);
        const Snarl* endSnarl = sm->into_which_snarl(chainEndIn.first, 
                                                     chainEndIn.second);

        ChainIndex& chainDists = chainDistances.at( chainStartIn.first);
        int64_t dsl = chainDists.chainDistance(chainStartIn, 
                              make_pair(nodeID1, nodeRev1), startSnarl, snarl1);
        int64_t dsr = chainDists.chainDistance(chainStartIn, 
                               make_pair(endID1, !endRev1), startSnarl, snarl1);
        int64_t der = chainDists.chainDistance(chainEndIn, 
                                 make_pair(endID1, !endRev1), endSnarl, snarl1);
        int64_t del = chainDists.chainDistance(chainEndIn,
                                make_pair(nodeID1, nodeRev1), endSnarl, snarl1);

        dsl = dsl == -1 || distL1 == -1? -1 : distL1 + dsl; 
        dsr =  dsr == -1 || distR1 == -1? -1 : distR1 + dsr; 
        der = der == -1 || distR1 == -1? -1 : distR1 + der; 
        del = del == -1 || distL1 == -1? -1 : distL1 + del; 
 
        distL1 = minPos({dsr, dsl});
        distR1 = minPos({der, del}); 



        nodeID1 =  chainStartIn.first;
        nodeRev1 = chainStartIn.second;
    }
    //Get dist from pos2 to ends of its chain 
    if (snarl2 != commonAncestor && sm->in_nontrivial_chain(snarl2)) {
        const Chain* chain = sm->chain_of(snarl2);
 
        Visit chainStart = get_start_of(*chain);
        Visit chainEnd = get_end_of(*chain);
        
     
        pair<id_t, bool> chainStartIn (chainStart.node_id(), chainStart.backward());
        pair<id_t, bool> chainEndIn (chainEnd.node_id(), !chainEnd.backward());

        const Snarl* startSnarl = sm->into_which_snarl(chainStart);
        const Snarl* endSnarl = sm->into_which_snarl(chainEndIn.first, 
                                                     chainEndIn.second);

        ChainIndex& chainDists = chainDistances.at( chainStartIn.first);


        int64_t dsl = chainDists.chainDistance(chainStartIn, 
                             make_pair(nodeID2, nodeRev2), startSnarl, snarl2);
        int64_t dsr = chainDists.chainDistance(chainStartIn,
                               make_pair(endID2, !endRev2), startSnarl, snarl2);
        int64_t der = chainDists.chainDistance(chainEndIn, 
                                 make_pair(endID2, !endRev2), endSnarl, snarl2);
        int64_t del = chainDists.chainDistance(chainEndIn, make_pair(nodeID2, 
                                             nodeRev2), endSnarl, snarl2);

        dsl = dsl == -1 || distL2 == -1? -1 : distL2 + dsl; 
        dsr = dsr == -1 || distR2 == -1? -1 : distR2 + dsr; 
        der = der == -1 || distR2 == -1? -1 : distR2 + der; 
        del = del == -1 || distL2 == -1? -1 : distL2 + del; 

        distL2 = minPos({dsr, dsl});
        distR2 = minPos({der, del});

        nodeID2 = chainStartIn.first; 
        nodeRev2 = get_start_of(*chain).backward();
    }
          
   
 
#ifdef printDistances
    cerr << "Distances to node in common ancestor: " << distL1 << ", " << distR1
              << "   " << distL2 << ", " << distR2 << endl;
#endif
    //Both nodes are nodes in common ancestor

    //Get distance between ends of nodes in common ancestor snarl
    NetGraph ng = NetGraph(commonAncestor->start(), 
                commonAncestor->end(),sm->chains_of(commonAncestor), graph);


    auto snarlDistsTmp = snarlDistances.find(make_pair(
                           commonAncestor->start().node_id(),
                                       commonAncestor->start().backward()));
    if (snarlDistsTmp == snarlDistances.end()) {
        snarlDistsTmp = snarlDistances.find(make_pair(
                           commonAncestor->end().node_id(),
                                       !commonAncestor->end().backward()));
    }
    SnarlIndex& snarlDists = snarlDistsTmp->second;


    int64_t d1 = snarlDists.snarlDistanceShort(
                make_pair(nodeID1, nodeRev1), make_pair(nodeID2, nodeRev2));
    d1 = (distR1 == -1 || distL2 == -1 || d1 == -1) ? -1 : 
                                                  distR1 + distL2 + d1; 

    int64_t d2 = snarlDists.snarlDistanceShort(
                   make_pair(nodeID1, nodeRev1), make_pair(nodeID2, !nodeRev2));

    d2 = (distR1 == -1 || distR2 == -1 || d2 == -1) ? -1 :
                                                   distR1 + distR2 + d2;
    int64_t d3 = snarlDists.snarlDistanceShort(
                make_pair(nodeID1, !nodeRev1), make_pair(nodeID2, nodeRev2));
    d3 = (distL1 == -1 || distL2 == -1 || d3 == -1) ? -1 : 
                                              distL1 + distL2 + d3; 
    int64_t d4 = snarlDists.snarlDistanceShort(
                make_pair(nodeID1, !nodeRev1), make_pair(nodeID2, !nodeRev2));
    d4 = (distL1 == -1 || distR2 == -1 || d4 == -1) ? -1 : 
                                                  distL1 + distR2 + d4; 

    shortestDistance =  minPos({d1, d2, d3, d4, chainDist, shortestDistance});

#ifdef printDistances
    cerr << "Distances within common ancestor: " << d1 << ", " << d2
                                        << ", " << d3 << ", " << d4 << endl;
    cerr << "Shortest dist only up to  common ancestor: " << shortestDistance
        << endl;
#endif
    
    
    //Find distances to the ends of the common ancestor snarl
    pair<int64_t, int64_t> endDists = snarlDists.distToEnds(graph, &ng, nodeID1, nodeRev1, distL1, distR1);
    distL1 = endDists.first;
    distR1 = endDists.second;

    endDists = snarlDists.distToEnds(graph, &ng, nodeID2, nodeRev2, distL2, distR2);
    distL2 = endDists.first;
    distR2 = endDists.second;
     
#ifdef printDistances
    cerr << "Distances to ends of common ancestor: " << distL1 << " " << distR1
         << " " << distL2 << " " << distR2 
        << endl;
#endif   

    const Snarl* currSnarl = commonAncestor;
    const Snarl* parentSnarl = sm->parent_of(currSnarl);
    id_t startID = currSnarl->start().node_id();
    id_t startRev = currSnarl->start().backward(); //pointing into snarl
    id_t endID = currSnarl->end().node_id();
    id_t endRev = currSnarl->end().backward();     //pointing out

    /*shortestDistance is now the shortest distance only traversing up to the 
      most recent common ancestor.   
      
      currSnarl is the common ancestor, start/end ID are a node in the 
      common ancestor, distances are up to a node in the common ancestor
      Traverse up to root and check for path at each level
    */
    
    while ( currSnarl != NULL) {

            
        if (sm->in_nontrivial_chain(currSnarl)) {
            //Find paths between ends of current chain

            const Chain* currChain= sm->chain_of(currSnarl);
            ChainIndex& chainDists = chainDistances.at(
                                            get_start_of(*currChain).node_id());

            //Distance from start (reverse) to start (forward)
            int64_t d1 = chainDists.chainDistanceShort(graph,
                  make_pair(startID, !startRev),
                   make_pair(startID, startRev), currSnarl, currSnarl);
            d1 = (distL1 == -1 || distL2 == -1 || d1 == -1) ? -1 : 
                                distL1 + distL2 + d1; 

            //Distance from start (reverse) to end (reverse)
            int64_t d = chainDists.chainDistanceShort(graph,
               make_pair(startID, !startRev), 
               make_pair(endID, !endRev), currSnarl, currSnarl);

            d2 = (distL1 == -1 || distR2 == -1 || d == -1) ? -1 : 
                                   distL1 + distR2 + d; 
            d3 = (distR1 == -1 || distL2 == -1 || d == -1) ? -1 : 
                                   distR1 + distL2 + d; 

           //Distance from end (fd) to end (rev)
            int64_t d4 =  chainDists.chainDistanceShort(graph,
                  make_pair(endID, endRev),
                  make_pair(endID, !endRev), currSnarl, currSnarl);
            d4 = (distR1 == -1 || distR2 == -1 || d4 == -1) ? -1 : 
                                   distR1 + distR2 + d4;
   
                  
            shortestDistance = minPos({shortestDistance, d1, d2, d3, d4});
           
            //Find distances to ends of the current chain

                    
            Visit chainStart = get_start_of(*currChain); //facing in
            Visit chainEnd = get_end_of(*currChain);//facing out

            pair<id_t, bool> chainStartIn (chainStart.node_id(), 
                                      chainStart.backward());
            pair<id_t, bool> chainEndIn (chainEnd.node_id(), 
                                        !chainEnd.backward());
            const Snarl* startSnarl = sm->into_which_snarl(chainStart);
            const Snarl* endSnarl = sm->into_which_snarl(chainEndIn.first, 
                                                         chainEndIn.second);


            int64_t dsl = chainDists.chainDistance(chainStartIn, 
                           make_pair(startID, startRev), startSnarl, currSnarl);
            int64_t dsr = chainDists.chainDistance(chainStartIn, 
                          make_pair(endID, !endRev), startSnarl, currSnarl);
            int64_t der = chainDists.chainDistance(chainEndIn, 
                          make_pair(endID, !endRev), endSnarl, currSnarl);
            int64_t del = chainDists.chainDistance(chainEndIn,
                            make_pair(startID, startRev), endSnarl, currSnarl);


            distL1 = minPos({dsr == -1 || distR1 == -1? -1 : distR1 + dsr, 
                             dsl == -1 || distL1 == -1? -1 : distL1 + dsl});
            distL2 = minPos({dsr == -1 || distR2 == -1? -1 : distR2 + dsr, 
                             dsl == -1 || distL2 == -1? -1 : distL2 + dsl});

            distR1 = minPos({der == -1 || distR1 == -1? -1 : distR1 + der, 
                             del == -1 || distL1 == -1? -1 : distL1 + del});
            distR2 = minPos({der == -1 || distR2 == -1? -1 : distR2 + der, 
                             del == -1 || distL2 == -1? -1 : distL2 + del});




            startID = chainStartIn.first;
            startRev = chainStartIn.second;
            endID = chainStartIn.first;
            endRev = !chainEndIn.second;

#ifdef printDistances
    cerr << "At chain " << startID << " dists to ends: " << distL1 << " " << 
         distR1 << " " << distL2 << " " << distR2 << endl;
    cerr << "distances: "  << d1 << " " << d2 << " " << d3 << " " << d4 << endl;
    cerr << " Shortest distance : "
          << shortestDistance << endl;
#endif  
        }
   
        if (parentSnarl == NULL) {break;}

        auto snarlDistsTmp = snarlDistances.find(
                                   make_pair(parentSnarl->start().node_id(),
                                             parentSnarl->start().backward()));
        if (snarlDistsTmp == snarlDistances.end()) {
            snarlDistsTmp = snarlDistances.find(
                                   make_pair(parentSnarl->end().node_id(),
                                           ! parentSnarl->end().backward()));
        }
        SnarlIndex& snarlDists = snarlDistsTmp->second;


        NetGraph ng = NetGraph(parentSnarl->start(), 
                parentSnarl->end(),sm->chains_of(parentSnarl), graph);

        //Find the shortest distance within the snarl

        //Dist from start to start
        d1 = snarlDists.snarlDistanceShort( 
                make_pair(startID, !startRev), make_pair(startID, startRev));
        d1 = (distL1 == -1 || distL2 == -1 || d1 == -1) ? -1 : 
                                                   distL1 + distL2 + d1; 
    
        //Dist from end to end
        d2 = snarlDists.snarlDistanceShort(
                   make_pair(startID, startRev), make_pair(startID, !startRev));

         d2 = (distR1 == -1 || distR2 == -1 || d2 == -1) ? -1 :
                                                   distR1 + distR2 + d2;
         //Dist from start to end
         int64_t dtemp = snarlDists.snarlDistanceShort(
                make_pair(startID, startRev), make_pair(startID, startRev));
         d3 = (distL1 == -1 || distR2 == -1 || dtemp == -1) ? -1 : 
                                              distL1 + distR2 + dtemp; 
         d4 = (distR1 == -1 || distL2 == -1 || dtemp == -1) ? -1 : 
                                                  distR1 + distL2 + dtemp; 


         shortestDistance =  minPos({d1, d2, d3, d4, shortestDistance});


        //Find the distances to ends of the snarl
        pair<int64_t, int64_t> endDists1 = snarlDists.distToEnds(graph, &ng, startID, 
                                                 startRev, distL1, distR1);
        distL1= endDists1.first; distR1= endDists1.second;

        pair<int64_t, int64_t> endDists2 = snarlDists.distToEnds(graph, &ng, startID,
                                                  startRev, distL2, distR2);
        distL2= endDists2.first; distR2= endDists2.second;

        startID = parentSnarl->start().node_id();
        startRev = parentSnarl->start().backward();
        endID = parentSnarl->end().node_id();
        endRev = parentSnarl->end().backward();


#ifdef printDistances
    cerr << "At snarl " << startID << " dists to ends: " << distL1 << " " << 
         distR1 << " " << distL2 << " " << distR2 << " Shortest distance : "
          << shortestDistance << endl;
#endif  
        currSnarl = parentSnarl;
        parentSnarl = sm->parent_of(currSnarl);
    }

    shortestDistance = shortestDistance == -1 ? -1 : shortestDistance - 1;
    return shortestDistance;

};


pair<pair<int64_t, int64_t>, const Snarl*> DistanceIndex::distToCommonAncestor(
          const Snarl* snarl, const Snarl* commonAncestor, pos_t& pos, bool rev){

    /* Find the distance from pos to either end of a snarl node in 
       commonAncestor. Doesn't find the distance to ends of a chain child of 
       common ancestor.
       Return the two distances and the Snarl whose parent is the is
       commonAncestor or commonAncestor if the position is on a node (not a
       snarl) in commonAncestor
    */

    int64_t distL; //Dist from pos1 to boundaries of curr snarl 
    int64_t distR; //To start and end of snarls, not necessarily left/right
    id_t nodeID = get_id(pos); 

    int64_t offset = get_offset(pos);
    #ifdef printDistances
    cerr << "Dist to common ancestor" << "node " << get_id(pos) << " offset " <<           offset <<" reversed " << is_rev(pos) 
              << " in snarl " << snarl->start().node_id();
    rev ? cerr << " end pos" << endl : cerr << " start pos" << endl;
    
    #endif
    if (is_rev(pos)) {
        distR = offset+1;
        distL = graph->get_length(graph->get_handle(get_id(pos), false))-offset;
    } else {
        distL = offset+1;
        distR = graph->get_length(graph->get_handle(get_id(pos), false))-offset;
    }
    if (rev == is_rev(pos)) {
        distL = -1;
    } else {
        distR = -1;
    }
    #ifdef printDistances
        cerr << "start pos: " << get_offset(pos) << "-> start: " << distL << 
               ", end: " << distR << endl;
    #endif
 
    if (commonAncestor != NULL &&
        snarl->start().node_id() == commonAncestor->start().node_id() &&
        snarl->start().backward() == commonAncestor->start().backward()) {
        /*If the node is a node in commonAncestor, return the distances to 
           the ends of the node
        */
        return make_pair(make_pair(distL, distR), snarl);
    }

    id_t startID = snarl->start().node_id(); 
    bool startRev = snarl->start().backward();

    auto snarlDistsTmp = snarlDistances.find(make_pair(startID, startRev));
    if (snarlDistsTmp == snarlDistances.end()) {
        snarlDistsTmp = snarlDistances.find(make_pair(startID, !startRev));
    }
    SnarlIndex& snarlDists = snarlDistsTmp->second;


    NetGraph ng (snarl->start(), snarl->end(), sm->chains_of(snarl), graph);
 
    pair<int64_t, int64_t> endDists = snarlDists.distToEnds(graph, &ng, 
                                             nodeID, false, distL, distR);
    distL = endDists.first;
    distR = endDists.second;

    #ifdef printDistances
    cerr << nodeID << "->" << startID << ": " << distL << ", " << distR << endl;
    #endif

    nodeID = startID;
    bool nodeRev = startRev;

    while ((sm->parent_of(snarl) != NULL && commonAncestor == NULL) || 
           (commonAncestor != NULL && 
           !(sm->parent_of(snarl)->start().node_id() == commonAncestor->start().node_id() &&
           sm->parent_of(snarl)->start().backward() == commonAncestor->start().backward()))) {
        //While snarl's parent doesn't equal common ancestor

        int64_t dsl; int64_t dsr; int64_t der; int64_t del;

        if (sm->in_nontrivial_chain(snarl)) {
            //Get distances to ends of chain
            id_t endID = snarl->end().node_id();
            bool endRev = snarl->end().backward();
            const Chain* chain = sm->chain_of(snarl);
 
            Visit chainStart = get_start_of(*chain);
            Visit chainEnd = get_end_of(*chain);
        
     
            pair<id_t, bool> chainStartIn (chainStart.node_id(), 
                                           chainStart.backward());
            pair<id_t, bool> chainEndIn (chainEnd.node_id(),
                                         !chainEnd.backward());

            const Snarl* startSnarl = sm->into_which_snarl(chainStart);
            const Snarl* endSnarl = sm->into_which_snarl(chainEndIn.first,
                                                         chainEndIn.second);

            ChainIndex& chainDists = chainDistances.at( chainStartIn.first);


            int64_t dsl = chainDists.chainDistance(chainStartIn, 
                              make_pair(nodeID, nodeRev), startSnarl, snarl);
            int64_t dsr = chainDists.chainDistance(chainStartIn, 
                            make_pair(endID, !endRev), startSnarl, snarl);
            int64_t der = chainDists.chainDistance(chainEndIn, 
                            make_pair(endID, !endRev), endSnarl, snarl);
            int64_t del = chainDists.chainDistance(chainEndIn,
                             make_pair(nodeID, nodeRev), endSnarl, snarl);

            dsl = dsl == -1 || distL == -1? -1 : distL + dsl; 
            dsr = dsr == -1 || distR == -1? -1 : distR + dsr; 
            der = der == -1 || distR == -1? -1 : distR + der; 
            del = del == -1 || distL == -1? -1 : distL + del; 
 
            distL = minPos({dsr, dsl});
            distR = minPos({der, del});


            nodeID = chainStartIn.first; 
            nodeRev = chainStartIn.second;
    #ifdef printDistances
        cerr << nodeID << "->" << chainStartIn.first << ": " << distL << ", " << distR
                                                                   << endl;
    #endif
        }
     
        //Get distances to ends of parent snarl
        snarl = sm->parent_of(snarl);
        id_t startNodeID = snarl->start().node_id();
        id_t startNodeRev = snarl->start().backward();
            
        auto snarlDistsTmp = snarlDistances.find(
                                          make_pair(startNodeID, startNodeRev));
        if (snarlDistsTmp == snarlDistances.end()) {
            snarlDistsTmp = snarlDistances.find(
                                         make_pair(startNodeID, !startNodeRev));
        }
        SnarlIndex& snarlDists = snarlDistsTmp->second;

        pair<int64_t, int64_t> endDists = snarlDists.distToEnds(
                                                graph, &ng, nodeID, nodeRev, distL, distR);
          
        distL = endDists.first;
        distR = endDists.second;
    #ifdef printDistances
        cerr << nodeID << "->" << startNodeID << ": " << distL << ", " << distR 
            << endl;
    #endif
        nodeID = startNodeID;
        nodeRev = startNodeRev;
    }
    return make_pair(make_pair(distL, distR), snarl);
};

int64_t DistanceIndex::minPos (vector<int64_t> vals) {
    /*return the minimum value in vals that is not -1, returns -1 if all
     values are -1 */
    return accumulate(vals.begin(), vals.end(), -1, 
          [](int x, int y) {if (x==-1) {return y;} 
                            else if (y == -1) {return x;}
                            else {return min(x, y);}} 
          ); 
   
};


const Snarl* DistanceIndex::snarlOf (id_t nodeID) {
    /*Given a node id, return the snarl that contains the node*/

    int64_t uintSID = nodeToSnarl[nodeID - minNodeID];
    const Snarl* s = sm->into_which_snarl(uintSID>>1, (uintSID % 2 == 1));
    return s;


}


void DistanceIndex::printSelf() {
    cerr << "Nodes : Snarls" << endl;
    for (size_t i = 0 ; i < nodeToSnarl.size() ; i++) {
        cerr << i << " " << nodeToSnarl[i] << endl;
    }
    cerr << "Snarls: " << endl;
    for (auto snarls : snarlDistances) {
        snarls.second.printSelf();
    }
    cerr << endl << "Chains:" << endl;
    for (auto chains : chainDistances) {
        chains.second.printSelf();
    }
    cerr << endl << "Maximum distances" << endl;
    maxIndex.printSelf();
}


DistanceIndex::SnarlIndex::SnarlIndex(DistanceIndex* di,
                      unordered_set<pair<id_t, bool>>& 
                      allNodes,  pair<id_t, bool> start, pair<id_t,bool> end) {
    /*Constructor for SnarlIndex object that stores distances between
        nodes in a snarl */
    distIndex = di;

    //Assign all nodes+direction in snarl to an index
    size_t snarlDistances = 0;
    for (pair<id_t, bool> node: allNodes) {
        visitToIndex[node] = snarlDistances++;
    }

    int size = visitToIndex.size();
    //Initialize all distances to 0 (representing -1)
    util::assign(distances, int_vector<>(((size+1)*size)/2, 0));
    snarlStart = start;
    snarlEnd = end;

}

DistanceIndex::SnarlIndex::SnarlIndex(DistanceIndex* di,
                      vector<int64_t> v) {
    /*Constructor for SnarlIndex object given vector from serialization */
    
    distIndex = di;
    int64_t numNodes = v[0];
    int64_t start = v[1];
    snarlStart = (start < 0) ? make_pair( (id_t) abs(start), true) : 
                               make_pair( (id_t) abs(start), false);

    int64_t end = v[2];
    snarlEnd = (end < 0) ? make_pair( (id_t) abs(end), true) : 
                               make_pair( (id_t) abs(end), false);

    //Get visitToIndex
    for (size_t i = 0; i < numNodes; i ++ ) {

        int64_t n = v[i + 3]; //Node
        pair<id_t, bool> node = (n < 0) ? make_pair( (id_t) abs(n), true) : 
                               make_pair( (id_t) abs(n), false);
       visitToIndex[node] = i; 
    }

    //Get distance vector
    distances.resize(((numNodes+1) *numNodes) / 2);
    size_t j = 0;
    for (size_t i = numNodes + 3; i < v.size(); i++) {

        distances[j++] = v[i];

    }
    util::bit_compress(distances);

}

 vector<int64_t>DistanceIndex::SnarlIndex::toVector() {
    /*Convert contents of object to vector for serialization
      Vector contains a header of four ints: #nodes, start node, end node
                  a vector representing visitToIndex [node1, node2, ...] where                          the nodes are ordered by the index they map to
                  a vector representing distances*/

    vector<int64_t> v;// v (1, 0, sizeof(int64_t));
    size_t numNodes = visitToIndex.size();//number of node+directions
    v.resize(numNodes + distances.size() + 3); //store map, distances, header

    v[0] = (int64_t) numNodes;
    v[1] = snarlStart.second ? -(int64_t) snarlStart.first :
                                                 (int64_t) snarlStart.first;
    v[2] =  snarlEnd.second ? -(int64_t) snarlEnd.first :
                                                 (int64_t) snarlEnd.first;

    for (pair<pair<id_t, bool>, size_t> p : visitToIndex) {
        pair<id_t, bool> node = p.first;
        int64_t index = (int64_t) p.second;
        v[3 + index] = node.second ? -(int64_t) node.first : 
                                      (int64_t) node.first;
    }
   
 
    size_t i = 3 + numNodes;   
    for (int64_t d : distances) {
        v[i++] = d;
    }
    return v;

}


size_t DistanceIndex::SnarlIndex::index(pair<id_t, bool> start, 
                                            pair<id_t, bool> end) {
    /*Get the index of dist from start to end in a snarl distance matrix
      given the node ids + direction */
    size_t length = visitToIndex.size();
    size_t i1 = visitToIndex.at(start);
    size_t i2 = visitToIndex.at(make_pair(end.first, !end.second));
    if (i1 > i2) {
        //Reverse order of nodes
        i1 = visitToIndex.at(make_pair(end.first, !end.second));
        i2 = visitToIndex.at(make_pair(start.first, start.second));
    }
    
    size_t k = length - i1;
    return ( ((length + 1) * length ) / 2 ) - ( ((k + 1) * k ) / 2 ) + i2-i1;
}

void DistanceIndex::SnarlIndex::insertDistance(pair<id_t, bool> start, 
                                           pair<id_t, bool> end, int64_t dist) {
    //Assign distance between start and end
    size_t i = index(start, end);

    distances[i] = dist + 1;
}
   
int64_t DistanceIndex::SnarlIndex::snarlDistance(HandleGraph* graph, 
           NetGraph* ng, pair<id_t, bool> start, pair<id_t, bool> end) {
    /*Distance between beginnings of two nodes n1 and n2 in snarl
    */
    size_t i = index(start, end);
    int64_t dist = int64_t(distances[i])-1;
    return dist == -1 ? -1 : dist + nodeLength(graph, ng,start.first); 
}

int64_t DistanceIndex::SnarlIndex::snarlDistanceShort(pair<id_t, bool> start, pair<id_t, bool> end) {
    /*Distance between end of node n1 and beginning of node n2 in snarl
    */
    size_t i = index(start, end);
    return int64_t(distances[i]) - 1; 
}
int64_t DistanceIndex::SnarlIndex::nodeLength(HandleGraph* graph, 
           NetGraph* ng, id_t node){

    //Get the length of the node. 
 
    handle_t handle = ng->get_handle(node, false);   
    SnarlManager* sm = distIndex->sm;
//TODO: Should be able to use is_child
    //Get the snarl that the node represents, if any
    const Snarl* tempSnarl = sm->into_which_snarl(
                                              node, false);
    const Snarl* currSnarl = tempSnarl == NULL ? 
                        sm->into_which_snarl(node, true) : tempSnarl; 

    if (node!= snarlStart.first && node!= snarlEnd.first && currSnarl != NULL) {
        //If node represents a chain or snarl
        auto chainDists = distIndex->chainDistances.find(node);

        if (chainDists != distIndex->chainDistances.end()) {
            //If chain
            return chainDists->second.chainLength();
        } else {
            //If snarl
            auto snarlDists = distIndex->snarlDistances.find(make_pair(node, false));
            auto snarlDists1 = distIndex->snarlDistances.find(make_pair(node,
                                                              true));
            if (snarlDists != distIndex->snarlDistances.end()) {
                return snarlDists->second.snarlLength(graph, ng);
            } else {
                return snarlDists1->second.snarlLength(graph, ng);
            }
        }
         
    } else {
        return graph->get_length(graph->get_handle(node, false));
    }

}

int64_t DistanceIndex::SnarlIndex::snarlLength(HandleGraph* graph, NetGraph* ng) {   
    //Return the length of the snarl- dist from beginning of start to end of end
    int64_t dist = snarlDistance(graph, ng, snarlStart, snarlEnd);
    
     //length of snarl
    if (dist == -1) {
        return -1;
    } else {
        int64_t nodeLen = graph->get_length(graph->get_handle(snarlEnd.first, 
				                          snarlEnd.second));
        return dist + nodeLen; 
    }
 
}

pair<int64_t, int64_t> DistanceIndex::SnarlIndex::distToEnds(HandleGraph* graph,
           NetGraph* ng, id_t node, bool rev, int64_t distL, int64_t distR) {
    /* Given the distances to either end of a node, find the distances to 
       either end of the snarl
       Rev is true if the node is reversed in the snarl
    */  
    if (rev) {
        int64_t temp = distL;
        distL = distR;
        distR = temp;
    }
    
    pair<id_t, bool> snarlEndRev = make_pair(snarlEnd.first, !snarlEnd.second);
    int64_t dsl = snarlDistance(graph, ng, snarlStart, make_pair(node, false)); 

    int64_t dsr = snarlDistance( graph, ng, snarlStart, make_pair(node, true));

    int64_t der = snarlDistance( graph, ng, snarlEndRev, make_pair(node, true));

    int64_t del = snarlDistance(graph, ng, snarlEndRev, make_pair(node, false));

    //If the current node is already the start or end position of the snarl
    //then there may be no path between them in the index but the distance is 0
    if (node == snarlStart.first) {
        if( rev == snarlStart.second) {
            dsl = 0;
        } else {
            dsr = 0;
        }
    }

    if (node == snarlEnd.first) {
        if (rev == !snarlEnd.second) {//node is snarl end pointing in
            del = 0;
        } else {
            der = 0;
        }
    }
 
    dsl = dsl == -1 || distL == -1? -1 : distL + dsl; 
    dsr =  dsr == -1 || distR == -1? -1 : distR + dsr; 
    der = der == -1 || distR == -1? -1 : distR + der; 
    del = del == -1 || distL == -1? -1 : distL + del; 

    int64_t distStart = minPos({dsr, dsl});

    int64_t distEnd = minPos({der, del});

    return make_pair(distStart, distEnd);
}

void DistanceIndex::SnarlIndex::printSelf() {
    //Print the nodes contained in SnarlDistance
    cerr << endl;
     
    cerr << "Snarl Distances for snarl starting at " << snarlStart.first;
    if (snarlStart.second) {cerr << " reverse and ending at ";} 
    else                   { cerr << " forward and ending at ";}
    cerr << snarlEnd.first;
    if (snarlEnd.second) {cerr << " reverse";} 
    else {cerr << " forward";}

    cerr << endl << "Indices:" << endl;
    
    for (auto n : visitToIndex) {
        cerr << n.first.first << ", " << n.first.second << ": " << n.second << endl;
    }
    cerr << "Distances:" << endl;
    cerr << "    ";
    for (auto n : visitToIndex) {
        cerr << n.first.first;
        if (n.first.second) {
           cerr << "r   "; 
         } else {
             cerr << "f    ";
         }
    }
    cerr << endl;
    for (auto n1 : visitToIndex) {
        if (n1.first.second) {
        cerr << n1.first.first << "r    ";
        } else {
        cerr << n1.first.first << "f    ";
        }
        for (auto n2 : visitToIndex) {
            size_t length = visitToIndex.size();
            size_t i1 = visitToIndex.at(n1.first);
            size_t i2 = visitToIndex.at(n2.first);
            size_t k = length - i1;
   
            size_t i =  ( ((length + 1) * length ) / 2 ) - ( ((k + 1) * k ) / 2 ) + i2-i1;
//            if (i1 <= i2) {
                cerr << snarlDistanceShort(n1.first, n2.first) << "   "; 
//            } else { 
//                cerr << "-   "; 
//            }
        }
        cerr << endl;
    }
    cerr << endl; 
}

//ChainDistance methods
DistanceIndex::ChainIndex::ChainIndex(DistanceIndex* di, id_t start, id_t end, 
                  hash_map<id_t, size_t> s, 
                 vector<int64_t> p, vector<int64_t> fd, vector<int64_t> rev) {
    
    distIndex = di;
    chainStartID = start;
    chainEndID = end;

    snarlToIndex = move(s); 
    util::assign(prefixSum, int_vector<>(p.size()));
    util::assign(loopFd, int_vector<>(fd.size()));
    util::assign(loopRev, int_vector<>(rev.size()));
      
    for (size_t i = 0; i < p.size(); i++) {
        prefixSum[i] = p[i] + 1;
    }

   
    for (size_t i = 0; i < fd.size(); i++) {
        loopFd[i] = fd[i] + 1;
    }

   
    for (size_t i = 0; i < rev.size(); i++) {
        loopRev[i] = rev[i] + 1;
    }
  
    util::bit_compress(prefixSum);
    util::bit_compress(loopFd);
    util::bit_compress(loopRev);

}

DistanceIndex::ChainIndex::ChainIndex(DistanceIndex* di, vector<int64_t> v) {
    //Constructor given vector of ints from serialization 
    
    distIndex = di;

    size_t numNodes = v.size()-2 / 4;
  
    chainStartID = v[0];
    chainEndID = v[1];

    prefixSum.resize(numNodes);
    loopFd.resize(numNodes);
    loopRev.resize(numNodes);

    for (size_t i = 0; i <  numNodes; i ++ ) {
        id_t node = (id_t) v[i*4+2];
        if (snarlToIndex.find(node) == snarlToIndex.end()) {
            snarlToIndex[node] = i;
        } 
        
        prefixSum[i] = v[i*4 + 3];
        loopFd[i] = v[i*4 + 4];
        loopRev[i] = v[i*4 + 5];
   
    }
}

vector<int64_t> DistanceIndex::ChainIndex::toVector() {
    /*Convert contents into vector of ints for serialization
     Stored as [node_id, prefix sum1, prefix sum2, loopfd,loop rev, node_id2...]
     */
   
    int64_t numNodes = snarlToIndex.size();
    bool loops = chainStartID == chainEndID;
    if (loops) { numNodes = numNodes + 1; } 

    vector<int64_t> v;
    v.resize(numNodes * 4 + 2);
    v[0] = chainStartID;
    v[1] = chainEndID;

    for (int i = 0 ; i < numNodes ; i++) {
        v[4*i + 3] = prefixSum[i];
        v[4*i + 4] = loopFd[i];
        v[4*i + 5] = loopRev[i];
    }
   
    for (pair<id_t, size_t> p : snarlToIndex) {
        v[p.second * 4 + 2] = (int64_t) p.first; 
    }
    if (loops) {v[(numNodes-1) * 4 + 2] = v[0];} //Last node id is first
   
    return v;

}
int64_t DistanceIndex::ChainIndex::chainDistance(pair<id_t, bool> start, 
         pair<id_t, bool> end, const Snarl* startSnarl, const Snarl* endSnarl,
         bool recurse) {

    /*
     * Return the distance between the given node sides, except node side is
     * specified relative to the reading orientation of the chain that the
     * nodes are in. 
     */
    size_t i1;
    size_t i2;
    if (!recurse) {
//TODO: This is a bad way of doing this, change it
        if (start.first == -1) {
            i1 =  snarlToIndex.size();
            start.first = chainEndID;
        } else {
            i1 = snarlToIndex.at(start.first);
        }
        if (end.first == -1) {
            i2 =  snarlToIndex.size();
            end.first = chainEndID;
        } else {
            i2 = snarlToIndex.at(end.first);
        }
    } else { 
        i1 = snarlToIndex.at(start.first);
        i2 = snarlToIndex.at(end.first);
    }
    
    //Orientation of snarl in chain
    bool chainRev1 = distIndex->sm->chain_orientation_of(startSnarl);
    bool chainRev2 = distIndex->sm->chain_orientation_of(endSnarl);

    //The orientation of the node in the snarl
    //TODO: Might not work for unary snarls??????
    bool snarlRev1 =
                  i1 == snarlToIndex[startSnarl->start().node_id()] ?
                                     startSnarl->start().backward() :
                                     startSnarl->end().backward();
    bool snarlRev2 = i2 == snarlToIndex[endSnarl->start().node_id()] ?
                                      endSnarl->start().backward() :
                                      endSnarl->end().backward();
    //If the snarl is reversed in the chain, the node is traversed reverse of snarl 
    //orientation of node in chain
    chainRev1 = chainRev1 ? !snarlRev1 : snarlRev1;
    chainRev2 = chainRev2 ? !snarlRev2 : snarlRev2;

    bool rev1 = chainRev1 ? !start.second : start.second;
    bool rev2 = chainRev2 ? !end.second : end.second;
    int64_t loopDist = -1;

    if (chainStartID == chainEndID && i1 != i2 && recurse) {
        //If the chain loops

         if (i1 == 0) {

             loopDist = chainDistance(make_pair(-1, start.second), 
                                              end, startSnarl, endSnarl, false); 

         } else if (i2 == 0) { 

             loopDist = chainDistance(start, 
                     make_pair(-1, end.second), startSnarl, endSnarl, false);

         } else if (i1 < i2  && start.second) {
             //If path could pass through first node in reverse

             loopDist = 
                  chainDistance(start, make_pair(chainStartID, start.second), 
                                                  startSnarl, endSnarl, false) 
               + chainDistance(make_pair(-1, start.second), end, startSnarl,
                                                         endSnarl, false);

         } else if (i1 > i2 && !rev1) {

             loopDist =
                chainDistance(start, make_pair(-1, start.second), startSnarl, 
                                                             endSnarl, false)
              + chainDistance(make_pair(chainStartID, start.second), end, startSnarl, 
                                                           endSnarl, false); 
         } 
          
    }
    HandleGraph* graph = distIndex->graph;

    if ((!rev1 && !rev2)) {
        //If start and end are facing forward relative to the start of the chain
        if (i1 <= i2) {
            int64_t dNoRev = prefixSum[i2] - prefixSum[i1] ; 
            return minPos({loopDist, dNoRev});
        } else {
            int64_t revID1 = loopFd[i1] - 1;
            int64_t revID2 = loopRev[i2] - 1;
            int64_t len1 = graph->get_length(graph->get_handle(start.first,
                                                                 start.second));
            int64_t len2 = graph->get_length(graph->get_handle(end.first,
                                                                 end.second));
            int64_t chainDist = (prefixSum[i1] + len1) - 
                                (prefixSum[i2] + len2); 
            return minPos({loopDist, (revID1 == -1 || revID2 == -1) ? -1 : 
                    chainDist + revID1 + revID2}); 
        }

    } else if (rev1 && rev2 ){
        //If start and end are both reversed relative to the start of the chain
        if (i1 >= i2) {
            int64_t len1 = graph->get_length(graph->get_handle(start.first,
                                                                 start.second));
            int64_t len2 = graph->get_length(graph->get_handle(end.first,
                                                                 end.second));
            int64_t dNoRev = (prefixSum[i1] + len1) - (prefixSum[i2]+len2);
            return minPos({loopDist, dNoRev});
        } else {
            int64_t revID1 = loopRev[i1] - 1;
            int64_t revID2 = loopFd[i2] - 1;
            int64_t chainDist = prefixSum[i2] - prefixSum[i1]; 
            return minPos({loopDist, ((revID1 == -1 || revID2 == -1) ? -1 : 
                     chainDist+ revID1 + revID2)}); 
        }
    } else if (!rev1 && rev2) {
        //Start is forward, end is reversed
        if (i1 <= i2) {
            int64_t rev = loopFd[i2] - 1;
            int64_t chainDist = prefixSum[i2]- prefixSum[i1];
            return minPos({loopDist, ((rev == -1) ? -1 : rev + chainDist )});
        } else {
            int64_t rev = loopFd[i1] - 1;
            int64_t len1 = graph->get_length(graph->get_handle(start.first,
                                                                 start.second));
            int64_t len2 = graph->get_length(graph->get_handle(end.first,
                                                                 end.second));
            int64_t chainDist = (prefixSum[i1]+len1) - (prefixSum[i2]+len2);
            return minPos({loopDist, ((rev == -1) ? -1 : rev + chainDist )});
        }
        
    } else {
        //start is reverse, end is forward
        if (i1 <= i2) {
            int64_t rev = loopRev[i1] - 1;
            int64_t chainDist = prefixSum[i2] - prefixSum[i1];
            return minPos({loopDist, (rev == -1 ? -1 : rev + chainDist )});

            
        } else {
            int64_t rev = loopRev[i2] - 1; 
            int64_t len1 = graph->get_length(graph->get_handle(start.first,
                                                                 start.second));
            int64_t len2 = graph->get_length(graph->get_handle(end.first,
                                                                 end.second));
            int64_t chainDist = (prefixSum[i1]+len1) - (prefixSum[i2]+len2);
            return minPos({loopDist, ((rev == -1) ? -1 : rev + chainDist )});
        }
    }
}

int64_t DistanceIndex::ChainIndex::chainDistanceShort(HandleGraph* graph, 
          pair<id_t, bool> start, pair<id_t, bool> end, const Snarl* startSnarl,
          const Snarl* endSnarl) {
    /*Distance between end of start node to beginning of end node in chain
      or the distance from the end of the end node to the start of the start 
      node
      If start and end are the same node, then return the length of that node
       because the length is needed for the distance calculation and a negative
       distance would indicate no path.
    */ 
    int64_t d1 = chainDistance(start, end, startSnarl, endSnarl);
    int64_t d2 = chainDistance(make_pair(end.first, !end.second), 
                               make_pair(start.first, !start.second), 
                                  startSnarl, endSnarl);
    if (start == end) {
        //If two positions are on different snarls that share a node
        return graph->get_length(graph->get_handle(start.first, start.second));
         
    }
    if (d1 == -1 && d2 == -1) {
        return -1;
    } else if (d2 == -1) {
        return d1 - graph->get_length(graph->get_handle(start.first, 
                                                        start.second));
    } else if (d1 == -1) {
        return d2 - graph->get_length(graph->get_handle(end.first, end.second));
    } else {
        return min(d1 - graph->get_length(graph->get_handle(start.first, 
                                                        start.second)), 
             d2 - graph->get_length(graph->get_handle(end.first, end.second)));
    }
}

int64_t DistanceIndex::ChainIndex::chainLength() {

    //Get the length of a chain including length of last node
    //TODO: if there is a unary snarl then this should be -1
    HandleGraph* graph = distIndex->graph;
    return prefixSum[prefixSum.size()-1] - 1 + 
                        graph->get_length(graph->get_handle(chainEndID, false));
}

void DistanceIndex::ChainIndex::printSelf() {
    //Print the contenst of ChainDistance
   
    cerr << "ChainDistance Indices:" << endl;
    
    for (auto n : snarlToIndex) {
        cerr << n.first  << ": " << n.second << endl;
    }
    cerr << "Distances:" << endl;
    cerr << endl;
    for (auto n : prefixSum) {
        cerr << n << " ";
    }
    cerr << endl; 
    cerr << "Loop Forward:" << endl;
    cerr << endl;
    for (auto n : loopFd) {
        cerr << n << " ";
    }
    cerr << endl; 
    cerr << "Loop Reverse:" << endl;
    cerr << endl;
    for (auto n : loopRev) {
        cerr << n << " ";
    }
    cerr << endl;
}




///////////////////////   MAXIMUM DISTANCE   ///////////////////////////////////


DistanceIndex::MaxDistanceIndex::MaxDistanceIndex() {
//TODO: Need this to compile for some reason?
}
DistanceIndex::MaxDistanceIndex::MaxDistanceIndex(DistanceIndex* di, const vector<const Snarl*> chain, uint64_t c) {

    //Calculate maximum distance index


    distIndex = di;
    cap = c;
    int64_t maxNodeID = distIndex->maxNodeID;
    int64_t minNodeID = distIndex->minNodeID;
    int_vector<> n(maxNodeID- minNodeID + 1, 0);
    nodeToComponent = n;

    int_vector<> max(maxNodeID - minNodeID + 1, 0);
    int_vector<> minFd(maxNodeID - minNodeID + 1, 0);
    int_vector<> minRev(maxNodeID - minNodeID + 1, 0);

    /////// DFS to get connected componpents that are in cycles
    numCycles= findComponents(nodeToComponent, max, minFd, minRev, 0, true);

    //Find connected components of nodes not in cycles

    findComponents(nodeToComponent, max, minFd, minRev, 
                                            numCycles, false);

    maxDistances = max;
    
    minDistances.resize(maxNodeID - minNodeID + 1);
    for (size_t i = 0; i < minDistances.size(); i++) {
        uint64_t d1 = minFd[i];
        uint64_t d2 = minRev[i];
        uint64_t d;
        if (d1 == 0) { 
            d = d2; 
        } else if (d2 == 0) {
            d = d1;
        } else {
            d = min(d1, d2);
        }
        minDistances[i] = d; 
    }
    
    util::bit_compress(nodeToComponent);
    util::bit_compress(minDistances);
    util::bit_compress(maxDistances);

}

int64_t DistanceIndex::MaxDistanceIndex::maxDistance(pos_t pos1, pos_t pos2) {
    //Upper bound of distance between two positions
    
    id_t node1 = get_id(pos1);
    HandleGraph* graph = distIndex->graph;
    int64_t len1 = max(get_offset(pos1),  
          graph->get_length(graph->get_handle(node1, false)) - 
          get_offset(pos1)) + 1;

    id_t node2 = get_id(pos2);
    int64_t len2 = max(get_offset(pos2),  
           graph->get_length(graph->get_handle(node2, false)) - 
           get_offset(pos2)) + 1;

    int64_t minNodeID = distIndex->minNodeID; 

    //Return the max distance between nodes plus maximum length of nodes
    uint64_t comp1 = nodeToComponent[node1-minNodeID];
    uint64_t comp2 = nodeToComponent[node2-minNodeID];
    if (comp1 != comp2 || comp1 <= numCycles) {
        //If they are in separate components or both in a cyclic component
        return cap;

    }

    uint64_t max1 = maxDistances[node1-minNodeID];
    uint64_t max2 = maxDistances[node2-minNodeID];
    uint64_t min1 = minDistances[node1-minNodeID];
    uint64_t min2 = minDistances[node2-minNodeID];

    uint64_t d1 = max1 > min2 ? max1-min2 : 0;
    uint64_t d2 = max2 > min1 ? max2-min1 : 0;

    return len1 + len2 + max(d1, d2);
}

uint64_t DistanceIndex::MaxDistanceIndex::findComponents( 
        int_vector<>& nodeToComponent, int_vector<>& maxDists, 
        int_vector<>& minDistsFd, int_vector<>& minDistsRev, 
        uint64_t currComponent, bool onlyCycles                       ){

    /*Assign nodes to a component
     *If onlyCycles, assign all nodes to a component of connected cycles 
               if in a cycle, 0 otherwise
      If not onlyCycles, assign all unassigned nodes to a connected component

      Returns the maximum component number, the number of connected components
    */

    int64_t minNodeID = distIndex->minNodeID;
    HandleGraph* graph = distIndex->graph;
    int64_t maxNodeID = distIndex->maxNodeID;
    hash_set<pair<id_t, bool>> seen;
    auto findComp = [&](const handle_t& h)-> bool { 
        id_t i = graph->get_id(h);

        if (nodeToComponent[i - minNodeID] == 0) {


            bool loops = distIndex->loopDistance(make_pair(i, false), 
                                                 make_pair(i, false)) > -1;
            if (onlyCycles == loops)  {
            //If this node hasn't been seen before and if only counting cycles,

                currComponent++;
                //Next nodes to look at; going forward at the end, remove from end 
                list<pair<pair<id_t, bool>, bool>> nextNodes;
            
                //Arbitrarily assign direction for DAG
                nextNodes.push_back(make_pair(make_pair(i, true), true));
                nextNodes.push_back(make_pair(make_pair(i, false), false));
                unordered_set<pair<id_t, bool>> sinkNodes;//Sink nodes of DAG
                unordered_set<pair<id_t, bool>> sourceNodes;//Source
                pair<id_t, bool> currNode; 

    
                while (nextNodes.size() > 0) {
                    //For each reachable node
 
                    //Traverse going forward first          
                    pair<pair<id_t, bool>, bool> next =  nextNodes.back();
                    nextNodes.pop_back();

                    currNode = next.first;
                    bool forward = next.second;

                    if (seen.count(currNode) == 0) {
                        //That hasn't been seen before
                    
                        seen.insert(currNode);
                        bool added = false;

                        auto addNextNodes = [&](const handle_t& h)-> bool {
                            //Helper fn to get adjacent nodes

                            pair<id_t, bool> node = make_pair(
                            graph->get_id(h), graph->get_is_reverse(h));
                            int64_t edgeLoop = distIndex->loopDistance(
                                                           currNode, node) > -1;
                            int64_t nodeLoop = distIndex->loopDistance(
                                                               node, node) > -1;
 
                            if ( 
                                 ((onlyCycles && edgeLoop && nodeLoop) || 
                                  (!onlyCycles && !edgeLoop && !nodeLoop)) ){
                                //Add nodes whose edges are in loops

                                added = true;
                            if (seen.count(node) == 0 ) {
                                if (forward) {
                                    nextNodes.push_back(make_pair(node, 
                                                                  forward));
                                } else {
                                    nextNodes.push_front(make_pair(node, 
                                                                   forward));
                                }

                                if (seen.count(make_pair(node.first, 
                                                         !node.second))
                                    == 0) { 
                                    if (forward) {
                                        nextNodes.push_front(make_pair(
                                          make_pair( node.first, !node.second), 
                                                                    !forward));
                                    } else {
                                        nextNodes.push_back(make_pair(
                                           make_pair( node.first, !node.second),
                                                                     !forward));
                                    }
                                }
                                }
                            }
                            return true;
                        };


                        nodeToComponent[currNode.first-minNodeID] = currComponent;


                        handle_t handle =graph->get_handle(currNode.first, 
                                                       currNode.second);
 
                        //Add nodes that are connected by edges in loops
                        graph->follow_edges(handle, false, addNextNodes);

                        if (!added && forward) {
                            //If there were no outgoing edges and this was a sink
                            sinkNodes.insert(currNode);
                        } else if (!added && !forward) {
                            //If there were no outgoing edges and this was a sink
                            sourceNodes.insert(currNode);
                        }
                     
                    }
                }
                //Found all nodes in current component 
                if (!onlyCycles) {
                    if (sinkNodes.size() == 0) {
                        calculateMaxDistances(sourceNodes, nodeToComponent, 
                                             maxDists, minDistsFd, minDistsRev);
                    } else {
                        calculateMaxDistances(sinkNodes, nodeToComponent, 
                                             maxDists, minDistsFd, minDistsRev);
                    }
                }
            }
        }
        return true;
    };

    graph->for_each_handle(findComp);
    return currComponent;
}


void DistanceIndex::MaxDistanceIndex::calculateMaxDistances(
                           unordered_set<pair<id_t, bool>>& sinkNodes,  
                           int_vector<>& nodeToComponent,
                           int_vector<>& maxDists, int_vector<>& minDistsFd, 
                           int_vector<>& minDistsRev){

    /*Given all nodes in a connected component and a set of source/sink nodes 
      (pointing out),get the max and min distances from each node to a sink node
    */

    HandleGraph* graph = distIndex->graph;
    int64_t minNodeID = distIndex->minNodeID; 
    
    list<pair<pair<id_t, bool>, pair<uint64_t, uint64_t>>> nextNodes;
    //Nodes that return to the component after leaving it
    list<pair<id_t, bool>> returnNodes;
    hash_map<pair<id_t, bool>, pair<uint64_t, uint64_t>> returnNodeVals;
    bool returned = false;
    uint64_t currComp; 
    
    for (pair<id_t, bool> sink : sinkNodes) {
        //Sink nodes are pointing out of the DAG

        pair<id_t, bool> currNode = make_pair(sink.first, !sink.second);
        currComp = nodeToComponent[currNode.first-minNodeID];
        nextNodes.push_back(make_pair(currNode, make_pair(1, 1))); 
        uint64_t len = graph->get_length(graph->get_handle(sink.first,
                                                         sink.second));
        returnNodes.push_back(sink);
        returnNodeVals[sink] = make_pair(0, len+1);
        //If a path leaves curr component, new min dist will never be a minimum
     }


    uint64_t maxMin = 0; // Largest min distance

    hash_set<pair<id_t, bool>> seenNodes;//Nodes that have been seen
    hash_set<pair<id_t, bool>> seenLoops;//Nodes in loops that have been seen- only traverse each loop at most once

    while (returnNodes.size() != 0) {
            //Traverse graph from one sink node
        
 
        if (nextNodes.size() == 0) {returned = true;}
        pair<id_t, bool> currNode;
        uint64_t minDist;
        uint64_t maxDist;
        if (returned ) {
            //If finished everything reachable without leaving component
            currNode = returnNodes.front();
            returnNodes.pop_front(); 
            pair<uint64_t, uint64_t> vals = returnNodeVals[currNode]; 
            returnNodeVals.erase(currNode);
            minDist = vals.first;
            maxDist = vals.second;
            
        } else {
            //Haven't left component
            pair<pair<id_t, bool>, pair<uint64_t, uint64_t>> next = 
                                                            nextNodes.front();
            nextNodes.pop_front(); 
            currNode = next.first;
            minDist = next.second.first;
            maxDist = next.second.second;
            seenNodes.insert(currNode);
        }
        
            
        uint64_t oldMin;
        uint64_t oldMax;
        if (nodeToComponent[currNode.first-minNodeID] == currComp) {
            //If in the same component - update distances

            ////// Update minimum distances depending on orientation of node

            if (currNode.second) {
                oldMin = minDistsFd[currNode.first-minNodeID]; 
                minDist = oldMin == 0 ? minDist :  min(oldMin, minDist);
                if (minDist != 0) {
                    minDistsFd[currNode.first-minNodeID] = minDist;
                }
            } else {
                oldMin = minDistsRev[currNode.first-minNodeID];
                minDist = oldMin == 0 ? minDist :  min(oldMin, minDist); 
                if (minDist != 0) {
                    minDistsRev[currNode.first-minNodeID] = minDist;
                }
            }

            //Update maximum distance
            oldMax = maxDists[currNode.first-minNodeID];
            maxDist = oldMax == 0 ? maxDist : max(oldMax, maxDist);
            maxDists[currNode.first-minNodeID] = maxDist;

            maxMin = max(maxMin, minDist);

        } else {
            seenLoops.insert(currNode);
        }

        int64_t nodeLen = graph->get_length(graph->get_handle(
                                          currNode.first, currNode.second));

        auto addNextNodes = [&](const handle_t& h)-> bool {
            //Helper fn to get adjacent nodes

            pair<id_t, bool> node = make_pair(
                   graph->get_id(h), graph->get_is_reverse(h));

            uint64_t nodeComp = nodeToComponent[node.first - minNodeID]; 

            if ( nodeComp == currComp) {
                //If next node is in the same component
                    
                if ( nodeToComponent[currNode.first - minNodeID] != currComp){
                    //If this is re-entering the current component

                    auto vals = returnNodeVals.find(node);
                    if (vals == returnNodeVals.end()) {
                        returnNodes.push_back(node);
                        returnNodeVals.insert(make_pair(node, 
                                            make_pair(0, maxDist+nodeLen+cap)));
                    } else {

                        returnNodeVals.erase(node);
                        returnNodeVals.insert(make_pair(node, make_pair(0, 
                               max(maxDist+nodeLen+cap, vals->second.second))));
                    }

                } else {
                    //If already in current component
                     
                       
                    if ((oldMin == 0 || oldMax == 0|| (maxMin+cap > oldMax))
                                && (maxDist > oldMax || minDist < oldMin)) {
                        //If this node hasn't been seen before
                        //If it has, then if the old distance was not already greater than the cap
                        //TODO: Arbitrarily breaking long loops
        
                        if (returned ) {

                           //If left component
                            auto vals = returnNodeVals.find(node);
                            if (vals == returnNodeVals.end()) {

                                returnNodes.push_back(node);
                                returnNodeVals.insert(make_pair( node,
                                                make_pair(0, maxDist+nodeLen)));
                            } else {

                                returnNodeVals.erase(node);
                                returnNodeVals.insert(make_pair(node,
                                          make_pair(0, max(maxDist+nodeLen, 
                                                        vals->second.second))));
                            }

                        } else { 
                            //If in current component and never left

                            bool add = true;
                            auto checkNodes = [&](const handle_t& p)-> bool {
                                //Check if incoming nodes have been seen
                                pair<id_t, bool> prev = make_pair(
                                    graph->get_id(p), graph->get_is_reverse(p));

                                uint64_t prevComp = nodeToComponent[
                                                        prev.first - minNodeID];
                                if (prevComp == currComp &&
                                    seenNodes.count(prev) == 0) {
                                    //If prev node in currComp and hasn't been 
                                    //seen yet
                                    add = false;

                                } 
                                return true;
                            };
                             
                            graph->follow_edges(h, true, checkNodes);

                            if (add) {

                                nextNodes.push_back(make_pair(node, 
                                                  make_pair(minDist + nodeLen, 
                                                             maxDist+nodeLen)));
                            }
                        }
                    }
               }

            } else if ( maxDist < maxMin + cap && seenLoops.count(node) == 0) {
                //The next node is in a different component
                //If the max distance that could be found is less than cap
                auto vals = returnNodeVals.find(node);

                if (vals == returnNodeVals.end()) {
                    
                    returnNodes.push_back(node);
                    returnNodeVals.insert(make_pair(node, 
                                          make_pair(0, maxDist+nodeLen)));
                } else { 

                    returnNodeVals.erase(node);
                    returnNodeVals.insert(make_pair(node, make_pair(0, 
                                   max(maxDist+nodeLen, vals->second.second))));
                }

            }
            return true;
        };
           

        handle_t handle = graph->get_handle(currNode.first, 
                                                currNode.second);

        graph->follow_edges(handle, false, addNextNodes);

            
    } 
    
    
   
}

void DistanceIndex::MaxDistanceIndex::printSelf() {

    cerr << "Number of cyclic components: " << numCycles << endl 
    << "Components: " << endl;
    for (auto x : nodeToComponent) {cerr << x << " ";}
    cerr << endl
    << "Min distances: " << endl;
    for (auto x : minDistances) {cerr << x << " " ;}
    cerr << endl
    << "Max distances: " << endl;
    for (auto x : maxDistances) {cerr << x << " " ;}
    cerr << endl << endl;
}


int64_t DistanceIndex::loopDistance(
                 pair<id_t, bool> node1, pair<id_t, bool> node2) {
    const Snarl* snarl1 = snarlOf(node1.first);
    const Snarl* snarl2 = snarlOf(node2.first); 
    return loopDistance(snarl1, snarl2, node1, node2);
}

int64_t DistanceIndex::loopDistance(const Snarl* snarl1,const Snarl* snarl2,
                 pair<id_t, bool> node1, pair<id_t, bool> node2) {
    /*Find the minimum distance to loop through the given edge or, if node1 and
      node2 are the same, to loop through that node    */
 
/*TODO: make a test using handle graphs
    if (node1 != node2 && !graph->has_edge(
                                      NodeSide(node1.first, !node1.second),
                                     NodeSide(node2.first, node2.second))){
        //Edge must exist
        throw runtime_error("Edge does not exist");       

    }
*/

#ifdef indexTraverse 
cerr << endl << " NEW LOOP CALCULATION: " << node1.first <<  " TO " << node2.first << endl;
#endif
              
    int64_t minLoop = -1;

    int64_t distSRev = 0; //Dist to start of snarl traversing node backward 
    int64_t distSFd = -1; // not including the length of the node
    int64_t distERev = -1;
    int64_t distEFd = 0;
    int64_t distERev1 = -1;
    int64_t distSFd2 = -1;

 
    const Snarl* snarl;
    

   //Length of current node passing through original node
    int64_t nodeLen;
    if (node1 == node2) { //Same node - look for loop through the node

        nodeLen = graph->get_length(graph->get_handle(node1.first, false));

    } else { //Look for loop that uses given edge

        nodeLen = graph->get_length(graph->get_handle(node1.first, false)) + 
                  graph->get_length(graph->get_handle(node2.first, false));

    }

    const Snarl* snarl1Rev = node1.first == snarl1->start().node_id() ?
               sm->into_which_snarl(node1.first, !snarl1->start().backward()) :
               sm->into_which_snarl(node1.first, snarl1->end().backward());

    const Snarl* snarl2Rev = node2.first == snarl2->start().node_id() ? 
                sm->into_which_snarl(node2.first, !snarl2->start().backward()) :
                sm->into_which_snarl(node2.first, snarl2->end().backward());
 
    if (snarl1 == snarl2) {

        snarl = snarl1;

    }  else if (sm->chain_of(snarl1) == sm->chain_of(snarl2))  {
        //If the two snarls are on the same chain
       
        const Chain* chain = sm->chain_of(snarl1);
        if ((node1.first == get_start_of(*chain).node_id() && 
                                    node2.first == get_end_of(*chain).node_id())
               ||
             (node2.first == get_start_of(*chain).node_id() && 
                             node1.first == get_end_of(*chain).node_id())) {
            /*If the nodes are on opposite sides of the chain, then the edge is
              part of a loop through the whole chain */
            auto chainDists = chainDistances.at(get_start_of(*chain).node_id());

            return chainDists.chainLength();
        }

        //At least one node must be the boundary node of a snarl
        if (node1.first == snarl1->start().node_id() || 
            node1.first == snarl1->end().node_id()) {

            snarl = sm->into_which_snarl(node1.first, node1.second);

        } else if (node2.first == snarl2->start().node_id() || 
                   node2.first == snarl2->end().node_id()){

            snarl = sm->into_which_snarl(node2.first, !node2.second);

        }


    } else if (sm->parent_of(snarl1) == sm->parent_of(snarl2)) { 
        //Snarls share a common parent snarl but aren't on the same chain

        int64_t length1 = 0; //Size of the snarl or chain of node1
        if (sm->in_nontrivial_chain(snarl1)) {
                //If chain, node is already a boundary node of snarl in chain 

            const Chain* chain = sm->chain_of(snarl1);
            Visit startVisit = get_start_of(*chain);
            id_t chainStartID = startVisit.node_id();

            ChainIndex& chainDists = chainDistances.at(chainStartID);


            pair<id_t, bool> bound;
            if (node1.first == chainStartID) {
                //if node is first in chain, bound is end
                Visit end = get_end_of(*chain);
                bound = make_pair(end.node_id(), !end.backward());
            } else {
                //if node is end of chain, bound is start
                 
                bound = make_pair(chainStartID, startVisit.backward());
            }
            const Snarl* 
                boundSnarl = sm->into_which_snarl(bound.first, bound.second);
          
            distSRev = chainDists.chainDistance(bound, node1, boundSnarl, snarl1);
            length1 = chainDists.chainLength();

            distERev = chainDists.chainDistance(
                              make_pair(node1.first, !node1.second),
                              node1, snarl1, snarl1);
            distERev1 = distERev;
            node1 = make_pair(chainStartID, node1.second);

#ifdef indexTraverse 
cerr << "DISTANCES TO ENDS OF CHAIN OF NODE 1: " << distSRev << " " << distSFd
     << " " << distERev << " " << distEFd << endl;
#endif  

        } else {
            //Node 1 is in a snarl
            SnarlIndex& snarlDists = snarlDistances.at(make_pair(
                         snarl1->start().node_id(),snarl1->start().backward()));

            NetGraph ng (snarl1->start(), snarl1->end(),
                                                  sm->chains_of(snarl1), graph);

            pair<id_t, bool> bound;
            if (node1.first == snarl1->start().node_id()) {
                bound = make_pair(snarl1->end().node_id(), 
                               !snarl1->end().backward());
            } else {
                bound = make_pair(snarl1->start().node_id(), 
                               snarl1->start().backward());
            }
            distSRev = snarlDists.snarlDistance(graph, &ng, bound, node1);
            length1 = snarlDists.snarlLength(graph, &ng);

            distERev = snarlDists.snarlDistance(graph, &ng, 
                          make_pair(node1.first, !node1.second), node1);
            distERev1 = distERev;

            node1 = make_pair(snarl1->start().node_id(), node1.second);

#ifdef indexTraverse 
cerr << "DISTANCES TO ENDS OF SNARL OF NODE 1: " << distSRev << " " << distSFd 
     << " " << distERev << " " << distEFd << endl;
#endif  


        }


        int64_t length2 = 0; //Size of the snarl or chain of node1
        if (sm->in_nontrivial_chain(snarl2)) {
                //If chain, node is already a boundary node of snarl in chain 

            const Chain* chain = sm->chain_of(snarl2);
            Visit startVisit = get_start_of(*chain);
            id_t chainStartID = startVisit.node_id();

            ChainIndex& chainDists = chainDistances.at(chainStartID);

            pair<id_t, bool> bound;
            if (node2.first == chainStartID) {
                Visit endVisit = get_end_of(*chain);
                bound = make_pair(endVisit.node_id(), !endVisit.backward());

            } else {
                bound = make_pair(chainStartID, startVisit.backward());
            }
            const Snarl* boundSnarl =
                                sm->into_which_snarl(bound.first, bound.second);


            distEFd = chainDists.chainDistance(bound, 
                 make_pair(node2.first, !node2.second), boundSnarl, snarl2);
            length2 = chainDists.chainLength();


            distSFd = chainDists.chainDistance(node2,
                    make_pair(node2.first, !node2.second), boundSnarl, snarl2);
            distSFd2 = distSFd;

            node2 = make_pair(chainStartID, node2.second);

#ifdef indexTraverse 
cerr << "DISTANCES TO ENDS OF CHAIN OF NODE 2: " << distSRev << " " << distSFd
     << " " << distERev << " " << distEFd << endl;
#endif  


        } else {
            //Node 2 is in a snarl
            SnarlIndex& snarlDists = snarlDistances.at(make_pair(
                         snarl2->start().node_id(),snarl2->start().backward()));

            NetGraph ng (snarl2->start(), snarl2->end(),
                                                  sm->chains_of(snarl2), graph);

            pair<id_t, bool> bound;
            if (node2.first == snarl2->start().node_id()) {
                bound = make_pair(snarl2->end().node_id(), 
                               !snarl2->end().backward());
            } else {
                bound = make_pair(snarl2->start().node_id(),
                               snarl2->start().backward());
            }
            distEFd = snarlDists.snarlDistance(graph, &ng, bound, 
                                        make_pair(node2.first, !node2.second));
            length2 = snarlDists.snarlLength(graph, &ng);

            distSFd = snarlDists.snarlDistance(graph, &ng, node2, 
                                     make_pair(node2.first, !node2.second));
            distSFd2 = distSFd;
            node2 = make_pair(snarl2->start().node_id(), node2.second);

#ifdef indexTraverse 
cerr << "DISTANCES TO ENDS OF SNARL OF NODE 2: " << distSRev << " " << distSFd
      << " " << distERev << " " << distEFd << endl;
#endif  


        }

        distSFd = distSFd == -1 ? -1 : distSFd + length1;
        distERev = distERev == -1 ? -1 : distERev + length2; 


        snarl = sm->parent_of(snarl1);

#ifdef indexTraverse 
cerr << "DISTANCES: " << distSRev << " " << distSFd << " " << distERev << " " << distEFd << endl;       
#endif

    } else {
        //One snarl must be the parent of the other

        if (snarl1Rev != NULL && sm->parent_of(snarl2) == snarl1Rev) {
            //Snarl1 is in a chain, adjacent snarl contains snarl 2
            snarl1 = snarl1Rev;
        } else if (snarl2Rev != NULL && sm->parent_of(snarl1) == snarl2Rev) {
            snarl2 = snarl2Rev;
        }
        if (sm->parent_of(snarl1) == snarl2) {

            //Snarl1 is start or end of child snarl in snarl2
            //Switch the orientation of the edge and continue to next condition
    
            pair<id_t, bool> node1Rev = make_pair(node1.first, !node1.second);
            pair<id_t, bool> node2Rev = make_pair(node2.first, !node2.second);
            node1 = node2Rev;
            node2 = node1Rev;
            const Snarl* temp = snarl1;
            snarl1 = snarl2;
            snarl2 = temp;




        } 
        if (sm->parent_of(snarl2) == snarl1) {
            //Snarl2 is start or end of child snarl in snarl1
            if (sm->in_nontrivial_chain(snarl2)) {
                //If chain, node is already a boundary node of snarl in chain 


                const Chain* chain = sm->chain_of(snarl2);

                Visit startVisit = get_start_of(*chain);
                Visit endVisit = get_end_of(*chain);
                id_t chainStartID = startVisit.node_id();
                pair<id_t, bool> chainStart;
                pair<id_t, bool> chainEnd;

                ChainIndex& chainDists = chainDistances.at(chainStartID);
                if (chainStartID == node2.first) {

                    chainStart = make_pair( startVisit.node_id(), 
                                                startVisit.backward() );
                    chainEnd = make_pair( endVisit.node_id(), 
                                             ! endVisit.backward());

                } else {

                    //Assume start of chain is the side node was on
                    chainEnd = make_pair( startVisit.node_id(), 
                                                startVisit.backward() );
                    chainStart = make_pair( endVisit.node_id(), 
                                             ! endVisit.backward());
                }


                const Snarl* chainStartSnarl = sm->into_which_snarl(
                                       chainStart.first, chainStart.second);

                const Snarl* chainEndSnarl = sm->into_which_snarl(
                                       chainEnd.first, chainEnd.second);





                pair<id_t, bool> node2Rev = make_pair(node2.first, !node2.second);
                distSFd = chainDists.chainDistance(chainStart, node2Rev, 
                                                   chainStartSnarl, snarl2);
                distERev = chainDists.chainDistance(node2Rev, chainEnd,
                                                    snarl2, chainEndSnarl);
                distEFd = chainDists.chainDistance(chainEnd, node2Rev,
                                                   chainEndSnarl, snarl2);
        
#ifdef indexTraverse 
cerr << "DISTANCES IN CHILD CHAIN: " << distSRev << " " << distSFd << " " << distERev << " " << distEFd << endl;       
#endif  
                node2 = make_pair(chainStartID, node2.second); 

            } else {
                //If only snarl

                SnarlIndex& snarlDists = snarlDistances.at(make_pair(
                         snarl2->start().node_id(),snarl2->start().backward()));
             
                pair<id_t, bool> snarlStart = snarlDists.snarlStart;
                pair<id_t, bool> snarlEnd = snarlDists.snarlEnd;
                snarlEnd = make_pair(snarlEnd.first, !snarlEnd.second);

                NetGraph ng (snarl2->start(), snarl2->end(),
                                                  sm->chains_of(snarl2), graph);
  
                if (node2.first != snarlStart.first) {
                    auto temp = snarlStart;
                    snarlStart = snarlEnd;
                    snarlEnd = temp;
                }

                pair<id_t, bool> node2Rev = make_pair(node2.first, 
                                                      !node2.second);

                distSFd = snarlDists.snarlDistance(graph, &ng,
                                                         snarlStart, node2Rev);

                distEFd = snarlDists.snarlDistance(graph, &ng, snarlEnd, 
                                                                      node2Rev);

node2 = node2.first == snarl2->start().node_id() ? 
        make_pair(snarl2->start().node_id(), snarl2->start().backward()) 
     :  make_pair(snarl2->start().node_id(), !snarl2->start().backward());

#ifdef indexTraverse 
cerr << "DISTANCES IN CHILD SNARL " << snarl2->start().node_id() << " : " << distSRev << " " << distSFd << " " << distERev << " " << distEFd << endl;       
#endif  
            }

            snarl = snarl1;

            auto snarlDists = snarlDistances.at(make_pair(
                         snarl->start().node_id(),snarl->start().backward()));

            NetGraph ng = NetGraph(snarl->start(), 
                               snarl->end(),sm->chains_of(snarl), graph);

            pair<id_t, bool> node1Rev = make_pair(node1.first, !node1.second);
            pair<id_t, bool> node2Rev = make_pair(node2.first, !node2.second);
            //Update snarl, node, and node length

            int64_t distSL = snarlDists.snarlDistanceShort(node2Rev, node1Rev);
            int64_t distEL = snarlDists.snarlDistanceShort(node2Rev, node1);
            int64_t distSR = snarlDists.snarlDistanceShort(node2, node1Rev);
            int64_t distER = snarlDists.snarlDistanceShort(node2, node1);

            int64_t distSFdTemp = minPos({ 
                (distSFd == -1 || distSL == -1) ? -1 : distSFd + distSL});

            int64_t distERevTemp = minPos({ 
                (distERev == -1 || distEL == -1) ? -1 : distERev + distEL});

            distSFd2 = distSFdTemp;
            distSRev = 0;

            distSFd = distSFdTemp == -1 ? -1 : distSFdTemp + 
                         snarlDists.nodeLength(graph, &ng, node1.first);
            distERev = distERevTemp == -1 ? -1 : distERevTemp + 
                          snarlDists.nodeLength(graph, &ng, node2.first);

#ifdef indexTraverse 
cerr << "DISTANCES: " << distSRev << " " << distSFd << " " << distERev << " " << distEFd << endl;       
#endif
        }

     
    }

   
    while (snarl != NULL) {
        //Check each ancestor snarl for a loop

#ifdef indexTraverse 
cerr << "SNARL: " << snarl->start() << endl;
#endif
        NetGraph ng = NetGraph(snarl->start(), 
                               snarl->end(),sm->chains_of(snarl), graph);

        SnarlIndex& snarlDists = snarlDistances.at(make_pair(
                           snarl->start().node_id(),snarl->start().backward()));

        pair<id_t, bool> node1Rev = make_pair(node1.first, !node1.second);
        pair<id_t, bool> node2Rev = make_pair(node2.first, !node2.second);

        int64_t loop = minPos({
                snarlDists.snarlDistanceShort(node2, node1),
                snarlDists.snarlDistanceShort(node1Rev, node2Rev)});

        int64_t loopL = snarlDists.snarlDistanceShort(node1Rev, node2);
        int64_t loopR = snarlDists.snarlDistanceShort(node2, node1Rev); 
#ifdef indexTraverse 
cerr << "SNARL LOOPS: " << loop << " " << loopL << " " << loopR << endl;
 #endif
        int64_t loop1 = loop == -1 || distSRev == -1 || distEFd == -1 ? -1 :
                                           loop + distSRev + distEFd + nodeLen;
        int64_t loop2 = loop == -1 || distSFd == -1 || distERev == -1 ? -1 : 
                                           loop + distSFd + distERev + nodeLen;
        int64_t loop3 = -1;
        if (node1 == node2) {

            loopL = loopL == -1 || distSFd == -1 || distSRev == -1 ? -1 : 
                                           loopL + distSFd + distSRev + nodeLen;
            loopR = loopR == -1 || distEFd == -1 || distERev == -1 ? -1 : 
                                           loopR + distEFd + distERev + nodeLen;
        } else {

            loopL = loopL == -1 || distSFd2 == -1 || distSRev == -1 ? -1 : 
                                        loopL + distSFd2 + distSRev + nodeLen;
            loopR = loopR == -1 || distEFd == -1 || distERev1 == -1 ? -1 : 
                                          loopR + distEFd + distERev1 + nodeLen;
            loop3 = distSFd2 == -1 || distERev1 == -1 ? -1 : distSFd2 + distERev1 + nodeLen;
        }
        

#ifdef indexTraverse 
cerr << "    LOOP DISTANCES: " << loop3 << " " << loop1 << " " << loop2 << " " << loopL << " " << loopR << endl;
#endif
        minLoop = minPos({minLoop, loop1, loop2, loop3, loopL, loopR});
          

        //Update snarl, node, and node length
        int64_t distSL = (node1 == snarlDists.snarlStart) ? 0 :
             snarlDists.snarlDistance(graph, &ng, make_pair(
                   snarl->start().node_id(), snarl->start().backward()), node1);
        int64_t distSR = (node2.first == snarlDists.snarlStart.first &&
                          node2.second != snarlDists.snarlStart.second) ? 0 : 
              snarlDists.snarlDistance(graph, &ng, make_pair(
                snarl->start().node_id(), snarl->start().backward()), node2Rev);
        int64_t distEL = (node1.first == snarlDists.snarlEnd.first &&
                          node1.second != snarlDists.snarlEnd.second) ? 0 : 
                 snarlDists.snarlDistance(graph, &ng, make_pair(
                      snarl->end().node_id(), !snarl->end().backward()), node1);
        int64_t distER = (node2 == snarlDists.snarlEnd) ? 0 :
                  snarlDists.snarlDistance(graph, &ng, make_pair(
                   snarl->end().node_id(), !snarl->end().backward()), node2Rev);
  
#ifdef indexTraverse 
cerr << "DISTANCES IN SNARL " << snarl->start().node_id() << " : " << distSL << " " << distSR << " " << distEL << " " << distER << endl;       
#endif
        int64_t distSRevTemp = minPos({ 
               ((distSRev == -1 || distSL == -1) ? -1 : distSRev + distSL), 
               ((distERev == -1 || distSR == -1) ? -1 : distERev + distSR)});

        int64_t distSFdTemp = minPos({
                ((distSFd == -1 || distSL == -1) ? -1 : distSFd + distSL),
                ((distEFd == -1 || distSR == -1) ? -1 : distEFd + distSR) });

        int64_t distERevTemp = minPos({
                ((distSRev == -1 || distEL == -1) ? -1 : distSRev + distEL),
                ((distERev == -1 || distER == -1) ? -1 : distERev + distER) });

        int64_t distEFdTemp = minPos({ 
                ((distSFd == -1 || distEL == -1) ? -1 : distSFd + distEL), 
                ((distEFd == -1 || distER == -1) ? -1 : distEFd + distER) });

        if (node1 != node2) {
            int64_t distSL2 =  snarlDists.snarlDistance(graph, &ng, make_pair(
                   snarl->start().node_id(), snarl->start().backward()), node2);

            int64_t distSR1 = snarlDists.snarlDistance(graph, &ng, make_pair(
                snarl->start().node_id(), snarl->start().backward()), node1Rev);

            int64_t distEL2 = snarlDists.snarlDistance(graph, &ng, make_pair(
                      snarl->end().node_id(), !snarl->end().backward()), node2);
            int64_t distER1 = snarlDists.snarlDistance(graph, &ng, make_pair(
                   snarl->end().node_id(), !snarl->end().backward()), node1Rev);

            distSRevTemp = minPos({distSRevTemp, 
              ((distERev1 == -1 || distSR1 == -1) ? -1 : distERev1 + distSR1)});

            distSFdTemp = minPos({distSFdTemp,
                ((distSFd2 == -1 || distSL2 == -1) ? -1 : distSFd2 + distSL2)});
 
            distERevTemp = minPos({distERevTemp,
             ((distERev1 == -1 || distER1 == -1) ? -1 : distERev1 + distER1) });

            distEFdTemp = minPos({distEFdTemp, 
               ((distSFd2 == -1 || distEL2 == -1) ? -1 : distSFd2 + distEL2) });
        }
        distSRev = distSRevTemp;
        distSFd = distSFdTemp;
        distERev = distERevTemp;
        distEFd = distEFdTemp;
  
#ifdef indexTraverse 
cerr << "DISTANCES AFTER SNARL: " << distSRev << " " << distSFd << " " << distERev << " " << distEFd << endl;
#endif
        node1 = snarlDists.snarlStart;

        node2 = node1; 

        if (sm->in_nontrivial_chain(snarl)) {
            //Loop distance through a chain

            node2 = snarlDists.snarlEnd;

            const Chain* chain = sm->chain_of(snarl);
           
            Visit startVisit = get_start_of(*chain);
            Visit endVisit = get_end_of(*chain);

            pair<id_t, bool> chainStartIn( startVisit.node_id(),
                                           startVisit.backward() );
            pair<id_t, bool> chainEndIn( endVisit.node_id(), 
                                           !endVisit.backward());

            ChainIndex& chainDists = chainDistances.at(chainStartIn.first);

            pair<id_t, bool> snarlStart = snarlDists.snarlStart;
            pair<id_t, bool> snarlEnd = snarlDists.snarlEnd;

            int64_t loopChain = chainDists.chainDistanceShort(graph,
                                    snarlEnd, snarlStart, snarl, snarl);

            int64_t loopL = chainDists.chainDistanceShort(graph, 
                               make_pair(snarlStart.first, !snarlStart.second),
                               snarlStart, snarl, snarl);
            int64_t loopR = chainDists.chainDistanceShort(graph, snarlEnd, 
                      make_pair(snarlEnd.first, !snarlEnd.second), snarl, snarl); 

#ifdef indexTraverse 
cerr << "LOOP DISTANCES IN CHAIN " << chainStartIn.first << " from node " << snarl->start().node_id() << " to " << snarl->end().node_id() << " : " << loopChain  <<  " " << loopL << " " << loopR << endl;
#endif

            int64_t loop1 = loopChain == -1 || distSRev == -1 || distEFd == -1 ? -1 :
                                     loopChain + distSRev + distEFd + nodeLen;
            int64_t loop2 = loopChain == -1 || distSFd == -1 || distERev == -1 ? -1 :
                                loopChain + distSFd + distERev + nodeLen;
            loopL = loopL == -1 || distSFd == -1 || distSRev == -1 ? -1 : 
                                           loopL + distSFd + distSRev + nodeLen;
            loopR = loopR == -1 || distEFd == -1 || distERev == -1 ? -1 : 
                                           loopR + distEFd + distERev + nodeLen;
            minLoop = minPos({minLoop, loop1, loop2, loopL, loopR });

#ifdef indexTraverse 
cerr << "   CHAIN LOOPS " << chainStartIn.first << " : " << loop1  <<  " " << loop2 << " " << loopL << " " << loopR << endl;
#endif

            pair<id_t, bool> node2Rev = make_pair(node2.first, !node2.second);
            const Snarl* startSnarl = sm->into_which_snarl(chainStartIn.first,
                                                          chainStartIn.second);
            const Snarl* endSnarl = sm->into_which_snarl(chainEndIn.first,
                                                         chainEndIn.second);
            //Get distance to ends of the chain
            int64_t distSL = chainDists.chainDistance(chainStartIn, node1, 
                                                        startSnarl, snarl);
            int64_t distSR = chainDists.chainDistance(chainStartIn, node2Rev, 
                                                        startSnarl, snarl);
            int64_t distEL = chainDists.chainDistance(chainEndIn, node1, 
                                                        endSnarl, snarl);
            int64_t distER = chainDists.chainDistance(chainEndIn, node2Rev, 
                                                        endSnarl, snarl);
        
       
            int64_t distSRevTemp = minPos({ 
               ((distSRev == -1 || distSL == -1) ? -1 : distSRev + distSL), 
               ((distERev == -1 || distSR == -1) ? -1 : distERev + distSR)});

            int64_t distSFdTemp = minPos({
                ((distSFd == -1 || distSL == -1) ? -1 : distSFd + distSL),
                ((distEFd == -1 || distSR == -1) ? -1 : distEFd + distSR) });

            int64_t distERevTemp = minPos({
                ((distSRev == -1 || distEL == -1) ? -1 : distSRev + distEL),
                ((distERev == -1 || distER == -1) ? -1 : distERev + distER) });

            int64_t distEFdTemp = minPos({ 
                ((distSFd == -1 || distEL == -1) ? -1 : distSFd + distEL), 
                ((distEFd == -1 || distER == -1) ? -1 : distEFd + distER) });

            distSRev = distSRevTemp;
            distSFd = distSFdTemp;
            distERev = distERevTemp;
            distEFd = distEFdTemp;     

#ifdef indexTraverse 

cerr << "DISTANCES chain? : " << distSL << " " << distSR << " " << distEL << " " << distER << endl;       
cerr << "DISTANCES TO ENDS OF CHAIN: " << distSRev << " " << distSFd << " " << distERev << " " << distEFd << endl;
#endif
            bool rev1 = node1.first == chainEndIn.first ? 
             !get_start_of(*chain).backward() : get_start_of(*chain).backward();
            node1 = make_pair(chainStartIn.first, rev1);
            node2 = node1;

  
        }
        snarl = sm->parent_of(snarl);

    }

    return minLoop;
    
}





//////////////////////Methods for testing
//


pair<int64_t, int64_t> DistanceIndex::sizeOf() {
    //Estimate of the size of the object in memory
   
    int64_t totalMin = 0;
   
    int64_t numSnarls = snarlDistances.size();

    int64_t snarlDists = 0;
    int64_t snarlNodes = 0; //# node ids + direction

    for (auto x : snarlDistances) {
        //Add size of each SnarlIndex object
        SnarlIndex sd = x.second;
        int64_t numNodes = sd.visitToIndex.size();

        snarlNodes += numNodes; 
        snarlDists += ((numNodes + 1) * numNodes) / 2;
  
        totalMin += numNodes * 17; //Add all elements in visitToIndex
        totalMin += sd.distances.capacity() / 8;
        
        totalMin += 3 * sizeof(pair<id_t, bool>);
        totalMin += sizeof(hash_map<pair<id_t, bool>, int64_t>);
        

    }
    
    int64_t chainDists = 0;
    int64_t chainNodes = 0;

    int64_t numChains = chainDistances.size();
  
    for (auto x : chainDistances) {
        ChainIndex cd = x.second;
        int64_t numNodes = cd.snarlToIndex.size();
        
        chainDists += numNodes*3;
        chainNodes += numNodes;

        totalMin += numNodes * 16; //Add all elements in snarlToIndex
        totalMin += cd.prefixSum.capacity() / 8;
        totalMin += cd.loopFd.capacity() / 8;
        totalMin += cd.loopRev.capacity() / 8;
        totalMin += sizeof(id_t) + sizeof(hash_map<id_t, int64_t>);
    }
 
    totalMin += nodeToSnarl.size() * 8;//TODO: ???
  
    int64_t totalMax = 0;
    totalMax += maxIndex.minDistances.capacity()/8;
    totalMax += maxIndex.maxDistances.capacity()/8;
    totalMax += maxIndex.nodeToComponent.capacity()/8;


    cerr << numSnarls << " snarls containing " << snarlNodes << " nodes" << endl;
    cerr << numChains << " chains containing " << chainNodes << " nodes" << endl;
    cerr << "Total for min index: " << totalMin << " bytes??" << endl;
    cerr << "Total for max index: " << totalMax << " bytes??" << endl;
    return make_pair(totalMin, totalMax); 
    

}


int64_t DistanceIndex::checkChainDist(id_t snarl, size_t index) {
    return chainDistances.at(snarl).prefixSum[index] - 1;
}
int64_t DistanceIndex::checkChainLoopFd(id_t snarl, size_t index) {
    return chainDistances.at(snarl).loopFd[index] - 1;
}
int64_t DistanceIndex::checkChainLoopRev(id_t snarl, size_t index) {
    return chainDistances.at(snarl).loopRev[index] - 1;
}
}
