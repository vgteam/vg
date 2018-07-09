//#define indexTraverse
//#define printDistances

#include "distance.hpp"

using namespace std;
namespace vg {

DistanceIndex::DistanceIndex::DistanceIndex(VG* vg, SnarlManager* snarlManager){
    /*Constructor for the distance index given a VG and snarl manager
    */
    graph = vg;
    sm = snarlManager;
    #ifdef indexTraverse
        cerr << endl << "Creating distance index"<< endl;
    #endif
    const vector<const Snarl*> topSnarls = sm->top_level_snarls();
    calculateIndex(&topSnarls);
};


int64_t DistanceIndex::calculateIndex(const Chain* chain) {
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
    Node* firstNode = graph->get_node(get_start_of(*chain).node_id());

    chainPrefixSum.push_back(firstNode->sequence().size());
    unordered_map<id_t, size_t> snarlToIndex;
    snarlToIndex[get_start_of(*chain).node_id()] = 0;

    ChainIterator chainEnd = chain_end(*chain);
    for (ChainIterator c = chain_begin(*chain); c != chainEnd; ++c) {
        //for each snarl in the chain 

        const Snarl* snarl = c->first;
        bool snarlRevInChain = c->second;

        id_t snarlStartID = snarl->start().node_id();
        bool snarlStartRev = snarl->start().backward(); //into snarl
        id_t snarlEndID = snarl->end().node_id();
        bool snarlEndRev = snarl->end().backward();   //pointing out

        if (snarlToIndex.find(snarlEndID) == snarlToIndex.end()){
//TODO Chain shouldn't loop
            //Store the index of the start of the snarl only if it hasn't
            //already been seen (if the chain loops)
            snarlToIndex[snarlEndID] = snarlToIndex.size()-1;
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
        snarlIndex.insert(make_pair( make_pair(snarlStartID,snarlStartRev),
            SnarlDistances(allNodes, make_pair(snarlStartID,snarlStartRev),
                 make_pair(snarlEndID,  snarlEndRev))));

        SnarlDistances& sd = snarlIndex.at(make_pair(snarlStartID, snarlStartRev));

        #ifdef indexTraverse
            cerr << "Snarl  at " << snarl->start().node_id() << endl;
            cerr << "    Contains nodes : ";
            for (pair<id_t, bool> node: allNodes) {
                cerr << node.first << " "; 
            }
            cerr << endl;
        #endif

        for (pair<id_t, bool> startID : allNodes){
            //Use each node in the snarl as start of djikstra search

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

            while (reachable.size() > 0) {
                pair<pair<id_t, bool>, int64_t> next = reachable.top();
                reachable.pop();
                pair<id_t, bool> currID= next.first;
                int64_t currDist = next.second;
                if ( sd.snarlDistance(startID, currID) == -1) {
                    //If node has not already been found:

                    //Save distance from start to current node 
                    if (currDist == 0 && firstLoop) { 
                        //if current node is the start node
                        sd.insertDistance(startID, currID, -1);
                    } else {
                        sd.insertDistance(startID, currID, currDist);
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
                            unordered_map<id_t, ChainDistances>::iterator 
                                 chainDists = chainIndex.find( get_start_of(
                                     *currChain).node_id());

                            if (chainDists != chainIndex.end()) {
                                //Length of chain has already been found
                               
                                //Get the length of the node (chain)
                                nodeLen = chainDists->second.chainLength();
       
                                //Get loop dist- enter and exit chain at same side
                                if (get_start_of(
                                     *currChain).backward() == currID.second) {
                                    //If traversing snarl forward in chain

                                    loopDist = chainDists->second.loopFd[0];

                                    if (loopDist != -1) {
                                       loopDist = loopDist  + graph->get_node(
                                       currID.first)->sequence().size();
                                     }

                                } else {
                                    loopDist = chainDists->second.loopRev[
                                           chainDists->second.loopRev.size()-1];

                                    if (loopDist != -1) { 
                                        loopDist = loopDist + 
                                           graph->get_node(get_end_of(*currChain
                                             ).node_id())->sequence().size();
                                     }
                                }

                            } else {//haven't recursed on this chain yet
                                #ifdef indexTraverse
                                    cerr << " recurse" << endl;
                                #endif
                                nodeLen = calculateIndex(currChain);

                                ChainDistances& currChainDists = 
                                         chainIndex.at( get_start_of(
                                                    *currChain).node_id());
                                if (get_start_of(
                                     *currChain).backward() == currID.second) {
                                    //If traversing snarl forward in chain

                                    loopDist = currChainDists.loopFd[0];

                                    if (loopDist != -1) {
                                        loopDist = loopDist + graph->get_node(
                                       currID.first)->sequence().size();
                                    }
                                } else {

                                    loopDist = currChainDists.loopRev[
                                            currChainDists.loopRev.size()-1];

                                    if (loopDist != -1) { 
                                        loopDist = loopDist +
                                           graph->get_node(get_end_of(*currChain
                                             ).node_id())->sequence().size();
                                     }
                                } 
                            }
                        } else {//Snarl

                            id_t snarlID = currSnarl->start().node_id();
                            bool snarlRev = currSnarl->start().backward();
                            id_t endID = currSnarl->end().node_id();
                            bool endRev = currSnarl->end().backward();
  
                            auto snarlDists = snarlIndex.find(make_pair(
                                        snarlID, snarlRev));

                            if (snarlDists != snarlIndex.end()) {//Already found
                                nodeLen = snarlDists->second.snarlLength();

                                //Find the distance to enter and exit snarl
                                //at the same side
                                if (currID.second == snarlRev) { 
                                    //If traversing snarl forward
                                    loopDist = snarlDists->second.snarlDistance(
                                           make_pair(snarlID, snarlRev), 
                                           make_pair(snarlID, !snarlRev));

                                     if (loopDist != -1) { 
                                         loopDist = loopDist+ graph->get_node(
                                           currID.first)->sequence().size();
                                     }
                                } else {
                                    loopDist = snarlDists->second.snarlDistance(
                                           make_pair(endID, !endRev), 
                                           make_pair(endID, endRev));

                                     if (loopDist != -1) { 
                                         loopDist = loopDist + 
                                              graph->get_node(currSnarl->end().
                                                node_id())->sequence().size();
                                     }
                                }
                            } else {//Haven't recursed on snarl yet
                                #ifdef indexTraverse
                                    cerr << " recurse" << endl;
                                #endif
                                
                                //Create chain to recurse on and recurse
                                Chain currChain;
                                currChain.push_back(currSnarl);
                                calculateIndex(&currChain);

                                SnarlDistances& currSnarlDists = 
                                     snarlIndex.at(make_pair(snarlID,snarlRev));

                                nodeLen = currSnarlDists.snarlLength(); 

                                //Find the distance to enter and exit snarl
                                //at the same side
                                if (currID.second == snarlRev) {

                                    loopDist = currSnarlDists.snarlDistance(
                                       make_pair(snarlID, snarlRev), 
                                       make_pair(snarlID, !snarlRev));

                                     if (loopDist != -1) { 
                                         loopDist = loopDist + graph->get_node(
                                       currID.first)->sequence().size();
                                     }

                                 } else {

                                     loopDist = currSnarlDists.snarlDistance(
                                        make_pair(endID, !endRev), 
                                        make_pair(endID, endRev));

                                     if (loopDist != -1) { 
                                         loopDist = loopDist+ 
                                             graph->get_node(currSnarl->end().
                                                node_id())->sequence().size(); 
                                      }
                                 }
                            }
                                        
                        }
                    } else { //Node is just a node
                        Node* n = graph->get_node(currID.first);
                        nodeLen = n->sequence().size();
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
                       reachable.push(make_pair(node, currDist + nodeLen));
                       
                         #ifdef indexTraverse
                             cerr << ng.get_id(h) << " " << ng.get_is_reverse(h) << ", ";
                         #endif
                         return true;
                    };
                    //Add reachable nodes to priority queue for unary snarl that doesn't loop - 0 distance
                    auto addHandle0 = [&](const handle_t& h)-> bool {
                         pair<id_t, bool> node = make_pair(
                                    ng.get_id(h), ng.get_is_reverse(h));
                       reachable.push(make_pair(node, 0));
                       
                         #ifdef indexTraverse
                             cerr << ng.get_id(h) << " " << ng.get_is_reverse(h) << ", ";
                         #endif
                         return true;
                    };


                    #ifdef indexTraverse
                         cerr << "    From start node " << startID.first << " " << startID.second<< " in snarl " << snarl->start().node_id() << " at " << ng.get_id(currHandle) << " " << ng.get_is_reverse(currHandle) << endl; 
                         cerr << "        Adding next nodes:  ";
                    #endif

                    if (nodeLen != -1) {

                        ng.follow_edges(currHandle, false, addHandle);
                    } else if (firstLoop) {
//TODO: This might be bad
                        //If the nodeLen is -1 then node is a unary snarl that doesn't have a path from start to end. If this is the start of the distance calculation then add subsequent nodes assuming that the node length was 0
                        ng.follow_edges(currHandle, false, addHandle0);
                    } 
                        

                    //Add edges between the boundary nodes that are not in 
                    //the net graph

                    if (currID.first == snarlStartID &&
                        currID.second != snarlStartRev ) {
                           
                        //If currently leaving start of snarl
                        NodeSide startSide = NodeSide(snarlStartID,
                                                      snarlStartRev);
                        NodeSide endSide = NodeSide(snarlEndID, !snarlEndRev);

                        if (graph->has_edge(startSide,startSide)) {
                             pair<id_t, bool> node = make_pair(
                                               snarlStartID, snarlStartRev);
                           reachable.push(make_pair(node, currDist + nodeLen));
                        } 
                        if (graph->has_edge(startSide, endSide)) {
                             pair<id_t, bool> node = make_pair(
                                               snarlEndID, !snarlEndRev);
                           reachable.push(make_pair(node, currDist + nodeLen));
                        }

                    } else if ( currID.first == snarlEndID &&
                                 currID.second == snarlEndRev ) {
                      
                        //If currently leaving end of snarl

                        NodeSide startSide = NodeSide(snarlStartID,
                                                      snarlStartRev);
                        NodeSide endSide = NodeSide(snarlEndID, !snarlEndRev);

                        if (graph->has_edge(endSide, endSide)) {
                             pair<id_t, bool> node = make_pair(
                                               snarlEndID, !snarlEndRev);
                           reachable.push(make_pair(node, currDist + nodeLen));
                        } 
                        if (graph->has_edge(startSide, endSide)) {
                             pair<id_t, bool> node = make_pair(
                                               snarlStartID, snarlStartRev);
                           reachable.push(make_pair(node, currDist + nodeLen));
                        }
                    }
                    #ifdef indexTraverse
                         cerr << "    prev dist: " << currDist << "+ new dist " << nodeLen << endl;
                    #endif
                } 
                firstLoop = false;
            }//End while loop
        }//End for loop over starting node/directions in a snarl

        #ifdef indexTraverse
            cerr << "End snarl " << snarl->start().node_id() << endl;
        #endif

        /*Add to prefix sum the distance to the beginning and end of the last
            node in the current snarl
        */

        int64_t dist;
        if (snarlRevInChain) {
            dist = sd.snarlDistance(
                 make_pair (snarlEndID, !snarlEndRev) , 
                 make_pair (snarlStartID, !snarlStartRev));

            chainPrefixSum.push_back(chainPrefixSum[chainPrefixSum.size()-2]+
                       dist);
            chainPrefixSum.push_back(chainPrefixSum[chainPrefixSum.size()-1] + 
                    graph->get_node(snarlStartID)->sequence().size());
        
        } else { 
            dist = sd.snarlDistance(
                  make_pair (snarlStartID, snarlStartRev),
                  make_pair (snarlEndID, snarlEndRev) );

            chainPrefixSum.push_back(chainPrefixSum[chainPrefixSum.size()-2]+
                       dist);
            chainPrefixSum.push_back(chainPrefixSum[chainPrefixSum.size()-1] + 
                    graph->get_node(snarlEndID)->sequence().size());
       }
        
         //length of snarl
        if (dist == -1) {
            sd.length = -1;
        } else {
            sd.length = dist + graph->get_node(snarlEndID)->sequence().size();
        }


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
        auto& sd =  snarlIndex.at(make_pair(snarlStartID, snarlStartRev)); 
        //Add reverse loop distances- from start node rev to start node fd



        if ( chainLoopRev.size() == 0) {
           int64_t firstRevDist;
           if (snarlRevInChain){ 
                //If this is the first snarl in the chain 
                firstRevDist = sd.snarlDistance(
                          make_pair(snarlEndID, snarlEndRev),
                          make_pair(snarlEndID, !snarlEndRev));
            } else {
                firstRevDist = sd.snarlDistance(
                   make_pair (snarlStartID, !snarlStartRev),
                   make_pair (snarlStartID, snarlStartRev));
            }


            if (snarlToIndex.size() == (chainPrefixSum.size()/2) -1) {
                //If the chain loops, might need distance from last snarl
                ChainIterator chainEndR = chain_rbegin(*chain);
                const Snarl* lastSnarl = chainEndR->first;
                bool lastRev = chainEndR->second;
                    
                id_t lastStartID = lastSnarl->start().node_id();
                bool lastStartRev = lastSnarl->start().backward();
                id_t lastEndID = lastSnarl->end().node_id();
                bool lastEndRev = lastSnarl->end().backward();
                auto& sdLast =  snarlIndex.at(make_pair(
                                              lastStartID, lastStartRev)); 

                if (lastRev) {
                    firstRevDist = minPos({firstRevDist, 
                             sdLast.snarlDistance(
                                  make_pair(lastStartID, lastStartRev), 
                                   make_pair(lastStartID, !lastStartRev)) });
                  
              
                } else { 
                    firstRevDist = minPos({firstRevDist, 
                             sdLast.snarlDistance(
                                  make_pair(lastEndID, !lastEndRev), 
                                   make_pair(lastEndID, lastEndRev)) });

                }
            }
            chainLoopRev.push_back(firstRevDist);
        }
        int64_t revLoopDist;

        if ( snarlRevInChain ) {

            revLoopDist = sd.snarlDistance(
              make_pair (snarlStartID, snarlStartRev),
             make_pair (snarlStartID, !snarlStartRev));
        } else {
            revLoopDist = sd.snarlDistance(
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
                 sd.snarlDistance(
                make_pair (snarlEndID, !snarlEndRev),
               make_pair (snarlStartID, !snarlStartRev))
            + 
                 sd.snarlDistance(
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
        auto& sd =  snarlIndex.at(make_pair(snarlStartID, snarlStartRev)); 


        if (c == chain_rbegin(*chain)) {
            //If this is the last snarl in the chain, push loop for last node

            int64_t loopDistLast; 
            if (snarlRevInChain) {
       
                loopDistLast = sd.snarlDistance(
                         make_pair(snarlStartID, !snarlStartRev), 
                         make_pair(snarlStartID, snarlStartRev));
            } else {

                loopDistLast = sd.snarlDistance(
                         make_pair(snarlEndID, snarlEndRev), 
                         make_pair(snarlEndID, !snarlEndRev));
            }

            if (snarlToIndex.size() == (chainPrefixSum.size()/2) -1) {
                //If the chain loops, might need distance from first snarl
                ChainIterator chainStart = chain_begin(*chain);
                const Snarl* firstSnarl = chainStart->first;
                bool firstSnarlRev = chainStart->second;
                    
                id_t firstStartID = firstSnarl->start().node_id();
                bool firstStartRev = firstSnarl->start().backward();
                id_t firstEndID = firstSnarl->end().node_id();
                bool firstEndRev = firstSnarl->end().backward();
                auto& sdFirst =  snarlIndex.at(make_pair(
                               firstStartID, firstStartRev)); 
                if (firstSnarlRev) {
                    loopDistLast = minPos({loopDistLast, 
                             sdFirst.snarlDistance(
                                  make_pair(firstEndID, !firstEndRev), 
                                   make_pair(firstEndID, firstEndRev)) });
                } else {
                    loopDistLast = minPos({loopDistLast, 
                             sdFirst.snarlDistance(
                                   make_pair(firstStartID, firstStartRev), 
                                   make_pair(firstStartID, !firstStartRev)) });

                }
              
            }
            chainLoopFd.push_back(loopDistLast);
        }

        int64_t fdLoopDist;


        if (snarlRevInChain) {
            //If the snarl is reversed in the chain
            fdLoopDist = sd.snarlDistance(
                  make_pair (snarlEndID, !snarlEndRev),
                  make_pair (snarlEndID, snarlEndRev));
        } else {
            fdLoopDist = sd.snarlDistance(
                  make_pair (snarlStartID, snarlStartRev),
                  make_pair (snarlStartID, !snarlStartRev));
        }

        int64_t lastLoop = chainLoopFd.back();

        if (lastLoop == -1) {

            chainLoopFd.push_back(fdLoopDist);

        } else {
        //push dist to end of snarl + loop dist + dist to start of snarl 

            int64_t loopDistance = minPos({fdLoopDist, lastLoop + 
                    sd.snarlDistance(make_pair(snarlEndID, !snarlEndRev),
                                     make_pair(snarlStartID, !snarlStartRev)) + 
                    sd.snarlDistance(make_pair(snarlStartID, snarlStartRev), 
                                     make_pair(snarlEndID, snarlEndRev))});
            chainLoopFd.push_back(loopDistance); 
        }           
      
    }
    reverse(chainLoopFd.begin(), chainLoopFd.end()); 
 
/*TODO

    if (snarlToIndex.size() == (chainPrefixSum.size()/2) -1) {
        //If the chain loops
        chainLoopFd[chainLoopFd.size()-1] = chainLoopFd[0];
        chainLoopRev[0] = chainLoopRev[chainLoopRev.size()-1];
    } 
*/
    if (chainPrefixSum.size() > 4) { //If chain and not just one snarl
        chainIndex.insert(make_pair(get_start_of(*chain).node_id(), 
                 ChainDistances(snarlToIndex, chainPrefixSum, chainLoopFd,
                                                               chainLoopRev)));
    }
    return chainPrefixSum.back();//return length of entire chain
};




//////////////////    Calculate distances

int64_t DistanceIndex::distance(const Snarl* snarl1, const Snarl* snarl2, 
                                                  pos_t& pos1, pos_t& pos2) {
    /*Find the distance between two positions
      pos1 and pos2 must be on nodes contained in snarl1/snarl2 */
    
    int64_t shortestDistance = -1; 

    if (get_id(pos1) == get_id(pos2)) { //if positions are on the same node
        int64_t nodeSize = graph->get_node(get_id(pos1))->sequence().size();
        int64_t offset1;
        if (is_rev(pos1)) {
            offset1 = nodeSize -get_offset(pos1) - 1;//Len of node - offset 
        } else {
            offset1 = get_offset(pos1);
        }

        int64_t offset2;
        if (is_rev(pos2)) {
            offset2 = nodeSize - get_offset(pos2) - 1;
        } else {
            offset2 = get_offset(pos2);
        }

        if (graph->has_edge(node_start(get_id(pos1)), node_end(get_id(pos1)))){
            //If there is an edge from start to end of node

            shortestDistance = min(   abs(offset1-offset2)+1,
                          nodeSize - abs(offset1-offset2) + 1  ); 

        } else {

            shortestDistance = abs(offset1-offset2)+1; //+1 to be consistent

        }
    }

    id_t nodeID1 = get_id(pos1);
    bool nodeRev1 = false;
    id_t nodeID2 = get_id(pos2); 
    bool nodeRev2 = false;


    const Snarl* commonAncestor = NULL; 


#ifdef printDistances
    cerr << endl << "Start distance calculation from " << nodeID1 << "->" <<
         nodeID2 << endl;

    cerr << "Shortes distance within same node: " << shortestDistance<<  endl;

    cerr << "Find common ancestor" << endl;
#endif


    //// Find common ancestor of the two snarls
    unordered_set<const Snarl*> ancestors1;//set of all ancestor snarls of node1
    const Snarl* ancestor1 = snarl1;

#ifdef printDistances
    cerr << "Ancestors of 1: ";
#endif


    while (ancestor1 != NULL) {
#ifdef printDistances
        cerr << ancestor1->start().node_id() << " ";
#endif
        ancestors1.emplace(ancestor1);
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


        if (ancestors1.count(ancestor2) > 0) { 
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
    distToCommonAncestor(snarl1, commonAncestor, pos1);
    pair<int64_t, int64_t> temp1 = p1.first; 
    snarl1 = p1.second;
    if (snarl1 != commonAncestor) {
        nodeID1 = snarl1->start().node_id();
        nodeRev1 = snarl1->start().backward();
    }
    int64_t distL1 = temp1.first; int64_t distR1 = temp1.second;
    
    pair<pair<int64_t, int64_t>, const Snarl*> p2 = 
    distToCommonAncestor(snarl2, commonAncestor, pos2);
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
        bool chainStartRev = get_start_of(*chain).backward(); 

        ChainDistances& chainDists = chainIndex.at( chainStartID); 
        bool snarlRev1 = chainDists.isReverse(snarl1, sm); 
        bool snarlRev2 = chainDists.isReverse(snarl2, sm);

        //Distance from left of s1 (reverse), left of s2 (forward)
        int64_t d1 = chainDists.chainDistanceShort(graph,
               make_pair(nodeID1, !snarlRev1), 
               make_pair(nodeID2, snarlRev2));
        d1 = (distL1 == -1 || distL2 == -1 || d1 == -1) ? -1 : 
                                distL1 + distL2 + d1; 

        //Distance from left of s1 (reverse) to right of s2 (reverse)
        int64_t d2 = chainDists.chainDistanceShort(graph,
               make_pair(nodeID1, !snarlRev1), 
               make_pair(endID2, !snarlRev2));
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
                  make_pair(endID1, snarlRev1), 
                  make_pair(nodeID2, snarlRev2));
        if (endID1 == nodeID2) {
            d3 = (distR1 == -1 || distL2 == -1 || d3 == -1) ? -1 : 
                                   distR1 + distL2 - d3; 
        } else {
            d3 = (distR1 == -1 || distL2 == -1 || d3 == -1) ? -1 : 
                                   distR1 + distL2 + d3; 
        }

        //Distance from right of s1 (fd) to right of s2 (rev)
        int64_t d4 =  chainDists.chainDistanceShort(graph,
              make_pair(endID1, snarlRev1),
              make_pair(endID2, !snarlRev2));
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
        return minPos({chainDist, shortestDistance});
    }  

    //Get dist from pos1 to ends of its chain 
    if (snarl1 != commonAncestor && sm->in_nontrivial_chain(snarl1)) {
        const Chain* chain = sm->chain_of(snarl1);

        id_t chainStartID = get_start_of(*chain).node_id();
        ChainDistances& chainDists = chainIndex.at( chainStartID);
        pair<int64_t, int64_t> endDists = chainDists.distToEnds(
                           make_pair(nodeID1, nodeRev1 != get_start_of(*chain).backward()), distL1, distR1);

        distL1 = endDists.first;
        distR1 = endDists.second;

        nodeID1 =  chainStartID;
        nodeRev1 = get_start_of(*chain).backward();
    }
    //Get dist from pos2 to ends of its chain 
    if (snarl2 != commonAncestor && sm->in_nontrivial_chain(snarl2)) {
        const Chain* chain = sm->chain_of(snarl2);
        id_t chainStartID = get_start_of(*chain).node_id();
        ChainDistances& chainDists = chainIndex.at( chainStartID);

        pair<int64_t, int64_t> endDists = chainDists.distToEnds(
                       make_pair(nodeID2, nodeRev2 != get_start_of(*chain).backward()), distL2, distR2);    
        distL2 = endDists.first;
        distR2 = endDists.second;

        nodeID2 = chainStartID; 
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
    SnarlDistances& snarlDists = snarlIndex.at(make_pair(
                           commonAncestor->start().node_id(),
                                       commonAncestor->start().backward()));

    int64_t d1 = snarlDists.snarlDistanceShort(graph, &ng,
                make_pair(nodeID1, nodeRev1), make_pair(nodeID2, nodeRev2));
    d1 = (distR1 == -1 || distL2 == -1 || d1 == -1) ? -1 : 
                                                  distR1 + distL2 + d1; 

    int64_t d2 = snarlDists.snarlDistanceShort(graph, &ng,
                   make_pair(nodeID1, nodeRev1), make_pair(nodeID2, !nodeRev2));

    d2 = (distR1 == -1 || distR2 == -1 || d2 == -1) ? -1 :
                                                   distR1 + distR2 + d2;
    int64_t d3 = snarlDists.snarlDistanceShort(graph, &ng,
                make_pair(nodeID1, !nodeRev1), make_pair(nodeID2, nodeRev2));
    d3 = (distL1 == -1 || distL2 == -1 || d3 == -1) ? -1 : 
                                              distL1 + distL2 + d3; 
    int64_t d4 = snarlDists.snarlDistanceShort(graph, &ng,
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
    pair<int64_t, int64_t> endDists = snarlDists.distToEnds(nodeID1, nodeRev1, distL1, 
                                                                        distR1);
    distL1 = endDists.first;
    distR1 = endDists.second;

    endDists = snarlDists.distToEnds(nodeID2, nodeRev2, distL2, distR2);
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
            ChainDistances& chainDists = chainIndex.at(
                                            get_start_of(*currChain).node_id());
            //TODO: If too slow, this is done twice        
            bool chainStartRev = get_start_of(*currChain).backward(); 

            bool snarlRev = chainDists.isReverse(currSnarl, sm); 

            //Distance from start (reverse) to start (forward)
            int64_t d1 = chainDists.chainDistanceShort(graph,
               make_pair(startID, !snarlRev), 
                   make_pair(startID, snarlRev));
            d1 = (distL1 == -1 || distL2 == -1 || d1 == -1) ? -1 : 
                                distL1 + distL2 + d1; 

            //Distance from start (reverse) to end (reverse)
            int64_t d = chainDists.chainDistanceShort(graph,
               make_pair(startID, !snarlRev), 
               make_pair(endID, !snarlRev));

            d2 = (distL1 == -1 || distR2 == -1 || d == -1) ? -1 : 
                                   distL1 + distR2 + d; 
            d3 = (distR1 == -1 || distL2 == -1 || d == -1) ? -1 : 
                                   distR1 + distL2 + d; 

           //Distance from end (fd) to end (rev)
            int64_t d4 =  chainDists.chainDistanceShort(graph,
                  make_pair(endID, snarlRev),
                  make_pair(endID, !snarlRev));
            d4 = (distR1 == -1 || distR2 == -1 || d4 == -1) ? -1 : 
                                   distR1 + distR2 + d4;
   
                  
            shortestDistance = minPos({shortestDistance, d1, d2, d3, d4});
           
            //Find distances to ends of the current chain

            pair<int64_t, int64_t> chainEnd1 = chainDists.distToEnds(
                      make_pair(startID, startRev != 
                          get_start_of(*currChain).backward()), distL1, distR1);
            distL1 = chainEnd1.first; distR1 = chainEnd1.second;

            pair<int64_t, int64_t> chainEnd2 = chainDists.distToEnds(
                          make_pair(startID, startRev!=
                          get_start_of(*currChain).backward()), distL2, distR2);
            distL2 = chainEnd2.first; distR2 = chainEnd2.second;


            startID = get_start_of(*currChain).node_id();
            startRev = get_start_of(*currChain).backward();
            endID = get_end_of(*currChain).node_id();
            endRev = get_end_of(*currChain).backward();

#ifdef printDistances
    cerr << "At chain " << startID << " dists to ends: " << distL1 << " " << 
         distR1 << " " << distL2 << " " << distR2 << endl;
    cerr << "distances: "  << d1 << " " << d2 << " " << d3 << " " << d4 << endl;
  << " Shortest distance : "
          << shortestDistance << endl;
#endif  
        }
   
        if (parentSnarl == NULL) {break;} //TODO: Make this cleaner

        SnarlDistances& snarlDists =snarlIndex.at(
                                   make_pair(parentSnarl->start().node_id(),
                                             parentSnarl->start().backward()));


        NetGraph ng = NetGraph(parentSnarl->start(), 
                parentSnarl->end(),sm->chains_of(parentSnarl), graph);

        //Find the shortest distance within the snarl

        //Dist from start to start
        d1 = snarlDists.snarlDistanceShort(graph, &ng, 
                make_pair(startID, !startRev), make_pair(startID, startRev));
        d1 = (distL1 == -1 || distL2 == -1 || d1 == -1) ? -1 : 
                                                   distL1 + distL2 + d1; 
    
        //Dist from end to end
        d2 = snarlDists.snarlDistanceShort(graph, &ng,
                   make_pair(startID, startRev), make_pair(startID, !startRev));

         d2 = (distR1 == -1 || distR2 == -1 || d2 == -1) ? -1 :
                                                   distR1 + distR2 + d2;
         //Dist from start to end
         int64_t dtemp = snarlDists.snarlDistanceShort(graph, &ng,
                make_pair(startID, startRev), make_pair(startID, startRev));
         d3 = (distL1 == -1 || distR2 == -1 || dtemp == -1) ? -1 : 
                                              distL1 + distR2 + dtemp; 
         d4 = (distR1 == -1 || distL2 == -1 || dtemp == -1) ? -1 : 
                                                  distR1 + distL2 + dtemp; 


         shortestDistance =  minPos({d1, d2, d3, d4, shortestDistance});


        //Find the distances to ends of the snarl
        pair<int64_t, int64_t> endDists1 = snarlDists.distToEnds(startID, 
                                                 startRev, distL1, distR1);
        distL1= endDists1.first; distR1= endDists1.second;

        pair<int64_t, int64_t> endDists2 = snarlDists.distToEnds(startID,
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

    return shortestDistance;

};


pair<pair<int64_t, int64_t>, const Snarl*> DistanceIndex::distToCommonAncestor(
          const Snarl* snarl, const Snarl* commonAncestor, pos_t& pos){

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
              << " in snarl " << snarl->start().node_id() << endl;
    
    #endif
    if (is_rev(pos)) {//Get distance to ends of current node
        distL = graph->get_node(get_id(pos))->sequence().size()-offset;
        distR = offset + 1;
    } else {
        distR = graph->get_node(get_id(pos))->sequence().size()-offset;
        distL = offset + 1;
    }

    #ifdef printDistances
        cerr << "start pos: " << get_offset(pos) << "-> start: " << distL << 
               ", end: " << distR << endl;
    #endif
 
    if (snarl == commonAncestor) {
        /*If the node is a node in commonAncestor, return the distances to 
           the ends of the node
        */
        return make_pair(make_pair(distL, distR), snarl);
    }

    id_t startID = snarl->start().node_id(); 
    bool startRev = snarl->start().backward();

    SnarlDistances& snarlDists = snarlIndex.at(make_pair(startID, startRev));
 
    pair<int64_t, int64_t> endDists = snarlDists.distToEnds(nodeID, false, 
                                                                 distL, distR);
    distL = endDists.first;
    distR = endDists.second;

    #ifdef printDistances
    cerr << nodeID << "->" << startID << ": " << distL << ", " << distR << endl;
    #endif

    nodeID = startID;
    bool nodeRev = startRev;

    while (sm->parent_of(snarl) != commonAncestor) {
        int64_t dsl; int64_t dsr; int64_t der; int64_t del;

        if (sm->in_nontrivial_chain(snarl)) {
            //Get distances to ends of chain

            const Chain* chain = sm->chain_of(snarl);
            id_t chainStartID = get_start_of(*chain).node_id();
            ChainDistances& chainDists =  chainIndex.at(chainStartID);

            pair<int64_t, int64_t> endDists = chainDists.distToEnds(
                          make_pair(nodeID, nodeRev != 
                               get_start_of(*chain).backward()), distL, distR);

            distL = endDists.first;   
            distR = endDists.second;

            nodeID = chainStartID; 
            nodeRev = get_start_of(*chain).backward();
    #ifdef printDistances
        cerr << nodeID << "->" << chainStartID << ": " << distL << ", " << distR
                                                                   << endl;
    #endif
        }
     
        //Get distances to ends of parent snarl
        snarl = sm->parent_of(snarl);
        id_t startNodeID = snarl->start().node_id();
        id_t startNodeRev = snarl->start().backward();
        SnarlDistances& snarlDists = 
                            snarlIndex.at(make_pair(startNodeID, startNodeRev));
            
        pair<int64_t, int64_t> endDists = snarlDists.distToEnds(
                                                nodeID, nodeRev, distL, distR);
          
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



void DistanceIndex::printSelf() {
    for (auto snarls : snarlIndex) {
        snarls.second.printSelf();
    }
    for (auto chains : chainIndex) {
        chains.second.printSelf();
    }
}


DistanceIndex::SnarlDistances::SnarlDistances(unordered_set<pair<id_t, bool>>& 
                      allNodes,  pair<id_t, bool> start, pair<id_t,bool> end) {
    /*Constructor for SnarlDistances object that stores distances between
        nodes in a snarl */

    //Assign all nodes+direction in snarl to an index
    size_t snarlIndex = 0;
    for (pair<id_t, bool> node: allNodes) {
        visitToIndex[node] = snarlIndex++;
    }

    int size = visitToIndex.size();
    //Initialize all distances to -1
    distances.assign(size*size, -1); 
    snarlStart = start;
    snarlEnd = end;

}


size_t DistanceIndex::SnarlDistances::index(pair<id_t, bool> start, 
                                            pair<id_t, bool> end) {
    /*Get the index of dist from start to end in a snarl distance matrix
      given the node ids + direction */
    size_t length = visitToIndex.size();
    size_t i1 = visitToIndex.at(start);
    size_t i2 = visitToIndex.at(end);
    return i1 * length + i2;
}

void DistanceIndex::SnarlDistances::insertDistance(pair<id_t, bool> start, 
                                           pair<id_t, bool> end, int64_t dist) {
    //Assign distance between start and end
    size_t i = index(start, end);
    distances.at(i) = dist;
}
   
int64_t DistanceIndex::SnarlDistances::snarlDistance(pair<id_t, bool> start,
                                                      pair<id_t, bool> end) {
    /*Distance between beginnings of two nodes n1 and n2 in snarl
    */
    size_t i = index(start, end);
    return distances.at(i); 
}

int64_t DistanceIndex::SnarlDistances::snarlDistanceShort(VG* graph, 
             NetGraph* ng,pair<id_t, bool> start, pair<id_t, bool> end) {
    /*Distance between end of node n1 and beginning of node n2 in snarl
    */
    
    //Check both directions and pick shortest one
    int64_t dist1 = snarlDistanceShortHelp(graph, ng, start, end);

    int64_t dist2 =  snarlDistanceShortHelp(graph, ng, 
      make_pair(end.first, !end.second), make_pair(start.first, !start.second));

    if (dist1 == -1) {
        return dist2;

    } else if (dist2 == -1){
        return dist1;

    } else {
        return min(dist1, dist2);
    }
}

int64_t DistanceIndex::SnarlDistances::snarlDistanceShortHelp(VG* graph, 
             NetGraph* ng,pair<id_t, bool> start, pair<id_t, bool> end) {
    /*Helper for snarlDistanceShort - does the actual calculation for one 
       direction*/

    size_t i = index(start, end);
    int64_t totalDist = distances.at(i); 

    if (totalDist == -1 ) {//No path between two nodes
        return -1;

    } else {
       /* If there is a path, return the distance between two nodes - length of
         first node */
       int64_t firstDist = nodeLength(graph, ng,start);
                        //Length of first node 
      
       return totalDist - firstDist;
       
    }
}

int64_t DistanceIndex::SnarlDistances::nodeLength(VG* graph, 
           NetGraph* ng, pair<id_t, bool> node){

    //Get the length of the node. 
    //Requires node has an edge to another node in the snarl

    handle_t handle = ng->get_handle(node.first, node.second);

    pair<id_t, bool> next (-1, -1); //index of a node after start
    auto getIndex = [&](const handle_t& h)-> bool {
    next.first = ng->get_id(h);
       next.second = ng->get_is_reverse(h); 
       return false;
    };
    ng->follow_edges(handle, false, getIndex);

    if (next.first == -1) { 
        return nodeLength(graph, ng, make_pair(node.first, !node.second));
    } else {

        size_t j = index(node, next); 
        return distances.at(j);
    }
}

int64_t DistanceIndex::SnarlDistances::snarlLength() {   
    //Return the length of the snarl- dist from beginning of start to end of end
//TODO: Might be worth it to calculate this each time rather than saving the length
    return length; 
 
}

pair<int64_t, int64_t> DistanceIndex::SnarlDistances::distToEnds(id_t node, 
                                     bool rev, int64_t distL, int64_t distR) {
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
    int64_t dsl = snarlDistance(snarlStart, make_pair(node, false)); 

    int64_t dsr = snarlDistance( snarlStart, make_pair(node, true));

    int64_t der = snarlDistance( snarlEndRev, make_pair(node, true));

    int64_t del = snarlDistance( snarlEndRev, make_pair(node, false));

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

    int64_t distStart;
    if (dsl == -1) {distStart = dsr;}
    else if (dsr ==-1) {distStart = dsl;}
    else {distStart = min(dsr, dsl);}

    int64_t distEnd;
    if (del == -1) {distEnd = der;}
    else if (der ==-1) {distEnd = del;}
    else {distEnd = min(der, del);}

    return make_pair(distStart, distEnd);
}

void DistanceIndex::SnarlDistances::printSelf() {
    //Print the nodes contained in SnarlDistance
    cerr << endl;
     
    cerr << "Snarl Distances for snarl starting at " << snarlStart.first;
    if (snarlStart.second) {cerr << " reverse and ending at ";} 
    else                   { cerr << " forward and ending at";}
    cerr << snarlEnd.first;
    if (snarlEnd.second) {cerr << " reverse";} 
    else {cerr << " forward";}
    cerr << "Indices:" << endl;
    
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
   
            size_t i = i1 * length + i2;;
            cerr << distances.at(i)  << "   "; 
        }
        cerr << endl;
    }
    cerr << endl; 
}

//ChainDistance methods
DistanceIndex::ChainDistances::ChainDistances(unordered_map<id_t, size_t> s, 
                 vector<int64_t> p, vector<int64_t> fd, vector<int64_t> rev) {
    
    snarlToIndex = move(s);
    prefixSum = move(p);
    loopFd = move(fd);
    loopRev = (rev);
}

int64_t DistanceIndex::ChainDistances::chainDistance(pair<id_t, bool> start, 
                                                      pair<id_t, bool> end) {
    /*Returns the distance between start of start and start of end in a chain
      Bools are true if traversed reverse relative to the start of the chain
    */
    size_t i1 = snarlToIndex.at(start.first);
    size_t i2 = snarlToIndex.at(end.first); 

    return chainDistanceHelper(make_pair(i1, start.second), 
                               make_pair(i2, end.second));
}

int64_t DistanceIndex::ChainDistances::chainDistanceHelper(
                             pair<size_t, bool> start, pair<size_t, bool> end) {
    /*Return the distance from the index start to end. Same as chainDistance
      but given the index of the node in the chain not node id */

    size_t i1 = start.first;
    size_t i2 = end.first;
    int64_t loopDist = -1;

    if (snarlToIndex.size() == (prefixSum.size()/2) -1 && i1 != i2) {
        //If the chain loops

         size_t size = snarlToIndex.size();
         if (i1 == 0) {
             loopDist = chainDistanceHelper(make_pair(size, start.second), 
                                                end); 
         } else if (i2 == 0) { 
             loopDist = chainDistanceHelper(start, 
                                            make_pair(size, end.second));
         }
    }

    if ((!start.second && !end.second)) {
        //If start and end are facing forward relative to the start of the chain
        if (i1 <= i2) {
            int64_t dNoRev = prefixSum[2*i2] - prefixSum[2*i1]; 
            return minPos({loopDist, dNoRev});
        } else {
            int64_t rev1 = loopFd[i1];
            int64_t rev2 = loopRev[i2];
            return minPos({loopDist, (rev1 == -1 || rev2 == -1) ? -1 : 
                     (prefixSum[2*i1+1] - prefixSum[2*i2+1]) + rev1 + rev2}); 
        }

    } else if (start.second && end.second ){
        //If start and end are both reversed relative to the start of the chain
        if (i1 >= i2) {
            int64_t dNoRev = prefixSum[2*i1+1] - prefixSum[2*i2+1]; 
            return minPos({loopDist, dNoRev});
        } else {
            int64_t rev1 = loopRev[i1];
            int64_t rev2 = loopFd[i2];
            return minPos({loopDist, ((rev1 == -1 || rev2 == -1) ? -1 : 
                     (prefixSum[2*i2] - prefixSum[2*i1]) + rev1 + rev2)}); 
        }
    } else if (!start.second && end.second) {
        //Start is forward, end is reversed
        if (i1 <= i2) {
            int64_t rev = loopFd[i2];
            return minPos({loopDist, ((rev == -1) ? -1 : rev + prefixSum[2*i2]
                                        - prefixSum[2*i1])});
        } else {
            int64_t rev = loopFd[i1];
            return minPos({loopDist, ((rev == -1) ? -1 : 
                               rev + prefixSum[2*i1+1] - prefixSum[2*i2+1])});
        }
        
    } else {
        //start is reverse, end is forward
        if (i1 <= i2) {
            int64_t rev = loopRev[i1];
            return minPos({loopDist, (rev == -1 ? -1 : rev + 
                                        (prefixSum[2*i2] - prefixSum[2*i1]))});

            
        } else {
            int64_t rev = loopRev[i2]; 
            return minPos({loopDist, ((rev == -1) ? -1 : rev +
                (prefixSum[2*i1+1] - prefixSum[2*i2+1]))});
        }
    }
}

int64_t DistanceIndex::ChainDistances::chainDistanceShort(VG* graph, 
                                 pair<id_t, bool> start, pair<id_t, bool> end) {
    /*Distance between end of start node to beginning of end node in chain
      or the distance from the end of the end node to the start of the start 
      node
      If start and end are the same node, then return the length of that node
       because the length is needed for the distance calculation and a negative
       distance would indicate no path.
    */ 
    int64_t d1 = chainDistance(start, end);
    int64_t d2 = chainDistance(make_pair(end.first, !end.second), 
                               make_pair(start.first, !start.second));
    if (start.first == end.first) {
        if (start.second == end.second) {
            //If two positions are on different snarls that share a node
            return graph->get_node(start.first)->sequence().size();
        } 
    }
    if (d1 == -1 && d2 == -1) {
        return -1;
    } else if (d2 == -1) {
        return d1 - graph->get_node(start.first)->sequence().size();
    } else if (d1 == -1) {
        return d2 - graph->get_node(end.first)->sequence().size();
    } else {
        return min(d1 - graph->get_node(start.first)->sequence().size(), 
                   d2 - graph->get_node(end.first)->sequence().size());
    }
}
bool DistanceIndex::ChainDistances::isReverse(const Snarl* snarl, SnarlManager* sm) {
    id_t start = snarl->start().node_id();
    id_t end = snarl->end().node_id();
    if (snarlToIndex.size() == (prefixSum.size()/2) -1 ){
        //If the chain loops
        //TODO there must be a better way of doing this 
        bool startRev = snarl->start().backward();
        bool endRev = snarl->end().backward();
        const Chain* chain = sm->chain_of(snarl);

        ChainIterator chainEnd = chain_end(*chain);
        for (ChainIterator c = chain_begin(*chain); c != chainEnd; ++c) {
            const Snarl* s = c->first; 
            if (start == s->start().node_id() && startRev == s->start().backward()) {
                return c->second;
            }
            
        }
         
cerr << " Shouldn't ever reach here" << endl;
        return true; //Should never reach here
        
    } else {
        //If the start of the snarl has a higher index than end, reversed 
        return snarlToIndex[start] > snarlToIndex[end];
    }   
}
int64_t DistanceIndex::ChainDistances::chainLength() {
    //Get the length of a chain including length of last node
    return prefixSum.back();
}

pair<int64_t, int64_t> DistanceIndex::ChainDistances::distToEnds(
              pair<id_t, bool> start, int64_t distL, int64_t distR) {
    /*Given the distance to either end of snarl starting at start, find the 
       distance to either end of the chain*/
     
    size_t startI = snarlToIndex[start.first];
    pair<size_t, bool> startFd = make_pair(startI, start.second); 
    pair<size_t, bool> startRev;
    if (start.second) {

        if (startI == 0) {
            startRev = make_pair(startI, !start.second);
        } else  { 
            startRev = make_pair(startI-1, !start.second); 
        }

    } else {

        startRev = make_pair(startI+1, !start.second);

    }

    int64_t dsl = chainDistanceHelper(make_pair(0,false), startFd);
    int64_t dsr = chainDistanceHelper(make_pair(0,false), startRev);
    int64_t der = chainDistanceHelper(make_pair(loopFd.size()-1,true),startRev);
    int64_t del = chainDistanceHelper(make_pair(loopFd.size()-1,true), startFd);

    dsl = dsl == -1 || distL == -1? -1 : distL + dsl; 
    dsr =  dsr == -1 || distR == -1? -1 : distR + dsr; 
    der = der == -1 || distR == -1? -1 : distR + der; 
    del = del == -1 || distL == -1? -1 : distL + del; 
 
    int64_t distStart;
    if (dsl == -1) {distStart = dsr;}
    else if (dsr ==-1) {distStart = dsl;}
    else {distStart = min(dsr, dsl);}

    int64_t distEnd;
    if (del == -1) {distEnd = der;}
    else if (der ==-1) {distEnd = del;}
    else {distEnd = min(der, del);}
    return make_pair(distStart, distEnd);
}

void DistanceIndex::ChainDistances::printSelf() {
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

//Methods for testing
int64_t DistanceIndex::checkChainDist(id_t snarl, size_t index) {
    return chainIndex.at(snarl).prefixSum[index];
}
int64_t DistanceIndex::checkChainLoopFd(id_t snarl, size_t index) {
    return chainIndex.at(snarl).loopFd[index];
}
int64_t DistanceIndex::checkChainLoopRev(id_t snarl, size_t index) {
    return chainIndex.at(snarl).loopRev[index];
}
}
