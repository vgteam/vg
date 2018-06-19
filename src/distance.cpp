#define indexTraverse
//#define debug
//#define printDistances

#include "distance.hpp"

using namespace std;
namespace vg {


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
    //return (length + 1) * (length / 2) + (i1+1)* (i1/2) + i2 - 1;
    /*TODO Currently store entire distance matrix since it's simpler
       dist from n1 -> n2 = -n2 -> -n1 - len(n1) + len(n2) */
    return i1 * length + i2;
}
   
int64_t DistanceIndex::SnarlDistances::snarlDistance(pair<id_t, bool> start,
                                                      pair<id_t, bool> end) {
    /*Distance between beginnings of two nodes n1 and n2 in snarl
    */
    size_t i = index(start, end);
    return distances.at(i); 
}

void DistanceIndex::SnarlDistances::insertDistance(pair<id_t, bool> start, 
                                           pair<id_t, bool> end, int64_t dist) {
    //Assign distance between start and end
    size_t i = index(start, end);
    distances.at(i) = dist;
}

int64_t DistanceIndex::SnarlDistances::snarlDistanceShort(NetGraph* graph, 
               SnarlManager* sm, pair<id_t, bool> start, pair<id_t, bool> end) {
    /*Distance between end of node n1 and beginning of node n2 in snarl
    */
    #ifdef debug
    cerr << start.first << start.second << "->" << end.first << end.second 
                  << endl;
    printSelf();
    #endif

    //If same node, return -1 since this distance will never be needed
    if (start.first == end.first) { return -1; }
    size_t i = index(start, end);
    int64_t totalDist = distances.at(i); 

    if (totalDist == -1) { //No path between two nodes
        return -1;

    } else {
       /* If there is a path, return the distance between two nodes - length of
         first node */
       handle_t handle = graph->get_handle(start.first, start.second);
       const Snarl* childSnarl= sm->into_which_snarl(start.first, start.second); 
       if (childSnarl != NULL && start.first != snarlStart.first) {
           handle = graph->get_handle(childSnarl->end().node_id(), 
                                         childSnarl->end().backward());
       }
       vector<const handle_t*> nextNodes;//vector of next nodes
       auto addHandle = [&](const handle_t& h)-> bool {
           nextNodes.push_back(&h);
           return true;
       };
       graph->follow_edges(handle, false, addHandle);
       handle_t nextHandle = *nextNodes[0];
     
       size_t j = index(start, make_pair(graph->get_id(nextHandle), 
                                         graph->get_is_reverse(nextHandle)));

       int64_t firstDist = distances.at(j);//Length of first node 
      
       return totalDist - firstDist;
    }
}

int64_t DistanceIndex::SnarlDistances::snarlLength(VG* graph, 
                                                        const Snarl* snarl) {   
    //Return the length of the snarl- dist from beginning of start to end of end
    return snarlDistance(
         make_pair(snarl->start().node_id(), snarl->start().backward()),
         make_pair(snarl->end().node_id(), snarl->end().backward()))
                                       + 
         graph->get_node(snarl->end().node_id())->sequence().size();
 
}

pair<int64_t, int64_t> DistanceIndex::SnarlDistances::distToEnds(id_t node, 
                                              int64_t distL, int64_t distR) {
    /* Given the distances to either end of a node, find the distances to 
       either end of the snarl
    */

    pair<id_t, bool> snarlEndRev = make_pair(snarlEnd.first, !snarlEnd.second);
    int64_t dsl = snarlDistance(snarlStart, make_pair(node, false)); 

    int64_t dsr = snarlDistance( snarlStart, make_pair(node, true));

    int64_t der = snarlDistance( snarlEndRev, make_pair(node, true));

    int64_t del = snarlDistance( snarlEndRev, make_pair(node, false));


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
    cerr << "Snarl Distances for " << snarlStart.first << endl;
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
      Bools are reverse relative to the start of the chain
    */
    //TODO: assumes that all snarls in a chain face the same direction
    size_t i1 = snarlToIndex.at(start.first);
    size_t i2 = snarlToIndex.at(end.first); 
    return chainDistanceHelper(make_pair(i1, start.second), 
                               make_pair(i2, end.second));
}
int64_t DistanceIndex::ChainDistances::chainDistanceHelper(
                             pair<size_t, bool> start, pair<size_t, bool> end) {
    /*Return the distance from the index start to end. Same as chainDistance
      but given the index not node id */
    size_t i1 = start.first;
    size_t i2 = end.first;
    if ((!start.second && !end.second)) {
        //If start and end are facing forward relative to the start of the chain
        if (i1 <= i2) {
            int64_t dNoRev = prefixSum[2*i2] - prefixSum[2*i1]; 
            return dNoRev;
        } else {
            int64_t rev1 = loopFd[i1];
            int64_t rev2 = loopRev[i2];
            return (rev1 == -1 || rev2 == -1) ? -1 : 
                     (prefixSum[2*i1+1] - prefixSum[2*i2+1]) + rev1 + rev2; 
        }
    } else if (start.second && end.second ){
        //If start and end are both reversed relative to the start of the chain
        if (i1 >= i2) {
            int64_t dNoRev = prefixSum[2*i1+1] - prefixSum[2*i2+1]; 
            return dNoRev;
        } else {
            int64_t rev1 = loopRev[i1];
            int64_t rev2 = loopFd[i2];
            return (rev1 == -1 || rev2 == -1) ? -1 : 
                     (prefixSum[2*i2] - prefixSum[2*i1]) + rev1 + rev2; 
        }
    } else if (!start.second && end.second) {
        //Start is forward, end is reversed
        if (i1 <= i2) {
            int64_t rev = loopFd[i2];
            return (rev == -1) ? -1 : rev + prefixSum[2*i2] - prefixSum[2*i1];
        } else {
            int64_t rev = loopFd[i1];
            return (rev == -1) ? -1 : rev + prefixSum[2*i1+1] - prefixSum[2*i2+1];
        }
        
    } else {
        //start is reverse, end is forward
        if (i1 <= i2) {
            int64_t rev = loopRev[i1];
            return rev == -1 ? -1 : rev + (prefixSum[2*i2] - prefixSum[2*i1]);

            
        } else {
            int64_t rev = loopRev[i2]; 
            return (rev == -1) ? -1 : rev + (prefixSum[2*i1+1] - 
                                                             prefixSum[2*i2+1]);
        }
    }
}

int64_t DistanceIndex::ChainDistances::chainDistanceShort(VG* graph, 
                                 pair<id_t, bool> start, pair<id_t, bool> end) {
    /*Distance between end of start node to beginning of end node in chain
      or the distance from the end of the end node to the start of the start 
      node
      If start and end are the same node, then return the length of that node
       because the length is needed for the distance calculation.
    */ 
    size_t i1 = 2*snarlToIndex.at(start.first);
    size_t i2 = 2*snarlToIndex.at(end.first); 
    int64_t d1 = chainDistance(start, end);
    int64_t d2 = chainDistance(make_pair(end.first, !end.second), 
                               make_pair(start.first, !start.second));
    if (i1 == i2) {
        return graph->get_node(start.first)->sequence().size();
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
    pair<size_t, bool> startRev = make_pair(startI+1, !start.second); 

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

int64_t DistanceIndex::calculateIndex(const Chain* chain) {
    /*Add chain to distance index and recursively add all distances in 
      snarls contained within chain
    */
    auto cmp = [] (pair<pair<id_t, bool>,int64_t> x,
                                           pair<pair<id_t, bool>,int64_t> y) {
        //Comparison function for the priority of a pair of handle, distance
        return (x.second > y.second);
    };

    //initialize chain prefix sum to [0, len of first node in chain]
    vector<int64_t> chainPrefixSum (1,0); 
    vector<int64_t> chainLoopFd; 
    vector<int64_t> chainLoopRev; 
    Node* firstNode = graph->get_node(get_start_of(*chain).node_id());
    chainPrefixSum.push_back(firstNode->sequence().size());
    unordered_map<id_t, size_t> snarlToIndex;
    snarlToIndex[get_start_of(*chain).node_id()] = snarlToIndex.size()-1;

    for (int i = 0; i < chain->size(); i++) {
        //for each snarl in the chain TODO ChainIterator?
        const Snarl* snarl = chain->at(i);
        snarlToIndex[snarl->end().node_id()] = snarlToIndex.size()-1;


        //Create a NetGraph for current snarl
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
        pair<id_t, bool> startPair = 
                  make_pair(snarl->start().node_id(),snarl->start().backward());
        //Create snarl distance obj for current snarl and add to distance index
        snarlIndex.insert(make_pair( startPair,
            SnarlDistances(allNodes, startPair,
                 make_pair(snarl->end().node_id(),  snarl->end().backward()))));

        SnarlDistances& sd = snarlIndex.at(make_pair(snarl->start().node_id(), 
                                                    snarl->start().backward()));

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

            //Priority queue of reachable nodes (handles)
            priority_queue<  pair<pair<id_t, bool>, int64_t>,  
                         vector<pair<pair<id_t, bool>, int64_t>>, 
                                     decltype(cmp)> reachable(cmp);
            reachable.push(make_pair(startID, 0));

            #ifdef indexTraverse
                cerr << "  Start Node: " << startID.first << "," 
                                                   << startID.second << endl;
            #endif
            while (reachable.size() > 0) {
                pair<pair<id_t, bool>, int64_t> next = reachable.top();
                reachable.pop();
                pair<id_t, bool> currID= next.first;
                int64_t currDist = next.second;
                if ( sd.snarlDistance(startID, currID) == -1) {
                    //If node has not already been found:
                    /*
                    TODO: Dist from node to itself is 0, maybe should be -1 
                          but would have to change distToCommonAncestor
                    if (currDist == 0) {
                        sd.insertDistance(startID, currID, -1);
                    } else {
                        sd.insertDistance(startID, currID, currDist);
                    }
                    */

                    //Save distance from start to current node 
                    sd.insertDistance(startID, currID, currDist);
               

                    const handle_t currHandle = ng.get_handle(currID.first, 
                                                               currID.second);
                    handle_t lastHandle;//To find next node in traversal
                #ifdef indexTraverse
                    if (ng.is_child(currHandle)) {
                        cerr << "    At child snarl " << currID.first << ","
                                                             << currID.second;
                    } else {
                        cerr << "    At node " << currID.first << "," 
                                                           << currID.second;
                    }
                #endif           

                    
                    int64_t nodeLen; //length of the current node

                    const Snarl* currSnarl = sm->into_which_snarl(
                                              currID.first, currID.second);
                    //TODO:If is_child, then this shouldn't be null?????
                    // assert(currSnarl != NULL);
                    if (currID.first != snarl->start().node_id() &&
                            currID.first != snarl->end().node_id() &&
                            currSnarl != NULL) {// ng.is_child(currHandle)){
                    //If current node is a child snarl/chain


                        if (sm->in_nontrivial_chain(currSnarl)) {//Chain

                            const Chain* currChain= sm->chain_of(currSnarl);
                            unordered_map<id_t, ChainDistances>::iterator 
                                 chainDists = chainIndex.find( get_start_of(
                                     *currChain).node_id());
                            if (chainDists != chainIndex.end()) {
                                //Length of chain has already been found
                                nodeLen =chainDists->second.chainLength();
                                //last element should be length of chain

                            } else {//haven't recursed on this chain yet
                                #ifdef indexTraverse
                                    cerr << " recurse" << endl;
                                #endif
                                nodeLen = calculateIndex(currChain);
                            }
                            /*handle to last node of chain at the other end
                              to find the next node in curr snarl*/
                            if (currID.first == 
                                            get_start_of(*currChain).node_id()
                                     && currID.second == 
                                         get_start_of(*currChain).backward()  ){

                                lastHandle = ng.get_handle(get_end_of(
                                          *currChain).node_id(), currID.second);
                            } else if (currID.first == 
                                               get_end_of(*currChain).node_id()
                                       && currID.second != 
                                          get_end_of(*currChain).backward()  ) {
 
                                lastHandle = ng.get_handle(get_start_of(
                                          *currChain).node_id(), currID.second);
                            } else {
                                lastHandle = currHandle;
                            } 
                        } else {//Snarl
                            auto snarlDists = snarlIndex.find(make_pair(
                                        currSnarl->start().node_id(),
                                        currSnarl->start().backward()));
                            if (snarlDists != snarlIndex.end()) {//Already found
                                nodeLen = snarlDists->second.snarlLength(
                                                         graph, currSnarl);

                            } else {//Haven't recursed on snarl yet
                                #ifdef indexTraverse
                                    cerr << " recurse" << endl;
                                #endif
                                Chain currChain;
                                currChain.push_back(currSnarl);
                                nodeLen = calculateIndex(&currChain);
                            }
                           //Handle to last node in snarl
                            if (currID.first == currSnarl->start().node_id()
                             && currID.second == currSnarl->start().backward()){
                                lastHandle = ng.get_handle(
                                        currSnarl->end().node_id(), 
                                        currSnarl->end().backward());
                            } else if (currID.first ==currSnarl->end().node_id()
                                     &&
                                 currID.second != currSnarl->end().backward()) {

                                lastHandle = ng.get_handle(
                                          currSnarl->start().node_id(), 
                                          !currSnarl->start().backward());
                            } else {
                                lastHandle = currHandle;
                            }
                                                         
                        }
                    } else { //Node is just a node
                        Node* n = graph->get_node(currID.first);
                        nodeLen = n->sequence().size();
                        lastHandle = currHandle;
                    }
       
                    //Add reachable nodes to priority queue
                    vector<pair<id_t, bool>> nextNodes;
                    auto addHandle = [&](const handle_t& h)-> bool {
                         nextNodes.push_back(make_pair(
                                    ng.get_id(h), ng.get_is_reverse(h)));
                         #ifdef indexTraverse
                             cerr << ng.get_id(h) << " ";
                         #endif
                         return true;
                    };

                    #ifdef indexTraverse
                         cerr << "    Adding next nodes: ";
                    #endif
                    ng.follow_edges(lastHandle, false, addHandle);
                        
                    #ifdef indexTraverse
                         cerr << endl;
                    #endif
                    for (pair<id_t, bool> node : nextNodes) { 
                       reachable.push(make_pair(node, currDist + nodeLen)); 
                    }
                    NodeSide nodeSide = NodeSide(ng.get_id(lastHandle), 
                                               !ng.get_is_reverse(lastHandle));
                    if (graph->has_edge(nodeSide,nodeSide)) {
                        //If there is a loop - may not be in netgraph
                        reachable.push(make_pair(make_pair(currID.first,
                                          !currID.second), currDist + nodeLen));
                    }
                } 
            }//End while loop
        }//End for loop over starting node/directions in a snarl

        #ifdef indexTraverse
            cerr << "End snarl " << snarl->start().node_id() << endl;
        #endif

        /*Add to prefix sum the distance to the beginning and end of the last
            node in the current snarl
        */
        int64_t dist = sd.snarlDistance(
                make_pair (snarl->start().node_id(), snarl->start().backward()),
               make_pair (snarl->end().node_id(), snarl->end().backward()) );
        chainPrefixSum.push_back(chainPrefixSum[chainPrefixSum.size()-2]+ dist);
        chainPrefixSum.push_back(chainPrefixSum[chainPrefixSum.size()-1] + 
                    graph->get_node(snarl->end().node_id())->sequence().size());
       
        //Add reverse loop distances- from start node back to start node
        int64_t revLoopDist = sd.snarlDistance(
              make_pair (snarl->start().node_id(), !snarl->start().backward()),
             make_pair (snarl->start().node_id(), snarl->start().backward()));
        if (chainLoopRev.size() == 0) { //if first snarl in chain
            chainLoopRev.push_back(revLoopDist);
        } else { 
            int64_t lastLoop = chainLoopRev.back();
            if (lastLoop == -1 || revLoopDist != -1) {
                chainLoopRev.push_back(revLoopDist);
            } else {
                //Push 2*length of prev snarl + loop dist of prev snarl
                chainLoopRev.push_back(lastLoop +  
                     (chainPrefixSum.at(chainPrefixSum.size()-1) -
                             chainPrefixSum.at(chainPrefixSum.size()-3)) +
                     (chainPrefixSum.at(chainPrefixSum.size()-2) -
                           chainPrefixSum.at(chainPrefixSum.size()-4)) );
            }
        }

    }//End for loop over snarls in chain

    //Add forward loop distances 
   
    //Check if there is an edge traversing last node in chain fd -> rev 
    id_t lastNodeID = get_end_of(*chain).node_id();
    bool firstRev=graph->has_edge(node_end(lastNodeID), node_end(lastNodeID));

    if (firstRev) {
        int64_t lastDist = graph->get_node(chain->back()->end().node_id())->
                    sequence().size();
        chainLoopFd.push_back(lastDist);
    } else {
        chainLoopFd.push_back(-1);
    }
    for (int i = chain->size()-1; i >= 0 ; i--) { 
        const Snarl* snarl = chain->at(i);
        auto sd =  snarlIndex.at(make_pair(
            snarl->start().node_id(), snarl->start().backward())); 
        int64_t fdLoopDist = sd.snarlDistance(
              make_pair (snarl->start().node_id(), snarl->start().backward()),
              make_pair (snarl->start().node_id(), !snarl->start().backward()));
        int64_t lastLoop = chainLoopFd.back();
        if (lastLoop == -1 || fdLoopDist != -1) {
            chainLoopFd.push_back(fdLoopDist);
        } else {
        //push dist to end of snarl + loop dist + dist to start of snarl 

            chainLoopFd.push_back(lastLoop + 
                   (chainPrefixSum.at(2*i+2) - chainPrefixSum.at(2*i)) + 
                       (chainPrefixSum.at(2*i+3) - chainPrefixSum.at(2*i+1))); 
        }           
      
    }
    reverse(chainLoopFd.begin(), chainLoopFd.end()); 
     
    if (chainPrefixSum.size() > 4) { //If chain and not just one snarl
        //i->cd[get_start_of(*chain).node_id()] = ChainDistances(snarlToIndex, chainPrefixSum);
        chainIndex.insert(make_pair(get_start_of(*chain).node_id(), 
                 ChainDistances(snarlToIndex, chainPrefixSum, chainLoopFd,
                                                               chainLoopRev)));
    }
    #ifdef indexTraverse
        cerr << "Return" << chainPrefixSum.back() << endl;
    #endif
    return chainPrefixSum.back();//return length of entire chain
};


DistanceIndex::DistanceIndex::DistanceIndex(VG* vg, SnarlManager* snarlManager){
    /*Wrapper for creating the distance index given a VG and snarl manager
    */
    graph = vg;
    sm = snarlManager;
    #ifdef indexTraverse
        cerr << endl << "Start distance calculation"<< endl;
    #endif
    const vector<const Snarl*> topSnarls = sm->top_level_snarls();
    calculateIndex(&topSnarls);
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


pair<pair<int64_t, int64_t>, const Snarl*> DistanceIndex::distToCommonAncestor(
          const Snarl* snarl, const Snarl* commonAncestor, pos_t& pos){

    /* Find the distance from pos to both boundary nodes of commonAncestor.
       Return the two distances and the node_id and Snarl of the ancestor snarl
       whose parent is the commonAncestor - distances to ends of a node in 
       the common ancestor
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
//TODO: Off by 1 error? - includes both positions in distance

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

    pair<int64_t, int64_t> endDists = snarlDists.distToEnds(nodeID, distL, 
                                                                        distR);
    distL = endDists.first;
    distR = endDists.second;
    #ifdef printDistances
      cerr << nodeID << "->" << startID << ": " << distL << ", " << distR << endl;

    #endif
    nodeID = startID;
    bool nodeRev = startRev;

    while (sm->parent_of(snarl) != commonAncestor) {
        int64_t dsl; int64_t dsr; int64_t der; int64_t del;
    assert(distL >= -1 && distR >= -1);//don't accumulate -1's
        if (sm->in_nontrivial_chain(snarl)) {//Get to ends of chain
            const Chain* chain = sm->chain_of(snarl);
            id_t chainStartID = get_start_of(*chain).node_id();
            ChainDistances& chainDists =  chainIndex.at(chainStartID);

            pair<int64_t, int64_t> endDists = chainDists.distToEnds(
                          make_pair(nodeID, nodeRev), distL, distR);
            distL = endDists.first;   
            distR = endDists.second;

            nodeID = chainStartID; 
            nodeRev = get_start_of(*chain).backward();
    #ifdef printDistances
        cerr << nodeID << "->" << chainStartID << ": " << distL << ", " << distR
                                                                   << endl;
    #endif
        }
        snarl = sm->parent_of(snarl);
        id_t startNodeID = snarl->start().node_id();
        id_t startNodeRev = snarl->start().backward();
        SnarlDistances& snarlDists = 
                            snarlIndex.at(make_pair(startNodeID, startNodeRev));
            
       pair<int64_t, int64_t> endDists = snarlDists.distToEnds(
                                                      nodeID, distL, distR);
          
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


int64_t DistanceIndex::distance(const Snarl* snarl1, const Snarl* snarl2, 
                                                  pos_t& pos1, pos_t& pos2) {
    /*Actual distance function to get the distance between two positions
    */
    
/*TODO: Assume given snarls - find a way to get snarl from just position*/
    if (get_id(pos1) == get_id(pos2)) { //if positions are on the same node
        int64_t offset1;
        if (is_rev(pos1)) {
            offset1 = graph->get_node(get_id(pos1))->sequence().size() -
                      get_offset(pos1) - 1;//Len of node - offset 
        } else {
            offset1 = get_offset(pos1);
        }
        int64_t offset2;
        if (is_rev(pos2)) {
            offset2 = graph->get_node(get_id(pos2))->sequence().size() - 
                      get_offset(pos2) - 1;
        } else {
            offset2 = get_offset(pos2);
        }
        return abs(offset1-offset2)+1; //+1 to be consistent


    } else { //Positions are on different nodes

        id_t nodeID1 = get_id(pos1);
        bool nodeRev1 = is_rev(pos1);
        id_t nodeID2 = get_id(pos2); 
        bool nodeRev2 = is_rev(pos2);
    #ifdef printDistances
        cerr << endl << "Start distance calculation from " << nodeID1 << "->" <<
             nodeID2 << endl;
    #endif
        //const Snarl* snarl1 = sm->into_which_snarl(node1, is_rev(pos1));
        //const Snarl* snarl2 = sm->into_which_snarl(node2, is_rev(pos2));
        const Snarl* commonAncestor = NULL; 

    #ifdef printDistances
        cerr << "Find common ancestor";
    #endif
        //// Find common ancestor
        unordered_set<const Snarl*> ancestors1;
                      //set of all ancestor snarls of node1
        const Snarl* ancestor1 = snarl1;

        while (ancestor1 != NULL) {
            ancestors1.emplace(ancestor1);
            ancestor1 = sm->parent_of(ancestor1);
        }

    #ifdef printDistances
        cerr << "Ancestors of 1: ";
        for (auto ancestorSnarl : ancestors1) {

            cerr << ancestorSnarl->start().node_id() << " ";
        }
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
        
        //Snarl1 and snarl2 are children of common ancestor or common ancestor

    #ifdef printDistances
        cerr << "Distances to snarl in common ancestor: " << distL1 << ", " <<
             distR1 << "   " << distL2 << ", " << distR2 << endl;
    #endif
        int64_t chainDist = -1; 
        //Find shortest distance between boundary nodes of snarls containing pos
        if (snarl1 != commonAncestor && snarl2 != commonAncestor && 
              sm->in_nontrivial_chain(snarl1) && sm->in_nontrivial_chain(snarl2)
               && sm->chain_of(snarl1) == sm->chain_of(snarl2)) {
            //If positions are in the same chain within common ancestor
            const Chain* chain = sm->chain_of(snarl1);
            id_t chainStartID = get_start_of(*chain).node_id();
            bool chainStartRev = get_start_of(*chain).backward();

            ChainDistances& chainDists = chainIndex.at( chainStartID); 

            //Distance from left of s1 (reverse), left of s2 (forward)
            int64_t d1 = chainDists.chainDistanceShort(graph,
                   make_pair(nodeID1, nodeRev1==chainStartRev), 
                   make_pair(nodeID2, nodeRev2!=chainStartRev));
            d1 = (distL1 == -1 || distL2 == -1 || d1 == -1) ? -1 : 
                                    distL1 + distL2 + d1; 
            //Distance from left of s1 (reverse) to right of s2 (reverse)
            int64_t d2 = chainDists.chainDistanceShort(graph,
                   make_pair(nodeID1, nodeRev1==chainStartRev), 
                   make_pair(snarl2->end().node_id(), nodeRev2==chainStartRev));
            if (nodeID1 == snarl2->end().node_id()) {
                d2 = distL1 + distR2 - d2; 
            } else {
                d2 = (distL1 == -1 || distR2 == -1 || d2 == -1) ? -1 : 
                                       distL1 + distR2 + d2; 
            }
            //Distance from right of s1 (fd) to left of s2 (fd)
            int64_t d3 = chainDists.chainDistanceShort(graph,
                   make_pair(snarl1->end().node_id(),nodeRev1!=chainStartRev), 
                           make_pair(nodeID2, nodeRev2!=chainStartRev));
            if(snarl1->end().node_id() == nodeID2) {
                d3 = (distR1 == -1 || distL2 == -1 || d3 == -1) ? -1 : 
                                       distR1 + distL2 - d3; 
            } else {
                d3 = (distR1 == -1 || distL2 == -1 || d3 == -1) ? -1 : 
                                       distR1 + distL2 + d3; 
            }
            //Distance from right of s1 (fd) to right of reverse (fd)
            int64_t d4 =  chainDists.chainDistanceShort(graph,
                  make_pair(snarl1->end().node_id(), nodeRev1!=chainStartRev),
                  make_pair(snarl2->end().node_id(), nodeRev2==chainStartRev));
            d4 = (distR1 == -1 || distR2 == -1 || d4 == -1) ? -1 : 
                                       distR1 + distR2 + d4;
            
            vector<int64_t> valsv {d1, d2, d3, d4};
                       
            chainDist = minPos(valsv);
    #ifdef printDistances
            cerr << "    Possible distances within chain: " << d1 << " " << d2
                 << " " << d3 << " " << d4 << endl;
            cerr << "Chain distance in common ancestor: " << chainDist << endl;
            
    #endif

        }
        if (commonAncestor == NULL) {
            return chainDist;
        } else { 
            //Get dist from pos1 to ends of chain 
            if (snarl1 != commonAncestor && sm->in_nontrivial_chain(snarl1)) {
                const Chain* chain = sm->chain_of(snarl1);

                id_t chainStartID = get_start_of(*chain).node_id();
                ChainDistances& chainDists = chainIndex.at( chainStartID);
                pair<int64_t, int64_t> endDists = chainDists.distToEnds(
                                make_pair(nodeID1, nodeRev1), distL1, distR1);

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
                           make_pair(nodeID2, nodeRev2), distL2, distR2);    
                distL2 = endDists.first;
                distR2 = endDists.second;

                nodeID2 = chainStartID; 
                nodeRev2 = get_start_of(*chain).backward();
            }
            
    
 
        #ifdef printDistances
            cerr << "Distances to common ancestor: " << distL1 << ", " << distR1
                  << "   " << distL2 << ", " << distR2 << endl;
        #endif
            //Both nodes are nodes in common ancestor
            //Get distance between end nodes in snarl
            NetGraph ng = NetGraph(commonAncestor->start(), 
                    commonAncestor->end(),sm->chains_of(commonAncestor), graph);
            SnarlDistances& snarlDists = snarlIndex.at(make_pair(
                               commonAncestor->start().node_id(),
                                           commonAncestor->start().backward()));
            int64_t d1 = snarlDists.snarlDistanceShort(&ng, sm,
                    make_pair(nodeID1, nodeRev1), make_pair(nodeID2, nodeRev2));
            d1 = (distR1 == -1 || distL2 == -1 || d1 == -1) ? -1 : 
                                                   distR1 + distL2 + d1; 
            int64_t d2 = snarlDists.snarlDistanceShort(&ng,sm,
                   make_pair(nodeID1, nodeRev1), make_pair(nodeID2, !nodeRev2));

            d2 = (distR1 == -1 || distR2 == -1 || d2 == -1) ? -1 :
                                                   distR1 + distR2 + d2;
            int64_t d3 = snarlDists.snarlDistanceShort(&ng,sm,
                  make_pair(nodeID1, !nodeRev1), make_pair(nodeID2, nodeRev2));
            d3 = (distL1 == -1 || distL2 == -1 || d3 == -1) ? -1 : 
                                                  distL1 + distL2 + d3; 
            int64_t d4 = snarlDists.snarlDistanceShort(&ng,sm,
                  make_pair(nodeID1, !nodeRev1), make_pair(nodeID2, !nodeRev2));
            d4 = (distL1 == -1 || distR2 == -1 || d4 == -1) ? -1 : 
                                                  distL1 + distR2 + d4; 
            vector<int64_t> vals {d1, d2, d3, d4, chainDist}; 

        #ifdef printDistances
            cerr << "Distances within common ancestor: " << d1 << ", " << d2
                                            << ", " << d3 << ", " << d4 << endl;
        #endif
            return minPos(vals);
        }
    }
};


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
