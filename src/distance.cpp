//#define indexTraverse
//#define debug
//#define printDistances

#include "distance.hpp"

using namespace std;
namespace vg {

int64_t minPos (vector<int64_t> vals) {
    /*return the minimum value in vals that is not -1, returns -1 if all
     values are -1 */
    return accumulate(vals.begin(), vals.end(), -1, 
          [](int x, int y) {if (x==-1) {return y;} 
                            else if (y == -1) {return x;}
                            else {return min(x, y);}} 
          ); 
   
};



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
    /*Distance between end of node n1 and beginning of node n2 in snarl
    */
    #ifdef debug
    cerr << start.first << start.second << "->" << end.first << end.second 
                  << endl;
    printSelf();
    #endif

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

//TODO: Can probably just get the length from the index
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
    Node* firstNode;
    firstNode = graph->get_node(chain->at(0)->start().node_id());
    chainPrefixSum.push_back(firstNode->sequence().size());
    unordered_map<id_t, size_t> snarlToIndex;
    snarlToIndex[get_start_of(*chain).node_id()] = snarlToIndex.size()-1;

    for (int i = 0; i < chain->size(); i++) {
        //for each snarl in the chain TODO ChainIterator?
        const Snarl* snarl = chain->at(i);
        if (snarlToIndex.find(snarl->end().node_id()) == snarlToIndex.end()){
            snarlToIndex[snarl->end().node_id()] = snarlToIndex.size()-1;
        }


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
            bool firstLoop = true;

            while (reachable.size() > 0) {
                pair<pair<id_t, bool>, int64_t> next = reachable.top();
                reachable.pop();
                pair<id_t, bool> currID= next.first;
                int64_t currDist = next.second;
                if ( sd.snarlDistance(startID, currID) == -1) {
                    //If node has not already been found:

                    //Save distance from start to current node 
                    if (currDist == 0) {
                        sd.insertDistance(startID, currID, -1);
                    } else {
                        sd.insertDistance(startID, currID, currDist);
                    }
               

                    const handle_t currHandle = ng.get_handle(currID.first, 
                                                               currID.second);

                    
                    int64_t nodeLen; //length of the current node
                       
                    int64_t loopDist = -1;
                         //Dist to enter curr node then exit at same side 

                    const Snarl* tempSnarl = sm->into_which_snarl(
                                              currID.first, currID.second);
                    const Snarl* currSnarl = tempSnarl == NULL ? 
                        sm->into_which_snarl(currID.first, !currID.second) :
                        tempSnarl; 

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
                                if (get_start_of(
                                     *currChain).backward() == currID.second) {
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
                                //last element should be length of chain

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

                                if (currID.second == snarlRev) { 
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
                                Chain currChain;
                                currChain.push_back(currSnarl);
                                calculateIndex(&currChain);

                                SnarlDistances& currSnarlDists = 
                                     snarlIndex.at(make_pair(snarlID,snarlRev));
                                nodeLen = currSnarlDists.snarlLength(); 
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
       
                    if (loopDist != -1 && !firstLoop) {
//TODO: I think this will add paths that exist but not in netgraph without having to use internal connectivity in netgraph
                    //If there is a reversing path from within the current node
                    //Add nodes reachable by reversing - leaving snarl at node entered


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

                    #ifdef indexTraverse
                         cerr << "    From start node " << startID.first << " " << startID.second<< " in snarl " << snarl->start().node_id() << " at " << ng.get_id(currHandle) << " " << ng.get_is_reverse(currHandle) << endl; 
                         cerr << "        Adding next nodes:  ";
                    #endif
/* TODO: Hopefully won't need this once path past common ancestor is found
         but might for a reversing edge at the end of a snarl???

                    NodeSide nodeSide = NodeSide(ng.get_id(currHandle), 
                                               !ng.get_is_reverse(currHandle));
                    if (graph->has_edge(nodeSide,nodeSide)) {
                        //If there is a loop - may not be in netgraph
                        reachable.push(make_pair(make_pair(currID.first,
                                          !currID.second), currDist + nodeLen));
                     
                    }
*/  
                    if (nodeLen != -1) {

                        ng.follow_edges(currHandle, false, addHandle);
                    }
                        

                    //Add edges between boundaries that leave snarl (not in ng)

                    id_t snarlStartID = snarl->start().node_id();
                    id_t snarlStartRev = snarl->start().backward();
                    id_t snarlEndID = snarl->end().node_id();
                    id_t snarlEndRev = snarl->end().backward();   //pointing out

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
        int64_t dist = sd.snarlDistance(
           make_pair (snarl->start().node_id(), snarl->start().backward()),
           make_pair (snarl->end().node_id(), snarl->end().backward()) );

        chainPrefixSum.push_back(chainPrefixSum[chainPrefixSum.size()-2]+
                       dist);
        chainPrefixSum.push_back(chainPrefixSum[chainPrefixSum.size()-1] + 
                    graph->get_node(snarl->end().node_id())->sequence().size());
        
         //length of snarl
        if (dist == -1) {
            sd.length = -1;
        } else {
            sd.length = dist  
               + 
               graph->get_node(snarl->end().node_id())->sequence().size();
        }

        //Add reverse loop distances- from start node rev to start node fd
//TODO: I think this assumes that the chain traverses all snarls forwards
        int64_t revLoopDist = sd.snarlDistance(
              make_pair (snarl->end().node_id(), !snarl->end().backward()),
             make_pair (snarl->end().node_id(), snarl->end().backward()));

        if (chainLoopRev.size() == 0) { 

            //if first snarl in chain, push dist for first node before
            // last node in snarl
            int64_t revDist = sd.snarlDistance(
              make_pair (snarl->start().node_id(), !snarl->start().backward()),
               make_pair (snarl->start().node_id(), snarl->start().backward()));

            chainLoopRev.push_back(revDist);

        } 

        int64_t lastLoop = chainLoopRev.back();

        if (lastLoop == -1) {

            chainLoopRev.push_back(revLoopDist);

        } else {

            //Push 2*length of prev snarl + loop dist of prev snarl
            int64_t loopDistance = minPos({revLoopDist, lastLoop +  
                 sd.snarlDistance(
                make_pair (snarl->end().node_id(), !snarl->end().backward()),
               make_pair (snarl->start().node_id(), !snarl->start().backward()))
            + 
                 sd.snarlDistance(
                make_pair (snarl->start().node_id(), snarl->start().backward()),
                 make_pair (snarl->end().node_id(), snarl->end().backward()))});

            chainLoopRev.push_back(loopDistance);
/* TODO: Fix chain loop distance
                 (chainPrefixSum.at(chainPrefixSum.size()-1) -
                         chainPrefixSum.at(chainPrefixSum.size()-3)) +
                 (chainPrefixSum.at(chainPrefixSum.size()-2) -
                       chainPrefixSum.at(chainPrefixSum.size()-4)) );
        
*/
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
        if (lastLoop == -1) {
            chainLoopFd.push_back(fdLoopDist);
        } else {
        //push dist to end of snarl + loop dist + dist to start of snarl 

            int64_t loopDistance = minPos({fdLoopDist, lastLoop + 
                   (chainPrefixSum.at(2*i+2) - chainPrefixSum.at(2*i)) + 
                       (chainPrefixSum.at(2*i+3) - chainPrefixSum.at(2*i+1))});
            chainLoopFd.push_back(loopDistance); 
        }           
      
    }
    reverse(chainLoopFd.begin(), chainLoopFd.end()); 
     
    if (chainPrefixSum.size() > 4) { //If chain and not just one snarl
        chainIndex.insert(make_pair(get_start_of(*chain).node_id(), 
                 ChainDistances(snarlToIndex, chainPrefixSum, chainLoopFd,
                                                               chainLoopRev)));
    }
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

    
    pair<int64_t, int64_t> endDists = snarlDists.distToEnds(nodeID, false, distL, 
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
        if (sm->in_nontrivial_chain(snarl)) {//Get to ends of chain
            const Chain* chain = sm->chain_of(snarl);
            id_t chainStartID = get_start_of(*chain).node_id();
            ChainDistances& chainDists =  chainIndex.at(chainStartID);

            pair<int64_t, int64_t> endDists = chainDists.distToEnds(
                          make_pair(nodeID, nodeRev!=get_start_of(*chain).backward()), distL, distR);
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


int64_t DistanceIndex::distance(const Snarl* snarl1, const Snarl* snarl2, 
                                                  pos_t& pos1, pos_t& pos2) {
    /*Actual distance function to get the distance between two positions
    */
    
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

    //// Find common ancestor
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
            //If snarls share a node, chainDistanceShort returns length of 
            //shared node
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

        //Distance from right of s1 (fd) to right of s2 (rev)
        int64_t d4 =  chainDists.chainDistanceShort(graph,
              make_pair(snarl1->end().node_id(), nodeRev1!=chainStartRev),
              make_pair(snarl2->end().node_id(), nodeRev2==chainStartRev));
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
    //Get dist from pos1 to ends of chain 
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
    
    
    //Find distances to the ends of the common ancestor
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

    //Paths between common ancestor's start and end nodes
    const Snarl* currSnarl = commonAncestor;
    const Snarl* parentSnarl = sm->parent_of(currSnarl);
    id_t startID = currSnarl->start().node_id();
    id_t startRev = currSnarl->start().backward(); //pointing into snarl
    id_t endID = currSnarl->end().node_id();
    id_t endRev = currSnarl->end().backward();     //pointing out

    NodeSide startSide = NodeSide(startID, startRev); //start of start
    NodeSide endSide = NodeSide(endID, !endRev);       //end of end
/* TODO: I think this is redundant
    //Dist from start to start
    snarlDists = snarlIndex.at(make_pair(
                           currSnarl->start().node_id(),
                                       currSnarl->start().backward()));


    d1 = snarlDists.snarlDistanceShort(graph, &ng, 
                make_pair(startID, !startRev), make_pair(startID, startRev));
    d1 = (distL1 == -1 || distL2 == -1 || d1 == -1) ? -1 : 
                                                  distL1 + distL2 + d1; 
    
    //Dist from end to end
    d2 = snarlDists.snarlDistanceShort(graph, &ng,
                   make_pair(endID, endRev), make_pair(endID, !endRev));

    d2 = (distR1 == -1 || distR2 == -1 || d2 == -1) ? -1 :
                                                   distR1 + distR2 + d2;
    //Dist from start to end
    d3 = snarlDists.snarlDistanceShort(graph, &ng,
                make_pair(endID, endRev), make_pair(startID, startRev));
    d3 = (distL1 == -1 || distR2 == -1 || d3 == -1) ? -1 : 
                                              distL1 + distR2 + d3; 
    d4 = (distR1 == -1 || distL2 == -1 || d4 == -1) ? -1 : 
                                                  distR1 + distL2 + d4; 


    shortestDistance =  minPos({d1, d2, d3, d4, chainDist, shortestDistance});

#ifdef printDistances
    cerr << "Distances within common ancestor: " << d1 << ", " << d2
                                        << ", " << d3 << ", " << d4 << endl;
    cerr << "Shortest dist only up to  common ancestor: " << shortestDistance
        << endl;
#endif
*/

    /*shortestDistance is now the shortest distance only traversing up to the 
      most recent common ancestor.   
      
      currSnarl is the common ancestor, start/end ID are a node in the 
      common ancestor, distances are up to a node in the common ancestor
      Traverse up to root and check for path at each level
    */
//TODO: Might be faster to store edges that loop in snarls/chains in the index rather than using has_edge
    
    while ( currSnarl != NULL) {

            
        if (sm->in_nontrivial_chain(currSnarl)) {
            //Find paths between ends of current chain 

            const Chain* currChain= sm->chain_of(currSnarl);
            ChainDistances& chainDists = chainIndex.at(
                                            get_start_of(*currChain).node_id());
            //TODO: If too slow, this is done twice        
            bool chainStartRev = get_start_of(*currChain).backward(); 

            //Distance from start (reverse) to start (forward)
            int64_t d1 = chainDists.chainDistanceShort(graph,
               make_pair(startID, startRev==chainStartRev), 
                   make_pair(startID, startRev!=chainStartRev));
            d1 = (distL1 == -1 || distL2 == -1 || d1 == -1) ? -1 : 
                                distL1 + distL2 + d1; 

            //Distance from start (reverse) to end (reverse)
            int64_t d = chainDists.chainDistanceShort(graph,
               make_pair(startID, startRev==chainStartRev), 
               make_pair(endID, endRev==chainStartRev));

            d2 = (distL1 == -1 || distR2 == -1 || d == -1) ? -1 : 
                                   distL1 + distR2 + d; 
            d3 = (distR1 == -1 || distL2 == -1 || d == -1) ? -1 : 
                                   distR1 + distL2 + d; 

           //Distance from end (fd) to end (rev)
            int64_t d4 =  chainDists.chainDistanceShort(graph,
                  make_pair(endID, endRev!=chainStartRev),
                  make_pair(endID, endRev==chainStartRev));
            d4 = (distR1 == -1 || distR2 == -1 || d4 == -1) ? -1 : 
                                   distR1 + distR2 + d4;
        
                   
            shortestDistance = minPos({shortestDistance, d1, d2, d3, d4});

            pair<int64_t, int64_t> chainEnd1 = chainDists.distToEnds(
                          make_pair(startID, startRev!=get_start_of(*currChain).backward()), distL1, distR1);
            distL1 = chainEnd1.first; distR1 = chainEnd1.second;

            pair<int64_t, int64_t> chainEnd2 = chainDists.distToEnds(
                          make_pair(startID, startRev!=get_start_of(*currChain).backward()), distL2, distR2);
            distL2 = chainEnd2.first; distR2 = chainEnd2.second;

            startID = get_start_of(*currChain).node_id();
            startRev = get_start_of(*currChain).backward();
            endID = get_end_of(*currChain).node_id();
            endRev = get_end_of(*currChain).backward();

            //If an edge connects ends of two paths, find shortest distance
            NodeSide startSide = NodeSide(startID, startRev); //start of start
            NodeSide endSide = NodeSide(endID, !endRev);       //end of end


            if (graph->has_edge(startSide, endSide)) {
                shortestDistance = minPos({shortestDistance, distL1+distR2,
                                        distL2+distR1});
           }
#ifdef printDistances
    cerr << "At chain " << startID << " dists to ends: " << distL1 << " " << 
         distR1 << " " << distL2 << " " << distR2 << " Shortest distance : "
          << shortestDistance << endl;
#endif  
        }
   
        if (parentSnarl == NULL) {break;} //TODO: Make this cleaner

        SnarlDistances& snarlDists =snarlIndex.at(make_pair(parentSnarl->start().node_id(), parentSnarl->start().backward()));


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
        pair<int64_t, int64_t> endDists1 = snarlDists.distToEnds(startID, startRev, distL1,
                                                                       distR1);
        distL1= endDists1.first; distR1= endDists1.second;

        pair<int64_t, int64_t> endDists2 = snarlDists.distToEnds(startID, startRev, distL2,
                                                                       distR2);
        distL2= endDists2.first; distR2= endDists2.second;

        startID = parentSnarl->start().node_id();
        startRev = parentSnarl->start().backward();
        endID = parentSnarl->end().node_id();
        endRev = parentSnarl->end().backward();

        startSide = NodeSide(startID, startRev); //start of start
        endSide = NodeSide(endID, !endRev);       //end of end

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

void DistanceIndex::printSelf() {
    for (auto snarls : snarlIndex) {
        snarls.second.printSelf();
    }
    for (auto chains : chainIndex) {
        chains.second.printSelf();
    }
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
