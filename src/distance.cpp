#include "distance.hpp"
using namespace std;
namespace vg {

    //////////  Helper Functions

    size_t snarlMatrixIndex (visit_to_index& indices, int start, int end) {
        //Get the index of dist from start to end in a snarl distance matrix
        size_t i1 = indices[start];
        size_t i2 = indices[end];
        size_t length = visit_to_index.size();
        return (length + 1) * (length / 2) + (i1+1)* (i1/2) + i2 - 1;
    };
   
    int64_t snarlDistance(SnarlDistance& sd,id_t& snarl, int n1, int n2){
        //Distance between two nodes n1 and n2 in snarl
        //nodes are node id- negative if reverse, positive if forward
        //TODO: the distance index is passed as an argument for now
        size_t i = snarlMatrixIndex(sd.visit_to_index, n1, n2);
        return sd.distances[i]; 
    };
   
    int64_t chainDistance(ChainDistances& cd, int v1, int v2) {
        //Distance between v1 and v2 in a chain
        size_t i1 = cd.visit_to_index(v1);
        size_t i2 = cd.visit_to_index(v2); 
        return abs(cd.prefix_sum[i1] - cd.prefix_sum[i2]);
    };

    int64_t chainLength(DistanceIndex& di, Chain& chain) {
        //Get the length of a chain including length of last node
        return di.cd.prefix_sum.back();
    };





    /////////   Create the distance index

    DistanceIndex makeDistanceIndex(VG graph, SnarlManager sm) {
        //Create the distance index given a VG and snarl manager
        DistanceIndex distances;
        return distances;
    };

    int64_t calculateIndex(DistanceIndex& di, VG& graph, SnarlManager& sm, 
                 Chain& chain) {
        bool cmp = [](pair<Node&, int64_t> x, pair<Node&, int64_t> y){
                 return x.second < y.second; };
        //Comparison of two pairs of <Node, distance> for priority queue
 
        vector<int64_t> chain_prefix_sum (1,0); //initialize to [0]
        unordered_map<int, size_t> chain_snarl_to_index;
        for (Snarl snarl = chain.begin(); snarl != chain.end(); ++snarl) { 
            //for each snarl in the chain TODO ChainIterator?
            chain_snarl_to_index.insert( make_pair<int, size_t>
                          (snarl.start().node_id(), chain_prefix_sum.size()-1));
            //Initialize components of SnarlDistances struct
            //TODO use visit as key to identify a node+direction? or node id? 
            unordered_map<int, size_t> snarl_visit_to_index;
            pair<unordered_set<Node*>, unordered_set<Edge*>> contents = 
                            sm.shallow_contents(snarl, graph, true);
            unordered_set<Node*> allNodes = contents.first;
            //all nodes in snarl includes both boundary nodes of child snarls
            /*TODO better way of getting all the nodes in the snarl? Currently
              don't include end boundary node of a child snarl in distance 
              matrix, maybe should?
            */
            size_t snarlIndex = 0;
            for (Node* node = allNodes.begin(); 
                    node != allNodes.end(); ++node) {
                snarl_visit_to_index[node->id] -> snarlIndex++;
                snarl_visit_to_index[-node->id] -> snarlIndex++;
            }
            int size = snarl_visit_to_index.size();
            vector<int64_t> snarl_distances ((size+1) * size / 2,-1);
            SnarlDistances sd;
            sd.visit_to_index = &snarl_visit_to_index;
            sd.distances = &snarl_distances;
            di.sd.insert(make_pair<id_t, snarlDistances*>
                                             (snarl.start().node_id(),&sd));
            
            //Create a NetGraph for current snarl
            NetGraph ng = NetGraph(snarl.start(), snarl.end(), chains??, graph);
            //TODO make a netgraph - need child chains??
            
            // For each node in snarl calculate distance to every reachable node
            for ( Node* startNode = allNodes.begin(); 
                    startNode != allNodes.end(); ++startNode){
                bool bools [2] = {true, false};
                for (bool rev : bools) {
                //TODO Better way to loop over all nodes/direction??
                if (rev) { int startID = -startNode->node_id(); }//reverse
                else     { int startID = startNode->node_id(); }
                handle_t handle = ng.get_handle(startNode->node_id,
                                                                is_reverse=rev);
                priority_queue<pair<handle_t&, int64_t>> reachable (cmp);
                pair<handle_t&, int64_t> init (handle, 0);
                reachable.push(init);
                while (reachable.size() > 0) {
                    pair<handle_t&, int64_t> next = reachable.pop();
                    handle_t* currHandle = next.first;
                    int64_t currDist = next.second;
                    if ng.get_is_reverse(currHandle) {
                        int currID = -currHandle.get_id();
                    } else {
                        int currID = currHandle.get_id();
                    }
 
                    if ( snarlDistance(sd, startID, startID, currID) == -1) {
                        //If node has not already been found:
   
                        //Get the length of the current node
                        index =  snarlMatrixIndex (&snarl_visit_to_index, 
                                                         startID, currID)
                        snarl_distances[index] = currDist;
                        if (abs(currID) != snarl.start().node_id() &&
                                                       ng.is_child(currHandle)){
                            //If a child snarl/chain begins at the current node
                            Snarl* currSnarl = sm.into_which_snarl(abs(currID), 
                                               (currID < 0));
                            if (sm.in_nontrivial_chain(currSnarl)) {//Chain
                                /*TODO assuming start and end are consistent 
                                  within the graph, not relative to the current
                                  orientation
                                */
                                Chain currChain = sm.chain_of(currSnarl);
                                unordered_map<id_t, chainDistances*>::iterator 
                                     chainDists = di.cd.find(
                                       get_start_of(chain).node_id());
                                if (chainDists != di.cd.end()) {
                                    //Length of chain has already been found
                                    size_t nodeLen = 
                                                  chainDists->prefix_sum.back();
                                    //last element should be length of chain
                                } else {//haven't recursed on this chain yet
                                    size_t nodeLen = calculateIndex( 
                                                     di, graph, sm, currChain);
                                }
                            } else {//Snarl
                                unordered_map<id_t, snarlDistances*>::iterator 
                                     snarlDists = di.sd.find(currID);
                                if (snarlDists != di.sd.end()) {//Already found
                                    size_t nodeLen = snarlDistance(snarlDists,
                                         currID, currSnarl.start().node_id(),
                                         currSnarl.end().node_id());
                                } else {//Haven't recursed on snarl yet
                                    Chain currChain;
                                    currChain.insert(&snarl);
                                    size_t nodeLen = calculateIndex( 
                                                     di, graph, sm, currChain);
                                }
                            }
                        } else { //Node is just a node
                            size_t nodeLen = ng.get_length(currHandle);
                            //TODO This doesn't work?
                        }

                        //Add next nodes to priority queue
                        vector<handle_t*> nextNodes;//vector of adjacent nodes
                        auto addHandle = [&](const handle_t& h)-> bool {
                            nextNodes.insert(&h);
                            return true;
                        };
                        ng.follow_edges(currHandle, false, addHandle);
                        for (handle_t* h = nextNodes.begin(); 
                                               h != nextNodes.end(); ++h) {
                            
                            pair<handle_t&, int64_t> p (h, currDist + nodeLen);
                            reachable.push(p);
                        }
                    } 
                }//End while loop
            }}//End for loop over starting node/directions
            /* add length of snarl (start of start node to start of end node) to
               the chain prefix sum vector */
            int64_t dist = snarlDistance(di, &snarl.start().node_id(), 
                                snarl.start().node_id(), snarl.end().node_id());
            chain_prefix_sum.push_back(chain_prefix_sum.back() + dist);
        }//End for loop over snarls in chain
        //Add the length of the last node in chain to get length of entire chain
        Visit lastVisit = get_end_of(chain).snarl.end;
        //TODO Visit points <-, snarl as well???
        Node lastNode = graph.node_by_id(lastVisit.node_id);
        
        chain_prefix_sum.push_back(chain_prefix_sum.back() + 
               lastNode.sequence.size());
 
        if (chain_prefix_sum.size() > 2) { //If chain and not just one snarl
            ChainDistances cd;
            cd.snarl_to_index = &chain_snarl_to_index;
            cd.prefix_sum = &chain_prefix_sum;
            di.cd.insert(make_pair<id_t, chainDistances*>
                                 (snarl.start().node_id(), &cd));
        }
        return chain_prefix_sum.back();//return length of entire chain
    };

    ////////    Calculate distances
}

