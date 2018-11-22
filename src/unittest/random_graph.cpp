#include "random_graph.hpp"

namespace vg {
namespace unittest {

VG randomGraph(int64_t seqSize, int64_t variantLen,
                    int64_t variantCount){
    //Create a random graph for a sequence of length seqSize
    //variantLen is the mean length of a larger variation and variationCount
    //is the number of variations in the graph

    VG graph;
    map<size_t, id_t> indexToNode;
                  //Index of original sequence to node that starts at that index

    random_device seed_source;
    default_random_engine generator(seed_source());
    
    uniform_int_distribution<int> baseDistribution(0, 3);
    poisson_distribution<size_t> lengthDistribution(variantLen);
    uniform_int_distribution<int> variantDistribution(0, 4);
    uniform_int_distribution<int> indexDistribution(0, seqSize-1);
    
    // Init the string with placeholder bases
    string seq(seqSize, 'A');
    // Set the string to random bases
    string alphabet = "ACGT";
    for (size_t i = 0; i < seq.size(); i++) {
        seq[i] = alphabet[baseDistribution(generator)];
    }
    
    Node* n = graph.create_node(seq);
    indexToNode.insert(make_pair(0, n->id()));

    //Get random number generator for lengths and variation types

    auto splitNode = [&] (size_t index) -> pair<Node, Node> {
        //Split graph at index. Split the node with the original sequence into
        //the original node and a new node and connect the two

        auto n = --indexToNode.upper_bound(index);//orig node containing pos
        size_t firstIndex = n->first;           //Index of first node
        size_t firstNodeID = n->second;         //Node ID of first node

        Node* firstNode = graph.get_node(firstNodeID);
        string origSeq = firstNode->sequence();
        if (index > firstIndex) {
            size_t nodeLength = index - firstIndex;

            firstNode->set_sequence(origSeq.substr(0, nodeLength));//replace seq
            Node* newNode = graph.create_node(origSeq.substr(nodeLength));
            indexToNode.insert(make_pair(index, newNode->id()));


            //Transfer outgoing edges from fist node to second node

            handle_t startFd = graph.get_handle(firstNodeID, false);
            handle_t endFd = graph.get_handle(newNode->id(), false);


            unordered_set<handle_t> nextHandles;
            auto addHandle = [&](const handle_t& h) ->bool {
                nextHandles.insert(h);
                return true;
            };
            graph.follow_edges(startFd, false, addHandle);

            for (handle_t h : nextHandles) {
                //for each edge from start node, delete and add to new node
                graph.destroy_edge(startFd, h);
                graph.create_edge(endFd, h);
            }
            graph.create_edge(firstNode, newNode);
            return make_pair(*firstNode, *newNode);
        } else {

            auto n = --indexToNode.lower_bound(index);
            Node* prevNode = graph.get_node(n->second);
            return make_pair(*prevNode, *firstNode);
        }

    };

    enum VariationType {SNP = 0, POINT_INDEL = 1, STRUCTURAL_INDEL = 2,
                        CNV = 3, INVERSION = 4};

    for (int j = 0; j < variantCount; j++) {
        //add variants
        int startIndex = indexDistribution(generator);
        VariationType variationType = 
                                 (VariationType) variantDistribution(generator);

        if (variationType == SNP) {

            if (startIndex == 0) {

                pair<Node, Node> endNodes = splitNode(startIndex+1);
                Node end = endNodes.second;
                Node* newNode = graph.create_node(string(1, alphabet[baseDistribution(generator)]));
                graph.create_edge(newNode, &end);

            } else if (startIndex < seqSize-2) {

                pair<Node, Node> startNodes = splitNode(startIndex);
                pair<Node, Node> endNodes = splitNode(startIndex+1);

                Node start = startNodes.first;
                Node end = endNodes.second;
                Node* newNode = graph.create_node(string(1, alphabet[baseDistribution(generator)]));
                graph.create_edge(&start, newNode);
                graph.create_edge(newNode, &end);

            } else if (startIndex == seqSize -2) {

                pair<Node, Node> startNodes = splitNode(startIndex);

                Node start = startNodes.first;
                Node* newNode = graph.create_node(string(1, alphabet[baseDistribution(generator)]));
                graph.create_edge(&start, newNode);

            }
        } else if (variationType == POINT_INDEL) {
            //Short indel - deletion of original

            if (startIndex > 0 && startIndex < seqSize-1) {
                pair<Node, Node> startNodes = splitNode(startIndex);
                pair<Node, Node> endNodes = splitNode(startIndex+1);

                Node start = startNodes.first;
                Node end = endNodes.second;
                graph.create_edge(&start, &end);

            }
        } else if (variationType == STRUCTURAL_INDEL) {
            //long indel
            size_t length = lengthDistribution(generator);

            if (length > 0 && startIndex > 0 &&
               length + startIndex < seqSize-1) {

                pair<Node, Node> startNodes = splitNode(startIndex);
                pair<Node, Node> endNodes = splitNode(startIndex+length);

                Node start = startNodes.first;
                Node end = endNodes.second;
                graph.create_edge(&start, &end);

            }
        } else if (variationType == CNV) {
            //Copy number variation
            size_t length = lengthDistribution(generator);

            if (length > 0 ) {
                if (startIndex == 0) {
                    pair<Node, Node> endNodes = splitNode(startIndex+length);
                    Node n = endNodes.first;

                    auto nodePair = --indexToNode.upper_bound(0);//first node
                    size_t firstNodeID = nodePair->second;

                    Node* firstNode = graph.get_node(firstNodeID);

                    graph.create_edge(&n, firstNode);

                } else if ( length + startIndex < seqSize-1) {

                    pair<Node, Node> startNodes = splitNode(startIndex);
                    pair<Node, Node> endNodes = splitNode(startIndex+length);

                    Node n1 = startNodes.second;
                    Node n2 = endNodes.first;
                    graph.create_edge(&n2, &n1);

                } else  if ( length + startIndex == seqSize -1) {
                    pair<Node, Node> startNodes = splitNode(startIndex);

                    auto nodePair = --indexToNode.upper_bound(seqSize);
                        //last node
                    size_t lastNodeID = nodePair->second;

                    Node* lastNode = graph.get_node(lastNodeID);

                    Node n = startNodes.second;
                    graph.create_edge(lastNode, &n);
                }
            }
        } else if (variationType == INVERSION){
            //Inversion
            size_t length = lengthDistribution(generator);
            if (length > 0) {
                if (startIndex == 0) {

                    pair<Node, Node> endNodes = splitNode(startIndex+length);

                    Node end = endNodes.second;
                    Node n = endNodes.first;
                    graph.create_edge(&n, &end, true, false);


                } else if ( length + startIndex < seqSize-1) {
                    pair<Node, Node> startNodes = splitNode(startIndex);
                    pair<Node, Node> endNodes = splitNode(startIndex+length);

                    Node start = startNodes.first;
                    Node end = endNodes.second;
                    Node n1 = startNodes.second;
                    Node n2 = endNodes.first;
                    graph.create_edge(&start, &n2, false, true);
                    graph.create_edge(&n1, &end, true, false);

               } else if (length + startIndex == seqSize - 1) {

                    pair<Node, Node> startNodes = splitNode(startIndex);

                    Node start = startNodes.first;
                    Node n = startNodes.second;
                    graph.create_edge(&start, &n, false, true);


               }
           }
        }
    }
    return graph;
};

}
}
