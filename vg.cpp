#include "vg.h"

using namespace std;
using namespace google::protobuf;
using namespace vg;


//VariantGraph::VariantGraph(void) { };
// construct from protobufs
VariantGraph::VariantGraph(ifstream& in) {
    ParseFromIstream(&in);
    // populate by-id node index
    for (int64_t i = 0; i < nodes_size(); ++i) {
        Node* n = mutable_nodes(i);
        node_index[n] = i;
        node_by_id[n->id()] = n;
    }
}

VariantGraph::VariantGraph(vector<Node>& nodesv) {
    for (vector<Node>::iterator n = nodesv.begin(); n != nodesv.end(); ++n) {
        Node& node = *n;
        int64_t id = node.id();
        if (current_id() < id) {
            set_current_id(id);
        }
        Node* new_node = add_nodes(); // add it to the graph
        new_node->set_sequence(node.sequence());
        new_node->set_id(node.id());
        //Node& new_node = nodes(nodes_size()-1); // get a reference to it
        node_by_id[new_node->id()] = new_node; // and insert into our id lookup table
        node_index[new_node] = nodes_size()-1;
    }
}

// construct from VCF records
// --------------------------
// algorithm
// maintain a core reference path upon which we add new variants as they come
// addition procedure is the following
// find reference node overlapping our start position
// if it is already the end of a node, add the new node
// if it is not the end of a node, break it, insert edges from old->new
// go to end position of alt allele (could be the same position)
// if it already has a break, just point to the next node in line
// if it is not broken, break it and point to the next node
// add new node for alt alleles, connect to start and end node in reference path
// store the ref mapping as a property of the edges and nodes (this allows deletion edges and insertion subpaths)
//
VariantGraph::VariantGraph(vcf::VariantCallFile& variantCallFile, FastaReference& reference) {
//// .... XXX

    //cerr << "target sequence " << targetSequence << endl;
    //cerr << "starts at " << offset << endl;

    for (vector<string>::iterator r = reference.index->sequenceNames.begin();
         r != reference.index->sequenceNames.end(); ++r) {

        string& seqName = *r;
        map<long, Node*> reference_path;
        //map<long, set<Node*> > nodes; // for maintaining a reference-sorted graph
        string seq = reference.getSequence(seqName);

        Node* ref_node = create_node(seq);
        reference_path[0] = ref_node;

        variantCallFile.setRegion(seqName);
        vcf::Variant var(variantCallFile);
        while (variantCallFile.getNextVariant(var)) {

            int current_pos = (long int) var.position - 1;
            // decompose the alt
            bool flat_input_vcf = false; // hack
            map<string, vector<vcf::VariantAllele> > alternates = (flat_input_vcf ? var.flatAlternates() : var.parsedAlternates());
            for (map<string, vector<vcf::VariantAllele> >::iterator va = alternates.begin(); va !=alternates.end(); ++va) {
                vector<vcf::VariantAllele>& alleles = va->second;

                for (vector<vcf::VariantAllele>::iterator a = alleles.begin(); a != alleles.end(); ++a) {
                    vcf::VariantAllele& allele = *a;

                    // reference alleles are provided naturally by the reference itself
                    if (allele.ref == allele.alt) {
                        continue;
                    }

                    long allele_start_pos = allele.position - 1;  // 0/1 based conversion... thanks vcflib!
                    long allele_end_pos = allele_start_pos + allele.ref.size();

                    Node* left_ref_node = NULL;
                    Node* middle_ref_node = NULL;
                    Node* right_ref_node = NULL;

                    // divide_path(map<long, Node*>& path, long pos, Node*& left, Node*& right) {
                    divide_path(reference_path,
                                allele_start_pos,
                                left_ref_node,
                                right_ref_node);

                    // if the ref portion of the allele is not empty, then we need to make another cut
                    if (!allele.ref.empty()) {
                        divide_path(reference_path,
                                    allele_start_pos,
                                    middle_ref_node,
                                    right_ref_node);
                    }

                    // create a new alt node and connect the pieces from before
                    if (!allele.alt.empty()) {
                        Node* alt_node = create_node(allele.alt);
                        //ref_map.add_node(alt_node, allele_start_pos, );
                        add_edge(left_ref_node, alt_node);
                        add_edge(alt_node, right_ref_node);
                    } else {// otherwise, we have a deletion
                        add_edge(left_ref_node, right_ref_node);
                    }
                }
            }
        }
    }
}

Edge* VariantGraph::create_edge(int64_t from, int64_t to) {
    Edge* edge = add_edges();
    edge->set_prev(from);
    edge->set_next(to);
    edge_by_ids[make_pair(from, to)] = edge;
    edge_index[edge] = edges_size()-1;
    return edge;
}

void VariantGraph::remove_edge_from_node(Node* node, Edge* edge) {
    int64_t prev_id = edge->prev();
    for (int i = 0; i < node->prev_size(); ++i) {
        if (node->prev(i) == prev_id) {
            node->mutable_prev()->SwapElements(i, node->prev_size()-1);
            node->mutable_prev()->RemoveLast();
            return;
        }
    }
    int64_t next_id = edge->next();
    for (int i = 0; i < node->next_size(); ++i) {
        if (node->next(i) == next_id) {
            node->mutable_next()->SwapElements(i, node->next_size()-1);
            node->mutable_next()->RemoveLast();
            return;
        }
    }
}

void VariantGraph::destroy_edge(Edge* edge) {
    int lei = edges_size()-1;
    int tei = edge_index[edge];
    remove_edge_from_node(node_by_id[edge->prev()], edge);
    remove_edge_from_node(node_by_id[edge->next()], edge);
    Edge* last = mutable_edges(lei);
    mutable_edges()->SwapElements(tei, lei);
    Edge* elast = mutable_edges(tei);
    edge_by_ids[make_pair(last->prev(), last->next())] = elast;
    edge_index.erase(last);
    edge_index[elast] = tei;
    edge_index.erase(edge);
    mutable_edges()->RemoveLast();
}

// use the VariantGraph class to generate ids
Node* VariantGraph::create_node(string seq) {
    // create the node
    Node* node = add_nodes();
    node->set_sequence(seq);
    node->set_id(current_id());
    set_current_id(current_id()+1);
    // copy it into the graph
    // and drop into our id index
    node_by_id[node->id()] = node;
    node_index[node] = nodes_size()-1;
    return node;
}

void VariantGraph::destroy_node(Node* node) {
    // swap node with the last in nodes
    // call RemoveLast() to drop the node
    int lni = nodes_size()-1;
    int tni = node_index[node];
    Node* last = mutable_nodes(lni);
    mutable_nodes()->SwapElements(tni, lni);
    Node* nlast = mutable_nodes(tni);
    node_by_id[last->id()] = nlast;
    node_index.erase(last);
    node_index[nlast] = tni;
    node_index.erase(node);
    mutable_nodes()->RemoveLast();
}

// utilities
void VariantGraph::divide_node(Node* node, int pos, Node*& left, Node*& right) {

    // make our left node
    left = create_node(node->sequence().substr(0,pos));

    // replace node connections to prev (left)
    for (int i = 0; i < node->prev_size(); ++i) {
        Node* p = node_by_id[node->prev(i)];
        node_replace_next(p, node, left);
    }

    // make our right node
    right = create_node(node->sequence().substr(pos,node->sequence().size()-1));

    // replace node connections to next (right)
    for (int i = 0; i < node->next_size(); ++i) {
        Node* n = node_by_id[node->next(i)];
        node_replace_prev(n, node, right);
    }

    // connect left to right
    add_edge(left, right);

    destroy_node(node);
}

Edge* VariantGraph::add_edge(Node* from, Node* to) {
    from->add_next(to->id());
    to->add_prev(from->id());
    return create_edge(from->id(), to->id());
}

// for dividing a path of nodes with an underlying coordinate system
void VariantGraph::divide_path(map<long, Node*>& path, long pos, Node*& left, Node*& right) {

    map<long, Node*>::iterator target = path.upper_bound(pos);
    --target; // we should now be pointing to the target ref node

    long node_pos = target->first;
    Node* old = target->second;
    
    // nothing to do
    if (node_pos == pos) {
        map<long, Node*>::iterator n = target; --n;
        left = n->second;
        right = target->second;

    } else {

        // divide the target node at our pos
        int diff = pos - node_pos;

        divide_node(old, diff, left, right);

        // left
        path[node_pos] = left;

        // right
        path[pos] = right;
    }
}

void VariantGraph::node_replace_prev(Node* node, Node* before, Node* after) {
    destroy_edge(edge_by_ids[make_pair(before->id(), node->id())]);
    add_edge(after, node);
}

void VariantGraph::node_replace_next(Node* node, Node* before, Node* after) {
    destroy_edge(edge_by_ids[make_pair(node->id(), before->id())]);
    add_edge(node, after);
}

//void align(Alignment& alignment);
