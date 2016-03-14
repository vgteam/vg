#include "vg.hpp"
#include "stream.hpp"

namespace vg {

using namespace std;
using namespace gfak;
using namespace supbub;


// construct from a stream of protobufs
VG::VG(istream& in, bool showp) {

    // set up uninitialized values
    init();
    show_progress = showp;
    // and if we should show progress
    function<void(uint64_t)> handle_count = [this](uint64_t count) {
        create_progress("loading graph", count);
    };

    // the graph is read in chunks, which are attached to this graph
    uint64_t i = 0;
    function<void(Graph&)> lambda = [this, &i](Graph& g) {
        update_progress(++i);
        // We expect these to not overlap in nodes or edges, so complain if they do.
        extend(g, true);
    };

    stream::for_each(in, lambda, handle_count);

    // store paths in graph
    paths.to_graph(graph);

    destroy_progress();

}

// construct from an arbitrary source of Graph protobuf messages
VG::VG(function<bool(Graph&)>& get_next_graph, bool showp) {
    // set up uninitialized values
    init();
    show_progress = showp;

    // We can't show loading progress since we don't know the total number of
    // subgraphs.

    // Try to load the first graph
    Graph subgraph;
    bool got_subgraph = get_next_graph(subgraph);
    while(got_subgraph) {
        // If there is a valid subgraph, add it to ourselves.
        // We expect these to not overlap in nodes or edges, so complain if they do.
        extend(subgraph, true);
        // Try and load the next subgraph, if it exists.
        got_subgraph = get_next_graph(subgraph);
    }

    // store paths in graph
    paths.to_graph(graph);
}

void VG::clear_paths(void) {
    paths.clear();
    graph.clear_path(); // paths.clear() should do this too
    sync_paths();
}

// synchronize the VG index and its backing store
void VG::sync_paths(void) {
    // ensure we can navigate paths correctly
    // by building paths.mapping_path_order
    paths.rebuild_mapping_aux();
}

void VG::serialize_to_ostream(ostream& out, id_t chunk_size) {

    sync_paths();

    // save the number of the messages to be serialized into the output file
    uint64_t count = graph.node_size() / chunk_size + 1;
    create_progress("saving graph", count);
    // partition the graph into a number of chunks (required by format)
    // constructing subgraphs and writing them to the stream
    function<Graph(uint64_t)> lambda =
        [this, chunk_size](uint64_t i) -> Graph {
        VG g;
        map<string, map<size_t, Mapping*> > sorted_paths;
        for (id_t j = i * chunk_size;
             j < (i+1)*chunk_size && j < graph.node_size();
             ++j) {
            Node* node = graph.mutable_node(j);
            // Grab the node and only the edges where it has the lower ID.
            // This prevents duplication of edges in the serialized output.
            nonoverlapping_node_context_without_paths(node, g);
            //set<pair<string, Mapping*> >& Paths::get_node_mapping(id_t id);
            auto& mappings = paths.get_node_mapping(node);
            //cerr << "getting node mappings for " << node->id() << endl;
            for (auto m : mappings) {
                auto& name = m.first;
                auto& mappings = m.second;
                for (auto& mapping : mappings) {
                    //cerr << "mapping " << name << pb2json(*mapping) << endl;
                    sorted_paths[name][paths.mapping_path_order[mapping]] = mapping;
                }
            }
        }
        // now get the paths for this chunk so that they are ordered correctly
        for (auto& p : sorted_paths) {
            auto& name = p.first;
            auto& path = p.second;
            // now sorted in ascending order
            // we could also assert that we have a contiguous path here
            for (auto& m : path) {
                g.paths.append_mapping(name, *m.second);
            }
        }

        // but this is broken as our paths have been reordered as
        // the nodes they cross are stored in graph.nodes
        g.paths.to_graph(g.graph);

        update_progress(i);
        return g.graph;
    };

    stream::write(out, count, lambda);

    destroy_progress();
}

void VG::serialize_to_file(const string& file_name, id_t chunk_size) {
    ofstream f(file_name);
    serialize_to_ostream(f);
    f.close();
}

VG::~VG(void) {
    destroy_alignable_graph();
}

VG::VG(void) {
    init();
}

void VG::init(void) {
    gssw_aligner = NULL;
    current_id = 1;
    show_progress = false;
    progress_message = "progress";
    progress = NULL;
}

VG::VG(set<Node*>& nodes, set<Edge*>& edges) {
    init();
    add_nodes(nodes);
    add_edges(edges);
    sort();
}


SB_Input VG::vg_to_sb_input(){
	//cout << this->edge_count() << endl;
  SB_Input sbi;
  sbi.num_vertices = this->edge_count();
	function<void(Edge*)> lambda = [&sbi](Edge* e){
		//cout << e->from() << " " << e->to() << endl;
    pair<id_t, id_t> dat = make_pair(e->from(), e->to() );
    sbi.edges.push_back(dat);
	};
	this->for_each_edge(lambda);
  return sbi;
}

vector<pair<id_t, id_t> > VG::get_superbubbles(SB_Input sbi){
    vector<pair<id_t, id_t> > ret;
    supbub::Graph sbg (sbi.num_vertices);
    supbub::DetectSuperBubble::SUPERBUBBLE_LIST superBubblesList{};
    supbub::DetectSuperBubble dsb;
    dsb.find(sbg, superBubblesList);
    supbub::DetectSuperBubble::SUPERBUBBLE_LIST::iterator it;
    for (it = superBubblesList.begin(); it != superBubblesList.end(); ++it) {
        ret.push_back(make_pair((*it).entrance, (*it).exit));
    }
    return ret;
}
vector<pair<id_t, id_t> > VG::get_superbubbles(){
    vector<pair<id_t, id_t> > ret;
    supbub::Graph sbg (this->edge_count());
    //load up the sbgraph with edges
    function<void(Edge*)> lambda = [&sbg](Edge* e){
            //cout << e->from() << " " << e->to() << endl;
        sbg.addEdge(e->from(), e->to());
    };

    this->for_each_edge(lambda);

    supbub::DetectSuperBubble::SUPERBUBBLE_LIST superBubblesList{};

    supbub::DetectSuperBubble dsb;
    dsb.find(sbg, superBubblesList);
    supbub::DetectSuperBubble::SUPERBUBBLE_LIST::iterator it;
    for (it = superBubblesList.begin(); it != superBubblesList.end(); ++it) {
        ret.push_back(make_pair((*it).entrance, (*it).exit));
    }
    return ret;
}
// check for conflict (duplicate nodes and edges) occurs within add_* functions

void VG::add_nodes(set<Node*>& nodes) {
    for (auto node : nodes) {
        add_node(*node);
    }
}

void VG::add_edges(set<Edge*>& edges) {
    for (auto edge : edges) {
        add_edge(*edge);
    }
}

void VG::add_nodes(vector<Node>& nodes) {
    for (auto& node : nodes) {
        add_node(node);
    }
}

void VG::add_edges(vector<Edge>& edges) {
    for (auto& edge : edges) {
        add_edge(edge);
    }
}

void VG::add_node(Node& node) {
    if (!has_node(node)) {
        Node* new_node = graph.add_node(); // add it to the graph
        *new_node = node; // overwrite it with the value of the given node
        node_by_id[new_node->id()] = new_node; // and insert into our id lookup table
        node_index[new_node] = graph.node_size()-1;
    }
}

void VG::add_edge(Edge& edge) {
    if (!has_edge(edge)) {
        Edge* new_edge = graph.add_edge(); // add it to the graph
        *new_edge = edge;
        set_edge(new_edge);
        edge_index[new_edge] = graph.edge_size()-1;
    }
}

id_t VG::node_count(void) {
    return graph.node_size();
}

id_t VG::edge_count(void) {
    return graph.edge_size();
}

vector<pair<id_t, bool>>& VG::edges_start(Node* node) {
    if(node == nullptr) {
        return empty_edge_ends;
    }
    return edges_start(node->id());
}

vector<pair<id_t, bool>>& VG::edges_start(id_t id) {
    if(edges_on_start.count(id) == 0) {
        return empty_edge_ends;
    }
    return edges_on_start[id];
}

vector<pair<id_t, bool>>& VG::edges_end(Node* node) {
    if(node == nullptr) {
        return empty_edge_ends;
    }
    return edges_end(node->id());
}

vector<pair<id_t, bool>>& VG::edges_end(id_t id) {
    if(edges_on_end.count(id) == 0) {
        return empty_edge_ends;
    }
    return edges_on_end[id];
}

int VG::start_degree(Node* node) {
    return edges_start(node).size();
}

int VG::end_degree(Node* node) {
    return edges_end(node).size();
}

int VG::left_degree(NodeTraversal node) {
    // If we're backward, the end is on the left. Otherwise, the start is.
    return node.backward ? end_degree(node.node) : start_degree(node.node);
}

int VG::right_degree(NodeTraversal node) {
    // If we're backward, the start is on the right. Otherwise, the end is.
    return node.backward ? start_degree(node.node) : end_degree(node.node);
}

void VG::edges_of_node(Node* node, vector<Edge*>& edges) {
    for(pair<id_t, bool>& off_start : edges_start(node)) {
        // Go through the edges on this node's start
        Edge* edge = edge_by_sides[NodeSide::pair_from_start_edge(node->id(), off_start)];
        if (!edge) {
            cerr << "error:[VG::edges_of_node] nonexistent start edge " << off_start.first << " start <-> "
                 << node->id() << (off_start.second ? " start" : " end") << endl;
            exit(1);
        }
        edges.push_back(edge);
    }

    for(pair<id_t, bool>& off_end : edges_end(node)) {
        // And on its end
        Edge* edge = edge_by_sides[NodeSide::pair_from_end_edge(node->id(), off_end)];
        if (!edge) {
            cerr << "error:[VG::edges_of_node] nonexistent end edge " << off_end.first << " end <-> "
                 << node->id() << (off_end.second ? " end" : " start") << endl;
            exit(1);
        }
        if(edge->from() == edge->to() && edge->from_start() == edge->to_end()) {
            // This edge touches both out start and our end, so we already
            // handled it on our start. Don't produce it twice.
            continue;
        }
        edges.push_back(edge);
    }
}

vector<Edge*> VG::edges_of(Node* node) {
    vector<Edge*> edges;
    edges_of_node(node, edges);
    return edges;
}

void VG::edges_of_nodes(set<Node*>& nodes, set<Edge*>& edges) {
    for (set<Node*>::iterator n = nodes.begin(); n != nodes.end(); ++n) {
        vector<Edge*> ev;
        edges_of_node(*n, ev);
        for (vector<Edge*>::iterator e = ev.begin(); e != ev.end(); ++e) {
            edges.insert(*e);
        }
    }
}

set<pair<NodeSide, bool>> VG::sides_context(id_t node_id) {
    // return the side we're going to and if we go from the start or end to get there
    set<pair<NodeSide, bool>> all;
    for (auto& s : sides_to(NodeSide(node_id, false))) {
        all.insert(make_pair(s, false));
    }
    for (auto& s : sides_to(NodeSide(node_id, true))) {
        all.insert(make_pair(s, true));
    }
    for (auto& s : sides_from(NodeSide(node_id, false))) {
        all.insert(make_pair(s, false));
    }
    for (auto& s : sides_from(NodeSide(node_id, true))) {
        all.insert(make_pair(s, true));
    }
    return all;
}

bool VG::same_context(id_t n1, id_t n2) {
    auto c1 = sides_context(n1);
    auto c2 = sides_context(n2);
    bool same = true;
    for (auto& s : c1) {
        if (!c2.count(s)) { same = false; break; }
    }
    return same;
}

bool VG::is_ancestor_prev(id_t node_id, id_t candidate_id) {
    set<id_t> seen;
    return is_ancestor_prev(node_id, candidate_id, seen);
}

bool VG::is_ancestor_prev(id_t node_id, id_t candidate_id, set<id_t>& seen, size_t steps) {
    if (node_id == candidate_id) return true;
    if (!steps) return false;
    for (auto& side : sides_to(NodeSide(node_id, false))) {
        if (seen.count(side.node)) continue;
        seen.insert(side.node);
        if (is_ancestor_prev(side.node, candidate_id, seen, steps-1)) return true;
    }
    return false;
}

bool VG::is_ancestor_next(id_t node_id, id_t candidate_id) {
    set<id_t> seen;
    return is_ancestor_next(node_id, candidate_id, seen);
}

bool VG::is_ancestor_next(id_t node_id, id_t candidate_id, set<id_t>& seen, size_t steps) {
    if (node_id == candidate_id) return true;
    if (!steps) return false;
    for (auto& side : sides_from(NodeSide(node_id, true))) {
        if (seen.count(side.node)) continue;
        seen.insert(side.node);
        if (is_ancestor_next(side.node, candidate_id, seen, steps-1)) return true;
    }
    return false;
}

id_t VG::common_ancestor_prev(id_t id1, id_t id2, size_t steps) {
    // arbitrarily step back from node 1 asking if we are prev-ancestral to node 2
    auto scan = [this](id_t id1, id_t id2, size_t steps) -> id_t {
        set<id_t> to_visit;
        to_visit.insert(id1);
        for (size_t i = 0; i < steps; ++i) {
            // collect nodes to visit
            set<id_t> to_visit_next;
            for (auto& id : to_visit) {
                if (is_ancestor_prev(id2, id)) return id;
                for (auto& side : sides_to(NodeSide(id, false))) {
                    to_visit_next.insert(side.node);
                }
            }
            to_visit = to_visit_next;
            if (to_visit.empty()) return -1; // we hit the end of the graph
        }
        return 0;
    };
    id_t id3 = scan(id1, id2, steps);
    if (id3) {
        return id3;
    } else {
        return scan(id2, id1, steps);
    }
}

id_t VG::common_ancestor_next(id_t id1, id_t id2, size_t steps) {
    // arbitrarily step forward from node 1 asking if we are next-ancestral to node 2
    auto scan = [this](id_t id1, id_t id2, size_t steps) -> id_t {
        set<id_t> to_visit;
        to_visit.insert(id1);
        for (size_t i = 0; i < steps; ++i) {
            // collect nodes to visit
            set<id_t> to_visit_next;
            for (auto& id : to_visit) {
                if (is_ancestor_next(id2, id)) return id;
                for (auto& side : sides_from(NodeSide(id, true))) {
                    to_visit_next.insert(side.node);
                }
            }
            to_visit = to_visit_next;
            if (to_visit.empty()) return -1; // we hit the end of the graph
        }
        return 0;
    };
    id_t id3 = scan(id1, id2, steps);
    if (id3) {
        return id3;
    } else {
        return scan(id2, id1, steps);
    }
}

set<NodeSide> VG::sides_of(NodeSide side) {
    set<NodeSide> v1 = sides_to(side);
    set<NodeSide> v2 = sides_from(side);
    for (auto s : v2) v1.insert(s);
    return v1;
}

set<NodeSide> VG::sides_to(NodeSide side) {
    set<NodeSide> other_sides;
    vector<Edge*> edges;
    edges_of_node(get_node(side.node), edges);
    for (auto* edge : edges) {
        if (edge->to() == side.node && edge->to_end() == side.is_end) {
            other_sides.insert(NodeSide(edge->from(), !edge->from_start()));
        }
    }
    return other_sides;
}

set<NodeSide> VG::sides_from(NodeSide side) {
    set<NodeSide> other_sides;
    vector<Edge*> edges;
    edges_of_node(get_node(side.node), edges);
    for (auto* edge : edges) {
        if (edge->from() == side.node && edge->from_start() != side.is_end) {
            other_sides.insert(NodeSide(edge->to(), edge->to_end()));
        }
    }
    return other_sides;
}

set<NodeSide> VG::sides_from(id_t id) {
    set<NodeSide> sides;
    for (auto side : sides_from(NodeSide(id, true))) {
        sides.insert(side);
    }
    for (auto side : sides_from(NodeSide(id, false))) {
        sides.insert(side);
    }
    return sides;
}

set<NodeSide> VG::sides_to(id_t id) {
    set<NodeSide> sides;
    for (auto side : sides_to(NodeSide(id, true))) {
        sides.insert(side);
    }
    for (auto side : sides_to(NodeSide(id, false))) {
        sides.insert(side);
    }
    return sides;
}

set<NodeTraversal> VG::siblings_to(const NodeTraversal& trav) {
    // find the sides to
    auto to_sides = sides_to(NodeSide(trav.node->id(), trav.backward));
    // and then find the traversals from them
    set<NodeTraversal> travs_from_to_sides;
    for (auto& s1 : to_sides) {
        // and the from-children of these
        for (auto& s2 : sides_from(s1)) {
            auto sib = NodeTraversal(get_node(s2.node), s2.is_end);
            // which are not this node
            if (sib != trav) {
                travs_from_to_sides.insert(sib);
            }
        }
    }
    return travs_from_to_sides;
}

set<NodeTraversal> VG::siblings_from(const NodeTraversal& trav) {
    // find the sides from
    auto from_sides = sides_from(NodeSide(trav.node->id(), !trav.backward));
    // and then find the traversals from them
    set<NodeTraversal> travs_to_from_sides;
    for (auto& s1 : from_sides) {
        // and the to-children of these
        for (auto& s2 : sides_to(s1)) {
            auto sib = NodeTraversal(get_node(s2.node), !s2.is_end);
            // which are not this node
            if (sib != trav) {
                travs_to_from_sides.insert(sib);
            }
        }
    }
    return travs_to_from_sides;
}

set<Node*> VG::siblings_of(Node* node) {
    set<Node*> sibs;
    for (auto& s : siblings_to(NodeTraversal(node, false))) {
        sibs.insert(s.node);
    }
    for (auto& s : siblings_to(NodeTraversal(node, true))) {
        sibs.insert(s.node);
    }
    for (auto& s : siblings_from(NodeTraversal(node, false))) {
        sibs.insert(s.node);
    }
    for (auto& s : siblings_from(NodeTraversal(node, true))) {
        sibs.insert(s.node);
    }
    return sibs;
}

set<NodeTraversal> VG::full_siblings_to(const NodeTraversal& trav) {
    // get the siblings of
    auto sibs_to = siblings_to(trav);
    // and filter them for nodes with the same inbound sides
    auto to_sides = sides_to(NodeSide(trav.node->id(), trav.backward));
    set<NodeTraversal> full_sibs_to;
    for (auto& sib : sibs_to) {
        auto sib_to_sides = sides_to(NodeSide(sib.node->id(), sib.backward));
        if (sib_to_sides == to_sides) {
            full_sibs_to.insert(sib);
        }
    }
    return full_sibs_to;
}

set<NodeTraversal> VG::full_siblings_from(const NodeTraversal& trav) {
    // get the siblings of
    auto sibs_from = siblings_from(trav);
    // and filter them for nodes with the same outbound sides
    auto from_sides = sides_from(NodeSide(trav.node->id(), !trav.backward));
    set<NodeTraversal> full_sibs_from;
    for (auto& sib : sibs_from) {
        auto sib_from_sides = sides_from(NodeSide(sib.node->id(), !sib.backward));
        if (sib_from_sides == from_sides) {
            full_sibs_from.insert(sib);
        }
    }
    return full_sibs_from;
}

// returns sets of sibling nodes that are only in one set of sibling nodes
set<set<NodeTraversal>> VG::transitive_sibling_sets(const set<set<NodeTraversal>>& sibs) {
    set<set<NodeTraversal>> trans_sibs;
    map<Node*, int> membership;
    // determine the number of sibling sets that each node is in
    for (auto& s : sibs) {
        for (auto& t : s) {
            if (membership.find(t.node) == membership.end()) {
                membership[t.node] = 1;
            } else {
                ++membership[t.node];
            }
        }
    }
    // now exclude components which are intransitive
    // by only keeping those sib sets whose members are in only one set
    for (auto& s : sibs) {
        // all members must only appear in this set
        bool is_transitive = true;
        for (auto& t : s) {
            if (membership[t.node] > 1) {
                is_transitive = false;
                break;
            }
        }
        if (is_transitive) {
            trans_sibs.insert(s);
        }
    }
    return trans_sibs;
}

set<set<NodeTraversal>> VG::identically_oriented_sibling_sets(const set<set<NodeTraversal>>& sibs) {
    set<set<NodeTraversal>> iosibs;
    for (auto& s : sibs) {
        int forward = 0;
        int reverse = 0;
        for (auto& t : s) {
            if (t.backward) {
                ++reverse;
            } else {
                ++forward;
            }
        }
        // if they are all forward or all reverse
        if (forward == 0 || reverse == 0) {
            iosibs.insert(s);
        }
    }
    return iosibs;
}

void VG::simplify_siblings(void) {
    // make a list of all the sets of siblings
    set<set<NodeTraversal>> to_sibs;
    for_each_node([this, &to_sibs](Node* n) {
            auto trav = NodeTraversal(n, false);
            auto tsibs = full_siblings_to(trav);
            tsibs.insert(trav);
            if (tsibs.size() > 1) {
                to_sibs.insert(tsibs);
            }
        });
    // make the sibling sets transitive
    // by removing any that are intransitive
    // then simplify
    simplify_to_siblings(
        identically_oriented_sibling_sets(
            transitive_sibling_sets(to_sibs)));
    // and remove any null nodes that result
    remove_null_nodes_forwarding_edges();

    // make a list of the from-siblings
    set<set<NodeTraversal>> from_sibs;
    for_each_node([this, &from_sibs](Node* n) {
            auto trav = NodeTraversal(n, false);
            auto fsibs = full_siblings_from(trav);
            fsibs.insert(trav);
            if (fsibs.size() > 1) {
                from_sibs.insert(fsibs);
            }
        });
    // then do the from direction
    simplify_from_siblings(
        identically_oriented_sibling_sets(
            transitive_sibling_sets(from_sibs)));
    // and remove any null nodes that result
    remove_null_nodes_forwarding_edges();

}

void VG::simplify_to_siblings(const set<set<NodeTraversal>>& to_sibs) {
    for (auto& sibs : to_sibs) {
        // determine the amount of sharing at the start
        // the to-sibs have the same parent(s) feeding into them
        // so we can safely make a single node out of the shared sequence
        // and link this to them and their parent to remove node level redundancy
        vector<string*> seqs;
        size_t min_seq_size = sibs.begin()->node->sequence().size();
        for (auto& sib : sibs) {
            auto seqp = sib.node->mutable_sequence();
            seqs.push_back(seqp);
            if (seqp->size() < min_seq_size) {
                min_seq_size = seqp->size();
            }
        }
        size_t i = 0;
        size_t j = 0;
        bool similar = true;
        for ( ; similar && i < min_seq_size; ++i) {
            //cerr << i << endl;
            char c = seqs.front()->at(i);
            for (auto s : seqs) {
                //cerr << "checking " << c << " vs " << s->at(i) << endl;
                if (c != s->at(i)) {
                    similar = false;
                    break;
                }
            }
            if (!similar) break;
            ++j;
        }
        size_t shared_start = j;
        //cerr << "sharing is " << shared_start << " for to-sibs of "
        //<< sibs.begin()->node->id() << endl;
        if (shared_start == 0) continue;

        // make a new node with the shared sequence
        string seq = seqs.front()->substr(0,shared_start);
        auto new_node = create_node(seq);

        // remove the sequence of the new node from the old nodes
        for (auto& sib : sibs) {
            //cerr << "to sib " << pb2json(*sib.node) << endl;
            *sib.node->mutable_sequence() = sib.node->sequence().substr(shared_start);
            // for each node mapping of the sibling
            // divide the mapping at the cut point

            // and then switch the node assignment for the cut nodes
            // for each mapping of the node
            for (auto& p : paths.get_node_mapping(sib.node)) {
                vector<Mapping*> v;
                for (auto& m : p.second) {
                    v.push_back(m);
                }
                for (auto m : v) {
                    auto mpts = paths.divide_mapping(m, shared_start);
                    // and then assign the first part of the mapping to the new node
                    auto o = mpts.first;
                    o->mutable_position()->set_offset(0);
                    auto n = mpts.second;
                    n->mutable_position()->set_offset(0);
                    paths.reassign_node(new_node->id(), n);
                    // note that the other part now maps to the correct (old) node
                }
            }
        }

        // connect the new node to the common *context* (the union of sides of the old nodes)

        // by definition we are only working with nodes that have exactly the same set of parents
        // so we just use the first node in the set to drive the reconnection
        auto new_left_side = NodeSide(new_node->id(), false);
        auto new_right_side = NodeSide(new_node->id(), true);
        for (auto side : sides_to(NodeSide(sibs.begin()->node->id(), sibs.begin()->backward))) {
            create_edge(side, new_left_side);
        }
        // disconnect the old nodes from their common parents
        for (auto& sib : sibs) {
            auto old_side = NodeSide(sib.node->id(), sib.backward);
            for (auto side : sides_to(old_side)) {
                destroy_edge(side, old_side);
            }
            // connect the new node to the old nodes
            create_edge(new_right_side, old_side);
        }
    }
}

void VG::simplify_from_siblings(const set<set<NodeTraversal>>& from_sibs) {
    for (auto& sibs : from_sibs) {
        // determine the amount of sharing at the end
        // the from-sibs have the same downstream nodes ("parents")
        // so we can safely make a single node out of the shared sequence at the end
        // and link this to them and their parent to remove node level redundancy
        vector<string*> seqs;
        size_t min_seq_size = sibs.begin()->node->sequence().size();
        for (auto& sib : sibs) {
            auto seqp = sib.node->mutable_sequence();
            seqs.push_back(seqp);
            if (seqp->size() < min_seq_size) {
                min_seq_size = seqp->size();
            }
        }
        size_t i = 0;
        size_t j = 0;
        bool similar = true;
        for ( ; similar && i < min_seq_size; ++i) {
            char c = seqs.front()->at(seqs.front()->size()-(i+1));
            for (auto s : seqs) {
                if (c != s->at(s->size()-(i+1))) {
                    similar = false;
                    break;
                }
            }
            if (!similar) break;
            ++j;
        }
        size_t shared_end = j;
        if (shared_end == 0) continue;
        // make a new node with the shared sequence
        string seq = seqs.front()->substr(seqs.front()->size()-shared_end);
        auto new_node = create_node(seq);
        // chop it off of the old nodes
        for (auto& sib : sibs) {
            *sib.node->mutable_sequence()
                = sib.node->sequence().substr(0, sib.node->sequence().size()-shared_end);

            // and then switch the node assignment for the cut nodes
            // for each mapping of the node
            for (auto& p : paths.get_node_mapping(sib.node)) {
                vector<Mapping*> v;
                for (auto& m : p.second) {
                    v.push_back(m);
                }
                for (auto m : v) {
                    auto mpts = paths.divide_mapping(m, sib.node->sequence().size());
                    // and then assign the second part of the mapping to the new node
                    auto o = mpts.first;
                    o->mutable_position()->set_offset(0);
                    paths.reassign_node(new_node->id(), o);
                    auto n = mpts.second;
                    n->mutable_position()->set_offset(0);

                    // note that the other part now maps to the correct (old) node
                }
            }
        }
        // connect the new node to the common downstream nodes
        // by definition we are only working with nodes that have exactly the same set of "children"
        // so we just use the first node in the set to drive the reconnection
        auto new_left_side = NodeSide(new_node->id(), false);
        auto new_right_side = NodeSide(new_node->id(), true);
        for (auto side : sides_from(NodeSide(sibs.begin()->node->id(), !sibs.begin()->backward))) {
            create_edge(new_right_side, side);
        }
        // disconnect the old nodes from their common "children"
        for (auto& sib : sibs) {
            auto old_side = NodeSide(sib.node->id(), !sib.backward);
            for (auto side : sides_from(old_side)) {
                destroy_edge(old_side, side);
            }
            // connect the new node to the old nodes
            create_edge(old_side, new_left_side);
        }
    }
}

// expand the context of the subgraph g by this many steps
// it's like a neighborhood function
void VG::expand_context(VG& g, size_t steps, bool add_paths) {
    set<id_t> to_visit;
    // start with the nodes in the subgraph
    g.for_each_node([&](Node* n) { to_visit.insert(n->id()); });
    g.for_each_edge([&](Edge* e) {
            to_visit.insert(e->from());
            to_visit.insert(e->to()); });
    // and expand
    for (size_t i = 0; i < steps; ++i) {
        set<id_t> to_visit_next;
        for (auto id : to_visit) {
            // build out the graph
            // if we have nodes we haven't seeen
            if (!g.has_node(id)) {
                g.add_node(*get_node(id));
            }
            for (auto& e : edges_of(get_node(id))) {
                g.add_edge(*e);
                if (e->from() == id) {
                    to_visit_next.insert(e->to());
                } else {
                    to_visit_next.insert(e->from());
                }
            }
        }
        to_visit = to_visit_next;
    }
    // then remove orphans
    g.remove_orphan_edges();
    // and add paths
    if (add_paths) {
        g.for_each_node([&](Node* n) {
                for (auto& path : paths.get_node_mapping(n)) {
                    for (auto& m : path.second) {
                        g.paths.append_mapping(path.first, *m);
                    }
                }
            });
        g.sync_paths();
    }
}

bool VG::adjacent(const Position& pos1, const Position& pos2) {
    // two positions are on the same node
    if (pos1.node_id() == pos2.node_id()) {
        if (pos1.offset() == pos1.offset()+1) {
            // and have adjacent offsets
            return true;
        } else {
            // if not, they aren't adjacent
            return false;
        }
    } else {
        // is the first at the end of its node
        // and the second at the start of its node
        // determine if the two nodes are connected
        auto* node1 = get_node(pos1.node_id());
        auto* node2 = get_node(pos2.node_id());
        if (pos1.offset() == node1->sequence().size()-1
            && pos2.offset() == 0) {
            // these are adjacent iff we have an edge
            return has_edge(NodeSide(pos1.node_id(), true),
                            NodeSide(pos2.node_id(), false));
        } else {
            // the offsets aren't at the end and start
            // so these positions can't be adjacent
            return false;
        }
    }
}

// edges which are both from_start and to_end can be represented naturally as
// a regular edge, from end to start, so we flip these as part of normalization
void VG::flip_doubly_reversed_edges(void) {
    for_each_edge([this](Edge* e) {
            if (e->from_start() && e->to_end()) {
                e->set_from_start(false);
                e->set_to_end(false);
                id_t f = e->to();
                id_t t = e->from();
                e->set_to(t);
                e->set_from(f);
            }
        });
    rebuild_edge_indexes();
}

// by definition, we can merge nodes that are a "simple component"
// without affecting the sequence or path space of the graph
// so we don't unchop nodes when they have mismatched path sets
void VG::unchop(void) {
    for (auto& comp : simple_multinode_components()) {
        concat_nodes(comp);
    }
    // rebuild path ranks, as these will be affected by mapping merging
    paths.compact_ranks();
}

void VG::normalize(void) {
    // convert edges that go from_start -> to_end to the equivalent "regular" edge
    flip_doubly_reversed_edges();
    // combine diced/chopped nodes (subpaths with no branching)
    unchop();
    // merge redundancy across multiple nodes into single nodes (requires flip_doubly_reversed_edges)
    simplify_siblings();
    // compact node ranks
    paths.compact_ranks();
    // there may now be some cut nodes that can be simplified
    unchop();
    // compact node ranks (again)
    paths.compact_ranks();
}

void VG::remove_non_path(void) {
    set<Edge*> path_edges;
    function<void(const Path&)> lambda = [this, &path_edges](const Path& path) {
        for (size_t i = 1; i < path.mapping_size(); ++i) {
            auto& m1 = path.mapping(i-1);
            auto& m2 = path.mapping(i);
            if (!adjacent_mappings(m1, m2)) continue; // the path is completely represented here
            auto s1 = NodeSide(m1.position().node_id(), (m1.position().is_reverse() ? false : true));
            auto s2 = NodeSide(m2.position().node_id(), (m2.position().is_reverse() ? true : false));
            // check that we always have an edge between the two nodes in the correct direction
            assert(has_edge(s1, s2));
            Edge* edge = get_edge(s1, s2);
            path_edges.insert(edge);
        }
    };
    paths.for_each(lambda);
    // now determine which edges aren't used
    set<Edge*> non_path_edges;
    for_each_edge([this, &path_edges, &non_path_edges](Edge* e) {
            if (!path_edges.count(e)) {
                non_path_edges.insert(e);
            }
        });
    // and destroy them
    for (auto* e : non_path_edges) {
        destroy_edge(e);
    }

    set<id_t> non_path_nodes;
    for_each_node([this, &non_path_nodes](Node* n) {
            if (!paths.has_node_mapping(n->id())) {
                non_path_nodes.insert(n->id());
            }
        });
    for (auto id : non_path_nodes) {
        destroy_node(id);
    }

    // re-compact ids if we have made changes to the graph
    if (!non_path_nodes.empty()) {
        sort();
        compact_ids();
    }
}

set<list<Node*>> VG::simple_multinode_components(void) {
    return simple_components(2);
}

// true if the mapping completely covers the node it maps to and is a perfect match
bool VG::mapping_is_total_match(const Mapping& m) {
    return mapping_is_simple_match(m)
        && mapping_from_length(m) == get_node(m.position().node_id())->sequence().size();
}

bool VG::nodes_are_perfect_path_neighbors(id_t id1, id_t id2) {
    // it is not possible for the nodes to be perfect neighbors if
    // they do not have exactly the same counts of paths
    if (paths.of_node(id1) != paths.of_node(id2)) return false;
    // now we know that the paths are identical in count and name between the two nodes

    // get the mappings for each node
    auto& m1 = paths.get_node_mapping(id1);
    auto& m2 = paths.get_node_mapping(id2);

    // verify that they are all perfect matches
    for (auto& p : m1) {
        for (auto* m : p.second) {
            if (!mapping_is_total_match(*m)) return false;
        }
    }
    for (auto& p : m2) {
        for (auto* m : p.second) {
            if (!mapping_is_total_match(*m)) return false;
        }
    }

    // it is still possible that we have the same path annotations
    // but the components of the paths we have are not contiguous across these nodes
    // to verify, we check that each mapping on the first immediately proceeds one on the second

    // order the mappings by rank so we can quickly check if everything is adjacent
    map<string, map<int, Mapping*>> r1, r2;
    for (auto& p : m1) {
        auto& name = p.first;
        auto& mp1 = p.second;
        auto& mp2 = m2[name];
        for (auto* m : mp1) r1[name][m->rank()] = m;
        for (auto* m : mp2) r2[name][m->rank()] = m;
    }
    // verify adjacency
    for (auto& p : r1) {
        auto& name = p.first;
        auto& ranked1 = p.second;
        map<int, Mapping*>& ranked2 = r2[name];
        for (auto& r : ranked1) {
            auto rank = r.first;
            auto& m = *r.second;
            auto f = ranked2.find(rank+(!m.position().is_reverse()? 1 : -1));
            if (f == ranked2.end()) return false;
            if (f->second->position().is_reverse() != m.position().is_reverse()) return false;
            ranked2.erase(f); // remove so we can verify that we have fully matched
        }
    }
    // verify that we fully matched the second node
    for (auto& p : r2) {
        if (!p.second.empty()) return false;
    }

    // we've passed all checks, so we have a node pair with mergable paths
    return true;
}

// the set of components that could be merged into single nodes without
// changing the path space of the graph
// respects stored paths
set<list<Node*>> VG::simple_components(int min_size) {

    // go around and establish groupings
    set<Node*> seen;
    set<list<Node*>> components;
    for_each_node([this, min_size, &components, &seen](Node* n) {
            if (seen.count(n)) return;
            seen.insert(n);
            // go left and right through each as far as we have only single edges connecting us
            // to nodes that have only single edges coming in or out
            // and these edges are "normal" in that they go from the tail to the head
            list<Node*> c;
            // go left
            {
                Node* l = n;
                auto sides = sides_to(NodeSide(l->id(), false));
                while (sides.size() == 1
                       && start_degree(l) == 1
                       && end_degree(get_node(sides.begin()->node)) == 1
                       && sides.begin()->is_end) {
                    id_t last_id = l->id();
                    l = get_node(sides.begin()->node);
                    seen.insert(l);
                    // avoid merging if it breaks stored paths
                    if (!nodes_are_perfect_path_neighbors(l->id(), last_id)) break;
                    sides = sides_to(NodeSide(l->id(), false));
                    c.push_front(l);
                }
            }
            // add the node (in the middle)
            c.push_back(n);
            // go right
            {
                Node* r = n;
                auto sides = sides_from(NodeSide(r->id(), true));
                while (sides.size() == 1
                       && end_degree(r) == 1
                       && start_degree(get_node(sides.begin()->node)) == 1
                       && !sides.begin()->is_end) {
                    id_t last_id = r->id();
                    seen.insert(r);
                    r = get_node(sides.begin()->node);
                    // avoid merging if it breaks stored paths
                    if (!nodes_are_perfect_path_neighbors(last_id, r->id())) break;
                    sides = sides_from(NodeSide(r->id(), true));
                    c.push_back(r);
                }
            }
            if (c.size() >= min_size) {
                components.insert(c);
            }
        });
    /*
    cerr << "components " << endl;
    for (auto& c : components) {
        for (auto x : c) {
            cerr << x->id() << " ";
        }
        cerr << endl;
    }
    */
    return components;
}

// merges right, so we take the rightmost rank as the new rank
map<string, map<int, Mapping>>
    VG::concat_mapping_groups(map<string, map<int, Mapping>>& r1,
                              map<string, map<int, Mapping>>& r2) {
    map<string, map<int, Mapping>> new_mappings;
    /*
    cerr << "merging mapping groups" << endl;
    cerr << "r1" << endl;
    for (auto& p : r1) {
        auto& name = p.first;
        auto& ranked1 = p.second;
        for (auto& r : ranked1) {
            cerr << pb2json(r.second) << endl;
        }
    }
    cerr << "r2" << endl;
    for (auto& p : r2) {
        auto& name = p.first;
        auto& ranked1 = p.second;
        for (auto& r : ranked1) {
            cerr << pb2json(r.second) << endl;
        }
    }
    cerr << "------------------" << endl;
    */
    // collect new mappings
    for (auto& p : r1) {
        auto& name = p.first;
        auto& ranked1 = p.second;
        map<int, Mapping>& ranked2 = r2[name];
        for (auto& r : ranked1) {
            auto rank = r.first;
            auto& m = r.second;
            auto f = ranked2.find(rank+(!m.position().is_reverse()? 1 : -1));
            //cerr << "seeking " << rank+(!m.position().is_reverse()? 1 : -1) << endl;
            assert(f != ranked2.end());
            auto& o = f->second;
            assert(m.position().is_reverse() == o.position().is_reverse());
            // make the new mapping for this pair of nodes
            Mapping n;
            if (!m.position().is_reverse()) {
                // in the forward orientation, we merge from left to right
                // and keep the right's rank
                n = concat_mappings(m, o);
                n.set_rank(o.rank());
            } else {
                // in the reverse orientation, we merge from left to right
                // but we keep the lower rank
                n = concat_mappings(o, m);
                n.set_rank(o.rank());
            }
            new_mappings[name][n.rank()] = n;
            ranked2.erase(f); // remove so we can verify that we have fully matched
        }
    }
    return new_mappings;
}

map<string, vector<Mapping>>
    VG::concat_mappings_for_nodes(const list<Node*>& nodes) {

    // determine the common paths that will apply to the new node
    // to do the ptahs right, we can only combine nodes if they also share all of their paths
    // and equal numbers of traversals
    set<map<string,int>> path_groups;
    for (auto n : nodes) {
        path_groups.insert(paths.of_node(n->id()));
    }

    if (path_groups.size() != 1) {
        cerr << "[VG::cat_nodes] error: cannot merge nodes with differing paths" << endl;
        exit(1); // we should be raising an error
    }

    auto ns = nodes; // to modify destructively
    auto np = nodes.front();
    ns.pop_front();
    // store the first base
    // we will use this to drive the merge
    auto base = paths.get_node_mapping_copies_by_rank(np->id());

    while (!ns.empty()) {
        // merge
        auto op = ns.front();
        ns.pop_front();
        // if this is our first batch, just keep them
        auto next = paths.get_node_mapping_copies_by_rank(op->id());
        // then merge the next in
        base = concat_mapping_groups(base, next);
    }

    // stores a merged mapping for each path traversal through the nodes we are merging
    map<string, vector<Mapping>> new_mappings;
    for (auto& p : base) {
        auto& name = p.first;
        for (auto& m : p.second) {
            new_mappings[name].push_back(m.second);
        }
    }

    return new_mappings;
}

Node* VG::concat_nodes(const list<Node*>& nodes) {

    // make the new mappings for the node
    map<string, vector<Mapping>> new_mappings = concat_mappings_for_nodes(nodes);

    // make a new node that concatenates the labels in the order they occur in the graph
    string seq;
    for (auto n : nodes) {
        seq += n->sequence();
    }
    auto node = create_node(seq);

    // remove the old mappings
    for (auto n : nodes) {
        set<Mapping*> to_remove;
        for (auto p : paths.get_node_mapping(n)) {
            for (auto* m : p.second) {
                to_remove.insert(m);
            }
        }
        for (auto m : to_remove) {
            paths.remove_mapping(m);
        }
    }

    // change the position of the new mappings to point to the new node
    // and store them in the path
    for (map<string, vector<Mapping>>::iterator nm = new_mappings.begin(); nm != new_mappings.end(); ++nm) {
        vector<Mapping>& ms = nm->second;
        for (vector<Mapping>::iterator m = ms.begin(); m != ms.end(); ++m) {
            m->mutable_position()->set_node_id(node->id());
            m->mutable_position()->set_offset(0); // uhhh
            if (m->position().is_reverse()) {
                paths.prepend_mapping(nm->first, *m);
            } else {
                paths.append_mapping(nm->first, *m);
            }
        }
    }

    // connect this node to the left and right connections of the set

    // do the left connections
    auto old_start = NodeSide(nodes.front()->id(), false);
    auto new_start = NodeSide(node->id(), false);
    // forward
    for (auto side : sides_to(old_start)) {
        create_edge(side, new_start);
    }
    // reverse
    for (auto side : sides_from(old_start)) {
        create_edge(new_start, side);
    }

    // do the right connections
    auto old_end = NodeSide(nodes.back()->id(), true);
    auto new_end = NodeSide(node->id(), true);
    // forward
    for (auto side : sides_from(old_end)) {
        create_edge(new_end, side);
    }
    // reverse
    for (auto side : sides_to(old_end)) {
        create_edge(side, new_end);
    }

    // remove the old nodes
    for (auto n : nodes) {
        destroy_node(n);
    }

    return node;
}

Node* VG::merge_nodes(const list<Node*>& nodes) {
    // make the new node (use the first one in the list)
    assert(!nodes.empty());
    Node* n = nodes.front();
    id_t nid = n->id();
    // create edges to the node
    for (auto& m : nodes) {
        if (m != n) { // skip first, which we're using
            //set<NodeSide> sides_of(NodeSide side);
            id_t id = m->id();
            for (auto& s : sides_to(NodeSide(id, false))) {
                create_edge(s, NodeSide(nid, false));
            }
            for (auto& s : sides_to(NodeSide(id, true))) {
                create_edge(s, NodeSide(nid, true));
            }
            for (auto& s : sides_from(NodeSide(id, false))) {
                create_edge(NodeSide(nid, false), s);
            }
            for (auto& s : sides_from(NodeSide(id, true))) {
                create_edge(NodeSide(nid, true), s);
            }
        }
    }
    // reassign mappings in paths to the new node
    hash_map<id_t, id_t> id_mapping;
    for (auto& m : nodes) {
        if (m != n) {
            id_mapping[m->id()] = nid;
        }
    }
    paths.swap_node_ids(id_mapping);
    // and erase the old nodes
    for (auto& m : nodes) {
        if (m != n) {
            destroy_node(m);
        }
    }
    // return the node we merged into
    return n;
}

id_t VG::total_length_of_nodes(void) {
    id_t length = 0;
    for (id_t i = 0; i < graph.node_size(); ++i) {
        Node* n = graph.mutable_node(i);
        length += n->sequence().size();
    }
    return length;
}

void VG::build_node_indexes(void) {
    for (id_t i = 0; i < graph.node_size(); ++i) {
        Node* n = graph.mutable_node(i);
        node_index[n] = i;
        node_by_id[n->id()] = n;
    }
}

void VG::build_edge_indexes(void) {
    for (id_t i = 0; i < graph.edge_size(); ++i) {
        Edge* e = graph.mutable_edge(i);
        edge_index[e] = i;
        set_edge(e);
    }
}

void VG::build_indexes(void) {
    build_node_indexes();
    build_edge_indexes();
}

void VG::clear_node_indexes(void) {
    node_index.clear();
    node_by_id.clear();
}

void VG::clear_node_indexes_no_resize(void) {
#ifdef USE_DENSE_HASH
    node_index.clear_no_resize();
    node_by_id.clear_no_resize();
#else
    clear_node_indexes();
#endif
}

void VG::clear_edge_indexes(void) {
    edge_by_sides.clear();
    edge_index.clear();
    edges_on_start.clear();
    edges_on_end.clear();
}

void VG::clear_edge_indexes_no_resize(void) {
#ifdef USE_DENSE_HASH
    edge_by_sides.clear_no_resize();
    edge_index.clear_no_resize();
    edges_on_start.clear_no_resize();
    edges_on_end.clear_no_resize();
#else
    clear_edge_indexes();
#endif
}

void VG::clear_indexes(void) {
    clear_node_indexes();
    clear_edge_indexes();
}

void VG::clear_indexes_no_resize(void) {
#ifdef USE_DENSE_HASH
    clear_node_indexes_no_resize();
    clear_edge_indexes_no_resize();
#else
    clear_indexes();
#endif
}

void VG::resize_indexes(void) {
    node_index.resize(graph.node_size());
    node_by_id.resize(graph.node_size());
    edge_by_sides.resize(graph.edge_size());
    edge_index.resize(graph.edge_size());
    edges_on_start.resize(graph.edge_size());
    edges_on_end.resize(graph.edge_size());
}

void VG::rebuild_indexes(void) {
    clear_indexes_no_resize();
    build_indexes();
    paths.rebuild_node_mapping();
}

void VG::rebuild_edge_indexes(void) {
    clear_edge_indexes_no_resize();
    build_edge_indexes();
}

bool VG::empty(void) {
    return graph.node_size() == 0 && graph.edge_size() == 0;
}

bool VG::has_node(Node* node) {
    return node && has_node(node->id());
}

bool VG::has_node(Node& node) {
    return has_node(node.id());
}

bool VG::has_node(id_t id) {
    return node_by_id.find(id) != node_by_id.end();
}

bool VG::has_edge(Edge* edge) {
    return edge && has_edge(*edge);
}

bool VG::has_edge(Edge& edge) {
    return edge_by_sides.find(NodeSide::pair_from_edge(edge)) != edge_by_sides.end();
}

bool VG::has_edge(const NodeSide& side1, const NodeSide& side2) {
    return edge_by_sides.find(minmax(side1, side2)) != edge_by_sides.end();
}

bool VG::has_edge(const pair<NodeSide, NodeSide>& sides) {
    return has_edge(sides.first, sides.second);
}

bool VG::has_inverting_edge(Node* n) {
    for (auto e : edges_of(n)) {
        if ((e->from_start() || e->to_end())
            && !(e->from_start() && e->to_end())) {
            return true;
        }
    }
    return false;
}

bool VG::has_inverting_edge_from(Node* n) {
    for (auto e : edges_of(n)) {
        if (e->from() == n->id()
            && (e->from_start() || e->to_end())
            && !(e->from_start() && e->to_end())) {
            return true;
        }
    }
    return false;
}

bool VG::has_inverting_edge_to(Node* n) {
    for (auto e : edges_of(n)) {
        if (e->to() == n->id()
            && (e->from_start() || e->to_end())
            && !(e->from_start() && e->to_end())) {
            return true;
        }
    }
    return false;
}

// remove duplicated nodes and edges that would occur if we merged the graphs
void VG::remove_duplicated_in(VG& g) {
    vector<Node*> nodes_to_destroy;
    for (id_t i = 0; i < graph.node_size(); ++i) {
        Node* n = graph.mutable_node(i);
        if (g.has_node(n)) {
            nodes_to_destroy.push_back(n);
        }
    }
    vector<Edge*> edges_to_destroy;
    for (id_t i = 0; i < graph.edge_size(); ++i) {
        Edge* e = graph.mutable_edge(i);
        if (g.has_edge(e)) {
            edges_to_destroy.push_back(e);
        }
    }
    for (vector<Node*>::iterator n = nodes_to_destroy.begin();
         n != nodes_to_destroy.end(); ++n) {
        g.destroy_node(g.get_node((*n)->id()));
    }
    for (vector<Edge*>::iterator e = edges_to_destroy.begin();
         e != edges_to_destroy.end(); ++e) {
        // Find and destroy the edge that does the same thing in g.
        g.destroy_edge(g.get_edge(NodeSide::pair_from_edge(*e)));
    }
}

void VG::merge_union(VG& g) {
    // remove duplicates, then merge
    remove_duplicated_in(g);
    if (g.graph.node_size() > 0) {
        merge(g.graph);
    }
}

void VG::merge(VG& g) {
    merge(g.graph);
}

// this merges without any validity checks
// this could be rather expensive if the graphs to merge are largely overlapping
void VG::merge(Graph& g) {
    graph.mutable_node()->MergeFrom(g.node());
    graph.mutable_edge()->MergeFrom(g.edge());
    rebuild_indexes();
}

// iterates over nodes and edges, adding them in when they don't already exist
void VG::extend(VG& g, bool warn_on_duplicates) {
    for (id_t i = 0; i < g.graph.node_size(); ++i) {
        Node* n = g.graph.mutable_node(i);
        if(n->id() == 0) {
            cerr << "[vg] warning: node ID 0 is not allowed. Skipping." << endl;
        } else if (!has_node(n)) {
            add_node(*n);
        } else if(warn_on_duplicates) {
            cerr << "[vg] warning: node ID " << n->id() << " appears multiple times. Skipping." << endl;
        }
    }
    for (id_t i = 0; i < g.graph.edge_size(); ++i) {
        Edge* e = g.graph.mutable_edge(i);
        if (!has_edge(e)) {
            add_edge(*e);
        } else if(warn_on_duplicates) {
            cerr << "[vg] warning: edge " << e->from() << (e->from_start() ? " start" : " end") << " <-> "
                 << e->to() << (e->to_end() ? " end" : " start") << " appears multiple times. Skipping." << endl;
        }
    }
    // todo breaks on cycles!
    paths.append(g.paths);
}

// TODO: unify with above. The only difference is what's done with the paths.
void VG::extend(Graph& graph, bool warn_on_duplicates) {
    for (id_t i = 0; i < graph.node_size(); ++i) {
        Node* n = graph.mutable_node(i);
        if(n->id() == 0) {
            cerr << "[vg] warning: node ID 0 is not allowed. Skipping." << endl;
        } else if (!has_node(n)) {
            add_node(*n);
        } else if(warn_on_duplicates) {
            cerr << "[vg] warning: node ID " << n->id() << " appears multiple times. Skipping." << endl;
        }
    }
    for (id_t i = 0; i < graph.edge_size(); ++i) {
        Edge* e = graph.mutable_edge(i);
        if (!has_edge(e)) {
            add_edge(*e);
        } else if(warn_on_duplicates) {
            cerr << "[vg] warning: edge " << e->from() << (e->from_start() ? " start" : " end") << " <-> "
                 << e->to() << (e->to_end() ? " end" : " start") << " appears multiple times. Skipping." << endl;
        }
    }
    paths.append(graph);
}

// extend this graph by g, connecting the tails of this graph to the heads of the other
// the ids of the second graph are modified for compact representation
void VG::append(VG& g) {

    // compact and increment the ids of g out of range of this graph
    //g.compact_ids();

    // assume we've already compacted the other, or that id compaction doesn't matter
    // just get out of the way
    g.increment_node_ids(max_node_id());

    // get the heads of the other graph, now that we've compacted the ids
    vector<Node*> heads = g.head_nodes();
    // The heads are guaranteed to be forward-oriented.
    vector<id_t> heads_ids;
    for (Node* n : heads) {
        heads_ids.push_back(n->id());
    }

    // get the current tails of this graph
    vector<Node*> tails = tail_nodes();
    // The tails are also guaranteed to be forward-oriented.
    vector<id_t> tails_ids;
    for (Node* n : tails) {
        tails_ids.push_back(n->id());
    }

    // add in the other graph
    // note that we don't use merge_union because we are ensured non-overlapping ids
    merge(g);

    /*
    cerr << "this graph size " << node_count() << " nodes " << edge_count() << " edges" << endl;
    cerr << "in append with " << heads.size() << " heads and " << tails.size() << " tails" << endl;
    */

    // now join the tails to heads
    for (id_t& tail : tails_ids) {
        for (id_t& head : heads_ids) {
            // Connect the tail to the head with a left to right edge.
            create_edge(tail, head);
        }
    }

    // wipe the ranks of the mappings, as these are destroyed in append
    // NB: append assumes that we are concatenating paths
    paths.clear_mapping_ranks();
    g.paths.clear_mapping_ranks();

    // and join paths that are embedded in the graph, where path names are the same
    paths.append(g.paths);
}

void VG::combine(VG& g) {
    // compact and increment the ids of g out of range of this graph
    //g.compact_ids();
    g.increment_node_ids(max_node_id());
    // now add it into the current graph, without connecting any nodes
    extend(g);
}

void VG::include(const Path& path) {
    for (size_t i = 0; i < path.mapping_size(); ++i) {
        if (!mapping_is_simple_match(path.mapping(i))) {
            cerr << "mapping " << pb2json(path.mapping(i)) << " cannot be included in the graph because it is not a simple match" << endl;
            //exit(1);
        }
    }
    paths.extend(path);
}

id_t VG::max_node_id(void) {
    id_t max_id = 0;
    for (int i = 0; i < graph.node_size(); ++i) {
        Node* n = graph.mutable_node(i);
        if (n->id() > max_id) {
            max_id = n->id();
        }
    }
    return max_id;
}

id_t VG::min_node_id(void) {
    id_t min_id = max_node_id();
    for (int i = 0; i < graph.node_size(); ++i) {
        Node* n = graph.mutable_node(i);
        if (n->id() < min_id) {
            min_id = n->id();
        }
    }
    return min_id;
}

void VG::compact_ids(void) {
    hash_map<id_t, id_t> new_id;
    id_t id = 1; // start at 1
    for_each_node([&id, &new_id](Node* n) {
            new_id[n->id()] = id++; });
//#pragma omp parallel for
    for_each_node([&new_id](Node* n) {
            n->set_id(new_id[n->id()]); });
//#pragma omp parallel for
    for_each_edge([&new_id](Edge* e) {
            e->set_from(new_id[e->from()]);
            e->set_to(new_id[e->to()]); });
    paths.swap_node_ids(new_id);
    rebuild_indexes();
}

void VG::increment_node_ids(id_t increment) {
    for_each_node_parallel([increment](Node* n) {
            n->set_id(n->id()+increment);
        });
    for_each_edge_parallel([increment](Edge* e) {
            e->set_from(e->from()+increment);
            e->set_to(e->to()+increment);
        });
    rebuild_indexes();
    paths.increment_node_ids(increment);
}

void VG::decrement_node_ids(id_t decrement) {
    increment_node_ids(-decrement);
}

void VG::swap_node_id(id_t node_id, id_t new_id) {
    swap_node_id(node_by_id[node_id], new_id);
}

void VG::swap_node_id(Node* node, id_t new_id) {

    int edge_n = edge_count();
    id_t old_id = node->id();
    node->set_id(new_id);
    node_by_id.erase(old_id);

    // we check if the old node exists, and bail out if we're not doing what we expect
    assert(node_by_id.find(new_id) == node_by_id.end());

    // otherwise move to a new id
    node_by_id[new_id] = node;

    // These are sets, so if we try to destroy and recreate the same edge from
    // both ends (i.e. if they both go to this node) we will only do it once.
    set<pair<NodeSide, NodeSide>> edges_to_destroy;
    set<pair<NodeSide, NodeSide>> edges_to_create;

    // Define a function that we will run on every edge this node is involved in
    auto fix_edge = [&](Edge* edge) {

        // Destroy that edge
        edges_to_destroy.emplace(NodeSide(edge->from(), !edge->from_start()), NodeSide(edge->to(), edge->to_end()));

        // Make a new edge with our new ID as from or to (or both), depending on which it was before.
        // TODO: Is there a cleaner way to do this?
        if(edge->from() == old_id) {
            if(edge->to() == old_id) {
                edges_to_create.emplace(NodeSide(new_id, !edge->from_start()), NodeSide(new_id, edge->to_end()));
            } else {
                edges_to_create.emplace(NodeSide(new_id, !edge->from_start()), NodeSide(edge->to(), edge->to_end()));
            }
        } else {
            edges_to_create.emplace(NodeSide(edge->from(), !edge->from_start()), NodeSide(new_id, edge->to_end()));
        }

    };

    for(pair<id_t, bool>& other : edges_start(old_id)) {
        // Get the actual Edge
        // We're at a start, so we go to the end of the other node normally, and the start if the other node is backward
        Edge* edge = edge_by_sides[minmax(NodeSide(old_id, false), NodeSide(other.first, !other.second))];

        // Plan to fix up its IDs.
        fix_edge(edge);
    }

    for(pair<id_t, bool>& other : edges_end(old_id)) {
        // Get the actual Edge
        // We're at an end, so we go to the start of the other node normally, and the end if the other node is backward
        Edge* edge = edge_by_sides[minmax(NodeSide(old_id, true), NodeSide(other.first, other.second))];

        // Plan to fix up its IDs.
        fix_edge(edge);
    }

    assert(edges_to_destroy.size() == edges_to_create.size());

    for (auto& e : edges_to_destroy) {
        // Destroy the edge (only one can exist between any two nodes)
        destroy_edge(e.first, e.second);
    }

    for (auto& e : edges_to_create) {
        // Make an edge with the appropriate start and end flags
        create_edge(e.first, e.second);
    }

    assert(edge_n == edge_count());

    // we maintain a valid graph
    // this an expensive check but should work (for testing only)
    //assert(is_valid());

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

void VG::vcf_records_to_alleles(vector<vcflib::Variant>& records,
                                map<long, set<vcflib::VariantAllele> >& altp,
                                int start_pos,
                                int stop_pos,
                                int max_node_size,
                                bool flat_input_vcf) {

    create_progress("parsing variants", records.size());

    for (int i = 0; i < records.size(); ++i) {
        vcflib::Variant& var = records.at(i);
        // decompose to alts
        map<string, vector<vcflib::VariantAllele> > alternates
            = (flat_input_vcf ? var.flatAlternates() : var.parsedAlternates());
        for (auto& alleles : alternates) {
            for (auto& allele : alleles.second) {
                altp[allele.position].insert(allele);
                if (i % 10000 == 0) {
                    update_progress(altp.size());
                }
            }
        }
    }
    destroy_progress();
}

void VG::slice_alleles(map<long, set<vcflib::VariantAllele> >& altp,
                       int start_pos,
                       int stop_pos,
                       int max_node_size) {

    auto enforce_node_size_limit =
        [this, max_node_size, &altp]
        (int curr_pos, int& last_pos) {
        int last_ref_size = curr_pos - last_pos;
        update_progress(last_pos);
        if (max_node_size && last_ref_size > max_node_size) {
            int div = 2;
            while (last_ref_size/div > max_node_size) {
                ++div;
            }
            int segment_size = last_ref_size/div;
            int i = 0;
            while (last_pos + i < curr_pos) {
                altp[last_pos+i];  // empty cut
                i += segment_size;
                update_progress(last_pos + i);
            }
        }
    };

    if (max_node_size > 0) {
        create_progress("enforcing node size limit ", (altp.empty()? 0 : altp.rbegin()->first));
        // break apart big nodes
        int last_pos = start_pos;
        for (auto& position : altp) {
            auto& alleles = position.second;
            enforce_node_size_limit(position.first, last_pos);
            for (auto& allele : alleles) {
                // cut the last reference sequence into bite-sized pieces
                last_pos = max(position.first + allele.ref.size(), (long unsigned int) last_pos);
            }
        }
        enforce_node_size_limit(stop_pos, last_pos);
        destroy_progress();
    }

}

void VG::dice_nodes(int max_node_size) {
    // We're going to chop up everything, so clear out the path ranks.
    paths.clear_mapping_ranks();

    if (max_node_size) {
        vector<Node*> nodes; nodes.reserve(size());
        for_each_node(
            [this, &nodes](Node* n) {
                nodes.push_back(n);
            });
        auto lambda =
            [this, max_node_size](Node* n) {
            int node_size = n->sequence().size();
            if (node_size > max_node_size) {
                Node* l = NULL;
                Node* r = NULL;
                int div = 2;
                while (node_size/div > max_node_size) {
                    ++div;
                }
                int segment_size = node_size/div;
                int i = 0;
                while (n->sequence().size() > segment_size) {
                    // For each division between the div new pieces,
                    // break off a piece of the determined size at the left
                    // TODO: for very large nodes, this is n * (n/max_node_size) = O(n^2).
                    // We could do n log n if we split nodes in their middles.
                    divide_node(n, segment_size, l, r);
                    // Keep dividing the right node
                    n = r;
                }
            }
        };
        for (int i = 0; i < nodes.size(); ++i) {
            lambda(nodes[i]);
        }
    }

    // Set the ranks again
    paths.rebuild_mapping_aux();
    paths.compact_ranks();
}

void VG::from_alleles(const map<long, set<vcflib::VariantAllele> >& altp,
                      string& seq,
                      string& name) {

    //init();
    this->name = name;

    int tid = omp_get_thread_num();
#ifdef debug
#pragma omp critical (cerr)
    {
        cerr << tid << ": in from_alleles" << endl;
        cerr << tid << ": with " << altp.size() << " vars" << endl;
        cerr << tid << ": and " << seq.size() << "bp" << endl;
        cerr << seq << endl;
    }
#endif


    // maintains the path of the seq in the graph
    map<long, id_t> seq_node_ids;
    // track the last nodes so that we can connect everything
    // completely when variants occur in succession
    map<long, set<Node*> > nodes_by_end_position;
    map<long, set<Node*> > nodes_by_start_position;


    Node* seq_node = create_node(seq);
    seq_node_ids[0] = seq_node->id();

    for (auto& va : altp) {

        const set<vcflib::VariantAllele>& alleles = va.second;

        // if alleles are empty, we just cut at this point
        if (alleles.empty()) {
            Node* l = NULL; Node* r = NULL;
            divide_path(seq_node_ids, va.first, l, r);
        }

        for (auto allele : alleles) {

            // skip ref-matching alleles; these are not informative
            if (allele.ref == allele.alt) {
                continue;
            }

#ifdef debug
#pragma omp critical (cerr)
            cerr << tid << ": " << allele << endl;
#endif

            // 0/1 based conversion happens in offset
            long allele_start_pos = allele.position;
            long allele_end_pos = allele_start_pos + allele.ref.size();
            // for ordering, set insertion start position at +1
            // otherwise insertions at the same position will loop infinitely
            //if (allele_start_pos == allele_end_pos) allele_end_pos++;

            if (allele_start_pos == 0) {
                // ensures that we can handle variation at first position
                // (important when aligning)
                Node* root = create_node("");
                seq_node_ids[-1] = root->id();
                nodes_by_start_position[-1].insert(root);
                nodes_by_end_position[0].insert(root);
            }

            Node* left_seq_node = NULL;
            Node* middle_seq_node = NULL;
            Node* right_seq_node = NULL;

            // make one cut at the ref-path relative start of the allele
            divide_path(seq_node_ids,
                        allele_start_pos,
                        left_seq_node,
                        right_seq_node);

            // if the ref portion of the allele is not empty, then we need to make another cut
            if (!allele.ref.empty()) {
                divide_path(seq_node_ids,
                            allele_end_pos,
                            middle_seq_node,
                            right_seq_node);
            }

            Node* alt_node = NULL;
            // create a new alt node and connect the pieces from before
            if (!allele.alt.empty() && !allele.ref.empty()) {
                //cerr << "both alt and ref have sequence" << endl;

                alt_node = create_node(allele.alt);
                create_edge(left_seq_node, alt_node);
                create_edge(alt_node, right_seq_node);

                nodes_by_end_position[allele_end_pos].insert(alt_node);
                nodes_by_end_position[allele_end_pos].insert(middle_seq_node);
                //nodes_by_end_position[allele_start_pos].insert(left_seq_node);
                nodes_by_start_position[allele_start_pos].insert(alt_node);
                nodes_by_start_position[allele_start_pos].insert(middle_seq_node);

            } else if (!allele.alt.empty()) { // insertion

                alt_node = create_node(allele.alt);
                create_edge(left_seq_node, alt_node);
                create_edge(alt_node, right_seq_node);
                nodes_by_end_position[allele_end_pos].insert(alt_node);
                nodes_by_end_position[allele_end_pos].insert(left_seq_node);
                nodes_by_start_position[allele_start_pos].insert(alt_node);

            } else {// otherwise, we have a deletion

                create_edge(left_seq_node, right_seq_node);
                nodes_by_end_position[allele_end_pos].insert(left_seq_node);
                nodes_by_start_position[allele_start_pos].insert(left_seq_node);

            }

#ifdef debug
#pragma omp critical (cerr)
            {
                if (left_seq_node) cerr << tid << ": left_ref " << left_seq_node->id()
                                        << " " << left_seq_node->sequence() << endl;
                if (middle_seq_node) cerr << tid << ": middle_ref " << middle_seq_node->id()
                                          << " " << middle_seq_node->sequence() << endl;
                if (alt_node) cerr << tid << ": alt_node " << alt_node->id()
                                   << " " << alt_node->sequence() << endl;
                if (right_seq_node) cerr << tid << ": right_ref " << right_seq_node->id()
                                         << " " << right_seq_node->sequence() << endl;
            }
#endif

            if (allele_end_pos == seq.size()) {
                // ensures that we can handle variation at last position (important when aligning)
                Node* end = create_node("");
                seq_node_ids[allele_end_pos] = end->id();
                // for consistency, this should be handled below in the start/end connections
                if (alt_node) {
                    create_edge(alt_node, end);
                }
                if (middle_seq_node) {
                    create_edge(middle_seq_node, end);
                }
            }

            //print_edges();
            /*
            if (!is_valid()) {
                cerr << "graph is invalid after variant " << *a << endl;
                std::ofstream out("fail.vg");
                serialize_to_ostream(out);
                out.close();
                exit(1);
            }
            */

        }

        map<long, set<Node*> >::iterator ep
            = nodes_by_end_position.find(va.first);
        map<long, set<Node*> >::iterator sp
            = nodes_by_start_position.find(va.first);
        if (ep != nodes_by_end_position.end()
            && sp != nodes_by_start_position.end()) {
            set<Node*>& previous_nodes = ep->second;
            set<Node*>& current_nodes = sp->second;
            for (set<Node*>::iterator n = previous_nodes.begin();
                 n != previous_nodes.end(); ++n) {
                for (set<Node*>::iterator m = current_nodes.begin();
                     m != current_nodes.end(); ++m) {
                    if (node_index.find(*n) != node_index.end()
                        && node_index.find(*m) != node_index.end()
                        && !(previous_nodes.count(*n) && current_nodes.count(*n)
                             && previous_nodes.count(*m) && current_nodes.count(*m))
                        ) {
                        /*
                        cerr << "connecting previous "
                             << (*n)->id() << " @end=" << ep->first << " to current "
                             << (*m)->id() << " @start=" << sp->first << endl;
                        */
                        create_edge(*n, *m);
                    }
                }
            }
        }

        // clean up previous
        while (!nodes_by_end_position.empty() && nodes_by_end_position.begin()->first < va.first) {
            nodes_by_end_position.erase(nodes_by_end_position.begin()->first);
        }

        while (!nodes_by_start_position.empty() && nodes_by_start_position.begin()->first < va.first) {
            nodes_by_start_position.erase(nodes_by_start_position.begin()->first);
        }

    }

    // serialize path
    for (auto& p : seq_node_ids) {
        paths.append_mapping(name, p.second);
    }
    // and set the mapping edits
    force_path_match();

    sort();
    compact_ids();

}

void VG::from_gfa(istream& in, bool showp) {
    // c++... split...
    // for line in stdin
    string line;
    auto too_many_fields = [&line]() {
        cerr << "[vg] error: too many fields in line " << endl << line << endl;
        exit(1);
    };

		GFAKluge gg;
		gg.parse_gfa_file(in);

		map<string, sequence_elem, custom_key> name_to_seq = gg.get_name_to_seq();
		map<std::string, vector<link_elem> > seq_to_link = gg.get_seq_to_link();
    map<string, vector<path_elem> > seq_to_paths = gg.get_seq_to_paths();
		map<string, sequence_elem>::iterator it;
		for (it = name_to_seq.begin(); it != name_to_seq.end(); it++){
			auto source_id = atol((it->second).name.c_str());
			//Make us some nodes
			Node n;
			n.set_sequence((it->second).sequence);
			n.set_id(source_id);
			add_node(n);
			// Now some edges. Since they're placed in this map
			// by their from_node, it's no big deal to just iterate
			// over them.
			for (link_elem l : seq_to_link[(it->second).name]){
				auto sink_id = atol(l.sink_name.c_str());
				Edge e;
				e.set_from(source_id);
				e.set_to(sink_id);
				e.set_from_start(!l.source_orientation_forward);
				e.set_to_end(!l.sink_orientation_forward);
				add_edge(e);
			}
      for (path_elem p: seq_to_paths[(it->second).name]){
        paths.append_mapping(p.name, source_id, p.rank ,p.is_reverse);
      }

		}

}

void VG::print_edges(void) {
    for (int i = 0; i < graph.edge_size(); ++i) {
        Edge* e = graph.mutable_edge(i);
        id_t f = e->from();
        id_t t = e->to();
        cerr << f << "->" << t << " ";
    }
    cerr << endl;
}

void VG::create_progress(const string& message, long count) {
    if (show_progress) {
        progress_message = message;
        create_progress(count);
    }
}

void VG::create_progress(long count) {
    if (show_progress) {
        progress_message.resize(30, ' ');
        progress_count = count;
        last_progress = 0;
        progress = new ProgressBar(progress_count, progress_message.c_str());
        progress->Progressed(0);
    }
}

void VG::update_progress(long i) {
    if (show_progress && progress) {
        if (i <= progress_count
            && (long double) (i - last_progress) / (long double) progress_count >= 0.001
            || i == progress_count) {
#pragma omp critical (progress)
            {
                progress->Progressed(i);
                last_progress = i;
            }
        }
    }
}

void VG::destroy_progress(void) {
    if (show_progress && progress) {
        update_progress(progress_count);
        cerr << endl;
        progress_message = "";
        progress_count = 0;
        delete progress;
        progress = NULL;
    }
}

VG::VG(vcflib::VariantCallFile& variantCallFile,
       FastaReference& reference,
       string& target_region,
       bool target_is_chrom,
       int vars_per_region,
       int max_node_size,
       bool flat_input_vcf,
       bool showprog) {

    init();

    omp_set_dynamic(1); // use dynamic scheduling

    show_progress = showprog;

    map<string, VG*> refseq_graph;

    vector<string> targets;
    if (!target_region.empty()) {
        targets.push_back(target_region);
    } else {
        for (vector<string>::iterator r = reference.index->sequenceNames.begin();
             r != reference.index->sequenceNames.end(); ++r) {
            targets.push_back(*r);
        }
    }

    // to scale up, we have to avoid big string memcpys
    // this could be accomplished by some deep surgery on the construction routines
    // however, that could be a silly thing to do,
    // because why break something that's conceptually clear
    // and anyway, we want to break the works into chunks
    //
    // there is a development that could be important
    // our chunk size isn't going to reach into the range where we'll have issues (>several megs)
    // so we'll run this for regions of moderate size, scaling up in the case that we run into a big deletion
    //
    //

    for (vector<string>::iterator t = targets.begin(); t != targets.end(); ++t) {

        //string& seq_name = *t;
        string seq_name;
        string target = *t;
        int start_pos = 0, stop_pos = 0;
        // nasty hack for handling single regions
        if (!target_is_chrom) {
            parse_region(target,
                         seq_name,
                         start_pos,
                         stop_pos);
            if (stop_pos > 0) {
                if (variantCallFile.is_open()) {
                    variantCallFile.setRegion(seq_name, start_pos, stop_pos);
                }
            } else {
                if (variantCallFile.is_open()) {
                    variantCallFile.setRegion(seq_name);
                }
                stop_pos = reference.sequenceLength(seq_name);
            }
        } else {
            // the user said the target is just a sequence name
            // and is unsafe to parse as it may contain ':' or '-'
            // for example "gi|568815592:29791752-29792749"
            if (variantCallFile.is_open()) {
                variantCallFile.setRegion(target);
            }
            stop_pos = reference.sequenceLength(target);
            seq_name = target;
        }
        vcflib::Variant var(variantCallFile);

        vector<vcflib::Variant>* region = NULL;

        // convert from 1-based input to 0-based internal format
        // and handle the case where we are already doing the whole chromosome
        id_t start = start_pos ? start_pos - 1 : 0;
        id_t end = start;

        create_progress("loading variants for " + target, stop_pos-start_pos);
        // get records
        vector<vcflib::Variant> records;
        int i = 0;
        while (variantCallFile.is_open() && variantCallFile.getNextVariant(var)) {
            // this ... maybe we should remove it as for when we have calls against N
            bool isDNA = allATGC(var.ref);
            for (vector<string>::iterator a = var.alt.begin(); a != var.alt.end(); ++a) {
                if (!allATGC(*a)) isDNA = false;
            }
            // only work with DNA sequences
            if (isDNA) {
                var.position -= 1; // convert to 0-based
                records.push_back(var);
            }
            if (++i % 1000 == 0) update_progress(var.position-start_pos);
        }
        destroy_progress();

        map<long,set<vcflib::VariantAllele> > alleles;
        // decompose records int alleles with offsets against our target sequence
        vcf_records_to_alleles(records, alleles, start_pos, stop_pos, max_node_size, flat_input_vcf);
        records.clear(); // clean up

        // enforce a maximum node size
        // by dividing nodes that are > than the max into the smallest number of
        // even pieces that would be smaller than the max
        slice_alleles(alleles, start_pos, stop_pos, max_node_size);

        // store our construction plans
        deque<Plan*> construction;
        // so we can check which graphs we can safely append
        set<VG*> graph_completed;
        // we add and remove from graph_completed, so track count for logging
        int graphs_completed = 0;
        int final_completed = -1; // hm
        // the construction queue
        list<VG*> graphq;
        int graphq_size = 0; // for efficiency
        // ^^^^ (we need to insert/remove things in the middle of the list,
        // but we also need to be able to quickly determine its size)
        // for tracking progress through the chromosome
        map<VG*, unsigned long> graph_end;

        create_progress("planning construction", stop_pos-start_pos);
        // break into chunks
        int chunk_start = start;
        bool invariant_graph = alleles.empty();
        while (invariant_graph || !alleles.empty()) {
            invariant_graph = false;
            auto* new_alleles = new map<long, set<vcflib::VariantAllele> >;
            // our start position is the "offset" we should subtract from the alleles
            // for correct construction
            //chunk_start = (!chunk_start ? 0 : alleles.begin()->first);
            int chunk_end = chunk_start;
            bool clean_end = true;
            for (int i = 0; (i < vars_per_region || !clean_end) && !alleles.empty(); ++i) {
                auto pos = alleles.begin()->first - chunk_start;
                chunk_end = max(chunk_end, (int)alleles.begin()->first);
                auto& pos_alleles = alleles.begin()->second;
                // apply offset when adding to the new alleles
                auto& curr_pos = (*new_alleles)[pos];
                for (auto& allele : pos_alleles) {
                    auto new_allele = allele;
                    int ref_end = new_allele.ref.size() + new_allele.position;
                    // look through the alleles to see if there is a longer chunk
                    if (ref_end > chunk_end) {
                        chunk_end = ref_end;
                    }
                    new_allele.position = pos;
                    curr_pos.insert(new_allele);
                }
                alleles.erase(alleles.begin());
                // TODO here we need to see if we are neighboring another variant
                // and if we are, keep constructing
                if (alleles.begin()->first <= chunk_end) {
                    clean_end = false;
                } else {
                    clean_end = true;
                }
            }
            // record end position, use target end in the case that we are at the end
            if (alleles.empty()) chunk_end = stop_pos;

            // we set the head graph to be this one, so we aren't obligated to copy the result into this object
            // make a construction plan
            Plan* plan = new Plan(graphq.empty() && targets.size() == 1 ? this : new VG,
                                  new_alleles,
                                  reference.getSubSequence(seq_name,
                                                           chunk_start,
                                                           chunk_end - chunk_start),
                                  seq_name);
            chunk_start = chunk_end;
#pragma omp critical (graphq)
            {
                graphq.push_back(plan->graph);
                construction.push_back(plan);
                if (show_progress) graph_end[plan->graph] = chunk_end;
                update_progress(chunk_end);
            }
        }
#ifdef debug
        cerr << omp_get_thread_num() << ": graphq size " << graphq.size() << endl;
#endif
        graphq_size = graphq.size();
        destroy_progress();

        // this system is not entirely general
        // there will be a problem when the regions of overlapping deletions become too large
        // then the inter-dependence of each region will make parallel construction in this way difficult
        // because the chunks will get too large

        // use this function to merge graphs both during and after the construction iteration
        auto merge_first_two_completed_graphs =
            [this, start_pos, &graph_completed, &graphq, &graphq_size, &graph_end, &final_completed](void) {
            // find the first two consecutive graphs which are completed
            VG* first = NULL;
            VG* second = NULL;
//#pragma omp critical (cerr)
//            cerr << omp_get_thread_num() << ": merging" << endl;
#pragma omp critical (graphq)
            {
                auto itp = graphq.begin(); // previous
                auto itn = itp; if (itp != graphq.end()) ++itn; // next
                // scan the graphq to find consecutive entries that are both completed
                while (itp != itn // there is > 1 entry
                       && itn != graphq.end() // we aren't yet at the end
                       && !(graph_completed.count(*itp) // the two we're looking at aren't completed
                            && graph_completed.count(*itn))) {
                    ++itp; ++itn;
                }

                if (itn != graphq.end()) {
                    // we have two consecutive graphs to merge!
                    first = *itp;
                    second = *itn;
                    // unset graph completed for both
                    graph_completed.erase(first);
                    graph_completed.erase(second);
                    graphq.erase(itn);
                    --graphq_size;
                }
            }

            if (first && second) {
                // combine graphs
                first->append(*second);
#pragma omp critical (graphq)
                {
                    if (final_completed != -1) update_progress(final_completed++);
                    graph_completed.insert(first);
                    graph_end.erase(second);
                }
                delete second;
            }
        };

        create_progress("constructing graph", construction.size());

        // (in parallel) construct each component of the graph
#pragma omp parallel for
        for (int i = 0; i < construction.size(); ++i) {

            int tid = omp_get_thread_num();
            Plan* plan = construction.at(i);
#ifdef debug
#pragma omp critical (cerr)
            cerr << tid << ": " << "constructing graph " << plan->graph << " over "
                 << plan->alleles->size() << " variants in " <<plan->seq.size() << "bp "
                 << plan->name << endl;
#endif

            plan->graph->from_alleles(*plan->alleles,
                                      plan->seq,
                                      plan->name);
#pragma omp critical (graphq)
            {
                update_progress(++graphs_completed);
                graph_completed.insert(plan->graph);
#ifdef debug
#pragma omp critical (cerr)
                cerr << tid << ": " << "constructed graph " << plan->graph << endl;
#endif
            }
            // clean up
            delete plan;

            // concatenate chunks of the result graph together
            merge_first_two_completed_graphs();

        }
        destroy_progress();

        // merge remaining graphs
        final_completed = 0;
        create_progress("merging remaining graphs", graphq.size());
#pragma omp parallel
        {
            bool more_to_merge = true;
            while (more_to_merge) {
                merge_first_two_completed_graphs();
                usleep(10);
#pragma omp critical (graphq)
                more_to_merge = graphq_size > 1;
            }
        }
        destroy_progress();

        // parallel end
        // finalize target

        // our target graph should be the only entry in the graphq
        assert(graphq.size() == 1);
        VG* target_graph = graphq.front();

        // store it in our results
        refseq_graph[target] = target_graph;

        create_progress("joining graphs", target_graph->size());
        // clean up "null" nodes that are used for maintaining structure between temporary subgraphs
        target_graph->remove_null_nodes_forwarding_edges();
        destroy_progress();

        // then use topological sorting and re-compression of the id space to make sure that
        create_progress("topologically sorting", target_graph->size());
        target_graph->sort();
        destroy_progress();

        create_progress("compacting ids", target_graph->size());
        // we get identical graphs no matter what the region size is
        target_graph->compact_ids();
        destroy_progress();

    }

    // hack for efficiency when constructing over a single chromosome
    if (refseq_graph.size() == 1) {
        // *this = *refseq_graph[targets.front()];
        // we have already done this because the first graph in the queue is this
    } else {
        // where we have multiple targets
        for (vector<string>::iterator t = targets.begin(); t != targets.end(); ++t) {
            // merge the variants into one graph
            VG& g = *refseq_graph[*t];
            combine(g);
        }
    }
    // rebuild the mapping ranks now that we've combined everything
    paths.clear_mapping_ranks();
    paths.rebuild_mapping_aux();
}

void VG::sort(void) {
    if (size() <= 1) return;
    // Topologically sort, which orders and orients all the nodes.
    deque<NodeTraversal> sorted_nodes;
    topological_sort(sorted_nodes);
    deque<NodeTraversal>::iterator n = sorted_nodes.begin();
    int i = 0;
    for ( ; i < graph.node_size() && n != sorted_nodes.end();
          ++i, ++n) {
        // Put the nodes in the order we got
        swap_nodes(graph.mutable_node(i), (*n).node);
    }
}

// Tarjan's strongly connected components algorithm
// transcribed from pseudocode
// https://en.wikipedia.org/wiki/Tarjan%27s_strongly_connected_components_algorithm
set<set<id_t> > VG::strongly_connected_components(void) {
    set<set<id_t> > components;
    deque<id_t> stack;
    set<id_t> on_stack;
    int64_t index = 0;
    map<id_t, int64_t> node_idx;
    map<id_t, int64_t> node_lowlink;

    //function<void(Node*)> strong_connect;
    function<void(id_t)> strong_connect = [&](id_t v) {
        node_idx[v] = index;
        node_lowlink[v] = index;
        ++index;
        stack.push_back(v);
        on_stack.insert(v);
        // Consider successors of v
        for (auto side : sides_from(v)) {
            id_t w = side.node;
            if (node_idx.find(w) == node_idx.end()) {
                // Successor w has not yet been visited; recurse on it
                strong_connect(w);
                node_lowlink[v] = min(node_lowlink[v], node_lowlink[w]);
            } else if (on_stack.count(w)) {
                // Successor w is in stack S and hence in the current SCC
                node_lowlink[v] = min(node_lowlink[v], node_idx[w]);
            }
        }
        // If v is a root node, pop the stack and generate an SCC
        if (node_lowlink[v] == node_idx[v]) {
            // start a new strongly connected component
            set<id_t> component;
            id_t w;
            do {
                // get w from the stack
                w = stack.back();
                stack.pop_back();
                on_stack.erase(w);
                // add w to current strongly connected component
                component.insert(w);
            } while (w != v);
            // output the current strongly connected component
            components.insert(component);
        }
    };

    // now run strong_connect across the graph
    for_each_node([&](Node* n) {
            id_t v = n->id();
            if (node_idx.find(v) == node_idx.end()) {
                strong_connect(v);
            }
        });

    return components;
}

vector<Edge> VG::break_cycles(void) {
    // ensure we are sorted
    sort();
    // remove any edge whose from has a higher index than its to
    vector<Edge*> to_remove;
    for_each_edge([&](Edge* e) {
            // if we cycle to this node or one before in the sort
            if (node_index[get_node(e->from())] >= node_index[get_node(e->to())]) {
                to_remove.push_back(e);
            }
        });
    vector<Edge> removed;
    for(Edge* edge : to_remove) {
        //cerr << "removing " << pb2json(*edge) << endl;
        removed.push_back(*edge);
        destroy_edge(edge);
    }
    return removed;
}

bool VG::is_acyclic(void) {
    return multinode_strongly_connected_components().empty();
}

set<set<id_t> > VG::multinode_strongly_connected_components(void) {
    set<set<id_t> > components;
    for (auto& c : strongly_connected_components()) {
        if (c.size() > 1) {
            components.insert(c);
        }
    }
    return components;
}

// keeping all components would be redundant, as every node is a self-component
void VG::keep_multinode_strongly_connected_components(void) {
    set<id_t> keep;
    for (auto& c : multinode_strongly_connected_components()) {
        for (auto& id : c) {
            keep.insert(id);
        }
    }
    set<Node*> remove;
    for_each_node([&](Node* n) {
            if (!keep.count(n->id())) {
                remove.insert(n);
            }
        });
    for (auto n : remove) {
        destroy_node(n);
    }
    remove_orphan_edges();
}

size_t VG::size(void) {
    return graph.node_size();
}

size_t VG::length(void) {
    size_t l = 0;
    for_each_node([&l](Node* n) { l+=n->sequence().size(); });
    return l;
}

void VG::swap_nodes(Node* a, Node* b) {
    int aidx = node_index[a];
    int bidx = node_index[b];
    graph.mutable_node()->SwapElements(aidx, bidx);
    node_index[a] = bidx;
    node_index[b] = aidx;
}

Edge* VG::create_edge(NodeTraversal left, NodeTraversal right) {
    // Connect to the start of the left node if it is backward, and the end of the right node if it is backward.
    return create_edge(left.node->id(), right.node->id(), left.backward, right.backward);
}

Edge* VG::create_edge(NodeSide side1, NodeSide side2) {
    // Connect to node 1 (from start if the first side isn't an end) to node 2 (to end if the second side is an end)
    return create_edge(side1.node, side2.node, !side1.is_end, side2.is_end);
}

Edge* VG::create_edge(Node* from, Node* to, bool from_start, bool to_end) {
    return create_edge(from->id(), to->id(), from_start, to_end);
}

Edge* VG::create_edge(id_t from, id_t to, bool from_start, bool to_end) {
    //cerr << "creating edge " << from << "->" << to << endl;
    // ensure the edge (or another between the same sides) does not already exist
    Edge* edge = get_edge(NodeSide(from, !from_start), NodeSide(to, to_end));
    if (edge) {
        // The edge we want to make exists.
        return edge;
    }
    // if not, create it
    edge = graph.add_edge();
    edge->set_from(from);
    edge->set_to(to);
    // Only set the backwardness fields if they are true.
    if(from_start) edge->set_from_start(from_start);
    if(to_end) edge->set_to_end(to_end);
    set_edge(edge);
    edge_index[edge] = graph.edge_size()-1;
    //cerr << "created edge " << edge->from() << "->" << edge->to() << endl;
    return edge;
}

Edge* VG::get_edge(const NodeSide& side1, const NodeSide& side2) {
    auto e = edge_by_sides.find(minmax(side1, side2));
    if (e != edge_by_sides.end()) {
        return e->second;
    } else {
        return NULL;
    }
}

Edge* VG::get_edge(const pair<NodeSide, NodeSide>& sides) {
    return get_edge(sides.first, sides.second);
}

Edge* VG::get_edge(const NodeTraversal& left, const NodeTraversal& right) {
    // We went from the right side of left to the left side of right.
    // We used the end of left if if isn't backward, and we used the end of right if it is.
    return get_edge(NodeSide(left.node->id(), !left.backward),
                    NodeSide(right.node->id(), right.backward));
}

void VG::set_edge(Edge* edge) {
    // Note: there must not be an edge between these sides of these nodes already.
    if (!has_edge(edge)) {
        // Note that we might add edges to nonexistent nodes (like in VG::node_context()). That's just fine.

        // Add the edge to the index by node side (edges_on_start, edges_on_end, and edge_by_sides)
        index_edge_by_node_sides(edge);
    }
}

void VG::for_each_edge_parallel(function<void(Edge*)> lambda) {
    create_progress(graph.edge_size());
    id_t completed = 0;
#pragma omp parallel for shared(completed)
    for (id_t i = 0; i < graph.edge_size(); ++i) {
        lambda(graph.mutable_edge(i));
        if (progress && completed++ % 1000 == 0) {
#pragma omp critical (progress_bar)
            update_progress(completed);
        }
    }
    destroy_progress();
}

void VG::for_each_edge(function<void(Edge*)> lambda) {
    for (id_t i = 0; i < graph.edge_size(); ++i) {
        lambda(graph.mutable_edge(i));
    }
}

void VG::destroy_edge(const NodeSide& side1, const NodeSide& side2) {
    destroy_edge(get_edge(side1, side2));
}

void VG::destroy_edge(const pair<NodeSide, NodeSide>& sides) {
    destroy_edge(sides.first, sides.second);
}


void VG::destroy_edge(Edge* edge) {
    //cerr << "destroying edge " << edge->from() << "->" << edge->to() << endl;

    // noop on NULL pointer or non-existent edge
    if (!has_edge(edge)) { return; }

    // first remove the edge from the edge-on-node-side indexes.
    unindex_edge_by_node_sides(edge);

    // get the last edge index (lei) and this edge index (tei)
    int lei = graph.edge_size()-1;
    int tei = edge_index[edge];

    // erase this edge from the index by node IDs.
    // we'll fix up below
    edge_index.erase(edge);

    // Why do we check that lei != tei?
    //
    // It seems, after an inordinate amount of testing and probing,
    // that if we call erase twice on the same entry, we'll end up corrupting the hash_map
    //
    // So, if the element is already at the end of the table,
    // take a little break and just remove the last edge in graph

    // if we need to move the element to the last position in the array...
    if (lei != tei) {

        // get a pointer to the last element
        Edge* last = graph.mutable_edge(lei);

        // erase from our index
        edge_index.erase(last);

        // swap
        graph.mutable_edge()->SwapElements(tei, lei);

        // point to new position
        Edge* nlast = graph.mutable_edge(tei);

        // insert the new edge index position
        edge_index[nlast] = tei;

        // and fix edge indexes for moved edge object
        set_edge(nlast);

    }

    // drop the last position, erasing the node
    // manually delete to free memory (RemoveLast does not free)
    Edge* last_edge = graph.mutable_edge()->ReleaseLast();
    delete last_edge;

    //if (!is_valid()) { cerr << "graph ain't valid" << endl; }

}

void VG::unindex_edge_by_node_sides(const NodeSide& side1, const NodeSide& side2) {
    unindex_edge_by_node_sides(get_edge(side1, side2));
}

void VG::unindex_edge_by_node_sides(Edge* edge) {
    // noop on NULL pointer or non-existent edge
    if (!has_edge(edge)) return;
    //if (!is_valid()) { cerr << "graph ain't valid" << endl; }
    // erase from indexes

    auto edge_pair = NodeSide::pair_from_edge(edge);

    //cerr << "erasing from indexes" << endl;

    //cerr << "Unindexing edge " << edge_pair.first << "<-> " << edge_pair.second << endl;

    // Remove from the edge by node side pair index
    edge_by_sides.erase(edge_pair);

    // Does this edge involve a change of relative orientation?
    bool relative_orientation = edge->from_start() != edge->to_end();

    // Un-index its from node, depending on whether it's attached to the start
    // or end.
    if(edge->from_start()) {
        // The edge is on the start of the from node, so remove it from the
        // start of the from node, with the correct relative orientation for the
        // to node.
        std::pair<id_t, bool> to_remove {edge->to(), relative_orientation};
        swap_remove(edges_start(edge->from()), to_remove);
        // removing the sub-indexes if they are now empty
        // we must do this to maintain a valid structure
        if (edges_on_start[edge->from()].empty()) edges_on_start.erase(edge->from());

        //cerr << "Removed " << edge->from() << "-start to " << edge->to() << " orientation " << relative_orientation << endl;
    } else {
        // The edge is on the end of the from node, do remove it form the end of the from node.
        std::pair<id_t, bool> to_remove {edge->to(), relative_orientation};
        swap_remove(edges_end(edge->from()), to_remove);
        if (edges_on_end[edge->from()].empty()) edges_on_end.erase(edge->from());

        //cerr << "Removed " << edge->from() << "-end to " << edge->to() << " orientation " << relative_orientation << endl;
    }

    if(edge->from() != edge->to() || edge->from_start() == edge->to_end()) {
        // Same for the to node, if we aren't just on the same node and side as with the from node.
        if(edge->to_end()) {
            std::pair<id_t, bool> to_remove {edge->from(), relative_orientation};
            swap_remove(edges_end(edge->to()), to_remove);
            if (edges_on_end[edge->to()].empty()) edges_on_end.erase(edge->to());

            //cerr << "Removed " << edge->to() << "-end to " << edge->from() << " orientation " << relative_orientation << endl;
        } else {
            std::pair<id_t, bool> to_remove {edge->from(), relative_orientation};
            swap_remove(edges_start(edge->to()), to_remove);
            if (edges_on_start[edge->to()].empty()) edges_on_start.erase(edge->to());

            //cerr << "Removed " << edge->to() << "-start to " << edge->from() << " orientation "
            //     << relative_orientation << endl;
        }
    }
}

void VG::index_edge_by_node_sides(Edge* edge) {

    // Generate sides, order them, and index the edge by them.
    edge_by_sides[NodeSide::pair_from_edge(edge)] = edge;

    // Index on ends appropriately depending on from_start and to_end.
    bool relative_orientation = edge->from_start() != edge->to_end();

    if(edge->from_start()) {
        edges_on_start[edge->from()].emplace_back(edge->to(), relative_orientation);
    } else {
        edges_on_end[edge->from()].emplace_back(edge->to(), relative_orientation);
    }

    if(edge->from() != edge->to() || edge->from_start() == edge->to_end()) {
        // Only index the other end of the edge if the edge isn't a self-loop on a single side.
        if(edge->to_end()) {
            edges_on_end[edge->to()].emplace_back(edge->from(), relative_orientation);
        } else {
            edges_on_start[edge->to()].emplace_back(edge->from(), relative_orientation);
        }
    }
}

Node* VG::get_node(id_t id) {
    hash_map<id_t, Node*>::iterator n = node_by_id.find(id);
    if (n != node_by_id.end()) {
        return n->second;
    } else {
        serialize_to_file("wtf.vg");
        throw runtime_error("No node " + to_string(id) + " in graph");
    }
}

Node* VG::create_node(const string& seq, id_t id) {
    // create the node
    Node* node = graph.add_node();
    node->set_sequence(seq);
    // ensure we properly update the current_id that's used to generate new ids
    // unless we have a specified id
    if (id == 0) {
        if (current_id == 1) current_id = max_node_id()+1;
        node->set_id(current_id++);
    } else {
        node->set_id(id);
    }
    // copy it into the graphnn
    // and drop into our id index
    node_by_id[node->id()] = node;
    node_index[node] = graph.node_size()-1;
    //if (!is_valid()) cerr << "graph invalid" << endl;
    return node;
}

void VG::for_each_node_parallel(function<void(Node*)> lambda) {
    create_progress(graph.node_size());
    id_t completed = 0;
    #pragma omp parallel for schedule(dynamic,1) shared(completed)
    for (id_t i = 0; i < graph.node_size(); ++i) {
        lambda(graph.mutable_node(i));
        if (progress && completed++ % 1000 == 0) {
            #pragma omp critical (progress_bar)
            update_progress(completed);
        }
    }
    destroy_progress();
}

void VG::for_each_node(function<void(Node*)> lambda) {
    for (id_t i = 0; i < graph.node_size(); ++i) {
        lambda(graph.mutable_node(i));
    }
}

void VG::for_each_connected_node(Node* node, function<void(Node*)> lambda) {
    // We keep track of nodes to visit.
    set<Node*> to_visit {node};
    // We mark all the nodes we have visited.
    set<Node*> visited;

    while(!to_visit.empty()) {
        // Grab some node to visit
        Node* visiting = *(to_visit.begin());
        to_visit.erase(to_visit.begin());

        // Visit it
        lambda(visiting);
        visited.insert(visiting);

        // Look at all its edges
        vector<Edge*> edges;
        edges_of_node(visiting, edges);
        for(Edge* edge : edges) {
            if(edge->from() != visiting->id() && !visited.count(get_node(edge->from()))) {
                // We found a new node on the from of the edge
                to_visit.insert(get_node(edge->from()));
            } else if(edge->to() != visiting->id() && !visited.count(get_node(edge->to()))) {
                // We found a new node on the to of the edge
                to_visit.insert(get_node(edge->to()));
            }
        }
    }
}

// a graph composed of this node and the edges that can be uniquely assigned to it
void VG::nonoverlapping_node_context_without_paths(Node* node, VG& g) {
    // add the node
    g.add_node(*node);

    auto grab_edge = [&](Edge* e) {
        // What node owns the edge?
        id_t owner_id = min(e->from(), e->to());
        if(node->id() == owner_id || !has_node(owner_id)) {
            // Either we are the owner, or the owner isn't in the graph to get serialized.
            g.add_edge(*e);
        }
    };

    // Go through all its edges
    vector<pair<id_t, bool>>& start = edges_start(node->id());
    for (auto& e : start) {
        grab_edge(get_edge(NodeSide::pair_from_start_edge(node->id(), e)));
    }
    vector<pair<id_t, bool>>& end = edges_end(node->id());
    for (auto& e : end) {
        grab_edge(get_edge(NodeSide::pair_from_end_edge(node->id(), e)));
    }
    // paths must be added externally
}

void VG::destroy_node(id_t id) {
    destroy_node(get_node(id));
}

void VG::destroy_node(Node* node) {
    //if (!is_valid()) cerr << "graph is invalid before destroy_node" << endl;
    //cerr << "destroying node " << node->id() << " degrees " << start_degree(node) << ", " << end_degree(node) << endl;
    // noop on NULL/nonexistent node
    if (!has_node(node)) { return; }
    // remove edges associated with node
    set<pair<NodeSide, NodeSide>> edges_to_destroy;

    for(auto& other_end : edges_start(node)) {
        // Destroy all the edges on its start
        edges_to_destroy.insert(NodeSide::pair_from_start_edge(node->id(), other_end));
    }

    for(auto& other_end : edges_end(node)) {
        // Destroy all the edges on its end
        edges_to_destroy.insert(NodeSide::pair_from_end_edge(node->id(), other_end));
    }

    for (auto& e : edges_to_destroy) {
        //cerr << "Destroying edge " << e.first << ", " << e.second << endl;
        destroy_edge(e.first, e.second);
        //cerr << "Edge destroyed" << endl;
    }

    // assert cleanup
    edges_on_start.erase(node->id());
    edges_on_end.erase(node->id());

    //assert(start_degree(node) == 0);
    //assert(end_degree(node) == 0);

    // swap node with the last in nodes
    // call RemoveLast() to drop the node
    int lni = graph.node_size()-1;
    int tni = node_index[node];

    if (lni != tni) {
        // swap this node with the last one
        Node* last = graph.mutable_node(lni);
        graph.mutable_node()->SwapElements(tni, lni);
        Node* nlast = graph.mutable_node(tni);

        // and fix up the indexes
        node_by_id[last->id()] = nlast;
        node_index.erase(last);
        node_index[nlast] = tni;
    }

    // remove this node (which is now the last one) and remove references from the indexes
    node_by_id.erase(node->id());
    node_index.erase(node);
    // manually delete to free memory (RemoveLast does not free)
    Node* last_node = graph.mutable_node()->ReleaseLast();
    delete last_node;
    //if (!is_valid()) { cerr << "graph is invalid after destroy_node" << endl; exit(1); }
}

void VG::remove_null_nodes(void) {
    vector<id_t> to_remove;
    for (int i = 0; i < graph.node_size(); ++i) {
        Node* node = graph.mutable_node(i);
        if (node->sequence().size() == 0) {
            to_remove.push_back(node->id());
        }
    }
    for (vector<id_t>::iterator n = to_remove.begin(); n != to_remove.end(); ++n) {
        destroy_node(*n);
    }
}

void VG::remove_null_nodes_forwarding_edges(void) {
    vector<Node*> to_remove;
    int i = 0;
    create_progress(graph.node_size()*2);
    for (i = 0; i < graph.node_size(); ++i) {
        Node* node = graph.mutable_node(i);
        if (node->sequence().size() == 0) {
            to_remove.push_back(node);
        }
        update_progress(i);
    }
    for (vector<Node*>::iterator n = to_remove.begin(); n != to_remove.end(); ++n, ++i) {
        remove_node_forwarding_edges(*n);
        update_progress(i);
    }
}

void VG::remove_node_forwarding_edges(Node* node) {
    // Grab all the nodes attached to our start, with true if the edge goes to their start
    vector<pair<id_t, bool>>& start = edges_start(node);
    // Grab all the nodes attached to our end, with true if the edge goes to their end
    vector<pair<id_t, bool>>& end = edges_end(node);

    // We instantiate the whole cross product first to avoid working on
    // references to the contents of containers we are modifying. This holds the
    // (node ID, relative orientation) pairs above.
    set<pair<pair<id_t, bool>, pair<id_t, bool>>> edges_to_create;

    // Make edges for the cross product of our start and end edges, making sure
    // to maintain relative orientation.
    for(auto& start_pair : start) {
        for(auto& end_pair : end) {
            // We already have the flags for from_start and to_end for the new edge.
            edges_to_create.emplace(start_pair, end_pair);
        }
    }

    for (auto& e : edges_to_create) {
        // make each edge we want to add
        create_edge(e.first.first, e.second.first, e.first.second, e.second.second);
    }

    // remove the node from paths
    if (paths.has_node_mapping(node)) {
        // We need to copy the set here because we're going to be throwing
        // things out of it while iterating over it.
        map<string, set<Mapping*>> node_mappings(paths.get_node_mapping(node));
        for (auto& p : node_mappings) {
            for (auto& m : p.second) {
                paths.remove_mapping(m);
            }
        }
    }
    // delete the actual node
    destroy_node(node);
}

void VG::remove_orphan_edges(void) {
    set<pair<NodeSide, NodeSide>> edges;
    for_each_edge([this,&edges](Edge* edge) {
            if (!has_node(edge->from())
                || !has_node(edge->to())) {
                edges.insert(NodeSide::pair_from_edge(edge));
            }
        });
    for (auto edge : edges) {
        destroy_edge(edge);
    }
}

void VG::keep_paths(set<string>& path_names, set<string>& kept_names) {
    // Strategy:

    // For each node in our graph
    // If it has node mappings, look at them
    // For each that occurs along a path we want to include
    // Mark that path as kept
    // Mark the node as kept
    // Peek in the path left and right from that mapping
    // Mark the edges traversed as kept
    // Mark all nodes and edges not marked for keeping as for removal
    // Remove all nodes and edges marked for removal
    // Drop all paths from the graph not requested to be kept

    // Previous code here would only keep edges between nodes with adjacent IDs
    // in the set of selected nodes, which is in general incorrect.

    vector<Node*> nodes_to_keep;
    vector<Node*> nodes_to_remove;
    // Holds node sides in smaller-to-larger order, so we can search for edges in it.
    set<pair<NodeSide, NodeSide>> edges_to_keep;

    for_each_node([&](Node* node) {
        // If we don't see anything about this node in our paths, throw it out.
        bool to_keep = false;

        if(paths.has_node_mapping(node)) {
            // This node appears on some paths. Look for appearances on paths we like.
            for(auto& appearance : paths.get_node_mapping(node)) {
                if(path_names.count(appearance.first)) {
                    // We found an appearance of this node on a path we are keeping. It comes with a Mapping*.

                    // Mark the path and node as kept
                    kept_names.insert(appearance.first);
                    to_keep = true;

                    // Walk left along the path and keep the edge(s) we traverse.
                    for (auto* m : appearance.second) {
                        Mapping* left_neighbor = paths.traverse_left(m);

                        if(left_neighbor != nullptr) {
                            // We aren't the first thing in the path, we want to keep the edge to our left.
                            // It may not exist, but we can still ask to keep it.

                            // What other node do we go to?
                            id_t neighbor_id = left_neighbor->position().node_id();

                            if(has_node(neighbor_id)) {
                                // Keep the edge if we actually have the other end.
                                // We know the other end is on this path and will be
                                // kept.

                                // We attach to the end of the previous node if it isn't
                                // backward along the path, and the end of this node if
                                // it is backward along the path.
                                edges_to_keep.insert(minmax(NodeSide(neighbor_id, !left_neighbor->position().is_reverse()),
                                                            NodeSide(node->id(), m->position().is_reverse())));
                            }
                        }
                    }

                    // We skip walking right along the path, because if the node
                    // we would find is in the graph, we're going to walk left
                    // from it eventually and find the same edge.
                }
            }
        }

        // Actually decide to keep or delete the node. Keep it if we ever saw it on a path.
        if (to_keep) {
            nodes_to_keep.push_back(node);
        } else {
            nodes_to_remove.push_back(node);
        }

    });

    // Mark any edges we don't keep for destruction
    set<pair<NodeSide, NodeSide>> edges_to_destroy;
    for_each_edge([this, &edges_to_keep, &edges_to_destroy](Edge* edge) {
            auto ep = NodeSide::pair_from_edge(edge);
            if (!edges_to_keep.count(ep)) {
                edges_to_destroy.insert(ep);
            }
        });

    // Destroy all the edges and nodes we don't want
    for (auto edge : edges_to_destroy) {
        destroy_edge(edge.first, edge.second);
    }
    for (auto node : nodes_to_remove) {
        destroy_node(node);
    }

    // Throw out all the paths data for paths we don't want to keep.
    paths.keep_paths(path_names);
}

void VG::keep_path(string& path_name) {
    set<string> s,k; s.insert(path_name);
    keep_paths(s, k);
}

// divide a node into two pieces at the given offset
void VG::divide_node(Node* node, int pos, Node*& left, Node*& right) {

#ifdef debug_divide
#pragma omp critical (cerr)
    cerr << "dividing node " << node->id() << " at " << pos << endl;
#endif

    if (pos < 0 || pos > node->sequence().size()) {
#pragma omp critical (cerr)
        {
            cerr << omp_get_thread_num() << ": cannot divide node " << node->id() << ":" << node->sequence()
                 << " -- position (" << pos << ") is less than 0 or greater than sequence length ("
                 << node->sequence().size() << ")" << endl;
            exit(1);
        }
    }

#ifdef debug_divide
#pragma omp critical (cerr)
    cerr << omp_get_thread_num() << ": in divide_node " << pos << " of " << node->sequence().size() << endl;
#endif

    // make our left node
    left = create_node(node->sequence().substr(0,pos));

    // make our right node
    right = create_node(node->sequence().substr(pos));

    // Create edges between the left node (optionally its start) and the right node (optionally its end)
    set<pair<pair<id_t, bool>, pair<id_t, bool>>> edges_to_create;

    // replace the connections to the node's start
    for(auto e : edges_start(node)) {
        // We have to check for self loops, as these will be clobbered by the division of the node
        if (e.first == node->id()) {
            // if it's a reversed edge, it would be from the start of this node
            if (e.second) {
                // so set it to the left node
                e.first = left->id();
            } else {
                // otherwise, it's from the end, so set it to the right side
                e.first = right->id();
            }
        }
        // Make an edge to the left node's start from wherever this edge went.
        edges_to_create.emplace(make_pair(e.first, e.second), make_pair(left->id(), false));
    }

    // replace the connections to the node's end
    for(auto e : edges_end(node)) {
        // We have to check for self loops, as these will be clobbered by the division of the node
        if (e.first == node->id()) {
            // if it's a reversed edge, it would be to the end of this node
            if (e.second) {
                // so set it to the right node
                e.first = right->id();
            } else {
                // otherwise, it's to the start, so set it to the left side
                e.first = left->id();
            }
        }
        // Make an edge from the right node's end to wherever this edge went.
        edges_to_create.emplace(make_pair(right->id(), false), make_pair(e.first, e.second));
    }

    // create the edges here as otherwise we will invalidate the iterators
    for (auto& e : edges_to_create) {
        // Swizzle the from_start and to_end bits to the right place.
        create_edge(e.first.first, e.second.first, e.first.second, e.second.second);
    }

    // connect left to right. This edge always goes from end to start.
    create_edge(left, right);

    //cerr << "Dividing into " << left->id() << " and " << right->id() << endl;

    // divide paths
    // note that we can't do this (yet) for non-exact matching paths
    if (paths.has_node_mapping(node)) {
        auto& node_path_mapping = paths.get_node_mapping(node);
        // apply to left and right
        vector<Mapping*> to_divide;
        for (auto& pm : node_path_mapping) {
            string path_name = pm.first;
            for (auto* m : pm.second) {
                to_divide.push_back(m);
            }
        }
        for (auto m : to_divide) {
            // we have to divide the mapping

#ifdef debug

#pragma omp critical (cerr)
            cerr << omp_get_thread_num() << ": dividing mapping " << pb2json(*m) << endl;
#endif

            string path_name = paths.mapping_path_name(m);
            // TODO: this only preserves perfect match paths.
            // TODO: warn if that precondition is violated
            //divide_invariant_mapping(*m, l, r, pos, node, left, right);
            auto n = cut_mapping(*m, pos);
            Mapping l=n.first, r=n.second;
            l.mutable_position()->set_node_id(left->id());
            l.mutable_position()->set_offset(0);
            r.mutable_position()->set_node_id(right->id());
            r.mutable_position()->set_offset(0);

            //divide_mapping
            // with the mapping divided, insert the pieces where the old one was
            auto mpit = paths.remove_mapping(m);
            if(m->position().is_reverse()) {
                // insert left then right in the path, snce we're going through
                // this node backward (insert puts *before* the iterator)
                mpit = paths.insert_mapping(mpit, path_name, l);
                mpit = paths.insert_mapping(mpit, path_name, r);
            } else {
                // insert right then left (insert puts *before* the iterator)
                mpit = paths.insert_mapping(mpit, path_name, r);
                mpit = paths.insert_mapping(mpit, path_name, l);
            }



#ifdef debug_divide
#pragma omp critical (cerr)
            cerr << omp_get_thread_num() << ": produced mappings " << pb2json(l) << " and " << pb2json(r) << endl;
#endif

        }
    }

    destroy_node(node);

}

// for dividing a path of nodes with an underlying coordinate system
void VG::divide_path(map<long, id_t>& path, long pos, Node*& left, Node*& right) {

    map<long, id_t>::iterator target = path.upper_bound(pos);
    --target; // we should now be pointing to the target ref node

    long node_pos = target->first;
    Node* old = get_node(target->second);

    // nothing to do
    if (node_pos == pos) {
        map<long, id_t>::iterator n = target; --n;
        left = get_node(n->second);
        right = get_node(target->second);
    } else {
        // divide the target node at our pos
        int diff = pos - node_pos;
        divide_node(old, diff, left, right);
        // left
        path[node_pos] = left->id();
        // right
        path[pos] = right->id();
    }
}

set<NodeTraversal> VG::travs_of(NodeTraversal node) {
    auto tos = travs_to(node);
    auto froms = travs_from(node);
    set<NodeTraversal> ofs;
    std::set_union(tos.begin(), tos.end(),
                   froms.begin(), froms.end(),
                   std::inserter(ofs, ofs.begin()));
    return ofs;
}

// traversals before this node on the same strand
set<NodeTraversal> VG::travs_to(NodeTraversal node) {
    set<NodeTraversal> tos;
    vector<NodeTraversal> tov;
    nodes_prev(node, tov);
    for (auto& t : tov) tos.insert(t);
    return tos;
}

// traversals after this node on the same strand
set<NodeTraversal> VG::travs_from(NodeTraversal node) {
    set<NodeTraversal> froms;
    vector<NodeTraversal> fromv;
    nodes_next(node, fromv);
    for (auto& t : fromv) froms.insert(t);
    return froms;
}

void VG::nodes_prev(NodeTraversal node, vector<NodeTraversal>& nodes) {
    // Get the node IDs that attach to the left of this node, and whether we are
    // attached relatively forward (false) or backward (true)
    vector<pair<id_t, bool>>& left_nodes = node.backward ? edges_end(node.node) : edges_start(node.node);
    for (auto& prev : left_nodes) {
        // Make a NodeTraversal that is an oriented description of the node attached to our relative left.
        // If we're backward, and it's in the same relative orientation as us, it needs to be backward too.
        nodes.emplace_back(get_node(prev.first), prev.second != node.backward);
    }
}

vector<NodeTraversal> VG::nodes_prev(NodeTraversal n) {
    vector<NodeTraversal> nodes;
    nodes_prev(n, nodes);
    return nodes;
}

void VG::nodes_next(NodeTraversal node, vector<NodeTraversal>& nodes) {
    // Get the node IDs that attach to the right of this node, and whether we
    // are attached relatively forward (false) or backward (true)
    vector<pair<id_t, bool>>& right_nodes = node.backward ? edges_start(node.node) : edges_end(node.node);
    for (auto& next : right_nodes) {
        // Make a NodeTraversal that is an oriented description of the node attached to our relative right.
        // If we're backward, and it's in the same relative orientation as us, it needs to be backward too.
        nodes.emplace_back(get_node(next.first), next.second != node.backward);
    }
}

vector<NodeTraversal> VG::nodes_next(NodeTraversal n) {
    vector<NodeTraversal> nodes;
    nodes_next(n, nodes);
    return nodes;
}

int VG::node_count_prev(NodeTraversal n) {
    vector<NodeTraversal> nodes;
    nodes_prev(n, nodes);
    return nodes.size();
}

int VG::node_count_next(NodeTraversal n) {
    vector<NodeTraversal> nodes;
    nodes_next(n, nodes);
    return nodes.size();
}

void VG::prev_kpaths_from_node(NodeTraversal node, int length, int edge_max, bool edge_bounding,
                               list<NodeTraversal> postfix, set<list<NodeTraversal> >& paths,
                               function<void(NodeTraversal)>& maxed_nodes) {

   // length gives the number of bases *off the end of the current node* to
   // look, and does not get the length of the current node chargged against it.
   // We keep the last node that length at all reaches into.

#ifdef debug
    cerr << "Looking left from " << node << " out to length " << length <<
        " with remaining edges " << edge_max << " on top of:" << endl;
    for(auto x : postfix) {
        cerr << "\t" << x << endl;
    }
#endif

    if(edge_bounding && edge_max < 0) {
        // (recursive) caller must check edge_max and call maxed_nodes for this node
        cerr << "Called prev_kpaths_from_node with negative edges left." << endl;
        exit(1);
    }

    // start at node
    // do a leftward DFS up to length limit to establish paths from the left of the node
    postfix.push_front(node);
    // Get all the nodes left of this one
    vector<NodeTraversal> prev_nodes;
    nodes_prev(node, prev_nodes);

    // If we can't find any valid extensions, we have to just emit up to here as
    // a path.
    bool valid_extensions = false;

    if(length > 0) {
        // We're allowed to look off our end

        for (NodeTraversal& prev : prev_nodes) {
#ifdef debug
            cerr << "Consider prev node " << prev << endl;
#endif

            if(edge_bounding && edge_max - (left_degree(node) > 1) < 0) {
                // We won't have been able to get anything out of this next node
                maxed_nodes(prev);

#ifdef debug
                cerr << "Out of edge-crossing range" << endl;
#endif

            } else {
#ifdef debug
                cerr << "Recursing..." << endl;
#endif
                prev_kpaths_from_node(prev,
                                      length - prev.node->sequence().size(),
                                      // Charge 1 against edge_max for every time we pass up alternative edges
                                      edge_max - (left_degree(node) > 1),
                                      // but only if we are using edge bounding
                                      edge_bounding,
                                      postfix, paths, maxed_nodes);

                // We found a valid extension of this node
                valid_extensions = true;

            }


        }
    } else {
#ifdef debug
        cerr << "No length remaining." << endl;
# endif
    }

    if(!valid_extensions) {
        // We didn't find an extension to do, either because we ran out of edge
        // crossings, or because our length will run out somewhere in this node.
        // Create a path for this node.
        paths.insert(postfix);
#ifdef debug
        cerr << "Reported path:" << endl;
        for(auto x : postfix) {
            cerr << "\t" << x << endl;
        }
#endif
    }
}

void VG::next_kpaths_from_node(NodeTraversal node, int length, int edge_max, bool edge_bounding,
                               list<NodeTraversal> prefix, set<list<NodeTraversal> >& paths,
                               function<void(NodeTraversal)>& maxed_nodes) {

    if(edge_bounding && edge_max < 0) {
        // (recursive) caller must check edge_max and call maxed_nodes for this node
        cerr << "Called next_kpaths_from_node with negative edges left." << endl;
        exit(1);
    }

    // start at node
    // do a leftward DFS up to length limit to establish paths from the left of the node
    prefix.push_back(node);
    vector<NodeTraversal> next_nodes;
    nodes_next(node, next_nodes);

    // If we can't find any valid extensions, we have to just emit up to here as
    // a path.
    bool valid_extensions = false;

    if(length > 0) {
        // We're allowed to look off our end

        for (NodeTraversal& next : next_nodes) {

            if(edge_bounding && edge_max - (right_degree(node) > 1) < 0) {
                // We won't have been able to get anything out of this next node
                maxed_nodes(next);
            } else {
                next_kpaths_from_node(next,
                                  length - next.node->sequence().size(),
                                  // Charge 1 against edge_max for every time we pass up alternative edges
                                  edge_max - (right_degree(node) > 1),
                                  // but only if we are using edge bounding
                                  edge_bounding,
                                  prefix, paths, maxed_nodes);

                // We found a valid extension of this node
                valid_extensions = true;

            }
        }
    }

    if(!valid_extensions) {
        // We didn't find an extension to do, either because we ran out of edge
        // crossings, or because our length will run out somewhere in this node.
        // Create a path for this node.
        paths.insert(prefix);
    }
}

// iterate over the kpaths in the graph, doing something

void VG::for_each_kpath(int k, int edge_max,
                        function<void(NodeTraversal)> prev_maxed,
                        function<void(NodeTraversal)> next_maxed,
                        function<void(list<NodeTraversal>::iterator,list<NodeTraversal>&)> lambda) {
    auto by_node = [k, edge_max, &lambda, &prev_maxed, &next_maxed, this](Node* node) {
        for_each_kpath_of_node(node, k, edge_max, prev_maxed, next_maxed, lambda);
    };
    for_each_node(by_node);
}

void VG::for_each_kpath(int k, int edge_max,
                        function<void(NodeTraversal)> prev_maxed,
                        function<void(NodeTraversal)> next_maxed,
                        function<void(size_t,Path&)> lambda) {
    auto by_node = [k, edge_max, &lambda, &prev_maxed, &next_maxed, this](Node* node) {
        for_each_kpath_of_node(node, k, edge_max, prev_maxed, next_maxed, lambda);
    };
    for_each_node(by_node);
}

// parallel versions of above
// this isn't by default because the lambda may have side effects
// that need to be guarded explicitly

void VG::for_each_kpath_parallel(int k, int edge_max,
                                 function<void(NodeTraversal)> prev_maxed, function<void(NodeTraversal)> next_maxed,
                                 function<void(list<NodeTraversal>::iterator,list<NodeTraversal>&)> lambda) {
    auto by_node = [k, edge_max, &prev_maxed, &next_maxed, &lambda, this](Node* node) {
        for_each_kpath_of_node(node, k, edge_max, prev_maxed, next_maxed, lambda);
    };
    for_each_node_parallel(by_node);
}

void VG::for_each_kpath_parallel(int k, int edge_max,
                                 function<void(NodeTraversal)> prev_maxed,
                                 function<void(NodeTraversal)> next_maxed,
                                 function<void(size_t,Path&)> lambda) {
    auto by_node = [k, edge_max, &prev_maxed, &next_maxed, &lambda, this](Node* node) {
        for_each_kpath_of_node(node, k, edge_max, prev_maxed, next_maxed, lambda);
    };
    for_each_node_parallel(by_node);
}

// per-node kpaths

void VG::for_each_kpath_of_node(Node* n, int k, int edge_max,
                                function<void(NodeTraversal)> prev_maxed,
                                function<void(NodeTraversal)> next_maxed,
                                function<void(size_t,Path&)> lambda) {
    auto apply_to_path = [&lambda, this](list<NodeTraversal>::iterator n, list<NodeTraversal>& p) {
        Path path = create_path(p);

        // Find the index of the node with a scan. We already did O(n) work to make the path.
        // TODO: is there a more efficient way?
        size_t index = 0;
        while(n != p.begin()) {
            --n;
            ++index;
        }

        lambda(index, path);
    };
    for_each_kpath_of_node(n, k, edge_max, prev_maxed, next_maxed, apply_to_path);
}

void VG::for_each_kpath_of_node(Node* node, int k, int edge_max,
                                function<void(NodeTraversal)> prev_maxed,
                                function<void(NodeTraversal)> next_maxed,
                                function<void(list<NodeTraversal>::iterator,list<NodeTraversal>&)> lambda) {
    // get left, then right
    set<list<NodeTraversal> > prev_paths;
    set<list<NodeTraversal> > next_paths;
    list<NodeTraversal> empty_list;

    prev_kpaths_from_node(NodeTraversal(node), k, edge_max, (edge_max != 0),
                          empty_list, prev_paths, prev_maxed);
    next_kpaths_from_node(NodeTraversal(node), k, edge_max, (edge_max != 0),
                          empty_list, next_paths, next_maxed);
    // now take the cross and give to the callback
    for (set<list<NodeTraversal> >::iterator p = prev_paths.begin(); p != prev_paths.end(); ++p) {
        for (set<list<NodeTraversal> >::iterator n = next_paths.begin(); n != next_paths.end(); ++n) {
            list<NodeTraversal> path = *p;

            // Find the iterator to this node in the list that will become the
            // path. We know it's the last thing in the prev kpath we made our path from.
            list<NodeTraversal>::iterator this_node = path.end(); this_node--;

            list<NodeTraversal>::const_iterator m = n->begin(); ++m; // skips current node, which is included in *p in the correct orientation
            while (m != n->end()) {
                // Add on all the other nodes from the next kpath.
                path.push_back(*m);
                ++m;
            }

            lambda(this_node, path);
        }
    }
}

void VG::kpaths_of_node(Node* node, set<list<NodeTraversal> >& paths,
                        int length, int edge_max,
                        function<void(NodeTraversal)> prev_maxed, function<void(NodeTraversal)> next_maxed) {
    auto collect_path = [&paths](list<NodeTraversal>::iterator n, list<NodeTraversal>& path) {
        paths.insert(path);
    };
    for_each_kpath_of_node(node, length, edge_max, prev_maxed, next_maxed, collect_path);
}

void VG::kpaths_of_node(Node* node, vector<Path>& paths,
                        int length, int edge_max,
                        function<void(NodeTraversal)> prev_maxed, function<void(NodeTraversal)> next_maxed) {
    set<list<NodeTraversal> > unique_paths;
    kpaths_of_node(node, unique_paths, length, edge_max, prev_maxed, next_maxed);
    for (auto& unique_path : unique_paths) {
        Path path = create_path(unique_path);
        paths.push_back(path);
    }
}

// aggregators, when a callback won't work

void VG::kpaths(set<list<NodeTraversal> >& paths, int length, int edge_max,
                function<void(NodeTraversal)> prev_maxed, function<void(NodeTraversal)> next_maxed) {
    for (int i = 0; i < graph.node_size(); ++i) {
        Node* node = graph.mutable_node(i);
        kpaths_of_node(node, paths, length, edge_max, prev_maxed, next_maxed);
    }
}

void VG::kpaths(vector<Path>& paths, int length, int edge_max,
                function<void(NodeTraversal)> prev_maxed, function<void(NodeTraversal)> next_maxed) {
    set<list<NodeTraversal> > unique_paths;
    kpaths(unique_paths, length, edge_max, prev_maxed, next_maxed);
    for (auto& unique_path : unique_paths) {
        Path path = create_path(unique_path);
        paths.push_back(path);
    }
}

// path utilities
// these are in this class because attributes of the path (such as its sequence) are a property of the graph

Path VG::create_path(const list<NodeTraversal>& nodes) {
    Path path;
    for (const NodeTraversal& n : nodes) {
        Mapping* mapping = path.add_mapping();
        mapping->mutable_position()->set_node_id(n.node->id());
        // If the node is backward along this path, note it.
        if(n.backward) mapping->mutable_position()->set_is_reverse(n.backward);
        // TODO: Is it OK if we say we're at a mapping at offset 0 of this node, backwards? Or should we offset ourselves to the end?
    }
    return path;
}

string VG::path_string(const list<NodeTraversal>& nodes) {
    string seq;
    for (const NodeTraversal& n : nodes) {
        if(n.backward) {
            seq.append(reverse_complement(n.node->sequence()));
        } else {
            seq.append(n.node->sequence());
        }
    }
    return seq;
}

string VG::path_string(const Path& path) {
    string seq;
    for (int i = 0; i < path.mapping_size(); ++i) {
        auto& m = path.mapping(i);
        Node* n = node_by_id[m.position().node_id()];
        seq.append(mapping_sequence(m, *n));
    }
    return seq;
}

void VG::expand_path(const list<NodeTraversal>& path, vector<NodeTraversal>& expanded) {
    for (list<NodeTraversal>::const_iterator n = path.begin(); n != path.end(); ++n) {
        NodeTraversal node = *n;
        int s = node.node->sequence().size();
        for (int i = 0; i < s; ++i) {
            expanded.push_back(node);
        }
    }
}

void VG::expand_path(list<NodeTraversal>& path, vector<list<NodeTraversal>::iterator>& expanded) {
    for (list<NodeTraversal>::iterator n = path.begin(); n != path.end(); ++n) {
        int s = (*n).node->sequence().size();
        for (int i = 0; i < s; ++i) {
            expanded.push_back(n);
        }
    }
}

// VG::edit_node takes a node and a set of mappings (consecutive sequences of edits starting at some
// offset in a node), modifies the node, and records its translation of the node space in cut_trans.
// We use cut_trans in VG::edit to handle multi-node events, such as deletions and jumps that
// traverse more than one node.
// Mappings are provided as a tuple of <offset, mapping, start/end indicator>, where the
// start/end indicator lets us know if the mapping may contain a soft clip from the start (-1) or
// end (1) of a path. (0) indicates full containment of the mapping in the path we are including.
// Noet that the caller will have to clear and regenerate the path ranks.

void VG::edit_node(id_t node_id,
                   const vector<tuple<Mapping, bool, bool> >& mappings,
                   map<pair<id_t, size_t>, pair<set<Node*>, set<Node*>>>& cut_trans) {

    //cerr << "editing node " << node_id << endl;

    // assume these are all for the same node
    //id_t node_id = mappings.front().node_id();
    // convert edits to a series of cut points
    set<int> cut_at;
    // map pairs of cut points to sequences of named paths joining the cut points
    // instead of just storing the sequence, we store the edit and whether or not
    // it lies at the start (-1) or end (1) of the read, or in the middle (0)
    map<pair<int, int>, map<string, tuple<Edit, bool, bool>>> cut_seqs;
    Node* node = get_node(node_id);
    string node_seq = node->sequence();
    for (auto& m : mappings) {
        auto& mapping = get<0>(m);
        bool at_start = get<1>(m);
        bool at_end = get<2>(m);
        //cerr << "editing on node " << node_id << " with mapping " << pb2json(mapping) << endl;
        // check that we're really working on one node
        assert(mapping.position().node_id() == node_id);
        size_t offset = mapping.position().offset();
        const auto& nitr = mapping.metadata().info().find("path_name");
        assert(nitr != mapping.metadata().info().end());
        const string& name = nitr->second.str();
        for (int i = 0; i < mapping.edit_size(); ++i) {
            const Edit& edit = mapping.edit(i);
            //edits[offset].push_back(edit);
            int end = offset + edit.from_length();
            //if (edit_is_match(edit) || edit_is_softclip(edit)) {
            if (edit.sequence().empty() && edit.from_length() == edit.to_length()) {
                // match, ignore
                // to handle soft clips
            } else if ((i == 0 || i == mapping.edit_size() - 1) && edit.from_length() == 0) {
                // end == offset by definition
                cut_at.insert(offset);
                cut_seqs[make_pair(offset, end)][name]
                    = make_tuple(edit,
                                 at_start && i == 0,
                                 at_end && i == mapping.edit_size()-1);
            } else if (edit_is_insertion(edit)) {
                //cerr << "found an insertion " << pb2json(edit) << endl;
                cut_at.insert(offset);
                cut_seqs[make_pair(offset, offset)][name]
                    = make_tuple(edit,
                                 at_start && i == 0,
                                 at_end && i == mapping.edit_size()-1);
            } else { //if (edit_is_deletion(edit)) {
                //cerr << "o/e " << offset << " " << end << endl;
                // should this be forced?
                cut_seqs[make_pair(offset, end)][name] = make_tuple(edit, false, false);
                cut_at.insert(offset);
                cut_at.insert(end);
            }
            offset = end;
        }
    }

    // here come the cuts
    // boom, dis is for debugging
    /*
    for (auto cut : cut_at) {
        cerr << "cut_at: " << cut << endl;
    }
    //     map<pair<int, int>, map<string, tuple<Edit, bool, bool>>> cut_seqs;
    for (auto cut : cut_seqs) {
        cerr << "cut_seq: " << cut.first.first << "," << cut.first.second << " ";
        for (auto s : cut.second) {
            cerr << " " << s.first << "," << pb2json(get<0>(s.second)) << " ";
        }
        cerr << endl;
    }
    */

    // tricky bit
    // we need to put the node pointers into structures to track the sides of each cut
    // however, as we make the subsequent cut we end up changing the nodes destructively

    Node* n = get_node(node_id);
    Node* l = NULL;
    Node* r = NULL;
    int offset = 0;

    map<int, tuple<set<Node*>, set<Node*>, set<Node*> > > cuts;

    // make the cuts
    set<Node*> emptyset;
    for (auto cut : cut_at) {
        //cerr << "cuttin at " << cut-offset << endl;
        // don't cut if it we already have the cut as the node end
        if (cut != 0 && cut != node_seq.size()) {
            divide_node(n, cut-offset, l, r);
            set<Node*> left{l};
            set<Node*> right{r};
            //cerr << "recording " << node_id << ":" << cut << " left " << pb2json(*l) << endl;
            //cerr << "recording " << node_id << ":" << cut << " right " << pb2json(*r) << endl;
            cuts[cut] = make_tuple(left, emptyset, right);
            // now that we've made the cuts record the mapping between the old node id and the new
            if (offset > 0) {
                // sets the last right to the current left
                get<2>(cuts[offset]) = left;
            }
            offset = cut;
            n = r;
        }
    }

    // now that we've made the cuts, store the translations for the node ends
    // we don't handle them below
    Node* start = nullptr;
    Node* end = nullptr;
    if (!cuts.empty()) {
        //cerr << "recording cut trans of " << node_id << ":" << 0 << " " << pb2json(*n) << " " << n << endl;
        auto& p1 = cut_trans[make_pair(node_id, 0)];
        start = *get<0>(cuts.begin()->second).begin();
        p1.second.insert(start);
        auto& p2 = cut_trans[make_pair(node_id, node_seq.size())];
        end = *get<2>(cuts.rbegin()->second).begin();
        p2.first.insert(end);
    } else {
        auto& p1 = cut_trans[make_pair(node_id, 0)];
        p1.second.insert(n);
        start = n;
        auto& p2 = cut_trans[make_pair(node_id, node_seq.size())];
        p2.first.insert(n);
        end = n;
    }

    // add the novel seqs
    // and join deletions
    for (auto& cs : cut_seqs) {
        // get the left and right
        int f = cs.first.first; // we go from after the left node at this cut
        int t = cs.first.second; // to before the second node in this cut
        auto& left = get<0>(cuts[f]);
        auto& center = get<1>(cuts[f]);
        auto& right = get<2>(cuts[t]);
        //cerr << f << "," << t << " " << left.size() << "," << center.size() << "," << right.size() << endl;
        for (auto& p : cs.second) {
            auto& name = get<0>(p);
            auto& edit = get<0>(get<1>(p));
            if (!edit.sequence().empty()) {
                Node* c = create_node(edit.sequence());
                // add to path
                center.insert(c);
            }
        }
    }

    // join in the novel seqs
    for (auto& cs : cut_seqs) {
        // get the left and right
        int f = cs.first.first; // we go from after the left node at this cut
        int t = cs.first.second; // to before the second node in this cut
        auto& left = get<0>(cuts[f]);
        auto& center = get<1>(cuts[f]);
        auto& right = get<2>(cuts[t]);
        //cerr << "cut seq " << f << ":" << t << " "
        //     << left.size() << "," << center.size() << "," << right.size() << endl;
        for (auto& e : cs.second) {
            const string& name = e.first;
            const Edit& edit = get<0>(e.second);
            bool at_start = get<1>(e.second);
            bool at_end = get<2>(e.second);
            // handle the case where we don't have any left or right nodes
            // such as when we are a soft clip dangling off the left or right end of the node
            if (at_start && start != nullptr && right.empty()) {
                for (auto& cp : center) {
                    create_edge(cp, start);
                }
            }
            // same for soft clips off end
            if (at_end && end != nullptr && left.empty()) {
                for (auto& cp : center) {
                    create_edge(end, cp);
                }
            }
            // otherwise, use the cut set to make the joins of novel sequences
            for (auto& lp : left) {
                for (auto& rp : right) {
                    if (!center.empty()) {
                        for (auto& cp : center) {
                            // make the edge across the left cut
                            if (!at_start) {
                                // we are not on a left-dangling node (from soft clip)
                                create_edge(lp, cp);
                            }
                            // record the cut translation (used for stitching deletions)
                            auto& p1 = cut_trans[make_pair(node_id, f)];
                            p1.first.insert(lp); p1.second.insert(cp);
                            //cerr << *p1.first.begin() << " " << *p1.second.begin() << endl;
                            // make the edge across the right cut
                            if (!at_end) {
                                // we are not on a right-dangling node (from soft clip)
                                create_edge(cp, rp);
                            }
                            // record the cut translation ...
                            auto& p2 = cut_trans[make_pair(node_id, t)];
                            p2.first.insert(cp); p2.second.insert(rp);
                            //cerr << *p2.first.begin() << " " << *p2.second.begin() << endl;
                        }
                    } else {
                        // simpler case when we don't have a center node
                        // make the edge across the single cut
                        // and record the translation
                        auto& p1 = cut_trans[make_pair(node_id, f)];
                        //cerr << "cut trans " << node_id << ":" << f << endl;
                        //cerr << pb2json(*lp) << " -> " << pb2json(*rp) << endl;
                        p1.first.insert(lp);
                        p1.second.insert(rp);
                        if (f != t) {
                            //cerr << "f != t" << endl;
                            //cerr << "cut trans " << node_id << ":" << t << endl;
                            auto& p2 = cut_trans[make_pair(node_id, t)];
                            //cerr << pb2json(*lp) << " -> " << pb2json(*rp) << endl;
                            p2.first.insert(lp);
                            p2.second.insert(rp);
                        }
                    }
                }
            }
        }
    }

    // now we should have incorporated the edits against a given node
    // deletions spanning multiple nodes must be handled externally using cut_trans
}

// edit the graph, so that the graph will include the paths which are provided
// beware:
// this doesn't work on the reverse strand
// to get this to work, we should be using NodeSides rather than pair<id_t, size_t>
// to track cuts and positions in the system
void VG::edit(const vector<Path>& paths) {
#ifdef debug
    // We should not add additional disconnected components
    list<VG> subgraphs;
    disjoint_subgraphs(subgraphs);
    int subgraph_count = subgraphs.size();
#endif


    //cerr << "editing graph" << endl;
    // deletions from this node position/offset on the forward strand to the second position/offset
    map<pair<id_t, size_t>, pair<id_t, size_t> > del_f;
    // reverse index of the above
    map<pair<id_t, size_t>, pair<id_t, size_t> > del_t;
    // mappings by node, each is stored as a tuple of <mapping, is_at_start, is_at_end>
    map<id_t, vector<tuple<Mapping, bool, bool>>> mappings; // by node
    // the cut_trans(lation) maps from positions in the old graph (pair<size_t, id_t>)
    // to tuples of nodes that sets of nodes that now lie on either side of these
    // to handle long jumps with sequence on the path
    map<pair<id_t, size_t>, pair<set<Node*>, set<Node*> > > cut_trans;
    bool in_del = false;
    pair<id_t, size_t> del_start;
    for (auto& p : paths) {
        auto path = simplify(p);
        for (int i = 0; i < path.mapping_size(); ++i) {
            Mapping mapping = path.mapping(i);
            //if (!mapping.has_position()) cerr << "woah mapping got no position " << pb2json(mapping) << endl;
            if (i > 0) {
                auto& last = path.mapping(i-1);
                // do we go from in-node to off-node between this mapping and the next
                // without reaching the end of the node we're on?
                // this would indicate a jump!
                id_t lid = last.position().node_id();
                id_t nid = mapping.position().node_id();
                //cerr << "there is a possible jump " << lid << " -> " << nid << endl;

                // now it's possible that we are cutting the two nodes as well with the jump
                // to make the cut, we have to find the mapping end offset on the first node
                size_t loff = last.position().offset() + mapping_from_length(last);
                // and the mapping start offset on the second
                size_t noff = mapping.position().offset();
                // should there be sequence on the path?
                // get the sequence after the loff and before the noff
                auto lafter = cut_mapping(last, loff).second;
                auto nbefore = cut_mapping(mapping, noff).first;
                string ins_seq;
                if (mapping_to_length(lafter) + mapping_to_length(nbefore)) {
                    //cerr << "there be an insertion!" << endl;
                    // path has been simplified, each insertion will be == 1 edit
                    ins_seq = (lafter.edit_size()?lafter.edit(lafter.edit_size()-1).sequence():"")
                        + (nbefore.edit_size()?nbefore.edit(0).sequence():"");
                }
                // record the jumps and cut trans in del_f and del_t
                // these they will be included at the end of this function
                if (ins_seq.empty()) {
                    auto dstart = make_pair(lid, loff);
                    auto dend =   make_pair(nid, noff);
                    del_f[dstart] = dend;
                    del_t[dend] = dstart;
                } else {
                    //cerr << "creating insertion node" << endl;
                    Node* ins = create_node(ins_seq);
                    auto dstart1 = make_pair(lid, loff);
                    auto dend1 = make_pair(ins->id(), 0);
                    del_f[dstart1] = dend1;
                    del_t[dend1] = dstart1;
                    auto& p1 = cut_trans[make_pair(lid, loff)];
                    p1.second.insert(ins);

                    // record the middle of the insertion in the trans
                    auto& p1m = cut_trans[make_pair(ins->id(), 0)];
                    p1m.second.insert(ins);
                    auto& p2m = cut_trans[make_pair(ins->id(), ins->sequence().size())];
                    p2m.first.insert(ins);

                    auto dstart2 = make_pair(ins->id(), ins_seq.size());
                    auto dend2 = make_pair(nid, noff);
                    del_f[dstart2] = dend2;
                    del_t[dend2] = dstart2;
                    auto& p2 = cut_trans[make_pair(nid, noff)];
                    p2.first.insert(ins);
                }
                // trigger the appropriate node-level edits by making two new mappings
                // against the nodes which include 1-bp no-seq insertions at the correct position
                // these won't be included but trigger the node cut
                Info info; info.set_str(path.name());
                Mapping lmap;
                {
                    lmap.mutable_position()->set_node_id(lid);
                    (*lmap.mutable_metadata()->mutable_info())["path_name"] = info;
                    Edit* edit = lmap.add_edit();
                    edit->set_to_length(loff);
                    edit->set_from_length(loff);
                    edit = lmap.add_edit();
                    edit->set_from_length(0);
                    edit->set_to_length(1);
                }
                Mapping nmap;
                {
                    nmap.mutable_position()->set_node_id(nid);
                    (*nmap.mutable_metadata()->mutable_info())["path_name"] = info;
                    Edit* edit = nmap.add_edit();
                    edit->set_to_length(noff);
                    edit->set_from_length(noff);
                    edit = nmap.add_edit();
                    edit->set_from_length(0);
                    edit->set_to_length(1);
                }
                //cerr << "lmap is " << pb2json(lmap) << " and nmap is " << pb2json(nmap) << endl;
                // and save them
                // note that we store them as being "from the start" only
                mappings[last.position().node_id()].push_back(make_tuple(lmap, true, false));
                mappings[mapping.position().node_id()].push_back(make_tuple(nmap, true, false));
            }
            Info info; info.set_str(path.name());
            (*mapping.mutable_metadata()->mutable_info())["path_name"] = info;
            // store in our mapping set, with flags for whether we're at the start or end
            mappings[mapping.position().node_id()].push_back(make_tuple(
                                                                 mapping,
                                                                 i==0,
                                                                 i==path.mapping_size()-1
                                                                 ));
            if (mapping_is_total_deletion(mapping)) {
                //&& path.mapping_size() > 1 && i != 0 && i != path.mapping_size()-1) {
                cerr << "[vg::VG] warning: mapping in " << path.name() << " is a total deletion" << endl;
            }
            // ^^ these should be stripped out using simplify_deletions(path)
            id_t node_id = mapping.position().node_id();

            // find internal deletions
            size_t ioff = 0;
            for (int j = 0; j < mapping.edit_size(); ++j) {
                auto& edit = mapping.edit(j);
                // is it an internal deletion?
                if (edit_is_deletion(edit) && j > 0 && j < mapping.edit_size()-1) {
                    auto dstart = make_pair(node_id, ioff);
                    auto dend =   make_pair(node_id, ioff+edit.from_length());
                    del_f[dstart] = dend;
                    del_t[dend] = dstart;
                }
                ioff += edit.from_length();
            }

            // not just ones that jump nodes
            if (!in_del) {
                if (mapping_starts_in_deletion(mapping) && i > 0) {
                    // this deletion should extend from the end of the last node to this one
                    id_t last_id = path.mapping(i-1).position().node_id();
                    del_start = make_pair(last_id, get_node(last_id)->sequence().size());
                }
                if (mapping_ends_in_deletion(mapping)) {
                    size_t del_off = get_node(node_id)->sequence().size()
                        - mapping.edit(mapping.edit_size()-1).from_length();
                    del_start = make_pair(node_id, del_off);
                    in_del = true;
                }
            } else {
                pair<id_t, size_t> del_end;
                if (mapping_starts_in_deletion(mapping)) {
                    size_t end_off = mapping.edit(0).from_length();
                    del_end = make_pair(node_id, end_off);
                    in_del = false;
                } else {
                    del_end = make_pair(node_id, 0);
                    in_del = false;
                }
                // record the deletion
                //cerr << "del_start " << del_start.first << ":" << del_start.second << " -> " << del_end.first << ":" << del_end.second << endl;
                del_f[del_start] = del_end;
                del_t[del_end] = del_start;
            }
        }
    }
    edit(mappings, cut_trans, del_f, del_t);

#ifdef debug
    // We should not add additional disconnected components
    subgraphs.clear();
    disjoint_subgraphs(subgraphs);
    int end_subgraph_count = subgraphs.size();
    assert(end_subgraph_count == subgraph_count);
#endif
}

// mappings sorted by node id
void VG::edit(const map<id_t, vector<tuple<Mapping, bool, bool> > >& mappings,
              map<pair<id_t, size_t>, pair<set<Node*>, set<Node*> > >& cut_trans,
              map<pair<id_t, size_t>, pair<id_t, size_t> >& del_f,
              map<pair<id_t, size_t>, pair<id_t, size_t> >& del_t) {

    map<id_t, pair<Node*, Node*> > node_end_map;
    for (auto& node : mappings) {
        //cerr << "edit: " << node.first << endl;
        edit_node(node.first, node.second, cut_trans);
    }

    // now resolve the node jumps (trans-insertions and deletions)
    //cerr << "resolving trans events" << endl;
    for (auto& d : del_f) {
        auto& f = d.first;
        auto& t = d.second;
        auto& l = cut_trans[f];
        auto& r = cut_trans[t];
        //cerr << f.first << ":" << f.second << " " << l.first.size() << " " << l.second.size() << endl;
        //cerr << t.first << ":" << t.second << " " << r.first.size() << " " << r.second.size() << endl;
        // connect all left to right
        for (auto& ln : l.first) {
            //cerr << "left node is " << pb2json(*ln) << " " << ln << endl;
            for (auto& rn : r.second) {
                //cerr << "right node is " << pb2json(*rn) << " " << rn << endl;
                create_edge(ln, rn);
            }
        }
    }
    remove_null_nodes_forwarding_edges();
}

// The correct way to edit the graph
void VG::edit_both_directions(const vector<Path>& paths_to_add) {
    // Collect the breakpoints
    map<id_t, set<pos_t>> breakpoints;

#ifdef debug
    for (auto& p : paths_to_add) {
        cerr << pb2json(p) << endl;
    }
#endif

    std::vector<Path> simplified_paths;

    for(auto path : paths_to_add) {
        // Simplify the path, just to eliminate adjacent match Edits in the same
        // Mapping (because we don't have or want a breakpoint there)
        simplified_paths.push_back(simplify(path));
    }

    for(auto path : simplified_paths) {
        // Add in breakpoints from each path
        find_breakpoints(path, breakpoints);
    }

    // Invert the breakpoints that are on the reverse strand
    breakpoints = forwardize_breakpoints(breakpoints);

    // Clear existing path ranks.
    paths.clear_mapping_ranks();

    // Break any nodes that need to be broken. Save the map we need to translate
    // from offsets on old nodes to new nodes. Note that this would mess up the
    // ranks of nodes in their existing paths, which is why we clear and rebuild
    // them.
    auto node_translation = ensure_breakpoints(breakpoints);

    for(auto path : simplified_paths) {
        // Now go through each new path again, and create new nodes/wire things up.
        add_nodes_and_edges(path, node_translation);
    }

    // TODO: add the new path to the graph, with perfect match mappings to all
    // the new and old stuff it visits.

    // Rebuild path ranks, aux mapping, etc. by compacting the path ranks
    paths.compact_ranks();

    // with the paths sorted, let's double-check that the edges are here
    paths.for_each([&](const Path& path) {
            for (size_t i = 1; i < path.mapping_size(); ++i) {
                auto& m1 = path.mapping(i-1);
                auto& m2 = path.mapping(i);
                if (!adjacent_mappings(m1, m2)) continue; // the path is completely represented here
                auto s1 = NodeSide(m1.position().node_id(), (m1.position().is_reverse() ? false : true));
                auto s2 = NodeSide(m2.position().node_id(), (m2.position().is_reverse() ? true : false));
                // check that we always have an edge between the two nodes in the correct direction
                if (!has_edge(s1, s2)) {
                    cerr << "graph path '" << path.name() << "' invalid: edge from "
                         << s1 << " to " << s2 << " does not exist" << endl;
                    cerr << "creating edge" << endl;
                    create_edge(s1, s2);
                }
            }
        });

    // execute a semi partial order sort on the nodes
    sort();
}

map<id_t, set<pos_t>> VG::forwardize_breakpoints(const map<id_t, set<pos_t>>& breakpoints) {
    map<id_t, set<pos_t>> fwd;
    for (auto& p : breakpoints) {
        id_t node_id = p.first;
        assert(has_node(node_id));
        size_t node_length = get_node(node_id)->sequence().size();
        auto bp = p.second;
        for (auto& pos : bp) {
            pos_t x = pos;
            if (offset(pos) == node_length) continue;
            if (offset(pos) > node_length) {
                cerr << "VG::forwardize_breakpoints error: failure, position " << pos << " is not inside node "
                     << pb2json(*get_node(node_id)) << endl;
                assert(false);
            }
            if (is_rev(pos)) {
                fwd[node_id].insert(reverse(pos, node_length));
            } else {
                fwd[node_id].insert(pos);
            }
        }
    }
    return fwd;
}

// returns breakpoints on the forward strand of the nodes
void VG::find_breakpoints(const Path& path, map<id_t, set<pos_t>>& breakpoints) {
    // We need to work out what offsets we will need to break each node at, if
    // we want to add in all the new material and edges in this path.

#ifdef debug
    cerr << "Processing path..." << endl;
#endif

    for (size_t i = 0; i < path.mapping_size(); ++i) {
        // For each Mapping in the path
        const Mapping& m = path.mapping(i);

        // What node are we on?
        id_t node_id = m.position().node_id();

        if(node_id == 0) {
            // Skip Mappings that aren't actually to nodes.
            continue;
        }

        // See where the next edit starts in the node. It is always included
        // (even when the edit runs backward), unless the edit has 0 length in
        // the reference.
        pos_t edit_first_position = make_pos_t(m.position());

#ifdef debug
        cerr << "Processing mapping " << pb2json(m) << endl;
#endif

        for(size_t j = 0; j < m.edit_size(); ++j) {
            // For each Edit in the mapping
            const Edit& e = m.edit(j);

            // We know where the mapping starts in its node. But where does it
            // end (inclusive)? Note that if the edit has 0 reference length,
            // this may not actually be included in the edit (and
            // edit_first_position will be further along than
            // edit_last_position).
            pos_t edit_last_position = edit_first_position;
            if (e.from_length()) {
                get_offset(edit_last_position) += e.from_length();
            }

#ifdef debug
            cerr << "Edit on " << node_id << " from " << edit_first_position << " to " << edit_last_position << endl;
            cerr << pb2json(e) << endl;
#endif

            if (!edit_is_match(e) || j == 0) {
                // If this edit is not a perfect match, or if this is the first
                // edit in this mapping and we had a previous mapping we may
                // need to connect to, we need to make sure we have a breakpoint
                // at the start of this edit.

#ifdef debug
                cerr << "Need to break " << node_id << " at edit lower end " <<
                    edit_first_position << endl;
#endif

                // We need to snip between edit_first_position and edit_first_position - direction.
                // Note that it doesn't matter if we put breakpoints at 0 and 1-past-the-end; those will be ignored.
                breakpoints[node_id].insert(edit_first_position);
            }

            if (!edit_is_match(e) || (j == m.edit_size() - 1)) {
                // If this edit is not a perfect match, or if it is the last
                // edit in a mapping and we have a subsequent mapping we might
                // need to connect to, make sure we have a breakpoint at the end
                // of this edit.

#ifdef debug
                cerr << "Need to break " << node_id << " at past edit upper end " <<
                    edit_last_position << endl;
#endif

                // We also need to snip between edit_last_position and edit_last_position + direction.
                breakpoints[node_id].insert(edit_last_position);
            }

            // TODO: for an insertion or substitution, note that we need a new
            // node and two new edges.

            // TODO: for a deletion, note that we need an edge. TODO: Catch
            // and complain about some things we can't handle (like a path with
            // a leading/trailing deletion)? Or just skip deletions when wiring.

            // Use up the portion of the node taken by this mapping, so we know
            // where the next mapping will start.
            edit_first_position = edit_last_position;
        }
    }

}

map<pos_t, Node*> VG::ensure_breakpoints(const map<id_t, set<pos_t>>& breakpoints) {
    // Set up the map we will fill in with the new node start positions in the
    // old nodes.
    map<pos_t, Node*> toReturn;

    for(auto& kv : breakpoints) {
        // Go through all the nodes we need to break up
        auto original_node_id = kv.first;

        // Save the original node length. We don;t want to break here (or later)
        // because that would be off the end.
        id_t original_node_length = get_node(original_node_id)->sequence().size();

        // We are going through the breakpoints left to right, so we need to
        // keep the node pointer for the right part that still needs further
        // dividing.
        Node* right_part = get_node(original_node_id);
        Node* left_part = nullptr;

        pos_t last_bp = make_pos_t(original_node_id, false, 0);
        // How far into the original node does our right part start?
        id_t current_offset = 0;

        for(auto breakpoint : kv.second) {
            // For every point at which we need to make a new node, in ascending
            // order (due to the way sets store ints)...

            // ensure that we're on the forward strand (should be the case due to forwardize_breakpoints)
            assert(!is_rev(breakpoint));

            // This breakpoint already exists, because the node starts or ends here
            if(offset(breakpoint) == 0
               || offset(breakpoint) == original_node_length) {
                continue;
            }

            // How far in do we need to break the remaining right part? And how
            // many bases will be in this new left part?
            id_t divide_offset = offset(breakpoint) - current_offset;


#ifdef debug
            cerr << "Need to divide original " << original_node_id << " at " << breakpoint << "/" <<

                original_node_length << endl;
            cerr << "Translates to " << right_part->id() << " at " << divide_offset << "/" <<
                right_part->sequence().size() << endl;
            cerr << "divide offset is " << divide_offset << endl;
#endif

            if (offset(breakpoint) <= 0) { cerr << "breakpoint is " << breakpoint << endl; }
            assert(offset(breakpoint) > 0);
            if (offset(breakpoint) >= original_node_length) { cerr << "breakpoint is " << breakpoint << endl; }
            assert(offset(breakpoint) < original_node_length);

            // Make a new left part and right part. This updates all the
            // existing perfect match paths in the graph.
            divide_node(right_part, divide_offset, left_part, right_part);

#ifdef debug

            cerr << "Produced " << left_part->id() << " (" << left_part->sequence().size() << " bp)" << endl;
            cerr << "Left " << right_part->id() << " (" << right_part->sequence().size() << " bp)" << endl;
#endif

            // The left part is now done. We know it started at current_offset
            // and ended before breakpoint, so record it by start position.

            // record forward and reverse
            toReturn[last_bp] = left_part;
            toReturn[reverse(breakpoint, original_node_length)] = left_part;

            // Record that more sequence has been consumed
            current_offset += divide_offset;
            last_bp = breakpoint;

        }

        // Now the right part is done too. It's going to be the part
        // corresponding to the remainder of the original node.
        toReturn[last_bp] = right_part;
        toReturn[make_pos_t(original_node_id, true, 0)] = right_part;

        // and record the start and end of the node
        toReturn[make_pos_t(original_node_id, true, original_node_length)] = nullptr;
        toReturn[make_pos_t(original_node_id, false, original_node_length)] = nullptr;

    }

    return toReturn;
}

void VG::add_nodes_and_edges(const Path& path, const map<pos_t, Node*>& node_translation) {
    // The basic algorithm is to traverse the path edit by edit, keeping track
    // of a NodeSide for the last piece of sequence we were on. If we hit an
    // edit that creates new sequence, we create that new sequence as a node,
    // and attach it to the dangling NodeSide, and leave its end dangling. If we
    // hit an edit that corresponds to a match, we know that there's a
    // breakpoint on each end (since it's bordered by a non-perfect-match or the
    // end of a node), so we can attach its start to the dangling NodeSide and
    // leave its end dangling.

    // We need node_translation to translate between node ID space, where the
    // paths are articulated, and new node ID space, where the edges are being
    // made.

    // We use this function to get the node that contains a position on an
    // original node.
    /*
    cerr << "node translation" << endl;
    for (auto& p : node_translation) {
        pos_t pos = p.first;
        cerr << pos << " " << (p.second != nullptr?pb2json(*p.second):"null") << endl;
    }
    */

    auto find_new_node = [&](pos_t old_pos) {
        if(node_translation.find(make_pos_t(id(old_pos), false, 0)) == node_translation.end()) {
            // The node is unchanged
            auto n = get_node(id(old_pos));
            assert(n != nullptr);
            return n;
        }
        // Otherwise, get the first new node starting after that position, and
        // then look left.
        auto found = node_translation.upper_bound(old_pos);
        assert(found != node_translation.end());
        if (id(found->first) != id(old_pos)
            || is_rev(found->first) != is_rev(old_pos)) {
            return (Node*)nullptr;
        }
        // Get the thing before that (last key <= the position we want
        --found;
        assert(found->second != nullptr);

        // Return the node we found.
        return found->second;
    };

    auto create_new_mappings = [&](pos_t p1, pos_t p2, bool is_rev) {
        vector<Mapping> mappings;
        vector<Node*> nodes;
        for (pos_t p = p1; p <= p2; ++get_offset(p)) {
            auto n = find_new_node(p);
            assert(n != nullptr);
            nodes.push_back(find_new_node(p));
        }
        auto np = nodes.begin();
        while (np != nodes.end()) {
            size_t c = 0;
            auto n1 = np;
            while (np != nodes.end() && *n1 == *np) {
                ++c;
                ++np; // we'll always increment once
            }
            assert(c);
            // set the mapping position
            Mapping m;
            m.mutable_position()->set_node_id((*n1)->id());
            m.mutable_position()->set_is_reverse(is_rev);
            // and the edit that says we match
            Edit* e = m.add_edit();
            e->set_from_length(c);
            e->set_to_length(c);
            mappings.push_back(m);
        }
        return mappings;
    };

    // What's dangling and waiting to be attached to? In current node ID space.
    // We use the default constructed one (id 0) as a placeholder.
    NodeSide dangling;

    for (size_t i = 0; i < path.mapping_size(); ++i) {
        // For each Mapping in the path
        const Mapping& m = path.mapping(i);

        // What node are we on? In old node ID space.
        id_t node_id = m.position().node_id();

        // See where the next edit starts in the node. It is always included
        // (even when the edit runs backward), unless the edit has 0 length in
        // the reference.
        pos_t edit_first_position = make_pos_t(m.position());

        for(size_t j = 0; j < m.edit_size(); ++j) {
            // For each Edit in the mapping
            const Edit& e = m.edit(j);

            // Work out where its end position on the original node is (inclusive)
            // We don't use this on insertions, so 0-from-length edits don't matter.
            pos_t edit_last_position = edit_first_position;
            //get_offset(edit_last_position) += (e.from_length()?e.from_length()-1:0);
            get_offset(edit_last_position) += (e.from_length()?e.from_length()-1:0);

//#define debug_edit true
#ifdef debug_edit
            cerr << "Edit on " << node_id << " from " << edit_first_position << " to " << edit_last_position << endl;
            cerr << pb2json(e) << endl;
#endif

            if(edit_is_insertion(e) || edit_is_sub(e)) {
                // This edit introduces new sequence.
#ifdef debug_edit
                cerr << "Handling ins/sub relative to " << node_id << endl;
#endif
                // Create the new node, reversing it if we are reversed
                Node* new_node = create_node(
                    m.position().is_reverse() ?
                    reverse_complement(e.sequence())
                    : e.sequence());

                if (!path.name().empty()) {
                    Mapping nm;
                    nm.mutable_position()->set_node_id(new_node->id());
                    nm.mutable_position()->set_is_reverse(m.position().is_reverse());
                    nm.set_rank(m.rank());
                    Edit* e = nm.add_edit();
                    size_t l = new_node->sequence().size();
                    e->set_from_length(l);
                    e->set_to_length(l);
                    // insert the mapping at the right place
                    paths.append_mapping(path.name(), nm);
                }

                if(dangling.node) {
                    // This actually referrs to a node.
#ifdef debug_edit
                    cerr << "Connecting " << dangling << " and " << NodeSide(new_node->id(), m.position().is_reverse()) << endl;
#endif
                    // Add an edge from the dangling NodeSide to the start of this new node
                    assert(create_edge(dangling, NodeSide(new_node->id(), m.position().is_reverse())));

                }

                // Dangle the end of this new node
                dangling = NodeSide(new_node->id(), !m.position().is_reverse());

                // save edit into translated path

            } else if(edit_is_match(e)) {
                // We're using existing sequence

                // We know we have breakpoints on both sides, but we also might
                // have additional breakpoints in the middle. So we need the
                // left node, that contains the first base of the match, and the
                // right node, that contains the last base of the match.
                Node* left_node = find_new_node(edit_first_position);
                Node* right_node = find_new_node(edit_last_position);

                // TODO: we just assume the outer edges of these nodes are in
                // the right places. They should be if we cut the breakpoints
                // right.

                // TODO include path
                // get the set of new nodes that we map to
                // and use the lengths of each to create new mappings
                // and append them to the path we are including
                if (!path.name().empty()) {
                    for (auto nm : create_new_mappings(edit_first_position,
                                                       edit_last_position,
                                                       m.position().is_reverse())) {
                        //cerr << "in match, adding " << pb2json(nm) << endl;
                        // no rank has been established, so get the next available
                        // otherwise we'll fail to include looping paths
                        nm.set_rank(m.rank());
                        paths.append_mapping(path.name(), nm);
                    }
                }

#ifdef debug_edit
                cerr << "Handling match relative to " << node_id << endl;
#endif

                if(dangling.node) {
#ifdef debug_edit
                    cerr << "Connecting " << dangling << " and " << NodeSide(left_node->id(), m.position().is_reverse()) << endl;
#endif

                    // Connect the left end of the left node we matched in the direction we matched it
                    assert(create_edge(dangling, NodeSide(left_node->id(), m.position().is_reverse())));
                }

                // Dangle the right end of the right node in the direction we matched it.
                if (right_node != nullptr) dangling = NodeSide(right_node->id(), !m.position().is_reverse());
            } else {
                // We don't need to deal with deletions since we'll deal with the actual match/insert edits on either side
                // Also, simplify() simplifies them out.
#ifdef debug_edit
                cerr << "Skipping other edit relative to " << node_id << endl;
#endif
            }

            // Advance in the right direction along the original node for this edit.
            // This way the next one will start at the right place.
            get_offset(edit_first_position) += e.from_length();


        }

    }

}

void VG::node_starts_in_path(const list<NodeTraversal>& path,
                             map<Node*, int>& node_start) {
    int i = 0;
    for (list<NodeTraversal>::const_iterator n = path.begin(); n != path.end(); ++n) {
        node_start[(*n).node] = i;
        int l = (*n).node->sequence().size();
        i += l;
    }
}

void VG::node_starts_in_path(list<NodeTraversal>& path,
                             map<NodeTraversal*, int>& node_start) {
    int i = 0;
    for (list<NodeTraversal>::iterator n = path.begin(); n != path.end(); ++n) {
        node_start[&(*n)] = i;
        int l = (*n).node->sequence().size();
        i += l;
    }
}

void VG::kpaths_of_node(id_t node_id, vector<Path>& paths, int length, int edge_max,
                        function<void(NodeTraversal)> prev_maxed, function<void(NodeTraversal)> next_maxed) {
    hash_map<id_t, Node*>::iterator n = node_by_id.find(node_id);
    if (n != node_by_id.end()) {
        Node* node = n->second;
        kpaths_of_node(node, paths, length, edge_max, prev_maxed, next_maxed);
    }
}

// todo record as an alignment rather than a string
Alignment VG::random_read(size_t read_len,
                                        mt19937& rng,
                                        id_t min_id,
                                        id_t max_id,
                                        bool either_strand) {
    // this is broken as it should be scaled by the sequence space
    // not node space
    // TODO BROKEN
    uniform_int_distribution<id_t> int64_dist(min_id, max_id);
    id_t id = int64_dist(rng);
    // We start at the node in its local forward orientation
    NodeTraversal node(get_node(id), false);
    int32_t start_pos = 0;
    if (node.node->sequence().size() > 1) {
        uniform_int_distribution<uint32_t> uint32_dist(0,node.node->sequence().size()-1);
        start_pos = uint32_dist(rng);
    }
    string read = node.node->sequence().substr(start_pos);
    Alignment aln;
    Path* path = aln.mutable_path();
    Mapping* mapping = path->add_mapping();
    Position* position = mapping->mutable_position();
    position->set_offset(start_pos);
    position->set_node_id(node.node->id());
    Edit* edit = mapping->add_edit();
    //edit->set_sequence(read);
    edit->set_from_length(read.size());
    edit->set_to_length(read.size());
    while (read.size() < read_len) {
        // pick a random downstream node
        vector<NodeTraversal> next_nodes;
        nodes_next(node, next_nodes);
        if (next_nodes.empty()) break;
        uniform_int_distribution<int> next_dist(0, next_nodes.size()-1);
        node = next_nodes.at(next_dist(rng));
        // Put in the node sequence in the correct relative orientation
        string addition = (node.backward
                           ? reverse_complement(node.node->sequence()) : node.node->sequence());
        read.append(addition);
        mapping = path->add_mapping();
        position = mapping->mutable_position();
        position->set_offset(0);
        position->set_node_id(node.node->id());
        edit = mapping->add_edit();
        //edit->set_sequence(addition);
        edit->set_from_length(addition.size());
        edit->set_to_length(addition.size());
    }
    aln.set_sequence(read);
    // fix up the length
    read = read.substr(0, read_len);
    size_t to_len = alignment_to_length(aln);
    if ((int)to_len - (int)read_len > 0) {
        aln = strip_from_end(aln, (int)to_len - (int)read_len);
    }
    uniform_int_distribution<int> binary_dist(0, 1);
    if (either_strand && binary_dist(rng) == 1) {
        // We can flip to the other strand (i.e. node's local reverse orientation).
        reverse_complement_alignment(aln,
                                     (function<id_t(id_t)>) ([this](id_t id) {
                                             return get_node(id)->sequence().size();
                                         }));
    }
    return aln;
}

bool VG::is_valid(bool check_nodes,
                  bool check_edges,
                  bool check_paths,
                  bool check_orphans) {

    if (check_nodes) {

        if (node_by_id.size() != graph.node_size()) {
            cerr << "graph invalid: node count is not equal to that found in node by-id index" << endl;
            return false;
        }

        for (int i = 0; i < graph.node_size(); ++i) {
            Node* n = graph.mutable_node(i);
            if (node_by_id.find(n->id()) == node_by_id.end()) {
                cerr << "graph invalid: node " << n->id() << " missing from by-id index" << endl;
                return false;
            }
        }
    }

    if (check_edges) {
        for (int i = 0; i < graph.edge_size(); ++i) {
            Edge* e = graph.mutable_edge(i);
            id_t f = e->from();
            id_t t = e->to();

            //cerr << "edge " << e << " " << e->from() << "->" << e->to() << endl;

            if (node_by_id.find(f) == node_by_id.end()) {
                cerr << "graph invalid: edge index=" << i
                     << " (" << f << "->" << t << ") cannot find node (from) " << f << endl;
                return false;
            }
            if (node_by_id.find(t) == node_by_id.end()) {
                cerr << "graph invalid: edge index=" << i
                     << " (" << f << "->" << t << ") cannot find node (to) " << t << endl;
                return false;
            }

            if (!edges_on_start.count(f) && !edges_on_end.count(f)) {
                // todo check if it's in the vector
                cerr << "graph invalid: edge index=" << i
                     << " could not find entry in either index for 'from' node " << f << endl;
                return false;
            }

            if (!edges_on_start.count(t) && !edges_on_end.count(t)) {
                // todo check if it's in the vector
                cerr << "graph invalid: edge index=" << i
                     << " could not find entry in either index for 'to' node " << t << endl;
                return false;
            }
        }

        for (pair<const id_t, vector<pair<id_t, bool>>>& start_and_edges : edges_on_start) {
            for (auto& edge_destination : start_and_edges.second) {
                // We're on the start, so we go to the end if we aren't a reversing edge
                Edge* e = get_edge(NodeSide::pair_from_start_edge(start_and_edges.first, edge_destination));
                if (!e) {
                    cerr << "graph invalid, edge is null" << endl;
                    return false;
                }
                if(start_and_edges.first != e->to() && start_and_edges.first != e->from()) {
                    // It needs to be attached to the node we looked up
                    cerr << "graph invalid: edge " << e->from() << "->" << e->to()
                         << " doesn't have start-indexed node in " << start_and_edges.first << "<->"
                         << edge_destination.first << endl;
                    return false;
                }
                if(edge_destination.first != e->to() && edge_destination.first != e->from()) {
                    // It also needs to be attached to the node it says it goes to
                    cerr << "graph invalid: edge " << e->from() << "->" << e->to()
                         << " doesn't have non-start-indexed node in " << start_and_edges.first << "<->"
                         << edge_destination.first << endl;
                    return false;
                }
                if(!((start_and_edges.first == e->to() && !e->to_end()) ||
                     (start_and_edges.first == e->from() && e->from_start()))) {

                    // The edge needs to actually attach to the start of the node we looked it up for.
                    // So at least one of its ends has to be to the start of the correct node.
                    // It may also be attached to the end.
                    cerr << "graph invalid: edge " << e->from() << "->" << e->to()
                         << " doesn't attach to start of " << start_and_edges.first << endl;
                    return false;
                }
                if (!has_node(e->from())) {
                    cerr << "graph invalid: edge from a non-existent node " << e->from() << "->" << e->to() << endl;
                    return false;
                }
                if (!has_node(e->to())) {
                    cerr << "graph invalid: edge to a non-existent node " << e->from() << "->" << e->to() << endl;
                    return false;
                }
            }
        }

        for (pair<const id_t, vector<pair<id_t, bool>>>& end_and_edges : edges_on_end) {
            for (auto& edge_destination : end_and_edges.second) {
                Edge* e = get_edge(NodeSide::pair_from_end_edge(end_and_edges.first, edge_destination));
                if (!e) {
                    cerr << "graph invalid, edge is null" << endl;
                    return false;
                }
                if(end_and_edges.first != e->to() && end_and_edges.first != e->from()) {
                    // It needs to be attached to the node we looked up
                    cerr << "graph invalid: edge " << e->from() << "->" << e->to()
                         << " doesn't have end-indexed node in " << end_and_edges.first << "<->"
                         << edge_destination.first << endl;
                    return false;
                }
                if(edge_destination.first != e->to() && edge_destination.first != e->from()) {
                    // It also needs to be attached to the node it says it goes to
                    cerr << "graph invalid: edge " << e->from() << "->" << e->to()
                         << " doesn't have non-end-indexed node in " << end_and_edges.first << "<->"
                         << edge_destination.first << endl;
                    return false;
                }
                if(!((end_and_edges.first == e->to() && e->to_end()) ||
                     (end_and_edges.first == e->from() && !e->from_start()))) {

                    // The edge needs to actually attach to the end of the node we looked it up for.
                    // So at least one of its ends has to be to the end of the correct node.
                    // It may also be attached to the start.
                    cerr << "graph invalid: edge " << e->from() << "->" << e->to()
                         << " doesn't attach to end of " << end_and_edges.first << endl;
                    return false;
                }
                if (!has_node(e->from())) {
                    cerr << "graph invalid: edge from a non-existent node " << e->from() << "->" << e->to() << endl;
                    return false;
                }
                if (!has_node(e->to())) {
                    cerr << "graph invalid: edge to a non-existent node " << e->from() << "->" << e->to() << endl;
                    return false;
                }
            }
        }
    }

    if (check_paths) {
        bool paths_ok = true;
        function<void(const Path&)> lambda =
            [this, &paths_ok]
            (const Path& path) {
            if (!paths_ok) {
                return;
            }
            if (!path.mapping_size()) {
                cerr << "graph invalid: path " << path.name() << " has no component mappings" << endl;
                paths_ok = false;
                return;
            }
            if (path.mapping_size() == 1) {
                // handle the single-entry case
                if (!path.mapping(0).has_position()) {
                    cerr << "graph path " << path.name() << " has no position in mapping "
                         << pb2json(path.mapping(0)) << endl;
                    paths_ok = false;
                    return;
                }
            }

            for (size_t i = 1; i < path.mapping_size(); ++i) {
                auto& m1 = path.mapping(i-1);
                auto& m2 = path.mapping(i);
                if (!m1.has_position()) {
                    cerr << "graph path " << path.name() << " has no position in mapping "
                         << pb2json(m1) << endl;
                    paths_ok = false;
                    return;
                }
                if (!m2.has_position()) {
                    cerr << "graph path " << path.name() << " has no position in mapping "
                         << pb2json(m2) << endl;
                    paths_ok = false;
                    return;
                }
                if (!adjacent_mappings(m1, m2)) continue; // the path is completely represented here
                auto s1 = NodeSide(m1.position().node_id(), (m1.position().is_reverse() ? false : true));
                auto s2 = NodeSide(m2.position().node_id(), (m2.position().is_reverse() ? true : false));
                // check that we always have an edge between the two nodes in the correct direction
                if (!has_edge(s1, s2)) {
                    cerr << "graph path '" << path.name() << "' invalid: edge from "
                         << s1 << " to " << s2 << " does not exist" << endl;
                    paths_ok = false;
                    //return;;
                }

                // in the four cases below, we check that edges always incident the tips of nodes
                // when edit length, offsets and strand flipping of mappings are taken into account:

                // NOTE: Because of the !adjacent_mappings check above, mappings that are out of order
                //       will be ignored.  If they are invalid, it won't be caught.  Solution is
                //       to sort by rank, but I'm not sure if any of this is by design or not...

                auto& p1 = m1.position();
                auto& n1 = *get_node(p1.node_id());
                auto& p2 = m2.position();
                auto& n2 = *get_node(p2.node_id());
                // count up how many bases of the node m1 covers.
                id_t m1_edit_length = m1.edit_size() == 0 ? n1.sequence().length() : 0;
                for (size_t edit_idx = 0; edit_idx < m1.edit_size(); ++edit_idx) {
                    m1_edit_length += m1.edit(edit_idx).from_length();
                }

                // verify that m1 ends at offset length-1 for forward mapping
                if (p1.offset() + m1_edit_length != n1.sequence().length()) {
                    cerr << "graph path '" << path.name() << "' has invalid mapping " << pb2json(m1)
                    << ": offset (" << p1.offset() << ") + from_length (" << m1_edit_length << ")"
                    << " != node length (" << n1.sequence().length() << ")" << endl;
                    paths_ok = false;
                    return;
                }
                // verify that m2 starts at offset 0 for forward mapping
                if (p2.offset() > 0) {
                    cerr << "graph path '" << path.name() << "' has invalid mapping " << pb2json(m2)
                    << ": offset=" << p2.offset() << " found when offset=0 expected" << endl;
                    paths_ok = false;
                        return;
                }
            }

            // check that the mappings have the right length
            for (size_t i = 0; i < path.mapping_size(); ++i) {
                auto& m = path.mapping(i);
                // get the node
                auto n = get_node(m.position().node_id());
                if (mapping_from_length(m) + m.position().offset() > n->sequence().size()) {
                    cerr << "graph path " << path.name() << " has a mapping which "
                         << "matches sequence outside of the node it maps to "
                         << pb2json(m)
                         << " vs "
                         << pb2json(*n) << endl;
                    paths_ok = false;
                    return;
                }
            }

            // check that the mappings all match the graph
            /*
            for (size_t i = 0; i < path.mapping_size(); ++i) {
                auto& m = path.mapping(i);
                if (!mapping_is_total_match(m)) {
                    cerr << "graph path " << path.name() << " has an imperfect mapping "
                         << pb2json(m) << endl;
                    paths_ok = false;
                    return;
                }
            }
            */

        };
        paths.for_each(lambda);
        if (!paths_ok) return false;
    }

    return true;
}

void VG::to_dot(ostream& out, vector<Alignment> alignments,
                bool show_paths,
                bool walk_paths,
                bool annotate_paths,
                bool show_mappings,
                bool simple_mode,
                bool invert_edge_ports,
                bool color_variants,
                int random_seed) {
    // setup graphviz output
    out << "digraph graphname {" << endl;
    out << "    node [shape=plaintext];" << endl;
    out << "    rankdir=LR;" << endl;
    //out << "    fontsize=22;" << endl;
    //out << "    colorscheme=paired12;" << endl;
    //out << "    splines=line;" << endl;
    //out << "    splines=true;" << endl;
    //out << "    smoothType=spring;" << endl;
    for (int i = 0; i < graph.node_size(); ++i) {
        Node* n = graph.mutable_node(i);
        auto node_paths = paths.of_node(n->id());
        if (!simple_mode) {
            out << "    " << n->id() << " [label=\"" << n->id() << ":" << n->sequence() << "\",shape=box,penwidth=2,";
        } else {
            out << "    " << n->id() << " [label=\"" << n->id() << "\",penwidth=2,shape=circle,";
        }
        // for neato output, which tends to randomly order the graph
        if (!simple_mode) {
            if (is_head_node(n)) {
                out << "pos=\"" << -graph.node_size()*100 << ", "<< -10 << "\"";
            } else if (is_tail_node(n)) {
                out << "pos=\"" << graph.node_size()*100 << ", "<< -10 << "\"";
            }
        }
        if (color_variants && node_paths.size() == 0){
            out << "color=red,";
        }
        out << "];" << endl;

  }

    // We're going to fill this in with all the path (symbol, color) label
    // pairs that each edge should get, by edge pointer. If a path takes an
    // edge multiple times, it will appear only once.
    map<Edge*, set<pair<string, string>>> symbols_for_edge;

    if(annotate_paths) {
        // We're going to annotate the paths, so we need to give them symbols and colors.
        Pictographs picts(random_seed);
        Colors colors(random_seed);
        // Work out what path symbols belong on what edges
        function<void(const Path&)> lambda = [this, &picts, &colors, &symbols_for_edge](const Path& path) {
            // Make up the path's label
            string path_label = picts.hashed(path.name());
            string color = colors.hashed(path.name());
            for (int i = 0; i < path.mapping_size(); ++i) {
                const Mapping& m1 = path.mapping(i);
                if (i < path.mapping_size()-1) {
                    const Mapping& m2 = path.mapping(i+1);
                    // skip if they are not contiguous
                    if (!adjacent_mappings(m1, m2)) continue;
                    // Find the Edge connecting the mappings in the order they occur in the path.
                    Edge* edge_used = get_edge(NodeTraversal(get_node(m1.position().node_id()), m1.position().is_reverse()),
                                               NodeTraversal(get_node(m2.position().node_id()), m2.position().is_reverse()));

                    // Say that edge should have this symbol
                    symbols_for_edge[edge_used].insert(make_pair(path_label, color));
                }
            }
        };
        paths.for_each(lambda);
    }

    for (int i = 0; i < graph.edge_size(); ++i) {
        Edge* e = graph.mutable_edge(i);
        auto from_paths = map_keys_to_set(paths.of_node(e->from()));
        auto to_paths = map_keys_to_set(paths.of_node(e->to()));
        set<string> both_paths;
        std::set_intersection(from_paths.begin(), from_paths.end(),
                              to_paths.begin(), to_paths.end(),
                              std::inserter(both_paths, both_paths.begin()));
        // are both nodes in the same path?
        /*
        bool in_path = !(!paths.empty()
                        && (both_paths.empty()
                            || !paths.are_consecutive_nodes_in_path(e->from(), e->to(),
                                                                    *both_paths.begin())));
                                                                    */

        // Grab the annotation symbols for this edge.
        auto annotations = symbols_for_edge.find(e);

        // Is the edge in the "wrong" direction for rank constraints?
        bool is_backward = e->from_start() && e->to_end();

        if(is_backward) {
            // Flip the edge around and write it forward.
            Edge* original = e;
            e = new Edge();

            e->set_from(original->to());
            e->set_from_start(!original->to_end());

            e->set_to(original->from());
            e->set_to_end(!original->from_start());
        }

        // display what kind of edge we have using different edge head and tail styles
        // depending on if the edge comes from the start or not
        if (!simple_mode) {
            out << "    " << e->from() << " -> " << e->to();
            out << " [dir=both,";
            if (!invert_edge_ports && e->from_start()
                || invert_edge_ports && !e->from_start()) {
                out << "arrowtail=none,";
                out << "tailport=sw,";
            } else {
                out << "arrowtail=none,";
                out << "tailport=ne,";
            }
            if (!invert_edge_ports && e->to_end()
                || invert_edge_ports && !e->to_end()) {
                out << "arrowhead=none,";
                out << "headport=se";
            } else {
                out << "arrowhead=none,";
                out << "headport=nw";
            }
            out << ",penwidth=2";

            if(annotations != symbols_for_edge.end()) {
                // We need to put a label on the edge with all the colored
                // characters for paths using it.
                out << ",label=<";

                for(auto& string_and_color : (*annotations).second) {
                    // Put every symbol in its font tag.
                    out << "<FONT COLOR=\"" << string_and_color.second << "\">" << string_and_color.first << "</FONT>";
                }

                out << ">";
            }
            out << "];" << endl;
        } else {
            out << "    " << e->from() << " -> " << e->to() << endl;
        }

        if(is_backward) {
            // We don't need this duplicate edge
            delete e;
        }
    }

    // add nodes for the alignments and link them to the nodes they match
    int alnid = max_node_id()+1;
    for (auto& aln : alignments) {
        // check direction
        if (!aln.has_path()) continue; // skip pathless alignments
        if (!aln.path().mapping(0).position().is_reverse()) {
            out << "    " << alnid << " [label=\"" << aln.name() << "\"];" << endl;
            out << "    " << alnid << " -> " << alnid+1 << " [dir=none,color=black];" << endl;
            alnid++;
            out << "    " << alnid << " [label=\"+\",fontcolor=green];" << endl;
            out << "    " << alnid << " -> " << alnid+1 << " [dir=none,color=green];" << endl;
        } else {
            out << "    " << alnid << " [label=\"-\",fontcolor=purple];" << endl;
            out << "    " << alnid << " -> " << alnid+1 << " [dir=none,color=purple];" << endl;
        }
        alnid++;
        for (int i = 0; i < aln.path().mapping_size(); ++i) {
            const Mapping& m = aln.path().mapping(i);
            //void mapping_cigar(const Mapping& mapping, vector<pair<int, char> >& cigar);
            //string cigar_string(vector<pair<int, char> >& cigar);
            //mapid << alnid << ":" << m.position().node_id() << ":" << cigar_string(cigar);
            //mapid << cigar_string(cigar) << ":"
            //      << (m.is_reverse() ? "-" : "+") << m.position().offset() << ":"
            string mstr;
            if (!simple_mode) {
                stringstream mapid;
                mapid << pb2json(m);
                mstr = mapid.str();
                mstr.erase(std::remove_if(mstr.begin(), mstr.end(),
                                          [](char c) { return c == '"'; }), mstr.end());
                mstr = wrap_text(mstr, 50);
            }
            // determine sequence of this portion of the alignment
            // set color based on cigar/mapping relationship
            // some mismatch, indicate with orange color
            string color;
            if (!simple_mode) {
                color = mapping_is_simple_match(m) ? "blue" : "orange";
            } else {
                color = "/rdylgn11/" + convert(round((1-divergence(m))*10)+1);
            }

            if (simple_mode) {
                out << "    "
                    << alnid << " [label=\""
                    << m.position().node_id() << "\"" << "shape=circle," //penwidth=2,"
                    << "style=filled,"
                    << "fillcolor=\"" << color << "\","
                    << "color=\"" << color << "\"];" << endl;

            } else {
                out << "    "
                    << alnid << " [label=\""
                    << mstr << "\",fontcolor=" << color << ",fontsize=10];" << endl;
            }
            if (i > 0) {
                out << "    "
                    << alnid-1 << " -> "
                    << alnid << "[dir=none,color=\"black\"];" << endl;
            }
            out << "    "
                << alnid << " -> " << m.position().node_id()
                << "[dir=none,style=invis];" << endl;
            out << "    { rank = same; " << alnid << "; " << m.position().node_id() << "; };" << endl;
            //out << "    " << m.position().node_id() << " -- " << alnid << "[color=" << color << ", style=invis];" << endl;
            alnid++;
        }
        if (!aln.path().mapping(aln.path().mapping_size()-1).position().is_reverse()) {
            out << "    " << alnid << " [label=\"-""\",fontcolor=purple];" << endl;
            out << "    " << alnid-1 << " -> " << alnid << " [dir=none,color=purple];" << endl;
        } else {
            out << "    " << alnid << " [label=\"+""\",fontcolor=green];" << endl;
            out << "    " << alnid-1 << " -> " << alnid << " [dir=none,color=green];" << endl;
            alnid++;
            out << "    " << alnid << " [label=\"" << aln.name() << "\"];" << endl;
            out << "    " << alnid << " -> " << alnid-1 << " [dir=none,color=black];" << endl;
        }
        alnid++;
    }

    // include paths
    if (show_paths || walk_paths) {
        int pathid = alnid;
        Pictographs picts(random_seed);
        Colors colors(random_seed);
        function<void(const Path&)> lambda =
            [this,&pathid,&out,&picts,&colors,show_paths,walk_paths,show_mappings]
            (const Path& path) {
            string path_label = picts.hashed(path.name());
            string color = colors.hashed(path.name());
            if (show_paths) {
                for (int i = 0; i < path.mapping_size(); ++i) {
                    const Mapping& m = path.mapping(i);
                    stringstream mapid;
                    mapid << path_label << " " << m.position().node_id();
                    stringstream mappings;
                    if (show_mappings) {
                        mappings << pb2json(m);
                    }
                    string mstr = mappings.str();
                    mstr.erase(std::remove_if(mstr.begin(), mstr.end(), [](char c) { return c == '"'; }), mstr.end());
                    mstr = wrap_text(mstr, 50);

                    if (i == 0) { // add the path name at the start
                        out << "    " << pathid << " [label=\"" << path_label << " "
                            << path.name() << "  " << m.position().node_id() << " "
                            << mstr << "\",fontcolor=\"" << color << "\"];" << endl;
                    } else {
                        out << "    " << pathid << " [label=\"" << mapid.str() << " "
                            << mstr
                            << "\",fontcolor=\"" << color << "\"];" << endl;
                    }
                    if (i > 0 && adjacent_mappings(path.mapping(i-1), m)) {
                        out << "    " << pathid-1 << " -> " << pathid << " [dir=none,color=\"" << color << "\"];" << endl;
                    }
                    out << "    " << pathid << " -> " << m.position().node_id() << " [dir=none,color=\"" << color << "\", style=invis];" << endl;
                    //out << "    " << pathid << " -> " << m.position().node_id() << "[headport=n,tailport=s,color=\"" << color << "\"];" << endl;
                    out << "    { rank = same; " << pathid << "; " << m.position().node_id() << "; };" << endl;
                    pathid++;
                }
            }
            if (walk_paths) {
                for (int i = 0; i < path.mapping_size(); ++i) {
                    const Mapping& m1 = path.mapping(i);
                    if (i < path.mapping_size()-1) {
                        const Mapping& m2 = path.mapping(i+1);
                        out << m1.position().node_id() << " -> " << m2.position().node_id() << " [dir=none,tailport=ne,headport=nw,color=\"" << color << "\",label=\"     " << path_label << "     \",fontcolor=\"" << color << "\"];" << endl;
                    }
                }
            }
        };
        paths.for_each(lambda);
    }

    out << "}" << endl;

}


void VG::to_gfa(ostream& out) {
  GFAKluge gg;
  gg.set_version();

    // TODO moving to GFAKluge
    // problem: protobuf longs don't easily go to strings....
    for (int i = 0; i < graph.node_size(); ++i){
        sequence_elem s_elem;
        Node* n = graph.mutable_node(i);
        // Fill seq element for a node
        s_elem.name = to_string(n->id());
        s_elem.sequence = n->sequence();
        gg.add_sequence(s_elem);

            auto& node_mapping = paths.get_node_mapping(n->id());
            set<Mapping*> seen;
            for (auto& p : node_mapping) {
                for (auto* m : p.second) {
                    if (seen.count(m)) continue;
                    else seen.insert(m);
                    const Mapping& mapping = *m;
                    string cigar;
                    if (mapping.edit_size() > 0) {
                        vector<pair<int, char> > cigarv;
                        mapping_cigar(mapping, cigarv);
                        cigar = cigar_string(cigarv);
                    } else {
                        // empty mapping edit implies perfect match
                        stringstream cigarss;
                        cigarss << n->sequence().size() << "M";
                        cigar = cigarss.str();
                    }
                    bool orientation = mapping.position().is_reverse();
                    path_elem p_elem;
                    p_elem.name = p.first;
                    p_elem.source_name = to_string(n->id());
                    p_elem.rank = mapping.rank();
                    p_elem.is_reverse = orientation;
                    p_elem.cigar = cigar;

                    gg.add_path(p_elem.source_name, p_elem);
                }
            }

    }

    for (int i = 0; i < graph.edge_size(); ++i) {
        Edge* e = graph.mutable_edge(i);
        link_elem l;
        l.source_name = to_string(e->from());
        l.sink_name = to_string(e->to());
        l.source_orientation_forward = ! e->from_start();
        l.sink_orientation_forward =  ! e->to_end();
        l.cigar = "0M";
        gg.add_link(l.source_name, l);
    }
    out << gg;


    // map<id_t, vector<string> > sorted_output;
    // out << "H" << "\t" << "HVN:Z:1.0" << endl;
    // for (int i = 0; i < graph.node_size(); ++i) {
    //     Node* n = graph.mutable_node(i);
    //     stringstream s;
    //     s << "S" << "\t" << n->id() << "\t" << n->sequence() << "\n";
    //     auto& node_mapping = paths.get_node_mapping(n->id());
    //     set<Mapping*> seen;
    //     for (auto& p : node_mapping) {
    //         for (auto* m : p.second) {
    //             if (seen.count(m)) continue;
    //             else seen.insert(m);
    //             const Mapping& mapping = *m;
    //             string cigar;
    //             if (mapping.edit_size() > 0) {
    //                 vector<pair<int, char> > cigarv;
    //                 mapping_cigar(mapping, cigarv);
    //                 cigar = cigar_string(cigarv);
    //             } else {
    //                 // empty mapping edit implies perfect match
    //                 stringstream cigarss;
    //                 cigarss << n->sequence().size() << "M";
    //                 cigar = cigarss.str();
    //             }
    //             string orientation = mapping.position().is_reverse() ? "-" : "+";
    //             s << "P" << "\t" << n->id() << "\t" << p.first << "\t"
    //               << mapping.rank() << "\t" << orientation << "\t" << cigar << "\n";
    //         }
    //     }
    //     sorted_output[n->id()].push_back(s.str());
    // }
    // for (int i = 0; i < graph.edge_size(); ++i) {
    //     Edge* e = graph.mutable_edge(i);
    //     stringstream s;
    //     s << "L" << "\t" << e->from() << "\t"
    //       << (e->from_start() ? "-" : "+") << "\t"
    //       << e->to() << "\t"
    //       << (e->to_end() ? "-" : "+") << "\t"
    //       << "0M" << endl;
    //     sorted_output[e->from()].push_back(s.str());
    // }
    // for (auto& chunk : sorted_output) {
    //     for (auto& line : chunk.second) {
    //         out << line;
    //     }
    // }
}

void VG::to_turtle(ostream& out, const string& rdf_base_uri) {
    out << "@base <http://example.org/vg/> . " << endl;
    out << "@prefix node: <" <<  rdf_base_uri <<"node/> . " << endl;
    out << "@prefix path: <" <<  rdf_base_uri <<"path/> . " << endl;
    out << "@prefix step: <" <<  rdf_base_uri <<"step/> . " << endl;
    out << "@prefix rdf: <http://www.w3.org/1999/02/22-rdf-syntax-ns#> . " << endl;
    //Ensure that mappings are sorted by ranks
    paths.sort_by_mapping_rank();
    for (int i = 0; i < graph.node_size(); ++i) {
        Node* n = graph.mutable_node(i);
        out << "node:" << n->id() << " rdf:value \"" << n->sequence() << "\" . " << endl ;

        auto& node_mapping = paths.get_node_mapping(n->id());
	    set<Mapping*> seen;
        for (auto& p : node_mapping) {
            for (auto* m : p.second) {
                if (seen.count(m)) continue;
                else seen.insert(m);
                const Mapping& mapping = *m;
                out << "s:" << p.first << "#" << mapping.rank() << "a <Step> ;" << endl ;
                out << " <rank> " << mapping.rank() << " ; "  << endl ;
                string orientation = mapping.position().is_reverse() ? "<reverseOfNode>" : "<node>";
                out << "\t" << orientation <<" n:" << n->id() << " ; " << endl;
                out << "\t<path> p:" << p.first << " . " << endl;
            }
        }
    }
    function<void(const Path&)> lambda = [&out]
        (const Path& path) {
            uint64_t offset=0; //We could have more than 2gigabases in a path
            for (auto &m : path.mapping()) {
                out << "step:" << path.name() << "#" << m.rank() << " <position> "<< offset<<" . " << endl;
                offset += mapping_to_length(m);
            }
        };
    paths.for_each(lambda);
    for (int i = 0; i < graph.edge_size(); ++i) {
        Edge* e = graph.mutable_edge(i);
        out << "node:" << e->from();
        if (e->from_start() && e->to_end()) {
          out << " <linksReverseToReverse> " ; // <--
        } else if (e->from_start() && !e->to_end()) {
          out << " <linksReverseToForward> " ; // -+
        } else if (e->to_end()) {
          out << " <linksForwardToReverse> " ; //+-
        } else {
          out << " <linksForwardToForward> " ; //++
        }
        out << "node:" << e->to() << " . " << endl;
    }
}
void VG::destroy_alignable_graph(void) {
    if (gssw_aligner != NULL) {
        delete gssw_aligner;
    }
}

void VG::connect_node_to_nodes(Node* node, vector<Node*>& nodes, bool from_start) {
    for (vector<Node*>::iterator n = nodes.begin(); n != nodes.end(); ++n) {
        // Connect them left to right, unless instructed otherwise
        create_edge(node, (*n), from_start, false);
    }
}

void VG::connect_nodes_to_node(vector<Node*>& nodes, Node* node, bool to_end) {
    for (vector<Node*>::iterator n = nodes.begin(); n != nodes.end(); ++n) {
        // Connect them left to right, unless instructed otherwise
        create_edge((*n), node, false, to_end);
    }
}

void VG::connect_node_to_nodes(NodeTraversal node, vector<NodeTraversal>& nodes) {
    for (vector<NodeTraversal>::iterator n = nodes.begin(); n != nodes.end(); ++n) {
        // Connect them left to right
        create_edge(node, (*n));
    }
}

void VG::connect_nodes_to_node(vector<NodeTraversal>& nodes, NodeTraversal node) {
    for (vector<NodeTraversal>::iterator n = nodes.begin(); n != nodes.end(); ++n) {
        // Connect them left to right
        create_edge((*n), node);
    }
}

// join all subgraphs together to a "null" head node
Node* VG::join_heads(void) {
    // Find the head nodes
    vector<Node*> heads;
    head_nodes(heads);

    // Then create the new node (so it isn't picked up as a head)
    current_id = max_node_id()+1;
    Node* root = create_node("N");

    // Wire it to all the heads and return
    connect_node_to_nodes(root, heads);
    return root;
}

void VG::join_heads(Node* node, bool from_start) {
    vector<Node*> heads;
    head_nodes(heads);

    // If the node we have been given shows up as a head, remove it.
    for(auto i = heads.begin(); i != heads.end(); ++i) {
        if(*i == node) {
            heads.erase(i);
            break;
        }
    }

    connect_node_to_nodes(node, heads, from_start);
}

void VG::join_tails(Node* node, bool to_end) {
    vector<Node*> tails;
    tail_nodes(tails);

    // If the node we have been given shows up as a tail, remove it.
    for(auto i = tails.begin(); i != tails.end(); ++i) {
        if(*i == node) {
            tails.erase(i);
            break;
        }
    }

    connect_nodes_to_node(tails, node, to_end);
}

void VG::add_start_end_markers(int length,
                               char start_char, char end_char,
                               Node*& start_node, Node*& end_node,
                               id_t start_id, id_t end_id) {
    // This set will hold all the nodes we haven't attached yet. But we don't
    // want it to hold the head_tail_node, so we fill it in now.
    set<Node*> unattached;
    for_each_node([&](Node* node) {
        unattached.insert(node);
    });

    // We handle the head and tail joining ourselves so we can do connected components.
    // We collect these before we add the new head/tail node so we don't have to filter it out later.
    vector<Node*> heads;
    head_nodes(heads);
    vector<Node*> tails;
    tail_nodes(tails);

    if(start_node == nullptr) {
        // We get to create the node. In its forward orientation it's the start node, so we use the start character.
        string start_string(length, start_char);
        start_node = create_node(start_string, start_id);
    } else {
        // We got a node to use
        add_node(*start_node);
    }

    if(end_node == nullptr) {
        // We get to create the node. In its forward orientation it's the end node, so we use the end character.
        string end_string(length, end_char);
        end_node = create_node(end_string, end_id);
    } else {
        // We got a node to use
        add_node(*end_node);
    }

#ifdef debug
    cerr << "Start node is " << start_node->id() << ", end node is " << end_node->id() << endl;
#endif

    for(Node* head : heads) {
        if(unattached.count(head)) {
            // This head is unconnected.

            // Mark everything it's attached to as attached
            for_each_connected_node(head, [&](Node* node) {
                unattached.erase(node);
            });
        }

        // Tie it to the start node
        create_edge(start_node, head);
#ifdef debug
    cerr << "Added edge " << start_node->id() << "->" << head->id() << endl;
#endif
    }

    for(Node* tail : tails) {
        if(unattached.count(tail)) {
            // This tail is unconnected.

            // Mark everything it's attached to as attached
            for_each_connected_node(tail, [&](Node* node) {
                unattached.erase(node);
            });
        }

        // Tie it to the end node
        create_edge(tail, end_node);
#ifdef debug
    cerr << "Added edge " << tail->id() << "->" << end_node->id() << endl;
#endif
    }

    // Find the connected components that aren't attached, if any.
    while(!unattached.empty()) {
        // Grab and attach some unattached node
        Node* to_attach = *(unattached.begin());

        // Mark everything it's attached to as attached
        for_each_connected_node(to_attach, [&](Node* node) {
            unattached.erase(node);
        });

        // Add the edge
        create_edge(start_node, to_attach);
#ifdef debug
        cerr << "Added cycle-breaking edge " << start_node->id() << "->" << to_attach->id() << endl;
#endif
        vector<Edge*> edges;
        edges_of_node(to_attach, edges);
        for (auto edge : edges) {
            //cerr << "edge of " << to_attach->id() << " " << edge->from() << " " << edge->to() << endl;
            if (edge->to() == to_attach->id() && edge->from() != start_node->id()) {
                //cerr << "creating edge" << endl;
                Edge* created = create_edge(edge->from(), end_node->id(), edge->from_start(), false);
#ifdef debug
                cerr << "Added edge " << pb2json(*created) << " in response to " << pb2json(*edge) << endl;
#endif
            }
        }
#ifdef debug
        cerr << "Broke into disconnected component at " << to_attach->id() << endl;
#endif
    }

    // Now we have no more disconnected stuff in our graph.

#ifdef debug
    cerr << "Start node edges: " << endl;
    vector<Edge*> edges;
    edges_of_node(start_node, edges);
    for(auto e : edges) {
        std::cerr << pb2json(*e) << std::endl;
    }

    cerr << "End node edges: " << endl;
    edges.clear();
    edges_of_node(end_node, edges);
    for(auto e : edges) {
        std::cerr << pb2json(*e) << std::endl;
    }

#endif

}

map<id_t, pair<id_t, bool> > VG::overlay_node_translations(const map<id_t, pair<id_t, bool> >& over,
                                                           const map<id_t, pair<id_t, bool> >& under) {
    map<id_t, pair<id_t, bool> > overlay = under;
    // for each over, check if we should map to the under
    // if so, adjust
    for (auto& o : over) {
        id_t new_id = o.first;
        id_t old_id = o.second.first;
        bool is_rev = o.second.second;
        auto u = under.find(old_id);
        if (u != under.end()) {
            id_t oldest_id = u->second.first;
            bool was_rev = u->second.second;
            overlay[new_id] = make_pair(oldest_id,
                                        is_rev ^ was_rev);
            /*
            cerr << "double trans "
                 << new_id << " -> " << old_id
                 << " -> " << oldest_id << endl;
            */
        } else {
            overlay[o.first] = o.second;
        }
    }
    /*
    for (auto& o : overlay) {
        cerr << o.first << " -> " << o.second.first
             << (o.second.second?"-":"+") << endl;
    }
    */
    return overlay;
}

Alignment VG::align(const Alignment& alignment) {

    //cerr << "aligning " << pb2json(alignment) << endl;
    // to be completely aligned, the graph's head nodes need to be fully-connected to a common root
    auto aln = alignment;

    map<id_t, pair<id_t, bool> > unfold_trans;
    map<id_t, pair<id_t, bool> > dagify_trans;
    uint32_t max_length = alignment.sequence().size();

    VG dag = unfold(max_length, unfold_trans).dagify(max_length, dagify_trans, max_length);
    // overlay the translations
    auto trans = overlay_node_translations(dagify_trans, unfold_trans);

    // Join to a common root, so alignment covers the entire graph
    // Put the nodes in sort order within the graph
    // and break any remaining cycles
    Node* root = dag.join_heads();
    dag.sort();

    gssw_aligner = new GSSWAligner(dag.graph);
    gssw_aligner->align(aln);

    delete gssw_aligner;
    gssw_aligner = NULL;

    /*
    auto check_aln = [&](VG& graph, const Alignment& a) {
        cerr << "checking alignment" << endl;
        if (a.has_path()) {
            auto seq = graph.path_string(a.path());
            //if (aln.sequence().find('N') == string::npos && seq != aln.sequence()) {
            if (seq != a.sequence()) {
                cerr << "alignment does not match graph " << endl
                     << pb2json(a) << endl
                     << "expect:\t" << a.sequence() << endl
                     << "got:\t" << seq << endl;
                write_alignment_to_file(a, "fail.gam");
                graph.serialize_to_file("fail.vg");
                assert(false);
            }
        }
    };
    */
    //check_aln(dag, aln);
    translate_nodes(aln, trans, [&](id_t node_id) {
            // We need to feed in the lengths of nodes, so the offsets in the alignment can be updated.
            return get_node(node_id)->sequence().size();
        });
    //check_aln(*this, aln);

    return aln;
}

Alignment VG::align(const string& sequence) {
    Alignment alignment;
    alignment.set_sequence(sequence);
    return align(alignment);
}

const string VG::hash(void) {
    stringstream s;
    serialize_to_ostream(s);
    return sha1sum(s.str());
}

void VG::for_each_kmer_parallel(int kmer_size,
                                int edge_max,
                                function<void(string&, list<NodeTraversal>::iterator, int, list<NodeTraversal>&, VG&)> lambda,
                                int stride,
                                bool allow_dups,
                                bool allow_negatives) {
    _for_each_kmer(kmer_size, edge_max, lambda, true, stride, allow_dups, allow_negatives);
}

void VG::for_each_kmer(int kmer_size,
                       int edge_max,
                       function<void(string&, list<NodeTraversal>::iterator, int, list<NodeTraversal>&, VG&)> lambda,
                       int stride,
                       bool allow_dups,
                       bool allow_negatives) {
    _for_each_kmer(kmer_size, edge_max, lambda, false, stride, allow_dups, allow_negatives);
}

void VG::for_each_kmer_of_node(Node* node,
                               int kmer_size,
                               int edge_max,
                               function<void(string&, list<NodeTraversal>::iterator, int, list<NodeTraversal>&, VG&)> lambda,
                               int stride,
                               bool allow_dups,
                               bool allow_negatives) {
    _for_each_kmer(kmer_size, edge_max, lambda, false, stride, allow_dups, allow_negatives, node);
}

void VG::_for_each_kmer(int kmer_size,
                        int edge_max,
                        function<void(string&, list<NodeTraversal>::iterator, int, list<NodeTraversal>&, VG&)> lambda,
                        bool parallel,
                        int stride,
                        bool allow_dups,
                        bool allow_negatives,
                        Node* node) {

#ifdef debug
    cerr << "Looking for kmers of size " << kmer_size << " over " << edge_max << " edges with node " << node << endl;
#endif

    // use an LRU cache to clean up duplicates over the last 1mb
    // use one per thread so as to avoid contention
    // If we aren't starting a parallel kmer iteration from here, just fill in 0.
    // TODO: How do we know this is big enough?
    map<int, LRUCache<string, bool>* > lru;
#pragma omp parallel
    {
#pragma omp single
        for (int i = 0; i < (parallel ? omp_get_num_threads() : 1); ++i) {
            lru[i] = new LRUCache<string, bool>(100000);
        }
    }
    // constructs the cache key
    // experiment -- use a struct here
    // We deduplicate kmers based on where they start, where they end (optionally), and where they are viewed from.
    auto make_cache_key = [](string& kmer,
                             id_t start_node, int start_pos,
                             id_t view_node, int view_pos,
                             id_t end_node, int end_pos) -> string {
        string cache_key = kmer;
        cache_key.resize(kmer.size() + sizeof(Node*) + 3*sizeof(id_t) + 3*sizeof(int));
        memcpy((char*)cache_key.c_str()+kmer.size(), &start_node, sizeof(id_t));
        memcpy((char*)cache_key.c_str()+kmer.size()+1*sizeof(id_t), &start_pos, sizeof(int));
        memcpy((char*)cache_key.c_str()+kmer.size()+1*sizeof(id_t)+1*sizeof(int), &view_node, sizeof(id_t));
        memcpy((char*)cache_key.c_str()+kmer.size()+2*sizeof(id_t)+1*sizeof(int), &view_pos, sizeof(int));
        memcpy((char*)cache_key.c_str()+kmer.size()+2*sizeof(id_t)+2*sizeof(int), &end_node, sizeof(id_t));
        memcpy((char*)cache_key.c_str()+kmer.size()+3*sizeof(id_t)+2*sizeof(int), &end_pos, sizeof(int));
        return cache_key;
    };

    auto handle_path = [this,
                        &lambda,
                        kmer_size,
                        stride,
                        allow_dups,
                        allow_negatives,
                        &lru,
                        &make_cache_key,
                        &parallel,
                        &node](list<NodeTraversal>::iterator forward_node, list<NodeTraversal>& forward_path) {
#ifdef debug
        cerr << "Handling path: " << endl;
        for(auto& traversal : forward_path) {
            cerr << "\t" << traversal.node->id() << " " << traversal.backward << endl;
        }
#endif

        // Reserve a place for our reversed path, if we find ourselves needing to flip it around for a negative offset kmer.
        list<NodeTraversal> reversed_path;

        // And one for the reversed version of this NodeTraversal on that path
        list<NodeTraversal>::iterator reversed_node;


        // expand the path into a vector :: 1,1,1,2,2,2,2,3,3 ... etc.
        // this makes it much easier to quickly get all the node matches of each kmer

        // We expand out as iterators in the path list, so we can get the
        // traversal but also distinguish different node instances in the path.
        vector<list<NodeTraversal>::iterator> node_by_path_position;
        expand_path(forward_path, node_by_path_position);

        // Go get the cache for this thread if _for_each_kmer launched threads,
        // and the only one we made (thread 0) if we're running this
        // _for_each_kmer call in a single thread. Remember that _for_each_kmer
        // itself may be called by many MPI threads in parallel.
        auto cache = lru[parallel ? omp_get_thread_num() : 0];
        assert(cache != nullptr);

        map<NodeTraversal*, int> node_start;
        node_starts_in_path(forward_path, node_start);

        // now process the kmers of this sequence
        // by first getting the sequence
        string seq = path_string(forward_path);

        // but bail out if the sequence is shorter than the kmer size
        if (seq.size() < kmer_size) return;

        // and then stepping across the path, finding the kmers, and then implied node overlaps
        for (int i = 0; i <= seq.size() - kmer_size; i+=stride) {

            // get the kmer
            string forward_kmer = seq.substr(i, kmer_size);
            // record when we get a kmer match

            // We'll fill this in if needed.
            string reversed_kmer;

            // execute our callback on each kmer/node/position
            // where node == node
            int j = 0;
            while (j < kmer_size) {
                if (forward_node == node_by_path_position[i+j]) {
                    // At this position along this possible kmer, we're in the
                    // instance of the node we're interested in the kmers of.


                    // Grab the node the kmer started at
                    list<NodeTraversal>::iterator start_node = node_by_path_position[i];
                    // And the one it's going to end at
                    list<NodeTraversal>::iterator end_node = node_by_path_position[i + kmer_size - 1];
                    // Work out how far into its actual starting node this kmer started.
                    size_t start_node_offset = i - node_start[&(*start_node)];

                    if(!allow_negatives && node == nullptr) {
                        // If we do allow negatives, we'll just articulate kmers
                        // from both sides whenever they cross edges. Otherwise,
                        // we only want edge-crossing kmers once, so we should
                        // only announce them to the callback from one of their
                        // ends. We arbitrarily choose the end with the lower
                        // node ID.

                        // We only do this when we aren't getting the kmers of a
                        // specific node.

                        if(forward_node == start_node &&
                           (*start_node).node->id() > (*end_node).node->id()) {
                            // We're on the start, but it's ID is larger than the end's.
                            // Announce the kmer from the end instead.
                            ++j;
                            continue;
                        }

                        if(forward_node == end_node &&
                           (*end_node).node->id() > (*start_node).node->id()) {
                            // We're on the end, but it's ID is larger than the start's.
                            // Announce the kmer from the start instead.
                            ++j;
                            continue;
                        }

                        if((*end_node).node->id() == (*start_node).node->id() &&
                            end_node != start_node &&
                            forward_node == end_node) {

                            // If this kmer starts and ends in different
                            // instances of the same node along the path, only
                            // announce it from the one it starts in. Skip the
                            // one it ends in.
                            ++j;
                            continue;
                        }
                    }

                    // We now know we should announce this kmer from this node.

                    // Work out where this node started along the path
                    int node_position = node_start[&(*forward_node)];
                    // And how far into the node this kmer started
                    int kmer_forward_relative_start = i - node_position;
                    // Negative-offset kmers will be processed, but will be
                    // corrected to the opposite strand if negative offsets are
                    // not allowed.
                    int kmer_reversed_relative_start;
                    // Did we flip?
                    bool reversed = false;
                    if(kmer_forward_relative_start < 0 && !allow_negatives) {
                        // This kmer starts at a negative offset from this node.
                        // We need to announce it with a positive offset.

                        size_t node_length = (*forward_node).node->sequence().size();

                        if(kmer_forward_relative_start + kmer_size > node_length) {
                            // If it doesn't start or end in this node, and is
                            // just passing through, there's no way to
                            // articulate it for this node with a positive
                            // offset, so we just skip to the next character in
                            // the kmer.
                            ++j;
                            continue;
                        }

                        // We know the kmer has its end in this node.

                        // If it ends in this node, we can reverse it and get a
                        // positive offset.

                        if(reversed_kmer.empty()) {
                            // Fill in the reversed kmer the first time we need it.
                            reversed_kmer = reverse_complement(forward_kmer);
                        }

                        if(reversed_path.empty()) {
                            // Only fill in the reversed path the first time we need it.
                            for(NodeTraversal traversal : forward_path) {
                                // Operate on copies here
                                traversal.backward = !traversal.backward;
                                reversed_path.push_front(traversal);
                            }

                            // Fill in the reversed iterator too
                            reversed_node = reversed_path.begin();
                            list<NodeTraversal>::iterator i(forward_node);

                            // If forward_node is the last non-end() item, we
                            // don't want to advance reversed_node at all from
                            // the first item of the reversed path.
                            ++i;

                            while(i != forward_path.end()) {
                                // Walk i towards the end of the forward path,
                                // and reversed_node in from the corresponding
                                // end of the reversed path.
                                ++i;
                                ++reversed_node;
                            }
                        }

                        // Flip the kmer start around to something that will be positive.
                        // We don't need a -1 here.
                        kmer_reversed_relative_start = node_length - (kmer_forward_relative_start + kmer_size);

                        // We flipped.
                        reversed = true;
                    }

                    // Set up some references so we don't need to make local
                    // copies unless we actually did need to reverse the kmer.
                    string& kmer = reversed ? reversed_kmer : forward_kmer;
                    list<NodeTraversal>::iterator& instance = reversed ? reversed_node : forward_node;
                    list<NodeTraversal>& path = reversed ? reversed_path : forward_path;
                    int& kmer_relative_start = reversed ? kmer_reversed_relative_start : kmer_forward_relative_start;

                    // Make sure we aren't disobeying instructions
                    assert(!(kmer_relative_start < 0 && !allow_negatives));

                    // What key in the cache do we say that we processed? We're
                    // going to cache kmers in their forward orientation, even
                    // if they are at negative relative offsets and we aren't
                    // supposed to be using those. Because for_each_kpath only
                    // ever presents a node to this function in its forward
                    // orientation, we'll never have a situation where we should
                    // have done the cache key on the opposite strand.
                    string cache_key;
                    if (allow_dups) {
                        // Duplicate kmers starting at the same place are allowed if the paths go to different places next.
                        // This is deduplicating by node ID so we don't have to worry about instances on the path.
                        // TODO: forward_node won't be passed in as a backward traversal from for_each_kpath, right?

                        // figure out past-the-end-of-the-kmer position and node
                        list<NodeTraversal>::iterator past_end = (i+kmer_size >= node_by_path_position.size())
                            ? path.end()
                            : node_by_path_position[i + kmer_size-1];
                        int node_past_end_position = (past_end == path.end()) ? 0 : i+kmer_size - node_start[&(*past_end)];

#ifdef debug
                        cerr << "Checking for duplicates of " << (*start_node).node->id() << "." << start_node_offset
                             << (reversed?"⍃":"⍄")
                             << "-" << forward_kmer  << "-" << (past_end == path.end() ? 0 : (*past_end).node->id()) << "."
                             << node_past_end_position << " viewed from "
                             << (*forward_node).node->id() << " offset " << kmer_forward_relative_start << endl;
#endif

                        cache_key = make_cache_key(forward_kmer, (*start_node).node->id(), start_node_offset,
                                                                 (*forward_node).node->id(), kmer_forward_relative_start,
                                                                 (past_end == path.end() ? 0 : (*past_end).node->id()),
                                                                 node_past_end_position);
                    } else {
                        // Duplicate kmers starting at the same place aren't allowed, no matter where they go after the end.
                        // This is deduplicating by node ID so we don't have to worry about instances on the path.
                        cache_key = make_cache_key(forward_kmer, (*start_node).node->id(), start_node_offset,
                                                                 (*forward_node).node->id(), kmer_forward_relative_start,
                                                                 0, 0);
                    }

                    // See if this kmer is mentioned in the cache already
                    pair<bool, bool> c = cache->retrieve(cache_key);
                    if (!c.second && (*instance).node != NULL) {
                        // TODO: how could we ever get a null node here?
                        // If not, put it in and run on it.
                        cache->put(cache_key, true);

                        lambda(kmer, instance, kmer_relative_start, path, *this);
                    } else {
#ifdef debug
                        cerr << "Skipped " << kmer << " because it was already done" << endl;
#endif
                    }

                }
                ++j;
            }
        }
    };

    auto noop = [](NodeTraversal) { };

    if(node == nullptr) {
        // Look at all the kpaths
        if (parallel) {
            for_each_kpath_parallel(kmer_size, edge_max, noop, noop, handle_path);
        } else {
            for_each_kpath(kmer_size, edge_max, noop, noop, handle_path);
        }
    } else {
        // Look only at kpaths of the specified node
        for_each_kpath_of_node(node, kmer_size, edge_max, noop, noop, handle_path);
    }

}

int VG::path_edge_count(list<NodeTraversal>& path, int32_t offset, int path_length) {
    int edge_count = 0;
    // starting from offset in the first node
    // how many edges do we cross?

    // This is the remaining path length
    int l = path_length;

    // This is the next node we are looking at.
    list<NodeTraversal>::iterator pitr = path.begin();

    // How many bases of the first node can we use?
    int available_in_first_node = (*pitr).node->sequence().size() - offset;

    if(available_in_first_node >= l) {
        // Cross no edges
        return 0;
    }

    l -= available_in_first_node;
    pitr++;
    while (l > 0) {
        // Now we can ignore node orientation
        ++edge_count;
        l -= (*pitr++).node->sequence().size();
    }
    return edge_count;
}

int VG::path_end_node_offset(list<NodeTraversal>& path, int32_t offset, int path_length) {
    // This is the remaining path length
    int l = path_length;

    // This is the next node we are looking at.
    list<NodeTraversal>::iterator pitr = path.begin();

    // How many bases of the first node can we use?
    int available_in_first_node = (*pitr).node->sequence().size() - offset;

    if(available_in_first_node >= l) {
        // Cross no edges
        return available_in_first_node - l;
    }

    l -= available_in_first_node;
    pitr++;
    while (l > 0) {
        l -= (*pitr++).node->sequence().size();
    }
    // Now back out the last node we just took.
    l += (*--pitr).node->sequence().size();

    // Measure form the far end of the last node.
    l = (*pitr).node->sequence().size() - l - 1;

    return l;
}

const vector<Alignment> VG::paths_as_alignments(void) {
    vector<Alignment> alns;
    paths.for_each([this,&alns](const Path& path) {
            alns.emplace_back();
            auto& aln = alns.back();
            *aln.mutable_path() = path; // copy the path
            // now reconstruct the sequence
            aln.set_sequence(this->path_sequence(path));
            aln.set_name(path.name());
        });
    return alns;
}

const string VG::path_sequence(const Path& path) {
    string seq;
    for (size_t i = 0; i < path.mapping_size(); ++i) {
        auto& m = path.mapping(i);
        seq.append(mapping_sequence(m, *get_node(m.position().node_id())));
    }
    return seq;
}

void VG::kmer_context(string& kmer,
                      int kmer_size,
                      int edge_max,
                      bool forward_only,
                      list<NodeTraversal>& path,
                      list<NodeTraversal>::iterator start_node,
                      int32_t start_offset,
                      list<NodeTraversal>::iterator& end_node,
                      int32_t& end_offset,
                      set<tuple<char, id_t, bool, int32_t>>& prev_positions,
                      set<tuple<char, id_t, bool, int32_t>>& next_positions) {

    // Say we couldn't find an and node. We'll replace this when we do.
    end_node = path.end();

    // Start our iterator at our kmer's starting node instance.
    list<NodeTraversal>::iterator np = start_node;

    if (start_offset == 0) {
        // for each node connected to this one
        // what's its last character?
        // add to prev_chars

        // TODO: do we have to check if we would cross too many edges here too?

        vector<NodeTraversal> prev_nodes;
        nodes_prev(*start_node, prev_nodes);
        for (auto n : prev_nodes) {
            const string& seq = n.node->sequence();
            // We have to find the last chartacter in either orientation.
            char c = n.backward ? reverse_complement(seq[0]) : seq[seq.size()-1];
            // Also note the previous position (which was always the last character in the orientation we'll be looking at it in)
            prev_positions.insert(
                make_tuple(c,
                           n.node->id(),
                           n.backward,
                           n.node->sequence().size() - 1));
        }
    } else {
        // Grab and point to the previous character in this orientation of this node.
        const string& seq = (*start_node).node->sequence();
        // If we're on the reverse strand, we go one character later in the
        // string than start_offset from its end. Otherwise we go one character
        // earlier than start_offset from its beginning.
        // TODO: Add some methods to get characters at offsets in NodeTraversals
        char c = (*start_node).backward ? reverse_complement(seq[seq.size() - start_offset]) : seq[start_offset - 1];
        prev_positions.insert(
            make_tuple(c,
                       (*start_node).node->id(),
                       (*start_node).backward,
                       start_offset - 1));
    }

    // find the kmer end
    int pos = start_offset; // point at start of kmer
    bool first_in_path = true;
    // while we're not through with the path
    while (np != path.end()) {
        NodeTraversal n = *np;
        // Every time we look into a new node, we add in the amount of sequence
        // in that node. So this is the number of bases after the start of our
        // kmer up through the node we're looking at.
        // TODO: rename this to something more descriptive.
        int newpos = pos + n.node->sequence().size();

        // QUESTION:
        // Would the edge_max constraint cause us to drop the next kmer?
        // ANSWER:
        // It will when the count of edges in the implied path would be >edge_max
        // So assemble these paths and answer that question.
        // --- you can assemble any one of them, or simpler,
        // the question to ask is 1) are we losing an edge crossing on the left?
        // 2) are we gaining a new edge crossing on the right?

        if (first_in_path) {
            // Start with newpos equal to the amount of sequence in the node the
            // kmer is on, minus that consumed by the offset to the start of the
            // kmer.
            newpos = n.node->sequence().size() - pos;
            first_in_path = false;
        }
        if (newpos == kmer.size()) {
            // The kmer ended right at the end of a node

            // Save the end node instance and offset for the kmer. Remember that
            // end_offset counts rightward.
            end_node = np;
            end_offset = 0;

            // we might lose the next kmer
            // if the current path crosses the edge max number of edges
            // and the next doesn't lose an edge on the left
            // TODO: implement that

            // Look at all the nodes that might come next.
            vector<NodeTraversal> next_nodes;
            nodes_next(n, next_nodes);
            for (auto m : next_nodes) {
                // How long is this next node?
                size_t node_length = m.node->sequence().size();
                // If the next node is backward, get the rc of its last character. Else get its first.
                char c = m.backward ? reverse_complement(m.node->sequence()[node_length - 1]) : m.node->sequence()[0];
                // We're going to the 0 offset on this node, no matter what orientation that actually is.
                next_positions.insert(
                    make_tuple(
                        c,
                        m.node->id(),
                        m.backward,
                        0));
            }
            break;
        } else if (newpos > kmer.size()) {
            // There is at least one character in this node after the kmer

            // How long is this node?
            size_t node_length = n.node->sequence().size();

            // How far in from the left of the oriented node is the first base not in the kmer?
            // We know this won't be 0, or we'd be in the if branch above.
            int off = node_length - (newpos - kmer.size());

            // Save the end node instance and offset for the kmer. Remember that
            // end_offset counts rightward, and we need to walk in 1 more to
            // actually be on the kmer (which we accomplish by omitting the -1
            // we would usually use when reversing).
            end_node = np;
            end_offset = node_length - off;

            // Fill in the next characters and positions. Remember that off
            // points to the first character after the end of the kmer.
            char c = n.backward ?
                reverse_complement(n.node->sequence()[node_length - off - 1]) :
                n.node->sequence()[off];
            next_positions.insert(
                make_tuple(c,
                           n.node->id(),
                           n.backward,
                           off));
            break;
        } else {
            pos = newpos;
            ++np;
        }
    }

    if(end_node == path.end()) {
        cerr << "Could not find end node for " << kmer << " at " << start_offset << " into "
             << (*start_node).node->id() << " " << (*start_node).backward << endl;
        for(auto t : path) {
            cerr << t.node->id() << " " << t.backward << endl;
        }
        assert(false);
    }
}

void VG::gcsa_handle_node_in_graph(Node* node, int kmer_size, int edge_max, int stride,
                                   bool forward_only,
                                   Node* head_node, Node* tail_node,
                                   function<void(KmerPosition&)> lambda) {

    // Go through all the kpaths of this node, and produce the GCSA2 kmers on both strands that start in this node.
#ifdef debug
    cerr << "Visiting node " << node->id() << endl;
#endif
    // This function runs in only one thread on a given node, so we can keep
    // our cache here. We gradually fill in each KmerPosition with all the
    // next positions and characters reachable with its string from its
    // orientation and offset along that strand in this node.
    map<tuple<string, bool, int32_t>, KmerPosition> cache;

    // We're going to visit every kmer of the node and run this:
    function<void(string&, list<NodeTraversal>::iterator, int, list<NodeTraversal>&, VG&)>
        visit_kmer = [&cache, &kmer_size, &edge_max, &node, &forward_only, &head_node, &tail_node, this]
        (string& kmer, list<NodeTraversal>::iterator start_node, int start_pos,
         list<NodeTraversal>& path, VG& graph) {

        // We should never see negative offset kmers; _for_each_kmer ought to
        // have turned them around for positive offsets on the opposite strand.
        assert(start_pos >= 0);

        // todo, handle edge bounding
        // we need to check if the previous or next kmer will be excluded based on
        // edge bounding
        // if so, we should connect to the source or sink node

        // Get the information from the graph about what's before and after
        // this kmer, and where it ends.
        list<NodeTraversal>::iterator end_node;
        int32_t end_pos; // This counts in from the right of end_node.

        set<tuple<char, id_t, bool, int32_t>> prev_positions;
        set<tuple<char, id_t, bool, int32_t>> next_positions;

        // Fill in prev_chars, next_chars, prev_positions, and next_positions for the kmer by walking the path.
        graph.kmer_context(kmer,
                           kmer_size,
                           edge_max,
                           forward_only,
                           path,
                           start_node,
                           start_pos,
                           end_node,
                           end_pos,
                           prev_positions,
                           next_positions);

        if(start_node->node == node) {
            if (forward_only && start_node->backward) return;
            // This kmer starts on the node we're currently processing.
            // Store the information about it's forward orientation.

            // Get the KmerPosition to fill, creating it if it doesn't exist already.
            auto cache_key = make_tuple(kmer, start_node->backward, start_pos);
#ifdef debug
            if(cache.count(cache_key)) {
                cerr << "F: Adding to " << kmer << " at " << start_node->node->id() << " " << start_node->backward
                     << " offset " << start_pos << endl;
            } else {
                cerr << "F: Creating " << kmer << " at " << start_node->node->id() << " " << start_node->backward
                     << " offset " << start_pos << endl;
            }
#endif
            KmerPosition& forward_kmer = cache[cache_key];

            // fix up the context by swapping reverse-complements of the head and tail node
            // with their forward versions
            auto next_positions_copy = next_positions;
            next_positions.clear();
            for (auto p : next_positions_copy) {
                char c = get<0>(p);
                id_t node_id = get<1>(p);
                bool is_backward = get<2>(p);
                int32_t pos = get<3>(p);
                if (node_id == tail_node->id() && is_backward) {
                    // fixe the char as well...
                    c = head_node->sequence()[0]; // let's hope it has sequence
                    node_id = head_node->id();
                    is_backward = false;
                } else if (node_id == head_node->id() && is_backward) {
                    c = tail_node->sequence()[0]; // let's hope it has sequence
                    node_id = tail_node->id();
                    is_backward = false;
                }
                next_positions.insert(make_tuple(c, node_id, is_backward, pos));
            }

            // Add in the kmer string
            if (forward_kmer.kmer.empty()) forward_kmer.kmer = kmer;

            // Add in the start position
            if (forward_kmer.pos.empty()) {
                // And the distance from the end of the kmer to the end of its ending node.
                stringstream ps;
                if (start_node->node->id() == tail_node->id() && start_node->backward) {
                    ps << head_node->id() << ":" << start_pos;
                } else if (start_node->node->id() == head_node->id() && start_node->backward) {
                    ps << tail_node->id() << ":" << start_pos;
                } else {
                    ps << start_node->node->id() << ":" << (start_node->backward?"-":"") << start_pos;
                }
                forward_kmer.pos = ps.str();
            }

            // Add in the prev and next characters.
            for (auto& t : prev_positions) {
                char c = get<0>(t);
                forward_kmer.prev_chars.insert(c);
            }
            for (auto& t : next_positions) {
                char c = get<0>(t);
                forward_kmer.next_chars.insert(c);
            }

            // Add in the next positions
            for (auto& p : next_positions) {
                // Figure out if the forward kmer should go next to the forward or reverse copy of the next node.
                id_t target_node = get<1>(p);
                bool target_node_backward = get<2>(p);
                int32_t target_off = get<3>(p);
                // Say we go to it at the correct offset
                stringstream ps;
                ps << target_node << ":" << (target_node_backward?"-":"") << target_off;
                forward_kmer.next_positions.insert(ps.str());
            }
        }

        if(end_node->node == node && !forward_only) {
            // This kmer ends on the node we're currently processing.
            // *OR* this kmer starts on the head node (in which case, we'd never see it from the other end)
            // Store the information about it's reverse orientation.

            // Get the KmerPosition to fill, creating it if it doesn't exist
            // already. We flip the backwardness because we look at the kmer
            // the other way, but since end_pos already counts from the end,
            // we don't touch it.
            auto cache_key = make_tuple(reverse_complement(kmer), !end_node->backward, end_pos);
#ifdef debug
            if(cache.count(cache_key)) {
                cerr << "R: Adding to " << reverse_complement(kmer) << " at " << end_node->node->id()
                     << " " << !(*end_node).backward << " offset " << end_pos << endl;
            } else {
                cerr << "R: Creating " << reverse_complement(kmer) << " at " << end_node->node->id()
                << " " << !(*end_node).backward << " offset " << end_pos << endl;
            }
#endif
            KmerPosition& reverse_kmer = cache[cache_key];

            // fix up the context by swapping reverse-complements of the head and tail node
            // with their forward versions
            auto prev_positions_copy = prev_positions;
            prev_positions.clear();
            for (auto p : prev_positions_copy) {
                char c = get<0>(p);
                id_t node_id = get<1>(p);
                bool is_backward = get<2>(p);
                int32_t pos = get<3>(p);
                if (node_id == tail_node->id() && !is_backward) {
                    node_id = head_node->id();
                    is_backward = true;
                    pos = graph.get_node(node_id)->sequence().size() - pos - 1;
                } else if (node_id == head_node->id() && !is_backward) {
                    node_id = tail_node->id();
                    is_backward = true;
                    pos = tail_node->sequence().size() - pos - 1;
                } else {
                    pos = graph.get_node(node_id)->sequence().size() - pos - 1;
                }
                prev_positions.insert(make_tuple(c, node_id, is_backward, pos));
            }

            // Add in the kmer string
            if (reverse_kmer.kmer.empty()) reverse_kmer.kmer = reverse_complement(kmer);

            // Add in the start position
            if (reverse_kmer.pos.empty()) {
                // Use the other node ID, facing the other way
                stringstream ps;
                if (end_node->node->id() == tail_node->id() && !end_node->backward) {
                    ps << head_node->id() << ":" << end_pos;
                } else if (end_node->node->id() == head_node->id() && !end_node->backward) {
                    ps << tail_node->id() << ":" << end_pos;
                } else {
                    ps << (*end_node).node->id() << ":" << (!end_node->backward?"-":"") << end_pos;
                }
                // And the distance from the end of the kmer to the end of its ending node.
                reverse_kmer.pos = ps.str();
            }

            // Add in the prev and next characters.
            // fixme ... reverse complements things that should be translated to the head or tail node
            for (auto& t : prev_positions) {
                char c = get<0>(t);
                reverse_kmer.next_chars.insert(reverse_complement(c));
            }
            for (auto& t : next_positions) {
                char c = get<0>(t);
                reverse_kmer.prev_chars.insert(reverse_complement(c));
            }

            // Add in the next positions (using the prev positions since we're reversing)
            for (auto& p : prev_positions) {
                // Figure out if the reverse kmer should go next to the forward or reverse copy of the next node.
                id_t target_node = get<1>(p);
                bool target_node_backward = get<2>(p);
                int32_t off = get<3>(p);

                // Say we go to it at the correct offset
                stringstream ps;
                ps << target_node << ":" << (!target_node_backward?"-":"") << off;
                reverse_kmer.next_positions.insert(ps.str());
            }
        }
    };

    // Now we visit every kmer of this node and fill in the cache. Don't
    // allow negative offsets; force them to be converted to positive
    // offsets on the reverse strand. But do allow different paths that
    // produce the same kmer, since GCSA2 needs those.
    for_each_kmer_of_node(node, kmer_size, edge_max, visit_kmer, stride, true, false);

    // Now that the cache is full and correct, containing each kmer starting
    // on either strand of this node, send out all its entries.
    for(auto& kv : cache) {
        lambda(kv.second);
    }

}

void VG::for_each_gcsa_kmer_position_parallel(int kmer_size, int edge_max, int stride,
                                              bool forward_only,
                                              id_t& head_id, id_t& tail_id,
                                              function<void(KmerPosition&)> lambda) {

    show_progress = show_progress;
    progress_message = "processing kmers of " + name;
    Node* head_node = nullptr, *tail_node = nullptr;
    if(head_id == 0) {
        assert(tail_id == 0); // they should be only set together
        // This is the first graph. We need to fill in the IDs.
        add_start_end_markers(kmer_size, '#', '$', head_node, tail_node, head_id, tail_id);
        // Save its ID
        head_id = head_node->id();
        tail_id = tail_node->id();
    } else {
        // Add the existing start/end node
        id_t maxid = max_node_id();
        if(head_id <= maxid || tail_id <= maxid) {
            // If the ID we got for the node when we made it in the
            // first graph is too small, we have to complain. It would be
            // nice if we could make a path through all the graphs, get the
            // max ID, and then use that to determine the new node ID.
            cerr << "error:[for_each_gcsa_kmer_position_parallel] created a start/end "
                 << "node in first graph with id used by later graph " << name
                 << ". Put the graph with the largest node id first and try again." << endl;
            exit(1);
        }
        add_start_end_markers(kmer_size, '#', '$', head_node, tail_node, head_id, tail_id);
    }

    // Copy the head and tail nodes
    Node local_head_node = *head_node;
    Node local_tail_node = *tail_node;

    // Now we have to drop unconnected start/end nodes, so we don't produce
    // kmers straight from # to $

    // Remember if we don't remove them from the graph here, because we'll need
    // to remove them from the graph later.
    bool head_node_in_graph = true;
    bool tail_node_in_graph = true;

    vector<Edge*> edges;
    edges_of_node(head_node, edges);
    if(edges.empty()) {
        // The head node is floating disconnected.
#ifdef debug
        cerr << "Head node wasn't used. Excluding from graph." << endl;
#endif
        destroy_node(head_node);
        // We still need one to represent the RC of the tail node, but it can't actually be in the graph.
        head_node = &local_head_node;
        head_node_in_graph = false;
    }
    edges.clear();
    edges_of_node(tail_node, edges);
    if(edges.empty()) {
        // The tail node is floating disconnected.
#ifdef debug
        cerr << "Tail node wasn't used. Excluding from graph." << endl;
#endif
        destroy_node(tail_node);
        tail_node = &local_tail_node;
        tail_node_in_graph = false;
    }

    if(forward_only && (!head_node_in_graph || !tail_node_in_graph)) {
        // TODO: break in arbitrarily if doing forward-only indexing and there's
        // no place to attach one of the start/end nodes.
        cerr << "error:[for_each_gcsa_kmer_position_parallel] attempted to forward-only index a graph "
            "that has only heads and no tails, or only tails and no heads. Only one of the start and "
            "end nodes could be attached." << endl;
        exit(1);
    }

    // Actually find the GCSA2 kmers. The head and tail node pointers point to
    // things, but the graph is only guaranteed to actually own one of those
    // things.
    for_each_node_parallel(
        [kmer_size, edge_max, stride, forward_only,
         head_node, tail_node, lambda, this](Node* node) {
            gcsa_handle_node_in_graph(node, kmer_size, edge_max, stride, forward_only,
                                      head_node, tail_node, lambda);
        });

    // cleanup
    if(head_node_in_graph) {
        destroy_node(head_node);
    }
    if(tail_node_in_graph) {
        destroy_node(tail_node);
    }
}

void VG::get_gcsa_kmers(int kmer_size, int edge_max, int stride,
                        bool forward_only,
                        vector<gcsa::KMer>& kmers_out,
                        id_t& head_id, id_t& tail_id) {

    // TODO: This function goes through an internal string format that should
    // really be replaced by making some API changes to gcsa2.

    // We need an alphabet to parse the internal string format
    const gcsa::Alphabet alpha;

    // Each thread is going to make its own KMers, then we'll concatenate these all together at the end.
    vector<vector<gcsa::KMer>> thread_outputs;

#pragma omp parallel
    {
#pragma omp single
        {
            // Become parallel, get our number of threads, and make one of them make the per-thread outputs big enough.
            thread_outputs.resize(omp_get_num_threads());
        }
    }

    auto convert_kmer = [&thread_outputs, &alpha, &head_id, &tail_id](KmerPosition& kp) {
        // Convert this KmerPosition to several gcsa::Kmers, and save them in thread_outputs

        // We need to make this kmer into a series of tokens
        vector<string> tokens;

        // First the kmer
        tokens.push_back(kp.kmer);

        // Then the node id:offset
        tokens.push_back(kp.pos);

        // Then the comma-separated preceeding characters. See <http://stackoverflow.com/a/18427254/402891>
        stringstream preceeding;
        copy(kp.prev_chars.begin(), kp.prev_chars.end(), ostream_iterator<char>(preceeding, ","));
        if(kp.prev_chars.empty()) {
            // If we don't have any previous characters, we come from "$"
            preceeding << "$";
        }
        tokens.push_back(preceeding.str());

        // And the comma-separated subsequent characters.
        stringstream subsequent;
        copy(kp.next_chars.begin(), kp.next_chars.end(), ostream_iterator<char>(subsequent, ","));
        if(kp.next_chars.empty()) {
            // If we don't have any next characters, we go to "#"
            subsequent << "#";
        }
        tokens.push_back(subsequent.str());

        // Finally, each of the node id:offset positions you can go to next (the successors).
        tokens.insert(tokens.end(), kp.next_positions.begin(), kp.next_positions.end());

        if (kp.next_positions.empty()) {
            // If we didn't have any successors, we have to say we go to the start of the start node
            tokens.push_back(to_string(tail_id) + ":0");
        }

        for(size_t successor_index = 4; successor_index < tokens.size(); successor_index++) {
            // Now make a GCSA KMer for each of those successors, by passing the
            // tokens, the alphabet, and the index in the tokens of the
            // successor.

            thread_outputs[omp_get_thread_num()].emplace_back(tokens, alpha, successor_index);

            // Kmers that go to the sink/have stop characters still need to be marked as sorted.
            if(kp.kmer.rfind('$') != string::npos) {
                //(*(thread_outputs[omp_get_thread_num()].rbegin())).makeSorted();
#pragma omp critical
                {
                   // cout << "Marked " << *(thread_outputs[omp_get_thread_num()].rbegin()) << " as sorted early" << endl;
                }
            }
        }

    };

    // Run on each KmerPosition. This populates start_end_id, if it was 0, before calling convert_kmer.
    for_each_gcsa_kmer_position_parallel(kmer_size, edge_max, stride,
                                         forward_only,
                                         head_id, tail_id,
                                         convert_kmer);

    for(auto& thread_output : thread_outputs) {
        // Now throw everything into the output vector
        kmers_out.insert(kmers_out.end(),
                         make_move_iterator(thread_output.begin()),
                         make_move_iterator(thread_output.end()));
    }

    for(auto& kmer : kmers_out) {
        // Mark kmers that go to the sink node as "sorted", since they have stop
        // characters in them and can't be extended.
        // If we don't do this GCSA will get unhappy and we'll see random segfalts and stack smashing errors
        if(gcsa::Node::id(kmer.to) == tail_id && gcsa::Node::offset(kmer.to) > 0) {
            kmer.makeSorted();
        }
    }

}

gcsa::GCSA* VG::build_gcsa_index(int kmer_size, bool forward_only,
                                 size_t doubling_steps,
                                 size_t size_limit) {
    vector<gcsa::KMer> kmers;
    id_t head_id=0, tail_id=0;
    get_gcsa_kmers(kmer_size, 0, 1, forward_only,
                   kmers, head_id, tail_id);
    gcsa::GCSA* result = new gcsa::GCSA(kmers, kmer_size, doubling_steps, size_limit);
    return result;
}

void VG::prune_complex_with_head_tail(int path_length, int edge_max) {
    Node* head_node = NULL;
    Node* tail_node = NULL;
    add_start_end_markers(path_length, '#', '$', head_node, tail_node);
    prune_complex(path_length, edge_max, head_node, tail_node);
    destroy_node(head_node);
    destroy_node(tail_node);
}

void VG::prune_complex(int path_length, int edge_max, Node* head_node, Node* tail_node) {
    vector<set<NodeTraversal> > prev_maxed_nodes;
    vector<set<NodeTraversal> > next_maxed_nodes;
#pragma omp parallel
    {
#pragma omp single
        {
            int threads = omp_get_num_threads();
            prev_maxed_nodes.resize(threads);
            next_maxed_nodes.resize(threads);
        }
    }
    auto prev_maxed = [this, &prev_maxed_nodes](NodeTraversal node) {
        int tid = omp_get_thread_num();
        prev_maxed_nodes[tid].insert(node);
    };
    auto next_maxed = [this, &next_maxed_nodes](NodeTraversal node) {
        int tid = omp_get_thread_num();
        next_maxed_nodes[tid].insert(node);
    };
    auto noop = [](list<NodeTraversal>::iterator node, list<NodeTraversal>& path) { };
    for_each_kpath_parallel(path_length, edge_max, prev_maxed, next_maxed, noop);

    // What nodes will we destroy because we got into them with too much complexity?
    set<Node*> to_destroy;

    set<NodeTraversal> prev;
    for (auto& p : prev_maxed_nodes) {
        for (auto node : p) {
            prev.insert(node);
        }
    }
    for (auto node : prev) {
        // This node was prev-maxed, meaning we tried to go left into it and couldn't.
        // We want to connect the end of the head node to everywhere we tried to come left from.

        // We drop any links going into the other side of this node.

        if(node.backward) {
            // Going left into it means coming to its start. True flags on the
            // edges mean connecting to the start of the other nodes.
            for (auto& e : edges_start(node.node)) {
                create_edge(e.first, head_node->id(), e.second, true);
            }
        } else {
            // Going left into it means coming to its end. True flags on the
            // edges mean connecting to the ends of the other nodes.
            for (auto& e : edges_end(node.node)) {
                create_edge(head_node->id(), e.first, false, e.second);
            }
        }

        // Remember to estroy the node. We might come to the same node in two directions.
        to_destroy.insert(node.node);
    }

    set<NodeTraversal> next;
    for (auto& n : next_maxed_nodes) {
        for (auto node : n) {
            next.insert(node);
        }
    }
    for (auto node : next) {
        // This node was next-maxed, meaning we tried to go right into it and couldn't.
        // We want to connect the start of the tail node to everywhere we tried to come right from.

        // We drop any links going into the other side of this node.

        if(node.backward) {
            // Going right into it means coming to its end. True flags on the
            // edges mean connecting to the end of the other nodes.
            for (auto& e : edges_end(node.node)) {
                create_edge(tail_node->id(), e.first, false, e.second);
            }
        } else {
            // Going right into it means coming to its start. True flags on the
            // edges mean connecting to the starts of the other nodes.
            for (auto& e : edges_start(node.node)) {
                create_edge(e.first, head_node->id(), e.second, true);
            }
        }

        // Remember to estroy the node. We might come to the same node in two directions.
        to_destroy.insert(node.node);
    }

    for(Node* n : to_destroy) {
        // Destroy all the nodes we wanted to destroy.
        if(n == head_node || n == tail_node) {
            // Keep these around
            continue;
        }

        // First delete any paths that touch it.
        // TODO: split the paths in two at this node somehow
        set<string> paths_to_remove;
        for(auto path_and_mapping : paths.get_node_mapping(n)) {
            paths_to_remove.insert(path_and_mapping.first);
        }
        paths.remove_paths(paths_to_remove);

        // Actually destroy the node
        destroy_node(n);
    }

    for (auto* n : head_nodes()) {
        if (n != head_node) {
            // Fix up multiple heads with a left-to-right edge
            create_edge(head_node, n);
        }
    }
    for (auto* n : tail_nodes()) {
        if (n != tail_node) {
            // Fix up multiple tails with a left-to-right edge
            create_edge(n, tail_node);
        }
    }
}

void VG::prune_short_subgraphs(size_t min_size) {
    list<VG> subgraphs;
    disjoint_subgraphs(subgraphs);
    for (auto& g : subgraphs) {
        vector<Node*> heads;
        g.head_nodes(heads);
        // calculate length
        // if < N
        if (g.total_length_of_nodes() < min_size) {
            //cerr << "removing" << endl;
            g.for_each_node([this](Node* n) {
                    // remove from this graph a node of the same id
                    //cerr << n->id() << endl;
                    this->destroy_node(n->id());
                });
        }
    }
}

/*
// todo
void VG::prune_complex_subgraphs(size_t ) {

}
*/

void VG::collect_subgraph(Node* start_node, set<Node*>& subgraph) {

    // add node to subgraph
    subgraph.insert(start_node);

    set<Node*> checked;
    set<Node*> to_check;
    to_check.insert(start_node);

    while (!to_check.empty()) {
        // for each predecessor of node
        set<Node*> curr_check = to_check;
        to_check.clear();
        for (auto* node : curr_check) {
            if (checked.count(node)) {
                continue;
            } else {
                checked.insert(node);
            }
            vector<NodeTraversal> prev;
            nodes_prev(node, prev);
            for (vector<NodeTraversal>::iterator p = prev.begin(); p != prev.end(); ++p) {
            // if it's not already been examined, collect its neighborhood
                if (!subgraph.count((*p).node)) {
                    subgraph.insert((*p).node);
                    to_check.insert((*p).node);
                }
            }
            // for each successor of node
            vector<NodeTraversal> next;
            nodes_next(node, next);
            for (vector<NodeTraversal>::iterator n = next.begin(); n != next.end(); ++n) {
                if (!subgraph.count((*n).node)) {
                    subgraph.insert((*n).node);
                    to_check.insert((*n).node);
                }
            }
        }
    }
    //cerr << "node " << start_node->id() << " subgraph size " << subgraph.size() << endl;
}

void VG::disjoint_subgraphs(list<VG>& subgraphs) {
    vector<Node*> heads;
    head_nodes(heads);
    map<Node*, set<Node*> > subgraph_by_head;
    map<Node*, set<Node*>* > subgraph_membership;
    // start at the heads, but keep in mind that we need to explore fully
    for (vector<Node*>::iterator h = heads.begin(); h != heads.end(); ++h) {
        if (subgraph_membership.find(*h) == subgraph_membership.end()) {
            set<Node*>& subgraph = subgraph_by_head[*h];
            collect_subgraph(*h, subgraph);
            for (set<Node*>::iterator n = subgraph.begin(); n != subgraph.end(); ++n) {
                subgraph_membership[*n] = &subgraph;
            }
        }
    }
    for (map<Node*, set<Node*> >::iterator g = subgraph_by_head.begin();
         g != subgraph_by_head.end(); ++ g) {
        set<Node*>& nodes = g->second;
        set<Edge*> edges;
        edges_of_nodes(nodes, edges);
        subgraphs.push_back(VG(nodes, edges));
    }
}

bool VG::is_head_node(id_t id) {
    return is_head_node(get_node(id));
}

bool VG::is_head_node(Node* node) {
    return start_degree(node) == 0;
}

void VG::head_nodes(vector<Node*>& nodes) {
    for (int i = 0; i < graph.node_size(); ++i) {
        Node* n = graph.mutable_node(i);
        if (is_head_node(n)) {
            nodes.push_back(n);
        }
    }
}

vector<Node*> VG::head_nodes(void) {
    vector<Node*> heads;
    head_nodes(heads);
    return heads;
}

/*
// TODO item
int32_t VG::distance_to_head(Position pos, int32_t limit) {
}
*/

int32_t VG::distance_to_head(NodeTraversal node, int32_t limit) {
    set<NodeTraversal> seen;
    return distance_to_head(node, limit, 0, seen);
}

int32_t VG::distance_to_head(NodeTraversal node, int32_t limit, int32_t dist, set<NodeTraversal>& seen) {
    NodeTraversal n = node;
    if (seen.count(n)) return -1;
    seen.insert(n);
    if (limit <= 0) {
        return -1;
    }
    if (is_head_node(n.node)) {
        return dist;
    }
    for (auto& trav : nodes_prev(n)) {
        size_t l = trav.node->sequence().size();
        size_t t = distance_to_head(trav, limit-l, dist+l, seen);
        if (t != -1) {
            return t;
        }
    }
    return -1;
}

bool VG::is_tail_node(id_t id) {
    return is_tail_node(get_node(id));
}

bool VG::is_tail_node(Node* node) {
    return end_degree(node) == 0;
}

void VG::tail_nodes(vector<Node*>& nodes) {
    for (int i = 0; i < graph.node_size(); ++i) {
        Node* n = graph.mutable_node(i);
        if (is_tail_node(n)) {
            nodes.push_back(n);
        }
    }
}

vector<Node*> VG::tail_nodes(void) {
    vector<Node*> tails;
    tail_nodes(tails);
    return tails;
}

int32_t VG::distance_to_tail(NodeTraversal node, int32_t limit) {
    set<NodeTraversal> seen;
    return distance_to_tail(node, limit, 0, seen);
}

int32_t VG::distance_to_tail(NodeTraversal node, int32_t limit, int32_t dist, set<NodeTraversal>& seen) {
    NodeTraversal n = node;
    if (seen.count(n)) return -1;
    seen.insert(n);
    if (limit <= 0) {
        return -1;
    }
    if (is_tail_node(n.node)) {
        return dist;
    }
    for (auto& trav : nodes_next(n)) {
        size_t l = trav.node->sequence().size();
        size_t t = distance_to_tail(trav, limit-l, dist+l, seen);
        if (t != -1) {
            return t;
        }
    }
    return -1;
}

void VG::wrap_with_null_nodes(void) {
    vector<Node*> heads;
    head_nodes(heads);
    Node* head = create_node("");
    for (vector<Node*>::iterator h = heads.begin(); h != heads.end(); ++h) {
        create_edge(head, *h);
    }

    vector<Node*> tails;
    tail_nodes(tails);
    Node* tail = create_node("");
    for (vector<Node*>::iterator t = tails.begin(); t != tails.end(); ++t) {
        create_edge(*t, tail);
    }
}

/**
Order and orient the nodes in the graph using a topological sort.

We use a bidirected adaptation of Kahn's topological sort (1962), which can handle components with no heads or tails.

L ← Empty list that will contain the sorted and oriented elements
S ← Set of nodes which have been oriented, but which have not had their downstream edges examined
N ← Set of all nodes that have not yet been put into S

while N is nonempty do
    remove a node from N, orient it arbitrarily, and add it to S
        (In practice, we use "seeds": the heads and any nodes we have seen that had too many incoming edges)
    while S is non-empty do
        remove an oriented node n from S
        add n to tail of L
        for each node m with an edge e from n to m do
            remove edge e from the graph
            if m has no other edges to that side then
                orient m such that the side the edge comes to is first
                remove m from N
                insert m into S
            otherwise
                put an oriented m on the list of arbitrary places to start when S is empty
                    (This helps start at natural entry points to cycles)
    return L (a topologically sorted order and orientation)
*/
void VG::topological_sort(deque<NodeTraversal>& l) {
    //assert(is_valid());

#ifdef debug
#pragma omp critical (cerr)
    cerr << "=====================STARTING SORT==========================" << endl;
#endif

    // using a map instead of a set ensures a stable sort across different systems
    map<id_t, NodeTraversal> s;

    // We find the head and tails, if there are any
    vector<Node*> heads;
    head_nodes(heads);
    vector<Node*> tails;
    tail_nodes(tails);

    // We'll fill this in with the heads, so we can orient things according to
    // them first, and then arbitrarily. We ignore tails since we only orient
    // right from nodes we pick.
    // Maps from node ID to first orientation we suggested for it.
    map<id_t, NodeTraversal> seeds;
    for(Node* head : heads) {
        seeds[head->id()] = NodeTraversal(head, false);
    }

    // We will use a std::map copy of node_by_id as our set of unvisited nodes,
    // and remove index entries from it when we visit nodes. It will be rebuilt
    // when we rebuild the indexes later. We know its order will be fixed across
    // systems.
    map<id_t, Node*> unvisited;
    // Fill it in, we can't use the copy constructor since std::map doesn't speak vg::hash_map
    // TODO: is the vg::hash_map order fixed across systems?
    for_each_node([&](Node* node) {
        unvisited[node->id()] = node;
    });

    // How many nodes have we ordered and oriented?
    id_t seen = 0;

    while(!unvisited.empty()) {

        // Put something in s. First go through seeds until we can find one
        // that's not already oriented.
        while(s.empty() && !seeds.empty()) {
            // Look at the first seed
            NodeTraversal first_seed = (*seeds.begin()).second;

            if(unvisited.count(first_seed.node->id())) {
                // We have an unvisited seed. Use it
#ifdef debug
#pragma omp critical (cerr)
                cerr << "Starting from seed " << first_seed.node->id() << " orientation " << first_seed.backward << endl;
#endif

                s[first_seed.node->id()] = first_seed;
                unvisited.erase(first_seed.node->id());
            }
            // Whether we used the seed or not, don't keep it around
            seeds.erase(seeds.begin());
        }

        if(s.empty()) {
            // If we couldn't find a seed, just grab any old node.
            // Since map order is stable across systems, we can take the first node by id and put it locally forward.
#ifdef debug
#pragma omp critical (cerr)
            cerr << "Starting from arbitrary node " << unvisited.begin()->first << " locally forward" << endl;
#endif

            s[unvisited.begin()->first] = NodeTraversal(unvisited.begin()->second, false);
            unvisited.erase(unvisited.begin()->first);
        }

        while (!s.empty()) {
            // Grab an oriented node
            NodeTraversal n = s.begin()->second;
            s.erase(n.node->id());
            l.push_back(n);
            ++seen;
#ifdef debug
#pragma omp critical (cerr)
            cerr << "Using oriented node " << n.node->id() << " orientation " << n.backward << endl;
#endif

            // See if it has an edge from its start to the start of some node
            // where both were picked as places to break into cycles. A
            // reversing self loop on a cycle entry point is a special case of
            // this.
            vector<NodeTraversal> prev;
            nodes_prev(n, prev);
            for(NodeTraversal& prev_node : prev) {
                if(!unvisited.count(prev_node.node->id())) {
#ifdef debug
#pragma omp critical (cerr)
                    cerr << "\tHas left-side edge to cycle entry point " << prev_node.node->id()
                         << " orientation " << prev_node.backward << endl;
#endif

                    // Unindex it
                    Edge* bad_edge = get_edge(prev_node, n);
#ifdef debug
#pragma omp critical (cerr)
                    cerr << "\t\tEdge: " << bad_edge << endl;
#endif
                    unindex_edge_by_node_sides(bad_edge);
                }
            }

            // All other connections and self loops are handled by looking off the right side.

            // See what all comes next, minus deleted edges.
            vector<NodeTraversal> next;
            nodes_next(n, next);

            for(NodeTraversal& next_node : next) {

#ifdef debug
#pragma omp critical (cerr)
                cerr << "\tHas edge to " << next_node.node->id() << " orientation " << next_node.backward << endl;
#endif

                // Unindex the edge connecting these nodes in this order and
                // relative orientation, so we can't traverse it with
                // nodes_next/nodes_prev.

                // Grab the edge
                Edge* edge = get_edge(n, next_node);

#ifdef debug
#pragma omp critical (cerr)
                cerr << "\t\tEdge: " << edge << endl;
#endif

                // Unindex it
                unindex_edge_by_node_sides(edge);

                if(unvisited.count(next_node.node->id())) {
                    // We haven't already started here as an arbitrary cycle entry point

#ifdef debug
#pragma omp critical (cerr)
                    cerr << "\t\tAnd node hasn't been visited yet" << endl;
#endif

                    if(node_count_prev(next_node) == 0) {

#ifdef debug
#pragma omp critical (cerr)
                        cerr << "\t\t\tIs last incoming edge" << endl;
#endif
                        // Keep this orientation and put it here
                        s[next_node.node->id()] = next_node;
                        // Remember that we've visited and oriented this node, so we
                        // don't need to use it as a seed.
                        unvisited.erase(next_node.node->id());

                    } else if(!seeds.count(next_node.node->id())) {
                        // We came to this node in this orientation; when we need a
                        // new node and orientation to start from (i.e. an entry
                        // point to the node's cycle), we might as well pick this
                        // one.
                        // Only take it if we don't already know of an orientation for this node.
                        seeds[next_node.node->id()] = next_node;

#ifdef debug
#pragma omp critical (cerr)
                        cerr << "\t\t\tSuggests seed " << next_node.node->id() << " orientation " << next_node.backward << endl;
#endif
                    }
                } else {
#ifdef debug
#pragma omp critical (cerr)
                    cerr << "\t\tAnd node was already visited (to break a cycle)" << endl;
#endif
                }
            }

            // The caller may put us in a progress context with the denominator
            // being the number of nodes in the graph.
            update_progress(seen);
        }
    }

    // There should be no edges left in the index
    if(!edges_on_start.empty() || !edges_on_end.empty()) {
#pragma omp critical (cerr)
        {
            cerr << "Error: edges remaining after topological sort and cycle breaking" << endl;

            // Dump the edges in question
            for(auto& on_start : edges_on_start) {
                cerr << "start: " << on_start.first << endl;
                for(auto& other_end : on_start.second) {
                    cerr << "\t" << other_end.first << " " << other_end.second << endl;
                }
            }
            for(auto& on_end : edges_on_end) {
                cerr << "end: " << on_end.first << endl;
                for(auto& other_end : on_end.second) {
                    cerr << "\t" << other_end.first << " " << other_end.second << endl;
                }
            }
            cerr << "By Sides:" << endl;
            for(auto& sides_and_edge : edge_by_sides) {
                cerr << sides_and_edge.first.first << "<->" << sides_and_edge.first.second << endl;
            }

            // Dump the whole graph if possible. May crash due to bad index.
            cerr << "Dumping to fail.vg" << endl;
            std::ofstream out("fail.vg");
            serialize_to_ostream(out);
            out.close();
        }
        exit(1);
    }

    // we have destroyed the graph's edge and node index to ensure its order
    // rebuild the indexes
    rebuild_indexes();
}

void VG::force_path_match(void) {
    for_each_node([&](Node* n) {
            Edit match;
            size_t seq_len = n->sequence().size();
            match.set_from_length(seq_len);
            match.set_to_length(seq_len);
            for (auto& p : paths.get_node_mapping(n)) {
                for (auto m : p.second) {
                    *m->add_edit() = match;
                }
            }
        });
}

void VG::fill_empty_path_mappings(void) {
    for_each_node([&](Node* n) {
            Edit match;
            size_t seq_len = n->sequence().size();
            match.set_from_length(seq_len);
            match.set_to_length(seq_len);
            for (auto& p : paths.get_node_mapping(n)) {
                for (auto m : p.second) {
                    if (m->edit_size() == 0) {
                        *m->add_edit() = match;
                    }
                }
            }
        });
}

// for each inverting edge
// we walk up to max_length across the inversion
// adding the forward complement of nodes we reach
// stopping if we reach an inversion back to the forward graph
// so as to ensure that sequences up to length max_length that cross the inversion
// would align to the forward component of the resulting graph
VG VG::unfold(uint32_t max_length,
              map<id_t, pair<id_t, bool> >& node_translation) {
    VG unfolded = *this;
    unfolded.flip_doubly_reversed_edges();
    if (!unfolded.has_inverting_edges()) return unfolded;
    // maps from entry id to the set of nodes
    // map from component root id to a translation
    // that maps the unrolled id to the original node and whether we've inverted or not
    // TODO
    // we need to first collect the components so we can ask quickly if a certain node is in one
    // then we need to
    set<NodeTraversal> travs_to_flip;
    //set<Edge*> edges_to_flip;
    set<pair<NodeTraversal, NodeTraversal> > edges_to_flip;
    //set<Edge*> edges_to_forward;
    set<pair<NodeTraversal, NodeTraversal> > edges_to_forward;
    //set<Edge*> edges_from_forward;
    set<pair<NodeTraversal, NodeTraversal> > edges_from_forward;

    // collect the set to invert
    // and we'll copy them over
    // then link back in according to how the inbound links work
    // so as to eliminate the reversing edges
    map<NodeTraversal, int32_t> seen;
    function<void(NodeTraversal,int32_t)> walk = [&](NodeTraversal curr,
                                                     int32_t length) {

        // check if we've passed our length limit
        // or if we've seen this node before at an earlier step
        // (in which case we're done b/c we will traverse the rest of the graph up to max_length)
        set<NodeTraversal> next;
        travs_to_flip.insert(curr);
        if (length <= 0 || seen.find(curr) != seen.end() && seen[curr] < length) {
            return;
        }
        seen[curr] = length;
        for (auto& trav : travs_from(curr)) {
            if (trav.backward) {
                edges_to_flip.insert(make_pair(curr, trav));
                walk(trav, length-trav.node->sequence().size());
            } else {
                // we would not continue, but we should retain this edge because it brings
                // us back into the forward strand
                edges_to_forward.insert(make_pair(curr, trav));
            }
        }
        for (auto& trav : travs_to(curr)) {
            if (trav.backward) {
                edges_to_flip.insert(make_pair(trav, curr));
                walk(trav, length-trav.node->sequence().size());
            } else {
                // we would not continue, but we should retain this edge because it brings
                // us back into the forward strand
                edges_from_forward.insert(make_pair(trav, curr));
            }
        }
    };

    // run over the inverting edges
    for_each_node([&](Node* node) {
            for (auto& trav : travs_of(NodeTraversal(node, false))) {
                if (trav.backward) {
                    walk(trav,  max_length);
                }
            }
        });
    // now build out the component of the graph that's reversed

    // our friend is map<id_t, pair<id_t, bool> >& node_translation;
    map<NodeTraversal, id_t> inv_translation;

    // first adding nodes that we've flipped
    for (auto t : travs_to_flip) {
        // make a new node, add it to the flattened version
        // record it in the translation
        //string seq = reverse_complement(t.node->sequence()) + ":" + convert(t.node->id());
        string seq = reverse_complement(t.node->sequence());
        id_t i = unfolded.create_node(seq)->id();
        node_translation[i] = make_pair(t.node->id(), t.backward);
        inv_translation[t] = i;
    }

    // then edges that we should translate into the reversed component
    for (auto e : edges_to_flip) {
        // both of the edges are now in nodes that have been flipped
        // we need to find the new nodes and add the natural edge
        // the edge will also be reversed
        Edge f;
        f.set_from(inv_translation[e.first]);
        f.set_to(inv_translation[e.second]);
        unfolded.add_edge(f);
    }

    // finally the edges that take us from the reversed component back to the original graph
    for (auto e : edges_to_forward) {
        Edge f;
        f.set_from(inv_translation[e.first]);//NodeTraversal(get_node(e->from()), true)]);
        f.set_to(e.second.node->id());
        unfolded.add_edge(f);
    }
    for (auto e : edges_from_forward) {
        Edge f;
        f.set_from(e.first.node->id());
        f.set_to(inv_translation[e.second]);
        unfolded.add_edge(f);
    }

    // now remove all inverting edges, so we have no more folds in the graph
    unfolded.remove_inverting_edges();

    return unfolded;
}

bool VG::has_inverting_edges(void) {
    for (id_t i = 0; i < graph.edge_size(); ++i) {
        auto& edge = graph.edge(i);
        if (!(edge.from_start() && edge.to_end())
            && (edge.from_start() || edge.to_end())) {
            return true;
        }
    }
    return false;
}

void VG::remove_inverting_edges(void) {
    set<pair<NodeSide, NodeSide>> edges;
    for_each_edge([this,&edges](Edge* edge) {
            if (!(edge->from_start() && edge->to_end())
                && (edge->from_start() || edge->to_end())) {
                edges.insert(NodeSide::pair_from_edge(edge));
            }
        });
    for (auto edge : edges) {
        destroy_edge(edge);
    }
}

bool VG::is_self_looping(NodeTraversal trav) {
    for (auto t : nodes_next(trav)) {
        if (t.node == trav.node) return true;
    }
    return false;
}

VG VG::dagify(uint32_t expand_scc_steps,
              map<id_t, pair<id_t, bool> >& node_translation,
              size_t target_min_walk_length) {

    VG dag;
    // Find the strongly connected components in the graph.
    set<set<id_t>> strong_components = strongly_connected_components();
    // map from component root id to a translation
    // that maps the unrolled id to the original node and whether we've inverted or not

    set<set<id_t>> strongly_connected_and_self_looping_components;
    set<id_t> weak_components;
    for (auto& component : strong_components) {
        // is this node a single component?
        // does this have an inversion as a child?
        // let's add in inversions
        if (component.size() == 1
            && !is_self_looping(NodeTraversal(get_node(*component.begin())))) {

            // not part of a SCC
            // copy into the new graph
            id_t id = *component.begin();
            Node* node = get_node(id);
            // this node translates to itself
            node_translation[id] = make_pair(node->id(), false);
            dag.add_node(*node);
            weak_components.insert(id);
        } else {
            strongly_connected_and_self_looping_components.insert(component);
        }
    }
    // add in the edges between the weak components
    for (auto& id : weak_components) {
        // get the edges from the graph that link it with other weak components
        for (auto e : edges_of(get_node(id))) {
            if (weak_components.count(e->from())
                && weak_components.count(e->to())) {
                dag.add_edge(*e);
            }
        }
    }

    for (auto& component : strongly_connected_and_self_looping_components) {

        // copy the SCC expand_scc_steps times, each time forwarding links from the old copy into the new
        // the result is a DAG even if the graph is initially cyclic

        // we need to record the minimum distances back to the root(s) of the graph
        // so we can (optionally) stop when we reach a given minimum minimum
        // we derive these using dynamic programming; the new min return length is
        // min(l_(i-1), \forall inbound links)
        size_t min_min_return_length = 0;
        map<Node*, size_t> min_return_length;
        // pointers to the last copy of the graph in the DAG
        map<id_t, Node*> last;
        for (uint32_t i = 0; i < expand_scc_steps+1; ++i) {
            map<id_t, Node*> curr;
            size_t curr_min_min_return_length = 0;
            // for each iteration, add in a copy of the nodes of the component
            for (auto id : component) {
                Node* node;
                if (last.empty() && !dag.has_node(id)) {
                    // create the node if we haven't already brought it over
                    node = dag.create_node(get_node(id)->sequence(), id);
                } else {
                    node = dag.create_node(get_node(id)->sequence());
                }
                curr[id] = node;
                node_translation[node->id()] = make_pair(id, false);
                // note that we don't record the lengths of these sequences
                // because we are measuring our min walk back length back to the first edge in the path
                // not the root of the graph
            }
            // preserve the edges that connect these nodes to the rest of the graph
            // And connect to the nodes in the previous component using the original edges as guide
            for (auto id : component) {
                for (auto e : edges_of(get_node(id))) {
                    if (e->from() == id && e->to() != id) {
                        // if other end is not in the component
                        if (!component.count(e->to())) {
                            // link the new node to the old one
                            Edge new_edge = *e;
                            new_edge.set_from(curr[id]->id());
                            dag.add_edge(new_edge);
                        }
                    } else if (e->to() == id && e->from() != id) {
                        // if other end is not in the component
                        if (!component.count(e->from())) {
                            // link the new node to the old one
                            Edge new_edge = *e;
                            new_edge.set_to(curr[id]->id());
                            dag.add_edge(new_edge);
                        } else if (!last.empty()) {
                            // but if we aren't in the first step
                            // and an edge is coming from a node in this component to this one
                            // add the edge that connects back to the previous node in the last copy
                            Edge new_edge = *e;
                            new_edge.set_to(curr[id]->id());
                            new_edge.set_from(last[e->from()]->id());
                            dag.add_edge(new_edge);
                            // update the min-min length
                            size_t& mm = min_return_length[curr[id]];
                            size_t inmm = curr[id]->sequence().size() + min_return_length[last[id]];
                            mm = (mm ? min(mm, inmm) : inmm);
                            curr_min_min_return_length = (curr_min_min_return_length ?
                                                          min(mm, curr_min_min_return_length)
                                                          : mm);
                        }
                    } else if (e->to() == id && e->from() == id) {
                        if (!last.empty()) {
                            // but if we aren't in the first step
                            // and an edge is coming from a node in this component to this one
                            // add the edge that connects back to the previous node in the last copy
                            Edge new_edge = *e;
                            new_edge.set_to(curr[id]->id());
                            new_edge.set_from(last[id]->id());
                            dag.add_edge(new_edge);
                            // update the min-min length
                            size_t& mm = min_return_length[curr[id]];
                            size_t inmm = curr[id]->sequence().size() + min_return_length[last[id]];
                            mm = (mm ? min(mm, inmm) : inmm);
                            curr_min_min_return_length = (curr_min_min_return_length ?
                                                          min(mm, curr_min_min_return_length)
                                                          : mm);
                        }
                    }
                }
            }
            min_min_return_length = curr_min_min_return_length;
            // finish if we've reached our target min walk length
            if (target_min_walk_length &&
                min_min_return_length >= target_min_walk_length) {
                break;
            }
            last = curr;
        }
    }

    // ensure normalized edges in output; we may introduced flipped/flipped edges
    // which are valid but can introduce problems for some algorithms
    dag.flip_doubly_reversed_edges();
    return dag;
}


// Unrolls the graph into a tree in which loops are "unrolled" into new nodes
// up to some max length away from the root node and orientations are flipped.
// A translation between the new nodes that are introduced and the old nodes and graph
// is returned so that we can map reads against this structure and resolve their mappings back
// to the original "rolled" and inverted graph.
// The graph which is returned can be seen as a tree rooted at the source node, so proper
// selection of the new root may be important for performance.
// Paths cannot be maintained provided their current implementation.
// Annotated collections of nodes, or subgraphs, may be a way to preserve the relationshp.
VG VG::backtracking_unroll(uint32_t max_length, uint32_t max_branch,
                           map<id_t, pair<id_t, bool> >& node_translation) {
    VG unrolled;
    // Find the strongly connected components in the graph.
    set<set<id_t>> strong_components = strongly_connected_components();
    // add in bits where we have inversions that we'd reach from the forward direction
    // we will "unroll" these regions as well to ensure that all is well
    // ....
    //
    map<id_t, VG> trees;
    // maps from entry id to the set of nodes
    map<id_t, set<id_t> > components;
    // map from component root id to a translation
    // that maps the unrolled id to the original node and whether we've inverted or not
    map<id_t, map<id_t, pair<id_t, bool> > > translations;
    map<id_t, map<pair<id_t, bool>, set<id_t> > > inv_translations;

    // ----------------------------------------------------
    // unroll the strong components of the graph into trees
    // ----------------------------------------------------
    //cerr << "unroll the strong components of the graph into trees" << endl;
    // Anything outside the strongly connected components can be copied directly into the DAG.
    // We can also reuse the entry nodes to the strongly connected components.
    for (auto& component : strong_components) {
        // is this node a single component?
        // does this have an inversion as a child?
        // let's add in inversions

        if (component.size() == 1) {
            // copy into the new graph
            id_t id = *component.begin();
            node_translation[id] = make_pair(id, false);
            // we will handle this node if it has an inversion originating from it
            //if (!has_inverting_edge(get_node(id))) {
            //nonoverlapping_node_context_without_paths(get_node(id), unrolled);
            unrolled.add_node(*get_node(id));
            continue;
            //}
            // otherwise we will consider it as a component
            // and it will be unrolled max_length bp
        }

        // we have a multi-node component
        // first find the entry points
        // entry points will be nodes that have connections outside of the component
        set<id_t> entries;
        set<id_t> exits;
        for (auto& n : component) {
            for (auto& e : edges_of(get_node(n))) {
                if (!component.count(e->from())) {
                    entries.insert(n);
                }
                if (!component.count(e->to())) {
                    exits.insert(n);
                }
            }
        }

        // Use backtracking search starting from each entry node of each strongly
        // connected component, keeping track of the nodes on the path from the
        // entry node to the current node:
        for (auto entrypoint : entries) {
            //cerr << "backtracking search from " << entrypoint << endl;
            // each entry point is going to make a tree
            // we can merge them later with some tricks
            // for now just make the tree
            VG& tree = trees[entrypoint];
            components[entrypoint] = component;
            // maps from new ids to old ones
            map<id_t, pair<id_t, bool> >& trans = translations[entrypoint];
            // maps from old to new ids
            map<pair<id_t, bool>, set<id_t> >& itrans = inv_translations[entrypoint];

            /*
            // backtracking search

            procedure bt(c)
            if reject(P,c) then return
            if accept(P,c) then output(P,c)
            s ← first(P,c)
            while s ≠ Λ do
            bt(s)
            s ← next(P,s)
            */

            function<void(pair<id_t,bool>,id_t,bool,uint32_t,uint32_t)>
                bt = [&](pair<id_t, bool> curr,
                         id_t parent,
                         bool in_cycle,
                         uint32_t length,
                         uint32_t branches) {
                //cerr << "bt: curr: " << curr.first << " parent: " << parent << " len: " << length << " " << branches << " " << (in_cycle?"loop":"acyc") << endl;

                // i.   If the current node is outside the component,
                //       terminate this branch and return to the previous branching point.
                if (!component.count(curr.first)) {
                    return;
                }
                // ii.  Create a new copy of the current node in the DAG
                //       and use that copy for this branch.
                Node* node = get_node(curr.first);
                string curr_node_seq = node->sequence();
                // if we've reversed, take the reverse complement of the node we're flipping
                if (curr.second) curr_node_seq = reverse_complement(curr_node_seq);
                // use the forward orientation by default
                id_t cn = tree.create_node(curr_node_seq)->id();
                // record the mapping from the new id to the old
                trans[cn] = curr;
                // and record the inverse mapping
                itrans[curr].insert(cn);
                // and build the tree by connecting to our parent
                if (parent) {
                    tree.create_edge(parent, cn);
                }

                // iii. If we have found the first cycle in the current branch
                //       (the current node is the first one occurring twice in the path),
                //       initialize path length to 1 for this branch.
                //       (We need to find a k-path starting from the last offset of the previous node.)
                // walk the path back to the root to determine if we are the first cycling node
                id_t p = cn;
                // check is borked
                while (!in_cycle) { // && trans[p] != entrypoint) {
                    auto parents = tree.sides_to(NodeSide(p, false));
                    if (parents.size() < 1) { break; }
                    assert(parents.size() == 1);
                    p = parents.begin()->node;
                    // this node cycles
                    if (trans[p] == trans[cn]) {
                        in_cycle = true;
                        //length = 1; // XXX should this be 1? or maybe the length of the node?
                        break;
                    }
                }

                // iv.  If we have found a cycle in the current branch,
                //       increment path length by the length of the label of the current node.
                if (in_cycle) {
                    length += curr_node_seq.length();
                } else {
                    // if we branch here so we need to record it
                    // so as to avoid branching too many times before reaching a loop
                    auto s = start_degree(node);
                    auto e = end_degree(node);
                    branches += max(s-1 + e-1, 0);
                }

                // v.   If path length >= k, terminate this branch
                //       and return to the previous branching point.
                if (length >= max_length || (max_branch > 0 && branches >= max_branch)) {
                    return;
                } else {
                    // for each next node
                    if (!curr.second) {
                        // forward direction -- sides from end
                        for (auto& side : sides_from(node_end(curr.first))) {
                            // we enter the next node in the forward direction if the side is
                            // the start, and the reverse if it is the end
                            bt(make_pair(side.node, side.is_end), cn, in_cycle, length, branches);
                        }
                        // this handles the case of doubly-inverted edges (from start to end)
                        // which we can follow as if they are forward
                        // we stay on the same strand
                        for (auto& side : sides_to(node_end(curr.first))) {
                            // we enter the next node in the forward direction if the side coming
                            // from the start, and the reverse if it is the end
                            bt(make_pair(side.node, !side.is_end), cn, in_cycle, length, branches);
                        }
                    } else { // we've inverted already
                        // so we look at the things that leave the node on the reverse strand
                        // inverted
                        for (auto& side : sides_from(node_start(curr.first))) {
                            // here, sides from the start to the end maintain the flip
                            // but if we go to the start of the next node, we invert back
                            bt(make_pair(side.node, side.is_end), cn, in_cycle, length, branches);
                        }
                        // odd, but important quirk
                        // we also follow the normal edges, but in the reverse direction
                        // because we store them in the forward orientation
                        for (auto& side : sides_to(node_start(curr.first))) {
                            bt(make_pair(side.node, side.is_end), cn, in_cycle, length, branches);
                        }
                    }
                }

            };

            // we start with the entrypoint and run backtracking search
            bt(make_pair(entrypoint, false), 0, false, 0, 0);
        }
    }

    // -------------------------
    // tree -> dag conversion
    // -------------------------
    // now simplify each tree into a dag
    // algorithm sketch:
    // for each tree, we'll make a dag
    //   1) we start by labeling each node with its rank from the root (we have a tree, then DAGs)
    //      among nodes with the same original identity
    //   2) then, we pick the set of nodes that forms the largest group with the same identity
    //      and merge them
    //      if there is no group >1, exit, else goto (1), relabeling
    map<id_t, VG> dags;
    for (auto& g : trees) {
        id_t entrypoint = g.first;
        VG& tree = trees[entrypoint];
        VG& dag = dags[entrypoint];
        dag = tree; // copy
        map<id_t, pair<id_t, bool> >& trans = translations[entrypoint];
        map<pair<id_t, bool>, set<id_t> >& itrans = inv_translations[entrypoint];
        // rank among nodes with same original identity labeling procedure
        map<pair<id_t, bool>, size_t> orig_off;
        size_t i = 0;
        for (auto& j : itrans) {
            orig_off[j.first] = i;
            ++i;
        }
        // set up the initial vector we'll use when we have no inputs
        vector<uint32_t> zeros(orig_off.size(), 0);
        bool stable = false;
        uint16_t iter = 0;
        do {
            // -------------------------
            // first establish the rank of each node among other nodes with the same original identity
            // -------------------------
            // we'll store our positional rank vectors in the current here
            map<id_t, vector<uint32_t> > rankmap;
            // the graph is sorted (and will stay so)
            // as such we can run across it in sorted order
            dag.for_each_node([&](Node* n) {
                    // collect inbound vectors
                    vector<vector<uint32_t> > iv;
                    for (auto& side : dag.sides_to(n->id())) {
                        id_t in = side.node;
                        // should be satisfied by partial order property of DAG
                        assert(rankmap.find(in) != rankmap.end());
                        assert(!rankmap[in].empty());
                        iv.push_back(rankmap[in]);
                    }
                    // take their maximum
                    vector<uint32_t> ranks = (iv.empty() ? zeros : vpmax(iv));
                    // and increment this node so that the rankmap shows our ranks for each
                    // node in the original set at this point in the graph
                    ++ranks[orig_off[trans[n->id()]]];
                    // then save it in the rankmap
                    rankmap[n->id()] = ranks;
                });
            /*
            for (auto& r : rankmap) {
                cerr << r.first << ":" << dag.get_node(r.first)->sequence() << " [ ";
                for (auto c : r.second) {
                    cerr << c << " ";
                }
                cerr << "]" << endl;
            }
            */

            // -------------------------
            // now establish the class relative ranks for each node
            // -------------------------
            // maps from node in the current graph to the original identitiy and its rank among
            // nodes that also are clones of that node
            map<id_t, pair<pair<id_t, bool>, uint32_t> > rank_among_same;
            dag.for_each_node([&](Node* n) {
                    id_t id = n->id();
                    rank_among_same[id] = make_pair(trans[id], rankmap[id][orig_off[trans[id]]]);
                });
            // dump
            /*
            for (auto& r : rank_among_same) {
                cerr << r.first << ":" << dag.get_node(r.first)->sequence()
                     << " " << r.second.first.first << (r.second.first.second?"-":"+")
                     << ":" << r.second.second << endl;
            }
            */
            // groups
            // populate group sizes
            // groups map from the
            map<pair<pair<id_t, bool>, uint32_t>, vector<id_t> > groups;
            for (auto& r : rank_among_same) {
                groups[r.second].push_back(r.first);
            }
            // and find the biggest one
            map<uint16_t, vector<pair<pair<id_t, bool>, uint32_t> > > groups_by_size;
            for (auto& g : groups) {
                groups_by_size[g.second.size()].push_back(g.first);
            }
            // do we have a group that can be merged?
            if (groups_by_size.rbegin()->first > 1) {
                // -----------------------
                // merge the nodes that are the same and in the largest group
                // -----------------------
                auto orig = groups_by_size.rbegin()->second.front();
                auto& group = groups[orig];
                list<Node*> to_merge;
                for (auto id : group) {
                    to_merge.push_back(dag.get_node(id));
                }
                auto merged = dag.merge_nodes(to_merge);
                // we've now merged the redundant nodes
                // now we need to update the translations to reflect the fact
                // that we no longer have certain nodes in the dag
                id_t new_id = merged->id();
                // remove all old translations from new to old
                // and insert the new one
                // do the same for the inverted translations
                // from old to new
                set<id_t>& inv = itrans[orig.first];
                for (auto id : group) {
                    trans.erase(id);
                    inv.erase(id);
                }
                // store the new new->old translation
                trans[new_id] = orig.first;
                // and its inverse
                inv.insert(new_id);
            } else {
                // we have no group with more than one member
                // we are done
                stable = true;
            }
            // sort the graph
            dag.sort();
        } while (!stable);
    }

    // recover all the edges that link the nodes in the acyclic components of the graph
    unrolled.for_each_node([&](Node* n) {
            // get its edges in the original graph, and check if their endpoints are now
            // included in the unrolled graph (in which case they must be acyclic)
            // if they are, add the edge to the graph
            for (auto e : this->edges_of(get_node(n->id()))) {
                if (unrolled.has_node(e->from()) && unrolled.has_node(e->to())) {
                    unrolled.add_edge(*e);
                }
            }
        });


    // -----------------------------------------------
    // connect unrolled components back into the graph
    // -----------------------------------------------
    // now that we've unrolled all of the components
    // what do we do
    // we link them back into the rest of the graph in two steps

    // for each component we've unrolled (actually, for each entrypoint)
    for (auto& g : dags) {
        id_t entrypoint = g.first;
        VG& dag = dags[entrypoint];
        set<id_t> component = components[entrypoint];
        map<id_t, pair<id_t, bool> >& trans = translations[entrypoint];
        map<pair<id_t, bool>, set<id_t> >& itrans = inv_translations[entrypoint];

        // 1) increment the node ids to not conflict with the rest of the graph
        //       while recording the changes to the translation
        id_t max_id = max_node_id();
        // update the node ids in the dag
        dag.increment_node_ids(max_id);
        map<id_t, pair<id_t, bool> > trans_incr;
        // increment the translation from new node to old
        for (auto t = trans.begin(); t != trans.end(); ++t) {
            trans_incr[t->first + max_id] = t->second;
        }
        // save the incremented translation
        trans = trans_incr;
        // and do the same for the reverse mapping from old ids to new ones
        for (auto t = itrans.begin(); t != itrans.end(); ++t) {
            set<id_t> n;
            for (auto i = t->second.begin(); i != t->second.end(); ++i) {
                n.insert(*i + max_id);
            }
            t->second = n;
        }
        // 2) now that we don't conflict, add the component to the graph
        unrolled.extend(dag);
        // and also record the translation into the external reference
        for (auto t : trans) {
            // we map from a node id in the graph to an original node id and its orientation
            node_translation[t.first] = t.second;
        }

        // 3) find all the links into the component
        //       then connect all the nodes they
        //       now relate to using itrans
        for (auto& i : itrans) {
            auto old_id = i.first.first;
            auto is_flipped = i.first.second;
            set<id_t>& new_ids = i.second;
            // collect connections to old id that aren't in the component
            // we need to take the reverse complement of these if we are flipped relatively
            // sides to forward
            // add reverse complement function for edges and nodesides
            // make the connections
            for (auto i : new_ids) {
                // sides to forward
                for (auto& s : sides_to(NodeSide(old_id, false))) {
                    // if we are not in the component
                    if (!component.count(s.node)) {
                        // if we aren't flipped, we just make the connection
                        if (!is_flipped) {
                            unrolled.create_edge(s, NodeSide(i, false));
                        } else {
                            // the new node is flipped in unrolled relative to the other graph
                            // so we need to connect from the opposite side
                            unrolled.create_edge(s, NodeSide(i, true));
                        }
                    }
                }
                // sides to reverse
                for (auto& s : sides_to(NodeSide(old_id, true))) {
                    // if we are not in the component
                    if (!component.count(s.node)) {
                        if (!is_flipped) {
                            unrolled.create_edge(s, NodeSide(i, true));
                        } else {
                            unrolled.create_edge(s, NodeSide(i, false));
                        }
                    }
                }
                // sides from forward
                for (auto& s : sides_from(NodeSide(old_id, true))) {
                    // if we are not in the component
                    if (!component.count(s.node)) {
                        if (!is_flipped) {
                            unrolled.create_edge(NodeSide(i, true), s);
                        } else {
                            unrolled.create_edge(NodeSide(i, false), s);
                        }
                    }
                }
                // sides from reverse
                for (auto& s : sides_from(NodeSide(old_id, false))) {
                    // if we are not in the component
                    if (!component.count(s.node)) {
                        if (!is_flipped) {
                            unrolled.create_edge(NodeSide(i, false), s);
                        } else {
                            unrolled.create_edge(NodeSide(i, true), s);
                        }
                    }
                }
            }
        }
    }

    return unrolled;
}

void VG::orient_nodes_forward(set<id_t>& nodes_flipped) {
    // TODO: update paths in the graph when you do this!

    // Clear the flipped nodes set.
    nodes_flipped.clear();

    // First do the topological sort to order and orient
    deque<NodeTraversal> order_and_orientation;
    topological_sort(order_and_orientation);
#ifdef debug
#pragma omp critical (cerr)
    cerr << "+++++++++++++++++++++DOING REORIENTATION+++++++++++++++++++++++" << endl;
#endif

    // These are the node IDs we've visited so far
    set<id_t> visited;

    for(auto& traversal : order_and_orientation) {
        // Say we visited this node
#ifdef debug
#pragma omp critical (cerr)
        cerr << "Visiting " << traversal.node->id() << endl;
#endif
        visited.insert(traversal.node->id());

        // Make sure this node is the "from" in all its edges with un-visited nodes.

        if(traversal.backward) {
            // We need to flip this node around.
#ifdef debug
#pragma omp critical (cerr)
            cerr << "Flipped node " << traversal.node->id() << endl;
#endif
            // Say we flipped it
            nodes_flipped.insert(traversal.node->id());

            // Flip the sequence
            traversal.node->set_sequence(reverse_complement(traversal.node->sequence()));

        }

        // Get all the edges
        vector<Edge*> node_edges;
        edges_of_node(traversal.node, node_edges);

        // We need to unindex all the edges we're going to change before any of
        // them get re-indexed. Otherwise, we might have to try to hold two of
        // the same edge in the index at the same time, if one edge is fixed up
        // to be equal to another original edge before the other original edge
        // is fixed.

        // Edges that go from unvisited things to here need to be flipped from/to.
        vector<Edge*> edges_to_flip;
        copy_if(node_edges.begin(), node_edges.end(), back_inserter(edges_to_flip), [&](Edge* edge) {
            // We only need to flip the edges that are from things we haven't visited yet to here.
            return edge->to() == traversal.node->id() && visited.count(edge->from()) == 0;
        });

        for(Edge* edge : (traversal.backward ? node_edges : edges_to_flip)) {
            // Unindex every edge if we flipped the node, or only the edges we are flipping otherwise
            unindex_edge_by_node_sides(edge);

#ifdef debug
#pragma omp critical (cerr)
            cerr << "Unindexed edge " << edge->from() << (edge->from_start() ? " start" : " end")
                 << " -> " << edge->to() << (edge->to_end() ? " end" : " start") << endl;
#endif
        }

        for(Edge* edge : edges_to_flip) {
            // Flip around all the edges that need flipping
            // Flip the nodes
            id_t temp_id = edge->from();
            edge->set_from(edge->to());
            edge->set_to(temp_id);

            // Move the directionality flags, but invert both.
            bool temp_orientation = !edge->from_start();
            edge->set_from_start(!edge->to_end());
            edge->set_to_end(temp_orientation);
#ifdef debug
#pragma omp critical (cerr)
            cerr << "Reversed edge direction to " << edge->from() << (edge->from_start() ? " start" : " end")
                 << " -> " << edge->to() << (edge->to_end() ? " end" : " start") << endl;
#endif
        }

        if(traversal.backward) {
            for(Edge* edge : node_edges) {
                // Now that all the edges have the correct to and from, flip the
                // appropriate from_start and to_end flags for the end(s) on this
                // node, since we flipped the node.

                if(edge->to() == traversal.node->id()) {
                    edge->set_to_end(!edge->to_end());
                }
                if(edge->from() == traversal.node->id()) {
                    edge->set_from_start(!edge->from_start());
                }
            }
        }

        for(Edge* edge : (traversal.backward ? node_edges : edges_to_flip)) {
            // Reindex exactly what was unindexed
            index_edge_by_node_sides(edge);

#ifdef debug
#pragma omp critical (cerr)
            cerr << "Reindexed edge " << edge->from() << (edge->from_start() ? " start" : " end")
                 << " -> " << edge->to() << (edge->to_end() ? " end" : " start") << endl;
#endif
        }

        // It should always work out that the edges are from end to
        // start when we are done, but right now they might not be,
        // because the nodes at the other ends may still need to be
        // flipped themselves.
    }

    // We now know there are no edges from later nodes to earlier nodes. But
    // edges that are reversing, or "from" earlier nodes and "to" later nodes
    // with both end flags set, will cause problems. So remove those.
    // This works out to just clearing any edges with from_start or to_end set.
    // We also need to clear otherwise-normal-looking self loops here.
    vector<Edge*> to_remove;
    for_each_edge([&](Edge* e) {
        if(e->from_start() || e->to_end() || e->from() == e->to()) {
            // gssw won't know how to read this. Get rid of it.
            to_remove.push_back(e);
        }
    });

    for(Edge* edge : to_remove) {
#ifdef debug
#pragma omp critical (cerr)
        cerr << "Removed cycle edge " << edge->from() << "->" << edge->to() << endl;
#endif
        destroy_edge(edge);
    }

}

NodeSide node_start(id_t id) {
    return NodeSide(id, false);
}

NodeSide node_end(id_t id) {
    return NodeSide(id, true);
}

} // end namespace
