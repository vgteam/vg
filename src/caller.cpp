#include <cstdlib>
#include <stdexcept>
#include "json2pb.h"
#include "caller.hpp"
#include "stream.hpp"

using namespace std;

namespace vg {

const double Caller::Log_zero = (double)-1e100;

const char Caller::Default_default_quality = 30;
const int Caller::Default_min_aug_support = 2;

Caller::Caller(VG* graph,
               int default_quality,
               int min_aug_support):
    _graph(graph),
    _default_quality(default_quality),
    _min_aug_support(min_aug_support) {
    assert(_min_aug_support > 0);
    _max_id = _graph->max_node_id();
    _node_divider._max_id = &_max_id;
}

// delete contents of table
Caller::~Caller() {
    clear();
}

void Caller::clear() {
    _node_calls.clear();
    _node_supports.clear();
    _insert_calls.clear();
    _insert_supports.clear();
    _augmented_graph.clear();
    _node_divider.clear();
    _visited_nodes.clear();
    _called_edges.clear();
    _augmented_edges.clear();
    _inserted_nodes.clear();
}

void Caller::write_augmented_graph(ostream& out, bool json) {
    if (json) {
        _augmented_graph.graph.paths.to_graph(_augmented_graph.graph.graph);
        out << pb2json(_augmented_graph.graph.graph);
    } else {
        _augmented_graph.graph.serialize_to_ostream(out);
    }
}

void Caller::call_node_pileup(const NodePileup& pileup) {

    _node = _graph->get_node(pileup.node_id());
    assert(_node != NULL);
    assert(_node->sequence().length() == pileup.base_pileup_size());
    
    _node_calls.clear();
    _insert_calls.clear();
    string def_char = "-";
    _node_calls.assign(_node->sequence().length(), Genotype(def_char, def_char));
    _insert_calls.assign(_node->sequence().length(), Genotype(def_char, def_char));
    _node_supports.clear();
    _node_supports.assign(_node->sequence().length(), make_pair(
                              StrandSupport(), StrandSupport()));
    _insert_supports.clear();
    _insert_supports.assign(_node->sequence().length(), make_pair(
                                StrandSupport(), StrandSupport()));

    // process each base in pileup individually
    #pragma omp parallel for
    for (int i = 0; i < pileup.base_pileup_size(); ++i) {
        int num_inserts = 0;
        for (auto b : pileup.base_pileup(i).bases()) {
            if (b == '+') {
                ++num_inserts;
            }
        }
        int pileup_depth = max(num_inserts, pileup.base_pileup(i).num_bases() - num_inserts);
        if (pileup_depth >= 1) {
            call_base_pileup(pileup, i, false);
            call_base_pileup(pileup, i, true);
        }
    }

    // add nodes and edges created when making calls to the output graph
    // (_side_map gets updated)
    create_node_calls(pileup);

    _visited_nodes.insert(_node->id());
}

void Caller::call_edge_pileup(const EdgePileup& pileup) {
    if (pileup.num_reads() >= 1) {

        double qual_sum = 0;
        
        for (int i = 0; i < pileup.num_reads(); ++i) {
            char qual = !pileup.qualities().empty() ? pileup.qualities().at(i) : _default_quality;
            qual_sum += (double)qual;
        }

        Edge edge = pileup.edge(); // gcc not happy about passing directly
        _called_edges[NodeSide::pair_from_edge(edge)] = StrandSupport(
            pileup.num_forward_reads(),
            pileup.num_reads() - pileup.num_forward_reads(),
            qual_sum);
    }
}

void Caller::update_augmented_graph() {
    
    // Add nodes we don't think necessarily exist.
    function<void(Node*)> add_node = [&](Node* node) {
        if (_visited_nodes.find(node->id()) == _visited_nodes.end()) {
            Node* call_node = _augmented_graph.graph.create_node(node->sequence(), node->id());
            _node_divider.add_fragment(node, 0, call_node, NodeDivider::EntryCat::Ref,
                                       vector<StrandSupport>());
        }
    };
    _graph->for_each_node(add_node);

    // map every edge in the original graph to equivalent sides
    // in the call graph. if both sides exist, make an edge in the call graph
    function<void(Edge*)> map_edge = [&](Edge* edge) {
        pair<NodeSide, NodeSide> sides = NodeSide::pair_from_edge(edge);
        // skip uncalled edges if not writing augmented graph
        auto called_it = _called_edges.find(sides);
        bool called = called_it != _called_edges.end();
        
        StrandSupport support = called ? called_it->second : StrandSupport();
        assert(support.fs >= 0 && support.rs >= 0);
        
        Node* side1 = _graph->get_node(sides.first.node);
        Node* side2 = _graph->get_node(sides.second.node);
        // find up to two nodes matching side1 in the call graph
        int from_offset = !sides.first.is_end ? 0 : side1->sequence().length() - 1;
        int to_offset = sides.second.is_end ? side2->sequence().length() - 1 : 0;
        char cat = called ? 'R' : 'U';
        create_augmented_edge(side1, from_offset, !sides.first.is_end, true,
                              side2, to_offset, !sides.second.is_end, true, cat,
                              support);
    };            

    function<void(bool)> process_augmented_edges = [&](bool pass1) {
        for (auto& i : _augmented_edges) {
            auto& sides = i.first;
            char cat = i.second;
            NodeOffSide os1 = sides.first;
            Node* node1;
            bool aug1;
            if (_graph->has_node(os1.first.node)) {
                node1 = _graph->get_node(os1.first.node);
                aug1 = true;
            } else {
                // snp or isnert node -- need to get from call grpah
                // note : that we should never break these as they aren't in
                // the divider structure (will be caught down the road)
                node1 = _augmented_graph.graph.get_node(os1.first.node);
                aug1 = false;
            }
            int from_offset = os1.second;
            bool left1 = !os1.first.is_end;
            NodeOffSide os2 = sides.second;
            Node* node2;
            bool aug2;
            if (_graph->has_node(os2.first.node)) {
                node2 = _graph->get_node(os2.first.node);
                aug2 = true;
            } else {
                // snp or insert node -- need to get from call graph
                node2 = _augmented_graph.graph.get_node(os2.first.node);
                aug2 = false;
            }
            // only need to pass support for here insertions, other cases handled elsewhere
            StrandSupport support(-1, -1);
            if (cat == 'I') {
                auto ins_it = _inserted_nodes.find(os1.first.node);
                if (ins_it != _inserted_nodes.end()) {
                    support = ins_it->second.sup;
                } else {
                    ins_it = _inserted_nodes.find(os2.first.node);
                    assert(ins_it != _inserted_nodes.end());
                    support = ins_it->second.sup;
                }
                assert(support.fs >= 0 && support.rs >= 0);
            }
            int to_offset = os2.second;
            bool left2 = !os2.first.is_end;
            // todo: clean this up
            if (!pass1) {
              create_augmented_edge(node1, from_offset, left1, aug1,  node2, to_offset, left2, aug2, cat, support);
            } else {
                _node_divider.break_end(node1, &_augmented_graph.graph, from_offset, left1);
                _node_divider.break_end(node2, &_augmented_graph.graph, to_offset, left2);
            }
        }
    };

    // two passes here is a hack to make sure break_end is called on all edge ends
    // before processing any of them. 
    process_augmented_edges(true);
    _graph->for_each_edge(map_edge);
    process_augmented_edges(false);

    // Annotate all the nodes in the divider structure in the AugmentedGraph
    annotate_augmented_nd();
    // add on the inserted nodes
    for (auto i : _inserted_nodes) {
        auto& n = i.second; 
        annotate_augmented_node(n.node, 'I', n.sup, n.orig_id, n.orig_offset);
    }
}


void Caller::map_paths() {
    // We don't remove any nodes, so paths always stay connected
    function<void(const Path&)> lambda = [&](const Path& path) {
        list<Mapping>& call_path = _augmented_graph.graph.paths.create_path(path.name());
        int64_t last_rank = -1;
        int64_t last_call_rank = 0;
        size_t running_len = 0;
        size_t path_len = 0;
        for (int i = 0; i < path.mapping_size(); ++i) {
            const Mapping& mapping = path.mapping(i);
            int64_t rank = mapping.rank() == 0 ? i+1 : mapping.rank();
            size_t len = mapping.edit_size() == 0 ? _graph->get_node(mapping.position().node_id())->sequence().length() :
                mapping.edit(0).from_length();
            size_t to_len = mapping.edit_size() == 0 ? len : mapping.edit(0).to_length();
            if (mapping.edit_size() > 1 || len != to_len || rank <= last_rank) {
                cerr << "rank " << rank;
                cerr << "mapping " << pb2json(mapping) << endl;
                cerr << "last_rank " << last_rank << endl;
                cerr << "i " << i << endl;
                cerr << "Skipping input path " << path.name()
                     << " because ranks out of order or non-trivial edits." << endl;
                set<string> s;
                s.insert(path.name());
                _augmented_graph.graph.paths.remove_paths(s);
                return;
            }
            int64_t node_id = mapping.position().node_id();
            Node* node = _graph->get_node(node_id);
            
            int start = mapping.position().offset();
            if (mapping.position().is_reverse()) {
                start = node->sequence().length() - 1 - start;
            }
            int end = mapping.position().is_reverse() ? start - len + 1 : start + len - 1;
            list<Mapping> call_mappings = _node_divider.map_node(node_id, start, len,
                                                                 mapping.position().is_reverse());
            for (auto& cm : call_mappings) {
                running_len += cm.edit(0).from_length();
                cm.set_rank(++last_call_rank);
                call_path.push_back(cm);
            }
            path_len += len;
            last_rank = rank; 
        }
        assert(running_len == path_len);
        verify_path(path, call_path);
    };
    _graph->paths.for_each(lambda);

    // make sure paths are saved
    _augmented_graph.graph.paths.rebuild_node_mapping();
    _augmented_graph.graph.paths.rebuild_mapping_aux();
    _augmented_graph.graph.paths.to_graph(_augmented_graph.graph.graph);    
}

void Caller::verify_path(const Path& in_path, const list<Mapping>& call_path) {
    function<string(VG*, const Mapping&)> lambda = [](VG* graph, const Mapping& mapping) {
        const Position& pos = mapping.position();
        const Node* node = graph->get_node(pos.node_id());
        string s = node->sequence();
        int64_t offset = pos.offset();
        if (mapping.position().is_reverse()) {
            s = reverse_complement(s);
        }
        if (offset > 0) {
            s = s.substr(pos.offset());
        }
        return s;
    };

    string in_string;
    for (int i = 0; i < in_path.mapping_size(); ++i) {
        in_string += lambda(_graph, in_path.mapping(i));
    }
    string call_string;
    for (auto& m : call_path) {
        call_string += lambda(&_augmented_graph.graph, m);
    }

    assert(in_string == call_string);

}

void Caller::create_augmented_edge(Node* node1, int from_offset, bool left_side1, bool aug1,
                                   Node* node2, int to_offset, bool left_side2, bool aug2, char cat,
                                   StrandSupport support) {
    NodeDivider::Entry call_sides1;
    NodeDivider::Entry call_sides2;

    if (aug1) {
        call_sides1 = _node_divider.break_end(node1, &_augmented_graph.graph, from_offset,
                                              left_side1);
    } else {
        call_sides1 = NodeDivider::Entry(node1, vector<StrandSupport>(1, support));
    }
    if (aug2) {
        call_sides2 = _node_divider.break_end(node2, &_augmented_graph.graph, to_offset,
                                              left_side2);
    } else {
        call_sides2 = NodeDivider::Entry(node2, vector<StrandSupport>(1, support));
    }
    
    // make up to 9 edges connecting them in the call graph
    for (int i = 0; i < (int)NodeDivider::EntryCat::Last; ++i) {
        for (int j = 0; j < (int)NodeDivider::EntryCat::Last; ++j) {
            if (call_sides1[i] != NULL && call_sides2[j] != NULL) {
                // always make links between alts and reference
                // (be more strict on deletion edges, only linking two reference)
                bool link_edge = ((i == (int)NodeDivider::EntryCat::Ref &&
                                  j == (int)NodeDivider::EntryCat::Ref) ||
                                  ((i == (int)NodeDivider::EntryCat::Ref || 
                                    j == (int)NodeDivider::EntryCat::Ref) &&
                                   cat != 'L'));

                NodeSide side1(call_sides1[i]->id(), !left_side1);
                NodeSide side2(call_sides2[j]->id(), !left_side2);
                if (!_augmented_graph.graph.has_edge(side1, side2)) {
                  StrandSupport edge_support = support >= StrandSupport() ? support :
                     min(avgSup(call_sides1.sup(i)), avgSup(call_sides2.sup(j)));

                  NodeOffSide no1(NodeSide(node1->id(), !left_side1), from_offset);
                  NodeOffSide no2(NodeSide(node2->id(), !left_side2), to_offset);
                  // take augmented deletion edge support from the pileup
                  if (cat == 'L') {
                    edge_support = _deletion_supports[minmax(no1, no2)];
                  }
                  // hack to decrease support for an edge that spans an insertion, by subtracting
                  // that insertion's copy number.  
                  auto is_it = _insertion_supports.find(minmax(no1, no2));
                  if (is_it != _insertion_supports.end()) {
                    edge_support = edge_support - is_it->second;
                  }

                  if (link_edge || edge_support.total() >= _min_aug_support) {
                        
                    Edge* edge = _augmented_graph.graph.create_edge(call_sides1[i], call_sides2[j],
                                                                    left_side1, !left_side2);

                    // TODO: can edges be annotated more than once with
                    // different cats? if so, last one will prevail. should
                    // check if this can impact vcf converter...
                    annotate_augmented_edge(edge, cat, edge_support);
                  }
                }
            }
        }
    }
}

void Caller::call_base_pileup(const NodePileup& np, int64_t offset, bool insertion) {
    const BasePileup& bp = np.base_pileup(offset);

    // parse the pilueup structure
    vector<pair<int64_t, int64_t> > base_offsets;
    Pileups::parse_base_offsets(bp, base_offsets);

    // compute top two most frequent bases and their counts
    string top_base;
    int top_count;
    int top_rev_count;
    string second_base;
    int second_count;
    int second_rev_count;
    int total_count;
    compute_top_frequencies(bp, base_offsets, top_base, top_count, top_rev_count,
                            second_base, second_count, second_rev_count, total_count,
                            insertion);

    // note first and second base will be upper case too
    string ref_base = string(1, ::toupper(bp.ref_base()));

    // get references to node-level members we want to update
    Genotype& base_call = insertion ? _insert_calls[offset] : _node_calls[offset];
    pair<StrandSupport, StrandSupport>& support = insertion ? _insert_supports[offset] : _node_supports[offset];

    if (top_count >= _min_aug_support || (top_base == ref_base && top_count > 0)) {
        base_call.first = top_base != ref_base ? top_base : ".";
        support.first.fs = top_count - top_rev_count;
        support.first.rs = top_rev_count;
        support.first.qual = total_base_quality(bp, base_offsets, top_base);
    }
    if (second_count >= _min_aug_support || (second_base == ref_base && second_count > 0)) { 
        base_call.second = second_base != ref_base ? second_base : ".";
        support.second.fs = second_count - second_rev_count;
        support.second.rs = second_rev_count;
        support.second.qual = total_base_quality(bp, base_offsets, second_base);
    }
}

void Caller::compute_top_frequencies(const BasePileup& bp,
                                     const vector<pair<int64_t, int64_t> >& base_offsets,
                                     string& top_base, int& top_count, int& top_rev_count,
                                     string& second_base, int& second_count, int& second_rev_count,
                                     int& total_count, bool inserts) {

    // histogram of pileup entries (base, indel)
    unordered_map<string, int> hist;
    // same thing but just reverse strand (used for strand bias filter)
    unordered_map<string, int> rev_hist;

    total_count = 0;
    const string& bases = bp.bases();
    string ref_base = string(1, ::toupper(bp.ref_base()));
    
    // compute histogram from pileup
    for (auto i : base_offsets) {
        string val = Pileups::extract(bp, i.first);

        if ((inserts && val[0] != '+') || (!inserts && val[0] == '+')) {
            // toggle inserts
            continue;
        }
        
        // We want to know if this pileup supports an N
        if (is_all_n(val)) {
            // N is not a real base, so we should never augment with it.
            continue;
        }
        
        ++total_count;
        
        // val will always be uppcase / forward strand.  we check
        // the original pileup to see if reversed
        bool reverse = bases[i.first] == ',' ||
            (bases[i.first] == '+' && ::islower(bases[i.first + val.length() - 1])) ||
            (bases[i.first] != '-' && ::islower(bases[i.first]));
        if (bases[i.first] == '-') {
            string tok = Pileups::extract(bp, i.first);
            bool is_reverse, from_start, to_end;
            int64_t from_id, from_offset, to_id, to_offset;
            Pileups::parse_delete(tok, is_reverse, from_id, from_offset, from_start, to_id, to_offset, to_end);
            reverse = is_reverse;
            // reset reverse to forward
            if (is_reverse) {
                Pileups::make_delete(val, false, from_id, from_offset, from_start, to_id, to_offset, to_end);
            }
        }

        if (hist.find(val) == hist.end()) {
            hist[val] = 0;
            rev_hist[val] = 0;
        }
        ++hist[val];

        if (reverse) {
            ++rev_hist[val];
        }
    }

    // tie-breaker heuristic:
    // reference > transition > transversion > delete > insert > N
    function<int(const string&)> base_priority = [&ref_base](const string& base) {
        size_t n = base.length();
        if(n == 0) {
            // base == '' -> Uncalled: 0 Points
            return 0;
        }
        else
        {
             char cbase = base[0];
            
             if(n == 1)
             {
                switch(cbase) {
                    case '.':// Ref: 6 Points
                        return 6;
                    break;
                    case '-': // Uncalled: 0 Points
                        return 0;
                    break;
                    
                    // Transition: 5 points.  Transversion: 4 points 
                    case 'A':
                    case 't':
                        return cbase == 'G' ? 5 : 4;
                    break;
                    
                    case 'C':
                    case 'g':
                        return cbase == 'T' ? 5 : 4;
                    break;
                    
                    case 'G':
                    case 'c':
                        return cbase == 'A' ? 5 : 4;
                    break;
                    
                    case 'T':
                    case 'a':
                        return cbase == 'C' ? 5 : 4;
                    break;
                }
            }
            
            // need to happen in any other case - i.e. also for n > 1
            switch(cbase) {
                case '-':
                    return 3; // Deletion: 3 Points
                break;
                case '+':
                    return 2; // Insertion: 2 Points
                break;
            }
        }
        return 1;
    };

    // compare to pileup entries, to see which has greater count, use tie breaker logic
    // if count is the same
    function<bool(const string&, int, const string&, int)> base_greater = [&base_priority] (
        const string& base1, int count1, const string& base2, int count2) {
        if (count1 == count2) {
            int p1 = base_priority(base1);
            int p2 = base_priority(base2);
            if (p1 == p2) {
                return base1 > base2;
            } else {
                return p1 > p2;
            }
        }
        return count1 > count2;
    };
        
    // find highest occurring string
    top_base.clear();
    top_count = 0;
    for (auto i : hist) {
        if (base_greater(i.first, i.second, top_base, top_count)) {
            top_base = i.first;
            top_count = i.second;
        }
    }

    // find second highest occurring string
    // todo: do it in same pass as above
    second_base.clear();
    second_count = 0;
    for (auto i : hist) {
        if (i.first != top_base &&
            base_greater(i.first, i.second, second_base, second_count)) {
            second_base = i.first;
            second_count = i.second;
        }
    }
    assert(top_base == "" || top_base != second_base);
    top_rev_count = rev_hist[top_base];
    second_rev_count = rev_hist[second_base];
}

double Caller::total_base_quality(const BasePileup& bp,
                                  const vector<pair<int64_t, int64_t> >& base_offsets,
                                  const string& val) {
    double qual_sum = 0;

    const string& bases = bp.bases();
    const string& quals = bp.qualities();

    for (int i = 0; i < base_offsets.size(); ++i) {
        string base = Pileups::extract(bp, base_offsets[i].first);

        // make sure deletes always compared without is_reverse flag
        if (base.length() > 1 && base[0] == '-') {
            bool is_reverse, from_start, to_end;
            int64_t from_id, from_offset, to_id, to_offset;
            Pileups::parse_delete(base, is_reverse, from_id, from_offset, from_start, to_id, to_offset, to_end);
            // reset reverse to forward
            if (is_reverse) {
                Pileups::make_delete(base, false, from_id, from_offset, from_start, to_id, to_offset, to_end);
            }
        }

        if (base == val) {
            char qual = base_offsets[i].second >= 0 ? quals[base_offsets[i].second] : _default_quality;
            qual_sum += (double)qual;
        }
    }
    
    return qual_sum;
}

// please refactor me! 
void Caller::create_node_calls(const NodePileup& np) {
    
    int n = _node->sequence().length();
    const string& seq = _node->sequence();
    int cur = 0;
    int cat = call_cat(_node_calls[cur]);

    // scan calls, merging contiguous reference calls.  only consider
    // ref / snp / inserts on first pass.  
    // scan contiguous chunks of a node with same call
    // (note: snps will always be 1-base -- never merged)
    for (int next = 1; next <= n; ++next) {
        int next_cat = next == n ? -1 : call_cat(_node_calls[next]);

        // for anything but case where we merge consec. ref/refs
        if (cat == 2 || cat != next_cat ||
            _insert_calls[next-1].first[0] == '+' || _insert_calls[next-1].second[0] == '+') {

            if (cat == 0 || cat == 1) {
                // add reference
                vector<StrandSupport> sup;
                if (_node_calls[cur].first == ".") {
                    for (int i = cur; i < next; ++i) {
                        sup.push_back(_node_supports[i].first);
                    }
                }
                if (_node_calls[cur].second == ".") {
                    assert (_node_calls[cur].first != ".");
                    for (int i = cur; i < next; ++i) {
                        sup.push_back(_node_supports[i].second);
                    }
                }
                string new_seq = seq.substr(cur, next - cur);
                Node* node = _augmented_graph.graph.create_node(new_seq, ++_max_id);
                _node_divider.add_fragment(_node, cur, node, NodeDivider::EntryCat::Ref, sup);
                // bridge to node
                NodeOffSide no1(NodeSide(_node->id(), true), cur-1);
                NodeOffSide no2(NodeSide(_node->id(), false), cur);
                _augmented_edges[make_pair(no1, no2)] = 'R';
                // bridge from node
                no1 = NodeOffSide(NodeSide(_node->id(), true), next-1);
                no2 = NodeOffSide(NodeSide(_node->id(), false), next);
                _augmented_edges[make_pair(no1, no2)] = 'R';
            }            
            else {
                // some mix of reference and alts
                assert(next == cur + 1);
                
                function<void(string&, StrandSupport, string&, NodeDivider::EntryCat)>  call_het =
                    [&](string& call1, StrandSupport support1, string& call2, NodeDivider::EntryCat altCat) {
                
                    if (call1 == "." || (altCat == NodeDivider::EntryCat::Alt1 && call2 != ".")) {
                        // reference base
                        StrandSupport sup = call1 == "." ? support1 : StrandSupport();
                        assert(call2 != "."); // should be handled above
                        string new_seq = seq.substr(cur, 1);
                        Node* node = _augmented_graph.graph.create_node(new_seq, ++_max_id);
                        _node_divider.add_fragment(_node, cur, node, NodeDivider::EntryCat::Ref,
                                                   vector<StrandSupport>(1, sup));
                        // bridge to node
                        NodeOffSide no1(NodeSide(_node->id(), true), cur-1);
                        NodeOffSide no2(NodeSide(_node->id(), false), cur);
                        _augmented_edges[make_pair(no1, no2)] = 'R';
                        // bridge from node
                        no1 = NodeOffSide(NodeSide(_node->id(), true), next-1);
                        no2 = NodeOffSide(NodeSide(_node->id(), false), next);
                        _augmented_edges[make_pair(no1, no2)] = 'R';
                    }
                    if (call1 != "." && call1[0] != '-' && call1[0] != '+' && (
                            // we only want to process a homozygous snp once:
                            call1 != call2 || altCat == NodeDivider::EntryCat::Alt1)) {
                        StrandSupport sup = support1;
                        // snp base
                        string new_seq = call1;
                        Node* node = _augmented_graph.graph.create_node(new_seq, ++_max_id);
                        _node_divider.add_fragment(_node, cur, node, altCat,
                                                   vector<StrandSupport>(1, sup));
                        // bridge to node
                        NodeOffSide no1(NodeSide(_node->id(), true), cur-1);
                        NodeOffSide no2(NodeSide(_node->id(), false), cur);
                        _augmented_edges[make_pair(no1, no2)] = 'S';
                        // bridge from node
                        no1 = NodeOffSide(NodeSide(_node->id(), true), next-1);
                        no2 = NodeOffSide(NodeSide(_node->id(), false), next);
                        _augmented_edges[make_pair(no1, no2)] = 'S';
                    }
                    else if (call1 != "." && call1[0] == '-' && call1.length() > 1 && (
                                 // we only want to process homozygous delete once
                                 call1 != call2 || altCat == NodeDivider::EntryCat::Alt1)) {
                        // delete
                        int64_t del_len;
                        bool from_start;
                        int64_t from_id;
                        int64_t from_offset;
                        int64_t to_id;
                        int64_t to_offset;
                        bool to_end;
                        bool reverse;
                        Pileups::parse_delete(call1, reverse, from_id, from_offset, from_start, to_id, to_offset, to_end);
                        NodeOffSide s1(NodeSide(from_id, !from_start), from_offset);
                        NodeOffSide s2(NodeSide(to_id, to_end), to_offset);
                        Node* node1 = _graph->get_node(from_id);
                        assert(from_offset >=0 && from_offset < node1->sequence().length());
                        Node* node2 = _graph->get_node(to_id);
                        assert(to_offset >=0 && to_offset < node2->sequence().length());
                        
                        // we're just going to update the divider here, since all
                        // edges get done at the end
                        _augmented_edges[make_pair(s1, s2)] = 'L';
                        // keep track of its support
                        _deletion_supports[minmax(s1, s2)] = support1;
                        
                        // also need to bridge any fragments created above
                        if ((from_start && from_offset > 0) ||
                            (!from_start && from_offset < node1->sequence().length() - 1)) {
                            NodeOffSide no1(NodeSide(from_id, !from_start), from_offset);
                            NodeOffSide no2(NodeSide(from_id, from_start),
                                            (from_start ? from_offset - 1 : from_offset + 1));
                            if (_augmented_edges.find(make_pair(no1, no2)) == _augmented_edges.end()) {
                                _augmented_edges[make_pair(no1, no2)] = 'R';
                            }
                        }
                        if ((!to_end && to_offset > 0) ||
                            (to_end && to_offset < node2->sequence().length() - 1)) {
                            NodeOffSide no1(NodeSide(to_id, to_end), to_offset);
                            NodeOffSide no2(NodeSide(to_id, !to_end), !to_end ? to_offset - 1 : to_offset + 1);
                            if (_augmented_edges.find(make_pair(no1, no2)) == _augmented_edges.end()) {
                                _augmented_edges[make_pair(no1, no2)] = 'R';
                            }

                        }
                    }
                };

                // apply same logic to both calls, updating opposite arrays
                call_het(_node_calls[cur].first, _node_supports[cur].first,
                         _node_calls[cur].second, NodeDivider::EntryCat::Alt1);
                call_het(_node_calls[cur].second, _node_supports[cur].second,
                         _node_calls[cur].first, NodeDivider::EntryCat::Alt2);                
            }

            // inserts done separate at end since they take start between cur and next
            function<void(string&, StrandSupport, string&, StrandSupport, NodeDivider::EntryCat)>  call_inserts =
                [&](string& ins_call1, StrandSupport ins_support1, string& ins_call2, StrandSupport ins_support2,
                    NodeDivider::EntryCat altCat) {
                if (ins_call1[0] == '+' && (
                        // we only want to process homozygous insert once
                        ins_call1 != ins_call2 || altCat == NodeDivider::EntryCat::Alt1)) {
                    int64_t ins_len;
                    string ins_seq;
                    bool ins_rev;
                    Pileups::parse_insert(ins_call1, ins_len, ins_seq, ins_rev);
                    // todo: check reverse?
                    Node* node = _augmented_graph.graph.create_node(ins_seq, ++_max_id);
                    StrandSupport sup = ins_support1;
                    InsertionRecord ins_rec = {node, sup, _node->id(), next-1};
                    _inserted_nodes[node->id()] = ins_rec;

                    // bridge to insert
                    NodeOffSide no1(NodeSide(_node->id(), true), next-1);
                    NodeOffSide no2(NodeSide(node->id(), false), 0);
                    _augmented_edges[make_pair(no1, no2)] = 'I';
                    // bridge from insert
                    if (next < _node->sequence().length()) {
                        NodeOffSide no3 = NodeOffSide(NodeSide(node->id(), true), node->sequence().length() - 1);
                        NodeOffSide no4 = NodeOffSide(NodeSide(_node->id(), false), next);
                        _augmented_edges[make_pair(no3, no4)] = 'I';
                        // bridge across insert
                        _augmented_edges[make_pair(no1, no4)] = 'R';
                        // remember support "lost" to insertion so we
                        // can subtract it from the bridge later on
                        if (_insertion_supports.count(minmax(no1, no4))) {
                          _insertion_supports[minmax(no1, no4)] += sup;
                        } else {
                          _insertion_supports[minmax(no1, no4)] = sup;
                        }                                                    
                    } else {
                        // we have to link all outgoing edges to our insert if
                        // we're at end of node (unlike snps, the fragment structure doesn't
                        // handle these cases)
                        vector<pair<id_t, bool>> next_nodes = _graph->edges_end(_node->id());
                        NodeOffSide no3 = NodeOffSide(NodeSide(node->id(), true), node->sequence().length() - 1);
                        for (auto nn : next_nodes) {
                            int64_t offset4 = !nn.second ? 0 : _graph->get_node(nn.first)->sequence().length() - 1;
                            NodeOffSide no4 = NodeOffSide(NodeSide(nn.first, nn.second), offset4);
                            _augmented_edges[make_pair(no3, no4)] = 'I';
                            // bridge across insert
                            _augmented_edges[make_pair(no1, no4)] = 'R';
                            // remember support "lost" to insertion so we
                            // can subtract it from the bridge later on
                            if (_insertion_supports.count(minmax(no1, no4))) {
                                _insertion_supports[minmax(no1, no4)] += sup;
                            } else {
                                _insertion_supports[minmax(no1, no4)] = sup;
                            }                            
                        }
                    }
                }
            };
            
            call_inserts(_insert_calls[next-1].first, _insert_supports[next-1].first,
                         _insert_calls[next-1].second, _insert_supports[next-1].second,
                         NodeDivider::EntryCat::Alt1);
            call_inserts(_insert_calls[next-1].second, _insert_supports[next-1].second,
                         _insert_calls[next-1].first, _insert_supports[next-1].first,
                         NodeDivider::EntryCat::Alt2);
                
            // shift right
            cur = next;
            cat = next_cat;
        }
    }
}

void AugmentedGraph::clear() {
    // Reset to default state
    *this = AugmentedGraph();
}

void Caller::annotate_augmented_node(Node* node, char call, StrandSupport support, int64_t orig_id, int orig_offset)
{
    _augmented_graph.node_calls[node] = (ElementCall) call;
    _augmented_graph.node_supports[node].set_forward(support.fs);
    _augmented_graph.node_supports[node].set_reverse(support.rs);
    _augmented_graph.node_supports[node].set_quality(support.qual);
    
    if (orig_id != 0 && call != 'S' && call != 'I') {
        // Add translations for preserved parts
        Translation trans;
        auto* new_mapping = trans.mutable_to()->add_mapping();
        new_mapping->mutable_position()->set_node_id(node->id());
        auto* new_edit = new_mapping->add_edit();
        new_edit->set_from_length(node->sequence().size());
        new_edit->set_to_length(node->sequence().size());
        auto* old_mapping = trans.mutable_from()->add_mapping();
        old_mapping->mutable_position()->set_node_id(orig_id);
        old_mapping->mutable_position()->set_offset(orig_offset);
        auto* old_edit = old_mapping->add_edit();
        old_edit->set_from_length(node->sequence().size());
        old_edit->set_to_length(node->sequence().size());
        
        _augmented_graph.translations.push_back(trans);
    }
}

void Caller::annotate_augmented_edge(Edge* edge, char call, StrandSupport support)
{
    _augmented_graph.edge_calls[edge] = (ElementCall) call;
    _augmented_graph.edge_supports[edge].set_forward(support.fs);
    _augmented_graph.edge_supports[edge].set_reverse(support.rs);
    _augmented_graph.edge_supports[edge].set_quality(support.qual);
}

void Caller::annotate_augmented_nd()
{
    for (auto& i : _node_divider.index) {
        int64_t orig_node_id = i.first;
        for (auto& j : i.second) {
            int64_t orig_node_offset = j.first;
            NodeDivider::Entry& entry = j.second;
            char call = entry.sup_ref.empty() || maxSup(entry.sup_ref) == StrandSupport() ? 'U' : 'R';
            annotate_augmented_node(entry.ref, call, maxSup(entry.sup_ref), orig_node_id, orig_node_offset);
            if (entry.alt1 != NULL) {
                annotate_augmented_node(entry.alt1, 'S', maxSup(entry.sup_alt1), orig_node_id, orig_node_offset);
            }
            if (entry.alt2 != NULL) {
                annotate_augmented_node(entry.alt2, 'S', maxSup(entry.sup_alt2), orig_node_id, orig_node_offset);
            }
        }
    }
}

void NodeDivider::add_fragment(const Node* orig_node, int offset, Node* fragment,
                               EntryCat cat, vector<StrandSupport> sup) {
    
    NodeHash::iterator i = index.find(orig_node->id());
    if (i == index.end()) {
        i = index.insert(make_pair(orig_node->id(), NodeMap())).first;
    }

    NodeMap& node_map = i->second;
    NodeMap::iterator j = node_map.find(offset);
    
    if (j != node_map.end()) {
        assert(j->second[cat] == NULL);
        j->second[cat] = fragment;
        j->second.sup(cat) = sup;
    } else {
        Entry ins_triple;
        ins_triple[cat] = fragment;
        ins_triple.sup(cat) = sup;
        j = node_map.insert(make_pair(offset, ins_triple)).first;
    }
    // sanity checks to make sure we don't introduce an overlap
    if (offset == 0) {
        assert(j == node_map.begin());
    } 
    if (offset + fragment->sequence().length() == orig_node->sequence().length()) {
        assert(j == --node_map.end());
    } else if (j != --node_map.end()) {
        NodeMap::iterator next = j;
        ++next;
        assert(offset + fragment->sequence().length() <= next->first);
    }
}
            
NodeDivider::Entry NodeDivider::break_end(const Node* orig_node, VG* graph, int offset, bool left_side) {
    NodeHash::iterator i = index.find(orig_node->id());
    if (i == index.end()) {
        return Entry();
    }
    NodeMap& node_map = i->second;
    NodeMap::iterator j = node_map.upper_bound(offset);
    if (j == node_map.begin()) {
        return Entry();
    }

    --j;
    int sub_offset = j->first;

    function<pair<Node*, vector<StrandSupport>>(Node*, EntryCat, vector<StrandSupport>& )>  lambda =
        [&](Node* fragment, EntryCat cat, vector<StrandSupport>& sup) {
        if (offset < sub_offset || offset >= sub_offset + fragment->sequence().length()) {
            return make_pair((Node*)NULL, vector<StrandSupport>());
        }

        // if our cut point is already the exact left or right side of the node, then
        // we don't have anything to do than return it.
        if (offset == sub_offset && left_side == true) {
            return make_pair(fragment, sup);
        }
        if (offset == sub_offset + fragment->sequence().length() - 1 && left_side == false) {
            return make_pair(fragment, sup);
        }
        
        // otherwise, we're somewhere in the middle, and have to subdivide the node
        // first, shorten the exsisting node
        int new_len = left_side ? offset - sub_offset : offset - sub_offset + 1;
        assert(new_len > 0 && new_len != fragment->sequence().length());
        string frag_seq = fragment->sequence();
        *fragment->mutable_sequence() = frag_seq.substr(0, new_len);

        // then make a new node for the right part
        Node* new_node = graph->create_node(frag_seq.substr(new_len, frag_seq.length() - new_len), ++(*_max_id));
        
        // now divide up the support, starting with the right bit
        vector<StrandSupport> new_sup;
        if (!sup.empty()) {
            new_sup = vector<StrandSupport>(sup.begin() + new_len, sup.end());
            // then cut the input (left bit) in place
            sup.resize(new_len);
        }

        // update the data structure with the new node
        add_fragment(orig_node, sub_offset + new_len, new_node, cat, new_sup);

        return make_pair(new_node, new_sup);
    };

    vector<StrandSupport>& sup_ref = j->second.sup_ref;
    Node* fragment_ref = j->second.ref;
    auto new_node_info = fragment_ref != NULL ? lambda(fragment_ref, Ref, sup_ref) :
        make_pair((Node*)NULL, vector<StrandSupport>());
    
    vector<StrandSupport>& sup_alt1 = j->second.sup_alt1;
    Node* fragment_alt1 = j->second.alt1;
    auto new_node_alt1_info = fragment_alt1 != NULL ? lambda(fragment_alt1, Alt1, sup_alt1) :
        make_pair((Node*)NULL, vector<StrandSupport>());
    
    vector<StrandSupport>& sup_alt2 = j->second.sup_alt2;
    Node* fragment_alt2 = j->second.alt2;
    auto new_node_alt2_info = fragment_alt2 != NULL ? lambda(fragment_alt2, Alt2, sup_alt2) :
        make_pair((Node*)NULL, vector<StrandSupport>());

    Entry ret = left_side ? Entry(new_node_info.first, new_node_info.second,
                                  new_node_alt1_info.first, new_node_alt1_info.second,
                                  new_node_alt2_info.first, new_node_alt2_info.second) :
        Entry(fragment_ref, sup_ref, fragment_alt1, sup_alt1, fragment_alt2, sup_alt2);

    return ret;
}

// this function only works if node is completely covered in divider structure,
list<Mapping> NodeDivider::map_node(int64_t node_id, int64_t start_offset, int64_t length, bool reverse){
    NodeHash::iterator i = index.find(node_id);
    assert(i != index.end());
    NodeMap& node_map = i->second;
    assert(!node_map.empty());
    list<Mapping> out_mappings;
    int cur_len = 0;
    if (!reverse) {
        for (auto i : node_map) {
            if (i.first >= start_offset && cur_len < length) {
                Node* call_node = i.second.ref;
                assert(call_node != NULL);
                Mapping mapping;
                mapping.mutable_position()->set_node_id(call_node->id());
                if (start_offset > i.first && out_mappings.empty()) {
                    mapping.mutable_position()->set_offset(start_offset - i.first);
                } else {
                    mapping.mutable_position()->set_offset(0);
                }
                int map_len = call_node->sequence().length() - mapping.position().offset();
                if (map_len + cur_len > length) {
                    map_len = length - cur_len;
                }
                Edit* edit = mapping.add_edit();
                edit->set_from_length(map_len);
                edit->set_to_length(map_len);
                cur_len += map_len;
                out_mappings.push_back(mapping);
            }
        }
    } else {
        // should fold into above when on less cold meds. 
        for (NodeMap::reverse_iterator i = node_map.rbegin(); i != node_map.rend(); ++i)
        {
            if (i->first <= start_offset && cur_len < length) {
                Node* call_node = i->second.ref;
                Mapping mapping;
                mapping.mutable_position()->set_is_reverse(true);
                mapping.mutable_position()->set_node_id(call_node->id());
                if (start_offset >= i->first && start_offset < i->first + call_node->sequence().length() - 1) {
                    assert(out_mappings.empty());
                    mapping.mutable_position()->set_offset(start_offset - i->first);
                } else {
                    mapping.mutable_position()->set_offset(call_node->sequence().length() - 1);
                }
                int map_len = mapping.position().offset() + 1;
                if (map_len + cur_len > length) {
                    map_len = length - cur_len;
                }
                // switch up to new-style offset (todo: revise whole function)
                mapping.mutable_position()->set_offset(call_node->sequence().length() - 1 -
                                                       mapping.position().offset());
                assert(map_len <= call_node->sequence().length());
                assert(mapping.position().offset() >= 0 &&
                       mapping.position().offset() < call_node->sequence().length());
                Edit* edit = mapping.add_edit();
                edit->set_from_length(map_len);
                edit->set_to_length(map_len);
                cur_len += map_len;
                out_mappings.push_back(mapping);
            }
        }
    }
                
    assert(cur_len == length);
    return out_mappings;
}

void NodeDivider::clear() {
    index.clear();
}

ostream& operator<<(ostream& os, const NodeDivider::NodeMap& nm) {
    for (auto& x : nm) {
        os << x.first << "[" << (x.second.ref ? x.second.ref->id() : -1) << ", "
           << (x.second.alt1 ? x.second.alt1->id() : -1) << ", "
           << (x.second.alt2 ? x.second.alt2->id() : -1) << "]" << endl;
    }
    return os;
}

ostream& operator<<(ostream& os, NodeDivider::Entry entry) {
    for (int i = 0; i < NodeDivider::EntryCat::Last; ++i) {
        if (entry[i] != NULL) {
            os << pb2json(*entry[i]);
        }
        else {
            os << "NULL";
        }
        os << ", ";
    }

    return os;
}

ostream& operator<<(ostream& os, const Caller::NodeOffSide& no) {
    os << "NOS(" << no.first.node << ":" << no.second << ",left=" << !no.first.is_end << ")";
    return os;
}

}
