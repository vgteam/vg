#include <cstdlib>
#include <stdexcept>
#include <regex>
#include "json2pb.h"
#include "pileup.hpp"
#include "stream.hpp"

using namespace std;

namespace vg {

void Pileups::clear() {
    for (auto& p : _node_pileups) {
        delete p.second;
    }
    _node_pileups.clear();

    for (auto& p : _edge_pileups) {
        delete p.second;
    }
    _edge_pileups.clear();
}

void Pileups::to_json(ostream& out) {
    out << "{\"node_pileups\": [";
    for (NodePileupHash::iterator i = _node_pileups.begin(); i != _node_pileups.end();) {
        out << pb2json(*i->second);
        ++i;
        if (i != _node_pileups.end()) {
            out << ",";
        }
    }
    out << "]," << endl << "\"edge_pileups\": [";
    for (EdgePileupHash::iterator i = _edge_pileups.begin(); i != _edge_pileups.end();) {
        out << pb2json(*i->second);
        ++i;
        if (i != _edge_pileups.end()) {
            out << ",";
        }
    }
    out << "]}" << endl;
}

void Pileups::load(istream& in) {
    function<void(Pileup&)> lambda = [this](Pileup& pileup) {
        extend(pileup);
    };
    stream::for_each(in, lambda);
}

void Pileups::write(ostream& out, uint64_t chunk_size) {

    int64_t count = max(_node_pileups.size(), _edge_pileups.size()) / chunk_size;
    if (max(_node_pileups.size(), _edge_pileups.size()) % chunk_size != 0) {
        ++count;
    }

    NodePileupHash::iterator node_it = _node_pileups.begin();
    EdgePileupHash::iterator edge_it = _edge_pileups.begin();
    Pileup pileup;

    // note: this won't work at all in parallel but presumably write
    // is single threaded...
    function<Pileup&(uint64_t)> lambda = [&](uint64_t i) -> Pileup& {
        pileup.clear_node_pileups();
        pileup.clear_edge_pileups();
        for (int j = 0; j < chunk_size && node_it != _node_pileups.end(); ++j, ++node_it) {
            NodePileup* np = pileup.add_node_pileups();
            *np = *node_it->second;
        }
        // unlike for Graph, we don't bother to try to group edges with nodes they attach
        for (int j = 0; j < chunk_size && edge_it != _edge_pileups.end(); ++j, ++edge_it) {
            EdgePileup* ep = pileup.add_edge_pileups();
            *ep = *edge_it->second;
        }

        return pileup;
    };

    stream::write(out, count, lambda);
}

void Pileups::for_each_node_pileup(const function<void(NodePileup&)>& lambda) {
    for (auto& p : _node_pileups) {
        lambda(*p.second);
    }
}

void Pileups::for_each_edge_pileup(const function<void(EdgePileup&)>& lambda) {
    for (auto& p : _edge_pileups) {
        lambda(*p.second);
    }
}

EdgePileup* Pileups::get_edge_pileup(pair<NodeSide, NodeSide> sides) {
    if (sides.second < sides.first) {
        swap(sides.first, sides.second);
    }
    auto p = _edge_pileups.find(sides);
    return p != _edge_pileups.end() ? p->second : NULL;
}
            
// get a pileup.  if it's null, create a new one and insert it.
EdgePileup* Pileups::get_create_edge_pileup(pair<NodeSide, NodeSide> sides) {
    if (sides.second < sides.first) {
        swap(sides.first, sides.second);
    }
    EdgePileup* p = get_edge_pileup(sides);
    if (p == NULL) {
        p = new EdgePileup();
        p->mutable_edge()->set_from(sides.first.node);
        p->mutable_edge()->set_from_start(!sides.first.is_end);
        p->mutable_edge()->set_to(sides.second.node);
        p->mutable_edge()->set_to_end(sides.second.is_end);
        _edge_pileups[sides] = p;
    }
    return p;
}


void Pileups::extend(Pileup& pileup) {
    for (int i = 0; i < pileup.node_pileups_size(); ++i) {
        insert_node_pileup(new NodePileup(pileup.node_pileups(i)));
    }
    for (int i = 0; i < pileup.edge_pileups_size(); ++i) {
        insert_edge_pileup(new EdgePileup(pileup.edge_pileups(i)));
    }
}

bool Pileups::insert_node_pileup(NodePileup* pileup) {
    NodePileup* existing = get_node_pileup(pileup->node_id());
    if (existing != NULL) {
        merge_node_pileups(*existing, *pileup);
        delete pileup;
    } else {
        _node_pileups[pileup->node_id()] = pileup;
    }
    return existing == NULL;
}

bool Pileups::insert_edge_pileup(EdgePileup* pileup) {
    EdgePileup* existing = get_edge_pileup(NodeSide::pair_from_edge(*pileup->mutable_edge()));
    if (existing != NULL) {
        merge_edge_pileups(*existing, *pileup);
        delete pileup;
    } else {
        _edge_pileups[NodeSide::pair_from_edge(*pileup->mutable_edge())] = pileup;
    }
    return existing == NULL;
}

void Pileups::compute_from_alignment(Alignment& alignment) {
    const Path& path = alignment.path();
    int64_t read_offset = 0;
    vector<int> mismatch_counts;
    count_mismatches(*_graph, path, mismatch_counts);
    // element i = location of rank i in the mapping array
    vector<int> ranks(path.mapping_size() + 1, -1);
    // keep track of read offset of mapping array element i
    vector<int64_t> in_read_offsets(path.mapping_size());
    vector<int64_t> out_read_offsets(path.mapping_size());
    // keep track of last mapping, offset of match, and open deletion for
    // calling deletion endpoints (which are beside, but not on the base offsets they get written to)
    pair<const Mapping*, int64_t> last_match(NULL, -1);
    pair<const Mapping*, int64_t> last_del(NULL, -1);
    pair<const Mapping*, int64_t> open_del(NULL, -1);
    for (int i = 0; i < path.mapping_size(); ++i) {
        const Mapping& mapping = path.mapping(i);
        int rank = mapping.rank() <= 0 ? i + 1 : mapping.rank(); 
        if (_graph->has_node(mapping.position().node_id())) {
            const Node* node = _graph->get_node(mapping.position().node_id());
            NodePileup* pileup = get_create_node_pileup(node);
            int64_t node_offset = mapping.position().offset();
            // utilize forward-relative node offset (old way), which
            // is not consistent with current protobuf.  conversion here.  
            if (mapping.position().is_reverse()) {
                node_offset = node->sequence().length() - 1 - node_offset;
            }
            in_read_offsets[i] = read_offset;
            for (int j = 0; j < mapping.edit_size(); ++j) {
                const Edit& edit = mapping.edit(j);
                // process all pileups in edit.
                // update the offsets as we go
                compute_from_edit(*pileup, node_offset, read_offset, *node,
                                  alignment, mapping, edit, mismatch_counts, last_match, last_del, open_del);
            }
            out_read_offsets[i] = read_offset - 1;

            if (rank <= 0 || rank >= ranks.size() || ranks[rank] != -1) {
            cerr << "Error determining rank of mapping " << i << " in path " << path.name() << ": "
                 << pb2json(mapping) << endl;
            }
            else {
                ranks[rank] = i;
            }
        } else {
            // node not in graph. that's okay, we do nothing but update the read_offset to
            // not trigger assert at end of this function
            for (int j = 0; j < mapping.edit_size(); ++j) {
                read_offset += mapping.edit(j).to_length();
            }
            ranks[rank] = -1;
        }
    }
    // loop again over all the edges crossed by the mapping alignment, using
    // the offsets and ranking information we got in the first pass
    for (int i = 2; i < ranks.size(); ++i) {
        int rank1_idx = ranks[i-1];
        int rank2_idx = ranks[i];
        if ((rank1_idx > 0 || rank2_idx > 0) && (rank1_idx >= 0 && rank2_idx >= 0)) {
            auto& m1 = path.mapping(rank1_idx);
            auto& m2 = path.mapping(rank2_idx);
            auto s1 = NodeSide(m1.position().node_id(), (m1.position().is_reverse() ? false : true));
            auto s2 = NodeSide(m2.position().node_id(), (m2.position().is_reverse() ? true : false));
            // no quality gives a free pass from quality filter
            char edge_qual = 127;
            if (!alignment.quality().empty()) {
                char from_qual = alignment.quality()[out_read_offsets[rank1_idx]];
                char to_qual = alignment.quality()[in_read_offsets[rank2_idx]];
                edge_qual = min(from_qual, to_qual);
            }
            if (edge_qual >= _min_quality) {
                EdgePileup* edge_pileup = get_create_edge_pileup(pair<NodeSide, NodeSide>(s1, s2));
                if (edge_pileup->num_reads() < _max_depth) {
                  edge_pileup->set_num_reads(edge_pileup->num_reads() + 1);
                  if (!alignment.quality().empty()) {
                    *edge_pileup->mutable_qualities() += edge_qual;
                  }
                }
            }
        }
    }
    
    assert(alignment.sequence().empty() ||
           alignment.path().mapping_size() == 0 ||
           read_offset == alignment.sequence().length());

}

void Pileups::compute_from_edit(NodePileup& pileup, int64_t& node_offset,
                                int64_t& read_offset,
                                const Node& node, const Alignment& alignment,
                                const Mapping& mapping, const Edit& edit,
                                const vector<int>& mismatch_counts,
                                pair<const Mapping*, int64_t>& last_match,
                                pair<const Mapping*, int64_t>& last_del,
                                pair<const Mapping*, int64_t>& open_del) {
    string seq = edit.sequence();
    // is the mapping reversed wrt read sequence? use for iterating
    bool map_reverse = mapping.position().is_reverse();
    
    // ***** MATCH *****
    if (edit.from_length() == edit.to_length()) {
        assert (edit.from_length() > 0);
        make_match(seq, edit.from_length(), map_reverse);
        assert(seq.length() == edit.from_length());            
        int64_t delta = map_reverse ? -1 : 1;
        for (int64_t i = 0; i < edit.from_length(); ++i) {
            if (pass_filter(alignment, read_offset, mismatch_counts)) {
                BasePileup* base_pileup = get_create_base_pileup(pileup, node_offset);
                if (base_pileup->num_bases() < _max_depth) {
                    // reference_base if empty
                    if (base_pileup->num_bases() == 0) {
                        base_pileup->set_ref_base(node.sequence()[node_offset]);
                    } else {
                        assert(base_pileup->ref_base() == node.sequence()[node_offset]);
                    }
                    // add base to bases field (converting to ,. if match)
                    char base = seq[i];
                    *base_pileup->mutable_bases() += base;
                    // add quality if there
                    if (!alignment.quality().empty()) {
                        *base_pileup->mutable_qualities() += alignment.quality()[read_offset];
                    }
                    // pileup size increases by 1
                    base_pileup->set_num_bases(base_pileup->num_bases() + 1);
                }
                // close off any open deletion
                if (open_del.first != NULL) {
                    string del_seq;
                    make_delete(del_seq, map_reverse, last_match, mapping, node_offset);
                    int64_t dp_node_id;
                    int64_t dp_node_offset;
                    // store in canonical position
                    if (make_pair(make_pair(last_del.first->position().node_id(), last_del.second),
                                  last_del.first->position().is_reverse()) <
                        make_pair(make_pair(open_del.first->position().node_id(), open_del.second),
                                  open_del.first->position().is_reverse())) {
                        dp_node_id = last_del.first->position().node_id();
                        dp_node_offset = last_del.second;
                    } else {
                        dp_node_id = open_del.first->position().node_id();
                        dp_node_offset = open_del.second;
                    }
                    Node* dp_node = _graph->get_node(dp_node_id);
                    NodePileup* dp_node_pileup = get_create_node_pileup(dp_node);
                    BasePileup* dp_base_pileup = get_create_base_pileup(*dp_node_pileup, dp_node_offset);
                    if (dp_base_pileup->num_bases() < _max_depth) {
                        // reference_base if empty
                        if (dp_base_pileup->num_bases() == 0) {
                            dp_base_pileup->set_ref_base(dp_node->sequence()[dp_node_offset]);
                        } else {
                            assert(dp_base_pileup->ref_base() == dp_node->sequence()[dp_node_offset]);
                        }
                        *dp_base_pileup->mutable_bases() += del_seq;
                        if (!alignment.quality().empty()) {
                            // we only use quality of one endpoint here.  should average
                            *dp_base_pileup->mutable_qualities() += alignment.quality()[read_offset];
                        }
                        dp_base_pileup->set_num_bases(dp_base_pileup->num_bases() + 1);
                    }
                    open_del = make_pair((Mapping*)NULL, -1);
                    last_del = make_pair((Mapping*)NULL, -1);
                }
                
                last_match = make_pair(&mapping, node_offset);
            }
            // move right along read, and left/right depending on strand on reference
            node_offset += delta;
            ++read_offset;
        }
    }
    // ***** INSERT *****
    else if (edit.from_length() < edit.to_length()) {
        if (pass_filter(alignment, read_offset, mismatch_counts)) {
            make_insert(seq, map_reverse);
            assert(edit.from_length() == 0);
            // we define insert (like sam) as insertion between current and next
            // position (on forward node coordinates). this means an insertion before
            // offset 0 is invalid! 
            int64_t insert_offset =  map_reverse ? node_offset : node_offset - 1;
            if (insert_offset >= 0) {        
                BasePileup* base_pileup = get_create_base_pileup(pileup, insert_offset);
                if (base_pileup->num_bases() < _max_depth) {
                    // reference_base if empty
                    if (base_pileup->num_bases() == 0) {
                        base_pileup->set_ref_base(node.sequence()[insert_offset]);
                    } else {
                        assert(base_pileup->ref_base() == node.sequence()[insert_offset]);
                    }
                    // add insertion string to bases field
                    base_pileup->mutable_bases()->append(seq);
                    if (!alignment.quality().empty()) {
                        *base_pileup->mutable_qualities() += alignment.quality()[read_offset];
                    }
                    // pileup size increases by 1
                    base_pileup->set_num_bases(base_pileup->num_bases() + 1);
                }
            }
            else {
                // need to check with aligner to make sure this doesn't happen, ie
                // inserts would hang off the end of previous node instead of start
                // of this node
                /*
                  stringstream ss;
                  ss << "Warning: pileup does not support insertions before 0th base in node."
                  << " Offending edit: " << pb2json(edit) << endl;
                  #pragma omp critical(cerr)
                  cerr << ss.str();
                */
            }
        }
        // move right along read (and stay put on reference)
        read_offset += edit.to_length();
    }
    // ***** DELETE *****
    else {
        if (pass_filter(alignment, read_offset, mismatch_counts)) {
            assert(edit.to_length() == 0);
            assert(edit.sequence().empty());

            // deltion will get written in the "Match" section
            // note: deletions will only get written if there's a match on either side
            // so deletions at beginning/end of read ignored in pileup
            if (open_del.first == NULL && last_match.first != NULL) {
                open_del = make_pair(&mapping, node_offset);
            }
            // open_del : first base deleted by deleltion
            // last_del : most recent base deleted by deletion
            // last_match : most recent base in a match
            // (most recent is in order we are scanning here)

            // a deletion will be an edge between two matches.
            // but in the pileup, it will be stored in either open_del or last_del
            // (which ever has lower coordinate).  
        }
        int64_t delta = map_reverse ? -edit.from_length() : edit.from_length();
        
        // stay put on read, move left/right depending on strand on reference
        node_offset += delta;

        last_del = make_pair(&mapping, map_reverse ? node_offset + 1 : node_offset - 1);
    }
}

void Pileups::count_mismatches(VG& graph, const Path& path,
                               vector<int>& mismatches,
                               bool skipIndels)
{
    mismatches.clear();
    int64_t read_offset = 0;
    for (int i = 0; i < path.mapping_size(); ++i) {
        const Mapping& mapping = path.mapping(i);
        if (graph.has_node(mapping.position().node_id())) {
            const Node* node = graph.get_node(mapping.position().node_id());
            int64_t node_offset = mapping.position().offset();
            // utilize forward-relative node offset (old way), which
            // is not consistent with current protobuf.  conversion here.  
            if (mapping.position().is_reverse()) {
                node_offset = node->sequence().length() - 1 - node_offset;
            }

            for (int j = 0; j < mapping.edit_size(); ++j) {
                const Edit& edit = mapping.edit(j);
                // process all pileups in edit.
                // update the offsets as we go
                string seq = edit.sequence();
                bool is_reverse = mapping.position().is_reverse();
                if (is_reverse) {
                    seq = reverse_complement(seq);
                }
    
                // ***** MATCH *****
                if (edit.from_length() == edit.to_length()) {
                    int64_t delta = is_reverse ? -1 : 1;
                    for (int64_t i = 0; i < edit.from_length(); ++i) {
                        if (!edit.sequence().empty() &&
                            !base_equal(seq[i], node->sequence()[node_offset], false)) {
                            mismatches.push_back(1);
                        }
                        else {
                            mismatches.push_back(0);
                        }
                        // move right along read, and left/right depending on strand on reference
                        node_offset += delta;
                        ++read_offset;
                    }
                }
                // ***** INSERT *****
                else if (edit.from_length() < edit.to_length()) {
                    if (skipIndels == false) {
                        mismatches.push_back(1);
                        for (int x = 1; x < edit.to_length(); ++x) {
                            mismatches.push_back(0);
                        }
                    }
                    // move right along read (and stay put on reference)
                    read_offset += edit.to_length();
                }
                // ***** DELETE *****
                else {
                    if (skipIndels == false) {
                        // since we're working in read coordinates, we count
                        // a single mismatch right before the delete.
                        if (mismatches.size() > 0) {
                            mismatches[mismatches.size() - 1] = 1;
                        }
                    }
                    int64_t delta = is_reverse ? -edit.from_length() : edit.from_length();
                    // stay put on read, move left/right depending on strand on reference
                    node_offset += delta;
                }
            }
        } else {
            // node not in graph: count 0 mismatches for each absent position
            for (int j = 0; j < mapping.edit_size(); ++j) {
                read_offset += mapping.edit(j).to_length();
                for (int k = 0; k < mapping.edit(j).to_length(); ++k) {
                    mismatches.push_back(0);
                }
            }
        }
    }
    assert(skipIndels || read_offset == mismatches.size());
    // too lazy to do full count inline.  sum up here
    int count = 0;
    for (int i = 0; i < mismatches.size(); ++i) {
        count += mismatches[i];
        mismatches[i] = count;
    }
}

bool Pileups::pass_filter(const Alignment& alignment, int64_t read_offset,
                          const vector<int>& mismatches) const
{
    bool passes = true;
    if (!alignment.quality().empty()) {
        passes = alignment.quality()[read_offset] >= _min_quality;
    }
    if (_window_size > 0 && passes) {
        // counts in left window
        int64_t left_point = max((int64_t)0, read_offset - _window_size / 2 - 1);
        int64_t right_point = max((int64_t)0, read_offset - 1);
        int64_t count = mismatches[right_point] - mismatches[left_point];
        // coutns in right window
        left_point = read_offset;
        right_point = min(read_offset + _window_size / 2, (int64_t)mismatches.size() - 1);
        count += mismatches[right_point] - mismatches[left_point];
        passes = passes && count <= _max_mismatches;
    }
    return passes;
}

Pileups& Pileups::merge(Pileups& other) {
    for (auto& p : other._node_pileups) {
        insert_node_pileup(p.second);
    }
    other._node_pileups.clear();
    for (auto& p : other._edge_pileups) {
        insert_edge_pileup(p.second);
    }
    other._edge_pileups.clear();
    return *this;
}

BasePileup& Pileups::merge_base_pileups(BasePileup& p1, BasePileup& p2) {
    assert(p1.num_bases() == 0 || p2.num_bases() == 0 ||
           p1.ref_base() == p2.ref_base());
    if (p1.num_bases() == 0) {
        p1.set_ref_base(p2.ref_base());
    }
    int merge_size = min(p2.num_bases(), _max_depth - p1.num_bases());
    p1.set_num_bases(p1.num_bases() + merge_size);
    if (merge_size == p2.num_bases()) {    
        p1.mutable_bases()->append(p2.bases());
        p1.mutable_qualities()->append(p2.qualities());
    } else if (merge_size > 0) {
        vector<pair<int64_t, int64_t> > offsets;
        parse_base_offsets(p2, offsets);
        int merge_length = offsets[merge_size].first;
        p1.mutable_bases()->append(p2.bases().substr(0, merge_length));
        if (!p2.qualities().empty()) {
            p1.mutable_qualities()->append(p2.qualities().substr(0, merge_size));
        }
    }
    p2.set_num_bases(0);
    p2.clear_bases();
    p2.clear_qualities();
    return p1;
}

NodePileup& Pileups::merge_node_pileups(NodePileup& p1, NodePileup& p2) {
    assert(p1.node_id() == p2.node_id());
    for (int i = 0; i < p2.base_pileup_size(); ++i) {
        BasePileup* bp1 = get_create_base_pileup(p1, i);
        BasePileup* bp2 = get_base_pileup(p2, i);
        merge_base_pileups(*bp1, *bp2);
    }
    p2.clear_base_pileup();
    return p1;
}

EdgePileup& Pileups::merge_edge_pileups(EdgePileup& p1, EdgePileup& p2) {
    assert(p1.edge().from() == p2.edge().from());
    assert(p1.edge().to() == p2.edge().to());
    assert(p1.edge().from_start() == p2.edge().from_start());
    assert(p1.edge().to_end() == p2.edge().to_end());
    int merge_size = min(p2.num_reads(), _max_depth - p1.num_reads());
    p1.set_num_reads(p1.num_reads() + merge_size);
    if (merge_size == p2.num_reads()) {
        p1.mutable_qualities()->append(p2.qualities());
    } else if (!p2.qualities().empty()) {
        p1.mutable_qualities()->append(p2.qualities().substr(0, merge_size));
    }
    p2.set_num_reads(0);
    p2.clear_qualities();
    return p1;
}

void Pileups::parse_base_offsets(const BasePileup& bp,
                                 vector<pair<int64_t, int64_t> >& offsets) {
    offsets.clear();
    
    const string& quals = bp.qualities();
    const string& bases = bp.bases();
    char ref_base = ::toupper(bp.ref_base());
    // we can use i to index the quality for the ith row of pileup, but
    // need base_offset to get position of appropriate token in bases string
    int64_t base_offset = 0;
    for (int i = 0; i < bp.num_bases(); ++i) {
        // insert
        if (bases[base_offset] == '+') {
            offsets.push_back(make_pair(base_offset, i < quals.length() ? i : -1));
            int64_t lf = base_offset + 1;
            int64_t rf = lf;
            while (rf < bases.length() && bases[rf] >= '0' && bases[rf] <= '9') {
                ++rf;
            }
            stringstream ss(bases.substr(lf, rf - lf + 1));
            int64_t indel_len;
            ss >> indel_len;
            // ex: +5aaaaa.  rf = lf = 1. indel_len = 5 -> increment 2+0+5=7
            base_offset += 1 + rf - lf + indel_len;
        // delete
        } else if (bases[base_offset] == '-') {
            offsets.push_back(make_pair(base_offset, i < quals.length() ? i : -1));
            int64_t lf = base_offset + 1;
            // eat up six semicolons
            for (int64_t sc_count = 0; sc_count < 6; ++lf) {
                if (bases[lf] == ';') {
                    ++sc_count;
                }
            }
            // and last number
            for (; bases[lf] >= '0' && bases[lf] <= '9'; ++lf);
            base_offset = lf;
        }
        // match / snp
        else {
            offsets.push_back(make_pair(base_offset, i < quals.length() ? i : -1));
            ++base_offset;
        }
    }
    assert(base_offset == bases.length());
}

// transform case of every character in string
void Pileups::casify(string& seq, bool is_reverse) {
    if (is_reverse) {
        transform(seq.begin(), seq.end(), seq.begin(), ::tolower);
    } else {
        transform(seq.begin(), seq.end(), seq.begin(), ::toupper);
    }
}

// make the sam pileup style token
void Pileups::make_match(string& seq, int64_t from_length, bool is_reverse) {
    if (seq.length() == 0) {
        seq = string(from_length, is_reverse ? ',' : '.');
    } else {
        casify(seq, is_reverse);
    }
}

void Pileups::make_insert(string& seq, bool is_reverse) {
    casify(seq, is_reverse);
    stringstream ss;
    ss << "+" << seq.length() << seq; 
    seq = ss.str();
}

void Pileups::make_delete(string& seq, bool is_reverse, const pair<const Mapping*, int64_t>& last_match,
                          const Mapping& mapping, int64_t node_offset){
    int64_t from_id = last_match.first->position().node_id();
    int64_t from_offset = last_match.second;
    bool from_start = last_match.first->position().is_reverse();
    int64_t to_id = mapping.position().node_id();
    int64_t to_offset = node_offset;
    bool to_end = mapping.position().is_reverse();

    // canonical order
    if (make_pair(make_pair(from_id, from_offset), from_start) >
        make_pair(make_pair(to_id, to_offset), to_end)) {
        swap(from_id, to_id);
        swap(from_offset, to_offset);
        swap(from_start, to_end);
        from_start = !from_start;
        to_end = !to_end;
    }
    
    make_delete(seq, is_reverse, from_id, from_offset, from_start, to_id, to_offset, to_end);
}

void Pileups::make_delete(string& seq, bool is_reverse,
                          int64_t from_id, int64_t from_offset, bool from_start,
                          int64_t to_id, int64_t to_offset, bool to_end) {
    // format : -is_reverse;from_id;from_offset;from_start;to_id;do_offset;to_end
    stringstream ss;
    ss << "-" << is_reverse << ";" << from_id << ";" << from_offset << ";" << from_start << ";"
       << to_id << ";" << to_offset << ";" << to_end;
    seq = ss.str();
}
        
void Pileups::parse_insert(const string& tok, int64_t& len, string& seq, bool& is_reverse) {
    assert(tok[0] == '+');
    int64_t i = 1;
    for (; tok[i] >= '0' && tok[i] <= '9'; ++i);
    stringstream ss;
    ss << tok.substr(1, i - 1);
    ss >> len;
    seq = tok.substr(i, tok.length() - i);
    is_reverse = ::islower(seq[0]);
}

void Pileups::parse_delete(const string& tok, bool& is_reverse,
                           int64_t& from_id, int64_t& from_offset, bool& from_start,
                           int64_t& to_id, int64_t& to_offset, bool& to_end) {
    assert(tok[0] == '-');
    vector<string> toks;
    regex sc_re(";");
    std::copy(sregex_token_iterator(tok.begin(), tok.end(), sc_re, -1),
              sregex_token_iterator(), back_inserter(toks));

    assert(toks.size() == 7);
    is_reverse = std::stoi(toks[0]) != 0;

    from_id = std::stoi(toks[1]);
    from_offset = std::stoi(toks[2]);
    from_start = std::stoi(toks[3]) != 0;
    
    to_id = std::stoi(toks[4]);
    to_offset = std::stoi(toks[5]);
    to_end = std::stoi(toks[6]) != 0;
}
    
bool Pileups::base_equal(char c1, char c2, bool is_reverse) {
    char t1 = ::toupper(c1);
    char t2 = ::toupper(c2);
    return is_reverse ? t1 == reverse_complement(t2) : t1 == t2;
}

char Pileups::extract_match(const BasePileup& bp, int64_t offset) {
    char v = bp.bases()[offset];
    assert(v != '+' && v != '-');
    if (v == ',' || v == '.') {
        return ::toupper(bp.ref_base());
    } else if (::islower(v)) {
        return reverse_complement(::toupper(v));
    }
    return v;
}

// get arbitrary value from offset on forward strand
string Pileups::extract(const BasePileup& bp, int64_t offset) {
    const string& bases = bp.bases();
    if (bases[offset] != '+' && bases[offset] != '-') {
        return string(1, extract_match(bp, offset));
    }
    else if (bases[offset] == '+') {
        string len_str;
        for (int64_t i = offset + 1; bases[i] >= '0' && bases[i] <= '9'; ++i) {
            len_str += bases[i];
        }
        int64_t len = atoi(len_str.c_str());
        // forward strand, return as is
        if (::isupper(bases[offset + 1 + len_str.length()])) {
            return bases.substr(offset, 1 + len_str.length() + len);
        }
        // reverse strand, flip the dna bases part and return upper case
        else {
            string dna = bases.substr(offset + 1 + len_str.length(), len);
            casify(dna, false);
            return string(1, bases[offset]) + len_str + reverse_complement(dna);
        }
    }
    else {
        assert(bases[offset] == '-');
        // todo : consolidate deletion parsing code better than this
        int64_t sc = 0;
        int64_t i = offset;
        for (; sc < 6; ++i) {
            if (bases[i] == ';') {
                ++sc;
            }
        }
        return bases.substr(offset, i - offset + 1);
    }
}

}
