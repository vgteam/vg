#include <cstdlib>
#include <stdexcept>
#include "json2pb.h"
#include "pileup.hpp"
#include "stream.hpp"

using namespace std;

namespace vg {

void Pileups::clear() {
    for (auto& p : _pileups) {
        delete p.second;
    }
    _pileups.clear();
}

void Pileups::to_json(ostream& out) {
    function<void(NodePileup&)> lambda = [this, &out](NodePileup& p) {
        out << pb2json(p) <<endl;
    };
    for_each(lambda);
}

void Pileups::load(istream& in) {
    function<void(NodePileup&)> lambda = [this](NodePileup& pileup) {
        insert(new NodePileup(pileup));
    };
    stream::for_each(in, lambda);
}

void Pileups::write(ostream& out, uint64_t buffer_size) {
    // maybe there's a more efficient way of getting
    // these into the stream from the hash table?  for now
    // we just buffer in a vector
    vector<NodePileup*> buf;
    function<NodePileup&(uint64_t)> lambda = [&](uint64_t i) -> NodePileup& {
        return *buf[i];
    };
    for (auto& p : _pileups) {
        buf.push_back(p.second);
        if (buf.size() >= buffer_size) {
            stream::write(out, buf.size(), lambda);
            buf.clear();
        }  
    }
    if (!buf.empty()) {
        stream::write(out, buf.size(), lambda);
    }
}

void Pileups::for_each(function<void(NodePileup&)>& lambda) {
    for (auto& p : _pileups) {
        lambda(*p.second);
    }
}

bool Pileups::insert(NodePileup* pileup) {
    NodePileup* existing = get(pileup->node_id());
    if (existing != NULL) {
        merge_node_pileups(*existing, *pileup);
        delete pileup;
    } else {
        _pileups[pileup->node_id()] = pileup;
    }
    return existing == NULL;
}

void Pileups::compute_from_alignment(VG& graph, Alignment& alignment) {
    if (alignment.is_reverse()) {
        throw runtime_error("pileup doesn\'t currently support reverse flag "
                            "in Alignment.  Reverse mappings must be specified "
                            "using the is_reverse falg in Mapping for now.");
    }
    const Path& path = alignment.path();
    int64_t read_offset = 0;
    for (int i = 0; i < path.mapping_size(); ++i) {
        const Mapping& mapping = path.mapping(i);
        if (graph.has_node(mapping.position().node_id())) {
            const Node* node = graph.get_node(mapping.position().node_id());
            NodePileup* pileup = get_create(node->id());
            int64_t node_offset = mapping.position().offset();
            for (int j = 0; j < mapping.edit_size(); ++j) {
                const Edit& edit = mapping.edit(j);
                // process all pileups in edit.
                // update the offsets as we go
                compute_from_edit(*pileup, node_offset, read_offset, *node,
                                  alignment, mapping, edit);
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
                                const Mapping& mapping, const Edit& edit) {
    string seq = edit.sequence();
    
    // ***** MATCH *****
    if (edit.from_length() == edit.to_length()) {
        assert (edit.from_length() > 0);
        make_match(seq, edit.from_length(), mapping.is_reverse());
        assert(seq.length() == edit.from_length());            
        int64_t delta = mapping.is_reverse() ? -1 : 1;
        for (int64_t i = 0; i < edit.from_length(); ++i) {
            BasePileup* base_pileup = get_create_base_pileup(pileup, node_offset);
            if (base_pileup->num_bases() == 0) {
                base_pileup->set_ref_base(node.sequence()[node_offset]);
            } else {
                assert(base_pileup->ref_base() == node.sequence()[node_offset]);
            }
            *base_pileup->mutable_bases() += seq[i];
            if (!alignment.quality().empty()) {
                *base_pileup->mutable_qualities() += alignment.quality()[read_offset];
            }
            base_pileup->set_num_bases(base_pileup->num_bases() + 1);
            node_offset += delta;
            ++read_offset;
        }
    }
    // ***** INSERT *****
    else if (edit.from_length() < edit.to_length()) {
        make_insert(seq, mapping.is_reverse());
        assert(edit.from_length() == 0);
        // we define insert (like sam) as insertion between current and next
        // position (on forward node coordinates). this means an insertion before
        // offset 0 is invalid! 
        int64_t insert_offset =  mapping.is_reverse() ? node_offset : node_offset - 1;
        if (insert_offset < 0) {
            stringstream ss;
            ss << "Error: pileup does not support insertions before 0th base in node."
               << " Offending edit: " << pb2json(edit) << "\nin alignment "
               << pb2json(alignment) << endl;
            throw runtime_error(ss.str());
        }
        BasePileup* base_pileup = get_create_base_pileup(pileup, insert_offset);
        if (base_pileup->num_bases() == 0) {
            base_pileup->set_ref_base(node.sequence()[insert_offset]);
        } else {
            assert(base_pileup->ref_base() == node.sequence()[insert_offset]);
        }
        base_pileup->mutable_bases()->append(seq);
        if (!alignment.quality().empty()) {
            *base_pileup->mutable_qualities() += alignment.quality()[read_offset];
        }
        base_pileup->set_num_bases(base_pileup->num_bases() + 1);
        read_offset += edit.to_length();
    }
    // ***** DELETE *****
    else {
        int64_t del_start = !mapping.is_reverse() ? node_offset :
            node_offset - edit.from_length();
        seq = node.sequence().substr(del_start, edit.from_length());
        make_delete(seq, mapping.is_reverse());        
        assert(edit.to_length() == 0);
        BasePileup* base_pileup = get_create_base_pileup(pileup, node_offset);
        if (base_pileup->num_bases() == 0) {
            base_pileup->set_ref_base(node.sequence()[node_offset]);
        } else {
            assert(base_pileup->ref_base() == node.sequence()[node_offset]);
        }
        base_pileup->mutable_bases()->append(seq);
        if (!alignment.quality().empty()) {
            *base_pileup->mutable_qualities() += alignment.quality()[read_offset];
        }
        base_pileup->set_num_bases(base_pileup->num_bases() + 1);
        int64_t delta = mapping.is_reverse() ? -edit.from_length() : edit.from_length();
        node_offset += delta;
        read_offset += edit.to_length();
    }
}

Pileups& Pileups::merge(Pileups& other) {
    for (auto& p : other._pileups) {
        insert(p.second);
    }
    other._pileups.clear();
    return *this;
}

BasePileup& Pileups::merge_base_pileups(BasePileup& p1, BasePileup& p2) {
    assert(p1.num_bases() == 0 || p2.num_bases() == 0 ||
           p1.ref_base() == p2.ref_base());
    if (p1.num_bases() == 0) {
        p1.set_ref_base(p2.ref_base());
    }
    p1.mutable_bases()->append(p2.bases());
    p1.mutable_qualities()->append(p2.qualities());
    p1.set_num_bases(p1.num_bases() + p2.num_bases());
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


}
