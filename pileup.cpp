#include <cstdlib>
#include <algorithm>
#include <stdexcept>
#include "json2pb.h"
#include "pileup.hpp"
#include "stream.hpp"

using namespace std;

namespace vg {

void Pileups::clear() {
    for (auto& p : _pileups) {
        delete p;
    }
    _pileups.clear();
}

void Pileups::to_json(ostream& out) {
    function<void(Pileup&)> lambda = [this, &out](Pileup& p) {
        out << pb2json(p) <<endl;
    };
    for_each(lambda);
}

void Pileups::load(istream& in) {
    function<void(Pileup&)> lambda = [this](Pileup& pileup) {
        insert(new Pileup(pileup));
    };
    stream::for_each(in, lambda);
}

void Pileups::write(ostream& out, uint64_t buffer_size) {
    // maybe there's a more efficient way of getting
    // these into the stream from the hash table?  for now
    // we just buffer in a vector
    vector<Pileup*> buf;
    function<Pileup&(uint64_t)> lambda = [&](uint64_t i) -> Pileup& {
        return *buf[i];
    };
    for (auto& p : _pileups) {
        buf.push_back(p);
        if (buf.size() >= buffer_size) {
            stream::write(out, buf.size(), lambda);
            buf.clear();
        }  
    }
    if (!buf.empty()) {
        stream::write(out, buf.size(), lambda);
    }
}

void Pileups::for_each(function<void(Pileup&)>& lambda) {
    for (auto& p : _pileups) {
        lambda(*p);
    }
}

bool Pileups::insert(Pileup* pileup) {
    Pileup* existing = get(*pileup);
    if (existing != NULL) {
        merge_pileup(*existing, *pileup);
        delete pileup;
    } else {
        assert(pileup->position().offset() >= 0);
        _pileups.insert(pileup);
    }
    return existing == NULL;
}

void Pileups::compute_from_alignment(VG& graph, Alignment& alignment) {
    Pileup query;
    const Path& path = alignment.path();
    int64_t read_offset = 0;
    for (int i = 0; i < path.mapping_size(); ++i) {
        const Mapping& mapping = path.mapping(i);
        if (graph.has_node(mapping.position().node_id())) {
            const Node* node = graph.get_node(mapping.position().node_id());
            query.mutable_position()->set_node_id(mapping.position().node_id());
            query.mutable_position()->set_offset(mapping.position().offset());
            for (int j = 0; j < mapping.edit_size(); ++j) {
                const Edit& edit = mapping.edit(j);
                // process all pileups in edit.
                // update the query position accordingly (to 1st coord of next edit)
                compute_from_edit(query, read_offset, *node, alignment, mapping, edit);
            }
        }
    }
}

void Pileups::compute_from_edit(Pileup& query, int64_t& read_offset,
                                const Node& node, const Alignment& alignment,
                                const Mapping& mapping, const Edit& edit) {
    
    string seq = edit.sequence();
    if (mapping.is_reverse()) {
        transform(seq.begin(), seq.end(), seq.begin(), ::tolower);
    } else {
        transform(seq.begin(), seq.end(), seq.begin(), ::toupper);
    }
    
    // ***** MATCH *****
    if (edit.from_length() == edit.to_length()) {
        assert (edit.from_length() > 0);
        if (seq.length() == 0) {
            seq = string(edit.from_length(), mapping.is_reverse() ? ',' : '.');
        }
        assert(seq.length() == edit.from_length());            
        int64_t delta = mapping.is_reverse() ? -1 : 1;
        for (int64_t i = 0; i < edit.from_length(); ++i) {
            query.set_ref_base(node.sequence()[query.position().offset()]);
            Pileup* pileup = get_create(query);
            *pileup->mutable_bases() += seq[i];
            if (!alignment.quality().empty()) {
                *pileup->mutable_qualities() += alignment.quality()[read_offset++];
            }
            pileup->set_num_bases(pileup->num_bases() + 1);
            query.mutable_position()->set_offset(query.position().offset() + delta);
        }
    }
    // ***** INSERT *****
    else if (edit.from_length() < edit.to_length()) {
        assert(edit.from_length() == 0);
        // clamp to endpoint, as query's position always set to the next
        // from_position (which may not exist). todo need to find if there's
        // some kind of convention that needs to be applied here!
        if (mapping.is_reverse() && query.position().offset() == -1) {
            query.mutable_position()->set_offset(0);
        } else if (!mapping.is_reverse() && query.position().offset() ==
                   node.sequence().length()) {
            query.mutable_position()->set_offset(node.sequence().length() - 1);
        }        
        query.set_ref_base(node.sequence()[query.position().offset()]);
        Pileup* pileup = get_create(query);
        pileup->mutable_bases()->append("+" + seq);
        if (!alignment.quality().empty()) {
            *pileup->mutable_qualities() += alignment.quality()[read_offset];
            read_offset += seq.length();
        }
        pileup->set_num_bases(pileup->num_bases() + 1);
    }
    // ***** DELETE *****
    else {
        assert(edit.to_length() == 0);
        // todo : edge cases, clipping ?
        query.set_ref_base(node.sequence()[query.position().offset()]);
        Pileup* pileup = get_create(query);
        pileup->mutable_bases()->append("-" + seq);
        if (!alignment.quality().empty()) {
            *pileup->mutable_qualities() += alignment.quality()[read_offset];
            read_offset += seq.length();
        }
        pileup->set_num_bases(pileup->num_bases() + 1);
        int64_t delta = mapping.is_reverse() ? -edit.from_length() : edit.from_length();
        query.mutable_position()->set_offset(query.position().offset() + delta);
    }
}

Pileups& Pileups::merge(Pileups& other) {
    for (auto& p : other._pileups) {
        insert(p);
    }
    other._pileups.clear();
    return *this;
}

Pileup& merge_pileup(Pileup& p1, const Pileup& p2) {
    assert(p1.position().node_id() == p2.position().node_id());
    assert(p1.position().offset() == p2.position().offset());
    assert(p1.ref_base() == p2.ref_base());
    p1.mutable_bases()->append(p2.bases());
    p1.mutable_qualities()->append(p2.qualities());
    p1.set_num_bases(p1.num_bases() + p2.num_bases());
    return p1;
}


}
