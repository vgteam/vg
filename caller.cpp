#include <cstdlib>
#include <stdexcept>
#include "json2pb.h"
#include "caller.hpp"
#include "stream.hpp"

using namespace std;

namespace vg {

void Caller::clear() {
    _alignments.clear();
}

void Caller::to_json(ostream& out) {
    for (auto aln : _alignments) {
        out << pb2json(aln);
    }
}

void Caller::write(ostream& out) {
    function<Alignment&(uint64_t)> lambda = [&](uint64_t i) -> Alignment& {
        return _alignments[i];
    };
    stream::write(out, _alignments.size(), lambda);
}

void Caller::flush_buffer(ostream& out, bool json) {    
#pragma omp critical (out)
    {
        if (json) {
            to_json(out);
        } else {
            write(out);
        }
    }
    _alignments.clear();
}

void Caller::call_node_pileup(const NodePileup& pileup, ostream& out, bool json) {
    for (int i = 0; i < pileup.base_pileup_size(); ++i) {
        call_base_pileup(pileup, i);
    }
    if (_alignments.size() >= _buffer_size) {
        flush_buffer(out, json);
    }
}

void Caller::call_base_pileup(const NodePileup& np, int64_t offset) {
    // todo: replace with something a little smarter that
    // takes into account mapping quality / error probability.
    const BasePileup& bp = np.base_pileup(offset);
    if (bp.num_bases() >= _min_depth) {
        // acgtn
        int hist[5] = {0};
        ++hist[nidx(bp.ref_base())];
        int base_offset = 0;
        for (int i = 0; i < bp.num_bases(); ++i) {
            // match
            if (bp.bases()[base_offset] == ',' ||
                bp.bases()[base_offset] == '.') {
                ++hist[nidx(bp.ref_base())];
                ++base_offset;
            }
            // snp
            else if (nidx(bp.bases()[base_offset]) < 4) {
                ++hist[nidx(bp.bases()[base_offset])];
                ++base_offset;
            }
            // indel
            else if (bp.bases()[base_offset] == '+' ||
                     bp.bases()[base_offset] == '-') {
                // todo
                int lf = base_offset + 1;
                int rf = bp.bases().find_last_of("0123456789", lf);
                stringstream ss(bp.bases().substr(lf, rf - lf + 1));
                int indel_len;
                ss >> indel_len;
                // ex: +5aaaaa.  rf = lf = 1. indel_len = 5 -> increment 2+0+5=7
                base_offset += 2 + rf - lf + indel_len;
            }
        }

        // find top two bases in pileup
        int first = -1;
        int second = -1;
        for (int i = 0; i < 4; ++i) {
            if (first < 0 || hist[i] > hist[first]) {
                first = i;
            } else if (second < 0 || hist[i] > hist[second]) {
                second = i;
            }
        }

        // make up to two snps
        if (idxn(first) != bp.ref_base() &&
            (double)hist[first] / (double)bp.num_bases() >= _min_frac) {
            create_snp(np, offset, idxn(first));
        }
        if (second != first &&
            idxn(second) != bp.ref_base() &&
            (double)hist[second] / (double)bp.num_bases() >= _min_frac) {
            create_snp(np, offset, idxn(second));
        }
    }   
}

void Caller::create_snp(const NodePileup& np, int64_t offset, char base) {
    _alignments.push_back(Alignment());
    Alignment& alignment = _alignments.back();
    string sequence;
    // note on context: since we don't have access to the Node object,
    // need to get sequence out of the pileups...
    // left context (run right-to-left so we can cut at empty base more easily
    for (int64_t i = offset - 1; i >= std::max((int64_t)0, offset - _context); --i) {
        const BasePileup* bp = Pileups::get_base_pileup(np, i);
        if (!bp || bp->num_bases() == 0) {
            break;
        } else {
            sequence.insert(0, 1, (char)bp->ref_base());
        }
    }
    // snp base
    int seq_offset = sequence.length();
    sequence += base;
    // right context
    for (int64_t i = offset + 1;
         i < std::min((int64_t)np.base_pileup_size(), offset + 1 + _context); ++i) {
        const BasePileup* bp = Pileups::get_base_pileup(np, i);
        if (!bp || bp->num_bases() == 0) {
            break;
        } else {
            sequence += (char)bp->ref_base();
        }
    }
    alignment.set_sequence(sequence);
    Path* path = alignment.mutable_path();
    Mapping* mapping = path->add_mapping();
    mapping->mutable_position()->set_node_id(np.node_id());
    mapping->mutable_position()->set_offset(offset - seq_offset);
    mapping->set_is_reverse(false);
    // create match edit from left context
    if (seq_offset > 0) {
        Edit* edit = mapping->add_edit();
        edit->set_to_length(seq_offset);
        edit->set_from_length(seq_offset);
    }
    // create edit from snp
    Edit* edit = mapping->add_edit();
    edit->set_to_length(1);
    edit->set_from_length(1);
    edit->set_sequence(string(1, base));
    // create edit from right context
    if (seq_offset < sequence.length() - 1) {
        Edit* edit = mapping->add_edit();
        edit->set_to_length(sequence.length() - 1 - seq_offset);
        edit->set_from_length(sequence.length() - 1 - seq_offset);
    }
}


}
