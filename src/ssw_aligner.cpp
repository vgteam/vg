#include "ssw_aligner.hpp"

namespace vg {

void SSWAligner::PrintAlignment(const StripedSmithWaterman::Alignment& alignment){
  cout << "===== SSW result =====" << endl;
  cout << "Best Smith-Waterman score:\t" << alignment.sw_score << endl
       << "Next-best Smith-Waterman score:\t" << alignment.sw_score_next_best << endl
       << "Reference start:\t" << alignment.ref_begin << endl
       << "Reference end:\t" << alignment.ref_end << endl
       << "Query start:\t" << alignment.query_begin << endl
       << "Query end:\t" << alignment.query_end << endl
       << "Next-best reference end:\t" << alignment.ref_end_next_best << endl
       << "Number of mismatches:\t" << alignment.mismatches << endl
       << "Cigar: " << alignment.cigar_string << endl;
  cout << "======================" << endl;
}

Alignment SSWAligner::align(const string& query, const string& ref) {
    StripedSmithWaterman::Aligner aligner(match,
                                          mismatch,
                                          gap_open,
                                          gap_extension);
    StripedSmithWaterman::Filter filter;
    StripedSmithWaterman::Alignment alignment;
    
    // We need to send out own mask length, recommended to be half the sequence length and at least 15.
    int32_t mask_len = min(max((size_t) 15, query.size()), (size_t) numeric_limits<int32_t>::max());
    
    assert(ref.size() <= numeric_limits<int>::max());
    
    aligner.Align(query.c_str(), ref.c_str(), (int) ref.size(), filter, &alignment, mask_len);
    return ssw_to_vg(alignment, query, ref);
}

Alignment SSWAligner::ssw_to_vg(const StripedSmithWaterman::Alignment& ssw_aln,
                                const string& query, const string& ref) {

    int from_pos = ssw_aln.ref_begin;
    int to_pos = 0;

    auto& to_seq = query;
    auto& from_seq = ref;

    Alignment vg_aln;
    Path* path = vg_aln.mutable_path();
    
    Mapping* mapping = path->add_mapping();
    mapping->mutable_position()->set_offset(from_pos);

    for (auto& elem : vcflib::splitCigar(ssw_aln.cigar_string)) {
        int32_t length = elem.first;
        string type(1, elem.second);
        Edit* edit;
        //cerr << e->length << e->type << endl;
        switch (type[0]) {
	case '=':
        case 'M':
        case 'X':
        case 'N': {
            // do the sequences match?
            // emit a stream of "SNPs" and matches
            int h = from_pos;
            int last_start = from_pos;
            int k = to_pos;
            for ( ; h < from_pos + length; ++h, ++k) {
                //cerr << h << ":" << k << " " << from_seq[h] << " " << to_seq[k] << endl;
                if (from_seq[h] != to_seq[k]) {
                    // emit the last "match" region
                    if (h-last_start > 0) {
                        edit = mapping->add_edit();
                        edit->set_from_length(h-last_start);
                        edit->set_to_length(h-last_start);
                    }
                    // set up the SNP
                    edit = mapping->add_edit();
                    edit->set_from_length(1);
                    edit->set_to_length(1);
                    edit->set_sequence(to_seq.substr(k,1));
                    last_start = h+1;
                }
            }
            // handles the match at the end or the case of no SNP
            if (h-last_start > 0) {
                edit = mapping->add_edit();
                edit->set_from_length(h-last_start);
                edit->set_to_length(h-last_start);
            }
            to_pos += length;
            from_pos += length;
        } break;
        case 'D':
            edit = mapping->add_edit();
            edit->set_from_length(length);
            edit->set_to_length(0);
            from_pos += length;
            break;
        case 'I':
            edit = mapping->add_edit();
            edit->set_from_length(0);
            edit->set_to_length(length);
            edit->set_sequence(to_seq.substr(to_pos, length));
            to_pos += length;
            break;
        case 'S':
            // note that soft clips and insertions are semantically equivalent
            // and can only be differentiated by their position in the read
            // with soft clips coming at the start or end
            edit = mapping->add_edit();
            edit->set_from_length(0);
            edit->set_to_length(length);
            edit->set_sequence(to_seq.substr(to_pos, length));
            to_pos += length;
            break;
        default:
            cerr << "error [GSSWAligner::gssw_mapping_to_alignment] "
                 << "unsupported cigar op type " << type << endl;
            exit(1);
            break;

        }

    }

    // set identity
    vg_aln.set_identity(identity(vg_aln.path()));

    // set score
    vg_aln.set_score(ssw_aln.sw_score);
    
    return vg_aln;

}

}
