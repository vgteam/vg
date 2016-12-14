#include "sampler.hpp"

namespace vg {

pos_t Sampler::position(void) {
    uniform_int_distribution<size_t> xdist(1, xgidx->seq_length);
    size_t offset = xdist(rng);
    id_t id = xgidx->node_at_seq_pos(offset);
    uniform_int_distribution<size_t> flip(0, 1);
    bool rev = forward_only ? false : flip(rng);
    // 1-0 base conversion
    size_t node_offset = offset - xgidx->node_start(id) - 1;
    return make_pos_t(id, rev, node_offset);
}

string Sampler::sequence(size_t length) {
    pos_t pos = position();
    cerr << pos << endl;
    string seq;
    while (seq.size() < length) {
        auto nextc = next_pos_chars(pos);
        if (nextc.empty()) break;
        vector<pos_t> nextp;
        for (auto& n : nextc) nextp.push_back(n.first);
        // pick one at random
        uniform_int_distribution<int> next_dist(0, nextc.size()-1);
        // update our position
        pos = nextp.at(next_dist(rng));
        // append to our sequence
        seq += nextc[pos];
    }
    return seq;
}


vector<Edit> Sampler::mutate_edit(const Edit& edit,
                                  const pos_t& position,
                                  double base_error,
                                  double indel_error,
                                  const string& bases,
                                  uniform_real_distribution<double>& rprob,
                                  uniform_int_distribution<int>& rbase) {

    // we will build up a mapping representing the modified edit
    Mapping new_mapping;
    //*new_mapping.mutable_position() = make_position(position);
    // determine to-length of edit
    size_t to_length = edit.to_length();
    // we will keep track of the current base using this
    pos_t curr_pos = position;
    /// TODO we should punt if we aren't a pure edit
    // as in, we are something with mixed to and from lengths; like a block sub with an indel
    if (edit_is_match(edit) || edit_is_sub(edit)
        || edit_is_insertion(edit)) {
        // distribute mutations across this length
        for (size_t k = 0; k < to_length; ++k) {
            char c = 'N'; // in the case that we are in an insertion
            if (!edit_is_insertion(edit)) {
                c = pos_char(curr_pos);
                ++get_offset(curr_pos);
            }
            if (rprob(rng) <= base_error) {
                // pick another base than what c is
                char n;
                do {
                    n = bases[rbase(rng)];
                } while (n == c);
                // make the edit for the sub
                Edit* e = new_mapping.add_edit();
                string s(1, n);
                e->set_sequence(s);
                e->set_from_length(1);
                e->set_to_length(1);
            // if we've got a indel
            // note that we're using a simple geometric indel dsitribution here
            } else if (rprob(rng) <= indel_error) {
                if (rprob(rng) < 0.5) {
                    char n = bases[rbase(rng)];
                    Edit* e = new_mapping.add_edit();
                    string s(1, c);
                    e->set_sequence(s);
                    e->set_to_length(1);
                } else {
                    Edit* e = new_mapping.add_edit();
                    e->set_from_length(1);
                }
            } else {
                // make the edit for the 1bp match
                Edit* e = new_mapping.add_edit();
                e->set_from_length(1);
                e->set_to_length(1);
            }

        }
    } else if (edit_is_deletion(edit)) {
        // special case: 0 (deletion)
        // maybe we do nothing; as there is no length in the read
    }
    // simplify the mapping
    new_mapping = simplify(new_mapping);
    // copy the new edits
    vector<Edit> new_edits;
    for (size_t i = 0; i < new_mapping.edit_size(); ++i) {
        new_edits.push_back(new_mapping.edit(i));
    }
    // and send them back
    return new_edits;
}

Alignment Sampler::mutate(const Alignment& aln,
                          double base_error,
                          double indel_error) {

    if (base_error == 0 && indel_error == 0) return aln;

    string bases = "ATGC";
    uniform_real_distribution<double> rprob(0, 1);
    uniform_int_distribution<int> rbase(0, 3);

    Alignment mutaln;
    for (size_t i = 0; i < aln.path().mapping_size(); ++i) {
        auto& orig_mapping = aln.path().mapping(i);
        Mapping* new_mapping = mutaln.mutable_path()->add_mapping();
        *new_mapping->mutable_position() = orig_mapping.position();
        // for each edit in the mapping
        for (size_t j = 0; j < orig_mapping.edit_size(); ++j) {
            auto& orig_edit = orig_mapping.edit(j);
            auto new_edits = mutate_edit(orig_edit, make_pos_t(orig_mapping.position()),
                                         base_error, indel_error,
                                         bases, rprob, rbase);
            for (auto& edit : new_edits) {
                *new_mapping->add_edit() = edit;
            }
        }
    }
    // re-derive the alignment's sequence
    mutaln = simplify(mutaln);
    mutaln.set_sequence(alignment_seq(mutaln));
    mutaln.set_name(aln.name());
    return mutaln;
}

string Sampler::alignment_seq(const Alignment& aln) {
    // get the graph corresponding to the alignment path
    Graph sub;
    for (int i = 0; i < aln.path().mapping_size(); ++ i) {
        auto& m = aln.path().mapping(i);
        if (m.has_position() && m.position().node_id()) {
            auto id = aln.path().mapping(i).position().node_id();
            xgidx->neighborhood(id, 2, sub);
        }
    }
    VG g; g.extend(sub);
    return g.path_string(aln.path());
}

vector<Alignment> Sampler::alignment_pair(size_t read_length, size_t fragment_length, double fragment_std_dev, double base_error, double indel_error) {
    // simulate forward/reverse pair by first simulating a long read
    normal_distribution<> norm_dist(fragment_length, fragment_std_dev);
    // bound at read length so we always get enough sequence
    int frag_len = max((int)read_length, (int)round(norm_dist(rng)));
    auto fragment = alignment_with_error(frag_len, base_error, indel_error);
    // then taking the ends
    auto fragments = alignment_ends(fragment, read_length, read_length);
    auto& aln1 = fragments.front();
    auto& aln2 = fragments.back();
    { // name the alignments
        string data;
        aln1.SerializeToString(&data);
        aln2.SerializeToString(&data);
        int n;
#pragma omp critical(nonce)
        n = nonce++;
        data += std::to_string(n);
        const string hash = sha1head(data, 16);
        aln1.set_name(hash + "_1");
        aln2.set_name(hash + "_2");
    }
    // set the appropriate flags for pairing
    aln1.mutable_fragment_next()->set_name(aln2.name());
    aln2.mutable_fragment_prev()->set_name(aln1.name());
    // reverse complement the back fragment
    fragments.back() = reverse_complement_alignment(fragments.back(),
                                                  (function<int64_t(int64_t)>) ([&](int64_t id) {
                                                          return (int64_t)node_length(id);
                                                      }));
    return fragments;
}

// generates a perfect alignment from the graph
Alignment Sampler::alignment(size_t length) {
    string seq;
    Alignment aln;
    Path* path = aln.mutable_path();
    pos_t pos = position();
    char c = pos_char(pos);
    // we do something wildly inefficient but conceptually clean
    // for each position in the mapping we add a mapping
    // at the end we will simplify the alignment, merging redundant mappings
    do {
        // add in the char for the current position
        seq += c;
        Mapping* mapping = path->add_mapping();
        *mapping->mutable_position() = make_position(pos);
        Edit* edit = mapping->add_edit();
        edit->set_from_length(1);
        edit->set_to_length(1);
        // decide the next position
        auto nextc = next_pos_chars(pos);
        // no new positions mean we are done; we've reached the end of the graph
        if (nextc.empty()) break;
        // what positions do we go to next?
        vector<pos_t> nextp;
        for (auto& n : nextc) nextp.push_back(n.first);
        // pick one at random
        uniform_int_distribution<int> next_dist(0, nextc.size()-1);
        // update our position
        pos = nextp.at(next_dist(rng));
        // update our char
        c = nextc[pos];
    } while (seq.size() < length);
    // save our sequence in the alignment
    aln.set_sequence(seq);
    aln = simplify(aln);
    { // name the alignment
        string data;
        aln.SerializeToString(&data);
        int n;
#pragma omp critical(nonce)
        n = nonce++;
        data += std::to_string(n);
        const string hash = sha1head(data, 16);
        aln.set_name(hash);
    }
    // and simplify it
    aln.set_identity(identity(aln.path()));
    return aln;
}

Alignment Sampler::alignment_with_error(size_t length,
                                        double base_error,
                                        double indel_error) {
    size_t maxiter = 100;
    Alignment aln;
    size_t iter = 0;
    if (base_error > 0 || indel_error > 0) {
        // sample a longer-than necessary alignment, then trim
        while (iter++ < maxiter) {
            aln = mutate(
                alignment(length + 2 * ((double) length * indel_error)),
                base_error, indel_error);
            if (aln.sequence().size() == length) {
                break;
            } else if (aln.sequence().size() > length) {
                aln = strip_from_end(aln, aln.sequence().size() - length);
                break;
            }
        }
    } else {
        size_t iter = 0;
        while (iter++ < maxiter) {
            aln = alignment(length);
            if (aln.sequence().size() == length) {
                break;
            }
        }
    }
    if (iter == maxiter) {
        cerr << "[vg::Sampler] Warning: could not generate alignment of sufficient length. "
             << "Graph may be too small, or indel rate too high." << endl;
    }
    aln.set_identity(identity(aln.path()));
    return aln;
}

size_t Sampler::node_length(id_t id) {
    return xg_cached_node_length(id, xgidx, node_cache);
}

char Sampler::pos_char(pos_t pos) {
    return xg_cached_pos_char(pos, xgidx, node_cache);
}

map<pos_t, char> Sampler::next_pos_chars(pos_t pos) {
    return xg_cached_next_pos_chars(pos, xgidx, node_cache);
}

}
