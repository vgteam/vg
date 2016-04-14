#include "sampler.hpp"

namespace vg {

pos_t Sampler::position(void) {
    uniform_int_distribution<size_t> xdist(1, xgidx->seq_length);
    size_t offset = xdist(rng);
    //cerr << "offset " << offset << endl;
    //cerr << "id at " << xgidx->node_at_seq_pos(offset) << endl;
    id_t id = xgidx->node_at_seq_pos(offset);
    size_t node_offset = offset - xgidx->node_start(id);
    return make_pos_t(id, false, node_offset);
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
                ++get_offset(curr_pos);
                c = pos_char(curr_pos);
            }
            if (rprob(rng) <= base_error) {
                // pick another base than what c is
                char n;
                do {
                    n = bases[rbase(rng)];
                } while (n == c);
                // make the edit for the sub
                Edit* e = new_mapping.add_edit();
                e->set_sequence(std::to_string(n));
                e->set_from_length(1);
                e->set_to_length(1);
            } else {
                // make the edit for the 1bp match
                Edit* e = new_mapping.add_edit();
                e->set_from_length(1);
                e->set_to_length(1);
            }
            // if we've got a indel
            // note that we're using a simple geometric indel dsitribution here
            if (rprob(rng) <= indel_error) {
                if (rprob(rng) < 0.5) {
                    char n = bases[rbase(rng)];
                    Edit* e = new_mapping.add_edit();
                    e->set_sequence(std::to_string(n));
                    e->set_to_length(1);
                } else {
                    Edit* e = new_mapping.add_edit();
                    e->set_from_length(1);
                }
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
    return mutaln;
}


Alignment Sampler::alignment(size_t length) {
    Alignment aln;
    return aln;
}

char Sampler::pos_char(pos_t pos) {
    return xgidx->pos_char(id(pos), is_rev(pos), offset(pos));
}

map<pos_t, char> Sampler::next_pos_chars(pos_t pos) {
    map<pos_t, char> nexts;
    Node node = xgidx->node(id(pos));
    // if we are still in the node, return the next position and character
    if (offset(pos) < node.sequence().size()-1) {
        ++get_offset(pos);
        nexts[pos] = pos_char(pos);
    } else {
        // look at the next positions we could reach
        if (!is_rev(pos)) {
            // we are on the forward strand, the next things from this node come off the end
            for (auto& edge : xgidx->edges_on_end(id(pos))) {
                if (edge.from() == id(pos)) {
                    pos_t p = make_pos_t(edge.to(), edge.to_end(), 0);
                    nexts[p] = pos_char(p);
                } else if (edge.from_start() && edge.to_end() && edge.to() == id(pos)) {
                    // doubly inverted, should be normalized to forward but we handle here for safety
                    pos_t p = make_pos_t(edge.from(), false, 0);
                    nexts[p] = pos_char(p);
                }
            }
        } else {
            // we are on the reverse strand, the next things from this node come off the start
            for (auto& edge : xgidx->edges_on_start(id(pos))) {
                if (edge.to() == id(pos)) {
                    pos_t p = make_pos_t(edge.from(), !edge.from_start(), 0);
                    nexts[p] = pos_char(p);
                } else if (edge.from_start() && edge.to_end() && edge.from() == id(pos)) {
                    // doubly inverted, should be normalized to forward but we handle here for safety
                    pos_t p = make_pos_t(edge.to(), true, 0);
                    nexts[p] = pos_char(p);
                }
            }
        }
    }
    return nexts;
}

}
