#include "sampler.hpp"
#include "path.hpp"

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
    
    // We punt if we aren't a well-defined kind of edit; from and to lengths
    // must be equal, or from length must be 0.
    if (edit_is_match(edit) || edit_is_sub(edit)
        || edit_is_insertion(edit)) {
        
#ifdef debug
        cerr << "Handle edit " << pb2json(edit) << endl;
#endif
        
        // distribute mutations across the to_length
        for (size_t k = 0; k < to_length;) {
            // Until we've consumed all the characters produced by the original edit
        
            // Get the character that the original edit has here
            char c;
            if (edit_is_match(edit)) {
                // It's not stored in the edit, so find it in the reference
                c = pos_char(curr_pos);
            } else {
                // It must be stored in the edit itself
                c = edit.sequence().at(k);
            }
            
#ifdef debug
            cerr << "At to_length position " << k << " of " << edit.to_length() << endl;
#endif
            
            // This is the edit we're going to create, to deside the fate of
            // this character from the old edit.
            Edit* e = nullptr;
            
            if (rprob(rng) <= base_error) {
                // We should do a substitution relative to the old edit.
                
                // pick another base than what c is
                char n;
                do {
                    n = bases[rbase(rng)];
                } while (n == c);
                // make the edit for the sub
                e = new_mapping.add_edit();
                string s(1, n);
                e->set_sequence(s);
                if (!edit_is_insertion(edit)) {
                    // We should stay aligned against whatever we were aligned
                    // against before.
                    e->set_from_length(1);
                }
                e->set_to_length(1);
                
#ifdef debug
                cerr << "Produced relative substitution " << pb2json(*e) << endl;
#endif
            } else if (rprob(rng) <= indel_error) {
                // We have an indel.
                // Note that we're using a simple geometric indel dsitribution here
                if (rprob(rng) < 0.5) {
                    // This should be an insertion relative to the original edit.
                    char n = bases[rbase(rng)];
                    e = new_mapping.add_edit();
                    string s(1, c);
                    e->set_sequence(s);
                    e->set_to_length(1);
                    
#ifdef debug
                    cerr << "Produced relative insertion " << pb2json(*e) << endl;
#endif
                    
                    // None of the graph sequence is consumed. But we are
                    // inserting relative to the old edit, so we want to give
                    // the base we just inserted before another shot. We'll
                    // continue now, before advancing our position in the edit
                    // we're modifying.
                    continue;
                } else {
                    // This should be a deletion relative to the edit we're
                    // modifying.
                    
                    // The old edit base isn't going to come out, but we need to
                    // consume it anyway.
                    k++;
                    
                    if (edit_is_insertion(edit)) {
                        // We just need to consume this base of the old edit and
                        // not emit any new edit.
#ifdef debug
                        cerr << "Skipped base for relative deletion" << endl;
#endif
                        continue;
                    } else {
                        // We have to delete the base that was substituted or
                        // matched against.
                        e = new_mapping.add_edit();
                        e->set_from_length(1);
#ifdef debug
                        cerr << "Produced relative deletion " << pb2json(*e) << endl;
#endif
                    }
                }
            } else {
                // make the edit for the 1bp match relative to the old edit
                // (which may actually be an insertion edit or a substitution
                // edit relative to the graph)
                e = new_mapping.add_edit();
                if (!edit_is_match(edit)) {
                    // We're not a match, so we need the sequence set.
                    string s(1, c);
                    e->set_sequence(s);
                }
                if (!edit_is_insertion(edit)) {
                    // We should stay aligned against whatever we were aligned
                    // against before.
                    e->set_from_length(1);
                }
                e->set_to_length(1);
                
#ifdef debug
                cerr << "Produced relative match " << pb2json(*e) << endl;
#endif
            }
            
            // Now advance in the old edit by the number of old edit bases used
            // in this new edit.
            k += e->to_length();
            // And in the graph by the number of graph bases consumed
            get_offset(curr_pos) += e->from_length();

        }
    } else if (edit_is_deletion(edit)) {
        // special case: 0 (deletion)
        // Just copy over the deletion edit
        *(new_mapping.add_edit()) = edit;
    }
    
#ifdef debug
    cerr << "Before merging adjacent edits: " << pb2json(new_mapping) << endl;
#endif

    // Merge adjacent edits, but don't get rid of leading or trailing deletions
    // (as with simplify), because we want a path that reflects the real
    // simulated history and because we don't output any modifications to the
    // position out of this function.
    new_mapping = merge_adjacent_edits(new_mapping);
    
#ifdef debug
    cerr << "Replacing " << pb2json(edit) << " with " << pb2json(new_mapping) << endl;
#endif
    
    assert(mapping_from_length(new_mapping) == edit.from_length());
    
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
    
    // Don't simplify the alignment, because it's nice to see the deletions as
    // originally generated. Also, simplified alignments won't validate.
    
    // re-derive the alignment's sequence.
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
            xgidx->get_id_range(id, id, sub);
        }
    }
    xgidx->expand_context(sub, 2, false);
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
    // Simplify the alignment to merge redundant mappings. There are no deletions to get removed.
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
    // And set its identity
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
    
    // Check the alignment to make sure we didn't mess it up
    assert(is_valid(aln));
    
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

bool Sampler::is_valid(const Alignment& aln) {
    for (auto i = 0; i + 1 < aln.path().mapping_size(); i++) {
        // For each mapping except the very last (which might not use its whole
        // node)
        auto& mapping = aln.path().mapping(i);
        
        // What's the number of bases it consumes?
        auto observed_from = mapping_from_length(mapping);
        
        // How many bases are accounted for?
        auto accounted_bases = observed_from + mapping.position().offset();
        
        // How many bases need to be accounted for?
        auto expected_bases = xgidx->node_length(mapping.position().node_id());
        
        if (accounted_bases != expected_bases) {
            cerr << "[vg::Sampler] Warning: alignment mapping " << i << " accounts for "
                << accounted_bases << " bases of graph sequence, but needs to account for "
                << expected_bases << endl;
            cerr << pb2json(aln) << endl;
            return false;
        }
    }
    
    // For now, we just say an alignment is valid if it accounts for all the
    // bases on its source nodes.
    return true;
}

}
