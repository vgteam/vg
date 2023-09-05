#include "sampler.hpp"

#include "path.hpp"
#include "utility.hpp"
#include "position.hpp"
#include "alignment.hpp"
#include "algorithms/path_string.hpp"
#include "algorithms/subgraph.hpp"
#include "algorithms/alignment_path_offsets.hpp"
#include "algorithms/next_pos_chars.hpp"

#include <sstream>

//#define debug_ngs_sim

using namespace vg::io;

namespace vg {

void AbstractReadSampler::annotate_with_path_positions(Alignment& aln) {
    // We need to annotate the alignment with the right kind of path positions
    // that we are configured for.
    if (multi_position_annotations) {
        algorithms::annotate_with_node_path_positions(graph, aln, 0, annotation_path_filter.get());
    } else {
        algorithms::annotate_with_initial_path_positions(graph, aln, 0, annotation_path_filter.get());
    }
}

void Sampler::set_source_paths(const vector<string>& source_paths,
                               const vector<double>& source_path_ploidies,
                               const vector<pair<string, double>>& transcript_expressions,
                               const vector<tuple<string, string, size_t>>& haplotype_transcripts) {
    if (!source_paths.empty() && !transcript_expressions.empty()) {
        cerr << "error:[Sampler] cannot sample from list of paths and from list of transcripts simultaneously" << endl;
        exit(1);
    }
    else if (!haplotype_transcripts.empty() && transcript_expressions.empty()) {
        cerr << "error:[Sampler] cannot sample from haplotype transcripts without an expression profile" << endl;
        exit(1);
    }
    if (!source_path_ploidies.empty() && source_path_ploidies.size() != source_paths.size()) {
        cerr << "error:[Sampler] cannot sample from list of paths with the wrong number of ploidy weights ("
            << source_path_ploidies.size() << " vs. " << source_paths.size() << ")" << endl;
        exit(1);
    }
    else if (!transcript_expressions.empty()) {
        this->source_paths.clear();
        vector<double> expression_values;
        if (haplotype_transcripts.empty()) {
            for (const pair<string, double>& transcript_expression : transcript_expressions) {
                this->source_paths.push_back(transcript_expression.first);
                size_t tx_len = graph.get_path_length(graph.get_path_handle(transcript_expression.first));
                expression_values.push_back(transcript_expression.second * tx_len);
            }
        }
        else {
            unordered_map<string, vector<size_t>> haplotypes_of_transcript;
            for (size_t i = 0; i < haplotype_transcripts.size(); ++i) {
                haplotypes_of_transcript[get<1>(haplotype_transcripts[i])].push_back(i);
            }
            for (const pair<string, double>& transcript_expression : transcript_expressions) {
                size_t total_haplotypes = 0;
                for (size_t i : haplotypes_of_transcript[transcript_expression.first]) {
                    total_haplotypes += get<2>(haplotype_transcripts[i]);
                }
                for (size_t i : haplotypes_of_transcript[transcript_expression.first]) {
                    double haplotype_expression = (transcript_expression.second * get<2>(haplotype_transcripts[i])) / total_haplotypes;
                    size_t hp_tx_len = graph.get_path_length(graph.get_path_handle(get<0>(haplotype_transcripts[i])));
                    expression_values.push_back(haplotype_expression * hp_tx_len);
                    this->source_paths.push_back(get<0>(haplotype_transcripts[i]));
                }
            }
        }
        path_sampler = vg::discrete_distribution<>(expression_values.begin(), expression_values.end());
    }
    else if (!source_paths.empty()) {
        this->source_paths = source_paths;
        vector<double> path_weights;
        path_weights.reserve(source_paths.size());
        for (size_t i = 0; i < source_paths.size(); i++) {
            // For each source path
            auto& source_path = source_paths[i];
            // Grab an applicable ploidy weight, or assume 1
            double ploidy = i >= source_path_ploidies.size() ? 1.0 : source_path_ploidies[i];
            
            // Add each path, weighted by ploidy and length, to the distribution for sampling paths
            path_weights.push_back(ploidy * graph.get_path_length(graph.get_path_handle(source_path)));
        }
        path_sampler = vg::discrete_distribution<>(path_weights.begin(), path_weights.end());
    }
    else {
        path_sampler = vg::discrete_distribution<>();
    }
}
    

/// We have a helper function to convert path positions and orientations to
/// pos_t values.
pos_t position_at(PathPositionHandleGraph* graph_ptr, const string& path_name, const size_t& path_offset, bool is_reverse) {
    path_handle_t path_handle = graph_ptr->get_path_handle(path_name);
    step_handle_t step = graph_ptr->get_step_at_position(path_handle, path_offset);
    handle_t handle = graph_ptr->get_handle_of_step(step);
    
    // Work out where in that mapping we should be.
    size_t node_offset = path_offset - graph_ptr->get_position_of_step(step);

    if (is_reverse) {
        // Flip the node offset around to be from the end and not the start
        node_offset = graph_ptr->get_length(handle) - node_offset - 1;
    }

    // Make a pos_t for where we are, on the appropriate strand
    pos_t pos = make_pos_t(graph_ptr->get_id(handle), graph_ptr->get_is_reverse(handle) != is_reverse, node_offset);
    
    return pos;
}

pos_t Sampler::position(void) {
    // We sample from the entire graph sequence, 1-based.
    vg::uniform_int_distribution<size_t> xdist(1, total_seq_length);
    size_t offset = xdist(rng);
    id_t id = dynamic_cast<VectorizableHandleGraph*>(&graph)->node_at_vector_offset(offset);
    vg::uniform_int_distribution<size_t> flip(0, 1);
    bool rev = forward_only ? false : flip(rng);
    // 1-0 base conversion
    size_t node_offset = offset - dynamic_cast<VectorizableHandleGraph*>(&graph)->node_vector_offset(id) - 1;
    // Ignore flipping the node offset because we're uniform over both strands
    return make_pos_t(id, rev, node_offset);
}

string Sampler::sequence(size_t length) {
    pos_t pos = position();

    string seq;
    while (seq.size() < length) {
        auto nextc = next_pos_chars(pos);
        if (nextc.empty()) break;
        vector<pos_t> nextp;
        for (auto& n : nextc) nextp.push_back(n.first);
        // pick one at random
        vg::uniform_int_distribution<int> next_dist(0, nextc.size()-1);
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
                                  vg::uniform_real_distribution<double>& rprob,
                                  vg::uniform_int_distribution<int>& rbase) {

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
                    string s(1, n);
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
    vg::uniform_real_distribution<double> rprob(0, 1);
    vg::uniform_int_distribution<int> rbase(0, 3);

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
    mutaln.clear_refpos();
    annotate_with_path_positions(mutaln);
    return mutaln;
}

string Sampler::alignment_seq(const Alignment& aln) {
    return algorithms::path_string(graph, aln.path());
}

vector<Alignment> Sampler::alignment_pair(size_t read_length, size_t fragment_length, double fragment_std_dev, double base_error, double indel_error) {
    // simulate forward/reverse pair by first simulating a long read
    vg::normal_distribution<> norm_dist(fragment_length, fragment_std_dev);
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
    // And annotate with true positions
    annotate_with_path_positions(aln1);
    annotate_with_path_positions(aln2);
    return fragments;
}

// generates a perfect alignment from the graph or the selected source path if one exists
Alignment Sampler::alignment(size_t length) {
    if (source_paths.empty()) {
        return alignment_to_graph(length);
    } else {
        return alignment_to_path(source_paths[path_sampler(rng)], length);
    }
}

// generates a perfect alignment from the graph
Alignment Sampler::alignment_to_path(const string& source_path, size_t length) {

    // Pick a starting point along the path and an orientation
    path_handle_t path_handle = graph.get_path_handle(source_path);
    uint64_t path_length = graph.get_path_length(path_handle);
    vg::uniform_int_distribution<size_t> xdist(0, path_length - 1);
    size_t path_offset = xdist(rng);
    vg::uniform_int_distribution<size_t> flip(0, 1);
    bool rev = forward_only ? false : flip(rng);
    
    // We will fill in this string
    string seq;
    // And this Alignment
    Alignment aln;
    Path* path = aln.mutable_path();
    
    while (seq.size() < length) {

        // Make a pos_t for where we are, on the appropriate strand
        pos_t pos = position_at(&graph, source_path, path_offset, rev);

        // Add that character to the sequence
        seq.push_back(pos_char(pos));
        
        // Add a perfect match edit for 1 base
        Mapping* mapping = path->add_mapping();
        *mapping->mutable_position() = make_position(pos);
        Edit* edit = mapping->add_edit();
        edit->set_from_length(1);
        edit->set_to_length(1);
        
        // Advance along the path in the appropriate direction
        if (rev) {
            if (path_offset == 0) {
                // Out of path!
                break;
            }
            path_offset--;
        } else {
            if (path_offset == path_length - 1) {
                // Out of path!
                break;
            }
            path_offset++;
        }
    }
    
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
    aln.clear_refpos();
    annotate_with_path_positions(aln);
    return aln;
}

// generates a perfect alignment from the graph
Alignment Sampler::alignment_to_graph(size_t length) {
    string seq;
    Alignment aln;
    Path* path = aln.mutable_path();
    pos_t pos = position();
    char c = pos_char(pos);
    // we do something wildly inefficient but conceptually clean
    // for each position in the mapping we add a mapping
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
        vg::uniform_int_distribution<int> next_dist(0, nextc.size()-1);
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
    annotate_with_path_positions(aln);
    return aln;
}

Alignment Sampler::alignment_with_error(size_t length,
                                        double base_error,
                                        double indel_error) {
    Alignment aln;
    size_t iter = 0;
    if (base_error > 0 || indel_error > 0) {
        // sample a longer-than necessary alignment, then trim
        while (iter++ < max_tries) {
            aln = mutate(
                alignment(length + 2 * ((double) length * indel_error)),
                base_error, indel_error);
            if (!(no_Ns && aln.sequence().find('N') != string::npos)) {
                if (aln.sequence().size() == length) {
                    break;
                } else if (aln.sequence().size() > length) {
                    aln = strip_from_end(aln, aln.sequence().size() - length);
                    break;
                }
            }
        }
    } else {
        size_t iter = 0;
        while (iter++ < max_tries) {
            aln = alignment(length);
            if (aln.sequence().size() == length
                && !(no_Ns && aln.sequence().find('N') != string::npos)) {
                break;
            }
        }
    }
    if (iter == max_tries) {
        cerr << "[vg::Sampler] Warning: could not generate alignment of sufficient length in "
             << max_tries << " tries. Graph may be too small, or indel rate too high." << endl;
    }
    aln.set_identity(identity(aln.path()));
    
    // Check the alignment to make sure we didn't mess it up
    assert(is_valid(aln));
    annotate_with_path_positions(aln);
    return aln;
}

size_t Sampler::node_length(id_t id) {
    return graph.get_length(graph.get_handle(id));
}

char Sampler::pos_char(pos_t pos) {
    return graph.get_base(graph.get_handle(id(pos), is_rev(pos)), offset(pos));
}

map<pos_t, char> Sampler::next_pos_chars(pos_t pos) {
    return algorithms::next_pos_chars(graph, pos);
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
        auto expected_bases = graph.get_length(graph.get_handle(mapping.position().node_id()));
        
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


const string NGSSimulator::alphabet = "ACGT";

NGSSimulator::NGSSimulator(PathPositionHandleGraph& graph,
                           const string& ngs_fastq_file,
                           const string& ngs_paired_fastq_file,
                           bool interleaved_fastq,
                           const vector<string>& source_paths_input,
                           const vector<double>& source_path_ploidies,
                           const vector<pair<string, double>>& transcript_expressions,
                           const vector<tuple<string, string, size_t>>& haplotype_transcripts,
                           double substition_polymorphism_rate,
                           double indel_polymorphism_rate,
                           double indel_error_proportion,
                           double fragment_length_mean,
                           double fragment_length_stdev,
                           double error_multiplier,
                           bool retry_on_Ns,
                           bool sample_unsheared_paths,
                           uint64_t manual_seed) :
      AbstractReadSampler(graph)
    , sub_poly_rate(substition_polymorphism_rate)
    , indel_poly_rate(indel_polymorphism_rate)
    , indel_error_prop(indel_error_proportion)
    , fragment_mean(fragment_length_mean)
    , fragment_sd(fragment_length_stdev)
    , retry_on_Ns(retry_on_Ns)
    , strand_sampler(0, 1)
    , background_sampler(0, alphabet.size() - 1)
    , mut_sampler(0, alphabet.size() - 2)
    , prob_sampler(0.0, 1.0)
    , seed(manual_seed)
    , source_paths(source_paths_input)
    , joint_initial_distr(manual_seed ? 1760681024122689423ull * manual_seed + 1107607255714504485ull : random_device()())
    , sample_unsheared_paths(sample_unsheared_paths)
{
    if (!ngs_paired_fastq_file.empty() && interleaved_fastq) {
        cerr << "error:[NGSSimulator] cannot indicate interleaved FASTQ and paired FASTQs simultaneously" << endl;
        exit(1);
    }
    
    if (!source_paths.empty() && !transcript_expressions.empty()) {
        cerr << "error:[NGSSimulator] cannot simultaneously limit sampling to paths and match an expression profile" << endl;
        exit(1);
    }
    
    if (!source_path_ploidies.empty() && source_path_ploidies.size() != source_paths.size()) {
        cerr << "error:[NGSSimulator] cannot sample from list of paths with the wrong number of ploidy weights ("
            << source_path_ploidies.size() << " vs. " << source_paths.size() << ")" << endl;
        exit(1);
    }
    
    if (!haplotype_transcripts.empty() && transcript_expressions.empty()) {
        cerr << "error:[NGSSimulator] cannot sample from haplotype transcripts without an expression profile" << endl;
        exit(1);
    }
    
    if (substition_polymorphism_rate < 0.0 || substition_polymorphism_rate > 1.0
        || indel_polymorphism_rate < 0.0 || indel_polymorphism_rate > 1.0
        || indel_error_proportion < 0.0 || indel_error_proportion > 1.0) {
        cerr << "error:[NGSSimulator] All proportions must be between 0.0 and 1.0" << endl;
        exit(1);
    }
    
    if (substition_polymorphism_rate + indel_polymorphism_rate > 1.0) {
        cerr << "error:[NGSSimulator] Indel polymorphism rate and substitution polymorphism rate cannot sum to greater than 1.0" << endl;
        exit(1);
    }
    
    if (fragment_length_mean <= 0.0) {
        cerr << "error:[NGSSimulator] Mean fragment length must be positive" << endl;
        exit(1);
    }
    
    if (fragment_length_stdev < 0.0) {
        cerr << "error:[NGSSimulator] Fragment length standard deviation must be positive" << endl;
        exit(1);
    }
    
    
    
    if (source_paths_input.empty() && transcript_expressions.empty()) {
        // we are sampling from all positions
        graph.for_each_handle([&](const handle_t& handle) {
            total_seq_length += graph.get_length(handle);
        });
        start_pos_samplers.emplace_back(1, total_seq_length);
    }
    else if (!source_paths_input.empty()) {
        // we are sampling from a given set of source paths
        // TODO: Deduplicate with Sampler's code that does almost exactly this.
        vector<double> path_weights;
        path_weights.reserve(source_paths.size());
        for (size_t i = 0; i < source_paths.size(); i++) {
            // For each source path
            auto& source_path = source_paths[i];
            
            size_t length = graph.get_path_length(graph.get_path_handle(source_path));
            
            // Always use accurate length for sampling start pos, even with sample_unsheared_paths
            start_pos_samplers.emplace_back(0, length - 1);
            
            if (length == 0) {
                path_weights.push_back(0.0);
            }
            else if (sample_unsheared_paths) {
                // sample uniformly between paths
                path_weights.push_back(1.0);
            }
            else {
                // Sample paths proportional to effective length and ploidy
                double eff_path_len;
                if (fragment_mean != numeric_limits<double>::max()) {
                    double trunc_mean = vg::truncated_normal_distribution<>(fragment_mean, fragment_sd, 1.0, length).mean();
                    eff_path_len = length - trunc_mean;
                }
                else {
                    eff_path_len = length;
                }
                
                // Grab an applicable ploidy weight, or assume 1 if not set or if using sample_unsheared_paths
                double ploidy = i >= source_path_ploidies.size() ? 1.0 : source_path_ploidies[i];
                
                // Add each path, weighted by ploidy and length, to the distribution for sampling paths
                path_weights.push_back(ploidy * eff_path_len);
            }
        }
        path_sampler = vg::discrete_distribution<>(path_weights.begin(), path_weights.end());
    }
    else {
        // we are sampling according to an expression profile
        vector<double> expression_values;
        
        if (haplotype_transcripts.empty()) {
            // no transcript name file provided, path names should match transcript names in the
            // expression file
            for (const pair<string, double>& transcript_expression : transcript_expressions) {
                size_t tx_len = graph.get_path_length(graph.get_path_handle(transcript_expression.first));
                if (tx_len == 0) {
                    continue;
                }
                source_paths.push_back(transcript_expression.first);
                start_pos_samplers.emplace_back(0, tx_len - 1);
                if (sample_unsheared_paths) {
                    expression_values.push_back(transcript_expression.second);
                }
                else {
                    double eff_tx_len;
                    if (fragment_mean != numeric_limits<double>::max()) {
                        double trunc_mean = vg::truncated_normal_distribution<>(fragment_mean, fragment_sd, 1.0, tx_len).mean();
                        eff_tx_len = tx_len - trunc_mean;
                    }
                    else {
                        eff_tx_len = tx_len;
                    }
                    expression_values.push_back(transcript_expression.second * eff_tx_len);
                }
            }
        }
        else {
            // map the transcript names to the haplotype transcript names
            unordered_map<string, vector<size_t>> haplotypes_of_transcript;
            for (size_t i = 0; i < haplotype_transcripts.size(); ++i) {
                haplotypes_of_transcript[get<1>(haplotype_transcripts[i])].push_back(i);
            }
            for (const pair<string, double>& transcript_expression : transcript_expressions) {
                // split the expression up among the haplotype transcripts according to their count
                size_t total_haplotypes = 0;
                for (size_t i : haplotypes_of_transcript[transcript_expression.first]) {
                    total_haplotypes += get<2>(haplotype_transcripts[i]);
                }
                for (size_t i : haplotypes_of_transcript[transcript_expression.first]) {
                    size_t hp_tx_len = graph.get_path_length(graph.get_path_handle(get<0>(haplotype_transcripts[i])));
                    if (hp_tx_len == 0) {
                        continue;
                    }
                    source_paths.push_back(get<0>(haplotype_transcripts[i]));
                    start_pos_samplers.emplace_back(0, hp_tx_len - 1);
                    double haplotype_expression = (transcript_expression.second * get<2>(haplotype_transcripts[i])) / total_haplotypes;
                    if (sample_unsheared_paths) {
                        expression_values.push_back(haplotype_expression);
                    }
                    else {
                        double eff_hp_tx_len;
                        if (fragment_mean != numeric_limits<double>::max()) {
                            double trunc_mean = vg::truncated_normal_distribution<>(fragment_mean, fragment_sd, 1.0, hp_tx_len).mean();
                            eff_hp_tx_len = hp_tx_len - trunc_mean;
                        }
                        else {
                            eff_hp_tx_len = hp_tx_len;
                        }
                        expression_values.push_back(haplotype_expression * eff_hp_tx_len);
                    }
                }
            }
        }
        
        path_sampler = vg::discrete_distribution<>(expression_values.begin(), expression_values.end());
    }
    
    // memoize phred conversions
    phred_prob.resize(256);
    for (int i = 1; i < phred_prob.size(); i++) {
        phred_prob[i] = error_multiplier * phred_to_prob((uint8_t)i);
    }
    
    for (size_t i = 0; i < alphabet.size(); i++) {
        mutation_alphabets[alphabet[i]] = string();
        for (size_t j = 0; j < alphabet.size(); j++) {
            if (j == i) {
                continue;
            }
            mutation_alphabets[alphabet[i]].push_back(alphabet[j]);
        }
    }
    
    // record read lengths and the empirical distribution of base qualities
    unordered_map<size_t, size_t> length_count;
    if (interleaved_fastq) {
        fastq_paired_interleaved_for_each(ngs_fastq_file, [&](const Alignment& aln_1, const Alignment& aln_2) {
            length_count[aln_1.quality().size()]++;
            length_count[aln_2.quality().size()]++;
            record_read_pair_quality(aln_1, aln_2);
        });
    }
    else if (!ngs_paired_fastq_file.empty()) {
        fastq_paired_two_files_for_each(ngs_fastq_file, ngs_paired_fastq_file, [&](const Alignment& aln_1, const Alignment& aln_2) {
            length_count[aln_1.quality().size()]++;
            length_count[aln_2.quality().size()]++;
            record_read_pair_quality(aln_1, aln_2);
        });
    }
    else {
        fastq_unpaired_for_each(ngs_fastq_file, [&](const Alignment& aln) {
            length_count[aln.quality().size()]++;
            record_read_quality(aln);
        });
    }
    
    // auto-detect the read length
    size_t modal_length = 0;
    size_t modal_length_count = 0;
    size_t total_reads = 0;
    for (const pair<size_t, size_t>& length_record : length_count) {
        if (length_record.second > modal_length_count) {
            modal_length_count = length_record.second;
            modal_length = length_record.first;
        }
        total_reads += length_record.second;
    }
    
    if (((double) modal_length_count) / total_reads < 0.5 && !sample_unsheared_paths) {
        cerr << "warning:[NGSSimulator] Auto-detected read length of " << modal_length << " encompasses less than half of training reads, NGSSimulator is optimized for training data in which most reads are the same length" << endl;
    }
    
    if (modal_length > fragment_length_mean - 2.0 * fragment_length_stdev && !sample_unsheared_paths) {
        cerr << "warning:[NGSSimulator] Auto-detected read length of " << modal_length << " is long compared to mean fragment length " << fragment_length_mean << " and standard deviation " << fragment_length_stdev << ". Sampling may take additional time and the statistical properties of the fragment length distribution may not reflect input parameters." << endl;
    }
    
    // shorten the quality string samplers until they are the modal length (this determines read length later)
    // if we're sampling unsheared paths, take the whole read
    while (transition_distrs_1.size() > modal_length && !sample_unsheared_paths) {
        transition_distrs_1.pop_back();
    }
    while (transition_distrs_2.size() > modal_length && !sample_unsheared_paths) {
        transition_distrs_2.pop_back();
    }
    
    if (transition_distrs_1.size() != transition_distrs_2.size() && transition_distrs_2.size() > 0) {
        cerr << "error:[NGSSimulator] One fragment end in training data has no reads at the modal length, cannot produce joint samples" << endl;
        exit(1);
    }
    
    finalize();
    
    
    uint64_t prng_seed = seed ? seed : random_device()();
    // engine with coding-time random coefficient to produce good seeds for each thread
    // from one seed
    linear_congruential_engine<uint64_t, 1094757125720465369ull, 10230831556735383564ull, 18446744073709551557ull>  seed_perturbor(prng_seed);
    // make a prng for each thread
    for (int i = 0, n = get_thread_count(); i < n; ++i) {
        prngs.emplace_back(seed_perturbor());
    }
    
#ifdef debug_ngs_sim
    cerr << "finished initializing simulator" << endl;
#endif
}

void NGSSimulator::connect_to_position_file(const string& filename) {
    if (source_paths.empty()) {
        cerr << "warning:[NGSSimulator] path position file will not be created because not simulating from paths" << endl;
        return;
    }
    position_file.open(filename);
    if (!position_file) {
        cerr << "error:[NGSSimulator] failed to open position file: " << filename << endl;
        exit(1);
    }
    position_file << "read\tpath\toffset\treverse" << endl;
}

mt19937_64& NGSSimulator::prng() {
    return prngs[omp_get_thread_num()];
}

void NGSSimulator::register_sampled_position(const Alignment& aln, const string& path_name,
                                             size_t offset, bool is_reverse) {
    if (position_file.is_open()) {
        // we're recording positions
        if (is_reverse) {
            // get the position of the end instead of the start
            offset -= path_from_length(aln.path());
        }
        string line = aln.name() + '\t' + path_name + '\t' + to_string(offset) + '\t' + to_string(is_reverse) + '\n';
#pragma omp critical
        position_file << line;
    }
}

Alignment NGSSimulator::sample_read() {
    
    Alignment aln;
    
    aln.set_name(get_read_name());
    
    // sample a quality string based on the trained distribution
    pair<string, vector<bool>> qual_and_masks = sample_read_quality();
    
#ifdef debug_ngs_sim
    cerr << "sampled qualities and N-mask:" << endl;
    cerr << string_quality_short_to_char(qual_and_masks.first) << endl;
    for (bool mask : qual_and_masks.second) {
        cerr << (mask ? "1" : "0");
    }
    cerr << endl;
#endif
    
    assert(qual_and_masks.first.size() == qual_and_masks.second.size());
    
    aln.set_quality(qual_and_masks.first);
    
    // We won't try indefinitely to find a place between Ns
    size_t failures_due_to_n_bases = 0;

    // Sample our path (if dealing with source_paths)
    size_t source_path_idx = sample_path();
    string source_path;
    if (source_path_idx != numeric_limits<size_t>::max()) {
        source_path = source_paths[source_path_idx];
    }
    
    // This is our offset along the source path, if in use
    int64_t sampled_offset;
    // And our direction to go along the source path, if in use
    bool sampled_is_reverse;
    
    // attempt samples until we get one that succeeds without walking
    // off the end of the graph
    while (!aln.has_path()) {
        // Populate the sample positions
        pos_t pos;
        sample_start_pos(source_path_idx, (int64_t) qual_and_masks.first.size(),
                         sampled_offset, sampled_is_reverse, pos);
        // copy the values so that we can change them without forgetting the start location
        int64_t offset = sampled_offset;
        bool is_reverse = sampled_is_reverse;
        // align the first end at this position on the source path or graph
        sample_read_internal(aln, offset, is_reverse, pos, source_path);
        
        // make sure we didn't sample sequence from
        if (retry_on_Ns) {
            if (aln.sequence().find('N') != string::npos) {
                aln.clear_path();
                aln.clear_sequence();
                failures_due_to_n_bases++;
                
                if (failures_due_to_n_bases >= max_tries) {
                    // We have hit Ns too many times in a row and need to bail out or give up.
                    stringstream ss;
                    ss << "Failed to sample a " << std::to_string(aln.quality().size()) << " bp sequence without Ns";
                    if (source_path_idx != numeric_limits<size_t>::max()) {
                        ss << " from path " << source_path;
                    }
                    ss << " for our maximum of " << max_tries << " tries. Is there such a sequence available?";
                    throw std::runtime_error(ss.str());
                }
            }
        }
    }
    
    // mask out any of the sequence that we sampled to be an 'N'
    apply_N_mask(*aln.mutable_sequence(), qual_and_masks.second);
    
    annotate_with_path_positions(aln);
    
    register_sampled_position(aln, source_path, sampled_offset + sampled_is_reverse, sampled_is_reverse);
    return aln;
}

pair<Alignment, Alignment> NGSSimulator::sample_read_pair() {
    pair<Alignment, Alignment> aln_pair;
    
    string name = get_read_name();
    aln_pair.first.set_name(name + "_1");
    aln_pair.second.set_name(name + "_2");
        
    pair<pair<string, vector<bool>>, pair<string, vector<bool>>> qual_and_mask_pair = sample_read_quality_pair();
    
#ifdef debug_ngs_sim
    cerr << "sampled qualities and N-masks:" << endl;
    cerr << string_quality_short_to_char(qual_and_mask_pair.first.first) << endl;
    for (bool mask : qual_and_mask_pair.first.second) {
        cerr << (mask ? "1" : "0");
    }
    cerr << endl;
    cerr << string_quality_short_to_char(qual_and_mask_pair.second.first) << endl;
    for (bool mask : qual_and_mask_pair.second.second) {
        cerr << (mask ? "1" : "0");
    }
    cerr << endl;
#endif

    
    assert(qual_and_mask_pair.first.first.size() == qual_and_mask_pair.first.second.size());
    assert(qual_and_mask_pair.second.first.size() == qual_and_mask_pair.second.second.size());
    
    // Sample our path (if dealing with source_paths)
    size_t source_path_idx = sample_path();
    string source_path;
    vg::truncated_normal_distribution<> fragment_sampler;
    if (source_path_idx != numeric_limits<size_t>::max()) {
        source_path = source_paths[source_path_idx];
#ifdef debug_ngs_sim
        cerr << "sampling from path " << source_path << " with length " << graph.get_path_length(graph.get_path_handle(source_path)) << endl;
#endif
        int64_t path_length = graph.get_path_length(graph.get_path_handle(source_path));
        fragment_sampler = vg::truncated_normal_distribution<>(fragment_mean, fragment_sd, 1.0, path_length);
    }
    else {
        fragment_sampler = vg::truncated_normal_distribution<>(fragment_mean, fragment_sd, 1.0);
    }
    int64_t fragment_length = round(fragment_sampler(prng()));
    
    
    // This is our offset along the source path, if in use
    int64_t sampled_offset;
    // And our direction to go along the source path, if in use
    bool sampled_is_reverse;
    
    int64_t walked_offset;
    
#ifdef debug_ngs_sim
    cerr << "sampled fragment length " << fragment_length << endl;
#endif
    
    if (fragment_length < transition_distrs_1.size()) {
        // the fragment is shorter than the sequencing length
        qual_and_mask_pair.first.first.resize(fragment_length);
        qual_and_mask_pair.first.second.resize(fragment_length);
        qual_and_mask_pair.second.first.resize(fragment_length);
        qual_and_mask_pair.second.second.resize(fragment_length);
#ifdef debug_ngs_sim
        cerr << "truncating reads to fragment length" << endl;
#endif
    }
    
    aln_pair.first.set_quality(qual_and_mask_pair.first.first);
    aln_pair.second.set_quality(qual_and_mask_pair.second.first);
    
    
    // reverse the quality string so that it acts like it's reading from the opposite end
    // when we walk forward from the beginning of the first read
    std::reverse(aln_pair.second.mutable_quality()->begin(),
                 aln_pair.second.mutable_quality()->end());
    
    while (!aln_pair.first.has_path() || !aln_pair.second.has_path()) {

        // Populate the sample positions
        pos_t pos;
        sample_start_pos(source_path_idx, fragment_length, sampled_offset, sampled_is_reverse, pos);
        // copy them so we can modify without losing the start pos info
        walked_offset = sampled_offset;
        bool is_reverse = sampled_is_reverse;
#ifdef debug_ngs_sim
        cerr << "sampled read 1 start pos " << pos << ", is reverse " << is_reverse << ", offset " << walked_offset << endl;
#endif
        
        // align the first end at this position on the source path or graph
        sample_read_internal(aln_pair.first, walked_offset, is_reverse, pos, source_path);
        
        if (retry_on_Ns) {
            if (aln_pair.first.sequence().find('N') != string::npos) {
#ifdef debug_ngs_sim
                cerr << "rejecting sample because of an N" << endl;
#endif
                aln_pair.first.clear_path();
                aln_pair.first.clear_sequence();
            }
        }
        
        if (!aln_pair.first.has_path()) {
            continue;
        }
        
#ifdef debug_ngs_sim
        cerr << "after first read, walked offset is " << walked_offset << endl;
#endif
        
        // walk out the unsequenced part of the fragment in the graph
        int64_t remaining_length = fragment_length - (aln_pair.first.quality().size() +
                                                      aln_pair.second.quality().size());
    
#ifdef debug_ngs_sim
        cerr << "walking " << remaining_length << " to start of next read in pair" << endl;
#endif
        if (remaining_length >= 0) {
            // we need to move forward from the end of the first read
            if (advance_by_distance(walked_offset, is_reverse, pos, remaining_length, source_path)) {
#ifdef debug_ngs_sim
                cerr << "rejecting sample because insert is off of path" << endl;
#endif
                // we hit the end of the graph trying to walk
                continue;
            }
        }
        else {
            // we need to walk backwards from the end of the first read
            if (walk_backwards(walked_offset, is_reverse, pos, -remaining_length, source_path, aln_pair.first.path())) {
#ifdef debug_ngs_sim
                cerr << "rejecting because backwards walk is off path/graph" << endl;
#endif
                continue;
            }

#ifdef debug_ngs_sim
            cerr << "walked backwards, walk length " << -remaining_length << ", is rev " << is_reverse << ", walked offset " << walked_offset << endl;
#endif
        }
#ifdef debug_ngs_sim
        cerr << "after moving for the fragment length walked offset is " << walked_offset << endl;
#endif
        // guard against running off the end of nodes
        // XXX this should not be happening
        // it seems to occur in some graphs due to the behavior of advance_by_distance
        if (vg::offset(pos) >= graph.get_length(graph.get_handle(id(pos)))) {
#ifdef debug_ngs_sim
            cerr << "rejecting sample because of invalid walked location" << endl;
#endif
            continue;
        }

        // align the second end starting at the walked position
        sample_read_internal(aln_pair.second, walked_offset, is_reverse, pos, source_path);
        
#ifdef debug_ngs_sim
        cerr << "after second read, walked offset is " << walked_offset << endl;
#endif
        
        if (retry_on_Ns) {
            if (aln_pair.second.sequence().find('N') != string::npos) {
                aln_pair.second.clear_path();
                aln_pair.second.clear_sequence();
            }
        }
    }
    
    // unreverse the second read in the pair
    reverse_complement_alignment_in_place(&aln_pair.second, [&](id_t node_id) {
        return graph.get_length(graph.get_handle(node_id));
    });
    
    for (auto aln : {&aln_pair.first, &aln_pair.second}) {
        for (const auto& mapping : aln->path().mapping()) {
            if (mapping.position().offset() < 0) {
                cerr << "error: invalid alignment!" << endl;
                cerr << pb2json(*aln) << endl;
                exit(1);
            }
        }
    }
    
    // mask out any of the sequence that we sampled to be an 'N'
    apply_N_mask(*aln_pair.first.mutable_sequence(), qual_and_mask_pair.first.second);
    apply_N_mask(*aln_pair.second.mutable_sequence(), qual_and_mask_pair.second.second);
        
    annotate_with_path_positions(aln_pair.first);
    annotate_with_path_positions(aln_pair.second);
    
    // take back the final base that we sampled
    register_sampled_position(aln_pair.first, source_path, sampled_offset + sampled_is_reverse, sampled_is_reverse);
    register_sampled_position(aln_pair.second, source_path, walked_offset + sampled_is_reverse, !sampled_is_reverse);
    
    return aln_pair;
}

void NGSSimulator::sample_read_internal(Alignment& aln, int64_t& offset, bool& is_reverse, pos_t& curr_pos,
                                        const string& source_path) {
    
    // we will accept a read that cannot be extended to the full read length if we're simulating from
    // a path that's too small or if we are sampling unsheared paths
    bool accept_partial = sample_unsheared_paths;
    if (!accept_partial && !source_path.empty()) {
        accept_partial = graph.get_path_length(graph.get_path_handle(source_path)) < transition_distrs_1.size();
    }
    
    // Make sure we are starting inside the node
    // XXX this is broken
    auto first_node_length = graph.get_length(graph.get_handle(id(curr_pos)));
    if (vg::offset(curr_pos) >= first_node_length) {
        cerr << "something wrong " << vg::offset(curr_pos) << " " << first_node_length << endl;
        cerr << vg::id(curr_pos) << ":" << vg::is_rev(curr_pos) << ":" << vg::offset(curr_pos) << endl;
    }
    assert(vg::offset(curr_pos) < first_node_length);
   
    aln.clear_path();
    aln.clear_sequence();
    
    char graph_char = graph.get_base(graph.get_handle(id(curr_pos), is_rev(curr_pos)), vg::offset(curr_pos));
    bool hit_end = false;
    
    // walk a path and generate a read sequence at the same time
    while (aln.sequence().size() < aln.quality().size() && !hit_end) {
        // sample insertion in the true graph path
        while (aln.sequence().size() < aln.quality().size() && prob_sampler(prng()) < indel_poly_rate * 0.5) {
            // TODO: no allowance for indel errors on inserted sequence
            
#ifdef debug_ngs_sim
            cerr << "insertion polymorphism at read idx " << aln.sequence().size() << ", graph pos " << curr_pos << endl;
#endif
            
            apply_insertion(aln, curr_pos);
        }
        if (aln.sequence().size() >= aln.quality().size() || hit_end) {
#ifdef debug_ngs_sim
            cerr << "break 1: ending sample with seq len " << aln.sequence().size() << ", qual len " << aln.quality().size() << ", hit end? " << hit_end << endl;
#endif
            break;
        }
        
        // sample errors
        double err_sample = prob_sampler(prng());
        double err_prob = phred_prob[aln.quality()[aln.sequence().size()]];
        while (err_sample < err_prob * indel_error_prop && !hit_end) {
            // indel errors
            if (prob_sampler(prng()) < 0.5) {
#ifdef debug_ngs_sim
                cerr << "insertion error at read idx " << aln.sequence().size() << ", graph pos " << curr_pos << endl;
#endif
                // insert error
                apply_insertion(aln, curr_pos);
                
                if (aln.sequence().size() >= aln.quality().size() || hit_end) {
#ifdef debug_ngs_sim
                    cerr << "break 2: ending sample with seq len " << aln.sequence().size() << ", qual len " << aln.quality().size() << ", hit end? " << hit_end << endl;
#endif
                    break;
                }
            }
            else {
#ifdef debug_ngs_sim
                cerr << "deletion error at read idx " << aln.sequence().size() << ", graph pos " << curr_pos << endl;
#endif
                // deletion error
                apply_deletion(aln, curr_pos);
                hit_end = advance(offset, is_reverse, curr_pos, graph_char, source_path);
            }
            
            err_sample = prob_sampler(prng());
            err_prob = phred_prob[aln.quality()[aln.sequence().size()]];
        }
        if (aln.sequence().size() >= aln.quality().size() || hit_end) {
#ifdef debug_ngs_sim
            cerr << "break 3: ending sample with seq len " << aln.sequence().size() << ", qual len " << aln.quality().size() << ", hit end? " << hit_end << endl;
#endif
            break;
        }
        
        // get the true graph char, possibly with a substitution polymorphism
        char poly_graph_char = graph_char;
        if (prob_sampler(prng()) < sub_poly_rate) {
            poly_graph_char = mutation_alphabets[poly_graph_char != 'N' ? poly_graph_char : alphabet[background_sampler(prng())]][mut_sampler(prng())];
        }
        
        // by default the read matches the true graph char
        char read_char = poly_graph_char;
        
        // sample substitution errors with the remaining err sample
        if (err_sample < err_prob) {
            // substitution error
            read_char = mutation_alphabets[read_char != 'N' ? read_char : alphabet[background_sampler(prng())]][mut_sampler(prng())];
        }
        
#ifdef debug_ngs_sim
        cerr << "aligned base at read idx " << aln.sequence().size() << ", graph pos " << curr_pos << endl;
#endif
        
        // add an aligned base (allowing errors to mask polymorphisms)
        apply_aligned_base(aln, curr_pos, graph_char, read_char);
        hit_end = advance(offset, is_reverse, curr_pos, graph_char, source_path);
        
        if (aln.sequence().size() >= aln.quality().size() || hit_end) {
#ifdef debug_ngs_sim
            cerr << "break 4: ending sample with seq len " << aln.sequence().size() << ", qual len " << aln.quality().size() << ", hit end? " << hit_end << endl;
#endif
            break;
        }
        
        // sample deletions in the true graph path
        while (prob_sampler(prng()) < indel_poly_rate * 0.5 && !hit_end) {
#ifdef debug_ngs_sim
            cerr << "deletion polymorphism at read idx " << aln.sequence().size() << ", graph pos " << curr_pos << endl;
#endif
            
            apply_deletion(aln, curr_pos);
            hit_end = advance(offset, is_reverse, curr_pos, graph_char, source_path);
        }
    }
    
    // remove the sequence and path if we hit the end the graph before finishing
    // the alignment
    if (aln.sequence().size() != aln.quality().size()) {
        if (accept_partial) {
            // we simulated the whole path, so we don't use the final quality values
            aln.mutable_quality()->resize(aln.sequence().size());
        }
        else {
            aln.clear_path();
            aln.clear_sequence();
        }
    }
    
#ifdef debug_ngs_sim
    cerr << "completed read: " << pb2json(aln) << endl;
#endif
}

bool NGSSimulator::advance(int64_t& offset, bool& is_reverse, pos_t& pos, char& graph_char, const string& source_path) {
    if (source_path.empty()) {
        return advance_on_graph(pos, graph_char);
    } else {
        return advance_on_path(offset, is_reverse, pos, graph_char, source_path);
    }
}

bool NGSSimulator::advance_on_graph(pos_t& pos, char& graph_char) {
    
    // choose a next position at random
    map<pos_t, char> next_pos_chars = algorithms::next_pos_chars(graph, pos);
    if (next_pos_chars.empty()) {
        return true;
    }
    
    vg::uniform_int_distribution<size_t> pos_distr(0, next_pos_chars.size() - 1);
    size_t next = pos_distr(prng());
    auto iter = next_pos_chars.begin();
    for (size_t i = 0; i != next; i++) {
        iter++;
    }
    pos = iter->first;
    graph_char = iter->second;
    
    return false;
}

bool NGSSimulator::advance_on_path(int64_t& offset, bool& is_reverse, pos_t& pos, char& graph_char, const string& source_path) {
    int64_t path_length = graph.get_path_length(graph.get_path_handle(source_path));
    if (is_reverse) {
        // Go left on the path
        offset--;
        if (offset < 0) {
            // We hit the end
            return true;
        }
    } else {
        // Go right on the path
        offset++;
        if (offset == path_length) {
            // We hit the end
            return true;
        }
    }
    
    // Set position according to position on path
    pos = position_at(&graph, source_path, offset, is_reverse);
    
    // And look up the character
    graph_char = graph.get_base(graph.get_handle(id(pos), is_rev(pos)), vg::offset(pos));
    
    return false;
}

bool NGSSimulator::advance_by_distance(int64_t& offset, bool& is_reverse, pos_t& pos, int64_t distance,
                                       const string& source_path) {
    if (source_path.empty()) {
        return advance_on_graph_by_distance(pos, distance);
    } else {
        return advance_on_path_by_distance(offset, is_reverse, pos, distance, source_path);
    }
}


bool NGSSimulator::advance_on_graph_by_distance(pos_t& pos, int64_t distance) {
    int64_t remaining = distance;
    handle_t handle = graph.get_handle(id(pos), is_rev(pos));
    int64_t node_length = graph.get_length(handle) - offset(pos);
    while (remaining >= node_length) {
        remaining -= node_length;
        vector<handle_t> nexts;
        graph.follow_edges(handle, false, [&](const handle_t& next) {
                nexts.push_back(next);
            });
        if (nexts.empty()) {
            return true;
        }
        size_t choice = vg::uniform_int_distribution<size_t>(0, nexts.size() - 1)(prng());
        handle = nexts[choice];
        node_length = graph.get_length(handle);
    }

    get_id(pos) = graph.get_id(handle);
    get_is_rev(pos) = graph.get_is_reverse(handle);
    get_offset(pos) += remaining;
    
    return false;
}

bool NGSSimulator::advance_on_path_by_distance(int64_t& offset, bool& is_reverse, pos_t& pos, int64_t distance,
                                               const string& source_path) {
    
    int64_t path_length = graph.get_path_length(graph.get_path_handle(source_path));
    if (is_reverse) {
        // Go left on the path
        offset -= distance;
    } else {
        // Go right on the path
        offset += distance;
    }
    if (offset < 0 || offset >= path_length) {
#ifdef debug_ngs_sim
        cerr << "walked offset of " << offset << " after advancing " << distance << ", rev ? " << is_reverse << " is outside of path of length " << path_length << endl;
#endif
        // We hit the end
        return true;
    }
    
    // Set position according to position on path
    pos = position_at(&graph, source_path, offset, is_reverse);
    
    return false;
}

bool NGSSimulator::walk_backwards_along_alignment(const Path& path, int64_t distance, pos_t& pos) {
    
    // convert to a distance forward, which is easier to implement
    int64_t remaining_to_walk = max<int64_t>(path_to_length(path) - distance, 0);
    
    for (size_t i = 0; i < path.mapping_size(); ++i) {
        const auto& mapping = path.mapping(i);
        int64_t walked_from_length = 0;
        for (size_t j = 0; j < mapping.edit_size(); ++j) {
            const auto& edit = mapping.edit(j);
            if (edit.to_length() < remaining_to_walk) {
                // we can continue to walk through this edit
                remaining_to_walk -= edit.to_length();
                walked_from_length += edit.from_length();
            }
            else {
                // this edit overlaps the place we want to walk to
                if (edit.to_length() == edit.from_length()) {
                    // this is a match/substitition, so we can walk the remaining distance
                    walked_from_length += remaining_to_walk;
                    remaining_to_walk = 0;
                }
                
                nid_t node_id = mapping.position().node_id();
                bool rev = mapping.position().is_reverse();
                size_t offset = mapping.position().offset() + walked_from_length;
                if (offset == graph.get_length(graph.get_handle(node_id))) {
                    // we're actually past-the-last on this node, which is not what we want
                    if (i + 1 < path.mapping_size()) {
                        // we can bump the position over to the next node
                        const auto& next_mapping = path.mapping(i + 1);
                        node_id = next_mapping.position().node_id();
                        rev = next_mapping.position().is_reverse();
                        offset = 0;
                    }
                    else {
                        // there's not really a "correct" option here, so we just adjust
                        // the offset by 1 and hope this doesn't lead to too much distortion
                        --offset;
                    }
                }
                pos = make_pos_t(node_id, rev, offset);
                return false;
            }
        }
    }
    return true;
}

bool NGSSimulator::walk_backwards(int64_t& offset, bool& is_reverse, pos_t& pos, int64_t distance,
                                  const string& source_path, const Path& path) {
    if (source_path.empty()) {
        return walk_backwards_along_alignment(path, distance, pos);
    }
    else {
        return advance_on_path_by_distance(offset, is_reverse, pos, -distance, source_path);
    }
}

void NGSSimulator::apply_aligned_base(Alignment& aln, const pos_t& pos, char graph_char,
                                      char read_char) {
    Path* path = aln.mutable_path();
    aln.mutable_sequence()->push_back(read_char);
    bool is_match = (graph_char == read_char);
    
    if (path->mapping_size() == 0) {
        Mapping* new_mapping = path->add_mapping();
        new_mapping->set_rank(1);
        
        Position* mapping_pos = new_mapping->mutable_position();
        mapping_pos->set_node_id(id(pos));
        mapping_pos->set_is_reverse(is_rev(pos));
        mapping_pos->set_offset(offset(pos));
        
        Edit* new_edit = new_mapping->add_edit();
        new_edit->set_from_length(1);
        new_edit->set_to_length(1);
        if (!is_match) {
            new_edit->mutable_sequence()->push_back(read_char);
        }
    }
    else {
        Mapping* last_mapping = path->mutable_mapping(path->mapping_size() - 1);
        if (last_mapping->position().node_id() == id(pos) &&
            last_mapping->position().is_reverse() == is_rev(pos)) {
            
            Edit* last_edit = last_mapping->mutable_edit(last_mapping->edit_size() - 1);
            if (last_edit->from_length() > 0 && last_edit->to_length() > 0) {
                if (last_edit->sequence().size() > 0 && !is_match) {
                    last_edit->set_from_length(last_edit->from_length() + 1);
                    last_edit->set_to_length(last_edit->to_length() + 1);
                    last_edit->mutable_sequence()->push_back(read_char);
                }
                else if (last_edit->sequence().size() == 0 && is_match) {
                    last_edit->set_from_length(last_edit->from_length() + 1);
                    last_edit->set_to_length(last_edit->to_length() + 1);
                }
                else {
                    Edit* new_edit = last_mapping->add_edit();
                    new_edit->set_from_length(1);
                    new_edit->set_to_length(1);
                    if (!is_match) {
                        new_edit->mutable_sequence()->push_back(read_char);
                    }
                }
            }
            else {
                Edit* new_edit = last_mapping->add_edit();
                new_edit->set_from_length(1);
                new_edit->set_to_length(1);
                if (!is_match) {
                    new_edit->mutable_sequence()->push_back(read_char);
                }
            }
        }
        else {
            Mapping* new_mapping = path->add_mapping();
            new_mapping->set_rank(last_mapping->rank() + 1);
            
            Position* mapping_pos = new_mapping->mutable_position();
            mapping_pos->set_node_id(id(pos));
            mapping_pos->set_is_reverse(is_rev(pos));
            
            Edit* new_edit = new_mapping->add_edit();
            new_edit->set_from_length(1);
            new_edit->set_to_length(1);
            if (!is_match) {
                new_edit->mutable_sequence()->push_back(read_char);
            }
        }
    }
}

void NGSSimulator::apply_deletion(Alignment& aln, const pos_t& pos) {
    Path* path = aln.mutable_path();
    if (path->mapping_size() == 0) {
        // don't introduce a deletion at the beginning of a read
        // TODO: check for deletions at the end of a read?
        return;
    }
    
    Mapping* last_mapping = path->mutable_mapping(path->mapping_size() - 1);
    if (last_mapping->position().node_id() == id(pos) &&
        last_mapping->position().is_reverse() == is_rev(pos)) {
        
        Edit* last_edit = last_mapping->mutable_edit(last_mapping->edit_size() - 1);
        if (last_edit->from_length() > 0 && last_edit->to_length() == 0) {
            last_edit->set_from_length(last_edit->from_length() + 1);
        }
        else {
            Edit* new_edit = last_mapping->add_edit();
            new_edit->set_from_length(1);
        }
    }
    else {
        Mapping* new_mapping = path->add_mapping();
        new_mapping->set_rank(last_mapping->rank() + 1);
        
        Position* mapping_pos = new_mapping->mutable_position();
        mapping_pos->set_node_id(id(pos));
        mapping_pos->set_is_reverse(is_rev(pos));
        
        Edit* new_edit = new_mapping->add_edit();
        new_edit->set_from_length(1);
    }
}

void NGSSimulator::apply_insertion(Alignment& aln, const pos_t& pos) {
    Path* path = aln.mutable_path();
    char insert_char = alphabet[background_sampler(prng())];
    aln.mutable_sequence()->push_back(insert_char);
    
    if (path->mapping_size() == 0) {
        Mapping* new_mapping = path->add_mapping();
        new_mapping->set_rank(1);
        
        Position* mapping_pos = new_mapping->mutable_position();
        mapping_pos->set_node_id(id(pos));
        mapping_pos->set_is_reverse(is_rev(pos));
        mapping_pos->set_offset(offset(pos));
        
        Edit* new_edit = new_mapping->add_edit();
        new_edit->set_to_length(1);
        new_edit->set_sequence(string(1, insert_char));
    }
    else {
        Mapping* last_mapping = path->mutable_mapping(path->mapping_size() - 1);
        Edit* last_edit = last_mapping->mutable_edit(last_mapping->edit_size() - 1);
        if (last_edit->from_length() == 0 && last_edit->to_length() > 0) {
            last_edit->set_to_length(last_edit->to_length() + 1);
            last_edit->mutable_sequence()->push_back(insert_char);
        }
        else {
            Edit* new_edit = last_mapping->add_edit();
            new_edit->set_to_length(1);
            new_edit->set_sequence(string(1, insert_char));
        }
    }
}

size_t NGSSimulator::sample_path() {
    if (source_paths.empty()) {
        return numeric_limits<size_t>::max();
    }
    else {
        size_t path_idx = path_sampler(prng());
        return path_idx;
    }
}

void NGSSimulator::sample_start_pos(const size_t& source_path_idx, const int64_t& fragment_length,
                                    int64_t& offset, bool& is_reverse, pos_t& pos) {
    if (source_paths.empty()) {
        pos = sample_start_graph_pos();
        offset = 0;
        is_reverse = false;
    }
    else {
        tie(offset, is_reverse, pos) = sample_start_path_pos(source_path_idx, fragment_length);
    }
}

pos_t NGSSimulator::sample_start_graph_pos() {
    // The start pos sampler has been set up in graph space, 1-based
    assert(start_pos_samplers.size() == 1);
    size_t idx = start_pos_samplers[0](prng());
    
    id_t id = dynamic_cast<VectorizableHandleGraph&>(graph).node_at_vector_offset(idx);
    bool rev = strand_sampler(prng());
    size_t node_offset = idx - dynamic_cast<VectorizableHandleGraph&>(graph).node_vector_offset(id) - 1;
    
    return make_pos_t(id, rev, node_offset);
}

tuple<int64_t, bool, pos_t> NGSSimulator::sample_start_path_pos(const size_t& source_path_idx,
                                                                const int64_t& fragment_length) {
    
    int64_t path_length = graph.get_path_length(graph.get_path_handle(source_paths[source_path_idx]));
    bool rev = strand_sampler(prng());
#ifdef debug_ngs_sim
    cerr << "sampling start position on path " << source_paths[source_path_idx] << ", strand " << rev << ", path length " << path_length << endl;
#endif
    int64_t offset;
    if (sample_unsheared_paths || path_length < transition_distrs_1.size() ||
        (fragment_length > 0 && fragment_length >= path_length)) {
        if (rev) {
            offset = path_length - 1;
        }
        else {
            offset = 0;
        }
    }
    else {
        
        // we'll not let it choose too unreasonable of a start position so that we can speed
        // up the process of continuously resampling impractically close to the end of paths
        bool feasible = false;
        int64_t shortened_fragment_length = fragment_length - ceil(sqrt(fragment_length));
        do {
            // The start pos sampler has been set up in path space, 0-based
            offset = start_pos_samplers[source_path_idx](prng());
            if (rev) {
                feasible = (offset - shortened_fragment_length >= -1);
            }
            else {
                feasible = (offset + shortened_fragment_length <= path_length);
            }
        } while (!feasible);
    }
    pos_t pos = position_at(&graph, source_paths[source_path_idx], offset, rev);
    
    return make_tuple(offset, rev, pos);
}

string NGSSimulator::get_read_name() {
    stringstream sstrm;
    size_t num;
#pragma omp atomic capture
    num = sample_counter++;
    sstrm << "seed_" << seed << "_fragment_" << num;
    return sstrm.str();
}

void NGSSimulator::record_read_quality(const Alignment& aln, bool read_2) {
    const string& quality = aln.quality();
    const string& sequence = aln.sequence();
    assert(sequence.size() == quality.size());
    auto& transition_distrs = read_2 ? transition_distrs_2 : transition_distrs_1;
    if (quality.empty()) {
        return;
    }
    while (transition_distrs.size() < quality.size()) {
        // coding-time random engine to perturb the seed for each position
        linear_congruential_engine<uint64_t, 16793141576979709161ull, 876936395080740889ull, 18446744073709551557ull>  seed_perturbor(seed + transition_distrs.size() + read_2);
        transition_distrs.emplace_back(seed ? seed_perturbor() : random_device()());
    }
    // record the initial quality and N-mask
    transition_distrs[0].record_transition(pair<uint8_t, bool>(0, false),
                                           pair<uint8_t, bool>(quality[0], sequence[0] == 'N'));
    // record the subsequent quality and N-mask transitions
    for (size_t i = 1; i < quality.size(); i++) {
        transition_distrs[i].record_transition(pair<uint8_t, bool>(quality[i - 1], sequence[i - 1] == 'N'),
                                               pair<uint8_t, bool>(quality[i], sequence[i] == 'N'));
    }
}
    
void NGSSimulator::record_read_pair_quality(const Alignment& aln_1, const Alignment& aln_2) {
    // record the transitions within the reads separates
    record_read_quality(aln_1, false);
    record_read_quality(aln_2, true);
    // record the joint distribution of the first quality and N-mask
    if (!aln_1.quality().empty() && !aln_2.quality().empty()) {
        joint_initial_distr.record_transition(pair<uint8_t, bool>(0, false),
                                              make_pair(pair<uint8_t, bool>(aln_1.quality()[0], aln_1.sequence()[0] == 'N'),
                                                        pair<uint8_t, bool>(aln_2.quality()[0], aln_2.sequence()[0] == 'N')));
    }
}

void NGSSimulator::finalize() {
    for (auto& markov_distr : transition_distrs_1) {
        markov_distr.finalize();
    }
    for (auto& markov_distr : transition_distrs_2) {
        markov_distr.finalize();
    }
    joint_initial_distr.finalize();
}

pair<string, vector<bool>> NGSSimulator::sample_read_quality() {
    // only use the first trained distribution (on the assumption that it better reflects the properties of
    // single-ended sequencing)
    return sample_read_quality_internal(transition_distrs_1[0].sample_transition(pair<uint8_t, bool>(0, false)),
                                        true);
}
    
pair<pair<string, vector<bool>>, pair<string, vector<bool>>> NGSSimulator::sample_read_quality_pair() {
    if (transition_distrs_2.empty()) {
        // no paired training data, sample qual strings independently
        return make_pair(sample_read_quality(), sample_read_quality());
    }
    else {
        // paired training data, sample the start quality jointly
        auto first_quals_and_masks = joint_initial_distr.sample_transition(pair<uint8_t, bool>(0, false));
        return make_pair(sample_read_quality_internal(first_quals_and_masks.first, true),
                         sample_read_quality_internal(first_quals_and_masks.second, false));
    }
}
    
                                
pair<string, vector<bool>> NGSSimulator::sample_read_quality_internal(pair<uint8_t, bool> first,
                                                                      bool transitions_1) {
    
    auto& transition_distrs = transitions_1 ? transition_distrs_1 : transition_distrs_2;
    string quality(transition_distrs.size(), first.first);
    vector<bool> n_masks(transition_distrs.size(), first.second);
    pair<uint8_t, bool> at = first;
    for (size_t i = 1; i < transition_distrs.size(); i++) {
        at = transition_distrs[i].sample_transition(at);
        quality[i] = at.first;
        n_masks[i] = at.second;
    }
    return make_pair(quality, n_masks);
}
                                              
void NGSSimulator::apply_N_mask(string& sequence, const vector<bool>& n_mask) {
    for (size_t i = 0; i < sequence.size(); i++) {
        if (n_mask[i]) {
            sequence[i] = 'N';
        }
    }
}

}
