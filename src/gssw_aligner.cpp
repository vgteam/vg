#include "gssw_aligner.hpp"
#include "json2pb.h"

// log(10)
static const double quality_scale_factor = 10.0 / log(10.0);
static const double exp_overflow_limit = log(std::numeric_limits<double>::max());

using namespace vg;
using namespace std;

BaseAligner::~BaseAligner(void) {
    free(nt_table);
    free(score_matrix);
}

gssw_graph* BaseAligner::create_gssw_graph(Graph& g) {
    
    // add a dummy sink node if we're pinning
    gssw_graph* graph = gssw_graph_create(g.node_size());
    unordered_map<int64_t, gssw_node*> nodes;
    
    for (int i = 0; i < g.node_size(); ++i) {
        Node* n = g.mutable_node(i);
        // switch any non-ATGCN characters from the node sequence to N
        auto cleaned_seq = nonATGCNtoN(n->sequence());
        gssw_node* node = (gssw_node*)gssw_node_create(n, n->id(),
                                                       cleaned_seq.c_str(),
                                                       nt_table,
                                                       score_matrix);
        nodes[n->id()] = node;
        gssw_graph_add_node(graph, node);
    }
    
    for (int i = 0; i < g.edge_size(); ++i) {
        // Convert all the edges
        Edge* e = g.mutable_edge(i);
        if(!e->from_start() && !e->to_end()) {
            // This is a normal end to start edge.
            gssw_nodes_add_edge(nodes[e->from()], nodes[e->to()]);
        } else if(e->from_start() && e->to_end()) {
            // This is a start to end edge, but isn't reversing and can be converted to a normal end to start edge.
            
            // Flip the start and end
            gssw_nodes_add_edge(nodes[e->to()], nodes[e->from()]);
        } else {
            // TODO: It's a reversing edge, which gssw doesn't support yet. What
            // we should really do is do a topological sort to break cycles, and
            // then flip everything at the lower-rank end of this edge around,
            // so we don't have to deal with its reversing-ness. But for now we
            // just die so we don't get nonsense into gssw.
#pragma omp critical
            {
                // We need the critical section so we don't throw uncaught
                // exceptions in multiple threads at once, leading to C++ trying
                // to run termiante in parallel. This doesn't make it safe, just
                // slightly safer.
                cerr << "Can't gssw over reversing edge " <<e->from() << (e->from_start() ? " start" : " end") << " -> "
                << e->to() << (e->to_end() ? " end" : " start")  << endl;
                // TODO: there's no safe way to kill the program without a way
                // to signal the master to do it, via a shared variable in the
                // clause that made us parallel.
            }
            exit(1);
        }
    }
    
    return graph;
    
}



void BaseAligner::gssw_mapping_to_alignment(gssw_graph* graph,
                                            gssw_graph_mapping* gm,
                                            Alignment& alignment,
                                            bool pinned,
                                            bool pin_left,
                                            bool print_score_matrices) {    
    alignment.clear_path();
    alignment.set_score(gm->score);
    alignment.set_query_position(0);
    Path* path = alignment.mutable_path();
    //alignment.set_cigar(graph_cigar(gm));
    
    gssw_graph_cigar* gc = &gm->cigar;
    gssw_node_cigar* ncs = gc->elements;
    //cerr << "gm->position " << gm->position << endl;
    string& to_seq = *alignment.mutable_sequence();
    //cerr << "-------------" << endl;
    
    if (print_score_matrices) {
        gssw_graph_print_score_matrices(graph, to_seq.c_str(), to_seq.size(), stderr);
        //cerr << alignment.DebugString() << endl;
    }
    
    int to_pos = 0;
    int from_pos = gm->position;
    
    for (int i = 0; i < gc->length; ++i) {
        // check that the current alignment has a non-zero length
        gssw_cigar* c = ncs[i].cigar;
        int l = c->length;
        if (l == 0) continue;
        gssw_cigar_element* e = c->elements;
        
        Node* from_node = (Node*) ncs[i].node->data;
        string& from_seq = *from_node->mutable_sequence();
        Mapping* mapping = path->add_mapping();
        
        if (i > 0) {
            // reset for each node after the first
            from_pos = 0;
        }
        
        mapping->mutable_position()->set_node_id(ncs[i].node->id);
        mapping->mutable_position()->set_offset(from_pos);
        mapping->set_rank(path->mapping_size());
        
        //cerr << from_node->id() << ":" << endl;
        
        for (int j=0; j < l; ++j, ++e) {
            int32_t length = e->length;
            //cerr << e->length << e->type << endl;
            
            Edit* edit;
            switch (e->type) {
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
                            if (h - last_start > 0) {
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
                    if (h - last_start > 0) {
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
                    cerr << "error:[Aligner::gssw_mapping_to_alignment] "
                    << "unsupported cigar op type " << e->type << endl;
                    exit(1);
                    break;
                    
            }
        }
    }
    
    // compute and set identity
    alignment.set_identity(identity(alignment.path()));
}

void BaseAligner::reverse_graph(Graph& g, Graph& reversed_graph_out) {
    if (reversed_graph_out.node_size()) {
        cerr << "error:[Aligner::reverse_graph] output graph is not empty" << endl;
        exit(EXIT_FAILURE);
    }
    
    // add reversed nodes in reverse order (Graphs come in topologically sorted and gssw
    // depends on this fact)
    for (int64_t i = g.node_size() - 1; i >= 0; i--) {
        const Node& original_node = g.node(i);
        
        Node* reversed_node = reversed_graph_out.add_node();
        
        // reverse the sequence
        string* reversed_node_seq = reversed_node->mutable_sequence();
        reversed_node_seq->resize(original_node.sequence().length());
        reverse_copy(original_node.sequence().begin(), original_node.sequence().end(), reversed_node_seq->begin());
        
        // preserve ids for easier translation
        reversed_node->set_id(original_node.id());
    }
    
    // add reversed edges
    for (int64_t i = 0; i < g.edge_size(); i++) {
        const Edge& original_edge = g.edge(i);
        
        Edge* reversed_edge = reversed_graph_out.add_edge();
        
        // reverse edge orientation
        reversed_edge->set_from(original_edge.to());
        reversed_edge->set_to(original_edge.from());
        
        // we will swap the 5'/3' labels of the node ends after reversing the sequence so that
        // an edge leaving an end now enters a beginning and an edge entering a beginning now
        // leaves an end
        reversed_edge->set_from_start(original_edge.to_end());
        reversed_edge->set_to_end(original_edge.from_start());
    }
    
}

void BaseAligner::unreverse_graph(Graph& graph) {
    // this is only for getting correct reference-relative edits, so we can get away with only
    // reversing the sequences and not paying attention to the edges
    for (int64_t i = 0; i < graph.node_size(); i++) {
        Node* node = graph.mutable_node(i);
        string* node_seq = node->mutable_sequence();
        reverse(node_seq->begin(), node_seq->end());
    }
}

void BaseAligner::unreverse_graph_mapping(gssw_graph_mapping* gm) {
    
    gssw_graph_cigar* graph_cigar = &(gm->cigar);
    gssw_node_cigar* node_cigars = graph_cigar->elements;
    
    // reverse the order of the node cigars
    int32_t num_switching_nodes = graph_cigar->length / 2;
    int32_t last_idx = graph_cigar->length - 1;
    for (int32_t i = 0; i < num_switching_nodes; i++) {
        swap(node_cigars[i], node_cigars[last_idx - i]);
    }
    
    // reverse the actual cigar string for each node cigar
    for (int32_t i = 0; i < graph_cigar->length; i++) {
        gssw_cigar* node_cigar = node_cigars[i].cigar;
        gssw_cigar_element* elements = node_cigar->elements;
        
        int32_t num_switching_elements = node_cigar->length / 2;
        last_idx = node_cigar->length - 1;
        for (int32_t j = 0; j < num_switching_elements; j++) {
            swap(elements[j], elements[last_idx - j]);
        }
    }
    
    // compute the position in the first node
    if (graph_cigar->length > 0) {
        gssw_cigar_element* first_node_elements = node_cigars[0].cigar->elements;
        int32_t num_first_node_elements = node_cigars[0].cigar->length;
        uint32_t num_ref_aligned = 0; // the number of characters on the node sequence that are aligned
        for (int32_t i = 0; i < num_first_node_elements; i++) {
            switch (first_node_elements[i].type) {
                case 'M':
                case 'X':
                case 'N':
                case 'D':
                    num_ref_aligned += first_node_elements[i].length;
                    break;
                    
            }
        }
        gm->position = node_cigars[0].node->len - num_ref_aligned - (graph_cigar->length == 1 ? gm->position : 0);
    }
    else {
        gm->position = 0;
    }
}

string BaseAligner::graph_cigar(gssw_graph_mapping* gm) {
    stringstream s;
    gssw_graph_cigar* gc = &gm->cigar;
    gssw_node_cigar* nc = gc->elements;
    int to_pos = 0;
    int from_pos = gm->position;
    //string& to_seq = *alignment.mutable_sequence();
    s << from_pos << '@';
    for (int i = 0; i < gc->length; ++i, ++nc) {
        if (i > 0) from_pos = 0; // reset for each node after the first
        Node* from_node = (Node*) nc->node->data;
        s << from_node->id() << ':';
        gssw_cigar* c = nc->cigar;
        int l = c->length;
        gssw_cigar_element* e = c->elements;
        for (int j=0; j < l; ++j, ++e) {
            s << e->length << e->type;
        }
        if (i + 1 < gc->length) {
            s << ",";
        }
    }
    return s.str();
}

void BaseAligner::init_mapping_quality(double gc_content) {
    log_base = gssw_dna_recover_log_base(match, mismatch, gc_content, 1e-12);
}

int32_t BaseAligner::score_gap(size_t gap_length) {
    return gap_length ? -gap_open - (gap_length - 1) * gap_extension : 0;
}

double BaseAligner::maximum_mapping_quality_exact(vector<double>& scaled_scores, size_t* max_idx_out) {
    size_t size = scaled_scores.size();
    
    // if necessary, assume a null alignment of 0.0 for comparison since this is local
    if (size == 1) {
        scaled_scores.push_back(0.0);
    }
    
    double max_score = scaled_scores[0];
    size_t max_idx = 0;
    for (size_t i = 1; i < size; i++) {
        if (scaled_scores[i] > max_score) {
            max_score = scaled_scores[i];
            max_idx = i;
        }
    }
    
    *max_idx_out = max_idx;

    double qual = 0;
    if (max_score * size < exp_overflow_limit) {
        // no risk of double overflow, sum exp directly (half as many transcendental function evals)
        double numer = 0.0;
        for (size_t i = 0; i < size; i++) {
            if (i == max_idx) {
                continue;
            }
            numer += exp(scaled_scores[i]);
        }
        qual = -10.0 * log10(numer / (numer + exp(scaled_scores[max_idx])));
    }
    else {
        // work in log transformed valued to avoid risk of overflow
        double log_sum_exp = scaled_scores[0];
        for (size_t i = 1; i < size; i++) {
            log_sum_exp = add_log(log_sum_exp, scaled_scores[i]);
        }
        qual = -10.0 * log10(1.0 - exp(scaled_scores[max_idx] - log_sum_exp));
    }

    return qual;
}

// TODO: this algorithm has numerical problems that would be difficult to solve without increasing the
// time complexity: adding the probability of the maximum likelihood tends to erase the contribution
// of the other terms so that when you subtract them off you get scores of 0 or infinity

//vector<double> Aligner::all_mapping_qualities_exact(vector<double>& scaled_scores) {
//
//    double max_score = *max_element(scaled_scores.begin(), scaled_scores.end());
//    size_t size = scaled_scores.size();
//
//    vector<double> mapping_qualities(size);
//
//    if (max_score * size < exp_overflow_limit) {
//        // no risk of double overflow, sum exp directly (half as many transcendental function evals)
//        vector<double> exp_scaled_scores(size);
//        for (size_t i = 0; i < size; i++) {
//            exp_scaled_scores[i] = exp(scaled_scores[i]);
//        }
//        double denom = std::accumulate(exp_scaled_scores.begin(), exp_scaled_scores.end(), 0.0);
//        for (size_t i = 0; i < size; i++) {
//            mapping_qualities[i] = -10.0 * log10((denom - exp_scaled_scores[i]) / denom);
//        }
//    }
//    else {
//        // work in log transformed valued to avoid risk of overflow
//        double log_sum_exp = scaled_scores[0];
//        for (size_t i = 1; i < size; i++) {
//            log_sum_exp = add_log(log_sum_exp, scaled_scores[i]);
//        }
//        for (size_t i = 0; i < size; i++) {
//            mapping_qualities[i] = -10.0 * log10(1.0 - exp(scaled_scores[i] - log_sum_exp));
//        }
//    }
//    return mapping_qualities;
//}

double BaseAligner::maximum_mapping_quality_approx(vector<double>& scaled_scores, size_t* max_idx_out) {
    
    // if necessary, assume a null alignment of 0.0 for comparison since this is local
    if (scaled_scores.size() == 1) {
        scaled_scores.push_back(0.0);
    }
    
    double max_score = scaled_scores[0];
    size_t max_idx = 0;
    
    double next_score = -std::numeric_limits<double>::max();
    int32_t next_count = 0;
    
    for (int32_t i = 1; i < scaled_scores.size(); i++) {
        double score = scaled_scores[i];
        if (score > max_score) {
            if (next_score == max_score) {
                next_count++;
            }
            else {
                next_score = max_score;
                next_count = 1;
            }
            max_score = score;
            max_idx = i;
        }
        else if (score > next_score) {
            next_score = score;
            next_count = 1;
        }
        else if (score == next_score) {
            next_count++;
        }
    }
    
    *max_idx_out = max_idx;

    return max(0.0, quality_scale_factor * (max_score - next_score - (next_count > 1 ? log(next_count) : 0.0)));
}

void BaseAligner::compute_mapping_quality(vector<Alignment>& alignments,
                                          int max_mapping_quality,
                                          bool fast_approximation,
                                          double cluster_mq,
                                          bool use_cluster_mq,
                                          int overlap_count,
                                          double mq_estimate,
                                          double identity_weight) {
    
    if (log_base <= 0.0) {
        cerr << "error:[Aligner] must call init_mapping_quality before computing mapping qualities" << endl;
        exit(EXIT_FAILURE);
    }
    
    if (alignments.empty()) {
        return;
    }
    
    vector<double> scaled_scores(alignments.size());
    for (size_t i = 0; i < alignments.size(); i++) {
        scaled_scores[i] = log_base * alignments[i].score();
    }

    double mapping_quality;
    size_t max_idx;
    if (!fast_approximation) {
        mapping_quality = maximum_mapping_quality_exact(scaled_scores, &max_idx);
    }
    else {
        mapping_quality = maximum_mapping_quality_approx(scaled_scores, &max_idx);
    }

    if (use_cluster_mq) {
        mapping_quality = prob_to_phred(sqrt(phred_to_prob(cluster_mq + mapping_quality)));
    }

    if (mq_estimate < mapping_quality) {
        mapping_quality = prob_to_phred(sqrt(phred_to_prob(mq_estimate + mapping_quality)));
    }

    auto& max_aln = alignments[max_idx];
    int l = max(alignment_to_length(max_aln), alignment_from_length(max_aln));
    double identity = 1. - (double)(l * match - max_aln.score()) / (match + mismatch) / l;

    mapping_quality *= pow(identity, identity_weight);

    if (mapping_quality > max_mapping_quality) {
        mapping_quality = max_mapping_quality;
    }

    if (alignments[max_idx].score() == 0) {
        mapping_quality = 0;
    }

    alignments[max_idx].set_mapping_quality(max(0, (int32_t) round(mapping_quality)));
    for (int i = 1; i < alignments.size(); ++i) {
        alignments[0].add_secondary_score(alignments[i].score());
    }
}

int32_t BaseAligner::compute_mapping_quality(vector<int32_t> scores, bool fast_approximation) {
    
    vector<double> scaled_scores(scores.size(), 0.0);
    for (size_t i = 0; i < scores.size(); i++) {
        scaled_scores[i] = log_base * scores[i];
    }
    size_t idx;
    return (int32_t) (fast_approximation ? maximum_mapping_quality_approx(scaled_scores, &idx)
                                         : maximum_mapping_quality_exact(scaled_scores, &idx));
}

void BaseAligner::compute_paired_mapping_quality(pair<vector<Alignment>, vector<Alignment>>& alignment_pairs,
                                                 const vector<double>& frag_weights,
                                                 int max_mapping_quality1,
                                                 int max_mapping_quality2,
                                                 bool fast_approximation,
                                                 double cluster_mq,
                                                 bool use_cluster_mq,
                                                 int overlap_count1,
                                                 int overlap_count2,
                                                 double mq_estimate1,
                                                 double mq_estimate2,
                                                 double identity_weight) {
    
    if (log_base <= 0.0) {
        cerr << "error:[Aligner] must call init_mapping_quality before computing mapping qualities" << endl;
        exit(EXIT_FAILURE);
    }
    
    size_t size = min(
                      alignment_pairs.first.size(),
                      alignment_pairs.second.size());
    
    if (size == 0) {
        return;
    }
    
    vector<double> scaled_scores(size);
    
    for (size_t i = 0; i < size; i++) {
        scaled_scores[i] = log_base * (alignment_pairs.first[i].score() + alignment_pairs.second[i].score());
        // + frag_weights[i]);
        // ^^^ we could also incorporate the fragment weights, but this does not seem to help performance
        // at least with the weights which we are using; todo explore this
    }

    size_t max_idx;
    double mapping_quality;
    if (!fast_approximation) {
        mapping_quality = maximum_mapping_quality_exact(scaled_scores, &max_idx);
    }
    else {
        mapping_quality = maximum_mapping_quality_approx(scaled_scores, &max_idx);
    }
    
    if (use_cluster_mq) {
        mapping_quality = prob_to_phred(sqrt(phred_to_prob(cluster_mq + mapping_quality)));
    }

    double mapping_quality1 = mapping_quality;
    double mapping_quality2 = mapping_quality;

    if (mq_estimate1 < mapping_quality2) {
        mapping_quality1 = prob_to_phred(sqrt(phred_to_prob(mq_estimate1 + mapping_quality1)));
    }
    if (mq_estimate2 < mapping_quality2) {
        mapping_quality2 = prob_to_phred(sqrt(phred_to_prob(mq_estimate2 + mapping_quality2)));
    }

    auto& max_aln1 = alignment_pairs.first[max_idx];
    int len1 = max(alignment_to_length(max_aln1), alignment_from_length(max_aln1));
    double identity1 = 1. - (double)(len1 * match - max_aln1.score()) / (match + mismatch) / len1;
    auto& max_aln2 = alignment_pairs.second[max_idx];
    int len2 = max(alignment_to_length(max_aln2), alignment_from_length(max_aln2));
    double identity2 = 1. - (double)(len2 * match - max_aln2.score()) / (match + mismatch) / len2;

    mapping_quality1 *= pow(identity1, identity_weight);
    mapping_quality2 *= pow(identity2, identity_weight);

    if (mapping_quality1 > max_mapping_quality1) {
        mapping_quality1 = max_mapping_quality1;
    }
    if (mapping_quality2 > max_mapping_quality2) {
        mapping_quality2 = max_mapping_quality2;
    }

    if (alignment_pairs.first[max_idx].score() == 0) {
        mapping_quality1 = 0;
    }
    if (alignment_pairs.second[max_idx].score() == 0) {
        mapping_quality2 = 0;
    }

    mapping_quality = max(0, (int32_t)round(min(mapping_quality1, mapping_quality2)));

    alignment_pairs.first[max_idx].set_mapping_quality(mapping_quality);
    alignment_pairs.second[max_idx].set_mapping_quality(mapping_quality);

    for (int i = 1; i < alignment_pairs.first.size(); ++i) {
        alignment_pairs.first[0].add_secondary_score(alignment_pairs.first[i].score());
    }
    for (int i = 1; i < alignment_pairs.second.size(); ++i) {
        alignment_pairs.second[0].add_secondary_score(alignment_pairs.second[i].score());
    }

}

double BaseAligner::mapping_quality_score_diff(double mapping_quality) const {
    return mapping_quality / (quality_scale_factor * log_base);
}

double BaseAligner::estimate_next_best_score(int length, double min_diffs) {
    return ((length - min_diffs) * match - min_diffs * mismatch);
}

double BaseAligner::estimate_max_possible_mapping_quality(int length, double min_diffs, double next_min_diffs) {
    double max_score = log_base * ((length - min_diffs) * match - min_diffs * mismatch);
    double next_max_score = log_base * ((length - next_min_diffs) * match - next_min_diffs * mismatch);
    vector<double> v = { max_score, next_max_score };
    size_t max_idx;
    return maximum_mapping_quality_approx(v, &max_idx);
}

double BaseAligner::score_to_unnormalized_likelihood_ln(double score) {
    // Log base needs to be set, or this can't work. It's set by default in
    // QualAdjAligner but needs to be set up manually in the normal Aligner.
    assert(log_base != 0);
    // Likelihood is proportional to e^(lambda * score), so ln is just the exponent.
    return log_base * score;
}

size_t BaseAligner::longest_detectable_gap(const Alignment& alignment, const string::const_iterator& read_pos) const {
    // algebraic solution for when score is > 0 assuming perfect match other than gap
    int64_t overhang_length = min(read_pos - alignment.sequence().begin(), alignment.sequence().end() - read_pos);
    int64_t numer = match * overhang_length + full_length_bonus;
    int64_t gap_length = (numer - gap_open) / gap_extension + 1;
    return gap_length >= 0 && overhang_length > 0 ? gap_length : 0;
}

size_t BaseAligner::longest_detectable_gap(const Alignment& alignment) const {
    // longest detectable gap across entire read is in the middle
    return longest_detectable_gap(alignment, alignment.sequence().begin() + (alignment.sequence().size() / 2));
    
}

int32_t BaseAligner::score_alignment(const Alignment& aln, const function<size_t(pos_t, pos_t, size_t)>& estimate_distance,
    bool strip_bonuses) {
    
    int score = 0;
    int read_offset = 0;
    auto& path = aln.path();
    for (int i = 0; i < path.mapping_size(); ++i) {
        // For each mapping
        auto& mapping = path.mapping(i);
        for (int j = 0; j < mapping.edit_size(); ++j) {
            // For each edit in the mapping
            auto& edit = mapping.edit(j);
            
            // Score the edit according to its type
            if (edit_is_match(edit)) {
                score += score_exact_match(aln, read_offset, edit.to_length());
            } else if (edit_is_sub(edit)) {
                score -= mismatch * edit.sequence().size();
            } else if (edit_is_deletion(edit)) {
                score -= edit.from_length() ? gap_open + (edit.from_length() - 1) * gap_extension : 0;
            } else if (edit_is_insertion(edit)
                       && !((i == 0 && j == 0)
                            || (i == path.mapping_size()-1
                                && j == mapping.edit_size()-1))) {
                // todo how do we score this qual adjusted?
                score -= edit.to_length() ? gap_open + (edit.to_length() - 1) * gap_extension : 0;
            }
            read_offset += edit.to_length();
        }
        // score any intervening gaps in mappings using approximate distances
        if (i+1 < path.mapping_size()) {
            // what is the distance between the last position of this mapping
            // and the first of the next
            Position last_pos = mapping.position();
            last_pos.set_offset(last_pos.offset() + mapping_from_length(mapping));
            Position next_pos = path.mapping(i+1).position();
            // Estimate the distance
            int dist = estimate_distance(make_pos_t(last_pos), make_pos_t(next_pos), aln.sequence().size());
            if (dist > 0) {
                // If it's nonzero, score it as a deletion gap
                score -= gap_open + (dist - 1) * gap_extension;
            }
        }
    }
    
    if (!strip_bonuses) {
        // We should report any bonuses used in the DP in the final score
        if (!softclip_start(aln)) {
            score += full_length_bonus;
        }
        if (!softclip_end(aln)) {
            score += full_length_bonus;
        }
    }
    
    return score;
}

int32_t BaseAligner::remove_bonuses(const Alignment& aln, bool pinned, bool pin_left) {
    int32_t score = aln.score();
    if (softclip_start(aln) == 0 && !(pinned && pin_left)) {
        // No softclip at the start, and a left end bonus was applied.
        score -= full_length_bonus;
    }
    if (softclip_end(aln) == 0 && !(pinned && !pin_left)) {
        // No softclip at the end, and a right end bonus was applied.
        score -= full_length_bonus;
    }
    return score;
}


Aligner::Aligner(int8_t _match,
                 int8_t _mismatch,
                 int8_t _gap_open,
                 int8_t _gap_extension,
                 int8_t _full_length_bonus,
                 double gc_content)
{
    match = _match;
    mismatch = _mismatch;
    gap_open = _gap_open;
    gap_extension = _gap_extension;
    full_length_bonus = _full_length_bonus;
    // these are used when setting up the nodes
    nt_table = gssw_create_nt_table();
    score_matrix = gssw_create_score_matrix(match, mismatch);
    BaseAligner::init_mapping_quality(gc_content);
}


void Aligner::align_internal(Alignment& alignment, vector<Alignment>* multi_alignments, Graph& g,
                             bool pinned, bool pin_left, int32_t max_alt_alns, bool print_score_matrices) {

    // check input integrity
    if (pin_left && !pinned) {
        cerr << "error:[Aligner] cannot choose pinned end in non-pinned alignment" << endl;
        exit(EXIT_FAILURE);
    }
    if (multi_alignments && !pinned) {
        cerr << "error:[Aligner] multiple traceback is not implemented in local alignment, only pinned and global" << endl;
        exit(EXIT_FAILURE);
    }
    if (!multi_alignments && max_alt_alns != 1) {
        cerr << "error:[Aligner] cannot specify maximum number of alignments in single alignment" << endl;
        exit(EXIT_FAILURE);
    }
    
    // alignment pinning algorithm is based on pinning in bottom right corner, if pinning in top
    // left we need to reverse all the sequences first and translate the alignment back later
    
    // create reversed graph if necessary
    Graph reversed_graph;
    if (pin_left) {
        reverse_graph(g, reversed_graph);
    }

    // choose forward or reversed objects
    // note: have to make a copy of the sequence because we will modify it to add a pinning point
    Graph* align_graph = &g;
    string align_sequence = alignment.sequence();
    if (pin_left) {
        align_graph = &reversed_graph;
        reverse(align_sequence.begin(), align_sequence.end());
    }
    
    // convert into gssw graph
    gssw_graph* graph = create_gssw_graph(*align_graph);
    
    // perform dynamic programming
    gssw_graph_fill_pinned(graph, align_sequence.c_str(),
                           nt_table, score_matrix,
                           gap_open, gap_extension, full_length_bonus,
                           pinned ? 0 : full_length_bonus, 15, 2);
    
    // traceback either from pinned position or optimal local alignment
    if (pinned) {
        // trace back pinned alignment
        gssw_graph_mapping** gms = gssw_graph_trace_back_pinned_multi (graph,
                                                                       max_alt_alns,
                                                                       true,
                                                                       align_sequence.c_str(),
                                                                       align_sequence.size(),
                                                                       nt_table,
                                                                       score_matrix,
                                                                       gap_open,
                                                                       gap_extension,
                                                                       full_length_bonus,
                                                                       0);
        
        if (pin_left) {
            // translate graph and mappings into original node space
            unreverse_graph(reversed_graph);
            for (int32_t i = 0; i < max_alt_alns; i++) {
                unreverse_graph_mapping(gms[i]);
            }
        }
        
        // convert optimal alignment and store it in the input Alignment object (in the multi alignment,
        // this will have been set to the first in the vector)
        if (gms[0]->score > 0) {
            // have a mapping, can just convert normally
            gssw_mapping_to_alignment(graph, gms[0], alignment, pinned, pin_left, print_score_matrices);
        }
        else if (g.node_size() > 0) {
            // gssw will not identify mappings with 0 score, infer location based on pinning
            
            Mapping* mapping = alignment.mutable_path()->add_mapping();
            mapping->set_rank(1);
            

            // locate at a beginning of an arbitrary source node or end of an arbitrary sink node as appropriate
            Position* position = mapping->mutable_position();
            if (pin_left) {
                position->set_node_id(g.node(0).id());
                position->set_offset(0);
            }
            else {
                position->set_node_id(g.node(g.node_size() - 1).id());
                position->set_offset(g.node(g.node_size() - 1).sequence().length());
            }
            
            // soft clip
            Edit* edit = mapping->add_edit();
            edit->set_to_length(alignment.sequence().length());
            edit->set_sequence(alignment.sequence());
        }
        
        if (multi_alignments) {
            // determine how many non-null alignments were returned
            int32_t num_non_null = max_alt_alns;
            for (int32_t i = 1; i < max_alt_alns; i++) {
                if (gms[i]->score <= 0) {
                    num_non_null = i;
                    break;
                }
            }
            
            // reserve to avoid illegal access errors that occur when the vector reallocates
            multi_alignments->reserve(num_non_null);
            
            // copy the primary alignment
            multi_alignments->emplace_back(alignment);
            
            // convert the alternate alignments and store them at the back of the vector (this will not
            // execute if we are doing single alignment)
            for (int32_t i = 1; i < num_non_null; i++) {
                gssw_graph_mapping* gm = gms[i];
                
                // make new alignment object
                multi_alignments->emplace_back();
                Alignment& next_alignment = multi_alignments->back();
                
                // copy over sequence information from the primary alignment
                next_alignment.set_sequence(alignment.sequence());
                next_alignment.set_quality(alignment.quality());
                
                // get path of the alternate alignment
                gssw_mapping_to_alignment(graph, gm, next_alignment, pinned, pin_left, print_score_matrices);
                
            }
            
        }
        
        for (int32_t i = 0; i < max_alt_alns; i++) {
            gssw_graph_mapping_destroy(gms[i]);
        }
        free(gms);
    }
    else {
        // trace back local alignment
        gssw_graph_mapping* gm = gssw_graph_trace_back (graph,
                                                        align_sequence.c_str(),
                                                        align_sequence.size(),
                                                        nt_table,
                                                        score_matrix,
                                                        gap_open,
                                                        gap_extension,
                                                        full_length_bonus,
                                                        full_length_bonus);
        
        gssw_mapping_to_alignment(graph, gm, alignment, pinned, pin_left, print_score_matrices);
        gssw_graph_mapping_destroy(gm);
    }
    
    //gssw_graph_print_score_matrices(graph, sequence.c_str(), sequence.size(), stderr);
    
    gssw_graph_destroy(graph);
}

void Aligner::align(Alignment& alignment, Graph& g, bool print_score_matrices) {
    
    align_internal(alignment, nullptr, g, false, false, 1, print_score_matrices);
}

void Aligner::align_pinned(Alignment& alignment, Graph& g, bool pin_left) {
    
    align_internal(alignment, nullptr, g, true, pin_left, 1, false);
}

void Aligner::align_pinned_multi(Alignment& alignment, vector<Alignment>& alt_alignments, Graph& g,
                                 bool pin_left, int32_t max_alt_alns) {
    
    if (alt_alignments.size() != 0) {
        cerr << "error:[Aligner::align_pinned_multi] output vector must be empty for pinned multi-aligning" << endl;
        exit(EXIT_FAILURE);
    }
    
    align_internal(alignment, &alt_alignments, g, true, pin_left, max_alt_alns, false);
}

void Aligner::align_global_banded(Alignment& alignment, Graph& g,
                                  int32_t band_padding, bool permissive_banding) {
    
    // We need to figure out what size ints we need to use.
    // Get upper and lower bounds on the scores. TODO: if these overflow int64 we're out of luck
    int64_t best_score = alignment.sequence().size() * match;
    size_t total_bases = 0;
    for(size_t i = 0; i < g.node_size(); i++) {
        total_bases += g.node(i).sequence().size();
    }
    int64_t worst_score = max(alignment.sequence().size(), total_bases) * -max(max(mismatch, gap_open), gap_extension);
    
    // TODO: put this all into another template somehow?
    
    if (best_score <= numeric_limits<int8_t>::max() && worst_score >= numeric_limits<int8_t>::min()) {
        // We'll fit in int8
        BandedGlobalAligner<int8_t> band_graph(alignment,
                                               g,
                                               band_padding,
                                               permissive_banding,
                                               false);
        
        band_graph.align(score_matrix, nt_table, gap_open, gap_extension);
    } else if (best_score <= numeric_limits<int16_t>::max() && worst_score >= numeric_limits<int16_t>::min()) {
        // We'll fit in int16
        BandedGlobalAligner<int16_t> band_graph(alignment,
                                                g,
                                                band_padding,
                                                permissive_banding,
                                                false);
        
        band_graph.align(score_matrix, nt_table, gap_open, gap_extension);
    } else if (best_score <= numeric_limits<int32_t>::max() && worst_score >= numeric_limits<int32_t>::min()) {
        // We'll fit in int32
        BandedGlobalAligner<int32_t> band_graph(alignment,
                                                g,
                                                band_padding,
                                                permissive_banding,
                                                false);
        
        band_graph.align(score_matrix, nt_table, gap_open, gap_extension);
    } else {
        // Fall back to int64
        BandedGlobalAligner<int64_t> band_graph(alignment,
                                                g,
                                                band_padding,
                                                permissive_banding,
                                                false);
        
        band_graph.align(score_matrix, nt_table, gap_open, gap_extension);
    }

}

void Aligner::align_global_banded_multi(Alignment& alignment, vector<Alignment>& alt_alignments, Graph& g,
                                        int32_t max_alt_alns, int32_t band_padding, bool permissive_banding) {
                                        
    // We need to figure out what size ints we need to use.
    // Get upper and lower bounds on the scores. TODO: if these overflow int64 we're out of luck
    int64_t best_score = alignment.sequence().size() * match;
    size_t total_bases = 0;
    for(size_t i = 0; i < g.node_size(); i++) {
        total_bases += g.node(i).sequence().size();
    }
    int64_t worst_score = max(alignment.sequence().size(), total_bases) * -max(max(mismatch, gap_open), gap_extension);
    
    if (best_score <= numeric_limits<int8_t>::max() && worst_score >= numeric_limits<int8_t>::min()) {
        // We'll fit in int8
        BandedGlobalAligner<int8_t> band_graph(alignment,
                                               g,
                                               alt_alignments,
                                               max_alt_alns,
                                               band_padding,
                                               permissive_banding,
                                               false);
        
        band_graph.align(score_matrix, nt_table, gap_open, gap_extension);
    } else if (best_score <= numeric_limits<int16_t>::max() && worst_score >= numeric_limits<int16_t>::min()) {
        // We'll fit in int16
        BandedGlobalAligner<int16_t> band_graph(alignment,
                                                g,
                                                alt_alignments,
                                                max_alt_alns,
                                                band_padding,
                                                permissive_banding,
                                                false);
        
        band_graph.align(score_matrix, nt_table, gap_open, gap_extension);
    } else if (best_score <= numeric_limits<int32_t>::max() && worst_score >= numeric_limits<int32_t>::min()) {
        // We'll fit in int32
        BandedGlobalAligner<int32_t> band_graph(alignment,
                                                g,
                                                alt_alignments,
                                                max_alt_alns,
                                                band_padding,
                                                permissive_banding,
                                                false);
        
        band_graph.align(score_matrix, nt_table, gap_open, gap_extension);
    } else {
        // Fall back to int64
        BandedGlobalAligner<int64_t> band_graph(alignment,
                                                g,
                                                alt_alignments,
                                                max_alt_alns,
                                                band_padding,
                                                permissive_banding,
                                                false);
        
        band_graph.align(score_matrix, nt_table, gap_open, gap_extension);
    }
}

// Scoring an exact match is very simple in an ordinary Aligner

int32_t Aligner::score_exact_match(const Alignment& aln, size_t read_offset, size_t length) {
    return match * length;
}

int32_t Aligner::score_exact_match(const string& sequence) const {
    return match * sequence.length();
}

int32_t Aligner::score_exact_match(string::const_iterator seq_begin, string::const_iterator seq_end) const {
    return match * (seq_end - seq_begin);
}

QualAdjAligner::QualAdjAligner(int8_t _match,
                               int8_t _mismatch,
                               int8_t _gap_open,
                               int8_t _gap_extension,
                               int8_t _full_length_bonus,
                               int8_t _max_scaled_score,
                               uint8_t _max_qual_score,
                               double gc_content)
{
    
    max_qual_score = _max_qual_score;
    match = _match;
    mismatch = _mismatch;
    gap_open = _gap_open;
    gap_extension = _gap_extension;
    full_length_bonus = _full_length_bonus;
    
    int8_t original_gap_open = gap_open;
    
    nt_table = gssw_create_nt_table();
    score_matrix = gssw_dna_scaled_adjusted_qual_matrix(_max_scaled_score, max_qual_score, &gap_open,
                                                        &gap_extension, match, mismatch,
                                                        gc_content, 1e-12);
    scale_factor = gap_open / original_gap_open;
    match *= scale_factor;
    mismatch *= scale_factor;
    full_length_bonus *= scale_factor;
    
    BaseAligner::init_mapping_quality(gc_content);
}

void QualAdjAligner::align_internal(Alignment& alignment, vector<Alignment>* multi_alignments, Graph& g,
                                    bool pinned, bool pin_left, int32_t max_alt_alns, bool print_score_matrices) {
    
    // check input integrity
    if (pin_left && !pinned) {
        cerr << "error:[Aligner] cannot choose pinned end in non-pinned alignment" << endl;
        exit(EXIT_FAILURE);
    }
    if (multi_alignments && !pinned) {
        cerr << "error:[Aligner] multiple traceback is not implemented in local alignment, only pinned and global" << endl;
        exit(EXIT_FAILURE);
    }
    if (!multi_alignments && max_alt_alns != 1) {
        cerr << "error:[Aligner] cannot specify maximum number of alignments in single alignment" << endl;
        exit(EXIT_FAILURE);
    }
    
    // alignment pinning algorithm is based on pinning in bottom right corner, if pinning in top
    // left we need to reverse all the sequences first and translate the alignment back later
    
    // create reversed graph if necessary
    Graph reversed_graph;
    if (pin_left) {
        reverse_graph(g, reversed_graph);
    }
    
    // choose forward or reversed objects
    // note: have to make copies of the strings because we will modify them to add a pinning point
    Graph* align_graph = &g;
    string align_sequence = alignment.sequence();
    string align_quality = alignment.quality();
    if (pin_left) {
        align_graph = &reversed_graph;
        reverse(align_sequence.begin(), align_sequence.end());
        reverse(align_quality.begin(), align_quality.end());
    }
    
    if (align_quality.length() != align_sequence.length()) {
        cerr << "error:[QualAdjAligner] Read " << alignment.name() << " has sequence and quality strings with different lengths. Cannot perform base quality adjusted alignment. Consider toggling off base quality adjusted alignment." << endl;
        exit(EXIT_FAILURE);
    }
    
    // convert into gssw graph and get dummy pinned node (if pinning)
    gssw_graph* graph = create_gssw_graph(*align_graph);
    
    // perform dynamic programming
    // offer a full length bonus on each end, or only on the left if the right end is pinned.
    gssw_graph_fill_pinned_qual_adj(graph, align_sequence.c_str(), align_quality.c_str(),
                                    nt_table, score_matrix,
                                    gap_open, gap_extension,
                                    full_length_bonus, pinned ? 0 : full_length_bonus, 15, 2);
    
    // traceback either from pinned position or optimal local alignment
    if (pinned) {
        // trace back pinned alignment
        gssw_graph_mapping** gms = gssw_graph_trace_back_pinned_qual_adj_multi (graph,
                                                                                max_alt_alns,
                                                                                true,
                                                                                align_sequence.c_str(),
                                                                                align_quality.c_str(),
                                                                                align_sequence.size(),
                                                                                nt_table,
                                                                                score_matrix,
                                                                                gap_open,
                                                                                gap_extension,
                                                                                full_length_bonus,
                                                                                0);
        
        if (pin_left) {
            // translate graph and mappings into original node space
            unreverse_graph(reversed_graph);
            for (int32_t i = 0; i < max_alt_alns; i++) {
                unreverse_graph_mapping(gms[i]);
            }
        }
        
        // convert optimal alignment and store it in the input Alignment object (in the multi alignment,
        // this will have been set to the first in the vector)
        if (gms[0]->score > 0) {
            // have a mapping, can just convert normally
            gssw_mapping_to_alignment(graph, gms[0], alignment, pinned, pin_left, print_score_matrices);
        }
        else if (g.node_size() > 0) {
            // gssw will not identify mappings with 0 score, infer location based on pinning
            
            Mapping* mapping = alignment.mutable_path()->add_mapping();
            mapping->set_rank(1);
            
            // locate at a beginning of a source node or end of a sink node as appropriate
            Position* position = mapping->mutable_position();
            if (pin_left) {
                position->set_node_id(g.node(0).id());
                position->set_offset(0);
            }
            else {
                position->set_node_id(g.node(g.node_size() - 1).id());
                position->set_offset(g.node(g.node_size() - 1).sequence().length());
            }
            
            // soft clip
            Edit* edit = mapping->add_edit();
            edit->set_to_length(alignment.sequence().length());
            edit->set_sequence(alignment.sequence());
        }
        
        
        if (multi_alignments) {
            // determine how many non-null alignments were returned
            int32_t num_non_null = max_alt_alns;
            for (int32_t i = 1; i < max_alt_alns; i++) {
                if (gms[i]->score <= 0) {
                    num_non_null = i;
                    break;
                }
            }
            
            // reserve to avoid illegal access errors that occur when the vector reallocates
            multi_alignments->reserve(num_non_null);
            
            // copy the primary alignment
            multi_alignments->emplace_back(alignment);
            
            // convert the alternate alignments and store them at the back of the vector (this will not
            // execute if we are doing single alignment)
            for (int32_t i = 1; i < num_non_null; i++) {
                gssw_graph_mapping* gm = gms[i];
                
                // make new alignment object
                multi_alignments->emplace_back();
                Alignment& next_alignment = multi_alignments->back();
                
                // copy over sequence information from the primary alignment
                next_alignment.set_sequence(alignment.sequence());
                next_alignment.set_quality(alignment.quality());
                
                // get path of the alternate alignment
                gssw_mapping_to_alignment(graph, gm, next_alignment, pinned, pin_left, print_score_matrices);
                
            }
        }
        
        for (int32_t i = 0; i < max_alt_alns; i++) {
            gssw_graph_mapping_destroy(gms[i]);
        }
        free(gms);
    }
    else {
        // trace back local alignment
        gssw_graph_mapping* gm = gssw_graph_trace_back_qual_adj (graph,
                                                                 align_sequence.c_str(),
                                                                 align_quality.c_str(),
                                                                 align_sequence.size(),
                                                                 nt_table,
                                                                 score_matrix,
                                                                 gap_open,
                                                                 gap_extension,
                                                                 full_length_bonus,
                                                                 full_length_bonus);
        
        gssw_mapping_to_alignment(graph, gm, alignment, pinned, pin_left, print_score_matrices);
        gssw_graph_mapping_destroy(gm);
    }
    
    //gssw_graph_print_score_matrices(graph, sequence.c_str(), sequence.size(), stderr);
    
    gssw_graph_destroy(graph);
    
}

void QualAdjAligner::align(Alignment& alignment, Graph& g, bool print_score_matrices) {
    
    align_internal(alignment, nullptr, g, false, false, 1, print_score_matrices);
}

void QualAdjAligner::align_pinned(Alignment& alignment, Graph& g, bool pin_left) {

    align_internal(alignment, nullptr, g, true, pin_left, 1, false);

}

void QualAdjAligner::align_pinned_multi(Alignment& alignment, vector<Alignment>& alt_alignments, Graph& g,
                                        bool pin_left, int32_t max_alt_alns) {
    align_internal(alignment, &alt_alignments, g, true, pin_left, max_alt_alns, false);
}

void QualAdjAligner::align_global_banded(Alignment& alignment, Graph& g,
                                         int32_t band_padding, bool permissive_banding) {
    
    BandedGlobalAligner<int16_t> band_graph = BandedGlobalAligner<int16_t>(alignment,
                                                                           g,
                                                                           band_padding,
                                                                           permissive_banding,
                                                                           true);
    
    band_graph.align(score_matrix, nt_table, gap_open, gap_extension);
}

void QualAdjAligner::align_global_banded_multi(Alignment& alignment, vector<Alignment>& alt_alignments, Graph& g,
                                               int32_t max_alt_alns, int32_t band_padding, bool permissive_banding) {
    
    BandedGlobalAligner<int16_t> band_graph = BandedGlobalAligner<int16_t>(alignment,
                                                                           g,
                                                                           alt_alignments,
                                                                           max_alt_alns,
                                                                           band_padding,
                                                                           permissive_banding,
                                                                           true);
    
    band_graph.align(score_matrix, nt_table, gap_open, gap_extension);
}

int32_t QualAdjAligner::score_exact_match(const Alignment& aln, size_t read_offset, size_t length) {
    auto& sequence = aln.sequence();
    auto& base_quality = aln.quality();
    int32_t score = 0;
    for (int32_t i = 0; i < length; i++) {
        // index 5 x 5 score matrices (ACGTN)
        // always have match so that row and column index are same and can combine algebraically
        score += score_matrix[25 * base_quality[read_offset + i] + 6 * nt_table[sequence[read_offset + i]]];
    }
    return score;
}

int32_t QualAdjAligner::score_exact_match(const string& sequence, const string& base_quality) const {
    int32_t score = 0;
    for (int32_t i = 0; i < sequence.length(); i++) {
        // index 5 x 5 score matrices (ACGTN)
        // always have match so that row and column index are same and can combine algebraically
        score += score_matrix[25 * base_quality[i] + 6 * nt_table[sequence[i]]];
    }
    return score;
}


int32_t QualAdjAligner::score_exact_match(string::const_iterator seq_begin, string::const_iterator seq_end,
                                          string::const_iterator base_qual_begin) const {
    int32_t score = 0;
    for (auto seq_iter = seq_begin, qual_iter = base_qual_begin; seq_iter != seq_end; seq_iter++) {
        // index 5 x 5 score matrices (ACGTN)
        // always have match so that row and column index are same and can combine algebraically
        score += score_matrix[25 * (*qual_iter) + 6 * nt_table[*seq_iter]];
        qual_iter++;
    }
    return score;
}
