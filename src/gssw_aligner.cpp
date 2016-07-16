#include <cmath>
#include "gssw_aligner.hpp"

// log(10)
static const double quality_scale_factor = 10.0 / log(10.0);
static const double exp_overflow_limit = log(std::numeric_limits<double>::max());

using namespace vg;

double add_log(double log_x, double log_y) {
    if (log_x > log_y) {
        return log_x + log(1.0 + exp(log_y - log_x));
    }
    else {
        return log_y + log(1.0 + exp(log_x - log_y));
    }
}

Aligner::~Aligner(void) {
    free(nt_table);
    free(score_matrix);
}

Aligner::Aligner(int32_t _match,
                 int32_t _mismatch,
                 int32_t _gap_open,
                 int32_t _gap_extension)
{
    match = _match;
    mismatch = _mismatch;
    gap_open = _gap_open;
    gap_extension = _gap_extension;
    
    // these are used when setting up the nodes
    nt_table = gssw_create_nt_table();
    score_matrix = gssw_create_score_matrix(match, mismatch);
    log_base = 0.0;
}


gssw_graph* Aligner::create_gssw_graph(Graph& g) {
    
    gssw_graph* graph = gssw_graph_create(g.node_size());

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

void Aligner::align(Alignment& alignment, Graph& g, bool print_score_matrices) {
    
    gssw_graph* graph = create_gssw_graph(g);
    const string& sequence = alignment.sequence();
    

    gssw_graph_fill(graph, sequence.c_str(),
                    nt_table, score_matrix,
                    gap_open, gap_extension, 15, 2);

    gssw_graph_mapping* gm = gssw_graph_trace_back (graph,
                                                    sequence.c_str(),
                                                    sequence.size(),
                                                    nt_table,
                                                    score_matrix,
                                                    gap_open,
                                                    gap_extension);

    //gssw_graph_print_score_matrices(graph, sequence.c_str(), sequence.size(), stderr);

    gssw_mapping_to_alignment(graph, gm, alignment, print_score_matrices);

#ifdef debug
    gssw_print_graph_mapping(gm, stderr);
#endif
    
    gssw_graph_mapping_destroy(gm);
    gssw_graph_destroy(graph);

}

void Aligner::gssw_mapping_to_alignment(gssw_graph* graph,
                                        gssw_graph_mapping* gm,
                                        Alignment& alignment,
                                        bool print_score_matrices) {

    alignment.clear_path();
    alignment.set_score(gm->score);
    alignment.set_query_position(0);
    Path* path = alignment.mutable_path();
    //alignment.set_cigar(graph_cigar(gm));

    gssw_graph_cigar* gc = &gm->cigar;
    gssw_node_cigar* nc = gc->elements;
    int to_pos = 0;
    int from_pos = gm->position;
    //cerr << "gm->position " << gm->position << endl;
    string& to_seq = *alignment.mutable_sequence();
    //cerr << "-------------" << endl;

    if (print_score_matrices) {
        gssw_graph_print_score_matrices(graph, to_seq.c_str(), to_seq.size(), stderr);
        //cerr << alignment.DebugString() << endl;
    }

    for (int i = 0; i < gc->length; ++i, ++nc) {
        if (i > 0) from_pos = 0; // reset for each node after the first
        // check that the current alignment has a non-zero length
        gssw_cigar* c = nc->cigar;
        int l = c->length;
        if (l == 0) continue;
        gssw_cigar_element* e = c->elements;

        Node* from_node = (Node*) nc->node->data;
        string& from_seq = *from_node->mutable_sequence();
        Mapping* mapping = path->add_mapping();
        mapping->mutable_position()->set_node_id(nc->node->id);
        mapping->mutable_position()->set_offset(from_pos);
        mapping->set_rank(path->mapping_size());

        //cerr << from_node->id() << ":" << endl;

        for (int j=0; j < l; ++j, ++e) {
            Edit* edit;
            int32_t length = e->length;
            //cerr << e->length << e->type << endl;
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
                cerr << "error [Aligner::gssw_mapping_to_alignment] "
                     << "unsupported cigar op type " << e->type << endl;
                exit(1);
                break;

            }

        }
        //cerr << "path to_length " << path_to_length(*path) << endl;
    }

    // set identity
    alignment.set_identity(identity(alignment.path()));

}

string Aligner::graph_cigar(gssw_graph_mapping* gm) {
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

void Aligner::init_mapping_quality(double gc_content) {
    log_base = gssw_dna_recover_log_base(match, mismatch, gc_content, 1e-12);
}

bool Aligner::is_mapping_quality_initialized() {
    return (log_base <= 0.0);
}

double Aligner::maximum_mapping_quality_exact(vector<double>& scaled_scores, size_t* max_idx_out) {
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
    
    if (max_score * size < exp_overflow_limit) {
        // no risk of double overflow, sum exp directly (half as many transcendental function evals)
        double numer = 0.0;
        for (size_t i = 0; i < size; i++) {
            if (i == max_idx) {
                continue;
            }
            numer += exp(scaled_scores[i]);
        }
        return -10.0 * log10(numer / (numer + exp(scaled_scores[max_idx])));
    }
    else {
        // work in log transformed valued to avoid risk of overflow
        double log_sum_exp = scaled_scores[0];
        for (size_t i = 1; i < size; i++) {
            log_sum_exp = add_log(log_sum_exp, scaled_scores[i]);
        }
        return -10.0 * log10(1.0 - exp(scaled_scores[max_idx] - log_sum_exp));
    }
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

double Aligner::maximum_mapping_quality_approx(vector<double>& scaled_scores, size_t* max_idx_out) {
    size_t size = scaled_scores.size();
    
    // if necessary, assume a null alignment of 0.0 for comparison since this is local
    if (size == 1) {
        scaled_scores.push_back(0.0);
    }

    double max_score = scaled_scores[0];
    size_t max_idx = 0;
    
    double next_score = -std::numeric_limits<double>::max();
    int32_t next_count = 0;
    
    for (int32_t i = 1; i < size; i++) {
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
    
    return quality_scale_factor * (max_score - next_count * next_score);
}

void Aligner::compute_mapping_quality(vector<Alignment>& alignments, bool fast_approximation) {
    if (log_base <= 0.0) {
        cerr << "error:[Aligner] must call init_mapping_quality before computing mapping qualities" << endl;
        exit(EXIT_FAILURE);
    }
    
    size_t size = alignments.size();
    
    if (size == 0) {
        return;
    }
    
    vector<double> scaled_scores(size);
    
    for (size_t i = 0; i < size; i++) {
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
    if (mapping_quality > std::numeric_limits<int32_t>::max()) {
        alignments[max_idx].set_mapping_quality(std::numeric_limits<int32_t>::max());
    }
    else {
        alignments[max_idx].set_mapping_quality((int32_t) mapping_quality);
    }
}

void Aligner::compute_paired_mapping_quality(pair<vector<Alignment>, vector<Alignment>>& alignment_pairs,
                                             bool fast_approximation) {
    if (log_base <= 0.0) {
        cerr << "error:[Aligner] must call init_mapping_quality before computing mapping qualities" << endl;
        exit(EXIT_FAILURE);
    }
    
    size_t size = alignment_pairs.first.size();
    
    if (size != alignment_pairs.second.size()) {
        cerr << "error:[Aligner] unpaired alignments included with pairs" << endl;
        exit(EXIT_FAILURE);
    }
    
    if (size == 0) {
        return;
    }
    
    vector<double> scaled_scores(size);
    
    for (size_t i = 0; i < size; i++) {
        scaled_scores[i] = log_base * (alignment_pairs.first[i].score() + alignment_pairs.second[i].score());
    }
    
    size_t max_idx;
    double mapping_quality;
    if (!fast_approximation) {
        mapping_quality = maximum_mapping_quality_exact(scaled_scores, &max_idx);
    }
    else {
        mapping_quality = maximum_mapping_quality_approx(scaled_scores, &max_idx);
    }
    
    if (mapping_quality > std::numeric_limits<int32_t>::max()) {
        alignment_pairs.first[max_idx].set_mapping_quality(std::numeric_limits<int32_t>::max());
        alignment_pairs.second[max_idx].set_mapping_quality(std::numeric_limits<int32_t>::max());
    }
    else {
        alignment_pairs.first[max_idx].set_mapping_quality((int32_t) mapping_quality);
        alignment_pairs.second[max_idx].set_mapping_quality((int32_t) mapping_quality);
    }

}

QualAdjAligner::QualAdjAligner(int8_t _match,
                               int8_t _mismatch,
                               int8_t _gap_open,
                               int8_t _gap_extension,
                               int8_t _max_scaled_score,
                               uint8_t _max_qual_score,
                               double gc_content) : Aligner(_match, _mismatch, _gap_open, _gap_extension) {
    
    init_quality_adjusted_scores(_max_scaled_score, _max_qual_score, gc_content);
}

void QualAdjAligner::init_quality_adjusted_scores(int8_t _max_scaled_score,
                                                  uint8_t _max_qual_score,
                                                  double gc_content) {
    max_qual_score = _max_qual_score;
    scaled_gap_open = gap_open;
    scaled_gap_extension = gap_extension;
    
    adjusted_score_matrix = gssw_dna_scaled_adjusted_qual_matrix(_max_scaled_score, max_qual_score, &scaled_gap_open,
                                                                &scaled_gap_extension, match, mismatch,
                                                                gc_content, 1e-12);
    init_mapping_quality(gc_content);
}

void QualAdjAligner::init_mapping_quality(double gc_content) {
    log_base = gssw_dna_recover_log_base(match, mismatch, gc_content, 1e-12);
    // adjust to scaled matrix (a bit hacky but correct)
    log_base /= (scaled_gap_open / gap_open);
}

QualAdjAligner::~QualAdjAligner(void) {
    free(adjusted_score_matrix);
}

void QualAdjAligner::align(Alignment& alignment, Graph& g, bool print_score_matrices) {
    
    gssw_graph* graph = create_gssw_graph(g);
    
    const string& sequence = alignment.sequence();
    const string& quality = alignment.quality();
    
    if (quality.length() != sequence.length()) {
        cerr << "error:[Aligner] sequence and quality strings different lengths, cannot perform base quality adjusted alignment" << endl;
    }

    gssw_graph_fill_qual_adj(graph, sequence.c_str(), quality.c_str(),
                             nt_table, adjusted_score_matrix,
                             scaled_gap_open, scaled_gap_extension, 15, 2);

    gssw_graph_mapping* gm = gssw_graph_trace_back_qual_adj (graph,
                                                             sequence.c_str(),
                                                             quality.c_str(),
                                                             sequence.size(),
                                                             nt_table,
                                                             adjusted_score_matrix,
                                                             gap_open,
                                                             gap_extension);

    gssw_mapping_to_alignment(graph, gm, alignment, print_score_matrices);

#ifdef debug
    gssw_print_graph_mapping(gm, stderr);
#endif

    gssw_graph_mapping_destroy(gm);
    gssw_graph_destroy(graph);
}
