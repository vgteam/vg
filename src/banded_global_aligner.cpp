//
//  banded_global_aligner.cpp
//  
//  Contains objects to support Aligner in performing banded global alignment against
//  a graph.
//

#include "banded_global_aligner.hpp"
#include "vg/io/json2pb.h"

#include <array>

//#define debug_banded_aligner_objects
//#define debug_banded_aligner_graph_processing
//#define debug_banded_aligner_fill_matrix
//#define debug_banded_aligner_traceback
//#define debug_banded_aligner_print_matrices
//#define debug_jemalloc

#ifdef debug_jemalloc
#include <jemalloc/jemalloc.h>
#endif

namespace vg {

template<class IntType>
BandedGlobalAligner<IntType>::BABuilder::BABuilder(Alignment& alignment) :
                                                   alignment(alignment),
                                                   matrix_state(Match),
                                                   matching(false),
                                                   current_node_id(0),
                                                   edit_length(0),
                                                   edit_read_end_idx(0)
{
    // nothing to do
}

template<class IntType>
BandedGlobalAligner<IntType>::BABuilder::~BABuilder() {
    // does not own any heap objects, nothing to do
}

template<class IntType>
void BandedGlobalAligner<IntType>::BABuilder::update_state(const HandleGraph& graph, matrix_t matrix,
                                                           const handle_t& node, int64_t read_idx, int64_t node_idx,
                                                           bool empty_node_seq) {
#ifdef debug_banded_aligner_traceback
    cerr << "[BABuilder::update_state] beginning " << (empty_node_seq ? "" : "non-") << "empty state update for read index " << read_idx << ", node seq index " << node_idx << endl;
#endif
    if (graph.get_id(node) != current_node_id) {
#ifdef debug_banded_aligner_traceback
        cerr << "[BABuilder::update_state] at new node " << graph.get_id(node) << " previously " << current_node_id << endl;
#endif
        // conclude current mapping and proceed to next node
        finish_current_node();
        current_node_id = graph.get_id(node);
        current_node_sequence = graph.get_sequence(node);
        matrix_state = matrix;
        if (matrix_state == Match) {
            matching = (alignment.sequence()[read_idx] == current_node_sequence[node_idx]);
        }
        edit_length = !empty_node_seq;
        edit_read_end_idx = read_idx;
    }
    else if (matrix != matrix_state) {
#ifdef debug_banded_aligner_traceback
        cerr << "[BABuilder::update_state] transitioning into a new matrix" << endl;
#endif
        // transitioned into another matrix, finish current edit and begin new one
        if (!empty_node_seq) {
            finish_current_edit();
        }
        matrix_state = matrix;
        if (matrix_state == Match) {
            matching = (alignment.sequence()[read_idx] == current_node_sequence[node_idx]);
        }
        edit_length = 1;
        edit_read_end_idx = read_idx;
    }
    else if (matrix == Match &&
             (alignment.sequence()[read_idx] == current_node_sequence[node_idx]) != matching) {
#ifdef debug_banded_aligner_traceback
        cerr << "[BABuilder::update_state] switching between match and mismatch" << endl;
#endif
        // switch from match to mismatch state or vice versa
        finish_current_edit();
        matching = !matching;
        edit_length = 1;
        edit_read_end_idx = read_idx;
    }
    else {
        // same edit, extend length
        edit_length++;
    }
    
#ifdef debug_banded_aligner_traceback
    cerr << "[BABuilder::update_state] finished updating state, matrix is " << (matrix_state == Match ? "match" : (matrix_state == InsertRow ? "insert row" : "insert column" )) << ", is matching? " << (matching ? "yes" : "no") << ", edit length " << edit_length << ", edit end index (on read) " << edit_read_end_idx << ", current node " << current_node_id << endl;
#endif
}

template <class IntType>
void BandedGlobalAligner<IntType>::BABuilder::finish_current_edit() {
    
#ifdef debug_banded_aligner_traceback
    cerr << "[BABuilder::finish_current_edit] finishing edit" << endl;
#endif
    
    mapping_edits.emplace_front();
    
    switch (matrix_state) {
        case Match:
            mapping_edits.front().set_from_length(edit_length);
            mapping_edits.front().set_to_length(edit_length);
            
            if (!matching) {
                mapping_edits.front().set_sequence(alignment.sequence().substr(edit_read_end_idx - edit_length + 1,
                                                                               edit_length));
            }
            
            break;
            
        case InsertRow:
            mapping_edits.front().set_from_length(0);
            mapping_edits.front().set_to_length(edit_length);
            mapping_edits.front().set_sequence(alignment.sequence().substr(edit_read_end_idx - edit_length + 1,
                                                                           edit_length));
            break;
            
        case InsertCol:
            mapping_edits.front().set_from_length(edit_length);
            mapping_edits.front().set_to_length(0);
            break;
            
        default:
            cerr << "error:[BandedGlobalAligner] unrecognized matrix type" << endl;
            assert(0);
            break;
    }
#ifdef debug_banded_aligner_traceback
    cerr << "[BABuilder::finish_current_edit] edit: " << pb2json(mapping_edits.front()) << endl;
#endif
}

template <class IntType>
void BandedGlobalAligner<IntType>::BABuilder::finish_current_node() {

    // sentinel for first iteration
    if (current_node_id == 0) {
#ifdef debug_banded_aligner_traceback
        cerr << "[BABuilder::finish_current_node] at beginning of traceback, not creating a mapping" << endl;
#endif
        return;
    }
    
    finish_current_edit();
    
#ifdef debug_banded_aligner_traceback
    cerr << "[BABuilder::finish_current_node] finishing mapping for node " << current_node_id << endl;
#endif
    
    node_mappings.emplace_front();
    for (Edit edit : mapping_edits) {
        *(node_mappings.front().add_edit()) = edit;
    }
    
#ifdef debug_banded_aligner_traceback
    cerr << "[BABuilder::finish_current_node] mapping: " << pb2json(node_mappings.front()) << endl;
#endif
    
    mapping_edits.clear();
    
    (*(node_mappings.front().mutable_position())).set_node_id(current_node_id);
    // note: global alignment always starts at beginning of node, default offset 0 is correct
}

template <class IntType>
void BandedGlobalAligner<IntType>::BABuilder::finalize_alignment(const list<int64_t>& empty_prefix) {
    
    finish_current_node();
    
#ifdef debug_banded_aligner_traceback
    cerr << "[BABuilder::finalize_alignment] finalizing alignment" << endl;
    cerr << "[BABuilder::finalize_alignment] empty prefix is size " << empty_prefix.size() << endl;
#endif
    
    alignment.clear_path();
    Path* path = alignment.mutable_path();
    
    int32_t mapping_rank = 1;
    for (Mapping& mapping : node_mappings) {
        mapping.set_rank(mapping_rank);
        *(path->add_mapping()) = mapping;
        mapping_rank++;
    }
    
    // if there were any empty nodes traversed before the traceback began, add those to the end
    for (int64_t empty_id : empty_prefix) {
        Mapping* mapping = path->add_mapping();
        mapping->add_edit();
        mapping->mutable_position()->set_node_id(empty_id);
        mapping->set_rank(mapping_rank);
        mapping_rank++;
    }
    
    node_mappings.clear();
    
#ifdef debug_banded_aligner_traceback
    cerr << "[BABuilder::finalize_alignment] alignment: " << pb2json(alignment) << endl;
#endif
}

template <class IntType>
BandedGlobalAligner<IntType>::BAMatrix::BAMatrix(Alignment& alignment, handle_t node, int64_t top_diag, int64_t bottom_diag,
                                                 const vector<BAMatrix*>& seeds, int64_t cumulative_seq_len) :
                                                 node(node),
                                                 top_diag(top_diag),
                                                 bottom_diag(bottom_diag),
                                                 seeds(seeds),
                                                 alignment(alignment),
                                                 cumulative_seq_len(cumulative_seq_len),
                                                 match(nullptr),
                                                 insert_col(nullptr),
                                                 insert_row(nullptr)
{
    // nothing to do
#ifdef debug_banded_aligner_objects
    cerr << "[BAMatrix]: constructor for node " << as_integer(node) << " and band from " << top_diag << " to " << bottom_diag << endl;;
#endif
}

template <class IntType>
BandedGlobalAligner<IntType>::BAMatrix::~BAMatrix() {
#ifdef debug_banded_aligner_objects
    cerr << "[BAMatrix::~BAMatrix] destructing matrix for handle " << handlegraph::as_integer(node) << endl;
#endif
    free(match);
    free(insert_row);
    free(insert_col);
}

template <class IntType>
void BandedGlobalAligner<IntType>::BAMatrix::fill_matrix(const HandleGraph& graph, int8_t* score_mat, int8_t* nt_table,
                                                         int8_t gap_open, int8_t gap_extend, bool qual_adjusted, IntType min_inf) {
    
#ifdef debug_banded_aligner_fill_matrix
    cerr << "[BAMatrix::fill_matrix] beginning DP on matrix for node " << as_integer(node) << endl;;
#endif
    
    // note: bottom has the higher index
    int64_t band_height = bottom_diag - top_diag + 1;
    int64_t ncols = graph.get_length(node);
    int64_t band_size = band_height * ncols;
    
#ifdef debug_banded_aligner_fill_matrix
    cerr << "[BAMatrix::fill_matrix]: allocating matrices of height " << band_height << " and width " << ncols << " for a total cell count of " << band_size << endl;
#endif
    
    string node_seq = graph.get_sequence(node);
    const string& read = alignment.sequence();
    const string& base_quality = alignment.quality();
    
    match = (IntType*) malloc(sizeof(IntType) * band_size);
    insert_col = (IntType*) malloc(sizeof(IntType) * band_size);
    insert_row = (IntType*) malloc(sizeof(IntType) * band_size);
    if (!match || !insert_col || !insert_row) {
        // An allocation has failed.
        // We may have run out of virtual memory.
        
#ifdef debug_jemalloc
        size_t requested_size = sizeof(IntType) * band_size;
        size_t usable_size[3] = {0, 0, 0};
#endif
        
        // Free up what we are holding, and also report how much usable mamoey jemalloc gave us for anything that passed.
        if (match) {
#ifdef debug_jemalloc
            usable_size[0] = malloc_usable_size(match);
#endif
            free(match);
        }
        if (insert_col) {
#ifdef debug_jemalloc
            usable_size[1] = malloc_usable_size(insert_col);
#endif
            free(insert_col);
        }
        if (insert_row) {
#ifdef debug_jemalloc
            usable_size[2] = malloc_usable_size(insert_row);
#endif
            free(insert_row);
        }
        
        cerr << "[BAMatrix::fill_matrix]: failed to allocate matrices of height " << band_height << " and width " << ncols << " for a total cell count of " << band_size << endl;
#ifdef debug_jemalloc
        cerr << "[BAMatrix::fill_matrix]: requested: " << requested_size << " actually obtained: " << usable_size[0] << " " << usable_size[1] << " " << usable_size[2] << endl;
#endif
        cerr << "[BAMatrix::fill_matrix]: is alignment problem too big for your virtual or physical memory?" << endl;
        
#ifdef debug_jemalloc
        // Dump the stats from the allocator.
        // TODO: skip when not building with jemalloc somehow.
        malloc_stats_print(nullptr, nullptr, "");
#endif
        
        // Bail out relatively safely
        throw std::bad_alloc();
    }
    
    /* these represent a band in a matrix, but we store it as a rectangle with chopped
     * corners
     *
     * 1XX2
     * XXXXX               2XXXXX4
     * XXXXXX             XXXXXXXX
     * 3XXXXXX           XXXXXXXXX
     *  XXXXXXX    -->  1XXXXXXXX6
     *   XXXXXXX        XXXXXXXXX
     *    XXXXXX4       XXXXXXXX
     *     XXXXXX       3XXXXX5
     *      XXXXX
     *       5XX6
     *
     * this way the dimensions are band_height x seq_len instead of read_len x seq_len
     * also note that the internal structure of each column is preserved and each row
     * in the rectangularized band corresponds to a diagonal in the original matrix
     *
     * the initial row and column can be reached via an implied row or column insertion
     * that is not represented in the matrix (this requires a number of edge cases)
     */
    
    if (!band_size) {
        return;
    }
    
    int64_t idx, up_idx, diag_idx, left_idx;
    
    // rows in the rectangularized band corresponding to diagonals
    int64_t iter_start = top_diag < 0 ? -top_diag : 0;
    int64_t iter_stop = bottom_diag >= (int64_t) read.length() ? band_height + (int64_t) read.length() - bottom_diag - 1 : band_height;
    
    // initialize with min infs (identity of max function)
    for (int64_t i = iter_start; i < iter_stop; i++) {
        idx = i * ncols;
        match[idx] = min_inf;
        insert_col[idx] = min_inf;
        // can skip insert row since it doesn't cross node boundaries
    }
    
    // make sure this one insert row value is there so we can use it for checking band boundaries
    // later
    insert_row[iter_start * ncols] = min_inf;
    
    // we will allow the alignment to treat this node as a source if it has no seeds or if it
    // is connected to a source node by a length 0 path (which we will check later)
    bool treat_as_source = seeds.empty();
    
    // initialize the queue with all of the predecessors
    vector<BAMatrix*> seed_queue = seeds;
    
    while (!seed_queue.empty()) {
        BAMatrix* seed = seed_queue.back();
        seed_queue.pop_back();
        
        if (seed == nullptr) {
#ifdef debug_banded_aligner_fill_matrix
            cerr << "[BAMatrix::fill_matrix]: seed is masked, skipping" << endl;
#endif
            continue;
        }
        
#ifdef debug_banded_aligner_fill_matrix
        cerr << "[BAMatrix::fill_matrix]: doing POA across boundary from seed node " << graph.get_id(seed->node) << " to node " << graph.get_id(node) << endl;
#endif
        
        int64_t seed_node_seq_len = graph.get_length(seed->node);
        
        if (seed_node_seq_len == 0) {
#ifdef debug_banded_aligner_fill_matrix
            cerr << "[BAMatrix::fill_matrix]: seed node " << graph.get_id(seed->node) << " has no sequence, adding its predecessors as seed nodes" << endl;
#endif
            // this is a length 0 node, so let this seed's predecessors seed into this one or
            // identify this node as a "source"
            treat_as_source = treat_as_source || seed->seeds.empty();
            for (auto& seed_of_seed : seed->seeds) {
                seed_queue.push_back(seed_of_seed);
            }
            continue;
        }
        
        // compute the interval of diagonals that this seed reaches
        int64_t seed_next_top_diag = seed->top_diag + seed_node_seq_len;
        int64_t seed_next_bottom_diag = seed->bottom_diag + seed_node_seq_len;
        
        // characterize the this interval
        bool beyond_top_of_matrix = seed_next_top_diag < 0;
        bool abutting_top_of_matrix = seed_next_top_diag == 0;
        bool beyond_bottom_of_matrix = seed_next_bottom_diag >= (int64_t) read.length();
        int64_t seed_next_top_diag_iter = beyond_top_of_matrix ? 0 : seed_next_top_diag;
        int64_t seed_next_bottom_diag_iter = beyond_bottom_of_matrix ? (int64_t) read.length() - 1 : seed_next_bottom_diag;
        
        // the shortest sequence path to any source that goes through this seed
        int64_t extended_cumulative_seq_len = seed->cumulative_seq_len + seed_node_seq_len;
        
#ifdef debug_banded_aligner_fill_matrix
        cerr << "[BAMatrix::fill_matrix]: this seed reaches diagonals " << seed_next_top_diag << " to " << seed_next_bottom_diag << " out of matrix range " << top_diag << " to " << bottom_diag << endl;
#endif
        // special logic for first row
        idx = (seed_next_top_diag_iter - top_diag) * ncols;
        
        IntType match_score;
        if (qual_adjusted) {
            match_score = score_mat[25 * base_quality[seed_next_top_diag_iter] + 5 * nt_table[node_seq[0]] + nt_table[read[seed_next_top_diag_iter]]];
        }
        else {
            match_score = score_mat[5 * nt_table[node_seq[0]] + nt_table[read[seed_next_top_diag_iter]]];
        }
        
        if (beyond_top_of_matrix) {
            // the implied cell above this cell is within the extended band form this seed, so we can extend from
            // paths through this node into both the match and insert row from a lead gap
            
            // match after implied gap along top edge
            match[idx] = max<IntType>(match_score - gap_open - (extended_cumulative_seq_len - 1) * gap_extend, match[idx]);
            // gap open after implied gap along top edge
            insert_row[idx] = max<IntType>(-2 * gap_open - extended_cumulative_seq_len * gap_extend, insert_row[idx]);
        }
        else if (abutting_top_of_matrix) {
            // the implied cell above this cell is not in the extended band, but the one diagonal is, so we can extend
            // into match from a lead gap but not insert row
            
            // match after implied gap along top edge
            match[idx] = max<IntType>(match_score - gap_open - (extended_cumulative_seq_len - 1) * gap_extend, match[idx]);
            
        }
        else {
#ifdef debug_banded_aligner_fill_matrix
            cerr << "[BAMatrix::fill_matrix]: top cell in match matrix is reachable without a lead gap" << endl;
#endif
            diag_idx = (seed_next_top_diag_iter - seed_next_top_diag) * seed_node_seq_len + seed_node_seq_len - 1;
            
            match[idx] = max<IntType>(match_score + max<IntType>(max<IntType>(seed->match[diag_idx],
                                                                              seed->insert_row[diag_idx]),
                                                                 seed->insert_col[diag_idx]), match[idx]);
        }
        
        if (seed_next_top_diag < seed_next_bottom_diag) {
#ifdef debug_banded_aligner_fill_matrix
            cerr << "[BAMatrix::fill_matrix]: seed band is greater than height 1, can extend column gap into first row" << endl;
#endif
            left_idx = (seed_next_top_diag_iter - seed_next_top_diag + 2) * seed_node_seq_len - 1;
            insert_col[idx] = max<IntType>(max<IntType>(max<IntType>(seed->match[left_idx] - gap_open,
                                                                     seed->insert_row[left_idx] - gap_open),
                                                        seed->insert_col[left_idx] - gap_extend), insert_col[idx]);
        }
        
        
        for (int64_t diag = seed_next_top_diag_iter + 1; diag < seed_next_bottom_diag_iter; diag++) {
            idx = (diag - top_diag) * ncols;
            
#ifdef debug_banded_aligner_fill_matrix
            cerr << "[BAMatrix::fill_matrix]: extending a match and column gap into matrix coord (" << diag << ", 0)" << ", rectangular coord coord (" << diag - top_diag << ", 0)" << endl;
#endif
            
            // extend a match
            diag_idx = (diag - seed_next_top_diag) * seed_node_seq_len + seed_node_seq_len - 1;
            if (qual_adjusted) {
                match_score = score_mat[25 * base_quality[diag] + 5 * nt_table[node_seq[0]] + nt_table[read[diag]]];
            }
            else {
                match_score = score_mat[5 * nt_table[node_seq[0]] + nt_table[read[diag]]];
            }
            
#ifdef debug_banded_aligner_fill_matrix
            cerr << "[BAMatrix::fill_matrix]: extending match from rectangular coord (" << diag - seed_next_top_diag << ", " << seed_node_seq_len - 1 << ")" << " with match score " << (int) match_score << ", scores are " << (int) seed->match[diag_idx] << " (M), " << (int) seed->insert_row[diag_idx] << " (Ir), and " << (int) seed->insert_col[diag_idx] << " (Ic), current score is " << (int) match[idx] << endl;
#endif
            
            match[idx] = max<IntType>(match_score + max<IntType>(max<IntType>(seed->match[diag_idx],
                                                                              seed->insert_row[diag_idx]),
                                                                 seed->insert_col[diag_idx]), match[idx]);
            
            // extend a column gap
            left_idx = (diag - seed_next_top_diag + 2) * seed_node_seq_len - 1;
            
#ifdef debug_banded_aligner_fill_matrix
            cerr << "[BAMatrix::fill_matrix]: extending match from rectangular coord (" << diag - seed_next_top_diag + 1 << ", " << seed_node_seq_len - 1 << ")" << ", scores are " << (int) seed->match[left_idx] << " (M), " << (int) seed->insert_row[left_idx] << " (Ir), and " << (int) seed->insert_col[left_idx] << " (Ic), current score is " << (int) insert_col[idx] << endl;
#endif
            insert_col[idx] = max<IntType>(max<IntType>(max<IntType>(seed->match[left_idx] - gap_open,
                                                                     seed->insert_row[left_idx] - gap_open),
                                                        seed->insert_col[left_idx] - gap_extend), insert_col[idx]);
            
#ifdef debug_banded_aligner_fill_matrix
            cerr << "[BAMatrix::fill_matrix]: score is now " << (int) insert_col[idx] << endl;
#endif
        }
        
        // don't handle final row edge case if we actually got it with the first row edge case
        if (seed_next_bottom_diag_iter != seed_next_top_diag_iter) {
#ifdef debug_banded_aligner_fill_matrix
            cerr << "[BAMatrix::fill_matrix]: edge case for final cell in column" << endl;
            cerr << "[BAMatrix::fill_matrix]: extending a match and column gap into matrix coord (" << seed_next_bottom_diag_iter << ", 0)" << ", rectangular coord coord (" << seed_next_bottom_diag_iter - top_diag << ", 0)" << endl;
#endif
            
            // may only be able to extend a match on last iteration
            idx = (seed_next_bottom_diag_iter - top_diag) * ncols;
            diag_idx = (seed_next_bottom_diag_iter - seed_next_top_diag) * seed_node_seq_len + seed_node_seq_len - 1;
            if (qual_adjusted) {
                match_score = score_mat[25 * base_quality[seed_next_bottom_diag_iter] + 5 * nt_table[node_seq[0]] + nt_table[read[seed_next_bottom_diag_iter]]];
            }
            else {
                match_score = score_mat[5 * nt_table[node_seq[0]] + nt_table[read[seed_next_bottom_diag_iter]]];
            }
            
#ifdef debug_banded_aligner_fill_matrix
            cerr << "[BAMatrix::fill_matrix]: extending match from rectangular coord (" << seed_next_bottom_diag_iter - seed_next_top_diag << ", " << seed_node_seq_len - 1 << ")" << " with match score " << (int) match_score << ", scores are " << (int) seed->match[diag_idx] << " (M), " << (int) seed->insert_row[diag_idx] << " (Ir), and " << (int) seed->insert_col[diag_idx] << " (Ic), current score is " << (int) match[idx] << endl;
#endif
            match[idx] = max<IntType>(match_score + max<IntType>(max<IntType>(seed->match[diag_idx],
                                                                              seed->insert_row[diag_idx]),
                                                                 seed->insert_col[diag_idx]), match[idx]);
            
            // can only extend column gap if the bottom of the matrix was hit in the last seed
            if (beyond_bottom_of_matrix) {
#ifdef debug_banded_aligner_fill_matrix
                cerr << "[BAMatrix::fill_matrix]: can also extend a column gap since already reached edge of matrix" << endl;
#endif
                left_idx = (seed_next_bottom_diag_iter - seed_next_top_diag + 2) * seed_node_seq_len - 1;
                insert_col[idx] = max<IntType>(max<IntType>(max<IntType>(seed->match[left_idx] - gap_open,
                                                                         seed->insert_row[left_idx] - gap_open),
                                                            seed->insert_col[left_idx] - gap_extend), insert_col[idx]);
            }
        }
    }
    
    // POA into left hand column from the seeds
    if (treat_as_source && ncols > 0) {
        if (cumulative_seq_len != 0) {
            cerr << "error:[BandedGlobalAligner] banded alignment has no node predecessor for node in middle of path" << endl;
            exit(1);
        }
        
#ifdef debug_banded_aligner_fill_matrix
        cerr << "[BAMatrix::fill_matrix]: this is a source node, computing implied edge conditions" << endl;
#endif
        
        // first node in alignment, fill out initial column based on implied lead gaps
        
        // find position of the first cell in the rectangularized band
        int64_t iter_start = -top_diag;
        idx = iter_start * ncols;
        
        // cap stop index if last diagonal is below bottom of matrix
        int64_t iter_stop = bottom_diag >= (int64_t) read.length() ? band_height + (int64_t) read.length() - bottom_diag - 1 : band_height;
        // match of first nucleotides
        if (qual_adjusted) {
            match[idx] = max<IntType>(score_mat[25 * base_quality[0] + 5 * nt_table[node_seq[0]] + nt_table[read[0]]], match[idx]);
            
#ifdef debug_banded_aligner_fill_matrix
            cerr << "[BAMatrix::fill_matrix]: set quality adjusted initial match cell to " << (int) match[idx] << " from node char " << node_seq[0] << ", read char " << read[0] << ", base qual " << (int) base_quality[0] << " for adjusted matrix index " << 25 * base_quality[0] + 5 * nt_table[node_seq[0]] + nt_table[read[0]] << " and score " << (int) score_mat[25 * base_quality[0] + 5 * nt_table[node_seq[0]] + nt_table[read[0]]] << endl;
#endif
        }
        else {
             match[idx] = max<IntType>(match[idx] = score_mat[5 * nt_table[node_seq[0]] + nt_table[read[0]]], match[idx]);
        }
        
        // only way to end an alignment in a gap here is to row and column gap
        insert_row[idx] = max<IntType>(-2 * gap_open, insert_row[idx]);
        insert_col[idx] = max<IntType>(-2 * gap_open, insert_col[idx]);
        
        for (int64_t i = iter_start + 1; i < iter_stop; i++) {
            idx = i * ncols;
            up_idx = idx - ncols;
            // score of a match in this cell
            IntType match_score;
            if (qual_adjusted) {
                match_score = score_mat[25 * base_quality[top_diag + i] + 5 * nt_table[node_seq[0]] + nt_table[read[top_diag + i]]];
            }
            else {
                match_score = score_mat[5 * nt_table[node_seq[0]] + nt_table[read[top_diag + i]]];
            }
            // must take one lead gap to get into first column
            match[idx] = max<IntType>(match_score - gap_open - (top_diag + i - 1) * gap_extend, match[idx]);
            // normal iteration along column
            insert_row[idx] = max<IntType>(max<IntType>(match[up_idx] - gap_open, insert_row[up_idx] - gap_extend),
                                           insert_col[up_idx] - gap_open);
            // must take two gaps to get into first column
            insert_col[idx] = max<IntType>(-2 * gap_open - (top_diag + i) * gap_extend, insert_col[idx]);
#ifdef debug_banded_aligner_fill_matrix
            cerr << "[BAMatrix::fill_matrix]: on left edge of matrix at rectangle coords (" << i << ", " << 0 << "), match score of node char " << 0 << " (" << node_seq[0] << ") and read char " << i + top_diag << " (" << read[i + top_diag] << ") is " << (int) match_score << ", leading gap length is " << top_diag + i << " for total match matrix score of " << (int) match[idx] << endl;
#endif
        }
        
        // fix the final insert column value, which actually would have been inserting from outside the band
        insert_col[idx] = min_inf;
    }
    else if (ncols > 0){
        // compute the insert row scores without any cases for lead gaps (these can be safely computed after
        // the POA iterations since they do not cross node boundaries)
        for (int64_t i = iter_start + 1; i < iter_stop; i++) {
            idx = i * ncols;
            up_idx = (i - 1) * ncols;
            
            insert_row[idx] = max<IntType>(max<IntType>(match[up_idx] - gap_open, insert_row[up_idx] - gap_extend),
                                           insert_col[up_idx] - gap_open);
        }
    }
    
#ifdef debug_banded_aligner_fill_matrix
    cerr << "[BAMatrix::fill_matrix]: seeding finished, moving to subsequent columns" << endl;
#endif
    
    // iterate through the rest of the columns
    for (int64_t j = 1; j < ncols; j++) {
        
        // are we clipping any diagonals because they are outside the range of the matrix in this column?
        bool bottom_diag_outside = bottom_diag + j >= int64_t(read.size());
        bool top_diag_outside = top_diag + j < 0;
        bool top_diag_abutting = top_diag + j == 0;
        
        int64_t iter_start = top_diag_outside ? -(top_diag + j) : 0;
        int64_t iter_stop = bottom_diag_outside ? band_height + int64_t(read.size()) - bottom_diag - j - 1 : band_height;
        
        idx = iter_start * ncols + j;
        
        IntType match_score;
        if (qual_adjusted) {
            match_score = score_mat[25 * base_quality[iter_start + top_diag + j] + 5 * nt_table[node_seq[j]] + nt_table[read[iter_start + top_diag + j]]];
        }
        else {
            match_score = score_mat[5 * nt_table[node_seq[j]] + nt_table[read[iter_start + top_diag + j]]];
        }
        if (top_diag_outside || top_diag_abutting) {
            // match after implied gap along top edge
            match[idx] = match_score - gap_open - (cumulative_seq_len + j - 1) * gap_extend;
            
#ifdef debug_banded_aligner_fill_matrix
            cerr << "[BAMatrix::fill_matrix]: on upper edge of matrix at rectangle coords (" << iter_start << ", " << j << "), match score of node char " << j << " (" << node_seq[j] << ") and read char " << iter_start + top_diag + j << " (" << read[iter_start + top_diag + j] << ") is " << (int) match_score << ", leading gap length is " << cumulative_seq_len + j << " for total match matrix score of " << (int) match[idx] << endl;
#endif
        }
        else {
            diag_idx = iter_start * ncols + (j - 1);
            // cells should be present to do normal diagonal iteration
            match[idx] = match_score + max(max(match[diag_idx], insert_row[diag_idx]), insert_col[diag_idx]);
        }
        
        if (top_diag_outside) {
            // gap open after implied gap along top edge
            insert_row[idx] = -2 * gap_open - (cumulative_seq_len + j) * gap_extend;
        }
        else {
            // cannot reach this node with row insert (outside the diagonal)
            insert_row[idx] = min_inf;
        }
        
        // normal iteration along row unless band height is 1
        if (band_height != 1) {
            int64_t left_idx = (iter_start + 1) * ncols + (j - 1);
            insert_col[idx] = max(max(match[left_idx] - gap_open, insert_row[left_idx] - gap_open),
                                  insert_col[left_idx] - gap_extend);
        }
        else {
            insert_col[idx] = min_inf;
        }
        
        
        for (int64_t i = iter_start + 1; i < iter_stop - 1; i++) {
            // indices of the current and previous cells in the rectangularized band
            idx = i * ncols + j;
            up_idx = (i - 1) * ncols + j;
            diag_idx = i * ncols + (j - 1);
            left_idx = (i + 1) * ncols + (j - 1);
            
            if (qual_adjusted) {
                match_score = score_mat[25 * base_quality[i + top_diag + j] + 5 * nt_table[node_seq[j]] + nt_table[read[i + top_diag + j]]];
            }
            else {
                match_score = score_mat[5 * nt_table[node_seq[j]] + nt_table[read[i + top_diag + j]]];
            }
            
            match[idx] = match_score + max(max(match[diag_idx], insert_row[diag_idx]), insert_col[diag_idx]);
            
            insert_row[idx] = max(max(match[up_idx] - gap_open, insert_row[up_idx] - gap_extend),
                                  insert_col[up_idx] - gap_open);
            
            insert_col[idx] = max(max(match[left_idx] - gap_open, insert_row[left_idx] - gap_open),
                                  insert_col[left_idx] - gap_extend);
            
#ifdef debug_banded_aligner_fill_matrix
            cerr << "[BAMatrix::fill_matrix]: in interior of matrix at rectangle coords (" << i << ", " << j << "), match score of node char " << j << " (" << node_seq[j] << ") and read char " << i + top_diag + j << " (" << read[i + top_diag + j] << ") is " << (int) match_score << ", leading gap length is " << cumulative_seq_len + j << " for total match matrix score of " << (int) match[idx] << endl;
#endif
        }
        
        // stop iteration one cell early to handle logic on bottom edge of band
        
        // skip this step in edge case where read length is 1
        if (iter_stop - 1 > iter_start) {
            idx = (iter_stop - 1) * ncols + j;
            up_idx = (iter_stop - 2) * ncols + j;
            diag_idx = (iter_stop - 1) * ncols + (j - 1);
            
            if (qual_adjusted) {
                match_score = score_mat[25 * base_quality[iter_stop + top_diag + j - 1] + 5 * nt_table[node_seq[j]] + nt_table[read[iter_stop + top_diag + j - 1]]];
            }
            else {
                match_score = score_mat[5 * nt_table[node_seq[j]] + nt_table[read[iter_stop + top_diag + j - 1]]];
            }
            
            match[idx] = match_score + max(max(match[diag_idx], insert_row[diag_idx]), insert_col[diag_idx]);
            
            insert_row[idx] = max(max(match[up_idx] - gap_open, insert_row[up_idx] - gap_extend),
                                  insert_col[up_idx] - gap_open);
            
            if (bottom_diag_outside) {
                // along the bottom edge of the matrix, so the cell to the right is still there
                left_idx = iter_stop * ncols + (j - 1);
                insert_col[idx] = max(max(match[left_idx] - gap_open, insert_row[left_idx] - gap_open),
                                      insert_col[left_idx] - gap_extend);
                
            }
            else {
                // cell to the right is outside the band
                insert_col[idx] = min_inf;
            }
        }
    }
    
#ifdef debug_banded_aligner_print_matrices
    print_full_matrices(graph);
    print_rectangularized_bands(graph);
#endif
}


template <class IntType>
void BandedGlobalAligner<IntType>::BAMatrix::init_traceback_indexes(const HandleGraph& graph, int64_t& i, int64_t& j) {
    
    // get coordinates of bottom right corner
    const string& read = alignment.sequence();
    int64_t ncols = graph.get_length(node);
    i = bottom_diag + ncols > (int64_t) read.length() ? (int64_t) read.length() - top_diag - ncols : bottom_diag - top_diag;
    j = ncols - 1;
}

template <class IntType>
void BandedGlobalAligner<IntType>::BAMatrix::traceback(const HandleGraph& graph, BABuilder& builder,
                                                       AltTracebackStack& traceback_stack,
                                                       int64_t& i, int64_t& j, matrix_t& mat, bool& in_lead_gap,
                                                       const int8_t* score_mat, const int8_t* nt_table,
                                                       const int8_t gap_open,  const int8_t gap_extend,
                                                       const bool qual_adjusted, const IntType min_inf) {
    
#ifdef debug_banded_aligner_traceback
    cerr << "[BAMatrix::traceback] starting traceback back through node " << graph.get_id(node) << " from rectangular coordinates (" << i << ", " << j << "), currently " << (in_lead_gap ? "" : "not ") << "in a lead gap" << endl;
#endif
    
    const string& read = alignment.sequence();
    const string& base_quality = alignment.quality();
    
    int64_t band_height = bottom_diag - top_diag + 1;
    string node_seq = graph.get_sequence(node);
    int64_t ncols = node_seq.size();
    int64_t node_id = graph.get_id(node);
    
    int64_t idx, next_idx;
    IntType curr_score;
    IntType source_score;
    IntType score_diff;
    IntType alt_score;
    IntType curr_traceback_score = traceback_stack.current_traceback_score();
    
    // do node traceback unless we are in the lead gap implied at the edge of the DP matrix or we
    // are already at a node boundary trying to get across
    while ((j > 0 || mat == InsertRow) && !in_lead_gap) {
        
#ifdef debug_banded_aligner_traceback
        cerr << "[BAMatrix::traceback] traceback coordinates (" << i << ", " << j << "), current matrix is " << (mat == Match ? "match" : (mat == InsertCol ? "insert column" : "insert row")) << endl;
#endif
        
        // add traceback step to alignment
        builder.update_state(graph, mat, node, i + top_diag + j, j);
        
        // check for a deflection
        if (traceback_stack.at_next_deflection(node_id, i, j)) {
            
            // move to the next position as dictated by current matrix
            switch (mat) {
                case Match:
                    j--;
                    break;
                    
                case InsertRow:
                    i--;
                    break;
                    
                case InsertCol:
                    j--;
                    i++;
                    break;
                    
                default:
                    cerr << "error:[BandedGlobalAligner] unrecognized matrix type at traceback deflection" << endl;
                    assert(0);
                    break;
            }
            
            // take the deflection and advance the deflection iterator
            mat = traceback_stack.deflect_to_matrix();
            
#ifdef debug_banded_aligner_traceback
            cerr << "[BAMatrix::traceback] taking inside matrix deflection to " << (mat == Match ? "match" : (mat == InsertCol ? "insert column" : "insert row")) << endl;
#endif
            
            continue;
        }
        
        array<matrix_t, 3> prev_mats{Match, InsertCol, InsertRow};
        
        // find optimal traceback
        idx = i * ncols + j;
        bool found_trace = false;
        switch (mat) {
            case Match:
            {
                if (i + j == -top_diag) {
#ifdef debug_banded_aligner_traceback
                    cerr << "[BAMatrix::traceback] next cell is outside matrix, opening implied lead gap" << endl;
#endif
                    // at top of matrix, move into implied lead gap along top edge
                    mat = InsertCol;
                    j--;
                    in_lead_gap = true;
                    // no where else to go, so break out of switch statement without checking alts
                    break;
                }
                
                curr_score = match[idx];
                next_idx = i * ncols + j - 1;
                
                IntType match_score;
                if (qual_adjusted) {
                    match_score = score_mat[25 * base_quality[i + top_diag + j] + 5 * nt_table[node_seq[j]] + nt_table[read[i + top_diag + j]]];
                }
                else {
                    match_score = score_mat[5 * nt_table[node_seq[j]] + nt_table[read[i + top_diag + j]]];
                }
                
#ifdef debug_banded_aligner_traceback
                cerr << "[BAMatrix::traceback] transitioning from match, current score " << (int) match[idx] << " match/mismatch score " << (int) match_score << " from node char " << j << " (" << node_seq[j] << ") and read char " << i + top_diag + j << " (" << read[i + top_diag + j] << ")" << endl;
#endif
                
                for (auto prev_mat : prev_mats) {
                    
                    switch (prev_mat) {
                        case Match:
                        {
                            source_score = match[next_idx];
                            score_diff = curr_score - (source_score + match_score);
                            if (score_diff == 0) {
#ifdef debug_banded_aligner_traceback
                                cerr << "[BAMatrix::traceback] found next cell in match matrix with score " << (int) match[next_idx] << endl;
#endif
                                mat = Match;
                                found_trace = true;
                            }
                            else if (source_score != min_inf) {
                                alt_score = curr_traceback_score - score_diff;
                                traceback_stack.propose_deflection(alt_score, node_id, i, j, node_id, Match);
                            }
                            break;
                        }
                        case InsertRow:
                        {
                            source_score = insert_row[next_idx];
                            if (source_score > min_inf) {
                                score_diff = curr_score - (source_score + match_score);
                                if (!found_trace && score_diff == 0) {
#ifdef debug_banded_aligner_traceback
                                    cerr << "[BAMatrix::traceback] found next cell in insert row matrix with score " << (int) insert_row[next_idx] << endl;
#endif
                                    mat = InsertRow;
                                    found_trace = true;
                                }
                                else {
                                    alt_score = curr_traceback_score - score_diff;
                                    traceback_stack.propose_deflection(alt_score, node_id, i, j, node_id, InsertRow);
                                }
                            }
                            break;
                        }
                        case InsertCol:
                        {
                            source_score = insert_col[next_idx];
                            if (source_score > min_inf) {
                                score_diff = curr_score - (source_score + match_score);
                                if (!found_trace && score_diff == 0) {
#ifdef debug_banded_aligner_traceback
                                    cerr << "[BAMatrix::traceback] found next cell in insert column matrix with score " << (int) insert_col[next_idx] << endl;
#endif
                                    mat = InsertCol;
                                    found_trace = true;
                                }
                                else if (source_score != min_inf) {
                                    alt_score = curr_traceback_score - score_diff;
                                    traceback_stack.propose_deflection(alt_score, node_id, i, j, node_id, InsertCol);
                                }
                                
                            }
                            break;
                        }
                        default:
                        {
                            cerr << "error: invalid previous matrix" << endl;
                            exit(1);
                        }
                    }
                }
                
                if (!found_trace) {
                    cerr << "error:[BandedGlobalAligner] traceback stuck in match matrix interior" << endl;
                    assert(0);
                }
                
                j--;
                
                break;
            }
                
            case InsertRow:
            {
                if (i == 0) {
                    // along top of band
                    cerr << "error:[BandedGlobalAligner] traceback attempted to leave band from top" << endl;
                    assert(0);
                }
                
                if (i + j == -top_diag) {
                    // at top of matrix, move into implied lead gap along top edge
                    in_lead_gap = true;
                    // no where else to go, so break out of switch statement without checking alts
                    i--;
                    break;
                }
                
                curr_score = insert_row[idx];
                next_idx = (i - 1) * ncols + j;

                for (auto prev_mat : prev_mats) {
                    
                    switch (prev_mat) {
                        case Match:
                        {
                            source_score = match[next_idx];
                            score_diff = curr_score - (source_score - gap_open);
                            if (score_diff == 0) {
                                mat = Match;
                                found_trace = true;
                            }
                            else if (source_score != min_inf) {
                                alt_score = curr_traceback_score - score_diff;
                                traceback_stack.propose_deflection(alt_score, node_id, i, j, node_id, Match);
                            }
                            break;
                        }
                        case InsertRow:
                        {
                            source_score = insert_row[next_idx];
                            if (source_score > min_inf) {
                                score_diff = curr_score - (source_score - gap_extend);
                                if (!found_trace && score_diff == 0) {
                                    mat = InsertRow;
                                    found_trace = true;
                                }
                                else if (source_score != min_inf) {
                                    alt_score = curr_traceback_score - score_diff;
                                    traceback_stack.propose_deflection(alt_score, node_id, i, j, node_id, InsertRow);
                                }
                            }
                            break;
                        }
                        case InsertCol:
                        {
                            source_score = insert_col[next_idx];
                            if (source_score > min_inf) {
                                score_diff = curr_score - (source_score - gap_open);
                                if (!found_trace && score_diff == 0) {
                                    mat = InsertCol;
                                    found_trace = true;
                                }
                                else if (source_score != min_inf) {
                                    alt_score = curr_traceback_score - score_diff;
                                    traceback_stack.propose_deflection(alt_score, node_id, i, j, node_id, InsertCol);
                                }
                            }
                            break;
                        }
                        default:
                        {
                            cerr << "error: invalid previous matrix" << endl;
                            exit(1);
                        }
                    }
                }
                
                if (!found_trace) {
                    cerr << "error:[BandedGlobalAligner] traceback stuck in insert row matrix interior" << endl;
                    assert(0);
                }
                
                i--;
                
                break;
            }
            
            case InsertCol:
            {
                if (i == band_height - 1) {
                    // along bottom of band
                    cerr << "error:[BandedGlobalAligner] traceback attempted to leave band from bottom" << endl;
                    assert(0);
                    
                }
                
                curr_score = insert_col[idx];
                next_idx = (i + 1) * ncols + j - 1;

                for (auto prev_mat : prev_mats) {
                    switch (prev_mat) {
                        case Match:
                        {
                            source_score = match[next_idx];
                            score_diff = curr_score - (source_score - gap_open);
                            if (score_diff == 0) {
                                mat = Match;
                                found_trace = true;
                            }
                            else if (source_score != min_inf) {
                                alt_score = curr_traceback_score - score_diff;
                                traceback_stack.propose_deflection(alt_score, node_id, i, j, node_id, Match);
                            }
                            break;
                        }
                        case InsertRow:
                        {
                            source_score = insert_row[next_idx];
                            if (source_score > min_inf) {
                                score_diff = curr_score - (source_score - gap_open);
                                if (!found_trace && score_diff == 0) {
                                    mat = InsertRow;
                                    found_trace = true;
                                }
                                else if (source_score != min_inf) {
                                    alt_score = curr_traceback_score - score_diff;
                                    traceback_stack.propose_deflection(alt_score, node_id, i, j, node_id, InsertRow);
                                }
                            }
                            break;
                        }
                        case InsertCol:
                        {
                            source_score = insert_col[next_idx];
                            if (source_score > min_inf) {
                                score_diff = curr_score - (source_score - gap_extend);
                                if (!found_trace && score_diff == 0) {
                                    mat = InsertCol;
                                    found_trace = true;
                                }
                                else if (source_score != min_inf) {
                                    alt_score = curr_traceback_score - score_diff;
                                    traceback_stack.propose_deflection(alt_score, node_id, i, j, node_id, InsertCol);
                                }
                            }
                            break;
                        }
                        default:
                        {
                            cerr << "error: invalid previous matrix" << endl;
                            exit(1);
                        }
                    }
                }
                
                if (!found_trace) {
                    cerr << "error:[BandedGlobalAligner] traceback stuck in insert column matrix interior" << endl;
                    assert(0);
                }
                
                i++;
                j--;
                
                break;
            }
                
            default:
            {
                cerr << "error:[BandedGlobalAligner] unrecognized matrix type given to traceback" << endl;
                assert(0);
                break;
            }
        }
    }
    
    if (in_lead_gap) {
#ifdef debug_banded_aligner_traceback
        cerr << "[BAMatrix::traceback] running through node sequence in a lead gap" << endl;
#endif
        // add lead column gaps until reaching edge of node
        mat = InsertCol;
        while (j > 0) {
            builder.update_state(graph, mat, node, -1, j);
            j--;
            i++;
        }
    }
    
}

template<class IntType>
void BandedGlobalAligner<IntType>::BAMatrix::traceback_over_edge(const HandleGraph& graph, BABuilder& builder,
                                                                 AltTracebackStack& traceback_stack,
                                                                 int64_t& i, int64_t& j, matrix_t& mat,
                                                                 bool& in_lead_gap, int64_t& node_id,
                                                                 const int8_t* score_mat, const int8_t* nt_table,
                                                                 const int8_t gap_open, const int8_t gap_extend,
                                                                 const bool qual_adjusted, IntType const min_inf) {
    
    // begin POA across the boundary
    
    bool treat_as_source = false;
    unordered_set<int64_t> traceback_source_nodes;
    if (seeds.empty()) {
        treat_as_source = true;
        traceback_source_nodes.insert(graph.get_id(node));
    }
    
    BAMatrix* traceback_seed = nullptr;
    int64_t traceback_seed_row = std::numeric_limits<int64_t>::min();
    int64_t traceback_seed_col = std::numeric_limits<int64_t>::min();
    matrix_t traceback_mat = Match;
    
    const string& read = alignment.sequence();
    const string& base_quality = alignment.quality();
    
    int64_t idx, next_idx;
    IntType score_diff;
    IntType alt_score;
    IntType curr_score;
    IntType source_score;
    IntType curr_traceback_score = traceback_stack.current_traceback_score();
    
    int64_t curr_diag = top_diag + i;
    string node_seq = graph.get_sequence(node);
    int64_t ncols = graph.get_length(node);
    
    if (traceback_stack.at_next_deflection(node_id, i, j)) {
#ifdef debug_banded_aligner_traceback
        cerr << "[BAMatrix::traceback_over_edge] at boundary, taking a deflection" << endl;
#endif
        
        builder.update_state(graph, mat, node, curr_diag, 0);
        
        // where to deflect to?
        int64_t deflect_node_id;
        matrix_t deflect_matrix = traceback_stack.deflect_to_matrix(deflect_node_id);
        
        // find which seed matrix to deflect to (don't have a better way of looking this up right now)
        vector<BAMatrix*> seed_path;
        for (auto initial_seed : seeds) {
            
            vector<BAMatrix*> seed_stack{initial_seed};
            
            while (!seed_stack.empty()) {
                BAMatrix* seed = seed_stack.back();
                seed_stack.pop_back();
                
                if (!seed) {
                    // pop off the path if we hit the stack marker
                    seed_path.pop_back();
                }
                else if (graph.get_id(seed->node) == deflect_node_id) {
                    // we found the traceback node
                    seed_path.push_back(seed);
                    break;
                }
                else if (graph.get_length(seed->node) == 0) {
                    // this is not the traceback node, but it is an empty node so the traceback
                    // might be on the other side of it
                    seed_path.push_back(seed);
                    seed_stack.push_back(nullptr);
                    for (auto seed_of_seed : seed->seeds) {
                        seed_stack.push_back(seed_of_seed);
                    }
                }
            }
            
            // stop looking if we found a path to the traceback seed
            if (!seed_path.empty()) {
                break;
            }
        }
        
        if (seed_path.empty()) {
            cerr << "error:[BandedGlobalAligner] traceback node boundary could not find node to deflect to" << endl;
            assert(0);
        }
        
        // if we traversed empty nodes on the way to the traceback node, add empty mapping for those
        if (seed_path.size() > 1) {
            auto end = seed_path.end();
            end--;
            for (auto iter = seed_path.begin(); iter != end; iter++) {
                builder.update_state(graph, mat, (*iter)->node, i, 0, true);
            }
        }
        
        BAMatrix* seed = seed_path.back();
        
        int64_t seed_ncols = graph.get_length(seed->node);
        traceback_seed_row = curr_diag - seed->top_diag - seed_ncols + (mat == InsertCol);
        traceback_seed_col = seed_ncols - 1;
        // check whether we're crossing into a lead gap
        in_lead_gap = (curr_diag == 0 && mat == Match) || in_lead_gap;
        
#ifdef debug_banded_aligner_traceback
        cerr << "[BAMatrix::traceback_over_edge] taking node boundary deflection to " << (deflect_matrix == Match ? "match" : (deflect_matrix == InsertCol ? "insert column" : "insert row")) << " in node " << deflect_node_id << ", will start at coordinates (" << traceback_seed_row << ", " << traceback_seed_col << ")" << endl;
#endif
        
        i = traceback_seed_row;
        j = traceback_seed_col;
        mat = deflect_matrix;
        node_id = deflect_node_id;
        return;
    }
    
    bool found_trace = false;
    // if we traverse through nodes with no sequence, we need to keep track of which ones
    vector<BAMatrix*> empty_intermediate_nodes;
    vector<BAMatrix*> empty_source_path;
    
    // a queue of seeds and their empty predecessors
    vector<pair<BAMatrix*, vector<BAMatrix*>>> seed_queue;
    for (auto seed : seeds) {
        seed_queue.emplace_back(seed, vector<BAMatrix*>());
    }
    
    if (in_lead_gap) {
#ifdef debug_banded_aligner_traceback
        cerr << "[BAMatrix::traceback_over_edge] at boundary, following seed backward from a lead gap" << endl;
#endif
        // we are in the implied lead gap along the top of the matrix
        // take the shortest path back to the origin of the global alignment
        
        // add final read deletion of node
        builder.update_state(graph, mat, node, -1, j);
        
        while (!seed_queue.empty()) {
            
            auto seed_record = seed_queue.back();
            BAMatrix* seed = seed_record.first;
            seed_queue.pop_back();
            
            // is seed masked?
            if (seed == nullptr) {
                continue;
            }
            
#ifdef debug_banded_aligner_traceback
            cerr << "[BAMatrix::traceback_over_edge] checking seed node " << graph.get_id(seed->node) << endl;
#endif
            
            // if this node is empty, add its predecessors to the queue
            if (graph.get_length(seed->node) == 0) {
#ifdef debug_banded_aligner_traceback
                cerr << "[BAMatrix::traceback_over_edge] seed node " << graph.get_id(seed->node) << " is empty, checking predecessors" << endl;
#endif
                // record that this seed comes before its predecessors in the traceback
                seed_record.second.push_back(seed);
                if (seed->seeds.empty()) {
#ifdef debug_banded_aligner_traceback
                    cerr << "[BAMatrix::traceback_over_edge] empty seed node " << graph.get_id(seed->node) << " is a source" << endl;
#endif
                    treat_as_source = true;
                    traceback_source_nodes.insert(graph.get_id(seed->node));
                    empty_source_path = seed_record.second;
                }
                
                for (auto seed_of_seed : seed->seeds) {
                    seed_queue.push_back(make_pair(seed_of_seed, seed_record.second));
                }
                continue;
            }
            
            score_diff = gap_extend * (seed->cumulative_seq_len + graph.get_length(seed->node) - cumulative_seq_len);
            if (score_diff == 0 && !found_trace) {
#ifdef debug_banded_aligner_traceback
                cerr << "[BAMatrix::traceback_over_edge] found a lead gap traceback to node " << graph.get_id(seed->node) << endl;
#endif
                traceback_seed = seed;
                empty_intermediate_nodes = seed_record.second;
                found_trace = true;
            }
            else {
                alt_score = curr_traceback_score - score_diff;
                traceback_stack.propose_deflection(alt_score, node_id, i, j, graph.get_id(seed->node), InsertCol);
            }
        }
        
        if (traceback_seed) {
            // where in the matrix is this?
            int64_t seed_ncols = graph.get_length(traceback_seed->node);
            int64_t seed_extended_top_diag = traceback_seed->top_diag + seed_ncols;
            traceback_seed_row = top_diag - seed_extended_top_diag + i + 1;
            traceback_seed_col = seed_ncols - 1;
        }
    }
    else {
        
        builder.update_state(graph, mat, node, curr_diag, 0);
        
        IntType match_score;
        switch (mat) {
            case Match:
            {
                curr_score = match[i * ncols];
                if (qual_adjusted) {
                    match_score = score_mat[25 * base_quality[i + top_diag] + 5 * nt_table[node_seq[j]] + nt_table[read[i + top_diag]]];
                }
                else {
                    match_score = score_mat[5 * nt_table[node_seq[j]] + nt_table[read[i + top_diag]]];
                }
                break;
            }
                
            case InsertCol:
            {
                curr_score = insert_col[i * ncols];
                break;
            }
                
            case InsertRow:
            {
                cerr << "error:[BandedGlobalAligner] traceback attempted to open row gap over node boundary" << endl;
                assert(0);
                break;
            }
                
            default:
            {
                cerr << "error:[BandedGlobalAligner] unrecognized matrix type given to traceback" << endl;
                assert(0);
                break;
            }
        }
        
#ifdef debug_banded_aligner_traceback
        cerr << "[BAMatrix::traceback_over_edge] at boundary, in node " << graph.get_id(node) << " following seed backward from " << (mat == Match ? "match" : "insert column") << " matrix with score " << (int) curr_score << endl;
#endif
        
        // matches stay on same diagonal, column insertions move over one diagonal
        // note that the indexing is relative to THIS matrix, not the seed
        
        // check traceback goes to each seed matrix
        while (!seed_queue.empty()) {
            auto seed_record = seed_queue.back();
            BAMatrix* seed = seed_record.first;
            seed_queue.pop_back();
            
            // is the matrix masked?
            if (seed == nullptr) {
#ifdef debug_banded_aligner_traceback
                cerr << "[BAMatrix::traceback_over_edge] seed is masked" << endl;
#endif
                continue;
            }
            
            if (graph.get_length(seed->node) == 0) {
                
                for (auto seed_of_seed : seed->seeds) {
                    seed_queue.push_back(make_pair(seed_of_seed, seed_record.second));
                    seed_queue.back().second.push_back(seed);
                }
                
                if (seed->seeds.empty()) {
                    treat_as_source = true;
                    // keep track of the path through empty nodes to a source
                    empty_source_path = seed_record.second;
                    empty_source_path.push_back(seed);
                }
                continue;
            }
            
            
#ifdef debug_banded_aligner_traceback
            cerr << "[BAMatrix::traceback_over_edge] checking seed node " << graph.get_id(seed->node) << endl;
#endif
            
            int64_t seed_node_id = graph.get_id(seed->node);
            int64_t seed_ncols = graph.get_length(seed->node);
            
            // the diagonals in the current matrix that this seed extends to
            int64_t seed_extended_top_diag = seed->top_diag + seed_ncols;
            int64_t seed_extended_bottom_diag = seed->bottom_diag + seed_ncols;
            
            // does the traceback diagonal extend backward to this matrix?
            if (curr_diag > seed_extended_bottom_diag - (mat == InsertCol) // col inserts hit 1 less of band
                || curr_diag < seed_extended_top_diag) {
#ifdef debug_banded_aligner_traceback
                cerr << "[BAMatrix::traceback_over_edge] seed extended diags are top: " << seed_extended_top_diag << ", bottom: " << seed_extended_bottom_diag << " and curr mat is " << (mat == InsertCol ? "" : "not") << " insert column, so we cannot extend from this seed to the current diag " << curr_diag << endl;
#endif
                continue;
            }
            
            int64_t seed_col = seed_ncols - 1;
            int64_t seed_row = -(seed_extended_top_diag - top_diag) + i + (mat == InsertCol);
            next_idx = seed_row * seed_ncols + seed_col;
            
#ifdef debug_banded_aligner_traceback
            cerr << "[BAMatrix::traceback_over_edge] checking seed rectangular coordinates (" << seed_row << ", " << seed_col << "), with indices calculated from current diagonal " << curr_diag << " (top diag " << top_diag << " + offset " << i << "), seed top diagonal " << seed->top_diag << ", seed seq length " << seed_ncols << " with insert column offset " << (mat == InsertCol) << endl;
#endif
            array<matrix_t, 3> prev_mats{Match, InsertCol, InsertRow};
            
            switch (mat) {
                case Match:
                {
#ifdef debug_banded_aligner_traceback
                    cerr << "[BAMatrix::traceback_over_edge] poa backwards from match, seed extended top diag " << seed_extended_top_diag << endl;
#endif
                    // does match lead into a lead row gap?
                    if (seed->top_diag + seed_row + seed_col == -1) {
#ifdef debug_banded_aligner_traceback
                        cerr << "[BAMatrix::traceback_over_edge] traceback points to a lead column gap of length " << seed->cumulative_seq_len + seed_ncols << " with score " << (int) -gap_open - (seed->cumulative_seq_len + seed_ncols - 1) * gap_extend << " extending to score here of " << (int) curr_score << " with match score " << (int) match_score << endl;
#endif
                        // score of implied column gap
                        source_score = -gap_open - (seed->cumulative_seq_len + seed_ncols - 1) * gap_extend;
                        score_diff = curr_score - (source_score + match_score);
                        if (score_diff == 0 && !found_trace) {
                            traceback_mat = InsertCol;
                            traceback_seed = seed;
                            traceback_seed_row = seed_row;
                            traceback_seed_col = seed_col;
                            in_lead_gap = true;
                            found_trace = true;
                            empty_intermediate_nodes = seed_record.second;
#ifdef debug_banded_aligner_traceback
                            cerr << "[BAMatrix::traceback_over_edge] hit found in lead gap with score " << -gap_open - (seed->cumulative_seq_len + seed_ncols - 1) * gap_extend << endl;
#endif
                        }
                        else {
                            alt_score = curr_traceback_score - score_diff;
                            traceback_stack.propose_deflection(alt_score, node_id, i, j, seed_node_id, InsertCol);
                        }
                        
                        // don't check any of the matrices because they will have garbage values in this position or seg fault
                        break;
                    }
                    
                    for (auto prev_mat : prev_mats) {
                        switch (prev_mat) {
                            case Match:
                            {
                                source_score = seed->match[next_idx];
                                // don't need to check edge condition because match does not have min inf
                                score_diff = curr_score - (source_score + match_score);
                                if (score_diff == 0 && !found_trace) {
                                    traceback_mat = Match;
                                    traceback_seed = seed;
                                    traceback_seed_row = seed_row;
                                    traceback_seed_col = seed_col;
                                    found_trace = true;
                                    empty_intermediate_nodes = seed_record.second;
#ifdef debug_banded_aligner_traceback
                                    cerr << "[BAMatrix::traceback_over_edge] hit found in match matrix with score " << (int) seed->match[next_idx] << endl;
#endif
                                }
                                else if (source_score != min_inf) {
                                    alt_score = curr_traceback_score - score_diff;
                                    traceback_stack.propose_deflection(alt_score, node_id, i, j, seed_node_id, Match);
                                }
                                break;
                            }
                            case InsertRow:
                            {
                                source_score = seed->insert_row[next_idx];
                                // check edge condition
                                if (source_score > min_inf) {
                                    score_diff = curr_score - (source_score + match_score);
                                    if (score_diff == 0 && !found_trace) {
#ifdef debug_banded_aligner_traceback
                                        cerr << "[BAMatrix::traceback_over_edge] hit found in insert row matrix  with score " << (int) seed->insert_row[next_idx] << endl;
#endif
                                        traceback_mat = InsertRow;
                                        traceback_seed = seed;
                                        traceback_seed_row = seed_row;
                                        traceback_seed_col = seed_col;
                                        found_trace = true;
                                        empty_intermediate_nodes = seed_record.second;
                                    }
                                    else {
                                        alt_score = curr_traceback_score - score_diff;
                                        traceback_stack.propose_deflection(alt_score, node_id, i, j, seed_node_id, InsertRow);
                                    }
                                }
                                break;
                            }
                            case InsertCol:
                            {
                                source_score = seed->insert_col[next_idx];
                                // check edge condition
                                if (source_score > min_inf) {
                                    score_diff = curr_score - (source_score + match_score);
                                    if (score_diff == 0 && !found_trace) {
#ifdef debug_banded_aligner_traceback
                                        cerr << "[BAMatrix::traceback_over_edge] hit found in insert column matrix  with score " << (int) seed->insert_col[next_idx] << endl;
#endif
                                        traceback_mat = InsertCol;
                                        traceback_seed = seed;
                                        traceback_seed_row = seed_row;
                                        traceback_seed_col = seed_col;
                                        found_trace = true;
                                        empty_intermediate_nodes = seed_record.second;
                                    }
                                    else {
                                        alt_score = curr_traceback_score - score_diff;
                                        traceback_stack.propose_deflection(alt_score, node_id, i, j, seed_node_id, InsertCol);
                                    }
                                }
                                break;
                            }
                            default:
                            {
                                cerr << "error: invalid matrix type" << endl;
                                exit(1);
                            }
                        }
                    }
                    
                    break;
                }
                    
                case InsertCol:
                {
                    for (auto prev_mat : prev_mats) {
                        switch (prev_mat) {
                            case Match:
                            {
                                source_score = seed->match[next_idx];
                                // don't need to check edge condition because match does not have min inf
                                score_diff = curr_score - (source_score - gap_open);
                                if (score_diff == 0 && !found_trace) {
#ifdef debug_banded_aligner_traceback
                                    cerr << "[BAMatrix::traceback_over_edge] hit found in match matrix with score " << (int) seed->match[next_idx] << endl;
#endif
                                    traceback_mat = Match;
                                    traceback_seed = seed;
                                    traceback_seed_row = seed_row;
                                    traceback_seed_col = seed_col;
                                    found_trace = true;
                                    empty_intermediate_nodes = seed_record.second;
                                }
                                else if (source_score != min_inf) {
                                    alt_score = curr_traceback_score - score_diff;
#ifdef debug_banded_aligner_traceback
                                    cerr << "[BAMatrix::traceback_over_edge] no hit in match matrix, proposing deflection with alt score " << (int) alt_score << " from current traceback score " << (int) curr_traceback_score << " and score diff " << (int) score_diff << endl;
#endif
                                    traceback_stack.propose_deflection(alt_score, node_id, i, j, seed_node_id, Match);
                                }
                                break;
                            }
                            case InsertRow:
                            {
                                source_score = seed->insert_row[next_idx];
                                // check edge condition
                                if (source_score > min_inf) {
                                    score_diff = curr_score - (source_score - gap_open);
                                    if (score_diff == 0 && !found_trace) {
#ifdef debug_banded_aligner_traceback
                                        cerr << "[BAMatrix::traceback_over_edge] hit found in insert row matrix with score " << (int) seed->match[next_idx] << endl;
#endif
                                        traceback_mat = InsertRow;
                                        traceback_seed = seed;
                                        traceback_seed_row = seed_row;
                                        traceback_seed_col = seed_col;
                                        found_trace = true;
                                        empty_intermediate_nodes = seed_record.second;
                                    }
                                    else {
                                        alt_score = curr_traceback_score - score_diff;
#ifdef debug_banded_aligner_traceback
                                        cerr << "[BAMatrix::traceback_over_edge] no hit in insert column matrix, proposing deflection with alt score " << (int) alt_score << " from current traceback score " << (int) curr_traceback_score << " and score diff " << (int) score_diff << endl;
#endif
                                        traceback_stack.propose_deflection(alt_score, node_id, i, j, seed_node_id, InsertRow);
                                    }
                                }
                                break;
                            }
                            case InsertCol:
                            {
                                source_score = seed->insert_col[next_idx];
                                // check edge condition
                                if (source_score > min_inf) {
                                    score_diff = curr_score - (source_score - gap_extend);
                                    if (score_diff == 0 && !found_trace) {
#ifdef debug_banded_aligner_traceback
                                        cerr << "[BAMatrix::traceback_over_edge] hit found in insert column matrix with score " << (int) seed->match[next_idx] << endl;
#endif
                                        traceback_mat = InsertCol;
                                        traceback_seed = seed;
                                        traceback_seed_row = seed_row;
                                        traceback_seed_col = seed_col;
                                        found_trace = true;
                                        empty_intermediate_nodes = seed_record.second;
                                    }
                                    else {
                                        alt_score = curr_traceback_score - score_diff;
#ifdef debug_banded_aligner_traceback
                                        cerr << "[BAMatrix::traceback_over_edge] no hit in insert row matrix, proposing deflection with alt score " << (int) alt_score << " from current traceback score " << (int) curr_traceback_score << " and score diff " << (int) score_diff << endl;
#endif
                                        traceback_stack.propose_deflection(alt_score, node_id, i, j, seed_node_id, InsertCol);
                                    }
                                }
                                break;
                            }
                            default:
                            {
                                cerr << "error: invalid matrix type" << endl;
                                exit(1);
                            }
                        }
                    }
                    
                    break;
                }
                    
                case InsertRow:
                {
                    cerr << "error:[BandedGlobalAligner] illegal matrix type for moving across node boundary" << endl;
                    assert(0);
                    break;
                }
                    
                default:
                {
                    cerr << "error:[BandedGlobalAligner] unrecognized matrix type given to traceback" << endl;
                    assert(0);
                    break;
                }
            }
        }
    }
    
    bool found_source_trace = false;
    
    if (treat_as_source) {
        // this is a source node, or it is connected to one by a zero-length path
        
#ifdef debug_banded_aligner_traceback
        cerr << "[BAMatrix::traceback_over_edge] at beginning of first node in alignment" << endl;
#endif
        if (in_lead_gap && !found_trace) {
            // this will always be the shortest gap
            found_source_trace = true;
            j--;
            i++;
        }
        else {
            switch (mat) {
                case Match:
                {
                    IntType match_score;
                    if (qual_adjusted) {
                        match_score = score_mat[25 * base_quality[i + top_diag + j] + 5 * nt_table[node_seq[j]] + nt_table[read[i + top_diag + j]]];
                    }
                    else {
                        match_score = score_mat[5 * nt_table[node_seq[j]] + nt_table[read[i + top_diag + j]]];
                    }
                    
                    source_score = curr_diag > 0 ? -gap_open - (curr_diag - 1) * gap_extend : 0;
                    score_diff = curr_score - (source_score + match_score);
                    if (score_diff == 0 && !found_trace) {
#ifdef debug_banded_aligner_traceback
                        cerr << "[BAMatrix::traceback_over_edge] alignment starts with match, adding read char " << top_diag + i << ": " << read[top_diag + i] << endl;
#endif
                        mat = InsertRow;
                        found_source_trace = true;
                    }
                    else {
                        alt_score = curr_traceback_score - score_diff;
                        for (int64_t source_node_id : traceback_source_nodes) {
                            traceback_stack.propose_deflection(alt_score, node_id, i, j, source_node_id, InsertRow);
                        }
                    }
                    j--;
                    break;
                }
                    
                case InsertCol:
                {
                    source_score = -gap_open - (i + top_diag) * gap_extend;
                    score_diff = curr_score - (source_score - gap_open);
                    if (score_diff == 0 && !found_trace) {
#ifdef debug_banded_aligner_traceback
                        cerr << "[BAMatrix::traceback_over_edge] alignment starts with column gap" << endl;
#endif
                        mat = InsertRow;
                        found_source_trace = true;
                    }
                    else {
                        alt_score = curr_traceback_score - score_diff;
                        for (int64_t source_node_id : traceback_source_nodes) {
                            traceback_stack.propose_deflection(alt_score, node_id, i, j, source_node_id, InsertRow);
                        }
                    }
                    j--;
                    i++;
                    break;
                }
                    
                default:
                {
                    cerr << "error:[BandedGlobalAligner] invalid matrix type for final traceback column" << endl;
                    assert(0);
                    break;
                }
            }
        }
    }
    
    
    if (found_source_trace) {
        // if we traversed any empty nodes before finding the traceback, add empty updates for them
        for (BAMatrix* seed_path_node : empty_source_path) {
            builder.update_state(graph, InsertCol, seed_path_node->node, i + top_diag, 0, true);
        }
        
        // add any lead row gaps necessary
        while (top_diag + i > 0) {
#ifdef debug_banded_aligner_traceback
            cerr << "[BAMatrix::traceback_over_edge] initial row gaps are present, adding read char " << top_diag + i - 1 << ": " << read[top_diag + i - 1] << endl;
#endif
            i--;
            const handle_t& end_node = empty_source_path.empty() ? node : empty_source_path.back()->node;
            builder.update_state(graph, InsertRow, end_node, i + top_diag, -1, graph.get_length(end_node) == 0);
        }
        
        // set the node ID to 0 to indicate completion
        node_id = 0;
        return;
    }
    else {
#ifdef debug_banded_aligner_traceback
        cerr << "[BAMatrix::traceback_over_edge] traversed " << empty_intermediate_nodes.size() << " empty nodes before finding trace" << endl;
#endif
        // if we traversed any empty nodes before finding the traceback, add empty updates for them
        for (BAMatrix* intermediate_node : empty_intermediate_nodes) {
            builder.update_state(graph, mat, intermediate_node->node, i + top_diag, 0, true);
        }
        
    }
    
    if (!found_trace) {
        cerr << "error:[BandedGlobalAligner] traceback stuck at node boundary" << endl;
        assert(0);
    }
    
    // set the traceback values for the next matrix
    i = traceback_seed_row;
    j = traceback_seed_col;
    mat = traceback_mat;
    node_id = graph.get_id(traceback_seed->node);
}

template <class IntType>
void BandedGlobalAligner<IntType>::BAMatrix::print_full_matrices(const HandleGraph& graph) {
    if (match == nullptr) {
        cerr << "error:[BandedGlobalAligner] cannot print matrix before performing dynamic programming" << endl;
        assert(0);
    }
    
    cerr << "matrices for node " << graph.get_id(node) << ":" << endl;
    
    for (matrix_t mat : {Match, InsertRow, InsertCol}) {
        print_matrix(graph, mat);
    }

}

template <class IntType>
void BandedGlobalAligner<IntType>::BAMatrix::print_rectangularized_bands(const HandleGraph& graph) {
    if (match == nullptr) {
        cerr << "error:[BandedGlobalAligner] cannot print band before performing dynamic programming" << endl;
        assert(0);
    }
    
    cerr << "rectangularized bands for node " << graph.get_id(node) << ":" << endl;
    
    for (matrix_t mat : {Match, InsertRow, InsertCol}) {
        print_band(graph, mat);
    }
}

template <class IntType>
void BandedGlobalAligner<IntType>::BAMatrix::print_matrix(const HandleGraph& graph, matrix_t which_mat) {

    const string& read = alignment.sequence();
    string node_seq = graph.get_sequence(node);
    
    IntType* band_rect;
    switch (which_mat) {
        case Match:
            cerr << "match:" << endl;
            band_rect = match;
            break;
            
        case InsertRow:
            cerr << "insert row:" << endl;
            band_rect = insert_row;
            break;
            
        case InsertCol:
            cerr << "insert column:" << endl;
            band_rect = insert_col;
            break;
            
        default:
            cerr << "error:[BandedGlobalAligner] unrecognized matrix type" << endl;
            assert(0);
            break;
    }
    
    for (auto iter = node_seq.begin(); iter != node_seq.end(); iter++) {
        cerr << "\t" << *iter;
    }
    cerr << endl;
    
    int64_t ncols = node_seq.length();
    
    for (int64_t i = 0; i < (int64_t) read.length(); i++) {
        cerr << read[i];
        for (int64_t j = 0; j < ncols; j++) {
            int64_t diag = i - j;
            if (diag < top_diag || diag > bottom_diag) {
                cerr << "\t.";
            }
            else {
                cerr << "\t" << (int) band_rect[(diag - top_diag) * ncols + j];
            }
        }
        cerr << endl;
    }
}

template <class IntType>
void BandedGlobalAligner<IntType>::BAMatrix::print_band(const HandleGraph& graph, matrix_t which_mat) {
    
    const string& read = alignment.sequence();
    string node_seq = graph.get_sequence(node);
    
    IntType* band_rect;
    switch (which_mat) {
        case Match:
            cerr << "match:" << endl;
            band_rect = match;
            break;
            
        case InsertRow:
            cerr << "insert row:" << endl;
            band_rect = insert_row;
            break;
            
        case InsertCol:
            cerr << "insert column:" << endl;
            band_rect = insert_col;
            break;
            
        default:
            cerr << "error:[BandedGlobalAligner] unrecognized matrix type" << endl;
            assert(0);
            break;
    }
    
    for (auto iter = node_seq.begin(); iter != node_seq.end(); iter++) {
        cerr << "\t" << *iter;
    }
    cerr << endl;
    
    int64_t band_height = bottom_diag - top_diag + 1;
    int64_t ncols = node_seq.length();
    
    for (int64_t i = 0; i < band_height; i++) {
        
        if (top_diag + i > 0) {
            cerr << read[top_diag + i - 1];
        }
        
        for (int64_t j = 0; j < ncols; j++) {
            if (i + j < -top_diag || top_diag + i + j >= (int64_t) read.length()) {
                cerr << "\t.";
            }
            else {
                cerr << "\t" << (int) band_rect[i * ncols + j];
            }
        }
        cerr << endl;
    }
    
    
    for (int64_t i = bottom_diag; i < min(bottom_diag + ncols, (int64_t) read.length()); i++) {
        cerr << read[i] << "\t";
    }
    cerr << endl;
}

template <class IntType>
BandedGlobalAligner<IntType>::BandedGlobalAligner(Alignment& alignment, const HandleGraph& g,
                                                  int64_t band_padding, bool permissive_banding,
                                                  bool adjust_for_base_quality) :
                                                  BandedGlobalAligner(alignment, g,
                                                                      nullptr, 1,
                                                                      band_padding,
                                                                      permissive_banding,
                                                                      adjust_for_base_quality)
{
    // nothing to do, just funnel into internal constructor
}

template <class IntType>
BandedGlobalAligner<IntType>::BandedGlobalAligner(Alignment& alignment, const HandleGraph& g,
                                                  vector<Alignment>& alt_alignments,
                                                  int64_t max_multi_alns, int64_t band_padding,
                                                  bool permissive_banding,
                                                  bool adjust_for_base_quality) :
                                                  BandedGlobalAligner(alignment, g,
                                                                      &alt_alignments,
                                                                      max_multi_alns,
                                                                      band_padding,
                                                                      permissive_banding,
                                                                      adjust_for_base_quality)
{
    // check data integrity and funnel into internal constructor
    if (!alt_alignments.empty()) {
        cerr << "error:[BandedGlobalAligner] alternate alignment vector must be empty before aligning" << endl;
        assert(0);
    }
}

template <class IntType>
BandedGlobalAligner<IntType>::BandedGlobalAligner(Alignment& alignment, const HandleGraph& g,
                                                  vector<Alignment>* alt_alignments,
                                                  int64_t max_multi_alns,
                                                  int64_t band_padding,
                                                  bool permissive_banding,
                                                  bool adjust_for_base_quality) :
                                                  graph(g),
                                                  alignment(alignment),
                                                  alt_alignments(alt_alignments),
                                                  max_multi_alns(max_multi_alns),
                                                  adjust_for_base_quality(adjust_for_base_quality),
                                                  // compute some graph features we will be frequently reusing
                                                  topological_order(handlealgs::lazier_topological_order(&g)),
                                                  source_nodes(handlealgs::head_nodes(&g)),
                                                  sink_nodes(handlealgs::tail_nodes(&g))
{
#ifdef debug_banded_aligner_objects
    cerr << "[BandedGlobalAligner]: constructing BandedBlobalAligner with " << band_padding << " padding, " << permissive_banding << " permissive, " << adjust_for_base_quality << " quality adjusted" << endl;
#endif
    if (adjust_for_base_quality) {
        if (alignment.quality().empty() && !alignment.sequence().empty()) {
            cerr << "error:[BandedGlobalAligner] alignment needs base quality to perform quality adjusted alignment" << endl;
            assert(0);
        }
    }
    
    if (topological_order.size() < graph.get_node_count()) {
        cerr << "error:[BandedGlobalAligner] alignment graph must be a DAG" << endl;
    }
    
    // TODO: this can waste memory, but reallocating the vector seems to throw an error in protobuf and
    // we won't know if there are fewer alignments than the max until the cycle is over
    if (alt_alignments) {
        alt_alignments->reserve(max_multi_alns);
    }
    
    // map node ids to indices
    for (int64_t i = 0; i < topological_order.size(); i++) {
        node_id_to_idx[graph.get_id(topological_order[i])] = i;
    }
    
#ifdef debug_banded_aligner_objects
    cerr << "[BandedGlobalAligner]: " << source_nodes.size() << " sources and " << sink_nodes.size() << " sinks" << endl;
    cerr << "[BandedGlobalAligner]: computing node bands" << endl;
#endif
    
    // figure out what the bands need to be for alignment and which nodes cannot complete a
    // global alignment within the band
    vector<bool> node_masked;
    vector<pair<int64_t, int64_t>> band_ends;
    find_banded_paths(permissive_banding, band_padding, node_masked, band_ends);
    
#ifdef debug_banded_aligner_objects
    cerr << "[BandedGlobalAligner]: identifying shortest paths" << endl;
#endif
    
    // find the shortest sequence leading to each node so we can infer the length
    // of lead deletions
    vector<int64_t> shortest_seqs;
    shortest_seq_paths(shortest_seqs);
    
#ifdef debug_banded_aligner_objects
    cerr << "[BandedGlobalAligner]: constructing banded matrix objects" << endl;
#endif
    
    // initialize DP matrices for each node
    banded_matrices.resize(graph.get_node_count(), nullptr);
    for (int64_t i = 0; i < topological_order.size(); i++) {
        
#ifdef debug_banded_aligner_objects
        cerr << "[BandedGlobalAligner]: creating matrix object for node " << graph.get_id(topological_order[i]) << " at index " << i << endl;
#endif
        
        if (!node_masked[i])  {
            
            const handle_t& node = topological_order[i];
            
            int64_t node_seq_len = graph.get_length(node);
            
#ifdef debug_banded_aligner_objects
            cerr << "[BandedGlobalAligner]: establishing seed list for node " << graph.get_id(node) << " at index " << i << endl;
#endif
            
            // POA predecessor matrices
            vector<BAMatrix*> seeds;
            graph.follow_edges(node, true, [&](const handle_t& prev) {
                seeds.push_back(banded_matrices[node_id_to_idx[graph.get_id(prev)]]);
            });
            
            banded_matrices[i] = new BAMatrix(alignment,
                                              node,
                                              band_ends[i].first,
                                              band_ends[i].second,
                                              std::move(seeds),
                                              shortest_seqs[i]);
            
        }
    }
    
    // check to see if we chose a banding that is too restrictive for making an
    // an alignment
    if (!permissive_banding) {
        bool sinks_masked = true;
        for (const handle_t& node : sink_nodes) {
            if (banded_matrices[node_id_to_idx[graph.get_id(node)]] != nullptr) {
                sinks_masked = false;
                break;
            }
            
        }
        if (sinks_masked) {
            // We couldn't find an alignment in this band. That's bad, but we
            // don't necessarily want to kill the whole program.
            throw NoAlignmentInBandException();
        }
    }
}

template <class IntType>
BandedGlobalAligner<IntType>::~BandedGlobalAligner() {
    
    for (BAMatrix* banded_matrix : banded_matrices) {
        if (banded_matrix != nullptr) {
            delete banded_matrix;
        }
    }
}

template <class IntType>
void BandedGlobalAligner<IntType>::path_lengths_to_sinks(vector<int64_t>& shortest_path_to_sink,
                                                         vector<int64_t>& longest_path_to_sink) {
#ifdef debug_banded_aligner_graph_processing
    cerr << "[BandedGlobalAligner::path_lengths_to_sinks]: finding longest and shortest paths to sink node" << endl;
#endif
    
    // find the longest path from the right side of each matrix to the end of the graph
    
    // set initial values
    longest_path_to_sink.resize(topological_order.size(), 0);
    shortest_path_to_sink.resize(topological_order.size(), numeric_limits<int64_t>::max());
    
    // set base case (longest path already set to 0)
    for (const handle_t& handle : sink_nodes) {
        shortest_path_to_sink[node_id_to_idx.at(graph.get_id(handle))] = 0;
    }
    
    // iterate in reverse order
    for (int64_t i = topological_order.size() - 1; i >= 0; i--) {
        
        int64_t node_seq_len = graph.get_length(topological_order[i]);
        // compute longest path through this node to right side of incoming matrices
        int64_t longest_path_length = longest_path_to_sink[i] + node_seq_len;
        int64_t shortest_path_length = shortest_path_to_sink[i] + node_seq_len;
        
        graph.follow_edges(topological_order[i], true, [&](const handle_t& prev) {
            
            int64_t prev_idx = node_id_to_idx.at(graph.get_id(prev));
            
            if (longest_path_to_sink[prev_idx] < longest_path_length) {
                
#ifdef debug_banded_aligner_graph_processing
                cerr << "[BandedGlobalAligner::path_lengths_to_sinks]: path through " << graph.get_id(prev) << " of length " << longest_path_length << " to node at index " << prev_idx << " is longer than current longest path " << longest_path_to_sink[prev_idx] << ", updating it now" << endl;
#endif
                
                longest_path_to_sink[prev_idx] = longest_path_length;
            }
            if (shortest_path_to_sink[prev_idx] > shortest_path_length) {
                
#ifdef debug_banded_aligner_graph_processing
                cerr << "[BandedGlobalAligner::path_lengths_to_sinks]: path through " << graph.get_id(prev) << " of length " << shortest_path_length << " to node at index " << prev_idx << " is shorter than current shortest path " << shortest_path_to_sink[prev_idx] << ", updating it now" << endl;
#endif
                
                shortest_path_to_sink[prev_idx] = shortest_path_length;
            }
        });
    }
}


// fills vectors with whether nodes are masked by the band width, and the band ends of each node
template <class IntType>
void BandedGlobalAligner<IntType>::find_banded_paths(bool permissive_banding, int64_t band_padding,
                                                     vector<bool>& node_masked,
                                                     vector<pair<int64_t, int64_t>>& band_ends) {
    
    // keeps track of which nodes cannot reach the bottom corner within the band
    node_masked.resize(topological_order.size(), false);
    
    // the bottom and top indices of the band in the rightmost column of each node's matrix
    band_ends.resize(topological_order.size(), make_pair(numeric_limits<int64_t>::max(),
                                                         numeric_limits<int64_t>::min()));
    
    // find the longest and shortest path from each node to any sink
    vector<int64_t> shortest_path_to_sink;
    vector<int64_t> longest_path_to_sink;
    path_lengths_to_sinks(shortest_path_to_sink, longest_path_to_sink);
    
    if (permissive_banding) {
        // initialize with wide enough bands that every source can hit every connected sink
        for (const handle_t& init_node : source_nodes) {
            int64_t init_node_idx = node_id_to_idx.at(graph.get_id(init_node));
            int64_t init_node_seq_len = graph.get_length(init_node);
            band_ends[init_node_idx].first = min(-band_padding,
                                                 int64_t(alignment.sequence().size())
                                                 - (init_node_seq_len + longest_path_to_sink[init_node_idx]) - band_padding);
            band_ends[init_node_idx].second = max(band_padding,
                                                  int64_t(alignment.sequence().size())
                                                  - (init_node_seq_len + shortest_path_to_sink[init_node_idx]) + band_padding);
#ifdef debug_banded_aligner_graph_processing
            cerr << "[BandedGlobalAligner::find_banded_paths]: initializing band path end at node " << graph.get_id(init_node) << " at index " << init_node_idx << " to top " << band_ends[init_node_idx].first << ", and bottom " << band_ends[init_node_idx].second << " from shortest and longest paths of length " << shortest_path_to_sink[init_node_idx] << " and " << longest_path_to_sink[init_node_idx] << " compared to read length " << alignment.sequence().size() << " with padding " << band_padding << endl;
#endif
        }
        
    }
    else {
        // initialize with band ends beginning with source nodes
        for (const handle_t& handle : source_nodes) {
            int64_t init_node_idx = node_id_to_idx.at(graph.get_id(handle));
            int64_t init_node_seq_len = graph.get_length(handle);
            band_ends[init_node_idx].first = -band_padding;
            band_ends[init_node_idx].second = band_padding;
#ifdef debug_banded_aligner_graph_processing
            cerr << "[BandedGlobalAligner::find_banded_paths]: initializing band path end at node " << graph.get_id(handle) << " at index " << init_node_idx << " to top " << band_ends[init_node_idx].first << ", and bottom " << band_ends[init_node_idx].second << endl;
#endif
        }
    }
    
    // iterate through the rest of the nodes in topological order
    for (int64_t i = 0; i < topological_order.size(); i++) {
        const handle_t& node = topological_order[i];
        int64_t node_seq_len = graph.get_length(node);
        
        int64_t extended_band_top = band_ends[i].first + node_seq_len;
        int64_t extended_band_bottom = band_ends[i].second + node_seq_len;
        
#ifdef debug_banded_aligner_graph_processing
        cerr << "[BandedGlobalAligner::find_banded_paths]: following edges out of node " << graph.get_id(node) << " at index " << i << " with sequence " << graph.get_sequence(node) << ", band of " << band_ends[i].first << ", " << band_ends[i].second << " extending to " << extended_band_top << ", " << extended_band_bottom << endl;
#endif
        // can alignments from this node reach the bottom right corner within the band?
        if (extended_band_top + shortest_path_to_sink[i] > int64_t(alignment.sequence().size())
            || extended_band_bottom + longest_path_to_sink[i] < int64_t(alignment.sequence().size())) {
            
            node_masked[i] = true;
            
#ifdef debug_banded_aligner_graph_processing
            cerr << "[BandedGlobalAligner::find_banded_paths]: cannot complete alignment to read of length " << alignment.sequence().size() << " along shortest path " << shortest_path_to_sink[i] << " or longest path " << longest_path_to_sink[i] << ", which reach range " << extended_band_top + shortest_path_to_sink[i] << ", " << extended_band_bottom + longest_path_to_sink[i] << endl;
#endif
            continue;
        }
        
        graph.follow_edges(node, false, [&](const handle_t& next) {
            
            int64_t node_out_idx = node_id_to_idx.at(graph.get_id(next));
            
#ifdef debug_banded_aligner_graph_processing
            cerr << "[BandedGlobalAligner::find_banded_paths]: extending band to node at index " << node_out_idx << endl;
#endif
            
            if (extended_band_top < band_ends[node_out_idx].first) {
#ifdef debug_banded_aligner_graph_processing
                cerr << "[BandedGlobalAligner::find_banded_paths]: updating top band limit from  " << band_ends[node_out_idx].first << " to " << extended_band_top << endl;
#endif
                band_ends[node_out_idx].first = extended_band_top;
            }
            
            if (extended_band_bottom > band_ends[node_out_idx].second) {
#ifdef debug_banded_aligner_graph_processing
                cerr << "[BandedGlobalAligner::find_banded_paths]: updating bottom band limit from  " << band_ends[node_out_idx].second << " to " << extended_band_bottom << endl;
#endif
                band_ends[node_out_idx].second = extended_band_bottom;
            }
        });
    }
}


// returns the shortest sequence from any source node to each node
template <class IntType>
void BandedGlobalAligner<IntType>::shortest_seq_paths(vector<int64_t>& seq_lens_out) {
    
    // initialize vector with min identity to store sequence lengths
    seq_lens_out.resize(topological_order.size(), numeric_limits<int64_t>::max());
    
    // base cases
    for (const handle_t& handle : source_nodes) {
        seq_lens_out[node_id_to_idx[graph.get_id(handle)]] = 0;
    }
    
    // dynamic programming to calculate sequence lengths for rest of nodes
    for (size_t i = 0; i < topological_order.size(); i++) {
        int64_t seq_len = graph.get_length(topological_order[i]) + seq_lens_out[i];

        graph.follow_edges(topological_order[i], false, [&](const handle_t& handle) {
            int64_t target_idx = node_id_to_idx.at(graph.get_id(handle));
            // find the shortest sequence that can reach the top left corner of the matrix
            if (seq_len < seq_lens_out[target_idx]) {
                seq_lens_out[target_idx] = seq_len;
            }
        });
    }
}

template <class IntType>
void BandedGlobalAligner<IntType>::align(int8_t* score_mat, int8_t* nt_table, int8_t gap_open, int8_t gap_extend) {
    
    // small enough number to never be accepted in alignment but also not trigger underflow
    IntType max_mismatch = numeric_limits<IntType>::max();
    for (int i = 0; i < 25; i++) {
        max_mismatch = min<IntType>(max_mismatch, score_mat[i]);
    }
    IntType min_inf = numeric_limits<IntType>::min() + max<IntType>((IntType) -max_mismatch, max<IntType>(gap_open, gap_extend));
    
    
    // fill each nodes matrix in topological order
    for (int64_t i = 0; i < banded_matrices.size(); i++) {
        BAMatrix* band_matrix = banded_matrices[i];
        
        // skip masked nodes
        if (band_matrix == nullptr) {
#ifdef debug_banded_aligner_fill_matrix
            cerr << "[BandedGlobalAligner::align] node " << graph.get_id(topological_order[i]) << " is masked, skipping" << endl;
#endif
            continue;
        }
        
#ifdef debug_banded_aligner_fill_matrix
        cerr << "[BandedGlobalAligner::align] at node " << graph.get_id(band_matrix->node) << " at index " << i << " with sequence " << graph.get_id(band_matrix->node) << endl;
        cerr << "[BandedGlobalAligner::align] node is not masked, filling matrix" << endl;
#endif
        band_matrix->fill_matrix(graph, score_mat, nt_table, gap_open, gap_extend, adjust_for_base_quality, min_inf);
    }
    
    traceback(score_mat, nt_table, gap_open, gap_extend, min_inf);
}

template <class IntType>
void BandedGlobalAligner<IntType>::traceback(int8_t* score_mat, int8_t* nt_table, int8_t gap_open, int8_t gap_extend, IntType min_inf) {
    
    // get the sink and source node matrices for alignment stack
    unordered_set<BAMatrix*> sink_node_matrices;
    unordered_set<BAMatrix*> source_node_matrices;
    for (const handle_t& node : sink_nodes) {
        sink_node_matrices.insert(banded_matrices[node_id_to_idx[graph.get_id(node)]]);
    }
    for (const handle_t& node : source_nodes) {
        source_node_matrices.insert(banded_matrices[node_id_to_idx[graph.get_id(node)]]);
    }
    
    int64_t read_length = alignment.sequence().length();
    int32_t empty_score = read_length > 0 ? -gap_open - (read_length - 1) * gap_extend : 0;
    
    // find the optimal alignment(s) and initialize stack
    AltTracebackStack traceback_stack(graph, max_multi_alns, empty_score, source_node_matrices, sink_node_matrices,
                                      gap_open, gap_extend, min_inf);
    
    while (traceback_stack.has_next()) {
        
        Alignment* next_alignment;
        if (!alt_alignments) {
            next_alignment = &alignment;
        }
        else if (alt_alignments->empty()){
            // do primary alignment first
            next_alignment = &alignment;
        }
        else {
            // create new alternate alignment object
            alt_alignments->emplace_back();
            next_alignment = &(alt_alignments->back());
            next_alignment->set_sequence(alignment.sequence());
            next_alignment->set_quality(alignment.quality());
        }
        
        if (traceback_stack.next_is_empty()) {
            traceback_stack.next_empty_alignment(*next_alignment);
#ifdef debug_banded_aligner_traceback
            cerr << "[BandedGlobalAligner::traceback] taking the next full empty alignment" << endl;
            cerr << pb2json(*next_alignment) << endl;
#endif
        }
        else {
            
            // what node does the alignment start at
            int64_t node_id;
            matrix_t mat;
            traceback_stack.get_alignment_start(node_id, mat);
            
            // start the row and column trackers
            int64_t i, j;
            banded_matrices[node_id_to_idx[node_id]]->init_traceback_indexes(graph, i, j);
            
            // we only start in a gap if the sequence if there is nothing to align
            bool in_lead_gap = alignment.sequence().empty();
            
#ifdef debug_banded_aligner_traceback
            cerr << "[BandedGlobalAligner::traceback] beginning traceback ending at node " << node_id << " in matrix " << (mat == Match ? "match" : (mat == InsertCol ? "insert column" : "insert row")) << endl;
#endif
            
            
            // do traceback
            BABuilder builder(*next_alignment);
            
            while (node_id != 0) {
                int64_t node_idx = node_id_to_idx[node_id];
                // trace through the matrix
                banded_matrices[node_idx]->traceback(graph, builder, traceback_stack, i, j, mat, in_lead_gap, score_mat,
                                                     nt_table, gap_open, gap_extend, adjust_for_base_quality, min_inf);
                // trace over edges
                banded_matrices[node_idx]->traceback_over_edge(graph, builder, traceback_stack, i, j, mat, in_lead_gap,
                                                               node_id, score_mat, nt_table, gap_open, gap_extend,
                                                               adjust_for_base_quality, min_inf);
            }
            
            // construct the alignment path
            builder.finalize_alignment(traceback_stack.current_empty_prefix());
            
            // add score to alignment
            next_alignment->set_score(traceback_stack.current_traceback_score());
            
            // advance to the next traceback
            traceback_stack.next_traceback_alignment();
        }
        
        if (alt_alignments) {
            if (alt_alignments->empty()) {
                // copy the primary into the alternates
                alt_alignments->emplace_back(*next_alignment);
            }
        }
    }
}

template <class IntType>
BandedGlobalAligner<IntType>::AltTracebackStack::AltTracebackStack(const HandleGraph& graph,
                                                                   int64_t max_multi_alns,
                                                                   int32_t empty_score,
                                                                   unordered_set<BAMatrix*>& source_node_matrices,
                                                                   unordered_set<BAMatrix*>& sink_node_matrices,
                                                                   int8_t gap_open,
                                                                   int8_t gap_extend,
                                                                   IntType min_inf) :
                                                                   graph(graph),
                                                                   empty_score(empty_score),
                                                                   max_multi_alns(max_multi_alns)
{
    // an empty trace back prefix for the initial tracebacks
    vector<Deflection> null_prefix;
    
    // check tracebacks for alignments ending in all sink nodes
    for (BAMatrix* sink_matrix : sink_node_matrices) {
        if (sink_matrix == nullptr) {
            // This is a masked sink node. Skip it.
            continue;
        }
    
        if (sink_matrix->match == nullptr) {
            cerr << "error:[BandedGlobalAligner] must fill dynamic programming matrices before finding optimal score" << endl;
            exit(1);
        }
        
        list<BAMatrix*> band_stack{sink_matrix};
        list<int64_t> path;
        while (!band_stack.empty()) {
            BAMatrix* band_matrix = band_stack.back();
            band_stack.pop_back();
            
            if (!band_matrix) {
#ifdef debug_banded_aligner_traceback
                cerr << "[BandedGlobalAligner::traceback] found stack marker, pulling " << path.front() << " from path" << endl;
#endif
                path.pop_front();
                continue;
            }
            
            if (graph.get_length(band_matrix->node) == 0) {
                path.push_front(graph.get_id(band_matrix->node));
                band_stack.push_back(nullptr);
                
#ifdef debug_banded_aligner_traceback
                cerr << "[BandedGlobalAligner::traceback] traversing initial empty path on " << graph.get_id(band_matrix->node) << endl;
#endif
                
                // we went all the way from a source to a sink using only nodes with
                // no sequence, keep track of these later so we can decide later
                // whether they are sufficiently high scoring alignments to yield
                if (source_node_matrices.count(band_matrix)) {
                    empty_full_paths.push_back(path);
#ifdef debug_banded_aligner_traceback
                    cerr << "[BandedGlobalAligner::traceback] found empty full path" << endl;
                    for (auto nid : path ) {
                        cerr << "\t" << nid << endl;
                    }
#endif
                    continue;
                }
                
                for (auto seed : band_matrix->seeds) {
                    if (seed) {
                        band_stack.push_back(seed);
                    }
                }
            }
            else {

                // get the coordinates of the bottom right corner
                const handle_t& node = band_matrix->node;
                int64_t node_id = graph.get_id(node);
                int64_t read_length = band_matrix->alignment.sequence().length();
                int64_t ncols = graph.get_length(node);
                
#ifdef debug_banded_aligner_traceback
                cerr << "[BandedGlobalAligner::traceback] initializing tracebacks on node " << node_id << endl;
#endif
                
                int64_t final_col = ncols - 1;
                int64_t final_row = band_matrix->bottom_diag + ncols > read_length ? read_length - band_matrix->top_diag - ncols : band_matrix->bottom_diag - band_matrix->top_diag;
                
                int64_t final_idx = final_row * ncols + final_col;
                
                if (band_matrix->alignment.sequence().empty()) {
                    // if the read sequence is empty then we can only insert relative to the graph
                    size_t graph_length = band_matrix->cumulative_seq_len + graph.get_length(band_matrix->node);
                    IntType insert_score = graph_length ? (graph_length - 1) * (-gap_extend) - gap_open : 0;
                    insert_traceback(null_prefix, insert_score, node_id, final_row, final_col, node_id, InsertCol, path);
                }
                else {
                    // let the insert routine figure out which one is the best and which ones to keep in the stack
                    array<matrix_t, 3> mats{Match, InsertCol, InsertRow};
                    for (auto mat : mats) {
                        switch (mat)
                        {
                            case Match:
                            {
                                if (band_matrix->match[final_idx] != min_inf) {
                                    insert_traceback(null_prefix, band_matrix->match[final_idx],
                                                     node_id, final_row, final_col, node_id, Match, path);
                                }
                                break;
                            }
                            case InsertCol:
                            {
                                if (band_matrix->insert_row[final_idx] != min_inf) {
                                    insert_traceback(null_prefix, band_matrix->insert_row[final_idx],
                                                     node_id, final_row, final_col, node_id, InsertRow, path);
                                }
                                break;
                            }
                            case InsertRow:
                            {
                                if (band_matrix->insert_col[final_idx] != min_inf) {
                                    insert_traceback(null_prefix, band_matrix->insert_col[final_idx],
                                                     node_id, final_row, final_col, node_id, InsertCol, path);
                                }
                                break;
                            }
                            default:
                            {
                                cerr << "error: invalid matrix type" << endl;
                                exit(1);
                            }
                        }
                    }
                }
            }
        }
    }
    

    // initialize the traceback trackers for the optimal traceback
    curr_traceback = alt_tracebacks.begin();
    if (curr_traceback != alt_tracebacks.end()) {
        curr_deflxn = get<0>(*curr_traceback).begin();
    }
}

template <class IntType>
BandedGlobalAligner<IntType>::AltTracebackStack::~AltTracebackStack() {
    // nothing to do
}

template <class IntType>
inline void BandedGlobalAligner<IntType>::AltTracebackStack::get_alignment_start(int64_t& node_id, matrix_t& matrix) {
    
    // move to next traceback
    if (curr_traceback != alt_tracebacks.end()) {
        // get the start node the first deflection
        curr_deflxn = get<0>(*curr_traceback).begin();
        node_id = (*curr_deflxn).from_node_id;
        // get the matrix and advance to the next deflection
        matrix = (*curr_deflxn).to_matrix;
        curr_deflxn++;
    }
    else {
        // no more tracebacks in stack
        node_id = -1;
        matrix = Match;
    }
}

template <class IntType>
inline bool BandedGlobalAligner<IntType>::AltTracebackStack::has_next() {
    return curr_traceback != alt_tracebacks.end() || !empty_full_paths.empty();
}

template <class IntType>
inline void BandedGlobalAligner<IntType>::AltTracebackStack::next_traceback_alignment() {
    if (curr_deflxn != get<0>(*curr_traceback).end()) {
        cerr << "warning:[BandedGlobalAligner] moving on to next alternate alignment without taking all deflections" << endl;
    }
    curr_traceback++;
    // the alt tracebacks list keeps track of how many alternates are left, so if we've
    // used them all up, dump the empty alignments list
    if (curr_traceback == alt_tracebacks.end()) {
        empty_full_paths.clear();
    }
}

template <class IntType>
inline bool BandedGlobalAligner<IntType>::AltTracebackStack::next_is_empty()  {
    return curr_traceback == alt_tracebacks.end() ? true : empty_score >= get<1>(*curr_traceback) && !empty_full_paths.empty();
}

template <class IntType>
inline void BandedGlobalAligner<IntType>::AltTracebackStack::next_empty_alignment(Alignment& alignment) {
#ifdef debug_banded_aligner_traceback
    cerr << "[BandedGlobalAligner::next_empty_alignment] creating empty alignment" << endl;
#endif
    
    // score to a full insertion
    alignment.set_score(empty_score);
    
    // add all the nodes in the path with empty edits
    Path* path = alignment.mutable_path();
    list<int64_t>& node_id_path = empty_full_paths.front();
    int32_t rank = 1;
    for (int64_t node_id : node_id_path) {
        Mapping* mapping = path->add_mapping();
        mapping->mutable_position()->set_node_id(node_id);
        mapping->set_rank(rank);
        mapping->add_edit();
        rank++;
    }
    
    // add the insertion onto the first mapping
    Edit* first_edit = path->mutable_mapping(0)->mutable_edit(0);
    first_edit->set_to_length(alignment.sequence().length());
    first_edit->set_sequence(alignment.sequence());
    
#ifdef debug_banded_aligner_traceback
    cerr << "[BandedGlobalAligner::next_empty_alignment] dequeueing empty alignment" << endl;
#endif
    
    // remove the empty path we just added
    empty_full_paths.pop_front();
    
    // we used up one of the alloted multi alignments on an empty path, so reduce
    // the max size of the stack and if necessary remove one from the back
    max_multi_alns--;
    if (alt_tracebacks.size() > max_multi_alns) {
#ifdef debug_banded_aligner_traceback
        cerr << "[BandedGlobalAligner::next_empty_alignment] removing alternate traceback alignments from stack after taking empty alignment" << endl;
#endif
        auto back = alt_tracebacks.end();
        back--;
        if (curr_traceback == back) {
            // we're moving the back of the stack past the current alternate traceback,
            // so we're done, dump the rest of the empty paths
            alt_tracebacks.pop_back();
            curr_traceback = alt_tracebacks.end();
            empty_full_paths.clear();
        }
        else {
            alt_tracebacks.pop_back();
        }
    }
}

template <class IntType>
inline void BandedGlobalAligner<IntType>::AltTracebackStack::propose_deflection(const IntType score, const int64_t from_node_id,
                                                                                const int64_t row_idx, const int64_t col_idx,
                                                                                const int64_t to_node_id, const matrix_t to_matrix) {
    // only propose deflections if we're going through a new untraversed section of the traceback
    if (curr_deflxn != get<0>(*curr_traceback).end()) {
        return;
    }
    
    // is the score good enough to be put on the stack?
    if (score <= get<1>(alt_tracebacks.back()) && alt_tracebacks.size() >= max_multi_alns) {
        return;
    }
#ifdef debug_banded_aligner_traceback
    cerr << "[AltTracebackStack::propose_deflection] inserting a traceback from deflection from " << from_node_id << " (" << row_idx << "," << col_idx << ")" << " to node " << to_node_id << " matrix " << (to_matrix == Match ? "match" : (to_matrix == InsertCol ? "insert col" : "insert row")) << endl;
#endif
    
    insert_traceback(get<0>(*curr_traceback), score, from_node_id, row_idx, col_idx, to_node_id, to_matrix, current_empty_prefix());
}

template <class IntType>
inline void BandedGlobalAligner<IntType>::AltTracebackStack::insert_traceback(const vector<Deflection>& traceback_prefix,
                                                                              const IntType score, const int64_t from_node_id,
                                                                              const int64_t row_idx, const int64_t col_idx,
                                                                              const int64_t to_node_id, const matrix_t to_matrix,
                                                                              const list<int64_t>& empty_node_prefix) {
#ifdef debug_banded_aligner_traceback
    cerr << "[AltTracebackStack::insert_traceback] adding traceback with score " << (int) score << ", new deflection at (" << row_idx << ", " << col_idx << ") on node " << from_node_id << " to " << (to_matrix == Match ? "match" : (to_matrix == InsertRow ? "insert row" : "insert column")) << " matrix on node " << to_node_id << endl;
#endif
    
    // find position in stack where this should go
    auto insert_after = alt_tracebacks.rbegin();
    while (score > get<1>(*insert_after)) {
        insert_after++;
        if (insert_after == alt_tracebacks.rend()) {
            break;
        }
    }
    
    // insert if score is high enough or stack is not full yet
    if (insert_after != alt_tracebacks.rbegin() || alt_tracebacks.size() < max_multi_alns) {
        
        // create a new traceback here
        auto new_traceback = alt_tracebacks.emplace(insert_after.base(), vector<Deflection>(), score, empty_node_prefix);
        
        vector<Deflection>& deflections = get<0>(*new_traceback);
        
        // add the deflections from the prefix
        deflections.reserve(traceback_prefix.size() + 1);
        for (const Deflection& deflection : traceback_prefix) {
            deflections.push_back(deflection);
        }
        
        // add the final deflection
        deflections.emplace_back(from_node_id, row_idx, col_idx, to_node_id, to_matrix);
    }
    
    // remove lowest scoring traceback if stack is over capacity
    if (alt_tracebacks.size() > max_multi_alns) {
        alt_tracebacks.pop_back();
    }
    
#ifdef debug_banded_aligner_traceback
    cerr << "[AltTracebackStack::insert_traceback] scores in alt traceback stack currently ";
    for (auto trace : alt_tracebacks) {
        cerr << (int) get<1>(trace) << " -> ";
    }
    cerr << endl;
#endif
}

template <class IntType>
inline IntType BandedGlobalAligner<IntType>::AltTracebackStack::current_traceback_score() {
    return get<1>(*curr_traceback);
}

template <class IntType>
inline const list<int64_t>& BandedGlobalAligner<IntType>::AltTracebackStack::current_empty_prefix() {
    return get<2>(*curr_traceback);
}

template <class IntType>
inline bool BandedGlobalAligner<IntType>::AltTracebackStack::at_next_deflection(int64_t node_id, int64_t row_idx,
                                                                                int64_t col_idx) {
    // taken all deflections already?
    if (curr_deflxn == get<0>(*curr_traceback).end()) {
        return false;
    }
    
    // at the correct coordinates?
    return node_id == (*curr_deflxn).from_node_id &&
           row_idx == (*curr_deflxn).row_idx &&
           col_idx == (*curr_deflxn).col_idx;
}

template <class IntType>
inline typename BandedGlobalAligner<IntType>::matrix_t BandedGlobalAligner<IntType>::AltTracebackStack::deflect_to_matrix() {
    matrix_t mat = (*curr_deflxn).to_matrix;
    curr_deflxn++;
    return mat;
}

template <class IntType>
inline typename BandedGlobalAligner<IntType>::matrix_t BandedGlobalAligner<IntType>::AltTracebackStack::deflect_to_matrix(int64_t& to_node_id) {
    matrix_t mat = (*curr_deflxn).to_matrix;
    to_node_id = (*curr_deflxn).to_node_id;
    curr_deflxn++;
    return mat;
}

template <class IntType>
BandedGlobalAligner<IntType>::AltTracebackStack::Deflection::Deflection(const int64_t from_node_id,
                                                                        const int64_t row_idx,
                                                                        const int64_t col_idx,
                                                                        const int64_t to_node_id,
                                                                        const matrix_t to_matrix) :
                                                                        from_node_id(from_node_id),
                                                                        row_idx(row_idx),
                                                                        col_idx(col_idx),
                                                                        to_node_id(to_node_id),
                                                                        to_matrix(to_matrix)
{
    // only set values
}

template <class IntType>
BandedGlobalAligner<IntType>::AltTracebackStack::Deflection::~Deflection() {
    // nothing to do
}

const string NoAlignmentInBandException::message = "error:[BandedGlobalAligner] cannot align to graph within band, consider permissive banding";

const char* NoAlignmentInBandException::what() const noexcept {
    return message.c_str();
}

// Define aligners for allowed integer types.
// This MUST occur after the actual members of the types are defined, according to the standard.
// Only members with definitions we have seen already are supposed to be created.
template class BandedGlobalAligner<int8_t>;
template class BandedGlobalAligner<int16_t>;
template class BandedGlobalAligner<int32_t>;
template class BandedGlobalAligner<int64_t>;

}


