//
//  banded_global_aligner.cpp
//  
//  Contains objects to support Aligner in performing banded global alignment against
//  a graph.
//

#include "banded_global_aligner.hpp"
#include "json2pb.h"

//#define debug_banded_aligner_objects
//#define debug_banded_aligner_graph_processing
//#define debug_banded_aligner_fill_matrix
//#define debug_banded_aligner_traceback
//#define debug_banded_aligner_print_matrices

using namespace vg;

BandedGlobalAligner::BABuilder::BABuilder(Alignment& alignment) :
                                          alignment(alignment),
                                          matrix_state(Match),
                                          matching(false),
                                          current_node(nullptr),
                                          edit_length(0),
                                          edit_read_end_idx(0)
{
    // nothing to do
}

BandedGlobalAligner::BABuilder::~BABuilder() {
    // does not own any heap objects, nothing to do
}

void BandedGlobalAligner::BABuilder::update_state(matrix_t matrix, Node* node,
                                                  int64_t read_idx, int64_t node_idx) {
#ifdef debug_banded_aligner_traceback
    cerr << "[BABuilder::update_state] beginning state update for read index " << read_idx << ", node seq index " << node_idx << endl;
#endif
    if (node != current_node) {
#ifdef debug_banded_aligner_traceback
        cerr << "[BABuilder::update_state] at new node" << endl;
#endif
        // conclude current mapping and proceed to next node
        finish_current_node();
        current_node = node;
        matrix_state = matrix;
        if (matrix_state == Match) {
            matching = (alignment.sequence()[read_idx] == current_node->sequence()[node_idx]);
        }
        edit_length = 1;
        edit_read_end_idx = read_idx;
    }
    else if (matrix != matrix_state) {
#ifdef debug_banded_aligner_traceback
        cerr << "[BABuilder::update_state] transitioning into a new matrix" << endl;
#endif
        // transitioned into another matrix, finish current edit and begin new one
        finish_current_edit();
        matrix_state = matrix;
        if (matrix_state == Match) {
            matching = (alignment.sequence()[read_idx] == current_node->sequence()[node_idx]);
        }
        edit_length = 1;
        edit_read_end_idx = read_idx;
    }
    else if (matrix == Match &&
             (alignment.sequence()[read_idx] == current_node->sequence()[node_idx]) != matching) {
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
    cerr << "[BABuilder::update_state] finished updating state, matrix is " << (matrix_state == Match ? "match" : (matrix_state == InsertRow ? "insert column" : "insert row" )) << ", is matching? " << (matching ? "yes" : "no") << ", edit length " << edit_length << ", edit end index (on read) " << edit_read_end_idx << ", current node " << current_node->id() << endl;
#endif
}

void BandedGlobalAligner::BABuilder::finish_current_edit() {
    
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
}

void BandedGlobalAligner::BABuilder::finish_current_node() {

    // sentinel for first iteration
    if (current_node == nullptr) {
#ifdef debug_banded_aligner_traceback
        cerr << "[BABuilder::finish_current_node] at beginning of traceback, not creating a mapping" << endl;
#endif
        return;
    }
    
    finish_current_edit();
    
#ifdef debug_banded_aligner_traceback
    cerr << "[BABuilder::finish_current_node] finishing mapping for node " << current_node->id() << endl;
#endif
    
    node_mappings.emplace_front();
    for (Edit edit : mapping_edits) {
        *(node_mappings.front().add_edit()) = edit;
    }
    mapping_edits.clear();
    
    (*(node_mappings.front().mutable_position())).set_node_id(current_node->id());
    // note: global alignment always starts at beginning of node, default offset 0 is correct
}

void BandedGlobalAligner::BABuilder::finalize_alignment() {
    
    finish_current_node();
    
#ifdef debug_banded_aligner_traceback
    cerr << "[BABuilder::finalize_alignment] finalizing alignment" << endl;
#endif
    
    alignment.clear_path();
    Path* path = alignment.mutable_path();
    
    int32_t mapping_rank = 1;
    for (Mapping& mapping : node_mappings) {
        mapping.set_rank(mapping_rank);
        *(path->add_mapping()) = mapping;
        mapping_rank++;
    }
    
    node_mappings.clear();
    
#ifdef debug_banded_aligner_traceback
    cerr << "[BABuilder::finalize_alignment] alignment: " << pb2json(alignment) << endl;
#endif
}

BandedGlobalAligner::BAMatrix::BAMatrix(Alignment& alignment, Node* node, int64_t top_diag,
                                        int64_t bottom_diag, BAMatrix** seeds, int64_t num_seeds,
                                        int64_t cumulative_seq_len) :
                                        node(node),
                                        top_diag(top_diag),
                                        bottom_diag(bottom_diag),
                                        seeds(seeds),
                                        alignment(alignment),
                                        num_seeds(num_seeds),
                                        cumulative_seq_len(cumulative_seq_len),
                                        match(nullptr),
                                        insert_col(nullptr),
                                        insert_row(nullptr)
{
    // nothing to do
#ifdef debug_banded_aligner_objects
    cerr << "[BAMatrix]: constructor for node " << node->id() << " with sequence " << node->sequence() << " and band from " << top_diag << " to " << bottom_diag << endl;;
#endif
}

BandedGlobalAligner::BAMatrix::~BAMatrix() {
#ifdef debug_banded_aligner_objects
    if (node != nullptr) {
        cerr << "[BAMatrix::~BAMatrix] destructing matrix for node " << node->id() << endl;
    }
    else {
        cerr << "[BAMatrix::~BAMatrix] destructing null matrix" << endl;
    }
#endif
    free(match);
    free(insert_row);
    free(insert_col);
    free(seeds);
}

void BandedGlobalAligner::BAMatrix::fill_matrix(int8_t* score_mat, int8_t* nt_table, int8_t gap_open,
                                                int8_t gap_extend, bool qual_adjusted, int8_t min_inf) {
    
#ifdef debug_banded_aligner_fill_matrix
    cerr << "[BAMatrix::fill_matrix] beginning DP on matrix for node " << node->id() << endl;;
#endif
    
    // note: bottom has the higher index
    int64_t band_height = bottom_diag - top_diag + 1;
    int64_t ncols = node->sequence().length();
    int64_t band_size = band_height * ncols;
    
    const string& node_seq = node->sequence();
    const string& read = alignment.sequence();
    const string& base_quality = alignment.quality();
    
    match = (int8_t*) malloc(sizeof(int8_t) * band_size);
    insert_col = (int8_t*) malloc(sizeof(int8_t) * band_size);
    insert_row = (int8_t*) malloc(sizeof(int8_t) * band_size);
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
    
    int64_t idx, up_idx, diag_idx, left_idx;
    
    // POA into left hand column from the seeds
    if (seeds == nullptr) {
        if (cumulative_seq_len != 0) {
            cerr << "error:[BandedGlobalAligner] banded alignment has no node predecessor for node in middle of path" << endl;
            assert(0);
        }
        
#ifdef debug_banded_aligner_fill_matrix
        cerr << "[BAMatrix::fill_matrix]: no previous seeds, computing implied edge conditions" << endl;
#endif
        
        // first node in alignment, fill out initial column based on implied lead gaps
        
        // find position of the first cell in the rectangularized band
        int64_t iter_start = -top_diag;
        idx = iter_start * ncols;
        
        // cap stop index if last diagonal is below bottom of matrix
        int64_t iter_stop = bottom_diag > (int64_t) read.length() ? band_height + (int64_t) read.length() - bottom_diag - 1 : band_height;
        

        
        // match of first nucleotides
        if (qual_adjusted) {
            match[idx] = score_mat[25 * base_quality[0] + 5 * nt_table[node_seq[0]] + nt_table[read[0]]];
#ifdef debug_banded_aligner_fill_matrix
            cerr << "[BAMatrix::fill_matrix]: set quality adjusted initial match cell to " << (int) match[idx] << " from node char " << node_seq[0] << ", read char " << read[0] << ", base qual " << (int) base_quality[0] << " for adjusted matrix index " << 25 * base_quality[0] + 5 * nt_table[node_seq[0]] + nt_table[read[0]] << " and score " << (int) score_mat[25 * base_quality[0] + 5 * nt_table[node_seq[0]] + nt_table[read[0]]] << endl;
#endif
        }
        else {
            match[idx] = score_mat[5 * nt_table[node_seq[0]] + nt_table[read[0]]];
        }
        

        
        // only way to end an alignment in a gap here is to row and column gap
        insert_col[idx] = -2 * gap_open;
        insert_row[idx] = -2 * gap_open;
        
        for (int64_t i = iter_start + 1; i < iter_stop; i++) {
            idx = i * ncols;
            up_idx = (i - 1) * ncols;
            // score of a match in this cell
            int8_t match_score;
            if (qual_adjusted) {
                match_score = score_mat[25 * base_quality[top_diag + i] + 5 * nt_table[node_seq[0]] + nt_table[read[top_diag + i]]];
            }
            else {
                match_score = score_mat[5 * nt_table[node_seq[0]] + nt_table[read[top_diag + i]]];
            }
            // must take one lead gap to get into first column
            match[idx] = match_score - gap_open - (top_diag + i - 1) * gap_extend;
            // normal iteration along column
            insert_row[idx] = max(max(match[up_idx] - gap_open, insert_row[up_idx] - gap_extend),
                                  insert_col[up_idx] - gap_open);
            // must take two gaps to get into first column
            insert_col[idx] = -2 * gap_open - (top_diag + i) * gap_extend;

#ifdef debug_banded_aligner_fill_matrix
            cerr << "[BAMatrix::fill_matrix]: on left edge of matrix at rectangle coords (" << i << ", " << 0 << "), match score of node char " << 0 << " (" << node_seq[0] << ") and read char " << i + top_diag << " (" << read[i + top_diag] << ") is " << (int) match_score << ", leading gap length is " << top_diag + i << " for total match matrix score of " << (int) match[idx] << endl;
#endif
        }
        
        // fix the final insert column value, which actually would have been inserting from outside the band
        insert_col[idx] = min_inf;
    }
    else {
        
#ifdef debug_banded_aligner_fill_matrix
        cerr << "[BAMatrix::fill_matrix] seed nodes are present, preparing to do POA" << endl;;
#endif
        
        // are we clipping any diagonals because they are outside the range of the matrix in this column?
        bool bottom_diag_outside = bottom_diag >= (int64_t) read.length();
        bool top_diag_outside = top_diag < 0;
        bool top_diag_abutting = top_diag == 0;
        
        // rows in the rectangularized band corresponding to diagonals
        int64_t iter_start = top_diag_outside ? -top_diag : 0;
        int64_t iter_stop = bottom_diag_outside ? band_height + (int64_t) read.length() - bottom_diag - 1 : band_height;
        
        // handle special logic for lead column insertions
        idx = iter_start * ncols;
        
        int8_t match_score;
        if (qual_adjusted) {
            match_score = score_mat[25 * base_quality[iter_start + top_diag] + 5 * nt_table[node_seq[0]] + nt_table[read[iter_start + top_diag]]];
        }
        else {
            match_score = score_mat[5 * nt_table[node_seq[0]] + nt_table[read[iter_start + top_diag]]];
        }
        if (top_diag_outside || top_diag_abutting) {
#ifdef debug_banded_aligner_fill_matrix
            cerr << "[BAMatrix::fill_matrix]: node still at top of matrix, adding implied gap to match of length " << cumulative_seq_len << "and initializing rest of column to -inf" << endl;
#endif
            // match after implied gap along top edge
            match[idx] = match_score - gap_open - (cumulative_seq_len - 1) * gap_extend;
        }
        else {
#ifdef debug_banded_aligner_fill_matrix
            cerr << "[BAMatrix::fill_matrix]: no lead gap possible, initializing first column to all -inf" << endl;
#endif
            // no lead gaps possible so seed with identity of max function to prepare for POA
            match[idx] = min_inf;
        }
        
        if (top_diag_outside) {
            // gap open after implied gap along top edge
            insert_row[idx] = -2 * gap_open - cumulative_seq_len * gap_extend;
        }
        else {
            // no lead gaps possible so seed with identity of max function to prepare for POA
            insert_row[idx] = min_inf;
        }
        
        insert_col[idx] = min_inf;
        
        // start with each cell in the leftmost column as identity of max function to prepare for POA
        for (int64_t i = iter_start + 1; i < iter_stop; i++) {
            idx = i * ncols;
            match[idx] = min_inf;
            insert_col[idx] = min_inf;
        }
        
        // the band size selection algorithm should guarantee that the band extension of each seed
        // is present in the current matrix, but not that every cell in the current matrix has a
        // predecessor from every seed
        
        // iterate through each seed and carry forward its subset of the current band
        for (int64_t seed_num = 0; seed_num < num_seeds; seed_num++) {
            
            BAMatrix* seed = seeds[seed_num];
            
            if (seed == nullptr) {
#ifdef debug_banded_aligner_fill_matrix
                cerr << "[BAMatrix::fill_matrix]: seed is masked, skipping" << endl;
#endif
                continue;
            }
            
#ifdef debug_banded_aligner_fill_matrix
            cerr << "[BAMatrix::fill_matrix]: doing POA across boundary from " << seed_num << "-th seed node " << seed->node->id() << endl;
#endif
            
            int64_t seed_node_seq_len = seed->node->sequence().length();
            
            int64_t seed_next_top_diag = seed->top_diag + seed_node_seq_len;
            int64_t seed_next_bottom_diag = seed->bottom_diag + seed_node_seq_len;
            
            bool beyond_top_of_matrix = seed_next_top_diag < 0;
            bool beyond_bottom_of_matrix = seed_next_bottom_diag >= (int64_t) read.length();
            int64_t seed_next_top_diag_iter = beyond_top_of_matrix ? 0 : seed_next_top_diag;
            int64_t seed_next_bottom_diag_iter = beyond_bottom_of_matrix ? (int64_t) read.length() - 1 : seed_next_bottom_diag;
            
#ifdef debug_banded_aligner_fill_matrix
            cerr << "[BAMatrix::fill_matrix]: this seed reaches diagonals " << seed_next_top_diag << " to " << seed_next_bottom_diag << " out of matrix range " << top_diag << " to " << bottom_diag << endl;
#endif
            // special logic for first row
            idx = (seed_next_top_diag_iter - top_diag) * ncols;
            
            // may not be able to extend a match if at top of matrix
            if (!beyond_top_of_matrix) {
                diag_idx = (seed_next_top_diag_iter - seed_next_top_diag) * seed_node_seq_len + seed_node_seq_len - 1;
                if (qual_adjusted) {
                    match_score = score_mat[25 * base_quality[seed_next_top_diag_iter] + 5 * nt_table[node_seq[0]] + nt_table[read[seed_next_top_diag_iter]]];
                }
                else {
                    match_score = score_mat[5 * nt_table[node_seq[0]] + nt_table[read[seed_next_top_diag_iter]]];
                }
                match[idx] = max<int8_t>(match_score + max<int8_t>(max<int8_t>(seed->match[diag_idx],
                                                                               seed->insert_row[diag_idx]),
                                                                   seed->insert_col[diag_idx]), match[idx]);
            }
            
            if (seed_next_top_diag_iter != seed_next_bottom_diag_iter) {
                left_idx = (seed_next_top_diag_iter - seed_next_top_diag + 2) * seed_node_seq_len - 1;
                insert_col[idx] = max<int8_t>(max<int8_t>(max<int8_t>(seed->match[left_idx] - gap_open,
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
                int8_t match_score;
                if (qual_adjusted) {
                    match_score = score_mat[25 * base_quality[diag] + 5 * nt_table[node_seq[0]] + nt_table[read[diag]]];
                }
                else {
                    match_score = score_mat[5 * nt_table[node_seq[0]] + nt_table[read[diag]]];
                }
                
#ifdef debug_banded_aligner_fill_matrix
                cerr << "[BAMatrix::fill_matrix]: extending match from rectangular coord (" << diag - seed_next_top_diag << ", " << seed_node_seq_len - 1 << ")" << " with match score " << (int) match_score << ", scores are " << (int) seed->match[diag_idx] << " (M), " << (int) seed->insert_row[diag_idx] << " (Ir), and " << (int) seed->insert_col[diag_idx] << " (Ic), current score is " << (int) match[idx] << endl;
#endif
                
                match[idx] = max<int8_t>(match_score + max<int8_t>(max<int8_t>(seed->match[diag_idx],
                                                                               seed->insert_row[diag_idx]),
                                                                   seed->insert_col[diag_idx]), match[idx]);
                
                
                // extend a column gap
                left_idx = (diag - seed_next_top_diag + 2) * seed_node_seq_len - 1;
                
#ifdef debug_banded_aligner_fill_matrix
                cerr << "[BAMatrix::fill_matrix]: extending match from rectangular coord (" << diag - seed_next_top_diag + 1 << ", " << seed_node_seq_len - 1 << ")" << ", scores are " << (int) seed->match[left_idx] << " (M), " << (int) seed->insert_row[left_idx] << " (Ir), and " << (int) seed->insert_col[left_idx] << " (Ic), current score is " << (int) insert_col[idx] << endl;
#endif
                insert_col[idx] = max<int8_t>(max<int8_t>(max<int8_t>(seed->match[left_idx] - gap_open,
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
                match[idx] = max<int8_t>(match_score + max<int8_t>(max<int8_t>(seed->match[diag_idx],
                                                                               seed->insert_row[diag_idx]),
                                                                   seed->insert_col[diag_idx]), match[idx]);

                
                // can only extend column gap if the bottom of the matrix was hit in the last seed
                if (beyond_bottom_of_matrix) {
#ifdef debug_banded_aligner_fill_matrix
                    cerr << "[BAMatrix::fill_matrix]: can also extend a column gap since already reached edge of matrix" << endl;
#endif
                    left_idx = (seed_next_bottom_diag_iter - seed_next_top_diag + 2) * seed_node_seq_len - 1;
                    insert_col[idx] = max<int8_t>(max<int8_t>(max<int8_t>(seed->match[left_idx] - gap_open,
                                                                          seed->insert_row[left_idx] - gap_open),
                                                              seed->insert_col[left_idx] - gap_extend), insert_col[idx]);
                }
            }

        }
#ifdef debug_banded_aligner_fill_matrix
        cerr << "[BAMatrix::fill_matrix]: filling out initial column of insert row matrix" << endl;
#endif
        // compute the insert row scores (they can be computed after the POA iterations since they do not
        // cross node boundaries)
        for (int64_t i = iter_start + 1; i < iter_stop; i++) {
            idx = i * ncols;
            up_idx = (i - 1) * ncols;
            
            insert_row[idx] = max(max(match[up_idx] - gap_open, insert_row[up_idx] - gap_extend),
                                  insert_col[up_idx] - gap_open);
        }
    }
    
#ifdef debug_banded_aligner_fill_matrix
    cerr << "[BAMatrix::fill_matrix]: seeding finished, moving to subsequent columns" << endl;
#endif
    
    // iterate through the rest of the columns
    for (int64_t j = 1; j < ncols; j++) {
        
        // are we clipping any diagonals because they are outside the range of the matrix in this column?
        bool bottom_diag_outside = bottom_diag + j >= (int64_t) read.length();
        bool top_diag_outside = top_diag + j < 0;
        bool top_diag_abutting = top_diag + j == 0;
        
        int64_t iter_start = top_diag_outside ? -(top_diag + j) : 0;
        int64_t iter_stop = bottom_diag_outside ? band_height + (int64_t) read.length() - bottom_diag - j - 1 : band_height;
        
        idx = iter_start * ncols + j;
        
        int8_t match_score;
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
    print_full_matrices();
    print_rectangularized_bands();
#endif
}

void BandedGlobalAligner::BAMatrix::traceback(BABuilder& builder, AltTracebackStack& traceback_stack, matrix_t start_mat,
                                              int8_t* score_mat, int8_t* nt_table, int8_t gap_open, int8_t gap_extend,
                                              bool qual_adjusted, int8_t min_inf) {
    
    // get coordinates of bottom right corner
    const string& read = alignment.sequence();
    int64_t ncols = node->sequence().length();
    int64_t row = bottom_diag + ncols > (int64_t) read.length() ? (int64_t) read.length() - top_diag - ncols : bottom_diag - top_diag;
    int64_t col = ncols - 1;
    int64_t idx = row * ncols + col;
    
    
#ifdef debug_banded_aligner_traceback
    cerr << "[BAMatrix::traceback] beginning traceback in matrices for node " << node->id() << " starting matrix is " << (start_mat == Match ? "match" : (start_mat == InsertCol ? "insert column" : "insert row")) << endl;
#endif
    
    traceback_internal(builder, traceback_stack, row, col, start_mat, false, score_mat, nt_table, gap_open, gap_extend,
                       qual_adjusted, min_inf);
}

void BandedGlobalAligner::BAMatrix::traceback_internal(BABuilder& builder, AltTracebackStack& traceback_stack, int64_t start_row,
                                                       int64_t start_col, matrix_t start_mat, bool in_lead_gap, int8_t* score_mat,
                                                       int8_t* nt_table, int8_t gap_open, int8_t gap_extend, bool qual_adjusted,
                                                       int8_t min_inf) {
    
#ifdef debug_banded_aligner_traceback
    cerr << "[BAMatrix::traceback_internal] starting traceback back through node " << node->id() << " from rectangular coordinates (" << start_row << ", " << start_col << "), currently " << (in_lead_gap ? "" : "not ") << "in a lead gap" << endl;
#endif
    
    //cerr << "node " << node->id() << endl;
    
    const string& read = alignment.sequence();
    const string& base_quality = alignment.quality();
    
    int64_t band_height = bottom_diag - top_diag + 1;
    const char* node_seq = node->sequence().c_str();
    int64_t ncols = node->sequence().length();
    int64_t node_id = node->id();
    
    int64_t idx, next_idx;
    int64_t i = start_row, j = start_col;
    matrix_t curr_mat = start_mat;
    int8_t curr_score;
    int8_t source_score;
    int8_t score_diff;
    int8_t alt_score;
    int8_t curr_traceback_score = traceback_stack.current_traceback_score();
    
    // do node traceback unless we are in the lead gap implied at the edge of the DP matrix or we
    // are already at a node boundary trying to get across
    while (j > 0 || curr_mat == InsertRow) {
        if (in_lead_gap) {
            break;
        }
        
#ifdef debug_banded_aligner_traceback
        cerr << "[BAMatrix::traceback_internal] traceback coordinates (" << i << ", " << j << "), current matrix is " << (curr_mat == Match ? "match" : (curr_mat == InsertCol ? "insert column" : "insert row")) << endl;
#endif
        
        // add traceback step to alignment
        builder.update_state(curr_mat, node, i + top_diag + j, j);
        
        // check for a deflection
        if (traceback_stack.at_next_deflection(node_id, i, j)) {
            
            // move to the next position as dictated by current matrix
            switch (curr_mat) {
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
            curr_mat = traceback_stack.deflect_to_matrix();
            
#ifdef debug_banded_aligner_traceback
            cerr << "[BAMatrix::traceback_internal] taking inside matrix deflection to " << (curr_mat == Match ? "match" : (curr_mat == InsertCol ? "insert column" : "insert row")) << endl;
#endif
            
            continue;
        }
        
        // find optimal traceback
        idx = i * ncols + j;
        bool found_trace = false;
        switch (curr_mat) {
            case Match:
            {
                if (i + j == -top_diag) {
#ifdef debug_banded_aligner_traceback
                    cerr << "[BAMatrix::traceback_internal] next cell is outside matrix, opening implied lead gap" << endl;
#endif
                    // at top of matrix, move into implied lead gap along top edge
                    j--;
                    in_lead_gap = true;
                    // no where else to go, so break out of switch statement without checking alts
                    break;
                }
                
                curr_score = match[idx];
                next_idx = i * ncols + j - 1;
                
                int8_t match_score;
                if (qual_adjusted) {
                    match_score = score_mat[25 * base_quality[i + top_diag + j] + 5 * nt_table[node_seq[j]] + nt_table[read[i + top_diag + j]]];
                }
                else {
                    match_score = score_mat[5 * nt_table[node_seq[j]] + nt_table[read[i + top_diag + j]]];
                }
                
#ifdef debug_banded_aligner_traceback
                cerr << "[BAMatrix::traceback_internal] transitioning from match, current score " << (int) match[idx] << " match/mismatch score " << (int) match_score << " from node char " << j << " (" << node_seq[j] << ") and read char " << i + top_diag + j << " (" << read[i + top_diag + j] << ")" << endl;
#endif
                
                source_score = match[next_idx];
                score_diff = curr_score - (source_score + match_score);
                if (score_diff == 0) {
#ifdef debug_banded_aligner_traceback
                    cerr << "[BAMatrix::traceback_internal] found next cell in match matrix with score " << (int) match[next_idx] << endl;
#endif
                    curr_mat = Match;
                    found_trace = true;
                }
                else {
                    alt_score = curr_traceback_score - score_diff;
                    traceback_stack.propose_deflection(alt_score, node_id, i, j, node_id, Match);
                }
                
                source_score = insert_row[next_idx];
                if (source_score > min_inf) {
                    score_diff = curr_score - (source_score + match_score);
                    if (!found_trace && score_diff == 0) {
#ifdef debug_banded_aligner_traceback
                        cerr << "[BAMatrix::traceback_internal] found next cell in insert row matrix with score " << (int) insert_row[next_idx] << endl;
#endif
                        curr_mat = InsertRow;
                        found_trace = true;
                    }
                    else {
                        alt_score = curr_traceback_score - score_diff;
                        traceback_stack.propose_deflection(alt_score, node_id, i, j, node_id, InsertRow);
                    }
                }
                
                source_score = insert_col[next_idx];
                if (source_score > min_inf) {
                    score_diff = curr_score - (source_score + match_score);
                    if (!found_trace && score_diff == 0) {
#ifdef debug_banded_aligner_traceback
                        cerr << "[BAMatrix::traceback_internal] found next cell in insert column matrix with score " << (int) insert_col[next_idx] << endl;
#endif
                        curr_mat = InsertCol;
                        found_trace = true;
                    }
                    else {
                        alt_score = curr_traceback_score - score_diff;
                        traceback_stack.propose_deflection(alt_score, node_id, i, j, node_id, InsertCol);
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
                    break;
                }
                
                curr_score = insert_row[idx];
                next_idx = (i - 1) * ncols + j;
                
                source_score = match[next_idx];
                score_diff = curr_score - (source_score - gap_open);
                if (score_diff == 0) {
                    curr_mat = Match;
                    found_trace = true;
                }
                else {
                    alt_score = curr_traceback_score - score_diff;
                    traceback_stack.propose_deflection(alt_score, node_id, i, j, node_id, Match);
                }
                
                source_score = insert_row[next_idx];
                if (source_score > min_inf) {
                    score_diff = curr_score - (source_score - gap_extend);
                    if (!found_trace && score_diff == 0) {
                        curr_mat = InsertRow;
                        found_trace = true;
                    }
                    else {
                        alt_score = curr_traceback_score - score_diff;
                        traceback_stack.propose_deflection(alt_score, node_id, i, j, node_id, InsertRow);
                    }
                }
                
                source_score = insert_col[next_idx];
                if (source_score > min_inf) {
                    score_diff = curr_score - (source_score - gap_open);
                    if (!found_trace && score_diff == 0) {
                        curr_mat = InsertCol;
                        found_trace = true;
                    }
                    else {
                        alt_score = curr_traceback_score - score_diff;
                        traceback_stack.propose_deflection(alt_score, node_id, i, j, node_id, InsertCol);
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

                source_score = match[next_idx];
                score_diff = curr_score - (source_score - gap_open);
                if (score_diff == 0) {
                    curr_mat = Match;
                    found_trace = true;
                }
                else {
                    alt_score = curr_traceback_score - score_diff;
                    traceback_stack.propose_deflection(alt_score, node_id, i, j, node_id, Match);
                }
                
                source_score = insert_row[next_idx];
                if (source_score > min_inf) {
                    score_diff = curr_score - (source_score - gap_open);
                    if (!found_trace && score_diff == 0) {
                        curr_mat = InsertRow;
                        found_trace = true;
                    }
                    else {
                        alt_score = curr_traceback_score - score_diff;
                        traceback_stack.propose_deflection(alt_score, node_id, i, j, node_id, InsertRow);
                    }
                }
                
                source_score = insert_col[next_idx];
                if (source_score > min_inf) {
                    score_diff = curr_score - (source_score - gap_extend);
                    if (!found_trace && score_diff == 0) {
                        curr_mat = InsertCol;
                        found_trace = true;
                    }
                    else {
                        alt_score = curr_traceback_score - score_diff;
                        traceback_stack.propose_deflection(alt_score, node_id, i, j, node_id, InsertCol);
                    }
                }
                
                if (!found_trace) {
                    cerr << "error:[BandedGlobalAligner] traceback stuck in insert row matrix interior" << endl;
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
        // add lead column gaps until reaching edge of node
        curr_mat = InsertCol;
        while (j > 0) {
            builder.update_state(curr_mat, node, -1, j);
            j--;
            i++;
        }
    }
    
    // are there are predecessor nodes?
    if (num_seeds > 0) {
        
        BAMatrix* traceback_seed;
        int64_t traceback_seed_row;
        int64_t traceback_seed_col;
        
        if (in_lead_gap) {
#ifdef debug_banded_aligner_traceback
            cerr << "[BAMatrix::traceback_internal] at boundary, following seed backward from a lead gap" << endl;
#endif
            // we are in the implied lead gap along the top of the matrix
            // take the shortest path back to the origin of the global alignment
            
            // add final read deletion of node
            builder.update_state(curr_mat, node, -1, j);
            
            int64_t shortest_seq_len = numeric_limits<int64_t>::max();
            for (int64_t k = 0; k < num_seeds; k++) {
                BAMatrix* seed = seeds[k];
                // is seed masked?
                if (seed == nullptr) {
                    continue;
                }
                // does a shorter path to origin go through this seed?
                if (seed->cumulative_seq_len < shortest_seq_len) {
                    traceback_seed = seed;
                    shortest_seq_len = seed->cumulative_seq_len;
                }
            }
            
            // where in the matrix is this?
            int64_t seed_ncols = traceback_seed->node->sequence().length();
            int64_t seed_extended_top_diag = traceback_seed->top_diag + seed_ncols;
            traceback_seed_row = top_diag - seed_extended_top_diag + i + 1;
            traceback_seed_col = seed_ncols - 1;
        }
        else {
            
            int64_t curr_diag = top_diag + i;
            builder.update_state(curr_mat, node, curr_diag, 0);
            
            if (traceback_stack.at_next_deflection(node_id, i, j)) {
#ifdef debug_banded_aligner_traceback
                cerr << "[BAMatrix::traceback_internal] at boundary, taking a deflection" << endl;
#endif
                
                // where to deflect to?
                int64_t deflect_node_id;
                matrix_t deflect_matrix = traceback_stack.deflect_to_matrix(deflect_node_id);
                

                
                // find which seed matrix to deflect to (don't have a better way of looking this up right now)
                BAMatrix* seed;
                for (int64_t k = 0; k < num_seeds; k++) {
                    seed = seeds[k];
                    if (seed->node->id() == deflect_node_id) {
                        break;
                    }
                }
                
                int64_t seed_ncols = seed->node->sequence().length();
                traceback_seed_row = curr_diag - seed->top_diag - seed_ncols + (curr_mat == InsertCol);
                traceback_seed_col = seed_ncols - 1;
                
#ifdef debug_banded_aligner_traceback
                cerr << "[BAMatrix::traceback_internal] taking node boundary deflection to " << (deflect_matrix == Match ? "match" : (deflect_matrix == InsertCol ? "insert column" : "insert row")) << " in node " << deflect_node_id << ", will start at coordinates (" << traceback_seed_row << ", " << traceback_seed_col << ")" << endl;
#endif
                
                // continue traceback in the next node
                seed->traceback_internal(builder, traceback_stack, traceback_seed_row, traceback_seed_col, deflect_matrix,
                                         in_lead_gap, score_mat, nt_table, gap_open, gap_extend, qual_adjusted, min_inf);
                return;
            }
            
#ifdef debug_banded_aligner_traceback
            cerr << "[BAMatrix::traceback_internal] at boundary, looking for predecessor node" << endl;
#endif
            
            int8_t match_score;
            switch (curr_mat) {
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
            cerr << "[BAMatrix::traceback_internal] at boundary, following seed backward from " << (curr_mat == Match ? "match" : "insert column") << " matrix with score " << (int) curr_score << endl;
#endif
            
            // matches stay on same diagonal, column insertions move over one diagonal
            // note that the indexing is relative to THIS matrix, not the seed
            
            bool found_trace = false;
            // check traceback goes to each seed matrix
            for (int64_t k = 0; k < num_seeds; k++) {
                BAMatrix* seed = seeds[k];

                // is the matrix masked?
                if (seed == nullptr) {
#ifdef debug_banded_aligner_traceback
                    cerr << "[BAMatrix::traceback_internal] seed " << k << " is masked" << endl;
#endif
                    continue;
                }
                
#ifdef debug_banded_aligner_traceback
                cerr << "[BAMatrix::traceback_internal] checking seed number " << k << " which corresponds to node " << seed->node->id() << endl;
#endif
                
                int64_t seed_node_id = seed->node->id();
                int64_t seed_ncols = seed->node->sequence().length();
                
                // the diagonals in the current matrix that this seed extends to
                int64_t seed_extended_top_diag = seed->top_diag + seed_ncols;
                int64_t seed_extended_bottom_diag = seed->bottom_diag + seed_ncols;
                
                // does the traceback diagonal extend backward to this matrix?
                if (curr_diag > seed_extended_bottom_diag - (curr_mat == InsertCol) // col inserts hit 1 less of band
                    || curr_diag < seed_extended_top_diag) {
                    continue;
                }
                
                int64_t seed_col = seed_ncols - 1;
                int64_t seed_row = curr_diag - seed_extended_top_diag + (curr_mat == InsertCol);
                next_idx = seed_row * seed_ncols + seed_col;
                
#ifdef debug_banded_aligner_traceback
                cerr << "[BAMatrix::traceback_internal] checking seed rectangular coordinates (" << seed_row << ", " << seed_col << ")" << endl;
#endif
                
                switch (curr_mat) {
                    case Match:
                    {
                        source_score = seed->match[next_idx];
                        score_diff = curr_score - (source_score + match_score);
                        if (score_diff == 0 && !found_trace) {
                            curr_mat = Match;
                            traceback_seed = seed;
                            traceback_seed_row = seed_row;
                            traceback_seed_col = seed_col;
                            found_trace = true;
#ifdef debug_banded_aligner_traceback
                            cerr << "[BAMatrix::traceback_internal] hit found in match matrix with score " << (int) seed->match[next_idx] << endl;
#endif
                        }
                        else {
                            alt_score = curr_traceback_score - score_diff;
                            traceback_stack.propose_deflection(alt_score, node_id, i, j, seed_node_id, Match);
                        }
                        
                        source_score = seed->insert_col[next_idx];
                        if (source_score > min_inf) {
                            score_diff = curr_score - (source_score + match_score);
                            if (score_diff == 0 && !found_trace) {
#ifdef debug_banded_aligner_traceback
                                cerr << "[BAMatrix::traceback_internal] hit found in insert column matrix  with score " << (int) seed->insert_col[next_idx] << endl;
#endif
                                curr_mat = InsertCol;
                                traceback_seed = seed;
                                traceback_seed_row = seed_row;
                                traceback_seed_col = seed_col;
                                found_trace = true;
                            }
                            else {
                                alt_score = curr_traceback_score - score_diff;
                                traceback_stack.propose_deflection(alt_score, node_id, i, j, seed_node_id, InsertCol);
                            }
                        }
                        
                        source_score = seed->insert_row[next_idx];
                        if (source_score > min_inf) {
                            score_diff = curr_score - (source_score + match_score);
                            if (score_diff == 0 && !found_trace) {
#ifdef debug_banded_aligner_traceback
                                cerr << "[BAMatrix::traceback_internal] hit found in insert row matrix  with score " << (int) seed->insert_row[next_idx] << endl;
#endif
                                curr_mat = InsertRow;
                                traceback_seed = seed;
                                traceback_seed_row = seed_row;
                                traceback_seed_col = seed_col;
                                found_trace = true;
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
                        source_score = seed->match[next_idx];
                        score_diff = curr_score - (source_score - gap_open);
                        if (score_diff == 0 && !found_trace) {
                            curr_mat = Match;
                            traceback_seed = seed;
                            traceback_seed_row = seed_row;
                            traceback_seed_col = seed_col;
                            found_trace = true;
                        }
                        else {
                            alt_score = curr_traceback_score - score_diff;
                            traceback_stack.propose_deflection(alt_score, node_id, i, j, seed_node_id, Match);
                        }
                        
                        source_score = seed->insert_col[next_idx];
                        if (source_score > min_inf) {
                            score_diff = curr_score - (source_score - gap_extend);
                            if (score_diff == 0 && !found_trace) {
                                curr_mat = InsertCol;
                                traceback_seed = seed;
                                traceback_seed_row = seed_row;
                                traceback_seed_col = seed_col;
                                found_trace = true;
                            }
                            else {
                                alt_score = curr_traceback_score - score_diff;
                                traceback_stack.propose_deflection(alt_score, node_id, i, j, seed_node_id, InsertCol);
                            }
                        }
                        
                        source_score = seed->insert_row[next_idx];
                        if (source_score > min_inf) {
                            score_diff = curr_score - (source_score - gap_open);
                            if (score_diff == 0 && !found_trace) {
                                curr_mat = InsertRow;
                                traceback_seed = seed;
                                traceback_seed_row = seed_row;
                                traceback_seed_col = seed_col;
                                found_trace = true;
                            }
                            else {
                                alt_score = curr_traceback_score - score_diff;
                                traceback_stack.propose_deflection(alt_score, node_id, i, j, seed_node_id, InsertRow);
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
            
            if (!found_trace) {
                cerr << "error:[BandedGlobalAligner] traceback stuck at node boundary" << endl;
                assert(0);
            }
        }
        
        // continue traceback in the next node
        traceback_seed->traceback_internal(builder, traceback_stack, traceback_seed_row, traceback_seed_col, curr_mat,
                                           in_lead_gap, score_mat, nt_table, gap_open, gap_extend, qual_adjusted, min_inf);
    }
    else {
        // this is the first node in the alignment, finish off the alignment and return
#ifdef debug_banded_aligner_traceback
        cerr << "[BAMatrix::traceback_internal] at beginning of first node in alignment" << endl;
#endif
        if (in_lead_gap) {
            builder.update_state(curr_mat, node, -1, 0);
        }
        else {
            
            builder.update_state(curr_mat, node, i + top_diag, 0);
            
            switch (curr_mat) {
                case Match:
                {
#ifdef debug_banded_aligner_traceback
                    cerr << "[BAMatrix::traceback_internal] alignment starts with match, adding read char " << top_diag + i << ": " << read[top_diag + i] << endl;
#endif
                    j--;
                    break;
                }
                    
                case InsertCol:
                {
#ifdef debug_banded_aligner_traceback
                    cerr << "[BAMatrix::traceback_internal] alignment starts with column gap" << endl;
#endif
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
            
            // add any lead row gaps necessary
            curr_mat = InsertRow;
            while (top_diag + i > 0) {
#ifdef debug_banded_aligner_traceback
                cerr << "[BAMatrix::traceback_internal] initial row gaps are present, adding read char " << top_diag + i - 1 << ": " << (char) tolower(read[top_diag + i - 1]) << " (lowercase of " << read[top_diag + i - 1]     << ")" << endl;
#endif
                builder.update_state(curr_mat, node, i + top_diag - 1, -1);
                i--;
            }
        }
    }
}

void BandedGlobalAligner::BAMatrix::print_full_matrices() {
    if (match == nullptr) {
        cerr << "error:[BandedGlobalAligner] cannot print matrix before performing dynamic programming" << endl;
        assert(0);
    }
    
    cout << "matrices for node " << node->id() << ":" << endl;
    
    for (matrix_t mat : {Match, InsertRow, InsertCol}) {
        print_matrix(mat);
    }

}

void BandedGlobalAligner::BAMatrix::print_rectangularized_bands() {
    if (match == nullptr) {
        cerr << "error:[BandedGlobalAligner] cannot print band before performing dynamic programming" << endl;
        assert(0);
    }
    
    cout << "rectangularized bands for node " << node->id() << ":" << endl;
    
    for (matrix_t mat : {Match, InsertRow, InsertCol}) {
        print_band(mat);
    }
}

void BandedGlobalAligner::BAMatrix::print_matrix(matrix_t which_mat) {

    const string& read = alignment.sequence();
    const string& node_seq = node->sequence();
    
    int8_t* band_rect;
    switch (which_mat) {
        case Match:
            cout << "match:" << endl;
            band_rect = match;
            break;
            
        case InsertRow:
            cout << "insert row:" << endl;
            band_rect = insert_row;
            break;
            
        case InsertCol:
            cout << "insert column:" << endl;
            band_rect = insert_col;
            break;
            
        default:
            cerr << "error:[BandedGlobalAligner] unrecognized matrix type" << endl;
            assert(0);
            break;
    }
    
    for (auto iter = node_seq.begin(); iter != node_seq.end(); iter++) {
        cout << "\t" << *iter;
    }
    cout << endl;
    
    int64_t ncols = node_seq.length();
    
    for (int64_t i = 0; i < (int64_t) read.length(); i++) {
        cout << read[i];
        for (int64_t j = 0; j < ncols; j++) {
            int64_t diag = i - j;
            if (diag < top_diag || diag > bottom_diag) {
                cout << "\t.";
            }
            else {
                cout << "\t" << (int) band_rect[(diag - top_diag) * ncols + j];
            }
        }
        cout << endl;
    }
}

void BandedGlobalAligner::BAMatrix::print_band(matrix_t which_mat) {
    
    const string& read = alignment.sequence();
    const string& node_seq = node->sequence();
    
    int8_t* band_rect;
    switch (which_mat) {
        case Match:
            cout << "match:" << endl;
            band_rect = match;
            break;
            
        case InsertRow:
            cout << "insert row:" << endl;
            band_rect = insert_row;
            break;
            
        case InsertCol:
            cout << "insert column:" << endl;
            band_rect = insert_col;
            break;
            
        default:
            cerr << "error:[BandedGlobalAligner] unrecognized matrix type" << endl;
            assert(0);
            break;
    }
    
    for (auto iter = node_seq.begin(); iter != node_seq.end(); iter++) {
        cout << "\t" << *iter;
    }
    cout << endl;
    
    int64_t band_height = bottom_diag - top_diag + 1;
    int64_t ncols = node_seq.length();
    
    for (int64_t i = 0; i < band_height; i++) {
        if (top_diag + i > 0) {
            cout << read[top_diag + i - 1];
        }
        
        for (int64_t j = 0; j < ncols; j++) {
            if (i + j < -top_diag || top_diag + i + j >= (int64_t) read.length()) {
                cout << "\t.";
            }
            else {
                cout << "\t" << (int) band_rect[i * ncols + j];
            }
        }
        cout << endl;
    }
    
    
    for (int64_t i = bottom_diag; i < min(bottom_diag + ncols, (int64_t) read.length()); i++) {
        cout << read[i] << "\t";
    }
    cout << endl;
}

BandedGlobalAligner::BandedGlobalAligner(Alignment& alignment, Graph& g,
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

BandedGlobalAligner::BandedGlobalAligner(Alignment& alignment, Graph& g,
                                         vector<Alignment>& alt_alignments,
                                         int64_t max_multi_alns, int64_t band_padding,
                                         bool permissive_banding,
                                         bool adjust_for_base_quality) :
                                         BandedGlobalAligner(alignment, g,
                                                             &alt_alignments,
                                                             max_multi_alns,
                                                             permissive_banding,
                                                             adjust_for_base_quality)
{
    // check data integrity and funnel into internal constructor
    if (!alt_alignments.empty()) {
        cerr << "error:[BandedGlobalAligner] alternate alignment vector must be empty before aligning" << endl;
        assert(0);
    }
}


BandedGlobalAligner::BandedGlobalAligner(Alignment& alignment, Graph& g,
                                         vector<Alignment>* alt_alignments,
                                         int64_t max_multi_alns,
                                         int64_t band_padding,
                                         bool permissive_banding,
                                         bool adjust_for_base_quality) :
                                         alignment(alignment),
                                         alt_alignments(alt_alignments),
                                         max_multi_alns(max_multi_alns),
                                         adjust_for_base_quality(adjust_for_base_quality)
{
    if (adjust_for_base_quality) {
        if (alignment.quality().empty()) {
            cerr << "error:[BandedGlobalAligner] alignment needs base quality to perform quality adjusted alignment" << endl;
            assert(0);
        }
    }
    
    // TODO: this can waste memory, but reallocating the vector seems to throw an error in protobuf and
    // we won't know if there are fewer alignments than the max until the cycle is over
    if (alt_alignments) {
        alt_alignments->reserve(max_multi_alns);
    }
    
    // map node ids to indices
    for (int64_t i = 0; i < g.node_size(); i++) {
        node_id_to_idx[g.node(i).id()] = i;
    }
    
#ifdef debug_banded_aligner_objects
    cerr << "[BandedGlobalAligner]: constructing edge lists by node" << endl;
#endif
    
    // convert the graph into lists of ids indicating incoming or outgoing edges
    vector<vector<int64_t>> node_edges_in;
    vector<vector<int64_t>> node_edges_out;
    graph_edge_lists(g, true, node_edges_out);
    graph_edge_lists(g, false, node_edges_in);
    
#ifdef debug_banded_aligner_objects
    cerr << "[BandedGlobalAligner]: performing topological sort" << endl;
#endif
    
    // compute topological ordering
    topological_sort(g, node_edges_out, topological_order);
    
#ifdef debug_banded_aligner_objects
    cerr << "[BandedGlobalAligner]: identifying source and sink nodes" << endl;
#endif
    
    // identify source and sink nodes in the graph
    for (int64_t i = 0; i < g.node_size(); i++) {
        if (node_edges_in[i].empty()) {
            source_nodes.insert(g.mutable_node(i));
        }
        if (node_edges_out[i].empty()) {
            sink_nodes.insert(g.mutable_node(i));
        }
    }
    
    if (source_nodes.empty() || sink_nodes.empty()) {
        cerr << "error:[BandedGlobalAligner] alignment graph must be a DAG" << endl;
    }
    
#ifdef debug_banded_aligner_objects
    cerr << "[BandedGlobalAligner]: computing node bands" << endl;
#endif
    
    // figure out what the bands need to be for alignment and which nodes cannot complete a
    // global alignment within the band
    vector<bool> node_masked;
    vector<pair<int64_t, int64_t>> band_ends;
    find_banded_paths(alignment.sequence(), permissive_banding, node_edges_in, node_edges_out, band_padding, node_masked, band_ends);
    
#ifdef debug_banded_aligner_objects
    cerr << "[BandedGlobalAligner]: identifying shortest paths" << endl;
#endif
    
    // find the shortest sequence leading to each node so we can infer the length
    // of lead deletions
    vector<int64_t> shortest_seqs;
    shortest_seq_paths(node_edges_out, shortest_seqs);
    
#ifdef debug_banded_aligner_objects
    cerr << "[BandedGlobalAligner]: constructing banded matrix objects" << endl;
#endif
    
    // initialize DP matrices for each node
    banded_matrices.resize(g.node_size());
    for (int64_t i = 0; i < g.node_size(); i++) {
        
#ifdef debug_banded_aligner_objects
        cerr << "[BandedGlobalAligner]: creating matrix object for node " << g.mutable_node(i)->id() << " at index " << i << endl;
#endif
        
        if (node_masked[i]) {
#ifdef debug_banded_aligner_objects
            cerr << "[BandedGlobalAligner]: node is masked, creating dummy matrix object" << endl;
#endif
            
            banded_matrices[i] = nullptr;
        }
        else {
            Node* node = g.mutable_node(i);
            int64_t node_seq_len = node->sequence().length();
            
#ifdef debug_banded_aligner_objects
            cerr << "[BandedGlobalAligner]: establishing seed list for node " << node->id() << " at index " << i << endl;
#endif
            
            // POA predecessor matrices
            BAMatrix** seeds;
            vector<int64_t>& edges_in = node_edges_in[i];
            if (edges_in.empty()) {
                
#ifdef debug_banded_aligner_objects
                cerr << "[BandedGlobalAligner]: no seeds, setting array to null" << endl;
#endif
                
                seeds = nullptr;
            }
            else {
                seeds = (BAMatrix**) malloc(sizeof(BAMatrix**) * edges_in.size());
                for (int64_t j = 0; j < edges_in.size(); j++) {
                    // add sentinel for masked nodes
                    seeds[j] = node_masked[edges_in[j]] ? nullptr : banded_matrices[edges_in[j]];
                }
            }
            
            banded_matrices[i] = new BAMatrix(alignment,
                                                           node,
                                                           band_ends[i].first,
                                                           band_ends[i].second,
                                                           seeds,
                                                           edges_in.size(),
                                                           shortest_seqs[i]);
            
        }
    }
    
    if (!permissive_banding) {
        bool sinks_masked = true;
        for (Node* node : sink_nodes) {
            if (banded_matrices[node_id_to_idx[node->id()]] != nullptr) {
                sinks_masked = false;
                break;
            }
            
        }
        if (sinks_masked) {
            cerr << "error:[BandedGlobalAligner] cannot align to graph within band, consider permissive banding" << endl;
            cerr << alignment.sequence() << endl;
            assert(0);
        }
    }
}

BandedGlobalAligner::~BandedGlobalAligner() {
    
    for (BAMatrix* banded_matrix : banded_matrices) {
        delete banded_matrix;
    }
    
}

// fills a vector with vectors ids that have edges to/from each node
void BandedGlobalAligner::graph_edge_lists(Graph& g, bool outgoing_edges, vector<vector<int64_t>>& out_edge_list) {
    if (outgoing_edges) {
        out_edge_list = vector<vector<int64_t>>(g.node_size());
        for (int64_t i = 0; i < g.edge_size(); i++) {
            Edge* edge = g.mutable_edge(i);
            out_edge_list[node_id_to_idx.at(edge->from())].push_back(node_id_to_idx.at(edge->to()));
        }
    }
    else {
        out_edge_list = vector<vector<int64_t>>(g.node_size());
        for (int64_t i = 0; i < g.edge_size(); i++) {
            Edge* edge = g.mutable_edge(i);
            out_edge_list[node_id_to_idx.at(edge->to())].push_back(node_id_to_idx.at(edge->from()));
        }
    }
}

// standard DFS-based topological sort algorithm
// NOTE: this is only valid if the Graph g has been dag-ified first and there are from_start
// or to_end edges.
void BandedGlobalAligner::topological_sort(Graph& g, vector<vector<int64_t>>& node_edges_out,
                                           vector<Node*>& out_topological_order) {
    if (g.node_size() == 0) {
        cerr << "warning:[BandedGlobalAligner] attempted to perform topological sort on empty graph" << endl;
        return;
    }
    
    // initialize return value
    out_topological_order = vector<Node*>(g.node_size());
    size_t order_index = g.node_size() - 1;
    
    // initialize iteration structures
    vector<bool> enqueued = vector<bool>(g.node_size());
    vector<int> edge_index = vector<int>(g.node_size());
    vector<int64_t> stack;
    
    // iterate through starting nodes
    for (int64_t init_node_id = 0; init_node_id < g.node_size(); init_node_id++) {
        if (enqueued[init_node_id]) {
            continue;
        }
        // navigate through graph with DFS
        stack.push_back(init_node_id);
        enqueued[init_node_id] = true;
        while (!stack.empty()) {
            int64_t node_id = stack[stack.size() - 1];
            if (edge_index[node_id] < node_edges_out[node_id].size()) {
                int64_t target_id = node_edges_out[node_id][edge_index[node_id]];
                if (enqueued[target_id]) {
                    edge_index[node_id]++;
                }
                else {
                    stack.push_back(target_id);
                    enqueued[target_id] = true;
                }
            }
            else {
                // add to topological order in reverse finishing order
                stack.pop_back();
                out_topological_order[order_index] = g.mutable_node(node_id);
                order_index--;
            }
        }
    }
}

void BandedGlobalAligner::path_lengths_to_sinks(const string& read, vector<vector<int64_t>>& node_edges_in,
                                                vector<int64_t>& shortest_path_to_sink,
                                                vector<int64_t>& longest_path_to_sink) {
#ifdef debug_banded_aligner_graph_processing
    cerr << "[BandedGlobalAligner::path_lengths_to_sinks]: finding longest and shortest paths to sink node" << endl;
#endif
    
    // find the longest path from the right side of each matrix to the end of the graph
    longest_path_to_sink = vector<int64_t>(topological_order.size());
    shortest_path_to_sink = vector<int64_t>(topological_order.size());
    
    // set initial values -- 0 default value is sufficient for longest path
    for (int64_t& initial_path_length :  shortest_path_to_sink) {
        initial_path_length = numeric_limits<int64_t>::max();
    }
    // set base case (longest path already set to 0)
    for (Node* node : sink_nodes) {
        shortest_path_to_sink[node_id_to_idx.at(node->id())] = 0;
    }
    
    // iterate in reverse order
    for (auto iter = topological_order.rbegin(); iter != topological_order.rend(); iter++) {
        Node* node = *iter;
        int64_t node_seq_len = node->sequence().length();
        int64_t node_idx = node_id_to_idx.at(node->id());
        // compute longest path through this node to right side of incoming matrices
        int64_t longest_path_length = longest_path_to_sink[node_idx] + node_seq_len;
        int64_t shortest_path_length = shortest_path_to_sink[node_idx] + node_seq_len;
        
#ifdef debug_banded_aligner_graph_processing
        cerr << "[BandedGlobalAligner::path_lengths_to_sinks]: processing node " << node->id() << " at index " << node_idx << " with longest/shortest distance to sink " << longest_path_to_sink[node_idx] << "/" << shortest_path_to_sink[node_idx] << " and sequence " << node->sequence() << " for total path length of " << longest_path_length << "/" << shortest_path_length << endl;
#endif
        
        for (int64_t node_in_idx : node_edges_in[node_idx] ) {
            if (longest_path_to_sink[node_in_idx] < longest_path_length) {
                
#ifdef debug_banded_aligner_graph_processing
                cerr << "[BandedGlobalAligner::path_lengths_to_sinks]: path through " << node->id() << " of length " << longest_path_length << " to node at index " << node_in_idx << " is longer than current longest path " << longest_path_to_sink[node_in_idx] << ", updating it now" << endl;
#endif
                
                longest_path_to_sink[node_in_idx] = longest_path_length;
            }
            
            if (shortest_path_to_sink[node_in_idx] > shortest_path_length) {
                
#ifdef debug_banded_aligner_graph_processing
                cerr << "[BandedGlobalAligner::path_lengths_to_sinks]: path through " << node->id() << " of length " << shortest_path_length << " to node at index " << node_in_idx << " is shorter than current shortest path " << shortest_path_to_sink[node_in_idx] << ", updating it now" << endl;
#endif
                
                shortest_path_to_sink[node_in_idx] = shortest_path_length;
            }
        }
    }
}


// fills vectors with whether nodes are masked by the band width, and the band ends of each node
void BandedGlobalAligner::find_banded_paths(const string& read, bool permissive_banding,
                                            vector<vector<int64_t>>& node_edges_in,
                                            vector<vector<int64_t>>& node_edges_out,
                                            int64_t band_padding, vector<bool>& node_masked,
                                            vector<pair<int64_t, int64_t>>& band_ends) {
    
    // find the longest and shortest path from each node to any sink
    vector<int64_t> shortest_path_to_sink;
    vector<int64_t> longest_path_to_sink;
    path_lengths_to_sinks(read, node_edges_in, shortest_path_to_sink, longest_path_to_sink);
    
    // keeps track of which nodes cannot reach the bottom corner within the band
    node_masked = vector<bool>(topological_order.size());
    
    // the bottom and top indices of the band in the rightmost column of each node's matrix
    band_ends = vector<pair<int64_t, int64_t>>(topological_order.size());
    
    // set band ends to identities of max / min functions
    for (int64_t i = 0; i < topological_order.size(); i++) {
        band_ends[i].first = numeric_limits<int64_t>::max();
        band_ends[i].second = numeric_limits<int64_t>::min();
    }
    
    
    if (permissive_banding) {
        // initialize with wide enough bands that every source can hit every connected sink
        for (Node* init_node : source_nodes) {
            int64_t init_node_idx = node_id_to_idx.at(init_node->id());
            int64_t init_node_seq_len = init_node->sequence().length();
            band_ends[init_node_idx].first = min(-band_padding,
                                                 (int64_t) read.length() - (init_node_seq_len + longest_path_to_sink[init_node_idx]) - band_padding);
            band_ends[init_node_idx].second = max(band_padding,
                                                  (int64_t) read.length() - (init_node_seq_len + shortest_path_to_sink[init_node_idx]) + band_padding);
#ifdef debug_banded_aligner_graph_processing
            cerr << "[BandedGlobalAligner::find_banded_paths]: initializing band path end at node " << init_node->id() << " at index " << init_node_idx << " to top " << band_ends[init_node_idx].first << ", and bottom " << band_ends[init_node_idx].second << " from shortest and longest paths of length " << shortest_path_to_sink[init_node_idx] << " and " << longest_path_to_sink[init_node_idx] << " compared to read length " << read.length() << " with padding " << band_padding << endl;
#endif
        }
        
    }
    else {
        // initialize with band ends beginning with source nodes
        for (Node* init_node : source_nodes) {
            int64_t init_node_idx = node_id_to_idx.at(init_node->id());
            int64_t init_node_seq_len = init_node->sequence().length();
            band_ends[init_node_idx].first = -band_padding;
            band_ends[init_node_idx].second = band_padding;
#ifdef debug_banded_aligner_graph_processing
            cerr << "[BandedGlobalAligner::find_banded_paths]: initializing band path end at node " << init_node->id() << " at index " << init_node_idx << " to top " << band_ends[init_node_idx].first << ", and bottom " << band_ends[init_node_idx].second << endl;
#endif
        }
    }
    
    
    // iterate through the rest of the nodes in topological order
    for (int64_t i = 0; i < topological_order.size(); i++) {
        Node* node = topological_order[i];
        int64_t node_idx = node_id_to_idx.at(node->id());
        int64_t node_seq_len = node->sequence().length();
        vector<int64_t>& edges_out = node_edges_out[node_idx];
        
        int64_t extended_band_top = band_ends[node_idx].first + node_seq_len;
        int64_t extended_band_bottom = band_ends[node_idx].second + node_seq_len;
        
#ifdef debug_banded_aligner_graph_processing
        cerr << "[BandedGlobalAligner::find_banded_paths]: following edges out of node " << node->id() << " at index " << node_idx << " with sequence " << node->sequence() << ", band of " << band_ends[node_idx].first << ", " << band_ends[node_idx].second << " extending to " << extended_band_top << ", " << extended_band_bottom << endl;
#endif
        // can alignments from this node reach the bottom right corner within the band?
        if (extended_band_top + shortest_path_to_sink[node_idx] > (int64_t) read.length()
            || extended_band_bottom + longest_path_to_sink[node_idx] < (int64_t) read.length()) {
            
            node_masked[node_idx] = true;
            
#ifdef debug_banded_aligner_graph_processing
            cerr << "[BandedGlobalAligner::find_banded_paths]: cannot complete alignment to read of length " << read.length() << " along shortest path " << shortest_path_to_sink[node_idx] << " or longest path " << longest_path_to_sink[node_idx] << ", which reach range " << extended_band_top + shortest_path_to_sink[node_idx] << ", " << extended_band_bottom + longest_path_to_sink[node_idx] << endl;
#endif
            continue;
        }
        
        // check if each edge out requires expanding the bands
        for (int64_t j = 0; j < edges_out.size(); j++) {
            int64_t node_out_idx = edges_out[j];
            
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
        }
    }
}


// returns the shortest sequence from any source node to each node
void BandedGlobalAligner::shortest_seq_paths(vector<vector<int64_t>>& node_edges_out,
                                             vector<int64_t>& seq_lens_out) {
    
    // initialize vector to store sequence lengths
    seq_lens_out = vector<int64_t>(topological_order.size());
    
    for (auto iter = topological_order.begin(); iter != topological_order.end(); iter++) {
        Node* node = *iter;
        int64_t node_idx = node_id_to_idx.at(node->id());
        int64_t seq_len = node->sequence().length() + seq_lens_out[node_idx];
        
        for (int64_t target_idx : node_edges_out[node_idx]) {
            if (seq_len < seq_lens_out[target_idx]) {
                seq_lens_out[target_idx] = seq_len;
            }
        }
    }
}

void BandedGlobalAligner::align(int8_t* score_mat, int8_t* nt_table, int8_t gap_open, int8_t gap_extend) {
    
    // small enough number to never be accepted in alignment but also not trigger underflow
    int8_t max_mismatch = numeric_limits<int8_t>::max();
    for (int i = 0; i < 25; i++) {
        max_mismatch = min<int8_t>(max_mismatch, score_mat[i]);
    }
    int8_t min_inf = numeric_limits<int8_t>::min() + max<int8_t>((int8_t) -max_mismatch, max<int8_t>(gap_open, gap_extend));
    
    
    // fill each nodes matrix in topological order
    //for (Node* node : topological_order) {
    for (int64_t i = 0; i < topological_order.size(); i++) {
        Node* node = topological_order[i];
        int64_t node_idx = node_id_to_idx.at(node->id());
        BAMatrix* band_matrix = banded_matrices[node_idx];
#ifdef debug_banded_aligner_fill_matrix
        cerr << "[BandedGlobalAligner::align] checking node " << node->id() << " at index " << node_idx << " with sequence " << node->sequence() << " and topological position " << i << endl;
#endif
        
        // skip masked nodes
        if (band_matrix == nullptr) {
#ifdef debug_banded_aligner_fill_matrix
            cerr << "[BandedGlobalAligner::align] node is masked, skipping" << endl;
#endif
            continue;
        }
#ifdef debug_banded_aligner_fill_matrix
        cerr << "[BandedGlobalAligner::align] node is not masked, filling matrix" << endl;
#endif
        band_matrix->fill_matrix(score_mat, nt_table, gap_open, gap_extend, adjust_for_base_quality, min_inf);
    }
    
    traceback(score_mat, nt_table, gap_open, gap_extend, min_inf);
}

void BandedGlobalAligner::traceback(int8_t* score_mat, int8_t* nt_table, int8_t gap_open, int8_t gap_extend, int8_t min_inf) {
    
    // get the sink node matrices for alignment stack
    vector<BAMatrix*> sink_node_matrices;
    sink_node_matrices.reserve(sink_nodes.size());
    for (Node* node : sink_nodes) {
        sink_node_matrices.push_back(banded_matrices[node_id_to_idx[node->id()]]);
    }
    
    // find the optimal alignment(s) and initialize stack
    AltTracebackStack traceback_stack(max_multi_alns, sink_node_matrices);
    
    for (; traceback_stack.has_next(); traceback_stack.next()) {
        int64_t end_node_id;
        matrix_t end_matrix;
        traceback_stack.get_alignment_start(end_node_id, end_matrix);
        int64_t end_node_idx = node_id_to_idx[end_node_id];
        
#ifdef debug_banded_aligner_traceback
        cerr << "[BandedGlobalAligner::traceback] beginning traceback ending at node " << end_node_id << " in matrix " << (end_matrix == Match ? "match" : (end_matrix == InsertCol ? "insert column" : "insert row")) << endl;
#endif
        
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
        
        // add score to alignment
        next_alignment->set_score(traceback_stack.current_traceback_score());
        
        // do traceback
        BABuilder builder(*next_alignment);
        banded_matrices[end_node_idx]->traceback(builder, traceback_stack, end_matrix, score_mat, nt_table,
                                                 gap_open, gap_extend, adjust_for_base_quality, min_inf);
        
        // construct the alignment path
        builder.finalize_alignment();
        
        if (alt_alignments) {
            if (alt_alignments->empty()) {
                // copy the primary into the alternates
                alt_alignments->emplace_back(*next_alignment);
            }
        }
    }
}

// TODO: should find a way to initialize it with the optimal alignments from the sink nodes in the constructor
BandedGlobalAligner::AltTracebackStack::AltTracebackStack(int64_t max_multi_alns,
                                                          vector<BAMatrix*> sink_node_matrices) :
                                                          max_multi_alns(max_multi_alns)
{
    
    // an empty trace back prefix for the initial tracebacks
    vector<Deflection> null_prefix;
    
    // check tracebacks for alignments ending in all sink nodes
    for (BAMatrix* band_matrix : sink_node_matrices) {
        if (band_matrix->match == nullptr) {
            cerr << "error:[BandedGlobalAligner] must fill dynamic programming matrices before finding optimal score" << endl;
            assert(0);
        }
        

        // get the coordinates of the bottom right corner
        Node* node = band_matrix->node;
        int64_t node_id = node->id();
        
        const string& read = band_matrix->alignment.sequence();
        int64_t ncols = node->sequence().length();
        
        int64_t final_col = ncols - 1;
        int64_t final_row = band_matrix->bottom_diag + ncols > (int64_t) read.length() ? (int64_t) read.length() - band_matrix->top_diag - ncols : band_matrix->bottom_diag - band_matrix->top_diag;
        
        int64_t final_idx = final_row * ncols + final_col;
        
        // let the insert routine figure out which one is the best and which ones to keep in the stack
        insert_traceback(null_prefix, band_matrix->match[final_idx],
                         node_id, final_row, final_col, node_id, Match);
        insert_traceback(null_prefix, band_matrix->insert_row[final_idx],
                         node_id, final_row, final_col, node_id, InsertRow);
        insert_traceback(null_prefix, band_matrix->insert_col[final_idx],
                         node_id, final_row, final_col, node_id, InsertCol);
    }
    
    // initialize the traceback trackers for the optimal traceback
    curr_traceback = alt_tracebacks.begin();
    curr_deflxn = (*curr_traceback).first.begin();
}

BandedGlobalAligner::AltTracebackStack::~AltTracebackStack() {
    // nothing to do
}

inline void BandedGlobalAligner::AltTracebackStack::get_alignment_start(int64_t& node_id, matrix_t& matrix) {
    
    // move to next traceback
    if (curr_traceback != alt_tracebacks.end()) {
        // get the start node the first deflection
        curr_deflxn = (*curr_traceback).first.begin();
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

inline void BandedGlobalAligner::AltTracebackStack::next() {
    if (curr_deflxn != (*curr_traceback).first.end()) {
        cerr << "warning:[BandedGlobalAligner] moving on to next alternate alignment without taking all deflections" << endl;
    }

    curr_traceback++;
}

inline bool BandedGlobalAligner::AltTracebackStack::has_next() {
    return (curr_traceback != alt_tracebacks.end());
}

inline void BandedGlobalAligner::AltTracebackStack::propose_deflection(const int8_t score, const int64_t from_node_id,
                                                                       const int64_t row_idx, const int64_t col_idx,
                                                                       const int64_t to_node_id, const matrix_t to_matrix) {
    // only propose deflections if we're going through a new untraversed section of the traceback
    if (curr_deflxn != (*curr_traceback).first.end()) {
        return;
    }
    
    // is the score good enough to be put on the stack?
    if (score <= alt_tracebacks.back().second && alt_tracebacks.size() >= max_multi_alns) {
        return;
    }
    
    insert_traceback((*curr_traceback).first, score, from_node_id, row_idx, col_idx, to_node_id, to_matrix);
}

inline void BandedGlobalAligner::AltTracebackStack::insert_traceback(const vector<Deflection>& traceback_prefix,
                                                                     const int8_t score, const int64_t from_node_id,
                                                                     const int64_t row_idx, const int64_t col_idx,
                                                                     const int64_t to_node_id, const matrix_t to_matrix) {
#ifdef debug_banded_aligner_traceback
    cerr << "[AltTracebackStack::insert_traceback] adding traceback with score " << (int) score << ", new deflection at (" << row_idx << ", " << col_idx << ") on node " << from_node_id << " to " << (to_matrix == Match ? "match" : (to_matrix == InsertRow ? "insert row" : "insert column")) << " matrix on node " << to_node_id << endl;
#endif
    
    // find position in stack where this should go
    auto insert_after = alt_tracebacks.rbegin();
    while (score > (*insert_after).second) {
        insert_after++;
        if (insert_after == alt_tracebacks.rend()) {
            break;
        }
    }
    
    // insert if score is high enough or stack is not full yet
    if (insert_after != alt_tracebacks.rbegin() || alt_tracebacks.size() < max_multi_alns) {
        
        // create a new traceback here
        auto new_traceback = alt_tracebacks.emplace(insert_after.base(), vector<Deflection>(), score);
        
        vector<Deflection>& deflections = (*new_traceback).first;
        
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
    for (pair<vector<Deflection>, int8_t> trace : alt_tracebacks) {
        cerr << (int) trace.second << " -> ";
    }
    cerr << endl;
#endif
}

inline int8_t BandedGlobalAligner::AltTracebackStack::current_traceback_score() {
    return (*curr_traceback).second;
}

inline bool BandedGlobalAligner::AltTracebackStack::at_next_deflection(int64_t node_id, int64_t row_idx,
                                                                       int64_t col_idx) {
    // taken all deflections already?
    if (curr_deflxn == (*curr_traceback).first.end()) {
        return false;
    }
    // correct coordinates?
    if (node_id != (*curr_deflxn).from_node_id) {
        return false;
    }
    if (row_idx != (*curr_deflxn).row_idx) {
        return false;
    }
    if (col_idx != (*curr_deflxn).col_idx) {
        return false;
    }
    return true;
}

inline BandedGlobalAligner::matrix_t BandedGlobalAligner::AltTracebackStack::deflect_to_matrix() {
    matrix_t mat = (*curr_deflxn).to_matrix;
    curr_deflxn++;
    return mat;
}

inline BandedGlobalAligner::matrix_t BandedGlobalAligner::AltTracebackStack::deflect_to_matrix(int64_t& to_node_id) {
    matrix_t mat = (*curr_deflxn).to_matrix;
    to_node_id = (*curr_deflxn).to_node_id;
    curr_deflxn++;
    return mat;
}

BandedGlobalAligner::AltTracebackStack::Deflection::Deflection(const int64_t from_node_id,
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

BandedGlobalAligner::AltTracebackStack::Deflection::~Deflection() {
    // nothing to do
}



