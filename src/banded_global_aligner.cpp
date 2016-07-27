//
//  banded_global_aligner.cpp
//  
//
//  Created by Jordan Eizenga on 7/25/16.
//
//

#include <ctype.h>
#include "banded_global_aligner.hpp"

using namespace vg;

BandedAlignmentMatrix::BandedAlignmentMatrix(string& read, Node* node,
                                             int64_t top_diag, int64_t bottom_diag,
                                             BandedAlignmentMatrix** seeds, int64_t num_seeds,
                                             int64_t cumulative_seq_len) :
                                             node(node),
                                             read_seq_len(read.length()),
                                             top_diag(top_diag),
                                             bottom_diag(bottom_diag),
                                             seeds(seeds),
                                             read(read.c_str()),
                                             num_seeds(num_seeds),
                                             cumulative_seq_len(cumulative_seq_len),
                                             match(nullptr),
                                             insert_col(nullptr),
                                             insert_row(nullptr)
{
    // nothing to do
}

BandedAlignmentMatrix::BandedAlignmentMatrix() : top_diag(0), bottom_diag(0), node(nullptr),
                                                 read(nullptr), read_seq_len(0), cumulative_seq_len(0),
                                                 seeds(nullptr), match(nullptr), insert_col(nullptr),
                                                 insert_row(nullptr)
{
    // default constructor
}

BandedAlignmentMatrix::~BandedAlignmentMatrix() {
    free(match);
    free(insert_row);
    free(insert_col);
    free(seeds);
    // note: BandedAlignmentMatrix doesn't own read sequence, don't free it
}

bool BandedAlignmentMatrix::is_masked() {
    // node will only be nullptr if masked
    return node == nullptr;
}

void BandedAlignmentMatrix::fill_matrix(int8_t* score_mat, int8_t* nt_table,
                                        int8_t gap_open, int8_t gap_extend) {
    // note: bottom has the higher index
    int64_t band_height = bottom_diag - top_diag + 1;
    int64_t ncols = node->sequence().length();
    int64_t band_size = band_height * ncols;
    
    const char* node_seq = node->sequence().c_str();
    
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
     *  this way the dimensions are band_height x seq_len instead of read_len x seq_len
     *  also note that the internal structure of each column is preserved and each row is 
     *  a diagonal in the original matrix
     */
    
    // TODO: calculating min_inf could be costly since we do it at every node and it looks through
    // the whole score matrix. need to move this out somewhere
    
    // find maximum mismatch penalty
    int8_t max_mismatch = numeric_limits<int8_t>::max();
    for (int i = 0; i < 25; i++) {
        max_mismatch = min<int8_t>(max_mismatch, score_mat[i]);
    }
    
    // small enough number to never be accepted in alignment but also not trigger underflow
    int8_t min_inf = numeric_limits<int8_t>::min() + max<int8_t>((int8_t) -max_mismatch, max<int8_t>(gap_open, gap_extend));
    
    int64_t idx, up_idx, diag_idx, left_idx;
    
    // POA into left hand column from the seeds
    if (seeds == nullptr) {
        if (cumulative_seq_len != 0) {
            cerr << "error:[BandedGlobalAligner] banded alignment has no node predecessor for node in middle of path" << endl;
            exit(EXIT_FAILURE);
        }
        
        // first node in alignment, fill out initial column based on implied lead gaps
        
        // find position of the first cell in the rectangularized band
        int64_t iter_start = -top_diag;
        idx = iter_start * ncols;
        
        // cap stop index if last diagonal is below bottom of matrix
        int64_t iter_stop = bottom_diag > read_seq_len ? band_height + read_seq_len - bottom_diag - 1 : band_height;
        
        // match of first nucleotides
        match[idx] = score_mat[5 * nt_table[node_seq[0]] + nt_table[read[0]]];
        // only way to end an alignment in a gap here is to row and column gap
        insert_col[idx] = -2 * gap_open;
        insert_row[idx] = -2 * gap_open;
        
        for (int64_t i = iter_start + 1; i < iter_stop; i++) {
            idx = i * ncols;
            up_idx = (i - 1) * ncols;
            // score of a match in this cell
            int8_t match_score = score_mat[5 * nt_table[node_seq[0]] + nt_table[read[top_diag + 1]]];
            // must take one lead gap to get into first column
            match[idx] = match_score - gap_open - (top_diag + i - 1) * gap_extend;
            // normal iteration along column
            insert_row[idx] = max(max(match[up_idx] - gap_open, insert_row[up_idx] - gap_extend),
                                  insert_col[up_idx] - gap_open);
            // must take two gaps to get into first column
            insert_col[idx] = -2 * gap_open - (top_diag + i) * gap_extend;
        }
    }
    else {
        
        // are we clipping any diagonals because they are outside the range of the matrix in this column?
        bool bottom_diag_outside = bottom_diag >= read_seq_len;
        bool top_diag_outside_or_abutting = top_diag <= 0;
        
        // rows in the rectangularized band corresponding to diagonals
        int64_t iter_start = top_diag_outside_or_abutting ? -top_diag : 0;
        int64_t iter_stop = bottom_diag_outside ? band_height + read_seq_len - bottom_diag - 1 : band_height;
        
        // handle special logic for lead column insertions
        idx = iter_start * ncols;
        
        int8_t match_score = score_mat[5 * nt_table[node_seq[0]] + nt_table[read[iter_start + top_diag]]];
        if (top_diag_outside_or_abutting) {
            // match after implied gap along top edge
            match[idx] = match_score - gap_open - (cumulative_seq_len - 1) * gap_extend;
            // gap open after implied gap along top edge
            insert_row[idx] = -2 * gap_open - cumulative_seq_len * gap_extend;
        }
        else {
            // no lead gaps possible so seed with identity of max function to prepare for POA
            insert_row[idx] = min_inf;
            match[idx] = min_inf;
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
            BandedAlignmentMatrix* seed = seeds[seed_num];
            if (seed == nullptr) {
                continue;
            }
            
            int64_t seed_seq_len = seed->node->sequence().length();
            
            int64_t seed_next_top_diag = seed->top_diag + seed_seq_len;
            int64_t seed_next_bottom_diag = seed->bottom_diag + seed_seq_len;
            
            int64_t seed_iter_start = seed_next_top_diag < 0 ? -seed_next_top_diag : 0;
            int64_t seed_iter_stop = seed_next_bottom_diag >= read_seq_len ? band_height + seed_next_bottom_diag - bottom_diag - 1 : band_height;
            
            for (int64_t i = seed_iter_start; i < seed_iter_stop - 1; i++) {
                idx = i * ncols;
                // are we at the first row in the band?
                if (i != iter_start) {
                    // can extend a match
                    diag_idx = i * ncols + (seed_seq_len - 1);
                    match_score = score_mat[5 * nt_table[node_seq[0]] + nt_table[read[i + top_diag]]];
                    
                    match[idx] = max<int8_t>(match_score + max<int8_t>(max<int8_t>(seed->match[diag_idx],
                                                                                   seed->insert_row[diag_idx]),
                                                                       seed->insert_col[diag_idx]), match[idx]);
                    
                }
                // are we at the last row in this seed?
                if (i != seed_iter_stop - 1) {
                    // can extend a column insertion
                    left_idx = (i + 1) * ncols + (seed_seq_len - 1);
                    
                    insert_col[idx] = max<int8_t>(max<int8_t>(max<int8_t>(match[left_idx] - gap_open,
                                                                          insert_row[left_idx] - gap_open),
                                                              insert_col[left_idx] - gap_extend), insert_col[idx]);
                }
            }
            
        }
        
        // compute the insert row scores (they can be computed after the POA iterations since they do not
        // cross node boundaries)
        for (int64_t i = iter_start + 1; i < iter_stop; i++) {
            idx = i * ncols;
            up_idx = (i - 1) * ncols;
            
            insert_row[idx] = max(max(match[up_idx] - gap_open, insert_row[up_idx] - gap_extend),
                                  insert_col[up_idx] - gap_open);
        }
    }
    
    // iterate through the rest of the columns
    for (int64_t j = 1; j < ncols; j++) {
        
        // are we clipping any diagonals because they are outside the range of the matrix in this column?
        bool bottom_diag_outside = bottom_diag + j >= read_seq_len;
        bool top_diag_outside_or_abutting = top_diag + j <= 0;
        
        int64_t iter_start = top_diag_outside_or_abutting ? -(top_diag + j) : 0;
        int64_t iter_stop = bottom_diag_outside ? band_height + read_seq_len - bottom_diag - j - 1 : band_height;
        
        idx = iter_start * ncols + j;
        
        int8_t match_score = score_mat[5 * nt_table[node_seq[j]] + nt_table[read[iter_start + top_diag]]];
        if (top_diag_outside_or_abutting) {
            // gap open after implied gap along top edge
            insert_row[idx] = -2 * gap_open - (cumulative_seq_len + j) * gap_extend;
            // match after implied gap along top edge
            match[idx] = match_score - gap_open - (cumulative_seq_len + j - 1) * gap_extend;
        }
        else {
            diag_idx = iter_start * ncols + (j - 1);
            // cannot reach this node with row insert (outside the diagonal)
            insert_row[idx] = min_inf;
            // cells should be present to do normal diagonal iteration
            match[idx] = match_score + max(max(match[diag_idx], insert_row[diag_idx]), insert_col[diag_idx]);
        }
        
        // normal iteration along row
        int64_t left_idx = (iter_start + 1) * ncols + (j - 1);
        insert_col[idx] = max(max(match[left_idx] - gap_open, insert_row[left_idx] - gap_open),
                              insert_col[left_idx] - gap_extend);
        
        for (int64_t i = iter_start + 1; i < iter_stop - 1; i++) {
            // indices of the current and previous cells in the rectangularized band
            idx = i * ncols + j;
            up_idx = (i - 1) * ncols + j;
            diag_idx = i * ncols + (j - 1);
            left_idx = (i + 1) * ncols + (j - 1);
            
            match_score = score_mat[5 * nt_table[node_seq[j]] + nt_table[read[i + top_diag]]];
            
            match[idx] = match_score + max(max(match[diag_idx], insert_row[diag_idx]), insert_col[diag_idx]);
            
            insert_row[idx] = max(max(match[up_idx] - gap_open, insert_row[up_idx] - gap_extend),
                                  insert_col[up_idx] - gap_open);
            
            insert_col[idx] = max(max(match[left_idx] - gap_open, insert_row[left_idx] - gap_open),
                                  insert_col[left_idx] - gap_extend);
        }
        
        // stop iteration one cell early to handle logic on bottom edge of band
        
        // skip this step in edge case where read length is 1
        if (iter_stop - 1 > iter_start) {
            idx = (iter_stop - 1) * ncols + j;
            up_idx = (iter_stop - 2) * ncols + j;
            diag_idx = (iter_stop - 1) * ncols + (j - 1);
            
            match_score = score_mat[5 * nt_table[node_seq[j]] + nt_table[read[iter_stop + top_diag - 1]]];
            
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
}

int8_t BandedAlignmentMatrix::final_score() {
    if (match == nullptr) {
        cerr << "error:[BandedGlobalAligner] must fill dynamic programming matrices before finding optimal score" << endl;
        exit(EXIT_FAILURE);
    }
    int64_t ncols = node->sequence().length();
    int64_t final_row = bottom_diag + ncols > read_seq_len ? 2 * read_seq_len - bottom_diag - ncols - 1 : read_seq_len - 1;
    int64_t final_idx = final_row * ncols + (ncols - 1);
    return max(max(match[final_idx], insert_row[final_idx]), insert_col[final_idx]);
}

void BandedAlignmentMatrix::traceback(stringstream& strm, int8_t* score_mat, int8_t* nt_table, int8_t gap_open,
                                      int8_t gap_extend) {
    
    // get coordinates of bottom right corner
    int64_t ncols = node->sequence().length();
    int64_t row = bottom_diag + ncols > read_seq_len ? 2 * read_seq_len - bottom_diag - ncols - 1 : read_seq_len - 1;
    int64_t col = ncols - 1;
    int64_t idx = row * ncols + col;
    
    // which matrix achieves the highest score?
    matrix start_mat = Match;
    int8_t highest_score = match[idx];
    if (insert_row[idx] > highest_score) {
        highest_score = insert_row[idx];
        start_mat = InsertRow;
    }
    if (insert_col[idx] > highest_score) {
        highest_score = insert_col[idx];
        start_mat = InsertCol;
    }
    
    traceback_internal(strm, row, col, start_mat, false, score_mat, nt_table, gap_open, gap_extend);
}

void BandedAlignmentMatrix::traceback_internal(stringstream& strm, int64_t start_row, int64_t start_col,
                                               matrix start_mat, bool in_lead_gap, int8_t* score_mat, int8_t* nt_table,
                                               int8_t gap_open, int8_t gap_extend) {
    
    
    int64_t band_height = bottom_diag - top_diag + 1;
    const char* node_seq = node->sequence().c_str();
    int64_t ncols = node->sequence().length();
    
    strm << "]";
    
    int64_t idx, next_idx;
    int64_t i = start_row, j = start_col;
    matrix curr_mat = start_mat;
    int8_t curr_score;
    
    // do node traceback unless we are in the lead gap implied at the edge of the DP matrix
    bool continue_node_traceback = !in_lead_gap;
    
    while (continue_node_traceback) {
        idx = i * ncols + j;
        
        switch (curr_mat) {
            case Match:
            {
                // add uppercase char for match
                strm << read[i + top_diag];
                
                curr_score = match[idx];
                next_idx = i * ncols + j - 1;
                int8_t match_score = score_mat[5 * nt_table[node_seq[j]] + nt_table[read[i + top_diag]]];
                
                if (i + j == -top_diag) {
                    // at top of matrix, move into implied lead gap along top edge
                    in_lead_gap = true;
                    continue_node_traceback = false;
                    
                }
                else if (match[idx] == match[next_idx] + match_score) {
                    curr_mat = Match;
                }
                else if (match[idx] == insert_row[next_idx] + match_score) {
                    curr_mat = InsertRow;
                }
                else if (match[idx] == insert_col[next_idx] + match_score) {
                    curr_mat = InsertCol;
                }
                else {
                    cerr << "error:[BandedGlobalAligner] traceback stuck in match matrix interior" << endl;
                    exit(EXIT_FAILURE);
                }
                
                j--;
                
                break;
            }
                
            case InsertRow:
            {
                // add lowercase char for read insert
                strm << tolower(read[i + top_diag]);
                
                curr_score = insert_row[idx];
                next_idx = (i - 1) * ncols + j;
                
                if (i + j == -top_diag) {
                    // at top of matrix, move into implied lead gap along top edge
                    in_lead_gap = true;
                    continue_node_traceback = false;
                    
                }
                else if (i == 0) {
                    // along top of band
                    cerr << "error:[BandedGlobalAligner] traceback attempted to leave band from top" << endl;
                    exit(EXIT_FAILURE);
                    
                }
                else {
                    // inside matrix interior
                    if (curr_score == match[next_idx] - gap_open) {
                        curr_mat = Match;
                    }
                    else if (curr_score == insert_row[next_idx] - gap_extend) {
                        curr_mat = InsertRow;
                    }
                    else if (curr_score == insert_col[next_idx] - gap_open) {
                        curr_mat = InsertCol;
                    }
                    else {
                        cerr << "error:[BandedGlobalAligner] traceback stuck in insert row matrix interior" << endl;
                        exit(EXIT_FAILURE);
                    }
                }
                
                i--;
                
                break;
            }
            
            case InsertCol:
            {
                // add read deletion
                strm << '-';
                
                curr_score = insert_col[idx];
                next_idx = (i + 1) * ncols + j - 1;
                
                if (i == band_height - 1) {
                    // along bottom of band
                    cerr << "error:[BandedGlobalAligner] traceback attempted to leave band from bottom" << endl;
                    exit(EXIT_FAILURE);
                    
                }
                else {
                    // inside matrix interior
                    if (curr_score == match[next_idx] - gap_open) {
                        curr_mat = Match;
                    }
                    else if (curr_score == insert_row[next_idx] - gap_open) {
                        curr_mat = InsertRow;
                    }
                    else if (curr_score == insert_col[next_idx] - gap_extend) {
                        curr_mat = InsertCol;
                    }
                    else {
                        cerr << "error:[BandedGlobalAligner] traceback stuck in insert row matrix interior" << endl;
                        exit(EXIT_FAILURE);
                    }
                }
                
                i++;
                j--;
                
                break;
            }
                
            default:
            {
                cerr << "error:[BandedGlobalAligner] unrecognized matrix type given to traceback" << endl;
                exit(EXIT_FAILURE);
                break;
            }
        }
        
        if (j == 0 && curr_mat != InsertRow) {
            // at a node boundary or the final iteration
            continue_node_traceback = false;
        }
    }
    
    if (in_lead_gap) {
        // add lead column gaps until reaching edge of node
        curr_mat = InsertCol;
        while (j > 0) {
            strm << '-';
            j--;
            i++;
        }
    }
    
    BandedAlignmentMatrix* traceback_seed = nullptr;
    int64_t seed_row;
    int64_t seed_col;
    
    // are there are predecessor nodes?
    if (num_seeds > 0) {
        if (in_lead_gap) {
            // we are in the implied lead gap along the top of the matrix
            // take the shortest path back to the origin of the global alignment
            
            // add final read deletion of node
            strm << '-';
            
            int64_t shortest_seq_len = numeric_limits<int64_t>::max();
            for (int64_t k = 0; k < num_seeds; k++) {
                BandedAlignmentMatrix* seed = seeds[k];
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
            seed_row = top_diag - seed_extended_top_diag + i + 1;
            seed_col = seed_ncols - 1;
        }
        else {
            int8_t match_score;
            switch (curr_mat) {
                case Match:
                {
                    // add final match of matrix
                    strm << read[top_diag + i];
                    curr_score = match[i * ncols];
                    match_score = score_mat[5 * nt_table[node_seq[j]] + nt_table[read[i + top_diag]]];
                    break;
                }
                    
                case InsertCol:
                {
                    // add final read deletion of matrix
                    strm << '-';
                    curr_score = insert_col[i * ncols];
                    break;
                }
                    
                case InsertRow:
                {
                    cerr << "error:[BandedGlobalAligner] traceback attempted to open row gap over node boundary" << endl;
                    exit(EXIT_FAILURE);
                    break;
                }
                    
                default:
                {
                    cerr << "error:[BandedGlobalAligner] unrecognized matrix type given to traceback" << endl;
                    exit(EXIT_FAILURE);
                    break;
                }
            }
            
            // matches stay on same diagonal, column insertions move over one diagonal
            // note that the indexing is relative to THIS matrix, not the seed
            int64_t curr_diag = top_diag + i;
            
            // check traceback goes to each seed matrix
            for (int64_t k = 0; k < num_seeds; k++) {
                BandedAlignmentMatrix* seed = seeds[k];
                // is the matrix masked?
                if (seed == nullptr) {
                    continue;
                }
                
                int64_t seed_ncols = seed->node->sequence().length();
                // the diagonals in the current matrix that this seed extends to
                int64_t seed_extended_top_diag = seed->top_diag + seed_ncols;
                int64_t seed_extended_bottom_diag = seed->bottom_diag + seed_ncols;
                
                // does the traceback diagonal extend backward to this matrix?
                if (curr_diag > seed_extended_bottom_diag - (curr_mat == InsertCol) // col inserts hit 1 less of band
                    || curr_diag < seed_extended_top_diag) {
                    continue;
                }
                
                seed_col = seed_ncols - 1;
                seed_row = curr_diag - seed_extended_top_diag + (curr_mat == InsertCol);
                next_idx = seed_row * seed_ncols + seed_col;
                
                // exit loop if we find a matching cell
                if (curr_mat == Match) {
                    if (curr_score == seed->match[next_idx] + match_score) {
                        curr_mat = Match;
                        traceback_seed = seed;
                        break;
                    }
                    else if (curr_score == seed->insert_col[next_idx] + match_score) {
                        curr_mat = InsertCol;
                        traceback_seed = seed;
                        break;
                    }
                    else if (curr_score == seed->insert_row[next_idx] + match_score) {
                        curr_mat = InsertRow;
                        traceback_seed = seed;
                        break;
                    }
                }
                else {
                    // only other matrix this could be is insert_col
                    if (curr_score == seed->match[next_idx] - gap_open) {
                        curr_mat = Match;
                        traceback_seed = seed;
                        break;
                    }
                    else if (curr_score == seed->insert_col[next_idx] - gap_extend) {
                        curr_mat = InsertCol;
                        traceback_seed = seed;
                        break;
                    }
                    else if (curr_score == seed->insert_row[next_idx] - gap_open) {
                        curr_mat = InsertRow;
                        traceback_seed = seed;
                        break;
                    }
                }
                // if no matches found, continue to next seed
            }
        }
    }
    else {
        // this is the first node in the alignment, finish off the alignment and return
        if (in_lead_gap) {
            strm << "-";
        }
        else {
            switch (curr_mat) {
                case Match:
                {
                    strm << read[i];
                    j--;
                    break;
                }
                    
                case InsertCol:
                {
                    strm << '-';
                    j--;
                    i++;
                    break;
                }
                
                default:
                {
                    cerr << "error:[BandedGlobalAligner] invalid matrix type for final traceback column" << endl;
                    exit(EXIT_FAILURE);
                    break;
                }
            }
            
            // add any lead row gaps necessary
            while (top_diag + i >= 0) {
                strm << tolower(read[top_diag + i]);
                i--;
            }
        }
    }
    
    strm << "[:";
    // reverse id string so that ends up forward when entire string is reversed
    string id_string = to_string(node->id());
    reverse(id_string.begin(), id_string.end());
    strm << id_string;
    
    // begin traceback through next node if necessary
    if (num_seeds > 0) {
        traceback_internal(strm, seed_row, seed_col, curr_mat, in_lead_gap, score_mat, nt_table,
                           gap_open, gap_extend);
    }
}

BandedGlobalAlignmentGraph::BandedGlobalAlignmentGraph(string& read, Graph& g,
                                                       int64_t band_width) : band_width(band_width) {
    // TODO: coordinate with Erik to make sure that the last and first character of the MEMs
    // are included in the cut up sub-graph
    
    // convert the graph into lists of ids indicating incoming or outgoing edges
    vector<vector<int64_t>> node_edges_in;
    vector<vector<int64_t>> node_edges_out;
    graph_edge_lists(g, true, node_edges_out);
    graph_edge_lists(g, false, node_edges_in);
    
    // compute topological ordering
    vector<Node*> topological_order;
    topological_sort(g, node_edges_out, topological_order);
    
    // identify source and sink nodes in the graph
    for (int64_t i = 0; i < g.node_size(); i++) {
        if (node_edges_in.empty()) {
            source_nodes.insert(g.mutable_node(i));
        }
        if (node_edges_out.empty()) {
            sink_nodes.insert(g.mutable_node(i));
        }
    }
    
    // figure out what the bands need to be for alignment and which nodes cannot complete a
    // global alignment within the band
    vector<bool> node_masked;
    vector<pair<int64_t, int64_t>> band_ends;
    find_banded_paths(read, node_edges_in, band_width, topological_order, node_masked, band_ends);
    
    // find the shortest sequence leading to each node so we can infer the length
    // of lead deletions
    vector<int64_t> shortest_seqs;
    shortest_seq_paths(topological_order, node_edges_out, shortest_seqs);
    
    // initialize DP matrices for each node
    banded_matrices.reserve(g.node_size());
    for (int64_t i = 0; i < g.node_size(); i++) {
        if (node_masked[i]) {
            banded_matrices.push_back(BandedAlignmentMatrix());
        }
        else {
            Node* node = g.mutable_node(i);
            int64_t node_seq_len = node->sequence().length();
            
            // POA predecessor matrices
            BandedAlignmentMatrix** seeds;
            vector<int64_t>& edges_in = node_edges_in[i];
            if (edges_in.empty()) {
                seeds = nullptr;
            }
            else {
                seeds = (BandedAlignmentMatrix**) malloc(sizeof(BandedAlignmentMatrix**) * edges_in.size());
                for (int64_t j = 0; j < edges_in.size(); j++) {
                    // add sentinel for masked nodes
                    seeds[j] = node_masked[edges_in[j]] ? nullptr : &banded_matrices[edges_in[j]];
                }
            }
            
            banded_matrices[i] = BandedAlignmentMatrix(read,
                                                       node,
                                                       band_ends[i].first - node_seq_len, // shift from last to first column
                                                       band_ends[i].second - node_seq_len,
                                                       seeds,
                                                       edges_in.size(),
                                                       shortest_seqs[i]);
        }
    }
}

BandedGlobalAlignmentGraph::~BandedGlobalAlignmentGraph() {
    // nothing to do
}

// fills a vector with vectors ids that have edges to/from each node
void BandedGlobalAlignmentGraph::graph_edge_lists(Graph& g, bool outgoing_edges, vector<vector<int64_t>>& edge_list) {
    if (outgoing_edges) {
        edge_list = vector<vector<int64_t>>(g.node_size());
        for (int64_t i = 0; i < g.edge_size(); i++) {
            Edge* edge = g.mutable_edge(i);
            edge_list[(int)edge->from()].push_back(edge->to());
        }
    }
    else {
        edge_list = vector<vector<int64_t>>(g.node_size());
        for (int64_t i = 0; i < g.edge_size(); i++) {
            Edge* edge = g.mutable_edge(i);
            edge_list[edge->to()].push_back(edge->from());
        }
    }
}

// standard DFS-based topological sort algorithm
// NOTE: this is only valid if the Graph g has been dag-ified first and there are from_start
// or to_end edges.
void BandedGlobalAlignmentGraph::topological_sort(Graph& g, vector<vector<int64_t>>& node_edges_out,
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

// fills vectors with whether nodes are masked by the band width, and the band ends of each node
void BandedGlobalAlignmentGraph::find_banded_paths(string& read, vector<vector<int64_t>>& node_edges_in,
                                                   int64_t band_width, vector<Node*>& topological_order,
                                                   vector<bool>& node_masked, vector<pair<int64_t, int64_t>>& band_ends) {
    
    // find the longest path from the right side of each matrix to the end of the graph
    // vector initialized with 0's -- sufficient base case for this algorithm
    vector<int64_t> longest_path_to_end = vector<int64_t>(topological_order.size());
    // iterate in reverse order
    for (auto iter = topological_order.rbegin(); iter != topological_order.rend(); iter++) {
        Node* node = *iter;
        int64_t node_seq_len = node->sequence().length();
        int64_t node_id = node->id();
        // compute longest path through this node to right side of incoming matrices
        int64_t path_length = longest_path_to_end[node_id] + node_seq_len;
        for (int64_t node_in_id : node_edges_in[node_id] ) {
            if (longest_path_to_end[node_in_id] < path_length) {
                longest_path_to_end[node_in_id] = path_length;
            }
        }
    }
    
    // keeps track of which nodes cannot reach the bottom corner within the band
    node_masked = vector<bool>(topological_order.size());
    
    // the bottom and top indices of the band in the rightmost column of each node's matrix
    band_ends = vector<pair<int64_t, int64_t>>(topological_order.size());
    
    // set band ends to identities of max / min functions
    for (int64_t i = 0; i < topological_order.size(); i++) {
        band_ends[i].first = numeric_limits<int64_t>::min();
        band_ends[i].second = numeric_limits<int64_t>::max();
    }
    
    int64_t read_len = read.length();
    
    // initial value: band ends after traversing first node
    Node* init_node = topological_order[0];
    int64_t init_node_id = init_node->id();
    int64_t init_node_seq_len = init_node->sequence().length();
    band_ends[init_node_id].first = init_node_seq_len + band_width;
    band_ends[init_node_id].second = init_node_seq_len - band_width;
    
    // iterate through the rest of the nodes in topological order
    for (int64_t i = 1; i < topological_order.size(); i++) {
        Node* node = topological_order[i];
        int64_t node_id = node->id();
        int64_t node_seq_len = node->sequence().length();
        vector<int64_t>& edges_in = node_edges_in[node_id];
        
        // check if each edge in requires extending the bands
        for (int64_t j = 0; j < edges_in.size(); j++) {
            int64_t node_in_id = edges_in[j];
            if (node_masked[node_in_id]) {
                // node is too distant from band to connect in global alignment, skip it
                continue;
            }
            
            // does extending alignments from these nodes require a wider band here?
            int64_t extended_band_bottom = band_ends[node_in_id].first + read_len;
            if (extended_band_bottom > band_ends[node_id].first) {
                band_ends[node_id].first = extended_band_bottom;
            }
            int64_t extended_band_top = band_ends[node_in_id].second + read_len;
            if (extended_band_top < band_ends[node_id].second) {
                band_ends[node_id].second = extended_band_top;
            }
        }
        
        // can alignments from this node reach the bottom right corner within the band?
        if (band_ends[node_id].first > read_len
            || band_ends[node_id].second + longest_path_to_end[node_id] < read_len) {
            node_masked[node_id] = true;
        }
    }
}

// returns the shortest sequence from any source node to each node
void BandedGlobalAlignmentGraph::shortest_seq_paths(vector<Node*>& topological_order,
                                                    vector<vector<int64_t>>& node_edges_out,
                                                    vector<int64_t>& seq_lens_out) {
    
    // initialize vector to store sequence lengths
    seq_lens_out = vector<int64_t>(topological_order.size());
    
    for (auto iter = topological_order.begin(); iter != topological_order.end(); iter++) {
        Node* node = *iter;
        int64_t node_id = node->id();
        int64_t seq_len = node->sequence().length() + seq_lens_out[node_id];
        
        for (int64_t target_id : node_edges_out[node_id]) {
            if (seq_len < seq_lens_out[target_id]) {
                seq_lens_out[target_id] = seq_len;
            }
        }
    }
}

void BandedGlobalAlignmentGraph::align(int8_t* score_mat, int8_t* nt_table, int8_t gap_open, int8_t gap_extend) {
    
    // fill each nodes matrix in topological order
    for (Node* node : topological_order) {
        BandedAlignmentMatrix& band_matrix = banded_matrices[node->id()];
        
        // skip masked nodes
        if (band_matrix.is_masked()) {
            continue;
        }
        
        band_matrix.fill_matrix(score_mat, nt_table, gap_open, gap_extend);
    }
}

string BandedGlobalAlignmentGraph::traceback(int8_t* score_mat, int8_t* nt_table, int8_t gap_open, int8_t gap_extend) {
    
    // find node at end of traceback with optimal score
    int64_t traceback_end_node_id;
    int8_t optimal_score = numeric_limits<int8_t>::min();
    for (Node* node : sink_nodes) {
        int64_t node_id = node->id();
        int8_t node_score = banded_matrices[node_id].final_score();
        if (node_score > optimal_score) {
            optimal_score = node_score;
            traceback_end_node_id = node_id;
        }
    }
    
    // get traceback string
    stringstream strm;
    banded_matrices[traceback_end_node_id].traceback(strm, score_mat, nt_table, gap_open, gap_extend);
    string traceback_string = strm.str();
    
    // string built in reverse, so flip it
    reverse(traceback_string.begin(), traceback_string.end());
    
    return traceback_string;
    
}







