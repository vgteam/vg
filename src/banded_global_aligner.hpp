//
//  banded_global_aligner.hpp
//  
//
//  Created by Jordan Eizenga on 7/25/16.
//
//

#ifndef banded_global_aligner_hpp
#define banded_global_aligner_hpp

#include <stdio.h>
#include <ctype.h>
#include <iostream>
#include <vector>
#include <unordered_set>
#include <unordered_map>
#include <sstream>
#include <algorithm>
#include "vg.pb.h"

#endif /* banded_global_aligner_hpp */

using namespace std;

namespace vg {
    class BandedAlignmentMatrix {
    private:
        
        enum matrix_t {Match, InsertCol, InsertRow};
        
        // these indicate the diagonals in this matrix that the band passes through
        // the bottom index is inclusive
        int64_t top_diag;
        int64_t bottom_diag;
        
        Node* node;
        
        const char* read;
        int64_t read_seq_len;
        
        int64_t cumulative_seq_len;
        
        BandedAlignmentMatrix** seeds;
        int64_t num_seeds;
        
        int8_t* match;
        int8_t* insert_col;
        int8_t* insert_row;
        
        void traceback_internal(stringstream& strm, int64_t start_row, int64_t start_col, matrix_t start_mat,
                                bool in_lead_gap, int8_t* score_mat, int8_t* nt_table, int8_t gap_open,
                                int8_t gap_extend);
        
        void print_matrix(matrix_t which_mat);
        void print_band(matrix_t which_mat);
        
    public:
        BandedAlignmentMatrix(string& read, Node* node, int64_t top_diag, int64_t bottom_diag,
                              BandedAlignmentMatrix** seeds, int64_t num_seeds, int64_t cumulative_seq_len);
        BandedAlignmentMatrix();
        ~BandedAlignmentMatrix();
        
        bool is_masked();
        
        void fill_matrix(int8_t* score_mat, int8_t* nt_table, int8_t gap_open, int8_t gap_extend);
        
        int8_t final_score();
        
        // TODO: coordinate with Erik about traceback semantics
        void traceback(stringstream& strm, int8_t* score_mat, int8_t* nt_table, int8_t gap_open, int8_t gap_extend);
        
        // debugging functions
        void print_full_matrices();
        void print_rectangularized_bands();
        
        friend class BandedAlignmentBuilder;
    };
    
    class BandedAlignmentBuilder {
    private:
        Alignment& alignment;
        const char* read;
        
        matrix_t matrix_state;
        Node* node;
        int64_t edit_length;
        int64_t edit_end;
        
    public:
        void update_state()
    };
    
    class BandedGlobalAlignmentGraph {
    public:
        
        BandedGlobalAlignmentGraph(string& read, Graph& g, int64_t band_padding, bool permissive_banding = false);
        ~BandedGlobalAlignmentGraph();
        
        // perform dynamic programming through the band
        void align(int8_t* score_mat, int8_t* nt_table, int8_t gap_open, int8_t gap_extend);
        // get traceback string after aligning
        string traceback(int8_t* score_mat, int8_t* nt_table, int8_t gap_open, int8_t gap_extend);
        
    private:
        vector<BandedAlignmentMatrix*> banded_matrices;
        
        unordered_map<int64_t, int64_t> node_id_to_idx;
        vector<Node*> topological_order;
        unordered_set<Node*> source_nodes;
        unordered_set<Node*> sink_nodes;
        
        // construction functions
        void graph_edge_lists(Graph& g, bool outgoing_edges, vector<vector<int64_t>>& out_edge_list);
        void topological_sort(Graph& g, vector<vector<int64_t>>& node_edges_out, vector<Node*>& out_topological_order);
        void path_lengths_to_sinks(string& read, vector<vector<int64_t>>& node_edges_in,
                                   vector<int64_t>& shortest_path_to_sink, vector<int64_t>& longest_path_to_sink);
        void find_banded_paths(string& read, bool permissive_banding, vector<vector<int64_t>>& node_edges_in,
                               vector<vector<int64_t>>& node_edges_out, int64_t band_padding,
                               vector<bool>& node_masked, vector<pair<int64_t, int64_t>>& band_ends);
//        void find_banded_paths_permissive(string& read, vector<vector<int64_t>>& node_edges_in,
//                                          int64_t band_padding, vector<pair<int64_t, int64_t>>& band_ends);
        void shortest_seq_paths(vector<vector<int64_t>>& node_edges_out, vector<int64_t>& seq_lens_out);
    };
}


