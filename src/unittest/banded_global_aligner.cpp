//
//  banded_global_aligner.cpp
//  
//
//  Created by Jordan Eizenga on 7/26/16.
//
//

#include <stdio.h>
#include "catch.hpp"
#include "gssw_aligner.hpp"
#include "vg.hpp"
#include "banded_global_aligner.hpp"

using namespace google::protobuf;

namespace vg {
    namespace unittest {
        
        
        TEST_CASE( "Banded global alignment works on a DAG with source and sink",
                   "[alignment][banded][DAG] ") {
            
            Aligner aligner = Aligner();
            
            VG graph;
            
            Node* n0 = graph.create_node("GGTG");
            Node* n1 = graph.create_node("C");
            Node* n2 = graph.create_node("A");
            Node* n3 = graph.create_node("TGAAGT");
            
            // bubble structure
            graph.create_edge(n0, n1);
            graph.create_edge(n0, n2);
            graph.create_edge(n1, n3);
            graph.create_edge(n2, n3);
            
            string read = string("GGTGCTGAAGT");
            
            int band_width = 1;
            BandedGlobalAlignmentGraph banded_aligner = BandedGlobalAlignmentGraph(read,
                                                                                   graph.graph,
                                                                                   band_width);
            
            banded_aligner.align(aligner.score_matrix, aligner.nt_table, aligner.gap_open,
                                 aligner.gap_extension);
            
            cerr << banded_aligner.traceback(aligner.score_matrix, aligner.nt_table, aligner.gap_open,
                                             aligner.gap_extension);
            
            REQUIRE(true);
        }
    }
}
