#ifndef VG_HOMOGENIZER
#define VG_HOMOGENIZEER
#include <iostream>
#include <vector>
#include "vg.hpp"
#include "translator.hpp"
#include "filter.hpp"
#include "mapper.hpp"
#include "bubbles.hpp"
#include "vg.pb.h"
#include "types.hpp"


using namespace vg;
using namespace std;

namespace vg{
    class Homogenizer{
        public:
            /* Locates tips in the graph
             * and tries to generate a single
             * edge / node to represent them. 
             * This edge is then added, the offending sequences
             * are remapped, and the process is repeated until the
             * graph becomes stable.
             */
            void homogenize(vg::VG* graph, xg::XG* xindex, gcsa::GCSA* gcsa_index, gcsa::LCPArray* lcp_index, Paths p, int kmer_size);
            void homogenize(vg::VG* graph, xg::XG* xindex, gcsa::GCSA* gcsa_index, gcsa::LCPArray* lcp_index, vg::Index reads_index);
        private:

            Translator translator;
            /** Find tips (nodes with an indegree/outdegree of 0 in the graph */
            vector<vg::id_t> find_tips(vg::VG* graph);

            /** Find non-ref tips */
            vector<vg::id_t> find_non_ref_tips(vg::VG* graph);

            /** remap a set of Alignments to the graph */
            int remap(vector<Alignment> reads, vg::VG graph);
            /** Remove all tips from the graph.
             * WARNING: may cut head/tail nodes.*/
            void cut_tips(vg::VG* graph);
            /** Remove specific nodes and their edges from the graph */
            void cut_tips(vector<id_t> tip_ids, vg::VG* graph);
            /** Remove non-reference tips from the graph. */
            void cut_nonref_tips(vg::VG* graph);
            

    };
}
#endif
