#ifndef DECON_HPP
#define DECON_HPP
#include <vector>
#include <array>
#include <list>
#include <string>
#include <iostream>
#include <unordered_map>
#include <climits>
#include <queue>
#include <fstream>
#include <cstdlib>
#include <sstream>
#include <stack>
#include "Variant.h"
#include "index.hpp"
#include "path.hpp"
#include "vg.hpp"
#include "vg.pb.h"
#include "Fasta.h"
#include "xg.hpp"
#include "position.hpp"
#include "vcfheader.hpp"
/**
* Deconstruct is getting rewritten.
* New functionality:
* -Detect superbubbles and bubbles
* -Fix command line interface.
* -harmonize on XG / raw graph (i.e. deprecate index)
* -Use unroll/DAGify if needed to avoid cycles

** Much of this is taken from Brankovic's
** "Linear-Time Superbubble Identification Algorithm for Genome Assembly"
*/
namespace vg{
    using namespace std;
    using namespace vcfh;
    struct SuperBubble{
      //A vector of topologically-sorted nodes in the superbubble.
      map<int, vector<id_t> > level_to_nodes;
      id_t start_node;
      id_t end_node;
      bool isNested;
    };

    class Deconstructor{
        public:

            Deconstructor();
            Deconstructor(VG* graph);
            ~Deconstructor();
            void unroll_my_graph(int steps);
            void dagify_my_graph(int steps);
            vg::VG* compact(vg::VG* graph);
            bool is_nested(SuperBubble sb);
            SuperBubble report_superbubble(int64_t start, int64_t end);
            vector<SuperBubble> get_all_superbubbles();
            void sb2vcf(vector<SuperBubble> sbs);
            

        private:

		  VG* my_vg;
          map<id_t, pair<id_t, bool> > my_translation;
          map<id_t, pair<id_t, bool> > my_unroll_trans;
          map<id_t, pair<id_t, bool> > my_dagify_trans;

          vector<SuperBubble> my_super_bubbles;
          size_t my_max_length;
          size_t my_max_component_length;

		  vector<int64_t> nt_to_ids(deque<NodeTraversal>& nt);
          id_t translate_id(id_t id);
		  void init();

    };
}
#endif
