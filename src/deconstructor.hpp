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
#include <fstream>
#include <cstdlib>

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
      vector<id_t> nodes;
      id_t start_node;
      id_t end_node;
      bool isNested;
    };

class CandList{


};
    class Deconstructor{
        public:

            Deconstructor();
            Deconstructor(xg::XG x);
			Deconstructor(VG v);
            ~Deconstructor();
            SuperBubble report_superbubble(int64_t start, int64_t end);
            vector<SuperBubble> get_all_superbubbles();
            

        private:
          xg::XG my_xg;
		  VG my_vg;
          vector<SuperBubble> my_super_bubbles;
		  vector<int64_t> nt_to_ids(deque<NodeTraversal>& nt);
		  void init();

    };
}
#endif
