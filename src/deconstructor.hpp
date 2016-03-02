#ifndef DECON_HPP
#define DECON_HPP
#include <vector>
#include <list>
#include <string>
#include <iostream>
#include <unordered_map>
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
    class Deconstructor{
        public:

            Deconstructor();
            Deconstructor(xg::XG x);
            ~Deconstructor();
            SuperBubble report_superbubble(int64_t start, int64_t end);
            vector<SuperBubble> get_all_superbubbles(VG v);
            int64_t validate_super_bubble(int64_t start, int64_t end);
            vector<int64_t> emit_nodes_on_path_through_superbubble(Path p, SuperBubble sb);
            vector<int64_t> emit_nodes_off_path_through_superbubble(SuperBubble sb);
            vector<vcflib::Variant> sb_to_variants(SuperBubble sb);

        private:
          xg::XG my_xg;
            int64_t rangemin(vector<int64_t> a, int i, int j);
            int64_t rangemax(vector<int64_t> a, int i, int j);

    };
}
#endif
