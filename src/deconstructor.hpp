#ifndef DECON_HPP
#define DECON_HPP
#include <vector>
#include <set>
#include <array>
#include <list>
#include <string>
#include <iostream>
#include <unordered_map>
#include <map>
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
/** \file
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
          class Deconstructor{
        public:

            Deconstructor();
            Deconstructor(VG* graph);
            ~Deconstructor();

            void deconstruct(string refpath);
            void deconstruct(vector<string> refpaths); 

        private:
		  VG* my_vg;




    };
}
#endif
