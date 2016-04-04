#include "filter.hpp"

using namespace std;

namespace vg{
    Alignment Filter::depth_filter(Alignment& aln){
      Path path = aln.path();

      for (int i = 0; i < path.mapping_size(); i++){
        Mapping mapping = path.mapping(i);
        Position start_pos = mapping.position();
        id_t start_node = start_pos.node_id();
        if (remove_failing_edits){
          return Alignment();
        }
        else {
          Alignment edited_aln = Alignment();
          return edited_aln;
        }
      }

      return aln;

    }
    Alignment& Filter::qual_filter(Alignment& aln){

    }
    Alignment& Filter::coverage_filter(Alignment& aln){

    }
    Alignment& Filter::avg_qual_filter(Alignment& aln){

    }
    Alignment& Filter::soft_clip_filter(Alignment& aln){

    }
    Alignment& Filter::percent_identity_filter(Alignment& aln){

    }
}
