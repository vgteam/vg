#include "filter.hpp"

using namespace std;

namespace vg{
  /**
  * CLI: vg filter -d 10 -q 40 -r -R
  * -r: track depth of both novel variants and those in the graph.
  * -R: remove edits that fail the filter (otherwise toss the whole alignment)
  */
    Alignment Filter::depth_filter(Alignment& aln){
      Path path = aln.path();
      //TODO handle reversing mappings

      for (int i = 0; i < path.mapping_size(); i++){
        Mapping mapping = path.mapping(i);
        Position start_pos = mapping.position();
        int64_t start_node = start_pos.node_id();
        int64_t start_offset = start_pos.offset();
        int64_t curr_offset_in_graph = 0;
        int64_t curr_offset_in_alignment = 0;

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
      double pctid = 0.0;
      graph_len = my_vg->total_length_of_nodes();

    }
}
