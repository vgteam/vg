#include <gbwtgraph/gbwtgraph.h>
#include "../handle.hpp"
#include "../subgraph.hpp"


namespace vg {
namespace algorithms {

class SnarlSequenceFinder {
  public:
    virtual ~SnarlSequenceFinder() = default;

    SnarlSequenceFinder(const SubHandleGraph &snarl,
                   const gbwtgraph::GBWTGraph &haploGraph, const id_t &source_id, 
                   const id_t &sink_id);

    tuple<vector<vector<handle_t>>, vector<vector<handle_t>>, unordered_set<handle_t>>
    find_gbwt_haps();

    pair<vector<string>, unordered_set<handle_t>> find_exhaustive_paths();

  protected:
    // member variables:
    // the handle graph with snarls to normalize
    const SubHandleGraph &_snarl;
    // GBWT graph with snarls to normalize, includes the embedded threads needed for the
    // GBWTPathFinder approach.
    const gbwtgraph::GBWTGraph &_haploGraph;
    const id_t _source_id; 
    const id_t _sink_id; 

    vector<vector<handle_t>> find_haplotypes_not_at_source(unordered_set<handle_t> &touched_handles);
};
}
}