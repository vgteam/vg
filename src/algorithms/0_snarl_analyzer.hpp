#include "../handle.hpp"
#include "../subgraph.hpp"

/**
 * Goal: use the output from this in python script, which will detect which snarls are 
 * split into multiple, where in graph snarls are shrinking or growing, and compare snarls
 * across multiple graphs.
 */
namespace vg {
namespace algorithms {

class SnarlAnalyzer {
  public:
    /// destructor/constructer
    virtual ~SnarlAnalyzer() = default;
    SnarlAnalyzer(const HandleGraph& graph, ifstream &snarl_stream, bool skip_source_sink);

    /// data on snarls:

    // // output mean snarl size; mode snarl size; number of snarls; etc. NOTE: this should be done in python.
    // void snarl_stats();
    
    // save three column tsv with snarl start-position, stop-position, and size of snarl.
    void output_snarl_sizes(string& file_name);

  protected:
    /// member variables:
    // the handle graph with snarls to normalize
    const HandleGraph &_graph;

    // Info on snarls
    vector<int> _snarl_sources;
    vector<int> _snarl_sinks;
    vector<int> _snarl_sizes;

    /// internal functions:
    // todo: un-steal this from SnarlNormalizer.
    // Extracts a subgraph representing a snarl.
    SubHandleGraph extract_subgraph(const HandleGraph &graph,
                                                 id_t source_id,
                                                 id_t sink_id,
                                                 const bool backwards);

};

}}