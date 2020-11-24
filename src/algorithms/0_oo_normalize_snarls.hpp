#include "../gbwt_helper.hpp"
#include "../handle.hpp"
#include "../subgraph.hpp"
#include "../vg.hpp"
#include "count_walks.hpp"
#include <string>
#include <gbwtgraph/gbwtgraph.h>


namespace vg {

class SnarlNormalizer {
  public:
    virtual ~SnarlNormalizer() = default;

    SnarlNormalizer(MutablePathDeletableHandleGraph &graph, const gbwtgraph::GBWTGraph &haploGraph,
                    const int &max_alignment_size = 200,
                    const string &path_finder = "GBWT" /*alternative is "exhaustive"*/);

    virtual void normalize_top_level_snarls(ifstream &snarl_stream);

    virtual vector<int> normalize_snarl(const id_t &source_id, const id_t &sink_id);

  protected:
    // member variables:
    // the handle graph with snarls to normalize
    MutablePathDeletableHandleGraph &_graph;
    // GBWT graph with snarls to normalize, includes the embedded threads needed for the
    // GBWTPathFinder approach.
    const gbwtgraph::GBWTGraph &_haploGraph;
    // the maximum number of threads allowed to align in a given snarl. If the number of
    // threads exceeds this threshold, the snarl is skipped.
    int _max_alignment_size;
    id_t _cur_source_id;
    id_t _cur_sink_id;
    const string &_path_finder;

    tuple<vector<vector<handle_t>>, vector<vector<handle_t>>, unordered_set<handle_t>>
    extract_gbwt_haplotypes(const SubHandleGraph &snarl,
                            const id_t &source_id, const id_t &sink_id);

    pair<vector<string>, unordered_set<handle_t>> source_to_sink_exhaustive_path_finder();

    vector<vector<handle_t>>
    find_haplotypes_not_at_source(unordered_set<handle_t> &touched_handles,
                                  const id_t &sink_id);

    vector<string> format_handle_haplotypes_to_strings(
        const vector<vector<handle_t>> &haplotype_handle_vectors);

    VG align_source_to_sink_haplotypes(vector<string> source_to_sink_haplotypes);

    void force_maximum_handle_size(MutableHandleGraph &graph, const size_t &max_size);

    vector<pair<step_handle_t, step_handle_t>>
    extract_embedded_paths_in_snarl(const PathHandleGraph &graph, const id_t &source_id,
                                    const id_t &sink_id);

    SubHandleGraph extract_subgraph(const HandleGraph &graph, const id_t &start_id,
                                    const id_t &end_id);

    void integrate_snarl(const HandleGraph &new_snarl,
                         const vector<pair<step_handle_t, step_handle_t>> embedded_paths);

    void move_path_to_snarl(const pair<step_handle_t, step_handle_t> &old_embedded_path,
                            vector<handle_t> &new_snarl_handles, id_t &new_source_id,
                            id_t &new_sink_id, const id_t &old_source_id,
                            const id_t &old_sink_id);

    bool source_and_sink_handles_map_properly(
        const HandleGraph &graph, const id_t &new_source_id, const id_t &new_sink_id,
        const bool &touching_source, const bool &touching_sink,
        const handle_t &potential_source, const handle_t &potential_sink);

    vector<int> check_handle_as_start_of_path_seq(const string &handle_seq,
                                                  const string &path_seq);

    // -------------------------------- DEBUG CODE BELOW:
    // ------------------------------------

    pair<vector<handle_t>, vector<handle_t>>
    debug_get_sources_and_sinks(const HandleGraph &graph);
};

} // namespace vg
