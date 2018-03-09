#include "build_index.hpp"

namespace vg {

void build_gcsa_lcp(VG& graph,
                    gcsa::GCSA*& gcsa,
                    gcsa::LCPArray*& lcp,
                    int kmer_size,
                    size_t doubling_steps,
                    size_t size_limit,
                    const string& base_file_name) {
    id_t max_id=0;
    graph.for_each_handle([&max_id,&graph](const handle_t& h) { max_id = max(graph.get_id(h), max_id); });
    id_t head_id = max_id+1;
    id_t tail_id = max_id+2;
    Node* head_node = nullptr; Node* tail_node = nullptr;
    // TODO add this for MutableHandleGraphs
    graph.add_start_end_markers(kmer_size, '#', '$', head_node, tail_node, head_id, tail_id);

    gcsa::ConstructionParameters params;
    params.setSteps(doubling_steps);
    params.setLimit(size_limit);

    // Generate the kmers and reduce the size limit by their size.
    size_t kmer_bytes = params.getLimitBytes();
    string tmpfile = write_gcsa_kmers_to_tmpfile(graph, kmer_size,
                                                 kmer_bytes,
                                                 head_id, tail_id,
                                                 base_file_name);
    params.reduceLimit(kmer_bytes);

    graph.destroy_node(head_node);
    graph.destroy_node(tail_node);
    // set up the input graph using the kmers
    gcsa::InputGraph input_graph({ tmpfile }, true);
    // run the GCSA construction
    gcsa = new gcsa::GCSA(input_graph, params);
    // and the LCP array construction
    lcp = new gcsa::LCPArray(input_graph, params);
    // delete the temporary debruijn graph file
    temp_file::remove(tmpfile);
    // results returned by reference
}

}
