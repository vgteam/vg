#include "build_index.hpp"
#include "source_sink_overlay.hpp"
#include "utility.hpp"

namespace vg {

void build_gcsa_lcp(const HandleGraph& graph,
                    gcsa::GCSA*& gcsa,
                    gcsa::LCPArray*& lcp,
                    int kmer_size,
                    size_t doubling_steps,
                    size_t size_limit,
                    const string& base_file_name) {
    
    // Add an overlay with the source and sink nodes for GCSA
    SourceSinkOverlay overlay(&graph, kmer_size);
    gcsa::ConstructionParameters params;
    params.setSteps(doubling_steps);
    params.setLimit(size_limit);

    // Generate the kmers and reduce the size limit by their size.
    size_t kmer_bytes = params.getLimitBytes();
    string tmpfile = write_gcsa_kmers_to_tmpfile(overlay, kmer_size,
                                                 kmer_bytes,
                                                 overlay.get_id(overlay.get_source_handle()),
                                                 overlay.get_id(overlay.get_sink_handle()),
                                                 base_file_name);
    params.reduceLimit(kmer_bytes);

    // set up the input graph using the kmers
    gcsa::InputGraph input_graph({ tmpfile }, true, params);
    // run the GCSA construction
    gcsa = new gcsa::GCSA(input_graph, params);
    // and the LCP array construction
    lcp = new gcsa::LCPArray(input_graph, params);
    // delete the temporary debruijn graph file
    temp_file::remove(tmpfile);
    // results returned by reference
}

}
