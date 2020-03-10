/**
 * \file index_manager.cpp: implementations of common indexing functionality
 */

#include <iostream>
#include <vector>
#include <string>

#include <vg/io/vpkg.hpp>
#include <vg/io/stream.hpp>
#include <bdsg/hash_graph.hpp>

#include "index_manager.hpp"
#include "utility.hpp"
#include "constructor.hpp"
#include "io/save_handle_graph.hpp"


using namespace std;

namespace vg {

IndexManager::IndexManager(const string& fasta_filename, const string& vcf_filename) : fasta_filename(fasta_filename), vcf_filename(vcf_filename) {
    // Work out the basename from the FASTA name, which may be .fa or .fa.gz
    pair<string, string> parts = split_ext(fasta_filename);
    if (parts.second == "gz") {
        // Split off whatever it was before the .gz
        parts = split_ext(parts.first);
    }
    basename = parts.first;
}

string IndexManager::get_filename(const string& extension) const {
    // Assume the indexes are all next to the FASTA, with the FASTA extension borken off and this one added on.
    return basename + "." + extension;
}

void IndexManager::ensure_graph() {
    
    if (graph) {
        // Already made
        return;
    }

    // Work out where to load from/save to
    string graph_filename = get_filename("vg");
    
    ifstream in(graph_filename);
    if (in) {
        // Load the graph
        auto loaded = vg::io::VPKG::load_one<handlegraph::PathHandleGraph>(in);
        // Make it owned by the shared_ptr
        graph.reset(loaded.release());
    } else {
        // Make the graph from the FASTA and VCF
        
        // Make sure we will be able to save the graph
        ofstream out(graph_filename);
        if (!out) {
            throw runtime_error("Cound not save graph to " + graph_filename);
        }
        
        // Make a graph and give ownership of it to the shared_ptr
        bdsg::HashGraph* mutable_graph = new bdsg::HashGraph();
        graph.reset(mutable_graph);
        
        Constructor constructor;
        constructor.alt_paths = true;
        constructor.max_node_size = 32;

        // Construct the graph.
        // TODO: We can't send a temporary vector to a constt reference for some reason.
        vector<string> fasta_filenames{fasta_filename};
        vector<string> vcf_filenames{vcf_filename};
        vector<string> insertion_filenames{};
        constructor.construct_graph(fasta_filenames, vcf_filenames, insertion_filenames, mutable_graph);
        
        // Save the graph
        vg::io::save_handle_graph(graph.get(), out);
    }
}

}
