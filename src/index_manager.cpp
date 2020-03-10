/**
 * \file index_manager.cpp: implementations of common indexing functionality
 */

#include <iostream>
#include <vector>
#include <string>

#include <vg/io/vpkg.hpp>
#include <vg/io/stream.hpp>
#include <bdsg/hash_graph.hpp>
#include <gbwtgraph/index.h> 

#include "index_manager.hpp"
#include "utility.hpp"
#include "constructor.hpp"
#include "haplotype_indexer.hpp"
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

template<typename IndexHolderType>
void IndexManager::ensure(IndexHolderType& member, const string& extension, const function<void(istream&)>& load, const function<void(ostream&)>& make_and_save) {
    if (member) {
        // Already made
        return;
    }

    // Work out where to load from/save to
    string filename = get_filename("extension");
    
    ifstream in(filename);
    if (in) {
        // Load the item
        load(in);
    } else {
        // Make the item and save it
        
        // Make sure we will be able to save
        ofstream out(filename);
        if (!out) {
            throw runtime_error("Cound not write to " + filename);
        }
        
        make_and_save(out);
    }
}

void IndexManager::ensure_graph() {
    ensure(graph, "vg", [&](istream& in) {
        // Load the graph
        auto loaded = vg::io::VPKG::load_one<handlegraph::PathHandleGraph>(in);
        // Make it owned by the shared_ptr
        graph.reset(loaded.release());
    }, [&](ostream& out) {
        // Make the graph from the FASTA and VCF
        
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
    });
}

void IndexManager::ensure_snarls() {
    ensure(snarls, "snarls", [&](istream& in) {
        // Load from the file
        snarls = make_shared<SnarlManager>(in);
    }, [&](ostream& out) {
        // Make the snarls and save them

        ensure_graph();

        // Make a snarl finder
        auto finder = make_unique<CactusSnarlFinder>(*graph);
        // Find the snarls and save them
        snarls = make_shared<SnarlManager>(std::move(finder->find_snarls_parallel()));
        // Delete the the snarl finder
        finder.reset();
        
        // Save the snarls
        snarls->serialize(out);
    });
}

void IndexManager::ensure_distance() {
    ensure(distance, "dist", [&](istream& in) {
        // Load distance index from the file
        auto loaded = vg::io::VPKG::load_one<MinimumDistanceIndex>(in);
        distance.reset(loaded.release());
    }, [&](ostream& out) {
        // Make and save

        ensure_graph();
        ensure_snarls();

        // Make it
        distance = make_shared<MinimumDistanceIndex>(graph.get(), snarls.get());
        
        // Save it
        vg::io::VPKG::save(*distance, out);
    });
}

void IndexManager::ensure_gbwt() {
    ensure(gbwt, "gbwt", [&](istream& in) {
        // Load GBWT from the file
        auto loaded = vg::io::VPKG::load_one<gbwt::DynamicGBWT>(in);
        gbwt.reset(loaded.release());
    }, [&](ostream& out) {
        // Make and save

        ensure_graph();
        assert(!vcf_filename.empty());

        HaplotypeIndexer indexer;
        indexer.show_progress = show_progress;

        // Make it
        map<string, Path> alt_paths;
        vector<string> insertion_filenames;
        auto built = indexer.build_gbwt(graph.get(), alt_paths, false, vcf_filename, insertion_filenames);
        gbwt.reset(built.release());

        // Save it
        vg::io::VPKG::save(*gbwt, out);
    });
}

void IndexManager::ensure_gbwtgraph() {
    ensure(gbwtgraph.first, "gg", [&](istream& in) {
        // Make sure GBWT is ready
        ensure_gbwt();

        // Load GBWTGraph from the file
        auto loaded = vg::io::VPKG::load_one<gbwtgraph::GBWTGraph>(in);
        loaded->set_gbwt(*gbwt);
        gbwtgraph.first.reset(loaded.release());
        gbwtgraph.second = gbwt;
    }, [&](ostream& out) {
        // Make and save

        ensure_graph();
        ensure_gbwt();

        // Make it
        gbwtgraph.first.reset(new gbwtgraph::GBWTGraph(*gbwt, *graph));
        gbwtgraph.second = gbwt;

        // Save it
        vg::io::VPKG::save(*gbwtgraph.first, out);
    });
}

void IndexManager::ensure_minimizer() {
    ensure(minimizer, "min", [&](istream& in) {
        // Load minimizer index from the file
        auto loaded = vg::io::VPKG::load_one<gbwtgraph::DefaultMinimizerIndex>(in);
        minimizer.reset(loaded.release());
    }, [&](ostream& out) {
        // Make and save

        ensure_gbwtgraph();

        // Make it
        minimizer.reset(new gbwtgraph::DefaultMinimizerIndex(minimizer_k, minimizer_w));
        gbwtgraph::index_haplotypes(*gbwtgraph.first, *minimizer);

        // Save it
        vg::io::VPKG::save(*minimizer, out);
    });
}

}
