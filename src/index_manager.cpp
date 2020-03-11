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

IndexManager::IndexManager(const string& fasta_filename, const string& vcf_filename) {
    set_fasta_filename(fasta_filename);
    set_vcf_filename(vcf_filename);
}

void IndexManager::set_fasta_filename(const string& filename) {
    fasta_filename = filename;
    if (!fasta_filename.empty()) {
        // Work out the basename from the FASTA name, which may be .fa or .fa.gz
        pair<string, string> parts = split_ext(fasta_filename);
        if (parts.second == "gz") {
            // Split off whatever it was before the .gz
            parts = split_ext(parts.first);
        }
        basename = parts.first;
    }
}

void IndexManager::set_vcf_filename(const string& filename) {
    vcf_filename = filename;
}

void IndexManager::set_minimizer_override(const string& filename) {
    minimizer_override = filename;
}

void IndexManager::set_gbwtgraph_override(const string& filename) {
    gbwtgraph_override = filename;
}

void IndexManager::set_gbwt_override(const string& filename) {
    gbwt_override = filename;
}

void IndexManager::set_distance_override(const string& filename) {
    distance_override = filename;
}

void IndexManager::set_snarls_override(const string& filename) {
    snarls_override = filename;
}

void IndexManager::set_graph_override(const string& filename) {
    graph_override = filename;
    
    if (!graph_override.empty()) {
        // Graph basename should be used instead
        // TODO: keep it if someone comes along and sets the FASTA afterward.
        pair<string, string> parts = split_ext(fasta_filename);
        if (!parts.first.empty()) {
            // We aren't working with ".xg" or something.
            basename = parts.first;
        }
    }
}

shared_ptr<gbwtgraph::DefaultMinimizerIndex> IndexManager::get_minimizer() {
    ensure_minimizer();
    return minimizer;
}

shared_ptr<gbwtgraph::GBWTGraph> IndexManager::get_gbwtgraph() {
    ensure_gbwtgraph();
    return gbwtgraph.first;
}

shared_ptr<gbwt::GBWT> IndexManager::get_gbwt() {
    ensure_gbwt();
    return gbwt;
}

shared_ptr<vg::MinimumDistanceIndex> IndexManager::get_distance() {
    ensure_distance();
    return distance;
}

shared_ptr<vg::SnarlManager> IndexManager::get_snarls() {
    ensure_snarls();
    return snarls;
}

shared_ptr<PathHandleGraph> IndexManager::get_graph() {
    ensure_graph();
    return graph;
}

string IndexManager::get_filename(const string& extension) const {
    // Assume the indexes are all arranged as basename.ext
    if (basename.empty()) {
        return ""; 
    }
    return basename + "." + extension;
}

template<typename IndexHolderType>
void IndexManager::ensure(IndexHolderType& member, const string& filename_override, const string& extension,
    const function<void(istream&)>& load, const function<void(ostream&)>& make_and_save) {
    if (member) {
        // Already made
        return;
    }


    // Work out where to try to load from
    string input_filename;

    if (!filename_override.empty()) {
        // Just use the override
        input_filename = filename_override;
    } else {
        // Try to get it based on a basename.
        input_filename = get_filename(extension);
    }

    ifstream in(input_filename);
    if (in) {
        // Load the item
        if (show_progress) {
            cerr << "Loading " << extension << " from " << input_filename << endl;
        }
        load(in);
    } else {
        // Make the item and save it

        ofstream out;

        string output_filename = get_filename(extension);
        if (!output_filename.empty()) {
            // User expects us to write
            out.open(output_filename);
            if (!out) {
                throw runtime_error("Cound not write to " + output_filename);
            }
        }
        
        if (show_progress) {
            cerr << "Building " << extension;
            if (out) {
                cerr << " to " << output_filename;
            }
            cerr << endl;
        }
        make_and_save(out);
    }
}

void IndexManager::ensure_graph() {
    ensure(graph, graph_override, "vg", [&](istream& in) {
        // Load the graph
        auto loaded = vg::io::VPKG::load_one<handlegraph::PathHandleGraph>(in);
        // Make it owned by the shared_ptr
        graph.reset(loaded.release());
    }, [&](ostream& out) {
        // Make the graph from the FASTA and VCF

        assert(!fasta_filename.empty());
        
        // Make a graph and give ownership of it to the shared_ptr
        bdsg::HashGraph* mutable_graph = new bdsg::HashGraph();
        graph.reset(mutable_graph);
        
        Constructor constructor;
        constructor.alt_paths = true;
        constructor.max_node_size = 32;

        // Construct the graph.
        // TODO: We can't send a temporary vector to a constt reference for some reason.
        vector<string> fasta_filenames{fasta_filename};
        vector<string> vcf_filenames{};
        if (!vcf_filename.empty()) {
            // We actually have variants
            vcf_filenames.push_back(vcf_filename);
        }
        vector<string> insertion_filenames{};
        constructor.construct_graph(fasta_filenames, vcf_filenames, insertion_filenames, mutable_graph);
        
        if (out) {
            // Save the graph
            vg::io::save_handle_graph(graph.get(), out);
        }
    });
}

void IndexManager::ensure_snarls() {
    ensure(snarls, snarls_override,  "snarls", [&](istream& in) {
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
        
        if (out) {
            // Save the snarls
            snarls->serialize(out);
        }
    });
}

void IndexManager::ensure_distance() {
    ensure(distance, distance_override, "dist", [&](istream& in) {
        // Load distance index from the file
        auto loaded = vg::io::VPKG::load_one<MinimumDistanceIndex>(in);
        distance.reset(loaded.release());
    }, [&](ostream& out) {
        // Make and save

        ensure_graph();
        ensure_snarls();

        // Make it
        distance = make_shared<MinimumDistanceIndex>(graph.get(), snarls.get());
        
        if (out) {
            // Save it
            vg::io::VPKG::save(*distance, out);
        }
    });
}

void IndexManager::ensure_gbwt() {
    ensure(gbwt, gbwt_override, "gbwt", [&](istream& in) {
        // Load GBWT from the file
        auto loaded = vg::io::VPKG::load_one<gbwt::GBWT>(in);
        gbwt.reset(loaded.release());
    }, [&](ostream& out) {
        // Make and save

        ensure_graph();
        assert(!vcf_filename.empty());

        HaplotypeIndexer indexer;
        indexer.show_progress = show_progress;

        // Make it, making sure to convert to non-dynamic GBWT
        map<string, Path> alt_paths;
        vector<string> insertion_filenames;
        auto built = indexer.build_gbwt(graph.get(), alt_paths, false, vcf_filename, insertion_filenames);
        gbwt = make_shared<gbwt::GBWT>(*built);
        built.reset();

        if (out) {
            // Save it
            vg::io::VPKG::save(*gbwt, out);
        }
    });
}

void IndexManager::ensure_gbwtgraph() {
    ensure(gbwtgraph.first, gbwtgraph_override, "gg", [&](istream& in) {
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

        if (out) {
            // Save it
            vg::io::VPKG::save(*gbwtgraph.first, out);
        }
    });
}

void IndexManager::ensure_minimizer() {
    ensure(minimizer, minimizer_override, "min", [&](istream& in) {
        // Load minimizer index from the file
        auto loaded = vg::io::VPKG::load_one<gbwtgraph::DefaultMinimizerIndex>(in);
        minimizer.reset(loaded.release());
    }, [&](ostream& out) {
        // Make and save

        ensure_gbwtgraph();

        // Make it
        minimizer.reset(new gbwtgraph::DefaultMinimizerIndex(minimizer_k, minimizer_w));
        gbwtgraph::index_haplotypes(*gbwtgraph.first, *minimizer);

        if (out) {
            // Save it
            vg::io::VPKG::save(*minimizer, out);
        }
    });
}

}
