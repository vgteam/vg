#include "rare_variant_simplifier.hpp"

namespace vg {

using namespace std;

RareVariantSimplifier::RareVariantSimplifier(MutablePathDeletableHandleGraph& graph, VcfBuffer& variant_source) : Progressive(), graph(graph), variant_source(variant_source) {
    // Nothing to do!
}

void RareVariantSimplifier::simplify() {
    // This holds the IDs of all the nodes we want to keep around
    unordered_set<id_t> to_keep;

    graph.for_each_path_handle([&](const path_handle_t& path) {
        // For each path
        
        if (!Paths::is_alt(graph.get_path_name(path))) {
            // If it isn't an alt path, we want to trace it

            // For each occurrence from start to end
            // Put the ID of the node we are visiting in the to-keep set
            for (handle_t handle : graph.scan_path(path)) {
                to_keep.insert(graph.get_id(handle));
            }
        }
    });

    variant_source.fill_buffer();
    while(variant_source.get() != nullptr) {
        // For each variant
        auto* variant = variant_source.get();

        // TODO: We will use getInfoValueFloat from vcflib, but that API
        // doesn't have support for sniffing the existence of fields, so we
        // have to grab stuff manually too.

        const map<string, vector<string> >& info = variant->info;

        // Count the AF, AC, and AN values
        size_t af_count = info.count("AF") ? info.at("AF").size() : 0;
        size_t ac_count = info.count("AC") ? info.at("AC").size() : 0;
        size_t an_count = info.count("AN") ? info.at("AN").size() : 0;

        // Default to keeping this variant
        bool keep = true;

        if (keep && min_frequency_to_keep != 0) {
            // A frequency condition is imposed, so we have to compute frequency

            // Determine the total frequency of alt alleles of this variant.
            // Sum that AF values if they exist, and sum the AC and AN values and divide otherwise.
            double frequency = 0;

            if (af_count == 0 && (ac_count == 0 || an_count == 0)) {
                // Bail out if any variant doesn't have a well-defined frequency
                cerr << "error[vg::RareVariantSimplifier]: variant at " << variant->sequenceName << ":"
                    << variant->position << " is missing sufficient AF, AC, and/or AN INFO tags to compute frequency" << endl;
                exit(1);
            }

            if (af_count >= ac_count && af_count >= an_count) {
                // AF is the best tag to use
                for (size_t i = 0; i < af_count; i++) {
                    frequency += variant->getInfoValueFloat("AF", i);
                }
            } else {
                // We have to use AC and AN

                if (ac_count != an_count) {
                    cerr << "error[vg::RareVariantSimplifier]: variant at " << variant->sequenceName << ":"
                        << variant->position << " has " << ac_count << " AC values but " << an_count << " AN values. "
                        << "Can't compute frequency!" << endl;
                    exit(1);
                }

                size_t ac_total = 0;
                size_t an_total = 0;

                for (size_t i = 0; i < ac_count; i++) {
                    // Sum up the AC and AN values.
                    // TODO: vcflib has no way to get an int except as a float.
                    ac_total += (size_t) variant->getInfoValueFloat("AC", i);
                    an_total += (size_t) variant->getInfoValueFloat("AN", i);
                }

                if (an_total == 0) {
                    // There are no calls so we can't compute a frequency
                    cerr << "error[vg::RareVariantSimplifier]: variant at " << variant->sequenceName << ":"
                        << variant->position << " has total AN of 0."
                        << "Can't compute frequency!" << endl;
                    exit(1);
                }

                // Compute the frequency
                frequency = (double) ac_total / (double) an_total;
            }

            // Keep only if frequency is sufficiently high
            keep = keep && (frequency >= min_frequency_to_keep);
        }

        if (keep && min_count_to_keep != 0) {
            // A count condition is imposed, so we need to check the raw count.

            if (ac_count == 0) {
                // If we have no AC info we can't do this threshold
                cerr << "error[vg::RareVariantSimplifier]: variant at " << variant->sequenceName << ":"
                    << variant->position << " has no AC INFO tag; can't apply count threshold." << endl;
                exit(1);
            }
            
            // Sum up the AC tags.
            // TODO: if using both conditions we do this loop twice!
            size_t ac_total = 0;
            for (size_t i = 0; i < ac_count; i++) {
                ac_total += (size_t) variant->getInfoValueFloat("AC", i);
            }

            // Keep the variant if it has a sufficient count
            keep = keep && (ac_total >= min_count_to_keep);
        }

        // Now we know if we are keeping the variant or not

        // Work out the variant's alt path names, to either trace or destroy
        vector<string> variant_alt_paths;

        // TODO: this should be factored out of here, construction, and GBWT generation into some central place.
        string var_name = make_variant_id(*variant);
        variant_alt_paths.push_back("_alt_" + var_name + "_0");
        for (size_t alt_index = 1; alt_index < variant->alleles.size(); alt_index++) {
            variant_alt_paths.push_back("_alt_" + var_name + "_" + to_string(alt_index));
        }

        if (keep) {
            // If it is sufficiently common, mark all its alt path nodes as to-keep
            for (auto& path_name : variant_alt_paths) {
                // For each alt path
                if (!graph.has_path(path_name)) {
                    // Skip those that do not exist
                    continue;
                }

                // For each occurrence from start to end
                // Put the ID of the node we are visiting in the to-keep set
                for (handle_t handle : graph.scan_path(graph.get_path_handle(path_name))) {
                    to_keep.insert(graph.get_id(handle));
                }
            }
        } else {
            // Otherwise delete all its alt paths and also its ref path
            for (auto& path_name : variant_alt_paths) {
                if (!graph.has_path(path_name)) {
                    // This path doesn't exist at all, so skip it
                    continue;
                }
                
                // For each alt path
                path_handle_t path = graph.get_path_handle(path_name);

                // Destroy it
                graph.destroy_path(path);

            }

            // The nodes will get destroyed if nothing else sufficiently frequent visited them.
        }

        // Load the next variant
        variant_source.handle_buffer();
        variant_source.fill_buffer();
    }
    
    graph.for_each_handle([&](const handle_t& handle) {
        // After going through all the variants, delete all nodes that aren't to-keep
        if (!to_keep.count(graph.get_id(handle))) {
            graph.destroy_handle(handle);
        }
    });
}

}
