/*
 * construct.cpp: contains implementations for vg construction functions.
 */


#include "vg.hpp"

namespace vg {

using namespace std;


// Implementations of VG functions. TODO: refactor out of VG class

VG::VG(vcflib::VariantCallFile& variantCallFile,
       FastaReference& reference,
       string& target_region,
       bool target_is_chrom,
       int vars_per_region,
       int max_node_size,
       bool flat_input_vcf,
       bool load_phasing_paths,
       bool load_variant_alt_paths,
       bool showprog,
       set<string>* allowed_variants) {

    init();

    omp_set_dynamic(1); // use dynamic scheduling

    show_progress = showprog;

    map<string, VG*> refseq_graph;

    vector<string> targets;
    if (!target_region.empty()) {
        targets.push_back(target_region);
    } else {
        for (vector<string>::iterator r = reference.index->sequenceNames.begin();
             r != reference.index->sequenceNames.end(); ++r) {
            targets.push_back(*r);
        }
    }

    // How many phase paths do we want to load?
    size_t num_phasings = load_phasing_paths ? variantCallFile.sampleNames.size() * 2 : 0;
    // We'll later split these where you would have to take an edge that doesn't exist.

    // to scale up, we have to avoid big string memcpys
    // this could be accomplished by some deep surgery on the construction routines
    // however, that could be a silly thing to do,
    // because why break something that's conceptually clear
    // and anyway, we want to break the works into chunks
    //
    // there is a development that could be important
    // our chunk size isn't going to reach into the range where we'll have issues (>several megs)
    // so we'll run this for regions of moderate size, scaling up in the case that we run into a big deletion
    //

    for (vector<string>::iterator t = targets.begin(); t != targets.end(); ++t) {

        //string& seq_name = *t;
        string seq_name;
        string target = *t;
        int start_pos = 0, stop_pos = 0;
        // nasty hack for handling single regions
        if (!target_is_chrom) {
            parse_region(target,
                         seq_name,
                         start_pos,
                         stop_pos);
            if (stop_pos > 0) {
                if (variantCallFile.is_open()) {
                    variantCallFile.setRegion(seq_name, start_pos, stop_pos);
                }
            } else {
                if (variantCallFile.is_open()) {
                    variantCallFile.setRegion(seq_name);
                }
                stop_pos = reference.sequenceLength(seq_name);
            }
        } else {
            // the user said the target is just a sequence name
            // and is unsafe to parse as it may contain ':' or '-'
            // for example "gi|568815592:29791752-29792749"
            if (variantCallFile.is_open()) {
                variantCallFile.setRegion(target);
            }
            stop_pos = reference.sequenceLength(target);
            seq_name = target;
        }
        vcflib::Variant var(variantCallFile);

        vector<vcflib::Variant>* region = NULL;

        // convert from 1-based input to 0-based internal format
        // and handle the case where we are already doing the whole chromosome
        id_t start = start_pos ? start_pos - 1 : 0;
        id_t end = start;

        create_progress("loading variants for " + target, stop_pos-start_pos);
        // get records
        vector<vcflib::Variant> records;

        // This is going to hold the alleles that occur at certain reference
        // positions, in addition to the reference allele. We keep them ordered
        // so we can refer to them by number.
        map<long,vector<vcflib::VariantAllele> > alleles;

        // This is going to hold, for each position, allele combination, a
        // vector of bools marking which phases of which samples visit that
        // allele. Each sample is stored at (sample number * 2) for phase 0 and
        // (sample number * 2 + 1) for phase 1. The reference may not always get
        // an allele, but if anything is reference it will show up as an
        // overlapping allele elsewhere.
        map<pair<long, int>, vector<bool>> phase_visits;

        // This is going to hold visits to VariantAlleles by the reference and
        // nonreference alts of variants. We map from VariantAllele index and
        // number to a list of the variant ID and alt number pairs that use the
        // VariantAllele.
        map<pair<long, int>, vector<pair<string, int>>> variant_alts;

        // We don't want to load all the vcf records into memory at once, since
        // the vcflib internal data structures are big compared to the info we
        // need.
        int64_t variant_chunk_size = 1000;

        auto parse_loaded_variants = [&]() {
            // Parse the variants we have loaded, and clear them out, so we can
            // go back and load a new batch of variants.

            // decompose records into alleles with offsets against our target
            // sequence Dump the collections of alleles (which are ref, alt
            // pairs) into the alleles map. Populate the phase visit map if
            // we're loading phasing paths, and the variant alt path map if
            // we're loading variant alts.
            vcf_records_to_alleles(records, alleles,
                load_phasing_paths ? &phase_visits : nullptr,
                load_variant_alt_paths ? &variant_alts : nullptr,
                flat_input_vcf);
            records.clear(); // clean up
        };

        int64_t i = 0;
        while (variantCallFile.is_open() && variantCallFile.getNextVariant(var)) {
            // this ... maybe we should remove it as for when we have calls against N
            bool isDNA = allATGC(var.ref);
            for (vector<string>::iterator a = var.alt.begin(); a != var.alt.end(); ++a) {
                if (!allATGC(*a)) isDNA = false;
            }
            // only work with DNA sequences
            if (isDNA) {
                string vrepr = var.vrepr();
                var.position -= 1; // convert to 0-based
                if (allowed_variants == nullptr
                    || allowed_variants->count(vrepr)) {
                    records.push_back(var);
                }
            }
            if (++i % 1000 == 0) update_progress(var.position-start_pos);
            // Periodically parse the records down to what we need and throw away the rest.
            if (i % variant_chunk_size == 0) parse_loaded_variants();
        }
        // Finish up any remaining unparsed variants
        parse_loaded_variants();

        destroy_progress();

        // store our construction plans
        deque<Plan*> construction;
        // so we can check which graphs we can safely append
        set<VG*> graph_completed;
        // we add and remove from graph_completed, so track count for logging
        int graphs_completed = 0;
        int final_completed = -1; // hm
        // the construction queue
        list<VG*> graphq;
        int graphq_size = 0; // for efficiency
        // ^^^^ (we need to insert/remove things in the middle of the list,
        // but we also need to be able to quickly determine its size)
        // for tracking progress through the chromosome
        map<VG*, unsigned long> graph_end;

        create_progress("planning construction", stop_pos-start_pos);
        // break into chunks
        int chunk_start = start;
        bool invariant_graph = alleles.empty();
        while (invariant_graph || !alleles.empty()) {
            invariant_graph = false;
            map<long, vector<vcflib::VariantAllele> > new_alleles;
            map<pair<long, int>, vector<bool>> new_phase_visits;
            map<pair<long, int>, vector<pair<string, int>>> new_variant_alts;
            // our start position is the "offset" we should subtract from the
            // alleles and the phase visits for correct construction
            //chunk_start = (!chunk_start ? 0 : alleles.begin()->first);
            int chunk_end = chunk_start;
            bool clean_end = true;
            for (int i = 0; (i < vars_per_region || !clean_end) && !alleles.empty(); ++i) {
                auto pos = alleles.begin()->first - chunk_start;
                chunk_end = max(chunk_end, (int)alleles.begin()->first);
                auto& pos_alleles = alleles.begin()->second;
                // apply offset when adding to the new alleles
                auto& curr_pos = new_alleles[pos];
                for (int j = 0; j < pos_alleles.size(); j++) {
                    // Go through every allele that occurs at this position, and
                    // update it to the offset position in new_alleles
                    auto& allele = pos_alleles[j];

                    // We'll clone and modify it.
                    auto new_allele = allele;
                    int ref_end = new_allele.ref.size() + new_allele.position;
                    // look through the alleles to see if there is a longer chunk
                    if (ref_end > chunk_end) {
                        chunk_end = ref_end;
                    }
                    new_allele.position = pos;
                    // Copy the modified allele over.
                    // No need to deduplicate.
                    curr_pos.push_back(new_allele);

                    // Also handle any visits to this allele
                    // We need the key, consisting of the old position and the allele number there.
                    auto old_allele_key = make_pair(alleles.begin()->first, j);
                    // Make the new key
                    auto new_allele_key = make_pair(pos, j);
                    if(phase_visits.count(old_allele_key)) {
                        // We have some usages of this allele for phase paths. We need to move them over.

                        // Move over the value and insert into the new map. See <http://stackoverflow.com/a/14816487/402891>
                        // TODO: would it be clearer with the braces instead?
                        new_phase_visits.insert(make_pair(new_allele_key, std::move(phase_visits.at(old_allele_key))));

                        // Now we've emptied out/made-undefined the old vector,
                        // so we probably should drop it from the old map.
                        phase_visits.erase(old_allele_key);
                    }

                    if(variant_alts.count(old_allele_key)) {
                        // We have some usages of this allele by variant alts. We need to move them over.

                        // Do a move operation
                        new_variant_alts.insert(make_pair(new_allele_key, std::move(variant_alts.at(old_allele_key))));
                        // Delete the olkd entry (just so we don't keep it around wasting time/space/being unspecified)
                        variant_alts.erase(old_allele_key);
                    }
                }
                alleles.erase(alleles.begin());
                // TODO here we need to see if we are neighboring another variant
                // and if we are, keep constructing
                if (alleles.begin()->first <= chunk_end) {
                    clean_end = false;
                } else {
                    clean_end = true;
                }
            }
            // record end position, use target end in the case that we are at the end
            if (alleles.empty()) chunk_end = stop_pos;

            // we set the head graph to be this one, so we aren't obligated to copy the result into this object
            // make a construction plan
            Plan* plan = new Plan(graphq.empty() && targets.size() == 1 ? this : new VG,
                                  std::move(new_alleles),
                                  std::move(new_phase_visits),
                                  std::move(new_variant_alts),
                                  reference.getSubSequence(seq_name,
                                                           chunk_start,
                                                           chunk_end - chunk_start),
                                  seq_name);
            chunk_start = chunk_end;
#pragma omp critical (graphq)
            {
                graphq.push_back(plan->graph);
                construction.push_back(plan);
                if (show_progress) graph_end[plan->graph] = chunk_end;
                update_progress(chunk_end);
            }
        }
#ifdef debug
        cerr << omp_get_thread_num() << ": graphq size " << graphq.size() << endl;
#endif
        graphq_size = graphq.size();
        destroy_progress();

        // this system is not entirely general
        // there will be a problem when the regions of overlapping deletions become too large
        // then the inter-dependence of each region will make parallel construction in this way difficult
        // because the chunks will get too large

        // use this function to merge graphs both during and after the construction iteration
        auto merge_first_two_completed_graphs =
            [this, start_pos, &graph_completed, &graphq, &graphq_size, &graph_end, &final_completed](void) {
            // find the first two consecutive graphs which are completed
            VG* first = NULL;
            VG* second = NULL;
//#pragma omp critical (cerr)
//            cerr << omp_get_thread_num() << ": merging" << endl;
#pragma omp critical (graphq)
            {
                auto itp = graphq.begin(); // previous
                auto itn = itp; if (itp != graphq.end()) ++itn; // next
                // scan the graphq to find consecutive entries that are both completed
                while (itp != itn // there is > 1 entry
                       && itn != graphq.end() // we aren't yet at the end
                       && !(graph_completed.count(*itp) // the two we're looking at aren't completed
                            && graph_completed.count(*itn))) {
                    ++itp; ++itn;
                }

                if (itn != graphq.end()) {
                    // we have two consecutive graphs to merge!
                    first = *itp;
                    second = *itn;
                    // unset graph completed for both
                    graph_completed.erase(first);
                    graph_completed.erase(second);
                    graphq.erase(itn);
                    --graphq_size;
                }
            }

            if (first && second) {
                // combine graphs
                first->append(*second);
#pragma omp critical (graphq)
                {
                    if (final_completed != -1) update_progress(final_completed++);
                    graph_completed.insert(first);
                    graph_end.erase(second);
                }
                delete second;
            }
        };

        create_progress("constructing graph", construction.size());

        // (in parallel) construct each component of the graph
#pragma omp parallel for
        for (int i = 0; i < construction.size(); ++i) {

            int tid = omp_get_thread_num();
            Plan* plan = construction.at(i);
#ifdef debug
#pragma omp critical (cerr)
            cerr << tid << ": " << "constructing graph " << plan->graph << " over "
                 << plan->alleles.size() << " variants in " <<plan->seq.size() << "bp "
                 << plan->name << endl;
#endif

            // Make the piece of graph, passing along the number of sample phases if we're making phase paths.
            plan->graph->from_alleles(plan->alleles,
                                      plan->phase_visits,
                                      num_phasings,
                                      plan->variant_alts,
                                      plan->seq,
                                      plan->name);

            // Break up the nodes ourselves
            if(max_node_size > 0) {
                plan->graph->dice_nodes(max_node_size);
            }

#pragma omp critical (graphq)
            {
                update_progress(++graphs_completed);
                graph_completed.insert(plan->graph);
#ifdef debug
#pragma omp critical (cerr)
                cerr << tid << ": " << "constructed graph " << plan->graph << endl;
#endif
            }
            // clean up
            delete plan;

            // concatenate chunks of the result graph together
            merge_first_two_completed_graphs();

        }
        destroy_progress();

        // merge remaining graphs
        final_completed = 0;
        create_progress("merging remaining graphs", graphq.size());
#pragma omp parallel
        {
            bool more_to_merge = true;
            while (more_to_merge) {
                merge_first_two_completed_graphs();
                usleep(10);
#pragma omp critical (graphq)
                more_to_merge = graphq_size > 1;
            }
        }
        destroy_progress();

        // parallel end
        // finalize target

        // our target graph should be the only entry in the graphq
        assert(graphq.size() == 1);
        VG* target_graph = graphq.front();

        // store it in our results
        refseq_graph[target] = target_graph;

        create_progress("joining graphs", target_graph->size());
        // clean up "null" nodes that are used for maintaining structure between temporary subgraphs
        target_graph->remove_null_nodes_forwarding_edges();
        destroy_progress();

        // then use topological sorting and re-compression of the id space to make sure that
        create_progress("topologically sorting", target_graph->size());
        target_graph->sort();
        destroy_progress();

        create_progress("compacting ids", target_graph->size());
        // we get identical graphs no matter what the region size is
        target_graph->compact_ids();
        destroy_progress();

    }

    // hack for efficiency when constructing over a single chromosome
    if (refseq_graph.size() == 1) {
        // *this = *refseq_graph[targets.front()];
        // we have already done this because the first graph in the queue is this
    } else {
        // where we have multiple targets
        for (vector<string>::iterator t = targets.begin(); t != targets.end(); ++t) {
            // merge the variants into one graph
            VG& g = *refseq_graph[*t];
            combine(g);
        }
    }
    // rebuild the mapping ranks now that we've combined everything
    paths.clear_mapping_ranks();
    paths.rebuild_mapping_aux();

    if(load_phasing_paths) {
        // Trace through all the phase paths, and, where they take edges that
        // don't exist, break them. TODO: we still might get spurious phasing
        // through a deletion where the two pahsed bits but up against each
        // other.

        create_progress("dividing phasing paths", num_phasings);
        for(size_t i = 0; i < num_phasings; i++) {
            // What's the path we want to trace?
            string original_path_name = "_phase" + to_string(i);

            list<Mapping>& path_mappings = paths.get_path(original_path_name);

            // What section of this phasing do we want to be outputting?
            size_t subpath = 0;
            // Make a name for it
            string subpath_name = "_phase" + to_string(i) + "_" + to_string(subpath);

            // For each mapping, we want to be able to look at the previous
            // mapping.
            list<Mapping>::iterator prev_mapping = path_mappings.end();
            for(list<Mapping>::iterator mapping = path_mappings.begin(); mapping != path_mappings.end(); ++mapping) {
                // For each mapping in the path
                if(prev_mapping != path_mappings.end()) {
                    // We have the previous mapping and this one

                    // Make the two sides of nodes that should be connected.
                    auto s1 = NodeSide(prev_mapping->position().node_id(),
                        (prev_mapping->position().is_reverse() ? false : true));
                    auto s2 = NodeSide(mapping->position().node_id(),
                        (mapping->position().is_reverse() ? true : false));
                    // check that we always have an edge between the two nodes in the correct direction
                    if (!has_edge(s1, s2)) {
                        // We need to split onto a new subpath;
                        subpath++;
                        subpath_name = "_phase" + to_string(i) + "_" + to_string(subpath);
                    }
                }

                // Now we just drop this node onto the current subpath
                paths.append_mapping(subpath_name, *mapping);

                // Save this mapping as the prev one
                prev_mapping = mapping;
            }

            // Now delete the original full phase path.
            // This invalidates the path_mappings reference!!!
            // We use the variant that actually unthreads the path from the indexes and doesn't erase and rebuild them.
            paths.remove_path(original_path_name);

            update_progress(i);
        }
        destroy_progress();


    }

    std::function<bool(string)> all_upper = [](string s){
        //GO until [size() - 1 ] to avoid the newline char
        for (int i = 0; i < s.size() - 1; i++){
            if (!isupper(s[i])){
                return false;
            }
        }
        return true;
    };

    for_each_node([&](Node* node) {
            if (!all_upper(node->sequence())){
                cerr << "WARNING: Lower case letters found during construction" << endl;
                cerr << "Sequences may not map to this graph." << endl;
                cerr << pb2json(*node) << endl;
            }
        });

}

map<id_t, vcflib::Variant> VG::get_node_id_to_variant(vcflib::VariantCallFile vfile){
    map<id_t, vcflib::Variant> ret;
    vcflib::Variant var;

    while(vfile.getNextVariant(var)){
        long nuc = var.position;
        id_t node_id = get_node_at_nucleotide(var.sequenceName, nuc);
        ret[node_id] = var;
    }

    return ret;
}

// construct from VCF records
// --------------------------
// algorithm
// maintain a core reference path upon which we add new variants as they come
// addition procedure is the following
// find reference node overlapping our start position
// if it is already the end of a node, add the new node
// if it is not the end of a node, break it, insert edges from old->new
// go to end position of alt allele (could be the same position)
// if it already has a break, just point to the next node in line
// if it is not broken, break it and point to the next node
// add new node for alt alleles, connect to start and end node in reference path
// store the ref mapping as a property of the edges and nodes (this allows deletion edges and insertion subpaths)
//

void VG::vcf_records_to_alleles(vector<vcflib::Variant>& records,
                                map<long, vector<vcflib::VariantAllele> >& altp,
                                map<pair<long, int>, vector<bool>>* phase_visits,
                                map<pair<long, int>, vector<pair<string, int>>>* alt_allele_visits,
                                bool flat_input_vcf) {



#ifdef debug
    cerr << "Processing " << records.size() << " vcf records..." << endl;
#endif

    for (int i = 0; i < records.size(); ++i) {
        vcflib::Variant& var = records.at(i);

        // What name should we use for the variant? We need to make sure it is
        // unique, even if there are multiple variant records at the same
        // position in the VCF. Also, we don't necessarily have every variant in
        // the VCF in our records vector.
        string var_name = make_variant_id(var);

        // decompose to alts
        // This holds a map from alt or ref allele sequence to a series of VariantAlleles describing an alignment.
        map<string, vector<vcflib::VariantAllele> > alternates
            = (flat_input_vcf ? var.flatAlternates() : var.parsedAlternates());

        if(!alternates.count(var.ref)) {
            // Ref is missing, as can happen with flat construction.
            // Stick the ref in, because we need to have ref.
            alternates[var.ref].push_back(vcflib::VariantAllele(var.ref, var.ref, var.position));
        }

        // This holds a map from alt index (0 for ref) to the phase sets
        // visiting it as a bool vector. No bit vector means no visits.
        map<int, vector<bool>> alt_usages;

        if(phase_visits != nullptr) {

            // Parse out what alleles each sample uses in its phase sets at this
            // VCF record.

            // Get all the sample names in order.
            auto& sample_names = var.vcf->sampleNames;

            for(int64_t j = 0; j < sample_names.size(); j++) {
                // For every sample, see if at this variant it uses this
                // allele in one or both phase sets.

                // Grab the genotypes
                string genotype = var.getGenotype(sample_names[j]);

                // Find the phasing bar
                auto bar_pos = genotype.find('|');

                if(bar_pos == string::npos || bar_pos == 0 || bar_pos + 1 >= genotype.size()) {
                    // Not phased here, or otherwise invalid
                    continue;
                }

                if(genotype.substr(0, bar_pos) == "." || genotype.substr(bar_pos + 1) == ".") {
                    // This site is uncalled
                    continue;
                }

                // Parse out the two alt indexes.
                // TODO: complain if there are more.
                int alt1index = stoi(genotype.substr(0, bar_pos));
                int alt2index = stoi(genotype.substr(bar_pos + 1));

                if(!alt_usages.count(alt1index)) {
                    // Make a new bit vector for the alt visited by 1
                    alt_usages[alt1index] = vector<bool>(var.getNumSamples() * 2, false);
                }
                // First phase of this phase set visits here.
                alt_usages[alt1index][j * 2] = true;

                if(!alt_usages.count(alt2index)) {
                    // Make a new bit vector for the alt visited by 2
                    alt_usages[alt2index] = vector<bool>(var.getNumSamples() * 2, false);
                }
                // Second phase of this phase set visits here.
                alt_usages[alt2index][j * 2 + 1] = true;
            }
        }

        for (auto& alleles : alternates) {

            // We'll point this to a vector flagging all the phase visits to
            // this alt (which may be the ref alt), if we want to record those.
            vector<bool>* visits = nullptr;

            // What alt number is this alt? (0 for ref)
            // -1 for nothing needs to visit it and we don't care.
            int alt_number = -1;

#ifdef debug
            cerr << "Considering alt " << alleles.first << " at " << var.position << endl;
            cerr << var << endl;
#endif

            if(phase_visits != nullptr || alt_allele_visits != nullptr) {
                // We actually have visits to look for. We need to know what
                // alt number we have here.

                // We need to copy out the alt sequence to appease the vcflib API
                string alt_sequence = alleles.first;

                // What alt number are we looking at
                if(alt_sequence == var.ref) {
                    // This is the ref allele
                    alt_number = 0;
                } else {
                    // This is an alternate allele
                    alt_number = var.getAltAlleleIndex(alt_sequence) + 1;
                }

#ifdef debug
                cerr << "Alt is number " << alt_number << endl;
#endif

                if(alt_usages.count(alt_number)) {
                    // Something did indeed visit. Point the pointer at the
                    // vector describing what visited.
                    visits = &alt_usages[alt_number];
                }
            }

            for (auto& allele : alleles.second) {
                // For each of the alignment bubbles or matches, add it in as something we'll need for the graph.
                // These may overlap between alleles, and not every allele will have one at all positions.
                // In general it has to be that way, because the alleles themselves can overlap.

                // TODO: we need these to be unique but also ordered by addition
                // order. For now we just check all previous entries before
                // adding and suffer being n^2 in vcf alts per variant. We
                // should use some kind of addition-ordered set.
                int found_at = -1;
                for(int j = 0; j < altp[allele.position].size(); j++) {
                    if(altp[allele.position][j].ref == allele.ref && altp[allele.position][j].alt == allele.alt) {
                        // TODO: no equality for VariantAlleles for some reason.
                        // We already have it at this index
                        found_at = j;
                        break;
                    }
                }
                if(found_at == -1) {
                    // We need to tack this on at the end.
                    found_at = altp[allele.position].size();
                    // Add the bubble made by this part of this alt at this
                    // position.
                    altp[allele.position].push_back(allele);
                }

                if(visits != nullptr && phase_visits != nullptr) {
                    // We have to record a phase visit

                    // What position, allele index pair are we visiting when we
                    // visit this alt?
                    auto visited = make_pair(allele.position, found_at);

                    if(!phase_visits->count(visited)) {
                        // Make sure we have a vector for visits to this allele, not
                        // just this alt. It needs an entry for each phase of each sample.
                        (*phase_visits)[visited] = vector<bool>(var.getNumSamples() * 2, false);
                    }

                    for(size_t j = 0; j < visits->size(); j++) {
                        // We need to toggle on all the phase sets that visited
                        // this alt as using this allele at this position.
                        if(visits->at(j) && !(*phase_visits)[visited].at(j)) {
                            // The bit needs to be set, because all the phases
                            // visiting this alt visit this allele that appears
                            // in it.
                            (*phase_visits)[visited][j] = true;
                        }

                    }
                }

                if(alt_allele_visits != nullptr && alt_number != -1) {
                    // We have to record a visit of this alt of this variant to
                    // this VariantAllele bubble/reference patch.

                    // What position, allele index pair are we visiting when we
                    // visit this alt?
                    auto visited = make_pair(allele.position, found_at);

#ifdef debug
                    cerr << var_name << " alt " << alt_number << " visits allele #" << found_at
                        << " at position " << allele.position << " of " << allele.ref << " -> " << allele.alt << endl;
#endif

                    // Say we visit this allele as part of this alt of this variant.
                    (*alt_allele_visits)[visited].push_back(make_pair(var_name, alt_number));
                }

            }
        }
    }
}

void VG::slice_alleles(map<long, vector<vcflib::VariantAllele> >& altp,
                       int start_pos,
                       int stop_pos,
                       int max_node_size) {

    // Slice up only the *reference*. Leaves the actual alt sequences alone.
    // Does *not* divide up the alt alleles into multiple pieces, despite its
    // name.

    auto enforce_node_size_limit =
        [this, max_node_size, &altp]
        (int curr_pos, int& last_pos) {
        int last_ref_size = curr_pos - last_pos;
        update_progress(last_pos);
        if (max_node_size && last_ref_size > max_node_size) {
            int div = 2;
            while (last_ref_size/div > max_node_size) {
                ++div;
            }
            int segment_size = last_ref_size/div;
            int i = 0;
            while (last_pos + i < curr_pos) {
                altp[last_pos+i];  // empty cut
                i += segment_size;
                update_progress(last_pos + i);
            }
        }
    };

    if (max_node_size > 0) {
        create_progress("enforcing node size limit ", (altp.empty()? 0 : altp.rbegin()->first));
        // break apart big nodes
        int last_pos = start_pos;
        for (auto& position : altp) {
            auto& alleles = position.second;
            enforce_node_size_limit(position.first, last_pos);
            for (auto& allele : alleles) {
                // cut the last reference sequence into bite-sized pieces
                last_pos = max(position.first + allele.ref.size(), (long unsigned int) last_pos);
            }
        }
        enforce_node_size_limit(stop_pos, last_pos);
        destroy_progress();
    }

}

}


