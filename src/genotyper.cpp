#include "genotyper.hpp"
#include "bubbles.hpp"

namespace vg {

using namespace std;

/**
 * Turn the given path (which must be a thread) into an allele. Drops the first
 * and last mappings and looks up the sequences for the nodes of the others.
 */
string allele_to_string(VG& graph, const Path& allele) {
    stringstream stream;
    
    for(size_t i = 1; i < allele.mapping_size() - 1; i++) {
        // Get the sequence for each node
        string node_string = graph.get_node(allele.mapping(i).position().node_id())->sequence();
        
        if(allele.mapping(i).position().is_reverse()) {
            // Flip it
            node_string = reverse_complement(node_string);
        }
        // Add it to the stream
        stream << node_string;
    }
    
    return stream.str();
}

double Genotyper::phred_to_prob(int phred) {
    return pow(10, -((double)phred) / 10);
}

int Genotyper::prob_to_phred(double prob) {
    return round(-10.0 * log10(prob));
}

int Genotyper::alignment_qual_score(const Alignment& alignment) {
    if(alignment.quality().empty()) {
        // Special case: qualities not given. Assume something vaguely sane so
        // we can genotype without quality.
#ifdef debug
        cerr << "No base qualities. Assuming default quality of " << default_sequence_quality << endl;
#endif
        return default_sequence_quality;
    }
    
    double total = 0;
    for(size_t i = 0; i < alignment.quality().size(); i++) {
#ifdef debug
        cerr << "Quality: " << (int)alignment.quality()[i] << endl;
#endif
        total += alignment.quality()[i];
    }
#ifdef debug
    cerr << "Count: " << alignment.quality().size() << endl;
#endif
    // Make the total now actually be an average
    total /= alignment.quality().size();
    return round(total);
}

map<pair<NodeTraversal, NodeTraversal>, set<id_t>> Genotyper::find_sites(VG& graph) {

    // Set up our output map
    map<pair<NodeTraversal, NodeTraversal>, set<id_t>> to_return;

    // Unfold the graph
    // Copy the graph and unfold the copy. We need to hold this translation
    // table from new node ID to old node and relative orientation.
    map<vg::id_t, pair<vg::id_t, bool>> unfold_translation;
    auto transformed = graph.unfold(unfold_max_length, unfold_translation);
    
    // Fix up any doubly reversed edges
    transformed.flip_doubly_reversed_edges();

    // Now dagify the graph. We need to hold this translation table from new
    // node ID to old node and relative orientation.
    map<vg::id_t, pair<vg::id_t, bool>> dag_translation;
    transformed = transformed.dagify(dagify_steps, dag_translation);
    
    // Compose the complete translation
    map<vg::id_t, pair<vg::id_t, bool>> overall_translation = transformed.overlay_node_translations(dag_translation, unfold_translation);
    dag_translation.clear();
    unfold_translation.clear();
    
    // Find the superbubbles in the DAG
    map<pair<id_t, id_t>, vector<id_t>> superbubbles = vg::superbubbles(transformed);
    
    for(auto& superbubble : superbubbles) {
        
        // Translate the superbubble coordinates into NodeTraversals
        auto& start_translation = overall_translation[superbubble.first.first];
        NodeTraversal start(graph.get_node(start_translation.first), start_translation.second);
        
        auto& end_translation = overall_translation[superbubble.first.second];
        NodeTraversal end(graph.get_node(end_translation.first), end_translation.second);
        
        // Find the vector we want all the nodes in
        auto& superbubble_nodes = to_return[make_pair(start, end)];
        
        for(auto id : superbubble.second) {
            // Translate each ID and put it in the set
            superbubble_nodes.insert(overall_translation[id].first);
        }
    }

    // Give back the collections of involved nodes by start and end
    return to_return;    
    
}

map<pair<NodeTraversal, NodeTraversal>, vector<id_t>> Genotyper::find_sites_with_cactus(VG& graph) {

    // Set up our output map
    map<pair<NodeTraversal, NodeTraversal>, vector<id_t>> to_return;

    // todo: use deomposition instead of converting tree into flat structure
    BubbleTree bubble_tree = cactusbubble_tree(graph);

    bubble_tree.for_each_preorder([&](Bubble& bubble) {
            // cut root to be consistent with superbubbles()
            if (bubble.start != bubble_tree.root->bubble.start ||
                bubble.end != bubble_tree.root->bubble.end) {
                vector<id_t> nodes = bubble.contents;
                // include endpoints to be consistent with superbubbles()
                nodes.push_back(bubble.start.node);
                nodes.push_back(bubble.end.node);
                // keep sorted by id to be consistent with superbubbles()
                sort(nodes.begin(), nodes.end());
                NodeTraversal start(graph.get_node(bubble.start.node), bubble.start.is_end);
                NodeTraversal end(graph.get_node(bubble.end.node), !bubble.start.is_end);
                to_return[minmax(start, end)] = nodes;
            }
        });

    return to_return;
}

vector<Path> Genotyper::get_paths_through_site(VG& graph, NodeTraversal start, NodeTraversal end) {
    // We're going to emit traversals supported by any paths in the graph.
    
    // Put all our subpaths in here to deduplicate them by sequence they spell out.
    map<string, list<NodeTraversal>> results;

#ifdef debug
    cerr << "Looking for paths between " << start << " and " << end << endl;
#endif
    
    if(graph.paths.has_node_mapping(start.node) && graph.paths.has_node_mapping(end.node)) {
        // If we have some paths that visit both ends (in some orientation)
        
        // Get all the mappings to the end node, by path name
        auto& endmappings_by_name = graph.paths.get_node_mapping(end.node);
        
        for(auto& name_and_mappings : graph.paths.get_node_mapping(start.node)) {
            // Go through the paths that visit the start node
            
            if(!endmappings_by_name.count(name_and_mappings.first)) {
                // No path by this name has any mappings to the end node. Skip
                // it early.
                continue;
            }
            
            for(auto* mapping : name_and_mappings.second) {
                // Start at each mapping in the appropriate orientation
                
#ifdef debug
                cerr << "Trying mapping of read " << name_and_mappings.first << endl;
#endif
                
                // How many times have we gone to the next mapping looking for a
                // mapping to the end node in the right orientation?
                size_t traversal_count = 0;
                
                // Do we want to go left (true) or right (false) from this
                // mapping? If start is a forward traversal and we found a
                // forward mapping, we go right. If either is backward we go
                // left, and if both are backward we go right again.
                bool traversal_direction = mapping->position().is_reverse() != start.backward;
                
                // What orientation would we want to find the end node in? If
                // we're traveling backward, we expect to find it in the
                // opposite direction to the one we were given.
                bool expected_end_orientation = end.backward != traversal_direction;
                
                // We're going to fill in this list with traversals.
                list<NodeTraversal> path_traversed;
                
                // And we're going to fill this with the sequence
                stringstream allele_stream;
                
                while(mapping != nullptr && traversal_count < max_path_search_steps) {
                    // Traverse along until we hit the end traversal or take too
                    // many steps
                    
#ifdef debug
                    cerr << "\tTraversing " << pb2json(*mapping) << endl;
#endif
                    
                    // Say we visit this node along the path, in this orientation
                    path_traversed.push_back(NodeTraversal(graph.get_node(mapping->position().node_id()), mapping->position().is_reverse() != traversal_direction));
                    
                    // Stick the sequence of the node (appropriately oriented) in the stream for the allele sequence
                    string seq = graph.get_node(mapping->position().node_id())->sequence();
                    allele_stream << (path_traversed.back().backward ? reverse_complement(seq) : seq);
                    
                    if(mapping->position().node_id() == end.node->id() && mapping->position().is_reverse() == expected_end_orientation) {
                        // We have stumbled upon the end node in the orientation we wanted it in.
                        
                        // Keep this subpath if it's for a new string
                        // Also ends up replacing anything that was there for the string before
                        results[allele_stream.str()] = path_traversed;
                        // Then try the next embedded path
                        break;
                    }
                    
                    // Otherwise just move to the right (or left)
                    if(traversal_direction) {
                        // We're going backwards
                        mapping = graph.paths.traverse_left(mapping);
                    } else {
                        // We're going forwards
                        mapping = graph.paths.traverse_right(mapping);
                    }
                    // Tick the counter so we don't go really far on long paths.
                    traversal_count++;
                    
                }
                
                
            }
        }
        
    }
    
    // Now convert the unique results into Paths
    vector<Path> toReturn;
    
    for(auto& string_and_traversals : results) {
        // Convert each list of node traversals to a proper path and add it to our collection of paths to return
        toReturn.push_back(path_from_node_traversals(string_and_traversals.second));
    }
    
    return toReturn;
}

map<Alignment*, vector<double>> Genotyper::get_affinities(VG& graph, const map<string, Alignment*>& reads_by_name,
    set<id_t>& superbubble_contents, vector<Path>& superbubble_paths) {
    
    // Grab our thread ID, which determines which aligner we get.
    int tid = omp_get_thread_num();
    
    // We're going to build this up gradually, appending to all the vectors.
    map<Alignment*, vector<double>> toReturn;
    
    // What reads are relevant to this superbubble?
    set<string> relevant_read_names;
    
#ifdef debug
    cerr << "Superbubble contains " << superbubble_contents.size() << " nodes" << endl;
#endif

    assert(!superbubble_paths.empty());
    assert(superbubble_paths.front().mapping_size() >= 2);

    // Find the start and end nodes of the superbubble
    id_t superbubble_start = superbubble_paths.front().mapping(0).position().node_id();
    id_t superbubble_end = superbubble_paths.front().mapping(superbubble_paths.front().mapping_size() - 1).position().node_id();
    
    for(auto id : superbubble_contents) {
        // For every node in the superbubble, what paths visit it?
        auto& mappings_by_name = graph.paths.get_node_mapping(id);
        for(auto& name_and_mappings : mappings_by_name) {
            // For each path visiting the node
            if(reads_by_name.count(name_and_mappings.first)) {
                // This path is a read, so add the name to the set if it's not
                // there already
                relevant_read_names.insert(name_and_mappings.first);
            }    
        }
    }
    
    // What IDs are visited by these reads?
    set<id_t> relevant_ids;
    
    for(auto& name : relevant_read_names) {
        // Get the mappings for each read
        auto& mappings = graph.paths.get_path(name);
        for(auto& mapping : mappings) {
            // Add in all the nodes that are visited
            relevant_ids.insert(mapping.position().node_id());
        }
    }
    
    for(auto id : superbubble_contents) {
        // Throw out all the IDs that are also used in the superbubble itself
        relevant_ids.erase(id);
    }
    
#ifdef debug
    cerr << relevant_read_names.size() << " reads visit an additional " << relevant_ids.size() << " nodes" << endl;
#endif
    
    // Make a vg graph with all the nodes used by the reads relevant to the
    // superbubble, but outside the superbubble itself.
    VG surrounding;
    for(auto id : relevant_ids) {
        // Add each node and its edges to the new graph. Ignore dangling edges.
        surrounding.add_node(*graph.get_node(id));
        surrounding.add_edges(graph.edges_of(graph.get_node(id)));
    }
    
    for(auto& path : superbubble_paths) {
        // Now for each superbubble path, make a copy of that graph with it in
        VG allele_graph(surrounding);
        
        for(size_t i = 0; i < path.mapping_size(); i++) {
            // Add in every node on the path to the new allele graph
            id_t id = path.mapping(i).position().node_id();
            allele_graph.add_node(*graph.get_node(id));
            // And its edges
            allele_graph.add_edges(graph.edges_of(graph.get_node(id)));
        }
        
        // Get rid of dangling edges
        allele_graph.remove_orphan_edges();
        
#ifdef debug
        #pragma omp critical (cerr)
        cerr << "Align to " << pb2json(allele_graph.graph) << endl;
#endif
        
        for(auto& name : relevant_read_names) {
            // For every read that touched the superbubble, grab its original
            // Alignment pointer.
            Alignment* read = reads_by_name.at(name);
            
            // Look to make sure it touches more than one node actually in the
            // superbubble, or a non-start, non-end node. If it just touches the
            // start or just touches the end, it can't be informative.
            set<id_t> touched_set;
            // Will this read be informative?
            bool informative = false;            
            for(size_t i = 0; i < read->path().mapping_size(); i++) {
                // Look at every node the read touches
                id_t touched = read->path().mapping(i).position().node_id();
                if(superbubble_contents.count(touched)) {
                    // If it's in the superbubble, keep it
                    touched_set.insert(touched);
                }
            }
            
            if(touched_set.size() >= 2) {
                // We touch both the start and end, or an internal node.
                informative = true;
            } else {
                // Throw out the start and end nodes, if we touched them.
                touched_set.erase(superbubble_start);
                touched_set.erase(superbubble_end);
                if(!touched_set.empty()) {
                    // We touch an internal node
                    informative = true;
                }
            }
            
            if(!informative) {
                // We only touch one of the start and end nodes, and can say nothing about the superbubble. Try the next read.
                continue;
            }
            
            // If we get here, we know this read is informative as to the internal status of this superbubble.
            Alignment aligned;
            if(read->sequence().size() == read->quality().size()) {
                // Re-align a copy to this graph (using quality-adjusted alignment).
                // TODO: actually use quality-adjusted alignment
                aligned = allele_graph.align(*read);
            } else {
                // If we don't have the right number of quality scores, use un-adjusted alignment instead.
                aligned = allele_graph.align(*read);
            }
            
#ifdef debug
            #pragma omp critical (cerr)
            cerr << "\t" << pb2json(aligned) << endl;
#endif
            
            // Grab the identity and save it for this read and superbubble path
            toReturn[read].push_back(aligned.identity());
            
        }
    }
    
    // After scoring all the reads against all the versions of the superbubble,
    // return the affinities
    return toReturn;
}

double Genotyper::get_genotype_probability(const vector<int>& genotype, const vector<pair<Alignment, vector<bool>>> alignment_consistency) {
    // For each genotype, calculate P(observed reads | genotype) as P(all reads
    // that don't support an allele from the genotype are mismapped or
    // miscalled) * P(all reads that do support alleles from the genotype ended
    // up apportioned across the alleles as they are)
    
    // This works out to the product over all reads that don't support either
    // alleles of 1 - ((1 - MAPQ) * (1 - P(bases called wrong)), times the
    // product over all the reads that do support one of the alleles of P(that
    // allele picked out of the one or two available).
    
    // TODO: handle contamination like Freebayes
    
    // TODO: split out a function to get P(obs | genotype) for a single genotype
    
    // This is the probability that all reads that don't support either allele in this genotype are wrong.
    double all_non_supporting_wrong = 1;
    
    // This is the probability that all the reads that do support alleles in this genotype were drawn from the genotype they support.
    double all_supporting_drawn = 1;
    
#ifdef debug
    if(genotype.size() > 1) {
        cerr << "Calculating P(a" << genotype[0] << "/a" << genotype[1] << ")" << endl;
    }
#endif
    
    for(auto& read_and_consistency : alignment_consistency) {
        // For each read, work out if it supports a genotype we have or not.
        
        // Split out the alignment from its consistency flags
        auto& read = read_and_consistency.first;
        auto& consistency = read_and_consistency.second;
        
        // How many of the alleles in our genotype is it consistent with?
        int consistent_alleles = 0;
        for(int allele : genotype) {
            if(consistency.at(allele)) {
                // We're consistent with this allele
                consistent_alleles++;
            }
        }
        
        auto read_qual = alignment_qual_score(read);
        
#ifdef debug
        cerr << "Read (qual score " << read_qual << ") consistent with " << consistent_alleles << " genotype alleles observed." << endl;
#endif
        
        if(consistent_alleles == 0) {
            // This read is inconsistent with all the alleles in the genotype,
            // so, given the genotype, the read must be sequenced or mapped
            // wrong.
            
            double prob_wrong;
            if(use_mapq) {
                // Compute P(mapped wrong or called wrong) = P(not (mapped right and called right)) = P(not (not mapped wrong and not called wrong))
                prob_wrong = (1.0 - (1.0 - phred_to_prob(read.mapping_quality())) * (1.0 - phred_to_prob(read_qual)));
            } else {
                // Compute P(called wrong).
                prob_wrong = phred_to_prob(read_qual);
            }
#ifdef debug
            cerr << "P(wrong) = " << prob_wrong << endl;
#endif
            all_non_supporting_wrong *= prob_wrong;
        } else {
            // This read is consistent with some of the alleles in the genotype,
            // so we must have drawn one of those alleles when sequencing.
            
            // So multiply in the probability that we hit one of those alleles
            double prob_drawn = ((double) consistent_alleles) / genotype.size();
            
#ifdef debug
            cerr << "P(drawn) = " << prob_drawn << endl;
#endif
            
            all_supporting_drawn *= prob_drawn;
        }
        
    }
    
    // Now we've looked at all the reads, so AND everything together
    return all_non_supporting_wrong * all_supporting_drawn;
}

/**
 * Get all the quality values in the alignment between the start and end nodes
 * of a site. Handles alignments that enter the site from the end, and
 * alignments that never make it through the site.
 *
 * If we run out of qualities, or qualities aren't present, returns no
 * qualities.
 *
 * If an alignment goes through the site multipe times, we get all the qualities from when it is in the site.
 */
string get_qualities_in_site(VG& graph, const Alignment& alignment, NodeTraversal start, NodeTraversal end) {
    // We'll fill this in.
    stringstream to_return;
    
    // Are we currently in the site?
    bool in_site = false;
    // What NodeTraversal do we need to see to leave?
    NodeTraversal expected;
    
    // Where are we in the quality string?
    size_t quality_pos = 0;
    
    for(size_t i = 0; i < alignment.path().mapping_size(); i++) {
        // For every mapping in the path in order
        auto& mapping = alignment.path().mapping(i);
        
        // What NodeTraversal is this?
        NodeTraversal traversal(graph.get_node(mapping.position().node_id()), mapping.position().is_reverse());
        
        if(!in_site) {
            // If we aren't in the site, we may be entering
            if(traversal == start) {
                // We entered through the start
                in_site = true;
                // We'll leave at the end
                expected = end;
            } else if(traversal.node == end.node && traversal.backward != end.backward) {
                // We entered through the end
                in_site = true;
                // We'll leave when we hit the start in reverse
                expected = start;
                expected.backward = !expected.backward;
            }
        }
        
        for(size_t j = 0; j < mapping.edit_size(); j++) {
            // For every edit
            auto& edit = mapping.edit(j);
        
            if(in_site) {
                for(size_t k = 0; k < edit.to_length(); k++) {
                    // Take each quality value from the edit and add it to our collection to return
                    if(quality_pos >= alignment.quality().size()) {
                        // If we've run out of quality values, give back no
                        // qualities, because base qualities aren't really being
                        // used.
                        return "";
                    }
                    to_return << (char)alignment.quality().at(quality_pos);
                    quality_pos++;
                }
            } else {
                // Skip this edit's qualities 
                quality_pos += edit.to_length();
            }
        }
        
        if(in_site && traversal == expected) {
            // This was the node we were supposed to leave the site at.
            in_site = false;
        }
    }
    
    return to_return.str();
    
}

#define debug
Genotype Genotyper::get_genotype(VG& graph, const vector<Path>& superbubble_paths, const map<Alignment*, vector<double>>& affinities) {
    
    // Freebayes way (improved with multi-support):
    
    if(superbubble_paths.empty()) {
        // No paths! Nothing to say!
        throw runtime_error("No paths through superbubble! Can't genotype!");
    }
    
    if(superbubble_paths.front().mapping_size() < 2) {
        throw runtime_error("Too few mappings in superbubble path!");
    }
    
    // Get Mappings bounding the superbubble
    const Mapping& start_mapping = superbubble_paths.front().mapping(0);
    const Mapping& end_mapping = superbubble_paths.front().mapping(superbubble_paths.front().mapping_size() - 1);
    
    // And NodeTraversals. TODO: pass them in instead?
    NodeTraversal start(graph.get_node(start_mapping.position().node_id()), start_mapping.position().is_reverse());
    NodeTraversal end(graph.get_node(end_mapping.position().node_id()), end_mapping.position().is_reverse());
    
#ifdef debug
    cerr << "Looking between " << start << " and " << end << endl;
#endif
    
    // We'll fill this in with the trimmed alignments and their consistency-with-alleles flags.
    vector<pair<Alignment, vector<bool>>> alignment_consistency;
    
    // We fill this in with totals
    vector<int> reads_consistent_with_allele(superbubble_paths.size(), 0);
    
    // We'll store affinities by read name and allele here, for printing later.
    map<string, vector<double>> debug_affinities;
    
    for(auto& alignment_and_affinities : affinities) {
        // We need to clip down to just the important quality values        
        Alignment trimmed = *alignment_and_affinities.first;
        
        // Trim down qualities to just those in the superbubble
        auto trimmed_qualities = get_qualities_in_site(graph, trimmed, start, end);
        
        // Stick them back in the now bogus-ified Alignment object.
        trimmed.set_quality(trimmed_qualities);
        
        // Hide all the affinities where we can pull them later
        debug_affinities[trimmed.name()] = alignment_and_affinities.second;
        
        // Calculate whether each affinity is enough to justify support.
        // We say only the maximal affinities are support.
        double max_affinity = *max_element(alignment_and_affinities.second.begin(), alignment_and_affinities.second.end());
        
        // Count up consistency
        size_t consistent = 0;
        
        // Flag alleles we're consistent with
        vector<bool> support_flags;
        for(auto affinity : alignment_and_affinities.second) {
            // Only maximal affinities represent consistency with an allele
            if(affinity == max_affinity) {
                // Read consistent with allele
                support_flags.push_back(true);
                // Count total alleles read is consistent with
                consistent++;
                // Count total reads consistent with allele
                reads_consistent_with_allele[support_flags.size() - 1]++;
            } else {
                // Read not consistent with allele
                support_flags.push_back(false);
            }
        }
        
        // Save the trimmed alignment and its consistency flags
        alignment_consistency.emplace_back(std::move(trimmed), std::move(support_flags));
        
    }
    
#ifdef debug
    for(int i = 0; i < reads_consistent_with_allele.size(); i++) {
        // Build a useful name for the allele
        stringstream allele_name;
        for(size_t j = 0; j < superbubble_paths[i].mapping_size(); j++) {
            allele_name << superbubble_paths[i].mapping(j).position().node_id() << ",";
        }
        cerr << "a" << i << "(" << allele_name.str() << "): " << reads_consistent_with_allele[i] << "/" << affinities.size() << " reads consistent" << endl;
        for(auto& read_and_consistency : alignment_consistency) {
            if(read_and_consistency.second[i]) {
                // Dump all the consistent reads
                cerr << "\t" << read_and_consistency.first.sequence() << " " << debug_affinities[read_and_consistency.first.name()][i] << endl;
            }
        }
    }
#endif

    // We'll go through all the genotypes and find the best
    vector<int> best_genotype;
    double best_probability = 0;

    for(int allele1 = 0; allele1 < superbubble_paths.size(); allele1++) {
        // For each first allele in the genotype
        for(int allele2 = 0; allele2 <= allele1; allele2++) {
            // For each second allele so we get all order-independent combinations
            
            // Make the combo
            vector<int> genotype = {allele1, allele2};
            
            // Compute the probability of the genotype
            double probability = get_genotype_probability(genotype, alignment_consistency);
            
#ifdef debug
            cerr << "P(obs | a" << allele1 << "/a" << allele2 << ") = " << probability << endl;
#endif
            
            if(probability > best_probability || best_genotype.empty()) {
                // This is the best genotype!
                best_probability = probability;
                best_genotype = genotype;
            }
        }
    }
    
    // Now compute the vector of ints into a real Genotype object
    Genotype genotype;
    
    // We only want to put each allele once
    set<int> used_alleles;
    
    for(int allele : best_genotype) {
        if(!used_alleles.count(allele)) {
            // For each allele we haven't reported, report it.
            used_alleles.insert(allele);
            *genotype.add_allele() = superbubble_paths.at(allele);
            Support* supoport = genotype.add_support();
            // TODO: turn overall genotype probability into quality and stick it in
        }
    }
    
    return genotype;
}
#undef debug

/**
 * Create the reference allele for an empty vcflib Variant, since apaprently
 * there's no method for that already. Must be called before any alt alleles are
 * added.
 */
void create_ref_allele(vcflib::Variant& variant, const std::string& allele) {
    // Set the ref allele
    variant.ref = allele;
    
    for(size_t i = 0; i < variant.ref.size(); i++) {
        // Look at all the bases
        if(variant.ref[i] != 'A' && variant.ref[i] != 'C' && variant.ref[i] != 'G' && variant.ref[i] != 'T') {
            // Correct anything bogus (like "X") to N
            variant.ref[i] = 'N';
        }
    }
    
    // Make it 0 in the alleles-by-index list
    variant.alleles.push_back(allele);
    // Build the reciprocal index-by-allele mapping
    variant.updateAlleleIndexes();
}

/**
 * Add a new alt allele to a vcflib Variant, since apaprently there's no method
 * for that already.
 *
 * If that allele already exists in the variant, does not add it again.
 *
 * Retuerns the allele number (0, 1, 2, etc.) corresponding to the given allele
 * string in the given variant. 
 */
int add_alt_allele(vcflib::Variant& variant, const std::string& allele) {
    // Copy the allele so we can throw out bad characters
    std::string fixed(allele);
    
    for(size_t i = 0; i < fixed.size(); i++) {
        // Look at all the bases
        if(fixed[i] != 'A' && fixed[i] != 'C' && fixed[i] != 'G' && fixed[i] != 'T') {
            // Correct anything bogus (like "X") to N
            fixed[i] = 'N';
        }
    }
    
    for(int i = 0; i < variant.alleles.size(); i++) {
        if(variant.alleles[i] == fixed) {
            // Already exists
            return i;
        }
    }

    // Add it as an alt
    variant.alt.push_back(fixed);
    // Make it next in the alleles-by-index list
    variant.alleles.push_back(fixed);
    // Build the reciprocal index-by-allele mapping
    variant.updateAlleleIndexes();

    // We added it in at the end
    return variant.alleles.size() - 1;
}

void Genotyper::write_vcf_header(std::ostream& stream, const std::string& sample_name, const std::string& contig_name, size_t contig_size) {
    stream << "##fileformat=VCFv4.2" << std::endl;
    stream << "##ALT=<ID=NON_REF,Description=\"Represents any possible alternative allele at this location\">" << std::endl;
    stream << "##INFO=<ID=XREF,Number=0,Type=Flag,Description=\"Present in original graph\">" << std::endl;
    stream << "##INFO=<ID=XSEE,Number=.,Type=String,Description=\"Original graph node:offset cross-references\">" << std::endl;
    stream << "##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Total Depth\">" << std::endl;
    stream << "##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Read Depth\">" << std::endl;
    stream << "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">" << std::endl;
    stream << "##FORMAT=<ID=AD,Number=.,Type=Integer,Description=\"Allelic depths for the ref and alt alleles in the order listed\">" << std::endl;
    stream << "##FORMAT=<ID=SB,Number=4,Type=Integer,Description=\"Forward and reverse support for ref and alt alleles.\">" << std::endl;
    // We need this field to stratify on for VCF comparison. The info is in SB but vcfeval can't pull it out
    stream << "##FORMAT=<ID=XAAD,Number=1,Type=Integer,Description=\"Alt allele read count.\">" << std::endl;
    if(!contig_name.empty()) {
        // Announce the contig as well.
        stream << "##contig=<ID=" << contig_name << ",length=" << contig_size << ">" << std::endl;
    }
    stream << "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" << sample_name << std::endl;
}

vcflib::VariantCallFile* Genotyper::start_vcf(std::ostream& stream, const ReferenceIndex& index, const string& ref_path_name) {
    // Generate a vcf header. We can't make Variant records without a
    // VariantCallFile, because the variants need to know which of their
    // available info fields or whatever are defined in the file's header, so
    // they know what to output.
    // Handle length override if specified.
    std::stringstream headerStream;
    write_vcf_header(headerStream, "SAMPLE", ref_path_name, index.sequence.size());
    
    // Load the headers into a new VCF file object
    vcflib::VariantCallFile* vcf = new vcflib::VariantCallFile();
    std::string headerString = headerStream.str();
    assert(vcf->openForOutput(headerString));
    
    // Spit out the header
    stream << headerStream.str();
    
    // Give back the created VCF
    return vcf;
}

vector<vcflib::Variant> Genotyper::genotype_to_variant(VG& graph, const ReferenceIndex& index, vcflib::VariantCallFile& vcf, const Genotype& genotype) {
    // Make a vector to fill in
    vector<vcflib::Variant> toReturn;
    
    // Make a new variant
    vcflib::Variant variant;
    // Attach it to the VCF
    variant.setVariantCallFile(vcf);
    // Crib the quality
    variant.quality = 0;
    
    // Make sure we have stuff
    if(genotype.allele_size() == 0) {
        throw runtime_error("Can't turn an empty genotype into VCF");
    }
    if(genotype.allele(0).mapping_size() == 0) {
        throw runtime_error("Can't turn an empty allele into VCF");
    }
    
    // All the paths start with the same beginning and end nodes, so we can just
    // look at those to get the reference interval.
    
    auto first_id = genotype.allele(0).mapping(0).position().node_id();
    auto last_id = genotype.allele(0).mapping(genotype.allele(0).mapping_size() - 1).position().node_id();
    
    if(!index.byId.count(first_id) || !index.byId.count(last_id)) {
        // We need to be anchored to the primary path to make a variant
        cerr << "Warning: Superbubble endpoints not on reference!" << endl;
        // If not return no variant
        return toReturn;
    }
    
    // The position we have stored for this start node is the first
    // position along the reference at which it occurs. Our bubble
    // goes forward in the reference, so we must come out of the
    // opposite end of the node from the one we have stored.
    auto referenceIntervalStart = index.byId.at(first_id).first + graph.get_node(first_id)->sequence().size();
    
    // The position we have stored for the end node is the first
    // position it occurs in the reference, and we know we go into
    // it in a reference-concordant direction, so we must have our
    // past-the-end position right there.
    auto referenceIntervalPastEnd = index.byId.at(last_id).first;
    
    // TODO: figure out how to handle superbubbles that come up backwards
    // relative to the primary reference.
    assert(referenceIntervalStart <= referenceIntervalPastEnd);
    
    // Get the string for the reference allele
    string refString = index.sequence.substr(referenceIntervalStart, referenceIntervalPastEnd - referenceIntervalStart);
    
    // And make strings for all the genotype alleles
    vector<string> allele_strings;
    
    for(size_t i = 0; i < genotype.allele_size(); i++) {
        // Get the string for each allele
        allele_strings.push_back(allele_to_string(graph, genotype.allele(i)));
    }
    
    // See if any alleles are empty
    bool empty_alleles = refString.empty();
    for(auto& allele : allele_strings) {
        if(allele == "") {
            empty_alleles = true;
        }
    }
    
    // Fix them up
    if(empty_alleles) {
        // Grab the character before our site
        string prefix = index.sequence.substr(referenceIntervalStart - 1, 1);
        for(auto& allele : allele_strings) {
            // Prepend it to every allele
            allele = prefix + allele;
        }
        // Also prepend it to the reference string
        refString = prefix + refString;
        
        // Budge the variant over
        referenceIntervalStart--;
    }
    
    // Make the ref allele
    create_ref_allele(variant, refString);
    
    // We're going to put all the allele indexes called as present in here.
    vector<int> called_alleles;
    
    for(size_t i = 0; i < genotype.allele_size(); i++) {
        // For each allele
        
        // Add it/find its number if it already exists (i.e. is the ref)
        int allele_number = add_alt_allele(variant, allele_strings[i]);
        
        // Remember we called a copy of it
        called_alleles.push_back(allele_number);
        
        if(i < genotype.support_size()) {
            // Add the quality
            variant.quality += genotype.support(i).quality();
        }
    }
    
    assert(!called_alleles.empty());
    
    while(called_alleles.size() < 2) {
        // Extend single-allele sites to be diploid.
        called_alleles.push_back(called_alleles.front());
    }
    
    // Compose the genotype
    variant.format.push_back("GT");
    auto& genotype_out = variant.samples["SAMPLE"]["GT"];
    genotype_out.push_back(to_string(called_alleles[0]) + "/" + to_string(called_alleles[1]));
    
    // Set the variant position (now that we have budged it left if necessary
    variant.position = referenceIntervalStart + 1;
    
    // Return the variant, since we managed to make it
    toReturn.push_back(variant);
    return toReturn;
    
}

ReferenceIndex::ReferenceIndex(vg::VG& vg, std::string refPathName) {
    // Make sure the reference path is present
    assert(vg.paths.has_path(refPathName));
    
    // We're also going to build the reference sequence string
    std::stringstream refSeqStream;
    
    // What base are we at in the reference
    size_t referenceBase = 0;
    
    // What was the last rank? Ranks must always go up.
    int64_t lastRank = -1;
    
    for(auto mapping : vg.paths.get_path(refPathName)) {
    
        if(!byId.count(mapping.position().node_id())) {
            // This is the first time we have visited this node in the reference
            // path.
            
            // Add in a mapping.
            byId[mapping.position().node_id()] = 
                std::make_pair(referenceBase, mapping.position().is_reverse());
#ifdef debug
            std::cerr << "Node " << mapping.position().node_id() << " rank " << mapping.rank()
                << " starts at base " << referenceBase << " with "
                << vg.get_node(mapping.position().node_id())->sequence() << std::endl;
#endif
            
            // Make sure ranks are monotonically increasing along the path.
            assert(mapping.rank() > lastRank);
            lastRank = mapping.rank();
        }
        
        // Find the node's sequence
        std::string sequence = vg.get_node(mapping.position().node_id())->sequence();
        
        while(referenceBase == 0 && sequence.size() > 0 &&
            (sequence[0] != 'A' && sequence[0] != 'T' && sequence[0] != 'C' &&
            sequence[0] != 'G' && sequence[0] != 'N')) {
            
            // If the path leads with invalid characters (like "X"), throw them
            // out when computing reference path positions.
            
            // TODO: this is a hack to deal with the debruijn-brca1-k63 graph,
            // which leads with an X.
            
            std::cerr << "Warning: dropping invalid leading character "
                << sequence[0] << " from node " << mapping.position().node_id()
                << std::endl;
                
            sequence.erase(sequence.begin());
        }
        
        if(mapping.position().is_reverse()) {
            // Put the reverse sequence in the reference path
            refSeqStream << vg::reverse_complement(sequence);
        } else {
            // Put the forward sequence in the reference path
            refSeqStream << sequence;
        }
            
        // Say that this node appears here along the reference in this
        // orientation.
        byStart[referenceBase] = vg::NodeTraversal(
            vg.get_node(mapping.position().node_id()), mapping.position().is_reverse()); 
            
        // Whether we found the right place for this node in the reference or
        // not, we still need to advance along the reference path. We assume the
        // whole node (except any leading bogus characters) is included in the
        // path (since it sort of has to be, syntactically, unless it's the
        // first or last node).
        referenceBase += sequence.size();
        
        // TODO: handle leading bogus characters in calls on the first node.
    }
    
    // Create the actual reference sequence we will use
    sequence = refSeqStream.str();
    
    // Announce progress.
    std::cerr << "Traced " << referenceBase << " bp reference path " << refPathName << "." << std::endl;
    
    if(sequence.size() < 100) {
        std::cerr << "Reference sequence: " << sequence << std::endl;
    }
}


}


