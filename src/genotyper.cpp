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

int Genotyper::alignment_score(const Alignment& alignment) {
    double total = 0;
    for(size_t i = 0; i < alignment.quality().size(); i++) {
        total += alignment.quality()[i];
    }
    // Make the total now actually be an average
    total /= alignment.quality().size();
    return round(total);
}

map<pair<NodeTraversal, NodeTraversal>, vector<id_t>> Genotyper::find_sites(VG& graph) {

    // Set up our output map
    map<pair<NodeTraversal, NodeTraversal>, vector<id_t>> to_return;

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
            // Translate each ID and put it in the vector
            superbubble_nodes.push_back(overall_translation[id].first);
        }
    }

    // Give back the collections of involved nodes by start and end
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
            
            // TODO: instead of walking, use cut_path or cut_alignment
            
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
                    allele_stream << path_traversed.back().backward ? reverse_complement(seq) : seq;
                    
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
    vector<id_t>& superbubble_contents, vector<Path>& superbubble_paths) {
    
    // Grab our thread ID, which determines which aligner we get.
    int tid = omp_get_thread_num();
    
    // We're going to build this up gradually, appending to all the vectors.
    map<Alignment*, vector<double>> toReturn;
    
    // What reads are relevant to this superbubble?
    set<string> relevant_read_names;
    
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
            // Add in every node on the path
            id_t id = path.mapping(i).position().node_id();
            surrounding.add_node(*graph.get_node(id));
            // And its edges
            surrounding.add_edges(graph.edges_of(graph.get_node(id)));
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
        
        if(consistent_alleles == 0) {
            // This read is inconsistent with all the alleles in the genotype,
            // so, given the genotype, the read must be sequenced or mapped
            // wrong.
            
            if(use_mapq) {
                // Compute P(mapped wrong or called wrong) = P(not (mapped right and called right)) = P(not (not mapped wrong and not called wrong))
                all_non_supporting_wrong *= (1.0 - (1.0 - phred_to_prob(read.mapping_quality())) * (1.0 - phred_to_prob(alignment_score(read))));
            } else {
                // Compute P(called wrong).
                all_non_supporting_wrong *= phred_to_prob(alignment_score(read));
            }
        } else {
            // This read is consistent with some of the alleles in the genotype,
            // so we must have drawn one of those alleles when sequencing.
            
            // So multiply in the probability that we hit one of those alleles
            all_supporting_drawn *= ((double) consistent_alleles) / genotype.size();
        }
        
    }
    
    // Now we've looked at all the reads, so AND everything together
    return all_non_supporting_wrong * all_supporting_drawn;
}

Genotype Genotyper::get_genotype(const vector<Path>& superbubble_paths, const map<Alignment*, vector<double>>& affinities) {
    
    // Freebayes way (improved with multi-support):
    
    // Decide that each read supports one or more alleles
    
    
    
    
    
    // Compute total affinity for each path
    vector<double> total_affinity(superbubble_paths.size(), 0);
    for(auto& alignment_and_affinities : affinities) {
        for(size_t i = 0; i < alignment_and_affinities.second.size(); i++) {
            // For each affinity of an alignment to a path
            double affinity = alignment_and_affinities.second.at(i);
            // Add it in to the total
            total_affinity.at(i) += affinity;
        }
    }
    
    Genotype genotype;
    
    if(superbubble_paths.empty()) {
        // No paths! Nothing to say!
        throw runtime_error("No paths through superbubble! Can't genotype!");
    }
    
    // If there's only one path, call hom ref
    if(superbubble_paths.size() < 2) {
        *genotype.add_allele() = superbubble_paths.front();
        Support* support = genotype.add_support();
        support->set_quality(total_affinity.front());
        return genotype;
    }
    
    // Find the two best paths
    int best_allele = -1;
    int second_best_allele = -1;
    
    for(int i = 0; i < superbubble_paths.size(); i++) {
        if(best_allele == -1 || total_affinity[best_allele] < total_affinity[i]) {
            // We have a new best allele
            second_best_allele = best_allele;
            best_allele = i;
        } else if(second_best_allele == -1 || total_affinity[second_best_allele] < total_affinity[i]) {
            // We have a new second best allele
            second_best_allele = i;
        }
    }
    
    // Add the best allele with the most support
    *genotype.add_allele() = superbubble_paths[best_allele];
    Support* best_support = genotype.add_support();
    best_support->set_quality(total_affinity[best_allele]);
    
    if(total_affinity[best_allele] > total_affinity[second_best_allele] * max_het_bias) {
        // If we're too biased to one side, call hom that
        return genotype;
    }
    
    // Else call het by adding the second best allele as well
    *genotype.add_allele() = superbubble_paths[second_best_allele];
    Support* second_best_support = genotype.add_support();
    second_best_support->set_quality(total_affinity[second_best_allele]);
    
    // TODO: use more support fields by investigating the Alignment associated
    // with the affinity, as it relates to the allele's path
    
    return genotype;
}

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


