#include <cstdint>
#include "genotyper.hpp"
#include "traversal_finder.hpp"
#include "cactus_snarl_finder.hpp"
#include "path_index.hpp"
#include "utility.hpp"
#include "translator.hpp"

//#define debug

namespace vg {

using namespace std;

// genotyper:
// use graph and reads to:
// - Augment graph
// - Find ultrabubbles or cactus branches to determine snarls
// - Generate proposals for paths through each snarl (from reads?)
// - Compute affinities of each read for each proposed path through a snarl
// - Compute diploid genotypes for each snarl
// - Output as vcf or as native format

void Genotyper::run(AugmentedGraph& augmented_graph,
                    ostream& out,
                    string ref_path_name,
                    string contig_name,
                    string sample_name,
                    string augmented_file_name,
                    bool output_vcf,
                    bool output_json,
                    int length_override,
                    int variant_offset) {

    // Unpack the graph
    VG& graph = augmented_graph.graph;
    // And all of the alignments as embedded in the graph
    const vector<const Alignment*>& alignments = augmented_graph.get_alignments();
    // And the translator that maps to/from the original unaugmented graph
    Translator& translator = augmented_graph.translator;

    normal_aligners.resize(get_thread_count());
    quality_aligners.resize(get_thread_count());

    if(output_vcf && show_progress) {
#pragma omp critical (cerr)
        cerr << "Calling against path " << ref_path_name << endl;
    }

    if(sample_name.empty()) {
        // Set a default sample name
        sample_name = "SAMPLE";
    }
    
    if(show_progress) {
#pragma omp critical (cerr)
        cerr << "Augmented graph; got " << translator.translations.size() << " translations" << endl;
    }

    // Make sure that we actually have an index for traversing along paths.
    graph.paths.rebuild_mapping_aux();

    if(!augmented_file_name.empty()) {
        ofstream augmented_stream(augmented_file_name);
        graph.serialize_to_ostream(augmented_stream);
        augmented_stream.close();
    }

    // store the reads that are embedded in the augmented graph, by their unique names
    map<string, const Alignment*> reads_by_name;
    for(const auto& alignment : alignments) {
        reads_by_name[alignment->name()] = alignment;
        // Each alignment is already as fully embedded in the graph as it should be.
    }
    
    if (show_progress) {
#pragma omp critical (cerr)
        cerr << "Converted " << alignments.size() << " alignments to embedded paths" << endl;
    }

    // Set up our genotypekit members.
    SnarlFinder* snarl_finder;

    // Find all the snarls in either the main graph or the subset, and put them
    // in this SnarlManager
    SnarlManager manager = CactusSnarlFinder(graph, ref_path_name).find_snarls();

    // TODO: I've removed this option because it has the chance to create incorrect Snarls from the subset
    // graph that are not actually snarls in the full graph. If we are going to retain this feature, we need
    // to use the subset graph throughout the genotyping process, but that's currently not possible because
    // the rest of the code uses an AugmentedGraph instead of a VG.
    
//    // We need to decide if we want to work on the full graph or just on the subgraph that has any support.
//    // We make the subset a local out here so it will stick around as long as we
//    // use the SnarlManager.
//    VG subset;
//    if (subset_graph) {
//        subset = make_subset_graph(graph, ref_path_name, reads_by_name); 
//        manager = subset.paths.has_path(ref_path_name) ?
//            CactusSnarlFinder(subset, ref_path_name).find_snarls()
//            : CactusSnarlFinder(subset).find_snarls();
//    } else {
//        algorithms::sort(&graph);
//        manager = CactusSnarlFinder(graph, ref_path_name).find_snarls();
//    }

    if(show_progress) {
        // Count up the ultrabubbles
        size_t ultrabubbles = 0;        
        manager.for_each_snarl_preorder([&](const Snarl* snarl) {
            if (snarl->type() == ULTRABUBBLE) {
                ultrabubbles++;
            }
        });    
#pragma omp critical (cerr)
            cerr << "Found " << ultrabubbles << " ultrabubbles" << endl;
    }

    // We're going to count up all the affinities we compute
    size_t total_affinities = 0;

    // We need a buffer for output
    vector<vector<Locus>> buffer;
    int thread_count = get_thread_count();
    buffer.resize(thread_count);

    // If we're doing VCF output we need a VCF header
    vcflib::VariantCallFile* vcf = nullptr;
    // And a reference index tracking the primary path
    PathIndex* reference_index = nullptr;
    if(output_vcf) {
        // Build a reference index on our reference path
        // Make sure to capture the sequence
        reference_index = new PathIndex(graph, ref_path_name, true);

        // Start up a VCF
        vcf = start_vcf(cout, *reference_index, sample_name, contig_name, length_override);
    }

    manager.for_each_snarl_parallel([&](const Snarl* snarl) {
        // For each snarl in parallel

        if (snarl->type() != ULTRABUBBLE) {
            // We only work on ultrabubbles right now
            cerr << "Skip snarl " << snarl->start() << " - " << snarl->end() << " due to not being an ultrabubble" << endl;
            return;
        }

        // Get the contents
        pair<unordered_set<id_t>, unordered_set<edge_t> > snarl_contents =
           manager.deep_contents(snarl, graph, true);

        // Test if the snarl can be longer than the reads
        bool read_bounded = is_snarl_smaller_than_reads(augmented_graph, snarl, snarl_contents, reads_by_name);
        TraversalAlg use_traversal_alg = traversal_alg;
        if (traversal_alg == TraversalAlg::Adaptive) {
            use_traversal_alg = read_bounded ? TraversalAlg::Reads : TraversalAlg::Representative;
        }

        if (use_traversal_alg == TraversalAlg::Exhaustive &&
            !manager.is_leaf(snarl)) {
            // The SupportRestrictedTraversalFinder we use in Exhaustive mode
            // can only handle leaf snarls.

            // Todo : support nesting hierarchy!
            if (show_progress) {
                cerr << "Skip snarl " << snarl->start() << " - " << snarl->end()
                    << " because it isn't a leaf and traversal algorithm "
                    << alg2name[use_traversal_alg] << " only works on leaves" << endl;
            }
            return;
        }

        if (use_traversal_alg == TraversalAlg::Representative && !manager.all_children_trivial(snarl, graph)) {
            // The RepresentativeTraversalFinder works for root and leaf
            // snarls, but unless we're in a leaf snarl, or a snarl with only
            // trivial children, it outputs traversals with child snarls in
            // them that the rest of genotype can't yet handle.

            // Todo : support nesting hierarchy!
            if (show_progress) {
                cerr << "Skip snarl " << snarl->start() << " - " << snarl->end()
                    << " because it has nontrivial children and traversal algorithm "
                    << alg2name[use_traversal_alg] << " will produce nested child snarl traversals" << endl;
            }
            return;
        }


        if (use_traversal_alg == TraversalAlg::Reads && !manager.is_root(snarl)) {
            // The ReadRestrictedTraversalFinder only works for root snarls.
            // TODO: How do we know this?
            
            // Todo : support nesting hierarchy!
            if (show_progress) {
                cerr << "Skip snarl " << snarl->start() << " - " << snarl->end()
                    << " because it isn't a root and traversal algorithm "
                    << alg2name[use_traversal_alg] << " only works on roots" << endl;
            }
            return;
        }
        
        // Report the snarl to our statistics code
        report_snarl(snarl, manager, reference_index, graph, reference_index);

        int tid = omp_get_thread_num();

        // Get the traverals
        vector<SnarlTraversal> paths = get_snarl_traversals(augmented_graph, manager, reads_by_name,
                                                            snarl, snarl_contents, reference_index,
                                                            use_traversal_alg);

        if(paths.empty()) {
            // Don't do anything for ultrabubbles with no routes through
            if(show_progress) {
#pragma omp critical (cerr)
                cerr << "Snarl " << snarl->start() << " - " << snarl->end() << " has " << paths.size() <<
                    " alleles: skipped for having no alleles" << endl;
            }
            return;
        }

        if(show_progress) {
#pragma omp critical (cerr)
            cerr << "Snarl " << snarl->start() << " - " << snarl->end() << " has " << paths.size() << " alleles" << endl;
            for(auto& path : paths) {
                // Announce each allele in turn
#pragma omp critical (cerr)
                cerr << "\t" << traversal_to_string(graph, path) << endl;
            }
        }

        // Compute the lengths of all the alleles
        set<size_t> allele_lengths;
        for(auto& path : paths) {
            allele_lengths.insert(traversal_to_string(graph, path).size());
        }

        // Get the affinities for all the paths
        map<const Alignment*, vector<Genotyper::Affinity>> affinities;

        if(allele_lengths.size() > 1 && (realign_indels || !read_bounded)) {
            // This is an indel, because we can change lengths. Use the slow route to do indel realignment.
            affinities = get_affinities(augmented_graph, reads_by_name, snarl, snarl_contents, manager, paths);
        } else {
            // Just use string comparison. Don't re-align when
            // length can't change, or when indle realignment is
            // off.
            affinities = get_affinities_fast(augmented_graph, reads_by_name, snarl, snarl_contents, manager, paths);
        }

        if(show_progress) {
            report_affinities(affinities, paths, graph);
            for(auto& alignment_and_affinities : affinities) {
#pragma omp critical (total_affinities)
                total_affinities += alignment_and_affinities.second.size();
            }
        }
        
        // Get a genotyped locus in the original frame
        Locus genotyped = genotype_snarl(graph, snarl, paths, affinities);

        if (output_vcf) {
            // Get 0 or more variants from the ultrabubble
            vector<vcflib::Variant> variants =
                locus_to_variant(graph, snarl, snarl_contents, manager, *reference_index,
                                 *vcf, genotyped, sample_name);
            for(auto& variant : variants) {
                // Fix up all the variants
                if(!contig_name.empty()) {
                    // Override path name
                    variant.sequenceName = contig_name;
                } else {
                    // Keep path name
                    variant.sequenceName = ref_path_name;
                }
                variant.position += variant_offset;

#pragma omp critical(cout)
                cout << variant << endl;
            }
        } else {
            // project into original graph (only need to do if we augmented with edit)
            if (!translator.translations.empty()) {
                genotyped = translator.translate(genotyped);
            }
            // record a consistent name based on the start and end position of the first allele
            stringstream name;
            if (genotyped.allele_size() && genotyped.allele(0).mapping_size()) {
                name << make_pos_t(genotyped.allele(0).mapping(0).position())
                     << "_"
                     << make_pos_t(genotyped
                                   .allele(0)
                                   .mapping(genotyped.allele(0).mapping_size()-1)
                                   .position());
            }
            genotyped.set_name(name.str());
            if (output_json) {
                // Dump in JSON
#pragma omp critical (cout)
                cout << pb2json(genotyped) << endl;
            } else {
                // Write out in Protobuf
                buffer[tid].push_back(genotyped);
                vg::io::write_buffered(cout, buffer[tid], 100);
            }
        }
    });

    if(!output_json && !output_vcf) {
        // Flush the protobuf output buffers
        for(int i = 0; i < buffer.size(); i++) {
            vg::io::write_buffered(cout, buffer[i], 0);
        }
    } 


    if(show_progress) {
#pragma omp critical (cerr)
        cerr << "Computed " << total_affinities << " affinities" << endl;
    }

    // Dump statistics before the snarls go away, so the pointers won't be dangling
    print_statistics(cerr);

    if(output_vcf) {
        delete vcf;
        delete reference_index;
    }

}


pair<pair<int64_t, int64_t>, bool> Genotyper::get_snarl_reference_bounds(const Snarl* snarl, const PathIndex& index,
    const HandleGraph* graph) {
    // Grab the start and end node IDs.
    auto first_id = snarl->start().node_id();
    auto last_id = snarl->end().node_id();

    if(!index.by_id.count(first_id) || !index.by_id.count(last_id)) {
        // Snarl isn't actually on the reference path so return a sentinel.
        return make_pair(make_pair(-1, -1), false);
    }

    // The position we have stored for this start node is the first
    // position along the reference at which it occurs. Our bubble
    // goes forward in the reference, so we must come out of the
    // opposnarl end of the node from the one we have stored.
    auto referenceIntervalStart = index.by_id.at(first_id).first + graph->get_length(graph->get_handle(snarl->start().node_id(),
        snarl->start().backward()));

    // The position we have stored for the end node is the first
    // position it occurs in the reference, and we know we go into
    // it in a reference-concordant direction, so we must have our
    // past-the-end position right there.
    auto referenceIntervalPastEnd = index.by_id.at(last_id).first;

    // Is this bubble articulated backwards relative to the reference?
    bool snarl_is_reverse = false;

    if(referenceIntervalStart > referenceIntervalPastEnd) {
        // Everything we know about the snarl is backwards relative to the reference. Flip it around frontways.
        snarl_is_reverse = true;
        std::swap(first_id, last_id);
        // Recalculate reference positions Use the end node, which we've now
        // made first_id, to get the length offset to the start of the actual
        // internal variable bit.
        referenceIntervalStart = index.by_id.at(first_id).first + graph->get_length(graph->get_handle(snarl->end().node_id(),
            snarl->end().backward()));
        referenceIntervalPastEnd = index.by_id.at(last_id).first;
    }

    return make_pair(make_pair(referenceIntervalStart, referenceIntervalPastEnd), snarl_is_reverse);
}

/**
 * Turn the given path into an allele. Drops the first
 * and last mappings and looks up the sequences for the nodes of the others.
 */
string allele_to_string(VG& graph, const Path& allele) {
    stringstream stream;

    for(size_t i = 1; i < allele.mapping_size() - 1; i++) {
        // Get the sequence for each node mapping
        const Node* node = graph.get_node(allele.mapping(i).position().node_id());
        stream << mapping_sequence(allele.mapping(i), *node);
    }

    return stream.str();
}

int Genotyper::alignment_qual_score(VG& graph, const Snarl* snarl, const Alignment& alignment) {
    if(alignment.quality().empty()) {
        // Special case: qualities not given. Assume something vaguely sane so
        // we can genotype without quality.
#ifdef debug
#pragma omp critical (cerr)
        cerr << "No base qualities. Assuming default quality of " << default_sequence_quality << endl;
#endif
        return default_sequence_quality;
    }

    // Go through all the qualities in the snarl
    // TODO: can we do this in place?
    string relevant_qualities = get_qualities_in_snarl(graph, snarl, alignment);

    if(relevant_qualities.empty()) {
        // No qualities available internal to the snarl for this read. Must be a
        // pure-deletion allele.
#ifdef debug
#pragma omp critical (cerr)
        cerr << "No internal qualities. Assuming default quality of " << default_sequence_quality << endl;
#endif
        // TODO: look at bases on either side of the deletion.
        return default_sequence_quality;
    }

    double total = 0;
    for(auto& quality : relevant_qualities) {
#ifdef debug
#pragma omp critical (cerr)
        cerr << "Quality: " << (int)quality << endl;
#endif
        total += quality;
    }
#ifdef debug
#pragma omp critical (cerr)
    cerr << "Count: " << relevant_qualities.size() << endl;
#endif
    // Make the total now actually be an average
    total /= relevant_qualities.size();
    return round(total);
}

bool Genotyper::is_snarl_smaller_than_reads(AugmentedGraph& augmented_graph,
                                            const Snarl* snarl,
                                            const pair<unordered_set<id_t>, unordered_set<edge_t> >& contents,
                                            map<string, const Alignment*>& reads_by_name) {
    size_t read_length = reads_by_name.empty() ? 50 : reads_by_name.begin()->second->sequence().length();
    size_t snarl_total_length = 0;
    for (auto snarl_id : contents.first) {
        
        if (snarl_id != snarl->start().node_id() &&
            snarl_id != snarl->end().node_id()) {
            snarl_total_length += augmented_graph.graph.get_length(augmented_graph.graph.get_handle(snarl_id));
        }
        if (snarl_total_length >= read_length) {
            return false;
        }
    }
    return true;
}

vector<SnarlTraversal> Genotyper::get_snarl_traversals(AugmentedGraph& augmented_graph, SnarlManager& manager,
                                                       map<string, const Alignment*>& reads_by_name,
                                                       const Snarl* snarl,
                                                       const pair<unordered_set<id_t>, unordered_set<edge_t> >& contents,
                                                       PathIndex* ref_path_index,
                                                       TraversalAlg use_traversal_alg) {
    vector<SnarlTraversal> paths;
    vector<SnarlTraversal> read_paths;

    // Get all the paths through the snarl supported by real named paths
    PathRestrictedTraversalFinder path_trav_finder(augmented_graph.graph, manager, reads_by_name);
    paths = path_trav_finder.find_traversals(*snarl);

    unique_ptr<TraversalFinder> read_trav_finder;
    if (use_traversal_alg == TraversalAlg::Reads) {
        // the original genotyper logic : paths for each read that spans site
        read_trav_finder = unique_ptr<TraversalFinder>(new ReadRestrictedTraversalFinder(augmented_graph, manager));
    } else if (use_traversal_alg == TraversalAlg::Exhaustive) {
        // exhaustive search, but only visit nodes and edges with given support
        read_trav_finder = unique_ptr<TraversalFinder>(new SupportRestrictedTraversalFinder(augmented_graph, manager, 2, 2));
    } else if (use_traversal_alg == TraversalAlg::Representative) {
        // representative search from vg call.  only current implementation that works for big sites
        // Now start looking for traversals of the sites.
        auto* finder = new RepresentativeTraversalFinder(
            augmented_graph.graph, manager, 1000, 1000,
            100, 1, 1, [&] (const Snarl& site) -> PathIndex* {
                return ref_path_index;
            },
            [&] (id_t node) -> Support {
                return augmented_graph.get_support(node);
            },
            [&] (edge_t edge) -> Support {
                return augmented_graph.get_support(edge);
            });

        // Since we can't sensibly handle any children, glom trivial children in.
        finder->eat_trivial_children = true;

        read_trav_finder = unique_ptr<TraversalFinder>(finder);
    } else {
        assert(false);
    }

    read_paths = read_trav_finder->find_traversals(*snarl);
    
    if (!read_paths.empty()) {
        // We want to log stats on reads that read all the
        // way through snarls. But since we may be called
        // multiple times we need to send the unique read
        // name too.            
        report_snarl_traversal(snarl, manager, augmented_graph.graph);
        }
    
    // Deduplicate into "paths"
    unordered_set<string> seen_sequences;
    for (auto& trav : paths) {
        // Make a set of the sequences already represented in paths        
        seen_sequences.insert(traversal_to_string(augmented_graph.graph, trav));
    }
    for (auto& trav : read_paths) {
        string traversal_seq = traversal_to_string(augmented_graph.graph, trav);

        if (traversal_seq.find_first_of("Nn") == string::npos &&
            !seen_sequences.count(traversal_seq)) {
            // If the sequence of a traversal from read_paths isn't shared with
            // something in paths, take it. We already know it's unique in
            // read_paths.  We also make sure it has no Ns
            paths.push_back(trav);
        }
    }

    return paths;
}

template<typename T> inline void set_intersection(const unordered_set<T>& set_1, const unordered_set<T>& set_2,
                                                  unordered_set<T>* out_intersection ) {
    bool set_1_smaller = set_1.size() < set_2.size();
    const unordered_set<T>& smaller_set = set_1_smaller ? set_1 : set_2;
    const unordered_set<T>& larger_set = set_1_smaller ? set_2 : set_1;

    *out_intersection = unordered_set<T>();
    unordered_set<T>& intersection = *out_intersection;
    for (T item : smaller_set) {
        if (larger_set.count(item)) {
            intersection.insert(item);
        }
    }
}


// TODO properly handle cycles inside ultrabubble by including multiplicity of an edge in a path
void Genotyper::edge_allele_labels(const VG& graph,
                                   const Snarl* snarl,
                                   const vector<list<NodeTraversal>>& snarl_paths,
                                   unordered_map<pair<NodeTraversal, NodeTraversal>,
                                   unordered_set<size_t>,
                                   hash_oriented_edge>* out_edge_allele_sets)
{
    // edges are indicated by the pair of node ids they connect and what orientation
    // they start and end in (true indicates forward)
    *out_edge_allele_sets = unordered_map<pair<NodeTraversal, NodeTraversal>, unordered_set<size_t>, hash_oriented_edge>();
    unordered_map<pair<NodeTraversal, NodeTraversal>, unordered_set<size_t>, hash_oriented_edge>& edge_allele_sets = *out_edge_allele_sets;

    for (size_t i = 0; i < snarl_paths.size(); i++) {
        list<NodeTraversal> path = snarl_paths[i];
        // start at second node so we can look at edge leading into it
        auto iter = path.begin();
        iter++;
        // can stop before last node because only interested in case where allele is ambiguous
        auto last_node = path.end();
        last_node--;
        for (; iter != last_node; iter++) {

            auto prev_iter = iter;
            prev_iter--;

            NodeTraversal node = *iter;
            NodeTraversal prev_node = *prev_iter;

            pair<NodeTraversal, NodeTraversal> edge = make_pair(prev_node, node);

            if (!edge_allele_sets.count(edge)) {
                edge_allele_sets.emplace(edge, unordered_set<size_t>());
            }

            // label the edge with this allele (indicated by its position in the vector)
            edge_allele_sets.at(edge).insert(i);
        }

    }
}

// find the log conditional probability of each ambiguous allele set given each true allele
void Genotyper::allele_ambiguity_log_probs(const VG& graph,
                                           const Snarl* snarl,
                                           const vector<list<NodeTraversal>>& snarl_paths,
                                           const unordered_map<pair<NodeTraversal, NodeTraversal>,
                                           unordered_set<size_t>,
                                           hash_oriented_edge>& edge_allele_sets,
                                           vector<unordered_map<vector<size_t>, double, hash_ambiguous_allele_set>>* out_allele_ambiguity_probs)
{
    *out_allele_ambiguity_probs = vector<unordered_map<vector<size_t>, double, hash_ambiguous_allele_set>>();
    vector<unordered_map<vector<size_t>, double, hash_ambiguous_allele_set>>& ambiguous_allele_probs = *out_allele_ambiguity_probs;
    ambiguous_allele_probs.reserve(snarl_paths.size());

    for (size_t i = 0; i < snarl_paths.size(); i++) {
        ambiguous_allele_probs[i] = unordered_map<vector<size_t>, double, hash_ambiguous_allele_set>();
        unordered_map<vector<size_t>, double, hash_ambiguous_allele_set>& allele_probs = ambiguous_allele_probs[i];
        list<NodeTraversal> path = snarl_paths[i];

        // consider both prefixes and suffixes that partially cross the ultrabubble
        for (bool forward : {true, false}) {
            // the set of alleles that this prefix of the current allele is consistent with
            unordered_set<size_t> prefix_allele_set;

            // find the first and last node to iterate through in this orientation
            // stop one node before last node because only interested in case where allele is ambiguous
            list<NodeTraversal>::iterator iter;
            list<NodeTraversal>::iterator final;
            if (forward) {
                // extend prefix forward through the ultrabubble
                iter = path.begin();
                final = path.end();
                final--;
            }
            else {
                // extend suffix backward through the ultrabubble
                iter = path.end();
                iter--;
                final = path.begin();
            }
            // iterate forwards or backwards along edges of path
            while (true) {
                // get the two nodes of the edge in the order they were entered into the allele label map
                NodeTraversal node;
                NodeTraversal next_node;
                auto next_iter = iter;
                if (forward) {
                    next_iter++;
                    node = *iter;
                    next_node = *next_iter;
                }
                else {
                    next_iter--;
                    node = *next_iter;
                    next_node = *iter;
                }
                if (next_iter == final) {
                    break;
                }

                pair<NodeTraversal, NodeTraversal> edge = make_pair(node, next_node);
                const unordered_set<size_t>& edge_allele_set = edge_allele_sets.at(edge);

                if (prefix_allele_set.empty()) {
                    // first edge in path, consistent with all alleles edge is labeled with
                    prefix_allele_set = edge_allele_set;
                }
                else {
                    // take set intersection of prefix alleles and the edge's allele labels
                    unordered_set<size_t> new_prefix_allele_set;
                    set_intersection<size_t>(prefix_allele_set, edge_allele_set, &new_prefix_allele_set);
                    prefix_allele_set = new_prefix_allele_set;
                }

                // convert unordered set into a sorted vector for consistent hash-key behavior
                vector<size_t> allele_set_key;
                allele_set_key.reserve(prefix_allele_set.size());
                for (size_t allele : prefix_allele_set) {
                    allele_set_key.emplace_back(allele);
                }
                sort(allele_set_key.begin(), allele_set_key.end());

                // add the length of the sequence to the probability of that allele set (will normalize later)
                if (allele_probs.count(allele_set_key)) {
                    allele_probs.at(allele_set_key) += node.node->sequence().length();
                }
                else {
                    allele_probs.at(allele_set_key) = node.node->sequence().length();
                }

                // iterate forward or backward through bubble
                if (forward) {
                    iter++;
                }
                else {
                    iter--;
                }
            }
        }

        // normalize lengths to probabilities (must sum to 1)
        size_t total_ambiguous_length = 0;
        for (auto& allele_length : allele_probs) {
            total_ambiguous_length += allele_length.second;
        }
        for (auto& allele_length : allele_probs) {
            allele_length.second = log(allele_length.second / total_ambiguous_length);
        }
    }
}




map<const Alignment*, vector<Genotyper::Affinity>>
    Genotyper::get_affinities(AugmentedGraph& aug,
                              const map<string, const Alignment*>& reads_by_name,
                              const Snarl* snarl,
                              const pair<unordered_set<id_t>, unordered_set<edge_t> >& contents,
                              const SnarlManager& manager,
                              const vector<SnarlTraversal>& snarl_paths) {

    // Grab our thread ID, which determines which aligner we get.
    int tid = omp_get_thread_num();
    
    Aligner& aligner = normal_aligners[tid];
    QualAdjAligner& qa_aligner = quality_aligners[tid];
    vector<MaximalExactMatch> mems; // to make sure we're using right align() interface.

    // We're going to build this up gradually, appending to all the vectors.
    map<const Alignment*, vector<Affinity>> to_return;

    // What reads are relevant to this ultrabubble?
    set<string> relevant_read_names;

    // What IDs are visited by these reads?
    unordered_set<id_t> relevant_ids;

#ifdef debug
#pragma omp critical (cerr)
    cerr << "Snarl contains " << contents.first.size() << " nodes" << endl;
#endif

    for(id_t node_id : contents.first) {
        // For every node in the ultrabubble, what reads visit it?
        for (const Alignment* aln : aug.get_alignments(node_id)) {
            // Each read that visits this node is relevant.
            relevant_read_names.insert(aln->name());

            for (size_t i = 0; i < aln->path().mapping_size(); ++i) {
                relevant_ids.insert(aln->path().mapping(i).position().node_id());
            }
        }
        // TODO: We populate relevant_read_names and then get the alignments by
        // name. Maybe we can just use alignment pointers?
    }

    for(id_t node_id : contents.first) {
        // Throw out all the IDs that are also used in the ultrabubble itself
        relevant_ids.erase(node_id);
    }

    // This is a temporary hack to break out of this function if there's too much to do. 
    // To fix properly will require a more general get_affinities_fast-type function, as well as,
    // probably, heuristics to reduce the search space
    if (relevant_read_names.size() * snarl_paths.size() > 1000) {
#pragma omp critical (cerr)
        cerr << "Skipping snarl " << pb2json(*snarl) << " with " << relevant_read_names.size() << " reads, "
             << snarl_paths.size() << " paths and " << contents.first.size()
             << " nodes as it is too complex to get affinities for (reads X paths > 1000)." << endl;
        return to_return;
    }

#ifdef debug
#pragma omp critical (cerr)
    cerr << relevant_read_names.size() << " reads visit an additional " << relevant_ids.size() << " nodes" << endl;
#endif

    // Make a vg graph with all the nodes used by the reads relevant to the
    // ultrabubble, but outside the ultrabubble itself.
    VG surrounding;
    for(auto id : relevant_ids) {
        // For all the IDs in the surrounding material
        if(min_recurrence != 0 &&
           aug.get_alignments(id).size() < min_recurrence) {
            // Skip nodes in the graph that have too little support. In practice
            // this means we'll restrict ourselves to supported, known nodes.
            // TODO: somehow do the same for edges.
            continue;
        }

        // Add each node and its edges to the new graph. Ignore dangling edges.
        // We'll keep edges dangling to the ultrabubble anchor nodes.
        surrounding.add_node(*aug.graph.get_node(id));
        surrounding.add_edges(aug.graph.edges_of(aug.graph.get_node(id)));
    }

    for(auto& path : snarl_paths) {
        // Now for each snarl path, make a copy of that graph with it in
        VG allele_graph(surrounding);

        for (size_t i = 0; i < path.visit_size(); i++) {
            // Add in every node on the path to the new allele graph
            allele_graph.add_node(*aug.graph.get_node(path.visit(i).node_id()));

            // Add in just the edge to the previous node on the path
            if(i != 0) {
                // There is something previous on the path.
                auto prev = i - 1;
                // Make an edge
                Edge path_edge;
                // And hook it to the correct side of the last node
                path_edge.set_from(path.visit(prev).node_id());
                path_edge.set_from_start(path.visit(prev).backward());
                // And the correct side of the next node
                path_edge.set_to(path.visit(i).node_id());
                path_edge.set_to_end(path.visit(i).backward());

                assert(aug.graph.has_edge(path_edge));

                // And add it in
                allele_graph.add_edge(path_edge);
            }
        }

        // Get rid of dangling edges
        allele_graph.remove_orphan_edges();

#ifdef debug_verbose
#pragma omp critical (cerr)
        cerr << "Align to " << pb2json(allele_graph.graph) << endl;
#endif

        // Grab the sequence of the path we are trying the reads against, so we
        // can check for identity across the snarl and not just globally for the
        // read.
        auto path_seq = traversal_to_string(aug.graph, path);

        for(auto& name : relevant_read_names) {
            // For every read that touched the ultrabubble, grab its original
            // Alignment pointer.
            const Alignment* read = reads_by_name.at(name);

            // Look to make sure it touches more than one node actually in the
            // ultrabubble, or a non-start, non-end node. If it just touches the
            // start or just touches the end, it can't be informative.
            set<id_t> touched_set;
            // Will this read be informative?
            bool informative = false;            
            for(size_t i = 0; i < read->path().mapping_size(); i++) {
                // Look at every node the read touches
                id_t touched = read->path().mapping(i).position().node_id();
                if(contents.first.count(touched)) {
                    // If it's in the ultrabubble, keep it
                    touched_set.insert(touched);
                }
            }

            if(touched_set.size() >= 2) {
                // We touch both the start and end, or an internal node.
                informative = true;
            } else {
                // Throw out the start and end nodes, if we touched them.
                touched_set.erase(snarl->start().node_id());
                touched_set.erase(snarl->end().node_id());
                if(!touched_set.empty()) {
                    // We touch an internal node
                    informative = true;
                }
            }

            if(!informative) {
                // We only touch one of the start and end nodes, and can say nothing about the ultrabubble. Try the next read.
                // TODO: mark these as ambiguous/consistent with everything (but strand?)
                continue;
            }

            // If we get here, we know this read is informative as to the internal status of this ultrabubble.
            Alignment aligned_fwd;
            Alignment aligned_rev;
            // We need a way to get graph node sizes to reverse these alignments
            auto get_node_size = [&](id_t id) {
                return aug.graph.get_node(id)->sequence().size();
            };

            if(read->sequence().size() == read->quality().size()) {
                // Re-align a copy to this graph (using quality-adjusted alignment).
                aligned_fwd = allele_graph.align_qual_adjusted(*read, &qa_aligner, mems);
                aligned_rev = allele_graph.align_qual_adjusted(reverse_complement_alignment(*read, get_node_size), &qa_aligner, mems);
            } else {
                // If we don't have the right number of quality scores, use un-adjusted alignment instead.
                aligned_fwd = allele_graph.align(*read, &aligner, mems);
                aligned_rev = allele_graph.align(reverse_complement_alignment(*read, get_node_size), &aligner, mems);
            }
            // Pick the best alignment, and emit in original orientation
            Alignment aligned = (aligned_rev.score() > aligned_fwd.score()) ? reverse_complement_alignment(aligned_rev, get_node_size) : aligned_fwd;

#ifdef debug
#pragma omp critical (cerr)
            cerr << path_seq << " vs " << aligned.sequence() << ": " << aligned.score() << endl;

#endif

#ifdef debug_verbose
#pragma omp critical (cerr)
            cerr << "\t" << pb2json(aligned) << endl;
#endif

            // Compute the score per base. TODO: is this at all comparable
            // between quality-adjusted and non-quality-adjusted reads?
            double score_per_base = (double)aligned.score() / aligned.sequence().size();

            // Save the score (normed per read base) and orientation
            // We'll normalize the affinities later to enforce the max of 1.0.
            Affinity affinity(score_per_base, aligned_rev.score() > aligned_fwd.score());

            // Compute the unnormalized likelihood of the read given the allele graph.
            if(read->sequence().size() == read->quality().size()) {
                // Use the quality-adjusted default scoring system
                affinity.likelihood_ln = qa_aligner.score_to_unnormalized_likelihood_ln(aligned.score());
            } else {
                // We will have aligned without quality adjustment, so interpret
                // score in terms of the normal scoring parameters.
                affinity.likelihood_ln = aligner.score_to_unnormalized_likelihood_ln(aligned.score());
            }

            // Get the NodeTraversals for the winning alignment through the snarl.
            auto read_traversal = get_traversal_of_snarl(aug.graph, snarl, manager, aligned.path());

            if(affinity.is_reverse) {
                // We really traversed this snarl backward. Flip it around.
                std::reverse(read_traversal.mutable_visit()->begin(), read_traversal.mutable_visit()->end());
                for (size_t i = 0; i < read_traversal.visit_size(); i++) {
                    // Flip around every traversal as well as reversing their order.
                    read_traversal.mutable_visit(i)->set_backward(!read_traversal.visit(i).backward());
                }

            }

            // Decide we're consistent if the alignment's string across the snarl
            // matches the string for the allele, anchored at the appropriate
            // ends.

            // Get the string this read spells out in its best alignment to this allele
            auto seq = traversal_to_string(aug.graph, read_traversal);

            // Now decide if the read's seq supports this path.
            if (read_traversal.visit_size() == 0) {
                // Aligned only to surrounding stuff, not consistent
                affinity.consistent = false;
            } else if(read_traversal.visit(0) == snarl->start() &&
               read_traversal.visit(read_traversal.visit_size() - 1) == snarl->end()) {
                // Anchored at both ends.
                // Need an exact match. Record if we have one or not.
                affinity.consistent = (seq == path_seq);
            } else if(read_traversal.visit(0) == snarl->start()) {
                // Anchored at start only.
                // seq needs to be a prefix of path_seq
                auto difference = std::mismatch(seq.begin(), seq.end(), path_seq.begin());
                // If the first difference is the past-the-end of the prefix, then it's a prefix
                affinity.consistent = (difference.first == seq.end());
            } else if(read_traversal.visit(read_traversal.visit_size() - 1) == snarl->end()) {
                // Anchored at end only.
                // seq needs to be a suffix of path_seq
                auto difference = std::mismatch(seq.rbegin(), seq.rend(), path_seq.rbegin());
                // If the first difference is the past-the-rend of the suffix, then it's a suffix
                affinity.consistent = (difference.first == seq.rend());
            } else {
                // This read doesn't touch either end. This might happen if the
                // snarl is very large. Just assume it's consistent and let
                // scoring work it out.
#pragma omp critical (cerr)
                cerr << "Warning: realigned read " << aligned.sequence() << " doesn't touch either end of its snarl!" << endl;
                affinity.consistent = true;
            }
            if(score_per_base < min_score_per_base) {
                // Say we can't really be consistent with this if we have such a
                // terrible score.
                affinity.consistent = false;
            }

            // Grab the identity and save it for this read and ultrabubble path
            to_return[read].push_back(affinity);

        }
    }

    for(auto& name : relevant_read_names) {
        // For every read that touched the ultrabubble, mark it consistent only
        // with its best-score alleles that don't mismatch in the allele.

        // So basically make everything that isn't normalized affinity 1.0
        // inconsistent if it wasn't already.

        // Grab its original Alignment pointer.
        const Alignment* read = reads_by_name.at(name);

        // Which is the best affinity we can get while being consistent?
        double best_consistent_affinity = 0;
        // And which is the best affinity we can get overall?
        double best_affinity = 0;
        for(auto& affinity : to_return[read]) {
            if(affinity.consistent && affinity.affinity > best_consistent_affinity) {
                // Look for the max affinity found on anything already consistent
                best_consistent_affinity = affinity.affinity;
            }

            if(affinity.affinity > best_affinity) {
                // And the max affinity overall, for normalizing to correct
                // affinity range of 0 to 1.
                best_affinity = affinity.affinity;
            }
        }

        for(auto& affinity : to_return[read]) {
            if(affinity.affinity < best_consistent_affinity) {
                // Mark all the Affinities that are worse than the best
                // consistent one as not actually being consistent.
                affinity.consistent = false;
            }

            if(best_affinity == 0) {
                affinity.affinity = 1.0;
            } else {
                // Normalize to the best affinity overall being 1.
                affinity.affinity /= best_affinity;
            }
        }

    }

    // After scoring all the reads against all the versions of the ultrabubble,
    // return the affinities
    return to_return;
}

map<const Alignment*, vector<Genotyper::Affinity> >
Genotyper::get_affinities_fast(AugmentedGraph& aug,
                               const map<string, const Alignment*>& reads_by_name,
                               const Snarl* snarl,
                               const pair<unordered_set<id_t>, unordered_set<edge_t> >& contents,
                               const SnarlManager& manager,
                               const vector<SnarlTraversal>& snarl_paths,
                               bool allow_internal_alignments) {

    // We're going to build this up gradually, appending to all the vectors.
    map<const Alignment*, vector<Affinity>> to_return;

    // What reads are relevant to this ultrabubble?
    set<string> relevant_read_names;

#ifdef debug
#pragma omp critical (cerr)
    cerr << "Ultrabubble contains " << contents.first.size() << " nodes" << endl;
#endif

    // Convert all the Paths used for alleles back to their strings.
    vector<string> allele_strings;
    for(auto& path : snarl_paths) {
        // Convert all the Paths used for alleles back to their strings.
        allele_strings.push_back(traversal_to_string(aug.graph, path));
    }

    for(id_t node_id : contents.first) {
        // For every node in the ultrabubble, what reads visit it?
        for (const Alignment* aln : aug.get_alignments(node_id)) {
            // Each read that visits this node is relevant.
            relevant_read_names.insert(aln->name());
        }
        // TODO: We populate relevant_read_names and then get the alignments by
        // name. Maybe we can just use alignment pointers?
    }

    for(auto name : relevant_read_names) {
        // For each relevant read, work out a string for the ultrabubble and whether
        // it's anchored on each end.

        // Make an Affinity to fill in
        Affinity base_affinity;

        // Get the NodeTraversals for this read through this snarl.
        auto read_traversal = get_traversal_of_snarl(aug.graph, snarl, manager, reads_by_name.at(name)->path());

        if(read_traversal.visit(0) == reverse(snarl->end()) ||
           read_traversal.visit(read_traversal.visit_size() - 1) == reverse(snarl->start())) {

            // We really traversed this snarl backward. Flip it around.
            
            std::reverse(read_traversal.mutable_visit()->begin(), read_traversal.mutable_visit()->end());
            for (size_t i = 0; i < read_traversal.visit_size(); i++) {
                // Flip around every traversal as well as reversing their order.
                read_traversal.mutable_visit(i)->set_backward(!read_traversal.visit(i).backward());
            }

            // We're on the reverse strand
            base_affinity.is_reverse = true;
        }

        if(read_traversal.visit_size() == 1 && (read_traversal.visit(0).node_id() == snarl->start().node_id() ||
            read_traversal.visit(read_traversal.visit_size() - 1).node_id() == snarl->end().node_id())) {
            
            // This read only touches the head or tail of the snarl, and so
            // cannot possibly be informative.
#ifdef debug
#pragma omp critical (cerr)
            cerr << "Non-informative read traversal being removed " << endl
                 << pb2json(read_traversal.visit(0)) << " to "
                 << pb2json(read_traversal.visit(read_traversal.visit_size() - 1)) << endl;  
#endif

            continue;
        }

        size_t total_supported = 0;

        // Get the string it spells out
        auto seq = traversal_to_string(aug.graph, read_traversal);

#ifdef debug
#pragma omp critical (cerr)
        cerr << "Consistency of " << reads_by_name.at(name)->sequence() << endl;
#endif

        // Now decide if the read's seq supports each path.
        for(auto& path_seq : allele_strings) {
            // We'll make an affinity for this allele
            Affinity affinity = base_affinity;
            if(read_traversal.visit(0) == snarl->start() &&
               read_traversal.visit(read_traversal.visit_size() - 1) == snarl->end()) {
                // Anchored at both ends.
                // Need an exact match. Record if we have one or not.
                affinity.consistent = (seq == path_seq);
            } else if(read_traversal.visit(0) == snarl->start()) {
                // Anchored at start only.
                // seq needs to be a prefix of path_seq
                auto difference = std::mismatch(seq.begin(), seq.end(), path_seq.begin());
                // If the first difference is the past-the-end of the prefix, then it's a prefix
                affinity.consistent = (difference.first == seq.end());
            } else if(read_traversal.visit(read_traversal.visit_size() - 1) == snarl->end()) {
                // Anchored at end only.
                // seq needs to be a suffix of path_seq
                auto difference = std::mismatch(seq.rbegin(), seq.rend(), path_seq.rbegin());
                // If the first difference is the past-the-rend of the suffix, then it's a suffix
                affinity.consistent = (difference.first == seq.rend());
            } else {
                // This read doesn't touch either end.
#pragma omp critical (cerr)
                cerr << "Warning: read " << read_traversal.visit(0) << " to " << read_traversal.visit(read_traversal.visit_size() - 1) << " doesn't touch either end of its snarl " << snarl->start() << " to " << snarl->end() << "!" << endl;
                if (allow_internal_alignments){

                }
            }

#ifdef debug
#pragma omp critical (cerr)
            cerr << "\t" << path_seq << " vs observed " << (read_traversal.visit(0).node_id() == snarl->start().node_id()) << " " << seq << " " << (read_traversal.visit(read_traversal.visit_size() - 1).node_id() == snarl->end().node_id()) << ": " << affinity.consistent << endl;
#endif

            // Fake a weight
            affinity.affinity = (double)affinity.consistent;
            to_return[reads_by_name.at(name)].push_back(affinity);

            // Add in to the total if it supports this
            total_supported += affinity.consistent;
        }

        if(total_supported == 0 && min_recurrence <= 1) {
            // This is weird. Doesn't match anything and we had no excuse to remove alleles.
#pragma omp critical (cerr)
            cerr << "Warning! Bubble sequence " << seq << " supports nothing!" << endl;
        }


    }


    // After scoring all the reads against all the versions of the ultrabubble,
    // return the affinities
    return to_return;
}

double Genotyper::get_genotype_log_likelihood(VG& graph, const Snarl* snarl, const vector<int>& genotype, const vector<pair<const Alignment*, vector<Affinity>>>& alignment_consistency) {
    // For each genotype, calculate P(observed reads | genotype) as P(all reads
    // that don't support an allele from the genotype are mismapped or
    // miscalled) * P(all reads that do support alleles from the genotype ended
    // up apportioned across the alleles as they are)

    // This works out to the product over all reads that don't support either
    // alleles of 1 - ((1 - MAPQ) * (1 - P(bases called wrong)), times the
    // likelihood of the observed (censored by multi-support) read counts coming
    // from the alleles they support, and the strands they are observed on.

    // TODO: handle contamination like Freebayes

    // This is the log probability that all reads that don't support either allele in this genotype are wrong.
    double all_non_supporting_wrong = prob_to_logprob(1);

    // This is the log probability that all the reads that do support alleles in this genotype were drawn from the genotype they support.
    double all_supporting_drawn = prob_to_logprob(1);

#ifdef debug
    if(genotype.size() == 2) {
#pragma omp critical (cerr)
        cerr << "Calculating P(a" << genotype[0] << "/a" << genotype[1] << ")" << endl;
    }
#endif

    // We want to keep count of how many reads support each allele in each
    // orientation. TODO: we might be sort of double-counting reads that support
    // multiple alleles, because we'll be treating them as having independent
    // orientations per allele, when really they're very likely to be oriented
    // the same way.

    // Maps from int allele number (as appears in Genotype) to total reads
    // forward and reverse giving support.
    map<int, pair<int, int>> strand_count_by_allele_and_orientation;

    for(auto& read_and_consistency : alignment_consistency) {
        // For each read, work out if it supports a genotype we have or not.

        // Split out the alignment from its consistency flags
        auto& read = *read_and_consistency.first;
        auto& consistency = read_and_consistency.second;

        // How many of the alleles in our genotype is it consistent with?
        int consistent_alleles = 0;
        // We only count each allele in the genotype once.
        set<int> alleles_seen;
        for(int allele : genotype) {
            if(alleles_seen.count(allele)) {
                // Counted up consistency with this allele already.
                continue;
            }
            alleles_seen.insert(allele);

            // For each unique allele in the genotype...

            if(consistency.size() > allele && consistency.at(allele).consistent) {
                // We're consistent with this allele
                consistent_alleles++;
            }
        }
        
        // Now we know how many alleles this read is consistent with.
        
        if (consistent_alleles == 1) {
            // If it supports exactly one allele, we add it to the per-strand support counts for that allele
            for(int allele : genotype) {
                if(consistency.size() > allele && consistency.at(allele).consistent) {
                    // We found the consustent allele again.
                    
                    if(consistency.at(allele).is_reverse) {
                        // Consistent with reverse
                        strand_count_by_allele_and_orientation[allele].second++;
                        break;
                    } else {
                        // Consistent with forward
                        strand_count_by_allele_and_orientation[allele].first++;
                        break;
                    }
                }
            }
        }

        auto read_qual = alignment_qual_score(graph, snarl, read);

#ifdef debug
#pragma omp critical (cerr)
        cerr << "Read (qual score " << read_qual << ") consistent with " << consistent_alleles
             << " genotype alleles observed." << endl;
#endif

        if(consistent_alleles == 0) {
            // This read is inconsistent with all the alleles in the genotype,
            // so, given the genotype, the read must be sequenced or mapped
            // wrong.

            double logprob_wrong;
            if(use_mapq) {
                // Compute P(mapped wrong or called wrong) = P(not (mapped right and called right)) = P(not (not mapped wrong and not called wrong))
                logprob_wrong = logprob_invert(logprob_invert(phred_to_logprob(read.mapping_quality())) +
                                               logprob_invert(phred_to_logprob(read_qual)));
            } else {
                // Compute P(called wrong).
                logprob_wrong = phred_to_logprob(read_qual);
            }

#ifdef debug
#pragma omp critical (cerr)
            cerr << "P(wrong) = " << logprob_to_prob(logprob_wrong) << endl;
#endif
            all_non_supporting_wrong += logprob_wrong;
        } else {
            // This read is consistent with some of the alleles in the genotype,
            // so we must have drawn one of those alleles when sequencing.

            // We account for this in a framework where we consider the reads
            // indistinguishable, so ignore it for now.
        }

    }

    // Multiply in in the probability that the supporting reads all came from
    // the strands they are on.
    double strands_as_specified = prob_to_logprob(1);
    // Each strand is equally likely
    vector<double> probs_by_orientation = {0.5, 0.5};
    for(auto& kv : strand_count_by_allele_and_orientation) {
        // For the forward and reverse strand counts for all the alleles
        auto& forward_count = kv.second.first;
        auto& reverse_count = kv.second.second;

        // Convert to a vector to satisfy the multinomial PMF function.
        vector<int> obs = {forward_count, reverse_count};

        // Get the log probability under multinomial.

        // TODO: we're counting oriented reads supporting multiple alleles
        // multiple times. Maybe we should look at orientation overall somehow?
        // Or treat orientations of alleles as the thing we do our master
        // censored multinomial likelihood over?
        double logprob = multinomial_sampling_prob_ln(probs_by_orientation, obs);

#ifdef debug
#pragma omp critical (cerr)
        cerr << "Allele "  << kv.first << " supported by " << forward_count << " forward, "
             << reverse_count << " reverse (P=" << logprob_to_prob(logprob) << ")" << endl;
#endif
        strands_as_specified += logprob;


    }

    // Multiply in probability that the reads came from alleles they support,
    // treating reads as indistinguishable and using a multinomial/binomial
    // model.
    double alleles_as_specified = prob_to_logprob(1);
    if(genotype.size() == 2 && genotype.at(0) != genotype.at(1)) {
        // For diploid heterozygotes, we can easily handle multi-support as
        // censorship. We know that at least the reads that only support allele
        // 0 are from allele 0, and that at most the reads that only support
        // allele 0 and those that support both alleles all are from allele 0.
        // We end up summing over a normalized choose (since the success
        // probability is known and fixed at 1/2), as described in Frey and
        // Marrero (2008) "A surprising MLE for Interval-Censored Binomial
        // Data". Diploid homozygotes don't need any of this logic, and keep the
        // default probability of 1 above for the support distribution across
        // alleles.

        // Work out how many reads support allele 0 only.
        int first_only_reads = 0;
        // And how many support both
        int ambiguous_reads = 0;
        // And how many total reads there are (# of trials).
        int total_reads = 0;

        for(auto& read_and_consistency : alignment_consistency) {
            auto& consistency = read_and_consistency.second;

            if(consistency.size() <= max(genotype.at(0), genotype.at(1))) {
                // No consistency info calculated for this useless uninformative
                // read.
                continue;
            }

            if(consistency.at(genotype.at(0)).consistent) {
                // Read is consistent with first allele
                if(consistency.at(genotype.at(1)).consistent) {
                    // And also second, so it's ambiguous
                    ambiguous_reads++;
#ifdef debug
#pragma omp critical (cerr)
                    cerr << "Ambiguous read: " << read_and_consistency.first->sequence() << endl;
                    for(int i = 0; i < consistency.size(); i++) {
#pragma omp critical (cerr)
                        cerr << "\t" << i << ": " << consistency[i].consistent << endl;
                    }
#endif
                } else {
                    // And only first allele
                    first_only_reads++;
                }
                total_reads++;
            } else if(consistency.at(genotype.at(1)).consistent) {
                // It's consistent with only the second allele.
                total_reads++;
            }
            // Don't count reads inconsistent with the genotype in this analysis.
        }

        // Now do the likelihood. We know each atom will be weighted by the same
        // factor (for assigning all the reads at 50% probability) so we can
        // pull it out of the sum.
        double log_atom_weight = prob_to_logprob(0.5) * total_reads;

        // We calculate the probability of each atom, weight, then sum.
        vector<double> weighted_atom_logprobs;

        for(int i = first_only_reads; i <= first_only_reads + ambiguous_reads; i++) {
            // For each possible actual number of reads from the first allele,
            // add in the probability of that atom.
            auto weighted_atom_logprob = choose_ln(total_reads, i) + log_atom_weight;
            weighted_atom_logprobs.push_back(weighted_atom_logprob);

#ifdef debug
#pragma omp critical (cerr)
            cerr << "P(" << i << " from first allele) = " << logprob_to_prob(weighted_atom_logprob) << endl;
#endif

        }

        // Atoms are already weighted above with the shared weight, and then multiply in (to
        // probability of 1) the probability of observing this range of possible
        // totals of reads from each of the two alleles.
        alleles_as_specified = logprob_sum(weighted_atom_logprobs);

#ifdef debug
#pragma omp critical (cerr)
        cerr << "P(" << first_only_reads << " to "  << (first_only_reads + ambiguous_reads) << " from first allele) = "
             << logprob_to_prob(logprob_sum(weighted_atom_logprobs)) 
             << " = " << logprob_to_prob(alleles_as_specified) << endl;
#endif

    } else if(genotype.size() > 2) {
        // This is tougher. We have to use a censored multinomial.
        
        // First we want to compress things down to distinct alleles in the
        // genotype, and weight them by the number of times they occur. This
        // saves us work because the multinomial will have fewer categories and
        // less ambiguity.
        double per_allele_prob = 1.0/genotype.size();
        
        // We'll put unique alleles in this set and store the total probability
        // of each under it.
        unordered_map<int, double> unique_alleles;
        
        for (auto& allele : genotype) {
            unique_alleles[allele] += per_allele_prob;
        }
        
        // Now convert to probs format (vector)
        vector<double> probs;
        for (auto& kv : unique_alleles) {
            // For each allele in whatever order the set gave them, put in the probability.
            probs.push_back(kv.second);
        }
        
        // Now we will assign reads to ambiguity classes and count them in here.
        unordered_map<vector<bool>, int> reads_by_class;
        
        for(auto& read_and_consistency : alignment_consistency) {
            // For each read, look at what it is consistent with
            auto& consistency = read_and_consistency.second;
            
            // Compute an ambiguity class for it
            vector<bool> ambiguity_class;

            for (auto& kv : unique_alleles) {
                // For each unique allele number in the genotype, in their assigned order
                const auto& allele_number = kv.first;
                
                // Add the consistency bit for this read against this allele to the class
                ambiguity_class.push_back(consistency.at(allele_number).consistent);
            }
            
            // Count the read as being in its class.
            reads_by_class[ambiguity_class]++;    
        }
        
        // Compute the censored multinomial probability of these potentially ambiguous reads given these class probabilities.
        alleles_as_specified = multinomial_censored_sampling_prob_ln(probs, reads_by_class);
        
#ifdef debug
#pragma omp critical (cerr)
        cerr << "P(reads drawn match specified ambiguity classes) = " << alleles_as_specified << endl;
#endif
        
    }
    // Haploid or 0-ploid genotypes will always have all the reads drawn from
    // their source alleles in the distribution observed, as only one is
    // possible.

    // Now we've looked at all the reads, so AND everything together
    double total_logprob = all_non_supporting_wrong + all_supporting_drawn + strands_as_specified + alleles_as_specified;

#ifdef debug
    if(genotype.size() == 2) {
#pragma omp critical (cerr)
        cerr << "logP(a" << genotype[0] << "/a" << genotype[1] << ") = " << all_non_supporting_wrong << " + "
             << all_supporting_drawn << " + " << strands_as_specified << " + " << alleles_as_specified << " = "
             << total_logprob << endl;
    }
#endif

    return total_logprob;


}

double Genotyper::get_genotype_log_prior(const vector<int>& genotype) {
    // Start with a prior probability of 100%
    double prior_logprob = prob_to_logprob(1);
    
    // The model we are workign under is:
    // We may be diploid. If so, we look at het and hom sites.
    // If not diploid, we may be haploid.
    // If not haploid, we may be deleted (0-ploid).
    // If not deleted, we will be polyploid (3+). We follow a geometric distribution on the extra copies.
    
    if (genotype.size() == 2) {
        // It's diploid
        prior_logprob += diploid_prior_logprob;
    
        // Priors are boring: certain amount for het, inverse of that for everyone else
        if(genotype[0] != genotype[1]) {
            // This is a het!
            prior_logprob += het_prior_logprob;
        } else {
            // This is a homozygote. Much more common.
            prior_logprob += logprob_invert(het_prior_logprob);
        }
    } else {
        // Not diploid
        prior_logprob += logprob_invert(diploid_prior_logprob);
        
        if (genotype.size() == 1) {
            // We're haploid
            prior_logprob += haploid_prior_logprob;
        } else {
            // Not haploid either
            prior_logprob += logprob_invert(haploid_prior_logprob);
            
            if (genotype.empty()) {
                // We're 0-ploid
                prior_logprob += deleted_prior_logprob;
            } else {
                // We're not 0-ploid either
                prior_logprob += logprob_invert(deleted_prior_logprob);
                
                // We must be polyploid
                
                // How much extra ploidy do we have
                auto extra_copies = genotype.size() - 2;
                
                // Charge for each additional copy except the last at the
                // failure price, and then the last at the success price.
                prior_logprob += geometric_sampling_prob_ln(polyploid_prior_success_logprob, extra_copies);
                
            }
        }
    }
    
    return prior_logprob;
}

string Genotyper::get_qualities_in_snarl(VG& graph, const Snarl* snarl, const Alignment& alignment) {
    // We'll fill this in.
    stringstream to_return;

    // Are we currently in the snarl?
    bool in_snarl = false;
    // What Visit do we need to see to leave?
    Visit expected;

    // Where are we in the quality string?
    size_t quality_pos = 0;

    for(size_t i = 0; i < alignment.path().mapping_size(); i++) {
        // For every mapping in the path in order
        auto& mapping = alignment.path().mapping(i);

        // What Visit is this?
        Visit traversal = to_visit(mapping.position().node_id(), mapping.position().is_reverse());

        if(!in_snarl) {
            // If we aren't in the snarl, we may be entering
            if(traversal == snarl->start()) {
                // We entered through the start
                in_snarl = true;
                // We'll leave at the end
                expected = snarl->end();
            } else if(traversal == reverse(snarl->end())) {
                // We entered through the end
                in_snarl = true;
                // We'll leave when we hit the start in reverse
                expected = reverse(snarl->start());
            }
        }

        for(size_t j = 0; j < mapping.edit_size(); j++) {
            // For every edit
            auto& edit = mapping.edit(j);

            if(in_snarl && mapping.position().node_id() != snarl->start().node_id()
               && mapping.position().node_id() != snarl->end().node_id()) {
                // We're in the snarl but not on the start or end nodes.
                // TODO: qualities for a deletion/insertion?
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

        if(in_snarl && traversal == expected) {
            // This was the node we were supposed to leave the snarl at.
            in_snarl = false;
        }
    }

    return to_return.str();

}

Locus Genotyper::genotype_snarl(VG& graph,
                               const Snarl* snarl,
                               const vector<SnarlTraversal>& snarl_paths,
                               const map<const Alignment*, vector<Affinity>>& affinities) {

    // Freebayes way (improved with multi-support)

    // We're going to populate this locus
    Locus to_return;

    for(auto& path : snarl_paths) {
        // Convert each allele to a Path and put it in the locus
        Path* allele_path = to_return.add_allele();
        for (size_t i = 0; i < path.visit_size(); i++) {
            // Convert each Visit to a Mapping and put it in the path.
            // TODO: we'll need a way to fill these in recursively somehow.
            *allele_path->add_mapping() = to_mapping(path.visit(i), graph);
        }
    }

#ifdef debug
#pragma omp critical (cerr)
    cerr << "Looking between " << snarl->start() << " and " << snarl->end() << endl;
#endif

    // We'll fill this in with the alignments for this snarl and their consistency-with-alleles flags.
    vector<pair<const Alignment*, vector<Affinity>>> alignment_consistency;

    // We fill this in with totals of reads supporting alleles
    vector<int> reads_consistent_with_allele(snarl_paths.size(), 0);
    // And this with the same thing split out by forward and reverse strand
    vector<pair<int, int>> strand_support_for_allele(snarl_paths.size(), make_pair(0, 0));

    // We'll store affinities by read name and allele here, for printing later.
    map<string, vector<Affinity>> debug_affinities;

    // We track overall forward and reverse strand reads, of reads that
    // support any allele.
    size_t overall_forward_reads = 0;
    size_t overall_reverse_reads = 0;

    for(auto& alignment_and_affinities : affinities) {
        // We need to clip down to just the important quality values        
        const Alignment& alignment = *alignment_and_affinities.first;

        // Hide all the affinities where we can pull them later
        debug_affinities[alignment.name()] = alignment_and_affinities.second;

        // Affinities already know whether they are consistent with an allele. Don't second-guess them.
        // Fine even with alignment; no read we embedded should ever have non-perfect identity.

        // We'll set these if the read supports anything in a forward or reverse
        // orientation.
        bool is_forward = false;
        bool is_reverse = false;

        // Of the alleles available, how many are consistent with this read?
        size_t consistent_alleles = 0;
        for(size_t i = 0; i < alignment_and_affinities.second.size(); i++) {
            consistent_alleles += alignment_and_affinities.second.at(i).consistent;
        }

        if (consistent_alleles == 1) {
            // This read is consistent with exactly one allele. Count it.
            
            for(size_t i = 0; i < alignment_and_affinities.second.size(); i++) {
                if(alignment_and_affinities.second.at(i).consistent) {
                    // We found the consistent allele again
                
                    // This read is consistent with this allele
                    reads_consistent_with_allele[i]++;
                    if(alignment_and_affinities.second.at(i).is_reverse) {
                        // It is on the reverse strand
                        strand_support_for_allele[i].second++;
                        is_reverse = true;
                    } else {
                        // It is on the forward strand
                        strand_support_for_allele[i].first++;
                        is_forward = true;
                    }
                }
            }
        }

        if(is_forward) {
            if(is_reverse) {
                // This is weird
#pragma omp critical (cerr)
                cerr << "Warning! Read supports alleles as both forward and reverse!" << endl;
                // Just call it forward
            }
            // This read supports an allele forward, so call it a forward read for the snarl
            overall_forward_reads++;
        } else if(is_reverse) {
            // This read supports an allele reverse, so call it a reverse read for the snarl
            overall_reverse_reads++;
        } else if(min_recurrence <= 1) {
            // Reads generally ought to support at least one allele, unless we
            // have weird softclips or they were for elided non-recurrent
            // alleles.
#pragma omp critical (cerr)
            cerr << "Warning! Read supports no alleles!" << endl;
        }

        // Save the alignment and its affinities, which we use to get GLs.
        alignment_consistency.push_back(alignment_and_affinities);

    }

#ifdef debug
    for(int i = 0; i < snarl_paths.size(); i++) {
        // Build a useful name for the allele
        stringstream allele_name;
        for (size_t j = 0; j < snarl_paths[i].visit_size(); j++) {
            allele_name << snarl_paths[i].visit(j).node_id() << ",";
        }
#pragma omp critical (cerr)
        {
            cerr << "a" << i << "(" << allele_name.str() << "): " << reads_consistent_with_allele[i] << "/" << affinities.size() << " reads consistent" << endl;
            for(auto& read_and_consistency : alignment_consistency) {
                if(read_and_consistency.second.size() > i && 
                   read_and_consistency.second[i].consistent &&
                   read_and_consistency.first->sequence().size() < 30) {
                    // Dump all the short consistent reads
                    cerr << "\t" << read_and_consistency.first->sequence() << " " << debug_affinities[read_and_consistency.first->name()][i].affinity << endl;
                }
            }
        }
    }
#endif

    // We'll go through all the genotypes, fill in their probabilities, put them
    // in here, and then sort them to find the best.
    vector<Genotype> genotypes_sorted;

    for(int allele1 = 0; allele1 < snarl_paths.size(); allele1++) {
        // For each first allele in the genotype
        for(int allele2 = 0; allele2 <= allele1; allele2++) {
            // For each second allele so we get all order-independent combinations

            // Make the combo
            vector<int> genotype_vector = {allele1, allele2};

            // Compute the log probability of the data given the genotype
            double log_likelihood = get_genotype_log_likelihood(graph, snarl, genotype_vector, alignment_consistency);

            // Compute the prior
            double log_prior = get_genotype_log_prior(genotype_vector);

            // Apply Bayes Rule
            double log_posterior_unnormalized = log_likelihood + log_prior;

#ifdef debug
#pragma omp critical (cerr)
            {
                cerr << "P(obs | a" << allele1 << "/a" << allele2 << ") = " << logprob_to_prob(log_likelihood) <<
                    " (" << log_likelihood << ")" << endl;
                cerr << "P(a" << allele1 << "/a" << allele2 << ") = " << logprob_to_prob(log_prior) <<
                    " (" << log_prior << ")" << endl;
                cerr << "P(a" << allele1 << "/a" << allele2 << " | obs) * P(obs) = " <<
                    logprob_to_prob(log_posterior_unnormalized) << " (" << log_posterior_unnormalized << ")" << endl;
            }
#endif

            // Fill in the actual Genotype object
            Genotype genotype;
            genotype.set_log_likelihood(log_likelihood);
            genotype.set_log_prior(log_prior);
            genotype.set_log_posterior(log_posterior_unnormalized);

            for(auto allele_id : genotype_vector) {
                // Copy over all the indexes of alleles in the genotype
                genotype.add_allele(allele_id);
            }

            // Put it in to sort
            genotypes_sorted.push_back(genotype);
        }
    }

    // Sort the genotypes in order of descending log posterior.
    sort(genotypes_sorted.begin(), genotypes_sorted.end(), [](const Genotype& a, const Genotype& b) {
            return a.log_posterior() > b.log_posterior();
        });

    for(size_t i = 0; i < snarl_paths.size(); i++) {
        // For each allele, make a support
        Support* support = to_return.add_support();
        // Set forward and reverse depth
        support->set_forward(strand_support_for_allele[i].first);
        support->set_reverse(strand_support_for_allele[i].second);
    }

    for(auto& genotype : genotypes_sorted) {
        // Add a genotype to the Locus for every one we looked at, in order by descending posterior
        *to_return.add_genotype() = genotype;
    }

    // Set up total support for overall depth
    Support* overall_support = to_return.mutable_overall_support();
    overall_support->set_forward(overall_forward_reads);
    overall_support->set_reverse(overall_reverse_reads);

    // Now we've populated the genotype so return it.
    return to_return;
}

void Genotyper::write_vcf_header(std::ostream& stream, const std::string& sample_name, const std::string& contig_name, size_t contig_size) {
    stream << "##fileformat=VCFv4.2" << std::endl;
    stream << "##ALT=<ID=NON_REF,Description=\"Represents any possible alternative allele at this location\">" << std::endl;
    stream << "##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Total Depth\">" << std::endl;
    stream << "##INFO=<ID=XSBB,Number=1,Type=Integer,Description=\"Ultrabubble Bases\">" << std::endl;
    stream << "##INFO=<ID=XSBN,Number=1,Type=Integer,Description=\"Ultrabubble Nodes\">" << std::endl;
    stream << "##FORMAT=<ID=DP,Number=1,Type=Float,Description=\"Read Depth\">" << std::endl;
    stream << "##FORMAT=<ID=GQ,Number=1,Type=Float,Description=\"Genotype Quality\">" << std::endl;
    stream << "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">" << std::endl;
    stream << "##FORMAT=<ID=AD,Number=.,Type=Float,Description=\"Allelic depths for the ref and alt alleles in the order listed\">" << std::endl;
    stream << "##FORMAT=<ID=SB,Number=4,Type=Float,Description=\"Forward and reverse support for ref and alt alleles.\">" << std::endl;
    stream << "##FORMAT=<ID=PL,Number=G,Type=String,Description=\"Log Likelihood\">" << std::endl;
    // We need this field to stratify on for VCF comparison. The info is in SB but vcfeval can't pull it out
    stream << "##FORMAT=<ID=XAAD,Number=1,Type=Float,Description=\"Alt allele read count.\">" << std::endl;
    if(!contig_name.empty()) {
        // Announce the contig as well.
        stream << "##contig=<ID=" << contig_name << ",length=" << contig_size << ">" << std::endl;
    }
    stream << "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" << sample_name << std::endl;
}

vcflib::VariantCallFile* Genotyper::start_vcf(std::ostream& stream, const PathIndex& index, const string& sample_name, const string& contig_name, size_t contig_size) {
    // Generate a vcf header. We can't make Variant records without a
    // VariantCallFile, because the variants need to know which of their
    // available info fields or whatever are defined in the file's header, so
    // they know what to output.
    // Handle length override if specified.
    std::stringstream headerStream;
    write_vcf_header(headerStream, sample_name, contig_name, contig_size > 0 ? contig_size : index.sequence.size());

    // Load the headers into a new VCF file object
    vcflib::VariantCallFile* vcf = new vcflib::VariantCallFile();
    std::string headerString = headerStream.str();
    assert(vcf->openForOutput(headerString));

    // Spit out the header
    stream << headerStream.str();

    // Give back the created VCF
    return vcf;
}

vector<vcflib::Variant>
Genotyper::locus_to_variant(VG& graph,
                            const Snarl* snarl,
                            const pair<unordered_set<id_t>, unordered_set<edge_t> >& contents,
                            const SnarlManager& manager,
                            const PathIndex& index,
                            vcflib::VariantCallFile& vcf,
                            const Locus& locus,
                            const string& sample_name) {

    // Make a vector to fill in
    vector<vcflib::Variant> to_return;

    // Make a new variant
    vcflib::Variant variant;
    // Attach it to the VCF
    variant.setVariantCallFile(vcf);
    // Fake the quality
    variant.quality = 0;

    // Make sure we have stuff
    if(locus.allele_size() == 0) {
        throw runtime_error("Can't turn an empty genotype into VCF");
    }
    if(locus.allele(0).mapping_size() == 0) {
        throw runtime_error("Can't turn an empty allele into VCF");
    }

    // Get the ultrabubble    
    auto first_id = snarl->start().node_id();
    auto last_id = snarl->end().node_id();

    if(!index.by_id.count(first_id) || !index.by_id.count(last_id)) {
        // We need to be anchored to the primary path to make a variant
#pragma omp critical (cerr)
        cerr << "Warning: Ultrabubble endpoints not on reference!" << endl;
        // If not return no variant
        return to_return;
    }

    // Compute the reference region occupied by the snarl, accounting for
    // orientation.
    auto bounds = get_snarl_reference_bounds(snarl, index, &graph);

    // Where does this bubble start and end in the reference?
    auto referenceIntervalStart = bounds.first.first;
    auto referenceIntervalPastEnd = bounds.first.second;

    // Is this bubble articulated backwards relative to the reference?
    bool snarl_is_reverse = bounds.second;
    if(snarl_is_reverse) {
        // Make sure our first and last IDs are actually accurate.
        std::swap(first_id, last_id);
    }

    // Get the string for the reference allele
    string ref_string = index.sequence.substr(referenceIntervalStart, referenceIntervalPastEnd - referenceIntervalStart);

    // And make strings for all the locus's alleles
    vector<string> allele_strings;

    for(size_t i = 0; i < locus.allele_size(); i++) {
        // Get the string for each allele
        string allele = allele_to_string(graph, locus.allele(i));
        if(snarl_is_reverse) {
            // Flip the alleles to match the reference orientation if necessary.
            allele = reverse_complement(allele);
        }
        allele_strings.push_back(allele);
    }

    // See if any alleles are empty
    bool empty_alleles = ref_string.empty();
    for(auto& allele : allele_strings) {
        if(allele == "") {
            empty_alleles = true;
        }
    }

    // Fix them up
    if(empty_alleles) {
        // Grab the character before our snarl
        string prefix = index.sequence.substr(referenceIntervalStart - 1, 1);
        for(auto& allele : allele_strings) {
            // Prepend it to every allele
            allele = prefix + allele;
        }
        // Also prepend it to the reference string
        ref_string = prefix + ref_string;

        // Budge the variant over
        referenceIntervalStart--;
    }

    // Trim fixed characters off the right
    while(ref_string.size() > 1) {
        // Grab the character at the end of the ref sequence
        char fixed = ref_string.back();

        // Set this to false if not all the alt strings also have this character
        // at the end, free to be chopped off.
        bool all_have_character = true;
        for(auto& allele_string : allele_strings) {
            if(allele_string.size() <= 1 || allele_string.back() != fixed) {
                // String is too short to trim or ends in the wrong character
                all_have_character = false;
                break;
            }
        }

        if(all_have_character) {
            // Trim it off
            ref_string.pop_back();
            for(auto& allele_string : allele_strings) {
                allele_string.pop_back();
            }

            // Record we budged the end of the interval left.
            referenceIntervalPastEnd--;
        } else {
            // Done trimming
            break;
        }
    }

    // Trim fixed characters off the left
    while(ref_string.size() > 1) {
        // Grab the character at the start of the ref sequence
        char fixed = ref_string.front();

        // Set this to false if not all the alt strings also have this character
        // at the start, free to be chopped off.
        bool all_have_character = true;
        for(auto& allele_string : allele_strings) {
            if(allele_string.size() <= 1 || allele_string.front() != fixed) {
                // String is too short to trim or starts with the wrong character
                all_have_character = false;
                break;
            }
        }

        if(all_have_character) {
            // Trim it off
            // TODO: this is going to be O(n^2)
            ref_string.erase(0, 1);
            for(auto& allele_string : allele_strings) {
                allele_string.erase(0, 1);
            }

            // Record that we budged the reference sequence start right.
            referenceIntervalStart++;
        } else {
            // Done trimming
            break;
        }
    }

    // Make the ref allele
    create_ref_allele(variant, ref_string);

    // Make a vector of supports by assigned VCF alt number
    vector<Support> support_by_alt;

    // This maps from locus allele index to VCF record alt number
    vector<int> allele_to_alt;

    // This records the max alt number used. We might have more alts than
    // alleles if the reference was never considered as an allele.
    int max_alt_number = 0;

    for(size_t i = 0; i < locus.allele_size(); i++) {
        // For each allele

        // Add it/find its number if it already exists (i.e. is the ref)
        int alt_number = add_alt_allele(variant, allele_strings[i]);

        // See if it's a new max
        max_alt_number = max(max_alt_number, alt_number);

        // Remember what VCF alt number it got
        allele_to_alt.push_back(alt_number);

        if(i < locus.support_size()) {
            if(support_by_alt.size() <= alt_number) {
                // Make sure we have a slot to put the support in
                support_by_alt.resize(alt_number + 1);
            }
            // Put it there
            support_by_alt[alt_number] = locus.support(i);
        }
    }

    // Get the best genotype
    assert(locus.genotype_size() > 0);
    const Genotype& best_genotype = locus.genotype(0);
    // And its support
    assert(locus.support_size() > 0);
    const Support& best_support = locus.support(0);
    // TODO: right now we only handle diploids
    assert(best_genotype.allele_size() == 2);

    // Compose the ML genotype
    variant.format.push_back("GT");
    auto& genotype_out = variant.samples[sample_name]["GT"];
    // Translate each allele to a VCF alt number, and put them in a string with the right separator
    genotype_out.push_back(to_string(allele_to_alt[best_genotype.allele(0)]) + (best_genotype.is_phased() ? "|" : "/")  +
                           to_string(allele_to_alt[best_genotype.allele(1)]));

    // Make sure that the called alleles have sufficient support on each strand
    for(size_t i = 0; i < best_genotype.allele_size(); i++) {
        // Check each allele marked present
        if(locus.support(best_genotype.allele(i)).forward() < min_unique_per_strand ||
           locus.support(best_genotype.allele(i)).reverse() < min_unique_per_strand) {
            // If there's not enough support for that allele in an orientation, skip the snarl. 

#ifdef debug_verbose
#pragma omp critical (cerr)
            cerr << "Warning: dropping locus from VCF due to insufficient per-strand unique support "
                 << locus.support(best_genotype.allele(i)).forward() << ", " 
                 << locus.support(best_genotype.allele(i)).reverse() << endl;
#endif
            return to_return;
        }
    }

    // Put the total depth overall (not double-counting)
    string depth_string = std::to_string(locus.overall_support().forward() + locus.overall_support().reverse());
    variant.format.push_back("DP");
    variant.samples[sample_name]["DP"].push_back(depth_string);
    variant.info["DP"].push_back(depth_string); // We only have one sample, so variant depth = sample depth

    // Also the snarl statistics
    // Ultrabubble bases
    size_t ultrabubble_bases = 0;
    for(id_t node_id : contents.first) {
        ultrabubble_bases += graph.get_length(graph.get_handle(node_id));
    }
    variant.info["XSBB"].push_back(to_string(ultrabubble_bases));
    // Ultrabubble nodes
    variant.info["XSBN"].push_back(to_string(contents.first.size()));

    variant.format.push_back("GQ");
    if(locus.genotype_size() > 1) {
        // Compute a quality from the difference between the best and second-
        // best genotype posteriors. Really this should be:

        // P(genotype is wrong) = sum(P(genotype) over other genotypes) / sum(P(genotype) over all genotypes)

        // When best genotype is much more probable than second best, which is
        // much more probable than all the rest, this approximation woks well.
        variant.samples[sample_name]["GQ"].push_back(to_string(
            logprob_to_phred(locus.genotype(1).log_posterior() - best_genotype.log_posterior())));
    } else {
        // This is very unlikely to be wrong. It can only be wrong if all the
        // reads in support of ref missed the haplotype on which an alt is.
        // TODO: this isn't exactly right; we should somehow account here for
        // reads we threw out for being on non-recurrent alleles...
        int total_reads = best_support.forward() + best_support.reverse();
        // Compute the likelihood that we missed everything, multiply it by a
        // prior of 5% for just being completely wrong, and treat that as the
        // posterior for the second best genotype.
        double all_missed_logprob = prob_to_logprob(0.5) * total_reads + prob_to_logprob(0.05);
        variant.samples[sample_name]["GQ"].push_back(to_string(
            logprob_to_phred(all_missed_logprob - best_genotype.log_posterior())));
    }

    // Compose the allele-specific depth
    variant.format.push_back("AD");
    for(auto& support : support_by_alt) {
        // Add the forward and reverse support together and use that for AD for the allele.
        variant.samples[sample_name]["AD"].push_back(to_string(support.forward() + support.reverse()));
    }

    // Do support for ref and alt alleles by strand
    variant.format.push_back("SB");
    for(auto& support : support_by_alt) {
        // Add the forward and reverse support in sequence, for ref and all the alts.
        // TODO: make this really only have the alt that's called.
        variant.samples[sample_name]["SB"].push_back(to_string(snarl_is_reverse ? support.reverse() : support.forward()));
        variant.samples[sample_name]["SB"].push_back(to_string(snarl_is_reverse ? support.forward() : support.reverse()));
    }


    // Work out genotype log likelihoods
    // Make a vector to shuffle them into VCF order. Fill it with inf for impossible genotypes.
    vector<double> log_likelihoods((max_alt_number * (max_alt_number + 1)) / 2 + max_alt_number + 1,
                                   numeric_limits<double>::infinity());
    for(size_t i = 0; i < locus.genotype_size(); i++) {
        // For every genotype, calculate its VCF-order index based on the VCF allele numbers

        // TODO: we can only do diploids
        assert(locus.genotype(i).allele_size() == 2);

        // We first need the low and high alt numbers
        size_t low_alt = allele_to_alt.at(locus.genotype(i).allele(0));
        size_t high_alt = allele_to_alt.at(locus.genotype(i).allele(1));
        if(low_alt > high_alt) {
            // Flip them to be the right way around
            std::swap(low_alt, high_alt);
        }

        // Compute the position as (sort of) specified in the VCF spec
        size_t index = (high_alt * (high_alt + 1)) / 2 + low_alt;
        // Store the log likelihood
        log_likelihoods.at(index) = locus.genotype(i).log_likelihood();
#ifdef debug
#pragma omp critical (cerr)
        cerr << high_alt << "/" << low_alt << ": " << index << " = " << pb2json(locus.genotype(i)) << endl;
#endif
    }

    variant.format.push_back("PL");
    for(auto& log_likelihood : log_likelihoods) {
        // Add all the likelihood strings, normalizing against the best
        // TODO: the best may not actually be the chosen genotype, because we genotype on posteriors.
        double pl_phred = 0.;
        if (!std::isinf((double)log_likelihood) && !(std::isinf((double)best_genotype.log_likelihood()))) {
            pl_phred = logprob_to_phred(log_likelihood - best_genotype.log_likelihood());
        }
        variant.samples[sample_name]["PL"].push_back(to_string(pl_phred));
    }

    // Set the variant position (now that we have budged it left if necessary)
    variant.position = referenceIntervalStart + 1;

    // Return the variant, since we managed to make it
    to_return.push_back(variant);
    return to_return;

}

void Genotyper::report_snarl(const Snarl* snarl, const SnarlManager& manager, const PathIndex* index,
                             VG& graph, PathIndex* reference_index) {
    // print some reference coordinate information
    if(show_progress &&
       reference_index != nullptr &&
       reference_index->by_id.count(snarl->start().node_id()) && 
       reference_index->by_id.count(snarl->end().node_id())) {
        // This snarl is on the reference (and we are indexing a reference because we are going to vcf)

        // Where do the start and end nodes fall in the reference?
        auto start_ref_appearance = reference_index->by_id.at(snarl->start().node_id());
        auto end_ref_appearance = reference_index->by_id.at(snarl->end().node_id());

        // Are the ends running with the reference (false) or against it (true)
        auto start_rel_orientation = (snarl->start().backward() != start_ref_appearance.second);
        auto end_rel_orientation = (snarl->end().backward() != end_ref_appearance.second);


        // Determine where the snarl starts and ends along the reference path
#pragma omp critical (cerr)
        cerr << "Snarl " << snarl->start() << " - " << snarl->end() << " runs reference " <<
            start_ref_appearance.first << " to " <<
            end_ref_appearance.first << endl;
                    
        if(!start_rel_orientation && !end_rel_orientation &&
           end_ref_appearance.first < start_ref_appearance.first) {
            // The snarl runs backward in the reference (but somewhat sensibly).
#pragma omp critical (cerr)
            cerr << "Warning! Snarl runs backwards!" << endl;
        }

    }
    
    // TODO: is there an easier way to detect trivial snarls?
    auto contents = manager.shallow_contents(snarl, graph, true);
    if(contents.first.size() == 2) {
        // Skip degenerate snarls
        return;
    }

    // Remember that we have the snarl
#pragma omp critical (all_snarls)
    all_snarls.insert(snarl);

    if(index != nullptr) {
        // We have an index of the reference.

        // Figure out its reference length and log it.
        auto bounds = get_snarl_reference_bounds(snarl, *index, &graph);

        if(bounds.first.first == -1) {
            // It's not really on that path
            return;
        }

        int64_t length = bounds.first.second - bounds.first.first;

#pragma omp critical (snarl_reference_length_histogram)
        snarl_reference_length_histogram[length]++;

    }

}

void Genotyper::report_snarl_traversal(const Snarl* snarl, const SnarlManager& manager, VG& graph) {
    // TODO: is there an easier way to detect trivial snarls?
    auto contents = manager.shallow_contents(snarl, graph, true);
    if(contents.first.size() == 2) {
        // Skip degenerate snarls
        return;
    }

    // Mark this read as traversing this snarl
#pragma omp critical (snarl_traversals)
    snarl_traversals.insert(snarl);
}

void Genotyper::report_affinities(map<const Alignment*, vector<Genotyper::Affinity>>& affinities,
                                  vector<SnarlTraversal>& paths, VG& graph) {
    // Sum up all the affinity counts by consistency flags
    map<string, size_t> consistency_combo_counts;

    // And average raw scores by alleles they are
    // consistent with, for things consistent with just
    // one allele.
    vector<double> score_totals(paths.size());
    vector<size_t> score_counts(paths.size());

    for(auto& alignment_and_affinities : affinities) {
        // For every alignment, make a string describing which alleles it is consistent with.
        string consistency;

        // How many alleles are we consstent with?
        size_t consistent_allele_count = 0;
        // And which one is it, if it's only one?
        int chosen = -1;

        for(size_t i = 0; i < alignment_and_affinities.second.size(); i++) {
            auto& affinity = alignment_and_affinities.second.at(i);
            if(affinity.consistent) {
                // Consistent alleles get marked with a 1
                consistency.push_back('1');

                // Say we're consistent with an allele
                chosen = i;
                consistent_allele_count++;

            } else {
                // Inconsistent ones get marked with a 0
                consistency.push_back('0');
            }
        }

        if(consistent_allele_count == 1) {
            // Add in the non-normalized score for the average
            score_totals.at(chosen) += alignment_and_affinities.second.at(chosen).score;
            score_counts.at(chosen)++;
        }

#ifdef debug
#pragma omp critical (cerr)
        cerr << consistency << ": " << alignment_and_affinities.first->sequence() << endl;
#endif


        // Increment the count for that pattern
        consistency_combo_counts[consistency]++;
    }

#pragma omp critical (cerr)
    {
        cerr << "Support patterns:" << endl;
        for(auto& combo_and_count : consistency_combo_counts) {
            // Spit out all the counts for all the combos
            cerr << "\t" << combo_and_count.first << ": " << combo_and_count.second << endl;
        }

        cerr << "Average scores for unique support:" << endl;
        for(size_t i = 0; i < score_totals.size(); i++) {
            // Spit out average scores of uniquely supporting reads for each allele that has them.
            if(score_counts.at(i) > 0) {
                cerr << "\t" << traversal_to_string(graph, paths.at(i)) << ": "
                     << score_totals.at(i) / score_counts.at(i) << endl;
            } else {
                cerr << "\t" << traversal_to_string(graph, paths.at(i)) << ": --" << endl;
            }
        }

    }
}

void Genotyper::print_statistics(ostream& out) {
    // Dump our stats to the given ostream.

    out << "Statistics:" << endl;
    out << "Number of Non-Degenerate Snarls: " << all_snarls.size() << endl;

    // How many snarls were actually traversed by reads?
    size_t snarls_traversed = 0;
    for(const Snarl* snarl : all_snarls) {
        // For every snarl
        if(snarl_traversals.count(snarl)) {
            // If it has a set of read names and the set is nonempty, it was traversed
            snarls_traversed++;
        }
    }
    out << "Snarls traversed by reads: " << snarls_traversed << endl;

    // How many snarls are on the reference? Only those that have defined lengths
    size_t snarls_on_reference = 0;
    for(auto& length_and_count : snarl_reference_length_histogram) {
        snarls_on_reference += length_and_count.second;
    }
    out << "Snarls on reference: " << snarls_on_reference << endl;

    // What's the length distribution?
    out << "Snarl length distribution: " << endl;
    for(auto& length_and_count : snarl_reference_length_histogram) {
        // Dump length and count as a TSV bit.
        out << length_and_count.first << "\t" << length_and_count.second << endl;
    }

}

VG Genotyper::make_subset_graph(VG& graph, const string& ref_path_name,
                                map<string, const Alignment*>& reads_by_name) {
    VG subset;
    
    // We'll collect the supported subset of the original graph
    set<Node*> supported_nodes;
    set<Edge*> supported_edges;

    for(auto& name_and_read : reads_by_name) {
        // Go through all the paths followed by reads
        auto& path = name_and_read.second->path();
        for(size_t i = 0; i < path.mapping_size(); i++) {
            // Look at all the nodes we visit along this read
            id_t node_id = path.mapping(i).position().node_id();
            // Make sure they are all supported
            supported_nodes.insert(graph.get_node(node_id));

            if(i > 0) {
                // We also need the edge from the last mapping to this one.
                // Make the two sides we connected moving from the last mapping to this one.
                NodeSide last(path.mapping(i - 1).position().node_id(), !path.mapping(i - 1).position().is_reverse());
                NodeSide here(node_id, path.mapping(i).position().is_reverse());

                Edge* edge = graph.get_edge(last, here);

                if(edge == nullptr) {
                    cerr << "Error! Edge " << last << " to " << here
                         << " from path " << name_and_read.first << " is missing!" << endl;
                    exit(1);
                }

                // We know the graph will have the edge
                supported_edges.insert(edge);
            }

        }
    }

    // We also want to support all nodes and edges used by the reference path.
    // TODO: once Cactus can root without hints, we can discard this
    if(graph.paths.has_path(ref_path_name)) {
        // We actually have a reference path, so get it for traversing.
        list<mapping_t>& ref_mappings = graph.paths.get_path(ref_path_name);
        // We need to remember the previous mapping for finding edges
        list<mapping_t>::iterator last_mapping = ref_mappings.end();
        for(list<mapping_t>::iterator mapping = ref_mappings.begin(); mapping != ref_mappings.end(); ++mapping) {
            // For each mapping along the reference path

            // What node is it on?
            id_t node_id = mapping->node_id();
            // Make sure it is supported
            supported_nodes.insert(graph.get_node(node_id));

            if(last_mapping != ref_mappings.end()) {
                // We're coming from another mapping and need to support the edge

                NodeSide last(last_mapping->node_id(), !last_mapping->is_reverse());
                NodeSide here(node_id, mapping->is_reverse());

                Edge* edge = graph.get_edge(last, here);

                if(edge == nullptr) {
                    cerr << "Error! Edge " << last << " to " << here
                         << " from path " << ref_path_name << " is missing!" << endl;
                    exit(1);
                }

                // We know the graph will have the edge
                supported_edges.insert(edge);

            }

            // Save the iterator so we can get the next edge
            last_mapping = mapping;

        }
    }

    // Make the subset graph of only supported nodes and edges (which will
    // internally contain copies of all of them).
    subset.add_nodes(supported_nodes);
    subset.add_edges(supported_edges);

    if(graph.paths.has_path(ref_path_name)) {
        // Copy over the reference path
        subset.paths.extend(graph.paths.path(ref_path_name));
    }


    if(show_progress) {
#pragma omp critical (cerr)
        cerr << "Looking at subset of " << subset.size() << " nodes" << endl;
    }

    return subset;
}

}


