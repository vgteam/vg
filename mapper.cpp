#include "mapper.hpp"

namespace vg {

Mapper::Mapper(Index* idex, gcsa::GCSA* g, xg::XG* xidex)
    : index(idex)
    , gcsa(g)
    , xindex(xidex)
    , best_clusters(0)
    , cluster_min(2)
    , hit_max(100)
    , hit_size_threshold(512)
    , kmer_min(0)
    , kmer_threshold(1)
    , kmer_sensitivity_step(3)
    , thread_extension(1)
    , thread_extension_max(80)
    , max_thread_gap(30)
    , context_depth(1)
    , max_multimaps(1)
    , max_attempts(7)
    , softclip_threshold(0)
    , prefer_forward(false)
    , greedy_accept(false)
    , target_score_per_bp(1.5)
    , min_score_per_bp(0)
    , min_kmer_entropy(0)
    , debug(false)
    , alignment_threads(1)
{
    // Nothing to do. We just hold the default parameter values.
}

Mapper::Mapper(Index* idex, gcsa::GCSA* g) : Mapper(idex, g, nullptr)
{
    if(idex == nullptr) {
        // With this constructor we need an index.
        cerr << "error:[vg::Mapper] cannot create a RocksDB-based Mapper with null index" << endl;
        exit(1);
    }
    
    kmer_sizes = index->stored_kmer_sizes();
    if (kmer_sizes.empty() && gcsa == NULL) {
        cerr << "error:[vg::Mapper] the index (" 
             << index->name << ") does not include kmers"
             << " and no GCSA index has been provided" << endl;
        exit(1);
    }
}

Mapper::Mapper(xg::XG* xidex, gcsa::GCSA* g) : Mapper(nullptr, g, xidex) {
    if(xidex == nullptr) {
        // With this constructor we need an XG graph.
        cerr << "error:[vg::Mapper] cannot create an xg-based Mapper with null xg index" << endl;
        exit(1);
    }
    
    if(g == nullptr) {
        // With this constructor we need a GCSA2 index too.
        cerr << "error:[vg::Mapper] cannot create an xg-based Mapper with null GCSA2 index" << endl;
        exit(1);
    }
}

Mapper::Mapper(void) : Mapper(nullptr, nullptr, nullptr) {
    // Nothing to do. Default constructed and can't really do anything.
}

Mapper::~Mapper(void) {
    // noop
}

Alignment Mapper::align(const string& seq, int kmer_size, int stride, int band_width) {
    Alignment aln;
    aln.set_sequence(seq);
    return align(aln, kmer_size, stride, band_width);
}

// align read2 near read1's mapping location
void Mapper::align_mate_in_window(Alignment& read1, Alignment& read2, int pair_window) {
    if (read1.score() == 0) return; // bail out if we haven't aligned the first
    // try to recover in region
    Path* path = read1.mutable_path();
    int64_t idf = path->mutable_mapping(0)->position().node_id();
    int64_t idl = path->mutable_mapping(path->mapping_size()-1)->position().node_id();
    // but which way should we expand? this will make things much easier
    // just use the whole "window" for now
    int64_t first = max((int64_t)0, idf - pair_window);
    int64_t last = idl + (int64_t) pair_window;
    VG* graph = new VG;
    
    // Now we need to get the neighborhood by ID and expand outward by actual
    // edges. How we do this depends on what indexing structures we have.
    if(xindex) {
        // should have callback here
        xindex->get_id_range(first, last, graph->graph);
        xindex->expand_context(graph->graph, context_depth);
        graph->rebuild_indexes();
    } else if(index) {
        index->get_range(first, last, *graph);
        index->expand_context(*graph, context_depth);
    } else {
        cerr << "error:[vg::Mapper] cannot align mate with no graph data" << endl;
        exit(1);
    }
    
    
    graph->remove_orphan_edges();
    read2.clear_path();
    read2.set_score(0);
    
    graph->align(read2);
    delete graph;
}

pair<vector<Alignment>, vector<Alignment>> Mapper::align_paired_multi(
    Alignment& read1, Alignment& read2, int kmer_size, int stride, int band_width, int pair_window) {

    // use paired-end resolution techniques
    //
    // map both reads independently
    // if both map, return the pair of all the mappings
    // if one maps but not the other, attempt to rescue by mapping the other nearby in both orientations
    //     for each mapping of the mapped read, pick the best mapping of the unmapped read near it
    
    // problem: need to develop model of pair orientations
    // solution: collect a buffer of alignments and then align them using unpaired approach
    //           detect read orientation and mean (and sd) of pair distance

    vector<Alignment> alignments1 = align_multi(read1, kmer_size, stride, band_width);
    vector<Alignment> alignments2 = align_multi(read2, kmer_size, stride, band_width);
    
    // We have some logic around align_mate_in_window to handle orientation
    // Since we now support reversing edges, we have to at least try opposing orientations for the reads.
    auto align_mate = [&](Alignment& read, Alignment& mate) {
        // Make an alignment to align in the same local orientation as the read
        Alignment aln_same = mate;
        // And one to align in the opposite local orientation
        Alignment aln_opposite = mate;
        
        // Is reverse going to be preferable in a tie?
        if (read.is_reverse()) {
            // If so, reverse the "same direction" Alignment sequence
            aln_same.set_sequence(reverse_complement(aln_same.sequence()));
            aln_same.set_is_reverse(true);
        } else {
            // Otherwise reverse the opposite direction sequence
            aln_opposite.set_sequence(reverse_complement(aln_opposite.sequence()));
            aln_opposite.set_is_reverse(true);
        }
        
        // Do both the alignments
        align_mate_in_window(read, aln_same, pair_window);
        align_mate_in_window(read, aln_opposite, pair_window);
        
        if(aln_same.score() >= aln_opposite.score()) {
            // The alignment in the same direction is best
            mate = aln_same;
        } else {
            // We have to say this pair has opposing local orientations
            mate = aln_opposite;
        }
    };
    
    // Try to rescue the unmapped end with each of the mapped end's  alignments.
    if(alignments1.empty() && !alignments2.empty()) {
        // We need to try aligning the mate near each of the places we aligned
        // the read. But we need to deduplicate those alignments. So we
        // serialize them into this set.
        set<string> serializedAlignmentsUsed;
    
        for(auto& aln2 : alignments2) {
            // Align the mate the best way for each alignment of read 2
            Alignment mate = read1;
            align_mate(aln2, mate);
            
            // Work out what it is as a string
            string serialized;
            mate.SerializeToString(&serialized);
            
            if(!serializedAlignmentsUsed.count(serialized)) {
                // It's not a duplicate
                alignments1.push_back(mate);
                serializedAlignmentsUsed.insert(serialized);
            }
        }
        
        // Sort these alignments by score, descending
        sort(alignments1.begin(), alignments1.end(), [](const Alignment& a, const Alignment& b) {
            return a.score() > b.score();
        });
        
        // Set the secondary bits on all but the first rescued alignment
        for(size_t i = 1; i < alignments1.size(); i++) {
            alignments1[i].set_is_secondary(true);
        }
    } else if(alignments2.empty() && !alignments1.empty()) {
        // We need to try aligning the mate near each of the places we aligned
        // the read. But we need to deduplicate those alignments. So we
        // serialize them into this set.
        set<string> serializedAlignmentsUsed;
    
        for(auto& aln1 : alignments1) {
            // Align the mate the best way for each alignment of read 1
            Alignment mate = read2;
            align_mate(aln1, mate);
            
            // Work out what it is as a string
            string serialized;
            mate.SerializeToString(&serialized);
            
            if(!serializedAlignmentsUsed.count(serialized)) {
                // It's not a duplicate
                alignments2.push_back(mate);
                serializedAlignmentsUsed.insert(serialized);
            }
        }
        
        // Sort these alignments by score, descending
        sort(alignments2.begin(), alignments2.end(), [](const Alignment& a, const Alignment& b) {
            return a.score() > b.score();
        });
        
        // Set the secondary bits on all but the first rescued alignment
        for(size_t i = 1; i < alignments2.size(); i++) {
            alignments2[i].set_is_secondary(true);
        }
    } 
    
    // link the fragments
    for(size_t i = 0; i < alignments1.size(); i++) {
        alignments1[i].mutable_fragment_next()->set_name(read2.name());
    }
    for(size_t i = 0; i < alignments2.size(); i++) {
        alignments2[i].mutable_fragment_prev()->set_name(read1.name());
    }
    
    // TODO: sort in order of best overall score, if the alignments happen to actually correspond?
    // TODO: pathfind between alignments?
    
    // TODO
    // mark them as discordant if there is an issue?
    // this needs to be detected with care using statistics built up from a bunch of reads
    return make_pair(alignments1, alignments2);

}

pair<Alignment, Alignment> Mapper::align_paired(Alignment& read1, Alignment& read2, int kmer_size, int stride, 
    int band_width, int pair_window) {
 
    pair<vector<Alignment>, vector<Alignment>> multimappings = align_paired_multi(read1, read2, 
        kmer_size, stride, band_width, pair_window);
        
    // Grab the input reads as alignments if nothing was found, and the found alignments otherwise
    Alignment aln1;
    if(multimappings.first.empty()) {
        aln1 = read1;
        aln1.clear_path();
        aln1.set_score(0);
        // Make sure to link up alignments even if they aren't mapped.
        aln1.mutable_fragment_next()->set_name(read2.name());
    } else {
        aln1 = multimappings.first[0];
    }
    Alignment aln2;
    if(multimappings.second.empty()) {
        aln2 = read2;
        aln2.clear_path();
        aln2.set_score(0);
        aln2.mutable_fragment_prev()->set_name(read1.name());
    } else {
        aln2 = multimappings.second[0];
    }
    
    // Stick the alignments together
    return make_pair(aln1, aln2);
}

Alignment Mapper::align_banded(Alignment& read, int kmer_size, int stride, int band_width) {
    // split the alignment up into overlapping chunks of band_width size
    list<Alignment> alignments;
    // force used bandwidth to be divisible by 4
    // round up so we have > band_width
    //cerr << "trying band width " << band_width << endl;
    if (band_width % 4) {
        band_width -= band_width % 4; band_width += 4;
    }
    assert(read.sequence().size() > band_width);
    int div = 2;
    while (read.sequence().size()/div > band_width) {
        ++div;
    }
    int segment_size = read.sequence().size()/div;
    // use segment sizes divisible by 4, as this simplifies math
    // round up as well
    // we'll divide the overlap by 2 and 2 and again when stripping from the start
    // and end of sub-reads
    if (segment_size % 4) {
        segment_size -= segment_size % 4; segment_size += 4;
    }
    //cerr << "Segment size be " << segment_size << endl;
    // and overlap them too
    size_t to_align = div * 2 - 1; // number of alignments we'll do
    vector<Alignment> alns; alns.resize(to_align);
    vector<size_t> overlaps; overlaps.resize(to_align);
    
    // We need a function to get the lengths of nodes, in case we need to
    // reverse an Alignment, including all its Mappings and Positions. TODO:
    // make this cache the node lengths for the nodes used in the actual
    // alignments somehow?
    std::function<int64_t(int64_t)> get_node_length = [&](int64_t node_id) {
        if(xindex) {
            // Grab the node sequence only from the XG index and get its size.
            return xindex->node_sequence(node_id).size();
        } else if(index) {
            // Get a 1-element range from the index and then use that.
            VG one_node_graph;
            index->get_range(node_id, node_id, one_node_graph);
            return one_node_graph.get_node(node_id)->sequence().size();
        } else {
            // Complain we don;t have the right indices.
            // This should be caught before here.
            throw runtime_error("No index to get nodes from.");
        }
    };
    
#pragma omp parallel for
    for (int i = 0; i < div; ++i) {
        {
            Alignment aln = read;
            if (i+1 == div) {
                // ensure we get all the sequence
                aln.set_sequence(read.sequence().substr(i*segment_size));
            } else {
                aln.set_sequence(read.sequence().substr(i*segment_size, segment_size));
            }
            // todo, possible problem
            // overlap can possibly go to 100% of the "last" read
            // if we aren't careful about this it might cause a problem in the merge
            size_t overlap = (i == 0? 0 : segment_size/2);
            //cerr << "overlap is " << overlap << endl;
            size_t idx = 2*i;
            overlaps[idx] = overlap;
            Alignment mapped_aln = align(aln, kmer_size, stride);
            if ((float) mapped_aln.score() / (float) mapped_aln.sequence().size()
                >= min_score_per_bp) {
            
                // We're actually going to use this alignment. We should make
                // sure to flip it if it's backward, though. We need all the
                // alignments to be a consistent orientation for merging.
                alns[idx] = mapped_aln.is_reverse() ? reverse_alignment(mapped_aln, get_node_length) : mapped_aln;
            } else {
                alns[idx] = aln; // unmapped
            }
            if (debug) {
#pragma omp critical (cerr)
                cerr << pb2json(alns[idx]) << endl;
            }
        }
        // step to next position
        // and the overlapped bit --- here we're using a hard-coded 50% overlap
        if (i != div-1) { // if we're not at the last sequence
            Alignment aln = read;
            aln.set_sequence(read.sequence().substr(i*segment_size+segment_size/2,
                                                    segment_size));
            size_t overlap = segment_size/2;
            size_t idx = 2*i+1;
            overlaps[idx] = overlap;
            Alignment mapped_aln = align(aln, kmer_size, stride);
            if ((float) mapped_aln.score() / (float) mapped_aln.sequence().size()
                >= min_score_per_bp) {
                alns[idx] = mapped_aln.is_reverse() ? reverse_alignment(mapped_aln, get_node_length) : mapped_aln;
            } else {
                alns[idx] = aln; // unmapped
            }
            if (debug) {
#pragma omp critical (cerr)
                cerr << pb2json(alns[idx]) << endl;
            }
        }
    }
    // by telling our merge the expected overlaps, it will correctly combine the alignments
    return merge_alignments(alns, overlaps, debug);
}

vector<Alignment> Mapper::align_multi(Alignment& aln, int kmer_size, int stride, int band_width) {

    if (aln.sequence().size() > band_width) {
        if (debug) cerr << "switching to banded alignment" << endl;
        return vector<Alignment>{align_banded(aln, kmer_size, stride, band_width)};
    }

    std::chrono::time_point<std::chrono::system_clock> start_both, end_both;
    if (debug) start_both = std::chrono::system_clock::now();
    const string& sequence = aln.sequence();

    // if kmer size is not specified, pick it up from the index
    // for simplicity, use the first available kmer size; this could change
    if (kmer_size == 0) kmer_size = *kmer_sizes.begin();
    // and start with stride such that we barely cover the read with kmers
    if (stride == 0)
        stride = sequence.size()
            / ceil((double)sequence.size() / kmer_size);

    int kmer_hit_count = 0;
    int kept_kmer_count = 0;

    if (debug) cerr << "aligning " << aln.sequence() << endl;

    // This will hold the best forward alignment (or an alignment with no apth and 0 score if no alignment is found).
    Alignment best_f = aln;
    // This will hold all of the forward alignments up to max_multimaps
    vector<Alignment> alignments_f;

    // This will similarly hold all the reverse alignments.
    // Right now we set it up to provide input to the actual alignment algorithm.
    Alignment best_r = aln;
    best_r.set_sequence(reverse_complement(aln.sequence()));
    best_r.set_is_reverse(true);
    // This will hold all of the reverse alignments up to max_multimaps
    vector<Alignment> alignments_r;

    auto increase_sensitivity = [this,
                                 &kmer_size,
                                 &stride,
                                 &sequence,
                                 &best_f,
                                 &best_r]() {
        kmer_size -= kmer_sensitivity_step;
        stride = sequence.size() / ceil((double)sequence.size() / kmer_size);
        if (debug) cerr << "realigning with " << kmer_size << " " << stride << endl;
        /*
        if ((double)stride/kmer_size < 0.5 && kmer_size -5 >= kmer_min) {
            kmer_size -= 5;
            stride = sequence.size() / ceil((double)sequence.size() / kmer_size);
            if (debug) cerr << "realigning with " << kmer_size << " " << stride << endl;
        } else if ((double)stride/kmer_size >= 0.5 && kmer_size >= kmer_min) {
            stride = max(1, stride/3);
            if (debug) cerr << "realigning with " << kmer_size << " " << stride << endl;
        }
        */
    };

    int attempt = 0;
    int kmer_count_f = 0;
    int kmer_count_r = 0;

    while (best_f.score() == 0 && best_r.score() == 0 && attempt < max_attempts) {

        {
            std::chrono::time_point<std::chrono::system_clock> start, end;
            if (debug) start = std::chrono::system_clock::now();
            // Go get all the forward alignments, putting the best one in best_f.
            alignments_f = align_threaded(best_f, kmer_count_f, kmer_size, stride, attempt);
            if (debug) {
                end = std::chrono::system_clock::now();
                std::chrono::duration<double> elapsed_seconds = end-start;
                cerr << elapsed_seconds.count() << "\t" << "+" << "\t" << best_f.sequence() << endl;
            }
        }

        if (!(prefer_forward && (float)best_f.score() / (float)sequence.size() >= target_score_per_bp))
        {
            // If we need to look on the reverse strand, do that too.
            std::chrono::time_point<std::chrono::system_clock> start, end;
            if (debug) start = std::chrono::system_clock::now();
            alignments_r = align_threaded(best_r, kmer_count_r, kmer_size, stride, attempt);
            if (debug) {
                end = std::chrono::system_clock::now();
                std::chrono::duration<double> elapsed_seconds = end-start;
                cerr << elapsed_seconds.count() << "\t" << "-" << "\t" << best_r.sequence() << endl;
            }
        }

        ++attempt;

        if (best_f.score() == 0 && best_r.score() == 0
            && kmer_size - kmer_sensitivity_step >= kmer_min) {
            // We couldn't find anything. Try harder.
            increase_sensitivity();
        } else {
            // We found at least one alignment
            break;
        }

    }

    if (debug) {
        end_both = std::chrono::system_clock::now();
        std::chrono::duration<double> elapsed_seconds = end_both-start_both;
        cerr << elapsed_seconds.count() << "\t" << "b" << "\t" << sequence << endl;
    }

    // Do a merge of the alignments up to the max limit
    vector<Alignment> merged;
    // What alignments are we looking at next in out in-order merge?
    size_t next_f = 0;
    size_t next_r = 0;
    
    // TODO: Apply a minimum score threshold?
    while(merged.size() < max_multimaps) {
        if(next_f < alignments_f.size()) {
            // We have an available forward alignment
            if(next_r < alignments_r.size()) {
                // We also have an available reverse alignment
                if(alignments_f[next_f].score() >= alignments_r[next_r].score()) {
                    // Take the forward alignment if it has a greater or equal score.
                    // TODO: this introduces a slight strand bias
                    merged.push_back(alignments_f[next_f]);
                    next_f++;
                } else {
                    // The reverse alignment is better
                    merged.push_back(alignments_r[next_r]);
                    next_r++;
                }
            } else {
                // Just take the forward alignments, since they are the only ones left
                merged.push_back(alignments_f[next_f]);
                next_f++;
            }
        } else if(next_r < alignments_r.size()) {
            // Just take the reverse alignments, since they are the only ones left
            merged.push_back(alignments_r[next_r]);
            next_r++;
        } else {
            // No alignments left, oh no!
            break;
        }
    }
    
    // Set all but the first alignment secondary.
    for(size_t i = 1; i < merged.size(); i++) {
        merged[i].set_is_secondary(true);
    }

    // Return the merged list of good alignments. Does not bother updating the input alignment.
    return merged;
}

Alignment Mapper::align(Alignment& aln, int kmer_size, int stride, int band_width) {
    // Do the multi-mapping
    vector<Alignment> best = align_multi(aln, kmer_size, stride, band_width);
    
    if(best.size() == 0) {
        // Spit back an alignment that says we failed, but make sure it has the right sequence in it.
        Alignment failed = aln;
        failed.clear_path();
        failed.set_score(0);
        return failed;
    }
    
    // Otherwise, just repoirt the best alignment, since we know one exists
    return best[0];
}

vector<Alignment> Mapper::align_threaded(Alignment& alignment, int& kmer_count, int kmer_size, int stride, int attempt) {

    // parameters, some of which should probably be modifiable
    // TODO -- move to Mapper object

    if (index == nullptr && (xindex == nullptr || gcsa == nullptr)) {
        cerr << "error:[vg::Mapper] index(es) missing, cannot map alignment!" << endl;
        exit(1);
    }

    const string& sequence = alignment.sequence();
    
    // Generate all the kmers we want to look up, with the correct stride.
    auto kmers = balanced_kmers(sequence, kmer_size, stride);

    //vector<uint64_t> sizes;
    //index->approx_sizes_of_kmer_matches(kmers, sizes);

    // Holds the map from node ID to collection of start offsets, one per kmer we're searching for.
    vector<map<int64_t, vector<int32_t> > > positions(kmers.size());
    int i = 0;
    for (auto& k : kmers) {
        if (!allATGC(k)) continue; // we can't handle Ns in this scheme
        //if (debug) cerr << "kmer " << k << " entropy = " << entropy(k) << endl;
        if (min_kmer_entropy > 0 && entropy(k) < min_kmer_entropy) continue;
        
        // We fill this in only once if we're using GCSA indexing
        gcsa::range_type gcsa_range;
        
        // Work out the number of *bytes* of matches for this kmer with the appropriate index.
        uint64_t approx_matches;
        if(gcsa) {
            // A little more complicated. We run the search and count the range size
            gcsa_range = gcsa->find(k);
            // Measure count and convert to bytes
            approx_matches = gcsa::Range::length(gcsa_range) * sizeof(gcsa::node_type);
        } else if(index) {
           approx_matches = index->approx_size_of_kmer_matches(k);
        } else {
            cerr << "error:[vg::Mapper] no search index present" << endl;
            exit(1);
        }
        
        // Report the approximate match byte size
        if (debug) cerr << k << "\t~" << approx_matches << endl;
        // if we have more than one block worth of kmers on disk, consider this kmer non-informative
        // we can do multiple mapping by relaxing this
        if (approx_matches > hit_size_threshold) {
            continue;
        }
        
        // Grab the map from node ID to kmer start positions for this particular kmer.
        auto& kmer_positions = positions.at(i);
        // Fill it in, since we know there won't be too many to work with.
        
        if(gcsa) {
            // We need to fill in this vector with the GCSA nodes and then convert.
            std::vector<gcsa::node_type> gcsa_nodes;
            gcsa->locate(gcsa_range, gcsa_nodes);
            
            for(gcsa::node_type gcsa_node : gcsa_nodes) {
                if(gcsa::Node::rc(gcsa_node)) {
                    // We found a kmer on the reverse strand. The old index
                    // didn't handle these, so we ignore them. TODO: figure out
                    // how to account for them.
                    continue;
                }
                // Decode the result's ID and offset and record it
                kmer_positions[gcsa::Node::id(gcsa_node)].push_back(gcsa::Node::offset(gcsa_node));
            }
            
        } else if(index) {
           index->get_kmer_positions(k, kmer_positions);
        } else {
            cerr << "error:[vg::Mapper] no search index present" << endl;
            exit(1);
        }
        
        
        // ignore this kmer if it has too many hits
        // typically this will be filtered out by the approximate matches filter
        if (kmer_positions.size() > hit_max) kmer_positions.clear();
        // Report the actual match count for the kmer
        if (debug) cerr << "\t=" << kmer_positions.size() << endl;
        kmer_count += kmer_positions.size();
        // break when we get more than a threshold number of kmers to seed further alignment
        //if (kmer_count >= kmer_threshold) break;
        ++i;
    }

    if (debug) cerr << "kept kmer hits " << kmer_count << endl;

    // make threads
    // these start whenever we have a kmer match which is outside of
    // one of the last positions (for the previous kmer) + the kmer stride % wobble (hmm)

    // For each node ID, holds the numbers of the kmers that we find on it, in
    // the order that they appear in the query. One would expect them to be
    // monotonically increasing.
    map<int64_t, vector<int> > node_kmer_order;
    
    // Maps from node ID and offset to a thread ending with the kmer that starts
    // there, if any such thread exists.
    map<pair<int64_t, int32_t>, vector<int64_t> > position_threads;
    
    // For each node, holds the last thread for that node. Because we only do
    // position wobble, threads only touch a single node.
    map<int64_t, vector<int64_t> > node_threads;
    
    //int node_wobble = 0; // turned off...
    
    // How far left or right from the "correct" position for the previous kmer
    // are we willing to search when looking for a thread to extend?
    int position_wobble = 2;

    int max_iter = sequence.size();
    int iter = 0;
    int64_t max_subgraph_size = 0;

    // This is basically the index on this loop over kmers and their position maps coming up
    i = 0;
    for (auto& p : positions) {
        // For every map from node ID to collection of kmer starts, for kmer i...
        
        // Grab the kmer and advance i for next loop iteration
        auto& kmer = kmers.at(i++);
        for (auto& x : p) {
            // For each node ID and the offsets on that node at which this kmer appears...        
            int64_t id = x.first;
            vector<int32_t>& pos = x.second;
            
            // Note that this kmer is the next kmer in the query to appear in that node.
            node_kmer_order[id].push_back(i-1);
            for (auto& y : pos) {
                // For every offset along the node at which this kmer appears, in order...
            
                //cerr << kmer << "\t" << i << "\t" << id << "\t" << y << endl;
                // thread rules
                // if we find the previous position
                
                // This holds the thread that this instance of this kmer on this node is involved in.
                vector<int64_t> thread;

                // If we can find a thread close enough to this kmer, we want to
                // continue it with this kmer. If nothing changed between the
                // query and the reference, we would expect to extend the thread
                // that has its last kmer starting exactly stride bases before
                // this kmer starts (i.e. at y - stride). However, due to indels
                // existing, we search with a "wobble" of up to position_wobble
                // in either direction, outwards from the center.
                
                // This holds the current wobble that we are searching (between
                // -position_wobble and +position_wobble).
                int m = 0;
                for (int j = 0; j < 2*position_wobble + 1; ++j) {
                    // For each of the 2 * position_wobble + 1 wobble values we
                    // need to try, calculate the jth wobble value out from the
                    // center.
                    if (j == 0) { // on point
                        // First we use the zero wobble, which we started with
                    } else if (j % 2 == 0) { // subtract
                        // Every even step except the first, we try the negative version of the positive wobble we just tried
                        m *= -1;
                    } else { // add
                        // Every odd step, we try the positive version of the
                        // negative (or 0) wobble we just tried, incremented by
                        // 1.
                        m *= -1; ++m;
                    }
                    
                    //cerr << "checking " << id << " " << y << " - " << kmer_size << " + " << m << endl;
                    
                    // See if we can find a thread at this wobbled position
                    auto previous = position_threads.find(make_pair(id, y - stride + m));
                    if (previous != position_threads.end()) {
                        // If we did find one, use it as our thread, remove it
                        // so it can't be extended by anything else, and stop
                        // searching more extreme wobbles.
                        
                        //length = position_threads[make_pair(id, y - stride + m)] + 1;
                        thread = previous->second;
                        position_threads.erase(previous);
                        //cerr << "thread is " << thread.size() << " long" << endl;
                        break;
                    }
                }

                // Now we either have the thread we are extending in thread, or we are starting a new thread.
                
                // Extend the thread with another kmer on this node ID. 
                thread.push_back(id);
                // Save the thread as ending with a kmer at this offset on this node.
                position_threads[make_pair(id, y)] = thread;
                
                // This is now the last thread for this node.
                node_threads[id] = thread;
            }
        }
    }

    // This maps from a thread length (in kmer instances) to all the threads of that length.
    map<int, vector<vector<int64_t> > > threads_by_length;
    for (auto& t : node_threads) {
        auto& thread = t.second;
        auto& threads = threads_by_length[thread.size()];
        threads.push_back(thread);
    }

    // now sort the threads and re-cluster them

    if (debug) {
        cerr << "initial threads" << endl;
        for (auto& t : threads_by_length) {
            auto& length = t.first;
            auto& threads = t.second;
            cerr << length << ":" << endl;
            for (auto& thread : threads) {
                cerr << "\t";
                for (auto& id : thread) {
                    cerr << id << " ";
                }
                cerr << endl;
            }
            cerr << endl;
        }
    }

    // sort threads by ids, taking advantage of vector comparison and how sets work
    set<vector<int64_t> > sorted_threads;
    auto tl = threads_by_length.rbegin();
    for (auto& t : node_threads) {
        auto& thread = t.second;
        sorted_threads.insert(thread);
    }
    threads_by_length.clear();

    // go back through and combine closely-linked threads
    // ... but only if their kmer order is proper
    
    // This holds threads by the last node ID they touch.
    map<int64_t, vector<int64_t> > threads_by_last;
    
    // go from threads that are longer to ones that are shorter
    for (auto& thread : sorted_threads) {
        //cerr << thread.front() << "-" << thread.back() << endl;
        
        // Find the earliest-ending thread that ends within max_thread_gap nodes of this thread's start
        auto prev = threads_by_last.upper_bound(thread.front()-max_thread_gap);
        //if (prev != threads_by_last.begin()) --prev;
        // now we should be at the highest thread within the bounds
        //cerr << prev->first << " " << thread.front() << endl;
        // todo: it may also make sense to check that the kmer order makes sense
        // what does this mean? it means that the previous 
        if (prev != threads_by_last.end()
            && prev->first > thread.front() - max_thread_gap) {
            // If we found such a thread, and it also *starts* within
            // max_thread_gap nodes of this thread's start, we want to add our
            // thread onto the end of it and keep only the combined longer
            // thread. TODO: this limits max thread length.
            vector<int64_t> new_thread;
            auto& prev_thread = prev->second;
            new_thread.reserve(prev_thread.size() + thread.size());
            new_thread.insert(new_thread.end(), prev_thread.begin(), prev_thread.end());
            new_thread.insert(new_thread.end(), thread.begin(), thread.end());
            threads_by_last.erase(prev);
            // this will clobber... not good
            // maybe overwrite only if longer?
            threads_by_last[new_thread.back()] = new_thread;
        } else {
            // We want to keep this thread since it couldn't attach to any other thread.
            threads_by_last[thread.back()] = thread;
        }
    }

    // debugging
    /*
    if (debug) {
        cerr << "threads by last" << endl;
        for (auto& t : threads_by_last) {
            auto& thread = t.second;
            cerr << t.first << "\t";
            for (auto& id : thread) {
                cerr << id << " ";
            }
            cerr << endl;
        }
    }
    */

    // rebuild our threads_by_length set
    for (auto& t : threads_by_last) {
        auto& thread = t.second;
        if (thread.size() >= cluster_min) {
            // Only keep threads if they have a sufficient number of kmer instances in them.
            auto& threads = threads_by_length[thread.size()];
            threads.push_back(thread);
        }
    }

    if (debug) {
        cerr << "threads ready for alignment" << endl;
        for (auto& t : threads_by_length) {
            auto& length = t.first;
            auto& threads = t.second;
            cerr << length << ":" << endl;
            for (auto& thread : threads) {
                cerr << "\t";
                for (auto& id : thread) {
                    cerr << id << " ";
                }
                cerr << endl;
            }
            cerr << endl;
        }
    }

    int thread_ex = thread_extension;
    map<vector<int64_t>*, Alignment> alignments;

    // collect the nodes from the best N threads by length
    // and expand subgraphs as before
    //cerr << "extending by " << thread_ex << endl;
    tl = threads_by_length.rbegin();
    bool accepted = false;
    for (int i = 0;
         !accepted
             && tl != threads_by_length.rend()
             && (best_clusters == 0 || i < best_clusters);
         ++i, ++tl) {
        auto& threads = tl->second;
        // by definition, our thread should construct a contiguous graph
        for (auto& thread : threads) {
            // Do an alignment to the subgraph for each thread.
        
            // thread extension should be determined during iteration
            // note that there is a problem and hits tend to be imbalanced
            int64_t first = max((int64_t)0, *thread.begin() - thread_ex);
            int64_t last = *thread.rbegin() + thread_ex;
            //int64_t first = *thread.begin();
            //int64_t last = *thread.rbegin();
            // so we can pick it up efficiently from the index by pulling the range from first to last
            if (debug) cerr << "getting node range " << first << "-" << last << endl;
            VG* graph = new VG;
            
            // Now we need to get the neighborhood by ID and expand outward by actual
            // edges. How we do this depends on what indexing structures we have.
            // TODO: We're repeating this code. Break it out into a function or something.
            if(xindex) {
                xindex->get_id_range(first, last, graph->graph);
                xindex->expand_context(graph->graph, context_depth);
            } else if(index) {
                index->get_range(first, last, *graph);
                index->expand_context(*graph, context_depth);
            } else {
                cerr << "error:[vg::Mapper] cannot align mate with no graph data" << endl;
                exit(1);
            }
            
            Alignment& ta = alignments[&thread];
            ta = alignment;
            // by default, expand the graph a bit so we are likely to map
            //index->get_connected_nodes(*graph);
            graph->remove_orphan_edges();

            if (debug) cerr << "got subgraph with " << graph->node_count() << " nodes, " 
                            << graph->edge_count() << " edges" << endl;
                            
            // Topologically sort the graph, breaking cycles and orienting all edges end to start.
            // This flips some nodes around, so we need to translate alignments back.
            set<int64_t> flipped_nodes;
            graph->orient_nodes_forward(flipped_nodes);
                            
            // align
            //graph->serialize_to_file("align2.vg");
            ta.clear_path();
            ta.set_score(0);
            graph->align(ta);

            // check if we start or end with soft clips
            // if so, try to expand the graph until we don't have any more (or we hit a threshold)
            // expand in the direction where there were soft clips

            if (!ta.has_path()) continue;
            
            int sc_start = softclip_start(ta);
            int sc_end = softclip_end(ta);

            if (sc_start > softclip_threshold || sc_end > softclip_threshold) {
                if (debug) cerr << "softclip handling " << sc_start << " " << sc_end << endl;
                Path* path = ta.mutable_path();
                int64_t idf = path->mutable_mapping(0)->position().node_id();
                int64_t idl = path->mutable_mapping(path->mapping_size()-1)->position().node_id();
                // step towards the side where there were soft clips
                // using 10x the thread_extension
                int64_t f = max((int64_t)0, idf - (int64_t) max(thread_ex, 1) * 10);
                int64_t l = idl + (int64_t) max(thread_ex, 1) * 10;
                
                // We're going to get three ranges. They need to be non-
                // overlapping since the xg range operation can't handle
                // repeated nodes.
                // The ranges are:
                // The original first-last range
                // The new range suggested by soft clip handling on the left (f to min(idf, first - 1))
                // The new range suggested by soft clip handling on the right (max(idl, last + 1) to l)
                // This way, if soft clip sends us off across an ID discontinuity, we don't try to get like half the graph.
                
                // The last two entries might be empty or backward, but the
                // range functions can just not get anything in those cases.
                
                if (debug) {
                    cerr << "getting node ranges: " << endl;
                    cerr << "\t" << first << "-" << last << endl;
                    cerr << "\t" << f << "-" << min(idf, first - 1) << endl;
                    cerr << "\t" << max(idl, last + 1) << "-" << l << endl;
                }

                // always rebuild the graph, since we messed it up by sorting it
                delete graph;
                graph = new VG;
                
                
                // Get the bigger range, but still go out to context depth afterwards.
                if(xindex) {
                    xindex->get_id_range(first, last, graph->graph);
                    xindex->get_id_range(f, min(idf, first - 1), graph->graph);
                    xindex->get_id_range(max(idl, last + 1), l, graph->graph);
                    xindex->expand_context(graph->graph, context_depth);
                    graph->rebuild_indexes();
                } else if(index) {
                    index->get_range(first, last, *graph);
                    index->get_range(f, min(idf, first - 1), *graph);
                    index->get_range(max(idl, last + 1), l, *graph);
                    index->expand_context(*graph, context_depth);
                } else {
                    cerr << "error:[vg::Mapper] cannot align mate with no graph data" << endl;
                    exit(1);
                }
                
                graph->remove_orphan_edges();
                
                if (debug) cerr << "got subgraph with " << graph->node_count() << " nodes, " 
                                << graph->edge_count() << " edges" << endl;

                //graph->serialize_to_file("align1.vg");
                ta.clear_path();
                ta.set_score(0);
                graph->align(ta);
                if (debug) cerr << "softclip after " << softclip_start(ta) << " " << softclip_end(ta) << endl;
            }

            delete graph;
            
            if (debug) cerr << "score per bp is " << (float)ta.score() / (float)ta.sequence().size() << endl;
            if (greedy_accept && (float)ta.score() / (float)ta.sequence().size() >= target_score_per_bp) {
                if (debug) cerr << "greedy accept" << endl;
                accepted = true;
                break;
            }
        }
    }

    // now find the best alignment
    int sum_score = 0;
    double mean_score = 0;
    map<int, set<Alignment*> > alignment_by_score;
    for (auto& ta : alignments) {
        Alignment* aln = &ta.second;
        alignment_by_score[aln->score()].insert(aln);
    }

    // Collect all the good alignments
    // TODO: Filter down subject to a minimum score per base or something?
    vector<Alignment> good;
    for(auto it = alignment_by_score.rbegin(); it != alignment_by_score.rend(); ++it) {
        // Copy over all the alignments in descending score order (following the pointers into the "alignments" vector)
        // Iterating through a set keyed on ints backward is in descending order.

        // This is going to let us deduplicate our alignments with this score, by storing them serialized to strings in this set.
        set<string> serializedAlignmentsUsed;
        
        for(Alignment* pointer : (*it).second) {
            // We serialize the alignment to a string
            string serialized;
            pointer->SerializeToString(&serialized);
        
            if(!serializedAlignmentsUsed.count(serialized)) {
                // This alignment hasn't been produced yet. Produce it. The
                // order in the alignment vector doesn't matter for things with
                // the same score.
                good.push_back(*pointer);
                
                // Save it so we can avoid putting it in the vector again
                serializedAlignmentsUsed.insert(serialized);
            }
            
            if(good.size() >= max_multimaps) {
                // Don't report too many mappings
                break;
            }
        }
        
        if(good.size() >= max_multimaps) {
            // Don't report too many mappings
            break;
        }
    }

    // get the best alignment
    if (!good.empty()) {
        alignment = good[0];
        if (debug) {
            cerr << "best alignment score " << alignment.score() << endl;
        }
    } else {
        alignment.clear_path();
        alignment.set_score(0);
        
        // We return an empty vector, and give an Alignment with no path and 0
        // score through our output parameter.
    }

    if (debug && alignment.score() == 0) cerr << "failed alignment" << endl;

    // Return all the multimappings
    return good;
}

Alignment& Mapper::align_simple(Alignment& alignment, int kmer_size, int stride) {

    if (index == NULL) {
        cerr << "error:[vg::Mapper] no index loaded, cannot map alignment!" << endl;
        exit(1);
    }

    // establish kmers
    const string& sequence = alignment.sequence();
    //  
    auto kmers = balanced_kmers(sequence, kmer_size, stride);

    map<string, int32_t> kmer_counts;
    vector<map<int64_t, vector<int32_t> > > positions(kmers.size());
    int i = 0;
    for (auto& k : kmers) {
        index->get_kmer_positions(k, positions.at(i++));
        kmer_counts[k] = positions.at(i-1).size();
    }
    positions.clear();
    VG* graph = new VG;
    for (auto& c : kmer_counts) {
        if (c.second < hit_max) {
            index->get_kmer_subgraph(c.first, *graph);
        }
    }

    int max_iter = sequence.size();
    int iter = 0;
    int context_step = 1;
    int64_t max_subgraph_size = 0;

    // use kmers which are informative
    // and build up the graph

    auto get_max_subgraph_size = [this, &max_subgraph_size, &graph]() {
        list<VG> subgraphs;
        graph->disjoint_subgraphs(subgraphs);
        for (auto& subgraph : subgraphs) {
            max_subgraph_size = max(subgraph.total_length_of_nodes(), max_subgraph_size);
        }
    };

    get_max_subgraph_size();

    while (max_subgraph_size < sequence.size()*2 && iter < max_iter) {
        index->expand_context(*graph, context_step);
        index->get_connected_nodes(*graph);
        get_max_subgraph_size();
        ++iter;
    }
    // ensure we have a complete graph prior to alignment
    index->get_connected_nodes(*graph);

    /*
    ofstream f("vg_align.vg");
    graph->serialize_to_ostream(f);
    f.close();
    */

    // Make sure the graph we're aligning to is all oriented
    graph->align(alignment);
    
    delete graph;

    return alignment;

}

const int balanced_stride(int read_length, int kmer_size, int stride) {
    double r = read_length;
    double k = kmer_size;
    double j = stride;
    int i = (r > j) ? round((r-k)/round((r-k)/j)) : j;
    return max(1, i);
}

const vector<string> balanced_kmers(const string& seq, const int kmer_size, const int stride) {
    // choose the closest stride that will generate balanced kmers
    vector<string> kmers;
    int b = balanced_stride(seq.size(), kmer_size, stride);
    if (!seq.empty()) {
        for (int i = 0; i+kmer_size <= seq.size(); i+=b) {
            kmers.push_back(seq.substr(i,kmer_size));
        }
    }
    return kmers;
}

}
