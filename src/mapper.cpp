#include "mapper.hpp"

namespace vg {

Mapper::Mapper(Index* idex,
               xg::XG* xidex,
               gcsa::GCSA* g,
               gcsa::LCPArray* a)
    : index(idex)
    , xindex(xidex)
    , gcsa(g)
    , lcp(a)
    , best_clusters(0)
    , cluster_min(2)
    , hit_max(0)
    , hit_size_threshold(512)
    , kmer_min(0)
    , kmer_sensitivity_step(3)
    , thread_extension(10)
    , max_thread_gap(30)
    , context_depth(1)
    , max_multimaps(1)
    , max_attempts(0)
    , softclip_threshold(0)
    , prefer_forward(false)
    , greedy_accept(false)
    , accept_norm_score(0.75)
    , min_norm_score(0)
    , min_kmer_entropy(0)
    , debug(false)
    , alignment_threads(1)
    , max_mem_length(0)
    , min_mem_length(0)
    , max_target_factor(100)
    , match(2)
    , mismatch(2)
    , gap_open(3)
    , gap_extension(1)
{
    // Nothing to do. We just hold the default parameter values.
}

Mapper::Mapper(Index* idex, gcsa::GCSA* g, gcsa::LCPArray* a) : Mapper(idex, nullptr, g, a)
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

Mapper::Mapper(xg::XG* xidex, gcsa::GCSA* g, gcsa::LCPArray* a) : Mapper(nullptr, xidex, g, a) {
    if(xidex == nullptr) {
        // With this constructor we need an XG graph.
        cerr << "error:[vg::Mapper] cannot create an xg-based Mapper with null xg index" << endl;
        exit(1);
    }

    if(g == nullptr || a == nullptr) {
        // With this constructor we need a GCSA2 index too.
        cerr << "error:[vg::Mapper] cannot create an xg-based Mapper with null GCSA2 index" << endl;
        exit(1);
    }
}

Mapper::Mapper(void) : Mapper(nullptr, nullptr, nullptr, nullptr) {
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
void Mapper::align_mate_in_window(const Alignment& read1, Alignment& read2, int pair_window) {
    if (read1.score() == 0) return; // bail out if we haven't aligned the first
    // try to recover in region
    auto& path = read1.path();
    int64_t idf = path.mapping(0).position().node_id();
    int64_t idl = path.mapping(path.mapping_size()-1).position().node_id();
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
        // don't get the paths (this isn't yet threadsafe in sdsl-lite)
        xindex->expand_context(graph->graph, context_depth, false);
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

    read2 = graph->align(read2, match, mismatch, gap_open, gap_extension);
    delete graph;
}

pair<vector<Alignment>, vector<Alignment>> Mapper::align_paired_multi(
    const Alignment& read1, const Alignment& read2, int kmer_size, int stride, int band_width, int pair_window) {

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
        if (read.has_path() && read.path().mapping(0).position().is_reverse()) {
            // If so, reverse the "same direction" Alignment sequence
            aln_same.set_sequence(reverse_complement(aln_same.sequence()));
        } else {
            // Otherwise reverse the opposite direction sequence
            aln_opposite.set_sequence(reverse_complement(aln_opposite.sequence()));
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

pair<Alignment, Alignment> Mapper::align_paired(const Alignment& read1, const Alignment& read2, int kmer_size, int stride,
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

// We need a function to get the lengths of nodes, in case we need to
// reverse an Alignment, including all its Mappings and Positions. TODO:
// make this cache the node lengths for the nodes used in the actual
// alignments somehow?
int64_t Mapper::get_node_length(int64_t node_id) {
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
}

Alignment Mapper::align_banded(const Alignment& read, int kmer_size, int stride, int band_width) {
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
    vector<pair<size_t, size_t>> to_strip; to_strip.resize(to_align);
    vector<Alignment> bands; bands.resize(to_align);

    // scan across the read choosing bands
    // these bands are hard coded to overlap by 50%
    // the last band is guaranteed to be segment_size long
    // overlap scheme example
    // read: ----------------------
    //       --------
    //           --------
    //               --------
    //                   --------
    //                     --------
    // Afterwards, we align each fragment, trim the overlaps implied by the layout
    // and concatenate the alignments. The result is a split read alignment that
    // can describe large indels, CNVs, and inversions natively, even though our
    // local alignment algorithm is only aware of alignment against DAGs.
    for (int i = 0; i < div; ++i) {
        size_t off = i*segment_size;
        auto aln = read;
        size_t addl_seq = 0;
        if (i+1 == div) {
            // ensure we have a full-length segment for the last alignment
            // otherwise we run the risk of trying to align a very tiny band
            size_t last_off = read.sequence().size() - segment_size;
            if (off > last_off) {
                // looks wrawng
                addl_seq = (off - last_off);
                aln.set_sequence(read.sequence().substr(last_off));
                //assert(aln.sequence().size() == segment_size);
            } else {
                aln.set_sequence(read.sequence().substr(off));
            }
        } else {
            aln.set_sequence(read.sequence().substr(off, segment_size));
        }
        size_t idx = 2*i;
        to_strip[idx].first = (i == 0 ? 0 : segment_size/4 + addl_seq);
        to_strip[idx].second = (i+1 == div ? 0 : segment_size/4);
        bands[idx] = aln;
        if (i != div-1) { // if we're not at the last sequence
            aln.set_sequence(read.sequence().substr(off+segment_size/2,
                                                    segment_size));
            idx = 2*i+1;
            to_strip[idx].first = segment_size/4;
            // record second but take account of case where we run off end
            to_strip[idx].second = segment_size/4 - (segment_size - aln.sequence().size());
            bands[idx] = aln;
        }
    }

    vector<vector<Alignment>> multi_alns;
    vector<Alignment> alns;
    if (max_multimaps > 1) multi_alns.resize(to_align);
    else alns.resize(to_align);

    auto check_alignment = [&](const Alignment& aln) {
        // use the graph to extract the sequence
        // assert that this == the alignment
        if (aln.path().mapping_size()) {
            // get the graph corresponding to the alignment path
            Graph sub;
            for (int i = 0; i < aln.path().mapping_size(); ++ i) {
                auto& m = aln.path().mapping(i);
                if (m.has_position() && m.position().node_id()) {
                    auto id = aln.path().mapping(i).position().node_id();
                    // XXXXXX this is single-threaded!
                    xindex->neighborhood(id, 2, sub);
                }
            }
            VG g; g.extend(sub);
            auto seq = g.path_string(aln.path());
            //if (aln.sequence().find('N') == string::npos && seq != aln.sequence()) {
            if (seq != aln.sequence()) {
                cerr << "alignment does not match graph " << endl
                     << pb2json(aln) << endl
                     << "expect:\t" << aln.sequence() << endl
                     << "got:\t" << seq << endl;
                // save alignment
                write_alignment_to_file(aln, "fail.gam");
                // save graph, bigger fragment
                xindex->expand_context(sub, 5, true);
                VG gn; gn.extend(sub);
                gn.serialize_to_file("fail.vg");
                assert(false);
            }
        }
    };

#pragma omp parallel for schedule(dynamic,1)
    for (int i = 0; i < bands.size(); ++i) {
        if (max_multimaps > 1) {
            vector<Alignment>& malns = multi_alns[i];
            malns = align_multi(bands[i], kmer_size, stride);
            // always include an unaligned mapping
            malns.push_back(bands[i]);
            for (vector<Alignment>::iterator a = malns.begin(); a != malns.end(); ++a) {
                Alignment& aln = *a;
                bool above_threshold = ((float) aln.score() / ((float) aln.sequence().size() * match)
                                        >= min_norm_score);
                if (!above_threshold) {
                    // treat as unmapped
                    aln = bands[i];
                }
                // strip overlaps
                aln = strip_from_start(aln, to_strip[i].first);
                aln = strip_from_end(aln, to_strip[i].second);
            }
        } else {
            Alignment& aln = alns[i];
            aln = align(bands[i], kmer_size, stride);
            bool above_threshold = ((float) aln.score() / ((float) aln.sequence().size() * match)
                                    >= min_norm_score);
            if (!above_threshold) {
                aln = bands[i]; // unmapped
            }
            // strip overlaps
            //cerr << "checking before strip" << endl;
            //check_alignment(aln);
            aln = strip_from_start(aln, to_strip[i].first);
            aln = strip_from_end(aln, to_strip[i].second);
            //cerr << "checking after strip" << endl;
            //check_alignment(aln);
            //cerr << "OK" << endl;
        }
    }

    // resolve the highest-scoring traversal of the multi-mappings
    if (max_multimaps > 1) {
        alns = resolve_banded_multi(multi_alns);
        multi_alns.clear(); // clean up
    }

    // check that the alignments are valid
    if (debug) {
        for (auto& aln : alns) {
            check_alignment(aln);
        }
    }

    // merge the resulting alignments
    Alignment merged = merge_alignments(alns, debug);
    // TODO RECALCULATE QUALITY BASED ON SCORING
    merged.set_quality(read.quality());
    merged.set_name(read.name());

    if(debug) {
        for(int i = 0; i < merged.path().mapping_size(); i++) {
            // Check each Mapping to make sure it doesn't go past the end of its
            // node.
            auto& mapping = merged.path().mapping(i);

            // What node is the mapping on
            int64_t node_id = mapping.position().node_id();
            if(node_id != 0) {
                // If it's actually on a node, get the node's sequence length
                int64_t node_length = get_node_length(node_id);

                // Make sure the mapping from length is shorter than the node length
                assert(node_length >= mapping_from_length(mapping));
            }
        }
    }

    return merged;
}

vector<Alignment> Mapper::resolve_banded_multi(vector<vector<Alignment>>& multi_alns) {
    // use a basic dynamic programming to score the path through the multi mapping
    // we add the score as long as our alignments overlap (we expect them to)
    // otherwise we add nothing
    // reads that are < the minimum alignment score threshold are dropped

    // a vector of
    // score, current alignment, parent alignment (direction)
    typedef tuple<int, Alignment*, size_t> score_t;
    vector<vector<score_t>> scores;
    scores.resize(multi_alns.size());
    // start with the scores for the first alignments
    if (debug) {
        cerr << "resolving banded multi over:" << endl;
        for (auto& alns : multi_alns) {
            for (auto& aln : alns) {
                if (aln.has_path()) {
                    cerr << aln.score() << "@ " << make_pos_t(aln.path().mapping(0).position()) <<", ";
                }
            }
            cerr << endl;
        }
    }
    for (auto& aln : multi_alns[0]) {
        scores.front().push_back(make_tuple(aln.score(), &aln, 0));
    }
    for (size_t i = 1; i < multi_alns.size(); ++i) {
        auto& curr_alns = multi_alns[i];
        vector<score_t>& curr_scores = scores[i];
        auto& prev_scores = scores[i-1];
        // find the best previous score
        score_t best_prev = prev_scores.front();
        size_t best_idx = 0;
        score_t unmapped_prev = prev_scores.front();
        size_t unmapped_idx = 0;
        size_t j = 0;
        for (auto& t : prev_scores) {
            if (get<0>(t) > get<0>(best_prev)) {
                best_prev = t;
                best_idx = j;
            }
            if (get<0>(t) == 0) {
                unmapped_idx = j;
                unmapped_prev = t;
            }
            ++j;
        }
        // for each alignment
        for (auto& aln : curr_alns) {
            // if it's not mapped, take the best previous score
            if (!aln.score()) {
                curr_scores.push_back(make_tuple(get<0>(best_prev),
                                                 &aln, best_idx));
            } else {
                // determine our start
                auto& curr_start = aln.path().mapping(0).position();
                // accumulate candidate alignments
                map<int, vector<pair<score_t, size_t>>> candidates;
                // for each previous alignment
                size_t k = 0;
                for (auto& score : prev_scores) {
                    auto old = get<1>(score);
                    if (!old->score()) continue; // unmapped
                    auto prev_end = path_end(old->path());
                    // save it as a candidate if the two are adjacent
                    // and in the same orientation
                    if (adjacent_positions(prev_end, curr_start)) {
                        candidates[get<0>(score)].push_back(make_pair(score,k));
                    }
                    ++k;
                }
                if (candidates.size()) {
                    // take the best one (at least the first best one we saw)
                    auto& opt = candidates.rbegin()->second.front();
                    // DP scoring step: add scores when we match head to tail
                    curr_scores.push_back(make_tuple(get<0>(opt.first) + aln.score(),
                                                     &aln, opt.second));
                } else {
                    // if there are no alignments matching our start
                    // just take the highest-scoring one
                    auto best_prev_aln = get<1>(prev_scores[best_idx]);
                    if (best_prev_aln->has_path()) {
                        curr_scores.push_back(make_tuple(get<0>(best_prev),
                                                         &aln, best_idx));
                    } else {
                        curr_scores.push_back(make_tuple(get<0>(unmapped_prev),
                                                         &aln, unmapped_idx));
                    }
                }
            }
        }
    }
    // find the best score at the end
    score_t best_last = scores.back().front();
    size_t best_last_idx = 0;
    size_t j = 0;
    for (auto& s : scores.back()) {
        if (get<0>(s) > get<0>(best_last)) {
            best_last = s;
            best_last_idx = j;
        }
        ++j;
    }
    // accumulate the alignments in the optimal path
    vector<Alignment> alns; alns.resize(multi_alns.size());
    size_t prev_best_idx = best_last_idx;
    for (int i = scores.size()-1; i >= 0; --i) {
        auto& score = scores[i][prev_best_idx];
        alns[i] = *get<1>(score); // save the alignment
        prev_best_idx = get<2>(score); // and where we go next
    }
    return alns;
}

bool Mapper::adjacent_positions(const Position& pos1, const Position& pos2) {
    // are they the same id, with offset differing by 1?
    if (pos1.node_id() == pos2.node_id()
        && pos1.offset() == pos2.offset()-1) {
        return true;
    }
    // otherwise, we're going to need to check via the index
    VG graph;
    // pick up a graph that's just the neighborhood of the start and end positions
    int64_t id1 = pos1.node_id();
    int64_t id2 = pos2.node_id();
    if(xindex) {
        // Grab the node sequence only from the XG index and get its size.
        xindex->get_id_range(id1, id1, graph.graph);
        xindex->get_id_range(id2, id2, graph.graph);
        xindex->expand_context(graph.graph, 1, false);
        graph.rebuild_indexes();
    } else if(index) {
        index->get_context(id1, graph);
        index->get_context(id2, graph);
        index->expand_context(graph, 1);
    } else {
        throw runtime_error("No index to get nodes from.");
    }
    // now look in the graph to figure out if we are adjacent
    return graph.adjacent(pos1, pos2);
}

vector<Alignment> Mapper::align_multi(const Alignment& aln, int kmer_size, int stride, int band_width) {
    // trigger a banded alignment if we need to
    // note that this will in turn call align_multi on fragments of the read
    if (aln.sequence().size() > band_width) {
        if (debug) cerr << "switching to banded alignment" << endl;
        return vector<Alignment>{align_banded(aln, kmer_size, stride, band_width)};
    }
    if (kmer_size || index != nullptr) {
        // if we've defined a kmer size, use the legacy style mapper
        return align_multi_kmers(aln, kmer_size, stride, band_width);
    } else {
        // otherwise use the mem mapper, which is a banded multi mapper by default
        return align_mem(aln);
    }
}

vector<Alignment> Mapper::align_multi_kmers(const Alignment& aln, int kmer_size, int stride, int band_width) {

    std::chrono::time_point<std::chrono::system_clock> start_both, end_both;
    if (debug) start_both = std::chrono::system_clock::now();
    const string& sequence = aln.sequence();

    // we assume a kmer size to be specified
    if (!kmer_size && !kmer_sizes.empty()) {
        // basically assumes one kmer size
        kmer_size = *kmer_sizes.begin();
    }
    assert(kmer_size);
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
    Alignment best_r = reverse_complement_alignment(aln,
                                                    (function<int64_t(int64_t)>) ([&](int64_t id) { return get_node_length(id); }));
    // This will hold all of the reverse alignments up to max_multimaps
    vector<Alignment> alignments_r;

    auto increase_sensitivity = [this,
                                 &kmer_size,
                                 &stride,
                                 &sequence,
                                 &best_f,
                                 &best_r]() {
        kmer_size -= kmer_sensitivity_step;
        stride = sequence.size() / ceil( (double)sequence.size() / kmer_size);
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

    while (!((float)best_f.score()/((float)sequence.size()*match) > min_norm_score
             || (float)best_r.score()/((float)sequence.size()*match) > min_norm_score)
           && attempt < max_attempts) {

        //cerr << "min score per bp " << min_score_per_bp << " " << (float)best_f.score()/(float)sequence.size()
        //<< " " << (float)best_r.score()/(float)sequence.size() << endl;
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

        if (!(prefer_forward && (float)best_f.score() / ((float)sequence.size()*match) >= accept_norm_score))
        {
            // If we need to look on the reverse strand, do that too.
            std::chrono::time_point<std::chrono::system_clock> start, end;
            if (debug) start = std::chrono::system_clock::now();
            auto alns =  align_threaded(best_r, kmer_count_r, kmer_size, stride, attempt);
            alignments_r = reverse_complement_alignments(alns,
                                                         (function<int64_t(int64_t)>) ([&](int64_t id) { return get_node_length(id); }));
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
		//TODO So merged contains qualities up to this point.
		// but once it returns to main_map, there will be no qualities...

    // Return the merged list of good alignments. Does not bother updating the input alignment.
    return merged;
}



Alignment Mapper::align(const Alignment& aln, int kmer_size, int stride, int band_width) {
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

set<pos_t> gcsa_nodes_to_positions(const vector<gcsa::node_type>& nodes) {
    set<pos_t> positions;
    for(gcsa::node_type node : nodes) {
        positions.insert(make_pos_t(gcsa::Node::id(node),
                                    gcsa::Node::rc(node),
                                    gcsa::Node::offset(node)));
    }
    return positions;    
}

set<pos_t> Mapper::sequence_positions(const string& seq) {
    gcsa::range_type gcsa_range = gcsa->find(seq);
    std::vector<gcsa::node_type> gcsa_nodes;
    gcsa->locate(gcsa_range, gcsa_nodes);
    return gcsa_nodes_to_positions(gcsa_nodes);
}

// Use the GCSA2 index to find super-maximal exact matches.
vector<MaximalExactMatch>
Mapper::find_smems(const string& seq) {

    string::const_iterator string_begin = seq.begin();
    string::const_iterator cursor = seq.end();
    vector<MaximalExactMatch> mems;

    // an empty sequence matches the entire bwt
    if (seq.empty()) {
        mems.emplace_back(
            MaximalExactMatch(seq.begin(), seq.end(),
                              gcsa::range_type(0, gcsa->size() - 1)));
        return mems;
    }

    // find MEMs using GCSA+LCP array
    // algorithm sketch:
    // set up a cursor pointing to the last position in the sequence
    // set up a structure to track our MEMs, and set it == "" and full range match
    // while our cursor is >= the beginning of the string
    //   try a step of backwards searching using LF mapping
    //   if our range goes to 0
    //       go back to the last non-empty range
    //       emit the MEM corresponding to this range
    //       start a new mem
    //           use the LCP array's parent function to cut off the end of the match
    //           (effectively, this steps up the suffix tree)
    //           and calculate the new end point using the LCP of the parent node
    // emit the final MEM, if we finished in a matching state

    // the temporary MEM we'll build up in this process
    MaximalExactMatch match(cursor, cursor, gcsa::range_type(0, gcsa->size() - 1));
    gcsa::range_type last_range = match.range;
    --cursor; // start off looking at the last character in the query
    while (cursor >= string_begin) {
        // hold onto our previous range
        last_range = match.range;
        // execute one step of LF mapping
        match.range = gcsa->LF(match.range, gcsa->alpha.char2comp[*cursor]);
        if (gcsa::Range::empty(match.range)
            || max_mem_length && match.end-cursor > max_mem_length
            || match.end-cursor > gcsa->order()) {
            // we've exhausted our BWT range, so the last match range was maximal
            // or: we have exceeded the order of the graph (FPs if we go further)
            //     we have run over our parameter-defined MEM limit
            // record the last MEM
            match.begin = cursor+1;
            match.range = last_range;
            mems.push_back(match);
            // length of last MEM, which we use to update our end pointer for the next MEM
            size_t last_mem_length = match.end - match.begin;
            // get the parent suffix tree node corresponding to the parent of the last MEM's STNode
            gcsa::STNode parent = lcp->parent(last_range);
            // change the end for the next mem to reflect our step size
            size_t step_size = last_mem_length - parent.lcp();
            match.end = mems.back().end-step_size;
            // and set up the next MEM using the parent node range
            match.range = parent.range();
        } else {
            // we are matching
            match.begin = cursor;
            // just step to the next position
            --cursor;
        }
    }
    // if we have a non-empty MEM at the end, record it
    if (match.end - match.begin > 0) mems.push_back(match);
    // return the matches in natural order
    std::reverse(mems.begin(), mems.end());
    // if debugging, verify the matches (costly)
    if (debug) { check_mems(mems); }
    return mems;
}

void Mapper::check_mems(const vector<MaximalExactMatch>& mems) {
    // for each copy of the mems
    for (auto mem : mems) {
        mem.fill_nodes(gcsa);
        if (sequence_positions(mem.sequence()) != gcsa_nodes_to_positions(mem.nodes)) {
            cerr << "mem failed! " << mem.sequence()
                 << " expected " << sequence_positions(mem.sequence()).size() << " hits "
                 << "but found " << gcsa_nodes_to_positions(mem.nodes).size() << endl;
        }
    }
}

// Use the GCSA2 index to find approximate maximal exact matches.
// When finding maximal exact matches (MEMs) we run backwards search in the index
// until we get a BWT range that's empty, or we exhaust the sequence or exceed the
// order of the index. The matches are approximately maximal due to these limitations.
// By converting the read to a series of MEMs we can leverage GCSA2 to determine the
// matching sequences in the read. Perfect matches generate MEMs. Interruptions due to
// errors and SNPs will tend to break MEMs, but may also generate shorter, spurious MEMs.
vector<MaximalExactMatch>
Mapper::find_forward_mems(const string& seq, size_t step) {
    string::const_iterator begin = seq.begin();
    string::const_iterator end = seq.end();
    vector<MaximalExactMatch> matches;
    // an empty sequence matches the entire bwt
    if(begin == end) {
        matches.emplace_back(MaximalExactMatch(begin, end, gcsa::range_type(0, gcsa->size() - 1)));
        return matches;
    }
    // otherwise, step through the sequence generating MEMs
    // running backward searching starting from end
    // and searching until we break, then emitting a match
    while (end != begin) {
        --end;
        matches.emplace_back(MaximalExactMatch(end, end+1,
                                               gcsa->charRange(gcsa->alpha.char2comp[*end])));
        MaximalExactMatch& match = matches.back();
        gcsa::range_type last_range = match.range;
        // run the backward search
        while (!gcsa::Range::empty(match.range) // until our BWT range is 0
               && end != begin                  // or we run out of sequence
               && (max_mem_length && match.end-end < max_mem_length || true)
               && match.end-end < gcsa->order()) { // or we exceed the order of the index
            --end;
            last_range = match.range;
            match.range = gcsa->LF(match.range, gcsa->alpha.char2comp[*end]);
        }
        // record the mem range in the sequence
        if (gcsa::Range::empty(match.range)) {
            match.range = last_range;
            match.begin = end+step;
        } else {
            // if we didn't exhaust the range, but ran out because our sequence ended
            // then we count this base in the MEM
            match.begin = end;
        }
    }
    // return the matches in natural order
    std::reverse(matches.begin(), matches.end());
    return matches;
}

const string mems_to_json(const vector<MaximalExactMatch>& mems) {
    stringstream s;
    s << "[";
    size_t j = 0;
    for (auto& mem : mems) {
        s << "[\"";
        s << mem.sequence();
        s << "\",[";
        size_t i = 0;
        for (auto& node : mem.nodes) {
            s << "\"" << gcsa::Node::decode(node) << "\"";
            if (++i < mem.nodes.size()) s << ",";
        }
        s << "]]";
        if (++j < mems.size()) s << ",";
    }
    s << "]";
    return s.str();
}

// these MEM walking functions are all rather inefficient
// for instance, we are getting the node each time; not good
// but my goal is to get them working before optimizing and this is clearer

// get the character at the given position in the graph
char Mapper::pos_char(pos_t pos) {
    Node node = xindex->node(id(pos)); // could we pass this in optionally to avoid reconstructing it?
    //cerr << "getting " << pos << " in " << pb2json(node) << endl;
    if (!is_rev(pos)) {
        return node.sequence().at(offset(pos));
    } else {
        return reverse_complement(node.sequence()).at(offset(pos));
    }
}

map<pos_t, char> Mapper::next_pos_chars(pos_t pos) {
    map<pos_t, char> nexts;
    Node node = xindex->node(id(pos));
    // if we are still in the node, return the next position and character
    if (offset(pos) < node.sequence().size()-1) {
        ++get_offset(pos);
        nexts[pos] = pos_char(pos);
    } else {
        // look at the next positions we could reach
        if (!is_rev(pos)) {
            // we are on the forward strand, the next things from this node come off the end
            for (auto& edge : xindex->edges_on_end(id(pos))) {
                if (edge.from() == id(pos)) {
                    pos_t p = make_pos_t(edge.to(), edge.to_end(), 0);
                    nexts[p] = pos_char(p);
                } else if (edge.from_start() && edge.to_end() && edge.to() == id(pos)) {
                    // doubly inverted, should be normalized to forward but we handle here for safety
                    pos_t p = make_pos_t(edge.from(), false, 0);
                    nexts[p] = pos_char(p);
                }
            }
        } else {
            // we are on the reverse strand, the next things from this node come off the start
            for (auto& edge : xindex->edges_on_start(id(pos))) {
                if (edge.to() == id(pos)) {
                    pos_t p = make_pos_t(edge.from(), !edge.from_start(), 0);
                    nexts[p] = pos_char(p);
                } else if (edge.from_start() && edge.to_end() && edge.from() == id(pos)) {
                    // doubly inverted, should be normalized to forward but we handle here for safety
                    pos_t p = make_pos_t(edge.to(), true, 0);
                    nexts[p] = pos_char(p);
                }
            }
        }
    }
    return nexts;
}

Alignment Mapper::walk_match(const string& seq, pos_t pos) {
    // go to the position in the xg index
    // and step in the direction given
    // until we exhaust our sequence
    // or hit another node
    Alignment aln; aln.set_sequence(seq);
    Path& path = *aln.mutable_path();
    Mapping* mapping = path.add_mapping();
    *mapping->mutable_position() = make_position(pos);
    // get the first node we match
    
    size_t match_len = 0;
    //cerr << "seq is " << seq << endl;
    for (size_t i = 0; i < seq.size()-1; ++i) {
        char c = seq[i];
        //cerr << pos_char(pos) << " == " << c << endl;
        auto nexts = next_pos_chars(pos);
        assert(nexts.size());
        if (nexts.size() == 1 && id(nexts.begin()->first) == id(pos)) {
            pos_t npos = nexts.begin()->first;
            //cerr << pos_char(npos) << " vs " << seq[i+1] << " i+1= " << i+1 << endl;
            assert(pos_char(npos) == seq[i+1]); // we must match
            ++match_len;
            ++get_offset(pos);
        } else { //if (nexts.size() > 1 || id(nexts.begin()->first) != id(pos)) {
            // we must be going into another node
            // emit the mapping for this node
            Edit* edit = mapping->add_edit();
            edit->set_from_length(match_len);
            edit->set_to_length(match_len);
            // find the next node that matches our MEM
            char n = seq[i+1];
            bool found_next = false;
            for (auto& p : nexts) {
                //cerr << "next : " << p.first << " " << p.second << endl;
                if (p.second == n) {
                    pos = p.first;
                    found_next = true;
                    break;
                }
            }
            assert(found_next);
            // set up a new mapping
            match_len = 0;
            mapping = path.add_mapping();
            *mapping->mutable_position() = make_position(pos);
        }
    }
    // add the last edit
    if (match_len) {
        Edit* edit = mapping->add_edit();
        edit->set_from_length(match_len);
        edit->set_to_length(match_len);
    }
    aln.set_score(aln.sequence().size()*match);
    return aln;
}

// convert one mem into a set of alignments, one for each exact match
vector<Alignment> Mapper::mem_to_alignments(MaximalExactMatch& mem) {
    vector<Alignment> alns;
    const string seq = mem.sequence();
    for (auto& node : mem.nodes) {
        pos_t pos = make_pos_t(gcsa::Node::id(node),
                               gcsa::Node::rc(node),
                               gcsa::Node::offset(node));
        alns.emplace_back(walk_match(seq, pos));
    }
    return alns;
}

vector<Alignment> Mapper::align_mem(const Alignment& alignment) {

    if (debug) cerr << "aligning " << pb2json(alignment) << endl;
    if (!gcsa || !xindex) {
        cerr << "error:[vg::Mapper] a GCSA2/xg index pair is required for MEM mapping" << endl;
        exit(1);
    }

    // use the GCSA index to determine the approximate MEMs of the read

    // remove those mems which are shorter than our kmer_min (could use a new parameter)
    // the short MEMs will have *many* hits, and we don't want to use them

    // for the moment, we can just get the positions in all the mems
    // but we may want to selectively use them
    // such as by trying to use them as seeds in order from shortest to longest
    vector<MaximalExactMatch> mems;
    for (auto& mem: find_smems(alignment.sequence())) {
        if (mem.end-mem.begin >= min_mem_length // require minimum length
            && allATGC(mem.sequence())) {       // and that the MEM does not contain Ns
            mems.push_back(mem);
        }
    }

    if (debug) cerr << "mems before filling " << mems_to_json(mems) << endl;
    for (auto& mem : mems) { get_mem_hits_if_under_max(mem); }
    if (debug) cerr << "mems after filtering " << mems_to_json(mems) << endl;

    return align_mem_multi(alignment, mems);

    /*
    // TODO
    if (max_multimaps > 1) {
        // align once per subgraph hit we get, take the best alignment
        // uses full-length alignment
        return align_mem_multi(alignment, mems);
    } else {
        // aligns only the regions of the read which aren't in the MEMs
        // stitching together the results using the resolve_banded_multi function
        // which presently only returns the optimal result
        return { align_mem_optimal(alignment, mems) };
    }
    */
}

bool Mapper::get_mem_hits_if_under_max(MaximalExactMatch& mem) {
    bool filled = false;
    // use the counting interface to determine the number of hits
    mem.fill_match_count(gcsa);
    // if we aren't filtering on hit count, or if we have up to the max allowed hits
    if (!hit_max || mem.match_count <= hit_max) {
        // extract the graph positions matching the range
        mem.fill_nodes(gcsa);
        filled = true;
    }
    return filled;
}


// Well, it's not really "optimal" yet
// but the idea here is that we can thread back through multi-mappings
// for each exact match
// use dynamic programming to fix up the intermediate portions
// and then thread back through with resolve_banded_multi
// VERY VERY EXPERIMENTAL
Alignment Mapper::align_mem_optimal(const Alignment& alignment, vector<MaximalExactMatch>& mems) {
    // run through the mems, generating a set of alignments for each
    vector<vector<Alignment> > alns;
    // run through the exact matches
    // if there are gaps between them, we need to fix them up using DP
    MaximalExactMatch* last_exact_match = nullptr;
    for (auto& mem : mems) {
        // get the >16mer centered on the non-matching portion
        if (last_exact_match) {
            // what is the gap?
            int gap_len = mem.begin - last_exact_match->end;
            if (gap_len > 0) {
                // stretch the gap
                // if in the middle of the read we get a 16-mer
                // if on the ends, we'll extend at least 8 opposite the end
                auto start = last_exact_match->end;
                auto end = mem.begin;
                int to_start = start - alignment.sequence().begin();
                int to_end = alignment.sequence().end() - end;
                int to_add = 16 - gap_len;
                int trim_from_start = 0;
                int trim_from_end = 0;
                if (to_add > 0) {
                    // how much can we do on the left
                    trim_from_start = min(to_add/2, to_start);
                    start -= trim_from_start;
                    // how much on the right
                    trim_from_end = min(to_add/2, to_end);
                    end += trim_from_end;
                }
                // get the sequence of the gap
                MaximalExactMatch gap(start, end, gcsa::range_type(1, 0));
                // now we need to align to all the possible paths from the ends of the last
                // set of alignments to the ends of the current set

                Alignment gapaln; gapaln.set_sequence(gap.sequence());
                vector<MaximalExactMatch> flanking_mems = { *last_exact_match, mem };
                int max_mem_gap_multi_mappings = 10;
                // then trim them using the trim counts
                auto malns = align_mem_multi(gapaln, flanking_mems, max_mem_gap_multi_mappings);
                vector<Alignment> trimmed;
                for (auto& aln : malns) {
                    //cerr << "aln " << pb2json(aln) << endl;
                    trimmed.push_back(simplify(strip_from_start(strip_from_end(aln, trim_from_end), trim_from_start)));
                    trimmed.back().set_score(0); // dunno
                    //cerr << "trimmed " << pb2json(trimmed.back()) << endl;
                }
                alns.emplace_back(trimmed);
            }
        }
        // add in this mem
        vector<Alignment> maln = mem_to_alignments(mem);
        alns.emplace_back(maln);
        last_exact_match = &mem;
    }
    auto optim_alns = resolve_banded_multi(alns);
    // merge the resulting alignments
    Alignment merged = merge_alignments(optim_alns, debug);
    merged.set_quality(alignment.quality());
    merged.set_name(alignment.name());
    return merged;
}

vector<Alignment> Mapper::align_mem_multi(const Alignment& alignment, vector<MaximalExactMatch>& mems, int max_multi) {

    struct StrandCounts {
        uint32_t forward;
        uint32_t reverse;
    };

    // we will use these to determine the alignment strand for each subgraph
    map<id_t, StrandCounts> node_strands;
    // records a mapping of id->MEMs, for cluster ranking
    map<id_t, vector<MaximalExactMatch> > id_to_mems;
    // for clustering
    set<id_t> ids;

    // run through the mems, generating a set of alignments for each
    for (auto& mem : mems) {
        size_t len = mem.begin - mem.end;
        // collect ids and orienations of hits to them on the forward mem
        for (auto& node : mem.nodes) {
            id_t id = gcsa::Node::id(node);
            id_to_mems[id].push_back(mem); // copy
            ids.insert(id);
            if (gcsa::Node::rc(node)) {
                node_strands[id].reverse++;
            } else {
                node_strands[id].forward++;
            }
        }

        // now handle the reverse complement mem to get the end position of the mem
        string rc_mem = reverse_complement(mem.sequence());
        auto rc_mems = find_smems(rc_mem);
        // also filter to remove spurious sub-hits
        // TODO: is it the case that the rc is definitively the same number as the forward?
        rc_mems.erase(std::remove_if(rc_mems.begin(), rc_mems.end(),
                                  [&](const MaximalExactMatch& m) { return m.end-m.begin < min_mem_length; }));
        for (auto& rc_mem : rc_mems) {
            get_mem_hits_if_under_max(rc_mem);
            for (auto& node : rc_mem.nodes) {
                id_t id = gcsa::Node::id(node);
                id_to_mems[id].push_back(rc_mem); // copy
                ids.insert(id);
                // we are looking at the position of the end of the mem on the reverse strand
                // so we flip it around when tabulating strands
                if (!gcsa::Node::rc(node)) {
                    node_strands[id].reverse++;
                } else {
                    node_strands[id].forward++;
                }
            }
        }
    }

    // collect the graph implied by the mems and their reverse complements
    // attempt to pick up ranges between successive nodes
    // when these are below our context_depth
    vector<vector<id_t> > clusters;
    for (auto& id : ids) {
        if (clusters.empty()) {
            clusters.emplace_back();
            auto& l = clusters.back();
            l.push_back(id);
        } else {
            auto& prev = clusters.back().back();
            if (id - prev <= thread_extension) {
                clusters.back().push_back(id);
            } else {
                clusters.emplace_back();
                auto& l = clusters.back();
                l.push_back(id);
            }
        }
    }

    // get target subgraphs
    vector<VG> subgraphs;
    int max_target_length = alignment.sequence().size() * max_target_factor;
    if (debug) cerr << clusters.size() << " clusters" << endl;
    for (auto& cluster : clusters) {
        subgraphs.emplace_back();
        auto& sub = subgraphs.back();
        auto id1 = cluster.front();
        auto id2 = cluster.back();
        xindex->get_id_range(id1, id2, sub.graph);
        xindex->expand_context(sub.graph, context_depth, false);
        sub.rebuild_indexes();
        if (max_target_factor && sub.length() > max_target_length) {
            subgraphs.pop_back();
        }
    }

    // order the subgraphs by number of hits
    // and go through them until we have enough multi-maps
    map<double, vector<VG*> > subgraphs_by_hits;;
    for (auto& subgraph : subgraphs) {
        double hit_length = 0;
        subgraph.for_each_node([&](Node* n) {
                auto m = id_to_mems.find(n->id());
                if (m != id_to_mems.end()) {
                    for (auto& mem : m->second) {
                        // we add weight for the hit length, but
                        // weight the hit by the number of matches
                        hit_length +=
                            (double)(mem.end - mem.begin)
                            /(double)mem.match_count;
                    }
                }
            });
        subgraphs_by_hits[hit_length].push_back(&subgraph);
    }

    // sort the ranked subgraphs and only keep the best N (attempts)
    vector<VG*> ranked_subgraphs;
    size_t i = 0;
    for (auto it = subgraphs_by_hits.rbegin();
         it != subgraphs_by_hits.rend()
             // if we've set max_attempts, only try up to that many graphs
             && (!max_attempts || i < max_attempts);
         ++it, ++i) {
        for (auto subgraph : it->second) {
            ranked_subgraphs.push_back(subgraph);
        }
    }

    // generate an alignment for each subgraph/orientation combination for which we have hits
    if (debug) cerr << "aligning to " << subgraphs.size() << " subgraphs" << endl;

    vector<Alignment> alns; // our alignments
    
    // set up our forward and reverse base alignments (these are just sequences in bare alignment objs)
    auto aln_fw = alignment;
    aln_fw.clear_path();
    aln_fw.set_score(0);
    auto aln_rc = aln_fw;
    aln_rc.set_sequence(reverse_complement(aln_fw.sequence()));

    size_t attempts = 0;
    for (auto s : ranked_subgraphs) {
        VG& subgraph = *s;
        // record our attempt count
        ++attempts;
        // determine the likely orientation
        uint32_t fw_mems = 0;
        uint32_t rc_mems = 0;
        subgraph.for_each_node([&](Node* n) {
                auto ns = node_strands.find(n->id());
                if (ns != node_strands.end()) {
                    fw_mems += ns->second.forward;
                    rc_mems += ns->second.reverse;
                }
            });
        if (debug) cerr << "got " << fw_mems << " forward and " << rc_mems << " reverse mems" << endl;
        if (fw_mems) {
            Alignment aln = subgraph.align(aln_fw, match, mismatch, gap_open, gap_extension);
            resolve_softclips(aln, subgraph);
            alns.push_back(aln);
            if (attempts >= max_multimaps &&
                greedy_accept && (float)aln.score() / ((float)aln.sequence().size()*match) >= accept_norm_score) {
                break;
            }
        }
        if (rc_mems) {
            Alignment aln = subgraph.align(aln_rc, match, mismatch, gap_open, gap_extension);
            resolve_softclips(aln, subgraph);
            alns.push_back(reverse_complement_alignment(aln,
                                                        (function<int64_t(int64_t)>)
                                                        ([&](int64_t id) { return get_node_length(id); })));
            if (attempts >= max_multimaps &&
                greedy_accept && (float)aln.score() / ((float)aln.sequence().size()*match) >= accept_norm_score) {
                break;
            }
        }
    }

    // free memory; we don't need the graphs for the last bit
    subgraphs.clear();

    int sum_score = 0;
    double mean_score = 0;
    map<int, set<Alignment*> > alignment_by_score;
    for (auto& ta : alns) {
        Alignment* aln = &ta;
        alignment_by_score[aln->score()].insert(aln);
    }

    // set the max multimap count from the mapper object if we haven't been given it
    if (!max_multi) max_multi = max_multimaps;

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

            if(good.size() >= max_multi) {
                // Don't report too many mappings
                break;
            }
        }

        if(good.size() >= max_multi) {
            // Don't report too many mappings
            break;
        }
    }

    // get the best alignment
    if (!good.empty()) {
        // TODO log best alignment score?
    } else {
        good.emplace_back();
        Alignment& aln = good.back();
        aln = alignment;
        aln.clear_path();
        aln.set_score(0);
    }

    // Return all the multimappings
    return good;
    
}

void Mapper::resolve_softclips(Alignment& aln, VG& graph) {

    if (!xindex) {
        cerr << "error:[vg::Mapper] xg index pair is required for dynamic softclip resolution" << endl;
        exit(1);
    }
    // we can't resolve softclips on a read without a mapping
    if (!aln.path().mapping_size()) return;
    // we can be more precise about our handling of softclips due to the low cost
    // of the fully in-memory xg index
    int sc_start = softclip_start(aln);
    int sc_end = softclip_end(aln);
    int last_score = aln.score();
    size_t itr = 0;
    Path* path = aln.mutable_path();
    int64_t idf = path->mutable_mapping(0)->position().node_id();
    int64_t idl = path->mutable_mapping(path->mapping_size()-1)->position().node_id();
    int32_t d_to_head = graph.distance_to_head(NodeTraversal(graph.get_node(idf), false), sc_start*3);
    int32_t d_to_tail = graph.distance_to_tail(NodeTraversal(graph.get_node(idl), false), sc_end*3);
    int max_target_length = aln.sequence().size() * max_target_factor;
    while (itr++ < 3
           && ((sc_start > softclip_threshold
                && d_to_head >= 0 && d_to_head < sc_start)
               || (sc_end > softclip_threshold
                   && d_to_tail >=0 && d_to_tail < sc_end))) {
        if (debug) {
            cerr << "softclip before " << sc_start << " " << sc_end << endl;
            cerr << "distance to head "
                 << graph.distance_to_head(NodeTraversal(graph.get_node(idf), false), sc_start*3)
                 << endl;
            cerr << "distance to tail "
                 << graph.distance_to_tail(NodeTraversal(graph.get_node(idl), false), sc_end*3)
                 << endl;
        }
        double avg_node_size = graph.length() / graph.size();
        if (debug) cerr << "average node size " << avg_node_size << endl;
        // step towards the side where there were soft clips
        if (sc_start) {
            Graph flank;
            xindex->get_id_range(idf-1, idf, flank);
            xindex->expand_context(flank,
                                   max(context_depth, (int)(sc_start/avg_node_size)),
                                   false);
            graph.extend(flank);
        }
        if (sc_end) {
            Graph flank;
            xindex->get_id_range(idl, idl+1, flank);
            xindex->expand_context(flank,
                                   max(context_depth, (int)(sc_end/avg_node_size)),
                                   false);
            graph.extend(flank);
        }
        graph.remove_orphan_edges();
        aln.clear_path();
        aln.set_score(0);

        // give up if the graph is too big
        if (max_target_factor && graph.length() >= max_target_length) break;

        // otherwise, align
        aln = graph.align(aln, match, mismatch, gap_open, gap_extension);

        sc_start = softclip_start(aln);
        sc_end = softclip_end(aln);
        if (debug) cerr << "softclip after " << sc_start << " " << sc_end << endl;
        // we are not improving, so increasing the window is unlikely to help
        if (last_score == aln.score()) break;
        // update tracking of path end
        last_score = aln.score();
        path = aln.mutable_path();
        idf = path->mutable_mapping(0)->position().node_id();
        idl = path->mutable_mapping(path->mapping_size()-1)->position().node_id();
        d_to_head = graph.distance_to_head(NodeTraversal(graph.get_node(idf), false), sc_start*3);
        d_to_tail = graph.distance_to_tail(NodeTraversal(graph.get_node(idl), false), sc_end*3);
    }
}

// core alignment algorithm that handles both kinds of sequence indexes
vector<Alignment> Mapper::align_threaded(const Alignment& alignment, int& kmer_count, int kmer_size, int stride, int attempt) {

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
        // we can do multiple mapping by rnelaxing this
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
            // due to the fact that we record the node position of the start of the kmer
            int64_t first = max((int64_t)0, *thread.begin());
            int64_t last = *thread.rbegin() + thread_ex;
            // so we can pick it up efficiently from the index by pulling the range from first to last
            if (debug) cerr << "getting node range " << first << "-" << last << endl;
            VG* graph = new VG;

            // Now we need to get the neighborhood by ID and expand outward by actual
            // edges. How we do this depends on what indexing structures we have.
            // TODO: We're repeating this code. Break it out into a function or something.
            if(xindex) {
                xindex->get_id_range(first, last, graph->graph);
                xindex->expand_context(graph->graph, context_depth, false);
                graph->rebuild_indexes();
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
            //serialize_to_file("init-" + alignment.sequence() + "-" + hash_alignment(alignment).substr(0,8) + "-" + hash().substr(0,8) + ".vg");

            // Topologically sort the graph, breaking cycles and orienting all edges end to start.
            // This flips some nodes around, so we need to translate alignments back.
            //set<int64_t> flipped_nodes;
            //graph->orient_nodes_forward(flipped_nodes);

            // align
            //graph->serialize_to_file("align2.vg");
            ta.clear_path();
            ta.set_score(0);

            ta = graph->align(ta, match, mismatch, gap_open, gap_extension);

            // check if we start or end with soft clips
            // if so, try to expand the graph until we don't have any more (or we hit a threshold)
            // expand in the direction where there were soft clips

            if (!ta.has_path()) continue;

            // we can be more precise about our handling of softclips due to the low cost
            // of the fully in-memory xg index
            int sc_start = softclip_start(ta);
            int sc_end = softclip_end(ta);
            int last_score = ta.score();
            size_t itr = 0;
            Path* path = ta.mutable_path();
            int64_t idf = path->mutable_mapping(0)->position().node_id();
            int64_t idl = path->mutable_mapping(path->mapping_size()-1)->position().node_id();
            int32_t d_to_head = graph->distance_to_head(NodeTraversal(graph->get_node(idf), false), sc_start*3);
            int32_t d_to_tail = graph->distance_to_tail(NodeTraversal(graph->get_node(idl), false), sc_end*3);
            while (itr++ < 3
                   && ((sc_start > softclip_threshold
                        && d_to_head >= 0 && d_to_head < sc_start)
                       || (sc_end > softclip_threshold
                           && d_to_tail >=0 && d_to_tail < sc_end))) {
                if (debug) {
                    cerr << "softclip before " << sc_start << " " << sc_end << endl;
                    cerr << "distance to head "
                         << graph->distance_to_head(NodeTraversal(graph->get_node(idf), false), sc_start*3)
                         << endl;
                    cerr << "distance to tail "
                         << graph->distance_to_tail(NodeTraversal(graph->get_node(idl), false), sc_end*3)
                         << endl;
                }
                double avg_node_size = graph->length() / graph->size();
                if (debug) cerr << "average node size " << avg_node_size << endl;
                // step towards the side where there were soft clips
                if (sc_start) {
                    if (xindex) {
                        Graph flank;
                        xindex->get_id_range(idf-1, idf, flank);
                        xindex->expand_context(flank,
                                               max(context_depth, (int)(sc_start/avg_node_size)),
                                               false);
                        graph->extend(flank);
                    } else if (index) {
                        VG flank;
                        index->get_range(max((int64_t)0, idf-thread_ex), idf, flank);
                        index->expand_context(flank, context_depth);
                        graph->extend(flank);
                    }
                }
                if (sc_end) {
                    if (xindex) {
                        Graph flank;
                        xindex->get_id_range(idl, idl+1, flank);
                        xindex->expand_context(flank,
                                               max(context_depth, (int)(sc_end/avg_node_size)),
                                               false);
                        graph->extend(flank);
                    } else if (index) {
                        VG flank;
                        index->get_range(idl, idl+thread_ex, flank);
                        index->expand_context(flank, context_depth);
                        graph->extend(flank);
                    }
                }
                graph->remove_orphan_edges();
                ta.clear_path();
                ta.set_score(0);

                ta = graph->align(ta, match, mismatch, gap_open, gap_extension);

                sc_start = softclip_start(ta);
                sc_end = softclip_end(ta);
                if (debug) cerr << "softclip after " << sc_start << " " << sc_end << endl;
                // we are not improving, so increasing the window is unlikely to help
                if (last_score == ta.score()) break;
                // update tracking of path end
                last_score = ta.score();
                path = ta.mutable_path();
                idf = path->mutable_mapping(0)->position().node_id();
                idl = path->mutable_mapping(path->mapping_size()-1)->position().node_id();
                d_to_head = graph->distance_to_head(NodeTraversal(graph->get_node(idf), false), sc_start*3);
                d_to_tail = graph->distance_to_tail(NodeTraversal(graph->get_node(idl), false), sc_end*3);
            }

            delete graph;

            if (debug) cerr << "normalized score is " << (float)ta.score() / ((float)ta.sequence().size()*match) << endl;
            if (greedy_accept && (float)ta.score() / ((float)ta.sequence().size()*match) >= accept_norm_score) {
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
        // TODO log best alignment score?
    } else {
        good.emplace_back();
        Alignment& aln = good.back();
        aln = alignment;
        aln.clear_path();
        aln.set_score(0);
    }

    // Return all the multimappings
    return good;
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
